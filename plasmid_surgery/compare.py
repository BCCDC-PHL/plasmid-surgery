"""
All-vs-all comparison of plasmid reconstructions (and references) for the cohort report.

Every member's contigs are aligned to every other's in one minimap2 call (asm20, all chains kept).
For each pair (A, B):

- **A in B**: the fraction of A's bases covered by its alignments to B (merged per contig; the inserted
  bases of an indel of at least MAX_INDEL bp within an alignment don't count), so a
  fragmentary short-read reconstruction contained in a complete plasmid shows A in B near 1 and B in A
  lower;
- **identity**: gap-compressed identity over their alignments (1 - minimap2's ``de``: an indel counts as
  one difference, not one per base, so a large insertion doesn't make two otherwise identical plasmids
  look divergent), weighted by alignment length.

Alignments shorter than MIN_ALIGNMENT bp or below MIN_IDENTITY gap-compressed identity don't count. Members are ordered
by average-linkage clustering on 1 - mean containment, so similar plasmids sit together.

**SNP distances** (pair_snps), from the same alignments:

- **one-to-one**: each stretch of either member is compared once (alignments kept best first, skipping any
  overlapping a kept one by more than OVERLAP_TOLERANCE bp on either side), and stretches aligning more
  than once between the pair (repeats, e.g. IS copies) are masked, so paralogous copies aren't compared;
- **trimmed**: differences within END_TRIM bp of an alignment's or a contig's end don't count (assembly
  errors concentrate at contig ends);
- substitutions (**SNPs**) and **indel** events are counted separately (long-read assembly errors are
  mostly homopolymer indels);
- SNPs in **clusters** (at least CLUSTER_MIN within CLUSTER_WINDOW bp) are counted apart: a block from a
  different lineage brings many SNPs at once;
- **compared bases** are the aligned positions that count, so a SNP count is always read with how much
  was compared.
"""
from __future__ import annotations

import bisect
import csv
import logging
import re
import subprocess

from pathlib import Path

import plasmid_surgery.alignment as alignment

logger = logging.getLogger(__name__)

MIN_ALIGNMENT = 500
MIN_IDENTITY = 0.9
MAX_INDEL = 50          # indels at least this long within an alignment aren't covered (e.g. a novel gene)
OVERLAP_TOLERANCE = 50  # one-to-one: alignments may overlap kept ones by this much
END_TRIM = 100          # differences this close to an alignment's or contig's end don't count
CLUSTER_WINDOW = 1000
CLUSTER_MIN = 3
MIN_COMPARED_FRACTION = 0.5   # a SNP distance is shown if at least this fraction of the shorter member was compared

CS_RE = re.compile(r'(:\d+|\*[a-z][a-z]|[+-][a-z]+)')

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

PAIR_FIELDNAMES = ['a', 'b', 'a_length', 'b_length', 'a_in_b', 'b_in_a', 'identity', 'aligned_bases_a', 'aligned_bases_b',
                   'compared_bases', 'compared_fraction_a', 'compared_fraction_b', 'snps', 'snps_in_clusters',
                   'snps_outside_clusters', 'snps_per_100kb', 'indels', 'snp_distance_reported']


def _merged_length(intervals) -> int:
    total, end = 0, None
    for s, e in sorted(intervals):
        if end is None or s > end:
            total += e - s
            end = e
        elif e > end:
            total += e - end
            end = e
    return total


def covered_blocks(f: list[str]) -> tuple[list, list]:
    """
    From a PAF line (with cg:Z), the query and target intervals the alignment covers, split at indels of at
    least MAX_INDEL bp (whose inserted bases aren't covered).
    """
    qs, qe, strand, ts = int(f[2]), int(f[3]), f[4], int(f[7])
    cigar = next((t[5:] for t in f[12:] if t.startswith('cg:Z:')), None)
    if cigar is None:
        return [(qs, qe)], [(ts, int(f[8]))]
    q_blocks, t_blocks = [], []
    q, t = 0, ts                       # q: bases of the query consumed, along the target's direction
    q_start, t_start = 0, ts
    for n, op in CIGAR_RE.findall(cigar):
        n = int(n)
        if op in 'M=X':
            q += n
            t += n
        elif op in 'IDN' and n >= MAX_INDEL:
            if q > q_start or t > t_start:
                q_blocks.append((q_start, q))
                t_blocks.append((t_start, t))
            if op == 'I':
                q += n
            else:
                t += n
            q_start, t_start = q, t
        elif op == 'I':
            q += n
        elif op in 'DN':
            t += n
    if q > q_start or t > t_start:
        q_blocks.append((q_start, q))
        t_blocks.append((t_start, t))
    if strand == '+':
        q_blocks = [(qs + a, qs + b) for a, b in q_blocks]
    else:
        q_blocks = [(qe - b, qe - a) for a, b in q_blocks]
    return q_blocks, t_blocks


def all_vs_all(members: list[dict], tmpdir, threads: int = 1) -> dict:
    """
    :param members: [{'name': str, 'contigs': [{'seq_id', 'seq'}, ...]}]
    :return: dict with 'names', 'lengths', 'containment' (containment[i][j]: fraction of member i in member j),
             'aligned' (aligned[i][j]: bases of i covered by j), 'identity' (identity[i][j], symmetric; None if
             no alignment) and 'order' (member indices, clustered)
    :raises alignment.AlignmentError: if minimap2 is missing or fails
    """
    tmpdir = Path(tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)
    fasta = tmpdir / 'members.fa'
    records = []
    for i, m in enumerate(members):
        for k, c in enumerate(m['contigs']):
            records.append((f'm{i}_{k}', c['seq']))
    alignment.write_fasta(fasta, records)
    n = len(members)
    lengths = [sum(len(c['seq']) for c in m['contigs']) for m in members]
    covered = {}                       # (i, j) -> {contig k: [(start, end)]}: parts of i aligning to j
    pair_alignments = {}               # (i, j), i < j -> alignments of i (query) to j (target), for SNPs
    matches = [[0] * n for _ in range(n)]
    blocks = [[0] * n for _ in range(n)]
    cmd = ['minimap2', '-c', '--cs', '-x', 'asm20', '-P', '-t', str(threads), str(fasta), str(fasta)]
    logger.info(f"Comparing {n} sequences all-vs-all")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError as e:
        raise alignment.AlignmentError("minimap2 not found on PATH") from e
    except subprocess.CalledProcessError as e:
        raise alignment.AlignmentError(f"minimap2 failed comparing the reconstructions:\n{e.stderr[-2000:]}")
    for line in result.stdout.splitlines():
        f = line.split('\t')
        qi, qk = map(int, f[0][1:].split('_'))
        ti, tk = map(int, f[5][1:].split('_'))
        if qi == ti:
            continue
        block = int(f[10])
        de = next((float(t[5:]) for t in f[12:] if t.startswith('de:f:')), None)
        ident = 1 - de if de is not None else int(f[9]) / block if block else 0
        if block < MIN_ALIGNMENT or ident < MIN_IDENTITY:
            continue
        q_blocks, t_blocks = covered_blocks(f)
        covered.setdefault((qi, ti), {}).setdefault(qk, []).extend(q_blocks)
        covered.setdefault((ti, qi), {}).setdefault(tk, []).extend(t_blocks)
        a, b = min(qi, ti), max(qi, ti)
        matches[a][b] += ident * block
        blocks[a][b] += block
        if qi < ti:
            cs = next((t[5:] for t in f[12:] if t.startswith('cs:Z:')), None)
            if cs is not None:
                pair_alignments.setdefault((qi, ti), []).append({
                    'qk': qk, 'qlen': int(f[1]), 'qs': int(f[2]), 'qe': int(f[3]), 'strand': f[4],
                    'tk': tk, 'tlen': int(f[6]), 'ts': int(f[7]), 'te': int(f[8]), 'block': block, 'cs': cs})
    aligned = [[0] * n for _ in range(n)]
    containment = [[0.0] * n for _ in range(n)]
    identity = [[None] * n for _ in range(n)]
    for i in range(n):
        containment[i][i] = 1.0
        for j in range(n):
            if i == j:
                continue
            aligned[i][j] = sum(_merged_length(ivs) for ivs in covered.get((i, j), {}).values())
            containment[i][j] = round(aligned[i][j] / lengths[i], 4) if lengths[i] else 0.0
            a, b = min(i, j), max(i, j)
            if blocks[a][b]:
                identity[i][j] = round(matches[a][b] / blocks[a][b], 5)
    distance = [[1 - (containment[i][j] + containment[j][i]) / 2 for j in range(n)] for i in range(n)]
    snps = {}
    for i in range(n):
        for j in range(i + 1, n):
            snps[(i, j)] = pair_snps(pair_alignments.get((i, j), []), lengths[i], lengths[j])
    return {'names': [m['name'] for m in members], 'lengths': lengths, 'containment': containment,
            'aligned': aligned, 'identity': identity, 'snps': snps, 'order': average_linkage_order(distance)}


# --- SNP distances ---

def _merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if e <= s:
            continue
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def _multiply_covered(intervals) -> list:
    """Where at least two of the intervals overlap."""
    events = sorted([(s, 1) for s, e in intervals] + [(e, -1) for s, e in intervals])
    out, depth, start = [], 0, None
    for pos, d in events:
        depth += d
        if depth >= 2 and start is None:
            start = pos
        elif depth < 2 and start is not None:
            if pos > start:
                out.append((start, pos))
            start = None
    return _merge(out)


def _overlap(a, b) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def one_to_one(alignments: list[dict]) -> list[dict]:
    """Alignments kept best (longest) first, skipping any overlapping a kept one by more than OVERLAP_TOLERANCE
    bp on the query or the target."""
    kept = []
    for a in sorted(alignments, key=lambda x: -x['block']):
        if all(not (a['qk'] == k['qk'] and _overlap((a['qs'], a['qe']), (k['qs'], k['qe'])) > OVERLAP_TOLERANCE)
               and not (a['tk'] == k['tk'] and _overlap((a['ts'], a['te']), (k['ts'], k['te'])) > OVERLAP_TOLERANCE)
               for k in kept):
            kept.append(a)
    return kept


class _Excluded:
    """Merged excluded intervals on one contig, with a fast membership test."""
    def __init__(self, intervals):
        self.ivs = _merge(intervals)
        self.starts = [s for s, _ in self.ivs]

    def __contains__(self, p):
        i = bisect.bisect_right(self.starts, p) - 1
        return i >= 0 and p < self.ivs[i][1]

    def within(self, start, end):
        """The excluded intervals clipped to [start, end)."""
        i = max(bisect.bisect_right(self.starts, start) - 1, 0)
        out = []
        while i < len(self.ivs) and self.ivs[i][0] < end:
            s, e = max(self.ivs[i][0], start), min(self.ivs[i][1], end)
            if e > s:
                out.append((s, e))
            i += 1
        return out


def clustered(positions: list[int]) -> set:
    """Indices of the SNP positions (sorted) lying in a window of CLUSTER_WINDOW bp holding at least CLUSTER_MIN."""
    out, j = set(), 0
    for i in range(len(positions)):
        while positions[i] - positions[j] >= CLUSTER_WINDOW:
            j += 1
        if i - j + 1 >= CLUSTER_MIN:
            out.update(range(j, i + 1))
    return out


def pair_snps(alignments: list[dict], length_a: int, length_b: int) -> dict:
    """
    SNPs and indels between members A (query) and B (target) from their alignments (with cs tags): one-to-one,
    repeats masked, ends trimmed (see the module docstring).

    :return: dict with 'compared_bases', 'snps', 'snps_in_clusters', 'indels'
    """
    # Repeats: stretches of either member aligned more than once between the pair
    q_mask, t_mask = {}, {}
    for a in alignments:
        q_mask.setdefault(a['qk'], []).append((a['qs'], a['qe']))
        t_mask.setdefault(a['tk'], []).append((a['ts'], a['te']))
    q_mask = {k: _multiply_covered(v) for k, v in q_mask.items()}
    t_mask = {k: _multiply_covered(v) for k, v in t_mask.items()}
    compared, indels = 0, 0
    snp_positions = {}
    for a in one_to_one(alignments):
        ts, te, qs, qe = a['ts'], a['te'], a['qs'], a['qe']
        t_excl = _Excluded(t_mask.get(a['tk'], []) + [(ts, ts + END_TRIM), (te - END_TRIM, te),
                                                      (0, END_TRIM), (a['tlen'] - END_TRIM, a['tlen'])])
        q_excl = _Excluded(q_mask.get(a['qk'], []) + [(0, END_TRIM), (a['qlen'] - END_TRIM, a['qlen'])])
        forward = a['strand'] == '+'
        t, qoff = ts, 0

        def qpos(off):
            return qs + off if forward else qe - 1 - off
        for op in CS_RE.findall(a['cs']):
            kind = op[0]
            if kind == ':':
                n = int(op[1:])
                # Excluded positions of this run, as offsets 0..n: the target's, and the query's mapped onto it
                parts = [(s - t, e - t) for s, e in t_excl.within(t, t + n)]
                if forward:
                    parts += [(s - qpos(qoff), e - qpos(qoff)) for s, e in q_excl.within(qpos(qoff), qpos(qoff) + n)]
                else:
                    lo, hi = qpos(qoff + n - 1), qpos(qoff) + 1
                    parts += [(qpos(qoff) - (e - 1), qpos(qoff) - s + 1) for s, e in q_excl.within(lo, hi)]
                compared += n - sum(e - s for s, e in _merge(parts))
                t += n
                qoff += n
            elif kind == '*':
                if t not in t_excl and qpos(qoff) not in q_excl:
                    compared += 1
                    snp_positions.setdefault(a['tk'], []).append(t)
                t += 1
                qoff += 1
            elif kind == '+':
                if t not in t_excl and qpos(qoff) not in q_excl:
                    indels += 1
                qoff += len(op) - 1
            else:   # '-'
                if t not in t_excl and qpos(qoff) not in q_excl:
                    indels += 1
                t += len(op) - 1
    snps = sum(len(v) for v in snp_positions.values())
    in_clusters = sum(len(clustered(sorted(v))) for v in snp_positions.values())
    return {'compared_bases': compared, 'snps': snps, 'snps_in_clusters': in_clusters, 'indels': indels}


def snp_distance_reported(compared: int, length_a: int, length_b: int) -> bool:
    """Whether enough was compared for the SNP count to stand as a distance."""
    shorter = min(length_a, length_b)
    return shorter > 0 and compared >= MIN_COMPARED_FRACTION * shorter


def average_linkage_order(distance: list[list[float]]) -> list[int]:
    """Leaf order of an average-linkage (UPGMA) clustering: merges the closest clusters first."""
    clusters = {i: [i] for i in range(len(distance))}
    while len(clusters) > 1:
        keys = sorted(clusters)
        best = None
        for x, a in enumerate(keys):
            for b in keys[x + 1:]:
                d = sum(distance[i][j] for i in clusters[a] for j in clusters[b]) / (len(clusters[a]) * len(clusters[b]))
                if best is None or d < best[0]:
                    best = (d, a, b)
        _, a, b = best
        clusters[a] = clusters[a] + clusters.pop(b)
    return next(iter(clusters.values()), [])


def pair_rows(result: dict) -> list[dict]:
    """One row per unordered pair (PAIR_FIELDNAMES)."""
    names, n = result['names'], len(result['names'])
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            rows.append({'a': names[i], 'b': names[j], 'a_length': result['lengths'][i], 'b_length': result['lengths'][j],
                         'a_in_b': result['containment'][i][j], 'b_in_a': result['containment'][j][i],
                         'identity': result['identity'][i][j] if result['identity'][i][j] is not None else 'NA',
                         'aligned_bases_a': result['aligned'][i][j], 'aligned_bases_b': result['aligned'][j][i],
                         **snp_fields(result, i, j)})
    return rows


def snp_fields(result: dict, i: int, j: int) -> dict:
    """The SNP columns of pair (i, j), i < j."""
    snp = result['snps'][(min(i, j), max(i, j))]
    la, lb = result['lengths'][i], result['lengths'][j]
    compared = snp['compared_bases']
    outside = snp['snps'] - snp['snps_in_clusters']
    return {'compared_bases': compared,
            'compared_fraction_a': round(compared / la, 4) if la else 0.0,
            'compared_fraction_b': round(compared / lb, 4) if lb else 0.0,
            'snps': snp['snps'], 'snps_in_clusters': snp['snps_in_clusters'], 'snps_outside_clusters': outside,
            'snps_per_100kb': round(outside / compared * 1e5, 2) if compared else 'NA',
            'indels': snp['indels'], 'snp_distance_reported': snp_distance_reported(compared, la, lb)}


def write_pairs(path, result: dict) -> None:
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=PAIR_FIELDNAMES, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(pair_rows(result))
