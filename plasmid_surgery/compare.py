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
"""
from __future__ import annotations

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

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

PAIR_FIELDNAMES = ['a', 'b', 'a_length', 'b_length', 'a_in_b', 'b_in_a', 'identity', 'aligned_bases_a', 'aligned_bases_b']


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
    matches = [[0] * n for _ in range(n)]
    blocks = [[0] * n for _ in range(n)]
    cmd = ['minimap2', '-c', '-x', 'asm20', '-P', '-t', str(threads), str(fasta), str(fasta)]
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
    return {'names': [m['name'] for m in members], 'lengths': lengths, 'containment': containment,
            'aligned': aligned, 'identity': identity, 'order': average_linkage_order(distance)}


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
                         'aligned_bases_a': result['aligned'][i][j], 'aligned_bases_b': result['aligned'][j][i]})
    return rows


def write_pairs(path, result: dict) -> None:
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=PAIR_FIELDNAMES, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(pair_rows(result))
