"""
Joins: where a reconstruction switches from one reference's sequence to another's (or to sequence no
reference explains), and whether long reads run across those points.

A reconstruction may join sequence from more than one plasmid: two plasmids in one bin, or one contig
joining two plasmids, which may be a real cointegrate or a misassembly (typically through a repeat both
carry, such as an IS element). The references assessed show where: the **composition** picks references
one at a time by how much new sequence each explains, and labels each stretch of each contig with the
first reference (in that order) that explains it, or none. A **join** is where the label changes
between stretches of at least ``min_block`` bp; where the two references overlap there (e.g. a shared
IS), the join covers the overlap. Each circular contig's own junction (end to start) is checked too.
Each join is then widened to cover any repeat of the reconstruction touching it (from aligning the
reconstruction to itself): a misjoin happens through such a repeat, and reads from either plasmid run
through the repeat, so only reads anchored in unique sequence on both sides show the join is real.

Read evidence at each join, from the reads aligned to the reconstruction's contigs (circular contigs
extended by their own start, so reads can cross the junction):

- long reads **spanning** it: one continuous alignment covering the join plus ``anchor`` bp each side;
- long reads **clipped** at it: an alignment ending (or starting) at the join with at least CLIP_MIN bp
  of the read unaligned beyond it, i.e. the read continues somewhere else;
- read depth in the flanks either side (short and long reads separately), to show a step.

Calls: ``supported`` (at least MIN_SPANNING spanning reads), ``contradicted`` (none spanning, at least
MIN_CLIPPED clipped), ``weak`` (otherwise), ``no_long_reads``.

Reads are aligned to the reconstruction only, so reads from elsewhere in the genome that share sequence
with it (e.g. the chromosome's copies of an IS) pile onto that sequence; they don't carry the unique
flanks, so they don't count as spanning, and the flank depths are taken outside the join.
"""
from __future__ import annotations

import csv
import gzip
import json
import logging
import os
import re
import subprocess
import tempfile

from array import array
from itertools import accumulate
from pathlib import Path
from statistics import median

import plasmid_surgery.alignment as alignment
import plasmid_surgery.assessment as assessment
import plasmid_surgery.reads as reads

logger = logging.getLogger(__name__)

DEFAULT_MIN_BLOCK = 1000      # shorter alignments / stretches don't make joins (e.g. scattered IS hits)
DEFAULT_ANCHOR = 1000         # bp a spanning read must align beyond each side of a join
CIRCULAR_EXTENSION = 30000    # circular contigs are extended by up to this much of their own start
MIN_SPANNING = 3
MIN_CLIPPED = 3
CLIP_MIN = 200                # bp of a read left unaligned beyond its alignment's end to count as clipped
CLIP_WINDOW = 100             # an alignment end within this of a join's edge counts as at the join
FLANK = 5000                  # bp each side of a join for its depth
CACHE_FORMAT_VERSION = 1
MIN_REPEAT = 200              # self-alignments at least this long count as repeats
REPEAT_TOLERANCE = 50         # a repeat within this of a join touches it

CIRCULAR_RE = re.compile(r'circular\s*=\s*(true|y|yes)\b|topology\s*=\s*circular', re.IGNORECASE)

COMPOSITION_FIELDNAMES = ['order', 'ref_id', 'new_bases', 'new_fraction', 'cumulative_fraction']

JOINS_FIELDNAMES = [
    'contig_id',
    'join',
    'type',
    'left',
    'right',
    'region_start',
    'region_end',
    'region_length',
    'long_reads_spanning',
    'long_reads_clipped',
    'short_depth_left',
    'short_depth_right',
    'long_depth_left',
    'long_depth_right',
    'call',
]


def is_circular(contig: dict) -> bool:
    return bool(CIRCULAR_RE.search(contig.get('description', '') or ''))


# --- composition and joins (from the references' placements) ---

def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def subtract(intervals, taken):
    """intervals minus taken (both merged, sorted)."""
    out = []
    for s, e in intervals:
        pos = s
        for ts, te in taken:
            if te <= pos or ts >= e:
                continue
            if ts > pos:
                out.append((pos, ts))
            pos = max(pos, te)
        if pos < e:
            out.append((pos, e))
    return out


def composition(contigs: list[dict], per_ref: list[tuple[str, dict]], min_block: int = DEFAULT_MIN_BLOCK) -> dict:
    """
    Pick references one at a time by the reconstruction sequence each newly explains.

    :param contigs: the reconstruction's contigs
    :param per_ref: (ref_id, intervals) per reference, intervals from contig_intervals: the contig
                    intervals its placements cover ('placements'; those shorter than min_block are
                    ignored here), and those its placements and alternative placements cover ('all')
    :return: dict with 'rows' (COMPOSITION_FIELDNAMES, plus a final 'none' row), 'segments'
             ({contig_id: [(start, end, ref_id or None)]}, tiling each contig) and 'blocks'
             ({ref_id: {contig_id: merged intervals of everything aligning to it}}, for the extent of a
             join: a repeat shared by two references, e.g. an IS at a join, is an alternative placement
             of one of them)
    """
    lengths = {c['seq_id']: len(c['seq']) for c in contigs}
    total = sum(lengths.values()) or 1
    labelling, blocks = {}, {}
    for ref_id, intervals in per_ref:
        labelling.setdefault(ref_id, {})
        blocks.setdefault(ref_id, {})
        for cid, ivs in intervals['placements'].items():
            long_ones = [(s, e) for s, e in ivs if e - s >= min_block]
            labelling[ref_id][cid] = merge(labelling[ref_id].get(cid, []) + long_ones)
        for cid, ivs in intervals['all'].items():
            blocks[ref_id][cid] = merge(blocks[ref_id].get(cid, []) + list(ivs))
    claimed = {cid: [] for cid in lengths}
    labelled = {cid: [] for cid in lengths}
    rows, cumulative, remaining = [], 0, dict(labelling)
    while remaining:
        best_id, best_new, best_parts = None, 0, None
        for ref_id, by_contig in remaining.items():
            parts = {cid: subtract(ivs, claimed.get(cid, [])) for cid, ivs in by_contig.items() if cid in lengths}
            new = sum(e - s for ivs in parts.values() for s, e in ivs)
            if new > best_new:
                best_id, best_new, best_parts = ref_id, new, parts
        if best_id is None or best_new < min_block:
            break
        del remaining[best_id]
        for cid, ivs in best_parts.items():
            claimed[cid] = merge(claimed[cid] + ivs)
            labelled[cid] += [(s, e, best_id) for s, e in ivs]
        cumulative += best_new
        rows.append({'order': len(rows) + 1, 'ref_id': best_id, 'new_bases': best_new,
                     'new_fraction': round(best_new / total, 4), 'cumulative_fraction': round(cumulative / total, 4)})
    unexplained = total - cumulative
    rows.append({'order': 'NA', 'ref_id': 'none', 'new_bases': unexplained,
                 'new_fraction': round(unexplained / total, 4), 'cumulative_fraction': 1.0 if total else 0.0})
    segments = {}
    for cid, n in lengths.items():
        segs = sorted(labelled[cid]) + [(s, e, None) for s, e in subtract([(0, n)], claimed[cid])]
        segs.sort()
        merged = []
        for s, e, lab in segs:
            if merged and merged[-1][2] == lab and merged[-1][1] == s:
                merged[-1] = (merged[-1][0], e, lab)
            else:
                merged.append((s, e, lab))
        segments[cid] = merged
    return {'rows': rows, 'segments': segments, 'blocks': blocks}


def _block_containing(blocks, pos):
    for s, e in blocks:
        if s <= pos < e or (pos == e and s < e):
            return s, e
    return None


def find_joins(contigs: list[dict], comp: dict, min_block: int = DEFAULT_MIN_BLOCK) -> list[dict]:
    """
    Where each contig's label changes between stretches of at least min_block bp, and each circular
    contig's own junction.

    :return: list of dicts: contig_id, type ('reference_switch', 'unexplained_boundary' or
             'circular_junction'), left, right (ref IDs or 'none'), region_start, region_end
    """
    joins = []
    for c in contigs:
        cid, n = c['seq_id'], len(c['seq'])
        major = [seg for seg in comp['segments'].get(cid, []) if seg[1] - seg[0] >= min_block]
        for (ls, le, la), (rs, re_, ra) in zip(major, major[1:]):
            if la == ra:
                continue
            start, end = le, rs
            if la is not None and ra is not None:
                a_block = _block_containing(comp['blocks'][la].get(cid, []), le - 1)
                b_block = _block_containing(comp['blocks'][ra].get(cid, []), rs)
                if a_block and b_block and b_block[0] < a_block[1]:
                    start, end = max(b_block[0], ls), min(a_block[1], re_)   # the overlap, e.g. a shared IS
            joins.append({'contig_id': cid,
                          'type': 'reference_switch' if la is not None and ra is not None else 'unexplained_boundary',
                          'left': la or 'none', 'right': ra or 'none',
                          'region_start': min(start, end), 'region_end': max(start, end)})
        if is_circular(c) and n:
            last = major[-1][2] if major else None
            first = major[0][2] if major else None
            joins.append({'contig_id': cid, 'type': 'circular_junction', 'left': last or 'none',
                          'right': first or 'none', 'region_start': n, 'region_end': n})
    return joins


def repeat_intervals(contigs: list[dict], tmpdir, threads: int = 1) -> dict:
    """
    Repeats within the reconstruction: {contig_id: [(start, end)]}, every copy of any sequence that occurs
    more than once in it (minimap2 self-alignment; at least MIN_REPEAT bp; the trivial self-hit excluded).
    """
    tmpdir = Path(tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)
    fasta = tmpdir / 'reconstruction_self.fa'
    alignment.write_fasta(fasta, [(f'c{i}', c['seq']) for i, c in enumerate(contigs)])
    # -D: skip only each sequence's trivial alignment to itself; repeats within a contig are kept
    cmd = ['minimap2', '-c', '-x', 'asm20', '-DP', '-t', str(threads), str(fasta), str(fasta)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError as e:
        raise alignment.AlignmentError("minimap2 not found on PATH") from e
    except subprocess.CalledProcessError as e:
        raise alignment.AlignmentError(f"minimap2 failed aligning the reconstruction to itself:\n{e.stderr}")
    out = {}
    for line in result.stdout.splitlines():
        f = line.split('\t')
        q, qs, qe, t, ts, te = int(f[0][1:]), int(f[2]), int(f[3]), int(f[5][1:]), int(f[7]), int(f[8])
        if q == t and qs == ts and qe == te:
            continue
        if qe - qs < MIN_REPEAT or te - ts < MIN_REPEAT:
            continue
        out.setdefault(contigs[q]['seq_id'], []).append((qs, qe))
        out.setdefault(contigs[t]['seq_id'], []).append((ts, te))
    return {cid: merge(v) for cid, v in out.items()}


def widen_to_repeats(found: list[dict], contigs: list[dict], repeats: dict) -> None:
    """
    Widen each join (in place) to cover any repeat touching it. A circular junction's region is in extended
    coordinates: from a repeat ending the contig to one starting it (region_end - length, past the end).
    """
    lengths = {c['seq_id']: len(c['seq']) for c in contigs}
    for j in found:
        n = lengths[j['contig_id']]
        reps = repeats.get(j['contig_id'], [])
        if j['type'] == 'circular_junction':
            a, b = n, n
            for s, e in reps:
                if e >= n - REPEAT_TOLERANCE:
                    a = min(a, s)
                if s <= REPEAT_TOLERANCE:
                    b = max(b, n + e)
        else:
            a, b = j['region_start'], j['region_end']
            changed = True
            while changed:
                changed = False
                for s, e in reps:
                    if s <= b + REPEAT_TOLERANCE and e >= a - REPEAT_TOLERANCE and (s < a or e > b):
                        a, b, changed = min(a, s), max(b, e), True
        j['region_start'], j['region_end'] = a, b


# --- reads aligned to the reconstruction ---

def extended(contig: dict) -> str:
    """A circular contig followed by up to CIRCULAR_EXTENSION bp of its own start (reads can cross the junction)."""
    seq = contig['seq']
    return seq + seq[:min(len(seq), CIRCULAR_EXTENSION)] if is_circular(contig) else seq


def align_reads_to_reconstruction(read_set: dict, contigs: list[dict], tmpdir, threads: int = 1) -> dict:
    """
    One read set aligned to the reconstruction's contigs (primary and supplementary alignments at least
    read_set['min_identity'] identical; PAF).

    :return: dict with 'depth' (per contig, depth along it, the extension folded back onto the start)
             and, for long reads, 'records': (contig index, start, end, clipped bp before start, clipped
             bp after end), in extended coordinates
    """
    tmpdir = Path(tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)
    fasta = tmpdir / f"reconstruction_for_{read_set['kind']}_reads.fa"
    seqs = [extended(c) for c in contigs]
    alignment.write_fasta(fasta, [(f'c{i}', s) for i, s in enumerate(seqs)])
    diffs = [[0] * (len(s) + 1) for s in seqs]
    records = [] if read_set['kind'] == 'long' else None
    cmd = ['minimap2', '-c', '-x', read_set['preset'], '--secondary=no', '-t', str(threads), str(fasta),
           *read_set['paths']]
    logger.info(f"Aligning {read_set['kind']} reads vs. the reconstruction ({len(contigs)} contig(s)) for joins")
    stderr_path = tmpdir / f"minimap2_{read_set['kind']}_reads_vs_reconstruction.log"
    try:
        with open(stderr_path, 'w') as stderr:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stderr, text=True)
    except FileNotFoundError as e:
        raise alignment.AlignmentError("minimap2 not found on PATH") from e
    for line in proc.stdout:
        f = line.split('\t', 12)
        matches, block = int(f[9]), int(f[10])
        if not block or matches / block < read_set['min_identity'] or '\ttp:A:S' in '\t' + f[12]:
            continue
        i = int(f[5][1:])
        qlen, qs, qe, ts, te = int(f[1]), int(f[2]), int(f[3]), int(f[7]), int(f[8])
        diffs[i][ts] += 1
        diffs[i][te] -= 1
        if records is not None:
            before, after = (qs, qlen - qe) if f[4] == '+' else (qlen - qe, qs)
            records.append((i, ts, te, before, after))
    if proc.wait() != 0:
        raise alignment.AlignmentError(f"minimap2 failed (exit code {proc.returncode}) aligning {read_set['kind']} "
                                       f"reads to the reconstruction:\n{stderr_path.read_text()[-2000:]}")
    depths = []
    for c, diff in zip(contigs, diffs):
        n = len(c['seq'])
        full = list(accumulate(diff))[:len(diff) - 1]
        d = full[:n]
        for p in range(n, len(full)):
            d[p - n] += full[p]
        depths.append(array('I', d))
    return {'depth': depths, 'records': records}


def _encode(result: dict) -> dict:
    return {'depth': [reads.run_length_encode(d) for d in result['depth']], 'records': result['records']}


def _decode(data: dict, contigs: list[dict]) -> dict:
    depths = []
    for c, segs in zip(contigs, data['depth']):
        d = array('I', bytes(4 * len(c['seq'])))
        for s, e, v in segs:
            if v:
                d[s:e] = array('I', [v]) * (e - s)
        depths.append(d)
    records = [tuple(r) for r in data['records']] if data['records'] is not None else None
    return {'depth': depths, 'records': records}


def reconstruction_reads(read_sets: list[dict], contigs: list[dict], tmpdir, threads: int = 1, cache_dir=None) -> dict:
    """
    Each read set aligned to the reconstruction (align_reads_to_reconstruction), cached in cache_dir by the
    reads and the reconstruction's sequences.

    :return: {kind ('short'/'long'): result}
    """
    out = {}
    target = 'reconstruction:' + reads.contigs_key(contigs) + ':' + ','.join(
        str(int(is_circular(c))) for c in contigs)
    for read_set in read_sets:
        path = (Path(cache_dir) / f"{reads.cache_key(read_set, target=target)}.reconstruction.json.gz"
                if cache_dir else None)
        result = None
        if path is not None and path.exists():
            try:
                with gzip.open(path, 'rt') as f:
                    data = json.load(f)
                if data.get('version') == CACHE_FORMAT_VERSION and data.get('ids') == [c['seq_id'] for c in contigs]:
                    result = _decode(data, contigs)
                    logger.info(f"Reads vs. the reconstruction ({read_set['kind']}) from the cache ({cache_dir})")
            except (OSError, ValueError, KeyError, EOFError):
                result = None
        if result is None:
            missing = [p for p in read_set['paths'] if not os.path.exists(p)]
            if missing:
                raise alignment.AlignmentError(f"Reads file(s) not found: {', '.join(missing)}")
            result = align_reads_to_reconstruction(read_set, contigs, tmpdir, threads)
            if path is not None:
                path.parent.mkdir(parents=True, exist_ok=True)
                fd, tmp = tempfile.mkstemp(dir=path.parent, prefix='.tmp_', suffix='.json.gz')
                os.close(fd)
                with gzip.open(tmp, 'wt') as f:
                    json.dump({'version': CACHE_FORMAT_VERSION, 'ids': [c['seq_id'] for c in contigs],
                               **_encode(result)}, f)
                os.replace(tmp, path)
        out[read_set['kind']] = result
    return out


# --- evidence at each join ---

def _flank_depth(depth, n, start, end, circular):
    positions = []
    for p in range(start, end):
        if circular:
            positions.append(depth[p % n])
        elif 0 <= p < n:
            positions.append(depth[p])
    return round(median(positions), 1) if positions else 'NA'


def join_evidence(join: dict, contig: dict, index: int, evidence: dict, anchor: int = DEFAULT_ANCHOR) -> dict:
    """Read evidence at one join (JOINS_FIELDNAMES, less the identifiers)."""
    n = len(contig['seq'])
    circular = is_circular(contig)
    a, b = join['region_start'], join['region_end']
    row = {}
    for kind in ('short', 'long'):
        result = evidence.get(kind)
        if result is None:
            row[f'{kind}_depth_left'] = row[f'{kind}_depth_right'] = 'NA'
            continue
        d = result['depth'][index]
        row[f'{kind}_depth_left'] = _flank_depth(d, n, a - FLANK, a, circular)
        row[f'{kind}_depth_right'] = _flank_depth(d, n, b, b + FLANK, circular)
    long = evidence.get('long')
    if long is None:
        row.update({'long_reads_spanning': 'NA', 'long_reads_clipped': 'NA', 'call': 'no_long_reads'})
        return row
    # The join's positions in extended coordinates: as is, and for a circular contig shifted by its length
    # (the copy of its start appended to it)
    ext = len(extended(contig))
    shifts = [0] + ([n] if circular and b + n <= ext else [])
    spanning, clipped = set(), set()
    for k, (i, ts, te, before, after) in enumerate(long['records']):
        if i != index:
            continue
        for shift in shifts:
            lo, hi = a + shift, b + shift
            if ts <= lo - anchor and te >= hi + anchor:
                spanning.add(k)
            if (lo - CLIP_WINDOW <= te <= hi + CLIP_WINDOW and after >= CLIP_MIN) or \
                    (lo - CLIP_WINDOW <= ts <= hi + CLIP_WINDOW and before >= CLIP_MIN):
                clipped.add(k)
    clipped -= spanning
    if len(spanning) >= MIN_SPANNING:
        call = 'supported'
    elif not spanning and len(clipped) >= MIN_CLIPPED:
        call = 'contradicted'
    else:
        call = 'weak'
    row.update({'long_reads_spanning': len(spanning), 'long_reads_clipped': len(clipped), 'call': call})
    return row


def assess_joins(contigs: list[dict], joins: list[dict], evidence: dict = None, anchor: int = DEFAULT_ANCHOR) -> list[dict]:
    """joins.tsv rows: each join with its read evidence (calls 'no_long_reads' without long reads)."""
    index = {c['seq_id']: i for i, c in enumerate(contigs)}
    by_id = {c['seq_id']: c for c in contigs}
    rows, counter = [], {}
    for j in joins:
        cid = j['contig_id']
        counter[cid] = counter.get(cid, 0) + 1
        row = {**j, 'join': f"join_{counter[cid]}", 'region_length': j['region_end'] - j['region_start']}
        if evidence:
            row.update(join_evidence(j, by_id[cid], index[cid], evidence, anchor))
        else:
            row.update({k: 'NA' for k in JOINS_FIELDNAMES if k not in row})
            row['call'] = 'no_long_reads'
        n = len(by_id[cid]['seq'])
        if row['region_end'] > n:      # a circular junction's region past the end: wraps to the start
            row['region_end'] -= n
        rows.append(row)
    return rows


def contig_intervals(placements: list[dict], alternatives: list[dict] = ()) -> dict:
    """
    The parts of each contig a reference's alignments cover: {'placements': {contig_id: [(start, end)]},
    'all': the same counting alternative placements too}.
    """
    out = {'placements': {}, 'all': {}}
    for group, alns in (('placements', placements), ('alternatives', alternatives)):
        for aln in alns:
            if aln.get('is_mapped') and aln.get('cigar_parsed'):
                s, e = assessment.query_interval_on_contig(aln)
                if group == 'placements':
                    out['placements'].setdefault(aln['contig_id'], []).append((s, e))
                out['all'].setdefault(aln['contig_id'], []).append((s, e))
    return out


def write_rows(path, rows: list[dict], fieldnames: list) -> None:
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def write_depths(path, contigs: list[dict], depths: list) -> None:
    """Depth along each reconstruction contig, as gzipped bedGraph."""
    with gzip.open(path, 'wt') as f:
        for c, d in zip(contigs, depths):
            for s, e, v in reads.run_length_encode(d):
                f.write(f"{c['seq_id']}\t{s}\t{e}\t{v}\n")
