"""
Reads as evidence: is the reference's sequence in the sample at all?

The sample's reads (short and/or long) are aligned to the references with minimap2 and turned into
per-base read depth along each reference. That separates what a reconstruction lacks into:

- not in the reads: the sample doesn't contain that sequence (it differs from the reference), and
- in the reads: the sample contains it, but the assembly (or, with --assembly, the binning) lost it.

Every alignment counts: primary, secondary (minimap2's other equally good placements, e.g. of a
repeat) and supplementary (e.g. a long read split across the reference's origin), as long as it is
at least ``min_identity`` identical (1 - NM / alignment columns). Reads are aligned to the references
only, so reads from the chromosome or another plasmid that share sequence with the reference (an IS
element, a transposon) count too: that sequence is in the sample, just not necessarily on this
plasmid. Repeats are flagged by the gap types ('repeat' vs. 'unique') rather than here.

Depth along each reference can be cached (``cache_dir``), keyed by the reads files (path and size),
the minimap2 preset, the identity threshold and the reference sequence, so assessing many bins of one
sample against the same references aligns the reads once.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import logging
import os
import re
import subprocess
import tempfile

from array import array
from itertools import accumulate
from pathlib import Path

import plasmid_surgery.alignment as alignment
import plasmid_surgery.assessment as assessment

logger = logging.getLogger(__name__)

SHORT_READ_PRESET = 'sr'
LONG_READ_PRESETS = ['map-ont', 'lr:hq', 'map-hifi', 'map-pb']
DEFAULT_LONG_READ_PRESET = 'map-ont'
DEFAULT_MIN_IDENTITY = {'short': 0.95, 'long': 0.85}
DEFAULT_MIN_DEPTH = 3
# Deletions within a read alignment shorter than this still count as covered (sequencing errors,
# e.g. nanopore homopolymer deletions, shouldn't punch holes in the depth)
MAX_COVERED_DELETION = 50
# A gap is 'in_reads' if at least this fraction of it has read depth >= min_depth, 'not_in_reads' if at
# most NOT_IN_READS_FRACTION does, and 'partly_in_reads' otherwise
IN_READS_FRACTION = 0.9
NOT_IN_READS_FRACTION = 0.1

CACHE_FORMAT_VERSION = 1

# Added to the assessment metrics ('NA' when no reads are given)
READS_METRICS_FIELDNAMES = [
    'reads_ref_fraction_covered',
    'reads_median_depth',
    'unique_gap_in_reads_length',
    'assembly_unique_gap_in_reads_length',
]

GAP_SUPPORT_FIELDNAMES = [
    'level',
    'gap',
    'type',
    'ref_start',
    'ref_end',
    'length',
    'reads_fraction_covered',
    'reads_median_depth',
    'call',
]

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def read_sets(short_reads=None, long_reads=None, long_read_preset=DEFAULT_LONG_READ_PRESET,
              min_identity=None) -> list[dict]:
    """
    The read sets to align: one per read type given.

    :param short_reads: 1 or 2 FASTQ files (single-end or interleaved; or R1 and R2)
    :param long_reads: FASTQ file(s)
    :param min_identity: minimum alignment identity for both types (default: DEFAULT_MIN_IDENTITY per type)
    """
    sets = []
    if short_reads:
        sets.append({'kind': 'short', 'paths': [str(p) for p in short_reads], 'preset': SHORT_READ_PRESET,
                     'min_identity': min_identity if min_identity is not None else DEFAULT_MIN_IDENTITY['short']})
    if long_reads:
        sets.append({'kind': 'long', 'paths': [str(p) for p in long_reads], 'preset': long_read_preset,
                     'min_identity': min_identity if min_identity is not None else DEFAULT_MIN_IDENTITY['long']})
    return sets


def cache_key(read_set: dict, ref_seq: str) -> str:
    """Cache key for one read set's depth along one reference: reads files by real path and size."""
    files = []
    for p in read_set['paths']:
        try:
            files.append([os.path.realpath(p), os.path.getsize(p)])
        except OSError:
            files.append([os.path.realpath(p), None])
    key = {
        'version': CACHE_FORMAT_VERSION,
        'kind': read_set['kind'],
        'preset': read_set['preset'],
        'min_identity': read_set['min_identity'],
        'files': files,
        'ref': hashlib.sha1(ref_seq.upper().encode()).hexdigest(),
    }
    return hashlib.sha1(json.dumps(key, sort_keys=True).encode()).hexdigest()


def run_length_encode(depth) -> list[tuple[int, int, int]]:
    """(start, end, depth) segments tiling the depth array, adjacent equal depths merged."""
    segments = []
    start = 0
    for i in range(1, len(depth) + 1):
        if i == len(depth) or depth[i] != depth[start]:
            segments.append((start, i, depth[start]))
            start = i
    return segments


def write_cached_depth(path: Path, depth) -> None:
    """Write a depth array (run-length encoded, gzipped), atomically: concurrent runs may share a cache."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix='.tmp_', suffix='.gz')
    os.close(fd)
    with gzip.open(tmp, 'wt') as f:
        f.write(f'{len(depth)}\n')
        for start, end, d in run_length_encode(depth):
            f.write(f'{start}\t{end}\t{d}\n')
    os.replace(tmp, path)


def read_cached_depth(path: Path, ref_length: int):
    """A cached depth array, or None if missing or not for a reference of this length."""
    try:
        with gzip.open(path, 'rt') as f:
            if int(f.readline()) != ref_length:
                return None
            depth = array('I', bytes(4 * ref_length))
            for line in f:
                start, end, d = map(int, line.split('\t'))
                if d:
                    depth[start:end] = array('I', [d]) * (end - start)
            return depth
    except (OSError, ValueError, EOFError):
        return None


def alignment_blocks(pos: int, cigar: str) -> tuple[list[tuple[int, int]], int]:
    """
    Reference intervals covered by an alignment (0-based start ``pos``), and its number of alignment
    columns (M/=/X/I/D, for identity). Deletions shorter than MAX_COVERED_DELETION are covered.
    """
    blocks = []
    columns = 0
    ref = pos
    for n, op in CIGAR_RE.findall(cigar):
        n = int(n)
        if op in 'M=X':
            blocks.append((ref, ref + n))
            ref += n
            columns += n
        elif op == 'D':
            if n < MAX_COVERED_DELETION:
                blocks.append((ref, ref + n))
            ref += n
            columns += n
        elif op == 'N':
            ref += n
        elif op == 'I':
            columns += n
    merged = []
    for s, e in blocks:
        if merged and merged[-1][1] == s:
            merged[-1] = (merged[-1][0], e)
        else:
            merged.append((s, e))
    return merged, columns


def align_read_set(read_set: dict, references: list[dict], tmpdir, threads: int = 1) -> list:
    """
    Align one read set to the references (one minimap2 call) and return the read depth along each.

    :return: list of depth arrays (array('I'), one per reference, in order)
    :raises alignment.AlignmentError: if minimap2 is missing or fails
    """
    tmpdir = Path(tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)
    refs_fasta = tmpdir / f"refs_for_{read_set['kind']}_reads.fa"
    # References get positional names: their IDs need not be unique
    alignment.write_fasta(refs_fasta, [(f'ref{i}', r['seq']) for i, r in enumerate(references)])
    lengths = [len(r['seq']) for r in references]
    diffs = [[0] * (n + 1) for n in lengths]
    cmd = ['minimap2', '-a', '-x', read_set['preset'], '--secondary=yes', '--sam-hit-only',
           '-t', str(threads), str(refs_fasta), *read_set['paths']]
    logger.info(f"Aligning {read_set['kind']} reads ({', '.join(read_set['paths'])}) vs. "
                f"{len(references)} reference(s)")
    counted = filtered = 0
    stderr_path = tmpdir / f"minimap2_{read_set['kind']}_reads.log"
    try:
        with open(stderr_path, 'w') as stderr:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stderr, text=True)
    except FileNotFoundError as e:
        raise alignment.AlignmentError("minimap2 not found on PATH") from e
    min_identity = read_set['min_identity']
    for line in proc.stdout:
        if line.startswith('@'):
            continue
        fields = line.rstrip('\n').split('\t')
        if int(fields[1]) & 4 or fields[5] == '*':
            continue
        i = int(fields[2][3:])
        blocks, columns = alignment_blocks(int(fields[3]) - 1, fields[5])
        nm = next((int(f[5:]) for f in fields[11:] if f.startswith('NM:i:')), None)
        if columns == 0 or nm is None or 1 - nm / columns < min_identity:
            filtered += 1
            continue
        counted += 1
        diff, n = diffs[i], lengths[i]
        for s, e in blocks:
            diff[max(s, 0)] += 1
            diff[min(e, n)] -= 1
    if proc.wait() != 0:
        raise alignment.AlignmentError(f"minimap2 failed (exit code {proc.returncode}) aligning "
                                       f"{read_set['kind']} reads:\n{stderr_path.read_text()[-2000:]}")
    logger.info(f"{counted} {read_set['kind']}-read alignment(s) counted, {filtered} below "
                f"{min_identity:.0%} identity ignored")
    return [array('I', list(accumulate(diff))[:n]) for diff, n in zip(diffs, lengths)]


def reads_depth(read_sets_: list[dict], references: list[dict], tmpdir, threads: int = 1,
                cache_dir=None) -> list:
    """
    Read depth along each reference, summed over the read sets; with ``cache_dir``, depth already
    cached for a (read set, reference) is reused and only the rest are aligned.

    :return: list of depth arrays (array('I')), one per reference, in order
    """
    totals = [array('I', bytes(4 * len(r['seq']))) for r in references]
    for read_set in read_sets_:
        depths = [None] * len(references)
        paths = [None] * len(references)
        if cache_dir:
            for i, r in enumerate(references):
                paths[i] = Path(cache_dir) / f"{cache_key(read_set, r['seq'])}.depth.gz"
                depths[i] = read_cached_depth(paths[i], len(r['seq']))
            cached = sum(d is not None for d in depths)
            if cached:
                logger.info(f"Read depth for {cached} of {len(references)} reference(s) from the cache ({cache_dir})")
        todo = [i for i, d in enumerate(depths) if d is None]
        if todo:
            missing = [p for p in read_set['paths'] if not os.path.exists(p)]
            if missing:
                raise alignment.AlignmentError(f"Reads file(s) not found: {', '.join(missing)}")
            aligned = align_read_set(read_set, [references[i] for i in todo], tmpdir, threads)
            for i, d in zip(todo, aligned):
                depths[i] = d
                if cache_dir:
                    write_cached_depth(paths[i], d)
        for total, d in zip(totals, depths):
            for j, v in enumerate(d):
                if v:
                    total[j] += v
    return totals


def _positions(start: int, end: int, ref_length: int):
    """Positions of an interval that may wrap the origin (end > ref_length)."""
    for s, e in assessment.split_circular_interval(start, end, ref_length):
        yield from range(s, e)


def fraction_covered(depth, start: int, end: int, min_depth: int) -> float:
    n = len(depth)
    positions = list(_positions(start, end, n))
    if not positions:
        return 0.0
    return round(sum(1 for p in positions if depth[p] >= min_depth) / len(positions), 4)


def median_depth(depth, start: int = 0, end: int = None):
    n = len(depth)
    values = sorted(depth[p] for p in _positions(start, n if end is None else end, n))
    if not values:
        return 0
    mid = len(values) // 2
    return values[mid] if len(values) % 2 else (values[mid - 1] + values[mid]) / 2


def gap_call(fraction: float) -> str:
    if fraction >= IN_READS_FRACTION:
        return 'in_reads'
    if fraction <= NOT_IN_READS_FRACTION:
        return 'not_in_reads'
    return 'partly_in_reads'


def gap_support(gaps: list[dict], depth, level: str, min_depth: int = DEFAULT_MIN_DEPTH) -> list[dict]:
    """
    Read support for each gap (as from assessment.find_typed_gaps).

    :param level: 'reconstruction' (named as in gaps.bed: gap_<n>_<type>) or 'assembly'
                  (assembly_gap_<n>_<type>: what the whole assembly lacks)
    """
    n = len(depth)
    rows = []
    for i, g in enumerate(gaps, start=1):
        fraction = fraction_covered(depth, g['start'], g['end'], min_depth)
        rows.append({
            'level': level,
            'gap': f"{'gap' if level == 'reconstruction' else 'assembly_gap'}_{i}_{g['type']}",
            'type': g['type'],
            'ref_start': g['start'],
            'ref_end': g['end'] - n if g['end'] > n else g['end'],
            'length': g['end'] - g['start'],
            'reads_fraction_covered': fraction,
            'reads_median_depth': median_depth(depth, g['start'], g['end']),
            'call': gap_call(fraction),
        })
    return rows


def _gap_bases_in_reads(gaps: list[dict], depth, min_depth: int) -> int:
    n = len(depth)
    return sum(1 for g in gaps if g['type'] == 'unique'
               for p in _positions(g['start'], g['end'], n) if depth[p] >= min_depth)


def reads_metrics(depth, gaps: list[dict], assembly_gaps: list[dict] = None,
                  min_depth: int = DEFAULT_MIN_DEPTH) -> dict:
    """READS_METRICS_FIELDNAMES for one reference."""
    n = len(depth)
    return {
        'reads_ref_fraction_covered': fraction_covered(depth, 0, n, min_depth) if n else 0.0,
        'reads_median_depth': median_depth(depth),
        'unique_gap_in_reads_length': _gap_bases_in_reads(gaps, depth, min_depth),
        'assembly_unique_gap_in_reads_length': (_gap_bases_in_reads(assembly_gaps, depth, min_depth)
                                                if assembly_gaps is not None else 'NA'),
    }


def write_bedgraph(path, ref_id: str, depth) -> None:
    """Read depth along the reference as bedGraph (zero-depth runs included)."""
    with open(path, 'w') as f:
        for start, end, d in run_length_encode(depth):
            f.write(f'{ref_id}\t{start}\t{end}\t{d}\n')


def write_gap_support(path, rows: list[dict]) -> None:
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=GAP_SUPPORT_FIELDNAMES, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
