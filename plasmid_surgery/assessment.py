"""
Assess how well a plasmid reconstruction matches a reference plasmid.

All functions take aligned segments as produced by
alignment.align_contig_vs_circular_ref: ref coordinates are relative to the
original (circular) reference, ref_start is in [0, ref_length) and ref_end may
exceed ref_length when an alignment wraps around the origin.
"""
from __future__ import annotations

import csv
import logging

from collections import defaultdict
from textwrap import wrap

logger = logging.getLogger(__name__)

QUERY_ALIGNED_OPS = ('alignment_match', 'sequence_match', 'sequence_mismatch', 'insertion_to_ref')
REF_ALIGNED_OPS = ('alignment_match', 'sequence_match', 'sequence_mismatch')
REF_ONLY_OPS = ('deletion_to_ref', 'skip_ref')
CLIP_OPS = ('soft_clip', 'hard_clip')

METRICS_FIELDNAMES = [
    'sample_id',
    'ref_id',
    'ref_length',
    'num_contigs',
    'total_contig_length',
    'ref_bases_covered',
    'ref_fraction_covered',
    'ref_bases_covered_multiple',
    'ref_fraction_covered_including_repeats',
    'num_gaps',
    'total_gap_length',
    'largest_gap_length',
    'num_repeat_gaps',
    'repeat_gap_length',
    'num_unique_gaps',
    'unique_gap_length',
    'largest_unique_gap_length',
    'contig_bases_aligned',
    'contig_fraction_aligned',
    'num_unaligned_contigs',
    'num_unaligned_regions',
    'unaligned_region_bases',
    'num_large_insertions',
    'num_large_deletions',
    'num_ambiguous_alignments',
    'num_ambiguous_contigs',
]

EVENTS_FIELDNAMES = [
    'contig_id',
    'event_type',
    'ref_start',
    'ref_end',
    'length',
    'contig_position',
]


def split_circular_interval(start: int, end: int, ref_length: int) -> list[tuple[int, int]]:
    """
    Split an interval that may wrap the origin (end > ref_length) into
    non-wrapping intervals within [0, ref_length].
    """
    if end <= ref_length:
        return [(start, end)]
    return [(start, ref_length), (0, end - ref_length)]


def coverage_segments(alignments: list[dict], ref_length: int) -> list[tuple[int, int, int]]:
    """
    Compute reference coverage depth.

    :return: Contiguous (start, end, depth) segments tiling [0, ref_length),
             with adjacent segments of equal depth merged.
    """
    deltas = defaultdict(int)
    deltas[0] += 0
    deltas[ref_length] += 0
    for aln in alignments:
        for start, end in split_circular_interval(aln['ref_start'], aln['ref_end'], ref_length):
            if end > start:
                deltas[start] += 1
                deltas[end] -= 1

    positions = sorted(deltas)
    segments = []
    depth = 0
    for pos, next_pos in zip(positions, positions[1:]):
        depth += deltas[pos]
        if next_pos <= pos:
            continue
        if segments and segments[-1][2] == depth and segments[-1][1] == pos:
            segments[-1] = (segments[-1][0], next_pos, depth)
        else:
            segments.append((pos, next_pos, depth))

    return segments


def find_typed_gaps(alignments: list[dict], alternative_placements: list[dict], ref_length: int) -> list[dict]:
    """
    Find uncovered regions of the reference, typed by whether an alternative
    placement of some contig covers them.

    'repeat' gaps are covered by an alternative placement (e.g. another copy of
    a collapsed repeat): the sequence is likely present in the reconstruction,
    just not placed there. 'unique' gaps are not covered by anything: the
    sequence is likely genuinely missing from the reconstruction. A gap that is
    partly covered is split into repeat and unique parts.

    :return: List of dicts with keys 'start', 'end' (> ref_length if the gap
             spans the origin) and 'type', ordered by start
    """
    alternative_segments = coverage_segments(alternative_placements, ref_length)
    typed = []
    for gap_start, gap_end, depth in coverage_segments(alignments, ref_length):
        if depth != 0:
            continue
        for seg_start, seg_end, alt_depth in alternative_segments:
            start = max(gap_start, seg_start)
            end = min(gap_end, seg_end)
            if end <= start:
                continue
            gap_type = 'repeat' if alt_depth > 0 else 'unique'
            if typed and typed[-1]['type'] == gap_type and typed[-1]['end'] == start:
                typed[-1]['end'] = end
            else:
                typed.append({'start': start, 'end': end, 'type': gap_type})

    if (len(typed) >= 2 and typed[0]['start'] == 0 and typed[-1]['end'] == ref_length
            and typed[0]['type'] == typed[-1]['type']):
        first = typed.pop(0)
        typed[-1]['end'] = ref_length + first['end']

    return typed


def leading_clip_length(aln: dict) -> int:
    """
    Number of query bases clipped (soft or hard) before the alignment starts.
    """
    clipped = 0
    for op, length in aln['cigar_parsed']:
        if op not in CLIP_OPS:
            break
        clipped += length
    return clipped


def to_contig_coords(aln: dict, start: int, end: int) -> tuple[int, int]:
    """
    Convert a query interval in alignment orientation into the original
    orientation of the contig.
    """
    if aln.get('is_reverse'):
        contig_length = aln['query_full_length']
        return contig_length - end, contig_length - start
    return start, end


def query_interval_on_contig(aln: dict) -> tuple[int, int]:
    """
    The interval of the contig (in its original orientation, 0-based half-open)
    covered by an alignment. Accounts for hard clipping on supplementary
    alignments, which query_start/query_end from pysam do not.
    """
    start = leading_clip_length(aln)
    aligned = sum(length for op, length in aln['cigar_parsed'] if op in QUERY_ALIGNED_OPS)
    return to_contig_coords(aln, start, start + aligned)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def complement_intervals(intervals: list[tuple[int, int]], length: int) -> list[tuple[int, int]]:
    """
    Regions of [0, length) not covered by any of the (merged) intervals.
    """
    regions = []
    pos = 0
    for start, end in merge_intervals(intervals):
        if start > pos:
            regions.append((pos, start))
        pos = max(pos, end)
    if pos < length:
        regions.append((pos, length))
    return regions


def indel_events(aln: dict, ref_length: int, min_length: int) -> list[dict]:
    """
    Insertions and deletions of at least min_length within an alignment.

    Insertions are sequence present in the contig but not the reference
    (ref_start == ref_end, the reference position the insertion follows).
    Deletions are reference sequence missing from the contig. Reference
    positions are folded into [0, ref_length); contig positions are in the
    contig's original orientation.
    """
    events = []
    r = aln['ref_start']
    q = leading_clip_length(aln)
    for op, length in aln['cigar_parsed']:
        if op in CLIP_OPS or op == 'padding':
            continue
        if op in REF_ALIGNED_OPS:
            r += length
            q += length
        elif op == 'insertion_to_ref':
            if length >= min_length:
                contig_start, _ = to_contig_coords(aln, q, q + length)
                events.append({
                    'contig_id': aln['contig_id'],
                    'event_type': 'insertion',
                    'ref_start': r % ref_length,
                    'ref_end': r % ref_length,
                    'length': length,
                    'contig_position': contig_start,
                })
            q += length
        elif op in REF_ONLY_OPS:
            if length >= min_length:
                contig_pos, _ = to_contig_coords(aln, q, q)
                events.append({
                    'contig_id': aln['contig_id'],
                    'event_type': 'deletion',
                    'ref_start': r % ref_length,
                    'ref_end': (r + length - 1) % ref_length + 1,
                    'length': length,
                    'contig_position': contig_pos,
                })
            r += length

    return events


def assess(contigs: list[dict], alignments: list[dict], ref_length: int,
           min_indel_length: int = 50, min_unaligned_length: int = 100,
           alternative_placements: list[dict] = ()) -> dict:
    """
    Assess a plasmid reconstruction against a reference.

    :param contigs: All contigs in the reconstruction, as parsed by parsers.parse_fasta
    :param alignments: Placement alignments (primary + supplementary) that passed QC
    :param alternative_placements: Secondary alignments that are alternative
                                   placements of a contig (see
                                   alignment.select_alternative_placements)
    :param ref_length: Length of the reference
    :param min_indel_length: Minimum insertion/deletion length to report as an event
    :param min_unaligned_length: Minimum length of an unaligned contig region to report
    :return: dict with keys:
             'metrics': dict of summary metrics (METRICS_FIELDNAMES, minus sample/ref IDs),
             'gaps': list of typed gap dicts (see find_typed_gaps),
             'events': list of indel event dicts,
             'unaligned_regions': list of dicts with contig_id, start, end, seq,
             'contig_fraction_aligned': dict of contig_id -> fraction
    """
    alignments = [a for a in alignments if a.get('is_mapped')]

    segments = coverage_segments(alignments, ref_length)
    ref_bases_covered = sum(end - start for start, end, depth in segments if depth >= 1)
    ref_bases_covered_multiple = sum(end - start for start, end, depth in segments if depth >= 2)
    alternative_placements = [a for a in alternative_placements if a.get('is_mapped')]
    gaps = find_typed_gaps(alignments, alternative_placements, ref_length)
    gap_lengths = [g['end'] - g['start'] for g in gaps]
    repeat_gap_lengths = [g['end'] - g['start'] for g in gaps if g['type'] == 'repeat']
    unique_gap_lengths = [g['end'] - g['start'] for g in gaps if g['type'] == 'unique']

    intervals_by_contig = defaultdict(list)
    for aln in alignments:
        intervals_by_contig[aln['contig_id']].append(query_interval_on_contig(aln))

    contig_fraction_aligned = {}
    unaligned_regions = []
    contig_bases_aligned = 0
    num_unaligned_contigs = 0
    for contig in contigs:
        contig_id = contig['seq_id']
        contig_length = len(contig['seq'])
        merged = merge_intervals(intervals_by_contig.get(contig_id, []))
        aligned = sum(end - start for start, end in merged)
        contig_bases_aligned += aligned
        contig_fraction_aligned[contig_id] = round(aligned / contig_length, 4) if contig_length else 0.0
        if aligned == 0:
            num_unaligned_contigs += 1
        for start, end in complement_intervals(merged, contig_length):
            if end - start >= min_unaligned_length:
                unaligned_regions.append({
                    'contig_id': contig_id,
                    'start': start,
                    'end': end,
                    'seq': contig['seq'][start:end],
                })

    events = []
    for aln in alignments:
        events += indel_events(aln, ref_length, min_indel_length)
    events.sort(key=lambda e: (e['ref_start'], e['contig_id']))

    total_contig_length = sum(len(c['seq']) for c in contigs)
    metrics = {
        'ref_length': ref_length,
        'num_contigs': len(contigs),
        'total_contig_length': total_contig_length,
        'ref_bases_covered': ref_bases_covered,
        'ref_fraction_covered': round(ref_bases_covered / ref_length, 4) if ref_length else 0.0,
        'ref_bases_covered_multiple': ref_bases_covered_multiple,
        'ref_fraction_covered_including_repeats': round((ref_bases_covered + sum(repeat_gap_lengths)) / ref_length, 4) if ref_length else 0.0,
        'num_gaps': len(gaps),
        'total_gap_length': sum(gap_lengths),
        'largest_gap_length': max(gap_lengths, default=0),
        'num_repeat_gaps': len(repeat_gap_lengths),
        'repeat_gap_length': sum(repeat_gap_lengths),
        'num_unique_gaps': len(unique_gap_lengths),
        'unique_gap_length': sum(unique_gap_lengths),
        'largest_unique_gap_length': max(unique_gap_lengths, default=0),
        'contig_bases_aligned': contig_bases_aligned,
        'contig_fraction_aligned': round(contig_bases_aligned / total_contig_length, 4) if total_contig_length else 0.0,
        'num_unaligned_contigs': num_unaligned_contigs,
        'num_unaligned_regions': len(unaligned_regions),
        'unaligned_region_bases': sum(r['end'] - r['start'] for r in unaligned_regions),
        'num_large_insertions': sum(1 for e in events if e['event_type'] == 'insertion'),
        'num_large_deletions': sum(1 for e in events if e['event_type'] == 'deletion'),
        'num_ambiguous_alignments': sum(1 for a in alignments if a.get('mapq') == 0),
        'num_ambiguous_contigs': len({a['contig_id'] for a in alternative_placements}),
    }

    return {
        'metrics': metrics,
        'gaps': gaps,
        'events': events,
        'unaligned_regions': unaligned_regions,
        'contig_fraction_aligned': contig_fraction_aligned,
    }


def write_metrics(path, metrics: dict) -> None:
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=METRICS_FIELDNAMES, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerow(metrics)


def write_gaps_bed(path, gaps: list[dict], ref_id: str, ref_length: int) -> None:
    """
    Write gaps as BED4, named 'gap_<n>_<type>' (type: repeat or unique). A gap
    spanning the origin is written as two records with the same name.
    """
    with open(path, 'w') as f:
        for i, gap in enumerate(gaps, start=1):
            for s, e in split_circular_interval(gap['start'], gap['end'], ref_length):
                if e > s:
                    f.write(f"{ref_id}\t{s}\t{e}\tgap_{i}_{gap['type']}\n")


def write_events(path, events: list[dict]) -> None:
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=EVENTS_FIELDNAMES, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for event in events:
            writer.writerow(event)


def write_unaligned_fasta(path, unaligned_regions: list[dict], bases_per_line: int = 70) -> None:
    """
    Write unaligned contig regions. Headers are '<contig_id>:<start>-<end>'
    (1-based, inclusive).
    """
    with open(path, 'w') as f:
        for region in unaligned_regions:
            length = region['end'] - region['start']
            f.write(f">{region['contig_id']}:{region['start'] + 1}-{region['end']} length={length}\n")
            for line in wrap(region['seq'], bases_per_line):
                f.write(f"{line}\n")
