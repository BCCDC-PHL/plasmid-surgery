from __future__ import annotations

import os
import logging
import subprocess

from pathlib import Path

import plasmid_surgery.parsers as parsers

logger = logging.getLogger(__name__)


class AlignmentError(Exception):
    pass


def align_contig_vs_ref(contig: dict, ref: dict, tmpdir: Path) -> list[dict]:
    """
    Align contig vs ref using minimap2

    :param contig: Contig, as parsed by parsers.parse_fasta
    :type contig: dict
    :param ref: Reference, as parsed by parsers.parse_fasta
    :type ref: dict
    :param tmpdir: Directory for intermediate files
    :type tmpdir: Path
    :return: Aligned segments, as parsed by parsers.parse_sam
    :rtype: list[dict]
    :raises AlignmentError: if minimap2 is missing or fails
    """
    os.makedirs(tmpdir, exist_ok=True)

    contig_id = contig['seq_id']
    ref_id = ref['seq_id']
    contig_path = os.path.join(tmpdir, f"{contig_id}.fa")

    with open(contig_path, 'w') as f:
        f.write(f">{contig['seq_id']}\n")
        f.write(f"{contig['seq']}\n")

    # Always (re)write the ref, so a different ref with the same ID
    # left over in tmpdir from a previous run is never reused.
    ref_path = os.path.join(tmpdir, f"{ref_id}.fa")
    with open(ref_path, 'w') as f:
        f.write(f">{ref['seq_id']}\n")
        f.write(f"{ref['seq']}\n")

    alignment_path = os.path.join(tmpdir, f"{contig_id}_vs_{ref_id}.sam")

    minimap2_cmd = [
        'minimap2',
        '-a',
        '-cx',
        'asm10',
        '--secondary=no',
        '--eqx',
        '-o', alignment_path,
        ref_path,
        contig_path,
    ]

    try:
        subprocess.run(minimap2_cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError as e:
        raise AlignmentError("minimap2 not found on PATH") from e
    except subprocess.CalledProcessError as e:
        raise AlignmentError(
            f"minimap2 failed (exit code {e.returncode}) aligning {contig_id} vs. {ref_id}:\n{e.stderr}"
        ) from e

    aligned_segments = parsers.parse_sam(alignment_path)
    for aligned_segment in aligned_segments:
        aligned_segment['query_full_length'] = len(contig['seq'])
        aligned_segment['ref_full_length'] = len(ref['seq'])

    return aligned_segments


def alignment_qc(aln, min_identity=80.0, min_aln_length=50, max_gap_frac=0.8):
    """
    Determine if alignment passes filters.

    :param aln: Aligned segment, as parsed by parsers.parse_sam
    :return: dict with keys 'qc_stats' and 'overall_qc_pass'
    """
    aln_length = aln.get('query_alignment_length') or 0
    gap_bases = sum(length for op, length in aln.get('cigar_parsed', []) if op in ('insertion_to_ref', 'deletion_to_ref'))
    gap_frac = gap_bases / aln_length if aln_length > 0 else 1.0

    alignment_qc = {
        'qc_stats': {
            'min_identity': {
                'min_value': min_identity,
                'value': aln['percent_identity'],
            },
            'aln_length': {
                'min_value': min_aln_length,
                'value': aln_length,
            },
            'gap_frac': {
                'max_value': max_gap_frac,
                'value': gap_frac,
            },
        },
        'overall_qc_pass': bool(aln.get('is_mapped', True)),
    }

    for qc_stat, qc_values in alignment_qc['qc_stats'].items():
        actual_value = qc_values['value']
        if 'min_value' in qc_values:
            if actual_value < qc_values['min_value']:
                alignment_qc['overall_qc_pass'] = False
        elif 'max_value' in qc_values:
            if actual_value > qc_values['max_value']:
                alignment_qc['overall_qc_pass'] = False

    return alignment_qc


def query_offset_at_ref_pos(aln: dict, ref_pos: int) -> int:
    """
    Find the offset into an alignment's query_alignment_sequence that corresponds
    to a given (0-based) position on the reference.

    If ref_pos falls inside a deletion, the offset of the next aligned query base
    is returned. Insertions immediately before ref_pos are skipped (i.e. they are
    treated as belonging to the preceding reference region).

    :param aln: Aligned segment, as parsed by parsers.parse_sam
    :param ref_pos: Position on the reference, ref_start <= ref_pos <= ref_end
    :return: Offset into aln['query_alignment_sequence']
    """
    if ref_pos <= aln['ref_start']:
        return 0
    if ref_pos >= aln['ref_end']:
        return len(aln['query_alignment_sequence'])

    r = aln['ref_start']
    q = 0
    for op, length in aln['cigar_parsed']:
        if op in ('alignment_match', 'sequence_match', 'sequence_mismatch'):
            if r + length > ref_pos:
                return q + (ref_pos - r)
            r += length
            q += length
        elif op in ('deletion_to_ref', 'skip_ref'):
            if r + length > ref_pos:
                return q
            r += length
        elif op == 'insertion_to_ref':
            q += length
        # soft/hard clips and padding are not part of query_alignment_sequence

    return q


def build_output_sequence(aligned_contigs: list[dict], ref_length: int) -> str:
    """
    Build a reference-ordered output sequence from a set of alignments.

    Alignments are placed in order of reference start position. Reference
    regions not covered by any alignment are filled with 'N' (one N per
    uncovered reference base). Where an alignment overlaps the region already
    covered by a previous one, its overlapping prefix is trimmed (using the
    CIGAR to find the matching query position). Alignments entirely contained
    in the already-covered region are skipped.

    Sets aln['included_in_output'] on every alignment.

    :param aligned_contigs: Aligned segments, as parsed by parsers.parse_sam
    :param ref_length: Length of the reference
    :return: Output sequence
    """
    for aln in aligned_contigs:
        aln['included_in_output'] = False

    ordered = sorted(aligned_contigs, key=lambda a: (a['ref_start'], -a['ref_end']))

    pieces = []
    covered_end = 0
    for aln in ordered:
        ref_start = aln['ref_start']
        ref_end = aln['ref_end']
        if ref_end <= covered_end:
            logger.info(f"Alignment of {aln['contig_id']} ({ref_start}-{ref_end}) is contained in already-covered region; excluding from output")
            continue
        if ref_start < covered_end:
            offset = query_offset_at_ref_pos(aln, covered_end)
            logger.info(f"Alignment of {aln['contig_id']} ({ref_start}-{ref_end}) overlaps previous alignment by {covered_end - ref_start} bp; trimming")
            seq = aln['query_alignment_sequence'][offset:]
        else:
            pieces.append('N' * (ref_start - covered_end))
            seq = aln['query_alignment_sequence']
        pieces.append(seq)
        aln['included_in_output'] = True
        covered_end = ref_end

    pieces.append('N' * max(0, ref_length - covered_end))

    return ''.join(pieces)
