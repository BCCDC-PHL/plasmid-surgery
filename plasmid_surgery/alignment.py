from __future__ import annotations

import json
import os
import logging
import subprocess

from pathlib import Path

import plasmid_surgery.parsers as parsers

logger = logging.getLogger(__name__)


MAX_SECONDARY_ALIGNMENTS = 10
# minimap2's own filter on secondary alignments (score relative to primary).
# Secondaries are filtered more strictly in select_alternative_placements.
MIN_SECONDARY_SCORE_RATIO_MINIMAP2 = 0.9

COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


class AlignmentError(Exception):
    pass


def revcomp(seq: str) -> str:
    return seq[::-1].translate(COMPLEMENT)


CLIP_OPS = ('soft_clip', 'hard_clip')
QUERY_ALIGNED_OPS = ('alignment_match', 'sequence_match', 'sequence_mismatch', 'insertion_to_ref')


def aligned_query_sequence(aln: dict, contig_seq: str) -> str:
    """
    The aligned part of the contig, in alignment orientation, derived from the
    CIGAR. Used for alignments whose SAM record has no SEQ (secondaries).
    """
    oriented = revcomp(contig_seq) if aln.get('is_reverse') else contig_seq
    start = 0
    for op, length in aln['cigar_parsed']:
        if op not in CLIP_OPS:
            break
        start += length
    aligned = sum(length for op, length in aln['cigar_parsed'] if op in QUERY_ALIGNED_OPS)
    return oriented[start:start + aligned]


def write_fasta(path, records: list[tuple[str, str]]) -> None:
    with open(path, 'w') as f:
        for seq_id, seq in records:
            f.write(f">{seq_id}\n{seq}\n")


def align_contigs_vs_ref(contigs: list[dict], ref: dict, tmpdir: Path, label: str = 'contigs',
                         threads: int = 1) -> list[list[dict]]:
    """
    Align contigs vs ref using minimap2, in a single call.

    minimap2 aligns each query independently, so this gives the same alignments
    as one call per contig. Queries keep their contig IDs: minimap2 breaks ties
    between equally good placements (e.g. repeat copies, MAPQ 0) with a hash of
    the query name, so renaming them would change which copy is primary. A
    repeated ID is made unique for the call ('<id>__<index>'); `contig_id` in the
    results is always the contig's own ID.

    :param contigs: Contigs, as parsed by parsers.parse_fasta
    :param ref: Reference, as parsed by parsers.parse_fasta
    :param tmpdir: Directory for intermediate files
    :param label: Name for this set of contigs' intermediate files (<label>.fa, <label>_vs_ref[_rotated].sam)
    :param threads: minimap2 threads
    :return: For each contig (in input order), its aligned segments, as parsed by parsers.parse_sam
    :raises AlignmentError: if minimap2 is missing or fails
    """
    if not contigs:
        return []
    os.makedirs(tmpdir, exist_ok=True)

    query_names = []
    index_by_name = {}
    for i, c in enumerate(contigs):
        name = c['seq_id']
        while name in index_by_name:
            name = f"{name}__{i}"
        index_by_name[name] = i
        query_names.append(name)
    contigs_path = os.path.join(tmpdir, f"{label}.fa")
    write_fasta(contigs_path, [(name, c['seq']) for name, c in zip(query_names, contigs)])

    # Always (re)write the ref, so a different ref left over in tmpdir from a
    # previous run is never reused.
    ref_name = 'ref_rotated' if 'rotation' in ref else 'ref'
    ref_path = os.path.join(tmpdir, f"{ref_name}.fa")
    write_fasta(ref_path, [(ref['seq_id'], ref['seq'])])

    alignment_path = os.path.join(tmpdir, f"{label}_vs_{ref_name}.sam")

    minimap2_cmd = [
        'minimap2',
        '-a',
        '-cx',
        'asm10',
        # Report secondary alignments, so that contigs which could be placed at
        # more than one location (e.g. collapsed repeats) can be identified.
        '-N', str(MAX_SECONDARY_ALIGNMENTS),
        '-p', str(MIN_SECONDARY_SCORE_RATIO_MINIMAP2),
        '--eqx',
        '-t', str(threads),
        '-o', alignment_path,
        ref_path,
        contigs_path,
    ]

    try:
        subprocess.run(minimap2_cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError as e:
        raise AlignmentError("minimap2 not found on PATH") from e
    except subprocess.CalledProcessError as e:
        raise AlignmentError(
            f"minimap2 failed (exit code {e.returncode}) aligning {label} vs. {ref['seq_id']}:\n{e.stderr}"
        ) from e

    by_contig = [[] for _ in contigs]
    ref_length = len(ref['seq'])
    for aligned_segment in parsers.parse_sam(alignment_path):
        i = index_by_name[aligned_segment['contig_id']]
        contig = contigs[i]
        aligned_segment['contig_id'] = contig['seq_id']
        aligned_segment['query_full_length'] = len(contig['seq'])
        aligned_segment['ref_full_length'] = ref_length
        if aligned_segment.get('is_mapped') and aligned_segment['query_alignment_sequence'] is None:
            # minimap2 omits SEQ for secondary alignments; recover it from the contig
            aligned_segment['query_alignment_sequence'] = aligned_query_sequence(aligned_segment, contig['seq'])
            aligned_segment['query_alignment_length'] = len(aligned_segment['query_alignment_sequence'])
        by_contig[i].append(aligned_segment)

    return by_contig


REF_ALIGNED_OPS = ('alignment_match', 'sequence_match', 'sequence_mismatch')
REF_ONLY_OPS = ('deletion_to_ref', 'skip_ref')


def rotate_ref(ref: dict, shift: int) -> dict:
    """
    Rotate a circular reference so that position `shift` becomes position 0.

    :param ref: Reference, as parsed by parsers.parse_fasta
    :param shift: Number of bases to rotate by
    :return: Rotated reference dict (with a distinct seq_id, and 'rotation' set to shift)
    """
    seq = ref['seq']
    return {
        'seq_id': f"{ref['seq_id']}__rotated_{shift}",
        'seq': seq[shift:] + seq[:shift],
        'rotation': shift,
    }


def unrotate_alignment(aln: dict, shift: int, ref_length: int, ref_id: str) -> None:
    """
    Convert an alignment against a reference rotated by `shift` back into the
    coordinates of the original reference, in place.

    After conversion ref_start is in [0, ref_length) and ref_end may exceed
    ref_length, meaning the alignment wraps around the origin.
    """
    aln['ref_id'] = ref_id
    if not aln.get('is_mapped'):
        return
    span = aln['ref_end'] - aln['ref_start']
    aln['ref_start'] = (aln['ref_start'] + shift) % ref_length
    aln['ref_end'] = aln['ref_start'] + span


def total_aligned_query_bases(aligned_segments: list[dict]) -> int:
    return sum(a['query_alignment_length'] for a in aligned_segments if a.get('is_mapped'))


def choose_circular_alignment(original: list[dict], rotated: list[dict], ref_length: int) -> list[dict]:
    """
    One contig's alignments against a circular reference: of its alignments
    against the reference and against a copy rotated by half its length (already
    converted back with unrotate_alignment), keep those covering more of the
    contig in fewer segments, preferring the unrotated alignment on ties.

    Sets 'used_rotated_ref' and 'spans_origin' on the chosen alignments.
    """
    def score(segments):
        # Primary + supplementary only: secondaries are alternative placements
        mapped = [a for a in segments if a.get('is_mapped') and not a.get('is_secondary')]
        return (total_aligned_query_bases(mapped), -len(mapped))

    use_rotated = score(rotated) > score(original)
    chosen = rotated if use_rotated else original
    if use_rotated and chosen:
        logger.info(f"Contig {chosen[0]['contig_id']} aligns better across the ref origin; using rotated-ref alignment")

    for aln in chosen:
        aln['used_rotated_ref'] = use_rotated
        aln['spans_origin'] = bool(aln.get('is_mapped')) and aln['ref_end'] > ref_length

    return chosen


def align_contigs_vs_circular_ref(contigs: list[dict], ref: dict, tmpdir: Path, label: str = 'contigs',
                                  threads: int = 1) -> list[list[dict]]:
    """
    Align contigs against a circular reference.

    minimap2 treats the reference as linear, so a contig spanning the
    reference origin is split into primary + supplementary alignments, and
    bases near the origin may be left unaligned. To avoid that, the contigs are
    aligned against both the reference and a copy rotated by half its length
    (one minimap2 call each), and for each contig the result covering more of
    it in fewer segments is kept (see choose_circular_alignment). Coordinates
    are always reported relative to the original reference; an alignment with
    ref_end > ref length wraps around the origin.

    The reference is not doubled, because a doubled reference would make every
    alignment ambiguous (MAPQ 0).

    :param contigs: Contigs, as parsed by parsers.parse_fasta
    :param ref: Circular reference, as parsed by parsers.parse_fasta
    :param tmpdir: Directory for intermediate files
    :param label: Name for this set of contigs' intermediate files
    :param threads: minimap2 threads
    :return: For each contig (in input order), its aligned segments, with
             'spans_origin' and 'used_rotated_ref' set
    """
    ref_length = len(ref['seq'])
    shift = ref_length // 2

    original = align_contigs_vs_ref(contigs, ref, tmpdir, label, threads)
    rotated = align_contigs_vs_ref(contigs, rotate_ref(ref, shift), tmpdir, label, threads)
    for alns in rotated:
        for aln in alns:
            unrotate_alignment(aln, shift, ref_length, ref['seq_id'])

    return [choose_circular_alignment(o, r, ref_length) for o, r in zip(original, rotated)]


def select_alternative_placements(aligned_segments: list[dict], min_score_ratio: float = 0.95) -> tuple[list[dict], list[dict], list[dict]]:
    """
    Separate one contig's alignments into its placement (primary + supplementary)
    and alternative placements (secondary alignments that are nearly as good as
    the primary, e.g. other copies of a repeat).

    A secondary alignment counts as an alternative placement only if its
    alignment score is at least min_score_ratio of the primary's, so that a
    unique contig which merely contains part of a repeat isn't treated as
    ambiguous.

    Sets aln['num_placements'] (1 + number of alternative placements) on all
    returned alignments.

    :param aligned_segments: All mapped alignments of a single contig
    :param min_score_ratio: Minimum secondary/primary alignment score ratio
    :return: (placement, alternatives, rejected) lists of alignments
    """
    placement = [a for a in aligned_segments if not a.get('is_secondary')]
    secondaries = [a for a in aligned_segments if a.get('is_secondary')]
    primary_scores = [a.get('alignment_score') for a in placement
                      if not a.get('is_supplementary') and a.get('alignment_score') is not None]
    alternatives = []
    rejected = []
    for aln in secondaries:
        score = aln.get('alignment_score')
        if primary_scores and score is not None and score >= min_score_ratio * primary_scores[0]:
            alternatives.append(aln)
        else:
            rejected.append(aln)

    for aln in placement + alternatives + rejected:
        aln['num_placements'] = 1 + len(alternatives)

    return placement, alternatives, rejected


def align_contigs(contigs: list[dict], ref: dict, tmpdir: Path, label: str = 'contigs', threads: int = 1) -> dict:
    """
    Align contigs against a circular reference, apply alignment QC, and
    separate each contig's placement from its alternative placements.

    :param contigs: Contigs, as parsed by parsers.parse_fasta
    :param ref: Circular reference, as parsed by parsers.parse_fasta
    :param tmpdir: Directory for intermediate files
    :param label: Name for this set of contigs' intermediate files
    :param threads: minimap2 threads
    :return: dict with lists of aligned segments:
             'placements' (primary + supplementary, passed QC),
             'alternatives' (secondaries that are alternative placements, passed QC),
             'failed' (failed QC, unmapped, or secondaries too weak to be alternatives)
    :raises AlignmentError: if minimap2 fails
    """
    placements = []
    alternatives = []
    failed = []
    logger.debug(f"Aligning {len(contigs)} contig(s) ({label}) vs. ref: {ref['seq_id']}")
    all_aligned_segments = align_contigs_vs_circular_ref(contigs, ref, tmpdir, label, threads)
    for contig, aligned_segments in zip(contigs, all_aligned_segments):
        contig_id = contig['seq_id']
        if not aligned_segments:
            logger.debug(f"Zero alignments for contig: {contig_id} vs. ref: {ref['seq_id']}")
            continue

        passed = []
        for aligned_segment in aligned_segments:
            qc = alignment_qc(aligned_segment)
            if not qc['overall_qc_pass']:
                if aligned_segment.get('is_mapped'):
                    logger.debug(f"Alignment of {contig_id} failed quality filters: " + json.dumps(qc['qc_stats']))
                aligned_segment['included_in_output'] = False
                failed.append(aligned_segment)
                continue
            passed.append(aligned_segment)

        placement, contig_alternatives, rejected = select_alternative_placements(passed)
        if contig_alternatives:
            logger.debug(f"Contig {contig_id} has {len(contig_alternatives)} alternative placement(s) on the ref (e.g. a repeat)")
        for aln in rejected:
            aln['included_in_output'] = False
        placements += placement
        alternatives += contig_alternatives
        failed += rejected

    return {'placements': placements, 'alternatives': alternatives, 'failed': failed}


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
        if op in REF_ALIGNED_OPS:
            if r + length > ref_pos:
                return q + (ref_pos - r)
            r += length
            q += length
        elif op in REF_ONLY_OPS:
            if r + length > ref_pos:
                return q
            r += length
        elif op == 'insertion_to_ref':
            q += length
        # soft/hard clips and padding are not part of query_alignment_sequence

    return q


def split_alignment_at_ref_pos(aln: dict, ref_pos: int) -> tuple[dict, dict]:
    """
    Split an alignment into two pieces at a reference position.

    The left piece covers [ref_start, ref_pos) and the right piece
    [ref_pos, ref_end). Insertions immediately before ref_pos go to the left
    piece (consistent with query_offset_at_ref_pos). Each piece gets its own
    CIGAR (without clipping) and query_alignment_sequence.

    :param aln: Aligned segment, with ref_start < ref_pos < ref_end
    :param ref_pos: Reference position to split at
    :return: (left, right) aligned-segment dicts
    """
    r = aln['ref_start']
    q = 0
    left_ops = []
    right_ops = []
    split_q = None
    for op, length in aln['cigar_parsed']:
        if op not in REF_ALIGNED_OPS + REF_ONLY_OPS + ('insertion_to_ref',):
            continue  # clipping / padding
        if split_q is not None:
            right_ops.append((op, length))
        elif op == 'insertion_to_ref':
            left_ops.append((op, length))
            q += length
        elif r + length <= ref_pos:
            left_ops.append((op, length))
            r += length
            if op in REF_ALIGNED_OPS:
                q += length
        else:
            k = ref_pos - r
            if k > 0:
                left_ops.append((op, k))
            right_ops.append((op, length - k))
            if op in REF_ALIGNED_OPS:
                q += k
            split_q = q
    if split_q is None:
        split_q = q

    seq = aln['query_alignment_sequence']
    left = dict(aln)
    left.update({
        'ref_end': ref_pos,
        'query_alignment_sequence': seq[:split_q],
        'query_alignment_length': split_q,
        'cigar_parsed': left_ops,
    })
    right = dict(aln)
    right.update({
        'ref_start': ref_pos,
        'query_alignment_sequence': seq[split_q:],
        'query_alignment_length': len(seq) - split_q,
        'cigar_parsed': right_ops,
    })

    return left, right


def _circular_pieces(alignments: list[dict], ref_length: int) -> list[tuple[dict, dict]]:
    """
    (piece, parent) pairs, with alignments wrapping the origin split there and
    the part past the origin shifted to the start of the reference.
    """
    pieces = []
    for aln in alignments:
        if aln['ref_end'] > ref_length:
            left, right = split_alignment_at_ref_pos(aln, ref_length)
            right['ref_start'] -= ref_length
            right['ref_end'] -= ref_length
            pieces += [(left, aln), (right, aln)]
        else:
            pieces.append((aln, aln))
    pieces.sort(key=lambda p: (p[0]['ref_start'], -p[0]['ref_end']))
    return pieces


def _fill_ref_interval(pieces, start: int, end: int, fill_gap, transform=lambda s: s) -> str:
    """
    Sequence for the reference interval [start, end), built from (piece, parent)
    pairs in reference order. The first piece to reach a region wins: later
    overlapping pieces are trimmed and contained pieces are skipped. Uncovered
    regions are filled by fill_gap(gap_start, gap_end).
    """
    out = []
    covered_end = start
    for piece, parent in pieces:
        piece_start = max(piece['ref_start'], covered_end)
        piece_end = min(piece['ref_end'], end)
        if piece_end <= piece_start:
            continue
        if piece['ref_start'] < covered_end:
            logger.debug(f"Alignment of {piece['contig_id']} ({piece['ref_start']}-{piece['ref_end']}) overlaps previous alignment by {covered_end - piece['ref_start']} bp; trimming")
        if piece_start > covered_end:
            out.append(fill_gap(covered_end, piece_start))
        seq = piece['query_alignment_sequence']
        out.append(transform(seq[query_offset_at_ref_pos(piece, piece_start):query_offset_at_ref_pos(piece, piece_end)]))
        parent['included_in_output'] = True
        covered_end = piece_end
    if covered_end < end:
        out.append(fill_gap(covered_end, end))

    return ''.join(out)


def build_output_sequence(aligned_contigs: list[dict], ref_length: int, alternative_placements: list[dict] = ()) -> str:
    """
    Build a reference-ordered output sequence from a set of alignments.

    The output starts at reference position 0. Alignments are placed in order
    of reference start position. Alignments wrapping around the origin
    (ref_end > ref_length) are split at the origin, with the part beyond the
    origin placed at the start of the output. Where an alignment overlaps the
    region already covered by a previous one, its overlapping prefix is trimmed
    (using the CIGAR to find the matching query position). Alignments entirely
    contained in the already-covered region are skipped.

    Reference regions not covered by any alignment are filled from
    alternative placements (e.g. other copies of a collapsed repeat), in
    lowercase to mark them as inferred, and otherwise with 'N' (one N per
    uncovered reference base).

    Sets aln['included_in_output'] on every alignment.

    :param aligned_contigs: Placement alignments (primary + supplementary)
    :param ref_length: Length of the reference
    :param alternative_placements: Secondary alignments used to fill gaps
    :return: Output sequence
    """
    for aln in list(aligned_contigs) + list(alternative_placements):
        aln['included_in_output'] = False

    primary_pieces = _circular_pieces(aligned_contigs, ref_length)
    alternative_pieces = _circular_pieces(alternative_placements, ref_length)

    def fill_with_ns(start, end):
        return 'N' * (end - start)

    def fill_from_alternatives(start, end):
        return _fill_ref_interval(alternative_pieces, start, end, fill_with_ns, transform=str.lower)

    return _fill_ref_interval(primary_pieces, 0, ref_length, fill_from_alternatives)
