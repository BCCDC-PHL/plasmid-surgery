from __future__ import annotations

import logging

import pysam
import re

from pathlib import Path

logger = logging.getLogger(__name__)

def parse_fasta(fasta_path: Path):
    """
    Parse a fasta file into a list of dicts.

    :param fasta_path: Path to the fasta file
    :type fasta_path: Path
    :return: The parsed fasta file. List of dicts with keys: 'seq_id', 'description', 'seq', 'length'
    :rtype: list[dict]
    """
    seqs = []
    seq = {'seq': ""}
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if 'seq_id' in seq:
                    seq['length'] = len(seq['seq'])
                    seqs.append(seq)
                seq = {'seq': ""}
                defline = line.lstrip('>')
                seq_id = defline.split()[0]
                seq['seq_id'] = seq_id
                description = ' '.join(defline[len(seq_id):].split())
                seq['description'] = description
            else:
                seq['seq'] += line

    if 'seq_id' in seq:
        seq['length'] = len(seq['seq'])
        seqs.append(seq)

    return seqs


def parse_cigar(cigar_str: str) -> list[tuple]:
    """
    Parse a CIGAR string into a list of (operation, length) tuples.
    Raises ValueError if the string is malformed.

    :param cigar_str:
    :type cigar_str:
    :return:
    :rtype:
    """
    cigar_regex = re.compile(r'(\d+)([MIDNSHP=X])')
    matches = cigar_regex.findall(cigar_str)
    total_len = sum(len(m[0]) + 1 for m in matches)
    if total_len != len(cigar_str):
        raise ValueError(f"Invalid CIGAR string: {cigar_str}")

    cigar_operation_codes = {
        'M': 'alignment_match',
        'I': 'insertion_to_ref',
        'D': 'deletion_to_ref',
        'N': 'skip_ref',
        'S': 'soft_clip',
        'H': 'hard_clip',
        'P': 'padding',
        '=': 'sequence_match',
        'X': 'sequence_mismatch',
        'B': 'back'
    }
    
    parsed_cigar = [(cigar_operation_codes[op], int(length)) for length, op in matches]

    return parsed_cigar


# pysam CIGAR operation codes
CIGAR_MATCH = 0       # M
CIGAR_INS = 1         # I
CIGAR_DEL = 2         # D
CIGAR_SEQ_MATCH = 7   # =
CIGAR_SEQ_MISMATCH = 8  # X


def calculate_percent_identity(cigartuples) -> float:
    """
    Gap-compressed percent identity: matches / (matches + mismatches + indel events).

    Each insertion or deletion counts once regardless of its length (as in
    minimap2's 'de' tag), so a single large structural indel doesn't swamp the
    identity of an otherwise near-identical alignment. Indel content is assessed
    separately by the gap fraction in alignment.alignment_qc.

    Requires a CIGAR with explicit '=' / 'X' operations (minimap2 --eqx).

    :param cigartuples: pysam-style list of (op_code, length) tuples
    :return: Percent identity, rounded to 3 decimal places (0.0 if no alignment)
    :rtype: float
    """
    if not cigartuples:
        return 0.0
    matches = sum(length for op, length in cigartuples if op == CIGAR_SEQ_MATCH)
    mismatches = sum(length for op, length in cigartuples if op in (CIGAR_SEQ_MISMATCH, CIGAR_MATCH))
    indel_events = sum(1 for op, _ in cigartuples if op in (CIGAR_INS, CIGAR_DEL))
    denominator = matches + mismatches + indel_events
    if denominator == 0:
        return 0.0

    return round(100 * matches / denominator, 3)


def parse_sam(sam_path) -> list[dict]:
    """
    Parse a sam file into a list of dicts.
    """
    samfile = pysam.AlignmentFile(sam_path, "r")
    aligned_segments = []
    for r in samfile:
        contig_id = r.query_name
        ref_id = r.reference_name  # None if unmapped
        ref_start = r.reference_start
        ref_end = r.reference_end
        query_start = r.query_alignment_start
        query_end = r.query_alignment_end
        query_alignment_length = r.query_alignment_length
        query_alignment_sequence = r.query_alignment_sequence
        is_reverse = r.is_reverse
        is_supplementary = r.is_supplementary
        is_secondary = r.is_secondary
        is_mapped = r.is_mapped
        cigar_str = r.cigarstring
        cigar_parsed = []
        if cigar_str:
            cigar_parsed = parse_cigar(cigar_str)
        tags = dict(r.tags)
        cg = tags.get("cg", None)
        edit_distance = tags.get("NM", None)

        percent_identity = calculate_percent_identity(r.cigartuples)

        aligned_segment = {
            'contig_id': contig_id,
            'ref_id': ref_id,
            'ref_start': ref_start,
            'ref_end': ref_end,
            'query_start': query_start,
            'query_end': query_end,
            'query_alignment_length': query_alignment_length,
            'query_alignment_sequence': query_alignment_sequence,
            'is_reverse': is_reverse,
            'is_supplementary': is_supplementary,
            'is_secondary': is_secondary,
            'alignment_score': tags.get('AS', None),
            'is_mapped': is_mapped,
            'mapq': r.mapping_quality,
            'cigar_str': cigar_str,
            'cigar_parsed': cigar_parsed,
            'tags': tags,
            'edit_distance': edit_distance,
            'percent_identity': percent_identity,
        }
        aligned_segments.append(aligned_segment)

    samfile.close()

    return aligned_segments
