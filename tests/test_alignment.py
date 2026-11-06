import pytest

from plasmid_surgery import alignment
from plasmid_surgery.parsers import parse_cigar


def make_aln(contig_id, ref_start, seq, cigar=None, ref_end=None):
    """Build a minimal aligned-segment dict. Default CIGAR is all matches."""
    if cigar is None:
        cigar = f"{len(seq)}="
    cigar_parsed = parse_cigar(cigar)
    if ref_end is None:
        ref_span = sum(l for op, l in cigar_parsed
                       if op in ('sequence_match', 'sequence_mismatch', 'alignment_match', 'deletion_to_ref'))
        ref_end = ref_start + ref_span
    return {
        'contig_id': contig_id,
        'ref_start': ref_start,
        'ref_end': ref_end,
        'query_alignment_sequence': seq,
        'query_alignment_length': len(seq),
        'cigar_parsed': cigar_parsed,
        'percent_identity': 100.0,
        'is_mapped': True,
    }


def test_build_output_fills_gaps_with_ref_length_ns():
    a = make_aln('a', 2, 'AAAA')
    b = make_aln('b', 8, 'CC')
    out = alignment.build_output_sequence([b, a], ref_length=12)
    assert out == 'NN' + 'AAAA' + 'NN' + 'CC' + 'NN'
    assert a['included_in_output'] and b['included_in_output']


def test_build_output_with_insertion_does_not_shift_downstream():
    # 'a' has a 3 bp insertion relative to the ref; 'b' must still appear intact
    # (the original splicing approach overwrote/displaced downstream contigs)
    a = make_aln('a', 0, 'AAAGGGAA', cigar='3=3I2=')   # ref span 5
    b = make_aln('b', 5, 'CCCCC')
    c = make_aln('c', 12, 'TTT')
    out = alignment.build_output_sequence([c, a, b], ref_length=15)
    assert out == 'AAAGGGAA' + 'CCCCC' + 'NN' + 'TTT'


def test_build_output_with_deletion():
    a = make_aln('a', 0, 'AAAAA', cigar='2=3D3=')   # ref span 8
    b = make_aln('b', 8, 'CC')
    out = alignment.build_output_sequence([a, b], ref_length=11)
    assert out == 'AAAAA' + 'CC' + 'N'


def test_build_output_trims_overlap():
    # 'b' overlaps the last 3 bp of 'a' on the reference
    a = make_aln('a', 0, 'AAAAAGGG')
    b = make_aln('b', 5, 'GGGCCCC')
    out = alignment.build_output_sequence([a, b], ref_length=12)
    assert out == 'AAAAAGGGCCCC'
    assert b['included_in_output']


def test_build_output_trims_overlap_with_indel_in_overlap():
    a = make_aln('a', 0, 'AAAAAGGG')                    # ref 0-8
    b = make_aln('b', 5, 'GGTGCCCC', cigar='2=1I5=')    # ref 5-12, insertion inside overlap
    out = alignment.build_output_sequence([a, b], ref_length=12)
    assert out == 'AAAAAGGG' + 'CCCC'


def test_build_output_skips_contained_alignment():
    a = make_aln('a', 0, 'AAAAAAAAAA')
    b = make_aln('b', 3, 'AAA')
    out = alignment.build_output_sequence([a, b], ref_length=10)
    assert out == 'AAAAAAAAAA'
    assert a['included_in_output']
    assert not b['included_in_output']


def test_build_output_empty():
    assert alignment.build_output_sequence([], ref_length=5) == 'NNNNN'


# query: 3 match, 3 ins, 2 match, 3 del, 2 match  => query len 10, ref 10..20
#   ref:   10 11 12  -  -  - 13 14 15 16 17 18 19
#   query:  0  1  2  3  4  5  6  7  -  -  -  8  9
@pytest.mark.parametrize("ref_pos,expected", [
    (10, 0),   # at ref_start
    (12, 2),   # inside first match block
    (13, 6),   # just after insertion -> insertion skipped
    (14, 7),
    (15, 8),   # inside deletion -> next aligned query base
    (17, 8),   # still inside deletion
    (18, 8),   # first base after deletion
    (19, 9),
    (20, 10),  # at ref_end
])
def test_query_offset_at_ref_pos(ref_pos, expected):
    aln = make_aln('a', 10, 'AAAGGGCCTT', cigar='3=3I2=3D2=')
    assert aln['ref_end'] == 20
    assert alignment.query_offset_at_ref_pos(aln, ref_pos) == expected


def test_alignment_qc_zero_length_does_not_crash():
    aln = {'query_alignment_length': 0, 'cigar_parsed': [], 'percent_identity': 0.0, 'is_mapped': False}
    qc = alignment.alignment_qc(aln)
    assert qc['overall_qc_pass'] is False


def test_alignment_qc_passes_good_alignment():
    aln = make_aln('a', 0, 'A' * 1000)
    assert alignment.alignment_qc(aln)['overall_qc_pass'] is True


def test_alignment_qc_fails_low_identity():
    aln = make_aln('a', 0, 'A' * 1000)
    aln['percent_identity'] = 50.0
    assert alignment.alignment_qc(aln)['overall_qc_pass'] is False
