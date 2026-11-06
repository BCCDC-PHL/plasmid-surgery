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


def test_rotate_ref():
    ref = {'seq_id': 'r', 'seq': 'AAAACCCC'}
    rotated = alignment.rotate_ref(ref, 3)
    assert rotated['seq'] == 'ACCCCAAA'
    assert rotated['seq_id'] != 'r'


def test_unrotate_alignment_wraps_origin():
    # ref length 100, rotated by 50: rotated position 40 == original position 90
    aln = make_aln('a', 40, 'A' * 20)
    alignment.unrotate_alignment(aln, shift=50, ref_length=100, ref_id='r')
    assert (aln['ref_start'], aln['ref_end']) == (90, 110)
    assert aln['ref_id'] == 'r'


def test_unrotate_alignment_no_wrap():
    aln = make_aln('a', 60, 'A' * 20)
    alignment.unrotate_alignment(aln, shift=50, ref_length=100, ref_id='r')
    assert (aln['ref_start'], aln['ref_end']) == (10, 30)


@pytest.mark.parametrize("ref_pos", range(11, 20))
def test_split_alignment_consistent_with_query_offset(ref_pos):
    aln = make_aln('a', 10, 'AAAGGGCCTT', cigar='3=3I2=3D2=')
    left, right = alignment.split_alignment_at_ref_pos(aln, ref_pos)
    offset = alignment.query_offset_at_ref_pos(aln, ref_pos)
    assert left['query_alignment_sequence'] == aln['query_alignment_sequence'][:offset]
    assert right['query_alignment_sequence'] == aln['query_alignment_sequence'][offset:]
    assert (left['ref_start'], left['ref_end']) == (10, ref_pos)
    assert (right['ref_start'], right['ref_end']) == (ref_pos, 20)
    # Each piece's CIGAR is consistent with its sequence and ref span
    for piece in (left, right):
        q_len = sum(l for op, l in piece['cigar_parsed'] if op in ('sequence_match', 'insertion_to_ref'))
        r_len = sum(l for op, l in piece['cigar_parsed'] if op in ('sequence_match', 'deletion_to_ref'))
        assert q_len == len(piece['query_alignment_sequence'])
        assert r_len == piece['ref_end'] - piece['ref_start']


def test_build_output_wrapping_alignment():
    # 'a' covers ref 8..12 on a length-10 ref, i.e. wraps the origin
    a = make_aln('a', 8, 'GGTT')
    b = make_aln('b', 4, 'CC')
    out = alignment.build_output_sequence([a, b], ref_length=10)
    assert out == 'TT' + 'NN' + 'CC' + 'NN' + 'GG'
    assert a['included_in_output'] and b['included_in_output']


def test_build_output_wrapping_alignment_overlaps_first_contig():
    # wrapped part of 'a' covers ref 0..2, 'b' starts at 1
    a = make_aln('a', 8, 'GGTT')
    b = make_aln('b', 1, 'TCCC')
    out = alignment.build_output_sequence([a, b], ref_length=10)
    assert out == 'TT' + 'CCC' + 'NNN' + 'GG'


# --- repeat awareness ---

def test_build_output_fills_gaps_from_alternative_placements_in_lowercase():
    a = make_aln('a', 0, 'AAAA')
    b = make_aln('b', 10, 'CCCC')
    repeat_alt = make_aln('r', 4, 'GGG')          # fills 4-7; 7-10 stays N
    out = alignment.build_output_sequence([a, b], ref_length=14, alternative_placements=[repeat_alt])
    assert out == 'AAAA' + 'ggg' + 'NNN' + 'CCCC'
    assert repeat_alt['included_in_output']


def test_build_output_alternative_trimmed_to_gap():
    a = make_aln('a', 0, 'AAAA')
    b = make_aln('b', 8, 'CCCC')
    repeat_alt = make_aln('r', 2, 'TTGGGGTT')      # ref 2-10, gap is 4-8
    out = alignment.build_output_sequence([a, b], ref_length=12, alternative_placements=[repeat_alt])
    assert out == 'AAAA' + 'gggg' + 'CCCC'


def test_build_output_alternative_not_used_where_covered():
    a = make_aln('a', 0, 'AAAAAAAA')
    repeat_alt = make_aln('r', 2, 'GGGG')
    out = alignment.build_output_sequence([a], ref_length=8, alternative_placements=[repeat_alt])
    assert out == 'AAAAAAAA'
    assert not repeat_alt['included_in_output']


def test_build_output_alternative_wrapping_origin():
    a = make_aln('a', 2, 'AAAAAA')               # ref 2-8
    repeat_alt = make_aln('r', 8, 'GGTT')        # ref 8-12 on length-10 ref: wraps
    out = alignment.build_output_sequence([a], ref_length=10, alternative_placements=[repeat_alt])
    assert out == 'tt' + 'AAAAAA' + 'gg'


def secondary(contig_id, score, **kwargs):
    aln = make_aln(contig_id, 0, 'A' * 10, **kwargs)
    aln.update({'is_secondary': True, 'is_supplementary': False, 'alignment_score': score})
    return aln


def test_select_alternative_placements():
    primary = make_aln('r', 0, 'A' * 10)
    primary.update({'is_secondary': False, 'is_supplementary': False, 'alignment_score': 1000})
    supplementary = make_aln('r', 50, 'A' * 10)
    supplementary.update({'is_secondary': False, 'is_supplementary': True, 'alignment_score': 200})
    good = secondary('r', 990)
    marginal = secondary('r', 950)
    poor = secondary('r', 600)   # e.g. a unique contig that just contains part of a repeat
    placement, alternatives, rejected = alignment.select_alternative_placements(
        [primary, supplementary, good, marginal, poor], min_score_ratio=0.95)
    assert placement == [primary, supplementary]
    assert alternatives == [good, marginal]
    assert rejected == [poor]
    assert all(a['num_placements'] == 3 for a in placement + alternatives + rejected)


def test_select_alternative_placements_no_secondaries():
    primary = make_aln('a', 0, 'A' * 10)
    primary.update({'is_secondary': False, 'is_supplementary': False, 'alignment_score': 10})
    placement, alternatives, rejected = alignment.select_alternative_placements([primary])
    assert (placement, alternatives, rejected) == ([primary], [], [])
    assert primary['num_placements'] == 1


@pytest.mark.parametrize("cigar,is_reverse,expected", [
    ('2S5=1S', False, 'TACCT'),
    ('3H5=', False, 'ACCTT'),
    # reverse: slice of revcomp(contig) = 'AAGGTACG'
    ('2S5=1S', True, 'GGTAC'),
])
def test_aligned_query_sequence(cigar, is_reverse, expected):
    contig_seq = 'CGTACCTT'
    aln = {'cigar_parsed': parse_cigar(cigar), 'is_reverse': is_reverse}
    assert alignment.aligned_query_sequence(aln, contig_seq) == expected
