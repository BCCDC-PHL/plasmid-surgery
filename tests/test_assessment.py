import pytest

from plasmid_surgery import assessment
from plasmid_surgery.parsers import parse_cigar


def make_aln(contig_id, ref_start, cigar, contig_length, is_reverse=False, mapq=60):
    ops = parse_cigar(cigar)
    ref_span = sum(l for op, l in ops if op in ('sequence_match', 'sequence_mismatch', 'deletion_to_ref'))
    return {
        'contig_id': contig_id,
        'ref_start': ref_start,
        'ref_end': ref_start + ref_span,
        'cigar_parsed': ops,
        'query_full_length': contig_length,
        'is_reverse': is_reverse,
        'is_mapped': True,
        'mapq': mapq,
    }


def contig(seq_id, length):
    return {'seq_id': seq_id, 'seq': 'A' * length}


# --- coverage and gaps ---

def test_coverage_segments_depth():
    alns = [make_aln('a', 10, '20=', 20), make_aln('b', 20, '20=', 20)]
    assert assessment.coverage_segments(alns, 50) == [
        (0, 10, 0), (10, 20, 1), (20, 30, 2), (30, 40, 1), (40, 50, 0),
    ]


def test_coverage_segments_merges_abutting_equal_depth():
    alns = [make_aln('a', 0, '10=', 10), make_aln('b', 10, '10=', 10)]
    assert assessment.coverage_segments(alns, 20) == [(0, 20, 1)]


def test_coverage_segments_wrapping_alignment():
    alns = [make_aln('a', 45, '10=', 10)]   # covers 45..50 and 0..5
    assert assessment.coverage_segments(alns, 50) == [(0, 5, 1), (5, 45, 0), (45, 50, 1)]


def test_find_gaps_merges_gap_across_origin():
    alns = [make_aln('a', 10, '20=', 20)]
    segments = assessment.coverage_segments(alns, 50)
    assert assessment.find_gaps(segments, 50) == [(30, 60)]


def test_find_gaps_no_coverage():
    segments = assessment.coverage_segments([], 50)
    assert assessment.find_gaps(segments, 50) == [(0, 50)]


def test_find_gaps_full_coverage():
    segments = assessment.coverage_segments([make_aln('a', 0, '50=', 50)], 50)
    assert assessment.find_gaps(segments, 50) == []


def test_write_gaps_bed_splits_wrapping_gap(tmp_path):
    path = tmp_path / 'gaps.bed'
    assessment.write_gaps_bed(path, [(5, 10), (30, 60)], 'ref', 50)
    assert path.read_text().splitlines() == [
        'ref\t5\t10\tgap_1',
        'ref\t30\t50\tgap_2',
        'ref\t0\t10\tgap_2',
    ]


# --- contig coordinates ---

@pytest.mark.parametrize("cigar,is_reverse,expected", [
    ('100=', False, (0, 100)),
    ('10S80=10S', False, (10, 90)),
    ('10S70=20S', True, (20, 90)),      # reverse: flip into contig orientation
    ('30H70=', False, (30, 100)),       # hard clip (supplementary) counted
    ('30H70=', True, (0, 70)),
    ('10S40=5I40=5S', False, (10, 95)),  # insertions are part of the aligned interval
])
def test_query_interval_on_contig(cigar, is_reverse, expected):
    aln = make_aln('a', 0, cigar, 100, is_reverse=is_reverse)
    assert assessment.query_interval_on_contig(aln) == expected


def test_complement_intervals():
    assert assessment.complement_intervals([(10, 20), (15, 30), (50, 60)], 70) == [(0, 10), (30, 50), (60, 70)]
    assert assessment.complement_intervals([], 5) == [(0, 5)]
    assert assessment.complement_intervals([(0, 5)], 5) == []


# --- indel events ---

def test_indel_events_forward():
    aln = make_aln('a', 100, '10S50=60I50=70D50=10S', 230)
    events = assessment.indel_events(aln, ref_length=1000, min_length=50)
    assert events == [
        {'contig_id': 'a', 'event_type': 'insertion', 'ref_start': 150, 'ref_end': 150, 'length': 60, 'contig_position': 60},
        {'contig_id': 'a', 'event_type': 'deletion', 'ref_start': 200, 'ref_end': 270, 'length': 70, 'contig_position': 170},
    ]


def test_indel_events_reverse_contig_position():
    # Same alignment, contig reverse-complemented. Insertion occupies alignment-
    # orientation query 60..120 -> contig 230-120..230-60 = 110..170
    aln = make_aln('a', 100, '10S50=60I50=70D50=10S', 230, is_reverse=True)
    events = assessment.indel_events(aln, ref_length=1000, min_length=50)
    assert events[0]['contig_position'] == 110
    assert events[1]['contig_position'] == 230 - 170


def test_indel_events_below_threshold_ignored():
    aln = make_aln('a', 0, '50=4D50=3I50=', 153)
    assert assessment.indel_events(aln, ref_length=1000, min_length=50) == []


def test_indel_events_fold_across_origin():
    # Deletion starting at ref 990 on a length-1000 ref, alignment wrapping origin
    aln = make_aln('a', 940, '50=60D50=', 100)
    events = assessment.indel_events(aln, ref_length=1000, min_length=50)
    assert events[0]['ref_start'] == 990
    assert events[0]['ref_end'] == 50


# --- assess ---

def test_assess():
    contigs = [contig('a', 100), contig('b', 100), contig('c', 300)]
    alns = [
        make_aln('a', 0, '100=', 100),
        make_aln('b', 150, '60=40S', 100, mapq=0),   # 40 bp unaligned tail
    ]
    result = assessment.assess(contigs, alns, ref_length=300, min_indel_length=50, min_unaligned_length=30)
    m = result['metrics']
    assert m['ref_bases_covered'] == 160
    assert m['ref_fraction_covered'] == round(160 / 300, 4)
    assert m['ref_bases_covered_multiple'] == 0
    assert result['gaps'] == [(100, 150), (210, 300)]
    assert m['num_gaps'] == 2
    assert m['total_gap_length'] == 140
    assert m['largest_gap_length'] == 90
    assert m['contig_bases_aligned'] == 160
    assert m['num_unaligned_contigs'] == 1
    assert [(r['contig_id'], r['start'], r['end']) for r in result['unaligned_regions']] == [('b', 60, 100), ('c', 0, 300)]
    assert m['unaligned_region_bases'] == 340
    assert m['num_ambiguous_alignments'] == 1
    assert result['contig_fraction_aligned'] == {'a': 1.0, 'b': 0.6, 'c': 0.0}


def test_assess_ignores_unmapped():
    alns = [{'contig_id': 'a', 'is_mapped': False}]
    result = assessment.assess([contig('a', 100)], alns, ref_length=200)
    assert result['metrics']['ref_bases_covered'] == 0
    assert result['metrics']['num_unaligned_contigs'] == 1
