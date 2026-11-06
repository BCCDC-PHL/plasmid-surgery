import pytest

from plasmid_surgery import context


def test_sequence_key_orientation_independent():
    assert context.sequence_key('AACGT') == context.sequence_key('ACGTT')   # revcomp
    assert context.sequence_key('aacgt') == context.sequence_key('AACGT')
    assert context.sequence_key('AACGT') != context.sequence_key('AACGA')


@pytest.mark.parametrize("seq_id,description,expected", [
    ('contig_1', 'length=1000 depth=2.34x', 2.34),                   # Unicycler
    ('NODE_3_length_5000_cov_12.5', '', 12.5),                       # SPAdes
    ('contig00001', 'len=5000 cov=31.2 corr=0 origname=NODE_1', 31.2),   # Shovill
    ('contig_1', 'length=1000', None),
])
def test_parse_depth(seq_id, description, expected):
    assert context.parse_depth({'seq_id': seq_id, 'description': description}) == expected


def test_length_weighted_mean_depth():
    contigs = [
        {'seq_id': 'a', 'description': 'depth=1x', 'seq': 'A' * 300},
        {'seq_id': 'b', 'description': 'depth=3x', 'seq': 'A' * 100},
        {'seq_id': 'c', 'description': '', 'seq': 'A' * 1000},   # no depth: ignored
    ]
    assert context.length_weighted_mean_depth(contigs) == pytest.approx(1.5)
    assert context.length_weighted_mean_depth([contigs[2]]) is None


def test_bin_labels(tmp_path):
    (tmp_path / 'plasmid_AA1.fasta').write_text(">c1\nAAAACCCC\n>c2\nGGGGTTTA\n")
    (tmp_path / 'chromosome.fasta').write_text(">c3\nACACACAC\n>c2_again\nTAAACCCC\n")   # revcomp of c2
    labels = context.bin_labels([tmp_path / 'plasmid_AA1.fasta', tmp_path / 'chromosome.fasta'])
    assert labels[context.sequence_key('AAAACCCC')] == 'plasmid_AA1'
    assert labels[context.sequence_key('ACACACAC')] == 'chromosome'
    assert labels[context.sequence_key('GGGGTTTA')] == 'chromosome,plasmid_AA1'


@pytest.mark.parametrize("num_placements,fraction,expected", [
    (1, 1.0, 'missed'),
    (1, 0.8, 'missed'),
    (1, 0.5, 'shared'),
    (3, 1.0, 'repeat'),
    (3, 0.1, 'shared'),   # e.g. a chromosome contig containing an IS element
])
def test_classify(num_placements, fraction, expected):
    assert context.classify(num_placements, fraction, min_missed_fraction=0.8) == expected


def test_ref_mask_wraps():
    alns = [{'ref_start': 8, 'ref_end': 12}, {'ref_start': 3, 'ref_end': 5}]
    assert list(context.ref_mask(alns, 10)) == [1, 1, 0, 1, 1, 0, 0, 0, 1, 1]
