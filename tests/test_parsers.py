import pytest

from plasmid_surgery import parsers


def test_parse_fasta_multiple_records(tmp_path):
    path = tmp_path / "seqs.fa"
    path.write_text(">seq1 some description\nACGT\nAC\n>seq2\nGGGG\n")
    seqs = parsers.parse_fasta(path)
    assert [s['seq_id'] for s in seqs] == ['seq1', 'seq2']
    assert seqs[0]['seq'] == 'ACGTAC'
    assert seqs[0]['length'] == 6
    assert seqs[0]['description'] == 'some description'
    assert seqs[1]['description'] == ''


def test_parse_fasta_description_not_stripped_by_char_set(tmp_path):
    # lstrip(seq_id) would have stripped the leading 'c', 'o', 'n'... from the description
    path = tmp_path / "seqs.fa"
    path.write_text(">contig_1 contig of length 4\nACGT\n")
    seqs = parsers.parse_fasta(path)
    assert seqs[0]['description'] == 'contig of length 4'


def test_parse_fasta_empty_file(tmp_path):
    path = tmp_path / "empty.fa"
    path.write_text("")
    assert parsers.parse_fasta(path) == []


def test_parse_cigar():
    assert parsers.parse_cigar("5S10=1X2I3D4=") == [
        ('soft_clip', 5),
        ('sequence_match', 10),
        ('sequence_mismatch', 1),
        ('insertion_to_ref', 2),
        ('deletion_to_ref', 3),
        ('sequence_match', 4),
    ]


def test_parse_cigar_invalid():
    with pytest.raises(ValueError):
        parsers.parse_cigar("10=Q")


@pytest.mark.parametrize("cigartuples,expected", [
    ([(7, 100)], 100.0),
    ([(7, 90), (8, 10)], 90.0),
    # 90 matches, one 5 bp insertion, one 5 bp deletion -> 90 / (90 + 2)
    ([(7, 45), (1, 5), (7, 45), (2, 5)], 97.826),
    # a single large insertion counts as one event
    ([(7, 999), (1, 3000)], 99.9),
    # soft clips don't count
    ([(4, 50), (7, 100), (4, 50)], 100.0),
    ([], 0.0),
    (None, 0.0),
])
def test_calculate_percent_identity(cigartuples, expected):
    assert parsers.calculate_percent_identity(cigartuples) == expected
