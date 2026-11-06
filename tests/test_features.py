import pytest

from plasmid_surgery import features
from plasmid_surgery.parsers import parse_cigar

HAMRONIZATION_HEADER = ['input_file_name', 'gene_symbol', 'gene_name', 'reference_database_name',
                        'reference_database_version', 'reference_accession', 'analysis_software_name',
                        'analysis_software_version', 'genetic_variation_type', 'drug_class',
                        'input_gene_start', 'input_gene_stop', 'input_sequence_id', 'strand_orientation']


def write_hamronization(path, rows):
    lines = ['\t'.join(HAMRONIZATION_HEADER)]
    for gene, start, stop, seq_id, strand, drug_class in rows:
        values = {'input_file_name': 'x.fa', 'gene_symbol': gene, 'gene_name': f"{gene} product",
                  'reference_database_name': 'ncbi', 'reference_database_version': '1',
                  'reference_accession': 'acc', 'analysis_software_name': 'abricate',
                  'analysis_software_version': '1.4.0', 'genetic_variation_type': 'gene_presence_detected',
                  'drug_class': drug_class, 'input_gene_start': str(start), 'input_gene_stop': str(stop),
                  'input_sequence_id': seq_id, 'strand_orientation': strand}
        lines.append('\t'.join(values[h] for h in HAMRONIZATION_HEADER))
    path.write_text('\n'.join(lines) + '\n')


# --- parsing ---

def test_parse_hamronization(tmp_path):
    path = tmp_path / 'h.tsv'
    write_hamronization(path, [('blaX', 101, 200, 'contig_1', '+', 'BETA-LACTAM'),
                               ('repA', 300, 251, 'contig_2', '-', '')])
    assert features.detect_format(path) == 'hamronization'
    [a, b] = features.parse_hamronization(path)
    assert (a['seq_id'], a['start'], a['end'], a['strand']) == ('contig_1', 100, 200, '+')
    assert (a['name'], a['product'], a['category'], a['source']) == ('blaX', 'blaX product', 'BETA-LACTAM', 'abricate:ncbi')
    assert (b['start'], b['end'], b['strand']) == (250, 300, '-')   # reversed coordinates normalized


def test_empty_file_has_no_features(tmp_path):
    # hAMRonization writes an empty file (no header) when the tool found nothing
    path = tmp_path / 'empty.tsv'
    path.write_text('')
    assert features.detect_format(path) == 'empty'
    assert features.read_features([path]) == ([], {})


def test_parse_gff3_drops_parent_genes_decodes_attributes_reads_fasta(tmp_path):
    path = tmp_path / 'a.gff3'
    path.write_text(
        "##gff-version 3\n"
        "contig_1\tBakta\tregion\t1\t100\t.\t+\t.\tID=contig_1\n"
        "contig_1\tBakta\tgene\t11\t40\t.\t+\t.\tID=g1;gene=traA\n"
        "contig_1\tBakta\tCDS\t11\t40\t.\t+\t0\tID=c1;Parent=g1;gene=traA;product=conjugal%20transfer%2C TraA\n"
        "contig_1\tBakta\tmobile_genetic_element\t51\t90\t.\t-\t.\tID=m1;Name=IS26\n"
        "##FASTA\n>contig_1\nACGT\nACGT\n"
    )
    assert features.detect_format(path) == 'gff3'
    found, sequences = features.parse_gff3(path)
    assert [(f['name'], f['category'], f['start'], f['end'], f['strand']) for f in found] == [
        ('traA', 'CDS', 10, 40, '+'), ('IS26', 'mobile_genetic_element', 50, 90, '-')]
    assert found[0]['product'] == 'conjugal transfer, TraA'
    assert sequences == {'contig_1': 'ACGTACGT'}


def test_parse_bed(tmp_path):
    path = tmp_path / 'f.bed'
    path.write_text("track name=x\nref\t10\t20\tgeneA\t0\t-\nref\t30\t40\n")
    assert features.detect_format(path) == 'bed'
    [a, b] = features.parse_bed(path)
    assert (a['name'], a['start'], a['end'], a['strand']) == ('geneA', 10, 20, '-')
    assert (b['name'], b['strand']) == ('ref:30-40', '.')


def test_unrecognized_format(tmp_path):
    path = tmp_path / 'x.txt'
    path.write_text("hello world\n")
    with pytest.raises(features.FeatureFormatError):
        features.detect_format(path)


# --- mapping ---

def feature(seq_id, start, end, strand='+', name='g'):
    return features._feature(seq_id, start, end, strand, name)


def test_map_features_by_id_and_unversioned_id():
    targets = [{'seq_id': 'CP0001.1', 'seq': 'ACGTACGTAA'}]
    mapped, unmapped = features.map_features([feature('CP0001.1', 0, 2), feature('CP0001', 2, 4), feature('other', 0, 1)], targets, {})
    assert [m['target_id'] for m in mapped] == ['CP0001.1', 'CP0001.1']
    assert unmapped == 1


def test_map_features_by_sequence_with_flip():
    # Feature annotated on 'asm_1' (e.g. by abricate on the assembly); the target
    # is the same contig renamed and reverse-complemented (e.g. by TaDReP)
    asm_seq = 'AAACCCGGGTTTA'
    rc = asm_seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    targets = [{'seq_id': 'tadrep-asm_1', 'seq': rc}]
    [m], unmapped = features.map_features([feature('asm_1', 3, 6, '+')], targets, {'asm_1': asm_seq})
    assert unmapped == 0
    assert (m['target_id'], m['start'], m['end'], m['strand']) == ('tadrep-asm_1', len(asm_seq) - 6, len(asm_seq) - 3, '-')
    assert rc[m['start']:m['end']] == asm_seq[3:6][::-1].translate(str.maketrans('ACGT', 'TGCA'))


# --- coordinates ---

def aln(ref_start, cigar, contig_length, is_reverse=False):
    ops = parse_cigar(cigar)
    span = sum(l for op, l in ops if op in ('sequence_match', 'sequence_mismatch', 'deletion_to_ref'))
    return {'contig_id': 'c', 'ref_start': ref_start, 'ref_end': ref_start + span, 'cigar_parsed': ops,
            'query_full_length': contig_length, 'is_reverse': is_reverse, 'is_mapped': True}


def test_matched_query_intervals_exclude_insertions():
    assert features.matched_query_intervals_on_contig(aln(0, '10S40=20I40=', 110)) == [(10, 50), (70, 110)]
    # reverse: flipped into contig orientation
    assert features.matched_query_intervals_on_contig(aln(0, '10S40=20I40=', 110, is_reverse=True)) == [(0, 40), (60, 100)]


@pytest.mark.parametrize("contig_pos,is_reverse,expected", [
    (10, False, 100),      # first aligned base
    (49, False, 139),
    (55, False, None),     # inside the insertion
    (70, False, 140),      # first base after the insertion
    (75, False, 150),      # after the 5 bp deletion following it
    (99, True, 100),       # reverse: the clip is at the contig's far end, so contig 99 is the first aligned base
    (100, True, None),     # clipped
])
def test_ref_pos_at_contig_pos(contig_pos, is_reverse, expected):
    # 10S 40= 20I 5= 5D 35=  on ref from 100
    a = aln(100, '10S40=20I5=5D35=', 110, is_reverse=is_reverse)
    assert features.ref_pos_at_contig_pos(a, contig_pos, ref_length=1000) == expected


def test_ref_pos_folds_across_origin():
    a = aln(990, '20=', 20)
    assert features.ref_pos_at_contig_pos(a, 15, ref_length=1000) == 5


# --- assessment ---

def placement(ref_start, length, contig_id='c'):
    a = aln(ref_start, f'{length}=', length)
    a['contig_id'] = contig_id
    return a


def test_assess_ref_features_statuses():
    reconstruction = {'placements': [placement(0, 100)], 'alternatives': [placement(200, 50)]}
    assembly = reconstruction['placements'] + [placement(300, 100, contig_id='chrom_1')]
    ref_features = [feature('r', 10, 50, name='present'), feature('r', 90, 110, name='partial'),
                    feature('r', 210, 240, name='inferred'), feature('r', 310, 350, name='missed'),
                    feature('r', 450, 480, name='absent')]
    rows, metrics = features.assess_ref_features(
        ref_features, 500, reconstruction, assembly_alignments=assembly,
        labels={'c': 'reconstruction', 'chrom_1': 'chromosome'})
    by_name = {r['name']: r for r in rows}
    assert {n: r['reconstruction_status'] for n, r in by_name.items()} == {
        'present': 'present', 'partial': 'partial', 'inferred': 'inferred', 'missed': 'absent', 'absent': 'absent'}
    assert by_name['partial']['reconstruction_fraction'] == 0.5
    assert by_name['inferred']['inferred_fraction'] == 1.0
    assert (by_name['missed']['assembly_status'], by_name['missed']['found_in']) == ('present', 'chromosome')
    assert (by_name['absent']['assembly_status'], by_name['absent']['found_in']) == ('absent', '-')
    assert metrics == {'num_ref_features': 5, 'num_ref_features_present': 1, 'num_ref_features_partial': 1,
                       'num_ref_features_absent': 2, 'num_ref_features_inferred': 1, 'num_ref_features_missed': 1}


def test_assess_ref_features_without_assembly():
    rows, metrics = features.assess_ref_features([feature('r', 10, 50)], 500, {'placements': [], 'alternatives': []})
    assert rows[0]['assembly_status'] == 'NA'
    assert metrics['num_ref_features_missed'] == 'NA'
