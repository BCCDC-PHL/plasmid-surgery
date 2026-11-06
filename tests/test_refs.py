import pytest

from plasmid_surgery import refs
from tests.scenarios import find_refs as scenario
from tests.scenarios.common import mash_available, minimap2_available, read_tsv, workspace


# --- unit ---

def test_estimate_containment():
    # A (100 k-mers) entirely inside B (400): J = 100/400, so all of A is in B, a quarter of B in A
    in_b, b_in = refs.estimate_containment(0.25, 100, 400)
    assert in_b == pytest.approx(1.0) and b_in == pytest.approx(0.25)
    # Identical sets
    assert refs.estimate_containment(1.0, 500, 500) == (pytest.approx(1.0), pytest.approx(1.0))
    # Half overlapping, equal sizes: |A∩B| = 50, |A∪B| = 150
    in_b, b_in = refs.estimate_containment(50 / 150, 100, 100)
    assert in_b == pytest.approx(0.5) and b_in == pytest.approx(0.5)
    assert refs.estimate_containment(0.0, 100, 100) == (0.0, 0.0)


def row(ref_id, distance, containment):
    return {'ref_id': ref_id, 'mash_distance': distance, 'est_reconstruction_containment': containment}


def test_select_candidates_takes_top_by_containment_and_by_distance():
    rows = [row('near_small', 0.001, 0.5), row('big', 0.05, 1.0), row('mid', 0.01, 0.9), row('far', 0.2, 0.1)]
    chosen = {r['ref_id'] for r in refs.select_candidates(rows, 1)}
    assert chosen == {'big', 'near_small'}   # best containment, and nearest
    assert len(refs.select_candidates(rows, 10)) == 4


def test_score_uses_assembly_coverage_when_present():
    assert refs.score({'reconstruction_fraction_aligned': 0.5, 'ref_fraction_covered': 0.8, 'assembly_ref_fraction_covered': 'NA'}) == 0.4
    assert refs.score({'reconstruction_fraction_aligned': 0.5, 'ref_fraction_covered': 0.8, 'assembly_ref_fraction_covered': 1.0}) == 0.5


def test_rank_ties_broken_by_identity_then_distance():
    rows = [{'ref_id': 'a', 'score': 0.9, 'identity': 99.0, 'mash_distance': 0.01},
            {'ref_id': 'b', 'score': 0.9, 'identity': 99.5, 'mash_distance': 0.02},
            {'ref_id': 'c', 'score': 0.9, 'identity': 99.5, 'mash_distance': 0.01},
            {'ref_id': 'd', 'score': 0.95, 'identity': 90.0, 'mash_distance': 0.1}]
    assert [r['ref_id'] for r in refs.rank(rows)] == ['d', 'c', 'b', 'a']
    assert [r['rank'] for r in rows] == [4, 3, 2, 1]


def test_target_intervals_exclude_large_deletions():
    # PAF columns: target start at 100, CIGAR with a 60 bp deletion
    cols = ['q', '100', '0', '100', '+', 't', '1000', '100', '260', '100', '160', '60', 'cg:Z:50M60D50M']
    assert refs._target_intervals(cols, 50) == [(100, 150), (210, 260)]
    assert refs._target_intervals(cols, 100) == [(100, 260)]


@pytest.fixture
def small_fasta(tmp_path):
    path = tmp_path / 'db.fa'
    path.write_text(">p1 first plasmid\nACGTACGT\nAC\n>p2\nGGGG\n>p3 third\nTTTT\n")
    return path


def test_extract_sequences_without_index(small_fasta):
    found = refs.extract_sequences(small_fasta, ['p1', 'p3', 'absent'])
    assert found == {'p1': {'seq': 'ACGTACGTAC', 'description': 'first plasmid'},
                     'p3': {'seq': 'TTTT', 'description': 'third'}}
    assert not (small_fasta.parent / 'db.fa.fai').exists()   # nothing written next to the FASTA


def test_extract_sequences_with_index(small_fasta, tmp_path):
    import pysam
    index = tmp_path / 'elsewhere.fai'
    pysam.faidx(str(small_fasta), '--fai-idx', str(index))
    assert refs.extract_sequences(small_fasta, ['p1', 'p3', 'absent'], index) == \
        refs.extract_sequences(small_fasta, ['p1', 'p3', 'absent'])


# --- with mash and minimap2 ---

needs_tools = pytest.mark.skipif(not (mash_available() and minimap2_available()), reason="needs mash and minimap2")


@pytest.fixture
def db(tmp_path):
    from tests.scenarios.common import Workspace, mash_sketch
    ws = Workspace(tmp_path)
    sequences, complete, fragment, assembly = scenario.build()
    db_fa = ws.write_fasta('plasmids.fa', sequences.items())
    return {'ws': ws, 'fasta': db_fa, 'sketch': mash_sketch(ws, db_fa, 'plasmids'),
            'complete': ws.write_fasta('complete.fa', complete), 'fragment': ws.write_fasta('fragment.fa', fragment)}


@needs_tools
def test_cli_outputs(db, tmp_path):
    out = tmp_path / 'out'
    code = refs.main(['--plasmid-reconstruction', str(db['complete']), '--db-sketch', str(db['sketch']),
                      '--db-fasta', str(db['fasta']), '-o', str(out), '-p', 's1', '--top', '2', '--log-level', 'error'])
    assert code == 0
    rows = read_tsv(out / 's1_candidates.tsv')
    assert list(rows[0]) == refs.CANDIDATE_FIELDNAMES
    assert [r['ref_id'] for r in rows[:3]] == ['pA', 'pA_variant', 'pA_core']
    assert [r['selected'] for r in rows[:3]] == ['True', 'True', 'False']
    written = sorted(p.name for p in (out / 's1_references').iterdir())
    assert written == ['1_pA.fa', '2_pA_variant.fa']
    assert (out / 's1_references' / '1_pA.fa').read_text().startswith('>pA\n')


@needs_tools
def test_index_gives_same_result(db, tmp_path):
    import pysam
    index = tmp_path / 'plasmids.fa.fai'
    args = ['--plasmid-reconstruction', str(db['fragment']), '--db-sketch', str(db['sketch']), '--db-fasta', str(db['fasta']),
            '--log-level', 'error']
    assert refs.main([*args, '-o', str(tmp_path / 'plain')]) == 0
    pysam.faidx(str(db['fasta']), '--fai-idx', str(index))
    assert refs.main([*args, '--db-fasta-index', str(index), '-o', str(tmp_path / 'indexed')]) == 0
    assert read_tsv(tmp_path / 'plain' / 'candidates.tsv') == read_tsv(tmp_path / 'indexed' / 'candidates.tsv')


@needs_tools
def test_no_shared_kmers(db, tmp_path):
    import random
    from tests.scenarios.common import random_seq
    unrelated = db['ws'].write_fasta('unrelated.fa', [('x', random_seq(random.Random(5), 20000))])
    out = tmp_path / 'none'
    assert refs.main(['--plasmid-reconstruction', str(unrelated), '--db-sketch', str(db['sketch']),
                      '--db-fasta', str(db['fasta']), '-o', str(out), '--log-level', 'error']) == 0
    assert read_tsv(out / 'candidates.tsv') == []
    assert list((out / 'references').iterdir()) == []


def test_missing_mash_is_an_error(monkeypatch, tmp_path):
    monkeypatch.setenv('PATH', str(tmp_path))   # nothing on PATH
    with pytest.raises(refs.RefSearchError):
        refs.sketch_info(tmp_path / 'x.msh')


@needs_tools
def test_plasmid_surgery_with_database(db, tmp_path):
    # End to end: find references in the database, then assess against the best two
    from plasmid_surgery.__main__ import main as plasmid_surgery
    out = tmp_path / 'e2e'
    code = plasmid_surgery(['--plasmid-reconstruction', str(db['complete']), '--db-sketch', str(db['sketch']),
                            '--db-fasta', str(db['fasta']), '--top', '2', '-o', str(out), '--log-level', 'error'])
    assert code == 0
    assert sorted(p.name for p in out.iterdir()) == ['assessment.tsv', 'candidates.tsv', 'pA', 'pA_variant']
    rows = read_tsv(out / 'assessment.tsv')
    assert [(r['ref_id'], r['rank']) for r in rows] == [('pA', '1'), ('pA_variant', '2')]
    assert [r['ref_id'] for r in read_tsv(out / 'candidates.tsv') if r['selected'] == 'True'] == ['pA', 'pA_variant']
    assert (out / 'pA' / 'reference.fa').read_text().startswith('>pA')
