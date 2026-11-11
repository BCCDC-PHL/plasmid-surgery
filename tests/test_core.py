"""
Core SNPs on the shared reference (plasmid_surgery.core).
"""
import random
import shutil

import pytest

from plasmid_surgery import core, cohort
from tests.simulate import random_seq

needs_minimap2 = pytest.mark.skipif(shutil.which('minimap2') is None, reason="minimap2 not on PATH")


def substitute(seq, positions):
    s = list(seq)
    for p in positions:
        s[p] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[s[p]]
    return ''.join(s)


def sample(name, *seqs):
    return {'name': name, 'contigs': [{'seq_id': f'{name}_{i}', 'seq': s} for i, s in enumerate(seqs)]}


@pytest.fixture(scope='module')
def cohort_seqs():
    rng = random.Random(13)
    insertion = random_seq(rng, 1500)
    ref = random_seq(rng, 15000) + insertion + random_seq(rng, 25000) + insertion + random_seq(rng, 18000)   # 61 kb
    shared = [5000, 40000]
    s1 = substitute(ref, shared + [50000])
    s2 = substitute(ref, shared + [8000, 9000, 33000, 45000, 55000])
    s3 = substitute(ref[:30000], [5000])                              # a fragment: the first half only
    s4 = ref[:20000] + ref[20500:]                                    # a 500 bp deletion
    return {'ref': ref, 'samples': [sample('S1', s1), sample('S2', s2), sample('S3', s3), sample('S4', s4)]}


@needs_minimap2
def test_strict_core(tmp_path, cohort_seqs):
    ref = {'seq_id': 'ref', 'seq': cohort_seqs['ref']}
    r = core.core_snps(ref, cohort_seqs['samples'], tmp_path, threads=2)
    n = len(cohort_seqs['ref'])
    core_set = set(r['core_positions'])
    assert r['names'] == ['ref', 'S1', 'S2', 'S3', 'S4']
    # Limited by the fragment S3 (first 30 kb, ends trimmed), without the reference's repeat copy at 15-16.5 kb,
    # nor S4's deletion
    assert max(core_set) < 30000 and 15500 not in core_set and 20200 not in core_set
    assert 25000 < len(core_set) + 1500 + 500 < 30000
    pair = {(p['a'], p['b']): p for p in r['pairs'].values()}
    assert pair[('S1', 'S2')]['snps'] == 2                            # 8,000 and 9,000 (S2 only); 5,000 shared
    assert pair[('ref', 'S1')]['snps'] == 1 and pair[('ref', 'S3')]['snps'] == 1
    assert pair[('S1', 'S4')]['snps'] == 1
    assert all(p['compared_positions'] == len(core_set) for p in r['pairs'].values())
    rows = {x['sample']: x for x in r['samples']}
    assert rows['S3']['called_fraction'] < 0.5 < rows['S1']['called_fraction']
    assert rows['S3']['core_without_it'] > 50000                      # the fragment is what shrinks the core
    assert rows['S1']['core_without_it'] == len(core_set)


@needs_minimap2
def test_soft_core(tmp_path, cohort_seqs):
    ref = {'seq_id': 'ref', 'seq': cohort_seqs['ref']}
    r = core.core_snps(ref, cohort_seqs['samples'], tmp_path, threads=2, min_fraction=0.75)
    assert max(r['core_positions']) > 50000                           # beyond the fragment: 3 of 4 samples call it
    pair = {(p['a'], p['b']): p for p in r['pairs'].values()}
    assert pair[('S1', 'S2')]['snps'] == 6                            # 8k, 9k, 33k, 45k, 50k (S1), 55k
    assert pair[('S1', 'S3')]['compared_positions'] < pair[('S1', 'S2')]['compared_positions']
    assert all(x['core_without_it'] == 'NA' for x in r['samples'])


@needs_minimap2
def test_outputs(tmp_path, cohort_seqs):
    r = core.core_snps({'seq_id': 'ref', 'seq': cohort_seqs['ref']}, cohort_seqs['samples'], tmp_path)
    paths = [tmp_path / x for x in ('core.fasta', 'pairs.tsv', 'samples.tsv')]
    core.write_outputs(r, *paths)
    seqs = paths[0].read_text().split('>')[1:]
    assert len(seqs) == 5 and len({len(''.join(s.splitlines()[1:])) for s in seqs}) == 1
    assert len(paths[1].read_text().splitlines()) == 1 + 10


def test_differing():
    assert core.differing(b'ACGT' * 40, b'ACGT' * 20 + b'TCGT' + b'AC-T' + b'ACGT' * 18) == {80}


@needs_minimap2
def test_cohort_core_section(tmp_path):
    from tests.scenarios.cohort_example import build_runs
    (tmp_path / 'runs').mkdir()
    runs = build_runs(tmp_path / 'runs')
    page = tmp_path / 'cohort.html'
    cohort.write_cohort(runs, page, threads=2)
    html = page.read_text()
    assert 'Core SNPs on the reference' in html
    for name in ('cohort_core.fasta', 'cohort_core_pairs.tsv', 'cohort_core_samples.tsv'):
        assert (tmp_path / name).exists()
