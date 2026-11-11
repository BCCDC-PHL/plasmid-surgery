"""
All-vs-all comparison of reconstructions (plasmid_surgery.compare) and its place in the cohort report.
"""
import csv
import random
import shutil

import pytest

from plasmid_surgery import compare, cohort
from tests.simulate import random_seq
from tests.scenarios.common import mutate

needs_minimap2 = pytest.mark.skipif(shutil.which('minimap2') is None, reason="minimap2 not on PATH")


def member(name, *seqs):
    return {'name': name, 'is_reference': False, 'contigs': [{'seq_id': f'{name}_{i}', 'seq': s} for i, s in enumerate(seqs)]}


def test_average_linkage_order_groups_close_members():
    d = [[0, .9, .1, .9], [.9, 0, .9, .2], [.1, .9, 0, .9], [.9, .2, .9, 0]]
    order = compare.average_linkage_order(d)
    assert sorted(order) == [0, 1, 2, 3]
    pos = {m: i for i, m in enumerate(order)}
    assert abs(pos[0] - pos[2]) == 1 and abs(pos[1] - pos[3]) == 1


@needs_minimap2
def test_all_vs_all_containment(tmp_path):
    rng = random.Random(3)
    a = random_seq(rng, 60000)
    members = [
        member('complete', a),
        member('fragments', a[2000:20000], a[30000:45000]),     # two contigs, both inside 'complete'
        member('variant', mutate(rng, a, 300)),                   # ~0.5% divergent: still contained
        member('other', random_seq(rng, 40000)),
    ]
    r = compare.all_vs_all(members, tmp_path, threads=2)
    c = r['containment']
    assert c[1][0] == pytest.approx(1.0, abs=0.01)                  # fragments in complete
    assert c[0][1] == pytest.approx(33000 / 60000, abs=0.01)        # complete in fragments
    assert c[2][0] > 0.98 and c[0][2] > 0.98
    assert 0.99 < r['identity'][0][2] < 0.997
    assert c[3][0] == 0 and c[0][3] == 0 and r['identity'][0][3] is None
    pos = {m: i for i, m in enumerate(r['order'])}
    assert pos[3] in (0, 3)                                          # the unrelated one at an end
    rows = compare.pair_rows(r)
    assert len(rows) == 6 and {(x['a'], x['b']) for x in rows} >= {('complete', 'fragments')}


@pytest.fixture(scope='module')
def runs(tmp_path_factory):
    from tests.scenarios.cohort_example import build_runs
    return build_runs(tmp_path_factory.mktemp('cohort_compare'))


@needs_minimap2
def test_cohort_report_compares_reconstructions(runs, tmp_path):
    page = tmp_path / 'cohort.html'
    cohort.write_cohort(runs, page, threads=2)
    html = page.read_text()
    assert 'Reconstructions compared' in html and 'cohort_pairs.tsv' in html
    with open(tmp_path / 'cohort_pairs.tsv') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    n = len(runs) + 1                                                # the samples and the reference
    assert len(rows) == n * (n - 1) // 2
    assert all(0 <= float(r['a_in_b']) <= 1 for r in rows)
    # Every sample's reconstruction is (mostly) contained in the shared reference
    ref_id = cohort.gather(runs)['ref_id']
    in_ref = [float(r['b_in_a']) for r in rows if r['a'] == ref_id]
    assert len(in_ref) == len(runs) and min(in_ref) > 0.5
    assert '<h3>SNP distances</h3>' in html
    # The example's samples differ from the reference by SNPs and indels, not by lineage: few SNPs per pair
    assert all(r['snp_distance_reported'] in ('True', 'False') and int(r['snps']) >= 0 for r in rows)
    s1 = next(r for r in rows if {r['a'], r['b']} == {ref_id, 'S1_complete'})
    assert s1['snp_distance_reported'] == 'True' and int(s1['compared_bases']) > 20000


@needs_minimap2
def test_cohort_without_reconstruction_fa(runs, tmp_path):
    copies = []
    for i, d in enumerate(runs):
        dest = tmp_path / f'run{i}'
        shutil.copytree(d, dest)
        if i == 0:
            next(dest.glob('*reconstruction.fa')).unlink()
        copies.append(dest)
    page = tmp_path / 'c.html'
    cohort.write_cohort(copies, page)
    html = page.read_text()
    assert 'Not included (no <code>reconstruction.fa</code>' in html
    cohort.write_cohort(copies, page, compare_reconstructions_=False)
    assert 'Reconstructions compared' not in page.read_text()


def test_covered_blocks_split_at_large_indels():
    # query 0-1100 (+), target 0-1000: 500M, 100I (a novel insertion), 500M
    f = ['q', '1100', '0', '1100', '+', 't', '1000', '0', '1000', '1000', '1100', '60', 'cg:Z:500M100I500M']
    assert compare.covered_blocks(f) == ([(0, 500), (600, 1100)], [(0, 500), (500, 1000)])
    # On the reverse strand, query coordinates run backwards from the end; small indels don't split
    f = ['q', '1000', '0', '1000', '-', 't', '1000', '0', '1000', '990', '1000', '60', 'cg:Z:400M5I200M200D395M']
    q, t = compare.covered_blocks(f)
    assert t == [(0, 600), (800, 1195)] and q == [(395, 1000), (0, 395)]


# --- SNP distances ---

def substitute(seq, positions):
    s = list(seq)
    for p in positions:
        s[p] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[s[p]]
    return ''.join(s)


def snp_row(members, tmp_path, i=0, j=1):
    r = compare.all_vs_all(members, tmp_path)
    row = next(x for x in compare.pair_rows(r) if x['a'] == members[i]['name'] and x['b'] == members[j]['name'])
    return row


@pytest.fixture(scope='module')
def base():
    return random_seq(random.Random(21), 50000)


@needs_minimap2
@pytest.mark.parametrize('reverse', [False, True])
def test_scattered_snps(tmp_path, base, reverse):
    variant = substitute(base, range(2000, 48000, 2300))             # 20 SNPs, 2.3 kb apart
    if reverse:
        from tests.simulate import revcomp
        variant = revcomp(variant)
    row = snp_row([member('a', base), member('b', variant)], tmp_path)
    assert (row['snps'], row['snps_in_clusters'], row['indels']) == (20, 0, 0)
    assert 49000 <= row['compared_bases'] <= 50000 and row['snp_distance_reported']
    assert row['snps_per_100kb'] == pytest.approx(40, abs=1)


@needs_minimap2
def test_clustered_snps_and_indels(tmp_path, base):
    variant = substitute(base, list(range(10000, 10900, 100)) + [20000, 30000])   # 9 in 900 bp + 2 apart
    variant = variant[:40000] + variant[40003:45000] + 'ACGT' + variant[45000:]    # a 3 bp deletion, a 4 bp insertion
    row = snp_row([member('a', base), member('b', variant)], tmp_path)
    assert (row['snps'], row['snps_in_clusters'], row['snps_outside_clusters'], row['indels']) == (11, 9, 2, 2)


@needs_minimap2
def test_contig_ends_trimmed(tmp_path, base):
    # B in two contigs; a SNP 30 bp from a contig end doesn't count, one well inside does
    variant = substitute(base, [24970, 35000])
    row = snp_row([member('a', base), member('b', variant[:25000], variant[25000:])], tmp_path)
    assert row['snps'] == 1


@needs_minimap2
def test_repeat_copies(tmp_path):
    rng = random.Random(5)
    insertion = random_seq(rng, 1500)
    diverged = substitute(insertion, range(100, 1500, 300))          # a second copy, 5 SNPs apart
    x, y, z = random_seq(rng, 20000), random_seq(rng, 20000), random_seq(rng, 20000)
    a = x + insertion + y + diverged + z
    # Whole plasmids: each copy is compared with its own counterpart
    row = snp_row([member('a', a), member('b', a)], tmp_path)
    assert row['snps'] == 0 and row['compared_bases'] >= len(a) - 400
    # B with the diverged copy as a contig of its own (as short-read assemblers do with repeats): it aligns to
    # both of A's copies, so it's masked rather than compared with the wrong one
    row = snp_row([member('a', a), member('b', x + insertion + y, diverged, z)], tmp_path / 'collapsed')
    assert row['snps'] == 0
    assert row['compared_bases'] <= len(a) - 1500
@needs_minimap2
def test_too_little_compared(tmp_path, base):
    row = snp_row([member('a', base), member('b', base[:20000] + random_seq(random.Random(9), 30000))], tmp_path)
    assert row['compared_bases'] < 25000 and not row['snp_distance_reported']


def test_clustered():
    assert compare.clustered([100, 500, 900, 5000, 5100]) == {0, 1, 2}
    assert compare.clustered([100, 1200, 2300]) == set()
