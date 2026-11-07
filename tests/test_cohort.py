import re
import shutil

import pytest

from plasmid_surgery import cohort
from tests.scenarios.common import minimap2_available


def test_runs_of():
    assert cohort.runs_of(bytearray([0, 0, 1, 1, 1, 0])) == [(0, 2, 0), (2, 5, 1), (5, 6, 0)]


def ref(gaps, sources=None):
    return {'gaps': gaps, 'contig_sources': sources}


def test_categories_without_assembly():
    cats = cohort.categories(ref([{'start': 10, 'end': 20, 'type': 'unique'}, {'start': 30, 'end': 40, 'type': 'repeat'}]), 50)
    assert cohort.runs_of(cats) == [(0, 10, cohort.COVERED), (10, 20, cohort.ABSENT), (20, 30, cohort.COVERED),
                                    (30, 40, cohort.REPEAT_GAP), (40, 50, cohort.COVERED)]


def test_categories_with_assembly_missed_takes_precedence():
    sources = [{'classification': 'shared', 'ref_start': '12', 'ref_end': '30'},
               {'classification': 'missed', 'ref_start': '15', 'ref_end': '18'},
               {'classification': 'in_reconstruction', 'ref_start': '0', 'ref_end': '10'}]
    cats = cohort.categories(ref([{'start': 10, 'end': 20, 'type': 'unique'}], sources), 50)
    assert cohort.runs_of(cats)[:5] == [(0, 10, cohort.COVERED), (10, 12, cohort.ABSENT), (12, 15, cohort.OTHER_ASSEMBLY),
                                        (15, 18, cohort.MISSED), (18, 20, cohort.OTHER_ASSEMBLY)]


def test_categories_wrapping_contig():
    sources = [{'classification': 'missed', 'ref_start': '45', 'ref_end': '5'}]   # spans the origin
    cats = cohort.categories(ref([{'start': 0, 'end': 10, 'type': 'unique'}, {'start': 40, 'end': 50, 'type': 'unique'}],
                                 sources), 50)
    assert cats[0] == cats[47] == cohort.MISSED and cats[7] == cats[42] == cohort.ABSENT


def test_no_runs_is_an_error(tmp_path):
    assert cohort.main([str(tmp_path), '-o', str(tmp_path / 'c.html'), '--log-level', 'error']) == 1


needs_minimap2 = pytest.mark.skipif(not minimap2_available(), reason="minimap2 not on PATH")


@pytest.fixture(scope='module')
def runs(tmp_path_factory):
    from tests.scenarios.cohort_example import build_runs
    return build_runs(tmp_path_factory.mktemp('cohort'))


def sample(c, name):
    return next(s for s in c['samples'] if s['sample_id'] == name)


@needs_minimap2
def test_gather_and_categories(runs):
    c = cohort.gather(runs)
    assert c['ref_id'] == 'pREF1' and len(c['samples']) == 5 and not c['skipped']
    length = c['ref_length']
    cats = {s['sample_id']: cohort.categories(s['ref'], length) for s in c['samples']}
    assert set(cats['S1_complete']) == {cohort.COVERED, cohort.REPEAT_GAP}
    # S2's third segment is in a missed contig (binning); S3's u2 second half is nowhere
    from tests.scenarios.cohort_example import R, U
    u3_start = 2 * U + 2 * R
    assert cats['S2_misbinned'][u3_start + 100] == cohort.MISSED
    assert cats['S3_incomplete'][U + R + 6000] == cohort.ABSENT
    # S4's deletion (u1 3000-4500) is uncovered and in no assembly contig
    assert set(cats['S4_indels'][3000:4500]) == {cohort.ABSENT}


@needs_minimap2
def test_cohort_page(runs, tmp_path):
    out = cohort.write_cohort(runs, tmp_path / 'cohort.html')
    page = out.read_text()
    assert len(re.findall(r'<div class="strip-row" data-order=', page)) == 5
    assert 'geneNovel' in page.split('Genes not on the reference')[1]
    assert page.count('class="binning"') >= 2          # S2's misbinned gene (plus the legend)
    assert 'href="results/S1_complete/report.html"' not in page   # links are relative to the page
    assert re.search(r'href="[^"]*S1_complete/report.html"', page)
    assert '<script>' in page and 'http' not in page.split('<script>')[1]


@needs_minimap2
def test_ref_choice_and_skipping(runs, tmp_path):
    with pytest.raises(cohort.CohortError):
        cohort.gather(runs, 'no_such_ref')
    # A run assessed against a different-length copy of the reference is skipped
    odd = tmp_path / 'odd'
    shutil.copytree(runs[0], odd)
    for p in odd.rglob('assessment.tsv'):
        p.write_text(p.read_text().replace('\t29400\t', '\t29401\t'))
    c = cohort.gather([*runs, odd])
    assert len(c['samples']) == 5 and any('29,401' in reason for _, reason in c['skipped'])
    assert cohort.main([*map(str, runs), '--ref', 'pREF1', '-o', str(tmp_path / 'c.html'), '--log-level', 'error']) == 0


def write_assessment(d, rows, prefix=''):
    d.mkdir(parents=True)
    lines = ['sample_id\tref_id\trank\tscore\tref_length\tref_fraction_covered']
    lines += ['\t'.join(map(str, (d.name, *r))) for r in rows]
    (d / f'{prefix}assessment.tsv').write_text('\n'.join(lines) + '\n')


@pytest.fixture
def assessments(tmp_path):
    #                         ref_id, rank, score, ref_length, ref_fraction_covered
    write_assessment(tmp_path / 's1', [('A', 1, 0.9, 5000, 0.95), ('B', 2, 0.8, 6000, 0.9)])
    write_assessment(tmp_path / 's2', [('B', 1, 0.85, 6000, 0.9), ('A', 2, 0.7, 5000, 0.8), ('C', 3, 0.5, 7000, 0.6)],
                     prefix='s2_')
    write_assessment(tmp_path / 's3', [('A', 1, 0.88, 5000, 0.93), ('C', 4, 0.1, 7000, 0.2)])
    (tmp_path / 'not_a_run').mkdir()
    return tmp_path


def test_reference_table(assessments):
    dirs = sorted(assessments.iterdir())
    table, runs = cohort.reference_table(dirs)
    assert runs == 3
    assert [r['ref_id'] for r in table] == ['A', 'B', 'C']
    a, b, c = table
    assert (a['runs_assessed'], a['rank_1'], a['in_top'], a['ref_length']) == (3, 2, 3, '5000')
    assert a['mean_score'] == round((0.9 + 0.7 + 0.88) / 3, 4)
    assert a['mean_ref_fraction_covered'] == round((0.95 + 0.8 + 0.93) / 3, 4)
    assert b['mean_score'] == round((0.8 + 0.85) / 3, 4)              # unassessed runs count as 0
    assert b['mean_ref_fraction_covered'] == 0.9                     # only over runs that assessed it
    assert (c['runs_assessed'], c['in_top']) == (2, 1)               # rank 4 is outside the top 3
    assert cohort.reference_table(dirs, top=4)[0][2]['in_top'] == 2


def test_list_refs_cli(assessments, capsys):
    dirs = [str(p) for p in sorted(assessments.iterdir())]
    assert cohort.main(dirs + ['--list-refs', '--top', '2', '--log-level', 'error']) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == 'ref_id\tref_length\truns_assessed\trank_1\tin_top_2\tmean_score\tmean_ref_fraction_covered'
    assert lines[1].split('\t')[:5] == ['A', '5000', '3/3', '2', '3']
    assert lines[3].split('\t')[:5] == ['C', '7000', '2/3', '0', '0']
    assert not (assessments / 'cohort.html').exists()


def test_list_refs_no_runs(tmp_path):
    assert cohort.main([str(tmp_path), '--list-refs', '--log-level', 'error']) == 1


@needs_minimap2
def test_reference_table_on_real_runs(runs):
    table, n = cohort.reference_table(runs)
    assert n == len(runs) and table[0]['ref_id'] == cohort.gather(runs)['ref_id']
