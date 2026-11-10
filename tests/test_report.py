import re
import xml.etree.ElementTree as ET

import pytest

from plasmid_surgery import report
from tests.scenarios.common import minimap2_available

SVG_NS = '{http://www.w3.org/2000/svg}'


def test_short_labels():
    assert report.short_labels(['SAMN1_12', 'SAMN1_7', 'SAMN1_103']) == {'SAMN1_12': '12', 'SAMN1_7': '7', 'SAMN1_103': '103'}
    assert report.short_labels(['contig_1', 'contig_10']) == {'contig_1': '1', 'contig_10': '10'}
    assert report.short_labels(['a1', 'b2']) == {'a1': 'a1', 'b2': 'b2'}      # no shared prefix ending in a separator
    assert report.short_labels(['only']) == {'only': 'only'}


def test_unwrap_and_wrap():
    assert report.unwrap(90, 10, 100) == (90, 110)
    assert report.unwrap(10, 90, 100) == (10, 90)
    assert report.wrap(90, 110, 100) == [(90, 100), (0, 10)]


def test_assign_lanes_avoids_overlap_including_wrapped():
    lanes = report.assign_lanes([(0, 50), (40, 60), (60, 100), (90, 110)], 100)
    assert lanes == [0, 1, 0, 1]   # (90, 110) wraps onto (0, 10), so it can't share lane 0


def test_missing_run_is_an_error(tmp_path):
    assert report.main([str(tmp_path), '--log-level', 'error']) == 1


needs_minimap2 = pytest.mark.skipif(not minimap2_available(), reason="minimap2 not on PATH")


@pytest.fixture(scope='module')
def example(tmp_path_factory):
    from tests.scenarios.report_example import build_run
    out = build_run(tmp_path_factory.mktemp('report'))
    path = report.write_report(out)
    return {'out': out, 'html': path.read_text(), 'run': report.load_run(out)}


@needs_minimap2
def test_load_run(example):
    run = example['run']
    assert [r['row']['ref_id'] for r in run['references']] == ['pREF1', 'pREF2_related']   # by rank
    assert run['prefix'] == '' and run['candidates'] is None
    ref = run['references'][0]
    assert ref['contig_sources'] is not None and ref['ref_features'] and ref['features']


@needs_minimap2
def test_map_svgs_are_valid_and_complete(example):
    for d in ('pREF1', 'pREF2_related'):
        root = ET.parse(example['out'] / d / 'map.svg').getroot()
        assert root.tag == f'{SVG_NS}svg'
        titles = [t.text for t in root.iter(f'{SVG_NS}title')]
        assert any(t.startswith('u2 (') for t in titles)              # a contig's tooltip
        assert any(' gap ' in t for t in titles)                       # a gap's tooltip
    titles = [t.text for t in ET.parse(example['out'] / 'pREF1' / 'map.svg').getroot().iter(f'{SVG_NS}title')]
    assert any(t.startswith('u3 (') and 'missed' in t for t in titles)        # the misbinned contig
    assert any(t.startswith('geneMisbinned') and 'absent' in t for t in titles)  # a ref feature's status
    assert any(t.startswith('insertion of') for t in titles)


@needs_minimap2
def test_report_page(example):
    page = example['html']
    # Panels in rank order, each with both scales
    assert page.index('id="ref-pREF1"') < page.index('id="ref-pREF2_related"')
    assert page.count('class="scroll map-shared"') == page.count('class="scroll map-fit"') == 2
    # A contig has the same colour in every panel
    fills = re.findall(r'<polygon [^>]*fill="(#[0-9a-f]{6})"[^>]*class="contig" data-contig="u2"', page)
    assert len(fills) >= 4 and len(set(fills)) == 1
    # The matrix has a row per reconstruction contig
    assert len(re.findall(r'<tr data-contig="', page)) == 3
    assert 'none annotated on this reference' in page          # pREF2 has no annotations
    assert 'Not on this reference: geneInInsertion' in page
    assert '<script>' in page and 'http' not in page.split('<script>')[1]   # no external resources


@needs_minimap2
def test_report_flag(tmp_path, example):
    from plasmid_surgery.__main__ import main
    out = tmp_path / 'flagged'
    ref = example['out'] / 'pREF1' / 'reference.fa'
    recon = example['out'].parent / 'plasmid_1.fasta'
    assert main(['--plasmid-reconstruction', str(recon), '--ref-plasmid', str(ref), '-o', str(out), '-p', 's1',
                 '--report', '--log-level', 'error']) == 0
    assert (out / 's1_report.html').exists() and (out / 'pREF1' / 's1_map.svg').exists()


def test_reconstruction_view_chimera(tmp_path):
    """A reconstruction that is 40 kb of one plasmid and 55 kb of another, against the first only:
    the view says how much each reference explains and marks the rest, and genes in it."""
    import random
    import shutil
    if shutil.which('minimap2') is None:
        pytest.skip("minimap2 not on PATH")
    from plasmid_surgery.__main__ import main
    from tests.simulate import write_fasta, random_seq
    rng = random.Random(7)
    a, b = random_seq(rng, 60000), random_seq(rng, 70000)
    write_fasta(tmp_path / 'recon.fa', [('hybrid_1', a[5000:45000] + b[10000:65000])])
    write_fasta(tmp_path / 'refs.fa', [('refA', a), ('unrelated', random_seq(rng, 30000))])
    (tmp_path / 'feats.bed').write_text("hybrid_1\t2000\t2900\tgeneA\t0\t+\nhybrid_1\t60000\t61200\tgeneB\t0\t+\n")
    assert main(['--plasmid-reconstruction', str(tmp_path / 'recon.fa'), '--ref-plasmid', str(tmp_path / 'refs.fa'),
                 '--features', str(tmp_path / 'feats.bed'), '--report', '-o', str(tmp_path / 'out'),
                 '--log-level', 'warning']) == 0
    html = (tmp_path / 'out' / 'report.html').read_text()
    assert '42% of the reconstruction (40,000 of 95,000 bp) aligns to at least one reference' in html
    assert 'The other 55,000 bp (1 stretch) align to none of them' in html
    assert 'hybrid_1:40,001–95,000 (55,000 bp): aligns to none of the references' in html
    assert 'Explaining none of it: unrelated' in html
    assert 'geneB' in html and 'in a part no reference explains' in html
    summary = {r['alignment_class']: r for r in report.read_tsv(tmp_path / 'out' / 'refA' / 'summary.tsv')}
    assert (summary['placement']['contig_start'], summary['placement']['contig_end']) == ('0', '40000')
