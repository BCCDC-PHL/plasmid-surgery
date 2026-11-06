"""Check every claim the Behaviour report makes (see tests/scenarios).

Each scenario runs the real CLI once; these tests assert what its text says, so the report can't
quietly drift from what the code does.
"""

import re
from functools import cache

import pytest

from tests.scenarios import MODULES, SCENARIOS, assembly, features, find_refs, indels, origin, overlaps, reference_choice, repeats
from tests.scenarios.common import mash_available, minimap2_available, render

pytestmark = pytest.mark.skipif(not minimap2_available(), reason="minimap2 not on PATH")


@cache
def outcome(module):
    return module.SCENARIO.run()


@pytest.mark.parametrize("module", MODULES, ids=lambda m: m.SCENARIO.slug)
def test_report_text_renders(module):
    if module is find_refs and not mash_available():
        pytest.skip("mash not on PATH")
    text = render(module.SCENARIO, outcome(module))
    assert not re.search(r"\{[a-z_0-9]+(:[^}]*)?\}", text), "unfilled placeholder"
    assert text.count("[[figure]]") <= 1
    assert module.SCENARIO.test.startswith("tests/test_scenarios.py::")
    assert len({s.slug for s in SCENARIOS}) == len(SCENARIOS)


def test_origin():
    f = outcome(origin).facts
    assert f['linear_segments'] == 2                        # minimap2 alone splits c1 at the origin
    assert f['c1_rows'] == 1 and f['c1_spans_origin']       # plasmid-surgery places it in one piece
    assert f['c1_end'] < f['c1_start']
    assert f['scaffold_is_sample'] and f['scaffold_starts_with_c1_tail']
    assert f['doubled_mapq'] == 0 and f['single_mapq'] == 60


def test_repeats():
    f = outcome(repeats).facts
    assert f['repeat_rows'] == 3 and f['repeat_placements'] == 3 and f['repeat_mapq'] == 0
    assert f['u3_placements'] == 1
    assert f['num_repeat_gaps'] == 2 and f['num_unique_gaps'] == 0
    assert f['covered'] < 1 and f['covered_with_repeats'] == 1.0
    assert f['lowercase'] == f['repeat_gap_length'] and f['scaffold_ns'] == 0
    assert f['scaffold_is_ref']


def test_overlaps():
    f = outcome(overlaps).facts
    assert f['scaffold_is_ref'] and f['scaffold_length'] == f['ref_length'] and f['scaffold_ns'] == 0
    assert not f['d_included'] and f['others_included']
    assert f['covered_multiple'] == f['expected_multiple']
    assert f['scaffold_base'] == f['ref_base'] != f['b_base']   # the first contig by position wins


def test_indels():
    f = outcome(indels).facts
    assert f['num_events'] == 2
    assert (f['event_deletion_start'], f['event_deletion_length']) == (f['deletion_start'], f['deletion_length'])
    assert (f['event_insertion_at'], f['event_insertion_length']) == (f['insertion_at'], f['insertion_length'])
    assert f['num_gaps'] == 1 and (f['gap_start'], f['gap_end']) == (f['deletion_start'], f['deletion_end'])
    assert f['y_fraction_aligned'] == 1.0
    assert f['y_identity'] > 99.5 and f['y_blast_identity'] < 85
    assert f['x_qc_pass']


def test_assembly():
    f = outcome(assembly).facts
    assert f['B_class'] == 'missed' and f['shared_class'] == 'shared' and f['D_class'] == 'shared'
    assert f['D_fraction'] < 0.8 <= 1.0
    assert not f['other_listed']
    assert f['B_depth'] == pytest.approx(1.05) and f['shared_depth'] == pytest.approx(0.5)
    assert f['covered'] == 0.5 and f['assembly_covered'] == 1.0
    assert f['unique_gap_length'] == 2 * f['segment'] and f['assembly_unique_gap_length'] == 0
    assert f['missed_ref_bases'] == f['segment'] and f['num_missed'] == 1 and f['num_shared'] == 2


def test_features():
    f = outcome(features).facts
    assert f['geneA_status'] == 'present'
    assert (f['deleted_status'], f['deleted_assembly_status']) == ('absent', 'absent')
    assert f['tnp_statuses'] == '`present` and `inferred`'
    assert f['tnp_assembly_statuses'] == ['inferred', 'present']
    assert (f['misbinned_status'], f['misbinned_assembly_status'], f['misbinned_found_in']) == ('absent', 'present', 'chromosome')
    assert f['sample_geneA_ref'] == ('1000', '1900')
    assert (f['insertion_gene_status'], f['insertion_gene_class']) == ('not_on_ref', 'in_reconstruction')
    assert (f['sample_misbinned_class'], f['sample_misbinned_bin']) == ('missed', 'chromosome')


def test_reference_choice():
    f = outcome(reference_choice).facts
    assert f['small_covered'] > 0.9 and f['small_aligned'] == 1.0 and f['small_missed'] == 0
    assert f['large_covered'] == pytest.approx(f['fragment_length'] / f['large'])
    assert f['large_missed'] == f['large'] - f['fragment_length']
    assert f['large_rank'] == 1 and f['small_rank'] == 2
    assert f['best_without_assembly'] == 'small_ref'


@pytest.mark.skipif(not mash_available(), reason="mash not on PATH")
def test_find_refs():
    f = outcome(find_refs).facts
    assert f['a_best'] == f['a_nearest'] == 'pA' and f['a_best_score'] >= 0.99
    assert f['a_variant_covered'] < 0.9 and f['a_variant_score'] < f['a_best_score']
    assert f['a_core_aligned'] < 0.7 and f['a_core_score'] < f['a_variant_score']
    assert f['a_unrelated_candidates'] == 0
    # Fragment: against the reconstruction alone, the small relative wins...
    assert f['b_best'] == 'pB_small' and f['b_small_score'] > f['b_large_score']
    assert f['b_large_aligned'] == 1.0 and f['b_large_covered'] < 0.4
    # ...and with the assembly, the large plasmid does
    assert f['b_asm_best'] == 'pB' and f['b_asm_large_covered'] == 1.0
    assert f['b_asm_large_score'] > f['b_asm_small_score']
