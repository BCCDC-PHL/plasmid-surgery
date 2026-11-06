"""The choice of reference: a fragment of a larger plasmid can look complete against a smaller one."""

from __future__ import annotations

import random

from tests.scenarios.common import Outcome, Scenario, outputs, plasmid_map, random_seq, read_tsv, workspace

LARGE = 40000
FRAGMENT = (10000, 25000)   # the part of the large plasmid that's in the reconstruction
EXTRA = 1000                # sequence in the small ref that isn't in the large plasmid


def build():
    rng = random.Random(83)
    large = random_seq(rng, LARGE)
    fragment = large[FRAGMENT[0]:FRAGMENT[1]]
    small = fragment + random_seq(rng, EXTRA)
    mid = (FRAGMENT[1] - FRAGMENT[0]) // 2
    reconstruction = [('r1', fragment[:mid]), ('r2', fragment[mid:])]
    elsewhere = [('e1', large[FRAGMENT[1]:]), ('e2', large[:FRAGMENT[0]])]   # the rest, in another bin
    return large, small, reconstruction, elsewhere


def run() -> Outcome:
    ws, tmp = workspace()
    large, small, reconstruction, elsewhere = build()
    with tmp:
        large_fa = ws.write_fasta('large_ref.fa', [('large_ref', large)])
        small_fa = ws.write_fasta('small_ref.fa', [('small_ref', small)])
        recon_fa = ws.write_fasta('plasmid_1.fasta', reconstruction)
        other_fa = ws.write_fasta('plasmid_2.fasta', elsewhere)
        assembly_fa = ws.write_fasta('assembly.fasta', reconstruction + elsewhere)
        out = ws.run('--plasmid-reconstruction', recon_fa, '--assembly', assembly_fa, '--bins', recon_fa, other_fa,
                     '--ref-plasmid', small_fa, large_fa, '-o', ws.path('out'))
        vs_small, vs_large = outputs(out, 'small_ref'), outputs(out, 'large_ref')
        ranked = read_tsv(out / 'assessment.tsv')
        without_assembly = read_tsv(ws.run('--plasmid-reconstruction', recon_fa, '--ref-plasmid', small_fa, large_fa,
                                           '-o', ws.path('out_without_assembly')) / 'assessment.tsv')

    s, l = vs_small['assessment'], vs_large['assessment']
    rows = [
        {'ref': f'small ({len(small):,} bp)', 'ref_fraction_covered': s['ref_fraction_covered'],
         'unique_gap_length': s['unique_gap_length'], 'contig_fraction_aligned': s['contig_fraction_aligned'],
         'missed_ref_bases': s['missed_ref_bases']},
        {'ref': f'large ({len(large):,} bp)', 'ref_fraction_covered': l['ref_fraction_covered'],
         'unique_gap_length': l['unique_gap_length'], 'contig_fraction_aligned': l['contig_fraction_aligned'],
         'missed_ref_bases': l['missed_ref_bases']},
    ]
    facts = {
        'large': LARGE, 'small': len(small), 'fragment_length': FRAGMENT[1] - FRAGMENT[0], 'extra': EXTRA,
        'small_covered': float(s['ref_fraction_covered']),
        'large_covered': float(l['ref_fraction_covered']),
        'small_aligned': float(s['contig_fraction_aligned']),
        'small_missed': int(s['missed_ref_bases']),
        'large_missed': int(l['missed_ref_bases']),
        'large_missed_contigs': int(l['num_missed_contigs']),
        'ranked': ', '.join(f"`{r['ref_id']}` (score {r['score']})" for r in ranked),
        'large_rank': int(l['rank']), 'small_rank': int(s['rank']),
        'best_without_assembly': without_assembly[0]['ref_id'],
    }
    columns = ['ref', 'ref_fraction_covered', 'unique_gap_length', 'contig_fraction_aligned', 'missed_ref_bases']
    data = {'tables': {'comparison': (columns, rows)},
            'maps': {'small': plasmid_map(vs_small, len(small)), 'large': plasmid_map(vs_large, LARGE)}}
    return Outcome(ws.commands, facts, data)


SCENARIO = Scenario(
    slug='reference-choice',
    title='The choice of reference',
    summary="Everything is measured against the reference given: a fragment of a larger plasmid can look "
            "nearly complete against a smaller plasmid that shares it.",
    text="""\
plasmid-surgery measures a reconstruction against the reference it's given. It can't tell whether
that's the right reference.

Here the reconstruction is two contigs, {fragment_length:,} bp in all, from a {large:,} bp plasmid;
the rest of that plasmid is in another bin. It's assessed against two references: the large plasmid,
and a smaller {small:,} bp plasmid that shares the same {fragment_length:,} bp (plus {extra:,} bp of
its own), in one run (`--ref-plasmid small_ref.fa large_ref.fa`), with the whole assembly and both
bins.

{comparison}

[[figure]]

Against the small reference, the reconstruction looks nearly complete: {small_covered:.1%} covered,
every contig aligned ({small_aligned:.0%}), and nothing missed elsewhere in the assembly
({small_missed} bp), because the rest of the large plasmid doesn't align to the small one. Against the
large reference it covers {large_covered:.1%}, and whole-assembly mode finds the rest
({large_missed:,} bp, in {large_missed_contigs} `missed` contigs).

Given several references, the top-level `assessment.tsv` ranks them by score: the fraction of the
reconstruction aligned × the fraction of the reference covered, by the whole assembly when there is
one (as in [Finding references](find-refs.md)). Here: {ranked}. Without the assembly, `{best_without_assembly}`
ranks first.

So a good result means "consistent with this reference", not "complete". Where the right reference is
uncertain, it's worth assessing against more than one candidate: the larger one reveals what the
smaller one can't. A reference much smaller than the plasmid's expected size is a warning sign.
""",
    run=run,
    test='tests/test_scenarios.py::test_reference_choice',
)
