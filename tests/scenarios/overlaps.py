"""Overlapping and contained contigs: how the scaffold is assembled from contigs that overlap on the ref."""

from __future__ import annotations

import random

from tests.scenarios.common import Outcome, Scenario, outputs, plasmid_map, random_seq, workspace

REF_LENGTH = 30000
OVERLAP = 127          # like the (k - 1) overlaps between adjacent de Bruijn graph contigs
SNP_AT = 10050         # inside the overlap of a and b
TRUTH = {'a': (0, 10000 + OVERLAP), 'b': (10000, 20000 + OVERLAP), 'c': (20000, 30000 + OVERLAP),
         'd': (12000, 13000)}


def build():
    rng = random.Random(31)
    ref = random_seq(rng, REF_LENGTH)
    doubled = ref + ref
    contigs = {name: doubled[s:e] for name, (s, e) in TRUTH.items()}
    # b has a different base where it overlaps a
    i = SNP_AT - TRUTH['b'][0]
    alt = 'A' if ref[SNP_AT] != 'A' else 'C'
    contigs['b'] = contigs['b'][:i] + alt + contigs['b'][i + 1:]
    return ref, contigs, alt


def run() -> Outcome:
    ws, tmp = workspace()
    ref, contigs, alt = build()
    with tmp:
        ref_fa = ws.write_fasta('ref.fa', [('ref', ref)])
        recon_fa = ws.write_fasta('reconstruction.fa', contigs.items())
        out = outputs(ws.run('--plasmid-reconstruction', recon_fa, '--ref-plasmid', ref_fa, '-o', ws.path('out')))

    included = {r['contig_id']: r['included_in_output'] == 'True' for r in out['summary']}
    scaffold = out['scaffold']
    facts = {
        'overlap': OVERLAP,
        'snp_at': SNP_AT,
        'ref_base': ref[SNP_AT],
        'b_base': alt,
        'scaffold_base': scaffold[SNP_AT],
        'scaffold_is_ref': scaffold == ref,
        'scaffold_length': len(scaffold),
        'ref_length': REF_LENGTH,
        'scaffold_ns': scaffold.count('N'),
        'd_included': included['d'],
        'others_included': all(included[c] for c in 'abc'),
        'covered_multiple': int(out['assessment']['ref_bases_covered_multiple']),
        'expected_multiple': 3 * OVERLAP + (TRUTH['d'][1] - TRUTH['d'][0]),
    }
    tables = {'included': (['contig', 'ref position', 'included_in_output'],
                           [{'contig': c, 'ref position': f"{s:,}–{e:,}", 'included_in_output': included[c]}
                            for c, (s, e) in TRUTH.items()])}
    return Outcome(ws.commands, facts, {'map': plasmid_map(out, REF_LENGTH, TRUTH), 'tables': tables})


SCENARIO = Scenario(
    slug='overlaps',
    title='Overlapping and contained contigs',
    summary="Contigs that overlap on the reference are trimmed so the overlap appears once (the first contig "
            "by position wins), and a contig entirely within another is left out of the scaffold.",
    text="""\
Adjacent contigs from a de Bruijn graph assembly typically overlap by a few dozen to a hundred or so
bases (k − 1). Here contigs a, b and c tile a {ref_length:,} bp reference, each overlapping the next by
{overlap} bp (c wraps the origin and overlaps a), and a short contig d lies entirely within b:

{included}

[[figure]]

The scaffold is built by walking along the reference in order of alignment start. Each contig is
placed from where the previous one ended: if it overlaps, its overlapping start is trimmed (using its
alignment to find the matching contig position, so indels in the overlap are handled). A contig
entirely within the region already covered, like d, is left out of the scaffold
(`included_in_output` = {d_included}), but still counts towards coverage:
`ref_bases_covered_multiple` = {covered_multiple:,} (3 × {overlap} + d's length).

The result is the plasmid with each overlap once: {scaffold_length:,} bp, {scaffold_ns} Ns, identical to
the reference: **{scaffold_is_ref}**.

**Which copy of an overlap wins?** The first contig by position. To show it, b has a different base
({b_base}) from a ({ref_base}) at {snp_at:,}, inside their overlap. The scaffold has a's base:
`{scaffold_base}`. This is simple and predictable, but not necessarily the better-supported base.
Overlaps between adjacent contigs are short and normally agree, so this rarely matters; where it
does, the contigs themselves disagree, which the scaffold can't resolve.
""",
    run=run,
    test='tests/test_scenarios.py::test_overlaps',
)
