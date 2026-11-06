"""Circular references: a contig spanning the reference's origin, and why the reference isn't doubled."""

from __future__ import annotations

import random

from plasmid_surgery import parsers
from tests.scenarios.common import Outcome, Scenario, outputs, plasmid_map, plural, random_seq, revcomp, workspace

REF_LENGTH = 30000
# Where each contig really belongs on the ref; c1 crosses the origin
TRUTH = {'c1': (24000, 36000), 'c2': (6000, 15000), 'c3': (15000, 24000)}


def build():
    rng = random.Random(11)
    ref = random_seq(rng, REF_LENGTH)
    doubled = ref + ref
    contigs = {name: doubled[start:end] for name, (start, end) in TRUTH.items()}
    contigs['c3'] = revcomp(contigs['c3'])
    return ref, contigs


def run() -> Outcome:
    ws, tmp = workspace()
    ref, contigs = build()
    with tmp:
        ref_fa = ws.write_fasta('ref.fa', [('ref', ref)])
        recon_fa = ws.write_fasta('reconstruction.fa', contigs.items())
        out = outputs(ws.run('--plasmid-reconstruction', recon_fa, '--ref-plasmid', ref_fa, '-o', ws.path('out')))

        # For comparison: minimap2 directly, against the ref as given (linear) and doubled
        c1_fa = ws.write_fasta('c1.fa', [('c1', contigs['c1'])])
        c2_fa = ws.write_fasta('c2.fa', [('c2', contigs['c2'])])
        doubled_fa = ws.write_fasta('ref_doubled.fa', [('ref_doubled', ref + ref)])
        linear = parsers.parse_sam(ws.minimap2('-a', '-cx', 'asm10', '--eqx', '-o', ws.path('c1_linear.sam'), ref_fa, c1_fa))
        doubled = parsers.parse_sam(ws.minimap2('-a', '-cx', 'asm10', '--eqx', '-o', ws.path('c2_doubled.sam'), doubled_fa, c2_fa))
        single = parsers.parse_sam(ws.minimap2('-a', '-cx', 'asm10', '--eqx', '-o', ws.path('c2_single.sam'), ref_fa, c2_fa))

    linear_mapped = [a for a in linear if a['is_mapped']]
    c1_rows = [r for r in out['summary'] if r['contig_id'] == 'c1']
    sample = ''.join(contigs[c] if c != 'c3' else revcomp(contigs[c]) for c in ('c2', 'c3', 'c1'))
    scaffold = out['scaffold']
    facts = {
        'ref_length': REF_LENGTH,
        'c1_length': len(contigs['c1']),
        'linear_segments': len(linear_mapped),
        'linear_segments_text': plural(len(linear_mapped), 'alignment'),
        'linear_pieces': ', '.join(f"{a['ref_start']:,}–{a['ref_end']:,}" for a in sorted(linear_mapped, key=lambda a: -a['ref_start'])),
        'c1_rows': len(c1_rows),
        'c1_spans_origin': c1_rows[0]['spans_origin'] == 'True',
        'c1_start': int(c1_rows[0]['alignment_start']),
        'c1_end': int(c1_rows[0]['alignment_end']),
        'scaffold_is_sample': len(scaffold) == len(sample) and scaffold.upper() in (sample + sample).upper(),
        'scaffold_starts_with_c1_tail': scaffold.upper().startswith(contigs['c1'][-6000:].upper()),
        'doubled_mapq': max(a['mapq'] for a in doubled if a['is_mapped']),
        'single_mapq': max(a['mapq'] for a in single if a['is_mapped']),
    }
    linear_map = {'ref_length': REF_LENGTH, 'arcs': [
        {'contig': 'c1', 'start': a['ref_start'], 'end': a['ref_end']} for a in linear_mapped]}
    return Outcome(ws.commands, facts, {
        'map': plasmid_map(out, REF_LENGTH, TRUTH),
        'linear': linear_map,
    })


SCENARIO = Scenario(
    slug='origin',
    title='Circular references',
    summary="A contig spanning the reference's origin is placed as one alignment, and the reference isn't "
            "doubled, so mapping qualities stay meaningful.",
    text="""\
A plasmid's sequence has to start somewhere, but the start (the *origin*) is arbitrary: nothing in
the molecule marks it. Here a {ref_length:,} bp reference is cut into three contigs, and contig **c1**
({c1_length:,} bp) crosses the origin: it covers the last 6,000 bp of the reference and then the
first 6,000.

minimap2, which plasmid-surgery uses for alignment, treats a reference as linear. Aligning c1 to the
reference as given produces {linear_segments_text} ({linear_pieces}): one piece on each side of the
origin. So plasmid-surgery aligns each contig twice, against the reference and against a copy
rotated by half its length, and keeps whichever alignment covers more of the contig in fewer pieces.
For c1 that's the rotated one: in `summary.tsv` it's a single row from {c1_start:,} to {c1_end:,}
(`alignment_end` < `alignment_start`, and `spans_origin` is **{c1_spans_origin}**). All coordinates
are converted back to the original reference.

[[figure]]

When the scaffold is built, an alignment that spans the origin is split there. The scaffold starts
at the reference's origin, so it begins with the end of c1 (**{scaffold_starts_with_c1_tail}**), and
it is the complete plasmid: **{scaffold_is_sample}**.

**Why not double the reference?** Concatenating the reference with itself would also let c1 align
in one piece. But then every contig aligns equally well to both copies, so minimap2 can't say which
copy is right and gives every alignment mapping quality 0. Here c2 aligns to the doubled reference
with mapping quality {doubled_mapq}, and to the single reference with {single_mapq}.
plasmid-surgery uses mapping quality, and alternative placements, to find contigs that really could
belong in more than one place (see [Repeats](repeats.md)), so it keeps the reference single.

The trade-off is two alignments per contig, and one remaining edge case: a contig spanning both the
origin and the halfway point (longer than half the plasmid) is split either way.
""",
    run=run,
    test='tests/test_scenarios.py::test_origin',
)
