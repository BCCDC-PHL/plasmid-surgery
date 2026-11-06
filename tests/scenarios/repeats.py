"""Repeats: a collapsed repeat contig, alternative placements, typed gaps, and inferred (lowercase) fill."""

from __future__ import annotations

import random

from tests.scenarios.common import Outcome, Scenario, outputs, plasmid_map, plural, random_seq, workspace

REPEAT_LENGTH = 1200
PARTIAL = 400   # how far contig u3 extends into the repeat copy after it


def build():
    rng = random.Random(23)
    repeat = random_seq(rng, REPEAT_LENGTH)
    u = [random_seq(rng, n) for n in (8000, 8000, 8000, 4000)]
    # ref: u1 IS u2 IS u3 IS u4
    ref = u[0] + repeat + u[1] + repeat + u[2] + repeat + u[3]
    starts = [len(u[0]), len(u[0]) + REPEAT_LENGTH + len(u[1]),
              len(u[0]) + 2 * (REPEAT_LENGTH) + len(u[1]) + len(u[2])]
    copies = [(s, s + REPEAT_LENGTH) for s in starts]
    contigs = {
        'u2': u[1],
        'u3': u[2] + repeat[:PARTIAL],          # runs a little way into the next repeat copy
        'u4_u1': u[3] + u[0],                    # wraps the origin
        'repeat': repeat,                        # all three copies collapsed into one contig
    }
    truth = {
        'u2': (copies[0][1], copies[1][0]),
        'u3': (copies[1][1], copies[2][0] + PARTIAL),
        'u4_u1': (copies[2][1], len(ref) + copies[0][0]),
    }
    return ref, contigs, copies, truth


def run() -> Outcome:
    ws, tmp = workspace()
    ref, contigs, copies, truth = build()
    with tmp:
        ref_fa = ws.write_fasta('ref.fa', [('ref', ref)])
        recon_fa = ws.write_fasta('reconstruction.fa', contigs.items())
        out = outputs(ws.run('--plasmid-reconstruction', recon_fa, '--ref-plasmid', ref_fa, '-o', ws.path('out')))

    rows = out['summary']
    repeat_rows = [r for r in rows if r['contig_id'] == 'repeat']
    primary = next(r for r in repeat_rows if r['is_secondary'] == 'False')
    primary_copy = next(i for i, (s, e) in enumerate(copies, 1) if s == int(primary['alignment_start']))
    a = out['assessment']
    scaffold = out['scaffold']
    lowercase = sum(c.islower() for c in scaffold)
    facts = {
        'repeat_length': REPEAT_LENGTH,
        'partial': PARTIAL,
        'repeat_rows': len(repeat_rows),
        'repeat_placements': int(primary['num_placements']),
        'repeat_mapq': int(primary['mapq']),
        'primary_copy': primary_copy,
        'u3_placements': int(next(r for r in rows if r['contig_id'] == 'u3')['num_placements']),
        'num_repeat_gaps': int(a['num_repeat_gaps']),
        'repeat_gap_length': int(a['repeat_gap_length']),
        'num_unique_gaps': int(a['num_unique_gaps']),
        'covered': float(a['ref_fraction_covered']),
        'covered_with_repeats': float(a['ref_fraction_covered_including_repeats']),
        'lowercase': lowercase,
        'scaffold_ns': scaffold.count('N'),
        'scaffold_is_ref': scaffold.upper() == ref.upper(),
        'num_ambiguous_contigs': int(a['num_ambiguous_contigs']),
        'repeat_gaps_text': plural(int(a['num_repeat_gaps']), 'repeat gap'),
    }
    return Outcome(ws.commands, facts, {'map': plasmid_map(out, len(ref), truth), 'copies': copies})


SCENARIO = Scenario(
    slug='repeats',
    title='Repeats',
    summary="A collapsed repeat contig is placed at every copy it fits: one copy in uppercase, the others "
            "inferred in lowercase, and the gaps it could fill are typed 'repeat', not 'unique'.",
    text="""\
Short-read assemblies usually collapse a repeat, such as an IS element, into one contig: the reads
can't say how many copies there are, or where. Here the reference has three identical copies of a
{repeat_length:,} bp repeat. The reconstruction has the unique sequence between them, and one
`repeat` contig. Contig **u3** also runs {partial} bp into the copy after it, as contigs often
partly extend into a repeat.

[[figure]]

**Alternative placements.** The repeat contig aligns equally well to all three copies. minimap2
picks one (copy {primary_copy}) as its primary alignment, with mapping quality {repeat_mapq}, and
reports the others as secondary alignments. plasmid-surgery asks minimap2 for secondary alignments
scoring at least 90% of the primary, and treats one as an *alternative placement* if it scores at
least 95%. So the repeat contig has {repeat_rows} rows in `summary.tsv`, and `num_placements` =
{repeat_placements}. u3's {partial} bp of repeat would also align to the other copies, but such an
alignment scores far below u3's own (which covers all of u3), so it isn't an alternative placement:
u3 isn't ambiguous, and `num_placements` = {u3_placements}.

**Typed gaps.** The repeat copies not covered by a primary alignment are gaps in the reconstruction,
but not the same kind of gap as missing unique sequence: the reconstruction does have that
sequence, it just can't say it's there. So each gap is typed `repeat` (an alternative placement
covers it) or `unique` (nothing could go there). Here there are {repeat_gaps_text}
({repeat_gap_length:,} bp) and {num_unique_gaps} unique gaps: the reconstruction covers {covered:.1%}
of the reference, or {covered_with_repeats:.1%} counting repeat gaps.

**Inferred fill.** In the scaffold, repeat gaps are filled from the alternative placements, in
lowercase to mark them as inferred: {lowercase:,} lowercase bases, {scaffold_ns} Ns. Uppercased, the
scaffold is the reference exactly: **{scaffold_is_ref}**.

The trade-off: which copy gets the primary alignment, and so stays uppercase, is minimap2's arbitrary
choice. The uppercase copy is no better supported than the lowercase ones, so for an ambiguous
contig (`num_placements` > 1), treat every copy's placement as inferred.
""",
    run=run,
    test='tests/test_scenarios.py::test_repeats',
)
