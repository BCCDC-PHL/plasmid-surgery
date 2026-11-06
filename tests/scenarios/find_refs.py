"""Finding references: screening a plasmid database with mash, confirming with minimap2, and ranking."""

from __future__ import annotations

import random

from tests.scenarios.common import (Outcome, Scenario, mash_sketch, mutate, random_seq, read_tsv, run_refs,
                                    workspace)


def build():
    rng = random.Random(97)
    pa = random_seq(rng, 30000)
    pb = random_seq(rng, 40000)
    db = {
        'pA': mutate(rng, pa, 30),                                              # the same plasmid (0.1% SNPs)
        'pA_variant': mutate(rng, pa[:15000] + random_seq(rng, 6000) + pa[15000:], 150),   # + 6 kb, 0.5% SNPs
        'pA_core': pa[5000:25000] + random_seq(rng, 2000),                       # shares 20 kb of pA
        'pB': pb,
        'pB_small': pb[10000:25000] + random_seq(rng, 1000),                     # shares the fragment of pB
    }
    for i in range(1, 6):
        db[f'unrelated_{i}'] = random_seq(rng, rng.randint(20000, 60000))
    complete = [('a1', pa[:12000]), ('a2', pa[12000:21000]), ('a3', pa[21000:])]
    fragment = [('b1', pb[10000:17500]), ('b2', pb[17500:25000])]
    assembly_b = fragment + [('b_rest1', pb[25000:]), ('b_rest2', pb[:10000]), ('chrom', random_seq(rng, 50000))]
    return db, complete, fragment, assembly_b


def table(rows, coverage='ref_fraction_covered'):
    return {r['ref_id']: {'rank': int(r['rank']), 'aligned': float(r['reconstruction_fraction_aligned']),
                          'covered': float(r[coverage]) if r[coverage] != 'NA' else None,
                          'ref_covered': float(r['ref_fraction_covered']),
                          'score': float(r['score']), 'mash': float(r['mash_distance']),
                          'est_in': float(r['est_reconstruction_containment']),
                          'est_ref': float(r['est_ref_containment']), 'length': int(r['ref_length'])}
            for r in rows}


def run() -> Outcome:
    ws, tmp = workspace()
    db, complete, fragment, assembly_b = build()
    with tmp:
        db_fa = ws.write_fasta('plasmids.fa', db.items())
        db_msh = mash_sketch(ws, db_fa, 'plasmids')
        a_fa = ws.write_fasta('reconstruction_a.fa', complete)
        b_fa = ws.write_fasta('reconstruction_b.fa', fragment)
        asm_fa = ws.write_fasta('assembly_b.fa', assembly_b)
        common = ['--db-sketch', db_msh, '--db-fasta', db_fa]
        a = read_tsv(run_refs(ws, '--plasmid-reconstruction', a_fa, *common, '-o', ws.path('refs_a')) / 'candidates.tsv')
        b = read_tsv(run_refs(ws, '--plasmid-reconstruction', b_fa, *common, '-o', ws.path('refs_b')) / 'candidates.tsv')
        b_asm = read_tsv(run_refs(ws, '--plasmid-reconstruction', b_fa, *common, '--assembly', asm_fa,
                                  '-o', ws.path('refs_b_assembly')) / 'candidates.tsv')

    ta, tb, tb_asm = table(a), table(b), table(b_asm, 'assembly_ref_fraction_covered')
    best = lambda t: min(t, key=lambda k: t[k]['rank'])  # noqa: E731
    nearest = lambda t: min(t, key=lambda k: (t[k]['mash'], t[k]['rank']))  # noqa: E731
    facts = {
        'a_best': best(ta), 'a_nearest': nearest(ta),
        'a_variant_covered': ta['pA_variant']['covered'], 'a_variant_score': ta['pA_variant']['score'],
        'a_core_aligned': ta['pA_core']['aligned'], 'a_core_score': ta['pA_core']['score'],
        'a_core_covered': ta['pA_core']['covered'],
        'a_best_score': ta[best(ta)]['score'],
        'a_core_est_ref': ta['pA_core']['est_ref'], 'a_core_mash': ta['pA_core']['mash'],
        'a_unrelated_candidates': sum(1 for k in ta if k.startswith('unrelated')),
        'b_best': best(tb), 'b_small_score': tb['pB_small']['score'], 'b_large_score': tb['pB']['score'],
        'b_large_aligned': tb['pB']['aligned'], 'b_large_covered': tb['pB']['covered'],
        'b_large_mash': tb['pB']['mash'], 'b_small_mash': tb['pB_small']['mash'],
        'b_large_est_in': tb['pB']['est_in'],
        'b_asm_best': best(tb_asm), 'b_asm_large_covered': tb_asm['pB']['covered'],
        'b_asm_large_score': tb_asm['pB']['score'], 'b_asm_small_score': tb_asm['pB_small']['score'],
        'b_asm_small_covered': tb_asm['pB_small']['covered'],
    }
    columns = ['ref_id', 'ref_length', 'mash_distance', 'est_reconstruction_containment', 'est_ref_containment',
               'reconstruction_fraction_aligned', 'ref_fraction_covered', 'score']
    pretty = lambda rows: [{c: r[c] for c in columns} for r in rows if not r['ref_id'].startswith('unrelated')]  # noqa: E731
    tables = {'a_table': (columns, pretty(a)), 'b_table': (columns + ['assembly_ref_fraction_covered'],
              [{**{c: r[c] for c in columns}, 'assembly_ref_fraction_covered': r['assembly_ref_fraction_covered']}
               for r in b_asm if not r['ref_id'].startswith('unrelated')])}
    return Outcome(ws.commands, facts, {'tables': tables, 'a': ta, 'b': tb, 'b_asm': tb_asm})


SCENARIO = Scenario(
    slug='find-refs',
    title='Finding references',
    summary="plasmid-surgery-refs screens a plasmid database with mash, confirms candidates with minimap2, and "
            "ranks them by how much of the reconstruction they explain times how much of them is covered.",
    text="""\
`plasmid-surgery-refs` looks for suitable references in a plasmid database: a mash sketch of it, and
the FASTA it was made from. Here the database has two plasmids and some of their relatives, and five
unrelated plasmids:

* `pA` (30 kb), `pA_variant` (pA with 6 kb more, and 0.5% SNPs), and `pA_core` (20 kb of pA, plus 2 kb of
  its own);
* `pB` (40 kb) and `pB_small` (15 kb of pB, plus 1 kb of its own).

**Screening.** The reconstruction is sketched as a whole and compared with every database sketch in
one `mash dist`. mash's distance comes from an estimate of the Jaccard index of the two k-mer sets:
it's symmetric, and high whenever the two differ in size, even if one contains the other. With the two
sizes (from the sketch), the same estimate gives *containment* in each direction: how much of the
reconstruction is in the reference, and how much of the reference is in the reconstruction. The
candidates are the top 20 by containment of the reconstruction, plus the top 20 by mash distance, so
that neither a size difference nor a shared fragment hides a good candidate. They're then extracted
from the FASTA (with its `.fai` index if there is one, otherwise in one pass) and aligned with minimap2
for exact figures.

**Ranking.** Each candidate's score is

    score = (fraction of the reconstruction aligned to it) × (fraction of it covered)

Both are near 1 only when the reference is the same plasmid as the reconstruction: one that
explains the reconstruction but has much more besides is penalised by the second factor, and one
that's covered but explains only part of the reconstruction by the first.

A complete reconstruction of pA:

{a_table}

`{a_best}` ranks first (score {a_best_score}), and is also the nearest by mash distance. The variant
explains the whole reconstruction, but its extra 6 kb is uncovered ({a_variant_covered:.0%} covered;
score {a_variant_score}). `pA_core` is {a_core_covered:.0%} covered but explains only
{a_core_aligned:.0%} of the reconstruction (score {a_core_score}). The unrelated plasmids share no
k-mers with the reconstruction, so none of them is a candidate.

**The trade-off: a fragment.** Reconstruction B is two contigs from the middle of pB, 15 kb of its
40 kb. Both `pB` and `pB_small` explain all of it, but the reconstruction covers almost all of
`pB_small` and only {b_large_covered:.0%} of `pB`, so against the reconstruction alone `{b_best}` ranks
first (score {b_small_score}, against {b_large_score} for pB). This is the limitation described in
[The choice of reference](reference-choice.md): from the reconstruction alone, a fragment of a larger
plasmid looks like a smaller one. With `--assembly`, the coverage factor is the whole assembly's
coverage of each candidate instead. The rest of pB is in the assembly, so pB is
{b_asm_large_covered:.0%} covered and ranks first (score {b_asm_large_score}); `pB_small`'s own 1 kb isn't in the assembly
(score {b_asm_small_score}):

{b_table}

[[figure]]

The candidates file reports each figure separately, so a different balance can be applied. If the best
candidate explains less than 80% of the reconstruction, a warning says the database may not contain a
closely related plasmid.
""",
    run=run,
    test='tests/test_scenarios.py::test_find_refs',
)
