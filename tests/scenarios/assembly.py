"""Whole-assembly mode: classifying the assembly's other contigs that align to the ref."""

from __future__ import annotations

import random

from tests.scenarios.common import Outcome, Scenario, outputs, plasmid_map, random_seq, workspace

SEGMENT = 7500
NAMES = ['A', 'B', 'C', 'D']   # the ref is A B C D, each SEGMENT long
FOREIGN = 2500                  # sequence not on the ref, at the end of contig 'D_junction'


def build():
    rng = random.Random(59)
    segments = {n: random_seq(rng, SEGMENT) for n in NAMES}
    ref = ''.join(segments[n] for n in NAMES)
    reconstruction = [('A depth=2.00x', segments['A']), ('C depth=2.00x', segments['C'])]
    chromosome = [
        ('B depth=2.10x', segments['B']),                                         # a misbinned plasmid contig
        ('D_junction depth=1.90x', segments['D'] + random_seq(rng, FOREIGN)),      # plasmid contig plus foreign sequence
        ('chrom_shared depth=1.00x', random_seq(rng, 9000) + segments['A'][3000:5000] + random_seq(rng, 9000)),
        ('chrom_other depth=1.00x', random_seq(rng, 20000)),                       # nothing to do with the plasmid
    ]
    return ref, reconstruction, chromosome


def run() -> Outcome:
    ws, tmp = workspace()
    ref, reconstruction, chromosome = build()
    with tmp:
        ref_fa = ws.write_fasta('ref.fa', [('ref', ref)])
        recon_fa = ws.write_fasta('plasmid_1.fasta', reconstruction)
        chrom_fa = ws.write_fasta('chromosome.fasta', chromosome)
        assembly_fa = ws.write_fasta('assembly.fasta', reconstruction + chromosome)
        out = outputs(ws.run('--plasmid-reconstruction', recon_fa, '--ref-plasmid', ref_fa, '--assembly', assembly_fa,
                             '--bins', recon_fa, chrom_fa, '-o', ws.path('out')))

    sources = {r['contig_id']: r for r in out['contig_sources']}
    a = out['assessment']
    truth = {n: (i * SEGMENT, (i + 1) * SEGMENT) for i, n in enumerate(NAMES)}
    columns = ['contig_id', 'bin', 'classification', 'contig_fraction_on_ref', 'depth_relative_to_reconstruction']
    facts = {
        'segment': SEGMENT, 'foreign': FOREIGN, 'ref_length': len(ref),
        'B_class': sources['B']['classification'],
        'D_class': sources['D_junction']['classification'],
        'D_fraction': float(sources['D_junction']['contig_fraction_on_ref']),
        'shared_class': sources['chrom_shared']['classification'],
        'shared_fraction': float(sources['chrom_shared']['contig_fraction_on_ref']),
        'other_listed': 'chrom_other' in sources,
        'B_depth': float(sources['B']['depth_relative_to_reconstruction']),
        'D_depth': float(sources['D_junction']['depth_relative_to_reconstruction']),
        'shared_depth': float(sources['chrom_shared']['depth_relative_to_reconstruction']),
        'covered': float(a['ref_fraction_covered']),
        'assembly_covered': float(a['assembly_ref_fraction_covered']),
        'unique_gap_length': int(a['unique_gap_length']),
        'assembly_unique_gap_length': int(a['assembly_unique_gap_length']),
        'missed_ref_bases': int(a['missed_ref_bases']),
        'num_missed': int(a['num_missed_contigs']),
        'num_shared': int(a['num_shared_contigs']),
    }
    tables = {'sources': (columns, [{c: sources[cid][c] for c in columns} for cid in sources])}
    data = {'map': plasmid_map(out, len(ref), truth), 'tables': tables,
            'sources': [{k: r[k] for k in ('contig_id', 'classification', 'ref_start', 'ref_end', 'bin',
                                           'contig_fraction_on_ref')} for r in out['contig_sources']]}
    return Outcome(ws.commands, facts, data)


SCENARIO = Scenario(
    slug='assembly',
    title='Whole-assembly mode',
    summary="With the full assembly, every contig aligning to the reference is classified: missed (likely in "
            "the wrong bin), shared (only partly on the reference) or repeat; read depth is supporting evidence.",
    text="""\
A reconstruction can lack part of a plasmid because it wasn't assembled, or because the contigs were
assembled but put in another bin, or in the chromosome. With `--assembly` (and optionally the other
bins, `--bins`), plasmid-surgery tells the two apart.

Here the reference is four {segment:,} bp segments, A–D. The reconstruction has A and C. The assembly's
chromosome bin has B, a contig with D plus {foreign:,} bp of other sequence, a chromosome contig that
shares 2,000 bp with A, and one unrelated contig. Contig headers carry read depth, as from Unicycler
(`depth=2.00x`). `contig_sources.tsv` lists every assembly contig that aligns to the reference:

{sources}

[[figure]]

* **B is `missed`**: all of it aligns to the reference, at one location, so it's very likely a plasmid
  contig binned as chromosome. Its depth relative to the reconstruction's contigs, {B_depth}, agrees.
* **The chromosome contig sharing part of A is `shared`**: only {shared_fraction:.0%} of it is on the
  reference, so it's more likely chromosome with some sequence in common (for example a transposon)
  than a binning error. Its relative depth is {shared_depth}.
* The unrelated contig doesn't align to the reference, so it isn't in the table.

The metrics separate the two kinds of loss: the reconstruction has {unique_gap_length:,} bp of
unique gaps, the assembly as a whole {assembly_unique_gap_length:,} bp; `missed_ref_bases` =
{missed_ref_bases:,}, the part of the reference only in `missed` contigs. B is written to
`missed_contigs.fa`.

**The trade-off is the threshold.** A contig is `missed` only if at least `--min-missed-fraction`
(default 0.8) of it aligns to the reference. The contig with D has {D_fraction:.0%} on the reference,
so it's `{D_class}`, though here it really is a plasmid contig: its other {foreign:,} bp are a junction
to sequence the reference doesn't have. A lower threshold would catch it, and also call more
chromosome contigs that just share a transposon `missed`. That's why the classification comes with
the evidence: `contig_fraction_on_ref` and relative depth ({D_depth} here, like a plasmid contig).
Depth from short contigs is noisy, and everything here is only as good as the reference.
""",
    run=run,
    test='tests/test_scenarios.py::test_assembly',
)
