"""An example run for the HTML report: one reconstruction against two references, with the assembly,
bins and annotations, built from the 'Feature annotations' scenario's plasmid. Used by the report's
tests and by the documentation's example report."""

from __future__ import annotations

import random
from pathlib import Path

from tests.scenarios import features as features_scenario
from tests.scenarios.common import Workspace, random_seq


def build_run(root: Path) -> Path:
    """Run plasmid-surgery in ``root``; returns the output directory."""
    ws = Workspace(Path(root))
    ref, reconstruction, chromosome, ref_features, sample_features = features_scenario.build()
    rng = random.Random(5)
    u = features_scenario.U
    # A second, related reference: shares the first two segments, but its third is different
    related = ref[:2 * u + 2 * features_scenario.REPEAT] + random_seq(rng, 6000)
    ref_fa = ws.write_fasta('refs.fa', [('pREF1', ref), ('pREF2_related', related)])
    recon_fa = ws.write_fasta('plasmid_1.fasta', reconstruction)
    chrom_fa = ws.write_fasta('chromosome.fasta', chromosome)
    assembly_fa = ws.write_fasta('assembly.fasta', reconstruction + chromosome)
    ref_gff = ws.write_text('refs.gff3', features_scenario.gff3(ref_features, 'pREF1'))
    sample_tsv = ws.write_text('assembly.hamronized.tsv', features_scenario.hamronization(sample_features))
    return ws.run('--plasmid-reconstruction', recon_fa, '--ref-plasmid', ref_fa, '--assembly', assembly_fa,
                  '--bins', recon_fa, chrom_fa, '--ref-features', ref_gff, '--features', sample_tsv,
                  '--sample-id', 'sample1', '-o', ws.path('out'))
