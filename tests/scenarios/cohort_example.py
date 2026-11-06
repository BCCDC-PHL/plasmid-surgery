"""An example cohort for the cohort report: five samples against one reference (the 'Feature
annotations' scenario's plasmid), each with a known difference. Used by the cohort report's tests and
by the documentation's example."""

from __future__ import annotations

import random
from pathlib import Path

from tests.scenarios import features as features_scenario
from tests.scenarios.common import Workspace, random_seq

U, R = features_scenario.U, features_scenario.REPEAT


def segments():
    ref, *_ = features_scenario.build()
    return ref, {'u1': ref[:U], 'repeat': ref[U:U + R], 'u2': ref[U + R:2 * U + R], 'u3': ref[2 * U + 2 * R:]}


def samples():
    """{sample: (reconstruction contigs, chromosome contigs, sample genes as (name, contig, start, end))}"""
    ref, seg = segments()
    rng = random.Random(41)
    novel = random_seq(rng, 2000)
    chrom = lambda n: ('chrom_' + n, random_seq(random.Random(len(n)), 30000))  # noqa: E731
    return {
        'S1_complete': ([('u1', seg['u1']), ('repeat', seg['repeat']), ('u2', seg['u2']), ('u3', seg['u3'])],
                        [chrom('a')], [('geneA', 'u1', 1000, 1900)]),
        'S2_misbinned': ([('u1', seg['u1']), ('repeat', seg['repeat']), ('u2', seg['u2'])],
                         [('u3', seg['u3']), chrom('bb')], [('geneA', 'u1', 1000, 1900), ('geneMisbinned', 'u3', 2000, 2900)]),
        'S3_incomplete': ([('u1', seg['u1']), ('repeat', seg['repeat']), ('u2a', seg['u2'][:4000]), ('u3', seg['u3'])],
                          [chrom('ccc')], [('geneA', 'u1', 1000, 1900)]),
        'S4_indels': ([('u1', seg['u1'][:3000] + seg['u1'][4500:]), ('repeat', seg['repeat']),
                       ('u2', seg['u2'][:4500] + random_seq(rng, 1500) + seg['u2'][4500:]), ('u3', seg['u3'])],
                      [chrom('dddd')], [('geneA', 'u1', 1000, 1900)]),
        'S5_novel_gene': ([('u1', seg['u1'][:6000] + novel + seg['u1'][6000:]), ('repeat', seg['repeat']), ('u2', seg['u2']),
                           ('u3', seg['u3'])], [chrom('eeeee')], [('geneA', 'u1', 1000, 1900), ('geneNovel', 'u1', 6300, 7500)]),
    }


def build_runs(root: Path) -> list[Path]:
    """Run plasmid-surgery (with --report) for each sample in ``root``; returns their output directories."""
    ws = Workspace(Path(root))
    ref, _ = segments()
    _, _, _, ref_features, _ = features_scenario.build()
    ref_fa = ws.write_fasta('pREF1.fa', [('pREF1', ref)])
    ref_gff = ws.write_text('pREF1.gff3', features_scenario.gff3(ref_features, 'pREF1'))
    runs = []
    for sample, (reconstruction, chromosome, genes) in samples().items():
        recon_fa = ws.write_fasta(f'{sample}_plasmid_1.fasta', reconstruction)
        chrom_fa = ws.write_fasta(f'{sample}_chromosome.fasta', chromosome)
        assembly_fa = ws.write_fasta(f'{sample}_assembly.fasta', reconstruction + chromosome)
        genes_tsv = ws.write_text(f'{sample}.hamronized.tsv', features_scenario.hamronization(genes))
        runs.append(ws.run('--plasmid-reconstruction', recon_fa, '--ref-plasmid', ref_fa, '--assembly', assembly_fa,
                           '--bins', recon_fa, chrom_fa, '--ref-features', ref_gff, '--features', genes_tsv,
                           '--sample-id', sample, '--report', '-o', ws.path(f'results/{sample}')))
    return runs
