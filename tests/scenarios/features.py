"""Feature annotations: reference genes' completeness, and where the sample's genes are."""

from __future__ import annotations

import random

from tests.scenarios.common import Outcome, Scenario, outputs, plasmid_map, random_seq, workspace

U = 9000
REPEAT = 1200
DELETED = (3000, 4500)         # within u1, missing from the sample
INSERTION_LENGTH = 1500        # novel sequence in the sample's u2, at its midpoint

HAMRONIZATION_COLUMNS = ['input_file_name', 'gene_symbol', 'gene_name', 'reference_database_name',
                         'reference_database_version', 'reference_accession', 'analysis_software_name',
                         'analysis_software_version', 'genetic_variation_type', 'drug_class',
                         'input_gene_start', 'input_gene_stop', 'input_sequence_id', 'strand_orientation']


def build():
    rng = random.Random(71)
    u1, u2, u3 = (random_seq(rng, U) for _ in range(3))
    repeat = random_seq(rng, REPEAT)
    ref = u1 + repeat + u2 + repeat + u3
    insertion = random_seq(rng, INSERTION_LENGTH)
    sample_u1 = u1[:DELETED[0]] + u1[DELETED[1]:]
    sample_u2 = u2[:U // 2] + insertion + u2[U // 2:]
    reconstruction = [('u1', sample_u1), ('repeat', repeat), ('u2', sample_u2)]
    chromosome = [('u3', u3)]
    copy1, copy2 = (U, U + REPEAT), (2 * U + REPEAT, 2 * U + 2 * REPEAT)
    u3_start = copy2[1]
    ref_features = [   # name, start, end (0-based), on the ref
        ('geneA', 1000, 1900),
        ('geneDeleted', DELETED[0] + 200, DELETED[1] - 200),
        ('tnpA_copy1', copy1[0] + 100, copy1[1] - 100),
        ('tnpA_copy2', copy2[0] + 100, copy2[1] - 100),
        ('geneMisbinned', u3_start + 2000, u3_start + 2900),
    ]
    sample_features = [   # name, contig, start, end (0-based)
        ('geneA', 'u1', 1000, 1900),
        ('geneInInsertion', 'u2', U // 2 + 300, U // 2 + 1200),
        ('geneMisbinned', 'u3', 2000, 2900),
    ]
    return ref, reconstruction, chromosome, ref_features, sample_features


def gff3(features, seq_id):
    lines = ['##gff-version 3']
    for i, (name, start, end) in enumerate(features, 1):
        lines.append('\t'.join([seq_id, 'annotator', 'CDS', str(start + 1), str(end), '.', '+', '0', f'ID=f{i};gene={name}']))
    return '\n'.join(lines) + '\n'


def hamronization(features):
    lines = ['\t'.join(HAMRONIZATION_COLUMNS)]
    for name, contig, start, end in features:
        row = {'input_file_name': 'assembly.fasta', 'gene_symbol': name, 'gene_name': name, 'reference_database_name': 'db',
               'reference_database_version': '1', 'reference_accession': '-', 'analysis_software_name': 'annotator',
               'analysis_software_version': '1', 'genetic_variation_type': 'gene_presence_detected', 'drug_class': '',
               'input_gene_start': str(start + 1), 'input_gene_stop': str(end), 'input_sequence_id': contig,
               'strand_orientation': '+'}
        lines.append('\t'.join(row[c] for c in HAMRONIZATION_COLUMNS))
    return '\n'.join(lines) + '\n'


def run() -> Outcome:
    ws, tmp = workspace()
    ref, reconstruction, chromosome, ref_features, sample_features = build()
    with tmp:
        ref_fa = ws.write_fasta('ref.fa', [('ref', ref)])
        recon_fa = ws.write_fasta('plasmid_1.fasta', reconstruction)
        chrom_fa = ws.write_fasta('chromosome.fasta', chromosome)
        assembly_fa = ws.write_fasta('assembly.fasta', reconstruction + chromosome)
        ref_gff = ws.write_text('ref.gff3', gff3(ref_features, 'ref'))
        sample_tsv = ws.write_text('assembly.hamronized.tsv', hamronization(sample_features))
        out = outputs(ws.run('--plasmid-reconstruction', recon_fa, '--ref-plasmid', ref_fa, '--assembly', assembly_fa,
                             '--bins', recon_fa, chrom_fa, '--ref-features', ref_gff, '--features', sample_tsv,
                             '-o', ws.path('out')))

    ref_rows = {r['name']: r for r in out['ref_features']}
    sample_rows = {r['name']: r for r in out['features']}
    ref_columns = ['name', 'reconstruction_fraction', 'inferred_fraction', 'reconstruction_status', 'assembly_status', 'found_in']
    sample_columns = ['name', 'contig_id', 'contig_classification', 'bin', 'fraction_on_ref', 'ref_status', 'ref_start', 'ref_end']
    tnp_statuses = sorted(ref_rows[n]['reconstruction_status'] for n in ('tnpA_copy1', 'tnpA_copy2'))
    facts = {
        'deleted_length': DELETED[1] - DELETED[0],
        'insertion_length': INSERTION_LENGTH,
        'repeat': REPEAT,
        'geneA_status': ref_rows['geneA']['reconstruction_status'],
        'deleted_status': ref_rows['geneDeleted']['reconstruction_status'],
        'deleted_assembly_status': ref_rows['geneDeleted']['assembly_status'],
        'tnp_statuses': ' and '.join(f'`{s}`' for s in reversed(tnp_statuses)),
        'tnp_assembly_statuses': sorted(ref_rows[n]['assembly_status'] for n in ('tnpA_copy1', 'tnpA_copy2')),
        'misbinned_status': ref_rows['geneMisbinned']['reconstruction_status'],
        'misbinned_assembly_status': ref_rows['geneMisbinned']['assembly_status'],
        'misbinned_found_in': ref_rows['geneMisbinned']['found_in'],
        'sample_geneA_ref': (sample_rows['geneA']['ref_start'], sample_rows['geneA']['ref_end']),
        'insertion_gene_status': sample_rows['geneInInsertion']['ref_status'],
        'insertion_gene_class': sample_rows['geneInInsertion']['contig_classification'],
        'sample_misbinned_class': sample_rows['geneMisbinned']['contig_classification'],
        'sample_misbinned_bin': sample_rows['geneMisbinned']['bin'],
        'metrics': {k: out['assessment'][k] for k in out['assessment'] if 'feature' in k},
    }
    tables = {
        'ref_table': (ref_columns, [{c: r[c] for c in ref_columns} for r in out['ref_features']]),
        'sample_table': (sample_columns, [{c: r[c] for c in sample_columns} for r in out['features']]),
    }
    return Outcome(ws.commands, facts, {'map': plasmid_map(out, len(ref)), 'tables': tables})


SCENARIO = Scenario(
    slug='features',
    title='Feature annotations',
    summary="Reference genes are reported present, inferred, partial or absent, and where else in the assembly "
            "they are; sample genes by their contig's classification and whether they're on the reference.",
    text="""\
Gene annotations (any kind: resistance, virulence, replicons, mobility, ...) can be given for the
reference (`--ref-features`) and for the sample's assembly (`--features`), as hAMRonization TSV, GFF3
or BED. No class of gene is treated specially.

Here the reference has three unique segments separated by two copies of a {repeat:,} bp repeat. The
sample is missing {deleted_length:,} bp of the first segment, has {insertion_length:,} bp of novel
sequence in the second, and its third segment is in the chromosome bin. The reference's annotations are
GFF3; the assembly's are hAMRonization TSV.

[[figure]]

**Reference features** (`ref_features.tsv`): how completely the reconstruction, and the assembly,
contain each one.

{ref_table}

* `geneA` is `{geneA_status}`.
* `geneDeleted` lies in the sequence the sample doesn't have: `{deleted_status}` in the reconstruction,
  and `{deleted_assembly_status}` in the assembly too.
* The repeat's gene, at both copies, is {tnp_statuses}: the copy where the collapsed repeat contig's
  primary alignment is, and the other copy, covered only by an alternative placement (see
  [Repeats](repeats.md)). The assembly status is the same, for the same reason.
* `geneMisbinned` is `{misbinned_status}` from the reconstruction but `{misbinned_assembly_status}` in the
  assembly, `found_in` the `{misbinned_found_in}` bin: a binning error, at gene level.

**Sample features** (`features.tsv`): for genes on the reconstruction, or on any assembly contig aligning
to the reference.

{sample_table}

* A gene's `fraction_on_ref` counts only bases matched to reference bases, so the gene inside the novel
  insertion is `{insertion_gene_status}`, though its contig is `{insertion_gene_class}` and aligns: it's
  sequence the reference doesn't have.
* The gene on the chromosome-binned contig is on a `{sample_misbinned_class}` contig in the
  `{sample_misbinned_bin}` bin, and on the reference.

The two views complement each other. Sample annotations can't report a gene that isn't assembled well
enough to be called, for example one split across short repeat contigs; the reference's annotation
still shows it as `partial` or `absent`. Reference annotations can't report genes the reference
doesn't have; the sample's annotations show them as `not_on_ref`.
""",
    run=run,
    test='tests/test_scenarios.py::test_features',
)
