# Getting started

## Install

In a dedicated environment (not conda's `base`). plasmid-surgery needs
[minimap2](https://github.com/lh3/minimap2) on `PATH`, and pysam:

```bash
conda env create -f environment.yml   # creates the 'plasmid-surgery' environment
conda activate plasmid-surgery
pip install -e .
plasmid-surgery --version
```

## A first run

You need a plasmid reconstruction (a FASTA of its contigs) and a reference plasmid (a FASTA; only the
first sequence is used):

```bash
plasmid-surgery \
  --plasmid-reconstruction sample1_plasmid.fasta \
  --ref-plasmid CP012345.fa \
  -o results/sample1 \
  -p sample1
```

This writes `results/sample1/sample1_assessment.tsv`, `sample1_summary.tsv`, `sample1_gaps.bed`,
`sample1_events.tsv`, `sample1_unaligned.fa` and `sample1_scaffold.fa`. Without `-p`, the files are
just `assessment.tsv`, and so on. Intermediate alignment files go to a temporary directory that's
removed afterwards, unless you give `--tmpdir`.

Start with `assessment.tsv`: one row of metrics, so rows from many samples can be concatenated into
one table. The most useful at first:

| Field | Meaning |
|---|---|
| `ref_fraction_covered` | How much of the reference the reconstruction's contigs cover |
| `ref_fraction_covered_including_repeats` | The same, counting repeat copies a collapsed repeat contig could fill |
| `unique_gap_length` | Reference sequence nothing in the reconstruction could cover: likely missing |
| `contig_fraction_aligned` | How much of the reconstruction is accounted for by the reference |
| `num_unaligned_contigs` | Contigs that don't align to the reference at all |
| `num_large_insertions`, `num_large_deletions` | Structural differences within aligned contigs |

[Outputs](reference/outputs.md) describes every file and field.

## With the full assembly

If the reconstruction came from an assembly, give the assembly too, and optionally the
reconstruction tool's other bins (its other plasmid FASTAs and its chromosome FASTA):

```bash
plasmid-surgery \
  --plasmid-reconstruction bins/plasmid_1.fasta \
  --ref-plasmid CP012345.fa \
  --assembly sample1.fasta \
  --bins bins/*.fasta \
  -o results/sample1 -p sample1
```

This adds `contig_sources.tsv` (every assembly contig aligning to the reference, classified, with
its bin) and `missed_contigs.fa`, and metrics separating sequence lost by binning
(`missed_ref_bases`) from sequence missing from the assembly itself
(`assembly_unique_gap_length`). See [Whole-assembly mode](concepts/whole-assembly.md).

## With gene annotations

Annotations of the reference, of the assembly, or both, as hAMRonization TSV, GFF3 or BED:

```bash
plasmid-surgery \
  --plasmid-reconstruction bins/plasmid_1.fasta \
  --ref-plasmid CP012345.fa \
  --assembly sample1.fasta --bins bins/*.fasta \
  --ref-features CP012345.amr.hamronized.tsv CP012345.replicons.hamronized.tsv \
  --features sample1.amr.hamronized.tsv sample1.replicons.hamronized.tsv \
  -o results/sample1 -p sample1
```

This adds `ref_features.tsv` and `features.tsv`. See [Feature annotations](concepts/features.md).

## Next

The [end-to-end guide](guides/end-to-end.md) puts these together for a typical workflow, and the
[Behaviour report](behaviour/index.md) shows what each output means on small examples.
