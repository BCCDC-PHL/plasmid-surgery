# Reads

Coverage against a reference says what a reconstruction lacks, but not why. A reference region no
contig covers may be:

- **not in the sample**: the sample's plasmid differs from the reference there (likely when the
  reference is only a relative), or
- **in the sample, but lost**: by the assembler, or, with [`--assembly`](whole-assembly.md), by the
  binning.

Without reads the two look the same. Given the sample's reads (`--short-reads` and/or
`--long-reads`), plasmid-surgery aligns them to each reference and reports read depth along it, so
each gap can be called one or the other.

```bash
plasmid-surgery --plasmid-reconstruction plasmid_AA001.fasta --ref-plasmid ref.fa \
  --assembly contigs.fa --bins plasmid_*.fasta chromosome.fasta \
  --short-reads sample_R1.fastq.gz sample_R2.fastq.gz -t 8 -o out
```

## Read depth

Reads are aligned with minimap2 (`sr` for short reads; `--long-read-preset`, default `map-ont`, for
long reads; `lr:hq` suits recent nanopore data), one call per read type for all references together.
Every alignment counts: primary, secondary (another equally good placement, e.g. of a repeat) and
supplementary (part of a read aligned separately, e.g. a long read crossing the reference's origin),
if it is at least `--min-read-identity` identical (default 0.95 for short reads, 0.85 for long reads).
Deletions shorter than 50 bp within a read's alignment count as covered, so sequencing errors don't
punch holes in the depth. Depth from short and long reads is added together.

A reference base is **in the reads** if its depth is at least `--min-read-depth` (default 3).

Reads are aligned to the references only, not to the rest of the genome. Reads from the chromosome
or another plasmid that share sequence with the reference (an IS element, a transposon) count as
well: that sequence is in the sample, though not necessarily on this plasmid. Such regions are
usually also [repeat gaps](coverage.md), so look at unique gaps first.

Near a structural difference between the sample and the reference (e.g. a deletion), reads crossing
it may align on one side only, so depth dips for about a read length on each side. At typical depths
this stays above `--min-read-depth`.

## Outputs

- `reads_depth.bedgraph`: depth along the reference.
- `gap_support.tsv`: one row per gap with the fraction of it in the reads and a call:
  `in_reads` (≥ 90% of the gap), `not_in_reads` (≤ 10%) or `partly_in_reads`. With `--assembly`, the
  assembly's own gaps (what no assembly contig covers) are listed too, as `assembly_gap_<n>_<type>`.
- In `assessment.tsv`: `reads_ref_fraction_covered`, `reads_median_depth`,
  `unique_gap_in_reads_length` and `assembly_unique_gap_in_reads_length`.
- In `ref_features.tsv`: `reads_fraction`, the fraction of each feature in the reads (is the gene in
  the sample at all?).

## Reading the numbers

With `--assembly`, a reference's unique sequence splits three ways:

| | in the reconstruction | not in the reconstruction |
|---|---|---|
| **in the assembly** | present | lost by binning (`missed_ref_bases`) |
| **not in the assembly, in the reads** | | lost by the assembler (`assembly_unique_gap_in_reads_length`) |
| **not in the reads** | | not in the sample: the sample differs from the reference |

`reads_ref_fraction_covered` is how much of the reference the sample contains at all, the ceiling
for any reconstruction against it. It is also a check on the inputs: if the reads cover little of a
reference that the reconstruction matches well, the reads and assembly probably aren't from the same
sample.

## Many bins, one sample

The reads are the same for every bin of a sample. With `--reads-cache DIR`, depth along each
reference is stored in `DIR`, keyed by the reads files (path and size), the minimap2 preset, the
identity threshold and the reference sequence, and later runs against the same references reuse it
instead of aligning the reads again. Runs may share a cache directory concurrently (files are
written atomically), but runs started at the same time each align the reads once before the cache
is filled; run one bin first, then the rest.
