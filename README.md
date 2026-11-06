# Plasmid Surgery

[![Tests](https://github.com/BCCDC-PHL/plasmid-surgery/actions/workflows/tests.yml/badge.svg)](https://github.com/BCCDC-PHL/plasmid-surgery/actions/workflows/tests.yml)

![Logo](docs/images/plasmid-surgery.png)

Tools like [mob-suite](https://github.com/phac-nml/mob-suite) produce plasmid reconstructions that consist of a set of contigs
without information about the order or orientation of the contigs in the complete plasmid.

This tool takes such a plasmid reconstruction along with a suitable reference plasmid, and attempts to identify the order,
orientation and (approximate) location of each contig in the plasmid reconstruction, based on where it best aligns on the
reference plasmid.

## Usage

```
usage: plasmid-surgery [-h] --plasmid-reconstruction PLASMID_RECONSTRUCTION
                       --ref-plasmid REF_PLASMID [-o OUTDIR] [-p PREFIX]
                       [--tmpdir TMPDIR] [--min-indel-length MIN_INDEL_LENGTH]
                       [--min-unaligned-length MIN_UNALIGNED_LENGTH]
                       [--sample-id SAMPLE_ID]
                       [--output-bases-per-line OUTPUT_BASES_PER_LINE]
                       [--log-level LOG_LEVEL]

options:
  -h, --help            show this help message and exit
  --plasmid-reconstruction PLASMID_RECONSTRUCTION
                        Plasmid reconstruction to be repaired.
  --ref-plasmid REF_PLASMID
                        Reference plasmid to be used for plasmid surgery.
  -o OUTDIR, --outdir OUTDIR
                        Output directory (default: current directory)
  -p PREFIX, --prefix PREFIX
                        Prefix for output filenames, e.g. '-p sample1' gives
                        'sample1_assessment.tsv' (default: no prefix)
  --tmpdir TMPDIR       Directory for intermediate minimap2 files, kept after
                        the run. If not supplied, a temporary directory is
                        used and removed afterwards.
  --min-indel-length MIN_INDEL_LENGTH
                        Minimum length of an insertion/deletion to report as
                        an event (default: 50)
  --min-unaligned-length MIN_UNALIGNED_LENGTH
                        Minimum length of an unaligned contig region to report
                        (default: 100)
  --sample-id SAMPLE_ID
                        Sample ID, used as the scaffold fasta header and the
                        sample_id in assessment output. If not supplied, uses
                        --prefix if given, otherwise the plasmid
                        reconstruction filename (minus file extension)
  --output-bases-per-line OUTPUT_BASES_PER_LINE
                        Number of bases per line in output fasta files
                        (default: 70)
  --log-level LOG_LEVEL
                        Log level (default: 'info')
```

Example:

```
plasmid-surgery \
  --plasmid-reconstruction sample1_plasmid_AA123.fasta \
  --ref-plasmid CP012345.fa \
  -o results/sample1 \
  -p sample1
```

This writes `results/sample1/sample1_assessment.tsv`, `results/sample1/sample1_gaps.bed`, etc.
Without `-p`, files are named `assessment.tsv`, `gaps.bed`, etc.

## Outputs

All outputs are written to `--outdir` (default: current directory), with filenames
prefixed by `<prefix>_` if `--prefix` is given.

| file             | description |
|------------------|-------------|
| `assessment.tsv` | One-row summary of how well the reconstruction matches the ref (see below). Rows from many samples can be concatenated. |
| `gaps.bed`       | BED4 of ref regions not covered by any contig. A gap spanning the origin is written as two records with the same name. |
| `events.tsv`     | Insertions and deletions of at least `--min-indel-length` bp within contig alignments. |
| `unaligned.fa`   | Contig regions of at least `--min-unaligned-length` bp that don't align to the ref (whole unaligned contigs, or unaligned contig ends). Headers are `<contig_id>:<start>-<end>`, 1-based inclusive. |
| `summary.tsv`    | One row per contig alignment (see below). |
| `scaffold.fa`    | Reference-ordered sequence built from the aligned contigs (see [Output sequence](#scaffold-sequence)). Header is the sample ID. |

Intermediate minimap2 files go to a temporary directory that is removed after the run,
unless `--tmpdir` is given, in which case they are kept there.

### Assessment metrics

| field                        | description |
|------------------------------|-------------|
| `sample_id`                  | `--sample-id` if given, else `--prefix` if given, else the reconstruction filename minus extension |
| `ref_id`, `ref_length`       | Ref plasmid ID and length |
| `num_contigs`, `total_contig_length` | Number and total length of contigs in the reconstruction |
| `ref_bases_covered`          | Ref bases covered by at least one contig alignment |
| `ref_fraction_covered`       | `ref_bases_covered / ref_length` |
| `ref_bases_covered_multiple` | Ref bases covered by more than one alignment (contig overlaps, or duplicated sequence) |
| `num_gaps`, `total_gap_length`, `largest_gap_length` | Uncovered ref regions |
| `contig_bases_aligned`       | Contig bases within an alignment to the ref (including inserted bases inside alignments; see `num_large_insertions`) |
| `contig_fraction_aligned`    | `contig_bases_aligned / total_contig_length`. Sequence not accounted for by the ref is in the unaligned output. |
| `num_unaligned_contigs`      | Contigs with no alignment to the ref that passed QC: possible chromosomal contamination or a different plasmid |
| `num_unaligned_regions`, `unaligned_region_bases` | Regions written to the unaligned output |
| `num_large_insertions`, `num_large_deletions` | Events written to the events output |
| `num_ambiguous_alignments`   | Alignments with mapping quality 0: the contig aligns equally well elsewhere on the ref (e.g. a repeated IS element), so its placement is arbitrary |

Assessment metrics use only alignments that passed QC.

### Events

| field             | description |
|-------------------|-------------|
| `contig_id`       | Contig containing the event |
| `event_type`      | `insertion` (sequence in the contig, not the ref) or `deletion` (ref sequence missing from the contig) |
| `ref_start`, `ref_end` | Ref position of the event. For insertions, both are the ref position the insertion follows. |
| `length`          | Event length (bp) |
| `contig_position` | 0-based position of the event in the contig, in its original orientation |

### Per-alignment summary

`summary.tsv` includes these fields:

| field                | description                                                             |
|----------------------|-------------------------------------------------------------------------|
| `contig_id`          | Contig ID from plasmid reconstruction                                   |
| `ref_id`             | Ref plasmid ID                                                          |
| `ref_length`         | Length of ref plasmid                                                   |
| `reverse_complement` | Aligned sequence is reverse-complement of original contig seq           |
| `is_mapped`          | Contig is mapped somewhere against the ref plasmid                      |
| `is_supplementary`   | Alignment is a supplementary (split) alignment of the contig            |
| `spans_origin`       | Alignment wraps around the origin of the (circular) ref plasmid. If so, `alignment_end` < `alignment_start` |
| `mapq`               | Mapping quality. 0 indicates the contig could equally well align elsewhere (e.g. a repeat such as an IS element) |
| `contig_length`      | Length of the contig                                                    |
| `contig_fraction_aligned` | Fraction of the contig covered by any of its alignments that passed QC |
| `alignment_start`    | Location on the ref plasmid where the contig alignment starts           |
| `alignment_end`      | Location on the ref plasmid where the contig alignment ends             |
| `alignment_length`   | Length of the alignment                                                 |
| `percent_identity`   | Percent identity of the alignment (matches / alignment columns, including indels) |
| `included_in_output` | The alignment is included in the output complete plasmid seq            |

## Circular references

The reference plasmid is treated as circular. Each contig is aligned against both
the reference and a copy rotated by half its length, and whichever alignment covers
more of the contig in fewer pieces is kept. This means a contig spanning the
reference origin is placed as a single alignment instead of being split by
minimap2 (which would otherwise treat the reference as linear and can leave bases
at the junction unaligned). All coordinates are reported relative to the original
reference.

The reference is deliberately not doubled: a doubled reference makes every
alignment ambiguous (mapping quality 0), hiding genuinely repetitive contigs.

## Scaffold sequence

`scaffold.fa` starts at the reference origin and is built by placing each
aligned contig segment in order of its start position on the reference. Alignments
that span the origin are split there, with the part past the origin placed at the
start of the output. Reference regions not covered by any contig are
filled with `N` (one per uncovered reference base). Where two alignments overlap on
the reference (e.g. the k-mer overlap between adjacent assembly contigs), the
overlapping part of the later one is trimmed. Alignments entirely contained within
an already-covered region are excluded (`included_in_output = False`).

## Development

```
conda env create -f environment.yml
conda activate plasmid-surgery
pip install -e .
pytest
```

