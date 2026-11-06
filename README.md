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
                       --ref-plasmid REF_PLASMID [--tmpdir TMPDIR]
                       [--summary-output SUMMARY_OUTPUT] [--output OUTPUT]
                       [--output-seq-id OUTPUT_SEQ_ID]
                       [--output-bases-per-line OUTPUT_BASES_PER_LINE]
                       [--log-level LOG_LEVEL]

options:
  -h, --help            show this help message and exit
  --plasmid-reconstruction PLASMID_RECONSTRUCTION
                        Plasmid reconstruction to be repaired.
  --ref-plasmid REF_PLASMID
                        Reference plasmid to be used for plasmid surgery.
  --tmpdir TMPDIR       temp dir to use for minimap2 alignments (default:
                        './plasmid-surgery_tmp')
  --summary-output SUMMARY_OUTPUT
                        Summary file for plasmid surgery process. (default:
                        './plasmid-surgery-summary.tsv')
  --output OUTPUT       Path to output file for repaired plasmid
                        reconstruction (default: './plasmid-surgery-
                        output.fa')
  --output-seq-id OUTPUT_SEQ_ID
                        Sequence ID to use in the output fasta header. If not
                        supplied, uses plasmid reconstruction filename (minus
                        file extension)
  --output-bases-per-line OUTPUT_BASES_PER_LINE
                        Number of bases per line in the output fasta (default:
                        70)
  --log-level LOG_LEVEL
                        Log level (default: 'info')
```

Example:

```
plasmid-surgery \
  --plasmid-reconstruction your_plasmid_reconstruction.fa \
  --ref-plasmid your_ref_plasmid.fa \
  --output fixed_plasmid_reconstruction.fa \
  --summary-output plasmid-surgery-summary.tsv 
```

## Outputs

The plasmid surgery summary tsv file includes these fields:

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

## Output sequence

The output sequence starts at the reference origin and is built by placing each
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

