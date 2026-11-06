# Feature annotations

Gene annotations, of any kind (resistance, virulence, replicons, mobility, ...), can be given for the
reference (`--ref-features`) and for the sample's assembly or reconstruction (`--features`). No class
of gene is treated specially: names, products and categories are passed through as given.

## Formats

Any number of files, in any of these formats (detected automatically):

* [hAMRonization](https://github.com/pha4ge/hAMRonization) TSV, one format for the output of many
  annotation tools. An empty file (no header) counts as no features.
* GFF3. Features that are only parents of other features (a `gene` with a `CDS` child) are reported
  once; a `##FASTA` section is read, if present.
* BED: `seq_id`, `start`, `end`, and optionally `name` (column 4) and `strand` (column 6).

Features are matched to sequences by ID, ignoring a `.N` version suffix, or by sequence where IDs
differ: through a GFF3 `##FASTA` section, or, for `--features`, through `--assembly` (so annotations
of the assembly work with a reconstruction whose contigs were renamed or reverse-complemented).

## Reference features

`ref_features.tsv` reports, for each feature on the reference, the fraction covered by the
reconstruction, and the further fraction covered only by alternative placements, with a status:

| Status | Meaning |
|---|---|
| `present` | Covered to at least `--min-feature-coverage` (default 0.9) |
| `inferred` | Present only counting alternative placements (e.g. a gene in a collapsed repeat) |
| `partial` | Partly covered |
| `absent` | Not covered |

With `--assembly`, the same for the whole assembly, and `found_in`: where the assembly's copy is.

## Sample features

`features.tsv` reports each feature on the reconstruction or (with `--assembly`) on any contig
aligning to the reference: its contig's classification and bin, and `fraction_on_ref`, the fraction
of its bases matched to reference bases. Bases in an insertion don't count, so a gene inside novel
sequence is `not_on_ref` even on a contig that otherwise aligns. Where both ends are on the reference,
`ref_start` and `ref_end` give its position there. Features on contigs that don't align to the
reference aren't reported.

## Two complementary views

Sample annotations can only report genes assembled well enough to be called. A gene split across
short contigs, for example, may not be called at all, but the reference's annotation still shows it
as `partial` or `absent`. Reference annotations can only report genes the reference has; the
sample's show the rest as `not_on_ref`. See [Feature annotations](../behaviour/features.md).
