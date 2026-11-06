# Outputs

All outputs are written to `--outdir` (default: the current directory), named `<prefix>_<file>` if
`--prefix` is given. Positions are 0-based, half-open (BED-style) unless stated otherwise.

| file             | description |
|------------------|-------------|
| `assessment.tsv` | One-row summary of how well the reconstruction matches the ref (see below). Rows from many samples can be concatenated. |
| `gaps.bed`       | BED4 of ref regions not covered by any contig, named `gap_<n>_<type>` where type is `repeat` or `unique` (see [Coverage, gaps and events](../concepts/coverage.md)). A gap spanning the origin is written as two records with the same name. |
| `events.tsv`     | Insertions and deletions of at least `--min-indel-length` bp within contig alignments. |
| `unaligned.fa`   | Contig regions of at least `--min-unaligned-length` bp that don't align to the ref (whole unaligned contigs, or unaligned contig ends). Headers are `<contig_id>:<start>-<end>`, 1-based inclusive. |
| `summary.tsv`    | One row per contig alignment (see below). |
| `scaffold.fa`    | Reference-ordered sequence built from the aligned contigs, with inferred repeat copies in lowercase (see [The scaffold](../concepts/scaffold.md)). Header is the sample ID. |
| `contig_sources.tsv` | With `--assembly` only: every assembly contig that aligns to the ref, and where it is (see [Whole-assembly mode](../concepts/whole-assembly.md)). |
| `missed_contigs.fa`  | With `--assembly` only: assembly contigs classified `missed`, with `bin=<label>` in the header when `--bins` is given. |
| `ref_features.tsv`   | With `--ref-features` only: each feature annotated on the ref, and how completely the reconstruction (and assembly) contains it. |
| `features.tsv`       | With `--features` only: each feature on the reconstruction's contigs (and with `--assembly`, on any contig aligning to the ref), with its contig's classification and bin, and whether the feature itself is on the ref. |

Intermediate alignment files go to a temporary directory that's removed after the run, unless
`--tmpdir` is given, in which case they're kept there.

## `assessment.tsv`

One row of metrics. Only alignments that passed the [quality filters](../concepts/alignment.md#quality-filters)
count. Columns for options not used are `NA`, so files from different runs can be concatenated.

| field                        | description |
|------------------------------|-------------|
| `sample_id`                  | `--sample-id` if given, else `--prefix` if given, else the reconstruction filename minus extension |
| `ref_id`, `ref_length`       | Ref plasmid ID and length |
| `num_contigs`, `total_contig_length` | Number and total length of contigs in the reconstruction |
| `ref_bases_covered`          | Ref bases covered by at least one contig alignment |
| `ref_fraction_covered`       | `ref_bases_covered / ref_length` |
| `ref_bases_covered_multiple` | Ref bases covered by more than one alignment (contig overlaps, or duplicated sequence) |
| `ref_fraction_covered_including_repeats` | As `ref_fraction_covered`, also counting repeat gaps as covered |
| `num_gaps`, `total_gap_length`, `largest_gap_length` | Uncovered ref regions (repeat and unique) |
| `num_repeat_gaps`, `repeat_gap_length` | Uncovered ref regions where an ambiguous contig could be placed: likely a collapsed repeat, present in the reconstruction but not placed there |
| `num_unique_gaps`, `unique_gap_length`, `largest_unique_gap_length` | Uncovered ref regions no contig could be placed at: likely genuinely missing from the reconstruction |
| `contig_bases_aligned`       | Contig bases within an alignment to the ref (including inserted bases inside alignments; see `num_large_insertions`) |
| `contig_fraction_aligned`    | `contig_bases_aligned / total_contig_length`. Sequence not accounted for by the ref is in the unaligned output. |
| `num_unaligned_contigs`      | Contigs with no alignment to the ref that passed QC: possible chromosomal contamination or a different plasmid |
| `num_unaligned_regions`, `unaligned_region_bases` | Regions written to the unaligned output |
| `num_large_insertions`, `num_large_deletions` | Events written to the events output |
| `num_ambiguous_alignments`   | Placement alignments with mapping quality 0: the contig aligns equally well elsewhere on the ref (e.g. a repeated IS element) |
| `num_ambiguous_contigs`      | Contigs with at least one alternative placement (see [Alignment and placement](../concepts/alignment.md)) |

With `--assembly` (otherwise `NA`):

| field                        | description |
|------------------------------|-------------|
| `assembly_ref_bases_covered`, `assembly_ref_fraction_covered` | Ref coverage by the whole assembly, as `ref_bases_covered` / `ref_fraction_covered` are for the reconstruction |
| `assembly_ref_fraction_covered_including_repeats` | As above, also counting repeat gaps as covered |
| `assembly_unique_gap_length` | Ref bases no assembly contig could be placed at: missing from the assembly itself, not just from the reconstruction |
| `missed_ref_bases`           | Ref bases covered by `missed` contigs but not by the reconstruction: what the reconstruction lost by binning |
| `num_missed_contigs`, `missed_contig_bases` | Contigs classified `missed` (written to `missed_contigs.fa`) |
| `num_shared_contigs`, `num_repeat_contigs` | Other contigs classified `shared` or `repeat` |
| `num_reconstruction_contigs_not_in_assembly` | Reconstruction contigs whose sequence isn't in `--assembly`. Should be 0; otherwise the reconstruction probably wasn't made from this assembly |

With `--ref-features` / `--features` (otherwise `NA`):

| field                        | description |
|------------------------------|-------------|
| `num_ref_features`           | Features annotated on the ref |
| `num_ref_features_present`, `_partial`, `_absent`, `_inferred` | By `reconstruction_status` in `ref_features.tsv` |
| `num_ref_features_missed`    | Ref features not present in the reconstruction but present in the assembly (with `--assembly`) |
| `num_features`               | Sample features reported in `features.tsv` |
| `num_features_not_on_ref`    | Of those, features not on the ref: sequence the ref doesn't have (e.g. an acquired gene, or a gene on a `shared` contig's other part) |
| `num_features_on_missed_contigs` | Of those, features on `missed` contigs |

## `summary.tsv`

One row per alignment, including those that failed the quality filters, and unmapped contigs.

| field                | description                                                             |
|----------------------|-------------------------------------------------------------------------|
| `contig_id`          | Contig ID from plasmid reconstruction                                   |
| `ref_id`             | Ref plasmid ID                                                          |
| `ref_length`         | Length of ref plasmid                                                   |
| `reverse_complement` | Aligned sequence is reverse-complement of original contig seq           |
| `is_mapped`          | Contig is mapped somewhere against the ref plasmid                      |
| `is_supplementary`   | Alignment is a supplementary (split) alignment of the contig            |
| `is_secondary`       | Alignment is a secondary alignment: an alternative placement of the contig if `included_in_output` is True, otherwise a weaker alignment that was not used |
| `num_placements`     | Number of places on the ref this contig could be placed (1 for unique contigs) |
| `spans_origin`       | Alignment wraps around the origin of the (circular) ref plasmid. If so, `alignment_end` < `alignment_start` |
| `mapq`               | Mapping quality. 0 indicates the contig could equally well align elsewhere (e.g. a repeat such as an IS element) |
| `contig_length`      | Length of the contig                                                    |
| `contig_fraction_aligned` | Fraction of the contig covered by any of its alignments that passed QC |
| `alignment_start`    | Where on the ref the alignment starts (0-based)                         |
| `alignment_end`      | Where on the ref the alignment ends (exclusive; < `alignment_start` if it spans the origin) |
| `alignment_length`   | Length of the alignment                                                 |
| `percent_identity`   | Gap-compressed percent identity: matches / (matches + mismatches + indel events) |
| `included_in_output` | The alignment was used in `scaffold.fa`                                 |

## `events.tsv`

| field             | description |
|-------------------|-------------|
| `contig_id`       | Contig containing the event |
| `event_type`      | `insertion` (sequence in the contig, not the ref) or `deletion` (ref sequence missing from the contig) |
| `ref_start`, `ref_end` | Ref position of the event. For insertions, both are the ref position the insertion follows. |
| `length`          | Event length (bp) |
| `contig_position` | 0-based position of the event in the contig, in its original orientation |

## `contig_sources.tsv`

With `--assembly`. One row per assembly contig that aligns to the reference.

| field                    | description |
|--------------------------|-------------|
| `contig_id`, `contig_length` | Assembly contig |
| `depth`                  | Read depth from the contig header (Unicycler `depth=`, SPAdes `_cov_`, `cov=`), if present |
| `depth_relative_to_reconstruction` | `depth` relative to the length-weighted mean depth of the reconstruction's contigs |
| `bin`                    | With `--bins`: the bin(s) containing the contig (FASTA filename without extension; comma-separated if more than one), or `unbinned`. `NA` without `--bins` |
| `classification`         | See below |
| `contig_fraction_on_ref` | Fraction of the contig aligned to the ref |
| `ref_bases_covered`      | Ref bases covered by the contig |
| `ref_start`, `ref_end`   | Position of the contig's primary alignment on the ref (`ref_end` < `ref_start` if it spans the origin) |
| `num_placements`, `mapq`, `percent_identity` | As in `summary.tsv` |

| classification       | meaning |
|----------------------|---------|
| `in_reconstruction`  | Part of the reconstruction being assessed |
| `shared`             | Less than `--min-missed-fraction` (default 0.8) of the contig aligns to the ref. Likely a contig from the chromosome or another plasmid sharing some sequence (e.g. a transposon) with this one: not a binning error |
| `repeat`             | Most of the contig aligns to the ref, but equally well at more than one place (e.g. an IS element), so it may belong on this plasmid, elsewhere, or both |
| `missed`             | Most of the contig aligns to the ref, at a single location: likely a plasmid contig the reconstruction tool put elsewhere |

## `ref_features.tsv`

With `--ref-features`.

| field                     | description |
|---------------------------|-------------|
| `name`, `product`, `category`, `source` | From the annotation |
| `ref_start`, `ref_end`, `strand`, `length` | Position on the ref |
| `reconstruction_fraction` | Fraction of the feature covered by the reconstruction |
| `inferred_fraction`       | Further fraction covered only by alternative placements (repeat copies; lowercase in the scaffold) |
| `reconstruction_status`   | `present` (≥ `--min-feature-coverage`, default 0.9), `inferred` (present only counting alternative placements), `partial`, or `absent` |
| `assembly_fraction`, `assembly_status` | The same for the whole assembly, counting its alternative placements for `inferred` (with `--assembly`) |
| `found_in`                | Where the assembly's copy is: `reconstruction` and/or bin labels (`-` if nowhere) |

## `features.tsv`

With `--features`.

| field                     | description |
|---------------------------|-------------|
| `name`, `product`, `category`, `source` | From the annotation |
| `contig_id`, `contig_start`, `contig_end`, `strand`, `length` | Position on the contig |
| `contig_classification`, `bin` | As in `contig_sources.tsv` (`in_reconstruction` without `--assembly`) |
| `fraction_on_ref`         | Fraction of the feature's bases matched to ref bases. Bases in an insertion relative to the ref don't count, so a gene inside novel sequence is `not_on_ref` even on a contig that otherwise aligns |
| `ref_status`              | `on_ref`, `partial` or `not_on_ref` |
| `ref_start`, `ref_end`    | Where the feature is on the ref, if both its ends are matched to it |
