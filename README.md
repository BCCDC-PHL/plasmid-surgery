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
                       --ref-plasmid REF_PLASMID [--assembly ASSEMBLY]
                       [--bins FASTA [FASTA ...]]
                       [--min-missed-fraction MIN_MISSED_FRACTION]
                       [--ref-features FILE [FILE ...]]
                       [--features FILE [FILE ...]]
                       [--min-feature-coverage MIN_FEATURE_COVERAGE]
                       [-o OUTDIR] [-p PREFIX] [--tmpdir TMPDIR]
                       [--min-indel-length MIN_INDEL_LENGTH]
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
  --assembly ASSEMBLY   Full assembly the reconstruction was made from.
                        Reports where the rest of the ref plasmid's sequence
                        is in the assembly (contig_sources.tsv,
                        missed_contigs.fa)
  --bins FASTA [FASTA ...]
                        Other bins of the reconstruction tool (e.g.
                        mob_recon's plasmid_*.fasta and chromosome.fasta),
                        used to label assembly contigs by bin. Requires
                        --assembly
  --min-missed-fraction MIN_MISSED_FRACTION
                        With --assembly: minimum fraction of an assembly
                        contig aligned to the ref (at a single location) to
                        call it 'missed' rather than 'shared' (default: 0.8)
  --ref-features FILE [FILE ...]
                        Feature (gene) annotations on the ref plasmid:
                        hAMRonization TSV, GFF3 or BED. Reports how completely
                        the reconstruction (and assembly) contains each
                        (ref_features.tsv)
  --features FILE [FILE ...]
                        Feature (gene) annotations on the sample's assembly or
                        reconstruction contigs: hAMRonization TSV, GFF3 or
                        BED. Reports, for features on the reconstruction (and
                        with --assembly, on any contig aligning to the ref),
                        where each is and whether it's on the ref
                        (features.tsv)
  --min-feature-coverage MIN_FEATURE_COVERAGE
                        Minimum fraction of a feature covered to call it
                        present / on the ref (default: 0.9)
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
                        an event. Deletions at least this long also count as
                        uncovered ref sequence (gaps) (default: 50)
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

With the full assembly, and mob_recon's other outputs as bins (see [Whole-assembly mode](#whole-assembly-mode)):

```
plasmid-surgery \
  --plasmid-reconstruction mob_recon/sample1/plasmid_AA123.fasta \
  --ref-plasmid CP012345.fa \
  --assembly sample1.fasta \
  --bins mob_recon/sample1/*.fasta \
  -o results/sample1 \
  -p sample1
```

With gene annotations of the ref plasmid and of the assembly (see [Feature annotations](#feature-annotations)):

```
plasmid-surgery \
  --plasmid-reconstruction mob_recon/sample1/plasmid_AA123.fasta \
  --ref-plasmid CP012345.fa \
  --assembly sample1.fasta \
  --bins mob_recon/sample1/*.fasta \
  --ref-features CP012345.ncbi.hamronized.tsv CP012345.plasmidfinder.hamronized.tsv \
  --features sample1.ncbi.hamronized.tsv sample1.plasmidfinder.hamronized.tsv \
  -o results/sample1 \
  -p sample1
```

## Outputs

All outputs are written to `--outdir` (default: current directory), with filenames
prefixed by `<prefix>_` if `--prefix` is given.

| file             | description |
|------------------|-------------|
| `assessment.tsv` | One-row summary of how well the reconstruction matches the ref (see below). Rows from many samples can be concatenated. |
| `gaps.bed`       | BED4 of ref regions not covered by any contig, named `gap_<n>_<type>` where type is `repeat` or `unique` (see [Repeats](#repeats)). A gap spanning the origin is written as two records with the same name. |
| `events.tsv`     | Insertions and deletions of at least `--min-indel-length` bp within contig alignments. |
| `unaligned.fa`   | Contig regions of at least `--min-unaligned-length` bp that don't align to the ref (whole unaligned contigs, or unaligned contig ends). Headers are `<contig_id>:<start>-<end>`, 1-based inclusive. |
| `summary.tsv`    | One row per contig alignment (see below). |
| `scaffold.fa`    | Reference-ordered sequence built from the aligned contigs, with inferred repeat copies in lowercase (see [Scaffold sequence](#scaffold-sequence)). Header is the sample ID. |
| `contig_sources.tsv` | With `--assembly` only: every assembly contig that aligns to the ref, and where it is (see [Whole-assembly mode](#whole-assembly-mode)). |
| `missed_contigs.fa`  | With `--assembly` only: assembly contigs classified `missed`, with `bin=<label>` in the header when `--bins` is given. |
| `ref_features.tsv`   | With `--ref-features` only: each feature annotated on the ref, and how completely the reconstruction (and assembly) contains it. |
| `features.tsv`       | With `--features` only: each feature on the reconstruction's contigs (and with `--assembly`, on any contig aligning to the ref), with its contig's classification and bin, and whether the feature itself is on the ref. |

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
| `num_ambiguous_contigs`      | Contigs with at least one alternative placement (see [Repeats](#repeats)) |

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
| `is_secondary`       | Alignment is a secondary alignment: an alternative placement of the contig if `included_in_output` is True, otherwise a weaker alignment that was not used |
| `num_placements`     | Number of places on the ref this contig could be placed (1 for unique contigs) |
| `spans_origin`       | Alignment wraps around the origin of the (circular) ref plasmid. If so, `alignment_end` < `alignment_start` |
| `mapq`               | Mapping quality. 0 indicates the contig could equally well align elsewhere (e.g. a repeat such as an IS element) |
| `contig_length`      | Length of the contig                                                    |
| `contig_fraction_aligned` | Fraction of the contig covered by any of its alignments that passed QC |
| `alignment_start`    | Location on the ref plasmid where the contig alignment starts           |
| `alignment_end`      | Location on the ref plasmid where the contig alignment ends             |
| `alignment_length`   | Length of the alignment                                                 |
| `percent_identity`   | Percent identity of the alignment (matches / alignment columns, including indels) |
| `included_in_output` | The alignment is included in the output complete plasmid seq            |

## Whole-assembly mode

A reconstruction can miss parts of a plasmid because they weren't assembled, or because the
reconstruction tool put contigs that were assembled in the wrong bin (another plasmid, or the
chromosome). Given the full assembly with `--assembly`, plasmid-surgery separates the two, and
lists every assembly contig that aligns to the ref in `contig_sources.tsv`:

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

`--bins` takes any number of FASTA files, e.g. mob_recon's `plasmid_*.fasta` and
`chromosome.fasta`, or another tool's per-plasmid outputs (e.g. TaDReP's `*-contigs.fna`). Contigs are matched between the
reconstruction, the assembly and the bins by sequence (in either orientation), not by name, so
tools that rename or reorient contigs work too.

A carbapenemase gene on a contig that mob_recon calls chromosomal is a common puzzle. If that
contig is classified `missed` against a plausible ref plasmid, its `depth_relative_to_reconstruction`
is close to 1, and it is at well above chromosomal depth, that is good evidence it's a misbinned
plasmid contig. If it's `shared` (only part of it matches the plasmid), a chromosomal location is
more plausible. This is evidence, not proof: depth is noisy for short contigs, and the
conclusion is only as good as the ref.

## Feature annotations

Gene (or any feature) annotations can be given for the ref (`--ref-features`) and for the sample's
assembly or reconstruction (`--features`), in any of these formats (detected automatically; any
number of files):

- [hAMRonization](https://github.com/pha4ge/hAMRonization) TSV: one format for many tools, e.g.
  abricate with any database (ncbi, plasmidfinder, vfdb, ...), AMRFinderPlus, RGI, ResFinder.
  Note that hAMRonization's AMRFinderPlus parser drops virulence hits.
- GFF3, e.g. from Bakta, Prokka or NCBI. Features that are only parents of other features (a
  `gene` with a `CDS` child) are reported once.
- BED: `seq_id`, `start`, `end`, and optionally `name` (column 4) and `strand` (column 6).

No feature class is treated specially: names, products and categories (drug class for
hAMRonization, feature type for GFF3) are passed through. Features are matched to sequences by ID
(ignoring a `.N` version suffix), or by sequence when IDs differ: from a GFF3 `##FASTA` section
(Bakta renames contigs by default), or, for `--features`, from `--assembly` (e.g. annotations of
the assembly, with a TaDReP reconstruction whose contigs are renamed and reoriented).

`ref_features.tsv` (one row per ref feature; coordinates 0-based, half-open):

| field                     | description |
|---------------------------|-------------|
| `name`, `product`, `category`, `source` | From the annotation |
| `ref_start`, `ref_end`, `strand`, `length` | Position on the ref |
| `reconstruction_fraction` | Fraction of the feature covered by the reconstruction |
| `inferred_fraction`       | Further fraction covered only by alternative placements (repeat copies; lowercase in the scaffold) |
| `reconstruction_status`   | `present` (≥ `--min-feature-coverage`, default 0.9), `inferred` (present only counting alternative placements), `partial`, or `absent` |
| `assembly_fraction`, `assembly_status` | The same for the whole assembly (with `--assembly`) |
| `found_in`                | Where the assembly's copy is: `reconstruction` and/or bin labels (`-` if nowhere) |

`features.tsv` (one row per sample feature on the reconstruction or, with `--assembly`, on a contig
aligning to the ref; features on other contigs aren't reported):

| field                     | description |
|---------------------------|-------------|
| `name`, `product`, `category`, `source` | From the annotation |
| `contig_id`, `contig_start`, `contig_end`, `strand`, `length` | Position on the contig |
| `contig_classification`, `bin` | As in `contig_sources.tsv` (`in_reconstruction` without `--assembly`) |
| `fraction_on_ref`         | Fraction of the feature's bases matched to ref bases. Bases in an insertion relative to the ref don't count, so a gene inside novel sequence is `not_on_ref` even on a contig that otherwise aligns |
| `ref_status`              | `on_ref`, `partial` or `not_on_ref` |
| `ref_start`, `ref_end`    | Where the feature is on the ref, if both its ends are matched to it |

Together they answer "is this gene where the ref has it?" from both sides. A ref gene `absent`
from the reconstruction but `found_in` another bin, or a sample gene on a `missed` contig, points
at a binning error. A ref gene `partial` or `absent` in the whole assembly means the assembly
itself doesn't have it complete (e.g. a repeat collapsed into fragments too short to call). A
sample gene `not_on_ref` is sequence the ref doesn't have.

## Reconstructions from other tools

Any FASTA of contigs works as `--plasmid-reconstruction`, e.g. a mob_recon `plasmid_*.fasta` or a
TaDReP `<sample>-<plasmid>-contigs.fna` (not `-pseudo.fna`, whose contigs are joined with fixed-length
N gaps). Contigs whose sequence (in either orientation) repeats an earlier one are ignored with a
warning: TaDReP, for example, writes a contig once per alignment hit.

## Coverage and deletions

Ref coverage (and so gaps, and feature coverage) counts the ref bases each alignment spans,
except deletions within an alignment of at least `--min-indel-length` (default 50) bp: ref
sequence the contig doesn't have. Such a deletion is reported both in `events.tsv` and as a
gap. Shorter deletions count as covered, so that small indels don't fragment coverage.

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

## Repeats

Short-read assemblies typically collapse repeats (e.g. IS26, Tn4401) into a single
contig, which aligns equally well to each copy on the ref. minimap2 picks one copy for
the primary alignment, and the others are reported as secondary alignments.

A secondary alignment is treated as an *alternative placement* of the contig if its
alignment score is at least 95% of the primary's. A unique contig that merely contains
part of a repeat has a much weaker secondary alignment, so isn't treated as
ambiguous. Alternative placements are used to:

- type each gap as `repeat` (an alternative placement covers it) or `unique`
  (nothing covers it). A gap partly covered by an alternative placement is split
  into a repeat part and a unique part. Unique gaps are the ones that indicate
  sequence genuinely missing from the reconstruction.
- fill repeat gaps in `scaffold.fa` (see below).

Alternative placements don't count towards `ref_bases_covered` or
`contig_bases_aligned`.

## Scaffold sequence

`scaffold.fa` starts at the reference origin and is built by placing each
aligned contig segment in order of its start position on the reference. Alignments
that span the origin are split there, with the part past the origin placed at the
start of the output. Reference regions not covered by any contig's placement are
filled from alternative placements where possible (i.e. repeat gaps), in
**lowercase** to mark them as inferred, and otherwise with `N` (one per uncovered
reference base). Where two alignments overlap on
the reference (e.g. the k-mer overlap between adjacent assembly contigs), the
overlapping part of the later one is trimmed. Alignments entirely contained within
an already-covered region are excluded (`included_in_output = False`).

Lowercase follows the common soft-masking convention. Note that some tools treat
lowercase as masked (e.g. BLAST with `-lcase_masking`) and others silently
uppercase input, so inferred regions may be ignored or indistinguishable
downstream.

## Development

```
conda env create -f environment.yml
conda activate plasmid-surgery
pip install -e .
pytest
```

