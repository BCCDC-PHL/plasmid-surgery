# The HTML reports

There are two: a report per run (one sample, against each of its references), and a
[cohort report](#the-cohort-report) (many samples, against one reference).

## The per-run report

`plasmid-surgery --report` (or `plasmid-surgery-plot <outdir>` on a finished run) writes
`report.html` to the output directory: one self-contained page with the results against every
reference, and a `map.svg` per reference in its directory, for slides and papers.

**[Open the example report](../example-report/report.html)**: one reconstruction against two
references, with the assembly, bins and gene annotations.

The page has no external dependencies, so it opens offline, and can be emailed or archived with the
results. It's drawn entirely from the files the run wrote.

### What's on it

1. **References**, ranked by score, with the headline metrics for each.
2. **Candidate references**, when they were found in a plasmid database (`--db-sketch`).
3. **Where each contig is**: a matrix of the reconstruction's contigs against the references, with
   each contig's position on each (or "not aligned"), and with reads, each contig's read depth.
4. **What the references explain**: the reconstruction on its own coordinates, its contigs end to end,
   with a row per reference showing which parts of each contig align to it, and a row for the parts
   that align to none of them (red). The maps below are drawn on each reference, where sequence that
   aligns to no reference has no place; this view shows how much of the reconstruction that is, and
   where. When the references together explain clearly more than the best one, the reconstruction
   may combine sequence from more than one plasmid (two plasmids in one bin, or a hybrid contig
   joining two plasmids). To find what the unexplained parts are, search a plasmid database with a
   reference's `unaligned.fa` (`plasmid-surgery-refs --plasmid-reconstruction <ref>/unaligned.fa ...`).
   With `--features`, the sample's genes are shown too, in red if mostly in an unexplained part.
   Joins between parts explained by different references (and circular contigs' own junctions) are
   marked, coloured by their long-read evidence (green supported, red contradicted, yellow weak, grey
   without long reads), with read-depth tracks along the contigs and the composition and joins tables
   below. See [Composition and joins](../concepts/joins.md).
5. **Maps**: one linear map per reference, in rank order.

### Reading the maps

<figure markdown="span">
  ![The map of the example's first reference](../example-report/map.svg)
  <figcaption>The example's first reference, as saved in <code>pREF1/map.svg</code>.</figcaption>
</figure>

Each map has up to six tracks, drawn on the reference's coordinates:

| Track | Shows |
|---|---|
| reference | The reference, with gaps: unique in red, repeat hatched. ▼ marks an insertion in a contig. With reads, each gap's tooltip says whether it's in the reads |
| read depth | With `--short-reads` / `--long-reads`: mean read depth along the reference (scaled to twice the median; taller bars are clipped and darker), with stretches below `--min-read-depth` in red underneath: sequence the sample doesn't have. See [Reads](../concepts/reads.md) |
| reconstruction | The reconstruction's contigs where they align, as arrows (pointing left on the reverse strand). Alternative placements (e.g. a collapsed repeat's other copies) are dashed outlines; contigs left out of the scaffold, dotted |
| other assembly | With `--assembly`: the assembly's other contigs aligning to the reference, coloured by classification (`missed`, `shared`, `repeat`) |
| ref features | With `--ref-features`: the reference's genes, coloured by status in the reconstruction (present, inferred, partial, absent) |
| sample features | With `--features`: the sample's genes where they are on the reference; genes not on it are listed below |

Hover over anything for its details (positions, identity, depth, bin, status, ...). Hovering over a
contig highlights it everywhere: in the matrix, and in every reference's map.

### Comparing references

* **Same colours throughout.** Each reconstruction contig has one colour in every map and in the
  matrix, so it can be followed from reference to reference.
* **One scale.** By default every map uses the same bp-per-pixel scale, so references' lengths, and
  the sizes of gaps and contigs, compare directly. "Fit each reference" scales each map to the page
  width instead, for detail.
* **Labels.** Contigs are labelled by their ID minus the prefix all of them share (e.g. `21` for
  `SAMN38646027_21`), where the label fits; the full ID is in the tooltip and the matrix.

### For slides and papers

Each `map.svg` is the reference's map on its own, fitted to its width, with a title. It's plain SVG,
so it opens in a browser, and in Inkscape or Illustrator for editing; tooltips are kept as SVG
`<title>` elements, which most viewers ignore.

## The cohort report

`plasmid-surgery-cohort` puts many samples' results against one shared reference in one page: for
example, every isolate carrying a plasmid of interest, each assessed against the same reference.

```bash
plasmid-surgery-cohort results/*/ --ref CP012345.1 -o cohort.html
```

It takes plasmid-surgery output directories (typically one per sample) and uses each one's results
against the reference given by `--ref`, or by default the reference assessed in the most runs. A run
without that reference, or with a copy of a different length, is skipped, and listed at the top.

Sample names link to each sample's own `report.html` by a relative path, so the links keep working if
the cohort page and the run directories are moved or copied together. To share the cohort page on its
own, add `--embed-reports`: each sample's report is stored inside the page (adding its size) and opens
in a new tab from there.

To choose the reference, `--list-refs` prints, instead of drawing the report, how each reference fared
across the runs (from each run's top-level `assessment.tsv`), the most often near the top first:

```bash
plasmid-surgery-cohort results/*/ --list-refs | column -t
```

```
ref_id      ref_length  runs_assessed  rank_1  in_top_3  mean_score  mean_ref_fraction_covered
CP012345.1  61250       12/12          9       12        0.912       0.947
CP067890.1  58731       10/12          3       9         0.701       0.903
```

| Column | Meaning |
|---|---|
| `runs_assessed` | Runs that assessed the reference, of all runs |
| `rank_1`, `in_top_N` | Runs ranking it first, or in their top N (`--top`, default 3) |
| `mean_score` | Its mean score over all runs, a run that didn't assess it counting as 0 |
| `mean_ref_fraction_covered` | Its mean coverage over the runs that assessed it |

A reference near the top on both `in_top_N` and `mean_score` suits most of the runs. Runs found
their references separately (`--db-sketch`) may not all have assessed it: re-run those with
`--ref-plasmid` and that reference to include them.

**[Open the example cohort](../example-cohort/cohort.html)**: five simulated samples against one
reference: complete, with a contig binned as chromosome, with a segment missing from the assembly,
with a deletion and an insertion, and with an extra gene.

### What's on it

1. **Samples**, with the headline metrics. Click a column to sort by it. Sample names link to their own
   reports, if they were made with `--report`.
2. **Coverage of the reference**: a strip per sample, along the reference (with its genes above, if
   annotated), coloured by what happened to each position:

    | Colour | Meaning |
    |---|---|
    | blue | covered by the sample's reconstruction |
    | hatched orange | a repeat gap: covered only by a repeat contig's alternative placement (inferred) |
    | orange | not in the reconstruction, but in a `missed` contig in the assembly: lost by binning |
    | lilac | not in the reconstruction; only in a `shared` or `repeat` contig |
    | red | not in the reconstruction, or the assembly (without `--assembly`: any unique gap) |

    ▼ marks an insertion. The assembly categories need the runs to have used `--assembly`.

3. **Reference genes**: each of the reference's genes (in reference order) by its status in each
   sample's reconstruction. ● marks a gene missing from the reconstruction but present in the assembly
   (binning). Needs `--ref-features`.
4. **Genes not on the reference**: the samples' genes that aren't on the reference, by sample. Needs
   `--features`.
5. **Reconstructions compared**: every sample's reconstruction and the reference, aligned to each
   other all-vs-all (one minimap2 call), independent of the reference's coordinates. Each cell is how
   much of the row's sequence the column contains; hover for both directions and identity. A
   fragmentary reconstruction contained in a complete plasmid is dark along its row and light down its
   column; sequence a sample has that others lack (e.g. an acquired gene) lightens its row. Ordered by
   average-linkage clustering, so similar plasmids sit together. Every pair is also written to
   `<report name>_pairs.tsv` beside the report (`a_in_b`, `b_in_a`, `identity`, aligned bases).
   It uses each run's `reconstruction.fa` (written by plasmid-surgery since this was added: older runs
   are listed and left out; rerun them to include them). Alignments count if at least 500 bp and 90%
   identical (gap-compressed); indels of 50 bp or more within an alignment don't count as contained.
   `--no-compare` skips it; `-t` sets minimap2's threads.

   Below it, **SNP distances** between every pair, from the same alignments:

    - **One-to-one**: each stretch of either reconstruction is compared once (longest alignments first),
      and stretches that align more than once between the two (repeats, e.g. an IS that one assembly
      put in a contig of its own) are left out, so repeat copies aren't compared with each other.
    - **Trimmed**: differences within 100 bp of an alignment's or a contig's end are ignored; assembly
      errors concentrate at contig ends.
    - **SNPs and indels separately**: long-read assembly errors are mostly indels (homopolymers), so
      the distance counts substitutions only; indels are in the tooltip and the table.
    - **Clusters**: SNPs with at least 2 others within 1 kb are counted apart, as a block from another
      lineage brings many at once. Cells show the SNPs outside clusters.
    - **How much was compared**: every count comes with the bases compared. A pair is grey, with no
      distance, when less than half of the shorter of the two was compared (e.g. a small fragment
      that shares only part of its sequence).

   In `<report name>_pairs.tsv`: `compared_bases`, `compared_fraction_a`/`_b`, `snps`,
   `snps_in_clusters`, `snps_outside_clusters`, `snps_per_100kb` (outside clusters), `indels`, and
   `snp_distance_reported`. Distances are pairwise: each pair is compared over the sequence the two
   share, so pairs aren't all compared over the same stretch.

The first four sections can be ordered together, by coverage (the default), score, sample or input
order. To compare samples, use coverage: the score ranks references within one sample's run, and with
`--assembly` it uses the assembly's coverage, so it doesn't reflect sequence lost by binning.
