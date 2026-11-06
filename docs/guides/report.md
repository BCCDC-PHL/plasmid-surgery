# The HTML report

`plasmid-surgery --report` (or `plasmid-surgery-plot <outdir>` on a finished run) writes
`report.html` to the output directory: one self-contained page with the results against every
reference, and a `map.svg` per reference in its directory, for slides and papers.

**[Open the example report](../example-report/report.html)**: one reconstruction against two
references, with the assembly, bins and gene annotations.

The page has no external dependencies, so it opens offline, and can be emailed or archived with the
results. It's drawn entirely from the files the run wrote.

## What's on it

1. **References**, ranked by score, with the headline metrics for each.
2. **Candidate references**, when they were found in a plasmid database (`--db-sketch`).
3. **Where each contig is**: a matrix of the reconstruction's contigs against the references, with
   each contig's position on each (or "not aligned").
4. **Maps**: one linear map per reference, in rank order.

## Reading the maps

<figure markdown="span">
  ![The map of the example's first reference](../example-report/map.svg)
  <figcaption>The example's first reference, as saved in <code>pREF1/map.svg</code>.</figcaption>
</figure>

Each map has up to five tracks, drawn on the reference's coordinates:

| Track | Shows |
|---|---|
| reference | The reference, with gaps: unique in red, repeat hatched. ▼ marks an insertion in a contig |
| reconstruction | The reconstruction's contigs where they align, as arrows (pointing left on the reverse strand). Alternative placements (e.g. a collapsed repeat's other copies) are dashed outlines; contigs left out of the scaffold, dotted |
| other assembly | With `--assembly`: the assembly's other contigs aligning to the reference, coloured by classification (`missed`, `shared`, `repeat`) |
| ref features | With `--ref-features`: the reference's genes, coloured by status in the reconstruction (present, inferred, partial, absent) |
| sample features | With `--features`: the sample's genes where they are on the reference; genes not on it are listed below |

Hover over anything for its details (positions, identity, depth, bin, status, ...). Hovering over a
contig highlights it everywhere: in the matrix, and in every reference's map.

## Comparing references

* **Same colours throughout.** Each reconstruction contig has one colour in every map and in the
  matrix, so it can be followed from reference to reference.
* **One scale.** By default every map uses the same bp-per-pixel scale, so references' lengths, and
  the sizes of gaps and contigs, compare directly. "Fit each reference" scales each map to the page
  width instead, for detail.
* **Labels.** Contigs are labelled by their ID minus the prefix all of them share (e.g. `21` for
  `SAMN38646027_21`), where the label fits; the full ID is in the tooltip and the matrix.

## For slides and papers

Each `map.svg` is the reference's map on its own, fitted to its width, with a title. It's plain SVG,
so it opens in a browser, and in Inkscape or Illustrator for editing; tooltips are kept as SVG
`<title>` elements, which most viewers ignore.
