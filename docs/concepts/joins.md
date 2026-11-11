# Composition and joins

A reconstruction can hold sequence from more than one plasmid: two plasmids in one bin, or one contig
joining two plasmids. A single contig like that may be a real **cointegrate** (two plasmids fused into
one molecule, often through an IS element such as IS26), or a **misassembly**: the assembler joined
two separate plasmids through a repeat both carry. The references assessed show where the switch from
one plasmid's sequence to another's happens, and the sample's long reads show whether the join is
real.

## Composition

`composition.tsv` picks the references one at a time by how much of the reconstruction each *newly*
explains: first the reference whose placements cover the most of it, then the one covering the most of
what's left, and so on, until no reference adds at least `--min-join-block` bp (default 1,000). The
last row is what none of them explain. Alignments shorter than `--min-join-block` don't count, so
scattered hits of shared IS elements don't make a reference look like part of the reconstruction.

For a reconstruction joining plasmid A and plasmid B, the first two rows are a relative of each, and
together explain nearly all of it, while each alone explains only its own part. The report's "What the
references explain" view shows the same thing along the reconstruction.

## Joins

Each stretch of the reconstruction is labelled with the first reference (in composition order) that
explains it, or none. A **join** is where the label changes between stretches of at least
`--min-join-block` bp:

| Type | Where |
|---|---|
| `reference_switch` | From one reference's sequence to another's |
| `unexplained_boundary` | Into or out of sequence no reference explains |
| `circular_junction` | A circular contig's own end-to-start junction (contigs with `circular=true`, `circular=Y` or `[topology=circular]` in the header). Checked whether or not the label changes there |

A join is a region, not a point. Where the two references' alignments overlap at the switch, as at an
IS element both plasmids carry, it covers the overlap. It is then widened to cover any repeat of the
reconstruction touching it, found by aligning the reconstruction to itself. That matters: a misjoin
happens through such a repeat, and reads from *either* plasmid run straight through the repeat into
the next sequence. Only reads anchored in unique sequence on both sides show that the join itself is
real.

## Long-read evidence

With `--long-reads`, the reads are aligned to the reconstruction's contigs. Circular contigs are
extended by their own start, so reads can run across the circular junction. At each join:

- **Spanning**: reads with one continuous alignment covering the join plus `--join-anchor` bp
  (default 1,000) of sequence on each side.
- **Clipped**: reads whose alignment ends (or starts) at the join, with at least 200 bp of the read
  left unaligned beyond it. The read continues somewhere else, which is what reads from two separate
  plasmids do at a misjoin.
- **Depth either side**: median depth in the 5 kb either side of the join, from short and long reads
  separately. Two plasmids at different copy numbers, joined by mistake, show a step.

| Call | Meaning |
|---|---|
| `supported` | At least 3 spanning reads: consistent with one molecule (e.g. a cointegrate) |
| `contradicted` | No spanning reads, and at least 3 clipped: suggests a misassembly, or a cointegrate that is rare in the culture |
| `weak` | Anything else, e.g. a repeat at the join longer than most reads |
| `no_long_reads` | No `--long-reads` given |

Cointegrates and their separate plasmids often coexist in one culture, so a join can have both spanning
and clipped reads. The counts are in `joins.tsv`.

Short reads can't span an IS-sized repeat, so they only contribute depth. Reads are aligned to the
reconstruction alone, so reads from elsewhere in the genome that share sequence with it (e.g. the
chromosome's copies of an IS) pile onto that sequence. They don't carry the unique flanks, so they don't
count as spanning, and the depth either side is taken outside the join.

## Outputs

- `composition.tsv` and `joins.tsv`, at the top level of each run (see [Outputs](../reference/outputs.md)).
- With reads: `reconstruction_short_reads_depth.bedgraph.gz` and `reconstruction_long_reads_depth.bedgraph.gz`,
  depth along each reconstruction contig.
- In the report's "What the references explain" view: a joins row (coloured by call; details in the
  tooltip), read-depth tracks along the contigs, and the composition and joins tables.
