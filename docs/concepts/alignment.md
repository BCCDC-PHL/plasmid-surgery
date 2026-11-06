# Alignment and placement

Each contig of the reconstruction is aligned to the reference with minimap2
(`-a -cx asm10 --eqx -N 10 -p 0.9`): the assembly-to-reference preset, with explicit match and
mismatch CIGAR operations, reporting up to 10 secondary alignments that score at least 90% of the
primary.

## A circular reference

A plasmid's origin is arbitrary, but minimap2 treats the reference as linear, so a contig crossing
the origin would come out as two alignments. Each contig is therefore aligned twice: against the
reference, and against a copy rotated by half its length. The result covering more of the contig, in
fewer pieces, is kept, and its coordinates are converted back to the original reference. An
alignment that wraps the origin has `alignment_end` < `alignment_start` and `spans_origin` = True in
`summary.tsv`.

The reference isn't doubled instead, because every contig would then align equally well to both
copies, with mapping quality 0. See [Circular references](../behaviour/origin.md).

## Quality filters

An alignment is used only if it passes all of:

| Filter | Default |
|---|---|
| Percent identity | ≥ 80% |
| Aligned query length | ≥ 50 bp |
| Fraction of the alignment that's insertions or deletions | ≤ 0.8 |

Alignments that fail are still listed in `summary.tsv`, with `included_in_output` = False, but
don't count towards any metric.

## Identity

`percent_identity` is *gap-compressed*: matches / (matches + mismatches + indel events). Each
insertion or deletion counts once, however long, so one structural difference doesn't swamp the
identity of an otherwise near-identical alignment, or push it below the identity filter. See
[Insertions, deletions and identity](../behaviour/indels.md).

## Placements and alternative placements

A contig's *placement* is its primary alignment, plus any supplementary alignments (other parts of
the contig aligning elsewhere, e.g. across a rearrangement). A secondary alignment whose score is at
least 95% of the primary's is an *alternative placement*: somewhere else the contig fits just as
well, as a collapsed repeat does at each of its copies. `num_placements` counts them (1 for a unique
contig). Weaker secondary alignments are ignored.

Alternative placements don't count as coverage, but they are used to type gaps and to fill the
scaffold; see [Repeats](../behaviour/repeats.md).

## Duplicate contigs

A reconstruction contig whose sequence (in either orientation) repeats an earlier contig's is
ignored, with a warning, so that a contig written twice isn't counted twice.
