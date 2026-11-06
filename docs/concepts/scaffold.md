# The scaffold

`scaffold.fa` is the reconstruction laid out in reference order: one sequence, starting at the
reference's origin.

It's built by walking along the reference in order of alignment start:

1. Each placement contributes its contig's aligned sequence (reverse-complemented if it aligned to the
   reverse strand).
2. An alignment overlapping the region already covered has its overlapping start trimmed, using its
   alignment to find the matching contig position, so each overlap appears once. The first contig by
   position wins where overlapping contigs disagree.
3. An alignment entirely within the region already covered is left out (`included_in_output` =
   False).
4. An alignment spanning the origin is split there; the part past the origin goes at the start.
5. Uncovered reference regions are filled from alternative placements (repeat gaps), in
   **lowercase** to mark them as inferred, and otherwise with one `N` per uncovered reference base.

Lowercase follows the common soft-masking convention. Some tools treat lowercase as masked, and others
uppercase their input, so inferred regions may be ignored, or indistinguishable, downstream.

The scaffold is built from the contigs' sequence, never from the reference's: where the sample
differs from the reference, the scaffold has the sample's sequence. Gaps filled with `N` are sized
by the reference, so they're estimates.

See [Overlapping and contained contigs](../behaviour/overlaps.md) and [Repeats](../behaviour/repeats.md).
