# Coverage, gaps and events

## Coverage

A reference base is *covered* if a placement alignment spans it. The exception is a deletion within
an alignment of at least `--min-indel-length` (default 50) bp: the contig doesn't have that sequence,
so it isn't covered. Shorter deletions count as covered, so that small indels don't break coverage
into pieces. Bases covered by more than one alignment (overlapping contigs, or duplicated sequence)
are counted in `ref_bases_covered_multiple`.

## Gaps

Uncovered regions are the gaps, written to `gaps.bed`. Each is typed:

* **`repeat`**: an alternative placement covers it. The sequence is in the reconstruction (in a
  collapsed repeat contig) but can't be placed there with certainty.
* **`unique`**: nothing covers it. The sequence is likely missing from the reconstruction.

A gap partly covered by an alternative placement is split into a repeat part and a unique part. A
gap spanning the origin is a single gap, written as two BED records with the same name.
`ref_fraction_covered_including_repeats` counts repeat gaps as covered.

## Events

Insertions and deletions of at least `--min-indel-length` bp within an alignment are written to
`events.tsv`, with their reference and contig positions. A large deletion is therefore both an event
and a unique gap. An insertion is sequence the reference doesn't have, inside an aligned contig: its
bases count towards `contig_fraction_aligned`, and are reported as an event rather than as unaligned
sequence.

## Unaligned sequence

Contig regions of at least `--min-unaligned-length` (default 100) bp outside any alignment (whole
contigs that don't align, or unaligned contig ends) are written to `unaligned.fa`.

See [Insertions, deletions and identity](../behaviour/indels.md) for an example.
