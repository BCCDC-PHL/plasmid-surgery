# Whole-assembly mode

A reconstruction can lack part of a plasmid for two different reasons: the sequence wasn't assembled,
or it was assembled but the contigs were put in another bin (another plasmid, or the chromosome).
Given the assembly the reconstruction was made from (`--assembly`), plasmid-surgery aligns every
assembly contig to the reference, and separates the two.

## Classification

Every assembly contig that aligns to the reference is listed in `contig_sources.tsv`, as one of:

| Classification | Meaning |
|---|---|
| `in_reconstruction` | Part of the reconstruction being assessed |
| `shared` | Less than `--min-missed-fraction` (default 0.8) of the contig aligns to the reference. Likely a contig from the chromosome or another plasmid sharing some sequence with this one (e.g. a transposon), rather than a binning error |
| `repeat` | Most of the contig aligns to the reference, but equally well at more than one place (e.g. an IS element), so it may belong on this plasmid, elsewhere, or both |
| `missed` | Most of the contig aligns to the reference, at a single location: likely a plasmid contig placed in another bin |

`missed` contigs are written to `missed_contigs.fa`.

## Bins

`--bins` takes any number of FASTA files: typically all of the reconstruction tool's outputs (its
plasmid FASTAs and chromosome FASTA). Each contig is labelled with the bin (file name without
extension) that contains it, or `unbinned`. Contigs are matched between the reconstruction, the
assembly and the bins by sequence, in either orientation, not by name, so tools that rename or
reorient contigs work too.

## Read depth

Where contig headers carry read depth (`depth=2.3x`, `_cov_12.5`, `cov=12.5`), it's reported, and
relative to the reconstruction's contigs (length-weighted mean). A `missed` contig at about the same
relative depth as the reconstruction is consistent with belonging to the same plasmid; one at
chromosomal depth less so. Depth is noisy for short contigs.

## Metrics

| Field | Meaning |
|---|---|
| `assembly_ref_fraction_covered` | Reference coverage by the whole assembly |
| `assembly_unique_gap_length` | Reference sequence no assembly contig covers: missing from the assembly itself |
| `missed_ref_bases` | Reference sequence covered by `missed` contigs and not by the reconstruction: lost by binning |

Interpreting a contig is a judgement, and the classification comes with the evidence for it:
`contig_fraction_on_ref`, `num_placements` and relative depth. See
[Whole-assembly mode](../behaviour/assembly.md) for an example, including a borderline contig, and
[The choice of reference](../behaviour/reference-choice.md) for what no single reference can show.
