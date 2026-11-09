# TODO

Open work for plasmid-surgery. Items marked *(PLSDB set)* come from exploring 25 public isolates
reassembled from their reads (Nov 2026; `exploration/plsdb-reassembly/`): with short reads, binning, not
assembly, accounts for most of what a reconstruction lacks. A formal validation comes later, once the
tool settles.

## Exploration (PLSDB set)

- [ ] Second pass with near-relative references: the closest PLSDB plasmids that aren't the isolate's
  own (closer to real use, where the database lacks the exact plasmid).
- [ ] Check read-truth concordance before using an isolate (one isolate's reads matched only 10-54% of
  its published plasmids in every assembler).
- [ ] Check what the bins matching no truth plasmid are (consistent across assemblers, so likely real
  replicons missing from the isolate's PLSDB set).
- [ ] Full mob_recon vs. TaDReP comparison against a full PLSDB/RefSeq database (cluster).
- [ ] Validate reference search on the cluster against a full PLSDB/RefSeq sketch.

## Assessment

- [ ] *(PLSDB set)* **Gaps at a circular contig's own junction.** Complete single-contig hybrid
  plasmids almost all show one small "unique gap" (2-947 bp, mostly 7-31) exactly where the contig's
  two ends meet (imprecise circularisation). Recognise these (adjacent ends of the same contig,
  `circular=Y`/`circular=true` in the header) and report them as a junction indel, not a plain gap.
- [ ] *(PLSDB set)* **Bins made of several references.** mob_recon put two complete truth plasmids in
  one bin (2 isolates, every assembler); the ranking lists both, but nothing says "this bin is A + B".
  Report a bin's composition: which contigs go to which reference, and the bin fraction explained by
  the references together.
- [ ] *(PLSDB set)* **No bin contains the reference.** Whole plasmids (a 46 kb NDM-5 plasmid in every
  assembler, 45 kb NDM-5, 2-4 kb plasmids) were binned as chromosome, so there's no reconstruction to
  point the tool at. Document passing `chromosome.fasta` (or the assembly) as the reconstruction, or
  add an assembly-only mode.
- [ ] **Structural consistency.** Walk placements around the circle and flag order breaks, strand
  flips and split alignments mapping to distant loci.
- [ ] Intersect gaps with the ref's IS/transposon features (from `--ref-features`, or GenBank): this
  likely explains most short-read gaps. Annotate gaps.bed / a notes column.
- [ ] Large insertions within alignments count as "aligned" contig bases; also write their sequence to
  the unaligned FASTA (or an insertions FASTA) so novel sequence is all in one place.
- [ ] Per-contig SNP and indel counts.
- [ ] Use contig read depth to estimate repeat copy number and check it against `num_placements`.
- [ ] A contig spanning *both* the origin and L/2 (longer than half the plasmid) is split either way;
  a third rotation, or merging split alignments adjacent on both query and ref, would cover it.
- [ ] Overlap resolution is "first by ref_start wins"; consider preferring higher MAPQ / identity /
  longer alignments.
- [ ] The primary copy of a repeat is minimap2's arbitrary choice (by a hash of the contig name):
  uppercase in the scaffold while the other copies are lowercase. Consider lowercasing all copies.

## Inputs and outputs

- [ ] GenBank input for refs (NCBI's usual download), if GFF3 isn't convenient.
- [ ] Write PAF, or SAM with `-Y`, so supplementary alignments keep correct query coordinates
  (currently hard-clipped).
- [ ] Scaffold output in standard formats: AGP + FASTA, gaps sized by estimate (or fixed 100 N when
  unknown), never filled from reference; unplaced contigs to a separate file.
- [ ] Optionally add `missed` contigs to the scaffold (marked), for a "best possible" reconstruction.
- [ ] Batch mode: sample sheet (or mob_recon output dirs + ref DB) -> one table across the
  back-catalogue. `plasmid-surgery-cohort` combines finished runs, but runs are started one at a time.
- [ ] `mafft` is in `environment.yml` but unused. Remove, or use it.

## Reports

- [ ] Cohort view across several references at once (e.g. one section per reference).
- [ ] Cluster samples in the cohort view by their coverage/gene profiles, rather than sorting.
- [ ] PNG export of maps (for tools that can't take SVG).
- [ ] More behaviour-report scenarios as behaviour is added (e.g. supplementary alignments from
  rearrangements, QC filter edge cases, unaligned contig ends).

## Later

- [ ] Validation: a set disjoint from every exploratory isolate, with questions and criteria fixed
  before looking at it. Optionally simulated reads (ART / InSilicoSeq) for exact truth.
- [ ] Assembly graph (GFA) input: check whether contigs on either side of a gap are connected by a
  graph path of about the right length (gap bridgeable vs. truly missing).
- [ ] `mash screen` as an alternative reference screen, if Jaccard-derived estimates prove too noisy
  for small reconstructions.
- [ ] Parse alignments into a typed record (dataclass) rather than dicts.
- [ ] Evaluate wrapping RagTag `scaffold` for ordering/orientation, focusing this tool on assessment.
- [ ] Prior art to compare against / borrow from: QUAST `-r` (genome fraction, misassemblies, Icarus),
  MUMmer `dnadiff`, RagTag `scaffold`.
