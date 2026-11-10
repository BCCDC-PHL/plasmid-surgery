# TODO

Open work for plasmid-surgery. Items marked *(PLSDB set)* come from exploring 25 public isolates
reassembled from their reads (Nov 2026; `exploration/plsdb-reassembly/`): with short reads, binning, not
assembly, accounts for most of what a reconstruction lacks. A formal validation comes later, once the
tool settles.

## Exploration (PLSDB set)

Next: both passes (truth, near-relative) again with reads (reads input done, below).

- [x] **Near-relative references** (2026-11-10, `near/FINDINGS.md`, not committed): stand-ins up to
  ~0.003 behave like the truth; drift grows from ~0.005; see the two Assessment items below. Repeated
  with a carbapenemase-plasmid database (`near-neighbor-with-carbapenemase/FINDINGS.md`): as close for
  carbapenemase plasmids, and more often carrying the gene; much worse for a sample's other plasmids.
- [ ] **Rerun both passes (truth, near-relative) with reads**: `scripts/rerun_with_reads.sh` (download ->
  plasmid-surgery for every reference set -> delete; ~7 min for SAMN40632425, the first isolate).
  That isolate: the 2.0 kb plasmid dragonflye lacks is in the reads (~6,700x), so the assembler dropped
  it; the short-read assemblies have it but no plasmid bin does.
- [ ] Check read-truth concordance before using an isolate (one isolate's reads matched only 10-54% of
  its published plasmids in every assembler). Covered by reads breadth on the reference, once added.
- [ ] Junction gaps (see Assessment): assemble the hybrid isolates' long reads with the production
  long-read pipeline (incl. dnaapler) and check whether they persist. dragonflye already runs dnaapler,
  and the gap sits at the rotated contig's ends (polishing a linearised circle?). With reads
  (SAMN40632425), the 7-13 bp junction gaps are all in the reads: an assembly error, not biology.
- [ ] Negatives for gene placement: assemblies where a carbapenemase is genuinely chromosomal (closed
  genome), from in-house hybrids or NCBI complete genomes. The PLSDB set only has positives.
- [ ] Check what the bins matching no truth plasmid are (consistent across assemblers, so likely real
  replicons missing from the isolate's PLSDB set).
- [ ] Full mob_recon vs. TaDReP comparison against a full PLSDB/RefSeq database (cluster).
- [ ] Validate reference search on the cluster against a full PLSDB/RefSeq sketch.

## Reads as input

Optional short and/or long reads (minimap2 `-ax sr` / `map-ont`; samtools). Keep v1 to mapping, depth
and spanning reads; leave variant calling to dedicated tools.

- [x] **Reads vs. reference:** breadth and depth over the reference (`--short-reads`, `--long-reads`;
  `reads_depth.bedgraph`, `reads_ref_fraction_covered`), so losses split three ways: not in the reads /
  not in the assembly / not in the bin. Also `reads_fraction` per ref feature.
- [x] **Gap adjudication:** `gap_support.tsv` calls each reconstruction and assembly gap `in_reads` /
  `not_in_reads` / `partly_in_reads`; `unique_gap_in_reads_length`, `assembly_unique_gap_in_reads_length`.
- [ ] **Circular junctions:** reads spanning a circular contig's join confirm circularity and the
  exact join (needs reads aligned to the contigs, not only the reference).
- [x] **Small plasmids missing from the assembly:** covered by reads breadth on the reference (with
  `--assembly`: `assembly_unique_gap_in_reads_length`).
- [x] Contig read depth from reads aligned to the assembly (`contig_depth.tsv`, `reads_depth` in
  `contig_sources.tsv`; consistent across assemblers). Next: use it for gene placement's copy-number
  evidence (relative to the chromosome bin).
- [x] Read depth in the HTML report: a depth track under each map, read support in gap tooltips, read
  columns in the tables, contig read depth in the contig table.
- [ ] Later: read support (allele depth) at the SNPs/indels the assessment already reports.

## Gene placement

- [ ] **Is a chromosome-binned resistance gene really chromosomal?** Starting from the contig carrying
  a gene (from `--features`), gather evidence for each placement: long reads spanning the gene and its
  flanking mobile element into unambiguous chromosome or plasmid sequence (strongest; reads input);
  depth relative to chromosome and plasmid bins (copy number); plasmid-DB and chromosome matches of the
  contig plus flanks (`plasmid-surgery-refs`; a plasmid hit alone isn't enough, these elements also
  insert into chromosomes); assembly-graph links (GFA, below). Check on the PLSDB positives (6-7
  short-read carbapenemases in the chromosome bin per assembler) and on true negatives.

## Assessment

- [ ] *(PLSDB set)* **Gaps at a circular contig's own junction.** Complete single-contig hybrid
  plasmids almost all show one small "unique gap" (2-947 bp, mostly 7-31) exactly where the contig's
  two ends meet (imprecise circularisation). Recognise these (adjacent ends of the same contig,
  `circular=Y`/`circular=true` in the header) and report them as a junction indel, not a plain gap.
  On hold until checked on the production long-read pipeline (Exploration); reads would settle it.
- [ ] *(PLSDB set)* **Bins made of several references.** mob_recon put two complete truth plasmids in
  one bin (2 isolates, every assembler); the ranking lists both, but nothing says "this bin is A + B".
  Report a bin's composition: which contigs go to which reference, and the bin fraction explained by
  the references together.
- [ ] *(PLSDB set)* **No bin contains the reference.** Whole plasmids (a 46 kb NDM-5 plasmid in every
  assembler, 45 kb NDM-5, 2-4 kb plasmids) were binned as chromosome, so there's no reconstruction to
  point the tool at. Document passing `chromosome.fasta` (or the assembly) as the reconstruction, or
  add an assembly-only mode.
- [ ] *(PLSDB near-relative pass)* **Report the reference's distance** (mash distance / identity to the
  reconstruction) and flag far references: from ~0.005 the reference's own differences dominate
  coverage and gaps (-6% at 0.005, -19% at 0.02), though bin choice and the bin-vs-assembly split
  stay robust.
- [ ] *(PLSDB near-relative pass)* **Genes on the reconstruction but not on the reference.** The nearest
  relative lacked the carbapenemase for 5 of 26 plasmids (10 of 18 at ~0.01). Report `--features`
  genes that have no counterpart on the reference (same gene family), so gene status doesn't depend
  on the reference carrying it.
- [ ] **Structural consistency.** Walk placements around the circle and flag order breaks, strand
  flips and split alignments mapping to distant loci.
- [ ] Intersect gaps with the ref's IS/transposon features (from `--ref-features`, or GenBank): this
  likely explains most short-read gaps. Annotate gaps.bed / a notes column.
- [ ] Large insertions within alignments count as "aligned" contig bases; also write their sequence to
  the unaligned FASTA (or an insertions FASTA) so novel sequence is all in one place.
- [ ] Per-contig SNP and indel counts.
- [ ] Use contig read depth to estimate repeat copy number and check it against `num_placements`
  (also feeds gene placement).
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
  graph path of about the right length (gap bridgeable vs. truly missing), and which contigs a gene
  contig links to (gene placement, short reads).
- [ ] `mash screen` as an alternative reference screen, if Jaccard-derived estimates prove too noisy
  for small reconstructions.
- [ ] Parse alignments into a typed record (dataclass) rather than dicts.
- [ ] Evaluate wrapping RagTag `scaffold` for ordering/orientation, focusing this tool on assessment.
- [ ] Prior art to compare against / borrow from: QUAST `-r` (genome fraction, misassemblies, Icarus),
  MUMmer `dnadiff`, RagTag `scaffold`.
