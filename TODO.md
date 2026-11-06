# TODO

Roadmap for plasmid-surgery, based on a review (Nov 2026) that ran the tool on a
synthetic plasmid with known ground truth: a 60.9 kb reference with 3 copies of a
1.3 kb IS-like repeat; a sample with 300 SNPs, a 2 kb deletion, a 4 bp deletion, a
3 kb novel insertion and a rotated origin; broken into 5 contigs at the repeats
(one contig spanning the reference origin, alternating orientations).

Findings in brief: per-contig placement (strand + reference position) was correct,
including the origin-spanning contig. The output FASTA construction was broken
whenever contigs had indels relative to the reference.

## Phase 1: Bug fixes and small issues

- [x] **Splicing bug.** Contigs were spliced into an N-string with
  `output[:ref_start] + contig + output[ref_end:]`. Because aligned query length
  differs from the ref span (indels), each splice shifted all downstream
  coordinates and later contigs were written at stale positions, overwriting
  earlier ones (in the test, contig_0 was overwritten after ~700 bp, others were
  displaced by 2 kb). Output is now built in a single left-to-right sweep over
  reference coordinates.
- [x] **No-op overlap check.** `ref_end <= ref_length + 1` was always true, so
  overlapping alignments overwrote each other. Overlaps (e.g. the k-1 bp overlap
  between adjacent SPAdes contigs) are now trimmed via the CIGAR; fully contained
  alignments are excluded from the output.
- [x] `__main__.py` had no `if __name__ == "__main__"` guard, so
  `python -m plasmid_surgery` did nothing.
- [x] `list[dict]` type hints broke Python 3.8; Python version now pinned.
- [x] `--output-bases-per-line` lacked `type=int` (crashed if passed).
- [x] `--output out.fa` (no directory part) crashed in `os.makedirs('')`.
- [x] `subprocess.run` lacked `check=True`, so minimap2 failures were silent.
- [x] Ref FASTA in tmpdir was only written if absent: a different ref with the
  same ID in the same tmpdir was silently reused.
- [x] `parse_fasta` used `lstrip(seq_id)`, which strips a character set, not a
  prefix (corrupted descriptions).
- [x] Percent identity ignored indels (`= / (= + X)`), overstating identity.
  Now gap-compressed: matches / (matches + mismatches + indel events), as in
  minimap2's `de` tag. (BLAST-style identity was tried first, but a single 3 kb
  insertion dropped an 18 kb contig to 83%, close to the 80% QC cutoff, so
  contigs with real structural indels would have been discarded. Indel content is
  already filtered separately by `gap_frac`.)
- [x] Possible divide-by-zero in `alignment_qc` for zero-length alignments.
- [x] `[circular=true]` in ref header had no effect on minimap2; removed.
- [x] `pyproject.toml` was empty; packaging metadata moved there from `setup.py`.
- [x] Summary sort used a magic `FAKE_ALIGNMENT_START`; replaced with a sort key.
- [x] Added `mapq` and `is_supplementary` to summary output.
- [x] README logo link was `[!Logo](...)` rather than `![Logo](...)`.
- [ ] `mafft` is in `environment.yml` but unused. Remove, or use it.
- [ ] Ref file must contain a single sequence; extra sequences are now warned
  about and ignored (see multi-reference support in Phase 5).

## Phase 2: Tests

- [x] pytest suite: unit tests for parsers, identity, CIGAR trimming, output
  construction.
- [x] Simulated-data integration test (needs `minimap2` on PATH): known ref, known
  edits, rotated origin, contigs broken at repeats, mixed orientations; assert
  every contig is recovered intact in the output.
- [ ] Larger validation: simulate 150 bp reads (ART / InSilicoSeq) from closed
  plasmids in our collection, run the real assembly -> mob_recon -> plasmid-surgery
  chain and measure against truth.
- [x] CI (GitHub Actions) running the test suite.

## Phase 3: Make it an assessment tool first

The goal is to assess how well a reconstruction matches a reference. A
ref-length, N-filled FASTA is a reference-guided consensus that hides exactly the
differences of interest. Make metrics the primary output and the scaffold an
optional by-product.

- [x] **Model the circle explicitly.** Each contig is aligned against the ref and
  a copy rotated by L/2; the alignment covering more of the contig in fewer
  segments wins, and coordinates are folded back to the original ref (an
  alignment with `ref_end > L` wraps the origin and is split there when building
  output). This fixed a case where minimap2 split an origin-spanning contig with
  7 bp aligned in neither piece. A doubled ref was rejected because it makes
  every alignment MAPQ 0.
  - [ ] Remaining edge case: a contig spanning *both* the origin and L/2 (i.e.
    longer than half the plasmid, which is rare for short-read contigs) is split
    either way. A third rotation, or merging split alignments that are adjacent
    on both query and ref, would cover this.
  - [ ] Two alignments per contig doubles minimap2 runtime; fine for now, revisit
    together with the "one minimap2 call per sample" item in Phase 4.
- [x] **Reference coverage.** Breadth covered 1x, covered >1x; BED of uncovered gaps
  with lengths.
- [x] **Reconstruction accounted for.** Per contig, fraction of contig length
  aligned. Unaligned contigs and long soft-clipped ends are candidate insertions,
  chromosomal contamination or a different plasmid. Report them, don't drop them.
- [x] **Deletions visible.** Large `D` operations inside an alignment (like the
  2 kb test deletion) are currently absorbed silently; report them as events.
  Done as `--assessment-output` (one row of metrics), `--gaps-output` (BED),
  `--events-output` (insertions/deletions >= `--min-indel-length`) and
  `--unaligned-output` (FASTA of contig regions not aligned to the ref).
  - [ ] Large insertions within alignments count as "aligned" contig bases;
    consider also writing their sequence to the unaligned FASTA (or a separate
    insertions FASTA) so novel sequence is all in one place.
  - [ ] Distinguish gaps at repeat copies (likely just collapsed repeats) from
    gaps in unique sequence (likely genuinely missing), once repeat awareness
    and/or ref annotation are in.
- [ ] **Structural consistency.** Walk placements around the circle and flag order
  breaks, strand flips and split alignments mapping to distant loci.
- [ ] **Repeat awareness.** Use MAPQ and multiple alignments (`-N 5` instead of
  `--secondary=no`) to flag contigs that place ambiguously. IS26 / Tn4401 etc. are
  ubiquitous in carbapenemase plasmids and will cause this constantly.
  Currently a MAPQ 0 repeat contig is placed at an arbitrary copy.
- [ ] If the ref is annotated (GenBank), intersect gaps with IS/transposon
  features: this likely explains most short-read gaps.
- [ ] Overlap resolution is currently "first by ref_start wins". Consider
  preferring higher MAPQ / higher identity / longer alignments.
- [ ] Per-contig SNP and indel counts.

## Phase 4: Architecture and outputs

- [ ] One minimap2 call per sample (all contigs) instead of one per contig; write
  PAF, or SAM with `-Y` so supplementary alignments keep correct query
  coordinates (currently hard-clipped).
- [ ] Avoid per-contig temp files named after contig IDs (unsafe characters).
- [x] Use `tempfile` and clean up by default (`--tmpdir` keeps intermediates).
- [x] `-o/--outdir` + optional `-p/--prefix` instead of one flag per output file.
- [ ] Parse alignments into a typed record (dataclass) rather than dicts.
- [ ] Scaffold output in standard formats: AGP + FASTA, gaps sized by estimate (or
  fixed 100 N when unknown), never filled from reference; unplaced contigs written
  to a separate file.
- [ ] Evaluate wrapping RagTag `scaffold` for ordering/orientation and focusing
  this tool on assessment and plasmid-specific logic.

## Phase 5: Scale and reporting

- [ ] Multiple candidate references per reconstruction; rank them.
- [ ] Batch mode: sample sheet (or mob_recon output dirs + ref DB) -> one table
  across the back-catalogue.
- [ ] Visualization: linear/circular plot of contigs on ref, coloured by strand,
  gaps and MAPQ 0 placements highlighted (self-contained HTML report, or
  pyCirclize).

## Prior art to compare against / borrow from

- QUAST with `-r` (genome fraction, misassemblies, unaligned contigs, Icarus viewer)
- MUMmer `dnadiff`
- RagTag `scaffold`
- Graph-aware checks against the assembly graph (GFA), to tell bridgeable gaps
  from truly missing sequence
