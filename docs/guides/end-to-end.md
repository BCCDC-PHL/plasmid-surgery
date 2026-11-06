# Assessing a reconstruction, end to end

A typical workflow, from an assembly with a plasmid reconstruction to an interpretation.

## 1. Inputs

* **The assembly**, as contigs (FASTA). Headers with read depth (e.g. `depth=2.3x`) are used as
  supporting evidence.
* **The reconstruction** of the plasmid of interest: a FASTA of its contigs. Use the contigs
  themselves, not a joined pseudomolecule with placeholder `N` gaps.
* **The reconstruction tool's other bins**, if it has them: its other plasmid FASTAs and its
  chromosome FASTA (e.g. MOB-recon's `plasmid_*.fasta` and `chromosome.fasta`).
* **One or more candidate reference plasmids** (FASTA): complete, closed plasmids believed to be
  related. Everything is measured against the reference, so see step 4.

## 2. Annotate (optional)

Annotate the assembly and the reference with the same tools, and convert to hAMRonization TSV. For
example, with abricate and two databases:

```bash
for db in ncbi plasmidfinder; do
  for f in sample1.fasta CP012345.fa; do
    abricate --db $db $f > ${f%.*}.$db.tsv
    hamronize abricate ${f%.*}.$db.tsv --analysis_software_version 1.4.0 \
      --reference_database_version unknown --format tsv --output ${f%.*}.$db.hamronized.tsv
  done
done
```

GFF3 (e.g. from a genome annotation pipeline) works too, and covers features that resistance and
replicon databases don't, such as transposases.

## 3. Run

```bash
plasmid-surgery \
  --plasmid-reconstruction bins/plasmid_1.fasta \
  --ref-plasmid CP012345.fa \
  --assembly sample1.fasta \
  --bins bins/*.fasta \
  --ref-features CP012345.*.hamronized.tsv \
  --features sample1.*.hamronized.tsv \
  -o results/sample1_vs_CP012345 -p sample1
```

## 4. Interpret

Work from the summary to the detail:

1. **Is this reference plausible?** `contig_fraction_aligned` near 1 and few `num_unaligned_contigs`
   mean the reconstruction is accounted for by the reference. A good result means "consistent with
   this reference", not "complete": a fragment of a larger plasmid can look complete against a
   smaller reference that shares it (see [The choice of reference](../behaviour/reference-choice.md)).
   Where the reference is uncertain, assess against more than one candidate.
2. **How complete is the reconstruction?** `ref_fraction_covered`, and
   `ref_fraction_covered_including_repeats`. The difference is collapsed repeats.
3. **What's missing?** `unique_gap_length`, and `gaps.bed` for where. Compare with
   `assembly_unique_gap_length`: if the assembly covers what the reconstruction doesn't, the loss is
   binning, and `contig_sources.tsv` shows which contigs, in which bins (`missed`).
4. **What's different?** `events.tsv` (insertions, deletions) and `unaligned.fa`.
5. **Genes.** In `ref_features.tsv`, reference genes not `present` in the reconstruction; if
   `found_in` another bin, binning; if `partial` or `absent` in the assembly too, the assembly
   doesn't have them complete. In `features.tsv`, sample genes on `missed` contigs (binned away from
   this plasmid), and `not_on_ref` genes (sequence the reference doesn't have).

### A gene on a contig binned as chromosome

If the contig is `missed` against a plausible reference (most of it aligns, at one place), and its
`depth_relative_to_reconstruction` is close to 1, it's most likely a plasmid contig in the wrong bin.
If it's `shared` (only part of it is on the reference), check whether the gene itself is on the
reference (`ref_status` in `features.tsv`): a gene in the contig's unaligned part, at chromosomal
depth, makes a chromosomal location more plausible. This is evidence, not proof.
