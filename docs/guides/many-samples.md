# Many samples

Each run's top-level `assessment.tsv` has a header and one row per reference, and its columns are the
same whichever options were used (options not used give `NA`), so the rows for many samples can be
concatenated into one table.

```bash
while read sample ref; do
  plasmid-surgery \
    --plasmid-reconstruction reconstructions/$sample.fasta \
    --ref-plasmid refs/$ref.fa \
    --assembly assemblies/$sample.fasta \
    --sample-id $sample \
    -o results/$sample
done < samples.tsv

# One table: the first file's header, then every file's rows (each run's top-level assessment.tsv)
first=$(ls results/*/assessment.tsv | head -1)
{ head -1 "$first"; for f in results/*/assessment.tsv; do tail -n +2 "$f"; done; } > assessment_all.tsv
```

Use `--sample-id` so the `sample_id` column identifies each row. Several references for one sample
go in one run (`--ref-plasmid` with several files, or `--db-sketch`), each in its own subdirectory.
Temporary files are per run, so runs can go in parallel.
