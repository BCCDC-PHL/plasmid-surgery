# Many samples

`assessment.tsv` is one row with a header, and its columns are the same whichever options were used
(options not used give `NA`), so the rows for many samples can be concatenated into one table.

```bash
while read sample ref; do
  plasmid-surgery \
    --plasmid-reconstruction reconstructions/$sample.fasta \
    --ref-plasmid refs/$ref.fa \
    --assembly assemblies/$sample.fasta \
    --sample-id $sample \
    -o results/$sample -p ${sample}_vs_$ref
done < samples.tsv

# One table: the first file's header, then every file's row
first=$(ls results/*/*_assessment.tsv | head -1)
{ head -1 "$first"; for f in results/*/*_assessment.tsv; do tail -n +2 "$f"; done; } > assessment_all.tsv
```

Use `--sample-id` so the `sample_id` column identifies each row, and `-p` so several comparisons
(e.g. one sample against several candidate references) can share an output directory without
overwriting each other. Temporary files are per run, so runs can go in parallel.
