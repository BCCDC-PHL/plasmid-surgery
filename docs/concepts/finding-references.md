# Finding references

`plasmid-surgery-refs` finds suitable reference plasmids for a reconstruction in a plasmid database:
a mash sketch of the database, and the FASTA it was made from (with sketch IDs matching the FASTA's
sequence IDs).

```bash
plasmid-surgery-refs \
  --plasmid-reconstruction bins/plasmid_1.fasta \
  --db-sketch plasmids.msh --db-fasta plasmids.fasta \
  --assembly sample1.fasta \
  -o refs/sample1 -p sample1 -t 8
```

It writes `candidates.tsv` (every candidate, ranked) and `references/`, one FASTA per selected
reference (`--top`, default 3), named `<rank>_<id>.fa`, ready to use as `--ref-plasmid`.

## How it works

1. **Screen** with mash. The reconstruction is sketched as a whole, with the database's k-mer and
   sketch sizes, and compared with every database sketch in one `mash dist`. mash's distance comes
   from an estimate of the Jaccard index of the two k-mer sets. That's symmetric, and penalises a
   difference in size even if one set contains the other, so on its own it would miss a large
   plasmid containing the reconstruction. With the two sizes (the database's from the sketch), the
   same estimate gives containment in each direction:

    ```text
    |A ∩ B| = J / (1 + J) × (|A| + |B|)
    ```

    The candidates are the top `-n` (default 20) by estimated containment of the reconstruction in
    the reference, plus the top `-n` by mash distance.

2. **Extract** the candidates from the FASTA: by random access with a samtools faidx index
   (`<fasta>.fai`, or `--db-fasta-index`), or else in one pass through the FASTA. Nothing is written
   next to the database.

3. **Confirm** with minimap2: the reconstruction's contigs (and with `--assembly`, the assembly's) are
   aligned to each candidate, giving the exact fraction of the reconstruction aligned, the fraction of
   the candidate covered (with deletions of 50 bp or more uncovered, as in the main assessment), and
   gap-compressed identity.

4. **Rank** by

    ```text
    score = (fraction of the reconstruction aligned) × (fraction of the reference covered)
    ```

    where the coverage is the reconstruction's, or with `--assembly`, the whole assembly's. Ties are
    broken by identity, then mash distance.

## Why this score, and its trade-off

Both factors are near 1 only when the reference is the same plasmid as the reconstruction. A
reference that explains the reconstruction but has much more sequence besides is penalised by the
second factor; one that's mostly covered but explains only part of the reconstruction, by the first.

But if the reconstruction is only a fragment of the real plasmid, the real plasmid is only partly
covered by it, and a smaller relative sharing the fragment can score higher. That's why `--assembly`
changes the coverage factor to the whole assembly's: the rest of the real plasmid is usually somewhere
in the assembly, even if not in the reconstruction. See [Finding references](../behaviour/find-refs.md)
in the Behaviour report.

`candidates.tsv` reports each figure separately (including mash's own estimates), so a different
balance can be applied. If the best candidate accounts for less than 80% of the reconstruction, a
warning says the database may not contain a closely related plasmid.

## Limitations

* Confirmation uses one linear alignment per contig, so a few bases at a candidate's origin can be
  missed. The full assessment (`plasmid-surgery`) handles the origin exactly.
* Containment estimates from mash are noisy for small reconstructions (few k-mers in the sketch); the
  confirmation step gives the exact figures for the candidates.
* The search can only find what's in the database.
