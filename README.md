# Plasmid Surgery

[![Tests](https://github.com/BCCDC-PHL/plasmid-surgery/actions/workflows/tests.yml/badge.svg)](https://github.com/BCCDC-PHL/plasmid-surgery/actions/workflows/tests.yml)
[![Docs](https://github.com/BCCDC-PHL/plasmid-surgery/actions/workflows/docs.yml/badge.svg)](https://bccdc-phl.github.io/plasmid-surgery/)

![Logo](docs/images/plasmid-surgery.png)

A plasmid reconstructed from a short-read assembly is usually a set of contigs, without their order,
orientation, or any indication of what's missing. **plasmid-surgery** assesses such a reconstruction
against a reference plasmid: where each contig belongs, which parts of the reference are covered,
which gaps are collapsed repeats and which are missing sequence, what the sample has that the
reference doesn't, and, given the full assembly, where the rest of the plasmid ended up. Gene
annotations of the reference and of the sample can be added, to report the same at gene level.

**Documentation: <https://bccdc-phl.github.io/plasmid-surgery/>**, including a
[Behaviour report](https://bccdc-phl.github.io/plasmid-surgery/behaviour/) that demonstrates each
design decision on a small simulated plasmid.

## Install

In a dedicated conda environment (needs minimap2 and pysam; `plasmid-surgery-refs` also needs mash):

```bash
conda env create -f environment.yml
conda activate plasmid-surgery
pip install -e .
```

## Usage

Three ways to work:

| To | Run |
|---|---|
| assess against references you have | `plasmid-surgery --ref-plasmid a.fa [b.fa …]` |
| find references in a plasmid database, then assess against the best | `plasmid-surgery --db-sketch db.msh --db-fasta db.fa` |
| only find references | `plasmid-surgery-refs --db-sketch db.msh --db-fasta db.fa` |

```bash
plasmid-surgery \
  --plasmid-reconstruction bins/plasmid_1.fasta \
  --ref-plasmid CP012345.fa \
  -o results/sample1 -p sample1
```

Each reference's outputs go to `<outdir>/<ref_id>/`, and `<outdir>/assessment.tsv` has one row per
reference, ranked. Add `--report` for a self-contained `report.html` with every reference's results on
linear maps (or run `plasmid-surgery-plot <outdir>` afterwards). `plasmid-surgery-cohort results/*/`
compares many samples against one shared reference in one page (coverage strips; genes by sample); its `--list-refs` shows which references rank well across the runs.

Optionally, with the full assembly and the reconstruction tool's other bins (`--assembly`,
`--bins`), and with gene annotations of the reference and the sample (`--ref-features`,
`--features`; hAMRonization TSV, GFF3 or BED):

```bash
plasmid-surgery \
  --plasmid-reconstruction bins/plasmid_1.fasta \
  --ref-plasmid CP012345.fa \
  --assembly sample1.fasta --bins bins/*.fasta \
  --ref-features CP012345.hamronized.tsv --features sample1.hamronized.tsv \
  -o results/sample1 -p sample1
```

### Finding references

`plasmid-surgery-refs` finds candidate references in a plasmid database (a mash sketch and the FASTA
it was made from): it screens every database plasmid with mash (distance, and estimated containment
in both directions), confirms the best candidates with minimap2, and ranks them. The best are written
out as FASTAs, ready for `--ref-plasmid`:

```bash
plasmid-surgery-refs \
  --plasmid-reconstruction bins/plasmid_1.fasta \
  --db-sketch plasmids.msh --db-fasta plasmids.fasta \
  --assembly sample1.fasta \
  -o refs/sample1 -t 8
```

## Outputs

| File | Contents |
|---|---|
| `assessment.tsv` | Top level: one row of metrics per reference, ranked: coverage, gaps (repeat vs. unique), unaligned sequence, events, and (with the options above) assembly-level and gene-level metrics. Rows from many samples concatenate into one table. |
| `candidates.tsv` | Top level, with `--db-sketch`: every candidate reference considered |
| *In `<ref_id>/`:* | |
| `summary.tsv` | One row per contig alignment: position, orientation, identity, placements |
| `gaps.bed` | Uncovered reference regions, typed `repeat` or `unique` |
| `events.tsv` | Large insertions and deletions within aligned contigs |
| `unaligned.fa` | Contig sequence that doesn't align to the reference |
| `scaffold.fa` | The reconstruction in reference order, with inferred repeat copies in lowercase |
| `contig_sources.tsv`, `missed_contigs.fa` | With `--assembly`: every assembly contig aligning to the reference, classified, with its bin |
| `ref_features.tsv`, `features.tsv` | With annotations: reference genes' completeness; sample genes' locations |
| `report.html`, `<ref_id>/map.svg` | With `--report`: the HTML report, and each reference's map for slides |

See [Outputs](https://bccdc-phl.github.io/plasmid-surgery/reference/outputs/) for every field.

## Development

```bash
conda activate plasmid-surgery
pip install -e '.[test,docs]'
pytest                # includes the behaviour report's scenarios (needs minimap2)
mkdocs serve          # live preview of the documentation
```

The behaviour report (`docs/behaviour/`) is generated at build time from `tests/scenarios/`, and
every claim in it is checked by `tests/test_scenarios.py`.
 
