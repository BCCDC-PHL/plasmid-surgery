# PLSDB reassembly set (exploratory)

**Phase: exploratory.** Gathering evidence about how plasmid-surgery behaves on typical public
assemblies; findings may change the tool. This set is "training data": a later test phase should use
isolates disjoint from these (`select_isolates.py --exclude`).

Only the selection, scripts and small provenance tables are committed (see `.gitignore`); the reads,
assemblies, mob_recon and plasmid-surgery outputs are regenerated with the scripts below.

25 isolates whose closed plasmids are in PLSDB (2024_05_31_v2), reassembled from their public reads with
contemporary, popular assemblers (not production versions). The published plasmids are the truth;
they were usually finished from these same reads.

## How it was made

1. `fetch_data.py select` (mash-cluster `validation/public/`): PLSDB plasmids that are complete,
   circular and carry a carbapenemase (AMRFinderPlus bla* with CARBAPENEM): 5,861 from 4,542 BioSamples
   -> `plasmids.tsv`.
2. `fetch_data.py runs`: every read run for those BioSamples (ENA), nanopore chemistry (NCBI)
   -> `runs.tsv`. 862 BioSamples have Illumina + nanopore; 384 Illumina only.
3. `scripts/select_isolates.py --hybrid 15 --short 10 --exclude <round-1 isolates>` -> `isolates.tsv`.
   Needs paired-end Illumina WGS, 0.3 Gb+ and up to 1.5 GB (run nearest 1 Gb chosen); nanopore up to 2.5 GB
   (largest chosen). Spread over carbapenemase families (quotas ~ sqrt(family size), families under 2%
   pooled as 'rare'), random within a family, at most 3 isolates per species. Seed 1.
4. `scripts/process_isolates.sh`: one isolate at a time (storage is limited), reads deleted after:

| Step | Tool | Settings |
|---|---|---|
| reads | fetch_data.py download | Illumina from NCBI's full-quality ("normalized") SRA copy, as ENA's FASTQ is often SRA-lite (placeholder qualities); nanopore from ENA, md5-checked |
| truth | NCBI efetch | every PLSDB plasmid of the BioSample (not only the carbapenemase ones) -> `truth/` |
| short (1) | shovill 1.1.0, SPAdes 4.3.0 | defaults (`--isolate` mode), `--cpus 16 --ram 32` |
| short (2) | fastp 1.4.0, then Unicycler 0.5.1 (SPAdes 4.3.0) | fastp `--cut_tail --trim_poly_g --overrepresentation_analysis --detect_adapter_for_pe`; Unicycler defaults |
| hybrid | dragonflye 1.2.1 (Flye 2.9.6) | `--racon 1 --polypolish 1`, `--nanohq` for R10.4.1 runs, `--gsize` from shovill's estimate (see below) |
| bins | mob_recon 3.1.9 | defaults, on each assembly |

Environments: `shovill-1.1.0`, `unicycler-0.5.1`, `dragonflye-1.2.1`, `mob-suite`, `sra-tools`.

## Running and resuming

```bash
# Once: the environments (mob-suite also needs its database: mob_init)
for f in scripts/envs/*.yml; do conda env create -f "$f"; done

# Selection (steps 1-3; fetch_data.py caches API responses in .cache/)
python scripts/fetch_data.py --cache .cache select --plsdb $PLSDB -o plasmids.tsv
python scripts/fetch_data.py --cache .cache runs -i plasmids.tsv -o runs.tsv --threads 8
python scripts/select_isolates.py --plsdb $PLSDB --plasmids plasmids.tsv --runs runs.tsv \
  --hybrid 15 --short 10 --exclude SAMN26420017 SAMN36514289 SAMN38646027 -o isolates.tsv

# Assembly (step 4), from this directory; threads, RAM in GB
PLSDB=$PLSDB nohup scripts/process_isolates.sh isolates.tsv 16 32 > logs/process.out 2>&1 &
tail -f logs/progress.log

# plasmid-surgery on every mob_recon bin vs. the isolate's truth plasmids (8 runs in parallel),
# then per-plasmid and per-bin tables (surgery/by_plasmid.tsv, surgery/by_bin.tsv)
conda activate plasmid-surgery
PLSDB=$PLSDB scripts/run_surgery.sh 8
python scripts/summarise_surgery.py
```

`$PLSDB` is a PLSDB 2024_05_31_v2 download (`nuccore.csv`, `biosample.csv`, `amr.tsv`, ...); the scripts
default to `../data/plsdb`.

Settings (paths, env names, download threads) are environment variables, listed at the top of
`process_isolates.sh`. `fetch_data.py` is a copy of mash-cluster's `validation/public/fetch_data.py`.

**Checkpoints.** Each step leaves `done/<biosample>.<step>` (`truth`, `shovill`, `unicycler`,
`dragonflye`, `mob_recon.<assembler>`), and `done/<biosample>` once all are done. After an interruption,
run the same command again: finished isolates and steps are skipped; an interrupted step is redone,
and the reads are downloaded again if an assembly is still missing. To redo a step, delete its marker.
To stop cleanly, kill the script's process group (`kill -- -<pid>`) so its children stop too.
To add isolates, append rows to `isolates.tsv` (or select a new file with `--exclude` of these).

Free space is checked on `/mnt/c` (WSL2: the Linux disk is a virtual disk file on C:, so `df .` shows
its nominal 1 TB; C: had ~155 GB free), and the loop stops before an isolate if under `MIN_FREE_GB`
(30). The virtual disk file doesn't shrink when files are deleted (space is reused inside WSL), so it
grows to the peak use. Peak disk use is one isolate's reads (up to ~4 GB gzipped, briefly more while fasterq-dump
converts) plus its assembly working files.

## Notes

- dragonflye estimates genome size from the long reads; on R9.4.1 reads it gave 1.7 Mb for a 5.9 Mb
  genome, then "subsampled" a 57x run as if it were 185x, discarding ~20% of reads. So `--gsize` comes
  from shovill's Illumina 21-mer estimate.
- Hybrid assemblies are not all complete: SAMN08133260 (57x R9.4.1) gave a linear chromosome and
  8 contigs. They're a better-but-imperfect tier, not a second truth.
- SPAdes 3.13.0 segfaults on this machine (also seen for cpo-genomics-db's fixtures).
- 2026-11-07: a WSL2 host-backed filesystem stopped answering flushes, so every `sync()` hung
  (SPAdes calls it: shovill stalled at K51 with no CPU use; `wchan` request_wait_answer, syscall 162).
  `process_isolates.sh` now preloads `scripts/nosync/libnosync.so` (a no-op `sync()`; `NOSYNC=0` to
  disable). Fixed for good by `wsl --shutdown`. It also reuses reads from a completed download
  after an interruption.
- 2026-11-07: Unicycler resumes from an existing output dir, and a run cut short (PC slept) left a
  half-made SPAdes dir that broke the resume (SAMN14640320: "SPAdes failed to produce an assembly
  graph"). The script now clears it first; SAMN14640320's Unicycler was redone after the main pass.

## Layout

```
isolates.tsv, plasmids.tsv, runs.tsv
scripts/            select_isolates.py, process_isolates.sh, fetch_data.py, run_surgery.sh,
                    summarise_surgery.py; envs/*.yml; nosync/ (build: gcc -shared -fPIC -o libnosync.so nosync.c)
truth/<bs>.fasta    published plasmids; truth.tsv.part (rows per plasmid)
assemblies/{shovill,unicycler,dragonflye}/<bs>/contigs.fa (+ .gfa, logs)
mob_recon/{shovill,unicycler,dragonflye}/<bs>/
done/               per-step markers (<bs>.<step>; <bs> = all done)
logs/               progress.log, per-step tool logs, reads_quality.tsv.part
surgery/            plasmid-surgery outputs per <assembler>/<bs>/<bin>/; by_plasmid.tsv, by_bin.tsv
```
