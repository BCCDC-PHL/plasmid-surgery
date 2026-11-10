#!/bin/bash
# Rerun plasmid-surgery with each isolate's reads, one isolate at a time (storage is limited):
# download the reads again, run every reference set (truth, near-relative tiers from all of PLSDB and
# from the carbapenemase-plasmid database) with --short-reads (and --long-reads for dragonflye's
# assemblies), then delete the reads. Assemblies and mob_recon outputs are reused, not redone.
# Outputs mirror the runs without reads, under with_reads/:
#   with_reads/surgery/, with_reads/surgery_near/<set>/, with_reads/near-neighbor-with-carbapenemase/surgery/<set>/
# Read depth is cached in with_reads/cache/ (so each isolate's reads are aligned once per reference).
# An isolate is skipped when done/<bs>.reads exists.
#
# Usage: scripts/rerun_with_reads.sh [jobs] [biosample ...]   (default: every isolate in isolates.tsv)
# Settings (environment): PLSDB, FETCH, DL_THREADS, ENV_SHOVILL (for fetch_data.py), ENV_SRA,
#   FREE_PATH, MIN_FREE_GB, PRIME_THREADS, as in process_isolates.sh / run_surgery.sh;
#   SETS (near-relative sets) [nearest d002 d005 d010 d020]
set -uo pipefail
JOBS=${1:-8}; shift || true
HERE=$(pwd)
SCRIPTS=$(dirname "$(readlink -f "$0")")
PLSDB=$(readlink -f "${PLSDB:-$HERE/../data/plsdb}")
FETCH=${FETCH:-$SCRIPTS/fetch_data.py}
DL_THREADS=${DL_THREADS:-8}
ENV_SHOVILL=${ENV_SHOVILL:-shovill-1.1.0}
ENV_SRA=${ENV_SRA:-sra-tools}
FASTERQ=${FASTERQ:-$(conda run -n "$ENV_SRA" which fasterq-dump)}
SETS=${SETS:-nearest d002 d005 d010 d020}
FREE_PATH=${FREE_PATH:-$( [ -d /mnt/c ] && echo /mnt/c || echo . )}
MIN_FREE_GB=${MIN_FREE_GB:-30}
export PLSDB PRIME_THREADS=${PRIME_THREADS:-16} READS_CACHE=with_reads/cache
mkdir -p done logs reads with_reads/cache

log() { echo "$(date +%FT%T) $*" | tee -a logs/reads_progress.log; }
col() {  # col <name> <biosample>: a field of isolates.tsv
  awk -F'\t' -v name="$1" -v bs="$2" 'NR==1 { for (i=1; i<=NF; i++) c[$i]=i; next } $1==bs { print $c[name] }' isolates.tsv
}
free_gb() { timeout 10 df -B1G --output=avail "$FREE_PATH" 2>/dev/null | tail -1 | tr -d ' '; }

isolates=("$@"); [ ${#isolates[@]} -gt 0 ] || mapfile -t isolates < <(tail -n +2 isolates.tsv | cut -f1)
for bs in "${isolates[@]}"; do
  [ -e "done/$bs.reads" ] && { log "$bs: done already"; continue; }
  free=$(free_gb)
  if [ -n "$free" ] && [ "$free" -lt "$MIN_FREE_GB" ]; then
    log "$bs: STOPPING: only ${free} GB free at $FREE_PATH (MIN_FREE_GB=$MIN_FREE_GB)"; break
  fi
  ill=$(col illumina_run "$bs"); ont=$(col ont_run "$bs"); chem=$(col ont_chemistry "$bs")
  t0=$(date +%s); dl=reads/$bs
  log "$bs: start (illumina $ill${ont:+; ont $ont $chem}); ${free:-?} GB free"
  if ! [ -s "$dl/reads/reads.tsv" ]; then
    rm -rf "$dl"; mkdir -p "$dl"
    awk -F'\t' -v a="$ill" -v b="$ont" 'NR==1 || $3==a || (b!="" && $3==b)' runs.tsv > "$dl/subset.tsv"
    if ! conda run -n "$ENV_SHOVILL" python -u "$FETCH" --cache .cache download -i "$dl/subset.tsv" -o "$dl" \
          --plsdb "$PLSDB" --plasmids plasmids.tsv --fasterq-dump "$FASTERQ" --threads "$DL_THREADS" \
          > "logs/$bs.reads_download.log" 2>&1; then
      log "$bs: download FAILED (logs/$bs.reads_download.log); skipping"; rm -rf "$dl"; continue
    fi
  fi
  r1=$(ls "$dl"/reads/$bs/${ill}_1.fastq.gz); r2=$(ls "$dl"/reads/$bs/${ill}_2.fastq.gz)
  lr=""; [ -n "$ont" ] && lr=$(ls "$dl"/reads/$bs/${ont}*.fastq.gz | head -1)
  preset=map-ont; [ "$chem" = "R10.4.1" ] && preset=lr:hq
  log "$bs: reads in $(( $(date +%s) - t0 )) s ($(du -sh "$dl" | cut -f1))"

  out=logs/$bs.reads_runs.out; : > "$out"
  export ONLY_BS=$bs READS_SHORT="$r1 $r2" READS_LONG="$lr" LONG_PRESET=$preset
  REFSET=truth OUT=with_reads/surgery "$SCRIPTS/run_surgery.sh" "$JOBS" >> "$out" 2>&1
  for s in $SETS; do
    [ -s "refsets/$s/$bs.fasta" ] && REFSET=$s OUT=with_reads/surgery_near/$s "$SCRIPTS/run_surgery.sh" "$JOBS" >> "$out" 2>&1
    [ -s "near-neighbor-with-carbapenemase/refsets/$s/$bs.fasta" ] && REFSETS=near-neighbor-with-carbapenemase/refsets \
      REFSET=$s OUT=with_reads/near-neighbor-with-carbapenemase/surgery/$s "$SCRIPTS/run_surgery.sh" "$JOBS" >> "$out" 2>&1
  done
  n_ok=$(grep -c '^ok' "$out"); n_fail=$(grep -c '^FAILED' "$out")
  rm -rf "$dl"
  if [ "$n_fail" = 0 ]; then
    touch "done/$bs.reads"; log "$bs: $n_ok runs ok in $(( $(date +%s) - t0 )) s; reads deleted"
  else
    log "$bs: $n_fail runs FAILED, $n_ok ok ($out); reads deleted; a rerun redoes the failed runs"
  fi
done
