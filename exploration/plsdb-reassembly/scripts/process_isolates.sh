#!/bin/bash
# Reassemble each isolate in isolates.tsv, one at a time. For each isolate:
#   1. download its reads (and the published sequences of all its PLSDB plasmids: the truth)
#   2. shovill 1.1.0 (SPAdes 4.3.0), defaults                                   -> assemblies/shovill/
#   3. fastp 1.4.0 (routine-assembly's options) + Unicycler 0.5.1 (SPAdes 4.3.0) -> assemblies/unicycler/
#   4. with a nanopore run: dragonflye 1.2.1 hybrid (Flye; racon x1; Polypolish x1) -> assemblies/dragonflye/
#   5. delete the reads
#   6. mob_recon 3.1.9 on each assembly                                         -> mob_recon/<assembler>/
# Each step leaves done/<biosample>.<step>, and is skipped when that exists, so a rerun only does
# what's missing (reads are downloaded only when an assembly is missing). done/<biosample> = all done.
# Run from the directory holding isolates.tsv, runs.tsv and plasmids.tsv; outputs go there too.
#
# Usage: scripts/process_isolates.sh [isolates.tsv] [threads] [ram_gb]
#
# Settings, as environment variables (defaults in brackets):
#   PLSDB          PLSDB metadata directory, with nuccore.csv and biosample.csv [../data/plsdb]
#   FETCH          fetch_data.py [next to this script]
#   FASTERQ        fasterq-dump [the one in env $ENV_SRA]
#   DL_THREADS     threads for downloads and fasterq-dump [8]
#   ENV_SHOVILL, ENV_UNICYCLER, ENV_DRAGONFLYE, ENV_MOB, ENV_SRA  conda env names
#                  [shovill-1.1.0, unicycler-0.5.1, dragonflye-1.2.1, mob-suite, sra-tools]
#                  (environment files in scripts/envs/)
#   FREE_PATH      where to check free space [/mnt/c if it exists, else .]: on WSL2 the Linux disk is a
#                  virtual disk file on C:, so df of . shows its nominal size, not C:'s real free space
#   MIN_FREE_GB    don't start an isolate with less than this free at FREE_PATH [30]
#   NOSYNC         1: preload scripts/nosync/libnosync.so, making sync() a no-op [1]. On WSL2 a
#                  host-backed filesystem can stop answering flushes, and then every sync() (SPAdes
#                  calls it) blocks forever. Build it with: gcc -shared -fPIC -O2 -o libnosync.so nosync.c
set -uo pipefail

ISOLATES=${1:-isolates.tsv}
THREADS=${2:-16}
RAM=${3:-32}
HERE=$(pwd)
SCRIPTS=$(dirname "$(readlink -f "$0")")
PLSDB=$(readlink -f "${PLSDB:-$HERE/../data/plsdb}")
FETCH=${FETCH:-$SCRIPTS/fetch_data.py}
DL_THREADS=${DL_THREADS:-8}
ENV_SHOVILL=${ENV_SHOVILL:-shovill-1.1.0}
ENV_UNICYCLER=${ENV_UNICYCLER:-unicycler-0.5.1}
ENV_DRAGONFLYE=${ENV_DRAGONFLYE:-dragonflye-1.2.1}
ENV_MOB=${ENV_MOB:-mob-suite}
ENV_SRA=${ENV_SRA:-sra-tools}
FASTERQ=${FASTERQ:-$(conda run -n "$ENV_SRA" which fasterq-dump)}
if [ "${NOSYNC:-1}" = 1 ]; then
  [ -e "$SCRIPTS/nosync/libnosync.so" ] || { echo "missing $SCRIPTS/nosync/libnosync.so (or set NOSYNC=0)" >&2; exit 1; }
  export LD_PRELOAD="$SCRIPTS/nosync/libnosync.so${LD_PRELOAD:+:$LD_PRELOAD}"
fi
for f in "$ISOLATES" runs.tsv plasmids.tsv "$PLSDB/nuccore.csv" "$FETCH" "$FASTERQ"; do
  [ -e "$f" ] || { echo "missing: $f" >&2; exit 1; }
done
mkdir -p done logs reads truth assemblies mob_recon

FREE_PATH=${FREE_PATH:-$( [ -d /mnt/c ] && echo /mnt/c || echo . )}
MIN_FREE_GB=${MIN_FREE_GB:-30}
free_gb() { timeout 10 df -B1G --output=avail "$FREE_PATH" 2>/dev/null | tail -1 | tr -d ' '; }

log() { echo "$(date +%FT%T) $*" | tee -a logs/progress.log; }
elapsed() { echo $(( $(date +%s) - $1 )); }
have() { [ -e "done/$bs.$1" ]; }
mark() { touch "done/$bs.$1"; }

col() {  # col <name> <tsv line>: a field of isolates.tsv by column name
  awk -F'\t' -v name="$1" -v line="$2" 'NR==1 { for (i=1; i<=NF; i++) c[$i]=i }
    END { split(line, f, "\t"); print f[c[name]] }' "$ISOLATES"
}

while IFS= read -r line <&3; do
  bs=$(col biosample "$line"); mode=$(col mode "$line")
  ill=$(col illumina_run "$line"); ont=$(col ont_run "$line"); chem=$(col ont_chemistry "$line")
  [ -e "done/$bs" ] && { log "$bs: done already"; continue; }
  assemblers="shovill unicycler"; [ -n "$ont" ] && assemblers="$assemblers dragonflye"
  missing=""; for a in $assemblers; do have "$a" || missing="$missing $a"; done
  free=$(free_gb)
  if [ -n "$free" ] && [ "$free" -lt "$MIN_FREE_GB" ]; then
    log "$bs: STOPPING: only ${free} GB free at $FREE_PATH (MIN_FREE_GB=$MIN_FREE_GB)"; break
  fi
  log "$bs: start ($mode; illumina $ill${ont:+; ont $ont $chem}; to assemble:${missing:- nothing}); ${free:-?} GB free at $FREE_PATH"
  t0=$(date +%s)
  dl=reads/$bs

  if [ -n "$missing" ] && [ -s "$dl/reads/reads.tsv" ] && [ -s "$dl/truth/$bs.fasta" ]; then
    log "$bs: reusing the reads already downloaded"   # a completed download from an interrupted run
    r1=$(ls "$dl"/reads/$bs/${ill}_1.fastq.gz); r2=$(ls "$dl"/reads/$bs/${ill}_2.fastq.gz)
  elif [ -n "$missing" ]; then
    # 1. Reads, and the truth
    rm -rf "$dl"; mkdir -p "$dl"
    awk -F'\t' -v a="$ill" -v b="$ont" 'NR==1 || $3==a || (b!="" && $3==b)' runs.tsv > "$dl/subset.tsv"
    if ! conda run -n "$ENV_SHOVILL" python -u "$FETCH" --cache .cache download -i "$dl/subset.tsv" -o "$dl" \
          --plsdb "$PLSDB" --plasmids plasmids.tsv --fasterq-dump "$FASTERQ" --threads "$DL_THREADS" \
          > "logs/$bs.download.log" 2>&1; then
      log "$bs: download FAILED (logs/$bs.download.log); skipping"; rm -rf "$dl"; continue
    fi
    r1=$(ls "$dl"/reads/$bs/${ill}_1.fastq.gz); r2=$(ls "$dl"/reads/$bs/${ill}_2.fastq.gz)
    log "$bs: reads in $(elapsed $t0) s ($(du -sh "$dl" | cut -f1))"
  fi
  if [ -n "$missing" ] && ! have truth; then
    cp "$dl/truth/$bs.fasta" truth/
    tail -n +2 "$dl/truth/truth.tsv" >> truth/truth.tsv.part
    tail -n +2 "$dl/reads/reads.tsv" >> logs/reads_quality.tsv.part
    mark truth
  fi

  # 2. shovill
  if ! have shovill; then
    t1=$(date +%s); out=assemblies/shovill/$bs
    if conda run -n "$ENV_SHOVILL" shovill --R1 "$r1" --R2 "$r2" --outdir "$out" \
          --cpus "$THREADS" --ram "$RAM" --force > "logs/$bs.shovill.log" 2>&1; then
      ( cd "$out" && rm -rf spades.fasta flash.* *.fq.gz spades contigs.fa.fai 2>/dev/null; true )
      mark shovill; log "$bs: shovill ok in $(elapsed $t1) s, $(grep -c '>' "$out/contigs.fa") contigs"
    else
      log "$bs: shovill FAILED (logs/$bs.shovill.log)"
    fi
  fi

  # 3. fastp (with routine-assembly's options) + Unicycler, defaults
  if ! have unicycler; then
    t1=$(date +%s); out=assemblies/unicycler/$bs; tmp=$dl/fastp
    # Unicycler resumes from an existing output dir; a run cut short leaves a half-made SPAdes dir that
    # breaks the resume (SPAdes failed to produce an assembly graph), so always start clean.
    rm -rf "$out/unicycler"; mkdir -p "$out" "$tmp"
    if conda run -n "$ENV_UNICYCLER" fastp -i "$r1" -I "$r2" -o "$tmp/R1.trim.fastq.gz" -O "$tmp/R2.trim.fastq.gz" \
          --cut_tail --trim_poly_g --overrepresentation_analysis --detect_adapter_for_pe -w "$THREADS" \
          -j "$out/fastp.json" -h "$out/fastp.html" > "logs/$bs.fastp.log" 2>&1 \
       && conda run -n "$ENV_UNICYCLER" unicycler --threads "$THREADS" -1 "$tmp/R1.trim.fastq.gz" -2 "$tmp/R2.trim.fastq.gz" \
          -o "$out/unicycler" > "logs/$bs.unicycler.log" 2>&1; then
      mv "$out/unicycler/assembly.fasta" "$out/contigs.fa"; mv "$out/unicycler/assembly.gfa" "$out/contigs.gfa"
      mv "$out/unicycler/unicycler.log" "$out/unicycler.log"; rm -rf "$out/unicycler"
      mark unicycler; log "$bs: unicycler ok in $(elapsed $t1) s, $(grep -c '>' "$out/contigs.fa") contigs"
    else
      log "$bs: unicycler FAILED (logs/$bs.fastp.log, logs/$bs.unicycler.log)"
    fi
  fi

  # 4. dragonflye hybrid
  if [ -n "$ont" ] && ! have dragonflye; then
    t1=$(date +%s); out=assemblies/dragonflye/$bs
    lr=$(ls "$dl"/reads/$bs/${ont}*.fastq.gz | head -1)
    hq=""; [ "$chem" = "R10.4.1" ] && hq="--nanohq"
    # Genome size from shovill's estimate (Illumina 21-mers): dragonflye's own, from noisy long reads,
    # can be far too small, which inflates its depth estimate and makes it discard reads.
    gs=$(grep -oP 'Using genome size \K[0-9]+' "logs/$bs.shovill.log" | head -1)
    gsize=""; [ -n "$gs" ] && gsize="--gsize $gs"
    if conda run -n "$ENV_DRAGONFLYE" dragonflye --reads "$lr" --R1 "$r1" --R2 "$r2" $hq $gsize --racon 1 --polypolish 1 \
          --outdir "$out" --cpus "$THREADS" --ram "$RAM" --force > "logs/$bs.dragonflye.log" 2>&1; then
      ( cd "$out" && rm -rf flye *.fq.gz *.bam* polypolish* racon* 2>/dev/null; true )
      mark dragonflye; log "$bs: dragonflye ok in $(elapsed $t1) s, $(grep -c '>' "$out/contigs.fa") contigs"
    else
      log "$bs: dragonflye FAILED (logs/$bs.dragonflye.log)"
    fi
  fi

  # 5. Reads no longer needed
  rm -rf "$dl"

  # 6. mob_recon on each assembly
  for a in $assemblers; do
    fa=assemblies/$a/$bs/contigs.fa
    { have "mob_recon.$a" || [ ! -s "$fa" ]; } && continue
    t3=$(date +%s)
    if conda run -n "$ENV_MOB" mob_recon -i "$fa" -o "mob_recon/$a/$bs" -n "$THREADS" --force \
          > "logs/$bs.mob_recon.$a.log" 2>&1; then
      mark "mob_recon.$a"
      log "$bs: mob_recon ($a) ok in $(elapsed $t3) s, $(ls mob_recon/$a/$bs/plasmid_*.fasta 2>/dev/null | wc -l) plasmid bins"
    else
      log "$bs: mob_recon ($a) FAILED (logs/$bs.mob_recon.$a.log)"
    fi
  done

  all=1; for a in $assemblers; do have "$a" && have "mob_recon.$a" || all=0; done
  if [ $all = 1 ]; then touch "done/$bs"; log "$bs: finished in $(elapsed $t0) s"
  else log "$bs: incomplete after $(elapsed $t0) s (see FAILED above); a rerun retries the missing steps"; fi
done 3< <(tail -n +2 "$ISOLATES")
n_done=$(tail -n +2 "$ISOLATES" | cut -f1 | while read -r b; do [ -e "done/$b" ] && echo; done | wc -l)
log "loop ended: $n_done of $(( $(wc -l < "$ISOLATES") - 1 )) isolates done"
