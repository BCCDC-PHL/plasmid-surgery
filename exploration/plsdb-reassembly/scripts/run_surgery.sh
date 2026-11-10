#!/bin/bash
# Run plasmid-surgery on every mob_recon plasmid bin of every assembly, against that isolate's own
# PLSDB plasmids (truth/<bs>.fasta; every plasmid a separate reference), with the whole assembly
# (--assembly) and all mob_recon bins incl. chromosome (--bins), and the truth plasmids' AMRFinderPlus
# genes from PLSDB amr.tsv (--ref-features, as BED). Outputs: surgery/<assembler>/<bs>/<bin>/.
# With REFSET=<set>, the references are refsets/<set>/<bs>.fasta instead (near relatives, from
# scripts/select_near_refs.py, with their features in refsets/<set>/<bs>.bed); outputs go to
# surgery_near/<set>/ and isolates without a reference in that set are skipped.
# A finished run (assessment.tsv present) is skipped. Usage: [REFSET=<set>] scripts/run_surgery.sh [jobs]
set -euo pipefail
cd "$(dirname "$0")/.."
JOBS=${1:-8}
PLSDB=${PLSDB:-../data/plsdb}
ENV_PS=${ENV_PS:-plasmid-surgery}
REFSET=${REFSET:-truth}  # truth, or a set in $REFSETS/
REFSETS=${REFSETS:-refsets}   # e.g. near-neighbor-with-carbapenemase/refsets (then set OUT too)
if [ "$REFSET" = truth ]; then
  REFDIR=truth; FEATDIR=surgery/features; OUT=${OUT:-surgery}   # features stay in surgery/features
else
  REFDIR=$REFSETS/$REFSET; FEATDIR=$REFDIR; OUT=${OUT:-surgery_near/$REFSET}
  [ -d "$REFDIR" ] || { echo "no reference set $REFDIR (scripts/select_near_refs.py)" >&2; exit 1; }
fi
THREADS=${THREADS:-1}   # plasmid-surgery -t
mkdir -p surgery/features "logs/$OUT"

# Ref features (truth set; select_near_refs.py writes them for the near-relative sets): AMRFinderPlus hits on each isolate's truth plasmids -> BED (name = gene symbol)
for f in truth/*.fasta; do
  bs=$(basename "$f" .fasta); bed=surgery/features/$bs.bed
  [ -s "$bed" ] && continue
  grep '^>' "$f" | sed 's/^>//; s/ .*//' > surgery/features/$bs.ids
  awk -F'\t' 'NR==FNR {ids[$1]; next} ($1 in ids) && $2=="amrfinderplus" {print $1"\t"$9-1"\t"$10"\t"$4"\t0\t"$12}' \
    surgery/features/$bs.ids "$PLSDB/amr.tsv" | sort -k1,1 -k2,2n -u > "$bed"
  rm surgery/features/$bs.ids
done

# One job per (assembler, isolate, bin)
jobs_list() {
  for d in mob_recon/*/*/; do
    a=$(basename "$(dirname "$d")"); bs=$(basename "$d")
    for bin in "$d"plasmid_*.fasta; do
      [ -e "$bin" ] || continue
      [ -s "$REFDIR/$bs.fasta" ] || continue
      echo "$a $bs $(basename "$bin" .fasta)"
    done
  done
}
run_one() {
  a=$1; bs=$2; b=$3; out=$OUT/$a/$bs/$b
  [ -s "$out/assessment.tsv" ] && return 0
  feat=(); [ -s "$FEATDIR/$bs.bed" ] && feat=(--ref-features "$FEATDIR/$bs.bed")
  if conda run -n "$ENV_PS" plasmid-surgery --plasmid-reconstruction "mob_recon/$a/$bs/$b.fasta" \
       --ref-plasmid "$REFDIR/$bs.fasta" --assembly "assemblies/$a/$bs/contigs.fa" \
       --bins mob_recon/$a/$bs/plasmid_*.fasta mob_recon/$a/$bs/chromosome.fasta "${feat[@]}" \
       --sample-id "${a}_${bs}_${b}" -t "$THREADS" --report -o "$out" > "logs/$OUT/${a}_${bs}_${b}.log" 2>&1; then
    echo "ok $a $bs $b"
  else
    echo "FAILED $a $bs $b (logs/$OUT/${a}_${bs}_${b}.log)"
  fi
}
export -f run_one; export ENV_PS OUT THREADS REFDIR FEATDIR
jobs_list | xargs -P "$JOBS" -L1 bash -c 'run_one "$@"' _
