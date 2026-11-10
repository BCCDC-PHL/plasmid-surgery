#!/usr/bin/env python3
"""Pick near-relative references for each truth plasmid: PLSDB plasmids that aren't the isolate's own.

Input: mash distances from every truth plasmid to the PLSDB sketch, e.g.
  cat truth/*.fasta > near/truth_all.fasta
  mash dist -p 16 -d 0.05 $PLSDB/plsdb_sketch.msh -i near/truth_all.fasta > near/mash_truth_vs_plsdb.tsv

Excluded as candidates: plasmids from the isolate's own BioSample (incl. the truth itself), and plasmids
PLSDB records as identical to a truth plasmid (nucc_identical.csv).

Reference sets (one reference per truth plasmid, where a candidate exists):
  nearest   the smallest mash distance (what a database search would most likely find)
  d002      distance in (0.001, 0.003]    }
  d005      (0.003, 0.0075]               }  the candidate nearest the tier's centre (geometric),
  d010      (0.0075, 0.015]               }  length within 0.5-2x of the truth plasmid
  d020      (0.015, 0.03]                 }

Writes:
  near/near_refs.tsv               one row per (set, truth plasmid): chosen reference, distance, length
  refsets/<set>/<bs>.fasta         the isolate's references for that set (deduplicated)
  refsets/<set>/<bs>.bed           their AMRFinderPlus genes from PLSDB amr.tsv (--ref-features)
Sequences come from PLSDB's sequences.fasta.bz2 (one streaming pass).
"""
import argparse
import bz2
import csv
import math
import os
import sys
from collections import defaultdict

TIERS = [('d002', 0.001, 0.003), ('d005', 0.003, 0.0075), ('d010', 0.0075, 0.015), ('d020', 0.015, 0.03)]
csv.field_size_limit(sys.maxsize)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--plsdb', required=True)
    p.add_argument('--mash', default='near/mash_truth_vs_plsdb.tsv')
    p.add_argument('--truth', default='truth/truth.tsv.part')
    p.add_argument('--out', default='refsets')
    p.add_argument('--table', default='near/near_refs.tsv')
    a = p.parse_args()

    truth = {}   # acc -> (bs, length, carbapenemases)
    with open(a.truth) as f:
        for row in csv.reader(f, delimiter='\t'):
            bs, acc, length, _, _, carb, _ = (row + [''] * 7)[:7]
            truth[acc] = (bs, int(length), carb)

    nuccore = {}  # acc -> (length, biosample_uid, description)
    with open(os.path.join(a.plsdb, 'nuccore.csv'), newline='') as f:
        for r in csv.DictReader(f):
            nuccore[r['NUCCORE_ACC']] = (int(r['NUCCORE_Length']), r['BIOSAMPLE_UID'], r['NUCCORE_Description'])
    # PLSDB drops exact duplicates, keeping one representative: a group = representative + its duplicates
    identical = defaultdict(set)
    with open(os.path.join(a.plsdb, 'nucc_identical.csv'), newline='') as f:
        for r in csv.DictReader(f):
            identical[r['NUCCORE_Identical']].add(r['NUCCORE_ACC'])
            identical[r['NUCCORE_ACC']].add(r['NUCCORE_Identical'])

    own_uids = defaultdict(set)   # bs -> BioSample UIDs of its truth plasmids
    for acc, (bs, _, _) in truth.items():
        if acc in nuccore and nuccore[acc][1] not in ('', '-1'):
            own_uids[bs].add(nuccore[acc][1])

    cands = defaultdict(list)     # truth acc -> [(dist, shared, ref acc)]
    with open(a.mash) as f:
        for ref, query, dist, _, shared in csv.reader(f, delimiter='\t'):
            if query not in truth or ref not in nuccore:
                continue
            bs = truth[query][0]
            if ref in truth and truth[ref][0] == bs:
                continue
            if nuccore[ref][1] in own_uids[bs]:
                continue
            if ref in identical[query] or any(ref in identical[x] for x in identical[query]):
                continue
            cands[query].append((float(dist), int(shared.split('/')[0]), ref))

    rows = []
    for acc, (bs, length, carb) in sorted(truth.items(), key=lambda x: (x[1][0], x[0])):
        cs = sorted(cands[acc], key=lambda c: (c[0], -c[1]))
        picks = []
        if cs:
            picks.append(('nearest', cs[0]))
        for name, lo, hi in TIERS:
            centre = math.sqrt(lo * hi)
            ok = [c for c in cs if lo < c[0] <= hi and 0.5 <= nuccore[c[2]][0] / length <= 2]
            if ok:
                picks.append((name, min(ok, key=lambda c: abs(math.log(max(c[0], 1e-6) / centre)))))
        for name, (dist, shared, ref) in picks:
            rows.append(dict(set=name, biosample=bs, truth_acc=acc, truth_length=length,
                             truth_carbapenemases=carb, ref_acc=ref, mash_dist=dist, shared_hashes=shared,
                             ref_length=nuccore[ref][0], ref_description=nuccore[ref][2]))
        if not cs:
            rows.append(dict(set='none', biosample=bs, truth_acc=acc, truth_length=length,
                             truth_carbapenemases=carb))

    os.makedirs(os.path.dirname(a.table), exist_ok=True)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with open(a.table, 'w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', restval='', lineterminator='\n')
        w.writeheader(); w.writerows(rows)

    # Sequences: one pass through PLSDB's FASTA
    wanted = {r['ref_acc'] for r in rows if r['set'] != 'none'}
    seqs, keep, name = {}, False, None
    with bz2.open(os.path.join(a.plsdb, 'sequences.fasta.bz2'), 'rt') as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()[0]
                keep = name in wanted
                if keep:
                    seqs[name] = [line]
            elif keep:
                seqs[name].append(line)
    missing = wanted - seqs.keys()
    if missing:
        print(f'warning: {len(missing)} references not in sequences.fasta.bz2: {sorted(missing)[:5]}', file=sys.stderr)

    by_set = defaultdict(lambda: defaultdict(list))   # set -> bs -> [ref acc] (deduplicated, in order)
    for r in rows:
        if r['set'] != 'none' and r['ref_acc'] in seqs and r['ref_acc'] not in by_set[r['set']][r['biosample']]:
            by_set[r['set']][r['biosample']].append(r['ref_acc'])
    beds = defaultdict(list)
    accs = {acc for s in by_set.values() for refs in s.values() for acc in refs}
    with open(os.path.join(a.plsdb, 'amr.tsv'), newline='') as f:
        for row in csv.reader(f, delimiter='\t'):
            if row[0] in accs and row[1] == 'amrfinderplus':
                beds[row[0]].append((row[0], int(row[8]) - 1, int(row[9]), row[3], '0', row[11]))
    for s, by_bs in by_set.items():
        os.makedirs(os.path.join(a.out, s), exist_ok=True)
        for bs, refs in by_bs.items():
            with open(os.path.join(a.out, s, f'{bs}.fasta'), 'w') as f:
                for acc in refs:
                    f.writelines(seqs[acc])
            with open(os.path.join(a.out, s, f'{bs}.bed'), 'w') as f:
                for b in sorted({b for acc in refs for b in beds[acc]}):
                    f.write('\t'.join(map(str, b)) + '\n')

    # Overview
    n_truth = len(truth)
    print(f'{n_truth} truth plasmids; {sum(1 for t in truth if not cands[t])} with no candidate within the mash cutoff')
    for s in ['nearest'] + [t[0] for t in TIERS]:
        rs = [r for r in rows if r['set'] == s]
        ds = sorted(r['mash_dist'] for r in rs)
        med = ds[len(ds) // 2] if ds else 'NA'
        print(f'{s:8s} {len(rs):3d} truth plasmids covered, {len(by_set[s])} isolates, median distance {med}')


if __name__ == '__main__':
    main()
