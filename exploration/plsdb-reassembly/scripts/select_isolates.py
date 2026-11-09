#!/usr/bin/env python3
"""
Choose isolates (BioSamples) to reassemble, from fetch_data.py's runs.tsv.

Each isolate needs a paired-end Illumina run; some also get a nanopore run (for a hybrid assembly).
Isolates are spread over carbapenemase families, with quotas proportional to the square root of
each family's number of candidates (rare families pooled), chosen at random within a family, with
at most --max-per-species isolates of one species.

Usage: select_isolates.py --plsdb DIR --plasmids plasmids.tsv --runs runs.tsv \
         --hybrid 15 --short 10 --exclude SAMN... -o isolates.tsv
"""
from __future__ import annotations

import argparse
import csv
import math
import random
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

csv.field_size_limit(sys.maxsize)

FIELDS = ['biosample', 'species', 'families', 'carbapenemases', 'carbapenemase_plasmids', 'plasmid_lengths',
          'n_plsdb_plasmids', 'mode', 'illumina_run', 'illumina_model', 'illumina_bases', 'illumina_bytes',
          'ont_run', 'ont_model', 'ont_chemistry', 'ont_bases', 'ont_bytes', 'other_long_read_platforms']


def read(path, delimiter=','):
    with open(path, newline='') as f:
        yield from csv.DictReader(f, delimiter=delimiter)


def family(gene: str) -> str:
    """blaKPC-2 -> KPC; OXA genes by group (OXA-48-like vs others), since they behave differently."""
    g = re.sub(r'^bla', '', gene)
    m = re.match(r'OXA-(\d+)', g)
    if m:
        n = int(m.group(1))
        return 'OXA-48-like' if n in {48, 162, 181, 204, 232, 244, 245, 370, 484, 515, 519} else 'OXA-other'
    return re.split(r'[-_]', g)[0]


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('--plsdb', required=True)
    p.add_argument('--plasmids', required=True)
    p.add_argument('--runs', required=True)
    p.add_argument('--hybrid', type=int, default=15, help="isolates with Illumina + nanopore")
    p.add_argument('--short', type=int, default=10, help="isolates with Illumina only (no nanopore run)")
    p.add_argument('--max-illumina-gb', type=float, default=1.5)
    p.add_argument('--max-ont-gb', type=float, default=2.5)
    p.add_argument('--min-illumina-bases', type=float, default=3e8, help="at least ~60x of a 5 Mb genome")
    p.add_argument('--max-per-species', type=int, default=3)
    p.add_argument('--exclude', nargs='*', default=[])
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('-o', '--output', required=True)
    args = p.parse_args(argv)
    plsdb = Path(args.plsdb)

    species = {r['TAXONOMY_UID']: r['TAXONOMY_taxon_name'] for r in read(plsdb / 'taxonomy.csv')}
    uid_to_acc = {r['BIOSAMPLE_UID']: r['BIOSAMPLE_ACC'] for r in read(plsdb / 'biosample.csv')}
    bs_species, n_plasmids = {}, Counter()
    for r in read(plsdb / 'nuccore.csv'):
        bs = uid_to_acc.get(r['BIOSAMPLE_UID'])
        if bs:
            n_plasmids[bs] += 1
            bs_species.setdefault(bs, species.get(r['TAXONOMY_UID'], ''))
    carrying = defaultdict(list)
    for r in read(args.plasmids, '\t'):
        if r['biosample']:
            carrying[r['biosample']].append(r)

    runs = defaultdict(list)
    for r in read(args.runs, '\t'):
        runs[r['biosample']].append(r)

    def size(r):
        return int(r['fastq_bytes'] or 0)

    candidates = {'hybrid': [], 'short': []}
    for bs, found in runs.items():
        if bs in args.exclude or bs not in carrying:
            continue
        ill = [r for r in found if r['platform'] == 'ILLUMINA' and r['library_layout'] == 'PAIRED'
               and r['library_strategy'] == 'WGS' and 0 < size(r) <= args.max_illumina_gb * 1e9
               and float(r['base_count'] or 0) >= args.min_illumina_bases]
        if not ill:
            continue
        ont = [r for r in found if r['platform'] == 'OXFORD_NANOPORE' and 0 < size(r) <= args.max_ont_gb * 1e9
               and float(r['base_count'] or 0) >= 1e8]
        other_long = sorted({r['platform'] for r in found if r['platform'] in ('PACBIO_SMRT',)})
        plasmids = carrying[bs]
        genes = sorted({g for r in plasmids for g in r['carbapenemases'].split(',') if g})
        row = {'biosample': bs, 'species': bs_species.get(bs, ''),
               'families': ','.join(sorted({family(g) for g in genes})), 'carbapenemases': ','.join(genes),
               'carbapenemase_plasmids': ','.join(r['nuccore_acc'] for r in plasmids),
               'plasmid_lengths': ','.join(r['length'] for r in plasmids), 'n_plsdb_plasmids': n_plasmids[bs],
               'other_long_read_platforms': ','.join(other_long)}
        i = min(ill, key=lambda r: abs(float(r['base_count']) - 1e9))   # nearest ~200x of 5 Mb
        row.update(illumina_run=i['run_accession'], illumina_model=i['instrument_model'],
                   illumina_bases=i['base_count'], illumina_bytes=size(i))
        if ont:
            o = max(ont, key=size)
            row.update(mode='hybrid', ont_run=o['run_accession'], ont_model=o['instrument_model'],
                       ont_chemistry=o['chemistry'], ont_bases=o['base_count'], ont_bytes=size(o))
            candidates['hybrid'].append(row)
        else:
            row.update(mode='short', ont_run='', ont_model='', ont_chemistry='', ont_bases='', ont_bytes='')
            candidates['short'].append(row)

    rng = random.Random(args.seed)
    chosen = []
    species_count = Counter()
    for mode, n in (('hybrid', args.hybrid), ('short', args.short)):
        pool = candidates[mode]
        rng.shuffle(pool)
        by_family = defaultdict(list)
        for r in pool:
            by_family[r['families'].split(',')[0]].append(r)
        # Families with under 2% of the candidates are pooled as 'rare'
        for f in [f for f in by_family if len(by_family[f]) < 0.02 * len(pool)]:
            by_family['rare'] += by_family.pop(f)
        # Quotas proportional to sqrt(family size) (largest remainder), at least one each
        weights = {f: math.sqrt(len(v)) for f, v in by_family.items()}
        total = sum(weights.values())
        quota = {f: max(1, int(n * w / total)) for f, w in weights.items()}
        for f in sorted(weights, key=lambda f: -(n * weights[f] / total % 1)):
            if sum(quota.values()) >= n:
                break
            quota[f] += 1
        picked = []
        for f in sorted(quota, key=lambda f: -len(by_family[f])):
            for r in by_family[f]:
                if sum(1 for x in picked if x['_family'] == f) >= quota[f] or len(picked) >= n:
                    break
                if species_count[r['species']] >= args.max_per_species:
                    continue
                picked.append({**r, '_family': f})
                species_count[r['species']] += 1
        chosen += picked
        print(f"{mode}: {len(pool)} candidates; quotas " +
              ', '.join(f"{f} {q}" for f, q in sorted(quota.items(), key=lambda kv: -kv[1])) +
              f"; chose {len(picked)}", file=sys.stderr)

    with open(args.output, 'w', newline='') as f:
        w = csv.DictWriter(f, FIELDS, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(chosen)


if __name__ == '__main__':
    main()
