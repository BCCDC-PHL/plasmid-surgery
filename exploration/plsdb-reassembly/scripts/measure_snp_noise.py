#!/usr/bin/env python3
"""False-SNP rate of plasmid-surgery's pairwise SNP distances (plasmid_surgery.compare), on the PLSDB set.

Per isolate, every truth plasmid (separately) and every mob_recon plasmid bin of every assembler that
matches a truth plasmid (surgery/by_bin.tsv) are compared all-vs-all, as plasmid-surgery-cohort does.
Pairs that are the same plasmid should differ by ~0 SNPs:

  <short|hybrid>_vs_truth   a bin vs. the truth plasmid it matches best
  <assembler>_vs_<assembler> two assemblers' bins matching the same truth plasmid

SAMN13951917 (reads don't match its truth) is left out of the *_vs_truth categories.

Writes snp_noise/pairs.tsv (one row per same-plasmid pair) and snp_noise/snps.tsv (each SNP, with its
distance to the nearest contig end on either member), and prints a summary per category.
Usage (from exploration/plsdb-reassembly/, env plasmid-surgery): python scripts/measure_snp_noise.py [threads]
"""
import csv
import os
import sys
import tempfile
from collections import defaultdict

from plasmid_surgery import compare, parsers

EXCLUDE_TRUTH = {'SAMN13951917'}
SHORT = {'shovill', 'unicycler'}


def main(threads=8):
    threads = int(threads)
    bins = defaultdict(list)          # bs -> [(assembler, bin, best_ref)]
    with open('surgery/by_bin.tsv') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['best_ref'] != 'none':
                bins[r['biosample']].append((r['assembler'], r['bin'], r['best_ref']))
    os.makedirs('snp_noise', exist_ok=True)
    pair_rows, snp_rows = [], []
    for bs in sorted(bins):
        members = []
        for rec in parsers.parse_fasta(f'truth/{bs}.fasta'):
            members.append({'name': rec['seq_id'], 'kind': 'truth', 'ref': rec['seq_id'], 'assembler': 'truth',
                            'contigs': [rec]})
        for a, b, ref in bins[bs]:
            contigs = parsers.parse_fasta(f'mob_recon/{a}/{bs}/{b}.fasta')
            members.append({'name': f'{a}:{b}', 'kind': 'bin', 'ref': ref, 'assembler': a, 'contigs': contigs})
        with tempfile.TemporaryDirectory() as tmp:
            result = compare.all_vs_all(members, tmp, threads)
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                mi, mj = members[i], members[j]
                if mi['ref'] != mj['ref'] or mi['assembler'] == mj['assembler']:
                    continue
                pair = sorted([mi['assembler'], mj['assembler']])
                if 'truth' in pair:
                    if bs in EXCLUDE_TRUTH:
                        continue
                    other = pair[0] if pair[1] == 'truth' else pair[1]
                    category = f"{'short' if other in SHORT else 'hybrid'}_vs_truth"
                else:
                    category = f'{pair[0]}_vs_{pair[1]}'
                f = compare.snp_fields(result, i, j)
                pair_rows.append({'biosample': bs, 'category': category, 'truth_plasmid': mi['ref'],
                                  'a': mi['name'], 'b': mj['name'], **f})
                snp = result['snps'][(i, j)]
                for tk, t, qk, q in snp['positions']:
                    ta, qa = mj['contigs'][tk], mi['contigs'][qk]
                    to_end = min(t, len(ta['seq']) - t, q, len(qa['seq']) - q)
                    snp_rows.append({'biosample': bs, 'category': category, 'a': mi['name'], 'a_contig': qa['seq_id'],
                                     'a_pos': q, 'b': mj['name'], 'b_contig': ta['seq_id'], 'b_pos': t,
                                     'bp_to_nearest_contig_end': to_end})
        print(f'{bs}: {len(members)} members', file=sys.stderr)

    for path, rows in (('snp_noise/pairs.tsv', pair_rows), ('snp_noise/snps.tsv', snp_rows)):
        with open(path, 'w', newline='') as f:
            if rows:
                w = csv.DictWriter(f, list(rows[0]), delimiter='\t', lineterminator='\n')
                w.writeheader(); w.writerows(rows)

    print('category\tpairs\treported\tcompared_Mb\tsnps_outside\tper_100kb\tpairs_0_snps\tsnps_in_clusters\t'
          'indels\tindels_per_100kb')
    for cat in sorted({r['category'] for r in pair_rows}):
        rs = [r for r in pair_rows if r['category'] == cat]
        rep = [r for r in rs if r['snp_distance_reported']]
        comp = sum(r['compared_bases'] for r in rs)
        out = sum(r['snps_outside_clusters'] for r in rs)
        ind = sum(r['indels'] for r in rs)
        print('\t'.join(map(str, [cat, len(rs), len(rep), round(comp / 1e6, 2), out,
                                  round(out / comp * 1e5, 3) if comp else 'NA',
                                  sum(r['snps_outside_clusters'] == 0 for r in rep),
                                  sum(r['snps_in_clusters'] for r in rs), ind,
                                  round(ind / comp * 1e5, 3) if comp else 'NA'])))
    print('\nSNPs by distance to the nearest contig end')
    edges = [(100, 500), (500, 2000), (2000, 10 ** 9)]
    for cat in sorted({r['category'] for r in snp_rows}):
        ds = [r['bp_to_nearest_contig_end'] for r in snp_rows if r['category'] == cat]
        print(cat, ' '.join(f"{lo}-{hi if hi < 10 ** 9 else ''}: {sum(lo <= d < hi for d in ds)}" for lo, hi in edges),
              f'(<100: {sum(d < 100 for d in ds)})')


if __name__ == '__main__':
    main(*sys.argv[1:])
