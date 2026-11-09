#!/usr/bin/env python3
"""Summarise plasmid-surgery runs (scripts/run_surgery.sh) over the PLSDB reassembly set.

Writes:
  surgery/by_plasmid.tsv  one row per (assembler, truth plasmid): best bin, its coverage and purity,
                          how many bins hold a piece, assembly-level coverage, binning vs assembly loss,
                          carbapenemase status
  surgery/by_bin.tsv      one row per (assembler, mob_recon bin): best truth plasmid, coverage, purity
and prints aggregates by assembler.
"""
import csv
import glob
import os
import statistics
import sys
from collections import defaultdict

PIECE = 0.05      # a bin holds a "piece" of a plasmid if it covers >= 5% of it
MATCH = 0.05      # a bin "matches no truth plasmid" if its best coverage is below this


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main(root='.'):
    os.chdir(root)
    truth = {}
    with open('truth/truth.tsv.part') as f:
        for row in csv.reader(f, delimiter='\t'):
            bs, acc, length, topo, compl, carb, is_carb = (row + [''] * 7)[:7]
            truth[(bs, acc)] = dict(length=int(length), carbapenemases=carb, carb_plasmid=is_carb)

    # assessments[(a, bs, acc)] = list of (bin, row)
    assessments = defaultdict(list)
    bins = []
    for path in sorted(glob.glob('surgery/*/*/*/assessment.tsv')):
        _, a, bs, b, _ = path.split('/')
        rows = read_tsv(path)
        bins.append((a, bs, b, rows))
        for r in rows:
            assessments[(a, bs, r['ref_id'])].append((b, r))

    assemblers = sorted({a for a, _, _ in assessments})
    by_plasmid = []
    for (bs, acc), t in sorted(truth.items()):
        for a in assemblers:
            runs = assessments.get((a, bs, acc))
            if not runs:
                if os.path.isdir(f'assemblies/{a}/{bs}'):
                    # no bins at all for this assembly: plasmid entirely unbinned
                    by_plasmid.append(dict(assembler=a, biosample=bs, ref_id=acc, ref_length=t['length'],
                                           carbapenemases=t['carbapenemases'], best_bin='none'))
                continue
            b, best = max(runs, key=lambda x: (float(x[1]['ref_fraction_covered']), float(x[1]['score'])))
            pieces = [x for x in runs if float(x[1]['ref_fraction_covered']) >= PIECE]
            bin_frac = float(best['ref_fraction_covered'])
            asm_frac = float(best['assembly_ref_fraction_covered'])
            carb_status = ''
            if t['carbapenemases']:
                feats = f'surgery/{a}/{bs}/{b}/{acc}/ref_features.tsv'
                if os.path.exists(feats):
                    wanted = set(t['carbapenemases'].split(','))
                    st = [f"{r['name']}:{r['reconstruction_status']}/{r['found_in']}"
                          for r in read_tsv(feats) if r['name'] in wanted]
                    carb_status = ';'.join(st) or 'not_annotated'
            by_plasmid.append(dict(
                assembler=a, biosample=bs, ref_id=acc, ref_length=t['length'],
                carbapenemases=t['carbapenemases'], best_bin=b,
                bin_fraction=bin_frac, bin_unique_gap=best['unique_gap_length'],
                bin_purity=best['contig_fraction_aligned'], bin_contigs=best['num_contigs'],
                n_bins_with_piece=len(pieces),
                assembly_fraction=asm_frac, assembly_unique_gap=best['assembly_unique_gap_length'],
                binning_loss=round(max(asm_frac - bin_frac, 0), 4), assembly_loss=round(1 - asm_frac, 4),
                missed_contigs=best['num_missed_contigs'], carbapenemase_status=carb_status))

    by_bin = []
    for a, bs, b, rows in bins:
        best = max(rows, key=lambda r: float(r['score']))
        cov = float(best['ref_fraction_covered'])
        by_bin.append(dict(assembler=a, biosample=bs, bin=b, bin_length=best['total_contig_length'],
                           bin_contigs=best['num_contigs'], best_ref=best['ref_id'] if cov >= MATCH else 'none',
                           ref_fraction=cov, purity=best['contig_fraction_aligned']))

    for name, rows in (('surgery/by_plasmid.tsv', by_plasmid), ('surgery/by_bin.tsv', by_bin)):
        fields = list(dict.fromkeys(k for r in rows for k in r))
        with open(name, 'w', newline='') as f:
            w = csv.DictWriter(f, fields, delimiter='\t', restval='NA', lineterminator='\n')
            w.writeheader(); w.writerows(rows)

    # Aggregates
    def med(xs):
        return round(statistics.median(xs), 3) if xs else 'NA'
    print('assembler\tplasmids\tno_bin\tbin>=0.95\tbin>=0.99\tmedian_bin\tmedian_asm\tsplit>1bin\t'
          'mean_binning_loss\tmean_assembly_loss\tcarb_present\tcarb_total')
    for a in assemblers:
        rs = [r for r in by_plasmid if r['assembler'] == a]
        ok = [r for r in rs if r['best_bin'] != 'none']
        carb = [r for r in rs if r['carbapenemases']]
        carb_ok = [r for r in carb if r.get('carbapenemase_status') and
                   all(':present/' in s for s in r['carbapenemase_status'].split(';'))]
        print('\t'.join(map(str, [
            a, len(rs), len(rs) - len(ok),
            sum(r['bin_fraction'] >= 0.95 for r in ok), sum(r['bin_fraction'] >= 0.99 for r in ok),
            med([r['bin_fraction'] for r in ok]), med([r['assembly_fraction'] for r in ok]),
            sum(r['n_bins_with_piece'] > 1 for r in ok),
            round(statistics.mean(r['binning_loss'] for r in ok), 4),
            round(statistics.mean(r['assembly_loss'] for r in ok), 4),
            len(carb_ok), len(carb)])))
    print()
    print('assembler\tbins\tno_truth_match\tmedian_purity\tpurity<0.9')
    for a in assemblers:
        bs_ = [r for r in by_bin if r['assembler'] == a]
        print('\t'.join(map(str, [a, len(bs_), sum(r['best_ref'] == 'none' for r in bs_),
                                  med([float(r['purity']) for r in bs_ if r['best_ref'] != 'none']),
                                  sum(float(r['purity']) < 0.9 for r in bs_ if r['best_ref'] != 'none')])))


if __name__ == '__main__':
    main(*sys.argv[1:])
