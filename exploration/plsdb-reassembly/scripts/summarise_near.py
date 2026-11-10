#!/usr/bin/env python3
"""Compare plasmid-surgery against near-relative references with the truth-reference pass.

For each reference set (surgery_near/<set>/, from REFSET=<set> scripts/run_surgery.sh) and each
(assembler, truth plasmid) with a stand-in reference in that set (near/near_refs.tsv): take the bin
that best matches the truth plasmid in the truth pass (surgery/by_plasmid.tsv, from
summarise_surgery.py), and compare its assessment against the truth plasmid and against the stand-in.

Writes near/by_plasmid.tsv and prints, per set and assembler: how often the stand-in picks the same
best bin, and the median drift (stand-in minus truth) in reference coverage, contig fraction aligned
and unique gap, and the carbapenemase gene's status.
"""
import csv
import os
import re
import statistics
import sys
from collections import defaultdict

SETS = ['nearest', 'd002', 'd005', 'd010', 'd020']


def family(gene):
    """blaKPC-3 -> blaKPC: the stand-in may carry another allele of the truth's gene."""
    return re.sub(r'-\d+$', '', gene)


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main(near_dir='near', surgery_dir='surgery_near', refsets_dir='refsets'):
    """Defaults: the PLSDB near-relative pass. Other passes: their own directories, e.g.
    summarise_near.py near-neighbor-with-carbapenemase{,/surgery,/refsets}"""
    truth_rows = {(r['assembler'], r['ref_id']): r for r in read_tsv('surgery/by_plasmid.tsv')}
    near = read_tsv(f'{near_dir}/near_refs.tsv')
    ref_genes = {}   # (set, bs) -> {ref acc: set of gene symbols}
    for r in near:
        if r['set'] == 'none':
            continue
        key = (r['set'], r['biosample'])
        if key not in ref_genes:
            genes = defaultdict(set)
            bed = f"{refsets_dir}/{r['set']}/{r['biosample']}.bed"
            if os.path.exists(bed):
                with open(bed) as f:
                    for line in f:
                        acc, _, _, name, *_ = line.rstrip('\n').split('\t')
                        genes[acc].add(name)
            ref_genes[key] = genes

    rows = []
    for n in near:
        s = n['set']
        if s == 'none' or not os.path.isdir(f'{surgery_dir}/{s}'):
            continue
        for a in ('shovill', 'unicycler', 'dragonflye'):
            t = truth_rows.get((a, n['truth_acc']))
            if not t or t.get('best_bin', 'none') == 'none':
                continue
            b = t['best_bin']
            path = f"{surgery_dir}/{s}/{a}/{n['biosample']}/{b}/assessment.tsv"
            if not os.path.exists(path):
                continue
            arows = read_tsv(path)
            stand = next((r for r in arows if r['ref_id'] == n['ref_acc']), None)
            if stand is None:
                continue
            # Which bin would the stand-in pick? (bin covering most of it, as for the truth)
            best_near_bin, best_cov = None, -1
            for p in sorted(os.listdir(f"{surgery_dir}/{s}/{a}/{n['biosample']}")):
                ap = f"{surgery_dir}/{s}/{a}/{n['biosample']}/{p}/assessment.tsv"
                if not os.path.exists(ap):
                    continue
                for r in read_tsv(ap):
                    if r['ref_id'] == n['ref_acc'] and float(r['ref_fraction_covered']) > best_cov:
                        best_near_bin, best_cov = p, float(r['ref_fraction_covered'])
            carbs = [c for c in n['truth_carbapenemases'].split(',') if c]
            on_ref = ref_genes.get((s, n['biosample']), {}).get(n['ref_acc'], set())
            carb_near = ''
            if carbs:
                feats = f"{surgery_dir}/{s}/{a}/{n['biosample']}/{b}/{n['ref_acc']}/ref_features.tsv"
                st = []
                if os.path.exists(feats):
                    fams = {family(c) for c in carbs}
                    st = [f"{r['name']}:{r['reconstruction_status']}/{r['found_in']}"
                          for r in read_tsv(feats) if family(r['name']) in fams]
                carb_near = ';'.join(st) or ('not_on_ref' if not ({family(c) for c in carbs} & {family(g) for g in on_ref})
                                             else 'not_annotated')
            rows.append(dict(
                set=s, assembler=a, biosample=n['biosample'], truth_acc=n['truth_acc'],
                ref_acc=n['ref_acc'], mash_dist=n['mash_dist'],
                truth_length=n['truth_length'], ref_length=n['ref_length'], bin=b,
                truth_cov=float(t['bin_fraction']), near_cov=float(stand['ref_fraction_covered']),
                truth_purity=float(t['bin_purity']), near_purity=float(stand['contig_fraction_aligned']),
                truth_gap=int(t['bin_unique_gap']), near_gap=int(stand['unique_gap_length']),
                truth_asm_cov=float(t['assembly_fraction']),
                near_asm_cov=float(stand['assembly_ref_fraction_covered']),
                same_best_bin=best_near_bin == b, near_best_bin=best_near_bin,
                carbapenemases=n['truth_carbapenemases'],
                carb_truth=t.get('carbapenemase_status', ''), carb_near=carb_near))

    fields = list(dict.fromkeys(k for r in rows for k in r))
    with open(f'{near_dir}/by_plasmid.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)

    def med(xs):
        return f'{statistics.median(xs):+.3f}' if xs else 'NA'
    print('set\tassembler\tn\tmedian_dist\tsame_best_bin\td_cov\td_purity\td_asm_cov\td_gap_bp\t'
          'carb_n\tcarb_not_on_ref\tcarb_status_same')
    for s in SETS:
        for a in ('shovill', 'unicycler', 'dragonflye'):
            rs = [r for r in rows if r['set'] == s and r['assembler'] == a]
            if not rs:
                continue
            carb = [r for r in rs if r['carbapenemases']]
            def status(x):   # 'present'/'partial'/'absent' per gene, ignoring found_in
                return sorted(p.split(':')[1].split('/')[0] for p in x.split(';')) if ':' in x else x
            print('\t'.join(map(str, [
                s, a, len(rs), f"{statistics.median(float(r['mash_dist']) for r in rs):.4f}",
                f"{sum(r['same_best_bin'] for r in rs)}/{len(rs)}",
                med([r['near_cov'] - r['truth_cov'] for r in rs]),
                med([r['near_purity'] - r['truth_purity'] for r in rs]),
                med([r['near_asm_cov'] - r['truth_asm_cov'] for r in rs]),
                statistics.median(r['near_gap'] - r['truth_gap'] for r in rs),
                len(carb), sum(r['carb_near'] == 'not_on_ref' for r in carb),
                f"{sum(status(r['carb_near']) == status(r['carb_truth']) for r in carb)}/{len(carb)}"])))


if __name__ == '__main__':
    main(*sys.argv[1:])
