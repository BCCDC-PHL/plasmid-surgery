#!/usr/bin/env python3
"""Build a carbapenemase-plasmid database from PLSDB: every plasmid with an AMRFinderPlus hit whose
antimicrobial_agent includes CARBAPENEM (amr.tsv).

Writes to OUTDIR:
  plasmids.fasta      their sequences (from sequences.fasta.bz2; 80 bp lines, so samtools faidx works)
  features.bed        all their AMRFinderPlus genes (name = gene symbol), for --ref-features
  plasmids.tsv        accession, length, carbapenemase gene(s), description
Then index and sketch it (same sketch settings as PLSDB's plsdb_sketch.msh: k=21, s=1000):
  samtools faidx OUTDIR/plasmids.fasta
  mash sketch -p 16 -k 21 -s 1000 -i -o OUTDIR/plasmids OUTDIR/plasmids.fasta
"""
import argparse
import bz2
import csv
import os
import sys
from collections import defaultdict

csv.field_size_limit(sys.maxsize)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--plsdb', required=True)
    p.add_argument('-o', '--outdir', required=True)
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    carb = defaultdict(set)
    beds = defaultdict(set)
    with open(os.path.join(a.plsdb, 'amr.tsv'), newline='') as f:
        rows = csv.reader(f, delimiter='\t')
        next(rows)
        for r in rows:
            if r[1] != 'amrfinderplus':
                continue
            beds[r[0]].add((r[0], int(r[8]) - 1, int(r[9]), r[3], '0', r[11]))
            if 'CARBAPENEM' in r[7]:
                carb[r[0]].add(r[3])

    desc = {}
    with open(os.path.join(a.plsdb, 'nuccore.csv'), newline='') as f:
        for r in csv.DictReader(f):
            if r['NUCCORE_ACC'] in carb:
                desc[r['NUCCORE_ACC']] = r['NUCCORE_Description']

    lengths = {}
    with bz2.open(os.path.join(a.plsdb, 'sequences.fasta.bz2'), 'rt') as f, \
            open(os.path.join(a.outdir, 'plasmids.fasta'), 'w') as out:
        name, chunks = None, []

        def flush():
            if name is not None and name in carb:
                seq = ''.join(chunks)
                lengths[name] = len(seq)
                out.write(f'>{name}\n')
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i + 80] + '\n')
        for line in f:
            if line.startswith('>'):
                flush()
                name, chunks = line[1:].split()[0], []
            elif name in carb:
                chunks.append(line.strip())
        flush()

    with open(os.path.join(a.outdir, 'features.bed'), 'w') as f:
        for acc in sorted(lengths):
            for b in sorted(beds[acc], key=lambda x: (x[1], x[2], x[3])):
                f.write('\t'.join(map(str, b)) + '\n')
    with open(os.path.join(a.outdir, 'plasmids.tsv'), 'w') as f:
        f.write('accession\tlength\tcarbapenemases\tdescription\n')
        for acc in sorted(lengths):
            f.write(f"{acc}\t{lengths[acc]}\t{','.join(sorted(carb[acc]))}\t{desc.get(acc, '')}\n")
    missing = carb.keys() - lengths.keys()
    print(f'{len(carb)} plasmids with a carbapenemase; {len(lengths)} written'
          + (f'; {len(missing)} not in sequences.fasta.bz2' if missing else ''))


if __name__ == '__main__':
    main()
