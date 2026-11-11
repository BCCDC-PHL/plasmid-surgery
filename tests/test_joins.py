"""
Joins (plasmid_surgery.joins): composition, where joins are, and long-read evidence at them.
"""
import csv
import random
import shutil

import pytest

from plasmid_surgery import joins
from plasmid_surgery.__main__ import main
from tests.simulate import write_fasta, random_seq, revcomp

needs_minimap2 = pytest.mark.skipif(shutil.which('minimap2') is None, reason="minimap2 not on PATH")


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def ivs(placements, alternatives=None):
    """contig_intervals-style input for one contig 'c'."""
    alternatives = alternatives or []
    return {'placements': {'c': placements}, 'all': {'c': placements + alternatives}}


# --- units ---

def test_composition_greedy_and_segments():
    contigs = [{'seq_id': 'c', 'seq': 'A' * 100000, 'description': ''}]
    per_ref = [('A1', ivs([(0, 40000)])), ('A2', ivs([(0, 38000)])), ('B', ivs([(40000, 90000)])),
               ('tiny', ivs([(95000, 95500)]))]          # shorter than min_block: ignored
    comp = joins.composition(contigs, per_ref, min_block=1000)
    assert [(r['ref_id'], r['new_bases']) for r in comp['rows']] == [('B', 50000), ('A1', 40000), ('none', 10000)]
    assert comp['rows'][1]['cumulative_fraction'] == 0.9
    assert comp['segments']['c'] == [(0, 40000, 'A1'), (40000, 90000, 'B'), (90000, 100000, None)]


def test_joins_cover_a_shared_repeat_and_the_circular_junction():
    contigs = [{'seq_id': 'c', 'seq': 'A' * 100000, 'description': 'len=100000 circular=true'}]
    # B's placement starts at 40,000; its IS copy at 38,800-40,000 is an alternative placement
    per_ref = [('A', ivs([(0, 40000)])), ('B', ivs([(40000, 100000)], [(38800, 40000)]))]
    comp = joins.composition(contigs, per_ref)
    found = joins.find_joins(contigs, comp)
    assert [(j['type'], j['left'], j['right'], j['region_start'], j['region_end']) for j in found] == [
        ('reference_switch', 'A', 'B', 38800, 40000),
        ('circular_junction', 'B', 'A', 100000, 100000),
    ]


def test_unexplained_boundaries():
    contigs = [{'seq_id': 'c', 'seq': 'A' * 50000, 'description': ''}]
    comp = joins.composition(contigs, [('A', ivs([(0, 20000)])), ('B', ivs([(30000, 50000)]))])
    found = joins.find_joins(contigs, comp)
    assert [(j['type'], j['left'], j['right'], j['region_start']) for j in found] == [
        ('unexplained_boundary', 'A', 'none', 20000), ('unexplained_boundary', 'none', 'B', 30000)]


def test_join_evidence_calls():
    contig = {'seq_id': 'c', 'seq': 'A' * 50000, 'description': ''}
    join = {'contig_id': 'c', 'region_start': 20000, 'region_end': 21000}
    depth = {'depth': [[10] * 50000]}

    def evidence(records):
        return {'long': {**depth, 'records': records}}
    spanning = [(0, 15000, 26000, 0, 0)] * 3
    clipped = [(0, 12000, 20990, 0, 3000)] * 4                 # ends at the join, the rest of the read elsewhere
    short_anchor = [(0, 19500, 26000, 0, 0)]                    # starts too close to the join to span it
    assert joins.join_evidence(join, contig, 0, evidence(spanning))['call'] == 'supported'
    row = joins.join_evidence(join, contig, 0, evidence(clipped + short_anchor))
    assert (row['long_reads_spanning'], row['long_reads_clipped'], row['call']) == (0, 4, 'contradicted')
    assert joins.join_evidence(join, contig, 0, evidence(spanning[:1] + clipped))['call'] == 'weak'
    assert joins.join_evidence(join, contig, 0, {})['call'] == 'no_long_reads'


def test_circular_contigs_extended():
    c = {'seq_id': 'c', 'seq': 'ACGT' * 10, 'description': 'circular=Y'}
    assert joins.extended(c) == c['seq'] * 2
    assert joins.extended({**c, 'description': ''}) == c['seq']
    assert joins.is_circular({'description': '[topology=circular]'})


# --- end to end: a cointegrate vs. two plasmids misassembled into one ---

def tile_circular(seq, length, step, rng, error_rate=0.03):
    doubled = seq + seq
    out = []
    for k, start in enumerate(range(0, len(seq), step)):
        r = doubled[start:start + length]
        r = ''.join(rng.choice('ACGT'.replace(b, '')) if rng.random() < error_rate else b for b in r)
        out.append(revcomp(r) if k % 2 else r)
    return out


def write_fastq(path, seqs):
    with open(path, 'w') as f:
        for i, s in enumerate(seqs):
            f.write(f"@read{i}\n{s}\n+\n{'I' * len(s)}\n")


@pytest.fixture(scope='module')
def plasmids():
    rng = random.Random(11)
    insertion = random_seq(rng, 1200)        # an IS on both plasmids
    a = random_seq(rng, 38800) + insertion
    b = random_seq(rng, 48800) + insertion
    # The cointegrate (one IS copy at each join), as one circular contig
    joined = a + b
    return {'a': a, 'b': b, 'joined': joined, 'rng': rng}


def run_joined(tmp_path, plasmids, reads_from):
    write_fasta(tmp_path / 'recon.fa', [('contig_1 len=100000 circular=true', plasmids['joined'])])
    write_fasta(tmp_path / 'refs.fa', [('pA', plasmids['a']), ('pB', plasmids['b'])])
    reads = []
    for seq in reads_from:
        reads += tile_circular(seq, 8000, 400, plasmids['rng'])
    write_fastq(tmp_path / 'long.fq', reads)
    out = tmp_path / 'out'
    assert main(['--plasmid-reconstruction', str(tmp_path / 'recon.fa'), '--ref-plasmid', str(tmp_path / 'refs.fa'),
                 '--long-reads', str(tmp_path / 'long.fq'), '-t', '2', '-o', str(out), '--report',
                 '--log-level', 'error']) == 0
    return out


@needs_minimap2
def test_cointegrate_joins_supported(tmp_path, plasmids):
    out = run_joined(tmp_path, plasmids, [plasmids['joined']])
    comp = read_tsv(out / 'composition.tsv')
    assert {r['ref_id'] for r in comp[:2]} == {'pA', 'pB'} and float(comp[1]['cumulative_fraction']) > 0.99
    rows = read_tsv(out / 'joins.tsv')
    assert {r['type'] for r in rows} == {'reference_switch', 'circular_junction'}
    switch = [r for r in rows if r['type'] == 'reference_switch']
    assert switch and all(int(r['region_length']) >= 1100 for r in switch)    # the shared IS
    assert all(r['call'] == 'supported' for r in rows), rows
    assert (out / 'reconstruction_long_reads_depth.bedgraph.gz').exists()


@needs_minimap2
def test_misjoined_plasmids_contradicted(tmp_path, plasmids):
    # The sample has the two plasmids separately; the contig joins them through the shared IS
    out = run_joined(tmp_path, plasmids, [plasmids['a'], plasmids['b']])
    rows = read_tsv(out / 'joins.tsv')
    assert rows and all(r['call'] == 'contradicted' for r in rows), rows
    assert all(int(r['long_reads_spanning']) == 0 and int(r['long_reads_clipped']) >= 3 for r in rows)
    html = (out / 'report.html').read_text()
    assert 'contradicted' in html


@needs_minimap2
def test_joins_without_reads(tmp_path, plasmids):
    write_fasta(tmp_path / 'recon.fa', [('contig_1 circular=true', plasmids['joined'])])
    write_fasta(tmp_path / 'refs.fa', [('pA', plasmids['a']), ('pB', plasmids['b'])])
    assert main(['--plasmid-reconstruction', str(tmp_path / 'recon.fa'), '--ref-plasmid', str(tmp_path / 'refs.fa'),
                 '-o', str(tmp_path / 'out'), '--log-level', 'error']) == 0
    rows = read_tsv(tmp_path / 'out' / 'joins.tsv')
    assert rows and all(r['call'] == 'no_long_reads' and r['long_reads_spanning'] == 'NA' for r in rows)


def test_widen_to_repeats():
    contigs = [{'seq_id': 'c', 'seq': 'A' * 100000, 'description': 'circular=true'}]
    found = [{'contig_id': 'c', 'type': 'reference_switch', 'region_start': 40000, 'region_end': 40000},
             {'contig_id': 'c', 'type': 'circular_junction', 'region_start': 100000, 'region_end': 100000}]
    # An IS copy ending at the join, and copies ending and starting the contig (around the circular junction)
    repeats = {'c': [(38800, 40000), (99000, 100000), (0, 300)]}
    joins.widen_to_repeats(found, contigs, repeats)
    assert (found[0]['region_start'], found[0]['region_end']) == (38800, 40000)
    assert (found[1]['region_start'], found[1]['region_end']) == (99000, 100300)    # past the end: wraps
    rows = joins.assess_joins(contigs, found)
    assert (rows[1]['region_start'], rows[1]['region_end'], rows[1]['region_length']) == (99000, 300, 1300)
