"""
Reads as input (plasmid_surgery.reads). The end-to-end tests need minimap2 on PATH.
"""
import csv
import random
import shutil

import pytest

from array import array

from plasmid_surgery import reads
from plasmid_surgery.__main__ import main
from tests.simulate import simulate, write_fasta, revcomp, random_seq

needs_minimap2 = pytest.mark.skipif(shutil.which('minimap2') is None, reason="minimap2 not on PATH")

# simulate() deletes ref[5000:7000] from the sample
DELETED_START, DELETED_END = 5000, 7000


def write_fastq(path, seqs):
    with open(path, 'w') as f:
        for i, s in enumerate(seqs):
            f.write(f"@read{i}\n{s}\n+\n{'I' * len(s)}\n")


def tile_reads(seq, length, step, circular=True, rng=None, error_rate=0.0):
    """Reads tiling seq (both strands alternately), wrapping the origin if circular."""
    doubled = seq + seq if circular else seq
    out = []
    last = len(seq) if circular else len(seq) - length + 1
    for k, start in enumerate(range(0, last, step)):
        r = doubled[start:start + length]
        if error_rate and rng:
            r = ''.join(rng.choice('ACGT'.replace(b, '')) if rng.random() < error_rate else b for b in r)
        out.append(revcomp(r) if k % 2 else r)
    return out


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


# --- units ---

def test_alignment_blocks():
    # 10S 50M 2D 30M 5I 20M 100D 10M: the 2 bp deletion is covered, the 100 bp one isn't
    blocks, columns = reads.alignment_blocks(100, '10S50M2D30M5I20M100D10M')
    assert blocks == [(100, 202), (302, 312)]
    assert columns == 50 + 2 + 30 + 5 + 20 + 100 + 10


def test_run_length_encode():
    assert reads.run_length_encode(array('I', [0, 0, 3, 3, 1])) == [(0, 2, 0), (2, 4, 3), (4, 5, 1)]


def test_cache_round_trip(tmp_path):
    depth = array('I', [0, 0, 5, 5, 5, 2, 0])
    reads.write_cached_depth(tmp_path / 'x.depth.gz', depth)
    assert reads.read_cached_depth(tmp_path / 'x.depth.gz', 7) == depth
    assert reads.read_cached_depth(tmp_path / 'x.depth.gz', 8) is None      # wrong reference length
    assert reads.read_cached_depth(tmp_path / 'missing.depth.gz', 7) is None


def test_cache_key_depends_on_reads_settings_and_ref(tmp_path):
    fq = tmp_path / 'r.fq'
    fq.write_text('@r\nACGT\n+\nIIII\n')
    rs = reads.read_sets(short_reads=[fq])[0]
    key = reads.cache_key(rs, 'ACGT')
    assert reads.cache_key(rs, 'acgt') == key
    assert reads.cache_key(rs, 'ACGA') != key
    assert reads.cache_key({**rs, 'min_identity': 0.9}, 'ACGT') != key
    fq.write_text('@r\nACGTA\n+\nIIIII\n')
    assert reads.cache_key(rs, 'ACGT') != key       # a different file at the same path


@pytest.mark.parametrize("fraction,call", [(1.0, 'in_reads'), (0.9, 'in_reads'), (0.5, 'partly_in_reads'),
                                           (0.1, 'not_in_reads'), (0.0, 'not_in_reads')])
def test_gap_call(fraction, call):
    assert reads.gap_call(fraction) == call


def test_gap_support_and_metrics_wrap_origin():
    depth = array('I', [5] * 10 + [0] * 80 + [5] * 10)
    gaps = [{'start': 95, 'end': 105, 'type': 'unique'},        # wraps the origin: all in reads
            {'start': 20, 'end': 60, 'type': 'unique'},         # not in reads
            {'start': 60, 'end': 70, 'type': 'repeat'}]
    rows = reads.gap_support(gaps, depth, 'reconstruction', min_depth=3)
    assert [(r['gap'], r['ref_start'], r['ref_end'], r['length'], r['call']) for r in rows] == [
        ('gap_1_unique', 95, 5, 10, 'in_reads'),
        ('gap_2_unique', 20, 60, 40, 'not_in_reads'),
        ('gap_3_repeat', 60, 70, 10, 'not_in_reads'),
    ]
    m = reads.reads_metrics(depth, gaps, None, min_depth=3)
    assert m['reads_ref_fraction_covered'] == 0.2
    assert m['unique_gap_in_reads_length'] == 10     # unique gaps only
    assert m['assembly_unique_gap_in_reads_length'] == 'NA'


# --- end to end ---

def run(tmp_path, sim, contigs, extra_args=(), name='out'):
    recon = tmp_path / 'recon.fa'
    ref = tmp_path / 'ref.fa'
    write_fasta(recon, contigs)
    write_fasta(ref, [('ref', sim['ref'])])
    out = tmp_path / name
    rc = main(['--plasmid-reconstruction', str(recon), '--ref-plasmid', str(ref), '-o', str(out),
               '--log-level', 'warning', *map(str, extra_args)])
    assert rc == 0
    return out


@pytest.fixture(scope='module')
def sim():
    return simulate()


@pytest.fixture(scope='module')
def short_reads(sim, tmp_path_factory):
    path = tmp_path_factory.mktemp('reads') / 'short.fq'
    write_fastq(path, tile_reads(sim['sample'], 150, 25))
    return path


def drop_longest(sim):
    """The reconstruction without its longest contig: a unique gap the reads do cover."""
    longest = max(sim['contigs'], key=lambda c: len(c[1]))
    return [c for c in sim['contigs'] if c[0] != longest[0]], longest


@needs_minimap2
def test_gaps_in_reads_vs_not(tmp_path, sim, short_reads):
    contigs, dropped = drop_longest(sim)
    out = run(tmp_path, sim, contigs, ['--short-reads', short_reads])
    m = read_tsv(out / 'assessment.tsv')[0]
    ref_length = len(sim['ref'])
    # Everything but the deleted 2 kb is in the reads
    assert float(m['reads_ref_fraction_covered']) == pytest.approx(1 - (DELETED_END - DELETED_START) / ref_length, abs=0.005)
    assert float(m['reads_median_depth']) >= 5

    rows = read_tsv(out / 'ref' / 'gap_support.tsv')
    assert {r['level'] for r in rows} == {'reconstruction'}
    # The sample's deletion: not in the reads
    deletion = [r for r in rows if int(r['ref_start']) <= DELETED_START + 50 and int(r['ref_end']) >= DELETED_END - 50]
    assert deletion and all(r['call'] == 'not_in_reads' for r in deletion)
    # The dropped contig's region: in the reads, and it's the bulk of unique_gap_in_reads_length
    # (it holds the sample's 3 kb novel insertion, which isn't on the reference)
    dropped_on_ref = len(dropped[1]) - 3000
    in_reads = [r for r in rows if r['call'] == 'in_reads' and r['type'] == 'unique']
    assert sum(int(r['length']) for r in in_reads) >= 0.9 * dropped_on_ref
    assert int(m['unique_gap_in_reads_length']) >= 0.9 * dropped_on_ref
    assert int(m['unique_gap_in_reads_length']) < int(m['unique_gap_length'])

    bedgraph = reads.read_bedgraph(out / 'ref' / 'reads_depth.bedgraph.gz')
    assert bedgraph[0][0] == 0
    assert bedgraph[-1][1] == ref_length


@needs_minimap2
def test_assembly_level_gaps(tmp_path, sim, short_reads):
    contigs, dropped = drop_longest(sim)
    assembly = tmp_path / 'assembly.fa'
    write_fasta(assembly, sim['contigs'])     # the assembly has the dropped contig
    out = run(tmp_path, sim, contigs, ['--short-reads', short_reads, '--assembly', assembly])
    m = read_tsv(out / 'assessment.tsv')[0]
    rows = read_tsv(out / 'ref' / 'gap_support.tsv')
    assembly_rows = [r for r in rows if r['level'] == 'assembly']
    assert assembly_rows and all(r['gap'].startswith('assembly_gap_') for r in assembly_rows)
    # The assembly lacks only what the sample lacks (and repeats), so little of its unique gap is in the reads
    assert int(m['assembly_unique_gap_in_reads_length']) < 0.1 * int(m['unique_gap_in_reads_length'])


@needs_minimap2
def test_long_reads_cover_origin(tmp_path, sim):
    rng = random.Random(3)
    fq = tmp_path / 'long.fq'
    # ~20x: near a structural difference (the deletion) reads align on one side only, so depth dips there
    write_fastq(fq, tile_reads(sim['sample'], 6000, 300, rng=rng, error_rate=0.04))
    out = run(tmp_path, sim, sim['contigs'], ['--long-reads', fq])
    m = read_tsv(out / 'assessment.tsv')[0]
    ref_length = len(sim['ref'])
    assert float(m['reads_ref_fraction_covered']) == pytest.approx(1 - (DELETED_END - DELETED_START) / ref_length, abs=0.01)
    depth = reads.read_bedgraph(out / 'ref' / 'reads_depth.bedgraph.gz')
    assert depth[0][2] >= 3 and depth[-1][2] >= 3      # both sides of the origin


@needs_minimap2
def test_divergent_reads_not_counted(tmp_path, sim):
    rng = random.Random(5)
    fq = tmp_path / 'other.fq'
    # ~7% divergent from the sample: below the 95% short-read identity threshold, but still aligned
    write_fastq(fq, tile_reads(sim['sample'], 150, 25, rng=rng, error_rate=0.07))
    out = run(tmp_path, sim, sim['contigs'], ['--short-reads', fq])
    m = read_tsv(out / 'assessment.tsv')[0]
    assert float(m['reads_ref_fraction_covered']) < 0.2     # a few reads are >= 95% identical by chance
    out = run(tmp_path, sim, sim['contigs'], ['--short-reads', fq, '--min-read-identity', '0.85'], name='out_loose')
    assert float(read_tsv(out / 'assessment.tsv')[0]['reads_ref_fraction_covered']) > 0.8


@needs_minimap2
def test_reads_cache_reused(tmp_path, sim, short_reads, monkeypatch):
    cache = tmp_path / 'cache'
    out1 = run(tmp_path, sim, sim['contigs'], ['--short-reads', short_reads, '--reads-cache', cache], name='out1')
    assert len(list(cache.glob('*.depth.gz'))) == 1

    assert len(list(cache.glob('*.contigs.json'))) == 1

    def fail(*args, **kwargs):
        raise AssertionError("reads aligned again despite the cache")
    monkeypatch.setattr(reads, 'align_read_set', fail)
    monkeypatch.setattr(reads, 'align_reads_to_contigs', fail)
    out2 = run(tmp_path, sim, sim['contigs'], ['--short-reads', short_reads, '--reads-cache', cache], name='out2')
    assert (reads.read_bedgraph(out1 / 'ref' / 'reads_depth.bedgraph.gz')
            == reads.read_bedgraph(out2 / 'ref' / 'reads_depth.bedgraph.gz'))
    assert (out1 / 'contig_depth.tsv').read_text() == (out2 / 'contig_depth.tsv').read_text()


@needs_minimap2
def test_contig_depth_from_assembly(tmp_path, sim, short_reads):
    contigs, dropped = drop_longest(sim)
    assembly = tmp_path / 'assembly.fa'
    write_fasta(assembly, sim['contigs'])
    out = run(tmp_path, sim, contigs, ['--short-reads', short_reads, '--assembly', assembly])
    rows = {r['contig_id']: r for r in read_tsv(out / 'contig_depth.tsv')}
    assert set(rows) == {cid for cid, _ in sim['contigs']}          # every assembly contig
    assert rows[dropped[0]]['in_reconstruction'] == 'False'
    assert all(rows[cid]['in_reconstruction'] == 'True' for cid, _ in contigs)
    # Reads tile the sample at 150/25 = 6x; the collapsed repeat contig gets all 3 copies' reads
    repeat_id = next(cid for cid, seq in sim['contigs'] if seq == sim['repeat'])
    unique = [float(r['reads_depth']) for cid, r in rows.items() if cid != repeat_id]
    assert all(5 <= d <= 7 for d in unique), unique
    assert 15 <= float(rows[repeat_id]['reads_depth']) <= 21
    m = read_tsv(out / 'assessment.tsv')[0]
    assert 5 <= float(m['reads_reconstruction_depth']) <= 9
    sources = {r['contig_id']: r for r in read_tsv(out / 'ref' / 'contig_sources.tsv')}
    assert float(sources[dropped[0]]['reads_depth']) == float(rows[dropped[0]]['reads_depth'])


@needs_minimap2
def test_contig_depth_without_assembly(tmp_path, sim, short_reads):
    out = run(tmp_path, sim, sim['contigs'], ['--short-reads', short_reads])
    rows = read_tsv(out / 'contig_depth.tsv')
    assert {r['contig_id'] for r in rows} == {cid for cid, _ in sim['contigs']}
    assert all(r['in_reconstruction'] == 'True' and r['bin'] == 'NA' for r in rows)


@needs_minimap2
def test_ref_features_reads_fraction(tmp_path, sim, short_reads):
    bed = tmp_path / 'feats.bed'
    bed.write_text(f"ref\t{DELETED_START + 100}\t{DELETED_START + 900}\tgone\t0\t+\n"
                   f"ref\t20000\t21000\tthere\t0\t+\n")
    out = run(tmp_path, sim, sim['contigs'], ['--short-reads', short_reads, '--ref-features', bed])
    rows = {r['name']: r for r in read_tsv(out / 'ref' / 'ref_features.tsv')}
    assert float(rows['gone']['reads_fraction']) == 0.0
    assert float(rows['there']['reads_fraction']) == 1.0


@needs_minimap2
def test_without_reads_fields_are_na(tmp_path, sim):
    out = run(tmp_path, sim, sim['contigs'])
    m = read_tsv(out / 'assessment.tsv')[0]
    assert all(m[f] == 'NA' for f in reads.READS_METRICS_FIELDNAMES)
    assert not (out / 'ref' / 'gap_support.tsv').exists()
    assert not (out / 'contig_depth.tsv').exists()


def test_short_reads_takes_at_most_two_files(tmp_path):
    with pytest.raises(SystemExit):
        main(['--plasmid-reconstruction', 'x.fa', '--ref-plasmid', 'r.fa', '--short-reads', 'a', 'b', 'c'])


@needs_minimap2
def test_report_shows_reads(tmp_path, sim, short_reads):
    contigs, dropped = drop_longest(sim)
    assembly = tmp_path / 'assembly.fa'
    write_fasta(assembly, sim['contigs'])
    out = run(tmp_path, sim, contigs, ['--short-reads', short_reads, '--assembly', assembly, '--report'])
    page = (out / 'report.html').read_text()
    assert 'read depth' in page and 'In the reads' in page and '<th>Read depth</th>' in page
    assert 'not in the reads' in page and 'in the reads (100%' in page      # gap tooltips
    assert 'below the minimum read depth' in page
