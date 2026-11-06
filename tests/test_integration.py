"""
End-to-end tests on simulated data. Require minimap2 on PATH.
"""
import csv
import random
import re
import shutil

import pytest

from plasmid_surgery.__main__ import main
from plasmid_surgery.parsers import parse_fasta
from tests.simulate import simulate, write_fasta, revcomp, random_seq

pytestmark = pytest.mark.skipif(shutil.which('minimap2') is None, reason="minimap2 not on PATH")


def read_summary(path):
    lines = path.read_text().splitlines()
    header = lines[0].split('\t')
    return [dict(zip(header, line.split('\t'))) for line in lines[1:]]


def run(tmp_path, sim, contigs=None, prefix=None, extra_args=()):
    recon = tmp_path / "recon.fa"
    ref = tmp_path / "ref.fa"
    write_fasta(recon, contigs if contigs is not None else sim['contigs'])
    write_fasta(ref, [('ref', sim['ref'])])
    out_dir = tmp_path / "out"
    rc = main([
        '--plasmid-reconstruction', str(recon),
        '--ref-plasmid', str(ref),
        '--tmpdir', str(tmp_path / "tmp"),
        '-o', str(out_dir),
        '--log-level', 'warning',
        *(['-p', prefix] if prefix else []),
        *extra_args,
    ])
    assert rc == 0
    name = (lambda f: f"{prefix}_{f}") if prefix else (lambda f: f)
    output = parse_fasta(out_dir / name("scaffold.fa"))
    assert len(output) == 1
    return output[0], read_summary(out_dir / name("summary.tsv"))


def contained_circularly(segment, seq):
    doubled = seq + seq
    return segment in doubled or revcomp(segment) in doubled


@pytest.fixture
def sim():
    return simulate(seed=1)


def test_every_unique_contig_recovered_intact(tmp_path, sim):
    output, _ = run(tmp_path, sim)
    for contig_id, segment in sim['truth_segments'][:-1]:  # last one is the repeat
        assert contained_circularly(segment, output['seq']), f"{contig_id} not intact in output"


def test_scaffold_reconstructs_sample_exactly(tmp_path, sim):
    # The two uncovered repeat copies are filled (in lowercase) from the repeat
    # contig's alternative placements, so the scaffold is the complete sample
    # plasmid, rotated to the ref origin.
    output, _ = run(tmp_path, sim)
    seq = output['seq']
    assert seq.count('N') == 0
    assert sum(c.islower() for c in seq) == 2 * len(sim['repeat'])
    assert len(seq) == len(sim['sample'])
    assert seq.upper() in sim['sample'] + sim['sample']


def test_inferred_lowercase_regions_are_repeat_copies(tmp_path, sim):
    output, _ = run(tmp_path, sim)
    runs = [m.group() for m in re.finditer('[acgtn]+', output['seq'])]
    assert len(runs) == 2
    assert all(r.upper() == sim['repeat'] for r in runs)


def test_summary(tmp_path, sim):
    _, rows = run(tmp_path, sim)
    by_id = {}
    for row in rows:
        by_id.setdefault(row['contig_id'], []).append(row)
    assert set(by_id) == {c for c, _ in sim['contigs']}
    # Odd-numbered contigs were reverse-complemented
    for contig_id, contig_rows in by_id.items():
        expected_rc = str(int(contig_id.split('_')[1]) % 2 == 1)
        assert all(r['reverse_complement'] == expected_rc for r in contig_rows)
    # The origin-spanning contig is placed as a single wrapping alignment
    # (rather than minimap2's primary + supplementary split on a linear ref)
    origin_rows = [r for r in rows if r['spans_origin'] == 'True']
    assert len(origin_rows) == 1
    origin_row = origin_rows[0]
    assert len(by_id[origin_row['contig_id']]) == 1
    assert origin_row['is_supplementary'] == 'False'
    assert int(origin_row['alignment_end']) < int(origin_row['alignment_start'])
    # The repeat contig maps ambiguously
    repeat_id = sim['truth_segments'][-1][0]
    assert by_id[repeat_id][0]['mapq'] == '0'
    assert all(r['included_in_output'] == 'True' for r in rows)
    # Rows sorted by alignment start
    starts = [int(r['alignment_start']) for r in rows]
    assert starts == sorted(starts)


def origin_spanning_contig_id(rows):
    spanning = {r['contig_id'] for r in rows if r['spans_origin'] == 'True'}
    assert len(spanning) == 1
    return spanning.pop()


def test_overlapping_contigs_not_duplicated(tmp_path):
    # Contigs overlap their neighbours by 127 bp, like a SPAdes k=128 assembly
    sim = simulate(seed=2, contig_overlap=127)
    output, rows = run(tmp_path, sim)
    for contig_id, segment in sim['truth_segments'][:-1]:
        # Each contig minus its trailing overlap should be present
        assert contained_circularly(segment[:-127], output['seq']), f"{contig_id} not intact in output"
    # No overlap sequence is duplicated and nothing is lost: with the repeat
    # copies filled from alternative placements, the scaffold is exactly the sample.
    seq = output['seq']
    assert seq.count('N') == 0
    assert seq.upper() in sim['sample'] + sim['sample']
    assert len(seq) == len(sim['sample'])
    # The contigs preceding two of the repeat copies extend 127 bp into them,
    # so less than the full two copies is inferred
    assert sum(c.islower() for c in seq) == 2 * len(sim['repeat']) - 2 * 127


def test_origin_spanning_contig_with_split_gap(tmp_path):
    # Regression: on a linear ref, minimap2 split this origin-spanning contig into
    # primary + supplementary alignments with 7 bp of query aligned in neither,
    # so those bases were lost from the output.
    sim = simulate(seed=2, contig_overlap=127)
    output, rows = run(tmp_path, sim)
    origin_id = origin_spanning_contig_id(rows)
    segment = dict(sim['truth_segments'])[origin_id]
    assert contained_circularly(segment[:-127], output['seq'])


EXPECTED_OUTPUTS = ['assessment.tsv', 'gaps.bed', 'events.tsv', 'unaligned.fa', 'summary.tsv', 'scaffold.fa']


def test_default_outdir_is_cwd_and_tmpdir_cleaned_up(tmp_path, sim, monkeypatch):
    monkeypatch.chdir(tmp_path)
    write_fasta('recon.fa', sim['contigs'])
    write_fasta('ref.fa', [('ref', sim['ref'])])
    rc = main(['--plasmid-reconstruction', 'recon.fa', '--ref-plasmid', 'ref.fa',
               '--output-bases-per-line', '60', '--log-level', 'warning'])
    assert rc == 0
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(EXPECTED_OUTPUTS + ['recon.fa', 'ref.fa'])
    lines = (tmp_path / 'scaffold.fa').read_text().splitlines()
    assert lines[0] == '>recon'
    assert max(len(l) for l in lines[1:]) == 60


def test_outdir_and_prefix(tmp_path, sim):
    recon = tmp_path / "recon.fa"
    ref = tmp_path / "ref.fa"
    write_fasta(recon, sim['contigs'])
    write_fasta(ref, [('ref', sim['ref'])])
    out_dir = tmp_path / "nested" / "out"
    rc = main(['--plasmid-reconstruction', str(recon), '--ref-plasmid', str(ref),
               '-o', str(out_dir), '-p', 'sample1', '--log-level', 'warning'])
    assert rc == 0
    assert sorted(p.name for p in out_dir.iterdir()) == sorted(f"sample1_{name}" for name in EXPECTED_OUTPUTS)
    # Prefix is the default sample ID
    assert (out_dir / 'sample1_scaffold.fa').read_text().startswith('>sample1\n')
    [m] = read_tsv(out_dir / 'sample1_assessment.tsv')
    assert m['sample_id'] == 'sample1'


def test_sample_id_overrides_prefix(tmp_path, sim):
    output, _ = run(tmp_path, sim, prefix='run42', extra_args=['--sample-id', 'sample1'])
    assert output['seq_id'] == 'sample1'
    out_dir = tmp_path / "out"
    [m] = read_tsv(out_dir / 'run42_assessment.tsv')
    assert m['sample_id'] == 'sample1'


def test_stale_ref_in_tmpdir_not_reused(tmp_path, sim):
    tmpdir = tmp_path / "tmp"
    tmpdir.mkdir()
    # A different sequence with the same ref ID, left over from a previous run
    write_fasta(tmpdir / "ref.fa", [('ref', 'ACGT' * 1000)])
    output, rows = run(tmp_path, sim)
    assert all(r['is_mapped'] == 'True' for r in rows)


# --- assessment outputs ---

def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def test_assessment_metrics(tmp_path, sim):
    run(tmp_path, sim)
    [m] = read_tsv(tmp_path / "out" / "assessment.tsv")
    ref_length = len(sim['ref'])
    repeat_length = len(sim['repeat'])
    assert m['sample_id'] == 'recon'
    assert int(m['ref_length']) == ref_length
    # The collapsed repeat contig covers one of three repeat copies
    assert int(m['num_gaps']) == 2
    assert int(m['total_gap_length']) == 2 * repeat_length
    assert int(m['ref_bases_covered']) == ref_length - 2 * repeat_length
    assert int(m['ref_bases_covered_multiple']) == 0
    assert m['contig_fraction_aligned'] == '1.0'
    assert int(m['num_unaligned_contigs']) == 0
    assert int(m['num_large_insertions']) == 1
    assert int(m['num_large_deletions']) == 1
    assert int(m['num_ambiguous_alignments']) == 1


def test_gaps_bed_are_repeat_copies(tmp_path, sim):
    run(tmp_path, sim)
    lines = (tmp_path / "out" / "gaps.bed").read_text().splitlines()
    assert len(lines) == 2
    for line in lines:
        _, start, end, _ = line.split('\t')
        assert sim['ref'][int(start):int(end)] == sim['repeat']


def test_events_match_simulated_indels(tmp_path, sim):
    # simulate() deletes ref[5000:7000] and inserts 3 kb of novel sequence at ref 40000
    run(tmp_path, sim)
    events = {e['event_type']: e for e in read_tsv(tmp_path / "out" / "events.tsv")}
    assert set(events) == {'insertion', 'deletion'}
    deletion = events['deletion']
    assert abs(int(deletion['ref_start']) - 5000) <= 10
    assert abs(int(deletion['length']) - 2000) <= 10
    insertion = events['insertion']
    assert abs(int(insertion['ref_start']) - 40000) <= 10
    assert abs(int(insertion['length']) - 3000) <= 10


def test_unaligned_sequence_reported(tmp_path, sim):
    rng = random.Random(99)
    foreign_tail = random_seq(rng, 2000)
    contaminant = random_seq(rng, 5000)
    contigs = list(sim['contigs'])
    tail_contig_id, tail_contig_seq = contigs[2]
    contigs[2] = (tail_contig_id, tail_contig_seq + foreign_tail)
    contigs.append(('contaminant', contaminant))
    _, rows = run(tmp_path, sim, contigs=contigs)

    [m] = read_tsv(tmp_path / "out" / "assessment.tsv")
    assert int(m['num_unaligned_contigs']) == 1
    assert int(m['num_unaligned_regions']) == 2

    unaligned = {r['seq_id']: r for r in parse_fasta(tmp_path / "out" / "unaligned.fa")}
    tail_region_id = f"{tail_contig_id}:{len(tail_contig_seq) + 1}-{len(tail_contig_seq) + 2000}"
    assert unaligned[tail_region_id]['seq'] == foreign_tail
    assert unaligned['contaminant:1-5000']['seq'] == contaminant

    fractions = {r['contig_id']: float(r['contig_fraction_aligned']) for r in rows}
    assert fractions['contaminant'] == 0.0
    assert fractions[tail_contig_id] == pytest.approx(len(tail_contig_seq) / (len(tail_contig_seq) + 2000), abs=1e-3)


# --- repeat awareness ---

def test_repeat_gaps_typed(tmp_path, sim):
    _, rows = run(tmp_path, sim)
    [m] = read_tsv(tmp_path / "out" / "assessment.tsv")
    assert int(m['num_repeat_gaps']) == 2
    assert int(m['repeat_gap_length']) == 2 * len(sim['repeat'])
    assert int(m['num_unique_gaps']) == 0
    assert m['ref_fraction_covered_including_repeats'] == '1.0'
    assert int(m['num_ambiguous_contigs']) == 1
    names = [line.split('\t')[3] for line in (tmp_path / "out" / "gaps.bed").read_text().splitlines()]
    assert names == ['gap_1_repeat', 'gap_2_repeat']
    repeat_id = sim['truth_segments'][-1][0]
    repeat_rows = [r for r in rows if r['contig_id'] == repeat_id]
    assert len(repeat_rows) == 3
    assert sorted(r['is_secondary'] for r in repeat_rows) == ['False', 'True', 'True']
    assert all(r['num_placements'] == '3' for r in repeat_rows)
    # Unique contigs are not ambiguous, even though they flank repeat copies
    assert all(r['num_placements'] == '1' for r in rows if r['contig_id'] != repeat_id)


def test_missing_contig_is_unique_gap(tmp_path, sim):
    # Drop a unique 9 kb contig: its region becomes a unique gap, left as N
    dropped_id, dropped_seq = sim['contigs'][2]
    contigs = [c for c in sim['contigs'] if c[0] != dropped_id]
    output, _ = run(tmp_path, sim, contigs=contigs)
    [m] = read_tsv(tmp_path / "out" / "assessment.tsv")
    assert int(m['num_unique_gaps']) == 1
    assert int(m['unique_gap_length']) == len(dropped_seq)
    assert int(m['num_repeat_gaps']) == 2
    assert output['seq'].count('N') == len(dropped_seq)
