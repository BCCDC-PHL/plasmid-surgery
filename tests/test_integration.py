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


# --- whole-assembly mode ---

@pytest.fixture
def misbinned(tmp_path, sim):
    """
    Assembly where the reconstruction tool put one unique plasmid contig in the
    chromosome bin. The chromosome bin also has a contig sharing 2 kb with the
    plasmid, one containing the IS repeat, and one unrelated contig.
    """
    rng = random.Random(7)
    plasmid_contigs = [(f"{cid} depth=2.00x", seq) for cid, seq in sim['contigs']]
    misbinned_id, misbinned_seq = plasmid_contigs[2]
    reconstruction = [c for c in plasmid_contigs if c[0] != misbinned_id]
    contig_0_seq = dict(sim['truth_segments'])['contig_0']
    chromosome = [
        (misbinned_id, misbinned_seq),
        ("chrom_shared depth=1.00x", random_seq(rng, 9000) + contig_0_seq[5000:7000] + random_seq(rng, 9000)),
        ("chrom_IS depth=1.00x", random_seq(rng, 15000) + sim['repeat'] + random_seq(rng, 5000)),
        ("chrom_plain depth=1.00x", random_seq(rng, 30000)),
    ]
    paths = {name: tmp_path / f"{name}.fasta" for name in ('reconstruction', 'chromosome', 'assembly', 'ref')}
    write_fasta(paths['reconstruction'], reconstruction)
    write_fasta(paths['chromosome'], chromosome)
    write_fasta(paths['assembly'], chromosome[::-1] + reconstruction)
    write_fasta(paths['ref'], [('ref', sim['ref'])])
    return paths, misbinned_id.split()[0], misbinned_seq


def run_with_assembly(tmp_path, paths, bins=True, extra_args=()):
    out_dir = tmp_path / "out"
    args = ['--plasmid-reconstruction', str(paths['reconstruction']), '--ref-plasmid', str(paths['ref']),
            '--assembly', str(paths['assembly']), '-o', str(out_dir), '--log-level', 'warning', *extra_args]
    if bins:
        args += ['--bins', str(paths['reconstruction']), str(paths['chromosome'])]
    assert main(args) == 0
    [m] = read_tsv(out_dir / "assessment.tsv")
    sources = {r['contig_id']: r for r in read_tsv(out_dir / "contig_sources.tsv")}
    return m, sources, out_dir


def test_assembly_mode_finds_misbinned_contig(tmp_path, misbinned):
    paths, misbinned_id, misbinned_seq = misbinned
    m, sources, out_dir = run_with_assembly(tmp_path, paths)

    missed = sources[misbinned_id]
    assert missed['classification'] == 'missed'
    assert missed['bin'] == 'chromosome'
    assert float(missed['contig_fraction_on_ref']) == 1.0
    assert float(missed['depth_relative_to_reconstruction']) == pytest.approx(1.0)

    assert int(m['unique_gap_length']) == len(misbinned_seq)          # reconstruction alone
    assert int(m['assembly_unique_gap_length']) == 0                    # assembly has everything
    assert m['assembly_ref_fraction_covered_including_repeats'] == '1.0'
    assert int(m['missed_ref_bases']) == len(misbinned_seq)
    assert int(m['num_missed_contigs']) == 1
    assert int(m['num_reconstruction_contigs_not_in_assembly']) == 0

    [missed_fasta] = parse_fasta(out_dir / "missed_contigs.fa")
    assert missed_fasta['seq_id'] == misbinned_id
    assert missed_fasta['description'] == 'bin=chromosome'
    assert missed_fasta['seq'] == misbinned_seq


def test_assembly_mode_classifies_shared_contigs(tmp_path, misbinned):
    paths, _, _ = misbinned
    m, sources, _ = run_with_assembly(tmp_path, paths)
    # Contigs that only partly match the plasmid aren't binning errors
    assert sources['chrom_shared']['classification'] == 'shared'
    assert sources['chrom_IS']['classification'] == 'shared'
    assert float(sources['chrom_IS']['depth_relative_to_reconstruction']) == pytest.approx(0.5)
    assert 'chrom_plain' not in sources
    assert int(m['num_shared_contigs']) == 2
    # Reconstruction contigs are listed too
    assert {r['classification'] for r in sources.values()} == {'in_reconstruction', 'missed', 'shared'}
    assert all(r['bin'] == 'reconstruction' for r in sources.values() if r['classification'] == 'in_reconstruction')


def test_assembly_mode_without_bins(tmp_path, misbinned):
    paths, misbinned_id, _ = misbinned
    _, sources, out_dir = run_with_assembly(tmp_path, paths, bins=False)
    assert sources[misbinned_id]['classification'] == 'missed'
    assert all(r['bin'] == 'NA' for r in sources.values())
    [missed_fasta] = parse_fasta(out_dir / "missed_contigs.fa")
    assert missed_fasta['description'] == ''


def test_assembly_mode_reconstruction_contig_not_in_assembly(tmp_path, misbinned):
    paths, _, _ = misbinned
    assembly = [(c['seq_id'], c['seq']) for c in parse_fasta(paths['assembly'])]
    write_fasta(paths['assembly'], [c for c in assembly if c[0] != 'contig_0'])
    m, _, _ = run_with_assembly(tmp_path, paths)
    assert int(m['num_reconstruction_contigs_not_in_assembly']) == 1


def test_assessment_has_na_context_columns_without_assembly(tmp_path, sim):
    run(tmp_path, sim)
    [m] = read_tsv(tmp_path / "out" / "assessment.tsv")
    assert m['missed_ref_bases'] == 'NA'
    assert m['assembly_ref_fraction_covered'] == 'NA'


def test_bins_require_assembly(tmp_path, sim):
    with pytest.raises(SystemExit):
        main(['--plasmid-reconstruction', 'r.fa', '--ref-plasmid', 'ref.fa', '--bins', 'b.fa'])
