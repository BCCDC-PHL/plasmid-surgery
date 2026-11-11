"""
Core SNPs: every sample's reconstruction placed on the shared reference, and SNPs counted at the same
reference positions for every pair.

Pairwise SNP distances (plasmid_surgery.compare) compare each pair over whatever the two share, so pairs
are compared over different stretches. Here, each reconstruction is aligned to the reference (one
minimap2 call, asm20, with cs tags) and **called** at each reference position: its base there, or
nothing. A position isn't called for a sample:

- outside its alignments, or in an alignment less than ``min_identity`` identical (a divergent homolog,
  as for pairwise distances);
- where its alignments cover the position more than once (which copy is which is ambiguous);
- within compare.END_TRIM bp of an alignment's end or of one of its contigs' ends;
- in a deletion relative to the reference;
- in a repeat of the reference itself (found by aligning the reference to itself): a sample's single
  copy may have aligned to either.

The **core** is the reference positions called in at least ``min_fraction`` of the samples (1.0: in all).
Each pair is compared at the core positions both called (with min_fraction 1.0, the whole core), and
SNPs in clusters (compare.CLUSTER_MIN within compare.CLUSTER_WINDOW bp along the reference) are counted
apart. The reference is a member too (called everywhere). The reference is aligned as linear, so the
END_TRIM bp either side of its origin are never called.
"""
from __future__ import annotations

import csv
import logging
import math
import re
import subprocess

from itertools import accumulate

from pathlib import Path

import plasmid_surgery.alignment as alignment
import plasmid_surgery.compare as compare
import plasmid_surgery.joins as joins

logger = logging.getLogger(__name__)

NOT_CALLED = ord('-')
UNCALLED_RE = re.compile(rb'-+')
CHUNK = 64

CORE_PAIR_FIELDNAMES = ['a', 'b', 'compared_positions', 'snps', 'snps_in_clusters', 'snps_outside_clusters']
CORE_SAMPLE_FIELDNAMES = ['sample', 'called_positions', 'called_fraction', 'core_without_it']


def call_samples(reference: dict, samples: list[dict], tmpdir, threads: int = 1,
                 min_identity: float = compare.SNP_MIN_IDENTITY) -> list[bytearray]:
    """
    Each sample's base at each reference position (upper case), or NOT_CALLED.

    :param reference: {'seq_id', 'seq'}
    :param samples: [{'name', 'contigs': [{'seq_id', 'seq'}]}]
    """
    tmpdir = Path(tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)
    ref_seq = reference['seq'].upper()
    n = len(ref_seq)
    ref_fasta, query_fasta = tmpdir / 'core_reference.fa', tmpdir / 'core_samples.fa'
    alignment.write_fasta(ref_fasta, [('ref', ref_seq)])
    alignment.write_fasta(query_fasta, [(f'm{i}_{k}', c['seq']) for i, s in enumerate(samples)
                                        for k, c in enumerate(s['contigs'])])
    cmd = ['minimap2', '-c', '--cs', '-x', 'asm20', '-P', '-t', str(threads), str(ref_fasta), str(query_fasta)]
    logger.info(f"Aligning {len(samples)} reconstructions to the reference {reference['seq_id']} for core SNPs")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError as e:
        raise alignment.AlignmentError("minimap2 not found on PATH") from e
    except subprocess.CalledProcessError as e:
        raise alignment.AlignmentError(f"minimap2 failed aligning the reconstructions to the reference:\n{e.stderr[-2000:]}")
    by_sample = {}
    for line in result.stdout.splitlines():
        f = line.split('\t')
        i, k = map(int, f[0][1:].split('_'))
        block = int(f[10])
        de = next((float(t[5:]) for t in f[12:] if t.startswith('de:f:')), None)
        ident = 1 - de if de is not None else (int(f[9]) / block if block else 0)
        cs = next((t[5:] for t in f[12:] if t.startswith('cs:Z:')), None)
        if block < compare.MIN_ALIGNMENT or cs is None:
            continue
        by_sample.setdefault(i, []).append({'qk': k, 'tk': 0, 'qlen': int(f[1]), 'qs': int(f[2]), 'qe': int(f[3]),
                                            'strand': f[4], 'ts': int(f[7]), 'te': int(f[8]), 'block': block,
                                            'cs': cs, 'identity': ident})
    ref_repeats = joins.repeat_intervals([{'seq_id': 'ref', 'seq': ref_seq}], tmpdir / 'reference_repeats', threads)
    ref_repeats = ref_repeats.get('ref', [])
    calls = []
    for i in range(len(samples)):
        calls.append(call_one(ref_seq, by_sample.get(i, []), ref_repeats, min_identity))
    return calls


def call_one(ref_seq: str, alns: list[dict], ref_repeats: list, min_identity: float) -> bytearray:
    """One sample's calls along the reference, from its alignments to it."""
    n = len(ref_seq)
    ref_bytes = ref_seq.encode()
    call = bytearray([NOT_CALLED]) * n
    # Reference positions this sample covers more than once (any alignment): ambiguous
    multi = compare._multiply_covered([(a['ts'], a['te']) for a in alns])
    q_mask = {}
    for a in alns:
        q_mask.setdefault(a['qk'], []).append((a['qs'], a['qe']))
    q_mask = {k: compare._multiply_covered(v) for k, v in q_mask.items()}
    kept = compare.one_to_one([a for a in alns if a['identity'] >= min_identity])
    T = compare.END_TRIM
    for a in kept:
        ts, te, qs, qe = a['ts'], a['te'], a['qs'], a['qe']
        t_excl = compare._Excluded(multi + list(ref_repeats) + [(ts, ts + T), (te - T, te), (0, T), (n - T, n)])
        q_excl = compare._Excluded(q_mask.get(a['qk'], []) + [(0, T), (a['qlen'] - T, a['qlen'])])
        forward = a['strand'] == '+'
        t, qoff = ts, 0

        def qpos(off):
            return qs + off if forward else qe - 1 - off
        for op in compare.CS_RE.findall(a['cs']):
            kind = op[0]
            if kind == ':':
                m = int(op[1:])
                parts = [(s - t, e - t) for s, e in t_excl.within(t, t + m)]
                if forward:
                    parts += [(s - qpos(qoff), e - qpos(qoff)) for s, e in q_excl.within(qpos(qoff), qpos(qoff) + m)]
                else:
                    lo, hi = qpos(qoff + m - 1), qpos(qoff) + 1
                    parts += [(qpos(qoff) - (e - 1), qpos(qoff) - s + 1) for s, e in q_excl.within(lo, hi)]
                pos = 0
                for s, e in compare._merge(parts) + [(m, m)]:
                    s, e = max(s, 0), min(e, m)
                    if s > pos:
                        call[t + pos:t + s] = ref_bytes[t + pos:t + s]
                    pos = max(pos, e)
                t += m
                qoff += m
            elif kind == '*':
                if t not in t_excl and qpos(qoff) not in q_excl:
                    call[t] = ord(op[2].upper())
                t += 1
                qoff += 1
            elif kind == '+':
                qoff += len(op) - 1
            else:   # '-': deleted in the sample: not called
                t += len(op) - 1
    return call


def core_snps(reference: dict, samples: list[dict], tmpdir, threads: int = 1, min_fraction: float = 1.0,
              min_identity: float = compare.SNP_MIN_IDENTITY) -> dict:
    """
    :return: dict with 'names' (the reference first, then the samples), 'core_positions' (reference positions
             in the core), 'calls' (per member, its bytes at the core positions), 'pairs' ({(i, j): CORE_PAIR
             fields}), 'samples' (CORE_SAMPLE_FIELDNAMES rows), 'reference_length', 'variable_sites'
    """
    sample_calls = call_samples(reference, samples, tmpdir, threads, min_identity)
    n = len(reference['seq'])
    ref_call = bytearray(reference['seq'].upper().encode())
    m = len(samples)
    need = max(1, math.ceil(round(min_fraction * m, 6))) if m else 0
    # Samples not called at each position, from runs of NOT_CALLED
    diff = [0] * (n + 1)
    for c in sample_calls:
        for run in UNCALLED_RE.finditer(c):
            diff[run.start()] += 1
            diff[run.end()] -= 1
    missing = list(accumulate(diff))[:n]
    core_runs = []
    for p in range(n):
        if m - missing[p] >= need:
            if core_runs and core_runs[-1][1] == p:
                core_runs[-1][1] = p + 1
            else:
                core_runs.append([p, p + 1])
    core = [p for a, b in core_runs for p in range(a, b)]
    members = [ref_call] + sample_calls
    names = [reference['seq_id']] + [s['name'] for s in samples]
    core_calls = [b''.join(bytes(c[a:b]) for a, b in core_runs) for c in members]
    # Called masks over the core as bitsets (popcount of two ANDed = positions both called)
    to_bits = bytes.maketrans(bytes(range(256)), b'1' * 256)
    to_bits = to_bits[:NOT_CALLED] + b'0' + to_bits[NOT_CALLED + 1:]
    masks = [int(c[::-1].translate(to_bits) or b'0', 2) for c in core_calls]
    # Variable sites: wherever a called member differs from the reference (any difference between two members
    # means at least one differs from the reference)
    variable = sorted(set().union(*(differing(core_calls[0], c) for c in core_calls[1:]))) if len(core_calls) > 1 else []
    pairs = {}
    for i in range(len(members)):
        for j in range(i + 1, len(members)):
            ci, cj = core_calls[i], core_calls[j]
            diffs = [core[x] for x in variable if ci[x] != NOT_CALLED and cj[x] != NOT_CALLED and ci[x] != cj[x]]
            in_clusters = len(compare.clustered(diffs))
            pairs[(i, j)] = {'a': names[i], 'b': names[j], 'compared_positions': (masks[i] & masks[j]).bit_count(),
                             'snps': len(diffs), 'snps_in_clusters': in_clusters,
                             'snps_outside_clusters': len(diffs) - in_clusters}
    # How big the strict core would be without each sample: the positions only it doesn't call
    only = list(accumulate(1 if x == 1 else 0 for x in missing))
    only_missing = []
    for c in sample_calls:
        total = 0
        for run in UNCALLED_RE.finditer(c):
            total += only[run.end() - 1] - (only[run.start() - 1] if run.start() else 0)
        only_missing.append(total)
    sample_rows = []
    for i, (s, c) in enumerate(zip(samples, sample_calls)):
        called = n - c.count(NOT_CALLED)
        sample_rows.append({'sample': s['name'], 'called_positions': called,
                            'called_fraction': round(called / n, 4) if n else 0.0,
                            'core_without_it': len(core) + only_missing[i] if need == m else 'NA'})
    return {'names': names, 'core_positions': core, 'calls': core_calls, 'pairs': pairs, 'samples': sample_rows,
            'reference_length': n, 'variable_sites': len(variable), 'min_fraction': min_fraction}


def differing(ref: bytes, other: bytes) -> set:
    """Positions where other is called and differs from ref (compared in chunks: identical stretches are skipped)."""
    out = set()
    for x in range(0, len(ref), CHUNK):
        a, b = ref[x:x + CHUNK], other[x:x + CHUNK]
        if a == b:
            continue
        for y, (u, v) in enumerate(zip(a, b)):
            if u != v and v != NOT_CALLED:
                out.add(x + y)
    return out


def write_outputs(result: dict, fasta_path, pairs_path, samples_path) -> None:
    """The core alignment (FASTA; '-' where a member isn't called), the pairwise core SNPs and the per-sample table."""
    with open(fasta_path, 'w') as f:
        for name, c in zip(result['names'], result['calls']):
            f.write(f">{name}\n")
            seq = c.decode()
            for x in range(0, len(seq), 80):
                f.write(seq[x:x + 80] + "\n")
    for path, rows, fields in ((pairs_path, list(result['pairs'].values()), CORE_PAIR_FIELDNAMES),
                               (samples_path, result['samples'], CORE_SAMPLE_FIELDNAMES)):
        with open(path, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(rows)
