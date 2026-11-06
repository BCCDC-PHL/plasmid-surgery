"""
Find suitable reference plasmids for a plasmid reconstruction, in a plasmid
database (a mash sketch, and the FASTA it was made from).

Two stages:

1. **Screen** (mash): one ``mash dist`` of the reconstruction (sketched as a
   whole) against every database sketch. mash's distance comes from an estimate
   of the Jaccard index J of the two k-mer sets, which is symmetric, and low
   when the sets differ in size even if one contains the other. With both
   sizes, J also gives estimated *containment* in each direction:

       |A ∩ B| = J / (1 + J) × (|A| + |B|)

   so the fraction of the reconstruction's k-mers in the reference, and the
   fraction of the reference's in the reconstruction. Candidates are the top N
   by containment of the reconstruction, and the top N by mash distance.

2. **Confirm** (minimap2): each candidate is aligned to the reconstruction's
   contigs (and optionally the whole assembly's), for exact fractions and
   identity, and ranked by

       score = (fraction of the reconstruction aligned to the reference)
             × (fraction of the reference covered)

   The coverage term is by the reconstruction, or with an assembly, by the
   whole assembly. Both factors are near 1 only for the same plasmid. Against
   the reconstruction alone, a fragment of a large plasmid scores better against
   a smaller relative that shares the fragment than against the large plasmid
   itself; the assembly's coverage reveals the large plasmid.
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import shutil
import subprocess
import sys
import tempfile

from collections import defaultdict
from pathlib import Path

import plasmid_surgery.assessment as assessment
import plasmid_surgery.parsers as parsers

logger = logging.getLogger(__name__)

DEFAULT_NUM_CANDIDATES = 20
DEFAULT_TOP = 3
# Warn when the best candidate accounts for less than this fraction of the reconstruction
WEAK_MATCH_FRACTION = 0.8

CANDIDATE_FIELDNAMES = [
    'rank',
    'ref_id',
    'ref_length',
    'description',
    'mash_distance',
    'mash_shared_hashes',
    'est_reconstruction_containment',
    'est_ref_containment',
    'reconstruction_fraction_aligned',
    'ref_fraction_covered',
    'assembly_ref_fraction_covered',
    'identity',
    'score',
    'selected',
]


class RefSearchError(Exception):
    pass


def run_command(cmd: list[str]) -> str:
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError as e:
        raise RefSearchError(f"{cmd[0]} not found on PATH") from e
    except subprocess.CalledProcessError as e:
        raise RefSearchError(f"`{' '.join(cmd)}` failed (exit code {e.returncode}):\n{e.stderr}") from e
    return result.stdout


# --- screening ---

def sketch_info(sketch_path) -> tuple[int, int, dict]:
    """
    :return: (k-mer size, sketch size, {sequence ID: length}) for a mash sketch
    """
    header = run_command(['mash', 'info', '-H', str(sketch_path)])
    k = s = None
    for line in header.splitlines():
        if 'K-mer size' in line:
            k = int(line.split(':')[1].split()[0])
        elif 'min-hashes per sketch' in line:
            s = int(line.split(':')[1].split()[0])
    if k is None or s is None:
        raise RefSearchError(f"Couldn't read k-mer and sketch size from `mash info -H {sketch_path}`")
    lengths = {}
    for line in run_command(['mash', 'info', '-t', str(sketch_path)]).splitlines():
        if line.startswith('#') or not line.strip():
            continue
        cols = line.split('\t')
        lengths[cols[2]] = int(cols[1])
    return k, s, lengths


def estimate_containment(jaccard: float, size_a: int, size_b: int) -> tuple[float, float]:
    """
    Estimated containment in each direction from a Jaccard index and the two
    set sizes: (fraction of A in B, fraction of B in A).
    """
    if jaccard <= 0 or size_a <= 0 or size_b <= 0:
        return 0.0, 0.0
    intersection = jaccard / (1 + jaccard) * (size_a + size_b)
    return min(1.0, intersection / size_a), min(1.0, intersection / size_b)


def screen(reconstruction_fasta, db_sketch, tmpdir, threads: int = 1) -> list[dict]:
    """
    mash dist of the reconstruction (as one sketch) against every database sketch.

    :return: one dict per database sequence sharing at least one hash, with
             mash_distance, mash_shared_hashes, est_reconstruction_containment,
             est_ref_containment and ref_length
    """
    k, s, lengths = sketch_info(db_sketch)
    contigs = parsers.parse_fasta(reconstruction_fasta)
    reconstruction_kmers = sum(max(0, len(c['seq']) - k + 1) for c in contigs)
    prefix = Path(tmpdir) / 'reconstruction'
    run_command(['mash', 'sketch', '-k', str(k), '-s', str(s), '-o', str(prefix), str(reconstruction_fasta)])
    output = run_command(['mash', 'dist', '-p', str(threads), str(db_sketch), f"{prefix}.msh"])

    rows = []
    for line in output.splitlines():
        ref_id, _, distance, _, shared = line.split('\t')
        shared_hashes, sketch_size = (int(x) for x in shared.split('/'))
        if shared_hashes == 0:
            continue
        ref_length = lengths.get(ref_id, 0)
        in_ref, ref_in = estimate_containment(shared_hashes / sketch_size, reconstruction_kmers, ref_length)
        rows.append({
            'ref_id': ref_id,
            'ref_length': ref_length,
            'mash_distance': float(distance),
            'mash_shared_hashes': f"{shared_hashes}/{sketch_size}",
            'est_reconstruction_containment': round(in_ref, 4),
            'est_ref_containment': round(ref_in, 4),
        })
    return rows


def select_candidates(rows: list[dict], n: int) -> list[dict]:
    """
    The top n by estimated containment of the reconstruction, and the top n by
    mash distance (a reference can be in both).
    """
    by_containment = sorted(rows, key=lambda r: (-r['est_reconstruction_containment'], r['mash_distance']))[:n]
    by_distance = sorted(rows, key=lambda r: (r['mash_distance'], -r['est_reconstruction_containment']))[:n]
    chosen = {}
    for r in by_containment + by_distance:
        chosen.setdefault(r['ref_id'], r)
    return list(chosen.values())


# --- extracting candidates ---

def extract_sequences(fasta_path, ids, index_path=None) -> dict:
    """
    Sequences (and header descriptions) for the given IDs from a FASTA:
    random access with a samtools faidx index if one exists (``<fasta>.fai``,
    or index_path), otherwise one streaming pass. Nothing is written next to
    the FASTA.

    :return: {id: {'seq': ..., 'description': ...}}
    """
    ids = set(ids)
    index_path = Path(index_path) if index_path else Path(f"{fasta_path}.fai")
    if index_path.exists():
        import pysam
        offsets = {}
        with open(index_path) as f:
            for line in f:
                cols = line.split('\t')
                if cols[0] in ids:
                    offsets[cols[0]] = int(cols[2])
        found = {}
        with pysam.FastaFile(str(fasta_path), filepath_index=str(index_path)) as fasta, open(fasta_path, 'rb') as raw:
            for seq_id, offset in offsets.items():
                found[seq_id] = {'seq': fasta.fetch(seq_id).upper(), 'description': _header_before(raw, offset, seq_id)}
        return found

    logger.info(f"No FASTA index for {fasta_path}; reading it through once (index it with `samtools faidx` "
                "for faster lookups)")
    found = {}
    current = None
    chunks = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if current:
                    found[current[0]] = {'seq': ''.join(chunks).upper(), 'description': current[1]}
                    if len(found) == len(ids):
                        return found
                defline = line[1:].strip()
                seq_id = defline.split()[0] if defline else ''
                current = (seq_id, defline[len(seq_id):].strip()) if seq_id in ids else None
                chunks = []
            elif current:
                chunks.append(line.strip())
    if current:
        found[current[0]] = {'seq': ''.join(chunks).upper(), 'description': current[1]}
    return found


def _header_before(raw, offset: int, seq_id: str) -> str:
    """The description on the header line just before byte offset (where a sequence starts, per the .fai)."""
    start = max(0, offset - 4096)
    raw.seek(start)
    chunk = raw.read(offset - start).decode(errors='replace').rstrip('\r\n')
    defline = chunk[chunk.rfind('>') + 1:].strip() if '>' in chunk else ''
    return defline[len(seq_id):].strip() if defline.startswith(seq_id) else ''


# --- confirming ---

def _paf(target_fasta, query_fasta, threads: int) -> list[list[str]]:
    output = run_command(['minimap2', '-c', '-x', 'asm10', '--secondary=no', '-t', str(threads),
                          str(target_fasta), str(query_fasta)])
    return [line.split('\t') for line in output.splitlines() if line]


def _union_length(intervals) -> int:
    return sum(e - s for s, e in assessment.merge_intervals(intervals))


def _target_intervals(cols: list[str], min_deletion_length: int) -> list[tuple[int, int]]:
    """
    Target (ref) intervals a PAF alignment covers, excluding deletions of at
    least min_deletion_length (ref sequence the query doesn't have), as in the
    full assessment. Uses the cg:Z: CIGAR (minimap2 -c).
    """
    cigar = next((t[5:] for t in cols[12:] if t.startswith('cg:Z:')), None)
    if cigar is None:
        return [(int(cols[7]), int(cols[8]))]
    aln = {'ref_start': int(cols[7]), 'cigar_parsed': parsers.parse_cigar(cigar)}
    return assessment.ref_intervals(aln, min_deletion_length)


def confirm(candidate_fasta, ref_length: int, reconstruction_fasta, reconstruction_length: int,
            assembly_fasta=None, threads: int = 1,
            min_deletion_length: int = assessment.DEFAULT_MIN_UNCOVERED_DELETION) -> dict:
    """
    Align the reconstruction (and optionally the assembly) to one candidate.

    Uses a single linear alignment per query (no rotation), so a few bases at
    the reference's origin can be missed; the full assessment
    (``plasmid-surgery``) handles the origin exactly. As there, deletions of
    at least min_deletion_length within an alignment don't count as covered.

    :return: dict with reconstruction_fraction_aligned, ref_fraction_covered,
             assembly_ref_fraction_covered ('NA' without an assembly) and identity
             (gap-compressed, alignment-length weighted)
    """
    query_intervals = defaultdict(list)
    target_intervals = []
    weighted_identity = 0.0
    total_block = 0
    for cols in _paf(candidate_fasta, reconstruction_fasta, threads):
        query_intervals[cols[0]].append((int(cols[2]), int(cols[3])))
        target_intervals += _target_intervals(cols, min_deletion_length)
        block = int(cols[10])
        divergence = next((float(t[5:]) for t in cols[12:] if t.startswith('de:f:')), None)
        if divergence is not None:
            weighted_identity += (1 - divergence) * block
            total_block += block
    aligned = sum(_union_length(v) for v in query_intervals.values())
    result = {
        'reconstruction_fraction_aligned': round(aligned / reconstruction_length, 4) if reconstruction_length else 0.0,
        'ref_fraction_covered': round(_union_length(target_intervals) / ref_length, 4) if ref_length else 0.0,
        'identity': round(100 * weighted_identity / total_block, 3) if total_block else 0.0,
        'assembly_ref_fraction_covered': 'NA',
    }
    if assembly_fasta:
        assembly_target = [iv for c in _paf(candidate_fasta, assembly_fasta, threads)
                           for iv in _target_intervals(c, min_deletion_length)]
        result['assembly_ref_fraction_covered'] = round(_union_length(assembly_target) / ref_length, 4) if ref_length else 0.0
    return result


def score(row: dict) -> float:
    coverage = row['assembly_ref_fraction_covered']
    if coverage == 'NA':
        coverage = row['ref_fraction_covered']
    return round(row['reconstruction_fraction_aligned'] * coverage, 4)


def rank(rows: list[dict]) -> list[dict]:
    """Sort by score, then identity, then mash distance; set 'rank' (1 = best)."""
    ranked = sorted(rows, key=lambda r: (-r['score'], -r['identity'], r['mash_distance']))
    for i, row in enumerate(ranked, start=1):
        row['rank'] = i
    return ranked


def find_refs(reconstruction_fasta, db_sketch, db_fasta, tmpdir, db_fasta_index=None, assembly_fasta=None,
              num_candidates: int = DEFAULT_NUM_CANDIDATES, threads: int = 1) -> tuple[list[dict], dict]:
    """
    Screen the database with mash, confirm candidates with minimap2, and rank them.

    :return: (ranked candidate rows, {ref_id: sequence} for the candidates)
    """
    contigs = parsers.parse_fasta(reconstruction_fasta)
    reconstruction_length = sum(len(c['seq']) for c in contigs)
    rows = screen(reconstruction_fasta, db_sketch, tmpdir, threads)
    logger.info(f"{len(rows)} database sequence(s) share k-mers with the reconstruction")
    candidates = select_candidates(rows, num_candidates)
    if not candidates:
        return [], {}
    sequences = extract_sequences(db_fasta, [c['ref_id'] for c in candidates], db_fasta_index)
    missing = [c['ref_id'] for c in candidates if c['ref_id'] not in sequences]
    if missing:
        logger.warning(f"{len(missing)} candidate(s) in the sketch but not in the FASTA, skipped: {', '.join(missing[:5])}")
    confirmed = []
    for candidate in candidates:
        entry = sequences.get(candidate['ref_id'])
        if entry is None:
            continue
        path = Path(tmpdir) / 'candidates' / f"{len(confirmed)}.fa"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f">{candidate['ref_id']}\n{entry['seq']}\n")
        row = {**candidate, 'ref_length': len(entry['seq']), 'description': entry['description']}
        row.update(confirm(path, len(entry['seq']), reconstruction_fasta, reconstruction_length, assembly_fasta, threads))
        row['score'] = score(row)
        confirmed.append(row)
    return rank(confirmed), {ref_id: e['seq'] for ref_id, e in sequences.items()}


# --- command line ---

def safe_name(seq_id: str) -> str:
    """A sequence ID made safe for use as a file or directory name."""
    return ''.join(c if c.isalnum() or c in '._-' else '_' for c in seq_id) or 'ref'


def add_database_arguments(parser) -> None:
    """
    The plasmid database options other than --db-sketch, shared by plasmid-surgery-refs and
    plasmid-surgery (which adds --db-sketch as an alternative to --ref-plasmid).
    """
    parser.add_argument('--db-fasta', help="The plasmid database's FASTA (IDs matching the sketch's)")
    parser.add_argument('--db-fasta-index', help="samtools faidx index for --db-fasta (default: <db-fasta>.fai, if it exists; "
                                                 "without an index the FASTA is read through once)")
    parser.add_argument('-n', '--num-candidates', type=int, default=DEFAULT_NUM_CANDIDATES,
                        help=f"Candidates to confirm with minimap2: the top N by estimated containment, and the top N by "
                             f"mash distance (default: {DEFAULT_NUM_CANDIDATES})")
    parser.add_argument('--top', type=int, default=DEFAULT_TOP, help=f"Number of best references to keep (default: {DEFAULT_TOP})")
    parser.add_argument('-t', '--threads', type=int, default=1, help="Threads for mash and minimap2 (default: 1)")


def write_candidates(path, ranked: list[dict], top: int) -> None:
    for row in ranked:
        row['selected'] = row['rank'] <= top
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CANDIDATE_FIELDNAMES, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(ranked)


def build_parser() -> argparse.ArgumentParser:
    from plasmid_surgery import __version__
    parser = argparse.ArgumentParser(
        prog='plasmid-surgery-refs',
        description="Find suitable reference plasmids for a plasmid reconstruction in a plasmid database "
                    "(a mash sketch and the FASTA it was made from): screen with mash, confirm with minimap2.")
    parser.add_argument('--plasmid-reconstruction', required=True, help="Plasmid reconstruction: a FASTA of its contigs")
    parser.add_argument('--db-sketch', required=True, help="mash sketch of the plasmid database (.msh)")
    add_database_arguments(parser)
    parser.add_argument('--assembly', help="Full assembly the reconstruction was made from. Ranks by the assembly's coverage "
                                           "of each candidate rather than the reconstruction's, so a reconstruction that's a "
                                           "fragment of a larger plasmid doesn't favour smaller relatives")
    parser.add_argument('-o', '--outdir', default='.', help="Output directory (default: current directory)")
    parser.add_argument('-p', '--prefix', help="Prefix for output filenames (default: no prefix)")
    parser.add_argument('--tmpdir', help="Directory for intermediate files, kept after the run (default: a temporary "
                                         "directory, removed afterwards)")
    parser.add_argument('--log-level', default='info', help="Log level (default: 'info')")
    parser.add_argument('--version', action='version', version=f"%(prog)s {__version__}")
    return parser


def write_outputs(ranked: list[dict], sequences: dict, outdir: Path, prefix: str, top: int) -> None:
    name = (lambda f: f"{prefix}_{f}") if prefix else (lambda f: f)
    outdir.mkdir(parents=True, exist_ok=True)
    write_candidates(outdir / name('candidates.tsv'), ranked, top)
    ref_dir = outdir / name('references')
    if ref_dir.exists():
        shutil.rmtree(ref_dir)
    ref_dir.mkdir()
    for row in ranked[:top]:
        seq = sequences[row['ref_id']]
        with open(ref_dir / f"{row['rank']}_{safe_name(row['ref_id'])}.fa", 'w') as f:
            header = f">{row['ref_id']} {row['description']}".rstrip()
            f.write(header + '\n')
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + '\n')


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.db_fasta:
        parser.error("--db-fasta is required")
    from plasmid_surgery.__main__ import setup_logging
    setup_logging(args.log_level)
    try:
        if args.tmpdir:
            os.makedirs(args.tmpdir, exist_ok=True)
            ranked, sequences = _run(args, args.tmpdir)
        else:
            with tempfile.TemporaryDirectory(prefix='plasmid-surgery-refs_') as tmpdir:
                ranked, sequences = _run(args, tmpdir)
    except RefSearchError as e:
        logger.error(str(e))
        return 1
    write_outputs(ranked, sequences, Path(args.outdir), args.prefix, args.top)
    if not ranked:
        logger.warning("No database sequence shares k-mers with the reconstruction")
    elif ranked[0]['reconstruction_fraction_aligned'] < WEAK_MATCH_FRACTION:
        logger.warning(f"The best candidate ({ranked[0]['ref_id']}) accounts for only "
                       f"{ranked[0]['reconstruction_fraction_aligned']:.0%} of the reconstruction: the database may not "
                       "contain a closely related plasmid")
    for row in ranked[:args.top]:
        logger.info(f"#{row['rank']} {row['ref_id']} ({row['ref_length']:,} bp): score {row['score']}, "
                    f"reconstruction aligned {row['reconstruction_fraction_aligned']}, "
                    f"ref covered {row['ref_fraction_covered']}, identity {row['identity']}%")
    return 0


def _run(args, tmpdir):
    return find_refs(args.plasmid_reconstruction, args.db_sketch, args.db_fasta, tmpdir,
                     db_fasta_index=args.db_fasta_index, assembly_fasta=args.assembly,
                     num_candidates=args.num_candidates, threads=args.threads)


if __name__ == '__main__':
    sys.exit(main())
