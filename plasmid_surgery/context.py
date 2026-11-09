"""
Whole-assembly context: where does the rest of the reference plasmid's sequence
sit in the assembly, if not in the plasmid reconstruction?

Given the full assembly the reconstruction was made from (and optionally the
other bins of the reconstruction tool, e.g. mob_recon's other plasmid_*.fasta
and chromosome.fasta), every assembly contig that aligns to the reference is
classified as:

- in_reconstruction: part of the reconstruction being assessed
- shared:  only part of the contig aligns to the reference. Likely a contig
           from the chromosome or another plasmid that shares some sequence
           with this one (e.g. a transposon), not a binning error.
- repeat:  most of the contig aligns to the reference, but equally well at
           more than one place (e.g. an IS element), so where it belongs is
           ambiguous: it may be on this plasmid, elsewhere, or both.
- missed:  most of the contig aligns to the reference, at a single location.
           Likely a plasmid contig the reconstruction tool put elsewhere.

Read depth from the assembler's contig headers, where available, is reported
relative to the reconstruction's contigs as supporting evidence: a missed
contig at a similar relative depth to the reconstruction is consistent with
belonging to the same plasmid.

Contigs are matched between the reconstruction, the assembly and the bins by
sequence (either orientation), not by name, so tools that rename or reorient
contigs are handled.
"""
from __future__ import annotations

import csv
import hashlib
import logging
import re

from collections import defaultdict
from pathlib import Path
from textwrap import wrap

import plasmid_surgery.alignment as alignment
import plasmid_surgery.assessment as assessment
import plasmid_surgery.parsers as parsers

logger = logging.getLogger(__name__)

DEFAULT_MIN_MISSED_FRACTION = 0.8

CONTIG_SOURCES_FIELDNAMES = [
    'contig_id',
    'contig_length',
    'depth',
    'depth_relative_to_reconstruction',
    'bin',
    'classification',
    'contig_fraction_on_ref',
    'ref_bases_covered',
    'ref_start',
    'ref_end',
    'num_placements',
    'mapq',
    'percent_identity',
]

# Assembly-context metrics, added to the assessment metrics ('NA' when no assembly is given)
CONTEXT_METRICS_FIELDNAMES = [
    'assembly_ref_bases_covered',
    'assembly_ref_fraction_covered',
    'assembly_ref_fraction_covered_including_repeats',
    'assembly_unique_gap_length',
    'missed_ref_bases',
    'num_missed_contigs',
    'missed_contig_bases',
    'num_shared_contigs',
    'num_repeat_contigs',
    'num_reconstruction_contigs_not_in_assembly',
]

DEPTH_PATTERNS = [
    re.compile(r'depth=([0-9.]+)x?'),        # Unicycler
    re.compile(r'_cov_([0-9.]+)'),           # SPAdes
    re.compile(r'\bcov=([0-9.]+)'),          # Shovill and others
]


def sequence_key(seq: str) -> str:
    """
    Orientation-independent identifier for a sequence.
    """
    seq = seq.upper()
    canonical = min(seq, alignment.revcomp(seq))
    return hashlib.md5(canonical.encode()).hexdigest()


def parse_depth(contig: dict):
    """
    Read depth from an assembler's contig header, if present (Unicycler
    'depth=1.2x', SPAdes '_cov_12.3', 'cov=12.3').

    :return: depth as float, or None
    """
    header = f"{contig['seq_id']} {contig.get('description', '')}"
    for pattern in DEPTH_PATTERNS:
        match = pattern.search(header)
        if match:
            try:
                return float(match.group(1).rstrip('.'))
            except ValueError:
                return None
    return None


def bin_labels(bin_paths: list) -> dict:
    """
    Map sequence key -> bin label (FASTA filename without extension) for every
    contig in the given bin FASTA files. A contig in more than one bin gets
    a comma-separated label.
    """
    labels = defaultdict(list)
    for path in bin_paths:
        label = Path(path).name
        for ext in ('.gz', '.fasta', '.fa', '.fna', '.fas'):
            if label.endswith(ext):
                label = label[:-len(ext)]
        for contig in parsers.parse_fasta(path):
            labels[sequence_key(contig['seq'])].append(label)
    return {key: ','.join(sorted(set(v))) for key, v in labels.items()}


def ref_mask(alignments: list[dict], ref_length: int,
             min_deletion_length: int = assessment.DEFAULT_MIN_UNCOVERED_DELETION) -> bytearray:
    """
    Per-base mask of the reference covered by the given alignments (excluding
    deletions of at least min_deletion_length; see assessment.ref_intervals).
    """
    mask = bytearray(ref_length)
    for aln in alignments:
        for start, end in assessment.covered_ref_intervals(aln, ref_length, min_deletion_length):
            mask[start:end] = b'\x01' * (end - start)
    return mask


def classify(num_placements: int, fraction_on_ref: float, min_missed_fraction: float) -> str:
    """
    Classify an assembly contig (not in the reconstruction) that aligns to the ref.
    """
    if fraction_on_ref < min_missed_fraction:
        return 'shared'
    if num_placements > 1:
        return 'repeat'
    return 'missed'


def length_weighted_mean_depth(contigs: list[dict]):
    """
    Length-weighted mean read depth of the contigs with a depth in their header,
    or None if none have one.
    """
    weighted = 0.0
    total_length = 0
    for contig in contigs:
        depth = parse_depth(contig)
        if depth is not None:
            weighted += depth * len(contig['seq'])
            total_length += len(contig['seq'])
    return weighted / total_length if total_length else None


def assess_assembly_context(reconstruction_contigs: list[dict], reconstruction_alignments: dict,
                            assembly_contigs: list[dict], ref: dict, tmpdir,
                            bins: dict = None, min_missed_fraction: float = DEFAULT_MIN_MISSED_FRACTION,
                            min_deletion_length: int = assessment.DEFAULT_MIN_UNCOVERED_DELETION,
                            threads: int = 1) -> dict:
    """
    Classify the assembly contigs that align to the reference, and summarize how
    much of the reference the assembly could have provided beyond the
    reconstruction.

    :param reconstruction_contigs: Contigs of the reconstruction being assessed
    :param reconstruction_alignments: alignment.align_contigs result for the reconstruction
    :param assembly_contigs: All contigs of the assembly the reconstruction came from
    :param ref: Reference plasmid
    :param tmpdir: Directory for intermediate files
    :param bins: Optional sequence key -> bin label (see bin_labels)
    :param min_missed_fraction: Minimum fraction of a contig aligned to the ref
                                (at a single location) to call it 'missed'
    :param threads: minimap2 threads
    :return: dict with 'metrics' (CONTEXT_METRICS_FIELDNAMES), 'contigs' (rows
             for CONTIG_SOURCES_FIELDNAMES), 'missed_contigs' (contig dicts),
             'placements' (placement alignments of all assembly contigs aligning
             to the ref, reconstruction included), 'alternatives' (their
             alternative placements), 'contigs_by_id' (those
             contigs), 'classifications' and 'locations' (contig ID ->
             classification, and -> bin label or 'reconstruction')
    """
    ref_length = len(ref['seq'])
    reconstruction_keys = {sequence_key(c['seq']) for c in reconstruction_contigs}
    assembly_keys = {sequence_key(c['seq']) for c in assembly_contigs}
    not_in_assembly = len(reconstruction_keys - assembly_keys)
    if not_in_assembly:
        logger.warning(f"{not_in_assembly} reconstruction contig(s) not found in the assembly: "
                       "was the reconstruction made from this assembly?")

    other_contigs = [c for c in assembly_contigs if sequence_key(c['seq']) not in reconstruction_keys]
    logger.info(f"Aligning {len(other_contigs)} other assembly contig(s) vs. ref plasmid: {ref['seq_id']}")
    other = alignment.align_contigs(other_contigs, ref, tmpdir, label='assembly', threads=threads)

    contigs_by_id = {c['seq_id']: c for c in reconstruction_contigs + other_contigs}
    reconstruction_depth = length_weighted_mean_depth(reconstruction_contigs)
    rows = []
    missed_contigs = []
    missed_alignments = []
    counts = defaultdict(int)
    for source, result in (('reconstruction', reconstruction_alignments), ('other', other)):
        by_contig = defaultdict(list)
        for aln in result['placements']:
            by_contig[aln['contig_id']].append(aln)
        for contig_id, alns in by_contig.items():
            contig = contigs_by_id[contig_id]
            contig_length = len(contig['seq'])
            merged = assessment.merge_intervals([assessment.query_interval_on_contig(a) for a in alns])
            fraction_on_ref = sum(e - s for s, e in merged) / contig_length if contig_length else 0.0
            num_placements = max(a.get('num_placements', 1) for a in alns)
            if source == 'reconstruction':
                classification = 'in_reconstruction'
            else:
                classification = classify(num_placements, fraction_on_ref, min_missed_fraction)
                counts[classification] += 1
                if classification == 'missed':
                    missed_contigs.append(contig)
                    missed_alignments += alns
            primary = next((a for a in alns if not a.get('is_supplementary')), alns[0])
            ref_end = primary['ref_end'] - ref_length if primary['ref_end'] > ref_length else primary['ref_end']
            depth = parse_depth(contig)
            rows.append({
                'contig_id': contig_id,
                'contig_length': contig_length,
                'depth': depth if depth is not None else 'NA',
                'depth_relative_to_reconstruction': (round(depth / reconstruction_depth, 3)
                                                     if depth is not None and reconstruction_depth else 'NA'),
                'bin': (bins.get(sequence_key(contig['seq']), 'unbinned') if bins is not None else 'NA'),
                'classification': classification,
                'contig_fraction_on_ref': round(fraction_on_ref, 4),
                'ref_bases_covered': sum(ref_mask(alns, ref_length, min_deletion_length)),
                'ref_start': primary['ref_start'],
                'ref_end': ref_end,
                'num_placements': num_placements,
                'mapq': primary.get('mapq'),
                'percent_identity': primary.get('percent_identity'),
            })

    rows.sort(key=lambda r: (r['ref_start'], r['contig_id']))

    all_placements = reconstruction_alignments['placements'] + other['placements']
    all_alternatives = reconstruction_alignments['alternatives'] + other['alternatives']
    assembly_mask = ref_mask(all_placements, ref_length, min_deletion_length)
    assembly_covered = sum(assembly_mask)
    gaps = assessment.find_typed_gaps(all_placements, all_alternatives, ref_length, min_deletion_length)
    repeat_gap_length = sum(g['end'] - g['start'] for g in gaps if g['type'] == 'repeat')
    unique_gap_length = sum(g['end'] - g['start'] for g in gaps if g['type'] == 'unique')

    reconstruction_mask = ref_mask(reconstruction_alignments['placements'], ref_length, min_deletion_length)
    missed_mask = ref_mask(missed_alignments, ref_length, min_deletion_length)
    missed_ref_bases = sum(1 for m, r in zip(missed_mask, reconstruction_mask) if m and not r)

    metrics = {
        'assembly_ref_bases_covered': assembly_covered,
        'assembly_ref_fraction_covered': round(assembly_covered / ref_length, 4) if ref_length else 0.0,
        'assembly_ref_fraction_covered_including_repeats': round((assembly_covered + repeat_gap_length) / ref_length, 4) if ref_length else 0.0,
        'assembly_unique_gap_length': unique_gap_length,
        'missed_ref_bases': missed_ref_bases,
        'num_missed_contigs': counts['missed'],
        'missed_contig_bases': sum(len(c['seq']) for c in missed_contigs),
        'num_shared_contigs': counts['shared'],
        'num_repeat_contigs': counts['repeat'],
        'num_reconstruction_contigs_not_in_assembly': not_in_assembly,
    }

    classifications = {r['contig_id']: r['classification'] for r in rows}
    locations = {}
    for r in rows:
        if r['classification'] == 'in_reconstruction':
            locations[r['contig_id']] = 'reconstruction'
        else:
            locations[r['contig_id']] = r['bin'] if r['bin'] != 'NA' else 'assembly'
    aligned_contigs = {cid: contigs_by_id[cid] for cid in classifications}

    return {
        'metrics': metrics,
        'contigs': rows,
        'missed_contigs': missed_contigs,
        'placements': all_placements,
        'alternatives': all_alternatives,
        'contigs_by_id': aligned_contigs,
        'classifications': classifications,
        'locations': locations,
    }


def write_contig_sources(path, rows: list[dict]) -> None:
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CONTIG_SOURCES_FIELDNAMES, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def write_missed_contigs(path, contigs: list[dict], bins: dict = None, bases_per_line: int = 70) -> None:
    with open(path, 'w') as f:
        for contig in contigs:
            label = bins.get(sequence_key(contig['seq']), 'unbinned') if bins is not None else None
            header = f">{contig['seq_id']}"
            if label:
                header += f" bin={label}"
            f.write(header + "\n")
            for line in wrap(contig['seq'], bases_per_line):
                f.write(line + "\n")
