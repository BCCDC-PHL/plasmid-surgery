"""
Feature (gene) annotations: on the reference plasmid (--ref-features), and on
the sample's contigs (--features).

Annotations are read into a common feature record, from any of:

- hAMRonization TSV (covers abricate with any database, AMRFinderPlus, RGI,
  ResFinder, ...)
- GFF3 (e.g. Bakta, Prokka, NCBI)
- BED (3+ columns; name in column 4 and strand in column 6 if present)

No gene class is treated specially: categories (e.g. drug class, feature type)
are passed through as given.

Feature records are dicts with keys: seq_id, start, end (0-based, half-open),
strand ('+', '-' or '.'), name, product, category, source.
"""
from __future__ import annotations

import csv
import logging
import re

from collections import defaultdict
from urllib.parse import unquote

import plasmid_surgery.assessment as assessment
import plasmid_surgery.context as context

logger = logging.getLogger(__name__)

DEFAULT_MIN_FEATURE_COVERAGE = 0.9

HAMRONIZATION_REQUIRED_COLUMNS = {'input_sequence_id', 'input_gene_start', 'input_gene_stop', 'gene_symbol'}

# GFF3 feature types that describe whole sequences rather than features on them
GFF3_SKIPPED_TYPES = {'region', 'source', 'databank_entry', 'chromosome', 'plasmid', 'contig'}

FEATURE_FIELDNAMES = ['name', 'product', 'category', 'source']

REF_FEATURES_FIELDNAMES = FEATURE_FIELDNAMES + [
    'ref_start',
    'ref_end',
    'strand',
    'length',
    'reconstruction_fraction',
    'inferred_fraction',
    'reconstruction_status',
    'assembly_fraction',
    'assembly_status',
    'found_in',
    'reads_fraction',   # with reads: fraction of the feature with read depth >= --min-read-depth
]

SAMPLE_FEATURES_FIELDNAMES = FEATURE_FIELDNAMES + [
    'contig_id',
    'contig_start',
    'contig_end',
    'strand',
    'length',
    'contig_classification',
    'bin',
    'fraction_on_ref',
    'ref_status',
    'ref_start',
    'ref_end',
]

FEATURE_METRICS_FIELDNAMES = [
    'num_ref_features',
    'num_ref_features_present',
    'num_ref_features_partial',
    'num_ref_features_absent',
    'num_ref_features_inferred',
    'num_ref_features_missed',
    'num_features',
    'num_features_not_on_ref',
    'num_features_on_missed_contigs',
]

STRANDS = {'+': '+', '-': '-', 'sense': '+', 'antisense': '-', '1': '+', '-1': '-'}


class FeatureFormatError(Exception):
    pass


# --- parsing ---

def _feature(seq_id, start, end, strand, name, product='', category='', source=''):
    start, end = int(start), int(end)
    if start > end:
        start, end = end, start
    return {
        'seq_id': seq_id,
        'start': start,
        'end': end,
        'strand': STRANDS.get(str(strand).strip(), '.'),
        'name': name or '',
        'product': product or '',
        'category': category or '',
        'source': source or '',
    }


def parse_hamronization(path) -> list[dict]:
    """
    hAMRonization TSV. Gene coordinates are taken as 1-based inclusive (as
    from abricate and AMRFinderPlus).
    """
    features = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        missing = HAMRONIZATION_REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise FeatureFormatError(f"{path}: not hAMRonization TSV (missing columns: {', '.join(sorted(missing))})")
        for row in reader:
            if not row['input_gene_start'] or not row['input_gene_stop']:
                continue
            source = ':'.join(x for x in (row.get('analysis_software_name'), row.get('reference_database_name')) if x)
            first, last = sorted((int(row['input_gene_start']), int(row['input_gene_stop'])))
            features.append(_feature(
                row['input_sequence_id'], first - 1, last,
                row.get('strand_orientation', ''), row['gene_symbol'],
                product=row.get('gene_name', ''), category=row.get('drug_class', ''), source=source,
            ))
    return features


def _gff3_attributes(column: str) -> dict:
    attributes = {}
    for item in column.strip().split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            attributes[key.strip()] = unquote(value.strip())
    return attributes


def parse_gff3(path) -> tuple[list[dict], dict]:
    """
    GFF3, including an optional ##FASTA section. Features that are only parents
    of other features (e.g. a 'gene' with a 'CDS' child) are dropped, so each
    gene is reported once.

    :return: (features, sequences from the ##FASTA section, by seq ID)
    """
    records = []
    sequences = {}
    with open(path) as f:
        in_fasta = False
        seq_id = None
        for line in f:
            line = line.rstrip('\n')
            if in_fasta:
                if line.startswith('>'):
                    seq_id = line[1:].split()[0]
                    sequences[seq_id] = []
                elif seq_id:
                    sequences[seq_id].append(line.strip())
                continue
            if line.startswith('##FASTA'):
                in_fasta = True
                continue
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) != 9:
                raise FeatureFormatError(f"{path}: not GFF3 (expected 9 tab-separated columns): {line[:80]}")
            if cols[2] in GFF3_SKIPPED_TYPES:
                continue
            records.append((cols, _gff3_attributes(cols[8])))

    parent_ids = {parent for _, attributes in records for parent in attributes.get('Parent', '').split(',') if parent}
    features = []
    for cols, attributes in records:
        if attributes.get('ID') in parent_ids:
            continue
        name = attributes.get('gene') or attributes.get('Name') or attributes.get('locus_tag') or attributes.get('ID', '')
        features.append(_feature(
            cols[0], int(cols[3]) - 1, int(cols[4]), cols[6], name,
            product=attributes.get('product', ''), category=cols[2], source=cols[1],
        ))
    return features, {k: ''.join(v) for k, v in sequences.items()}


def parse_bed(path) -> list[dict]:
    features = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 3:
                raise FeatureFormatError(f"{path}: not BED (fewer than 3 columns): {line[:80]}")
            features.append(_feature(
                cols[0], cols[1], cols[2], cols[5] if len(cols) > 5 else '.',
                cols[3] if len(cols) > 3 else f"{cols[0]}:{cols[1]}-{cols[2]}",
                source='bed',
            ))
    return features


def detect_format(path) -> str:
    """
    :return: 'hamronization', 'gff3', 'bed', or 'empty' (no content, e.g.
             hAMRonization output for a tool that found nothing)
    """
    with open(path) as f:
        if not any(line.strip() for line in f):
            return 'empty'
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('##gff-version'):
                return 'gff3'
            cols = line.rstrip('\n').split('\t')
            if HAMRONIZATION_REQUIRED_COLUMNS <= set(cols):
                return 'hamronization'
            if line.startswith('#') or line.startswith(('track', 'browser')):
                continue
            if len(cols) == 9 and cols[3].isdigit() and cols[4].isdigit():
                return 'gff3'
            if len(cols) >= 3 and cols[1].isdigit() and cols[2].isdigit():
                return 'bed'
            break
    raise FeatureFormatError(f"{path}: unrecognized feature format (expected hAMRonization TSV, GFF3 or BED)")


def read_features(paths) -> tuple[list[dict], dict]:
    """
    Read features from any number of files, detecting each file's format.

    :return: (features, sequences embedded in the files (GFF3 ##FASTA), by seq ID)
    """
    features = []
    sequences = {}
    for path in paths:
        fmt = detect_format(path)
        if fmt == 'empty':
            logger.info(f"No features in {path} (empty file)")
            continue
        if fmt == 'hamronization':
            file_features = parse_hamronization(path)
        elif fmt == 'gff3':
            file_features, file_sequences = parse_gff3(path)
            sequences.update(file_sequences)
        else:
            file_features = parse_bed(path)
        logger.info(f"Read {len(file_features)} feature(s) from {path} ({fmt})")
        features += file_features
    return features, sequences


# --- mapping features onto sequences ---

def _strip_version(seq_id: str) -> str:
    return re.sub(r'\.\d+$', '', seq_id)


def map_features(features: list[dict], targets: list[dict], known_sequences: dict) -> tuple[list[dict], int]:
    """
    Map features onto target sequences (contigs, or the ref).

    A feature's sequence is identified by its seq_id: directly, if a target has
    that ID (or the same ID without a '.N' version suffix); otherwise through
    known_sequences (e.g. a GFF3 ##FASTA section, or the assembly), matched to
    a target by sequence in either orientation. Features on a target that is
    the reverse complement of the annotated sequence are flipped.

    :param targets: Sequence dicts (seq_id, seq) to map onto
    :param known_sequences: seq ID -> sequence for annotated sequences
    :return: (features with 'target_id' set, coordinates on that target;
              number of features that could not be mapped)
    """
    by_id = {t['seq_id']: t for t in targets}
    by_unversioned_id = {_strip_version(t['seq_id']): t for t in targets}
    by_key = {context.sequence_key(t['seq']): t for t in targets}
    mapped = []
    unmapped = 0
    for feature in features:
        target = by_id.get(feature['seq_id']) or by_unversioned_id.get(_strip_version(feature['seq_id']))
        flip = False
        source_seq = known_sequences.get(feature['seq_id'])
        if target is not None and source_seq is not None and source_seq.upper() != target['seq'].upper():
            target = None   # same ID but different sequence: match by sequence instead
        if target is None and source_seq is not None:
            target = by_key.get(context.sequence_key(source_seq))
            if target is not None:
                flip = source_seq.upper() != target['seq'].upper()
        if target is None:
            unmapped += 1
            continue
        mapped_feature = dict(feature, target_id=target['seq_id'])
        if flip:
            length = len(target['seq'])
            mapped_feature.update({
                'start': length - feature['end'],
                'end': length - feature['start'],
                'strand': {'+': '-', '-': '+'}.get(feature['strand'], '.'),
            })
        mapped.append(mapped_feature)
    return mapped, unmapped


# --- reference features ---

def _fraction(mask: bytearray, start: int, end: int) -> float:
    length = end - start
    return round(sum(mask[start:end]) / length, 4) if length > 0 else 0.0


def _status(fraction: float, min_coverage: float) -> str:
    if fraction >= min_coverage:
        return 'present'
    if fraction > 0:
        return 'partial'
    return 'absent'


def _with_inferred(fraction: float, inferred_fraction: float, min_coverage: float) -> str:
    """present / partial / absent, or 'inferred' if present only counting alternative placements."""
    status = _status(fraction, min_coverage)
    if status != 'present' and _status(fraction + inferred_fraction, min_coverage) == 'present':
        return 'inferred'
    return status


def assess_ref_features(ref_features: list[dict], ref_length: int, reconstruction: dict,
                        assembly_alignments: list[dict] = None, labels: dict = None,
                        assembly_alternatives: list[dict] = (),
                        min_coverage: float = DEFAULT_MIN_FEATURE_COVERAGE,
                        min_deletion_length: int = assessment.DEFAULT_MIN_UNCOVERED_DELETION) -> tuple[list[dict], dict]:
    """
    For each feature annotated on the ref: how much of it the reconstruction
    covers, how much more is inferred from alternative placements (repeat
    copies), and, with the assembly, how much the whole assembly covers and
    which bins those contigs are in.

    :param ref_features: Features mapped onto the ref (map_features)
    :param reconstruction: alignment.align_contigs result for the reconstruction
    :param assembly_alignments: Placement alignments of all assembly contigs (or None)
    :param assembly_alternatives: Their alternative placements (repeat copies)
    :param labels: contig ID -> location label for assembly contigs ('reconstruction', bin, ...)
    :return: (rows for REF_FEATURES_FIELDNAMES, metrics)
    """
    reconstruction_mask = context.ref_mask(reconstruction['placements'], ref_length, min_deletion_length)
    inferred_mask = context.ref_mask(reconstruction['alternatives'], ref_length, min_deletion_length)
    assembly_mask = (context.ref_mask(assembly_alignments, ref_length, min_deletion_length)
                     if assembly_alignments is not None else None)
    assembly_inferred_mask = context.ref_mask(list(assembly_alternatives), ref_length, min_deletion_length)

    rows = []
    counts = defaultdict(int)
    for feature in sorted(ref_features, key=lambda f: (f['start'], f['end'])):
        start, end = feature['start'], min(feature['end'], ref_length)
        reconstruction_fraction = _fraction(reconstruction_mask, start, end)
        either = bytearray(a | b for a, b in zip(reconstruction_mask[start:end], inferred_mask[start:end]))
        inferred_fraction = round(_fraction(either, 0, end - start) - reconstruction_fraction, 4)
        status = _with_inferred(reconstruction_fraction, inferred_fraction, min_coverage)
        counts[status] += 1
        row = {
            **{k: feature[k] for k in FEATURE_FIELDNAMES},
            'ref_start': start,
            'ref_end': end,
            'strand': feature['strand'],
            'length': end - start,
            'reconstruction_fraction': reconstruction_fraction,
            'inferred_fraction': inferred_fraction,
            'reconstruction_status': status,
            'assembly_fraction': 'NA',
            'assembly_status': 'NA',
            'found_in': 'NA',
        }
        if assembly_mask is not None:
            assembly_fraction = _fraction(assembly_mask, start, end)
            assembly_either = bytearray(a | b for a, b in zip(assembly_mask[start:end], assembly_inferred_mask[start:end]))
            assembly_inferred_fraction = round(_fraction(assembly_either, 0, end - start) - assembly_fraction, 4)
            assembly_status = _with_inferred(assembly_fraction, assembly_inferred_fraction, min_coverage)
            found_in = sorted({
                labels.get(a['contig_id'], 'assembly') for a in list(assembly_alignments) + list(assembly_alternatives)
                if any(s < end and e > start for s, e in assessment.covered_ref_intervals(a, ref_length, min_deletion_length))
            })
            row.update({
                'assembly_fraction': assembly_fraction,
                'assembly_status': assembly_status,
                'found_in': ','.join(found_in) if found_in else '-',
            })
            if status != 'present' and assembly_status == 'present':
                counts['missed'] += 1
        rows.append(row)

    metrics = {
        'num_ref_features': len(rows),
        'num_ref_features_present': counts['present'],
        'num_ref_features_partial': counts['partial'],
        'num_ref_features_absent': counts['absent'],
        'num_ref_features_inferred': counts['inferred'],
        'num_ref_features_missed': counts['missed'] if assembly_mask is not None else 'NA',
    }
    return rows, metrics


# --- sample features ---

def matched_query_intervals_on_contig(aln: dict) -> list[tuple[int, int]]:
    """
    Intervals of the contig (original orientation) whose bases are matched to
    ref bases (=, X, M), excluding insertions relative to the ref.
    """
    intervals = []
    q = assessment.leading_clip_length(aln)
    for op, length in aln['cigar_parsed']:
        if op in assessment.REF_ALIGNED_OPS:
            intervals.append(assessment.to_contig_coords(aln, q, q + length))
            q += length
        elif op == 'insertion_to_ref':
            q += length
    return assessment.merge_intervals(intervals)


def ref_pos_at_contig_pos(aln: dict, contig_pos: int, ref_length: int):
    """
    Ref position aligned to a contig position (original contig orientation),
    folded into [0, ref_length), or None if the position isn't matched to the ref.
    """
    q_target = aln['query_full_length'] - 1 - contig_pos if aln.get('is_reverse') else contig_pos
    q = assessment.leading_clip_length(aln)
    r = aln['ref_start']
    for op, length in aln['cigar_parsed']:
        if op in assessment.REF_ALIGNED_OPS:
            if q <= q_target < q + length:
                return (r + q_target - q) % ref_length
            q += length
            r += length
        elif op == 'insertion_to_ref':
            q += length
        elif op in assessment.REF_ONLY_OPS:
            r += length
    return None


def assess_sample_features(features: list[dict], contigs: dict, placements: list[dict], ref_length: int,
                           classifications: dict, bins: dict = None,
                           min_coverage: float = DEFAULT_MIN_FEATURE_COVERAGE) -> tuple[list[dict], dict]:
    """
    For each feature on a contig of interest (the reconstruction's contigs, and
    with the assembly, any contig aligning to the ref): the contig's
    classification and bin, how much of the feature is matched to the ref, and
    where on the ref.

    :param features: Features mapped onto contigs (map_features)
    :param contigs: contig ID -> contig dict, for the contigs of interest
    :param placements: Placement alignments of those contigs
    :param classifications: contig ID -> classification ('in_reconstruction', 'missed', ...)
    :param bins: Optional sequence key -> bin label
    :return: (rows for SAMPLE_FEATURES_FIELDNAMES, metrics)
    """
    alignments_by_contig = defaultdict(list)
    for aln in placements:
        alignments_by_contig[aln['contig_id']].append(aln)

    rows = []
    for feature in sorted(features, key=lambda f: (f['target_id'], f['start'])):
        contig_id = feature['target_id']
        if contig_id not in contigs:
            continue
        start, end = feature['start'], feature['end']
        length = end - start
        alns = alignments_by_contig.get(contig_id, [])
        on_ref = 0
        best_aln, best_overlap = None, 0
        for aln in alns:
            overlap = sum(max(0, min(e, end) - max(s, start)) for s, e in matched_query_intervals_on_contig(aln))
            if overlap > best_overlap:
                best_aln, best_overlap = aln, overlap
        merged = assessment.merge_intervals([iv for aln in alns for iv in matched_query_intervals_on_contig(aln)])
        on_ref = sum(max(0, min(e, end) - max(s, start)) for s, e in merged)
        fraction = round(on_ref / length, 4) if length > 0 else 0.0
        ref_start = ref_end = 'NA'
        if best_aln is not None:
            positions = [p for p in (ref_pos_at_contig_pos(best_aln, x, ref_length)
                                     for x in (start, end - 1)) if p is not None]
            if len(positions) == 2:
                ref_start, ref_end = min(positions), max(positions) + 1
        contig = contigs[contig_id]
        rows.append({
            **{k: feature[k] for k in FEATURE_FIELDNAMES},
            'contig_id': contig_id,
            'contig_start': start,
            'contig_end': end,
            'strand': feature['strand'],
            'length': length,
            'contig_classification': classifications.get(contig_id, 'NA'),
            'bin': (bins.get(context.sequence_key(contig['seq']), 'unbinned') if bins is not None else 'NA'),
            'fraction_on_ref': fraction,
            'ref_status': {'present': 'on_ref', 'partial': 'partial', 'absent': 'not_on_ref'}[_status(fraction, min_coverage)],
            'ref_start': ref_start,
            'ref_end': ref_end,
        })

    metrics = {
        'num_features': len(rows),
        'num_features_not_on_ref': sum(1 for r in rows if r['ref_status'] == 'not_on_ref'),
        'num_features_on_missed_contigs': sum(1 for r in rows if r['contig_classification'] == 'missed'),
    }
    return rows, metrics


def write_rows(path, rows: list[dict], fieldnames: list) -> None:
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
