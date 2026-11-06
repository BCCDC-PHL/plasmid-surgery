import argparse
import csv
import json
import logging
import os
import sys
import tempfile

from pathlib import Path
from textwrap import wrap

import plasmid_surgery.parsers as parsers
import plasmid_surgery.alignment as alignment
import plasmid_surgery.assessment as assessment
import plasmid_surgery.context as context
import plasmid_surgery.features as features

logger = logging.getLogger(__name__)

# Output files written to --outdir, named '<prefix>_<filename>' if --prefix is given
OUTPUT_FILENAMES = {
    'assessment': 'assessment.tsv',
    'gaps': 'gaps.bed',
    'events': 'events.tsv',
    'unaligned': 'unaligned.fa',
    'summary': 'summary.tsv',
    'scaffold': 'scaffold.fa',
    # Only with --assembly
    'contig_sources': 'contig_sources.tsv',
    'missed_contigs': 'missed_contigs.fa',
    # Only with --ref-features / --features
    'ref_features': 'ref_features.tsv',
    'features': 'features.tsv',
}

SUMMARY_OUTPUT_FIELDNAMES = [
    'contig_id',
    'ref_id',
    'ref_length',
    'reverse_complement',
    'is_mapped',
    'is_supplementary',
    'is_secondary',
    'num_placements',
    'spans_origin',
    'mapq',
    'contig_length',
    'contig_fraction_aligned',
    'alignment_start',
    'alignment_end',
    'alignment_length',
    'percent_identity',
    'included_in_output',
]


def setup_logging(log_level):
    """
    Configure logging at the given level name (e.g. 'info', 'debug').
    """
    numeric_log_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_log_level, int):
        raise ValueError(f"Invalid log level: {log_level}")

    logging.basicConfig(format='%(asctime)s :: %(levelname)s :: %(message)s', datefmt='%Y-%m-%dT%H:%M:%S', level=numeric_log_level)

    return None


def parse_args(argv=None):
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(prog='plasmid-surgery')
    parser.add_argument('--plasmid-reconstruction', required=True, help="Plasmid reconstruction to be repaired.")
    parser.add_argument('--ref-plasmid', required=True, help="Reference plasmid to be used for plasmid surgery.")
    parser.add_argument('--assembly', help="Full assembly the reconstruction was made from. Reports where the rest of the ref plasmid's sequence is in the assembly (contig_sources.tsv, missed_contigs.fa)")
    parser.add_argument('--bins', nargs='+', metavar='FASTA', help="Other bins of the reconstruction tool (e.g. mob_recon's plasmid_*.fasta and chromosome.fasta), used to label assembly contigs by bin. Requires --assembly")
    parser.add_argument('--min-missed-fraction', type=float, default=context.DEFAULT_MIN_MISSED_FRACTION, help=f"With --assembly: minimum fraction of an assembly contig aligned to the ref (at a single location) to call it 'missed' rather than 'shared' (default: {context.DEFAULT_MIN_MISSED_FRACTION})")
    parser.add_argument('--ref-features', nargs='+', metavar='FILE', help="Feature (gene) annotations on the ref plasmid: hAMRonization TSV, GFF3 or BED. Reports how completely the reconstruction (and assembly) contains each (ref_features.tsv)")
    parser.add_argument('--features', nargs='+', metavar='FILE', help="Feature (gene) annotations on the sample's assembly or reconstruction contigs: hAMRonization TSV, GFF3 or BED. Reports, for features on the reconstruction (and with --assembly, on any contig aligning to the ref), where each is and whether it's on the ref (features.tsv)")
    parser.add_argument('--min-feature-coverage', type=float, default=features.DEFAULT_MIN_FEATURE_COVERAGE, help=f"Minimum fraction of a feature covered to call it present / on the ref (default: {features.DEFAULT_MIN_FEATURE_COVERAGE})")
    parser.add_argument('-o', '--outdir', default='.', help="Output directory (default: current directory)")
    parser.add_argument('-p', '--prefix', help="Prefix for output filenames, e.g. '-p sample1' gives 'sample1_assessment.tsv' (default: no prefix)")
    parser.add_argument('--tmpdir', help="Directory for intermediate minimap2 files, kept after the run. If not supplied, a temporary directory is used and removed afterwards.")
    parser.add_argument('--min-indel-length', type=int, default=50, help="Minimum length of an insertion/deletion to report as an event. Deletions at least this long also count as uncovered ref sequence (gaps) (default: 50)")
    parser.add_argument('--min-unaligned-length', type=int, default=100, help="Minimum length of an unaligned contig region to report (default: 100)")
    parser.add_argument('--sample-id', help="Sample ID, used as the scaffold fasta header and the sample_id in assessment output. If not supplied, uses --prefix if given, otherwise the plasmid reconstruction filename (minus file extension)")
    parser.add_argument('--output-bases-per-line', type=int, default=70, help="Number of bases per line in output fasta files (default: 70)")
    parser.add_argument('--log-level', default='info', help="Log level (default: 'info')")
    args = parser.parse_args(argv)
    if args.bins and not args.assembly:
        parser.error("--bins requires --assembly")

    return args


def output_paths(outdir, prefix=None) -> dict:
    """
    Paths of all output files.

    :param outdir: Output directory
    :param prefix: Optional filename prefix
    :return: dict of output name (keys of OUTPUT_FILENAMES) -> Path
    """
    outdir = Path(outdir)
    return {
        name: outdir / (f"{prefix}_{filename}" if prefix else filename)
        for name, filename in OUTPUT_FILENAMES.items()
    }


def summary_row(aligned_contig: dict, ref_id: str, ref_length: int, contig_fraction_aligned: dict) -> dict:
    """
    Build a summary output row for one aligned segment.
    """
    is_mapped = aligned_contig.get('is_mapped', None)
    row = {
        'contig_id': aligned_contig.get('contig_id', None),
        'ref_id': ref_id,
        'ref_length': ref_length,
        'reverse_complement': aligned_contig.get('is_reverse', None),
        'is_mapped': is_mapped,
        'is_supplementary': aligned_contig.get('is_supplementary', None),
        'is_secondary': aligned_contig.get('is_secondary', None),
        'num_placements': aligned_contig.get('num_placements', None),
        'spans_origin': aligned_contig.get('spans_origin', False),
        'mapq': aligned_contig.get('mapq', None),
        'contig_length': aligned_contig.get('query_full_length', None),
        'contig_fraction_aligned': contig_fraction_aligned.get(aligned_contig.get('contig_id'), 0.0),
        'alignment_start': aligned_contig.get('ref_start', None),
        'alignment_end': aligned_contig.get('ref_end', None),
        'alignment_length': aligned_contig.get('query_alignment_length', None),
        'percent_identity': aligned_contig.get('percent_identity', None),
        'included_in_output': aligned_contig.get('included_in_output', False),
    }
    if is_mapped and row['alignment_end'] > ref_length:
        # Alignment wraps around the origin: alignment_end < alignment_start
        row['alignment_end'] -= ref_length
    if not is_mapped:
        for field in ['reverse_complement', 'is_supplementary', 'is_secondary', 'num_placements', 'spans_origin', 'mapq', 'alignment_start', 'alignment_end', 'alignment_length', 'percent_identity']:
            row[field] = 'NA'

    return row


def deduplicate_contigs(contigs: list[dict]) -> list[dict]:
    """
    Drop contigs whose sequence (in either orientation) repeats an earlier one.
    Some tools write a contig once per alignment hit (e.g. TaDReP, for a contig
    split across the reference origin), which would otherwise be counted twice.
    """
    seen = set()
    unique = []
    for contig in contigs:
        key = context.sequence_key(contig['seq'])
        if key in seen:
            logger.warning(f"Reconstruction contig {contig['seq_id']} duplicates the sequence of an earlier contig; ignoring it")
            continue
        seen.add(key)
        unique.append(contig)
    return unique


def summary_sort_key(row: dict):
    """
    Sort mapped rows by alignment start; unmapped rows last.
    """
    if row['is_mapped']:
        return (0, row['alignment_start'])
    return (1, 0)


def main(argv=None):
    args = parse_args(argv)
    setup_logging(args.log_level)
    outputs = output_paths(args.outdir, args.prefix)
    os.makedirs(args.outdir, exist_ok=True)

    if args.tmpdir:
        os.makedirs(args.tmpdir, exist_ok=True)
        return run(args, outputs, args.tmpdir)
    with tempfile.TemporaryDirectory(prefix='plasmid-surgery_') as tmpdir:
        return run(args, outputs, tmpdir)


def run(args, outputs: dict, tmpdir) -> int:
    """
    Align, assess and write outputs.

    :param args: Parsed command-line arguments
    :param outputs: Output paths, from output_paths
    :param tmpdir: Directory for intermediate files
    :return: Exit code
    """

    plasmid_reconstruction_contigs = deduplicate_contigs(parsers.parse_fasta(args.plasmid_reconstruction))
    ref_plasmids = parsers.parse_fasta(args.ref_plasmid)
    if not ref_plasmids:
        logger.error(f"No sequences found in ref plasmid file: {args.ref_plasmid}")
        return 1
    if len(ref_plasmids) > 1:
        logger.warning(f"Ref plasmid file contains {len(ref_plasmids)} sequences; using only the first: {ref_plasmids[0]['seq_id']}")
    ref_plasmid = ref_plasmids[0]
    ref_plasmid_length = len(ref_plasmid['seq'])
    ref_plasmid_id = ref_plasmid['seq_id']

    if args.sample_id:
        sample_id = args.sample_id
    elif args.prefix:
        sample_id = args.prefix
    else:
        plasmid_reconstruction_basename = os.path.basename(args.plasmid_reconstruction)
        sample_id = plasmid_reconstruction_basename.split('.')[0]

    logger.info(f"Aligning {len(plasmid_reconstruction_contigs)} reconstruction contig(s) vs. ref plasmid: {ref_plasmid_id}")
    try:
        recon_alignments = alignment.align_contigs(plasmid_reconstruction_contigs, ref_plasmid, tmpdir)
    except alignment.AlignmentError as e:
        logger.error(str(e))
        return 1
    aligned_contigs = recon_alignments['placements']
    alternative_placements = recon_alignments['alternatives']
    failed_alignments = recon_alignments['failed']

    output_seq = alignment.build_output_sequence(aligned_contigs, ref_plasmid_length, alternative_placements)

    assessment_result = assessment.assess(
        plasmid_reconstruction_contigs,
        aligned_contigs,
        ref_plasmid_length,
        min_indel_length=args.min_indel_length,
        min_unaligned_length=args.min_unaligned_length,
        alternative_placements=alternative_placements,
    )
    metrics = {'sample_id': sample_id, 'ref_id': ref_plasmid_id, **assessment_result['metrics']}

    assembly_contigs = None
    context_result = None
    if args.assembly:
        assembly_contigs = parsers.parse_fasta(args.assembly)
        bins = context.bin_labels(args.bins) if args.bins else None
        try:
            context_result = context.assess_assembly_context(
                plasmid_reconstruction_contigs, recon_alignments, assembly_contigs, ref_plasmid, tmpdir,
                bins=bins, min_missed_fraction=args.min_missed_fraction,
                min_deletion_length=args.min_indel_length,
            )
        except alignment.AlignmentError as e:
            logger.error(str(e))
            return 1
        metrics.update(context_result['metrics'])
        context.write_contig_sources(outputs['contig_sources'], context_result['contigs'])
        context.write_missed_contigs(outputs['missed_contigs'], context_result['missed_contigs'], bins, args.output_bases_per_line)
    else:
        bins = None

    if args.ref_features:
        ref_feature_list, ref_feature_sequences = features.read_features(args.ref_features)
        mapped, unmapped = features.map_features(ref_feature_list, [ref_plasmid], ref_feature_sequences)
        if unmapped:
            log = logger.warning if not mapped else logger.info
            log(f"{unmapped} ref feature(s) are not on the ref plasmid ({ref_plasmid_id}) and were ignored")
        rows, feature_metrics = features.assess_ref_features(
            mapped, ref_plasmid_length, recon_alignments,
            assembly_alignments=context_result['placements'] if context_result else None,
            labels=context_result['locations'] if context_result else None,
            min_coverage=args.min_feature_coverage,
            min_deletion_length=args.min_indel_length,
        )
        features.write_rows(outputs['ref_features'], rows, features.REF_FEATURES_FIELDNAMES)
        metrics.update(feature_metrics)

    if args.features:
        sample_feature_list, sample_feature_sequences = features.read_features(args.features)
        contigs_of_interest = {c['seq_id']: c for c in plasmid_reconstruction_contigs}
        classifications = {cid: 'in_reconstruction' for cid in contigs_of_interest}
        placements = recon_alignments['placements']
        known_sequences = dict(sample_feature_sequences)
        if context_result:
            contigs_of_interest.update(context_result['contigs_by_id'])
            classifications.update(context_result['classifications'])
            placements = context_result['placements']
            known_sequences = {**{c['seq_id']: c['seq'] for c in assembly_contigs}, **known_sequences}
        mapped, unmapped = features.map_features(sample_feature_list, list(contigs_of_interest.values()), known_sequences)
        logger.info(f"{len(mapped)} feature(s) on the reconstruction{' or contigs aligning to the ref' if context_result else ''}; "
                    f"{unmapped} elsewhere (not reported)")
        if not mapped and not context_result:
            logger.warning("No features matched the reconstruction's contigs. If they were annotated on the assembly, "
                           "and the reconstruction tool renamed contigs, add --assembly so they can be matched by sequence")
        rows, feature_metrics = features.assess_sample_features(
            mapped, contigs_of_interest, placements, ref_plasmid_length, classifications,
            bins=bins, min_coverage=args.min_feature_coverage,
        )
        features.write_rows(outputs['features'], rows, features.SAMPLE_FEATURES_FIELDNAMES)
        metrics.update(feature_metrics)

    logger.info("Assessment: " + json.dumps(metrics))
    assessment.write_metrics(outputs['assessment'], metrics, assessment.METRICS_FIELDNAMES + context.CONTEXT_METRICS_FIELDNAMES + features.FEATURE_METRICS_FIELDNAMES)
    assessment.write_gaps_bed(outputs['gaps'], assessment_result['gaps'], ref_plasmid_id, ref_plasmid_length)
    assessment.write_events(outputs['events'], assessment_result['events'])
    assessment.write_unaligned_fasta(outputs['unaligned'], assessment_result['unaligned_regions'], args.output_bases_per_line)

    contig_fraction_aligned = assessment_result['contig_fraction_aligned']
    summary_output = [summary_row(a, ref_plasmid_id, ref_plasmid_length, contig_fraction_aligned) for a in aligned_contigs + alternative_placements + failed_alignments]
    summary_output_sorted = sorted(summary_output, key=summary_sort_key)
    with open(outputs['summary'], 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_OUTPUT_FIELDNAMES, delimiter='\t', extrasaction='ignore', quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        for row in summary_output_sorted:
            writer.writerow(row)

    output_seq_wrapped = wrap(output_seq, args.output_bases_per_line)
    with open(outputs['scaffold'], 'w') as f:
        f.write(f">{sample_id}\n")
        for line in output_seq_wrapped:
            f.write(f"{line}\n")

    return 0


if __name__ == '__main__':
    sys.exit(main())
