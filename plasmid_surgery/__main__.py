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

logger = logging.getLogger(__name__)

# Output files written to --outdir, named '<prefix>_<filename>' if --prefix is given
OUTPUT_FILENAMES = {
    'assessment': 'assessment.tsv',
    'gaps': 'gaps.bed',
    'events': 'events.tsv',
    'unaligned': 'unaligned.fa',
    'summary': 'summary.tsv',
    'scaffold': 'scaffold.fa',
}

SUMMARY_OUTPUT_FIELDNAMES = [
    'contig_id',
    'ref_id',
    'ref_length',
    'reverse_complement',
    'is_mapped',
    'is_supplementary',
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
    parser.add_argument('-o', '--outdir', default='.', help="Output directory (default: current directory)")
    parser.add_argument('-p', '--prefix', help="Prefix for output filenames, e.g. '-p sample1' gives 'sample1_assessment.tsv' (default: no prefix)")
    parser.add_argument('--tmpdir', help="Directory for intermediate minimap2 files, kept after the run. If not supplied, a temporary directory is used and removed afterwards.")
    parser.add_argument('--min-indel-length', type=int, default=50, help="Minimum length of an insertion/deletion to report as an event (default: 50)")
    parser.add_argument('--min-unaligned-length', type=int, default=100, help="Minimum length of an unaligned contig region to report (default: 100)")
    parser.add_argument('--sample-id', help="Sample ID, used as the scaffold fasta header and the sample_id in assessment output. If not supplied, uses --prefix if given, otherwise the plasmid reconstruction filename (minus file extension)")
    parser.add_argument('--output-bases-per-line', type=int, default=70, help="Number of bases per line in output fasta files (default: 70)")
    parser.add_argument('--log-level', default='info', help="Log level (default: 'info')")
    args = parser.parse_args(argv)

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
        for field in ['reverse_complement', 'is_supplementary', 'spans_origin', 'mapq', 'alignment_start', 'alignment_end', 'alignment_length', 'percent_identity']:
            row[field] = 'NA'

    return row


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

    plasmid_reconstruction_contigs = parsers.parse_fasta(args.plasmid_reconstruction)
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

    aligned_contigs = []
    failed_alignments = []
    for contig in plasmid_reconstruction_contigs:
        contig_id = contig['seq_id']
        logger.info(f"Aligning contig: {contig_id} vs. ref plasmid: {ref_plasmid_id}")
        try:
            aligned_segments = alignment.align_contig_vs_circular_ref(contig, ref_plasmid, tmpdir)
        except alignment.AlignmentError as e:
            logger.error(str(e))
            return 1

        if len(aligned_segments) == 0:
            logger.info(f"Zero alignments for contig: {contig_id} vs. ref_plasmid: {ref_plasmid_id}")
            continue
        if any(s['is_supplementary'] for s in aligned_segments):
            logger.info(f"Supplementary alignment found for contig: {contig_id}")

        for aligned_segment in aligned_segments:
            aligned_segment_qc_stats = alignment.alignment_qc(aligned_segment)
            if not aligned_segment_qc_stats['overall_qc_pass']:
                logger.info(f"Alignment failed quality filters: {contig_id}")
                logger.info("Alignment QC stats: " + json.dumps(aligned_segment_qc_stats['qc_stats']))
                aligned_segment['included_in_output'] = False
                failed_alignments.append(aligned_segment)
                continue

            aligned_contigs.append(aligned_segment)

    output_seq = alignment.build_output_sequence(aligned_contigs, ref_plasmid_length)

    assessment_result = assessment.assess(
        plasmid_reconstruction_contigs,
        aligned_contigs,
        ref_plasmid_length,
        min_indel_length=args.min_indel_length,
        min_unaligned_length=args.min_unaligned_length,
    )
    metrics = {'sample_id': sample_id, 'ref_id': ref_plasmid_id, **assessment_result['metrics']}
    logger.info("Assessment: " + json.dumps(metrics))
    assessment.write_metrics(outputs['assessment'], metrics)
    assessment.write_gaps_bed(outputs['gaps'], assessment_result['gaps'], ref_plasmid_id, ref_plasmid_length)
    assessment.write_events(outputs['events'], assessment_result['events'])
    assessment.write_unaligned_fasta(outputs['unaligned'], assessment_result['unaligned_regions'], args.output_bases_per_line)

    contig_fraction_aligned = assessment_result['contig_fraction_aligned']
    summary_output = [summary_row(a, ref_plasmid_id, ref_plasmid_length, contig_fraction_aligned) for a in aligned_contigs + failed_alignments]
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
