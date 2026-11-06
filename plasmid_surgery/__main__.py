import argparse
import csv
import json
import logging
import os
import sys

from textwrap import wrap

import plasmid_surgery.parsers as parsers
import plasmid_surgery.alignment as alignment

logger = logging.getLogger(__name__)

SUMMARY_OUTPUT_FIELDNAMES = [
    'contig_id',
    'ref_id',
    'ref_length',
    'reverse_complement',
    'is_mapped',
    'is_supplementary',
    'mapq',
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
    parser.add_argument('--tmpdir', default="./plasmid-surgery_tmp", help="temp dir to use for minimap2 alignments (default: './plasmid-surgery_tmp')")
    parser.add_argument('--summary-output', default='./plasmid-surgery-summary.tsv', help="Summary file for plasmid surgery process. (default: './plasmid-surgery-summary.tsv')")
    parser.add_argument('--output', default='./plasmid-surgery-output.fa', help="Path to output file for repaired plasmid reconstruction (default: './plasmid-surgery-output.fa')")
    parser.add_argument('--output-seq-id', help="Sequence ID to use in the output fasta header. If not supplied, uses plasmid reconstruction filename (minus file extension)")
    parser.add_argument('--output-bases-per-line', type=int, default=70, help="Number of bases per line in the output fasta (default: 70)")
    parser.add_argument('--log-level', default='info', help="Log level (default: 'info')")
    args = parser.parse_args(argv)

    return args


def ensure_parent_dir(path):
    """
    Create the parent directory of path, if it has one and it doesn't exist.
    """
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def summary_row(aligned_contig: dict, ref_id: str, ref_length: int) -> dict:
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
        'mapq': aligned_contig.get('mapq', None),
        'alignment_start': aligned_contig.get('ref_start', None),
        'alignment_end': aligned_contig.get('ref_end', None),
        'alignment_length': aligned_contig.get('query_alignment_length', None),
        'percent_identity': aligned_contig.get('percent_identity', None),
        'included_in_output': aligned_contig.get('included_in_output', False),
    }
    if not is_mapped:
        for field in ['reverse_complement', 'is_supplementary', 'mapq', 'alignment_start', 'alignment_end', 'alignment_length', 'percent_identity']:
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
    ensure_parent_dir(args.output)
    ensure_parent_dir(args.summary_output)

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

    if args.output_seq_id:
        output_seq_id = args.output_seq_id
    else:
        plasmid_reconstruction_basename = os.path.basename(args.plasmid_reconstruction)
        output_seq_id = plasmid_reconstruction_basename.split('.')[0]

    aligned_contigs = []
    failed_alignments = []
    for contig in plasmid_reconstruction_contigs:
        contig_id = contig['seq_id']
        logger.info(f"Aligning contig: {contig_id} vs. ref plasmid: {ref_plasmid_id}")
        try:
            aligned_segments = alignment.align_contig_vs_ref(contig, ref_plasmid, args.tmpdir)
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

    summary_output = [summary_row(a, ref_plasmid_id, ref_plasmid_length) for a in aligned_contigs + failed_alignments]
    summary_output_sorted = sorted(summary_output, key=summary_sort_key)
    with open(args.summary_output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_OUTPUT_FIELDNAMES, delimiter='\t', extrasaction='ignore', quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        for row in summary_output_sorted:
            writer.writerow(row)

    output_seq_wrapped = wrap(output_seq, args.output_bases_per_line)
    with open(args.output, 'w') as f:
        f.write(f">{output_seq_id}\n")
        for line in output_seq_wrapped:
            f.write(f"{line}\n")

    return 0


if __name__ == '__main__':
    sys.exit(main())
