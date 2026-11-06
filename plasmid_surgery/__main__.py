import argparse
import csv
import json
import logging
import os
import sys
import tempfile

from pathlib import Path
from textwrap import wrap

from plasmid_surgery import __version__
import plasmid_surgery.parsers as parsers
import plasmid_surgery.alignment as alignment
import plasmid_surgery.assessment as assessment
import plasmid_surgery.context as context
import plasmid_surgery.features as features
import plasmid_surgery.refs as refs

logger = logging.getLogger(__name__)

# Output files written to <outdir>/<ref_id>/ for each reference, named '<prefix>_<filename>' if --prefix is given
OUTPUT_FILENAMES = {
    'reference': 'reference.fa',
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

# Written to <outdir>/ itself
TOP_LEVEL_FILENAMES = {
    'assessment': 'assessment.tsv',   # one row per reference, ranked
    'candidates': 'candidates.tsv',   # with --db-sketch only
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


def build_parser() -> argparse.ArgumentParser:
    """
    The command-line parser (also used to generate the documentation's CLI reference).
    """
    parser = argparse.ArgumentParser(
        prog='plasmid-surgery',
        description="Assess a plasmid reconstruction (a set of contigs) against one or more reference plasmids: where each "
                    "contig belongs, what's missing or extra, and optionally where the rest of the plasmid is in the assembly. "
                    "References are given (--ref-plasmid), or found in a plasmid database (--db-sketch, --db-fasta). Each "
                    "reference's outputs go to <outdir>/<ref_id>/, and <outdir>/assessment.tsv has one row per reference, ranked.")
    parser.add_argument('--plasmid-reconstruction', required=True, help="Plasmid reconstruction to assess: a FASTA of its contigs")
    source = parser.add_argument_group("references (one of)")
    exclusive = source.add_mutually_exclusive_group(required=True)
    exclusive.add_argument('--ref-plasmid', nargs='+', metavar='FASTA', help="Reference plasmid(s) to assess against: FASTA file(s); every sequence in them is a separate reference")
    exclusive.add_argument('--db-sketch', help="mash sketch of a plasmid database to find references in (with --db-fasta; see plasmid-surgery-refs)")
    refs.add_database_arguments(parser.add_argument_group("plasmid database (with --db-sketch)"))
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
    parser.add_argument('--report', action='store_true', help="Also write report.html (all references on linear maps) and a map.svg per reference; the same as running plasmid-surgery-plot on the output directory")
    parser.add_argument('--log-level', default='info', help="Log level (default: 'info')")
    parser.add_argument('--version', action='version', version=f"%(prog)s {__version__}")
    return parser


def parse_args(argv=None):
    """
    Parse command-line arguments.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.bins and not args.assembly:
        parser.error("--bins requires --assembly")
    if args.db_sketch and not args.db_fasta:
        parser.error("--db-sketch requires --db-fasta")

    return args


def output_paths(outdir, prefix=None, filenames=OUTPUT_FILENAMES) -> dict:
    """
    Paths of output files.

    :param outdir: Output directory (for a reference's outputs, <outdir>/<ref_id>)
    :param prefix: Optional filename prefix
    :return: dict of output name (keys of filenames) -> Path
    """
    outdir = Path(outdir)
    return {
        name: outdir / (f"{prefix}_{filename}" if prefix else filename)
        for name, filename in filenames.items()
    }


def unique_dirnames(ref_ids: list[str]) -> list[str]:
    """A directory name per reference: its ID made safe, with _2, _3, ... for repeats."""
    names, seen = [], {}
    for ref_id in ref_ids:
        base = refs.safe_name(ref_id)
        seen[base] = seen.get(base, 0) + 1
        names.append(base if seen[base] == 1 else f"{base}_{seen[base]}")
    return names


def reference_score(metrics: dict) -> float:
    """
    (fraction of the reconstruction aligned) x (fraction of the reference covered), the coverage
    being the whole assembly's when there is one: as plasmid-surgery-refs ranks candidates.
    """
    coverage = metrics.get('assembly_ref_fraction_covered', 'NA')
    if coverage in ('NA', None):
        coverage = metrics['ref_fraction_covered']
    return round(float(metrics['contig_fraction_aligned']) * float(coverage), 4)


def rank_references(results: list[dict]) -> None:
    """Set 'score' and 'rank' (1 = best) on each reference's metrics, in place."""
    for m in results:
        m['score'] = reference_score(m)
    order = sorted(range(len(results)), key=lambda i: (-results[i]['score'],
                                                        -float(results[i]['ref_fraction_covered_including_repeats']), i))
    for rank, i in enumerate(order, start=1):
        results[i]['rank'] = rank


def assessment_fieldnames() -> list:
    fields = assessment.METRICS_FIELDNAMES + context.CONTEXT_METRICS_FIELDNAMES + features.FEATURE_METRICS_FIELDNAMES
    return fields[:2] + ['rank', 'score'] + fields[2:]


def write_fasta(path, seq_id: str, seq: str, bases_per_line: int, description: str = '') -> None:
    with open(path, 'w') as f:
        f.write(f">{seq_id} {description}".rstrip() + "\n")
        for line in wrap(seq, bases_per_line):
            f.write(f"{line}\n")


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
    os.makedirs(args.outdir, exist_ok=True)

    try:
        if args.tmpdir:
            os.makedirs(args.tmpdir, exist_ok=True)
            return run(args, args.tmpdir)
        with tempfile.TemporaryDirectory(prefix='plasmid-surgery_') as tmpdir:
            return run(args, tmpdir)
    except (alignment.AlignmentError, refs.RefSearchError) as e:
        logger.error(str(e))
        return 1


def load_references(paths) -> list[dict]:
    """Every sequence in the given FASTA files, each a separate reference."""
    references = []
    for path in paths:
        records = parsers.parse_fasta(path)
        if not records:
            logger.warning(f"No sequences in ref plasmid file: {path}")
        references += records
    return references


def find_references(args, tmpdir, outdir_files: dict) -> list[dict]:
    """References from the plasmid database (plasmid-surgery-refs' search); writes candidates.tsv."""
    search_tmp = Path(tmpdir) / 'reference_search'
    search_tmp.mkdir(parents=True, exist_ok=True)
    ranked, sequences = refs.find_refs(
        args.plasmid_reconstruction, args.db_sketch, args.db_fasta, search_tmp,
        db_fasta_index=args.db_fasta_index, assembly_fasta=args.assembly,
        num_candidates=args.num_candidates, threads=args.threads,
    )
    refs.write_candidates(outdir_files['candidates'], ranked, args.top)
    if ranked and ranked[0]['reconstruction_fraction_aligned'] < refs.WEAK_MATCH_FRACTION:
        logger.warning(f"The best candidate reference ({ranked[0]['ref_id']}) accounts for only "
                       f"{ranked[0]['reconstruction_fraction_aligned']:.0%} of the reconstruction: the database may not "
                       "contain a closely related plasmid")
    return [{'seq_id': r['ref_id'], 'seq': sequences[r['ref_id']], 'description': r['description']}
            for r in ranked[:args.top]]


def run(args, tmpdir) -> int:
    """
    Assess the reconstruction against each reference, rank the references, and write outputs.

    :param args: Parsed command-line arguments
    :param tmpdir: Directory for intermediate files
    :return: Exit code
    """
    outdir = Path(args.outdir)
    top_level = output_paths(outdir, args.prefix, TOP_LEVEL_FILENAMES)
    plasmid_reconstruction_contigs = deduplicate_contigs(parsers.parse_fasta(args.plasmid_reconstruction))

    if args.sample_id:
        sample_id = args.sample_id
    elif args.prefix:
        sample_id = args.prefix
    else:
        plasmid_reconstruction_basename = os.path.basename(args.plasmid_reconstruction)
        sample_id = plasmid_reconstruction_basename.split('.')[0]

    # Inputs shared by every reference, read once
    assembly_contigs = parsers.parse_fasta(args.assembly) if args.assembly else None
    bins = context.bin_labels(args.bins) if args.bins else None
    ref_feature_input = features.read_features(args.ref_features) if args.ref_features else None
    sample_feature_input = features.read_features(args.features) if args.features else None

    if args.ref_plasmid:
        references = load_references(args.ref_plasmid)
    else:
        references = find_references(args, tmpdir, top_level)
    if not references:
        logger.error("No references to assess against")
        with open(top_level['assessment'], 'w') as f:
            f.write('\t'.join(assessment_fieldnames()) + '\n')
        return 1

    results = []
    per_reference_outputs = []
    for ref_plasmid, dirname in zip(references, unique_dirnames([r['seq_id'] for r in references])):
        outputs = output_paths(outdir / dirname, args.prefix)
        os.makedirs(outdir / dirname, exist_ok=True)
        ref_tmpdir = Path(tmpdir) / dirname
        ref_tmpdir.mkdir(parents=True, exist_ok=True)
        write_fasta(outputs['reference'], ref_plasmid['seq_id'], ref_plasmid['seq'], args.output_bases_per_line,
                    ref_plasmid.get('description', ''))
        results.append(assess_reference(args, plasmid_reconstruction_contigs, ref_plasmid, outputs, ref_tmpdir, sample_id,
                                        assembly_contigs, bins, ref_feature_input, sample_feature_input))
        per_reference_outputs.append(outputs)

    rank_references(results)
    fieldnames = assessment_fieldnames()
    for metrics, outputs in zip(results, per_reference_outputs):
        assessment.write_metrics(outputs['assessment'], metrics, fieldnames)
    with open(top_level['assessment'], 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore', restval='NA')
        writer.writeheader()
        writer.writerows(sorted(results, key=lambda m: m['rank']))
    if len(results) > 1:
        best = min(results, key=lambda m: m['rank'])
        logger.info(f"Best reference: {best['ref_id']} (score {best['score']})")
    if args.report:
        from plasmid_surgery.report import write_report
        logger.info(f"Wrote {write_report(outdir)}")
    return 0


def assess_reference(args, plasmid_reconstruction_contigs, ref_plasmid, outputs, tmpdir, sample_id,
                     assembly_contigs, bins, ref_feature_input, sample_feature_input) -> dict:
    """
    Assess the reconstruction against one reference, writing that reference's outputs
    (except assessment.tsv, written once all references are ranked).

    :return: the assessment metrics
    """
    ref_plasmid_length = len(ref_plasmid['seq'])
    ref_plasmid_id = ref_plasmid['seq_id']

    logger.info(f"Aligning {len(plasmid_reconstruction_contigs)} reconstruction contig(s) vs. ref plasmid: {ref_plasmid_id}")
    try:
        recon_alignments = alignment.align_contigs(plasmid_reconstruction_contigs, ref_plasmid, tmpdir)
    except alignment.AlignmentError:
        raise
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

    context_result = None
    if assembly_contigs is not None:
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

    if ref_feature_input is not None:
        ref_feature_list, ref_feature_sequences = ref_feature_input
        mapped, unmapped = features.map_features(ref_feature_list, [ref_plasmid], ref_feature_sequences)
        if unmapped:
            log = logger.warning if not mapped else logger.info
            log(f"{unmapped} ref feature(s) are not on the ref plasmid ({ref_plasmid_id}) and were ignored")
        rows, feature_metrics = features.assess_ref_features(
            mapped, ref_plasmid_length, recon_alignments,
            assembly_alignments=context_result['placements'] if context_result else None,
            labels=context_result['locations'] if context_result else None,
            assembly_alternatives=context_result['alternatives'] if context_result else (),
            min_coverage=args.min_feature_coverage,
            min_deletion_length=args.min_indel_length,
        )
        features.write_rows(outputs['ref_features'], rows, features.REF_FEATURES_FIELDNAMES)
        metrics.update(feature_metrics)

    if sample_feature_input is not None:
        sample_feature_list, sample_feature_sequences = sample_feature_input
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

    logger.info(f"Assessment vs. {ref_plasmid_id}: " + json.dumps(metrics))
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

    return metrics


if __name__ == '__main__':
    sys.exit(main())
