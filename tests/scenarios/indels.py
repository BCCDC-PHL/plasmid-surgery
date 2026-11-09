"""Insertions, deletions and identity: what counts as covered, what's an event, and how identity is measured."""

from __future__ import annotations

import random

from plasmid_surgery import parsers
from tests.scenarios.common import Outcome, Scenario, mutate, outputs, plasmid_map, random_seq, workspace

REF_LENGTH = 30000
DELETION = (6000, 8000)       # ref sequence missing from the sample
SMALL_DELETION = (12000, 12004)
INSERTION_AT, INSERTION_LENGTH = 22000, 3000
NUM_SNPS = 20


def build():
    rng = random.Random(47)
    ref = random_seq(rng, REF_LENGTH)
    x = mutate(rng, ref[:15000], NUM_SNPS // 2)
    y = mutate(rng, ref[15000:], NUM_SNPS // 2)
    x = x[:DELETION[0]] + x[DELETION[1]:SMALL_DELETION[0]] + x[SMALL_DELETION[1]:]
    i = INSERTION_AT - 15000
    insertion = random_seq(rng, INSERTION_LENGTH)
    y = y[:i] + insertion + y[i:]
    return ref, {'x': x, 'y': y}


def blast_identity(cigar_parsed) -> float:
    """matches / alignment columns, every inserted and deleted base counting as a column."""
    count = {op: sum(n for o, n in cigar_parsed if o == op) for op in
             ('sequence_match', 'sequence_mismatch', 'insertion_to_ref', 'deletion_to_ref')}
    return 100 * count['sequence_match'] / sum(count.values())


def run() -> Outcome:
    ws, tmp = workspace()
    ref, contigs = build()
    with tmp:
        ref_fa = ws.write_fasta('ref.fa', [('ref', ref)])
        recon_fa = ws.write_fasta('reconstruction.fa', contigs.items())
        out = outputs(ws.run('--plasmid-reconstruction', recon_fa, '--ref-plasmid', ref_fa,
                             '-o', ws.path('out'), '--tmpdir', ws.path('tmp')))
        y_sam = [a for a in parsers.parse_sam(ws.path('tmp') / 'ref' / 'reconstruction_vs_ref.sam') if a['contig_id'] == 'y' and a['is_mapped'] and not a['is_supplementary']]

    rows = {r['contig_id']: r for r in out['summary']}
    events = {e['event_type']: e for e in out['events']}
    a = out['assessment']
    facts = {
        'ref_length': REF_LENGTH,
        'deletion_start': DELETION[0], 'deletion_end': DELETION[1], 'deletion_length': DELETION[1] - DELETION[0],
        'small_deletion_length': SMALL_DELETION[1] - SMALL_DELETION[0], 'small_deletion_at': SMALL_DELETION[0],
        'insertion_at': INSERTION_AT, 'insertion_length': INSERTION_LENGTH,
        'num_snps': NUM_SNPS,
        'num_events': len(out['events']),
        'event_deletion_start': int(events['deletion']['ref_start']),
        'event_deletion_length': int(events['deletion']['length']),
        'event_insertion_at': int(events['insertion']['ref_start']),
        'event_insertion_length': int(events['insertion']['length']),
        'num_gaps': int(a['num_gaps']),
        'gap': f"{int(out['gaps'][0]['start']):,}–{int(out['gaps'][0]['end']):,}" if out['gaps'] else '-',
        'gap_start': int(out['gaps'][0]['start']), 'gap_end': int(out['gaps'][0]['end']),
        'covered': float(a['ref_fraction_covered']),
        'y_length': len(contigs['y']),
        'y_fraction_aligned': float(rows['y']['contig_fraction_aligned']),
        'y_identity': float(rows['y']['percent_identity']),
        'y_blast_identity': round(blast_identity(y_sam[0]['cigar_parsed']), 3),
        'x_identity': float(rows['x']['percent_identity']),
        'x_qc_pass': rows['x']['included_in_output'] == 'True' and rows['y']['included_in_output'] == 'True',
    }
    return Outcome(ws.commands, facts, {'map': plasmid_map(out, REF_LENGTH)})


SCENARIO = Scenario(
    slug='indels',
    title='Insertions, deletions and identity',
    summary="Large deletions are uncovered reference (gaps); insertions are sequence the reference lacks; "
            "both are reported as events, and identity counts each one once rather than per base.",
    text="""\
Two contigs, x and y, cover a {ref_length:,} bp reference with {num_snps} SNPs and three differences
from it:

* a {deletion_length:,} bp **deletion** (reference {deletion_start:,}–{deletion_end:,} is missing from x),
* a {small_deletion_length} bp deletion in x, at {small_deletion_at:,},
* a {insertion_length:,} bp **insertion** of sequence the reference doesn't have, in y at {insertion_at:,}.

Each contig still aligns in one piece, with the differences inside its alignment.

[[figure]]

**Events.** Insertions and deletions of at least `--min-indel-length` (default 50) bp within an
alignment are reported in `events.tsv`: here {num_events}, the deletion at {event_deletion_start:,}
({event_deletion_length:,} bp) and the insertion at {event_insertion_at:,} ({event_insertion_length:,} bp).
The {small_deletion_length} bp deletion is below the threshold.

**Coverage.** x's alignment spans the deleted region, but x doesn't have that sequence, so a deletion
of at least `--min-indel-length` bp counts as *uncovered*, and appears as a gap: {num_gaps} gap, at
{gap}. Shorter deletions count as covered, so that small indels don't break coverage into pieces. The
reconstruction covers {covered:.1%} of the reference.

**Insertions count as aligned contig sequence.** y's inserted bases are inside its alignment, so
`contig_fraction_aligned` is {y_fraction_aligned} for y, even though {insertion_length:,} of its
{y_length:,} bp aren't on the reference. The insertion is reported as an event rather than as
unaligned sequence. (Features are stricter: a gene inside the insertion is `not_on_ref`; see
[Feature annotations](features.md).)

**Identity** is *gap-compressed*: matches / (matches + mismatches + indel events), so each insertion
or deletion counts once, however long. y's identity is {y_identity}%. Counting every inserted base as
a mismatched column instead (BLAST-style identity) would give {y_blast_identity}%: one structural difference
would swamp the near-perfect identity of the rest of the contig, and with the default 80% identity
filter, a slightly larger insertion would discard the contig entirely. Indel content is filtered
separately (by the fraction of the alignment that is gaps). Both contigs pass QC here:
**{x_qc_pass}**.
""",
    run=run,
    test='tests/test_scenarios.py::test_indels',
)
