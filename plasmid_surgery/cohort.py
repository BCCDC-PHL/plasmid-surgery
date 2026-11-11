"""
A cohort report: many samples' plasmid-surgery results against one shared reference, in one page.

Takes several run output directories (typically one per sample) and, for the reference they share,
shows per sample: the headline metrics, a coverage strip along the reference (covered; repeat gap;
missing from the reconstruction but in a 'missed' contig, i.e. lost by binning; only in a shared or
repeat contig; missing from the assembly), and the reference's genes by status, plus genes the
samples have that the reference doesn't. Rows can be sorted together, and every element has a
tooltip. Like the per-run report it's self-contained: plain SVG and a little inline JavaScript.
"""
from __future__ import annotations

import argparse
import base64
import json
import logging
import os
import sys
import tempfile

from collections import Counter
from datetime import date
from pathlib import Path

from plasmid_surgery import report
import plasmid_surgery.alignment as alignment
import plasmid_surgery.compare as compare
import plasmid_surgery.core as core
import plasmid_surgery.parsers as parsers

logger = logging.getLogger(__name__)

STRIP_WIDTH = 900
STRIP_HEIGHT = 12

# Per-position categories along the reference, in drawing order
COVERED, REPEAT_GAP, MISSED, OTHER_ASSEMBLY, ABSENT = range(5)
CATEGORY_NAMES = {
    COVERED: 'covered by the reconstruction',
    REPEAT_GAP: 'repeat gap (inferred)',
    MISSED: 'not in the reconstruction; in a missed contig (binning)',
    OTHER_ASSEMBLY: 'not in the reconstruction; only in a shared or repeat contig',
    ABSENT: 'not in the reconstruction or the assembly',
}
CATEGORY_COLOURS = {
    COVERED: '#2f6db3',
    REPEAT_GAP: 'url(#cohort-hatch-repeat)',
    MISSED: '#d1663a',
    OTHER_ASSEMBLY: '#c9a0dc',
    ABSENT: '#c0392b',
}
LEGEND_COLOURS = {COVERED: '#2f6db3', REPEAT_GAP: report.GAP_COLOURS['repeat'], MISSED: '#d1663a',
                  OTHER_ASSEMBLY: '#c9a0dc', ABSENT: '#c0392b'}


class CohortError(Exception):
    pass


# --- gathering ---

def gather(run_dirs, ref_id=None) -> dict:
    """
    Each run's results against one reference.

    :param run_dirs: plasmid-surgery output directories
    :param ref_id: the shared reference; if None, the one assessed in the most runs
    :return: dict with 'ref_id', 'ref_length', 'samples' (one per run with that reference, in input
             order: run, reference results), and 'skipped' (run directory, reason)
    """
    runs = []
    skipped = []
    for d in run_dirs:
        try:
            runs.append(report.load_run(d))
        except FileNotFoundError as e:
            skipped.append((str(d), str(e)))
    if not runs:
        raise CohortError("No plasmid-surgery runs found")
    counts = Counter(ref['row']['ref_id'] for run in runs for ref in run['references'])
    if ref_id is None:
        if not counts:
            raise CohortError("The runs have no references")
        top = counts.most_common()
        if len(top) > 1 and top[0][1] == top[1][1]:
            tied = ', '.join(r for r, n in top if n == top[0][1])
            raise CohortError(f"Several references are assessed in the same number of runs ({tied}); choose one with --ref")
        ref_id = top[0][0]
    lengths = Counter(ref['length'] for run in runs for ref in run['references'] if ref['row']['ref_id'] == ref_id)
    if not lengths:
        raise CohortError(f"No run has results against {ref_id}")
    ref_length = lengths.most_common(1)[0][0]
    samples = []
    for run in runs:
        ref = next((r for r in run['references'] if r['row']['ref_id'] == ref_id), None)
        if ref is None:
            skipped.append((str(run['outdir']), f"not assessed against {ref_id}"))
        elif ref['length'] != ref_length:
            skipped.append((str(run['outdir']), f"its {ref_id} is {ref['length']:,} bp, not {ref_length:,} bp"))
        else:
            samples.append({'run': run, 'ref': ref, 'sample_id': ref['row']['sample_id']})
    for d, reason in skipped:
        logger.warning(f"Skipped {d}: {reason}")
    return {'ref_id': ref_id, 'ref_length': ref_length, 'samples': samples, 'skipped': skipped}


def categories(ref: dict, length: int) -> bytearray:
    """Per reference position, the category (COVERED, ...) for one sample."""
    cat = bytearray(length)   # COVERED
    sources = ref['contig_sources'] or []
    gap_mask = bytearray(length)
    for g in ref['gaps']:
        for a, b in report.wrap(g['start'], g['end'], length):
            code = REPEAT_GAP if g['type'] == 'repeat' else ABSENT
            cat[a:b] = bytes([code]) * (b - a)
            if g['type'] == 'unique':
                gap_mask[a:b] = b'\x01' * (b - a)
    if ref['contig_sources'] is not None:
        # Unique gaps covered by other assembly contigs: missed contigs take precedence
        for wanted, code in (({'shared', 'repeat'}, OTHER_ASSEMBLY), ({'missed'}, MISSED)):
            for r in sources:
                if r['classification'] not in wanted:
                    continue
                start, end = report.unwrap(r['ref_start'], r['ref_end'], length)
                for a, b in report.wrap(start, end, length):
                    for i in range(a, b):
                        if gap_mask[i]:
                            cat[i] = code
    return cat


def runs_of(values: bytearray) -> list[tuple[int, int, int]]:
    """(start, end, value) runs of equal values."""
    out = []
    start = 0
    for i in range(1, len(values) + 1):
        if i == len(values) or values[i] != values[start]:
            out.append((start, i, values[start]))
            start = i
    return out


# --- drawing ---

def strip_svg(sample: dict, length: int, px_per_bp: float) -> str:
    """One sample's coverage strip, with insertion ticks."""
    ref = sample['ref']
    s = report.Svg()
    top = 5
    for start, end, code in runs_of(categories(ref, length)):
        s.rect(start * px_per_bp, top, (end - start) * px_per_bp, STRIP_HEIGHT, CATEGORY_COLOURS[code],
               title=f"{sample['sample_id']}: {start:,}–{end:,} ({end - start:,} bp): {CATEGORY_NAMES[code]}")
    for e in ref['events']:
        x = int(e['ref_start']) * px_per_bp
        if e['event_type'] == 'insertion':
            s.polygon([(x - 4, 0), (x + 4, 0), (x, top + 1)], report.INK,
                      title=f"{sample['sample_id']}: insertion of {int(e['length']):,} bp at {int(e['ref_start']):,} ({e['contig_id']})")
    width = length * px_per_bp
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{top + STRIP_HEIGHT + 2}" '
            f'class="strip">{PATTERNS}{"".join(s.parts)}</svg>')


PATTERNS = (
    '<defs><pattern id="cohort-hatch-repeat" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
    f'<rect width="6" height="6" fill="{report.GAP_COLOURS["repeat"]}" opacity="0.35"/>'
    f'<line x1="0" y1="0" x2="0" y2="6" stroke="{report.GAP_COLOURS["repeat"]}" stroke-width="3"/></pattern></defs>'
)


def axis_svg(length: int, px_per_bp: float, ref_features) -> str:
    """The reference's ruler, with its features (names above them) if annotated."""
    s = report.Svg()
    y = 4
    feats = sorted(ref_features or [], key=lambda f: int(f['ref_start']))
    if feats:
        lanes = report.assign_lanes(
            [report.labelled_span(int(f['ref_start']), int(f['ref_end']), f['name'], px_per_bp) for f in feats], length)
        for f, lane in zip(feats, lanes):
            a, b = int(f['ref_start']), int(f['ref_end'])
            yy = y + 10 + lane * report.FEATURE_LANE
            s.rect(a * px_per_bp, yy, (b - a) * px_per_bp, 8, '#555555',
                   title=f"{f['name']}{': ' + f['product'] if f['product'] and f['product'] != f['name'] else ''} ({a:,}–{b:,})")
            s.text(a * px_per_bp, yy - 2, f['name'], size=9)
        y += (max(lanes) + 1) * report.FEATURE_LANE + 6
    step = report.nice_tick(length, px_per_bp)
    s.line(0, y + 14, length * px_per_bp, y + 14)
    for t in range(0, length + 1, step):
        s.line(t * px_per_bp, y + 10, t * px_per_bp, y + 14)
        s.text(min(max(t * px_per_bp, 12), length * px_per_bp - 12), y + 8,
               f"{t / 1000:g} kb" if step >= 1000 else f"{t:,}", size=9, fill=report.MUTED, anchor='middle')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{length * px_per_bp:.0f}" height="{y + 18}" '
            f'class="axis">{"".join(s.parts)}</svg>')


def report_path(sample: dict):
    """The sample's own report.html, if there is one."""
    run = sample['run']
    path = run['outdir'] / f"{run['prefix']}report.html"
    return path if path.exists() else None


def report_href(path: Path, page_dir: Path) -> str:
    """A link from the cohort page to a report: always relative (going up with ../ if needed), so the
    cohort page and the run directories can be moved or copied together."""
    return Path(os.path.relpath(Path(path).resolve(), Path(page_dir).resolve())).as_posix()


def report_link(sample: dict, page_dir: Path) -> str:
    path = report_path(sample)
    if path is None:
        return report.esc(sample['sample_id'])
    href = report.esc(report_href(path, page_dir))
    return f'<a href="{href}" data-report="{href}" target="_blank">{report.esc(sample["sample_id"])}</a>'


def sort_attrs(sample: dict, index: int) -> str:
    row = sample['ref']['row']
    return (f'data-order="{index}" data-sample="{report.esc(sample["sample_id"])}" '
            f'data-score="{report.esc(row.get("score", "0"))}" data-covered="{report.esc(row["ref_fraction_covered"])}"')


def summary_table(cohort: dict, page_dir: Path) -> str:
    samples = cohort['samples']
    has_assembly = any(s['ref']['row'].get('assembly_ref_fraction_covered', 'NA') != 'NA' for s in samples)
    has_features = any(s['ref']['row'].get('num_ref_features', 'NA') != 'NA' for s in samples)
    head = ['Sample', 'Rank', 'Score', 'Covered', 'Covered incl. repeats', 'Unique gaps', 'Reconstruction aligned',
            'Insertions', 'Deletions']
    if has_assembly:
        head += ['Missed (binning)', 'Unique gaps in assembly']
    if has_features:
        head += ['Ref features present']
    out = ['<table class="summary sortable"><thead><tr>'
           + ''.join(f'<th data-col="{i}">{h}</th>' for i, h in enumerate(head)) + '</tr></thead><tbody class="sortable-rows">']
    for i, s in enumerate(samples):
        r = s['ref']['row']
        cells = [(report_link(s, page_dir), s['sample_id']), (r.get('rank', ''), r.get('rank', '')),
                 (r.get('score', ''), r.get('score', '')),
                 (report.fraction(r['ref_fraction_covered']), r['ref_fraction_covered']),
                 (report.fraction(r['ref_fraction_covered_including_repeats']), r['ref_fraction_covered_including_repeats']),
                 (f"{report.bp(r['unique_gap_length'])} bp", r['unique_gap_length']),
                 (report.fraction(r['contig_fraction_aligned']), r['contig_fraction_aligned']),
                 (r['num_large_insertions'], r['num_large_insertions']), (r['num_large_deletions'], r['num_large_deletions'])]
        if has_assembly:
            cells += [(f"{report.bp(r['missed_ref_bases'])} bp", r['missed_ref_bases']),
                      (f"{report.bp(r['assembly_unique_gap_length'])} bp", r['assembly_unique_gap_length'])]
        if has_features:
            cells += [(f"{r['num_ref_features_present']} of {r['num_ref_features']}", r['num_ref_features_present'])]
        out.append(f'<tr {sort_attrs(s, i)}>' + ''.join(f'<td data-value="{report.esc(v)}">{c}</td>' for c, v in cells) + '</tr>')
    out.append('</tbody></table>')
    return ''.join(out)


def strips(cohort: dict, page_dir: Path) -> str:
    length = cohort['ref_length']
    px = STRIP_WIDTH / length
    ref_features = next((s['ref']['ref_features'] for s in cohort['samples'] if s['ref']['ref_features']), None)
    rows = [f'<div class="strip-row header"><div class="strip-label"></div><div>{axis_svg(length, px, ref_features)}</div></div>',
            '<div class="sortable-rows">']
    for i, s in enumerate(cohort['samples']):
        rows.append(f'<div class="strip-row" {sort_attrs(s, i)}><div class="strip-label">{report_link(s, page_dir)}</div>'
                    f'<div>{strip_svg(s, length, px)}</div></div>')
    rows.append('</div>')
    return ''.join(rows)


def gene_matrix(cohort: dict, page_dir: Path) -> str:
    """Reference genes (by position) x samples, coloured by status in the reconstruction."""
    samples = [s for s in cohort['samples'] if s['ref']['ref_features'] is not None]
    if not samples:
        return '<p class="sub">No reference annotations (<code>--ref-features</code>) in these runs.</p>'
    genes = {}
    for s in samples:
        for f in s['ref']['ref_features']:
            genes.setdefault((int(f['ref_start']), int(f['ref_end']), f['name']), f)
    columns = sorted(genes)
    head = ''.join(f'<th class="gene" title="{report.esc(n)} ({a:,}–{b:,})"><span>{report.esc(n)}</span></th>' for a, b, n in columns)
    out = [f'<table class="genes"><thead><tr><th>Sample</th>{head}</tr></thead><tbody class="sortable-rows">']
    for i, s in enumerate(cohort['samples']):
        if s['ref']['ref_features'] is None:
            continue
        by_key = {(int(f['ref_start']), int(f['ref_end']), f['name']): f for f in s['ref']['ref_features']}
        cells = []
        for key in columns:
            f = by_key.get(key)
            if f is None:
                cells.append('<td class="na">·</td>')
                continue
            status = f['reconstruction_status']
            binning = status in ('absent', 'partial') and f['assembly_status'] in ('present', 'inferred')
            colour = report.STATUS_COLOURS.get(status, report.MUTED)
            tip = (f"{s['sample_id']}: {f['name']} {status} ({float(f['reconstruction_fraction']):.0%})"
                   + (f"; assembly {f['assembly_status']}, found in {f['found_in']}" if f['assembly_status'] != 'NA' else ''))
            mark = '<span class="binning">●</span>' if binning else ''
            style = (f'background:repeating-linear-gradient(45deg,{colour} 0 3px,#fff 3px 6px)' if status == 'inferred'
                     else f'background:{colour}')
            cells.append(f'<td class="cell" style="{style}" title="{report.esc(tip)}">{mark}</td>')
        out.append(f'<tr {sort_attrs(s, i)}><td>{report_link(s, page_dir)}</td>{"".join(cells)}</tr>')
    out.append('</tbody></table>')
    return ''.join(out)


def novel_genes(cohort: dict, page_dir: Path) -> str:
    """Sample genes not on the reference (features.tsv, ref_status not_on_ref) x samples."""
    samples = [s for s in cohort['samples'] if s['ref']['features'] is not None]
    if not samples:
        return '<p class="sub">No sample annotations (<code>--features</code>) in these runs.</p>'
    names = sorted({f['name'] for s in samples for f in s['ref']['features'] if f['ref_status'] == 'not_on_ref'})
    if not names:
        return '<p class="sub">Every annotated sample gene on these contigs is on the reference.</p>'
    head = ''.join(f'<th class="gene"><span>{report.esc(n)}</span></th>' for n in names)
    out = [f'<table class="genes"><thead><tr><th>Sample</th>{head}</tr></thead><tbody class="sortable-rows">']
    for i, s in enumerate(cohort['samples']):
        if s['ref']['features'] is None:
            continue
        found = {}
        for f in s['ref']['features']:
            if f['ref_status'] == 'not_on_ref':
                found.setdefault(f['name'], []).append(f)
        cells = []
        for n in names:
            if n in found:
                tip = '; '.join(f"{f['contig_id']} ({f['contig_classification']}, bin {f['bin']})" for f in found[n])
                cells.append(f'<td class="cell" style="background:#5c6bc0" title="{report.esc(s["sample_id"] + ": " + n + " on " + tip)}"></td>')
            else:
                cells.append('<td class="na">·</td>')
        out.append(f'<tr {sort_attrs(s, i)}><td>{report_link(s, page_dir)}</td>{"".join(cells)}</tr>')
    out.append('</tbody></table>')
    return ''.join(out)


def legend(has_assembly: bool) -> str:
    codes = [COVERED, REPEAT_GAP] + ([MISSED, OTHER_ASSEMBLY, ABSENT] if has_assembly else [ABSENT])
    items = []
    for code in codes:
        label = CATEGORY_NAMES[code] if has_assembly or code != ABSENT else 'unique gap'
        kind = 'hatch' if code == REPEAT_GAP else 'swatch'
        items.append(f'<span class="legend-item"><span class="key {kind}" style="--c:{LEGEND_COLOURS[code]}"></span>{report.esc(label)}</span>')
    items.append(f'<span class="legend-item"><span class="key tri" style="--c:{report.INK}"></span>insertion</span>')
    return ''.join(items)


def gene_legend() -> str:
    items = [f'<span class="legend-item"><span class="key swatch" style="--c:{report.STATUS_COLOURS[k]}"></span>{k}</span>'
             for k in ('present', 'partial', 'absent')]
    items.append(f'<span class="legend-item"><span class="key hatch" style="--c:{report.STATUS_COLOURS["present"]}"></span>inferred</span>')
    items.append('<span class="legend-item"><span class="binning">●</span> in the assembly, not the reconstruction (binning)</span>')
    return ''.join(items)


CSS = report.CSS + """
.strip-row { display: grid; grid-template-columns: 200px auto; align-items: center; gap: 8px; }
.strip-row .strip-label { font-size: 12.5px; text-align: right; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.strip-row svg { display: block; }
.strip-row.header { align-items: end; }
.strips { overflow-x: auto; padding-bottom: 6px; }
table.genes th.gene { height: 110px; white-space: nowrap; vertical-align: bottom; padding: 0 2px; background: none; position: static; }
table.genes th.gene span { display: inline-block; transform: rotate(-60deg); transform-origin: left bottom; width: 16px; font-weight: 500; font-size: 12px; }
table.genes td { padding: 0; height: 20px; min-width: 20px; text-align: center; border: 1px solid #fff; }
table.genes td:first-child { padding: 2px 10px 2px 0; text-align: right; border: none; white-space: nowrap; }
table.genes td.na { color: #bbb; }
.binning { color: #fff; font-size: 11px; text-shadow: 0 0 2px #000; }
table.sortable th { cursor: pointer; } table.sortable th:hover { text-decoration: underline; }
.skipped { color: var(--muted); font-size: 12.5px; }
"""

JS = """
(function () {
  var current = 'covered';
  function key(el, by) {
    if (by === 'sample') { return el.getAttribute('data-sample'); }
    if (by === 'order') { return parseFloat(el.getAttribute('data-order')); }
    return -parseFloat(el.getAttribute('data-' + by)) || 0;
  }
  function sortAll(by) {
    current = by;
    document.querySelectorAll('.sortable-rows').forEach(function (container) {
      var rows = Array.prototype.slice.call(container.children);
      rows.sort(function (a, b) {
        var ka = key(a, by), kb = key(b, by);
        if (ka < kb) { return -1; } if (ka > kb) { return 1; }
        return parseFloat(a.getAttribute('data-order')) - parseFloat(b.getAttribute('data-order'));
      });
      rows.forEach(function (r) { container.appendChild(r); });
    });
  }
  document.querySelectorAll('input[name="sort"]').forEach(function (input) {
    input.addEventListener('change', function () { if (input.checked) { sortAll(input.value); } });
  });
  // Clicking a summary column sorts that table only, by that column (again to reverse)
  document.querySelectorAll('table.sortable th').forEach(function (th) {
    th.addEventListener('click', function () {
      var col = parseInt(th.getAttribute('data-col'), 10);
      var body = th.closest('table').querySelector('tbody');
      var rows = Array.prototype.slice.call(body.children);
      var dir = th.getAttribute('data-dir') === 'asc' ? -1 : 1;
      th.setAttribute('data-dir', dir === 1 ? 'asc' : 'desc');
      rows.sort(function (a, b) {
        var va = a.children[col].getAttribute('data-value'), vb = b.children[col].getAttribute('data-value');
        var na = parseFloat(va), nb = parseFloat(vb);
        if (!isNaN(na) && !isNaN(nb)) { return dir * (na - nb); }
        return dir * va.localeCompare(vb);
      });
      rows.forEach(function (r) { body.appendChild(r); });
    });
  });
  sortAll(current);
  // With --embed-reports: open a sample's report from the copy embedded in this page
  var store = document.getElementById('embedded-reports');
  if (store) {
    var reports = JSON.parse(store.textContent);
    document.querySelectorAll('a[data-report]').forEach(function (a) {
      var encoded = reports[a.getAttribute('data-report')];
      if (!encoded) { return; }
      a.addEventListener('click', function (event) {
        event.preventDefault();
        var bytes = Uint8Array.from(atob(encoded), function (c) { return c.charCodeAt(0); });
        var html = new TextDecoder('utf-8').decode(bytes);
        // A new blank window, written into (works from file:// pages, unlike blob: or data: URLs)
        var w = window.open('', '_blank');
        if (!w) { return; }
        w.document.open(); w.document.write(html); w.document.close();
      });
    });
  }
})();
"""


# --- comparing the reconstructions with each other ---

def comparison_members(cohort: dict) -> tuple[list[dict], list[str]]:
    """The shared reference and each sample's reconstruction (reconstruction.fa), and the samples without one."""
    members, missing = [], []
    for s in cohort['samples']:
        run = s['run']
        ref_fasta = s['ref']['dir'] / f"{run['prefix']}reference.fa"
        if ref_fasta.exists() and not members:
            members.append({'name': cohort['ref_id'], 'is_reference': True, 'contigs': parsers.parse_fasta(ref_fasta)})
    for s in cohort['samples']:
        run = s['run']
        path = run['outdir'] / f"{run['prefix']}reconstruction.fa"
        if path.exists():
            members.append({'name': s['sample_id'], 'is_reference': False, 'contigs': parsers.parse_fasta(path),
                            'sample': s})
        else:
            missing.append(s['sample_id'])
    return members, missing


def compare_reconstructions(cohort: dict, threads: int = 1, min_identity: float = compare.SNP_MIN_IDENTITY,
                            core_min_fraction: float = 1.0) -> dict:
    """All-vs-all comparison of the reconstructions and the reference (compare.all_vs_all), and core SNPs on the
    reference (core.core_snps), or why not."""
    members, missing = comparison_members(cohort)
    if missing:
        logger.warning(f"{len(missing)} run(s) have no reconstruction.fa (made before plasmid-surgery kept it; rerun to "
                       f"include them): {', '.join(missing[:5])}{' ...' if len(missing) > 5 else ''}")
    out = {'members': members, 'missing': missing, 'result': None, 'error': None, 'core': None}
    if sum(not m['is_reference'] for m in members) < 2:
        out['error'] = 'fewer than two reconstructions to compare'
        return out
    try:
        with tempfile.TemporaryDirectory(prefix='plasmid-surgery-cohort_') as tmp:
            out['result'] = compare.all_vs_all(members, tmp, threads, min_identity)
            reference = next((m for m in members if m['is_reference']), None)
            if reference is not None:
                samples = [m for m in members if not m['is_reference']]
                out['core'] = core.core_snps({'seq_id': reference['name'], 'seq': reference['contigs'][0]['seq']},
                                             samples, Path(tmp) / 'core', threads, core_min_fraction, min_identity)
    except alignment.AlignmentError as e:
        logger.error(str(e))
        out['error'] = str(e).splitlines()[0]
    return out


def containment_colour(v: float) -> str:
    """White (0) to dark blue (1)."""
    lo, hi = (255, 255, 255), (31, 78, 140)
    t = max(0.0, min(1.0, v))
    return '#' + ''.join(f'{round(a + (b - a) * t):02x}' for a, b in zip(lo, hi))


def _heatmap(comparison: dict, cell_fn) -> str:
    """Members against each other in clustered order; cell_fn(i, j) -> (fill, tooltip, text or '')."""
    result, members = comparison['result'], comparison['members']
    order = result['order']
    n = len(order)
    cell = max(6, min(22, 640 // max(n, 1)))
    names = result['names']
    shown = [names[i] + (' (reference)' if members[i]['is_reference'] else '') for i in range(len(names))]
    label_w = min(260, max(len(x) for x in shown) * 6.6 + 14)
    left, top = label_w, label_w
    width, height = left + n * cell + 10, top + n * cell + 10
    s = report.Svg()
    for r, i in enumerate(order):
        weight = 'bold' if members[i]['is_reference'] else 'normal'
        s.text(left - 6, top + r * cell + cell / 2 + 4, shown[i], size=10, anchor='end', weight=weight)
        x = left + r * cell + cell / 2 + 4
        s.text(x, top - 6, shown[i], size=10, weight=weight, extra=f'transform="rotate(-90 {x:.1f} {top - 6:.1f})"')
    for r, i in enumerate(order):
        for c, j in enumerate(order):
            if i == j:
                s.rect(left + c * cell, top + r * cell, cell - 1, cell - 1, report.FAINT)
                continue
            fill, tip, text = cell_fn(i, j)
            s.rect(left + c * cell, top + r * cell, cell - 1, cell - 1, fill, title=tip)
            if text and cell >= 16:
                dark = fill.startswith("#") and sum(int(fill[k:k + 2], 16) for k in (1, 3, 5)) < 480
                s.text(left + c * cell + (cell - 1) / 2, top + r * cell + cell / 2 + 3, text, size=8,
                       fill='white' if dark else report.INK, anchor='middle', extra='pointer-events="none"')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
            f'viewBox="0 0 {width:.0f} {height:.0f}" class="map">{"".join(s.parts)}</svg>')


def comparison_heatmap(comparison: dict, page_dir: Path) -> str:
    """The members against each other: cell (row, column) = how much of the row's sequence the column contains."""
    result = comparison['result']
    names = result['names']

    def cell(i, j):
        v = result['containment'][i][j]
        ident = result['identity'][i][j]
        tip = (f"{names[i]} in {names[j]}: {v:.1%} ({result['aligned'][i][j]:,} of {result['lengths'][i]:,} bp)\n"
               f"{names[j]} in {names[i]}: {result['containment'][j][i]:.1%} "
               f"({result['aligned'][j][i]:,} of {result['lengths'][j]:,} bp)"
               + (f"\nidentity {ident:.2%}" if ident is not None else ''))
        return containment_colour(v), tip, ''
    scale = ''.join(f'<span class="key" style="--c:{containment_colour(v)};border:1px solid #ccc"></span>'
                    for v in (0, 0.25, 0.5, 0.75, 1))
    return (f'<div class="legend"><span class="legend-item">contained: 0% {scale} 100%</span></div>'
            f'<div class="scroll">{_heatmap(comparison, cell)}</div>')


SNP_SCALE_MAX = 50     # SNPs (outside clusters) at which the SNP heatmap's colour is lightest


def snp_colour(snps: int) -> str:
    """Dark (0 SNPs) to light (SNP_SCALE_MAX or more), on a square-root scale."""
    return containment_colour(1 - min(snps, SNP_SCALE_MAX) ** 0.5 / SNP_SCALE_MAX ** 0.5)


def snp_heatmap(comparison: dict) -> str:
    """SNP distances (outside clusters), symmetric, in the same order; pairs compared over too little in grey."""
    result = comparison['result']
    names = result['names']

    def cell(i, j):
        f = compare.snp_fields(result, i, j)
        tip = (f"{names[i]} vs. {names[j]}: {f['snps']} SNPs ({f['snps_in_clusters']} in clusters), {f['indels']} indels\n"
               f"compared {f['compared_bases']:,} bp ({f['compared_fraction_a']:.0%} of {names[i]}, "
               f"{f['compared_fraction_b']:.0%} of {names[j]})"
               + (f"\n{f['snps_per_100kb']} SNPs per 100 kb outside clusters" if f['snps_per_100kb'] != 'NA' else ''))
        if not f['snp_distance_reported']:
            return '#eeeeee', tip + "\ntoo little compared for a distance", ''
        return snp_colour(f['snps_outside_clusters']), tip, str(f['snps_outside_clusters'])
    scale = ''.join(f'<span class="key" style="--c:{snp_colour(v)};border:1px solid #ccc"></span>'
                    for v in (0, 2, 10, 25, 50))
    return (f'<div class="legend"><span class="legend-item">SNPs: 0 {scale} {SNP_SCALE_MAX}+</span>'
            f'<span class="legend-item"><span class="key" style="--c:#eeeeee;border:1px solid #ccc"></span>'
            f'less than {compare.MIN_COMPARED_FRACTION:.0%} of the shorter compared</span></div>'
            f'<div class="scroll">{_heatmap(comparison, cell)}</div>')


def core_section(comparison: dict) -> str:
    """Core SNPs on the shared reference: a summary, the samples' calls, and a heatmap in the same order."""
    result = comparison.get('core')
    if result is None:
        return ''
    n, core_len = result['reference_length'], len(result['core_positions'])
    strict = result['min_fraction'] >= 1
    files = comparison.get('core_files') or []
    head = ('<h3>Core SNPs on the reference</h3><p class="sub">Every reconstruction placed on the shared reference and '
            'called at each of its positions (the same filters as above: close alignments only, one-to-one, repeats '
            'of the sample or of the reference left out, ends trimmed, deletions not called). The <b>core</b> is the '
            + ('positions called in every sample' if strict else
               f"positions called in at least {result['min_fraction']:.0%} of the samples; each pair is compared at the "
               'core positions both called')
            + f", so every pair is compared over the same positions. Core: <b>{core_len:,} bp</b> "
            f"({core_len / n:.0%} of the reference), {result['variable_sites']:,} variable sites"
            + (f"; alignment in <code>{report.esc(files[0].name)}</code>, pairs in <code>{report.esc(files[1].name)}</code>"
               if files else '') + '.</p>')
    rows = ''.join(f"<tr><td>{report.esc(r['sample'])}</td><td>{report.bp(r['called_positions'])}</td>"
                   f"<td>{report.fraction(r['called_fraction'])}</td>"
                   f"<td>{report.bp(r['core_without_it']) if r['core_without_it'] != 'NA' else 'NA'}</td></tr>"
                   for r in sorted(result['samples'], key=lambda r: r['called_fraction']))
    table = ('<p class="sub">How much of the reference each sample is called at, fewest first. '
             + ('"Core without it" is the core\'s size if that sample were left out: a fragmentary reconstruction '
                'shrinks the core for everyone.' if strict else '') + '</p>'
             '<div class="matrix-wrap"><table><thead><tr><th>Sample</th><th>Called</th><th>Of the reference</th>'
             '<th>Core without it</th></tr></thead><tbody>' + rows + '</tbody></table></div>')
    names = result['names']

    def cell(i, j):
        p = result['pairs'][(min(i, j), max(i, j))]
        tip = (f"{names[i]} vs. {names[j]}: {p['snps']} core SNPs ({p['snps_in_clusters']} in clusters)\n"
               f"compared at {p['compared_positions']:,} core positions")
        return snp_colour(p['snps_outside_clusters']), tip, str(p['snps_outside_clusters'])
    if core_len == 0:
        return head + table + '<p class="sub na">No core: no reference position is called in enough samples.</p>'
    return head + table + f'<div class="scroll">{_heatmap(comparison, cell)}</div>'


def comparison_section(comparison: dict, page_dir: Path, pairs_name: str = None) -> str:
    if comparison is None:
        return ''
    head = ('<h2>Reconstructions compared</h2><p class="sub">Every reconstruction (and the reference) aligned to '
            'every other, independent of the reference. Each cell is how much of the <b>row</b>\'s sequence the '
            '<b>column</b> contains (alignments of at least '
            f'{compare.MIN_ALIGNMENT} bp and {compare.MIN_IDENTITY:.0%} identity); hover for both directions and '
            'identity. A fragmentary reconstruction contained in a complete plasmid is dark in its row and lighter '
            'in its column. Ordered by clustering, so similar plasmids sit together'
            + (f'; all pairs in <code>{report.esc(pairs_name)}</code>' if pairs_name and comparison['result'] else '')
            + '.</p>')
    notes = ''
    if comparison['missing']:
        notes += (f'<p class="sub na">Not included (no <code>reconstruction.fa</code>: made before plasmid-surgery kept '
                  f'it; rerun to include): {report.esc(", ".join(comparison["missing"]))}</p>')
    if comparison['result'] is None:
        return head + notes + f'<p class="sub na">Not drawn: {report.esc(comparison["error"] or "")}</p>'
    snp_head = ('<h3>SNP distances</h3><p class="sub">SNPs between each pair, from one-to-one alignments: each '
                'stretch compared once, stretches aligning more than once between the pair (repeats) left out, and '
                f'differences within {compare.END_TRIM} bp of an alignment or contig end ignored. Cells show SNPs '
                f'outside clusters ({compare.CLUSTER_MIN} or more within {compare.CLUSTER_WINDOW:,} bp, e.g. a block '
                'from another lineage); hover for all SNPs, indels (counted separately) and how much was compared. '
                'Read every count with how much was compared: a pair is grey when less than '
                f'{compare.MIN_COMPARED_FRACTION:.0%} of the shorter of the two was.</p>')
    return (head + notes + comparison_heatmap(comparison, page_dir) + snp_head + snp_heatmap(comparison)
            + core_section(comparison))


def embedded_reports(cohort: dict, page_dir: Path) -> str:
    """The samples' own reports, base64-encoded in a JSON script block keyed by their link, for --embed-reports."""
    store = {}
    for s in cohort['samples']:
        path = report_path(s)
        if path is not None:
            store[report_href(path, page_dir)] = base64.b64encode(path.read_bytes()).decode('ascii')
    return f'<script type="application/json" id="embedded-reports">{json.dumps(store)}</script>'


def render_cohort(cohort: dict, page_dir: Path, title: str = None, embed_reports: bool = False,
                  comparison: dict = None, pairs_name: str = None) -> str:
    from plasmid_surgery import __version__
    samples = cohort['samples']
    has_assembly = any(s['ref']['contig_sources'] is not None for s in samples)
    title = title or f"plasmid-surgery cohort: {cohort['ref_id']}"
    skipped = ''
    if cohort['skipped']:
        skipped = ('<p class="skipped">Skipped: ' + '; '.join(f"{report.esc(d)} ({report.esc(r)})" for d, r in cohort['skipped'])
                   + '</p>')
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{report.esc(title)}</title><style>{CSS}</style></head>
<body><main>
<h1>{report.esc(title)}</h1>
<p class="sub">{len(samples)} sample{'s' if len(samples) != 1 else ''} against <code>{report.esc(cohort['ref_id'])}</code>
({cohort['ref_length']:,} bp). Hover over anything for details; sample names link to their own reports, where available.</p>
{skipped}
<div class="controls"><span>Order samples by:</span>
<label><input type="radio" name="sort" value="covered" checked> coverage</label>
<label><input type="radio" name="sort" value="score"> score</label>
<label><input type="radio" name="sort" value="sample"> sample</label>
<label><input type="radio" name="sort" value="order"> input order</label></div>
<h2>Samples</h2>
<p class="sub">Click a column to sort this table by it. "Covered" is how much of the reference each sample's
reconstruction covers. "Score" ranks references within one sample's run, and with the assembly it uses the
assembly's coverage, so it doesn't reflect sequence lost by binning: compare samples by coverage.</p>
<div class="scroll">{summary_table(cohort, page_dir)}</div>
<h2>Coverage of the reference</h2>
<div class="legend">{legend(has_assembly)}</div>
<div class="strips">{strips(cohort, page_dir)}</div>
<h2>Reference genes</h2>
<p class="sub">Each reference gene's status in each sample's reconstruction, in reference order.</p>
<div class="legend">{gene_legend()}</div>
<div class="scroll">{gene_matrix(cohort, page_dir)}</div>
<h2>Genes not on the reference</h2>
<p class="sub">Sample genes on the reconstruction (or, with the assembly, on contigs aligning to the reference) that
aren't on the reference itself.</p>
<div class="scroll">{novel_genes(cohort, page_dir)}</div>
{comparison_section(comparison, page_dir, pairs_name)}
<footer>Generated by plasmid-surgery {report.esc(__version__)} on {date.today().isoformat()}.</footer>
</main>{embedded_reports(cohort, page_dir) if embed_reports else ''}<script>{JS}</script></body></html>
"""


# --- choosing a reference ---

LIST_FIELDS = ['ref_id', 'ref_length', 'runs_assessed', 'rank_1', 'in_top', 'mean_score', 'mean_ref_fraction_covered']


def reference_table(run_dirs, top: int = 3) -> tuple[list[dict], int]:
    """
    How each reference fared across runs, to choose one for the cohort report. Reads only each run's
    top-level assessment.tsv.

    :param top: count the runs ranking a reference this high or better ('in_top')
    :return: (rows, number of runs read): one row per reference (LIST_FIELDS), the most often in the top first, then by mean score.
             'mean_score' is over all runs (a run that didn't assess the reference counts as 0);
             'mean_ref_fraction_covered' is over the runs that did. 'ref_length' is the most common.
    """
    runs = 0
    stats = {}
    for d in run_dirs:
        d = Path(d)
        try:
            rows = report.read_tsv(d / f"{report.find_prefix(d)}assessment.tsv")
        except FileNotFoundError as e:
            logger.warning(f"Skipped {d}: {e}")
            continue
        runs += 1
        for row in rows:
            s = stats.setdefault(row['ref_id'], {'lengths': Counter(), 'n': 0, 'rank_1': 0, 'in_top': 0,
                                                 'score': 0.0, 'covered': 0.0})
            rank = int(row['rank'])
            s['lengths'][row['ref_length']] += 1
            s['n'] += 1
            s['rank_1'] += rank == 1
            s['in_top'] += rank <= top
            s['score'] += float(row['score'])
            s['covered'] += float(row['ref_fraction_covered'])
    if not runs:
        raise CohortError("No plasmid-surgery runs found")
    table = [{'ref_id': ref_id, 'ref_length': s['lengths'].most_common(1)[0][0],
              'runs_assessed': s['n'], 'rank_1': s['rank_1'], 'in_top': s['in_top'],
              'mean_score': round(s['score'] / runs, 4), 'mean_ref_fraction_covered': round(s['covered'] / s['n'], 4)}
             for ref_id, s in stats.items()]
    table.sort(key=lambda r: (-r['in_top'], -r['mean_score'], r['ref_id']))
    return table, runs


def write_reference_table(run_dirs, out, top: int = 3) -> None:
    """Write reference_table() as a TSV to ``out`` (runs_assessed as 'n/total', in_top as 'in_top_<top>')."""
    table, runs = reference_table(run_dirs, top)
    header = [f'in_top_{top}' if f == 'in_top' else f for f in LIST_FIELDS]
    out.write('\t'.join(header) + '\n')
    for r in table:
        values = {**r, 'runs_assessed': f"{r['runs_assessed']}/{runs}"}
        out.write('\t'.join(str(values[f]) for f in LIST_FIELDS) + '\n')


def write_cohort(run_dirs, output, ref_id=None, title=None, embed_reports=False, compare_reconstructions_=True,
                 threads: int = 1, min_identity: float = compare.SNP_MIN_IDENTITY, core_min_fraction: float = 1.0) -> Path:
    """Write the cohort report; with the comparison, also <output stem>_pairs.tsv beside it, and the core SNP
    outputs (<output stem>_core.fasta, _core_pairs.tsv, _core_samples.tsv)."""
    cohort = gather([Path(d) for d in run_dirs], ref_id)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    comparison, pairs = None, None
    if compare_reconstructions_:
        comparison = compare_reconstructions(cohort, threads, min_identity, core_min_fraction)
        if comparison['result'] is not None:
            pairs = output.with_name(f"{output.stem}_pairs.tsv")
            compare.write_pairs(pairs, comparison['result'])
        if comparison['core'] is not None:
            comparison['core_files'] = [output.with_name(f"{output.stem}_core{x}") for x in
                                        ('.fasta', '_pairs.tsv', '_samples.tsv')]
            core.write_outputs(comparison['core'], *comparison['core_files'])
    output.write_text(render_cohort(cohort, output.parent, title, embed_reports, comparison,
                                    pairs.name if pairs else None))
    return output


def build_parser() -> argparse.ArgumentParser:
    from plasmid_surgery import __version__
    parser = argparse.ArgumentParser(
        prog='plasmid-surgery-cohort',
        description="Draw a self-contained HTML report of many samples' plasmid-surgery results against one shared "
                    "reference: coverage strips, and the reference's genes by sample.")
    parser.add_argument('runs', nargs='+', help="plasmid-surgery output directories (e.g. one per sample)")
    parser.add_argument('--ref', help="The reference to compare against (default: the one assessed in the most runs; see --list-refs to choose)")
    parser.add_argument('-o', '--output', default='cohort.html', help="Report path (default: cohort.html)")
    parser.add_argument('--list-refs', action='store_true',
                        help="Instead of drawing the report, print a TSV of how each reference fared across the runs "
                             "(how many runs assessed it, ranked it first, ranked it in the top --top; mean score and "
                             "coverage), to choose one for --ref")
    parser.add_argument('--top', type=int, default=3, help="With --list-refs: the rank cutoff for in_top_N (default: 3)")
    parser.add_argument('--title', help="Report title (default: 'plasmid-surgery cohort: <ref_id>')")
    parser.add_argument('--embed-reports', action='store_true', help="Embed each sample's own report.html in the cohort page, so its links work when the page is moved or shared on its own (otherwise links are relative paths to the run directories). Adds each report's size to the page")
    parser.add_argument('--no-compare', action='store_true', help="Don't compare the reconstructions with each other (by default, every run's reconstruction.fa and the reference are aligned all-vs-all with minimap2: a clustered matrix in the report, and <output stem>_pairs.tsv)")
    parser.add_argument('--snp-min-identity', type=float, default=compare.SNP_MIN_IDENTITY, help=f"SNP distances (pairwise and core): alignments less identical than this aren't compared, as they're likely a divergent homolog rather than the same locus (default: {compare.SNP_MIN_IDENTITY})")
    parser.add_argument('--core-min-fraction', type=float, default=1.0, help="Core SNPs: the fraction of samples a reference position must be called in to be in the core (default: 1.0, all). Below 1, each pair is compared at the core positions both called")
    parser.add_argument('-t', '--threads', type=int, default=1, help="Threads for minimap2 (the comparison) (default: 1)")
    parser.add_argument('--log-level', default='info', help="Log level (default: 'info')")
    parser.add_argument('--version', action='version', version=f"%(prog)s {__version__}")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    from plasmid_surgery.__main__ import setup_logging
    setup_logging(args.log_level)
    try:
        if args.list_refs:
            write_reference_table(args.runs, sys.stdout, args.top)
            return 0
        path = write_cohort(args.runs, args.output, args.ref, args.title, args.embed_reports,
                            not args.no_compare, args.threads, args.snp_min_identity, args.core_min_fraction)
    except CohortError as e:
        logger.error(str(e))
        return 1
    logger.info(f"Wrote {path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
