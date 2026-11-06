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
import logging
import sys

from collections import Counter
from datetime import date
from pathlib import Path

from plasmid_surgery import report

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


def report_link(sample: dict, page_dir: Path) -> str:
    run = sample['run']
    path = run['outdir'] / f"{run['prefix']}report.html"
    if not path.exists():
        return report.esc(sample['sample_id'])
    try:
        href = Path(path).resolve().relative_to(page_dir.resolve())
    except ValueError:
        href = Path(path).resolve()
    return f'<a href="{report.esc(href.as_posix())}">{report.esc(sample["sample_id"])}</a>'


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
})();
"""


def render_cohort(cohort: dict, page_dir: Path, title: str = None) -> str:
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
<footer>Generated by plasmid-surgery {report.esc(__version__)} on {date.today().isoformat()}.</footer>
</main><script>{JS}</script></body></html>
"""


def write_cohort(run_dirs, output, ref_id=None, title=None) -> Path:
    cohort = gather([Path(d) for d in run_dirs], ref_id)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_cohort(cohort, output.parent, title))
    return output


def build_parser() -> argparse.ArgumentParser:
    from plasmid_surgery import __version__
    parser = argparse.ArgumentParser(
        prog='plasmid-surgery-cohort',
        description="Draw a self-contained HTML report of many samples' plasmid-surgery results against one shared "
                    "reference: coverage strips, and the reference's genes by sample.")
    parser.add_argument('runs', nargs='+', help="plasmid-surgery output directories (e.g. one per sample)")
    parser.add_argument('--ref', help="The reference to compare against (default: the one assessed in the most runs)")
    parser.add_argument('-o', '--output', default='cohort.html', help="Report path (default: cohort.html)")
    parser.add_argument('--title', help="Report title (default: 'plasmid-surgery cohort: <ref_id>')")
    parser.add_argument('--log-level', default='info', help="Log level (default: 'info')")
    parser.add_argument('--version', action='version', version=f"%(prog)s {__version__}")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    from plasmid_surgery.__main__ import setup_logging
    setup_logging(args.log_level)
    try:
        path = write_cohort(args.runs, args.output, args.ref, args.title)
    except CohortError as e:
        logger.error(str(e))
        return 1
    logger.info(f"Wrote {path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
