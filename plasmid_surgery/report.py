"""
A self-contained HTML report of a plasmid-surgery run: every reference's results in one page, drawn on
linear maps with a shared scale so they can be compared, plus an SVG map per reference for slides.

The report is drawn entirely from the files a run wrote (no re-analysis), as plain SVG with native
tooltips and a few lines of inline JavaScript (highlighting a contig across all panels, and switching
between a shared and a per-reference scale). It has no external dependencies and works offline.
"""
from __future__ import annotations

import argparse
import csv
import html
import logging
import sys

from datetime import date
from pathlib import Path

from plasmid_surgery.reads import read_bedgraph

logger = logging.getLogger(__name__)

# Reconstruction contig colours (cycled when there are more contigs than colours; labels and tooltips
# disambiguate). They avoid the colours that carry meaning elsewhere: gap red, missed orange, shared
# grey, repeat purple, and the feature-status greens and yellow.
PALETTE = ['#2f6db3', '#4fa3b8', '#b8963a', '#c97fb8', '#2e7d7d', '#5c6bc0', '#7a6a4f', '#1f4e8c',
           '#9e7bb5', '#6fa0c8', '#a0785a', '#3c9c9c', '#8c6d1f', '#b06090', '#4a5a8a', '#77a8a8']
INK, MUTED, FAINT, REF = '#222222', '#666666', '#dddddd', '#c8c8c8'
GAP_COLOURS = {'unique': '#c0392b', 'repeat': '#e69f00'}
CLASS_COLOURS = {'missed': '#d1663a', 'shared': '#9a9a9a', 'repeat': '#8a5cb8'}
STATUS_COLOURS = {'present': '#3a9a5b', 'inferred': '#9ccc9c', 'partial': '#e6b800', 'absent': '#c0392b'}
REF_STATUS_COLOURS = {'on_ref': '#3a9a5b', 'partial': '#e6b800', 'not_on_ref': '#c0392b'}
DEPTH_COLOUR, LOW_DEPTH_COLOUR = '#7f98b5', '#c0392b'
DEPTH_TRACK = 34       # height of the read depth track
CALL_LABELS = {'in_reads': 'in the reads', 'not_in_reads': 'not in the reads', 'partly_in_reads': 'partly in the reads'}

DRAW_WIDTH = 1000      # px for the longest reference (shared scale) or each reference (fitted)
LEFT, RIGHT = 120, 20  # margins: track labels on the left
LANE = 18              # vertical spacing of contig lanes
FEATURE_LANE = 26      # feature lanes: a label above each feature
CHAR_PX = 5.6          # approximate width of a character at the label font size
ARROW_HEIGHT = 12


# --- reading a run ---

def read_tsv(path: Path) -> list[dict]:
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def read_bed(path: Path) -> list[dict]:
    rows = []
    for line in path.read_text().splitlines():
        cols = line.split('\t')
        if len(cols) >= 4:
            rows.append({'start': int(cols[1]), 'end': int(cols[2]), 'name': cols[3],
                         'type': cols[3].rsplit('_', 1)[-1]})
    return rows


def find_prefix(outdir: Path) -> str:
    """The filename prefix of a run ('' if none), from its top-level assessment.tsv."""
    candidates = sorted(p.name for p in outdir.glob('*assessment.tsv') if p.is_file())
    if not candidates:
        raise FileNotFoundError(f"No assessment.tsv in {outdir}: is this a plasmid-surgery output directory?")
    if len(candidates) > 1:
        logger.warning(f"Several assessment files in {outdir}; using {candidates[0]}")
    return candidates[0][:-len('assessment.tsv')]


def load_run(outdir) -> dict:
    """
    Everything the report needs from a run's output directory.

    :return: dict with 'prefix', 'assessment' (top-level rows), 'candidates' (or None), and
             'references': per reference, its row, directory and parsed outputs, by rank
    """
    outdir = Path(outdir)
    prefix = find_prefix(outdir)
    name = lambda f: f"{prefix}{f}"  # noqa: E731
    top = read_tsv(outdir / name('assessment.tsv'))
    candidates_path = outdir / name('candidates.tsv')
    references = []
    for d in sorted(p for p in outdir.iterdir() if p.is_dir()):
        if not (d / name('assessment.tsv')).exists():
            continue
        [row] = read_tsv(d / name('assessment.tsv'))
        ref = {'row': row, 'dir': d, 'length': int(row['ref_length']),
               'summary': read_tsv(d / name('summary.tsv')),
               'gaps': read_bed(d / name('gaps.bed')),
               'events': read_tsv(d / name('events.tsv'))}
        for optional in ('contig_sources', 'ref_features', 'features', 'gap_support'):
            path = d / name(f'{optional}.tsv')
            ref[optional] = read_tsv(path) if path.exists() else None
        depth_path = d / name('reads_depth.bedgraph.gz')
        ref['reads_depth'] = read_bedgraph(depth_path) if depth_path.exists() else None
        references.append(ref)
    references.sort(key=lambda r: int(r['row'].get('rank') or 0))
    contig_depth_path = outdir / name('contig_depth.tsv')
    return {'outdir': outdir, 'prefix': prefix, 'assessment': top,
            'candidates': read_tsv(candidates_path) if candidates_path.exists() else None,
            'contig_depth': ({r['contig_id']: r for r in read_tsv(contig_depth_path)}
                             if contig_depth_path.exists() else None),
            'references': references}


def reconstruction_contigs(run: dict) -> list[tuple[str, int]]:
    """(contig_id, length) for every reconstruction contig, longest first."""
    lengths = {}
    for ref in run['references']:
        for r in ref['summary']:
            if r['contig_length'] not in ('', 'NA'):
                lengths[r['contig_id']] = int(r['contig_length'])
    return sorted(lengths.items(), key=lambda kv: (-kv[1], kv[0]))


def short_labels(ids) -> dict:
    """
    Labels for the maps: each ID minus the prefix all the IDs share, up to a separator (e.g.
    'SAMN123_17' -> '17'), so they fit inside contig arrows. Full IDs stay in tooltips and tables.
    """
    ids = sorted(set(ids))
    if len(ids) < 2:
        return {i: i for i in ids}
    first, last = ids[0], ids[-1]
    common = 0
    while common < min(len(first), len(last)) and first[common] == last[common]:
        common += 1
    cut = max((first.rfind(sep, 0, common) for sep in '_.-|:'), default=-1) + 1
    return {i: (i[cut:] or i) for i in ids}


def contig_colours(contigs) -> dict:
    return {cid: PALETTE[i % len(PALETTE)] for i, (cid, _) in enumerate(contigs)}


# --- drawing helpers ---

def esc(text) -> str:
    return html.escape(str(text), quote=True)


def fmt_bp(n) -> str:
    n = int(n)
    return f"{n / 1000:.1f} kb" if n >= 10000 else f"{n:,} bp"


def wrap(start, end, length):
    """Split an interval that may run past the origin into pieces within [0, length]."""
    if end <= length:
        return [(start, end)]
    return [(start, length), (0, end - length)]


def unwrap(start, end, length):
    """Coordinates from the outputs, where end < start means the interval spans the origin."""
    start, end = int(start), int(end)
    return (start, end + length) if end < start else (start, end)


def assign_lanes(intervals, length):
    """Greedy lane assignment so intervals in one lane don't overlap."""
    occupied, out = [], []
    for start, end in intervals:
        pieces = wrap(start, end, length)
        for lane, used in enumerate(occupied):
            if all(e <= s2 or s >= e2 for s, e in pieces for s2, e2 in used):
                used.extend(pieces)
                out.append(lane)
                break
        else:
            occupied.append(list(pieces))
            out.append(len(occupied) - 1)
    return out


def labelled_span(start, end, label, px_per_bp):
    """A feature's extent on the ref, widened to cover its label (drawn from its start), for lane assignment."""
    return start, max(end, start + int((len(label) * CHAR_PX + 6) / px_per_bp))


def nice_tick(length, px_per_bp):
    """A tick interval at least ~70 px apart."""
    for step in (100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000):
        if step * px_per_bp >= 70:
            return step
    return 1000000


class Svg:
    """A minimal SVG builder."""

    def __init__(self):
        self.parts = []

    def add(self, s):
        self.parts.append(s)

    def rect(self, x, y, w, h, fill, *, title=None, cls=None, data=None, stroke='none', extra=''):
        attrs = f'x="{x:.1f}" y="{y:.1f}" width="{max(w, 0.8):.1f}" height="{h:.1f}" fill="{fill}" stroke="{stroke}" {extra}'
        self._shape('rect', attrs, title, cls, data)

    def polygon(self, points, fill, *, title=None, cls=None, data=None, stroke='none', extra=''):
        pts = ' '.join(f'{x:.1f},{y:.1f}' for x, y in points)
        self._shape('polygon', f'points="{pts}" fill="{fill}" stroke="{stroke}" {extra}', title, cls, data)

    def text(self, x, y, s, *, size=10, fill=INK, anchor='start', weight='normal', extra=''):
        self.add(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" text-anchor="{anchor}" '
                 f'font-weight="{weight}" {extra}>{esc(s)}</text>')

    def line(self, x1, y1, x2, y2, stroke=MUTED, width=1):
        self.add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{width}"/>')

    def _shape(self, tag, attrs, title, cls, data):
        cls_attr = f' class="{cls}"' if cls else ''
        data_attr = ''.join(f' data-{k}="{esc(v)}"' for k, v in (data or {}).items())
        if title:
            self.add(f'<{tag} {attrs}{cls_attr}{data_attr}><title>{esc(title)}</title></{tag}>')
        else:
            self.add(f'<{tag} {attrs}{cls_attr}{data_attr}/>')


PATTERNS = (
    '<defs>'
    '<pattern id="hatch-repeat" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
    f'<rect width="6" height="6" fill="{GAP_COLOURS["repeat"]}" opacity="0.35"/>'
    f'<line x1="0" y1="0" x2="0" y2="6" stroke="{GAP_COLOURS["repeat"]}" stroke-width="3"/></pattern>'
    '<pattern id="hatch-inferred" width="5" height="5" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
    f'<rect width="5" height="5" fill="white"/><line x1="0" y1="0" x2="0" y2="5" stroke="{STATUS_COLOURS["present"]}" stroke-width="2"/></pattern>'
    '</defs>'
)


def arrow_points(x0, x1, y, reverse, height=ARROW_HEIGHT):
    tip = min((x1 - x0) * 0.4, 7)
    top, bottom, mid = y, y + height, y + height / 2
    if reverse:
        return [(x0, mid), (x0 + tip, top), (x1, top), (x1, bottom), (x0 + tip, bottom)]
    return [(x0, top), (x1 - tip, top), (x1, mid), (x1 - tip, bottom), (x0, bottom)]


# --- one reference's map ---

def depth_bins(segments, length, n_bins) -> list[float]:
    """Mean depth in each of n_bins equal bins along the reference, from (start, end, depth) segments."""
    sums = [0.0] * n_bins
    width = length / n_bins
    for start, end, d in segments:
        if not d:
            continue
        b0, b1 = int(start / width), min(int((end - 1) / width), n_bins - 1)
        for b in range(b0, b1 + 1):
            lo, hi = max(start, b * width), min(end, (b + 1) * width)
            if hi > lo:
                sums[b] += d * (hi - lo)
    return [v / width for v in sums]


def render_depth_track(s, ref: dict, x, px_per_bp: float, y: float) -> float:
    """Read depth along the reference: mean depth per ~2 px bar, scaled to 2x the median (higher bars are
    clipped, shown darker), and the stretches below the minimum read depth in red underneath."""
    length = ref['length']
    segments = ref['reads_depth']
    row = ref['row']
    try:
        min_depth = int(row.get('reads_min_depth', 'NA'))
    except ValueError:
        min_depth = 3
    n_bins = max(1, int(length * px_per_bp / 2))
    means = depth_bins(segments, length, n_bins)
    ordered = sorted(means)
    median = ordered[len(ordered) // 2] if ordered else 0
    top = max(2 * median, max(means, default=0) if median == 0 else 0, 1)
    s.text(LEFT - 8, y + DEPTH_TRACK / 2 + 4, 'read depth', size=10, fill=MUTED, anchor='end')
    s.text(LEFT - 8, y + 8, f"{top:,.0f}", size=8, fill=MUTED, anchor='end')
    s.line(x(0), y + DEPTH_TRACK, x(length), y + DEPTH_TRACK, FAINT)
    bin_bp = length / n_bins
    for i, m in enumerate(means):
        if m <= 0:
            continue
        h = min(m / top, 1) * DEPTH_TRACK
        a = i * bin_bp
        s.rect(x(a), y + DEPTH_TRACK - h, max(bin_bp * px_per_bp, 0.5), h,
               DEPTH_COLOUR if m <= top else '#4a6585',
               title=f"{a / 1000:.1f}–{(a + bin_bp) / 1000:.1f} kb: mean read depth {m:,.1f}",
               extra='shape-rendering="crispEdges"')
    low = []
    for start, end, d in segments:
        if d < min_depth:
            if low and low[-1][1] == start:
                low[-1] = (low[-1][0], end)
            else:
                low.append((start, end))
    for start, end in low:
        s.rect(x(start), y + DEPTH_TRACK + 1, max((end - start) * px_per_bp, 1), 3, LOW_DEPTH_COLOUR,
               title=f"{start:,}–{end:,} ({end - start:,} bp): read depth below {min_depth}")
    return y + DEPTH_TRACK + 12


def render_map(ref: dict, colours: dict, px_per_bp: float, *, standalone=False, title=None, labels=None) -> str:
    """The linear map of one reference, as an SVG string. labels: contig ID -> short label for the map."""
    labels = labels or {}
    length = ref['length']
    x = lambda pos: LEFT + pos * px_per_bp  # noqa: E731
    width = LEFT + length * px_per_bp + RIGHT
    s = Svg()
    y = 8
    if title:
        s.text(LEFT, y + 12, title, size=13, weight='bold')
        y += 24

    # Ruler
    step = nice_tick(length, px_per_bp)
    s.line(x(0), y + 14, x(length), y + 14, MUTED)
    for t in range(0, length + 1, step):
        s.line(x(t), y + 10, x(t), y + 14, MUTED)
        s.text(x(t), y + 8, f"{t / 1000:g} kb" if step >= 1000 else f"{t:,}", size=9, fill=MUTED, anchor='middle')
    y += 22

    # Reference, with gaps and insertions
    events = [e for e in ref['events'] if e['event_type'] == 'insertion']
    if events:
        y += 12
    s.text(LEFT - 8, y + 11, 'reference', size=10, fill=MUTED, anchor='end')
    s.rect(x(0), y, length * px_per_bp, 14, REF,
           title=f"{ref['row']['ref_id']}: {length:,} bp")
    support = {r['gap']: r for r in (ref['gap_support'] or []) if r['level'] == 'reconstruction'}
    for g in ref['gaps']:
        sup = support.get(g['name'])
        reads_note = (f"\n{CALL_LABELS.get(sup['call'], sup['call'])} ({float(sup['reads_fraction_covered']):.0%}, "
                      f"median depth {sup['reads_median_depth']})" if sup else '')
        for a, b in wrap(g['start'], g['end'], length):
            fill = 'url(#hatch-repeat)' if g['type'] == 'repeat' else GAP_COLOURS['unique']
            s.rect(x(a), y, (b - a) * px_per_bp, 14, fill,
                   title=f"{g['type']} gap {g['name']}: {g['start']:,}–{g['end'] % length if g['end'] > length else g['end']:,} "
                         f"({g['end'] - g['start']:,} bp){reads_note}")
    for e in events:
        px = x(int(e['ref_start']))
        s.polygon([(px - 5, y - 9), (px + 5, y - 9), (px, y - 1)], INK,
                  title=f"insertion of {int(e['length']):,} bp in {e['contig_id']} at {int(e['ref_start']):,}")
    y += 22

    # Read depth along the reference (with reads)
    if ref['reads_depth'] is not None:
        y = render_depth_track(s, ref, x, px_per_bp, y)

    # Reconstruction contigs
    rows = [r for r in ref['summary'] if r['is_mapped'] == 'True']
    kinds = {'placement': 0, 'excluded': 1, 'alternative': 2}
    items = []
    for r in rows:
        kind = 'alternative' if r['is_secondary'] == 'True' else 'placement'
        if r['included_in_output'] != 'True':
            if kind == 'alternative':
                continue   # weak secondaries: not alternative placements
            kind = 'excluded'
        items.append((kinds[kind], kind, r, unwrap(r['alignment_start'], r['alignment_end'], length)))
    items.sort(key=lambda t: (t[0], t[3][0]))
    lanes = assign_lanes([iv for *_, iv in items], length)
    s.text(LEFT - 8, y + 10, 'reconstruction', size=10, fill=MUTED, anchor='end')
    for (_, kind, r, (start, end)), lane in zip(items, lanes):
        cid = r['contig_id']
        c = colours.get(cid, MUTED)
        yy = y + lane * LANE
        tip = (f"{cid} ({int(r['contig_length']):,} bp){' ' + kind if kind != 'placement' else ''}\n"
               f"ref {int(r['alignment_start']):,}–{int(r['alignment_end']):,} ({'−' if r['reverse_complement'] == 'True' else '+'})"
               f"{', spans origin' if r['spans_origin'] == 'True' else ''}\n"
               f"identity {r['percent_identity']}%, MAPQ {r['mapq']}, placements {r['num_placements']}")
        for a, b in wrap(start, end, length):
            pts = arrow_points(x(a), x(b), yy, r['reverse_complement'] == 'True')
            if kind == 'placement':
                s.polygon(pts, c, title=tip, cls='contig', data={'contig': cid}, stroke='white', extra='stroke-width="0.5"')
            elif kind == 'alternative':
                s.polygon(pts, 'white', title=tip, cls='contig', data={'contig': cid}, stroke=c,
                          extra='stroke-width="1.5" stroke-dasharray="3,2"')
            else:
                s.polygon(pts, 'none', title=tip, cls='contig', data={'contig': cid}, stroke=MUTED,
                          extra='stroke-width="1" stroke-dasharray="2,2"')
        a, b = max(wrap(start, end, length), key=lambda p: p[1] - p[0])
        label = labels.get(cid, cid)
        if (b - a) * px_per_bp >= len(label) * 6 + 8 and kind == 'placement':
            s.text((x(a) + x(b)) / 2, yy + 9.5, label, size=9, fill='white', anchor='middle', extra='pointer-events="none"')
    y += (max(lanes, default=-1) + 1) * LANE + 8

    # Other assembly contigs (with --assembly)
    if ref['contig_sources'] is not None:
        others = [r for r in ref['contig_sources'] if r['classification'] != 'in_reconstruction']
        spans = [unwrap(r['ref_start'], r['ref_end'], length) for r in others]
        lanes = assign_lanes(spans, length)
        s.text(LEFT - 8, y + 10, 'other assembly', size=10, fill=MUTED, anchor='end')
        if not others:
            s.text(LEFT, y + 10, 'none align to this reference', size=9, fill=MUTED)
            lanes = [0]
        for r, (start, end), lane in zip(others, spans, lanes):
            c = CLASS_COLOURS.get(r['classification'], MUTED)
            yy = y + lane * LANE
            tip = (f"{r['contig_id']} ({int(r['contig_length']):,} bp): {r['classification']}, bin {r['bin']}\n"
                   f"{float(r['contig_fraction_on_ref']):.0%} on ref, depth {r['depth']} "
                   f"({r['depth_relative_to_reconstruction']}× reconstruction), placements {r['num_placements']}"
                   + (f"\nread depth {r['reads_depth']}" if r.get('reads_depth', 'NA') not in ('NA', '') else ''))
            for a, b in wrap(start, end, length):
                s.polygon(arrow_points(x(a), x(b), yy, False), c, title=tip, cls='assembly-contig',
                          data={'class': r['classification']}, stroke='white', extra='stroke-width="0.5"')
            a, b = max(wrap(start, end, length), key=lambda p: p[1] - p[0])
            label = labels.get(r['contig_id'], r['contig_id'])
            if (b - a) * px_per_bp >= len(label) * 6 + 8:
                s.text((x(a) + x(b)) / 2, yy + 9.5, label, size=9, fill='white', anchor='middle', extra='pointer-events="none"')
        y += (max(lanes, default=0) + 1) * LANE + 8

    # Reference features (with --ref-features)
    if ref['ref_features'] is not None:
        s.text(LEFT - 8, y + 10, 'ref features', size=10, fill=MUTED, anchor='end')
        if not ref['ref_features']:
            s.text(LEFT, y + 10, 'none annotated on this reference', size=9, fill=MUTED)
        feats = sorted(ref['ref_features'], key=lambda f: int(f['ref_start']))
        lanes = assign_lanes([labelled_span(int(f['ref_start']), int(f['ref_end']), f['name'], px_per_bp) for f in feats],
                             length)
        y += 10   # room for the labels
        for f, lane in zip(feats, lanes):
            status = f['reconstruction_status']
            fill = 'url(#hatch-inferred)' if status == 'inferred' else STATUS_COLOURS.get(status, MUTED)
            a, b = int(f['ref_start']), int(f['ref_end'])
            tip = (f"{f['name']}{': ' + f['product'] if f['product'] and f['product'] != f['name'] else ''}\n"
                   f"{f['category'] + ', ' if f['category'] else ''}{f['source']}\n"
                   f"reconstruction: {status} ({float(f['reconstruction_fraction']):.0%})"
                   + (f"\nassembly: {f['assembly_status']}, found in {f['found_in']}" if f['assembly_status'] != 'NA' else ''))
            yy = y + lane * FEATURE_LANE
            s.rect(x(a), yy, (b - a) * px_per_bp, 12, fill, title=tip, cls='feature',
                   data={'status': status}, stroke=INK, extra='stroke-width="0.4"')
            s.text(x(a), yy - 3, f['name'], size=9, fill=INK)
        y += max(max(lanes, default=-1) + 1, 1) * FEATURE_LANE

    # Sample features on the reference (with --features)
    not_on_ref = []
    if ref['features'] is not None:
        s.text(LEFT - 8, y + 10, 'sample features', size=10, fill=MUTED, anchor='end')
        placed = [f for f in ref['features'] if f['ref_start'] != 'NA']
        if not placed:
            s.text(LEFT, y + 10, 'none on this reference', size=9, fill=MUTED)
        not_on_ref = [f for f in ref['features'] if f['ref_start'] == 'NA']
        lanes = assign_lanes([labelled_span(int(f['ref_start']), int(f['ref_end']), f['name'], px_per_bp) for f in placed],
                             length)
        y += 10
        for f, lane in zip(placed, lanes):
            a, b = int(f['ref_start']), int(f['ref_end'])
            c = REF_STATUS_COLOURS.get(f['ref_status'], MUTED)
            cx, cy = (x(a) + x(b)) / 2, y + lane * FEATURE_LANE + 6
            s.text(x(a), cy - 8, f['name'], size=9, fill=INK)
            tip = (f"{f['name']} on {f['contig_id']} ({f['contig_classification']}, bin {f['bin']})\n"
                   f"{f['ref_status']} ({float(f['fraction_on_ref']):.0%} on ref), ref {a:,}–{b:,}")
            s.rect(x(a), cy - 2, (b - a) * px_per_bp, 4, c, cls='feature', data={'contig': f['contig_id']})
            s.polygon([(cx - 5, cy), (cx, cy - 5), (cx + 5, cy), (cx, cy + 5)], c, title=tip, cls='feature',
                      data={'contig': f['contig_id']}, stroke=INK, extra='stroke-width="0.4"')
        y += max(max(lanes, default=-1) + 1, 1) * FEATURE_LANE
        if not_on_ref:
            names = ', '.join(sorted({f['name'] for f in not_on_ref}))
            s.text(LEFT, y + 6, f"Not on this reference: {names}", size=9, fill=MUTED)
            y += 16

    height = y + 4
    body = ''.join(s.parts)
    style = (f'<style>text{{font-family:-apple-system,"Segoe UI",Helvetica,Arial,sans-serif}}</style>' if standalone else '')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
            f'viewBox="0 0 {width:.0f} {height:.0f}" class="map">{style}{PATTERNS}{body}</svg>')


# --- the reconstruction's own coordinates ---

EXPLAINED_COLOUR, ALTERNATIVE_COLOUR, UNEXPLAINED_COLOUR = '#4f6f8f', '#a9bccf', '#c0392b'
CONTIG_GAP_PX = 6


def explained_intervals(ref: dict) -> dict:
    """Per reconstruction contig, the (start, end, row) intervals of it that align to the reference: placements
    and alternative placements (from summary.tsv's contig_start/contig_end)."""
    out = {}
    for r in ref['summary']:
        if (r['is_mapped'] == 'True' and r.get('alignment_class') in ('placement', 'alternative')
                and r.get('contig_start', 'NA') not in ('NA', '')):
            out.setdefault(r['contig_id'], []).append((int(r['contig_start']), int(r['contig_end']), r))
    return out


def merged_length(intervals) -> int:
    total, end = 0, None
    for s, e in sorted(intervals):
        if end is None or s > end:
            total += e - s
            end = e
        elif e > end:
            total += e - end
            end = e
    return total


def uncovered(intervals, length) -> list[tuple[int, int]]:
    gaps, pos = [], 0
    for s, e in sorted(intervals):
        if s > pos:
            gaps.append((pos, s))
        pos = max(pos, e)
    if pos < length:
        gaps.append((pos, length))
    return gaps


def reconstruction_explained(run: dict, contigs) -> dict:
    """How much of the reconstruction each reference explains, and all of them together."""
    total = sum(n for _, n in contigs) or 1
    per_ref = []
    union = {cid: [] for cid, _ in contigs}
    for ref in run['references']:
        intervals = explained_intervals(ref)
        bases = 0
        for cid, n in contigs:
            ivs = [(s, e) for s, e, _ in intervals.get(cid, [])]
            bases += merged_length(ivs)
            union[cid] += ivs
        per_ref.append((ref, intervals, bases / total))
    unexplained = {cid: uncovered(union[cid], n) for cid, n in contigs}
    unexplained_bp = sum(e - s for gaps in unexplained.values() for s, e in gaps)
    return {'total': total, 'per_ref': per_ref, 'unexplained': unexplained,
            'any_fraction': 1 - unexplained_bp / total, 'unexplained_bp': unexplained_bp}


def render_reconstruction(run: dict, contigs, colours: dict, labels: dict, explained: dict) -> str:
    """
    The reconstruction on its own coordinates (contigs end to end, longest first): which parts of it each
    reference explains (aligns to), and which parts no reference does. The reference maps can't show the
    latter: sequence that aligns to no reference has no place on them.
    """
    total = explained['total']
    right = 56
    usable = DRAW_WIDTH - CONTIG_GAP_PX * max(len(contigs) - 1, 0)
    px = max(usable, 200) / total
    offsets, pos = {}, LEFT
    for cid, n in contigs:
        offsets[cid] = pos
        pos += n * px + CONTIG_GAP_PX
    width = pos - CONTIG_GAP_PX + right
    xs = lambda cid, p: offsets[cid] + p * px  # noqa: E731
    s = Svg()
    y = 8

    # Ruler: positions within each contig wide enough to carry ticks
    step = nice_tick(max(n for _, n in contigs), px) if contigs else 1
    for cid, n in contigs:
        if n * px < 60:
            continue
        s.line(offsets[cid], y + 14, offsets[cid] + n * px, y + 14, MUTED)
        for t in range(0, n + 1, step):
            s.line(xs(cid, t), y + 10, xs(cid, t), y + 14, MUTED)
            if t and (n - t) * px < 28:
                continue   # too close to the contig's end (the next contig's 0)
            s.text(xs(cid, t), y + 8, f"{t / 1000:g} kb" if step >= 1000 else f"{t:,}", size=9, fill=MUTED, anchor='middle')
    y += 22

    # Contigs
    s.text(LEFT - 8, y + 11, 'contigs', size=10, fill=MUTED, anchor='end')
    for cid, n in contigs:
        s.rect(offsets[cid], y, max(n * px, 1), 14, colours.get(cid, MUTED), title=f"{cid} ({n:,} bp)",
               cls='contig', data={'contig': cid})
        label = labels.get(cid, cid)
        if n * px >= len(label) * 6 + 8:
            s.text(offsets[cid] + n * px / 2, y + 10.5, label, size=9, fill='white', anchor='middle',
                   extra='pointer-events="none"')
    y += 24

    # One row per reference that explains any of it; the rest in one line
    none = [ref['row']['ref_id'] for ref, _, frac in explained['per_ref'] if frac == 0]
    for ref, intervals, frac in explained['per_ref']:
        if frac == 0:
            continue
        row = ref['row']
        name = f"#{row.get('rank', '')} {row['ref_id']}"
        s.text(LEFT - 8, y + 10, name if len(name) <= 20 else name[:19] + '…', size=10, fill=MUTED, anchor='end')
        s.line(LEFT, y + 6, pos - CONTIG_GAP_PX, y + 6, FAINT)
        for cid, n in contigs:
            for a, b, r in sorted(intervals.get(cid, []), key=lambda t: t[2].get('alignment_class') != 'alternative'):
                alt = r.get('alignment_class') == 'alternative'
                tip = (f"{cid}:{a + 1:,}–{b:,} ({b - a:,} bp) aligns to {row['ref_id']} "
                       f"{int(r['alignment_start']):,}–{int(r['alignment_end']):,} "
                       f"({'−' if r['reverse_complement'] == 'True' else '+'})\n"
                       f"identity {r['percent_identity']}%{', alternative placement (repeat)' if alt else ''}")
                s.rect(xs(cid, a), y, max((b - a) * px, 1), 12, ALTERNATIVE_COLOUR if alt else EXPLAINED_COLOUR,
                       title=tip, cls='explained', data={'contig': cid})
        s.text(pos - CONTIG_GAP_PX + 6, y + 10, f"{frac:.0%}", size=10, fill=INK)
        y += LANE

    # No reference
    s.text(LEFT - 8, y + 10, 'no reference', size=10, fill=MUTED, anchor='end')
    s.line(LEFT, y + 6, pos - CONTIG_GAP_PX, y + 6, FAINT)
    for cid, n in contigs:
        for a, b in explained['unexplained'][cid]:
            s.rect(xs(cid, a), y, max((b - a) * px, 1), 12, UNEXPLAINED_COLOUR,
                   title=f"{cid}:{a + 1:,}–{b:,} ({b - a:,} bp): aligns to none of the references",
                   cls='unexplained', data={'contig': cid})
    s.text(pos - CONTIG_GAP_PX + 6, y + 10, f"{1 - explained['any_fraction']:.0%}", size=10, fill=UNEXPLAINED_COLOUR)
    y += LANE + 4
    if none:
        listed = ', '.join(none[:8]) + (f' and {len(none) - 8} more' if len(none) > 8 else '')
        s.text(LEFT, y + 6, f"Explaining none of it: {listed}", size=9, fill=MUTED)
        y += 16

    # The sample's features on the reconstruction (with --features)
    feats = next((ref['features'] for ref in run['references'] if ref['features'] is not None), None)
    if feats is not None:
        placed = [f for f in feats if f['contig_id'] in offsets and f['contig_start'] not in ('NA', '')]
        seen, unique = set(), []
        for f in placed:
            key = (f['contig_id'], f['contig_start'], f['contig_end'], f['name'])
            if key not in seen:
                seen.add(key)
                unique.append(f)
        unique.sort(key=lambda f: xs(f['contig_id'], int(f['contig_start'])))
        s.text(LEFT - 8, y + 10 + 10, 'sample features', size=10, fill=MUTED, anchor='end')
        spans = [(xs(f['contig_id'], int(f['contig_start'])),
                  max(xs(f['contig_id'], int(f['contig_end'])), xs(f['contig_id'], int(f['contig_start'])) + len(f['name']) * CHAR_PX + 6))
                 for f in unique]
        lanes = assign_lanes([(int(a), int(b)) for a, b in spans], int(width) + 1)
        y += 10
        for f, lane in zip(unique, lanes):
            a, b = int(f['contig_start']), int(f['contig_end'])
            yy = y + lane * FEATURE_LANE
            outside = sum(max(0, min(b, e) - max(a, st)) for st, e in explained['unexplained'][f['contig_id']])
            off_refs = b > a and outside / (b - a) > 0.5
            tip = (f"{f['name']}{': ' + f['product'] if f.get('product') and f['product'] != f['name'] else ''}\n"
                   f"{f['contig_id']}:{a + 1:,}–{b:,}" + ("\nin a part no reference explains" if off_refs else ''))
            s.rect(xs(f['contig_id'], a), yy, max((b - a) * px, 1.5), 10, UNEXPLAINED_COLOUR if off_refs else INK,
                   title=tip, cls='feature', data={'contig': f['contig_id']})
            s.text(xs(f['contig_id'], a), yy - 3, f['name'], size=9, fill=INK)
        y += max(max(lanes, default=-1) + 1, 1) * FEATURE_LANE
        if not unique:
            s.text(LEFT, y, 'none on the reconstruction', size=9, fill=MUTED)
            y += 12

    height = y + 4
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
            f'viewBox="0 0 {width:.0f} {height:.0f}" class="map">{"".join(s.parts)}</svg>')


def reconstruction_summary(run: dict, explained: dict) -> str:
    if not explained['per_ref']:
        return ''
    best_ref, _, best = max(explained['per_ref'], key=lambda t: t[2])
    n_stretches = sum(len(g) for g in explained['unexplained'].values())
    text = (f"{explained['any_fraction']:.0%} of the reconstruction ({explained['total'] - explained['unexplained_bp']:,} of "
            f"{explained['total']:,} bp) aligns to at least one reference; the best single reference, "
            f"{best_ref['row']['ref_id']}, explains {best:.0%}.")
    if explained['unexplained_bp']:
        text += (f" The other {explained['unexplained_bp']:,} bp ({n_stretches} stretch{'es' if n_stretches != 1 else ''}) "
                 "align to none of them: to find what they are, search a plasmid database with a reference's "
                 "<code>unaligned.fa</code> (e.g. <code>plasmid-surgery-refs</code>).")
    if explained['any_fraction'] - best > 0.1:
        text += (" The references together explain clearly more than any one of them: the reconstruction may combine "
                 "sequence from more than one plasmid.")
    return text


def reconstruction_section(run, contigs, colours, labels, explained) -> str:
    if not any('contig_start' in r for ref in run['references'] for r in ref['summary']):
        return ('<p class="sub na">Not available: this run predates the contig coordinates in <code>summary.tsv</code>; '
                'rerun plasmid-surgery to draw it.</p>')
    return (f"""<p class="sub">The reconstruction's contigs end to end (longest first), with the parts of each that align to each
reference (dark; alternative placements of repeats lighter), and the parts that align to none (red). Percentages
are of the whole reconstruction. With <code>--features</code>, genes mostly in parts no reference explains are red.
{reconstruction_summary(run, explained)}</p>
<div class="scroll">{render_reconstruction(run, contigs, colours, labels, explained)}</div>""")


# --- the page ---

def fraction(v) -> str:
    try:
        return f"{float(v):.1%}"
    except (TypeError, ValueError):
        return 'NA'


def bp(v) -> str:
    try:
        return f"{int(v):,}"
    except (TypeError, ValueError):
        return 'NA'


def summary_table(run: dict) -> str:
    rows = sorted(run['assessment'], key=lambda r: int(r.get('rank') or 0))
    has_assembly = any(r.get('assembly_ref_fraction_covered', 'NA') != 'NA' for r in rows)
    has_features = any(r.get('num_ref_features', 'NA') != 'NA' for r in rows)
    has_reads = any(r.get('reads_ref_fraction_covered', 'NA') != 'NA' for r in rows)
    head = ['Rank', 'Reference', 'Length', 'Score', 'Covered', 'Covered incl. repeats', 'Unique gaps',
            'Reconstruction aligned']
    if has_assembly:
        head += ['Assembly covers', 'Missed (binning)', 'Unique gaps in assembly']
    if has_reads:
        head += ['In the reads', 'Read depth (median)', 'Unique gaps in the reads']
    if has_features:
        head += ['Ref features present']
    out = ['<table class="summary"><thead><tr>' + ''.join(f'<th>{h}</th>' for h in head) + '</tr></thead><tbody>']
    for r in rows:
        cells = [r.get('rank', ''), f'<a href="#ref-{esc(r["ref_id"])}">{esc(r["ref_id"])}</a>', bp(r['ref_length']),
                 r.get('score', ''), fraction(r['ref_fraction_covered']), fraction(r['ref_fraction_covered_including_repeats']),
                 f"{bp(r['unique_gap_length'])} bp ({r['num_unique_gaps']})", fraction(r['contig_fraction_aligned'])]
        if has_assembly:
            cells += [fraction(r['assembly_ref_fraction_covered']), f"{bp(r['missed_ref_bases'])} bp",
                      f"{bp(r['assembly_unique_gap_length'])} bp"]
        if has_reads:
            cells += [fraction(r['reads_ref_fraction_covered']), r.get('reads_median_depth', 'NA'),
                      f"{bp(r['unique_gap_in_reads_length'])} bp"]
        if has_features:
            cells += [f"{r['num_ref_features_present']} of {r['num_ref_features']}"
                      + (f" (+{r['num_ref_features_inferred']} inferred)" if r['num_ref_features_inferred'] not in ('0', 'NA') else '')]
        out.append('<tr>' + ''.join(f'<td>{c}</td>' for c in cells) + '</tr>')
    out.append('</tbody></table>')
    return ''.join(out)


def contig_matrix(run: dict, contigs, colours) -> str:
    """Rows: reconstruction contigs; columns: references; cells: where each contig is on each."""
    refs = run['references']
    depth = run.get('contig_depth')
    head = (['Contig', 'Length'] + (['Read depth'] if depth is not None else [])
            + [f'#{r["row"]["rank"]} {esc(r["row"]["ref_id"])}' for r in refs])
    out = ['<table class="matrix"><thead><tr>' + ''.join(f'<th>{h}</th>' for h in head) + '</tr></thead><tbody>']
    for cid, length in contigs:
        cells = [f'<span class="swatch" style="background:{colours[cid]}"></span>{esc(cid)}', f'{length:,}']
        if depth is not None:
            cells.append(esc(depth[cid]['reads_depth']) if cid in depth else '<span class="na">NA</span>')
        for ref in refs:
            rows = [r for r in ref['summary'] if r['contig_id'] == cid and r['is_mapped'] == 'True'
                    and r['is_secondary'] != 'True']
            alts = [r for r in ref['summary'] if r['contig_id'] == cid and r['is_secondary'] == 'True'
                    and r['included_in_output'] == 'True']
            if not rows:
                cells.append('<span class="na">not aligned</span>')
                continue
            frac = float(rows[0]['contig_fraction_aligned'])
            where = ', '.join(f"{int(r['alignment_start']) / 1000:.1f}–{int(r['alignment_end']) / 1000:.1f} kb"
                              f"{' (−)' if r['reverse_complement'] == 'True' else ''}" for r in rows)
            extra = f' <span class="na">+{len(alts)} alt.</span>' if alts else ''
            cells.append(f'<span title="{frac:.0%} of the contig aligned">{where}</span>{extra}'
                         + ('' if frac >= 0.99 else f' <span class="na">({frac:.0%})</span>'))
        out.append(f'<tr data-contig="{esc(cid)}">' + ''.join(f'<td>{c}</td>' for c in cells) + '</tr>')
    out.append('</tbody></table>')
    return ''.join(out)


def legend(run: dict) -> str:
    items = [('swatch', REF, 'reference'), ('swatch', GAP_COLOURS['unique'], 'unique gap'),
             ('hatch', GAP_COLOURS['repeat'], 'repeat gap'), ('tri', INK, 'insertion'),
             ('outline', MUTED, 'alternative placement (dashed) / not in scaffold (dotted)')]
    if any(r['contig_sources'] is not None for r in run['references']):
        items += [('swatch', CLASS_COLOURS['missed'], 'missed contig'), ('swatch', CLASS_COLOURS['shared'], 'shared contig'),
                  ('swatch', CLASS_COLOURS['repeat'], 'repeat contig')]
    if any(r['reads_depth'] is not None for r in run['references']):
        items += [('swatch', DEPTH_COLOUR, 'read depth'), ('swatch', LOW_DEPTH_COLOUR, 'below the minimum read depth')]
    if any(r['ref_features'] is not None for r in run['references']):
        items += [('swatch', STATUS_COLOURS[k], f'ref feature {k}') for k in ('present', 'partial', 'absent')]
        items += [('hatch', STATUS_COLOURS['present'], 'ref feature inferred')]
    if any(r['features'] is not None for r in run['references']):
        items += [('diamond', REF_STATUS_COLOURS['on_ref'], 'sample feature on ref'),
                  ('diamond', REF_STATUS_COLOURS['partial'], 'partly')]
    out = []
    for kind, colour, label in items:
        out.append(f'<span class="legend-item"><span class="key {kind}" style="--c:{colour}"></span>{esc(label)}</span>')
    return ''.join(out)


CSS = """
:root { --ink: #222; --muted: #666; --line: #ddd; --bg: #fff; --panel: #fafafa; }
* { box-sizing: border-box; }
body { margin: 0; font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; color: var(--ink);
       background: var(--bg); font-size: 14px; line-height: 1.45; }
main { max-width: 1240px; margin: 0 auto; padding: 24px 16px 64px; }
h1 { font-size: 22px; margin: 0 0 4px; } h2 { font-size: 17px; margin: 32px 0 8px; }
.sub { color: var(--muted); margin: 0 0 16px; }
table { border-collapse: collapse; font-size: 13px; margin: 8px 0; }
th, td { border-bottom: 1px solid var(--line); padding: 5px 10px; text-align: left; vertical-align: top; }
th { font-weight: 600; background: var(--panel); position: sticky; top: 0; }
td:not(:first-child) { font-variant-numeric: tabular-nums; }
.scroll { overflow-x: auto; max-width: 100%; }
.matrix-wrap { max-height: 360px; overflow: auto; border: 1px solid var(--line); }
.na { color: var(--muted); }
.swatch { display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 6px; vertical-align: -1px; }
.panel { border: 1px solid var(--line); border-radius: 6px; padding: 12px 12px 4px; margin: 14px 0; background: var(--bg); }
.panel h3 { font-size: 15px; margin: 0 0 2px; } .panel .stats { color: var(--muted); font-size: 12.5px; margin: 0 0 6px; }
.panel .scroll svg { display: block; }
.controls { display: flex; gap: 18px; align-items: center; flex-wrap: wrap; margin: 10px 0; color: var(--muted); }
.legend { display: flex; flex-wrap: wrap; gap: 6px 16px; font-size: 12px; color: var(--muted); margin: 6px 0 4px; }
.legend-item { display: inline-flex; align-items: center; gap: 6px; }
.key { display: inline-block; width: 16px; height: 10px; border-radius: 2px; background: var(--c); }
.key.hatch { background: repeating-linear-gradient(45deg, var(--c) 0 2px, transparent 2px 5px); border: 1px solid var(--c); }
.key.outline { background: none; border: 1.5px dashed var(--c); }
.key.tri { width: 0; height: 0; border-left: 6px solid transparent; border-right: 6px solid transparent;
           border-top: 9px solid var(--c); background: none; border-radius: 0; }
.key.diamond { width: 9px; height: 9px; transform: rotate(45deg); border-radius: 1px; }
body.highlighting .contig:not(.hl), body.highlighting tr[data-contig]:not(.hl) { opacity: 0.25; }
.contig.hl { stroke: #000 !important; stroke-width: 1.8 !important; }
tr[data-contig] { cursor: default; } tr.hl { background: #fff6dd; }
.scale-fit .map-fit, .scale-shared .map-shared { display: block; } .map-fit, .map-shared { display: none; }
footer { color: var(--muted); font-size: 12px; margin-top: 40px; }
code { font-size: 12.5px; }
"""

JS = """
(function () {
  function setHighlight(contig, on) {
    document.querySelectorAll('[data-contig]').forEach(function (el) {
      if (el.getAttribute('data-contig') === contig) { el.classList.toggle('hl', on); }
    });
    document.body.classList.toggle('highlighting', on);
  }
  document.querySelectorAll('[data-contig]').forEach(function (el) {
    el.addEventListener('mouseenter', function () { setHighlight(el.getAttribute('data-contig'), true); });
    el.addEventListener('mouseleave', function () { setHighlight(el.getAttribute('data-contig'), false); });
  });
  document.querySelectorAll('input[name="scale"]').forEach(function (input) {
    input.addEventListener('change', function () {
      document.body.classList.toggle('scale-fit', input.value === 'fit' && input.checked);
      document.body.classList.toggle('scale-shared', input.value === 'shared' && input.checked);
    });
  });
})();
"""


def panel_stats(row: dict) -> str:
    parts = [f"score {row.get('score', 'NA')}",
             f"covered {fraction(row['ref_fraction_covered'])} ({fraction(row['ref_fraction_covered_including_repeats'])} with repeats)",
             f"unique gaps {bp(row['unique_gap_length'])} bp",
             f"reconstruction aligned {fraction(row['contig_fraction_aligned'])}"]
    if row.get('num_large_insertions') not in (None, 'NA', '0') or row.get('num_large_deletions') not in (None, 'NA', '0'):
        parts.append(f"{row['num_large_insertions']} insertion(s), {row['num_large_deletions']} deletion(s)")
    if row.get('assembly_ref_fraction_covered', 'NA') != 'NA':
        parts.append(f"assembly covers {fraction(row['assembly_ref_fraction_covered'])}, "
                     f"{bp(row['missed_ref_bases'])} bp only in missed contigs")
    if row.get('reads_ref_fraction_covered', 'NA') != 'NA':
        parts.append(f"in the reads {fraction(row['reads_ref_fraction_covered'])} (median depth {row['reads_median_depth']}), "
                     f"{bp(row['unique_gap_in_reads_length'])} bp of unique gaps in the reads")
    return ' · '.join(parts)


def render_report(run: dict, title: str = None) -> tuple[str, dict]:
    """
    :return: (the HTML page, {reference directory: standalone SVG map})
    """
    from plasmid_surgery import __version__
    refs = run['references']
    contigs = reconstruction_contigs(run)
    colours = contig_colours(contigs)
    all_ids = [c for c, _ in contigs] + [r['contig_id'] for ref in refs for r in (ref['contig_sources'] or [])]
    labels = short_labels(all_ids)
    explained = reconstruction_explained(run, contigs)
    longest = max((r['length'] for r in refs), default=1)
    shared_scale = DRAW_WIDTH / longest
    sample = run['assessment'][0]['sample_id'] if run['assessment'] else ''
    title = title or f"plasmid-surgery: {sample}"

    panels, maps = [], {}
    for ref in refs:
        row = ref['row']
        fit_scale = DRAW_WIDTH / ref['length']
        maps[ref['dir']] = render_map(ref, colours, fit_scale, standalone=True, labels=labels,
                                      title=f"{sample} vs. {row['ref_id']} ({ref['length']:,} bp)")
        panels.append(
            f'<section class="panel" id="ref-{esc(row["ref_id"])}">'
            f'<h3>#{esc(row.get("rank", ""))} {esc(row["ref_id"])} <span class="na">({ref["length"]:,} bp)</span></h3>'
            f'<p class="stats">{esc(panel_stats(row))}</p>'
            f'<div class="scroll map-shared">{render_map(ref, colours, shared_scale, labels=labels)}</div>'
            f'<div class="scroll map-fit">{render_map(ref, colours, fit_scale, labels=labels)}</div>'
            '</section>')

    candidates = ''
    if run['candidates']:
        top = run['candidates'][:10]
        candidates = ('<h2>Candidate references</h2><p class="sub">From the plasmid database search '
                      f'(top {len(top)} of {len(run["candidates"])}; all in <code>candidates.tsv</code>).</p>'
                      '<div class="scroll"><table><thead><tr><th>Rank</th><th>Reference</th><th>Length</th><th>mash distance</th>'
                      '<th>Reconstruction aligned</th><th>Covered</th><th>Score</th><th>Assessed</th></tr></thead><tbody>'
                      + ''.join(f"<tr><td>{c['rank']}</td><td>{esc(c['ref_id'])}</td><td>{bp(c['ref_length'])}</td>"
                                f"<td>{float(c['mash_distance']):.4f}</td><td>{fraction(c['reconstruction_fraction_aligned'])}</td>"
                                f"<td>{fraction(c['assembly_ref_fraction_covered'] if c['assembly_ref_fraction_covered'] != 'NA' else c['ref_fraction_covered'])}</td>"
                                f"<td>{c['score']}</td><td>{'yes' if c['selected'] == 'True' else ''}</td></tr>" for c in top)
                      + '</tbody></table></div>')

    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)}</title><style>{CSS}</style></head>
<body class="scale-shared"><main>
<h1>{esc(title)}</h1>
<p class="sub">Reconstruction <code>{esc(sample)}</code> ({len(contigs)} contig{'s' if len(contigs) != 1 else ''},
{sum(l for _, l in contigs):,} bp) against {len(refs)} reference{'s' if len(refs) != 1 else ''}, ranked by score
(fraction of the reconstruction aligned × fraction of the reference covered). Hover over anything for details;
hovering over a contig highlights it everywhere.</p>
<h2>References</h2>
<div class="scroll">{summary_table(run)}</div>
{candidates}
<h2>Where each contig is</h2>
<p class="sub">Positions on each reference (kb; − for reverse strand). "alt." counts alternative placements
(e.g. repeat copies); a percentage is the fraction of the contig aligned, if not all of it.</p>
<div class="matrix-wrap">{contig_matrix(run, contigs, colours)}</div>
<h2>What the references explain</h2>
{reconstruction_section(run, contigs, colours, labels, explained)}
<h2>Maps</h2>
<div class="controls"><span>Scale:</span>
<label><input type="radio" name="scale" value="shared" checked> same for every reference</label>
<label><input type="radio" name="scale" value="fit"> fit each reference</label></div>
<div class="legend">{legend(run)}</div>
{''.join(panels)}
<footer>Generated by plasmid-surgery {esc(__version__)} on {date.today().isoformat()} from <code>{esc(run['outdir'])}</code>.
Each reference's map is also saved as <code>map.svg</code> in its directory.</footer>
</main><script>{JS}</script></body></html>
"""
    return page, maps


def write_report(outdir, output=None, title=None) -> Path:
    """Write <outdir>/[prefix]report.html and a map.svg per reference; returns the report's path."""
    run = load_run(outdir)
    page, maps = render_report(run, title)
    prefix = run['prefix']
    path = Path(output) if output else Path(outdir) / f"{prefix}report.html"
    path.write_text(page)
    for directory, svg in maps.items():
        (directory / f"{prefix}map.svg").write_text(svg)
    return path


def build_parser() -> argparse.ArgumentParser:
    from plasmid_surgery import __version__
    parser = argparse.ArgumentParser(
        prog='plasmid-surgery-plot',
        description="Draw a self-contained HTML report of a plasmid-surgery run (all references, on linear maps), "
                    "and an SVG map per reference.")
    parser.add_argument('outdir', help="A plasmid-surgery output directory")
    parser.add_argument('-o', '--output', help="Report path (default: <outdir>/[prefix_]report.html)")
    parser.add_argument('--title', help="Report title (default: 'plasmid-surgery: <sample_id>')")
    parser.add_argument('--log-level', default='info', help="Log level (default: 'info')")
    parser.add_argument('--version', action='version', version=f"%(prog)s {__version__}")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    from plasmid_surgery.__main__ import setup_logging
    setup_logging(args.log_level)
    try:
        path = write_report(args.outdir, args.output, args.title)
    except FileNotFoundError as e:
        logger.error(str(e))
        return 1
    logger.info(f"Wrote {path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
