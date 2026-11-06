"""Figures for the Behaviour report, drawn from each scenario's outcome (see tests/scenarios).

Every figure is drawn from the files plasmid-surgery wrote in that scenario's run (summary.tsv,
gaps.bed, events.tsv, ...), so it shows what the code actually did.
"""

from __future__ import annotations

import io
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch, Polygon, Rectangle, Wedge  # noqa: E402

PALETTE = ["#2f6db3", "#d1663a", "#3a9a5b", "#8a5cb8", "#b8963a", "#5c8a8a", "#b84a6a", "#4a6ab8"]
INK, MUTED, FAINT = "#222222", "#666666", "#cccccc"
GAP_COLOURS = {"unique": "#c0392b", "repeat": "#e69f00"}
CLASS_COLOURS = {"in_reconstruction": "#2f6db3", "missed": "#d1663a", "shared": "#999999", "repeat": "#8a5cb8"}
STATUS_COLOURS = {"present": "#3a9a5b", "inferred": "#9ccc9c", "partial": "#e6b800", "absent": "#c0392b"}

plt.rcParams.update({
    "font.size": 9,
    "font.family": "sans-serif",
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "svg.fonttype": "none",
    "svg.hashsalt": "plasmid-surgery",  # stable ids, so regenerated figures only differ when they change
})


def svg(fig) -> str:
    buf = io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight", facecolor="white", metadata={"Date": None})
    plt.close(fig)
    return buf.getvalue()


def colours(names):
    return {n: PALETTE[i % len(PALETTE)] for i, n in enumerate(sorted(set(names)))}


def wrap(start, end, length):
    """Split an interval that may run past the origin into pieces within [0, length]."""
    if end <= length:
        return [(start, end)]
    return [(start, length), (0, end - length)]


def lanes(intervals, length):
    """Greedy lane assignment so intervals in the same lane don't overlap. Returns lane per interval."""
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


def arrow(ax, start, end, y, height, colour, *, reverse=False, length=1, hatch=None, dashed=False, fill=True):
    """A contig as an arrow from start to end at height y (pointing left if reverse)."""
    tip = min((end - start) * 0.35, length * 0.012)
    top, bottom = y + height / 2, y - height / 2
    if reverse:
        xy = [(start, y), (start + tip, top), (end, top), (end, bottom), (start + tip, bottom)]
    else:
        xy = [(start, bottom), (start, top), (end - tip, top), (end, y), (end - tip, bottom)]
    ax.add_patch(Polygon(xy, closed=True, fc=colour if fill else "white", alpha=0.85 if fill and not hatch else 1,
                         ec=colour, lw=1.2, ls="--" if dashed else "-", hatch=hatch))


def ref_track(ax, m, y):
    """The reference as a bar, with gaps coloured by type."""
    length = m["ref_length"]
    ax.add_patch(Rectangle((0, y - 0.15), length, 0.3, fc=FAINT, ec="none"))
    for g in m["gaps"]:
        for s, e in wrap(g["start"], g["end"], length):
            ax.add_patch(Rectangle((s, y - 0.15), e - s, 0.3, fc=GAP_COLOURS.get(g["type"], MUTED), ec="none",
                                   hatch="////" if g["type"] == "repeat" else None, alpha=0.9))
    for e in m["events"]:
        if e["type"] == "insertion":
            ax.plot(e["start"], y + 0.32, "v", color=INK, ms=7, zorder=4)
            ax.annotate(f"+{e['length'] / 1000:.1f} kb", (e["start"], y + 0.32), xytext=(0, 6),
                        textcoords="offset points", ha="center", fontsize=7.5, color=INK)
        elif e["type"] == "deletion":
            ax.annotate(f"−{e['length'] / 1000:.1f} kb", ((e["start"] + e["end"]) / 2, y - 0.15), xytext=(0, -11),
                        textcoords="offset points", ha="center", fontsize=7.5, color=GAP_COLOURS["unique"])
    ax.text(-0.01 * length, y, "reference", ha="right", va="center", fontsize=8, color=MUTED)


def contig_colours(m):
    """One colour per contig, shared by the placement and truth tracks."""
    return colours([p["contig"] for p in m["placements"]] + list(m["truth"]))


def label(ax, start, end, y, text, length, *, on_fill):
    """A contig's name: inside its arrow if it fits, otherwise just above it."""
    if (end - start) >= (len(text) * 0.0105 + 0.01) * length:
        ax.text((start + end) / 2, y, text, ha="center", va="center", fontsize=7.5,
                color="white" if on_fill else INK, zorder=5)
    else:
        ax.text((start + end) / 2, y + 0.2, text, ha="center", va="bottom", fontsize=7, color=INK, zorder=5)


def placements_track(ax, m, y0):
    """Contig placements below the reference, one lane per non-overlapping set. Returns the lowest y used."""
    length = m["ref_length"]
    order = sorted(m["placements"], key=lambda p: ({"placement": 0, "excluded": 1, "alternative": 2, "unused": 3}[p["kind"]], p["start"]))
    colour = contig_colours(m)
    lane_of = lanes([(p["start"], p["end"]) for p in order], length)
    for p, lane in zip(order, lane_of):
        y = y0 - 0.5 * lane
        c = colour[p["contig"]]
        for s, e in wrap(p["start"], p["end"], length):
            if p["kind"] == "placement":
                arrow(ax, s, e, y, 0.3, c, reverse=p["reverse"], length=length)
            elif p["kind"] == "alternative":
                arrow(ax, s, e, y, 0.3, c, reverse=p["reverse"], length=length, hatch="....", fill=False)
            else:
                arrow(ax, s, e, y, 0.3, MUTED, reverse=p["reverse"], length=length, dashed=True, fill=False)
        s, e = max(wrap(p["start"], p["end"], length), key=lambda piece: piece[1] - piece[0])
        label(ax, s, e, y, p["contig"], length, on_fill=p["kind"] == "placement")
    return y0 - 0.5 * (max(lane_of, default=0))


def truth_track(ax, m, y):
    length = m["ref_length"]
    colour = contig_colours(m)
    for name, (start, end) in m["truth"].items():
        for s, e in wrap(start, end, length):
            ax.add_patch(Rectangle((s, y - 0.05), e - s, 0.1, fc=colour[name], ec="none", alpha=0.6))
    ax.text(-0.01 * length, y, "truth", ha="right", va="center", fontsize=8, color=MUTED)


def features_track(ax, m, y):
    length = m["ref_length"]
    for f in m["features"]:
        c = STATUS_COLOURS.get(f["status"], MUTED)
        ax.add_patch(Rectangle((f["start"], y - 0.12), f["end"] - f["start"], 0.24, fc=c, ec=INK, lw=0.5,
                               hatch="////" if f["status"] == "inferred" else None))
        ax.annotate(f["name"], ((f["start"] + f["end"]) / 2, y + 0.12), xytext=(0, 3), textcoords="offset points",
                    ha="center", fontsize=7, color=INK, rotation=20)
    ax.text(-0.01 * length, y, "ref features", ha="right", va="center", fontsize=8, color=MUTED)


def finish(ax, m, top, bottom, title):
    length = m["ref_length"]
    ax.set_xlim(-0.005 * length, length * 1.005)
    ax.set_ylim(bottom - 0.4, top + 0.6)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ticks = range(0, length + 1, 5000 if length <= 40000 else 10000)
    ax.set_xticks(list(ticks))
    ax.set_xticklabels([f"{t / 1000:g} kb" for t in ticks])
    ax.set_title(title, loc="left", fontsize=10, color=INK)


def legend(fig, *, gaps=True, alternative=False, excluded=False, classes=None, statuses=None, ncol=None):
    handles = []
    if gaps:
        handles += [Patch(fc=GAP_COLOURS["unique"], label="unique gap"),
                    Patch(fc=GAP_COLOURS["repeat"], hatch="////", label="repeat gap")]
    if alternative:
        handles.append(Patch(fc="white", ec=MUTED, hatch="....", label="alternative placement"))
    if excluded:
        handles.append(Patch(fc="white", ec=MUTED, ls="--", label="not in scaffold"))
    for name in classes or []:
        handles.append(Patch(fc=CLASS_COLOURS[name], label=name))
    for name in statuses or []:
        handles.append(Patch(fc=STATUS_COLOURS[name], hatch="////" if name == "inferred" else None, label=name))
    fig.legend(handles=handles, loc="upper center", ncol=ncol or len(handles), frameon=False, fontsize=8,
               bbox_to_anchor=(0.5, 0.0))


def linear_map(m, title, *, truth=False, features=False, alternative=False, excluded=False, height=2.6):
    fig, ax = plt.subplots(figsize=(9, height))
    y = 0.0
    top = y + (0.5 if any(e["type"] == "insertion" for e in m["events"]) else 0)
    if features and m["features"]:
        features_track(ax, m, y + 0.9)
        top = y + 1.4
    ref_track(ax, m, y)
    bottom = placements_track(ax, m, y - 0.6)
    if truth and m["truth"]:
        bottom -= 0.5
        truth_track(ax, m, bottom)
    finish(ax, m, top, bottom, title)
    legend(fig, alternative=alternative, excluded=excluded, statuses=list(STATUS_COLOURS) if features else None)
    return fig


# --- one function per scenario: outcome -> [(caption, svg), ...] ---

def origin(outcome):
    linear, m = outcome.data["linear"], outcome.data["map"]
    length = m["ref_length"]
    colour = contig_colours(m)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.4))

    def ring(ax, title):
        ax.add_patch(Wedge((0, 0), 1.0, 0, 360, width=0.06, fc=FAINT, ec="none"))
        ax.plot([0, 0], [0.7, 1.08], color=INK, lw=1.2)
        ax.text(0, 0.64, "origin", ha="center", va="top", fontsize=8, color=INK)
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(title, fontsize=10, color=INK)

    def arc(ax, start, end, radius, c, label=None, hatch=None):
        theta1, theta2 = 90 - 360 * end / length, 90 - 360 * start / length
        ax.add_patch(Wedge((0, 0), radius, theta1, theta2, width=0.14, fc=c, ec="white", lw=1, alpha=0.9, hatch=hatch))
        if label:
            mid = math.radians((theta1 + theta2) / 2)
            ax.text(1.28 * math.cos(mid), 1.28 * math.sin(mid), label, ha="center", va="center", fontsize=8.5,
                    color=c, weight="bold")

    ring(axes[0], "c1, aligned to the reference as given")
    for k, a in enumerate(linear["arcs"], 1):
        arc(axes[0], a["start"], a["end"], 0.88, colour["c1"], f"c1 ({k})")
    ring(axes[1], "plasmid-surgery's placements")
    for p in m["placements"]:
        arc(axes[1], p["start"], p["end"], 0.88, colour[p["contig"]], p["contig"])
    fig.text(0.5, 0.02, f"{length / 1000:g} kb reference; positions run clockwise from the origin",
             ha="center", fontsize=8, color=MUTED)
    return [("Left: minimap2 alone splits c1 into two alignments at the origin. Right: plasmid-surgery places "
             "c1 as one alignment across the origin, alongside c2 and c3.", svg(fig))]


def repeats(outcome):
    fig = linear_map(outcome.data["map"], "Placements on the reference", truth=True, alternative=True, height=3.0)
    return [("The collapsed repeat contig's primary placement (solid) and its two alternative placements (dotted), "
             "at the three repeat copies. The copies not covered by a primary placement are repeat gaps. "
             "'truth' shows where each unique contig really belongs.", svg(fig))]


def overlaps(outcome):
    fig = linear_map(outcome.data["map"], "Placements on the reference", truth=True, excluded=True, height=2.8)
    return [("a, b and c overlap their neighbours slightly, and c wraps the origin. d lies entirely within b, so "
             "it's left out of the scaffold (dashed).", svg(fig))]


def indels(outcome):
    fig = linear_map(outcome.data["map"], "Placements on the reference", height=2.4)
    return [("x spans the 2 kb deletion, which is a unique gap; the insertion in y (▼) is sequence the reference "
             "doesn't have. The 4 bp deletion is too small to show, or to count as a gap.", svg(fig))]


def assembly(outcome):
    m = outcome.data["map"]
    length = m["ref_length"]
    fig, ax = plt.subplots(figsize=(9, 2.8))
    ref_track(ax, m, 0)
    rows = [r for r in outcome.data["sources"]]
    spans = []
    for r in rows:
        s, e = int(r["ref_start"]), int(r["ref_end"])
        spans.append((s, e + length if e < s else e))
    for r, (s, e), lane in zip(rows, spans, lanes(spans, length)):
        y = -0.6 - 0.5 * lane
        c = CLASS_COLOURS[r["classification"]]
        arrow(ax, s, e, y, 0.3, c, length=length)
        text = r["contig_id"] if r["classification"] == "in_reconstruction" else f"{r['contig_id']} ({r['bin']})"
        label(ax, s, e, y, text, length, on_fill=True)
    finish(ax, m, 0, -0.6 - 0.5 * max(lanes(spans, length)), "Assembly contigs on the reference")
    legend(fig, gaps=True, classes=list(CLASS_COLOURS)[:3])
    return [("The reconstruction (blue) covers A and C, leaving unique gaps at B and D. In the assembly, B is a "
             "missed contig and the contig with D is shared (only 75% on the reference); a chromosome contig "
             "shares a small part of A.", svg(fig))]


def features(outcome):
    fig = linear_map(outcome.data["map"], "Reference features and placements", features=True, alternative=True, height=3.2)
    return [("Reference features coloured by their status in the reconstruction. geneDeleted is in sequence the "
             "sample lacks; geneMisbinned is in the region left uncovered by the reconstruction, whose contig is in "
             "the chromosome bin.", svg(fig))]


def reference_choice(outcome):
    maps = outcome.data["maps"]
    fig, axes = plt.subplots(2, 1, figsize=(9, 3.6), gridspec_kw={"height_ratios": [1, 1]})
    for ax, (name, title) in zip(axes, (("small", "Against the small reference"), ("large", "Against the large reference"))):
        m = maps[name]
        ref_track(ax, m, 0)
        bottom = placements_track(ax, m, -0.6)
        finish(ax, m, 0, bottom, title)
        ax.set_xlim(-0.005 * maps["large"]["ref_length"], maps["large"]["ref_length"] * 1.005)
    fig.tight_layout()
    legend(fig, gaps=True)
    return [("The same two contigs against the two references, drawn to the same scale.", svg(fig))]


def find_refs(outcome):
    import numpy as np
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), sharey=True)
    panels = (("Reconstruction A (complete)", outcome.data["a"], None),
              ("Reconstruction B (a fragment)", outcome.data["b"], outcome.data["b_asm"]))
    x = np.linspace(0.05, 1, 200)
    for ax, (title, rows, with_assembly) in zip(axes, panels):
        for level in (0.25, 0.5, 0.75, 0.9):
            ax.plot(x, level / x, color=FAINT, lw=0.8, zorder=0)
            ax.text(level / 1.06, 1.06, f"{level:g}", fontsize=6.5, color=MUTED, ha="right", va="bottom")
        ax.text(0.02, 1.06, "score", fontsize=6.5, color=MUTED, ha="left", va="bottom")
        colour = colours(rows)
        # Points something else will be drawn just above get their label below
        arrived = [(a["aligned"], a["covered"]) for a in (with_assembly or {}).values()]
        for name, r in rows.items():
            crowded = any(abs(ax_ - r["aligned"]) < 0.05 and 0 < ay - r["covered"] < 0.12 for ax_, ay in arrived)
            ax.plot(r["aligned"], r["covered"], "o", color=colour[name], ms=8, zorder=3)
            ax.annotate(f"{name} (#{r['rank']})", (r["aligned"], r["covered"]), xytext=(-8, -14 if crowded else 6),
                        textcoords="offset points", ha="right", fontsize=7.5, color=INK)
            if with_assembly and name in with_assembly:
                after = with_assembly[name]
                if abs(after["covered"] - r["covered"]) > 0.01:
                    ax.annotate("", xy=(after["aligned"], after["covered"]), xytext=(r["aligned"], r["covered"]),
                                arrowprops=dict(arrowstyle="->", color=colour[name], lw=1.3, ls="--"))
                    ax.annotate(f"{name} with --assembly (#{after['rank']})", (after["aligned"], after["covered"]),
                                xytext=(-8, 4), textcoords="offset points", ha="right", fontsize=7.5,
                                color=colour[name])
        ax.set_xlim(0, 1.05)
        ax.set_ylim(0, 1.08)
        ax.set_xlabel("fraction of the reconstruction aligned")
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        for side in ("right", "top"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("fraction of the reference covered")
    fig.tight_layout()
    return [("Each candidate's two fractions, with curves of equal score (their product). For the fragment, "
             "--assembly replaces the reconstruction's coverage of each candidate with the whole assembly's, "
             "and the large plasmid moves to the top.", svg(fig))]


FIGURES = {
    "find-refs": find_refs,
    "origin": origin,
    "repeats": repeats,
    "overlaps": overlaps,
    "indels": indels,
    "assembly": assembly,
    "features": features,
    "reference-choice": reference_choice,
}
