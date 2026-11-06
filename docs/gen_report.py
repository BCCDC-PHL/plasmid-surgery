"""Write the Behaviour report (docs/behaviour/) at build time (mkdocs-gen-files).

Runs every scenario in tests/scenarios with the installed code, and writes one page per scenario,
with its figures, the exact commands, and the test that checks it.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

import mkdocs_gen_files

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT), str(ROOT / "docs")]

from plasmid_surgery import __version__  # noqa: E402
from report_figures import FIGURES  # noqa: E402
from tests.scenarios import SCENARIOS  # noqa: E402
from tests.scenarios.common import render  # noqa: E402

OUT = "behaviour"
FIGURE_MARKER = "[[figure]]"


def page(scenario, outcome) -> str:
    lines = [f"# {scenario.title}", "", '!!! abstract "In short"', f"    {scenario.summary}", ""]
    figures = []
    for k, (caption, svg) in enumerate(FIGURES.get(scenario.slug, lambda o: [])(outcome), start=1):
        name = f"{scenario.slug}-{k}.svg"
        with mkdocs_gen_files.open(f"{OUT}/{name}", "w") as fh:
            fh.write(svg)
        figures.append("\n".join(['<figure markdown="span">', f"  ![{caption}]({name})",
                                   f"  <figcaption>{caption}</figcaption>", "</figure>"]))
    # Figures go where the text has a [[figure]] marker, or else at the end
    body = render(scenario, outcome).rstrip()
    if FIGURE_MARKER in body:
        body = body.replace(FIGURE_MARKER, "\n\n".join(figures), 1)
    else:
        body = "\n\n".join([body, *figures])
    lines += [body, ""]
    lines += ['??? note "Commands"', "    Run in a scratch directory, on the scenario's simulated FASTA files:", "", "    ```bash"]
    lines += [f"    {c}" for c in outcome.commands]
    lines += ["    ```", "", f"Checked by `{scenario.test}`.", ""]
    return "\n".join(lines)


index = [
    "# Behaviour report",
    "",
    "Each page here demonstrates one of plasmid-surgery's design decisions, and its trade-offs, on a",
    "small simulated plasmid. Each scenario builds a reference and a reconstruction from random",
    "sequence, so where every contig really belongs is known, then runs the real `plasmid-surgery`",
    "command (shown on each page) and reports what it found. The figures are drawn from the files that",
    "run wrote.",
    "",
    "The report is generated when the documentation is built, and every claim in it is also a test",
    "(`tests/test_scenarios.py`), so it always describes what the code actually does.",
    "",
    "| Scenario | In short |",
    "|---|---|",
]
for scenario in SCENARIOS:
    outcome = scenario.run()
    with mkdocs_gen_files.open(f"{OUT}/{scenario.slug}.md", "w") as fh:
        fh.write(page(scenario, outcome))
    index.append(f"| [{scenario.title}]({scenario.slug}.md) | {scenario.summary} |")
index += ["", f"*Generated with plasmid-surgery {__version__} on {date.today().isoformat()}.*", ""]
with mkdocs_gen_files.open(f"{OUT}/index.md", "w") as fh:
    fh.write("\n".join(index))
