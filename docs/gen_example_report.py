"""Write an example HTML report into the site (example-report/) at build time (mkdocs-gen-files).

Runs plasmid-surgery on the example in tests/scenarios/report_example.py (one reconstruction against
two references, with the assembly, bins and annotations) and copies the report and a map into the docs.
"""

import sys
import tempfile
from pathlib import Path

import mkdocs_gen_files

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT)]

from plasmid_surgery.report import write_report  # noqa: E402
from tests.scenarios.report_example import build_run  # noqa: E402

with tempfile.TemporaryDirectory(prefix="plasmid-surgery-example-") as tmp:
    out = build_run(Path(tmp))
    page = write_report(out, title="plasmid-surgery: example report").read_text()
    page = page.replace(str(out), "results/sample1")   # don't publish the temporary path
    with mkdocs_gen_files.open("example-report/report.html", "w") as fh:
        fh.write(page)
    with mkdocs_gen_files.open("example-report/map.svg", "w") as fh:
        fh.write((out / "pREF1" / "map.svg").read_text())

# An example cohort: five samples against one reference, each with a known difference
from plasmid_surgery.cohort import write_cohort  # noqa: E402
from tests.scenarios.cohort_example import build_runs  # noqa: E402

with tempfile.TemporaryDirectory(prefix="plasmid-surgery-cohort-") as tmp:
    root = Path(tmp)
    runs = build_runs(root)
    page = write_cohort(runs, root / "cohort.html", title="plasmid-surgery: example cohort").read_text()
    with mkdocs_gen_files.open("example-cohort/cohort.html", "w") as fh:
        fh.write(page)
    for run in runs:   # the per-sample reports the cohort links to
        sample_page = (run / "report.html").read_text().replace(str(run), f"results/{run.name}")
        with mkdocs_gen_files.open(f"example-cohort/results/{run.name}/report.html", "w") as fh:
            fh.write(sample_page)
