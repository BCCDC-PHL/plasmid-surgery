"""Write docs/reference/cli.md from the argparse parser, at build time (mkdocs-gen-files).

The page is generated, so it always matches the installed code.
"""

import argparse

import mkdocs_gen_files

from plasmid_surgery.__main__ import build_parser

parser = build_parser()
parser.formatter_class = lambda prog: argparse.HelpFormatter(prog, width=100)
lines = [
    "# Command line",
    "",
    "Generated from `plasmid-surgery --help`.",
    "",
    "```text",
    parser.format_help().rstrip(),
    "```",
    "",
]
with mkdocs_gen_files.open("reference/cli.md", "w") as fh:
    fh.write("\n".join(lines))
