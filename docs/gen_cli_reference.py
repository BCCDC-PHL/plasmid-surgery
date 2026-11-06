"""Write docs/reference/cli.md from the argparse parser, at build time (mkdocs-gen-files).

The page is generated, so it always matches the installed code.
"""

import argparse

import mkdocs_gen_files

from plasmid_surgery.__main__ import build_parser
from plasmid_surgery.refs import build_parser as build_refs_parser

lines = ["# Command line", "", "Generated from each tool's `--help`.", ""]
for tool, make_parser in (("plasmid-surgery", build_parser), ("plasmid-surgery-refs", build_refs_parser)):
    parser = make_parser()
    parser.formatter_class = lambda prog: argparse.HelpFormatter(prog, width=100)
    lines += [f"## `{tool}`", "", "```text", parser.format_help().rstrip(), "```", ""]
with mkdocs_gen_files.open("reference/cli.md", "w") as fh:
    fh.write("\n".join(lines))
