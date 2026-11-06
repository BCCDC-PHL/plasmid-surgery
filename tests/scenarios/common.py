"""Shared helpers for behaviour scenarios: synthetic plasmids, FASTA files, and CLI runs.

Each scenario builds a small plasmid from random sequence (seeded, so it's the same every time), cuts
it into contigs the way a short-read assembly would, and runs the real ``plasmid-surgery`` command on
them. Because the scenario made the contigs, it knows exactly where each one belongs, so the report
can say what the right answer is and show that the tool finds it.
"""

from __future__ import annotations

import csv
import logging
import random
import shlex
import shutil
import subprocess
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from plasmid_surgery.__main__ import main as plasmid_surgery

COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def revcomp(seq: str) -> str:
    return seq[::-1].translate(COMPLEMENT)


def random_seq(rng: random.Random, n: int) -> str:
    return ''.join(rng.choice('ACGT') for _ in range(n))


def mutate(rng: random.Random, seq: str, n: int) -> str:
    """``seq`` with ``n`` random substitutions."""
    s = list(seq)
    for i in rng.sample(range(len(s)), n):
        s[i] = rng.choice('ACGT'.replace(s[i].upper(), ''))
    return ''.join(s)


def plural(n: int, word: str, words: str | None = None) -> str:
    return f"{n} {word if n == 1 else (words or word + 's')}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def read_fasta(path: Path) -> dict[str, str]:
    seqs, seq_id = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                seq_id = line[1:].split()[0]
                seqs[seq_id] = []
            elif seq_id:
                seqs[seq_id].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


def minimap2_available() -> bool:
    return shutil.which('minimap2') is not None


@dataclass
class Workspace:
    """A temporary directory where a scenario writes FASTA files and runs commands.

    Every command is recorded with paths relative to the workspace, as it would be typed there.
    """

    root: Path
    commands: list[str] = field(default_factory=list)

    def path(self, name: str) -> Path:
        return self.root / name

    def write_fasta(self, name: str, records: Iterable[tuple[str, str]]) -> Path:
        """``records`` are (header, sequence); the header may include a description after the ID."""
        path = self.path(name)
        with open(path, 'w') as fh:
            for header, seq in records:
                fh.write(f">{header}\n")
                for i in range(0, len(seq), 70):
                    fh.write(seq[i:i + 70] + '\n')
        return path

    def write_text(self, name: str, text: str) -> Path:
        path = self.path(name)
        path.write_text(text)
        return path

    def _shown(self, args: Sequence[str]) -> list[str]:
        return [str(Path(a).relative_to(self.root)) if a.startswith(str(self.root)) else a for a in args]

    def run(self, *args: str | Path) -> Path:
        """Run ``plasmid-surgery`` with ``args``; returns the output directory (``-o``)."""
        argv = [str(a) for a in args]
        self.commands.append(shlex.join(['plasmid-surgery', *self._shown(argv)]))
        logging.disable(logging.WARNING)   # expected warnings are noise here
        try:
            code = plasmid_surgery([*argv, '--log-level', 'error'])
        finally:
            logging.disable(logging.NOTSET)
        if code != 0:
            raise RuntimeError(f"`{self.commands[-1]}` failed")
        return Path(argv[argv.index('-o') + 1])

    def minimap2(self, *args: str | Path) -> Path:
        """Run minimap2 directly (for comparisons); ``args`` must include ``-o <out>``."""
        argv = [str(a) for a in args]
        self.commands.append(shlex.join(['minimap2', *self._shown(argv)]))
        subprocess.run(['minimap2', *argv], check=True, capture_output=True)
        return Path(argv[argv.index('-o') + 1])


def workspace() -> tuple[Workspace, tempfile.TemporaryDirectory]:
    tmp = tempfile.TemporaryDirectory(prefix='plasmid-surgery-scenario-')
    return Workspace(Path(tmp.name)), tmp


def outputs(outdir: Path, ref_id: str = None) -> dict:
    """Everything a run wrote for one reference (``<outdir>/<ref_id>/``; the only one if not given),
    parsed (missing optional outputs are absent)."""
    if ref_id is None:
        [outdir] = [d for d in Path(outdir).iterdir() if d.is_dir()]
    else:
        outdir = Path(outdir) / ref_id
    out = {
        'assessment': read_tsv(outdir / 'assessment.tsv')[0],
        'summary': read_tsv(outdir / 'summary.tsv'),
        'events': read_tsv(outdir / 'events.tsv'),
        'gaps': [dict(zip(('ref', 'start', 'end', 'name'), line.split('\t')))
                 for line in (outdir / 'gaps.bed').read_text().splitlines()],
        'scaffold': next(iter(read_fasta(outdir / 'scaffold.fa').values()), ''),
    }
    for name in ('contig_sources', 'ref_features', 'features'):
        if (outdir / f'{name}.tsv').exists():
            out[name] = read_tsv(outdir / f'{name}.tsv')
    return out


def plasmid_map(out: dict, ref_length: int, truth: Mapping[str, tuple[int, int]] = None) -> dict:
    """Figure data for a linear map of the ref: contig placements, gaps, events and features.

    :param truth: optional contig ID -> (start, end) on the ref, where each contig really belongs
    """
    placements = []
    for row in out['summary']:
        if row['is_mapped'] != 'True':
            continue
        start, end = int(row['alignment_start']), int(row['alignment_end'])
        if end < start:
            end += ref_length
        kind = ('alternative' if row['is_secondary'] == 'True' else 'placement')
        if row['included_in_output'] != 'True':
            kind = 'unused' if kind == 'alternative' else 'excluded'
        placements.append({'contig': row['contig_id'], 'start': start, 'end': end, 'kind': kind,
                           'reverse': row['reverse_complement'] == 'True'})
    gaps = []
    for g in out['gaps']:
        gaps.append({'start': int(g['start']), 'end': int(g['end']), 'type': g['name'].rsplit('_', 1)[-1]})
    events = [{'type': e['event_type'], 'start': int(e['ref_start']), 'end': int(e['ref_end']), 'length': int(e['length'])}
              for e in out['events']]
    features = [{'name': f['name'], 'start': int(f['ref_start']), 'end': int(f['ref_end']), 'status': f['reconstruction_status']}
                for f in out.get('ref_features', [])]
    return {'ref_length': ref_length, 'placements': placements, 'gaps': gaps, 'events': events,
            'features': features, 'truth': dict(truth or {})}


@dataclass
class Outcome:
    """What a scenario did and found.

    :ivar facts: named results, checked by the tests and quoted in the report's text.
    :ivar data: anything else the report's figures and tables need.
    """

    commands: list[str]
    facts: dict
    data: dict = field(default_factory=dict)


@dataclass(frozen=True)
class Scenario:
    """One behaviour to demonstrate: its title, explanation and how to run it.

    ``text`` is Markdown with ``{fact}`` placeholders filled from the outcome's facts, and
    ``{table_name}`` placeholders filled from ``outcome.data["tables"]``.
    """

    slug: str
    title: str
    summary: str
    text: str
    run: object  # () -> Outcome
    test: str  # the test that checks it, for the report


def markdown_table(columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> str:
    lines = ['| ' + ' | '.join(f'`{c}`' for c in columns) + ' |', '|' + '---|' * len(columns)]
    lines += ['| ' + ' | '.join(str(r.get(c, '')) for c in columns) + ' |' for r in rows]
    return '\n'.join(lines)


def render(scenario: Scenario, outcome: Outcome) -> str:
    """The scenario's text, with facts and tables filled in."""
    tables = {name: markdown_table(cols, rows) for name, (cols, rows) in outcome.data.get('tables', {}).items()}
    return scenario.text.format(**outcome.facts, **tables)


def run_refs(ws: Workspace, *args: str | Path) -> Path:
    """Run ``plasmid-surgery-refs`` in a workspace; returns the output directory (``-o``)."""
    from plasmid_surgery.refs import main as plasmid_surgery_refs
    argv = [str(a) for a in args]
    ws.commands.append(shlex.join(['plasmid-surgery-refs', *ws._shown(argv)]))
    logging.disable(logging.WARNING)
    try:
        code = plasmid_surgery_refs([*argv, '--log-level', 'error'])
    finally:
        logging.disable(logging.NOTSET)
    if code != 0:
        raise RuntimeError(f"`{ws.commands[-1]}` failed")
    return Path(argv[argv.index('-o') + 1])


def mash_sketch(ws: Workspace, fasta: Path, out_prefix: str) -> Path:
    """``mash sketch -i`` (one sketch per sequence) of a FASTA, recorded like the other commands."""
    argv = ['sketch', '-i', '-k', '21', '-s', '1000', '-o', str(ws.path(out_prefix)), str(fasta)]
    ws.commands.append(shlex.join(['mash', *ws._shown(argv)]))
    subprocess.run(['mash', *argv], check=True, capture_output=True)
    return ws.path(f"{out_prefix}.msh")


def mash_available() -> bool:
    return shutil.which('mash') is not None
