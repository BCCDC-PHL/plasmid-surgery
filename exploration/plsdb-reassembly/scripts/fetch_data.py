#!/usr/bin/env python3
"""Gather the public data for the validation reports: PLSDB plasmids, their read runs, and the reads.

Standard library only, so it runs anywhere with Python >= 3.10. Steps, each a subcommand:

  select    PLSDB plasmids to validate with (carbapenemase-carrying, plus optional random others),
            complete and circular, from PLSDB's metadata files.
  runs      Every read run for each plasmid's BioSample (ENA portal API), and for nanopore runs the
            flow cell chemistry (R9.4.1, R10.4.1, ...) from NCBI SRA's experiment metadata.
  subset    Choose BioSamples to download: those with both Illumina and nanopore runs, optionally
            only R10.4.1, within a size budget.
  download  The reads for a subset (from ENA, md5-checked, resumable), and the published sequences
            of every PLSDB plasmid from those BioSamples (the truth set).

API responses are cached under --cache (default .cache/), so reruns are fast and repeatable. Set
NCBI_API_KEY to raise NCBI's rate limit from 3 to 10 requests a second.

Example (small run):

  python fetch_data.py select --plsdb ../../test_input/plsdb -o plasmids.tsv
  python fetch_data.py runs -i plasmids.tsv -o runs.tsv --sample 200 --seed 1
  python fetch_data.py subset -i runs.tsv -o subset.tsv --chemistry R10.4.1 --max-biosamples 3 --max-gb 2
  python fetch_data.py download -i subset.tsv --plasmids plasmids.tsv --plsdb ../../test_input/plsdb -o data
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import logging
import os
import math
import random
import statistics
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

log = logging.getLogger("fetch_data")
csv.field_size_limit(sys.maxsize)

ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "cluster-muster-validation/0.1 (python urllib)"

ENA_FIELDS = (
    "run_accession", "sample_accession", "experiment_accession", "study_accession",
    "instrument_platform", "instrument_model", "library_layout", "library_strategy", "library_source",
    "read_count", "base_count", "fastq_ftp", "fastq_bytes", "fastq_md5", "submitted_ftp",
    "first_public", "library_construction_protocol", "experiment_title", "run_alias", "library_name",
)

RUN_FIELDS = (
    "biosample", "plasmids", "run_accession", "platform", "instrument_model", "library_layout",
    "library_strategy", "read_count", "base_count", "fastq_bytes", "first_public",
    "chemistry", "chemistry_evidence", "median_read_q", "fastq_ftp", "fastq_md5",
)

CARBAPENEM = "CARBAPENEM"

# Nanopore chemistry from free text in run metadata. Order matters: the first rule that matches wins,
# and the most specific come first. Kit numbers: kit 14 (SQK-xxx114) runs on R10.4.1; kit 12
# (SQK-xxx112) on R10.4; kits 9 and 10 (SQK-xxx109/110, SQK-RBK004, ...) on R9.4.1, so they're only
# used when no flow cell is named. Barcode expansions (EXP-NBD114 etc.) belong to kit 9 and say
# nothing about the flow cell, so they aren't used.
CHEMISTRY_RULES: tuple[tuple[str, str], ...] = (
    ("R10.4.1", r"r10\.4\.1|flo-(min|pro|flg)114|e8\.2|dna_r10\.4\.1"),
    ("R10.4", r"r10\.4(?!\.)|flo-(min|pro|flg)112"),
    ("R10.3", r"r10\.3|flo-(min|pro)111"),
    ("R9.4.1", r"r9\.4\.1|r9\.4(?!\.)|flo-(min|pro)10[67]|flo-pro002|flo-flg001|dna_r9\.4"),
    ("R9.5", r"r9\.5|flo-min107"),
    ("R10.4.1", r"sqk-[a-z]{3}114|\bkit ?14\b|\bv14\b"),
    ("R10.4", r"sqk-[a-z]{3}112|\bkit ?12\b"),
    ("R9.4.1", r"sqk-(lsk|rbk|nbd|rad|rpb|pcb|lsk)1(08|09|10)|sqk-(rbk|rad|rpb|rlk)00[0-9]"),
)


# ---------------------------------------------------------------- HTTP, with a cache and rate limits


class Client:
    """GET with an on-disk cache, retries, and a per-host minimum interval between requests."""

    def __init__(self, cache: Path, ncbi_key: str | None = None):
        self.cache = cache
        self.cache.mkdir(parents=True, exist_ok=True)
        self.ncbi_key = ncbi_key
        self._lock = threading.Lock()
        self._next = defaultdict(float)
        self._interval = {"eutils.ncbi.nlm.nih.gov": 0.11 if ncbi_key else 0.35, "www.ebi.ac.uk": 0.05}

    def get(self, url: str, params: dict, *, cache: bool = True) -> str:
        host = urllib.parse.urlsplit(url).hostname or ""
        key = hashlib.sha1((url + "?" + urllib.parse.urlencode(sorted(params.items()))).encode()).hexdigest()
        path = self.cache / host / key[:2] / key
        if cache and path.exists():
            return path.read_text()
        if host.endswith("ncbi.nlm.nih.gov") and self.ncbi_key:
            params = {**params, "api_key": self.ncbi_key}
        full = url + "?" + urllib.parse.urlencode(params)
        for attempt in range(5):
            self._wait(host)
            try:
                req = urllib.request.Request(full, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(req, timeout=60) as r:
                    text = r.read().decode()
                break
            except urllib.error.HTTPError as e:
                if e.code in (429, 500, 502, 503, 504) and attempt < 4:
                    time.sleep(2 ** attempt)
                    continue
                raise
            except (urllib.error.URLError, TimeoutError):
                if attempt == 4:
                    raise
                time.sleep(2 ** attempt)
        if cache:
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(".tmp")
            tmp.write_text(text)
            tmp.replace(path)
        return text

    def _wait(self, host: str) -> None:
        with self._lock:
            now = time.monotonic()
            start = max(now, self._next[host])
            self._next[host] = start + self._interval.get(host, 0.1)
        time.sleep(max(0.0, start - now))


# ---------------------------------------------------------------- select


def read_csv(path: Path, delimiter=","):
    with open(path, newline="") as fh:
        yield from csv.DictReader(fh, delimiter=delimiter)


def cmd_select(args) -> None:
    plsdb = Path(args.plsdb)
    genes: dict[str, set[str]] = defaultdict(set)
    for row in read_csv(plsdb / "amr.tsv", "\t"):
        if row["gene_symbol"].startswith("bla") and CARBAPENEM in row["antimicrobial_agent"].upper():
            genes[row["NUCCORE_ACC"]].add(row["gene_symbol"])
    biosample_acc = {r["BIOSAMPLE_UID"]: r["BIOSAMPLE_ACC"] for r in read_csv(plsdb / "biosample.csv")}

    carrying, others = [], []
    for row in read_csv(plsdb / "nuccore.csv"):
        if args.complete_circular and (row["NUCCORE_Topology"] != "circular"
                                       or row["NUCCORE_Completeness"] != "complete"):
            continue
        rec = {
            "nuccore_acc": row["NUCCORE_ACC"],
            "biosample": biosample_acc.get(row["BIOSAMPLE_UID"], ""),
            "length": row["NUCCORE_Length"],
            "topology": row["NUCCORE_Topology"],
            "completeness": row["NUCCORE_Completeness"],
            "create_date": row["NUCCORE_CreateDate"][:10],
            "carbapenemases": ",".join(sorted(genes.get(row["NUCCORE_ACC"], ()))),
            "set": "carbapenemase" if row["NUCCORE_ACC"] in genes else "other",
        }
        (carrying if rec["set"] == "carbapenemase" else others).append(rec)
    rng = random.Random(args.seed)
    chosen = carrying + (rng.sample(others, min(args.others, len(others))) if args.others else [])
    write_tsv(args.output, chosen, list(chosen[0]) if chosen else ["nuccore_acc"])
    n_bs = len({r["biosample"] for r in chosen if r["biosample"]})
    log.info("%d carbapenemase plasmids and %d others, from %d BioSamples, written to %s",
             len(carrying), len(chosen) - len(carrying), n_bs, args.output)


# ---------------------------------------------------------------- runs


def ena_runs(client: Client, biosample: str) -> list[dict]:
    text = client.get(ENA_FILEREPORT, {"accession": biosample, "result": "read_run",
                                       "fields": ",".join(ENA_FIELDS), "format": "tsv"})
    lines = [l for l in text.splitlines() if l.strip()]
    if len(lines) < 2:
        return []
    return list(csv.DictReader(lines, delimiter="\t"))


def ncbi_experiment_text(client: Client, runs: list[str]) -> dict[str, str]:
    """Free text describing each run's library and experiment, from NCBI SRA (in batches)."""
    out: dict[str, str] = {}
    for i in range(0, len(runs), 100):
        batch = runs[i:i + 100]
        xml = client.get(f"{NCBI_EUTILS}/efetch.fcgi", {"db": "sra", "id": ",".join(batch)})
        try:
            root = ET.fromstring(xml)
        except ET.ParseError:
            log.warning("unparseable SRA XML for %s", ",".join(batch))
            continue
        for pkg in root.iter("EXPERIMENT_PACKAGE"):
            texts = []
            for tag in ("TITLE", "DESIGN_DESCRIPTION", "LIBRARY_NAME", "LIBRARY_CONSTRUCTION_PROTOCOL"):
                texts += [e.text for e in pkg.iter(tag) if e.text]
            for attr in pkg.iter("EXPERIMENT_ATTRIBUTE"):
                texts.append(" ".join(e.text or "" for e in attr))
            for run in pkg.iter("RUN"):
                run_texts = list(texts)
                run_texts += [f.get("filename", "") for f in run.iter("SRAFile")]
                run_texts += [" ".join(e.text or "" for e in a) for a in run.iter("RUN_ATTRIBUTE")]
                out[run.get("accession")] = " | ".join(t.strip() for t in run_texts if t and t.strip())
    return out


def classify_chemistry(text: str) -> tuple[str, str]:
    """(chemistry, the matching text), or ("unknown", "") when nothing in ``text`` says."""
    low = text.lower()
    for chemistry, pattern in CHEMISTRY_RULES:
        m = re.search(pattern, low)
        if m:
            a, b = max(0, m.start() - 30), min(len(text), m.end() + 30)
            return chemistry, text[a:b].replace("\t", " ").replace("\n", " ")
    return "unknown", ""


def median_read_q(url: str, reads: int, max_bytes: int = 50_000_000) -> str:
    """Median per-read mean quality of the first ``reads`` reads, streamed (nothing is saved).

    A weak hint at chemistry when the metadata says nothing: R10.4.1 reads basecalled with recent
    hac/sup models are typically Q16 and up, R9.4.1 reads Q9 to Q14. Returns "constant" when every
    base has the same score (SRA "lite" files replace qualities with a placeholder), and "" on failure.
    """
    qs, seen = [], set()
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            limited = _Limited(r, max_bytes)
            with gzip.open(limited, "rt") as fh:
                for i, line in enumerate(fh):
                    if i % 4 != 3:
                        continue
                    line = line.strip()
                    if not line:
                        continue
                    seen.update(set(line[:200]))
                    errors = sum(10 ** (-(ord(c) - 33) / 10) for c in line) / len(line)
                    qs.append(-10 * math.log10(errors))
                    if len(qs) >= reads:
                        break
    except (OSError, EOFError, urllib.error.URLError, TimeoutError) as e:
        if not qs:
            log.warning("  couldn't read %s: %s", url, e)
            return ""
    if not qs:
        return ""
    return "constant" if len(seen) == 1 else f"{statistics.median(qs):.1f}"


class _Limited:
    """A read-only stream that stops after ``limit`` bytes."""

    def __init__(self, raw, limit):
        self.raw, self.left = raw, limit

    def read(self, n=-1):
        if self.left <= 0:
            return b""
        n = self.left if n is None or n < 0 else min(n, self.left)
        data = self.raw.read(n)
        self.left -= len(data)
        return data


def cmd_runs(args) -> None:
    client = Client(Path(args.cache), os.environ.get("NCBI_API_KEY"))
    plasmids: dict[str, list[str]] = defaultdict(list)
    for row in read_csv(Path(args.input), "\t"):
        if row.get("biosample"):
            plasmids[row["biosample"]].append(row["nuccore_acc"])
    biosamples = sorted(plasmids)
    if args.sample and args.sample < len(biosamples):
        biosamples = sorted(random.Random(args.seed).sample(biosamples, args.sample))
    log.info("looking up runs for %d BioSamples at ENA", len(biosamples))

    runs: dict[str, list[dict]] = {}
    with concurrent.futures.ThreadPoolExecutor(args.threads) as pool:
        for k, (bs, found) in enumerate(zip(biosamples, pool.map(lambda b: ena_runs(client, b), biosamples)), 1):
            runs[bs] = found
            if k % 250 == 0:
                log.info("  %d / %d", k, len(biosamples))

    ont = [r["run_accession"] for found in runs.values() for r in found if r["instrument_platform"] == "OXFORD_NANOPORE"]
    log.info("looking up chemistry for %d nanopore runs at NCBI", len(ont))
    ncbi_text = ncbi_experiment_text(client, ont)

    rows = []
    for bs in biosamples:
        if not runs[bs]:
            rows.append({"biosample": bs, "plasmids": ",".join(plasmids[bs]), "platform": "NONE"})
            continue
        for r in runs[bs]:
            chemistry, evidence = "", ""
            if r["instrument_platform"] == "OXFORD_NANOPORE":
                text = " | ".join(t for t in (ncbi_text.get(r["run_accession"], ""), r["library_construction_protocol"],
                                              r["experiment_title"], r["run_alias"], r["library_name"],
                                              r["submitted_ftp"]) if t)
                chemistry, evidence = classify_chemistry(text)
            rows.append({
                "biosample": bs, "plasmids": ",".join(plasmids[bs]), "run_accession": r["run_accession"],
                "platform": r["instrument_platform"], "instrument_model": r["instrument_model"],
                "library_layout": r["library_layout"], "library_strategy": r["library_strategy"],
                "read_count": r["read_count"], "base_count": r["base_count"],
                "fastq_bytes": sum(int(b) for b in r["fastq_bytes"].split(";") if b),
                "first_public": r["first_public"], "chemistry": chemistry, "chemistry_evidence": evidence,
                "fastq_ftp": r["fastq_ftp"], "fastq_md5": r["fastq_md5"],
            })
    if args.peek_quality:
        todo = [r for r in rows if r["platform"] == "OXFORD_NANOPORE" and r["fastq_ftp"]
                and (args.peek_all or r["chemistry"] == "unknown")]
        log.info("estimating read quality for %d nanopore runs (first %d reads each)", len(todo), args.peek_quality)
        with concurrent.futures.ThreadPoolExecutor(args.threads) as pool:
            for r, q in zip(todo, pool.map(lambda r: median_read_q("https://" + r["fastq_ftp"].split(";")[0],
                                                                   args.peek_quality), todo)):
                r["median_read_q"] = q
    write_tsv(args.output, rows, RUN_FIELDS)
    summarise_runs(rows)


def biosample_platforms(rows) -> dict[str, set[str]]:
    by_bs: dict[str, set[str]] = defaultdict(set)
    for r in rows:
        by_bs[r["biosample"]].add(r["platform"])
    return by_bs


def summarise_runs(rows) -> None:
    by_bs = biosample_platforms(rows)
    combos = defaultdict(int)
    for platforms in by_bs.values():
        combos[" + ".join(sorted(p for p in platforms if p != "NONE")) or "no runs"] += 1
    log.info("BioSamples by platforms:")
    for combo, n in sorted(combos.items(), key=lambda kv: -kv[1]):
        log.info("  %5d  %s", n, combo)
    chem = defaultdict(int)
    for r in rows:
        if r["platform"] == "OXFORD_NANOPORE":
            chem[r["chemistry"]] += 1
    log.info("nanopore runs by chemistry: %s", ", ".join(f"{k} {v}" for k, v in sorted(chem.items(), key=lambda kv: -kv[1])))


# ---------------------------------------------------------------- subset


def cmd_subset(args) -> None:
    rows = list(read_csv(Path(args.input), "\t"))
    by_bs: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        if r["platform"] != "NONE":
            by_bs[r["biosample"]].append(r)

    def wanted_ont(r):
        if r["platform"] != "OXFORD_NANOPORE" or int(r["fastq_bytes"] or 0) == 0:
            return False
        if args.chemistry == "any" or r["chemistry"] == args.chemistry:
            return True
        q = r.get("median_read_q", "")
        return bool(args.unknown_min_q and r["chemistry"] == "unknown" and re.fullmatch(r"[0-9.]+", q)
                    and float(q) >= args.unknown_min_q)

    def wanted_illumina(r):
        return r["platform"] == "ILLUMINA" and r["library_layout"] == "PAIRED" and int(r["fastq_bytes"] or 0) > 0

    candidates = []
    for bs, found in by_bs.items():
        ont = [r for r in found if wanted_ont(r)]
        ill = [r for r in found if wanted_illumina(r)]
        if not ont or (args.require_illumina and not ill):
            continue
        # One run of each kind per BioSample: the largest, which is usually the main sequencing run.
        pick = [max(ont, key=lambda r: int(r["fastq_bytes"]))]
        if ill:
            pick.append(max(ill, key=lambda r: int(r["fastq_bytes"])))
        size = sum(int(r["fastq_bytes"]) for r in pick)
        if args.max_run_gb and size > args.max_run_gb * 1e9:
            continue
        candidates.append((bs, pick, size))
    rng = random.Random(args.seed)
    rng.shuffle(candidates)

    chosen, total = [], 0
    for bs, pick, size in candidates:
        if args.max_biosamples and len(chosen) >= args.max_biosamples:
            break
        if args.max_gb and total + size > args.max_gb * 1e9:
            continue
        chosen.append((bs, pick, size))
        total += size
    out = [r for _, pick, _ in chosen for r in pick]
    write_tsv(args.output, out, RUN_FIELDS)
    log.info("%d candidate BioSamples; chose %d, %.2f GB of reads, written to %s",
             len(candidates), len(chosen), total / 1e9, args.output)


# ---------------------------------------------------------------- download


def download(url: str, dest: Path, md5: str | None, *, retries: int = 3) -> None:
    if dest.exists() and (not md5 or file_md5(dest) == md5):
        log.info("  have %s", dest.name)
        return
    part = dest.with_name(dest.name + ".part")
    for attempt in range(retries):
        start = part.stat().st_size if part.exists() else 0
        headers = {"User-Agent": USER_AGENT, **({"Range": f"bytes={start}-"} if start else {})}
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=120) as r, \
                    open(part, "ab" if start and r.status == 206 else "wb") as fh:
                while chunk := r.read(1 << 20):
                    fh.write(chunk)
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            log.warning("  %s: %s (attempt %d)", dest.name, e, attempt + 1)
            time.sleep(2 ** attempt)
            continue
        if md5 and file_md5(part) != md5:
            log.warning("  %s: md5 mismatch, starting again", dest.name)
            part.unlink()
            continue
        part.replace(dest)
        log.info("  got %s", dest.name)
        return
    raise RuntimeError(f"couldn't download {url}")


def file_md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def cmd_download(args) -> None:
    out = Path(args.outdir)
    rows = list(read_csv(Path(args.input), "\t"))
    biosamples = sorted({r["biosample"] for r in rows})

    if not args.skip_reads:
        client = Client(Path(args.cache), os.environ.get("NCBI_API_KEY"))
        checks = []
        for r in rows:
            dest_dir = out / "reads" / r["biosample"]
            dest_dir.mkdir(parents=True, exist_ok=True)
            log.info("%s %s (%s, %.0f MB)", r["biosample"], r["run_accession"], r["platform"], int(r["fastq_bytes"]) / 1e6)
            names = [url.rsplit("/", 1)[1] for url in r["fastq_ftp"].split(";") if url]
            source = "ena"
            if r["platform"] == "ILLUMINA" and args.illumina_source == "ncbi":
                if all((dest_dir / n).exists() and quality_values(dest_dir / n) > 1 for n in names):
                    log.info("  have %s", ", ".join(names))
                    source = "ncbi"
                elif fetch_normalized(client, args, r["run_accession"], dest_dir, names):
                    source = "ncbi"
                else:
                    log.warning("  no full-quality copy at NCBI; using ENA's")
            if source == "ena":
                for url, md5 in zip(r["fastq_ftp"].split(";"), r["fastq_md5"].split(";")):
                    download("https://" + url, dest_dir / url.rsplit("/", 1)[1], md5 or None)
            for n in names:
                q = quality_values(dest_dir / n)
                if q <= 1:
                    log.warning("  %s: every base has the same quality score (an SRA lite copy); "
                                "quality trimming and read correction won't work on it", n)
                checks.append({"biosample": r["biosample"], "run_accession": r["run_accession"],
                               "platform": r["platform"], "file": n, "source": source,
                               "quality_values": q, "real_qualities": "yes" if q > 1 else "no"})
        write_tsv(out / "reads" / "reads.tsv", checks, ["biosample", "run_accession", "platform", "file", "source",
                                                       "quality_values", "real_qualities"])

    # The truth set: every PLSDB plasmid from these BioSamples, not just the carbapenemase ones.
    accessions = truth_accessions(args, set(biosamples))
    if accessions:
        write_truth(args, out / "truth", accessions)


def quality_values(path: Path, reads: int = 2000) -> int:
    """How many distinct quality characters the first ``reads`` reads use (1 means placeholders)."""
    seen: set[str] = set()
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 3:
                seen.update(line.strip())
            if i >= reads * 4:
                break
    return len(seen)


def normalized_sra_url(client: Client, run: str) -> tuple[str, str] | None:
    """(url, md5) of NCBI's "SRA Normalized" copy of ``run`` (original qualities), free to download."""
    xml = client.get(f"{NCBI_EUTILS}/efetch.fcgi", {"db": "sra", "id": run})
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return None
    for f in root.iter("SRAFile"):
        if f.get("semantic_name") != "SRA Normalized":
            continue
        for alt in f.iter("Alternatives"):
            if alt.get("free_egress") == "worldwide" and alt.get("access_type") == "anonymous" \
                    and alt.get("url", "").startswith("https://"):
                return alt.get("url"), f.get("md5", "")
        if f.get("url", "").startswith("https://"):
            return f.get("url"), f.get("md5", "")
    return None


def fetch_normalized(client: Client, args, run: str, dest_dir: Path, names: list[str]) -> bool:
    """Download NCBI's full-quality copy of ``run`` and convert it to gzipped FASTQ as ``names``.

    ENA's FASTQ for many SRA runs is made from SRA "lite" files, whose qualities are all Q30
    placeholders. NCBI keeps the original ("normalized") copy, which fasterq-dump converts.
    """
    found = normalized_sra_url(client, run)
    if not found:
        return False
    fasterq = shutil.which(args.fasterq_dump)
    if not fasterq:
        sys.exit(f"{args.fasterq_dump} isn't on the PATH: install sra-tools, use --fasterq-dump, "
                 "or --illumina-source ena")
    url, md5 = found
    with tempfile.TemporaryDirectory(prefix=f".{run}.", dir=dest_dir) as tmp:
        sra = Path(tmp) / run
        download(url, sra, md5 or None)
        subprocess.run([fasterq, "--split-files", "--threads", str(args.threads), "--temp", tmp,
                        "--outdir", tmp, str(sra)], check=True, stdout=subprocess.DEVNULL)
        produced = sorted(Path(tmp).glob(f"{run}_*.fastq"))
        expected = sorted(n.removesuffix(".gz") for n in names)
        if [p.name for p in produced] != expected:
            log.warning("  fasterq-dump made %s, expected %s", [p.name for p in produced], expected)
            return False
        zipper = [shutil.which("pigz"), "-p", str(args.threads)] if shutil.which("pigz") else ["gzip"]
        for fq in produced:
            subprocess.run([*zipper, str(fq)], check=True)
            Path(str(fq) + ".gz").replace(dest_dir / (fq.name + ".gz"))
    log.info("  got %s (full quality, from NCBI)", ", ".join(names))
    return True


def truth_accessions(args, biosamples: set[str]) -> dict[str, list[dict]]:
    """Every PLSDB record from ``biosamples``, with its length, topology and completeness."""
    plsdb = Path(args.plsdb)
    uid_to_acc = {r["BIOSAMPLE_UID"]: r["BIOSAMPLE_ACC"] for r in read_csv(plsdb / "biosample.csv")}
    found: dict[str, list[dict]] = defaultdict(list)
    for row in read_csv(plsdb / "nuccore.csv"):
        bs = uid_to_acc.get(row["BIOSAMPLE_UID"])
        if bs in biosamples:
            found[bs].append({"nuccore_acc": row["NUCCORE_ACC"], "length": row["NUCCORE_Length"],
                              "topology": row["NUCCORE_Topology"], "completeness": row["NUCCORE_Completeness"]})
    return found


def write_truth(args, outdir: Path, accessions: dict[str, list[dict]]) -> None:
    client = Client(Path(args.cache), os.environ.get("NCBI_API_KEY"))
    outdir.mkdir(parents=True, exist_ok=True)
    meta = {r["nuccore_acc"]: r for r in read_csv(Path(args.plasmids), "\t")} if args.plasmids else {}
    rows = []
    for bs, records in sorted(accessions.items()):
        accs = [r["nuccore_acc"] for r in records]
        fasta = client.get(f"{NCBI_EUTILS}/efetch.fcgi", {"db": "nuccore", "id": ",".join(sorted(accs)),
                                                          "rettype": "fasta", "retmode": "text"})
        if not fasta.startswith(">"):
            log.warning("no sequences returned for %s", bs)
            continue
        (outdir / f"{bs}.fasta").write_text(fasta)
        for rec in sorted(records, key=lambda r: r["nuccore_acc"]):
            carbapenemases = meta.get(rec["nuccore_acc"], {}).get("carbapenemases", "")
            rows.append({"biosample": bs, **rec, "carbapenemases": carbapenemases,
                         "in_selection": "yes" if rec["nuccore_acc"] in meta else "no"})
    write_tsv(outdir / "truth.tsv", rows, ["biosample", "nuccore_acc", "length", "topology", "completeness",
                                           "carbapenemases", "in_selection"])
    log.info("published plasmid sequences for %d BioSamples (%d plasmids) in %s",
             len(accessions), sum(map(len, accessions.values())), outdir)


# ---------------------------------------------------------------- common


def write_tsv(path, rows, fields) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fields), delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--cache", default=".cache", help="where API responses are cached (default .cache)")
    p.add_argument("-v", "--verbose", action="store_true")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("select", help="choose PLSDB plasmids")
    s.add_argument("--plsdb", required=True, help="directory with PLSDB's metadata (nuccore.csv, amr.tsv, biosample.csv)")
    s.add_argument("-o", "--output", required=True)
    s.add_argument("--others", type=int, default=0, help="also include this many random plasmids without a carbapenemase")
    s.add_argument("--all-topologies", dest="complete_circular", action="store_false",
                   help="keep linear and incomplete records too (default: complete circular only)")
    s.add_argument("--seed", type=int, default=1)
    s.set_defaults(func=cmd_select)

    s = sub.add_parser("runs", help="find read runs and nanopore chemistry for each BioSample")
    s.add_argument("-i", "--input", required=True, help="plasmids.tsv from select")
    s.add_argument("-o", "--output", required=True)
    s.add_argument("--sample", type=int, default=0, help="only this many random BioSamples (default: all)")
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--threads", type=int, default=4)
    s.add_argument("--peek-quality", type=int, default=0, metavar="READS",
                   help="for nanopore runs of unknown chemistry, stream this many reads (e.g. 2000) and record their "
                        "median quality as a hint (default: off)")
    s.add_argument("--peek-all", action="store_true", help="with --peek-quality, do it for every nanopore run")
    s.set_defaults(func=cmd_runs)

    s = sub.add_parser("subset", help="choose BioSamples to download")
    s.add_argument("-i", "--input", required=True, help="runs.tsv from runs")
    s.add_argument("-o", "--output", required=True)
    s.add_argument("--chemistry", default="R10.4.1", help="nanopore chemistry to require, or 'any' (default R10.4.1)")
    s.add_argument("--no-illumina", dest="require_illumina", action="store_false",
                   help="don't require paired Illumina reads as well")
    s.add_argument("--unknown-min-q", type=float, default=0, metavar="Q",
                   help="also accept nanopore runs of unknown chemistry whose median read quality (runs "
                        "--peek-quality) is at least Q, e.g. 16 (default: off)")
    s.add_argument("--max-biosamples", type=int, default=0)
    s.add_argument("--max-gb", type=float, default=0, help="total size budget for the reads (gzipped)")
    s.add_argument("--max-run-gb", type=float, default=0, help="skip BioSamples whose reads are bigger than this")
    s.add_argument("--seed", type=int, default=1)
    s.set_defaults(func=cmd_subset)

    s = sub.add_parser("download", help="download reads and published sequences for a subset")
    s.add_argument("-i", "--input", required=True, help="subset.tsv from subset")
    s.add_argument("-o", "--outdir", required=True)
    s.add_argument("--plsdb", required=True, help="PLSDB metadata directory, to find every plasmid of each BioSample")
    s.add_argument("--plasmids", help="plasmids.tsv from select, to annotate the truth set")
    s.add_argument("--skip-reads", action="store_true", help="only fetch the published sequences")
    s.add_argument("--illumina-source", choices=("ncbi", "ena"), default="ncbi",
                   help="ncbi (default): NCBI's full-quality copy, converted with fasterq-dump; ena: ENA's "
                        "FASTQ, which for many SRA runs has placeholder qualities")
    s.add_argument("--fasterq-dump", default="fasterq-dump", help="fasterq-dump to use (from sra-tools)")
    s.add_argument("--threads", type=int, default=4)
    s.set_defaults(func=cmd_download)

    args = p.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")
    args.func(args)


if __name__ == "__main__":
    main()
