#!/usr/bin/env python3
"""GEO (Gene Expression Omnibus) retrieval skill.

NCBI GEO is the largest public repository of high-throughput functional
genomics data. This skill wraps the entry points researchers actually
need:

  1. **Search** GEO for studies (Series / GSE) matching keywords,
     organism, platform, or study type via NCBI E-utilities.
  2. **Inspect** a single GSE — pulling its SOFT-format metadata,
     sample sheet (GSM accessions + characteristics), platform info,
     and raw / processed file inventory.
  3. **Download** processed expression matrices, SOFT files, and
     supplementary archives from GEO's FTP mirror without leaving the
     terminal.
  4. **Hand off** to the RNA-seq analysis skill: emit a sample-sheet
     CSV + counts-matrix-path so `igvfagent rnaseq pipeline` can pick
     up where this skill left off.

Subcommands

  search          Keyword / organism / platform search (E-utilities)
  series          Pull metadata + sample sheet for one GSE
  download        Pull supplementary files for a GSE
  list-files      List downloadable artifacts under a GSE without
                  fetching them
  sample-sheet    Convert a GSE's GSM list into a CSV ready to feed
                  into the rnaseq skill
  write-playbook  Emit Docs/Skills/GEO_RETRIEVAL_SKILLS.md

Endpoints used (resolved through ``Scripts/_endpoints.py``):
  - NCBI E-utilities: esearch / esummary / efetch on the ``gds``
    database (GEO datasets/series).
  - GEO FTP: ``ftp.ncbi.nlm.nih.gov/geo/series/GSEnnnnnn/<GSE>/...`` —
    sigh-aware, GEO splits series into thousand-buckets:
    ``geo/series/GSE9nnn/GSE9574/`` for GSE9574.

Outputs land under:
  Data/IGVF/GEO/Metadata/<GSE>.json     hydrated metadata
  Data/IGVF/GEO/Downloads/<GSE>/        downloaded files
  Data/Manifests/GEO/                   manifests + sample sheets
  Docs/GEO/                             reports + plots
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import logging
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import OrderedDict
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "GEO"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
MANIFEST_DIR = DATA_DIR / "Manifests" / "GEO"
METADATA_DIR = DATA_DIR / "IGVF" / "GEO" / "Metadata"
DOWNLOAD_DIR = DATA_DIR / "IGVF" / "GEO" / "Downloads"

EUTILS_BASE = _resolve_endpoint("pubmed_eutils", "PUBMED_EUTILS_BASE")
GEO_FTP = _resolve_endpoint("geo_ftp", "GEO_FTP_BASE")

USER_AGENT = "IGVFdataAgent-GEO/0.1"

# GSE files live in thousand-buckets at GEO's FTP mirror — e.g. GSE9574
# is under series/GSE9nnn/GSE9574/. Build the prefix from the accession.
def _gse_ftp_prefix(gse: str) -> str:
    """e.g. 'GSE9574' -> 'series/GSE9nnn/GSE9574'."""
    if not gse.startswith("GSE"):
        raise SystemExit(f"Bad GSE accession: {gse!r}")
    digits = gse[3:]
    if not digits.isdigit():
        raise SystemExit(f"Bad GSE accession: {gse!r}")
    bucket = "GSE" + digits[:-3] + "nnn" if len(digits) > 3 else "GSE0nnn"
    return f"series/{bucket}/{gse}"


# --------------------------- Project plumbing ------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"geo_retrieval_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    logging.info("Log file: %s", log_path)
    return log_path


def mkdirs() -> None:
    for d in (REPORT_DIR, MANIFEST_DIR, METADATA_DIR, DOWNLOAD_DIR,
              SKILL_DOC_DIR):
        d.mkdir(parents=True, exist_ok=True)


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


# --------------------------- HTTP helpers ----------------------------------

def http_get(url: str, params: dict | None = None,
             timeout: int = 60) -> "tuple[int, str]":
    if params:
        sep = "&" if "?" in url else "?"
        url = (url + sep
                + urllib.parse.urlencode({k: v for k, v in params.items()
                                            if v is not None}))
    logging.info("GET %s", url)
    req = urllib.request.Request(url,
                                    headers={"User-Agent": USER_AGENT,
                                             "Accept": "*/*"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status, resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        return e.code, (e.read().decode("utf-8", errors="replace")
                          if e.fp else "")
    except Exception as e:
        return 0, f"network_error: {e}"


def download_file(url: str, dest: Path) -> int:
    logging.info("Download %s -> %s", url, dest)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    written = 0
    with urllib.request.urlopen(req, timeout=300) as resp, \
         dest.open("wb") as f:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            f.write(chunk); written += len(chunk)
    return written


def write_csv(path: Path, rows: "list[dict]",
                cols: "Optional[list[str]]" = None) -> None:
    if not rows:
        path.write_text("")
        return
    cols = cols or sorted({k for r in rows for k in r})
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# --------------------------- E-utilities -----------------------------------

def search_geo(query: str, *, organism: str | None = None,
                  platform: str | None = None,
                  study_type: str | None = None,
                  limit: int = 25) -> "list[dict]":
    """Search NCBI gds (GEO DataSets / Series) and return summaries."""
    parts = [f"({query})"]
    if organism:
        parts.append(f'"{organism}"[Organism]')
    if platform:
        parts.append(f'"{platform}"[Platform]')
    if study_type:
        parts.append(f'"{study_type}"[DataSet Type]')
    parts.append('"GSE"[Entry Type]')
    full = " AND ".join(parts)

    sc, body = http_get(f"{EUTILS_BASE}/esearch.fcgi", params={
        "db": "gds", "term": full, "retmax": limit,
        "retmode": "json", "sort": "relevance",
    })
    if sc != 200:
        logging.warning("GEO esearch %s", sc)
        return []
    try:
        ids = (json.loads(body).get("esearchresult", {})
                 .get("idlist", []) or [])
    except Exception:
        return []
    if not ids:
        return []
    sc, body = http_get(f"{EUTILS_BASE}/esummary.fcgi", params={
        "db": "gds", "id": ",".join(ids), "retmode": "json",
    })
    if sc != 200:
        return []
    try:
        data = json.loads(body).get("result", {})
    except Exception:
        return []
    out: "list[dict]" = []
    for uid in ids:
        s = data.get(uid)
        if not isinstance(s, dict):
            continue
        accession = s.get("accession", "") or s.get("entrytype", "")
        if not accession.startswith("GSE"):
            # Other entry types (GDS, GPL) — skip; we want Series here
            continue
        out.append({
            "uid":          uid,
            "accession":    accession,
            "title":        s.get("title", ""),
            "summary":      (s.get("summary") or "")[:280],
            "n_samples":    s.get("n_samples", ""),
            "platform":     ", ".join(p.get("accession", "")
                                        for p in (s.get("platforms") or [])),
            "organism":     ", ".join(s.get("taxon", []) or [])
                              if isinstance(s.get("taxon"), list)
                              else str(s.get("taxon", "")),
            "pdat":         s.get("pdat", ""),
            "url":          f"https://www.ncbi.nlm.nih.gov/geo/query/"
                              f"acc.cgi?acc={accession}",
        })
    return out


# --------------------------- SOFT parser -----------------------------------

_SOFT_KEY_RE = re.compile(r"^!([A-Za-z0-9_]+)\s*=\s*(.*)$")


def _parse_soft(text: str) -> "dict[str, dict[str, list[str]]]":
    """Parse a SOFT-format response into per-entity attribute maps.

    Returns a dict keyed by entity ID (e.g. ``GSE9574``, ``GSM242540``,
    ``GPL96``) with values being maps from field name (``title``,
    ``summary``, ``characteristics_ch1`` …) to a list of values."""
    entities: "dict[str, dict[str, list[str]]]" = OrderedDict()
    current: "Optional[str]" = None
    for raw_line in text.splitlines():
        line = raw_line.strip("\r\n")
        if not line:
            continue
        if line.startswith("^"):
            try:
                _, _, value = line.partition("=")
            except ValueError:
                continue
            current = value.strip().split()[0] if value.strip() else None
            if current:
                entities.setdefault(current, OrderedDict())
            continue
        if not current:
            continue
        m = _SOFT_KEY_RE.match(line)
        if not m:
            continue
        key, val = m.group(1).lower(), m.group(2).strip()
        entities[current].setdefault(key, []).append(val)
    return entities


def fetch_soft(gse: str, view: str = "quick") -> "dict[str, dict[str, list[str]]]":
    """Pull the GEO SOFT (text) record for a Series. ``view='quick'``
    returns just the GSE record; ``view='full'`` includes all GSM and
    GPL sub-entities (much larger payload)."""
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    params = {"acc": gse, "targ": "self" if view == "quick" else "all",
                "form": "text", "view": view}
    sc, body = http_get(url, params=params, timeout=120)
    if sc != 200:
        logging.warning("GEO SOFT %s for %s", sc, gse)
        return {}
    return _parse_soft(body)


def summarize_series(soft_quick: "dict[str, dict]",
                       gse: str) -> dict:
    """Pull a one-page summary out of a quick-view SOFT record."""
    rec = soft_quick.get(gse, {})
    return {
        "accession":  gse,
        "title":      "; ".join(rec.get("series_title", [])),
        "summary":    " ".join(rec.get("series_summary", []))[:600],
        "overall_design":
            " ".join(rec.get("series_overall_design", []))[:600],
        "type":       "; ".join(rec.get("series_type", [])),
        "submission_date":
            "; ".join(rec.get("series_submission_date", [])),
        "last_update_date":
            "; ".join(rec.get("series_last_update_date", [])),
        "platforms":  "; ".join(rec.get("series_platform_id", [])),
        "samples":    "; ".join(rec.get("series_sample_id", [])),
        "n_samples":  len(rec.get("series_sample_id", [])),
        "pubmed_id":  "; ".join(rec.get("series_pubmed_id", [])),
        "contact":    "; ".join(
            rec.get("series_contact_name", [])
            + rec.get("series_contact_institute", [])
        ),
        "url":        f"https://www.ncbi.nlm.nih.gov/geo/query/"
                        f"acc.cgi?acc={gse}",
    }


def parse_sample_sheet(soft_full: "dict[str, dict]") -> "list[dict]":
    """Extract one row per GSM with its characteristics_ch1 fields
    expanded into columns (e.g. ``cell line: GM12878`` →
    ``cell_line=GM12878``)."""
    rows: "list[dict]" = []
    for entity_id, fields in soft_full.items():
        if not entity_id.startswith("GSM"):
            continue
        row = {
            "gsm":            entity_id,
            "title":          "; ".join(fields.get("sample_title", [])),
            "source_name":    "; ".join(fields.get("sample_source_name_ch1", [])),
            "organism":       "; ".join(fields.get("sample_organism_ch1", [])),
            "molecule":       "; ".join(fields.get("sample_molecule_ch1", [])),
            "library_strategy":
                "; ".join(fields.get("sample_library_strategy", [])),
            "library_source":
                "; ".join(fields.get("sample_library_source", [])),
            "library_selection":
                "; ".join(fields.get("sample_library_selection", [])),
            "instrument":     "; ".join(fields.get("sample_instrument_model", [])),
            "platform_id":    "; ".join(fields.get("sample_platform_id", [])),
        }
        for c in fields.get("sample_characteristics_ch1", []):
            if ":" not in c:
                continue
            k, _, v = c.partition(":")
            key = "char_" + safe_label(k.strip().lower())
            row[key] = v.strip()
        # SRA / supplementary URLs
        sra = [u for u in fields.get("sample_supplementary_file", [])
               if "sra" in u.lower() or u.endswith((".sra", ".fastq",
                                                       ".fastq.gz"))]
        if sra:
            row["sra_files"] = "; ".join(sra)
        rows.append(row)
    return rows


# --------------------------- FTP file listing ------------------------------

# GEO's FTP server runs an Apache index at the HTTPS mirror, so we can
# scrape filenames from a directory listing's HTML. Cheap and robust.
_FILE_LINK_RE = re.compile(r'<a\s+href="([^"?][^"]*)"', re.I)


def list_ftp_dir(rel_path: str) -> "list[dict]":
    """List files in a GEO FTP subdirectory (e.g.
    'series/GSE9nnn/GSE9574/suppl/'). Returns name + size (best effort)."""
    url = f"{GEO_FTP}/{rel_path.strip('/')}/"
    sc, body = http_get(url, timeout=60)
    if sc != 200:
        return []
    out: "list[dict]" = []
    seen: "set[str]" = set()
    for m in _FILE_LINK_RE.finditer(body):
        name = m.group(1).rstrip("/")
        if not name or name in seen or name.startswith("?")\
                or name == ".." or "/" in name:
            continue
        seen.add(name)
        if name in (".", ".."):
            continue
        out.append({
            "name": name,
            "url":  url + name,
            "kind": "dir" if m.group(1).endswith("/") else "file",
        })
    return out


def list_gse_files(gse: str) -> "list[dict]":
    """Scrape the supplementary + matrix subfolders of a GSE."""
    prefix = _gse_ftp_prefix(gse)
    out: "list[dict]" = []
    for sub in ("matrix", "suppl", "soft"):
        rows = list_ftp_dir(f"{prefix}/{sub}")
        for r in rows:
            r["category"] = sub
            r["gse"] = gse
            out.append(r)
    return out


# --------------------------- Subcommands -----------------------------------

def cmd_search(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    rows = search_geo(args.query, organism=args.organism,
                         platform=args.platform,
                         study_type=args.study_type, limit=args.limit)
    ts = timestamp()
    label = safe_label(args.label or args.query)[:60]
    out = MANIFEST_DIR / f"{ts}_search_{label}.csv"
    write_csv(out, rows,
              cols=["accession", "title", "n_samples", "organism",
                     "platform", "pdat", "summary", "url"])
    print(f"Found {len(rows)} GSE matches.")
    print(f"Manifest: {out}")
    if rows:
        print()
        for r in rows[: min(10, len(rows))]:
            print(f"  {r['accession']:>9}  {r['title'][:80]}")
    return out


def cmd_series(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    gse = args.gse.strip()
    label = safe_label(args.label or gse)
    full = fetch_soft(gse, view="full" if args.full_samples else "quick")
    if not full:
        raise SystemExit(f"No SOFT response for {gse}.")
    quick = full if args.full_samples else fetch_soft(gse, view="quick")
    summary = summarize_series(quick, gse)
    (METADATA_DIR / f"{gse}.json").write_text(
        json.dumps({"summary": summary, "soft_keys":
                       {k: list(v.keys()) for k, v in full.items()}},
                     indent=2, default=str))

    samples = parse_sample_sheet(full) if args.full_samples else []
    if samples:
        sheet_path = MANIFEST_DIR / f"{ts}_{label}_samples.csv"
        write_csv(sheet_path, samples)
        print(f"Sample sheet: {sheet_path} ({len(samples)} GSMs)")

    files = list_gse_files(gse)
    files_path = MANIFEST_DIR / f"{ts}_{label}_files.csv"
    write_csv(files_path, files,
              cols=["gse", "category", "name", "url", "kind"])

    report = REPORT_DIR / f"{ts}_{label}_geo_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# GEO Series — `{gse}`", "",
             f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}", ""]
    for k, v in summary.items():
        lines.append(f"- **{k}**: {v}")
    lines += ["", f"## Files ({len(files)})", "",
              "| category | name | url |", "|---|---|---|"]
    for f in files[:50]:
        lines.append(f"| {f['category']} | {f['name']} | {f['url']} |")
    if len(files) > 50:
        lines.append(f"| _… and {len(files) - 50} more in the CSV_ | | |")
    if samples:
        lines += ["", f"## Samples (first 10 of {len(samples)})", "",
                   "| GSM | title | source | strategy |",
                   "|---|---|---|---|"]
        for s in samples[:10]:
            lines.append(f"| {s['gsm']} | {s['title'][:60]} | "
                          f"{s['source_name'][:60]} | "
                          f"{s.get('library_strategy','')} |")
    report.write_text("\n".join(lines))
    print(f"Report:  {report}")
    print(f"Files:   {files_path}")
    print(f"Summary: {summary['title']}")
    return report


def cmd_list_files(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    files = list_gse_files(args.gse.strip())
    out = MANIFEST_DIR / f"{ts}_{safe_label(args.gse)}_files.csv"
    write_csv(out, files, cols=["gse", "category", "name", "url", "kind"])
    print(f"Found {len(files)} files for {args.gse}")
    print(f"Manifest: {out}")
    return out


def cmd_download(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    gse = args.gse.strip()
    out_dir = DOWNLOAD_DIR / gse
    out_dir.mkdir(parents=True, exist_ok=True)
    files = list_gse_files(gse)
    if args.only:
        only = {x.lower() for x in args.only}
        files = [f for f in files if f["category"] in only]
    if args.pattern:
        rx = re.compile(args.pattern, re.I)
        files = [f for f in files if rx.search(f["name"])]
    files = [f for f in files if f["kind"] != "dir"]
    cap = int(args.max_download_gb * (1024 ** 3)) if args.max_download_gb else None
    used = 0
    log: "list[dict]" = []
    for f in files:
        dest = out_dir / f["name"]
        if dest.exists() and dest.stat().st_size > 0:
            log.append({**f, "downloaded_path": str(dest),
                          "downloaded_bytes": dest.stat().st_size,
                          "status": "cached"})
            continue
        try:
            n = download_file(f["url"], dest)
        except Exception as e:
            log.append({**f, "status": f"failed: {e}"})
            continue
        used += n
        log.append({**f, "downloaded_path": str(dest),
                      "downloaded_bytes": n, "status": "ok"})
        if cap and used > cap:
            logging.info("Download cap reached.")
            break
    log_path = MANIFEST_DIR / f"{timestamp()}_{safe_label(gse)}_download_log.csv"
    write_csv(log_path, log)
    print(f"Downloaded {sum(1 for r in log if r.get('status') == 'ok')} new "
          f"files / {used / (1024**3):.2f} GB into {out_dir}")
    print(f"Log: {log_path}")
    return log_path


def cmd_sample_sheet(args: argparse.Namespace) -> Path:
    """Fetch full SOFT for a GSE, then emit a clean sample-sheet CSV
    with one row per GSM expanded by characteristics."""
    setup_logging(); mkdirs()
    full = fetch_soft(args.gse.strip(), view="full")
    samples = parse_sample_sheet(full)
    if not samples:
        raise SystemExit(f"No samples parsed from SOFT for {args.gse}")
    ts = timestamp()
    label = safe_label(args.label or args.gse)
    out = MANIFEST_DIR / f"{ts}_{label}_sample_sheet.csv"
    write_csv(out, samples)
    print(f"Sample sheet: {out}  ({len(samples)} GSMs)")
    return out


def cmd_write_playbook(_a) -> Path:
    mkdirs()
    path = SKILL_DOC_DIR / "GEO_RETRIEVAL_SKILLS.md"
    lines = [
        "# Skill: GEO retrieval",
        "",
        "Wraps NCBI E-utilities + the GEO FTP mirror so the agent can "
        "search GEO, pull a Series' metadata + sample sheet, and "
        "download supplementary expression matrices from the terminal.",
        "",
        "## Subcommands",
        "",
        "### `search` — keyword / organism / platform search",
        "",
        "```bash",
        "igvfagent geo search --query 'GM12878 RNA-seq' \\",
        "    --organism 'Homo sapiens' --limit 25",
        "igvfagent geo search --query 'lymphoblastoid H3K27ac' --study-type 'Expression profiling by high throughput sequencing'",
        "```",
        "",
        "### `series` — pull metadata + file inventory for a GSE",
        "",
        "```bash",
        "igvfagent geo series --gse GSE9574 --full-samples",
        "```",
        "",
        "Writes a markdown report under `Docs/GEO/<ts>_<label>_geo_report.md`, "
        "a sample sheet (`<label>_samples.csv`), and a file inventory "
        "(`<label>_files.csv`).",
        "",
        "### `list-files` — what's available on GEO FTP for a GSE",
        "",
        "```bash",
        "igvfagent geo list-files --gse GSE9574",
        "```",
        "",
        "### `download` — pull supplementary / matrix files",
        "",
        "```bash",
        "igvfagent geo download --gse GSE9574 --only matrix --max-download-gb 1",
        "igvfagent geo download --gse GSE9574 --pattern 'counts.*tsv|tpm.*tsv' --max-download-gb 2",
        "```",
        "",
        "By default the ``download`` step pulls everything under "
        "``series/GSEXnnn/GSEX/{matrix,suppl,soft}/``; ``--only`` "
        "restricts to one of those categories and ``--pattern`` is a "
        "case-insensitive regex applied to the filename.",
        "",
        "### `sample-sheet` — clean sample CSV for the rnaseq skill",
        "",
        "```bash",
        "igvfagent geo sample-sheet --gse GSE9574 --label gse9574_sheet",
        "```",
        "",
        "Reads the SOFT record's full sample block and produces a CSV "
        "with one row per GSM, columns expanded from "
        "``characteristics_ch1`` (cell line, treatment, etc.). The output "
        "feeds straight into ``igvfagent rnaseq pipeline --sample-sheet "
        "<csv>``.",
        "",
        "## How this chains with other skills",
        "",
        "- After `series` + `download`, hand the sample sheet to "
        "`igvfagent rnaseq pipeline` for QC + DEG + plotting.",
        "- The DEG list from rnaseq can then be cross-referenced with "
        "`igvfagent kg variant` / `igvfagent encode integrate-ccre` for "
        "the controlling regulatory elements.",
        "- For literature corroboration of any GSE record, "
        "`igvfagent ref validate --input <degs.csv>`.",
        "",
        "## Outputs",
        "",
        "- Metadata JSON: `Data/IGVF/GEO/Metadata/<GSE>.json`",
        "- Downloads:     `Data/IGVF/GEO/Downloads/<GSE>/`",
        "- Manifests:     `Data/Manifests/GEO/<ts>_<label>_*.csv`",
        "- Reports:       `Docs/GEO/<ts>_<label>_geo_report.md`",
    ]
    path.write_text("\n".join(lines))
    print(f"Playbook: {path}")
    return path


# --------------------------------- CLI -------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(
        description="GEO (Gene Expression Omnibus) retrieval skill — "
                    "search NCBI for studies, pull metadata + sample "
                    "sheet, download supplementary expression files.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("search", help="Keyword search of GEO Series.")
    s.add_argument("--query", required=True)
    s.add_argument("--organism", default=None,
                    help='e.g. "Homo sapiens", "Mus musculus".')
    s.add_argument("--platform", default=None, help='e.g. "GPL11154".')
    s.add_argument("--study-type", default=None,
                    help='e.g. "Expression profiling by high throughput '
                         'sequencing".')
    s.add_argument("--limit", type=int, default=25)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_search)

    s = sub.add_parser("series",
                        help="Pull metadata + samples for one GSE.")
    s.add_argument("--gse", required=True)
    s.add_argument("--full-samples", action="store_true",
                    help="Also fetch every GSM's SOFT block (slower).")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_series)

    s = sub.add_parser("list-files",
                        help="Scrape the FTP listing for a GSE.")
    s.add_argument("--gse", required=True)
    s.set_defaults(func=cmd_list_files)

    s = sub.add_parser("download",
                        help="Download supplementary / matrix files.")
    s.add_argument("--gse", required=True)
    s.add_argument("--only", nargs="*", default=None,
                    choices=["matrix", "suppl", "soft"])
    s.add_argument("--pattern", default=None,
                    help="Regex over filenames (case-insensitive).")
    s.add_argument("--max-download-gb", type=float, default=None)
    s.set_defaults(func=cmd_download)

    s = sub.add_parser("sample-sheet",
                        help="Emit a clean sample CSV for the rnaseq skill.")
    s.add_argument("--gse", required=True)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_sample_sheet)

    s = sub.add_parser("write-playbook",
                        help="Emit Docs/Skills/GEO_RETRIEVAL_SKILLS.md.")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
