#!/usr/bin/env python3
"""Cross-source survey + downloader for single-cell multiome data.

Surveys IGVF Portal, ENCODE Portal, and GEO for **10x Multiome / SHARE-seq /
single-nucleus multiome** datasets, normalizes their metadata into a unified
manifest, and lets you download the training-relevant files (cell-by-gene /
cell-by-peak matrices, ATAC fragments, cell annotations) with a size cap.

Subcommands
-----------

  survey-igvf      Pull IGVF Portal AnalysisSets/MeasurementSets tagged
                   ``10x multiome`` (and SHARE-seq).
  survey-encode    Search ENCODE for multiome / SHARE-seq experiments.
  survey-geo       Keyword-search GEO (E-utilities) for multiome studies.
  survey-all       Run all three, write per-source + unified reports.
  manifest         Materialize a unified downloadable-file CSV from the
                   most recent surveys, classified by file kind.
  download         Download files from a manifest with size cap.
  inventory        Scan ``Data/MultiomeSurvey/`` and write a fresh local
                   inventory CSV.
  write-playbook   Emit ``Docs/Skills/MULTIOME_SURVEY_SKILLS.md``.

All endpoints are resolved via ``Scripts/_endpoints.py``; no URLs or
credentials appear in source.  Outputs use repo-relative paths and
runtime caches live under ``Data/MultiomeSurvey/`` which is gitignored.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Iterable

# Resolve service endpoints via the shared helper (kept out of source).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # type: ignore  # noqa: E402

# ---------------------------------------------------------------------------
# Paths and endpoints
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"

SURVEY_DIR = DATA_DIR / "MultiomeSurvey"
DOWNLOAD_DIR = SURVEY_DIR
MANIFEST_DIR = DATA_DIR / "Manifests" / "MultiomeSurvey"
REPORT_DIR = DOCS_DIR / "MultiomeSurvey"
SKILL_DOC = DOCS_DIR / "Skills" / "MULTIOME_SURVEY_SKILLS.md"
INVENTORY_CSV = SURVEY_DIR / "inventory.csv"

PORTAL_API_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_API_BASE")
ENCODE_BASE = _resolve_endpoint("encode", "ENCODE_BASE")
EUTILS_BASE = _resolve_endpoint("pubmed_eutils", "PUBMED_EUTILS_BASE")
GEO_FTP_BASE = _resolve_endpoint("geo_ftp", "GEO_FTP_BASE")

USER_AGENT = "IGVFagent-multiome-survey/0.1"

# ---------------------------------------------------------------------------
# Survey-time defaults
# ---------------------------------------------------------------------------

# Assay-title variants we treat as "multiome" across the three sources.
IGVF_PREFERRED_ASSAY_TITLES = [
    "10x multiome",
    "10x multiome with MULTI-seq",
]
IGVF_SHARESEQ_TITLES = [
    "SHARE-seq",
]
ENCODE_MULTIOME_SEARCH_TERMS = [
    "10x multiome",
    "multiome ATAC + Gene Expression",
    "SHARE-seq",
]
GEO_MULTIOME_QUERIES = [
    '"10x multiome"',
    '"10x Genomics multiome"',
    '"single-nucleus multiome"',
    '"single-cell multiome"',
    '"SHARE-seq"',
    '"Multiome ATAC + Gene Expression"',
]

# Map content_type / file_format → training-relevant "kind" used in the
# unified manifest.  Anything not matched is left as ``other``.
KIND_RULES: list[tuple[str, re.Pattern[str]]] = [
    ("matrix_rna", re.compile(r"sparse gene count matrix|gene quantifications|filtered feature barcode matrix|raw feature barcode matrix", re.I)),
    ("matrix_atac", re.compile(r"annotated sparse peak count matrix|peak.*matrix|cell.by.peak", re.I)),
    ("fragments", re.compile(r"fragments?\b", re.I)),
    ("peaks", re.compile(r"^peaks?$|peak calls?", re.I)),
    ("annotations", re.compile(r"cell annotations?|cell.metadata|sample.sheet|metadata\.tsv", re.I)),
    ("alignments", re.compile(r"alignments?|bam\b", re.I)),
    ("index", re.compile(r"\bindex\b|\.bai|\.tbi|\.crai", re.I)),
    ("raw_reads", re.compile(r"reads|fastq", re.I)),
]


def classify_kind(content_type: str | None, file_format: str | None, filename: str | None = None) -> str:
    """Map a (content_type, file_format, filename) to a kind label."""
    haystack = " ".join(str(x) for x in (content_type, file_format, filename) if x)
    for kind, rx in KIND_RULES:
        if rx.search(haystack):
            return kind
    return "other"


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"multiome_survey_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)
    return log_path


def mkdirs() -> None:
    for d in (SURVEY_DIR, MANIFEST_DIR, REPORT_DIR, LOG_DIR):
        d.mkdir(parents=True, exist_ok=True)


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in s)


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def rel_path(p: Path | str) -> str:
    pp = Path(p).resolve()
    try:
        return pp.relative_to(ROOT).as_posix()
    except ValueError:
        return pp.as_posix()


def abs_from_rel(rel: str) -> Path:
    p = Path(rel)
    return p if p.is_absolute() else (ROOT / p)


def http_get(url: str, params: dict[str, Any] | None = None, *,
             headers: dict[str, str] | None = None,
             timeout: int = 60, retries: int = 3) -> tuple[int, bytes]:
    if params:
        url = f"{url}?{urllib.parse.urlencode(params, doseq=True)}"
    base_headers = {"User-Agent": USER_AGENT, "Accept": "application/json,*/*"}
    if headers:
        base_headers.update(headers)
    last_exc: Exception | None = None
    for attempt in range(retries):
        req = urllib.request.Request(url, headers=base_headers, method="GET")
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as exc:
            body = exc.read() if hasattr(exc, "read") else b""
            return exc.code, body
        except urllib.error.URLError as exc:
            last_exc = exc
            wait = 1.5 * (attempt + 1)
            logging.warning("network error on %s: %s (retry in %.1fs)", url, exc.reason, wait)
            time.sleep(wait)
    if last_exc is not None:
        logging.error("giving up on %s after %d retries: %s", url, retries, last_exc)
    return 0, b""


def http_get_json(url: str, params: dict[str, Any] | None = None, *,
                  headers: dict[str, str] | None = None, timeout: int = 60) -> Any:
    status, body = http_get(url, params, headers=headers, timeout=timeout)
    if status == 0 or status >= 400 or not body:
        return None
    try:
        return json.loads(body)
    except json.JSONDecodeError:
        return None


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def stringify(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ";".join(str(v) for v in value if v is not None)
    return str(value)


# ---------------------------------------------------------------------------
# Unified record schema
# ---------------------------------------------------------------------------

DATASET_FIELDS = [
    "source",            # igvf | encode | geo
    "dataset_accession", # IGVFDS… / ENCSR… / GSE…
    "title",
    "summary",
    "organism",
    "assay",
    "samples",           # joined sample/tissue/cell-line summary
    "n_samples",
    "publications",
    "url",               # human-browsable URL
    "lab_or_pi",
    "creation_date",
]

FILE_FIELDS = [
    "source",
    "dataset_accession",
    "file_accession",    # IGVFFI… / ENCFF… / GSM… (or filename for GEO suppl)
    "kind",              # classify_kind() output
    "content_type",
    "file_format",
    "file_size_bytes",
    "download_url",
    "local_path",        # repo-relative target
    "md5sum",
    "notes",
]


# ---------------------------------------------------------------------------
# IGVF Portal survey
# ---------------------------------------------------------------------------


def igvf_search(params: dict[str, Any], *, limit: int = 200) -> list[dict[str, Any]]:
    """Query IGVF /search/ returning @graph rows.  Uses frame=object for rich fields."""
    p = dict(params)
    p.setdefault("format", "json")
    p.setdefault("frame", "object")
    p["limit"] = str(limit)
    payload = http_get_json(f"{PORTAL_API_BASE}/search/", p)
    if not isinstance(payload, dict):
        return []
    rows = payload.get("@graph") or []
    return [r for r in rows if isinstance(r, dict)]


def igvf_object(rel_path_value: str) -> dict[str, Any] | None:
    """Fetch a single IGVF object by its @id path."""
    if not rel_path_value:
        return None
    path = rel_path_value if rel_path_value.startswith("/") else f"/{rel_path_value}"
    payload = http_get_json(f"{PORTAL_API_BASE}{path}", {"format": "json", "frame": "object"})
    return payload if isinstance(payload, dict) else None


def survey_igvf(*, limit: int = 200, label: str = "igvf_multiome",
                fetch_files: bool = False,
                include_shareseq: bool = True) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Return (datasets, files) for IGVF multiome AnalysisSets."""
    titles = list(IGVF_PREFERRED_ASSAY_TITLES)
    if include_shareseq:
        titles += IGVF_SHARESEQ_TITLES

    seen: dict[str, dict[str, Any]] = {}
    all_files: list[dict[str, Any]] = []

    for title in titles:
        rows = igvf_search(
            {"type": "AnalysisSet", "preferred_assay_titles": title}, limit=limit,
        )
        logging.info("IGVF AnalysisSet  preferred_assay_titles=%-30s rows=%d", title, len(rows))
        for row in rows:
            acc = row.get("accession")
            if not acc or acc in seen:
                continue
            samples = row.get("samples") or []
            sample_summaries = []
            for s in samples[:3]:
                # samples can be either dicts or @id strings depending on frame
                if isinstance(s, dict):
                    sample_summaries.append(s.get("summary") or s.get("accession") or "")
                else:
                    sample_summaries.append(str(s))
            rec = {
                "source": "igvf",
                "dataset_accession": acc,
                "title": row.get("summary") or row.get("description") or "",
                "summary": row.get("summary") or "",
                "organism": stringify(row.get("organism") or row.get("taxa")),
                "assay": stringify(row.get("preferred_assay_titles") or row.get("assay_titles")),
                "samples": stringify(sample_summaries),
                "n_samples": len(samples),
                "publications": stringify(row.get("publications")),
                "url": f"https://data.igvf.org/analysis-sets/{acc}/",
                "lab_or_pi": stringify(row.get("lab")),
                "creation_date": row.get("creation_timestamp", "") or row.get("release_timestamp", ""),
            }
            seen[acc] = rec

            if fetch_files:
                file_paths = row.get("files") or []
                for fp in file_paths:
                    # file references in @graph are @id strings
                    fpath = fp if isinstance(fp, str) else fp.get("@id", "")
                    obj = igvf_object(fpath)
                    if not obj:
                        continue
                    href = obj.get("href") or ""
                    download_url = f"{PORTAL_API_BASE}{href}" if href.startswith("/") else href
                    file_acc = obj.get("accession") or Path(href).stem
                    content_type = obj.get("content_type") or ""
                    file_format = obj.get("file_format") or ""
                    file_record = {
                        "source": "igvf",
                        "dataset_accession": acc,
                        "file_accession": file_acc,
                        "kind": classify_kind(content_type, file_format, href),
                        "content_type": content_type,
                        "file_format": file_format,
                        "file_size_bytes": obj.get("file_size") or 0,
                        "download_url": download_url,
                        "local_path": rel_path(DOWNLOAD_DIR / "igvf" / acc / f"{file_acc}.{file_format}"),
                        "md5sum": obj.get("md5sum", ""),
                        "notes": stringify(obj.get("content_summary")),
                    }
                    all_files.append(file_record)

    datasets = list(seen.values())
    logging.info("IGVF survey: %d unique AnalysisSets, %d files (fetch_files=%s)",
                 len(datasets), len(all_files), fetch_files)
    return datasets, all_files


# ---------------------------------------------------------------------------
# ENCODE survey
# ---------------------------------------------------------------------------


def encode_search(params: dict[str, Any], *, limit: int = 200) -> list[dict[str, Any]]:
    p = dict(params)
    p.setdefault("format", "json")
    p["limit"] = str(limit)
    payload = http_get_json(f"{ENCODE_BASE}/search/", p)
    if not isinstance(payload, dict):
        return []
    rows = payload.get("@graph") or []
    return [r for r in rows if isinstance(r, dict)]


def encode_object(at_id: str) -> dict[str, Any] | None:
    if not at_id:
        return None
    path = at_id if at_id.startswith("/") else f"/{at_id}"
    payload = http_get_json(f"{ENCODE_BASE}{path}", {"format": "json"})
    return payload if isinstance(payload, dict) else None


def survey_encode(*, limit: int = 200, label: str = "encode_multiome",
                  fetch_files: bool = False) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Return (datasets, files) for ENCODE multiome / SHARE-seq experiments + Series."""
    seen: dict[str, dict[str, Any]] = {}
    all_files: list[dict[str, Any]] = []

    # Multiome data on ENCODE shows up under both Experiment (single-cell
    # paired RNA + ATAC) and SingleCellRnaSeries / FunctionalCharacterizationSeries.
    for type_filter in ("Experiment", "SingleCellRnaSeries", "FunctionalCharacterizationSeries"):
        for term in ENCODE_MULTIOME_SEARCH_TERMS:
            rows = encode_search(
                {"type": type_filter, "searchTerm": term, "status": "released"},
                limit=limit,
            )
            logging.info("ENCODE type=%-30s searchTerm=%-40s rows=%d", type_filter, term, len(rows))
            for row in rows:
                acc = row.get("accession")
                if not acc or acc in seen:
                    continue
                bios = row.get("biosample_summary") or row.get("biosample_ontology") or ""
                rec = {
                    "source": "encode",
                    "dataset_accession": acc,
                    "title": row.get("description") or row.get("assay_title") or "",
                    "summary": row.get("biosample_summary") or "",
                    "organism": stringify(row.get("replicates_organism")
                                          or row.get("organism") or row.get("biosample_organism")),
                    "assay": stringify(row.get("assay_title") or row.get("assay_term_name")),
                    "samples": stringify(bios),
                    "n_samples": int(row.get("replication_type") is not None) + 1,
                    "publications": stringify(row.get("references")),
                    "url": f"{ENCODE_BASE}/experiments/{acc}/",
                    "lab_or_pi": stringify(row.get("lab")),
                    "creation_date": row.get("date_released") or row.get("date_created", ""),
                }
                seen[acc] = rec

                if fetch_files:
                    file_paths = row.get("files") or []
                    for fp in file_paths:
                        fid = fp if isinstance(fp, str) else fp.get("@id", "")
                        obj = encode_object(fid)
                        if not obj:
                            continue
                        href = obj.get("href") or ""
                        download_url = (
                            f"{ENCODE_BASE}{href}" if href.startswith("/") else href
                        )
                        file_acc = obj.get("accession") or Path(href).stem
                        content_type = obj.get("output_type") or obj.get("content_type") or ""
                        file_format = obj.get("file_format") or ""
                        file_record = {
                            "source": "encode",
                            "dataset_accession": acc,
                            "file_accession": file_acc,
                            "kind": classify_kind(content_type, file_format, href),
                            "content_type": content_type,
                            "file_format": file_format,
                            "file_size_bytes": obj.get("file_size") or 0,
                            "download_url": download_url,
                            "local_path": rel_path(DOWNLOAD_DIR / "encode" / acc / f"{file_acc}.{file_format}"),
                            "md5sum": obj.get("md5sum", ""),
                            "notes": stringify(obj.get("output_category")),
                        }
                        all_files.append(file_record)

    datasets = list(seen.values())
    logging.info("ENCODE survey: %d unique experiments, %d files (fetch_files=%s)",
                 len(datasets), len(all_files), fetch_files)
    return datasets, all_files


# ---------------------------------------------------------------------------
# GEO survey via NCBI E-utilities
# ---------------------------------------------------------------------------


def eutils_search_gds(query: str, *, retmax: int = 200,
                     organism: str | None = None) -> list[str]:
    """Run esearch.fcgi on the gds database; return GDS UIDs (not GSE IDs yet)."""
    term = query
    if organism:
        term += f' AND "{organism}"[Organism]'
    # Restrict to GEO Series records (entry_type=gse) — gives the parent study UIDs.
    term += " AND gse[Entry Type]"
    params = {
        "db": "gds",
        "term": term,
        "retmode": "json",
        "retmax": str(retmax),
    }
    payload = http_get_json(f"{EUTILS_BASE}/esearch.fcgi", params)
    if not isinstance(payload, dict):
        return []
    ids = payload.get("esearchresult", {}).get("idlist") or []
    return [str(i) for i in ids]


def eutils_summary_gds(uids: list[str]) -> list[dict[str, Any]]:
    """esummary in batches of 200 — returns GDS summary records."""
    out: list[dict[str, Any]] = []
    for i in range(0, len(uids), 200):
        chunk = uids[i:i + 200]
        params = {
            "db": "gds",
            "id": ",".join(chunk),
            "retmode": "json",
        }
        payload = http_get_json(f"{EUTILS_BASE}/esummary.fcgi", params)
        if not isinstance(payload, dict):
            continue
        result = payload.get("result", {})
        for uid in result.get("uids", []):
            rec = result.get(uid)
            if isinstance(rec, dict):
                out.append(rec)
    return out


def geo_ftp_listing(gse: str) -> list[dict[str, Any]]:
    """Scrape the GEO FTP HTTPS listing for a GSE and return supplementary files."""
    if not gse.upper().startswith("GSE"):
        return []
    gse_id = gse.upper()
    n = gse_id.replace("GSE", "")
    if not n.isdigit():
        return []
    # GEO uses GSE9nnn buckets for GSE9574 etc.
    bucket = f"GSE{n[:-3]}nnn" if len(n) >= 4 else "GSEnnn"
    base = f"{GEO_FTP_BASE}/series/{bucket}/{gse_id}/suppl/"
    # FTP listing via HTTPS mirror is an HTML index; parse hrefs of files.
    status, body = http_get(base, headers={"Accept": "text/html,*/*"})
    if status == 0 or status >= 400 or not body:
        return []
    text = body.decode("utf-8", errors="replace")
    rows: list[dict[str, Any]] = []
    for m in re.finditer(r'href=\"([^\"?]+?)\"', text):
        name = m.group(1)
        if name in ("../", "/", ""):
            continue
        if name.endswith("/"):
            continue
        url = base + name
        # try to grep the size from the same line; the HTML mirror prints a size next to the link
        rows.append({
            "filename": name,
            "url": url,
            "kind": classify_kind(None, None, name),
            "size_bytes": 0,  # FTP HTML listing's size cell isn't always reliable; fill from HEAD on download.
        })
    return rows


def survey_geo(*, limit: int = 50, label: str = "geo_multiome",
               extra_queries: list[str] | None = None,
               organism: str | None = None,
               fetch_files: bool = False) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Return (datasets, files) for GEO multiome studies."""
    queries = list(GEO_MULTIOME_QUERIES)
    if extra_queries:
        queries.extend(extra_queries)

    seen_uids: set[str] = set()
    seen_gse: set[str] = set()
    all_summaries: list[dict[str, Any]] = []

    for q in queries:
        uids = eutils_search_gds(q, retmax=limit, organism=organism)
        new_uids = [u for u in uids if u not in seen_uids]
        seen_uids.update(new_uids)
        logging.info("GEO esearch  q=%-50s -> %d uids (new: %d)", q[:50], len(uids), len(new_uids))
        if new_uids:
            all_summaries.extend(eutils_summary_gds(new_uids))
        # E-utilities asks for ≥3 requests/sec spacing without an API key; play nice.
        time.sleep(0.4)

    datasets: list[dict[str, Any]] = []
    files: list[dict[str, Any]] = []
    for rec in all_summaries:
        # gds returns ``GSE####`` in the ``accession`` field on series records.
        gse = rec.get("accession") or rec.get("entrytype", "")
        if not gse or not gse.upper().startswith("GSE"):
            # Some result types are GDS datasets, not series — skip.
            continue
        if gse in seen_gse:
            continue
        seen_gse.add(gse)
        samples = rec.get("samples") or []
        sample_summaries = [s.get("title") if isinstance(s, dict) else str(s) for s in samples[:3]]
        ds = {
            "source": "geo",
            "dataset_accession": gse,
            "title": rec.get("title") or "",
            "summary": rec.get("summary") or "",
            "organism": stringify(rec.get("taxon")),
            "assay": stringify(rec.get("gdstype") or rec.get("gpl")),
            "samples": stringify(sample_summaries),
            "n_samples": int(rec.get("n_samples") or len(samples) or 0),
            "publications": stringify(rec.get("pubmedids")),
            "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}",
            "lab_or_pi": stringify(rec.get("contact_name") or rec.get("contact")),
            "creation_date": rec.get("pdat") or rec.get("gdate", ""),
        }
        datasets.append(ds)

        if fetch_files:
            for f in geo_ftp_listing(gse):
                files.append({
                    "source": "geo",
                    "dataset_accession": gse,
                    "file_accession": f["filename"],
                    "kind": f["kind"],
                    "content_type": "",
                    "file_format": Path(f["filename"]).suffix.lstrip("."),
                    "file_size_bytes": f.get("size_bytes", 0),
                    "download_url": f["url"],
                    "local_path": rel_path(DOWNLOAD_DIR / "geo" / gse / f["filename"]),
                    "md5sum": "",
                    "notes": "GEO supplementary",
                })

    logging.info("GEO survey: %d unique GSEs, %d supplementary files (fetch_files=%s)",
                 len(datasets), len(files), fetch_files)
    return datasets, files


# ---------------------------------------------------------------------------
# Per-source writers
# ---------------------------------------------------------------------------


def write_survey_outputs(source: str, label: str,
                        datasets: list[dict[str, Any]],
                        files: list[dict[str, Any]]) -> dict[str, Path]:
    mkdirs()
    ts = timestamp()
    ds_csv = MANIFEST_DIR / f"{ts}_{safe_label(label)}_{source}_datasets.csv"
    fl_csv = MANIFEST_DIR / f"{ts}_{safe_label(label)}_{source}_files.csv"
    write_csv(ds_csv, datasets, DATASET_FIELDS)
    write_csv(fl_csv, files, FILE_FIELDS)
    report = REPORT_DIR / f"{ts}_{safe_label(label)}_{source}_report.md"
    write_source_report(report, source, datasets, files, ds_csv, fl_csv)
    logging.info("[%s] wrote datasets=%s files=%s report=%s", source,
                 rel_path(ds_csv), rel_path(fl_csv), rel_path(report))
    return {"datasets_csv": ds_csv, "files_csv": fl_csv, "report": report}


def write_source_report(path: Path, source: str,
                       datasets: list[dict[str, Any]],
                       files: list[dict[str, Any]],
                       ds_csv: Path, fl_csv: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append(f"# {source.upper()} multiome survey")
    lines.append("")
    lines.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    lines.append("")
    lines.append(f"- Datasets: **{len(datasets)}**")
    lines.append(f"- Files indexed: **{len(files)}**")
    lines.append(f"- Datasets CSV: `{rel_path(ds_csv)}`")
    lines.append(f"- Files CSV: `{rel_path(fl_csv)}`")
    lines.append("")

    if files:
        # Group sizes by kind
        from collections import Counter, defaultdict
        n_by_kind = Counter(f.get("kind", "other") for f in files)
        sz_by_kind: dict[str, int] = defaultdict(int)
        for f in files:
            sz_by_kind[f.get("kind", "other")] += int(f.get("file_size_bytes") or 0)
        lines.append("## File inventory by kind")
        lines.append("")
        lines.append("| kind | count | total GB |")
        lines.append("|---|---:|---:|")
        for k, n in n_by_kind.most_common():
            gb = sz_by_kind[k] / 1024 / 1024 / 1024
            lines.append(f"| {k} | {n} | {gb:.2f} |")
        lines.append("")

    lines.append("## First 25 datasets")
    lines.append("")
    lines.append("| accession | title | organism | assay | samples |")
    lines.append("|---|---|---|---|---|")
    for d in datasets[:25]:
        title = (d.get("title") or "")[:60].replace("|", "\\|")
        org = (d.get("organism") or "")[:20]
        assay = (d.get("assay") or "")[:30].replace("|", "\\|")
        samples = (d.get("samples") or "")[:40].replace("|", "\\|")
        lines.append(f"| [{d['dataset_accession']}]({d.get('url','')}) | {title} | {org} | {assay} | {samples} |")
    lines.append("")
    if len(datasets) > 25:
        lines.append(f"_… plus {len(datasets) - 25} more datasets in the CSV._")
        lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Unified manifest
# ---------------------------------------------------------------------------


def newest_files_csv(source: str) -> Path | None:
    cands = sorted(MANIFEST_DIR.glob(f"*_{source}_files.csv"))
    return cands[-1] if cands else None


def build_unified_manifest(label: str = "unified") -> Path:
    rows: list[dict[str, Any]] = []
    for source in ("igvf", "encode", "geo"):
        csv_path = newest_files_csv(source)
        if csv_path is None:
            logging.info("no files manifest yet for source=%s; skipping in unified", source)
            continue
        with csv_path.open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                rows.append(r)
    out = MANIFEST_DIR / f"{timestamp()}_{safe_label(label)}_unified_manifest.csv"
    write_csv(out, rows, FILE_FIELDS)
    logging.info("Wrote unified manifest %s (%d rows)", rel_path(out), len(rows))
    return out


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------


def newest_unified() -> Path | None:
    cands = sorted(MANIFEST_DIR.glob("*_unified_manifest.csv"))
    return cands[-1] if cands else None


def download_file(url: str, target: Path) -> dict[str, Any]:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0:
        return {"status": "exists", "bytes": target.stat().st_size}
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"})
    started = time.time()
    try:
        with urllib.request.urlopen(req, timeout=600) as r, target.open("wb") as fh:
            total = 0
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
                total += len(chunk)
        logging.info("downloaded %s (%.2f MB in %.1fs)",
                     target.name, total / 1024 / 1024, time.time() - started)
        return {"status": "downloaded", "bytes": total}
    except urllib.error.HTTPError as exc:
        return {"status": f"http_error_{exc.code}", "bytes": 0}
    except urllib.error.URLError as exc:
        return {"status": f"network_error_{exc.reason}", "bytes": 0}


def download_from_manifest(manifest: Path | None, *,
                           only: set[str] | None = None,
                           pattern: re.Pattern[str] | None = None,
                           max_gb: float = 5.0,
                           dry_run: bool = False) -> Path:
    if manifest is None:
        manifest = newest_unified()
    if manifest is None:
        raise SystemExit("No unified manifest found. Run `manifest` first.")
    rows = list(csv.DictReader(manifest.open()))
    total_bytes = 0
    cap_bytes = int(max_gb * 1024 ** 3)
    log_rows: list[dict[str, Any]] = []
    for r in rows:
        kind = r.get("kind", "")
        if only and kind not in only:
            continue
        if pattern and not pattern.search(r.get("file_accession", "") + " " + r.get("download_url", "")):
            continue
        size = int(float(r.get("file_size_bytes") or 0))
        if size and total_bytes + size > cap_bytes:
            logging.info("hit --max-download-gb cap (%.2f); skipping %s/%s",
                         max_gb, r.get("dataset_accession"), r.get("file_accession"))
            continue
        url = r.get("download_url", "")
        target = abs_from_rel(r.get("local_path", ""))
        if not url or not target:
            continue
        if dry_run:
            log_rows.append({**r, "status": "dry_run"})
            continue
        result = download_file(url, target)
        total_bytes += int(result.get("bytes") or 0)
        log_rows.append({**r, "status": result["status"], "bytes_actual": result["bytes"]})
        print(f"  {r.get('source','?'):6s} {r.get('dataset_accession',''):14s} "
              f"{r.get('file_accession','')[:18]:18s} {kind:14s} "
              f"{result['status']:>14s} {result['bytes']:>13,}")
    log_path = manifest.with_name(manifest.stem + "_download_log.csv")
    fields = list(FILE_FIELDS) + ["status", "bytes_actual"]
    write_csv(log_path, log_rows, fields)
    logging.info("Wrote download log: %s", rel_path(log_path))
    return log_path


# ---------------------------------------------------------------------------
# Local inventory
# ---------------------------------------------------------------------------


def inventory_local() -> Path:
    SURVEY_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for source_dir in SURVEY_DIR.iterdir():
        if not source_dir.is_dir() or source_dir.name.startswith("."):
            continue
        for ds_dir in source_dir.iterdir():
            if not ds_dir.is_dir() or ds_dir.name.startswith("."):
                continue
            for f in ds_dir.rglob("*"):
                if not f.is_file() or f.name.startswith("._"):
                    continue
                rows.append({
                    "source": source_dir.name,
                    "dataset_accession": ds_dir.name,
                    "file_accession": f.name,
                    "kind": classify_kind(None, None, f.name),
                    "content_type": "",
                    "file_format": f.suffix.lstrip("."),
                    "file_size_bytes": f.stat().st_size,
                    "download_url": "",
                    "local_path": rel_path(f),
                    "md5sum": "",
                    "notes": "on-disk",
                })
    rows.sort(key=lambda r: (r["source"], r["dataset_accession"], r["file_accession"]))
    write_csv(INVENTORY_CSV, rows, FILE_FIELDS)
    logging.info("Wrote inventory %s (%d files, %.2f GB)",
                 rel_path(INVENTORY_CSV), len(rows),
                 sum(r["file_size_bytes"] for r in rows) / 1024 / 1024 / 1024)
    return INVENTORY_CSV


# ---------------------------------------------------------------------------
# Skill playbook
# ---------------------------------------------------------------------------

SKILL_DOC_TEMPLATE = """# Skill: Multiome cross-source survey

Surveys **IGVF Portal**, **ENCODE**, and **GEO** for single-cell multiome
datasets (10x Multiome / SHARE-seq / single-nucleus multiome), produces a
unified file manifest, and downloads the training-relevant files into
``Data/MultiomeSurvey/``.

## Subcommands

### `survey-igvf`

```bash
igvfagent multiome-survey survey-igvf --limit 200 --fetch-files
```

Queries ``preferred_assay_titles=10x multiome`` (and SHARE-seq) on the
IGVF Portal.  Writes:

- `Data/Manifests/MultiomeSurvey/<ts>_<label>_igvf_datasets.csv`
- `Data/Manifests/MultiomeSurvey/<ts>_<label>_igvf_files.csv`
- `Docs/MultiomeSurvey/<ts>_<label>_igvf_report.md`

### `survey-encode`

```bash
igvfagent multiome-survey survey-encode --limit 200 --fetch-files
```

Searches ENCODE Experiments and Series for the multiome / SHARE-seq keywords.

### `survey-geo`

```bash
igvfagent multiome-survey survey-geo --limit 50 --organism 'Homo sapiens' --fetch-files
```

Runs NCBI E-utilities `esearch+esummary` on the `gds` database for the
canonical multiome queries.  Optionally scrapes the GEO FTP listing for
each GSE to expose supplementary files.

### `survey-all`

```bash
igvfagent multiome-survey survey-all --limit 100 --fetch-files
```

Runs all three.

### `manifest`

```bash
igvfagent multiome-survey manifest --label unified_v1
```

Builds a unified file manifest from the **most recent** per-source files
CSVs.  Use this as the input to `download`.

### `download`

```bash
igvfagent multiome-survey download \\
    --manifest Data/Manifests/MultiomeSurvey/<ts>_<label>_unified_manifest.csv \\
    --only matrix_rna,matrix_atac,fragments,annotations \\
    --max-download-gb 20
```

Cap defaults to 5 GB.  ``--only`` accepts a comma-separated list of
kinds: ``matrix_rna``, ``matrix_atac``, ``fragments``, ``peaks``,
``annotations``, ``alignments``, ``index``, ``raw_reads``, ``other``.
``--pattern`` is a case-insensitive regex applied to the filename + URL.

### `inventory`

```bash
igvfagent multiome-survey inventory
```

Walks ``Data/MultiomeSurvey/`` and writes a fresh inventory CSV.

## Output layout

```
Data/
  MultiomeSurvey/
    igvf/<IGVFDS...>/<file>
    encode/<ENCSR...>/<file>
    geo/<GSE...>/<file>
    inventory.csv
  Manifests/MultiomeSurvey/
    <ts>_<label>_<source>_datasets.csv
    <ts>_<label>_<source>_files.csv
    <ts>_<label>_unified_manifest.csv
    <ts>_<label>_unified_manifest_download_log.csv
Docs/
  MultiomeSurvey/
    <ts>_<label>_<source>_report.md
  Skills/MULTIOME_SURVEY_SKILLS.md   (this file)
```

## Privacy

All endpoint URLs are resolved through `Scripts/_endpoints.py`; no URLs,
cookies, or credentials are written to source.  `Data/MultiomeSurvey/`
matches `Data/*` in the repo `.gitignore`, so downloaded payload never
accidentally lands in commits.
"""


def write_playbook() -> Path:
    SKILL_DOC.parent.mkdir(parents=True, exist_ok=True)
    SKILL_DOC.write_text(SKILL_DOC_TEMPLATE, encoding="utf-8")
    logging.info("Wrote skill doc: %s", rel_path(SKILL_DOC))
    return SKILL_DOC


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_only(value: str | None) -> set[str] | None:
    if not value:
        return None
    return {v.strip() for v in value.split(",") if v.strip()}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    s_igvf = sub.add_parser("survey-igvf", help="Survey IGVF Portal multiome datasets.")
    s_igvf.add_argument("--limit", type=int, default=200)
    s_igvf.add_argument("--label", default="multiome_survey")
    s_igvf.add_argument("--fetch-files", action="store_true")
    s_igvf.add_argument("--no-shareseq", action="store_true",
                        help="Exclude SHARE-seq from the IGVF query.")

    s_enc = sub.add_parser("survey-encode", help="Survey ENCODE multiome datasets.")
    s_enc.add_argument("--limit", type=int, default=200)
    s_enc.add_argument("--label", default="multiome_survey")
    s_enc.add_argument("--fetch-files", action="store_true")

    s_geo = sub.add_parser("survey-geo", help="Survey GEO multiome studies.")
    s_geo.add_argument("--limit", type=int, default=50)
    s_geo.add_argument("--label", default="multiome_survey")
    s_geo.add_argument("--organism", default=None)
    s_geo.add_argument("--extra-query", action="append", default=None,
                       help="Additional GEO query strings (repeatable).")
    s_geo.add_argument("--fetch-files", action="store_true")

    s_all = sub.add_parser("survey-all", help="Run survey-igvf + survey-encode + survey-geo.")
    s_all.add_argument("--limit", type=int, default=100)
    s_all.add_argument("--label", default="multiome_survey")
    s_all.add_argument("--organism", default=None)
    s_all.add_argument("--fetch-files", action="store_true")

    s_man = sub.add_parser("manifest", help="Build unified manifest from latest surveys.")
    s_man.add_argument("--label", default="unified")

    s_dl = sub.add_parser("download", help="Download files from a unified manifest.")
    s_dl.add_argument("--manifest", default=None, help="Path to manifest CSV (default: newest).")
    s_dl.add_argument("--only", default=None,
                      help="Comma-separated kinds (matrix_rna,matrix_atac,fragments,annotations,…)")
    s_dl.add_argument("--pattern", default=None,
                      help="Case-insensitive regex applied to filename + URL.")
    s_dl.add_argument("--max-download-gb", type=float, default=5.0)
    s_dl.add_argument("--dry-run", action="store_true")

    sub.add_parser("inventory", help="Scan Data/MultiomeSurvey and emit inventory CSV.")
    sub.add_parser("write-playbook", help="Emit Docs/Skills/MULTIOME_SURVEY_SKILLS.md.")

    args = parser.parse_args(argv)
    setup_logging()
    mkdirs()

    if args.command == "survey-igvf":
        datasets, files = survey_igvf(
            limit=args.limit, label=args.label,
            fetch_files=args.fetch_files,
            include_shareseq=not args.no_shareseq,
        )
        out = write_survey_outputs("igvf", args.label, datasets, files)
        print(f"IGVF: {len(datasets)} datasets, {len(files)} files")
        print(f"  datasets: {rel_path(out['datasets_csv'])}")
        print(f"  files:    {rel_path(out['files_csv'])}")
        print(f"  report:   {rel_path(out['report'])}")
        return 0

    if args.command == "survey-encode":
        datasets, files = survey_encode(
            limit=args.limit, label=args.label, fetch_files=args.fetch_files,
        )
        out = write_survey_outputs("encode", args.label, datasets, files)
        print(f"ENCODE: {len(datasets)} datasets, {len(files)} files")
        print(f"  datasets: {rel_path(out['datasets_csv'])}")
        print(f"  files:    {rel_path(out['files_csv'])}")
        print(f"  report:   {rel_path(out['report'])}")
        return 0

    if args.command == "survey-geo":
        datasets, files = survey_geo(
            limit=args.limit, label=args.label, organism=args.organism,
            extra_queries=args.extra_query, fetch_files=args.fetch_files,
        )
        out = write_survey_outputs("geo", args.label, datasets, files)
        print(f"GEO: {len(datasets)} datasets, {len(files)} files")
        print(f"  datasets: {rel_path(out['datasets_csv'])}")
        print(f"  files:    {rel_path(out['files_csv'])}")
        print(f"  report:   {rel_path(out['report'])}")
        return 0

    if args.command == "survey-all":
        ds_i, fl_i = survey_igvf(limit=args.limit, label=args.label, fetch_files=args.fetch_files)
        write_survey_outputs("igvf", args.label, ds_i, fl_i)
        ds_e, fl_e = survey_encode(limit=args.limit, label=args.label, fetch_files=args.fetch_files)
        write_survey_outputs("encode", args.label, ds_e, fl_e)
        ds_g, fl_g = survey_geo(limit=args.limit, label=args.label,
                                organism=args.organism, fetch_files=args.fetch_files)
        write_survey_outputs("geo", args.label, ds_g, fl_g)
        build_unified_manifest(label=args.label)
        print(f"survey-all: igvf={len(ds_i)} encode={len(ds_e)} geo={len(ds_g)} datasets")
        return 0

    if args.command == "manifest":
        out = build_unified_manifest(label=args.label)
        print(f"unified manifest: {rel_path(out)}")
        return 0

    if args.command == "download":
        mpath = Path(args.manifest).resolve() if args.manifest else None
        only = _parse_only(args.only)
        pattern = re.compile(args.pattern, re.I) if args.pattern else None
        log = download_from_manifest(mpath, only=only, pattern=pattern,
                                     max_gb=args.max_download_gb, dry_run=args.dry_run)
        print(f"download log: {rel_path(log)}")
        return 0

    if args.command == "inventory":
        out = inventory_local()
        print(f"inventory: {rel_path(out)}")
        return 0

    if args.command == "write-playbook":
        out = write_playbook()
        print(f"wrote {rel_path(out)}")
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
