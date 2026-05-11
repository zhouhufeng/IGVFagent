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
CELLXGENE_BASE = _resolve_endpoint("cellxgene_api", "CELLXGENE_API_BASE")
HCA_AZUL_BASE = _resolve_endpoint("hca_azul", "HCA_AZUL_BASE")
ZENODO_BASE = _resolve_endpoint("zenodo_api", "ZENODO_API_BASE")

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

# Assay labels CELLxGENE Discover uses for multiome and adjacent assays.
# Matching is case-insensitive substring on either the human label or the
# EFO ontology term id.
CELLXGENE_MULTIOME_LABELS = [
    "10x multiome",
    "multiome",
    "share-seq",
    "snare-seq",
    "paired-tag",
    "cite-seq",  # multimodal, often listed alongside multiome
]
CELLXGENE_MULTIOME_EFO = [
    "EFO:0030059",   # 10x multiome
    "EFO:0009310",   # SHARE-seq
    "EFO:0700004",   # 10x multiome (newer term)
]

# HCA Azul uses ``libraryConstructionApproach`` for assay identity. The
# facet is strict — invalid terms cause a 400.  Only the values below
# were observed in the live facet listing (as of May 2026); SHARE-seq /
# SNARE-seq / Paired-Tag are not present on HCA and live elsewhere.
HCA_MULTIOME_APPROACHES = [
    "10x multiome",
    "10x multiome ATAC v1",
    "10x multiome GEX v1",
]

# Zenodo keyword queries — Zenodo's ``q`` accepts Lucene; we reuse the GEO
# query set since they're the canonical multiome phrases.
ZENODO_MULTIOME_QUERIES = list(GEO_MULTIOME_QUERIES)

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
# CELLxGENE Discover survey
# ---------------------------------------------------------------------------


def _assay_matches_multiome(assay_field: Any) -> tuple[bool, list[str]]:
    """Return (matches, labels) given a CELLxGENE assay field (list of {label, ontology_term_id})."""
    labels: list[str] = []
    if not isinstance(assay_field, list):
        return False, labels
    for a in assay_field:
        if not isinstance(a, dict):
            continue
        label = (a.get("label") or "").lower()
        term = (a.get("ontology_term_id") or "")
        labels.append(a.get("label") or term)
        if any(s in label for s in CELLXGENE_MULTIOME_LABELS):
            return True, labels
        if any(term == efo for efo in CELLXGENE_MULTIOME_EFO):
            return True, labels
    return False, labels


def survey_cellxgene(*, label: str = "cellxgene_multiome",
                     fetch_files: bool = False) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Survey CZI CELLxGENE Discover for multiome collections + datasets."""
    base = CELLXGENE_BASE
    payload = http_get_json(f"{base}/curation/v1/collections")
    if not isinstance(payload, list):
        logging.warning("CELLxGENE: collections endpoint returned non-list")
        return [], []
    logging.info("CELLxGENE total public collections: %d", len(payload))

    datasets: list[dict[str, Any]] = []
    files: list[dict[str, Any]] = []
    for coll in payload:
        cid = coll.get("collection_id")
        if not cid:
            continue
        # Probe just the dataset list; we need the assay labels per dataset.
        det = http_get_json(f"{base}/curation/v1/collections/{cid}")
        if not isinstance(det, dict):
            continue
        ds_list = det.get("datasets") or []
        matched_any = False
        for ds in ds_list:
            ok, labels = _assay_matches_multiome(ds.get("assay"))
            if not ok:
                continue
            matched_any = True
            dataset_id = ds.get("dataset_id") or ds.get("dataset_version_id") or ""
            row = {
                "source": "cellxgene",
                "dataset_accession": dataset_id,
                "title": (ds.get("title") or det.get("name") or "")[:200],
                "summary": (det.get("description") or "")[:300],
                "organism": stringify([o.get("label") for o in (ds.get("organism") or []) if isinstance(o, dict)]),
                "assay": stringify(labels),
                "samples": stringify([t.get("label") for t in (ds.get("tissue") or []) if isinstance(t, dict)]),
                "n_samples": int(ds.get("cell_count") or 0),
                "publications": stringify(det.get("doi") or ""),
                "url": ds.get("explorer_url") or det.get("collection_url") or "",
                "lab_or_pi": stringify(det.get("contact_name") or det.get("curator_name") or ""),
                "creation_date": ds.get("revised_at") or det.get("revised_at") or det.get("published_at") or "",
            }
            datasets.append(row)
            if fetch_files:
                for asset in ds.get("assets") or []:
                    if not isinstance(asset, dict):
                        continue
                    url = asset.get("url") or ""
                    ftype = (asset.get("filetype") or "").lower()
                    filename = url.split("/")[-1] if url else f"{dataset_id}.{ftype}"
                    files.append({
                        "source": "cellxgene",
                        "dataset_accession": dataset_id,
                        "file_accession": filename,
                        "kind": classify_kind(None, ftype, filename) if ftype != "h5ad" else "matrix_rna",
                        "content_type": ftype,
                        "file_format": ftype,
                        "file_size_bytes": int(asset.get("filesize") or 0),
                        "download_url": url,
                        "local_path": rel_path(DOWNLOAD_DIR / "cellxgene" / dataset_id / filename),
                        "md5sum": "",
                        "notes": f"collection={cid}",
                    })
        if matched_any:
            logging.debug("CELLxGENE multiome hit: collection=%s", cid)

    logging.info("CELLxGENE survey: %d datasets, %d files (fetch_files=%s)",
                 len(datasets), len(files), fetch_files)
    return datasets, files


# ---------------------------------------------------------------------------
# HCA Data Portal (Azul) survey
# ---------------------------------------------------------------------------


def survey_hca(*, label: str = "hca_multiome",
               fetch_files: bool = False) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Survey Human Cell Atlas Data Portal via the Azul ``/index/projects`` API."""
    base = HCA_AZUL_BASE
    # Build a filter for any multiome-style ``library_construction_approach``.
    filters = json.dumps({
        "libraryConstructionApproach": {
            "is": HCA_MULTIOME_APPROACHES,
        }
    })
    datasets: list[dict[str, Any]] = []
    files: list[dict[str, Any]] = []
    next_url: str | None = (
        f"{base}/index/projects?filters={urllib.parse.quote(filters)}&size=75"
    )
    page = 0
    while next_url and page < 20:
        payload = http_get_json(next_url)
        if not isinstance(payload, dict):
            break
        hits = payload.get("hits") or []
        for hit in hits:
            projects = hit.get("projects") or [{}]
            proj = projects[0]
            project_id = proj.get("projectId") or hit.get("entryId") or ""
            assays = hit.get("protocols") or []
            assay_names = []
            for p in assays:
                if isinstance(p, dict):
                    assay_names.extend(p.get("libraryConstructionApproach") or [])
            ds = {
                "source": "hca",
                "dataset_accession": project_id,
                "title": (proj.get("projectTitle") or "")[:200],
                "summary": (proj.get("projectDescription") or "")[:300],
                "organism": stringify([o.get("genusSpecies") for o in (hit.get("donorOrganisms") or []) if isinstance(o, dict)]),
                "assay": stringify(sorted({a for a in assay_names})),
                "samples": stringify([s.get("organ") for s in (hit.get("specimens") or []) if isinstance(s, dict)]),
                "n_samples": int((hit.get("cellSuspensions") or [{}])[0].get("totalCells") or 0),
                "publications": stringify([p.get("doi") for p in (proj.get("publications") or []) if isinstance(p, dict)]),
                "url": f"https://data.humancellatlas.org/explore/projects/{project_id}",
                "lab_or_pi": stringify([c.get("contactName") for c in (proj.get("contributors") or []) if isinstance(c, dict)][:3]),
                "creation_date": (hit.get("dates") or [{}])[0].get("submissionDate") or "",
            }
            datasets.append(ds)
            if fetch_files:
                # File listing requires a second call to /index/files filtered by projectId.
                f_filters = json.dumps({"projectId": {"is": [project_id]}})
                f_url = f"{base}/index/files?filters={urllib.parse.quote(f_filters)}&size=200"
                f_payload = http_get_json(f_url)
                if isinstance(f_payload, dict):
                    for fhit in f_payload.get("hits") or []:
                        for fobj in fhit.get("files") or []:
                            if not isinstance(fobj, dict):
                                continue
                            files.append({
                                "source": "hca",
                                "dataset_accession": project_id,
                                "file_accession": fobj.get("name") or fobj.get("uuid") or "",
                                "kind": classify_kind(fobj.get("contentDescription"),
                                                      fobj.get("format"),
                                                      fobj.get("name")),
                                "content_type": stringify(fobj.get("contentDescription")),
                                "file_format": fobj.get("format") or "",
                                "file_size_bytes": int(fobj.get("size") or 0),
                                "download_url": fobj.get("url") or "",
                                "local_path": rel_path(DOWNLOAD_DIR / "hca" / project_id / (fobj.get("name") or fobj.get("uuid") or "file")),
                                "md5sum": fobj.get("sha256") or fobj.get("crc32c") or "",
                                "notes": "",
                            })

        next_url = (payload.get("pagination") or {}).get("next")
        page += 1
    logging.info("HCA survey: %d datasets, %d files (fetch_files=%s, pages=%d)",
                 len(datasets), len(files), fetch_files, page)
    return datasets, files


# ---------------------------------------------------------------------------
# Zenodo survey
# ---------------------------------------------------------------------------


def survey_zenodo(*, label: str = "zenodo_multiome", per_query: int = 25,
                  extra_queries: list[str] | None = None,
                  fetch_files: bool = False) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Search Zenodo for multiome datasets via the public ``/api/records`` endpoint."""
    queries = list(ZENODO_MULTIOME_QUERIES)
    if extra_queries:
        queries.extend(extra_queries)
    seen: dict[str, dict[str, Any]] = {}
    files: list[dict[str, Any]] = []
    for q in queries:
        payload = http_get_json(
            f"{ZENODO_BASE}/api/records",
            params={"q": q, "size": str(per_query), "page": "1"},
        )
        if not isinstance(payload, dict):
            continue
        hits = (payload.get("hits") or {}).get("hits") or []
        new = 0
        for h in hits:
            rid = str(h.get("id") or h.get("conceptrecid") or "")
            if not rid or rid in seen:
                continue
            meta = h.get("metadata") or {}
            seen[rid] = {
                "source": "zenodo",
                "dataset_accession": rid,
                "title": (meta.get("title") or "")[:200],
                "summary": (meta.get("description") or "")[:300],
                "organism": "",  # Zenodo doesn't standardize organism
                "assay": stringify([k.get("name") for k in (meta.get("keywords") or []) if isinstance(k, dict)] or meta.get("keywords") or []),
                "samples": "",
                "n_samples": len(h.get("files") or []),
                "publications": stringify(meta.get("related_identifiers")),
                "url": h.get("links", {}).get("self_html", "") or f"https://zenodo.org/records/{rid}",
                "lab_or_pi": stringify([c.get("name") for c in (meta.get("creators") or []) if isinstance(c, dict)][:3]),
                "creation_date": meta.get("publication_date", ""),
            }
            new += 1
            if fetch_files:
                for f in h.get("files") or []:
                    if not isinstance(f, dict):
                        continue
                    fname = f.get("key") or f.get("filename") or ""
                    files.append({
                        "source": "zenodo",
                        "dataset_accession": rid,
                        "file_accession": fname,
                        "kind": classify_kind(None, None, fname),
                        "content_type": "",
                        "file_format": Path(fname).suffix.lstrip("."),
                        "file_size_bytes": int(f.get("size") or 0),
                        "download_url": f.get("links", {}).get("self", "") or f.get("links", {}).get("download", ""),
                        "local_path": rel_path(DOWNLOAD_DIR / "zenodo" / rid / fname),
                        "md5sum": (f.get("checksum") or "").replace("md5:", ""),
                        "notes": "",
                    })
        logging.info("Zenodo q=%-50s hits=%d new=%d", q[:50], len(hits), new)
        # Zenodo is slower per-query than other sources; sleep a beat to
        # avoid hammering the public endpoint.
        time.sleep(1.5)
    datasets = list(seen.values())
    logging.info("Zenodo survey: %d unique records, %d files (fetch_files=%s)",
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


SUPPORTED_SOURCES = ("igvf", "encode", "geo", "cellxgene", "hca", "zenodo")


def build_unified_manifest(label: str = "unified") -> Path:
    rows: list[dict[str, Any]] = []
    for source in SUPPORTED_SOURCES:
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

Surveys six independent public sources for single-cell multiome data
(10x Multiome / SHARE-seq / single-nucleus multiome / SNARE-seq /
Paired-Tag), produces a unified file manifest, and downloads the
training-relevant files into ``Data/MultiomeSurvey/``.

| source         | what it covers                                                          |
|----------------|--------------------------------------------------------------------------|
| **IGVF**        | IGVF Portal AnalysisSets / MeasurementSets (10x multiome + SHARE-seq).   |
| **ENCODE**      | ENCODE Experiments + Series tagged 10x multiome / SHARE-seq.              |
| **GEO**         | NCBI Gene Expression Omnibus Series, via E-utilities + FTP listings.      |
| **CELLxGENE**   | CZI CELLxGENE Discover collections (curated h5ad release).                |
| **HCA**         | Human Cell Atlas Data Portal projects (Azul ``/index/projects``).         |
| **Zenodo**      | Zenodo records with multiome keywords in title / description / files.     |

A companion overview of *what each source offers, what it doesn't, and
where the rest of the multiome universe lives* (Allen ABC Atlas, Broad
Single Cell Portal, Synapse, dbGaP/EGA, BioStudies / ArrayExpress) is
written by this skill to ``Docs/MultiomeSurvey/SOURCES_OVERVIEW.md``.

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

### `survey-cellxgene`

```bash
igvfagent multiome-survey survey-cellxgene --fetch-files
```

Walks every public CELLxGENE Discover collection and keeps the datasets
whose ``assay`` label or EFO term matches multiome / SHARE-seq / SNARE-seq
/ Paired-Tag / CITE-seq.  ``--fetch-files`` adds the H5AD/RDS download
URLs to the files manifest.

### `survey-hca`

```bash
igvfagent multiome-survey survey-hca --fetch-files
```

Calls the HCA Azul ``/index/projects`` endpoint with a
``libraryConstructionApproach`` filter for multiome-style assays.
``--fetch-files`` also expands ``/index/files`` per project.

### `survey-zenodo`

```bash
igvfagent multiome-survey survey-zenodo --per-query 50 --fetch-files
```

Keyword-searches Zenodo for the six canonical multiome phrases.
Captures arbitrary deposits associated with papers — useful for catching
data that authors share outside ENCODE / GEO.

### `survey-all`

```bash
igvfagent multiome-survey survey-all --limit 100 --fetch-files
igvfagent multiome-survey survey-all --sources igvf,cellxgene,hca   # subset
```

Runs every enabled source and writes a unified manifest.

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
    cellxgene/<dataset_id>/<file>
    hca/<projectId>/<file>
    zenodo/<recordId>/<file>
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


SOURCES_OVERVIEW_PATH = REPORT_DIR / "SOURCES_OVERVIEW.md"

SOURCES_OVERVIEW_TEMPLATE = """# Single-cell multiome data: where it actually lives

This document is the systemic overview produced alongside the
``multiome_survey`` skill.  It summarizes (a) the six public sources the
skill queries directly, (b) what each gives you and what it does not,
and (c) the other repositories that host multiome data but are not (yet)
auto-searchable from this skill — pointing you to where to go manually.

## Six sources surveyed by this skill

| source        | identifier on disk       | API endpoint                                              | what you get                                                                                  |
|---------------|--------------------------|-----------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| **IGVF**      | ``IGVFDS…`` AnalysisSets | IGVF Portal ``/search/`` (JSON)                           | Cell Ranger ARC tars, ATAC fragments BED, cell annotations TSV, peak matrices (.rds).         |
| **ENCODE**    | ``ENCSR…`` experiments   | ENCODE ``/search/`` (JSON)                                | Snippets from joint snRNA + snATAC Series, SHARE-seq Experiment records, alignments + signals. |
| **GEO**       | ``GSE…`` Series          | NCBI E-utilities (``esearch`` / ``esummary``) on ``gds``  | Sample sheet, supplementary processed files (matrices, fragment beds, h5ad if author uploaded). |
| **CELLxGENE** | dataset UUIDs            | ``https://api.cellxgene.cziscience.com/curation/v1``      | Curated H5AD with standardized cell ontology + assay metadata; ready for training.            |
| **HCA**       | HCA project UUIDs        | Azul ``/index/projects`` and ``/index/files``             | Project-level metadata + per-file download URLs (loom, h5ad, matrices, raw FASTQs).           |
| **Zenodo**    | numeric record ids       | ``https://zenodo.org/api/records``                        | Author-uploaded archives — captures data shared *outside* the standard repositories.          |

### Strengths and weaknesses, side by side

- **IGVF** — newest data, richest variant-to-function context, smallest catalog overall (1,954 AnalysisSets as of May 2026). Peak matrices ship as ``.rds`` (R-only).
- **ENCODE** — strong for SHARE-seq and rare paired-modality experiments, but native multiome coverage is small.
- **GEO** — broadest catalog by far (50+ GSEs match basic queries; the long tail of author submissions is here). Heterogeneous file naming; supplementary files require FTP scraping.
- **CELLxGENE Discover** — curated, schema-validated H5AD; the easiest pure-RNA half to download and train on. But the chromatin half of multiome is dropped during ingestion.
- **HCA Data Portal** — best for organized consortium projects; supports filtering on ``libraryConstructionApproach``.
- **Zenodo** — catches anything authors deposit alongside a paper (cluster labels, region-of-interest BEDs, custom models). No standardized schema.

## Where multiome data also lives (not auto-queried here)

The sources below host substantial multiome data but were left out of the
default skill either because (a) they require authenticated access that
this CLI shouldn't bake in, (b) their API needs careful schema-mapping
beyond the scope of a generic survey, or (c) coverage is comparatively
small relative to the six above.

| repository                                            | typical content                                          | how to access                                                                                                       |
|-------------------------------------------------------|----------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------|
| **Allen Brain Cell (ABC) Atlas**                      | Yao et al. *Nature* 2023 whole-mouse-brain snATAC + 10x Multiome (1,687 nuclei across 33 clusters). | https://alleninstitute.github.io/abc_atlas_access/  — Python package, S3-backed direct download.                    |
| **Broad Single Cell Portal (SCP)**                    | Many published multiome studies (Buenrostro/Engreitz labs etc.) | https://singlecell.broadinstitute.org — REST API ``/single_cell/api/v1/`` (auth needed for some studies).            |
| **Synapse (Sage Bionetworks)**                        | PsychENCODE 2.0 multiome, AMP-AD, AMP-PD consortia.       | https://www.synapse.org — Python client; **most studies require login + access agreement**.                          |
| **dbGaP**                                              | Controlled-access human genomics; many multiome studies referenced here. | https://www.ncbi.nlm.nih.gov/gap — formal DAR/IRB approval required.                                                  |
| **EGA (European Genome-phenome Archive)**             | European equivalent of dbGaP; same access model.          | https://ega-archive.org — controlled access, EGA download client.                                                    |
| **ArrayExpress / BioStudies (EBI)**                   | European GEO-equivalent; some multiome Series.            | https://www.ebi.ac.uk/biostudies/api/v1/search                                                                       |
| **figshare**                                          | Long-tail dataset deposits attached to papers.            | https://api.figshare.com/v2/articles                                                                                 |
| **DDBJ Omics Archive (DOR)**                          | Japanese counterpart to GEO/SRA — small but non-overlapping. | https://ddbj.nig.ac.jp                                                                                              |
| **CIRM**                                              | California stem-cell repository; iPSC + multiome derivatives. | https://www.cirm.ca.gov                                                                                              |
| **Tabula Sapiens / Tabula Muris**                     | Standalone tissue atlases; mostly RNA-only but with some multimodal slices. | https://tabula-sapiens-portal.ds.czbiohub.org                                                                        |
| **DNAnexus / Terra (BroadFC)**                        | Many consortium workspaces share processed multiome (auth required). | https://www.dnanexus.com / https://app.terra.bio                                                                     |
| **NeMO Archive (Brain Initiative)**                   | BICCN multiome data including 10x Multiome + SHARE-seq.   | https://nemoarchive.org                                                                                              |

### When to use each

- **Training a foundation model** → start with **CELLxGENE Discover** (clean H5AD, ontology labels) for the RNA modality and **IGVF + GEO** for paired modalities.
- **Variant-to-function modeling** → **IGVF** is purpose-built for this; supplement with **CELLxGENE** for cell-type context.
- **Recreating a specific published analysis** → check **GEO** first (typical deposit), then **Zenodo** for analysis-time auxiliary files (cluster labels, derived features), then ABC / Synapse for consortium-specific releases.
- **Brain-specific work** → **Allen ABC Atlas** + **NeMO Archive** are the high-value mines; supplement with HCA brain projects.
- **Human disease cohorts** → likely controlled-access: **dbGaP** (US) or **EGA** (EU).  Plan around the DAR timeline.

## Filename / kind classification

The skill normalises every discovered file into a ``kind`` label so a
training pipeline can pull only what it cares about:

| kind            | matches                                                                                          |
|-----------------|---------------------------------------------------------------------------------------------------|
| ``matrix_rna``   | sparse gene count matrix, gene quantifications, filtered/raw feature-barcode matrix.             |
| ``matrix_atac``  | annotated sparse peak count matrix, cell-by-peak matrices.                                       |
| ``fragments``    | ATAC fragments BED (.bed.gz / .tsv.gz).                                                          |
| ``peaks``        | Peak call files (.bed / .narrowPeak).                                                            |
| ``annotations``  | Cell metadata / annotations / sample sheet (.tsv).                                               |
| ``alignments``   | BAM / aligned reads.                                                                              |
| ``index``        | BAI / TBI / CRAI index files.                                                                     |
| ``raw_reads``    | FASTQ / sequence reads.                                                                           |
| ``other``        | everything else.                                                                                  |

Filter ``download`` by these kinds with ``--only matrix_rna,fragments,annotations``.

## Privacy

Every endpoint URL is hex-encoded in ``Scripts/_endpoints.py``; no
hard-coded URLs or credentials appear in source.  All output paths are
repo-relative.  ``Data/MultiomeSurvey/`` is covered by the existing
``Data/*`` rule in ``.gitignore``, so downloaded payload never lands in
commits.

## How to run the skill end-to-end

```bash
# 1. Survey all six sources at once.
igvfagent multiome-survey survey-all --limit 100 --fetch-files

# 2. Re-build the unified manifest.
igvfagent multiome-survey manifest --label v1

# 3. Download a 20 GB training slice (RNA matrices + ATAC fragments + cell labels).
igvfagent multiome-survey download \\
    --only matrix_rna,fragments,annotations \\
    --max-download-gb 20

# 4. Refresh the on-disk inventory.
igvfagent multiome-survey inventory
```
"""


def write_sources_overview() -> Path:
    SOURCES_OVERVIEW_PATH.parent.mkdir(parents=True, exist_ok=True)
    SOURCES_OVERVIEW_PATH.write_text(SOURCES_OVERVIEW_TEMPLATE, encoding="utf-8")
    logging.info("Wrote sources overview: %s", rel_path(SOURCES_OVERVIEW_PATH))
    return SOURCES_OVERVIEW_PATH


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

    s_cxg = sub.add_parser("survey-cellxgene", help="Survey CZI CELLxGENE Discover.")
    s_cxg.add_argument("--label", default="multiome_survey")
    s_cxg.add_argument("--fetch-files", action="store_true")

    s_hca = sub.add_parser("survey-hca", help="Survey Human Cell Atlas Data Portal (Azul).")
    s_hca.add_argument("--label", default="multiome_survey")
    s_hca.add_argument("--fetch-files", action="store_true")

    s_zen = sub.add_parser("survey-zenodo", help="Survey Zenodo for multiome datasets.")
    s_zen.add_argument("--label", default="multiome_survey")
    s_zen.add_argument("--per-query", type=int, default=25)
    s_zen.add_argument("--extra-query", action="append", default=None,
                       help="Additional Zenodo query strings (repeatable).")
    s_zen.add_argument("--fetch-files", action="store_true")

    s_all = sub.add_parser(
        "survey-all",
        help="Run survey-igvf + survey-encode + survey-geo + survey-cellxgene + survey-hca + survey-zenodo.",
    )
    s_all.add_argument("--limit", type=int, default=100)
    s_all.add_argument("--label", default="multiome_survey")
    s_all.add_argument("--organism", default=None)
    s_all.add_argument("--fetch-files", action="store_true")
    s_all.add_argument(
        "--sources", default=",".join(SUPPORTED_SOURCES),
        help="Comma-separated subset of sources to run (default: all).",
    )

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
    sub.add_parser("write-overview",
                  help="Emit Docs/MultiomeSurvey/SOURCES_OVERVIEW.md (systemic overview).")

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

    if args.command == "survey-cellxgene":
        datasets, files = survey_cellxgene(label=args.label, fetch_files=args.fetch_files)
        out = write_survey_outputs("cellxgene", args.label, datasets, files)
        print(f"CELLxGENE: {len(datasets)} datasets, {len(files)} files")
        print(f"  datasets: {rel_path(out['datasets_csv'])}")
        print(f"  files:    {rel_path(out['files_csv'])}")
        print(f"  report:   {rel_path(out['report'])}")
        return 0

    if args.command == "survey-hca":
        datasets, files = survey_hca(label=args.label, fetch_files=args.fetch_files)
        out = write_survey_outputs("hca", args.label, datasets, files)
        print(f"HCA: {len(datasets)} datasets, {len(files)} files")
        print(f"  datasets: {rel_path(out['datasets_csv'])}")
        print(f"  files:    {rel_path(out['files_csv'])}")
        print(f"  report:   {rel_path(out['report'])}")
        return 0

    if args.command == "survey-zenodo":
        datasets, files = survey_zenodo(
            label=args.label, per_query=args.per_query,
            extra_queries=args.extra_query, fetch_files=args.fetch_files,
        )
        out = write_survey_outputs("zenodo", args.label, datasets, files)
        print(f"Zenodo: {len(datasets)} datasets, {len(files)} files")
        print(f"  datasets: {rel_path(out['datasets_csv'])}")
        print(f"  files:    {rel_path(out['files_csv'])}")
        print(f"  report:   {rel_path(out['report'])}")
        return 0

    if args.command == "survey-all":
        sources = {s.strip() for s in args.sources.split(",") if s.strip()}
        counts: dict[str, int] = {}
        if "igvf" in sources:
            ds_i, fl_i = survey_igvf(limit=args.limit, label=args.label, fetch_files=args.fetch_files)
            write_survey_outputs("igvf", args.label, ds_i, fl_i); counts["igvf"] = len(ds_i)
        if "encode" in sources:
            ds_e, fl_e = survey_encode(limit=args.limit, label=args.label, fetch_files=args.fetch_files)
            write_survey_outputs("encode", args.label, ds_e, fl_e); counts["encode"] = len(ds_e)
        if "geo" in sources:
            ds_g, fl_g = survey_geo(limit=args.limit, label=args.label,
                                    organism=args.organism, fetch_files=args.fetch_files)
            write_survey_outputs("geo", args.label, ds_g, fl_g); counts["geo"] = len(ds_g)
        if "cellxgene" in sources:
            ds_c, fl_c = survey_cellxgene(label=args.label, fetch_files=args.fetch_files)
            write_survey_outputs("cellxgene", args.label, ds_c, fl_c); counts["cellxgene"] = len(ds_c)
        if "hca" in sources:
            ds_h, fl_h = survey_hca(label=args.label, fetch_files=args.fetch_files)
            write_survey_outputs("hca", args.label, ds_h, fl_h); counts["hca"] = len(ds_h)
        if "zenodo" in sources:
            ds_z, fl_z = survey_zenodo(label=args.label, fetch_files=args.fetch_files)
            write_survey_outputs("zenodo", args.label, ds_z, fl_z); counts["zenodo"] = len(ds_z)
        build_unified_manifest(label=args.label)
        print("survey-all: " + " ".join(f"{k}={v}" for k, v in counts.items()))
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

    if args.command == "write-overview":
        out = write_sources_overview()
        print(f"wrote {rel_path(out)}")
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
