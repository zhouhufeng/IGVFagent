#!/usr/bin/env python3
"""Retrieve and summarize IGVF 10x multiome datasets.

The pipeline targets IGVF Portal AnalysisSet records with
``preferred_assay_titles=10x multiome``. It preserves complete JSON metadata,
writes analysis/file/sample manifests, and optionally downloads manageable files
with a total-size cap.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import logging
import os
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "Multiome10x"
PLOT_DIR = REPORT_DIR / "Plots"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
MANIFEST_DIR = DATA_DIR / "Manifests" / "Multiome10x"
METADATA_DIR = DATA_DIR / "IGVF" / "10xMultiome" / "Metadata"
DOWNLOAD_DIR = DATA_DIR / "IGVF" / "10xMultiome" / "Downloads"

PORTAL_API_BASE = os.environ.get("IGVF_PORTAL_API_BASE", "https://api.data.igvf.org").rstrip("/")
PORTAL_PUBLIC_BASE = os.environ.get("IGVF_PORTAL_PUBLIC_BASE", "https://data.igvf.org").rstrip("/")

PROCESSED_CONTENT_TYPES = {
    "annotated sparse peak count matrix",
    "sparse gene count matrix",
    "filtered feature barcode matrix",
    "cell annotations",
    "fragments",
    "peaks",
}


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"multiome_10x_pipeline_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)
    return log_path


def safe_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in label)


def mkdirs() -> None:
    for path in (DATA_DIR, REPORT_DIR, PLOT_DIR, MANIFEST_DIR, METADATA_DIR, DOWNLOAD_DIR, SKILL_DOC_DIR):
        path.mkdir(parents=True, exist_ok=True)


def build_url(base: str, path: str, params: dict[str, Any] | None = None) -> str:
    normalized = path if path.startswith("/") else f"/{path}"
    url = f"{base}{normalized}"
    if params:
        url = f"{url}?{urllib.parse.urlencode(params, doseq=True, quote_via=urllib.parse.quote)}"
    return url


def request_headers() -> dict[str, str]:
    headers = {
        "Accept": "application/json,*/*",
        "User-Agent": "IGVFdataAgent/0.1 10x-multiome",
    }
    if os.environ.get("IGVF_PORTAL_COOKIE"):
        headers["Cookie"] = os.environ["IGVF_PORTAL_COOKIE"]
    return headers


def fetch_json(url: str) -> tuple[int, Any]:
    logging.info("Request: GET %s", url)
    request = urllib.request.Request(url, headers=request_headers(), method="GET")
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            content = response.read()
            if "json" in response.headers.get("Content-Type", ""):
                return response.status, json.loads(content)
            return response.status, {"text_response": content.decode(errors="replace"), "url": url}
    except urllib.error.HTTPError as exc:
        content = exc.read()
        try:
            data = json.loads(content)
        except json.JSONDecodeError:
            data = {"http_error_body": content.decode(errors="replace"), "url": url}
        return exc.code, data
    except urllib.error.URLError as exc:
        return 0, {"network_error": str(exc.reason), "url": url}


def rows_from_response(data: Any) -> list[dict[str, Any]]:
    if isinstance(data, list):
        return [item for item in data if isinstance(item, dict)]
    if not isinstance(data, dict):
        return []
    for key in ("@graph", "graph", "results", "result", "data", "items"):
        value = data.get(key)
        if isinstance(value, list):
            return [item for item in value if isinstance(item, dict)]
    return []


def scalar_strings(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (str, int, float, bool)):
        return [str(value)]
    if isinstance(value, list):
        out: list[str] = []
        for item in value:
            out.extend(scalar_strings(item))
        return out
    if isinstance(value, dict):
        for key in ("title", "term_name", "accession", "@id", "uuid", "href", "summary", "name"):
            if key in value:
                return scalar_strings(value[key])
    return []


def text_value(value: Any, limit: int = 6) -> str:
    return "; ".join(scalar_strings(value)[:limit])


def counter_for(rows: list[dict[str, Any]], fields: tuple[str, ...], limit: int = 12) -> list[tuple[str, int]]:
    counter: collections.Counter[str] = collections.Counter()
    for row in rows:
        for field in fields:
            for value in scalar_strings(row.get(field)):
                if value:
                    counter[value] += 1
    return counter.most_common(limit)


def bytes_to_gb(size: int | float | None) -> float:
    return float(size or 0) / 1_000_000_000


def download_url(href: str | None) -> str:
    if not href:
        return ""
    if href.startswith("http://") or href.startswith("https://"):
        return href
    return f"{PORTAL_API_BASE}{href}"


def search_analysis_sets(limit: int) -> tuple[int, Any, list[dict[str, Any]], Path]:
    params = {
        "type": "AnalysisSet",
        "preferred_assay_titles": "10x multiome",
        "file_set_type": "principal analysis",
        "status": "released",
        "controlled_access": "false",
        "format": "json",
        "limit": str(limit),
    }
    url = build_url(PORTAL_API_BASE, "/search/", params)
    status, data = fetch_json(url)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    path = METADATA_DIR / f"{stamp}_10x_multiome_analysis_set_search.json"
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    rows = rows_from_response(data)
    logging.info("Search returned %s rows; saved %s", len(rows), path)
    return status, data, rows, path


def select_analysis_sets(rows: list[dict[str, Any]], count: int) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in rows:
        accession = row.get("accession")
        assays = set(scalar_strings(row.get("assay_titles")))
        searchable_text = " ".join(
            scalar_strings(row.get("summary"))
            + scalar_strings(row.get("description"))
            + scalar_strings(row.get("preferred_assay_titles"))
            + list(assays)
        ).lower()
        has_rna = any("rna" in assay.lower() for assay in assays) or "rna" in searchable_text
        has_atac = any("atac" in assay.lower() for assay in assays) or "atac" in searchable_text
        if not accession or accession in seen:
            continue
        if row.get("status") != "released" or row.get("controlled_access") is True:
            continue
        if row.get("file_set_type") != "principal analysis":
            continue
        if not (has_rna and has_atac):
            continue
        selected.append(row)
        seen.add(accession)
        if len(selected) >= count:
            break
    return selected


def fetch_item_detail(item_id: str) -> tuple[int, Any]:
    return fetch_json(build_url(PORTAL_API_BASE, item_id, {"format": "json"}))


def fetch_full_metadata(selected_rows: list[dict[str, Any]], *, fetch_files: bool) -> list[dict[str, Any]]:
    details: list[dict[str, Any]] = []
    for row in selected_rows:
        item_id = row.get("@id")
        if not item_id:
            continue
        status, detail = fetch_item_detail(item_id)
        if status < 200 or status >= 400 or not isinstance(detail, dict):
            logging.warning("Could not fetch detail for %s: HTTP %s", item_id, status)
            continue
        if fetch_files:
            hydrated_files = []
            for file_row in detail.get("files", []):
                file_id = file_row.get("@id") if isinstance(file_row, dict) else None
                if not file_id:
                    hydrated_files.append(file_row)
                    continue
                file_status, file_detail = fetch_item_detail(file_id)
                if 200 <= file_status < 400 and isinstance(file_detail, dict):
                    hydrated_files.append(file_detail)
                else:
                    hydrated_files.append(file_row)
            detail["files"] = hydrated_files
        details.append(detail)
        logging.info("Fetched %s with %s files", detail.get("accession"), len(detail.get("files", [])))
    return details


def save_full_metadata(details: list[dict[str, Any]], label: str) -> Path:
    path = METADATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_full_metadata.json"
    path.write_text(json.dumps(details, indent=2, sort_keys=True), encoding="utf-8")
    logging.info("Wrote full metadata: %s", path)
    return path


def dataset_manifest_rows(details: list[dict[str, Any]]) -> list[dict[str, str]]:
    rows = []
    for detail in details:
        samples = detail.get("samples", [])
        donors = detail.get("donors", [])
        files = [file_row for file_row in detail.get("files", []) if isinstance(file_row, dict)]
        rows.append(
            {
                "analysis_set_accession": text_value(detail.get("accession")),
                "analysis_set_id": text_value(detail.get("@id")),
                "status": text_value(detail.get("status")),
                "controlled_access": text_value(detail.get("controlled_access")),
                "file_set_type": text_value(detail.get("file_set_type")),
                "assay_titles": text_value(detail.get("assay_titles")),
                "preferred_assay_titles": text_value(detail.get("preferred_assay_titles")),
                "sample_summary": text_value(detail.get("sample_summary")),
                "sample_terms": text_value([sample.get("sample_terms") for sample in samples if isinstance(sample, dict)]),
                "donors": text_value(donors),
                "lab": text_value(detail.get("lab")),
                "award": text_value(detail.get("award")),
                "file_count": str(len(files)),
                "downloadable_file_count": str(sum(1 for file_row in files if file_row.get("href"))),
                "total_file_size_gb": f"{sum(bytes_to_gb(file_row.get('file_size')) for file_row in files):.3f}",
                "summary": text_value(detail.get("summary")),
            }
        )
    return rows


def file_manifest_rows(details: list[dict[str, Any]]) -> list[dict[str, str]]:
    rows = []
    for detail in details:
        set_accession = text_value(detail.get("accession"))
        sample_summary = text_value(detail.get("sample_summary"))
        for file_row in detail.get("files", []):
            if not isinstance(file_row, dict):
                continue
            href = file_row.get("href")
            rows.append(
                {
                    "analysis_set_accession": set_accession,
                    "file_accession": text_value(file_row.get("accession")),
                    "file_id": text_value(file_row.get("@id")),
                    "status": text_value(file_row.get("status")),
                    "content_type": text_value(file_row.get("content_type") or file_row.get("content_types")),
                    "content_summary": text_value(file_row.get("content_summary")),
                    "summary": text_value(file_row.get("summary")),
                    "file_format": text_value(file_row.get("file_format")),
                    "assembly": text_value(file_row.get("assembly")),
                    "transcriptome_annotation": text_value(file_row.get("transcriptome_annotation")),
                    "file_size_bytes": str(file_row.get("file_size") or ""),
                    "file_size_gb": f"{bytes_to_gb(file_row.get('file_size')):.3f}",
                    "md5sum": text_value(file_row.get("md5sum")),
                    "href": text_value(href),
                    "download_url": download_url(href),
                    "s3_uri": text_value(file_row.get("s3_uri")),
                    "sample_summary": sample_summary,
                }
            )
    return rows


def sample_manifest_rows(details: list[dict[str, Any]]) -> list[dict[str, str]]:
    rows = []
    for detail in details:
        set_accession = text_value(detail.get("accession"))
        for sample in detail.get("samples", []):
            if not isinstance(sample, dict):
                continue
            rows.append(
                {
                    "analysis_set_accession": set_accession,
                    "sample_accession": text_value(sample.get("accession")),
                    "sample_id": text_value(sample.get("@id")),
                    "status": text_value(sample.get("status")),
                    "taxa": text_value(sample.get("taxa")),
                    "classifications": text_value(sample.get("classifications")),
                    "sample_terms": text_value(sample.get("sample_terms")),
                    "summary": text_value(sample.get("summary")),
                    "institutional_certificates": text_value(sample.get("institutional_certificates")),
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote CSV: %s", path)


def write_manifests(details: list[dict[str, Any]], label: str) -> dict[str, Path]:
    stamp = time.strftime("%Y%m%d_%H%M%S")
    paths = {
        "analysis_sets": MANIFEST_DIR / f"{stamp}_{safe_label(label)}_analysis_sets.csv",
        "files": MANIFEST_DIR / f"{stamp}_{safe_label(label)}_files.csv",
        "samples": MANIFEST_DIR / f"{stamp}_{safe_label(label)}_samples.csv",
    }
    write_csv(paths["analysis_sets"], dataset_manifest_rows(details))
    write_csv(paths["files"], file_manifest_rows(details))
    write_csv(paths["samples"], sample_manifest_rows(details))
    return paths


def svg_bar_plot(items: list[tuple[str, int | float]], title: str, path: Path, *, width: int = 980, height: int = 440) -> None:
    top = items[:14]
    margin_left = 260
    margin_right = 40
    margin_top = 54
    row_height = 25
    chart_width = width - margin_left - margin_right
    max_value = max((float(value) for _, value in top), default=1.0)
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfbf8"/>',
        f'<text x="24" y="32" font-family="Arial, sans-serif" font-size="20" font-weight="700" fill="#1f2933">{escape_xml(title)}</text>',
    ]
    for index, (name, value) in enumerate(top):
        y = margin_top + index * row_height
        bar_width = 0 if max_value == 0 else int(chart_width * float(value) / max_value)
        lines.extend(
            [
                f'<text x="24" y="{y + 16}" font-family="Arial, sans-serif" font-size="12" fill="#28323c">{escape_xml(str(name)[:42])}</text>',
                f'<rect x="{margin_left}" y="{y + 3}" width="{bar_width}" height="16" rx="3" fill="#287c71"/>',
                f'<text x="{margin_left + bar_width + 8}" y="{y + 16}" font-family="Arial, sans-serif" font-size="12" fill="#28323c">{value}</text>',
            ]
        )
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")
    logging.info("Wrote plot: %s", path)


def escape_xml(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


def counter_table(items: list[tuple[str, int | float]]) -> str:
    if not items:
        return "- none"
    return "\n".join(f"- {name}: {value}" for name, value in items)


def summarize(details: list[dict[str, Any]]) -> dict[str, Any]:
    files = [file_row for detail in details for file_row in detail.get("files", []) if isinstance(file_row, dict)]
    samples = [sample for detail in details for sample in detail.get("samples", []) if isinstance(sample, dict)]
    return {
        "dataset_count": len(details),
        "file_count": len(files),
        "sample_count": len(samples),
        "downloadable_file_count": sum(1 for file_row in files if file_row.get("href")),
        "total_file_size_gb": round(sum(bytes_to_gb(file_row.get("file_size")) for file_row in files), 3),
        "assays": counter_for(details, ("assay_titles", "preferred_assay_titles")),
        "sample_terms": counter_for(samples, ("sample_terms", "summary")),
        "labs": counter_for(details, ("lab",)),
        "file_formats": counter_for(files, ("file_format",)),
        "content_types": counter_for(files, ("content_type", "content_types")),
        "assemblies": counter_for(files, ("assembly",)),
        "statuses": counter_for(details, ("status",)),
        "controlled_access": counter_for(details, ("controlled_access",)),
    }


def write_plots(summary: dict[str, Any], label: str) -> list[Path]:
    stamp = time.strftime("%Y%m%d_%H%M%S")
    plots = [
        PLOT_DIR / f"{stamp}_{safe_label(label)}_content_types.svg",
        PLOT_DIR / f"{stamp}_{safe_label(label)}_sample_terms.svg",
        PLOT_DIR / f"{stamp}_{safe_label(label)}_file_formats.svg",
    ]
    svg_bar_plot(summary["content_types"], "10x Multiome File Content Types", plots[0])
    svg_bar_plot(summary["sample_terms"], "10x Multiome Sample Terms", plots[1])
    svg_bar_plot(summary["file_formats"], "10x Multiome File Formats", plots[2])
    return plots


def should_download(file_row: dict[str, Any], policy: str) -> bool:
    if policy == "none":
        return False
    if policy == "all":
        return bool(file_row.get("href"))
    content = str(file_row.get("content_type") or "").lower()
    return bool(file_row.get("href")) and content in PROCESSED_CONTENT_TYPES


def download_files(details: list[dict[str, Any]], policy: str, max_download_gb: float, label: str) -> list[dict[str, str]]:
    if policy == "none":
        return []
    destination_root = DOWNLOAD_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    destination_root.mkdir(parents=True, exist_ok=True)
    downloaded: list[dict[str, str]] = []
    spent = 0.0
    for detail in details:
        dataset_dir = destination_root / safe_label(text_value(detail.get("accession")))
        dataset_dir.mkdir(parents=True, exist_ok=True)
        for file_row in detail.get("files", []):
            if not isinstance(file_row, dict) or not should_download(file_row, policy):
                continue
            file_size_gb = bytes_to_gb(file_row.get("file_size"))
            if spent + file_size_gb > max_download_gb:
                downloaded.append(
                    {
                        "analysis_set_accession": text_value(detail.get("accession")),
                        "file_accession": text_value(file_row.get("accession")),
                        "status": "skipped_size_cap",
                        "file_size_gb": f"{file_size_gb:.3f}",
                        "path": "",
                    }
                )
                continue
            url = download_url(file_row.get("href"))
            suffix = Path(url.split("?")[0]).name or text_value(file_row.get("accession"))
            output_path = dataset_dir / safe_label(suffix)
            logging.info("Downloading %s to %s", url, output_path)
            request = urllib.request.Request(url, headers=request_headers(), method="GET")
            try:
                with urllib.request.urlopen(request, timeout=120) as response, output_path.open("wb") as handle:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        handle.write(chunk)
                spent += file_size_gb
                downloaded.append(
                    {
                        "analysis_set_accession": text_value(detail.get("accession")),
                        "file_accession": text_value(file_row.get("accession")),
                        "status": "downloaded",
                        "file_size_gb": f"{file_size_gb:.3f}",
                        "path": str(output_path),
                    }
                )
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
                downloaded.append(
                    {
                        "analysis_set_accession": text_value(detail.get("accession")),
                        "file_accession": text_value(file_row.get("accession")),
                        "status": f"download_error: {exc}",
                        "file_size_gb": f"{file_size_gb:.3f}",
                        "path": str(output_path),
                    }
                )
    if downloaded:
        write_csv(MANIFEST_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_download_manifest.csv", downloaded)
    return downloaded


def write_report(
    details: list[dict[str, Any]],
    search_status: int,
    search_total: int | None,
    search_path: Path,
    metadata_path: Path,
    manifest_paths: dict[str, Path],
    summary: dict[str, Any],
    plots: list[Path],
    downloads: list[dict[str, str]],
    label: str,
    download_policy: str,
    max_download_gb: float,
) -> Path:
    path = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_report.md"
    lines = [
        "# IGVF 10x Multiome Retrieval And Processing Report",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        "## Retrieval",
        "",
        f"- Portal API base: `{PORTAL_API_BASE}`",
        "- Query: `type=AnalysisSet&preferred_assay_titles=10x multiome&file_set_type=principal analysis&status=released&controlled_access=false`",
        f"- Search HTTP status: `{search_status}`",
        f"- Portal reported total matching principal public AnalysisSets: `{search_total}`",
        f"- Selected datasets: `{summary['dataset_count']}`",
        f"- Full metadata JSON: `{metadata_path}`",
        f"- Search response JSON: `{search_path}`",
        f"- Analysis set manifest: `{manifest_paths['analysis_sets']}`",
        f"- File manifest: `{manifest_paths['files']}`",
        f"- Sample manifest: `{manifest_paths['samples']}`",
        "",
        "## Processed Inventory",
        "",
        f"- Files listed: `{summary['file_count']}`",
        f"- Downloadable files with href: `{summary['downloadable_file_count']}`",
        f"- Total listed file size: `{summary['total_file_size_gb']}` GB",
        f"- Download policy used in this run: `{download_policy}` with max `{max_download_gb}` GB",
        "",
        "Content types:",
        counter_table(summary["content_types"]),
        "",
        "File formats:",
        counter_table(summary["file_formats"]),
        "",
        "Sample terms and summaries:",
        counter_table(summary["sample_terms"]),
        "",
        "Labs:",
        counter_table(summary["labs"]),
        "",
        "Assemblies:",
        counter_table(summary["assemblies"]),
        "",
        "## Selected Datasets",
        "",
        "| AnalysisSet | Sample | Files | Size GB | Lab |",
        "|---|---|---:|---:|---|",
    ]
    for detail in details:
        files = [file_row for file_row in detail.get("files", []) if isinstance(file_row, dict)]
        size_gb = sum(bytes_to_gb(file_row.get("file_size")) for file_row in files)
        lines.append(
            "| "
            + " | ".join(
                [
                    text_value(detail.get("accession")),
                    text_value(detail.get("sample_summary")).replace("|", "/"),
                    str(len(files)),
                    f"{size_gb:.3f}",
                    text_value(detail.get("lab")).replace("|", "/"),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Plots",
            "",
            *[f"- `{plot}`" for plot in plots],
            "",
            "## Download Results",
            "",
        ]
    )
    if downloads:
        downloaded = sum(1 for row in downloads if row.get("status") == "downloaded")
        skipped = sum(1 for row in downloads if row.get("status") == "skipped_size_cap")
        errors = len(downloads) - downloaded - skipped
        lines.extend(
            [
                f"- Downloaded files: `{downloaded}`",
                f"- Skipped by size cap: `{skipped}`",
                f"- Download errors: `{errors}`",
            ]
        )
    else:
        lines.extend(
            [
                "- No file payloads were downloaded in this run.",
                "- The file manifest contains full Portal download URLs and S3 URIs for each selected dataset file.",
            ]
        )
    lines.extend(
        [
            "",
            "## Analysis Notes",
            "",
            "- Each selected AnalysisSet is a released, public principal analysis containing paired single-nucleus ATAC and RNA 10x multiome outputs.",
            "- The recurrent processed outputs are ATAC peak matrices, RNA gene count matrices, cell annotations, and ATAC fragments.",
            "- Use the cell annotations as the join key when integrating RNA and ATAC modalities; use `assembly` and `transcriptome_annotation` fields to keep reference builds aligned.",
            "- For full downstream processing, load RNA sparse gene matrices into AnnData/Scanpy and ATAC peak matrices/fragments into ArchR, Signac, or SnapATAC-style workflows, then integrate by barcode/cell annotation metadata.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")
    logging.info("Wrote report: %s", path)
    return path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def process_local_manifests(file_manifest: Path, download_manifest: Path, label: str) -> dict[str, Path]:
    file_rows = read_csv(file_manifest)
    download_rows = read_csv(download_manifest)
    path_by_accession = {
        row["file_accession"]: Path(row["path"])
        for row in download_rows
        if row.get("file_accession") and row.get("path") and row.get("status") == "downloaded"
    }
    annotation_rows: list[dict[str, str]] = []
    archive_rows: list[dict[str, str]] = []
    cell_terms: collections.Counter[str] = collections.Counter()
    brain_regions: collections.Counter[str] = collections.Counter()
    total_cells = 0
    cells_in_both = 0
    per_dataset_overlap: list[dict[str, Any]] = []

    for row in file_rows:
        path = path_by_accession.get(row.get("file_accession", ""))
        if not path or not path.exists():
            continue
        if row.get("content_type") == "cell annotations":
            with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fields = reader.fieldnames or []
                row_count = 0
                both_count = 0
                local_terms: collections.Counter[str] = collections.Counter()
                for record in reader:
                    row_count += 1
                    term = (record.get("CL_term_name") or record.get("cell_description") or "").strip()
                    brain_region = (record.get("BrainRegion") or "").strip()
                    if term:
                        cell_terms[term] += 1
                        local_terms[term] += 1
                    if brain_region:
                        brain_regions[brain_region] += 1
                    in_gex = str(record.get("in_GEX", "")).lower() in {"true", "1", "yes"}
                    in_atac = str(record.get("in_ATAC", "")).lower() in {"true", "1", "yes"}
                    if in_gex and in_atac:
                        both_count += 1
                total_cells += row_count
                cells_in_both += both_count
                per_dataset_overlap.append(
                    {
                        "analysis_set_accession": row.get("analysis_set_accession", ""),
                        "cells": row_count,
                        "cells_in_both_modalities": both_count,
                        "pct_both": round(100 * both_count / row_count, 2) if row_count else 0,
                    }
                )
                annotation_rows.append(
                    {
                        "analysis_set_accession": row.get("analysis_set_accession", ""),
                        "file_accession": row.get("file_accession", ""),
                        "path": str(path),
                        "rows_cells": str(row_count),
                        "cells_in_both_modalities": str(both_count),
                        "pct_both": f"{round(100 * both_count / row_count, 2) if row_count else 0}",
                        "columns": ";".join(fields),
                        "top_cell_terms": json.dumps(local_terms.most_common(10)),
                    }
                )
        elif row.get("content_type") == "sparse gene count matrix":
            members: list[dict[str, Any]] = []
            try:
                with tarfile.open(path, "r:gz") as tar:
                    for member in tar.getmembers()[:30]:
                        members.append(
                            {
                                "name": member.name,
                                "size": member.size,
                                "type": "dir" if member.isdir() else "file",
                            }
                        )
            except (tarfile.TarError, OSError) as exc:
                members.append({"error": str(exc)})
            archive_rows.append(
                {
                    "analysis_set_accession": row.get("analysis_set_accession", ""),
                    "file_accession": row.get("file_accession", ""),
                    "path": str(path),
                    "members_json": json.dumps(members),
                }
            )

    stamp = time.strftime("%Y%m%d_%H%M%S")
    annotation_path = REPORT_DIR / f"{stamp}_{safe_label(label)}_cell_annotation_summary.csv"
    archive_path = REPORT_DIR / f"{stamp}_{safe_label(label)}_rna_archive_summary.csv"
    summary_path = REPORT_DIR / f"{stamp}_{safe_label(label)}_processing_summary.json"
    report_path = REPORT_DIR / f"{stamp}_{safe_label(label)}_processing_report.md"
    write_csv(annotation_path, annotation_rows)
    write_csv(archive_path, archive_rows)
    summary = {
        "total_cells": total_cells,
        "cells_in_both_modalities": cells_in_both,
        "pct_both_modalities": round(100 * cells_in_both / total_cells, 2) if total_cells else 0,
        "annotation_files": len(annotation_rows),
        "rna_archives": len(archive_rows),
        "top_cell_terms": cell_terms.most_common(25),
        "brain_regions": brain_regions.most_common(),
        "per_dataset_modality_overlap": per_dataset_overlap,
        "annotation_summary_csv": str(annotation_path),
        "rna_archive_summary_csv": str(archive_path),
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    lines = [
        "# IGVF 10x Multiome Local Processing Report",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"- File manifest: `{file_manifest}`",
        f"- Download manifest: `{download_manifest}`",
        f"- Annotation summary CSV: `{annotation_path}`",
        f"- RNA archive summary CSV: `{archive_path}`",
        f"- Processing summary JSON: `{summary_path}`",
        "",
        "## Cell Annotation Summary",
        "",
        f"- Annotation files parsed: `{len(annotation_rows)}`",
        f"- Total annotated cells: `{total_cells}`",
        f"- Cells present in both RNA and ATAC flags: `{cells_in_both}` ({summary['pct_both_modalities']}%)",
        "",
        "Top cell terms:",
        counter_table(cell_terms.most_common(15)),
        "",
        "Brain regions:",
        counter_table(brain_regions.most_common()),
        "",
        "## Per-Dataset Modality Overlap",
        "",
        "| AnalysisSet | Cells | Cells in both modalities | Percent |",
        "|---|---:|---:|---:|",
    ]
    for item in per_dataset_overlap:
        lines.append(
            f"| {item['analysis_set_accession']} | {item['cells']} | {item['cells_in_both_modalities']} | {item['pct_both']} |"
        )
    report_path.write_text("\n".join(lines), encoding="utf-8")
    logging.info("Wrote local processing report: %s", report_path)
    return {
        "annotation_summary": annotation_path,
        "rna_archive_summary": archive_path,
        "processing_summary": summary_path,
        "processing_report": report_path,
    }


def write_playbook() -> Path:
    path = SKILL_DOC_DIR / "10X_MULTIOME_SKILLS.md"
    lines = [
        "# Skill: IGVF 10x Multiome Retrieval And Analysis",
        "",
        "Use this skill when the agent needs to retrieve, summarize, or prepare paired 10x multiome RNA/ATAC datasets from the IGVF Portal.",
        "",
        "## Run",
        "",
        "```bash",
        "python3 Scripts/multiome_10x_pipeline.py retrieve --count 20 --fetch-file-details",
        "python3 Scripts/multiome_10x_pipeline.py process-local --file-manifest Data/Manifests/Multiome10x/<files.csv> --download-manifest Data/Manifests/Multiome10x/<download_manifest.csv>",
        "```",
        "",
        "The script searches public released IGVF Portal `AnalysisSet` records where `preferred_assay_titles=10x multiome`, preserves full JSON metadata, writes dataset/file/sample manifests, and creates a processing report.",
        "",
        "## Inputs To Prefer",
        "",
        "- RNA sparse gene count matrices.",
        "- ATAC annotated sparse peak count matrices.",
        "- ATAC fragments.",
        "- Cell annotations for barcode-level modality integration.",
        "",
        "## Processing Pattern",
        "",
        "1. Build file and sample manifests before downloading payloads.",
        "2. Verify each file is released, public, and aligned to the expected assembly and transcriptome annotation.",
        "3. Load RNA matrices into AnnData/Scanpy and ATAC matrices or fragments into ArchR, Signac, or SnapATAC-style workflows.",
        "4. Join modalities with cell annotations and barcodes.",
        "5. Summarize QC metrics, cell states, accessible peaks, gene expression, and peak-gene or variant-gene context for IGVF interpretation.",
        "",
        "## Download Policy",
        "",
        "Use `--download-policy processed --max-download-gb N` only after reviewing the file manifest. Raw and processed multiome payloads can be large.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    logging.info("Wrote playbook: %s", path)
    return path


def run_retrieve(args: argparse.Namespace) -> int:
    mkdirs()
    search_status, search_data, rows, search_path = search_analysis_sets(max(args.count * 3, args.limit))
    selected = select_analysis_sets(rows, args.count)
    if len(selected) < args.count and len(rows) < args.limit:
        logging.warning("Selected only %s datasets from %s search rows.", len(selected), len(rows))
    details = fetch_full_metadata(selected, fetch_files=args.fetch_file_details)
    metadata_path = save_full_metadata(details, args.label)
    manifest_paths = write_manifests(details, args.label)
    summary = summarize(details)
    plots = write_plots(summary, args.label)
    downloads = download_files(details, args.download_policy, args.max_download_gb, args.label)
    report_path = write_report(
        details,
        search_status,
        search_data.get("total") if isinstance(search_data, dict) else None,
        search_path,
        metadata_path,
        manifest_paths,
        summary,
        plots,
        downloads,
        args.label,
        args.download_policy,
        args.max_download_gb,
    )
    print(f"Selected datasets: {summary['dataset_count']}")
    print(f"Files listed: {summary['file_count']}")
    print(f"Downloadable files: {summary['downloadable_file_count']}")
    print(f"Total listed size GB: {summary['total_file_size_gb']}")
    print(f"Metadata: {metadata_path}")
    print(f"Analysis manifest: {manifest_paths['analysis_sets']}")
    print(f"File manifest: {manifest_paths['files']}")
    print(f"Sample manifest: {manifest_paths['samples']}")
    print(f"Report: {report_path}")
    return 0 if summary["dataset_count"] >= min(args.count, 10) and summary["file_count"] else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Retrieve and summarize IGVF 10x multiome datasets.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    retrieve = subparsers.add_parser("retrieve", help="Retrieve 10x multiome AnalysisSet metadata and manifests.")
    retrieve.add_argument("--count", type=int, default=20, help="Number of public principal AnalysisSets to select.")
    retrieve.add_argument("--limit", type=int, default=80, help="Portal search rows to request.")
    retrieve.add_argument("--label", default="igvf_10x_multiome", help="Output label.")
    retrieve.add_argument("--fetch-file-details", action="store_true", help="Fetch each file object for complete file metadata.")
    retrieve.add_argument(
        "--download-policy",
        choices=["none", "processed", "all"],
        default="none",
        help="Optional file payload download policy.",
    )
    retrieve.add_argument("--max-download-gb", type=float, default=2.0, help="Maximum total payload download size.")

    process = subparsers.add_parser("process-local", help="Process downloaded 10x multiome payload manifests.")
    process.add_argument("--file-manifest", required=True, type=Path)
    process.add_argument("--download-manifest", required=True, type=Path)
    process.add_argument("--label", default="igvf_10x_multiome_local")

    subparsers.add_parser("write-playbook", help="Write the reusable 10x multiome skill document.")

    args = parser.parse_args(argv)
    setup_logging()
    if args.command == "write-playbook":
        mkdirs()
        path = write_playbook()
        print(f"Wrote {path}")
        return 0
    if args.command == "retrieve":
        return run_retrieve(args)
    if args.command == "process-local":
        mkdirs()
        paths = process_local_manifests(args.file_manifest, args.download_manifest, args.label)
        for name, path in paths.items():
            print(f"{name}: {path}")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
