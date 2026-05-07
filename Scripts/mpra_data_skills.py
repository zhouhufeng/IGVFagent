#!/usr/bin/env python3
"""MPRA metadata retrieval and lightweight analysis skills.

This script can:
- Pull MPRA/STARR/BlueSTARR metadata from IGVF Catalog, IGVF Portal, and ENCODE.
- Summarize local MPRA result tables.
- Generate dependency-free SVG plots for quick review.
"""

from __future__ import annotations

import argparse
import collections
import csv
import html
import json
import logging
import math
import os
import statistics
import sys
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
REPORT_DIR = DOCS_DIR / "MPRA"
PLOT_DIR = REPORT_DIR / "Plots"
MANIFEST_DIR = DATA_DIR / "Manifests" / "MPRA"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

CATALOG_API_BASE = os.environ.get(
    "IGVF_CATALOG_API_BASE", "https://api.catalogkg.igvf.org"
).rstrip("/")
PORTAL_BASE = os.environ.get("IGVF_PORTAL_BASE", "https://api.data.igvf.org").rstrip("/")
ENCODE_BASE = os.environ.get("ENCODE_BASE", "https://www.encodeproject.org").rstrip("/")


CATALOG_MPRA_QUERIES = [
    (
        "catalog_mpra_variant_biosample",
        "/api/variants/biosamples",
        {"method": "MPRA", "limit": "25", "page": "0", "verbose": "false"},
        "Catalog variant-biosample MPRA evidence.",
    ),
    (
        "catalog_starr_variant_biosample",
        "/api/variants/biosamples",
        {"method": "STARR-seq", "limit": "25", "page": "0", "verbose": "false"},
        "Catalog variant-biosample STARR-seq evidence.",
    ),
    (
        "catalog_bluestarr_variant_biosample",
        "/api/variants/biosamples",
        {"method": "BlueSTARR", "limit": "25", "page": "0", "verbose": "false"},
        "Catalog variant-biosample BlueSTARR evidence.",
    ),
]

PORTAL_MPRA_QUERIES = [
    ("portal_mpra_measurement_sets", {"type": "MeasurementSet", "assay_titles": "MPRA", "limit": "25"}, "IGVF MPRA measurement sets."),
    ("portal_starr_measurement_sets", {"type": "MeasurementSet", "assay_titles": "STARR-seq", "limit": "25"}, "IGVF STARR-seq measurement sets."),
    ("portal_mpra_files", {"type": "File", "searchTerm": "MPRA", "limit": "25"}, "IGVF Portal files matching MPRA."),
]

ENCODE_MPRA_QUERIES = [
    ("encode_mpra_experiments", {"type": "Experiment", "searchTerm": "MPRA", "limit": "25"}, "ENCODE experiments matching MPRA."),
    ("encode_starr_experiments", {"type": "Experiment", "searchTerm": "STARR", "limit": "25"}, "ENCODE experiments matching STARR."),
    ("encode_reporter_files", {"type": "File", "searchTerm": "reporter", "limit": "25"}, "ENCODE reporter-assay-related files."),
]

PORTAL_MANY_MPRA_QUERIES = [
    ("portal_mpra_measurement_sets", {"type": "MeasurementSet", "searchTerm": "MPRA", "limit": "100"}, "IGVF Portal measurement sets matching MPRA."),
    ("portal_starr_measurement_sets", {"type": "MeasurementSet", "searchTerm": "STARR", "limit": "100"}, "IGVF Portal measurement sets matching STARR/BlueSTARR."),
    ("portal_mpra_analysis_sets", {"type": "AnalysisSet", "searchTerm": "MPRA", "limit": "100"}, "IGVF Portal MPRA analysis sets."),
    ("portal_mpra_files", {"type": "File", "searchTerm": "MPRA", "limit": "100"}, "IGVF Portal files matching MPRA."),
    ("portal_starr_files", {"type": "File", "searchTerm": "STARR", "limit": "100"}, "IGVF Portal files matching STARR/BlueSTARR."),
    ("portal_reporter_files", {"type": "File", "searchTerm": "reporter", "limit": "100"}, "IGVF Portal reporter-assay files."),
]

KNOWN_IGVF_MPRA_DATASETS = [
    ("IGVFDS1589ODOW", "MPRAsnakeflow paper: 240K-HepG2 and 240K-HEK293T"),
    ("IGVFDS1419ZPHD", "MPRAsnakeflow paper: 80K-neurons and 80K-WTC11"),
    ("IGVFDS5307RCQG", "MPRAsnakeflow paper: 20K-HepG2_P"),
    ("IGVFDS8114HMGD", "MPRAsnakeflow paper: 12K-cardiop and 12K-cardiom"),
    ("IGVFDS9668QVZD", "MPRAsnakeflow paper: 8K-neurons"),
]

MPRA_LITERATURE_USE_CASES = [
    {
        "paper": "Tewhey et al.",
        "year": "2016",
        "venue": "Cell / PubMed",
        "url": "https://pubmed.ncbi.nlm.nih.gov/27259153/",
        "application": "causal regulatory variant prioritization",
        "plot_templates": "allelic effect volcano, significant variant table, GWAS/eQTL overlap summary",
        "lesson": "Use MPRA log2 allelic effects and significance to prioritize noncoding variants that directly modulate expression.",
    },
    {
        "paper": "Klein et al.",
        "year": "2020",
        "venue": "Nature Methods",
        "url": "https://www.nature.com/articles/s41592-020-0965-y",
        "application": "assay design and context benchmarking",
        "plot_templates": "design-by-activity bars, input/output count QC, activity dynamic-range comparison",
        "lesson": "Compare construct design, sequence length, and assay context before over-interpreting activity differences.",
    },
    {
        "paper": "Zhao et al.",
        "year": "2023",
        "venue": "Nature Genetics",
        "url": "https://www.nature.com/articles/s41588-022-01278-7",
        "application": "cell-type-specific regulatory activity",
        "plot_templates": "cell-type activity heatmap, cluster-colored activity plot, element-by-cell-type specificity summary",
        "lesson": "Single-cell MPRA extends reporter readouts to heterogeneous cell mixtures and highlights cell-type-specific enhancers.",
    },
    {
        "paper": "Melnikov et al.",
        "year": "2012",
        "venue": "Nature Biotechnology",
        "url": "https://www.nature.com/articles/nbt.2137",
        "application": "regulatory grammar and enhancer optimization",
        "plot_templates": "saturation-mutagenesis effect map, sequence-position activity heatmap, model-vs-observed scatter",
        "lesson": "Dense MPRA tiling/mutagenesis can map functional bases and train sequence-activity models.",
    },
    {
        "paper": "ENCODE-style lentiMPRA consortium study",
        "year": "2024",
        "venue": "Nature",
        "url": "https://www.nature.com/articles/s41586-024-08430-9",
        "application": "large-scale cCRE functional annotation",
        "plot_templates": "element activity distribution, cCRE-class enrichment, active-element genomic browser view",
        "lesson": "MPRA can convert candidate regulatory catalogs into functional activity maps across biological contexts.",
    },
    {
        "paper": "Oliveros et al.",
        "year": "2023",
        "venue": "Cell Genomics / Cell Press",
        "url": "https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00134-3",
        "application": "disease-locus regulatory variant characterization",
        "plot_templates": "variant effect volcano, locus summary, target-gene interpretation table",
        "lesson": "MPRA variant effects are most useful when integrated with GWAS fine-mapping, expression, chromatin, and gene-linked evidence.",
    },
    {
        "paper": "MPRAsnakeflow / IGVF MPRA uniform processing",
        "year": "2025",
        "venue": "PubMed Central / IGVF workflow paper",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12621732/",
        "application": "uniform processing of IGVF and ENCODE MPRA datasets",
        "plot_templates": "dataset inventory bars, barcode-count QC, replicate/activity summary panels",
        "lesson": "For many IGVF MPRA datasets, start with a manifest, then compare processed activity outputs with the same QC and visualization templates.",
    },
]


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"mpra_data_skills_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)
    return log_path


def safe_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in label)


def save_json(label: str, data: Any) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}.json"
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    logging.info("Saved JSON: %s", path)
    return path


def build_url(base: str, path: str, params: dict[str, Any]) -> str:
    return f"{base}{path}?{urllib.parse.urlencode(params, doseq=True)}"


def fetch_json(
    base: str,
    path: str,
    label: str,
    params: dict[str, Any],
    *,
    portal_auth: bool = False,
) -> tuple[int, Any, Path]:
    query = dict(params)
    if base in (PORTAL_BASE, ENCODE_BASE):
        query["format"] = "json"
    headers = {"Accept": "application/json,*/*", "User-Agent": "IGVFdataAgent/0.1"}
    if portal_auth and os.environ.get("IGVF_PORTAL_COOKIE"):
        headers["Cookie"] = os.environ["IGVF_PORTAL_COOKIE"]
    url = build_url(base, path, query)
    logging.info("Request: GET %s", url)
    request = urllib.request.Request(url, headers=headers, method="GET")
    try:
        with urllib.request.urlopen(request, timeout=45) as response:
            content = response.read()
            data = json.loads(content) if "json" in response.headers.get("Content-Type", "") else content.decode(errors="replace")
            saved = save_json(label, data)
            return response.status, data, saved
    except urllib.error.HTTPError as exc:
        content = exc.read()
        try:
            data = json.loads(content)
        except json.JSONDecodeError:
            data = {"http_error_body": content.decode(errors="replace")}
        saved = save_json(label, data)
        return exc.code, data, saved
    except urllib.error.URLError as exc:
        data = {"network_error": str(exc.reason), "url": url}
        saved = save_json(label, data)
        return 0, data, saved


def rows_from_response(data: Any) -> list[dict[str, Any]]:
    if isinstance(data, list):
        return [item for item in data if isinstance(item, dict)]
    if not isinstance(data, dict):
        return []
    for key in ("@graph", "graph", "results", "result", "data", "items"):
        value = data.get(key)
        if isinstance(value, list):
            return [item for item in value if isinstance(item, dict)]
    if data.get("@id") or data.get("accession"):
        return [data]
    return []


def total_from_response(data: Any, rows: list[dict[str, Any]]) -> int:
    if isinstance(data, dict):
        for key in ("total", "@total", "total_results", "count"):
            if isinstance(data.get(key), int):
                return int(data[key])
    return len(rows)


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
        for key in ("accession", "@id", "href", "method", "source", "name", "title", "biological_context"):
            if key in value:
                return scalar_strings(value[key])
    return []


def text_value(value: Any) -> str:
    return ";".join(scalar_strings(value)[:6])


def count_field(rows: list[dict[str, Any]], fields: tuple[str, ...], limit: int = 12) -> list[tuple[str, int]]:
    counter: collections.Counter[str] = collections.Counter()
    for row in rows:
        for field in fields:
            for value in scalar_strings(row.get(field)):
                if value:
                    counter[value] += 1
    return counter.most_common(limit)


def write_manifest(rows: list[dict[str, str]], label: str) -> Path:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    path = MANIFEST_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_manifest.csv"
    fields = ["label", "description", "accession", "id", "status", "method", "source", "context", "file_format", "output_type", "href", "source_url"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


def manifest_rows(label: str, description: str, data: Any, limit: int | None = 10) -> list[dict[str, str]]:
    rows = []
    source_rows = rows_from_response(data)
    if limit is not None:
        source_rows = source_rows[:limit]
    for row in source_rows:
        rows.append(
            {
                "label": label,
                "description": description,
                "accession": text_value(row.get("accession")),
                "id": text_value(row.get("@id") or row.get("uuid")),
                "status": text_value(row.get("status")),
                "method": text_value(row.get("method") or row.get("assay_title") or row.get("assay_titles")),
                "source": text_value(row.get("source")),
                "context": text_value(row.get("biological_context") or row.get("biosample_summary")),
                "file_format": text_value(row.get("file_format")),
                "output_type": text_value(row.get("output_type") or row.get("content_type")),
                "href": text_value(row.get("href")),
                "source_url": text_value(row.get("source_url")),
            }
        )
    return rows


def write_portal_manifest_report(summaries: list[dict[str, Any]], manifest_path: Path, plots: list[Path], label: str) -> Path:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    path = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_portal_many_report.md"
    lines = [
        "# IGVF Portal MPRA Dataset Pull",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"Portal base: `{PORTAL_BASE}`",
        f"Manifest: `{manifest_path}`",
        "",
        "## Query Results",
        "",
    ]
    for summary in summaries:
        lines.extend(
            [
                f"### {summary['label']}",
                "",
                summary["description"],
                "",
                f"- HTTP status: {summary['http_status']}",
                f"- Returned rows: {summary['returned_rows']}",
                f"- Reported total: {summary['total']}",
                f"- Saved response: `{summary['saved_response']}`",
                "",
                "Top file/output types:",
                format_counter(summary["file_formats"]),
                "",
                "Top contexts:",
                format_counter(summary["contexts"]),
                "",
            ]
        )
    lines.extend(["## Plots", "", *[f"- `{plot}`" for plot in plots], ""])
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def row_text(row: dict[str, Any]) -> str:
    try:
        return json.dumps(row, sort_keys=True).lower()
    except TypeError:
        return str(row).lower()


def filter_rows_by_terms(data: Any, terms: tuple[str, ...]) -> Any:
    rows = rows_from_response(data)
    matched = [row for row in rows if any(term.lower() in row_text(row) for term in terms)]
    if isinstance(data, dict):
        filtered = dict(data)
        for key in ("@graph", "graph", "results", "result", "data", "items"):
            if isinstance(filtered.get(key), list):
                filtered[key] = matched
                filtered["total"] = len(matched)
                return filtered
    return matched


def summarize_remote(label: str, description: str, status: int, data: Any, saved: Path) -> dict[str, Any]:
    rows = rows_from_response(data)
    return {
        "label": label,
        "description": description,
        "http_status": status,
        "saved_response": str(saved),
        "returned_rows": len(rows),
        "total": total_from_response(data, rows),
        "methods": count_field(rows, ("method", "assay_title", "assay_titles")),
        "sources": count_field(rows, ("source", "lab", "award")),
        "contexts": count_field(rows, ("biological_context", "biosample_summary")),
        "file_formats": count_field(rows, ("file_format", "output_type", "content_type")),
    }


def format_counter(items: list[tuple[str, int]]) -> str:
    return "\n".join(f"- {name}: {count}" for name, count in items) if items else "- none observed"


def run_pull(args: argparse.Namespace) -> int:
    summaries: list[dict[str, Any]] = []
    manifest: list[dict[str, str]] = []
    sources = ["catalog", "portal", "encode"] if args.source == "all" else [args.source]
    for source in sources:
        if source == "catalog":
            for label, path, params, description in CATALOG_MPRA_QUERIES:
                query = dict(params)
                query["limit"] = str(args.limit)
                status, data, saved = fetch_json(CATALOG_API_BASE, path, label, query)
                summaries.append(summarize_remote(label, description, status, data, saved))
                manifest.extend(manifest_rows(label, description, data))
                print(f"{label}: HTTP {status}, rows={len(rows_from_response(data))}")
        elif source == "portal":
            for label, params, description in PORTAL_MPRA_QUERIES:
                query = dict(params)
                query["limit"] = str(args.limit)
                status, data, saved = fetch_json(PORTAL_BASE, "/search/", label, query, portal_auth=True)
                summaries.append(summarize_remote(label, description, status, data, saved))
                manifest.extend(manifest_rows(label, description, data))
                print(f"{label}: HTTP {status}, rows={len(rows_from_response(data))}")
        elif source == "encode":
            for label, params, description in ENCODE_MPRA_QUERIES:
                query = dict(params)
                query["limit"] = str(args.limit)
                status, data, saved = fetch_json(ENCODE_BASE, "/search/", label, query)
                summaries.append(summarize_remote(label, description, status, data, saved))
                manifest.extend(manifest_rows(label, description, data))
                print(f"{label}: HTTP {status}, rows={len(rows_from_response(data))}")
    save_json(f"mpra_{args.source}_summary", summaries)
    manifest_path = write_manifest(manifest, f"mpra_{args.source}")
    report_path = write_remote_report(summaries, manifest_path, f"mpra_{args.source}")
    print(f"Wrote manifest: {manifest_path}")
    print(f"Wrote report: {report_path}")
    return 0 if any(200 <= item["http_status"] < 400 and item["returned_rows"] > 0 for item in summaries) else 1


def write_remote_report(summaries: list[dict[str, Any]], manifest_path: Path, label: str) -> Path:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    path = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_metadata_report.md"
    lines = ["# MPRA Metadata Retrieval Report", "", f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}", "", f"Manifest: `{manifest_path}`", ""]
    for summary in summaries:
        lines.extend(
            [
                f"## {summary['label']}",
                "",
                summary["description"],
                "",
                f"HTTP status: {summary['http_status']}",
                f"Returned rows: {summary['returned_rows']}",
                f"Reported total: {summary['total']}",
                f"Saved response: `{summary['saved_response']}`",
                "",
                "Top methods/assays:",
                format_counter(summary["methods"]),
                "",
                "Top sources:",
                format_counter(summary["sources"]),
                "",
                "Top contexts:",
                format_counter(summary["contexts"]),
                "",
                "Top file/output types:",
                format_counter(summary["file_formats"]),
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    value = value.strip()
    if value in {"", "NA", "NaN", "#N/A", "FALSE", "TRUE"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_local_mpra(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def numeric_column(rows: list[dict[str, str]], column: str) -> list[float]:
    return [value for value in (parse_float(row.get(column)) for row in rows) if value is not None]


def stats(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"n": 0, "mean": None, "median": None, "min": None, "max": None}
    return {
        "n": len(values),
        "mean": statistics.fmean(values),
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
    }


def svg_header(width: int, height: int) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n'
        '<rect width="100%" height="100%" fill="white"/>\n'
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#1f2933}.title{font-size:18px;font-weight:700}.axis{stroke:#344054;stroke-width:1.2}.tick{font-size:11px;fill:#4f5b66}.label{font-size:13px;font-weight:600}</style>\n'
    )


def write_histogram(values: list[float], path: Path, title: str, bins: int = 30, x_label: str = "value", y_label: str = "count") -> None:
    width, height = 760, 420
    left, right, top, bottom = 70, 30, 50, 60
    if not values:
        path.write_text(svg_header(width, height) + f'<text x="20" y="40">{title}: no data</text></svg>\n', encoding="utf-8")
        return
    lo, hi = min(values), max(values)
    if lo == hi:
        lo -= 0.5
        hi += 0.5
    counts = [0] * bins
    for value in values:
        idx = min(bins - 1, int((value - lo) / (hi - lo) * bins))
        counts[idx] += 1
    max_count = max(counts) or 1
    plot_w = width - left - right
    plot_h = height - top - bottom
    bar_w = plot_w / bins
    parts = [svg_header(width, height), f'<text class="title" x="{left}" y="30">{html.escape(title)}</text>']
    parts.append(f'<line class="axis" x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}"/>')
    parts.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}"/>')
    for tick in range(5):
        y = height - bottom - tick * plot_h / 4
        count = max_count * tick / 4
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}" stroke="#d0d5dd" stroke-width=".8"/>')
        parts.append(f'<text class="tick" x="{left - 8}" y="{y + 4:.1f}" text-anchor="end">{count:.0f}</text>')
    for i, count in enumerate(counts):
        bar_h = count / max_count * plot_h
        x = left + i * bar_w
        y = height - bottom - bar_h
        parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{max(bar_w-1, 1):.2f}" height="{bar_h:.2f}" fill="#4C78A8"/>')
    parts.append(f'<text class="tick" x="{left}" y="{height-20}">{lo:.3g}</text>')
    parts.append(f'<text class="tick" x="{width-right-60}" y="{height-20}">{hi:.3g}</text>')
    parts.append(f'<text class="label" x="{width/2:.0f}" y="{height-12}" text-anchor="middle">{html.escape(x_label)}</text>')
    parts.append(f'<text class="label" transform="translate(16 {height/2:.0f}) rotate(-90)" text-anchor="middle">{html.escape(y_label)}</text>')
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_scatter(x_values: list[float], y_values: list[float], path: Path, title: str, x_label: str, y_label: str) -> None:
    width, height = 760, 480
    left, right, top, bottom = 80, 40, 55, 70
    pairs = [(x, y) for x, y in zip(x_values, y_values) if x is not None and y is not None]
    if not pairs:
        path.write_text(svg_header(width, height) + f'<text x="20" y="40">{title}: no data</text></svg>\n', encoding="utf-8")
        return
    xs, ys = [p[0] for p in pairs], [p[1] for p in pairs]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if xmin == xmax:
        xmin -= 0.5
        xmax += 0.5
    if ymin == ymax:
        ymin -= 0.5
        ymax += 0.5
    plot_w = width - left - right
    plot_h = height - top - bottom
    parts = [svg_header(width, height), f'<text class="title" x="{left}" y="30">{html.escape(title)}</text>']
    parts.append(f'<line class="axis" x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}"/>')
    parts.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}"/>')
    for tick in range(6):
        x = left + tick * plot_w / 5
        y = height - bottom - tick * plot_h / 5
        x_value = xmin + tick * (xmax - xmin) / 5
        y_value = ymin + tick * (ymax - ymin) / 5
        parts.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{height-bottom}" stroke="#d0d5dd" stroke-width=".8"/>')
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}" stroke="#d0d5dd" stroke-width=".8"/>')
        parts.append(f'<text class="tick" x="{x:.1f}" y="{height-bottom+18}" text-anchor="middle">{x_value:.2g}</text>')
        parts.append(f'<text class="tick" x="{left-8}" y="{y+4:.1f}" text-anchor="end">{y_value:.2g}</text>')
    for x, y in pairs[:5000]:
        px = left + (x - xmin) / (xmax - xmin) * plot_w
        py = height - bottom - (y - ymin) / (ymax - ymin) * plot_h
        parts.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="2" fill="#E45756" opacity="0.5"/>')
    parts.append(f'<text class="label" x="{width/2:.0f}" y="{height-18}" text-anchor="middle">{html.escape(x_label)}</text>')
    parts.append(f'<text class="label" transform="translate(16 {height/2:.0f}) rotate(-90)" text-anchor="middle">{html.escape(y_label)}</text>')
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_bar_plot(items: list[tuple[str, int]], path: Path, title: str, x_label: str = "category", y_label: str = "count", limit: int = 20) -> None:
    width, height = 900, 520
    left, right, top, bottom = 78, 28, 54, 128
    plot_w = width - left - right
    plot_h = height - top - bottom
    data = items[:limit]
    max_count = max([count for _, count in data], default=1)
    parts = [svg_header(width, height), f'<text class="title" x="28" y="30">{html.escape(title)}</text>']
    parts.append(f'<line class="axis" x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}"/>')
    parts.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}"/>')
    for tick in range(5):
        y = top + plot_h - tick * plot_h / 4
        value = max_count * tick / 4
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" stroke="#d0d5dd" stroke-width=".8"/>')
        parts.append(f'<text class="tick" x="{left - 8}" y="{y + 4:.1f}" text-anchor="end">{value:.0f}</text>')
    if data:
        slot = plot_w / len(data)
        bar_w = max(8, min(44, slot * 0.72))
        colors = ["#4C78A8", "#E45756", "#54A24B", "#B279A2", "#F58518", "#72B7B2", "#9D755D"]
        for index, (name, count) in enumerate(data):
            x = left + index * slot + (slot - bar_w) / 2
            h = count / max_count * plot_h
            y = top + plot_h - h
            parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{colors[index % len(colors)]}"/>')
            parts.append(f'<text class="tick" x="{x + bar_w/2:.1f}" y="{y - 5:.1f}" text-anchor="middle">{count}</text>')
            parts.append(f'<text class="tick" transform="translate({x + bar_w/2:.1f} {top + plot_h + 12}) rotate(55)" text-anchor="start">{html.escape(name[:32])}</text>')
    parts.append(f'<text class="label" x="{width/2:.0f}" y="{height-12}" text-anchor="middle">{html.escape(x_label)}</text>')
    parts.append(f'<text class="label" transform="translate(16 {height/2:.0f}) rotate(-90)" text-anchor="middle">{html.escape(y_label)}</text>')
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_workflow_svg(path: Path) -> None:
    width, height = 1120, 420
    steps = [
        ("Select variants/elements", "GWAS, eQTL, cCRE, saturation mutagenesis"),
        ("Build MPRA library", "alleles, barcodes, controls, assay design"),
        ("QC barcode counts", "DNA input, RNA output, replicate balance"),
        ("Estimate activity/effects", "activity, log2FC, p/q values"),
        ("Integrate evidence", "chromatin, eQTL, CRISPRi, linked genes"),
        ("Prioritize biology", "causal variants, enhancers, target genes"),
    ]
    parts = [svg_header(width, height), '<text class="title" x="28" y="32">Literature-Informed MPRA Analysis Workflow</text>']
    box_w, box_h = 160, 92
    y = 130
    for index, (title, detail) in enumerate(steps):
        x = 38 + index * 178
        parts.append(f'<rect x="{x}" y="{y}" width="{box_w}" height="{box_h}" rx="6" fill="#f7f9fc" stroke="#98a2b3"/>')
        parts.append(f'<text x="{x+12}" y="{y+26}" font-size="13" font-weight="700">{html.escape(title)}</text>')
        parts.append(f'<foreignObject x="{x+12}" y="{y+38}" width="{box_w-24}" height="48"><div xmlns="http://www.w3.org/1999/xhtml" style="font-family:Arial;font-size:11px;color:#475467;line-height:1.25">{html.escape(detail)}</div></foreignObject>')
        if index < len(steps) - 1:
            parts.append(f'<path d="M {x+box_w+4} {y+box_h/2} L {x+174} {y+box_h/2}" stroke="#344054" stroke-width="1.5" marker-end="url(#arrow)"/>')
    parts.insert(2, '<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#344054"/></marker></defs>')
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_analysis_report(summary: dict[str, Any], plots: list[Path], output_csv: Path, label: str) -> Path:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    path = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_analysis_report.md"
    lines = ["# MPRA Analysis Report", "", f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}", "", f"Summary CSV: `{output_csv}`", ""]
    lines.append("## Summary Statistics")
    lines.append("")
    for column, values in summary["numeric_stats"].items():
        lines.append(f"- `{column}`: {values}")
    lines.append("")
    lines.append("## Counts")
    lines.append("")
    for key, value in summary["counts"].items():
        lines.append(f"- {key}: {value}")
    lines.append("")
    lines.append("## Plots")
    lines.append("")
    for plot in plots:
        lines.append(f"- `{plot}`")
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def run_portal_manifest(args: argparse.Namespace) -> int:
    summaries: list[dict[str, Any]] = []
    manifest: list[dict[str, str]] = []
    for accession, description in KNOWN_IGVF_MPRA_DATASETS:
        status, data, saved = fetch_json(PORTAL_BASE, f"/measurement-sets/{accession}/", f"portal_known_{accession}", {})
        summary = summarize_remote(f"portal_known_{accession}", description, status, data, saved)
        summaries.append(summary)
        manifest.extend(manifest_rows(f"portal_known_{accession}", description, data, limit=None))
        print(f"portal_known_{accession}: HTTP {status}, rows={len(rows_from_response(data))}")
    for label, params, description in PORTAL_MANY_MPRA_QUERIES:
        query = dict(params)
        query["limit"] = str(args.limit)
        status, data, saved = fetch_json(PORTAL_BASE, "/search/", label, query, portal_auth=True)
        if "starr" in label:
            filtered_data = filter_rows_by_terms(data, ("starr", "bluestarr"))
        elif "reporter" in label:
            filtered_data = filter_rows_by_terms(data, ("reporter", "mpra", "starr"))
        else:
            filtered_data = filter_rows_by_terms(data, ("mpra", "massively parallel reporter", "reporter"))
        summary = summarize_remote(label, description, status, filtered_data, saved)
        summaries.append(summary)
        manifest.extend(manifest_rows(label, description, filtered_data, limit=None))
        print(f"{label}: HTTP {status}, matched_rows={len(rows_from_response(filtered_data))}, reported_total={summary['total']}")
    manifest_path = write_manifest(manifest, f"{args.label}_portal_many")
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    query_counts = collections.Counter({summary["label"]: int(summary["returned_rows"]) for summary in summaries})
    status_counts = collections.Counter(str(row.get("status") or "unknown") for row in manifest)
    format_counts = collections.Counter(str(row.get("file_format") or row.get("output_type") or "metadata") for row in manifest)
    plots = [
        PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_portal_query_rows.svg",
        PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_portal_status_counts.svg",
        PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_portal_format_counts.svg",
    ]
    write_bar_plot(query_counts.most_common(), plots[0], "IGVF Portal MPRA Query Rows", "Portal query", "rows")
    write_bar_plot(status_counts.most_common(), plots[1], "IGVF Portal MPRA Status Counts", "status", "rows")
    write_bar_plot(format_counts.most_common(), plots[2], "IGVF Portal MPRA File/Output Types", "file or output type", "rows")
    save_json(f"mpra_{args.label}_portal_many_summary", summaries)
    report_path = write_portal_manifest_report(summaries, manifest_path, plots, args.label)
    print(f"Wrote manifest: {manifest_path}")
    print(f"Wrote report: {report_path}")
    for plot in plots:
        print(f"Wrote plot: {plot}")
    return 0 if any(200 <= item["http_status"] < 400 and item["returned_rows"] > 0 for item in summaries) else 1


def top_mpra_rows(rows: list[dict[str, str]], limit: int = 12) -> list[dict[str, str]]:
    def score(row: dict[str, str]) -> float:
        effect = abs(parse_float(row.get("MPRA.log2FoldChange")) or 0)
        sig = parse_float(row.get("MPRA.minusLog10PValue")) or 0
        activity = abs(parse_float(row.get("MPRA.max_activity")) or 0)
        return effect * 3 + sig + activity

    return sorted(rows, key=score, reverse=True)[:limit]


def markdown_table(rows: list[dict[str, str]], columns: list[str]) -> list[str]:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows:
        values = []
        for column in columns:
            value = row.get(column, "")
            if len(value) > 42:
                value = value[:39] + "..."
            values.append(html.escape(value))
        lines.append("| " + " | ".join(values) + " |")
    return lines


def run_literature_demo(args: argparse.Namespace) -> int:
    input_path = Path(args.input).resolve()
    rows = load_local_mpra(input_path)
    if not rows:
        print(f"No rows found in {input_path}", file=sys.stderr)
        return 1
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    label = safe_label(args.label)
    plots: list[Path] = []
    class_counts = collections.Counter(row.get("MPRA.class") or "unknown" for row in rows)
    sublib_counts = collections.Counter(row.get("MPRA.sublib") or row.get("libraryName") or "unknown" for row in rows)
    source_counts = collections.Counter(row.get("MPRA.source") or "unknown" for row in rows)
    active_sig_counts = collections.Counter()
    for row in rows:
        active = (row.get("MPRA.active") or "").upper() == "TRUE"
        signif = (row.get("MPRA.signif") or "").upper() == "TRUE"
        if active and signif:
            key = "active + significant"
        elif active:
            key = "active only"
        elif signif:
            key = "significant only"
        else:
            key = "neither"
        active_sig_counts[key] += 1
    plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}_mpra_class_counts.svg"
    write_bar_plot(class_counts.most_common(), plot, "MPRA Rows By Activity Class", "MPRA class", "variants")
    plots.append(plot)
    plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}_mpra_active_significant.svg"
    write_bar_plot(active_sig_counts.most_common(), plot, "MPRA Active And Significant Calls", "call category", "variants")
    plots.append(plot)
    plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}_mpra_sublib_counts.svg"
    write_bar_plot(sublib_counts.most_common(), plot, "MPRA Library/Sub-Library Counts", "library", "variants")
    plots.append(plot)
    plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}_mpra_source_counts.svg"
    write_bar_plot(source_counts.most_common(), plot, "MPRA Source Counts", "source", "variants")
    plots.append(plot)
    if "MPRA.max_activity" in rows[0] and "MPRA.log2FoldChange" in rows[0]:
        plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}_activity_vs_allelic_effect.svg"
        write_scatter(
            [parse_float(row.get("MPRA.max_activity")) for row in rows],
            [parse_float(row.get("MPRA.log2FoldChange")) for row in rows],
            plot,
            "MPRA Activity vs Allelic Effect",
            "max activity",
            "log2 fold-change",
        )
        plots.append(plot)
    if "MPRA.delta_log10_input" in rows[0] and "MPRA.log2FoldChange" in rows[0]:
        plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}_input_balance_vs_effect.svg"
        write_scatter(
            [parse_float(row.get("MPRA.delta_log10_input")) for row in rows],
            [abs(parse_float(row.get("MPRA.log2FoldChange")) or 0) for row in rows],
            plot,
            "MPRA Input Balance vs Absolute Allelic Effect",
            "delta log10 input",
            "absolute log2 fold-change",
        )
        plots.append(plot)
    use_case_counts = collections.Counter(item["application"] for item in MPRA_LITERATURE_USE_CASES)
    plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}_literature_use_cases.svg"
    write_bar_plot(use_case_counts.most_common(), plot, "MPRA Literature Use-Case Templates", "research application", "papers/examples")
    plots.append(plot)
    workflow = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}_literature_workflow.svg"
    write_workflow_svg(workflow)
    plots.append(workflow)
    top_rows = top_mpra_rows(rows, limit=12)
    summary = {
        "input": str(input_path),
        "rows": len(rows),
        "active_significant_counts": active_sig_counts.most_common(),
        "class_counts": class_counts.most_common(20),
        "sublib_counts": sublib_counts.most_common(20),
        "source_counts": source_counts.most_common(20),
        "literature_use_cases": MPRA_LITERATURE_USE_CASES,
        "plots": [str(plot) for plot in plots],
    }
    json_path = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}_literature_demo_summary.json"
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    report = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}_literature_informed_demo.md"
    lines = [
        "# Literature-Informed MPRA Research Demo",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"Input table: `{input_path}`",
        f"Summary JSON: `{json_path}`",
        "",
        "## What Recent MPRA Studies Use The Data For",
        "",
    ]
    for item in MPRA_LITERATURE_USE_CASES:
        lines.extend(
            [
                f"### {item['paper']} ({item['venue']}, {item['year']})",
                "",
                f"- Link: {item['url']}",
                f"- Research use: {item['application']}",
                f"- Plot templates to mimic: {item['plot_templates']}",
                f"- Skill rule: {item['lesson']}",
                "",
            ]
        )
    lines.extend(
        [
            "## Local Example Summary",
            "",
            f"- Rows: {len(rows):,}",
            f"- Active/significant categories: {', '.join(f'{k}: {v:,}' for k, v in active_sig_counts.most_common())}",
            f"- Top MPRA classes: {', '.join(f'{k}: {v:,}' for k, v in class_counts.most_common(8))}",
            "",
            "Top variants by combined activity/effect/significance score:",
            *markdown_table(top_rows, ["rsID", "SPDI", "MPRA.log2FoldChange", "MPRA.minusLog10PValue", "MPRA.class", "MPRA.nearestGene"]),
            "",
            "## Plots",
            "",
            *[f"- `{plot}`" for plot in plots],
            "",
            "## Interpretation Pattern",
            "",
            "Treat MPRA as direct regulatory activity evidence, then ask whether the active or allelically biased sequence overlaps cCREs/peaks, has a linked gene through rE2G/eQTL/CRISPRi, and falls in the right cell or tissue context. For disease loci, prioritize variants with strong MPRA effect, adequate barcode/input support, fine-mapping support, and concordant gene-linkage evidence.",
            "",
        ]
    )
    report.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote report: {report}")
    print(f"Wrote JSON: {json_path}")
    for plot in plots:
        print(f"Wrote plot: {plot}")
    return 0


def analyze_local(args: argparse.Namespace) -> int:
    input_path = Path(args.input).resolve()
    rows = load_local_mpra(input_path)
    mpra_columns = [column for column in rows[0].keys() if column.startswith("MPRA.")] if rows else []
    numeric_columns = [
        "MPRA.log2FoldChange",
        "MPRA.minusLog10PValue",
        "MPRA.minusLog10QValue",
        "MPRA.inputCountRef",
        "MPRA.outputCountRef",
        "MPRA.inputCountAlt",
        "MPRA.outputCountAlt",
        "MPRA.meanInput",
        "MPRA.max_activity",
        "MPRA.delta_log10_input",
    ]
    numeric_stats = {column: stats(numeric_column(rows, column)) for column in numeric_columns if column in rows[0]}
    counts = {
        "rows": len(rows),
        "rows_with_mpra_log2fc": len(numeric_column(rows, "MPRA.log2FoldChange")),
        "mpra_signif_TRUE": sum(1 for row in rows if (row.get("MPRA.signif") or "").upper() == "TRUE"),
        "mpra_active_TRUE": sum(1 for row in rows if (row.get("MPRA.active") or "").upper() == "TRUE"),
        "mpra_columns": len(mpra_columns),
    }
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    output_csv = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_summary_stats.csv"
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["column", "n", "mean", "median", "min", "max"])
        for column, values in numeric_stats.items():
            writer.writerow([column, values["n"], values["mean"], values["median"], values["min"], values["max"]])
    plots: list[Path] = []
    if "MPRA.log2FoldChange" in rows[0]:
        plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_log2fc_hist.svg"
        write_histogram(numeric_column(rows, "MPRA.log2FoldChange"), plot, "MPRA log2 fold-change")
        plots.append(plot)
    if "MPRA.minusLog10PValue" in rows[0]:
        plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_minuslog10p_hist.svg"
        write_histogram(numeric_column(rows, "MPRA.minusLog10PValue"), plot, "MPRA -log10 P-value")
        plots.append(plot)
    if "MPRA.log2FoldChange" in rows[0] and "MPRA.minusLog10PValue" in rows[0]:
        plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_volcano.svg"
        x = [parse_float(row.get("MPRA.log2FoldChange")) for row in rows]
        y = [parse_float(row.get("MPRA.minusLog10PValue")) for row in rows]
        write_scatter(x, y, plot, "MPRA volcano-style plot", "log2 fold-change", "-log10 P-value")
        plots.append(plot)
    if "MPRA.inputCountRef" in rows[0] and "MPRA.outputCountRef" in rows[0]:
        plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_ref_counts_scatter.svg"
        x = [parse_float(row.get("MPRA.inputCountRef")) for row in rows]
        y = [parse_float(row.get("MPRA.outputCountRef")) for row in rows]
        write_scatter(x, y, plot, "MPRA reference input/output counts", "input ref", "output ref")
        plots.append(plot)
    summary = {"counts": counts, "numeric_stats": numeric_stats, "input": str(input_path)}
    save_json(f"mpra_{args.label}_analysis_summary", summary)
    report = write_analysis_report(summary, plots, output_csv, args.label)
    print(f"Wrote stats CSV: {output_csv}")
    print(f"Wrote report: {report}")
    for plot in plots:
        print(f"Wrote plot: {plot}")
    return 0


def write_playbook() -> Path:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "MPRA_ANALYSIS_SKILLS.md"
    path.write_text(
        """# Skill: MPRA Data Retrieval And Analysis

Use this skill when the agent needs to retrieve MPRA/STARR/BlueSTARR evidence from IGVF Catalog, IGVF Portal, or ENCODE, or analyze local MPRA result tables.

## Retrieval

```bash
python3 Scripts/mpra_data_skills.py pull --source catalog --limit 25
python3 Scripts/mpra_data_skills.py pull --source all --limit 25
python3 Scripts/mpra_data_skills.py portal-manifest --limit 100 --label igvf_portal_mpra_many
```

Use `IGVF_PORTAL_COOKIE` for unreleased Portal datasets.

## Local Analysis

```bash
python3 Scripts/mpra_data_skills.py analyze-local --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra
python3 Scripts/mpra_data_skills.py literature-demo --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra_literature_demo
```

The analysis writes summary stats, SVG plots, and a Markdown report.

## Literature-Informed Templates

- Cell / PubMed variant-effect MPRA: use volcano plots, significant allelic-effect tables, and integration with GWAS/eQTL/fine-mapping.
- Nature Methods MPRA design benchmarking: use count QC, assay-design comparisons, activity dynamic range, and sequence-context summaries.
- Nature Genetics single-cell MPRA: use cell-type activity heatmaps and cluster/cell-type specificity summaries.
- Nature Biotechnology regulatory grammar MPRA: use saturation-mutagenesis maps, position-effect heatmaps, and model-vs-observed plots.
- Nature large-scale cCRE MPRA: use activity distributions, cCRE class enrichment, and genome-browser style views.

## Reuse Rules

- Inspect metadata before downloading large MPRA files.
- Preserve variant identifiers, allele orientation, library/source, biosample, activity class, input/output counts, log2 fold-change, P/Q values, and significance calls.
- Check input balance and low-count variants before interpreting effect sizes.
- Compare MPRA effects with IGVF Catalog variant-gene, variant-biosample, regulatory-element, and enhancer-gene prediction evidence.
- Use `portal-manifest` to pull many IGVF Portal MPRA/STARR/reporter datasets before choosing files to download.
- Use `literature-demo` to create paper-style interpretation plots and a research-use report from a local MPRA table.
""",
        encoding="utf-8",
    )
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="MPRA retrieval and analysis skills.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    pull = subparsers.add_parser("pull", help="Pull MPRA metadata/evidence.")
    pull.add_argument("--source", default="catalog", choices=["catalog", "portal", "encode", "all"])
    pull.add_argument("--limit", type=int, default=25)
    portal = subparsers.add_parser("portal-manifest", help="Pull many IGVF Portal MPRA/STARR/reporter datasets into a manifest with plots.")
    portal.add_argument("--limit", type=int, default=100)
    portal.add_argument("--label", default="igvf_portal_mpra_many")
    analyze = subparsers.add_parser("analyze-local", help="Analyze a local MPRA result table.")
    analyze.add_argument("--input", default=str(DATA_DIR / "Input" / "VariantList" / "example_variants.csv"))
    analyze.add_argument("--label", default="mpra_local")
    literature = subparsers.add_parser("literature-demo", help="Create literature-informed MPRA plots and interpretation report.")
    literature.add_argument("--input", default=str(DATA_DIR / "Input" / "VariantList" / "example_variants.csv"))
    literature.add_argument("--label", default="mpra_literature_demo")
    subparsers.add_parser("write-playbook", help="Write MPRA skill documentation.")

    args = parser.parse_args(argv)
    setup_logging()
    if args.command == "pull":
        return run_pull(args)
    if args.command == "portal-manifest":
        return run_portal_manifest(args)
    if args.command == "analyze-local":
        return analyze_local(args)
    if args.command == "literature-demo":
        return run_literature_demo(args)
    if args.command == "write-playbook":
        path = write_playbook()
        print(f"Wrote {path}")
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
