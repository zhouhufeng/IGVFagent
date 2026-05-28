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


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "MPRA"
PLOT_DIR = REPORT_DIR / "Plots"
MANIFEST_DIR = DATA_DIR / "Manifests" / "MPRA"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint

CATALOG_API_BASE = _resolve_endpoint("catalog_api", "IGVF_CATALOG_API_BASE")
PORTAL_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_BASE")
ENCODE_BASE = _resolve_endpoint("encode", "ENCODE_BASE")


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


# ---------------------------------------------------------------------------
# MPRASuite / esMPRA-style core analytics
# ---------------------------------------------------------------------------
# Clean-room reimplementations of the algorithms used in the canonical MPRA
# pipelines (Tewhey-lab MPRASuite, Apache-2.0; WangLabTHU esMPRA). No code is
# copied -- algorithms paraphrased from the published descriptions:
#   * Per-oligo activity = log2(RNA/DNA), modeled as a negative-binomial GLM
#     (DESeq2 Wald test). After fitting, size factors for the RNA samples are
#     shifted so the mode of the log2-fold-change density lies at 0
#     ("summit-shift" normalization, MPRAmodel).
#   * Allelic skew = paired t-test of the per-replicate log2-ratio between
#     reference and alternate alleles of the same SNP-window-strand-haplotype
#     element. BH-FDR controls multiple-testing.
#   * Replicate concordance = pairwise Pearson/Spearman correlation between
#     normalized RNA counts across replicates.
#
# Heavy numerical deps (numpy, pandas, scipy, pydeseq2, statsmodels,
# matplotlib) are imported lazily so the lightweight pull / portal-manifest /
# write-playbook commands keep working when those packages are absent.


def _require_pkg(name: str, hint: str) -> Any:
    try:
        return __import__(name)
    except Exception as exc:  # pragma: no cover - defensive
        raise SystemExit(
            f"Missing dependency '{name}'. {hint}\nInstall with: pip install {name}"
        ) from exc


def _detect_sample_columns(columns: list[str]) -> dict[str, list[str]]:
    """Group columns by DNA/RNA condition.

    Accepts columns named ``DNA_rep1``, ``DNA.r1``, ``DNA1``, ``RNA_HepG2_rep2``,
    etc. Returns ``{"DNA": [...], "RNA": [...]}``.
    """
    import re
    pat = re.compile(r"^(DNA|RNA)([._-].*)?$", re.IGNORECASE)
    buckets: dict[str, list[str]] = {"DNA": [], "RNA": []}
    for column in columns:
        match = pat.match(column)
        if match:
            buckets[match.group(1).upper()].append(column)
    return buckets


def _load_counts_table(path: Path) -> tuple[Any, dict[str, list[str]]]:
    """Load a counts table into a pandas DataFrame and detect sample columns."""
    pd = _require_pkg("pandas", "Needed for MPRA analytics.")
    df = pd.read_csv(path, sep=None, engine="python")
    if "Oligo" not in df.columns:
        raise SystemExit(
            "Counts table must contain an 'Oligo' column. "
            "Expected schema: Oligo,Allele,...,DNA_rep1,...,RNA_rep1,..."
        )
    samples = _detect_sample_columns(list(df.columns))
    if not samples["DNA"] or not samples["RNA"]:
        raise SystemExit(
            "Could not find DNA_* and RNA_* sample columns. "
            f"Detected: {samples}. Rename columns to DNA_rep1, RNA_rep1, etc."
        )
    if "Barcode" in df.columns:
        keep_meta = [c for c in df.columns if c not in samples["DNA"] + samples["RNA"] + ["Barcode"]]
        meta = df.groupby("Oligo", as_index=False)[keep_meta].first()
        counts = df.groupby("Oligo", as_index=False)[samples["DNA"] + samples["RNA"]].sum()
        df = meta.merge(counts, on="Oligo", how="inner")
    return df, samples


def _summit_shift(log2fc: Any) -> float:
    """Return the shift (in log2 units) that places the mode of ``log2fc`` at 0."""
    np = _require_pkg("numpy", "Needed for MPRA analytics.")
    scipy_stats = __import__("scipy.stats", fromlist=["stats"])
    values = np.asarray([v for v in log2fc if v is not None and not np.isnan(v)], dtype=float)
    if values.size < 8:
        return 0.0
    kde = scipy_stats.gaussian_kde(values)
    grid = np.linspace(float(np.percentile(values, 1)), float(np.percentile(values, 99)), 2048)
    return float(grid[int(np.argmax(kde(grid)))])


def _run_deseq_activity(df: Any, samples: dict[str, list[str]]) -> Any:
    """Run per-oligo NB GLM Wald (RNA vs DNA) with summit-shift normalization."""
    pd = _require_pkg("pandas", "Needed for MPRA analytics.")
    np = _require_pkg("numpy", "Needed for MPRA analytics.")
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    sample_cols = samples["DNA"] + samples["RNA"]
    counts = df.set_index("Oligo")[sample_cols].T.astype(int)
    metadata = pd.DataFrame(
        {"condition": ["DNA"] * len(samples["DNA"]) + ["RNA"] * len(samples["RNA"])},
        index=sample_cols,
    )
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design_factors="condition",
        refit_cooks=False,
    )
    dds.deseq2()
    res = DeseqStats(dds, contrast=["condition", "RNA", "DNA"])
    res.summary()
    res_df = res.results_df.copy()
    shift = _summit_shift(res_df["log2FoldChange"].to_numpy())
    if abs(shift) > 1e-6:
        sf = dds.obsm["size_factors"].copy()
        mask = (metadata["condition"] == "RNA").to_numpy()
        sf[mask] *= 2 ** shift
        dds.obsm["size_factors"] = sf
        try:
            dds.fit_genewise_dispersions()
            dds.fit_dispersion_trend()
            dds.fit_dispersion_prior()
            dds.fit_MAP_dispersions()
            dds.fit_LFC()
            dds.calculate_cooks()
        except Exception:
            pass
        res = DeseqStats(dds, contrast=["condition", "RNA", "DNA"])
        res.summary()
        res_df = res.results_df.copy()
    res_df.attrs["summit_shift"] = shift
    return res_df


def cmd_activity(args: argparse.Namespace) -> int:
    """Compute per-oligo MPRA activity (log2FC RNA/DNA) with NB GLM + summit-shift."""
    pd = _require_pkg("pandas", "Needed for MPRA analytics.")
    setup_logging()
    input_path = Path(args.input).resolve()
    df, samples = _load_counts_table(input_path)
    logging.info(
        "Loaded counts: %d oligos, %d DNA + %d RNA samples",
        len(df), len(samples["DNA"]), len(samples["RNA"]),
    )
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    res_df = _run_deseq_activity(df, samples)
    out_csv = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_activity.out"
    res_df.to_csv(out_csv, sep="\t", index_label="Oligo")
    padj = pd.to_numeric(res_df["padj"], errors="coerce")
    summary = {
        "input": str(input_path),
        "oligos": int(len(res_df)),
        "summit_shift_log2": float(res_df.attrs.get("summit_shift", 0.0)),
        "active_fdr_0.05": int((padj < 0.05).sum()),
        "active_fdr_0.10": int((padj < 0.10).sum()),
        "median_log2fc": float(pd.to_numeric(res_df["log2FoldChange"], errors="coerce").median()),
    }
    save_json(f"mpra_{args.label}_activity_summary", summary)
    print(f"Wrote activity table: {out_csv}")
    print(f"Active @ FDR 0.05: {summary['active_fdr_0.05']} / {summary['oligos']}")
    print(f"Active @ FDR 0.10: {summary['active_fdr_0.10']} / {summary['oligos']}")
    print(f"Summit-shift applied (log2): {summary['summit_shift_log2']:+.4f}")
    return 0


def _pair_alleles(df: Any) -> tuple[Any, Any, list[str]]:
    """Pair ref/alt oligos by SNP_window_strand_haplotype."""
    cols = [c for c in ("SNP", "Window", "Strand", "Haplotype") if c in df.columns]
    if not cols:
        raise SystemExit(
            "Skew test needs at least one of SNP/Window/Strand/Haplotype columns "
            "to pair ref vs alt oligos."
        )
    if "Allele" not in df.columns:
        raise SystemExit("Skew test needs an 'Allele' column (REF/ALT or A/B).")
    df = df.copy()
    df["__key"] = df[cols].astype(str).agg("_".join, axis=1)
    allele_norm = df["Allele"].astype(str).str.upper().map(
        {"REF": "REF", "ALT": "ALT", "A": "REF", "B": "ALT"}
    ).fillna(df["Allele"].astype(str).str.upper())
    df["__allele"] = allele_norm
    ref = df[df["__allele"] == "REF"].set_index("__key")
    alt = df[df["__allele"] == "ALT"].set_index("__key")
    common = sorted(set(ref.index) & set(alt.index))
    return ref.loc[common], alt.loc[common], common


def cmd_skew(args: argparse.Namespace) -> int:
    """Compute allelic skew via paired t-test on log2(RNA/DNA) ratios."""
    pd = _require_pkg("pandas", "Needed for MPRA analytics.")
    np = _require_pkg("numpy", "Needed for MPRA analytics.")
    scipy_stats = __import__("scipy.stats", fromlist=["stats"])
    from statsmodels.stats.multitest import multipletests
    setup_logging()
    input_path = Path(args.input).resolve()
    df, samples = _load_counts_table(input_path)
    ref, alt, keys = _pair_alleles(df)
    logging.info("Paired %d ref/alt oligos", len(keys))
    if not keys:
        raise SystemExit("No ref/alt pairs found. Check Allele and pairing columns.")
    dna_cols = samples["DNA"]
    rna_cols = samples["RNA"]
    pseudo = 1.0
    dna_mean_ref = ref[dna_cols].astype(float).mean(axis=1).replace(0, np.nan)
    dna_mean_alt = alt[dna_cols].astype(float).mean(axis=1).replace(0, np.nan)
    n_rep = len(rna_cols)
    log2_ref = np.log2(
        (ref[rna_cols].astype(float).to_numpy() + pseudo)
        / (dna_mean_ref.to_numpy()[:, None] + pseudo)
    )
    log2_alt = np.log2(
        (alt[rna_cols].astype(float).to_numpy() + pseudo)
        / (dna_mean_alt.to_numpy()[:, None] + pseudo)
    )
    log2_skew = log2_alt - log2_ref
    if n_rep > 1:
        se = log2_skew.std(axis=1, ddof=1) / np.sqrt(n_rep)
        tstat, pvals = scipy_stats.ttest_rel(log2_alt, log2_ref, axis=1, nan_policy="omit")
    else:
        se = np.full(len(keys), np.nan)
        tstat = np.full(len(keys), np.nan)
        pvals = np.full(len(keys), np.nan)
    pvals_arr = np.asarray(pvals, dtype=float)
    valid = ~np.isnan(pvals_arr)
    padj = np.full_like(pvals_arr, np.nan)
    if valid.any():
        _, padj_valid, _, _ = multipletests(pvals_arr[valid], method="fdr_bh")
        padj[valid] = padj_valid
    out_df = pd.DataFrame({
        "Element": keys,
        "Log2Skew": log2_skew.mean(axis=1),
        "LogSkew_SE": se,
        "tstat": tstat,
        "pvalue": pvals_arr,
        "padj": padj,
        "RefOligo": ref["Oligo"].to_numpy() if "Oligo" in ref.columns else ref.index,
        "AltOligo": alt["Oligo"].to_numpy() if "Oligo" in alt.columns else alt.index,
    })
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out_csv = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_skew.out"
    out_df.to_csv(out_csv, sep="\t", index=False)
    padj_series = pd.to_numeric(out_df["padj"], errors="coerce")
    summary = {
        "input": str(input_path),
        "pairs": int(len(out_df)),
        "skew_fdr_0.05": int((padj_series < 0.05).sum()),
        "skew_fdr_0.10": int((padj_series < 0.10).sum()),
        "median_abs_log2skew": float(out_df["Log2Skew"].abs().median()),
    }
    save_json(f"mpra_{args.label}_skew_summary", summary)
    print(f"Wrote skew table: {out_csv}")
    print(f"Allelic skew @ FDR 0.05: {summary['skew_fdr_0.05']} / {summary['pairs']}")
    print(f"Allelic skew @ FDR 0.10: {summary['skew_fdr_0.10']} / {summary['pairs']}")
    return 0


def cmd_qc(args: argparse.Namespace) -> int:
    """Replicate concordance + barcodes/oligo + counts/oligo histograms."""
    pd = _require_pkg("pandas", "Needed for MPRA analytics.")
    np = _require_pkg("numpy", "Needed for MPRA analytics.")
    setup_logging()
    input_path = Path(args.input).resolve()
    raw = pd.read_csv(input_path, sep=None, engine="python")
    samples = _detect_sample_columns(list(raw.columns))
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    plots: list[Path] = []
    if "Oligo" not in raw.columns:
        raise SystemExit("QC needs an 'Oligo' column.")
    if "Barcode" in raw.columns:
        bc_per = raw.groupby("Oligo")["Barcode"].nunique().to_numpy()
        bc_path = PLOT_DIR / f"{ts}_{safe_label(args.label)}_qc_barcodes_per_oligo.svg"
        write_histogram(
            [float(x) for x in bc_per],
            bc_path,
            "Unique barcodes per oligo",
            x_label="barcodes / oligo",
        )
        plots.append(bc_path)
    sample_cols = samples["DNA"] + samples["RNA"]
    per_oligo = raw.groupby("Oligo")[sample_cols].sum() if sample_cols else None
    if per_oligo is not None:
        total_counts = per_oligo.sum(axis=1).to_numpy()
        counts_path = PLOT_DIR / f"{ts}_{safe_label(args.label)}_qc_counts_per_oligo.svg"
        write_histogram(
            [float(np.log10(x + 1)) for x in total_counts],
            counts_path,
            "Total counts per oligo (log10)",
            x_label="log10 (counts + 1)",
        )
        plots.append(counts_path)
    corr_summary: dict[str, Any] = {}
    if per_oligo is not None:
        for cond in ("DNA", "RNA"):
            cols = samples[cond]
            if len(cols) < 2:
                continue
            log_counts = np.log10(per_oligo[cols].to_numpy() + 1.0)
            corr_matrix = np.corrcoef(log_counts.T)
            tri = corr_matrix[np.triu_indices_from(corr_matrix, k=1)]
            corr_summary[cond] = {
                "samples": cols,
                "min_pearson": float(tri.min()),
                "median_pearson": float(np.median(tri)),
                "max_pearson": float(tri.max()),
            }
            heat_path = PLOT_DIR / f"{ts}_{safe_label(args.label)}_qc_repcorr_{cond}.svg"
            _write_corr_heatmap(corr_matrix, cols, heat_path, f"{cond} replicate Pearson r")
            plots.append(heat_path)
    summary = {
        "input": str(input_path),
        "n_oligos": int(raw["Oligo"].nunique()),
        "samples": samples,
        "replicate_concordance": corr_summary,
    }
    save_json(f"mpra_{args.label}_qc_summary", summary)
    report = REPORT_DIR / f"{ts}_{safe_label(args.label)}_qc_report.md"
    lines = [
        "# MPRA QC Report", "",
        f"Input: `{input_path}`",
        f"Oligos: {summary['n_oligos']}",
        f"DNA samples: {len(samples['DNA'])}",
        f"RNA samples: {len(samples['RNA'])}",
        "",
        "## Replicate concordance (Pearson r on log10 counts)", "",
    ]
    for cond, info in corr_summary.items():
        lines.append(
            f"- {cond}: median r = {info['median_pearson']:.3f}, "
            f"min = {info['min_pearson']:.3f}, max = {info['max_pearson']:.3f} "
            f"(n = {len(info['samples'])})"
        )
    lines.append("")
    lines.append("## Plots")
    lines.append("")
    for plot in plots:
        lines.append(f"- `{plot}`")
    report.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote QC report: {report}")
    for plot in plots:
        print(f"Wrote plot: {plot}")
    return 0


def _write_corr_heatmap(matrix: Any, labels: list[str], path: Path, title: str) -> None:
    """Tiny dependency-free correlation heatmap (SVG, blue ramp)."""
    n = len(labels)
    cell = max(38, min(56, 420 // max(n, 1)))
    margin_left, margin_top = 110, 60
    width = margin_left + cell * n + 80
    height = margin_top + cell * n + 90
    parts = [svg_header(width, height)]
    parts.append(f'<text class="title" x="20" y="30">{html.escape(title)}</text>')
    for i in range(n):
        for j in range(n):
            value = float(matrix[i][j])
            intensity = max(0.0, min(1.0, (value + 1) / 2))
            red = int(247 - 200 * intensity)
            green = int(251 - 180 * intensity)
            blue = int(255 - 100 * intensity)
            x = margin_left + j * cell
            y = margin_top + i * cell
            parts.append(
                f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="rgb({red},{green},{blue})" stroke="#cbd2da" stroke-width=".5"/>'
            )
            parts.append(
                f'<text x="{x + cell/2:.1f}" y="{y + cell/2 + 4:.1f}" text-anchor="middle" font-size="10" fill="#1f2933">{value:.2f}</text>'
            )
    for i, label in enumerate(labels):
        parts.append(
            f'<text class="tick" x="{margin_left - 8}" y="{margin_top + i * cell + cell/2 + 4:.1f}" text-anchor="end">{html.escape(label)}</text>'
        )
        parts.append(
            f'<text class="tick" transform="translate({margin_left + i * cell + cell/2:.1f} {margin_top + n*cell + 14}) rotate(45)" text-anchor="start">{html.escape(label)}</text>'
        )
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def _write_volcano_4panel(
    activity: list[tuple],
    skew: list[tuple],
    path: Path,
    title: str,
    fdr: float = 0.05,
) -> None:
    """4-panel volcano: activity volcano + skew volcano + activity MA + skew MA."""
    width, height = 1200, 920
    parts = [svg_header(width, height)]
    parts.append(f'<text class="title" x="28" y="30">{html.escape(title)}</text>')
    panels = [
        ("Activity volcano", "log2(RNA/DNA)", "-log10 padj", activity, 30, 70),
        ("Allelic-skew volcano", "Log2Skew (alt - ref)", "-log10 padj", skew, 620, 70),
        ("Activity MA", "log10 mean RNA+DNA", "log2(RNA/DNA)", activity, 30, 500),
        ("Skew MA", "log10 mean RNA", "Log2Skew", skew, 620, 500),
    ]
    plot_w, plot_h = 540, 360
    for sub_title, x_label, y_label, points, ox, oy in panels:
        pts = [(x, y, q) for x, y, q in points if x is not None and y is not None]
        parts.append(f'<rect x="{ox-10}" y="{oy-30}" width="{plot_w+30}" height="{plot_h+70}" fill="none" stroke="#cbd2da"/>')
        parts.append(f'<text x="{ox}" y="{oy - 6}" font-size="13" font-weight="700">{html.escape(sub_title)}</text>')
        if not pts:
            parts.append(f'<text x="{ox + 20}" y="{oy + 40}" fill="#666">no data</text>')
            continue
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        if xmin == xmax: xmin, xmax = xmin - 0.5, xmax + 0.5
        if ymin == ymax: ymin, ymax = ymin - 0.5, ymax + 0.5
        parts.append(f'<line class="axis" x1="{ox}" y1="{oy+plot_h}" x2="{ox+plot_w}" y2="{oy+plot_h}"/>')
        parts.append(f'<line class="axis" x1="{ox}" y1="{oy}" x2="{ox}" y2="{oy+plot_h}"/>')
        for tick in range(5):
            tx = ox + tick * plot_w / 4
            ty = oy + plot_h - tick * plot_h / 4
            xv = xmin + tick * (xmax - xmin) / 4
            yv = ymin + tick * (ymax - ymin) / 4
            parts.append(f'<line x1="{tx:.1f}" y1="{oy}" x2="{tx:.1f}" y2="{oy+plot_h}" stroke="#e6eaf0" stroke-width=".7"/>')
            parts.append(f'<line x1="{ox}" y1="{ty:.1f}" x2="{ox+plot_w}" y2="{ty:.1f}" stroke="#e6eaf0" stroke-width=".7"/>')
            parts.append(f'<text class="tick" x="{tx:.1f}" y="{oy+plot_h+16}" text-anchor="middle">{xv:.2g}</text>')
            parts.append(f'<text class="tick" x="{ox-8}" y="{ty+4:.1f}" text-anchor="end">{yv:.2g}</text>')
        for x, y, q in pts[:6000]:
            px = ox + (x - xmin) / (xmax - xmin) * plot_w
            py = oy + plot_h - (y - ymin) / (ymax - ymin) * plot_h
            color = "#E45756" if (q is not None and q < fdr) else "#8c98a7"
            opacity = "0.85" if (q is not None and q < fdr) else "0.35"
            parts.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="2.4" fill="{color}" opacity="{opacity}"/>')
        parts.append(f'<text class="label" x="{ox + plot_w/2:.0f}" y="{oy+plot_h+34}" text-anchor="middle">{html.escape(x_label)}</text>')
        parts.append(f'<text class="label" transform="translate({ox-44} {oy+plot_h/2:.0f}) rotate(-90)" text-anchor="middle">{html.escape(y_label)}</text>')
    parts.append(f'<text class="tick" x="{width - 220}" y="{height - 18}">red = padj &lt; {fdr:g}</text>')
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def cmd_volcano(args: argparse.Namespace) -> int:
    """4-panel volcano (activity + skew + MA plots) from existing .out tables."""
    pd = _require_pkg("pandas", "Needed for MPRA analytics.")
    np = _require_pkg("numpy", "Needed for MPRA analytics.")
    setup_logging()
    activity_pts: list[tuple] = []
    skew_pts: list[tuple] = []
    if args.activity:
        act = pd.read_csv(args.activity, sep="\t")
        for _, row in act.iterrows():
            lfc = row.get("log2FoldChange")
            padj = row.get("padj")
            if pd.notna(lfc) and pd.notna(padj):
                activity_pts.append((float(lfc), float(-np.log10(max(float(padj), 1e-300))), float(padj)))
    if args.skew:
        sk = pd.read_csv(args.skew, sep="\t")
        for _, row in sk.iterrows():
            log2s = row.get("Log2Skew"); padj = row.get("padj")
            if pd.notna(log2s) and pd.notna(padj):
                skew_pts.append((float(log2s), float(-np.log10(max(float(padj), 1e-300))), float(padj)))
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    out = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_volcano4.svg"
    _write_volcano_4panel(activity_pts, skew_pts, out, args.title or "MPRA activity & allelic skew", fdr=args.fdr)
    print(f"Wrote 4-panel volcano: {out}")
    return 0


def write_playbook() -> Path:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "MPRA_ANALYSIS_SKILLS.md"
    path.write_text(
        """# Skill: MPRA Data Retrieval And Analysis

Use this skill when the agent needs to retrieve MPRA/STARR/BlueSTARR evidence
from IGVF Catalog, IGVF Portal, or ENCODE, or to run end-to-end MPRA
analytics (per-element activity, allelic skew, QC) on a local counts table.

The analytical methods are clean-room reimplementations of the canonical
public pipelines (algorithms only, no source copied):

- **Tewhey-lab MPRASuite** (Apache-2.0): `MPRAmodel` per-element NB GLM
  Wald test (DESeq2) of RNA vs DNA with **summit-shift** size-factor
  normalization, then paired-t allelic skew across replicates with BH-FDR.
- **WangLabTHU esMPRA**: count-based activity ratios, replicate concordance,
  barcode-per-oligo QC.

## Retrieval

```bash
python3 Scripts/mpra_data_skills.py pull --source catalog --limit 25
python3 Scripts/mpra_data_skills.py pull --source all --limit 25
python3 Scripts/mpra_data_skills.py portal-manifest --limit 100 --label igvf_portal_mpra_many
```

Use `IGVF_PORTAL_COOKIE` for unreleased Portal datasets.

## Local Result-Table Summaries

```bash
python3 Scripts/mpra_data_skills.py analyze-local --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra
python3 Scripts/mpra_data_skills.py literature-demo --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra_literature_demo
```

These take an existing MPRA result table (with columns like
`MPRA.log2FoldChange`, `MPRA.minusLog10PValue`) and emit summary stats,
SVG plots and a Markdown report.

## End-to-End Counts Analytics (NEW)

Run these on a counts table with columns `Oligo,Allele,SNP,Window,Strand,
Haplotype,DNA_rep1,DNA_rep2,...,RNA_rep1,RNA_rep2,...` (barcode-level
tables with a `Barcode` column are summed per oligo automatically):

```bash
# QC: replicate Pearson r heatmaps, BC/oligo, log10 counts/oligo
python3 Scripts/mpra_data_skills.py qc       --input counts.tsv --label myrun

# Per-oligo activity (NB GLM Wald + summit-shift) -> *_activity.out
python3 Scripts/mpra_data_skills.py activity --input counts.tsv --label myrun

# Allelic skew (paired t-test of ALT vs REF on log2(RNA/DNA)) -> *_skew.out
python3 Scripts/mpra_data_skills.py skew     --input counts.tsv --label myrun

# 4-panel volcano: activity, skew, MA-activity, MA-skew (red = padj < FDR)
python3 Scripts/mpra_data_skills.py volcano  \
    --activity Docs/MPRA/<ts>_myrun_activity.out \
    --skew     Docs/MPRA/<ts>_myrun_skew.out      \
    --label myrun --fdr 0.05
```

Dependencies (installed once into the project venv): `pandas`, `numpy`,
`scipy`, `statsmodels`, `pydeseq2` (BSD-3).

### Output schemas

- **`*_activity.out`** (TSV) -- one row per Oligo, columns:
  `baseMean, log2FoldChange, lfcSE, stat, pvalue, padj`. The
  `summit_shift_log2` value used to recenter size factors is stored in
  the JSON summary saved next to the table.
- **`*_skew.out`** (TSV) -- one row per paired element, columns:
  `Element, Log2Skew, LogSkew_SE, tstat, pvalue, padj, RefOligo, AltOligo`.
- **`*_qc_report.md`** -- per-condition replicate Pearson stats (median,
  min, max r) and links to BC/oligo + counts/oligo + Pearson heatmap SVGs.

## Literature-Informed Templates

- Cell / PubMed variant-effect MPRA: volcano plots, significant allelic-
  effect tables, integration with GWAS/eQTL/fine-mapping.
- Nature Methods MPRA design benchmarking: count QC, assay-design
  comparisons, activity dynamic range, sequence-context summaries.
- Nature Genetics single-cell MPRA: cell-type activity heatmaps and
  cluster/cell-type specificity summaries.
- Nature Biotechnology regulatory grammar MPRA: saturation-mutagenesis
  maps, position-effect heatmaps, and model-vs-observed plots.
- Nature large-scale cCRE MPRA: activity distributions, cCRE class
  enrichment, genome-browser style views.

## Reuse Rules

- Inspect metadata before downloading large MPRA files.
- For raw counts: always run `qc` first to check replicate concordance
  before trusting activity / skew calls.
- Use `summit-shift` (default in `activity`) to center the log2FC mode
  on zero before declaring "active" elements.
- Preserve variant identifiers, allele orientation, library/source,
  biosample, activity class, counts, log2FC, P/Q values, significance.
- Cross-check active or skewed elements against IGVF Catalog variant-gene,
  variant-biosample, regulatory-element, and enhancer-gene predictions.
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
    pull.add_argument("--label", default="mpra_pull",
                       help="Output label; defaults to 'mpra_pull'. "
                             "Accepted (and used by save_json) so the tool "
                             "dispatcher's --label flag isn't a hard error.")
    portal = subparsers.add_parser("portal-manifest", help="Pull many IGVF Portal MPRA/STARR/reporter datasets into a manifest with plots.")
    portal.add_argument("--limit", type=int, default=100)
    portal.add_argument("--label", default="igvf_portal_mpra_many")
    analyze = subparsers.add_parser("analyze-local", help="Analyze a local MPRA result table.")
    analyze.add_argument("--input", default=str(DATA_DIR / "Input" / "VariantList" / "example_variants.csv"))
    analyze.add_argument("--label", default="mpra_local")
    literature = subparsers.add_parser("literature-demo", help="Create literature-informed MPRA plots and interpretation report.")
    literature.add_argument("--input", default=str(DATA_DIR / "Input" / "VariantList" / "example_variants.csv"))
    literature.add_argument("--label", default="mpra_literature_demo")

    activity = subparsers.add_parser(
        "activity",
        help="Per-oligo MPRA activity via NB GLM Wald test (DESeq2) with summit-shift normalization.",
    )
    activity.add_argument("--input", required=True, help="Counts table: Oligo,Allele,...,DNA_rep1,...,RNA_rep1,...")
    activity.add_argument("--label", default="mpra_activity")

    skew = subparsers.add_parser(
        "skew",
        help="Allelic-skew paired t-test on log2(RNA/DNA) per replicate, with BH-FDR.",
    )
    skew.add_argument("--input", required=True, help="Counts table with Allele + at least one of SNP/Window/Strand/Haplotype.")
    skew.add_argument("--label", default="mpra_skew")

    qc = subparsers.add_parser(
        "qc",
        help="Replicate concordance, barcodes/oligo, counts/oligo histograms.",
    )
    qc.add_argument("--input", required=True, help="Counts table (barcode-level allowed via 'Barcode' column).")
    qc.add_argument("--label", default="mpra_qc")

    volcano = subparsers.add_parser(
        "volcano",
        help="4-panel volcano (activity + allelic skew + MA plots) from .out tables.",
    )
    volcano.add_argument("--activity", help="Activity .out table from `mpra activity`.")
    volcano.add_argument("--skew", help="Skew .out table from `mpra skew`.")
    volcano.add_argument("--label", default="mpra_volcano")
    volcano.add_argument("--title", default=None)
    volcano.add_argument("--fdr", type=float, default=0.05)

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
    if args.command == "activity":
        return cmd_activity(args)
    if args.command == "skew":
        return cmd_skew(args)
    if args.command == "qc":
        return cmd_qc(args)
    if args.command == "volcano":
        return cmd_volcano(args)
    if args.command == "write-playbook":
        path = write_playbook()
        print(f"Wrote {path}")
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
