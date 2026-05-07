#!/usr/bin/env python3
"""CRISPRi data retrieval and functional-annotation integration skills."""

from __future__ import annotations

import argparse
import collections
import csv
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
REPORT_DIR = DOCS_DIR / "CRISPRi"
PLOT_DIR = REPORT_DIR / "Plots"
MANIFEST_DIR = DATA_DIR / "Manifests" / "CRISPRi"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

CATALOG_API_BASE = os.environ.get(
    "IGVF_CATALOG_API_BASE", "https://api.catalogkg.igvf.org"
).rstrip("/")
PORTAL_BASE = os.environ.get("IGVF_PORTAL_BASE", "https://data.igvf.org").rstrip("/")
ENCODE_BASE = os.environ.get("ENCODE_BASE", "https://www.encodeproject.org").rstrip("/")


CATALOG_CRISPRI_QUERIES = [
    (
        "catalog_crispr_facs_element_gene",
        "/api/genomic-elements/genes",
        {"method": "CRISPR FACS screen", "limit": "25", "page": "0", "verbose": "false"},
        "Catalog genomic element-gene links from CRISPR FACS screen evidence.",
    ),
    (
        "catalog_perturb_seq_gene_element",
        "/api/genes/genomic-elements",
        {"method": "Perturb-seq", "limit": "25", "page": "0", "verbose": "false"},
        "Catalog gene-element links from Perturb-seq evidence.",
    ),
    (
        "catalog_crispr_region_element_gene",
        "/api/genomic-elements/genes",
        {"region": "chr1:903900-904900", "limit": "25", "page": "0", "verbose": "false"},
        "Region-based regulatory element-gene evidence useful for CRISPRi linkage review.",
    ),
]

PORTAL_CRISPRI_QUERIES = [
    ("portal_crispri_files", {"type": "File", "searchTerm": "CRISPRi", "limit": "25"}, "IGVF Portal files matching CRISPRi."),
    ("portal_crispri_measurement_sets", {"type": "MeasurementSet", "searchTerm": "CRISPRi", "limit": "25"}, "IGVF Portal measurement sets matching CRISPRi."),
    ("portal_crispri_analysis_sets", {"type": "AnalysisSet", "searchTerm": "CRISPRi", "limit": "25"}, "IGVF Portal analysis sets matching CRISPRi."),
    ("portal_crispr_facs_files", {"type": "File", "searchTerm": "CRISPR FACS", "limit": "25"}, "IGVF Portal files matching CRISPR FACS."),
]

ENCODE_CRISPRI_QUERIES = [
    ("encode_crispri_experiments", {"type": "Experiment", "searchTerm": "CRISPRi", "limit": "25"}, "ENCODE experiments matching CRISPRi."),
    ("encode_crispr_screen_files", {"type": "File", "searchTerm": "CRISPR screen", "limit": "25"}, "ENCODE files matching CRISPR screen."),
    ("encode_guide_quantifications", {"type": "File", "output_type": "guide quantifications", "limit": "25"}, "Guide quantification files relevant to perturbation screens."),
]


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"crispri_data_skills_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
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
        for key in ("accession", "@id", "href", "method", "source", "name", "title", "biological_context", "gene_name", "gene_id"):
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
        "genes": count_field(rows, ("gene", "target_gene", "gene_name")),
        "file_formats": count_field(rows, ("file_format", "output_type", "content_type")),
    }


def manifest_rows(label: str, description: str, data: Any) -> list[dict[str, str]]:
    rows = []
    for row in rows_from_response(data)[:10]:
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
                "gene": text_value(row.get("gene") or row.get("target_gene") or row.get("gene_name")),
                "file_format": text_value(row.get("file_format")),
                "output_type": text_value(row.get("output_type") or row.get("content_type")),
                "href": text_value(row.get("href")),
                "source_url": text_value(row.get("source_url")),
            }
        )
    return rows


def write_manifest(rows: list[dict[str, str]], label: str) -> Path:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    path = MANIFEST_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_manifest.csv"
    fields = ["label", "description", "accession", "id", "status", "method", "source", "context", "gene", "file_format", "output_type", "href", "source_url"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


def format_counter(items: list[tuple[str, int]]) -> str:
    return "\n".join(f"- {name}: {count}" for name, count in items) if items else "- none observed"


def write_remote_report(summaries: list[dict[str, Any]], manifest_path: Path, label: str) -> Path:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    path = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_metadata_report.md"
    lines = ["# CRISPRi Metadata Retrieval Report", "", f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}", "", f"Manifest: `{manifest_path}`", ""]
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
                "Top genes:",
                format_counter(summary["genes"]),
                "",
                "Top file/output types:",
                format_counter(summary["file_formats"]),
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def run_pull(args: argparse.Namespace) -> int:
    summaries: list[dict[str, Any]] = []
    manifest: list[dict[str, str]] = []
    sources = ["catalog", "portal", "encode"] if args.source == "all" else [args.source]
    for source in sources:
        if source == "catalog":
            for label, path, params, description in CATALOG_CRISPRI_QUERIES:
                query = dict(params)
                query["limit"] = str(args.limit)
                status, data, saved = fetch_json(CATALOG_API_BASE, path, label, query)
                summaries.append(summarize_remote(label, description, status, data, saved))
                manifest.extend(manifest_rows(label, description, data))
                print(f"{label}: HTTP {status}, rows={len(rows_from_response(data))}")
        elif source == "portal":
            for label, params, description in PORTAL_CRISPRI_QUERIES:
                query = dict(params)
                query["limit"] = str(args.limit)
                status, data, saved = fetch_json(PORTAL_BASE, "/search/", label, query, portal_auth=True)
                summaries.append(summarize_remote(label, description, status, data, saved))
                manifest.extend(manifest_rows(label, description, data))
                print(f"{label}: HTTP {status}, rows={len(rows_from_response(data))}")
        elif source == "encode":
            for label, params, description in ENCODE_CRISPRI_QUERIES:
                query = dict(params)
                query["limit"] = str(args.limit)
                status, data, saved = fetch_json(ENCODE_BASE, "/search/", label, query)
                summaries.append(summarize_remote(label, description, status, data, saved))
                manifest.extend(manifest_rows(label, description, data))
                print(f"{label}: HTTP {status}, rows={len(rows_from_response(data))}")
    save_json(f"crispri_{args.source}_summary", summaries)
    manifest_path = write_manifest(manifest, f"crispri_{args.source}")
    report_path = write_remote_report(summaries, manifest_path, f"crispri_{args.source}")
    print(f"Wrote manifest: {manifest_path}")
    print(f"Wrote report: {report_path}")
    return 0 if any(200 <= item["http_status"] < 400 and item["returned_rows"] > 0 for item in summaries) else 1


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


def numeric_column(rows: list[dict[str, str]], column: str) -> list[float]:
    return [value for value in (parse_float(row.get(column)) for row in rows) if value is not None]


def stats(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"n": 0, "mean": None, "median": None, "min": None, "max": None}
    return {"n": len(values), "mean": statistics.fmean(values), "median": statistics.median(values), "min": min(values), "max": max(values)}


def boolish(row: dict[str, str], column: str) -> bool:
    return (row.get(column) or "").strip().upper() in {"1", "TRUE", "YES"}


def svg_header(width: int, height: int) -> str:
    return f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n<rect width="100%" height="100%" fill="#fbfaf7"/>\n'


def escape_xml(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


def tick_values(lo: float, hi: float, count: int = 5) -> list[float]:
    if hi == lo:
        return [lo]
    step = (hi - lo) / count
    return [lo + step * index for index in range(count + 1)]


def write_histogram(values: list[float], path: Path, title: str, bins: int = 30) -> None:
    width, height = 920, 520
    left, right, top, bottom = 92, 42, 72, 86
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
    axis_y = height - bottom
    parts = [
        svg_header(width, height),
        f'<text x="{left}" y="32" font-size="20" font-weight="700" font-family="Arial" fill="#1f2933">{escape_xml(title)}</text>',
        f'<text x="{left}" y="54" font-size="12" font-family="Arial" fill="#5b6470">n={len(values)}, min={lo:.3g}, max={hi:.3g}; x axis is binned value, y axis is variant count.</text>',
        f'<line x1="{left}" y1="{axis_y}" x2="{width-right}" y2="{axis_y}" stroke="#222"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{axis_y}" stroke="#222"/>',
        f'<text x="{width/2-58:.0f}" y="{height-28}" font-size="13" font-family="Arial" fill="#1f2933">binned value</text>',
        f'<text x="20" y="{height/2+45:.0f}" font-size="13" font-family="Arial" fill="#1f2933" transform="rotate(-90 20 {height/2+45:.0f})">variant count</text>',
    ]
    for tick in tick_values(lo, hi):
        x = left + (tick - lo) / (hi - lo) * plot_w
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{axis_y}" stroke="#e0e4e8"/>')
        parts.append(f'<line x1="{x:.2f}" y1="{axis_y}" x2="{x:.2f}" y2="{axis_y+5}" stroke="#222"/>')
        parts.append(f'<text x="{x-18:.2f}" y="{axis_y+22}" font-size="11" font-family="Arial" fill="#47515c">{tick:.2g}</text>')
    for tick in tick_values(0, max_count):
        y = axis_y - (tick / max_count * plot_h if max_count else 0)
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width-right}" y2="{y:.2f}" stroke="#eef1f3"/>')
        parts.append(f'<line x1="{left-5}" y1="{y:.2f}" x2="{left}" y2="{y:.2f}" stroke="#222"/>')
        parts.append(f'<text x="{left-44}" y="{y+4:.2f}" font-size="11" font-family="Arial" fill="#47515c">{tick:.0f}</text>')
    for i, count in enumerate(counts):
        bar_h = count / max_count * plot_h
        x = left + i * bar_w
        y = height - bottom - bar_h
        parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{max(bar_w-1, 1):.2f}" height="{bar_h:.2f}" fill="#3d7a9c" opacity="0.9"/>')
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_scatter(x_values: list[float | None], y_values: list[float | None], path: Path, title: str, x_label: str, y_label: str) -> None:
    width, height = 920, 560
    left, right, top, bottom = 96, 44, 76, 92
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
    axis_y = height - bottom
    parts = [
        svg_header(width, height),
        f'<text x="{left}" y="32" font-size="20" font-weight="700" font-family="Arial" fill="#1f2933">{escape_xml(title)}</text>',
        f'<text x="{left}" y="54" font-size="12" font-family="Arial" fill="#5b6470">n={len(pairs)} paired variants; each point is one variant with both measurements.</text>',
        f'<line x1="{left}" y1="{axis_y}" x2="{width-right}" y2="{axis_y}" stroke="#222"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{axis_y}" stroke="#222"/>',
    ]
    for tick in tick_values(xmin, xmax):
        px = left + (tick - xmin) / (xmax - xmin) * plot_w
        parts.append(f'<line x1="{px:.2f}" y1="{top}" x2="{px:.2f}" y2="{axis_y}" stroke="#e0e4e8"/>')
        parts.append(f'<line x1="{px:.2f}" y1="{axis_y}" x2="{px:.2f}" y2="{axis_y+5}" stroke="#222"/>')
        parts.append(f'<text x="{px-18:.2f}" y="{axis_y+22}" font-size="11" font-family="Arial" fill="#47515c">{tick:.2g}</text>')
    for tick in tick_values(ymin, ymax):
        py = axis_y - (tick - ymin) / (ymax - ymin) * plot_h
        parts.append(f'<line x1="{left}" y1="{py:.2f}" x2="{width-right}" y2="{py:.2f}" stroke="#eef1f3"/>')
        parts.append(f'<line x1="{left-5}" y1="{py:.2f}" x2="{left}" y2="{py:.2f}" stroke="#222"/>')
        parts.append(f'<text x="{left-56}" y="{py+4:.2f}" font-size="11" font-family="Arial" fill="#47515c">{tick:.2g}</text>')
    if xmin < 0 < xmax:
        x0 = left + (0 - xmin) / (xmax - xmin) * plot_w
        parts.append(f'<line x1="{x0:.2f}" y1="{top}" x2="{x0:.2f}" y2="{axis_y}" stroke="#9099a1" stroke-dasharray="4 4"/>')
    if ymin < 0 < ymax:
        y0 = axis_y - (0 - ymin) / (ymax - ymin) * plot_h
        parts.append(f'<line x1="{left}" y1="{y0:.2f}" x2="{width-right}" y2="{y0:.2f}" stroke="#9099a1" stroke-dasharray="4 4"/>')
    for x, y in pairs[:5000]:
        px = left + (x - xmin) / (xmax - xmin) * plot_w
        py = height - bottom - (y - ymin) / (ymax - ymin) * plot_h
        color = "#c84444" if abs(x) >= 1 or abs(y) >= 1 else "#4d7898"
        parts.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="2.3" fill="{color}" opacity="0.55"/>')
    parts.append(f'<text x="{width/2-90:.0f}" y="{height-30}" font-size="13" font-family="Arial" fill="#1f2933">{escape_xml(x_label)}</text>')
    parts.append(f'<text x="22" y="{height/2+55:.0f}" font-size="13" font-family="Arial" fill="#1f2933" transform="rotate(-90 22 {height/2+55:.0f})">{escape_xml(y_label)}</text>')
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_count_bar(items: list[tuple[str, int]], path: Path, title: str, x_label: str = "variant count") -> None:
    width = 920
    top, right, bottom, left = 72, 110, 76, 285
    row_h = 34
    height = max(300, top + bottom + row_h * max(1, len(items)))
    max_count = max((count for _, count in items), default=1)
    plot_w = width - left - right
    axis_y = height - bottom
    parts = [
        svg_header(width, height),
        f'<text x="{left}" y="32" font-size="20" font-weight="700" font-family="Arial" fill="#1f2933">{escape_xml(title)}</text>',
        f'<text x="{left}" y="54" font-size="12" font-family="Arial" fill="#5b6470">Bars show how many variants satisfy each evidence flag.</text>',
        f'<line x1="{left}" y1="{axis_y}" x2="{width-right}" y2="{axis_y}" stroke="#222"/>',
        f'<line x1="{left}" y1="{top-8}" x2="{left}" y2="{axis_y}" stroke="#222"/>',
        f'<text x="{width/2-45:.0f}" y="{height-24}" font-size="13" font-family="Arial" fill="#1f2933">{escape_xml(x_label)}</text>',
    ]
    for tick in tick_values(0, max_count):
        x = left + (tick / max_count if max_count else 0) * plot_w
        parts.append(f'<line x1="{x:.2f}" y1="{top-8}" x2="{x:.2f}" y2="{axis_y}" stroke="#e0e4e8"/>')
        parts.append(f'<text x="{x-14:.2f}" y="{axis_y+22}" font-size="11" font-family="Arial" fill="#47515c">{tick:.0f}</text>')
    for index, (name, count) in enumerate(items):
        y = top + index * row_h
        bar_w = count / max_count * plot_w if max_count else 0
        parts.append(f'<text x="24" y="{y+21}" font-size="12" font-family="Arial" fill="#28323c">{escape_xml(name[:38])}</text>')
        parts.append(f'<rect x="{left}" y="{y+5}" width="{bar_w:.2f}" height="22" rx="2" fill="#2f776d"/>')
        parts.append(f'<text x="{left + bar_w + 8:.2f}" y="{y+21}" font-size="12" font-family="Arial" fill="#28323c">{count}</text>')
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def integrated_score(row: dict[str, str]) -> float:
    score = 0.0
    if boolish(row, "CRISPRi.LDL.uptake.significant..1.TRUE..0.tested.and.FALSE."):
        score += 3
    crispri_abs = parse_float(row.get("CRISPRi.LDL.uptake.abs.mu_Z"))
    if crispri_abs is not None:
        score += min(crispri_abs, 5) / 5 * 2
    if boolish(row, "MPRA.signif"):
        score += 1.5
    if boolish(row, "in_HepG2_ATAC.seq_peak"):
        score += 1
    if (row.get("SCREEN.cCRE..Nov.2025.") or "").strip() not in {"", "NA", "not_within_cCRE"}:
        score += 1
    if boolish(row, "FAVOR_Predicted_Functional"):
        score += 1
    if boolish(row, "BE.targets.ref..1.TRUE..0.FALSE."):
        score += 1
    return score


def analyze_local(args: argparse.Namespace) -> int:
    input_path = Path(args.input).resolve()
    rows = load_csv(input_path)
    if not rows:
        print(f"No rows found in {input_path}")
        return 2
    for row in rows:
        row["CRISPRi.functional_annotation_score"] = f"{integrated_score(row):.4f}"
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    annotated_path = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_crispri_functional_annotation.csv"
    fields = list(rows[0].keys())
    with annotated_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    numeric_columns = [
        "CRISPRi.LDL.uptake.mu_Z",
        "CRISPRi.LDL.uptake.mu_Z.synced.for.MPRA..flipped.if.LOF.increases.LDL.uptake.",
        "CRISPRi.LDL.uptake.abs.mu_Z",
        "BE.LDL.uptake.abs.mu_Z_adj",
        "MPRA.log2FoldChange",
        "MPRA.minusLog10PValue",
        "CRISPRi.functional_annotation_score",
    ]
    numeric_stats = {column: stats(numeric_column(rows, column)) for column in numeric_columns if column in rows[0]}
    counts = {
        "rows": len(rows),
        "crispri_tested": len(numeric_column(rows, "CRISPRi.LDL.uptake.mu_Z")),
        "crispri_significant": sum(1 for row in rows if boolish(row, "CRISPRi.LDL.uptake.significant..1.TRUE..0.tested.and.FALSE.")),
        "mpra_significant": sum(1 for row in rows if boolish(row, "MPRA.signif")),
        "in_hepg2_atac_peak": sum(1 for row in rows if boolish(row, "in_HepG2_ATAC.seq_peak")),
        "in_ccre": sum(1 for row in rows if (row.get("SCREEN.cCRE..Nov.2025.") or "").strip() not in {"", "NA", "not_within_cCRE"}),
        "favor_predicted_functional": sum(1 for row in rows if boolish(row, "FAVOR_Predicted_Functional")),
    }
    plots: list[Path] = []
    plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_crispri_abs_mu_z_hist.svg"
    write_histogram(numeric_column(rows, "CRISPRi.LDL.uptake.abs.mu_Z"), plot, "CRISPRi absolute LDL uptake mu Z")
    plots.append(plot)
    plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_functional_score_hist.svg"
    write_histogram(numeric_column(rows, "CRISPRi.functional_annotation_score"), plot, "CRISPRi integrated functional annotation score")
    plots.append(plot)
    plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_crispri_vs_mpra.svg"
    write_scatter(
        [parse_float(row.get("CRISPRi.LDL.uptake.mu_Z")) for row in rows],
        [parse_float(row.get("MPRA.log2FoldChange")) for row in rows],
        plot,
        "CRISPRi effect vs MPRA log2FC",
        "CRISPRi LDL uptake mu Z",
        "MPRA log2FC",
    )
    plots.append(plot)
    plot = PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_evidence_counts.svg"
    write_count_bar(
        [
            ("CRISPRi tested", counts["crispri_tested"]),
            ("CRISPRi significant", counts["crispri_significant"]),
            ("MPRA significant", counts["mpra_significant"]),
            ("HepG2 ATAC peak", counts["in_hepg2_atac_peak"]),
            ("SCREEN cCRE", counts["in_ccre"]),
            ("FAVOR predicted functional", counts["favor_predicted_functional"]),
        ],
        plot,
        "Functional Evidence Overlap Counts",
    )
    plots.append(plot)
    stats_path = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_summary_stats.csv"
    with stats_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["column", "n", "mean", "median", "min", "max"])
        for column, values in numeric_stats.items():
            writer.writerow([column, values["n"], values["mean"], values["median"], values["min"], values["max"]])
    summary = {"input": str(input_path), "annotated_csv": str(annotated_path), "counts": counts, "numeric_stats": numeric_stats, "plots": [str(plot) for plot in plots]}
    save_json(f"crispri_{args.label}_analysis_summary", summary)
    report_path = write_analysis_report(summary, stats_path, args.label)
    print(f"Wrote annotated CSV: {annotated_path}")
    print(f"Wrote stats CSV: {stats_path}")
    print(f"Wrote report: {report_path}")
    for plot_path in plots:
        print(f"Wrote plot: {plot_path}")
    return 0


def write_analysis_report(summary: dict[str, Any], stats_path: Path, label: str) -> Path:
    path = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_analysis_report.md"
    lines = ["# CRISPRi Functional Annotation Integration Report", "", f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}", "", f"Stats CSV: `{stats_path}`", f"Annotated CSV: `{summary['annotated_csv']}`", ""]
    lines.append("## Counts")
    lines.append("")
    for key, value in summary["counts"].items():
        lines.append(f"- {key}: {value}")
    lines.append("")
    lines.append("## Summary Statistics")
    lines.append("")
    for column, values in summary["numeric_stats"].items():
        lines.append(f"- `{column}`: {values}")
    lines.append("")
    lines.append("## Plots")
    lines.append("")
    for plot in summary["plots"]:
        lines.append(f"- `{plot}`")
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_playbook() -> Path:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "CRISPRI_ANALYSIS_SKILLS.md"
    path.write_text(
        """# Skill: CRISPRi Data Access, Download Planning, And Functional Annotation Integration

Use this skill when the agent needs CRISPRi, CRISPR FACS screen, or Perturb-seq regulatory perturbation evidence from IGVF Portal, IGVF Catalog, or ENCODE.

## Bioinformatics Research Use Cases

CRISPRi data are most useful when the question needs endogenous perturbation evidence rather than correlation-only regulatory annotation.

1. **Functional cCRE discovery and benchmarking.** Noncoding CRISPRi screens can directly test candidate cis-regulatory elements and compare hits against cCRE, ATAC/DNase, H3K27ac, TF-binding, MPRA, and rE2G annotations. The ENCODE multicenter analysis assembled 108 noncoding CRISPR screens, >540,000 perturbations, and 24.85 Mb of tested sequence, then benchmarked analysis tools and screen design rules.
2. **Enhancer-to-gene mapping.** CRISPRi-FlowFISH, CRISPRi growth screens, and TAP-seq/Perturb-seq can connect perturbed enhancers to target genes. TAP-seq demonstrated perturbation-based enhancer-target maps for 1,778 enhancers.
3. **Variant-to-gene interpretation.** Variants inside CRISPRi-supported regulatory elements should be prioritized when the same locus also has cCRE overlap, accessibility, MPRA allelic effect, rE2G/single-cell linkage, eQTL/QTL, or disease fine-mapping evidence.
4. **Single-cell perturbation programs and GRNs.** Perturb-seq and compressed Perturb-seq use CRISPRi/CRISPR perturbations plus single-cell RNA-seq readouts to infer regulatory circuits, cell states, pathway programs, and genetic interactions.
5. **Multiome perturbation interpretation.** Multiome Perturb-seq extends CRISPRi screens to paired gene expression and chromatin accessibility, enabling analysis of how perturbations alter both transcription and regulatory DNA state.
6. **Screen QC and method comparison.** CRISPRi screen analysis should inspect guide efficiency, low-specificity guides, strand/orientation artifacts, replicate concordance, negative controls, local chromatin context, and whether effects are expression, growth, reporter, or cell-state readouts.

## Literature Anchors

- Multicenter integrated analysis of noncoding CRISPRi screens, Nature Methods 2024: https://www.nature.com/articles/s41592-024-02216-7
- Targeted Perturb-seq enables genome-scale genetic screens in single cells, Nature Methods 2020: https://www.nature.com/articles/s41592-020-0837-5
- Scalable genetic screening for regulatory circuits using compressed Perturb-seq, Nature Biotechnology 2024: https://www.nature.com/articles/s41587-023-01964-9
- Multiome Perturb-seq unlocks scalable discovery of integrated perturbation effects on the transcriptome and epigenome, Cell Systems 2025 / PubMed: https://pubmed.ncbi.nlm.nih.gov/39091800/
- Enhlink infers distal and context-specific enhancer-promoter linkages, Genome Biology 2024: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03374-9

## Retrieval

```bash
python3 Scripts/crispri_data_skills.py pull --source catalog --limit 25
python3 Scripts/crispri_data_skills.py pull --source portal --limit 25
python3 Scripts/crispri_data_skills.py pull --source all --limit 25
```

Use `IGVF_PORTAL_COOKIE` for unreleased IGVF Portal datasets.

## Functional Annotation Integration

```bash
python3 Scripts/crispri_data_skills.py analyze-local --input Data/Input/VariantList/example_variants.csv --label my_locus_crispri
```

The integration score combines CRISPRi significance/effect size, MPRA significance, HepG2 ATAC peak overlap, cCRE annotation, FAVOR predicted function, and base-editing target evidence.

## Suggested IGVFdataAgent Workflow

1. Pull CRISPRi/Perturb-seq metadata from Catalog, Portal, and ENCODE.
2. Build a download manifest before fetching large count, guide, or processed result files.
3. Normalize guide-level or element-level output to stable columns: element interval, guide ID, target gene, biosample, method, effect size, p/q value, direction, source file, and screen readout.
4. Join with cCRE classes, ENCODE-rE2G links, single-cell linkage files, MPRA effects, eQTL/QTL, ATAC/DNase/H3K27ac peaks, and FAVOR/IGVF variant annotations.
5. Summarize hit rates by cCRE class, cell type, target gene, distance-to-TSS, and evidence overlap.
6. Plot effect-size distributions, guide concordance, CRISPRi-vs-MPRA agreement, evidence overlap counts, and IGV-like browser views around prioritized loci.
7. Report research interpretation: candidate causal CREs, likely target genes, cell contexts, phenotype direction, and follow-up experiments.

## Reuse Rules

- Preserve guide/element, target gene, biosample, assay, source fileset, effect size, significance, and direction fields.
- Check whether CRISPRi direction has been synchronized to MPRA or gene loss-of-function assumptions before combining evidence.
- Integrate CRISPRi with MPRA, ATAC/cCRE, eQTL/QTL, enhancer-gene prediction, and base-editing evidence.
- Use manifests before downloading large files.
- Treat CRISPRi evidence as context-specific: a negative result in one cell type does not rule out regulatory function in another context.
""",
        encoding="utf-8",
    )
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="CRISPRi retrieval and analysis skills.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    pull = subparsers.add_parser("pull", help="Pull CRISPRi metadata/evidence.")
    pull.add_argument("--source", default="catalog", choices=["catalog", "portal", "encode", "all"])
    pull.add_argument("--limit", type=int, default=25)
    analyze = subparsers.add_parser("analyze-local", help="Integrate local CRISPRi columns with functional annotation evidence.")
    analyze.add_argument("--input", default=str(DATA_DIR / "Input" / "VariantList" / "example_variants.csv"))
    analyze.add_argument("--label", default="crispri_local")
    subparsers.add_parser("write-playbook", help="Write CRISPRi skill documentation.")

    args = parser.parse_args(argv)
    setup_logging()
    if args.command == "pull":
        return run_pull(args)
    if args.command == "analyze-local":
        return analyze_local(args)
    if args.command == "write-playbook":
        path = write_playbook()
        print(f"Wrote {path}")
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
