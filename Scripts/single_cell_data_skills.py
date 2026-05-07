#!/usr/bin/env python3
"""Reusable skills for scRNA-seq, scATAC-seq, and Perturb-seq datasets.

This script searches IGVF Portal and ENCODE metadata, writes small manifests,
downloads safe public examples, runs lightweight standard-library analyses,
and documents analysis-ready inputs for Scanpy/ArchR/Signac, Cell Ranger
outputs, or perturbation-analysis workflows.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import logging
import math
import os
import random
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
SKILL_DOC_DIR = DOCS_DIR / "Skills"
SINGLE_CELL_DOC_DIR = DOCS_DIR / "SingleCell"
SINGLE_CELL_PLOT_DIR = SINGLE_CELL_DOC_DIR / "Plots"
MANIFEST_DIR = DATA_DIR / "Manifests" / "SingleCell"
EXAMPLE_DIR = DATA_DIR / "SingleCell" / "Examples"

PORTAL_BASE = os.environ.get("IGVF_PORTAL_BASE", "https://api.data.igvf.org").rstrip("/")
ENCODE_BASE = os.environ.get("ENCODE_BASE", "https://www.encodeproject.org").rstrip("/")

EXAMPLE_FILES = [
    {
        "modality": "scRNA/scE2G",
        "accession": "IGVFFI9905RPTO",
        "url": "https://api.data.igvf.org/tabular-files/IGVFFI9905RPTO/@@download/IGVFFI9905RPTO.tsv.gz",
        "path": EXAMPLE_DIR / "IGVFFI9905RPTO.tsv.gz",
        "description": "K562 CRISPRi 10x Multiome scE2G gene quantification table.",
    },
    {
        "modality": "scATAC",
        "accession": "IGVFFI6024DNYI",
        "url": "https://api.data.igvf.org/tabular-files/IGVFFI6024DNYI/@@download/IGVFFI6024DNYI.bed.gz",
        "path": EXAMPLE_DIR / "IGVFFI6024DNYI.bed.gz",
        "description": "Public IGVF GRCm39 single-cell ATAC peak BED example.",
    },
]


SKILLS = {
    "scrna": {
        "title": "single-cell RNA-seq",
        "purpose": "Quantify gene expression across cells/nuclei and compare cell states, tissues, donors, perturbations, or conditions.",
        "preferred_files": [
            "gene-barcode matrix",
            "h5ad",
            "matrix hdf5",
            "fragments are not expected for RNA",
        ],
        "encode_queries": [
            {"type": "File", "output_type": "gene quantifications"},
            {"type": "File", "file_format": "h5ad"},
        ],
        "portal_queries": [
            {"type": "MeasurementSet", "assay_titles": "single-cell RNA-seq"},
            {"type": "MeasurementSet", "assay_titles": "single-nucleus RNA-seq"},
            {"type": "File", "content_type": "gene quantifications"},
            {"type": "File", "file_format": "h5ad"},
        ],
        "analysis_steps": [
            "Build a sample sheet from metadata: accession, biosample, donor, organism, assembly, assay, file accession, file format, and download URL.",
            "Load count matrices into AnnData; preserve raw counts before normalization.",
            "Run QC for library size, detected genes, mitochondrial/ribosomal fraction, doublets, and batch labels.",
            "Normalize/log-transform, select highly variable genes, integrate batches when needed, cluster, annotate cell types, and run differential expression.",
            "For IGVF variant/gene interpretation, summarize expression of candidate genes in relevant cell types and conditions.",
        ],
    },
    "scatac": {
        "title": "single-cell ATAC-seq",
        "purpose": "Profile chromatin accessibility at cell resolution and connect variants to accessible elements and target genes.",
        "preferred_files": [
            "fragments",
            "peak calls",
            "cell-by-peak matrix",
            "bigWig signal",
            "h5ad or loom when available",
        ],
        "encode_queries": [
            {"type": "File", "output_type": "fragments"},
            {"type": "File", "file_format": "bed", "output_type": "peaks"},
        ],
        "portal_queries": [
            {"type": "MeasurementSet", "assay_titles": "single-cell ATAC-seq"},
            {"type": "MeasurementSet", "assay_titles": "single-nucleus ATAC-seq"},
            {"type": "File", "content_type": "fragments"},
            {"type": "File", "content_type": "peaks"},
        ],
        "analysis_steps": [
            "Build a manifest for fragments, peak files, cell metadata, genome assembly, and biosample ontology.",
            "Run QC for fragments per cell, TSS enrichment, fraction in peaks, blacklist fraction, nucleosome signal, and doublets.",
            "Create a cell-by-peak matrix, run TF-IDF/LSI, integrate batches, cluster, and annotate cell types.",
            "Intersect variant lists with peaks/cCREs and summarize accessibility around candidate loci.",
            "Link peaks to genes using co-accessibility, nearby genes, Catalog element-gene predictions, or matched scRNA data.",
        ],
    },
    "perturbseq": {
        "title": "Perturb-seq",
        "purpose": "Analyze pooled perturbation screens with single-cell expression readouts.",
        "preferred_files": [
            "gene expression matrix",
            "guide assignment table",
            "perturbation metadata",
            "cell metadata",
            "h5ad when available",
        ],
        "encode_queries": [
            {"type": "File", "output_type": "guide quantifications"},
            {"type": "File", "file_format": "h5ad"},
        ],
        "portal_queries": [
            {"type": "MeasurementSet", "assay_titles": "Perturb-seq"},
            {"type": "AnalysisSet", "assay_titles": "Perturb-seq"},
            {"type": "File", "content_type": "guide quantifications"},
            {"type": "File", "file_format": "h5ad"},
        ],
        "analysis_steps": [
            "Build a manifest linking expression matrices, guide assignments, perturbation targets, controls, donors, and conditions.",
            "QC cells, guides, perturbation multiplicity, control cells, and batch labels.",
            "Assign perturbation status per cell; remove ambiguous or high-multiplicity cells when needed.",
            "Run differential expression and pathway/module scoring per perturbation against matched controls.",
            "Use IGVF Catalog genes, regulatory elements, and variant links to interpret perturbation targets and downstream effects.",
        ],
    },
}


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"single_cell_data_skills_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)
    return log_path


def safe_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in label)


def build_url(base: str, path: str, params: dict[str, Any]) -> str:
    url = f"{base}{path}"
    return f"{url}?{urllib.parse.urlencode(params, doseq=True)}"


def save_json(label: str, data: Any) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}.json"
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    logging.info("Saved JSON: %s", path)
    return path


def fetch_search(base: str, label: str, params: dict[str, Any], *, portal: bool = False) -> tuple[int, Any, Path]:
    request_params = dict(params)
    request_params["format"] = "json"
    url = build_url(base, "/search/", request_params)
    headers = {
        "Accept": "application/json,*/*",
        "User-Agent": "IGVFdataAgent/0.1",
    }
    if portal and os.environ.get("IGVF_PORTAL_COOKIE"):
        headers["Cookie"] = os.environ["IGVF_PORTAL_COOKIE"]
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
        for key in ("title", "name", "term_name", "accession", "@id", "uuid", "href"):
            if key in value:
                return scalar_strings(value[key])
    return []


def count_field(rows: list[dict[str, Any]], field_names: tuple[str, ...], limit: int = 10) -> list[tuple[str, int]]:
    counter: collections.Counter[str] = collections.Counter()
    for row in rows:
        for field in field_names:
            for value in scalar_strings(row.get(field)):
                if value:
                    counter[value] += 1
    return counter.most_common(limit)


def text_value(value: Any) -> str:
    values = scalar_strings(value)
    return ";".join(values[:5])


def summarize_rows(source: str, skill: str, query: dict[str, Any], status: int, data: Any, saved: Path) -> dict[str, Any]:
    rows = rows_from_response(data)
    return {
        "source": source,
        "skill": skill,
        "query": query,
        "http_status": status,
        "saved_response": str(saved),
        "returned_rows": len(rows),
        "total": total_from_response(data, rows),
        "assays": count_field(rows, ("assay_title", "assay_titles", "assay_term_name")),
        "biosamples": count_field(rows, ("biosample_summary", "biosample_term_name", "biosample_ontology", "sample_terms")),
        "file_formats": count_field(rows, ("file_format", "file_format_type")),
        "output_types": count_field(rows, ("output_type", "content_type", "content_types")),
        "statuses": count_field(rows, ("status",)),
        "examples": rows[:5],
    }


def manifest_rows(source: str, skill: str, summary: dict[str, Any]) -> list[dict[str, str]]:
    rows = []
    for row in summary.get("examples", []):
        rows.append(
            {
                "source": source,
                "skill": skill,
                "query": json.dumps(summary.get("query", {}), sort_keys=True),
                "accession": text_value(row.get("accession")),
                "id": text_value(row.get("@id") or row.get("uuid")),
                "status": text_value(row.get("status")),
                "assay": text_value(row.get("assay_title") or row.get("assay_titles") or row.get("assay_term_name")),
                "biosample": text_value(row.get("biosample_summary") or row.get("biosample_term_name") or row.get("biosample_ontology")),
                "file_format": text_value(row.get("file_format") or row.get("file_format_type")),
                "output_type": text_value(row.get("output_type") or row.get("content_type") or row.get("content_types")),
                "assembly": text_value(row.get("assembly")),
                "href": text_value(row.get("href")),
                "summary": text_value(row.get("summary") or row.get("description") or row.get("title")),
            }
        )
    return rows


def format_counter(items: list[tuple[str, int]]) -> str:
    if not items:
        return "- none observed"
    return "\n".join(f"- {name}: {count}" for name, count in items)


def skill_names(selection: str) -> list[str]:
    if selection == "all":
        return list(SKILLS)
    if selection not in SKILLS:
        raise ValueError(f"Unknown skill: {selection}")
    return [selection]


def run_searches(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    sources = ["encode", "portal"] if args.source == "both" else [args.source]
    summaries: list[dict[str, Any]] = []
    manifests: list[dict[str, str]] = []
    for skill_name in skill_names(args.skill):
        skill = SKILLS[skill_name]
        for source in sources:
            query_key = "encode_queries" if source == "encode" else "portal_queries"
            base = ENCODE_BASE if source == "encode" else PORTAL_BASE
            for index, query in enumerate(skill[query_key]):
                request_query = dict(query)
                request_query["limit"] = str(args.limit)
                label = f"{source}_{skill_name}_{index}_{safe_label('_'.join(f'{k}_{v}' for k, v in request_query.items()))}"
                status, data, saved = fetch_search(base, label, request_query, portal=(source == "portal"))
                summary = summarize_rows(source, skill_name, request_query, status, data, saved)
                summaries.append(summary)
                manifests.extend(manifest_rows(source, skill_name, summary))
                print(
                    f"{source} {skill_name} query {index}: HTTP {status}, "
                    f"rows={summary['returned_rows']}, total={summary['total']}"
                )
    return summaries, manifests


def write_manifest(manifests: list[dict[str, str]], label: str) -> Path:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    path = MANIFEST_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_manifest.csv"
    fields = [
        "source",
        "skill",
        "query",
        "accession",
        "id",
        "status",
        "assay",
        "biosample",
        "file_format",
        "output_type",
        "assembly",
        "href",
        "summary",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in manifests:
            writer.writerow(row)
    logging.info("Wrote manifest: %s", path)
    return path


def write_report(summaries: list[dict[str, Any]], manifest_path: Path, label: str) -> Path:
    SINGLE_CELL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SINGLE_CELL_DOC_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_report.md"
    sections = [
        "# Single-Cell Dataset Search Report",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"Manifest: `{manifest_path}`",
        "",
    ]
    for summary in summaries:
        skill = SKILLS[summary["skill"]]
        sections.extend(
            [
                f"## {summary['source']} {skill['title']}",
                "",
                f"Query: `{summary['query']}`",
                f"HTTP status: {summary['http_status']}",
                f"Returned rows: {summary['returned_rows']}",
                f"Reported total: {summary['total']}",
                f"Saved response: `{summary['saved_response']}`",
                "",
                "Top assays:",
                format_counter(summary["assays"]),
                "",
                "Top biosamples:",
                format_counter(summary["biosamples"]),
                "",
                "Top file/output types:",
                format_counter(summary["file_formats"] + summary["output_types"]),
                "",
            ]
        )
    path.write_text("\n".join(sections), encoding="utf-8")
    logging.info("Wrote report: %s", path)
    return path


def maybe_float(value: str | None) -> float | None:
    if value is None:
        return None
    value = value.strip()
    if not value or value.upper() in {"NA", "NAN", "NULL", "NONE"}:
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    if math.isnan(number) or math.isinf(number):
        return None
    return number


def median(values: list[float]) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def quantile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    index = (len(ordered) - 1) * q
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[int(index)]
    return ordered[lower] * (upper - index) + ordered[upper] * (index - lower)


def summarize_numeric(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"n": 0, "min": None, "q25": None, "median": None, "q75": None, "max": None, "mean": None}
    return {
        "n": len(values),
        "min": min(values),
        "q25": quantile(values, 0.25),
        "median": median(values),
        "q75": quantile(values, 0.75),
        "max": max(values),
        "mean": sum(values) / len(values),
    }


def fmt_number(value: Any, digits: int = 3) -> str:
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if abs(value) >= 1000:
            return f"{value:,.0f}"
        return f"{value:.{digits}f}".rstrip("0").rstrip(".")
    return str(value)


def escape_xml(value: Any) -> str:
    text = str(value)
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


def palette(index: int) -> str:
    colors = [
        "#2f6f9f",
        "#c44e52",
        "#55a868",
        "#8172b2",
        "#ccb974",
        "#64b5cd",
        "#8c613c",
        "#dc7ec0",
        "#4c72b0",
        "#dd8452",
        "#937860",
        "#6acc64",
    ]
    return colors[index % len(colors)]


def scale(value: float, source_min: float, source_max: float, target_min: float, target_max: float) -> float:
    if source_max == source_min:
        return (target_min + target_max) / 2
    return target_min + (value - source_min) * (target_max - target_min) / (source_max - source_min)


def svg_frame(width: int, height: int, title: str, x_label: str, y_label: str, body: list[str]) -> str:
    return "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            "<style>",
            "text { font-family: Arial, Helvetica, sans-serif; fill: #202124; }",
            ".title { font-size: 18px; font-weight: 700; }",
            ".axis { stroke: #333; stroke-width: 1.2; }",
            ".grid { stroke: #d7dce2; stroke-width: 0.8; }",
            ".tick { font-size: 11px; fill: #4f5b66; }",
            ".label { font-size: 13px; font-weight: 600; }",
            ".legend { font-size: 11px; }",
            "</style>",
            f'<rect width="{width}" height="{height}" fill="#ffffff"/>',
            f'<text class="title" x="28" y="30">{escape_xml(title)}</text>',
            *body,
            f'<text class="label" x="{width / 2}" y="{height - 14}" text-anchor="middle">{escape_xml(x_label)}</text>',
            f'<text class="label" transform="translate(15 {height / 2}) rotate(-90)" text-anchor="middle">{escape_xml(y_label)}</text>',
            "</svg>",
        ]
    )


def write_histogram(values: list[float], path: Path, title: str, x_label: str, bins: int = 25) -> Path:
    SINGLE_CELL_PLOT_DIR.mkdir(parents=True, exist_ok=True)
    width, height = 840, 520
    left, right, top, bottom = 78, 26, 52, 72
    plot_width = width - left - right
    plot_height = height - top - bottom
    if not values:
        path.write_text(svg_frame(width, height, title, x_label, "Count", []), encoding="utf-8")
        return path
    low, high = min(values), max(values)
    if low == high:
        low -= 0.5
        high += 0.5
    counts = [0 for _ in range(bins)]
    for value in values:
        index = int((value - low) / (high - low) * bins)
        counts[min(index, bins - 1)] += 1
    max_count = max(counts) or 1
    body = [
        f'<line class="axis" x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}"/>',
    ]
    for tick in range(5):
        y = top + plot_height - tick * plot_height / 4
        count = max_count * tick / 4
        body.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{left + plot_width}" y2="{y:.1f}"/>')
        body.append(f'<text class="tick" x="{left - 9}" y="{y + 4:.1f}" text-anchor="end">{fmt_number(count, 1)}</text>')
    bar_width = plot_width / bins
    for index, count in enumerate(counts):
        x = left + index * bar_width
        y = top + plot_height - (count / max_count) * plot_height
        body.append(
            f'<rect x="{x + 1:.1f}" y="{y:.1f}" width="{max(bar_width - 2, 1):.1f}" '
            f'height="{top + plot_height - y:.1f}" fill="#3b7ea1"/>'
        )
    for tick in range(6):
        value = low + (high - low) * tick / 5
        x = left + plot_width * tick / 5
        body.append(f'<line class="axis" x1="{x:.1f}" y1="{top + plot_height}" x2="{x:.1f}" y2="{top + plot_height + 5}"/>')
        body.append(f'<text class="tick" x="{x:.1f}" y="{top + plot_height + 20}" text-anchor="middle">{fmt_number(value, 2)}</text>')
    path.write_text(svg_frame(width, height, title, x_label, "Count", body), encoding="utf-8")
    return path


def write_bar_plot(items: list[tuple[str, int]], path: Path, title: str, x_label: str, y_label: str, limit: int = 20) -> Path:
    SINGLE_CELL_PLOT_DIR.mkdir(parents=True, exist_ok=True)
    width, height = 960, 560
    left, right, top, bottom = 82, 28, 54, 132
    plot_width = width - left - right
    plot_height = height - top - bottom
    data = items[:limit]
    max_count = max([count for _, count in data], default=1)
    body = [
        f'<line class="axis" x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}"/>',
    ]
    for tick in range(5):
        y = top + plot_height - tick * plot_height / 4
        count = max_count * tick / 4
        body.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{left + plot_width}" y2="{y:.1f}"/>')
        body.append(f'<text class="tick" x="{left - 9}" y="{y + 4:.1f}" text-anchor="end">{fmt_number(count, 1)}</text>')
    if data:
        slot = plot_width / len(data)
        bar_width = min(slot * 0.74, 42)
        for index, (label, count) in enumerate(data):
            x = left + index * slot + (slot - bar_width) / 2
            bar_height = count / max_count * plot_height
            y = top + plot_height - bar_height
            body.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{bar_height:.1f}" '
                f'fill="{palette(index)}"/>'
            )
            body.append(f'<text class="tick" x="{x + bar_width / 2:.1f}" y="{y - 5:.1f}" text-anchor="middle">{count}</text>')
            body.append(
                f'<text class="tick" transform="translate({x + bar_width / 2:.1f} {top + plot_height + 12}) rotate(55)" '
                f'text-anchor="start">{escape_xml(label[:28])}</text>'
            )
    path.write_text(svg_frame(width, height, title, x_label, y_label, body), encoding="utf-8")
    return path


def write_scatter(
    points: list[dict[str, Any]],
    path: Path,
    title: str,
    x_label: str,
    y_label: str,
    *,
    label_key: str | None = None,
    color_key: str | None = None,
    max_points: int = 7000,
) -> Path:
    SINGLE_CELL_PLOT_DIR.mkdir(parents=True, exist_ok=True)
    width, height = 900, 620
    left, right, top, bottom = 82, 180 if color_key else 32, 56, 74
    plot_width = width - left - right
    plot_height = height - top - bottom
    usable = [p for p in points if p.get("x") is not None and p.get("y") is not None]
    if len(usable) > max_points:
        random.seed(17)
        usable = random.sample(usable, max_points)
    xs = [float(p["x"]) for p in usable] or [0.0]
    ys = [float(p["y"]) for p in usable] or [0.0]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    if x_min == x_max:
        x_min -= 0.5
        x_max += 0.5
    if y_min == y_max:
        y_min -= 0.5
        y_max += 0.5
    x_pad = (x_max - x_min) * 0.04
    y_pad = (y_max - y_min) * 0.04
    x_min -= x_pad
    x_max += x_pad
    y_min -= y_pad
    y_max += y_pad
    categories: dict[str, str] = {}
    if color_key:
        for point in usable:
            category = str(point.get(color_key) or "unknown")
            if category not in categories:
                categories[category] = palette(len(categories))
    body = [
        f'<line class="axis" x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}"/>',
    ]
    for tick in range(6):
        x = left + tick * plot_width / 5
        value = x_min + tick * (x_max - x_min) / 5
        body.append(f'<line class="grid" x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + plot_height}"/>')
        body.append(f'<text class="tick" x="{x:.1f}" y="{top + plot_height + 19}" text-anchor="middle">{fmt_number(value, 2)}</text>')
        y = top + plot_height - tick * plot_height / 5
        y_value = y_min + tick * (y_max - y_min) / 5
        body.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{left + plot_width}" y2="{y:.1f}"/>')
        body.append(f'<text class="tick" x="{left - 9}" y="{y + 4:.1f}" text-anchor="end">{fmt_number(y_value, 2)}</text>')
    for point in usable:
        x = scale(float(point["x"]), x_min, x_max, left, left + plot_width)
        y = scale(float(point["y"]), y_min, y_max, top + plot_height, top)
        color = categories.get(str(point.get(color_key)), "#3b7ea1") if color_key else "#3b7ea1"
        opacity = 0.55 if len(usable) > 1000 else 0.78
        body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.8" fill="{color}" fill-opacity="{opacity}"/>')
    if label_key:
        for point in sorted(usable, key=lambda item: float(item.get("label_score") or item.get("y") or 0), reverse=True)[:14]:
            x = scale(float(point["x"]), x_min, x_max, left, left + plot_width)
            y = scale(float(point["y"]), y_min, y_max, top + plot_height, top)
            label = escape_xml(point.get(label_key, ""))
            if label:
                body.append(f'<text class="tick" x="{x + 5:.1f}" y="{y - 5:.1f}">{label}</text>')
    if color_key and categories:
        legend_x = left + plot_width + 20
        body.append(f'<text class="label" x="{legend_x}" y="{top + 5}">{escape_xml(color_key)}</text>')
        for index, (category, color) in enumerate(list(categories.items())[:18]):
            y = top + 25 + index * 21
            body.append(f'<circle cx="{legend_x + 5}" cy="{y - 4}" r="5" fill="{color}"/>')
            body.append(f'<text class="legend" x="{legend_x + 17}" y="{y}">{escape_xml(category[:24])}</text>')
    path.write_text(svg_frame(width, height, title, x_label, y_label, body), encoding="utf-8")
    return path


def download_file(url: str, path: Path, force: bool = False) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not force:
        return {"path": str(path), "url": url, "status": "exists", "bytes": path.stat().st_size}
    headers = {"User-Agent": "IGVFdataAgent/0.1", "Accept": "*/*"}
    request = urllib.request.Request(url, headers=headers, method="GET")
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            path.write_bytes(response.read())
        return {"path": str(path), "url": url, "status": "downloaded", "bytes": path.stat().st_size}
    except urllib.error.HTTPError as exc:
        return {"path": str(path), "url": url, "status": f"http_error_{exc.code}", "bytes": 0}
    except urllib.error.URLError as exc:
        return {"path": str(path), "url": url, "status": f"network_error_{exc.reason}", "bytes": 0}


def download_examples(force: bool = False) -> tuple[list[dict[str, Any]], Path]:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for item in EXAMPLE_FILES:
        result = download_file(item["url"], Path(item["path"]), force=force)
        result.update(
            {
                "accession": item["accession"],
                "modality": item["modality"],
                "description": item["description"],
            }
        )
        rows.append(result)
    path = MANIFEST_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_single_cell_example_downloads.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["accession", "modality", "description", "url", "path", "status", "bytes"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    logging.info("Wrote example download manifest: %s", path)
    return rows, path


def analyze_scrna_example(path: Path) -> dict[str, Any]:
    genes: list[dict[str, Any]] = []
    metadata: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header: list[str] | None = None
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                metadata.append(line[1:].strip())
                continue
            if header is None:
                header = line.split("\t")
                continue
            parts = line.split("\t")
            row = dict(zip(header, parts))
            record = {
                "gene": row.get("GeneSymbol", ""),
                "chrom": row.get("TSSChr", ""),
                "promoter_atac": maybe_float(row.get("normalizedATAC_prom")),
                "rna_log": maybe_float(row.get("RNA_meanLogNorm")),
                "tpm": maybe_float(row.get("RNA_pseudobulkTPM")),
                "percent_detected": maybe_float(row.get("RNA_percentCellsDetected")),
                "ubiq_expressed": str(row.get("ubiqExpressed", "")).upper() == "TRUE",
            }
            genes.append(record)
    expressed = [gene for gene in genes if gene["rna_log"] is not None]
    promoter_values = [float(gene["promoter_atac"]) for gene in genes if gene["promoter_atac"] is not None]
    rna_values = [float(gene["rna_log"]) for gene in expressed]
    tpm_values = [float(gene["tpm"]) for gene in genes if gene["tpm"] is not None]
    pct_values = [float(gene["percent_detected"]) for gene in genes if gene["percent_detected"] is not None]
    scatter_points = [
        {
            "x": gene["promoter_atac"],
            "y": gene["rna_log"],
            "gene": gene["gene"],
            "label_score": (gene["rna_log"] or 0) + math.log10((gene["promoter_atac"] or 0) + 1),
        }
        for gene in genes
        if gene["promoter_atac"] is not None and gene["rna_log"] is not None
    ]
    scatter_path = SINGLE_CELL_PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_scrna_promoter_atac_vs_rna.svg"
    write_scatter(
        scatter_points,
        scatter_path,
        "K562 CRISPRi 10x Multiome Gene Expression vs Promoter Accessibility",
        "normalized promoter ATAC",
        "RNA mean log-normalized expression",
        label_key="gene",
    )
    pct_hist_path = SINGLE_CELL_PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_scrna_percent_cells_detected.svg"
    write_histogram(pct_values, pct_hist_path, "RNA Detection Fraction Across Genes", "percent cells detected")
    paired_promoter_values = [float(gene["promoter_atac"]) for gene in expressed if gene["promoter_atac"] is not None]
    if paired_promoter_values:
        cutoff = quantile(paired_promoter_values, 0.75) or 0
        low_cutoff = quantile(paired_promoter_values, 0.25) or 0
    else:
        cutoff = low_cutoff = 0
    high_promoter_rna = [float(g["rna_log"]) for g in expressed if g["promoter_atac"] is not None and g["promoter_atac"] >= cutoff]
    low_promoter_rna = [float(g["rna_log"]) for g in expressed if g["promoter_atac"] is not None and g["promoter_atac"] <= low_cutoff]
    top_by_expression = sorted(expressed, key=lambda gene: gene["rna_log"] or -1, reverse=True)[:15]
    top_by_detection = sorted(
        [gene for gene in genes if gene["percent_detected"] is not None],
        key=lambda gene: gene["percent_detected"] or -1,
        reverse=True,
    )[:15]
    return {
        "input": str(path),
        "metadata": metadata[:12],
        "genes_total": len(genes),
        "genes_with_rna": len(expressed),
        "ubiquitously_expressed": sum(1 for gene in genes if gene["ubiq_expressed"]),
        "rna_mean_log_norm": summarize_numeric(rna_values),
        "rna_pseudobulk_tpm": summarize_numeric(tpm_values),
        "rna_percent_cells_detected": summarize_numeric(pct_values),
        "promoter_atac": summarize_numeric(promoter_values),
        "high_vs_low_promoter_accessibility": {
            "high_promoter_cutoff_q75": cutoff,
            "low_promoter_cutoff_q25": low_cutoff,
            "high_promoter_gene_count": len(high_promoter_rna),
            "low_promoter_gene_count": len(low_promoter_rna),
            "median_rna_high_promoter": median(high_promoter_rna),
            "median_rna_low_promoter": median(low_promoter_rna),
        },
        "top_genes_by_expression": top_by_expression,
        "top_genes_by_detection": top_by_detection,
        "plots": [str(scatter_path), str(pct_hist_path)],
    }


def analyze_scatac_example(path: Path, max_rows: int | None = None) -> dict[str, Any]:
    chrom_counts: collections.Counter[str] = collections.Counter()
    lengths: list[float] = []
    signal_values: list[float] = []
    points: list[dict[str, Any]] = []
    top_peaks: list[dict[str, Any]] = []
    rows = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            start = maybe_float(parts[1])
            end = maybe_float(parts[2])
            if start is None or end is None:
                continue
            rows += 1
            chrom = parts[0]
            length = max(0.0, end - start)
            signal = maybe_float(parts[6] if len(parts) > 6 else None)
            score = maybe_float(parts[7] if len(parts) > 7 else None)
            chrom_counts[chrom] += 1
            lengths.append(length)
            if signal is not None:
                signal_values.append(signal)
                points.append({"x": length, "y": signal, "peak": parts[3] if len(parts) > 3 else chrom, "label_score": signal})
                top_peaks.append(
                    {
                        "chrom": chrom,
                        "start": int(start),
                        "end": int(end),
                        "peak": parts[3] if len(parts) > 3 else "",
                        "signal": signal,
                        "score": score,
                    }
                )
            if max_rows and rows >= max_rows:
                break
    top_peaks = sorted(top_peaks, key=lambda peak: peak.get("signal") or -1, reverse=True)[:15]
    length_hist_path = SINGLE_CELL_PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_scatac_peak_length_histogram.svg"
    write_histogram(lengths, length_hist_path, "Single-Cell ATAC Peak Length Distribution", "peak length (bp)")
    chrom_bar_path = SINGLE_CELL_PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_scatac_peaks_by_chromosome.svg"
    write_bar_plot(chrom_counts.most_common(), chrom_bar_path, "Single-Cell ATAC Peaks By Chromosome", "chromosome", "peak count", limit=25)
    signal_scatter_path = SINGLE_CELL_PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_scatac_peak_length_vs_signal.svg"
    write_scatter(points, signal_scatter_path, "Single-Cell ATAC Peak Length vs Signal", "peak length (bp)", "signal value", label_key="peak")
    return {
        "input": str(path),
        "peaks_total": rows,
        "peak_length_bp": summarize_numeric(lengths),
        "signal_value": summarize_numeric(signal_values),
        "chromosomes": chrom_counts.most_common(25),
        "top_peaks_by_signal": top_peaks,
        "plots": [str(length_hist_path), str(chrom_bar_path), str(signal_scatter_path)],
    }


def find_multiome_cell_annotation_files() -> list[Path]:
    base = DATA_DIR / "IGVF" / "10xMultiome" / "Downloads"
    if not base.exists():
        return []
    candidates: list[Path] = []
    for path in sorted(base.glob("*/*/*.tsv.gz")):
        try:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
                header = handle.readline().rstrip("\n").split("\t")
            if {"cell_barcode", "CL_term_name", "BrainRegion"}.issubset(set(header)):
                candidates.append(path)
        except OSError:
            continue
    return candidates


def analyze_multiome_cell_annotations(max_cells: int = 12000) -> dict[str, Any]:
    files = find_multiome_cell_annotation_files()
    cell_type_counts: collections.Counter[str] = collections.Counter()
    region_counts: collections.Counter[str] = collections.Counter()
    joint_counts: collections.Counter[tuple[str, str]] = collections.Counter()
    modality_counts: collections.Counter[str] = collections.Counter()
    points: list[dict[str, Any]] = []
    rows = 0
    top_cell_type_order: list[str] = []
    random.seed(42)
    for path in files:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                cell_type = row.get("CL_term_name") or row.get("cell_name") or "unknown"
                region = row.get("BrainRegion") or "unknown"
                in_gex = str(row.get("in_GEX", "")).upper() == "TRUE"
                in_atac = str(row.get("in_ATAC", "")).upper() == "TRUE"
                if in_gex and in_atac:
                    modality = "RNA+ATAC"
                elif in_gex:
                    modality = "RNA only"
                elif in_atac:
                    modality = "ATAC only"
                else:
                    modality = "metadata only"
                cell_type_counts[cell_type] += 1
                region_counts[region] += 1
                joint_counts[(cell_type, region)] += 1
                modality_counts[modality] += 1
                rows += 1
                if rows <= max_cells:
                    if cell_type not in top_cell_type_order:
                        top_cell_type_order.append(cell_type)
                    cell_index = top_cell_type_order.index(cell_type)
                    region_index = list(region_counts.keys()).index(region)
                    angle = 2 * math.pi * (cell_index % 12) / 12
                    radius = 5 + (cell_index // 12) * 1.2
                    x = math.cos(angle) * radius + (region_index % 5 - 2) * 0.35 + random.gauss(0, 0.33)
                    y = math.sin(angle) * radius + (region_index // 5) * 0.35 + random.gauss(0, 0.33)
                    points.append({"x": x, "y": y, "cell_type": cell_type, "region": region, "modality": modality})
    top_cell_types = cell_type_counts.most_common(20)
    top_regions = region_counts.most_common(20)
    cell_plot_path = SINGLE_CELL_PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_multiome_annotation_cell_embedding.svg"
    write_scatter(
        points,
        cell_plot_path,
        "10x Multiome Cell Landscape From IGVF Cell Annotations",
        "embedding 1 (annotation-derived demo)",
        "embedding 2 (annotation-derived demo)",
        color_key="cell_type",
        max_points=max_cells,
    )
    cell_type_bar_path = SINGLE_CELL_PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_multiome_cell_type_counts.svg"
    write_bar_plot(top_cell_types, cell_type_bar_path, "10x Multiome Cell Type Counts", "cell type", "cell count", limit=18)
    region_bar_path = SINGLE_CELL_PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_multiome_brain_region_counts.svg"
    write_bar_plot(top_regions, region_bar_path, "10x Multiome Brain Region Counts", "brain region", "cell count", limit=18)
    return {
        "input_files": [str(path) for path in files],
        "cells_total": rows,
        "cells_plotted": len(points),
        "cell_types": top_cell_types,
        "brain_regions": top_regions,
        "modality_flags": modality_counts.most_common(),
        "top_cell_type_region_pairs": [(f"{cell_type} / {region}", count) for (cell_type, region), count in joint_counts.most_common(20)],
        "plots": [str(cell_plot_path), str(cell_type_bar_path), str(region_bar_path)],
        "note": "Coordinates are deterministic annotation-derived demo coordinates for visualization. True tSNE/UMAP requires a cell-by-gene or cell-by-peak matrix plus Scanpy/ArchR/Signac.",
    }


def read_tar_text_lines(archive_path: Path, suffix: str) -> list[str]:
    with tarfile.open(archive_path, "r:gz") as archive:
        member = next((item for item in archive.getmembers() if item.name.endswith(suffix)), None)
        if member is None:
            return []
        extracted = archive.extractfile(member)
        if extracted is None:
            return []
        raw = extracted.read()
    if suffix.endswith(".gz"):
        return gzip.decompress(raw).decode("utf-8", errors="replace").splitlines()
    return raw.decode("utf-8", errors="replace").splitlines()


def cell_annotations_by_barcode(directory: Path) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    for path in sorted(directory.glob("*.tsv.gz")):
        try:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if not reader.fieldnames or "cell_barcode" not in reader.fieldnames:
                    continue
                for row in reader:
                    barcode = row.get("cell_barcode", "")
                    if "|" in barcode:
                        barcode = barcode.split("|", 1)[1]
                    if barcode:
                        mapping[barcode] = row
        except OSError:
            continue
    return mapping


def analyze_expression_matrix_embedding(max_cells: int = 3500) -> dict[str, Any]:
    base = DATA_DIR / "IGVF" / "10xMultiome" / "Downloads"
    archives = sorted(base.glob("*/*/*.tar.gz"))
    if not archives:
        return {"available": False, "reason": "No local 10x multiome Matrix Market tarballs found.", "plots": []}
    archive_path = archives[0]
    directory = archive_path.parent
    barcodes = read_tar_text_lines(archive_path, "_barcodes.tsv.gz")
    features = read_tar_text_lines(archive_path, "_features.tsv.gz")
    if not barcodes or not features:
        return {"available": False, "reason": f"Missing barcodes/features in {archive_path}", "plots": []}
    feature_to_index: dict[str, int] = {}
    for index, line in enumerate(features, start=1):
        parts = line.split("\t")
        names = [part for part in parts[:2] if part]
        for name in names:
            feature_to_index.setdefault(name.upper(), index)
    marker_modules = {
        "neuron": ["SNAP25", "RBFOX3", "SYT1", "TUBB3"],
        "excitatory": ["SLC17A7", "SLC17A6", "CAMK2A"],
        "inhibitory": ["GAD1", "GAD2", "PVALB", "SST", "VIP"],
        "oligodendrocyte": ["MBP", "MOG", "PLP1", "MAG"],
        "astrocyte": ["GFAP", "AQP4", "SLC1A3", "ALDH1L1"],
        "microglia": ["C1QA", "CSF1R", "CX3CR1", "TYROBP"],
        "opc": ["PDGFRA", "CSPG4", "VCAN"],
    }
    index_to_modules: dict[int, list[str]] = collections.defaultdict(list)
    found_markers: dict[str, list[str]] = {}
    for module, genes in marker_modules.items():
        for gene in genes:
            feature_index = feature_to_index.get(gene.upper())
            if feature_index:
                index_to_modules[feature_index].append(module)
                found_markers.setdefault(module, []).append(gene)
    if not index_to_modules:
        return {"available": False, "reason": f"No marker genes found in {archive_path}", "plots": []}
    selected_cells = set(range(1, min(len(barcodes), max_cells) + 1))
    module_scores: dict[int, collections.Counter[str]] = {
        cell_index: collections.Counter() for cell_index in selected_cells
    }
    library_sizes: collections.Counter[int] = collections.Counter()
    with tarfile.open(archive_path, "r:gz") as archive:
        member = next((item for item in archive.getmembers() if item.name.endswith("_counts.mtx")), None)
        if member is None:
            return {"available": False, "reason": f"Missing counts.mtx in {archive_path}", "plots": []}
        extracted = archive.extractfile(member)
        if extracted is None:
            return {"available": False, "reason": f"Could not read counts.mtx in {archive_path}", "plots": []}
        header_seen = False
        for raw in extracted:
            line = raw.decode("utf-8", errors="replace").strip()
            if not line or line.startswith("%"):
                continue
            if not header_seen:
                header_seen = True
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            feature_index = int(parts[0])
            cell_index = int(parts[1])
            value = float(parts[2])
            if cell_index not in selected_cells:
                continue
            library_sizes[cell_index] += int(value)
            modules = index_to_modules.get(feature_index)
            if modules:
                for module in modules:
                    module_scores[cell_index][module] += value
    annotations = cell_annotations_by_barcode(directory)
    points: list[dict[str, Any]] = []
    umap_points: list[dict[str, Any]] = []
    tsne_points: list[dict[str, Any]] = []
    cluster_counts: collections.Counter[str] = collections.Counter()
    module_names = list(marker_modules)
    random.seed(71)
    for cell_index in sorted(selected_cells):
        barcode = barcodes[cell_index - 1]
        annotation = annotations.get(barcode, {})
        library = max(library_sizes[cell_index], 1)
        scores = {
            module: math.log1p((module_scores[cell_index][module] / library) * 10000)
            for module in marker_modules
        }
        neuronal_axis = scores["neuron"] + scores["excitatory"] + scores["inhibitory"] - scores["oligodendrocyte"] - scores["astrocyte"] - scores["microglia"]
        glia_axis = scores["oligodendrocyte"] + scores["astrocyte"] + scores["microglia"] + scores["opc"] - scores["neuron"]
        dominant_module = max(scores, key=lambda module: scores[module])
        cell_type = annotation.get("CL_term_name") or annotation.get("cell_name") or "unlabeled"
        cluster = cell_type if cell_type != "unlabeled" else dominant_module
        cluster_counts[cluster] += 1
        module_index = module_names.index(dominant_module)
        module_angle = 2 * math.pi * module_index / len(module_names)
        module_strength = scores[dominant_module]
        umap_x = neuronal_axis + math.cos(module_angle) * 0.9 + random.gauss(0, 0.16)
        umap_y = glia_axis + math.sin(module_angle) * 0.9 + random.gauss(0, 0.16)
        tsne_radius = 4.6 + min(module_strength, 3.5) * 0.32
        tsne_x = math.cos(module_angle) * tsne_radius + random.gauss(0, 0.42)
        tsne_y = math.sin(module_angle) * tsne_radius + random.gauss(0, 0.42)
        point_common = {
            "cell_type": cell_type,
            "cluster": cluster,
            "dominant_marker_module": dominant_module,
            "barcode": barcode,
            "label_score": abs(neuronal_axis) + abs(glia_axis),
        }
        points.append(
            {
                "x": neuronal_axis,
                "y": glia_axis,
                **point_common,
            }
        )
        umap_points.append(
            {
                "x": umap_x,
                "y": umap_y,
                **point_common,
            }
        )
        tsne_points.append(
            {
                "x": tsne_x,
                "y": tsne_y,
                **point_common,
            }
        )
    coordinates_path = SINGLE_CELL_DOC_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_scrna_umap_tsne_demo_coordinates.csv"
    with coordinates_path.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "cell_barcode",
            "cell_type",
            "cluster",
            "dominant_marker_module",
            "umap_1",
            "umap_2",
            "tsne_1",
            "tsne_2",
            "marker_axis_1",
            "marker_axis_2",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for base_point, umap_point, tsne_point in zip(points, umap_points, tsne_points):
            writer.writerow(
                {
                    "cell_barcode": base_point["barcode"],
                    "cell_type": base_point["cell_type"],
                    "cluster": base_point["cluster"],
                    "dominant_marker_module": base_point["dominant_marker_module"],
                    "umap_1": f"{umap_point['x']:.5f}",
                    "umap_2": f"{umap_point['y']:.5f}",
                    "tsne_1": f"{tsne_point['x']:.5f}",
                    "tsne_2": f"{tsne_point['y']:.5f}",
                    "marker_axis_1": f"{base_point['x']:.5f}",
                    "marker_axis_2": f"{base_point['y']:.5f}",
                }
            )
    plot_path = SINGLE_CELL_PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_multiome_expression_marker_embedding.svg"
    write_scatter(
        points,
        plot_path,
        "10x Multiome Expression Marker-Score Embedding",
        "neuronal marker score minus glial marker score",
        "glial marker score minus neuronal marker score",
        color_key="cell_type",
        max_points=max_cells,
    )
    umap_plot_path = SINGLE_CELL_PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_scrna_umap_style_clusters.svg"
    write_scatter(
        umap_points,
        umap_plot_path,
        "scRNA-seq UMAP-Style Clusters Colored By Cell Type",
        "UMAP 1 (marker-score demo)",
        "UMAP 2 (marker-score demo)",
        color_key="cluster",
        max_points=max_cells,
    )
    tsne_plot_path = SINGLE_CELL_PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_scrna_tsne_style_clusters.svg"
    write_scatter(
        tsne_points,
        tsne_plot_path,
        "scRNA-seq tSNE-Style Clusters Colored By Cell Type",
        "tSNE 1 (marker-score demo)",
        "tSNE 2 (marker-score demo)",
        color_key="cluster",
        max_points=max_cells,
    )
    return {
        "available": True,
        "input_archive": str(archive_path),
        "cells_plotted": len(points),
        "barcodes_total": len(barcodes),
        "features_total": len(features),
        "found_markers": found_markers,
        "clusters": cluster_counts.most_common(20),
        "coordinates": str(coordinates_path),
        "plots": [str(plot_path), str(umap_plot_path), str(tsne_plot_path)],
        "note": "These are lightweight expression-matrix marker-score embeddings from Matrix Market counts. The UMAP/tSNE-style plots are demo coordinates, not Scanpy/Seurat UMAP or tSNE, but they are useful paper-style cluster illustrations before the heavier analysis stack is installed.",
    }


def markdown_table(rows: list[dict[str, Any]], columns: list[tuple[str, str]], limit: int = 10) -> list[str]:
    lines = ["| " + " | ".join(header for header, _ in columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows[:limit]:
        lines.append("| " + " | ".join(escape_xml(fmt_number(row.get(key))) for _, key in columns) + " |")
    return lines


def write_example_analysis_report(summary: dict[str, Any], label: str) -> Path:
    SINGLE_CELL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SINGLE_CELL_DOC_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_single_cell_example_analysis_report.md"
    scrna = summary["scrna"]
    scatac = summary["scatac"]
    multiome = summary["multiome_annotations"]
    expression = summary.get("expression_embedding", {})
    lines = [
        "# Single-Cell RNA/ATAC Example Analysis",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        "This report uses small public IGVF examples for smoke analysis and the already-downloaded IGVF 10x multiome cell annotation files for a cell-landscape demo. The plots are intentionally standard-library SVGs so the skill can run on a clean workstation.",
        "",
        "## Inputs",
        "",
        f"- scRNA/scE2G example: `{scrna['input']}`",
        f"- scATAC peaks example: `{scatac['input']}`",
        f"- 10x multiome cell annotation files: {len(multiome['input_files'])}",
        "",
        "## scRNA/scE2G Summary",
        "",
        f"- Genes total: {scrna['genes_total']:,}",
        f"- Genes with RNA values: {scrna['genes_with_rna']:,}",
        f"- Ubiquitously expressed genes: {scrna['ubiquitously_expressed']:,}",
        f"- Median RNA mean log-normalized expression: {fmt_number(scrna['rna_mean_log_norm']['median'])}",
        f"- Median RNA percent cells detected: {fmt_number(scrna['rna_percent_cells_detected']['median'])}",
        f"- Median RNA in high promoter-accessibility genes: {fmt_number(scrna['high_vs_low_promoter_accessibility']['median_rna_high_promoter'])}",
        f"- Median RNA in low promoter-accessibility genes: {fmt_number(scrna['high_vs_low_promoter_accessibility']['median_rna_low_promoter'])}",
        "",
        "Top genes by expression:",
        *markdown_table(
            scrna["top_genes_by_expression"],
            [("Gene", "gene"), ("RNA log", "rna_log"), ("TPM", "tpm"), ("% cells", "percent_detected"), ("Promoter ATAC", "promoter_atac")],
            limit=10,
        ),
        "",
        "Plots:",
        *[f"- `{plot}`" for plot in scrna["plots"]],
        "",
        "## scATAC Peak Summary",
        "",
        f"- Peaks parsed: {scatac['peaks_total']:,}",
        f"- Median peak length: {fmt_number(scatac['peak_length_bp']['median'])} bp",
        f"- Median signal value: {fmt_number(scatac['signal_value']['median'])}",
        "",
        "Top peaks by signal:",
        *markdown_table(
            scatac["top_peaks_by_signal"],
            [("Peak", "peak"), ("Chrom", "chrom"), ("Start", "start"), ("End", "end"), ("Signal", "signal"), ("Score", "score")],
            limit=10,
        ),
        "",
        "Plots:",
        *[f"- `{plot}`" for plot in scatac["plots"]],
        "",
        "## 10x Multiome Cell Landscape",
        "",
        f"- Cells summarized: {multiome['cells_total']:,}",
        f"- Cells plotted: {multiome['cells_plotted']:,}",
        f"- Modality flags: {', '.join(f'{name}: {count:,}' for name, count in multiome['modality_flags'])}",
        "",
        "Top cell types:",
        *[f"- {name}: {count:,}" for name, count in multiome["cell_types"][:12]],
        "",
        "Top brain regions:",
        *[f"- {name}: {count:,}" for name, count in multiome["brain_regions"][:12]],
        "",
        "Plots:",
        *[f"- `{plot}`" for plot in multiome["plots"]],
        "",
        "Important interpretation note: the multiome embedding is an annotation-driven demo of cell-type separation, not a matrix-derived tSNE/UMAP. For a publication-grade figure, use the optional Scanpy/ArchR/Signac workflow in the skill document on h5ad, Cell Ranger, fragments, and peak matrices.",
        "",
        "## Expression-Matrix Cell Embedding",
        "",
    ]
    if expression.get("available"):
        lines.extend(
            [
                f"- Matrix archive: `{expression['input_archive']}`",
                f"- Cells plotted: {expression['cells_plotted']:,}",
                f"- Features in matrix: {expression['features_total']:,}",
                f"- Marker modules found: {', '.join(f'{module} ({len(genes)})' for module, genes in expression['found_markers'].items())}",
                f"- UMAP/tSNE-style coordinate table: `{expression['coordinates']}`",
                "",
                "Top plotted clusters:",
                *[f"- {name}: {count:,}" for name, count in expression.get("clusters", [])[:12]],
                "",
                "Plots:",
                *[f"- `{plot}`" for plot in expression["plots"]],
                "",
                expression["note"],
                "",
            ]
        )
    else:
        lines.extend([f"- Not available: {expression.get('reason', 'unknown reason')}", ""])
    lines.extend(
        [
        "## Full Matrix Workflow",
        "",
        "For a true tSNE/UMAP differential cell-state plot, install Scanpy/AnnData or Seurat/Signac and run the full matrix workflow from the skill document on h5ad or Cell Ranger outputs. The lightweight plots here are meant to make IGVF files understandable before the heavier analysis stack is available.",
        "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")
    logging.info("Wrote example analysis report: %s", path)
    return path


def analyze_examples(max_cells: int = 12000) -> tuple[dict[str, Any], Path, Path]:
    SINGLE_CELL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    SINGLE_CELL_PLOT_DIR.mkdir(parents=True, exist_ok=True)
    scrna = analyze_scrna_example(EXAMPLE_DIR / "IGVFFI9905RPTO.tsv.gz")
    scatac = analyze_scatac_example(EXAMPLE_DIR / "IGVFFI6024DNYI.bed.gz")
    multiome = analyze_multiome_cell_annotations(max_cells=max_cells)
    expression = analyze_expression_matrix_embedding(max_cells=min(max_cells, 3500))
    summary = {"scrna": scrna, "scatac": scatac, "multiome_annotations": multiome, "expression_embedding": expression}
    json_path = SINGLE_CELL_DOC_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_single_cell_example_analysis_summary.json"
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    report_path = write_example_analysis_report(summary, "igvf_examples")
    return summary, report_path, json_path


def write_playbook() -> Path:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "SINGLE_CELL_ANALYSIS_SKILLS.md"
    lines = [
        "# Skill: Single-Cell RNA-seq, Single-Cell ATAC-seq, And Perturb-seq Analysis",
        "",
        "Use this skill when the agent needs to find, summarize, download-plan, or analyze single-cell datasets from IGVF Portal or ENCODE.",
        "",
        "## Metadata First",
        "",
        "Always build a manifest before downloading data. Preserve source database, accession, assay, biosample, file format, output type, assembly, status, and href/download URL.",
        "",
        "Run:",
        "",
        "```bash",
        "python3 Scripts/single_cell_data_skills.py smoke --skill all --source encode --limit 5",
        "python3 Scripts/single_cell_data_skills.py manifest --skill scrna --source both --limit 25",
        "python3 Scripts/single_cell_data_skills.py download-examples",
        "python3 Scripts/single_cell_data_skills.py analyze-examples --max-cells 12000",
        "```",
        "",
        "## Nature/Cell-Style Analysis Pattern",
        "",
        "Recent high-impact single-cell RNA/ATAC papers usually follow the same spine: careful cell/sample QC, latent-space construction, cell-state annotation, differential testing, and biological interpretation through regulatory elements, genes, and perturbations. The command above runs a lightweight smoke version with public IGVF examples; full publication workflows should use Scanpy/AnnData, Seurat/Signac, ArchR, or SnapATAC on full matrices.",
        "",
        "The built-in example analysis includes a real Matrix Market expression-count pass when local 10x multiome tarballs are available. It scores marker modules and plots cells in a matrix-derived marker-score embedding; this is useful for quick interpretation before installing the heavier tSNE/UMAP stack.",
        "",
        "For scRNA-seq:",
        "",
        "1. Load raw counts into AnnData/Seurat and preserve raw counts.",
        "2. QC cells by total UMIs, detected genes, mitochondrial/ribosomal fraction, doublet score, donor, batch, and chemistry.",
        "3. Normalize, log-transform, select highly variable genes, regress only justified covariates, run PCA, neighbors, UMAP/tSNE, clustering, marker discovery, and differential expression.",
        "4. Annotate clusters with marker genes and reference mapping; report pseudobulk differential expression for condition contrasts when donors/replicates exist.",
        "",
        "For scATAC-seq:",
        "",
        "1. Start from fragments, peaks, cell metadata, and genome assembly.",
        "2. QC by fragments per cell, TSS enrichment, FRiP, nucleosome signal, blacklist fraction, and doublets.",
        "3. Build a cell-by-peak matrix, run TF-IDF/LSI, neighbors, UMAP/tSNE, clustering, differential accessibility, motif enrichment, and gene-activity estimates.",
        "4. Link peaks to genes with co-accessibility, multiome paired RNA, eQTL/rE2G evidence, enhancer perturbation evidence, and IGVF Catalog links.",
        "",
        "For multiome integration:",
        "",
        "1. Use paired cells to jointly inspect RNA expression and ATAC accessibility by cell type.",
        "2. Integrate with weighted nearest neighbor, MultiVI/scVI-style latent models, ArchR gene scores, or Signac bridge integration.",
        "3. Interpret variant loci by asking whether the variant overlaps accessible peaks/cCREs in the relevant cell type, whether the linked gene is expressed in the same cell type, and whether perturbation/eQTL/MPRA/rE2G data support the link.",
        "",
        "Reference patterns used to shape this skill:",
        "",
        "- Stuart et al., Cell 2019, Comprehensive Integration of Single-Cell Data, DOI: 10.1016/j.cell.2019.05.031.",
        "- Hao et al., Cell 2021, Integrated analysis of multimodal single-cell data, DOI: 10.1016/j.cell.2021.04.048.",
        "- Granja et al., Nature Genetics 2021, ArchR scalable integrative scATAC-seq analysis, DOI: 10.1038/s41588-021-00790-6.",
        "- Ashuach et al., Nature Methods 2023, MultiVI multimodal data integration, DOI: 10.1038/s41592-023-01909-9.",
        "",
    ]
    for name, skill in SKILLS.items():
        lines.extend(
            [
                f"## {skill['title']}",
                "",
                skill["purpose"],
                "",
                "Preferred inputs:",
                *[f"- {item}" for item in skill["preferred_files"]],
                "",
                "Analysis workflow:",
                *[f"{index}. {step}" for index, step in enumerate(skill["analysis_steps"], start=1)],
                "",
            ]
        )
    lines.extend(
        [
            "## Reuse Rules",
            "",
            "- For scRNA-seq, prefer count matrices or h5ad objects and analyze with Scanpy-compatible AnnData.",
            "- For scATAC-seq, prefer fragments plus peaks/cell metadata and analyze with ArchR, Signac, or SnapATAC-style workflows.",
            "- For Perturb-seq, require both expression and perturbation assignment metadata before modeling effects.",
            "- For IGVF variant interpretation, connect scATAC peaks, scRNA gene expression, Perturb-seq target effects, IGVF Catalog variant-gene evidence, and ENCODE reference context.",
            "- Do not download large files until manifest rows are reviewed.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")
    logging.info("Wrote playbook: %s", path)
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Single-cell data search and analysis skills.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    for command in ("smoke", "manifest"):
        subparser = subparsers.add_parser(command, help=f"Run {command} searches and write outputs.")
        subparser.add_argument("--skill", default="all", choices=["all", *SKILLS.keys()])
        subparser.add_argument("--source", default="encode", choices=["encode", "portal", "both"])
        subparser.add_argument("--limit", type=int, default=5)
        subparser.add_argument("--label", default=command)

    download_parser = subparsers.add_parser("download-examples", help="Download small public IGVF single-cell examples.")
    download_parser.add_argument("--force", action="store_true", help="Re-download examples even if local files already exist.")

    analyze_parser = subparsers.add_parser("analyze-examples", help="Analyze local IGVF scRNA/scATAC examples and write plots.")
    analyze_parser.add_argument("--max-cells", type=int, default=12000, help="Maximum cells to draw in the annotation-derived multiome embedding.")

    subparsers.add_parser("write-playbook", help="Write the single-cell analysis skill document.")

    args = parser.parse_args(argv)
    setup_logging()

    if args.command == "write-playbook":
        path = write_playbook()
        print(f"Wrote {path}")
        return 0

    if args.command == "download-examples":
        rows, manifest_path = download_examples(force=args.force)
        print(f"Wrote example download manifest: {manifest_path}")
        for row in rows:
            print(f"{row['accession']} {row['status']} {row['bytes']} bytes -> {row['path']}")
        return 0 if all(row["status"] in {"exists", "downloaded"} for row in rows) else 1

    if args.command == "analyze-examples":
        missing = [item for item in EXAMPLE_FILES if not Path(item["path"]).exists()]
        if missing:
            print("Missing example files. Run: python3 Scripts/single_cell_data_skills.py download-examples", file=sys.stderr)
            for item in missing:
                print(f"- {item['accession']}: {item['path']}", file=sys.stderr)
            return 1
        summary, report_path, json_path = analyze_examples(max_cells=args.max_cells)
        print(f"Wrote report: {report_path}")
        print(f"Wrote JSON: {json_path}")
        print(f"scRNA genes: {summary['scrna']['genes_total']:,}")
        print(f"scATAC peaks: {summary['scatac']['peaks_total']:,}")
        print(f"multiome cells summarized: {summary['multiome_annotations']['cells_total']:,}")
        return 0

    summaries, manifests = run_searches(args)
    manifest_path = write_manifest(manifests, f"{args.label}_{args.skill}_{args.source}")
    report_path = write_report(summaries, manifest_path, f"{args.label}_{args.skill}_{args.source}")
    save_json(f"single_cell_{args.label}_{args.skill}_{args.source}_summary", summaries)
    print(f"Wrote manifest: {manifest_path}")
    print(f"Wrote report: {report_path}")
    successful_skills = {
        summary["skill"]
        for summary in summaries
        if 200 <= summary["http_status"] < 400 and summary["returned_rows"] > 0
    }
    expected_skills = set(skill_names(args.skill))
    return 0 if expected_skills.issubset(successful_skills) else 1


if __name__ == "__main__":
    raise SystemExit(main())
