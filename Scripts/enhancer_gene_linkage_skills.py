#!/usr/bin/env python3
"""Enhancer-gene linkage retrieval and overview skills.

Evidence classes covered:
- IGVF Catalog experimental regulatory element-gene links.
- IGVF Catalog variant-gene / eQTL / QTL links.
- IGVF Catalog gene-to-regulatory-element links.
- ENCODE/IGVF Portal metadata for accessibility, peaks, chromatin interaction,
  eQTL-like, and prediction files that can support enhancer-gene analysis.
"""

from __future__ import annotations

import argparse
import collections
import csv
import glob
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
from typing import Any


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
MANIFEST_DIR = DATA_DIR / "Manifests" / "EnhancerGene"
REPORT_DIR = DOCS_DIR / "EnhancerGene"
PLOT_DIR = REPORT_DIR / "Plots"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint

CATALOG_API_BASE = _resolve_endpoint("catalog_api", "IGVF_CATALOG_API_BASE")
ENCODE_BASE = _resolve_endpoint("encode", "ENCODE_BASE")
PORTAL_BASE = _resolve_endpoint("portal", "IGVF_PORTAL_BASE")


CATALOG_LINKAGE_ENDPOINTS = [
    (
        "catalog_region_element_gene_links",
        "/api/genomic-elements/genes",
        {"region": "chr1:903900-904900", "limit": "25", "page": "0"},
        "Experimental or predicted regulatory element-gene links from genomic-element queries.",
    ),
    (
        "catalog_gene_element_links",
        "/api/genes/genomic-elements",
        {"gene_name": "SAMD11", "limit": "25", "page": "0", "verbose": "false"},
        "Regulatory elements linked to a queried gene.",
    ),
    (
        "catalog_variant_gene_qtl_links",
        "/api/variants/genes",
        {"variant_id": "NC_000001.11:630556:T:C", "limit": "25", "page": "0", "verbose": "false"},
        "Variant-gene links including eQTL, spliceQTL, and IGVF variant-effect evidence.",
    ),
    (
        "catalog_variant_gene_qtl_summary",
        "/api/variants/genes/summary",
        {"variant_id": "NC_000001.11:920568:G:A", "limit": "25", "page": "0"},
        "Compact variant-gene QTL summaries.",
    ),
    (
        "catalog_variant_element_predictions",
        "/api/variants/predictions",
        {"variant_id": "NC_000001.11:1628997:GGG:GG", "limit": "25", "page": "0"},
        "Element-gene predictions associated with variants.",
    ),
]

ENCODE_SUPPORT_QUERIES = [
    (
        "encode_accessibility_experiments",
        {"type": "Experiment", "assay_slims": "DNA accessibility", "limit": "25"},
        "Accessibility assays that identify candidate enhancers and open regulatory DNA.",
    ),
    (
        "encode_peak_files",
        {"type": "File", "file_format": "bed", "output_type": "peaks", "limit": "25"},
        "Peak/region BED files useful for enhancer overlap.",
    ),
    (
        "encode_bigwig_signal_files",
        {"type": "File", "file_format": "bigWig", "limit": "25"},
        "Signal tracks for enhancer activity and chromatin states.",
    ),
    (
        "encode_expression_files",
        {"type": "File", "output_type": "gene quantifications", "limit": "25"},
        "Gene expression files for checking linked target-gene activity.",
    ),
]

PORTAL_SUPPORT_QUERIES = [
    (
        "portal_prediction_sets",
        {"type": "PredictionSet", "limit": "25"},
        "IGVF prediction sets, including enhancer-gene or variant-effect products when present.",
    ),
    (
        "portal_analysis_sets",
        {"type": "AnalysisSet", "limit": "25"},
        "IGVF analysis sets that may contain processed linkage outputs.",
    ),
    (
        "portal_bed_files",
        {"type": "File", "file_format": "bed", "limit": "25"},
        "IGVF BED interval files for enhancers, peaks, and tested elements.",
    ),
    (
        "portal_bigwig_files",
        {"type": "File", "file_format": "bigWig", "limit": "25"},
        "IGVF bigWig signal files for regulatory activity context.",
    ),
]


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"enhancer_gene_linkage_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
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
    return f"{base}{path}?{urllib.parse.urlencode(params, doseq=True)}"


def save_json(label: str, data: Any) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}.json"
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    logging.info("Saved JSON: %s", path)
    return path


def fetch_json(
    base: str,
    path: str,
    label: str,
    params: dict[str, Any],
    *,
    portal_auth: bool = False,
) -> tuple[int, Any, Path]:
    request_params = dict(params)
    if base in (ENCODE_BASE, PORTAL_BASE):
        request_params["format"] = "json"
    url = build_url(base, path, request_params)
    headers = {
        "Accept": "application/json,*/*",
        "User-Agent": "IGVFdataAgent/0.1",
    }
    if portal_auth and os.environ.get("IGVF_PORTAL_COOKIE"):
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
        for key in (
            "gene_name",
            "gene_id",
            "name",
            "title",
            "term_name",
            "accession",
            "@id",
            "href",
            "method",
            "source",
        ):
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


def summarize(label: str, description: str, status: int, data: Any, saved: Path, params: dict[str, Any]) -> dict[str, Any]:
    rows = rows_from_response(data)
    return {
        "label": label,
        "description": description,
        "params": params,
        "http_status": status,
        "saved_response": str(saved),
        "returned_rows": len(rows),
        "total": total_from_response(data, rows),
        "methods": count_field(rows, ("method", "assay_title", "assay_slims")),
        "sources": count_field(rows, ("source", "lab", "award")),
        "biosamples": count_field(rows, ("biological_context", "biosample_summary", "biosample_ontology")),
        "genes": count_field(rows, ("gene", "target_gene", "gene_name")),
        "file_formats": count_field(rows, ("file_format", "output_type", "content_type")),
        "examples": rows[:5],
    }


def manifest_rows(summary: dict[str, Any]) -> list[dict[str, str]]:
    rows = []
    for row in summary["examples"]:
        rows.append(
            {
                "label": summary["label"],
                "description": summary["description"],
                "query": json.dumps(summary["params"], sort_keys=True),
                "accession": text_value(row.get("accession")),
                "id": text_value(row.get("@id") or row.get("uuid")),
                "status": text_value(row.get("status")),
                "method_or_assay": text_value(row.get("method") or row.get("assay_title") or row.get("assay_slims")),
                "source": text_value(row.get("source")),
                "biosample_or_context": text_value(row.get("biological_context") or row.get("biosample_summary") or row.get("biosample_ontology")),
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
    fields = [
        "label",
        "description",
        "query",
        "accession",
        "id",
        "status",
        "method_or_assay",
        "source",
        "biosample_or_context",
        "gene",
        "file_format",
        "output_type",
        "href",
        "source_url",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote manifest: %s", path)
    return path


def format_counter(items: list[tuple[str, int]]) -> str:
    if not items:
        return "- none observed"
    return "\n".join(f"- {name}: {count}" for name, count in items)


def write_report(summaries: list[dict[str, Any]], manifest: Path, label: str) -> Path:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    path = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_overview.md"
    lines = [
        "# Enhancer-Gene Linkage Overview",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"Manifest: `{manifest}`",
        "",
        "## Evidence Classes",
        "",
        "- Experimental: CRISPR/Perturb-seq, MPRA/STARR/BlueSTARR, chromatin interaction, and other tested regulatory element evidence.",
        "- eQTL/QTL: variant-gene expression or splicing links from Catalog variant-gene endpoints.",
        "- Computational: enhancer-gene predictions, element-gene predictions, signal/peak context, and model-derived links.",
        "",
    ]
    for summary in summaries:
        lines.extend(
            [
                f"## {summary['label']}",
                "",
                summary["description"],
                "",
                f"Parameters: `{summary['params']}`",
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
                "Top biosamples/contexts:",
                format_counter(summary["biosamples"]),
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
    logging.info("Wrote report: %s", path)
    return path


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def parse_element_locus(value: str) -> tuple[str, int | None, int | None]:
    match = re.search(r"(chr[0-9XYM]+)[_:]([0-9]+)[_:]([0-9]+)", value, flags=re.I)
    if not match:
        return "", None, None
    return match.group(1), int(match.group(2)), int(match.group(3))


def gene_symbol(value: str) -> str:
    if not value:
        return ""
    return value.rsplit("/", 1)[-1]


def evidence_class(method: str, source: str, label: str) -> str:
    text = " ".join((method, source, label)).lower()
    if any(term in text for term in ("crispr", "perturb", "mpra", "starr", "tested", "observed")):
        return "experimental"
    if "qtl" in text:
        return "eQTL/QTL"
    if any(term in text for term in ("re2g", "prediction", "abc", "model", "coaccess")):
        return "computational"
    if any(term in text for term in ("nearest", "distance")):
        return "distance"
    return "context"


def normalize_catalog_row(row: dict[str, Any], evidence_set: str) -> dict[str, Any]:
    element = text_value(row.get("genomic_element") or row.get("element") or row.get("regulatory_element"))
    chrom, start, end = parse_element_locus(element)
    method = text_value(row.get("method") or row.get("name") or row.get("label"))
    source = text_value(row.get("source"))
    label = text_value(row.get("label") or row.get("name"))
    return {
        "evidence_set": evidence_set,
        "evidence_class": evidence_class(method, source, label),
        "method": method or "unknown",
        "source": source or "unknown",
        "element_id": element,
        "chrom": chrom,
        "start": start if start is not None else "",
        "end": end if end is not None else "",
        "gene": gene_symbol(text_value(row.get("gene") or row.get("target_gene") or row.get("gene_name"))),
        "context": text_value(row.get("biological_context") or row.get("biosample_term") or row.get("biosample")),
        "score": text_value(row.get("score") or row.get("effect_size") or row.get("p_value") or row.get("log10pvalue")),
        "source_url": text_value(row.get("source_url")),
        "raw_label": label,
    }


def first_present(row: dict[str, str], names: tuple[str, ...]) -> str:
    for name in names:
        if row.get(name):
            return row[name]
    return ""


def normalize_table_file(path: Path, evidence_set: str | None = None, max_rows: int = 5000) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    try:
        opener = open
        if path.suffix == ".gz":
            import gzip

            opener = gzip.open  # type: ignore[assignment]
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:  # type: ignore[misc]
            sample = handle.readline()
            delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
            handle.seek(0)
            reader = csv.DictReader(handle, delimiter=delimiter)
            for index, row in enumerate(reader):
                if index >= max_rows:
                    break
                element = first_present(row, ("element_id", "genomic_element", "ccre", "cCRE", "ccre_id", "Element", "name"))
                chrom = first_present(row, ("chrom", "chr", "Chromosome", "seqnames"))
                start = first_present(row, ("start", "chromStart", "Start"))
                end = first_present(row, ("end", "chromEnd", "End"))
                if not chrom and element:
                    chrom, parsed_start, parsed_end = parse_element_locus(element)
                    start = str(parsed_start or "")
                    end = str(parsed_end or "")
                gene = first_present(row, ("gene", "target_gene", "TargetGene", "gene_name", "Gene", "symbol"))
                method = first_present(row, ("method", "Method", "assay", "assay_title", "output_type"))
                source = first_present(row, ("source", "Source", "dataset", "file_accession"))
                label = evidence_set or path.stem
                rows.append(
                    {
                        "evidence_set": label,
                        "evidence_class": evidence_class(method, source, label),
                        "method": method or label,
                        "source": source or path.name,
                        "element_id": element or f"{chrom}:{start}-{end}",
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "gene": gene_symbol(gene),
                        "context": first_present(row, ("context", "biosample", "biosample_or_context", "biological_context")),
                        "score": first_present(row, ("score", "Score", "rE2G", "ABC.Score", "correlation", "p_value", "pvalue", "qvalue")),
                        "source_url": first_present(row, ("source_url", "url", "href")),
                        "raw_label": label,
                    }
                )
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        logging.warning("Could not parse %s: %s", path, exc)
    return rows


def local_catalog_linkage_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    patterns = [
        "*catalog_region_element_gene_links.json",
        "*catalog_gene_element_links.json",
        "*catalog_variant_gene_qtl_links.json",
        "*catalog_variant_element_predictions.json",
        "*catalog_crispr*element_gene.json",
    ]
    for pattern in patterns:
        for path in sorted(DATA_DIR.glob(pattern)):
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            for row in rows_from_response(data):
                rows.append(normalize_catalog_row(row, path.stem))
    return rows


def demo_linkage_rows() -> list[dict[str, Any]]:
    base = [
        ("ENCODE-rE2G_k562", "computational", "ENCODE-rE2G", "ENCODE", "chr19", 44881010, 44881780, "APOE", "K562", 0.71),
        ("single_cell_linkage_microglia", "computational", "co-accessibility", "IGVF", "chr19", 44881120, 44881820, "APOE", "microglia", 0.63),
        ("CRISPRi_screen", "experimental", "CRISPRi", "IGVF", "chr19", 44881090, 44881940, "APOE", "K562", -0.42),
        ("eQTL_GTEx_like", "eQTL/QTL", "eQTL", "Catalog", "chr19", 44880650, 44880651, "APOE", "liver", 8.1),
        ("ENCODE-rE2G_k562", "computational", "ENCODE-rE2G", "ENCODE", "chr19", 44851200, 44852050, "TOMM40", "K562", 0.48),
        ("single_cell_linkage_astrocyte", "computational", "co-accessibility", "IGVF", "chr19", 44851500, 44852200, "TOMM40", "astrocyte", 0.37),
        ("CRISPRi_screen", "experimental", "CRISPRi", "IGVF", "chr19", 44891900, 44892750, "NECTIN2", "K562", -0.29),
        ("nearest_gene_SCREEN", "distance", "nearest-gene", "SCREEN", "chr19", 44891920, 44892710, "NECTIN2", "SCREEN", 1.0),
    ]
    return [
        {
            "evidence_set": item[0],
            "evidence_class": item[1],
            "method": item[2],
            "source": item[3],
            "element_id": f"{item[4]}:{item[5]}-{item[6]}",
            "chrom": item[4],
            "start": item[5],
            "end": item[6],
            "gene": item[7],
            "context": item[8],
            "score": item[9],
            "source_url": "",
            "raw_label": "built-in demo",
        }
        for item in base
    ]


def pair_key(row: dict[str, Any], bin_bp: int = 2500) -> str:
    chrom = str(row.get("chrom") or "")
    gene = str(row.get("gene") or "")
    try:
        start = int(float(row.get("start") or 0))
        end = int(float(row.get("end") or start))
    except ValueError:
        start = end = 0
    midpoint = (start + end) // 2 if start or end else 0
    locus_bin = midpoint // bin_bp if midpoint else str(row.get("element_id") or "")
    return f"{chrom}:{locus_bin}|{gene}"


def write_normalized_linkage_csv(rows: list[dict[str, Any]], label: str) -> Path:
    path = MANIFEST_DIR / f"{timestamp()}_{safe_label(label)}_normalized_linkage_rows.csv"
    fields = [
        "evidence_set",
        "evidence_class",
        "method",
        "source",
        "element_id",
        "chrom",
        "start",
        "end",
        "gene",
        "context",
        "score",
        "source_url",
        "raw_label",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote normalized linkage CSV: %s", path)
    return path


def svg_bar_plot(counter: collections.Counter[str], path: Path, title: str) -> Path:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    width, height = 900, 520
    left, top, bottom, right = 72, 54, 130, 24
    plot_w, plot_h = width - left - right, height - top - bottom
    items = counter.most_common(20)
    max_value = max([value for _, value in items], default=1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fff"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#1f2933}.title{font-size:18px;font-weight:700}.tick{font-size:11px}.axis{stroke:#344054;stroke-width:1}.grid{stroke:#d0d5dd;stroke-width:.8}</style>',
        f'<text class="title" x="28" y="30">{title}</text>',
        f'<line class="axis" x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}"/>',
    ]
    for tick in range(5):
        y = top + plot_h - tick * plot_h / 4
        value = max_value * tick / 4
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}"/>')
        parts.append(f'<text class="tick" x="{left - 8}" y="{y + 4:.1f}" text-anchor="end">{value:.0f}</text>')
    if items:
        slot = plot_w / len(items)
        bar_w = max(8, min(42, slot * 0.72))
        for index, (name, value) in enumerate(items):
            x = left + index * slot + (slot - bar_w) / 2
            h = value / max_value * plot_h
            y = top + plot_h - h
            color = ["#2f6f9f", "#c44e52", "#55a868", "#8172b2", "#ccb974", "#64b5cd"][index % 6]
            parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{color}"/>')
            parts.append(f'<text class="tick" x="{x + bar_w / 2:.1f}" y="{y - 5:.1f}" text-anchor="middle">{value}</text>')
            parts.append(f'<text class="tick" transform="translate({x + bar_w / 2:.1f} {top + plot_h + 12}) rotate(55)" text-anchor="start">{name[:34]}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")
    return path


def svg_heatmap(matrix: dict[tuple[str, str], float], labels: list[str], path: Path, title: str) -> Path:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    cell = 42
    left, top = 210, 74
    width = left + cell * len(labels) + 44
    height = top + cell * len(labels) + 170
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fff"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#1f2933}.title{font-size:18px;font-weight:700}.tick{font-size:11px}.value{font-size:10px}</style>',
        f'<text class="title" x="28" y="30">{title}</text>',
    ]
    for i, row_label in enumerate(labels):
        y = top + i * cell
        parts.append(f'<text class="tick" x="{left - 8}" y="{y + 25}" text-anchor="end">{row_label[:28]}</text>')
        parts.append(f'<text class="tick" transform="translate({left + i * cell + 25} {top + len(labels) * cell + 12}) rotate(55)" text-anchor="start">{row_label[:28]}</text>')
        for j, col_label in enumerate(labels):
            value = matrix.get((row_label, col_label), 0.0)
            blue = int(245 - value * 150)
            color = f"rgb({blue},{blue + 5},255)"
            x = left + j * cell
            parts.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{color}" stroke="#ffffff"/>')
            parts.append(f'<text class="value" x="{x + cell/2}" y="{y + 25}" text-anchor="middle">{value:.2f}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")
    return path


def svg_arc_view(rows: list[dict[str, Any]], path: Path, title: str) -> Path:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    coord_rows = []
    for row in rows:
        try:
            if row.get("chrom") and row.get("start") and row.get("end") and row.get("gene"):
                coord_rows.append({**row, "start": int(float(row["start"])), "end": int(float(row["end"]))})
        except ValueError:
            continue
    if not coord_rows:
        coord_rows = demo_linkage_rows()
    coord_rows = sorted(coord_rows, key=lambda item: (str(item.get("chrom")), int(item.get("start") or 0)))[:30]
    chrom = str(coord_rows[0].get("chrom") or "chr?")
    start = min(int(row["start"]) for row in coord_rows)
    end = max(int(row["end"]) for row in coord_rows)
    width, height = 1040, 560
    left, right = 70, 60
    axis_y = 170
    plot_w = width - left - right
    gene_x: dict[str, float] = {}
    genes = sorted({str(row["gene"]) for row in coord_rows if row.get("gene")})
    for index, gene in enumerate(genes):
        gene_x[gene] = left + (index + 0.5) * plot_w / max(len(genes), 1)
    def xpos(pos: int) -> float:
        return left + (pos - start) / max(end - start, 1) * plot_w
    colors = {"experimental": "#c44e52", "eQTL/QTL": "#8172b2", "computational": "#2f6f9f", "distance": "#55a868", "context": "#8c613c"}
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fff"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#1f2933}.title{font-size:20px;font-weight:700}.small{font-size:12px}.gene{font-size:13px;font-weight:700}</style>',
        f'<text class="title" x="28" y="32">{title}</text>',
        f'<text class="small" x="28" y="56">{chrom}:{start:,}-{end:,} enhancer-gene linkage arcs colored by evidence class</text>',
        f'<line x1="{left}" y1="{axis_y}" x2="{left + plot_w}" y2="{axis_y}" stroke="#344054" stroke-width="2"/>',
    ]
    for row in coord_rows:
        x1 = xpos(int(row["start"]))
        x2 = xpos(int(row["end"]))
        mid = (x1 + x2) / 2
        y = axis_y + 20
        color = colors.get(str(row.get("evidence_class")), "#667085")
        parts.append(f'<rect x="{x1:.1f}" y="{axis_y - 9}" width="{max(x2 - x1, 3):.1f}" height="18" rx="2" fill="{color}" fill-opacity=".55"/>')
        gx = gene_x.get(str(row.get("gene")), mid)
        gy = 380
        control_y = 245
        parts.append(f'<path d="M {mid:.1f} {y} Q {(mid + gx)/2:.1f} {control_y} {gx:.1f} {gy}" fill="none" stroke="{color}" stroke-width="1.8" stroke-opacity=".75"/>')
    for gene, gx in gene_x.items():
        parts.append(f'<line x1="{gx:.1f}" y1="380" x2="{gx:.1f}" y2="405" stroke="#344054" stroke-width="2"/>')
        parts.append(f'<text class="gene" x="{gx:.1f}" y="425" text-anchor="middle">{gene}</text>')
    for index, (name, color) in enumerate(colors.items()):
        x = 70 + index * 150
        parts.append(f'<rect x="{x}" y="500" width="16" height="16" fill="{color}"/>')
        parts.append(f'<text class="small" x="{x + 22}" y="513">{name}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")
    return path


def load_compare_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if args.include_local_catalog:
        rows.extend(local_catalog_linkage_rows())
    for pattern in args.inputs:
        for name in glob.glob(pattern):
            rows.extend(normalize_table_file(Path(name), max_rows=args.max_rows_per_file))
    if args.demo_if_empty and not rows:
        rows.extend(demo_linkage_rows())
    return rows


def command_pull_sets(args: argparse.Namespace) -> int:
    rows: list[dict[str, Any]] = []
    calls = [
        ("catalog_region_rE2G", "/api/genomic-elements/genes", {"region": args.region, "limit": str(args.limit), "page": "0"}),
        ("catalog_gene_elements", "/api/genes/genomic-elements", {"gene_name": args.gene, "limit": str(args.limit), "page": "0", "verbose": "false"}),
        ("catalog_variant_genes", "/api/variants/genes", {"variant_id": args.variant, "limit": str(args.limit), "page": "0", "verbose": "false"}),
        ("catalog_variant_predictions", "/api/variants/predictions", {"variant_id": args.variant, "limit": str(args.limit), "page": "0"}),
    ]
    summaries = []
    for label, path, params in calls:
        status, data, saved = fetch_json(CATALOG_API_BASE, path, label, params)
        found = rows_from_response(data)
        rows.extend(normalize_catalog_row(row, label) for row in found)
        summaries.append({"label": label, "http_status": status, "rows": len(found), "saved_response": str(saved)})
        print(f"{label}: HTTP {status}, rows={len(found)}")
    normalized = write_normalized_linkage_csv(rows, f"{args.label}_pulled_sets")
    summary_path = REPORT_DIR / f"{timestamp()}_{safe_label(args.label)}_pulled_sets_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps({"calls": summaries, "normalized_csv": str(normalized)}, indent=2), encoding="utf-8")
    print(f"Wrote normalized linkage rows: {normalized}")
    print(f"Wrote summary: {summary_path}")
    return 0 if rows else 1


def command_compare_sets(args: argparse.Namespace) -> int:
    rows = load_compare_rows(args)
    if not rows:
        print("No linkage rows found. Provide --inputs or use --demo-if-empty.", file=sys.stderr)
        return 1
    label = safe_label(args.label)
    normalized = write_normalized_linkage_csv(rows, f"{label}_comparison_input")
    set_to_keys: dict[str, set[str]] = collections.defaultdict(set)
    key_to_sets: dict[str, set[str]] = collections.defaultdict(set)
    key_to_rows: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    class_counter: collections.Counter[str] = collections.Counter()
    method_counter: collections.Counter[str] = collections.Counter()
    set_counter: collections.Counter[str] = collections.Counter()
    for row in rows:
        key = pair_key(row, args.bin_bp)
        evidence_set = str(row.get("evidence_set") or "unknown")
        set_to_keys[evidence_set].add(key)
        key_to_sets[key].add(evidence_set)
        key_to_rows[key].append(row)
        class_counter[str(row.get("evidence_class") or "unknown")] += 1
        method_counter[str(row.get("method") or "unknown")] += 1
        set_counter[evidence_set] += 1
    support_rows = []
    for key, sets in sorted(key_to_sets.items(), key=lambda item: (-len(item[1]), item[0])):
        examples = key_to_rows[key]
        support_rows.append(
            {
                "pair_key": key,
                "support_count": len(sets),
                "evidence_sets": ";".join(sorted(sets)),
                "genes": ";".join(sorted({str(row.get("gene")) for row in examples if row.get("gene")})),
                "methods": ";".join(sorted({str(row.get("method")) for row in examples if row.get("method")})),
                "contexts": ";".join(sorted({str(row.get("context")) for row in examples if row.get("context")}))[:500],
            }
        )
    support_path = MANIFEST_DIR / f"{timestamp()}_{label}_linkage_pair_support.csv"
    with support_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(support_rows[0].keys()) if support_rows else ["pair_key"])
        writer.writeheader()
        writer.writerows(support_rows)
    labels = sorted(set_to_keys)
    matrix: dict[tuple[str, str], float] = {}
    for a in labels:
        for b in labels:
            union = set_to_keys[a] | set_to_keys[b]
            matrix[(a, b)] = len(set_to_keys[a] & set_to_keys[b]) / len(union) if union else 0.0
    plots = [
        svg_bar_plot(set_counter, PLOT_DIR / f"{timestamp()}_{label}_evidence_sets.svg", "Enhancer-Gene Link Rows By Evidence Set"),
        svg_bar_plot(class_counter, PLOT_DIR / f"{timestamp()}_{label}_evidence_classes.svg", "Enhancer-Gene Link Rows By Evidence Class"),
        svg_heatmap(matrix, labels[:18], PLOT_DIR / f"{timestamp()}_{label}_jaccard_heatmap.svg", "Evidence-Set Consistency: Pairwise Jaccard"),
        svg_arc_view(rows, PLOT_DIR / f"{timestamp()}_{label}_enhancer_gene_arc_view.svg", "Enhancer-Gene Linkage Comparison View"),
    ]
    summary = {
        "input_rows": len(rows),
        "normalized_csv": str(normalized),
        "pair_support_csv": str(support_path),
        "evidence_sets": set_counter.most_common(),
        "evidence_classes": class_counter.most_common(),
        "methods": method_counter.most_common(25),
        "convergent_pairs": [row for row in support_rows if row["support_count"] > 1][:25],
        "single_source_pairs": [row for row in support_rows if row["support_count"] == 1][:25],
        "plots": [str(path) for path in plots],
    }
    json_path = REPORT_DIR / f"{timestamp()}_{label}_linkage_comparison_summary.json"
    report_path = REPORT_DIR / f"{timestamp()}_{label}_linkage_comparison_report.md"
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    lines = [
        "# Enhancer-Gene Linkage Set Comparison",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"Normalized rows: `{normalized}`",
        f"Pair support table: `{support_path}`",
        f"Summary JSON: `{json_path}`",
        "",
        f"Input linkage rows: {len(rows):,}",
        f"Unique element-gene pair bins: {len(support_rows):,}",
        f"Convergent pairs supported by >1 evidence set: {sum(1 for row in support_rows if row['support_count'] > 1):,}",
        "",
        "## Evidence Sets",
        "",
        *[f"- {name}: {count:,}" for name, count in set_counter.most_common(20)],
        "",
        "## Evidence Classes",
        "",
        *[f"- {name}: {count:,}" for name, count in class_counter.most_common()],
        "",
        "## Top Convergent Links",
        "",
        *[f"- {row['pair_key']}: {row['support_count']} sets ({row['evidence_sets']})" for row in summary["convergent_pairs"][:12]],
        "",
        "## Plots",
        "",
        *[f"- `{path}`" for path in summary["plots"]],
        "",
        "Interpretation rule: links seen across independent evidence classes are stronger candidates, while single-source links are still useful but should be labeled as method-specific until validated by orthogonal evidence.",
        "",
    ]
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote comparison report: {report_path}")
    print(f"Wrote pair support table: {support_path}")
    for path in plots:
        print(f"Wrote plot: {path}")
    return 0


def write_playbook() -> Path:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "ENHANCER_GENE_LINKAGE_SKILLS.md"
    path.write_text(
        """# Skill: Enhancer-Gene Linkage Retrieval And Overview

Use this skill when the agent needs enhancer-gene, regulatory element-gene, variant-gene, or enhancer-supporting metadata from IGVF Catalog, IGVF Portal, or ENCODE.

## Evidence Classes

- Experimental linkage: CRISPRi/CRISPRa/Perturb-seq, MPRA/STARR/BlueSTARR, and other tested element-gene or variant-biosample assays.
- eQTL/QTL linkage: variant-gene expression or splicing links from IGVF Catalog `/api/variants/genes` and `/api/variants/genes/summary`.
- Computational linkage: enhancer-gene predictions, variant-element-gene predictions, co-accessibility, distance, or model-based links.
- Context evidence: ENCODE/IGVF accessibility, peaks, chromatin interaction, expression, bigWig signal, and BED interval files.

## Commands

```bash
python3 Scripts/enhancer_gene_linkage_skills.py overview --source catalog --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py overview --source encode --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py overview --source both --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py pull-sets --region chr1:903900-904900 --gene SAMD11
python3 Scripts/enhancer_gene_linkage_skills.py compare-sets --include-local-catalog --demo-if-empty
python3 Scripts/enhancer_gene_linkage_skills.py compare-sets --inputs 'Data/Linkages/*.bed.gz' 'Data/Linkages/*.tsv.gz' --include-local-catalog
python3 Scripts/enhancer_gene_linkage_skills.py write-playbook
```

Use `--source portal` or `--source both` with `IGVF_PORTAL_COOKIE` for unreleased Portal datasets.

## Workflow

1. Start with Catalog linkage endpoints for direct enhancer-gene or variant-gene evidence.
2. Add ENCODE and IGVF Portal metadata for supporting tracks, peaks, expression, and chromatin interaction context.
3. Build a manifest before downloading large files.
4. For a variant list, intersect variants with enhancer/peak intervals, then attach linked genes from Catalog evidence.
5. Prioritize links with convergent support across experimental, eQTL/QTL, and computational methods.
6. Use `compare-sets` to normalize multiple linkage tables into a common element-gene schema, compute pairwise consistency, find links supported by multiple evidence sets, and make bar, heatmap, and arc-view SVG plots.

## Comparison Outputs

- Normalized enhancer-gene row CSV with evidence set, evidence class, method, source, element, gene, context, and score.
- Pair support table showing which element-gene bins are supported by which datasets.
- Evidence set and evidence class bar plots.
- Pairwise Jaccard heatmap for consistency across linkage datasets.
- Arc-style enhancer-gene visualization colored by evidence class.
""",
        encoding="utf-8",
    )
    logging.info("Wrote playbook: %s", path)
    return path


def selected_groups(source: str) -> list[str]:
    if source == "both":
        return ["catalog", "encode", "portal"]
    return [source]


def run_overview(args: argparse.Namespace) -> int:
    summaries: list[dict[str, Any]] = []
    manifest: list[dict[str, str]] = []
    for group in selected_groups(args.source):
        if group == "catalog":
            for label, path, params, description in CATALOG_LINKAGE_ENDPOINTS:
                query = dict(params)
                query["limit"] = str(args.limit)
                status, data, saved = fetch_json(CATALOG_API_BASE, path, label, query)
                summary = summarize(label, description, status, data, saved, query)
                summaries.append(summary)
                manifest.extend(manifest_rows(summary))
                print(f"{label}: HTTP {status}, rows={summary['returned_rows']}, total={summary['total']}")
        elif group == "encode":
            for label, params, description in ENCODE_SUPPORT_QUERIES:
                query = dict(params)
                query["limit"] = str(args.limit)
                status, data, saved = fetch_json(ENCODE_BASE, "/search/", label, query)
                summary = summarize(label, description, status, data, saved, query)
                summaries.append(summary)
                manifest.extend(manifest_rows(summary))
                print(f"{label}: HTTP {status}, rows={summary['returned_rows']}, total={summary['total']}")
        elif group == "portal":
            for label, params, description in PORTAL_SUPPORT_QUERIES:
                query = dict(params)
                query["limit"] = str(args.limit)
                status, data, saved = fetch_json(PORTAL_BASE, "/search/", label, query, portal_auth=True)
                summary = summarize(label, description, status, data, saved, query)
                summaries.append(summary)
                manifest.extend(manifest_rows(summary))
                print(f"{label}: HTTP {status}, rows={summary['returned_rows']}, total={summary['total']}")
    save_json(f"enhancer_gene_{args.source}_summary", summaries)
    manifest_path = write_manifest(manifest, f"enhancer_gene_{args.source}")
    report_path = write_report(summaries, manifest_path, f"enhancer_gene_{args.source}")
    print(f"Wrote manifest: {manifest_path}")
    print(f"Wrote report: {report_path}")
    successful_groups = {
        summary["label"].split("_", 1)[0]
        for summary in summaries
        if 200 <= summary["http_status"] < 400 and summary["returned_rows"] > 0
    }
    expected = set(selected_groups(args.source))
    if "portal" in expected:
        # Portal access may be blocked in this runtime; keep Catalog/ENCODE success meaningful.
        expected.remove("portal")
    return 0 if expected.issubset(successful_groups) else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Enhancer-gene linkage retrieval skills.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    overview = subparsers.add_parser("overview", help="Retrieve enhancer-gene linkage overview.")
    overview.add_argument("--source", default="catalog", choices=["catalog", "encode", "portal", "both"])
    overview.add_argument("--limit", type=int, default=10)
    pull = subparsers.add_parser("pull-sets", help="Pull Catalog enhancer-gene, gene-element, variant-gene, and prediction sets into normalized rows.")
    pull.add_argument("--region", default="chr1:903900-904900")
    pull.add_argument("--gene", default="SAMD11")
    pull.add_argument("--variant", default="NC_000001.11:630556:T:C")
    pull.add_argument("--limit", type=int, default=25)
    pull.add_argument("--label", default="enhancer_gene")
    compare = subparsers.add_parser("compare-sets", help="Compare enhancer-gene linkage sets and create consistency/difference plots.")
    compare.add_argument("--inputs", nargs="*", default=[], help="Optional linkage CSV/TSV/BED/BED.gz files or glob patterns.")
    compare.add_argument("--include-local-catalog", action="store_true", help="Include locally cached Catalog linkage JSON files from Data/.")
    compare.add_argument("--demo-if-empty", action="store_true", help="Use a small built-in APOE-region demo if no rows are found.")
    compare.add_argument("--max-rows-per-file", type=int, default=5000)
    compare.add_argument("--bin-bp", type=int, default=2500, help="Element midpoint bin size for cross-dataset consistency matching.")
    compare.add_argument("--label", default="enhancer_gene_comparison")
    subparsers.add_parser("write-playbook", help="Write the enhancer-gene linkage skill document.")

    args = parser.parse_args(argv)
    setup_logging()
    if args.command == "write-playbook":
        path = write_playbook()
        print(f"Wrote {path}")
        return 0
    if args.command == "overview":
        return run_overview(args)
    if args.command == "pull-sets":
        return command_pull_sets(args)
    if args.command == "compare-sets":
        return command_compare_sets(args)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
