#!/usr/bin/env python3
"""IGVF Portal and Catalog data overview and smoke-analysis skills.

The portal side is authentication-aware. Released portal data should be public;
set IGVF_PORTAL_COOKIE from a local logged-in browser session when unreleased
or protected portal metadata is needed.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import logging
import os
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

PORTAL_BASE = os.environ.get("IGVF_PORTAL_BASE", "https://data.igvf.org").rstrip("/")
CATALOG_API_BASE = os.environ.get(
    "IGVF_CATALOG_API_BASE", "https://api.catalogkg.igvf.org"
).rstrip("/")
ENCODE_BASE = os.environ.get("ENCODE_BASE", "https://www.encodeproject.org").rstrip("/")

PORTAL_OBJECT_TYPES = [
    "File",
    "MeasurementSet",
    "AnalysisSet",
    "PredictionSet",
    "ModelSet",
    "ConstructLibrarySet",
    "Sample",
    "Donor",
    "Software",
    "Document",
]

CATALOG_SMOKE_ENDPOINTS = [
    (
        "catalog_files_filesets",
        "/api/files-filesets",
        {"limit": "25", "offset": "0"},
        "Catalog file/fileset inventory used to connect portal filesets to KG-facing rows.",
    ),
    (
        "catalog_mpra_variants",
        "/api/variants/biosamples",
        {"method": "MPRA", "limit": "10", "page": "0"},
        "Measured noncoding variant activity evidence.",
    ),
    (
        "catalog_coding_variant_scores",
        "/api/genes/coding-variants/scores",
        {"gene_name": "TP53", "limit": "10", "page": "0"},
        "Coding variant score evidence for one well-known gene.",
    ),
    (
        "catalog_enhancer_gene_predictions",
        "/api/genomic-elements/genes",
        {"region": "chr1:903900-904900", "limit": "10", "page": "0"},
        "Regulatory element-gene evidence rows for a small genomic region.",
    ),
]

ENCODE_OBJECT_TYPES = [
    "Experiment",
    "File",
    "Biosample",
    "Dataset",
    "Annotation",
    "Reference",
    "Software",
    "Publication",
]

ENCODE_SMOKE_QUERIES = [
    (
        "encode_atac_seq_experiments",
        {"type": "Experiment", "assay_title": "ATAC-seq", "limit": "10"},
        "Chromatin accessibility experiments.",
    ),
    (
        "encode_chip_seq_experiments",
        {"type": "Experiment", "assay_slims": "DNA binding", "limit": "10"},
        "Transcription factor or histone ChIP-seq experiments.",
    ),
    (
        "encode_rna_seq_experiments",
        {"type": "Experiment", "assay_slims": "Transcription", "limit": "10"},
        "Transcriptome profiling experiments.",
    ),
    (
        "encode_dna_accessibility_experiments",
        {"type": "Experiment", "assay_slims": "DNA accessibility", "limit": "10"},
        "DNA accessibility experiments useful for regulatory variant interpretation.",
    ),
    (
        "encode_bigwig_files",
        {"type": "File", "file_format": "bigWig", "limit": "10"},
        "Signal tracks commonly used for regulatory annotation.",
    ),
    (
        "encode_bed_files",
        {"type": "File", "file_format": "bed", "limit": "10"},
        "Region/peak calls and genomic interval files.",
    ),
]


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"igvf_data_skills_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)
    return log_path


def safe_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in label)


def save_bytes(label: str, content: bytes, content_type: str) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    suffix = "json" if "json" in content_type else "txt"
    path = DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}.{suffix}"
    path.write_bytes(content)
    logging.info("Saved response: %s", path)
    return path


def save_json(label: str, data: Any) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}.json"
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    logging.info("Saved JSON: %s", path)
    return path


def build_url(base: str, path: str, params: dict[str, Any] | None = None) -> str:
    normalized = path if path.startswith("/") else f"/{path}"
    url = f"{base}{normalized}"
    if params:
        url = f"{url}?{urllib.parse.urlencode(params, doseq=True)}"
    return url


def fetch_json(
    base: str,
    path: str,
    label: str,
    params: dict[str, Any] | None = None,
    *,
    portal_auth: bool = False,
) -> tuple[int, Any, Path]:
    headers = {
        "Accept": "application/json,*/*",
        "User-Agent": "IGVFdataAgent/0.1",
    }
    if portal_auth and os.environ.get("IGVF_PORTAL_COOKIE"):
        headers["Cookie"] = os.environ["IGVF_PORTAL_COOKIE"]
    url = build_url(base, path, params)
    logging.info("Request: GET %s", url)
    req = urllib.request.Request(url, headers=headers, method="GET")
    try:
        with urllib.request.urlopen(req, timeout=45) as response:
            content = response.read()
            content_type = response.headers.get("Content-Type", "")
            saved = save_bytes(label, content, content_type)
            if "json" in content_type:
                return response.status, json.loads(content), saved
            return response.status, content.decode(errors="replace"), saved
    except urllib.error.HTTPError as exc:
        content = exc.read()
        saved = save_bytes(label, content, exc.headers.get("Content-Type", "text/plain"))
        return exc.code, content.decode(errors="replace"), saved
    except urllib.error.URLError as exc:
        message = {"network_error": str(exc.reason), "url": url}
        saved = save_json(label, message)
        return 0, message, saved


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
        for key in ("title", "name", "term_name", "accession", "@id", "uuid"):
            if key in value:
                return scalar_strings(value[key])
    return []


def count_field(rows: list[dict[str, Any]], field_names: tuple[str, ...], limit: int = 12) -> list[tuple[str, int]]:
    counter: collections.Counter[str] = collections.Counter()
    for row in rows:
        for field in field_names:
            for value in scalar_strings(row.get(field)):
                if value:
                    counter[value] += 1
    return counter.most_common(limit)


def summarize_rows(object_type: str, status: int, data: Any, saved: Path) -> dict[str, Any]:
    rows = rows_from_response(data)
    summary = {
        "object_type": object_type,
        "http_status": status,
        "saved_response": str(saved),
        "returned_rows": len(rows),
        "total": total_from_response(data, rows),
        "statuses": count_field(rows, ("status",)),
        "assays": count_field(rows, ("assay_title", "assay_titles", "assay_term_name")),
        "file_formats": count_field(rows, ("file_format", "file_format_type")),
        "content_types": count_field(rows, ("content_type", "content_types")),
        "file_sets": count_field(rows, ("file_set", "file_sets", "fileset", "filesets")),
        "labs": count_field(rows, ("lab", "award")),
        "examples": [
            {
                "accession": row.get("accession"),
                "id": row.get("@id") or row.get("uuid"),
                "title": row.get("title") or row.get("summary"),
                "status": row.get("status"),
            }
            for row in rows[:5]
        ],
    }
    return summary


def portal_search(object_type: str, limit: int) -> tuple[int, Any, Path]:
    return fetch_json(
        PORTAL_BASE,
        "/search/",
        f"portal_{object_type}_search",
        {"type": object_type, "format": "json", "limit": str(limit)},
        portal_auth=True,
    )


def encode_search(params: dict[str, Any], label: str) -> tuple[int, Any, Path]:
    request_params = dict(params)
    request_params["format"] = "json"
    return fetch_json(ENCODE_BASE, "/search/", label, request_params)


def catalog_fetch(path: str, label: str, params: dict[str, Any]) -> tuple[int, Any, Path]:
    return fetch_json(CATALOG_API_BASE, path, label, params)


def write_markdown_report(path: Path, title: str, sections: list[str]) -> None:
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    body = [f"# {title}", "", f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}", ""]
    for section in sections:
        body.append(section.rstrip())
        body.append("")
    path.write_text("\n".join(body), encoding="utf-8")
    logging.info("Wrote report: %s", path)


def format_counter(items: list[tuple[str, int]]) -> str:
    if not items:
        return "- none observed"
    return "\n".join(f"- {name}: {count}" for name, count in items)


def section_for_summary(summary: dict[str, Any]) -> str:
    return f"""## {summary['object_type']}

HTTP status: {summary['http_status']}
Returned rows: {summary['returned_rows']}
Reported total: {summary['total']}
Saved response: `{summary['saved_response']}`

Top statuses:
{format_counter(summary['statuses'])}

Top assays:
{format_counter(summary['assays'])}

Top file formats:
{format_counter(summary['file_formats'])}

Top content types:
{format_counter(summary['content_types'])}
"""


def section_for_encode_summary(summary: dict[str, Any]) -> str:
    return f"""## {summary['object_type']}

HTTP status: {summary['http_status']}
Returned rows: {summary['returned_rows']}
Reported total: {summary['total']}
Saved response: `{summary['saved_response']}`

Top statuses:
{format_counter(summary['statuses'])}

Top assays:
{format_counter(summary['assays'])}

Top biosample terms:
{format_counter(summary['biosamples'])}

Top organisms:
{format_counter(summary['organisms'])}

Top file formats:
{format_counter(summary['file_formats'])}

Top output/content types:
{format_counter(summary['content_types'])}
"""


def summarize_encode_rows(object_type: str, status: int, data: Any, saved: Path) -> dict[str, Any]:
    summary = summarize_rows(object_type, status, data, saved)
    rows = rows_from_response(data)
    summary["biosamples"] = count_field(
        rows,
        (
            "biosample_ontology",
            "biosample_summary",
            "biosample_term_name",
            "biosample_type",
        ),
    )
    summary["organisms"] = count_field(rows, ("organism", "replicates.library.biosample.organism"))
    summary["targets"] = count_field(rows, ("target", "target.label", "target.title"))
    summary["award_rfas"] = count_field(rows, ("award.rfa", "award", "lab"))
    return summary


def overview(args: argparse.Namespace) -> int:
    object_types = args.type or PORTAL_OBJECT_TYPES
    summaries = []
    for object_type in object_types:
        status, data, saved = portal_search(object_type, args.limit)
        summary = summarize_rows(object_type, status, data, saved)
        summaries.append(summary)
        print(f"{object_type}: HTTP {status}, rows={summary['returned_rows']}, total={summary['total']}")
    save_json("portal_overview_summary", summaries)
    report_path = DOCS_DIR / "IGVF_PORTAL_DATA_OVERVIEW.md"
    sections = [
        "## Scope\n\nThis overview samples IGVF Portal metadata object types through `/search/?type=...&format=json`. Released data should be public; unreleased datasets require a local authenticated session via `IGVF_PORTAL_COOKIE`.",
        *[section_for_summary(summary) for summary in summaries],
    ]
    write_markdown_report(report_path, "IGVF Portal Data Overview", sections)
    print(f"Wrote {report_path}")
    return 0 if all(200 <= item["http_status"] < 400 for item in summaries) else 1


def smoke(args: argparse.Namespace) -> int:
    object_types = args.type or ["File", "MeasurementSet", "AnalysisSet", "PredictionSet", "ModelSet"]
    summaries = []
    for object_type in object_types:
        status, data, saved = portal_search(object_type, args.limit)
        summaries.append(summarize_rows(object_type, status, data, saved))
    save_json("portal_smoke_summary", summaries)
    report_path = DOCS_DIR / "IGVF_PORTAL_SMOKE_ANALYSIS.md"
    intro = """## Major Smoke Classes

- Files: check file formats, content types, statuses, and file-to-fileset links.
- Measurement sets: check assay titles and biological contexts for raw/processed experimental inputs.
- Analysis and prediction sets: check computational outputs, model/pipeline-linked products, and downstream Catalog inputs.
- Model sets: check predictive model products and metadata needed to reproduce computational outputs.

The smoke analysis is intentionally shallow: it verifies access, metadata shape, counts, and representative examples before any large data download."""
    sections = [intro, *[section_for_summary(summary) for summary in summaries]]
    write_markdown_report(report_path, "IGVF Portal Smoke Analysis", sections)
    print(f"Wrote {report_path}")
    return 0 if all(200 <= item["http_status"] < 400 for item in summaries) else 1


def catalog_smoke(args: argparse.Namespace) -> int:
    summaries = []
    for label, path, params, description in CATALOG_SMOKE_ENDPOINTS:
        request_params = dict(params)
        if args.limit is not None:
            request_params["limit"] = str(args.limit)
        status, data, saved = catalog_fetch(path, label, request_params)
        rows = rows_from_response(data)
        summary = {
            "label": label,
            "description": description,
            "path": path,
            "params": request_params,
            "http_status": status,
            "saved_response": str(saved),
            "returned_rows": len(rows),
            "total": total_from_response(data, rows),
            "examples": rows[:3],
        }
        summaries.append(summary)
        print(f"{label}: HTTP {status}, rows={summary['returned_rows']}")
    save_json("catalog_smoke_summary", summaries)
    report_path = DOCS_DIR / "IGVF_CATALOG_SMOKE_ANALYSIS.md"
    sections = [
        "## Purpose\n\nThese public Catalog API checks exercise data classes that are derived from or linked back to IGVF Portal submissions: file/fileset inventory, measured noncoding activity, coding variant scores, and enhancer-gene predictions.",
    ]
    for summary in summaries:
        sections.append(
            f"""## {summary['label']}

{summary['description']}

Endpoint: `{summary['path']}`
Parameters: `{summary['params']}`
HTTP status: {summary['http_status']}
Returned rows: {summary['returned_rows']}
Saved response: `{summary['saved_response']}`
"""
        )
    write_markdown_report(report_path, "IGVF Catalog Smoke Analysis", sections)
    print(f"Wrote {report_path}")
    return 0 if all(200 <= item["http_status"] < 400 for item in summaries) else 1


def encode_overview(args: argparse.Namespace) -> int:
    object_types = args.type or ENCODE_OBJECT_TYPES
    summaries = []
    for object_type in object_types:
        status, data, saved = encode_search(
            {"type": object_type, "limit": str(args.limit)},
            f"encode_{object_type}_search",
        )
        summary = summarize_encode_rows(object_type, status, data, saved)
        summaries.append(summary)
        print(f"ENCODE {object_type}: HTTP {status}, rows={summary['returned_rows']}, total={summary['total']}")
    save_json("encode_overview_summary", summaries)
    report_path = DOCS_DIR / "ENCODE_DATA_OVERVIEW.md"
    sections = [
        "## Scope\n\nThis overview samples public ENCODE metadata through `/search/?type=...&format=json`. Use it to understand available experiments, files, biosamples, annotations, references, software, and publications before downloading data.",
        *[section_for_encode_summary(summary) for summary in summaries],
    ]
    write_markdown_report(report_path, "ENCODE Data Overview", sections)
    print(f"Wrote {report_path}")
    return 0 if all(200 <= item["http_status"] < 400 for item in summaries) else 1


def encode_smoke(args: argparse.Namespace) -> int:
    summaries = []
    for label, params, description in ENCODE_SMOKE_QUERIES:
        request_params = dict(params)
        if args.limit is not None:
            request_params["limit"] = str(args.limit)
        status, data, saved = encode_search(request_params, label)
        rows = rows_from_response(data)
        summary = summarize_encode_rows(label, status, data, saved)
        summary.update(
            {
                "label": label,
                "description": description,
                "params": request_params,
                "examples": rows[:3],
            }
        )
        summaries.append(summary)
        print(f"{label}: HTTP {status}, rows={summary['returned_rows']}, total={summary['total']}")
    save_json("encode_smoke_summary", summaries)
    report_path = DOCS_DIR / "ENCODE_SMOKE_ANALYSIS.md"
    sections = [
        "## Purpose\n\nThese public ENCODE checks exercise major input data classes for IGVF interpretation work: accessibility, DNA binding, transcription, signal tracks, and interval/peak files.",
    ]
    for summary in summaries:
        sections.append(
            f"""## {summary['label']}

{summary['description']}

Parameters: `{summary['params']}`
HTTP status: {summary['http_status']}
Returned rows: {summary['returned_rows']}
Reported total: {summary['total']}
Saved response: `{summary['saved_response']}`

Top assays:
{format_counter(summary['assays'])}

Top biosamples:
{format_counter(summary['biosamples'])}

Top file formats:
{format_counter(summary['file_formats'])}
"""
        )
    write_markdown_report(report_path, "ENCODE Smoke Analysis", sections)
    print(f"Wrote {report_path}")
    return 0 if all(200 <= item["http_status"] < 400 for item in summaries) else 1


def encode_export_csv(args: argparse.Namespace) -> int:
    params = {"type": args.type, "limit": str(args.limit)}
    for param in args.param:
        if "=" not in param:
            print(f"Expected KEY=VALUE for --param, got: {param}")
            return 2
        key, value = param.split("=", 1)
        params.setdefault(key, [])
        if isinstance(params[key], list):
            params[key].append(value)
        else:
            params[key] = [params[key], value]
    status, data, saved = encode_search(params, f"encode_{args.type}_export")
    rows = rows_from_response(data)
    output = DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_encode_{safe_label(args.type)}_metadata.csv"
    fields = [
        "accession",
        "@id",
        "uuid",
        "status",
        "assay_title",
        "biosample_summary",
        "file_format",
        "output_type",
        "assembly",
        "href",
        "summary",
        "description",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            flat = {field: ", ".join(scalar_strings(row.get(field))) for field in fields}
            writer.writerow(flat)
    print(f"Fetched HTTP {status} from {saved}")
    print(f"Wrote {output}")
    return 0 if 200 <= status < 400 else 1


def export_csv(args: argparse.Namespace) -> int:
    status, data, saved = portal_search(args.type, args.limit)
    rows = rows_from_response(data)
    output = DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_portal_{safe_label(args.type)}_metadata.csv"
    fields = ["accession", "@id", "uuid", "status", "file_format", "content_type", "assay_title", "summary"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            flat = {field: ", ".join(scalar_strings(row.get(field))) for field in fields}
            writer.writerow(flat)
    print(f"Fetched HTTP {status} from {saved}")
    print(f"Wrote {output}")
    return 0 if 200 <= status < 400 else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="IGVF Portal and Catalog data skills")
    subparsers = parser.add_subparsers(dest="command", required=True)

    overview_parser = subparsers.add_parser("overview", help="Sample major IGVF Portal object types.")
    overview_parser.add_argument("--type", action="append", help="Portal object type. May be repeated.")
    overview_parser.add_argument("--limit", type=int, default=25)

    smoke_parser = subparsers.add_parser("smoke", help="Run authenticated portal smoke analysis.")
    smoke_parser.add_argument("--type", action="append", help="Portal object type. May be repeated.")
    smoke_parser.add_argument("--limit", type=int, default=25)

    catalog_parser = subparsers.add_parser("catalog-smoke", help="Run public Catalog smoke analysis.")
    catalog_parser.add_argument("--limit", type=int, default=10)

    encode_overview_parser = subparsers.add_parser("encode-overview", help="Sample major ENCODE object types.")
    encode_overview_parser.add_argument("--type", action="append", help="ENCODE object type. May be repeated.")
    encode_overview_parser.add_argument("--limit", type=int, default=25)

    encode_smoke_parser = subparsers.add_parser("encode-smoke", help="Run public ENCODE smoke analysis.")
    encode_smoke_parser.add_argument("--limit", type=int, default=10)

    encode_csv_parser = subparsers.add_parser("encode-export-csv", help="Export ENCODE metadata as CSV.")
    encode_csv_parser.add_argument("--type", default="Experiment")
    encode_csv_parser.add_argument("--limit", type=int, default=100)
    encode_csv_parser.add_argument(
        "--param",
        action="append",
        default=[],
        help="Additional ENCODE search parameter as KEY=VALUE. May be repeated.",
    )

    csv_parser = subparsers.add_parser("export-csv", help="Export a small metadata sample as CSV.")
    csv_parser.add_argument("--type", default="File")
    csv_parser.add_argument("--limit", type=int, default=100)

    args = parser.parse_args(argv)
    setup_logging()

    if args.command == "overview":
        return overview(args)
    if args.command == "smoke":
        return smoke(args)
    if args.command == "catalog-smoke":
        return catalog_smoke(args)
    if args.command == "encode-overview":
        return encode_overview(args)
    if args.command == "encode-smoke":
        return encode_smoke(args)
    if args.command == "encode-export-csv":
        return encode_export_csv(args)
    if args.command == "export-csv":
        return export_csv(args)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
