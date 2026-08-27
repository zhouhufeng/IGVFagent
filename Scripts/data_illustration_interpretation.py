#!/usr/bin/env python3
"""Explain, illustrate, and download-plan IGVF Portal or ENCODE data IDs/URLs.

Users can paste an IGVF/ENCODE accession, object URL, or search URL. The script
fetches metadata, builds a file inventory, creates small SVG overview plots,
and writes a plain-language report about what the data are and how to use them.
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


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "DataIllustration"
PLOT_DIR = REPORT_DIR / "Plots"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
MANIFEST_DIR = DATA_DIR / "Manifests" / "DataIllustration"
RAW_DIR = DATA_DIR / "Interpreted" / "Metadata"
DOWNLOAD_DIR = DATA_DIR / "Interpreted" / "Downloads"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint

IGVF_API_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_API_BASE")
ENCODE_BASE = _resolve_endpoint("encode", "ENCODE_BASE")


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"data_illustration_interpretation_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)
    return log_path


def mkdirs() -> None:
    for path in (REPORT_DIR, PLOT_DIR, MANIFEST_DIR, RAW_DIR, DOWNLOAD_DIR, SKILL_DOC_DIR):
        path.mkdir(parents=True, exist_ok=True)


def safe_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in label)[:160]


def request_headers(source: str) -> dict[str, str]:
    headers = {"Accept": "application/json,*/*", "User-Agent": "IGVFdataAgent/0.1 data-illustration"}
    if source == "igvf" and os.environ.get("IGVF_PORTAL_COOKIE"):
        headers["Cookie"] = os.environ["IGVF_PORTAL_COOKIE"]
    return headers


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
        for key in ("accession", "@id", "uuid", "title", "name", "term_name", "summary", "href"):
            if key in value:
                return scalar_strings(value[key])
    return []


def text_value(value: Any, limit: int = 8) -> str:
    return "; ".join(scalar_strings(value)[:limit])


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


def build_json_url(target: str) -> tuple[str, str, str]:
    """Return source, normalized URL, and label seed."""
    parsed = urllib.parse.urlparse(target)
    if parsed.scheme and parsed.netloc:
        host = parsed.netloc.lower()
        source = "encode" if "encodeproject.org" in host else "igvf"
        if source == "igvf":
            path = parsed.path if parsed.path else "/"
            base = IGVF_API_BASE
        else:
            path = parsed.path if parsed.path else "/"
            base = ENCODE_BASE
        query_pairs = urllib.parse.parse_qsl(parsed.query, keep_blank_values=True)
        query = dict(query_pairs)
        query["format"] = "json"
        url = f"{base}{path}?{urllib.parse.urlencode(query, doseq=True, quote_via=urllib.parse.quote)}"
        label = safe_label(path.strip("/") or "home")
        return source, url, label

    accession = target.strip().strip("/")
    source = "encode" if accession.startswith(("ENC", "ENCSR", "ENCDO", "ENCA")) else "igvf"
    if source == "igvf":
        if accession.startswith("IGVFDS"):
            path = f"/search/?type=FileSet&accession={urllib.parse.quote(accession)}&format=json"
        elif accession.startswith("IGVFFI"):
            path = f"/search/?type=File&accession={urllib.parse.quote(accession)}&format=json"
        else:
            path = f"/search/?searchTerm={urllib.parse.quote(accession)}&format=json"
        return source, f"{IGVF_API_BASE}{path}", safe_label(accession)
    # Fetch the OBJECT, not a full-text search. `searchTerm=<accession>`
    # returns whatever mentions the accession — for ENCSR000EMT that is 25
    # related files and annotations, and not the experiment itself — so a
    # summary built from those hits describes the wrong thing. ENCODE serves a
    # bare accession as a 301 to its typed path, and urlopen follows redirects,
    # so the object endpoint resolves directly.
    path = f"/{urllib.parse.quote(accession)}/?format=json"
    return source, f"{ENCODE_BASE}{path}", safe_label(accession)


def fetch_json_with_fallback(source: str, url: str,
                              accession: str) -> "tuple[int, Any, str]":
    """Object endpoint first, search as fallback.

    Returns (status, data, url_used). An accession that is not a resolvable
    object (a free-text term, a retired id) still needs the search path, so
    the fallback is kept rather than replaced.
    """
    status, data = fetch_json(source, url)
    if 200 <= status < 400 and isinstance(data, dict) and "@graph" not in data:
        return status, data, url
    base = ENCODE_BASE if source == "encode" else IGVF_API_BASE
    alt = f"{base}/search/?searchTerm={urllib.parse.quote(accession)}&format=json"
    if alt == url:
        return status, data, url
    logging.info("Object endpoint returned %s; falling back to search", status)
    s2, d2 = fetch_json(source, alt)
    return s2, d2, alt


def fetch_json(source: str, url: str) -> tuple[int, Any]:
    logging.info("Request: GET %s", url)
    request = urllib.request.Request(url, headers=request_headers(source), method="GET")
    try:
        with urllib.request.urlopen(request, timeout=90) as response:
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


def save_raw(label: str, data: Any) -> Path:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    path = RAW_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}.json"
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    logging.info("Wrote raw metadata: %s", path)
    return path


def hydrate_rows(source: str, rows: list[dict[str, Any]], limit: int) -> list[dict[str, Any]]:
    hydrated: list[dict[str, Any]] = []
    for row in rows[:limit]:
        item_id = text_value(row.get("@id"))
        if not item_id.startswith("/"):
            hydrated.append(row)
            continue
        base = ENCODE_BASE if source == "encode" else IGVF_API_BASE
        url = f"{base}{item_id}?format=json"
        status, detail = fetch_json(source, url)
        if 200 <= status < 400 and isinstance(detail, dict):
            hydrated.append(detail)
        else:
            hydrated.append(row)
    return hydrated


def absolute_download_url(source: str, href: str) -> str:
    if not href:
        return ""
    if href.startswith("http://") or href.startswith("https://"):
        return href
    base = ENCODE_BASE if source == "encode" else IGVF_API_BASE
    return f"{base}{href}"


def is_file_object(row: dict[str, Any]) -> bool:
    types = set(scalar_strings(row.get("@type")))
    return "File" in types or "File" in str(row.get("@id", "")) or bool(row.get("file_format") and row.get("href"))


def extract_files(source: str, data: Any, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    file_rows: list[dict[str, Any]] = []
    candidates = rows if rows else ([data] if isinstance(data, dict) else [])
    for row in candidates:
        if is_file_object(row):
            file_rows.append(row)
        for file_row in row.get("files", []) if isinstance(row, dict) else []:
            if isinstance(file_row, dict):
                file_rows.append(file_row)
    unique: dict[str, dict[str, Any]] = {}
    for file_row in file_rows:
        key = text_value(file_row.get("accession")) or text_value(file_row.get("@id")) or json.dumps(file_row, sort_keys=True)[:80]
        unique[key] = file_row
    return list(unique.values())


def counter_for(rows: list[dict[str, Any]], fields: tuple[str, ...], limit: int = 15) -> list[tuple[str, int]]:
    counter: collections.Counter[str] = collections.Counter()
    for row in rows:
        for field in fields:
            for value in scalar_strings(row.get(field)):
                if value:
                    counter[value] += 1
    return counter.most_common(limit)


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fields = sorted({key for row in rows for key in row}) or ["empty"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote CSV: %s", path)


def manifest_rows(source: str, files: list[dict[str, Any]]) -> list[dict[str, str]]:
    rows = []
    for file_row in files:
        href = text_value(file_row.get("href"))
        size = file_row.get("file_size") or file_row.get("size")
        rows.append(
            {
                "source": source,
                "accession": text_value(file_row.get("accession")),
                "id": text_value(file_row.get("@id") or file_row.get("uuid")),
                "status": text_value(file_row.get("status")),
                "controlled_access": text_value(file_row.get("controlled_access")),
                "file_format": text_value(file_row.get("file_format")),
                "content_type": text_value(file_row.get("content_type") or file_row.get("output_type")),
                "assembly": text_value(file_row.get("assembly")),
                "file_size_bytes": str(size or ""),
                "file_size_gb": f"{float(size or 0) / 1_000_000_000:.3f}" if size else "",
                "md5sum": text_value(file_row.get("md5sum")),
                "href": href,
                "download_url": absolute_download_url(source, href),
                "s3_uri": text_value(file_row.get("s3_uri") or file_row.get("cloud_metadata")),
                "summary": text_value(file_row.get("summary") or file_row.get("description") or file_row.get("title")),
            }
        )
    return rows


def escape_xml(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


def write_bar_plot(items: list[tuple[str, int | float]], title: str, path: Path, x_label: str = "value") -> None:
    width = 1080
    top = 70
    right = 120
    bottom = 78
    left = 310
    row_height = 32
    shown = items[:14]
    height = max(300, top + bottom + row_height * max(len(shown), 1))
    chart = width - left - right
    max_value = max((float(value) for _, value in shown), default=1.0)
    tick_count = 5
    tick_step = max_value / tick_count if max_value else 1.0
    axis_y = height - bottom
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf7"/>',
        f'<text x="24" y="34" font-family="Arial, sans-serif" font-size="21" font-weight="700" fill="#1f2933">{escape_xml(title)}</text>',
        f'<text x="24" y="56" font-family="Arial, sans-serif" font-size="12" fill="#5b6470">Horizontal axis shows count or total size; vertical axis lists metadata categories.</text>',
        f'<line x1="{left}" y1="{axis_y}" x2="{width-right}" y2="{axis_y}" stroke="#1f2933" stroke-width="1"/>',
        f'<line x1="{left}" y1="{top-10}" x2="{left}" y2="{axis_y}" stroke="#1f2933" stroke-width="1"/>',
        f'<text x="{left + chart / 2 - 65:.0f}" y="{height-22}" font-family="Arial, sans-serif" font-size="13" fill="#1f2933">{escape_xml(x_label)}</text>',
        f'<text x="18" y="{top + 10}" font-family="Arial, sans-serif" font-size="13" fill="#1f2933">category</text>',
    ]
    for tick in range(tick_count + 1):
        value = tick * tick_step
        x = left + chart * tick / tick_count
        label = f"{value:.1f}" if max_value < 10 and any(isinstance(v, float) for _, v in shown) else f"{value:.0f}"
        lines.extend(
            [
                f'<line x1="{x:.2f}" y1="{top-10}" x2="{x:.2f}" y2="{axis_y}" stroke="#d9dee2" stroke-width="1"/>',
                f'<line x1="{x:.2f}" y1="{axis_y}" x2="{x:.2f}" y2="{axis_y+5}" stroke="#1f2933" stroke-width="1"/>',
                f'<text x="{x-12:.2f}" y="{axis_y+22}" font-family="Arial, sans-serif" font-size="11" fill="#47515c">{escape_xml(label)}</text>',
            ]
        )
    if not shown:
        lines.append(f'<text x="{left}" y="{top+30}" font-family="Arial, sans-serif" font-size="14" fill="#5b6470">No data available for this plot.</text>')
    for index, (name, value) in enumerate(shown):
        y = top + index * row_height
        bar_width = int(chart * float(value) / max_value) if max_value else 0
        value_label = f"{float(value):.3f}" if isinstance(value, float) and float(value) < 10 else f"{value}"
        lines.append(f'<text x="24" y="{y + 20}" font-family="Arial, sans-serif" font-size="12" fill="#28323c">{escape_xml(str(name)[:48])}</text>')
        lines.append(f'<rect x="{left}" y="{y + 5}" width="{bar_width}" height="20" rx="2" fill="#2f776d"/>')
        lines.append(f'<text x="{min(left + bar_width + 8, width-right+8)}" y="{y + 20}" font-family="Arial, sans-serif" font-size="12" fill="#28323c">{escape_xml(value_label)}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")
    logging.info("Wrote plot: %s", path)


def sum_file_size_by(files: list[dict[str, Any]], fields: tuple[str, ...]) -> list[tuple[str, float]]:
    totals: collections.Counter[str] = collections.Counter()
    for file_row in files:
        names: list[str] = []
        for field in fields:
            names = scalar_strings(file_row.get(field))
            if names:
                break
        if not names:
            names = ["unknown"]
        size = float(file_row.get("file_size") or file_row.get("size") or 0) / 1_000_000_000
        for name in names[:1]:
            totals[name] += size
    return [(name, round(value, 3)) for name, value in totals.most_common(15)]


def infer_research_uses(source: str, data: Any, rows: list[dict[str, Any]], files: list[dict[str, Any]]) -> list[str]:
    text = " ".join(
        scalar_strings(data)
        + [text_value(row.get("summary") or row.get("description") or row.get("title")) for row in rows[:20]]
        + [text_value(file_row.get("content_type") or file_row.get("output_type") or file_row.get("summary")) for file_row in files[:50]]
    ).lower()
    uses = []
    if "re2g" in text or "element to gene" in text or "regulatory element" in text:
        uses.append("Use as enhancer-to-gene or regulatory-element-to-gene evidence for variant-to-gene nomination, GWAS locus interpretation, and target-gene prioritization.")
    if "curated" in text:
        uses.append("Use as a curated entry point: first inventory the linked datasets/files, then choose the assay-specific processed files rather than raw reads when possible.")
    if "fragments" in text or "atac" in text or "chromatin" in text:
        uses.append("Use for chromatin accessibility analysis: QC fragments, call or reuse peaks, cluster cell states, and intersect variants or cCREs with accessible regions.")
    if "gene count" in text or "rna" in text or "expression" in text:
        uses.append("Use for expression analysis: load count matrices into AnnData/Scanpy, normalize, annotate cell types, and test candidate genes across tissues or conditions.")
    if "binding effect" in text or "variant binding" in text:
        uses.append("Use for regulatory variant interpretation: rank variants by predicted TF binding gain/loss and integrate with MPRA, CRISPRi, eQTL, and accessibility evidence.")
    if "sge" in text or "saturation genome editing" in text or "variant effect" in text:
        uses.append("Use as experimental variant-function evidence: standardize alleles and score direction, then combine with gene, disease, and regulatory annotations.")
    if "annotation" in text and source == "encode":
        uses.append("Use as reference annotation; check the biosample, genome assembly, annotation method, and file format before joining to IGVF variant or gene evidence.")
    if not uses:
        uses.append("Use the report as a triage map: inspect file formats, content types, sample metadata, and status before downloading large payloads.")
    return uses


def plain_summary(source: str, data: Any, rows: list[dict[str, Any]], files: list[dict[str, Any]], total: int) -> list[str]:
    obj = data if isinstance(data, dict) else {}
    types = counter_for(rows or ([obj] if obj else []), ("@type",), limit=8)
    assays = counter_for(rows or ([obj] if obj else []), ("assay_title", "assay_titles", "preferred_assay_titles", "assay_term_name"), limit=8)
    file_formats = counter_for(files, ("file_format",), limit=8)
    content_types = counter_for(files, ("content_type", "output_type"), limit=8)
    lines = [
        f"Source: `{source}`",
        f"Returned/represented items: `{len(rows) if rows else 1}`; reported total: `{total}`",
        f"Files discovered directly in metadata: `{len(files)}`",
        f"Primary summary: {text_value(obj.get('summary') or obj.get('description') or obj.get('title')) or 'not provided in top-level metadata'}",
        "",
        "Observed object types:",
        format_counter(types),
        "",
        "Observed assays:",
        format_counter(assays),
        "",
        "File formats:",
        format_counter(file_formats),
        "",
        "File/content types:",
        format_counter(content_types),
    ]
    return lines


def format_counter(items: list[tuple[str, int | float]]) -> str:
    if not items:
        return "- none observed"
    return "\n".join(f"- {name}: {value}" for name, value in items)


def write_report(
    source: str,
    target: str,
    url: str,
    status: int,
    data: Any,
    raw_path: Path,
    manifest_path: Path,
    plots: list[Path],
    rows: list[dict[str, Any]],
    files: list[dict[str, Any]],
    label: str,
) -> Path:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    total = total_from_response(data, rows)
    uses = infer_research_uses(source, data, rows, files)
    path = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_interpretation_report.md"
    lines = [
        "# Data Illustration And Interpretation Report",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"Input: `{target}`",
        f"Resolved JSON URL: `{url}`",
        f"HTTP status: `{status}`",
        f"Raw metadata JSON: `{raw_path}`",
        f"File/download manifest: `{manifest_path}`",
        "",
        "## What This Data Appears To Be",
        "",
        *plain_summary(source, data, rows, files, total),
        "",
        "## How To Use It",
        "",
        *[f"- {item}" for item in uses],
        "",
        "## Practical Next Steps",
        "",
        "- Start with processed files when available: matrices, BED/TSV annotations, bigWig/bedGraph signal, or variant score tables.",
        "- Confirm genome assembly, biosample/tissue/cell type, assay, status, and controlled-access fields before analysis.",
        "- Keep raw reads for reprocessing only when the processed outputs do not answer the question.",
        "- For variant/gene interpretation, join this evidence to IGVF Catalog/KG genes, variants, enhancer-gene links, MPRA, CRISPRi, eQTL, and ENCODE references.",
        "",
        "## Plots",
        "",
        *[f"- `{plot}`" for plot in plots],
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    logging.info("Wrote report: %s", path)
    return path


def download_files(source: str, manifest: list[dict[str, str]], label: str, max_download_gb: float) -> Path:
    destination = DOWNLOAD_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    destination.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    spent = 0.0
    for row in manifest:
        url = row.get("download_url", "")
        if not url:
            continue
        size_gb = float(row.get("file_size_gb") or 0)
        out = dict(row)
        if spent + size_gb > max_download_gb:
            out["download_status"] = "skipped_size_cap"
            out["local_path"] = ""
            rows.append(out)
            continue
        filename = Path(urllib.parse.urlparse(url).path).name or row.get("accession") or "downloaded_file"
        local_path = destination / safe_label(filename)
        logging.info("Downloading %s to %s", url, local_path)
        request = urllib.request.Request(url, headers=request_headers(source), method="GET")
        try:
            with urllib.request.urlopen(request, timeout=180) as response, local_path.open("wb") as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
            spent += size_gb
            out["download_status"] = "downloaded"
            out["local_path"] = str(local_path)
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            out["download_status"] = f"download_error: {exc}"
            out["local_path"] = str(local_path)
        rows.append(out)
    path = MANIFEST_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_download_results.csv"
    write_csv(path, rows)
    return path


def run_explain(args: argparse.Namespace) -> int:
    mkdirs()
    source, url, label_seed = build_json_url(args.target)
    label = args.label or label_seed
    status, data, url = fetch_json_with_fallback(source, url, args.target.strip().strip("/"))
    raw_path = save_raw(label, data)
    rows = rows_from_response(data)
    hydrated_rows = hydrate_rows(source, rows, args.hydrate_limit) if rows and args.hydrate_limit else []
    hydrated_path = None
    if hydrated_rows:
        hydrated_path = save_raw(f"{label}_hydrated_items", hydrated_rows)
    files = extract_files(source, data, rows)
    if hydrated_rows:
        files.extend(extract_files(source, hydrated_rows, hydrated_rows))
        unique_files: dict[str, dict[str, Any]] = {}
        for file_row in files:
            key = text_value(file_row.get("accession")) or text_value(file_row.get("@id"))
            if key:
                unique_files[key] = file_row
        files = list(unique_files.values())
    manifest = manifest_rows(source, files)
    manifest_path = MANIFEST_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_file_manifest.csv"
    write_csv(manifest_path, manifest)
    plots = [
        PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_file_formats.svg",
        PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_content_types.svg",
        PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_statuses.svg",
        PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_content_type_size_gb.svg",
    ]
    write_bar_plot(counter_for(files, ("file_format",)), "File Formats: Number Of Files", plots[0], "file count")
    write_bar_plot(counter_for(files, ("content_type", "output_type")), "File / Content Types: Number Of Files", plots[1], "file count")
    write_bar_plot(counter_for(rows or files, ("status",)), "Item Statuses: Number Of Items", plots[2], "item count")
    write_bar_plot(sum_file_size_by(files, ("content_type", "output_type")), "Download Size By Content Type (GB)", plots[3], "total GB")
    report_path = write_report(source, args.target, url, status, data, raw_path, manifest_path, plots, hydrated_rows or rows, files, label)
    print(f"Source: {source}")
    print(f"HTTP status: {status}")
    print(f"Rows/items: {len(rows) if rows else 1}")
    print(f"Files discovered: {len(files)}")
    # Structured identity BEFORE the paths, so the model summarises what was
    # actually retrieved instead of guessing from the accession alone.
    # rows_from_response only unwraps @graph-style envelopes, so a direct
    # object fetch yields no rows — fall back to the object itself, which is
    # precisely the case that carries the identifying metadata.
    identity_rows = hydrated_rows or rows
    if not identity_rows and isinstance(data, dict) and data.get("accession"):
        identity_rows = [data]
    for line in identity_lines(identity_rows, args.target):
        print(line)
    print(f"Raw metadata: {raw_path}")
    if hydrated_path:
        print(f"Hydrated item metadata: {hydrated_path}")
    print(f"Manifest: {manifest_path}")
    print(f"Report: {report_path}")
    if args.download:
        download_results = download_files(source, manifest, label, args.max_download_gb)
        print(f"Download results: {download_results}")
    return 0 if 200 <= status < 400 else 1


def _scalar(value):
    """ENCODE/IGVF fields are inconsistently scalar, list, or linked object."""
    if isinstance(value, dict):
        return (value.get("term_name") or value.get("label")
                or value.get("title") or value.get("name"))
    if isinstance(value, list):
        vals = [_scalar(v) for v in value]
        vals = [v for v in vals if v]
        return "; ".join(dict.fromkeys(str(v) for v in vals)) or None
    return value


def identity_lines(rows, target: str) -> "list[str]":
    """Human- and model-readable identity of what was actually retrieved.

    Without this the tool printed only counts and file paths, so a model asked
    what an experiment *is* had nothing but the accession to go on and would
    confabulate a plausible assay and cell line. (Observed: ENCSR000EMT, a
    DNase-seq experiment in GM12878, summarised as RNA-seq in K562.) The
    identifying fields are already in hand here; they simply were never
    surfaced. Printing them is what makes the final answer checkable against
    the source rather than invented.
    """
    if not rows:
        return []
    fields = [("assay_term_name", "Assay"), ("assay_title", "Assay title"),
              ("biosample_ontology", "Biosample"), ("target", "Target"),
              ("lab", "Lab"), ("status", "Status"), ("assembly", "Assembly")]

    # The queried accession is not necessarily the first row: an accession
    # search also returns related annotations and derived items.
    t = (target or "").strip().upper()
    primary = next((r for r in rows
                    if str(r.get("accession", "")).upper() == t), None)

    out = ["Identity:"]
    if primary is not None:
        out.append(f"  accession: {primary.get('accession')}  (queried)")
        for key, lbl in fields:
            v = _scalar(primary.get(key))
            if v:
                out.append(f"  {lbl}: {v}")
        desc = primary.get("description")
        if desc:
            out.append(f"  Description: {str(desc)[:200]}")
    else:
        out.append(f"  NOTE: no returned record has accession {target!r}; "
                   f"the {len(rows)} item(s) below are related records.")

    if len(rows) > 1:
        out.append(f"  Retrieved {len(rows)} items in total. Distribution:")
        for key, lbl in (("assay_term_name", "Assay"),
                          ("biosample_ontology", "Biosample"),
                          ("@type", "Type")):
            counts = {}
            for r in rows:
                v = _scalar(r.get(key)) or "(none)"
                counts[str(v)] = counts.get(str(v), 0) + 1
            if len(counts) > 1 or "(none)" not in counts:
                top = sorted(counts.items(), key=lambda kv: -kv[1])[:5]
                out.append("    " + lbl + ": "
                           + ", ".join(f"{k} ({n})" for k, n in top))
    return out


def write_playbook() -> Path:
    mkdirs()
    path = SKILL_DOC_DIR / "DATA_ILLUSTRATION_INTERPRETATION_SKILLS.md"
    lines = [
        "# Skill: Data Illustration And Interpretation For IGVF/ENCODE",
        "",
        "Use this skill when a user pastes an IGVF Portal or ENCODE accession, object URL, or search URL and needs to understand what the data actually are, which files matter, how to download them, and how to apply them in research.",
        "",
        "## Commands",
        "",
        "```bash",
        "python3 Scripts/data_illustration_interpretation.py explain '<igvf-portal-url>/curated-sets/IGVFDS2544COZH/'",
        "python3 Scripts/data_illustration_interpretation.py explain '<encode-portal-url>/search/?type=Annotation&searchTerm=encode-re2g&status!=archived'",
        "python3 Scripts/data_illustration_interpretation.py explain IGVFDS2544COZH --download --max-download-gb 2",
        "```",
        "",
        "## Workflow",
        "",
        "1. Resolve the pasted ID or URL to JSON metadata using the configured IGVF API base or ENCODE `format=json`.",
        "2. Preserve raw JSON under `Data/Interpreted/Metadata/`.",
        "3. Extract directly linked files and write a download manifest under `Data/Manifests/DataIllustration/`.",
        "4. Generate SVG plots for file formats, content types, and statuses under `Docs/DataIllustration/Plots/`.",
        "5. Write a plain-language report explaining what the object/search appears to represent and how to use it.",
        "6. Only download payloads when explicitly requested with `--download` and a size cap.",
        "",
        "## Interpretation Rules",
        "",
        "- Translate vague metadata into concrete analysis questions: expression, accessibility, enhancer-gene linkage, variant effect, binding effect, SGE, annotation, or curated collection.",
        "- Prefer processed files for first-pass research use; raw reads are for reprocessing or method development.",
        "- Always check assembly, biosample, assay, status, controlled access, file size, and provenance before joining data sources.",
        "- Suggest IGVF Catalog/KG and ENCODE integrations when the data can support variant-to-gene, enhancer-to-gene, regulatory element, or disease-locus interpretation.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    logging.info("Wrote playbook: %s", path)
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Explain and illustrate IGVF/ENCODE data IDs or URLs.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    explain = subparsers.add_parser("explain", help="Fetch, summarize, plot, and optionally download an IGVF/ENCODE target.")
    explain.add_argument("target", help="IGVF/ENCODE accession, object URL, or search URL.")
    explain.add_argument("--label", default="", help="Output label.")
    explain.add_argument("--download", action="store_true", help="Download files found in the manifest.")
    explain.add_argument("--max-download-gb", type=float, default=2.0, help="Maximum total payload download size.")
    explain.add_argument("--hydrate-limit", type=int, default=25, help="Fetch detail JSON for this many search-result rows.")

    subparsers.add_parser("write-playbook", help="Write the data illustration skill document.")

    args = parser.parse_args(argv)
    setup_logging()
    if args.command == "write-playbook":
        path = write_playbook()
        print(f"Wrote {path}")
        return 0
    if args.command == "explain":
        return run_explain(args)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
