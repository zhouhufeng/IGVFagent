#!/usr/bin/env python3
"""cCRE, linkage, FAVOR, and IGVF Catalog annotation skills.

This script keeps large-data operations manifest-first. It can discover SCREEN
cCRE download links, plan or download ENCODE rE2G/single-cell linkage files,
annotate variant lists with FAVOR plus IGVF Catalog evidence, and summarize
local cCRE/linkage files with simple plots.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import html
import json
import logging
import mimetypes
import os
import re
import shutil
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
SKILL_DOC_DIR = DOCS_DIR / "Skills"
CCRE_DATA_DIR = DATA_DIR / "cCRE"
LINKAGE_DATA_DIR = DATA_DIR / "Linkages"
MANIFEST_DIR = DATA_DIR / "Manifests" / "cCRELinkage"
CACHE_DIR = DATA_DIR / "Cache" / "cCRELinkage"
REPORT_DIR = DOCS_DIR / "cCRELinkage"
PLOT_DIR = REPORT_DIR / "Plots"

SCREEN_DOWNLOAD_PAGES = [
    "https://screen.wenglab.org/downloads",
    "https://screen.beta.wenglab.org/downloads",
]
ENCODE_BASE = os.environ.get("ENCODE_BASE", "https://www.encodeproject.org").rstrip("/")
IGVF_PORTAL_BASE = os.environ.get("IGVF_PORTAL_BASE", "https://data.igvf.org").rstrip("/")
CATALOG_API_BASE = os.environ.get(
    "IGVF_CATALOG_API_BASE", "https://api.catalogkg.igvf.org"
).rstrip("/")
FAVOR_API_BASE = os.environ.get("FAVOR_API_BASE", "https://api.genohub.org").rstrip("/")

FALLBACK_SCREEN_DOWNLOADS = [
    {
        "label": "all_human_ccres",
        "category": "cCREs by Class",
        "description": "All Human cCREs, SCREEN V4, GRCh38/hg38",
        "url": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed",
        "expected_count": "2348854",
        "size_mb": "129.1",
    },
    {
        "label": "promoter_like_pls",
        "category": "cCREs by Class",
        "description": "Promoter-like cCREs, PLS",
        "url": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.PLS.bed",
        "expected_count": "47532",
        "size_mb": "2.6",
    },
    {
        "label": "candidate_enhancers_pels_dels",
        "category": "cCREs by Class",
        "description": "All candidate enhancers, pELS and dELS",
        "url": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.ELS.bed",
        "expected_count": "1718669",
        "size_mb": "94.4",
    },
    {
        "label": "proximal_enhancer_like_pels",
        "category": "cCREs by Class",
        "description": "Proximal enhancer-like cCREs, pELS",
        "url": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.pELS.bed",
        "expected_count": "249464",
        "size_mb": "13.7",
    },
    {
        "label": "distal_enhancer_like_dels",
        "category": "cCREs by Class",
        "description": "Distal enhancer-like cCREs, dELS",
        "url": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.dELS.bed",
        "expected_count": "1469205",
        "size_mb": "80.7",
    },
    {
        "label": "ctcf_bound_ccres",
        "category": "cCREs by Class",
        "description": "CTCF-bound cCREs",
        "url": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.CTCF.bed",
        "expected_count": "948642",
        "size_mb": "63.0",
    },
]

CATALOG_VARIANT_ENDPOINTS = {
    "summary": ("/api/variants/summary", {}),
    "qtl_genes": ("/api/variants/genes/summary", {"limit": "20"}),
    "genomic_elements": ("/api/variants/genomic-elements", {"limit": "20"}),
    "predictions": ("/api/variants/predictions", {"limit": "20"}),
    "phenotypes": ("/api/variants/phenotypes", {"limit": "20", "verbose": "false"}),
}

CATALOG_LINKAGE_ENDPOINTS = [
    (
        "region_element_gene_links",
        "/api/genomic-elements/genes",
        {"region": "chr1:903900-904900", "limit": "25", "page": "0"},
    ),
    (
        "gene_element_links",
        "/api/genes/genomic-elements",
        {"gene_name": "SAMD11", "limit": "25", "page": "0", "verbose": "false"},
    ),
    (
        "variant_gene_links",
        "/api/variants/genes",
        {"variant_id": "NC_000001.11:630556:T:C", "limit": "25", "page": "0"},
    ),
    (
        "variant_element_predictions",
        "/api/variants/predictions",
        {"variant_id": "NC_000001.11:1628997:GGG:GG", "limit": "25", "page": "0"},
    ),
]


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"ccre_linkage_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)
    return log_path


def ensure_dirs() -> None:
    for directory in (
        CCRE_DATA_DIR,
        LINKAGE_DATA_DIR,
        MANIFEST_DIR,
        CACHE_DIR,
        REPORT_DIR,
        PLOT_DIR,
        SKILL_DOC_DIR,
    ):
        directory.mkdir(parents=True, exist_ok=True)


def safe_label(value: str) -> str:
    safe = "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in value)
    return safe.strip("_") or "item"


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def request_json(url: str, *, headers: dict[str, str] | None = None, timeout: int = 60) -> tuple[int, Any]:
    req_headers = {
        "Accept": "application/json,*/*",
        "User-Agent": "IGVFdataAgent/0.1",
    }
    if headers:
        req_headers.update(headers)
    request = urllib.request.Request(url, headers=req_headers, method="GET")
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            content = response.read()
            content_type = response.headers.get("Content-Type", "")
            if "json" in content_type:
                return response.status, json.loads(content)
            return response.status, content.decode(errors="replace")
    except urllib.error.HTTPError as exc:
        body = exc.read()
        try:
            return exc.code, json.loads(body)
        except json.JSONDecodeError:
            return exc.code, body.decode(errors="replace")
    except urllib.error.URLError as exc:
        return 0, {"network_error": str(exc.reason), "url": url}


def request_text(url: str, timeout: int = 60) -> tuple[int, str]:
    request = urllib.request.Request(
        url,
        headers={"Accept": "text/html,text/plain,*/*", "User-Agent": "IGVFdataAgent/0.1"},
        method="GET",
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.status, response.read().decode(errors="replace")
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read().decode(errors="replace")
    except urllib.error.URLError as exc:
        return 0, f"Network error: {exc.reason}"


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        field_order: list[str] = []
        for row in rows:
            for key in row:
                if key not in field_order:
                    field_order.append(key)
        fields = field_order or ["empty"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote CSV: %s", path)
    return path


def write_json(path: Path, payload: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    logging.info("Wrote JSON: %s", path)
    return path


def rows_from_response(data: Any) -> list[dict[str, Any]]:
    if isinstance(data, list):
        return [row for row in data if isinstance(row, dict)]
    if isinstance(data, dict):
        for key in ("@graph", "graph", "results", "result", "data", "items"):
            value = data.get(key)
            if isinstance(value, list):
                return [row for row in value if isinstance(row, dict)]
    return []


def scalar(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (str, int, float, bool)):
        return str(value)
    if isinstance(value, list):
        return ";".join(scalar(item) for item in value[:8] if scalar(item))
    if isinstance(value, dict):
        for key in (
            "name",
            "title",
            "accession",
            "gene_name",
            "gene_id",
            "target_gene",
            "method",
            "source",
            "@id",
            "href",
        ):
            if value.get(key):
                return scalar(value[key])
        return json.dumps(value, sort_keys=True)[:240]
    return str(value)


def build_url(base: str, path: str, params: dict[str, Any] | None = None) -> str:
    url = f"{base}{path if path.startswith('/') else '/' + path}"
    if params:
        url = f"{url}?{urllib.parse.urlencode(params, doseq=True)}"
    return url


def infer_size_mb(text: str) -> str:
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(GB|MB|KB)", text, flags=re.I)
    if not match:
        return ""
    value = float(match.group(1))
    unit = match.group(2).upper()
    if unit == "GB":
        value *= 1024
    elif unit == "KB":
        value /= 1024
    return f"{value:.3f}".rstrip("0").rstrip(".")


def parse_screen_downloads(page_url: str, html_text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for match in re.finditer(r'<a[^>]+href="([^"]+)"[^>]*>(.*?)</a>', html_text, flags=re.I | re.S):
        href = html.unescape(match.group(1))
        label_text = re.sub(r"<[^>]+>", " ", match.group(2))
        label_text = html.unescape(re.sub(r"\s+", " ", label_text)).strip()
        if "downloads.wenglab.org" not in href:
            continue
        url = urllib.parse.urljoin(page_url, href)
        lower = label_text.lower()
        if "nearest gene" in lower:
            category = "Nearest Gene by cCRE"
        elif "link" in lower or "association" in lower or "eqtl" in lower or "crispr" in lower:
            category = "cCRE-Gene Links"
        elif "silencer" in lower or "maff" in lower:
            category = "Functional Characterization"
        else:
            category = "cCREs by Class"
        rows.append(
            {
                "label": safe_label(label_text.lower())[:80],
                "category": category,
                "description": label_text,
                "url": url,
                "expected_count": "".join(re.findall(r"\(([0-9,]+)\)", label_text)[:1]).replace(",", ""),
                "size_mb": infer_size_mb(label_text),
                "source_page": page_url,
            }
        )
    return rows


def discover_screen_downloads() -> list[dict[str, str]]:
    discovered: list[dict[str, str]] = []
    for page in SCREEN_DOWNLOAD_PAGES:
        status, text = request_text(page)
        cache = CACHE_DIR / f"{timestamp()}_{safe_label(page)}.html"
        cache.write_text(text, encoding="utf-8")
        logging.info("SCREEN downloads page %s: HTTP %s", page, status)
        if status and 200 <= status < 400:
            discovered.extend(parse_screen_downloads(page, text))
    seen: set[str] = set()
    unique: list[dict[str, str]] = []
    for row in discovered + FALLBACK_SCREEN_DOWNLOADS:
        url = row["url"]
        if url not in seen:
            seen.add(url)
            unique.append(row)
    return unique


def command_screen_manifest(args: argparse.Namespace) -> int:
    rows = discover_screen_downloads()
    if args.category != "all":
        needle = args.category.lower()
        rows = [row for row in rows if needle in row.get("category", "").lower() or needle in row.get("label", "").lower()]
    path = MANIFEST_DIR / f"{timestamp()}_screen_ccre_download_manifest.csv"
    write_csv(path, rows)
    report = REPORT_DIR / f"{timestamp()}_screen_ccre_download_manifest.md"
    lines = [
        "# SCREEN cCRE Download Manifest",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        "Source: https://screen.wenglab.org/downloads",
        "",
        "The manifest is intentionally separate from download execution because full cCRE, cCRE-gene, and matrix files can be large.",
        "",
        f"Manifest CSV: `{path}`",
        "",
        "## Discovered Downloads",
        "",
    ]
    for row in rows:
        size = f", {row.get('size_mb')} MB" if row.get("size_mb") else ""
        lines.append(f"- {row.get('description')} [{row.get('category')}{size}]: {row.get('url')}")
    report.write_text("\n".join(lines), encoding="utf-8")
    print(f"Discovered SCREEN rows: {len(rows)}")
    print(f"Wrote manifest: {path}")
    print(f"Wrote report: {report}")
    return 0


def load_manifest_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def resolve_download_manifest(args: argparse.Namespace) -> Path:
    if args.manifest:
        return Path(args.manifest)
    candidates = sorted(MANIFEST_DIR.glob("*screen_ccre_download_manifest.csv"))
    if not candidates:
        raise FileNotFoundError("No SCREEN manifest found. Run screen-manifest first.")
    return candidates[-1]


def estimate_manifest_gb(rows: list[dict[str, str]]) -> float:
    total_mb = 0.0
    for row in rows:
        try:
            total_mb += float(row.get("size_mb") or 0)
        except ValueError:
            continue
    return total_mb / 1024


def download_url(url: str, dest: Path) -> tuple[str, int]:
    dest.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "IGVFdataAgent/0.1"}, method="GET")
    try:
        with urllib.request.urlopen(request, timeout=120) as response, dest.open("wb") as handle:
            shutil.copyfileobj(response, handle)
            return "downloaded", int(response.headers.get("Content-Length") or dest.stat().st_size)
    except urllib.error.HTTPError as exc:
        return f"http_{exc.code}", 0
    except urllib.error.URLError as exc:
        return f"network_error_{safe_label(str(exc.reason))}", 0


def command_screen_download(args: argparse.Namespace) -> int:
    manifest = resolve_download_manifest(args)
    rows = load_manifest_rows(manifest)
    if args.only:
        needles = [item.strip().lower() for item in args.only.split(",") if item.strip()]
        rows = [
            row
            for row in rows
            if any(needle in (row.get("label", "") + row.get("description", "") + row.get("category", "")).lower() for needle in needles)
        ]
    estimated_gb = estimate_manifest_gb(rows)
    if not args.download:
        print(f"Planned downloads: {len(rows)}")
        print(f"Estimated size from manifest: {estimated_gb:.2f} GB")
        print("Add --download to fetch files.")
        return 0
    if estimated_gb > args.max_download_gb:
        print(f"Refusing download: estimated {estimated_gb:.2f} GB exceeds --max-download-gb {args.max_download_gb}.")
        return 2
    out_rows: list[dict[str, Any]] = []
    for row in rows:
        url = row.get("url", "")
        name = Path(urllib.parse.urlparse(url).path).name or f"{safe_label(row.get('label', 'download'))}.dat"
        dest = CCRE_DATA_DIR / name
        status, bytes_written = download_url(url, dest)
        out_row = dict(row)
        out_row.update({"status": status, "path": str(dest), "bytes": bytes_written})
        out_rows.append(out_row)
        print(f"{status}: {url}")
    path = MANIFEST_DIR / f"{timestamp()}_screen_ccre_download_results.csv"
    write_csv(path, out_rows)
    print(f"Wrote download results: {path}")
    return 0 if all(row["status"] == "downloaded" for row in out_rows) else 1


def encode_search(search_term: str, *, limit: str = "all") -> tuple[int, Any]:
    params = {
        "type": "Annotation",
        "searchTerm": search_term,
        "status!": "archived",
        "limit": limit,
        "format": "json",
    }
    return request_json(build_url(ENCODE_BASE, "/search/", params))


def collect_file_rows_from_encode(rows: list[dict[str, Any]], source_label: str) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for item in rows:
        files = item.get("files") or []
        if isinstance(files, dict):
            files = [files]
        for file_obj in files:
            if not isinstance(file_obj, dict):
                continue
            href = scalar(file_obj.get("href"))
            url = f"{ENCODE_BASE}{href}" if href.startswith("/") else href
            out.append(
                {
                    "source": source_label,
                    "dataset_accession": scalar(item.get("accession")),
                    "dataset_id": scalar(item.get("@id")),
                    "dataset_title": scalar(item.get("description") or item.get("summary") or item.get("title")),
                    "file_accession": scalar(file_obj.get("accession")),
                    "file_id": scalar(file_obj.get("@id")),
                    "file_format": scalar(file_obj.get("file_format")),
                    "output_type": scalar(file_obj.get("output_type")),
                    "content_type": scalar(file_obj.get("content_type")),
                    "assembly": scalar(file_obj.get("assembly")),
                    "href": href,
                    "url": url,
                    "status": scalar(file_obj.get("status")),
                }
            )
    return out


def hydrate_encode_annotation_files(
    rows: list[dict[str, Any]], source_label: str, hydrate_limit: int
) -> tuple[list[dict[str, str]], list[dict[str, Any]]]:
    file_rows: list[dict[str, str]] = []
    hydrated: list[dict[str, Any]] = []
    for index, item in enumerate(rows):
        if hydrate_limit >= 0 and index >= hydrate_limit:
            break
        item_id = scalar(item.get("@id"))
        if not item_id:
            continue
        status, data = request_json(build_url(ENCODE_BASE, item_id, {"format": "json"}))
        hydrated.append({"annotation": item_id, "http_status": status})
        if not (200 <= status < 400) or not isinstance(data, dict):
            continue
        write_json(CACHE_DIR / f"{timestamp()}_{safe_label(source_label)}_{safe_label(item_id)}.json", data)
        file_rows.extend(collect_file_rows_from_encode([data], source_label))
    return file_rows, hydrated


def portal_search(search_term: str, *, limit: str = "100") -> tuple[int, Any]:
    params = {
        "type": "File",
        "searchTerm": search_term,
        "status!": "archived",
        "limit": limit,
        "format": "json",
    }
    headers = {}
    if os.environ.get("IGVF_PORTAL_COOKIE"):
        headers["Cookie"] = os.environ["IGVF_PORTAL_COOKIE"]
    return request_json(build_url(IGVF_PORTAL_BASE, "/search/", params), headers=headers)


def command_linkage_manifest(args: argparse.Namespace) -> int:
    manifests: list[dict[str, str]] = []
    summaries: list[dict[str, Any]] = []
    if args.source in {"encode", "all"}:
        for term, label in (("encode-re2g", "encode_rE2G"), ("single cell linkage", "encode_single_cell_linkage")):
            status, data = encode_search(term, limit=args.limit)
            rows = rows_from_response(data)
            file_rows = collect_file_rows_from_encode(rows, label)
            hydrated: list[dict[str, Any]] = []
            if args.hydrate_limit != 0 and not file_rows:
                hydrate_limit = len(rows) if args.hydrate_limit < 0 else args.hydrate_limit
                file_rows, hydrated = hydrate_encode_annotation_files(rows, label, hydrate_limit)
            summaries.append(
                {
                    "source": label,
                    "http_status": status,
                    "rows": len(rows),
                    "files": len(file_rows),
                    "hydrated": hydrated[:20],
                    "example": rows[:3],
                }
            )
            manifests.extend(file_rows)
            write_json(CACHE_DIR / f"{timestamp()}_{label}_search.json", data)
            print(f"{label}: HTTP {status}, datasets={len(rows)}, files={len(file_rows)}")
    if args.source in {"portal", "all"}:
        for term, label in (("rE2G", "igvf_portal_rE2G"), ("single cell linkage", "igvf_portal_single_cell_linkage")):
            status, data = portal_search(term, limit=args.limit if args.limit != "all" else "100")
            rows = rows_from_response(data)
            summaries.append({"source": label, "http_status": status, "rows": len(rows), "example": rows[:3]})
            for row in rows:
                href = scalar(row.get("href"))
                manifests.append(
                    {
                        "source": label,
                        "dataset_accession": scalar(row.get("dataset") or row.get("accession")),
                        "dataset_id": scalar(row.get("@id")),
                        "dataset_title": scalar(row.get("summary") or row.get("description")),
                        "file_accession": scalar(row.get("accession")),
                        "file_id": scalar(row.get("@id")),
                        "file_format": scalar(row.get("file_format")),
                        "output_type": scalar(row.get("output_type")),
                        "content_type": scalar(row.get("content_type")),
                        "assembly": scalar(row.get("assembly")),
                        "href": href,
                        "url": f"{IGVF_PORTAL_BASE}{href}" if href.startswith("/") else href,
                        "status": scalar(row.get("status")),
                    }
                )
            write_json(CACHE_DIR / f"{timestamp()}_{label}_search.json", data)
            print(f"{label}: HTTP {status}, files={len(rows)}")
    if args.source in {"catalog", "all"}:
        for label, path, params in CATALOG_LINKAGE_ENDPOINTS:
            status, data = request_json(build_url(CATALOG_API_BASE, path, params))
            rows = rows_from_response(data)
            summaries.append({"source": f"catalog_{label}", "http_status": status, "rows": len(rows), "example": rows[:3]})
            write_json(CACHE_DIR / f"{timestamp()}_catalog_{label}.json", data)
            print(f"catalog_{label}: HTTP {status}, rows={len(rows)}")
    manifest_path = MANIFEST_DIR / f"{timestamp()}_{safe_label(args.source)}_linkage_file_manifest.csv"
    write_csv(manifest_path, manifests)
    summary_path = REPORT_DIR / f"{timestamp()}_{safe_label(args.source)}_linkage_manifest_summary.json"
    write_json(summary_path, summaries)
    report_path = REPORT_DIR / f"{timestamp()}_{safe_label(args.source)}_linkage_manifest_report.md"
    lines = [
        "# rE2G And Single-Cell Linkage Manifest",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"Manifest CSV: `{manifest_path}`",
        f"Summary JSON: `{summary_path}`",
        "",
        "Use `linkage-download --manifest ... --download` to fetch the files after reviewing size and file type.",
        "",
        "## Sources Queried",
        "",
    ]
    lines.extend(
        f"- {item['source']}: HTTP {item['http_status']}, rows {item['rows']}, files {item.get('files', 0)}"
        for item in summaries
    )
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote manifest: {manifest_path}")
    print(f"Wrote report: {report_path}")
    return 0


def command_linkage_download(args: argparse.Namespace) -> int:
    rows = load_manifest_rows(Path(args.manifest))
    if args.only:
        needles = [item.strip().lower() for item in args.only.split(",") if item.strip()]
        rows = [row for row in rows if any(needle in json.dumps(row).lower() for needle in needles)]
    if not args.download:
        print(f"Planned linkage downloads: {len(rows)}")
        print("Add --download to fetch files.")
        return 0
    results: list[dict[str, Any]] = []
    for row in rows:
        url = row.get("url", "")
        if not url:
            continue
        name = Path(urllib.parse.urlparse(url).path).name or f"{safe_label(row.get('file_accession') or row.get('source') or 'linkage')}.dat"
        dest = LINKAGE_DATA_DIR / name
        status, bytes_written = download_url(url, dest)
        result = dict(row)
        result.update({"status": status, "path": str(dest), "bytes": bytes_written})
        results.append(result)
        print(f"{status}: {url}")
    out = MANIFEST_DIR / f"{timestamp()}_linkage_download_results.csv"
    write_csv(out, results)
    print(f"Wrote download results: {out}")
    return 0 if all(row.get("status") == "downloaded" for row in results) else 1


def infer_variant_query(row: dict[str, str]) -> tuple[str, str]:
    for key in ("variant_vcf", "favor_variant_vcf", "hg38_vcf"):
        value = (row.get(key) or "").strip()
        if re.match(r"^(chr)?[0-9XYM]+-[0-9]+-[ACGTN]+-[ACGTN]+$", value, flags=re.I):
            return "variant", value.replace("chr", "")
    rsid = (row.get("rsID") or row.get("rsid") or "").strip()
    if rsid.lower().startswith("rs"):
        return "rsid", rsid.lower()
    chrom = (row.get("hg38Chromosome") or row.get("chrom") or row.get("chr") or "").strip()
    pos = (row.get("hg38Position") or row.get("position") or row.get("pos") or "").strip()
    ref = (row.get("hg38Ref") or row.get("ref") or "").strip()
    alt = (row.get("hg38Alt") or row.get("alt") or "").strip()
    if chrom and pos and ref and alt:
        return "variant", f"{chrom.replace('chr', '')}-{pos}-{ref}-{alt}"
    return "", ""


def catalog_variant_params(row: dict[str, str]) -> dict[str, str]:
    spdi = (row.get("SPDI") or row.get("spdi") or "").strip()
    if spdi:
        return {"variant_id": spdi}
    rsid = (row.get("rsID") or row.get("rsid") or "").strip()
    if rsid.lower().startswith("rs"):
        return {"rsid": rsid.lower()}
    kind, value = infer_variant_query(row)
    if kind == "variant":
        chrom, pos, ref, alt = value.split("-", 3)
        return {"region": f"chr{chrom}:{pos}-{pos}"}
    return {}


def fetch_favor_variant(kind: str, value: str) -> tuple[int, Any]:
    if kind == "rsid":
        return request_json(build_url(FAVOR_API_BASE, f"/v1/rsids/{urllib.parse.quote(value)}"))
    if kind == "variant":
        return request_json(build_url(FAVOR_API_BASE, f"/v1/variants/{urllib.parse.quote(value)}"))
    return 0, {"query_error": "missing FAVOR variant identifier"}


def flatten_favor(data: Any) -> dict[str, str]:
    if isinstance(data, list):
        data = data[0] if data else {}
    if not isinstance(data, dict):
        return {}
    keys = [
        "variant_vcf",
        "rsid",
        "chromosome",
        "position",
        "filter_status",
        "bravo_af",
        "genecode_comprehensive_category",
        "genecode_comprehensive_info",
        "genecode_comprehensive_exonic_category",
        "geneinfo",
        "cage_enhancer",
        "cage_promoter",
        "genehancer",
        "super_enhancer",
        "linsight",
        "fathmm_xf",
        "encode_dnase_sum",
        "encodeh3k27ac_sum",
        "funseq_value",
        "funseq_description",
        "clnsig",
        "clndn",
        "clndisdb",
    ]
    return {f"favor_{key}": scalar(data.get(key)) for key in keys}


def fetch_catalog_variant(row: dict[str, str]) -> dict[str, str]:
    params = catalog_variant_params(row)
    out = {"igvf_catalog_query": json.dumps(params, sort_keys=True)}
    if not params:
        out["igvf_catalog_status"] = "missing_query"
        return out
    for name, (path, extras) in CATALOG_VARIANT_ENDPOINTS.items():
        query = dict(params)
        query.update(extras)
        status, data = request_json(build_url(CATALOG_API_BASE, path, query))
        rows = rows_from_response(data)
        out[f"igvf_catalog_{name}_http_status"] = str(status)
        out[f"igvf_catalog_{name}_count"] = str(len(rows) if rows else (1 if isinstance(data, dict) and data else 0))
        out[f"igvf_catalog_{name}_genes"] = ";".join(
            sorted(
                {
                    scalar(item.get("gene") or item.get("target_gene") or item.get("gene_name"))
                    for item in rows
                    if scalar(item.get("gene") or item.get("target_gene") or item.get("gene_name"))
                }
            )[:10]
        )
    return out


def command_annotate_variants(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    label = safe_label(args.label or input_path.stem)
    rows_out: list[dict[str, str]] = []
    evidence_path = DATA_DIR / "Annotated" / "VariantList" / f"{timestamp()}_{label}_favor_igvf_evidence.jsonl"
    evidence_path.parent.mkdir(parents=True, exist_ok=True)
    with input_path.open(newline="", encoding="utf-8-sig") as handle, evidence_path.open("w", encoding="utf-8") as evidence_handle:
        reader = csv.DictReader(handle)
        for index, row in enumerate(reader):
            if args.max_rows is not None and index >= args.max_rows:
                break
            kind, value = infer_variant_query(row)
            favor_status, favor_data = fetch_favor_variant(kind, value)
            catalog_data = fetch_catalog_variant(row)
            annotation = dict(row)
            annotation["favor_query_type"] = kind
            annotation["favor_query"] = value
            annotation["favor_http_status"] = str(favor_status)
            annotation.update(flatten_favor(favor_data))
            annotation.update(catalog_data)
            rows_out.append(annotation)
            evidence_handle.write(
                json.dumps(
                    {
                        "input": row,
                        "favor_query_type": kind,
                        "favor_query": value,
                        "favor_http_status": favor_status,
                        "favor_data": favor_data,
                        "catalog_annotation": catalog_data,
                    },
                    sort_keys=True,
                )
                + "\n"
            )
            if args.sleep:
                time.sleep(args.sleep)
    csv_path = DATA_DIR / "Annotated" / "VariantList" / f"{timestamp()}_{label}_favor_igvf_annotated.csv"
    write_csv(csv_path, rows_out)
    report_path = REPORT_DIR / f"{timestamp()}_{label}_favor_igvf_variant_annotation_report.md"
    report_path.write_text(
        "\n".join(
            [
                "# FAVOR And IGVF Catalog Variant Annotation",
                "",
                f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
                "",
                f"Input: `{input_path}`",
                f"Annotated CSV: `{csv_path}`",
                f"Evidence JSONL: `{evidence_path}`",
                f"Rows annotated: {len(rows_out)}",
                "",
                "FAVOR annotations include coding/noncoding category, nearby genes, regulatory annotations such as CAGE/GeneHancer/super-enhancer fields, ENCODE DNase/H3K27ac summary scores, FunSeq, and ClinVar-style disease fields when present.",
                "IGVF Catalog annotations add variant-gene, variant-genomic-element, predictions, phenotypes, and summary evidence.",
            ]
        ),
        encoding="utf-8",
    )
    print(f"Annotated rows: {len(rows_out)}")
    print(f"Wrote CSV: {csv_path}")
    print(f"Wrote evidence: {evidence_path}")
    print(f"Wrote report: {report_path}")
    return 0


def contains_cosmic(row: dict[str, Any]) -> bool:
    text = json.dumps(row, sort_keys=True).lower()
    return "cosmic" in text


def command_cosmic_from_favor(args: argparse.Namespace) -> int:
    region = args.region.replace("chr", "").replace(":", "-").replace(",", "")
    if "-" not in region:
        print("Use FAVOR region format chr-start-end or chr:start-end.")
        return 2
    endpoint_results: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    status = 0
    data: Any = {}
    used_endpoint = ""
    for suffix in ("", "/snv", "/indel"):
        endpoint = f"/v1/regions/{urllib.parse.quote(region)}{suffix}"
        status, data = request_json(build_url(FAVOR_API_BASE, endpoint, {"limit": args.limit, "page": args.page}))
        candidate_rows = rows_from_response(data)
        endpoint_results.append({"endpoint": endpoint, "http_status": status, "rows": len(candidate_rows)})
        if 200 <= status < 400:
            rows.extend(candidate_rows)
            used_endpoint = endpoint if not used_endpoint else f"{used_endpoint};{endpoint}"
        if suffix == "" and 200 <= status < 400 and candidate_rows:
            break
    cosmic_rows = [row for row in rows if contains_cosmic(row)]
    label = safe_label(args.label or region)
    all_path = DATA_DIR / "FAVOR" / f"{timestamp()}_{label}_favor_region_variants.csv"
    cosmic_path = DATA_DIR / "FAVOR" / f"{timestamp()}_{label}_favor_cosmic_like_variants.csv"
    write_csv(all_path, [{key: scalar(value) for key, value in row.items()} for row in rows])
    write_csv(cosmic_path, [{key: scalar(value) for key, value in row.items()} for row in cosmic_rows])
    report = REPORT_DIR / f"{timestamp()}_{label}_cosmic_from_favor_report.md"
    report.write_text(
        "\n".join(
            [
                "# COSMIC-Like Variant Retrieval From FAVOR",
                "",
                f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
                "",
                f"FAVOR region query: `{region}`",
                f"Endpoint(s) used: `{used_endpoint or 'none succeeded'}`",
                f"Endpoint attempts: `{endpoint_results}`",
                f"Rows returned on queried page(s): {len(rows)}",
                f"Rows with COSMIC text in FAVOR fields: {len(cosmic_rows)}",
                "",
                f"All FAVOR rows: `{all_path}`",
                f"COSMIC-like filtered rows: `{cosmic_path}`",
                "",
                "Note: the public FAVOR documentation describes variant, rsID, region, gene, and batch endpoints. It does not document a dedicated COSMIC-only endpoint, so this skill retrieves FAVOR variants and filters any returned fields containing COSMIC.",
            ]
        ),
        encoding="utf-8",
    )
    print(f"FAVOR endpoint attempts: {endpoint_results}; rows={len(rows)}; cosmic_like={len(cosmic_rows)}")
    print(f"Wrote report: {report}")
    return 0 if rows or any(200 <= item["http_status"] < 400 for item in endpoint_results) else 1


def open_text_auto(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def ccre_class_from_row(row: list[str]) -> str:
    text = "\t".join(row)
    for klass in ("PLS", "pELS", "dELS", "CA-CTCF", "CA-H3K4me3", "CA-TF", "CA", "TF", "CTCF"):
        if klass in text:
            return klass
    return "unknown"


def summarize_ccre_bed(path: Path, limit: int) -> tuple[list[dict[str, str]], collections.Counter[str]]:
    examples: list[dict[str, str]] = []
    counts: collections.Counter[str] = collections.Counter()
    with open_text_auto(path) as handle:
        for index, line in enumerate(handle):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            klass = ccre_class_from_row(parts)
            counts[klass] += 1
            if len(examples) < limit:
                examples.append(
                    {
                        "chrom": parts[0],
                        "start": parts[1],
                        "end": parts[2],
                        "ccre_id": parts[3] if len(parts) > 3 else "",
                        "class": klass,
                    }
                )
    return examples, counts


def write_bar_plot(counter: collections.Counter[str], title: str, x_label: str, path: Path) -> Path:
    items = counter.most_common(20)
    width = 980
    left, right, top, bottom = 260, 100, 72, 72
    row_h = 30
    height = max(280, top + bottom + row_h * max(1, len(items)))
    max_count = max((count for _, count in items), default=1)
    plot_w = width - left - right
    axis_y = height - bottom
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf7"/>',
        f'<text x="{left}" y="32" font-size="21" font-weight="700" font-family="Arial" fill="#1f2933">{html.escape(title)}</text>',
        f'<line x1="{left}" y1="{axis_y}" x2="{width-right}" y2="{axis_y}" stroke="#222"/>',
        f'<line x1="{left}" y1="{top-10}" x2="{left}" y2="{axis_y}" stroke="#222"/>',
        f'<text x="{width/2-60:.0f}" y="{height-24}" font-size="13" font-family="Arial" fill="#1f2933">{html.escape(x_label)}</text>',
    ]
    for tick in range(6):
        value = max_count * tick / 5
        x = left + tick / 5 * plot_w
        parts.append(f'<line x1="{x:.2f}" y1="{top-10}" x2="{x:.2f}" y2="{axis_y}" stroke="#e0e4e8"/>')
        parts.append(f'<text x="{x-18:.2f}" y="{axis_y+20}" font-size="11" font-family="Arial" fill="#47515c">{value:.0f}</text>')
    for index, (name, count) in enumerate(items):
        y = top + index * row_h
        bar_w = count / max_count * plot_w if max_count else 0
        parts.append(f'<text x="22" y="{y+20}" font-size="12" font-family="Arial" fill="#28323c">{html.escape(name[:36])}</text>')
        parts.append(f'<rect x="{left}" y="{y+5}" width="{bar_w:.2f}" height="20" rx="2" fill="#2f776d"/>')
        parts.append(f'<text x="{left + bar_w + 8:.2f}" y="{y+20}" font-size="12" font-family="Arial" fill="#28323c">{count}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")
    return path


def infer_linkage_kind(path: Path, header: list[str]) -> str:
    text = " ".join([path.name] + header).lower()
    if "crispr" in text:
        return "CRISPR"
    if "eqtl" in text or "qtl" in text:
        return "eQTL/QTL"
    if "re2g" in text or "e2g" in text:
        return "ENCODE-rE2G"
    if "single" in text or "scrna" in text or "scatac" in text:
        return "single-cell linkage"
    if "nearest" in text:
        return "nearest-gene"
    return "other"


def summarize_linkage_file(path: Path, sample_rows: int) -> dict[str, Any]:
    with open_text_auto(path) as handle:
        sample = handle.readline()
        delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter=delimiter)
        header = reader.fieldnames or []
        kind = infer_linkage_kind(path, header)
        gene_counter: collections.Counter[str] = collections.Counter()
        ccre_counter: collections.Counter[str] = collections.Counter()
        rows = 0
        for row in reader:
            rows += 1
            for key in ("TargetGene", "target_gene", "gene", "gene_name", "Gene", "symbol"):
                if row.get(key):
                    gene_counter[row[key]] += 1
                    break
            for key in ("ccre", "cCRE", "ccre_id", "cCRE_ID", "name", "element", "Element"):
                if row.get(key):
                    ccre_counter[row[key]] += 1
                    break
            if rows >= sample_rows:
                break
    return {
        "path": str(path),
        "kind": kind,
        "sampled_rows": rows,
        "top_genes": gene_counter.most_common(10),
        "top_ccres": ccre_counter.most_common(10),
    }


def parse_region(region: str) -> tuple[str, int, int]:
    clean = region.replace(",", "").strip()
    match = re.match(r"^(chr)?([0-9XYM]+)[:-]([0-9]+)-([0-9]+)$", clean, flags=re.I)
    if not match:
        raise ValueError("Region must look like chr19:44850000-44910000.")
    chrom = f"chr{match.group(2)}"
    start = int(match.group(3))
    end = int(match.group(4))
    if end <= start:
        raise ValueError("Region end must be greater than start.")
    return chrom, start, end


def overlap(start: int, end: int, region_start: int, region_end: int) -> bool:
    return start < region_end and end > region_start


def load_ccre_browser_records(path: str, chrom: str, start: int, end: int, max_records: int) -> list[dict[str, Any]]:
    if not path:
        return []
    records: list[dict[str, Any]] = []
    with open_text_auto(Path(path)) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] != chrom:
                continue
            try:
                item_start = int(parts[1])
                item_end = int(parts[2])
            except ValueError:
                continue
            if not overlap(item_start, item_end, start, end):
                continue
            records.append(
                {
                    "chrom": parts[0],
                    "start": item_start,
                    "end": item_end,
                    "name": parts[3] if len(parts) > 3 else f"{parts[0]}:{parts[1]}-{parts[2]}",
                    "class": ccre_class_from_row(parts),
                    "source": Path(path).name,
                }
            )
            if len(records) >= max_records:
                break
    return records


def first_int(row: dict[str, str], keys: tuple[str, ...]) -> int | None:
    for key in keys:
        value = (row.get(key) or "").replace(",", "").strip()
        if not value:
            continue
        try:
            return int(float(value))
        except ValueError:
            continue
    return None


def first_text(row: dict[str, str], keys: tuple[str, ...]) -> str:
    for key in keys:
        value = (row.get(key) or "").strip()
        if value:
            return value
    return ""


def load_linkage_browser_records(
    paths: list[str], chrom: str, start: int, end: int, max_records: int
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for pattern in paths:
        candidates = sorted(ROOT.glob(pattern) if not Path(pattern).is_absolute() else Path("/").glob(str(Path(pattern)).lstrip("/")))
        for path in candidates:
            if not path.is_file():
                continue
            with open_text_auto(path) as handle:
                sample = handle.readline()
                delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
                handle.seek(0)
                reader = csv.DictReader(handle, delimiter=delimiter)
                for row in reader:
                    row_chrom = first_text(row, ("chrom", "chr", "Chromosome", "seqnames", "#chrom", "element_chrom"))
                    if row_chrom and not row_chrom.startswith("chr"):
                        row_chrom = f"chr{row_chrom}"
                    if row_chrom and row_chrom != chrom:
                        continue
                    item_start = first_int(row, ("start", "chromStart", "Start", "element_start", "ElementStart"))
                    item_end = first_int(row, ("end", "chromEnd", "End", "element_end", "ElementEnd"))
                    if item_start is None or item_end is None:
                        continue
                    if not overlap(item_start, item_end, start, end):
                        continue
                    gene = first_text(row, ("TargetGene", "target_gene", "gene", "gene_name", "Gene", "symbol"))
                    score = first_text(row, ("score", "Score", "rE2G", "ABC.Score", "correlation", "qvalue", "pvalue"))
                    records.append(
                        {
                            "start": item_start,
                            "end": item_end,
                            "gene": gene or "linked_gene",
                            "score": score,
                            "method": infer_linkage_kind(path, reader.fieldnames or []),
                            "source": path.name,
                        }
                    )
                    if len(records) >= max_records:
                        return records
    return records


def example_browser_data(region: tuple[str, int, int]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    chrom, start, _end = region
    genes = [
        {"name": "TOMM40", "start": 44890500, "end": 44896000, "strand": "+", "tss": 44890500},
        {"name": "APOE", "start": 44905790, "end": 44909390, "strand": "+", "tss": 44905790},
        {"name": "APOC1", "start": 44912000, "end": 44916000, "strand": "-", "tss": 44916000},
    ]
    ccres = [
        {"chrom": chrom, "start": start + 4500, "end": start + 5200, "name": "demo-cCRE-PLS-TOMM40", "class": "PLS", "source": "demo"},
        {"chrom": chrom, "start": start + 17200, "end": start + 18050, "name": "demo-cCRE-pELS-1", "class": "pELS", "source": "demo"},
        {"chrom": chrom, "start": start + 27500, "end": start + 28650, "name": "demo-cCRE-dELS-1", "class": "dELS", "source": "demo"},
        {"chrom": chrom, "start": start + 40600, "end": start + 41450, "name": "demo-cCRE-CTCF", "class": "CTCF", "source": "demo"},
        {"chrom": chrom, "start": start + 55200, "end": start + 56100, "name": "demo-cCRE-PLS-APOE", "class": "PLS", "source": "demo"},
    ]
    links = [
        {"start": ccres[1]["start"], "end": ccres[1]["end"], "gene": "TOMM40", "score": "0.42", "method": "ENCODE-rE2G", "source": "demo"},
        {"start": ccres[2]["start"], "end": ccres[2]["end"], "gene": "APOE", "score": "0.68", "method": "ENCODE-rE2G", "source": "demo"},
        {"start": ccres[3]["start"], "end": ccres[3]["end"], "gene": "APOE", "score": "0.31", "method": "single-cell linkage", "source": "demo"},
        {"start": ccres[4]["start"], "end": ccres[4]["end"], "gene": "APOE", "score": "promoter", "method": "cCRE-gene", "source": "demo"},
    ]
    return ccres, links, genes


def gene_models_for_region(chrom: str, start: int, end: int) -> list[dict[str, Any]]:
    if chrom == "chr19" and overlap(44850000, 44920000, start, end):
        genes = [
            {"name": "TOMM40", "start": 44890500, "end": 44896000, "strand": "+", "tss": 44890500},
            {"name": "APOE", "start": 44905790, "end": 44909390, "strand": "+", "tss": 44905790},
            {"name": "APOC1", "start": 44912000, "end": 44916000, "strand": "-", "tss": 44916000},
        ]
        return [gene for gene in genes if overlap(int(gene["start"]), int(gene["end"]), start, end)]
    return []


def ccre_color(ccre_class: str) -> str:
    palette = {
        "PLS": "#b7342f",
        "pELS": "#d08b18",
        "dELS": "#2f776d",
        "CTCF": "#4f64a3",
        "CA-CTCF": "#4f64a3",
        "CA-H3K4me3": "#7b5aa6",
        "CA-TF": "#5f8d3b",
        "CA": "#667085",
        "unknown": "#88919b",
    }
    return palette.get(ccre_class, "#667085")


def method_color(method: str) -> str:
    lower = method.lower()
    if "re2g" in lower:
        return "#2563a6"
    if "single" in lower:
        return "#8a5a15"
    if "crispr" in lower:
        return "#9b2f5f"
    if "qtl" in lower:
        return "#6b4aa0"
    return "#315f56"


def write_browser_svg(
    region: tuple[str, int, int],
    ccres: list[dict[str, Any]],
    links: list[dict[str, Any]],
    genes: list[dict[str, Any]],
    path: Path,
    title: str,
) -> Path:
    chrom, start, end = region
    width, height = 1240, 620
    left, right = 100, 60
    axis_y = 88
    ccre_y = 170
    link_y = 315
    gene_y = 455
    plot_w = width - left - right

    def x_pos(pos: int) -> float:
        clipped = min(max(pos, start), end)
        return left + (clipped - start) / (end - start) * plot_w

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf7"/>',
        f'<text x="32" y="34" font-size="22" font-weight="700" font-family="Arial" fill="#1f2933">{html.escape(title)}</text>',
        f'<text x="32" y="58" font-size="13" font-family="Arial" fill="#5b6470">{html.escape(chrom)}:{start:,}-{end:,}  |  IGV-like static cCRE plus enhancer-gene linkage view</text>',
        f'<line x1="{left}" y1="{axis_y}" x2="{width-right}" y2="{axis_y}" stroke="#1f2933" stroke-width="1.5"/>',
    ]
    for tick in range(6):
        pos = int(start + (end - start) * tick / 5)
        x = x_pos(pos)
        parts.append(f'<line x1="{x:.2f}" y1="{axis_y-7}" x2="{x:.2f}" y2="{axis_y+7}" stroke="#1f2933"/>')
        parts.append(f'<text x="{x-42:.2f}" y="{axis_y+26}" font-size="11" font-family="Arial" fill="#47515c">{chrom}:{pos:,}</text>')

    parts.extend(
        [
            '<text x="32" y="154" font-size="13" font-weight="700" font-family="Arial" fill="#1f2933">SCREEN cCRE track</text>',
            f'<line x1="{left}" y1="{ccre_y}" x2="{width-right}" y2="{ccre_y}" stroke="#d8dde3"/>',
        ]
    )
    lanes: dict[str, int] = {}
    for record in ccres:
        klass = str(record.get("class") or "unknown")
        lane = lanes.setdefault(klass, len(lanes))
        y = ccre_y + lane * 28
        x1 = x_pos(int(record["start"]))
        x2 = max(x1 + 4, x_pos(int(record["end"])))
        parts.append(f'<rect x="{x1:.2f}" y="{y-10}" width="{x2-x1:.2f}" height="20" rx="2" fill="{ccre_color(klass)}"/>')
        parts.append(f'<text x="{x2+5:.2f}" y="{y+4}" font-size="10" font-family="Arial" fill="#28323c">{html.escape(klass)}</text>')
    if not ccres:
        parts.append(f'<text x="{left}" y="{ccre_y+5}" font-size="12" font-family="Arial" fill="#667085">No cCRE records found for this region.</text>')

    gene_tss = {gene["name"]: int(gene.get("tss") or gene["start"]) for gene in genes}
    parts.extend(
        [
            '<text x="32" y="286" font-size="13" font-weight="700" font-family="Arial" fill="#1f2933">rE2G / cCRE-gene linkage arcs</text>',
            f'<line x1="{left}" y1="{link_y}" x2="{width-right}" y2="{link_y}" stroke="#d8dde3"/>',
        ]
    )
    for index, link in enumerate(links):
        element_center = int((int(link["start"]) + int(link["end"])) / 2)
        target = gene_tss.get(str(link.get("gene")), element_center)
        x1 = x_pos(element_center)
        x2 = x_pos(target)
        arc_h = 45 + 18 * (index % 4)
        color = method_color(str(link.get("method") or "linkage"))
        parts.append(
            f'<path d="M {x1:.2f} {link_y} C {x1:.2f} {link_y-arc_h}, {x2:.2f} {link_y-arc_h}, {x2:.2f} {link_y}" fill="none" stroke="{color}" stroke-width="2.2" opacity="0.82"/>'
        )
        label_x = min(x1, x2) + abs(x2 - x1) / 2 - 35
        label_y = link_y - arc_h - 4
        label = f"{link.get('gene')} {link.get('score')}".strip()
        parts.append(f'<text x="{label_x:.2f}" y="{label_y:.2f}" font-size="10" font-family="Arial" fill="{color}">{html.escape(label[:28])}</text>')
    if not links:
        parts.append(f'<text x="{left}" y="{link_y-10}" font-size="12" font-family="Arial" fill="#667085">No linkage records found for this region.</text>')

    parts.extend(
        [
            '<text x="32" y="438" font-size="13" font-weight="700" font-family="Arial" fill="#1f2933">Gene models</text>',
            f'<line x1="{left}" y1="{gene_y}" x2="{width-right}" y2="{gene_y}" stroke="#d8dde3"/>',
        ]
    )
    for idx, gene in enumerate(genes):
        y = gene_y + idx * 34
        x1 = x_pos(int(gene["start"]))
        x2 = max(x1 + 12, x_pos(int(gene["end"])))
        tss_x = x_pos(int(gene.get("tss") or gene["start"]))
        label_x = x2 + 8 if x2 < width - right - 140 else max(left, x1 - 138)
        parts.append(f'<line x1="{x1:.2f}" y1="{y}" x2="{x2:.2f}" y2="{y}" stroke="#26323f" stroke-width="3"/>')
        parts.append(f'<rect x="{tss_x-3:.2f}" y="{y-8}" width="6" height="16" fill="#26323f"/>')
        parts.append(f'<text x="{label_x:.2f}" y="{y+4}" font-size="12" font-family="Arial" fill="#26323f">{html.escape(gene["name"])} ({html.escape(gene.get("strand", "."))})</text>')
    if not genes:
        parts.append(f'<text x="{left}" y="{gene_y+5}" font-size="12" font-family="Arial" fill="#667085">No built-in gene model available. Add a gene BED/GTF parser in a downstream extension.</text>')

    legend_y = 572
    legend = [("PLS", ccre_color("PLS")), ("pELS", ccre_color("pELS")), ("dELS", ccre_color("dELS")), ("CTCF", ccre_color("CTCF")), ("ENCODE-rE2G", method_color("ENCODE-rE2G")), ("single-cell", method_color("single-cell linkage"))]
    x = 32
    for name, color in legend:
        parts.append(f'<rect x="{x}" y="{legend_y-12}" width="16" height="12" fill="{color}"/>')
        parts.append(f'<text x="{x+22}" y="{legend_y-2}" font-size="11" font-family="Arial" fill="#47515c">{html.escape(name)}</text>')
        x += 118
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")
    return path


def command_browser_demo(args: argparse.Namespace) -> int:
    region = parse_region(args.region)
    chrom, start, end = region
    ccres = load_ccre_browser_records(args.ccre_bed, chrom, start, end, args.max_records)
    links = load_linkage_browser_records(args.linkage_files, chrom, start, end, args.max_records)
    genes = gene_models_for_region(chrom, start, end)
    used_demo_data = False
    if not ccres and not links and args.demo_fallback:
        ccres, links, genes = example_browser_data(region)
        used_demo_data = True
    genes = [gene for gene in genes if overlap(int(gene["start"]), int(gene["end"]), start, end)]
    label = safe_label(args.label or args.region)
    svg_path = PLOT_DIR / f"{timestamp()}_{label}_igv_like_ccre_re2g_browser.svg"
    write_browser_svg(region, ccres, links, genes, svg_path, args.title)
    summary = {
        "region": {"chrom": chrom, "start": start, "end": end},
        "used_demo_data": used_demo_data,
        "ccre_count": len(ccres),
        "linkage_count": len(links),
        "gene_count": len(genes),
        "svg": str(svg_path),
        "ccre_examples": ccres[:20],
        "linkage_examples": links[:20],
    }
    json_path = REPORT_DIR / f"{timestamp()}_{label}_igv_like_browser_summary.json"
    write_json(json_path, summary)
    report = REPORT_DIR / f"{timestamp()}_{label}_igv_like_browser_report.md"
    lines = [
        "# IGV-Like cCRE And rE2G Browser Demo",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"Region: `{chrom}:{start}-{end}`",
        f"SVG browser view: `{svg_path}`",
        f"Summary JSON: `{json_path}`",
        f"Used built-in demo data: `{used_demo_data}`",
        "",
        "## Tracks",
        "",
        f"- SCREEN cCRE intervals: {len(ccres)}",
        f"- rE2G / cCRE-gene linkage arcs: {len(links)}",
        f"- Gene models: {len(genes)}",
        "",
        "## How To Use",
        "",
        "Open the SVG in a browser or Markdown preview. The top track is genomic coordinate position, the middle tracks show cCRE intervals and enhancer-gene arcs, and the bottom track shows gene models/TSS markers.",
        "",
        "For real data, provide `--ccre-bed Data/cCRE/<file.bed>` and one or more `--linkage-files 'Data/Linkages/*'` paths after downloading cCRE and rE2G files.",
    ]
    report.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote browser SVG: {svg_path}")
    print(f"Wrote report: {report}")
    print(f"Wrote summary JSON: {json_path}")
    return 0


def command_summarize_local(args: argparse.Namespace) -> int:
    ccre_examples: list[dict[str, str]] = []
    ccre_counts: collections.Counter[str] = collections.Counter()
    if args.ccre_bed:
        examples, counts = summarize_ccre_bed(Path(args.ccre_bed), args.example_rows)
        ccre_examples = examples
        ccre_counts = counts
    linkage_summaries: list[dict[str, Any]] = []
    linkage_kind_counts: collections.Counter[str] = collections.Counter()
    for pattern in args.linkage_files:
        for file_path in sorted(ROOT.glob(pattern) if not Path(pattern).is_absolute() else Path("/").glob(str(Path(pattern)).lstrip("/"))):
            if file_path.is_file():
                summary = summarize_linkage_file(file_path, args.sample_rows)
                linkage_summaries.append(summary)
                linkage_kind_counts[summary["kind"]] += summary["sampled_rows"]
    label = safe_label(args.label)
    ccre_plot = write_bar_plot(ccre_counts, "cCRE Class Counts", "cCRE count", PLOT_DIR / f"{timestamp()}_{label}_ccre_classes.svg")
    linkage_plot = write_bar_plot(linkage_kind_counts, "Linkage Evidence Rows By Type", "sampled rows", PLOT_DIR / f"{timestamp()}_{label}_linkage_types.svg")
    summary = {
        "ccre_bed": args.ccre_bed,
        "ccre_counts": ccre_counts.most_common(),
        "ccre_examples": ccre_examples,
        "linkage_summaries": linkage_summaries,
        "plots": [str(ccre_plot), str(linkage_plot)],
    }
    json_path = REPORT_DIR / f"{timestamp()}_{label}_ccre_linkage_summary.json"
    write_json(json_path, summary)
    report = REPORT_DIR / f"{timestamp()}_{label}_ccre_linkage_summary.md"
    lines = [
        "# cCRE And Linkage Summary",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"Machine-readable summary: `{json_path}`",
        "",
        "## cCRE Classes",
        "",
    ]
    lines.extend(f"- {name}: {count}" for name, count in ccre_counts.most_common()[:20])
    lines.extend(["", "## Linkage Evidence Types", ""])
    lines.extend(f"- {name}: {count}" for name, count in linkage_kind_counts.most_common())
    lines.extend(["", "## Plots", "", f"- `{ccre_plot}`", f"- `{linkage_plot}`"])
    report.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote report: {report}")
    print(f"Wrote summary JSON: {json_path}")
    return 0


def write_playbook() -> Path:
    path = SKILL_DOC_DIR / "CCRE_LINKAGE_FAVOR_SKILLS.md"
    path.write_text(
        """# Skill: cCRE, Linkage, FAVOR, And COSMIC-Aware Variant Annotation

Use this skill when users need cCRE downloads from SCREEN, rE2G/single-cell linkage files, FAVOR functional annotations, IGVF Catalog regulatory evidence, or COSMIC-like variant filtering through FAVOR outputs.

## Main Commands

```bash
python3 Scripts/ccre_linkage_annotation_skills.py screen-manifest
python3 Scripts/ccre_linkage_annotation_skills.py screen-download --only PLS --download --max-download-gb 1
python3 Scripts/ccre_linkage_annotation_skills.py linkage-manifest --source all --limit all --hydrate-limit -1
python3 Scripts/ccre_linkage_annotation_skills.py linkage-download --manifest Data/Manifests/cCRELinkage/<manifest.csv> --download
python3 Scripts/ccre_linkage_annotation_skills.py annotate-variants --input Data/Input/VariantList/example_variants.csv --max-rows 10
python3 Scripts/ccre_linkage_annotation_skills.py cosmic-from-favor --region chr19:44851820-44908922
python3 Scripts/ccre_linkage_annotation_skills.py summarize-local --ccre-bed Data/cCRE/GRCh38-cCREs.bed --linkage-files 'Data/Linkages/*'
python3 Scripts/ccre_linkage_annotation_skills.py browser-demo --region chr19:44850000-44910000
```

## Workflow

1. Run `screen-manifest` first. Review sizes and categories before downloading.
2. Download cCRE class files or cCRE-gene association files with `screen-download`.
3. Run `linkage-manifest --source all --limit all --hydrate-limit -1` to discover ENCODE-rE2G, single-cell linkage, IGVF Portal, and IGVF Catalog linkage evidence. Use a smaller `--hydrate-limit` for smoke tests.
4. Download selected linkage files only after checking the manifest.
5. Annotate user variants with `annotate-variants`, which joins FAVOR fields and IGVF Catalog evidence.
6. Use `cosmic-from-favor` for region-level FAVOR retrieval and COSMIC text filtering. FAVOR docs do not expose a COSMIC-only endpoint, so the skill filters returned variant records.
7. Use `summarize-local` to make cCRE class and linkage evidence plots from downloaded files.
8. Use `browser-demo` to generate an IGV-like static SVG with coordinate axis, cCRE tracks, rE2G/cCRE-gene arcs, and gene models.

## Interpretation

- cCRE classes answer what regulatory element class a variant or linkage lands in.
- rE2G links provide computational enhancer-to-gene predictions across ENCODE biosamples.
- Single-cell linkage files support cell-type-specific enhancer-gene interpretation when available.
- IGVF Catalog edges add experimental, QTL, and prediction evidence around variants, genes, and genomic elements.
- FAVOR adds broad functional annotation, including coding class, nearby genes, population frequency, epigenetic scores, GeneHancer/CAGE/super-enhancer fields, and disease fields.
- IGV-like browser SVGs are for static reporting and quick interpretation; use IGV/JBrowse/UCSC for interactive inspection of very large tracks.
""",
        encoding="utf-8",
    )
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="cCRE/linkage/FAVOR/IGVF Catalog skill pack.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    screen_manifest = subparsers.add_parser("screen-manifest", help="Discover SCREEN cCRE download URLs.")
    screen_manifest.add_argument("--category", default="all", help="Filter by category or label substring.")

    screen_download = subparsers.add_parser("screen-download", help="Plan or download SCREEN cCRE files.")
    screen_download.add_argument("--manifest", default="", help="SCREEN manifest CSV. Defaults to latest.")
    screen_download.add_argument("--only", default="", help="Comma-separated label/category substrings to download.")
    screen_download.add_argument("--download", action="store_true", help="Actually download files.")
    screen_download.add_argument("--max-download-gb", type=float, default=1.0)

    linkage_manifest = subparsers.add_parser("linkage-manifest", help="Discover rE2G, single-cell, and Catalog linkage metadata.")
    linkage_manifest.add_argument("--source", choices=["encode", "portal", "catalog", "all"], default="all")
    linkage_manifest.add_argument("--limit", default="100", help="ENCODE/Portal limit, or all.")
    linkage_manifest.add_argument(
        "--hydrate-limit",
        type=int,
        default=50,
        help="Number of ENCODE annotations to hydrate for files. Use -1 for all, 0 to disable.",
    )

    linkage_download = subparsers.add_parser("linkage-download", help="Plan or download linkage files from a manifest.")
    linkage_download.add_argument("--manifest", required=True)
    linkage_download.add_argument("--only", default="", help="Comma-separated substring filter.")
    linkage_download.add_argument("--download", action="store_true")

    annotate = subparsers.add_parser("annotate-variants", help="Annotate variant CSV with FAVOR and IGVF Catalog.")
    annotate.add_argument("--input", required=True)
    annotate.add_argument("--label", default="")
    annotate.add_argument("--max-rows", type=int)
    annotate.add_argument("--sleep", type=float, default=0.05)

    cosmic = subparsers.add_parser("cosmic-from-favor", help="Retrieve FAVOR region variants and filter COSMIC-like fields.")
    cosmic.add_argument("--region", required=True, help="Region such as chr19:44851820-44908922.")
    cosmic.add_argument("--limit", type=int, default=100)
    cosmic.add_argument("--page", type=int, default=1)
    cosmic.add_argument("--label", default="")

    summarize = subparsers.add_parser("summarize-local", help="Summarize local cCRE and linkage files with plots.")
    summarize.add_argument("--ccre-bed", default="")
    summarize.add_argument("--linkage-files", nargs="*", default=[])
    summarize.add_argument("--sample-rows", type=int, default=200000)
    summarize.add_argument("--example-rows", type=int, default=20)
    summarize.add_argument("--label", default="local")

    browser = subparsers.add_parser("browser-demo", help="Create an IGV-like cCRE/rE2G SVG browser view.")
    browser.add_argument("--region", default="chr19:44850000-44910000")
    browser.add_argument("--ccre-bed", default="", help="Optional cCRE BED/BED.GZ file.")
    browser.add_argument("--linkage-files", nargs="*", default=[], help="Optional linkage BED/CSV/TSV files or glob patterns.")
    browser.add_argument("--max-records", type=int, default=80)
    browser.add_argument("--label", default="")
    browser.add_argument("--title", default="cCRE And rE2G Linkage Browser Demo")
    browser.add_argument(
        "--no-demo-fallback",
        dest="demo_fallback",
        action="store_false",
        help="Do not use built-in illustrative records if no local records overlap the region.",
    )
    browser.set_defaults(demo_fallback=True)

    subparsers.add_parser("write-playbook", help="Write skill documentation.")

    args = parser.parse_args(argv)
    ensure_dirs()
    setup_logging()

    if args.command == "screen-manifest":
        return command_screen_manifest(args)
    if args.command == "screen-download":
        return command_screen_download(args)
    if args.command == "linkage-manifest":
        return command_linkage_manifest(args)
    if args.command == "linkage-download":
        return command_linkage_download(args)
    if args.command == "annotate-variants":
        return command_annotate_variants(args)
    if args.command == "cosmic-from-favor":
        return command_cosmic_from_favor(args)
    if args.command == "summarize-local":
        return command_summarize_local(args)
    if args.command == "browser-demo":
        return command_browser_demo(args)
    if args.command == "write-playbook":
        path = write_playbook()
        print(f"Wrote {path}")
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
