#!/usr/bin/env python3
"""Refresh front-page IGVF Portal and Knowledge Graph summary statistics."""

from __future__ import annotations

import argparse
import base64
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
SUMMARY_DIR = DATA_DIR / "Summaries"
README = ROOT / "README.md"

PORTAL_BASE = os.environ.get("IGVF_PORTAL_BASE", "https://api.data.igvf.org").rstrip("/")
CATALOG_API_BASE = os.environ.get("IGVF_CATALOG_API_BASE", "https://api.catalogkg.igvf.org").rstrip("/")
ARANGO_BASE = os.environ.get("IGVF_ARANGO_BASE", "https://db.catalog.igvf.org/_db/igvf").rstrip("/")

README_START = "<!-- IGVF_FRONT_PAGE_STATS_START -->"
README_END = "<!-- IGVF_FRONT_PAGE_STATS_END -->"

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

CATALOG_ENDPOINTS = [
    ("Files / filesets", "/api/files-filesets", {"limit": "1", "offset": "0"}),
    ("Gene-element links", "/api/genes/genomic-elements", {"gene_name": "SAMD11", "limit": "1", "page": "0"}),
    ("Region-gene links", "/api/genomic-elements/genes", {"region": "chr1:903900-904900", "limit": "1", "page": "0"}),
    ("Variant-gene links", "/api/variants/genes", {"variant_id": "NC_000001.11:630556:T:C", "limit": "1", "page": "0"}),
    ("Variant predictions", "/api/variants/predictions", {"variant_id": "NC_000001.11:1628997:GGG:GG", "limit": "1", "page": "0"}),
    ("MPRA variant evidence", "/api/variants/biosamples", {"method": "MPRA", "limit": "1", "page": "0"}),
]


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"igvf_frontpage_summary_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    return log_path


def ensure_dirs() -> None:
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    DOCS_DIR.mkdir(parents=True, exist_ok=True)


def build_url(base: str, path: str, params: dict[str, Any] | None = None) -> str:
    url = f"{base}{path if path.startswith('/') else '/' + path}"
    if params:
        url = f"{url}?{urllib.parse.urlencode(params, doseq=True)}"
    return url


def request_json(
    url: str,
    *,
    headers: dict[str, str] | None = None,
    basic_auth: tuple[str, str] | None = None,
    timeout: int = 60,
) -> tuple[int, Any]:
    request_headers = {"Accept": "application/json,*/*", "User-Agent": "IGVFdataAgent/0.1"}
    if headers:
        request_headers.update(headers)
    if basic_auth:
        token = base64.b64encode(f"{basic_auth[0]}:{basic_auth[1]}".encode()).decode()
        request_headers["Authorization"] = f"Basic {token}"
    request = urllib.request.Request(url, headers=request_headers, method="GET")
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


def rows_from_response(data: Any) -> list[dict[str, Any]]:
    if isinstance(data, list):
        return [item for item in data if isinstance(item, dict)]
    if isinstance(data, dict):
        for key in ("@graph", "graph", "results", "result", "data", "items"):
            value = data.get(key)
            if isinstance(value, list):
                return [item for item in value if isinstance(item, dict)]
    return []


def total_from_response(data: Any, rows: list[dict[str, Any]]) -> int:
    if isinstance(data, dict):
        for key in ("total", "@total", "total_results", "count"):
            value = data.get(key)
            if isinstance(value, int):
                return value
    return len(rows)


def portal_headers() -> dict[str, str]:
    headers: dict[str, str] = {}
    if os.environ.get("IGVF_PORTAL_COOKIE"):
        headers["Cookie"] = os.environ["IGVF_PORTAL_COOKIE"]
    return headers


def summarize_portal() -> dict[str, Any]:
    object_summaries = []
    headers = portal_headers()
    for object_type in PORTAL_OBJECT_TYPES:
        params = {"type": object_type, "format": "json", "limit": "1"}
        status, data = request_json(build_url(PORTAL_BASE, "/search/", params), headers=headers)
        rows = rows_from_response(data)
        object_summaries.append(
            {
                "object_type": object_type,
                "http_status": status,
                "total": total_from_response(data, rows),
                "returned_rows": len(rows),
            }
        )
    return {
        "base_url": PORTAL_BASE,
        "used_cookie": bool(headers.get("Cookie")),
        "objects": object_summaries,
        "note": "README front-page stats use robust object-type totals. Data-type-specific discovery is handled by dedicated manifest/search skills.",
    }


def summarize_catalog() -> dict[str, Any]:
    endpoint_summaries = []
    for label, path, params in CATALOG_ENDPOINTS:
        status, data = request_json(build_url(CATALOG_API_BASE, path, params))
        rows = rows_from_response(data)
        endpoint_summaries.append(
            {
                "label": label,
                "path": path,
                "params": params,
                "http_status": status,
                "total": total_from_response(data, rows),
                "returned_rows": len(rows),
            }
        )
    return {"base_url": CATALOG_API_BASE, "endpoints": endpoint_summaries}


def arango_auth() -> tuple[str, str] | None:
    password = os.environ.get("IGVF_ARANGO_PASSWORD")
    if not password:
        return None
    return os.environ.get("IGVF_ARANGO_USER", "guest"), password


def summarize_arango(max_collections: int) -> dict[str, Any]:
    auth = arango_auth()
    if auth is None:
        return {
            "base_url": ARANGO_BASE,
            "queried": False,
            "reason": "IGVF_ARANGO_PASSWORD not set",
            "collections": [],
        }
    status, data = request_json(build_url(ARANGO_BASE, "/_api/collection"), basic_auth=auth)
    rows = []
    if isinstance(data, dict):
        rows = [item for item in data.get("result", []) if isinstance(item, dict)]
    collections = []
    for row in rows[:max_collections]:
        name = row.get("name")
        if not name or str(name).startswith("_"):
            continue
        count_status, count_data = request_json(
            build_url(ARANGO_BASE, f"/_api/collection/{urllib.parse.quote(str(name))}/count"),
            basic_auth=auth,
        )
        count = count_data.get("count") if isinstance(count_data, dict) else None
        collections.append({"name": name, "http_status": count_status, "count": count})
    return {
        "base_url": ARANGO_BASE,
        "queried": True,
        "http_status": status,
        "collection_count": len(rows),
        "collections": collections,
    }


def fmt_int(value: Any) -> str:
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def status_note(items: list[dict[str, Any]]) -> str:
    failures = [item for item in items if not (200 <= int(item.get("http_status") or 0) < 400)]
    if not failures:
        return "all queried endpoints responded successfully"
    return f"{len(failures)} queried endpoint(s) did not return HTTP 2xx/3xx"


def top_items(items: list[dict[str, Any]], key: str = "total", limit: int = 6) -> list[dict[str, Any]]:
    return sorted(items, key=lambda item: int(item.get(key) or 0), reverse=True)[:limit]


def readme_block(summary: dict[str, Any], report_path: Path, json_path: Path) -> str:
    portal = summary["portal"]
    catalog = summary["catalog"]
    arango = summary["knowledge_graph"]
    generated = summary["generated_at"]
    portal_objects = top_items(portal["objects"])
    catalog_endpoints = catalog["endpoints"]
    lines = [
        README_START,
        "## Current IGVF Data Overview",
        "",
        f"Last refreshed: `{generated}`",
        "",
        "This block is generated by:",
        "",
        "```bash",
        "python3 Scripts/igvf_frontpage_summary.py refresh --update-readme",
        "```",
        "",
        "### IGVF Portal Snapshot",
        "",
        f"- Portal API base: `{portal['base_url']}`",
        f"- Authenticated portal cookie used: `{portal['used_cookie']}`",
        f"- Status: {status_note(portal['objects'])}",
        "",
        "| Object type | Reported total | HTTP |",
        "| --- | ---: | ---: |",
    ]
    lines.extend(
        f"| {item['object_type']} | {fmt_int(item['total'])} | {item['http_status']} |"
        for item in portal_objects
    )
    lines.extend(["", f"Note: {portal.get('note', '')}"])
    lines.extend(
        [
            "",
            "### IGVF Catalog / Knowledge Graph Snapshot",
            "",
            f"- Catalog API base: `{catalog['base_url']}`",
            f"- Catalog API status: {status_note(catalog_endpoints)}",
            "",
            "| Catalog evidence class | Returned rows in smoke query | HTTP |",
            "| --- | ---: | ---: |",
        ]
    )
    lines.extend(
        f"| {item['label']} | {fmt_int(item['returned_rows'])} | {item['http_status']} |"
        for item in catalog_endpoints
    )
    lines.extend(["", "Knowledge Graph / ArangoDB:"])
    if arango.get("queried"):
        lines.append(f"- Collections visible: `{fmt_int(arango.get('collection_count'))}`")
        if arango.get("collections"):
            lines.extend(["", "| Collection | Count | HTTP |", "| --- | ---: | ---: |"])
            lines.extend(
                f"| {item['name']} | {fmt_int(item.get('count'))} | {item.get('http_status')} |"
                for item in arango["collections"][:8]
            )
    else:
        lines.append(f"- Not queried in this refresh: {arango.get('reason')}")
    lines.extend(
        [
            "",
            f"Full report: `{display_path(report_path)}`",
            f"Machine-readable summary: `{display_path(json_path)}`",
            README_END,
        ]
    )
    return "\n".join(lines)


def write_report(summary: dict[str, Any], path: Path, json_path: Path) -> None:
    block = readme_block(summary, path, json_path)
    path.write_text(block.replace(README_START + "\n", "").replace("\n" + README_END, "") + "\n", encoding="utf-8")


def update_readme(block: str) -> None:
    text = README.read_text(encoding="utf-8")
    if README_START in text and README_END in text:
        before = text.split(README_START, 1)[0].rstrip()
        after = text.split(README_END, 1)[1].lstrip()
        updated = f"{before}\n\n{block}\n\n{after}"
    else:
        insert_after = "## What It Can Do"
        if insert_after in text:
            head, tail = text.split(insert_after, 1)
            updated = f"{head}{block}\n\n{insert_after}{tail}"
        else:
            updated = f"{block}\n\n{text}"
    README.write_text(updated, encoding="utf-8")


def refresh(args: argparse.Namespace) -> int:
    ensure_dirs()
    summary = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "portal": summarize_portal(),
        "catalog": summarize_catalog(),
        "knowledge_graph": summarize_arango(args.max_collections),
    }
    json_path = SUMMARY_DIR / "igvf_frontpage_summary.json"
    timestamped_json = SUMMARY_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_igvf_frontpage_summary.json"
    report_path = DOCS_DIR / "IGVF_FRONT_PAGE_DATA_SUMMARY.md"
    json_text = json.dumps(summary, indent=2, sort_keys=True)
    json_path.write_text(json_text, encoding="utf-8")
    timestamped_json.write_text(json_text, encoding="utf-8")
    write_report(summary, report_path, json_path)
    block = readme_block(summary, report_path, json_path)
    if args.update_readme:
        update_readme(block)
    print(f"Wrote JSON: {json_path}")
    print(f"Wrote timestamped JSON: {timestamped_json}")
    print(f"Wrote report: {report_path}")
    if args.update_readme:
        print(f"Updated README: {README}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Refresh IGVF front-page summary statistics.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    refresh_parser = subparsers.add_parser("refresh", help="Refresh Portal/Catalog/KG summary stats.")
    refresh_parser.add_argument("--update-readme", action="store_true", help="Replace the generated README stats block.")
    refresh_parser.add_argument("--max-collections", type=int, default=12, help="Arango collections to count when KG credentials are available.")
    args = parser.parse_args(argv)
    setup_logging()
    if args.command == "refresh":
        return refresh(args)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
