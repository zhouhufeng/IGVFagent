#!/usr/bin/env python3
"""Starter client for IGVF Portal, Catalog docs, and Knowledge Graph access."""

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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
LOG_DIR = ROOT / "Docs" / "Logs"

PORTAL_BASE = _resolve_endpoint("portal", "IGVF_PORTAL_BASE")
# Data plane. Serves search/object queries anonymously, unlike the
# browser host above.
PORTAL_API_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_API_BASE")
CATALOG_DOCS_BASE = _resolve_endpoint("catalog_docs", "IGVF_CATALOG_DOCS_BASE")
CATALOG_API_BASE = _resolve_endpoint("catalog_api", "IGVF_CATALOG_API_BASE")
ARANGO_BASE = _resolve_endpoint("arango", "IGVF_ARANGO_BASE")
ENCODE_BASE = _resolve_endpoint("encode", "ENCODE_BASE")


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"igvf_client_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)
    return log_path


def output_path(label: str, suffix: str = "json") -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    safe = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in label)
    return DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe}.{suffix}"


def build_request(
    url: str,
    *,
    method: str = "GET",
    body: bytes | None = None,
    headers: dict[str, str] | None = None,
    basic_auth: tuple[str, str] | None = None,
) -> urllib.request.Request:
    req_headers = {
        "Accept": "application/json,text/plain,*/*",
        "User-Agent": "IGVFdataAgent/0.1",
    }
    if headers:
        req_headers.update(headers)
    if basic_auth:
        token = base64.b64encode(f"{basic_auth[0]}:{basic_auth[1]}".encode()).decode()
        req_headers["Authorization"] = f"Basic {token}"
    return urllib.request.Request(url, data=body, headers=req_headers, method=method)


def fetch(
    url: str,
    *,
    method: str = "GET",
    payload: dict[str, Any] | None = None,
    headers: dict[str, str] | None = None,
    basic_auth: tuple[str, str] | None = None,
) -> tuple[int, bytes, str]:
    body = None
    req_headers = dict(headers or {})
    if payload is not None:
        body = json.dumps(payload).encode()
        req_headers["Content-Type"] = "application/json"
    request = build_request(
        url, method=method, body=body, headers=req_headers, basic_auth=basic_auth
    )
    logging.info("Request: %s %s", method, url)
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            content = response.read()
            content_type = response.headers.get("Content-Type", "")
            logging.info("Response: %s %s bytes", response.status, len(content))
            return response.status, content, content_type
    except urllib.error.HTTPError as exc:
        content = exc.read()
        logging.warning("HTTP error: %s %s bytes", exc.code, len(content))
        return exc.code, content, exc.headers.get("Content-Type", "")
    except urllib.error.URLError as exc:
        message = f"Network error while fetching {url}: {exc.reason}\n".encode()
        logging.error(message.decode().strip())
        return 0, message, "text/plain"


def save_response(label: str, content: bytes, content_type: str) -> Path:
    suffix = "json" if "json" in content_type else "txt"
    path = output_path(label, suffix)
    path.write_bytes(content)
    logging.info("Saved response: %s", path)
    return path


def check_catalog_docs() -> int:
    status, content, content_type = fetch(f"{CATALOG_DOCS_BASE}/introduction")
    save_response("catalog_introduction", content, content_type)
    print(f"Catalog docs introduction: HTTP {status}")
    return 0 if 200 <= status < 400 else 1


def check_catalog_index() -> int:
    status, content, content_type = fetch(f"{CATALOG_DOCS_BASE}/llms.txt")
    save_response("catalog_llms_index", content, content_type)
    print(f"Catalog docs llms.txt: HTTP {status}")
    return 0 if 200 <= status < 400 else 1


# ---------------------------------------------------------------------------
# Access tiers
# ---------------------------------------------------------------------------
#
# "Public" does not mean "freely accessible". A resource can be public and
# still require a session cookie, an approved data-use agreement, or be under
# embargo. Recording the tier per resource lets error messages, manifests and
# tool descriptions say which of those applies instead of surfacing a bare 403.
#
# Verified by probing each endpoint anonymously; see Docs/THREAT_MODEL.md.
ACCESS_TIERS = {
    "igvf_portal_api": {"tier": "anonymous-api",
                        "host": "api.data.igvf.org",
                        "note": "search, objects and file metadata need no credential"},
    "igvf_portal_web": {"tier": "cookie",
                        "host": "data.igvf.org",
                        "note": "browser host declines anonymous API requests (403); "
                                "not used by any skill",
                        "env": "IGVF_PORTAL_COOKIE"},
    "igvf_catalog":    {"tier": "anonymous-api", "host": "api.catalogkg.igvf.org"},
    "encode":          {"tier": "anonymous-api", "host": "www.encodeproject.org"},
    "favor":           {"tier": "anonymous-api", "host": "api.genohub.org"},
    "geo":             {"tier": "anonymous-api", "host": "eutils.ncbi.nlm.nih.gov"},
    "cellxgene":       {"tier": "anonymous-api", "host": "api.cellxgene.cziscience.com"},
    "zenodo":          {"tier": "anonymous-api", "host": "zenodo.org"},
    "chipatlas":       {"tier": "anonymous-api", "host": "chip-atlas.dbcls.jp"},
    "igvf_arango":     {"tier": "cookie", "note": "guest read-only by default",
                        "env": "IGVF_ARANGO_PASSWORD"},
    "synapse":         {"tier": "dua-gated",
                        "note": "per-study terms; needs an authenticated client"},
    "dbgap":           {"tier": "dua-gated",
                        "note": "controlled access via an approved application"},
    "ega":             {"tier": "dua-gated", "note": "controlled access"},
}


def access_tier(resource: str) -> dict:
    """Machine-readable access tier for a resource key."""
    return ACCESS_TIERS.get(resource, {"tier": "unknown"})


def cmd_auth_check(_args=None) -> int:
    """Probe every resource anonymously and report what actually works.

    Answers "which resources need credentials, and do mine work" with
    measurements rather than assumptions — the question a bare 403 leaves a
    user unable to answer.
    """
    import json as _json
    import urllib.error
    import urllib.request

    probes = [
        ("IGVF Portal (data API)", "igvf_portal_api",
         f"{PORTAL_API_BASE}/search/?type=MeasurementSet&limit=1&format=json"),
        ("IGVF Portal (web host)", "igvf_portal_web", f"{PORTAL_BASE}/"),
        ("IGVF Catalog", "igvf_catalog",
         f"{CATALOG_API_BASE}/api/genes?name=APOE&limit=1"),
        ("ENCODE", "encode", f"{ENCODE_BASE}/ENCSR000EMT/?format=json"),
    ]
    cookie = os.environ.get("IGVF_PORTAL_COOKIE")
    rows = []
    for label, key, url in probes:
        meta = access_tier(key)
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json",
                              **({"Cookie": cookie} if cookie else {})})
            with urllib.request.urlopen(req, timeout=30) as r:
                status = r.status
        except urllib.error.HTTPError as e:
            status = e.code
        except Exception:
            status = 0
        ok = 200 <= status < 400
        rows.append({"resource": label, "key": key, "tier": meta["tier"],
                     "status": status, "reachable": ok,
                     "note": meta.get("note", "")})
        mark = "ok " if ok else "FAIL"
        print(f"  [{mark}] {label:26s} tier={meta['tier']:14s} HTTP {status}")
        if not ok and meta.get("env"):
            print(f"         needs {meta['env']} — export a session cookie "
                  f"from a logged-in browser")
        elif not ok and meta["tier"] == "dua-gated":
            print(f"         controlled access: requires an approved "
                  f"data-use agreement")

    api_ok = next((r["reachable"] for r in rows
                   if r["key"] == "igvf_portal_api"), False)
    web_ok = next((r["reachable"] for r in rows
                   if r["key"] == "igvf_portal_web"), False)
    print()
    if api_ok and not web_ok:
        print("  IGVF Portal: anonymous data access WORKS. The web host's 403 "
              "is expected and affects nothing — every Portal skill uses the "
              "data API.")
    elif not api_ok:
        print("  IGVF Portal data API unreachable. Set IGVF_PORTAL_COOKIE, or "
              "check network egress.")
    print(f"\n  Cookie present: {'yes' if cookie else 'no'} "
          f"(needed only for unreleased or restricted records)")
    # save_response takes bytes, not str.
    save_response("auth_check", _json.dumps(rows, indent=2).encode(), "application/json")
    return 0 if api_ok else 1


def check_portal() -> int:
    """Report Portal reachability for BOTH hosts, separately.

    The Portal is served from two hosts with different auth behaviour, and
    conflating them is misleading: the browser host (data.igvf.org) refuses
    anonymous requests with 403, while the data API (api.data.igvf.org) serves
    search and object queries anonymously. Probing only the browser host and
    concluding "authenticated access is likely required" tells a user their
    setup is broken when every data-access skill in fact works — which is what
    happened in review.
    """
    session_cookie = os.environ.get("IGVF_PORTAL_COOKIE")
    headers = {"Cookie": session_cookie} if session_cookie else {}

    web_status, content, content_type = fetch(f"{PORTAL_BASE}/", headers=headers)
    save_response("portal_home", content, content_type)
    print(f"IGVF Portal (web, {PORTAL_BASE}): HTTP {web_status}")

    api_url = f"{PORTAL_API_BASE}/search/?type=MeasurementSet&limit=1&format=json"
    api_status, api_content, api_ct = fetch(api_url, headers=headers)
    save_response("portal_api_probe", api_content, api_ct)
    print(f"IGVF Portal (data API, {PORTAL_API_BASE}): HTTP {api_status}")

    if 200 <= api_status < 400:
        print("Anonymous data access works: search, object and file metadata "
              "queries need no cookie. Every Portal skill uses this host.")
        if web_status == 403:
            print(f"The {web_status} above is the browser host only, which "
                  "declines anonymous requests. It does not affect analysis.")
        return 0

    if not session_cookie:
        print("The data API is not reachable anonymously from here. Set "
              "IGVF_PORTAL_COOKIE for authenticated access, or check network "
              "egress to the Portal.")
    return 1


def catalog_api_url(path: str, params: dict[str, Any] | None = None) -> str:
    normalized = path if path.startswith("/") else f"/{path}"
    url = f"{CATALOG_API_BASE}{normalized}"
    if params:
        url = f"{url}?{urllib.parse.urlencode(params, doseq=True)}"
    return url


def fetch_catalog_api(path: str, label: str, params: dict[str, Any] | None = None) -> int:
    url = catalog_api_url(path, params)
    status, content, content_type = fetch(url)
    save_response(label, content, content_type)
    print(f"Catalog API {path}: HTTP {status}")
    return 0 if 200 <= status < 400 else 1


def fetch_catalog_files(limit: int, offset: int) -> int:
    params = {"limit": limit, "offset": offset}
    return fetch_catalog_api("/api/files-filesets", "catalog_files_filesets", params)


def catalog_get(
    path: str, label: str, params: dict[str, Any] | None = None
) -> tuple[int, Any, Path]:
    url = catalog_api_url(path, params)
    status, content, content_type = fetch(url)
    saved = save_response(label, content, content_type)
    if "json" not in content_type:
        return status, content.decode(errors="replace"), saved
    try:
        return status, json.loads(content), saved
    except json.JSONDecodeError:
        return status, content.decode(errors="replace"), saved


def infer_gene_params(query: str) -> tuple[dict[str, str], dict[str, str]]:
    if query.startswith("ENSG"):
        return {"gene_id": query}, {"gene_id": query}
    if query.startswith("HGNC:"):
        return {"hgnc_id": query}, {"hgnc_id": query}
    if query.startswith("ENTREZ:"):
        return {"entrez": query}, {"gene_name": query}
    return {"name": query}, {"gene_name": query}


def infer_variant_params(query: str) -> dict[str, str]:
    lower = query.lower()
    if lower.startswith("rs"):
        return {"rsid": lower}
    if query.startswith("CA"):
        return {"ca_id": query}
    if query.startswith("NC_") and ":" in query:
        return {"variant_id": query}
    if query.startswith("chr") and ":" in query:
        return {"region": query}
    return {"variant_id": query}


def result_count(value: Any) -> int:
    if isinstance(value, list):
        return len(value)
    if isinstance(value, dict):
        for key in ("results", "result", "data", "items"):
            if isinstance(value.get(key), list):
                return len(value[key])
        return 1 if value else 0
    return 0


def first_items(value: Any, limit: int = 5) -> list[Any]:
    if isinstance(value, list):
        return value[:limit]
    if isinstance(value, dict):
        for key in ("results", "result", "data", "items"):
            if isinstance(value.get(key), list):
                return value[key][:limit]
    return []


def compact_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        if value and isinstance(value[0], dict) and "variant" in value[0]:
            first_variant = value[0].get("variant", {})
            first_scores = value[0].get("scores", [])
            scores = ", ".join(
                f"{score.get('method')}={score.get('score')}"
                for score in first_scores[:3]
                if isinstance(score, dict)
            )
            suffix = f"; {scores}" if scores else ""
            return f"{len(value)} variant(s), first {compact_value(first_variant)}{suffix}"
        return ", ".join(compact_value(item) for item in value[:4])
    if isinstance(value, dict):
        if "coding_variant_id" in value:
            change = value.get("hgvsp") or value.get("coding_variant_id")
            protein = value.get("protein_name") or value.get("protein_id")
            return " ".join(str(part) for part in (change, protein) if part)
        if "variant" in value and isinstance(value.get("variant"), dict):
            return compact_value(value["variant"])
        if "scores" in value and isinstance(value.get("scores"), list):
            return ", ".join(
                f"{score.get('method')}={score.get('score')}"
                for score in value["scores"][:3]
                if isinstance(score, dict)
            )
        for key in ("name", "gene_name", "gene_id", "_id", "spdi", "hgvs", "label"):
            if value.get(key):
                return str(value[key])
        return json.dumps(value, sort_keys=True)[:160]
    return str(value)


def print_gene_records(records: Any) -> None:
    items = first_items(records)
    if not items:
        print("No gene records returned.")
        return
    print("Gene records:")
    for item in items:
        if not isinstance(item, dict):
            print(f"- {compact_value(item)}")
            continue
        coords = ""
        if item.get("chr") and item.get("start") is not None and item.get("end") is not None:
            coords = f" {item['chr']}:{item['start']}-{item['end']}"
        print(
            f"- {item.get('name', item.get('_id', 'unknown'))}"
            f" [{item.get('_id', 'no id')}]{coords}"
            f" {item.get('gene_type', '')}".rstrip()
        )


def print_variant_summary(summary: Any) -> None:
    if not isinstance(summary, dict):
        print(compact_value(summary))
        return
    variant = summary.get("summary", summary)
    if isinstance(variant, dict):
        rsid = compact_value(variant.get("rsid"))
        print("Variant summary:")
        print(f"- SPDI: {variant.get('spdi') or variant.get('variant_id') or ''}")
        print(f"- HGVS: {variant.get('hgvs') or ''}")
        print(f"- rsID: {rsid}")
        print(f"- ClinGen Allele ID: {variant.get('ca_id') or ''}")
        print(f"- Alleles: {variant.get('ref') or ''}>{variant.get('alt') or ''}")
    cadd = summary.get("cadd_scores")
    if isinstance(cadd, dict):
        print(f"- CADD: raw={cadd.get('raw')} phred={cadd.get('phread')}")
    nearest = summary.get("nearest_genes")
    if isinstance(nearest, dict):
        nearest_parts = []
        for key in ("nearestGene", "nearestCodingGene"):
            entry = nearest.get(key)
            if isinstance(entry, dict) and entry.get("gene_name"):
                distance = entry.get("distance")
                distance_text = f", {distance} bp" if distance is not None else ""
                nearest_parts.append(f"{key}={entry['gene_name']} ({entry.get('id')}{distance_text})")
        if nearest_parts:
            print(f"- Nearest genes: {'; '.join(nearest_parts)}")
    elif nearest:
        print(f"- Nearest genes: {compact_value(nearest)}")


def print_association_examples(title: str, records: Any, keys: tuple[str, ...]) -> None:
    print(f"{title}: {result_count(records)} result(s)")
    for item in first_items(records, 5):
        if not isinstance(item, dict):
            print(f"- {compact_value(item)}")
            continue
        parts = [compact_value(item.get(key)) for key in keys if compact_value(item.get(key))]
        if not parts:
            parts = [compact_value(item)]
        print(f"- {' | '.join(parts)}")


def explore_gene(query: str, limit: int) -> int:
    gene_params, related_params = infer_gene_params(query)
    gene_params["limit"] = str(limit)
    related_params["limit"] = str(limit)
    related_params["verbose"] = "false"

    status, genes, gene_path = catalog_get("/api/genes", "gene_lookup", gene_params)
    print(f"Gene lookup: HTTP {status} saved={gene_path}")
    print_gene_records(genes)

    status, coding, coding_path = catalog_get(
        "/api/genes/coding-variants/scores", "gene_coding_variant_scores", related_params
    )
    print(f"Coding variant scores: HTTP {status} saved={coding_path}")
    print_association_examples(
        "Coding variant evidence",
        coding,
        ("protein_change", "variants", "scores"),
    )

    status, elements, elements_path = catalog_get(
        "/api/genes/genomic-elements", "gene_genomic_elements", related_params
    )
    print(f"Gene-genomic element links: HTTP {status} saved={elements_path}")
    print_association_examples(
        "Regulatory element evidence",
        elements,
        ("gene", "genomic_element", "biological_context", "method", "source"),
    )
    return 0 if 200 <= status < 400 else 1


def explore_variant(query: str, limit: int) -> int:
    params = infer_variant_params(query)
    params["limit"] = str(limit)
    params["verbose"] = "false"

    status, summary, summary_path = catalog_get("/api/variants/summary", "variant_summary", params)
    print(f"Variant summary endpoint: HTTP {status} saved={summary_path}")
    print_variant_summary(summary)

    status, genes, genes_path = catalog_get("/api/variants/genes/summary", "variant_gene_summary", params)
    print(f"Variant-gene QTL summary: HTTP {status} saved={genes_path}")
    print_association_examples(
        "Gene modulation evidence",
        genes,
        ("qtl_type", "gene", "effect_size", "log10pvalue", "biological_context", "name"),
    )

    status, phenotypes, phenotypes_path = catalog_get(
        "/api/variants/phenotypes", "variant_phenotypes", params
    )
    print(f"Variant-phenotype links: HTTP {status} saved={phenotypes_path}")
    print_association_examples(
        "Phenotype evidence",
        phenotypes,
        (
            "phenotype_term",
            "phenotype_id",
            "biological_context",
            "score",
            "method",
            "class",
            "source",
        ),
    )

    status, ld, ld_path = catalog_get("/api/variants/variant-ld/summary", "variant_ld_summary", params)
    print(f"Variant LD summary: HTTP {status} saved={ld_path}")
    print_association_examples("LD evidence", ld, ("sequence variant", "r2", "d_prime", "ancestry"))
    return 0 if 200 <= status < 400 else 1


def parse_repeated_params(params: list[str]) -> dict[str, list[str]]:
    parsed: dict[str, list[str]] = {}
    for param in params:
        if "=" not in param:
            raise ValueError(f"Expected KEY=VALUE for --param, got: {param}")
        key, value = param.split("=", 1)
        parsed.setdefault(key, []).append(value)
    return parsed


def encode_search_url(object_type: str, params: list[str]) -> str:
    query: dict[str, list[str]] = {"type": [object_type], "format": ["json"]}
    for key, values in parse_repeated_params(params).items():
        query.setdefault(key, []).extend(values)
    return f"{ENCODE_BASE}/search/?{urllib.parse.urlencode(query, doseq=True)}"


def encode_search(object_type: str, params: list[str]) -> int:
    url = encode_search_url(object_type, params)
    status, content, content_type = fetch(url)
    save_response(f"encode_search_{object_type}", content, content_type)
    print(f"ENCODE search {object_type}: HTTP {status}")
    if status < 400:
        try:
            parsed = json.loads(content)
            total = parsed.get("total", "unknown")
            print(f"Total results: {total}")
        except json.JSONDecodeError:
            print("Response was not JSON.")
    return 0 if 200 <= status < 400 else 1


def arango_auth() -> tuple[str, str] | None:
    user = os.environ.get("IGVF_ARANGO_USER", "guest")
    password = os.environ.get("IGVF_ARANGO_PASSWORD")
    if not password:
        print("Set IGVF_ARANGO_PASSWORD locally before using Knowledge Graph commands.")
        return None
    return user, password


def list_collections() -> int:
    auth = arango_auth()
    if auth is None:
        return 2
    status, content, content_type = fetch(f"{ARANGO_BASE}/_api/collection", basic_auth=auth)
    save_response("arango_collections", content, content_type)
    print(f"Knowledge Graph collections: HTTP {status}")
    if status < 400:
        try:
            parsed = json.loads(content)
            names = [item.get("name") for item in parsed.get("result", []) if item.get("name")]
            for name in sorted(names):
                print(name)
        except json.JSONDecodeError:
            print("Response was not JSON.")
    return 0 if status < 400 else 1


def run_aql(query: str, bind_vars: dict[str, Any] | None = None) -> int:
    auth = arango_auth()
    if auth is None:
        return 2
    payload = {"query": query, "bindVars": bind_vars or {}}
    status, content, content_type = fetch(
        f"{ARANGO_BASE}/_api/cursor", method="POST", payload=payload, basic_auth=auth
    )
    save_response("arango_aql", content, content_type)
    print(f"AQL query: HTTP {status}")
    if status < 400:
        try:
            parsed = json.loads(content)
            print(json.dumps(parsed.get("result", parsed), indent=2)[:4000])
        except json.JSONDecodeError:
            print(content.decode(errors="replace")[:4000])
    return 0 if status < 400 else 1


def fetch_url(url: str, label: str) -> int:
    status, content, content_type = fetch(url)
    save_response(label, content, content_type)
    print(f"{url}: HTTP {status}")
    return 0 if status < 400 else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="IGVF starter client")
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("check", help="Check portal, docs, docs index, and KG collections.")
    subparsers.add_parser("auth-check",
                    help="Probe every resource anonymously and report access tiers")
    subparsers.add_parser("docs", help="Fetch IGVF Catalog introduction.")
    subparsers.add_parser("docs-index", help="Fetch IGVF Catalog llms.txt index.")
    subparsers.add_parser("portal", help="Fetch IGVF Portal home page.")
    catalog_parser = subparsers.add_parser("catalog-api", help="Fetch a Catalog API path.")
    catalog_parser.add_argument("path", help="API path, for example / or /api/files-filesets")
    catalog_parser.add_argument(
        "--param",
        action="append",
        default=[],
        help="Query parameter as KEY=VALUE. May be repeated.",
    )
    files_parser = subparsers.add_parser("catalog-files", help="Fetch Catalog files/filesets table data.")
    files_parser.add_argument("--limit", type=int, default=25)
    files_parser.add_argument("--offset", type=int, default=0)
    gene_parser = subparsers.add_parser("gene", help="Explore and summarize an IGVF Catalog gene.")
    gene_parser.add_argument("query", help="Gene symbol, Ensembl ID, HGNC ID, or other supported gene id.")
    gene_parser.add_argument("--limit", type=int, default=10)
    variant_parser = subparsers.add_parser(
        "variant", help="Explore and summarize an IGVF Catalog variant."
    )
    variant_parser.add_argument("query", help="rsID, SPDI/variant_id, CA id, or small chr region.")
    variant_parser.add_argument("--limit", type=int, default=10)
    encode_parser = subparsers.add_parser("encode-search", help="Search ENCODE metadata.")
    encode_parser.add_argument("--type", default="Experiment", help="ENCODE object type.")
    encode_parser.add_argument(
        "--param",
        action="append",
        default=[],
        help="ENCODE search parameter as KEY=VALUE. May be repeated.",
    )
    subparsers.add_parser("collections", help="List Knowledge Graph collections.")
    aql_parser = subparsers.add_parser("aql", help="Run a read-oriented AQL query.")
    aql_parser.add_argument("query", help="AQL query string.")
    fetch_parser = subparsers.add_parser("fetch-url", help="Fetch a URL and save the response.")
    fetch_parser.add_argument("url")
    fetch_parser.add_argument("--label", default="manual_fetch")

    args = parser.parse_args(argv)
    setup_logging()

    if args.command == "auth-check":
        return cmd_auth_check(args)
    if args.command == "check":
        statuses = [
            check_catalog_docs(),
            check_catalog_index(),
            fetch_catalog_api("/", "catalog_api_root"),
            check_portal(),
            list_collections(),
        ]
        return 0 if all(status == 0 for status in statuses) else 1
    if args.command == "docs":
        return check_catalog_docs()
    if args.command == "docs-index":
        return check_catalog_index()
    if args.command == "portal":
        return check_portal()
    if args.command == "catalog-api":
        try:
            params = parse_repeated_params(args.param)
            return fetch_catalog_api(args.path, "catalog_api_fetch", params)
        except ValueError as exc:
            print(exc)
            return 2
    if args.command == "catalog-files":
        return fetch_catalog_files(args.limit, args.offset)
    if args.command == "gene":
        return explore_gene(args.query, args.limit)
    if args.command == "variant":
        return explore_variant(args.query, args.limit)
    if args.command == "encode-search":
        try:
            return encode_search(args.type, args.param)
        except ValueError as exc:
            print(exc)
            return 2
    if args.command == "collections":
        return list_collections()
    if args.command == "aql":
        return run_aql(args.query)
    if args.command == "fetch-url":
        return fetch_url(args.url, args.label)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
