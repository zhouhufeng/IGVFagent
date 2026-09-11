#!/usr/bin/env python3
"""IGVF Knowledge Graph traversal skill.

Iteratively walks the IGVF Catalog Knowledge Graph (db.catalog.igvf.org +
api.catalogkg.igvf.org) starting from a single entity (gene, variant, or
genomic region) and assembles a unified evidence pack across:

  • the entity itself (gene metadata / variant summary / region info)
  • directly linked entities  (variants, transcripts, proteins, regulatory
                                elements / cCREs, diseases, pathways)
  • second-degree relations   (per-variant phenotypes, biosamples, predictions;
                                per-cCRE method/biosample; per-protein partners)
  • cross-skill enrichment    (FAVOR annotations, enhancer-gene linkage,
                                IGVF single-cell datasets, prior literature)

Designed to be the orchestrator-friendly entry point: a single CLI call
yields a per-relation manifest set + a comprehensive markdown report + a
machine-readable JSON evidence pack ready to feed into downstream analyses
or another skill.

Subcommands

  gene <symbol>        Comprehensive gene-centric traversal.
  variant <id>         Variant-centric traversal (rsID / SPDI / HGVS).
  region <chr:start-end>
                       Region-centric traversal: genes + cCREs + variants in
                       window.
  aql <query>          Pass-through to direct ArangoDB AQL.
  write-playbook       Emit Docs/Skills/IGVF_KG_TRAVERSAL_SKILLS.md.
"""

from __future__ import annotations

import argparse
import csv
import io
import hashlib
import json
import logging
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Any, Iterable

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # _endpoints applies the IPv4-preferred DNS fix on import

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "KGTraversal"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
MANIFEST_DIR = DATA_DIR / "Manifests" / "KGTraversal"
CACHE_DIR = DATA_DIR / "Cache" / "KGTraversal"

CATALOG_API_BASE = _resolve_endpoint("catalog_api", "IGVF_CATALOG_API_BASE")
ARANGO_BASE = _resolve_endpoint("arango", "IGVF_ARANGO_BASE")
PORTAL_API_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_API_BASE")
FAVOR_API_BASE = _resolve_endpoint("favor", "FAVOR_API_BASE")

USER_AGENT = "IGVFdataAgent-KGTraversal/0.1"


# ----------------------------- Project plumbing ------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"kg_traversal_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    logging.info("Log file: %s", log_path)
    return log_path


def mkdirs() -> None:
    for d in (REPORT_DIR, MANIFEST_DIR, CACHE_DIR, SKILL_DOC_DIR):
        d.mkdir(parents=True, exist_ok=True)


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


# ----------------------------- HTTP helpers ----------------------------------

def request_headers(json_only: bool = True) -> dict[str, str]:
    h = {"User-Agent": USER_AGENT,
         "Accept": "application/json,*/*"}
    if not json_only:
        h["Accept"] = "*/*"
    if os.environ.get("IGVF_PORTAL_COOKIE"):
        h["Cookie"] = os.environ["IGVF_PORTAL_COOKIE"]
    return h


def fetch_json(url: str, timeout: int | None = None,
                 retries: int = 2, retry_backoff: float = 2.0
                 ) -> tuple[int, Any]:
    """GET a URL with a hard total timeout, JSON-or-text parsing, and
    auto-retry with exponential backoff on transient failures.

    The previous implementation passed ``timeout=60`` straight to
    ``urlopen``, but Python's stdlib treats that as a *per-read*
    inactivity timeout — if the server trickles even a few bytes
    before stalling, the call hangs forever. We hit this against
    api.catalogkg.igvf.org during a server-side slowdown: a ``kg gene
    APOE --call-* …`` run sat on one socket for > 30 min without
    raising. The fix:

      * `timeout` here is a *connection + per-read* ceiling, defaults
        to 30 s (env-tunable via ``IGVF_KG_HTTP_TIMEOUT``). Lower than
        the previous 60 s on purpose: a healthy Catalog response is
        sub-second; anything > 30 s means upstream trouble and we
        should fail-fast rather than hang.
      * `retries` controls re-attempts on timeout / network errors
        (default 2 retries → 3 total attempts).
      * Exponential backoff between retries (2 s, 4 s, 8 s).
      * Returns ``(0, {"network_error": ...})`` after the last retry
        fails, so callers see a clean error code instead of hanging.

    Override via env:
      IGVF_KG_HTTP_TIMEOUT   default 30 (seconds)
      IGVF_KG_HTTP_RETRIES   default 2
    """
    if timeout is None:
        try:
            timeout = int(os.environ.get("IGVF_KG_HTTP_TIMEOUT", "30"))
        except ValueError:
            timeout = 30
    try:
        retries = int(os.environ.get("IGVF_KG_HTTP_RETRIES",
                                        str(retries)))
    except ValueError:
        pass

    logging.info("GET %s", url)
    last_err: Exception | None = None
    for attempt in range(retries + 1):
        req = urllib.request.Request(url, headers=request_headers(),
                                       method="GET")
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                content = resp.read()
                try:
                    return resp.status, json.loads(content)
                except json.JSONDecodeError:
                    return resp.status, {
                        "text_response": content.decode(errors="replace"),
                        "url": url}
        except urllib.error.HTTPError as e:
            # 5xx + 408 + 429 are retryable; 4xx (except 408/429) are not.
            if e.code in (408, 429) or (500 <= e.code < 600):
                last_err = e
                if attempt < retries:
                    wait = retry_backoff * (2 ** attempt)
                    logging.warning("  HTTP %d on attempt %d/%d; "
                                      "retrying in %.0fs",
                                      e.code, attempt + 1, retries + 1, wait)
                    time.sleep(wait)
                    continue
                # Fall through after final retry to return the body
            body = e.read().decode(errors="replace")
            try:
                data = json.loads(body)
            except json.JSONDecodeError:
                data = {"http_error_body": body, "url": url}
            return e.code, data
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last_err = e
            reason = (getattr(e, "reason", None) or str(e))
            if attempt < retries:
                wait = retry_backoff * (2 ** attempt)
                logging.warning("  network error %r on attempt %d/%d; "
                                  "retrying in %.0fs",
                                  reason, attempt + 1, retries + 1, wait)
                time.sleep(wait)
                continue
            logging.error("  network error %r after %d attempts; giving up",
                            reason, retries + 1)
            return 0, {"network_error": str(reason), "url": url,
                        "attempts": retries + 1}
    # Unreachable, but keep mypy/pyright happy
    return 0, {"network_error": str(last_err), "url": url}


def catalog_get(path: str, **params) -> tuple[int, Any]:
    url = path if path.startswith("http") else CATALOG_API_BASE + path
    if params:
        url = f"{url}{'&' if '?' in url else '?'}" + urllib.parse.urlencode(
            {k: v for k, v in params.items() if v is not None}
        )
    return fetch_json(url)


# Measured against /api/genomic-elements/genes on 2026-09-09:
#   limit=1000 returns 500, so 500 is a server-side cap;
#   skip= and offset= are SILENTLY IGNORED -- they return page 1 again, so a
#     skip-based pager loops on the first page forever;
#   page= works, is 0-based, and strides by `limit` exactly (page 0 -> rows
#     0-4, page 1 -> rows 5-9, ... verified against a single limit=20 request
#     with no gaps and no overlaps).
# Anything that pages this API must therefore use `page`.
CATALOG_PAGE_MAX = 500


def catalog_paged(path: str, *, page_limit: int = CATALOG_PAGE_MAX,
                   max_pages: int = 40, **params) -> tuple[list[dict], dict]:
    """Read every page of a Catalog endpoint, or say why it stopped.

    Returns (rows, meta) with meta carrying `pages`, `returned_count`,
    `truncated` and `stopped_because`. The caller needs `truncated` to know
    whether "no matching record" means "none exist" or "none in what we
    looked at" -- reporting the second as the first is how this workflow
    claimed WT1 and HNF4A had no kidney regulatory records when WT1 has 390
    target-gene rows, 28 of them kidney.
    """
    rows: list[dict] = []
    seen: set = set()
    page = 0
    stopped = "exhausted"
    while page < max_pages:
        status, data = catalog_get(path, limit=page_limit, page=page, **params)
        if status != 200:
            stopped = f"http {status} on page {page}"
            break
        batch = listify(data)
        if not batch:
            break
        fresh = 0
        for r in batch:
            k = json.dumps(r, sort_keys=True, default=str)
            if k not in seen:
                seen.add(k)
                rows.append(r)
                fresh += 1
        if fresh == 0:
            # The endpoint is repeating itself: stop rather than spin.
            stopped = f"page {page} returned no new rows"
            break
        if len(batch) < page_limit:
            break
        page += 1
    else:
        stopped = f"hit max_pages={max_pages}"
    return rows, {"pages": page + 1, "returned_count": len(rows),
                   "truncated": stopped != "exhausted",
                   "stopped_because": stopped}


def portal_get(path: str, **params) -> tuple[int, Any]:
    url = PORTAL_API_BASE + path
    if params:
        url = f"{url}{'&' if '?' in url else '?'}" + urllib.parse.urlencode(
            {k: v for k, v in params.items() if v is not None}
        )
    return fetch_json(url)


def favor_get(path: str, **params) -> tuple[int, Any]:
    url = FAVOR_API_BASE + path
    if params:
        url = f"{url}{'&' if '?' in url else '?'}" + urllib.parse.urlencode(params)
    return fetch_json(url)


# --------------------------- KG record normalization -------------------------

# An error body is a dict too. {"message": "Not found", "code": "NOT_FOUND"}
# used to come back from listify() as ONE DATA ROW, so a dead endpoint looked
# like a database with one record in it.
_ERROR_KEYS = ({"message", "code"}, {"error"}, {"detail"})


def is_error_body(data: Any) -> bool:
    if not isinstance(data, dict):
        return False
    keys = set(data)
    return any(sig <= keys for sig in _ERROR_KEYS)


def listify(data: Any) -> list[dict]:
    if isinstance(data, list):
        return [d for d in data if isinstance(d, dict)]
    if isinstance(data, dict):
        for k in ("results", "items", "data", "@graph", "edges"):
            if k in data and isinstance(data[k], list):
                return [d for d in data[k] if isinstance(d, dict)]
        if is_error_body(data):
            return []
        return [data]
    return []


# ----------------------------- Manifest writing ------------------------------

def write_csv(path: Path, rows: list[dict], cols: list[str] | None = None) -> None:
    if not rows:
        path.write_text("")
        return
    cols = cols or sorted({k for r in rows for k in r})
    flat: list[dict] = []
    for r in rows:
        flat.append({k: _flatten_cell(r.get(k)) for k in cols})
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(flat)


def _flatten_cell(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, (list, tuple)):
        if all(isinstance(x, (str, int, float)) for x in v):
            return "; ".join(str(x) for x in v)
        return json.dumps(v, default=str)
    if isinstance(v, dict):
        return json.dumps(v, default=str)
    return str(v)


# --------------------------- Gene-centric traversal --------------------------

GENE_RELATIONS = {
    # label -> (api path, params for default lookup, extra kwargs to keep)
    "metadata":            ("/api/genes",                              {"name": "{symbol}"}),
    "variants":            ("/api/genes/variants",                      {"gene_name": "{symbol}"}),
    "transcripts":         ("/api/genes/transcripts",                   {"gene_name": "{symbol}"}),
    "proteins":            ("/api/genes/proteins",                      {"gene_name": "{symbol}"}),
    "regulatory_elements": ("/api/genes/genomic-elements",              {"gene_name": "{symbol}"}),
    "diseases":            ("/api/genes/diseases",                      {"gene_name": "{symbol}"}),
    # `verbose=true` inlines the full pathway node on each edge, so the
    # manifest carries the human-readable name directly. Without it the
    # edge only yields `pathways/R-HSA-n` refs and every name needs a
    # separate node lookup.
    "pathways":            ("/api/genes/pathways",                      {"gene_name": "{symbol}", "verbose": "true"}),
    "coding_variant_scores": ("/api/genes/coding-variants/scores",      {"gene_name": "{symbol}"}),
}

# Variant-side fan-out (called per variant when --depth>=2)
VARIANT_RELATIONS = {
    "summary":           ("/api/variants/summary",                     {"variant_id": "{vid}"}),
    "qtl_genes":         ("/api/variants/genes/summary",               {"variant_id": "{vid}", "verbose": "false"}),
    "phenotypes":        ("/api/variants/phenotypes",                  {"variant_id": "{vid}", "verbose": "false"}),
    "biosamples":        ("/api/variants/biosamples",                  {"variant_id": "{vid}", "verbose": "false"}),
    "genomic_elements":  ("/api/variants/genomic-elements",            {"variant_id": "{vid}"}),
    "predictions":       ("/api/variants/predictions",                 {"variant_id": "{vid}"}),
}


def _format_params(template: dict, **values) -> dict:
    out = {}
    for k, v in template.items():
        if isinstance(v, str) and "{" in v:
            try:
                out[k] = v.format(**values)
            except KeyError:
                out[k] = v
        else:
            out[k] = v
    return out


def gene_metadata(symbol: str) -> dict:
    status, data = catalog_get("/api/genes", name=symbol, limit=1)
    if status != 200:
        return {"_status": status, "_payload": data, "name": symbol}
    rows = listify(data)
    return rows[0] if rows else {"name": symbol, "_not_found": True}


# Nested-node key -> ArangoDB collection, so a lifted node can be collapsed
# back to the exact `<collection>/<id>` ref the non-verbose response gives.
_NODE_COLLECTION = {
    "gene": "genes", "pathway": "pathways", "parent_pathway": "pathways",
    "child_pathway": "pathways", "protein": "proteins",
    "transcript": "transcripts", "variant": "variants",
    "ontology_term": "ontology_terms", "drug": "drugs",
    "complex": "complexes", "genomic_element": "genomic_elements",
    "study": "studies",
}


def lift_nested_nodes(rows: list[dict]) -> list[dict]:
    """Flatten `verbose=true` edge rows for CSV.

    An edge fetched with ``verbose=true`` inlines the endpoint node as a
    nested dict (``{"pathway": {"_id": "R-HSA-174824", "name": "..."}}``),
    which ``_flatten_cell`` would JSON-dump into a single unreadable
    cell. Lift each nested node's id and name into ``<key>_id`` /
    ``<key>_name`` columns and collapse the dict back to the plain
    ``pathways/R-HSA-n`` ref the non-verbose response would have given —
    so the row shape stays comparable across runs, plus two useful
    columns. Rows without nested nodes pass through untouched.
    """
    for row in rows:
        if not isinstance(row, dict):
            continue
        for key, val in list(row.items()):
            if not isinstance(val, dict):
                continue
            node_id, node_name = val.get("_id"), val.get("name")
            if node_id is None and node_name is None:
                continue          # not a node record — leave it alone
            if node_id is not None:
                row.setdefault(f"{key}_id", node_id)
            if node_name is not None:
                row.setdefault(f"{key}_name", node_name)
            # Collapse back to the ref string (`pathways/R-HSA-174824`).
            # Never guess the collection by pluralising the key — that turns
            # `parent_pathway` into `parent_pathways/`. Unknown keys keep the
            # bare id.
            if node_id is not None:
                coll = _NODE_COLLECTION.get(key)
                row[key] = f"{coll}/{node_id}" if coll else node_id
    return rows


def fetch_gene_relations(symbol: str, limit: int) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for rel, (path, tpl) in GENE_RELATIONS.items():
        params = _format_params(tpl, symbol=symbol)
        params["limit"] = limit
        status, data = catalog_get(path, **params)
        if status == 200:
            rows = lift_nested_nodes(listify(data))
            out[rel] = rows
            logging.info("[%s] %d rows", rel, len(rows))
        else:
            out[rel] = []
            logging.warning("[%s] HTTP %s -> %r", rel, status,
                             str(data)[:120])
        time.sleep(0.1)
    return out


_SPDI_RE = re.compile(r"^NC_\d+\.\d+:\d+:[ACGTN]*:[ACGTN]*$")
_RSID_RE = re.compile(r"^rs\d+$")
_CA_RE   = re.compile(r"^CA\d+$")
_HGVS_RE = re.compile(r"^(?:NM|NR|NC|NP|ENS[TPG])[A-Z]?\d+(?:\.\d+)?:[gcpnmrn]\..+$",
                        re.IGNORECASE)


def resolve_to_spdi(vid: str) -> tuple[str, str]:
    """Resolve any variant identifier to its canonical SPDI form.

    Returns ``(spdi, input_form)``. The IGVF Catalog edge endpoints
    (``/api/variants/summary``, ``/predictions``, etc.) only reliably
    accept the canonical SPDI as the ``variant_id`` parameter — passing
    an rsID may return a *different* variant from ``/summary`` (the API
    treats it as a substring/fuzzy match), and passing an rsID to
    ``/predictions`` 400s. So we always pre-resolve before the fan-out.

    Heuristic:
      1. If ``vid`` already matches the SPDI pattern, return it as-is.
      2. Else call ``/api/variants?<rsid|ca_id|hgvs>=vid`` and pull
         ``_id`` from the first record.
      3. If the lookup yields no record, fall back to the raw input
         (caller will see empty results, which is informative).
    """
    s = vid.strip()
    if _SPDI_RE.match(s):
        return (s, "spdi")
    if _RSID_RE.match(s):
        param = "rsid"
    elif _CA_RE.match(s):
        param = "ca_id"
    elif _HGVS_RE.match(s):
        param = "hgvs"
    else:
        # Unknown form — best-effort, let the API decide.
        return (s, "unknown")
    status, data = catalog_get("/api/variants", **{param: s, "limit": 1})
    rows = listify(data) if status == 200 else []
    if not rows:
        logging.warning("resolve_to_spdi(%s): no record via %s",
                          s, param)
        return (s, param)
    spdi = rows[0].get("_id") or s
    logging.info("resolved %s (%s) -> %s", s, param, spdi)
    return (spdi, param)


def fetch_variant_relations(vid: str, limit: int = 25,
                              skip: tuple[str, ...] = ()) -> dict[str, list[dict]]:
    # Resolve to canonical SPDI before any edge call. This single change
    # turns every VARIANT_RELATIONS endpoint from "Variant not found"
    # into real data for any rsID / CA-ID / HGVS input.
    spdi, _form = resolve_to_spdi(vid)
    out: dict[str, list[dict]] = {}
    for rel, (path, tpl) in VARIANT_RELATIONS.items():
        if rel in skip:
            continue
        params = _format_params(tpl, vid=spdi)
        params["limit"] = limit
        status, data = catalog_get(path, **params)
        out[rel] = listify(data) if status == 200 else []
        time.sleep(0.05)
    return out


def variant_id_from_record(v: dict) -> str:
    """Extract a stable variant identifier from a KG variant *edge* record.

    Gene→variant edges from /api/genes/variants come back as edge documents
    like ``{"sequence_variant": "variants/NC_000019.10:44941484:A:G",
    "gene": "genes/ENSG00000130203", "method": "eQTL", ...}``. The actual
    variant identifier lives in ``sequence_variant`` (or in ``_from``/``_to``
    when AQL is used). We strip the collection prefix.
    """
    for k in ("sequence_variant", "variant", "variant_id", "spdi", "hgvs",
              "rsid", "_id", "_from"):
        if k in v and isinstance(v[k], str) and v[k]:
            val = v[k]
            if val.startswith("variants/"):
                return val.split("/", 1)[1]
            return val
    return ""


def biological_contexts_from_edges(edges: Iterable[dict]) -> list[str]:
    """Aggregate distinct biosample / cell-type contexts from gene-edge rows."""
    out: list[str] = []
    for e in edges or []:
        for k in ("biological_context", "biosample", "biosample_term_name",
                  "cell_type", "tissue"):
            v = e.get(k)
            if isinstance(v, str) and v and v not in out:
                out.append(v)
    return out


# --------------------------- Region helpers ----------------------------------

REGION_RE = re.compile(r"^(chr\w+):(\d+)-(\d+)$", re.I)


def parse_region(region: str) -> tuple[str, int, int]:
    m = REGION_RE.match(region.strip())
    if not m:
        raise SystemExit(f"Bad region: {region}. Use chr19:44900000-44910000")
    return m.group(1), int(m.group(2)), int(m.group(3))


def gene_region_string(gene: dict) -> str | None:
    chrom = gene.get("chr") or gene.get("chromosome") or gene.get("seqid")
    start = gene.get("start") or gene.get("start_position")
    end = gene.get("end") or gene.get("end_position")
    if chrom and start and end:
        return f"{chrom}:{start}-{end}"
    return None


# --------------------------- FAVOR side-call ---------------------------------

def favor_query_region(region: str, max_variants: int = 50) -> list[dict]:
    """Tries a small set of known FAVOR endpoints. Returns a list of variant
    annotation rows; empty if FAVOR is unreachable from this network."""
    chrom, start, end = parse_region(region)
    chrom_short = chrom.replace("chr", "")
    candidates = [
        ("/region", {"chr": chrom_short, "start": start, "end": end,
                      "limit": max_variants}),
        ("/api/range", {"chr": chrom, "start": start, "end": end,
                          "limit": max_variants}),
    ]
    for path, params in candidates:
        status, data = favor_get(path, **params)
        if status == 200:
            return listify(data)
    return []


# --------------------------- Single-cell side-call ---------------------------

def search_singlecell_for_gene(symbol: str, contexts: list[str] | None = None,
                                  limit: int = 25) -> list[dict]:
    """Find IGVF Portal AnalysisSets relevant to single-cell expression of a gene.

    Strategy:

    1. Use the biosample / cell-type contexts surfaced by the gene's eQTL and
       regulatory-element edges to drive a tissue-aware search of single-cell
       AnalysisSets — these are the cell types where the gene is *known* to
       be expressed.
    2. Always perform a broad gene-symbol fallback search.
    3. Filter to single-cell / single-nucleus assays (RNA / ATAC / multiome /
       SPLiT-seq / Perturb-seq).

    Returns a deduplicated list of candidate AnalysisSets; downstream skills
    (`single_cell_data_skills.py`, `splitseq_pipeline.py`,
    `multiome_10x_pipeline.py`) can drill into them.
    """
    sc_assays = (
        "Parse SPLiT-seq", "10x multiome",
        "single-nucleus RNA sequencing assay",
        "single-nucleus ATAC-seq",
        "single-cell RNA sequencing assay",
        "single-cell ATAC-seq",
        "Perturb-seq",
    )
    seen: dict[str, dict] = {}
    used_terms: set[str] = set()

    def _normalize(ctx: str) -> str:
        t = re.sub(r"[(),]", " ", ctx).strip()
        return re.sub(r"\s+", " ", t)

    for ctx in (contexts or [])[:6]:
        term = _normalize(ctx)
        if not term or term.lower() in used_terms:
            continue
        used_terms.add(term.lower())
        for assay in sc_assays:
            status, data = portal_get(
                "/search/", type="AnalysisSet", format="json",
                limit=limit, preferred_assay_titles=assay, searchTerm=term,
            )
            if status != 200 or not isinstance(data, dict):
                continue
            for g in (data.get("@graph") or [])[:limit]:
                acc = g.get("accession") or ""
                if not acc or acc in seen:
                    continue
                seen[acc] = {
                    "accession":              acc,
                    "preferred_assay_titles": assay,
                    "description":            (g.get("description") or "")[:280],
                    "summary":                g.get("summary", ""),
                    "lab":                    (g.get("lab") or {}).get("title", "")
                                                if isinstance(g.get("lab"), dict)
                                                else "",
                    "status":                 g.get("status", ""),
                    "matched_via":            f"{ctx}|{assay}",
                }

    # Symbol fallback: pure searchTerm=symbol with each assay
    for assay in sc_assays:
        status, data = portal_get(
            "/search/", type="AnalysisSet", format="json",
            limit=limit, preferred_assay_titles=assay, searchTerm=symbol,
        )
        if status != 200 or not isinstance(data, dict):
            continue
        for g in (data.get("@graph") or [])[:limit]:
            acc = g.get("accession") or ""
            if not acc or acc in seen:
                continue
            seen[acc] = {
                "accession":              acc,
                "preferred_assay_titles": assay,
                "description":            (g.get("description") or "")[:280],
                "summary":                g.get("summary", ""),
                "lab":                    (g.get("lab") or {}).get("title", "")
                                            if isinstance(g.get("lab"), dict)
                                            else "",
                "status":                 g.get("status", ""),
                "matched_via":            f"symbol={symbol}|{assay}",
            }
    return list(seen.values())


# --------------------------- Linkage side-call -------------------------------

def target_gene_id(row: dict) -> str:
    """The Ensembl gene this row is evidence FOR, bare of any prefix.

    The field arrives as "genes/ENSG00000184937", and may be an embedded
    object instead of a string. Splitting it on "." rather than "/" silently
    matched nothing, which reads exactly like "this gene has no records".
    """
    g = row.get("gene")
    if isinstance(g, dict):
        g = g.get("_id") or g.get("gene_id") or g.get("id") or g.get("name") or ""
    return str(g or "").split("/")[-1].split(".")[0].strip()


def fetch_linkage_for_region(region: str, limit: int = 25,
                              gene_id: "Optional[str]" = None,
                              exhaustive: bool = True,
                              max_pages: int = 40) -> dict:
    """Element-to-gene edges over a region, optionally for ONE target gene.

    THE TRAP THIS EXISTS TO CLOSE. The endpoint returns edges for regulatory
    elements that OVERLAP the submitted region. A row is evidence for the
    gene named in its own `gene` field -- not for the gene whose locus was
    used to build the region. Those are mostly different genes: of the first
    25 rows over the WT1 locus, 23 distinct target genes appear and NONE of
    them is WT1. Reporting them under WT1 asserts element-to-WT1 links that
    the source does not claim.

    Pass `gene_id` to keep only rows whose target gene matches, which is what
    a gene-centric question means. Without it the rows are returned unfiltered
    and `target_gene_id` is None in the meta, so a caller cannot mistake them
    for gene-specific evidence.

    /api/regulatory-regions/genes is NOT tried any more: it answers 404 for
    every region, so every call paid for a failed request and then fell
    through, and any future non-404 error would have been masked the same way.
    """
    out: dict = {"region_predictions": [], "qtl_links": [],
                  "meta": {"endpoint": "/api/genomic-elements/genes",
                            "region": region, "target_gene_id": gene_id}}
    if exhaustive:
        rows, meta = catalog_paged("/api/genomic-elements/genes",
                                    region=region, max_pages=max_pages)
    else:
        status, data = catalog_get("/api/genomic-elements/genes",
                                    region=region, limit=limit)
        rows = listify(data) if status == 200 else []
        meta = {"pages": 1, "returned_count": len(rows),
                 "truncated": len(rows) >= limit,
                 "stopped_because": f"single request, limit={limit}"}
    out["meta"].update(meta)
    out["meta"]["rows_before_gene_filter"] = len(rows)
    out["meta"]["distinct_target_genes"] = len({target_gene_id(r) for r in rows})

    if gene_id:
        want = str(gene_id).split("/")[-1].split(".")[0].strip()
        kept = [r for r in rows if target_gene_id(r) == want]
        out["meta"]["rows_for_target_gene"] = len(kept)
        out["meta"]["rows_dropped_other_genes"] = len(rows) - len(kept)
        rows = kept

    # Predictions and observed measurements are different kinds of claim and
    # must not be pooled: a prediction carries a model score with p_value_adj
    # and significant both null, so calling it significant or non-significant
    # invents a status the source never gave.
    obs = [r for r in rows if str(r.get("class", "")).lower().startswith("observed")]
    pred = [r for r in rows if str(r.get("class", "")).lower() == "prediction"]
    out["meta"]["observed_count"] = len(obs)
    out["meta"]["prediction_count"] = len(pred)
    out["meta"]["other_class_count"] = len(rows) - len(obs) - len(pred)
    out["region_predictions"] = rows
    out["observed"] = obs
    out["predictions"] = pred
    return out


def evidence_note(meta: dict) -> str:
    """One line a report can print instead of an unqualified absence claim."""
    n = meta.get("rows_for_target_gene", meta.get("returned_count", 0))
    parts = [f"{n} row(s) for the target gene"]
    if meta.get("rows_dropped_other_genes"):
        parts.append(f"{meta['rows_dropped_other_genes']} dropped as other "
                      f"target genes")
    parts.append(f"{meta.get('observed_count', 0)} observed / "
                  f"{meta.get('prediction_count', 0)} predicted")
    if meta.get("truncated"):
        parts.append(f"RETRIEVAL WAS TRUNCATED ({meta.get('stopped_because')}) "
                      f"-- absence cannot be concluded from this")
    else:
        parts.append(f"retrieval exhausted over {meta.get('pages', 1)} page(s)")
    return "; ".join(parts)



# --------------------------- Answer verification -----------------------------

# The fields an answer is allowed to quote about a regulatory link. Each one
# must come from the SAME source record: a row is a measurement, and taking
# the biosample from one row, the score from another and the p-value from a
# third produces a sentence that is false about every record it was built
# from while every individual value is real.
CITABLE_FIELDS = ("gene", "genomic_element", "biological_context", "method",
                   "class", "source", "source_url", "score", "effect_size",
                   "log2FC", "p_value", "p_value_adj", "significant",
                   "crispr_modality", "files_filesets")


def record_id(row: dict) -> str:
    """Stable id for one source record, from its identifying fields."""
    basis = json.dumps(
        {k: row.get(k) for k in ("gene", "genomic_element", "source_url",
                                  "biological_context", "method", "class",
                                  "score", "p_value_adj")},
        sort_keys=True, default=str)
    return hashlib.sha256(basis.encode()).hexdigest()[:12]


def evidence_records(rows: "list[dict]") -> "list[dict]":
    """Citable, record-atomic view of retrieved rows, each with a record_id."""
    out = []
    for r in rows:
        rec = {"record_id": record_id(r)}
        for k in CITABLE_FIELDS:
            rec[k] = r.get(k)
        out.append(rec)
    return out


def _same(a: Any, b: Any) -> bool:
    """Field comparison that tolerates presentation, not substance.

    Numbers are compared with a relative tolerance because an answer rounds
    (0.7911 for 0.79107...); strings are compared case-insensitively after
    stripping any "genes/" style prefix. None matches only None -- a null
    adjusted p-value quoted as a number is exactly the confusion between
    prediction scores and significance that this is here to catch.
    """
    if a is None or b is None:
        return a is None and b is None
    if isinstance(a, bool) or isinstance(b, bool):
        return bool(a) == bool(b)
    try:
        fa, fb = float(a), float(b)
        if fa == fb:
            return True
        scale = max(abs(fa), abs(fb), 1e-12)
        return abs(fa - fb) / scale < 1e-3
    except (TypeError, ValueError):
        pass
    sa = str(a).strip().lower().split("/")[-1]
    sb = str(b).strip().lower().split("/")[-1]
    return sa == sb


def verify_claim(claim: "dict", records: "list[dict]") -> dict:
    """Does one claimed row correspond to a single retrieved record?

    Returns {ok, record_id, reason, fields_found_elsewhere}. When no single
    record carries every claimed field but each value does occur somewhere,
    the verdict is CROSS_RECORD -- the assembly failure the evaluation asked
    us to reject, and the one that is invisible to spot-checking because
    every value in the sentence is genuine.
    """
    fields = {k: v for k, v in (claim or {}).items()
              if k in CITABLE_FIELDS and v is not None}
    if not fields:
        return {"ok": False, "record_id": None, "reason": "NO_CITABLE_FIELDS",
                 "fields_found_elsewhere": []}
    for rec in records:
        if all(_same(rec.get(k), v) for k, v in fields.items()):
            return {"ok": True, "record_id": rec["record_id"], "reason": "",
                     "fields_found_elsewhere": []}
    elsewhere = [k for k, v in fields.items()
                 if any(_same(rec.get(k), v) for rec in records)]
    if len(elsewhere) == len(fields):
        return {"ok": False, "record_id": None, "reason": "CROSS_RECORD",
                 "fields_found_elsewhere": elsewhere}
    missing = [k for k in fields if k not in elsewhere]
    return {"ok": False, "record_id": None, "reason": "NOT_IN_SOURCE",
             "fields_found_elsewhere": elsewhere, "unmatched_fields": missing}


def verify_table(claims: "list[dict]", records: "list[dict]") -> dict:
    """Verify every claimed row. Returns a verdict plus per-row detail."""
    rows = [dict(verify_claim(c, records), claim_index=i)
            for i, c in enumerate(claims or [])]
    bad = [r for r in rows if not r["ok"]]
    return {"n_claims": len(rows), "n_verified": len(rows) - len(bad),
             "n_rejected": len(bad),
             "cross_record": sum(1 for r in bad if r["reason"] == "CROSS_RECORD"),
             "not_in_source": sum(1 for r in bad if r["reason"] == "NOT_IN_SOURCE"),
             "ok": not bad, "rows": rows}

# --------------------------- Literature side-call ----------------------------

def call_literature_validate(symbol: str, context: list[str], top: int = 10) -> list[dict]:
    """Lazily import the reference skill's PubMed search to keep this skill
    standalone if the reference module is unavailable."""
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import importlib
        ref = importlib.import_module("reference_skill")
    except Exception as e:
        logging.warning("reference_skill unavailable: %s", e)
        return []
    q = symbol
    if context:
        q = f"{symbol} AND ({' AND '.join(context)})"
    try:
        recs = ref.pubmed_search(q, limit=top)
        recs += ref.semanticscholar_search(q, limit=top)
        deduped = ref.dedup_records(recs)
        deduped.sort(key=lambda r: ref.score_relevance(r, [symbol] + (context or [])),
                     reverse=True)
        return deduped[:top]
    except Exception as e:
        logging.warning("Literature retrieval failed: %s", e)
        return []


# --------------------------- Reporting ---------------------------------------

def summarize_relation(rows: list[dict]) -> str:
    if not rows:
        return "_no records returned_"
    lines = []
    sample = rows[:5]
    for r in sample:
        lines.append("- " + _row_oneline(r))
    if len(rows) > 5:
        lines.append(f"- _… and {len(rows) - 5} more_")
    return "\n".join(lines)


def _row_oneline(r: dict) -> str:
    """Best-effort one-line printable summary for a KG row."""
    if not isinstance(r, dict):
        return str(r)
    # Try the most useful fields first
    parts = []
    # `<key>_name` (lifted from a verbose edge) comes first: on an edge row
    # the bare `name` is the edge label ("belongs to"), not the node's name.
    for k in ("pathway_name", "name", "symbol", "term_name", "label", "id",
              "_key", "_id", "spdi", "rsid", "hgvs", "variant_id",
              "uniprot_id", "accession"):
        if r.get(k):
            parts.append(str(r[k])[:80])
            break
    for k in ("biotype", "consequence", "type", "method", "tissue",
              "biosample", "phenotype", "score", "p_value", "log2fc",
              "category"):
        if r.get(k):
            parts.append(f"{k}={str(r[k])[:60]}")
    return " | ".join(parts) if parts else json.dumps(r, default=str)[:160]


def render_gene_report(symbol: str, meta: dict, rels: dict[str, list[dict]],
                        deep_var: dict[str, dict[str, list[dict]]],
                        favor_rows: list[dict],
                        singlecell_hits: list[dict],
                        linkage: dict[str, list[dict]],
                        literature: list[dict],
                        manifest_paths: dict[str, Path],
                        out_path: Path) -> Path:
    region = gene_region_string(meta)
    coords = (f"`{region}` ({meta.get('strand','?')} strand)"
              if region else "_no region in metadata_")
    biotype = meta.get("biotype") or meta.get("gene_type") or ""
    ensembl_id = meta.get("ensembl_id") or meta.get("gene_id") or meta.get("_key", "")
    n_variants = len(rels.get("variants", []))
    n_re = len(rels.get("regulatory_elements", []))
    n_dis = len(rels.get("diseases", []))
    n_path = len(rels.get("pathways", []))
    n_prot = len(rels.get("proteins", []))
    n_cv = len(rels.get("coding_variant_scores", []))

    lines = [
        f"# IGVF KG Traversal: gene `{symbol}`",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        "## Gene metadata",
        "",
        f"- Ensembl ID: `{ensembl_id}`",
        f"- Coordinates: {coords}",
        (f"- Biotype: `{biotype}`" if biotype else ""),
        (f"- Description: {meta.get('description','')[:280]}"
         if meta.get("description") else ""),
        "",
        "## Direct neighbors (Catalog API)",
        "",
        f"| Relation | n | manifest |",
        "|---|---:|---|",
    ]
    for label, key in [
        ("Variants",                  "variants"),
        ("Coding-variant scores",     "coding_variant_scores"),
        ("Regulatory elements (cCREs)", "regulatory_elements"),
        ("Transcripts",               "transcripts"),
        ("Proteins",                  "proteins"),
        ("Diseases",                  "diseases"),
        ("Pathways",                  "pathways"),
    ]:
        n = len(rels.get(key, []))
        path = manifest_paths.get(key)
        rel_path = path.relative_to(ROOT) if path else None
        lines.append(f"| {label} | {n} | "
                      f"{'`'+str(rel_path)+'`' if rel_path else ''} |")
    lines += ["", "### Variants (sample)", "",
               summarize_relation(rels.get("variants", [])), "",
               "### Coding variant scores (sample)", "",
               summarize_relation(rels.get("coding_variant_scores", [])), "",
               "### Regulatory elements / cCRE links (sample)", "",
               summarize_relation(rels.get("regulatory_elements", [])), "",
               "### Diseases (sample)", "",
               summarize_relation(rels.get("diseases", [])), "",
               "### Pathways (sample)", "",
               summarize_relation(rels.get("pathways", [])), "",
               "### Proteins (sample)", "",
               summarize_relation(rels.get("proteins", [])), ""]

    if deep_var:
        lines += ["## Per-variant fan-out (depth ≥ 2)", ""]
        for vid, sub in list(deep_var.items())[:10]:
            lines += [f"### Variant `{vid}`", ""]
            for k in ("summary", "qtl_genes", "phenotypes", "biosamples",
                       "genomic_elements", "predictions"):
                if not sub.get(k):
                    continue
                lines.append(f"- **{k}** ({len(sub[k])}): "
                             + _row_oneline(sub[k][0]))
            lines.append("")
        if len(deep_var) > 10:
            lines.append(f"_… and {len(deep_var) - 10} more variants in the manifest._")
            lines.append("")

    if linkage:
        lmeta = linkage.get("meta") or {}
        tgt = lmeta.get("target_gene_id")
        lines += ["## Enhancer-gene linkage"
                   + (f" — edges whose TARGET GENE is `{tgt}`" if tgt
                      else " — edges over this region, TARGET GENES VARY"),
                   "",
                   f"_Retrieval: {evidence_note(lmeta)}._", ""]
        if not tgt:
            others = {}
            for r in linkage.get("region_predictions", []):
                g = target_gene_id(r)
                others[g] = others.get(g, 0) + 1
            top = ", ".join(f"`{g}` ({n})" for g, n in
                             sorted(others.items(), key=lambda kv: -kv[1])[:8])
            lines += [f"Target genes present: {top or 'none'}", ""]
        if lmeta.get("prediction_count") and lmeta.get("observed_count") == 0:
            lines += ["> Every row below is `class=prediction`: a model score, "
                      "with `p_value_adj` and `significant` both null. These "
                      "are not statistically significant or non-significant "
                      "findings and must not be described as either.", ""]
        lines += [summarize_relation(linkage.get("region_predictions", [])), ""]
    if favor_rows:
        lines += ["## FAVOR functional annotation (region)",
                   "",
                   summarize_relation(favor_rows), ""]
    if singlecell_hits:
        lines += ["## IGVF single-cell datasets mentioning this gene", "",
                   "| Accession | Assay | Lab | Description |",
                   "|---|---|---|---|"]
        for h in singlecell_hits[:20]:
            lines.append(
                f"| {h.get('accession','')} | "
                f"{h.get('preferred_assay_titles','')} | {h.get('lab','')} | "
                f"{(h.get('description') or '')[:120]} |"
            )
        lines.append("")
    if literature:
        lines += ["## Literature corroboration", ""]
        for r in literature[:10]:
            lines.append(
                f"- {r.get('title','(untitled)')} — "
                f"**{r.get('journal','')}** ({r.get('year','')})"
                + (f" · https://doi.org/{r['doi']}" if r.get('doi') else "")
            )
        lines.append("")

    lines += ["## Manifests", ""]
    for k, p in manifest_paths.items():
        if p:
            lines.append(f"- **{k}**: `{p.relative_to(ROOT)}`")
    out_path.write_text("\n".join(l for l in lines if l is not None))
    return out_path


# --------------------------- Subcommands -------------------------------------

def cmd_genes(args: argparse.Namespace) -> Path:
    """Several genes' regulatory evidence in ONE tool call.

    WHY THIS EXISTS. The six-gene prompt took over 25 minutes and never
    finished, where a single gene takes 86-350 seconds. Measurement showed
    retrieval is not the cost: an exhaustive paged pull is 6-8 seconds per
    gene (WT1 8,926 rows over 18 pages in 8.1s; GATA3 11,345 over 23 in
    6.3s), so six genes is about 45 seconds of API time. The cost is the
    AGENT LOOP -- 5-12 seconds per iteration of model latency, and two
    single-gene runs hit the 25-iteration cap on their own. Six genes
    multiplies the orchestration, not the fetching.

    So this does the whole job in one call: for each gene, resolve the
    target gene id, pull element-to-gene edges exhaustively, keep only rows
    whose OWN target gene matches, split observations from predictions, and
    apply the tissue filter -- then write one combined manifest and one
    summary. The agent reads a single result instead of driving forty
    iterations by hand.

    Tissue matching takes a list and an exclusion list for the reason the
    independent audit found: `kidney` alone misses "renal cortical epithelial
    cell", and `renal` alone also matches "ADRENAL gland". The auditors hit
    the second themselves and had to redo their first pass.
    """
    setup_logging(); mkdirs()
    ts = timestamp()
    symbols = [s.strip() for s in args.symbols.split(",") if s.strip()]
    label = safe_label(args.label or f"genes_{len(symbols)}")
    out_dir = REPORT_DIR / f"{ts}_{label}"
    (out_dir / "Manifests").mkdir(parents=True, exist_ok=True)

    include = [t.strip().lower() for t in (args.tissue or "").split(",") if t.strip()]
    exclude = [t.strip().lower() for t in (args.exclude_tissue or "").split(",") if t.strip()]

    per_gene, combined = {}, []
    for sym in symbols:
        # The gene id and the region both come from the gene's own metadata.
        # target_gene_id() takes a ROW and says which gene that row is
        # evidence for -- it does not resolve a symbol, and using it that way
        # raised AttributeError on a str.
        gmeta = gene_metadata(sym) or {}
        gid = gmeta.get("_id") or gmeta.get("gene_id") or gmeta.get("_key")
        region = gene_region_string(gmeta)
        entry = {"symbol": sym, "target_gene_id": gid, "region": region}
        if not gid or not region:
            entry["error"] = ("could not resolve this gene's id or region, so "
                               "NO linkage was retrieved. That is not the same "
                               "as the gene having no evidence.")
            per_gene[sym] = entry
            continue
        link = fetch_linkage_for_region(region, gene_id=gid,
                                         exhaustive=not args.no_exhaustive,
                                         max_pages=args.max_pages)
        meta = link.get("meta") or {}
        obs = link.get("observed") or []
        pred = link.get("predictions") or []
        rows = link.get("region_predictions") or []

        def tissue_ok(r):
            hay = " ".join(str(v) for v in r.values()).lower()
            if include and not any(t in hay for t in include):
                return False
            if exclude and any(t in hay for t in exclude):
                return False
            return True

        t_obs = [r for r in obs if tissue_ok(r)]
        t_pred = [r for r in pred if tissue_ok(r)]
        entry.update({
            "rows_retrieved": meta.get("rows_before_gene_filter", len(rows)),
            "rows_on_target_gene": meta.get("rows_for_target_gene", len(rows)),
            "rows_dropped_other_genes": meta.get("rows_dropped_other_genes", 0),
            "observations": len(obs), "predictions": len(pred),
            "tissue_filtered_observations": len(t_obs),
            "tissue_filtered_predictions": len(t_pred),
            # Two DIFFERENT statuses, kept apart on purpose: whether the
            # Catalog traversal saw everything, and whether this report shows
            # everything it saw. Conflating them produced an answer that
            # called an exhaustive retrieval "TRUNCATED".
            "catalog_retrieval": ("exhausted" if not meta.get("truncated")
                                   else "truncated"),
            "catalog_stopped_because": meta.get("stopped_because"),
            "pages": meta.get("pages"),
        })
        per_gene[sym] = entry
        for r in t_obs + t_pred:
            combined.append({"query_gene": sym, **r})

    man = out_dir / "Manifests" / "linkage_all_genes.csv"
    if combined:
        cols = sorted({k for r in combined for k in r})
        with open(man, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in combined:
                w.writerow(r)

    report = {
        "generated": ts, "symbols": symbols,
        "tissue_include": include, "tissue_exclude": exclude,
        "genes": per_gene,
        "combined_manifest": str(man) if combined else None,
        "combined_rows": len(combined),
        "note": ("`catalog_retrieval` says whether the Catalog traversal was "
                  "exhaustive. It is NOT a statement about how much of the "
                  "manifest any later read covered -- rank the manifest with "
                  "`igvfagent artifact top` rather than reading an excerpt "
                  "if the question asks for the highest-scoring rows."),
        # A hosted answer suggested a follow-up using `--filter` and
        # `--export_full`, neither of which exists on this skill (the model
        # appears to have borrowed `--filters` from catalog_query_skill).
        # Printing the real commands removes the need to guess one.
        "next_commands": [
            (f"igvfagent artifact top --path {man} --column score --n 10 "
              f"--where kidney,renal --exclude adrenal"),
            (f"igvfagent kg genes {args.symbols} --tissue kidney,renal "
              f"--exclude-tissue adrenal"),
        ],
    }
    path = out_dir / "report.json"
    path.write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps({k: v for k, v in report.items() if k != "genes"}, indent=2,
                      default=str))
    for sym, e in per_gene.items():
        print(f"  {sym:8} on-target {e.get('rows_on_target_gene', 0):>6,}  "
              f"obs {e.get('observations', 0):>4}  pred {e.get('predictions', 0):>6,}  "
              f"tissue obs/pred {e.get('tissue_filtered_observations', 0)}/"
              f"{e.get('tissue_filtered_predictions', 0)}  "
              f"[{e.get('catalog_retrieval', e.get('error', '?'))}]")
    print(f"\nReport: {path}")
    return path


def cmd_gene(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    symbol = args.symbol
    label = safe_label(args.label or f"gene_{symbol}")
    out_dir = REPORT_DIR / f"{ts}_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifests_dir = out_dir / "Manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)

    meta = gene_metadata(symbol)
    rels = fetch_gene_relations(symbol, limit=args.limit)

    # Save per-relation manifests
    manifest_paths: dict[str, Path] = {}
    for rel, rows in rels.items():
        path = manifests_dir / f"{rel}.csv"
        write_csv(path, rows)
        manifest_paths[rel] = path

    # Per-variant fan-out at depth ≥ 2
    deep_var: dict[str, dict[str, list[dict]]] = {}
    if args.depth >= 2 and rels.get("variants"):
        deep_rows: list[dict] = []
        for v in rels["variants"][: args.max_variants]:
            vid = variant_id_from_record(v)
            if not vid:
                continue
            sub = fetch_variant_relations(vid, limit=args.subvariant_limit)
            deep_var[vid] = sub
            for k, vrows in sub.items():
                for vr in vrows:
                    deep_rows.append({"_variant": vid, "_relation": k, **vr})
        if deep_rows:
            depth_path = manifests_dir / "variant_depth2.csv"
            write_csv(depth_path, deep_rows)
            manifest_paths["variant_depth2"] = depth_path

    # Side-calls
    region = gene_region_string(meta)
    favor_rows: list[dict] = []
    linkage: dict[str, list[dict]] = {}
    if region and args.call_favor:
        favor_rows = favor_query_region(region, max_variants=args.favor_max)
        if favor_rows:
            p = manifests_dir / "favor.csv"
            write_csv(p, favor_rows); manifest_paths["favor"] = p
    if region and args.call_linkage:
        # gene_id is REQUIRED here. Without it these rows are edges for
        # whatever genes happen to sit near this locus, and 23 of the first
        # 25 over the WT1 region target other genes.
        linkage = fetch_linkage_for_region(
            region, limit=args.limit, gene_id=(meta or {}).get("_id"),
            exhaustive=not getattr(args, "no_exhaustive_linkage", False))
        if linkage.get("region_predictions"):
            p = manifests_dir / "linkage_region_predictions.csv"
            write_csv(p, linkage["region_predictions"])
            manifest_paths["linkage_region_predictions"] = p
        for k in ("observed", "predictions"):
            if linkage.get(k):
                p = manifests_dir / f"linkage_{k}.csv"
                write_csv(p, linkage[k])
                manifest_paths[f"linkage_{k}"] = p
        # Record-atomic view: one row per source record, each with a
        # record_id. An answer should quote a row of THIS table, so every
        # field it states came from one measurement and can be checked back
        # against it with verify_claim().
        if linkage.get("region_predictions"):
            recs = evidence_records(linkage["region_predictions"])
            p = manifests_dir / "linkage_evidence_records.csv"
            write_csv(p, recs)
            manifest_paths["linkage_evidence_records"] = p

    singlecell_hits: list[dict] = []
    if args.call_singlecell:
        contexts = (biological_contexts_from_edges(rels.get("variants", []))
                    + biological_contexts_from_edges(rels.get("regulatory_elements", [])))
        singlecell_hits = search_singlecell_for_gene(
            symbol, contexts=contexts, limit=args.limit)
        if singlecell_hits:
            p = manifests_dir / "single_cell_datasets.csv"
            write_csv(p, singlecell_hits); manifest_paths["single_cell_datasets"] = p

    literature: list[dict] = []
    if args.call_literature:
        literature = call_literature_validate(
            symbol, args.literature_context or [], top=args.literature_top,
        )
        if literature:
            p = manifests_dir / "literature.csv"
            write_csv(p, literature); manifest_paths["literature"] = p

    # Evidence pack: full JSON dump
    pack = {
        "symbol": symbol,
        "metadata": meta,
        "relations": rels,
        "variant_depth2": deep_var,
        "favor": favor_rows,
        "linkage": linkage,
        "single_cell_datasets": singlecell_hits,
        "literature": literature,
    }
    pack_path = out_dir / "evidence_pack.json"
    pack_path.write_text(json.dumps(pack, indent=2, default=str))

    report = render_gene_report(symbol, meta, rels, deep_var, favor_rows,
                                  singlecell_hits, linkage, literature,
                                  manifest_paths,
                                  out_dir / f"gene_{safe_label(symbol)}_report.md")
    print(f"Report:        {report}")
    print(f"Evidence pack: {pack_path}")
    print(f"Manifests:     {manifests_dir}")
    return report


def cmd_variant(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    vid = args.variant
    label = safe_label(args.label or f"variant_{vid}")
    out_dir = REPORT_DIR / f"{ts}_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifests_dir = out_dir / "Manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)
    rels = fetch_variant_relations(vid, limit=args.limit)
    manifest_paths = {}
    for k, rows in rels.items():
        p = manifests_dir / f"{k}.csv"
        write_csv(p, rows); manifest_paths[k] = p

    favor_rows: list[dict] = []
    if args.call_favor:
        # Best-effort: parse SPDI -> region
        m = re.match(r"NC_(\d+)\.\d+:(\d+):", vid)
        if m:
            chrom = f"chr{int(m.group(1))}"
            pos = int(m.group(2))
            favor_rows = favor_query_region(f"{chrom}:{pos-1}-{pos+1}",
                                             max_variants=10)
    literature: list[dict] = []
    if args.call_literature:
        literature = call_literature_validate(vid, args.literature_context or [],
                                                top=args.literature_top)

    pack = {"variant": vid, "relations": rels, "favor": favor_rows,
             "literature": literature}
    (out_dir / "evidence_pack.json").write_text(json.dumps(pack, indent=2,
                                                            default=str))
    lines = [f"# IGVF KG Traversal: variant `{vid}`",
             f"\nGenerated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n",
             "## Direct relations\n",
             "| Relation | n |\n|---|---:|"]
    for k, rows in rels.items():
        lines.append(f"| {k} | {len(rows)} |")
    for k, rows in rels.items():
        lines += [f"\n### {k} (sample)\n", summarize_relation(rows)]
    if favor_rows:
        lines += ["\n## FAVOR\n", summarize_relation(favor_rows)]
    if literature:
        lines += ["\n## Literature\n"]
        for r in literature[:10]:
            lines.append(f"- {r.get('title','')} — **{r.get('journal','')}** "
                          f"({r.get('year','')})")
    report = out_dir / f"variant_{safe_label(vid)}_report.md"
    report.write_text("\n".join(lines))
    print(f"Report:        {report}")
    return report


def cmd_region(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    region = args.region
    chrom, start, end = parse_region(region)
    label = safe_label(args.label or f"region_{chrom}_{start}_{end}")
    out_dir = REPORT_DIR / f"{ts}_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifests_dir = out_dir / "Manifests"; manifests_dir.mkdir(parents=True, exist_ok=True)

    # Genes in region
    s, d = catalog_get("/api/genes", region=region, limit=args.limit)
    genes = listify(d) if s == 200 else []
    write_csv(manifests_dir / "genes_in_region.csv", genes)

    # cCREs / regulatory elements in region
    s, d = catalog_get("/api/genomic-elements", region=region, limit=args.limit)
    ccres = listify(d) if s == 200 else []
    write_csv(manifests_dir / "regulatory_elements.csv", ccres)

    # Linkage for region. Deliberately NOT gene-filtered: this IS a region
    # question, so edges to any target gene in the window are the answer.
    # The report says which genes they target so they are not read as
    # evidence for one gene.
    linkage = fetch_linkage_for_region(region, limit=args.limit,
                                        exhaustive=not getattr(args, "no_exhaustive_linkage", False))
    write_csv(manifests_dir / "linkage_region_predictions.csv",
               linkage.get("region_predictions", []))

    favor_rows: list[dict] = []
    if args.call_favor:
        favor_rows = favor_query_region(region, max_variants=args.favor_max)
        write_csv(manifests_dir / "favor.csv", favor_rows)

    pack = {"region": region, "genes": genes, "regulatory_elements": ccres,
             "linkage": linkage, "favor": favor_rows}
    (out_dir / "evidence_pack.json").write_text(json.dumps(pack, indent=2,
                                                            default=str))
    lines = [f"# IGVF KG Traversal: region `{region}`",
             f"\nGenerated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n",
             f"Genes: **{len(genes)}**, regulatory elements: **{len(ccres)}**, "
             f"linkage rows: **{len(linkage.get('region_predictions', []))}**, "
             f"FAVOR variants: **{len(favor_rows)}**\n"]
    lines += ["\n## Genes\n", summarize_relation(genes)]
    lines += ["\n## Regulatory elements (cCREs)\n", summarize_relation(ccres)]
    lines += ["\n## Linkage region predictions\n",
               summarize_relation(linkage.get("region_predictions", []))]
    if favor_rows:
        lines += ["\n## FAVOR\n", summarize_relation(favor_rows)]
    report = out_dir / f"region_{safe_label(region)}_report.md"
    report.write_text("\n".join(lines))
    print(f"Report: {report}")
    return report


def cmd_aql(args: argparse.Namespace) -> Path:
    """Pass-through to the existing igvf_client AQL helper for direct
    ArangoDB access (same as `python3 Scripts/igvf_client.py aql ...`)."""
    setup_logging(); mkdirs()
    import importlib
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    client = importlib.import_module("igvf_client")
    # The igvf_client AQL helper writes to Data/; just print its output here.
    rc = client.run_aql(args.query, limit=args.limit, log_only=False)
    return Path(rc) if isinstance(rc, str) else Path("")


def cmd_write_playbook(_args) -> Path:
    mkdirs()
    path = SKILL_DOC_DIR / "IGVF_KG_TRAVERSAL_SKILLS.md"
    lines = [
        "# Skill: IGVF Knowledge Graph traversal",
        "",
        "Iteratively walks the IGVF Catalog Knowledge Graph "
        "(`api.catalogkg.igvf.org` + the underlying ArangoDB at "
        "`db.catalog.igvf.org`) starting from a single entity and assembles "
        "a unified evidence pack across direct neighbors, second-degree "
        "relations, and optional cross-skill enrichment (FAVOR, "
        "enhancer-gene linkage, IGVF single-cell datasets, prior literature).",
        "",
        "Designed as the orchestrator-friendly **comprehensive context** "
        "tool: one CLI call → per-relation manifests + JSON evidence pack "
        "+ markdown report.",
        "",
        "## Subcommands",
        "",
        "### 1. `gene <symbol>` — comprehensive gene-centric traversal",
        "",
        "```bash",
        "python3 Scripts/kg_traversal_skill.py gene APOE \\",
        "    --depth 2 --limit 50 \\",
        "    --max-variants 25 --subvariant-limit 10 \\",
        "    --call-favor --call-linkage --call-singlecell --call-literature \\",
        "    --literature-context Alzheimer cardiovascular \\",
        "    --label apoe_full",
        "```",
        "",
        "Default direct relations (each saved to its own manifest CSV):",
        "",
        "- `variants` — variants on the gene",
        "- `coding_variant_scores` — MutPred2 / ESM-1v predictions per coding change",
        "- `regulatory_elements` — gene→cCRE links (CRISPRi / Perturb-seq / etc.)",
        "- `transcripts` — gene model isoforms",
        "- `proteins` — protein records linked to the gene",
        "- `diseases` — gene→disease/phenotype links",
        "- `pathways` — gene→pathway membership",
        "",
        "At `--depth 2`, each variant additionally fans out into its own "
        "summary, QTL genes, phenotypes, biosamples (CRISPRi / MPRA), "
        "genomic-element overlaps, and prediction sets.",
        "",
        "Optional side-calls:",
        "",
        "- `--call-favor` — pulls FAVOR functional annotations for the "
        "gene region.",
        "- `--call-linkage` — adds enhancer-gene linkage predictions for "
        "the gene region (rE2G / catalog regulatory-region links).",
        "- `--call-singlecell` — searches the IGVF Portal for single-cell "
        "AnalysisSets that mention the gene; surfaces candidate datasets "
        "for downstream expression analysis with "
        "`Scripts/single_cell_data_skills.py` or "
        "`Scripts/splitseq_pipeline.py`.",
        "- `--call-literature` — runs `Scripts/reference_skill.py "
        "validate` on the gene + your context terms.",
        "",
        "### 2. `variant <id>` — variant-centric traversal",
        "",
        "```bash",
        "python3 Scripts/kg_traversal_skill.py variant rs429358 \\",
        "    --call-favor --call-literature --label apoe_e4_variant",
        "```",
        "",
        "Variant ID accepted as rsID, SPDI, HGVS, or chr:pos:ref:alt where "
        "the Catalog API understands it.",
        "",
        "### 3. `region <chr:start-end>`",
        "",
        "```bash",
        "python3 Scripts/kg_traversal_skill.py region chr19:44903000-44912000 \\",
        "    --call-favor --label apoe_locus",
        "```",
        "",
        "Returns: genes overlapping the region, regulatory elements (cCREs) "
        "in the region, region-predictor enhancer-gene linkage rows, and "
        "(optional) FAVOR variant annotations.",
        "",
        "### 4. `aql` — direct ArangoDB AQL pass-through",
        "",
        "```bash",
        "python3 Scripts/kg_traversal_skill.py aql \\",
        "    'FOR g IN genes FILTER g.name == \"APOE\" RETURN g'",
        "```",
        "",
        "## Outputs",
        "",
        "Each run writes a timestamped folder under `Docs/KGTraversal/`:",
        "",
        "  Docs/KGTraversal/<timestamp>_<label>/",
        "    ├─ <entity>_<key>_report.md   # full markdown report",
        "    ├─ evidence_pack.json         # complete JSON with every relation",
        "    └─ Manifests/                  # one CSV per relation",
        "",
        "## How this chains with other skills",
        "",
        "- The `evidence_pack.json` is designed to be consumed by other "
        "skills — variant manifests can be fed to "
        "`Scripts/advanced_variant_analysis.py` or "
        "`Scripts/annotate_variant_list.py`; cCRE manifests can be fed to "
        "`Scripts/enhancer_gene_linkage_skills.py compare-sets`.",
        "- The single-cell dataset hits can be drilled into with "
        "`Scripts/single_cell_data_skills.py manifest` or "
        "`Scripts/splitseq_pipeline.py manifest`.",
        "- The literature manifest can be cross-checked with "
        "`Scripts/reference_skill.py validate`.",
        "- The internal orchestrator (Plan → Action → Results → Evaluation) "
        "uses this skill as the primary 'comprehensive context' action: "
        "given a gene of interest from the planning step, this skill "
        "supplies all the multi-omic evidence the analysis and evaluation "
        "steps need.",
    ]
    path.write_text("\n".join(lines))
    print(f"Playbook: {path}")
    return path


# --------------------------------- CLI ---------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(
        description="IGVF Knowledge Graph traversal: gene / variant / region "
                    "centric multi-hop evidence retrieval."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    g2 = sub.add_parser("genes", help="SEVERAL genes' regulatory evidence in "
                                       "ONE call — use this for multi-gene "
                                       "questions instead of looping `gene`.")
    g2.add_argument("symbols", help="Comma-separated, e.g. PAX2,LHX1,WT1")
    g2.add_argument("--tissue", help="Comma-separated include terms; a row "
                                      "matching ANY is kept (kidney,renal).")
    g2.add_argument("--exclude-tissue", default="adrenal",
                    help="Comma-separated terms to drop (default: adrenal, "
                         "which 'renal' otherwise matches).")
    g2.add_argument("--max-pages", type=int, default=40)
    g2.add_argument("--no-exhaustive", action="store_true")
    g2.add_argument("--label", default="")
    g2.set_defaults(func=cmd_genes)

    s = sub.add_parser("gene", help="Comprehensive gene-centric traversal.")
    s.add_argument("symbol")
    s.add_argument("--limit", type=int, default=25)
    s.add_argument("--depth", type=int, default=1,
                    help="1 = gene -> direct relations; 2 = also fan out from variants.")
    s.add_argument("--max-variants", type=int, default=25,
                    help="Cap on per-variant fan-out at depth>=2.")
    s.add_argument("--subvariant-limit", type=int, default=10)
    s.add_argument("--call-favor", action="store_true")
    s.add_argument("--favor-max", type=int, default=50)
    s.add_argument("--call-linkage", action="store_true")
    s.add_argument("--no-exhaustive-linkage", action="store_true",
                    help="Single capped request instead of paging every "
                         "page. Faster, but then absence cannot be "
                         "concluded from the result.")
    s.add_argument("--call-singlecell", action="store_true")
    s.add_argument("--call-literature", action="store_true")
    s.add_argument("--literature-context", nargs="*", default=None)
    s.add_argument("--literature-top", type=int, default=10)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_gene)

    s = sub.add_parser("variant", help="Variant-centric traversal.")
    s.add_argument("variant")
    s.add_argument("--limit", type=int, default=25)
    s.add_argument("--call-favor", action="store_true")
    s.add_argument("--call-literature", action="store_true")
    s.add_argument("--literature-context", nargs="*", default=None)
    s.add_argument("--literature-top", type=int, default=10)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_variant)

    s = sub.add_parser("region", help="Region-centric traversal.")
    s.add_argument("region", help="chr19:44903000-44912000")
    s.add_argument("--limit", type=int, default=50)
    s.add_argument("--call-favor", action="store_true")
    s.add_argument("--favor-max", type=int, default=100)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_region)

    s = sub.add_parser("aql", help="Direct AQL pass-through.")
    s.add_argument("query")
    s.add_argument("--limit", type=int, default=25)
    s.set_defaults(func=cmd_aql)

    s = sub.add_parser("write-playbook",
                        help="Write Docs/Skills/IGVF_KG_TRAVERSAL_SKILLS.md.")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
