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

def listify(data: Any) -> list[dict]:
    if isinstance(data, list):
        return [d for d in data if isinstance(d, dict)]
    if isinstance(data, dict):
        for k in ("results", "items", "data", "@graph", "edges"):
            if k in data and isinstance(data[k], list):
                return [d for d in data[k] if isinstance(d, dict)]
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
    "pathways":            ("/api/genes/pathways",                      {"gene_name": "{symbol}"}),
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


def fetch_gene_relations(symbol: str, limit: int) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for rel, (path, tpl) in GENE_RELATIONS.items():
        params = _format_params(tpl, symbol=symbol)
        params["limit"] = limit
        status, data = catalog_get(path, **params)
        if status == 200:
            rows = listify(data)
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

def fetch_linkage_for_region(region: str, limit: int = 25) -> dict[str, list[dict]]:
    """Pull IGVF Catalog enhancer-gene linkage evidence overlapping a region.

    Uses /api/regulatory-regions/genes (region predictor predictions)."""
    out: dict[str, list[dict]] = {"region_predictions": [],
                                    "qtl_links": []}
    status, data = catalog_get(
        "/api/regulatory-regions/genes", region=region, limit=limit)
    if status == 200:
        out["region_predictions"] = listify(data)
    else:
        # second pass: search via genomic-element predictions endpoint
        status2, data2 = catalog_get(
            "/api/genomic-elements/genes", region=region, limit=limit)
        if status2 == 200:
            out["region_predictions"] = listify(data2)
    return out


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
    for k in ("name", "symbol", "term_name", "label", "id", "_key", "_id",
              "spdi", "rsid", "hgvs", "variant_id", "uniprot_id", "accession"):
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

    if linkage and any(linkage.values()):
        lines += ["## Enhancer-gene linkage in this region",
                   "",
                   summarize_relation(linkage.get("region_predictions", [])),
                   ""]
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
        linkage = fetch_linkage_for_region(region, limit=args.limit)
        if linkage.get("region_predictions"):
            p = manifests_dir / "linkage_region_predictions.csv"
            write_csv(p, linkage["region_predictions"])
            manifest_paths["linkage_region_predictions"] = p

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

    # Linkage for region
    linkage = fetch_linkage_for_region(region, limit=args.limit)
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
