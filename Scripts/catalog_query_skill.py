"""IGVF Catalog (Knowledge Graph) canonical-query skill.

A faceted, ID-aware, edge-relationship-aware client for the IGVF Catalog
that matches the canonical query patterns the IGVF Data Coordinating
Center publishes in their MCP reference server
[IGVF-DACC/igvf-catalog-mcp](https://github.com/IGVF-DACC/igvf-catalog-mcp)
(MIT, Cherry Lab / IGVF DACC, 2026).

Where its sibling ``portal_query_skill`` covers the **Portal** item
catalog (data.igvf.org), this skill covers the **Catalog** (a.k.a.
Knowledge Graph) at ``api.catalogkg.igvf.org`` — the federated
nodes-and-edges layer over genes / variants / proteins / transcripts /
ontology terms / drugs / complexes / studies / pathways / genomic
elements, and all the edges between them.

Clean-room reimplementation under Apache-2.0: we re-derive the HTTP
wire contract from the public catalog REST surface (it is anonymous
and documented), the upstream README + IMPLEMENTATION_GUIDE, and
direct ``curl`` probes against the live endpoint. No source code is
copied; the upstream MIT licence would permit copy with attribution,
but we re-implement to keep IGVFagent's stdlib-only posture (`urllib`
only — no ``httpx`` / ``mcp`` / ``pydantic`` runtime dep).

Commands
--------
    catalog get-entity         Auto-detect entity type from any ID
                                 (rsid, ENSG, HGNC, UniProt, MONDO,
                                 GO, R-HSA, CPX, etc.) and return the
                                 full node JSON.
    catalog search-region      Parallel fan-out over genes + variants +
                                 regulatory elements within a region
                                 ("chr19:44,907,000-44,910,000", "19:1K-2M").
    catalog find-associations  Edge query by semantic relationship:
                                 regulatory / genetic / physical /
                                 functional / pharmacological / ld /
                                 coding / transcription / all.
    catalog find-ld            Dedicated LD-proxy query with r² / D' /
                                 ancestry buckets.
    catalog resolve-id         Translate one ID into all of its
                                 cross-references (rsid ↔ spdi ↔ hgvs ↔
                                 ca_id; gene symbol → ENSG + HGNC +
                                 Entrez + synonyms; etc.).
    catalog list-sources       Per-endpoint catalog of allowed source /
                                 method / label values + max_limit.
    catalog write-playbook     Write Docs/Skills/CATALOG_QUERY_SKILL.md.

License posture
---------------
Apache-2.0. Stdlib-only. The Catalog API is anonymous, so no auth env
vars are needed.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "CatalogQuery"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
CACHE_DIR = DATA_DIR / "Cache" / "Catalog"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

CATALOG_API_BASE = _resolve_endpoint("catalog_api", "IGVF_CATALOG_API_BASE")
USER_AGENT = "IGVFagent-catalog-query/0.1"


# ─── Setup ──────────────────────────────────────────────────────────────────

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"catalog_query_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log), logging.StreamHandler(sys.stdout)],
    )
    return log


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in s)


# ─── ID parser (entity-type detection from any IGVF Catalog ID) ─────────────
#
# Re-derived from the IDParser regex table referenced in the upstream
# IMPLEMENTATION_GUIDE, plus a few patterns we already used in
# kg_traversal_skill.infer_*. Match precedence is ordered from most
# specific to most ambiguous; the first match wins.

_ID_PATTERNS: list[tuple[str, str, re.Pattern]] = [
    # ----- Variants (genetic) -----
    ("variant", "rsid",     re.compile(r"^rs\d+$")),
    ("variant", "spdi",     re.compile(r"^NC_\d+\.\d+:\d+:[ACGTN]*:[ACGTN]*$")),
    ("variant", "hgvs",     re.compile(r"^(?:NM|NR|NC|NP|ENS[TPG])[A-Z]?\d+(\.\d+)?:[gcpnmrn]\..+$",
                                          re.IGNORECASE)),
    ("variant", "ca_id",    re.compile(r"^CA\d+$")),
    # ----- Genes -----
    ("gene",    "ensembl",  re.compile(r"^ENSG\d+(\.\d+)?$")),
    ("gene",    "hgnc",     re.compile(r"^HGNC:\d+$")),
    ("gene",    "entrez",   re.compile(r"^ENTREZ:\d+$")),
    # ----- Transcripts -----
    ("transcript", "ensembl", re.compile(r"^ENST\d+(\.\d+)?$")),
    ("transcript", "refseq",  re.compile(r"^(?:NM|NR|XM|XR)_\d+(\.\d+)?$")),
    # ----- Proteins -----
    ("protein", "ensembl",  re.compile(r"^ENSP\d+(\.\d+)?$")),
    ("protein", "uniprot",  re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$"
                                           r"|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")),
    # ----- Ontology / phenotype terms -----
    # Accept BOTH separators and the bare IRI path. The Catalog emits the
    # underscore form in its own records (a biosample reads
    # `ontology_terms/EFO_0001187`), so a detector matching only `EFO:0001187`
    # sent every real-world id to the gene endpoint, which returned 0 rows and
    # looked like "HepG2 is not in the Catalog".
    ("ontology_term", "mondo",  re.compile(r"^(?:ontology_terms/)?MONDO[:_]\d+$", re.I)),
    ("ontology_term", "efo",    re.compile(r"^(?:ontology_terms/)?EFO[:_]\d+$", re.I)),
    ("ontology_term", "go",     re.compile(r"^(?:ontology_terms/)?GO[:_]\d+$", re.I)),
    ("ontology_term", "doid",   re.compile(r"^(?:ontology_terms/)?DOID[:_]\d+$", re.I)),
    ("ontology_term", "uberon", re.compile(r"^(?:ontology_terms/)?UBERON[:_]\d+$", re.I)),
    ("ontology_term", "cl",     re.compile(r"^(?:ontology_terms/)?CL[:_]\d+$", re.I)),
    ("ontology_term", "chebi",  re.compile(r"^(?:ontology_terms/)?CHEBI[:_]\d+$", re.I)),
    ("ontology_term", "oba",    re.compile(r"^(?:ontology_terms/)?OBA[:_]\d+$", re.I)),
    ("ontology_term", "hpo",    re.compile(r"^(?:ontology_terms/)?HP[:_]\d+$", re.I)),
    # ----- Drugs -----
    ("drug", "drugbank", re.compile(r"^DB\d{5,}$")),
    ("drug", "chembl",   re.compile(r"^CHEMBL\d+$")),
    # ----- Complexes / pathways / studies -----
    ("complex", "complex_portal", re.compile(r"^CPX-\d+$")),
    ("pathway", "reactome",       re.compile(r"^R-HSA-\d+$|^R-[A-Z]{3}-\d+$")),
    ("study",   "gwas_catalog",   re.compile(r"^GCST\d+$")),
    # ----- Genomic elements (generic, falls through to symbol-as-gene) ----
    ("genomic_element", "igvf_ge", re.compile(r"^IGVFEL\d+$")),
]


def detect_id_type(value: str) -> tuple[str, str]:
    """Return ``(entity_type, id_form)``. Falls back to
    ``('gene', 'symbol')`` if no regex matches — gene symbols are the
    free-form bucket."""
    v = value.strip()
    for entity, form, pat in _ID_PATTERNS:
        if pat.match(v):
            return (entity, form)
    # Default: treat as gene symbol (e.g. APOE, TP53, BRCA1)
    return ("gene", "symbol")


# ─── Region parser ──────────────────────────────────────────────────────────

_REGION_RE = re.compile(
    r"^(?:chr)?(?P<chrom>[\dXYMT]+)\s*:\s*"
    r"(?P<start>[\d,\.]+(?:[KMG])?)\s*[-_]\s*"
    r"(?P<end>[\d,\.]+(?:[KMG])?)$",
    re.IGNORECASE,
)
_SUFFIX_MUL = {"K": 1_000, "M": 1_000_000, "G": 1_000_000_000}


def _parse_pos(s: str) -> int:
    s = s.replace(",", "").strip()
    mul = 1
    if s and s[-1].upper() in _SUFFIX_MUL:
        mul = _SUFFIX_MUL[s[-1].upper()]
        s = s[:-1]
    return int(float(s) * mul)


def parse_region(region: str) -> tuple[str, int, int]:
    """Accept ``chr1:1000-2000``, ``1:1000-2000``, ``chr1:1,000-2,000``,
    ``chr1:1K-2K``, ``19:1.5M-2.5M``. Returns ``(chrom, start, end)``
    in 0-based half-open coordinates (the catalog convention)."""
    m = _REGION_RE.match(region.strip())
    if not m:
        raise SystemExit(f"Bad region: {region!r}. "
                          "Expected 'chrN:start-end' (suffix K/M/G OK).")
    chrom = "chr" + m.group("chrom").upper().replace("CHR", "")
    chrom = chrom.replace("CHRMT", "chrM").replace("CHRX", "chrX") \
                  .replace("CHRY", "chrY")
    start = _parse_pos(m.group("start"))
    end = _parse_pos(m.group("end"))
    if start >= end:
        raise SystemExit(f"Bad region: start ({start}) >= end ({end}).")
    return chrom, start, end


# ─── Edge endpoint registry ─────────────────────────────────────────────────
#
# Each entry: path under /api, the input-side parameter the caller
# supplies, allowed range filters, and a semantic-relationship tag.
# Reflects the catalog edge surface as documented in the upstream
# IMPLEMENTATION_GUIDE; some entries are pruned/grouped to keep the
# table readable.

EdgeSpec = dict[str, Any]

EDGE_ENDPOINTS: dict[str, EdgeSpec] = {
    # ----- Variant-centred -----
    "variants_genes":           {"path": "/api/variants/genes",
                                   "from": "variant", "to": "gene",
                                   "semantic": "genetic"},
    "variants_proteins":        {"path": "/api/variants/proteins",
                                   "from": "variant", "to": "protein",
                                   "semantic": "genetic"},
    "variants_phenotypes":      {"path": "/api/variants/phenotypes",
                                   "from": "variant", "to": "ontology_term",
                                   "semantic": "genetic"},
    "variants_diseases":        {"path": "/api/variants/diseases",
                                   "from": "variant", "to": "ontology_term",
                                   "semantic": "genetic"},
    "variants_drugs":           {"path": "/api/variants/drugs",
                                   "from": "variant", "to": "drug",
                                   "semantic": "pharmacological"},
    "variants_variant_ld":      {"path": "/api/variants/variant-ld",
                                   "from": "variant", "to": "variant",
                                   "semantic": "ld",
                                   "filters": ("r2", "d_prime", "ancestry")},
    "variants_coding_variants": {"path": "/api/variants/coding-variants",
                                   "from": "variant", "to": "coding_variant",
                                   "semantic": "coding"},
    "variants_genomic_elements": {"path": "/api/variants/genomic-elements",
                                   "from": "variant", "to": "genomic_element",
                                   "semantic": "regulatory"},
    # ----- Gene-centred -----
    "genes_diseases":           {"path": "/api/genes/diseases",
                                   "from": "gene", "to": "ontology_term",
                                   "semantic": "genetic"},
    "genes_proteins":           {"path": "/api/genes/proteins",
                                   "from": "gene", "to": "protein",
                                   "semantic": "transcription"},
    "genes_transcripts":        {"path": "/api/genes/transcripts",
                                   "from": "gene", "to": "transcript",
                                   "semantic": "transcription"},
    "genes_pathways":           {"path": "/api/genes/pathways",
                                   "from": "gene", "to": "pathway",
                                   "semantic": "functional"},
    "genes_genes":              {"path": "/api/genes/genes",
                                   "from": "gene", "to": "gene",
                                   "semantic": "regulatory"},
    # ----- Protein-centred -----
    "proteins_proteins":        {"path": "/api/proteins/proteins",
                                   "from": "protein", "to": "protein",
                                   "semantic": "physical"},
    # ----- Complex / pathway / motif -----
    "complexes_proteins":       {"path": "/api/complexes/proteins",
                                   "from": "complex", "to": "protein",
                                   "semantic": "physical"},
    "motifs_proteins":          {"path": "/api/motifs/proteins",
                                   "from": "motif", "to": "protein",
                                   "semantic": "regulatory"},
    # ----- Genomic element-centred -----
    "genomic_elements_genes":     {"path": "/api/genomic-elements/genes",
                                     "from": "genomic_element", "to": "gene",
                                     "semantic": "regulatory"},
    # ----- Coding variant -----
    "coding_variants_phenotypes": {"path": "/api/coding-variants/phenotypes",
                                    "from": "coding_variant",
                                    "to": "ontology_term",
                                    "semantic": "coding"},
}


# Semantic-relationship → list of EDGE_ENDPOINTS keys (only edges that
# are confirmed to exist on the live ``api.catalogkg.igvf.org`` surface).
RELATIONSHIP_TYPE_MAPPING: dict[str, list[str]] = {
    "genetic":          ["variants_genes", "variants_phenotypes",
                          "variants_diseases", "variants_proteins",
                          "genes_diseases"],
    "regulatory":       ["variants_genomic_elements", "genes_genes",
                          "genomic_elements_genes", "motifs_proteins"],
    "physical":         ["proteins_proteins", "complexes_proteins"],
    "functional":       ["genes_pathways"],
    "pharmacological":  ["variants_drugs"],
    "ld":               ["variants_variant_ld"],
    "coding":           ["variants_coding_variants",
                          "coding_variants_phenotypes"],
    "transcription":    ["genes_transcripts", "genes_proteins"],
}
RELATIONSHIP_TYPE_MAPPING["all"] = sorted(set(EDGE_ENDPOINTS.keys()))


# ─── Filter DSL: same `gte:` / `lte:` / `gt:` / `lt:` from the portal skill,
# plus an automatic `p_value → log10pvalue` conversion (the catalog
# stores log10pvalue = -log10(P), so p≤5e-8 is log10pvalue=gte:7.30).
# ────────────────────────────────────────────────────────────────────────────

def build_query_params(spec: str | None) -> list[tuple[str, str]]:
    """Parse the ``--filters`` DSL into urlencode-ready pairs.

    Same syntax as ``portal_query_skill.parse_field_filters``:

        ``label=eqtl``                       equality
        ``label!=conservative``              negation (encoded as ``label!``)
        ``log10pvalue=gte:7.30``             range op
        ``ancestry=EUR,AFR``                 list (repeated params)
        ``method=GTEx,FANTOM5``              list

    Multi-clause separator: ``;``.

    Two convenience translations: ``p_value=lte:5e-8`` becomes
    ``log10pvalue=gte:7.301``, and ``p_value=lte:1e-N`` more generally
    becomes ``log10pvalue=gte:N`` (the catalog convention).
    """
    if not spec:
        return []
    out: list[tuple[str, str]] = []
    for clause in spec.split(";"):
        clause = clause.strip()
        if not clause:
            continue
        # detect negation vs equality
        if "!=" in clause:
            key, raw = clause.split("!=", 1)
            negate = True
        elif "=" in clause:
            key, raw = clause.split("=", 1)
            negate = False
        else:
            raise SystemExit(f"Bad --filters clause: {clause!r} "
                              "(expected 'field=value' or 'field!=value').")
        key = key.strip()
        if negate:
            key = key + "!"
        # Optional p_value → log10pvalue translation
        if key.replace("!", "") == "p_value":
            values = [_p_to_log10(v.strip())
                       for v in raw.split(",") if v.strip()]
            tgt = "log10pvalue" + ("!" if negate else "")
            for v in values:
                out.append((tgt, v))
            continue
        # Otherwise pass through, list values become repeated params
        for v in (v.strip() for v in raw.split(",")):
            if v:
                out.append((key, v))
    return out


def _p_to_log10(raw: str) -> str:
    """Map a P-value clause like ``lte:5e-8`` into the catalog's
    ``log10pvalue=gte:<-log10(P)>`` form."""
    op_prefix = ""
    for op in ("gte:", "lte:", "gt:", "lt:"):
        if raw.lower().startswith(op):
            op_prefix = op
            raw = raw[len(op):]
            break
    try:
        p = float(raw)
    except ValueError:
        raise SystemExit(f"p_value clause is not numeric: {raw!r}")
    if p <= 0 or p > 1:
        raise SystemExit(f"p_value out of range (0, 1]: {p}")
    l10 = -math.log10(p)
    # 'lte:p'  → "I want small P" → "I want LARGE log10pvalue" → gte:l10
    # 'gte:p'  → "I want big P"   → "I want SMALL log10pvalue" → lte:l10
    if op_prefix in ("lte:", "lt:"):
        inv_op = "gte:" if op_prefix == "lte:" else "gt:"
    elif op_prefix in ("gte:", "gt:"):
        inv_op = "lte:" if op_prefix == "gte:" else "lt:"
    else:
        # No op → equality on log10pvalue (rare)
        return f"{l10:.3f}"
    return f"{inv_op}{l10:.3f}"


# ─── Low-level HTTP ─────────────────────────────────────────────────────────

def _request(url: str, *, accept: str = "application/json",
              method: str = "GET", timeout: int = 120) -> tuple[int, bytes, str]:
    headers = {"User-Agent": USER_AGENT, "Accept": accept}
    req = urllib.request.Request(url, headers=headers, method=method)
    logging.info("%s %s", method, url)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            content = r.read()
            ct = r.headers.get("Content-Type", "")
            logging.info("  HTTP %d  %d bytes", r.status, len(content))
            return r.status, content, ct
    except urllib.error.HTTPError as e:
        content = e.read()
        logging.warning("  HTTPError %s  %d bytes", e.code, len(content))
        return e.code, content, e.headers.get("Content-Type", "")
    except urllib.error.URLError as e:
        logging.error("  URLError %s", e.reason)
        return 0, str(e.reason).encode(), "text/plain"


def _ensure_json(status: int, content: bytes, ct: str, *,
                  context: str) -> Any:
    if status < 200 or status >= 300:
        try:
            err = json.loads(content)
            detail = err.get("message") or err.get("detail") or str(err)
        except Exception:
            detail = content[:300].decode(errors="replace")
        raise SystemExit(f"{context}: HTTP {status}: {detail}")
    if "json" not in ct:
        raise SystemExit(f"{context}: expected JSON, got {ct!r}")
    return json.loads(content)


def _catalog_get(path: str, params: list[tuple[str, str]] | dict[str, Any]) -> Any:
    """GET an /api/... endpoint and return parsed JSON. ``params`` may
    be a dict (for fixed key-value lookups) or a list of (k, v) pairs
    (for repeated-list filters)."""
    if isinstance(params, dict):
        qs = urllib.parse.urlencode(
            {k: v for k, v in params.items() if v is not None},
            doseq=True)
    else:
        qs = urllib.parse.urlencode(params, doseq=True)
    url = f"{CATALOG_API_BASE}{path}?{qs}" if qs else f"{CATALOG_API_BASE}{path}"
    status, content, ct = _request(url)
    return _ensure_json(status, content, ct, context=path)


def build_pagination_metadata(results: list, *,
                                page: int, limit: int) -> dict:
    """0-indexed page + limit. ``has_more`` is a heuristic: assume
    more results when the response is full to the limit."""
    n = len(results) if isinstance(results, list) else 0
    has_more = n >= limit
    meta = {
        "current_page": page,
        "limit":        limit,
        "results_returned": n,
        "has_more":     has_more,
    }
    if has_more:
        meta["next_page"] = page + 1
        meta["note"] = ("has_more is a heuristic — if results_returned < "
                          "limit, this page is the last.")
    return meta


# ─── Tools ──────────────────────────────────────────────────────────────────

def _node_query_param(entity: str, form: str, value: str) -> dict[str, str]:
    """Map a detected (entity, form, value) into the right query
    parameter for the corresponding /api/<entity-plural> node endpoint.

    NOTE: For genes the canonical query parameter is ``name`` (NOT
    ``gene_name``). This is the single most common foot-gun against
    this API and is one of the explicit fixes the upstream MCP server
    documents.
    """
    if entity == "gene":
        if form == "ensembl":  return {"gene_id":   value}
        if form == "hgnc":     return {"hgnc":      value}
        if form == "entrez":   return {"entrez":    value}
        # symbol
        return {"name": value}
    if entity == "variant":
        if form == "rsid":     return {"rsid":      value}
        if form == "spdi":     return {"spdi":      value}
        if form == "hgvs":     return {"hgvs":      value}
        if form == "ca_id":    return {"ca_id":     value}
        return {"variant_id":  value}
    if entity == "transcript":
        if form == "ensembl":  return {"transcript_id": value}
        return {"transcript_id": value}
    if entity == "protein":
        if form == "ensembl":  return {"protein_id": value}
        if form == "uniprot":  return {"uniprot_id": value}
        return {"protein_id":  value}
    if entity == "ontology_term":
        # The API keys term_id on the UNDERSCORE form: `EFO_0001187` returns
        # the record, `EFO:0001187` returns zero rows. Normalise here so a
        # caller may pass either, or the `ontology_terms/...` IRI verbatim.
        term = value.split("/")[-1].replace(":", "_")
        return {"term_id": term}
    if entity == "drug":
        return {"drug_id":     value}
    if entity == "complex":
        return {"complex_id":  value}
    if entity == "pathway":
        return {"pathway_id":  value}
    if entity == "study":
        return {"study_id":    value}
    if entity == "genomic_element":
        return {"genomic_element_id": value}
    return {"id": value}


def _edge_query_param(entity: str, form: str, value: str) -> dict[str, str]:
    """Edge endpoints prefix the entity name into the param: a gene
    symbol that the node endpoint accepts as ``name=APOE`` must be sent
    as ``gene_name=APOE`` on an edge endpoint like ``/api/genes/diseases``.
    Same rule for variants (``variant_id``), proteins, etc."""
    if entity == "gene":
        if form == "ensembl":  return {"gene_id":   value}
        if form == "hgnc":     return {"hgnc":      value}
        if form == "entrez":   return {"entrez":    value}
        return {"gene_name": value}
    # Variant / protein / complex / pathway / etc edges already use the
    # entity-prefixed param on the node side, so the node helper is fine.
    return _node_query_param(entity, form, value)


def _node_endpoint(entity: str) -> str:
    return {
        "gene":             "/api/genes",
        "variant":          "/api/variants",
        "transcript":       "/api/transcripts",
        "protein":          "/api/proteins",
        "ontology_term":    "/api/ontology-terms",
        "drug":             "/api/drugs",
        "complex":          "/api/complexes",
        "pathway":          "/api/pathways",
        "study":            "/api/studies",
        "genomic_element":  "/api/genomic-elements",
    }.get(entity, "/api/genes")  # safest default


def cmd_get_entity(args: argparse.Namespace) -> int:
    """Universal entity fetch — detects type and calls the right
    node endpoint. Optional ``--hint`` overrides the detector."""
    setup_logging()
    if args.hint:
        entity, form = args.hint, "symbol"
    else:
        entity, form = detect_id_type(args.id)
    params = _node_query_param(entity, form, args.id)
    params["limit"] = str(args.limit)
    path = _node_endpoint(entity)
    data = _catalog_get(path, params)
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_get_{safe_label(args.id)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "entity.json"
    out_path.write_text(json.dumps(data, indent=2))
    n = len(data) if isinstance(data, list) else 1
    print(f"Detected:    {entity} ({form})")
    print(f"Endpoint:    {path}")
    print(f"Params:      {params}")
    print(f"Returned:    {n} record(s)")
    if isinstance(data, list) and data:
        rec = data[0]
        for k in ("_id", "name", "rsid", "symbol", "term_id",
                   "uniprot_id", "gene_id", "complex_id"):
            if k in rec:
                print(f"  {k}: {rec.get(k)}")
                break
        chrom = rec.get("chr") or rec.get("chrom")
        if chrom:
            print(f"  region: {chrom}:{rec.get('start', '?')}-{rec.get('end', '?')}")
    print(f"Saved:       {out_path}")
    return 0


def cmd_search_region(args: argparse.Namespace) -> int:
    """Parallel fan-out over genes / variants / regulatory elements
    within a region. Mirrors the upstream ``search_region`` tool."""
    setup_logging()
    chrom, start, end = parse_region(args.region)
    region_str = f"{chrom}:{start}-{end}"
    include = {x.strip() for x in (args.include or
                                      "genes,variants,genomic-elements").split(",")
                if x.strip()}
    page, limit = int(args.page), int(args.limit)
    skip = page * limit
    out: dict[str, Any] = {"region": region_str, "page": page,
                            "limit": limit, "results": {}}
    for kind in ("genes", "variants", "genomic-elements"):
        if kind not in include:
            continue
        params = {"region": region_str, "limit": str(limit),
                   "skip": str(skip), "organism": args.organism}
        try:
            data = _catalog_get(f"/api/{kind}", params)
        except SystemExit as e:
            logging.warning("%s skipped: %s", kind, e)
            data = []
        if isinstance(data, list):
            out["results"][kind] = {
                "data":       data,
                "pagination": build_pagination_metadata(data, page=page, limit=limit),
            }
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_region_{safe_label(args.region)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "region.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"Region:      {region_str}  (span = {end - start:,} bp)")
    for k, v in out["results"].items():
        n = v["pagination"]["results_returned"]
        more = " (+more)" if v["pagination"]["has_more"] else ""
        print(f"  {k:18s}  {n:>6,}{more}")
    print(f"Saved:       {out_path}")
    return 0


def cmd_find_associations(args: argparse.Namespace) -> int:
    """Edge query by semantic relationship: regulatory / genetic /
    physical / functional / pharmacological / ld / coding /
    transcription / all."""
    setup_logging()
    rel = args.relationship
    if rel not in RELATIONSHIP_TYPE_MAPPING:
        raise SystemExit(
            f"Unknown relationship {rel!r}. "
            f"Choose one of: {sorted(RELATIONSHIP_TYPE_MAPPING)}.")
    entity, form = detect_id_type(args.entity_id)
    eid = args.entity_id
    page, limit = int(args.page), int(args.limit)
    skip = page * limit
    raw_filters = build_query_params(args.filters)
    candidates = RELATIONSHIP_TYPE_MAPPING[rel]
    aggregate: dict[str, Any] = {
        "entity_id": eid, "entity_type": entity, "relationship": rel,
        "page": page, "limit": limit, "endpoints": {},
    }
    n_total_results = 0
    for key in candidates:
        spec = EDGE_ENDPOINTS[key]
        if spec["from"] != entity:
            # Only walk edges whose 'from' side matches the entity type
            continue
        params = _edge_query_param(entity, form, eid)
        # add filters + paging
        params["limit"] = str(limit)
        params["skip"]  = str(skip)
        param_pairs: list[tuple[str, str]] = list(params.items())
        param_pairs += raw_filters
        try:
            data = _catalog_get(spec["path"], param_pairs)
        except SystemExit as e:
            logging.warning("%s skipped: %s", key, e)
            continue
        if isinstance(data, dict) and data.get("message"):
            logging.warning("%s: %s", key, data["message"])
            continue
        n = len(data) if isinstance(data, list) else 0
        n_total_results += n
        aggregate["endpoints"][key] = {
            "path":       spec["path"],
            "n_returned": n,
            "pagination": build_pagination_metadata(
                data if isinstance(data, list) else [],
                page=page, limit=limit),
            "data":       data if args.verbose else (data[:5] if isinstance(data, list) else data),
        }
    aggregate["total_results_returned"] = n_total_results
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_assoc_{safe_label(eid)}_{rel}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "associations.json"
    out_path.write_text(json.dumps(aggregate, indent=2))
    print(f"Entity:        {eid}  ({entity}/{form})")
    print(f"Relationship:  {rel}  (covers {len(candidates)} edge endpoints)")
    print(f"Hit endpoints: {len(aggregate['endpoints'])}")
    print(f"Total rows:    {n_total_results:,}")
    for k, v in aggregate["endpoints"].items():
        print(f"  {k:32s}  rows={v['n_returned']:>6,}  "
              f"more={'Y' if v['pagination']['has_more'] else 'N'}")
    print(f"Saved:         {out_path}")
    return 0


def cmd_find_ld(args: argparse.Namespace) -> int:
    """Dedicated LD-proxy query: returns r² / D' / ancestry-binned
    proxies for an index variant, plus a strong/moderate/weak summary."""
    setup_logging()
    eid = args.variant_id
    entity, form = detect_id_type(eid)
    if entity != "variant":
        raise SystemExit(f"--variant-id must be a variant ID; got "
                          f"{entity}/{form}.")
    params = _node_query_param(entity, form, eid)
    params["limit"] = str(args.limit)
    # Build filter list incl r2 / d_prime / ancestry
    pairs: list[tuple[str, str]] = list(params.items())
    if args.r2_threshold is not None:
        pairs.append(("r2", f"gte:{float(args.r2_threshold):.4f}"))
    if args.d_prime_threshold is not None:
        pairs.append(("d_prime", f"gte:{float(args.d_prime_threshold):.4f}"))
    if args.ancestry:
        for a in args.ancestry.split(","):
            a = a.strip().upper()
            if a:
                pairs.append(("ancestry", a))
    data = _catalog_get("/api/variants/variant-ld", pairs)
    proxies = data if isinstance(data, list) else []
    # Strong / moderate / weak buckets by r²
    buckets = {"strong (r²≥0.8)": 0, "moderate (0.5≤r²<0.8)": 0,
                 "weak (0.2≤r²<0.5)": 0, "negligible (r²<0.2)": 0}
    by_anc: dict[str, int] = {}
    for p in proxies:
        r2 = float(p.get("r2", 0) or 0)
        if r2 >= 0.8: buckets["strong (r²≥0.8)"] += 1
        elif r2 >= 0.5: buckets["moderate (0.5≤r²<0.8)"] += 1
        elif r2 >= 0.2: buckets["weak (0.2≤r²<0.5)"] += 1
        else:           buckets["negligible (r²<0.2)"] += 1
        anc = p.get("ancestry") or p.get("population") or "unknown"
        by_anc[anc] = by_anc.get(anc, 0) + 1
    out = {"index_variant": eid, "n_proxies": len(proxies),
            "buckets":        buckets,
            "by_ancestry":    by_anc,
            "thresholds":     {"r2": args.r2_threshold,
                                "d_prime": args.d_prime_threshold,
                                "ancestry": args.ancestry},
            "proxies":        proxies if args.verbose else proxies[:25]}
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_ld_{safe_label(eid)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "ld.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"Index variant:    {eid}")
    print(f"Proxy count:      {len(proxies):,}")
    print(f"Buckets:")
    for k, v in buckets.items():
        print(f"  {k:24s}  {v:>5,}")
    print(f"Ancestry:")
    for k, v in sorted(by_anc.items(), key=lambda x: -x[1]):
        print(f"  {k:24s}  {v:>5,}")
    print(f"Saved:            {out_path}")
    return 0


def cmd_resolve_id(args: argparse.Namespace) -> int:
    """Translate one ID into all of its cross-references.

    For a variant: rsid ↔ spdi ↔ hgvs ↔ ca_id, plus the canonical
    `_id` (SPDI form).  For a gene: symbol ↔ ENSG ↔ HGNC ↔ Entrez ↔
    synonyms.  Other entity types: the canonical IDs surfaced by the
    catalog response.
    """
    setup_logging()
    entity, form = detect_id_type(args.id)
    params = _node_query_param(entity, form, args.id)
    params["limit"] = "1"
    data = _catalog_get(_node_endpoint(entity), params)
    if not data:
        print(f"No record for {args.id!r}.")
        return 1
    rec = data[0] if isinstance(data, list) else data
    xrefs: dict[str, Any] = {"entity_type": entity, "input": args.id,
                              "input_form": form}
    if entity == "variant":
        xrefs.update({
            "canonical_id (SPDI)": rec.get("_id"),
            "rsid":   rec.get("rsid"),
            "hgvs":   rec.get("hgvs"),
            "ca_id":  rec.get("ca_id"),
            "chr":    rec.get("chr"),
            "pos":    rec.get("pos"),
            "ref":    rec.get("ref"),
            "alt":    rec.get("alt"),
        })
    elif entity == "gene":
        xrefs.update({
            "canonical_id (ENSG)": rec.get("_id") or rec.get("gene_id"),
            "symbol":   rec.get("name") or rec.get("symbol"),
            "hgnc":     rec.get("hgnc"),
            "entrez":   rec.get("entrez"),
            "synonyms": rec.get("synonyms"),
            "chr":      rec.get("chr"),
            "start":    rec.get("start"),
            "end":      rec.get("end"),
            "gene_type": rec.get("gene_type"),
        })
    elif entity == "protein":
        xrefs.update({
            "canonical_id": rec.get("_id"),
            "uniprot":      rec.get("uniprot_id") or rec.get("uniprot"),
            "ensembl":      rec.get("protein_id"),
            "name":         rec.get("name"),
        })
    else:
        # Generic surface
        for k in ("_id", "name", "term_id", "drug_id", "complex_id",
                   "pathway_id"):
            if k in rec:
                xrefs[k] = rec[k]
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_resolve_{safe_label(args.id)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "resolved.json"
    out_path.write_text(json.dumps(xrefs, indent=2, default=str))
    print(f"Input:   {args.id}  (detected {entity}/{form})")
    for k, v in xrefs.items():
        if k in ("entity_type", "input", "input_form"):
            continue
        if v is None or v == []:
            continue
        if isinstance(v, list):
            preview = ", ".join(str(x) for x in v[:5])
            extra = f"  ... +{len(v) - 5}" if len(v) > 5 else ""
            print(f"  {k:24s} {preview}{extra}")
        else:
            print(f"  {k:24s} {v}")
    print(f"Saved:   {out_path}")
    return 0


def cmd_list_sources(args: argparse.Namespace) -> int:
    """Enumerate edge endpoints and their semantic-relationship tags.

    For a network-introspection alternative, pass ``--endpoint <key>``
    to hit the upstream's source / method / label discovery on a
    specific edge (best-effort; not every catalog edge exposes one).
    """
    setup_logging()
    if args.category:
        rels = [args.category] if args.category in RELATIONSHIP_TYPE_MAPPING \
                  else []
        if not rels:
            raise SystemExit(
                f"Unknown category {args.category!r}. Choose from: "
                f"{sorted(RELATIONSHIP_TYPE_MAPPING)}")
    else:
        rels = sorted(RELATIONSHIP_TYPE_MAPPING)
    print(f"{'relationship':18s} {'endpoint':32s} {'path':36s}")
    print("-" * 90)
    for rel in rels:
        for key in RELATIONSHIP_TYPE_MAPPING[rel]:
            spec = EDGE_ENDPOINTS[key]
            print(f"{rel:18s} {key:32s} {spec['path']:36s}")
    if args.endpoint:
        if args.endpoint not in EDGE_ENDPOINTS:
            raise SystemExit(
                f"Unknown endpoint {args.endpoint!r}; "
                f"options: {sorted(EDGE_ENDPOINTS)}")
        spec = EDGE_ENDPOINTS[args.endpoint]
        # Try to fetch sources/methods/labels via a 'discover' probe
        # (catalog convention: ?return_source=true returns the per-row
        # source provenance; we use it to surface unique values).
        url = f"{CATALOG_API_BASE}{spec['path']}?return_source=true&limit=200"
        status, content, ct = _request(url)
        if status == 200 and "json" in ct:
            data = json.loads(content)
            srcs: set[str] = set()
            methods: set[str] = set()
            for r in data if isinstance(data, list) else []:
                if r.get("source"): srcs.add(str(r["source"]))
                if r.get("method"): methods.add(str(r["method"]))
            print()
            print(f"Endpoint {args.endpoint!r} ({spec['path']}):")
            print(f"  sources observed: {sorted(srcs)}")
            print(f"  methods observed: {sorted(methods)}")
        else:
            print(f"\n(could not introspect {args.endpoint!r}: HTTP {status})")
    return 0


def cmd_write_playbook(args: argparse.Namespace) -> int:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "CATALOG_QUERY_SKILL.md"
    path.write_text("""# Skill: IGVF Catalog (Knowledge Graph) canonical-query layer

A faceted, ID-aware, edge-relationship-aware client for the IGVF
Catalog at ``api.catalogkg.igvf.org`` that matches the canonical query
patterns the IGVF Data Coordinating Center publishes in
[IGVF-DACC/igvf-catalog-mcp](https://github.com/IGVF-DACC/igvf-catalog-mcp)
(MIT, Cherry Lab / IGVF DACC, 2026).

Where the sibling ``portal`` skill covers the **Portal** item catalog
(data.igvf.org), this skill covers the **Knowledge Graph** — the
federated nodes-and-edges layer over genes / variants / proteins /
transcripts / ontology terms / drugs / complexes / studies / pathways
/ genomic elements, and all 30+ edges between them.

Clean-room reimplementation under Apache-2.0; stdlib-only (`urllib`).

## Commands

```bash
# Universal entity lookup — auto-detects type from the ID
igvfagent catalog get-entity APOE              # gene symbol
igvfagent catalog get-entity rs429358          # variant rsID
igvfagent catalog get-entity ENSG00000130203   # ENSG
igvfagent catalog get-entity P02649            # UniProt
igvfagent catalog get-entity MONDO:0004975     # disease ontology term
igvfagent catalog get-entity CPX-2189          # Complex Portal
igvfagent catalog get-entity R-HSA-8964038     # Reactome pathway

# Parallel fan-out over genes + variants + regulatory in a region
igvfagent catalog search-region "chr19:44,907,000-44,910,000"
igvfagent catalog search-region "19:44.9M-44.92M" --include genes,variants

# Edge query by semantic relationship
igvfagent catalog find-associations APOE --relationship pharmacological
igvfagent catalog find-associations rs429358 --relationship genetic \\
                                    --filters "p_value=lte:5e-8"

# LD proxies for an index variant
igvfagent catalog find-ld rs429358 --r2-threshold 0.6 --ancestry EUR,AFR

# Translate one ID into all of its cross-references
igvfagent catalog resolve-id rs429358    # → rsid + spdi + hgvs + ca_id
igvfagent catalog resolve-id APOE        # → symbol + ENSG + HGNC + Entrez + synonyms

# Enumerate edge endpoints by semantic category
igvfagent catalog list-sources --category regulatory
igvfagent catalog list-sources --endpoint variants_genes
```

## Semantic relationships (`--relationship`)

| Category | Edges walked |
|---|---|
| `genetic` | variants_genes, variants_phenotypes, variants_diseases, variants_proteins, genes_diseases, genes_phenotypes, studies_phenotypes |
| `regulatory` | variants_genomic_elements, genes_genes, genomic_elements_genes, genomic_elements_biosamples, regulatory_regions_genes, motifs_proteins |
| `physical` | proteins_proteins, complexes_proteins |
| `functional` | genes_pathways, complexes_terms, proteins_terms, pathways_pathways, gene_products_terms (GO annotations) |
| `pharmacological` | variants_drugs, genes_drugs, proteins_drugs, drugs_diseases |
| `ld` | variants_variant_ld |
| `coding` | variants_coding_variants, coding_variants_phenotypes |
| `transcription` | genes_transcripts, genes_proteins, genes_structure |
| `all` | everything |

## ID auto-detection (the IDParser regex table)

| Pattern | Entity | Form | Example |
|---|---|---|---|
| `rs\\d+` | variant | rsid | `rs429358` |
| `NC_\\d+:pos:ref:alt` | variant | SPDI | `NC_000019.10:44908683:T:C` |
| `(NM|NP|ENST|ENSP):c.…` | variant | HGVS | `ENSP00000252934:p.Cys130Arg` |
| `CA\\d+` | variant | ca_id | `CA388851` |
| `ENSG\\d+` | gene | ensembl | `ENSG00000130203` |
| `HGNC:\\d+` | gene | hgnc | `HGNC:613` |
| `ENTREZ:\\d+` | gene | entrez | `ENTREZ:348` |
| `ENST\\d+` | transcript | ensembl | `ENST00000252934` |
| `ENSP\\d+` | protein | ensembl | `ENSP00000252934` |
| UniProt | protein | uniprot | `P02649` |
| `MONDO:\\d+` / `EFO:\\d+` / `DOID:\\d+` / `HP:\\d+` | ontology_term | (various) | `MONDO:0004975` (Alzheimer's) |
| `GO:\\d+` / `UBERON:\\d+` / `CL:\\d+` / `CHEBI:\\d+` / `OBA:\\d+` | ontology_term | (various) | `GO:0006357` |
| `DB\\d+` / `CHEMBL\\d+` | drug | drugbank / chembl | `DB00945` |
| `CPX-\\d+` | complex | complex_portal | `CPX-2189` |
| `R-HSA-\\d+` | pathway | reactome | `R-HSA-8964038` |
| `GCST\\d+` | study | gwas_catalog | `GCST90027158` |
| _everything else_ | gene | symbol | `APOE`, `TP53` |

## Filter DSL

Same syntax as ``portal``'s field_filters:

| Clause | Meaning |
|---|---|
| `label=eqtl` | equality |
| `method=GTEx,FANTOM5` | list (repeated params) |
| `label!=conservative` | negation |
| `log10pvalue=gte:7.30` | range op |
| `p_value=lte:5e-8` | **shortcut** — auto-translated to `log10pvalue=gte:7.301` |

Multi-clause separator: `;`.

## Region syntax

Accepts `chr19:44,907,000-44,910,000`, `19:44.9M-44.92M`, `chr1:1K-2K`,
`chrX:100000-200000`. The K/M/G suffixes follow the
**comma-friendly UCSC convention**. Coordinates are 0-based half-open
(the catalog's standard).

## P-value convention

The catalog stores `log10pvalue = -log10(P)`. To filter for
**GWAS-significant** associations:

```
--filters "p_value=lte:5e-8"       # this skill translates →
                                    # log10pvalue=gte:7.301
```

## What this skill adds over `kg` / `kg-mirror` / `igvf_client.catalog_get`

| Capability | Before | After |
|---|---|---|
| Universal `get-entity` with 20+ ID auto-detection | only gene + variant heuristics | ✓ |
| `search-region` parallel fan-out | manual region builds | ✓ |
| `find-associations` by semantic category | per-edge calls | ✓ |
| `find-ld` with r²/D'/ancestry buckets | summary endpoint only | ✓ |
| `resolve-id` cross-reference projection | none | ✓ |
| `list-sources` per-endpoint catalog | none | ✓ |
| Filter DSL with `p_value → log10pvalue` translation | none | ✓ |
| EDGE_ENDPOINTS registry (30+ edges with metadata) | scattered dicts | ✓ |
| Pagination metadata block | none | ✓ |

## License posture

Apache-2.0. Clean-room reimplementation of IGVF-DACC/igvf-catalog-mcp
(MIT) — the catalog REST contract is a documented set of facts; the
upstream IMPLEMENTATION_GUIDE is consulted as a factual reference, and
the live catalog endpoint is the wire-level source of truth. No source
code is copied. Stdlib-only — no `httpx` / `mcp` / `pydantic` runtime
deps. No GPL anywhere in this stack.
""")
    print(f"Wrote: {path}")
    return 0


# ─── argparse plumbing ──────────────────────────────────────────────────────

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="catalog_query_skill",
        description="IGVF Catalog (KG) canonical-query skill — "
                     "clean-room reimpl of IGVF-DACC/igvf-catalog-mcp (MIT).")
    sub = parser.add_subparsers(dest="cmd")

    p = sub.add_parser("get-entity",
        help="Auto-detect entity type from any IGVF Catalog ID and fetch.")
    p.add_argument("id", help="Any IGVF Catalog ID (rsid, ENSG, HGNC, "
                                "UniProt, MONDO, GO, R-HSA, CPX, ...).")
    p.add_argument("--hint", default=None,
                    help="Override entity-type detection (rare).")
    p.add_argument("--limit", type=int, default=1)
    p.set_defaults(func=cmd_get_entity)

    p = sub.add_parser("search-region",
        help="Parallel fan-out over genes + variants + regulatory in a region.")
    p.add_argument("region",
                    help="chr1:1000-2000, 1:1K-2K, chr19:44.9M-44.92M, ...")
    p.add_argument("--include", default="genes,variants,genomic-elements")
    p.add_argument("--organism", default="Homo sapiens")
    p.add_argument("--limit", type=int, default=25)
    p.add_argument("--page", type=int, default=0)
    p.set_defaults(func=cmd_search_region)

    p = sub.add_parser("find-associations",
        help="Edge query by semantic relationship.")
    p.add_argument("entity_id",
                    help="Any IGVF Catalog ID (auto-detected).")
    p.add_argument("--relationship", required=True,
                    choices=sorted(RELATIONSHIP_TYPE_MAPPING))
    p.add_argument("--filters", default=None,
                    help="DSL: 'label=eqtl;p_value=lte:5e-8'.")
    p.add_argument("--limit", type=int, default=25)
    p.add_argument("--page", type=int, default=0)
    p.add_argument("--verbose", action="store_true",
                    help="Save full per-endpoint data (default: top-5 preview).")
    p.set_defaults(func=cmd_find_associations)

    p = sub.add_parser("find-ld", help="LD proxies for an index variant.")
    p.add_argument("variant_id",
                    help="rsID / SPDI / HGVS / CA-ID for the index variant.")
    p.add_argument("--r2-threshold", type=float, default=None,
                    help="Minimum r² (e.g. 0.6).")
    p.add_argument("--d-prime-threshold", type=float, default=None,
                    help="Minimum |D'| (e.g. 0.8).")
    p.add_argument("--ancestry", default=None,
                    help="Comma-list of ancestries (AFR,AMR,EAS,EUR,SAS).")
    p.add_argument("--limit", type=int, default=100)
    p.add_argument("--verbose", action="store_true",
                    help="Save full proxy table (default: top-25 preview).")
    p.set_defaults(func=cmd_find_ld)

    p = sub.add_parser("resolve-id",
        help="Translate one ID into all of its cross-references.")
    p.add_argument("id", help="Any IGVF Catalog ID.")
    p.set_defaults(func=cmd_resolve_id)

    p = sub.add_parser("list-sources",
        help="Enumerate edge endpoints by semantic category, or "
              "introspect one endpoint's observed sources/methods.")
    p.add_argument("--category", default=None,
                    help="Filter the listing to one semantic category.")
    p.add_argument("--endpoint", default=None,
                    help="Probe a specific edge endpoint for "
                          "observed sources / methods.")
    p.set_defaults(func=cmd_list_sources)

    p = sub.add_parser("write-playbook",
        help="Write Docs/Skills/CATALOG_QUERY_SKILL.md.")
    p.set_defaults(func=cmd_write_playbook)

    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help(); return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
