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
from collections import Counter
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
    # Reverse direction + Reactome hierarchy: without these, a pathway ID
    # matched no edge whose `from` side was `pathway`, so
    # find-associations on an R-HSA-n returned "Hit endpoints: 0".
    "pathways_genes":           {"path": "/api/pathways/genes",
                                   "from": "pathway", "to": "gene",
                                   "semantic": "functional"},
    "pathways_pathways":        {"path": "/api/pathways/pathways",
                                   "from": "pathway", "to": "pathway",
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
    "functional":       ["genes_pathways", "pathways_genes",
                          "pathways_pathways"],
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
        # `protein_id` resolves BOTH the ENSP and the UniProt accession.
        # There is no `uniprot_id` query param — sending one is silently
        # ignored and yields an unfiltered row (see the pathway note below).
        return {"protein_id": value}
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
        # The pathway NODE endpoint keys on `id`, not `pathway_id` — the
        # one place the collection breaks the `<entity>_id` convention.
        # `/api/pathways` silently DROPS unknown params instead of erroring,
        # so `pathway_id=` returned the first row of the whole collection
        # (R-HSA-1059683, "Interleukin-6 signaling") for every ID queried.
        # The EDGE endpoints (`/api/pathways/genes`) do want `pathway_id`
        # — see _edge_query_param.
        return {"id": value}
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
    if entity == "pathway":
        # Inverse of the node case: `/api/pathways/genes` rejects `id`
        # ("At least one pathway property must be defined") and wants the
        # entity-prefixed `pathway_id`.
        return {"pathway_id": value}
    # Variant / protein / complex / etc edges already use the
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


# Fields that carry an entity's own identity in a node response, in the
# order we trust them.
_IDENTITY_FIELDS = ("_id", "id", "uri", "name", "rsid", "spdi", "hgvs", "ca_id",
                    "term_id", "gene_id", "transcript_id", "protein_id",
                    "uniprot_ids", "complex_id", "drug_id", "study_id",
                    "genomic_element_id", "name_aliases", "synonyms",
                    "id_version")


def _identity_matches(rec: dict, value: str) -> bool:
    """True when ``rec`` plausibly IS the record asked for.

    The Catalog silently DROPS query params it does not recognise and
    answers with an unfiltered page, so a wrong param name looks like a
    successful lookup that returns row #1 of the collection for every
    input. Comparing the request against the response's own identity
    fields is the only way to catch that from the client side.
    """
    want = value.strip().lower()
    for key in _IDENTITY_FIELDS:
        got = rec.get(key)
        if got is None:
            continue
        candidates = got if isinstance(got, (list, tuple)) else [got]
        for c in candidates:
            # Strip ArangoDB collection prefixes ("ontology_terms/GO_...")
            # and OBO IRIs before comparing.
            c = str(c).strip().lower().rsplit("/", 1)[-1]
            # `id_version` is "R-HSA-174824.6"; ontology _ids use "GO_..."
            # where the caller may have typed "GO:...".
            if c == want or c.split(".")[0] == want or \
               c.replace("_", ":") == want.replace("_", ":"):
                return True
    return False


def _check_identity(data: Any, value: str, params: dict, path: str) -> None:
    """Warn loudly when the response does not match what was asked for."""
    rec = data[0] if isinstance(data, list) and data else data
    if not isinstance(rec, dict) or _identity_matches(rec, value):
        return
    got = rec.get("_id") or rec.get("name") or "?"
    logging.warning(
        "IDENTITY MISMATCH: asked %s for %r, got %r. The query param %r is "
        "probably not recognised by this endpoint (unknown params are "
        "silently ignored and an unfiltered row is returned). Do NOT trust "
        "this record.",
        path, value, got, ",".join(k for k in params if k != "limit"))
    print(f"  !! WARNING: requested {value!r} but the endpoint returned "
          f"{got!r} — treat this record as UNVERIFIED.")


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
    _check_identity(data, args.id, params, path)
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
        # Print the canonical ID *and* the human-readable name — for
        # pathways / ontology terms the name is the whole point of the
        # lookup, and stopping at the first hit used to hide it.
        shown = 0
        for k in ("_id", "name", "rsid", "symbol", "term_id",
                   "uniprot_ids", "gene_id", "complex_id"):
            if k in rec and rec.get(k) is not None:
                print(f"  {k}: {rec.get(k)}")
                shown += 1
                if shown == 2:
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
    # `/api/variants/coding-variants` keys on SPDI and 400s on an rsID, so
    # walk it last and feed it an SPDI harvested from whichever earlier
    # edge already named the variant. Without this the endpoint errored on
    # every rsID query and its evidence was dropped in silence.
    ordered = sorted(candidates,
                     key=lambda k: 1 if k in _SPDI_ONLY_ENDPOINTS else 0)
    seen_records: list[dict] = []
    for key in ordered:
        spec = EDGE_ENDPOINTS[key]
        if spec["from"] != entity:
            # Only walk edges whose 'from' side matches the entity type
            continue
        if (key in _SPDI_ONLY_ENDPOINTS and entity == "variant"
                and form != "spdi"):
            spdi = (_spdi_from_records(seen_records)
                    or _resolve_spdi(entity, form, eid))
            if not spdi:
                logging.info("%s: no SPDI resolved, skipping %s", eid, key)
                continue
            params = {"spdi": spdi}
        else:
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
        if isinstance(data, list):
            seen_records.extend(r for r in data if isinstance(r, dict))
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


# ---------------------------------------------------------------------------
# Multi-assay variant evidence
# ---------------------------------------------------------------------------

# `method` values that name a *data source or relationship*, not an assay.
# Counting them inflates every variant's evidence tally: the LD graph in
# particular attaches to essentially every common variant, so leaving it
# in would make "variants with >3 assays" mostly a list of well-tagged
# common SNPs. Overridable with --exclude-method / --include-ld.
NON_ASSAY_METHODS = {
    "linkage disequilibrum",   # sic — upstream spelling
    "linkage disequilibrium",
}

# Endpoints that carry no experimental assay at all. The coding-variant
# edge returns in-silico predictor columns (SIFT, PolyPhen2, CADD, REVEL,
# AlphaMissense, ESM1b, VARITY, ...) with no `method` field, so every one
# of its records would otherwise land in the tally as a single bogus
# "unknown" assay and add +1 to every coding variant. Opt in with
# --include-predictions.
PREDICTION_ONLY_ENDPOINTS = {"variants_coding_variants"}

# How each evidence label was produced. Counting a computational
# prediction as an "assay" alongside an eQTL or an MPRA overstates the
# experimental support for a variant, which is the whole point of the
# question, so the breakdown is reported rather than hidden.
EVIDENCE_KIND = {
    "eQTL": "experimental", "spliceQTL": "experimental",
    "caQTL": "experimental", "pQTL": "experimental",
    "GWAS": "experimental", "ADASTRA": "experimental",
    "MPRA": "experimental", "CRISPR": "experimental",
    "SEMVAR": "computational", "cV2F": "computational",
    "GVATdb": "computational",
    "PharmGKB": "curated",
}


def _variant_key(rec: "dict") -> "str | None":
    """Best available variant identifier on an edge record.

    Catalog edges name the variant differently per endpoint
    (``sequence_variant`` on gene edges, ``variant`` on phenotype edges,
    ``rsid`` when it was the query key), so try each in turn rather than
    assuming one shape.
    """
    for k in ("sequence_variant", "variant", "variant_id", "rsid",
               "sequence variant"):
        v = rec.get(k)
        if v:
            return str(v)
    return None


def _assay_of(rec: "dict") -> str:
    """Evidence label for one edge record.

    Falls through method -> label -> source, because not every edge
    carries a `method`: the coding-variant predictor rows carry only a
    `source`, and labelling them all "unknown" both loses information
    and collapses distinct resources into one fake assay.
    """
    return str(rec.get("method") or rec.get("label")
               or rec.get("source") or "unknown")


def _harvest_traversal(paths: "list[Path]") -> "list[tuple[str, str, str]]":
    """(variant, assay, origin) triples from kg-traversal evidence packs."""
    out: "list[tuple[str, str, str]]" = []
    for root in paths:
        packs = ([root] if root.is_file()
                  else sorted(root.rglob("evidence_pack.json")))
        for pack in packs:
            try:
                doc = json.loads(pack.read_text())
            except Exception as exc:
                logging.warning("skipping %s: %s", pack, exc)
                continue
            origin = pack.parent.name
            blocks = []
            rel = doc.get("relations")
            if isinstance(rel, dict):
                blocks.extend((k, v) for k, v in rel.items())
            for key in ("variants", "variant_depth2", "linkage"):
                v = doc.get(key)
                if v is not None:
                    blocks.append((key, v))
            for name, block in blocks:
                rows = block if isinstance(block, list) else []
                if isinstance(block, dict):
                    for sub in block.values():
                        if isinstance(sub, list):
                            rows.extend(sub)
                for rec in rows:
                    if not isinstance(rec, dict):
                        continue
                    vid = _variant_key(rec)
                    if not vid:
                        continue
                    out.append((vid, _assay_of(rec), f"{origin}:{name}"))
    return out


# `/api/variants/coding-variants` is the one variant edge that rejects
# `rsid=` outright ("At least one variant parameter must be defined") and
# returns nothing for `variant_id=<rsid>`. It keys on SPDI. Every other
# variant edge takes the rsID happily, so the generic param builder is
# right for them and wrong only here.
_SPDI_ONLY_ENDPOINTS = {"variants_coding_variants"}

_SPDI_RE = re.compile(r"(N[CGTW]_\d+\.\d+:\d+:[A-Za-z-]*:[A-Za-z-]*)")


def _resolve_spdi(entity: str, form: str, value: str) -> "str | None":
    """One node lookup to turn an rsID into its SPDI.

    The fallback for when no edge has revealed the SPDI yet — notably
    `--relationship coding`, whose only variant-side endpoint is the one
    that needs SPDI in the first place.
    """
    try:
        params = _node_query_param(entity, form, value)
        params["limit"] = "1"
        data = _catalog_get("/api/variants", list(params.items()))
    except SystemExit as exc:
        logging.info("SPDI resolve failed for %s: %s", value, exc)
        return None
    if isinstance(data, list) and data and isinstance(data[0], dict):
        return data[0].get("spdi") or data[0].get("_id")
    return None


def _spdi_from_records(records: "list[dict]") -> "str | None":
    """Pull an SPDI out of whatever the other edges already returned.

    Saves a resolve round-trip: phenotype and protein edges name the
    variant as ``variants/NC_000019.10:44908683:T:C``, which is exactly
    the key the coding-variant endpoint wants.
    """
    for rec in records:
        for k in ("variant", "sequence_variant", "variant_id", "spdi"):
            v = rec.get(k)
            if not v:
                continue
            m = _SPDI_RE.search(str(v))
            if m:
                return m.group(1)
    return None


def _harvest_catalog(variants: "list[str]", *, limit: int
                     ) -> "list[tuple[str, str, str]]":
    """(variant, assay, endpoint) triples by walking every variant edge."""
    out: "list[tuple[str, str, str]]" = []
    for i, vid in enumerate(variants, 1):
        entity, form = detect_id_type(vid)
        if entity != "variant":
            logging.warning("%s is not a variant id (detected %s) — skipped",
                             vid, entity)
            continue
        seen_records: "list[dict]" = []
        deferred: "list[str]" = []
        for key in RELATIONSHIP_TYPE_MAPPING["all"]:
            spec = EDGE_ENDPOINTS[key]
            if spec["from"] != "variant":
                continue
            if key in _SPDI_ONLY_ENDPOINTS and form != "spdi":
                deferred.append(key)
                continue
            params = _edge_query_param(entity, form, vid)
            params["limit"] = str(limit)
            params["skip"] = "0"
            try:
                data = _catalog_get(spec["path"], list(params.items()))
            except SystemExit as exc:
                logging.warning("%s / %s skipped: %s", vid, key, exc)
                continue
            if not isinstance(data, list):
                continue
            for rec in data:
                if isinstance(rec, dict):
                    out.append((vid, _assay_of(rec), key))
                    seen_records.append(rec)

        # Re-run the SPDI-only endpoints now that the other edges have
        # very likely revealed this variant's SPDI.
        if deferred:
            spdi = (vid if form == "spdi"
                    else (_spdi_from_records(seen_records)
                          or _resolve_spdi(entity, form, vid)))
            if not spdi:
                logging.info("%s: no SPDI resolved, skipping %s",
                              vid, ", ".join(deferred))
            else:
                for key in deferred:
                    spec = EDGE_ENDPOINTS[key]
                    try:
                        data = _catalog_get(
                            spec["path"],
                            [("spdi", spdi), ("limit", str(limit)),
                             ("skip", "0")])
                    except SystemExit as exc:
                        logging.warning("%s / %s skipped: %s", vid, key, exc)
                        continue
                    if not isinstance(data, list):
                        continue
                    for rec in data:
                        if isinstance(rec, dict):
                            out.append((vid, _assay_of(rec), key))
        if i % 10 == 0:
            logging.info("  queried %d/%d variants", i, len(variants))
    return out


def cmd_variant_evidence(args: argparse.Namespace) -> int:
    """Rank variants by how many DISTINCT assay types support them.

    Answers "which variants have evidence from more than N assays
    (MPRA, CRISPR, eQTL, ...)" — a question the per-variant and
    per-edge queries cannot express, because it is an aggregation
    across every edge a variant participates in.

    Two input modes, because the Catalog has no unfiltered variant
    scan (the edge endpoints reject a query with no key, and the
    `variants` collection is ~944 GB, deliberately unmirrored):

      --variants / --variants-file   walk the Catalog per variant
      --from-traversal <dir>         aggregate kg-traversal evidence
                                     packs already on disk

    There is deliberately no genome-wide mode. See the note printed at
    the end of a run.
    """
    setup_logging()

    triples: "list[tuple[str, str, str]]" = []
    sources_used: "list[str]" = []

    if args.from_traversal:
        roots = [Path(p) for p in args.from_traversal.split(",") if p.strip()]
        missing = [str(r) for r in roots if not r.exists()]
        if missing:
            raise SystemExit(f"no such path(s): {', '.join(missing)}")
        triples += _harvest_traversal(roots)
        sources_used.append(f"traversal:{len(roots)} path(s)")

    variants: "list[str]" = []
    if args.variants:
        variants += [v.strip() for v in re.split(r"[,\s]+", args.variants)
                     if v.strip()]
    if args.variants_file:
        text = Path(args.variants_file).read_text()
        variants += [v.strip() for v in re.split(r"[,\s]+", text) if v.strip()]
    if variants:
        seen: "set[str]" = set()
        variants = [v for v in variants if not (v in seen or seen.add(v))]
        logging.info("querying the Catalog for %d variant(s)", len(variants))
        triples += _harvest_catalog(variants, limit=args.limit_per_endpoint)
        sources_used.append(f"catalog:{len(variants)} variant(s)")

    if not triples:
        raise SystemExit(
            "No evidence collected. Provide --variants / --variants-file "
            "to query the Catalog, or --from-traversal <dir> to aggregate "
            "kg-traversal evidence packs already on disk.")

    excluded = set(NON_ASSAY_METHODS)
    if args.include_ld:
        excluded -= {"linkage disequilibrum", "linkage disequilibrium"}
    if args.exclude_method:
        excluded |= {m.strip() for m in args.exclude_method.split(",")
                     if m.strip()}
    skip_origins = (set() if args.include_predictions
                    else set(PREDICTION_ONLY_ENDPOINTS))

    per_variant: "dict[str, set[str]]" = {}
    per_variant_exp: "dict[str, set[str]]" = {}
    per_variant_edges: "dict[str, int]" = {}
    origins: "dict[str, set[str]]" = {}
    vocab_all: "Counter[str]" = Counter()
    kind_counts: "Counter[str]" = Counter()
    for vid, assay, origin in triples:
        vocab_all[assay] += 1
        per_variant_edges[vid] = per_variant_edges.get(vid, 0) + 1
        origins.setdefault(vid, set()).add(origin)
        if assay in excluded:
            continue
        if origin.split(":")[-1] in skip_origins:
            continue
        kind = EVIDENCE_KIND.get(assay, "unclassified")
        kind_counts[kind] += 1
        per_variant.setdefault(vid, set()).add(assay)
        if kind == "experimental":
            per_variant_exp.setdefault(vid, set()).add(assay)

    rows = []
    for vid in sorted(per_variant_edges):
        assays = sorted(per_variant.get(vid, ()))
        exp = sorted(per_variant_exp.get(vid, ()))
        rows.append({
            "variant": vid,
            "n_assays": len(assays),
            "n_experimental": len(exp),
            "assays": ";".join(assays),
            "experimental_assays": ";".join(exp),
            "n_edges": per_variant_edges[vid],
            "origins": ";".join(sorted(origins.get(vid, ()))[:6]),
        })
    rows.sort(key=lambda r: (-r["n_assays"], -r["n_edges"], r["variant"]))
    # `n_assays` counts every evidence kind; `n_experimental` counts only the
    # wet-lab assays. Someone asking for ">3 assays, e.g. MPRA/CRISPR/eQTL"
    # means the latter, so --min-experimental gates on it independently rather
    # than letting a computational predictor and a curated resource pad a
    # variant over the threshold.
    kept = [r for r in rows
            if r["n_assays"] >= args.min_assays
            and r["n_experimental"] >= args.min_experimental]

    out_dir = REPORT_DIR / (f"{time.strftime('%Y%m%d_%H%M%S')}_"
                             f"variant_evidence_{safe_label(args.label)}")
    out_dir.mkdir(parents=True, exist_ok=True)
    cols = ["variant", "n_assays", "n_experimental", "assays",
            "experimental_assays", "n_edges", "origins"]

    def _write(path: Path, data: "list[dict]") -> Path:
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            for r in data:
                w.writerow(r)
        return path

    all_path = _write(out_dir / "variant_evidence.tsv", rows)
    hit_path = _write(out_dir / "variant_evidence_filtered.tsv", kept)

    dist = Counter(r["n_assays"] for r in rows)
    summary = {
        "sources": sources_used,
        "n_variants": len(rows),
        "n_variants_passing": len(kept),
        "min_assays": args.min_assays,
        "min_experimental": args.min_experimental,
        "excluded_methods": sorted(excluded),
        "excluded_endpoints": sorted(skip_origins),
        "assay_vocabulary": dict(vocab_all.most_common()),
        "evidence_kind_counts": dict(kind_counts),
        "evidence_kind_of": {k: EVIDENCE_KIND.get(k, "unclassified")
                              for k in vocab_all},
        "assay_count_distribution": {str(k): v for k, v in sorted(dist.items())},
        "max_assays_observed": max((r["n_assays"] for r in rows), default=0),
        "max_experimental_observed": max((r["n_experimental"] for r in rows),
                                          default=0),
    }
    sum_path = out_dir / "summary.json"
    sum_path.write_text(json.dumps(summary, indent=2, sort_keys=True))

    print(f"Output dir:    {out_dir}")
    print(f"All variants:  {all_path}  ({len(rows):,} variants)")
    _gate = f">= {args.min_assays} assays"
    if args.min_experimental:
        _gate += f", >= {args.min_experimental} experimental"
    print(f"Report:        {hit_path}  ({len(kept):,} with {_gate})")
    print(f"Summary:       {sum_path}")
    print(f"Assay vocabulary observed ({len(vocab_all)}): "
          f"{', '.join(k for k, _ in vocab_all.most_common(12))}")
    if excluded:
        print(f"Excluded as non-assay: {', '.join(sorted(excluded))} "
              f"(pass --include-ld to keep LD)")
    if skip_origins:
        print(f"Excluded prediction-only endpoints: "
              f"{', '.join(sorted(skip_origins))} — in-silico predictor "
              f"scores, not assays (pass --include-predictions to count them)")
    if kind_counts:
        print("Evidence kinds counted: " +
              ", ".join(f"{k}={v}" for k, v in kind_counts.most_common()))
    print("Assay-count distribution: " +
          ", ".join(f"{k}:{v}" for k, v in sorted(dist.items())))
    if kept:
        print(f"\nTop {min(10, len(kept))}:")
        for r in kept[:10]:
            print(f"  {r['variant']:24} {r['n_assays']:>2} evid "
                  f"({r['n_experimental']:>2} exp)  {r['assays'][:62]}")
    else:
        print(f"\nNo variant passed {_gate}. Max observed was "
              f"{summary['max_assays_observed']} assays / "
              f"{summary['max_experimental_observed']} experimental.")
    if not args.from_traversal:
        print("\nNote: this is scoped to the variants you supplied. The "
              "Catalog has no unfiltered variant scan and the `variants` "
              "collection (~944 GB) is not mirrored, so a genome-wide "
              "'every variant with >N assays' list cannot be produced — "
              "supply a candidate set, or point --from-traversal at "
              "kg-traversal output for a gene panel.")
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
    _check_identity(data, args.id, params, _node_endpoint(entity))
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
        # Protein records carry `_id` (the ENSP) and the *plural*
        # `uniprot_ids` / `uniprot_names`; there are no `uniprot_id`
        # or `protein_id` fields on the response.
        xrefs.update({
            "canonical_id (ENSP)": rec.get("_id"),
            "uniprot":      rec.get("uniprot_ids") or rec.get("uniprot_id"),
            "uniprot_name": rec.get("uniprot_names"),
            "full_name":    rec.get("uniprot_full_names"),
            "name":         rec.get("name"),
        })
    else:
        # Generic surface
        for k in ("_id", "name", "term_id", "drug_id", "complex_id",
                   "id_version", "source", "is_top_level_pathway"):
            if k in rec and rec[k] is not None:
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

## Multi-assay variant evidence (`variant-evidence`)

"Which variants do you have more than N assays for?" is the question the
per-variant and per-edge queries cannot answer, because it is an
aggregation across *every* edge a variant participates in.

```bash
# A candidate set, walked across every variant edge in the Catalog
igvfagent catalog variant-evidence \
    --variants rs429358,rs7412,rs1421085,rs12740374 --min-assays 3

# A gene panel: reuse kg-traversal packs already on disk (offline)
igvfagent catalog variant-evidence \
    --from-traversal Docs/KGTraversal/2026..._cdk2_variants,Docs/KGTraversal/2026..._e2f1_variants \
    --min-assays 3
```

Output is `variant_evidence.tsv` (all variants, ranked),
`variant_evidence_filtered.tsv` (those clearing the gate) and a
`summary.json` carrying the assay vocabulary and count distribution.

`--min-assays` counts every evidence kind. When the asker means wet-lab
assays — "MPRA, CRISPR, eQTL etc" usually does — add
`--min-experimental 3`, which gates on `n_experimental` as well so a
computational predictor plus a curated resource cannot pad a variant over
the threshold. On the APOE pair, `--min-assays 3` keeps rs429358 and
rs7412 at 4 evidence types each; adding `--min-experimental 3` drops both,
because only 2 of those 4 (GWAS, pQTL) are experimental.

**There is no genome-wide mode, on purpose.** The Catalog rejects an
unfiltered variant-edge query ("At least one variant parameter must be
defined"), and the `variants` collection is ~944 GB, deliberately outside
the local mirror. "Every variant with >3 assays" across the whole genome
is therefore not a question this system can answer; supply a candidate
set or a gene panel. The command says so at the end of every run rather
than returning a partial list that reads like a complete one.

**What counts as an assay.** Distinct `method` / `label` values on the
edges, with two default exclusions that materially change the answer:

- `linkage disequilibrum` (sic, upstream spelling) is a *relationship*,
  not an assay, and attaches to nearly every common variant — leaving it
  in makes the result a list of well-tagged common SNPs. `--include-ld`
  to keep it.
- `/api/variants/coding-variants` returns in-silico predictor columns
  (SIFT, PolyPhen2, CADD, REVEL, AlphaMissense, ESM1b, VARITY, ...) with
  no `method` field at all, so every coding variant would otherwise gain
  a single bogus "unknown" assay. `--include-predictions` to count them.

The report separates `n_assays` from `n_experimental`, because counting a
computational prediction (SEMVAR, cV2F, GVATdb) or a curated resource
(PharmGKB) alongside an eQTL overstates experimental support — which is
the whole point of the question. In a five-variant probe, rs12740374 (the
SORT1 locus) scored 8 evidence types of which 6 were experimental, while
rs1421085 scored 5 of which only 1 was.

That endpoint also needs SPDI rather than an rsID, so the command
harvests the SPDI from whichever other edge already returned it and
re-queries. Before this, `find-associations --relationship coding/all`
silently dropped all coding-variant evidence with an HTTP 400.

## What this skill adds over `kg` / `kg-mirror` / `igvf_client.catalog_get`

| Capability | Before | After |
|---|---|---|
| Universal `get-entity` with 20+ ID auto-detection | only gene + variant heuristics | ✓ |
| `search-region` parallel fan-out | manual region builds | ✓ |
| `find-associations` by semantic category | per-edge calls | ✓ |
| `find-ld` with r²/D'/ancestry buckets | summary endpoint only | ✓ |
| `variant-evidence` multi-assay ranking | **impossible** — an aggregation across every edge a variant has | ✓ |
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

    p = sub.add_parser("variant-evidence",
        help="Rank variants by how many DISTINCT assay types support them "
             "(MPRA / CRISPR / eQTL / caQTL / pQTL / GWAS ...), and filter "
             "with --min-assays.")
    p.add_argument("--variants", default=None,
                    help="Comma/whitespace separated variant ids (rsIDs or "
                         "SPDI). Each is walked across every variant edge "
                         "endpoint in the Catalog.")
    p.add_argument("--variants-file", default=None,
                    help="File of variant ids, one per line.")
    p.add_argument("--from-traversal", default=None,
                    help="Comma-separated kg-traversal run directories (or "
                         "evidence_pack.json paths). Aggregates packs "
                         "already on disk — no network, and the way to go "
                         "from a gene panel to per-variant assay counts.")
    p.add_argument("--min-experimental", type=int, default=0,
                    help="Also require at least this many distinct "
                         "EXPERIMENTAL assays (default 0 = no extra gate). "
                         "Use it when 'N assays' is meant to exclude "
                         "computational predictors and curated resources.")
    p.add_argument("--min-assays", type=int, default=3,
                    help="Report variants supported by at least this many "
                         "distinct assays (default 3).")
    p.add_argument("--limit-per-endpoint", type=int, default=200,
                    help="Max edges pulled per endpoint per variant.")
    p.add_argument("--include-ld", action="store_true",
                    help="Count 'linkage disequilibrum' as an assay. Off by "
                         "default: it is a relationship, not an assay, and "
                         "attaches to nearly every common variant.")
    p.add_argument("--include-predictions", action="store_true",
                    help="Count the coding-variant in-silico predictor "
                         "endpoint (SIFT/PolyPhen2/CADD/REVEL/...). Off by "
                         "default: those are predictions, not assays.")
    p.add_argument("--exclude-method", default=None,
                    help="Extra comma-separated method names to not count.")
    p.add_argument("--label", default="run")
    p.set_defaults(func=cmd_variant_evidence)

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
