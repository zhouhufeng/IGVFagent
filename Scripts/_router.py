"""Deterministic query router for IGVFagent.

Detects the *shape* of a user query (gene symbol, IGVF/ENCODE accession, URL,
variant rsID, genomic region) and maps it to a fixed first tool call — the same
one regardless of which LLM drives the loop. This removes the biggest source of
cross-backend divergence: weaker models picking a different tool (or different
arguments) than Claude for an unambiguous query.

The router is intentionally conservative: it only fires when the query shape is
unmistakable. Everything else returns ``[]`` and the LLM plans freely. The
runtime (``_agent.run``) executes any routed calls first, injects their results,
then lets the LLM continue — so routing *seeds* the plan rather than replacing it.

Pure functions only (no I/O), so the consistency test harness can assert
``route(q)`` is stable without touching a backend.
"""
from __future__ import annotations

import re
from typing import Optional

# rsID: rs followed by >=2 digits (rs429358)
_RE_RSID = re.compile(r"\brs\d{2,}\b", re.IGNORECASE)
# IGVF accession: IGVF + 2 letters + alnum (IGVFDS3222WCZH); ENCODE: ENCSR..., ENCFF...
_RE_IGVF = re.compile(r"\bIGVF[A-Z]{2}[0-9A-Z]{6,}\b")
_RE_ENCODE = re.compile(r"\bENC[A-Z]{2}[0-9A-Z]{6,}\b")
# URL to a known portal
_RE_URL = re.compile(r"https?://\S*(?:igvf|encodeproject|encode)\S*", re.IGNORECASE)
# region: chr1:1,000-2,000 (commas allowed)
_RE_REGION = re.compile(r"\bchr[0-9XYMT]{1,2}:[0-9,]+-[0-9,]+\b", re.IGNORECASE)
# a single bare token that looks like an HGNC gene symbol
_RE_GENE_TOKEN = re.compile(r"^[A-Z][A-Z0-9]{0,9}(?:-[A-Z0-9]{1,4})?$")
# Intent to DO something with the data, not merely learn what it is.
# "raw data" counts on its own: nobody mentions raw reads to be told a title.
_RE_ANALYSE_INTENT = re.compile(
    r"\b(analy[sz]e|analy[sz]ing|analysis|process(?:ing)?|re-?process\w*|"
    r"pipeline|quantif\w+|align\w*|map\s+the\s+reads|count\s+matrix|"
    r"cluster\w*|umap|demultiplex\w*|raw\s+(?:data|reads|fastq))\b",
    re.IGNORECASE)

# All-caps acronyms that are NOT genes — keep bare-token gene routing safe.
_NOT_GENES = {
    "IGVF", "ENCODE", "ATAC", "RNA", "DNA", "SCRNA", "SNRNA", "CRISPR",
    "CRISPRI", "MPRA", "STARR", "GWAS", "QTL", "EQTL", "SQTL", "CCRE", "CRE",
    "TSS", "UMAP", "TSNE", "PCA", "QC", "TF", "TFS", "KG", "LD", "SNP", "SNV",
    "VCF", "BED", "GEO", "SGE", "VEP", "FAVOR", "ABC", "E2G", "WNN", "LSI",
    "HELP", "TODO", "API", "CLI", "UI", "PBMC", "HTO", "GRN", "DEG", "ORA",
    "GSEA", "GO", "PPI", "AQL", "HIC", "PDF", "CSV", "TSV", "JSON",
}


def route(query: str) -> "list[dict]":
    """Return a deterministic list of first tool calls for ``query``.

    Each element is ``{"tool": name, "arguments": {...}, "shape": label}``.
    Returns ``[]`` when no shape matches unambiguously (LLM plans freely).
    Precedence: URL > IGVF/ENCODE accession > rsID > region > bare gene symbol.
    """
    q = (query or "").strip()
    if not q:
        return []

    m = _RE_URL.search(q)
    if m:
        return [{"tool": "explain_dataset",
                 "arguments": {"accession_or_url": m.group(0)},
                 "shape": "url"}]

    m = _RE_IGVF.search(q) or _RE_ENCODE.search(q)
    if m:
        acc = m.group(0)
        # "Analyse / process this data" is a different request from "what is
        # this". Seeding only explain_dataset for both meant the analysis ask
        # was answered with a description plus a list of commands for the
        # user to run themselves. Seed the planner too, so the very first
        # thing on the transcript is the route, the download size, and
        # whether an aligner is present -- and the model continues from a
        # plan instead of from a report that already reads like an answer.
        calls = [{"tool": "explain_dataset",
                  "arguments": {"accession_or_url": acc},
                  "shape": "accession"}]
        if _RE_ANALYSE_INTENT.search(q) and acc.upper().startswith("IGVF"):
            calls.append({"tool": "raw_pipeline_plan",
                          "arguments": {"accession": acc},
                          "shape": "accession_analyse"})
        return calls

    m = _RE_RSID.search(q)
    if m:
        return [{"tool": "kg_variant",
                 "arguments": {"variant": m.group(0)},
                 "shape": "variant_rsid"}]

    m = _RE_REGION.search(q)
    if m:
        return [{"tool": "kg_region",
                 "arguments": {"region": m.group(0).replace(",", "")},
                 "shape": "region"}]

    # bare gene symbol: the WHOLE query is one gene-like token
    tok = q.rstrip("?.!").strip()
    if (" " not in tok and _RE_GENE_TOKEN.match(tok)
            and tok.upper() not in _NOT_GENES
            and any(c.isalpha() for c in tok)
            and len(tok) >= 2):
        return [{"tool": "kg_gene",
                 "arguments": {"symbol": tok.upper()},
                 "shape": "gene_symbol"}]

    return []


def describe() -> str:
    """Human-readable summary of the routing table (for docs / --help)."""
    return (
        "Deterministic routes (same on every backend):\n"
        "  URL to igvf/encode         -> explain_dataset\n"
        "  IGVF*/ENC* accession       -> explain_dataset\n"
        "  ...+ analyse/process intent -> + raw_pipeline_plan\n"
        "  rsID (rs\\d+)               -> kg_variant\n"
        "  chrN:start-end region      -> kg_region\n"
        "  bare gene symbol token     -> kg_gene\n"
        "  anything else              -> (LLM plans freely)"
    )


__all__ = ["route", "describe"]
