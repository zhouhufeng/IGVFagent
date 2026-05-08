#!/usr/bin/env python3
"""IGVF Portal skills for specialized assay and functional data types.

This script searches IGVF Portal metadata for data types that need recurring
agent workflows, writes auditable manifests, and documents download/processing
analysis playbooks. It keeps payload downloads as an explicit follow-up step.
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
REPORT_DIR = DOCS_DIR / "SpecializedIGVF"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
MANIFEST_DIR = DATA_DIR / "Manifests" / "SpecializedIGVF"
RAW_DIR = DATA_DIR / "IGVF" / "SpecializedIGVF"

PORTAL_API_BASE = os.environ.get("IGVF_PORTAL_API_BASE", "https://api.data.igvf.org").rstrip("/")


SKILLS: dict[str, dict[str, Any]] = {
    "parse_split_seq": {
        "title": "Parse SPLiT-seq",
        "purpose": "Access and analyze Parse Biosciences SPLiT-seq single-cell RNA/nuclei datasets.",
        "portal_queries": [
            {"type": "MeasurementSet", "preferred_assay_titles": "Parse SPLiT-seq"},
            {"type": "AnalysisSet", "preferred_assay_titles": "Parse SPLiT-seq"},
            {"type": "File", "preferred_assay_titles": "Parse SPLiT-seq"},
            {"type": "FileSet", "preferred_assay_titles": "Parse SPLiT-seq"},
        ],
        "preferred_files": [
            "gene count matrices",
            "cell annotations",
            "h5ad or sparse matrix archives",
            "barcode/sample demultiplexing metadata",
        ],
        "process_steps": [
            "Build file/sample manifests from MeasurementSet and AnalysisSet records.",
            "Download matrix and cell annotation payloads only after checking file size and access status.",
            "Load RNA counts into AnnData, preserve raw counts, and run mitochondrial/ribosomal and library-size QC.",
            "Normalize, integrate donors/batches, cluster, annotate cell types, and test gene expression across biological groups.",
        ],
    },
    "snatac_scale_preindexing": {
        "title": "10x snATAC-seq With Scale Pre-indexing",
        "purpose": "Retrieve and analyze single-nucleus ATAC-seq libraries that use Scale-style pre-indexing before 10x capture.",
        "portal_queries": [
            {"type": "MeasurementSet", "preferred_assay_titles": "10x snATAC-seq with Scale pre-indexing"},
            {"type": "AnalysisSet", "preferred_assay_titles": "10x snATAC-seq with Scale pre-indexing"},
            {"type": "File", "preferred_assay_titles": "10x snATAC-seq with Scale pre-indexing"},
            {"type": "File", "preferred_assay_titles": "10x snATAC-seq with Scale pre-indexing", "content_type": "fragments"},
        ],
        "preferred_files": [
            "fragments",
            "peak calls or peak matrix",
            "cell annotations",
            "barcode-to-sample mapping and pre-index metadata",
        ],
        "process_steps": [
            "Confirm the pre-index/sample mapping before merging fragments.",
            "Compute ATAC QC metrics: fragments per nucleus, TSS enrichment, fraction in peaks, blacklist fraction, and doublets.",
            "Create or reuse a cell-by-peak matrix; run TF-IDF, LSI, clustering, and cell-type annotation.",
            "Intersect peaks/fragments with variants, cCREs, and candidate regulatory elements.",
        ],
    },
    "multiome_multi_seq": {
        "title": "10x Multiome With MULTI-seq",
        "purpose": "Retrieve paired RNA/ATAC multiome datasets that use MULTI-seq sample multiplexing.",
        "portal_queries": [
            {"type": "AnalysisSet", "preferred_assay_titles": "10x multiome with MULTI-seq"},
            {"type": "MeasurementSet", "preferred_assay_titles": "10x multiome with MULTI-seq"},
            {"type": "File", "preferred_assay_titles": "10x multiome with MULTI-seq"},
            {"type": "File", "searchTerm": "MULTI-seq"},
        ],
        "preferred_files": [
            "RNA sparse gene count matrices",
            "ATAC peak matrices",
            "ATAC fragments",
            "cell annotations with MULTI-seq demultiplexing labels",
        ],
        "process_steps": [
            "Build paired RNA/ATAC manifests and verify matched cell barcodes and sample labels.",
            "Parse MULTI-seq/demultiplexing metadata and remove ambiguous or multiplet assignments.",
            "Process RNA with AnnData/Scanpy and ATAC with Signac, ArchR, or SnapATAC-style tools.",
            "Integrate modalities by barcode and summarize gene expression, accessibility, cell types, and variant-overlap signals.",
        ],
    },
    "share_seq_rna": {
        "title": "SHARE-seq RNA Modality",
        "purpose": "Access the RNA/gene-expression side of SHARE-seq paired chromatin-expression datasets.",
        "portal_queries": [
            {"type": "MeasurementSet", "preferred_assay_titles": "SHARE-seq", "assay_titles": "single-cell RNA sequencing assay"},
            {"type": "AnalysisSet", "preferred_assay_titles": "SHARE-seq"},
            {"type": "File", "preferred_assay_titles": "SHARE-seq", "assay_titles": "single-cell RNA sequencing assay"},
            {"type": "File", "preferred_assay_titles": "SHARE-seq", "content_type": "sparse gene count matrix"},
        ],
        "preferred_files": [
            "gene count matrix",
            "cell annotations",
            "paired ATAC cell identifiers",
            "sample/donor metadata",
        ],
        "process_steps": [
            "Confirm SHARE-seq RNA matrices are paired to matching ATAC/chromatin files.",
            "Load gene counts into AnnData and run standard RNA QC, normalization, clustering, and differential expression.",
            "Keep shared cell identifiers for modality integration.",
            "Use expression results to prioritize genes linked to accessible regulatory elements or variants.",
        ],
    },
    "binding_effect": {
        "title": "Binding Effect",
        "purpose": "Retrieve variant, motif, model, or prediction data that estimate regulatory binding effects.",
        "portal_queries": [
            {"type": "PredictionSet", "file_set_type": "binding effect"},
            {"type": "FileSet", "file_set_type": "binding effect"},
            {"type": "ModelSet", "searchTerm": "binding effect"},
            {"type": "File", "searchTerm": "binding effect"},
        ],
        "preferred_files": [
            "variant-level prediction tables",
            "model output files",
            "motif or TF metadata",
            "reference genome/assembly fields",
        ],
        "process_steps": [
            "Build variant-level manifests with model, assay, biosample, assembly, and score columns.",
            "Normalize score direction and effect allele conventions before joining with variant lists.",
            "Annotate variants with predicted TF binding gain/loss, model provenance, and evidence confidence.",
            "Integrate binding effect evidence with MPRA, CRISPRi, eQTL, accessibility, and enhancer-gene linkage evidence.",
        ],
    },
    "share_seq_atac": {
        "title": "SHARE-seq ATAC/Chromatin Modality",
        "purpose": "Access the chromatin-accessibility side of SHARE-seq paired expression-accessibility datasets.",
        "portal_queries": [
            {"type": "MeasurementSet", "preferred_assay_titles": "SHARE-seq", "assay_titles": "single-cell ATAC-seq"},
            {"type": "AnalysisSet", "preferred_assay_titles": "SHARE-seq"},
            {"type": "File", "preferred_assay_titles": "SHARE-seq", "assay_titles": "single-cell ATAC-seq"},
            {"type": "File", "preferred_assay_titles": "SHARE-seq", "content_type": "fragments"},
        ],
        "preferred_files": [
            "fragments",
            "peak calls",
            "cell-by-peak matrix",
            "cell annotations with paired RNA identifiers",
        ],
        "process_steps": [
            "Verify that ATAC barcodes map to the SHARE-seq RNA modality.",
            "Run fragments/peak QC and construct cell-by-peak matrices when not already provided.",
            "Run TF-IDF/LSI, clustering, motif enrichment, and peak-to-gene or co-accessibility analysis.",
            "Intersect peaks with variants and connect accessibility changes to paired RNA gene expression.",
        ],
    },
    "sge": {
        "title": "SGE / Saturation Genome Editing",
        "purpose": "Retrieve saturation genome editing and other dense variant-function maps from IGVF Portal.",
        "portal_queries": [
            {"type": "MeasurementSet", "preferred_assay_titles": "SGE"},
            {"type": "AnalysisSet", "preferred_assay_titles": "SGE"},
            {"type": "File", "preferred_assay_titles": "SGE"},
            {"type": "FileSet", "preferred_assay_titles": "SGE"},
        ],
        "preferred_files": [
            "variant effect tables",
            "allele/count tables",
            "guide or edit design tables",
            "quality-control summaries",
        ],
        "process_steps": [
            "Build a variant effect manifest with genomic coordinates, alleles, edited sequence, target gene, condition, and score.",
            "Normalize variant identifiers to chr-pos-ref-alt and rsID/CAid where available.",
            "Summarize score distributions, replicate concordance, controls, and significant/deleterious calls.",
            "Join SGE scores to IGVF Catalog variants, genes, enhancer-gene links, MPRA, CRISPRi, and binding-effect evidence.",
        ],
    },
    "sge_variant_annotation": {
        "title": "SGE Variant Annotation And Integration",
        "purpose": "Analyze SGE variant-function tables in the context of IGVF Catalog/KG variant and gene evidence.",
        "portal_queries": [
            {"type": "File", "searchTerm": "variant effect"},
            {"type": "File", "searchTerm": "functional score"},
            {"type": "PredictionSet", "searchTerm": "variant effect"},
            {"type": "AnalysisSet", "preferred_assay_titles": "SGE"},
        ],
        "preferred_files": [
            "scored variant tables",
            "functional annotation tables",
            "gene or target annotations",
            "links to model/prediction provenance",
        ],
        "process_steps": [
            "Parse SGE/variant-effect tables and standardize coordinates, alleles, score names, and score directions.",
            "Flag coding, splice, promoter, enhancer, and UTR contexts with Catalog/KG and ENCODE annotations.",
            "Rank variants by experimental score, predicted binding effect, eQTL/CRISPRi/MPRA evidence, and gene relevance.",
            "Write an annotated variant table plus plots for score distribution, locus tracks, and evidence overlap.",
        ],
    },
}


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"igvf_specialized_data_skills_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)
    return log_path


def safe_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in label)


def build_url(path: str, params: dict[str, Any]) -> str:
    url = f"{PORTAL_API_BASE}{path}"
    return f"{url}?{urllib.parse.urlencode(params, doseq=True, quote_via=urllib.parse.quote)}"


def request_headers() -> dict[str, str]:
    headers = {
        "Accept": "application/json,*/*",
        "User-Agent": "IGVFdataAgent/0.1 specialized-skills",
    }
    if os.environ.get("IGVF_PORTAL_COOKIE"):
        headers["Cookie"] = os.environ["IGVF_PORTAL_COOKIE"]
    return headers


def save_json(label: str, data: Any) -> Path:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    path = RAW_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}.json"
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    logging.info("Saved JSON: %s", path)
    return path


def fetch_search(label: str, params: dict[str, Any]) -> tuple[int, Any, Path, str]:
    query = dict(params)
    query["format"] = "json"
    url = build_url("/search/", query)
    logging.info("Request: GET %s", url)
    request = urllib.request.Request(url, headers=request_headers(), method="GET")
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            content = response.read()
            data = json.loads(content) if "json" in response.headers.get("Content-Type", "") else content.decode(errors="replace")
            saved = save_json(label, data)
            return response.status, data, saved, url
    except urllib.error.HTTPError as exc:
        content = exc.read()
        try:
            data = json.loads(content)
        except json.JSONDecodeError:
            data = {"http_error_body": content.decode(errors="replace")}
        saved = save_json(label, data)
        return exc.code, data, saved, url
    except urllib.error.URLError as exc:
        data = {"network_error": str(exc.reason), "url": url}
        saved = save_json(label, data)
        return 0, data, saved, url


def fetch_item_detail(item_id: str) -> tuple[int, Any]:
    url = build_url(item_id, {"format": "json"})
    logging.info("Request: GET %s", url)
    request = urllib.request.Request(url, headers=request_headers(), method="GET")
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            content = response.read()
            data = json.loads(content) if "json" in response.headers.get("Content-Type", "") else content.decode(errors="replace")
            return response.status, data
    except urllib.error.HTTPError as exc:
        content = exc.read()
        try:
            data = json.loads(content)
        except json.JSONDecodeError:
            data = {"http_error_body": content.decode(errors="replace")}
        return exc.code, data
    except urllib.error.URLError as exc:
        return 0, {"network_error": str(exc.reason), "url": url}


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
        for key in ("accession", "@id", "uuid", "href", "title", "name", "term_name", "summary", "file_format"):
            if key in value:
                return scalar_strings(value[key])
    return []


def text_value(value: Any, limit: int = 5) -> str:
    return "; ".join(scalar_strings(value)[:limit])


def count_field(rows: list[dict[str, Any]], fields: tuple[str, ...], limit: int = 10) -> list[tuple[str, int]]:
    counter: collections.Counter[str] = collections.Counter()
    for row in rows:
        for field in fields:
            for value in scalar_strings(row.get(field)):
                if value:
                    counter[value] += 1
    return counter.most_common(limit)


def download_url(href: str) -> str:
    if not href:
        return ""
    if href.startswith("http://") or href.startswith("https://"):
        return href
    return f"{PORTAL_API_BASE}{href}"


def summarize_query(skill_name: str, query: dict[str, Any], status: int, data: Any, saved: Path, url: str) -> dict[str, Any]:
    rows = rows_from_response(data)
    return {
        "skill": skill_name,
        "query": query,
        "url": url,
        "http_status": status,
        "saved_response": str(saved),
        "returned_rows": len(rows),
        "total": total_from_response(data, rows),
        "assays": count_field(rows, ("assay_title", "assay_titles", "preferred_assay_titles", "assay_term_name")),
        "statuses": count_field(rows, ("status",)),
        "file_formats": count_field(rows, ("file_format", "file_format_type")),
        "content_types": count_field(rows, ("content_type", "content_types", "output_type")),
        "biosamples": count_field(rows, ("sample_summary", "biosample_summary", "sample_terms", "biosample_ontology")),
        "examples": rows[: min(10, len(rows))],
    }


def manifest_rows(summary: dict[str, Any]) -> list[dict[str, str]]:
    rows = []
    for row in summary["examples"]:
        href = text_value(row.get("href"))
        rows.append(
            {
                "skill": summary["skill"],
                "query": json.dumps(summary["query"], sort_keys=True),
                "source_url": summary["url"],
                "accession": text_value(row.get("accession")),
                "id": text_value(row.get("@id") or row.get("uuid")),
                "item_api_url": download_url(text_value(row.get("@id"))),
                "status": text_value(row.get("status")),
                "controlled_access": text_value(row.get("controlled_access")),
                "assay": text_value(row.get("assay_title") or row.get("assay_titles") or row.get("preferred_assay_titles")),
                "sample_summary": text_value(row.get("sample_summary") or row.get("biosample_summary") or row.get("sample_terms")),
                "file_format": text_value(row.get("file_format") or row.get("file_format_type")),
                "content_type": text_value(row.get("content_type") or row.get("content_types") or row.get("output_type")),
                "file_size": text_value(row.get("file_size")),
                "href": href,
                "download_url": download_url(href),
                "summary": text_value(row.get("summary") or row.get("description") or row.get("title")),
            }
        )
    return rows


def format_counter(items: list[tuple[str, int]]) -> str:
    if not items:
        return "- none observed"
    return "\n".join(f"- {name}: {count}" for name, count in items)


def skill_names(selection: str) -> list[str]:
    if selection == "all":
        return list(SKILLS)
    if selection not in SKILLS:
        raise ValueError(f"Unknown skill: {selection}")
    return [selection]


def run_queries(selection: str, limit: int, label: str) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    summaries: list[dict[str, Any]] = []
    manifest: list[dict[str, str]] = []
    for skill_name in skill_names(selection):
        skill = SKILLS[skill_name]
        for index, query in enumerate(skill["portal_queries"]):
            request_query = dict(query)
            request_query["limit"] = str(limit)
            query_label = f"{label}_{skill_name}_{index}_{safe_label('_'.join(f'{k}_{v}' for k, v in request_query.items()))}"
            status, data, saved, url = fetch_search(query_label, request_query)
            summary = summarize_query(skill_name, request_query, status, data, saved, url)
            summaries.append(summary)
            manifest.extend(manifest_rows(summary))
            print(f"{skill_name} query {index}: HTTP {status}, rows={summary['returned_rows']}, total={summary['total']}")
    return summaries, manifest


def write_csv(rows: list[dict[str, str]], label: str, suffix: str) -> Path:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    path = MANIFEST_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_{suffix}.csv"
    fields = sorted({key for row in rows for key in row}) or ["empty"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote CSV: %s", path)
    return path


def write_report(summaries: list[dict[str, Any]], manifest_path: Path, label: str) -> Path:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    path = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_report.md"
    lines = [
        "# IGVF Specialized Data Skills Search Report",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"Manifest: `{manifest_path}`",
        "",
    ]
    for summary in summaries:
        skill = SKILLS[summary["skill"]]
        lines.extend(
            [
                f"## {skill['title']}",
                "",
                f"Query: `{summary['query']}`",
                f"HTTP status: `{summary['http_status']}`",
                f"Returned rows: `{summary['returned_rows']}`",
                f"Reported total: `{summary['total']}`",
                f"Saved response: `{summary['saved_response']}`",
                "",
                "Top assays:",
                format_counter(summary["assays"]),
                "",
                "Top file/content types:",
                format_counter(summary["file_formats"] + summary["content_types"]),
                "",
                "Top samples:",
                format_counter(summary["biosamples"]),
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")
    logging.info("Wrote report: %s", path)
    return path


def write_download_plan(manifest: list[dict[str, str]], label: str) -> Path:
    planned: list[dict[str, str]] = []
    seen: set[str] = set()
    seen_files: set[str] = set()
    for row in manifest:
        query = json.loads(row.get("query") or "{}")
        row_type = query.get("type")
        item_id = row.get("id", "")
        key = row.get("accession") or item_id
        if key in seen:
            continue
        seen.add(key)
        planned_row = dict(row)
        if not planned_row.get("download_url") and row_type == "File" and item_id.startswith("/"):
            status, detail = fetch_item_detail(item_id)
            planned_row["detail_http_status"] = str(status)
            if isinstance(detail, dict):
                href = text_value(detail.get("href"))
                planned_row["href"] = href
                planned_row["download_url"] = download_url(href)
                planned_row["file_size"] = text_value(detail.get("file_size")) or planned_row.get("file_size", "")
                planned_row["file_size_gb"] = (
                    f"{float(detail.get('file_size') or 0) / 1_000_000_000:.3f}"
                    if detail.get("file_size") is not None
                    else ""
                )
                planned_row["md5sum"] = text_value(detail.get("md5sum"))
                planned_row["s3_uri"] = text_value(detail.get("s3_uri"))
                planned_row["content_type"] = text_value(detail.get("content_type")) or planned_row.get("content_type", "")
                planned_row["file_format"] = text_value(detail.get("file_format")) or planned_row.get("file_format", "")
        if planned_row.get("download_url") or row_type == "File":
            planned.append(planned_row)
        if row_type in {"FileSet", "AnalysisSet", "MeasurementSet", "PredictionSet"} and item_id.startswith("/"):
            status, detail = fetch_item_detail(item_id)
            if not isinstance(detail, dict):
                continue
            for file_detail in detail.get("files", []):
                if not isinstance(file_detail, dict):
                    continue
                file_id = text_value(file_detail.get("@id"))
                file_key = text_value(file_detail.get("accession")) or file_id
                if file_key in seen_files:
                    continue
                seen_files.add(file_key)
                href = text_value(file_detail.get("href"))
                file_http_status = status
                if not href:
                    file_status, hydrated = fetch_item_detail(file_id) if file_id.startswith("/") else (0, {})
                    if isinstance(hydrated, dict):
                        file_detail = hydrated
                        href = text_value(file_detail.get("href"))
                    file_http_status = file_status
                if href:
                    planned.append(
                        {
                            "skill": row.get("skill", ""),
                            "query": row.get("query", ""),
                            "source_url": row.get("source_url", ""),
                            "parent_accession": row.get("accession", ""),
                            "parent_id": item_id,
                            "parent_type": row_type,
                            "detail_http_status": str(file_http_status),
                            "accession": text_value(file_detail.get("accession")),
                            "id": text_value(file_detail.get("@id")),
                            "item_api_url": download_url(text_value(file_detail.get("@id"))),
                            "status": text_value(file_detail.get("status")),
                            "controlled_access": text_value(file_detail.get("controlled_access")),
                            "file_format": text_value(file_detail.get("file_format")),
                            "content_type": text_value(file_detail.get("content_type")),
                            "file_size": text_value(file_detail.get("file_size")),
                            "file_size_gb": (
                                f"{float(file_detail.get('file_size') or 0) / 1_000_000_000:.3f}"
                                if file_detail.get("file_size") is not None
                                else ""
                            ),
                            "md5sum": text_value(file_detail.get("md5sum")),
                            "href": href,
                            "download_url": download_url(href),
                            "s3_uri": text_value(file_detail.get("s3_uri")),
                            "summary": text_value(file_detail.get("summary")),
                        }
                    )
    return write_csv(planned, label, "download_plan")


def write_playbook() -> Path:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "IGVF_SPECIALIZED_DATA_SKILLS.md"
    lines = [
        "# Skill: IGVF Specialized Data Access, Download, Processing, And Analysis",
        "",
        "Use these skills when the agent needs to find and work with specialized IGVF Portal datasets for Parse SPLiT-seq, 10x snATAC-seq with Scale pre-indexing, 10x multiome with MULTI-seq, SHARE-seq, binding effect, and SGE workflows.",
        "",
        "## Commands",
        "",
        "```bash",
        "python3 Scripts/igvf_specialized_data_skills.py smoke --skill all --limit 5",
        "python3 Scripts/igvf_specialized_data_skills.py manifest --skill all --limit 25",
        "python3 Scripts/igvf_specialized_data_skills.py download-plan --skill all --limit 25",
        "python3 Scripts/igvf_specialized_data_skills.py write-playbook",
        "```",
        "",
        "## Shared Rules",
        "",
        "- Always run a metadata manifest before downloading payloads.",
        "- Preserve accession, `@id`, assay, sample summary, file format, content type, file size, href, download URL, status, and controlled-access fields.",
        "- Prefer released public files for automated smoke runs. Use `IGVF_PORTAL_COOKIE` only for unreleased data that the local logged-in account is authorized to access.",
        "- Keep large payload downloads explicit and size-capped.",
        "- Store raw API JSON under `Data/IGVF/SpecializedIGVF/`, manifests under `Data/Manifests/SpecializedIGVF/`, reports under `Docs/SpecializedIGVF/`, and logs under `Docs/Logs/`.",
        "",
    ]
    for index, (name, skill) in enumerate(SKILLS.items(), start=1):
        lines.extend(
            [
                f"## {index}. {skill['title']}",
                "",
                f"Skill key: `{name}`",
                "",
                skill["purpose"],
                "",
                "Data access:",
                *[f"- Query IGVF Portal with `{query}`." for query in skill["portal_queries"]],
                "",
                "Download inputs to prefer:",
                *[f"- {item}" for item in skill["preferred_files"]],
                "",
                "Processing and analysis:",
                *[f"{step_index}. {step}" for step_index, step in enumerate(skill["process_steps"], start=1)],
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")
    logging.info("Wrote playbook: %s", path)
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="IGVF specialized data access/download/process/analysis skills.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    for command in ("smoke", "manifest", "download-plan"):
        subparser = subparsers.add_parser(command, help=f"Run {command} searches.")
        subparser.add_argument("--skill", default="all", choices=["all", *SKILLS.keys()])
        subparser.add_argument("--limit", type=int, default=5 if command == "smoke" else 25)
        subparser.add_argument("--label", default=command)

    subparsers.add_parser("write-playbook", help="Write the specialized IGVF data skills document.")

    args = parser.parse_args(argv)
    setup_logging()

    if args.command == "write-playbook":
        path = write_playbook()
        print(f"Wrote {path}")
        return 0

    summaries, manifest = run_queries(args.skill, args.limit, args.label)
    manifest_path = write_csv(manifest, f"{args.label}_{args.skill}", "manifest")
    report_path = write_report(summaries, manifest_path, f"{args.label}_{args.skill}")
    print(f"Wrote manifest: {manifest_path}")
    print(f"Wrote report: {report_path}")
    if args.command == "download-plan":
        download_plan_path = write_download_plan(manifest, f"{args.label}_{args.skill}")
        print(f"Wrote download plan: {download_plan_path}")
    successful = {
        summary["skill"]
        for summary in summaries
        if 200 <= summary["http_status"] < 400 and summary["returned_rows"] > 0
    }
    expected = set(skill_names(args.skill))
    return 0 if expected.intersection(successful) else 1


if __name__ == "__main__":
    raise SystemExit(main())
