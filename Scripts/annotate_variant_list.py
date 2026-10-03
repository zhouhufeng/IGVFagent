#!/usr/bin/env python3
"""Annotate a variant list with IGVF Catalog / Knowledge Graph evidence.

The script reads CSV variant lists, queries the public IGVF Catalog API, caches
per-variant responses, and writes an augmented CSV plus a JSONL evidence file.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
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
LOG_DIR = ROOT / "Docs" / "Logs"
REPORT_DIR = ROOT / "Docs" / "VariantAnnotation"
OUTPUT_DIR = DATA_DIR / "Annotated" / "VariantList"
CACHE_DIR = DATA_DIR / "Cache" / "IGVFCatalogVariantAnnotations"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint

CATALOG_API_BASE = _resolve_endpoint("catalog_api", "IGVF_CATALOG_API_BASE")

DEFAULT_ENDPOINTS = {
    "summary": ("/api/variants/summary", {}),
    "qtl_genes": ("/api/variants/genes/summary", {"limit": "10"}),
    "phenotypes": ("/api/variants/phenotypes", {"limit": "10", "verbose": "false"}),
    "biosamples": ("/api/variants/biosamples", {"limit": "10", "verbose": "false"}),
    "genomic_elements": ("/api/variants/genomic-elements", {"limit": "10"}),
    "predictions": ("/api/variants/predictions", {"limit": "10"}),
}


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"annotate_variant_list_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)
    return log_path


def safe_label(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in value)


def cache_key(endpoint: str, params: dict[str, Any]) -> str:
    payload = json.dumps({"endpoint": endpoint, "params": params}, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()


def build_url(path: str, params: dict[str, Any]) -> str:
    url = f"{CATALOG_API_BASE}{path}"
    if params:
        url = f"{url}?{urllib.parse.urlencode(params, doseq=True)}"
    return url


def fetch_catalog_json(
    endpoint_name: str,
    path: str,
    params: dict[str, Any],
    *,
    use_cache: bool,
    sleep_seconds: float,
) -> tuple[int, Any, Path]:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    key = cache_key(path, params)
    cache_path = CACHE_DIR / f"{endpoint_name}_{key}.json"
    if use_cache and cache_path.exists():
        try:
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            return int(cached["http_status"]), cached["data"], cache_path
        except (json.JSONDecodeError, KeyError, ValueError):
            logging.warning("Ignoring unreadable cache file: %s", cache_path)

    url = build_url(path, params)
    logging.info("Request: GET %s", url)
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json,*/*",
            "User-Agent": "IGVFdataAgent/0.1",
        },
        method="GET",
    )
    try:
        with urllib.request.urlopen(request, timeout=45) as response:
            content = response.read()
            content_type = response.headers.get("Content-Type", "")
            if "json" in content_type:
                data = json.loads(content)
            else:
                data = content.decode(errors="replace")
            status = response.status
    except urllib.error.HTTPError as exc:
        content = exc.read()
        status = exc.code
        try:
            data = json.loads(content)
        except json.JSONDecodeError:
            data = content.decode(errors="replace")
    except urllib.error.URLError as exc:
        status = 0
        data = {"network_error": str(exc.reason), "url": url}

    cache_path.write_text(
        json.dumps(
            {
                "endpoint": path,
                "params": params,
                "url": url,
                "http_status": status,
                "data": data,
                "cached_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    if sleep_seconds:
        time.sleep(sleep_seconds)
    return status, data, cache_path


def rows_from_response(data: Any) -> list[dict[str, Any]]:
    if isinstance(data, list):
        return [item for item in data if isinstance(item, dict)]
    if isinstance(data, dict):
        for key in ("@graph", "results", "result", "data", "items"):
            value = data.get(key)
            if isinstance(value, list):
                return [item for item in value if isinstance(item, dict)]
    return []


def text_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (str, int, float, bool)):
        return str(value)
    if isinstance(value, list):
        return ";".join(text_value(item) for item in value[:8] if text_value(item))
    if isinstance(value, dict):
        for key in (
            "gene_name",
            "gene_id",
            "phenotype_term",
            "phenotype_id",
            "biological_context",
            "name",
            "method",
            "source",
            "spdi",
            "hgvs",
            "_id",
            "id",
        ):
            if value.get(key):
                return text_value(value[key])
    return json.dumps(value, sort_keys=True)[:200]


def unique_join(values: list[Any], limit: int = 8) -> str:
    seen: list[str] = []
    for value in values:
        text = text_value(value)
        if text and text not in seen:
            seen.append(text)
        if len(seen) >= limit:
            break
    return ";".join(seen)


def infer_variant_params(row: dict[str, str]) -> dict[str, str]:
    spdi = (row.get("SPDI") or "").strip()
    if spdi and spdi.upper() not in {"NA", "N/A", "NONE"}:
        return {"variant_id": spdi}
    rsid = (row.get("rsID") or row.get("rsid") or "").strip()
    if rsid and rsid.lower().startswith("rs"):
        return {"rsid": rsid.lower()}
    ca_id = (row.get("ca_id") or row.get("CA_ID") or "").strip()
    if ca_id:
        return {"ca_id": ca_id}
    chrom = (row.get("hg38Chromosome") or row.get("chrom") or row.get("chr") or "").strip()
    pos = (row.get("hg38Position") or row.get("position") or row.get("pos") or "").strip()
    ref = (row.get("hg38Ref") or row.get("ref") or "").strip()
    alt = (row.get("hg38Alt") or row.get("alt") or "").strip()
    if chrom and pos and ref and alt:
        try:
            zero_based = int(pos) - 1
        except ValueError:
            zero_based = pos
        return {"variant_id": f"chr{chrom}:{zero_based}:{ref}:{alt}"}
    return {}


def summary_fields(data: Any) -> dict[str, str]:
    out = {
        "igvf_kg_spdi": "",
        "igvf_kg_hgvs": "",
        "igvf_kg_rsids": "",
        "igvf_kg_ca_id": "",
        "igvf_kg_cadd_raw": "",
        "igvf_kg_cadd_phred": "",
        "igvf_kg_nearest_gene": "",
        "igvf_kg_nearest_gene_id": "",
        "igvf_kg_nearest_gene_distance": "",
        "igvf_kg_nearest_coding_gene": "",
        "igvf_kg_nearest_coding_gene_id": "",
        "igvf_kg_nearest_coding_gene_distance": "",
    }
    if not isinstance(data, dict):
        return out
    summary = data.get("summary", data)
    if isinstance(summary, dict):
        out["igvf_kg_spdi"] = text_value(summary.get("spdi") or summary.get("variant_id"))
        out["igvf_kg_hgvs"] = text_value(summary.get("hgvs"))
        out["igvf_kg_rsids"] = text_value(summary.get("rsid"))
        out["igvf_kg_ca_id"] = text_value(summary.get("ca_id"))
    cadd = data.get("cadd_scores")
    if isinstance(cadd, dict):
        out["igvf_kg_cadd_raw"] = text_value(cadd.get("raw"))
        out["igvf_kg_cadd_phred"] = text_value(cadd.get("phread"))
    nearest = data.get("nearest_genes")
    if isinstance(nearest, dict):
        nearest_gene = nearest.get("nearestGene")
        nearest_coding_gene = nearest.get("nearestCodingGene")
        if isinstance(nearest_gene, dict):
            out["igvf_kg_nearest_gene"] = text_value(nearest_gene.get("gene_name"))
            out["igvf_kg_nearest_gene_id"] = text_value(nearest_gene.get("id"))
            out["igvf_kg_nearest_gene_distance"] = text_value(nearest_gene.get("distance"))
        if isinstance(nearest_coding_gene, dict):
            out["igvf_kg_nearest_coding_gene"] = text_value(nearest_coding_gene.get("gene_name"))
            out["igvf_kg_nearest_coding_gene_id"] = text_value(nearest_coding_gene.get("id"))
            out["igvf_kg_nearest_coding_gene_distance"] = text_value(nearest_coding_gene.get("distance"))
    return out


def qtl_fields(data: Any) -> dict[str, str]:
    rows = rows_from_response(data)
    return {
        "igvf_kg_qtl_gene_count": str(len(rows)),
        "igvf_kg_qtl_top_genes": unique_join([row.get("gene") for row in rows]),
        "igvf_kg_qtl_contexts": unique_join([row.get("biological_context") for row in rows]),
        "igvf_kg_qtl_types": unique_join([row.get("qtl_type") or row.get("name") for row in rows]),
        "igvf_kg_qtl_max_log10p": text_value(max([row.get("log10pvalue") for row in rows if isinstance(row.get("log10pvalue"), (int, float))], default="")),
    }


def phenotype_fields(data: Any) -> dict[str, str]:
    rows = rows_from_response(data)
    scores = [row.get("score") for row in rows if isinstance(row.get("score"), (int, float))]
    return {
        "igvf_kg_phenotype_count": str(len(rows)),
        "igvf_kg_phenotype_terms": unique_join([row.get("phenotype_term") or row.get("phenotype") for row in rows]),
        "igvf_kg_phenotype_contexts": unique_join([row.get("biological_context") for row in rows]),
        "igvf_kg_phenotype_methods": unique_join([row.get("method") for row in rows]),
        "igvf_kg_phenotype_max_score": text_value(max(scores, default="")),
    }


def evidence_fields(prefix: str, data: Any) -> dict[str, str]:
    rows = rows_from_response(data)
    return {
        f"{prefix}_count": str(len(rows)),
        f"{prefix}_methods": unique_join([row.get("method") for row in rows]),
        f"{prefix}_contexts": unique_join([row.get("biological_context") for row in rows]),
        f"{prefix}_sources": unique_join([row.get("source") for row in rows]),
        f"{prefix}_source_urls": unique_join([row.get("source_url") for row in rows], limit=4),
    }


def prediction_fields(data: Any) -> dict[str, str]:
    rows = rows_from_response(data)
    return {
        "igvf_kg_prediction_count": str(len(rows)),
        "igvf_kg_prediction_target_genes": unique_join([row.get("target_gene") or row.get("gene") for row in rows]),
        "igvf_kg_prediction_contexts": unique_join([row.get("biological_context") for row in rows]),
        "igvf_kg_prediction_sources": unique_join([row.get("source") for row in rows]),
    }


def annotate_row(
    row: dict[str, str],
    *,
    endpoints: dict[str, tuple[str, dict[str, str]]],
    use_cache: bool,
    sleep_seconds: float,
) -> tuple[dict[str, str], dict[str, Any]]:
    params = infer_variant_params(row)
    evidence: dict[str, Any] = {
        "input_variant": row.get("SPDI") or row.get("rsID") or row.get("prioritizedVariantID"),
        "query_params": params,
        "endpoints": {},
    }
    annotation = {
        "igvf_kg_query": json.dumps(params, sort_keys=True),
        "igvf_kg_query_status": "missing_variant_identifier" if not params else "queried",
    }
    if not params:
        return annotation, evidence

    for endpoint_name, (path, extra_params) in endpoints.items():
        request_params = dict(params)
        request_params.update(extra_params)
        status, data, cache_path = fetch_catalog_json(
            endpoint_name,
            path,
            request_params,
            use_cache=use_cache,
            sleep_seconds=sleep_seconds,
        )
        evidence["endpoints"][endpoint_name] = {
            "path": path,
            "params": request_params,
            "http_status": status,
            "cache_path": str(cache_path),
            "data": data,
        }
        annotation[f"igvf_kg_{endpoint_name}_http_status"] = str(status)
        if endpoint_name == "summary":
            annotation.update(summary_fields(data))
        elif endpoint_name == "qtl_genes":
            annotation.update(qtl_fields(data))
        elif endpoint_name == "phenotypes":
            annotation.update(phenotype_fields(data))
        elif endpoint_name == "biosamples":
            annotation.update(evidence_fields("igvf_kg_biosample_evidence", data))
        elif endpoint_name == "genomic_elements":
            annotation.update(evidence_fields("igvf_kg_genomic_element_evidence", data))
        elif endpoint_name == "predictions":
            annotation.update(prediction_fields(data))
    return annotation, evidence


def output_paths(input_path: Path, label: str | None) -> tuple[Path, Path, Path]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    stem = safe_label(label or input_path.stem)
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    return (
        OUTPUT_DIR / f"{timestamp}_{stem}_igvf_annotated.csv",
        OUTPUT_DIR / f"{timestamp}_{stem}_igvf_evidence.jsonl",
        REPORT_DIR / f"{timestamp}_{stem}_annotation_report.md",
    )


def write_report(
    report_path: Path,
    input_path: Path,
    output_csv: Path,
    output_jsonl: Path,
    processed: int,
    missing: int,
    status_counts: dict[str, int],
) -> None:
    lines = [
        "# IGVF Variant Annotation Report",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        f"Input: `{input_path}`",
        f"Annotated CSV: `{output_csv}`",
        f"Evidence JSONL: `{output_jsonl}`",
        "",
        f"Processed rows: {processed}",
        f"Rows missing variant identifiers: {missing}",
        "",
        "## Summary Endpoint Status",
        "",
    ]
    if status_counts:
        lines.extend(f"- HTTP {status}: {count}" for status, count in sorted(status_counts.items()))
    else:
        lines.append("- none")
    lines.extend(
        [
            "",
            "## Annotation Columns",
            "",
            "- `igvf_kg_*_http_status`: HTTP status for each Catalog endpoint.",
            "- `igvf_kg_spdi`, `igvf_kg_hgvs`, `igvf_kg_rsids`, `igvf_kg_ca_id`: variant identity fields.",
            "- `igvf_kg_cadd_*`: CADD scores from the Catalog summary endpoint.",
            "- `igvf_kg_nearest_*`: nearest gene annotations.",
            "- `igvf_kg_qtl_*`: variant-gene / QTL evidence.",
            "- `igvf_kg_phenotype_*`: variant-phenotype evidence.",
            "- `igvf_kg_biosample_evidence_*`: MPRA/STARR/other biosample-linked evidence.",
            "- `igvf_kg_genomic_element_evidence_*`: regulatory element links.",
            "- `igvf_kg_prediction_*`: element-gene prediction links.",
        ]
    )
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Annotate variant CSV with IGVF Catalog/KG evidence.")
    parser.add_argument(
        "--input",
        default=str(DATA_DIR / "Input" / "VariantList" / "example_variants.csv"),
        help="Input CSV path. Provide your own variant list (CSV with rsID/chr-pos-ref-alt columns).",
    )
    parser.add_argument("--label", help="Output label. Defaults to input stem.")
    parser.add_argument("--max-rows", type=int, help="Limit rows for smoke runs.")
    parser.add_argument("--start-row", type=int, default=0, help="Zero-based data row offset.")
    parser.add_argument("--sleep", type=float, default=0.05, help="Seconds to sleep after uncached API requests.")
    parser.add_argument("--no-cache", action="store_true", help="Ignore existing cache files.")
    parser.add_argument(
        "--endpoints",
        default="summary,qtl_genes,phenotypes,biosamples,genomic_elements,predictions",
        help="Comma-separated endpoint set.",
    )
    args = parser.parse_args(argv)

    setup_logging()
    input_path = Path(args.input).resolve()
    selected_names = [name.strip() for name in args.endpoints.split(",") if name.strip()]
    unknown = [name for name in selected_names if name not in DEFAULT_ENDPOINTS]
    if unknown:
        print(f"Unknown endpoint(s): {', '.join(unknown)}")
        return 2
    endpoints = {name: DEFAULT_ENDPOINTS[name] for name in selected_names}
    output_csv, output_jsonl, report_path = output_paths(input_path, args.label)

    processed = 0
    missing = 0
    status_counts: dict[str, int] = {}
    annotation_fields: list[str] = []

    with input_path.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames is None:
            print(f"No CSV header found in {input_path}")
            return 2
        source_fields = list(reader.fieldnames)
        with output_jsonl.open("w", encoding="utf-8") as evidence_handle:
            writer: csv.DictWriter[str] | None = None
            csv_handle = output_csv.open("w", newline="", encoding="utf-8")
            try:
                for row_index, row in enumerate(reader):
                    if row_index < args.start_row:
                        continue
                    if args.max_rows is not None and processed >= args.max_rows:
                        break
                    annotation, evidence = annotate_row(
                        row,
                        endpoints=endpoints,
                        use_cache=not args.no_cache,
                        sleep_seconds=args.sleep,
                    )
                    if annotation.get("igvf_kg_query_status") == "missing_variant_identifier":
                        missing += 1
                    status = annotation.get("igvf_kg_summary_http_status")
                    if status:
                        status_counts[status] = status_counts.get(status, 0) + 1
                    for field in annotation:
                        if field not in annotation_fields:
                            annotation_fields.append(field)
                    if writer is None:
                        writer = csv.DictWriter(
                            csv_handle,
                            fieldnames=source_fields + annotation_fields,
                            extrasaction="ignore",
                        )
                        writer.writeheader()
                    else:
                        missing_fields = [field for field in annotation_fields if field not in writer.fieldnames]
                        if missing_fields:
                            # Keep the initial smoke outputs simple and deterministic. A full run should see
                            # the same annotation fields after the first row because keys are endpoint-driven.
                            logging.warning("New annotation fields after header will be omitted: %s", missing_fields)
                    output_row = dict(row)
                    output_row.update(annotation)
                    writer.writerow(output_row)
                    evidence_handle.write(json.dumps(evidence, sort_keys=True) + "\n")
                    processed += 1
                    if processed % 25 == 0:
                        logging.info("Annotated %s rows", processed)
            finally:
                csv_handle.close()

    write_report(report_path, input_path, output_csv, output_jsonl, processed, missing, status_counts)
    print(f"Annotated rows: {processed}")
    print(f"Wrote CSV: {output_csv}")
    print(f"Wrote evidence JSONL: {output_jsonl}")
    print(f"Wrote report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

