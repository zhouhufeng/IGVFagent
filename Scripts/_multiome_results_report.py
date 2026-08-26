#!/usr/bin/env python3
"""Internal helper: build the run-results report for multiome-survey.

Reads the most recent unified manifest under Data/Manifests/MultiomeSurvey/,
aggregates by source / kind, and writes
``Docs/MultiomeSurvey/LATEST_RUN_SUMMARY.md`` plus a per-source ``how to
download`` recipe.

Run after ``multiome_survey.py survey-all --fetch-files``.
"""

from __future__ import annotations

import collections
import csv
import sys
import os
from pathlib import Path

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
MANIFEST_DIR = ROOT / "Data" / "Manifests" / "MultiomeSurvey"
REPORT_DIR = ROOT / "Docs" / "MultiomeSurvey"
OUT_REPORT = REPORT_DIR / "LATEST_RUN_SUMMARY.md"


def fmt_size(b: int) -> str:
    b = float(b or 0)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if b < 1024:
            return f"{b:.2f} {unit}"
        b /= 1024
    return f"{b:.2f} PB"


def latest_unified() -> Path | None:
    cands = sorted(MANIFEST_DIR.glob("*_unified_manifest.csv"))
    return cands[-1] if cands else None


def latest_per_source(source: str) -> Path | None:
    cands = sorted(MANIFEST_DIR.glob(f"*_{source}_files.csv"))
    return cands[-1] if cands else None


def latest_datasets(source: str) -> Path | None:
    cands = sorted(MANIFEST_DIR.glob(f"*_{source}_datasets.csv"))
    return cands[-1] if cands else None


def read_csv(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.exists():
        return []
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main() -> int:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    unified = latest_unified()
    if unified is None:
        print("no unified manifest yet")
        return 2
    rows = read_csv(unified)
    if not rows:
        print("unified manifest is empty")
        return 2

    # Aggregate
    sources_in_order = ["igvf", "encode", "geo", "cellxgene", "hca", "zenodo"]
    by_source: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    for r in rows:
        by_source[r.get("source", "")].append(r)

    src_stats: dict[str, dict[str, object]] = {}
    for src in sources_in_order:
        sr = by_source.get(src, [])
        kind_counts = collections.Counter(r.get("kind", "other") for r in sr)
        kind_sizes: dict[str, int] = collections.defaultdict(int)
        total = 0
        for r in sr:
            sz = int(float(r.get("file_size_bytes") or 0))
            total += sz
            kind_sizes[r.get("kind", "other")] += sz
        n_datasets = len({r.get("dataset_accession") for r in sr})
        src_stats[src] = {
            "n_files": len(sr),
            "n_datasets": n_datasets,
            "total_bytes": total,
            "kinds": dict(kind_counts),
            "kind_sizes": dict(kind_sizes),
        }

    grand_total = sum(int(s["total_bytes"]) for s in src_stats.values())
    grand_files = sum(int(s["n_files"]) for s in src_stats.values())
    grand_datasets = sum(int(s["n_datasets"]) for s in src_stats.values())

    # Top 10 largest files (downloadable)
    rows_sized = [r for r in rows if int(float(r.get("file_size_bytes") or 0)) > 0]
    rows_sized.sort(key=lambda r: -int(float(r["file_size_bytes"])))
    top10 = rows_sized[:10]

    lines: list[str] = []
    lines.append("# Multiome survey — run results and download recipes")
    lines.append("")
    lines.append(f"Source: `{unified.relative_to(ROOT).as_posix()}`  ")
    lines.append(f"Generated: {Path(unified).stat().st_mtime}")
    lines.append("")
    lines.append("## Totals")
    lines.append("")
    lines.append(f"- **{grand_datasets:,} multiome datasets** across the six surveyed sources")
    lines.append(f"- **{grand_files:,} downloadable files indexed**")
    lines.append(f"- **{fmt_size(grand_total)}** total file size if every indexed file were downloaded")
    lines.append("")

    lines.append("## Per-source breakdown")
    lines.append("")
    lines.append("| source | datasets | files | total size | top kinds (count) |")
    lines.append("|---|---:|---:|---:|---|")
    for src in sources_in_order:
        s = src_stats[src]
        kinds = ", ".join(f"{k}={v}" for k, v in sorted(s["kinds"].items(), key=lambda kv: -kv[1])[:4])
        lines.append(
            f"| {src} | {s['n_datasets']:,} | {s['n_files']:,} | {fmt_size(int(s['total_bytes']))} | {kinds} |"
        )
    lines.append("")

    lines.append("## Size by kind, across all sources")
    lines.append("")
    kind_totals: dict[str, int] = collections.defaultdict(int)
    kind_counts: dict[str, int] = collections.defaultdict(int)
    for r in rows:
        kind_totals[r.get("kind", "other")] += int(float(r.get("file_size_bytes") or 0))
        kind_counts[r.get("kind", "other")] += 1
    lines.append("| kind | files | total size |")
    lines.append("|---|---:|---:|")
    for k in sorted(kind_totals, key=lambda x: -kind_totals[x]):
        lines.append(f"| {k} | {kind_counts[k]:,} | {fmt_size(kind_totals[k])} |")
    lines.append("")

    lines.append("## Top-10 largest single files")
    lines.append("")
    lines.append("| source | dataset | file | kind | size |")
    lines.append("|---|---|---|---|---:|")
    for r in top10:
        sz = int(float(r.get("file_size_bytes") or 0))
        lines.append(
            f"| {r.get('source','')} | {r.get('dataset_accession','')[:18]} | "
            f"{(r.get('file_accession') or '')[:36]} | {r.get('kind','')} | {fmt_size(sz)} |"
        )
    lines.append("")

    lines.append("## How to download")
    lines.append("")
    lines.append("All commands resolve their endpoints through ``Scripts/_endpoints.py``; "
                 "no URLs are hard-coded.  Data lands under ``Data/MultiomeSurvey/<source>/<accession>/`` "
                 "which is gitignored.")
    lines.append("")
    lines.append("### Quickstart — RNA + ATAC + annotations only, 20 GB cap")
    lines.append("")
    lines.append("```bash")
    lines.append("igvfagent multiome-survey download \\")
    lines.append("    --only matrix_rna,matrix_atac,fragments,annotations \\")
    lines.append("    --max-download-gb 20")
    lines.append("```")
    lines.append("")
    lines.append("### Just RNA H5AD (CELLxGENE + IGVF + GEO matrices), 10 GB cap")
    lines.append("")
    lines.append("```bash")
    lines.append("igvfagent multiome-survey download \\")
    lines.append("    --only matrix_rna --max-download-gb 10")
    lines.append("```")
    lines.append("")
    lines.append("### Per-source filtered downloads with a regex")
    lines.append("")
    lines.append("```bash")
    lines.append("# Just IGVF h5ad RNA + ATAC fragments")
    lines.append("igvfagent multiome-survey download --pattern 'IGVFFI.*\\.(h5ad|bed\\.gz)$' --max-download-gb 5")
    lines.append("")
    lines.append("# CELLxGENE H5AD only")
    lines.append("igvfagent multiome-survey download --pattern '\\.h5ad$' --max-download-gb 5")
    lines.append("```")
    lines.append("")
    lines.append("### Dry-run first to see what you'd pull")
    lines.append("")
    lines.append("```bash")
    lines.append("igvfagent multiome-survey download \\")
    lines.append("    --only matrix_rna,annotations --max-download-gb 5 --dry-run")
    lines.append("```")
    lines.append("")

    # Per-source recipes
    lines.append("## Per-source download notes")
    lines.append("")
    for src in sources_in_order:
        s = src_stats[src]
        if s["n_files"] == 0:
            continue
        ds_csv = latest_datasets(src)
        fl_csv = latest_per_source(src)
        lines.append(f"### {src.upper()}")
        lines.append("")
        lines.append(f"- Datasets indexed: **{s['n_datasets']:,}**, files indexed: **{s['n_files']:,}**, "
                     f"total: **{fmt_size(int(s['total_bytes']))}**.")
        if ds_csv: lines.append(f"- Datasets manifest: `{ds_csv.relative_to(ROOT).as_posix()}`")
        if fl_csv: lines.append(f"- Files manifest:    `{fl_csv.relative_to(ROOT).as_posix()}`")
        # Top 3 example files
        examples = [r for r in by_source[src] if int(float(r.get("file_size_bytes") or 0)) > 0]
        examples.sort(key=lambda r: -int(float(r["file_size_bytes"])))
        if examples:
            lines.append("- Example largest files:")
            for r in examples[:3]:
                sz = int(float(r.get("file_size_bytes") or 0))
                lines.append(f"    - {r.get('file_accession','')[:48]} ({r.get('kind')}) — {fmt_size(sz)}  → `{r.get('download_url','')[:90]}`")
        lines.append("")

    OUT_REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {OUT_REPORT.relative_to(ROOT).as_posix()}")
    print(f"totals: {grand_datasets} datasets / {grand_files} files / {fmt_size(grand_total)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
