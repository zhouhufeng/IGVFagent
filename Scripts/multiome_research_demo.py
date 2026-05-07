#!/usr/bin/env python3
"""Create a 10x multiome research-use-case demo from local IGVF downloads."""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DOC_DIR = ROOT / "Docs" / "Demos" / "10xMultiomeResearchDemo"
PLOT_DIR = DOC_DIR / "Plots"
DATA_DIR = ROOT / "Data"


PAPER_USE_CASES = [
    {
        "theme": "Developmental cell atlas and trajectories",
        "paper": "Single-cell multiomics decodes regulatory programs for mouse secondary palate development",
        "published": "2024-01-27",
        "url": "https://www.nature.com/articles/s41467-024-45199-x",
        "pattern": "Joint RNA and ATAC from the same cells, cell-type annotation, gene activity, marker expression, and developmental-stage trajectories.",
        "demo": "Use IGVF brain-region samples as a mini atlas: summarize cell states by region and modality availability.",
    },
    {
        "theme": "Disease or tissue-state regulatory programs",
        "paper": "Multiomic analysis of human kidney disease identifies a tractable inflammatory and pro-fibrotic tubular cell phenotype",
        "published": "2025-05-22",
        "url": "https://www.nature.com/articles/s41467-025-59997-4",
        "pattern": "Integrate expression and chromatin accessibility to detect disease-associated cell phenotypes and regulatory modules.",
        "demo": "Use brain-region metadata and cell annotations to nominate cell classes and regions for disease/locus follow-up.",
    },
    {
        "theme": "Gene regulatory network inference",
        "paper": "Inferring gene regulatory networks from single-cell multiome data using atlas-scale external data",
        "published": "2024-04-12",
        "url": "https://www.nature.com/articles/s41587-024-02182-7",
        "pattern": "Use paired chromatin accessibility and gene expression to infer TF-gene and enhancer-gene regulatory mechanisms.",
        "demo": "Build a local scaffold: cell-type abundance, modality overlap, and regulatory file inventory for GRN-ready inputs.",
    },
    {
        "theme": "QC and real-cell calling for multiome droplets",
        "paper": "EmptyDropsMultiome discriminates real cells from background in single-cell multiomics assays",
        "published": "2024-05-13",
        "url": "https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03259-x",
        "pattern": "Distinguish true nuclei-containing droplets from background by modeling both RNA and ATAC modalities.",
        "demo": "Summarize cells with RNA-only, ATAC-only, and both-modality flags from IGVF cell annotations.",
    },
    {
        "theme": "Organogenesis regulatory programs",
        "paper": "Single-cell multiomics analysis reveals CTCF as a key regulator of lung morphogenesis and progenitor maintenance",
        "published": "2025-11-28",
        "url": "https://www.nature.com/articles/s41467-025-65757-1",
        "pattern": "Use 10x multiome to connect cell-type specification, accessibility, TF regulation, and gene expression in organ development.",
        "demo": "Mimic the interpretation layer with brain cell classes: top cell states, region enrichment, and data products needed for TF/gene analysis.",
    },
]


def safe_label(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in text)


def find_latest_download_manifest() -> Path:
    manifests = sorted((DATA_DIR / "Manifests" / "Multiome10x").glob("*full_payload_download_manifest.csv"))
    if not manifests:
        raise FileNotFoundError("No full-payload 10x multiome download manifest found.")
    return manifests[-1]


def read_download_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def find_cell_annotation_files(rows: list[dict[str, str]]) -> list[tuple[str, Path]]:
    out = []
    for row in rows:
        path = Path(row.get("path") or "")
        if row.get("status") == "downloaded" and path.exists() and path.name.endswith(".tsv.gz"):
            # Cell annotation files are the tiny TSVs with annotation columns. Peek at the header.
            try:
                with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
                    header = handle.readline()
                if "CL_term_name" in header and "cell_barcode" in header:
                    out.append((row.get("analysis_set_accession", ""), path))
            except OSError:
                continue
    return out


def load_cells(annotation_files: list[tuple[str, Path]]) -> list[dict[str, str]]:
    cells: list[dict[str, str]] = []
    for analysis_set, path in annotation_files:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                row["analysis_set_accession"] = analysis_set
                cells.append(row)
    return cells


def boolish(value: str | None) -> bool:
    return str(value or "").strip().lower() in {"true", "1", "yes"}


def counter_rows(cells: list[dict[str, str]], field: str, limit: int = 20) -> list[tuple[str, int]]:
    counter = collections.Counter((cell.get(field) or "unknown").strip() or "unknown" for cell in cells)
    return counter.most_common(limit)


def region_by_cell_type(cells: list[dict[str, str]], top_cell_terms: list[str]) -> dict[str, dict[str, int]]:
    table: dict[str, dict[str, int]] = collections.defaultdict(lambda: collections.Counter())
    for cell in cells:
        region = (cell.get("BrainRegion") or "unknown").strip() or "unknown"
        term = (cell.get("CL_term_name") or "unknown").strip() or "unknown"
        if term in top_cell_terms:
            table[region][term] += 1
    return {region: dict(counter) for region, counter in table.items()}


def modality_counts(cells: list[dict[str, str]]) -> list[tuple[str, int]]:
    counts = collections.Counter()
    for cell in cells:
        gex = boolish(cell.get("in_GEX"))
        atac = boolish(cell.get("in_ATAC"))
        if gex and atac:
            counts["RNA + ATAC"] += 1
        elif gex:
            counts["RNA only"] += 1
        elif atac:
            counts["ATAC only"] += 1
        else:
            counts["neither flag"] += 1
    return counts.most_common()


def escape_xml(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


def write_bar_plot(items: list[tuple[str, int]], title: str, x_label: str, path: Path) -> None:
    width = 1080
    top, right, bottom, left = 74, 120, 78, 315
    row_h = 31
    height = max(320, top + bottom + row_h * max(1, len(items)))
    max_count = max((value for _, value in items), default=1)
    plot_w = width - left - right
    axis_y = height - bottom
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf7"/>',
        f'<text x="{left}" y="32" font-size="21" font-weight="700" font-family="Arial" fill="#1f2933">{escape_xml(title)}</text>',
        f'<text x="{left}" y="55" font-size="12" font-family="Arial" fill="#5b6470">Horizontal axis: {escape_xml(x_label)}. Vertical axis: category.</text>',
        f'<line x1="{left}" y1="{axis_y}" x2="{width-right}" y2="{axis_y}" stroke="#222"/>',
        f'<line x1="{left}" y1="{top-10}" x2="{left}" y2="{axis_y}" stroke="#222"/>',
        f'<text x="{width/2-58:.0f}" y="{height-24}" font-size="13" font-family="Arial" fill="#1f2933">{escape_xml(x_label)}</text>',
    ]
    for tick in range(6):
        value = max_count * tick / 5
        x = left + (tick / 5) * plot_w
        parts.append(f'<line x1="{x:.2f}" y1="{top-10}" x2="{x:.2f}" y2="{axis_y}" stroke="#e0e4e8"/>')
        parts.append(f'<text x="{x-18:.2f}" y="{axis_y+22}" font-size="11" font-family="Arial" fill="#47515c">{value:.0f}</text>')
    for index, (name, value) in enumerate(items):
        y = top + index * row_h
        bar_w = value / max_count * plot_w if max_count else 0
        parts.append(f'<text x="24" y="{y+20}" font-size="12" font-family="Arial" fill="#28323c">{escape_xml(name[:45])}</text>')
        parts.append(f'<rect x="{left}" y="{y+5}" width="{bar_w:.2f}" height="20" rx="2" fill="#2f776d"/>')
        parts.append(f'<text x="{left + bar_w + 8:.2f}" y="{y+20}" font-size="12" font-family="Arial" fill="#28323c">{value}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_heatmap(table: dict[str, dict[str, int]], cell_terms: list[str], path: Path) -> None:
    regions = sorted(table)
    cell_w, cell_h = 92, 32
    left, top = 215, 92
    width = left + cell_w * len(cell_terms) + 60
    height = top + cell_h * len(regions) + 90
    max_value = max((table[r].get(t, 0) for r in regions for t in cell_terms), default=1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf7"/>',
        '<text x="24" y="32" font-size="21" font-weight="700" font-family="Arial" fill="#1f2933">Cell Type Counts By Brain Region</text>',
        '<text x="24" y="55" font-size="12" font-family="Arial" fill="#5b6470">Darker cells indicate more annotated cells in that region/type combination.</text>',
    ]
    for j, term in enumerate(cell_terms):
        x = left + j * cell_w + 10
        parts.append(f'<text x="{x}" y="82" font-size="10" font-family="Arial" fill="#28323c" transform="rotate(-35 {x} 82)">{escape_xml(term[:20])}</text>')
    for i, region in enumerate(regions):
        y = top + i * cell_h
        parts.append(f'<text x="24" y="{y+21}" font-size="12" font-family="Arial" fill="#28323c">{escape_xml(region)}</text>')
        for j, term in enumerate(cell_terms):
            value = table[region].get(term, 0)
            opacity = 0.12 + 0.88 * (value / max_value if max_value else 0)
            x = left + j * cell_w
            parts.append(f'<rect x="{x}" y="{y}" width="{cell_w-3}" height="{cell_h-3}" fill="#2f776d" opacity="{opacity:.3f}"/>')
            if value:
                parts.append(f'<text x="{x+8}" y="{y+20}" font-size="10" font-family="Arial" fill="#10201d">{value}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_report(summary: dict[str, Any], plots: list[Path], output_json: Path) -> Path:
    report = DOC_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_10x_multiome_research_demo.md"
    lines = [
        "# 10x Multiome Research Use Cases And Local IGVF Demo",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        "",
        "## Recent Research Patterns",
        "",
    ]
    for item in PAPER_USE_CASES:
        lines.extend(
            [
                f"### {item['theme']}",
                "",
                f"Paper: [{item['paper']}]({item['url']})",
                f"Published: {item['published']}",
                "",
                f"Pattern to mimic: {item['pattern']}",
                "",
                f"Local demo: {item['demo']}",
                "",
            ]
        )
    lines.extend(
        [
            "## Local Demo Built From IGVF 10x Multiome Downloads",
            "",
            "This demo mimics the first-pass interpretation layer common to recent 10x Multiome papers: identify the dataset shape, check paired RNA/ATAC coverage, describe cell-state composition, and choose cell types/regions for downstream regulatory analysis.",
            "",
            f"Download manifest: `{summary['download_manifest']}`",
            f"Annotation files parsed: `{summary['annotation_files']}`",
            f"Cells parsed: `{summary['total_cells']}`",
            f"Cells with both RNA and ATAC flags: `{summary['both_modalities']}` ({summary['pct_both']}%)",
            "",
            "Top cell classes:",
            *[f"- {name}: {count}" for name, count in summary["top_cell_terms"]],
            "",
            "Brain regions:",
            *[f"- {name}: {count}" for name, count in summary["brain_regions"]],
            "",
            "## Demo Research Questions Users Can Run",
            "",
            "1. Which brain cell classes are most represented in this IGVF multiome bundle?",
            "2. Which regions have strong paired RNA+ATAC modality coverage?",
            "3. Which cell classes should be prioritized for variant-to-gene or enhancer-to-gene interpretation?",
            "4. Which files are needed for a full GRN analysis: RNA matrices, ATAC peak matrices, fragments, and cell annotations?",
            "5. How would disease-locus interpretation differ if a variant lies in a cell-type-specific accessible element?",
            "",
            "## Plots",
            "",
            *[f"- `{plot}`" for plot in plots],
            "",
            f"Machine-readable summary: `{output_json}`",
            "",
        ]
    )
    report.write_text("\n".join(lines), encoding="utf-8")
    return report


def run_demo(args: argparse.Namespace) -> int:
    DOC_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = Path(args.download_manifest) if args.download_manifest else find_latest_download_manifest()
    rows = read_download_manifest(manifest)
    annotation_files = find_cell_annotation_files(rows)
    cells = load_cells(annotation_files)
    top_terms = [name for name, _ in counter_rows(cells, "CL_term_name", args.top_n)]
    modality = modality_counts(cells)
    regions = counter_rows(cells, "BrainRegion", 20)
    both = dict(modality).get("RNA + ATAC", 0)
    summary = {
        "download_manifest": str(manifest),
        "annotation_files": len(annotation_files),
        "total_cells": len(cells),
        "both_modalities": both,
        "pct_both": round(100 * both / len(cells), 2) if cells else 0,
        "top_cell_terms": counter_rows(cells, "CL_term_name", args.top_n),
        "brain_regions": regions,
        "modality_counts": modality,
    }
    output_json = DOC_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_10x_multiome_research_demo_summary.json"
    output_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    plots = [
        PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_top_cell_types.svg",
        PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_modality_coverage.svg",
        PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_brain_regions.svg",
        PLOT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_region_cell_type_heatmap.svg",
    ]
    write_bar_plot(summary["top_cell_terms"], "Top Cell Classes In IGVF 10x Multiome Demo", "cell count", plots[0])
    write_bar_plot(summary["modality_counts"], "RNA/ATAC Modality Coverage", "cell count", plots[1])
    write_bar_plot(summary["brain_regions"], "Cells By Brain Region", "cell count", plots[2])
    write_heatmap(region_by_cell_type(cells, top_terms[:8]), top_terms[:8], plots[3])
    report = write_report(summary, plots, output_json)
    print(f"Wrote report: {report}")
    print(f"Wrote summary: {output_json}")
    for plot in plots:
        print(f"Wrote plot: {plot}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a 10x multiome research-use-case demo.")
    parser.add_argument("--download-manifest", default="", help="10x multiome full-payload download manifest.")
    parser.add_argument("--top-n", type=int, default=15)
    args = parser.parse_args()
    return run_demo(args)


if __name__ == "__main__":
    raise SystemExit(main())
