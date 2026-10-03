"""Perturbation Catalogue retrieval skill.

Pulls metadata and per-row perturbation effects from the Perturbation
Catalogue Search API
(``https://perturbation-catalogue-be-…europe-west2.run.app``). The
catalogue currently indexes **~1,222 perturbation datasets** across
three modalities:

  * **MAVE**           — multiplex assays of variant effect (e.g. DMS,
                         VAMP-seq family, deep mutational scanning)
  * **CRISPR screen**  — pooled-cell CRISPR knock-out / KD / activation
                         screens with bulk or FACS readouts
  * **Perturb-seq**    — single-cell CRISPR perturbation with scRNA-seq
                         readout (incl. precomputed GSEA tables)

Subcommands
-----------
    summary                Landing-page statistics (n_datasets, top tissues, …)
    search                 Global gene/target search across all modalities
    search-modality        Modality-scoped search (mave|crispr-screen|perturb-seq)
    dataset                Fetch one dataset record by id
    dataset-rows           Fetch perturbation rows within one dataset
    gsea                   Query the perturb-seq GSEA endpoint
    download               Download bulk modality or dataset data
    pipeline               One-shot: summary + per-modality search +
                           sample downloads for a perturbed gene
    write-playbook         Write Docs/Skills/PERTURBATION_CATALOG_SKILLS.md
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import logging
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Optional


# ---------------------------------------------------------------------------
# Project paths + endpoint
# ---------------------------------------------------------------------------

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data" / "Perturbation"
SEARCHES_DIR = DATA_DIR / "Searches"
DATASETS_DIR = DATA_DIR / "Datasets"
DOWNLOADS_DIR = DATA_DIR / "Downloads"
DOCS_DIR = ROOT / "Docs" / "Perturbation"
LOG_DIR = ROOT / "Docs" / "Logs"
PLAYBOOK_PATH = ROOT / "Docs" / "Skills" / "PERTURBATION_CATALOG_SKILLS.md"

BASE = _resolve_endpoint("perturb_cat", "PERTURB_CAT_BASE")

MODALITIES = ("mave", "crispr-screen", "perturb-seq")

logger = logging.getLogger("perturbation_catalog")


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_") else "_"
                    for c in (s or "run"))


def mkdirs() -> None:
    for d in (DATA_DIR, SEARCHES_DIR, DATASETS_DIR, DOWNLOADS_DIR,
              DOCS_DIR, LOG_DIR):
        d.mkdir(parents=True, exist_ok=True)


def setup_logging(label: str) -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"perturbation_catalog_{timestamp()}_{safe_label(label)}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    return log_path


def _http_get(path: str, *, params: Optional[dict] = None,
               timeout: int = 60) -> "tuple[int, Any]":
    """GET to the catalogue. Returns (status, parsed_json or bytes)."""
    url = BASE + path
    if params:
        clean = {k: v for k, v in params.items()
                  if v is not None and v != ""}
        if clean:
            url += ("&" if "?" in url else "?") + urllib.parse.urlencode(
                clean, doseq=True)
    req = urllib.request.Request(url, headers={
        "Accept": "application/json,*/*",
        "User-Agent": "IGVFagent-PerturbationCatalog/0.1",
    }, method="GET")
    logger.info("GET %s", url)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = resp.read()
            ctype = resp.headers.get("Content-Type", "")
            if "json" in ctype:
                return resp.status, json.loads(body)
            return resp.status, body
    except urllib.error.HTTPError as e:
        body = b""
        try:
            body = e.read()
        except Exception:
            pass
        try:
            return e.code, json.loads(body)
        except Exception:
            return e.code, {"http_error": e.code, "body": body[:400].decode(
                "utf-8", "replace")}
    except urllib.error.URLError as e:
        return 0, {"network_error": str(e.reason)}


def _http_download(path: str, dest: Path, *,
                    params: Optional[dict] = None,
                    timeout: int = 600) -> int:
    """Stream a binary/file response to disk."""
    url = BASE + path
    if params:
        clean = {k: v for k, v in params.items()
                  if v is not None and v != ""}
        if clean:
            url += ("&" if "?" in url else "?") + urllib.parse.urlencode(
                clean, doseq=True)
    req = urllib.request.Request(url, headers={
        "User-Agent": "IGVFagent-PerturbationCatalog/0.1",
    })
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    n = 0
    with urllib.request.urlopen(req, timeout=timeout) as resp, tmp.open("wb") as f:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            n += len(chunk)
    tmp.replace(dest)
    logger.info("Downloaded %s (%.1f MB)", dest.name, n / (1 << 20))
    return n


def _save_json(obj: Any, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str))
    return path


# ---------------------------------------------------------------------------
# Core API wrappers
# ---------------------------------------------------------------------------


def get_summary() -> dict:
    s, d = _http_get("/summary")
    if s != 200 or not isinstance(d, dict):
        raise RuntimeError(f"/summary HTTP {s}: {d}")
    return d


def global_search(query: str, *, page: int = 1, size: int = 25,
                    facets: Optional[str] = None) -> dict:
    params = {"query": query, "page": max(1, page), "size": size}
    if facets:
        params["facets"] = facets
    s, d = _http_get("/search", params=params)
    if s != 200 or not isinstance(d, dict):
        raise RuntimeError(f"/search HTTP {s}: {d}")
    return d


def modality_search(modality: str, *,
                      query: Optional[str] = None,
                      perturbation_gene_name: Optional[str] = None,
                      perturbation_position: Optional[str] = None,
                      effect_score_name: Optional[str] = None,
                      effect_score_value: Optional[str] = None,
                      dataset_limit: int = 25,
                      dataset_offset: int = 0,
                      rows_per_dataset_limit: int = 5,
                      sort: Optional[str] = None,
                      **extra) -> dict:
    """Modality-scoped search: mave | crispr-screen | perturb-seq."""
    if modality not in MODALITIES:
        raise ValueError(f"modality must be one of {MODALITIES}")
    params = {
        "query": query,
        "perturbation_gene_name": perturbation_gene_name,
        "perturbation_position": perturbation_position,
        "effect_score_name": effect_score_name,
        "effect_score_value": effect_score_value,
        "dataset_limit": dataset_limit,
        "dataset_offset": dataset_offset,
        "rows_per_dataset_limit": rows_per_dataset_limit,
        "sort": sort,
    }
    params.update(extra)
    s, d = _http_get(f"/v1/{modality}/search", params=params)
    if s != 200 or not isinstance(d, dict):
        raise RuntimeError(f"/v1/{modality}/search HTTP {s}: {d}")
    return d


def dataset_search(modality: str, dataset_id: str, *,
                     limit: int = 100, offset: int = 0,
                     **extra) -> dict:
    """Search rows *within* one dataset."""
    if modality not in MODALITIES:
        raise ValueError(f"modality must be one of {MODALITIES}")
    params = {"limit": limit, "offset": offset}
    params.update(extra)
    s, d = _http_get(f"/v1/{modality}/{dataset_id}/search", params=params)
    if s != 200 or not isinstance(d, dict):
        raise RuntimeError(f"/v1/{modality}/{dataset_id}/search HTTP {s}: {d}")
    return d


def get_dataset(dataset_id: str) -> dict:
    """Full dataset record."""
    s, d = _http_get(f"/dataset/{dataset_id}")
    if s != 200 or not isinstance(d, dict):
        raise RuntimeError(f"/dataset/{dataset_id} HTTP {s}: {d}")
    return d


def get_gsea(*, query: Optional[str] = None,
              dataset_id: Optional[str] = None,
              page: int = 1, size: int = 50,
              **extra) -> dict:
    params = {"query": query, "dataset_id": dataset_id,
               "page": max(1, page), "size": size}
    params.update(extra)
    s, d = _http_get("/v1/perturb-seq-gsea", params=params)
    if s != 200 or not isinstance(d, dict):
        raise RuntimeError(f"/v1/perturb-seq-gsea HTTP {s}: {d}")
    return d


def download_data(modality: str, *,
                    dataset_id: Optional[str] = None,
                    dest_dir: Optional[Path] = None,
                    **extra) -> Path:
    """Download bulk modality data (optionally scoped to one dataset)."""
    if modality not in MODALITIES:
        raise ValueError(f"modality must be one of {MODALITIES}")
    dest_dir = dest_dir or (DOWNLOADS_DIR / modality)
    dest_dir.mkdir(parents=True, exist_ok=True)
    if dataset_id:
        path = f"/v1/{modality}/{dataset_id}/download"
        fname = f"{safe_label(dataset_id)}.bin"
    else:
        path = f"/v1/{modality}/download"
        fname = f"{modality}_{timestamp()}.bin"
    dest = dest_dir / fname
    n = _http_download(path, dest, params=extra or None)
    # Heuristically rename .bin to a more useful suffix based on the bytes
    new_suffix = _sniff_suffix(dest)
    if new_suffix and new_suffix != ".bin":
        new_dest = dest.with_suffix(new_suffix)
        dest.rename(new_dest)
        return new_dest
    return dest


def download_gsea(*, dataset_id: Optional[str] = None,
                    dest_dir: Optional[Path] = None,
                    **extra) -> Path:
    dest_dir = dest_dir or (DOWNLOADS_DIR / "perturb-seq-gsea")
    dest_dir.mkdir(parents=True, exist_ok=True)
    fname = (f"{safe_label(dataset_id)}.bin" if dataset_id
              else f"gsea_{timestamp()}.bin")
    dest = dest_dir / fname
    params = dict(extra or {})
    if dataset_id:
        params["dataset_id"] = dataset_id
    n = _http_download("/v1/perturb-seq-gsea/download", dest, params=params)
    new_suffix = _sniff_suffix(dest)
    if new_suffix and new_suffix != ".bin":
        new_dest = dest.with_suffix(new_suffix)
        dest.rename(new_dest)
        return new_dest
    return dest


def _sniff_suffix(path: Path) -> Optional[str]:
    """Inspect the first few bytes to guess the right file extension."""
    try:
        with path.open("rb") as f:
            head = f.read(64)
    except Exception:
        return None
    if head.startswith(b"\x1f\x8b"):
        return ".csv.gz"
    if head.startswith(b"PK"):
        return ".zip"
    if head.startswith((b"[", b"{")):
        return ".json"
    if head.startswith(b"\x89PNG"):
        return ".png"
    # Heuristic: CSV vs TSV vs plain TSV from header line
    try:
        line = head.split(b"\n", 1)[0].decode("utf-8", "replace")
    except Exception:
        return None
    if "," in line and "\t" not in line:
        return ".csv"
    if "\t" in line:
        return ".tsv"
    return None


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------


def cmd_summary(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("summary")
    d = get_summary()
    out = DOCS_DIR / f"{timestamp()}_summary"
    out.mkdir(parents=True, exist_ok=True)
    _save_json(d, out / "summary.json")

    md = ["# Perturbation Catalogue summary",
           "",
           f"- Datasets: **{d.get('n_datasets', '?'):,}**  ·  "
           f"experiments: **{d.get('n_experiments', '?'):,}**",
           f"- Year range: **{d.get('min_year', '?')} – {d.get('max_year', '?')}**",
           f"- Targets: **{d.get('n_targets', '?'):,}**  ·  "
           f"tissues: **{d.get('n_tissues', '?')}**  ·  "
           f"cell types: **{d.get('n_cell_types', '?')}**  ·  "
           f"cell lines: **{d.get('n_cell_lines', '?'):,}**  ·  "
           f"diseases: **{d.get('n_diseases', '?')}**",
           ""]
    for fld, title in [
        ("top_modalities", "Top modalities"),
        ("top_tissues", "Top tissues"),
        ("top_cell_types", "Top cell types"),
        ("top_cell_lines", "Top cell lines"),
        ("top_diseases", "Top diseases"),
        ("top_perturbation_types", "Top perturbation types"),
        ("top_sexes", "Top sexes"),
        ("top_dev_stages", "Top developmental stages"),
    ]:
        rows = d.get(fld) or []
        if not rows:
            continue
        md += [f"## {title}", "", "| Value | Datasets |", "|---|---|"]
        for r in rows[:20]:
            md.append(f"| {r.get('value','?')} | {r.get('n_datasets','?'):,} |")
        md.append("")
    (out / "report.md").write_text("\n".join(md))
    print(f"Output: {out}")
    print("\n".join(md[:14]))
    return 0


def cmd_search(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("search_" + safe_label(args.query))
    d = global_search(args.query, page=args.page, size=args.size,
                        facets=args.facets)
    out = SEARCHES_DIR / f"{timestamp()}_global_{safe_label(args.query)}.json"
    _save_json(d, out)
    n_total = d.get("total", 0)
    print(f"Total: {n_total}  ·  saved: {out}")
    for r in (d.get("results") or [])[:args.size]:
        gene = r.get("perturbed_target_symbol", "?")
        line = (f"  {gene:18s}  "
                f"crispr={r.get('n_crispr',0):4d}  "
                f"mave={r.get('n_mave',0):3d}  "
                f"perturb_seq={r.get('n_perturb_seq',0):3d}")
        print(line)
    return 0


def cmd_search_modality(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging(f"search_{args.modality}_"
                              + safe_label(args.query or
                                            args.perturbation_gene_name or "all"))
    d = modality_search(
        args.modality,
        query=args.query,
        perturbation_gene_name=args.perturbation_gene_name,
        perturbation_position=args.perturbation_position,
        effect_score_name=args.effect_score_name,
        effect_score_value=args.effect_score_value,
        dataset_limit=args.dataset_limit,
        dataset_offset=args.dataset_offset,
        rows_per_dataset_limit=args.rows_per_dataset_limit,
        sort=args.sort,
    )
    label = safe_label(args.modality + "_"
                         + (args.query or args.perturbation_gene_name or "all"))
    out = SEARCHES_DIR / f"{timestamp()}_modality_{label}.json"
    _save_json(d, out)
    total = d.get("total_datasets_count") or d.get("total", "?")
    print(f"{args.modality}: total datasets = {total}  ·  saved: {out}")
    for ds in (d.get("datasets") or [])[:10]:
        # Modality search wraps each row as {dataset: {...}, results: [...]}.
        meta = (ds.get("dataset") or ds.get("metadata") or {}) \
            if isinstance(ds, dict) else {}
        ds_id = (meta.get("dataset_id") or ds.get("dataset_id")
                  or ds.get("_id") or "?")
        title = (meta.get("dataset_study_title")
                  or meta.get("dataset_experiment_title")
                  or ds.get("title") or "")[:80]
        print(f"  {ds_id:30s}  {title}")
    return 0


def cmd_dataset(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("dataset_" + safe_label(args.dataset_id))
    d = get_dataset(args.dataset_id)
    out = DATASETS_DIR / f"{safe_label(args.dataset_id)}.json"
    _save_json(d, out)
    print(f"Saved: {out}")
    # Print key metadata fields when present
    for k in ("dataset_id", "dataset_study_title", "dataset_first_author",
               "dataset_study_year", "data_modality", "n_perturbed_targets",
               "dataset_tissues", "dataset_cell_lines", "dataset_diseases"):
        if k in d:
            print(f"  {k:30s}  {str(d[k])[:120]}")
    return 0


def cmd_dataset_rows(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("dataset_rows_"
                              + safe_label(args.modality + "_" + args.dataset_id))
    d = dataset_search(args.modality, args.dataset_id,
                         limit=args.limit, offset=args.offset)
    out = (DATASETS_DIR
            / f"{safe_label(args.dataset_id)}_rows_{args.offset}.json")
    _save_json(d, out)
    print(f"Saved: {out}")
    rows = d.get("results") or d.get("rows") or []
    print(f"Rows: {len(rows)}  (total {d.get('total','?')})")
    return 0


def cmd_gsea(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("gsea_" + safe_label(args.query or
                                                    args.dataset_id or "all"))
    d = get_gsea(query=args.query, dataset_id=args.dataset_id,
                  page=args.page, size=args.size)
    label = safe_label(args.query or args.dataset_id or "all")
    out = SEARCHES_DIR / f"{timestamp()}_gsea_{label}.json"
    _save_json(d, out)
    print(f"Total: {d.get('total','?')}  ·  saved: {out}")
    for r in (d.get("results") or d.get("rows") or [])[:10]:
        term = r.get("term") or r.get("pathway") or r.get("name") or "?"
        nes = r.get("nes") or r.get("score") or "?"
        p = r.get("padj") or r.get("p_value") or "?"
        print(f"  {term[:60]:60s}  NES={nes}  padj={p}")
    return 0


def cmd_download(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("download_"
                              + safe_label(args.modality + "_"
                                            + (args.dataset_id or "all")))
    if args.gsea:
        dest = download_gsea(dataset_id=args.dataset_id)
    else:
        dest = download_data(args.modality, dataset_id=args.dataset_id)
    print(f"Downloaded: {dest}  ({dest.stat().st_size:,} bytes)")
    return 0


def cmd_pipeline(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("pipeline_" + safe_label(args.gene or "run"))
    out = DOCS_DIR / f"{timestamp()}_{safe_label(args.label or args.gene)}"
    out.mkdir(parents=True, exist_ok=True)
    summary = get_summary()
    _save_json(summary, out / "summary.json")
    gene_hit = global_search(args.gene, size=5)
    _save_json(gene_hit, out / "global_search.json")
    modality_hits: "dict[str, Any]" = {}
    for m in MODALITIES:
        try:
            modality_hits[m] = modality_search(
                m, perturbation_gene_name=args.gene,
                dataset_limit=args.dataset_limit,
                rows_per_dataset_limit=args.rows_per_dataset_limit,
            )
            _save_json(modality_hits[m], out / f"{m}_search.json")
        except Exception as e:
            logger.warning("modality %s: %s", m, e)
            modality_hits[m] = {"error": str(e)}

    # Per-modality summary
    md = ["# Perturbation Catalogue — pipeline run for "
          f"`{args.gene}`", "",
          f"- Catalogue total datasets: {summary.get('n_datasets', '?'):,}",
          ""]
    if gene_hit.get("results"):
        r = gene_hit["results"][0]
        md += [f"## Global gene record — {r.get('perturbed_target_symbol','?')}",
                "",
                f"- CRISPR-screen datasets: **{r.get('n_crispr', 0)}** "
                f"(significant: {r.get('n_sig_crispr', 0)})",
                f"- MAVE datasets: **{r.get('n_mave', 0)}**",
                f"- Perturb-seq datasets: **{r.get('n_perturb_seq', 0)}** "
                f"(up: {r.get('n_sig_perturb_pairs_up', 0)}, "
                f"down: {r.get('n_sig_perturb_pairs_down', 0)})",
                f"- Top GSEA terms: {', '.join((r.get('top_gsea_terms') or [])[:5])}",
                ""]

    for m, h in modality_hits.items():
        if isinstance(h, dict) and "error" in h:
            md += [f"## {m}", "", f"_(query failed: {h['error']})_", ""]
            continue
        n = h.get("total_datasets_count") or h.get("total") or 0
        md += [f"## {m}  ({n} datasets)", "",
                "| Dataset id | Study title |", "|---|---|"]
        for ds in (h.get("datasets") or [])[:10]:
            # Modality search wraps each row as {dataset: {...}, results: [...]}.
            meta = ((ds.get("dataset") or ds.get("metadata") or {})
                    if isinstance(ds, dict) else {})
            ds_id = (meta.get("dataset_id") or ds.get("dataset_id")
                      or ds.get("_id") or "?")
            title = (meta.get("dataset_study_title")
                      or meta.get("dataset_experiment_title")
                      or "")[:100]
            md.append(f"| `{ds_id}` | {title} |")
        md.append("")

    (out / "report.md").write_text("\n".join(md))
    print(f"Output: {out}")
    return 0


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------


PLAYBOOK_TEXT = """\
# Skill: Perturbation Catalogue retrieval

Pulls metadata and per-row perturbation effects from the Perturbation
Catalogue Search API. The catalogue indexes **~1,222 perturbation
datasets** across MAVE (DMS / VAMP-seq family), CRISPR screens, and
Perturb-seq.

## Subcommands

### summary
```
igvfagent perturb-catalog summary
```
Landing-page stats: total datasets, top modalities / tissues / cell
types / diseases / perturbation types. Writes
`Docs/Perturbation/<ts>_summary/{summary.json, report.md}`.

### search (global, gene-level)
```
igvfagent perturb-catalog search --query BRCA1 --size 10
```
Returns one row per perturbed gene with counts in each modality, top
GSEA terms, and faceted tissue/cell-type info. Use this when the
user asks "what does the catalogue have on gene X?".

### search-modality
```
igvfagent perturb-catalog search-modality --modality crispr-screen \\
    --perturbation-gene-name BRCA1 --dataset-limit 25
igvfagent perturb-catalog search-modality --modality mave \\
    --perturbation-gene-name TP53 \\
    --perturbation-position 100_300 \\
    --effect-score-name vamp_score
igvfagent perturb-catalog search-modality --modality perturb-seq \\
    --query "lung cancer" --dataset-limit 10
```
Modality-scoped search returns datasets + sample rows per dataset.
Supports filters on perturbation gene name, position range (e.g.
`100_300`), score name and value range (e.g. `0.5_1.0`), tissue,
cell line, disease, study year, library generation type, sequencing
platform, and many more (full filter list mirrors the catalogue's
faceted UI).

### dataset
```
igvfagent perturb-catalog dataset --dataset-id <id>
```
Full dataset record (single JSON document).

### dataset-rows
```
igvfagent perturb-catalog dataset-rows --modality mave \\
    --dataset-id <id> --limit 500 --offset 0
```
Paginate through the **per-perturbation rows** inside one dataset
(variant- or gRNA-level effect scores).

### gsea
```
igvfagent perturb-catalog gsea --query BRCA1 --size 50
igvfagent perturb-catalog gsea --dataset-id <id>
```
Query the Perturb-seq GSEA endpoint for hallmark/pathway enrichment
tables.

### download
```
igvfagent perturb-catalog download --modality crispr-screen --dataset-id <id>
igvfagent perturb-catalog download --gsea --dataset-id <id>
```
Bulk download of one dataset (or all rows for a modality). Output
extension is auto-sniffed (`.csv.gz`, `.json`, `.zip`, …) and saved
under `Data/Perturbation/Downloads/<modality>/`.

### pipeline (one-shot)
```
igvfagent perturb-catalog pipeline --gene BRCA1
```
Runs summary + global gene search + modality-scoped searches for
each of MAVE / CRISPR-screen / Perturb-seq and writes a markdown
report under `Docs/Perturbation/<ts>_<gene>/`. The headline command
to call when the user asks "what perturbation data exists for
gene X?".

### write-playbook
```
igvfagent perturb-catalog write-playbook
```

## Cross-skill chaining

- `proteomics vampseq-analyze` (canonical MAVE / VAMP-seq scoresets
  from MaveDB) → `perturb-catalog search-modality --modality mave`
  for related Perturbation Catalogue datasets.
- `rnaseq deg` → `perturb-catalog search` for each significant gene
  to pull CRISPR / Perturb-seq evidence.
- `kg gene <SYM>` (IGVF Catalog) → `perturb-catalog pipeline --gene <SYM>`
  to attach perturbation evidence to a gene-centric workflow.
- `proteomics kg-visualize --gene <SYM>` → uses top hubs from the
  Perturb-seq GSEA terms to suggest functional cofactors.

## Output layout

```
Data/Perturbation/
    Searches/<ts>_<scope>.json
    Datasets/<dataset_id>.json
    Downloads/<modality>/<dataset_id>.<auto-sniffed-ext>
    Downloads/perturb-seq-gsea/<id>.csv.gz
Docs/Perturbation/<ts>_<label>/
    summary.json
    global_search.json
    {mave,crispr-screen,perturb-seq}_search.json
    report.md
```
"""


def cmd_write_playbook(_a) -> int:
    PLAYBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    PLAYBOOK_PATH.write_text(PLAYBOOK_TEXT)
    print(f"Wrote: {PLAYBOOK_PATH}")
    return 0


# ---------------------------------------------------------------------------
# argparse
# ---------------------------------------------------------------------------


def main(argv: "Optional[list[str]]" = None) -> int:
    p = argparse.ArgumentParser(
        prog="perturb-catalog",
        description="Perturbation Catalogue retrieval: MAVE / CRISPR-screen "
                    "/ Perturb-seq datasets, perturbation rows, GSEA tables, "
                    "and bulk downloads.")
    sub = p.add_subparsers(dest="cmd")

    s = sub.add_parser("summary",
                        help="Landing-page stats (n_datasets, top facets).")
    s.set_defaults(func=cmd_summary)

    s = sub.add_parser("search",
                        help="Global gene/term search across all modalities.")
    s.add_argument("--query", required=True)
    s.add_argument("--page", type=int, default=1)
    s.add_argument("--size", type=int, default=25)
    s.add_argument("--facets", default=None,
                    help="Comma-separated facet names to include.")
    s.set_defaults(func=cmd_search)

    s = sub.add_parser("search-modality",
                        help="Modality-scoped search "
                              "(mave|crispr-screen|perturb-seq).")
    s.add_argument("--modality", required=True, choices=MODALITIES)
    s.add_argument("--query", default=None)
    s.add_argument("--perturbation-gene-name", default=None)
    s.add_argument("--perturbation-position", default=None,
                    help="Position or range (e.g. 100 or 100_300).")
    s.add_argument("--effect-score-name", default=None)
    s.add_argument("--effect-score-value", default=None,
                    help="Score range (e.g. 0.5_1.0).")
    s.add_argument("--dataset-limit", type=int, default=25)
    s.add_argument("--dataset-offset", type=int, default=0)
    s.add_argument("--rows-per-dataset-limit", type=int, default=5)
    s.add_argument("--sort", default=None)
    s.set_defaults(func=cmd_search_modality)

    s = sub.add_parser("dataset", help="Fetch one dataset record by id.")
    s.add_argument("--dataset-id", required=True)
    s.set_defaults(func=cmd_dataset)

    s = sub.add_parser("dataset-rows",
                        help="Fetch perturbation rows within one dataset.")
    s.add_argument("--modality", required=True, choices=MODALITIES)
    s.add_argument("--dataset-id", required=True)
    s.add_argument("--limit", type=int, default=100)
    s.add_argument("--offset", type=int, default=0)
    s.set_defaults(func=cmd_dataset_rows)

    s = sub.add_parser("gsea",
                        help="Query the Perturb-seq GSEA endpoint.")
    s.add_argument("--query", default=None)
    s.add_argument("--dataset-id", default=None)
    s.add_argument("--page", type=int, default=1)
    s.add_argument("--size", type=int, default=50)
    s.set_defaults(func=cmd_gsea)

    s = sub.add_parser("download",
                        help="Download bulk modality / dataset data.")
    s.add_argument("--modality", default=MODALITIES[1], choices=MODALITIES)
    s.add_argument("--dataset-id", default=None)
    s.add_argument("--gsea", action="store_true",
                    help="Use the perturb-seq-gsea/download endpoint.")
    s.set_defaults(func=cmd_download)

    s = sub.add_parser("pipeline",
                        help="One-shot: summary + per-modality search for "
                              "a perturbed gene.")
    s.add_argument("--gene", required=True)
    s.add_argument("--label", default=None)
    s.add_argument("--dataset-limit", type=int, default=10)
    s.add_argument("--rows-per-dataset-limit", type=int, default=3)
    s.set_defaults(func=cmd_pipeline)

    s = sub.add_parser("write-playbook",
                        help="Write Docs/Skills/PERTURBATION_CATALOG_SKILLS.md")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args(argv)
    if not getattr(args, "func", None):
        p.print_help()
        return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
