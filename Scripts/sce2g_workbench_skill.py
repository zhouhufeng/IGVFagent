#!/usr/bin/env python3
"""scE2G workbench: set up, configure, check, run and benchmark scE2G model training.

Absorbs the Engreitz-lab scE2G training workflow (https://github.com/EngreitzLab/scE2G,
MIT) and the CRISPR benchmarking pipeline it is evaluated with
(https://github.com/EngreitzLab/CRISPR_comparison, MIT), following the group's
"Quick Start on Building New E2G Models" walkthrough for adding crowdsourced
features to the multiome model. The feature table that the main branch lacks
comes from https://github.com/kaybrand/scE2G (rebase/multiple-cell-types).

Relationship to the upstreams: scE2G itself is NOT reimplemented -- it is a
Snakemake + R + Python workflow with a Singularity container and a Slurm
profile, and the model's scientific asset is its trained weights. This module
WRAPS it: it prepares a checkout the way the walkthrough says, writes every
config file in the exact column order upstream reads, checks that all the
paths line up before hours of compute are spent, and runs (or prints) the
snakemake command. The CRISPR benchmark's evaluation is a clean-room
reimplementation of the CRISPR_comparison overlap-and-aggregate scoring, PR
curves, AUPRC, precision at 70% recall and bootstrap intervals, and it also
writes a `pred_config.txt` + `config.yml` so the upstream pipeline can be run
on the same inputs for a byte-level comparison.

What the walkthrough says, and what each subcommand does about it:

  setup      clone scE2G (default branch fix/dag-staleness-integration, which
             stops Snakemake re-running rules needlessly), init the ENCODE-rE2G
             submodule, then apply the four patches: drop `conda: "mamba"` from
             both Snakefile_training files, add SCRIPTS_DIR to the top-level
             one, add RNA_matrix_filtered / max_cell_count to
             config_training.yaml, and put the missing
             resources/feature_tables/multiome_arc_n6.tsv in place. Idempotent;
             `--check` only reports.
  features   turn a crowdsourced E2G feature table (ElementChr/ElementStart/
             ElementEnd/GeneSymbol + feature columns, as shared on Synapse) into
             the `source_file` scE2G merges (chr/start/end/TargetGene, spaces in
             feature names -> underscores, .tsv.gz), plus the five-column
             external_features_config_<name>.tsv (aggregate_function=mean,
             join_by=overlap) and a feature_table_<name>.tsv extending
             multiome_arc_n6 with one row per feature (aggregate_function=max,
             fill_value=0, nice_name with the spaces put back).
  configure  add or replace the cluster row in config_cell_clusters.tsv
             (rna_matrix_file, atac_frag_file, model_dir, and the extra
             external_features_config column) and the model row in
             config_models.tsv (dataset == cluster, ABC_directory blank,
             polynomial False), and write run_training_<model>.sh with the
             Slurm-profile snakemake command and the memory/runtime overrides.
  check      every path referenced by the configs exists; cluster == dataset;
             the feature_table's extra rows match the external config's
             input_cols; the source_file has chr/start/end/TargetGene and the
             source_cols; no feature name contains a space; reference files in
             config_training.yaml exist. Exit 1 with the list of problems.
  run        `check`, then execute snakemake (dry-run by default) with the
             profile / job count you ask for. Needs snakemake on PATH.
  predictions  describe one or more scE2G prediction tables (*.e2g.tsv): links,
             genes, elements, score quantiles, element classes, links above the
             threshold, and pairwise overlap between tables.
  benchmark  CRISPR_comparison-style evaluation of prediction tables against an
             EPCrisprBenchmark CRISPR file: per predictor aggregate_function,
             fill_value, inverse_predictor and boolean semantics; AUPRC,
             precision at 70% recall, bootstrap 95% intervals, PR-curve figure;
             writes pred_config.txt and config.yml for the upstream pipeline.
  selftest   builds a fake scE2G checkout and synthetic feature / prediction /
             CRISPR files in a temp dir, runs setup --check + patches, features,
             configure, check and benchmark, and asserts on every artefact.

Validated against the upstream pipeline itself (2026-09-22, CRISPR_comparison
@ 5058742 run in a conda container on the same inputs): the per-pair merged
scores were identical for all 10,356 K562 CRISPR pairs, precision at 70%
recall matched to three decimals, and `auprc_crispr_comparison` -- the
trapezoid over yardstick's tie-aware PR points with the first and last rows
dropped, which is what upstream's performance_summary.txt reports -- matched
to four decimals (rE2G base 0.6326, scE2G 0.5190, Pinloop 0.3311). Two
earlier defects were found by that comparison and fixed: inverse predictors
were negated after fill/aggregation (distance baseline read 0.106 instead of
~0.42), and tied scores were ranked one at a time (Pinloop moved by 0.015
with the tie order). `auprc` (step rule, scikit-learn convention) is kept for
continuity with benchmark #12.

Standard library for setup/configure/check/run; pandas for features,
predictions and benchmark; matplotlib optional (figures skipped without it).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT")
            or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "scE2G"

SCE2G_REPO = "https://github.com/EngreitzLab/scE2G.git"
SCE2G_DEFAULT_BRANCH = "fix/dag-staleness-integration"
FEATURE_TABLE_URL = ("https://raw.githubusercontent.com/kaybrand/scE2G/"
                     "rebase/multiple-cell-types/resources/feature_tables/multiome_arc_n6.tsv")
# The file as published in that fork on 2026-09-22, so `setup` works offline
# and a later upstream change is visible as a diff rather than a silent swap.
MULTIOME_ARC_N6 = """feature\tinput_col\tsecond_input\taggregate_function\tfill_value\tnice_name
numTSSEnhGene\tnumTSSEnhGene\tNA\tmax\t0\t# TSSs between E and P
normalizedATAC_prom\tnormalized_atac_prom\tNA\tmean\t0\tATAC signal at P
numNearbyEnhancers\tnumNearbyEnhancers\tNA\tmax\t0\t# peaks within 5Kb of E
ubiqExpressed\tis_ubiquitous_uniform\tNA\tmax\t0\tUbiquitous expression
numCandidateEnhGene\tnumCandidateEnhGene\tNA\tmax\t0\t# peaks between E and P
ARC.E2G.Score\tARC.E2G.Score\tNA\tmean\t0\tARC-E2G score
"""

CLUSTER_COLS = ["cluster", "rna_matrix_file", "atac_frag_file", "HiC_file", "HiC_type",
                "HiC_resolution", "alt_TSS", "alt_genes", "model_dir"]
EXT_COL = "external_features_config"
MODEL_COLS = ["model", "dataset", "ABC_directory", "feature_table", "polynomial", "override_params"]
EXT_CONFIG_COLS = ["input_col", "source_col", "aggregate_function", "join_by", "source_file"]
FEATURE_TABLE_COLS = ["feature", "input_col", "second_input", "aggregate_function", "fill_value", "nice_name"]
SOURCE_RENAME = {"ElementChr": "chr", "ElementStart": "start", "ElementEnd": "end", "GeneSymbol": "TargetGene",
                 "chrom": "chr", "TargetGene": "TargetGene"}
SOURCE_REQUIRED = ["chr", "start", "end", "TargetGene"]
PRED_CONFIG_COLS = ["pred_id", "pred_col", "boolean", "alpha", "aggregate_function", "fill_value",
                    "inverse_predictor", "pred_name_long", "color"]
TRAINING_YAML_ADDITION = """
# --- added by igvfagent sce2g setup (crowdsourced-features walkthrough) ---
# set this to True if the RNA matrix contains the exact set of cells as the ATAC fragment file, and False if it contains more cells (default if not specified: True)
RNA_matrix_filtered: True
# If the cell count in a cell type exceeds the max_cell_count, randomly extract max_cell_count cells.
max_cell_count: 20000
"""
SNAKEMAKE_TRAINING = ("snakemake -s workflow/Snakefile_training --profile {profile} "
                      "--configfile config/config_training.yaml --rerun-incomplete "
                      "--set-resources make_external_features_config:mem_mb=64000 "
                      "make_external_features_config:runtime=1440 add_external_features:mem_mb=256000")
SNAKEMAKE_TRAINING_LOCAL = ("snakemake -s workflow/Snakefile_training -j {jobs} --use-conda "
                            "--configfile config/config_training.yaml --rerun-incomplete")

BLUE, ORANGE, INK, INK2, AXIS, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"sce2g_workbench_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def run_dir(label: str) -> Path:
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _open_text(path: Path, mode: str = "rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_tsv_rows(path: Path) -> "tuple[list[str], list[dict]]":
    with _open_text(path) as fh:
        rows = [ln for ln in fh if not ln.startswith("#")]
    if not rows:
        return [], []
    reader = csv.DictReader(rows, delimiter="\t")
    header = list(reader.fieldnames or [])
    return header, [dict(r) for r in reader]


def write_tsv_rows(path: Path, header: "list[str]", rows: "Iterable[dict]") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", extrasaction="ignore",
                           lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in header})
    return path


def _pd():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit("this subcommand needs pandas: pip install 'igvfagent[analysis]'") from e


def _plt():
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
        return plt
    except Exception:
        return None


def _style(ax, title: str, xlabel: str = "", ylabel: str = "") -> None:
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold", color=INK)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.set_facecolor(SURFACE)


def md_table(headers: "list[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 60) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for i, r in enumerate(rows):
        if i >= max_rows:
            lines.append("| … |" + " |" * (len(headers) - 1))
            break
        lines.append("| " + " | ".join(str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(lines)


def relpath(p: "str | Path", base: Path) -> str:
    try:
        return str(Path(p).resolve().relative_to(base.resolve()))
    except ValueError:
        return str(p)


# ---------------------------------------------------------------------------
# setup: checkout + the walkthrough's patches
# ---------------------------------------------------------------------------

_MAMBA_RE = re.compile(r'^\s*conda:\s*["\']mamba["\']\s*$')
SCRIPTS_DIR_LINE = 'SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "workflow", "scripts")'


def patch_snakefile(path: Path, add_scripts_dir: bool, apply: bool) -> "list[str]":
    """Drop `conda: "mamba"`; add SCRIPTS_DIR after WORKFLOW_DIR. Returns notes."""
    notes = []
    if not path.is_file():
        return [f"missing: {path}"]
    lines = path.read_text().splitlines()
    kept = [ln for ln in lines if not _MAMBA_RE.match(ln)]
    n_dropped = len(lines) - len(kept)
    if n_dropped:
        notes.append(f"{path.name}: {'dropped' if apply else 'would drop'} {n_dropped} `conda: \"mamba\"` line(s)")
    else:
        notes.append(f"{path.name}: no `conda: \"mamba\"` line (already patched)")
    if add_scripts_dir:
        if any(ln.strip().startswith("SCRIPTS_DIR =") for ln in kept):
            notes.append(f"{path.name}: SCRIPTS_DIR already defined")
        else:
            idx = next((i for i, ln in enumerate(kept) if ln.strip().startswith("WORKFLOW_DIR =")), None)
            if idx is None:
                notes.append(f"{path.name}: WARNING no `WORKFLOW_DIR =` line; SCRIPTS_DIR not added")
            else:
                kept.insert(idx + 1, SCRIPTS_DIR_LINE)
                notes.append(f"{path.name}: {'added' if apply else 'would add'} SCRIPTS_DIR after WORKFLOW_DIR")
    if apply and kept != lines:
        path.write_text("\n".join(kept) + "\n")
    return notes


def patch_training_yaml(path: Path, apply: bool) -> "list[str]":
    if not path.is_file():
        return [f"missing: {path}"]
    text = path.read_text()
    have = [k for k in ("RNA_matrix_filtered", "max_cell_count") if re.search(rf"(?m)^{k}\s*:", text)]
    if len(have) == 2:
        return [f"{path.name}: RNA_matrix_filtered and max_cell_count already set"]
    if apply:
        add = TRAINING_YAML_ADDITION
        for k in have:                      # keep whichever one is already there
            add = "\n".join(ln for ln in add.splitlines() if not ln.startswith(f"{k}:")) + "\n"
        path.write_text(text.rstrip("\n") + "\n" + add)
    return [f"{path.name}: {'added' if apply else 'would add'} " + ", ".join(
        k for k in ("RNA_matrix_filtered", "max_cell_count") if k not in have)]


def ensure_feature_table(repo: Path, apply: bool, offline: bool = False) -> "list[str]":
    dest = repo / "resources" / "feature_tables" / "multiome_arc_n6.tsv"
    if dest.is_file():
        return [f"{relpath(dest, repo)}: present"]
    if not apply:
        return [f"{relpath(dest, repo)}: MISSING (setup would fetch it from kaybrand/scE2G)"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    text = None
    if not offline:
        try:
            with urllib.request.urlopen(FEATURE_TABLE_URL, timeout=30) as r:
                text = r.read().decode()
            if not text.startswith("feature\tinput_col"):
                text = None
        except Exception as exc:
            logging.warning("feature table fetch failed: %s", exc)
    src = "fetched from kaybrand/scE2G"
    if text is None:
        text, src = MULTIOME_ARC_N6, "written from the copy embedded in igvfagent"
    dest.write_text(text)
    return [f"{relpath(dest, repo)}: {src}"]


def git(repo: Path, *args: str) -> "tuple[int, str]":
    try:
        p = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, timeout=600)
        return p.returncode, (p.stdout + p.stderr).strip()
    except (OSError, subprocess.TimeoutExpired) as exc:
        return 1, str(exc)


def cmd_setup(args: argparse.Namespace) -> int:
    log_path = setup_logging()
    repo = Path(args.repo_dir).expanduser().resolve()
    notes: "list[str]" = []
    apply = not args.check
    if not (repo / "workflow").is_dir():
        if args.check:
            print(f"ERROR: {repo} is not an scE2G checkout (no workflow/ directory)", file=sys.stderr)
            return 2
        if shutil.which("git") is None:
            print("ERROR: git is not on PATH; clone https://github.com/EngreitzLab/scE2G yourself", file=sys.stderr)
            return 2
        print(f"Cloning scE2G ({args.branch}) into {repo} ...")
        p = subprocess.run(["git", "clone", "--recurse-submodules", "--branch", args.branch, SCE2G_REPO, str(repo)],
                           capture_output=True, text=True)
        if p.returncode != 0:
            print(f"ERROR: git clone failed:\n{p.stderr.strip()}", file=sys.stderr)
            return 1
        notes.append(f"cloned {SCE2G_REPO} @ {args.branch}")
    else:
        rc, head = git(repo, "rev-parse", "--abbrev-ref", "HEAD")
        rc2, sha = git(repo, "rev-parse", "--short", "HEAD")
        notes.append(f"checkout at {repo} on {head if rc == 0 else '?'} ({sha if rc2 == 0 else '?'})")
        if apply and head != args.branch and rc == 0 and not args.keep_branch:
            rc3, out = git(repo, "checkout", args.branch)
            notes.append(f"checkout {args.branch}: {'ok' if rc3 == 0 else 'FAILED: ' + out.splitlines()[-1]}")
        if apply and not (repo / "ENCODE-rE2G" / "workflow").is_dir():
            rc4, out = git(repo, "submodule", "update", "--init", "--recursive")
            notes.append(f"submodule init: {'ok' if rc4 == 0 else 'FAILED: ' + out.splitlines()[-1] if out else 'FAILED'}")
    notes += patch_snakefile(repo / "workflow" / "Snakefile_training", add_scripts_dir=True, apply=apply)
    notes += patch_snakefile(repo / "ENCODE-rE2G" / "workflow" / "Snakefile_training", add_scripts_dir=False, apply=apply)
    notes += patch_training_yaml(repo / "config" / "config_training.yaml", apply=apply)
    notes += ensure_feature_table(repo, apply=apply, offline=args.offline)
    for c in ("config/config_cell_clusters.tsv", "config/config_models.tsv", "models/multiome_powerlaw_v3"):
        notes.append(f"{c}: {'present' if (repo / c).exists() else 'MISSING'}")
    for n in notes:
        print(("  " if not n.startswith("cloned") else "") + n)
    problems = [n for n in notes if "MISSING" in n or "FAILED" in n or "WARNING" in n or n.startswith("missing")]
    print(f"Setup {'check' if args.check else 'done'}: {len(problems)} problem(s)")
    print(f"Log: {log_path}")
    return 1 if problems and not args.check else 0


# ---------------------------------------------------------------------------
# features: Synapse feature table -> source_file + external config + feature table
# ---------------------------------------------------------------------------

def sanitize_feature(name: str) -> str:
    return re.sub(r"\s+", "_", name.strip())


def convert_feature_table(src: Path, dest: Path, features: "Optional[list[str]]" = None) -> "dict":
    """ElementChr/Start/End/GeneSymbol + feature columns -> chr/start/end/TargetGene + renamed features."""
    pd = _pd()
    df = pd.read_csv(src, sep="\t", comment="#", low_memory=False)
    ren = {c: SOURCE_RENAME[c] for c in df.columns if c in SOURCE_RENAME}
    df = df.rename(columns=ren)
    missing = [c for c in SOURCE_REQUIRED if c not in df.columns]
    if missing:
        raise SystemExit(f"{src}: cannot find {missing} (have {list(df.columns)[:12]}...). Expected "
                         f"ElementChr/ElementStart/ElementEnd/GeneSymbol or chr/start/end/TargetGene.")
    drop = {"ElementName", "ElementClass", "isSelfPromoter", "TargetGeneEnsemblID", "TargetGeneTSS",
            "CellType", "Score", "DistanceToTSS", "class"}
    candidates = [c for c in df.columns if c not in SOURCE_REQUIRED and c not in drop]
    if features:
        want = set(features) | {sanitize_feature(f) for f in features}
        candidates = [c for c in candidates if c in want or sanitize_feature(c) in want]
        unknown = [f for f in features if f not in df.columns and sanitize_feature(f) not in
                   {sanitize_feature(c) for c in df.columns}]
        if unknown:
            raise SystemExit(f"features not in {src.name}: {unknown}")
    numeric = [c for c in candidates if pd.api.types.is_numeric_dtype(df[c]) or pd.api.types.is_bool_dtype(df[c])]
    skipped = [c for c in candidates if c not in numeric]
    mapping = {c: sanitize_feature(c) for c in numeric}
    out = df[SOURCE_REQUIRED + numeric].rename(columns=mapping)
    out["start"] = out["start"].astype("int64")
    out["end"] = out["end"].astype("int64")
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(dest, sep="\t", index=False, compression="gzip" if str(dest).endswith(".gz") else None)
    return {"rows": int(len(out)), "features": mapping, "skipped_non_numeric": skipped,
            "renamed_coords": ren, "genes": int(out["TargetGene"].nunique())}


def write_external_config(path: Path, source_file: str, features: "dict[str, str]",
                          aggregate: str = "mean", join_by: str = "overlap") -> Path:
    rows = [{"input_col": new, "source_col": new, "aggregate_function": aggregate, "join_by": join_by,
             "source_file": source_file} for _orig, new in features.items()]
    return write_tsv_rows(path, EXT_CONFIG_COLS, rows)


def write_feature_table(path: Path, base_table: "Optional[Path]", features: "dict[str, str]",
                        aggregate: str = "max", fill_value: str = "0") -> Path:
    if base_table is not None and base_table.is_file():
        header, rows = read_tsv_rows(base_table)
    else:
        rows = list(csv.DictReader(MULTIOME_ARC_N6.splitlines(), delimiter="\t"))
    existing = {r["feature"] for r in rows}
    for orig, new in features.items():
        if new in existing:
            continue
        rows.append({"feature": new, "input_col": new, "second_input": "NA", "aggregate_function": aggregate,
                     "fill_value": fill_value, "nice_name": orig.strip()})
    return write_tsv_rows(path, FEATURE_TABLE_COLS, rows)


def cmd_features(args: argparse.Namespace) -> int:
    log_path = setup_logging()
    repo = Path(args.repo_dir).expanduser().resolve() if args.repo_dir else None
    name = safe_label(args.name)
    src = Path(args.features).expanduser()
    if not src.is_file():
        print(f"ERROR: feature table not found: {src}", file=sys.stderr)
        return 2
    base = repo if repo else Path(args.out_dir or ".").resolve()
    source_file = base / "resources" / "external_features" / f"{name}.tsv.gz"
    ext_config = base / "config" / f"external_features_config_{name}.tsv"
    feat_table = base / "resources" / "feature_tables" / f"multiome_arc_n6_{name}.tsv"
    base_table = (repo / "resources" / "feature_tables" / "multiome_arc_n6.tsv") if repo else None
    info = convert_feature_table(src, source_file, [f for f in (args.select or "").split(",") if f] or None)
    if not info["features"]:
        print("ERROR: no numeric feature columns left after selection", file=sys.stderr)
        return 1
    write_external_config(ext_config, relpath(source_file, base) if repo else str(source_file), info["features"],
                          aggregate=args.merge_aggregate)
    write_feature_table(feat_table, base_table, info["features"], aggregate=args.benchmark_aggregate,
                        fill_value=args.fill_value)
    print(f"Source file: {source_file}  ({info['rows']:,} rows, {info['genes']:,} genes)")
    print("Features: " + ", ".join(f"{o} -> {n}" if o != n else n for o, n in info["features"].items()))
    if info["skipped_non_numeric"]:
        print("Skipped non-numeric columns: " + ", ".join(info["skipped_non_numeric"]))
    if info["renamed_coords"]:
        print("Renamed: " + ", ".join(f"{a} -> {b}" for a, b in info["renamed_coords"].items()))
    print(f"External features config: {ext_config}")
    print(f"Feature table: {feat_table}")
    print("Next: igvfagent sce2g configure --repo-dir <repo> --cluster <name> --rna ... --atac-frag ... "
          f"--model <model> --name {name}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# configure: cluster + model rows, run script
# ---------------------------------------------------------------------------

def upsert_row(path: Path, key_col: str, key: str, row: dict, base_cols: "list[str]",
               extra_cols: "Iterable[str]" = ()) -> "tuple[Path, str]":
    header, rows = read_tsv_rows(path) if path.is_file() else (list(base_cols), [])
    header = list(header) or list(base_cols)
    for c in list(base_cols) + list(extra_cols):
        if c not in header:
            header.append(c)
    action = "added"
    out = []
    for r in rows:
        if r.get(key_col) == key:
            action = "replaced"
            continue
        out.append(r)
    out.append(row)
    write_tsv_rows(path, header, out)
    return path, action


def cmd_configure(args: argparse.Namespace) -> int:
    log_path = setup_logging()
    repo = Path(args.repo_dir).expanduser().resolve()
    name = safe_label(args.name) if args.name else ""
    model = args.model or (f"multiome_{name}" if name else f"model_{args.cluster}")
    ext_cfg = (f"config/external_features_config_{name}.tsv" if name else "")
    feat_table = args.feature_table or (f"resources/feature_tables/multiome_arc_n6_{name}.tsv" if name
                                        else "resources/feature_tables/multiome_arc_n6.tsv")
    cluster_row = {"cluster": args.cluster, "rna_matrix_file": args.rna or "", "atac_frag_file": args.atac_frag,
                   "HiC_file": args.hic or "", "HiC_type": args.hic_type or "", "HiC_resolution": args.hic_resolution or "",
                   "alt_TSS": "", "alt_genes": "", "model_dir": args.model_dir, EXT_COL: ext_cfg}
    p1, a1 = upsert_row(repo / "config" / "config_cell_clusters.tsv", "cluster", args.cluster, cluster_row,
                        CLUSTER_COLS, [EXT_COL] if ext_cfg else [])
    model_row = {"model": model, "dataset": args.cluster, "ABC_directory": "", "feature_table": feat_table,
                 "polynomial": "False", "override_params": ""}
    p2, a2 = upsert_row(repo / "config" / "config_models.tsv", "model", model, model_row, MODEL_COLS)
    cmd = SNAKEMAKE_TRAINING.format(profile=args.profile) if args.profile else \
        SNAKEMAKE_TRAINING_LOCAL.format(jobs=args.jobs)
    script = repo / f"run_training_{safe_label(model)}.sh"
    script.write_text("#!/usr/bin/env bash\n# written by igvfagent sce2g configure\nset -euo pipefail\n"
                      f"cd \"$(dirname \"$0\")\"\n{cmd} \"$@\"\n")
    script.chmod(0o755)
    print(f"Cluster row {a1}: {args.cluster} -> {p1}")
    print(f"Model row {a2}: {model} (dataset {args.cluster}, feature_table {feat_table}) -> {p2}")
    if ext_cfg:
        print(f"External features config: {ext_cfg}")
    print(f"Run script: {script}")
    print(f"Command: {cmd}")
    print(f"Next: igvfagent sce2g check --repo-dir {repo} --model {model}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# check: do the paths and names line up?
# ---------------------------------------------------------------------------

def _yaml_scalars(path: Path) -> "dict[str, str]":
    """Top-level `key: value` pairs of a simple YAML file, without PyYAML."""
    out = {}
    if not path.is_file():
        return out
    for ln in path.read_text().splitlines():
        m = re.match(r'^([A-Za-z_][\w]*)\s*:\s*(.*?)\s*$', ln)
        if m and not m.group(2).startswith("["):
            out[m.group(1)] = m.group(2).strip().strip('"').strip("'")
    return out


def check_setup(repo: Path, model: "Optional[str]" = None) -> "tuple[list[str], list[str]]":
    """(problems, notes) for a training setup."""
    problems, notes = [], []
    cfg = repo / "config"
    yaml = _yaml_scalars(cfg / "config_training.yaml")
    if not yaml:
        return [f"config/config_training.yaml missing or unreadable"], notes
    for k in ("RNA_matrix_filtered", "max_cell_count"):
        (notes if k in yaml else problems).append(f"config_training.yaml {k}: {yaml.get(k, 'NOT SET')}")
    for k in ("gene_annotations", "gene_TSS500", "genes", "crispr_dataset"):
        v = yaml.get(k)
        if v and not (repo / v).exists():
            problems.append(f"config_training.yaml {k} -> {v}: file missing")
    sf = repo / "workflow" / "Snakefile_training"
    if sf.is_file():
        txt = sf.read_text()
        if re.search(r'(?m)^\s*conda:\s*["\']mamba["\']', txt):
            problems.append("workflow/Snakefile_training still has `conda: \"mamba\"` (run sce2g setup)")
        if "SCRIPTS_DIR" not in txt:
            problems.append("workflow/Snakefile_training lacks SCRIPTS_DIR (run sce2g setup)")
    else:
        problems.append("workflow/Snakefile_training missing")
    sf2 = repo / "ENCODE-rE2G" / "workflow" / "Snakefile_training"
    if sf2.is_file() and re.search(r'(?m)^\s*conda:\s*["\']mamba["\']', sf2.read_text()):
        problems.append("ENCODE-rE2G/workflow/Snakefile_training still has `conda: \"mamba\"`")
    if not sf2.is_file():
        problems.append("ENCODE-rE2G submodule not initialised (git submodule update --init --recursive)")

    ch, clusters = read_tsv_rows(cfg / "config_cell_clusters.tsv") if (cfg / "config_cell_clusters.tsv").is_file() else ([], [])
    mh, models = read_tsv_rows(cfg / "config_models.tsv") if (cfg / "config_models.tsv").is_file() else ([], [])
    if not clusters:
        problems.append("config/config_cell_clusters.tsv has no rows")
    if not models:
        problems.append("config/config_models.tsv has no rows")
    if models and model:
        models = [m for m in models if m.get("model") == model] or problems.append(f"model {model!r} not in config_models.tsv") or []
    by_cluster = {c.get("cluster"): c for c in clusters}
    for m in models:
        ds = m.get("dataset", "")
        c = by_cluster.get(ds)
        if c is None:
            problems.append(f"model {m.get('model')}: dataset {ds!r} has no cluster row")
            continue
        notes.append(f"model {m.get('model')} <- cluster {ds}")
        for col in ("rna_matrix_file", "atac_frag_file"):
            v = c.get(col, "")
            if not v:
                (problems if col == "atac_frag_file" else notes).append(f"cluster {ds}: {col} empty")
            elif not (repo / v).exists() and not Path(v).exists():
                problems.append(f"cluster {ds}: {col} -> {v}: file missing")
        for md in (c.get("model_dir") or "").split(","):
            md = md.strip()
            if md and not (repo / md).is_dir():
                problems.append(f"cluster {ds}: model_dir {md} missing")
        if str(m.get("polynomial", "")).strip() not in ("False", "TRUE", "True", "FALSE", ""):
            problems.append(f"model {m.get('model')}: polynomial should be True/False")
        if m.get("ABC_directory"):
            notes.append(f"model {m.get('model')}: ABC_directory set ({m['ABC_directory']}); the walkthrough leaves it blank")
        ft = m.get("feature_table", "")
        ftp = repo / ft
        if not ft or not ftp.is_file():
            problems.append(f"model {m.get('model')}: feature_table {ft!r} missing")
            continue
        fth, ftrows = read_tsv_rows(ftp)
        if fth != FEATURE_TABLE_COLS:
            problems.append(f"{ft}: columns {fth} != {FEATURE_TABLE_COLS}")
        spaced = [r["feature"] for r in ftrows if " " in (r.get("feature") or "") or " " in (r.get("input_col") or "")]
        if spaced:
            problems.append(f"{ft}: feature names with spaces: {spaced}")
        ext = c.get(EXT_COL, "")
        if not ext:
            notes.append(f"cluster {ds}: no {EXT_COL} (plain multiome model, no crowdsourced features)")
            continue
        extp = repo / ext
        if not extp.is_file():
            problems.append(f"cluster {ds}: {EXT_COL} -> {ext}: file missing")
            continue
        eh, erows = read_tsv_rows(extp)
        if eh != EXT_CONFIG_COLS:
            problems.append(f"{ext}: columns {eh} != {EXT_CONFIG_COLS}")
        ft_inputs = {r.get("input_col") for r in ftrows} | {r.get("second_input") for r in ftrows}
        for r in erows:
            ic = r.get("input_col", "")
            if ic not in ft_inputs:
                problems.append(f"{ext}: input_col {ic!r} is not in {ft} (add a feature row or the model never sees it)")
            if r.get("join_by") not in ("overlap", "TargetGene"):
                problems.append(f"{ext}: join_by {r.get('join_by')!r} must be overlap or TargetGene")
            sfile = r.get("source_file", "")
            sp = repo / sfile
            if not sfile or not (sp.is_file() or Path(sfile).is_file()):
                problems.append(f"{ext}: source_file {sfile!r} missing")
                continue
            sp = sp if sp.is_file() else Path(sfile)
            with _open_text(sp) as fh:
                head = fh.readline().rstrip("\n").split("\t")
            need = SOURCE_REQUIRED if r.get("join_by") == "overlap" else ["TargetGene"]
            for col in need + [r.get("source_col", "")]:
                if col not in head:
                    problems.append(f"{sfile}: column {col!r} missing (has {head[:8]}...)")
            if any(" " in h for h in head):
                problems.append(f"{sfile}: header has spaces: {[h for h in head if ' ' in h]}")
        notes.append(f"cluster {ds}: {len(erows)} external feature(s) via {ext}")
    return problems, notes


def cmd_check(args: argparse.Namespace) -> int:
    log_path = setup_logging()
    repo = Path(args.repo_dir).expanduser().resolve()
    problems, notes = check_setup(repo, args.model)
    for n in notes:
        print(f"  ok    {n}")
    for p in problems:
        print(f"  FAIL  {p}")
    print(f"Check: {len(problems)} problem(s)")
    print("Command: " + (SNAKEMAKE_TRAINING.format(profile=args.profile) if args.profile else
                         SNAKEMAKE_TRAINING_LOCAL.format(jobs=args.jobs)))
    print(f"Log: {log_path}")
    return 1 if problems else 0


# ---------------------------------------------------------------------------
# run: snakemake, dry-run by default
# ---------------------------------------------------------------------------

def cmd_run(args: argparse.Namespace) -> int:
    log_path = setup_logging()
    repo = Path(args.repo_dir).expanduser().resolve()
    problems, _ = check_setup(repo, args.model)
    if problems and not args.force:
        for p in problems:
            print(f"  FAIL  {p}")
        print("Refusing to run with an inconsistent setup (use --force to override).", file=sys.stderr)
        return 1
    if shutil.which("snakemake") is None:
        print("ERROR: snakemake is not on PATH. Activate the scE2G environment "
              "(conda/mamba env with snakemake >= 7) and re-run.", file=sys.stderr)
        return 2
    cmd = (SNAKEMAKE_TRAINING.format(profile=args.profile) if args.profile else
           SNAKEMAKE_TRAINING_LOCAL.format(jobs=args.jobs)).split()
    if args.dry_run:
        cmd.append("-n")
    if args.extra:
        cmd += args.extra
    print("Running: " + " ".join(cmd))
    out_dir = run_dir(f"train_{args.model or 'all'}")
    log_file = out_dir / "snakemake.log"
    with open(log_file, "w") as fh:
        p = subprocess.Popen(cmd, cwd=str(repo), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        assert p.stdout is not None
        tail: "list[str]" = []
        for line in p.stdout:
            fh.write(line)
            tail.append(line.rstrip("\n"))
            tail = tail[-25:]
        rc = p.wait()
    for ln in tail:
        print(ln)
    print(f"snakemake exit {rc} ({'dry run' if args.dry_run else 'real run'})")
    print(f"Log: {log_file}")
    print(f"Log: {log_path}")
    return rc


# ---------------------------------------------------------------------------
# predictions: describe / compare scE2G output tables
# ---------------------------------------------------------------------------

def _pred_label(spec: str) -> "tuple[str, Path]":
    label, _, p = spec.partition("=")
    if not p:
        return Path(spec).name.split(".")[0][:24], Path(spec)
    return label, Path(p)


def _gene_col(df) -> str:
    for c in ("GeneSymbol", "TargetGene", "gene", "measuredGeneSymbol"):
        if c in df.columns:
            return c
    raise SystemExit(f"no gene column among {list(df.columns)[:12]}")


def _coord_cols(df) -> "tuple[str, str, str]":
    for trio in (("ElementChr", "ElementStart", "ElementEnd"), ("chr", "start", "end"), ("chrom", "chromStart", "chromEnd")):
        if all(c in df.columns for c in trio):
            return trio
    raise SystemExit(f"no element coordinate columns among {list(df.columns)[:12]}")


def cmd_predictions(args: argparse.Namespace) -> int:
    pd = _pd()
    log_path = setup_logging()
    out_dir = run_dir(args.label or "predictions")
    tables = {}
    rows = []
    for spec in args.predictions:
        label, p = _pred_label(spec)
        df = pd.read_csv(p, sep="\t", comment="#", low_memory=False)
        g = _gene_col(df)
        cc = _coord_cols(df)
        score = args.score_col if args.score_col in df.columns else next(
            (c for c in ("E2G.Score.qnorm", "Score", "E2G.Score", "ABC.Score") if c in df.columns), None)
        df["_pair"] = df[cc[0]].astype(str) + ":" + df[cc[1]].astype(str) + "-" + df[cc[2]].astype(str) + "|" + df[g].astype(str)
        tables[label] = (df, score)
        r = {"table": label, "links": len(df), "genes": int(df[g].nunique()),
             "elements": int((df[cc[0]].astype(str) + ":" + df[cc[1]].astype(str) + "-" + df[cc[2]].astype(str)).nunique()),
             "score_col": score or ""}
        if score:
            s = pd.to_numeric(df[score], errors="coerce")
            r.update({"score_min": float(s.min()), "score_median": float(s.median()), "score_max": float(s.max()),
                      f"links_ge_{args.threshold}": int((s >= args.threshold).sum())})
        if "ElementClass" in df.columns:
            r["element_classes"] = "; ".join(f"{k}={v}" for k, v in df["ElementClass"].value_counts().items())
        if "isSelfPromoter" in df.columns:
            r["self_promoter_links"] = int(df["isSelfPromoter"].astype(str).str.lower().eq("true").sum())
        per_gene = df.groupby(g).size()
        r["links_per_gene_median"] = float(per_gene.median())
        rows.append(r)
    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "predictions_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    labels = list(tables)
    if len(labels) >= 2:
        cmp_rows = []
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                a, sa = tables[labels[i]]
                b, sb = tables[labels[j]]
                pa, pb = set(a["_pair"]), set(b["_pair"])
                shared = pa & pb
                row = {"a": labels[i], "b": labels[j], "pairs_a": len(pa), "pairs_b": len(pb), "shared_pairs": len(shared),
                       "jaccard": round(len(shared) / max(1, len(pa | pb)), 4),
                       "shared_genes": len(set(a[_gene_col(a)]) & set(b[_gene_col(b)]))}
                if sa and sb and shared:
                    m = a.drop_duplicates("_pair").set_index("_pair")[sa].to_frame("x").join(
                        b.drop_duplicates("_pair").set_index("_pair")[sb].to_frame("y"), how="inner")
                    row["spearman_shared"] = round(float(m["x"].rank().corr(m["y"].rank())), 4) if len(m) > 2 else float("nan")
                cmp_rows.append(row)
        cmp = pd.DataFrame(cmp_rows)
        cmp.to_csv(out_dir / "predictions_pairwise.tsv", sep="\t", index=False)
        print(cmp.to_string(index=False))
        print(f"CSV: {out_dir / 'predictions_pairwise.tsv'}")
    plt = _plt()
    if plt and not args.no_plots:
        fig, ax = plt.subplots(figsize=(8, 4))
        for k, label in enumerate(labels):
            df, score = tables[label]
            if score:
                ax.hist(pd.to_numeric(df[score], errors="coerce").dropna(), bins=50, histtype="step", lw=1.5,
                        color=SERIES[k % len(SERIES)], label=label)
        ax.axvline(args.threshold, color=INK2, ls="--", lw=0.8)
        _style(ax, "scE2G score distributions", "score", "links")
        if len(labels) > 1:
            ax.legend(fontsize=8, frameon=False)
        fig.patch.set_facecolor(SURFACE)
        fig.tight_layout()
        fig.savefig(out_dir / "score_distributions.png", dpi=140)
        plt.close(fig)
        print(f"Figure: {out_dir / 'score_distributions.png'}")
    print(f"CSV: {out_dir / 'predictions_summary.tsv'}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# benchmark: CRISPR_comparison semantics, clean-room
# ---------------------------------------------------------------------------

def read_crispr(path: Path) -> "list[dict]":
    """EPCrisprBenchmark rows as chrom/start/end/gene/regulated (+ CellType when present)."""
    from e2g_benchmark_eval import read_ground_truth  # noqa: E402
    truth = read_ground_truth(str(path))
    with _open_text(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    if rows and "CellType" in rows[0] and len(rows) == len(truth):
        for t, r in zip(truth, rows):
            t["CellType"] = r.get("CellType", "")
    return truth


def parse_pred_config(path: "Optional[Path]", predictions: "list[str]", default_col: str) -> "list[dict]":
    """pred_config.txt rows, or defaults (max / fill 0 / not inverse) per table."""
    if path and path.is_file():
        _h, rows = read_tsv_rows(path)
        out = []
        for r in rows:
            if str(r.get("include", "TRUE")).upper() in ("FALSE", "0"):
                continue
            out.append({"pred_id": r["pred_id"], "pred_col": r["pred_col"],
                        "boolean": str(r.get("boolean", "FALSE")).upper() == "TRUE",
                        "alpha": float(r["alpha"]) if r.get("alpha") not in (None, "", "NA") else None,
                        "aggregate_function": (r.get("aggregate_function") or "max").lower(),
                        "fill_value": float(r.get("fill_value") or 0),
                        "inverse_predictor": str(r.get("inverse_predictor", "FALSE")).upper() == "TRUE",
                        "pred_name_long": r.get("pred_name_long") or r["pred_id"],
                        "color": r.get("color") or ""})
        return out
    out = []
    for k, spec in enumerate(predictions):
        label, _ = _pred_label(spec)
        out.append({"pred_id": label, "pred_col": default_col, "boolean": False, "alpha": None,
                    "aggregate_function": "max", "fill_value": 0.0, "inverse_predictor": False,
                    "pred_name_long": label, "color": SERIES[k % len(SERIES)]})
    return out


META_COLS = {"ElementChr", "ElementStart", "ElementEnd", "ElementName", "ElementClass", "GeneTSS", "GeneSymbol",
             "GeneEnsemblID", "SampleSummaryShort", "isSelfPromoter", "chr", "start", "end", "TargetGene", "chrom",
             "chromStart", "chromEnd", "TargetGeneEnsemblID", "TargetGeneTSS", "CellType", "name", "class", "gene"}
SHARED_UNIVERSE_COLS = {"E2G_Distance", "GeneTSS"}
_INVERSE_RE = re.compile(r"(pval|p_value|pvalue|fdr|qval|distance|loeuf|rank)", re.I)
_PRED_CACHE: "dict[str, Any]" = {}


def load_predictions(path: Path, genes: "Optional[set]" = None, chunksize: int = 1_000_000):
    """A prediction table restricted to CRISPR-tested genes, read in chunks.

    The crowdsourced K562 feature tables are 11 million element-gene rows
    each (EPCOT's is 2 GB compressed with ~100 columns); only rows whose
    gene is in the CRISPR set can ever be scored, and that is ~5%.
    """
    pd = _pd()
    key = f"{path.resolve()}|{len(genes) if genes else 'all'}"
    if key in _PRED_CACHE:
        return _PRED_CACHE[key]
    frames = []
    for chunk in pd.read_csv(path, sep="\t", comment="#", low_memory=False, chunksize=chunksize):
        if genes is not None:
            g = _gene_col(chunk)
            chunk = chunk[chunk[g].astype(str).isin(genes)]
        if len(chunk):
            frames.append(chunk)
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    _PRED_CACHE[key] = df
    return df


def feature_columns(df) -> "list[str]":
    """Numeric / boolean columns that are features rather than coordinates or ids."""
    pd = _pd()
    out = []
    for c in df.columns:
        if c in META_COLS or c.startswith("_"):
            continue
        if pd.api.types.is_numeric_dtype(df[c]) or pd.api.types.is_bool_dtype(df[c]):
            out.append(c)
        elif df[c].dtype == object:
            low = df[c].dropna().astype(str).str.lower().unique()[:5]
            if len(low) and set(low) <= {"true", "false", "0", "1"}:
                out.append(c)
    return out


def score_crispr_pairs(truth: "list[dict]", pred_df, pred_col: str, aggregate: str, fill_value: float,
                       inverse: bool, boolean: bool) -> "tuple[list[float], int]":
    """One score per CRISPR pair: aggregate of overlapping predicted elements for the same gene."""
    from bisect import bisect_left
    pd = _pd()
    g = _gene_col(pred_df)
    cc = _coord_cols(pred_df)
    if pred_col not in pred_df.columns:
        raise SystemExit(f"pred_col {pred_col!r} not in prediction table (have {list(pred_df.columns)[:12]}...)")
    vals = pd.to_numeric(pred_df[pred_col], errors="coerce")
    if boolean:
        vals = pred_df[pred_col].astype(str).str.lower().isin(("true", "1", "t", "yes")).astype(float)
    idx: "dict[tuple[str, str], list]" = {}
    for gene, chrom, st, en, v in zip(pred_df[g].astype(str), pred_df[cc[0]].astype(str), pred_df[cc[1]],
                                       pred_df[cc[2]], vals):
        if v != v:
            continue
        idx.setdefault((gene, chrom.replace("chr", "")), []).append((int(st), int(en), float(v)))
    for k in idx:
        idx[k].sort()
    # Inverse predictors (distance, p-values) are negated BEFORE aggregation,
    # so `max` picks the closest / most significant overlapping element, and
    # a tested pair with no overlapping element gets the WORST score, not
    # the fill value negated. The first version negated after aggregating
    # and filling: every unmatched pair scored -0, which out-ranked every
    # real (negative) distance and put the distance baseline at AUPRC 0.106
    # where CRISPR_comparison's own distance baseline scores 0.436.
    sign = -1.0 if inverse else 1.0
    scores, matched, unmatched_idx = [], 0, []
    for t in truth:
        ivs = idx.get((t["gene"], str(t["chrom"]).replace("chr", "")), [])
        starts = [i[0] for i in ivs]
        hi = bisect_left(starts, t["end"] + 1)
        hits = [sign * v for st, en, v in ivs[:hi] if en >= t["start"] and st <= t["end"]]
        if not hits:
            unmatched_idx.append(len(scores))
            scores.append(sign * fill_value)
            continue
        matched += 1
        if aggregate == "mean":
            s = sum(hits) / len(hits)
        elif aggregate == "sum":
            s = sum(hits)
        elif aggregate == "min":
            s = min(hits)
        else:
            s = max(hits)
        scores.append(s)
    if inverse and unmatched_idx:
        worst = min(scores[i] for i in range(len(scores)) if i not in set(unmatched_idx)) if matched else 0.0
        for i in unmatched_idx:
            scores[i] = worst - 1.0
    return scores, matched


def pr_curve(labels: "list[bool]", scores: "list[float]", min_recall: float = 0.7
             ) -> "tuple[float, float, dict]":
    """Tie-aware precision-recall: (AUPRC, precision at first recall >= min_recall, curve).

    One threshold per DISTINCT score, as ROCR (CRISPR_comparison) and
    scikit-learn do. Walking tied pairs one at a time in arbitrary order --
    what e2g_benchmark_eval._manual_pr does -- is fine for a model score
    with 10,000 distinct values and wrong for a feature like Pinloop where
    3,207 of 10,356 pairs share the value 0: the order chosen inside the tie
    moved its AUPRC between 0.351 and 0.365 while the tie-aware value is 0.357.
    AUPRC is the step integral sum(precision_k * (recall_k - recall_{k-1})),
    scikit-learn's average_precision convention.
    """
    n_pos = sum(1 for y in labels if y)
    if n_pos == 0 or not scores:
        return 0.0, 0.0, {"precision": [], "recall": []}
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    tp = fp = 0
    prev_recall = 0.0
    ap = 0.0
    p_at = None
    pcs, rcs = [], []
    i = 0
    while i < len(order):
        j = i
        while j < len(order) and scores[order[j]] == scores[order[i]]:
            if labels[order[j]]:
                tp += 1
            else:
                fp += 1
            j += 1
        prec = tp / (tp + fp)
        rec = tp / n_pos
        ap += prec * (rec - prev_recall)
        prev_recall = rec
        pcs.append(prec)
        rcs.append(rec)
        if p_at is None and rec >= min_recall:
            p_at = prec
        i = j
    # CRISPR_comparison integrates the same tie-aware points with the
    # trapezoid rule (caTools::trapz over ROCR's curve, NaN start dropped);
    # scikit-learn and e2g_benchmark_eval use the step rule. Both are kept:
    # `ap` (step) for continuity with benchmark #12, `trapz` for the
    # number-for-number comparison with the upstream pipeline.
    trapz = 0.0
    for k in range(1, len(pcs)):
        trapz += (rcs[k] - rcs[k - 1]) * (pcs[k] + pcs[k - 1]) / 2.0
    # CRISPR_comparison's exact definition (crisprComparisonBootstrapFunctions.R,
    # calculate_auprc): yardstick::pr_curve, drop the first row (threshold Inf,
    # recall 0) and the LAST row -- "ROCR always includes a Recall = 100% point
    # even if the predictor cannot achieve it" -- then caTools::trapz. The last
    # block is where every fill-value pair sits, so for a genome-wide table this
    # removes the final low-precision rise to recall 1 and reads a little lower
    # than `trapz`; per-pair scores were verified identical to the upstream
    # merged table, so this is the number to compare with its
    # performance_summary.txt.
    trapz_cc = 0.0
    for k in range(1, len(pcs) - 1):
        trapz_cc += (rcs[k] - rcs[k - 1]) * (pcs[k] + pcs[k - 1]) / 2.0
    if len(pcs) < 3:
        # A binary feature has two thresholds; dropping first and last leaves
        # nothing to integrate. Upstream excludes boolean predictors from its
        # PR summary for the same reason; report NaN, not 0.
        trapz_cc = float("nan")
    step = max(1, len(pcs) // 400)
    return ap, (p_at if p_at is not None else 0.0), {"precision": pcs[::step], "recall": rcs[::step],
                                                     "auprc_trapz": trapz, "auprc_crispr_comparison": trapz_cc}


def bootstrap_auprc(labels: "list[bool]", scores: "list[float]", n: int, seed: int = 1) -> "tuple[float, float]":
    _manual_pr = pr_curve
    rng = random.Random(seed)
    N = len(labels)
    vals = []
    for _ in range(n):
        pick = [rng.randrange(N) for _ in range(N)]
        ap, _p70, _c = _manual_pr([labels[i] for i in pick], [scores[i] for i in pick])
        vals.append(ap)
    vals.sort()
    return vals[int(0.025 * (n - 1))], vals[int(0.975 * (n - 1))]



# ---------------------------------------------------------------------------
# CRISPR_comparison extras: baseline predictors, TSS filter, distance bins,
# delta AUPRC (crisprComparisonSimplePredictors.R, crisprComparisonLoadInputData.R
# filterPredictionsTSS, crisprComparisonBootstrapFunctions.R bootstrapDeltaPerformance)
# ---------------------------------------------------------------------------

BASELINES_ALL = ("distToTSS", "distToGene", "nearestTSS", "nearestGene", "within100kbTSS", "within100kbGene",
                 "nearestExprTSS", "nearestExprGene", "within100kbExprTSS", "within100kbExprGene")
BASELINE_INVERSE = {"distToTSS", "distToGene"}


def read_bed_annot(path: Path) -> "list[tuple[str, int, int, str]]":
    """chr, start, end, name from a BED(6); header / track lines skipped."""
    out = []
    with _open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser", "chr\t")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                out.append((f[0], int(f[1]), int(f[2]), f[3]))
            except ValueError:
                continue
    return out


def _center_1bp(start: int, end: int) -> int:
    """GenomicRanges::resize(width = 1, fix = "center") on [start, end] read as 1-based closed."""
    return start + (end - start) // 2


def _gr_distance(a0: int, a1: int, b0: int, b1: int) -> int:
    """GenomicRanges::distance between closed ranges: gap in bp, 0 when overlapping or adjacent."""
    return max(b0 - a1 - 1, a0 - b1 - 1, 0)


def compute_baseline(name: str, truth: "list[dict]", tss: "list[tuple]", genes: "Optional[list[tuple]]",
                     expressed: "Optional[set]" = None) -> "list[float]":
    """One value per CRISPR pair, as computeBaselinePreds; NaN when the gene has no annotation."""
    use_gene = name.endswith("Gene")
    annot = genes if use_gene else tss
    if annot is None:
        raise SystemExit(f"baseline {name} needs --gene-bed")
    if "Expr" in name:
        if expressed is None:
            raise SystemExit(f"baseline {name} needs --expressed-genes")
        annot = [a for a in annot if a[3] in expressed]
    by_name: "dict[str, tuple]" = {}
    for a in annot:
        by_name.setdefault(a[3], a)
    by_chr: "dict[str, list]" = {}
    for a in annot:
        by_chr.setdefault(a[0], []).append(a)
    for k in by_chr:
        by_chr[k].sort(key=lambda a: a[1])
    out = []
    for t in truth:
        chrom = t["chrom"] if str(t["chrom"]).startswith("chr") else f"chr{t['chrom']}"
        if name.startswith("distTo"):
            a = by_name.get(t["gene"])
            if a is None:
                out.append(float("nan"))
                continue
            e = _center_1bp(t["start"], t["end"])
            if use_gene:
                out.append(float(_gr_distance(e, e, a[1], a[2])))
            else:
                c = _center_1bp(a[1], a[2])
                out.append(float(_gr_distance(e, e, c, c)))
        elif name.startswith("nearest"):
            best, bd = None, None
            for a in by_chr.get(chrom, []):
                d = _gr_distance(t["start"], t["end"], a[1], a[2])
                if bd is None or d < bd:
                    best, bd = a, d
            out.append(1.0 if best is not None and best[3] == t["gene"] else 0.0)
        else:  # within100kb*: CRE resized to 200 kb around its centre, any overlap with the gene's feature
            c = _center_1bp(t["start"], t["end"])
            lo, hi = c - 100_000 + 1, c + 100_000
            hit = any(a[3] == t["gene"] and a[1] <= hi and a[2] >= lo for a in by_chr.get(chrom, []))
            out.append(1.0 if hit else 0.0)
    return out


def filter_predictions_tss(df, tss: "list[tuple]"):
    """Drop predicted elements overlapping any TSS window (filterPredictionsTSS, 0-based TSS starts)."""
    from bisect import bisect_left
    cc = _coord_cols(df)
    by_chr: "dict[str, list]" = {}
    for a in tss:
        by_chr.setdefault(a[0].replace("chr", ""), []).append((a[1] + 1, a[2]))
    for k in by_chr:
        by_chr[k].sort()
    starts = {k: [x[0] for x in v] for k, v in by_chr.items()}
    maxlen = {k: max((x[1] - x[0] for x in v), default=0) for k, v in by_chr.items()}
    keep = []
    for chrom, st, en in zip(df[cc[0]].astype(str), df[cc[1]], df[cc[2]]):
        k = chrom.replace("chr", "")
        ivs = by_chr.get(k)
        if not ivs:
            keep.append(True)
            continue
        st1, en1 = int(st) + 1, int(en)
        i = bisect_left(starts[k], st1 - maxlen[k])
        hit = False
        while i < len(ivs) and ivs[i][0] <= en1:
            if ivs[i][1] >= st1:
                hit = True
                break
            i += 1
        keep.append(not hit)
    return df[keep]


def bootstrap_delta_auprc(labels: "list[bool]", s1: "list[float]", s2: "list[float]", n: int, seed: int = 1
                          ) -> "tuple[float, float, float, float]":
    """delta = AUPRC(s1) - AUPRC(s2) on the same resampled pairs; percentile CI and two-sided p-value
    (smallest alpha whose 1-alpha percentile interval excludes 0, as boot.pval type="perc")."""
    def auc(lab, sc):
        return pr_curve(lab, sc)[2]["auprc_crispr_comparison"]
    d0 = auc(labels, s1) - auc(labels, s2)
    rng = random.Random(seed)
    N = len(labels)
    ds = []
    for _ in range(n):
        pick = [rng.randrange(N) for _ in range(N)]
        lab = [labels[i] for i in pick]
        ds.append(auc(lab, [s1[i] for i in pick]) - auc(lab, [s2[i] for i in pick]))
    ds = sorted(x for x in ds if x == x)
    if not ds:
        return d0, float("nan"), float("nan"), float("nan")
    lo, hi = ds[int(0.025 * (len(ds) - 1))], ds[int(0.975 * (len(ds) - 1))]
    below = sum(1 for x in ds if x <= 0) / len(ds)
    above = sum(1 for x in ds if x >= 0) / len(ds)
    p = min(1.0, 2 * min(below, above))
    return d0, lo, hi, max(p, 1.0 / len(ds))


def write_crispr_comparison_configs(out_dir: Path, preds: "list[dict]", pred_paths: "dict[str, Path]",
                                    crispr: Path, name: str) -> "tuple[Path, Path]":
    rows = []
    for k, p in enumerate(preds):
        rows.append({"pred_id": p["pred_id"], "pred_col": p["pred_col"], "boolean": str(p["boolean"]).upper(),
                     "alpha": "" if p["alpha"] is None else p["alpha"], "aggregate_function": p["aggregate_function"],
                     "fill_value": p["fill_value"], "inverse_predictor": str(p["inverse_predictor"]).upper(),
                     "pred_name_long": p["pred_name_long"], "color": p["color"] or SERIES[k % len(SERIES)]})
    pc = write_tsv_rows(out_dir / "pred_config.txt", PRED_CONFIG_COLS, rows)
    lines = ["# written by igvfagent sce2g benchmark -- drop into CRISPR_comparison/config/config.yml",
             "comparisons:", f"  {name}:", "    pred:"]
    for p in preds:
        lines.append(f"      {p['pred_id']}: \"{pred_paths[p['pred_id']]}\"")
    lines += [f"    expt:", f"      crispr: \"{crispr}\"", f"    pred_config: \"{pc}\""]
    cy = out_dir / "config.yml"
    cy.write_text("\n".join(lines) + "\n")
    return pc, cy


def cmd_benchmark(args: argparse.Namespace) -> int:
    pd = _pd()
    _manual_pr = pr_curve
    log_path = setup_logging()
    out_dir = run_dir(args.label or "benchmark")
    truth = read_crispr(Path(args.crispr).expanduser())
    if args.cell_type:
        before = len(truth)
        truth = [t for t in truth if t.get("CellType") == args.cell_type]
        print(f"CellType {args.cell_type}: {len(truth)} / {before} CRISPR pairs kept")
    labels = [t["regulated"] for t in truth]
    genes = {t["gene"] for t in truth}
    pred_paths = {_pred_label(s)[0]: _pred_label(s)[1].expanduser().resolve() for s in args.predictions}
    if args.all_features:
        # One predictor per feature column of every table: the crowdsourced-
        # feature benchmark. Direction is unknown, so p-value / FDR / distance
        # style names are inverted and every row also reports the AUPRC of
        # the negated score.
        preds = []
        seen_shared: "set[str]" = set()
        for label, path in pred_paths.items():
            df = load_predictions(path, genes)
            for k, col in enumerate(feature_columns(df)):
                if col in SHARED_UNIVERSE_COLS:
                    # E2G_Distance is carried by every crowdsourced table; one
                    # copy is the distance baseline, the others are noise.
                    if col in seen_shared:
                        continue
                    seen_shared.add(col)
                    label_here = "baseline"
                else:
                    label_here = label
                preds.append({"pred_id": f"{label_here}:{col}", "pred_col": col, "boolean": False, "alpha": None,
                              "aggregate_function": "max", "fill_value": 0.0,
                              "inverse_predictor": bool(_INVERSE_RE.search(col)), "pred_name_long": col,
                              "color": SERIES[k % len(SERIES)], "_table": label})
        print(f"--all-features: {len(preds)} feature column(s) across {len(pred_paths)} table(s)")
    else:
        preds = parse_pred_config(Path(args.pred_config) if args.pred_config else None, args.predictions, args.score_col)
    results, curves = [], {}
    tss_filter = read_bed_annot(Path(args.filter_pred_tss).expanduser()) if getattr(args, "filter_pred_tss", None) else None
    _tss_cache: "dict[str, Any]" = {}
    pair_scores: "dict[str, dict]" = {}
    for p in preds:
        path = pred_paths.get(p.get("_table") or p["pred_id"])
        if path is None:
            print(f"WARNING: pred_config row {p['pred_id']} has no matching --predictions LABEL=PATH; skipped")
            continue
        df = load_predictions(path, genes)
        if tss_filter is not None and len(df):
            key = str(path)
            if key not in _tss_cache:
                before = len(df)
                _tss_cache[key] = filter_predictions_tss(df, tss_filter)
                print(f"filter-pred-tss: {path.name}: {before - len(_tss_cache[key]):,} of {before:,} rows overlap a TSS and are dropped")
            df = _tss_cache[key]
        truth_p, labels_p, n_missing = truth, labels, 0
        if args.gene_universe_filter:
            # CRISPR_comparison's filterExptGeneUniverse(): tested pairs whose
            # gene is absent from the predictor's gene universe are DROPPED
            # before scoring, not filled. Without this flag they score the
            # fill value and count against the predictor.
            universe = set(df[_gene_col(df)].astype(str)) if len(df) else set()
            truth_p = [t for t in truth if t["gene"] in universe]
            labels_p = [t["regulated"] for t in truth_p]
            n_missing = len(truth) - len(truth_p)
        scores, matched = score_crispr_pairs(truth_p, df, p["pred_col"], p["aggregate_function"], p["fill_value"],
                                             p["inverse_predictor"], p["boolean"])
        ap, p70, curve = _manual_pr(labels_p, scores)
        ap_inv = _manual_pr(labels_p, [-x for x in scores])[0]
        lo, hi = bootstrap_auprc(labels_p, scores, args.bootstrap) if args.bootstrap > 0 else (float("nan"), float("nan"))
        n_pos = sum(labels_p)
        results.append({"pred_id": p["pred_id"], "pred_col": p["pred_col"], "aggregate": p["aggregate_function"],
                        "fill_value": p["fill_value"], "inverse": p["inverse_predictor"],
                        "n_pairs": len(labels_p), "n_positive": int(n_pos), "baseline_precision": round(n_pos / max(1, len(labels_p)), 4),
                        "n_pairs_dropped_missing_gene": int(n_missing),
                        "pairs_overlapping_prediction": matched,
                        "frac_overlapping": round(matched / max(1, len(labels_p)), 4),
                        "auprc": round(ap, 4), "auprc_ci95_low": round(lo, 4), "auprc_ci95_high": round(hi, 4),
                        "auprc_trapz": round(curve.get("auprc_trapz", float("nan")), 4),
                        "auprc_crispr_comparison": round(curve.get("auprc_crispr_comparison", float("nan")), 4),
                        "precision_at_70_recall": round(p70, 4), "auprc_negated": round(ap_inv, 4),
                        "n_nonmissing_scores": int(sum(1 for x in scores if x != p["fill_value"] and x != -p["fill_value"]))})
        curves[p["pred_id"]] = curve
        pair_scores[p["pred_id"]] = {f"{t['chrom']}:{t['start']}-{t['end']}|{t['gene']}": sc for t, sc in zip(truth_p, scores)}
        if not args.all_features or args.keep_scored_pairs:
            pd.DataFrame({"pair": [f"{t['chrom']}:{t['start']}-{t['end']}|{t['gene']}" for t in truth_p],
                          "regulated": labels_p, "score": scores}).to_csv(out_dir / f"scored_pairs_{safe_label(p['pred_id'])}.tsv",
                                                                           sep="\t", index=False)
    # ---- baseline predictors, computed per CRISPR pair (no overlap, no fill)
    tss_annot = read_bed_annot(Path(args.tss_bed).expanduser()) if getattr(args, "tss_bed", None) else None
    gene_annot = read_bed_annot(Path(args.gene_bed).expanduser()) if getattr(args, "gene_bed", None) else None
    expressed = None
    if getattr(args, "expressed_genes", None):
        _h, erows = read_tsv_rows(Path(args.expressed_genes).expanduser())
        expressed = {r.get("gene") for r in erows if str(r.get("expressed", "TRUE")).upper() == "TRUE"
                     and (not args.cell_type or r.get("cell_type", args.cell_type) == args.cell_type)}
    key_of = lambda t: f"{t['chrom']}:{t['start']}-{t['end']}|{t['gene']}"
    for b in (getattr(args, "baselines", None) or []):
        if b not in BASELINES_ALL:
            raise SystemExit(f"unknown baseline {b}; choose from {', '.join(BASELINES_ALL)}")
        if tss_annot is None:
            raise SystemExit("--baselines needs --tss-bed (the TSS universe, e.g. CollapsedGeneBounds.hg38.TSS500bp.bed)")
        tss_genes = {a[3] for a in tss_annot}
        truth_b = [t for t in truth if t["gene"] in tss_genes]        # filterExptGeneUniverse
        vals = compute_baseline(b, truth_b, tss_annot, gene_annot, expressed)
        keep = [i for i, v in enumerate(vals) if v == v]
        truth_b = [truth_b[i] for i in keep]
        vals = [vals[i] for i in keep]
        inv = b in BASELINE_INVERSE
        scores = [-v for v in vals] if inv else vals
        labels_b = [t["regulated"] for t in truth_b]
        ap, p70, curve = pr_curve(labels_b, scores)
        lo, hi = bootstrap_auprc(labels_b, scores, args.bootstrap) if args.bootstrap > 0 else (float("nan"), float("nan"))
        n_pos = sum(labels_b)
        pid = f"baseline.{b}"
        results.append({"pred_id": pid, "pred_col": b, "aggregate": "none", "fill_value": float("nan"), "inverse": inv,
                        "n_pairs": len(labels_b), "n_positive": int(n_pos), "baseline_precision": round(n_pos / max(1, len(labels_b)), 4),
                        "n_pairs_dropped_missing_gene": len(truth) - len(truth_b), "pairs_overlapping_prediction": len(truth_b),
                        "frac_overlapping": 1.0, "auprc": round(ap, 4), "auprc_ci95_low": round(lo, 4), "auprc_ci95_high": round(hi, 4),
                        "auprc_trapz": round(curve.get("auprc_trapz", float("nan")), 4),
                        "auprc_crispr_comparison": round(curve.get("auprc_crispr_comparison", float("nan")), 4),
                        "precision_at_70_recall": round(p70, 4), "auprc_negated": round(pr_curve(labels_b, [-x for x in scores])[0], 4),
                        "n_nonmissing_scores": len(scores)})
        curves[pid] = curve
        pair_scores[pid] = {key_of(t): sc for t, sc in zip(truth_b, scores)}
        print(f"baseline {b}: {len(truth_b):,} pairs, AUPRC(CRISPR_comparison) {curve.get('auprc_crispr_comparison', float('nan')):.4f}, "
              f"P@70 {p70:.4f}")

    # ---- AUPRC by distance-to-TSS bin (comparePredictionsToExperiment.Rmd distanceBins)
    bins_rows = []
    if getattr(args, "dist_bins_kb", None) is not None and tss_annot is not None:
        dist = dict(zip((key_of(t) for t in truth), compute_baseline("distToTSS", truth, tss_annot, None)))
        dvals = [v / 1000.0 for v in dist.values() if v == v]
        if len(args.dist_bins_kb) == 0:
            lo_d, hi_d = min(dvals), max(dvals)                        # cut(breaks = 4): 4 equal-width bins
            edges = [lo_d + (hi_d - lo_d) * k / 4 for k in range(5)]
            edges[-1] += 1e-9
        else:
            edges = [float(x) for x in args.dist_bins_kb]
        lab_of = {key_of(t): t["regulated"] for t in truth}
        for a, b in zip(edges[:-1], edges[1:]):
            in_bin = {k for k, v in dist.items() if v == v and a <= v / 1000.0 < b}
            for pid, sc in pair_scores.items():
                ks = [k for k in sc if k in in_bin]
                if not ks:
                    continue
                labs = [lab_of[k] for k in ks]
                ap_b, p70_b, cv_b = pr_curve(labs, [sc[k] for k in ks])
                bins_rows.append({"pred_id": pid, "dist_bin_kb": f"[{a:g},{b:g})", "n_pairs": len(ks), "n_positive": sum(labs),
                                  "auprc": round(ap_b, 4), "auprc_crispr_comparison": round(cv_b.get("auprc_crispr_comparison", float("nan")), 4),
                                  "precision_at_70_recall": round(p70_b, 4)})
        if bins_rows:
            pd.DataFrame(bins_rows).to_csv(out_dir / "auprc_by_distance_bin.tsv", sep="\t", index=False)
            print(f"CSV: {out_dir / 'auprc_by_distance_bin.tsv'}")

    # ---- delta AUPRC between predictor pairs on their shared CRISPR pairs
    delta_rows = []
    for spec in (getattr(args, "delta", None) or []):
        if ":" not in spec and "," not in spec:
            raise SystemExit(f"--delta expects PRED1,PRED2 (got {spec})")
        a_id, b_id = spec.split(",", 1) if "," in spec else spec.split(":", 1)
        if a_id not in pair_scores or b_id not in pair_scores:
            print(f"WARNING: --delta {spec}: unknown predictor (have {', '.join(pair_scores)})")
            continue
        lab_of = {key_of(t): t["regulated"] for t in truth}
        ks = sorted(set(pair_scores[a_id]) & set(pair_scores[b_id]))
        d0, dlo, dhi, pv = bootstrap_delta_auprc([lab_of[k] for k in ks], [pair_scores[a_id][k] for k in ks],
                                                 [pair_scores[b_id][k] for k in ks], max(args.bootstrap, 100))
        delta_rows.append({"pred1": a_id, "pred2": b_id, "n_pairs": len(ks), "delta_auprc": round(d0, 4),
                           "delta_ci95_low": round(dlo, 4), "delta_ci95_high": round(dhi, 4), "pvalue": pv})
        print(f"delta AUPRC {a_id} - {b_id}: {d0:+.4f} [{dlo:+.4f}, {dhi:+.4f}] p={pv:.3g} ({len(ks):,} shared pairs)")
    if delta_rows:
        pd.DataFrame(delta_rows).to_csv(out_dir / "delta_auprc.tsv", sep="\t", index=False)
        print(f"CSV: {out_dir / 'delta_auprc.tsv'}")

    res = pd.DataFrame(results).sort_values("auprc", ascending=False) if results else pd.DataFrame()
    res.to_csv(out_dir / "benchmark_summary.tsv", sep="\t", index=False)
    cc_preds = [dict(p, pred_id=p["pred_id"].replace(":", "_")) for p in preds if (p.get("_table") or p["pred_id"]) in pred_paths]
    cc_paths = {p["pred_id"].replace(":", "_"): pred_paths[p.get("_table") or p["pred_id"]] for p in preds
                if (p.get("_table") or p["pred_id"]) in pred_paths}
    pc, cy = write_crispr_comparison_configs(out_dir, cc_preds, cc_paths,
                                             Path(args.crispr).expanduser().resolve(), args.label or "sce2g")
    curve_rows = [{"pred_id": k, "recall": r, "precision": pr} for k, c in curves.items()
                  for r, pr in zip(c["recall"], c["precision"])]
    pd.DataFrame(curve_rows).to_csv(out_dir / "pr_curves.tsv", sep="\t", index=False)
    plt = _plt()
    fig_path = None
    if plt and curves and not args.no_plots and len(res):
        top_ids = list(res["pred_id"][:8])          # at most eight curves, the palette's limit
        fig, ax = plt.subplots(figsize=(7, 5))
        for k, (pid, c) in enumerate((pid, curves[pid]) for pid in top_ids):
            r = res[res["pred_id"] == pid].iloc[0]
            ax.plot(c["recall"], c["precision"], lw=1.6, color=SERIES[k % len(SERIES)],
                    label=f"{pid} (AUPRC {r['auprc']:.3f})")
        base = sum(labels) / max(1, len(labels))
        ax.axhline(base, color=AXIS, ls="--", lw=0.8)
        ax.axvline(0.7, color=INK2, ls=":", lw=0.8)
        _style(ax, "Precision-recall on CRISPR element-gene pairs", "recall", "precision")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.legend(fontsize=8, frameon=False, loc="upper right")
        fig.patch.set_facecolor(SURFACE)
        fig.tight_layout()
        fig_path = out_dir / "pr_curves.png"
        fig.savefig(fig_path, dpi=140)
        plt.close(fig)
    if len(res) and not args.no_plots and plt and len(res) > 1:
        fig, ax = plt.subplots(figsize=(9, max(2.5, 0.28 * len(res) + 1)))
        names = list(res["pred_id"])[::-1]
        vals = list(res["auprc"])[::-1]
        ax.barh(names, vals, color=BLUE, height=0.72)
        ax.errorbar(vals, range(len(vals)), xerr=[[max(0, v - lo) for v, lo in zip(vals, list(res["auprc_ci95_low"])[::-1])],
                                                    [max(0, hi - v) for v, hi in zip(vals, list(res["auprc_ci95_high"])[::-1])]],
                    fmt="none", ecolor=INK2, elinewidth=0.8, capsize=2)
        ax.axvline(sum(labels) / max(1, len(labels)), color=AXIS, ls="--", lw=0.8)
        _style(ax, "AUPRC per predictor (bootstrap 95% interval; dashed = random baseline)", "AUPRC")
        fig.patch.set_facecolor(SURFACE)
        fig.tight_layout()
        fig.savefig(out_dir / "auprc_by_predictor.png", dpi=140)
        plt.close(fig)
        print(f"Figure: {out_dir / 'auprc_by_predictor.png'}")
    if len(res):
        cols = ["pred_id", "n_pairs", "n_positive", "frac_overlapping", "auprc", "auprc_ci95_low",
                "auprc_ci95_high", "precision_at_70_recall", "auprc_negated"]
        print(res[cols].head(40).to_string(index=False))
        if len(res) > 40:
            print(f"... {len(res) - 40} more rows in benchmark_summary.tsv")
    print(f"CSV: {out_dir / 'benchmark_summary.tsv'}")
    print(f"CSV: {out_dir / 'pr_curves.tsv'}")
    if fig_path:
        print(f"Figure: {fig_path}")
    print(f"CRISPR_comparison pred_config: {pc}")
    print(f"CRISPR_comparison config: {cy}")
    print(f"Log: {log_path}")
    return 0 if results else 1


# ---------------------------------------------------------------------------
# merge: several benchmark runs -> one ranked table and one figure
# ---------------------------------------------------------------------------

def cmd_merge(args: argparse.Namespace) -> int:
    pd = _pd()
    log_path = setup_logging()
    out_dir = run_dir(args.label or "benchmark_merged")
    frames = []
    for r in args.runs:
        rp = Path(r).expanduser()
        f = rp / "benchmark_summary.tsv" if rp.is_dir() else rp
        if not f.is_file():
            print(f"WARNING: no benchmark_summary.tsv in {r}; skipped")
            continue
        df = pd.read_csv(f, sep="\t")
        df["run"] = rp.name if rp.is_dir() else rp.parent.name
        frames.append(df)
    if not frames:
        print("ERROR: nothing to merge", file=sys.stderr)
        return 1
    res = pd.concat(frames, ignore_index=True)
    # Shared-universe columns (E2G_Distance) appear once per table in runs made
    # before they were collapsed; they are one baseline predictor.
    col = res["pred_id"].str.split(":").str[-1]
    res.loc[col.isin(SHARED_UNIVERSE_COLS), "pred_id"] = "baseline:" + col[col.isin(SHARED_UNIVERSE_COLS)]
    res = res.drop_duplicates("pred_id").sort_values("auprc", ascending=False)
    baseline = float(res["n_positive"].iloc[0]) / float(res["n_pairs"].iloc[0]) if len(res) else float("nan")
    res["auprc_over_random"] = (res["auprc"] / baseline).round(2)
    res["best_orientation_auprc"] = res[["auprc", "auprc_negated"]].max(axis=1)
    res.to_csv(out_dir / "benchmark_summary.tsv", sep="\t", index=False)
    cols = ["pred_id", "auprc", "auprc_crispr_comparison", "auprc_ci95_low", "auprc_ci95_high", "precision_at_70_recall",
            "auprc_negated", "auprc_over_random", "frac_overlapping"]
    show = res[[c for c in cols if c in res.columns]]
    print(show.head(args.top).to_string(index=False))
    if len(res) > args.top:
        print(f"... {len(res) - args.top} more in benchmark_summary.tsv")
    print(f"random baseline precision {baseline:.4f} ({int(res['n_positive'].iloc[0])} / {int(res['n_pairs'].iloc[0])} pairs)")
    plt = _plt()
    if plt and not args.no_plots:
        top = res.head(args.top)
        fig, ax = plt.subplots(figsize=(9, max(3, 0.3 * len(top) + 1)))
        names, vals = list(top["pred_id"])[::-1], list(top["auprc"])[::-1]
        colors = [ORANGE if str(n).startswith("baseline:") else BLUE for n in names]
        ax.barh(names, vals, color=colors, height=0.72)
        lo = [max(0, v - l) for v, l in zip(vals, list(top["auprc_ci95_low"])[::-1])]
        hi = [max(0, h - v) for v, h in zip(vals, list(top["auprc_ci95_high"])[::-1])]
        ax.errorbar(vals, range(len(vals)), xerr=[lo, hi], fmt="none", ecolor=INK2, elinewidth=0.8, capsize=2)
        ax.axvline(baseline, color=AXIS, ls="--", lw=0.8)
        for i, v in enumerate(vals):
            ax.text(v, i, f" {v:.3f}", va="center", fontsize=8, color=INK)
        _style(ax, f"AUPRC on K562 CRISPR pairs, top {len(top)} predictors (dashed = random; orange = baseline)", "AUPRC")
        fig.patch.set_facecolor(SURFACE)
        fig.tight_layout()
        fig.savefig(out_dir / "auprc_ranked.png", dpi=140)
        plt.close(fig)
        print(f"Figure: {out_dir / 'auprc_ranked.png'}")
    # report.md: the merged ranking with the method stated once
    L = [f"# E2G predictor benchmark on CRISPR element-gene pairs: {args.label or 'merged'}", "",
         f"Generated {time.strftime('%Y-%m-%d %H:%M')} by `igvfagent sce2g merge` from {len(frames)} benchmark run(s): "
         + ", ".join(sorted(set(res['run']))) + ".", "",
         f"CRISPR pairs: {int(res['n_pairs'].iloc[0]):,}, of which {int(res['n_positive'].iloc[0]):,} regulated "
         f"(random-baseline precision {baseline:.4f}). Each predictor's score for a CRISPR pair is the aggregate of the "
         "predicted elements overlapping the tested element for the same gene (max unless stated), `fill_value` when "
         "none overlaps; p-value / FDR / distance style columns are inverted so that higher means more likely. "
         "`auprc_negated` is the AUPRC of the negated score, so a feature that works in the other direction is visible.", "",
         "## Ranked predictors", "",
         md_table(["rank", "predictor", "AUPRC", "95% CI", "precision @ 70% recall", "AUPRC negated", "x random", "pairs overlapped"],
                  [[i + 1, r.pred_id, f"{r.auprc:.4f}", f"{r.auprc_ci95_low:.3f}-{r.auprc_ci95_high:.3f}",
                    f"{r.precision_at_70_recall:.3f}", f"{r.auprc_negated:.4f}", f"{r.auprc_over_random:.1f}",
                    f"{r.frac_overlapping:.3f}"] for i, r in enumerate(res.itertuples())], max_rows=400), "",
         "## Reading the table", "",
         "- A predictor at the random baseline carries no information about which tested elements regulate their gene.",
         "- `pairs overlapped` below 1.0 means some CRISPR elements had no predicted element for that gene in the table; "
         "those pairs scored `fill_value` and count against the predictor, exactly as in CRISPR_comparison.",
         "- Bootstrap intervals resample CRISPR pairs with replacement; overlapping intervals mean the ranking between two "
         "predictors is not settled by this dataset.", "",
         "## Files", "", "- `benchmark_summary.tsv`", "- `auprc_ranked.png`" if plt and not args.no_plots else "", ""]
    (out_dir / "report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    # summary.json: what Benchmarks/concordance.py scores. Keys are the pred_ids
    # with '.' replaced, because the scorer walks dotted paths.
    def _key(pid: str) -> str:
        return str(pid).replace(".", "_")
    feats = res[~res["pred_id"].astype(str).str.startswith("baseline:") & res["pred_id"].astype(str).str.contains(":")]
    refs = res[~res["pred_id"].astype(str).str.contains(":")]
    summary = {
        "label": args.label, "runs": sorted(set(res["run"])),
        "n_pairs": int(res["n_pairs"].iloc[0]), "n_positive": int(res["n_positive"].iloc[0]),
        "baseline_precision": round(baseline, 5), "n_predictors": int(len(res)),
        "n_single_features": int(len(feats)),
        "best_single_feature": ({"pred_id": str(feats.iloc[0]["pred_id"]), "auprc": float(feats.iloc[0]["auprc"])}
                                if len(feats) else None),
        "best_reference": ({"pred_id": str(refs.iloc[0]["pred_id"]), "auprc": float(refs.iloc[0]["auprc"])}
                           if len(refs) else None),
        "n_features_above_distance_baseline": int((feats["auprc"] > float(
            res.loc[res["pred_id"] == "baseline:E2G_Distance", "auprc"].iloc[0]) if (res["pred_id"] == "baseline:E2G_Distance").any() else 0).sum()),
        "predictors": {_key(r.pred_id): {"auprc": float(r.auprc), "auprc_ci95_low": float(r.auprc_ci95_low),
                                          "auprc_ci95_high": float(r.auprc_ci95_high),
                                          "auprc_trapz": float(getattr(r, "auprc_trapz", float("nan"))),
                                          "auprc_crispr_comparison": float(getattr(r, "auprc_crispr_comparison", float("nan"))),
                                          "precision_at_70_recall": float(r.precision_at_70_recall),
                                          "auprc_negated": float(r.auprc_negated),
                                          "frac_overlapping": float(r.frac_overlapping)}
                       for r in res.itertuples()},
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"JSON: {out_dir / 'summary.json'}")
    print(f"Report: {out_dir / 'report.md'}")
    print(f"CSV: {out_dir / 'benchmark_summary.tsv'}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# inventory: what is in a folder of crowdsourced feature tables?
# ---------------------------------------------------------------------------

def cmd_inventory(args: argparse.Namespace) -> int:
    pd = _pd()
    log_path = setup_logging()
    out_dir = run_dir(args.label or "feature_inventory")
    files = sorted(p for p in Path(args.dir).expanduser().glob("*") if p.suffix in (".gz", ".tsv", ".txt"))
    rows, feat_rows = [], []
    ref_pairs = None
    for f in files:
        n = 0
        nan: "dict[str, int]" = {}
        feats: "list[str]" = []
        sample_pairs = []
        header = None
        try:
            for chunk in pd.read_csv(f, sep="\t", comment="#", low_memory=False, chunksize=1_000_000):
                if header is None:
                    header = list(chunk.columns)
                    feats = feature_columns(chunk)
                    for c in feats:
                        nan[c] = 0
                for c in feats:
                    nan[c] += int(chunk[c].isna().sum())
                if {"ElementName", "GeneSymbol"} <= set(chunk.columns):
                    sub = chunk.iloc[::1000]
                    sample_pairs.extend((sub["ElementName"].astype(str) + "|" + sub["GeneSymbol"].astype(str)).tolist())
                n += len(chunk)
                if args.max_rows and n >= args.max_rows:
                    break
        except Exception as exc:
            rows.append({"file": f.name, "rows": n, "error": str(exc)[:120]})
            continue
        pairs = set(sample_pairs)
        if ref_pairs is None:
            ref_pairs = pairs
            same = "reference"
        else:
            same = f"{len(pairs & ref_pairs) / max(1, len(ref_pairs)):.3f}"
        rows.append({"file": f.name, "rows": n, "columns": len(header or []), "features": len(feats),
                     "feature_names": "; ".join(feats)[:300], "sampled_pair_overlap_vs_first": same,
                     "size_mb": round(f.stat().st_size / 1e6, 1)})
        for c in feats:
            feat_rows.append({"file": f.name, "feature": c, "clean_name": sanitize_feature(c),
                              "missing_fraction": round(nan[c] / max(1, n), 4)})
    inv = pd.DataFrame(rows)
    fr = pd.DataFrame(feat_rows)
    inv.to_csv(out_dir / "feature_tables.tsv", sep="\t", index=False)
    fr.to_csv(out_dir / "features.tsv", sep="\t", index=False)
    print(inv.drop(columns=[c for c in ("feature_names",) if c in inv.columns]).to_string(index=False))
    print(f"{len(fr)} feature column(s) in {len(inv)} table(s)")
    print(f"CSV: {out_dir / 'feature_tables.tsv'}")
    print(f"CSV: {out_dir / 'features.tsv'}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------

def _fake_repo(root: Path) -> Path:
    repo = root / "scE2G"
    (repo / "workflow").mkdir(parents=True)
    (repo / "ENCODE-rE2G" / "workflow").mkdir(parents=True)
    (repo / "config").mkdir()
    (repo / "models" / "multiome_powerlaw_v3").mkdir(parents=True)
    (repo / "resources" / "genome_annotations").mkdir(parents=True)
    (repo / "workflow" / "Snakefile_training").write_text(
        'import os\nconfigfile: "config/config_training.yaml"\nconda: "mamba"\n\n'
        "WORKFLOW_DIR = os.path.dirname(str(workflow.current_basedir))\n"
        'include: os.path.join("rules", "utils.smk")\n')
    (repo / "ENCODE-rE2G" / "workflow" / "Snakefile_training").write_text(
        'configfile: "config/config_training.yaml"\nconda: "mamba"\nrule all:\n    input: []\n')
    (repo / "config" / "config_training.yaml").write_text(
        'cell_clusters: "config/config_cell_clusters.tsv"\nmodel_config: "config/config_models.tsv"\n'
        'results_dir: "results/"\ngene_annotations: "resources/genome_annotations/genes.gtf.gz"\n'
        'genes: "resources/genome_annotations/genes.bed"\n')
    for f in ("genes.gtf.gz", "genes.bed"):
        (repo / "resources" / "genome_annotations" / f).write_bytes(b"")
    write_tsv_rows(repo / "config" / "config_cell_clusters.tsv", CLUSTER_COLS, [
        {"cluster": "K562_chr22_cluster1", "rna_matrix_file": "resources/example/rna.csv.gz",
         "atac_frag_file": "resources/example/atac.tsv.gz", "model_dir": "models/multiome_powerlaw_v3"}])
    write_tsv_rows(repo / "config" / "config_models.tsv", MODEL_COLS, [
        {"model": "sample_model_multiome", "dataset": "K562_chr22_cluster1", "ABC_directory": "",
         "feature_table": "resources/feature_tables/multiome_arc_n6.tsv", "polynomial": "False", "override_params": ""}])
    (repo / "resources" / "example").mkdir()
    for f in ("rna.csv.gz", "atac.tsv.gz"):
        (repo / "resources" / "example" / f).write_bytes(b"")
    return repo


def _synthetic_inputs(root: Path) -> dict:
    rng = random.Random(3)
    genes = [f"G{i}" for i in range(60)]
    feat_rows, pred_a, pred_b, crispr = [], [], [], []
    for gi, g in enumerate(genes):
        gene_tss = 1_000_000 + gi * 100_000
        for k in range(6):
            st = gene_tss + (k - 3) * 8_000
            en = st + 500
            truth_reg = k == 1                        # element 1 of each gene truly regulates it
            feat_rows.append({"ElementChr": "chr1", "ElementStart": st, "ElementEnd": en, "ElementName": f"e{gi}_{k}",
                              "GeneSymbol": g, "H3K27ac mean signal": rng.uniform(0, 5) + (3 if truth_reg else 0),
                              "CTCF binding": rng.uniform(0, 1), "notes": "text"})
            good = (0.7 + rng.uniform(0, 0.3)) if truth_reg else rng.uniform(0, 0.4)
            pred_a.append({"ElementChr": "chr1", "ElementStart": st, "ElementEnd": en, "GeneSymbol": g,
                           "E2G.Score.qnorm": round(good, 4), "ElementClass": "intergenic" if k != 3 else "promoter",
                           "isSelfPromoter": k == 3})
            pred_b.append({"chr": "chr1", "start": st, "end": en, "TargetGene": g, "ABC.Score": round(rng.uniform(0, 1), 4)})
            if k in (0, 1, 2):
                crispr.append({"chrom": "chr1", "chromStart": st + 50, "chromEnd": en - 50, "measuredGeneSymbol": g,
                               "EffectSize": -0.5 if truth_reg else 0.01, "Significant": str(truth_reg).upper(),
                               "Regulated": str(truth_reg).upper(), "CellType": "K562"})
    paths = {}
    for name, rows in (("features.tsv", feat_rows), ("pred_sce2g.e2g.tsv", pred_a), ("pred_abc.tsv", pred_b),
                       ("crispr.tsv", crispr)):
        p = root / name
        with open(p, "w", newline="") as fh:
            if name.endswith(".e2g.tsv"):
                fh.write("# synthetic scE2G output\n")
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        paths[name] = p
    return paths


def cmd_selftest(args: argparse.Namespace) -> int:
    pd = _pd()
    checks: "list[tuple[bool, str]]" = []

    def check(cond: bool, msg: str) -> None:
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        repo = _fake_repo(root)
        inputs = _synthetic_inputs(root)
        # setup --check reports, changes nothing
        rc = cmd_setup(argparse.Namespace(repo_dir=str(repo), branch=SCE2G_DEFAULT_BRANCH, check=True,
                                          keep_branch=True, offline=True))
        check(rc == 0 and 'conda: "mamba"' in (repo / "workflow" / "Snakefile_training").read_text(),
              "setup --check leaves the checkout untouched")
        rc = cmd_setup(argparse.Namespace(repo_dir=str(repo), branch=SCE2G_DEFAULT_BRANCH, check=False,
                                          keep_branch=True, offline=True))
        sf = (repo / "workflow" / "Snakefile_training").read_text()
        check(rc == 0, "setup applies without problems")
        check('conda: "mamba"' not in sf and 'conda: "mamba"' not in (repo / "ENCODE-rE2G" / "workflow" / "Snakefile_training").read_text(),
              "both `conda: \"mamba\"` lines removed")
        check(sf.splitlines()[sf.splitlines().index(next(l for l in sf.splitlines() if l.startswith("WORKFLOW_DIR"))) + 1] == SCRIPTS_DIR_LINE,
              "SCRIPTS_DIR inserted right after WORKFLOW_DIR")
        y = (repo / "config" / "config_training.yaml").read_text()
        check("RNA_matrix_filtered: True" in y and "max_cell_count: 20000" in y, "training yaml gained the two keys")
        check((repo / "resources" / "feature_tables" / "multiome_arc_n6.tsv").read_text() == MULTIOME_ARC_N6,
              "missing multiome_arc_n6.tsv written (offline copy)")
        rc2 = cmd_setup(argparse.Namespace(repo_dir=str(repo), branch=SCE2G_DEFAULT_BRANCH, check=False,
                                           keep_branch=True, offline=True))
        check(rc2 == 0 and y == (repo / "config" / "config_training.yaml").read_text(), "setup is idempotent")

        # features
        rc = cmd_features(argparse.Namespace(repo_dir=str(repo), out_dir=None, name="crowd", features=str(inputs["features.tsv"]),
                                             select=None, merge_aggregate="mean", benchmark_aggregate="max", fill_value="0"))
        check(rc == 0, "features conversion runs")
        src = repo / "resources" / "external_features" / "crowd.tsv.gz"
        with gzip.open(src, "rt") as fh:
            head = fh.readline().rstrip("\n").split("\t")
        check(head == ["chr", "start", "end", "TargetGene", "H3K27ac_mean_signal", "CTCF_binding"],
              f"source_file header renamed and de-spaced ({head})")
        eh, erows = read_tsv_rows(repo / "config" / "external_features_config_crowd.tsv")
        check(eh == EXT_CONFIG_COLS and [r["input_col"] for r in erows] == ["H3K27ac_mean_signal", "CTCF_binding"]
              and all(r["aggregate_function"] == "mean" and r["join_by"] == "overlap" for r in erows)
              and erows[0]["source_file"] == "resources/external_features/crowd.tsv.gz",
              "external_features_config has the five columns, mean/overlap, relative source_file")
        fh_, frows = read_tsv_rows(repo / "resources" / "feature_tables" / "multiome_arc_n6_crowd.tsv")
        new = [r for r in frows if r["feature"] in ("H3K27ac_mean_signal", "CTCF_binding")]
        check(fh_ == FEATURE_TABLE_COLS and len(frows) == 8 and len(new) == 2
              and new[0]["nice_name"] == "H3K27ac mean signal" and new[0]["aggregate_function"] == "max"
              and new[0]["fill_value"] == "0" and new[0]["second_input"] == "NA",
              "feature table = 6 base rows + 2 new rows with max / 0 / NA / nice_name with spaces")

        # configure + check
        rc = cmd_configure(argparse.Namespace(repo_dir=str(repo), cluster="K562_crowd", rna="resources/example/rna.csv.gz",
                                              atac_frag="resources/example/atac.tsv.gz", hic=None, hic_type=None,
                                              hic_resolution=None, model_dir="models/multiome_powerlaw_v3", name="crowd",
                                              model="multiome_crowd", feature_table=None, profile="profiles/slurm", jobs=4))
        ch, crows = read_tsv_rows(repo / "config" / "config_cell_clusters.tsv")
        mh, mrows = read_tsv_rows(repo / "config" / "config_models.tsv")
        check(rc == 0 and ch == CLUSTER_COLS + [EXT_COL] and len(crows) == 2
              and crows[1]["external_features_config"] == "config/external_features_config_crowd.tsv",
              "cluster row appended with the external_features_config column")
        check(mh == MODEL_COLS and mrows[1] == {"model": "multiome_crowd", "dataset": "K562_crowd", "ABC_directory": "",
                                                 "feature_table": "resources/feature_tables/multiome_arc_n6_crowd.tsv",
                                                 "polynomial": "False", "override_params": ""},
              "model row matches the walkthrough")
        script = (repo / "run_training_multiome_crowd.sh").read_text()
        check("--profile profiles/slurm" in script and "add_external_features:mem_mb=256000" in script,
              "run script carries the Slurm profile and resource overrides")
        problems, notes = check_setup(repo, "multiome_crowd")
        check(not problems, f"check passes on a consistent setup ({problems})")
        # break it: a feature in the external config that the feature table does not list
        erows.append({"input_col": "Orphan", "source_col": "Orphan", "aggregate_function": "mean", "join_by": "overlap",
                      "source_file": "resources/external_features/crowd.tsv.gz"})
        write_tsv_rows(repo / "config" / "external_features_config_crowd.tsv", EXT_CONFIG_COLS, erows)
        problems, _ = check_setup(repo, "multiome_crowd")
        check(any("Orphan" in p and "feature row" in p for p in problems) and any("column 'Orphan' missing" in p for p in problems),
              "check catches an external feature missing from the feature table and the source file")
        # re-configure replaces, not duplicates
        cmd_configure(argparse.Namespace(repo_dir=str(repo), cluster="K562_crowd", rna="resources/example/rna.csv.gz",
                                         atac_frag="resources/example/atac.tsv.gz", hic=None, hic_type=None, hic_resolution=None,
                                         model_dir="models/multiome_powerlaw_v3", name="crowd", model="multiome_crowd",
                                         feature_table=None, profile=None, jobs=2))
        _, crows2 = read_tsv_rows(repo / "config" / "config_cell_clusters.tsv")
        check(len(crows2) == 2, "re-running configure replaces the row instead of duplicating it")

        # predictions + benchmark
        rc = cmd_predictions(argparse.Namespace(predictions=[f"scE2G={inputs['pred_sce2g.e2g.tsv']}", f"ABC={inputs['pred_abc.tsv']}"],
                                                score_col="E2G.Score.qnorm", threshold=0.177, label="selftest_pred", no_plots=args.no_plots))
        out_p = sorted(OUT_ROOT.glob("*_selftest_pred"))[-1]
        summ = pd.read_csv(out_p / "predictions_summary.tsv", sep="\t")
        pair = pd.read_csv(out_p / "predictions_pairwise.tsv", sep="\t")
        check(rc == 0 and list(summ["links"]) == [360, 360] and int(pair["shared_pairs"][0]) == 360 and float(pair["jaccard"][0]) == 1.0,
              "predictions summary and pairwise overlap (same pairs, different score columns)")
        pc_path = root / "pred_config.txt"
        write_tsv_rows(pc_path, PRED_CONFIG_COLS + ["include"], [
            {"pred_id": "scE2G", "pred_col": "E2G.Score.qnorm", "boolean": "FALSE", "alpha": 0.177, "aggregate_function": "max",
             "fill_value": 0, "inverse_predictor": "FALSE", "pred_name_long": "scE2G multiome", "color": "blue", "include": "TRUE"},
            {"pred_id": "ABC", "pred_col": "ABC.Score", "boolean": "FALSE", "alpha": 0.02, "aggregate_function": "max",
             "fill_value": 0, "inverse_predictor": "FALSE", "pred_name_long": "ABC", "color": "grey", "include": "TRUE"}])
        rc = cmd_benchmark(argparse.Namespace(predictions=[f"scE2G={inputs['pred_sce2g.e2g.tsv']}", f"ABC={inputs['pred_abc.tsv']}"],
                                              crispr=str(inputs["crispr.tsv"]), pred_config=str(pc_path), score_col="E2G.Score.qnorm",
                                              cell_type=None, bootstrap=50, label="selftest_bench", no_plots=args.no_plots,
                                              all_features=False, keep_scored_pairs=False,
                                              gene_universe_filter=False))
        out_b = sorted(OUT_ROOT.glob("*_selftest_bench"))[-1]
        res = pd.read_csv(out_b / "benchmark_summary.tsv", sep="\t").set_index("pred_id")
        check(rc == 0 and res.loc["scE2G", "auprc"] > 0.9 and res.loc["ABC", "auprc"] < 0.6
              and res.loc["scE2G", "auprc_ci95_low"] <= res.loc["scE2G", "auprc"] <= res.loc["scE2G", "auprc_ci95_high"],
              f"planted predictor wins the benchmark (scE2G {res.loc['scE2G', 'auprc']:.3f} vs ABC {res.loc['ABC', 'auprc']:.3f}), CI brackets it")
        check(int(res.loc["scE2G", "n_pairs"]) == 180 and int(res.loc["scE2G", "n_positive"]) == 60
              and float(res.loc["scE2G", "frac_overlapping"]) == 1.0, "all 180 CRISPR pairs scored, 60 positives, all overlapped")
        pch, pcrows = read_tsv_rows(out_b / "pred_config.txt")
        check(pch == PRED_CONFIG_COLS and [r["pred_id"] for r in pcrows] == ["scE2G", "ABC"] and "comparisons:" in (out_b / "config.yml").read_text(),
              "CRISPR_comparison pred_config.txt and config.yml written")
        # inverse predictor semantics: negating the good score should flip the ranking
        df = pd.read_csv(inputs["pred_sce2g.e2g.tsv"], sep="\t", comment="#")
        from e2g_benchmark_eval import read_ground_truth
        _manual_pr = pr_curve
        truth = read_ground_truth(str(inputs["crispr.tsv"]))
        lab = [t["regulated"] for t in truth]
        s_norm, _ = score_crispr_pairs(truth, df, "E2G.Score.qnorm", "max", 0.0, False, False)
        s_inv, _ = score_crispr_pairs(truth, df, "E2G.Score.qnorm", "max", 0.0, True, False)
        check(_manual_pr(lab, s_norm)[0] > 0.9 and _manual_pr(lab, s_inv)[0] < 0.5, "inverse_predictor flips the ranking")
        s_mean, _ = score_crispr_pairs(truth, df, "E2G.Score.qnorm", "mean", 0.0, False, False)
        check(s_mean == s_norm, "mean == max when exactly one element overlaps each pair")
        # tie-aware PR: two orderings of the same tied data must give one AUPRC, and it must match the
        # closed form for a single tie block (all pairs at one score -> precision = prevalence).
        ap_a = pr_curve([True, False, True, False], [1, 1, 1, 1])[0]
        ap_b = pr_curve([False, True, False, True], [1, 1, 1, 1])[0]
        check(abs(ap_a - 0.5) < 1e-12 and abs(ap_b - 0.5) < 1e-12, "tie-aware AUPRC is order-independent (single tie block = prevalence)")
        ap_untied = pr_curve([True, True, False, False], [4, 3, 2, 1])[0]
        check(abs(ap_untied - 1.0) < 1e-12, "perfect ranking gives AUPRC 1")
        c4 = pr_curve([True, False, True, False], [4, 3, 2, 1])[2]
        # points: (r=.5,p=1) (.5,.5) (1,.667) (1,.5): trapz = 0 + .5*(.5+.667)/2 + 0 = 0.2917
        check(abs(c4["auprc_trapz"] - 0.2916667) < 1e-6, f"trapezoid AUPRC matches hand calculation ({c4['auprc_trapz']:.4f})")
        # CRISPR_comparison baselines on a toy locus: TSS of G1 at 1000, G2 at 5000; gene bodies 800-2000 / 4800-9000
        tss_t = [("chr1", 999, 1000, "G1"), ("chr1", 4999, 5000, "G2")]
        gene_t = [("chr1", 800, 2000, "G1"), ("chr1", 4800, 9000, "G2")]
        truth_t = [{"chrom": "chr1", "start": 1100, "end": 1300, "gene": "G1", "regulated": True},
                   {"chrom": "chr1", "start": 1100, "end": 1300, "gene": "G2", "regulated": False},
                   {"chrom": "chr1", "start": 300000, "end": 300200, "gene": "G1", "regulated": False}]
        d_tss = compute_baseline("distToTSS", truth_t, tss_t, gene_t)
        check(d_tss[:2] == [200.0, 3798.0], f"distToTSS: centre-to-centre GenomicRanges distance ({d_tss[:2]})")
        check(compute_baseline("distToGene", truth_t, tss_t, gene_t)[:2] == [0.0, 3599.0], "distToGene: 0 inside the gene body")
        check(compute_baseline("nearestTSS", truth_t, tss_t, gene_t) == [1.0, 0.0, 0.0], "nearestTSS: only the closest gene scores 1")
        check(compute_baseline("within100kbTSS", truth_t, tss_t, gene_t) == [1.0, 1.0, 0.0], "within100kbTSS: +/-100 kb window")
        pd_ = _pd()
        pr_t = pd_.DataFrame({"chr": ["chr1", "chr1"], "start": [900, 3000], "end": [1100, 3200], "TargetGene": ["G1", "G1"]})
        check(len(filter_predictions_tss(pr_t, tss_t)) == 1, "filter-pred-tss drops the element overlapping a TSS")
        d0, dlo, dhi, pv = bootstrap_delta_auprc([True, False, True, False, True, False] * 5, list(range(30, 0, -1)),
                                                 [((i * 7) % 11) for i in range(30)], 200)
        check(d0 > 0 and dlo <= d0 <= dhi and 0 < pv <= 1, f"delta AUPRC bootstrap: {d0:+.3f} [{dlo:+.3f}, {dhi:+.3f}] p={pv:.3g}")
        if not args.keep:
            shutil.rmtree(out_p, ignore_errors=True)
            shutil.rmtree(out_b, ignore_errors=True)
    ok = all(c for c, _ in checks)
    print("selftest: all checks pass" if ok else f"selftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[list[str]]" = None) -> int:
    p = argparse.ArgumentParser(prog="igvfagent sce2g",
                                description="scE2G workbench: set up, configure, check, run and benchmark scE2G "
                                            "model training with crowdsourced features.")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("setup", help="Clone/patch an scE2G checkout the way the walkthrough says.")
    s.add_argument("--repo-dir", required=True, help="Where the scE2G checkout is / should go.")
    s.add_argument("--branch", default=SCE2G_DEFAULT_BRANCH)
    s.add_argument("--check", action="store_true", help="Report only; change nothing.")
    s.add_argument("--keep-branch", action="store_true", help="Do not switch an existing checkout's branch.")
    s.add_argument("--offline", action="store_true", help="Never fetch; use the embedded feature table.")
    s.set_defaults(func=cmd_setup)

    s = sub.add_parser("features", help="Synapse feature table -> source_file + external config + feature table.")
    s.add_argument("--features", required=True, help="Crowdsourced E2G feature table (TSV; ElementChr/.../GeneSymbol).")
    s.add_argument("--name", required=True, help="Feature-set name, e.g. h3k27ac (used in every file name).")
    s.add_argument("--repo-dir", help="scE2G checkout; files land under its config/ and resources/.")
    s.add_argument("--out-dir", help="Alternative output root when no --repo-dir.")
    s.add_argument("--select", help="Comma list of feature columns to keep (default: every numeric column).")
    s.add_argument("--merge-aggregate", default="mean", help="external_features_config aggregate_function.")
    s.add_argument("--benchmark-aggregate", default="max", help="feature_table aggregate_function for CRISPR benchmark.")
    s.add_argument("--fill-value", default="0")
    s.set_defaults(func=cmd_features)

    s = sub.add_parser("configure", help="Write the cluster + model rows and the run script.")
    s.add_argument("--repo-dir", required=True)
    s.add_argument("--cluster", required=True, help="Arbitrary cluster name (also the model's dataset).")
    s.add_argument("--rna", help="rna_matrix_file (.csv.gz); omit for ATAC-only.")
    s.add_argument("--atac-frag", required=True, help="atac_frag_file (.tsv.gz).")
    s.add_argument("--hic", help="HiC_file (optional).")
    s.add_argument("--hic-type", help="HiC_type (optional).")
    s.add_argument("--hic-resolution", help="HiC_resolution (optional).")
    s.add_argument("--model-dir", default="models/multiome_powerlaw_v3")
    s.add_argument("--name", help="Feature-set name from `features` (adds external_features_config).")
    s.add_argument("--model", help="Model name (default multiome_<name>).")
    s.add_argument("--feature-table", help="Override feature_table path.")
    s.add_argument("--profile", help="Snakemake profile dir, e.g. profiles/slurm (default: local -j).")
    s.add_argument("--jobs", type=int, default=4)
    s.set_defaults(func=cmd_configure)

    s = sub.add_parser("check", help="Validate paths and names before spending compute.")
    s.add_argument("--repo-dir", required=True)
    s.add_argument("--model", help="Check one model (default all).")
    s.add_argument("--profile")
    s.add_argument("--jobs", type=int, default=4)
    s.set_defaults(func=cmd_check)

    s = sub.add_parser("run", help="check, then snakemake (dry-run unless --execute).")
    s.add_argument("--repo-dir", required=True)
    s.add_argument("--model")
    s.add_argument("--profile")
    s.add_argument("--jobs", type=int, default=4)
    s.add_argument("--execute", dest="dry_run", action="store_false", help="Really run (default dry-run).")
    s.add_argument("--force", action="store_true", help="Run even if check fails.")
    s.add_argument("extra", nargs="*", help="Extra snakemake arguments.")
    s.set_defaults(func=cmd_run, dry_run=True)

    s = sub.add_parser("predictions", help="Describe / compare scE2G prediction tables.")
    s.add_argument("--predictions", action="append", required=True, metavar="LABEL=PATH")
    s.add_argument("--score-col", default="E2G.Score.qnorm")
    s.add_argument("--threshold", type=float, default=0.177)
    s.add_argument("--label", default="")
    s.add_argument("--no-plots", action="store_true")
    s.set_defaults(func=cmd_predictions)

    s = sub.add_parser("benchmark", help="CRISPR_comparison-style evaluation against CRISPR element-gene pairs.")
    s.add_argument("--predictions", action="append", required=True, metavar="LABEL=PATH")
    s.add_argument("--crispr", required=True, help="EPCrisprBenchmark TSV (chrom, chromStart, chromEnd, measuredGeneSymbol, Regulated).")
    s.add_argument("--pred-config", help="CRISPR_comparison pred_config.txt (else defaults: max / fill 0).")
    s.add_argument("--score-col", default="E2G.Score.qnorm", help="Default pred_col when no pred_config.")
    s.add_argument("--cell-type", help="Restrict CRISPR pairs to this CellType.")
    s.add_argument("--bootstrap", type=int, default=200, help="Bootstrap resamples for the AUPRC interval (0 = off).")
    s.add_argument("--gene-universe-filter", action="store_true",
                   help="CRISPR_comparison semantics: drop tested pairs whose gene has no row in the prediction "
                        "table instead of scoring them as fill_value.")
    s.add_argument("--all-features", action="store_true",
                   help="Benchmark EVERY feature column of each table as its own predictor (crowdsourced-feature benchmark).")
    s.add_argument("--keep-scored-pairs", action="store_true", help="With --all-features, still write per-predictor scored pairs.")
    s.add_argument("--tss-bed", help="TSS universe BED (CRISPR_comparison tss_universe), needed for --baselines / --dist-bins-kb.")
    s.add_argument("--gene-bed", help="Gene-body BED (gene_universe), needed for distToGene / nearestGene / within100kbGene.")
    s.add_argument("--expressed-genes", help="TSV cell_type, gene, expressed (for the *Expr* baselines).")
    s.add_argument("--baselines", nargs="+", action="extend", choices=list(BASELINES_ALL),
                   help="CRISPR_comparison baseline predictors to add (e.g. distToTSS nearestTSS within100kbTSS).")
    s.add_argument("--filter-pred-tss", help="TSS BED: drop predicted elements overlapping a gene TSS (upstream filter_pred_tss: True).")
    s.add_argument("--dist-bins-kb", nargs="*", action="extend", type=float,
                   help="AUPRC by distance-to-TSS bin; give bin edges in kb (e.g. 0 20 100 2500) or no values for 4 equal-width bins.")
    s.add_argument("--delta", action="append", metavar="PRED1,PRED2",
                   help="Bootstrap delta AUPRC (upstream definition) between two predictors on shared pairs (repeatable).")
    s.add_argument("--label", default="")
    s.add_argument("--no-plots", action="store_true")
    s.set_defaults(func=cmd_benchmark)

    s = sub.add_parser("merge", help="Merge several benchmark runs into one ranked table and figure.")
    s.add_argument("--runs", action="append", required=True, help="Benchmark run directory (repeatable).")
    s.add_argument("--top", type=int, default=30)
    s.add_argument("--label", default="")
    s.add_argument("--no-plots", action="store_true")
    s.set_defaults(func=cmd_merge)

    s = sub.add_parser("inventory", help="Describe a folder of crowdsourced feature tables (rows, features, missingness, universe agreement).")
    s.add_argument("--dir", required=True)
    s.add_argument("--max-rows", type=int, default=0, help="Stop after this many rows per file (0 = all).")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_inventory)

    s = sub.add_parser("selftest", help="Fake checkout + synthetic inputs; asserts every artefact.")
    s.add_argument("--no-plots", action="store_true")
    s.add_argument("--keep", action="store_true")
    s.set_defaults(func=cmd_selftest)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
