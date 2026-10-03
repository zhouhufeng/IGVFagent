#!/usr/bin/env python3
"""Bulk RNA-seq analysis skill.

Picks up where the GEO retrieval skill (or any local pipeline) leaves
off: takes a counts matrix + sample sheet and runs:

  1. **QC**: per-sample library size, gene complexity, top-genes table.
  2. **Normalisation**: TPM / CPM / VST-style log scaling.
  3. **PCA + sample correlation heatmap** for batch sanity-checks.
  4. **Differential expression**:
       - If ``pydeseq2`` is installed, runs a proper DESeq2 fit.
       - Otherwise falls back to per-gene Welch's t-test on log-CPM
         + Benjamini-Hochberg FDR. Coarser than DESeq2 but doesn't
         pull a heavy dep.
  5. **Volcano + MA plots, top-genes heatmap**.
  6. **DEG → controlling cCRE / regulatory-element linkage** by
     querying the IGVF Catalog enhancer-gene endpoints for each
     significant DEG.

Subcommands

  qc            QC report for a counts matrix.
  pca           PCA + sample correlation heatmap.
  deg           Differential expression call.
  link-cre      For a DEG list, query the IGVF Catalog for the
                regulatory elements that control each gene.
  pipeline      Full QC -> PCA -> DEG -> link-cre in one command.
  write-playbook  Emit Docs/Skills/RNASEQ_ANALYSIS_SKILLS.md.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
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
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "RNAseq"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
MANIFEST_DIR = DATA_DIR / "Manifests" / "RNAseq"

CATALOG_API_BASE = _resolve_endpoint("catalog_api", "IGVF_CATALOG_API_BASE")

USER_AGENT = "IGVFdataAgent-RNAseq/0.1"


# --------------------------- Project plumbing ------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"rnaseq_analysis_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    logging.info("Log file: %s", log_path)
    return log_path


def mkdirs() -> None:
    for d in (REPORT_DIR, MANIFEST_DIR, SKILL_DOC_DIR):
        d.mkdir(parents=True, exist_ok=True)


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


# --------------------------- IO --------------------------------------------

def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def load_counts_matrix(path: Path) -> "tuple[Any, list[str], list[str]]":
    """Load a counts matrix (TSV/CSV, optionally gzipped) into a numpy
    array of shape (genes, samples). The first column is taken as the
    gene identifier and the first row as sample names.

    Returns (matrix, gene_ids, sample_ids).
    """
    try:
        import numpy as np
        import pandas as pd
    except ImportError:
        raise RuntimeError(
            "RNA-seq analysis needs `pandas` + `numpy`. "
            "Install with: pip install 'igvfagent[analysis]'"
        )
    sep = "\t" if str(path).endswith((".tsv", ".tsv.gz")) else None
    df = pd.read_csv(_open_text(path), sep=sep, engine="python",
                       index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    return df.to_numpy(dtype=float), list(df.index.astype(str)), list(df.columns.astype(str))


def load_sample_sheet(path: Path) -> "list[dict]":
    sep = "\t" if str(path).endswith((".tsv", ".tsv.gz")) else ","
    rows: "list[dict]" = []
    with _open_text(path) as f:
        rdr = csv.DictReader(f, delimiter=sep)
        for r in rdr:
            rows.append({k: (v or "").strip() for k, v in r.items()})
    return rows


def write_csv(path: Path, rows: "list[dict]",
                cols: "Optional[list[str]]" = None) -> None:
    if not rows:
        path.write_text("")
        return
    cols = cols or sorted({k for r in rows for k in r})
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# --------------------------- QC --------------------------------------------

def run_qc(counts, gene_ids: "list[str]",
              sample_ids: "list[str]") -> "list[dict]":
    """Per-sample library size, genes detected, top-fraction stats."""
    import numpy as np
    libsize = counts.sum(axis=0)
    detected = (counts > 0).sum(axis=0)
    top10 = np.zeros(counts.shape[1])
    for j in range(counts.shape[1]):
        col = counts[:, j]
        s = col.sum()
        if s <= 0:
            continue
        top10[j] = np.sort(col)[-10:].sum() / s
    rows = []
    for i, sid in enumerate(sample_ids):
        rows.append({
            "sample":         sid,
            "library_size":   int(libsize[i]),
            "genes_detected": int(detected[i]),
            "top10_frac":     float(top10[i]),
        })
    return rows


def cpm_log(counts) -> "Any":
    import numpy as np
    libsize = counts.sum(axis=0).astype(float)
    libsize[libsize == 0] = 1.0
    cpm = counts / libsize * 1e6
    return np.log2(cpm + 1)


# --------------------------- DEG -------------------------------------------

def deg_welch(counts, gene_ids, sample_ids, sample_sheet,
                 condition_col: str, group_a: str, group_b: str
                 ) -> "list[dict]":
    """Welch's t-test on log-CPM across samples, BH-adjusted."""
    import numpy as np
    from scipy import stats

    # Map sample -> condition
    by_sample = {r.get("sample") or r.get("gsm") or r.get("sample_id") or "":
                  r.get(condition_col, "")
                  for r in sample_sheet}
    a_idx = [i for i, s in enumerate(sample_ids)
             if by_sample.get(s) == group_a]
    b_idx = [i for i, s in enumerate(sample_ids)
             if by_sample.get(s) == group_b]
    if len(a_idx) < 2 or len(b_idx) < 2:
        raise SystemExit(
            f"Need at least 2 samples per group. Got "
            f"{len(a_idx)} {group_a!r} / {len(b_idx)} {group_b!r}. "
            f"Sample-sheet column was {condition_col!r}; the "
            f"matrix sample IDs are {sample_ids[:5]}{'...' if len(sample_ids)>5 else ''}."
        )

    log_cpm = cpm_log(counts)
    a = log_cpm[:, a_idx]
    b = log_cpm[:, b_idx]
    mean_a = a.mean(axis=1)
    mean_b = b.mean(axis=1)
    log2fc = mean_b - mean_a   # positive = higher in group_b

    t, p = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    p = np.where(np.isnan(p), 1.0, p)

    # BH correction
    n = len(p)
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    padj = np.minimum(1.0, p * n / ranks)
    # enforce monotonicity from largest p down
    sorted_idx = order[::-1]
    running_min = 1.0
    for i in sorted_idx:
        running_min = min(running_min, padj[i])
        padj[i] = running_min

    rows = []
    for i, gid in enumerate(gene_ids):
        rows.append({
            "gene":     gid,
            "log2FC":   float(log2fc[i]),
            "mean_a":   float(mean_a[i]),
            "mean_b":   float(mean_b[i]),
            "t_stat":   float(t[i]) if not np.isnan(t[i]) else 0.0,
            "p_value":  float(p[i]),
            "padj":     float(padj[i]),
        })
    return rows


def deg_pydeseq2(counts, gene_ids, sample_ids, sample_sheet,
                    condition_col: str, group_a: str, group_b: str
                    ) -> "list[dict]":
    """Optional pyDESeq2 path when the package is installed.
    Otherwise raises ImportError so caller falls back to Welch."""
    import numpy as np
    import pandas as pd
    from pydeseq2.dds import DeseqDataSet  # type: ignore
    from pydeseq2.ds import DeseqStats     # type: ignore

    by_sample = {r.get("sample") or r.get("gsm") or r.get("sample_id") or "":
                  r.get(condition_col, "")
                  for r in sample_sheet}
    keep = [i for i, s in enumerate(sample_ids)
             if by_sample.get(s) in (group_a, group_b)]
    if not keep:
        raise SystemExit(f"No samples match groups {group_a!r}/{group_b!r}.")
    sub = counts[:, keep]
    cols = [sample_ids[i] for i in keep]
    metadata = pd.DataFrame({"condition":
                                [by_sample.get(c) for c in cols]},
                                index=cols)
    expr = pd.DataFrame(sub.T, index=cols, columns=gene_ids).round().astype(int)
    dds = DeseqDataSet(counts=expr, metadata=metadata,
                          design_factors="condition", quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=("condition", group_b, group_a),
                      quiet=True)
    ds.summary()
    res = ds.results_df.reset_index().rename(
        columns={"baseMean": "base_mean", "log2FoldChange": "log2FC"})
    rows = res.to_dict(orient="records")
    return rows


def call_deg(counts, gene_ids, sample_ids, sample_sheet,
                condition_col: str, group_a: str, group_b: str,
                use_pydeseq2: bool = True) -> "list[dict]":
    if use_pydeseq2:
        try:
            return deg_pydeseq2(counts, gene_ids, sample_ids, sample_sheet,
                                  condition_col, group_a, group_b)
        except ImportError:
            logging.info("pydeseq2 not installed, falling back to Welch's t-test")
        except Exception as e:
            logging.warning("pydeseq2 failed (%s); falling back to Welch", e)
    return deg_welch(counts, gene_ids, sample_ids, sample_sheet,
                       condition_col, group_a, group_b)


# --------------------------- Plots -----------------------------------------

def plot_volcano(deg_rows: "list[dict]", out_path: Path,
                    label: str = "DEG",
                    fc_cut: float = 1.0, padj_cut: float = 0.05) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return
    if not deg_rows:
        return
    fc = np.array([r.get("log2FC", 0) for r in deg_rows], dtype=float)
    p = np.array([r.get("padj", 1.0) for r in deg_rows], dtype=float)
    p = np.clip(p, 1e-300, 1.0)
    nlp = -np.log10(p)
    plt.figure(figsize=(7, 6))
    sig_up = (fc >= fc_cut) & (p <= padj_cut)
    sig_dn = (fc <= -fc_cut) & (p <= padj_cut)
    other = ~(sig_up | sig_dn)
    plt.scatter(fc[other], nlp[other], s=4, c="#bbbbbb", alpha=0.5)
    plt.scatter(fc[sig_up], nlp[sig_up], s=8, c="#D55E00",
                  label=f"up (padj<{padj_cut}, |log2FC|>={fc_cut})")
    plt.scatter(fc[sig_dn], nlp[sig_dn], s=8, c="#0072B2",
                  label="down")
    plt.axvline(fc_cut, ls="--", color="#888", lw=0.5)
    plt.axvline(-fc_cut, ls="--", color="#888", lw=0.5)
    plt.axhline(-np.log10(padj_cut), ls="--", color="#888", lw=0.5)
    plt.xlabel("log2 fold change")
    plt.ylabel("-log10(padj)")
    plt.title(f"Volcano — {label}")
    plt.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    for ext in ("png", "svg"):
        plt.savefig(out_path.with_suffix(f".{ext}"), dpi=200,
                      bbox_inches="tight", facecolor="white")
    plt.close()


def plot_ma(deg_rows: "list[dict]", out_path: Path,
              label: str = "DEG", fc_cut: float = 1.0,
              padj_cut: float = 0.05) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return
    if not deg_rows:
        return
    fc = np.array([r.get("log2FC", 0) for r in deg_rows], dtype=float)
    a = np.array([(r.get("mean_a", 0) + r.get("mean_b", 0)) / 2
                    for r in deg_rows], dtype=float)
    if not (a > 0).any():
        a = np.array([r.get("base_mean", 0) for r in deg_rows], dtype=float)
        a = np.log2(a + 1)
    p = np.array([r.get("padj", 1.0) for r in deg_rows], dtype=float)
    sig = p <= padj_cut
    plt.figure(figsize=(7, 5))
    plt.scatter(a[~sig], fc[~sig], s=4, c="#bbbbbb", alpha=0.5)
    plt.scatter(a[sig], fc[sig], s=8, c="#D55E00",
                  label=f"padj<{padj_cut}")
    plt.axhline(0, color="#888", lw=0.5)
    plt.axhline(fc_cut, ls="--", color="#888", lw=0.5)
    plt.axhline(-fc_cut, ls="--", color="#888", lw=0.5)
    plt.xlabel("mean log2 expression")
    plt.ylabel("log2 fold change")
    plt.title(f"MA plot — {label}")
    plt.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    for ext in ("png", "svg"):
        plt.savefig(out_path.with_suffix(f".{ext}"), dpi=200,
                      bbox_inches="tight", facecolor="white")
    plt.close()


def plot_pca(counts, sample_ids, sample_sheet, condition_col,
                out_path: Path, label: str = "PCA") -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return
    log_cpm = cpm_log(counts).T   # (samples, genes)
    # Top variable genes
    var = log_cpm.var(axis=0)
    top = np.argsort(-var)[:1000]
    X = log_cpm[:, top]
    X = X - X.mean(axis=0, keepdims=True)
    # SVD-based PCA
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    pc = U * S
    by_sample = {r.get("sample") or r.get("gsm") or r.get("sample_id") or "":
                  r.get(condition_col, "(none)")
                  for r in sample_sheet}
    groups = [by_sample.get(s, "(none)") for s in sample_ids]
    uniq = sorted(set(groups))
    palette = {g: c for g, c in zip(uniq,
        ["#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9",
          "#CC79A7", "#7B5BA6"] * 4)}
    plt.figure(figsize=(7, 6))
    for g in uniq:
        mask = [x == g for x in groups]
        plt.scatter(pc[mask, 0], pc[mask, 1], s=60, color=palette[g],
                      label=g, edgecolor="black", linewidth=0.4)
    for i, sid in enumerate(sample_ids):
        plt.annotate(sid, (pc[i, 0], pc[i, 1]),
                       textcoords="offset points", xytext=(4, 4),
                       fontsize=7)
    var_explained = (S ** 2) / (S ** 2).sum()
    plt.xlabel(f"PC1 ({var_explained[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({var_explained[1]*100:.1f}%)")
    plt.title(f"PCA — {label}")
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    for ext in ("png", "svg"):
        plt.savefig(out_path.with_suffix(f".{ext}"), dpi=200,
                      bbox_inches="tight", facecolor="white")
    plt.close()


def plot_heatmap_top(counts, gene_ids, sample_ids, deg_rows,
                        out_path: Path, top_n: int = 50,
                        label: str = "Top DEGs") -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return
    if not deg_rows:
        return
    sorted_rows = sorted(deg_rows,
                            key=lambda r: r.get("padj", 1.0))[:top_n]
    keep_genes = [r["gene"] for r in sorted_rows]
    gene_to_idx = {g: i for i, g in enumerate(gene_ids)}
    idx = [gene_to_idx[g] for g in keep_genes if g in gene_to_idx]
    if not idx:
        return
    log_cpm = cpm_log(counts)
    sub = log_cpm[idx]
    sub = (sub - sub.mean(axis=1, keepdims=True))
    sub = sub / (sub.std(axis=1, keepdims=True) + 1e-9)
    plt.figure(figsize=(max(5.0, 0.4 * counts.shape[1]),
                          max(4.0, 0.18 * len(idx))))
    plt.imshow(sub, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    plt.colorbar(label="z-score (log2 CPM)")
    plt.yticks(range(len(idx)), [gene_ids[i] for i in idx], fontsize=6)
    plt.xticks(range(len(sample_ids)), sample_ids, rotation=70, fontsize=8)
    plt.title(f"{label} (top {len(idx)} by padj)")
    plt.tight_layout()
    for ext in ("png", "svg"):
        plt.savefig(out_path.with_suffix(f".{ext}"), dpi=200,
                      bbox_inches="tight", facecolor="white")
    plt.close()


# --------------------------- DEG → cCRE linkage ----------------------------

def http_get_json(url: str, timeout: int = 30) -> "tuple[int, Any]":
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT,
                                                  "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status, json.loads(resp.read())
    except urllib.error.HTTPError as e:
        return e.code, {}
    except Exception:
        return 0, {}


def link_deg_to_cre(deg_rows: "list[dict]",
                       padj_cut: float = 0.05,
                       fc_cut: float = 1.0,
                       limit_per_gene: int = 10,
                       max_genes: int = 50,
                       max_consecutive_failures: int = 3
                       ) -> "list[dict]":
    """For each significant DEG, query the IGVF Catalog
    /api/genes/genomic-elements endpoint to find regulatory elements
    that target it. Bounded by max_genes and consecutive-failure cap
    so a slow Catalog can never hang the analysis."""
    sig = [r for r in deg_rows
            if (r.get("padj", 1.0) <= padj_cut and
                abs(r.get("log2FC", 0.0)) >= fc_cut)]
    sig.sort(key=lambda r: r.get("padj", 1.0))
    sig = sig[:max_genes]
    out: "list[dict]" = []
    consec = 0
    for r in sig:
        gene = r["gene"]
        url = (CATALOG_API_BASE + "/api/genes/genomic-elements?"
                + urllib.parse.urlencode({"gene_name": gene,
                                            "limit": limit_per_gene,
                                            "verbose": "false"}))
        sc, data = http_get_json(url, timeout=20)
        if sc != 200:
            consec += 1
            if consec >= max_consecutive_failures:
                logging.warning(
                    "Catalog DEG-to-cCRE linkage aborted after %d "
                    "consecutive failures.", consec,
                )
                break
            continue
        consec = 0
        rows = (data if isinstance(data, list)
                else (data.get("results") or data.get("@graph") or []))
        for row in rows:
            if not isinstance(row, dict):
                continue
            ge = row.get("genomic_element") or ""
            element_id = (ge.split("/", 1)[1]
                            if isinstance(ge, str) and "/" in ge else ge)
            out.append({
                "gene":      gene,
                "log2FC":    r.get("log2FC", 0),
                "padj":      r.get("padj", 1.0),
                "element":   element_id,
                "method":    row.get("method", ""),
                "score":     row.get("score", ""),
                "source":    row.get("source", ""),
                "biological_context":
                    row.get("biological_context", ""),
            })
        time.sleep(0.05)
    return out


# --------------------------- Subcommands -----------------------------------

def cmd_qc(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    counts, gene_ids, sample_ids = load_counts_matrix(Path(args.counts))
    rows = run_qc(counts, gene_ids, sample_ids)
    ts = timestamp()
    label = safe_label(args.label or Path(args.counts).stem)
    out = REPORT_DIR / f"{ts}_qc_{label}"
    out.mkdir(parents=True, exist_ok=True)
    csv_out = out / f"{label}_qc.csv"
    write_csv(csv_out, rows)
    print(f"QC: {len(sample_ids)} samples × {len(gene_ids)} genes")
    print(f"CSV: {csv_out}")
    return csv_out


def cmd_pca(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    counts, gene_ids, sample_ids = load_counts_matrix(Path(args.counts))
    sheet = load_sample_sheet(Path(args.sample_sheet))
    ts = timestamp()
    label = safe_label(args.label or Path(args.counts).stem)
    out_dir = REPORT_DIR / f"{ts}_pca_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_pca(counts, sample_ids, sheet, args.condition_col,
                out_path=out_dir / f"{label}_pca.png", label=label)
    print(f"PCA plot: {out_dir / f'{label}_pca.png'}")
    return out_dir


def cmd_deg(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    counts, gene_ids, sample_ids = load_counts_matrix(Path(args.counts))
    sheet = load_sample_sheet(Path(args.sample_sheet))
    rows = call_deg(counts, gene_ids, sample_ids, sheet,
                       args.condition_col, args.group_a, args.group_b,
                       use_pydeseq2=not args.no_pydeseq2)
    ts = timestamp()
    label = safe_label(args.label or
                          f"deg_{args.group_b}_vs_{args.group_a}")
    out_dir = REPORT_DIR / f"{ts}_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_out = out_dir / f"{label}_deg.csv"
    write_csv(csv_out, rows)
    plot_volcano(rows, out_dir / f"{label}_volcano.png", label=label,
                    fc_cut=args.fc_cut, padj_cut=args.padj_cut)
    plot_ma(rows, out_dir / f"{label}_ma.png", label=label,
              fc_cut=args.fc_cut, padj_cut=args.padj_cut)
    plot_heatmap_top(counts, gene_ids, sample_ids, rows,
                          out_dir / f"{label}_heatmap_top.png",
                          top_n=args.heatmap_top, label=label)
    sig = [r for r in rows
            if r.get("padj", 1.0) <= args.padj_cut
            and abs(r.get("log2FC", 0.0)) >= args.fc_cut]
    print(f"DEGs found ({len(sig):,} significant of {len(rows):,} tested)")
    print(f"Table: {csv_out}")
    return csv_out


def cmd_link_cre(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    rows = list(csv.DictReader(open(args.deg)))
    for r in rows:
        for k in ("log2FC", "padj"):
            try:
                r[k] = float(r.get(k, 1.0))
            except ValueError:
                r[k] = 0.0 if k == "log2FC" else 1.0
    links = link_deg_to_cre(rows, padj_cut=args.padj_cut,
                                fc_cut=args.fc_cut,
                                limit_per_gene=args.limit_per_gene,
                                max_genes=args.max_genes)
    ts = timestamp()
    label = safe_label(args.label or "deg_to_cre")
    out_dir = REPORT_DIR / f"{ts}_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_out = out_dir / f"{label}_deg_to_cre.csv"
    write_csv(csv_out, links)
    print(f"DEG→cCRE rows: {len(links)} from "
          f"{len({r['gene'] for r in links})} genes")
    print(f"Table: {csv_out}")
    return csv_out


def cmd_pipeline(args: argparse.Namespace) -> Path:
    setup_logging(); mkdirs()
    ts = timestamp()
    label = safe_label(args.label or
                          f"rnaseq_{args.group_b}_vs_{args.group_a}")
    out_dir = REPORT_DIR / f"{ts}_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)

    counts, gene_ids, sample_ids = load_counts_matrix(Path(args.counts))
    sheet = load_sample_sheet(Path(args.sample_sheet))

    print(f"▶ QC")
    qc_rows = run_qc(counts, gene_ids, sample_ids)
    write_csv(out_dir / f"{label}_qc.csv", qc_rows)

    print(f"▶ PCA + correlation")
    plot_pca(counts, sample_ids, sheet, args.condition_col,
                out_path=out_dir / f"{label}_pca.png", label=label)

    print(f"▶ DEG ({args.group_b} vs {args.group_a})")
    rows = call_deg(counts, gene_ids, sample_ids, sheet,
                       args.condition_col, args.group_a, args.group_b,
                       use_pydeseq2=not args.no_pydeseq2)
    write_csv(out_dir / f"{label}_deg.csv", rows)
    plot_volcano(rows, out_dir / f"{label}_volcano.png", label=label,
                    fc_cut=args.fc_cut, padj_cut=args.padj_cut)
    plot_ma(rows, out_dir / f"{label}_ma.png", label=label,
              fc_cut=args.fc_cut, padj_cut=args.padj_cut)
    plot_heatmap_top(counts, gene_ids, sample_ids, rows,
                          out_dir / f"{label}_heatmap_top.png",
                          top_n=args.heatmap_top, label=label)

    sig = [r for r in rows
            if r.get("padj", 1.0) <= args.padj_cut
            and abs(r.get("log2FC", 0.0)) >= args.fc_cut]
    sig_up = [r for r in sig if r.get("log2FC", 0) >= args.fc_cut]
    sig_dn = [r for r in sig if r.get("log2FC", 0) <= -args.fc_cut]

    print(f"▶ Link DEGs to controlling cCREs (Catalog)")
    links: "list[dict]" = []
    if not args.skip_link_cre:
        links = link_deg_to_cre(rows, padj_cut=args.padj_cut,
                                    fc_cut=args.fc_cut,
                                    limit_per_gene=args.limit_per_gene,
                                    max_genes=args.max_link_genes)
        write_csv(out_dir / f"{label}_deg_to_cre.csv", links)

    report = out_dir / f"{label}_report.md"
    lines = [f"# RNA-seq pipeline — `{label}`",
             f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
             "",
             "## Inputs",
             f"- Counts: `{args.counts}`",
             f"- Sample sheet: `{args.sample_sheet}`",
             f"- Condition column: `{args.condition_col}`",
             f"- Compared: `{args.group_b}` vs `{args.group_a}`",
             "",
             "## Sizes",
             f"- Genes (rows): **{len(gene_ids):,}**",
             f"- Samples (cols): **{len(sample_ids):,}**",
             "",
             "## DEG summary",
             f"- Tested: **{len(rows):,}**",
             f"- Significant (padj < {args.padj_cut}, "
             f"|log2FC| >= {args.fc_cut}): **{len(sig):,}** "
             f"({len(sig_up)} up, {len(sig_dn)} down)",
             "",
             "## Top 25 up-regulated", "",
             "| gene | log2FC | padj |", "|---|---:|---:|"]
    for r in sorted(sig_up, key=lambda r: r["padj"])[:25]:
        lines.append(f"| {r['gene']} | {r['log2FC']:+.2f} | {r['padj']:.2e} |")
    lines += ["", "## Top 25 down-regulated", "",
               "| gene | log2FC | padj |", "|---|---:|---:|"]
    for r in sorted(sig_dn, key=lambda r: r["padj"])[:25]:
        lines.append(f"| {r['gene']} | {r['log2FC']:+.2f} | {r['padj']:.2e} |")
    if links:
        lines += ["", "## Top DEG → controlling cCRE",
                   f"_(Top {min(40, len(links))} of {len(links)} rows)_", "",
                   "| gene | element | method | score | bio_context |",
                   "|---|---|---|---:|---|"]
        for L in links[:40]:
            lines.append(
                f"| {L['gene']} | {L['element']} | {L['method']} "
                f"| {L['score']} | {L.get('biological_context','')} |"
            )
    lines += ["", "## Plots",
               f"- Volcano: `Plots/{label}_volcano.png`",
               f"- MA plot: `Plots/{label}_ma.png`",
               f"- Top-DEG heatmap: `Plots/{label}_heatmap_top.png`",
               f"- PCA: `Plots/{label}_pca.png`",
               ""]
    report.write_text("\n".join(lines))
    print(f"\n✓ Run dir:    {out_dir}")
    print(f"  Report:     {report}")
    return report


def cmd_write_playbook(_a) -> Path:
    mkdirs()
    path = SKILL_DOC_DIR / "RNASEQ_ANALYSIS_SKILLS.md"
    lines = [
        "# Skill: bulk RNA-seq analysis + DEG → controlling cCRE",
        "",
        "Picks up where the GEO retrieval skill leaves off: takes a "
        "counts matrix + sample sheet, runs QC, PCA, differential "
        "expression, and (optionally) links each significant DEG back "
        "to the regulatory elements that control it via the IGVF "
        "Catalog.",
        "",
        "DEG engine: pyDESeq2 if installed (preferred); otherwise "
        "Welch's t-test on log-CPM with Benjamini-Hochberg FDR. The "
        "fallback is fast and dependency-free; for publication-quality "
        "calls, install pyDESeq2 or use proper tools (DESeq2/edgeR/limma "
        "in R).",
        "",
        "## Subcommands",
        "",
        "### `pipeline` — one-command end-to-end",
        "",
        "```bash",
        "igvfagent rnaseq pipeline \\",
        "    --counts counts.tsv --sample-sheet samples.csv \\",
        "    --condition-col condition --group-a control --group-b treated \\",
        "    --label demo_run",
        "```",
        "",
        "Writes QC CSV + PCA + volcano + MA + heatmap + DEG CSV + (with "
        "Catalog reachable) DEG → controlling cCRE table + a summary "
        "markdown report under `Docs/RNAseq/<ts>_<label>/`.",
        "",
        "### `qc`, `pca`, `deg`, `link-cre` — individual steps",
        "",
        "```bash",
        "igvfagent rnaseq qc --counts counts.tsv",
        "igvfagent rnaseq pca --counts counts.tsv --sample-sheet samples.csv \\",
        "    --condition-col condition",
        "igvfagent rnaseq deg --counts counts.tsv --sample-sheet samples.csv \\",
        "    --condition-col condition --group-a control --group-b treated",
        "igvfagent rnaseq link-cre --deg <deg.csv> --padj-cut 0.05 --fc-cut 1.0",
        "```",
        "",
        "## How this chains with other skills",
        "",
        "1. `igvfagent geo series --gse GSE9574 --full-samples` →",
        "2. `igvfagent geo download --gse GSE9574 --only matrix` →",
        "3. `igvfagent geo sample-sheet --gse GSE9574` →",
        "4. `igvfagent rnaseq pipeline --counts <matrix> "
            "--sample-sheet <sheet>` →",
        "5. `igvfagent kg gene <SYMBOL> --depth 2 --call-literature` "
            "for any DEG you'd like to drill into.",
        "",
        "## Outputs",
        "",
        "- Per-sample QC: `Docs/RNAseq/<ts>_<label>/<label>_qc.csv`",
        "- DEG table:     `<label>_deg.csv`",
        "- Plots:         volcano / MA / heatmap / PCA "
        "(.png + .svg)",
        "- DEG → cCRE:    `<label>_deg_to_cre.csv` (when Catalog "
        "reachable)",
        "- Report:        `<label>_report.md`",
    ]
    path.write_text("\n".join(lines))
    print(f"Playbook: {path}")
    return path


# --------------------------------- CLI -------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(
        description="Bulk RNA-seq QC + DEG + DEG→cCRE linkage skill.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("qc", help="Counts QC report.")
    s.add_argument("--counts", required=True)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_qc)

    s = sub.add_parser("pca", help="Sample PCA + correlation heatmap.")
    s.add_argument("--counts", required=True)
    s.add_argument("--sample-sheet", required=True)
    s.add_argument("--condition-col", default="condition")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_pca)

    s = sub.add_parser("deg", help="Differential expression call.")
    s.add_argument("--counts", required=True)
    s.add_argument("--sample-sheet", required=True)
    s.add_argument("--condition-col", default="condition")
    s.add_argument("--group-a", required=True,
                    help="Control / reference group label.")
    s.add_argument("--group-b", required=True,
                    help="Treated / test group label.")
    s.add_argument("--padj-cut", type=float, default=0.05)
    s.add_argument("--fc-cut",   type=float, default=1.0)
    s.add_argument("--heatmap-top", type=int, default=50)
    s.add_argument("--no-pydeseq2", action="store_true",
                    help="Force the Welch's t-test fallback even when "
                         "pyDESeq2 is installed.")
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_deg)

    s = sub.add_parser("link-cre",
                        help="For a DEG CSV, query the IGVF Catalog "
                             "for the regulatory elements that control "
                             "each significant gene.")
    s.add_argument("--deg", required=True)
    s.add_argument("--padj-cut", type=float, default=0.05)
    s.add_argument("--fc-cut",   type=float, default=1.0)
    s.add_argument("--limit-per-gene", type=int, default=10)
    s.add_argument("--max-genes", type=int, default=50)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_link_cre)

    s = sub.add_parser("pipeline",
                        help="QC + PCA + DEG + DEG→cCRE in one run.")
    s.add_argument("--counts", required=True)
    s.add_argument("--sample-sheet", required=True)
    s.add_argument("--condition-col", default="condition")
    s.add_argument("--group-a", required=True)
    s.add_argument("--group-b", required=True)
    s.add_argument("--padj-cut", type=float, default=0.05)
    s.add_argument("--fc-cut",   type=float, default=1.0)
    s.add_argument("--heatmap-top", type=int, default=50)
    s.add_argument("--no-pydeseq2", action="store_true")
    s.add_argument("--skip-link-cre", action="store_true")
    s.add_argument("--limit-per-gene", type=int, default=10)
    s.add_argument("--max-link-genes", type=int, default=50)
    s.add_argument("--label", default="")
    s.set_defaults(func=cmd_pipeline)

    s = sub.add_parser("write-playbook",
                        help="Emit Docs/Skills/RNASEQ_ANALYSIS_SKILLS.md.")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
