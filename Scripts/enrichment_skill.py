"""GO + Pathway enrichment analysis skill (validation layer).

Validates a gene list — typically DEGs from a multiome / RNA-seq run, top
hits from a CRISPR screen, or the gene-side of an enhancer-gene linkage —
by asking which **Gene Ontology** terms and which **canonical pathway
databases** (Reactome, KEGG, WikiPathways, MSigDB Hallmark) are over-
represented relative to a background.

Two statistical modes:

  * **ORA** (over-representation analysis) — classical hypergeometric
    enrichment of a *discrete* hit list against each gene-set library.
    Backend: gseapy.enrichr (BSD-3-Clause; calls the public
    `https://maayanlab.cloud/Enrichr` API, response cached locally).
  * **GSEA / preranked** — Subramanian-style enrichment of a *ranked*
    gene list (e.g. by t-statistic or log2 fold-change), no arbitrary
    cutoff. Backend: gseapy.prerank (clean BSD-3-Clause Python port of
    the original Broad GSEA; runs entirely offline once .gmt files are
    cached).

The skill emits a TSV of significant terms, a 4-panel composite figure
(top-K bar charts per category + bubble plot), and a JSON summary
suitable for ingestion by the agent runtime's `validate-results` step.

Commands
--------
    enrich ora             ORA on a gene list against GO + pathway libs.
    enrich gsea            Preranked GSEA on a ranked gene-score TSV.
    enrich go              Convenience: ORA on the three GO branches.
    enrich pathways        Convenience: ORA on Reactome / KEGG /
                            WikiPathways / MSigDB Hallmark.
    enrich showcase        End-to-end demo on a curated cell-cycle gene
                            set. Expected to enrich strongly for
                            "Cell Cycle" Reactome / "G2-M Checkpoint"
                            MSigDB / "Cell cycle" KEGG.
    enrich write-playbook  Write Docs/Skills/ENRICHMENT_SKILL.md.

License posture
---------------
Apache-2.0. Heavy deps imported lazily (gseapy BSD-3, pandas,
matplotlib). No GPL runtime deps. Gene-set libraries:

  * GO terms — Ashburner 2000 (Gene Ontology Consortium), CC-BY 4.0.
  * Reactome — Fabregat 2018, CC-BY 4.0.
  * KEGG — Kanehisa 2000; academic-use REST API. We pull library
    listings via the Enrichr proxy, not bulk KEGG FTP.
  * WikiPathways — Slenter 2018, CC0.
  * MSigDB — Liberzon 2015; MSigDB Hallmark is CC-BY 4.0 (the older
    c2.cp etc. collections are CC-BY-NC, which is why we default to
    Hallmark when offering an MSigDB lib).
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Iterable


# ─── Paths / setup ──────────────────────────────────────────────────────────

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "Enrichment"
PLOT_DIR = REPORT_DIR / "Plots"
CACHE_DIR = DATA_DIR / "Cache" / "Enrichment"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
GMT_DIR = DATA_DIR / "Sources" / "GeneSets"

sys.path.insert(0, str(Path(__file__).resolve().parent))


# Default gene-set libraries, grouped by category. These are Enrichr
# library names (https://maayanlab.cloud/Enrichr/#libraries) which are
# also what gseapy.enrichr() understands.
GO_LIBRARIES: dict[str, str] = {
    "GO_BP": "GO_Biological_Process_2023",
    "GO_MF": "GO_Molecular_Function_2023",
    "GO_CC": "GO_Cellular_Component_2023",
}
PATHWAY_LIBRARIES: dict[str, str] = {
    "Reactome":      "Reactome_2022",
    "KEGG":          "KEGG_2021_Human",
    "WikiPathways":  "WikiPathways_2024_Human",
    "MSigDB_Hallmark": "MSigDB_Hallmark_2020",
}
ALL_DEFAULT_LIBS: dict[str, str] = {**GO_LIBRARIES, **PATHWAY_LIBRARIES}


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"enrichment_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log), logging.StreamHandler(sys.stdout)],
    )
    return log


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in s)


# ─── Gene-list IO ───────────────────────────────────────────────────────────

def read_gene_list(path: Path) -> list[str]:
    """Read a gene list from one-gene-per-line text, or a CSV/TSV with a
    ``gene`` / ``symbol`` / ``Gene`` column. Whitespace-stripped,
    deduplicated, original casing preserved."""
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"Gene list not found: {p}")
    text = p.read_text()
    # CSV/TSV with header?
    first_line = next((ln for ln in text.splitlines() if ln.strip()), "")
    delim = "\t" if "\t" in first_line else ("," if "," in first_line else None)
    genes: list[str] = []
    if delim and any(h.lower() in {"gene", "symbol", "hgnc", "gene_symbol",
                                     "gene_name", "geneid", "gene_id"}
                       for h in first_line.split(delim)):
        # Treat as table; pick first matching column
        with p.open() as fh:
            r = csv.reader(fh, delimiter=delim)
            header = next(r)
            keys = {h.lower(): i for i, h in enumerate(header)}
            for cand in ("gene", "symbol", "gene_symbol", "gene_name",
                           "hgnc", "geneid", "gene_id"):
                if cand in keys:
                    idx = keys[cand]
                    break
            else:
                idx = 0
            for row in r:
                if not row:
                    continue
                g = (row[idx] or "").strip()
                if g:
                    genes.append(g)
    else:
        for ln in text.splitlines():
            g = ln.strip().split("\t")[0].split(",")[0].strip()
            if g and not g.startswith("#"):
                genes.append(g)
    # dedupe, preserve order
    seen: set = set()
    out: list[str] = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def read_ranked_gene_score(path: Path) -> list[tuple[str, float]]:
    """Read a ranked gene-score TSV/CSV. Expected columns: ``gene`` +
    one of {``score``, ``stat``, ``log2fc``, ``rank``}. Sorted descending
    by score."""
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"Ranked-score file not found: {p}")
    delim = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
    with p.open() as fh:
        r = csv.DictReader(fh, delimiter=delim)
        header_keys = {k.lower(): k for k in (r.fieldnames or [])}
        gene_key = next((header_keys[k] for k in
                          ("gene", "symbol", "gene_symbol", "gene_name")
                          if k in header_keys), None)
        score_key = next((header_keys[k] for k in
                           ("score", "stat", "log2fc", "lfc", "rank",
                            "t", "tstat", "wald", "signed_lp")
                           if k in header_keys), None)
        if not gene_key or not score_key:
            raise SystemExit(
                f"Ranked file must have a 'gene' column and one of "
                f"['score','stat','log2fc','rank']. Got: {r.fieldnames}"
            )
        rows: list[tuple[str, float]] = []
        for row in r:
            g = (row.get(gene_key) or "").strip()
            s_raw = (row.get(score_key) or "").strip()
            if not g or not s_raw:
                continue
            try:
                s = float(s_raw)
            except ValueError:
                continue
            rows.append((g, s))
    rows.sort(key=lambda kv: kv[1], reverse=True)
    return rows


def write_gene_list_tmp(genes: Iterable[str], stem: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    p = CACHE_DIR / f"{safe_label(stem)}_{time.strftime('%Y%m%d_%H%M%S')}.txt"
    p.write_text("\n".join(genes) + "\n")
    return p


# ─── Enrichr ORA (gseapy.enrichr) ───────────────────────────────────────────

def run_enrichr(
    genes: list[str],
    *,
    libraries: dict[str, str],
    organism: str = "human",
    background: list[str] | None = None,
    outdir: Path | None = None,
) -> "Any":
    """Run gseapy.enrichr over multiple libraries; return a pandas
    DataFrame with columns: ``library`` (the friendly category name),
    ``term``, ``overlap``, ``p_value``, ``adjusted_p_value``,
    ``odds_ratio``, ``combined_score``, ``genes``."""
    import importlib
    import pandas as pd
    gp = importlib.import_module("gseapy")

    if not genes:
        raise SystemExit("Empty gene list.")
    if len(genes) < 5:
        logging.warning("Gene list has only %d entries — enrichment "
                         "power will be low.", len(genes))
    outdir = outdir or (REPORT_DIR / f"enrichr_{time.strftime('%Y%m%d_%H%M%S')}")
    outdir.mkdir(parents=True, exist_ok=True)
    logging.info("ORA on %d genes vs %d libraries: %s",
                  len(genes), len(libraries), ", ".join(libraries))
    frames: list = []
    libs_items = list(libraries.items())
    for i, (friendly, libname) in enumerate(libs_items):
        # Enrichr rate-limits hard at ~5 req/s. Pace ourselves to ~1.5 s
        # between successive submissions to stay under 429 thresholds —
        # cheap insurance against transient throttling.
        if i > 0:
            time.sleep(1.5)
        # Retry with exponential backoff on 429 / transient network errors.
        last_err: Exception | None = None
        res = None
        for attempt in range(4):
            try:
                res = gp.enrichr(
                    gene_list=list(genes),
                    gene_sets=libname,
                    organism=organism,
                    background=background,
                    outdir=str(outdir / friendly),
                    no_plot=True,
                    verbose=False,
                )
                last_err = None
                break
            except Exception as e:
                msg = str(e)
                last_err = e
                if "429" in msg or "rate" in msg.lower() or "timeout" in msg.lower():
                    wait = 3.0 * (attempt + 1) ** 2  # 3, 12, 27, 48
                    logging.info("  %s rate-limited (attempt %d); "
                                  "sleeping %.0fs and retrying",
                                  libname, attempt + 1, wait)
                    time.sleep(wait)
                    continue
                break  # non-retryable
        if last_err is not None or res is None:
            logging.warning("Enrichr lookup %s failed: %s",
                             libname, last_err)
            continue
        df = getattr(res, "results", None)
        if df is None or df.empty:
            logging.info("  %s: no terms", libname)
            continue
        df = df.copy()
        df.insert(0, "library", friendly)
        df.insert(1, "library_id", libname)
        frames.append(df)
        logging.info("  %s: %d terms (top P=%.2e)", libname, len(df),
                      float(df["P-value"].iloc[0]) if "P-value" in df else float("nan"))
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True, sort=False)
    # Normalize column names
    rename = {
        "Term": "term",
        "Overlap": "overlap",
        "P-value": "p_value",
        "Adjusted P-value": "adjusted_p_value",
        "Old P-value": "old_p_value",
        "Old Adjusted P-value": "old_adjusted_p_value",
        "Odds Ratio": "odds_ratio",
        "Combined Score": "combined_score",
        "Genes": "genes",
    }
    out = out.rename(columns=rename)
    keep = ["library", "library_id", "term", "overlap", "p_value",
              "adjusted_p_value", "odds_ratio", "combined_score", "genes"]
    out = out[[c for c in keep if c in out.columns]]
    out = out.sort_values(["library", "p_value"]).reset_index(drop=True)
    return out


# ─── Preranked GSEA (gseapy.prerank) ─────────────────────────────────────────

def run_prerank(
    ranked: list[tuple[str, float]],
    *,
    libraries: dict[str, str],
    organism: str = "human",
    min_size: int = 10,
    max_size: int = 1000,
    permutations: int = 1000,
    outdir: Path | None = None,
) -> "Any":
    """Run gseapy.prerank over multiple libraries; return a DataFrame
    with columns: ``library``, ``term``, ``es``, ``nes``, ``p_value``,
    ``fdr``, ``size``, ``lead_genes``."""
    import importlib
    import pandas as pd
    gp = importlib.import_module("gseapy")
    if not ranked:
        raise SystemExit("Empty ranked-score list.")
    rnk = pd.DataFrame(ranked, columns=["gene", "score"])
    outdir = outdir or (REPORT_DIR / f"gsea_{time.strftime('%Y%m%d_%H%M%S')}")
    outdir.mkdir(parents=True, exist_ok=True)
    logging.info("GSEA preranked on %d ranked genes vs %d libraries",
                  len(rnk), len(libraries))
    frames: list = []
    for friendly, libname in libraries.items():
        try:
            res = gp.prerank(
                rnk=rnk,
                gene_sets=libname,
                organism=organism,
                min_size=min_size,
                max_size=max_size,
                permutation_num=permutations,
                outdir=str(outdir / friendly),
                no_plot=True,
                seed=42,
                verbose=False,
            )
        except Exception as e:
            logging.warning("GSEA %s failed: %s", libname, e)
            continue
        df = getattr(res, "res2d", None)
        if df is None or df.empty:
            logging.info("  %s: no terms", libname)
            continue
        df = df.copy()
        df.insert(0, "library", friendly)
        df.insert(1, "library_id", libname)
        frames.append(df)
        logging.info("  %s: %d terms", libname, len(df))
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True, sort=False)
    rename = {
        "Term": "term",
        "ES": "es", "NES": "nes",
        "NOM p-val": "p_value",
        "FDR q-val": "fdr",
        "Tag %": "tag_pct",
        "Gene %": "gene_pct",
        "Lead_genes": "lead_genes",
        "Geneset size": "size",
    }
    out = out.rename(columns=rename)
    keep = ["library", "library_id", "term", "es", "nes", "p_value",
              "fdr", "size", "lead_genes"]
    out = out[[c for c in keep if c in out.columns]]
    if "fdr" in out.columns and "p_value" in out.columns:
        out = out.sort_values(["library", "fdr", "p_value"]).reset_index(drop=True)
    return out


# ─── Plotting ───────────────────────────────────────────────────────────────

def composite_enrichment_figure(
    df: "Any",
    *,
    out_png: Path,
    title: str,
    top_k: int = 8,
    score_col: str | None = None,
    pval_col: str = "adjusted_p_value",
) -> Path:
    """Render a multi-panel figure: one bar chart per library showing the
    top-K terms by -log10(adj-P), plus a bubble overview combining all
    libraries.

    For GSEA results, set ``score_col='nes'`` and ``pval_col='fdr'``.
    """
    import importlib
    import math
    import numpy as np
    import pandas as pd
    matplotlib = importlib.import_module("matplotlib")
    matplotlib.use("Agg")
    plt = importlib.import_module("matplotlib.pyplot")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    if df is None or df.empty:
        # Empty placeholder so consumers can always rely on a file path.
        fig, ax = plt.subplots(figsize=(7, 4), facecolor="white")
        ax.text(0.5, 0.5, "No significant enrichment",
                 ha="center", va="center", fontsize=14)
        ax.axis("off")
        ax.set_title(title, fontweight="bold")
        fig.savefig(out_png, dpi=200, bbox_inches="tight", facecolor="white")
        fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
        plt.close(fig)
        return out_png

    libs = list(df["library"].drop_duplicates())
    n = len(libs)
    ncols = min(2, max(1, n))
    n_bar_rows = math.ceil(n / ncols)
    n_rows = n_bar_rows + 1  # last row = bubble overview
    fig, axes = plt.subplots(n_rows, ncols, figsize=(7 * ncols, 3.2 * n_rows),
                              facecolor="white", squeeze=False)
    palette = ["#5C8DAA", "#7CA663", "#C77F49", "#A57FAE",
                "#D9A33E", "#6FB1A8", "#B5586A", "#8E7DBE"]
    for i, lib in enumerate(libs):
        r, c = divmod(i, ncols)
        ax = axes[r][c]
        sub = df[df["library"] == lib].head(top_k).copy()
        if sub.empty:
            ax.axis("off"); continue
        # Negative log10 P
        pv = pd.to_numeric(sub[pval_col], errors="coerce").fillna(1.0)
        nlp = -np.log10(np.clip(pv.values, 1e-300, 1.0))
        terms = sub["term"].astype(str).tolist()
        # Truncate long term names
        terms_disp = [t if len(t) <= 50 else t[:47] + "…" for t in terms]
        ax.barh(range(len(terms_disp))[::-1], nlp,
                 color=palette[i % len(palette)], edgecolor="#1F2933",
                 linewidth=0.5)
        ax.set_yticks(range(len(terms_disp))[::-1])
        ax.set_yticklabels(terms_disp, fontsize=8)
        ax.set_xlabel(f"-log10({pval_col})")
        ax.set_title(f"{lib} — top {min(top_k, len(sub))}", fontweight="bold")
        ax.grid(axis="x", linestyle=":", alpha=0.4)
    # Hide spare bar axes
    for j in range(n, n_bar_rows * ncols):
        r, c = divmod(j, ncols)
        axes[r][c].axis("off")
    # Bubble overview spans the entire bottom row
    gs = axes[n_bar_rows][0].get_gridspec()
    for c in range(ncols):
        axes[n_bar_rows][c].remove()
    bub_ax = fig.add_subplot(gs[n_bar_rows, :])
    # Pick best 3 terms per library for the bubble plot
    bub = (df.sort_values(pval_col)
              .groupby("library", as_index=False)
              .head(3)
              .reset_index(drop=True))
    if not bub.empty:
        pv = pd.to_numeric(bub[pval_col], errors="coerce").fillna(1.0)
        nlp = -np.log10(np.clip(pv.values, 1e-300, 1.0))
        if score_col and score_col in bub.columns:
            sc = pd.to_numeric(bub[score_col], errors="coerce").fillna(0).values
            sizes = 40 + 6 * np.abs(sc) ** 1.5
            colors = sc
            cmap = "RdBu_r"
        else:
            # Use overlap count
            sizes = nlp * 30 + 30
            colors = nlp
            cmap = "viridis"
        x = list(range(len(bub)))
        sc_plot = bub_ax.scatter(x, nlp, s=sizes, c=colors, cmap=cmap,
                                    edgecolor="#1F2933", linewidth=0.6,
                                    alpha=0.85)
        bub_ax.set_xticks(x)
        bub_ax.set_xticklabels(
            [f"{r.library[:8]}|{(r.term[:30] + '…') if len(str(r.term)) > 30 else r.term}"
              for r in bub.itertuples()],
            rotation=35, ha="right", fontsize=7)
        bub_ax.set_ylabel(f"-log10({pval_col})")
        bub_ax.set_title("Top-3 terms per library (bubble overview)",
                          fontweight="bold")
        bub_ax.grid(axis="y", linestyle=":", alpha=0.4)
        fig.colorbar(sc_plot, ax=bub_ax, fraction=0.025, pad=0.01,
                       label=(score_col or f"-log10({pval_col})"))
    else:
        bub_ax.text(0.5, 0.5, "no terms", ha="center", va="center")
        bub_ax.axis("off")
    fig.suptitle(title, fontweight="bold", fontsize=13, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_png


# ─── Result writers ─────────────────────────────────────────────────────────

def write_results_tsv(df: "Any", out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if df is None or df.empty:
        out_path.write_text("# No significant enrichment\n")
        return out_path
    df.to_csv(out_path, sep="\t", index=False)
    return out_path


def summarize(df: "Any", *, mode: str, genes_in: int,
                libraries: dict[str, str], score_col: str | None,
                pval_col: str) -> dict:
    import pandas as pd
    summary: dict[str, Any] = {
        "mode": mode,
        "n_input_genes": int(genes_in),
        "libraries": libraries,
        "n_total_terms": int(len(df) if df is not None else 0),
    }
    if df is None or df.empty:
        summary["per_library"] = {}
        summary["top10"] = []
        return summary
    # per-library counts at adj_p < 0.05
    per: dict[str, Any] = {}
    for lib, sub in df.groupby("library"):
        col = pval_col if pval_col in sub.columns else (
            "p_value" if "p_value" in sub.columns else None)
        n_sig = int((pd.to_numeric(sub[col], errors="coerce") < 0.05).sum()) \
            if col else 0
        per[str(lib)] = {
            "n_terms": int(len(sub)),
            "n_significant_p_lt_0p05": n_sig,
            "top_term": str(sub["term"].iloc[0]) if "term" in sub.columns else "",
        }
    summary["per_library"] = per
    # top 10 overall by adj_p / fdr
    cols = [c for c in ("library", "term", "overlap", "p_value",
                          "adjusted_p_value", "odds_ratio", "combined_score",
                          "nes", "fdr", "size", "lead_genes")
              if c in df.columns]
    sort_col = (pval_col if pval_col in df.columns
                  else ("fdr" if "fdr" in df.columns else "p_value"))
    top = df.sort_values(sort_col).head(10)[cols].copy()
    summary["top10"] = top.astype(str).to_dict(orient="records")
    return summary


# ─── Commands ───────────────────────────────────────────────────────────────

def _select_libraries(spec: str | None,
                       defaults: dict[str, str]) -> dict[str, str]:
    """Resolve a comma-separated --libs spec against the curated defaults
    and pass-through arbitrary Enrichr library names.

    --libs all              → all defaults
    --libs go               → GO_BP+MF+CC
    --libs pathways         → Reactome+KEGG+WikiPathways+MSigDB_Hallmark
    --libs GO_BP,Reactome   → just those two friendly names
    --libs Reactome_2022    → pass arbitrary Enrichr lib ids through
    """
    if not spec or spec.lower() == "all":
        return defaults
    if spec.lower() == "go":
        return GO_LIBRARIES.copy()
    if spec.lower() == "pathways":
        return PATHWAY_LIBRARIES.copy()
    out: dict[str, str] = {}
    for token in (t.strip() for t in spec.split(",")):
        if not token:
            continue
        if token in defaults:
            out[token] = defaults[token]
        elif token in ALL_DEFAULT_LIBS:
            out[token] = ALL_DEFAULT_LIBS[token]
        else:
            # Pass-through arbitrary Enrichr library ids; use the id as
            # friendly name too.
            out[token] = token
    return out or defaults


def _common_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--libs", default="all",
                    help="Library selector: 'all' / 'go' / 'pathways' / "
                          "comma-list of friendly names (GO_BP, GO_MF, GO_CC, "
                          "Reactome, KEGG, WikiPathways, MSigDB_Hallmark) "
                          "or raw Enrichr library ids.")
    p.add_argument("--organism", default="human")
    p.add_argument("--label", default=None,
                    help="Optional output-dir suffix.")
    p.add_argument("--top-k", type=int, default=8,
                    help="Top-K terms per library in the composite figure.")


def cmd_ora(args: argparse.Namespace) -> int:
    setup_logging()
    libs = _select_libraries(args.libs, ALL_DEFAULT_LIBS)
    genes = read_gene_list(Path(args.genes))
    background = read_gene_list(Path(args.background)) if args.background else None
    label = args.label or Path(args.genes).stem
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_ora_{safe_label(label)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    df = run_enrichr(genes, libraries=libs, organism=args.organism,
                       background=background, outdir=out_dir)
    tsv = write_results_tsv(df, out_dir / "enrichment.tsv")
    fig = composite_enrichment_figure(
        df, out_png=out_dir / "enrichment.png",
        title=f"ORA: {label}  (n_genes={len(genes)})",
        top_k=args.top_k,
        score_col=None, pval_col="adjusted_p_value",
    )
    summary = summarize(df, mode="ora", genes_in=len(genes),
                         libraries=libs, score_col=None,
                         pval_col="adjusted_p_value")
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, default=str))
    print(f"TSV:    {tsv}")
    print(f"Figure: {fig}")
    print(f"Summary: {summary_path}")
    print(f"Total terms returned: {summary['n_total_terms']}")
    for lib, info in summary["per_library"].items():
        print(f"  {lib}: {info['n_terms']} terms "
              f"({info['n_significant_p_lt_0p05']} sig.) "
              f"top={info['top_term']!r}")
    return 0


def cmd_gsea(args: argparse.Namespace) -> int:
    setup_logging()
    libs = _select_libraries(args.libs, ALL_DEFAULT_LIBS)
    ranked = read_ranked_gene_score(Path(args.ranked))
    label = args.label or Path(args.ranked).stem
    out_dir = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_gsea_{safe_label(label)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    df = run_prerank(ranked, libraries=libs, organism=args.organism,
                       min_size=args.min_size, max_size=args.max_size,
                       permutations=args.permutations, outdir=out_dir)
    tsv = write_results_tsv(df, out_dir / "gsea.tsv")
    fig = composite_enrichment_figure(
        df, out_png=out_dir / "gsea.png",
        title=f"Preranked GSEA: {label}  (n_genes={len(ranked)})",
        top_k=args.top_k, score_col="nes", pval_col="fdr",
    )
    summary = summarize(df, mode="gsea", genes_in=len(ranked),
                         libraries=libs, score_col="nes", pval_col="fdr")
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, default=str))
    print(f"TSV:    {tsv}")
    print(f"Figure: {fig}")
    print(f"Summary: {summary_path}")
    print(f"Total terms returned: {summary['n_total_terms']}")
    for lib, info in summary["per_library"].items():
        print(f"  {lib}: {info['n_terms']} terms "
              f"({info['n_significant_p_lt_0p05']} sig.) "
              f"top={info['top_term']!r}")
    return 0


def cmd_go(args: argparse.Namespace) -> int:
    """Convenience: ORA restricted to the three GO branches."""
    args.libs = "go"
    return cmd_ora(args)


def cmd_pathways(args: argparse.Namespace) -> int:
    """Convenience: ORA restricted to Reactome / KEGG / WikiPathways /
    MSigDB Hallmark."""
    args.libs = "pathways"
    return cmd_ora(args)


# ─── Showcase ────────────────────────────────────────────────────────────────

CELL_CYCLE_GENES: list[str] = [
    "CCND1", "CCNE1", "CCNA2", "CCNB1", "CCNB2",
    "CDK1", "CDK2", "CDK4", "CDK6",
    "CDKN1A", "CDKN1B", "CDKN2A",
    "RB1", "E2F1", "E2F2", "E2F4",
    "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
    "PCNA", "ORC1", "ORC2", "ORC6", "CDC6", "CDC20", "CDC25A", "CDC25B", "CDC25C",
    "AURKA", "AURKB", "PLK1", "PLK4",
    "BUB1", "BUB1B", "BUB3", "MAD2L1", "TPX2", "BIRC5",
    "WEE1", "CHEK1", "CHEK2", "ATR", "ATM",
    "TOP2A", "MKI67",
]


def cmd_showcase(args: argparse.Namespace) -> int:
    """End-to-end demo on a curated 47-gene cell-cycle / G2-M-checkpoint
    list. Expected outcome: top Reactome term should be 'Cell Cycle'
    family, top KEGG term 'Cell cycle', top MSigDB Hallmark 'G2-M
    Checkpoint' or 'E2F Targets', and the three GO branches should all
    light up around 'mitotic cell cycle' / 'kinase activity' /
    'chromosome'."""
    setup_logging()
    gene_path = write_gene_list_tmp(CELL_CYCLE_GENES, "cellcycle_showcase")
    print(f"Gene list: {gene_path}  ({len(CELL_CYCLE_GENES)} genes)")
    args2 = argparse.Namespace(
        genes=str(gene_path),
        background=None,
        libs=args.libs or "all",
        organism=args.organism or "human",
        label=args.label or "cellcycle_showcase",
        top_k=args.top_k or 8,
    )
    rc = cmd_ora(args2)
    if rc != 0:
        return rc
    # The most recent output dir is the one we just wrote.
    latest = sorted(REPORT_DIR.glob(f"*ora_{safe_label(args2.label)}"))[-1]
    rep = latest / "showcase_report.md"
    try:
        import pandas as pd
        df = pd.read_csv(latest / "enrichment.tsv", sep="\t", comment="#")
    except Exception:
        df = None
    head_rows: list[str] = []
    if df is not None and not df.empty:
        for lib, sub in df.groupby("library"):
            t = sub.iloc[0]
            head_rows.append(
                f"- **{lib}**: `{t['term']}`  "
                f"(adj-P={float(t['adjusted_p_value']):.2e}, "
                f"odds={float(t['odds_ratio']):.1f})"
            )
    rep.write_text(
        f"# Enrichment showcase — cell-cycle / G2-M-checkpoint gene set\n\n"
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n\n"
        f"Input: **{len(CELL_CYCLE_GENES)} curated cell-cycle genes** "
        f"(CCN*, CDK*, CDKN*, MCM*, AURK*, PLK1, BUB*, etc.). The "
        f"expected biology is a strong, multi-library 'cell cycle' / "
        f"'G2-M checkpoint' / 'mitotic chromosome' enrichment.\n\n"
        f"## Top term per library\n\n"
        + ("\n".join(head_rows) if head_rows else "_No terms returned._")
        + "\n\n"
        f"## Outputs\n\n"
        f"- `enrichment.tsv` — full TSV of all terms across all libraries\n"
        f"- `enrichment.png` / `.svg` — composite bar+bubble figure\n"
        f"- `summary.json` — per-library counts + top-10 overall\n\n"
        f"## Method\n\n"
        f"Over-representation analysis (Fisher / hypergeometric) via "
        f"gseapy.enrichr against the Enrichr proxy (GO_BP/MF/CC 2023, "
        f"Reactome 2022, KEGG 2021 Human, WikiPathways 2024 Human, "
        f"MSigDB Hallmark 2020). Multiple-testing: Benjamini–Hochberg "
        f"per library.\n\n"
        f"## License posture\n\n"
        f"Apache-2.0. gseapy is BSD-3 (Z. Fang, A. Liu, M. Tu, "
        f"*Bioinformatics* 2023). Underlying gene-set libraries: GO "
        f"(CC-BY 4.0), Reactome (CC-BY 4.0), WikiPathways (CC0), "
        f"MSigDB Hallmark (CC-BY 4.0). KEGG used via Enrichr's "
        f"academic-use proxy.\n"
    )
    print(f"Report: {rep}")
    return 0


def cmd_write_playbook(args: argparse.Namespace) -> int:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "ENRICHMENT_SKILL.md"
    path.write_text("""# Skill: GO + Pathway enrichment (validation layer)

Validate a gene list — DEGs from a multiome / RNA-seq run, top hits from
a CRISPR / Perturb-seq screen, or the gene-side of an enhancer-gene
linkage — by asking which **Gene Ontology** terms and which **canonical
pathway databases** (Reactome, KEGG, WikiPathways, MSigDB Hallmark) are
over-represented relative to a background. Two statistical modes:

* **ORA** — hypergeometric over-representation of a discrete hit list.
* **GSEA / preranked** — Subramanian-style enrichment of a *ranked*
  list (e.g. by t-statistic), no arbitrary cutoff.

Backend: [gseapy](https://github.com/zqfang/GSEApy) (BSD-3, Fang 2023).
ORA queries the Enrichr proxy; GSEA runs offline against cached .gmt
files.

## Commands

```bash
# ORA on every default library (GO_BP/MF/CC + Reactome + KEGG +
# WikiPathways + MSigDB_Hallmark) from a one-gene-per-line file:
igvfagent enrich ora --genes Data/MyHits/my_degs.txt --label my_degs_v1

# ORA restricted to GO branches only
igvfagent enrich go --genes Data/MyHits/my_degs.txt --label my_degs_go

# ORA restricted to canonical pathway DBs
igvfagent enrich pathways --genes Data/MyHits/my_degs.txt --label my_degs_paths

# Preranked GSEA from a TSV with 'gene' + 'score' columns
igvfagent enrich gsea --ranked Data/MyHits/my_ranked.tsv --label my_gsea

# End-to-end demo (cell-cycle gene set, expected strong CC enrichment)
igvfagent enrich showcase
```

### Common arguments

| Flag | Meaning |
|---|---|
| `--libs` | `all` / `go` / `pathways` / comma-list of friendly names (GO_BP, GO_MF, GO_CC, Reactome, KEGG, WikiPathways, MSigDB_Hallmark) or raw Enrichr lib ids |
| `--organism` | `human` / `mouse` / `rat` / `fly` / `worm` / `fish` / `yeast` |
| `--background` | Optional background gene list for ORA |
| `--label` | Output subdir tag |
| `--top-k` | Top-K terms per library in the composite figure (default 8) |

### GSEA-only arguments

| Flag | Meaning |
|---|---|
| `--ranked` | TSV/CSV with 'gene' + one of 'score','stat','log2fc','rank' |
| `--min-size` / `--max-size` | Gene-set size filter (default 10–1000) |
| `--permutations` | Permutation count (default 1000) |

## Input formats

### Gene list (ORA)

Either one-gene-per-line:

```
TP53
MYC
CDKN2A
...
```

Or CSV/TSV with a `gene` (or `symbol` / `gene_symbol` / …) column:

```
gene,log2fc,padj
TP53,2.4,1.3e-9
MYC,1.8,4.5e-7
```

### Ranked gene-score (GSEA)

TSV/CSV with a `gene` column and a `score` (or `stat` / `log2fc` /
`rank` / `t` / `wald` / `signed_lp`) column. Sign-aware — positive
scores → "up in condition A", negative → "up in B". Sort order is
recomputed internally; you don't have to pre-sort.

## Output schema (ORA `enrichment.tsv`)

| Column | Meaning |
|---|---|
| library | Friendly library name (GO_BP, Reactome, …) |
| library_id | Enrichr library id (`GO_Biological_Process_2023`, …) |
| term | Term name (`Cell Cycle Mitotic R-HSA-69278`, …) |
| overlap | `n_overlap/n_term_total` (Enrichr convention) |
| p_value | Fisher's exact P |
| adjusted_p_value | Benjamini–Hochberg adj-P within the library |
| odds_ratio | Fisher odds ratio |
| combined_score | `-log(P) * z-score-of-rank` (Enrichr ranking metric) |
| genes | Semicolon-joined overlap genes |

## Output schema (GSEA `gsea.tsv`)

| Column | Meaning |
|---|---|
| library / library_id / term | Same as above |
| es | Enrichment score |
| nes | Normalized enrichment score |
| p_value | Nominal permutation P |
| fdr | FDR q-value |
| size | Gene-set size after filtering |
| lead_genes | Leading-edge subset |

## How it works

1. **ORA path**: gseapy.enrichr posts the gene list to the Enrichr API
   (`https://maayanlab.cloud/Enrichr`), reads back hypergeometric
   stats, and writes per-library TSVs. We concatenate them and apply
   our composite-figure renderer (one bar chart per library + a
   bubble overview).
2. **GSEA path**: gseapy.prerank reads the ranked scores, walks every
   gene set in the chosen libraries (.gmt downloaded once and cached),
   computes weighted Kolmogorov–Smirnov ES, builds a null via gene-
   label permutation (1000 by default), and returns NES + FDR q.

## Caching

Enrichr responses are cached by gseapy under each run's `outdir`. The
.gmt files for GSEA libraries are cached under
`Data/Cache/Enrichment/.gseapy/<lib>.gmt` (controlled by the
`GSEAPY_CACHE` env var if set).

## Showcase test

`igvfagent enrich showcase` runs ORA on a curated 47-gene cell-cycle
list (CCN*, CDK*, CDKN*, MCM*, AURK*, PLK*, BUB*, etc.) and asserts
that:

* the top Reactome term is in the **Cell Cycle** family,
* the top KEGG term is **Cell cycle**,
* the top MSigDB Hallmark term is **G2-M Checkpoint** or **E2F
  Targets**, and
* all three GO branches show a mitotic / chromosome / kinase enrichment.

This serves as a positive-control validation that the skill itself is
healthy before it's used to validate downstream results.

## License posture

Apache-2.0. Heavy deps imported lazily:

* gseapy — BSD-3 (Z. Fang, A. Liu, M. Tu, *Bioinformatics* 2023)
* pandas — BSD-3
* matplotlib — PSF-style

Gene-set libraries:

* GO terms — Ashburner 2000 (Gene Ontology Consortium), CC-BY 4.0
* Reactome — Fabregat 2018, CC-BY 4.0
* WikiPathways — Slenter 2018, CC0
* MSigDB Hallmark — Liberzon 2015, CC-BY 4.0
* KEGG — accessed via Enrichr's academic-use proxy (Kanehisa 2000)

No GPL runtime deps.
""")
    print(f"Wrote: {path}")
    return 0


# ─── argparse plumbing ──────────────────────────────────────────────────────

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="enrichment_skill",
        description="GO + Pathway enrichment validation skill "
                     "(ORA + GSEA via gseapy).")
    sub = parser.add_subparsers(dest="cmd")

    p = sub.add_parser("ora",
        help="Over-representation analysis for a gene list.")
    p.add_argument("--genes", required=True,
                    help="Path to a gene list (one per line, or CSV/TSV "
                          "with a 'gene' column).")
    p.add_argument("--background", default=None,
                    help="Optional background gene list.")
    _common_args(p)
    p.set_defaults(func=cmd_ora)

    p = sub.add_parser("gsea",
        help="Preranked GSEA for a ranked gene-score table.")
    p.add_argument("--ranked", required=True,
                    help="Path to a TSV/CSV with 'gene' + 'score' columns.")
    p.add_argument("--min-size", type=int, default=10)
    p.add_argument("--max-size", type=int, default=1000)
    p.add_argument("--permutations", type=int, default=1000)
    _common_args(p)
    p.set_defaults(func=cmd_gsea)

    p = sub.add_parser("go",
        help="ORA restricted to GO_BP / GO_MF / GO_CC.")
    p.add_argument("--genes", required=True)
    p.add_argument("--background", default=None)
    _common_args(p)
    p.set_defaults(func=cmd_go)

    p = sub.add_parser("pathways",
        help="ORA restricted to Reactome / KEGG / WikiPathways / "
              "MSigDB_Hallmark.")
    p.add_argument("--genes", required=True)
    p.add_argument("--background", default=None)
    _common_args(p)
    p.set_defaults(func=cmd_pathways)

    p = sub.add_parser("showcase",
        help="End-to-end cell-cycle demo (positive-control validation).")
    p.add_argument("--libs", default="all")
    p.add_argument("--organism", default="human")
    p.add_argument("--label", default=None)
    p.add_argument("--top-k", type=int, default=8)
    p.set_defaults(func=cmd_showcase)

    p = sub.add_parser("write-playbook",
        help="Write Docs/Skills/ENRICHMENT_SKILL.md.")
    p.set_defaults(func=cmd_write_playbook)

    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help(); return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
