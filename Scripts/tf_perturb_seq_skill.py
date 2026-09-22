#!/usr/bin/env python3
"""IGVF TF Perturb-seq: calibrated cis / trans effects -> disease and GWAS overlay.

Port of the IGVF TF Perturb-seq consortium's Working Group 3 (disease & GWAS)
jamboree notebook, ``perturb_seq_analysis_v4.ipynb`` in
https://github.com/IGVF/tf_perturb_seq (docs/jamborees/2026_UTSW/
working_groups/wg3_disease_gwas). The notebook was consulted line by line and
its analysis translated into a skill; no notebook code was copied verbatim.
Relationship: port.

The upstream analysis, as run on the Huangfu-lab HUES8 definitive-endoderm
CRISPRi TF Perturb-seq (2,161 targets, 269,491 cells over eight IGVF
measurement sets), in order:

  1. guide-library overview from ``inference_mudata.h5mu`` (guide types,
     guides per cell, guides per target, batches);
  2. the IGVF CRISPR pipeline's *calibrated* inference tables -- direct
     target (self-knockdown), cis (nearby genes), trans (genome-wide) -- with
     empirical, multiple-testing-adjusted p-values;
  3. significance filters (adjusted p < 0.05; |log2FC| > 0.2 cis / > 1.0
     trans; targeting elements only) and distribution plots;
  4. top trans regulators and their strongest up / down targets;
  5. GWAS Catalog positional overlap: SNPs within a window of the perturbed
     TF elements and of the trans target genes, then Fisher trait enrichment
     against the whole catalog with BH correction;
  6. scE2G enhancer-gene predictions (ESC and definitive-endoderm models):
     which perturbed TFs and trans targets carry E2G links, and which E2G
     elements contain GWAS SNPs -- including E2G assignments that disagree
     with the catalog's nearest-gene mapping;
  7. optional ChIP-seq support: bigWig signal at each GWAS position, called
     present above a fraction of the track maximum.

WHERE THIS DIFFERS FROM THE NOTEBOOK, so nobody mistakes the two:

  * Positional overlap is vectorised (sorted positions + searchsorted per
    chromosome) instead of a Python loop over every gene against the whole
    catalog. Same answer; minutes become seconds.
  * The 22-million-row trans table is filtered in chunks and the filtered
    rows cached under Data/Cache/TFPerturbSeq, keyed on the file and the
    thresholds, so a re-run with the same cutoffs does not re-read 3 GB.
  * BH adjustment is computed here (the same procedure statsmodels'
    ``fdr_bh`` runs), so statsmodels is not a dependency.
  * The Fisher enrichment counts catalog ASSOCIATIONS, exactly as upstream
    does. Several SNPs in one locus, or one SNP reported by several studies,
    each count once, so odds ratios for traits mapped to a handful of loci
    (the memory / Alzheimer's block around APOE in the upstream run) are
    inflated. That is a property of the method, reproduced faithfully and
    stated in the report rather than silently changed.
  * ``mudata`` is optional: the .h5mu is read through h5py + anndata when
    it is not installed. ``pyBigWig`` is optional and only needed for the
    ChIP step.

    igvfagent tf-perturb run --calibrated-prefix <dir>/<run>_calibrated_ \\
        [--mudata inference_mudata.h5mu] [--gwas gwas_catalog_associations.tsv] \\
        [--e2g ESC=<h7.e2g.tsv> --e2g DE=<de.e2g.tsv>] [--bigwig <a.bigwig> ...] \\
        [--gene-coords genes.tsv] [--label huangfu_de]
    igvfagent tf-perturb selftest        # synthetic inputs, planted signals

Writes ``Docs/TFPerturbSeq/<timestamp>_<label>/`` with report.md, the filtered
result tables, every overlap and enrichment table as CSV, and PNG figures.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import re
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT")
            or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "TFPerturbSeq"
CACHE_DIR = ROOT / "Data" / "Cache" / "TFPerturbSeq"

# Palette: the dataviz reference instance. Down-regulation blue, up orange.
BLUE, ORANGE, INK, INK2, AXIS, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"

CAL_SUFFIXES = {"direct": "direct_target_results.tsv",
                "cis": "cis_results.tsv",
                "trans": "trans_results.tsv"}
CAL_COLUMNS = ["element_id", "element_symbol", "element_label", "tested_gene_id",
               "tested_gene_symbol", "n_cells", "log2fc", "log2fc_se", "is_cis",
               "is_direct_target", "posterior_pval", "empirical_pval",
               "empirical_pval_adj"]
GWAS_COLUMNS = ["CHR_ID", "CHR_POS", "DISEASE/TRAIT", "SNPS", "P-VALUE", "MAPPED_GENE"]
E2G_COLUMNS = ["ElementChr", "ElementStart", "ElementEnd", "ElementName",
               "GeneSymbol", "Score", "ElementClass", "isSelfPromoter"]

_ELEMENT_RE = re.compile(r"^(?P<gene>[^|]+)\|(?P<chr>chr[\w]+):(?P<start>\d+)-(?P<end>\d+)$")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"tf_perturb_seq_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def _pd():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit("tf-perturb needs pandas/numpy/scipy/matplotlib: "
                         "pip install 'igvfagent[analysis]'") from e


def read_table(path: "str | Path", usecols: "Optional[list[str]]" = None, **kw):
    pd = _pd()
    p = Path(path)
    if p.suffix == ".parquet":
        return pd.read_parquet(p, columns=usecols)
    sep = "," if p.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(p, sep=sep, usecols=usecols, comment=kw.pop("comment", None),
                       low_memory=False, **kw)


def write_csv(df, path: Path) -> Path:
    df.to_csv(path, index=False)
    return path


def bh_adjust(pvals) -> Any:
    """Benjamini-Hochberg, the same monotone step-up statsmodels' fdr_bh uses."""
    import numpy as np
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def md_table(headers: "list[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 30) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for i, r in enumerate(rows):
        if i >= max_rows:
            lines.append(f"| … | {'|'.join(' ' for _ in headers[1:])} |")
            break
        lines.append("| " + " | ".join(_fmt(c) for c in r) + " |")
    return "\n".join(lines)


def _fmt(v: Any) -> str:
    if isinstance(v, float):
        if v != v:
            return ""
        if abs(v) < 1e-3 and v != 0:
            return f"{v:.2e}"
        return f"{v:,.3f}".rstrip("0").rstrip(".") if abs(v) < 1e6 else f"{v:,.0f}"
    if isinstance(v, int):
        return f"{v:,}"
    return str(v).replace("|", "\\|")


# ---------------------------------------------------------------------------
# 1. MuData: guide library overview (optional input)
# ---------------------------------------------------------------------------

def load_mudata_meta(path: Path) -> dict:
    """Guide/gene metadata and the guide-assignment layer, with or without mudata.

    An .h5mu is HDF5 with each modality stored AnnData-style under ``mod/``,
    so anndata's element reader gets everything the overview needs when the
    mudata package is absent.
    """
    import numpy as np
    out: dict = {"path": str(path)}
    try:
        import mudata as md  # type: ignore
        m = md.read_h5mu(str(path))
        gene_var, guide_var = m.mod["gene"].var, m.mod["guide"].var
        obs = m.obs
        assignment = m.mod["guide"].layers.get("guide_assignment")
        n_cells = m.n_obs
    except ImportError:
        import h5py  # type: ignore
        from anndata.experimental import read_elem  # type: ignore
        with h5py.File(str(path), "r") as h:
            gene_var = read_elem(h["mod/gene/var"])
            guide_var = read_elem(h["mod/guide/var"])
            obs = read_elem(h["obs"]) if "obs" in h else read_elem(h["mod/gene/obs"])
            assignment = (read_elem(h["mod/guide/layers/guide_assignment"])
                          if "mod/guide/layers/guide_assignment" in h else None)
            n_cells = obs.shape[0]
    out["n_cells"], out["n_genes"], out["n_guides"] = int(n_cells), int(gene_var.shape[0]), int(guide_var.shape[0])
    out["gene_var"], out["guide_var"] = gene_var, guide_var
    out["batches"] = obs["batch"].value_counts().to_dict() if "batch" in obs.columns else {}
    out["guide_types"] = guide_var["type"].value_counts().to_dict() if "type" in guide_var.columns else {}
    if "intended_target_name" in guide_var.columns:
        out["n_targets"] = int(guide_var["intended_target_name"].nunique())
        targeting = guide_var[guide_var.get("targeting", True) == True]  # noqa: E712
        gpt = targeting.groupby("intended_target_name", observed=True).size()
        out["guides_per_target_median"] = float(gpt.median()) if len(gpt) else float("nan")
        out["guides_per_target_mean"] = float(gpt.mean()) if len(gpt) else float("nan")
    if assignment is not None:
        gpc = np.asarray(assignment.sum(axis=1)).ravel()
        out["cells_with_guide"] = int((gpc > 0).sum())
        out["guides_per_cell_median"] = float(np.median(gpc))
        out["guides_per_cell_mean"] = float(gpc.mean())
    # Ensembl id -> symbol, from both modalities, as the notebook builds it.
    id_to_symbol = {}
    if "symbol" in gene_var.columns:
        id_to_symbol.update(gene_var["symbol"].to_dict())
    if {"intended_target_name", "gene_name"} <= set(guide_var.columns):
        id_to_symbol.update(guide_var[["intended_target_name", "gene_name"]].drop_duplicates()
                            .set_index("intended_target_name")["gene_name"].to_dict())
    out["id_to_symbol"] = id_to_symbol
    return out


def gene_coords_from(meta: "Optional[dict]", coords_path: "Optional[Path]"):
    """Gene coordinates: gene_id, symbol, chr (no 'chr' prefix), start, end."""
    pd = _pd()
    if coords_path is not None:
        df = read_table(coords_path)
        cols = {c.lower(): c for c in df.columns}
        pick = {k: cols[k] for k in ("gene_id", "symbol", "chr", "start", "end") if k in cols}
        for alt, k in (("gene_chr", "chr"), ("gene_start", "start"), ("gene_end", "end"),
                       ("chrom", "chr"), ("ensembl_id", "gene_id"), ("gene_symbol", "symbol")):
            if k not in pick and alt in cols:
                pick[k] = cols[alt]
        missing = {"symbol", "chr", "start", "end"} - set(pick)
        if missing:
            raise SystemExit(f"--gene-coords lacks columns {sorted(missing)}")
        out = pd.DataFrame({k: df[v] for k, v in pick.items()})
        if "gene_id" not in out.columns:
            out["gene_id"] = out["symbol"]
    elif meta is not None and {"gene_chr", "gene_start", "gene_end"} <= set(meta["gene_var"].columns):
        gv = meta["gene_var"]
        out = pd.DataFrame({"gene_id": gv.index.astype(str), "symbol": gv["symbol"].astype(str),
                            "chr": gv["gene_chr"], "start": gv["gene_start"], "end": gv["gene_end"]})
    else:
        return None
    out = out[out["chr"].notna() & out["start"].notna() & out["end"].notna()].copy()
    out["chr"] = out["chr"].astype(str).str.replace("^chr", "", regex=True)
    out["start"] = out["start"].astype(float).astype(int)
    out["end"] = out["end"].astype(float).astype(int)
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# 2-3. Calibrated tables and significance filters
# ---------------------------------------------------------------------------

def _cache_key(path: Path, *parts: Any) -> str:
    st = path.stat()
    raw = "|".join([str(path.resolve()), str(st.st_size), str(int(st.st_mtime)), *map(str, parts)])
    return hashlib.sha1(raw.encode()).hexdigest()[:16]


def load_trans_significant(path: Path, padj: float, lfc_min: float,
                           chunksize: int = 2_000_000, use_cache: bool = True):
    """Rows of the trans table passing the cutoffs, read in chunks and cached."""
    pd = _pd()
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache = CACHE_DIR / f"trans_sig_{_cache_key(path, padj, lfc_min)}.parquet"
    if use_cache and cache.is_file():
        logging.info("trans: cache hit %s", cache)
        return pd.read_parquet(cache), True
    keep = []
    n_total = 0
    if path.suffix == ".parquet":
        chunks = [pd.read_parquet(path)]
    else:
        chunks = pd.read_csv(path, sep="\t", chunksize=chunksize, low_memory=False)
    for chunk in chunks:
        n_total += len(chunk)
        m = (chunk["empirical_pval_adj"] < padj) & (chunk["log2fc"].abs() > lfc_min)
        keep.append(chunk[m])
    sig = pd.concat(keep, ignore_index=True) if keep else pd.DataFrame(columns=CAL_COLUMNS)
    sig.attrs["n_total"] = n_total
    try:
        sig.to_parquet(cache, index=False)
    except Exception as exc:  # pragma: no cover - pyarrow missing
        logging.warning("trans cache not written: %s", exc)
    return sig, False


def filter_targeting(df, targeting_symbols: "Optional[set]" = None):
    """Targeting elements only: drops non-targeting / negative / positive controls."""
    out = df[df["element_label"].astype(str) == "targeting"]
    if targeting_symbols:
        out = out[out["element_symbol"].isin(targeting_symbols)]
    return out.copy()


def parse_element_coords(df):
    """element_id 'ENSG…|chr6:41072917-41073528' -> el_chr / el_start / el_end."""
    parts = df["element_id"].astype(str).str.extract(r"^(?:[^|]+)\|chr(\w+):(\d+)-(\d+)$")
    out = df.copy()
    out["el_chr"] = parts[0].values
    out["el_start"] = _pd().to_numeric(parts[1], errors="coerce").values
    out["el_end"] = _pd().to_numeric(parts[2], errors="coerce").values
    return out[out["el_chr"].notna()]


# ---------------------------------------------------------------------------
# 5. GWAS Catalog overlap and trait enrichment
# ---------------------------------------------------------------------------

def load_gwas(path: Path):
    pd = _pd()
    g = read_table(path, usecols=lambda c: c in GWAS_COLUMNS)
    g = g[g["CHR_ID"].notna() & g["CHR_POS"].notna()].copy()
    g["CHR_POS"] = pd.to_numeric(g["CHR_POS"], errors="coerce")
    g = g[g["CHR_POS"].notna()].copy()
    g["CHR_POS"] = g["CHR_POS"].astype("int64")
    g["CHR_ID"] = g["CHR_ID"].astype(str).str.replace("^chr", "", regex=True)
    g["P-VALUE"] = pd.to_numeric(g["P-VALUE"], errors="coerce")
    return g.reset_index(drop=True)


class PositionIndex:
    """Per-chromosome sorted SNP positions for interval queries."""

    def __init__(self, gwas) -> None:
        import numpy as np
        self.by_chr: dict = {}
        for chrom, sub in gwas.groupby("CHR_ID", sort=False):
            pos = sub["CHR_POS"].to_numpy()
            order = np.argsort(pos, kind="stable")
            self.by_chr[str(chrom)] = (pos[order], sub.index.to_numpy()[order])

    def query(self, chrom: str, start: int, end: int, window: int = 0):
        """Row indices of SNPs with start-window <= pos <= end+window."""
        import numpy as np
        entry = self.by_chr.get(str(chrom).replace("chr", ""))
        if entry is None:
            return np.empty(0, dtype=int)
        pos, idx = entry
        lo = np.searchsorted(pos, start - window, side="left")
        hi = np.searchsorted(pos, end + window, side="right")
        return idx[lo:hi]


def overlap_intervals(index: PositionIndex, gwas, intervals, key_cols: "list[str]",
                      chr_col: str, start_col: str, end_col: str, window: int):
    """Every (interval, SNP) pair within `window` bp; one row per pair."""
    import numpy as np
    pd = _pd()
    frames = []
    for _, row in intervals.iterrows():
        hits = index.query(row[chr_col], int(row[start_col]), int(row[end_col]), window)
        if hits.size == 0:
            continue
        sub = gwas.loc[hits, ["SNPS", "DISEASE/TRAIT", "P-VALUE", "MAPPED_GENE", "CHR_POS"]].copy()
        for k in key_cols:
            sub[k] = row[k]
        frames.append(sub)
    if not frames:
        return pd.DataFrame(columns=key_cols + ["snp", "trait", "pvalue", "mapped_gene", "pos"])
    out = pd.concat(frames, ignore_index=True)
    out = out.rename(columns={"SNPS": "snp", "DISEASE/TRAIT": "trait", "P-VALUE": "pvalue",
                              "MAPPED_GENE": "mapped_gene", "CHR_POS": "pos"})
    return out[key_cols + ["snp", "trait", "pvalue", "mapped_gene", "pos"]]


def trait_enrichment(hits, gwas, min_hits: int = 5):
    """Fisher exact (greater) per trait, hits vs the whole catalog; BH FDR.

    The contingency table is the notebook's: observed associations near the
    query set vs. the trait's catalog total, against all other associations.
    """
    from scipy.stats import fisher_exact  # type: ignore
    pd = _pd()
    trait_counts = gwas["DISEASE/TRAIT"].value_counts()
    total = int(len(gwas))
    n_hits = int(len(hits))
    hit_traits = hits["trait"].value_counts()
    hit_traits = hit_traits[hit_traits >= min_hits]
    rows = []
    for trait, observed in hit_traits.items():
        catalog_n = int(trait_counts.get(trait, 1))
        observed = int(observed)
        expected = catalog_n * (n_hits / total) if total else float("nan")
        a, b = observed, max(catalog_n - observed, 0)
        c = max(n_hits - observed, 0)
        d = max(total - catalog_n - c, 0)
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({"trait": trait, "observed": observed, "catalog_total": catalog_n,
                     "expected": round(expected, 2),
                     "enrichment": round(observed / max(expected, 0.01), 2),
                     "odds_ratio": (round(float(odds), 2) if odds != float("inf") else float("inf")),
                     "pvalue": float(p)})
    edf = pd.DataFrame(rows, columns=["trait", "observed", "catalog_total", "expected",
                                      "enrichment", "odds_ratio", "pvalue"])
    if len(edf):
        edf["fdr"] = bh_adjust(edf["pvalue"].to_numpy())
        edf = edf.sort_values(["fdr", "pvalue"]).reset_index(drop=True)
    else:
        edf["fdr"] = []
    return edf


# ---------------------------------------------------------------------------
# 6. scE2G predictions
# ---------------------------------------------------------------------------

def load_e2g(specs: "list[str]", score_min: float):
    """`LABEL=path` specs -> one combined frame with a `source` column."""
    pd = _pd()
    frames = []
    for spec in specs:
        label, _, p = spec.partition("=")
        if not p:
            label, p = Path(spec).stem.split(".")[0][:12], spec
        df = read_table(p, comment="#")
        missing = {"GeneSymbol", "Score", "ElementChr", "ElementStart", "ElementEnd"} - set(df.columns)
        if missing:
            raise SystemExit(f"E2G file {p} lacks columns {sorted(missing)}")
        if "ElementName" not in df.columns:
            df["ElementName"] = (df["ElementChr"].astype(str) + ":" + df["ElementStart"].astype(str)
                                 + "-" + df["ElementEnd"].astype(str))
        df["source"] = label
        df["n_total"] = len(df)
        frames.append(df[df["Score"] >= score_min].copy())
    if not frames:
        return None
    out = pd.concat(frames, ignore_index=True)
    out["chr_clean"] = out["ElementChr"].astype(str).str.replace("^chr", "", regex=True)
    return out


# ---------------------------------------------------------------------------
# 7. bigWig ChIP support (optional)
# ---------------------------------------------------------------------------

def bigwig_support(loci, bigwigs: "list[str]", frac_of_max: float):
    """Signal at each (chrom, pos) in every bigWig; pass = >= frac * track max."""
    try:
        import pyBigWig  # type: ignore
    except ImportError:
        raise RuntimeError("--bigwig needs pyBigWig: pip install 'igvfagent[hic]'")
    res = loci.copy()
    summary = []
    for bwf in bigwigs:
        bw = pyBigWig.open(str(bwf))
        bw_max = float(bw.header().get("maxVal") or 0.0)
        thr = bw_max * frac_of_max
        label = re.sub(r"\.(bigwig|bw)$", "", Path(bwf).name, flags=re.I)
        vals = []
        for chrom, pos in zip(res["chrom"], res["pos"]):
            try:
                v = bw.values(str(chrom), int(pos), int(pos) + 1)[0]
            except Exception:
                v = None
            vals.append(0.0 if v is None or v != v else float(v))
        bw.close()
        res[f"{label}_signal"] = vals
        res[f"{label}_pass"] = [v >= thr for v in vals]
        summary.append({"track": label, "max": bw_max, "threshold": thr,
                        "passing": int(sum(res[f"{label}_pass"])), "n": len(res)})
    pass_cols = [c for c in res.columns if c.endswith("_pass")]
    res["n_chipseq_passing"] = res[pass_cols].sum(axis=1) if pass_cols else 0
    return res, summary


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _plt():
    import matplotlib  # type: ignore
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore
    return plt


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


def fig_distributions(cis, trans_sig, padj: float, lfc_cis: float, out_dir: Path) -> "list[Path]":
    import numpy as np
    plt = _plt()
    floor = 1e-300
    paths = []
    if cis is not None and len(cis):
        fig, axes = plt.subplots(2, 3, figsize=(15, 8))
        fig.suptitle("Cis results", fontsize=13, fontweight="bold", color=INK, x=0.02, ha="left")
        fin = cis[np.isfinite(cis["log2fc"])]
        axes[0, 0].hist(-np.log10(cis["empirical_pval"].clip(lower=floor)), bins=50, color=BLUE)
        _style(axes[0, 0], "Empirical p-values", "-log10(empirical p)", "count")
        axes[0, 1].hist(-np.log10(cis["empirical_pval_adj"].clip(lower=floor)), bins=50, color=BLUE)
        _style(axes[0, 1], "Adjusted p-values", "-log10(adj p)", "count")
        axes[0, 2].hist(fin["log2fc"], bins=50, color=BLUE)
        _style(axes[0, 2], "Effect sizes", "log2FC", "count")
        ps = np.sort(cis["empirical_pval_adj"].to_numpy())
        axes[1, 0].plot(np.arange(1, len(ps) + 1), -np.log10(np.clip(ps, floor, None)), lw=0.8, color=BLUE)
        axes[1, 0].axhline(-np.log10(padj), color=ORANGE, ls="--", lw=0.8)
        axes[1, 0].text(0.98, 0.95, f"{int((ps < padj).sum()):,} significant", transform=axes[1, 0].transAxes,
                        ha="right", va="top", fontsize=8, color=INK)
        _style(axes[1, 0], "Adj p-value rank", "rank", "-log10(adj p)")
        pp = np.sort(cis["posterior_pval"].to_numpy())
        axes[1, 1].plot(np.arange(1, len(pp) + 1), -np.log10(np.clip(pp, floor, None)), lw=0.8, color=BLUE)
        _style(axes[1, 1], "Posterior p-value rank", "rank", "-log10(posterior p)")
        lfc = np.sort(fin["log2fc"].to_numpy())
        axes[1, 2].scatter(np.arange(len(lfc)), lfc, s=2, c=np.where(lfc < 0, BLUE, ORANGE), alpha=0.6, lw=0)
        for y in (lfc_cis, -lfc_cis):
            axes[1, 2].axhline(y, color=INK2, ls="--", lw=0.6)
        _style(axes[1, 2], "log2FC rank (blue down, orange up)", "rank", "log2FC")
        fig.patch.set_facecolor(SURFACE)
        fig.tight_layout()
        p = out_dir / "cis_distributions.png"
        fig.savefig(p, dpi=140)
        plt.close(fig)
        paths.append(p)
    if trans_sig is not None and len(trans_sig):
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        fig.suptitle("Trans results (significant subset)", fontsize=13, fontweight="bold", color=INK, x=0.02, ha="left")
        axes[0].hist(-np.log10(trans_sig["empirical_pval_adj"].clip(lower=floor)), bins=50, color=BLUE)
        _style(axes[0], "Adjusted p-values", "-log10(adj p)", "count")
        axes[1].hist(trans_sig["log2fc"], bins=50, color=BLUE)
        _style(axes[1], "Effect sizes", "log2FC", "count")
        lfc = np.sort(trans_sig["log2fc"].to_numpy())
        axes[2].scatter(np.arange(len(lfc)), lfc, s=2, c=np.where(lfc < 0, BLUE, ORANGE), alpha=0.6, lw=0)
        axes[2].axhline(0, color=AXIS, lw=0.6)
        _style(axes[2], "log2FC rank (blue down, orange up)", "rank", "log2FC")
        fig.patch.set_facecolor(SURFACE)
        fig.tight_layout()
        p = out_dir / "trans_distributions.png"
        fig.savefig(p, dpi=140)
        plt.close(fig)
        paths.append(p)
    return paths


def fig_hbar(items: "list[tuple[str, int | float]]", title: str, xlabel: str, path: Path,
             color: str = BLUE) -> "Optional[Path]":
    if not items:
        return None
    plt = _plt()
    names = [str(n)[:40] for n, _ in items][::-1]
    vals = [float(v) for _, v in items][::-1]
    fig, ax = plt.subplots(figsize=(9, max(2.5, 0.3 * len(items) + 1)))
    ax.barh(names, vals, color=color, height=0.72)
    for i, v in enumerate(vals):
        ax.text(v, i, f" {v:,.0f}" if float(v).is_integer() else f" {v:,.2f}", va="center", fontsize=8, color=INK)
    _style(ax, title, xlabel)
    fig.patch.set_facecolor(SURFACE)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def fig_e2g(e2g, score_min: float, out_dir: Path) -> "list[Path]":
    import numpy as np
    plt = _plt()
    paths = []
    sources = list(dict.fromkeys(e2g["source"]))
    fig, axes = plt.subplots(len(sources), 3, figsize=(15, 4 * len(sources)), squeeze=False)
    for r, src in enumerate(sources):
        sub = e2g[e2g["source"] == src]
        axes[r, 0].hist(sub["Score"], bins=50, color=BLUE)
        axes[r, 0].axvline(score_min, color=ORANGE, ls="--", lw=0.8)
        _style(axes[r, 0], f"{src}: score distribution", "E2G score", "count")
        axes[r, 1].hist(sub["Score"], bins=50, cumulative=True, density=True, color=BLUE)
        _style(axes[r, 1], f"{src}: cumulative", "E2G score", "fraction")
        if "ElementClass" in sub.columns:
            vc = sub["ElementClass"].value_counts()
            axes[r, 2].bar(vc.index.astype(str), vc.values, color=BLUE, width=0.72)
            _style(axes[r, 2], f"{src}: element classes", "", "count")
    fig.patch.set_facecolor(SURFACE)
    fig.tight_layout()
    p = out_dir / "e2g_distributions.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)
    fig, axes = plt.subplots(1, len(sources), figsize=(6 * len(sources), 4), squeeze=False)
    for i, src in enumerate(sources):
        s = np.sort(e2g.loc[e2g["source"] == src, "Score"].to_numpy())[::-1]
        axes[0, i].plot(np.arange(len(s)), s, lw=0.8, color=BLUE)
        axes[0, i].axhline(score_min, color=ORANGE, ls="--", lw=0.8)
        axes[0, i].text(0.98, 0.95, f"{len(s):,} links >= {score_min}", transform=axes[0, i].transAxes,
                        ha="right", va="top", fontsize=8, color=INK)
        _style(axes[0, i], f"{src}: score rank", "rank", "E2G score")
    fig.patch.set_facecolor(SURFACE)
    fig.tight_layout()
    p = out_dir / "e2g_score_ranks.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)
    return paths


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run(args: argparse.Namespace) -> int:
    import numpy as np
    pd = _pd()
    log_path = setup_logging()
    ts = time.strftime("%Y%m%d_%H%M%S")
    out_dir = OUT_ROOT / f"{ts}_{safe_label(args.label or 'tf_perturb')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    R: "list[str]" = []           # report lines
    tables: "list[Path]" = []
    figures: "list[Path]" = []
    summary: dict = {"label": args.label, "thresholds": {
        "padj": args.padj, "lfc_cis": args.lfc_cis, "lfc_trans": args.lfc_trans,
        "gwas_window": args.gwas_window, "e2g_score_min": args.e2g_score_min,
        "signal_frac_of_max": args.signal_frac}}

    R.append("# IGVF TF Perturb-seq: cis / trans effects with disease and GWAS overlay")
    R.append("")
    R.append(f"Generated {time.strftime('%Y-%m-%d %H:%M')} by `igvfagent tf-perturb run`; port of the IGVF "
             f"tf_perturb_seq WG3 notebook `perturb_seq_analysis_v4.ipynb`. Thresholds: adjusted p < {args.padj}, "
             f"|log2FC| > {args.lfc_cis} (cis / direct) and > {args.lfc_trans} (trans), GWAS window "
             f"{args.gwas_window:,} bp, E2G score >= {args.e2g_score_min}.")
    R.append("")

    # ---- 1. MuData overview -------------------------------------------------
    meta = None
    if args.mudata:
        meta = load_mudata_meta(Path(args.mudata))
        summary["mudata"] = {k: v for k, v in meta.items() if k not in ("gene_var", "guide_var", "id_to_symbol")}
        R.append("## 1. Guide library and cells")
        R.append("")
        R.append(md_table(["metric", "value"], [
            ["cells", meta["n_cells"]], ["genes (expression modality)", meta["n_genes"]],
            ["guides", meta["n_guides"]], ["unique intended targets", meta.get("n_targets", "")],
            ["cells with >= 1 guide", meta.get("cells_with_guide", "")],
            ["guides per cell (median / mean)",
             f"{meta.get('guides_per_cell_median', float('nan')):.0f} / {meta.get('guides_per_cell_mean', float('nan')):.2f}"],
            ["guides per target (median / mean)",
             f"{meta.get('guides_per_target_median', float('nan')):.0f} / {meta.get('guides_per_target_mean', float('nan')):.1f}"]]))
        R.append("")
        if meta["guide_types"]:
            R.append(md_table(["guide type", "n"], sorted(meta["guide_types"].items(), key=lambda kv: -kv[1])))
            R.append("")
        if meta["batches"]:
            R.append(md_table(["batch (IGVF measurement set)", "cells"],
                              sorted(meta["batches"].items(), key=lambda kv: -kv[1])))
            R.append("")
    targeting_symbols = None
    if meta is not None and {"targeting", "intended_target_name"} <= set(meta["guide_var"].columns):
        gv = meta["guide_var"]
        elems = set(gv.loc[gv["targeting"] == True, "intended_target_name"].astype(str))  # noqa: E712
        targeting_symbols = {meta["id_to_symbol"].get(e, e) for e in elems}

    # ---- 2-3. Calibrated tables ---------------------------------------------
    prefix = args.calibrated_prefix
    paths = {k: Path(prefix + v) for k, v in CAL_SUFFIXES.items()}
    for k, p in paths.items():
        if not p.is_file():
            alt = p.with_suffix(".parquet")
            if alt.is_file():
                paths[k] = alt
            else:
                print(f"ERROR: calibrated {k} table not found: {p}", file=sys.stderr)
                return 2
    direct = read_table(paths["direct"])
    cis = read_table(paths["cis"])
    trans_sig, cached = load_trans_significant(paths["trans"], args.padj, args.lfc_trans,
                                               use_cache=not args.no_cache)
    n_trans_total = trans_sig.attrs.get("n_total") if hasattr(trans_sig, "attrs") else None

    direct_sig = direct[(direct["empirical_pval_adj"] < args.padj) & (direct["log2fc"].abs() > args.lfc_cis)].copy()
    cis_sig = cis[(cis["empirical_pval_adj"] < args.padj) & (cis["log2fc"].abs() > args.lfc_cis)].copy()
    n_trans_before = len(trans_sig)
    trans_sig = filter_targeting(trans_sig, targeting_symbols)
    summary["counts"] = {
        "direct_total": int(len(direct)), "direct_sig": int(len(direct_sig)),
        "cis_total": int(len(cis)), "cis_sig": int(len(cis_sig)),
        "trans_total": n_trans_total, "trans_sig_before_targeting_filter": int(n_trans_before),
        "trans_sig": int(len(trans_sig)),
        "trans_regulators": int(trans_sig["element_symbol"].nunique()),
        "trans_target_genes": int(trans_sig["tested_gene_symbol"].nunique())}
    for name, df in (("direct_sig", direct_sig), ("cis_sig", cis_sig), ("trans_sig", trans_sig)):
        tables.append(write_csv(df, out_dir / f"{name}.csv"))

    R.append("## 2. Calibrated results and significance filters")
    R.append("")
    R.append("Tables are the IGVF CRISPR pipeline's calibrated inference output; `empirical_pval_adj` is "
             "already corrected for multiple testing. Trans rows are restricted to targeting elements "
             "(non-targeting, negative and positive controls dropped).")
    R.append("")
    R.append(md_table(["table", "rows", "significant", "unique elements", "unique genes"], [
        ["direct target (self-knockdown)", len(direct), len(direct_sig),
         direct_sig["element_symbol"].nunique(), direct_sig["tested_gene_symbol"].nunique()],
        ["cis", len(cis), len(cis_sig), cis_sig["element_symbol"].nunique(), cis_sig["tested_gene_symbol"].nunique()],
        ["trans", n_trans_total if n_trans_total else "(cached filter)", len(trans_sig),
         trans_sig["element_symbol"].nunique(), trans_sig["tested_gene_symbol"].nunique()]]))
    R.append("")
    R.append("### Strongest direct-target knockdowns")
    R.append("")
    top_kd = direct_sig.sort_values("empirical_pval_adj").head(20)
    R.append(md_table(["element", "gene", "log2FC", "adj p", "cells"],
                      top_kd[["element_symbol", "tested_gene_symbol", "log2fc", "empirical_pval_adj", "n_cells"]].values.tolist()))
    R.append("")
    R.append("### Strongest cis effects")
    R.append("")
    top_cis = cis_sig.sort_values("empirical_pval_adj").head(15)
    R.append(md_table(["element", "gene", "log2FC", "adj p", "cells"],
                      top_cis[["element_symbol", "tested_gene_symbol", "log2fc", "empirical_pval_adj", "n_cells"]].values.tolist()))
    R.append("")
    if not args.no_plots:
        figures += fig_distributions(cis, trans_sig, args.padj, args.lfc_cis, out_dir)

    # ---- 4. Top trans regulators --------------------------------------------
    reg_counts = trans_sig.groupby("element_symbol").size().sort_values(ascending=False)
    top_reg = reg_counts.head(args.top_regulators)
    tables.append(write_csv(reg_counts.rename("n_trans_targets").reset_index(), out_dir / "trans_regulators.csv"))
    R.append("## 3. Top trans regulators")
    R.append("")
    R.append(md_table(["TF element", "significant trans targets"], list(top_reg.items())))
    R.append("")
    for sym in list(top_reg.index[:3]):
        sub = trans_sig[trans_sig["element_symbol"] == sym]
        up = sub[sub["log2fc"] > 0].nlargest(5, "log2fc")
        down = sub[sub["log2fc"] < 0].nsmallest(5, "log2fc")
        R.append(f"**{sym}** ({len(sub)} trans targets). Up: "
                 + ", ".join(f"{g} ({v:+.2f})" for g, v in zip(up["tested_gene_symbol"], up["log2fc"]))
                 + ". Down: " + ", ".join(f"{g} ({v:+.2f})" for g, v in zip(down["tested_gene_symbol"], down["log2fc"])) + ".")
        R.append("")
    if not args.no_plots:
        p = fig_hbar(list(top_reg.items()), "Top trans regulators (significant targets)", "trans targets",
                     out_dir / "top_trans_regulators.png")
        if p:
            figures.append(p)
    if args.tf_list:
        tfs = set(read_table(args.tf_list, header=None).iloc[:, 0].astype(str))
        perturbed = set(trans_sig["element_symbol"]) if targeting_symbols is None else targeting_symbols
        summary["tf_list"] = {"n_tfs": len(tfs), "perturbed_in_list": len(perturbed & tfs),
                              "with_trans_hits_in_list": len(set(trans_sig["element_symbol"]) & tfs)}
        R.append(f"TF list `{Path(args.tf_list).name}`: {len(tfs):,} names; {len(perturbed & tfs):,} perturbed "
                 f"elements are on it, {len(set(trans_sig['element_symbol']) & tfs):,} of those have trans hits.")
        R.append("")

    # ---- 5. GWAS ------------------------------------------------------------
    gwas = index = None
    gene_coords = gene_coords_from(meta, Path(args.gene_coords) if args.gene_coords else None)
    if args.gwas:
        gwas = load_gwas(Path(args.gwas))
        index = PositionIndex(gwas)
        summary["gwas"] = {"associations": int(len(gwas)), "traits": int(gwas["DISEASE/TRAIT"].nunique())}
        R.append("## 4. GWAS Catalog overlap")
        R.append("")
        R.append(f"Catalog: {len(gwas):,} mapped associations, {gwas['DISEASE/TRAIT'].nunique():,} traits. "
                 f"Window: {args.gwas_window:,} bp either side.")
        R.append("")
        # 5a. SNPs near the top TF elements
        top_syms = list(reg_counts.head(args.top_regulators).index)
        el = parse_element_coords(trans_sig[trans_sig["element_symbol"].isin(top_syms)])
        el = el.drop_duplicates("element_symbol")
        tf_gwas = overlap_intervals(index, gwas, el, ["element_symbol"], "el_chr", "el_start", "el_end", args.gwas_window)
        tf_gwas = tf_gwas.rename(columns={"element_symbol": "tf_symbol"})
        tables.append(write_csv(tf_gwas, out_dir / "gwas_near_tf_elements.csv"))
        R.append(f"**SNPs within {args.gwas_window // 1000} kb of the top {len(el)} TF elements:** "
                 f"{len(tf_gwas):,} associations, {tf_gwas['tf_symbol'].nunique()} TFs, {tf_gwas['trait'].nunique():,} traits.")
        R.append("")
        # 5b. SNPs near trans target genes
        gene_gwas = None
        if gene_coords is not None:
            tg = gene_coords[gene_coords["gene_id"].isin(set(trans_sig["tested_gene_id"].astype(str)))
                             | gene_coords["symbol"].isin(set(trans_sig["tested_gene_symbol"].astype(str)))]
            gene_gwas = overlap_intervals(index, gwas, tg, ["gene_id", "symbol"], "chr", "start", "end", args.gwas_window)
            gene_gwas = gene_gwas.rename(columns={"symbol": "gene_symbol"})
            tables.append(write_csv(gene_gwas, out_dir / "gwas_near_target_genes.csv"))
            R.append(f"**SNPs within {args.gwas_window // 1000} kb of {len(tg):,} trans target genes:** "
                     f"{len(gene_gwas):,} associations, {gene_gwas['gene_symbol'].nunique():,} genes, "
                     f"{gene_gwas['trait'].nunique():,} traits.")
            R.append("")
        else:
            R.append("_Gene coordinates unavailable (no --mudata with gene_chr/start/end and no --gene-coords), "
                     "so the target-gene overlap was skipped._")
            R.append("")
        # 5c. Trait enrichment
        R.append("### Trait enrichment (Fisher exact, BH-adjusted)")
        R.append("")
        R.append("Counts are catalog associations, as in the notebook: several SNPs in one locus, or one SNP in "
                 "several studies, each count once, so traits mapped to few loci show inflated odds ratios. "
                 "Read `observed` alongside `odds_ratio`.")
        R.append("")
        enr_tf = trait_enrichment(tf_gwas, gwas, args.min_hits)
        tables.append(write_csv(enr_tf, out_dir / "trait_enrichment_near_tf_elements.csv"))
        summary["enrichment_top_tf"] = enr_tf.head(5).to_dict("records")
        R.append(f"**Near TF elements** (traits with >= {args.min_hits} hits):")
        R.append("")
        R.append(md_table(["trait", "observed", "catalog", "expected", "enrichment", "odds", "FDR"],
                          enr_tf.head(20)[["trait", "observed", "catalog_total", "expected", "enrichment", "odds_ratio", "fdr"]].values.tolist()))
        R.append("")
        if gene_gwas is not None:
            enr_gene = trait_enrichment(gene_gwas, gwas, args.min_hits)
            tables.append(write_csv(enr_gene, out_dir / "trait_enrichment_near_target_genes.csv"))
            summary["enrichment_top_gene"] = enr_gene.head(5).to_dict("records")
            R.append("**Near trans target genes:**")
            R.append("")
            R.append(md_table(["trait", "observed", "catalog", "expected", "enrichment", "odds", "FDR"],
                              enr_gene.head(20)[["trait", "observed", "catalog_total", "expected", "enrichment", "odds_ratio", "fdr"]].values.tolist()))
            R.append("")
            if not args.no_plots and len(enr_gene):
                p = fig_hbar([(t, e) for t, e in zip(enr_gene.head(15)["trait"], enr_gene.head(15)["enrichment"])],
                             "GWAS traits enriched near trans target genes (top 15 by FDR)", "observed / expected",
                             out_dir / "trait_enrichment_target_genes.png", color=ORANGE)
                if p:
                    figures.append(p)

    # ---- 6. E2G -------------------------------------------------------------
    e2g = load_e2g(args.e2g, args.e2g_score_min) if args.e2g else None
    e2g_gwas = None
    if e2g is not None:
        R.append("## 5. scE2G enhancer-gene predictions")
        R.append("")
        per_src = e2g.groupby("source").agg(links=("Score", "size"), genes=("GeneSymbol", "nunique"),
                                            elements=("ElementName", "nunique")).reset_index()
        R.append(md_table(["model", "links >= threshold", "genes", "elements"], per_src.values.tolist()))
        R.append("")
        srcs = list(dict.fromkeys(e2g["source"]))
        if len(srcs) >= 2:
            shared = set.intersection(*(set(e2g.loc[e2g["source"] == s, "GeneSymbol"]) for s in srcs))
            R.append(f"Genes with links in every model: {len(shared):,}.")
            R.append("")
        e2g_genes = set(e2g["GeneSymbol"].astype(str))
        tfs = set(trans_sig["element_symbol"].astype(str))
        targets = set(trans_sig["tested_gene_symbol"].astype(str))
        tf_in, tg_in = tfs & e2g_genes, targets & e2g_genes
        a = trans_sig[trans_sig["element_symbol"].isin(tf_in)].copy()
        a["e2g_match"] = "perturbed_TF"
        b = trans_sig[trans_sig["tested_gene_symbol"].isin(tg_in)].copy()
        b["e2g_match"] = "target_gene"
        trans_e2g = pd.concat([a, b], ignore_index=True).drop_duplicates(["element_id", "tested_gene_id", "e2g_match"])
        gene_source = e2g.groupby("GeneSymbol")["source"].apply(lambda x: ",".join(sorted(set(x)))).to_dict()
        trans_e2g["e2g_source"] = [gene_source.get(r.element_symbol if r.e2g_match == "perturbed_TF" else r.tested_gene_symbol, "")
                                   for r in trans_e2g.itertuples()]
        tables.append(write_csv(trans_e2g, out_dir / "trans_sig_e2g.csv"))
        summary["e2g"] = {"perturbed_tfs_with_e2g": len(tf_in), "perturbed_tfs": len(tfs),
                          "targets_with_e2g": len(tg_in), "targets": len(targets),
                          "trans_rows_with_e2g": int(len(trans_e2g))}
        R.append(f"Perturbed TFs with E2G links: {len(tf_in)} / {len(tfs)}; trans target genes with E2G links: "
                 f"{len(tg_in)} / {len(targets)}. Trans rows touching an E2G gene: {len(trans_e2g):,}.")
        R.append("")
        R.append(md_table(["E2G model(s)", "trans rows"], trans_e2g["e2g_source"].value_counts().items()))
        R.append("")
        if not args.no_plots:
            figures += fig_e2g(e2g, args.e2g_score_min, out_dir)
        # 6b. GWAS SNPs inside E2G elements of those genes
        if gwas is not None:
            genes_in = tf_in | tg_in
            els = e2g[e2g["GeneSymbol"].isin(genes_in)][["chr_clean", "ElementStart", "ElementEnd", "ElementName",
                                                          "GeneSymbol", "Score", "source"]].drop_duplicates()
            e2g_gwas = overlap_intervals(index, gwas, els, ["ElementName", "GeneSymbol", "Score", "source"],
                                         "chr_clean", "ElementStart", "ElementEnd", 0)
            e2g_gwas = e2g_gwas.rename(columns={"ElementName": "element_name", "GeneSymbol": "e2g_gene",
                                                "Score": "e2g_score", "source": "e2g_source", "pvalue": "gwas_pvalue"})
            e2g_gwas["gene_matches_gwas"] = e2g_gwas["e2g_gene"].astype(str) == e2g_gwas["mapped_gene"].astype(str)
            tables.append(write_csv(e2g_gwas, out_dir / "e2g_elements_gwas_overlap.csv"))
            strong = e2g_gwas[(e2g_gwas["e2g_score"] > 0.5) & (e2g_gwas["gwas_pvalue"] < 5e-8)]
            novel = e2g_gwas[~e2g_gwas["gene_matches_gwas"]]
            summary["e2g_gwas"] = {"overlaps": int(len(e2g_gwas)), "elements": int(e2g_gwas["element_name"].nunique()),
                                   "genes": int(e2g_gwas["e2g_gene"].nunique()), "snps": int(e2g_gwas["snp"].nunique()),
                                   "traits": int(e2g_gwas["trait"].nunique()),
                                   "frac_gene_matches_mapped": float(e2g_gwas["gene_matches_gwas"].mean()) if len(e2g_gwas) else float("nan"),
                                   "novel_assignments": int(len(novel)), "strong_hits": int(len(strong))}
            R.append("### GWAS SNPs inside E2G elements")
            R.append("")
            R.append(f"{len(e2g_gwas):,} SNP-element overlaps across {e2g_gwas['element_name'].nunique():,} elements, "
                     f"{e2g_gwas['e2g_gene'].nunique():,} genes, {e2g_gwas['snp'].nunique():,} SNPs and "
                     f"{e2g_gwas['trait'].nunique():,} traits (E2G elements of {len(els):,} linked to {len(genes_in):,} genes).")
            if len(e2g_gwas):
                R.append(f"E2G gene equals the catalog's mapped gene in {e2g_gwas['gene_matches_gwas'].mean():.1%} of overlaps; "
                         f"{len(novel):,} overlaps ({novel['e2g_gene'].nunique():,} genes) point at a different gene than "
                         f"proximity. Strong hits (E2G > 0.5 and GWAS p < 5e-8): {len(strong):,} across "
                         f"{strong['e2g_gene'].nunique():,} genes and {strong['trait'].nunique():,} traits.")
            R.append("")
            if len(e2g_gwas):
                R.append("Genes with most GWAS traits via their E2G elements:")
                R.append("")
                gt = e2g_gwas.groupby("e2g_gene")["trait"].nunique().sort_values(ascending=False).head(20)
                R.append(md_table(["gene", "traits"], list(gt.items())))
                R.append("")
                R.append("Top traits in E2G elements:")
                R.append("")
                R.append(md_table(["trait", "overlaps"], list(e2g_gwas["trait"].value_counts().head(20).items())))
                R.append("")
                if not args.no_plots:
                    p = fig_hbar(list(gt.items()), "Genes with most GWAS traits in their E2G elements", "distinct traits",
                                 out_dir / "e2g_gwas_genes.png")
                    if p:
                        figures.append(p)

    # ---- 7. ChIP support ----------------------------------------------------
    if args.bigwig and e2g_gwas is not None and len(e2g_gwas):
        loci = e2g_gwas[["snp", "pos"]].copy()
        chrom_of = gwas.drop_duplicates("SNPS").set_index("SNPS")["CHR_ID"]
        loci["chrom"] = "chr" + loci["snp"].map(chrom_of).astype(str)
        loci = loci.drop_duplicates(["chrom", "pos"]).reset_index(drop=True)
        try:
            chip, chip_summary = bigwig_support(loci, args.bigwig, args.signal_frac)
        except RuntimeError as exc:
            R.append(f"_ChIP step skipped: {exc}_")
            R.append("")
            chip = None
        if chip is not None:
            merged = e2g_gwas.merge(chip.drop(columns=["snp"]), left_on=["pos"], right_on=["pos"], how="left")
            tables.append(write_csv(merged, out_dir / "e2g_gwas_chipseq.csv"))
            supported = merged[merged["n_chipseq_passing"] > 0]
            tables.append(write_csv(supported, out_dir / "e2g_gwas_chipseq_supported.csv"))
            summary["chip"] = {"tracks": chip_summary, "positions": int(len(loci)),
                               "positions_with_signal": int((chip["n_chipseq_passing"] > 0).sum()),
                               "rows_supported": int(len(supported)), "rows": int(len(merged))}
            R.append("## 6. ChIP-seq support at GWAS positions")
            R.append("")
            R.append(f"{len(loci):,} distinct positions queried; signal called present at >= {args.signal_frac:.0%} of each "
                     f"track's maximum. {int((chip['n_chipseq_passing'] > 0).sum()):,} positions have signal in at least one "
                     f"track; {len(supported):,} / {len(merged):,} E2G-GWAS rows are supported.")
            R.append("")
            R.append(md_table(["track", "max", "threshold", "positions passing"],
                              [[s["track"], s["max"], s["threshold"], f"{s['passing']:,} / {s['n']:,}"] for s in chip_summary]))
            R.append("")

    # ---- report footer -----------------------------------------------------
    n = R and sum(1 for ln in R if ln.startswith("## ")) + 1
    R.append(f"## {n}. Files")
    R.append("")
    for f in figures:
        R.append(f"- `{f.name}`  ![{f.stem}]({f.name})")
    for t in tables:
        R.append(f"- `{t.name}`")
    R.append("")
    R.append("## Method notes")
    R.append("")
    R.append("- Port of IGVF/tf_perturb_seq `docs/jamborees/2026_UTSW/working_groups/wg3_disease_gwas/perturb_seq_analysis_v4.ipynb`.")
    R.append("- Positional overlaps use sorted per-chromosome positions (identical results to the notebook's loops).")
    R.append("- Trans significance filtering is chunked; the filtered rows are cached under `Data/Cache/TFPerturbSeq` keyed on file and thresholds"
             + (" (this run used the cache)." if cached else "."))
    R.append("- BH adjustment is the standard step-up procedure; enrichment counts associations, not loci.")
    R.append(f"- Inputs: calibrated prefix `{prefix}`" + (f", mudata `{args.mudata}`" if args.mudata else "")
             + (f", GWAS `{args.gwas}`" if args.gwas else "") + (f", E2G {args.e2g}" if args.e2g else "")
             + (f", bigWig {args.bigwig}" if args.bigwig else "") + ".")
    report = out_dir / "report.md"
    report.write_text("\n".join(R) + "\n", encoding="utf-8")
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    c = summary["counts"]
    print(f"Direct targets significant: {c['direct_sig']:,} / {c['direct_total']:,}")
    print(f"Cis significant: {c['cis_sig']:,} / {c['cis_total']:,}")
    print(f"Trans significant (targeting): {c['trans_sig']:,} rows, {c['trans_regulators']:,} regulators, "
          f"{c['trans_target_genes']:,} target genes")
    print("Top trans regulators: " + ", ".join(f"{s} ({n})" for s, n in list(top_reg.items())[:8]))
    if "gwas" in summary:
        top = summary.get("enrichment_top_gene") or summary.get("enrichment_top_tf") or []
        if top:
            print("Top GWAS trait enrichment: " + "; ".join(f"{t['trait']} (obs {t['observed']}, FDR {t['fdr']:.1e})" for t in top[:3]))
    if "e2g_gwas" in summary:
        eg = summary["e2g_gwas"]
        print(f"E2G x GWAS: {eg['overlaps']:,} overlaps, {eg['genes']:,} genes, {eg['strong_hits']:,} strong hits")
    if "chip" in summary:
        print(f"ChIP support: {summary['chip']['rows_supported']:,} / {summary['chip']['rows']:,} rows")
    print(f"Report: {report}")
    for t in tables[:6]:
        print(f"CSV: {t}")
    if len(tables) > 6:
        print(f"CSV: ... {len(tables) - 6} more in {out_dir}")
    for f in figures:
        print(f"Figure: {f}")
    print(f"JSON: {out_dir / 'summary.json'}")
    print(f"Log: {log_path}")
    return 0


# ---------------------------------------------------------------------------
# Self-test on synthetic data with planted signals
# ---------------------------------------------------------------------------

def make_synthetic(dirpath: Path, seed: int = 7) -> dict:
    """Small inputs in the upstream formats, with signals to recover.

    Planted: FOXH1 has 60 strong trans targets and SOX17 40; a 'Type 2
    diabetes' block of GWAS SNPs sits inside FOXH1's target genes; one E2G
    element for a FOXH1 target contains one of those SNPs; a non-targeting
    element carries 'significant' rows that must be dropped.
    """
    import numpy as np
    pd = _pd()
    rng = np.random.default_rng(seed)
    n_genes, n_tfs = 300, 40
    genes = [f"ENSG{100000 + i:08d}" for i in range(n_genes)]
    symbols = [f"G{i}" for i in range(n_genes)]
    chroms = ["1", "2", "3"]
    gene_chr = [chroms[i % 3] for i in range(n_genes)]
    gene_start = [1_000_000 + (i // 3) * 200_000 for i in range(n_genes)]
    gene_end = [s + 20_000 for s in gene_start]
    coords = pd.DataFrame({"gene_id": genes, "symbol": symbols, "chr": gene_chr,
                           "start": gene_start, "end": gene_end})
    coords.to_csv(dirpath / "genes.tsv", sep="\t", index=False)

    tf_syms = ["FOXH1", "SOX17"] + [f"TF{i}" for i in range(n_tfs - 2)]
    tf_ids = [f"ENSG{200000 + i:08d}" for i in range(n_tfs)]
    el_ids = [f"{tid}|chr{chroms[i % 3]}:{5_000_000 + i * 50_000}-{5_000_000 + i * 50_000 + 400}"
              for i, tid in enumerate(tf_ids)]

    def row(el, sym, label, gid, gsym, lfc, se, padj, cis=False, direct=False, n=800):
        return {"element_id": el, "element_symbol": sym, "element_label": label, "tested_gene_id": gid,
                "tested_gene_symbol": gsym, "n_cells": n, "log2fc": lfc, "log2fc_se": se, "is_cis": cis,
                "is_direct_target": direct, "posterior_pval": padj / 10, "empirical_pval": padj / 5,
                "empirical_pval_adj": padj}

    direct_rows, cis_rows, trans_rows = [], [], []
    planted = {"FOXH1": set(range(0, 60)), "SOX17": set(range(100, 140))}
    for i, (el, sym, tid) in enumerate(zip(el_ids, tf_syms, tf_ids)):
        strong = i < 30
        direct_rows.append(row(el, sym, "targeting", tid, sym, -1.8 if strong else -0.05, 0.2,
                               1e-8 if strong else 0.6, cis=True, direct=True))
        cis_rows.append(direct_rows[-1])
        cis_rows.append(row(el, sym, "targeting", genes[i], symbols[i], 0.4 if i < 10 else 0.02, 0.1,
                            1e-3 if i < 10 else 0.9, cis=True))
        for g in range(n_genes):
            if g in planted.get(sym, set()):
                lfc = float(rng.choice([-1, 1]) * rng.uniform(1.5, 3.0))
                padj = 10 ** -rng.uniform(6, 30)
            else:
                lfc = float(rng.normal(0, 0.15))
                padj = float(rng.uniform(0.06, 1.0))
            trans_rows.append(row(el, sym, "targeting", genes[g], symbols[g], lfc, 0.1, padj))
    # A non-targeting element with rows that pass the cutoffs: must be excluded.
    for g in range(10):
        trans_rows.append(row("NT|chr1:1-2", "non-targeting", "non-targeting", genes[g], symbols[g], 2.5, 0.1, 1e-9))
    prefix = str(dirpath / "synthetic_calibrated_")
    pd.DataFrame(direct_rows).to_csv(prefix + "direct_target_results.tsv", sep="\t", index=False)
    pd.DataFrame(cis_rows).to_csv(prefix + "cis_results.tsv", sep="\t", index=False)
    pd.DataFrame(trans_rows).to_csv(prefix + "trans_results.tsv", sep="\t", index=False)

    # GWAS catalog: planted T2D SNPs inside FOXH1's targets, background elsewhere.
    gw = []
    for k, g in enumerate(sorted(planted["FOXH1"])[:25]):
        gw.append({"CHR_ID": gene_chr[g], "CHR_POS": gene_start[g] + 5_000 + k, "DISEASE/TRAIT": "Type 2 diabetes",
                   "SNPS": f"rs{1000 + k}", "P-VALUE": 1e-12, "MAPPED_GENE": symbols[g]})
    for k in range(600):
        chrom = chroms[k % 3]
        gw.append({"CHR_ID": chrom, "CHR_POS": int(rng.integers(30_000_000, 60_000_000)),
                   "DISEASE/TRAIT": f"Trait {k % 40}", "SNPS": f"rs{50000 + k}", "P-VALUE": 1e-9, "MAPPED_GENE": "X"})
    # One SNP inside an E2G element for FOXH1 target G0, mapped by proximity to another gene.
    gw.append({"CHR_ID": gene_chr[0], "CHR_POS": gene_start[0] - 30_000 + 150, "DISEASE/TRAIT": "Type 2 diabetes",
               "SNPS": "rs777", "P-VALUE": 1e-20, "MAPPED_GENE": "OTHER"})
    pd.DataFrame(gw).to_csv(dirpath / "gwas.tsv", sep="\t", index=False)

    # E2G: two models; an enhancer 30 kb upstream of G0 -> G0 (contains rs777), promoters for a few genes.
    e2g_rows = []
    for src in ("ESC", "DE"):
        e2g_rows.append({"ElementChr": f"chr{gene_chr[0]}", "ElementStart": gene_start[0] - 30_000,
                         "ElementEnd": gene_start[0] - 29_500, "ElementName": f"chr{gene_chr[0]}:e0",
                         "GeneSymbol": symbols[0], "Score": 0.9, "ElementClass": "intergenic", "isSelfPromoter": False,
                         "source_note": src})
        for g in range(1, 50):
            e2g_rows.append({"ElementChr": f"chr{gene_chr[g]}", "ElementStart": gene_start[g] - 500,
                             "ElementEnd": gene_start[g] + 500, "ElementName": f"chr{gene_chr[g]}:p{g}",
                             "GeneSymbol": symbols[g], "Score": float(rng.uniform(0.05, 1.0)),
                             "ElementClass": "promoter", "isSelfPromoter": True, "source_note": src})
        e2g_rows.append({"ElementChr": "chr1", "ElementStart": 10, "ElementEnd": 20, "ElementName": "chr1:tf",
                         "GeneSymbol": "FOXH1", "Score": 0.8, "ElementClass": "genic", "isSelfPromoter": False,
                         "source_note": src})
        df = pd.DataFrame([r for r in e2g_rows if r["source_note"] == src]).drop(columns=["source_note"])
        with open(dirpath / f"{src.lower()}.e2g.tsv", "w") as fh:
            fh.write("# synthetic scE2G\n")
            df.to_csv(fh, sep="\t", index=False)
    mudata_path = None
    try:
        import h5py  # type: ignore
        import scipy.sparse as sp  # type: ignore
        from anndata.experimental import write_elem  # type: ignore
        n_cells = 500
        gene_var = pd.DataFrame({"symbol": symbols, "gene_chr": ["chr" + c for c in gene_chr],
                                 "gene_start": gene_start, "gene_end": gene_end}, index=genes)
        guide_var = pd.DataFrame({
            "guide_id": [f"g{i}" for i in range(n_tfs * 2 + 4)],
            "targeting": [True] * (n_tfs * 2) + [False] * 4,
            "type": ["targeting"] * (n_tfs * 2) + ["non-targeting"] * 4,
            "intended_target_name": [tf_ids[i // 2] for i in range(n_tfs * 2)] + ["non-targeting"] * 4,
            "gene_name": [tf_syms[i // 2] for i in range(n_tfs * 2)] + ["non-targeting"] * 4,
        }).set_index("guide_id")
        assign = sp.random(n_cells, len(guide_var), density=0.05, format="csr", random_state=1)
        assign.data[:] = 1
        obs = pd.DataFrame({"batch": [f"IGVFDS{i % 2}" for i in range(n_cells)]},
                           index=[f"cell{i}" for i in range(n_cells)])
        mudata_path = dirpath / "inference_mudata.h5mu"
        with h5py.File(str(mudata_path), "w") as h:
            write_elem(h, "obs", obs)
            write_elem(h, "mod/gene/var", gene_var)
            write_elem(h, "mod/guide/var", guide_var)
            write_elem(h, "mod/guide/layers/guide_assignment", assign)
    except Exception as exc:  # pragma: no cover - h5py/anndata absent
        logging.warning("synthetic h5mu not written: %s", exc)
    return {"prefix": prefix, "gwas": str(dirpath / "gwas.tsv"), "genes": str(dirpath / "genes.tsv"),
            "mudata": str(mudata_path) if mudata_path else None,
            "e2g": [f"ESC={dirpath / 'esc.e2g.tsv'}", f"DE={dirpath / 'de.e2g.tsv'}"],
            "planted": {k: len(v) for k, v in planted.items()}}


def selftest(args: argparse.Namespace) -> int:
    pd = _pd()
    with tempfile.TemporaryDirectory() as td:
        syn = make_synthetic(Path(td))
        ns = argparse.Namespace(calibrated_prefix=syn["prefix"], mudata=syn["mudata"], gwas=syn["gwas"],
                                e2g=syn["e2g"], bigwig=[], gene_coords=None if syn["mudata"] else syn["genes"],
                                tf_list=None,
                                label="selftest", padj=0.05, lfc_cis=0.2, lfc_trans=1.0, top_regulators=20,
                                gwas_window=50_000, e2g_score_min=0.177, min_hits=5, signal_frac=0.10,
                                no_plots=args.no_plots, no_cache=True)
        rc = run(ns)
        if rc != 0:
            print("selftest: run() failed", file=sys.stderr)
            return 1
        out_dir = sorted(OUT_ROOT.glob("*_selftest"))[-1]
        s = json.loads((out_dir / "summary.json").read_text())
        checks = []
        c = s["counts"]
        checks.append((c["trans_sig"] == 100, f"trans_sig == 100 planted (got {c['trans_sig']}; non-targeting rows excluded)"))
        checks.append((c["direct_sig"] == 30, f"direct_sig == 30 (got {c['direct_sig']})"))
        regs = pd.read_csv(out_dir / "trans_regulators.csv")
        checks.append((list(regs["element_symbol"][:2]) == ["FOXH1", "SOX17"], f"top regulators FOXH1, SOX17 (got {list(regs['element_symbol'][:2])})"))
        enr = pd.read_csv(out_dir / "trait_enrichment_near_target_genes.csv")
        checks.append((len(enr) and enr.iloc[0]["trait"] == "Type 2 diabetes" and enr.iloc[0]["fdr"] < 1e-6,
                       f"planted trait is the top enrichment (got {enr.iloc[0]['trait'] if len(enr) else None})"))
        eg = pd.read_csv(out_dir / "e2g_elements_gwas_overlap.csv")
        checks.append((bool((eg["snp"] == "rs777").any()) and bool((~eg.loc[eg["snp"] == "rs777", "gene_matches_gwas"]).all()),
                       "rs777 found inside the G0 enhancer and flagged as a non-nearest-gene assignment"))
        te = pd.read_csv(out_dir / "trans_sig_e2g.csv")
        checks.append(("perturbed_TF" in set(te["e2g_match"]) and "target_gene" in set(te["e2g_match"]),
                       "E2G intersection found both perturbed-TF and target-gene rows"))
        # Vectorised overlap agrees with brute force on the synthetic catalog.
        import numpy as np
        gwas = load_gwas(Path(syn["gwas"]))
        idx = PositionIndex(gwas)
        coords = read_table(syn["genes"])
        coords["chr"] = coords["chr"].astype(str)
        fast = overlap_intervals(idx, gwas, coords, ["symbol"], "chr", "start", "end", 50_000)
        brute = 0
        for _, g in coords.iterrows():
            brute += int(((gwas["CHR_ID"] == g["chr"]) & (gwas["CHR_POS"] >= g["start"] - 50_000)
                          & (gwas["CHR_POS"] <= g["end"] + 50_000)).sum())
        checks.append((len(fast) == brute, f"vectorised overlap == brute force ({len(fast)} vs {brute})"))
        checks.append((np.allclose(bh_adjust([0.01, 0.04, 0.03, 0.2]), [0.04, 0.0533333, 0.0533333, 0.2]),
                       "BH adjustment matches the step-up procedure"))
        if syn["mudata"]:
            m = s.get("mudata", {})
            checks.append((m.get("n_cells") == 500 and m.get("n_targets") == 41 and m.get("cells_with_guide", 0) > 0
                           and set(m.get("batches", {})) == {"IGVFDS0", "IGVFDS1"},
                           f"h5mu read without the mudata package (cells {m.get('n_cells')}, targets {m.get('n_targets')})"))
            checks.append((c["trans_sig"] == 100, "targeting filter from the guide table keeps the planted rows"))
        ok = True
        for passed, msg in checks:
            print(("  ok    " if passed else "  FAIL  ") + msg)
            ok &= bool(passed)
        if not args.keep:
            import shutil
            shutil.rmtree(out_dir, ignore_errors=True)
        print("selftest: all checks pass" if ok else "selftest: FAILED")
        return 0 if ok else 1


def main(argv: "Optional[list[str]]" = None) -> int:
    p = argparse.ArgumentParser(prog="igvfagent tf-perturb",
                                description="IGVF TF Perturb-seq calibrated results -> disease / GWAS overlay "
                                            "(port of the tf_perturb_seq WG3 notebook).")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="Run the analysis on calibrated result tables.")
    r.add_argument("--calibrated-prefix", required=True,
                   help="Path prefix of the calibrated TSVs, e.g. data/<run>_calibrated_ "
                        "(expects <prefix>direct_target_results.tsv, cis_results.tsv, trans_results.tsv).")
    r.add_argument("--mudata", help="inference_mudata.h5mu for the guide-library overview and gene coordinates.")
    r.add_argument("--gene-coords", help="TSV with gene_id, symbol, chr, start, end (used when no mudata).")
    r.add_argument("--gwas", help="GWAS Catalog associations TSV (full download).")
    r.add_argument("--e2g", action="append", default=[], metavar="LABEL=PATH",
                   help="scE2G predictions TSV; repeat for several models (ESC=..., DE=...).")
    r.add_argument("--bigwig", action="append", default=[], help="ChIP-seq bigWig; repeatable (needs pyBigWig).")
    r.add_argument("--tf-list", help="Optional text file of TF names, one per line.")
    r.add_argument("--label", default="")
    r.add_argument("--padj", type=float, default=0.05)
    r.add_argument("--lfc-cis", type=float, default=0.2)
    r.add_argument("--lfc-trans", type=float, default=1.0)
    r.add_argument("--top-regulators", type=int, default=20)
    r.add_argument("--gwas-window", type=int, default=50_000)
    r.add_argument("--e2g-score-min", type=float, default=0.177)
    r.add_argument("--min-hits", type=int, default=5, help="Minimum hits for a trait to be tested.")
    r.add_argument("--signal-frac", type=float, default=0.10, help="ChIP pass = signal >= frac * track max.")
    r.add_argument("--no-plots", action="store_true")
    r.add_argument("--no-cache", action="store_true", help="Re-filter the trans table even if cached.")
    r.set_defaults(func=run)
    s = sub.add_parser("selftest", help="Synthetic inputs with planted signals; checks they are recovered.")
    s.add_argument("--no-plots", action="store_true")
    s.add_argument("--keep", action="store_true", help="Keep the selftest output directory.")
    s.set_defaults(func=selftest)
    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
