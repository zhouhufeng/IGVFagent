#!/usr/bin/env python3
"""Pooled CRISPR screen toolkit: screen object, normalisation, LFC, sorting-screen stats, QC, guide design (port of IGVF-CRISPR/perturb-tools).

Port of https://github.com/IGVF-CRISPR/perturb-tools (MIT License, Copyright (c)
2021 Michael Vinyard; authors M. E. Vinyard and J. K. Ryu, Pinello lab).  The
upstream is a small Python package ("pt") for bulk pooled CRISPR screens built
on AnnData: a screen is samples x guides (`X`), with `obs` (sample annotation),
`var` (guide annotation), `layers`, `obsm`, `obsp` and `uns`.  Every module of
the default branch (main @ 81eaae0) was read and re-derived here in
numpy / pandas / scipy; nothing was vendored.  Two files that main's
`pt.io.__init__` imports but that are missing on main (so `import perturb_tools`
fails at that commit) were taken from the `fix_import` branch (d0148ce):
`_to_from_csv.py` (pt.io.read_csv / write_csv) and `_export_as_input.py`
(pt.io.to_mageck_input).  The guide-library design code (`_experimental_design`),
the replicate-consistency / guide-enrichment plots and the single-cell
`sccrispr_tools.calculate_distance` that exist only on the `dev` branch
(5fc786b) are ported too.  The three bulk tutorials (notebooks/bulk/
basic_api_demo, anndata_demo, sample_quality_report) are the `run` pipeline.
Relationship: port (algorithms translated, source consulted, code rewritten).

Definitions reproduced from the upstream
  screen layout     X is samples x guides; obs index = sample names, var index =
                    guide ids.  pt.io.read_csv reads a guides x samples count
                    table (first column = guide id) and transposes it.
  TKO sample names  (notebooks) replicate = last character, time = int(name[1:-1])
                    or -1 when empty ("T0" -> replicate "0", time -1).
  add               X1 + X2 when obs and var indices are identical, else error.
  log_norm          layer "lognorm_counts" = log2(X / rowsum(X) * 1e6 + 1)
                    (reads per million per sample, then log2(x + 1)); with a
                    read-count layer L the output layer is "lognorm_L".  This is
                    PoolQ's lognormalized-counts formula (checked against PoolQ
                    3.3.2 test output to 4e-15).
  log_fold_change   var["{s1}_{s2}.lfc"] = lognorm[s1] - lognorm[s2].
  fold_change       var["{c1}_{c2}.fc"], same arithmetic.
  lfc reps          per replicate r (obs[rep_col]) the unique sample with
                    obs[compare_col] == cond1 minus the one with cond2 ->
                    columns "{r}.{cond1}_{cond2}.lfc"; aggregate over replicates
                    with mean, median or sd (np.std, ddof = 0) ->
                    "{cond1}_{cond2}.lfc.{fn}".
  sorting-screen LFC (pt.calc.log_fold_change / LogFoldChangeModule)
                    columns of the log-normalised guides x samples table matched
                    by regex (pandas .filter(regex=...)) for condit_1, condit_2,
                    control; guides filtered to var.target in targets; baseline
                    subtraction c - control (element-wise, per replicate column);
                    lfc.mean = mean_r(c1 - ctrl) - mean_r(c2 - ctrl); lfc.stdev =
                    std_r((c1 - ctrl) - (c2 - ctrl)) (ddof = 0); lfc.p_value =
                    two-sided equal-variance t-test (scipy ttest_ind) of the
                    replicate baseline-LFCs of c1 vs c2 per guide;
                    "-log10(pval)".
  sample correlation  pearson or spearman between samples over guides ->
                    obsm["{prefix}corr_{layer}"], obs["{prefix}mean_corr_{layer}"]
                    (column mean including the diagonal) and
                    obs["{prefix}median_corr_{layer}"].
  Gini (_G)         bins = linspace(0, 1, 50); y(b) = sum(v[v <= quantile(v, b)]) /
                    sum(v); gini = (A_eq - A_lorenz) / A_eq with trapezoid areas
                    (A_eq = 0.5).  Ties at the minimum count towards y(0), so a
                    constant vector gives -1 (upstream behaviour, kept).
  count dist        per-sample histogram on 10 ** linspace(0, log10(max + 1), 100)
                    bins (log x) with the per-sample median.
  LFC correlation   lfc reps on the (optionally positive-control) guide subset,
                    then DataFrame.corr(method) between replicate LFC columns;
                    replicate consistency = lower-triangle scatter matrix with
                    Pearson r (NaN pairs dropped).
  outlier guides    RPM layer "{layer}_RPM" = X / rowsum * 1e6; within each
                    obs[cond_col] condition a guide is an outlier in sample i when
                    RPM_i > mad_z_thres * median_condition(RPM of that guide) AND
                    RPM_i > abs_RPM_thres (defaults 5 and 10000).  Output columns:
                    <var index name>, sample, RPM.
  MAGeCK input      columns sgRNA, gene, <sample_prefix + sample...>; counts
                    NaN -> 0, int; spaces in sgRNA -> "_"; rows with a null or
                    empty gene dropped; tab-separated.
  Excel export      sheets X (guides x samples), each layer, guides, samples,
                    obsm.*, obsp.*, optional uns.*; the experiment report appends
                    extra tables as sheets.
  PoolQ reader      counts.txt (fallback expected-counts.txt) -> X, var
                    [barcode, barcode_id] from "Row Barcode"/"Row Barcode IDs"
                    (or PoolQ 2 "Construct Barcode"/"Construct IDs");
                    lognormalized-counts.txt -> layer lognorm_counts;
                    barcode-counts.txt, unexpected-sequences.txt -> uns tables;
                    correlation.txt -> obsp; quality.txt parsed line-wise
                    (":" in field 1 -> metadata [metric, statistic]; > 2 fields ->
                    bc_readcounts, first such line is the header; 2 fields ->
                    common_barcodes; "Read counts..." lines skipped); runinfo.txt.
  guide annotation  protospacer = barcode when every barcode is 20 nt; target =
                    direct barcode_id -> target map, else the first annotation
                    string contained in barcode_id, else barcode_id; position =
                    first regex match of the protospacer in the region's + strand
                    sequence (Start = region start + span0, End = + span1), else in
                    its reverse complement (Start = region end - span0, End =
                    region end - span1, i.e. Start is the 5' end on "-");
                    center = (Start + End) / 2.
  guide library     (dev branch) PAM (default NGG, N = any base) found with
                    overlapping matches on both strands; forward PAM_loci = first
                    PAM base, reverse PAM_loci = L - (start in the reverse
                    complement); a PAM belongs to region (Start, End] (pd.cut
                    right-closed bins); + guides: protospacer = seq[p-20:p],
                    PAM = seq[p:p+3]; - guides: protospacer = rc(seq[p:p+20]),
                    PAM = rc(seq[p-3:p]); region_center = round(mean(left, right)),
                    PAM_distance_to_center = round(|p - center|); guide_context =
                    seq[p-30:p+13] (+) or seq[p-13:p+30] (-), flank 10; BsmbI flag =
                    "TCTC" or "GAGAC" in the context or its reverse complement,
                    flagged guides filtered by default; duplicates by protospacer
                    dropped; exons from the GTF, exact duplicates and duplicate
                    end/start coordinates removed, fragments that totally engulf
                    another dropped, partial overlaps reported; optional widening.
  guide biochemistry  polyT/polyG/polyC/polyA = contains TTTT/GGGG/CCCC/AAAA;
                    GC_content = (#G + #C) * 100 / 20.
  pair distance     (dev sccrispr_tools) same-chromosome pairs between two
                    feature tables; distance = |start1 - start2| (what upstream
                    computes, although its docstring names "midpoint").

Deliberate deviations (documented, upstream behaviour available where it runs)
  * get_outlier_guides indexes the samples x guides RPM matrix as if it were
    guides x samples (median over axis 1, column i = "sample i"), which mixes
    guides and samples.  The port applies the documented rule;
    --upstream-compat reproduces the transposed indexing.
  * _scan_for_BsmbI tests `seq.find(motif)` for truthiness, so a motif at
    position 0 is missed; corrected, --upstream-compat reproduces it.
  * to_mageck_input does not transpose a non-X count layer (shape error) and
    _write_to_csv flattens X with ndarray.tofile; both fixed (the CSV matrix is
    the guides x samples table pt.io.read_csv reads back).
  * _add_guide_target_metadata only assigns targets when a DirectPairDict is
    given and hard-codes "3 annotations"; the port applies the intended rule.
  * Exon selection uses an exact gene_name match (upstream: prefix match on
    the 11th attribute token, GENCODE-specific).  Sequences are upper-cased.
    Region coordinates are bp by default (upstream metadata used Mb: --units mb).
  * Sheet names are "obsm.<key>" (upstream literally writes "df_dict.<key>")
    and truncated to Excel's 31 characters.
  * Upstream plots are seaborn / plotly; here matplotlib PNGs.

Subcommands
  make-screen      counts table (+ guide / sample tables, TKO name parsing) ->
                   screen (CSV dir + .h5ad when anndata is installed)
  add-screens      add technical replicate screens (pt.add)
  log-norm         RPM + log2(x + 1) layer (pt.pp.log_norm)
  lfc              one sample vs another (pt.pp.log_fold_change / fold_change)
  lfc-reps         per-replicate LFCs + mean / median / sd aggregate
  sort-lfc         sorting-screen delta-LFC vs control with t-test p-values and
                   the guide-enrichment position plot (pt.calc.log_fold_change)
  qc               sample quality report: count distribution, Gini, count
                   correlation, LFC correlation (all / positive controls),
                   replicate consistency, outlier guides
  outlier-guides   pt.qc.get_outlier_guides
  to-mageck        MAGeCK count table (pt.io.to_mageck_input)
  to-excel         workbook of every screen table (+ experiment report sheets)
  to-csv           screen CSV directory (pt.io.write_csv); read back by --screen
  read-poolq       PoolQ output directory -> screen (pt.io PoolQ reader)
  annotate-guides  protospacer, target and genomic position of each guide
  design-library   PAM scan of exons / regions -> sgRNA library (dev branch)
  guide-biochem    homopolymer + GC content of protospacers
  pair-distance    same-chromosome distances between two feature tables
  run              the bulk tutorial pipeline end to end
  selftest         synthetic screens with planted effects; asserts every step

Output: Docs/PerturbTools/<timestamp>_<label>/ with report.md and summary.json.
numpy, pandas, scipy required; anndata (h5ad), openpyxl (Excel), matplotlib
(figures) and pysam (indexed FASTA) optional.

Usage:
    igvfagent perturb-tools make-screen --counts readcount-HeLa-lib1 --parse-tko-names --label hela
    igvfagent perturb-tools lfc-reps --screen screen.h5ad --cond1 18 --cond2 8 --rep-col replicate --compare-col time --exclude replicate=0
    igvfagent perturb-tools sort-lfc --screen screen.h5ad --condit-1 high --condit-2 low --control presort --targets MYC_enh1 MYC_enh2
    igvfagent perturb-tools qc --counts readcount-HeLa-lib1 --parse-tko-names --cond1 8 --cond2 18 --compare-col time --pos-ctrl core-essential-genes-sym_HGNCID --target-col GENE
    igvfagent perturb-tools to-mageck --screen screen.h5ad --target-col GENE
    igvfagent perturb-tools read-poolq --poolq-dir poolq_out/ --sample-metadata conditions.csv
    igvfagent perturb-tools design-library --fasta hg38.fa --gtf gencode.gtf --gene KMT2C --chrom chr7
    igvfagent perturb-tools run --counts readcount-HeLa-lib1 --parse-tko-names --cond1 18 --cond2 8 --compare-col time --target-col GENE
    igvfagent perturb-tools selftest --no-plots
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import math
import os
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "PerturbTools"

UPSTREAM_REPO = "IGVF-CRISPR/perturb-tools"
UPSTREAM_COMMIT = "81eaae070be97396a62f87654cf5fc5132474af2"

BLUE, ORANGE, INK, INK2, AXIS, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948",
          "#792374", "#006479", "#0096a0", "#d3a9ce", "#96a0b3", "#5e5948", "#916953", "#888364"]

FIVE_PRIME = "5'"
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "N": "[ACGT]", "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]",
         "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]"}
IUPAC_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M",
              "M": "K", "B": "V", "V": "B", "D": "H", "H": "D"}

log = logging.getLogger("perturb_tools")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"perturb_tools_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in str(label))[:80]


def run_dir(label: str) -> Path:
    base = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d, i = base, 1
    while d.exists():
        i += 1
        d = Path(f"{base}_{i}")
    d.mkdir(parents=True, exist_ok=True)
    return d


def _pd():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit("perturb-tools needs pandas + numpy + scipy: pip install 'igvfagent[analysis]'") from e


def _np():
    import numpy as np  # type: ignore
    return np


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


def _save(fig, path: Path, figs: Optional[list] = None) -> Path:
    fig.patch.set_facecolor(SURFACE)
    fig.savefig(str(path), bbox_inches="tight", dpi=130)
    import matplotlib.pyplot as plt  # type: ignore
    plt.close(fig)
    print(f"Figure: {path}")
    if figs is not None:
        figs.append(str(path))
    return path


def md_table(headers: Sequence[str], rows: Iterable[Iterable[Any]], max_rows: int = 40) -> str:
    lines = ["| " + " | ".join(str(h) for h in headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for i, r in enumerate(rows):
        if i >= max_rows:
            lines.append("| … |" + " |" * (len(headers) - 1))
            break
        lines.append("| " + " | ".join(_fmt(c) if isinstance(c, float) else str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(lines)


def _fmt(x, nd: int = 3) -> str:
    try:
        if x is None or (isinstance(x, float) and not math.isfinite(x)):
            return "nan"
        if isinstance(x, float) and x != 0 and abs(x) < 10 ** -nd:
            return f"{x:.2e}"
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return str(x)


def write_tsv(df, path: Path, index: bool = False) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=index, compression="gzip" if str(path).endswith(".gz") else None)
    print(f"TSV: {path}")
    return path


def _jsonable(o):
    np = _np()
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, float) and not math.isfinite(o):
        return None
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    return o


def finish(out: Path, title: str, sections: List[str], summary: dict) -> None:
    summary = dict(summary)
    summary.setdefault("upstream", {"repo": UPSTREAM_REPO, "ref": UPSTREAM_COMMIT})
    (out / "summary.json").write_text(json.dumps(_jsonable(summary), indent=2, default=str))
    print(f"JSON: {out / 'summary.json'}")
    body = [f"# {title}", "", f"perturb-tools port ({UPSTREAM_REPO} @ {UPSTREAM_COMMIT[:7]}); run dir `{out}`.", ""]
    body.extend(sections)
    (out / "report.md").write_text("\n".join(body) + "\n")
    print(f"Report: {out / 'report.md'}")


def _trapz(y, x) -> float:
    np = _np()
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    return float(np.sum((x[1:] - x[:-1]) * (y[1:] + y[:-1]) / 2.0))


# ---------------------------------------------------------------------------
# Screen object (AnnData layout: samples x guides)
# ---------------------------------------------------------------------------

class Screen:
    """Minimal AnnData-shaped container: X (n_samples x n_guides), obs, var, layers, obsm, obsp, uns."""

    def __init__(self, X, obs, var, layers=None, obsm=None, obsp=None, uns=None):
        np = _np()
        self.X = np.asarray(X, dtype=float)
        self.obs = obs
        self.var = var
        self.layers: Dict[str, Any] = dict(layers or {})
        self.obsm: Dict[str, Any] = dict(obsm or {})
        self.obsp: Dict[str, Any] = dict(obsp or {})
        self.uns: Dict[str, Any] = dict(uns or {})
        if self.X.shape != (len(obs), len(var)):
            raise ValueError(f"X shape {self.X.shape} != (n_obs {len(obs)}, n_vars {len(var)})")

    @property
    def shape(self) -> Tuple[int, int]:
        return self.X.shape

    def copy(self) -> "Screen":
        return Screen(self.X.copy(), self.obs.copy(), self.var.copy(), {k: v.copy() for k, v in self.layers.items()},
                      {k: (v.copy() if hasattr(v, "copy") else v) for k, v in self.obsm.items()},
                      {k: (v.copy() if hasattr(v, "copy") else v) for k, v in self.obsp.items()}, dict(self.uns))

    def subset(self, samples=None, guides=None) -> "Screen":
        np = _np()
        si = np.arange(self.shape[0]) if samples is None else np.asarray(samples)
        gi = np.arange(self.shape[1]) if guides is None else np.asarray(guides)
        if si.dtype == bool:
            si = np.where(si)[0]
        if gi.dtype == bool:
            gi = np.where(gi)[0]
        obsm = {}
        for k, v in self.obsm.items():
            try:
                obsm[k] = v[si][:, si] if getattr(v, "shape", (0, 0))[:2] == (self.shape[0], self.shape[0]) else v[si]
            except Exception:
                pass
        return Screen(self.X[np.ix_(si, gi)], self.obs.iloc[si].copy(), self.var.iloc[gi].copy(),
                      {k: v[np.ix_(si, gi)] for k, v in self.layers.items()}, obsm, {}, dict(self.uns))

    def matrix(self, layer: Optional[str] = None):
        if layer is None or layer == "X":
            return self.X
        if layer not in self.layers:
            raise ValueError(f"Invalid layer key {layer}; layers: {list(self.layers)}")
        return self.layers[layer]

    def __repr__(self) -> str:
        return (f"Screen n_samples x n_guides = {self.shape[0]} x {self.shape[1]}; obs: {list(self.obs.columns)}; "
                f"var: {list(self.var.columns)}; layers: {list(self.layers)}")


def parse_tko_sample_names(obs):
    """Notebook convention: replicate = last char, time = int(name[1:-1]) or -1."""
    obs = obs.copy()
    names = obs.index.astype(str)
    obs["replicate"] = [n[-1] if n else "" for n in names]
    obs["time"] = [int(n[1:-1]) if n[1:-1].strip().lstrip("-").isdigit() else -1 for n in names]
    return obs


def _read_any_table(path: Path, sep: Optional[str] = None):
    pd = _pd()
    p = str(path)
    if sep is None:
        low = p.lower().replace(".gz", "")
        sep = "," if low.endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, low_memory=False)


def screen_from_counts_table(df, guide_id_col: Optional[str] = None, annotation_cols: Optional[List[str]] = None) -> Screen:
    """guides x samples table: guide id column, optional annotation columns (auto: every non-numeric column), sample columns."""
    pd, np = _pd(), _np()
    gid = guide_id_col or df.columns[0]
    if gid not in df.columns:
        raise ValueError(f"guide id column {gid} not in {list(df.columns)}")
    ann = list(annotation_cols or [])
    for c in df.columns:
        if c == gid or c in ann:
            continue
        if not pd.api.types.is_numeric_dtype(df[c]):
            ann.append(c)
    samples = [c for c in df.columns if c != gid and c not in ann]
    var = df[[gid] + ann].copy().set_index(gid)
    var.index = var.index.astype(str)
    obs = pd.DataFrame(index=pd.Index([str(s) for s in samples], name="sample"))
    X = df[samples].to_numpy(dtype=float).T
    return Screen(X, obs, var)


def read_csv_screen(X_path=None, guide_path=None, sample_path=None, sep: str = ",") -> Screen:
    """pt.io.read_csv: X (guides x samples, first column = guide id) transposed; guide and sample tables."""
    pd = _pd()
    if X_path is None:
        raise ValueError("X_path is required")
    X_df = pd.read_csv(X_path, sep=sep, header=0, index_col=0)
    X_df.index = X_df.index.astype(str)
    var = pd.DataFrame(index=X_df.index.copy())
    if guide_path is not None:
        g = pd.read_csv(guide_path, sep=sep)
        key = g.columns[0]
        if set(g[key].astype(str)) >= set(X_df.index):
            g[key] = g[key].astype(str)
            g = g.drop_duplicates(key).set_index(key).reindex(X_df.index)
        elif len(g) == len(X_df):
            g.index = X_df.index
        else:
            raise ValueError("guide table does not match the count table")
        var = g
    obs = pd.DataFrame(index=pd.Index([str(c) for c in X_df.columns], name="sample"))
    if sample_path is not None:
        obs = attach_sample_table(obs, pd.read_csv(sample_path, sep=sep))
    return Screen(X_df.to_numpy(dtype=float).T, obs, var)


def attach_sample_table(obs, s):
    key = "sample" if "sample" in s.columns else s.columns[0]
    s = s.copy()
    s[key] = s[key].astype(str)
    s = s.drop_duplicates(key).set_index(key)
    missing = [i for i in obs.index if i not in s.index]
    if missing:
        raise ValueError(f"sample table lacks {missing[:5]}")
    out = obs.join(s.loc[list(obs.index)], how="left", rsuffix="_meta")
    out.index.name = obs.index.name or "sample"
    return out


def write_csv_screen(screen: Screen, out_path: Path) -> Path:
    """pt.io.write_csv (fixed): screen.X.csv as a guides x samples matrix, plus var / obs / layers / obsm / obsp."""
    pd = _pd()
    out_path.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(screen.X.T, index=screen.var.index, columns=screen.obs.index).to_csv(out_path / "screen.X.csv")
    screen.var.to_csv(out_path / "screen.var.csv")
    screen.obs.to_csv(out_path / "screen.obs.csv")
    for k, v in screen.layers.items():
        pd.DataFrame(v.T, index=screen.var.index, columns=screen.obs.index).to_csv(out_path / f"screen.layers.{safe_label(k)}.csv")
    for kind, dct in (("obsm", screen.obsm), ("obsp", screen.obsp)):
        for k, v in dct.items():
            df = v if hasattr(v, "to_csv") else pd.DataFrame(v, index=screen.obs.index if len(v) == screen.shape[0] else None)
            df.to_csv(out_path / f"screen.{kind}.{safe_label(k)}.csv")
    return out_path


def read_screen_dir(d: Path) -> Screen:
    pd = _pd()
    X_df = pd.read_csv(d / "screen.X.csv", index_col=0)
    X_df.index = X_df.index.astype(str)
    var = pd.read_csv(d / "screen.var.csv", index_col=0) if (d / "screen.var.csv").is_file() else pd.DataFrame(index=X_df.index)
    var.index = var.index.astype(str)
    obs = pd.read_csv(d / "screen.obs.csv", index_col=0) if (d / "screen.obs.csv").is_file() else pd.DataFrame(index=X_df.columns)
    obs.index = obs.index.astype(str)
    layers = {}
    for f in sorted(d.glob("screen.layers.*.csv")):
        key = f.name[len("screen.layers."):-4]
        layers[key] = pd.read_csv(f, index_col=0).to_numpy(dtype=float).T
    obsm = {}
    for f in sorted(d.glob("screen.obsm.*.csv")):
        obsm[f.name[len("screen.obsm."):-4]] = pd.read_csv(f, index_col=0).to_numpy()
    return Screen(X_df.to_numpy(dtype=float).T, obs, var.reindex(X_df.index), layers, obsm)


def _clean_frame_for_h5ad(df):
    df = df.copy()
    df.index = df.index.astype(str)
    df.columns = [str(c) for c in df.columns]
    if df.index.name in df.columns:
        df.index.name = None
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].map(lambda v: "" if v is None or (isinstance(v, float) and math.isnan(v)) else str(v)).astype("category")
    return df


def write_h5ad(screen: Screen, path: Path) -> Optional[Path]:
    try:
        import anndata as ad  # type: ignore
    except Exception:
        print("anndata not installed: .h5ad not written (CSV screen directory written instead)")
        return None
    a = ad.AnnData(X=screen.X.astype("float32"), obs=_clean_frame_for_h5ad(screen.obs), var=_clean_frame_for_h5ad(screen.var))
    for k, v in screen.layers.items():
        a.layers[k] = v
    for k, v in screen.obsm.items():
        try:
            a.obsm[k] = v.to_numpy() if hasattr(v, "to_numpy") else v
        except Exception:
            pass
    a.write_h5ad(str(path))
    print(f"Wrote: {path}")
    return path


def read_h5ad(path: Path) -> Screen:
    try:
        import anndata as ad  # type: ignore
    except Exception as e:
        raise SystemExit("reading .h5ad needs anndata (pip install anndata), or pass --counts / a screen CSV directory") from e
    np = _np()
    a = ad.read_h5ad(str(path))
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    layers = {k: (v.toarray() if hasattr(v, "toarray") else np.asarray(v)) for k, v in a.layers.items()}
    obsm = {k: np.asarray(v) for k, v in a.obsm.items()}
    return Screen(X, a.obs.copy(), a.var.copy(), layers, obsm, {}, {})


def save_screen(screen: Screen, out: Path, name: str = "screen") -> dict:
    d = write_csv_screen(screen, out / name)
    print(f"Wrote: {d}")
    h = write_h5ad(screen, out / f"{name}.h5ad")
    return {"screen_dir": str(d), "h5ad": str(h) if h else None}


def add(screen1: Screen, screen2: Screen) -> Screen:
    """pt.add: sum two screens whose guides and samples match exactly."""
    if (len(screen1.var) == len(screen2.var) and len(screen1.obs) == len(screen2.obs)
            and all(screen1.var.index == screen2.var.index) and all(screen1.obs.index == screen2.obs.index)):
        return Screen(screen1.X + screen2.X, screen1.obs.copy(), screen1.var.copy())
    raise ValueError("Guides/sample description mismatch")


# ---------------------------------------------------------------------------
# Normalisation and log fold change (pt.pp)
# ---------------------------------------------------------------------------

def read_count_normalize(X):
    """Reads per million per sample (rows)."""
    np = _np()
    X = np.asarray(X, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return X / X.sum(axis=1)[:, None] * 1e6


def log_transform_read_count(X):
    return _np().log2(X + 1)


def log_normalize_read_count(X):
    """log2(RPM + 1) (PLoS One 2017 e0170445 methods; identical to PoolQ lognormalized-counts)."""
    return log_transform_read_count(read_count_normalize(X))


def log_norm(screen: Screen, output_layer: str = "lognorm_counts", read_count_layer: Optional[str] = None) -> str:
    if read_count_layer is None:
        screen.layers[output_layer] = log_normalize_read_count(screen.X)
        return output_layer
    output_layer = f"lognorm_{read_count_layer}"
    screen.layers[output_layer] = log_normalize_read_count(screen.layers[read_count_layer])
    return output_layer


def _sample_index(screen: Screen, name) -> int:
    np = _np()
    idx = np.where(screen.obs.index.astype(str) == str(name))[0]
    if len(idx) == 0:
        raise ValueError(f"No sample named {name} in Screen object.")
    if len(idx) > 1:
        raise ValueError(f"Duplicate sample name {name} in Screen object")
    return int(idx[0])


def log_fold_change(screen: Screen, sample1, sample2, lognorm_counts_key: str = "lognorm_counts",
                    out_guides_suffix: str = "lfc", return_result: bool = False, name: Optional[str] = None):
    if "lognorm" not in lognorm_counts_key:
        log.warning("The layer specified must be log-normalized values using log_norm().")
    if lognorm_counts_key not in screen.layers:
        raise ValueError("Specified normalized count isn't in your layer. First run log_norm().")
    i, j = _sample_index(screen, sample1), _sample_index(screen, sample2)
    M = screen.layers[lognorm_counts_key]
    lfc = M[i, :] - M[j, :]
    if return_result:
        return lfc
    col = name or f"{sample1}_{sample2}.{out_guides_suffix}"
    screen.var[col] = lfc
    return col


def fold_change(screen: Screen, cond1, cond2, lognorm_counts_key: str = "lognorm_counts", return_result: bool = False):
    """pt.pp.fold_change: same difference, stored as "{c1}_{c2}.fc"."""
    return log_fold_change(screen, cond1, cond2, lognorm_counts_key, out_guides_suffix="fc", return_result=return_result)


def _match(series, value):
    """obs value equality tolerant of int/str (CLI passes strings)."""
    s = series.astype(str)
    v = str(value)
    try:
        if float(v) == int(float(v)):
            v2 = str(int(float(v)))
            return (s == v) | (s == v2)
    except ValueError:
        pass
    return s == v


def log_fold_change_reps(screen: Screen, cond1, cond2, lognorm_counts_key: str = "lognorm_counts",
                         rep_col: str = "replicate", compare_col: str = "sort", out_guides_suffix: str = "lfc",
                         keep_result: bool = False):
    pd, np = _pd(), _np()
    if rep_col not in screen.obs.columns:
        raise ValueError(f"{rep_col} not in sample features")
    if compare_col not in screen.obs.columns:
        raise ValueError(f"{compare_col} not in sample features")
    reps = list(pd.unique(screen.obs[rep_col]))
    lfcs = []
    for rep in reps:
        r = (screen.obs[rep_col] == rep).to_numpy()
        i1 = np.where(r & _match(screen.obs[compare_col], cond1).to_numpy())[0]
        i2 = np.where(r & _match(screen.obs[compare_col], cond2).to_numpy())[0]
        if len(i1) != 1 or len(i2) != 1:
            raise ValueError(f"Conditions are not unique for each replicate to be aggregated (replicate {rep}: "
                             f"{len(i1)} x {cond1}, {len(i2)} x {cond2}).")
        lfcs.append(log_fold_change(screen, screen.obs.index[i1[0]], screen.obs.index[i2[0]],
                                    lognorm_counts_key=lognorm_counts_key, return_result=True))
    cols = [f"{r}.{cond1}_{cond2}.{out_guides_suffix}" for r in reps]
    df = pd.DataFrame(np.vstack(lfcs).T, index=screen.var.index, columns=cols)
    if keep_result:
        for c in cols:
            screen.var[c] = df[c].to_numpy()
    return df


def aggregate_lfcs(lfcs_df, aggregate_fn: str = "median"):
    np = _np()
    M = lfcs_df.to_numpy(dtype=float)
    if aggregate_fn == "mean":
        v = np.mean(M, axis=1)
    elif aggregate_fn == "median":
        v = np.median(M, axis=1)
    elif aggregate_fn == "sd":
        v = np.std(M, axis=1)
    else:
        raise ValueError("Only 'mean', 'median', and 'sd' are supported for aggregating LFCs.")
    return v


def log_fold_change_aggregate(screen: Screen, cond1, cond2, lognorm_counts_key: str = "lognorm_counts",
                              aggregate_col: str = "replicate", compare_col: str = "sort", out_guides_suffix: str = "lfc",
                              aggregate_fn: str = "median", name: Optional[str] = None, return_result: bool = False,
                              keep_per_replicate: bool = False):
    df = log_fold_change_reps(screen, cond1, cond2, lognorm_counts_key, aggregate_col, compare_col, out_guides_suffix,
                              keep_per_replicate)
    v = aggregate_lfcs(df, aggregate_fn)
    if return_result:
        return v
    col = name or f"{cond1}_{cond2}.{out_guides_suffix}.{aggregate_fn}"
    screen.var[col] = v
    return col


# ---------------------------------------------------------------------------
# Sorting-screen delta LFC (pt.calc.log_fold_change / LogFoldChangeModule)
# ---------------------------------------------------------------------------

def sort_screen_lfc(counts_df, guides, condit_1: str, condit_2: str, control: str,
                    targets: Optional[List[str]] = None, target_col: str = "target"):
    """counts_df: log-normalised guides x samples DataFrame (index aligned with guides)."""
    pd, np = _pd(), _np()
    from scipy import stats  # type: ignore
    keep = np.ones(len(counts_df), dtype=bool)
    if targets:
        if target_col not in guides.columns:
            raise ValueError(f"{target_col} not in guide features")
        keep = guides[target_col].astype(str).isin([str(t) for t in targets]).to_numpy()
    cdf = counts_df.iloc[np.where(keep)[0]]
    g = guides.iloc[np.where(keep)[0]].copy()
    c1, c2, ct = cdf.filter(regex=condit_1), cdf.filter(regex=condit_2), cdf.filter(regex=control)
    for nm, part, rx in (("condit_1", c1, condit_1), ("condit_2", c2, condit_2), ("control", ct, control)):
        if part.shape[1] == 0:
            raise ValueError(f"no sample column matches {nm} regex {rx!r}; samples: {list(counts_df.columns)}")
    overlap = (set(c1.columns) & set(c2.columns)) | (set(c1.columns) & set(ct.columns)) | (set(c2.columns) & set(ct.columns))
    if overlap:
        raise ValueError(f"regexes match the same sample(s) {sorted(overlap)}; make them specific")
    ctrl = ct.to_numpy(dtype=float)
    if ct.shape[1] not in (1, c1.shape[1]) or c1.shape[1] != c2.shape[1]:
        log.warning("control/condition replicate counts differ; baseline = control mean (upstream requires equal shapes)")
        ctrl = ct.to_numpy(dtype=float).mean(axis=1, keepdims=True)
    a = c1.to_numpy(dtype=float) - ctrl
    b = c2.to_numpy(dtype=float) - ctrl
    lfc_mean = a.mean(axis=1) - b.mean(axis=1)
    lfc_sd = (a - b).std(axis=1) if a.shape == b.shape else np.full(len(a), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.array([stats.ttest_ind(a[i], b[i])[1] for i in range(len(a))], dtype=float)
    g["lfc.mean"] = lfc_mean
    g["lfc.stdev"] = lfc_sd
    g["lfc.p_value"] = p
    with np.errstate(divide="ignore"):
        g["-log10(pval)"] = -np.log10(p)
    return g, {"condit_1": list(c1.columns), "condit_2": list(c2.columns), "control": list(ct.columns)}


# ---------------------------------------------------------------------------
# QC (pt.qc, pt.pl)
# ---------------------------------------------------------------------------

def set_sample_correlation(screen: Screen, method: str = "pearson", count_layer: Optional[str] = None,
                           guide_idx=None, prefix: str = "") -> str:
    pd, np = _pd(), _np()
    if method not in ("pearson", "spearman"):
        raise ValueError(f"method must be pearson or spearman, not {method}")
    sub = screen if guide_idx is None else screen.subset(guides=np.asarray(guide_idx))
    M = sub.matrix(count_layer)
    label = count_layer or "X"
    pre = f"{prefix}_" if prefix else ""
    M = np.where(np.isfinite(M), M, np.nan)
    C = pd.DataFrame(M.T, columns=screen.obs.index).corr(method=method).to_numpy()
    key = f"{pre}corr_{label}"
    screen.obsm[key] = C
    screen.obs[f"{pre}mean_corr_{label}"] = C.mean(axis=0)
    screen.obs[f"{pre}median_corr_{label}"] = np.median(C, axis=1)
    return key


def gini_curve(v) -> Tuple[Any, List[float], float]:
    """_G: Lorenz curve on 50 quantile bins and the Gini coefficient (upstream definition)."""
    np = _np()
    v = np.asarray(v, dtype=float)
    bins = np.linspace(0.0, 1.0, 50)
    total = float(np.sum(v))
    yvals = [float(np.sum(v[v <= np.quantile(v, b)]) / total) if total > 0 else float("nan") for b in bins]
    pe_area = _trapz(bins, bins)
    lorenz = _trapz(yvals, bins)
    return bins, yvals, (pe_area - lorenz) / pe_area


def get_outlier_guides(screen: Screen, cond_col: str, count_layer: Optional[str] = None, mad_z_thres: float = 5,
                       abs_RPM_thres: float = 10000, upstream_compat: bool = False):
    pd, np = _pd(), _np()
    label = count_layer or "X"
    M = screen.matrix(count_layer)
    key = f"{label}_RPM"
    if key not in screen.layers:
        screen.layers[key] = read_count_normalize(M)
    R = screen.layers[key]
    if cond_col not in screen.obs.columns:
        raise ValueError(f"{cond_col} not in sample features")
    idx_name = screen.var.index.name or "index"
    rows = []
    for cond in pd.unique(screen.obs[cond_col]):
        si = np.where((screen.obs[cond_col] == cond).to_numpy())[0]
        Rs = R[si, :]
        if upstream_compat:
            med = np.nanmedian(Rs, axis=1)
            for i in range(len(si)):
                if i >= Rs.shape[1]:
                    break
                col = Rs[:, i]
                oi = np.where((col > med * mad_z_thres) & (col > abs_RPM_thres))[0]
                for o in oi:
                    rows.append({idx_name: screen.var.index[o], "sample": screen.obs.index[si[i]], "RPM": Rs[o, i],
                                 cond_col: cond})
            continue
        med = np.nanmedian(Rs, axis=0)
        for i in range(len(si)):
            oi = np.where((Rs[i] > med * mad_z_thres) & (Rs[i] > abs_RPM_thres))[0]
            for o in oi:
                rows.append({idx_name: screen.var.index[o], "sample": screen.obs.index[si[i]], "RPM": float(Rs[i, o]),
                             cond_col: cond, "condition_median_RPM": float(med[o])})
    cols = [idx_name, "sample", "RPM", cond_col] + ([] if upstream_compat else ["condition_median_RPM"])
    return pd.DataFrame(rows, columns=cols)


def sample_qc_table(screen: Screen, count_layer: Optional[str] = None):
    pd, np = _pd(), _np()
    M = screen.matrix(count_layer)
    rows = []
    for i, s in enumerate(screen.obs.index):
        v = M[i, :]
        _, _, g = gini_curve(v)
        rows.append({"sample": s, "total_reads": float(np.nansum(v)), "median_count": float(np.nanmedian(v)),
                     "zero_count_guides": int((v == 0).sum()), "gini": g})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Export (pt.io)
# ---------------------------------------------------------------------------

def to_mageck_input(screen: Screen, out_path: Optional[Path] = None, count_layer: Optional[str] = None,
                    sgrna_column: Optional[str] = None, target_column: str = "target_id", sample_prefix: str = ""):
    pd = _pd()
    M = screen.matrix(count_layer)
    if target_column not in screen.var.columns:
        raise ValueError(f"{target_column} not in guide features {list(screen.var.columns)}")
    df = pd.DataFrame(M.T, columns=[sample_prefix + str(s) for s in screen.obs.index], index=screen.var.index)
    df = df.fillna(0).astype(int)
    if sgrna_column is None or screen.var.index.name == sgrna_column:
        sg = [str(s) for s in screen.var.index]
    elif sgrna_column in screen.var.columns:
        sg = screen.var[sgrna_column].astype(str).tolist()
    else:
        raise ValueError(f"{sgrna_column} not found in guide features.")
    df.insert(0, "sgRNA", [s.replace(" ", "_") for s in sg])
    df.insert(1, "gene", screen.var[target_column].to_numpy())
    keep = df["gene"].map(lambda o: not pd.isnull(o) and bool(o) and str(o) != "")
    df = df.loc[keep.to_numpy()]
    if out_path is not None:
        df.to_csv(out_path, sep="\t", index=False)
    return df


def check_fix_file_extension(filepath: str, extension: str, silent: bool = False) -> str:
    if not str(filepath).endswith(extension):
        if not silent:
            print(f"filepath: {filepath} does not have the correct extension: {extension}; appending it")
        return str(filepath) + extension
    return str(filepath)


def collect_screen_dfs(screen: Screen, include_uns: bool = False, guide_rows: bool = True):
    pd = _pd()

    def mat(m):
        return (pd.DataFrame(m.T, index=screen.var.index, columns=screen.obs.index) if guide_rows
                else pd.DataFrame(m, index=screen.obs.index, columns=screen.var.index))

    dfs, names = [mat(screen.X)], ["X"]
    for k, m in screen.layers.items():
        dfs.append(mat(m)); names.append(k)
    dfs.extend([screen.var, screen.obs]); names.extend(["guides", "samples"])
    for kind, dct in (("obsm", screen.obsm), ("obsp", screen.obsp)):
        for k, v in dct.items():
            dfs.append(v if hasattr(v, "to_excel") else pd.DataFrame(v)); names.append(f"{kind}.{k}")
    if include_uns:
        for k, v in screen.uns.items():
            try:
                dfs.append(v if hasattr(v, "to_excel") else pd.DataFrame(v)); names.append(f"uns.{k}")
            except Exception:
                pass
    return dfs, names


def write_screen_to_excel(screen: Screen, workbook_path: str, index: bool = True, include_uns: bool = False,
                          guide_rows: bool = True, extra: Optional[List[Tuple[str, Any]]] = None) -> Optional[str]:
    pd = _pd()
    try:
        import openpyxl  # type: ignore  # noqa: F401
    except Exception:
        print("openpyxl not installed: Excel workbook not written (pip install openpyxl)")
        return None
    dfs, names = collect_screen_dfs(screen, include_uns, guide_rows)
    for nm, df in extra or []:
        dfs.append(df); names.append(nm)
    workbook_path = check_fix_file_extension(workbook_path, ".xlsx", silent=True)
    used = set()
    with pd.ExcelWriter(workbook_path) as w:
        for n, (df, nm) in enumerate(zip(dfs, names)):
            sheet = re.sub(r"[\[\]\:\*\?\/\\]", "_", nm)[:31] or f"Sheet_{n}"
            base, k = sheet, 1
            while sheet in used:
                k += 1
                sheet = f"{base[:28]}_{k}"
            used.add(sheet)
            try:
                df.to_excel(w, sheet_name=sheet, index=index)
            except ValueError as e:
                print(f"While writing sheet {nm}: {e}; skipped")
    print(f"Wrote: {workbook_path}")
    return workbook_path


# ---------------------------------------------------------------------------
# PoolQ reader
# ---------------------------------------------------------------------------

def _read_lines(path: Path) -> List[str]:
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return fh.readlines()


def parse_poolq_quality(path: Path) -> Dict[str, Any]:
    pd = _pd()
    meta, bc, common = [], [], []
    bc_cols, common_cols = None, None
    for line in _read_lines(path):
        if line == "\n" or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if f[0].startswith("Read counts"):
            continue
        if ":" in f[0]:
            parts = f[0].split(":", 1)
            meta.append([parts[0].strip(), parts[1].strip()])
        elif len(f) > 2:
            if bc_cols is None:
                bc_cols = f
            else:
                bc.append(f)
        elif len(f) == 2:
            if common_cols is None:
                common_cols = f
            else:
                common.append(f)
    return {"metadata": pd.DataFrame(meta, columns=["metric", "statistic"]),
            "bc_readcounts": pd.DataFrame(bc, columns=bc_cols) if bc_cols else pd.DataFrame(),
            "common_barcodes": pd.DataFrame(common, columns=["barcode", "count"]) if common_cols or common else pd.DataFrame(columns=["barcode", "count"])}


def read_poolq(poolq_dir: Path, sample_metadata: Optional[Path] = None, merge_on: str = "Condition") -> Screen:
    pd, np = _pd(), _np()
    d = Path(poolq_dir)
    counts_p = d / "counts.txt" if (d / "counts.txt").is_file() else d / "expected-counts.txt"
    if not counts_p.is_file():
        raise SystemExit(f"{d}: no counts.txt / expected-counts.txt")
    try:
        cdf = pd.read_csv(counts_p, sep="\t", encoding_errors="replace")
    except TypeError:  # pandas < 1.3
        cdf = pd.read_csv(counts_p, sep="\t")
    idc = list(cdf.columns[:2])
    samples = [str(c) for c in cdf.columns[2:]]
    var = pd.DataFrame({"barcode": cdf[idc[0]].astype(str).to_numpy(), "barcode_id": cdf[idc[1]].astype(str).to_numpy()})
    var.index = pd.Index(var["barcode_id"].astype(str), name="barcode_id")
    if var.index.duplicated().any():
        var.index = pd.Index(var["barcode"].astype(str), name="barcode")
    obs = pd.DataFrame({"Condition": samples}, index=pd.Index(samples, name="sample"))
    X = cdf[cdf.columns[2:]].to_numpy(dtype=float).T
    layers, uns, obsp = {}, {}, {}
    ln = d / "lognormalized-counts.txt"
    if ln.is_file():
        ldf = pd.read_csv(ln, sep="\t")
        layers["lognorm_counts"] = ldf[ldf.columns[2:]].to_numpy(dtype=float).T
    for nm in ("barcode-counts", "unexpected-sequences"):
        p = d / f"{nm}.txt"
        if p.is_file():
            uns[nm] = pd.read_csv(p, sep="\t")
    cp = d / "correlation.txt"
    if cp.is_file():
        rows = [l.rstrip("\n").split("\t") for l in _read_lines(cp) if l.strip()]
        hdr = rows[0][1:]
        vals = []
        for r in rows[1:]:
            vals.append([_to_float(x) for x in r[1:]])
        C = pd.DataFrame(vals, index=[r[0] for r in rows[1:]], columns=hdr)
        obsp["poolq_correlation"] = C.reindex(index=samples, columns=samples).to_numpy()
    qp = d / "quality.txt"
    if qp.is_file():
        for k, v in parse_poolq_quality(qp).items():
            uns[f"quality_{k}"] = v
    rp = d / "runinfo.txt"
    if rp.is_file():
        uns["run_info"] = "".join(_read_lines(rp))
    lp = d / "poolq3.log"
    if lp.is_file():
        uns["poolq3"] = "".join(_read_lines(lp))
    if sample_metadata:
        meta = _read_any_table(Path(sample_metadata))
        if merge_on not in meta.columns:
            meta = meta.rename(columns={meta.columns[0]: merge_on})
        meta[merge_on] = meta[merge_on].astype(str)
        merged = obs.reset_index().merge(meta.drop_duplicates(merge_on), on=merge_on, how="left").set_index("sample")
        obs = merged
    s = Screen(X, obs, var, layers, {}, obsp, uns)
    return s


def _to_float(x: str) -> float:
    try:
        return float(x)
    except ValueError:
        return float("nan")


# ---------------------------------------------------------------------------
# Sequences, FASTA, GTF
# ---------------------------------------------------------------------------

def complement(seq: str) -> str:
    return "".join(COMPLEMENT.get(c, "N") for c in seq.upper())


def reverse_complement(seq: str) -> str:
    return complement(seq)[::-1]


def read_fasta(path: Path, chroms: Optional[Iterable[str]] = None) -> Dict[str, str]:
    """Chromosome sequences (upper case) by the first word of the header; pysam used for an indexed FASTA."""
    want = set(chroms) if chroms is not None else None
    p = str(path)
    if want is not None and (Path(p + ".fai").is_file()):
        try:
            import pysam  # type: ignore
            fa = pysam.FastaFile(p)
            return {c: fa.fetch(c).upper() for c in want if c in fa.references}
        except Exception:
            pass
    opener = gzip.open if p.endswith(".gz") else open
    out: Dict[str, List[str]] = {}
    cur = None
    with opener(p, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0] if line[1:].strip() else ""
                cur = name if (want is None or name in want) else None
                if cur is not None:
                    out[cur] = []
                if want is not None and cur is None and want <= set(out):
                    break
            elif cur is not None:
                out[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in out.items()}


def gtf_gene_features(gtf_path: Path, gene: str, chromosome: Optional[str] = None, feature: str = "exon",
                      source: Optional[str] = None):
    """Rows of a GTF for gene_name == gene (exact) and the requested feature -> Chromosome, Start, End (1-based GTF)."""
    pd = _pd()
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    rows = []
    pat = re.compile(r'gene_name "([^"]+)"')
    with opener(str(gtf_path), "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != feature:
                continue
            if chromosome and f[0] != chromosome:
                continue
            if source and f[1] != source:
                continue
            m = pat.search(f[8])
            if m and m.group(1) == gene:
                rows.append({"Chromosome": f[0], "Start": int(f[3]), "End": int(f[4]), "strand": f[6]})
    return pd.DataFrame(rows, columns=["Chromosome", "Start", "End", "strand"])


def identify_unique_fragments(df, start_key: str = "Start", end_key: str = "End"):
    df = df.reset_index(drop=True)
    a = df.sort_values(start_key, kind="mergesort").reset_index(drop=True).drop_duplicates(subset=end_key, keep="first")
    b = a.sort_values(end_key, kind="mergesort").reset_index(drop=True).drop_duplicates(subset=start_key, keep="first")
    return b.reset_index(drop=True)


def total_overlapping_fragments(df, start_key: str = "Start", end_key: str = "End") -> List[int]:
    """Indices of fragments that totally engulf another (strict on both ends); these are dropped upstream."""
    s = df[start_key].to_numpy(); e = df[end_key].to_numpy()
    out = set()
    for i in range(len(df)):
        for j in range(len(df)):
            if s[i] < s[j] and e[i] > e[j]:
                out.add(i)
    return sorted(out)


def partial_overlapping_fragments(df, start_key: str = "Start", end_key: str = "End") -> Dict[str, List[List[int]]]:
    s = df[start_key].to_numpy(); e = df[end_key].to_numpy()
    res: Dict[str, List[List[int]]] = {"5'": [], "3'": []}
    for i in range(len(df)):
        for j in range(len(df)):
            if s[j] < e[i] and e[j] > e[i]:
                res["5'"].append([i, j])
            if s[j] < s[i] and e[j] > s[i]:
                res["3'"].append([i, j])
    return res


def resolve_fragments(df):
    """_OverlappingFragments.find_all: dedupe, unique ends, drop engulfing fragments, report partial overlaps."""
    df = identify_unique_fragments(df.drop_duplicates())
    total = total_overlapping_fragments(df)
    df = df.drop(total).reset_index(drop=True)
    return df, {"dropped_engulfing": len(total), "partial": partial_overlapping_fragments(df)}


def widen_regions(regions, amount: int):
    regions = regions.copy()
    regions["Start"] = regions["Start"] - amount
    regions["End"] = regions["End"] + amount
    return regions


def pam_regex(pam: str) -> "re.Pattern":
    return re.compile("(?=(" + "".join(IUPAC.get(c, c) for c in pam.upper()) + "))")


def find_pam_loci(seq: str, pam: str = "NGG", offset: int = 0, chrom_len: Optional[int] = None):
    """Forward PAM starts and reverse-strand PAM_loci (upstream L - start in the reverse complement) in forward coords."""
    fwd = [m.start() + offset for m in pam_regex(pam).finditer(seq)]
    rc_pam = "".join(IUPAC_COMP.get(c, c) for c in pam.upper())[::-1]
    rev = [m.start() + offset + len(pam) for m in pam_regex(rc_pam).finditer(seq)]
    return fwd, rev


def annotate_guide_biochemistry(df, col: str = "protospacer"):
    df = df.copy()
    s = df[col].astype(str).str.upper()
    for nm, motif in (("polyT", "TTTT"), ("polyG", "GGGG"), ("polyC", "CCCC"), ("polyA", "AAAA")):
        df[nm] = s.str.contains(motif, regex=False)
    df["GC_content"] = (s.str.count("G") + s.str.count("C")) * 100 / 20
    return df


def scan_for_bsmbi(df, motifs: Sequence[str] = ("TCTC", "GAGAC"), upstream_compat: bool = False):
    flags = []
    for ctx in df["guide_context"].astype(str):
        rc = reverse_complement(ctx)
        hit = False
        for s in (ctx, rc):
            for m in motifs:
                pos = s.find(m)
                if (pos > 0) if upstream_compat else (pos >= 0):
                    hit = True
        flags.append(hit)
    df = df.copy()
    df["BsmbI_flag"] = flags
    return df


def make_guide_library(regions, seqs: Dict[str, str], pam: str = "NGG", flank: int = 10, protospacer_len: int = 20):
    """PAM scan of each region (Start, End] -> one row per (PAM, strand) with protospacer, PAM and context."""
    pd, np = _pd(), _np()
    rows = []
    plen = len(pam)
    for _, r in regions.iterrows():
        chrom = str(r["Chromosome"])
        if chrom not in seqs:
            raise ValueError(f"chromosome {chrom} not in the FASTA")
        seq = seqs[chrom]
        L = len(seq)
        s, e = int(r["Start"]), int(r["End"])
        lo = max(0, s + 1 - plen)
        win = seq[lo:min(L, e + plen)]
        fwd, rev = find_pam_loci(win, pam, offset=lo)
        center = int(np.round((s + e) / 2.0))
        for strand, loci in (("+", fwd), ("-", rev)):
            for p in loci:
                if not (s < p <= e):
                    continue
                if strand == "+":
                    gb, ge = p - protospacer_len, p
                    proto = seq[max(0, gb):ge]
                    pam_seq = seq[p:p + plen]
                    ctx = seq[max(0, p - protospacer_len - flank):p + plen + flank]
                else:
                    gb, ge = p + protospacer_len, p
                    proto = reverse_complement(seq[p:gb])
                    pam_seq = reverse_complement(seq[max(0, p - plen):p])
                    ctx = seq[max(0, p - plen - flank):p + protospacer_len + flank]
                rows.append({"Chromosome": chrom, "PAM_loci": int(p), "region_left": s, "region_right": e, "strand": strand,
                             "region_center": center, "PAM_distance_to_center": int(round(abs(p - center))),
                             "guide_begin": int(gb), "guide_end": int(ge), "PAM": pam_seq, "protospacer": proto,
                             "guide_context": ctx})
    cols = ["Chromosome", "PAM_loci", "region_left", "region_right", "strand", "region_center", "PAM_distance_to_center",
            "guide_begin", "guide_end", "PAM", "protospacer", "guide_context"]
    df = pd.DataFrame(rows, columns=cols)
    df = df.sort_values(["Chromosome", "PAM_loci"], kind="mergesort").drop_duplicates("protospacer").reset_index(drop=True)
    return df


def annotate_protospacer(guide_df, barcode_col: str = "barcode", protospacer_length: int = 20):
    guide_df = guide_df.copy()
    if barcode_col in guide_df.columns:
        lens = guide_df[barcode_col].astype(str).str.len().unique()
        if len(lens) and all(int(x) == protospacer_length for x in lens):
            guide_df["protospacer"] = guide_df[barcode_col].astype(str).str.upper()
    return guide_df


def add_guide_target_metadata(guide_df, annotations: Optional[List[str]] = None, direct_pairs: Optional[Dict[str, str]] = None,
                              id_col: str = "barcode_id"):
    guide_df = guide_df.copy()
    anns = list(annotations or [])
    targets = []
    for gid in guide_df[id_col].astype(str):
        if direct_pairs and gid in direct_pairs:
            targets.append(direct_pairs[gid])
            continue
        hit = next((a for a in anns if a in gid), None)
        targets.append(hit if hit is not None else gid)
    guide_df["target"] = targets
    return guide_df


def annotate_guide_position(guide_df, regions, seqs: Dict[str, str], units: str = "bp"):
    """Position each protospacer inside its target region: + strand first, else the reverse complement."""
    pd, np = _pd(), _np()
    scale = 1e6 if units == "mb" else 1.0
    reg = {}
    for _, r in regions.iterrows():
        name = str(r["target"] if "target" in regions.columns else r[regions.columns[0]])
        s, e = int(round(float(r["Start"]) * scale)), int(round(float(r["End"]) * scale))
        seq = seqs[str(r["Chromosome"])][s:e]
        reg[name] = (str(r["Chromosome"]), s, e, seq, reverse_complement(seq))
    strand, chrom, start, end = [], [], [], []
    n_miss = 0
    for t, sp in zip(guide_df["target"].astype(str), guide_df["protospacer"].astype(str)):
        res = (None, None, float("nan"), float("nan"))
        if t in reg:
            c, s, e, plus, minus = reg[t]
            m = re.search(re.escape(sp), plus)
            if m:
                res = ("+", c, s + m.span()[0], s + m.span()[1])
            else:
                m = re.search(re.escape(sp), minus)
                if m:
                    res = ("-", c, e - m.span()[0], e - m.span()[1])
        if res[0] is None:
            n_miss += 1
        strand.append(res[0]); chrom.append(res[1]); start.append(res[2]); end.append(res[3])
    out = guide_df.copy()
    out["strand"] = strand
    out["Chromosome"] = chrom
    out["Start"] = start
    out["End"] = end
    out["center"] = np.mean([np.asarray(start, dtype=float), np.asarray(end, dtype=float)], axis=0)
    return out, n_miss


def pair_distance(df1, df2, chrom_col: str = "chr", start_col: str = "start", end_col: str = "end", how: str = "start",
                  id1: Optional[str] = None, id2: Optional[str] = None, max_distance: Optional[float] = None):
    """Same-chromosome pairs; distance |start1 - start2| (upstream), midpoint, or shortest gap (0 if overlapping)."""
    pd, np = _pd(), _np()
    for nm, df in (("first", df1), ("second", df2)):
        for c in (chrom_col, start_col, end_col):
            if c not in df.columns:
                raise ValueError(f"{c} not in the {nm} table columns {list(df.columns)}")
    a = df1.copy(); b = df2.copy()
    a["_id1"] = a[id1].astype(str) if id1 else a.index.astype(str)
    b["_id2"] = b[id2].astype(str) if id2 else b.index.astype(str)
    a["_chr"] = a[chrom_col].astype(str); b["_chr"] = b[chrom_col].astype(str)
    out = []
    for c, sa in a.groupby("_chr", sort=False):
        sb = b[b["_chr"] == c]
        if sb.empty:
            continue
        s1 = sa[start_col].to_numpy(float)[:, None]; e1 = sa[end_col].to_numpy(float)[:, None]
        s2 = sb[start_col].to_numpy(float)[None, :]; e2 = sb[end_col].to_numpy(float)[None, :]
        if how == "start":
            D = np.abs(s1 - s2)
        elif how == "midpoint":
            D = np.abs((s1 + e1) / 2 - (s2 + e2) / 2)
        elif how == "shortest":
            D = np.maximum(0, np.maximum(s2 - e1, s1 - e2))
        else:
            raise ValueError("how must be start, midpoint or shortest")
        ii, jj = np.nonzero(np.ones_like(D, dtype=bool) if max_distance is None else D <= max_distance)
        out.append(pd.DataFrame({"id1": sa["_id1"].to_numpy()[ii], "id2": sb["_id2"].to_numpy()[jj], "chr": c,
                                 "distance": D[ii, jj]}))
    if not out:
        return pd.DataFrame(columns=["id1", "id2", "chr", "distance"])
    return pd.concat(out, ignore_index=True)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def fig_count_dist(screen: Screen, path: Path, figs: list, count_layer: Optional[str] = None, n_bins: int = 100) -> None:
    plt, np = _plt(), _np()
    if plt is None:
        return
    M = screen.matrix(count_layer)
    bins = 10 ** np.linspace(0, np.log10(np.nanmax(M) + 1), n_bins)
    fig, ax = plt.subplots(figsize=(7, 4))
    for i, s in enumerate(screen.obs.index):
        ax.hist(M[i], bins, histtype="step", color=SERIES[i % len(SERIES)], label=f"{s} (median={np.nanmedian(M[i]):.3g})")
    ax.set_xscale("log")
    _style(ax, "Per-sample guide count distribution", "read count per guide", "# guides")
    ax.legend(bbox_to_anchor=(1.02, 1), fontsize=7, frameon=False)
    _save(fig, path, figs)


def fig_gini(screen: Screen, path: Path, figs: list, count_layer: Optional[str] = None) -> None:
    plt = _plt()
    if plt is None:
        return
    M = screen.matrix(count_layer)
    fig, ax = plt.subplots(figsize=(7, 4))
    for i, s in enumerate(screen.obs.index):
        b, y, g = gini_curve(M[i])
        ax.plot(b, y, color=SERIES[i % len(SERIES)], label=f"{s} ({g:.3f})")
    ax.plot((0, 1), (0, 1), "--", color=INK2, label="perfect eq.")
    _style(ax, "Gini index (guide coverage uniformity)", "cumulative fraction of guides, lowest to highest count",
           "cumulative fraction of reads")
    ax.legend(bbox_to_anchor=(1.02, 1), fontsize=7, frameon=False)
    _save(fig, path, figs)


def fig_heatmap(C, labels, title: str, path: Path, figs: list, annot: bool = True) -> None:
    plt, np = _plt(), _np()
    if plt is None:
        return
    n = len(labels)
    fig, ax = plt.subplots(figsize=(max(4, 0.5 * n + 2), max(3.5, 0.5 * n + 1.5)))
    im = ax.imshow(np.asarray(C, dtype=float), cmap="viridis", vmin=min(0, float(np.nanmin(C))), vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(n)); ax.set_yticklabels(labels, fontsize=7)
    if annot and n <= 20:
        for i in range(n):
            for j in range(n):
                ax.text(j, i, f"{C[i][j]:.2f}", ha="center", va="center", fontsize=6, color="white" if C[i][j] < 0.7 else INK)
    fig.colorbar(im, ax=ax, shrink=0.7)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold", color=INK)
    _save(fig, path, figs)


def fig_rep_consistency(lfcs, path: Path, figs: list) -> None:
    plt, np = _plt(), _np()
    if plt is None or lfcs.shape[1] < 2:
        return
    from scipy import stats  # type: ignore
    cols = list(lfcs.columns)
    k = len(cols)
    fig, axes = plt.subplots(k, k, figsize=(2.2 * k, 2.2 * k), squeeze=False)
    for i in range(k):
        for j in range(k):
            ax = axes[i][j]
            if j > i:
                ax.axis("off"); continue
            x = lfcs[cols[j]].to_numpy(float); y = lfcs[cols[i]].to_numpy(float)
            if i == j:
                ax.hist(x[np.isfinite(x)], bins=40, color=BLUE)
            else:
                ok = np.isfinite(x) & np.isfinite(y)
                ax.scatter(x[ok], y[ok], s=3, color=BLUE, alpha=0.4)
                if ok.sum() > 2:
                    ax.annotate(f"r={stats.pearsonr(x[ok], y[ok])[0]:.2f}", xy=(0.05, 0.88), xycoords="axes fraction", fontsize=8)
            _style(ax, "", cols[j] if i == k - 1 else "", cols[i] if j == 0 else "")
    fig.suptitle("Replicate consistency of LFC", fontsize=10, fontweight="bold")
    _save(fig, path, figs)


def fig_sort_lfc(res, path_volcano: Path, path_pos: Path, figs: list, exons=None, gene: Optional[str] = None) -> None:
    plt, np = _plt(), _np()
    if plt is None:
        return
    fig, ax = plt.subplots(figsize=(5.5, 4))
    y = res["-log10(pval)"].replace([np.inf], np.nan)
    ax.scatter(res["lfc.mean"], y, s=6, color=BLUE, alpha=0.6)
    ax.axvline(0, color=AXIS, lw=0.8)
    _style(ax, "Sorting screen: delta LFC vs control", "Δlog₂ fold change (condit_1 - condit_2)", "-log10(p)")
    _save(fig, path_volcano, figs)
    if "center" in res.columns and res["center"].notna().any():
        r = res[res["center"].notna()]
        fig, ax = plt.subplots(figsize=(9, 3.6))
        for _, row in r.iterrows():
            lo, hi = sorted([row["Start"], row["End"]])
            ax.fill_between([lo, hi], 0, row["lfc.mean"], color="navy", alpha=0.1, lw=0)
        ax.scatter(r["center"], r["lfc.mean"], s=8, color=BLUE)
        ax.axhline(0, color=ORANGE, lw=0.6, ls="--")
        if exons is not None and len(exons):
            ymin = r["lfc.mean"].min() - (r["lfc.mean"].max() - r["lfc.mean"].min()) / 10
            ax.hlines(ymin - 0.15, exons["Start"].min(), exons["End"].max(), color=ORANGE, lw=2)
            for _, e in exons.iterrows():
                ax.fill_between([e["Start"], e["End"]], ymin - 0.3, ymin, color=ORANGE, alpha=0.35, lw=0)
            if gene:
                ax.text((exons["Start"].min() + exons["End"].max()) / 2, ymin + 0.05, gene, ha="center", fontsize=8)
        chrom = r["Chromosome"].dropna().unique()
        _style(ax, "Guide enrichment along the locus", f"genomic coordinate ({chrom[0] if len(chrom) else ''})",
               "Δlog₂ fold change")
        _save(fig, path_pos, figs)


def fig_lfc_hist(values, groups, title: str, path: Path, figs: list) -> None:
    plt, np = _plt(), _np()
    if plt is None:
        return
    fig, ax = plt.subplots(figsize=(6, 3.6))
    if groups is None:
        ax.hist(values[np.isfinite(values)], bins=60, color=BLUE)
    else:
        for k, (nm, m) in enumerate(groups.items()):
            v = values[m & np.isfinite(values)]
            if len(v):
                ax.hist(v, bins=60, histtype="step", color=SERIES[k % len(SERIES)], label=f"{nm} (n={len(v)})", density=True)
        ax.legend(fontsize=7, frameon=False)
    _style(ax, title, "log₂ fold change", "guides")
    _save(fig, path, figs)


# ---------------------------------------------------------------------------
# Screen loading from CLI arguments
# ---------------------------------------------------------------------------

def add_screen_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("screen input (one of --screen / --counts)")
    g.add_argument("--screen", help=".h5ad or a screen CSV directory written by this skill")
    g.add_argument("--counts", help="guides x samples count table (TSV/CSV; first column guide id; text columns become guide annotations)")
    g.add_argument("--guides", help="guide annotation table (first column = guide id)")
    g.add_argument("--samples", help="sample annotation table (column 'sample' or first column = sample name)")
    g.add_argument("--guide-id-col", help="guide id column of --counts (default: first column)")
    g.add_argument("--annotation-cols", nargs="+", help="extra columns of --counts to treat as guide annotation")
    g.add_argument("--parse-tko-names", action="store_true", help="derive replicate (last char) and time (int of name[1:-1]) from sample names")
    g.add_argument("--exclude", nargs="+", default=[], help="drop samples by obs value, e.g. replicate=0")


def load_screen(args: argparse.Namespace) -> Screen:
    pd = _pd()
    if getattr(args, "screen", None):
        p = Path(args.screen)
        s = read_screen_dir(p) if p.is_dir() else read_h5ad(p)
    elif getattr(args, "counts", None):
        s = screen_from_counts_table(_read_any_table(Path(args.counts)), args.guide_id_col, args.annotation_cols)
    else:
        raise SystemExit("pass --screen or --counts")
    if getattr(args, "guides", None):
        g = _read_any_table(Path(args.guides))
        key = g.columns[0]
        g[key] = g[key].astype(str)
        g = g.drop_duplicates(key).set_index(key)
        extra = [c for c in g.columns if c not in s.var.columns]
        s.var = s.var.join(g[extra], how="left")
    if getattr(args, "samples", None):
        s.obs = attach_sample_table(s.obs, _read_any_table(Path(args.samples)))
    if getattr(args, "parse_tko_names", False):
        s.obs = parse_tko_sample_names(s.obs)
    for ex in getattr(args, "exclude", None) or []:
        if "=" not in ex:
            raise SystemExit(f"--exclude expects col=value, got {ex}")
        col, val = ex.split("=", 1)
        if col not in s.obs.columns:
            raise SystemExit(f"--exclude: {col} not in sample features {list(s.obs.columns)}")
        keep = ~_match(s.obs[col], val).to_numpy()
        s = s.subset(samples=keep)
    if s.obs.index.name is None:
        s.obs.index.name = "sample"
    _ = pd
    return s


def _screen_section(s: Screen) -> str:
    return (f"Screen: {s.shape[0]} samples x {s.shape[1]} guides; sample features {list(s.obs.columns)}; "
            f"guide features {list(s.var.columns)}; layers {list(s.layers)}.")


def _ensure_lognorm(s: Screen, key: str, count_layer: Optional[str] = None) -> str:
    if key in s.layers:
        return key
    return log_norm(s, output_layer=key, read_count_layer=count_layer)


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_make_screen(args: argparse.Namespace) -> int:
    s = load_screen(args)
    out = run_dir(args.label)
    if args.log_norm:
        log_norm(s)
    files = save_screen(s, out)
    write_tsv(s.obs.reset_index(), out / "samples.tsv")
    write_tsv(s.var.reset_index(), out / "guides.tsv")
    finish(out, "Screen object", [_screen_section(s), "", md_table(["sample"] + list(s.obs.columns),
                                                                     ([i] + list(r) for i, r in s.obs.iterrows()))],
           {"n_samples": s.shape[0], "n_guides": s.shape[1], "obs_columns": list(s.obs.columns),
            "var_columns": list(s.var.columns), "layers": list(s.layers), **files})
    return 0


def cmd_add_screens(args: argparse.Namespace) -> int:
    screens = []
    for p in args.screens:
        pp = Path(p)
        screens.append(read_screen_dir(pp) if pp.is_dir() else read_h5ad(pp))
    total = screens[0]
    for s in screens[1:]:
        total = add(total, s)
    out = run_dir(args.label)
    files = save_screen(total, out)
    finish(out, "Added screens", [f"Summed {len(screens)} screens (technical replicates).", _screen_section(total)],
           {"n_screens": len(screens), "n_samples": total.shape[0], "n_guides": total.shape[1],
            "total_reads": float(total.X.sum()), **files})
    return 0


def cmd_log_norm(args: argparse.Namespace) -> int:
    s = load_screen(args)
    key = log_norm(s, read_count_layer=args.count_layer)
    out = run_dir(args.label)
    files = save_screen(s, out)
    pd = _pd()
    write_tsv(pd.DataFrame(s.layers[key].T, index=s.var.index, columns=s.obs.index).reset_index(), out / f"{key}.tsv")
    finish(out, "Log-normalised counts", [f"Layer `{key}` = log2(reads per million + 1).", _screen_section(s)],
           {"layer": key, "n_samples": s.shape[0], "n_guides": s.shape[1], **files})
    return 0


def cmd_lfc(args: argparse.Namespace) -> int:
    s = load_screen(args)
    key = _ensure_lognorm(s, args.lognorm_key, args.count_layer)
    col = log_fold_change(s, args.sample1, args.sample2, key, out_guides_suffix=args.suffix)
    out = run_dir(args.label)
    write_tsv(s.var.reset_index(), out / "guides_lfc.tsv")
    files = save_screen(s, out)
    np = _np()
    v = s.var[col].to_numpy(float)
    figs: list = []
    if not args.no_plots:
        fig_lfc_hist(v, None, f"{col}", out / "lfc_hist.png", figs)
    finish(out, "Log fold change", [f"`{col}` = lognorm[{args.sample1}] - lognorm[{args.sample2}] over {len(v)} guides; "
                                    f"median {np.nanmedian(v):.3f}, sd {np.nanstd(v):.3f}."],
           {"column": col, "median": float(np.nanmedian(v)), "sd": float(np.nanstd(v)), "figures": figs, **files})
    return 0


def _group_masks(s: Screen, target_col: Optional[str], pos_ctrl: Optional[set], neg_ctrl: Optional[set]):
    np = _np()
    if not target_col or target_col not in s.var.columns:
        return None
    t = s.var[target_col].astype(str)
    groups = {}
    if pos_ctrl:
        groups["positive controls"] = t.isin(pos_ctrl).to_numpy()
    if neg_ctrl:
        groups["negative controls"] = t.isin(neg_ctrl).to_numpy()
    if groups:
        other = np.ones(len(t), dtype=bool)
        for m in groups.values():
            other &= ~m
        groups["other"] = other
        return groups
    return None


def _read_gene_list(p: Optional[str]) -> Optional[set]:
    if not p:
        return None
    pd = _pd()
    path = Path(p)
    if path.is_file():
        df = pd.read_csv(path, sep="\t", header=None, comment="#")
        return set(df.iloc[:, 0].astype(str))
    return set(x for x in re.split(r"[,\s]+", p) if x)


def cmd_lfc_reps(args: argparse.Namespace) -> int:
    np = _np()
    s = load_screen(args)
    key = _ensure_lognorm(s, args.lognorm_key, args.count_layer)
    reps = log_fold_change_reps(s, args.cond1, args.cond2, key, args.rep_col, args.compare_col, args.suffix, keep_result=True)
    agg_cols = []
    for fn in args.aggregate:
        agg_cols.append(log_fold_change_aggregate(s, args.cond1, args.cond2, key, args.rep_col, args.compare_col, args.suffix, fn))
    out = run_dir(args.label)
    write_tsv(s.var.reset_index(), out / "guides_lfc.tsv")
    files = save_screen(s, out)
    corr = reps.corr(method=args.method)
    write_tsv(corr.reset_index().rename(columns={"index": "lfc"}), out / "lfc_replicate_correlation.tsv")
    figs: list = []
    pos = _read_gene_list(args.pos_ctrl)
    neg = _read_gene_list(args.neg_ctrl)
    groups = _group_masks(s, args.target_col, pos, neg)
    med_col = next((c for c in agg_cols if c.endswith(".median")), agg_cols[0] if agg_cols else None)
    grp_rows = []
    if groups and med_col:
        v = s.var[med_col].to_numpy(float)
        for nm, m in groups.items():
            grp_rows.append([nm, int(m.sum()), float(np.nanmedian(v[m])) if m.any() else float("nan")])
    if not args.no_plots:
        fig_rep_consistency(reps, out / "replicate_consistency.png", figs)
        fig_heatmap(corr.to_numpy(), list(corr.columns), f"LFC {args.method} correlation", out / "lfc_correlation.png", figs)
        if med_col:
            fig_lfc_hist(s.var[med_col].to_numpy(float), groups, med_col, out / "lfc_aggregate_hist.png", figs)
    sections = [f"Per-replicate LFC `{args.cond1}` vs `{args.cond2}` (replicates from `{args.rep_col}`, conditions from "
                f"`{args.compare_col}`): {list(reps.columns)}; aggregates {agg_cols}.", "",
                "## Replicate LFC correlation", "", md_table(["lfc"] + list(corr.columns), ([i] + list(r) for i, r in corr.iterrows()))]
    if grp_rows:
        sections += ["", f"## `{med_col}` by guide group", "", md_table(["group", "guides", "median LFC"], grp_rows)]
    finish(out, "Replicate log fold changes", sections,
           {"replicate_columns": list(reps.columns), "aggregate_columns": agg_cols,
            "replicate_correlation": corr.round(4).to_dict(), "groups": {r[0]: {"n": r[1], "median_lfc": r[2]} for r in grp_rows},
            "figures": figs, **files})
    return 0


def cmd_sort_lfc(args: argparse.Namespace) -> int:
    pd, np = _pd(), _np()
    s = load_screen(args)
    key = _ensure_lognorm(s, args.lognorm_key, args.count_layer)
    counts_df = pd.DataFrame(s.layers[key].T, index=s.var.index, columns=[str(c) for c in s.obs.index])
    res, used = sort_screen_lfc(counts_df, s.var, args.condit_1, args.condit_2, args.control, args.targets, args.target_col)
    out = run_dir(args.label)
    res = res.reset_index()
    write_tsv(res, out / "guide_lfc.tsv")
    figs: list = []
    exons = None
    if args.gtf and args.gene:
        exons = gtf_gene_features(Path(args.gtf), args.gene, None, "exon", args.exon_source)
    if not args.no_plots:
        fig_sort_lfc(res, out / "volcano.png", out / "guide_enrichment.png", figs, exons, args.gene)
    p = res["lfc.p_value"].to_numpy(float)
    top = res.sort_values("lfc.p_value").head(20)
    idc = res.columns[0]
    sections = [f"Samples matched: condit_1 {used['condit_1']}, condit_2 {used['condit_2']}, control {used['control']}.",
                f"{len(res)} guides" + (f" targeting {args.targets}" if args.targets else "") +
                f"; {int((p < 0.05).sum())} with p < 0.05.", "", "## Top guides", "",
                md_table([idc, "lfc.mean", "lfc.stdev", "lfc.p_value"], top[[idc, "lfc.mean", "lfc.stdev", "lfc.p_value"]].itertuples(index=False))]
    finish(out, "Sorting-screen log fold change", sections,
           {"samples": used, "n_guides": len(res), "n_p_lt_0.05": int((p < 0.05).sum()),
            "median_lfc_mean": float(np.nanmedian(res["lfc.mean"])), "figures": figs})
    return 0


def run_qc(s: Screen, out: Path, args: argparse.Namespace, figs: list) -> Tuple[List[str], dict]:
    pd, np = _pd(), _np()
    cl = args.count_layer
    label = cl or "X"
    qc = sample_qc_table(s, cl)
    set_sample_correlation(s, args.method, cl)
    qc[f"mean_corr_{label}"] = s.obs[f"mean_corr_{label}"].to_numpy()
    qc[f"median_corr_{label}"] = s.obs[f"median_corr_{label}"].to_numpy()
    write_tsv(qc, out / "sample_qc.tsv")
    C = s.obsm[f"corr_{label}"]
    write_tsv(pd.DataFrame(C, index=s.obs.index, columns=s.obs.index).reset_index(), out / f"sample_count_correlation_{args.method}.tsv")
    summary: dict = {"samples": qc.to_dict(orient="records")}
    sections = ["## Sample quality", "", md_table(list(qc.columns), qc.itertuples(index=False))]
    if not args.no_plots:
        fig_count_dist(s, out / "sample_count_dist.png", figs, cl)
        fig_gini(s, out / "sample_count_gini.png", figs, cl)
        fig_heatmap(C, list(s.obs.index), f"Guide count {args.method} correlation", out / "sample_count_correlation.png", figs)
    if args.cond1 is not None and args.cond2 is not None:
        key = _ensure_lognorm(s, "lognorm_counts" if cl is None else f"lognorm_{cl}", cl)
        lf = log_fold_change_reps(s, args.cond1, args.cond2, key, args.rep_col, args.compare_col)
        corr = lf.corr(method=args.method)
        write_tsv(corr.reset_index().rename(columns={"index": "lfc"}), out / "lfc_correlation_all.tsv")
        summary["lfc_correlation_all"] = corr.round(4).to_dict()
        sections += ["", f"## LFC {args.method} correlation ({args.cond1} vs {args.cond2}), all guides", "",
                     md_table(["lfc"] + list(corr.columns), ([i] + list(r) for i, r in corr.iterrows()))]
        if not args.no_plots:
            fig_heatmap(corr.to_numpy(), list(corr.columns), "LFC correlation (all guides)", out / "lfc_correlation_all.png", figs)
            fig_rep_consistency(lf, out / "replicate_consistency.png", figs)
        pos = _read_gene_list(args.pos_ctrl)
        if pos and args.target_col in s.var.columns:
            m = s.var[args.target_col].astype(str).isin(pos).to_numpy()
            if m.sum() >= 3:
                sub = s.subset(guides=m)
                lf2 = log_fold_change_reps(sub, args.cond1, args.cond2, key, args.rep_col, args.compare_col)
                c2 = lf2.corr(method=args.method)
                write_tsv(c2.reset_index().rename(columns={"index": "lfc"}), out / "lfc_correlation_pos_ctrl.tsv")
                summary["lfc_correlation_pos_ctrl"] = c2.round(4).to_dict()
                summary["n_pos_ctrl_guides"] = int(m.sum())
                sections += ["", f"## LFC correlation, {int(m.sum())} positive-control guides", "",
                             md_table(["lfc"] + list(c2.columns), ([i] + list(r) for i, r in c2.iterrows()))]
                if not args.no_plots:
                    fig_heatmap(c2.to_numpy(), list(c2.columns), "LFC correlation (positive controls)", out / "lfc_correlation_pos_ctrl.png", figs)
    if args.outlier_cond_col:
        og = get_outlier_guides(s, args.outlier_cond_col, cl, args.mad_z_thres, args.abs_rpm_thres, args.upstream_compat)
        write_tsv(og, out / "outlier_guides.tsv")
        summary["n_outlier_guides"] = int(len(og))
        sections += ["", f"## Outlier guides (RPM > {args.mad_z_thres} x condition median and > {args.abs_rpm_thres})", "",
                     f"{len(og)} guide-sample outliers." if len(og) else "No outlier guides.", "",
                     md_table(list(og.columns), og.itertuples(index=False)) if len(og) else ""]
    return sections, summary


def cmd_qc(args: argparse.Namespace) -> int:
    s = load_screen(args)
    out = run_dir(args.label)
    figs: list = []
    sections, summary = run_qc(s, out, args, figs)
    summary["figures"] = figs
    finish(out, "Screen sample / guide quality report", [_screen_section(s), ""] + sections, summary)
    return 0


def cmd_outliers(args: argparse.Namespace) -> int:
    s = load_screen(args)
    og = get_outlier_guides(s, args.cond_col, args.count_layer, args.mad_z_thres, args.abs_rpm_thres, args.upstream_compat)
    out = run_dir(args.label)
    write_tsv(og, out / "outlier_guides.tsv")
    finish(out, "Outlier guides", [f"Rule: RPM > {args.mad_z_thres} x median RPM of the guide within its `{args.cond_col}` "
                                   f"condition and RPM > {args.abs_rpm_thres}" + (" (upstream-compat indexing)" if args.upstream_compat else "") + ".",
                                   "", md_table(list(og.columns), og.itertuples(index=False)) if len(og) else "No outlier guides."],
           {"n_outliers": int(len(og)), "outliers": og.head(200).to_dict(orient="records"), "upstream_compat": args.upstream_compat})
    return 0


def cmd_to_mageck(args: argparse.Namespace) -> int:
    s = load_screen(args)
    out = run_dir(args.label)
    path = out / "mageck_count.txt"
    df = to_mageck_input(s, path, args.count_layer, args.sgrna_col, args.target_col, args.sample_prefix)
    print(f"TSV: {path}")
    cmd = f"mageck test -k {path} -t <treatment samples> -c <control samples> -n {out / 'mageck'}"
    finish(out, "MAGeCK input", [f"{len(df)} sgRNAs x {df.shape[1] - 2} samples (rows without a gene dropped: {s.shape[1] - len(df)}).",
                                 "", "Next step (MAGeCK not run here):", "", f"    {cmd}"],
           {"path": str(path), "n_sgrna": int(len(df)), "n_dropped": int(s.shape[1] - len(df)), "columns": list(df.columns),
            "mageck_command": cmd})
    return 0


def cmd_to_excel(args: argparse.Namespace) -> int:
    pd = _pd()
    s = load_screen(args)
    if args.log_norm:
        log_norm(s)
    extra = []
    for spec in args.extra_sheet or []:
        if "=" not in spec:
            raise SystemExit("--extra-sheet expects NAME=path.tsv")
        nm, p = spec.split("=", 1)
        extra.append((nm, _read_any_table(Path(p))))
    out = run_dir(args.label)
    wb = write_screen_to_excel(s, str(out / args.workbook), index=True, include_uns=args.include_uns,
                               guide_rows=not args.sample_rows, extra=extra)
    names = collect_screen_dfs(s, args.include_uns, not args.sample_rows)[1] + [e[0] for e in extra]
    finish(out, "Excel export", [f"Workbook: {wb or 'not written (openpyxl missing)'}", "", "Sheets: " + ", ".join(names)],
           {"workbook": wb, "sheets": names})
    _ = pd
    return 0


def cmd_to_csv(args: argparse.Namespace) -> int:
    s = load_screen(args)
    out = run_dir(args.label)
    d = write_csv_screen(s, out / "screen")
    print(f"Wrote: {d}")
    finish(out, "Screen CSV directory", [f"{d}: screen.X.csv (guides x samples), screen.var.csv, screen.obs.csv, layers/obsm/obsp.",
                                         "Read back with `--screen <dir>` or pt.io.read_csv(X_path=screen.X.csv, ...)."],
           {"screen_dir": str(d), "files": sorted(p.name for p in d.iterdir())})
    return 0


def cmd_read_poolq(args: argparse.Namespace) -> int:
    np = _np()
    s = read_poolq(Path(args.poolq_dir), Path(args.sample_metadata) if args.sample_metadata else None, args.merge_on)
    out = run_dir(args.label)
    files = save_screen(s, out)
    write_tsv(s.var.reset_index(drop=True), out / "guides.tsv")
    sections = [_screen_section(s), ""]
    summ: dict = {"n_samples": s.shape[0], "n_guides": s.shape[1], "total_reads": float(s.X.sum()), **files}
    for k, v in s.uns.items():
        if hasattr(v, "to_csv"):
            write_tsv(v, out / f"poolq_{safe_label(k)}.tsv")
    if "quality_metadata" in s.uns:
        q = s.uns["quality_metadata"]
        sections += ["## PoolQ quality", "", md_table(["metric", "statistic"], q.itertuples(index=False))]
        summ["quality"] = dict(zip(q["metric"], q["statistic"]))
    if "lognorm_counts" in s.layers:
        ours = log_normalize_read_count(s.X)
        ok = np.isfinite(ours) & np.isfinite(s.layers["lognorm_counts"])
        dev = float(np.max(np.abs(ours[ok] - s.layers["lognorm_counts"][ok]))) if ok.any() else float("nan")
        summ["lognorm_max_abs_diff_vs_log2_rpm_plus_1"] = dev
        sections += ["", f"PoolQ lognormalized-counts vs log2(RPM + 1): max |diff| = {dev:.2e}."]
    finish(out, "PoolQ screen", sections, summ)
    return 0


def cmd_annotate_guides(args: argparse.Namespace) -> int:
    pd = _pd()
    g = _read_any_table(Path(args.guides_table))
    if args.id_col not in g.columns:
        raise SystemExit(f"{args.id_col} not in {list(g.columns)}")
    g = annotate_protospacer(g, args.barcode_col)
    direct = None
    if args.direct_pairs:
        dp = _read_any_table(Path(args.direct_pairs))
        direct = dict(zip(dp.iloc[:, 0].astype(str), dp.iloc[:, 1].astype(str)))
    regions = _read_any_table(Path(args.regions)) if args.regions else None
    anns = args.annotations or (list(regions["target"].astype(str)) if regions is not None and "target" in regions.columns else [])
    if "target" not in g.columns or args.force_targets:
        g = add_guide_target_metadata(g, anns, direct, args.id_col)
    n_miss = None
    if regions is not None and args.fasta:
        if "protospacer" not in g.columns:
            raise SystemExit("guides lack 20-nt barcodes / a protospacer column; cannot position them")
        seqs = read_fasta(Path(args.fasta), set(regions["Chromosome"].astype(str)))
        g, n_miss = annotate_guide_position(g, regions, seqs, args.units)
    if "protospacer" in g.columns and args.biochem:
        g = annotate_guide_biochemistry(g)
    out = run_dir(args.label)
    write_tsv(g, out / "guides_annotated.tsv")
    summ = {"n_guides": len(g), "columns": list(g.columns), "n_non_mapping": n_miss,
            "targets": g["target"].value_counts().to_dict() if "target" in g.columns else {}}
    finish(out, "Guide annotation", [f"{len(g)} guides; columns {list(g.columns)}." +
                                     (f" Non-mapping guides: {n_miss}." if n_miss is not None else "")], summ)
    _ = pd
    return 0


def cmd_design_library(args: argparse.Namespace) -> int:
    pd = _pd()
    notes = []
    if args.regions:
        regions = _read_any_table(Path(args.regions))
        ren = {c: c2 for c, c2 in (("chr", "Chromosome"), ("chrom", "Chromosome"), ("start", "Start"), ("end", "End")) if c in regions.columns}
        regions = regions.rename(columns=ren)
        if args.units == "mb":
            regions["Start"] = (regions["Start"].astype(float) * 1e6).round().astype(int)
            regions["End"] = (regions["End"].astype(float) * 1e6).round().astype(int)
        regions = regions[["Chromosome", "Start", "End"]]
    elif args.gtf and args.gene:
        regions = gtf_gene_features(Path(args.gtf), args.gene, args.chrom, args.feature, args.exon_source)
        if regions.empty:
            raise SystemExit(f"no {args.feature} rows for gene_name {args.gene} in {args.gtf}")
        regions = regions[["Chromosome", "Start", "End"]]
    else:
        raise SystemExit("pass --regions or --gtf with --gene")
    regions, ov = resolve_fragments(regions.astype({"Start": int, "End": int}))
    notes.append(f"{len(regions)} regions after removing duplicates and {ov['dropped_engulfing']} engulfing fragment(s); "
                 f"partial overlaps: {len(ov['partial'][FIVE_PRIME])}.")
    if args.widen:
        regions = widen_regions(regions, args.widen)
        notes.append(f"regions widened by {args.widen} bp on both sides.")
    seqs = read_fasta(Path(args.fasta), set(regions["Chromosome"].astype(str)))
    lib = make_guide_library(regions, seqs, args.pam, args.flank)
    lib = scan_for_bsmbi(lib, args.bsmbi_motifs, args.upstream_compat)
    lib = annotate_guide_biochemistry(lib)
    flagged = lib[lib["BsmbI_flag"]]
    out = run_dir(args.label)
    write_tsv(regions, out / "regions.tsv")
    if args.keep_bsmbi:
        final = lib
    else:
        final = lib[~lib["BsmbI_flag"]].reset_index(drop=True)
        write_tsv(flagged.reset_index(drop=True), out / "bsmbi_flagged_guides.tsv")
    write_tsv(final, out / "guide_library.tsv")
    summ = {"n_regions": len(regions), "n_candidates": len(lib), "n_bsmbi_flagged": int(len(flagged)), "n_library": len(final),
            "pam": args.pam, "strand_counts": final["strand"].value_counts().to_dict(),
            "median_gc": float(final["GC_content"].median()) if len(final) else None,
            "n_polyT": int(final["polyT"].sum()) if len(final) else 0}
    finish(out, "sgRNA library design", notes + ["", f"{len(lib)} candidate guides ({args.pam}), {len(flagged)} with BsmbI motifs; "
                                                 f"library {len(final)} guides.", "",
                                                 md_table(["Chromosome", "PAM_loci", "strand", "protospacer", "PAM", "GC_content"],
                                                          final[["Chromosome", "PAM_loci", "strand", "protospacer", "PAM", "GC_content"]].itertuples(index=False), 20)],
           summ)
    _ = pd
    return 0


def cmd_guide_biochem(args: argparse.Namespace) -> int:
    g = _read_any_table(Path(args.guides_table))
    if args.protospacer_col not in g.columns:
        raise SystemExit(f"{args.protospacer_col} not in {list(g.columns)}")
    g = annotate_guide_biochemistry(g, args.protospacer_col)
    out = run_dir(args.label)
    write_tsv(g, out / "guides_biochem.tsv")
    summ = {"n_guides": len(g), "median_gc": float(g["GC_content"].median()),
            **{k: int(g[k].sum()) for k in ("polyT", "polyG", "polyC", "polyA")}}
    finish(out, "Guide biochemistry", [md_table(["metric", "value"], summ.items())], summ)
    return 0


def cmd_pair_distance(args: argparse.Namespace) -> int:
    a = _read_any_table(Path(args.table1)); b = _read_any_table(Path(args.table2))
    res = pair_distance(a, b, args.chrom_col, args.start_col, args.end_col, args.how, args.id1, args.id2, args.max_distance)
    out = run_dir(args.label)
    write_tsv(res, out / "pair_distance.tsv.gz" if len(res) > 200000 else out / "pair_distance.tsv")
    finish(out, "Pairwise feature distance", [f"{len(res)} same-chromosome pairs (how = {args.how})."],
           {"n_pairs": int(len(res)), "how": args.how, "median_distance": float(res["distance"].median()) if len(res) else None})
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    np = _np()
    s = load_screen(args)
    out = run_dir(args.label)
    figs: list = []
    key = log_norm(s)
    sections = [_screen_section(s), "", "Steps: log_norm -> LFC (replicates + aggregate) -> sorting-screen LFC (optional) -> QC -> exports.", ""]
    summ: dict = {"n_samples": s.shape[0], "n_guides": s.shape[1]}
    if args.cond1 is not None and args.cond2 is not None:
        reps = log_fold_change_reps(s, args.cond1, args.cond2, key, args.rep_col, args.compare_col, "lfc", keep_result=True)
        aggs = [log_fold_change_aggregate(s, args.cond1, args.cond2, key, args.rep_col, args.compare_col, "lfc", fn)
                for fn in ("mean", "median", "sd")]
        summ["lfc_replicate_columns"] = list(reps.columns)
        summ["lfc_aggregate_columns"] = aggs
        pos = _read_gene_list(args.pos_ctrl)
        groups = _group_masks(s, args.target_col, pos, _read_gene_list(args.neg_ctrl))
        if groups:
            v = s.var[aggs[1]].to_numpy(float)
            summ["median_lfc_by_group"] = {nm: float(np.nanmedian(v[m])) for nm, m in groups.items() if m.any()}
            sections += ["## LFC by guide group", "", md_table(["group", "guides", "median LFC"],
                                                               [[nm, int(m.sum()), float(np.nanmedian(v[m]))] for nm, m in groups.items() if m.any()]), ""]
        if not args.no_plots:
            fig_lfc_hist(s.var[aggs[1]].to_numpy(float), groups, aggs[1], out / "lfc_aggregate_hist.png", figs)
    if args.condit_1 and args.condit_2 and args.control:
        import pandas as pd  # type: ignore
        cdf = pd.DataFrame(s.layers[key].T, index=s.var.index, columns=[str(c) for c in s.obs.index])
        res, used = sort_screen_lfc(cdf, s.var, args.condit_1, args.condit_2, args.control, args.targets, args.target_col or "target")
        write_tsv(res.reset_index(), out / "sort_guide_lfc.tsv")
        summ["sort_lfc"] = {"samples": used, "n_p_lt_0.05": int((res["lfc.p_value"] < 0.05).sum())}
        if not args.no_plots:
            fig_sort_lfc(res.reset_index(), out / "sort_volcano.png", out / "sort_guide_enrichment.png", figs)
    qc_args = argparse.Namespace(count_layer=None, method=args.method, cond1=args.cond1, cond2=args.cond2, rep_col=args.rep_col,
                                 compare_col=args.compare_col, pos_ctrl=args.pos_ctrl, target_col=args.target_col or "target",
                                 outlier_cond_col=args.outlier_cond_col or (args.compare_col if args.compare_col in s.obs.columns else None),
                                 mad_z_thres=args.mad_z_thres, abs_rpm_thres=args.abs_rpm_thres, upstream_compat=False,
                                 no_plots=args.no_plots)
    qsec, qsum = run_qc(s, out, qc_args, figs)
    sections += qsec
    summ["qc"] = qsum
    write_tsv(s.var.reset_index(), out / "guides_lfc.tsv")
    if args.target_col and args.target_col in s.var.columns:
        mp = out / "mageck_count.txt"
        df = to_mageck_input(s, mp, None, None, args.target_col)
        print(f"TSV: {mp}")
        summ["mageck_input"] = {"path": str(mp), "n_sgrna": int(len(df))}
    if args.excel:
        summ["workbook"] = write_screen_to_excel(s, str(out / "screen.xlsx"))
    summ.update(save_screen(s, out))
    summ["figures"] = figs
    finish(out, "Bulk CRISPR screen analysis", sections, summ)
    return 0


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def _random_seq(rng, n: int) -> str:
    return "".join(rng.choice(list("ACGT"), size=n))


def synthetic_screen(d: Path, seed: int = 3) -> dict:
    """TKO-style dropout screen: 60 genes x 5 guides, T0 + T08A-C + T18A-C, 8 essential genes at LFC -2, one jackpot guide."""
    pd, np = _pd(), _np()
    rng = np.random.default_rng(seed)
    genes = [f"G{i:02d}" for i in range(60)]
    essential = genes[:8]
    gids, gcol, seqs = [], [], []
    for g in genes:
        for k in range(5):
            sq = _random_seq(rng, 20)
            gids.append(f"{g}_{sq}"); gcol.append(g); seqs.append(sq)
    gids[-1] = "odd guide " + seqs[-1]       # space -> underscore in MAGeCK export
    gcol[-2] = ""                            # empty gene -> dropped from MAGeCK export
    n = len(gids)
    base = rng.lognormal(mean=5.5, sigma=0.5, size=n)
    samples = ["T0", "T08A", "T08B", "T08C", "T18A", "T18B", "T18C"]
    ess = np.isin(np.array(gcol), essential)
    ess_lfc = -1.0 - 2.0 * rng.random(n)          # per-guide effect, mean -2
    cols = {}
    for sname in samples:
        mu = base.copy()
        if sname.startswith("T18"):
            mu = mu * np.where(ess, 2.0 ** ess_lfc, 1.0)
        cols[sname] = rng.poisson(mu * 1.0)
    jack = 17
    cols["T08B"][jack] = int(cols["T08B"].sum() * 0.05)   # ~5 % of reads -> RPM ~ 50,000
    df = pd.DataFrame({"GENE_CLONE": gids, "GENE": gcol, **cols})
    p = d / "readcount-synthetic"
    df.to_csv(p, sep="\t", index=False)
    return {"path": p, "df": df, "essential": essential, "jack": gids[jack], "ess_mask": ess, "samples": samples}


def synthetic_sort(d: Path, seed: int = 5) -> dict:
    """Sorting screen: high / low bins and presort, 2 replicates; guides in enh1 enriched +1.5 log2 in high vs low."""
    pd, np = _pd(), _np()
    rng = np.random.default_rng(seed)
    targets = ["enh1"] * 30 + ["enh2"] * 30 + ["NT"] * 40
    n = len(targets)
    base = rng.lognormal(6.0, 0.4, n)
    eff = np.where(np.array(targets) == "enh1", 0.75, 0.0)
    cols = {}
    for r in (1, 2):
        cols[f"presort_r{r}"] = rng.poisson(base)
        cols[f"high_r{r}"] = rng.poisson(base * 2.0 ** eff)
        cols[f"low_r{r}"] = rng.poisson(base * 2.0 ** -eff)
    ids = [f"sg{i:03d}" for i in range(n)]
    df = pd.DataFrame({"guide": ids, "target": targets, **cols})
    p = d / "sort_counts.tsv"
    df.to_csv(p, sep="\t", index=False)
    return {"path": p, "targets": targets}


def synthetic_poolq(d: Path, seed: int = 11) -> dict:
    pd, np = _pd(), _np()
    rng = np.random.default_rng(seed)
    pq = d / "poolq"
    pq.mkdir()
    bcs = [_random_seq(rng, 20) for _ in range(12)]
    ids = [f"BRDN{i:010d}" for i in range(12)]
    conds = ["A03", "A05", "A10"]
    C = rng.poisson(50, size=(12, 3))
    cdf = pd.DataFrame(C, columns=conds)
    cdf.insert(0, "Row Barcode IDs", ids)
    cdf.insert(0, "Row Barcode", bcs)
    cdf.to_csv(pq / "counts.txt", sep="\t", index=False)
    L = np.log2(C / C.sum(0) * 1e6 + 1)
    ldf = pd.DataFrame(L, columns=conds)
    ldf.insert(0, "Row Barcode IDs", ids)
    ldf.insert(0, "Row Barcode", bcs)
    ldf.to_csv(pq / "lognormalized-counts.txt", sep="\t", index=False)
    bc = cdf.rename(columns={"A03": "CGTTTCAA", "A05": "GAGAGGTT", "A10": "GCGGCCAT"})
    bc.to_csv(pq / "barcode-counts.txt", sep="\t", index=False)
    (pq / "unexpected-sequences.txt").write_text("Sequence\tTotal\tCGTTTCAA\tGAGAGGTT\tGCGGCCAT\tPotential IDs\nAACG\t1\t0\t1\t0\t\n")
    with open(pq / "correlation.txt", "wb") as fh:
        fh.write(b"\tA03\tA05\tA10\nA03\t1.00\t0.10\t\xef\xbf\xbd\nA05\t0.10\t1.00\t0.20\nA10\t\xef\xbf\xbd\t0.20\t1.00\n")
    (pq / "quality.txt").write_text(
        "Total reads: 2000\nMatching reads: 1800\n1-base mismatch reads: 10\n\nOverall % match: 90.00\n\n"
        "Read counts for sample barcodes with associated conditions:\n"
        "Barcode\tCondition\tMatched (Construct+Sample Barcode)\tMatched Sample Barcode\t% Match\tNormalized Match\n"
        "CGTTTCAA\tA03\t600\t650\t92.31\t1.000\nGAGAGGTT\tA05\t600\t640\t93.75\t1.000\nGCGGCCAT\tA10\t600\t660\t90.91\t1.000\n\n"
        "Read counts for most common sample barcodes without associated conditions:\nBarcode\tCount\nTTTTTTTT\t12\nGGGGGGGG\t7\n")
    (pq / "runinfo.txt").write_text("PoolQ version: 3.3.2\nPoolQ command-line settings:\n  --counts counts.txt \\\n")
    meta = d / "poolq_conditions.csv"
    pd.DataFrame({"Condition": conds, "time": [0, 7, 14]}).to_csv(meta, index=False)
    return {"dir": pq, "C": C, "L": L, "meta": meta}


def synthetic_genome(d: Path, seed: int = 21) -> dict:
    """One chromosome with a 3-exon gene (plus an engulfing duplicate exon and a planted BsmbI site) and a GTF."""
    np = _np()
    rng = np.random.default_rng(seed)
    seq = list(_random_seq(rng, 6000))
    s = "".join(seq)
    fa = d / "genome.fa"
    with open(fa, "w") as fh:
        fh.write(">chrS test chromosome\n")
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60].lower() if i // 60 % 7 == 0 else s[i:i + 60])  # soft-masked lines
            fh.write("\n")
        fh.write(">chrOther\n" + _random_seq(rng, 300) + "\n")
    exons = [(1000, 1200), (2000, 2300), (4000, 4150)]
    gtf = d / "genes.gtf"
    lines = []
    for (a, b) in exons + [(1990, 2310)]:   # engulfing duplicate of exon 2 (dropped upstream: engulfs another)
        lines.append(f"chrS\tHAVANA\texon\t{a}\t{b}\t.\t+\t.\tgene_id \"ENSG1\"; transcript_id \"T1\"; gene_type \"protein_coding\"; gene_name \"GENEX\";")
    lines.append("chrS\tHAVANA\texon\t5000\t5100\t.\t+\t.\tgene_id \"ENSG2\"; gene_name \"GENEXL\";")  # prefix-match trap
    gtf.write_text("\n".join(lines) + "\n")
    return {"fa": fa, "seq": s, "gtf": gtf, "exons": exons}


def cmd_selftest(args: argparse.Namespace) -> int:
    import shutil
    import tempfile
    pd, np = _pd(), _np()
    from scipy import stats  # type: ignore
    checks = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    def last_run(label):
        return sorted(OUT_ROOT.glob(f"*_{label}*"))[-1]

    made: List[Path] = []
    nop = args.no_plots
    t0 = time.time()
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        W = synthetic_screen(d)

        print("\nscreen object + normalisation")
        s = screen_from_counts_table(pd.read_csv(W["path"], sep="\t"))
        check(s.shape == (7, 300) and list(s.var.columns) == ["GENE"], "make screen: 7 samples x 300 guides; text column GENE -> guide annotation")
        obs = parse_tko_sample_names(s.obs)
        check(list(obs["replicate"]) == ["0", "A", "B", "C", "A", "B", "C"] and list(obs["time"]) == [-1, 8, 8, 8, 18, 18, 18],
              "TKO names: T0 -> replicate '0', time -1; T18B -> B, 18")
        s.obs = obs
        log_norm(s)
        L = s.layers["lognorm_counts"]
        check(np.allclose((2 ** L - 1).sum(axis=1), 1e6), "log_norm: 2^x - 1 sums to 1e6 per sample (RPM)")
        check(np.allclose(L[1], np.log2(s.X[1] / s.X[1].sum() * 1e6 + 1)), "log_norm: log2(X / rowsum * 1e6 + 1)")
        k2 = log_norm(s, read_count_layer="lognorm_counts")
        check(k2 == "lognorm_lognorm_counts", "log_norm with read_count_layer L writes lognorm_L")
        del s.layers[k2]
        col = log_fold_change(s, "T18A", "T08A")
        check(col == "T18A_T08A.lfc" and np.allclose(s.var[col], L[4] - L[1]), "log_fold_change: var['T18A_T08A.lfc'] = lognorm[T18A] - lognorm[T08A]")
        check(fold_change(s, "T18A", "T08A") == "T18A_T08A.fc", "fold_change: '.fc' suffix")
        try:
            log_fold_change(s, "T99", "T08A"); ok = False
        except ValueError:
            ok = True
        check(ok, "log_fold_change: unknown sample raises")
        sub = s.subset(samples=~_match(s.obs["replicate"], 0).to_numpy())
        reps = log_fold_change_reps(sub, 18, 8, rep_col="replicate", compare_col="time")
        check(list(reps.columns) == ["A.18_8.lfc", "B.18_8.lfc", "C.18_8.lfc"], "lfc reps: columns '{rep}.{cond1}_{cond2}.lfc'")
        mc = log_fold_change_aggregate(sub, 18, 8, aggregate_col="replicate", compare_col="time", aggregate_fn="median")
        sd = log_fold_change_aggregate(sub, 18, 8, aggregate_col="replicate", compare_col="time", aggregate_fn="sd")
        check(mc == "18_8.lfc.median" and np.allclose(sub.var[mc], np.median(reps.to_numpy(), axis=1))
              and np.allclose(sub.var[sd], reps.to_numpy().std(axis=1, ddof=0)), "aggregate: median column, sd = np.std(ddof=0)")
        ess = W["ess_mask"]
        med_e, med_o = float(np.median(sub.var[mc][ess])), float(np.median(sub.var[mc][~ess]))
        check(med_e < -1.5 and abs(med_o) < 0.35, f"planted dropout recovered: essential median LFC {med_e:.2f} (planted -2 before renormalisation), others {med_o:.2f}")
        try:
            log_fold_change_reps(s, 18, 8, rep_col="replicate", compare_col="time"); ok = False
        except ValueError:
            ok = True
        check(ok, "lfc reps: replicate '0' with no T18/T08 sample raises (conditions not unique per replicate)")
        s2 = add(s, s)
        check(np.allclose(s2.X, 2 * s.X), "add: technical replicates summed")
        try:
            add(s, s.subset(guides=np.arange(10))); ok = False
        except ValueError:
            ok = True
        check(ok, "add: guide mismatch raises")

        print("\nQC")
        _, _, g = gini_curve(np.array([0, 0, 0, 1.0]))
        check(abs(g - (1 - 1 / 49)) < 1e-12, f"gini: [0,0,0,1] -> 1 - 1/49 = {g:.5f} (50-bin Lorenz curve)")
        _, _, g1 = gini_curve(np.array([5.0] * 10))
        check(abs(g1 + 1) < 1e-12, "gini: constant vector -> -1 (upstream tie handling kept)")
        key = set_sample_correlation(s, "spearman")
        C = s.obsm[key]
        check(key == "corr_X" and np.allclose(np.diag(C), 1) and np.allclose(C, C.T) and "mean_corr_X" in s.obs and "median_corr_X" in s.obs,
              "sample correlation: obsm corr_X symmetric, unit diagonal; obs mean_corr_X / median_corr_X")
        check(abs(s.obs["mean_corr_X"].iloc[0] - C[:, 0].mean()) < 1e-12, "mean_corr includes the diagonal (column mean)")
        set_sample_correlation(s, "pearson", guide_idx=np.arange(50), prefix="sub")
        check("sub_corr_X" in s.obsm, "sample correlation: guide subset + prefix -> sub_corr_X")
        og = get_outlier_guides(s, "time")
        check(len(og) == 1 and og.iloc[0]["GENE_CLONE"] == W["jack"] and og.iloc[0]["sample"] == "T08B" and og.iloc[0]["RPM"] > 10000,
              f"outlier guides: exactly the planted jackpot ({W['jack'][:12]}… in T08B, RPM {og.iloc[0]['RPM']:.0f})" if len(og) else "outlier guides: jackpot found")
        check("X_RPM" in s.layers, "outlier guides: X_RPM layer stored")
        og5 = get_outlier_guides(s, "time", abs_RPM_thres=1e9)
        check(len(og5) == 0, "outlier guides: abs RPM threshold respected")
        ogc = get_outlier_guides(s, "time", upstream_compat=True)
        check(list(ogc.columns) == ["GENE_CLONE", "sample", "RPM", "time"] and not (len(ogc) == 1 and ogc.iloc[0]["GENE_CLONE"] == W["jack"]),
              f"--upstream-compat reproduces the transposed indexing (does not isolate the jackpot; {len(ogc)} rows)")
        qc = sample_qc_table(s)
        check(list(qc.columns) == ["sample", "total_reads", "median_count", "zero_count_guides", "gini"] and qc["gini"].between(-1, 1).all(),
              "sample QC table")

        print("\nexports")
        mg = to_mageck_input(s, d / "m.txt", target_column="GENE")
        check(list(mg.columns[:2]) == ["sgRNA", "gene"] and list(mg.columns[2:]) == W["samples"] and len(mg) == 299
              and any(x.startswith("odd_guide_") for x in mg["sgRNA"]) and str(mg.dtypes.iloc[2]).startswith("int"),
              "to_mageck_input: sgRNA, gene, samples; ints; empty gene row dropped; spaces -> '_'")
        mg2 = to_mageck_input(s, None, count_layer="lognorm_counts", target_column="GENE", sample_prefix="S_")
        check(mg2.shape == (299, 9) and mg2.columns[2] == "S_T0", "to_mageck_input: count layer transposed (upstream shape bug fixed), sample prefix")
        cd = write_csv_screen(s, d / "csvscreen")
        back = read_screen_dir(cd)
        check(np.allclose(back.X, s.X) and np.allclose(back.layers["lognorm_counts"], L) and list(back.obs["time"]) == list(s.obs["time"]),
              "to-csv / read back: X, layers and obs round-trip")
        rc = read_csv_screen(cd / "screen.X.csv", cd / "screen.var.csv", cd / "screen.obs.csv")
        check(np.allclose(rc.X, s.X) and "GENE" in rc.var.columns and "time" in rc.obs.columns, "pt.io.read_csv on the CSV matrix (guides x samples, transposed)")
        wb = write_screen_to_excel(s, str(d / "screen"))
        if wb:
            import openpyxl  # type: ignore
            names = openpyxl.load_workbook(wb, read_only=True).sheetnames
            check(wb.endswith(".xlsx") and names[:5] == ["X", "lognorm_counts", "X_RPM", "guides", "samples"] and "obsm.corr_X" in names,
                  f"to_Excel: sheets {names[:6]}…")
        else:
            check(True, "to_Excel skipped (openpyxl missing)")
        h = write_h5ad(s, d / "s.h5ad")
        if h:
            hs = read_h5ad(h)
            check(np.allclose(hs.X, s.X, rtol=1e-6) and "lognorm_counts" in hs.layers, ".h5ad round-trip")

        print("\nsorting-screen LFC")
        S = synthetic_sort(d)
        ss = screen_from_counts_table(pd.read_csv(S["path"], sep="\t"))
        log_norm(ss)
        cdf = pd.DataFrame(ss.layers["lognorm_counts"].T, index=ss.var.index, columns=ss.obs.index)
        res, used = sort_screen_lfc(cdf, ss.var, "high", "low", "presort")
        check(list(res.columns[-4:]) == ["lfc.mean", "lfc.stdev", "lfc.p_value", "-log10(pval)"] and used["control"] == ["presort_r1", "presort_r2"],
              "sort LFC: columns lfc.mean, lfc.stdev, lfc.p_value, -log10(pval); regex-selected samples")
        e1 = res["target"] == "enh1"
        m1, m0 = float(res.loc[e1, "lfc.mean"].median()), float(res.loc[~e1, "lfc.mean"].median())
        check(abs((m1 - m0) - 1.5) < 0.3 and m0 < 0, f"sort LFC: planted +1.5 recovered relative to null guides ({m1:.2f} vs {m0:.2f}; "
              "RPM normalisation shifts the nulls because 30 % of guides are enriched)")
        a = cdf[["high_r1", "high_r2"]].to_numpy() - cdf[["presort_r1", "presort_r2"]].to_numpy()
        b = cdf[["low_r1", "low_r2"]].to_numpy() - cdf[["presort_r1", "presort_r2"]].to_numpy()
        check(np.allclose(res["lfc.mean"], a.mean(1) - b.mean(1)) and np.allclose(res["lfc.stdev"], (a - b).std(1))
              and np.allclose(res["lfc.p_value"].iloc[:5], [stats.ttest_ind(a[i], b[i])[1] for i in range(5)]),
              "sort LFC: baseline subtraction per replicate, delta of means, sd(ddof=0), ttest_ind p")
        check(float(np.median(res.loc[e1, "lfc.p_value"])) < float(np.median(res.loc[~e1, "lfc.p_value"])),
              "sort LFC: planted guides have smaller p than null guides")
        rt, _ = sort_screen_lfc(cdf, ss.var, "high", "low", "presort", targets=["enh1", "enh2"])
        check(len(rt) == 60, "sort LFC: targets filter")

        print("\nPoolQ reader")
        P = synthetic_poolq(d)
        ps = read_poolq(P["dir"], P["meta"])
        check(ps.shape == (3, 12) and list(ps.var.columns) == ["barcode", "barcode_id"] and np.allclose(ps.X, P["C"].T),
              "read_poolq: counts.txt -> X (samples x guides), var barcode / barcode_id")
        check(np.allclose(ps.layers["lognorm_counts"], log_normalize_read_count(ps.X)), "read_poolq: PoolQ lognormalized-counts == log2(RPM + 1)")
        qm = ps.uns["quality_metadata"]
        check(dict(zip(qm["metric"], qm["statistic"]))["Total reads"] == "2000" and len(ps.uns["quality_bc_readcounts"]) == 3
              and list(ps.uns["quality_common_barcodes"]["barcode"]) == ["TTTTTTTT", "GGGGGGGG"],
              "read_poolq: quality.txt -> metadata / bc_readcounts / common_barcodes")
        check(math.isnan(ps.obsp["poolq_correlation"][0, 2]) and ps.obsp["poolq_correlation"][1, 2] == 0.2 and list(ps.obs["time"]) == [0, 7, 14],
              "read_poolq: correlation.txt with undecodable NaN bytes; sample metadata merged on Condition")

        print("\nguide annotation + library design")
        G = synthetic_genome(d)
        seq = G["seq"]
        region_s, region_e = 1000, 1400
        plus_sp, minus_sp = seq[1100:1120], reverse_complement(seq[1300:1320])
        gt = pd.DataFrame({"barcode": [plus_sp, minus_sp, "A" * 20], "barcode_id": ["enhA_g1", "g2_custom", "enhA_bad"]})
        gt.to_csv(d / "guides_ann.tsv", sep="\t", index=False)
        reg = pd.DataFrame({"target": ["enhA"], "Chromosome": ["chrS"], "Start": [region_s], "End": [region_e]})
        reg.to_csv(d / "regions.tsv", sep="\t", index=False)
        ga = annotate_protospacer(gt)
        ga = add_guide_target_metadata(ga, ["enhA"], {"g2_custom": "enhA"})
        check(list(ga["protospacer"]) == list(gt["barcode"]) and list(ga["target"]) == ["enhA", "enhA", "enhA"],
              "annotate: protospacer = 20-nt barcode; target by direct pair and substring")
        seqs = read_fasta(G["fa"], {"chrS"})
        check(len(seqs["chrS"]) == 6000 and seqs["chrS"] == seq, "read_fasta: soft-masked lines upper-cased, other contigs skipped")
        gp, miss = annotate_guide_position(ga, reg, seqs)
        r0, r1 = gp.iloc[0], gp.iloc[1]
        check(r0["strand"] == "+" and r0["Start"] == 1100 and r0["End"] == 1120 and r0["center"] == 1110,
              "annotate: + strand guide at [1100, 1120), center 1110")
        check(r1["strand"] == "-" and r1["Start"] == 1320 and r1["End"] == 1300 and miss == 1,
              "annotate: - strand guide Start = region end - span0 (1320), End = 1300; poly-A guide non-mapping")
        rg = gtf_gene_features(G["gtf"], "GENEX")
        check(len(rg) == 4, "GTF: exact gene_name match (GENEXL not taken)")
        rr, ov = resolve_fragments(rg[["Chromosome", "Start", "End"]])
        check(len(rr) == 3 and ov["dropped_engulfing"] == 1 and sorted(rr["Start"]) == [1000, 2000, 4000],
              "overlapping fragments: the exon engulfing another is dropped")
        lib = make_guide_library(rr, seqs, "NGG", 10)
        # independent brute force: forward PAMs p in (s, e] with seq[p+1:p+3] == GG, reverse with seq[p-3:p-1] == CC
        exp_f = sum(1 for (a_, b_) in G["exons"] for p in range(a_ + 1, b_ + 1) if seq[p + 1:p + 3] == "GG")
        exp_r = sum(1 for (a_, b_) in G["exons"] for p in range(a_ + 1, b_ + 1) if seq[p - 3:p - 1] == "CC")
        nf, nr = int((lib["strand"] == "+").sum()), int((lib["strand"] == "-").sum())
        check(abs(nf - exp_f) <= 2 and abs(nr - exp_r) <= 2 and nf + nr > 50,
              f"PAM scan: {nf} + / {nr} - guides vs brute force {exp_f} / {exp_r} (duplicates by protospacer dropped)")
        okf = all(seq[int(r.guide_begin):int(r.guide_end)] == r.protospacer and r.PAM[1:] == "GG" and seq[int(r.PAM_loci):int(r.PAM_loci) + 3] == r.PAM
                  for r in lib[lib.strand == "+"].itertuples())
        okr = all(reverse_complement(seq[int(r.guide_end):int(r.guide_begin)]) == r.protospacer and r.PAM[1:] == "GG"
                  for r in lib[lib.strand == "-"].itertuples())
        check(okf and okr, "protospacer + PAM: + strand seq[p-20:p] / seq[p:p+3]; - strand reverse complements; every PAM is NGG")
        okc = all(len(r.guide_context) == 43 for r in lib.itertuples())
        check(okc and set(lib.columns[:12]) >= {"PAM_distance_to_center", "region_center", "guide_context"}, "guide context: 20 + 3 + 2 x 10 flank = 43 nt")
        ctx = pd.DataFrame({"guide_context": ["TCTC" + "A" * 39, "AAATCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "A" * 43]})
        f = scan_for_bsmbi(ctx)
        fc = scan_for_bsmbi(ctx, upstream_compat=True)
        check(list(f["BsmbI_flag"]) == [True, True, False] and list(fc["BsmbI_flag"]) == [False, True, False],
              "BsmbI: TCTC / GAGAC on either strand; --upstream-compat misses a motif at position 0")
        bio = annotate_guide_biochemistry(pd.DataFrame({"protospacer": ["GGGGCCCCAAAATTTTACGT", "ACGTACGTACGTACGTACGT"]}))
        check(list(bio["GC_content"]) == [50.0, 50.0] and list(bio["polyT"]) == [True, False] and bool(bio["polyG"].iloc[0]),
              "guide biochemistry: GC% = (G + C) x 100 / 20; homopolymers")
        pdist = pair_distance(pd.DataFrame({"chr": ["chr1", "chr1", "chr2"], "start": [100, 1000, 5], "end": [200, 1100, 10]}),
                              pd.DataFrame({"chr": ["chr1", "chr3"], "start": [150, 1], "end": [160, 2]}))
        check(sorted(pdist["distance"].tolist()) == [50.0, 850.0] and set(pdist["chr"]) == {"chr1"},
              "pair distance: |start1 - start2| for same-chromosome pairs only")
        ps2 = pair_distance(pd.DataFrame({"chr": ["chr1"], "start": [100], "end": [200]}), pd.DataFrame({"chr": ["chr1"], "start": [150], "end": [300]}), how="shortest")
        check(ps2["distance"].iloc[0] == 0, "pair distance: shortest gap 0 for overlapping features")

        print("\nCLI subcommands")
        ca = ["--counts", str(W["path"]), "--parse-tko-names"]
        rc = main(["make-screen"] + ca + ["--label", "st_pt_make"]); r = last_run("st_pt_make"); made.append(r)
        check(rc == 0 and (r / "screen" / "screen.X.csv").is_file() and (r / "report.md").is_file(), "make-screen: CSV screen dir + report")
        sdir = str(r / "screen")
        rc = main(["add-screens", "--screens", sdir, sdir, "--label", "st_pt_add"]); r = last_run("st_pt_add"); made.append(r)
        check(rc == 0 and json.loads((r / "summary.json").read_text())["total_reads"] == 2 * float(s.X.sum()), "add-screens: totals doubled")
        rc = main(["log-norm", "--screen", sdir, "--label", "st_pt_ln"]); r = last_run("st_pt_ln"); made.append(r)
        check(rc == 0 and (r / "lognorm_counts.tsv").is_file(), "log-norm")
        rc = main(["lfc", "--screen", sdir, "--sample1", "T18A", "--sample2", "T08A", "--label", "st_pt_lfc"] + (["--no-plots"] if nop else []))
        r = last_run("st_pt_lfc"); made.append(r)
        check(rc == 0 and "T18A_T08A.lfc" in pd.read_csv(r / "guides_lfc.tsv", sep="\t").columns, "lfc: T18A_T08A.lfc column")
        rc = main(["lfc-reps", "--screen", sdir, "--cond1", "18", "--cond2", "8", "--rep-col", "replicate", "--compare-col", "time",
                   "--exclude", "replicate=0", "--target-col", "GENE", "--pos-ctrl", ",".join(W["essential"]), "--label", "st_pt_reps"] + (["--no-plots"] if nop else []))
        r = last_run("st_pt_reps"); made.append(r)
        sm = json.loads((r / "summary.json").read_text())
        check(rc == 0 and sm["aggregate_columns"] == ["18_8.lfc.median"] and sm["groups"]["positive controls"]["median_lfc"] < -1.5,
              f"lfc-reps: aggregate + positive-control median {sm['groups'].get('positive controls', {}).get('median_lfc', float('nan')):.2f}")
        rc = main(["sort-lfc", "--counts", str(S["path"]), "--condit-1", "high", "--condit-2", "low", "--control", "presort",
                   "--label", "st_pt_sort"] + (["--no-plots"] if nop else []))
        r = last_run("st_pt_sort"); made.append(r)
        check(rc == 0 and len(pd.read_csv(r / "guide_lfc.tsv", sep="\t")) == 100, "sort-lfc: guide_lfc.tsv")
        rc = main(["qc"] + ca + ["--cond1", "8", "--cond2", "18", "--compare-col", "time", "--exclude", "replicate=0", "--target-col", "GENE",
                                 "--pos-ctrl", ",".join(W["essential"]), "--outlier-cond-col", "time", "--label", "st_pt_qc"] + (["--no-plots"] if nop else []))
        r = last_run("st_pt_qc"); made.append(r)
        sm = json.loads((r / "summary.json").read_text())
        pc = pd.read_csv(r / "lfc_correlation_pos_ctrl.tsv", sep="\t").set_index("lfc")
        al = pd.read_csv(r / "lfc_correlation_all.tsv", sep="\t").set_index("lfc")
        check(rc == 0 and sm["n_outlier_guides"] == 1 and sm["n_pos_ctrl_guides"] == 40 and len(sm["samples"]) == 6,
              "qc: sample table, outlier jackpot, 40 positive-control guides")
        check(pc.iloc[0, 1] > al.iloc[0, 1], f"qc: positive-control LFCs correlate across replicates ({pc.iloc[0, 1]:.2f}) more than all guides ({al.iloc[0, 1]:.2f})")
        if not nop:
            check(len(sm["figures"]) >= 5, f"qc: {len(sm['figures'])} figures")
        rc = main(["outlier-guides", "--screen", sdir, "--cond-col", "time", "--label", "st_pt_out"]); r = last_run("st_pt_out"); made.append(r)
        check(rc == 0 and json.loads((r / "summary.json").read_text())["n_outliers"] == 1, "outlier-guides")
        rc = main(["to-mageck", "--screen", sdir, "--target-col", "GENE", "--label", "st_pt_mageck"]); r = last_run("st_pt_mageck"); made.append(r)
        check(rc == 0 and (r / "mageck_count.txt").read_text().startswith("sgRNA\tgene\tT0"), "to-mageck: MAGeCK count table")
        rc = main(["to-excel", "--screen", sdir, "--log-norm", "--extra-sheet", f"sortlfc={d / 'sort_counts.tsv'}", "--label", "st_pt_xlsx"])
        r = last_run("st_pt_xlsx"); made.append(r)
        check(rc == 0 and "sortlfc" in json.loads((r / "summary.json").read_text())["sheets"], "to-excel: experiment report sheet appended")
        rc = main(["to-csv", "--screen", sdir, "--label", "st_pt_csv"]); r = last_run("st_pt_csv"); made.append(r)
        check(rc == 0 and (r / "screen" / "screen.X.csv").is_file(), "to-csv")
        rc = main(["read-poolq", "--poolq-dir", str(P["dir"]), "--sample-metadata", str(P["meta"]), "--label", "st_pt_poolq"])
        r = last_run("st_pt_poolq"); made.append(r)
        sm = json.loads((r / "summary.json").read_text())
        check(rc == 0 and sm["lognorm_max_abs_diff_vs_log2_rpm_plus_1"] < 1e-9 and sm["quality"]["Matching reads"] == "1800", "read-poolq")
        pd.DataFrame({"barcode_id": ["g2_custom"], "target": ["enhA"]}).to_csv(d / "pairs.tsv", sep="\t", index=False)
        rc = main(["annotate-guides", "--guides-table", str(d / "guides_ann.tsv"), "--regions", str(d / "regions.tsv"), "--fasta", str(G["fa"]),
                   "--direct-pairs", str(d / "pairs.tsv"), "--biochem", "--label", "st_pt_ann"]); r = last_run("st_pt_ann"); made.append(r)
        an = pd.read_csv(r / "guides_annotated.tsv", sep="\t")
        check(rc == 0 and list(an["strand"].fillna("NA")) == ["+", "-", "NA"] and "GC_content" in an.columns, "annotate-guides")
        rc = main(["design-library", "--fasta", str(G["fa"]), "--gtf", str(G["gtf"]), "--gene", "GENEX", "--label", "st_pt_design"])
        r = last_run("st_pt_design"); made.append(r)
        sm = json.loads((r / "summary.json").read_text())
        libt = pd.read_csv(r / "guide_library.tsv", sep="\t")
        check(rc == 0 and sm["n_regions"] == 3 and sm["n_library"] == sm["n_candidates"] - sm["n_bsmbi_flagged"] and not libt["BsmbI_flag"].any()
              and len(libt) == len(lib[~scan_for_bsmbi(lib)["BsmbI_flag"]]), f"design-library: {sm['n_library']} guides after BsmbI filter ({sm['n_bsmbi_flagged']} flagged)")
        rc = main(["guide-biochem", "--guides-table", str(r / "guide_library.tsv"), "--label", "st_pt_bio"]); rb = last_run("st_pt_bio"); made.append(rb)
        check(rc == 0 and (rb / "guides_biochem.tsv").is_file(), "guide-biochem")
        pd.DataFrame({"name": ["a", "b"], "chr": ["chr1", "chr1"], "start": [100, 900], "end": [120, 950]}).to_csv(d / "f1.tsv", sep="\t", index=False)
        pd.DataFrame({"name": ["x"], "chr": ["chr1"], "start": [500], "end": [510]}).to_csv(d / "f2.tsv", sep="\t", index=False)
        rc = main(["pair-distance", "--table1", str(d / "f1.tsv"), "--table2", str(d / "f2.tsv"), "--id1", "name", "--id2", "name",
                   "--how", "midpoint", "--label", "st_pt_dist"]); r = last_run("st_pt_dist"); made.append(r)
        pdd = pd.read_csv(r / "pair_distance.tsv", sep="\t")
        check(rc == 0 and sorted(pdd["distance"]) == [395.0, 420.0], "pair-distance --how midpoint")
        rc = main(["run"] + ca + ["--cond1", "18", "--cond2", "8", "--compare-col", "time", "--exclude", "replicate=0", "--target-col", "GENE",
                                  "--pos-ctrl", ",".join(W["essential"]), "--label", "st_pt_run"] + (["--no-plots"] if nop else []))
        r = last_run("st_pt_run"); made.append(r)
        sm = json.loads((r / "summary.json").read_text())
        check(rc == 0 and sm["median_lfc_by_group"]["positive controls"] < -1.5 and abs(sm["median_lfc_by_group"]["other"]) < 0.35
              and sm["mageck_input"]["n_sgrna"] == 299 and sm["qc"]["n_outlier_guides"] == 1,
              "run: dropout recovered, MAGeCK table, QC with the jackpot outlier")
        rc = main(["run", "--counts", str(S["path"]), "--condit-1", "high", "--condit-2", "low", "--control", "presort",
                   "--target-col", "target", "--label", "st_pt_runsort", "--no-plots"])
        r = last_run("st_pt_runsort"); made.append(r)
        check(rc == 0 and (r / "sort_guide_lfc.tsv").is_file(), "run: sorting-screen branch")
    for p in made:
        shutil.rmtree(p, ignore_errors=True)
    try:
        OUT_ROOT.rmdir()  # only succeeds when the selftest left it empty
    except OSError:
        pass
    ok = all(c for c, _ in checks)
    print(f"\n({time.time() - t0:.1f} s)")
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent perturb-tools", description=__doc__.split("\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def lfc_common(p):
        p.add_argument("--lognorm-key", default="lognorm_counts", help="log-normalised layer (computed if absent)")
        p.add_argument("--count-layer", help="raw count layer to normalise (default X)")

    p = sub.add_parser("make-screen", help="counts table (+ guide/sample tables) -> screen")
    add_screen_args(p)
    p.add_argument("--log-norm", action="store_true", help="also add the lognorm_counts layer")
    p.add_argument("--label", default="screen")
    p.set_defaults(func=cmd_make_screen)

    p = sub.add_parser("add-screens", help="sum technical-replicate screens (identical guides and samples)")
    p.add_argument("--screens", nargs="+", required=True, help=".h5ad files or screen CSV dirs")
    p.add_argument("--label", default="added")
    p.set_defaults(func=cmd_add_screens)

    p = sub.add_parser("log-norm", help="log2(RPM + 1) layer")
    add_screen_args(p)
    p.add_argument("--count-layer", help="normalise this layer instead of X (output lognorm_<layer>)")
    p.add_argument("--label", default="lognorm")
    p.set_defaults(func=cmd_log_norm)

    p = sub.add_parser("lfc", help="LFC of one sample vs another")
    add_screen_args(p)
    lfc_common(p)
    p.add_argument("--sample1", required=True)
    p.add_argument("--sample2", required=True)
    p.add_argument("--suffix", default="lfc", help="column suffix (fc for pt.pp.fold_change)")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="lfc")
    p.set_defaults(func=cmd_lfc)

    p = sub.add_parser("lfc-reps", help="per-replicate LFC + mean/median/sd aggregate")
    add_screen_args(p)
    lfc_common(p)
    p.add_argument("--cond1", required=True)
    p.add_argument("--cond2", required=True)
    p.add_argument("--rep-col", default="replicate")
    p.add_argument("--compare-col", default="sort")
    p.add_argument("--aggregate", nargs="+", default=["median"], choices=["mean", "median", "sd"])
    p.add_argument("--suffix", default="lfc")
    p.add_argument("--method", default="spearman", choices=["pearson", "spearman"], help="replicate LFC correlation")
    p.add_argument("--target-col", help="guide feature with the target gene (for control groups)")
    p.add_argument("--pos-ctrl", help="positive-control genes: file (first column) or comma list")
    p.add_argument("--neg-ctrl", help="negative-control genes: file or comma list")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="lfc_reps")
    p.set_defaults(func=cmd_lfc_reps)

    p = sub.add_parser("sort-lfc", help="sorting-screen delta LFC vs control with t-test")
    add_screen_args(p)
    lfc_common(p)
    p.add_argument("--condit-1", required=True, help="regex selecting condition-1 sample columns (e.g. high)")
    p.add_argument("--condit-2", required=True, help="regex selecting condition-2 sample columns (e.g. low)")
    p.add_argument("--control", required=True, help="regex selecting control columns (e.g. presort)")
    p.add_argument("--targets", nargs="+", help="keep guides whose target is in this list")
    p.add_argument("--target-col", default="target")
    p.add_argument("--gtf", help="GTF for the exon track of the position plot")
    p.add_argument("--gene", help="gene_name for the exon track")
    p.add_argument("--exon-source", default="HAVANA", help="GTF source column filter for the exon track (upstream HAVANA)")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="sort_lfc")
    p.set_defaults(func=cmd_sort_lfc)

    def qc_common(p):
        p.add_argument("--count-layer")
        p.add_argument("--method", default="spearman", choices=["pearson", "spearman"])
        p.add_argument("--cond1", help="LFC correlation: condition 1 value of --compare-col")
        p.add_argument("--cond2")
        p.add_argument("--rep-col", default="replicate")
        p.add_argument("--compare-col", default="time")
        p.add_argument("--target-col", default="target")
        p.add_argument("--pos-ctrl", help="positive-control genes: file (first column) or comma list")
        p.add_argument("--outlier-cond-col", help="obs column defining conditions for outlier guides")
        p.add_argument("--mad-z-thres", type=float, default=5.0)
        p.add_argument("--abs-rpm-thres", type=float, default=10000.0)
        p.add_argument("--upstream-compat", action="store_true", help="outlier guides with the upstream transposed indexing")

    p = sub.add_parser("qc", help="sample quality report")
    add_screen_args(p)
    qc_common(p)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="qc")
    p.set_defaults(func=cmd_qc)

    p = sub.add_parser("outlier-guides", help="jackpot guides per condition")
    add_screen_args(p)
    p.add_argument("--cond-col", required=True)
    p.add_argument("--count-layer")
    p.add_argument("--mad-z-thres", type=float, default=5.0)
    p.add_argument("--abs-rpm-thres", type=float, default=10000.0)
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--label", default="outliers")
    p.set_defaults(func=cmd_outliers)

    p = sub.add_parser("to-mageck", help="MAGeCK count table")
    add_screen_args(p)
    p.add_argument("--count-layer")
    p.add_argument("--sgrna-col", help="guide feature with sgRNA names (default: guide index)")
    p.add_argument("--target-col", default="target_id")
    p.add_argument("--sample-prefix", default="")
    p.add_argument("--label", default="mageck")
    p.set_defaults(func=cmd_to_mageck)

    p = sub.add_parser("to-excel", help="Excel workbook of the screen (+ extra report sheets)")
    add_screen_args(p)
    p.add_argument("--workbook", default="CRISPR_screen.xlsx")
    p.add_argument("--log-norm", action="store_true")
    p.add_argument("--include-uns", action="store_true")
    p.add_argument("--sample-rows", action="store_true", help="samples as rows (default guides as rows)")
    p.add_argument("--extra-sheet", nargs="+", help="NAME=table.tsv sheets (experiment report)")
    p.add_argument("--label", default="excel")
    p.set_defaults(func=cmd_to_excel)

    p = sub.add_parser("to-csv", help="screen CSV directory")
    add_screen_args(p)
    p.add_argument("--label", default="csv")
    p.set_defaults(func=cmd_to_csv)

    p = sub.add_parser("read-poolq", help="PoolQ output directory -> screen")
    p.add_argument("--poolq-dir", required=True)
    p.add_argument("--sample-metadata", help="CSV/TSV merged on --merge-on")
    p.add_argument("--merge-on", default="Condition")
    p.add_argument("--label", default="poolq")
    p.set_defaults(func=cmd_read_poolq)

    p = sub.add_parser("annotate-guides", help="protospacer, target and genomic position of guides")
    p.add_argument("--guides-table", required=True, help="table with barcode and barcode_id columns")
    p.add_argument("--barcode-col", default="barcode")
    p.add_argument("--id-col", default="barcode_id")
    p.add_argument("--regions", help="target regions: target, Chromosome, Start, End")
    p.add_argument("--units", default="bp", choices=["bp", "mb"], help="region coordinate units (upstream metadata: mb)")
    p.add_argument("--fasta", help="reference FASTA (positions)")
    p.add_argument("--annotations", nargs="+", help="target names searched in barcode_id (default: region targets)")
    p.add_argument("--direct-pairs", help="two-column table barcode_id -> target")
    p.add_argument("--force-targets", action="store_true", help="recompute target even if the table has one")
    p.add_argument("--biochem", action="store_true", help="add GC content and homopolymer flags")
    p.add_argument("--label", default="annotate")
    p.set_defaults(func=cmd_annotate_guides)

    p = sub.add_parser("design-library", help="PAM scan of exons / regions -> sgRNA library")
    p.add_argument("--fasta", required=True)
    p.add_argument("--regions", help="BED-like table Chromosome/chr, Start/start, End/end")
    p.add_argument("--units", default="bp", choices=["bp", "mb"])
    p.add_argument("--gtf")
    p.add_argument("--gene")
    p.add_argument("--chrom")
    p.add_argument("--feature", default="exon")
    p.add_argument("--exon-source", help="GTF source filter (e.g. HAVANA)")
    p.add_argument("--pam", default="NGG")
    p.add_argument("--flank", type=int, default=10)
    p.add_argument("--widen", type=int, default=0)
    p.add_argument("--bsmbi-motifs", nargs="+", default=["TCTC", "GAGAC"], help="upstream motifs (the full BsmbI site is CGTCTC)")
    p.add_argument("--keep-bsmbi", action="store_true", help="keep flagged guides (upstream filter_BsmbI=False)")
    p.add_argument("--upstream-compat", action="store_true", help="miss motifs at position 0 like upstream")
    p.add_argument("--label", default="library")
    p.set_defaults(func=cmd_design_library)

    p = sub.add_parser("guide-biochem", help="GC content and homopolymers of protospacers")
    p.add_argument("--guides-table", required=True)
    p.add_argument("--protospacer-col", default="protospacer")
    p.add_argument("--label", default="biochem")
    p.set_defaults(func=cmd_guide_biochem)

    p = sub.add_parser("pair-distance", help="same-chromosome distances between two feature tables")
    p.add_argument("--table1", required=True)
    p.add_argument("--table2", required=True)
    p.add_argument("--id1"); p.add_argument("--id2")
    p.add_argument("--chrom-col", default="chr"); p.add_argument("--start-col", default="start"); p.add_argument("--end-col", default="end")
    p.add_argument("--how", default="start", choices=["start", "midpoint", "shortest"])
    p.add_argument("--max-distance", type=float)
    p.add_argument("--label", default="distance")
    p.set_defaults(func=cmd_pair_distance)

    p = sub.add_parser("run", help="bulk tutorial pipeline end to end")
    add_screen_args(p)
    p.add_argument("--cond1"); p.add_argument("--cond2")
    p.add_argument("--rep-col", default="replicate")
    p.add_argument("--compare-col", default="time")
    p.add_argument("--condit-1"); p.add_argument("--condit-2"); p.add_argument("--control")
    p.add_argument("--targets", nargs="+")
    p.add_argument("--target-col", help="guide feature with target gene (MAGeCK export, control groups)")
    p.add_argument("--pos-ctrl"); p.add_argument("--neg-ctrl")
    p.add_argument("--method", default="spearman", choices=["pearson", "spearman"])
    p.add_argument("--outlier-cond-col")
    p.add_argument("--mad-z-thres", type=float, default=5.0)
    p.add_argument("--abs-rpm-thres", type=float, default=10000.0)
    p.add_argument("--excel", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="bulk_screen")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("selftest", help="synthetic screens, all subcommands, assertions")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    if not any(isinstance(h, logging.FileHandler) for h in logging.getLogger().handlers):
        logp = setup_logging()
        print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
