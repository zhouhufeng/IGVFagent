#!/usr/bin/env python3
"""4th IGVF CRISPR Jamboree (2025) inference benchmark task (port of IGVF-CRISPR/CRISPR-Jamboree_2025).

Port of https://github.com/IGVF-CRISPR/CRISPR-Jamboree_2025 (No LICENSE file;
README + Setup.ipynb + RunTask.ipynb).  The repository drives a checkpointed
run of the IGVF single-cell CRISPR pipeline on two datasets (H9, WTC11):
prepare_user_guide_inference -> inference_sceptre + inference_perturbo ->
mergedResults, then `evaluation_curve.py -i pipeline_outputs/inference_mudata
.h5mu` reports the AUPRC and AUROC with which each method separates the control
pairs, as a 2 x 2 grid of precision-recall and ROC panels.  The task bundle
(jamboree_task.zip: main.nf, input_h9/wtc11.config, evaluation_curve.py, the
MuDatas and the PerTurbo wheel) was hosted on Dropbox and has been deleted, so
the pipeline steps were re-derived from the IGVF CRISPR_Pipeline revision the
bundle ran (github.com/IGVF/CRISPR_Pipeline 1bb2ac2, 2025-02-21: processes
prepare_user_guide_inference / inference_sceptre / inference_perturbo /
mergedResults, bin/prepare_inference.py, inference_sceptre.R,
perturbo_inference.py, export_output_multiple.py / export_output_single.py)
and evaluation_curve.py from the notebook's recorded output (log lines
"Assigning true labels... Performing binary evaluation... Results for
sceptre_p_value / perturbo_p_value", the figure) plus the jamboree-2
perform_binary_evaluation it generalises.  Every step is re-implemented in
Python; no code was copied.  Relationship: port (evaluation_curve.py:
clean-room reconstruction from its outputs).  Not vendored: the H9 / WTC11
inference MuDatas and pairs-to-test sheets (Google Sheets, not public).

Definitions
  prepare         prepare_inference.py: user pairs CSV (guide_id, gene_name,
                  intended_target_name, pair_type) restricted to gene_name in
                  gene.var_names and guide_id in guide.var.guide_id, then
                  gene_name -> gene_id; stored as uns['pairs_to_test'].
  sceptre         inference_sceptre.R: formula ~ log(response_n_nonzero) +
                  log(response_n_umis) (+ non-redundant covariates of the
                  covariate string when formula_object != 'default'),
                  assign_grnas(thresholding, threshold = 1), QC off, side
                  both, union integration, skew-normal approximation; results
                  per (intended_target_name, gene_id), distinct, left-joined to
                  pairs_to_test.  Re-implemented (NB score statistic,
                  conditional randomisation / permutation null, skew-normal
                  fit), as in crispr-jamboree2 / crispr-jamboree3.
  negbinom        extra in-module comparator (the jamboree-2 R module):
                  MASS::glm.nb(y ~ element + offset(log total UMIs)), Wald p.
  perturbo        perturbo_inference.py (package required; tested guides and
                  genes subset, batch_key 'batch', library size
                  total_gene_umis, covariate total_guide_umis, 20 epochs,
                  lr 0.01, batch 128; q_value -> p_value, loc / ln 2 ->
                  log2_fc; merged back to pairs_to_test).
  merge           export_output_multiple.py: sceptre results renamed
                  sceptre_log2_fc / sceptre_p_value, PerTurbo results
                  perturbo_log2_fc / perturbo_p_value, inner-joined on
                  (guide_id, gene_id); per_guide_output.tsv = that joined with
                  guide.var on (intended_target_name, guide_id), guide_id ->
                  guide_id(s), cell_number = cells in the MuData,
                  avg_gene_expression = mean of gene.X over all cells for
                  genes whose symbol == intended_target_name, columns
                  intended_target_name, guide_id(s), intended_target_chr,
                  intended_target_start, intended_target_end, gene_id,
                  sceptre_log2_fc, sceptre_p_value, perturbo_log2_fc,
                  perturbo_p_value, cell_number, avg_gene_expression;
                  per_element_output.tsv = guide ids comma-joined per element.
                  Upstream groups the per-element table by gene_id alone
                  (two elements tested against one gene collapse into one
                  row); the port groups by (intended_target_name, gene_id)
                  and reproduces the gene_id grouping with --upstream-compat.
                  Single-method runs follow export_output_single.py (the
                  other method's columns empty).
  evaluate        evaluation_curve.py: true label from pair_type (positive:
                  positive / direct-targeting controls; negative: negative /
                  non-targeting controls; other pairs dropped), score
                  -log10(p) per <method>_p_value column; sklearn
                  precision_recall_curve + roc_curve, AUPRC = auc(recall,
                  precision) (trapezoid), AUROC; average precision reported
                  too; figure: one row per method, PR left, ROC right.

Subcommands
  prepare    user pairs -> uns['pairs_to_test'] (mudata_inference_input.h5mu)
  infer      sceptre / negbinom / perturbo on the pairs (test_results per method)
  merge      method results -> inference_mudata.h5mu, per_guide_output.tsv,
             per_element_output.tsv
  evaluate   AUPRC / AUROC of every *_p_value column (one or more result
             files, e.g. H9 and WTC11), curves figure, metrics table
  run        prepare -> infer -> merge -> evaluate
  selftest   synthetic inference MuData + pairs sheet with planted
             positive / negative controls; all subcommands asserted

Output: Docs/CRISPRJamboree4/<timestamp>_<label>/.

Usage:
    igvfagent crispr-jamboree4 prepare --mudata mudata_guide_assignment.h5mu --pairs pairs_to_test_h9.csv
    igvfagent crispr-jamboree4 infer --mudata mudata_inference_input.h5mu --methods sceptre negbinom
    igvfagent crispr-jamboree4 merge --mudata mudata_inference_input.h5mu --sceptre-results test_results.sceptre.tsv --perturbo-results perturbo.tsv
    igvfagent crispr-jamboree4 evaluate --results H9=per_guide_output.tsv WTC11=wtc11/per_guide_output.tsv
    igvfagent crispr-jamboree4 run --mudata mudata_guide_assignment.h5mu --pairs pairs_to_test_h9.csv --label h9
    igvfagent crispr-jamboree4 selftest --no-plots
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import os
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "CRISPRJamboree4"

UPSTREAM_REPO = "IGVF-CRISPR/CRISPR-Jamboree_2025"
UPSTREAM_COMMIT = "ca56b3249373cba508c32f8269d884fb16c768b7"
PIPELINE_REVISION = "IGVF/CRISPR_Pipeline@1bb2ac269bf4b69fdeb8e7e03118fcf0f23ccbee"
METHODS = ["sceptre", "negbinom", "perturbo"]
POS_TYPES = {"positive_control", "positive control", "poscontrol", "positive", "direct targeting", "direct_targeting", "cis"}
NEG_TYPES = {"negative_control", "negative control", "negcontrol", "negative", "non-targeting", "non_targeting",
             "non-targeting_negative_control", "targeting_negative_control", "safe-targeting", "safe_targeting"}
PER_GUIDE_COLS = ["intended_target_name", "guide_id(s)", "intended_target_chr", "intended_target_start", "intended_target_end", "gene_id",
                  "sceptre_log2_fc", "sceptre_p_value", "perturbo_log2_fc", "perturbo_p_value", "cell_number", "avg_gene_expression"]

INK, INK2, AXIS, SURFACE = "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"
BLUE, ORANGE, GREEN, RED, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#e34948", "#96a0b3"

log = logging.getLogger("crispr_jamboree4")


# ---------------------------------------------------------------------------
# Plumbing and MuData I/O (h5py + anndata; the `mudata` package is optional)
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"crispr_jamboree4_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def run_dir(label: str) -> Path:
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _pd():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit("this subcommand needs pandas + numpy + scipy + anndata + h5py") from e


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


def _save(fig, path: Path) -> Path:
    fig.patch.set_facecolor(SURFACE)
    fig.savefig(str(path), bbox_inches="tight", dpi=130)
    import matplotlib.pyplot as plt  # type: ignore
    plt.close(fig)
    print(f"Figure: {path}")
    return path


def md_table(headers: "list[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 60) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for i, r in enumerate(rows):
        if i >= max_rows:
            lines.append("| … |" + " |" * (len(headers) - 1))
            break
        lines.append("| " + " | ".join(str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(lines)


def _fmt(x, nd: int = 3) -> str:
    try:
        if x is None or (isinstance(x, float) and not math.isfinite(x)):
            return "nan"
        v = float(x)
        if v != 0 and abs(v) < 10 ** -nd:
            return f"{v:.2e}"
        return f"{v:.{nd}f}"
    except (TypeError, ValueError):
        return str(x)


def write_tsv(df, path: Path, label: str = "TSV") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    sep = "," if str(path).endswith(".csv") else "\t"
    df.to_csv(path, sep=sep, index=False, compression="gzip" if str(path).endswith(".gz") else None)
    print(f"{'CSV' if sep == ',' else label}: {path}")
    return path


def write_json(obj, path: Path) -> Path:
    path.write_text(json.dumps(obj, indent=2, default=_json_default))
    print(f"JSON: {path}")
    return path


def _json_default(o):
    np = _np()
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    return str(o)


def write_report(d: Path, title: str, sections: "list[str]", summary: dict) -> None:
    rp = d / "report.md"
    rp.write_text(f"# {title}\n\n" + "\n\n".join(sections) + "\n")
    print(f"Report: {rp}")
    write_json(summary, d / "summary.json")


def read_table(path: Path):
    """TSV / CSV / XLSX guide metadata or pairs tables."""
    pd = _pd()
    s = str(path).lower()
    if s.endswith((".xlsx", ".xls")):
        return pd.read_excel(path)
    if s.endswith((".csv", ".csv.gz")):
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def _elem_io():
    try:
        from anndata.io import read_elem, write_elem  # type: ignore
    except ImportError:
        from anndata.experimental import read_elem, write_elem  # type: ignore
    return read_elem, write_elem


class MuDataLite:
    """Modalities (name -> AnnData) + global obs + uns, read / written as .h5mu."""

    def __init__(self, mods: dict, obs=None, uns: Optional[dict] = None) -> None:
        pd = _pd()
        self.mods = dict(mods)
        first = next(iter(self.mods.values()))
        self.obs = obs if obs is not None else pd.DataFrame(index=first.obs_names.copy())
        self.uns = dict(uns or {})

    def __getitem__(self, k: str):
        return self.mods[k]

    def __contains__(self, k: str) -> bool:
        return k in self.mods

    @property
    def gene(self):
        return self.mods["gene"]

    @property
    def guide(self):
        return self.mods["guide"]

    def subset_cells(self, names) -> "MuDataLite":
        names = list(names)
        mods = {k: v[names].copy() for k, v in self.mods.items()}
        obs = self.obs.loc[[n for n in names if n in self.obs.index]] if len(self.obs.columns) else None
        return MuDataLite(mods, obs, self.uns)


def read_h5mu(path: Path) -> MuDataLite:
    import h5py  # type: ignore
    read_elem, _ = _elem_io()
    with h5py.File(str(path), "r") as h:
        mods = {k: read_elem(h[f"mod/{k}"]) for k in h["mod"].keys()}
        obs = read_elem(h["obs"]) if "obs" in h else None
        uns = read_elem(h["uns"]) if "uns" in h else {}
    if "gene" not in mods and "rna" in mods:
        mods["gene"] = mods.pop("rna")
    if "guide" not in mods and "gRNA" in mods:
        mods["guide"] = mods.pop("gRNA")
    return MuDataLite(mods, obs, uns)


def write_h5mu(md: MuDataLite, path: Path) -> Path:
    import h5py  # type: ignore
    np, pd = _np(), _pd()
    _, write_elem = _elem_io()
    path.parent.mkdir(parents=True, exist_ok=True)
    obs = md.obs.copy() if md.obs is not None else pd.DataFrame(index=md.gene.obs_names)
    with h5py.File(str(path), "w") as h:
        h.attrs["encoding-type"] = "MuData"
        h.attrs["encoding-version"] = "0.1.0"
        for k, a in md.mods.items():
            write_elem(h, f"mod/{k}", a)
        h["mod"].attrs["mod-order"] = np.array(list(md.mods), dtype=object).astype("S")
        write_elem(h, "obs", obs)
        var_index = pd.Index(np.concatenate([np.asarray(a.var_names) for a in md.mods.values()]))
        write_elem(h, "var", pd.DataFrame(index=var_index.astype(str)))
        write_elem(h, "obsm", {k: np.isin(np.asarray(obs.index), np.asarray(a.obs_names)).astype(np.int32) for k, a in md.mods.items()})
        write_elem(h, "varm", {k: np.isin(np.asarray(var_index), np.asarray(a.var_names)).astype(np.int32) for k, a in md.mods.items()})
        uns = {}
        for k, v in md.uns.items():
            uns[k] = v.reset_index(drop=True) if isinstance(v, pd.DataFrame) else v
        write_elem(h, "uns", uns)
    print(f"Wrote: {path}")
    return path


def uns_df(md: MuDataLite, key: str, required: bool = True):
    pd = _pd()
    v = md.uns.get(key)
    if v is None:
        if required:
            raise SystemExit(f"'{key}' not found in mdata.uns (available: {sorted(md.uns)})")
        return None
    if isinstance(v, pd.DataFrame):
        return v.reset_index(drop=True).copy()
    if isinstance(v, dict):
        return pd.DataFrame({k: (list(x) if not hasattr(x, "shape") else x) for k, x in v.items()})
    return pd.DataFrame(v)


def _uns_scalar(adata, key: str, default: str = "") -> str:
    v = adata.uns.get(key, default)
    if hasattr(v, "__len__") and not isinstance(v, str):
        v = v[0] if len(v) else default
    return str(v)


def moi_of(md: MuDataLite, override: Optional[str] = None) -> str:
    return (override or _uns_scalar(md.guide, "moi", "high")).lower()


def truthy_str(v) -> bool:
    return str(v).strip().upper() in {"TRUE", "T", "1", "YES"}


def is_targeting(guide_var):
    np = _np()
    if "targeting" in guide_var.columns:
        return np.array([truthy_str(v) for v in guide_var["targeting"]])
    names = guide_var["intended_target_name"].astype(str).str.lower()
    return ~(names.str.contains("non-targeting") | names.str.contains("non_targeting") | names.str.contains("safe"))


def dense(m):
    np = _np()
    return np.asarray(m.todense()) if hasattr(m, "todense") else np.asarray(m)


def assignment_matrix(guide, layer: str = "guide_assignment"):
    np = _np()
    if layer in guide.layers:
        return dense(guide.layers[layer]) > 0
    return dense(guide.X) > 0


def element_indicator(guide, target: str, layer: str = "guide_assignment"):
    """max over the element's gRNAs of the assignment layer (the union strategy)."""
    np = _np()
    cols = np.where(guide.var["intended_target_name"].astype(str).values == str(target))[0]
    if cols.size == 0:
        return np.zeros(guide.n_obs, dtype=bool)
    A = guide.layers[layer] if layer in guide.layers else guide.X
    sub = A[:, cols]
    return np.asarray((sub.max(axis=1).todense() if hasattr(sub, "todense") else sub.max(axis=1))).ravel() > 0


def nt_cells(guide, layer: str = "guide_assignment"):
    np = _np()
    nt = ~is_targeting(guide.var)
    A = assignment_matrix(guide, layer)
    return A[:, nt].any(axis=1) if nt.any() else np.zeros(guide.n_obs, dtype=bool)


def gene_counts(gene, gene_ids):
    """cells x len(gene_ids) dense count matrix (float)."""
    np = _np()
    idx = gene.var_names.get_indexer(list(gene_ids))
    if (idx < 0).any():
        missing = [g for g, i in zip(gene_ids, idx) if i < 0][:5]
        raise SystemExit(f"gene_id(s) not in gene modality: {missing}")
    return dense(gene.X[:, idx]).astype(float)


def total_umis(gene):
    np = _np()
    if "total_gene_umis" in gene.obs.columns:
        return np.asarray(gene.obs["total_gene_umis"], dtype=float)
    return np.asarray(gene.X.sum(axis=1)).ravel().astype(float)


# ---------------------------------------------------------------------------
# Statistics: NB GLM, sceptre-style test, metrics (shared with crispr-jamboree2/3)
# ---------------------------------------------------------------------------

def glm_log_irls(y, X, offset=None, theta: float = float("inf"), weights=None, maxit: int = 50, tol: float = 1e-8, beta=None):
    """Log-link GLM (Poisson if theta = inf, NB2 otherwise) by IRLS. Returns (beta, cov, mu)."""
    np = _np()
    n, p = X.shape
    off = np.zeros(n) if offset is None else np.asarray(offset, dtype=float)
    pw = np.ones(n) if weights is None else np.asarray(weights, dtype=float)
    if beta is None:
        mu = y + 0.1
        eta = np.log(mu)
    else:
        eta = X @ beta + off
        mu = np.exp(np.clip(eta, -30, 30))
    dev_old = np.inf
    b = np.zeros(p) if beta is None else np.asarray(beta, float)
    for _ in range(maxit):
        vr = 1.0 + (mu / theta if math.isfinite(theta) else 0.0)
        w = pw * mu / vr
        z = eta - off + (y - mu) / mu
        WX = X * w[:, None]
        A = X.T @ WX
        try:
            b = np.linalg.solve(A + 1e-10 * np.eye(p), WX.T @ z)
        except np.linalg.LinAlgError:
            b = np.linalg.lstsq(A, WX.T @ z, rcond=None)[0]
        eta = X @ b + off
        mu = np.exp(np.clip(eta, -30, 30))
        if math.isfinite(theta):
            dev = float(np.sum(pw * (-(y * np.log(mu + 1e-300)) + (y + theta) * np.log(mu + theta))))
        else:
            dev = float(np.sum(pw * (mu - y * np.log(mu + 1e-300))))
        if abs(dev - dev_old) < tol * (abs(dev) + 0.1):
            break
        dev_old = dev
    vr = 1.0 + (mu / theta if math.isfinite(theta) else 0.0)
    w = pw * mu / vr
    try:
        cov = np.linalg.inv(X.T @ (X * w[:, None]))
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(X.T @ (X * w[:, None]))
    return b, cov, mu


def theta_ml(y, mu, limit: int = 25, eps: float = 1e-8) -> float:
    """MASS::theta.ml — Newton iterations on the NB2 profile score for theta."""
    np = _np()
    from scipy.special import digamma, polygamma  # type: ignore
    n = len(y)
    denom = np.sum((y / mu - 1.0) ** 2)
    t = n / denom if denom > 0 else 1e4
    t = min(max(t, 1e-3), 1e6)
    for _ in range(limit):
        score = np.sum(digamma(t + y) - digamma(t) + np.log(t) + 1 - np.log(t + mu) - (y + t) / (mu + t))
        info = np.sum(-polygamma(1, t + y) + polygamma(1, t) - 1 / t + 2 / (mu + t) - (y + t) / (mu + t) ** 2)
        if not np.isfinite(score) or not np.isfinite(info) or info <= 0:
            break
        delta = score / info
        t_new = t + delta
        if t_new <= 0:
            t_new = t / 10
        if abs(t_new - t) < eps * t:
            t = t_new
            break
        t = min(t_new, 1e8)
    return float(t)


def glm_nb(y, X, offset=None, maxit: int = 25):
    """MASS::glm.nb: alternate theta.ml and IRLS. Returns dict(beta, se, z, p, theta, mu)."""
    np = _np()
    from scipy.stats import norm  # type: ignore
    b, cov, mu = glm_log_irls(y, X, offset)
    theta = theta_ml(y, mu)
    for _ in range(maxit):
        b_new, cov, mu = glm_log_irls(y, X, offset, theta=theta, beta=b)
        th_new = theta_ml(y, mu)
        done = np.max(np.abs(b_new - b)) < 1e-7 and abs(th_new - theta) < 1e-5 * theta
        b, theta = b_new, th_new
        if done:
            break
    b, cov, mu = glm_log_irls(y, X, offset, theta=theta, beta=b)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    z = b / np.where(se > 0, se, np.nan)
    return {"beta": b, "se": se, "z": z, "p": 2 * norm.sf(np.abs(z)), "theta": theta, "mu": mu}


def logistic_irls(x, Z, maxit: int = 50):
    np = _np()
    n, p = Z.shape
    b = np.zeros(p)
    xm = min(max(x.mean(), 1e-4), 1 - 1e-4)
    b[0] = math.log(xm / (1 - xm))
    for _ in range(maxit):
        eta = Z @ b
        pr = 1 / (1 + np.exp(-np.clip(eta, -30, 30)))
        w = np.clip(pr * (1 - pr), 1e-10, None)
        zz = eta + (x - pr) / w
        b_new = np.linalg.solve((Z * w[:, None]).T @ Z + 1e-8 * np.eye(p), (Z * w[:, None]).T @ zz)
        if np.max(np.abs(b_new - b)) < 1e-8:
            b = b_new
            break
        b = b_new
    return 1 / (1 + np.exp(-np.clip(Z @ b, -30, 30)))


def binary_metrics(labels, scores) -> dict:
    """sklearn precision_recall_curve / roc_curve + auc(), as the jamboree notebooks."""
    np = _np()
    labels = np.asarray(labels, int)
    scores = np.asarray(scores, float)
    ok = np.isfinite(scores)
    labels, scores = labels[ok], scores[ok]
    out = {"n": int(len(labels)), "n_pos": int(labels.sum()), "n_neg": int((1 - labels).sum()),
           "auprc": float("nan"), "auroc": float("nan"), "average_precision": float("nan")}
    if out["n_pos"] == 0 or out["n_neg"] == 0:
        return out
    from sklearn.metrics import auc, average_precision_score, precision_recall_curve, roc_curve  # type: ignore
    pre, rec, _ = precision_recall_curve(labels, scores)
    fpr, tpr, _ = roc_curve(labels, scores)
    out.update({"auprc": float(auc(rec, pre)), "auroc": float(auc(fpr, tpr)),
                "average_precision": float(average_precision_score(labels, scores)),
                "_curves": {"precision": pre, "recall": rec, "fpr": fpr, "tpr": tpr}})
    return out


def sceptre_covariates(gene, obs_extra=None):
    """Design matrix for ~ log(response_n_nonzero) + log(response_n_umis) [+ factors]."""
    np, pd = _np(), _pd()
    X = gene.X
    n_umis = np.asarray(X.sum(axis=1)).ravel().astype(float)
    n_nonzero = np.asarray((X > 0).sum(axis=1)).ravel().astype(float)
    cols = [np.ones(gene.n_obs), np.log(np.maximum(n_nonzero, 1)), np.log(np.maximum(n_umis, 1))]
    names = ["(Intercept)", "log(response_n_nonzero)", "log(response_n_umis)"]
    if obs_extra is not None and len(obs_extra.columns):
        for c in obs_extra.columns:
            dm = pd.get_dummies(obs_extra[c].astype(str), prefix=c, drop_first=True).astype(float)
            for cc in dm.columns:
                cols.append(dm[cc].values)
                names.append(cc)
    Z = np.column_stack(cols)
    Z[:, 1:3] = Z[:, 1:3] - Z[:, 1:3].mean(axis=0)
    return Z, names


def nb_score_stats(Xt, y, mu, theta, Z):
    """Score z-statistics for adding each row of Xt (B x n, 0/1) to the NB null model."""
    np = _np()
    w = mu / (1 + mu / theta)
    r = (y - mu) / (1 + mu / theta)
    U = Xt @ r
    ZW = Z * w[:, None]
    A_inv = np.linalg.pinv(Z.T @ ZW)
    XW = Xt * w[None, :]
    xwx = XW.sum(axis=1) if set(np.unique(Xt)) <= {0.0, 1.0} else np.einsum("ij,ij->i", XW, Xt)
    XWZ = XW @ Z
    info = xwx - np.einsum("ij,jk,ik->i", XWZ, A_inv, XWZ)
    return U / np.sqrt(np.clip(info, 1e-12, None))


def skewnorm_pvalue(z_obs: float, null, side: str) -> "tuple[float, float, bool]":
    """(p_skew_normal, p_empirical, fit_ok) as sceptre's resampling_approximation = skew_normal."""
    np = _np()
    from scipy import stats  # type: ignore
    null = np.asarray(null, float)
    B = len(null)
    emp_l = (1 + np.sum(null <= z_obs)) / (B + 1)
    emp_r = (1 + np.sum(null >= z_obs)) / (B + 1)
    emp = {"left": emp_l, "right": emp_r, "both": min(1.0, 2 * min(emp_l, emp_r))}[side]
    try:
        a, loc, scale = stats.skewnorm.fit(null)
        ok = np.isfinite([a, loc, scale]).all() and scale > 0
        # goodness of fit: KS on the null draws (sceptre falls back to the empirical p when the fit is poor)
        ks = stats.kstest(null, "skewnorm", args=(a, loc, scale)).pvalue
        ok = ok and ks > 0.01
        cdf = stats.skewnorm.cdf(z_obs, a, loc, scale)
        sf = stats.skewnorm.sf(z_obs, a, loc, scale)
        p = {"left": cdf, "right": sf, "both": min(1.0, 2 * min(cdf, sf))}[side]
        return (float(p) if ok else float(emp)), float(emp), bool(ok)
    except Exception:
        return float(emp), float(emp), False


def sceptre_test_pairs(md: MuDataLite, pairs, side: str = "both", mechanism: str = "default", B: int = 499,
                       moi: Optional[str] = None, seed: int = 4, covariate_obs=None, layer: str = "guide_assignment"):
    """p_value + log2_fc per (intended_target_name, gene_id) pair."""
    np, pd = _np(), _pd()
    rng = np.random.default_rng(seed)
    gene, guide = md.gene, md.guide
    moi = moi_of(md, moi)
    mech = mechanism if mechanism != "default" else ("crt" if moi == "high" else "permutations")
    Z_all, _ = sceptre_covariates(gene, covariate_obs)
    ntc = nt_cells(guide, layer)
    out = []
    null_cache: dict = {}
    targets = pairs["intended_target_name"].astype(str).unique()
    for target in targets:
        trt = element_indicator(guide, target, layer)
        if moi == "low":
            keep = trt | (ntc & ~trt)
        else:
            keep = np.ones(gene.n_obs, dtype=bool)
        x = trt[keep].astype(float)
        Z = Z_all[keep]
        sub = pairs[pairs["intended_target_name"].astype(str) == target]
        if x.sum() == 0 or x.sum() == len(x):
            for _, r in sub.iterrows():
                out.append({"intended_target_name": target, "gene_id": r["gene_id"], "p_value": np.nan, "log2_fc": np.nan,
                            "n_treatment": int(x.sum()), "n_control": int(len(x) - x.sum())})
            continue
        if mech == "crt":
            pr = logistic_irls(x, Z)
            Xt = (rng.random((B, len(x))) < pr[None, :]).astype(float)
        else:
            Xt = np.array([rng.permutation(x) for _ in range(B)])
        Y = gene_counts(gene, sub["gene_id"].tolist())[keep]
        for j, (_, r) in enumerate(sub.iterrows()):
            y = Y[:, j]
            key = (r["gene_id"], moi == "low" and target)
            if key not in null_cache:
                if y.sum() == 0:
                    null_cache[key] = None
                else:
                    b, _, mu = glm_log_irls(y, Z)
                    th = theta_ml(y, mu)
                    b, _, mu = glm_log_irls(y, Z, theta=th, beta=b)
                    null_cache[key] = (mu, th, b)
            fit = null_cache[key]
            if fit is None:
                out.append({"intended_target_name": target, "gene_id": r["gene_id"], "p_value": np.nan, "log2_fc": np.nan,
                            "n_treatment": int(x.sum()), "n_control": int(len(x) - x.sum())})
                continue
            mu, th, b = fit
            z_obs = float(nb_score_stats(x[None, :], y, mu, th, Z)[0])
            z_null = nb_score_stats(Xt, y, mu, th, Z)
            p, p_emp, ok = skewnorm_pvalue(z_obs, z_null, side)
            b1, _, _ = glm_log_irls(y, np.column_stack([Z, x]), theta=th, beta=np.concatenate([b, [0.0]]))
            out.append({"intended_target_name": target, "gene_id": r["gene_id"], "p_value": p, "log2_fc": float(b1[-1] / math.log(2)),
                        "z": z_obs, "p_value_empirical": p_emp, "skew_normal_fit": ok,
                        "n_treatment": int(x.sum()), "n_control": int(len(x) - x.sum())})
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Pipeline helpers shared with crispr-jamboree3 (prepare_inference.py, inference_sceptre.R, perturbo_inference.py)
# ---------------------------------------------------------------------------

def strip_version(ids):
    return [str(x).split(".")[0] for x in ids]


def restrict_pairs(pairs, md: MuDataLite):
    genes = set(map(str, md.gene.var_names))
    guides = set(map(str, md.guide.var["guide_id"])) if "guide_id" in md.guide.var.columns else set(map(str, md.guide.var_names))
    inc1 = set(pairs["gene_name"].astype(str)) & genes
    inc2 = set(pairs["guide_id"].astype(str)) & guides
    print(f"Number of genes in common: {len(inc1)}")
    print(f"Number of guides in common: {len(inc2)}")
    sub = pairs[pairs["gene_name"].astype(str).isin(inc1) & pairs["guide_id"].astype(str).isin(inc2)].copy()
    if sub.empty:
        raise SystemExit("The subset of guide_inference is empty after filtering. Please check your input data.")
    return sub.rename(columns={"gene_name": "gene_id"}).reset_index(drop=True)


def identify_non_redundant_covariates(data, cov_string: str) -> str:
    covs = [c for c in re.split(r"\s*\+\s*", re.sub(r"\s*\+\s*$", "", cov_string.strip())) if c]
    if len(covs) <= 1:
        return cov_string

    def same_levels(a, b):
        if a not in data.columns or b not in data.columns:
            raise SystemExit(f"Columns '{a}' or '{b}' do not exist in the data")
        la, lb = data[a].unique(), data[b].unique()
        if len(la) == 0 or len(lb) == 0:
            return False
        return len(la) == len(lb) and all(len(set(data.loc[data[a] == l, b])) == 1 for l in la)
    groups: list = []
    for i in range(len(covs) - 1):
        for j in range(i + 1, len(covs)):
            if same_levels(covs[i], covs[j]):
                grp = [covs[i], covs[j]]
                for gg in groups:
                    if set(grp) & set(gg):
                        gg.extend(x for x in grp if x not in gg)
                        break
                else:
                    groups.append(grp)
    keep = [g[0] for g in groups] + [c for c in covs if not any(c in g for g in groups)]
    keep = [c for c in keep if data[c].nunique() > 1]
    return " + ".join(keep)


def sceptre_inference(md: MuDataLite, side: str = "both", formula: str = "default", cov_string: str = "", B: int = 499,
                      seed: int = 4, moi: Optional[str] = None):
    pairs = uns_df(md, "pairs_to_test")
    pairs["gene_id"] = pairs["gene_id"].astype(str)
    pairs["intended_target_name"] = pairs["intended_target_name"].astype(str)
    cov_obs = None
    used = ""
    if formula != "default" and cov_string.strip():
        used = identify_non_redundant_covariates(md.obs, cov_string)
        cols = [c for c in re.split(r"\s*\+\s*", used) if c]
        cov_obs = md.obs[cols] if cols else None
    uniq = pairs[["intended_target_name", "gene_id"]].drop_duplicates()
    res = sceptre_test_pairs(md, uniq, side=side, B=B, moi=moi, seed=seed, covariate_obs=cov_obs)
    res = res[["intended_target_name", "gene_id", "p_value", "log2_fc"]].drop_duplicates(["gene_id", "intended_target_name"])
    return pairs.merge(res, how="left", on=["intended_target_name", "gene_id"]), used


def infer_perturbo(md_path: Path, out_path: Path):
    try:
        import mudata  # type: ignore  # noqa: F401
        import perturbo  # type: ignore  # noqa: F401
    except ImportError:
        print("perturbo: package not importable (pip install perturbo-0.0.1-py3-none-any.whl mudata); not re-implemented here.")
        print(f"Command: perturbo_inference.py {md_path} {out_path}")
        return None
    import mudata as mdmod  # type: ignore
    import numpy as np  # type: ignore
    import pandas as pd  # type: ignore
    import perturbo  # type: ignore
    m = mdmod.read_h5mu(str(md_path))
    m["gene"].obs = (m.obs.join(m["gene"].obs.drop(columns=m.obs.columns, errors="ignore"))
                     .join(m["guide"].obs.drop(columns=m.obs.columns.union(m["gene"].obs.columns), errors="ignore"))
                     .assign(log1p_total_guide_umis=lambda x: np.log1p(x["total_guide_umis"])))
    m["guide"].X = m["guide"].layers["guide_assignment"]
    ptt = pd.DataFrame(m.uns["pairs_to_test"])
    names = sorted(pd.unique(ptt["intended_target_name"]))
    m.uns["intended_target_names"] = names
    agg = ptt.assign(value=1).groupby(["gene_id", "intended_target_name"]).agg(value=("value", "max")).reset_index()
    m["gene"].varm["intended_targets"] = agg.pivot(index="gene_id", columns="intended_target_name", values="value").reindex(m["gene"].var_names).fillna(0)
    m["guide"].varm["intended_targets"] = pd.get_dummies(m["guide"].var["intended_target_name"]).astype(float)[names]
    perturbo.PERTURBO.setup_mudata(m, batch_key="batch", library_size_key="total_gene_umis", continuous_covariates_keys=["total_guide_umis"],
                                   guide_by_element_key="intended_targets", gene_by_element_key="intended_targets",
                                   modalities={"rna_layer": "gene", "perturbation_layer": "guide"})
    model = perturbo.PERTURBO(m, likelihood="nb")
    model.train(20, lr=0.01, batch_size=128)
    eff = (model.get_element_effects().rename(columns={"element": "intended_target_name", "gene": "gene_id", "q_value": "p_value"})
           .assign(log2_fc=lambda x: x["loc"] / np.log(2)).merge(ptt))
    return eff[[c for c in ["gene_id", "guide_id", "intended_target_name", "log2_fc", "p_value", "pair_type"] if c in eff.columns]]


# ---------------------------------------------------------------------------
# Inference (inference_sceptre.R, perturbo_inference.py, glm.nb comparator)
# ---------------------------------------------------------------------------

def negbinom_inference(md: MuDataLite, moi: Optional[str] = None):
    np, pd = _np(), _pd()
    pairs = uns_df(md, "pairs_to_test")
    pairs["gene_id"] = pairs["gene_id"].astype(str)
    pairs["intended_target_name"] = pairs["intended_target_name"].astype(str)
    moi = moi_of(md, moi)
    ntc = nt_cells(md.guide)
    umis = total_umis(md.gene)
    rows = []
    for (t, g), _ in pairs.groupby(["intended_target_name", "gene_id"]):
        trt = element_indicator(md.guide, t)
        keep = ((trt | ntc) if moi == "low" else np.ones(md.gene.n_obs, bool)) & (umis > 0)
        y = gene_counts(md.gene, [g])[keep, 0]
        x = trt[keep].astype(float)
        if x.sum() == 0 or x.sum() == len(x) or y.sum() == 0:
            rows.append({"intended_target_name": t, "gene_id": g, "p_value": np.nan, "log2_fc": np.nan})
            continue
        fit = glm_nb(y, np.column_stack([np.ones(len(y)), x]), offset=np.log(umis[keep]))
        rows.append({"intended_target_name": t, "gene_id": g, "p_value": float(fit["p"][1]), "log2_fc": float(fit["beta"][1] / math.log(2))})
    return pairs.merge(pd.DataFrame(rows), how="left", on=["intended_target_name", "gene_id"])


def run_method(md: MuDataLite, method: str, md_path: Optional[Path] = None, work: Optional[Path] = None, side: str = "both", B: int = 499,
               seed: int = 4, formula: str = "default", cov_string: str = ""):
    if method == "sceptre":
        return sceptre_inference(md, side, formula, cov_string, B, seed)[0]
    if method == "negbinom":
        return negbinom_inference(md)
    if method == "perturbo":
        if md_path is None:
            md_path = (work or Path(".")) / "perturbo_input.h5mu"
            write_h5mu(md, md_path)
        return infer_perturbo(md_path, (work or Path(".")) / "perturbo_inference_mudata.h5mu")
    raise SystemExit(f"unknown method {method}; choose from {METHODS}")


def cmd_prepare(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    md = read_h5mu(Path(args.mudata))
    pairs = read_table(Path(args.pairs))
    if "guide_id" not in md.guide.var.columns:
        md.guide.var["guide_id"] = [str(v) for v in md.guide.var_names]
    sub = restrict_pairs(pairs, md)
    md.uns["pairs_to_test"] = sub
    out = Path(args.out) if args.out else d / "mudata_inference_input.h5mu"
    write_h5mu(md, out)
    counts = sub["pair_type"].astype(str).value_counts().to_dict() if "pair_type" in sub.columns else {}
    write_report(d, "CRISPR Jamboree 4 — prepare_user_guide_inference",
                 [f"{len(pairs)} pairs in the sheet -> {len(sub)} with the gene and guide in the MuData.", md_table(["pair_type", "pairs"], counts.items())],
                 {"n_input": len(pairs), "n_pairs": len(sub), "pair_type": counts, "out": out})
    return 0


def cmd_infer(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    md = read_h5mu(Path(args.mudata))
    rows, done = [], {}
    for m in args.methods:
        t0 = time.time()
        tr = run_method(md, m, Path(args.mudata), d, args.side, args.resamples, args.seed, args.formula, args.cov_string or "")
        if tr is None:
            rows.append([m, "skipped (not installed)", "", ""])
            continue
        write_tsv(tr, d / f"test_results.{m}.tsv")
        done[m] = str(d / f"test_results.{m}.tsv")
        rows.append([m, len(tr), int((tr["p_value"] < 0.05).sum()), f"{time.time() - t0:.1f}s"])
    write_report(d, "CRISPR Jamboree 4 — inference", [md_table(["method", "pairs", "p < 0.05", "time"], rows)],
                 {"results": done, "methods": {r[0]: r[1] for r in rows}})
    return 0


# ---------------------------------------------------------------------------
# mergedResults (export_output_multiple.py / export_output_single.py)
# ---------------------------------------------------------------------------

def read_results(path: Path):
    if str(path).endswith(".h5mu"):
        return uns_df(read_h5mu(path), "test_results")
    return read_table(path)


def merge_method_results(results: dict):
    """{method: test_results} -> one table with <method>_log2_fc / <method>_p_value, inner-joined on (guide_id, gene_id)."""
    keys = None
    merged = None
    for m, tr in results.items():
        t = tr.copy()
        t["gene_id"] = t["gene_id"].astype(str)
        k = ["guide_id", "gene_id"] if "guide_id" in t.columns else ["intended_target_name", "gene_id"]
        keys = keys or k
        t = t.rename(columns={"log2_fc": f"{m}_log2_fc", "p_value": f"{m}_p_value"})
        if merged is None:
            merged = t
        else:
            cols = [c for c in keys if c in t.columns] + [f"{m}_log2_fc", f"{m}_p_value"]
            merged = merged.merge(t[cols], on=[c for c in keys if c in t.columns])
    return merged


def per_guide_table(merged, md: MuDataLite):
    np, pd = _np(), _pd()
    gv = md.guide.var.reset_index(drop=True).copy()
    if "guide_id" not in gv.columns:
        gv["guide_id"] = [str(v) for v in md.guide.var_names]
    extra = [c for c in gv.columns if c not in merged.columns or c in ("intended_target_name", "guide_id")]
    if "guide_id" in merged.columns:
        pg = merged.merge(gv[extra].drop_duplicates(["intended_target_name", "guide_id"]), how="left", on=["intended_target_name", "guide_id"])
    else:
        pg = merged.copy()
        pg["guide_id"] = ""
    pg = pg.rename(columns={"guide_id": "guide_id(s)"})
    sym = md.gene.var["symbol"].astype(str).values if "symbol" in md.gene.var.columns else np.asarray(md.gene.var_names).astype(str)
    cache: dict = {}
    avg = []
    for t in pg["intended_target_name"].astype(str):
        if t not in cache:
            mask = sym == t
            cache[t] = float(dense(md.gene.X[:, np.where(mask)[0]]).mean()) if mask.any() else np.nan
        avg.append(cache[t])
    pg["cell_number"] = md.gene.n_obs
    pg["avg_gene_expression"] = avg
    for c in ("sceptre_log2_fc", "sceptre_p_value", "perturbo_log2_fc", "perturbo_p_value"):
        if c not in pg.columns:
            pg[c] = None
    other = [c for c in pg.columns if c.endswith(("_log2_fc", "_p_value")) and c not in PER_GUIDE_COLS]
    return pg.reindex(columns=PER_GUIDE_COLS + other)


def per_element_table(pg, compat: bool = False):
    keys = ["gene_id"] if compat else ["intended_target_name", "gene_id"]
    ids = pg.groupby(keys, as_index=False).agg({"guide_id(s)": lambda x: ",".join(x.dropna().astype(str))})
    first = pg.drop_duplicates(keys).drop(columns=["guide_id(s)"])
    return ids.merge(first, on=keys, how="left")[list(pg.columns)]


def cmd_merge(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    md = read_h5mu(Path(args.mudata))
    results = {}
    if args.sceptre_results:
        results["sceptre"] = read_results(Path(args.sceptre_results))
    if args.perturbo_results:
        results["perturbo"] = read_results(Path(args.perturbo_results))
    for spec in args.extra_results or []:
        name, _, p = spec.partition("=")
        results[name] = read_results(Path(p))
    if not results:
        raise SystemExit("give --sceptre-results and/or --perturbo-results (or --extra-results name=path)")
    merged = merge_method_results(results)
    md.uns["test_results"] = merged
    write_h5mu(md, d / "inference_mudata.h5mu")
    pg = per_guide_table(merged, md)
    pe = per_element_table(pg, args.upstream_compat)
    write_tsv(pg, d / "per_guide_output.tsv")
    write_tsv(pe, d / "per_element_output.tsv")
    write_report(d, "CRISPR Jamboree 4 — merged results",
                 [f"Methods {', '.join(results)}: {len(merged)} merged pair rows, {len(pe)} per-element rows"
                  f"{' (grouped by gene_id, as upstream)' if args.upstream_compat else ''}."],
                 {"methods": list(results), "n_merged": len(merged), "n_per_element": len(pe), "per_guide": d / "per_guide_output.tsv"})
    return 0


# ---------------------------------------------------------------------------
# evaluation_curve.py
# ---------------------------------------------------------------------------

def control_labels(tr):
    np = _np()
    if "pair_type" not in tr.columns:
        raise SystemExit("results need a pair_type column to assign true labels")
    pt = tr["pair_type"].astype(str).str.strip().str.lower()
    pos, neg = pt.isin(POS_TYPES), pt.isin(NEG_TYPES)
    return pos.astype(int).values, (pos | neg).values


def evaluate_results_table(tr) -> dict:
    np = _np()
    lab, keep = control_labels(tr)
    out = {}
    for c in [c for c in tr.columns if str(c).endswith("_p_value")] or (["p_value"] if "p_value" in tr.columns else []):
        p = __import__("pandas").to_numeric(tr[c], errors="coerce").values.astype(float)
        score = -np.log10(np.clip(p, 1e-300, 1.0))
        out[c] = binary_metrics(lab[keep], score[keep])
    return out


def plot_curves(evals: dict, path: Path):
    plt = _plt()
    if plt is None or not evals:
        return None
    items = [(k, e) for k, e in evals.items() if "_curves" in e]
    if not items:
        return None
    fig, axes = plt.subplots(len(items), 2, figsize=(10, 4 * len(items)), squeeze=False)
    for i, (k, e) in enumerate(items):
        c = e["_curves"]
        axes[i, 0].plot(c["recall"], c["precision"], lw=1.2, color=BLUE, label=f"AUPRC={e['auprc']:.3f}")
        axes[i, 1].plot(c["fpr"], c["tpr"], lw=1.2, color=BLUE, label=f"AUROC={e['auroc']:.3f}")
        _style(axes[i, 0], f"{k} - Precision-Recall", "Recall", "Precision")
        _style(axes[i, 1], f"{k} - ROC", "False Positive Rate", "True Positive Rate")
        axes[i, 0].legend(frameon=False, fontsize=8)
        axes[i, 1].legend(frameon=False, fontsize=8)
    return _save(fig, path)


def cmd_evaluate(args: argparse.Namespace) -> int:
    pd = _pd()
    d = run_dir(args.label)
    sources = []
    for spec in args.results:
        name, sep, p = spec.partition("=")
        if not sep:
            name, p = Path(spec).parent.name or "results", spec
        sources.append((name, Path(p)))
    evals, rows = {}, []
    for name, p in sources:
        print(f"Reading data from {p}")
        tr = read_results(p)
        if "pair_type" not in tr.columns and args.pairs:
            lab = read_table(Path(args.pairs)).rename(columns={"gene_name": "gene_id"})
            gcol = "guide_id(s)" if "guide_id(s)" in tr.columns else "guide_id"
            lab = lab.rename(columns={"guide_id": gcol})[[gcol, "gene_id", "pair_type"]].drop_duplicates([gcol, "gene_id"])
            tr = tr.assign(gene_id=tr["gene_id"].astype(str)).merge(lab.assign(gene_id=lab["gene_id"].astype(str)), on=[gcol, "gene_id"], how="left")
        print("Assigning true labels...")
        print("Performing binary evaluation...")
        for col, e in evaluate_results_table(tr).items():
            key = f"{name}: {col}" if len(sources) > 1 else col
            evals[key] = e
            print(f"\nResults for {col}:\nArea under Precision-Recall Curve : {e['auprc']:.3f}\nArea under Receiver-Operating Curve: {e['auroc']:.3f}")
            rows.append({"dataset": name, "method": col.replace("_p_value", ""), "n_pairs": e["n"], "n_positive": e["n_pos"],
                         "n_negative": e["n_neg"], "auprc": e["auprc"], "auroc": e["auroc"], "average_precision": e["average_precision"]})
    tab = pd.DataFrame(rows)
    write_tsv(tab, d / "evaluation_metrics.tsv")
    figs = []
    if not args.no_plots:
        f = plot_curves(evals, d / "evaluation_curves.png")
        if f:
            figs.append(f)
    write_report(d, "CRISPR Jamboree 4 — control-set evaluation",
                 [md_table(["dataset", "method", "pairs", "positives", "negatives", "AUPRC", "AUROC", "AP"],
                           [[r["dataset"], r["method"], r["n_pairs"], r["n_positive"], r["n_negative"], _fmt(r["auprc"]), _fmt(r["auroc"]),
                             _fmt(r["average_precision"])] for r in rows])],
                 {"metrics": rows, "figures": figs})
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    md = read_h5mu(Path(args.mudata))
    if "guide_id" not in md.guide.var.columns:
        md.guide.var["guide_id"] = [str(v) for v in md.guide.var_names]
    md.uns["pairs_to_test"] = restrict_pairs(read_table(Path(args.pairs)), md)
    p_in = d / "mudata_inference_input.h5mu"
    write_h5mu(md, p_in)
    results = {}
    for m in args.methods:
        tr = run_method(md, m, p_in, d, args.side, args.resamples, args.seed, args.formula, args.cov_string or "")
        if tr is not None:
            results[m] = tr
            write_tsv(tr, d / f"test_results.{m}.tsv")
    if not results:
        raise SystemExit("no method produced results")
    merged = merge_method_results(results)
    md.uns["test_results"] = merged
    write_h5mu(md, d / "inference_mudata.h5mu")
    pg = per_guide_table(merged, md)
    write_tsv(pg, d / "per_guide_output.tsv")
    write_tsv(per_element_table(pg, args.upstream_compat), d / "per_element_output.tsv")
    evals = evaluate_results_table(merged)
    figs = []
    if not args.no_plots:
        f = plot_curves(evals, d / "evaluation_curves.png")
        if f:
            figs.append(f)
    rows = [[k.replace("_p_value", ""), e["n"], e["n_pos"], e["n_neg"], _fmt(e["auprc"]), _fmt(e["auroc"])] for k, e in evals.items()]
    write_report(d, f"CRISPR Jamboree 4 — {args.label}",
                 [f"{len(md.uns['pairs_to_test'])} pairs; methods {', '.join(results)}.",
                  md_table(["method", "pairs", "positives", "negatives", "AUPRC", "AUROC"], rows)],
                 {"methods": {k.replace("_p_value", ""): {kk: v for kk, v in e.items() if kk != "_curves"} for k, e in evals.items()},
                  "n_pairs": len(md.uns["pairs_to_test"]), "figures": figs})
    return 0


# ---------------------------------------------------------------------------
# Synthetic checkpoint + self-test
# ---------------------------------------------------------------------------

def synthetic_checkpoint(d: Path, seed: int = 21, n_cells: int = 1500) -> dict:
    """Guide-assigned MuData (inference checkpoint) + a pairs sheet with direct-targeting and negative controls."""
    np, pd = _np(), _pd()
    import anndata as ad  # type: ignore
    import scipy.sparse as sp  # type: ignore
    rng = np.random.default_rng(seed)
    n_genes = 200
    targets = [f"TF{i}" for i in range(1, 11)]
    symbols = targets + [f"G{i}" for i in range(n_genes - len(targets))]
    genes = [f"ENSG{800000 + i:08d}" for i in range(n_genes)]
    guide_ids = [f"{t}_sg{k}" for t in targets for k in (1, 2)] + [f"non-targeting_{k:05d}_sg1" for k in range(6)]
    tname = [g.rsplit("_sg", 1)[0] for g in guide_ids]
    ng = len(guide_ids)
    A = np.zeros((n_cells, ng), bool)
    for c in range(n_cells):
        A[c, rng.choice(ng, size=1 + rng.poisson(1.2), replace=False)] = True
    kd = {t: f for t, f in zip(targets, [0.2, 0.25, 0.3, 0.3, 0.35, 0.4, 0.45, 0.5, 1.0, 1.0])}
    lib = rng.lognormal(0, 0.3, n_cells)
    base = rng.gamma(2, 1.5, n_genes) + 0.5
    mu = lib[:, None] * base[None, :]
    for t, f in kd.items():
        cells = A[:, [i for i, x in enumerate(tname) if x == t]].any(axis=1)
        mu[cells, targets.index(t)] *= f
    X = rng.negative_binomial(5, 5 / (5 + mu)).astype(np.float64)
    gumi = np.where(A, rng.negative_binomial(3, 3 / 33.0, A.shape) + 2, rng.poisson(0.2, A.shape)).astype(np.float64)
    cells = [f"BC{c:06d}" for c in range(n_cells)]
    batch = np.where(np.arange(n_cells) % 2 == 0, "b1", "b2")
    gene = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame({"total_gene_umis": X.sum(axis=1), "num_expressed_genes": (X > 0).sum(axis=1),
                                                          "batch": batch}, index=cells),
                      var=pd.DataFrame({"symbol": symbols}, index=genes))
    gv = pd.DataFrame({"guide_id": guide_ids, "intended_target_name": tname,
                       "targeting": ["FALSE" if t.startswith("non-targeting") else "TRUE" for t in tname],
                       "intended_target_chr": ["chr2" if not t.startswith("non") else "non_targeting" for t in tname],
                       "intended_target_start": [2_000_000 + 100_000 * targets.index(t) if t in targets else 0 for t in tname],
                       "intended_target_end": [2_000_020 + 100_000 * targets.index(t) if t in targets else 0 for t in tname]}, index=guide_ids)
    guide = ad.AnnData(X=sp.csr_matrix(gumi), obs=pd.DataFrame({"total_guide_umis": gumi.sum(axis=1), "batch": batch}, index=cells), var=gv)
    guide.layers["guide_assignment"] = sp.csr_matrix(A.astype(np.float64))
    guide.uns["moi"] = np.array(["high"], dtype=object)
    guide.uns["capture_method"] = np.array(["CROP-seq"], dtype=object)
    md = MuDataLite({"gene": gene, "guide": guide}, pd.DataFrame({"batch": batch}, index=cells), {})
    p = d / "mudata_guide_assignment.h5mu"
    write_h5mu(md, p)
    rows = []
    for g, t in zip(guide_ids, tname):
        if t in targets:
            rows.append({"guide_id": g, "gene_name": genes[targets.index(t)], "intended_target_name": t, "pair_type": "Direct targeting"})
        else:
            for k in range(2):
                rows.append({"guide_id": g, "gene_name": genes[(guide_ids.index(g) + k) % len(targets)], "intended_target_name": t,
                             "pair_type": "Non-targeting_negative_control"})
    rows.append({"guide_id": "TF1_sg1", "gene_name": genes[150], "intended_target_name": "TF1", "pair_type": "discovery"})
    rows.append({"guide_id": "absent_sg9", "gene_name": genes[0], "intended_target_name": "TF1", "pair_type": "Direct targeting"})
    pairs = pd.DataFrame(rows)
    pp = d / "pairs_to_test.csv"
    pairs.to_csv(pp, index=False)
    return {"h5mu": p, "pairs": pp, "genes": genes, "targets": targets, "kd": kd, "n_pairs": len(rows) - 1, "n_cells": n_cells}


def cmd_selftest(args: argparse.Namespace) -> int:
    import tempfile
    np, pd = _np(), _pd()
    checks = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    def last(label):
        return sorted(OUT_ROOT.glob(f"*_{label}"))[-1]

    try:
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            W = synthetic_checkpoint(d)

            print("\nmetrics")
            m = binary_metrics([1, 1, 1, 0, 0], [5, 4, 1, 2, 0.5])
            check(abs(m["auroc"] - 5 / 6) < 1e-12 and 0 < m["auprc"] <= 1, f"binary_metrics: AUROC {m['auroc']:.3f} = 5/6; AUPRC by trapezoid {m['auprc']:.3f}")
            lab, keep = control_labels(pd.DataFrame({"pair_type": ["Direct targeting", "Positive_control", "Non-targeting_negative_control",
                                                                   "Targeting_negative_control", "discovery"]}))
            check(list(lab[keep]) == [1, 1, 0, 0] and not keep[-1], "labels: direct-targeting / positive = 1, non-targeting / targeting negatives = 0, discovery dropped")

            print("\nprepare")
            rc = cmd_prepare(argparse.Namespace(mudata=str(W["h5mu"]), pairs=str(W["pairs"]), out=str(d / "inference_input.h5mu"), label="st_j4_prep"))
            last("st_j4_prep")
            mdi = read_h5mu(d / "inference_input.h5mu")
            ptt = uns_df(mdi, "pairs_to_test")
            check(rc == 0 and len(ptt) == W["n_pairs"] and "gene_id" in ptt.columns, f"prepare: {len(ptt)} pairs kept (absent guide dropped); gene_name -> gene_id")

            print("\ninfer")
            sc_ = run_method(mdi, "sceptre", B=299)
            nb = run_method(mdi, "negbinom")
            for name, tr in (("sceptre", sc_), ("negbinom", nb)):
                dt = tr[tr["pair_type"] == "Direct targeting"].groupby("intended_target_name")["p_value"].max()
                check(all(dt[t] < 1e-3 for t in W["targets"][:8]) and dt["TF9"] > 1e-3 and dt["TF10"] > 1e-3 and len(tr) == len(ptt),
                      f"{name}: 8 knocked-down targets p <= {max(dt[t] for t in W['targets'][:8]):.1e}; unperturbed TF9 / TF10 p {dt['TF9']:.2f} / {dt['TF10']:.2f}")
            rc = cmd_infer(argparse.Namespace(mudata=str(d / "inference_input.h5mu"), methods=["sceptre", "negbinom", "perturbo"], side="both",
                                              resamples=199, seed=4, formula="default", cov_string=None, label="st_j4_infer"))
            r = last("st_j4_infer")
            s = json.loads((r / "summary.json").read_text())
            check(rc == 0 and set(s["results"]) == {"sceptre", "negbinom"} and "skipped" in str(s["methods"]["perturbo"]),
                  "infer: sceptre + negbinom tables; perturbo skipped when the package is absent")

            print("\nmerge")
            fake = nb.copy()
            rng = np.random.default_rng(3)
            fake["p_value"] = np.clip(fake["p_value"] * 10 ** rng.normal(0, 0.5, len(fake)), 1e-300, 1)
            fake = fake[fake["pair_type"] != "discovery"]
            fp = d / "perturbo_results.tsv"
            fake.to_csv(fp, sep="\t", index=False)
            rc = cmd_merge(argparse.Namespace(mudata=str(d / "inference_input.h5mu"), sceptre_results=str(r / "test_results.sceptre.tsv"),
                                              perturbo_results=str(fp), extra_results=None, upstream_compat=False, label="st_j4_merge"))
            rm = last("st_j4_merge")
            pg = pd.read_csv(rm / "per_guide_output.tsv", sep="\t")
            pe = pd.read_csv(rm / "per_element_output.tsv", sep="\t")
            tf1 = pg[(pg["intended_target_name"] == "TF1") & (pg["gene_id"] == W["genes"][0])]
            x_tf1 = dense(mdi.gene.X[:, 0]).mean()
            check(rc == 0 and list(pg.columns) == PER_GUIDE_COLS and len(pg) == len(fake) and (pg["cell_number"] == W["n_cells"]).all()
                  and abs(tf1["avg_gene_expression"].iloc[0] - x_tf1) < 1e-9 and tf1["intended_target_chr"].iloc[0] == "chr2",
                  "merge: export_output_multiple columns, inner join on (guide_id, gene_id), cell_number, avg_gene_expression of the target's own gene")
            e1 = pe[(pe["intended_target_name"] == "TF1") & (pe["gene_id"] == W["genes"][0])]
            check(set(e1["guide_id(s)"].iloc[0].split(",")) == {"TF1_sg1", "TF1_sg2"} and len(pe) == len(pg.drop_duplicates(["intended_target_name", "gene_id"])),
                  "merge: per-element table joins the element's guide ids, one row per (element, gene)")
            rc = cmd_merge(argparse.Namespace(mudata=str(d / "inference_input.h5mu"), sceptre_results=str(r / "test_results.sceptre.tsv"),
                                              perturbo_results=str(fp), extra_results=None, upstream_compat=True, label="st_j4_mergec"))
            pec = pd.read_csv(last("st_j4_mergec") / "per_element_output.tsv", sep="\t")
            check(rc == 0 and len(pec) == pg["gene_id"].nunique() < len(pe), f"merge --upstream-compat: grouped by gene_id -> {len(pec)} rows (elements sharing a gene collapse)")

            print("\nevaluate")
            rc = cmd_evaluate(argparse.Namespace(results=[f"H9={rm / 'inference_mudata.h5mu'}", f"per_guide={rm / 'per_guide_output.tsv'}"],
                                                 pairs=str(W["pairs"]), no_plots=args.no_plots, label="st_j4_eval"))
            re_ = last("st_j4_eval")
            met = pd.read_csv(re_ / "evaluation_metrics.tsv", sep="\t")
            sce = met[(met["dataset"] == "H9") & (met["method"] == "sceptre")].iloc[0]
            pgm = met[(met["dataset"] == "per_guide") & (met["method"] == "sceptre")].iloc[0]
            check(rc == 0 and set(met["method"]) == {"sceptre", "perturbo"} and sce["auroc"] > 0.9 and sce["auprc"] > 0.9
                  and sce["n_positive"] == 20 and sce["n_negative"] == 12 and abs(pgm["auroc"] - sce["auroc"]) < 1e-12,
                  f"evaluate: sceptre AUROC {sce['auroc']:.3f}, AUPRC {sce['auprc']:.3f} on 20 direct-targeting vs 12 non-targeting controls; inference_mudata and per_guide_output + pairs sheet agree")
            rand = pd.read_csv(rm / "per_guide_output.tsv", sep="\t")
            rand["pair_type"] = uns_df(read_h5mu(rm / "inference_mudata.h5mu"), "test_results")["pair_type"].values
            rand["sceptre_p_value"] = np.random.default_rng(0).uniform(size=len(rand))
            ev = evaluate_results_table(rand)
            check(ev["sceptre_p_value"]["auroc"] < 0.8, f"evaluate: random p-values -> AUROC {ev['sceptre_p_value']['auroc']:.2f}")

            print("\nrun")
            rc = cmd_run(argparse.Namespace(mudata=str(W["h5mu"]), pairs=str(W["pairs"]), methods=["sceptre", "negbinom"], side="both", resamples=199, seed=4,
                                            formula="default", cov_string=None, upstream_compat=False, no_plots=args.no_plots, label="st_j4_run"))
            rr = last("st_j4_run")
            s = json.loads((rr / "summary.json").read_text())
            check(rc == 0 and s["methods"]["sceptre"]["auroc"] > 0.9 and s["methods"]["negbinom"]["auroc"] > 0.9 and (rr / "per_element_output.tsv").is_file(),
                  f"run: prepare -> sceptre + negbinom -> merge -> evaluate (AUROC {s['methods']['sceptre']['auroc']:.3f} / {s['methods']['negbinom']['auroc']:.3f})")
            if not args.no_plots:
                check(len(s["figures"]) == 1, "run: evaluation_curves.png")
    finally:
        for p in sorted(OUT_ROOT.glob("*_st_j4_*")):
            shutil.rmtree(p, ignore_errors=True)
        if OUT_ROOT.is_dir() and not any(OUT_ROOT.iterdir()):
            OUT_ROOT.rmdir()
    ok = all(c for c, _ in checks)
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[list[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent crispr-jamboree4", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add_infer_opts(p):
        p.add_argument("--side", choices=["left", "right", "both"], default="both")
        p.add_argument("--resamples", type=int, default=499, help="sceptre resamples")
        p.add_argument("--seed", type=int, default=4)
        p.add_argument("--formula", default="default", help="'default' or any other value to add --cov-string covariates")
        p.add_argument("--cov-string", help="covariates, e.g. 'batch + lane'")

    p = sub.add_parser("prepare", help="user pairs -> uns['pairs_to_test']")
    p.add_argument("--mudata", required=True, help="guide-assigned MuData (the inference checkpoint)")
    p.add_argument("--pairs", required=True, help="CSV guide_id, gene_name, intended_target_name, pair_type")
    p.add_argument("--out")
    p.add_argument("--label", default="prepare")
    p.set_defaults(func=cmd_prepare)

    p = sub.add_parser("infer", help="sceptre / negbinom / perturbo inference")
    p.add_argument("--mudata", required=True)
    p.add_argument("--methods", nargs="+", choices=METHODS, default=["sceptre", "perturbo"])
    add_infer_opts(p)
    p.add_argument("--label", default="infer")
    p.set_defaults(func=cmd_infer)

    p = sub.add_parser("merge", help="method results -> inference_mudata + per_guide / per_element outputs")
    p.add_argument("--mudata", required=True)
    p.add_argument("--sceptre-results")
    p.add_argument("--perturbo-results")
    p.add_argument("--extra-results", nargs="+", help="name=path for other methods")
    p.add_argument("--upstream-compat", action="store_true", help="per-element table grouped by gene_id only")
    p.add_argument("--label", default="merge")
    p.set_defaults(func=cmd_merge)

    p = sub.add_parser("evaluate", help="AUPRC / AUROC of each *_p_value column against the control pairs")
    p.add_argument("--results", nargs="+", required=True, help="[name=]path to inference_mudata.h5mu / test_results / per_guide_output.tsv")
    p.add_argument("--pairs", help="pairs sheet with pair_type, joined when a results table lacks it (per_guide_output.tsv)")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="evaluate")
    p.set_defaults(func=cmd_evaluate)

    p = sub.add_parser("run", help="prepare -> infer -> merge -> evaluate")
    p.add_argument("--mudata", required=True)
    p.add_argument("--pairs", required=True)
    p.add_argument("--methods", nargs="+", choices=METHODS, default=["sceptre", "perturbo"])
    add_infer_opts(p)
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="jamboree4_run")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("selftest", help="synthetic checkpoint MuData + pairs sheet; all subcommands asserted")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
