#!/usr/bin/env python3
"""Fishash Table 2 reproduction: CLEANSER + SCEPTRE-mixture guide assignment scored on human/mouse barnyards (port of IGVF-CRISPR/fishash-table2-reproduction).

Port of https://github.com/IGVF-CRISPR/fishash-table2-reproduction (No LICENSE
file; Python + R + bash), which reproduces the SCEPTRE and CLEANSER rows of the
four-sample species-assignment benchmark in Table 2 of the Fishash preprint
(doi:10.64898/2026.01.22.701179) on the Liu et al. GSE272457 human (293T) /
mouse (NIH3T3) barnyard mixtures.  Every script was read (build_inputs.py,
run_sceptre_table2.R, score_table2.R, run_reproduction.sh) and the methods were
re-derived in Python; the upstream has no licence, so the DEFINITIONS were
ported, not the code.  The two guide callers the benchmark scores were also
re-implemented from their definitions:
  * CLEANSER 1.2.1 (Gersbachlab-Bioinformatics/CLEANSER @ d28c31b5, MIT): the
    zero-truncated Poisson/NB (CROP-seq, `cs`) and NB/NB (direct capture,
    `dc`) mixtures of cs-guide-mixture.stan / dc-guide-mixture.stan, with the
    same priors, bounds, library-size normalisation and "median posterior
    PZi" output.  Stan/NUTS is replaced by an adaptive random-walk Metropolis
    sampler on Stan's unconstrained scale (Jacobians included), started at
    the MAP; `--method map` gives the plug-in posterior at the MAP.
  * SCEPTRE 0.10.3 assign_grnas(method = "mixture"): per-gRNA Poisson GLM on
    the cell covariates, then the reduced two-component EM (pi and the
    perturbation effect g) with the GLM fit as offset, 5 random restarts,
    posterior >= 0.8; gRNAs with < 10 nonzero cells fall back to count >= 5.
    This is a Python approximation; `sceptre-mixture --engine r` writes and
    (when Rscript + sceptre 0.10.3 are installed) runs the upstream-equivalent
    R call instead.
Relationship: port.

Definitions, exactly as upstream computes them
  inputs (build_inputs.py)  features typed "CRISPR Guide Capture" are guides,
                   "Gene Expression" are genes; symbols prefixed GRCh38 are
                   human, mm10 mouse; homo_sum_gex / mus_sum_gex = UMI sums;
                   mito_sum = UMIs of GRCh38_MT-* + mm10___mt-* genes;
                   features_gex = number of genes with count > 0; guide
                   nt_<k>: k <= 100 -> homo_guide, else mus_guide.
  cohort (score_table2.R)  sum_gex = homo + mus; frac_human = homo / sum_gex;
                   selected = mito_sum/sum_gex < 0.15 & 1500 <= features_gex
                   <= 6000 & 3500 <= sum_gex <= 20000 & (frac_human < 0.1 |
                   frac_human > 0.9); truth (human) = frac_human >= 0.1 &
                   sum_gex > 100.
  prediction       has_human = any human guide assigned; has_mouse likewise;
                   prediction = has_human, NA when both or neither.
  accuracy         (tp + tn) / total over the selected cells, NA counted as
                   wrong; stderr = sqrt((T acc_T (1-acc_T) + F acc_F (1-acc_F))
                   / total^2) with T = tp + fn + NA_true, F = tn + fp +
                   NA_false, acc_T = tp/T, acc_F = tn/F; both_or_neither =
                   number of NA predictions.
  methods          cleanser cs (Table 2, 0.80), cleanser dc (Table 2, 0.50),
                   sceptre mixture, and the public-code 0.95 cut-offs for cs
                   and dc; a CLEANSER call is posterior >= cut-off.
  comparison       published Table 2 accuracy / SE per method x sample;
                   accuracy_delta, stderr_delta, *_match_4dp (round(x, 4) ==
                   published), accuracy_within_0_001.
  CLEANSER         L_i = (sum of cell i's guide counts <= lpf [2]; 0 -> 1) /
                   mean over cells; only nonzero guide counts are modelled
                   (zero-truncated likelihood).
                   cs: (1-r) Pois(x | lambda) + r NB2(x | mu L_i, phi), divided
                   by 1 - [(1-r) e^-lambda + r (phi/(phi+mu L_i))^phi];
                   lambda in (1e-6, 5), r in (1e-6, 0.1), mu > lambda,
                   phi > 1e-6; r ~ Beta(1,10), lambda ~ LN(log 1.1, log 1.1),
                   mu ~ LN(log 100, log 100), phi ~ Beta(1,10).
                   dc: (1-r) NB2(x | m0 L_i, d0) + r NB2(x | mu L_i, phi),
                   truncated likewise; m0 > 1e-6, d0 > 0.1, r in (1e-6, 0.1),
                   mu > 5, phi > 0.1; r ~ Beta(1,10), phi ~ LN(log 3, log 2),
                   m0 ~ LN(log 1.075, log 1.1), d0 ~ LN(log 1.1, log 1.1),
                   mu ~ LN(log 100, log 3).
                   PZi = exp(LP1 - logsumexp(LP0, LP1)); output = median over
                   draws; defaults chains 4, warmup 300, samples 1000,
                   seed + guide index per guide.
  SCEPTRE mixture  covariates log(response_n_nonzero), log(response_n_umis),
                   log(grna_n_nonzero + 1), log(grna_n_umis + 1) (+ batch);
                   probability_threshold 0.8, n_em_rep 5, pi_guess_range
                   (1e-5, 0.1), g_pert_guess_range log(10)..log(5000),
                   n_nonzero_cells_cutoff 10, backup_threshold 5.

Subcommands
  build-inputs     raw GSE272457-style 10x triplets -> per-sample guide count
                   .mtx (CLEANSER input), .h5ad, _meta.csv and _guides.csv.
  cleanser         CLEANSER cs / dc posterior per guide x cell (Python MCMC or
                   MAP; --engine cleanser runs the real binary if on PATH).
  sceptre-mixture  SCEPTRE-style mixture assignment (Python; --engine r
                   writes/runs the sceptre 0.10.3 R script).
  score            Table 2 cohort + species accuracy for every method x sample,
                   compared with the published Table 2 -> the upstream
                   extended_table2_reproduction.csv schema.
  compare          a scored CSV vs the published Table 2 and vs the upstream
                   repository's own reproduced accuracies.
  run              build-inputs -> cleanser cs + dc -> sceptre-mixture -> score.
  selftest         synthetic barnyard with planted species guides and ambient
                   noise; runs every subcommand and asserts the definitions.

Output: Docs/FishashTable2/<timestamp>_<label>/.  numpy, pandas, scipy
required; anndata optional (.h5ad skipped); matplotlib optional.

Usage:
    igvfagent fishash-table2 build-inputs --raw-dir GSE272457 --work-dir work
    igvfagent fishash-table2 cleanser --input work/mix0hr_Cropseq_grna_counts.mtx --mode cs --out work/mix0hr_Cropseq_cleanser_cs_posterior.mtx
    igvfagent fishash-table2 sceptre-mixture --work-dir work --sample mix0hr_Cropseq --raw-prefix GSE272457/GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix
    igvfagent fishash-table2 score --work-dir work --label table2
    igvfagent fishash-table2 compare --scored Docs/FishashTable2/<run>/extended_table2_reproduction.csv
    igvfagent fishash-table2 run --raw-dir GSE272457 --work-dir work --label table2
    igvfagent fishash-table2 selftest --no-plots
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "FishashTable2"

UPSTREAM_REPO = "IGVF-CRISPR/fishash-table2-reproduction"
UPSTREAM_COMMIT = "7d5d9763368682214952719931e86229a075a0be"
CLEANSER_COMMIT = "d28c31b52cdaece66605896852abf502eb595d3a"
MAX_SEED_INT = 4_294_967_295

SAMPLES = {
    "mix0hr_Cropseq": "GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix",
    "mix72hr_Cropseq": "GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_72hr_mix",
    "mix0hr_DirectCapture": "GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_0hr_mix",
    "mix72hr_DirectCapture": "GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_72hr_mix",
}
SAMPLE_ORDER = list(SAMPLES)

# method -> (kind, mode, cutoff)
SPECS = [
    ("cleanser cs (Table 2, 0.80)", "cleanser", "cs", 0.80),
    ("cleanser dc (Table 2, 0.50)", "cleanser", "dc", 0.50),
    ("sceptre mixture", "mtx", None, None),
    ("cleanser cs (public-code 0.95)", "cleanser", "cs", 0.95),
    ("cleanser dc (public-code 0.95)", "cleanser", "dc", 0.95),
]
PUBLISHED = {
    "cleanser cs (Table 2, 0.80)": [0.9414, 0.9410, 0.8612, 0.9221],
    "cleanser dc (Table 2, 0.50)": [0.9178, 0.8459, 0.8614, 0.9308],
    "sceptre mixture": [0.9491, 0.9472, 0.8953, 0.9408],
}
PUBLISHED_SE = {
    "cleanser cs (Table 2, 0.80)": [0.0030, 0.0027, 0.0039, 0.0028],
    "cleanser dc (Table 2, 0.50)": [0.0034, 0.0041, 0.0039, 0.0027],
    "sceptre mixture": [0.0028, 0.0026, 0.0036, 0.0025],
}
# accuracies the upstream repository itself reproduced (results/extended_table2_reproduction.csv @ UPSTREAM_COMMIT)
UPSTREAM_REPRODUCED = {
    "cleanser cs (Table 2, 0.80)": [0.941223703227879, 0.940993368520774, 0.86123101518785, 0.922124380847551],
    "cleanser dc (Table 2, 0.50)": [0.917777420908945, 0.845851942076059, 0.860911270983213, 0.930792515134838],
    "sceptre mixture": [0.949092660992452, 0.947218838814454, 0.895283772981615, 0.94083654375344],
    "cleanser cs (public-code 0.95)": [0.933194154488518, 0.933143862498308, 0.855635491606715, 0.918409466153],
    "cleanser dc (public-code 0.95)": [0.847759755901718, 0.694275274056029, 0.84220623501199, 0.92294991744634],
}
UPSTREAM_N_QC = [6227, 7389, 6255, 7268]
SCORED_COLS = ["sample", "method", "cutoff", "accuracy", "stderr", "n_qc", "both_or_neither",
               "nonfatal_sampler_warnings", "convergence_warnings", "divergent_chain_warnings",
               "published_accuracy", "published_stderr", "accuracy_delta", "stderr_delta",
               "accuracy_match_4dp", "stderr_match_4dp", "accuracy_within_0_001"]

INK, INK2, AXIS, SURFACE = "#1f2328", "#57606a", "#8c959f", "#ffffff"
PALETTE = ["#3b6fb6", "#d9822b", "#3a9d5d", "#9b59b6", "#c0392b"]


# ---------------------------------------------------------------------------
# plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"fishash_table2_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
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
        raise SystemExit("this subcommand needs pandas + numpy + scipy: pip install 'igvfagent[analysis]'") from e


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


def md_table(headers: "list[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 80) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for i, r in enumerate(rows):
        if i >= max_rows:
            lines.append("| ... |" + " |" * (len(headers) - 1))
            break
        lines.append("| " + " | ".join(str(c).replace("|", "\\|") for c in r) + " |")
    return "\n".join(lines)


def _fmt(x, nd: int = 4) -> str:
    try:
        if x is None or (isinstance(x, float) and not math.isfinite(x)):
            return "NA"
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return str(x)


def _jsonable(o):
    np = _np()
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return f if math.isfinite(f) else None
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
    return o


def write_json(obj, path: Path) -> Path:
    path.write_text(json.dumps(_jsonable(obj), indent=2))
    print(f"JSON: {path}")
    return path


def write_csv(df, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, na_rep="")
    print(f"CSV: {path}")
    return path


def write_report(path: Path, text: str) -> Path:
    path.write_text(text)
    print(f"Report: {path}")
    return path


# ---------------------------------------------------------------------------
# build-inputs (build_inputs.py)
# ---------------------------------------------------------------------------

def read_features(path: Path):
    pd = _pd()
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        rows = [line.rstrip("\n").split("\t")[:3] for line in fh if line.strip()]
    return pd.DataFrame(rows, columns=["ID", "Symbol", "type"])


def read_barcodes(path: Path) -> "list[str]":
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        return [line.strip() for line in fh if line.strip()]


def read_mm(path: Path):
    import scipy.io  # type: ignore
    if str(path).endswith(".gz"):
        with gzip.open(path, "rb") as fh:
            return scipy.io.mmread(fh).tocsr()
    return scipy.io.mmread(str(path)).tocsr()


def guide_index(name: str) -> int:
    return int(str(name).replace("nt_", ""))


def build_sample(prefix: str, sample: str, out_dir: Path, human_max_index: int = 100) -> dict:
    """One raw 10x triplet -> the four upstream per-sample inputs."""
    pd, np = _pd(), _np()
    import scipy.io  # type: ignore
    features = read_features(Path(f"{prefix}_features.tsv.gz"))
    barcodes = read_barcodes(Path(f"{prefix}_barcodes.tsv.gz"))
    counts = read_mm(Path(f"{prefix}_matrix.mtx.gz"))
    if counts.shape != (len(features), len(barcodes)):
        raise SystemExit(f"{prefix}: matrix {counts.shape} vs {len(features)} features x {len(barcodes)} barcodes")
    is_guide = (features["type"] == "CRISPR Guide Capture").to_numpy()
    is_gex = (features["type"] == "Gene Expression").to_numpy()
    if not np.all(is_guide | is_gex):
        raise SystemExit(f"{prefix}: feature types other than Gene Expression / CRISPR Guide Capture present")
    sym = features["Symbol"].astype(str)
    refs = np.where(sym.str.startswith("GRCh38"), "homo", np.where(sym.str.startswith("mm10"), "mus", "unknown"))
    if np.any(refs[is_gex] == "unknown"):
        raise SystemExit(f"{prefix}: gene symbols must be prefixed GRCh38 / mm10 (barnyard reference)")
    gex, grna = counts[is_gex], counts[is_guide]
    gex_refs = refs[is_gex]
    symbols = sym[is_gex].to_numpy()
    homo_sum = np.asarray(gex[gex_refs == "homo"].sum(axis=0)).ravel()
    mus_sum = np.asarray(gex[gex_refs == "mus"].sum(axis=0)).ravel()
    is_mito = np.array([s.startswith("GRCh38_MT-") or s.startswith("mm10___mt-") for s in symbols])
    mito_sum = np.asarray(gex[is_mito].sum(axis=0)).ravel() if is_mito.any() else np.zeros(len(barcodes))
    features_gex = np.asarray((gex > 0).sum(axis=0)).ravel()
    guide_names = sym[is_guide].to_numpy()
    gidx = np.array([guide_index(x) for x in guide_names])
    gtype = np.where(gidx <= human_max_index, "homo_guide", "mus_guide")
    out_dir.mkdir(parents=True, exist_ok=True)
    mtx_path = out_dir / f"{sample}_grna_counts.mtx"
    scipy.io.mmwrite(str(mtx_path), grna.tocoo().astype(np.int64))
    print(f"Wrote: {mtx_path}")
    try:
        import anndata as ad  # type: ignore
        a = ad.AnnData(grna.T.tocsr().astype(np.int64))
        a.obs_names = barcodes
        a.var_names = list(guide_names)
        a.obs["batch"] = 0
        h5 = out_dir / f"{sample}_grna.h5ad"
        a.write_h5ad(str(h5))
        print(f"Wrote: {h5}")
    except ImportError:
        logging.warning("anndata missing: %s_grna.h5ad not written", sample)
    write_csv(pd.DataFrame({"barcode": barcodes, "homo_sum_gex": homo_sum, "mus_sum_gex": mus_sum,
                            "mito_sum": mito_sum, "features_gex": features_gex}), out_dir / f"{sample}_meta.csv")
    write_csv(pd.DataFrame({"guide": guide_names, "gidx": gidx, "guide_type": gtype}), out_dir / f"{sample}_guides.csv")
    info = {"sample": sample, "cells": len(barcodes), "guides": int(is_guide.sum()), "gex": int(is_gex.sum())}
    print(f"[{sample}] cells={info['cells']} guides={info['guides']} gex={info['gex']}")
    return info


def sample_prefixes(args) -> "list[tuple[str, str]]":
    if getattr(args, "raw_prefix", None):
        return [(args.sample or Path(args.raw_prefix).name, args.raw_prefix)]
    samples = args.samples or SAMPLE_ORDER
    return [(s, str(Path(args.raw_dir) / SAMPLES.get(s, s))) for s in samples]


def cmd_build_inputs(args: argparse.Namespace) -> int:
    out = Path(args.work_dir)
    infos = [build_sample(prefix, s, out, args.human_max_index) for s, prefix in sample_prefixes(args)]
    write_json({"samples": infos}, out / "build_inputs.json")
    return 0


# ---------------------------------------------------------------------------
# CLEANSER (cs / dc mixture models)
# ---------------------------------------------------------------------------

def read_cleanser_mtx(path: Path) -> "tuple[str, list[tuple[int, int, int]]]":
    """CLEANSER's reader: skip '%' lines, keep the dims line, then guide cell count triples."""
    header, rows = "", []
    with open(path) as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            header = line
            break
        for line in fh:
            if not line.strip():
                continue
            g, c, n = line.split()[:3]
            rows.append((int(g), int(c), int(float(n))))
    return header, rows


def cleanser_normalise(rows: "list[tuple[int, int, int]]", lpf: int = 2) -> "dict[int, float]":
    """L_i: per-cell sum of guide counts <= lpf (all counts when lpf == 0), zeros -> 1, divided by the mean."""
    tot: "dict[int, int]" = {}
    for _, c, n in rows:
        if c not in tot:
            tot[c] = 0
        if lpf:
            if n <= lpf:
                tot[c] += n
        else:
            tot[c] += n
    for k, v in tot.items():
        if v == 0:
            tot[k] = 1
    avg = sum(tot.values()) / max(len(tot), 1)
    return {k: v / avg for k, v in tot.items()}


def _nb2_lpmf(x, mu, phi):
    from scipy.special import gammaln  # type: ignore
    np = _np()
    return (gammaln(x + phi) - gammaln(phi) - gammaln(x + 1.0) + phi * (np.log(phi) - np.log(phi + mu))
            + x * (np.log(mu) - np.log(phi + mu)))


def _pois_lpmf(x, lam):
    from scipy.special import gammaln  # type: ignore
    np = _np()
    return x * np.log(lam) - lam - gammaln(x + 1.0)


def _lognormal_lpdf(y, m, s):
    return -math.log(y) - math.log(s) - 0.5 * math.log(2 * math.pi) - (math.log(y) - m) ** 2 / (2 * s * s)


def _beta_1_10_lpdf(y):
    if not (0.0 < y < 1.0):
        return -math.inf
    return math.log(10.0) + 9.0 * math.log1p(-y)


def _log1mexp(a):
    """log(1 - exp(a)) for a < 0."""
    np = _np()
    a = np.minimum(a, -1e-300)
    return np.where(a > -0.693, np.log(-np.expm1(a)), np.log1p(-np.exp(a)))


def _inv_logit(u):
    return 1.0 / (1.0 + math.exp(-u)) if u >= 0 else math.exp(u) / (1.0 + math.exp(u))


def _bounded(u, lo, hi):
    """Stan lower-upper transform and its log-Jacobian."""
    p = _inv_logit(u)
    val = lo + (hi - lo) * p
    lj = math.log(hi - lo) + (math.log(p) if p > 0 else -745.0) + (math.log1p(-p) if p < 1 else -745.0)
    return val, lj


class CleanserModel:
    """One guide's truncated mixture on Stan's unconstrained scale."""

    def __init__(self, mode: str, X, L):
        np = _np()
        self.mode = mode
        self.X = np.asarray(X, dtype=float)
        self.L = np.asarray(L, dtype=float)
        self.dim = 4 if mode == "cs" else 5

    def constrain(self, u) -> "tuple[dict, float]":
        if self.mode == "cs":
            lam, j1 = _bounded(u[0], 1e-6, 5.0)
            r, j2 = _bounded(u[1], 1e-6, 0.1)
            mu = lam + math.exp(min(u[2], 700)); j3 = u[2]
            phi = 1e-6 + math.exp(min(u[3], 700)); j4 = u[3]
            return {"lambda": lam, "r": r, "nbMean": mu, "nbDisp": phi}, j1 + j2 + j3 + j4
        m0 = 1e-6 + math.exp(min(u[0], 700)); j1 = u[0]
        d0 = 0.1 + math.exp(min(u[1], 700)); j2 = u[1]
        r, j3 = _bounded(u[2], 1e-6, 0.1)
        mu = 5.0 + math.exp(min(u[3], 700)); j4 = u[3]
        phi = 0.1 + math.exp(min(u[4], 700)); j5 = u[4]
        return {"n_nbMean": m0, "n_nbDisp": d0, "r": r, "nbMean": mu, "nbDisp": phi}, j1 + j2 + j3 + j4 + j5

    def components(self, p: dict):
        """LP0, LP1 per cell (Stan's generated quantities) and the truncation log-denominator."""
        np = _np()
        X, L, r = self.X, self.L, p["r"]
        mu = p["nbMean"] * L
        phi = p["nbDisp"]
        lp1 = math.log(r) + _nb2_lpmf(X, mu, phi)
        lz1 = math.log(r) + phi * (math.log(phi) - np.log(phi + mu))
        if self.mode == "cs":
            lam = p["lambda"]
            lp0 = math.log1p(-r) + _pois_lpmf(X, lam)
            lz0 = np.full_like(X, math.log1p(-r) - lam)
        else:
            m0 = p["n_nbMean"] * L
            d0 = p["n_nbDisp"]
            lp0 = math.log1p(-r) + _nb2_lpmf(X, m0, d0)
            lz0 = math.log1p(-r) + d0 * (math.log(d0) - np.log(d0 + m0))
        lz = np.logaddexp(lz0, lz1)
        return lp0, lp1, _log1mexp(lz)

    def log_prior(self, p: dict) -> float:
        if self.mode == "cs":
            return (_beta_1_10_lpdf(p["r"]) + _lognormal_lpdf(p["lambda"], math.log(1.1), math.log(1.1))
                    + _lognormal_lpdf(p["nbMean"], math.log(100), math.log(100)) + _beta_1_10_lpdf(p["nbDisp"]))
        return (_beta_1_10_lpdf(p["r"]) + _lognormal_lpdf(p["nbDisp"], math.log(3), math.log(2))
                + _lognormal_lpdf(p["n_nbMean"], math.log(1.075), math.log(1.1))
                + _lognormal_lpdf(p["n_nbDisp"], math.log(1.1), math.log(1.1))
                + _lognormal_lpdf(p["nbMean"], math.log(100), math.log(3)))

    def log_post(self, u) -> float:
        np = _np()
        try:
            p, lj = self.constrain(u)
        except (OverflowError, ValueError):
            return -math.inf
        lpr = self.log_prior(p)
        if not math.isfinite(lpr):
            return -math.inf
        with np.errstate(all="ignore"):
            lp0, lp1, ldenom = self.components(p)
            ll = float(np.sum(np.logaddexp(lp0, lp1) - ldenom))
        if not math.isfinite(ll):
            return -math.inf
        return ll + lpr + lj

    def pzi(self, u):
        np = _np()
        p, _ = self.constrain(u)
        with np.errstate(all="ignore"):
            lp0, lp1, _ = self.components(p)
            return np.exp(lp1 - np.logaddexp(lp0, lp1))

    def init_points(self) -> "list[list[float]]":
        np = _np()
        X = self.X
        hi = X[X > 5]
        mu0 = float(np.median(hi / np.maximum(self.L[X > 5], 1e-3))) if hi.size else 50.0
        mu0 = min(max(mu0, 6.0), 1e5)
        lg = lambda v: math.log(max(v, 1e-9))  # noqa: E731
        logit = lambda q: math.log(q / (1 - q))  # noqa: E731
        pts = []
        for rr in (0.02, 0.05, 0.09):
            if self.mode == "cs":
                lam = 0.5
                pts.append([logit((lam - 1e-6) / (5 - 1e-6)), logit((rr - 1e-6) / (0.1 - 1e-6)), lg(mu0 - lam), lg(0.3)])
            else:
                pts.append([lg(1.0), lg(1.0), logit((rr - 1e-6) / (0.1 - 1e-6)), lg(mu0 - 5.0), lg(2.0)])
        return pts


def cleanser_map(model: CleanserModel):
    from scipy.optimize import minimize  # type: ignore
    np = _np()
    best = None
    for x0 in model.init_points():
        f = lambda u: -model.log_post(u) if math.isfinite(model.log_post(u)) else 1e300  # noqa: E731
        res = minimize(f, np.asarray(x0, dtype=float), method="Nelder-Mead",
                       options={"maxiter": 4000, "xatol": 1e-6, "fatol": 1e-8})
        if best is None or res.fun < best.fun:
            best = res
    return np.asarray(best.x, dtype=float)


def cleanser_mcmc(model: CleanserModel, u0, n_warmup: int, n_samples: int, rng) -> "tuple[Any, dict]":
    """Adaptive random-walk Metropolis (Haario) on the unconstrained scale; returns the post-warmup draws."""
    np = _np()
    d = model.dim
    u = np.asarray(u0, dtype=float).copy()
    lp = model.log_post(u)
    scale = 0.1
    cov = np.eye(d) * 0.05
    hist = []
    draws = np.empty((n_samples, d))
    acc_w = acc_s = 0
    for it in range(n_warmup + n_samples):
        if it < n_warmup and it >= 50 and it % 25 == 0 and len(hist) > 20:
            H = np.asarray(hist)
            cov = np.cov(H.T) + np.eye(d) * 1e-6
        try:
            step = rng.multivariate_normal(np.zeros(d), (scale ** 2) * cov * (2.38 ** 2 / d))
        except np.linalg.LinAlgError:
            step = rng.normal(0, 0.05, size=d)
        prop = u + step
        lpp = model.log_post(prop)
        if math.isfinite(lpp) and math.log(rng.uniform()) < lpp - lp:
            u, lp = prop, lpp
            if it < n_warmup:
                acc_w += 1
            else:
                acc_s += 1
        if it < n_warmup:
            hist.append(u.copy())
            if it >= 20 and it % 20 == 0:  # tune the global scale toward ~0.3 acceptance
                rate = acc_w / (it + 1)
                scale *= math.exp(rate - 0.3)
                scale = min(max(scale, 0.05), 5.0)
        else:
            draws[it - n_warmup] = u
    return draws, {"acceptance": acc_s / max(n_samples, 1), "warmup_acceptance": acc_w / max(n_warmup, 1)}


def run_cleanser(input_path: Path, mode: str, out_path: Path, method: str = "mcmc", chains: int = 4,
                 n_warmup: int = 300, n_samples: int = 1000, seed: int = 20260810, lpf: int = 2,
                 samples_out: "Optional[Path]" = None, log_path: "Optional[Path]" = None) -> dict:
    """CLEANSER in Python: per-guide posteriors written in CLEANSER's banner-less coordinate format."""
    np = _np()
    header, rows = read_cleanser_mtx(input_path)
    rows = sorted(rows, key=lambda t: (t[0], t[1]))
    L = cleanser_normalise(rows, lpf)
    per_guide: "dict[int, list[tuple[int, int]]]" = {}
    for g, c, n in rows:
        per_guide.setdefault(g, []).append((c, n))
    lines, stats, sample_rows = [], [], []
    logf = open(log_path, "w") if log_path else None
    try:
        for g in sorted(per_guide):
            cells = per_guide[g]
            X = [n for _, n in cells]
            Lg = [L[c] for c, _ in cells]
            model = CleanserModel(mode, X, Lg)
            u_map = cleanser_map(model)
            if method == "map":
                post = model.pzi(u_map)
                pmed = model.constrain(u_map)[0]
                diag = {"acceptance": None}
            else:
                pz_draws, par_draws, accs = [], [], []
                for ch in range(max(chains, 1)):
                    rng = np.random.default_rng((seed + g + 7919 * ch) % MAX_SEED_INT)
                    start = u_map + rng.normal(0, 0.05, size=model.dim)
                    draws, dg = cleanser_mcmc(model, start, n_warmup, n_samples, rng)
                    accs.append(dg["acceptance"])
                    thin = draws[:: max(1, n_samples // 250)]
                    pz_draws.extend(model.pzi(u) for u in thin)
                    par_draws.extend(model.constrain(u)[0] for u in draws)
                post = np.median(np.asarray(pz_draws), axis=0)
                pmed = {k: float(np.median([p[k] for p in par_draws])) for k in par_draws[0]}
                diag = {"acceptance": float(np.mean(accs))}
                if samples_out is not None:
                    for p in par_draws:
                        sample_rows.append((g, p))
                if logf and diag["acceptance"] is not None and diag["acceptance"] < 0.05:
                    logf.write(f"guide {g}: Some chains may have failed to converge (acceptance {diag['acceptance']:.3f})\n")
            for (c, _), pz in zip(cells, post):
                lines.append(f"{g}\t{c}\t{float(pz)}\n")
            st = {"guide": g, "n_cells": len(cells), **{k: float(v) for k, v in pmed.items()}, **diag}
            stats.append(st)
            if logf:
                logf.write("\t".join(f"{k}={v}" for k, v in st.items()) + "\n")
        if logf:
            logf.write(f"Random seed: {seed}\n")
    finally:
        if logf:
            logf.close()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(header if header.endswith("\n") else header + "\n")
        fh.writelines(lines)
    print(f"Wrote: {out_path}")
    if samples_out is not None and sample_rows:
        keys = list(sample_rows[0][1])
        with open(samples_out, "w") as fh:
            fh.write("guide id\t" + "\t".join(keys) + "\n")
            for g, p in sample_rows:
                fh.write(f"{g}\t" + "\t".join(str(p[k]) for k in keys) + "\n")
        print(f"Wrote: {samples_out}")
    return {"mode": mode, "method": method, "n_guides": len(per_guide), "n_entries": len(lines), "guides": stats}


def cleanser_binary_cmd(input_path: str, out_path: str, mode: str, cpus: int, seed: int) -> "list[str]":
    return ["cleanser", "-i", input_path, "-o", out_path, f"--{mode}", "-p", str(cpus), "-s", str(seed)]


def cmd_cleanser(args: argparse.Namespace) -> int:
    inp, out = Path(args.input), Path(args.out)
    log_path = Path(args.log) if args.log else out.with_name(out.name.replace("_posterior.mtx", "") + ".log")
    if args.engine == "cleanser":
        cmd = cleanser_binary_cmd(str(inp), str(out), args.mode, args.cpus, args.seed)
        if shutil.which("cleanser") is None:
            print("cleanser not on PATH; the upstream command would be:\n  " + " ".join(cmd) + f" > {log_path} 2>&1")
            print("Install: pip install 'CLEANSER @ git+https://github.com/Gersbachlab-Bioinformatics/CLEANSER.git@" + CLEANSER_COMMIT + "'")
            return 0
        with open(log_path, "w") as lf:
            rc = subprocess.call(cmd, stdout=lf, stderr=subprocess.STDOUT)
        print(f"Wrote: {out}")
        return rc
    res = run_cleanser(inp, args.mode, out, method=args.method, chains=args.chains, n_warmup=args.num_warmup,
                       n_samples=args.num_samples, seed=args.seed, lpf=args.lpf,
                       samples_out=Path(args.samples_out) if args.samples_out else None, log_path=log_path)
    print(f"Wrote: {log_path}")
    if args.threshold is not None:
        n = sum(1 for line in open(out).readlines()[1:] if float(line.split()[2]) >= args.threshold)
        print(f"calls at posterior >= {args.threshold}: {n}")
    write_json(res, out.with_suffix(".stats.json"))
    return 0


# ---------------------------------------------------------------------------
# SCEPTRE mixture assignment (Python approximation + R runner)
# ---------------------------------------------------------------------------

def poisson_glm(y, X, offset=None, max_iter: int = 50, tol: float = 1e-8):
    """IRLS Poisson regression with log link; returns (beta, linear predictor)."""
    np = _np()
    n, p = X.shape
    off = np.zeros(n) if offset is None else offset
    beta = np.zeros(p)
    ybar = max(float(np.mean(y)), 1e-8)
    beta[0] = math.log(ybar)
    eta = X @ beta + off
    for _ in range(max_iter):
        mu = np.exp(np.clip(eta, -30, 30))
        z = eta - off + (y - mu) / mu
        W = mu
        XtW = X.T * W
        try:
            new = np.linalg.solve(XtW @ X + 1e-8 * np.eye(p), XtW @ z)
        except np.linalg.LinAlgError:
            break
        if np.max(np.abs(new - beta)) < tol:
            beta = new
            break
        beta = new
        eta = X @ beta + off
    eta = X @ beta + off
    return beta, eta


def reduced_em(y, offsets, pi0: float, g0: float, max_iter: int = 200, tol: float = 1e-6):
    """Two-component Poisson mixture with fixed offsets: y ~ (1-pi) Pois(e^o) + pi Pois(e^(o+g))."""
    np = _np()
    from scipy.special import gammaln  # type: ignore
    base = np.exp(np.clip(offsets, -30, 30))
    lgy = gammaln(y + 1.0)
    pi, g = pi0, g0
    prev = -math.inf
    for _ in range(max_iter):
        l0 = math.log1p(-pi) + y * offsets - base - lgy
        l1 = math.log(pi) + y * (offsets + g) - base * math.exp(g) - lgy
        lz = np.logaddexp(l0, l1)
        T = np.exp(l1 - lz)
        ll = float(np.sum(lz))
        pi = min(max(float(np.mean(T)), 1e-10), 1 - 1e-10)
        num, den = float(np.sum(T * y)), float(np.sum(T * base))
        if num <= 0 or den <= 0:
            break
        g = math.log(num / den)
        if abs(ll - prev) < tol * (1 + abs(ll)):
            break
        prev = ll
    l0 = math.log1p(-pi) + y * offsets - base - lgy
    l1 = math.log(pi) + y * (offsets + g) - base * math.exp(g) - lgy
    lz = np.logaddexp(l0, l1)
    return {"pi": pi, "g": g, "loglik": float(np.sum(lz)), "posterior": np.exp(l1 - lz)}


def sceptre_covariates(gex_counts, grna_counts, batch=None, grna_covariates: bool = True):
    """Design matrix of sceptre's default high-MOI formula (cells x p), intercept first."""
    np = _np()
    resp_umis = np.asarray(gex_counts.sum(axis=0)).ravel().astype(float)
    resp_nz = np.asarray((gex_counts > 0).sum(axis=0)).ravel().astype(float)
    cols = [np.ones_like(resp_umis), np.log(np.maximum(resp_nz, 1)), np.log(np.maximum(resp_umis, 1))]
    if grna_covariates:
        g_umis = np.asarray(grna_counts.sum(axis=0)).ravel().astype(float)
        g_nz = np.asarray((grna_counts > 0).sum(axis=0)).ravel().astype(float)
        cols += [np.log(g_nz + 1), np.log(g_umis + 1)]
    if batch is not None and len(set(batch)) > 1:
        levels = sorted(set(batch))
        for lv in levels[1:]:
            cols.append((np.asarray(batch) == lv).astype(float))
    X = np.column_stack(cols)
    # standardise non-intercept columns for IRLS stability (same fitted values)
    mu, sd = X[:, 1:].mean(axis=0), X[:, 1:].std(axis=0)
    sd[sd == 0] = 1.0
    X[:, 1:] = (X[:, 1:] - mu) / sd
    return X


def sceptre_mixture_assign(grna_counts, X, probability_threshold: float = 0.8, n_em_rep: int = 5,
                           pi_guess_range=(1e-5, 0.1), g_pert_guess_range=(math.log(10), math.log(5000)),
                           n_nonzero_cells_cutoff: int = 10, backup_threshold: int = 5, seed: int = 4):
    """guides x cells boolean assignment + per-guide fit table."""
    np, pd = _np(), _pd()
    import scipy.sparse as sp  # type: ignore
    rng = np.random.default_rng(seed)
    G = sp.csr_matrix(grna_counts)
    n_g, n_c = G.shape
    rows, cols, fits = [], [], []
    for gi in range(n_g):
        y = np.asarray(G[gi].todense()).ravel().astype(float)
        nnz = int((y > 0).sum())
        if nnz < n_nonzero_cells_cutoff:
            idx = np.where(y >= backup_threshold)[0]
            fits.append({"grna_index": gi, "method": "backup_threshold", "n_nonzero": nnz, "pi": np.nan, "g_pert": np.nan, "n_assigned": len(idx)})
        else:
            _, eta = poisson_glm(y, X)
            best = None
            for _ in range(n_em_rep):
                pi0 = rng.uniform(*pi_guess_range)
                g0 = rng.uniform(*g_pert_guess_range)
                res = reduced_em(y, eta, pi0, g0)
                if best is None or res["loglik"] > best["loglik"]:
                    best = res
            if best["g"] <= 0:  # the "perturbed" component must be the high one
                idx = np.where(y >= backup_threshold)[0]
                meth = "backup_threshold"
            else:
                idx = np.where(best["posterior"] >= probability_threshold)[0]
                meth = "mixture"
            fits.append({"grna_index": gi, "method": meth, "n_nonzero": nnz, "pi": best["pi"], "g_pert": best["g"], "n_assigned": len(idx)})
        rows.extend([gi] * len(idx))
        cols.extend(idx.tolist())
    A = sp.csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n_g, n_c))
    return A, pd.DataFrame(fits)


SCEPTRE_R_TEMPLATE = r"""#!/usr/bin/env Rscript
# Written by igvfagent fishash-table2 (equivalent to upstream scripts/run_sceptre_table2.R).
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(sceptre))
if (as.character(packageVersion("sceptre")) != "0.10.3") stop("Table 2 requires sceptre 0.10.3; found ", packageVersion("sceptre"))
prefix <- "{prefix}"
features <- read.delim(gzfile(paste0(prefix, "_features.tsv.gz")), header = FALSE, col.names = c("ID", "Symbol", "type"), check.names = FALSE)
barcodes <- scan(gzfile(paste0(prefix, "_barcodes.tsv.gz")), what = character(), quiet = TRUE)
counts <- as(readMM(gzfile(paste0(prefix, "_matrix.mtx.gz"))), "CsparseMatrix")
rownames(counts) <- make.unique(features$Symbol); colnames(counts) <- barcodes
guide_counts <- counts[features$type == "CRISPR Guide Capture", , drop = FALSE]
gex_counts <- counts[features$type == "Gene Expression", , drop = FALSE]
batch <- sub("^.*/GSE272457_", "", prefix)
seed_key <- paste0("run_sceptre_mixture.R_", batch, ".Rds")
set.seed(abs(digest::digest2int(seed_key)))
so <- sceptre::import_data(response_matrix = gex_counts, grna_matrix = guide_counts,
  grna_target_data_frame = data.frame(grna_id = rownames(guide_counts), grna_target = "non-targeting"), moi = "high")
so <- sceptre::set_analysis_parameters(so)
cpus <- {cpus}L
so <- if (cpus == 1L) sceptre::assign_grnas(so, method = "mixture", parallel = FALSE) else sceptre::assign_grnas(so, method = "mixture", parallel = TRUE, n_processors = cpus)
a <- sceptre::get_grna_assignments(so)[rownames(guide_counts), ]
colnames(a) <- colnames(guide_counts)
writeMM(a, "{out}")
cat(sprintf("sceptre=%s sample=%s guides=%d cells=%d calls=%d seed_key=%s\n", packageVersion("sceptre"), batch, nrow(a), ncol(a), sum(a), seed_key))
"""


def cmd_sceptre_mixture(args: argparse.Namespace) -> int:
    np = _np()
    import scipy.io  # type: ignore
    work = Path(args.work_dir)
    out = Path(args.out) if args.out else work / f"{args.sample}_sceptre_mixture.mtx"
    prefix = args.raw_prefix or (str(Path(args.raw_dir) / SAMPLES.get(args.sample, args.sample)) if args.raw_dir else None)
    if args.engine == "r":
        if not prefix:
            raise SystemExit("--engine r needs --raw-prefix or --raw-dir")
        script = out.with_suffix(".R")
        script.write_text(SCEPTRE_R_TEMPLATE.replace("{prefix}", prefix).replace("{out}", str(out)).replace("{cpus}", str(args.cpus)))
        print(f"Wrote: {script}")
        cmd = ["Rscript", str(script)]
        if shutil.which("Rscript") is None:
            print("Rscript not on PATH; run with an R library holding sceptre 0.10.3:\n  R_LIBS_USER=/path/to/sceptre-lib " + " ".join(cmd))
            return 0
        return subprocess.call(cmd)
    if not prefix:
        raise SystemExit("sceptre-mixture needs the raw gene-expression counts: --raw-prefix or --raw-dir")
    features = read_features(Path(f"{prefix}_features.tsv.gz"))
    counts = read_mm(Path(f"{prefix}_matrix.mtx.gz"))
    is_guide = (features["type"] == "CRISPR Guide Capture").to_numpy()
    is_gex = (features["type"] == "Gene Expression").to_numpy()
    X = sceptre_covariates(counts[is_gex], counts[is_guide], grna_covariates=not args.no_grna_covariates)
    A, fits = sceptre_mixture_assign(counts[is_guide], X, probability_threshold=args.probability_threshold,
                                     n_em_rep=args.n_em_rep, seed=args.seed)
    out.parent.mkdir(parents=True, exist_ok=True)
    scipy.io.mmwrite(str(out), A.tocoo().astype(np.int64))
    print(f"Wrote: {out}")
    fits["guide"] = features.loc[is_guide, "Symbol"].to_numpy()
    fpath = out.with_suffix(".fits.csv")
    write_csv(fits, fpath)
    print(f"sceptre-mixture (python): guides={A.shape[0]} cells={A.shape[1]} calls={int(A.sum())} "
          f"mixture={int((fits['method'] == 'mixture').sum())} backup={int((fits['method'] == 'backup_threshold').sum())}")
    return 0


# ---------------------------------------------------------------------------
# scoring (score_table2.R)
# ---------------------------------------------------------------------------

def binary_accuracy_metrics(pred, truth) -> dict:
    """pred: float array with 1/0/nan; truth: bool array (same length)."""
    np = _np()
    pred = np.asarray(pred, dtype=float)
    truth = np.asarray(truth, dtype=bool)
    na = np.isnan(pred)
    p1 = (pred == 1) & ~na
    p0 = (pred == 0) & ~na
    tp, fp = int(np.sum(p1 & truth)), int(np.sum(p1 & ~truth))
    tn, fn = int(np.sum(p0 & ~truth)), int(np.sum(p0 & truth))
    na_t, na_f = int(np.sum(na & truth)), int(np.sum(na & ~truth))
    total = tp + fp + tn + fn + na_t + na_f
    tt, tf = tp + fn + na_t, tn + fp + na_f
    acc_t = tp / tt if tt else float("nan")
    acc_f = tn / tf if tf else float("nan")
    se = math.sqrt((tt * acc_t * (1 - acc_t) + tf * acc_f * (1 - acc_f)) / total ** 2) if total else float("nan")
    return {"accuracy": (tp + tn) / total if total else float("nan"), "stderr": se, "n": total,
            "both_or_neither": na_t + na_f, "tp": tp, "fp": fp, "tn": tn, "fn": fn}


def table2_cohort(meta) -> "tuple[Any, Any]":
    np = _np()
    sum_gex = (meta["homo_sum_gex"] + meta["mus_sum_gex"]).to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = meta["homo_sum_gex"].to_numpy(dtype=float) / sum_gex
        mito = meta["mito_sum"].to_numpy(dtype=float) / sum_gex
    fg = meta["features_gex"].to_numpy(dtype=float)
    selected = ((mito < 0.15) & (fg <= 6000) & (sum_gex <= 20000) & (fg >= 1500) & (sum_gex >= 3500)
                & ((frac < 0.1) | (frac > 0.9)))
    truth = (frac >= 0.1) & (sum_gex > 100)
    return np.nan_to_num(selected, nan=False).astype(bool), np.nan_to_num(truth, nan=False).astype(bool)


def read_mtx_assignment(path: Path):
    return (read_mm(path) > 0).tocsr()


def read_cleanser_assignment(path: Path, cutoff: float):
    """Banner-less CLEANSER coordinate output: first line dims, then guide cell posterior; keep >= cutoff."""
    np = _np()
    import scipy.sparse as sp  # type: ignore
    rows, cols = [], []
    dims = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("%") or not line.strip():
                continue
            parts = line.split()
            if dims is None:
                dims = (int(parts[0]), int(parts[1]))
                continue
            if float(parts[2]) >= cutoff:
                rows.append(int(parts[0]) - 1)
                cols.append(int(parts[1]) - 1)
    if dims is None:
        raise SystemExit(f"empty CLEANSER output: {path}")
    return sp.csr_matrix((np.ones(len(rows), dtype=bool), (rows, cols)), shape=dims)


def score_assignment(assignment, guide_type, truth, selected) -> dict:
    np = _np()
    gt = np.asarray(guide_type)
    if assignment.shape != (len(gt), len(truth)):
        raise SystemExit(f"assignment {assignment.shape} vs {len(gt)} guides x {len(truth)} cells")
    has_h = np.asarray(assignment[gt == "homo_guide"].sum(axis=0)).ravel() > 0
    has_m = np.asarray(assignment[gt == "mus_guide"].sum(axis=0)).ravel() > 0
    pred = has_h.astype(float)
    pred[(has_h & has_m) | (~has_h & ~has_m)] = np.nan
    return binary_accuracy_metrics(pred[selected], truth[selected])


def log_warning_counts(log_path: Path) -> "tuple[bool, int, int, int]":
    lines = log_path.read_text(errors="replace").splitlines() if log_path.is_file() else []
    complete = any("Random seed:" in x for x in lines)
    return (complete, sum("Non-fatal error during sampling" in x for x in lines),
            sum("Some chains may have failed to converge" in x for x in lines),
            sum("divergent transitions" in x for x in lines))


def score_work_dir(work: Path, samples: "list[str]", strict: bool = False) -> "tuple[Any, list[str]]":
    pd, np = _pd(), _np()
    out, notes = [], []
    for si, sample in enumerate(samples):
        meta = pd.read_csv(work / f"{sample}_meta.csv")
        guides = pd.read_csv(work / f"{sample}_guides.csv")
        selected, truth = table2_cohort(meta)
        col = SAMPLE_ORDER.index(sample) if sample in SAMPLE_ORDER else None
        for method, kind, mode, cutoff in SPECS:
            if kind == "cleanser":
                path = work / f"{sample}_cleanser_{mode}_posterior.mtx"
                logp = work / f"{sample}_cleanser_{mode}.log"
                if not path.is_file() or path.stat().st_size == 0:
                    msg = f"missing or empty CLEANSER output: {path}"
                    if strict:
                        raise SystemExit(msg)
                    notes.append(msg)
                    continue
                complete, nonfatal, conv, div = log_warning_counts(logp)
                if not complete:
                    msg = f"CLEANSER log lacks completion marker: {logp}"
                    if strict:
                        raise SystemExit(msg)
                    notes.append(msg)
                assignment = read_cleanser_assignment(path, cutoff)
            else:
                path = work / f"{sample}_sceptre_mixture.mtx"
                if not path.is_file() or path.stat().st_size == 0:
                    msg = f"missing or empty SCEPTRE output: {path}"
                    if strict:
                        raise SystemExit(msg)
                    notes.append(msg)
                    continue
                assignment = read_mtx_assignment(path)
                nonfatal = conv = div = np.nan
            r = score_assignment(assignment, guides["guide_type"].to_numpy(), truth, selected)
            exp = PUBLISHED[method][col] if (col is not None and method in PUBLISHED) else np.nan
            exp_se = PUBLISHED_SE[method][col] if (col is not None and method in PUBLISHED_SE) else np.nan
            acc, se = r["accuracy"], r["stderr"]
            has = not (isinstance(exp, float) and math.isnan(exp))
            out.append({"sample": sample, "method": method, "cutoff": cutoff if cutoff is not None else np.nan,
                        "accuracy": acc, "stderr": se, "n_qc": r["n"], "both_or_neither": r["both_or_neither"],
                        "nonfatal_sampler_warnings": nonfatal, "convergence_warnings": conv, "divergent_chain_warnings": div,
                        "published_accuracy": exp, "published_stderr": exp_se,
                        "accuracy_delta": acc - exp if has else np.nan, "stderr_delta": se - exp_se if has else np.nan,
                        "accuracy_match_4dp": (round(acc, 4) == exp) if has else None,
                        "stderr_match_4dp": (round(se, 4) == exp_se) if has else None,
                        "accuracy_within_0_001": (abs(acc - exp) <= 0.001) if has else None})
    return pd.DataFrame(out, columns=SCORED_COLS), notes


def detect_samples(work: Path) -> "list[str]":
    found = sorted(p.name[: -len("_meta.csv")] for p in work.glob("*_meta.csv"))
    ordered = [s for s in SAMPLE_ORDER if s in found]
    return ordered + [s for s in found if s not in ordered]


def plot_scores(df, path: Path) -> "Optional[Path]":
    plt = _plt()
    if plt is None or df.empty:
        return None
    np = _np()
    samples = list(dict.fromkeys(df["sample"]))
    methods = list(dict.fromkeys(df["method"]))
    fig, ax = plt.subplots(figsize=(max(6, 1.6 * len(samples)), 3.6))
    w = 0.8 / max(len(methods), 1)
    for i, m in enumerate(methods):
        sub = df[df["method"] == m].set_index("sample").reindex(samples)
        x = np.arange(len(samples)) + i * w
        ax.bar(x, sub["accuracy"], w, yerr=sub["stderr"], color=PALETTE[i % len(PALETTE)], label=m)
        pub = sub["published_accuracy"].to_numpy(dtype=float)
        ok = np.isfinite(pub)
        ax.scatter(x[ok], pub[ok], marker="_", s=200, color=INK, zorder=3)
    ax.set_xticks(np.arange(len(samples)) + 0.4 - w / 2)
    ax.set_xticklabels(samples, fontsize=8)
    ax.set_ylim(max(0, float(df["accuracy"].min()) - 0.1), 1.0)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    _style(ax, "Species-assignment accuracy (bars) vs Fishash Table 2 (ticks)", "", "accuracy")
    return _save(fig, path)


def cmd_score(args: argparse.Namespace) -> int:
    work = Path(args.work_dir)
    samples = args.samples or detect_samples(work)
    if not samples:
        raise SystemExit(f"no <sample>_meta.csv in {work}; run build-inputs first")
    df, notes = score_work_dir(work, samples, strict=args.strict)
    rd = run_dir(args.label)
    csvp = write_csv(df, rd / "extended_table2_reproduction.csv")
    figs = []
    if not args.no_plots:
        f = plot_scores(df, rd / "accuracy_vs_table2.png")
        if f:
            figs.append(str(f))
    pubrows = df[df["published_accuracy"].notna()]
    within = int(sum(v is True or v == True for v in pubrows["accuracy_within_0_001"].tolist()))  # noqa: E712
    summ = {"work_dir": str(work), "samples": samples, "n_rows": len(df), "n_published_rows": len(pubrows),
            "n_within_0_001": within, "notes": notes, "figures": figs, "csv": str(csvp),
            "rows": df.to_dict(orient="records")}
    write_json(summ, rd / "summary.json")
    rows = [(r["sample"], r["method"], _fmt(r["accuracy"]), _fmt(r["stderr"]), int(r["n_qc"]), int(r["both_or_neither"]),
             _fmt(r["published_accuracy"]), _fmt(r["accuracy_delta"], 5)) for _, r in df.iterrows()]
    text = (f"# Fishash Table 2 reproduction — {args.label}\n\n"
            f"Port of [{UPSTREAM_REPO}](https://github.com/{UPSTREAM_REPO}) @ `{UPSTREAM_COMMIT[:8]}`.\n\n"
            f"Work dir `{work}`; samples: {', '.join(samples)}.\n\n"
            "Cohort: mito/sum < 0.15, 1500 <= genes <= 6000, 3500 <= UMIs <= 20000, human fraction < 0.1 or > 0.9. "
            "Accuracy counts both-or-neither calls as wrong.\n\n"
            + md_table(["sample", "method", "accuracy", "stderr", "n_qc", "both/neither", "Table 2", "delta"], rows)
            + f"\n\n{within} of {len(pubrows)} rows with a published value are within 0.001 of Table 2.\n")
    if notes:
        text += "\nNotes:\n" + "\n".join(f"- {n}" for n in notes) + "\n"
    write_report(rd / "report.md", text)
    return 0


def cmd_compare(args: argparse.Namespace) -> int:
    pd, np = _pd(), _np()
    df = pd.read_csv(args.scored)
    rows = []
    for _, r in df.iterrows():
        col = SAMPLE_ORDER.index(r["sample"]) if r["sample"] in SAMPLE_ORDER else None
        pub = PUBLISHED.get(r["method"], [np.nan] * 4)[col] if col is not None else np.nan
        upr = UPSTREAM_REPRODUCED.get(r["method"], [np.nan] * 4)[col] if col is not None else np.nan
        rows.append({"sample": r["sample"], "method": r["method"], "accuracy": r["accuracy"],
                     "published_accuracy": pub, "upstream_reproduced_accuracy": upr,
                     "delta_vs_published": r["accuracy"] - pub, "delta_vs_upstream": r["accuracy"] - upr,
                     "within_0_001_published": abs(r["accuracy"] - pub) <= args.tolerance if np.isfinite(pub) else None,
                     "within_0_001_upstream": abs(r["accuracy"] - upr) <= args.tolerance if np.isfinite(upr) else None})
    cmp_df = pd.DataFrame(rows)
    rd = run_dir(args.label)
    write_csv(cmp_df, rd / "comparison.csv")
    n_pub = int(cmp_df["within_0_001_published"].dropna().astype(bool).sum())
    n_up = int(cmp_df["within_0_001_upstream"].dropna().astype(bool).sum())
    summ = {"scored": args.scored, "tolerance": args.tolerance, "n_rows": len(cmp_df),
            "within_published": n_pub, "n_published": int(cmp_df["within_0_001_published"].notna().sum()),
            "within_upstream": n_up, "n_upstream": int(cmp_df["within_0_001_upstream"].notna().sum())}
    write_json(summ, rd / "summary.json")
    text = (f"# Fishash Table 2 comparison — {args.label}\n\nScored table `{args.scored}`; tolerance {args.tolerance}.\n\n"
            + md_table(["sample", "method", "accuracy", "Table 2", "upstream repro", "within (T2)", "within (upstream)"],
                       [(r["sample"], r["method"], _fmt(r["accuracy"]), _fmt(r["published_accuracy"]),
                         _fmt(r["upstream_reproduced_accuracy"]), r["within_0_001_published"], r["within_0_001_upstream"])
                        for _, r in cmp_df.iterrows()])
            + f"\n\nWithin tolerance of Table 2: {n_pub}/{summ['n_published']}; of the upstream reproduction: {n_up}/{summ['n_upstream']}.\n")
    write_report(rd / "report.md", text)
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    work = Path(args.work_dir)
    pairs = sample_prefixes(args)
    for s, prefix in pairs:
        if not (work / f"{s}_meta.csv").is_file() or args.force:
            build_sample(prefix, s, work, args.human_max_index)
        if not args.skip_sceptre:
            out = work / f"{s}_sceptre_mixture.mtx"
            if not out.is_file() or args.force:
                cmd_sceptre_mixture(argparse.Namespace(work_dir=str(work), sample=s, raw_prefix=prefix, raw_dir=None, out=str(out),
                                                       engine=args.sceptre_engine, cpus=args.cpus, no_grna_covariates=False,
                                                       probability_threshold=0.8, n_em_rep=5, seed=args.seed))
        for mode in ("cs", "dc"):
            post = work / f"{s}_cleanser_{mode}_posterior.mtx"
            logp = work / f"{s}_cleanser_{mode}.log"
            if post.is_file() and post.stat().st_size and log_warning_counts(logp)[0] and not args.force:
                continue
            print(f"Running CLEANSER {mode} for {s}; sampler output goes to {logp}.")
            cmd_cleanser(argparse.Namespace(input=str(work / f"{s}_grna_counts.mtx"), mode=mode, out=str(post), log=str(logp),
                                            engine=args.cleanser_engine, method=args.cleanser_method, chains=args.chains,
                                            num_warmup=args.num_warmup, num_samples=args.num_samples, seed=args.seed,
                                            lpf=2, samples_out=None, threshold=None, cpus=args.cpus))
    return cmd_score(argparse.Namespace(work_dir=str(work), samples=[s for s, _ in pairs], strict=False,
                                        label=args.label, no_plots=args.no_plots))


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------

def synthetic_barnyard(d: Path, sample: str, n_cells: int = 700, n_genes_each: int = 2500, n_guides: int = 20,
                       direct_capture: bool = False, seed: int = 3) -> dict:
    """Human/mouse barnyard: species genes, one species-matched guide per cell with NB counts, ambient noise."""
    np = _np()
    import scipy.io  # type: ignore
    import scipy.sparse as sp  # type: ignore
    rng = np.random.default_rng(seed)
    species = rng.choice(["homo", "mus", "doublet"], size=n_cells, p=[0.47, 0.47, 0.06])
    genes = [f"GRCh38_G{i}" for i in range(n_genes_each - 5)] + [f"GRCh38_MT-{i}" for i in range(5)] + \
            [f"mm10___g{i}" for i in range(n_genes_each - 5)] + [f"mm10___mt-{i}" for i in range(5)]
    is_h = np.array([g.startswith("GRCh38") for g in genes])
    depth = rng.lognormal(math.log(9000), 0.35, size=n_cells)
    base = rng.gamma(2.0, 1.0, size=len(genes))
    base[["_MT-" in g or "_mt-" in g for g in genes]] *= 3
    gex = np.zeros((len(genes), n_cells))
    for c in range(n_cells):
        w = base * np.where(is_h, 1.0 if species[c] == "homo" else (0.02 if species[c] == "mus" else 0.5),
                            1.0 if species[c] == "mus" else (0.02 if species[c] == "homo" else 0.5))
        gex[:, c] = rng.poisson(depth[c] * w / w.sum())
    # guides: nt_1..nt_{n/2} human, the rest mouse (upstream: 100 + 100); one species-matched guide per cell
    h = n_guides // 2
    true_guide = np.array([rng.integers(1, h + 1) if s == "homo" else (rng.integers(h + 1, n_guides + 1) if s == "mus" else rng.integers(1, n_guides + 1))
                           for s in species])
    lib = rng.lognormal(0, 0.3, size=n_cells)
    grna = np.zeros((n_guides, n_cells))
    for c in range(n_cells):
        mu = 60 * lib[c]
        phi = 1.5
        grna[true_guide[c] - 1, c] = rng.negative_binomial(phi, phi / (phi + mu)) + 1
    noise_mean = 0.7  # ambient guide UMIs per guide per cell: most nonzero entries are noise, as in the real libraries
    amb = rng.poisson(noise_mean * (1.4 if direct_capture else 1.0) * lib[None, :], size=grna.shape)
    grna += amb
    names = [f"nt_{k}" for k in range(1, n_guides + 1)]
    feats = [(f"ENSG{i}", g, "Gene Expression") for i, g in enumerate(genes)] + [(n, n, "CRISPR Guide Capture") for n in names]
    barcodes = [f"CELL{c:05d}-1" for c in range(n_cells)]
    prefix = d / f"GSE272457_{sample}"
    with gzip.open(f"{prefix}_features.tsv.gz", "wt") as fh:
        fh.write("".join("\t".join(f) + "\n" for f in feats))
    with gzip.open(f"{prefix}_barcodes.tsv.gz", "wt") as fh:
        fh.write("".join(b + "\n" for b in barcodes))
    M = sp.vstack([sp.csr_matrix(gex), sp.csr_matrix(grna)]).tocoo().astype(np.int64)
    with gzip.open(f"{prefix}_matrix.mtx.gz", "wb") as fh:
        scipy.io.mmwrite(fh, M)
    return {"prefix": str(prefix), "species": species, "true_guide": true_guide, "grna": grna, "names": names}


def cmd_selftest(args: argparse.Namespace) -> int:
    pd, np = _pd(), _np()
    import tempfile
    checks = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    made = []
    t0 = time.time()
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        print("\nmetric definitions")
        pred = np.array([1, 1, 0, 0, np.nan, np.nan, 1, 0])
        truth = np.array([True, False, False, True, True, False, True, False])
        m = binary_accuracy_metrics(pred, truth)
        # tp=2 fp=1 tn=2 fn=1 na_t=1 na_f=1; T=4 F=4
        exp_se = math.sqrt((4 * 0.5 * 0.5 + 4 * 0.5 * 0.5) / 64)
        check(m["tp"] == 2 and m["tn"] == 2 and m["both_or_neither"] == 2 and abs(m["accuracy"] - 0.5) < 1e-12 and abs(m["stderr"] - exp_se) < 1e-12,
              "binary_accuracy_metrics: NA counted wrong, stderr = sqrt((T aT(1-aT) + F aF(1-aF)) / n^2)")
        meta = pd.DataFrame({"homo_sum_gex": [9000, 500, 4000, 9000, 20000, 9000], "mus_sum_gex": [100, 8000, 4000, 100, 5000, 100],
                             "mito_sum": [100, 100, 100, 2000, 100, 100], "features_gex": [3000, 3000, 3000, 3000, 3000, 1000]})
        sel, tr = table2_cohort(meta)
        check(list(sel) == [True, True, False, False, False, False] and list(tr) == [True, False, True, True, True, True],
              "table2_cohort: purity, mito < 0.15, UMI and gene windows; truth = human fraction >= 0.1")
        import scipy.sparse as sp  # type: ignore
        A = sp.csr_matrix(np.array([[1, 0, 1, 0], [0, 1, 1, 0]]))
        r = score_assignment(A, np.array(["homo_guide", "mus_guide"]), np.array([True, False, True, False]), np.array([True] * 4))
        check(r["tp"] == 1 and r["tn"] == 1 and r["both_or_neither"] == 2, "score_assignment: human-only -> human, mouse-only -> mouse, both/neither -> NA")
        rows = [(1, 1, 1), (2, 1, 5), (1, 2, 2), (3, 3, 30), (2, 3, 1)]
        L = cleanser_normalise(rows, 2)
        check(abs(L[1] - 1 / (4 / 3)) < 1e-12 and abs(L[2] - 2 / (4 / 3)) < 1e-12 and abs(L[3] - 1 / (4 / 3)) < 1e-12,
              "cleanser_normalise: sums of counts <= lpf per cell, divided by the mean over cells")

        print("\nCLEANSER model")
        rng = np.random.default_rng(1)
        n = 600
        Lc = rng.lognormal(0, 0.3, n)
        z = rng.uniform(size=n) < 0.06
        x = np.where(z, rng.negative_binomial(0.8, 0.8 / (0.8 + 80 * Lc)), rng.poisson(0.6, n))
        keep = x > 0
        mdl = CleanserModel("cs", x[keep], Lc[keep])
        u = cleanser_map(mdl)
        p_map = mdl.constrain(u)[0]
        pz = mdl.pzi(u)
        check(0.0 < p_map["lambda"] < 2.0 and p_map["nbMean"] > 20, f"cs MAP: lambda {p_map['lambda']:.2f}, nbMean {p_map['nbMean']:.1f}, r {p_map['r']:.3f}")
        hi_x = mdl.X >= 15
        check(pz[hi_x].min() > 0.9 and pz[mdl.X == 1].max() < 0.2, "cs posterior: counts >= 15 -> PZi > 0.9; count 1 -> PZi < 0.2")
        draws, dg = cleanser_mcmc(mdl, u, 150, 200, np.random.default_rng(2))
        check(draws.shape == (200, 4) and 0.05 < dg["acceptance"] < 0.8, f"cs MCMC: adaptive Metropolis acceptance {dg['acceptance']:.2f}")
        # truncated likelihood normalises: sum over x >= 1 of exp(loglik) == 1
        grid = CleanserModel("dc", np.arange(1, 4000), np.ones(3999))
        pdc = {"n_nbMean": 1.2, "n_nbDisp": 1.1, "r": 0.05, "nbMean": 100.0, "nbDisp": 3.0}
        lp0, lp1, ld = grid.components(pdc)
        check(abs(float(np.sum(np.exp(np.logaddexp(lp0, lp1) - ld))) - 1.0) < 1e-6, "dc likelihood: zero-truncated mixture sums to 1 over x >= 1")
        gcs = CleanserModel("cs", np.arange(1, 4000), np.ones(3999))
        lp0, lp1, ld = gcs.components({"lambda": 0.8, "r": 0.05, "nbMean": 100.0, "nbDisp": 0.5})
        check(abs(float(np.sum(np.exp(np.logaddexp(lp0, lp1) - ld))) - 1.0) < 1e-5, "cs likelihood: zero-truncated mixture sums to 1 over x >= 1")

        print("\nsynthetic barnyards: build-inputs")
        raw = d / "raw"
        raw.mkdir()
        work = d / "work"
        syn = {}
        for s, dc in (("mix0hr_Cropseq", False), ("mix0hr_DirectCapture", True)):
            syn[s] = synthetic_barnyard(raw, s, direct_capture=dc, seed=5 if dc else 3)
        rc = cmd_build_inputs(argparse.Namespace(raw_dir=None, raw_prefix=syn["mix0hr_Cropseq"]["prefix"], sample="mix0hr_Cropseq",
                                                 samples=None, work_dir=str(work), human_max_index=10))
        rc2 = cmd_build_inputs(argparse.Namespace(raw_dir=None, raw_prefix=syn["mix0hr_DirectCapture"]["prefix"], sample="mix0hr_DirectCapture",
                                                  samples=None, work_dir=str(work), human_max_index=10))
        meta = pd.read_csv(work / "mix0hr_Cropseq_meta.csv")
        gdf = pd.read_csv(work / "mix0hr_Cropseq_guides.csv")
        check(rc == 0 and rc2 == 0 and list(meta.columns) == ["barcode", "homo_sum_gex", "mus_sum_gex", "mito_sum", "features_gex"]
              and (meta["mito_sum"] > 0).all(), "build-inputs: meta.csv columns; mito from GRCh38_MT- / mm10___mt- genes")
        check(list(gdf.columns) == ["guide", "gidx", "guide_type"] and (gdf["guide_type"] == "homo_guide").sum() == 10,
              "build-inputs: guides.csv, nt_k with k <= human-max-index -> homo_guide")
        hdr, mrows = read_cleanser_mtx(work / "mix0hr_Cropseq_grna_counts.mtx")
        check(hdr.split()[:2] == ["20", "700"] and len(mrows) == int((syn["mix0hr_Cropseq"]["grna"] > 0).sum()),
              "build-inputs: guide x cell Matrix Market for CLEANSER (1-based guide, cell, count)")

        print("\ncleanser (python)")
        for s in syn:
            for mode in ("cs", "dc"):
                rc = cmd_cleanser(argparse.Namespace(input=str(work / f"{s}_grna_counts.mtx"), mode=mode,
                                                     out=str(work / f"{s}_cleanser_{mode}_posterior.mtx"), log=None,
                                                     engine="python", method="mcmc" if mode == "cs" and s == "mix0hr_Cropseq" else "map",
                                                     chains=1, num_warmup=100, num_samples=150, seed=20260810, lpf=2,
                                                     samples_out=None, threshold=0.8, cpus=1))
                check(rc == 0, f"cleanser {mode} {s}: posterior written")
        comp, *_ = log_warning_counts(work / "mix0hr_Cropseq_cleanser_cs.log")
        check(comp, "cleanser: log carries the 'Random seed:' completion marker score_table2.R checks")
        post = read_cleanser_assignment(work / "mix0hr_Cropseq_cleanser_cs_posterior.mtx", 0.8)
        tg = syn["mix0hr_Cropseq"]["true_guide"]
        hit = np.mean([post[tg[c] - 1, c] for c in range(len(tg))])
        check(hit > 0.9, f"cleanser cs: planted guide called in {hit:.1%} of cells at posterior >= 0.8")
        amb_calls = post.sum() - sum(post[tg[c] - 1, c] for c in range(len(tg)))
        check(amb_calls < 0.05 * len(tg), f"cleanser cs: ambient noise mostly rejected ({int(amb_calls)} off-target calls)")

        print("\nsceptre-mixture (python)")
        for s in syn:
            rc = cmd_sceptre_mixture(argparse.Namespace(work_dir=str(work), sample=s, raw_prefix=syn[s]["prefix"], raw_dir=None, out=None,
                                                        engine="python", cpus=1, no_grna_covariates=False,
                                                        probability_threshold=0.8, n_em_rep=5, seed=4))
            check(rc == 0 and (work / f"{s}_sceptre_mixture.mtx").is_file(), f"sceptre-mixture {s}: assignment written")
        sa = read_mtx_assignment(work / "mix0hr_Cropseq_sceptre_mixture.mtx")
        hit_s = np.mean([sa[tg[c] - 1, c] for c in range(len(tg))])
        check(hit_s > 0.85, f"sceptre-mixture: planted guide assigned in {hit_s:.1%} of cells")
        fits = pd.read_csv(work / "mix0hr_Cropseq_sceptre_mixture.fits.csv")
        check((fits["method"] == "mixture").mean() > 0.8 and (fits.loc[fits["method"] == "mixture", "g_pert"] > 2).all(),
              "sceptre-mixture: EM used for guides with >= 10 nonzero cells; perturbation effect g > 0")
        rc = cmd_sceptre_mixture(argparse.Namespace(work_dir=str(work), sample="mix0hr_Cropseq", raw_prefix=syn["mix0hr_Cropseq"]["prefix"],
                                                    raw_dir=None, out=str(d / "r_out.mtx"), engine="r", cpus=2, no_grna_covariates=False,
                                                    probability_threshold=0.8, n_em_rep=5, seed=4))
        rtxt = (d / "r_out.R").read_text()
        check(rc == 0 and 'method = "mixture"' in rtxt and "digest2int" in rtxt, "sceptre-mixture --engine r: sceptre 0.10.3 script written (runs only with Rscript)")

        print("\nscore")
        rc = cmd_score(argparse.Namespace(work_dir=str(work), samples=None, strict=True, label="st_fishash", no_plots=args.no_plots))
        rd = sorted(OUT_ROOT.glob("*_st_fishash"))[-1]; made.append(rd)
        sc = pd.read_csv(rd / "extended_table2_reproduction.csv")
        check(rc == 0 and list(sc.columns) == SCORED_COLS and len(sc) == 10, "score: 5 methods x 2 samples in the upstream CSV schema")
        cs = sc[(sc["sample"] == "mix0hr_Cropseq")].set_index("method")
        check(cs.loc["cleanser cs (Table 2, 0.80)", "accuracy"] > 0.9 and cs.loc["sceptre mixture", "accuracy"] > 0.9,
              f"score: planted barnyard accuracies cs {cs.loc['cleanser cs (Table 2, 0.80)', 'accuracy']:.3f}, sceptre {cs.loc['sceptre mixture', 'accuracy']:.3f}")
        n80 = int(read_cleanser_assignment(work / "mix0hr_Cropseq_cleanser_cs_posterior.mtx", 0.80).sum())
        n95 = int(read_cleanser_assignment(work / "mix0hr_Cropseq_cleanser_cs_posterior.mtx", 0.95).sum())
        check(n95 <= n80 and cs.loc["cleanser cs (public-code 0.95)", "cutoff"] == 0.95,
              f"score: public-code 0.95 cut-off scored alongside Table 2's 0.80 ({n95} vs {n80} calls)")
        n_sel = int(table2_cohort(pd.read_csv(work / "mix0hr_Cropseq_meta.csv"))[0].sum())
        check((cs["n_qc"] == n_sel).all() and n_sel < 700, f"score: n_qc = cohort size {n_sel} (doublets and impure cells removed)")
        check(abs(cs.loc["sceptre mixture", "published_accuracy"] - 0.9491) < 1e-12 and
              abs(cs.loc["sceptre mixture", "accuracy_delta"] - (cs.loc["sceptre mixture", "accuracy"] - 0.9491)) < 1e-12,
              "score: published Table 2 value and delta attached per method x sample")
        summ = json.loads((rd / "summary.json").read_text())
        check((rd / "report.md").is_file() and summ["n_rows"] == 10, "score: report.md + summary.json")

        print("\ncompare")
        rc = cmd_compare(argparse.Namespace(scored=str(rd / "extended_table2_reproduction.csv"), tolerance=0.001, label="st_fishash_cmp"))
        rc2 = sorted(OUT_ROOT.glob("*_st_fishash_cmp"))[-1]; made.append(rc2)
        up = pd.DataFrame({"sample": [s for s in SAMPLE_ORDER for _ in SPECS], "method": [m for _ in SAMPLE_ORDER for m, *_ in SPECS]})
        up["accuracy"] = [UPSTREAM_REPRODUCED[m][SAMPLE_ORDER.index(s)] for s, m in zip(up["sample"], up["method"])]
        up_path = d / "upstream.csv"
        up.to_csv(up_path, index=False)
        rc3 = cmd_compare(argparse.Namespace(scored=str(up_path), tolerance=0.001, label="st_fishash_up"))
        rd3 = sorted(OUT_ROOT.glob("*_st_fishash_up"))[-1]; made.append(rd3)
        s3 = json.loads((rd3 / "summary.json").read_text())
        check(rc == 0 and rc3 == 0 and s3["within_published"] == 12 and s3["within_upstream"] == 20,
              "compare: the upstream repository's own 20 accuracies -> all 12 published rows within 0.001 (as its README states)")

        print("\nrun (end to end)")
        rc = cmd_run(argparse.Namespace(raw_dir=None, raw_prefix=syn["mix0hr_Cropseq"]["prefix"], sample="mix0hr_Cropseq", samples=None,
                                        work_dir=str(d / "work_run"), human_max_index=10, force=False, skip_sceptre=False,
                                        sceptre_engine="python", cleanser_engine="python", cleanser_method="map", chains=1,
                                        num_warmup=100, num_samples=100, seed=1, cpus=1, label="st_fishash_run", no_plots=True))
        rr = sorted(OUT_ROOT.glob("*_st_fishash_run"))[-1]; made.append(rr)
        check(rc == 0 and len(pd.read_csv(rr / "extended_table2_reproduction.csv")) == 5, "run: build -> sceptre -> cleanser cs/dc -> score")
        rc = cmd_cleanser(argparse.Namespace(input=str(work / "mix0hr_Cropseq_grna_counts.mtx"), mode="cs", out=str(d / "bin_post.mtx"), log=None,
                                             engine="cleanser", method="map", chains=1, num_warmup=1, num_samples=1, seed=1, lpf=2,
                                             samples_out=None, threshold=None, cpus=1))
        check(rc == 0, "cleanser --engine cleanser: runs the binary or prints the exact upstream command")
    for p in made:
        shutil.rmtree(p, ignore_errors=True)
    try:
        OUT_ROOT.rmdir()  # only if the selftest left it empty
    except OSError:
        pass
    ok = all(c for c, _ in checks)
    print(f"\n({time.time() - t0:.0f} s)")
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[list[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent fishash-table2", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add_samples(p):
        p.add_argument("--raw-dir", help="directory with the GSE272457 <prefix>_{features.tsv,barcodes.tsv,matrix.mtx}.gz files")
        p.add_argument("--raw-prefix", help="one raw 10x triplet prefix (instead of --raw-dir)")
        p.add_argument("--sample", help="sample name for --raw-prefix")
        p.add_argument("--samples", nargs="+", help=f"samples to process with --raw-dir (default: {' '.join(SAMPLE_ORDER)})")
        p.add_argument("--human-max-index", type=int, default=100, help="nt_<k> guides with k <= this are human")

    p = sub.add_parser("build-inputs", help="raw barnyard 10x files -> guide mtx / h5ad / meta.csv / guides.csv")
    add_samples(p)
    p.add_argument("--work-dir", required=True)
    p.set_defaults(func=cmd_build_inputs)

    p = sub.add_parser("cleanser", help="CLEANSER cs/dc mixture posteriors (python MCMC/MAP, or the real binary)")
    p.add_argument("--input", required=True, help="guide x cell Matrix Market (CLEANSER input)")
    p.add_argument("--mode", choices=["cs", "dc"], required=True, help="cs = CROP-seq model, dc = direct-capture model")
    p.add_argument("--out", required=True, help="posterior output (CLEANSER coordinate format)")
    p.add_argument("--log", help="log file (default: <out without _posterior.mtx>.log)")
    p.add_argument("--engine", choices=["python", "cleanser"], default="python")
    p.add_argument("--method", choices=["mcmc", "map"], default="mcmc")
    p.add_argument("--chains", type=int, default=4)
    p.add_argument("--num-warmup", type=int, default=300)
    p.add_argument("--num-samples", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20260810)
    p.add_argument("--lpf", type=int, default=2, help="normalisation low-pass: only counts <= lpf enter the library size (0 = all)")
    p.add_argument("--samples-output", dest="samples_out", help="write the parameter draws (CLEANSER --so)")
    p.add_argument("--threshold", type=float, help="also report the number of calls at this posterior")
    p.add_argument("--cpus", type=int, default=8)
    p.set_defaults(func=cmd_cleanser)

    p = sub.add_parser("sceptre-mixture", help="SCEPTRE-style mixture guide assignment (python, or the sceptre R package)")
    p.add_argument("--work-dir", default=".")
    p.add_argument("--sample", default="sample")
    p.add_argument("--raw-prefix")
    p.add_argument("--raw-dir")
    p.add_argument("--out", help="assignment .mtx (default: <work-dir>/<sample>_sceptre_mixture.mtx)")
    p.add_argument("--engine", choices=["python", "r"], default="python")
    p.add_argument("--cpus", type=int, default=8)
    p.add_argument("--no-grna-covariates", action="store_true", help="drop log(grna_n_nonzero+1), log(grna_n_umis+1) from the formula")
    p.add_argument("--probability-threshold", type=float, default=0.8)
    p.add_argument("--n-em-rep", type=int, default=5)
    p.add_argument("--seed", type=int, default=4)
    p.set_defaults(func=cmd_sceptre_mixture)

    p = sub.add_parser("score", help="Table 2 cohort + species accuracy vs the published Table 2")
    p.add_argument("--work-dir", required=True)
    p.add_argument("--samples", nargs="+")
    p.add_argument("--strict", action="store_true", help="stop on missing outputs / logs like score_table2.R")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="table2")
    p.set_defaults(func=cmd_score)

    p = sub.add_parser("compare", help="scored CSV vs published Table 2 and the upstream repository's reproduction")
    p.add_argument("--scored", required=True, help="CSV with sample, method, accuracy")
    p.add_argument("--tolerance", type=float, default=0.001)
    p.add_argument("--label", default="table2_compare")
    p.set_defaults(func=cmd_compare)

    p = sub.add_parser("run", help="build-inputs -> sceptre-mixture -> cleanser cs/dc -> score")
    add_samples(p)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--force", action="store_true")
    p.add_argument("--skip-sceptre", action="store_true")
    p.add_argument("--sceptre-engine", choices=["python", "r"], default="python")
    p.add_argument("--cleanser-engine", choices=["python", "cleanser"], default="python")
    p.add_argument("--cleanser-method", choices=["mcmc", "map"], default="mcmc")
    p.add_argument("--chains", type=int, default=4)
    p.add_argument("--num-warmup", type=int, default=300)
    p.add_argument("--num-samples", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20260810)
    p.add_argument("--cpus", type=int, default=8)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="table2")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("selftest", help="synthetic barnyards, all subcommands, assertions")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
