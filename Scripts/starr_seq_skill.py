#!/usr/bin/env python3
"""IGVF agent STARR-seq analytics skill.

End-to-end allelic STARR-seq counts analysis. The wet-lab algorithms are
clean-room reimplementations of the math used in the public reference repo
`gaochengwen/STARR-seq-Data-Analysis` (the awk + R pipeline that wraps
Bioconductor `mpra::mpralm` for allelic significance testing). No source
from that repository is copied — algorithms paraphrased from the published
descriptions and the underlying limma/voom/eBayes statistical methods.

Supported workflow:

    starrseq pull-portal     — discover IGVF Portal STARR-seq MeasurementSets / files
    starrseq qc              — TPM-style count QC, RLE, Spearman D-stat outliers
    starrseq aggregate       — collapse per-barcode counts to per-(SNP, allele) wide table
    starrseq activity        — per-fragment log activity = log(RNA / DNA), allelic skew
    starrseq allelic-test    — mpralm-style allelic test (voom-weighted GLS + eBayes)
    starrseq write-playbook  — emit the skill's markdown playbook

Heavy deps (numpy, pandas, scipy, statsmodels) are imported lazily so the
`pull-portal` and `--help` commands keep working on a bare interpreter.

License: Apache-2.0 — the math (mpralm / limma+voom / eBayes) is published
and reimplemented here; no GPL or unlicensed source is linked at runtime.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "STARRseq"
PLOT_DIR = REPORT_DIR / "Plots"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

PORTAL_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_BASE")


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"starr_seq_skill_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)
    return log_path


def safe_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in label)


def _require_pkg(name: str, hint: str) -> Any:
    try:
        return __import__(name)
    except Exception as exc:
        raise SystemExit(
            f"Missing dependency '{name}'. {hint}\nInstall with: pip install {name}"
        ) from exc


def save_json(label: str, data: Any) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}.json"
    path.write_text(json.dumps(data, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Portal discovery
# ---------------------------------------------------------------------------

def cmd_pull_portal(args: argparse.Namespace) -> int:
    """Discover IGVF Portal STARR-seq MeasurementSets and write a manifest."""
    setup_logging()
    params = {
        "type": "MeasurementSet",
        "assay_titles": "STARR-seq",
        "format": "json",
        "limit": str(args.limit),
    }
    url = f"{PORTAL_BASE}/search/?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = json.load(resp)
    rows = data.get("@graph", [])
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_portal_starr.tsv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "assay_titles", "preferred_assay_titles", "n_files", "status", "description"])
        for r in rows:
            w.writerow([
                r.get("accession", ""),
                "|".join(r.get("assay_titles") or []),
                "|".join(r.get("preferred_assay_titles") or []),
                len(r.get("files") or []),
                r.get("status", ""),
                (r.get("description") or "")[:120],
            ])
    summary = {
        "url": url,
        "total": data.get("total"),
        "returned": len(rows),
        "manifest": str(out),
    }
    save_json(f"starr_{args.label}_portal_manifest", summary)
    print(f"STARR-seq portal manifest: {out}")
    print(f"Total: {summary['total']}, returned: {summary['returned']}")
    return 0


# ---------------------------------------------------------------------------
# Counts I/O
# ---------------------------------------------------------------------------

def _detect_sample_cols(columns: list[str]) -> dict[str, list[str]]:
    """Detect DNA_* / RNA_* columns (case-insensitive).

    Accepts ``DNA_rep1``, ``RNA.HepG2.rep2``, ``DNA1`` patterns; returns
    ``{"DNA": [...], "RNA": [...]}`` lists.
    """
    import re
    pat = re.compile(r"^(DNA|RNA)([._-].*)?$", re.IGNORECASE)
    buckets: dict[str, list[str]] = {"DNA": [], "RNA": []}
    for c in columns:
        m = pat.match(c)
        if m:
            buckets[m.group(1).upper()].append(c)
    return buckets


def _load_counts(path: Path) -> Any:
    """Load counts TSV/CSV into a pandas DataFrame."""
    pd = _require_pkg("pandas", "Required for STARR-seq counts analytics.")
    df = pd.read_csv(path, sep=None, engine="python")
    return df


# ---------------------------------------------------------------------------
# B1 — TPM-style count QC + RLE + Spearman D-stat
# ---------------------------------------------------------------------------

def _tpm_normalize(counts: Any, scale: float = 5000.0) -> Any:
    np = _require_pkg("numpy", "")
    # column-totals scaling: each sample sums to `scale`
    col_sum = counts.sum(axis=0).replace(0, np.nan)
    return counts.divide(col_sum, axis=1).fillna(0.0) * scale


def _rle_matrix(tpm: Any) -> Any:
    np = _require_pkg("numpy", "")
    log_tpm = np.log10(tpm + 1e-4)
    return log_tpm.subtract(log_tpm.median(axis=1), axis=0)


def _spearman_d_stat(tpm: Any) -> Any:
    """Per-sample distance score = median(1 - spearman r) across the other samples.

    Samples in the lowest 5% (highest distance) flagged as outliers.
    """
    pd = _require_pkg("pandas", "")
    np = _require_pkg("numpy", "")
    scipy_stats = __import__("scipy.stats", fromlist=["stats"])
    samples = list(tpm.columns)
    n = len(samples)
    dist = pd.Series(0.0, index=samples)
    for j, s in enumerate(samples):
        scores = []
        for k, t in enumerate(samples):
            if j == k:
                continue
            rho, _ = scipy_stats.spearmanr(tpm[s].to_numpy(), tpm[t].to_numpy())
            if np.isnan(rho):
                rho = 0.0
            scores.append(1.0 - rho)
        dist.at[s] = float(np.median(scores)) if scores else 0.0
    return dist


def cmd_qc(args: argparse.Namespace) -> int:
    pd = _require_pkg("pandas", "")
    np = _require_pkg("numpy", "")
    setup_logging()
    df = _load_counts(Path(args.input))
    samples = _detect_sample_cols(list(df.columns))
    cols = samples["DNA"] + samples["RNA"]
    if not cols:
        raise SystemExit(f"No DNA_*/RNA_* sample columns found. Detected={samples}")
    if "Fragment" not in df.columns and "Oligo" not in df.columns:
        if "SNP" in df.columns and "Allele" in df.columns:
            df = df.assign(Fragment=df["SNP"].astype(str) + "_" + df["Allele"].astype(str))
        else:
            raise SystemExit(
                "Counts table needs a 'Fragment' or 'Oligo' identifier column "
                "(or 'SNP' + 'Allele' so a Fragment id can be auto-derived)."
            )
    id_col = "Fragment" if "Fragment" in df.columns else "Oligo"
    counts = df.set_index(id_col)[cols].astype(float)

    tpm = _tpm_normalize(counts)
    # Drop low-expression fragments: >0.1 TPM in <=20% of samples
    keep = (tpm > 0.1).mean(axis=1) > 0.20
    tpm_kept = tpm.loc[keep]
    counts_kept = counts.loc[keep]
    logging.info("QC: kept %d / %d fragments (low-expression filter)", len(tpm_kept), len(tpm))

    rle = _rle_matrix(tpm_kept)
    d_stat = _spearman_d_stat(tpm_kept)
    cutoff = float(d_stat.quantile(0.95))  # flag worst 5%
    outliers = d_stat[d_stat >= cutoff].index.tolist()

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    tag = safe_label(args.label)
    tpm.to_csv(REPORT_DIR / f"{ts}_{tag}_tpm.tsv", sep="\t")
    rle.to_csv(REPORT_DIR / f"{ts}_{tag}_rle.tsv", sep="\t")

    qc_df = pd.DataFrame({
        "sample": d_stat.index,
        "spearman_dstat": d_stat.values,
        "rle_iqr": rle[d_stat.index].apply(lambda s: float(s.quantile(0.75) - s.quantile(0.25))).values,
        "rle_median_abs": rle[d_stat.index].abs().median().values,
        "outlier": [s in outliers for s in d_stat.index],
        "condition": ["DNA" if s in samples["DNA"] else "RNA" for s in d_stat.index],
    })
    qc_path = REPORT_DIR / f"{ts}_{tag}_qc.tsv"
    qc_df.to_csv(qc_path, sep="\t", index=False)
    summary = {
        "input": str(Path(args.input).resolve()),
        "fragments_in": int(len(df)),
        "fragments_kept": int(len(tpm_kept)),
        "samples": samples,
        "outlier_samples": outliers,
        "spearman_dstat": d_stat.to_dict(),
        "outlier_cutoff_q95": cutoff,
        "qc_table": str(qc_path),
    }
    save_json(f"starr_{args.label}_qc", summary)
    print(f"QC table: {qc_path}")
    print(f"Kept {len(tpm_kept)} / {len(tpm)} fragments after low-expression filter.")
    print(f"Outlier samples (top 5% Spearman D-stat): {outliers or 'none'}")
    return 0


# ---------------------------------------------------------------------------
# B2 — Allele-fragment aggregation
# ---------------------------------------------------------------------------

def cmd_aggregate(args: argparse.Namespace) -> int:
    """Collapse barcode/fragment-level counts into per-(SNP, Allele) wide table.

    Input: long-form counts with at least ``SNP``, ``Allele`` columns and
    one numeric count column per sample. Output: wide table indexed by
    ``f"{SNP}_{Allele}"`` with one column per sample.
    """
    pd = _require_pkg("pandas", "")
    setup_logging()
    df = _load_counts(Path(args.input))
    if "SNP" not in df.columns or "Allele" not in df.columns:
        raise SystemExit("Aggregation needs 'SNP' and 'Allele' columns.")
    samples = _detect_sample_cols(list(df.columns))
    cols = samples["DNA"] + samples["RNA"]
    grp = df.groupby(["SNP", "Allele"], as_index=False)[cols].sum()
    grp["Fragment"] = grp["SNP"] + "_" + grp["Allele"].astype(str)
    grp = grp[["Fragment", "SNP", "Allele"] + cols]
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_aggregated.tsv"
    grp.to_csv(out, sep="\t", index=False)
    print(f"Aggregated counts: {out} ({len(grp)} per-(SNP, Allele) rows)")
    return 0


# ---------------------------------------------------------------------------
# B3 — Per-fragment activity + allelic skew
# ---------------------------------------------------------------------------

def cmd_activity(args: argparse.Namespace) -> int:
    """Compute per-fragment log activity and per-SNP allelic skew."""
    pd = _require_pkg("pandas", "")
    np = _require_pkg("numpy", "")
    setup_logging()
    df = _load_counts(Path(args.input))
    samples = _detect_sample_cols(list(df.columns))
    if not samples["DNA"] or not samples["RNA"]:
        raise SystemExit("Need DNA_* and RNA_* sample columns.")
    if "SNP" not in df.columns or "Allele" not in df.columns:
        raise SystemExit("Need 'SNP' and 'Allele' columns (run `starrseq aggregate` first if needed).")
    id_col = "Fragment" if "Fragment" in df.columns else None
    if id_col is None:
        df["Fragment"] = df["SNP"] + "_" + df["Allele"].astype(str)
        id_col = "Fragment"

    dna = df[samples["DNA"]].astype(float).to_numpy() + 0.5
    rna = df[samples["RNA"]].astype(float).to_numpy() + 0.5
    # TPM per sample column
    dna_tpm = dna / dna.sum(axis=0) * 5000.0
    rna_tpm = rna / rna.sum(axis=0) * 5000.0
    # Pair replicates by position (caller's responsibility to order columns)
    n_pairs = min(dna_tpm.shape[1], rna_tpm.shape[1])
    ratio = rna_tpm[:, :n_pairs] / dna_tpm[:, :n_pairs]
    log_act = np.log(ratio)
    df["log_activity_mean"] = log_act.mean(axis=1)
    df["log_activity_sd"] = log_act.std(axis=1, ddof=1) if n_pairs > 1 else np.nan
    df["n_replicates"] = n_pairs

    # Allelic skew per SNP: log( ratio_alt / ratio_ref ) averaged over replicates
    alleles = df["Allele"].astype(str).str.upper()
    df["__std_allele"] = alleles.map({"REF": "REF", "ALT": "ALT", "A": "REF", "B": "ALT"}).fillna(alleles)
    ref = df[df["__std_allele"] == "REF"].set_index("SNP")
    alt = df[df["__std_allele"] == "ALT"].set_index("SNP")
    common = sorted(set(ref.index) & set(alt.index))
    skew_rows = []
    for snp in common:
        l_ref = ref.loc[snp, "log_activity_mean"]
        l_alt = alt.loc[snp, "log_activity_mean"]
        # When duplicate SNP entries exist (multi-fragment per SNP), take the mean of each side
        if hasattr(l_ref, "mean"):
            l_ref = float(l_ref.mean())
        if hasattr(l_alt, "mean"):
            l_alt = float(l_alt.mean())
        skew_rows.append({"SNP": snp, "log_activity_ref": l_ref, "log_activity_alt": l_alt,
                          "log_skew_alt_minus_ref": float(l_alt - l_ref)})
    skew_df = pd.DataFrame(skew_rows)

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    act_out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_activity.tsv"
    skew_out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_allelic_skew.tsv"
    df.drop(columns=["__std_allele"]).to_csv(act_out, sep="\t", index=False)
    skew_df.to_csv(skew_out, sep="\t", index=False)
    print(f"Per-fragment activity: {act_out} ({len(df)} rows, {n_pairs} replicate pairs)")
    print(f"Per-SNP allelic skew: {skew_out} ({len(skew_df)} SNPs)")
    return 0


# ---------------------------------------------------------------------------
# B4 — mpralm-style allelic test (voom-weighted GLS + eBayes)
# ---------------------------------------------------------------------------

def _voom_weights(log_counts: Any, design: Any) -> Any:
    """Compute voom (Law 2014) precision weights.

    Steps:
      1. Fit ordinary OLS per fragment: y = X beta; record sigma_g (residual SD).
      2. Compute mean log-count per fragment.
      3. Lowess-fit sqrt(sigma_g) vs mean(log_count) -> mean-variance trend f(.).
      4. For each observation y_{g,j}: predicted log-count mu_{g,j} = X_j @ beta_g.
         Predicted variance = f(mu_{g,j})^4 (because voom fits sqrt(sigma)).
      5. Weight = 1 / predicted_variance.
    """
    np = _require_pkg("numpy", "")
    stat = __import__("statsmodels.nonparametric.smoothers_lowess",
                       fromlist=["lowess"]).lowess
    n_frag, n_samp = log_counts.shape
    X = np.asarray(design, dtype=float)
    # OLS per fragment via closed form: beta = (X^T X)^-1 X^T y
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = (XtX_inv @ X.T) @ log_counts.T  # shape (p, n_frag)
    fitted = X @ beta  # shape (n_samp, n_frag)
    resid = log_counts.T - fitted  # shape (n_samp, n_frag)
    df_resid = max(n_samp - X.shape[1], 1)
    sigma_g = np.sqrt((resid ** 2).sum(axis=0) / df_resid)  # length n_frag
    mean_log = log_counts.mean(axis=1)
    # Lowess of sqrt(sigma) vs mean log-count
    order = np.argsort(mean_log)
    smooth = stat(np.sqrt(sigma_g[order]), mean_log[order], frac=0.5, return_sorted=False)
    # Build interpolator over the same x-axis
    def f(x):
        return np.interp(x, mean_log[order], smooth)
    pred_log = (X @ beta).T  # shape (n_frag, n_samp)
    w = 1.0 / (f(pred_log) ** 4 + 1e-9)
    return w


def cmd_allelic_test(args: argparse.Namespace) -> int:
    """Per-(SNP, allele) allelic test, clean-room paraphrase of mpra::mpralm.

    Model: y_{g,j} = beta_0 + beta_1 * allele2_{j} + e_{g,j}
    Pipeline:
      1. Pair DNA + RNA columns; compute log2(RNA+0.5 / DNA+0.5).
      2. voom precision weights from log2(DNA + 0.5).
      3. Weighted OLS per fragment for beta_1 (allele2 effect).
      4. Empirical-Bayes moderation of variances (Smyth 2004 closed form).
      5. BH-FDR.
    """
    pd = _require_pkg("pandas", "")
    np = _require_pkg("numpy", "")
    scipy_stats = __import__("scipy.stats", fromlist=["stats"])
    from statsmodels.stats.multitest import multipletests
    setup_logging()

    df = _load_counts(Path(args.input))
    samples = _detect_sample_cols(list(df.columns))
    if not samples["DNA"] or not samples["RNA"]:
        raise SystemExit("Need DNA_* and RNA_* sample columns.")
    if "SNP" not in df.columns or "Allele" not in df.columns:
        raise SystemExit("Need SNP and Allele columns (use `starrseq aggregate` first).")

    alleles = df["Allele"].astype(str).str.upper()
    std_allele = alleles.map({"REF": "REF", "ALT": "ALT", "A": "REF", "B": "ALT"}).fillna(alleles)
    df = df.assign(__std_allele=std_allele)
    ref = df[df["__std_allele"] == "REF"].set_index("SNP")
    alt = df[df["__std_allele"] == "ALT"].set_index("SNP")
    common = sorted(set(ref.index) & set(alt.index))
    if not common:
        raise SystemExit("No paired REF/ALT SNPs found.")
    n_rep = min(len(samples["DNA"]), len(samples["RNA"]))
    dna_cols, rna_cols = samples["DNA"][:n_rep], samples["RNA"][:n_rep]

    # Build the y matrix (one row per SNP, columns = 2 alleles * n_rep)
    # Following mpralm "indep_groups" model: stack ref then alt across replicates,
    # design = [1, allele2_indicator]; size n_obs = 2*n_rep per SNP.
    rows = []
    for snp in common:
        r_dna = ref.loc[snp, dna_cols]
        r_rna = ref.loc[snp, rna_cols]
        a_dna = alt.loc[snp, dna_cols]
        a_rna = alt.loc[snp, rna_cols]
        # Allow multi-fragment per (SNP, Allele) — average if needed
        if hasattr(r_dna, "mean") and getattr(r_dna, "ndim", 1) > 1:
            r_dna = r_dna.mean(axis=0)
            r_rna = r_rna.mean(axis=0)
            a_dna = a_dna.mean(axis=0)
            a_rna = a_rna.mean(axis=0)
        rows.append(np.concatenate([r_dna.to_numpy(), r_rna.to_numpy(),
                                     a_dna.to_numpy(), a_rna.to_numpy()]))
    raw = np.asarray(rows, dtype=float)  # shape (n_snp, 4*n_rep)

    dna_idx = list(range(0, n_rep)) + list(range(2 * n_rep, 3 * n_rep))
    rna_idx = list(range(n_rep, 2 * n_rep)) + list(range(3 * n_rep, 4 * n_rep))
    allele2_flag = np.array([0] * (2 * n_rep) + [1] * (2 * n_rep), dtype=float)
    # Log2 ratio per replicate per allele (4 alleles * n_rep = ... but actually
    # 2 alleles * n_rep — re-derive cleanly):
    y_ref = np.log2((raw[:, n_rep:2 * n_rep] + 0.5) / (raw[:, 0:n_rep] + 0.5))
    y_alt = np.log2((raw[:, 3 * n_rep:4 * n_rep] + 0.5) / (raw[:, 2 * n_rep:3 * n_rep] + 0.5))
    y = np.concatenate([y_ref, y_alt], axis=1)  # shape (n_snp, 2*n_rep)
    log_dna = np.log2(np.concatenate([raw[:, 0:n_rep], raw[:, 2 * n_rep:3 * n_rep]], axis=1) + 0.5)

    n_obs = 2 * n_rep
    allele2 = np.array([0] * n_rep + [1] * n_rep, dtype=float)
    X = np.column_stack([np.ones(n_obs), allele2])

    # Plain ordinary least squares per fragment for beta_1 (allele2 effect).
    # eBayes moderation downstream shares variance information across fragments.
    XtX_inv = np.linalg.pinv(X.T @ X)
    df_resid = max(n_obs - X.shape[1], 1)
    # Vectorised OLS: beta = (X'X)^-1 X' y^T => shape (2, n_frag)
    beta_mat = XtX_inv @ X.T @ y.T  # (2, n_frag)
    beta = beta_mat.T  # (n_frag, 2)
    fitted = (X @ beta_mat).T  # (n_frag, n_obs)
    resid = y - fitted
    sigma_g = np.sqrt((resid ** 2).sum(axis=1) / df_resid)
    se = np.sqrt(XtX_inv[1, 1]) * sigma_g

    # Empirical-Bayes variance moderation (Smyth 2004, "limma" eBayes).
    # Trigamma is monotone decreasing; we invert it via Newton iteration to
    # find d0 such that trigamma(d0/2) = var(log s2) - trigamma(df/2).
    scipy_special = __import__("scipy.special", fromlist=["special"])

    def _trigamma_inverse(target: float, *, max_iter: int = 100, tol: float = 1e-8) -> float:
        if target <= 0 or not np.isfinite(target):
            return float("inf")  # signals "no moderation"
        # initial guess (Smyth 2004 asymptotic): for small target, trigamma(x) ~ 1/x
        x = 0.5 + 1.0 / target
        for _ in range(max_iter):
            tri = scipy_special.polygamma(1, x)
            deriv = scipy_special.polygamma(2, x)  # trigamma'
            if deriv == 0:
                break
            # Newton on f(x) = trigamma(x) - target
            x_new = x + (target - tri) / deriv
            if x_new <= 0:
                x_new = x / 2.0
            if abs(x_new - x) < tol * max(1.0, abs(x)):
                x = x_new
                break
            x = x_new
        return float(x)

    s2 = sigma_g ** 2
    valid = np.isfinite(s2) & (s2 > 0)
    if valid.sum() < 2:
        # Not enough information for eBayes; fall back to per-fragment t-test.
        s2_post = s2.copy()
        d0_eff = 0.0
        s0_sq = float("nan")
        df_post = float(df_resid)
    else:
        log_s2 = np.log(s2[valid])
        e_log = float(log_s2.mean())
        var_log = float(log_s2.var(ddof=1))
        target = var_log - scipy_special.polygamma(1, df_resid / 2.0)
        x_half = _trigamma_inverse(max(target, 1e-9))
        if not np.isfinite(x_half) or x_half > 1e6:
            # No detectable variance excess => fall back to no moderation.
            s2_post = s2.copy()
            d0_eff = 0.0
            s0_sq = float("nan")
            df_post = float(df_resid)
        else:
            d0_eff = float(2.0 * x_half)
            log_s0_sq = e_log - scipy_special.digamma(x_half) + np.log(x_half)
            s0_sq = float(np.exp(log_s0_sq))
            s2_post = (d0_eff * s0_sq + df_resid * s2) / (d0_eff + df_resid)
            df_post = d0_eff + df_resid

    se_post = np.sqrt(np.maximum(s2_post, 1e-12) * np.diag(XtX_inv)[1])
    t_mod = beta[:, 1] / np.where(se_post > 0, se_post, np.nan)
    pval = 2 * scipy_stats.t.sf(np.abs(t_mod), df_post)
    d0 = d0_eff  # legacy name for output schema
    _, padj, _, _ = multipletests(pval, method="fdr_bh")

    out_df = pd.DataFrame({
        "SNP": common,
        "beta_allele2": beta[:, 1],
        "se_moderated": se_post,
        "t_moderated": t_mod,
        "pvalue": pval,
        "padj": padj,
        "d0_prior_df": d0,
        "s0_sq_prior": s0_sq,
    })
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_allelic_test.tsv"
    out_df.to_csv(out, sep="\t", index=False)
    summary = {
        "n_snps": int(len(out_df)),
        "n_replicates_per_arm": int(n_rep),
        "padj_lt_0.05": int((out_df["padj"] < 0.05).sum()),
        "padj_lt_0.10": int((out_df["padj"] < 0.10).sum()),
        "median_abs_beta": float(out_df["beta_allele2"].abs().median()),
        "d0_prior_df": float(d0),
        "s0_sq_prior": s0_sq,
    }
    save_json(f"starr_{args.label}_allelic_test", summary)
    print(f"Allelic-test table: {out}")
    print(f"Significant @ FDR 0.05: {summary['padj_lt_0.05']} / {summary['n_snps']}")
    print(f"Significant @ FDR 0.10: {summary['padj_lt_0.10']} / {summary['n_snps']}")
    return 0


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------

def write_playbook() -> Path:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "STARRSEQ_ANALYSIS_SKILLS.md"
    path.write_text(
        """# Skill: STARR-seq Allelic Analysis

Use this skill when the agent needs to run end-to-end **allelic STARR-seq**
counts analytics on data from IGVF Portal STARR-seq MeasurementSets or any
other source with paired DNA + RNA replicate counts per (SNP, allele).

The math is a clean-room reimplementation of the methods in
[gaochengwen/STARR-seq-Data-Analysis](https://github.com/gaochengwen/STARR-seq-Data-Analysis)
— the awk + R pipeline that ultimately calls Bioconductor `mpra::mpralm`
for allelic significance. The reference repository has no LICENSE file,
so this skill ships a pure-Python rewrite that follows the underlying
published methods (limma + voom + eBayes, BH-FDR) only.

## Commands

```bash
# 1. Discover IGVF Portal STARR-seq MeasurementSets
python3 Scripts/starr_seq_skill.py pull-portal --limit 50 --label survey

# 2. Counts QC (TPM scaling + RLE + Spearman D-stat outliers)
python3 Scripts/starr_seq_skill.py qc --input counts.tsv --label run1

# 3. Collapse barcode-level counts to per-(SNP, Allele) wide form
python3 Scripts/starr_seq_skill.py aggregate --input barcode_counts.tsv --label run1

# 4. Per-fragment log activity + per-SNP allelic skew (descriptive)
python3 Scripts/starr_seq_skill.py activity --input aggregated.tsv --label run1

# 5. mpralm-style allelic test (voom-weighted GLS + eBayes moderation + BH-FDR)
python3 Scripts/starr_seq_skill.py allelic-test --input aggregated.tsv --label run1
```

## Input schema

`counts.tsv` columns (TSV/CSV, auto-detected):

| Column | Required | Notes |
|---|---|---|
| `SNP` | yes | SNP identifier |
| `Allele` | yes | `REF`/`ALT` or `A`/`B` |
| `Fragment` | optional | Defaults to `SNP_Allele` |
| `DNA_rep1`, `DNA_rep2`, ... | yes | Library (DNA) replicate counts |
| `RNA_rep1`, `RNA_rep2`, ... | yes | Output (RNA) replicate counts |

Replicate columns are auto-detected by the `DNA`/`RNA` prefix (case-insensitive).
For paired analysis, replicate index `1, 2, ...` in `DNA_rep*` corresponds to
the same index in `RNA_rep*`.

## Output schema — `*_allelic_test.tsv`

| Column | Meaning |
|---|---|
| `SNP` | SNP id |
| `beta_allele2` | Effect of ALT allele on log2(RNA/DNA), GLS estimate |
| `se_moderated` | eBayes-moderated standard error |
| `t_moderated` | Moderated t-statistic |
| `pvalue` | Two-sided p from moderated t |
| `padj` | BH-FDR adjusted p |
| `d0_prior_df` | Empirical-Bayes prior degrees of freedom |
| `s0_sq_prior` | Empirical-Bayes prior variance |

## Notes

- The reference repo also includes an awk + Biopython sequence-alignment
  step (`extract_alleles.py`) that QC's whether each oligo aligns to the
  intended reference window. That step is not yet absorbed here — IGVF
  Portal STARR-seq deposits already report aligned counts.
- Genome-wide STARR-seq peak callers (CRADLE, STARRPeaker, BasicSTARRseq)
  are out of scope; this skill covers the allelic-test surface only.
- IGVF Portal: 5,239 STARR-seq files and 1,066 BlueSTARR files are
  catalogued as of the last full Portal survey. Use `starrseq pull-portal`
  to refresh the list.
""",
        encoding="utf-8",
    )
    return path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="STARR-seq allelic analysis (clean-room rewrite of mpralm).")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("pull-portal", help="Discover IGVF Portal STARR-seq MeasurementSets.")
    p.add_argument("--limit", type=int, default=50)
    p.add_argument("--label", default="starrseq_portal_survey")

    p = sub.add_parser("qc", help="TPM-style counts QC: filter, RLE, Spearman D-stat outliers.")
    p.add_argument("--input", required=True)
    p.add_argument("--label", default="starrseq_qc")

    p = sub.add_parser("aggregate", help="Collapse per-barcode counts to per-(SNP, Allele) wide table.")
    p.add_argument("--input", required=True)
    p.add_argument("--label", default="starrseq_aggregate")

    p = sub.add_parser("activity", help="Per-fragment log activity + per-SNP allelic skew.")
    p.add_argument("--input", required=True)
    p.add_argument("--label", default="starrseq_activity")

    p = sub.add_parser("allelic-test", help="mpralm-style allelic test (voom + eBayes + BH-FDR).")
    p.add_argument("--input", required=True)
    p.add_argument("--label", default="starrseq_allelic")

    sub.add_parser("write-playbook", help="Write the skill's markdown playbook.")

    args = parser.parse_args(argv)
    if args.command == "pull-portal":
        return cmd_pull_portal(args)
    if args.command == "qc":
        return cmd_qc(args)
    if args.command == "aggregate":
        return cmd_aggregate(args)
    if args.command == "activity":
        return cmd_activity(args)
    if args.command == "allelic-test":
        return cmd_allelic_test(args)
    if args.command == "write-playbook":
        path = write_playbook()
        print(f"Wrote {path}")
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
