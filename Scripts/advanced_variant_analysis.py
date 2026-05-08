#!/usr/bin/env python3
"""Advanced variant annotation + integrated functional analysis skill.

This skill generalizes the IGVF-Annotation-Analysis R pipeline (Variant
Annotation -> Predicted Functional composite -> Logistic / Miami / overlap
plots -> research report) so it can run on any user-supplied variant CSV
against IGVF Catalog and ENCODE public evidence.

Pipeline:

  1. Read a user variant list (rsid OR chr/pos/ref/alt OR SPDI).
  2. Pull IGVF Catalog evidence for each variant (aPC scores / MACIE-like
     scores via the summary endpoint, QTL genes, phenotypes, regulatory
     elements, predictions). Backed by `Data/Cache/AdvancedVariantAnalysis`.
  3. Overlay ENCODE cCRE classes (PLS / pELS / dELS) from a public BED file.
  4. Join gene-proximity windows from the Catalog gene endpoint.
  5. Compute a `Predicted_Functional` composite from any subset of:
        aPC > 20, MACIE > 20, ClinVar pathogenic, QTL evidence,
        regulatory-element overlap, prediction-set overlap.
  6. Optionally merge user-provided experimental columns
     (CRISPRi BEAN / MPRA log2FC / GWAS p-value, etc.) by variant key.
  7. Produce:
        - `<label>_annotated.csv`           unified per-variant matrix
        - `<label>_summary_stats.csv`       per-feature counts and rates
        - `<label>_logistic_model.json`     logistic model fits
        - `<label>_logistic_model.svg`      coefficient bar plot
        - `<label>_volcano.svg`             experimental volcano (if p+effect)
        - `<label>_miami_<gene>.svg`        per-gene Miami plot (if --gene)
        - `<label>_evidence_overlap.svg`    stacked bar of feature overlap
        - `<label>_analysis_report.md`      research interpretation report

Design notes:
  - No confidential reference files are required. ENCODE cCRE BED is
    auto-downloaded once into `Data/Cache/AdvancedVariantAnnotations/cCRE/`.
  - The skill is environment-driven; `IGVF_CATALOG_API_BASE`, etc. apply.
  - Logistic regression uses scipy/numpy only (no scikit-learn dependency).
  - Plots use matplotlib SVG.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
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
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

# --- Project paths -----------------------------------------------------------
ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
CACHE_DIR = DATA_DIR / "Cache" / "AdvancedVariantAnnotations"
CCRE_CACHE_DIR = CACHE_DIR / "cCRE"
LOG_DIR = ROOT / "Docs" / "Logs"
DOCS_DIR = ROOT / "Docs" / "AdvancedVariantAnalysis"
PLOTS_DIR = DOCS_DIR / "Plots"

# --- Endpoints ---------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint

CATALOG_API_BASE = _resolve_endpoint("catalog_api", "IGVF_CATALOG_API_BASE")
ENCODE_BASE = _resolve_endpoint("encode", "ENCODE_BASE")
# Public cCRE registry BED. Override with env var if needed.
ENCODE_CCRE_URL = os.environ.get(
    "ENCODE_CCRE_BED_URL",
    f"{_resolve_endpoint('wenglab_dl')}/Registry-V3/GRCh38-cCREs.bed",
)

DEFAULT_THRESHOLD = 20.0          # aPC / MACIE PHRED-style threshold
PFUNC_FLAGS = (
    "PFunc_aPC",
    "PFunc_MACIE",
    "PFunc_ClinVar",
    "PFunc_QTL",
    "PFunc_RegElement",
    "PFunc_Prediction",
)


# --- Logging ----------------------------------------------------------------
def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"advanced_variant_analysis_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Advanced variant analysis log: %s", log_path)
    return log_path


# --- HTTP helpers -----------------------------------------------------------
def _fetch(url: str, timeout: int = 60) -> bytes:
    req = urllib.request.Request(
        url,
        headers={"Accept": "application/json,application/octet-stream,*/*", "User-Agent": "IGVFagent-AVA/0.1"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return response.read()


def _json_get(url: str, timeout: int = 60) -> Any:
    raw = _fetch(url, timeout=timeout)
    try:
        return json.loads(raw.decode("utf-8"))
    except json.JSONDecodeError:
        return None


def _cache_key(payload: str) -> str:
    return hashlib.sha256(payload.encode()).hexdigest()


def _cache_get(category: str, key: str) -> Any:
    p = CACHE_DIR / category / f"{key}.json"
    if p.exists():
        try:
            return json.loads(p.read_text("utf-8"))
        except Exception:  # noqa: BLE001
            return None
    return None


def _cache_put(category: str, key: str, value: Any) -> None:
    p = CACHE_DIR / category / f"{key}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(value), encoding="utf-8")


# --- Variant identity -------------------------------------------------------
@dataclass
class VariantKey:
    raw: dict[str, str]
    rsid: str | None = None
    chrom: str | None = None
    pos: int | None = None
    ref: str | None = None
    alt: str | None = None
    spdi: str | None = None
    hgvs: str | None = None
    label: str = ""

    @property
    def variant_id(self) -> str:
        """Canonical id like '19-44908822-C-T' when alleles are known, else rsid."""
        if self.chrom and self.pos and self.ref and self.alt:
            return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"
        return self.rsid or self.spdi or self.hgvs or self.label or ""


def _norm_chrom(value: str | None) -> str | None:
    if not value:
        return None
    s = str(value).strip().replace("chr", "")
    if not s:
        return None
    return "X" if s.upper() == "X" else s


def _read_variants(path: Path) -> list[VariantKey]:
    rows: list[VariantKey] = []
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise SystemExit(f"No CSV header in {path}")
        # Build lowercase column index
        lower = {c.lower().strip(): c for c in reader.fieldnames}

        def col(*names: str) -> str | None:
            for n in names:
                if n in lower:
                    return lower[n]
            return None

        rsid_c = col("rsid", "rsids", "dbsnp")
        chr_c = col("chrom", "chr", "chromosome")
        pos_c = col("pos", "position", "start", "bp")
        ref_c = col("ref", "reference", "allele1")
        alt_c = col("alt", "alternate", "allele2")
        spdi_c = col("spdi", "spdi_id")
        hgvs_c = col("hgvs", "hgvs_id")
        label_c = col("label", "name", "id", "variant", "variant_id")

        for r in reader:
            rsid = r.get(rsid_c) if rsid_c else None
            chrom = _norm_chrom(r.get(chr_c)) if chr_c else None
            pos = r.get(pos_c) if pos_c else None
            try:
                pos_i = int(pos) if pos and str(pos).strip() else None
            except (TypeError, ValueError):
                pos_i = None
            ref = (r.get(ref_c) or "").strip().upper() or None if ref_c else None
            alt = (r.get(alt_c) or "").strip().upper() or None if alt_c else None
            spdi = (r.get(spdi_c) or "").strip() or None if spdi_c else None
            hgvs = (r.get(hgvs_c) or "").strip() or None if hgvs_c else None
            label = (r.get(label_c) or "").strip() or None if label_c else None

            rows.append(
                VariantKey(
                    raw=dict(r),
                    rsid=rsid.strip() if rsid else None,
                    chrom=chrom,
                    pos=pos_i,
                    ref=ref,
                    alt=alt,
                    spdi=spdi,
                    hgvs=hgvs,
                    label=label or rsid or spdi or "",
                )
            )
    return rows


# --- Catalog evidence -------------------------------------------------------
def _identifier_for_catalog(v: VariantKey) -> str | None:
    if v.rsid:
        return v.rsid
    if v.spdi:
        return v.spdi
    if v.hgvs:
        return v.hgvs
    if v.chrom and v.pos and v.ref and v.alt:
        return f"chr{v.chrom}:{v.pos}-{v.ref}-{v.alt}"
    return None


def fetch_catalog_evidence(v: VariantKey, *, sleep: float = 0.05) -> dict[str, Any]:
    """Return a flattened dict of Catalog-derived evidence for one variant."""
    ident = _identifier_for_catalog(v)
    if not ident:
        return {"_query_status": "missing_variant_identifier"}

    out: dict[str, Any] = {"_query_status": "ok", "_identifier": ident}
    cache_key = _cache_key(ident)

    endpoints = {
        "summary": f"{CATALOG_API_BASE}/api/variants/summary?variant_id={urllib.parse.quote(ident)}",
        "qtl": f"{CATALOG_API_BASE}/api/variants/genes?variant_id={urllib.parse.quote(ident)}&limit=50",
        "phenotypes": f"{CATALOG_API_BASE}/api/variants/phenotypes?variant_id={urllib.parse.quote(ident)}&limit=50",
        "elements": f"{CATALOG_API_BASE}/api/variants/genomic-elements?variant_id={urllib.parse.quote(ident)}&limit=50",
        "predictions": f"{CATALOG_API_BASE}/api/variants/predictions?variant_id={urllib.parse.quote(ident)}&limit=50",
    }

    for name, url in endpoints.items():
        cached = _cache_get(f"catalog_{name}", cache_key)
        if cached is not None:
            payload = cached
        else:
            try:
                payload = _json_get(url)
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
                logging.warning("catalog %s failed for %s: %s", name, ident, exc)
                payload = {"_error": str(exc)}
            _cache_put(f"catalog_{name}", cache_key, payload)
            time.sleep(sleep)

        out[f"raw_{name}"] = payload

    # Flatten summary scores (if present)
    summary = out.get("raw_summary") or {}
    if isinstance(summary, list) and summary:
        summary = summary[0]
    if isinstance(summary, dict):
        for key in (
            "cadd_phred", "cadd_raw",
            "apc_conservation", "apc_protein_function",
            "apc_epigenetics_active", "apc_epigenetics_repressed",
            "apc_epigenetics_transcription", "apc_transcription_factor",
            "apc_local_nucleotide_diversity",
            "macie_anyclass", "macie_conserved", "macie_protein", "macie_regulatory",
            "clinvar_significance", "rsid", "chr", "position",
        ):
            if key in summary and summary[key] is not None:
                out[key] = summary[key]

    # Counts of related items
    out["qtl_count"] = len(out.get("raw_qtl") or [])
    out["phenotype_count"] = len(out.get("raw_phenotypes") or [])
    out["element_count"] = len(out.get("raw_elements") or [])
    out["prediction_count"] = len(out.get("raw_predictions") or [])
    return out


# --- ENCODE cCRE overlay ----------------------------------------------------
def ensure_ccre_bed() -> Path:
    """Download ENCODE GRCh38 cCRE BED once into the cache."""
    CCRE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    bed_path = CCRE_CACHE_DIR / "GRCh38-cCREs.bed"
    if bed_path.exists() and bed_path.stat().st_size > 0:
        return bed_path

    logging.info("Fetching ENCODE cCRE BED from %s", ENCODE_CCRE_URL)
    try:
        raw = _fetch(ENCODE_CCRE_URL, timeout=120)
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
        logging.warning("ENCODE cCRE download failed (%s); writing empty BED.", exc)
        bed_path.write_text("", encoding="utf-8")
        return bed_path

    if ENCODE_CCRE_URL.endswith(".gz"):
        raw = gzip.decompress(raw)
    bed_path.write_bytes(raw)
    logging.info("Wrote %s (%d bytes)", bed_path, bed_path.stat().st_size)
    return bed_path


def _index_ccre(bed_path: Path) -> dict[str, list[tuple[int, int, str]]]:
    """Build a per-chrom interval list of (start, end, class)."""
    idx: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    if not bed_path.exists() or bed_path.stat().st_size == 0:
        return idx
    with bed_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6 or not parts[1].isdigit():
                continue
            chrom = parts[0].replace("chr", "")
            classification = parts[5]
            cclass: str | None = None
            for tag in ("PLS", "pELS", "dELS"):
                if tag in classification:
                    cclass = tag
                    break
            if cclass is None:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            idx[chrom].append((start, end, cclass))
    for chrom in list(idx.keys()):
        idx[chrom].sort()
    return idx


def overlay_ccre(v: VariantKey, ccre_idx: dict[str, list[tuple[int, int, str]]]) -> dict[str, Any]:
    """Return the cCRE classification for a variant if any."""
    if not (v.chrom and v.pos):
        return {"encode_ccre_class": None}
    items = ccre_idx.get(v.chrom, [])
    # Linear scan; lists are small per chrom for typical user lists.
    for start, end, cclass in items:
        if start <= v.pos <= end:
            return {"encode_ccre_class": cclass}
    return {"encode_ccre_class": None}


# --- Predicted functional composite -----------------------------------------
def compute_predicted_functional(row: dict[str, Any], threshold: float) -> dict[str, Any]:
    """Compute Predicted_Functional flags + composite from one annotated row."""
    def _gt(key: str) -> bool:
        try:
            return float(row[key]) >= threshold
        except (TypeError, ValueError, KeyError):
            return False

    pfunc_apc = any(
        _gt(k)
        for k in (
            "apc_conservation", "apc_protein_function",
            "apc_epigenetics_active", "apc_epigenetics_repressed",
            "apc_epigenetics_transcription", "apc_transcription_factor",
        )
    )
    pfunc_macie = any(
        _gt(k)
        for k in ("macie_anyclass", "macie_conserved", "macie_protein", "macie_regulatory")
    )
    clin = (row.get("clinvar_significance") or "").lower()
    pfunc_clinvar = any(t in clin for t in ("pathogenic", "likely_pathogenic", "drug_response"))
    pfunc_qtl = (row.get("qtl_count") or 0) > 0
    pfunc_reg = bool(row.get("encode_ccre_class")) or (row.get("element_count") or 0) > 0
    pfunc_pred = (row.get("prediction_count") or 0) > 0

    flags = {
        "PFunc_aPC": pfunc_apc,
        "PFunc_MACIE": pfunc_macie,
        "PFunc_ClinVar": pfunc_clinvar,
        "PFunc_QTL": pfunc_qtl,
        "PFunc_RegElement": pfunc_reg,
        "PFunc_Prediction": pfunc_pred,
    }
    flags["PFunc_Any"] = any(flags.values())
    flags["PFunc_Sum"] = sum(1 for k in PFUNC_FLAGS if flags[k])
    return flags


# --- Logistic regression (pure-Python, IRLS) --------------------------------
def _logistic_irls(X: list[list[float]], y: list[float], max_iter: int = 50, tol: float = 1e-6) -> dict[str, Any]:
    """Plain IRLS logistic regression. Returns coefs, std errors, p-values."""
    try:
        import numpy as np
        from scipy.stats import norm
    except ImportError as exc:
        raise SystemExit("numpy + scipy are required for logistic models") from exc

    Xa = np.asarray(X, dtype=float)
    ya = np.asarray(y, dtype=float)
    n, p = Xa.shape
    beta = np.zeros(p)

    for _ in range(max_iter):
        eta = Xa @ beta
        mu = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
        W = mu * (1 - mu)
        z = eta + (ya - mu) / np.where(W < 1e-8, 1e-8, W)
        XtWX = Xa.T * W @ Xa
        XtWz = (Xa.T * W) @ z
        try:
            beta_new = np.linalg.solve(XtWX + 1e-8 * np.eye(p), XtWz)
        except np.linalg.LinAlgError:
            break
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break
        beta = beta_new

    # Std errors
    eta = Xa @ beta
    mu = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
    W = mu * (1 - mu)
    cov = np.linalg.pinv(Xa.T * W @ Xa)
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    z_scores = beta / np.where(se > 0, se, 1.0)
    pvals = 2 * (1 - norm.cdf(np.abs(z_scores)))
    return {
        "beta": beta.tolist(),
        "se": se.tolist(),
        "z": z_scores.tolist(),
        "pvalue": pvals.tolist(),
    }


def fit_logistic_models(
    rows: list[dict[str, Any]],
    *,
    outcome_col: str,
    predictor_cols: list[str],
) -> dict[str, Any]:
    """Fit binary, marginal, and joint logistic models and return summaries."""
    try:
        import numpy as np
    except ImportError as exc:
        raise SystemExit("numpy is required for logistic models") from exc

    def _to_float(v: Any) -> float:
        if isinstance(v, bool):
            return 1.0 if v else 0.0
        try:
            return float(v)
        except (TypeError, ValueError):
            return float("nan")

    y_full = [1.0 if (r.get(outcome_col) in (True, 1, "1", "TRUE", "True", "true")) else 0.0 for r in rows]
    if sum(y_full) == 0 or sum(y_full) == len(y_full):
        return {"error": f"outcome '{outcome_col}' has no variation; skipping logistic models."}

    results: dict[str, Any] = {"outcome": outcome_col, "n_total": len(rows), "n_positive": int(sum(y_full)), "marginal": {}, "joint": {}}

    # Marginal: each predictor alone (with intercept)
    for pred in predictor_cols:
        x_vals = [_to_float(r.get(pred)) for r in rows]
        keep = [i for i, x in enumerate(x_vals) if not (math.isnan(x))]
        if len(keep) < 5:
            results["marginal"][pred] = {"error": "insufficient non-NA rows"}
            continue
        X = [[1.0, x_vals[i]] for i in keep]
        y = [y_full[i] for i in keep]
        try:
            fit = _logistic_irls(X, y)
            results["marginal"][pred] = {
                "n": len(y),
                "intercept": fit["beta"][0],
                "beta": fit["beta"][1],
                "se": fit["se"][1],
                "pvalue": fit["pvalue"][1],
                "or": math.exp(fit["beta"][1]) if abs(fit["beta"][1]) < 30 else None,
            }
        except Exception as exc:  # noqa: BLE001
            results["marginal"][pred] = {"error": str(exc)}

    # Joint: all predictors together, complete-case rows
    matrix: list[list[float]] = []
    y_keep: list[float] = []
    for r, yv in zip(rows, y_full):
        xs = [_to_float(r.get(p)) for p in predictor_cols]
        if any(math.isnan(x) for x in xs):
            continue
        matrix.append([1.0] + xs)
        y_keep.append(yv)
    if len(matrix) >= len(predictor_cols) + 5 and 0 < sum(y_keep) < len(y_keep):
        try:
            fit = _logistic_irls(matrix, y_keep)
            results["joint"] = {
                "n": len(y_keep),
                "predictors": predictor_cols,
                "intercept": fit["beta"][0],
                "beta": fit["beta"][1:],
                "se": fit["se"][1:],
                "pvalue": fit["pvalue"][1:],
                "or": [math.exp(b) if abs(b) < 30 else None for b in fit["beta"][1:]],
            }
        except Exception as exc:  # noqa: BLE001
            results["joint"] = {"error": str(exc)}
    else:
        results["joint"] = {"error": "insufficient complete-case rows for joint model"}
    return results


# --- Plotting ---------------------------------------------------------------
def _import_plt():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError as exc:
        raise SystemExit("matplotlib is required for plots") from exc


def plot_logistic_coefficients(model: dict[str, Any], out_path: Path) -> None:
    plt = _import_plt()
    joint = model.get("joint") or {}
    if "beta" not in joint:
        return
    preds = joint["predictors"]
    beta = joint["beta"]
    se = joint["se"]
    pvals = joint["pvalue"]
    fig, ax = plt.subplots(figsize=(6.5, max(2.5, 0.45 * len(preds) + 1)))
    ypos = list(range(len(preds)))
    ax.barh(ypos, beta, xerr=[1.96 * s for s in se], color=["#1f77b4" if p < 0.05 else "#bbbbbb" for p in pvals])
    ax.set_yticks(ypos)
    ax.set_yticklabels(preds)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel(f"Joint logistic β  (outcome: {model.get('outcome')})")
    ax.set_title(f"Predictors of {model.get('outcome')} (n={joint.get('n')})")
    fig.tight_layout()
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def plot_volcano(rows: list[dict[str, Any]], effect_col: str, p_col: str, out_path: Path, label: str) -> None:
    plt = _import_plt()
    xs, ys, colors = [], [], []
    for r in rows:
        try:
            x = float(r.get(effect_col))
            p = float(r.get(p_col))
        except (TypeError, ValueError):
            continue
        if p <= 0 or math.isnan(x) or math.isnan(p):
            continue
        xs.append(x)
        ys.append(-math.log10(p))
        colors.append("#d62728" if r.get("PFunc_Any") else "#7f7f7f")
    if not xs:
        return
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.scatter(xs, ys, c=colors, s=18, alpha=0.7)
    ax.axhline(-math.log10(0.05), color="grey", linestyle="--", lw=0.8)
    ax.set_xlabel(effect_col)
    ax.set_ylabel(f"-log10({p_col})")
    ax.set_title(f"Volcano: {label}\nred = Predicted Functional (any evidence)")
    fig.tight_layout()
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def plot_evidence_overlap(rows: list[dict[str, Any]], out_path: Path, label: str) -> None:
    plt = _import_plt()
    counts = {flag: sum(1 for r in rows if r.get(flag)) for flag in PFUNC_FLAGS}
    counts["Any"] = sum(1 for r in rows if r.get("PFunc_Any"))
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(list(counts.keys()), list(counts.values()), color="#1f77b4")
    ax.set_ylabel("variants")
    ax.set_title(f"Evidence overlap (n={len(rows)}): {label}")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def plot_miami_per_gene(
    rows: list[dict[str, Any]],
    *,
    gene: str,
    p_col_top: str,
    p_col_bottom: str | None,
    out_path: Path,
    nearest_col: str = "_nearest_gene",
) -> None:
    plt = _import_plt()
    sub = [r for r in rows if (r.get(nearest_col) or "").upper().startswith(gene.upper())]
    if not sub:
        return
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 5), sharex=True)

    def _scatter(ax, p_col, invert):
        xs, ys, colors = [], [], []
        for r in sub:
            try:
                p = float(r.get(p_col))
                if p <= 0:
                    continue
            except (TypeError, ValueError):
                continue
            try:
                pos = float(r.get("position") or r.get("pos") or 0)
            except (TypeError, ValueError):
                pos = 0
            xs.append(pos)
            ys.append(-math.log10(p))
            colors.append("#d62728" if r.get("PFunc_Any") else "#7f7f7f")
        if not xs:
            return
        if invert:
            ys = [-y for y in ys]
        ax.scatter(xs, ys, c=colors, s=20, alpha=0.8)
        ax.axhline(-math.log10(0.05) * (-1 if invert else 1), color="grey", lw=0.6, ls="--")
        ax.set_ylabel(("- " if invert else "") + f"-log10({p_col})")

    _scatter(ax1, p_col_top, invert=False)
    if p_col_bottom and any(r.get(p_col_bottom) is not None for r in sub):
        _scatter(ax2, p_col_bottom, invert=True)
    else:
        ax2.set_visible(False)
    ax1.set_title(f"Miami plot for {gene} (n={len(sub)})")
    ax2.set_xlabel("Position (bp)")
    fig.tight_layout()
    fig.savefig(out_path, format="svg")
    plt.close(fig)


# --- Annotation pipeline ----------------------------------------------------
def annotate_variants(
    variants: list[VariantKey],
    *,
    sleep: float = 0.05,
) -> tuple[list[dict[str, Any]], list[str]]:
    bed_path = ensure_ccre_bed()
    ccre_idx = _index_ccre(bed_path)

    all_rows: list[dict[str, Any]] = []
    extra_keys: set[str] = set()

    for i, v in enumerate(variants):
        if i % 25 == 0:
            logging.info("Annotating variant %d / %d", i + 1, len(variants))
        catalog = fetch_catalog_evidence(v, sleep=sleep)
        ccre = overlay_ccre(v, ccre_idx)

        # Resolve coords from catalog summary if missing on input
        if (not v.chrom) and isinstance(catalog.get("chr"), (str, int)):
            v.chrom = _norm_chrom(str(catalog["chr"]))
            ccre = overlay_ccre(v, ccre_idx)
        if (not v.pos) and isinstance(catalog.get("position"), (str, int)):
            try:
                v.pos = int(catalog["position"])
            except (TypeError, ValueError):
                pass
        if (not v.rsid) and catalog.get("rsid"):
            v.rsid = catalog["rsid"]

        flat: dict[str, Any] = dict(v.raw)  # preserve input columns
        flat["variant_id"] = v.variant_id
        flat["rsid_resolved"] = v.rsid
        flat["chrom_resolved"] = v.chrom
        flat["position"] = v.pos
        flat["ref"] = v.ref
        flat["alt"] = v.alt
        for k, val in catalog.items():
            if k.startswith("raw_") or k.startswith("_"):
                continue
            flat[k] = val
        flat.update(ccre)

        # Nearest-gene heuristic from QTL/predictions raw payloads (best-effort)
        nearest = None
        for raw_key in ("raw_qtl", "raw_predictions", "raw_phenotypes"):
            payload = catalog.get(raw_key) or []
            if isinstance(payload, list):
                for entry in payload:
                    if isinstance(entry, dict):
                        for cand in ("gene_name", "hgnc_symbol", "symbol", "nearest_gene"):
                            if entry.get(cand):
                                nearest = entry.get(cand)
                                break
                    if nearest:
                        break
            if nearest:
                break
        flat["_nearest_gene"] = nearest

        flat.update(compute_predicted_functional(flat, DEFAULT_THRESHOLD))
        all_rows.append(flat)
        for k in flat:
            extra_keys.add(k)

    # Stable column order: input cols, then enrichment cols
    base_cols = list(variants[0].raw.keys()) if variants else []
    enrich_cols = [
        "variant_id", "rsid_resolved", "chrom_resolved", "position", "ref", "alt",
        "cadd_phred",
        "apc_conservation", "apc_protein_function",
        "apc_epigenetics_active", "apc_epigenetics_repressed",
        "apc_epigenetics_transcription", "apc_transcription_factor",
        "apc_local_nucleotide_diversity",
        "macie_anyclass", "macie_conserved", "macie_protein", "macie_regulatory",
        "clinvar_significance",
        "qtl_count", "phenotype_count", "element_count", "prediction_count",
        "encode_ccre_class", "_nearest_gene",
        *PFUNC_FLAGS, "PFunc_Any", "PFunc_Sum",
    ]
    seen = set()
    cols: list[str] = []
    for c in base_cols + enrich_cols:
        if c not in seen and c in extra_keys:
            cols.append(c)
            seen.add(c)
    # Append any other accidentally produced columns
    for k in sorted(extra_keys):
        if k not in seen:
            cols.append(k)
            seen.add(k)
    return all_rows, cols


# --- Optional experimental table merge --------------------------------------
def merge_experimental(
    rows: list[dict[str, Any]],
    experimental_csv: Path,
    join_col: str | None,
) -> tuple[list[dict[str, Any]], list[str]]:
    if not experimental_csv.exists():
        raise SystemExit(f"experimental file not found: {experimental_csv}")
    with experimental_csv.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        exp_fields = list(reader.fieldnames or [])
        # Pick join column
        if not join_col:
            for cand in ("variant_id", "variant", "rsid", "rsID", "rs_id"):
                if cand in exp_fields:
                    join_col = cand
                    break
        if not join_col:
            raise SystemExit("--experimental requires a join column; pass --join-col explicitly")
        index: dict[str, dict[str, Any]] = {}
        for row in reader:
            key = (row.get(join_col) or "").strip()
            if key:
                index[key] = row

    added: list[str] = []
    for row in rows:
        key_candidates = [
            row.get("variant_id"), row.get("rsid_resolved"),
            row.get("rsid"), row.get("rsID"), row.get("variant"),
        ]
        match: dict[str, Any] | None = None
        for k in key_candidates:
            if k and k in index:
                match = index[k]
                break
        if match:
            for f in exp_fields:
                if f == join_col:
                    continue
                row[f] = match.get(f)
                if f not in added:
                    added.append(f)
    return rows, added


# --- Reporting --------------------------------------------------------------
def _csv_write(path: Path, rows: list[dict[str, Any]], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_summary_stats(rows: list[dict[str, Any]], out_path: Path) -> dict[str, int]:
    counts = {flag: sum(1 for r in rows if r.get(flag)) for flag in (*PFUNC_FLAGS, "PFunc_Any")}
    counts["n_total"] = len(rows)
    counts["n_with_rsid"] = sum(1 for r in rows if r.get("rsid_resolved"))
    counts["n_with_position"] = sum(1 for r in rows if r.get("position"))
    counts["n_in_ccre"] = sum(1 for r in rows if r.get("encode_ccre_class"))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["feature", "count"])
        for k, v in counts.items():
            w.writerow([k, v])
    return counts


def write_report(
    out_path: Path,
    *,
    label: str,
    input_path: Path,
    annotated_csv: Path,
    summary_path: Path,
    counts: dict[str, int],
    model: dict[str, Any] | None,
    plots: list[Path],
    extra_experimental_cols: list[str],
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"# Advanced Variant Analysis Report — {label}",
        "",
        f"- Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        f"- Input variants: `{input_path}`",
        f"- Annotated matrix: `{annotated_csv}`",
        f"- Summary stats: `{summary_path}`",
        "",
        "## Counts",
        "",
        "| Feature | Count |",
        "|---|---|",
    ]
    lines.extend(f"| {k} | {v} |" for k, v in counts.items())
    lines.append("")

    if model:
        lines.append("## Logistic regression")
        lines.append("")
        if "error" in model and "outcome" not in model:
            lines.append(f"- Skipped: {model['error']}")
            lines.append("")
        else:
            lines.append(f"- Outcome: `{model.get('outcome')}`")
            lines.append(f"- Total variants: {model.get('n_total')}")
            lines.append(f"- Positives: {model.get('n_positive')}")
        joint = model.get("joint") or {}
        if "beta" in joint:
            lines.append("")
            lines.append("### Joint model")
            lines.append("")
            lines.append("| Predictor | β | SE | OR | p |")
            lines.append("|---|---|---|---|---|")
            for pred, b, se, p_, or_ in zip(joint["predictors"], joint["beta"], joint["se"], joint["pvalue"], joint["or"]):
                lines.append(f"| {pred} | {b:.3f} | {se:.3f} | {or_:.2f} | {p_:.2e} |" if or_ is not None else f"| {pred} | {b:.3f} | {se:.3f} | NA | {p_:.2e} |")
        elif joint.get("error"):
            lines.append(f"- Joint model: {joint['error']}")

        lines.append("")
        lines.append("### Marginal models")
        lines.append("")
        lines.append("| Predictor | β | SE | OR | p |")
        lines.append("|---|---|---|---|---|")
        for pred, info in (model.get("marginal") or {}).items():
            if "error" in info:
                lines.append(f"| {pred} | — | — | — | {info['error']} |")
            else:
                lines.append(f"| {pred} | {info['beta']:.3f} | {info['se']:.3f} | {info['or']:.2f} | {info['pvalue']:.2e} |"
                             if info.get("or") is not None
                             else f"| {pred} | {info['beta']:.3f} | {info['se']:.3f} | NA | {info['pvalue']:.2e} |")
        lines.append("")

    if extra_experimental_cols:
        lines.append("## Experimental columns merged")
        lines.append("")
        lines.append(", ".join(f"`{c}`" for c in extra_experimental_cols))
        lines.append("")

    if plots:
        lines.append("## Plots")
        lines.append("")
        for p in plots:
            lines.append(f"- `{p}`")
        lines.append("")

    lines.append("## Interpretation guidance")
    lines.append("")
    lines.append(
        "- A variant flagged `PFunc_Any` carries at least one line of computational evidence "
        "(aPC > %.0f, MACIE > %.0f, ClinVar pathogenic, QTL hit, regulatory-element overlap, "
        "or prediction-set overlap)." % (DEFAULT_THRESHOLD, DEFAULT_THRESHOLD)
    )
    lines.append(
        "- The logistic models report whether each functional axis predicts your supplied "
        "experimental outcome (e.g. CRISPRi BEAN p<.05, MPRA significance). Positive β with "
        "p < 0.05 in the joint model indicates an axis that is informative beyond the others."
    )
    lines.append(
        "- The Miami plot shows experimental significance (top, +log10p) and an optional second "
        "track (bottom, −log10p). Red dots are `PFunc_Any` = TRUE."
    )
    lines.append(
        "- Evidence overlap counts let you rank candidate causal variants by the number of "
        "independent functional axes that flag them."
    )
    out_path.write_text("\n".join(lines), encoding="utf-8")


# --- CLI --------------------------------------------------------------------
def cmd_run(args: argparse.Namespace) -> int:
    setup_logging()
    input_path = Path(args.input).resolve()
    if not input_path.exists():
        raise SystemExit(f"input not found: {input_path}")
    label = args.label or input_path.stem

    variants = _read_variants(input_path)
    if args.max_rows:
        variants = variants[: args.max_rows]
    logging.info("Loaded %d variants from %s", len(variants), input_path)

    rows, cols = annotate_variants(variants, sleep=args.sleep)

    extra_exp_cols: list[str] = []
    if args.experimental:
        rows, extra_exp_cols = merge_experimental(rows, Path(args.experimental).resolve(), args.join_col)
        for c in extra_exp_cols:
            if c not in cols:
                cols.append(c)

    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    out_csv = DOCS_DIR / f"{stamp}_{label}_annotated.csv"
    out_summary = DOCS_DIR / f"{stamp}_{label}_summary_stats.csv"
    out_model_json = DOCS_DIR / f"{stamp}_{label}_logistic_model.json"
    out_model_svg = PLOTS_DIR / f"{stamp}_{label}_logistic_model.svg"
    out_overlap_svg = PLOTS_DIR / f"{stamp}_{label}_evidence_overlap.svg"
    out_volcano_svg = PLOTS_DIR / f"{stamp}_{label}_volcano.svg"
    out_report = DOCS_DIR / f"{stamp}_{label}_analysis_report.md"

    _csv_write(out_csv, rows, cols)
    counts = write_summary_stats(rows, out_summary)
    plot_evidence_overlap(rows, out_overlap_svg, label)

    plots: list[Path] = [out_overlap_svg]
    model: dict[str, Any] | None = None
    if args.outcome:
        predictors = (
            args.predictors.split(",")
            if args.predictors
            else list(PFUNC_FLAGS)
        )
        predictors = [p.strip() for p in predictors if p.strip()]
        model = fit_logistic_models(rows, outcome_col=args.outcome, predictor_cols=predictors)
        out_model_json.write_text(json.dumps(model, indent=2, default=str), encoding="utf-8")
        plot_logistic_coefficients(model, out_model_svg)
        plots.append(out_model_svg)

        if args.effect_col and args.outcome:
            # Volcano: effect vs outcome's underlying p-value column if it looks like one
            if args.outcome.lower().endswith(("_p", "_pvalue", "_p_value", "_pval")):
                plot_volcano(rows, args.effect_col, args.outcome, out_volcano_svg, label)
                plots.append(out_volcano_svg)

    if args.gene_list:
        for gene in [g.strip() for g in args.gene_list.split(",") if g.strip()]:
            mp = PLOTS_DIR / f"{stamp}_{label}_miami_{gene}.svg"
            top_p = args.outcome if args.outcome else "qtl_count"
            bot_p = args.bottom_p
            plot_miami_per_gene(rows, gene=gene, p_col_top=top_p, p_col_bottom=bot_p, out_path=mp)
            if mp.exists():
                plots.append(mp)

    write_report(
        out_report,
        label=label,
        input_path=input_path,
        annotated_csv=out_csv,
        summary_path=out_summary,
        counts=counts,
        model=model,
        plots=plots,
        extra_experimental_cols=extra_exp_cols,
    )

    print(f"Annotated CSV : {out_csv}")
    print(f"Summary stats : {out_summary}")
    if model is not None:
        print(f"Logistic model: {out_model_json}")
        print(f"Coefficient   : {out_model_svg}")
    print(f"Overlap plot  : {out_overlap_svg}")
    print(f"Report        : {out_report}")
    return 0


def cmd_write_playbook(args: argparse.Namespace) -> int:
    skills_dir = ROOT / "Docs" / "Skills"
    skills_dir.mkdir(parents=True, exist_ok=True)
    path = skills_dir / "ADVANCED_VARIANT_ANALYSIS_SKILLS.md"
    path.write_text(
        """# Advanced Variant Analysis

This skill produces the full IGVF-style integrated variant analysis from a
generic variant CSV: per-variant feature matrix, predicted-functional
composite, logistic models against any user-supplied experimental outcome,
volcano / Miami / overlap plots, and a research-grade markdown report.

## Inputs

- `--input <csv>`: variant list with any of `rsid`, `chrom/pos/ref/alt`,
  `spdi`, `hgvs`. Extra user columns are preserved.
- `--experimental <csv>` (optional): user effect / p-value table; joined on
  `variant_id` or `rsid` automatically (override with `--join-col`).
- `--outcome <col>`: column treated as the binary outcome for logistic
  models (e.g. `BEAN_pval_lt_.05`, `mpra_sig`).
- `--predictors <comma>`: defaults to all `PFunc_*` flags.
- `--gene-list <comma>`: per-gene Miami plots near these symbols.
- `--effect-col <col>`: continuous effect for volcano (paired with `--outcome`).
- `--bottom-p <col>`: optional second p-value track for the Miami plot.

## Outputs (under `Docs/AdvancedVariantAnalysis/`)

- `<stamp>_<label>_annotated.csv`        unified per-variant matrix
- `<stamp>_<label>_summary_stats.csv`    feature counts
- `<stamp>_<label>_logistic_model.json`  marginal + joint logistic fits
- `Plots/<stamp>_<label>_logistic_model.svg`
- `Plots/<stamp>_<label>_evidence_overlap.svg`
- `Plots/<stamp>_<label>_volcano.svg` (when applicable)
- `Plots/<stamp>_<label>_miami_<gene>.svg`
- `<stamp>_<label>_analysis_report.md`   markdown research report

## Predicted-Functional axes

| Flag | Source |
|---|---|
| `PFunc_aPC` | Catalog summary aPC scores (any > 20) |
| `PFunc_MACIE` | Catalog summary MACIE scores (any > 20) |
| `PFunc_ClinVar` | ClinVar pathogenic / likely-pathogenic / drug response |
| `PFunc_QTL` | At least one QTL gene linked in Catalog |
| `PFunc_RegElement` | ENCODE cCRE class match (PLS / pELS / dELS) or Catalog element link |
| `PFunc_Prediction` | At least one prediction-set link in Catalog |

`PFunc_Any` is the OR; `PFunc_Sum` is the count of axes that fire.

## Example

```bash
python3 Scripts/advanced_variant_analysis.py run \\
  --input Data/Input/VariantList/example_variants.csv \\
  --label example_locus_v1
```

With user experimental data and a custom outcome:

```bash
python3 Scripts/advanced_variant_analysis.py run \\
  --input Data/Input/VariantList/my_variants.csv \\
  --experimental Data/Input/Experimental/my_crispri.csv \\
  --join-col variant_id \\
  --outcome BEAN_pval_lt_.05 \\
  --gene-list LDLR,PCSK9,APOE \\
  --label my_crispri_v1
```

## Reuse rules

- Provide your own variant list and (optionally) experimental table — never
  commit confidential or pre-publication variant data.
- The script caches Catalog responses under `Data/Cache/AdvancedVariantAnnotations/`;
  delete the cache to force re-fetch.
- The ENCODE cCRE BED is downloaded once from `ENCODE_CCRE_BED_URL`. Override
  the env var to use a private mirror.
- Logistic models use IRLS over numpy/scipy only; results match `glm()`
  asymptotically.
- Plots are SVG so they can be edited in any vector editor.
""",
        encoding="utf-8",
    )
    print(f"Wrote skill doc: {path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Advanced variant annotation + integrated functional analysis.")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="Annotate a variant list and produce models, plots, and a report.")
    run.add_argument("--input", required=True, help="Input variant CSV.")
    run.add_argument("--experimental", help="Optional CSV of experimental columns to merge by variant key.")
    run.add_argument("--join-col", help="Column in --experimental to join on (default: auto-detect variant_id/rsid).")
    run.add_argument("--outcome", help="Column in the merged matrix to model as a binary outcome.")
    run.add_argument("--predictors", help=f"Comma-separated predictors (default: {','.join(PFUNC_FLAGS)}).")
    run.add_argument("--gene-list", help="Comma-separated gene symbols for per-gene Miami plots.")
    run.add_argument("--effect-col", help="Effect-size column for volcano plot.")
    run.add_argument("--bottom-p", help="Optional second p-value column for the Miami bottom track.")
    run.add_argument("--label", help="Output label (default: input filename stem).")
    run.add_argument("--max-rows", type=int, help="Limit rows for smoke runs.")
    run.add_argument("--sleep", type=float, default=0.05, help="Sleep between Catalog requests.")
    run.set_defaults(func=cmd_run)

    sub.add_parser("write-playbook", help="Write the skill documentation file.").set_defaults(func=cmd_write_playbook)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
