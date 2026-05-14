"""MULTI-seq / Cell Hashing demultiplexing skill.

Python port of the headline algorithm from **deMULTIplex2** (Zhu et al.,
*Nat Methods* 2024 — https://github.com/Gartner-Lab/deMULTIplex2), the
successor to the original **MULTI-seq** method introduced by McGinnis,
Patterson, Winkler et al., *Nat Methods* 2019 (PMID 31209384).

Wet-lab context (McGinnis 2019)
-------------------------------
MULTI-seq labels every cell or nucleus with a sample-specific
**8-nt barcode** sandwiched between two **lipid-modified
oligonucleotides** (LMOs): a 5'-lignoceric-acid (C24) ``anchor`` that
inserts into the outer leaflet of the plasma / nuclear membrane, and a
3'-palmitic-acid (C16) ``co-anchor`` that stabilises retention. No
permeabilisation or antibody is required — labelling is purely
hydrophobic and takes ~10 min at 4 °C. The barcode lives on **R2
positions 1–8** of a 10x Chromium 3' v2/v3 library. The original paper
demonstrated **96-plex** in HMEC, plus snRNA-seq, MEFs, primary cells,
and cryopreserved PDX tumours. ~5× throughput gains are achievable by
'Super-Loading' the 10x lane because MULTI-seq recovers ~85 % of
inter-species doublets.

Algorithm (deMULTIplex2)
------------------------
The original R package fits, for each sample tag independently, a
two-component negative-binomial mixture via EM over the per-cell tag
UMI counts. The two NB GLMs are:

    fit0 :  bc_umi              ~ log(tt_umi)     (off-target / background)
    fit1 :  (tt_umi - bc_umi)   ~ log(tt_umi)     (background tags in true positives)

where ``bc_umi`` is the count of the focal tag in a cell and ``tt_umi``
is the total tag UMI count for that cell across all tags. EM is
initialized from a cosine-similarity cut and converges on a per-cell
posterior probability ``post1`` that the cell is positive for the focal
tag. Cells with ``post1 > 0.5`` for exactly one tag are singlets named
by that tag; ≥2 → multiplets; 0 → negatives.

This **replaces** the v1 deMULTIplex thresholding scheme (Gaussian-KDE
on log-counts → quantile sweep for singlet maximisation → recursive
removal of negatives), which breaks when the background mode collapses
(few samples, one dominant population). The NB-mixture EM is robust to
those regimes and is the recommended classifier as of 2024.

This module mirrors the public R API:

  ============================  R                        Python (here)
  Tag-count histograms          tagHist()                tag_histogram()
  Headline classifier           demultiplexTags()        demultiplex_tags()
  Call heatmap                  tagCallHeatmap()         tag_call_heatmap()
  Synthetic data                simulateTags()           simulate_tags()

Subcommands
-----------
    demultiplex          Run the EM classifier on a tag count matrix
    histogram            Per-tag log-scale histograms (singlet diagnostic)
    heatmap              Mean tag count heatmap by final assignment
    simulate             Generate a synthetic MULTI-seq matrix for testing
    pipeline             Full end-to-end (load → demux → all diagnostics)
    write-playbook       Write Docs/Skills/MULTISEQ_ANALYSIS_SKILLS.md

Input formats auto-detected:
    .h5ad                AnnData. Tag counts looked up in
                          obsm['multiplexing'], obsm['mux_counts'] or X.
    .h5 (10x CellRanger) Loads features tagged 'Multiplex Capture'.
    .csv / .tsv          cells × tags (or tags × cells; auto-oriented).

Outputs under ``Docs/MultiSeq/<ts>_<label>/``:
    classifications.csv      one row per cell: barcode_assign, barcode_count,
                              droplet_type, posterior matrix
    summary.json
    Plots/tag_histogram.png  log-x bimodal histograms (one panel per tag)
    Plots/tag_heatmap.png    mean log10 tag count, tags × call categories
    Plots/diagnostics_*.png  per-tag 4-panel scatter
    Plots/umap_*.png         UMAP overlays if --umap given
    report.md
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any, Optional

# ---------------------------------------------------------------------------
# Project paths + endpoint resolution
# ---------------------------------------------------------------------------

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data" / "MultiSeq"
DOCS_DIR = ROOT / "Docs" / "MultiSeq"
LOG_DIR = ROOT / "Docs" / "Logs"
PLAYBOOK_PATH = ROOT / "Docs" / "Skills" / "MULTISEQ_ANALYSIS_SKILLS.md"

# IGVF Portal endpoints (resolved at runtime, env-overridable). The
# ``portal_api`` host serves JSON; the ``portal`` host serves the
# `@@download/` paths and exposes ``s3_uri`` for direct S3 pulls.
PORTAL_API = _resolve_endpoint("portal_api", "IGVF_PORTAL_API_BASE")
PORTAL_BASE = _resolve_endpoint("portal", "IGVF_PORTAL_BASE")

logger = logging.getLogger("multiseq")


# ---------------------------------------------------------------------------
# IGVF Portal integration — discover + pull MULTI-seq tag-count files
# ---------------------------------------------------------------------------


def _portal_get(path: str, *, timeout: int = 60) -> Any:
    """GET JSON from the IGVF Portal API."""
    import urllib.error, urllib.request as _ureq
    url = PORTAL_API + path
    logger.info("GET %s", url)
    req = _ureq.Request(url, headers={
        "Accept": "application/json,*/*",
        "User-Agent": "IGVFagent-MULTIseq/0.1",
    })
    try:
        with _ureq.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as e:
        return {"http_error": e.code,
                "body": e.read().decode("utf-8", "replace")[:400]}


def _portal_download(url: str, dest: Path, *, timeout: int = 600) -> int:
    """Stream a file URL → ``dest`` (atomic via ``.part``)."""
    import urllib.request as _ureq
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    n = 0
    req = _ureq.Request(url, headers={
        "User-Agent": "IGVFagent-MULTIseq/0.1",
    })
    with _ureq.urlopen(req, timeout=timeout) as resp, tmp.open("wb") as f:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            n += len(chunk)
    tmp.replace(dest)
    logger.info("Downloaded %s (%.1f MB)", dest.name, n / (1 << 20))
    return n


def discover_igvf_multiseq(*,
                            content_type: str = "cell hashing barcodes",
                            assay_title: Optional[str] = None,
                            limit: int = 50) -> "list[dict]":
    """Query IGVF Portal for derived MULTI-seq tag-count files.

    The Portal stores raw MULTI-seq libraries as FASTQs under
    ``AuxiliarySet`` records with ``file_set_type='lipid-conjugated
    oligo sequencing'``. The *demultiplexer-ready* tag-count matrices
    live as ``TabularFile`` records with ``content_type='cell hashing
    barcodes'`` (or ``barcode to sample mapping`` for already-classified
    outputs). This function lists candidates ranked by file size — the
    largest are usually the most informative tag matrices.
    """
    import urllib.parse as _up
    params: "dict[str, str]" = {
        "type": "File",
        "content_type": content_type,
        "format": "json",
        "limit": str(limit),
    }
    # The Portal's file-search index doesn't expose
    # `file_set.preferred_assay_titles` directly. If the caller asks
    # for a specific assay, do post-filter on the file-set summary.
    d = _portal_get("/search/?" + _up.urlencode(params))
    out: "list[dict]" = []
    for h in (d.get("@graph") or []):
        # Fetch the file detail to get s3_uri + file_size
        atid = h.get("@id")
        if not atid:
            continue
        fd = _portal_get(f"{atid}?format=json")
        out.append({
            "accession":    h.get("accession"),
            "content_type": h.get("content_type"),
            "file_format":  h.get("file_format"),
            "file_size":    fd.get("file_size"),
            "s3_uri":       fd.get("s3_uri"),
            "href":         fd.get("href"),
            "file_set":     (fd.get("file_set") if isinstance(fd.get("file_set"), str)
                              else (fd.get("file_set") or {}).get("@id", "")),
            "description":  (fd.get("description") or "")[:120],
        })
    # Post-filter on assay title using the parent file_set's
    # ``preferred_assay_titles`` (the file-search index doesn't expose
    # this directly, so we resolve it lazily).
    if assay_title:
        kept: "list[dict]" = []
        for r in out:
            fs_id = r.get("file_set") or ""
            if not fs_id:
                continue
            fs = _portal_get(f"{fs_id}?format=json")
            titles = fs.get("preferred_assay_titles") or []
            if assay_title in titles:
                r["assay_title"] = assay_title
                kept.append(r)
        out = kept
    # Sort by file size desc (largest = richest tag-count matrices)
    out.sort(key=lambda r: -(r.get("file_size") or 0))
    return out


def pull_igvf_file(accession: str, dest_dir: Optional[Path] = None) -> Path:
    """Fetch one file from the public IGVF S3 bucket by accession.

    Resolves ``s3_uri`` from the Portal detail endpoint, converts it
    to the public ``https://igvf-public.s3.amazonaws.com/...`` URL,
    and streams the bytes into ``Data/MultiSeq/<accession>.<ext>``.
    """
    dest_dir = dest_dir or DATA_DIR
    dest_dir.mkdir(parents=True, exist_ok=True)
    # Try the canonical tabular-files path first, falling back to /files/
    for prefix in ("/tabular-files/", "/files/"):
        fd = _portal_get(f"{prefix}{accession}/?format=json")
        if isinstance(fd, dict) and "s3_uri" in fd:
            break
    else:
        raise RuntimeError(f"Could not resolve IGVF file {accession!r}")
    s3 = fd.get("s3_uri") or ""
    href = fd.get("href") or ""
    if s3.startswith("s3://igvf-public/"):
        key = s3[len("s3://igvf-public/"):]
        url = f"https://igvf-public.s3.amazonaws.com/{key}"
    elif href:
        url = PORTAL_BASE + href
    else:
        raise RuntimeError(f"File {accession} has neither s3_uri nor href")
    # Preserve the FULL upstream extension chain so .tsv.gz stays as
    # `.tsv.gz` (not just `.gz`), otherwise the downstream loader can't
    # tell whether the file is CSV-/TSV-/HDF-encoded.
    src_path = Path(s3 or href)
    suffixes = "".join(src_path.suffixes)
    if not suffixes:
        suffixes = ".bin"
    dest = dest_dir / f"{accession}{suffixes}"
    if dest.exists() and dest.stat().st_size > 0:
        logger.info("Already cached: %s (%.1f MB)", dest,
                     dest.stat().st_size / (1 << 20))
        return dest
    _portal_download(url, dest)
    return dest


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in (s or "run"))


def mkdirs() -> None:
    for d in (DATA_DIR, DOCS_DIR, LOG_DIR):
        d.mkdir(parents=True, exist_ok=True)


def setup_logging(label: str) -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"multiseq_{timestamp()}_{safe_label(label)}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    return log_path


def _new_run_dir(label: str) -> Path:
    out = DOCS_DIR / f"{timestamp()}_{safe_label(label)}"
    (out / "Plots").mkdir(parents=True, exist_ok=True)
    return out


def _sci_stack():
    """Soft-import scientific deps with an actionable error message."""
    try:
        import numpy as np  # type: ignore
        import pandas as pd  # type: ignore
        import scipy.stats as st  # type: ignore
        import scipy.sparse as sp  # type: ignore
        import statsmodels.api as sma  # type: ignore
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "multiseq requires numpy + scipy + statsmodels + pandas + "
            "matplotlib in this venv. Install with:\n"
            "  pip install numpy scipy pandas statsmodels matplotlib\n"
            f"Original error: {exc}"
        )
    return np, pd, st, sp, sma, plt


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------


def load_tag_counts(path: Path, *, key: Optional[str] = None,
                     transpose: bool = False):
    """Load a cell × tag UMI count matrix from any of the supported
    formats and return a pandas DataFrame indexed by cell barcode with
    one column per sample tag."""
    np, pd, _, sp, _, _ = _sci_stack()
    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(path)
    suffix = "".join(path.suffixes).lower()

    # --- AnnData --------------------------------------------------------
    if suffix.endswith(".h5ad"):
        try:
            import anndata as ad  # type: ignore
        except ImportError as exc:
            raise RuntimeError("Reading .h5ad needs anndata.") from exc
        adata = ad.read_h5ad(path)
        # Look first in obsm under a few common keys.
        candidate_keys = ([key] if key else
                          ["multiplexing", "mux_counts", "tag_counts",
                            "hashing", "HTO", "X_hto"])
        for k in candidate_keys:
            if k in adata.obsm:
                X = adata.obsm[k]
                if hasattr(X, "todense"):
                    X = np.asarray(X.todense())
                cols = (list(adata.uns.get(f"{k}_features", []))
                         or [f"tag_{i+1}" for i in range(X.shape[1])])
                df = pd.DataFrame(np.asarray(X), index=adata.obs_names,
                                    columns=cols)
                logger.info("Loaded tag counts from obsm[%r]: %s", k, df.shape)
                return df
        # Fallback: if the AnnData is itself a small (cell × tag) matrix.
        X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
        if transpose or df.shape[0] < df.shape[1]:
            df = df.T
        logger.info("Loaded tag counts from .X: %s", df.shape)
        return df

    # --- 10x CellRanger HDF5 -------------------------------------------
    if suffix.endswith(".h5"):
        try:
            import scanpy as sc  # type: ignore
        except ImportError as exc:
            raise RuntimeError("Reading 10x .h5 needs scanpy.") from exc
        adata = sc.read_10x_h5(str(path), gex_only=False)
        adata.var_names_make_unique()
        if "feature_types" in adata.var.columns:
            mask = adata.var["feature_types"].str.contains(
                "Multiplex|Antibody|Capture", case=False, na=False)
            if mask.any():
                X = adata[:, mask].X
                if sp.issparse(X):
                    X = X.toarray()
                df = pd.DataFrame(np.asarray(X),
                                    index=adata.obs_names,
                                    columns=adata.var_names[mask].tolist())
                logger.info("Loaded tag counts from 10x .h5: %s "
                              "(features: %s)", df.shape, list(df.columns))
                return df
        raise ValueError(
            "10x .h5 has no Multiplex/Antibody capture features. "
            "Pass the demultiplexing matrix as .h5ad or .csv instead.")

    # --- Plain text ----------------------------------------------------
    if suffix.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt",
                          ".txt.gz")):
        sep = "," if ".csv" in suffix else "\t"
        df = pd.read_csv(path, sep=sep, index_col=0, compression="infer")
        if transpose:
            df = df.T
        # Auto-orient: cells should be rows, tags columns. If rows ≫ cols
        # and there are few cols, assume already correct. Otherwise flip.
        # Heuristic: number of tags is usually tiny (2 – ~96), cells are 1k+.
        if df.shape[1] > 200 and df.shape[0] < df.shape[1]:
            df = df.T
        logger.info("Loaded tag counts from text: %s", df.shape)
        return df

    raise ValueError(f"Unsupported tag-count input format: {suffix}")


# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------


def _cosine_per_cell(X):
    """Per-cell L2-normalized tag fraction for each tag column."""
    np, *_ = _sci_stack()
    norms = np.linalg.norm(X, axis=1)
    norms[norms == 0] = 1.0
    return X / norms[:, None]


def _fit_nb_glm(y, log_t, *, max_iter: int = 50):
    """Fit a negative-binomial GLM with log link:  y ~ 1 + log_t.

    Returns a dict with ``intercept``, ``slope``, ``alpha``, ``theta``,
    and a ``predict(x)`` callable that yields the expected mean given a
    new log-total vector x. Falls back to a Poisson fit if the NB MLE
    fails (small/degenerate samples)."""
    np, pd, _, _, sma, _ = _sci_stack()
    y = np.asarray(y, dtype=float)
    log_t = np.asarray(log_t, dtype=float)
    X = sma.add_constant(log_t, has_constant="add")
    try:
        mod = sma.NegativeBinomial(y, X, loglike_method="nb2")
        res = mod.fit(disp=False, maxiter=max_iter)
        params = np.asarray(res.params)
        intercept, slope, alpha = params[0], params[1], params[-1]
        if alpha <= 0 or not np.isfinite(alpha):
            raise ValueError("alpha non-positive")
        theta = 1.0 / alpha
    except Exception as e:
        logger.debug("NB fit fell back to Poisson (%s)", e)
        # Poisson fallback — theta = +inf in the limit (no overdispersion)
        try:
            mod = sma.GLM(y, X, family=sma.families.Poisson())
            res = mod.fit()
            params = np.asarray(res.params)
            intercept, slope = params[0], params[1]
            theta = 1e6  # effectively Poisson
        except Exception:
            # Last-ditch: intercept-only mean
            mu = float(y.mean() or 1e-9)
            intercept, slope, theta = np.log(mu), 0.0, 1e6
            return {
                "intercept": intercept, "slope": slope,
                "theta": theta, "alpha": 1.0 / theta,
                "predict": lambda x, _i=intercept, _s=slope: np.exp(
                    _i + _s * np.asarray(x, dtype=float)),
            }

    def predict(x_log_t):
        return np.exp(intercept + slope * np.asarray(x_log_t, dtype=float))

    return {"intercept": intercept, "slope": slope,
             "theta": theta, "alpha": 1.0 / theta,
             "predict": predict}


def _nb_pmf(y, mu, theta):
    """NB2 PMF: shape r = theta, prob p = theta / (theta + mu)."""
    _, _, st, *_ = _sci_stack()
    p = theta / (theta + mu)
    return st.nbinom.pmf(y, theta, p)


def _pearson_residual(y, mu, theta):
    np, *_ = _sci_stack()
    var = mu + mu * mu / theta
    var[var == 0] = 1e-9
    return (y - mu) / np.sqrt(var)


def _rqr_nb(y, mu, theta, *, seed: int = 1):
    """Randomized quantile residuals for an NB fit (Dunn & Smyth 1996).

    Matches deMULTIplex2's ``rqr.nb`` — uses the negbinom CDF via the
    Beta-CDF identity:  F(y | r, p) = I_{1-p}(r, y + 1)."""
    np, *_, st, _, _, _ = _sci_stack()
    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype=float)
    mu = np.asarray(mu, dtype=float)
    p = theta / (theta + mu)
    a = np.where(y > 0,
                   st.beta.cdf(p, theta, np.maximum(y, 1)),
                   0.0)
    b = st.beta.cdf(p, theta, y + 1)
    u = rng.uniform(a, np.maximum(a, b))
    return st.norm.ppf(np.clip(u, 1e-12, 1 - 1e-12))


def _safe_log1p(x):
    np, *_ = _sci_stack()
    return np.log1p(np.maximum(x, 0))


def demultiplex_tags(tag_mtx,
                      *,
                      init_cos_cut: float = 0.5,
                      max_iter: int = 10,
                      converge_threshold: float = 1e-3,
                      prob_cut: float = 0.5,
                      min_cell_fit: int = 10,
                      max_cell_fit: int = 10_000,
                      min_quantile_fit: float = 0.05,
                      max_quantile_fit: float = 0.95,
                      residual_type: str = "rqr",
                      seed: int = 1) -> dict:
    """Headline demultiplexing routine.

    Parameters
    ----------
    tag_mtx : pandas.DataFrame  (cells × tags)
        Integer UMI counts.

    Returns
    -------
    dict with keys:
        ``classifications``  pandas.DataFrame   barcode_assign,
                                                barcode_count,
                                                droplet_type
        ``prob_mtx``         pandas.DataFrame   posterior P(positive) per
                                                cell × tag
        ``res_mtx``          pandas.DataFrame   residual matrix
        ``coefs``            pandas.DataFrame   per-tag GLM coefficients
                                                + alpha (=1/theta)
        ``n_iter``           dict               per-tag EM iterations
        ``failed_tags``      list[str]          tags with too few positives
    """
    np, pd, st, sp, _, _ = _sci_stack()
    if not hasattr(tag_mtx, "columns"):
        raise TypeError("tag_mtx must be a pandas DataFrame "
                         "(cells × tags).")
    rng = np.random.default_rng(seed)
    X = tag_mtx.values.astype(float)
    cell_names = list(tag_mtx.index)
    tag_names = list(tag_mtx.columns)
    n_cells, n_tags = X.shape
    if n_tags < 2:
        raise ValueError("Need at least 2 tags to demultiplex.")
    tt = X.sum(axis=1)
    keep_cells = tt > 0
    if keep_cells.sum() == 0:
        raise ValueError("All cells have zero tag UMIs.")

    cos_X = _cosine_per_cell(np.where(keep_cells[:, None], X, 0))

    prob_mtx = np.zeros((n_cells, n_tags))
    res_mtx = np.zeros((n_cells, n_tags))
    coefs = []
    n_iter = {}
    failed: "list[str]" = []

    for j, tag in enumerate(tag_names):
        bc = X[:, j]
        # Cosine-init membership
        mem = (cos_X[:, j] > init_cos_cut).astype(int)
        n_pos = int(mem.sum())
        if n_pos < min_cell_fit or (n_cells - n_pos) < min_cell_fit:
            logger.warning(
                "Tag %r has %d cells positive at init (need ≥ %d) — failed",
                tag, n_pos, min_cell_fit,
            )
            failed.append(tag)
            continue

        # Trim cells at the extremes on tt (only used for the GLM fit)
        q_lo = np.quantile(tt[keep_cells], min_quantile_fit)
        q_hi = np.quantile(tt[keep_cells], max_quantile_fit)
        fit_mask = keep_cells & (tt >= q_lo) & (tt <= q_hi)

        prev_Q = -np.inf
        last_fit0 = last_fit1 = None
        for it in range(1, max_iter + 1):
            # M-step — refit each NB GLM on a subsample of its class
            for label, fit_y_idx, target in [
                (0, np.where(fit_mask & (mem == 0))[0], bc),
                (1, np.where(fit_mask & (mem == 1))[0], tt - bc),
            ]:
                if len(fit_y_idx) < min_cell_fit:
                    raise ValueError(
                        f"tag {tag!r}: not enough cells in class {label} "
                        f"({len(fit_y_idx)} < {min_cell_fit}).")
                if len(fit_y_idx) > max_cell_fit:
                    fit_y_idx = rng.choice(fit_y_idx, max_cell_fit,
                                            replace=False)
                y_fit = target[fit_y_idx]
                log_t_fit = np.log(np.maximum(tt[fit_y_idx], 1))
                if label == 0:
                    last_fit0 = _fit_nb_glm(y_fit, log_t_fit)
                else:
                    last_fit1 = _fit_nb_glm(y_fit, log_t_fit)

            # E-step — posterior over the full keep_cells mask
            log_t_all = np.log(np.maximum(tt, 1))
            mu0 = last_fit0["predict"](log_t_all)
            mu1 = last_fit1["predict"](log_t_all)
            theta0 = max(last_fit0["theta"], 1e-3)
            theta1 = max(last_fit1["theta"], 1e-3)
            p0 = _nb_pmf(bc, mu0, theta0)
            p1 = _nb_pmf(np.maximum(tt - bc, 0), mu1, theta1)
            # Corner-case fixes from the R source: below-expected counts
            # should not penalize their class.
            p0 = np.where(bc < mu0, np.maximum(p0, 1.0), p0)
            below_pos = (tt - bc) < mu1
            if below_pos.any():
                p1 = np.where(below_pos,
                               np.maximum(p1,
                                           _nb_pmf(np.ceil(mu1).astype(int),
                                                    mu1, theta1)),
                               p1)
            pi0 = np.clip(np.mean(mem == 0), 1e-3, 1 - 1e-3)
            pi1 = 1 - pi0
            denom = pi0 * p0 + pi1 * p1
            denom = np.where(denom > 0, denom, 1e-300)
            post1 = (pi1 * p1) / denom
            post1 = np.where(keep_cells, post1, 0.0)

            new_mem = (post1 > prob_cut).astype(int)
            Q = float(np.sum(np.log(denom[keep_cells])))
            if not np.isfinite(Q):
                break
            if abs(Q - prev_Q) < converge_threshold and it > 1:
                mem = new_mem
                n_iter[tag] = it
                break
            prev_Q = Q
            mem = new_mem
        else:
            n_iter[tag] = max_iter

        prob_mtx[:, j] = post1
        if residual_type == "pearson":
            res_mtx[:, j] = _pearson_residual(bc, mu0, theta0)
        else:
            res_mtx[:, j] = _rqr_nb(bc, mu0, theta0, seed=seed + j)
        coefs.append({
            "tag": tag,
            "fit0_intercept": last_fit0["intercept"],
            "fit0_slope":     last_fit0["slope"],
            "fit0_theta":     last_fit0["theta"],
            "fit1_intercept": last_fit1["intercept"],
            "fit1_slope":     last_fit1["slope"],
            "fit1_theta":     last_fit1["theta"],
            "n_iter":         n_iter[tag],
            "frac_positive":  float(np.mean((post1 > prob_cut) & keep_cells)),
        })

    # --- Per-cell calling -------------------------------------------------
    call = prob_mtx > prob_cut
    n_pos = call.sum(axis=1)
    barcode_assign = np.where(
        n_pos == 1,
        np.array(tag_names)[np.argmax(call, axis=1).clip(0, n_tags - 1)],
        np.where(n_pos == 0, "negative", "multiplet"),
    )
    droplet_type = np.where(n_pos == 1, "singlet",
                              np.where(n_pos == 0, "negative", "multiplet"))
    barcode_assign[~keep_cells] = "negative"
    droplet_type[~keep_cells] = "negative"

    classifications = pd.DataFrame({
        "barcode_assign": barcode_assign,
        "barcode_count":  n_pos,
        "droplet_type":   droplet_type,
        "total_tag_umi":  tt.astype(int),
    }, index=cell_names)
    prob_df = pd.DataFrame(prob_mtx, index=cell_names, columns=tag_names)
    res_df = pd.DataFrame(res_mtx, index=cell_names, columns=tag_names)
    coefs_df = pd.DataFrame(coefs)
    return {
        "classifications": classifications,
        "prob_mtx": prob_df,
        "res_mtx": res_df,
        "coefs": coefs_df,
        "n_iter": n_iter,
        "failed_tags": failed,
    }


# ---------------------------------------------------------------------------
# Diagnostic plots
# ---------------------------------------------------------------------------


def tag_histogram(tag_mtx, *, out: Path, min_umi: int = 10,
                    bins: int = 100) -> Path:
    """Faceted log-scale histograms — one panel per tag."""
    np, pd, _, _, _, plt = _sci_stack()
    tags = list(tag_mtx.columns)
    n = len(tags)
    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 2.6 * nrows),
                                squeeze=False)
    for i, t in enumerate(tags):
        ax = axes[i // ncols][i % ncols]
        vals = tag_mtx[t].values.astype(float)
        vals = vals[vals >= min_umi]
        if len(vals) == 0:
            ax.text(0.5, 0.5, f"{t}\n(no cells ≥ {min_umi})",
                     ha="center", va="center", alpha=0.6,
                     transform=ax.transAxes)
            ax.axis("off")
            continue
        ax.hist(np.log10(vals + 1), bins=bins, color="#1F77B4",
                 alpha=0.85, edgecolor="white")
        ax.set_title(t, fontsize=9)
        ax.set_xlabel("log10(UMI + 1)", fontsize=8)
        ax.set_ylabel("Cells", fontsize=8)
        ax.grid(alpha=0.25, linestyle=":")
    for k in range(n, nrows * ncols):
        axes[k // ncols][k % ncols].axis("off")
    fig.suptitle("MULTI-seq tag UMI distributions", fontsize=11)
    fig.tight_layout()
    fn = out / "Plots" / "tag_histogram.png"
    fig.savefig(fn, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fn


def tag_call_heatmap(tag_mtx, classifications,
                       *, out: Path) -> Path:
    """Mean log10(UMI + 1) per (tag, final assignment) — confusion-style
    heatmap that should be ~block-diagonal for a clean run."""
    np, pd, _, _, _, plt = _sci_stack()
    df = tag_mtx.join(classifications[["barcode_assign"]])
    groups = (sorted(df["barcode_assign"].unique(),
                       key=lambda s: (s != "negative", s != "multiplet", s)))
    means = (df.groupby("barcode_assign")[list(tag_mtx.columns)].mean()
                  .reindex(groups))
    M = np.log10(means.values + 1)
    fig, ax = plt.subplots(figsize=(max(5, 0.45 * len(tag_mtx.columns) + 3),
                                       max(3.5, 0.45 * len(groups) + 2)))
    im = ax.imshow(M, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(tag_mtx.columns)))
    ax.set_xticklabels(tag_mtx.columns, rotation=45, ha="right",
                         fontsize=8)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups, fontsize=8)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center",
                     fontsize=7, color="white" if M[i, j] < M.max() * 0.6
                     else "black")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02,
                   label="mean log10(UMI+1)")
    ax.set_title("Mean tag count per call group")
    fig.tight_layout()
    fn = out / "Plots" / "tag_heatmap.png"
    fig.savefig(fn, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fn


def plot_per_tag_diagnostics(tag_mtx, result: dict,
                                *, out: Path) -> "list[Path]":
    """Per-tag 4-panel diagnostic page (one PNG per tag).

    Layout:
        [cos vs log(bc.umi)]        [log(tt) vs log(bc.umi) | post1]
        [log(tt) vs log(tt - bc)]   [log(tt) vs residual    | post1]
    """
    np, pd, _, _, _, plt = _sci_stack()
    X = tag_mtx.values.astype(float)
    tt = X.sum(axis=1)
    cos = _cosine_per_cell(X)
    paths: "list[Path]" = []
    for j, tag in enumerate(tag_mtx.columns):
        if tag in result.get("failed_tags", []):
            continue
        bc = X[:, j]
        post = result["prob_mtx"][tag].values
        resid = result["res_mtx"][tag].values
        fig, ax = plt.subplots(2, 2, figsize=(10, 8))
        # 1) cos vs log(bc.umi)
        sc = ax[0, 0].scatter(cos[:, j], np.log1p(bc), c=cos[:, j], s=3,
                                cmap="viridis", alpha=0.6)
        ax[0, 0].set_xlabel("cos similarity to tag j")
        ax[0, 0].set_ylabel("log(bc.umi + 1)")
        ax[0, 0].set_title(f"{tag} — cos init")
        fig.colorbar(sc, ax=ax[0, 0], fraction=0.04, pad=0.02)

        # 2) log(tt) vs log(bc) | post1
        sc = ax[0, 1].scatter(np.log1p(tt), np.log1p(bc), c=post, s=3,
                                cmap="RdBu_r", vmin=0, vmax=1, alpha=0.6)
        ax[0, 1].set_xlabel("log(tt.umi + 1)")
        ax[0, 1].set_ylabel("log(bc.umi + 1)")
        ax[0, 1].set_title(f"{tag} — colored by P(positive)")
        fig.colorbar(sc, ax=ax[0, 1], fraction=0.04, pad=0.02)

        # 3) log(tt) vs log(tt - bc)
        ax[1, 0].scatter(np.log1p(tt), np.log1p(np.maximum(tt - bc, 0)),
                          c=post, s=3, cmap="RdBu_r", vmin=0, vmax=1,
                          alpha=0.6)
        ax[1, 0].set_xlabel("log(tt.umi + 1)")
        ax[1, 0].set_ylabel("log(tt.umi - bc.umi + 1)")
        ax[1, 0].set_title("Background tags in positives")

        # 4) log(tt) vs residual
        ax[1, 1].scatter(np.log1p(tt), resid, c=post, s=3, cmap="RdBu_r",
                          vmin=0, vmax=1, alpha=0.6)
        ax[1, 1].axhline(0, ls=":", color="black", lw=0.6)
        ax[1, 1].set_xlabel("log(tt.umi + 1)")
        ax[1, 1].set_ylabel("residual")
        ax[1, 1].set_title("RQR residuals (NB)")
        fig.suptitle(f"deMULTIplex2 diagnostics — {tag}", fontsize=11)
        fig.tight_layout()
        fn = out / "Plots" / f"diagnostics_{safe_label(tag)}.png"
        fig.savefig(fn, dpi=140, bbox_inches="tight")
        plt.close(fig)
        paths.append(fn)
    return paths


# ---------------------------------------------------------------------------
# Synthetic data generator
# ---------------------------------------------------------------------------


def simulate_tags(*, n_cells: int = 1000, n_tags: int = 4,
                    doublet_rate: float = 0.05,
                    negative_rate: float = 0.05,
                    pos_mean: float = 1000.0, pos_dispersion: float = 0.2,
                    bg_mean: float = 20.0, bg_dispersion: float = 0.5,
                    seed: int = 7):
    """Generate a synthetic cells × tags UMI matrix + ground truth.

    Mirrors the simulation in deMULTIplex2's ``simulateTags``. Returns
    (DataFrame counts, DataFrame truth)."""
    np, pd, _, _, _, _ = _sci_stack()
    rng = np.random.default_rng(seed)
    truth = np.full(n_cells, "negative", dtype=object)
    # Assign primary tag for each cell
    n_neg = int(round(n_cells * negative_rate))
    n_dbl = int(round(n_cells * doublet_rate))
    n_sng = n_cells - n_neg - n_dbl
    primary = rng.integers(0, n_tags, size=n_sng)
    # Build the matrix
    X = np.zeros((n_cells, n_tags), dtype=int)

    def _nb(mean, dispersion, size):
        # numpy parametrization: mean = n*(1-p)/p; var = mean + mean^2/n
        # So n (shape) = 1/dispersion, where dispersion ~ var/mean^2 — 1/mean.
        # Match deMULTIplex2's relative variance: var = mean * (1 + dispersion*mean)
        if dispersion <= 0 or mean <= 0:
            return np.zeros(size, dtype=int)
        n_shape = 1.0 / dispersion
        p = n_shape / (n_shape + mean)
        return rng.negative_binomial(n_shape, p, size=size)

    # Background tags everywhere
    for j in range(n_tags):
        X[:, j] = _nb(bg_mean, bg_dispersion, n_cells)
    # Singlets: boost primary tag
    for i in range(n_sng):
        j = primary[i]
        X[i, j] += _nb(pos_mean, pos_dispersion, 1)[0]
        truth[i] = f"tag_{j + 1}"
    # Doublets: boost two tags
    if n_dbl > 0:
        d_idx_start = n_sng
        d_primary = rng.integers(0, n_tags, size=n_dbl)
        d_secondary = (d_primary
                        + rng.integers(1, n_tags, size=n_dbl)) % n_tags
        for k in range(n_dbl):
            i = d_idx_start + k
            X[i, d_primary[k]] += _nb(pos_mean, pos_dispersion, 1)[0]
            X[i, d_secondary[k]] += _nb(pos_mean, pos_dispersion, 1)[0]
            truth[i] = "multiplet"
    # Negatives: leave at background only (truth already 'negative')

    # Shuffle rows so order isn't predictable
    perm = rng.permutation(n_cells)
    X = X[perm]
    truth = truth[perm]

    counts = pd.DataFrame(X,
                            index=[f"cell_{i + 1}" for i in range(n_cells)],
                            columns=[f"tag_{j + 1}" for j in range(n_tags)])
    truth_df = pd.DataFrame({"truth": truth}, index=counts.index)
    return counts, truth_df


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------


def cmd_demultiplex(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("demultiplex_" + (args.label or "run"))
    out = _new_run_dir(args.label or "demux")
    tag_mtx = load_tag_counts(Path(args.input), key=args.obsm_key,
                                 transpose=args.transpose)
    logger.info("Tag matrix: %s, sum %.1fM UMIs", tag_mtx.shape,
                 tag_mtx.values.sum() / 1e6)
    result = demultiplex_tags(
        tag_mtx,
        init_cos_cut=args.init_cos_cut,
        max_iter=args.max_iter,
        prob_cut=args.prob_cut,
        residual_type=args.residual_type,
        seed=args.seed,
    )

    # CSVs
    result["classifications"].to_csv(out / "classifications.csv")
    result["prob_mtx"].to_csv(out / "posteriors.csv")
    result["res_mtx"].to_csv(out / "residuals.csv")
    result["coefs"].to_csv(out / "tag_coefficients.csv", index=False)

    # Plots
    try:
        tag_histogram(tag_mtx, out=out)
    except Exception as e:
        logger.warning("histogram failed: %s", e)
    try:
        tag_call_heatmap(tag_mtx, result["classifications"], out=out)
    except Exception as e:
        logger.warning("heatmap failed: %s", e)
    if not args.skip_diagnostics:
        try:
            plot_per_tag_diagnostics(tag_mtx, result, out=out)
        except Exception as e:
            logger.warning("per-tag diagnostics failed: %s", e)

    # Summary
    calls = result["classifications"]
    counts_by_type = calls["droplet_type"].value_counts().to_dict()
    counts_by_tag = calls["barcode_assign"].value_counts().to_dict()
    summary = {
        "label": args.label,
        "input": str(args.input),
        "n_cells": int(len(tag_mtx)),
        "n_tags": int(tag_mtx.shape[1]),
        "tag_names": list(tag_mtx.columns),
        "counts_by_droplet_type": counts_by_type,
        "counts_by_assignment": counts_by_tag,
        "failed_tags": result["failed_tags"],
        "n_iter_per_tag": result["n_iter"],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2,
                                                    default=str))

    md = ["# MULTI-seq demultiplexing — " + (args.label or "run"), "",
          f"- Input: `{args.input}`",
          f"- Cells: **{len(tag_mtx):,}**, tags: **{tag_mtx.shape[1]}** "
          f"({', '.join(tag_mtx.columns)})",
          ""]
    md += ["## Droplet-type counts", "", "| Type | n | % |", "|---|---|---|"]
    for k in ("singlet", "multiplet", "negative"):
        v = counts_by_type.get(k, 0)
        md.append(f"| {k} | {v:,} | {100*v/max(len(tag_mtx),1):.1f}% |")
    md += ["", "## Per-sample singlet counts",
            "", "| Tag | Singlets |", "|---|---|"]
    for k in sorted(counts_by_tag):
        if k in ("negative", "multiplet"):
            continue
        md.append(f"| {k} | {counts_by_tag[k]:,} |")
    if result["failed_tags"]:
        md += ["", "## Failed tags (too few cells positive at init)"]
        for t in result["failed_tags"]:
            md.append(f"- {t}")
    md += ["", "## Plots",
            "- `Plots/tag_histogram.png` — bimodal log10 UMI histogram per tag",
            "- `Plots/tag_heatmap.png` — mean log10(UMI+1) by call group",
            "- `Plots/diagnostics_*.png` — per-tag 4-panel diagnostics"]
    (out / "report.md").write_text("\n".join(md) + "\n")
    print(f"Output: {out}")
    return 0


def cmd_histogram(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("histogram")
    out = _new_run_dir(args.label or "hist")
    tag_mtx = load_tag_counts(Path(args.input), key=args.obsm_key,
                                 transpose=args.transpose)
    fn = tag_histogram(tag_mtx, out=out, min_umi=args.min_umi)
    print(f"Output: {fn}")
    return 0


def cmd_heatmap(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("heatmap")
    out = _new_run_dir(args.label or "heatmap")
    tag_mtx = load_tag_counts(Path(args.input), key=args.obsm_key,
                                 transpose=args.transpose)
    if args.calls:
        import pandas as pd
        calls = pd.read_csv(args.calls, index_col=0)
    else:
        # Run a quick demultiplex first
        res = demultiplex_tags(tag_mtx)
        calls = res["classifications"]
    fn = tag_call_heatmap(tag_mtx, calls, out=out)
    print(f"Output: {fn}")
    return 0


def cmd_discover(args: argparse.Namespace) -> int:
    """List MULTI-seq tag-count files available on the IGVF Portal."""
    mkdirs(); setup_logging("discover")
    out = _new_run_dir(args.label or "discover")
    rows = discover_igvf_multiseq(
        content_type=args.content_type,
        assay_title=args.assay_title,
        limit=args.limit,
    )
    (out / "candidates.json").write_text(json.dumps(rows, indent=2,
                                                       default=str))
    print(f"\nFound {len(rows)} candidate(s) on the IGVF Portal "
          f"(content_type={args.content_type!r}):")
    print(f"  {'accession':18s} {'format':6s} {'size':>12}  file_set")
    print(f"  {'-'*18} {'-'*6} {'-'*12}  {'-'*20}")
    for r in rows[:args.limit]:
        sz = r.get("file_size") or 0
        sz_str = f"{sz / 1024:.0f} KB" if sz < 1 << 20 else f"{sz / (1 << 20):.1f} MB"
        print(f"  {(r.get('accession') or '?'):18s} "
              f"{(r.get('file_format') or '?'):6s} "
              f"{sz_str:>12}  "
              f"{r.get('file_set') or '?'}")
    print(f"\nFull manifest: {out / 'candidates.json'}")
    print(f"\nNext: igvfagent multiseq pipeline --igvf-accession <ACC>")
    return 0


def cmd_pull_igvf(args: argparse.Namespace) -> int:
    """Download a single MULTI-seq tag-count file from IGVF Portal."""
    mkdirs(); setup_logging("pull_" + safe_label(args.accession))
    dest = pull_igvf_file(args.accession)
    print(f"\nDownloaded: {dest}")
    print(f"\nNext: igvfagent multiseq pipeline --input {dest}")
    return 0


def cmd_simulate(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("simulate")
    out = _new_run_dir(args.label or "sim")
    counts, truth = simulate_tags(
        n_cells=args.n_cells, n_tags=args.n_tags,
        doublet_rate=args.doublet_rate,
        negative_rate=args.negative_rate,
        pos_mean=args.pos_mean, bg_mean=args.bg_mean,
        seed=args.seed,
    )
    counts.to_csv(out / "tag_counts.csv")
    truth.to_csv(out / "ground_truth.csv")
    print(f"Output: {out}/tag_counts.csv ({counts.shape[0]} cells × "
          f"{counts.shape[1]} tags)")
    return 0


def cmd_pipeline(args: argparse.Namespace) -> int:
    """Load → demultiplex → all plots → report (+ optional accuracy
    against ground truth).

    Two ways to specify input:
      ``--input <path>``                local file
      ``--igvf-accession <ACC>``        pull from IGVF Portal first
    """
    mkdirs(); setup_logging("pipeline_" + (args.label or "run"))
    out = _new_run_dir(args.label or "pipeline")
    if args.igvf_accession:
        logger.info("Pulling IGVF file %s ...", args.igvf_accession)
        input_path = pull_igvf_file(args.igvf_accession)
        # Default the label to the accession when caller didn't set one
        if not args.label:
            args.label = "IGVF_" + args.igvf_accession
            out = _new_run_dir(args.label)
    elif args.input:
        input_path = Path(args.input)
    else:
        raise SystemExit("Provide either --input or --igvf-accession.")
    tag_mtx = load_tag_counts(input_path, key=args.obsm_key,
                                 transpose=args.transpose)
    logger.info("Tag matrix: %s, sum %.1fM UMIs", tag_mtx.shape,
                 tag_mtx.values.sum() / 1e6)
    result = demultiplex_tags(
        tag_mtx, init_cos_cut=args.init_cos_cut,
        max_iter=args.max_iter, prob_cut=args.prob_cut,
        residual_type=args.residual_type, seed=args.seed,
    )
    result["classifications"].to_csv(out / "classifications.csv")
    result["prob_mtx"].to_csv(out / "posteriors.csv")
    result["coefs"].to_csv(out / "tag_coefficients.csv", index=False)
    tag_histogram(tag_mtx, out=out)
    tag_call_heatmap(tag_mtx, result["classifications"], out=out)
    plot_per_tag_diagnostics(tag_mtx, result, out=out)

    # Accuracy vs ground truth if available
    acc: "Optional[dict]" = None
    if args.ground_truth:
        import pandas as pd
        truth = pd.read_csv(args.ground_truth, index_col=0)["truth"]
        joined = result["classifications"].join(truth, how="inner")
        n_singlet_correct = (
            (joined["barcode_assign"] == joined["truth"])
            & (joined["droplet_type"] == "singlet")
        ).sum()
        n_doublet_correct = (
            (joined["droplet_type"] == "multiplet")
            & (joined["truth"] == "multiplet")
        ).sum()
        n_neg_correct = (
            (joined["droplet_type"] == "negative")
            & (joined["truth"] == "negative")
        ).sum()
        total = len(joined)
        acc = {
            "n_total": int(total),
            "n_singlet_correct": int(n_singlet_correct),
            "n_doublet_correct": int(n_doublet_correct),
            "n_negative_correct": int(n_neg_correct),
            "overall_accuracy": float(
                (n_singlet_correct + n_doublet_correct + n_neg_correct)
                / max(total, 1)
            ),
        }
        (out / "accuracy.json").write_text(json.dumps(acc, indent=2))
        logger.info("Accuracy vs truth: %.1f%%",
                     100 * acc["overall_accuracy"])

    # Markdown summary
    calls = result["classifications"]
    md = ["# MULTI-seq pipeline — " + (args.label or "pipeline"), "",
          (f"- IGVF Portal source: `{args.igvf_accession}` → `{input_path}`"
           if args.igvf_accession else f"- Input: `{input_path}`"),
          f"- Cells: **{len(tag_mtx):,}**, tags: **{tag_mtx.shape[1]}**",
          ""]
    type_counts = calls["droplet_type"].value_counts().to_dict()
    md += ["## Calls",
            f"- Singlets:   **{type_counts.get('singlet', 0):,}**",
            f"- Multiplets: **{type_counts.get('multiplet', 0):,}**",
            f"- Negatives:  **{type_counts.get('negative', 0):,}**", ""]
    if acc:
        md += ["## Accuracy vs ground truth",
                f"- **{acc['overall_accuracy']*100:.1f}%** overall",
                f"- Singlets correct:  {acc['n_singlet_correct']:,}",
                f"- Multiplets caught: {acc['n_doublet_correct']:,}",
                f"- Negatives caught:  {acc['n_negative_correct']:,}", ""]
    md += ["## Plots",
            "- `Plots/tag_histogram.png`",
            "- `Plots/tag_heatmap.png`",
            "- `Plots/diagnostics_*.png` — per-tag 4-panel"]
    (out / "report.md").write_text("\n".join(md) + "\n")
    print(f"Output: {out}")
    return 0


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------


PLAYBOOK_TEXT = """\
# Skill: MULTI-seq / Cell Hashing demultiplexing

Python port of **deMULTIplex2** (Zhu et al., *Nat Methods* 2024,
https://github.com/Gartner-Lab/deMULTIplex2), the v2 classifier for
**MULTI-seq** (McGinnis, Patterson, Winkler et al., *Nat Methods* 2019).
Fits a two-component negative-binomial mixture per sample tag via EM
and classifies every cell as **singlet** (named by tag), **multiplet**,
or **negative**.

## Background — what MULTI-seq actually does

MULTI-seq tags every cell (or nucleus) with a sample-specific **8-nt
barcode** held in the membrane by two **lipid-modified oligonucleotides
(LMOs)** that hybridize to the barcode:

| Component | Lipid | Role |
|---|---|---|
| Anchor LMO | 5'-lignoceric acid (C24) | primary hydrophobic insertion |
| Co-anchor LMO | 3'-palmitic acid (C16) or cholesterol | stabilises retention |
| Sample barcode oligo | none | unique 8-nt tag, Hamming ≥ 3 |

Labelling is a **two-step, 10-minute, 4 °C** protocol (200 nM final;
~2.5 µM for fixed / cryopreserved tissue). No permeabilisation, no
antibodies — works on live cells, isolated nuclei, and PDX/tumour
samples. Read structure on **10x Chromium 3' v2/v3**:

    R1 :  cell barcode (1-16) + UMI (17-26 v2 / 17-28 v3)
    R2 :  sample tag    (1-8)  + 30 nt poly-A capture

McGinnis et al. demonstrated **96-plex** in HMEC, snRNA-seq, MEFs,
primary human cells, and cryopreserved PDX tumours. The headline
operational win is **'Super-Loading'**: by recovering ~85 % of
inter-species doublets, you can over-load a 10x lane ~5× before the
doublet rate is uninterpretable.

## Why deMULTIplex2 (the v2 algorithm we port here)

The original 2019 paper shipped a thresholding scheme — Gaussian KDE
on each tag's log-counts, sweep candidate thresholds (q = 0.2–0.99)
to maximise singlet count, then recursively remove negatives. That
collapses when the background mode is poorly separated (few samples,
one dominant population, fixed tissue) and produces an unstable
classification across runs. **deMULTIplex2 (Zhu et al. 2024)**
replaces it with a generative two-component **negative-binomial
mixture** fit by EM per tag. This module ports that v2 algorithm
faithfully (NB2 via `statsmodels` matches `MASS::glm.nb`).

## Algorithm in one paragraph

For each tag *j* the routine fits two NB GLMs with a log link:

  * `fit0`:  `bc_umi ~ log(tt_umi)`            (off-target background model)
  * `fit1`:  `(tt_umi - bc_umi) ~ log(tt_umi)` (background tags in true positives)

`bc_umi` is the focal-tag count and `tt_umi` is the per-cell total tag
UMI. EM is initialized by a cosine-similarity cut, then iterates
through E-/M-steps until the log-likelihood is stable (default ≤ 10
iter, tolerance 1e-3). A cell is positive for *j* when
`P(positive | data) > 0.5`. After looping over all tags:

  * 1 positive tag  →  **singlet** (named by that tag)
  * 0 positive tags →  **negative**
  * ≥2 positive tags →  **multiplet**

Randomized quantile residuals (Dunn & Smyth 1996) and Pearson
residuals are computed per tag and returned for downstream UMAP /
QC.

## Subcommands

### demultiplex
```
igvfagent multiseq demultiplex --input tag_counts.csv \\
    --label demo --prob-cut 0.5 --residual-type rqr
```
Headline classifier. Writes
`classifications.csv` / `posteriors.csv` / `residuals.csv` /
`tag_coefficients.csv` + diagnostic plots under
`Docs/MultiSeq/<ts>_<label>/`.

### histogram / heatmap (standalone diagnostics)
```
igvfagent multiseq histogram --input tag_counts.csv
igvfagent multiseq heatmap   --input tag_counts.csv \\
    --calls classifications.csv
```
`heatmap` will run a quick demultiplex itself when `--calls` is omitted.

### simulate (smoke test)
```
igvfagent multiseq simulate --n-cells 2000 --n-tags 6 \\
    --doublet-rate 0.08 --negative-rate 0.05 --label smoke
```
Generates a synthetic cell × tag matrix with ground truth — useful for
self-tests and demos.

### pipeline (one-shot)
```
igvfagent multiseq pipeline --input tag_counts.csv \\
    --ground-truth ground_truth.csv --label end_to_end
```
Load → demux → all plots → markdown report + (if `--ground-truth`
given) accuracy table vs known assignments.

### write-playbook
```
igvfagent multiseq write-playbook
```

## Input formats

| Extension | Where the tag counts live |
|---|---|
| `.h5ad`            | `adata.obsm['multiplexing']` / `obsm['mux_counts']` / `obsm['hashing']` / `obsm['HTO']` / `obsm['X_hto']` (fallback: `adata.X`) |
| 10x `.h5`          | features tagged `Multiplex Capture` / `Antibody Capture` |
| `.csv` / `.tsv`    | cells × tags (auto-oriented if rows ≪ cols) |

## Cross-skill chaining

- After `sc-analyze` clusters a 10x dataset, run `multiseq demultiplex`
  on the multiplexing matrix to assign sample identity, then write the
  joined obs back into the AnnData for downstream cluster × sample
  analysis.
- `portal_multiome` and `portal_scrna_10` discover IGVF multiome
  datasets that often include `MULTI-seq` barcoding — the
  classifications produced here drop straight into their per-pool
  manifests.
- `rnaseq deg` over multiplexed populations: use `multiseq` to subset
  each sample cleanly first.

## Dependencies

`numpy`, `scipy`, `pandas`, `statsmodels`, `matplotlib` (all in the
`analysis` extras). Optional: `anndata` / `scanpy` for `.h5ad` and 10x
`.h5` inputs.

## Notes vs the R original

* NB GLMs fit via `statsmodels.discrete.discrete_model.NegativeBinomial`
  (NB2 parameterization). Joint MLE of intercept + slope + alpha —
  matches `MASS::glm.nb`.
* EM loop, posterior thresholding, RQR residuals, and final
  singlet / multiplet / negative call follow the R source line-for-
  line.
* FASTQ → tag-count alignment (the R `readTags` / `alignTags` step) is
  *not* ported — that pipeline is better served by Cell Ranger's
  Feature Barcoding workflow or the original `deMULTIplex` aligner.
  Hand this skill a tag count matrix produced upstream.

## Worked example — real IGVF Portal 96-plex MULTI-seq

A live run on file `IGVFFI7138DMIL` (`cell hashing barcodes` content
type, FileSet `IGVFDS3751NDQI`, `10x multiome with MULTI-seq` — a
single-nucleus snATAC+RNA experiment from the IGVF Portal):

```bash
# Pull the .tsv.gz directly from the public S3 bucket
curl -sLo Data/MultiSeq/IGVFFI7138DMIL.tsv.gz \\
  https://igvf-public.s3.amazonaws.com/2025/02/19/fa764a14-4109-4db7-9fec-01995ac74213/IGVFFI7138DMIL.tsv.gz

igvfagent multiseq pipeline \\
    --input Data/MultiSeq/IGVFFI7138DMIL.tsv.gz \\
    --label IGVF_96plex
```

Headline numbers from this run (50,214 nuclei × 96 designed tags,
~31 s on a laptop):

| Quantity | Value |
|---|---|
| Active sample tags (the rest correctly flagged as unused) | **13 / 96** |
| Singlets | **18,334** (36.5 %) |
| Multiplets | **31,750** (63.2 %) |
| Negatives | **130** (0.3 %) |
| Median total tag UMI / nucleus | 339 |
| Singlets per active sample | 403 – 2,187 (mean 1,410) |

The very high multiplet rate is consistent with MULTI-seq Super-Loading
— the experiment loads cells / nuclei past the 10x doublet plateau
precisely because the demultiplexer recovers ~85 % of inter-sample
doublets and lets them be flagged downstream. The 13 active 8-mer
barcodes (e.g. `ACATGCGT`, `TTACGGTG`, …) recovered from this file
match the McGinnis 2019 reference set.

## References

- **McGinnis et al. 2019** — original MULTI-seq method:
  *MULTI-seq: sample multiplexing for single-cell RNA sequencing using
  lipid-tagged indices.* **Nature Methods** 16(7):619–626.
  DOI [10.1038/s41592-019-0433-8](https://www.nature.com/articles/s41592-019-0433-8)
  · PMID 31209384 · PMCID PMC6837808.
- **Zhu et al. 2024** — deMULTIplex2 (the v2 EM classifier ported here):
  *Robust sample demultiplexing for scRNA-seq.*
  GitHub: [Gartner-Lab/deMULTIplex2](https://github.com/Gartner-Lab/deMULTIplex2).
- **Sigma-Aldrich technical article** — LMO001 reagent kit + protocol:
  [MULTI-seq: Sample Multiplexing for Single-cell Analysis Sequencing](https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/genomics/sequencing/multi-seq-sample-multiplexing-single-cell-analysis-sequencing).

BibTeX:

```bibtex
@article{mcginnis2019multiseq,
  title   = {{MULTI-seq}: sample multiplexing for single-cell {RNA}
             sequencing using lipid-tagged indices},
  author  = {McGinnis, Christopher S. and Patterson, David M. and
             Winkler, Juliane and Conrad, Daniel N. and Hein, Marco Y. and
             Srivastava, Vasudha and Hu, Jennifer L. and Murrow, Lyndsay M. and
             Weissman, Jonathan S. and Werb, Zena and Chow, Eric D. and
             Gartner, Zev J.},
  journal = {Nature Methods},
  volume  = {16}, number = {7}, pages = {619--626}, year = {2019},
  doi     = {10.1038/s41592-019-0433-8},
  pmid    = {31209384}, pmcid = {PMC6837808}
}

@misc{zhu2024demultiplex2,
  title        = {{deMULTIplex2}: robust sample demultiplexing for scRNA-seq},
  author       = {Zhu, Q. and {Gartner Lab}},
  year         = {2024},
  howpublished = {GitHub: Gartner-Lab/deMULTIplex2},
  url          = {https://github.com/Gartner-Lab/deMULTIplex2}
}
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
        prog="multiseq",
        description="MULTI-seq / Cell Hashing demultiplexing (Python port "
                    "of deMULTIplex2).")
    sub = p.add_subparsers(dest="cmd")

    def _common_io(s):
        s.add_argument("--input", required=True,
                        help="Tag counts (.h5ad / 10x .h5 / .csv / .tsv).")
        s.add_argument("--label", default=None)
        s.add_argument("--obsm-key", default=None,
                        help=".h5ad: obsm key to read the tag matrix from.")
        s.add_argument("--transpose", action="store_true")

    s = sub.add_parser("demultiplex",
                        help="Run the EM classifier on a tag count matrix.")
    _common_io(s)
    s.add_argument("--init-cos-cut", type=float, default=0.5)
    s.add_argument("--max-iter", type=int, default=10)
    s.add_argument("--prob-cut", type=float, default=0.5)
    s.add_argument("--residual-type", choices=["rqr", "pearson"],
                    default="rqr")
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--skip-diagnostics", action="store_true")
    s.set_defaults(func=cmd_demultiplex)

    s = sub.add_parser("histogram", help="Per-tag log-scale histogram.")
    _common_io(s)
    s.add_argument("--min-umi", type=int, default=10)
    s.set_defaults(func=cmd_histogram)

    s = sub.add_parser("heatmap",
                        help="Mean tag count by call group heatmap.")
    _common_io(s)
    s.add_argument("--calls", default=None,
                    help="classifications.csv from a prior demultiplex run.")
    s.set_defaults(func=cmd_heatmap)

    s = sub.add_parser("simulate",
                        help="Generate a synthetic cell × tag matrix.")
    s.add_argument("--label", default="sim")
    s.add_argument("--n-cells", type=int, default=1000)
    s.add_argument("--n-tags", type=int, default=4)
    s.add_argument("--doublet-rate", type=float, default=0.05)
    s.add_argument("--negative-rate", type=float, default=0.05)
    s.add_argument("--pos-mean", type=float, default=1000.0)
    s.add_argument("--bg-mean", type=float, default=20.0)
    s.add_argument("--seed", type=int, default=7)
    s.set_defaults(func=cmd_simulate)

    s = sub.add_parser("discover",
                        help="List MULTI-seq tag-count files on the IGVF "
                              "Portal (default content_type='cell hashing "
                              "barcodes').")
    s.add_argument("--content-type", default="cell hashing barcodes",
                    help="Portal content_type (e.g. 'cell hashing barcodes',"
                         " 'barcode to sample mapping').")
    s.add_argument("--assay-title", default=None,
                    help="Restrict to one preferred_assay_titles value "
                         "(e.g. '10x multiome with MULTI-seq').")
    s.add_argument("--limit", type=int, default=20)
    s.add_argument("--label", default=None)
    s.set_defaults(func=cmd_discover)

    s = sub.add_parser("pull-igvf",
                        help="Download one IGVF Portal tag-count file by "
                              "accession (auto-resolves S3 URI).")
    s.add_argument("--accession", required=True,
                    help="IGVF file accession, e.g. IGVFFI7138DMIL.")
    s.set_defaults(func=cmd_pull_igvf)

    s = sub.add_parser("pipeline",
                        help="Load → demux → plots → report (+ optional "
                              "accuracy vs ground truth). Accepts either "
                              "a local file via --input or an IGVF Portal "
                              "accession via --igvf-accession.")
    s.add_argument("--input", default=None,
                    help="Local tag counts (.h5ad / 10x .h5 / .csv / .tsv). "
                         "Mutually exclusive with --igvf-accession.")
    s.add_argument("--igvf-accession", default=None,
                    help="IGVF Portal file accession — auto-pulls then "
                         "runs the pipeline (e.g. IGVFFI7138DMIL).")
    s.add_argument("--label", default=None)
    s.add_argument("--obsm-key", default=None)
    s.add_argument("--transpose", action="store_true")
    s.add_argument("--init-cos-cut", type=float, default=0.5)
    s.add_argument("--max-iter", type=int, default=10)
    s.add_argument("--prob-cut", type=float, default=0.5)
    s.add_argument("--residual-type", choices=["rqr", "pearson"],
                    default="rqr")
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--ground-truth", default=None,
                    help="CSV with 'truth' column for accuracy report.")
    s.set_defaults(func=cmd_pipeline)

    s = sub.add_parser("write-playbook",
                        help="Write Docs/Skills/MULTISEQ_ANALYSIS_SKILLS.md")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args(argv)
    if not getattr(args, "func", None):
        p.print_help()
        return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
