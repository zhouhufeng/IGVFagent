#!/usr/bin/env python3
"""SCEPTRE for IGVF CRISPR MuData: guide assignment, calibration, power and discovery (port of IGVF-CRISPR/sceptreIGVF).

Port of https://github.com/IGVF-CRISPR/sceptreIGVF (MIT, "sceptreIGVF authors",
Eugene Katsevich), the R wrapper that runs SCEPTRE (Barry, Roeder & Katsevich)
on the IGVF CRISPR single-cell MuData schema: a `gene` modality (UMI counts),
a `guide` modality (gRNA UMI counts, optional `guide_assignment` layer,
`intended_target_name` / `targeting` in var, `moi` in uns), cell covariates in
the top-level obs and `pairs_to_test` in the top-level uns.  Every R function
of the package (R/mudata.R, R/grna-assignment.R, R/inference.R, the tests and
the sample-data script) was read, and the functions were re-derived in Python;
no code was copied.  Relationship: port.  The SCEPTRE statistics themselves
live in the `sceptre` R package, which sceptreIGVF calls; they are
re-implemented here in numpy + statsmodels as a documented approximation (see
"SCEPTRE statistics").  `r-runner` calls the real R packages when Rscript,
MuData and sceptreIGVF/sceptre are installed.

The IGVF CRISPR_Pipeline evolved the wrapper (bin/inference_sceptre.R,
bin/assign_grnas_sceptre.R, bin/merge_cis_trans_results.py): greedy collinear
covariate removal, the complement control group, union (per-element) and
singleton (per-guide) analyses written to uns `per_element_results` /
`per_guide_results`, and their cis_/trans_ prefixed forms.  `inference
--mode pipeline` reproduces that schema too.

Definitions reproduced from sceptreIGVF
  convert (convert_mudata_to_sceptre_object)
      response matrix = mod/gene counts (first assay), gRNA matrix =
      mod/guide `guide_assignment` layer when present, else the counts;
      grna_target = guide var `intended_target_name`; moi = guide uns `moi`;
      extra covariates = top-level obs.  remove_collinear_covariates: build
      model.matrix(~ .) of the covariates; if rank < ncol, drop ALL extra
      covariates (sceptreIGVF, --collinear-mode all-or-nothing) -- or, as the
      pipeline does, drop invariant columns and keep covariates greedily while
      the design stays full rank (--collinear-mode greedy).
  assign-guides (assign_grnas_sceptre)
      convert with remove_collinear_covariates = TRUE, empty discovery set,
      sceptre::assign_grnas(method = "mixture"); the 0/1 matrix is written
      back as layer `guide_assignment` of mod/guide (and guide_assignment.mtx,
      guides x cells, as the pipeline does).  Defaults probability_threshold
      0.8, n_em_rep 5.
  inference (inference_sceptre)
      discovery pairs = uns pairs_to_test renamed intended_target_name ->
      grna_target, gene_id -> response_id; formula_object default
      ~ log(response_n_nonzero) + log(response_n_umis); gRNA assignment by
      thresholding at 1 UMI; run_qc(n_nonzero_trt_thresh = 0,
      n_nonzero_cntrl_thresh = 0, p_mito_threshold = 1) (the cell-wise
      response_n_umis / response_n_nonzero 1%-99% range filters still apply);
      run_discovery_analysis; test_results = pairs_to_test LEFT JOIN
      (p_value, log2_fc = log_2_fold_change) on (intended_target_name,
      gene_id), stored in uns["test_results"] with the columns of
      pairs_to_test followed by p_value, log2_fc.
  run (sceptre workflow + sceptre_object_to_mudata)
      assign (sceptre default: mixture in high MOI, maximum in low MOI),
      default QC, calibration check, power check on positive controls,
      discovery; output MuData: gene obs num_expressed_genes /
      total_gene_umis, guide obs num_expressed_guides / total_guide_umis,
      guide var targeting ("TRUE"/"FALSE", grna_target != "non-targeting"),
      intended_target_name, intended_target_chr ("" if NA),
      intended_target_start / _end (-9 if NA), guide uns moi, top-level obs
      = covariates without the computed ones, uns pairs_to_test
      (intended_target_name, gene_id, pair_type; positive controls first) and
      test_results (+ p_value, log2_fc).
  export-mudata (sceptre_object_to_mudata_inputs_outputs + save_mudata_list)
      positive controls = power-check pairs with a result; discovery pairs =
      min(#significant, round(N/2)) significant + the rest non-significant,
      sampled; rerun with default assignment and QC; eight objects
      inference_output / inference_input (no test_results) /
      guide_assignment_output (no pairs_to_test) / guide_assignment_input (no
      guide_assignment layer), each also *_minimal (no covariates, no gene
      var/obs, guide var = targeting + intended_target_name, pairs_to_test =
      intended_target_name + gene_id, test_results + p_value); written to
      <out>/guide_assignment/<prefix><name>.h5mu and <out>/inference/...

SCEPTRE statistics (re-implemented; defaults as in sceptre 0.10)
  covariates      response_n_nonzero, response_n_umis, response_p_mito (genes
                  ^MT-), grna_n_nonzero, grna_n_umis, + extra covariates.
                  auto formula: log(x) of the count covariates (log(x + 1) if
                  any zero), response_p_mito, gRNA covariates only in high
                  MOI, extra covariates; invariant terms dropped.
  assignment      thresholding: count >= threshold (default 5; 1 in
                  inference).  maximum (low MOI): the gRNA with most UMIs if
                  it holds >= umi_fraction_threshold (0.5) of the cell's gRNA
                  UMIs and the cell has >= min_grna_n_umis_threshold (5)
                  UMIs.  mixture: per gRNA, Poisson GLM of its counts on the
                  response-side covariates gives background means mu_i; EM on
                  y_i ~ pi Pois(f mu_i) + (1 - pi) Pois(mu_i), n_em_rep random
                  starts (pi in [1e-5, 0.1], log f in [0.5, 10]); assigned if
                  posterior >= probability_threshold and y_i > 0; gRNAs with
                  < 10 non-zero cells fall back to thresholding at 5.
  cell-wise QC    drop cells outside the [1%, 99%] quantiles of
                  response_n_umis and response_n_nonzero, response_p_mito >
                  p_mito_threshold (0.2), and in low MOI cells with zero or
                  2+ assigned gRNAs.
  pairwise QC     n_nonzero_trt / n_nonzero_cntrl = treated / control cells
                  with a non-zero response; pass_qc if both >= 7 (0 in
                  inference); failing pairs get p_value NA.
  grouping        union (per element: any gRNA of the target), singleton
                  (each gRNA), bonferroni (singleton p-values, min p x k).
  control group   complement (high MOI default; all other cells) or nt_cells
                  (low MOI default; cells carrying a non-targeting gRNA).
  test            NB GLM of the response on the formula covariates (theta by
                  ML on the Poisson fit, then statsmodels NB GLM with that
                  theta), fitted on all cells (complement) or on the NT cells
                  (nt_cells); GLM score statistic of the treatment indicator
                  z = sum x_i r_i / sqrt(x'Wx - x'WZ(Z'WZ)^-1Z'Wx),
                  r_i = (y_i - mu_i) / (1 + mu_i / theta),
                  w_i = mu_i / (1 + mu_i / theta).
  null            permutations (low MOI): the treatment labels permuted within
                  treatment + control cells; crt (high MOI): conditional
                  resampling x_i ~ Bernoulli(pi_i), pi_i from a logistic
                  regression of the indicator on the covariates.  The same
                  resamples are reused across responses of a gRNA group.
  p-value         B1 = 499 resamples; if the empirical p <= 0.02, B2 = 4999
                  resamples and a skew-normal fitted by the method of moments
                  (sceptre fits a skew-normal too); if the fit is poor (KS
                  distance > 0.03) or resampling_approximation =
                  no_approximation, B3 = 24999 resamples, empirical.
                  Empirical p = (1 + #{null at least as extreme}) / (1 + B);
                  side left / right / both (2 min(left, right)).
  fold change     log_2_fold_change = MLE of b in y_i ~ NB(mu_i e^b, theta)
                  over treated cells, / log 2.
  significance    BH over pairs passing QC, significant = q < 0.1.
  calibration     "undercover" NT groups (NT gRNAs partitioned into groups of
                  the median targeting-group size) x random responses,
                  n_calibration_pairs = number of discovery pairs passing QC;
                  in nt_cells mode the rest of the NT cells are the control.
  cis / trans     construct_cis_pairs: targets vs genes whose TSS is within
                  500 kb of the target; construct_trans_pairs: every
                  targeting group x every response.

Deviations (documented; not bugs in the port)
  * sceptre's null resampling fits a skew-normal inside C++ with its own
    check; ours uses a moment fit + KS check.  p-values agree in distribution
    (calibrated nulls, same order of magnitude on signals), not digit-for-digit.
  * sceptreIGVF's inference_sceptre tests `formula_object == "default"`
    with a non-short-circuit `|` on an undefined variable, so passing no
    formula_object errors in R unless the name exists; here `--formula
    default` (or omitted) means the documented default formula.
  * mixture assignment uses response-side covariates only (the gRNA
    covariates are functions of the counts being modelled).
  * log_2_fold_change uses the null-model means as offset; with the
    complement control those means are fitted on all cells, so strong effects
    in large treatment groups are shrunk towards 0 (planted log2 -1.73 is
    estimated about -1.5 in the selftest).
  * assign-guides always models the raw gRNA counts (X); upstream would
    re-assign an existing guide_assignment layer if the input carried one.

Subcommands
  convert             MuData -> sceptre object directory (matrices,
                      grna_target_data_frame.tsv, covariate_data_frame.tsv,
                      formula.txt) with the collinearity check.
  assign-guides       assign_grnas_sceptre: mixture / thresholding / maximum,
                      writes the guide_assignment layer + .mtx.
  qc                  cell-wise and pairwise QC summary for a MuData.
  inference           inference_sceptre (sceptreIGVF) or the pipeline's
                      per-element/per-guide variant (--mode pipeline, --scope
                      cis|trans).
  calibration-check   negative-control (NT) pairs; false discoveries, QQ.
  power-check         positive-control pairs.
  run                 the whole SCEPTRE workflow; output MuData in the
                      sceptre_object_to_mudata schema.
  export-mudata       sceptre_object_to_mudata_inputs_outputs +
                      save_mudata_list (the eight benchmark MuData files).
  r-runner            run sceptreIGVF in R (Rscript + MuData + sceptreIGVF);
                      prints the command when R is absent.
  compare-results     two test_results tables (e.g. Python vs R): p-value
                      rank agreement and significance concordance.
  selftest            synthetic high- and low-MOI screens with planted
                      effects; every subcommand asserted.

Output: Docs/SceptreIGVF/<timestamp>_<label>/.  numpy, pandas, scipy,
statsmodels, patsy, anndata and h5py required; mudata optional (h5mu is
read and written through h5py + anndata when it is absent); matplotlib
optional.

Usage:
    igvfagent sceptre-igvf convert --mudata inference_input.h5mu --remove-collinear-covariates
    igvfagent sceptre-igvf assign-guides --mudata guide_assignment_input.h5mu --label gasperini
    igvfagent sceptre-igvf inference --mudata inference_input.h5mu --side left --formula "~ prep_batch + log(response_n_nonzero) + log(response_n_umis)"
    igvfagent sceptre-igvf inference --mudata inference_input.h5mu --mode pipeline --scope cis
    igvfagent sceptre-igvf calibration-check --mudata inference_input.h5mu
    igvfagent sceptre-igvf power-check --mudata inference_input.h5mu
    igvfagent sceptre-igvf run --mudata screen.h5mu --label screen
    igvfagent sceptre-igvf export-mudata --mudata screen.h5mu --num-discovery-pairs 100 --prefix gasperini_
    igvfagent sceptre-igvf r-runner --step inference --mudata inference_input.h5mu
    igvfagent sceptre-igvf compare-results --a py_test_results.tsv --b r_test_results.tsv
    igvfagent sceptre-igvf selftest --no-plots
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import time
import warnings
import zlib
from collections import OrderedDict
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "SceptreIGVF"

UPSTREAM_REPO = "IGVF-CRISPR/sceptreIGVF"
UPSTREAM_COMMIT = "87a35d6736a5ddbf9c2975323d55e4423bd88219"

DEFAULT_INFERENCE_FORMULA = "~ log(response_n_nonzero) + log(response_n_umis)"
COMPUTED_COVARIATES = ["grna_n_nonzero", "grna_n_umis", "response_n_nonzero", "response_n_umis", "response_p_mito"]
RESULT_COLS = ["response_id", "grna_target", "n_nonzero_trt", "n_nonzero_cntrl", "pass_qc", "p_value",
               "log_2_fold_change", "significant"]
NT_LABEL = "non-targeting"

BLUE, ORANGE, INK, INK2, AXIS, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"

log = logging.getLogger("sceptre_igvf")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"sceptre_igvf_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(path)], force=True)
    return path


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in label)[:80]


def run_dir(label: str) -> Path:
    d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}"
    i = 1
    while d.exists():
        d = OUT_ROOT / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}_{i}"
        i += 1
    d.mkdir(parents=True, exist_ok=True)
    return d


def _np():
    import numpy as np  # type: ignore
    return np


def _pd():
    try:
        import pandas as pd  # type: ignore
        return pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit("sceptre-igvf needs numpy, pandas, scipy, statsmodels: pip install 'igvfagent[analysis]'") from e


def _sp():
    import scipy.sparse as sp  # type: ignore
    return sp


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
    fig.tight_layout()
    fig.savefig(str(path), dpi=130)
    import matplotlib.pyplot as plt  # type: ignore
    plt.close(fig)
    print(f"Figure: {path}")
    return path


def md_table(headers: "list[str]", rows: "Iterable[Iterable[Any]]", max_rows: int = 40) -> str:
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
            return "NA"
        v = float(x)
        if v != 0 and abs(v) < 10 ** (-nd):
            return f"{v:.2e}"
        return f"{v:.{nd}f}"
    except (TypeError, ValueError):
        return str(x)


def write_tsv(df, path: Path, announce: bool = True) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip" if str(path).endswith(".gz") else None)
    if announce:
        print(f"TSV: {path}")
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
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
    return str(o)


def write_report(d: Path, title: str, sections: "list[str]") -> Path:
    p = d / "report.md"
    p.write_text(f"# {title}\n\n" + "\n\n".join(sections) + "\n")
    print(f"Report: {p}")
    return p


def bh_adjust(pvals):
    np = _np()
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    n = q.size
    if n == 0:
        return out
    order = np.argsort(q)
    ranked = q[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


def _seed_for(*parts) -> int:
    return zlib.crc32("|".join(str(p) for p in parts).encode()) & 0x7FFFFFFF


# ---------------------------------------------------------------------------
# MuData I/O (mudata when installed, else h5py + anndata)
# ---------------------------------------------------------------------------

class Modalities:
    def __init__(self, gene, guide, obs, uns: dict, source: Optional[Path] = None, backend: str = "memory") -> None:
        self.gene, self.guide, self.obs, self.uns, self.source, self.backend = gene, guide, obs, uns, source, backend


def read_h5mu(path: Path) -> Modalities:
    try:
        from tf_perturb_stages import read_mudata  # type: ignore
        m = read_mudata(Path(path), need_uns=True)
        return Modalities(m.gene, m.guide, m.obs, dict(m.uns), Path(path), m.backend)
    except ImportError:
        pass
    import h5py  # type: ignore
    try:
        from anndata.io import read_elem  # type: ignore
    except ImportError:
        from anndata.experimental import read_elem  # type: ignore
    with h5py.File(str(path), "r") as h:
        mods = list(h["mod"].keys())
        gk = "gene" if "gene" in mods else mods[0]
        uk = "guide" if "guide" in mods else mods[-1]
        gene, guide = read_elem(h[f"mod/{gk}"]), read_elem(h[f"mod/{uk}"])
        obs = read_elem(h["obs"]) if "obs" in h else gene.obs.iloc[:, :0]
        uns = read_elem(h["uns"]) if "uns" in h else {}
    return Modalities(gene, guide, obs, dict(uns), Path(path), "h5py")


def write_h5mu(path: Path, gene, guide, obs, uns: dict, announce: bool = True) -> Path:
    """Write a two-modality MuData (gene, guide) with top-level obs and uns."""
    pd, np = _pd(), _np()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    obs = obs.copy() if obs is not None else pd.DataFrame(index=gene.obs_names)
    obs.index = obs.index.astype(str)
    try:
        import mudata as md  # type: ignore
        m = md.MuData({"gene": gene, "guide": guide})
        for c in obs.columns:
            m.obs[c] = obs.reindex(m.obs_names)[c].values
        for k, v in uns.items():
            m.uns[k] = v
        m.write(str(path))
    except ImportError:
        import h5py  # type: ignore
        try:
            from anndata.io import write_elem  # type: ignore
        except ImportError:
            from anndata.experimental import write_elem  # type: ignore
        with h5py.File(str(path), "w") as f:
            f.attrs["encoding-type"] = "MuData"
            f.attrs["encoding-version"] = "0.1.0"
            write_elem(f, "mod/gene", gene)
            write_elem(f, "mod/guide", guide)
            f["mod"].attrs["mod-order"] = np.array(["gene", "guide"], dtype=object).astype("S")
            write_elem(f, "obs", obs)
            var_index = list(gene.var_names.astype(str)) + list(guide.var_names.astype(str))
            write_elem(f, "var", pd.DataFrame(index=pd.Index(var_index).astype(str)))
            n_obs, n_g, n_u = len(obs), gene.n_vars, guide.n_vars
            gpos = {c: i + 1 for i, c in enumerate(gene.obs_names.astype(str))}
            upos = {c: i + 1 for i, c in enumerate(guide.obs_names.astype(str))}
            write_elem(f, "obsmap", {"gene": np.array([gpos.get(c, 0) for c in obs.index], dtype=np.uint32),
                                     "guide": np.array([upos.get(c, 0) for c in obs.index], dtype=np.uint32)})
            vg = np.zeros(n_g + n_u, dtype=np.uint32); vg[:n_g] = np.arange(1, n_g + 1)
            vu = np.zeros(n_g + n_u, dtype=np.uint32); vu[n_g:] = np.arange(1, n_u + 1)
            write_elem(f, "varmap", {"gene": vg, "guide": vu})
            write_elem(f, "obsm", {})
            write_elem(f, "varm", {})
            write_elem(f, "uns", _uns_clean(uns))
        _ = n_obs
    if announce:
        print(f"Wrote: {path}")
    return path


def _uns_clean(uns: dict) -> dict:
    """anndata's writer wants plain types; DataFrames get string/bool/float columns."""
    pd = _pd()
    out = {}
    for k, v in uns.items():
        if isinstance(v, pd.DataFrame):
            df = v.copy().reset_index(drop=True)
            for c in df.columns:
                if df[c].dtype == object:
                    df[c] = df[c].map(lambda x: "" if x is None or (isinstance(x, float) and math.isnan(x)) else str(x))
            df.index = df.index.astype(str)
            out[k] = df
        elif isinstance(v, dict):
            out[k] = _uns_clean(v)
        else:
            out[k] = v
    return out


def uns_table(uns: dict, key: str):
    pd = _pd()
    v = uns.get(key)
    if v is None:
        return None
    if isinstance(v, pd.DataFrame):
        return v.reset_index(drop=True)
    if isinstance(v, dict):
        return pd.DataFrame({k: (list(x) if not hasattr(x, "shape") else x) for k, x in v.items()})
    return pd.DataFrame(v)


def _csr(x):
    sp, np = _sp(), _np()
    if sp.issparse(x):
        return x.tocsr().astype(np.float64)
    return sp.csr_matrix(np.asarray(x, dtype=np.float64))


def _align(mods: Modalities) -> Modalities:
    """Subset gene and guide to the cells they share, in gene order."""
    g, u = mods.gene, mods.guide
    gn, un = list(g.obs_names.astype(str)), list(u.obs_names.astype(str))
    if gn == un:
        return mods
    common = [c for c in gn if c in set(un)]
    if not common:
        raise SystemExit("gene and guide modalities share no cell barcodes")
    log.warning("aligning modalities to %d shared cells", len(common))
    return Modalities(g[common].copy(), u[common].copy(), mods.obs.reindex(common), mods.uns, mods.source, mods.backend)


# ---------------------------------------------------------------------------
# convert_mudata_to_sceptre_object
# ---------------------------------------------------------------------------

def is_nt_target(t) -> bool:
    s = str(t).strip().lower().replace("_", "-").replace(" ", "-")
    return s.startswith("non-targeting") or s.startswith("nontargeting")


class SceptreData:
    """The Python analogue of a sceptre_object after import_data()."""

    def __init__(self, response, response_ids, response_names, grna, grna_ids, grna_var, grna_target, moi,
                 covariates, extra_covariate_names, cell_ids, grna_is_assignment: bool):
        self.response = response          # cells x genes, csc
        self.response_ids = response_ids
        self.response_names = response_names
        self.grna = grna                  # cells x guides, csc (counts or 0/1 assignment)
        self.grna_ids = grna_ids
        self.grna_var = grna_var
        self.grna_target = grna_target    # pandas Series grna_id -> grna_target
        self.moi = moi
        self.covariates = covariates
        self.extra_covariate_names = extra_covariate_names
        self.cell_ids = cell_ids
        self.grna_is_assignment = grna_is_assignment
        np = _np()
        tv = grna_var["targeting"].astype(str).str.upper().isin(["FALSE", "F", "0"]) if "targeting" in grna_var.columns else None
        nt = np.array([is_nt_target(t) for t in grna_target.values])
        if tv is not None:
            nt = nt | tv.values
        self.grna_is_nt = nt

    @property
    def low_moi(self) -> bool:
        return self.moi == "low"

    @property
    def n_cells(self) -> int:
        return self.response.shape[0]


def model_matrix_rank_ok(df) -> bool:
    """R's model.matrix(~ ., data) full-rank check (characters -> treatment-coded factors, intercept)."""
    np, pd = _np(), _pd()
    if df.shape[1] == 0:
        return True
    parts = [np.ones((len(df), 1))]
    for c in df.columns:
        col = df[c]
        if col.dtype == object or str(col.dtype) in ("category", "bool") or col.dtype == bool:
            d = pd.get_dummies(col.astype(str), drop_first=True).values.astype(float)
            if d.size:
                parts.append(d)
        else:
            v = pd.to_numeric(col, errors="coerce").astype(float)
            parts.append(v.fillna(v.median() if v.notna().any() else 0.0).values[:, None])
    M = np.hstack(parts)
    return int(np.linalg.matrix_rank(M)) == M.shape[1]


def prepare_extra_covariates(obs, remove_collinear: bool, mode: str = "all-or-nothing"):
    """sceptreIGVF: drop all extra covariates if model.matrix(~ .) is rank deficient.
    greedy (CRISPR_Pipeline): characters/logicals -> factors, drop invariant columns, keep covariates in order while
    the design remains full rank."""
    pd = _pd()
    if obs is None or obs.shape[1] == 0:
        return pd.DataFrame(index=obs.index if obs is not None else None), []
    cov = obs[[c for c in obs.columns if ":" not in str(c) and c not in COMPUTED_COVARIATES]].copy()
    for c in cov.columns:
        if str(cov[c].dtype) == "category":
            cov[c] = cov[c].astype(str)
        elif cov[c].dtype == bool:
            cov[c] = cov[c].astype(str)
    dropped: "list[str]" = []
    if mode == "greedy":
        varying = [c for c in cov.columns if cov[c].dropna().nunique() > 1]
        dropped += [c for c in cov.columns if c not in varying]
        cov = cov[varying]
        if remove_collinear and cov.shape[1] > 1:
            kept: "list[str]" = []
            for c in cov.columns:
                if model_matrix_rank_ok(cov[kept + [c]]):
                    kept.append(c)
                else:
                    dropped.append(c)
            cov = cov[kept]
    elif remove_collinear and cov.shape[1] > 0 and not model_matrix_rank_ok(cov):
        print("Removing multicollinear covariates")
        dropped = list(cov.columns)
        cov = cov.iloc[:, :0]
    return cov, dropped


def _gene_names(gene) -> "list[str]":
    for c in ("gene_name", "symbol", "gene_symbols", "feature_name", "gene_names"):
        if c in gene.var.columns:
            return [str(x) for x in gene.var[c].values]
    return [str(x) for x in gene.var_names]


def convert_mudata_to_sceptre(mods: Modalities, remove_collinear_covariates: bool = False,
                              collinear_mode: str = "all-or-nothing", target_col: str = "intended_target_name",
                              moi: Optional[str] = None, prefer_counts: bool = False) -> "tuple[SceptreData, list[str]]":
    np, pd = _np(), _pd()
    mods = _align(mods)
    gene, guide = mods.gene, mods.guide
    moi = moi or guide.uns.get("moi")
    if isinstance(moi, (list, tuple)) or hasattr(moi, "shape"):
        moi = str(list(np.ravel(moi))[0])
    if moi is None:
        log.warning("guide.uns['moi'] missing; assuming high MOI")
        moi = "high"
    moi = str(moi).lower()
    if moi not in ("low", "high"):
        raise SystemExit(f"moi must be 'low' or 'high', got {moi!r}")
    response = _csr(gene.X).tocsc()
    has_assign = ("guide_assignment" in guide.layers) and not prefer_counts
    grna = _csr(guide.layers["guide_assignment"] if has_assign else guide.X).tocsc()
    gv = guide.var.copy()
    if target_col not in gv.columns:
        if target_col == "intended_target_key" and "intended_target_name" in gv.columns:
            target_col = "intended_target_name"
        else:
            raise SystemExit(f"guide var lacks '{target_col}' (columns: {list(gv.columns)})")
    grna_ids = [str(x) for x in guide.var_names]
    targets = pd.Series([str(x) for x in gv[target_col].values], index=grna_ids, name="grna_target")
    extra, dropped = prepare_extra_covariates(mods.obs.reindex(gene.obs_names) if mods.obs is not None else None,
                                              remove_collinear_covariates, collinear_mode)
    names = _gene_names(gene)
    cov = compute_covariates(response, grna, names)
    cov.index = gene.obs_names.astype(str)
    for c in extra.columns:
        cov[c] = extra[c].values
    data = SceptreData(response, [str(x) for x in gene.var_names], names, grna, grna_ids, gv, targets, moi, cov,
                       list(extra.columns), list(gene.obs_names.astype(str)), has_assign)
    return data, dropped


def compute_covariates(response, grna, response_names):
    np, pd = _np(), _pd()
    r = response.tocsr()
    g = grna.tocsr()
    cov = pd.DataFrame({
        "response_n_nonzero": np.asarray((r != 0).sum(1)).ravel().astype(float),
        "response_n_umis": np.asarray(r.sum(1)).ravel(),
    })
    mt = np.array([bool(re.match(r"^MT-", n, flags=re.I)) for n in response_names])
    if mt.any():
        mt_umis = np.asarray(r[:, np.where(mt)[0]].sum(1)).ravel()
        with np.errstate(divide="ignore", invalid="ignore"):
            cov["response_p_mito"] = np.where(cov["response_n_umis"] > 0, mt_umis / cov["response_n_umis"], 0.0)
    cov["grna_n_nonzero"] = np.asarray((g != 0).sum(1)).ravel().astype(float)
    cov["grna_n_umis"] = np.asarray(g.sum(1)).ravel()
    return cov


def auto_formula(cov, low_moi: bool, include_grna_covariates: Optional[bool] = None,
                 extra: Optional["list[str]"] = None) -> str:
    """sceptre:::auto_construct_formula_object."""
    if include_grna_covariates is None:
        include_grna_covariates = not low_moi
    terms = []
    counts = ["response_n_nonzero", "response_n_umis"] + (["grna_n_nonzero", "grna_n_umis"] if include_grna_covariates else [])
    for c in counts:
        if c in cov.columns and cov[c].nunique() > 1:
            terms.append(f"log({c} + 1)" if (cov[c] == 0).any() else f"log({c})")
    if "response_p_mito" in cov.columns and cov["response_p_mito"].nunique() > 1:
        terms.append("response_p_mito")
    for c in (extra or []):
        if c in cov.columns and cov[c].nunique() > 1:
            terms.append(c if re.match(r"^[A-Za-z_][A-Za-z0-9_.]*$", c) and "." not in c else f'Q("{c}")')
    return "~ " + (" + ".join(terms) if terms else "1")


def design_matrix(formula: str, cov):
    """patsy design with R-style log(); returns (Z ndarray, column names, finite-row mask)."""
    np = _np()
    import patsy  # type: ignore
    f = formula.strip()
    if f.lower() == "default":
        f = DEFAULT_INFERENCE_FORMULA
    if not f.startswith("~"):
        f = "~ " + f
    with np.errstate(divide="ignore", invalid="ignore"):
        env = patsy.EvalEnvironment([{"log": np.log, "log1p": np.log1p, "exp": np.exp}])
        dm = patsy.dmatrix(f, cov, return_type="dataframe", eval_env=env, NA_action=patsy.NAAction(NA_types=[]))
    Z = dm.values.astype(float)
    ok = np.isfinite(Z).all(axis=1)
    Z[~ok] = 0.0
    keep = _independent_columns(Z[ok])
    return Z[:, keep], [dm.columns[i] for i in keep], ok


def _independent_columns(Z) -> "list[int]":
    np = _np()
    keep: "list[int]" = []
    for j in range(Z.shape[1]):
        cand = keep + [j]
        if np.linalg.matrix_rank(Z[:, cand]) == len(cand):
            keep.append(j)
    return keep


# ---------------------------------------------------------------------------
# gRNA assignment
# ---------------------------------------------------------------------------

def assign_thresholding(grna, threshold: float):
    sp, np = _sp(), _np()
    g = grna.tocsc().copy()
    g.data = (g.data >= threshold).astype(np.float64)
    g.eliminate_zeros()
    return g.astype(bool).tocsc()


def assign_maximum(grna, umi_fraction_threshold: float = 0.5, min_grna_n_umis_threshold: float = 5):
    sp, np = _sp(), _np()
    g = grna.tocsr()
    n, k = g.shape
    tot = np.asarray(g.sum(1)).ravel()
    rows, cols = [], []
    for i in range(n):
        s, e = g.indptr[i], g.indptr[i + 1]
        if e == s or tot[i] < min_grna_n_umis_threshold:
            continue
        d = g.data[s:e]
        j = int(np.argmax(d))
        if d[j] / tot[i] >= umi_fraction_threshold:
            rows.append(i)
            cols.append(g.indices[s + j])
    m = sp.csr_matrix((np.ones(len(rows), dtype=bool), (rows, cols)), shape=(n, k))
    return m.tocsc()


def _poisson_background(y, Z):
    np = _np()
    import statsmodels.api as sm  # type: ignore
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = sm.GLM(y, Z, family=sm.families.Poisson()).fit(maxiter=50)
        mu = np.asarray(fit.fittedvalues, dtype=float)
    except Exception:
        mu = np.full(len(y), max(y.mean(), 1e-8))
    return np.clip(mu, 1e-10, None)


def em_mixture(y, mu, n_em_rep: int = 5, pi_guess_range=(1e-5, 0.1), effect_size_guess_range=(0.5, 10.0),
               max_iter: int = 100, tol: float = 1e-6, seed: int = 4):
    """Reduced EM for y ~ pi Pois(f mu) + (1 - pi) Pois(mu); returns (posterior, pi, log f, loglik)."""
    np = _np()
    rng = np.random.default_rng(seed)
    logmu = np.log(mu)
    best = None
    for _ in range(n_em_rep):
        pi = rng.uniform(*pi_guess_range)
        b = rng.uniform(*effect_size_guess_range)
        prev = -np.inf
        for _it in range(max_iter):
            l1 = y * (logmu + b) - mu * np.exp(b)
            l0 = y * logmu - mu
            a1 = np.log(pi) + l1
            a0 = np.log1p(-pi) + l0
            m = np.maximum(a1, a0)
            lse = m + np.log(np.exp(a1 - m) + np.exp(a0 - m))
            T = np.exp(a1 - lse)
            ll = float(lse.sum())
            pi = float(np.clip(T.mean(), 1e-8, 1 - 1e-8))
            num, den = float((T * y).sum()), float((T * mu).sum())
            b = math.log(max(num, 1e-12) / max(den, 1e-12))
            b = max(b, 0.0)
            if abs(ll - prev) < tol * max(1.0, abs(ll)):
                break
            prev = ll
        if best is None or ll > best[3]:
            best = (T, pi, b, ll)
    return best


def assign_mixture(grna, Z, probability_threshold: float = 0.8, n_em_rep: int = 5, n_nonzero_cells_cutoff: int = 10,
                   backup_threshold: float = 5, seed: int = 4):
    sp, np = _sp(), _np()
    g = grna.tocsc()
    n, k = g.shape
    rows, cols = [], []
    info = []
    for j in range(k):
        y = np.asarray(g[:, j].todense()).ravel()
        nz = int((y > 0).sum())
        if nz < n_nonzero_cells_cutoff:
            idx = np.where(y >= backup_threshold)[0]
            method = "backup_thresholding"
            pi_hat, b_hat = float("nan"), float("nan")
        else:
            mu = _poisson_background(y, Z)
            T, pi_hat, b_hat, _ll = em_mixture(y, mu, n_em_rep=n_em_rep, seed=seed + j)
            idx = np.where((T >= probability_threshold) & (y > 0))[0]
            method = "mixture"
        rows += list(idx)
        cols += [j] * len(idx)
        info.append({"grna_index": j, "method": method, "n_nonzero": nz, "n_assigned": len(idx), "pi": pi_hat,
                     "log_fold": b_hat})
    m = sp.csr_matrix((np.ones(len(rows), dtype=bool), (rows, cols)), shape=(n, k)).tocsc()
    return m, info


def assign_grnas(data: SceptreData, method: str = "default", threshold: Optional[float] = None,
                 probability_threshold: float = 0.8, n_em_rep: int = 5, formula: Optional[str] = None,
                 umi_fraction_threshold: float = 0.5, min_grna_n_umis_threshold: float = 5):
    if method == "default":
        method = "maximum" if data.low_moi else "mixture"
    info: list = []
    if method == "thresholding":
        A = assign_thresholding(data.grna, 5 if threshold is None else threshold)
    elif method == "maximum":
        A = assign_maximum(data.grna, umi_fraction_threshold, min_grna_n_umis_threshold)
    elif method == "mixture":
        f = formula or auto_formula(data.covariates, data.low_moi, include_grna_covariates=False,
                                    extra=data.extra_covariate_names)
        Z, _, ok = design_matrix(f, data.covariates)
        A, info = assign_mixture(data.grna, Z, probability_threshold=probability_threshold, n_em_rep=n_em_rep)
    else:
        raise SystemExit(f"unknown assignment method {method}")
    return A, method, info


# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------

def run_qc(data: SceptreData, assignment, response_n_umis_range=(0.01, 0.99), response_n_nonzero_range=(0.01, 0.99),
           p_mito_threshold: float = 0.2, remove_cells_w_zero_or_twoplus_grnas: Optional[bool] = None,
           additional_cells_to_remove: Optional["list[int]"] = None):
    np = _np()
    cov = data.covariates
    keep = np.ones(data.n_cells, dtype=bool)
    reasons = {}
    for col, rng_ in (("response_n_umis", response_n_umis_range), ("response_n_nonzero", response_n_nonzero_range)):
        v = cov[col].values.astype(float)
        lo, hi = np.quantile(v, rng_[0]), np.quantile(v, rng_[1])
        bad = (v < lo) | (v > hi)
        reasons[f"{col}_range"] = int(bad.sum())
        keep &= ~bad
    if "response_p_mito" in cov.columns:
        bad = cov["response_p_mito"].values > p_mito_threshold
        reasons["response_p_mito"] = int(bad.sum())
        keep &= ~bad
    if remove_cells_w_zero_or_twoplus_grnas is None:
        remove_cells_w_zero_or_twoplus_grnas = data.low_moi
    n_assigned = np.asarray(assignment.sum(1)).ravel()
    if remove_cells_w_zero_or_twoplus_grnas:
        bad = n_assigned != 1
        reasons["zero_or_twoplus_grnas"] = int(bad.sum())
        keep &= ~bad
    if additional_cells_to_remove:
        keep[list(additional_cells_to_remove)] = False
        reasons["additional"] = len(additional_cells_to_remove)
    return keep, reasons


# ---------------------------------------------------------------------------
# SCEPTRE test engine
# ---------------------------------------------------------------------------

def nb_loglik_theta(y, mu, theta):
    from scipy.special import gammaln  # type: ignore
    np = _np()
    return float(np.sum(gammaln(y + theta) - gammaln(theta) - gammaln(y + 1) + theta * np.log(theta / (theta + mu))
                        + y * np.log(mu / (theta + mu))))


def estimate_theta(y, mu) -> float:
    from scipy.optimize import minimize_scalar  # type: ignore
    np = _np()
    r = minimize_scalar(lambda lt: -nb_loglik_theta(y, mu, math.exp(lt)), bounds=(-6.0, 9.0), method="bounded",
                        options={"xatol": 1e-4})
    return float(math.exp(r.x)) if np.isfinite(r.fun) else 100.0


def fit_nb_glm(y, Z):
    """Poisson GLM -> theta ML -> NB GLM with that theta. Returns (beta, theta) or None."""
    np = _np()
    import statsmodels.api as sm  # type: ignore
    if y.sum() == 0:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            pfit = sm.GLM(y, Z, family=sm.families.Poisson()).fit(maxiter=50)
            mu0 = np.clip(np.asarray(pfit.fittedvalues), 1e-10, None)
            theta = estimate_theta(y, mu0)
            nfit = sm.GLM(y, Z, family=sm.families.NegativeBinomial(alpha=1.0 / theta)).fit(
                start_params=pfit.params, maxiter=50)
            beta = np.asarray(nfit.params, dtype=float)
            if not np.all(np.isfinite(beta)):
                beta = np.asarray(pfit.params, dtype=float)
        except Exception:
            return None
    return beta, theta


def skewnorm_mom(x):
    """Method-of-moments skew-normal fit: returns (xi, omega, alpha)."""
    np = _np()
    from scipy import stats  # type: ignore
    m, s = float(np.mean(x)), float(np.std(x, ddof=1))
    g = float(stats.skew(x))
    gmax = 0.99
    g = max(min(g, gmax), -gmax)
    a = abs(g) ** (2.0 / 3.0)
    delta = math.copysign(math.sqrt((math.pi / 2) * a / (a + ((4 - math.pi) / 2) ** (2.0 / 3.0))), g) if g != 0 else 0.0
    delta = max(min(delta, 0.995), -0.995)
    alpha = delta / math.sqrt(1 - delta ** 2)
    omega = s / math.sqrt(1 - 2 * delta ** 2 / math.pi)
    xi = m - omega * delta * math.sqrt(2 / math.pi)
    return xi, omega, alpha


def empirical_p(z: float, null, side: str) -> float:
    np = _np()
    B = len(null)
    left = (1 + np.sum(null <= z)) / (1 + B)
    right = (1 + np.sum(null >= z)) / (1 + B)
    if side == "left":
        return float(left)
    if side == "right":
        return float(right)
    return float(min(1.0, 2 * min(left, right)))


def skewnorm_p(z: float, params, side: str) -> float:
    from scipy import stats  # type: ignore
    xi, omega, alpha = params
    left = float(stats.skewnorm.cdf(z, alpha, loc=xi, scale=omega))
    right = float(stats.skewnorm.sf(z, alpha, loc=xi, scale=omega))
    if side == "left":
        return left
    if side == "right":
        return right
    return min(1.0, 2 * min(left, right))


def skewnorm_fit_ok(null, params, max_ks: float = 0.03) -> bool:
    from scipy import stats  # type: ignore
    xi, omega, alpha = params
    if not (math.isfinite(xi) and math.isfinite(omega) and omega > 0):
        return False
    ks = stats.kstest(null, lambda v: stats.skewnorm.cdf(v, alpha, loc=xi, scale=omega)).statistic
    return ks <= max_ks


class Tester:
    """Runs SCEPTRE pair tests over the cells that passed QC."""

    def __init__(self, data: SceptreData, assignment, cells_in_use, formula: str, side: str = "both",
                 control_group: str = "default", resampling_mechanism: str = "default",
                 resampling_approximation: str = "skew_normal", B1: int = 499, B2: int = 4999, B3: int = 24999,
                 p_thresh: float = 0.02, n_nonzero_trt_thresh: int = 7, n_nonzero_cntrl_thresh: int = 7,
                 seed: int = 4):
        np = _np()
        self.data = data
        self.used = np.where(cells_in_use)[0]
        self.A = assignment.tocsc()[self.used, :]
        cov = data.covariates.iloc[self.used].reset_index(drop=True)
        Z, cols, ok = design_matrix(formula, cov)
        if not ok.all():
            log.warning("%d cells with non-finite covariates dropped from testing", int((~ok).sum()))
            self.used = self.used[ok]
            self.A = self.A[np.where(ok)[0], :]
            Z = Z[ok]
        self.Z, self.Z_cols, self.formula = Z, cols, formula
        self.Y = data.response.tocsc()[self.used, :]
        self.side = side
        self.control_group = control_group if control_group != "default" else ("nt_cells" if data.low_moi else "complement")
        self.mechanism = resampling_mechanism if resampling_mechanism != "default" else ("permutations" if data.low_moi else "crt")
        self.approximation = resampling_approximation
        self.B1, self.B2, self.B3, self.p_thresh = B1, B2, B3, p_thresh
        self.trt_thresh, self.cntrl_thresh = n_nonzero_trt_thresh, n_nonzero_cntrl_thresh
        self.seed = seed
        self.nt_idx = [j for j in range(len(data.grna_ids)) if data.grna_is_nt[j]]
        self.nt_mask = (np.asarray(self.A[:, self.nt_idx].sum(1)).ravel() > 0) if self.nt_idx else np.zeros(len(self.used), bool)
        self._fit_cache: "OrderedDict" = OrderedDict()
        self._null_cache: "OrderedDict" = OrderedDict()
        self._pi_cache: "OrderedDict" = OrderedDict()
        self.target_to_idx: "dict[str, list[int]]" = {}
        for j, t in enumerate(data.grna_target.values):
            self.target_to_idx.setdefault(t, []).append(j)
        self.gene_pos = {g: i for i, g in enumerate(data.response_ids)}
        for i, nm in enumerate(data.response_names):
            self.gene_pos.setdefault(nm, i)

    # -- indicators -------------------------------------------------------
    def indicator(self, grna_idx: "list[int]"):
        np = _np()
        return np.asarray(self.A[:, grna_idx].sum(1)).ravel() > 0

    def control_mask(self, x, exclude_idx: "Optional[list[int]]" = None):
        np = _np()
        if self.control_group == "complement":
            return ~x
        nt = self.nt_idx if exclude_idx is None else [j for j in self.nt_idx if j not in set(exclude_idx)]
        m = (np.asarray(self.A[:, nt].sum(1)).ravel() > 0) if nt else np.zeros(len(x), bool)
        return m & ~x

    # -- response model ---------------------------------------------------
    def response_fit(self, gi: int):
        np = _np()
        key = gi
        if key in self._fit_cache:
            self._fit_cache.move_to_end(key)
            return self._fit_cache[key]
        y = np.asarray(self.Y[:, gi].todense()).ravel()
        fit_rows = np.arange(len(y)) if self.control_group == "complement" else np.where(self.nt_mask)[0]
        res = fit_nb_glm(y[fit_rows], self.Z[fit_rows]) if len(fit_rows) > self.Z.shape[1] else None
        out = None
        if res is not None:
            beta, theta = res
            mu = np.exp(np.clip(self.Z @ beta, -30, 30))
            out = (y, mu, theta)
        self._fit_cache[key] = out
        if len(self._fit_cache) > 256:
            self._fit_cache.popitem(last=False)
        return out

    # -- null resamples ---------------------------------------------------
    def _crt_pi(self, key, x, cells):
        np = _np()
        if key in self._pi_cache:
            return self._pi_cache[key]
        import statsmodels.api as sm  # type: ignore
        xa = x[cells].astype(float)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = sm.GLM(xa, self.Z[cells], family=sm.families.Binomial()).fit(maxiter=50)
            pi = np.clip(np.asarray(fit.fittedvalues), 1e-6, 1 - 1e-6)
        except Exception:
            pi = np.full(len(cells), max(xa.mean(), 1e-6))
        self._pi_cache[key] = pi
        if len(self._pi_cache) > 64:
            self._pi_cache.popitem(last=False)
        return pi

    def null_matrix(self, key, x, cells, B: int, stage: int):
        """B x n_A CSR matrix of resampled treatment indicators (analysis-set coordinates)."""
        sp, np = _sp(), _np()
        ck = (key, stage)
        if ck in self._null_cache:
            return self._null_cache[ck]
        rng = np.random.default_rng(_seed_for(self.seed, key, stage))
        nA = len(cells)
        if self.mechanism == "permutations":
            k = int(x[cells].sum())
            idx = np.empty((B, k), dtype=np.int64)
            for b in range(B):
                idx[b] = rng.choice(nA, size=k, replace=False)
            M = sp.csr_matrix((np.ones(B * k), idx.ravel(), np.arange(0, B * k + 1, k) if k else np.zeros(B + 1, int)),
                              shape=(B, nA))
        else:
            pi = self._crt_pi(key, x, cells)
            blocks = []
            chunk = max(1, int(2_000_000 // max(nA, 1)))
            for s in range(0, B, chunk):
                b = min(chunk, B - s)
                blocks.append(sp.csr_matrix(rng.random((b, nA)) < pi[None, :]).astype(np.float64))
            M = sp.vstack(blocks).tocsr()
        if stage == 1:
            self._null_cache[ck] = M
            if len(self._null_cache) > 32:
                self._null_cache.popitem(last=False)
        return M

    # -- one pair ---------------------------------------------------------
    def test(self, key: str, grna_idx: "list[int]", gi: int, control_exclude: "Optional[list[int]]" = None) -> dict:
        np = _np()
        x = self.indicator(grna_idx)
        ctrl = self.control_mask(x, control_exclude)
        fit = self.response_fit(gi)
        y = fit[0] if fit is not None else np.asarray(self.Y[:, gi].todense()).ravel()
        n_nz_trt = int(((y > 0) & x).sum())
        n_nz_ctl = int(((y > 0) & ctrl).sum())
        out = {"n_nonzero_trt": n_nz_trt, "n_nonzero_cntrl": n_nz_ctl,
               "pass_qc": bool(n_nz_trt >= self.trt_thresh and n_nz_ctl >= self.cntrl_thresh and x.sum() > 0 and ctrl.sum() > 0),
               "p_value": float("nan"), "log_2_fold_change": float("nan"), "z": float("nan"), "p_method": ""}
        if not out["pass_qc"] or fit is None:
            return out
        _, mu, theta = fit
        cells = np.where(x | ctrl)[0]
        ya, mua, xa = y[cells], mu[cells], x[cells].astype(float)
        Za = self.Z[cells]
        denom = 1 + mua / theta
        r = (ya - mua) / denom
        w = mua / denom
        ZtWZ = Za.T @ (w[:, None] * Za)
        try:
            L = np.linalg.cholesky(ZtWZ + 1e-10 * np.eye(ZtWZ.shape[0]))
        except np.linalg.LinAlgError:
            return out
        V = np.linalg.solve(L, (w[:, None] * Za).T).T  # n_A x p

        def stats_for(M):
            U = M @ r
            I = M @ w - np.sum(np.asarray(M @ V) ** 2, axis=1)
            return U / np.sqrt(np.clip(I, 1e-12, None))

        z_obs = float(stats_for(xa[None, :])[0])
        out["z"] = z_obs
        out["log_2_fold_change"] = log2_fold_change(ya[xa > 0], mua[xa > 0], theta)
        nkey = key
        null1 = stats_for(self.null_matrix(nkey, x, cells, self.B1, 1))
        p = empirical_p(z_obs, null1, self.side)
        method = f"empirical_B{self.B1}"
        if p <= self.p_thresh:
            done = False
            if self.approximation == "skew_normal":
                null2 = stats_for(self.null_matrix(nkey, x, cells, self.B2, 2))
                prm = skewnorm_mom(null2)
                if skewnorm_fit_ok(null2, prm):
                    p = skewnorm_p(z_obs, prm, self.side)
                    method = "skew_normal"
                    done = True
            if not done:
                null3 = stats_for(self.null_matrix(nkey, x, cells, self.B3, 3))
                p = empirical_p(z_obs, null3, self.side)
                method = f"empirical_B{self.B3}"
        out["p_value"] = float(p)
        out["p_method"] = method
        return out


def log2_fold_change(y, mu, theta) -> float:
    """MLE of b in y ~ NB(mu e^b, theta) over treated cells, in log2 units."""
    np = _np()
    from scipy.optimize import brentq  # type: ignore
    if len(y) == 0:
        return float("nan")

    def g(b):
        m = mu * math.exp(b)
        return float(np.sum((y - m) / (1 + m / theta)))

    lo, hi = -10.0, 10.0
    if g(lo) <= 0:
        return lo / math.log(2)
    if g(hi) >= 0:
        return hi / math.log(2)
    return brentq(g, lo, hi, xtol=1e-8) / math.log(2)


# ---------------------------------------------------------------------------
# Pair sets and analyses
# ---------------------------------------------------------------------------

def construct_trans_pairs(data: SceptreData):
    pd = _pd()
    targets = sorted({t for t, nt in zip(data.grna_target.values, data.grna_is_nt) if not nt})
    return pd.DataFrame([(t, g) for t in targets for g in data.response_ids], columns=["grna_target", "response_id"])


def construct_cis_pairs(data: SceptreData, gene_var, distance_threshold: int = 500_000):
    """Targets vs genes whose TSS (gene start; end on '-' strand) lies within distance_threshold of the target."""
    np, pd = _np(), _pd()
    gv = data.grna_var
    need = ["intended_target_chr", "intended_target_start", "intended_target_end"]
    if not all(c in gv.columns for c in need):
        raise SystemExit("cis pairs need intended_target_chr/start/end in guide var")
    chr_col = next((c for c in ("chr", "chromosome", "gene_chr", "seqname") if c in gene_var.columns), None)
    start_col = next((c for c in ("start", "gene_start", "tss") if c in gene_var.columns), None)
    if chr_col is None or start_col is None:
        raise SystemExit("cis pairs need gene var chr + start (or tss) columns")
    tss = gene_var[start_col].astype(float).values
    if "strand" in gene_var.columns and "end" in gene_var.columns:
        tss = np.where(gene_var["strand"].astype(str).values == "-", gene_var["end"].astype(float).values, tss)
    gchr = gene_var[chr_col].astype(str).values
    rows = []
    tgt = pd.DataFrame({"t": data.grna_target.values, "chr": gv["intended_target_chr"].astype(str).values,
                        "s": pd.to_numeric(gv["intended_target_start"], errors="coerce").values,
                        "e": pd.to_numeric(gv["intended_target_end"], errors="coerce").values, "nt": data.grna_is_nt})
    tgt = tgt[~tgt["nt"]].groupby("t").agg(chr=("chr", "first"), s=("s", "min"), e=("e", "max")).reset_index()
    for _, t in tgt.iterrows():
        if not np.isfinite(t["s"]):
            continue
        d = np.where(tss < t["s"], t["s"] - tss, np.where(tss > t["e"], tss - t["e"], 0))
        hit = np.where((gchr == t["chr"]) & (d <= distance_threshold))[0]
        rows += [(t["t"], data.response_ids[i]) for i in hit]
    return pd.DataFrame(rows, columns=["grna_target", "response_id"])


def groups_for(tester: Tester, strategy: str, target: str) -> "list[tuple[str, list[int]]]":
    idx = tester.target_to_idx.get(target, [])
    if strategy == "union":
        return [(target, idx)]
    return [(tester.data.grna_ids[j], [j]) for j in idx]


def run_pairs(tester: Tester, pairs, strategy: str = "union", alpha: float = 0.1, correct: bool = True,
              control_exclude_map: "Optional[dict]" = None, group_idx_map: "Optional[dict]" = None):
    """Test grna_target x response_id pairs; returns a sceptre-schema results DataFrame."""
    pd, np = _pd(), _np()
    rows = []
    pairs = pairs.reset_index(drop=True)
    order = pairs.sort_values(["grna_target", "response_id"]).index
    for i in order:
        t, gname = str(pairs.at[i, "grna_target"]), str(pairs.at[i, "response_id"])
        gi = tester.gene_pos.get(gname)
        if group_idx_map is not None and t in group_idx_map:
            groups = [(t, group_idx_map[t])]
        else:
            groups = groups_for(tester, "singleton" if strategy in ("singleton", "bonferroni") else "union", t)
        if gi is None or not groups or not any(g[1] for g in groups):
            base = {"response_id": gname, "grna_target": t, "n_nonzero_trt": 0, "n_nonzero_cntrl": 0, "pass_qc": False,
                    "p_value": np.nan, "log_2_fold_change": np.nan, "z": np.nan, "p_method": "missing"}
            if strategy == "singleton":
                base["grna_id"] = ""
            rows.append(base)
            continue
        excl = control_exclude_map.get(t) if control_exclude_map else None
        res = [(gid, tester.test(gid, idx, gi, excl)) for gid, idx in groups]
        if strategy == "singleton":
            for gid, r in res:
                rows.append({"response_id": gname, "grna_id": gid, "grna_target": t, **r})
        elif strategy == "bonferroni":
            ps = [r["p_value"] for _, r in res if np.isfinite(r["p_value"])]
            u = tester.test(t, tester.target_to_idx.get(t, []), gi, excl)
            if ps:
                u["p_value"] = float(min(1.0, len(res) * min(ps)))
                u["p_method"] = "bonferroni"
            else:
                u["p_value"] = np.nan
            rows.append({"response_id": gname, "grna_target": t, **u})
        else:
            rows.append({"response_id": gname, "grna_target": t, **res[0][1]})
    df = pd.DataFrame(rows)
    if df.empty:
        cols = RESULT_COLS + (["grna_id"] if strategy == "singleton" else [])
        return pd.DataFrame(columns=cols + ["z", "p_method"])
    if correct:
        q = bh_adjust(df["p_value"].where(df["pass_qc"]).values)
        df["significant"] = np.where(np.isfinite(q), q < alpha, False)
    lead = ["response_id"] + (["grna_id"] if "grna_id" in df.columns else []) + ["grna_target"]
    rest = [c for c in RESULT_COLS if c not in lead and c in df.columns]
    return df[lead + rest + [c for c in df.columns if c not in lead + rest]]


def calibration_pairs(tester: Tester, n_pairs: int, group_size: int, seed: int = 4):
    """Undercover NT groups x random responses (sceptre run_calibration_check)."""
    np, pd = _np(), _pd()
    rng = np.random.default_rng(seed)
    nt = list(tester.nt_idx)
    if not nt:
        raise SystemExit("calibration check needs non-targeting gRNAs")
    rng.shuffle(nt)
    group_size = max(1, min(group_size, len(nt) - (1 if tester.control_group == "nt_cells" else 0)))
    groups = {}
    for s in range(0, len(nt) - group_size + 1, group_size):
        idx = nt[s:s + group_size]
        groups["&".join(tester.data.grna_ids[j] for j in idx)] = idx
    names = list(groups)
    genes = list(tester.data.response_ids)
    pairs = []
    seen = set()
    tries = 0
    while len(pairs) < n_pairs and tries < n_pairs * 50:
        tries += 1
        g = names[int(rng.integers(len(names)))]
        r = genes[int(rng.integers(len(genes)))]
        if (g, r) in seen:
            continue
        seen.add((g, r))
        pairs.append((g, r))
    return pd.DataFrame(pairs, columns=["grna_target", "response_id"]), groups


def run_calibration_check(tester: Tester, n_calibration_pairs: int, alpha: float = 0.1, group_size: Optional[int] = None,
                          seed: int = 4):
    np = _np()
    if group_size is None:
        sizes = [len(v) for t, v in tester.target_to_idx.items() if not is_nt_target(t)
                 and not all(tester.data.grna_is_nt[j] for j in v)]
        group_size = int(np.median(sizes)) if sizes else 1
    pairs, groups = calibration_pairs(tester, n_calibration_pairs, group_size, seed=seed)
    excl = {k: v for k, v in groups.items()}
    res = run_pairs(tester, pairs, "union", alpha=alpha, control_exclude_map=excl, group_idx_map=groups)
    ok = res["pass_qc"] & res["p_value"].notna()
    summ = {"n_calibration_pairs": int(len(res)), "n_passing_qc": int(ok.sum()),
            "n_false_discoveries": int(res["significant"].sum()), "alpha": alpha,
            "median_p": float(res.loc[ok, "p_value"].median()) if ok.any() else float("nan"),
            "mean_log2_fc": float(res.loc[ok, "log_2_fold_change"].mean()) if ok.any() else float("nan"),
            "undercover_group_size": group_size}
    if ok.sum() >= 5:
        from scipy import stats  # type: ignore
        summ["ks_uniform_p"] = float(stats.kstest(res.loc[ok, "p_value"].values, "uniform").pvalue)
    return res, summ


def pairs_from_uns(mods: Modalities, data: SceptreData, target_col: str = "intended_target_name"):
    """pairs_to_test (any schema) -> DataFrame with grna_target, response_id (+ original columns)."""
    pd = _pd()
    p = uns_table(mods.uns, "pairs_to_test")
    if p is None:
        return None
    p = p.copy()
    if "gene_id" not in p.columns and "gene_name" in p.columns:
        p = p.rename(columns={"gene_name": "gene_id"})
    if target_col not in p.columns:
        if "intended_target_name" in p.columns:
            key = dict(zip(data.grna_var["intended_target_name"].astype(str), data.grna_target.values))
            p[target_col] = p["intended_target_name"].astype(str).map(key)
        elif "guide_id" in p.columns:
            gid = data.grna_var["guide_id"].astype(str) if "guide_id" in data.grna_var.columns else pd.Series(data.grna_ids)
            key = dict(zip(gid.values, data.grna_target.values))
            p[target_col] = p["guide_id"].astype(str).map(key)
        else:
            raise SystemExit("pairs_to_test needs intended_target_name (or intended_target_key / guide_id) and gene_id")
    p["grna_target"] = p[target_col].astype(str)
    p["response_id"] = p["gene_id"].astype(str)
    return p


def read_pairs_tsv(path: str):
    pd = _pd()
    p = pd.read_csv(path, sep=None, engine="python")
    ren = {"intended_target_name": "grna_target", "gene_id": "response_id"}
    p = p.rename(columns={k: v for k, v in ren.items() if k in p.columns and v not in p.columns})
    if not {"grna_target", "response_id"} <= set(p.columns):
        raise SystemExit(f"{path}: need grna_target/response_id (or intended_target_name/gene_id)")
    return p


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def qq_figure(d: Path, series: "dict[str, Any]", title: str, name: str = "qq_plot.png"):
    plt = _plt()
    np = _np()
    if plt is None:
        return None
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    colors = [BLUE, ORANGE, "#1baf7a", "#4a3aa7"]
    mx = 1.0
    for (lab, p), c in zip(series.items(), colors):
        p = np.sort(np.asarray(p, float)[np.isfinite(p)])
        if len(p) == 0:
            continue
        e = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
        o = -np.log10(np.clip(p, 1e-300, 1))
        mx = max(mx, e.max())
        ax.scatter(e, o, s=10, color=c, label=f"{lab} (n={len(p)})", alpha=0.8, edgecolors="none")
    ax.plot([0, mx], [0, mx], color=INK2, lw=0.8, ls="--")
    ax.legend(fontsize=7, frameon=False)
    _style(ax, title, "expected -log10 p", "observed -log10 p")
    return _save(fig, d / name)


def volcano_figure(d: Path, res, title: str):
    plt = _plt()
    np = _np()
    if plt is None or res.empty:
        return None
    ok = res["p_value"].notna()
    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    sig = res["significant"].fillna(False).astype(bool) if "significant" in res.columns else ok & False
    ax.scatter(res.loc[ok & ~sig, "log_2_fold_change"], -np.log10(res.loc[ok & ~sig, "p_value"].clip(1e-300)), s=10,
               color=AXIS, label="not significant", edgecolors="none")
    ax.scatter(res.loc[ok & sig, "log_2_fold_change"], -np.log10(res.loc[ok & sig, "p_value"].clip(1e-300)), s=12,
               color=ORANGE, label="significant (BH q < 0.1)", edgecolors="none")
    ax.legend(fontsize=7, frameon=False)
    _style(ax, title, "log2 fold change", "-log10 p")
    return _save(fig, d / "volcano.png")


def assignment_figure(d: Path, A):
    plt = _plt()
    np = _np()
    if plt is None:
        return None
    per_cell = np.asarray(A.sum(1)).ravel()
    per_grna = np.asarray(A.sum(0)).ravel()
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.3))
    axes[0].hist(per_cell, bins=np.arange(per_cell.max() + 2) - 0.5, color=BLUE)
    _style(axes[0], "gRNAs per cell", "assigned gRNAs", "cells")
    axes[1].hist(per_grna, bins=30, color=ORANGE)
    _style(axes[1], "cells per gRNA", "assigned cells", "gRNAs")
    return _save(fig, d / "guide_assignment.png")


# ---------------------------------------------------------------------------
# Shared workflow pieces
# ---------------------------------------------------------------------------

def _load(args, remove_collinear: bool, collinear_mode: str = "all-or-nothing", target_col: str = "intended_target_name",
          prefer_counts: bool = False):
    mods = read_h5mu(Path(args.mudata))
    data, dropped = convert_mudata_to_sceptre(mods, remove_collinear, collinear_mode, target_col,
                                              moi=getattr(args, "moi", None), prefer_counts=prefer_counts)
    return mods, data, dropped


def _tester_kwargs(args) -> dict:
    return dict(side=args.side, control_group=args.control_group, resampling_mechanism=args.resampling_mechanism,
                resampling_approximation=args.resampling_approximation, B1=args.b1, B2=args.b2, B3=args.b3,
                seed=args.seed)


def _result_md(res, n: int = 25) -> str:
    cols = [c for c in ["grna_target", "grna_id", "response_id", "n_nonzero_trt", "n_nonzero_cntrl", "pass_qc", "p_value",
                        "log_2_fold_change", "significant"] if c in res.columns]
    rows = res.sort_values("p_value").head(n)
    return md_table(cols, [[_fmt(r[c]) if c in ("p_value", "log_2_fold_change") else r[c] for c in cols]
                           for _, r in rows.iterrows()])


def sceptre_to_mudata(data: SceptreData, cells_in_use, assignment, positive_control_pairs, discovery_pairs,
                      power_result, discovery_result, gene_var=None):
    """sceptre_object_to_mudata: returns (gene AnnData, guide AnnData, obs, uns)."""
    import anndata as ad  # type: ignore
    np, pd = _np(), _pd()
    used = np.where(cells_in_use)[0]
    cell_ids = [data.cell_ids[i] for i in used]
    cov = data.covariates.iloc[used]
    sample_df = cov[[c for c in cov.columns if c not in COMPUTED_COVARIATES]].copy()
    sample_df.index = cell_ids
    gv = data.grna_var.copy()
    rowdata = pd.DataFrame(index=pd.Index(data.grna_ids, name=None))
    rowdata["targeting"] = np.where(data.grna_target.values != NT_LABEL, "TRUE", "FALSE")
    rowdata["intended_target_name"] = data.grna_target.values
    for src, dst, fill in (("intended_target_chr", "intended_target_chr", ""), ("chr", "intended_target_chr", ""),
                           ("intended_target_start", "intended_target_start", -9), ("start", "intended_target_start", -9),
                           ("intended_target_end", "intended_target_end", -9), ("end", "intended_target_end", -9)):
        if src in gv.columns and dst not in rowdata.columns:
            v = gv[src].values
            if fill == "":
                rowdata[dst] = ["" if (x is None or (isinstance(x, float) and math.isnan(x)) or str(x) == "nan") else str(x) for x in v]
            else:
                rowdata[dst] = pd.to_numeric(pd.Series(v), errors="coerce").fillna(fill).values
    gene_obs = pd.DataFrame({"num_expressed_genes": cov["response_n_nonzero"].values,
                             "total_gene_umis": cov["response_n_umis"].values}, index=cell_ids)
    guide_obs = pd.DataFrame({"num_expressed_guides": cov["grna_n_nonzero"].values,
                              "total_guide_umis": cov["grna_n_umis"].values}, index=cell_ids)
    gvar = gene_var.copy() if gene_var is not None else pd.DataFrame(index=data.response_ids)
    gvar.index = gvar.index.astype(str)
    gene = ad.AnnData(X=data.response.tocsr()[used, :], obs=gene_obs, var=gvar)
    guide = ad.AnnData(X=data.grna.tocsr()[used, :], obs=guide_obs, var=rowdata)
    guide.layers["guide_assignment"] = assignment.tocsr()[used, :].astype(np.float64)
    guide.uns["moi"] = data.moi
    pairs = pd.concat([positive_control_pairs[["grna_target", "response_id"]].assign(pair_type="positive_control"),
                       discovery_pairs[["grna_target", "response_id"]].assign(pair_type="discovery")], ignore_index=True)
    pt = pairs.rename(columns={"grna_target": "intended_target_name", "response_id": "gene_id"})
    both = pd.concat([r[["grna_target", "response_id", "p_value", "log_2_fold_change"]] for r in (power_result, discovery_result)
                      if r is not None and len(r)], ignore_index=True) if any(r is not None and len(r) for r in (power_result, discovery_result)) \
        else pd.DataFrame(columns=["grna_target", "response_id", "p_value", "log_2_fold_change"])
    tr = pairs.merge(both, on=["grna_target", "response_id"], how="left").rename(
        columns={"grna_target": "intended_target_name", "response_id": "gene_id", "log_2_fold_change": "log2_fc"})
    uns = {"pairs_to_test": pt, "test_results": tr}
    return gene, guide, sample_df, uns


def full_workflow(data: SceptreData, positive_pairs, discovery_pairs, formula: Optional[str], assign_method: str,
                  tkw: dict, qc_kw: Optional[dict] = None, strategy: str = "union", alpha: float = 0.1,
                  n_calibration_pairs: Optional[int] = None, run_calibration: bool = True, assign_kw: Optional[dict] = None):
    A, method, ainfo = assign_grnas(data, assign_method, **(assign_kw or {}))
    qc_kw = qc_kw or {}
    trt, ctl = qc_kw.pop("n_nonzero_trt_thresh", 7), qc_kw.pop("n_nonzero_cntrl_thresh", 7)
    keep, reasons = run_qc(data, A, **qc_kw)
    f = formula or auto_formula(data.covariates, data.low_moi, extra=data.extra_covariate_names)
    tester = Tester(data, A, keep, f, n_nonzero_trt_thresh=trt, n_nonzero_cntrl_thresh=ctl, **tkw)
    disc = run_pairs(tester, discovery_pairs, strategy, alpha=alpha)
    power = run_pairs(tester, positive_pairs, strategy, correct=False) if len(positive_pairs) else None
    if power is not None and "significant" in power.columns:
        power = power.drop(columns=["significant"])
    calib, csumm = (None, {})
    if run_calibration and tester.nt_idx:
        n_cal = n_calibration_pairs or max(1, int(disc["pass_qc"].sum()))
        calib, csumm = run_calibration_check(tester, n_cal, alpha=alpha, seed=tkw.get("seed", 4))
    return {"assignment": A, "assign_method": method, "assign_info": ainfo, "cells_in_use": keep, "qc_reasons": reasons,
            "formula": f, "tester": tester, "discovery": disc, "power": power, "calibration": calib,
            "calibration_summary": csumm}


def split_pairs(mods: Modalities, data: SceptreData, args):
    pd = _pd()
    pos = pd.DataFrame(columns=["grna_target", "response_id"])
    disc = None
    p = pairs_from_uns(mods, data)
    if p is not None:
        if "pair_type" in p.columns:
            pos = p[p["pair_type"].astype(str) == "positive_control"][["grna_target", "response_id"]]
            disc = p[p["pair_type"].astype(str) != "positive_control"][["grna_target", "response_id"]]
        else:
            disc = p[["grna_target", "response_id"]]
    pc_uns = uns_table(mods.uns, "positive_control_pairs")
    if pc_uns is not None and len(pos) == 0:
        pc_uns = pc_uns.rename(columns={"intended_target_name": "grna_target", "gene_id": "response_id"})
        pos = pc_uns[["grna_target", "response_id"]]
    if getattr(args, "positive_control_pairs", None):
        pos = read_pairs_tsv(args.positive_control_pairs)[["grna_target", "response_id"]]
    if getattr(args, "discovery_pairs", None):
        disc = read_pairs_tsv(args.discovery_pairs)[["grna_target", "response_id"]]
    pair_set = getattr(args, "pair_set", None)
    if pair_set == "trans" or (disc is None and pair_set != "cis"):
        disc = construct_trans_pairs(data)
    elif pair_set == "cis":
        disc = construct_cis_pairs(data, mods.gene.var, getattr(args, "cis_distance", 500_000))
    return pos.reset_index(drop=True), disc.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_convert(args) -> int:
    np, sp = _np(), _sp()
    mods, data, dropped = _load(args, args.remove_collinear_covariates, args.collinear_mode, args.target_col)
    d = run_dir(args.label)
    sp.save_npz(d / "response_matrix.npz", data.response.T.tocsc())
    sp.save_npz(d / "grna_matrix.npz", data.grna.T.tocsc())
    print(f"Wrote: {d / 'response_matrix.npz'}")
    print(f"Wrote: {d / 'grna_matrix.npz'}")
    pd = _pd()
    gtdf = pd.DataFrame({"grna_id": data.grna_ids, "grna_target": data.grna_target.values})
    write_tsv(gtdf, d / "grna_target_data_frame.tsv")
    cov = data.covariates.copy()
    cov.insert(0, "cell_id", data.cell_ids)
    write_tsv(cov, d / "covariate_data_frame.tsv")
    f = auto_formula(data.covariates, data.low_moi, extra=data.extra_covariate_names)
    (d / "formula.txt").write_text(f + "\n")
    print(f"Wrote: {d / 'formula.txt'}")
    summ = {"n_cells": data.n_cells, "n_responses": len(data.response_ids), "n_grnas": len(data.grna_ids),
            "n_nt_grnas": int(data.grna_is_nt.sum()), "n_targets": int(data.grna_target[~data.grna_is_nt].nunique()),
            "moi": data.moi, "grna_matrix_source": "guide_assignment" if data.grna_is_assignment else "counts",
            "extra_covariates": data.extra_covariate_names, "dropped_covariates": dropped, "auto_formula": f,
            "remove_collinear_covariates": args.remove_collinear_covariates, "collinear_mode": args.collinear_mode}
    write_json(summ, d / "summary.json")
    write_report(d, "sceptreIGVF: MuData -> sceptre object", [
        f"Source: `{args.mudata}`",
        md_table(["field", "value"], [[k, v] for k, v in summ.items()]),
        "Files: response_matrix.npz and grna_matrix.npz (features x cells, as sceptre stores them), "
        "grna_target_data_frame.tsv, covariate_data_frame.tsv, formula.txt."])
    return 0


def cmd_assign_guides(args) -> int:
    np, pd = _np(), _pd()
    mods, data, dropped = _load(args, True, args.collinear_mode, prefer_counts=True)
    d = run_dir(args.label)
    A, method, info = assign_grnas(data, args.method, threshold=args.threshold,
                                   probability_threshold=args.probability_threshold, n_em_rep=args.n_em_rep)
    guide = mods.guide.copy() if mods.guide.obs_names.equals(mods.gene.obs_names) else _align(mods).guide.copy()
    guide.layers["guide_assignment"] = A.tocsr().astype(np.float64)
    out = Path(args.out) if args.out else d / "guide_assignment_output.h5mu"
    write_h5mu(out, _align(mods).gene, guide, _align(mods).obs, mods.uns)
    try:
        from scipy.io import mmwrite  # type: ignore
        mmwrite(str(d / "guide_assignment.mtx"), A.T.astype(np.float64).tocoo())
        print(f"Wrote: {d / 'guide_assignment.mtx'}")
    except Exception as e:  # pragma: no cover
        log.warning("mtx write failed: %s", e)
    per_cell = np.asarray(A.sum(1)).ravel()
    per_grna = np.asarray(A.sum(0)).ravel()
    tab = pd.DataFrame({"grna_id": data.grna_ids, "grna_target": data.grna_target.values,
                        "n_cells_assigned": per_grna.astype(int),
                        "n_nonzero": np.asarray((data.grna != 0).sum(0)).ravel().astype(int)})
    if info:
        tab["assignment_method"] = [i["method"] for i in info]
        tab["pi"] = [i["pi"] for i in info]
        tab["log_fold"] = [i["log_fold"] for i in info]
    write_tsv(tab, d / "per_grna_assignment.tsv")
    if not args.no_plots:
        assignment_figure(d, A)
    summ = {"method": method, "moi": data.moi, "n_cells": data.n_cells, "n_grnas": len(data.grna_ids),
            "mean_grnas_per_cell": float(per_cell.mean()), "frac_cells_zero_grnas": float((per_cell == 0).mean()),
            "frac_cells_one_grna": float((per_cell == 1).mean()), "median_cells_per_grna": float(np.median(per_grna)),
            "dropped_covariates": dropped, "output": str(out)}
    write_json(summ, d / "summary.json")
    write_report(d, "sceptreIGVF: gRNA assignment", [
        f"Source: `{args.mudata}`; method **{method}** (sceptre::assign_grnas); output MuData `{out}` with layer "
        "`guide_assignment` in mod/guide; guide_assignment.mtx is guides x cells.",
        md_table(["metric", "value"], [[k, _fmt(v) if isinstance(v, float) else v] for k, v in summ.items()])])
    return 0


def cmd_qc(args) -> int:
    np, pd = _np(), _pd()
    mods, data, _ = _load(args, False)
    d = run_dir(args.label)
    A, method, _ = assign_grnas(data, args.method, threshold=args.threshold)
    keep, reasons = run_qc(data, A, p_mito_threshold=args.p_mito_threshold,
                           response_n_umis_range=tuple(args.response_n_umis_range),
                           response_n_nonzero_range=tuple(args.response_n_nonzero_range))
    write_tsv(pd.DataFrame({"cell_id": data.cell_ids, "pass_qc": keep}), d / "cells_in_use.tsv")
    summ = {"assign_method": method, "n_cells": data.n_cells, "n_cells_in_use": int(keep.sum()), "removed": reasons}
    write_json(summ, d / "summary.json")
    write_report(d, "sceptreIGVF: cell-wise QC", [
        md_table(["filter", "cells flagged"], [[k, v] for k, v in reasons.items()]),
        f"{int(keep.sum())} of {data.n_cells} cells pass (sceptre::run_qc defaults unless overridden)."])
    return 0


def _pipeline_tables(data: SceptreData, union_res, single_res):
    pd = _pd()
    gv = data.grna_var
    look = pd.DataFrame({"intended_target_key": data.grna_target.values})
    for c in ("intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end"):
        look[c] = gv[c].values if c in gv.columns else (data.grna_target.values if c == "intended_target_name" else "")
    look = look.drop_duplicates("intended_target_key")
    pe = union_res.rename(columns={"response_id": "gene_id", "grna_target": "intended_target_key",
                                   "log_2_fold_change": "log2_fc"}).merge(look, on="intended_target_key", how="left")
    pe = pe[["gene_id", "intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end",
             "log2_fc", "p_value"]]
    gid = gv["guide_id"].astype(str).values if "guide_id" in gv.columns else data.grna_ids
    gmap = dict(zip(data.grna_ids, gid))
    pg = single_res.rename(columns={"response_id": "gene_id", "log_2_fold_change": "log2_fc"})
    pg["guide_id"] = pg["grna_id"].map(gmap)
    pg = pg[["gene_id", "guide_id", "p_value", "log2_fc"]]
    return pe, pg


def cmd_inference(args) -> int:
    np, pd = _np(), _pd()
    pipeline = args.mode == "pipeline"
    tcol = "intended_target_key" if pipeline else "intended_target_name"
    mods, data, dropped = _load(args, pipeline, "greedy" if pipeline else "all-or-nothing", tcol)
    d = run_dir(args.label)
    pairs = pairs_from_uns(mods, data, tcol)
    if pairs is None:
        if not pipeline:
            raise SystemExit("uns['pairs_to_test'] missing (sceptreIGVF inference needs it; --mode pipeline falls back to trans pairs)")
        disc = construct_trans_pairs(data)
        pairs = disc.assign(intended_target_name=disc["grna_target"], gene_id=disc["response_id"])
    if pipeline:
        formula = args.formula if args.formula not in (None, "default") else auto_formula(
            data.covariates, data.low_moi, include_grna_covariates=False, extra=data.extra_covariate_names)
        control = "complement"
        if args.control_group not in ("default", "complement"):
            print(f"Overriding control_group='{args.control_group}' with 'complement' for SCEPTRE DE analysis.")
    else:
        formula = args.formula if args.formula not in (None, "default") else DEFAULT_INFERENCE_FORMULA
        control = args.control_group
    A = assign_thresholding(data.grna, 1)
    keep, reasons = run_qc(data, A, p_mito_threshold=1.0)
    tkw = _tester_kwargs(args)
    tkw["control_group"] = control
    tester = Tester(data, A, keep, formula, n_nonzero_trt_thresh=0, n_nonzero_cntrl_thresh=0, **tkw)
    disc_pairs = pairs[["grna_target", "response_id"]].dropna().drop_duplicates()
    disc_pairs = disc_pairs[disc_pairs["grna_target"].str.lower() != "nan"]
    uns = dict(mods.uns)
    outputs = {}
    if pipeline:
        union = run_pairs(tester, disc_pairs, "union", alpha=args.alpha)
        single = run_pairs(tester, disc_pairs, "singleton", alpha=args.alpha)
        pe, pg = _pipeline_tables(data, union, single)
        pre = f"{args.scope}_" if args.scope else ""
        uns[f"{pre}per_element_results"] = pe
        uns[f"{pre}per_guide_results"] = pg
        for k in ("per_guide_results", "per_element_results"):
            if pre and k in uns:
                del uns[k]
        outputs["per_element"] = write_tsv(pe, d / f"{pre}per_element_output.tsv.gz")
        outputs["per_guide"] = write_tsv(pg, d / f"{pre}per_guide_output.tsv.gz")
        write_tsv(union, d / "sceptre_union_results.tsv")
        write_tsv(single, d / "sceptre_singleton_results.tsv")
        main_res = union
        tr_desc = f"uns `{pre}per_element_results` ({len(pe)} rows) and `{pre}per_guide_results` ({len(pg)} rows)"
    else:
        disc = run_pairs(tester, disc_pairs, args.grna_integration_strategy, alpha=args.alpha)
        write_tsv(disc, d / "discovery_results.tsv")
        dr = disc.rename(columns={"response_id": "gene_id", "grna_target": "intended_target_name",
                                  "log_2_fold_change": "log2_fc"})
        dr = dr[["gene_id", "intended_target_name", "p_value", "log2_fc"]].drop_duplicates(["intended_target_name", "gene_id"])
        orig = [c for c in pairs.columns if c not in ("grna_target", "response_id")]
        test_results = pairs[orig].merge(dr, on=["intended_target_name", "gene_id"], how="left")
        uns["test_results"] = test_results
        outputs["test_results"] = write_tsv(test_results, d / "test_results.tsv")
        main_res = disc
        tr_desc = f"uns `test_results` ({len(test_results)} rows; columns {list(test_results.columns)})"
    al = _align(mods)
    out = Path(args.out) if args.out else d / "inference_mudata.h5mu"
    write_h5mu(out, al.gene, al.guide, al.obs, uns)
    if not args.no_plots:
        qq_figure(d, {"discovery": main_res["p_value"].values}, "Discovery p-values")
        volcano_figure(d, main_res, "Discovery pairs")
    summ = {"mode": args.mode, "scope": args.scope, "formula": formula, "control_group": tester.control_group,
            "resampling_mechanism": tester.mechanism, "side": args.side, "n_cells_in_use": int(len(tester.used)),
            "qc_removed": reasons, "n_pairs": int(len(disc_pairs)), "n_tested": int(main_res["p_value"].notna().sum()),
            "n_significant": int(main_res["significant"].sum()) if "significant" in main_res.columns else 0,
            "dropped_covariates": dropped, "output": str(out)}
    write_json(summ, d / "summary.json")
    write_report(d, "sceptreIGVF: inference", [
        f"Source: `{args.mudata}`.  Mode **{args.mode}** — gRNA assignment by thresholding at 1 UMI, relaxed QC "
        f"(n_nonzero thresholds 0, p_mito 1), formula `{formula}`, control group {tester.control_group}, "
        f"{tester.mechanism} null, side {args.side}.  Wrote {tr_desc} to `{out}`.",
        md_table(["metric", "value"], [[k, v] for k, v in summ.items() if k not in ("qc_removed",)]),
        "Top pairs:\n\n" + _result_md(main_res)])
    return 0


def _positive_or_die(pos):
    if pos is None or len(pos) == 0:
        raise SystemExit("no positive-control pairs: give --positive-control-pairs, or uns pairs_to_test with "
                         "pair_type == 'positive_control', or uns positive_control_pairs")


def _workflow_tester(args, data):
    A, method, _ = assign_grnas(data, args.method, threshold=args.threshold)
    keep, reasons = run_qc(data, A)
    f = args.formula if args.formula not in (None, "default") else auto_formula(data.covariates, data.low_moi,
                                                                                 extra=data.extra_covariate_names)
    tester = Tester(data, A, keep, f, n_nonzero_trt_thresh=args.n_nonzero_trt_thresh,
                    n_nonzero_cntrl_thresh=args.n_nonzero_cntrl_thresh, **_tester_kwargs(args))
    return tester, method, reasons


def cmd_calibration_check(args) -> int:
    mods, data, _ = _load(args, True)
    d = run_dir(args.label)
    tester, method, reasons = _workflow_tester(args, data)
    n = args.n_calibration_pairs
    if n is None:
        _, disc = split_pairs(mods, data, args)
        n = max(1, len(disc))
    res, summ = run_calibration_check(tester, n, alpha=args.alpha, group_size=args.group_size, seed=args.seed)
    write_tsv(res, d / "calibration_check_results.tsv")
    summ.update({"assign_method": method, "control_group": tester.control_group, "resampling_mechanism": tester.mechanism,
                 "formula": tester.formula})
    write_json(summ, d / "summary.json")
    if not args.no_plots:
        qq_figure(d, {"negative control": res["p_value"].values}, "Calibration check (NT pairs)")
    write_report(d, "sceptreIGVF: calibration check", [
        f"{summ['n_calibration_pairs']} undercover non-targeting pairs (groups of {summ['undercover_group_size']} NT "
        f"gRNAs x random responses).  A calibrated test makes ~0 discoveries at BH {args.alpha}: here "
        f"**{summ['n_false_discoveries']}**; median p {_fmt(summ['median_p'])}; mean log2 FC {_fmt(summ['mean_log2_fc'])}.",
        md_table(["metric", "value"], [[k, _fmt(v) if isinstance(v, float) else v] for k, v in summ.items()])])
    return 0


def cmd_power_check(args) -> int:
    mods, data, _ = _load(args, True)
    d = run_dir(args.label)
    pos, _ = split_pairs(mods, data, args)
    _positive_or_die(pos)
    tester, method, reasons = _workflow_tester(args, data)
    res = run_pairs(tester, pos, args.grna_integration_strategy, correct=False)
    write_tsv(res, d / "power_check_results.tsv")
    ok = res["p_value"].notna()
    summ = {"n_positive_control_pairs": int(len(res)), "n_tested": int(ok.sum()),
            "median_p": float(res.loc[ok, "p_value"].median()) if ok.any() else float("nan"),
            "frac_p_below_1e-5": float((res.loc[ok, "p_value"] < 1e-5).mean()) if ok.any() else float("nan"),
            "median_log2_fc": float(res.loc[ok, "log_2_fold_change"].median()) if ok.any() else float("nan"),
            "assign_method": method, "control_group": tester.control_group}
    write_json(summ, d / "summary.json")
    write_report(d, "sceptreIGVF: power check", [
        f"{len(res)} positive-control pairs; median p {_fmt(summ['median_p'])}; median log2 FC {_fmt(summ['median_log2_fc'])}.",
        _result_md(res)])
    return 0


def cmd_run(args) -> int:
    np, pd = _np(), _pd()
    mods, data, dropped = _load(args, True)
    d = run_dir(args.label)
    pos, disc = split_pairs(mods, data, args)
    qc_kw = {"n_nonzero_trt_thresh": args.n_nonzero_trt_thresh, "n_nonzero_cntrl_thresh": args.n_nonzero_cntrl_thresh,
             "p_mito_threshold": args.p_mito_threshold}
    formula = args.formula if args.formula not in (None, "default") else None
    W = full_workflow(data, pos, disc, formula, args.method, _tester_kwargs(args), qc_kw, args.grna_integration_strategy,
                      args.alpha, args.n_calibration_pairs, not args.skip_calibration,
                      assign_kw={"threshold": args.threshold})
    write_tsv(W["discovery"], d / "discovery_results.tsv")
    if W["power"] is not None:
        write_tsv(W["power"], d / "power_check_results.tsv")
    if W["calibration"] is not None:
        write_tsv(W["calibration"], d / "calibration_check_results.tsv")
    gene, guide, obs, uns = sceptre_to_mudata(data, W["cells_in_use"], W["assignment"], pos, disc,
                                              W["power"], W["discovery"], gene_var=_align(mods).gene.var)
    out = Path(args.out) if args.out else d / "sceptre_mudata.h5mu"
    write_h5mu(out, gene, guide, obs, uns)
    write_tsv(uns["test_results"], d / "test_results.tsv")
    if not args.no_plots:
        series = {"discovery": W["discovery"]["p_value"].values}
        if W["calibration"] is not None:
            series["negative control"] = W["calibration"]["p_value"].values
        if W["power"] is not None:
            series["positive control"] = W["power"]["p_value"].values
        qq_figure(d, series, "SCEPTRE p-values")
        volcano_figure(d, W["discovery"], "Discovery pairs")
        assignment_figure(d, W["assignment"])
    pw = W["power"]
    summ = {"moi": data.moi, "assign_method": W["assign_method"], "formula": W["formula"],
            "control_group": W["tester"].control_group, "resampling_mechanism": W["tester"].mechanism,
            "n_cells": data.n_cells, "n_cells_in_use": int(W["cells_in_use"].sum()), "qc_removed": W["qc_reasons"],
            "n_discovery_pairs": int(len(disc)), "n_discovery_tested": int(W["discovery"]["p_value"].notna().sum()),
            "n_discoveries": int(W["discovery"]["significant"].sum()),
            "n_positive_control_pairs": int(len(pos)),
            "power_median_p": float(pw["p_value"].median()) if pw is not None and len(pw) else None,
            "calibration": W["calibration_summary"], "dropped_covariates": dropped, "output": str(out)}
    write_json(summ, d / "summary.json")
    cs = W["calibration_summary"]
    write_report(d, "sceptreIGVF: SCEPTRE workflow", [
        f"Source: `{args.mudata}` ({data.moi} MOI).  Assignment **{W['assign_method']}**; "
        f"{summ['n_cells_in_use']}/{data.n_cells} cells pass QC; formula `{W['formula']}`; "
        f"{W['tester'].control_group} control, {W['tester'].mechanism} null.",
        f"Calibration: {cs.get('n_false_discoveries', 'NA')} false discoveries in {cs.get('n_calibration_pairs', 0)} NT pairs "
        f"(median p {_fmt(cs.get('median_p'))}).  Power: median positive-control p {_fmt(summ['power_median_p'])}.  "
        f"Discovery: {summ['n_discoveries']} of {summ['n_discovery_tested']} pairs significant at BH {args.alpha}.",
        "Top discovery pairs:\n\n" + _result_md(W["discovery"]),
        f"Output MuData (sceptre_object_to_mudata schema): `{out}`."])
    return 0


def _minimal(gene, guide, obs, uns):
    import anndata as ad  # type: ignore
    pd = _pd()
    g2 = ad.AnnData(X=gene.X, obs=pd.DataFrame(index=gene.obs_names), var=pd.DataFrame(index=gene.var_names))
    gv = guide.var[["targeting", "intended_target_name"]].copy()
    u2 = ad.AnnData(X=guide.X, obs=pd.DataFrame(index=guide.obs_names), var=gv)
    for k, v in guide.layers.items():
        u2.layers[k] = v
    u2.uns = dict(guide.uns)
    o2 = pd.DataFrame(index=obs.index)
    un2 = {}
    if "pairs_to_test" in uns:
        un2["pairs_to_test"] = uns["pairs_to_test"][["intended_target_name", "gene_id"]]
    if "test_results" in uns:
        un2["test_results"] = uns["test_results"][["intended_target_name", "gene_id", "p_value"]]
    return g2, u2, o2, un2


def build_mudata_list(gene, guide, obs, uns) -> "dict[str, tuple]":
    """sceptre_object_to_mudata_inputs_outputs: the eight objects keyed by name."""
    def strip(g, u, o, un, drop_tr=False, drop_pairs=False, drop_layer=False):
        un = dict(un)
        if drop_tr:
            un.pop("test_results", None)
        if drop_pairs:
            un.pop("pairs_to_test", None)
        if drop_layer:
            u = u.copy()
            if "guide_assignment" in u.layers:
                del u.layers["guide_assignment"]
        return g, u, o, un
    out = {}
    io = (gene, guide, obs, uns)
    out["inference_output"] = io
    out["inference_input"] = strip(*io, drop_tr=True)
    out["guide_assignment_output"] = strip(*out["inference_input"], drop_pairs=True)
    out["guide_assignment_input"] = strip(*out["guide_assignment_output"], drop_layer=True)
    mn = _minimal(*io)
    out["inference_output_minimal"] = mn
    out["inference_input_minimal"] = strip(*mn, drop_tr=True)
    out["guide_assignment_output_minimal"] = strip(*out["inference_input_minimal"], drop_pairs=True)
    out["guide_assignment_input_minimal"] = strip(*out["guide_assignment_output_minimal"], drop_layer=True)
    return out


def save_mudata_list(objs: dict, path: Path, prefix: str = "") -> "list[Path]":
    written = []
    for name, (g, u, o, un) in objs.items():
        if "guide_assignment" in name:
            sub = path / "guide_assignment"
        elif "inference" in name:
            sub = path / "inference"
        else:
            raise SystemExit("mudata_name must contain either 'guide_assignment' or 'inference'")
        written.append(write_h5mu(sub / f"{prefix}{name}.h5mu", g, u, o, un))
    return written


def cmd_export_mudata(args) -> int:
    np, pd = _np(), _pd()
    mods, data, _ = _load(args, True)
    d = run_dir(args.label)
    pos, disc = split_pairs(mods, data, args)
    tkw = _tester_kwargs(args)
    W = full_workflow(data, pos, disc, None, "default", tkw, run_calibration=False)
    pw = W["power"]
    pos2 = pw.dropna(subset=["p_value"])[["response_id", "grna_target"]] if pw is not None else pd.DataFrame(columns=["response_id", "grna_target"])
    dr = W["discovery"].dropna(subset=["p_value"])
    n_sig, n_non = int(dr["significant"].sum()), int((~dr["significant"].astype(bool)).sum())
    k_sig = min(n_sig, int(round(args.num_discovery_pairs / 2)))
    k_non = min(n_non, args.num_discovery_pairs - k_sig)
    rng = np.random.default_rng(args.seed)
    sig_rows = dr[dr["significant"].astype(bool)]
    non_rows = dr[~dr["significant"].astype(bool)]
    disc2 = pd.concat([sig_rows.iloc[rng.choice(len(sig_rows), k_sig, replace=False)] if k_sig else sig_rows.iloc[:0],
                       non_rows.iloc[rng.choice(len(non_rows), k_non, replace=False)] if k_non else non_rows.iloc[:0]])[["response_id", "grna_target"]]
    W2 = full_workflow(data, pos2, disc2, W["formula"], "default", tkw, run_calibration=False)
    gene_var = pd.read_csv(args.gene_info, sep=None, engine="python", index_col=0) if args.gene_info else None
    gene, guide, obs, uns = sceptre_to_mudata(data, W2["cells_in_use"], W2["assignment"], pos2, disc2, W2["power"],
                                              W2["discovery"], gene_var=gene_var)
    if args.guide_capture_method:
        guide.uns["capture_method"] = args.guide_capture_method
    out = Path(args.out_dir) if args.out_dir else d
    written = save_mudata_list(build_mudata_list(gene, guide, obs, uns), out, args.prefix)
    summ = {"n_positive_control_pairs": int(len(pos2)), "n_discovery_pairs": int(len(disc2)),
            "n_significant_kept": k_sig, "n_non_significant_kept": k_non, "files": [str(p) for p in written]}
    write_json(summ, d / "summary.json")
    write_report(d, "sceptreIGVF: benchmark MuData export", [
        f"{len(pos2)} positive controls and {len(disc2)} discovery pairs ({k_sig} significant + {k_non} not), "
        "re-analysed with default assignment and QC, exported as the eight inputs/outputs objects "
        "(save_mudata_list layout).",
        "\n".join(f"- `{p}`" for p in written)])
    return 0


R_TEMPLATE = """suppressPackageStartupMessages({{library(MuData); library(sceptreIGVF)}})
m <- MuData::readH5MU("{inp}")
{call}
MuData::writeH5MU(object = out, file = "{out}")
tr <- MultiAssayExperiment::metadata(out)$test_results
if (!is.null(tr)) write.table(as.data.frame(tr), file = "{tsv}", sep = "\\t", row.names = FALSE, quote = FALSE)
"""


def cmd_r_runner(args) -> int:
    d = run_dir(args.label)
    out = Path(args.out) if args.out else d / f"r_{args.step}_output.h5mu"
    extra = []
    if args.side:
        extra.append(f'side = "{args.side}"')
    if args.formula and args.formula != "default":
        extra.append(f"formula_object = stats::formula({args.formula})")
    if args.step == "assign":
        call = "out <- sceptreIGVF::assign_grnas_sceptre(m)"
    else:
        call = "out <- sceptreIGVF::inference_sceptre(m" + (", " + ", ".join(extra) if extra else "") + ")"
    script = d / f"run_sceptreIGVF_{args.step}.R"
    script.write_text(R_TEMPLATE.format(inp=Path(args.mudata).resolve(), call=call, out=out, tsv=d / "r_test_results.tsv"))
    print(f"Wrote: {script}")
    cmd = ["Rscript", str(script)]
    rs = shutil.which("Rscript")
    status = "not_run"
    msg = ""
    if rs is None:
        msg = ("Rscript not on PATH: install R, then Bioconductor MuData and "
               f"remotes::install_github('{UPSTREAM_REPO}@{UPSTREAM_COMMIT}') (pulls in sceptre); command that would run:")
    else:
        chk = subprocess.run([rs, "-e", "suppressPackageStartupMessages({library(MuData); library(sceptreIGVF)})"],
                             capture_output=True, text=True)
        if chk.returncode != 0:
            msg = "R is present but MuData/sceptreIGVF are not installed; command that would run:"
        else:
            p = subprocess.run(cmd, capture_output=True, text=True)
            (d / "r_stdout.txt").write_text(p.stdout + "\n" + p.stderr)
            status = "ok" if p.returncode == 0 else f"failed ({p.returncode})"
            if p.returncode == 0:
                print(f"Wrote: {out}")
    if msg:
        print(msg)
    print("  " + " ".join(cmd))
    summ = {"step": args.step, "status": status, "command": cmd, "script": str(script), "output": str(out)}
    write_json(summ, d / "summary.json")
    write_report(d, "sceptreIGVF: R runner", [f"Status: **{status}**. {msg}", "```\n" + " ".join(cmd) + "\n```",
                                              "```r\n" + script.read_text() + "```"])
    return 0


def cmd_compare_results(args) -> int:
    np, pd = _np(), _pd()
    from scipy import stats  # type: ignore
    d = run_dir(args.label)
    a = pd.read_csv(args.a, sep=None, engine="python")
    b = pd.read_csv(args.b, sep=None, engine="python")
    keys = [k for k in args.keys if k in a.columns and k in b.columns]
    if not keys:
        raise SystemExit(f"no shared key columns among {args.keys}")
    m = a.merge(b, on=keys, suffixes=("_a", "_b"))
    ok = m["p_value_a"].notna() & m["p_value_b"].notna()
    la = -np.log10(m.loc[ok, "p_value_a"].clip(1e-300))
    lb = -np.log10(m.loc[ok, "p_value_b"].clip(1e-300))
    rho = float(stats.spearmanr(la, lb).correlation) if ok.sum() > 2 else float("nan")
    qa, qb = bh_adjust(m["p_value_a"].values), bh_adjust(m["p_value_b"].values)
    sa, sb = np.nan_to_num(qa, nan=1) < args.alpha, np.nan_to_num(qb, nan=1) < args.alpha
    summ = {"n_a": int(len(a)), "n_b": int(len(b)), "n_matched": int(len(m)), "n_both_tested": int(ok.sum()),
            "spearman_neglog10p": rho, "n_sig_a": int(sa.sum()), "n_sig_b": int(sb.sum()),
            "n_sig_both": int((sa & sb).sum()), "significance_agreement": float((sa == sb).mean()) if len(m) else float("nan")}
    fc = [c for c in ("log2_fc", "log_2_fold_change") if f"{c}_a" in m.columns and f"{c}_b" in m.columns]
    if fc:
        f_ok = m[f"{fc[0]}_a"].notna() & m[f"{fc[0]}_b"].notna()
        summ["pearson_log2_fc"] = float(np.corrcoef(m.loc[f_ok, f"{fc[0]}_a"], m.loc[f_ok, f"{fc[0]}_b"])[0, 1]) if f_ok.sum() > 2 else float("nan")
    write_tsv(m, d / "merged_results.tsv")
    write_json(summ, d / "summary.json")
    write_report(d, "sceptreIGVF: results comparison", [f"A = `{args.a}`, B = `{args.b}`, keys {keys}.",
                                                       md_table(["metric", "value"], [[k, _fmt(v) if isinstance(v, float) else v] for k, v in summ.items()])])
    return 0


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------

def make_synthetic(path: Path, moi: str = "high", n_cells: int = 1500, seed: int = 11) -> dict:
    """IGVF-schema MuData with planted effects: T1 -> G000 (log -1.2), T2 -> G001 (-1.0), T3 -> G002 (+0.8)."""
    import anndata as ad  # type: ignore
    np, pd, sp = _np(), _pd(), _sp()
    rng = np.random.default_rng(seed)
    targets = [f"T{i}" for i in range(1, 7)]
    guides, gtargets = [], []
    for t in targets:
        for k in range(2):
            guides.append(f"{t}_g{k + 1}")
            gtargets.append(t)
    n_nt = 8
    for k in range(n_nt):
        guides.append(f"NT_g{k + 1}")
        gtargets.append(NT_LABEL)
    K = len(guides)
    if moi == "high":
        present = rng.random((n_cells, K)) < 0.08
    else:
        present = np.zeros((n_cells, K), bool)
        pick = rng.integers(0, K, n_cells)
        pick = np.where(rng.random(n_cells) < 0.3, rng.integers(12, K, n_cells), pick)  # enrich NT cells
        present[np.arange(n_cells), pick] = True
        extra = rng.random(n_cells) < 0.03
        present[np.where(extra)[0], rng.integers(0, K, int(extra.sum()))] = True
    s = rng.lognormal(0, 0.3, n_cells)
    gcounts = np.where(present, rng.poisson(25 * s[:, None], (n_cells, K)), rng.poisson(0.2, (n_cells, K)))
    G = 30
    gene_ids = [f"ENSG{i:011d}" for i in range(G - 2)] + ["ENSG_MT1", "ENSG_MT2"]
    gene_names = [f"G{i:03d}" for i in range(G - 2)] + ["MT-CO1", "MT-ND1"]
    batch = rng.integers(0, 2, n_cells)
    L = rng.lognormal(0, 0.35, n_cells)
    base = rng.uniform(-0.5, 1.5, G)
    base[-2:] = 0.3
    eta = base[None, :] + np.log(L)[:, None] + 0.3 * batch[:, None]
    tmask = {t: present[:, [i for i, gt in enumerate(gtargets) if gt == t]].any(1) for t in targets}
    effects = {("T1", 0): -1.2, ("T2", 1): -1.0, ("T3", 2): 0.8}
    for (t, g), e in effects.items():
        eta[:, g] += e * tmask[t]
    mu = np.exp(eta)
    theta = 5.0
    Y = rng.poisson(rng.gamma(theta, mu / theta))
    obs = pd.DataFrame({"batch": np.where(batch == 1, "b2", "b1")}, index=[f"cell{i}" for i in range(n_cells)])
    gvar = pd.DataFrame({"gene_name": gene_names, "chr": ["chr1"] * 3 + ["chr2"] * (G - 3),
                         "start": [1_000_000, 5_000_000, 9_000_000] + list(range(1_000_000, 1_000_000 + 100_000 * (G - 3), 100_000)),
                         "end": [1_010_000, 5_010_000, 9_010_000] + list(range(1_010_000, 1_010_000 + 100_000 * (G - 3), 100_000)),
                         "strand": "+"}, index=gene_ids)
    tpos = {"T1": ("chr1", 1_200_000), "T2": ("chr1", 5_100_000), "T3": ("chr1", 9_300_000), "T4": ("chr3", 100),
            "T5": ("chr3", 50_000_000), "T6": ("chr3", 90_000_000)}
    uvar = pd.DataFrame({
        "guide_id": guides, "intended_target_name": gtargets,
        "targeting": ["FALSE" if t == NT_LABEL else "TRUE" for t in gtargets],
        "type": ["non-targeting" if t == NT_LABEL else "targeting" for t in gtargets],
        "intended_target_chr": [tpos[t][0] if t in tpos else "" for t in gtargets],
        "intended_target_start": [tpos[t][1] if t in tpos else -9 for t in gtargets],
        "intended_target_end": [tpos[t][1] + 500 if t in tpos else -9 for t in gtargets]}, index=guides)
    uvar["intended_target_key"] = uvar["intended_target_name"]
    gene = ad.AnnData(X=sp.csr_matrix(Y.astype(np.float32)), obs=pd.DataFrame(index=obs.index), var=gvar)
    guide = ad.AnnData(X=sp.csr_matrix(gcounts.astype(np.float32)), obs=pd.DataFrame(index=obs.index), var=uvar)
    guide.uns["moi"] = moi
    pc = [("T1", gene_ids[0]), ("T2", gene_ids[1])]
    disc = [(t, gene_ids[g]) for t in ["T3", "T4", "T5", "T6"] for g in range(2, 8)]
    ptt = pd.DataFrame(pc + disc, columns=["intended_target_name", "gene_id"])
    ptt["pair_type"] = ["positive_control"] * len(pc) + ["discovery"] * len(disc)
    write_h5mu(path, gene, guide, obs, {"pairs_to_test": ptt}, announce=False)
    return {"path": path, "present": present, "gene_ids": gene_ids, "planted_disc": ("T3", gene_ids[2]),
            "pc": pc, "disc": disc, "effects": effects}


def cmd_selftest(args) -> int:
    import tempfile
    np, pd = _np(), _pd()
    checks: list = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    def ns(**kw):
        base = dict(moi=None, side="both", control_group="default", resampling_mechanism="default",
                    resampling_approximation="skew_normal", b1=499, b2=4999, b3=24999, seed=4, no_plots=args.no_plots,
                    alpha=0.1, formula=None, out=None, method="default", threshold=None,
                    n_nonzero_trt_thresh=7, n_nonzero_cntrl_thresh=7, positive_control_pairs=None, discovery_pairs=None,
                    pair_set=None, cis_distance=500_000, n_calibration_pairs=None, group_size=None,
                    grna_integration_strategy="union", p_mito_threshold=0.2, skip_calibration=False)
        base.update(kw)
        return argparse.Namespace(**base)

    before = set(OUT_ROOT.glob("*")) if OUT_ROOT.exists() else set()

    def latest(label):
        c = sorted(p for p in OUT_ROOT.glob(f"*_{label}*") if p not in before)
        return c[-1] if c else None

    t0 = time.time()
    try:
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            print("\nstatistics primitives")
            rng = np.random.default_rng(0)
            from scipy import stats  # type: ignore
            x = stats.skewnorm.rvs(4, loc=1, scale=2, size=20000, random_state=1)
            xi, om, al = skewnorm_mom(x)
            check(abs(xi - 1) < 0.15 and abs(om - 2) < 0.15 and 2 < al < 8, f"skewnorm_mom recovers (1, 2, 4): ({xi:.2f}, {om:.2f}, {al:.2f})")
            null = np.arange(1, 100, dtype=float)
            check(abs(empirical_p(99.5, null, "right") - 1 / 100) < 1e-12 and abs(empirical_p(0.0, null, "left") - 1 / 100) < 1e-12
                  and empirical_p(50, null, "both") <= 1.0, "empirical_p = (1 + #extreme) / (1 + B), sides")
            check(np.allclose(bh_adjust([0.01, 0.04, 0.03, 0.5]), [0.04, 0.16 / 3, 0.16 / 3, 0.5]), "bh_adjust matches p.adjust('BH')")
            y = rng.poisson(rng.gamma(4.0, 3.0 * 2 ** -1.0 / 4.0, 4000))
            check(abs(log2_fold_change(y, np.full(4000, 3.0), 4.0) + 1.0) < 0.08, "log2_fold_change: NB MLE with offset recovers -1")
            cov = pd.DataFrame({"a": ["x", "y"] * 5, "b": ["x", "y"] * 5, "c": np.arange(10.0)})
            e1, d1 = prepare_extra_covariates(cov, True, "all-or-nothing")
            e2, d2 = prepare_extra_covariates(cov, True, "greedy")
            check(e1.shape[1] == 0 and list(e2.columns) == ["a", "c"] and d2 == ["b"],
                  "collinear covariates: sceptreIGVF drops all, pipeline greedy drops only the duplicate")
            f = auto_formula(pd.DataFrame({"response_n_nonzero": [1.0, 2], "response_n_umis": [0.0, 5], "grna_n_nonzero": [1.0, 2],
                                           "grna_n_umis": [3.0, 4], "batch": ["a", "b"]}), False, extra=["batch"])
            check(f == "~ log(response_n_nonzero) + log(response_n_umis + 1) + log(grna_n_nonzero) + log(grna_n_umis) + batch",
                  f"auto_formula (high MOI, zero-guard, extra covariates): {f}")

            print("\nhigh-MOI synthetic screen")
            H = make_synthetic(td / "high.h5mu", "high")
            mods = read_h5mu(H["path"])
            data, _ = convert_mudata_to_sceptre(mods, True)
            check(data.moi == "high" and data.n_cells == 1500 and "response_p_mito" in data.covariates.columns
                  and data.extra_covariate_names == ["batch"] and int(data.grna_is_nt.sum()) == 8,
                  "convert: moi, cells, computed covariates (p_mito from MT- genes), extra covariate, NT gRNAs")
            rc = cmd_convert(ns(mudata=str(H["path"]), remove_collinear_covariates=True, collinear_mode="all-or-nothing",
                                target_col="intended_target_name", label="st_sceptre_convert"))
            cd = latest("st_sceptre_convert")
            check(rc == 0 and cd is not None and (cd / "covariate_data_frame.tsv").exists() and (cd / "grna_matrix.npz").exists(),
                  "convert subcommand writes the sceptre object directory")

            thr = assign_thresholding(data.grna, 5)
            check(np.array_equal(thr.toarray(), data.grna.toarray() >= 5), "thresholding assignment = counts >= threshold")
            rc = cmd_assign_guides(ns(mudata=str(H["path"]), method="mixture", probability_threshold=0.8, n_em_rep=5,
                                      collinear_mode="all-or-nothing", label="st_sceptre_assign"))
            ad_ = latest("st_sceptre_assign")
            am = read_h5mu(ad_ / "guide_assignment_output.h5mu")
            G_ = am.guide.layers["guide_assignment"]
            G_ = G_.toarray() if hasattr(G_, "toarray") else np.asarray(G_)
            acc = float((G_.astype(bool) == H["present"]).mean())
            sens = float(G_.astype(bool)[H["present"]].mean())
            check(rc == 0 and G_.shape == H["present"].shape and set(np.unique(G_)) <= {0, 1} and acc > 0.98 and sens > 0.95,
                  f"assign-guides (mixture): layer guide_assignment written, 0/1, accuracy {acc:.3f}, sensitivity {sens:.3f}")
            check((ad_ / "guide_assignment.mtx").exists(), "assign-guides writes guide_assignment.mtx (guides x cells)")

            rc = cmd_qc(ns(mudata=str(H["path"]), method="thresholding", threshold=5, response_n_umis_range=[0.01, 0.99],
                           response_n_nonzero_range=[0.01, 0.99], label="st_sceptre_qc"))
            qs = json.loads((latest("st_sceptre_qc") / "summary.json").read_text())
            check(rc == 0 and 1400 <= qs["n_cells_in_use"] < 1500, f"qc: 1%-99% range filters keep {qs['n_cells_in_use']} of 1500 cells")

            print("\ninference (sceptreIGVF inference_sceptre on the assign-guides output, as test-integration.R)")
            inf_in = ad_ / "guide_assignment_output.h5mu"
            check("pairs_to_test" in am.uns, "assign-guides output keeps uns pairs_to_test")
            rc = cmd_inference(ns(mudata=str(inf_in), mode="sceptreIGVF", scope=None, label="st_sceptre_inf"))
            idir = latest("st_sceptre_inf")
            im = read_h5mu(idir / "inference_mudata.h5mu")
            tr = uns_table(im.uns, "test_results")
            check(rc == 0 and list(tr.columns) == ["intended_target_name", "gene_id", "pair_type", "p_value", "log2_fc"]
                  and len(tr) == 26, f"test_results columns = pairs_to_test + p_value, log2_fc; {len(tr)} rows (upstream test)")
            pcp = tr[tr["pair_type"] == "positive_control"]["p_value"].astype(float)
            check(pcp.median() < 1e-10, f"median positive-control p < 1e-10 (upstream test): {pcp.median():.2e}")
            t1 = tr[(tr["intended_target_name"] == "T1") & (tr["gene_id"] == H["gene_ids"][0])].iloc[0]
            check(abs(float(t1["log2_fc"]) - (-1.2 / math.log(2))) < 0.35, f"log2_fc of T1->G000 {float(t1['log2_fc']):.2f} ~ planted {-1.2 / math.log(2):.2f}")
            pl = tr[(tr["intended_target_name"] == "T3") & (tr["gene_id"] == H["gene_ids"][2])].iloc[0]
            nulls = tr[(tr["pair_type"] == "discovery") & ~((tr["intended_target_name"] == "T3") & (tr["gene_id"] == H["gene_ids"][2]))]
            check(float(pl["p_value"]) < 1e-6 and float(pl["log2_fc"]) > 0.6, f"planted discovery T3->G002 found: p {float(pl['p_value']):.1e}")
            nfp = int((nulls["p_value"].astype(float) < 0.01).sum())
            check(nfp <= 1 and nulls["p_value"].astype(float).median() > 0.2, f"null discovery pairs calibrated: {nfp} of {len(nulls)} with p < 0.01, median p {nulls['p_value'].astype(float).median():.2f}")

            rc = cmd_inference(ns(mudata=str(inf_in), mode="sceptreIGVF", scope=None, side="left",
                                  formula="~ batch + log(response_n_nonzero) + log(response_n_umis)", label="st_sceptre_inf_left"))
            tr2 = pd.read_csv(latest("st_sceptre_inf_left") / "test_results.tsv", sep="\t")
            p3 = float(tr2[(tr2["intended_target_name"] == "T3") & (tr2["gene_id"] == H["gene_ids"][2])]["p_value"].iloc[0])
            check(rc == 0 and tr2[tr2["pair_type"] == "positive_control"]["p_value"].median() < 1e-10 and p3 > 0.5,
                  f"inference side=left + batch formula (upstream test): repression found, up-regulated T3 not (p {p3:.2f})")

            rc = cmd_inference(ns(mudata=str(inf_in), mode="pipeline", scope="cis", label="st_sceptre_inf_pipe"))
            pm = read_h5mu(latest("st_sceptre_inf_pipe") / "inference_mudata.h5mu")
            pe, pg = uns_table(pm.uns, "cis_per_element_results"), uns_table(pm.uns, "cis_per_guide_results")
            check(rc == 0 and pe is not None and pg is not None and list(pe.columns) == ["gene_id", "intended_target_name", "intended_target_chr",
                  "intended_target_start", "intended_target_end", "log2_fc", "p_value"] and list(pg.columns) == ["gene_id", "guide_id", "p_value", "log2_fc"]
                  and len(pe) == 26 and len(pg) == 52, "inference --mode pipeline --scope cis: cis_per_element_results / cis_per_guide_results schema")
            g1 = pg[(pg["guide_id"].isin(["T1_g1", "T1_g2"])) & (pg["gene_id"] == H["gene_ids"][0])]["p_value"].astype(float)
            check(len(g1) == 2 and (g1 < 1e-4).all(), "per-guide (singleton) results: both T1 gRNAs repress G000")

            print("\ncalibration / power / run")
            rc = cmd_calibration_check(ns(mudata=str(H["path"]), n_calibration_pairs=40, label="st_sceptre_calib"))
            cs = json.loads((latest("st_sceptre_calib") / "summary.json").read_text())
            check(rc == 0 and cs["n_false_discoveries"] <= 1 and cs["median_p"] > 0.2 and cs.get("ks_uniform_p", 0) > 0.01,
                  f"calibration check (undercover NT groups of {cs['undercover_group_size']}): {cs['n_false_discoveries']} false discoveries, "
                  f"median p {cs['median_p']:.2f}, KS p {cs.get('ks_uniform_p', float('nan')):.2f}")
            rc = cmd_power_check(ns(mudata=str(H["path"]), label="st_sceptre_power"))
            ps = json.loads((latest("st_sceptre_power") / "summary.json").read_text())
            check(rc == 0 and ps["median_p"] < 1e-8 and ps["median_log2_fc"] < -1, f"power check: median p {ps['median_p']:.1e}, log2 FC {ps['median_log2_fc']:.2f}")
            rc = cmd_run(ns(mudata=str(H["path"]), label="st_sceptre_run"))
            rd = latest("st_sceptre_run")
            rm = read_h5mu(rd / "sceptre_mudata.h5mu")
            rpt = uns_table(rm.uns, "pairs_to_test")
            rtr = uns_table(rm.uns, "test_results")
            check(rc == 0 and list(rpt.columns) == ["intended_target_name", "gene_id", "pair_type"] and list(rpt["pair_type"][:2]) == ["positive_control"] * 2
                  and list(rtr.columns) == ["intended_target_name", "gene_id", "pair_type", "p_value", "log2_fc"],
                  "run: sceptre_object_to_mudata uns schema (positive controls first)")
            check("guide_assignment" in rm.guide.layers and list(rm.guide.var.columns[:2]) == ["targeting", "intended_target_name"]
                  and set(rm.guide.var["targeting"]) == {"TRUE", "FALSE"} and {"num_expressed_guides", "total_guide_umis"} <= set(rm.guide.obs.columns)
                  and {"num_expressed_genes", "total_gene_umis"} <= set(rm.gene.obs.columns) and "batch" in rm.obs.columns,
                  "run: guide var targeting/intended_target_*, obs renames, covariates in top-level obs")
            rs = json.loads((rd / "summary.json").read_text())
            check(rs["n_discoveries"] >= 1 and rs["calibration"]["n_false_discoveries"] <= 1 and rs["power_median_p"] < 1e-8,
                  f"run: {rs['n_discoveries']} discoveries, {rs['calibration']['n_false_discoveries']} calibration false discoveries")

            cp = construct_cis_pairs(data, mods.gene.var, 500_000)
            check(("T1", H["gene_ids"][0]) in set(map(tuple, cp.values)) and not any(cp["grna_target"] == "T4"),
                  f"construct_cis_pairs: T1-G000 within 500 kb kept, distant T4 has none ({len(cp)} pairs)")
            check(len(construct_trans_pairs(data)) == 6 * 30, "construct_trans_pairs: 6 targets x 30 responses")

            print("\nexport-mudata (sceptre_object_to_mudata_inputs_outputs)")
            ed = td / "export"
            rc = cmd_export_mudata(ns(mudata=str(H["path"]), num_discovery_pairs=6, out_dir=str(ed), prefix="syn_",
                                      gene_info=None, guide_capture_method="direct capture", label="st_sceptre_export"))
            names = ["inference_output", "inference_input", "guide_assignment_output", "guide_assignment_input"]
            files = [ed / ("inference" if "inference" in n else "guide_assignment") / f"syn_{n}{s}.h5mu" for n in names for s in ("", "_minimal")]
            check(rc == 0 and all(p.exists() for p in files), "export-mudata: eight h5mu files in guide_assignment/ and inference/")
            gi = read_h5mu(ed / "guide_assignment" / "syn_guide_assignment_input.h5mu")
            ii = read_h5mu(ed / "inference" / "syn_inference_input.h5mu")
            om = read_h5mu(ed / "inference" / "syn_inference_output_minimal.h5mu")
            check("guide_assignment" not in gi.guide.layers and "pairs_to_test" not in gi.uns and "test_results" not in ii.uns
                  and "pairs_to_test" in ii.uns, "export-mudata: input objects drop test_results / pairs_to_test / guide_assignment layer")
            check(list(om.guide.var.columns) == ["targeting", "intended_target_name"] and list(uns_table(om.uns, "test_results").columns) == ["intended_target_name", "gene_id", "p_value"]
                  and om.obs.shape[1] == 0 and str(np.ravel(gi.guide.uns.get("capture_method"))[0]) == "direct capture",
                  "export-mudata: minimal schema; capture_method in guide uns")
            es = json.loads((latest("st_sceptre_export") / "summary.json").read_text())
            check(es["n_discovery_pairs"] == 6 and es["n_significant_kept"] <= 3, f"export-mudata: {es['n_significant_kept']} significant + {es['n_non_significant_kept']} non-significant discovery pairs")

            print("\nlow-MOI synthetic screen")
            Lw = make_synthetic(td / "low.h5mu", "low", seed=5)
            ldata, _ = convert_mudata_to_sceptre(read_h5mu(Lw["path"]), True)
            A_max = assign_maximum(ldata.grna)
            one = Lw["present"].sum(1) == 1
            accm = float((A_max.toarray()[one] == Lw["present"][one]).all(1).mean())
            check(ldata.low_moi and accm > 0.97, f"maximum assignment recovers the single gRNA in {accm:.3f} of single-gRNA cells")
            keep, reasons = run_qc(ldata, A_max)
            check(reasons.get("zero_or_twoplus_grnas", 0) > 0, f"low-MOI QC removes cells with 0 or 2+ gRNAs ({reasons.get('zero_or_twoplus_grnas')})")
            rc = cmd_run(ns(mudata=str(Lw["path"]), label="st_sceptre_low"))
            ls = json.loads((latest("st_sceptre_low") / "summary.json").read_text())
            check(rc == 0 and ls["control_group"] == "nt_cells" and ls["resampling_mechanism"] == "permutations"
                  and ls["assign_method"] == "maximum" and ls["power_median_p"] < 1e-6 and ls["calibration"]["n_false_discoveries"] <= 1,
                  f"low-MOI run: nt_cells control, permutations, power p {ls['power_median_p']:.1e}, "
                  f"{ls['calibration']['n_false_discoveries']} calibration false discoveries")

            print("\nR runner + comparison")
            rc = cmd_r_runner(argparse.Namespace(step="inference", mudata=str(H["path"]), out=None, side="left",
                                                 formula="~ log(response_n_nonzero) + log(response_n_umis)", label="st_sceptre_r"))
            rs_ = json.loads((latest("st_sceptre_r") / "summary.json").read_text())
            check(rc == 0 and Path(rs_["script"]).exists()
                  and "inference_sceptre(m, side = \"left\"" in Path(rs_["script"]).read_text(),
                  f"r-runner writes the sceptreIGVF script and exits 0 (status {rs_['status']})")
            a = idir / "test_results.tsv"
            rc = cmd_compare_results(argparse.Namespace(a=str(a), b=str(a), keys=["intended_target_name", "gene_id"], alpha=0.1,
                                                        label="st_sceptre_cmp"))
            cs2 = json.loads((latest("st_sceptre_cmp") / "summary.json").read_text())
            check(rc == 0 and abs(cs2["spearman_neglog10p"] - 1) < 1e-9 and cs2["significance_agreement"] == 1.0,
                  "compare-results: identical tables agree perfectly")
    finally:
        for p in (set(OUT_ROOT.glob("*")) if OUT_ROOT.exists() else set()) - before:
            shutil.rmtree(p, ignore_errors=True)
    ok = all(c for c, _ in checks)
    print(f"\n({time.time() - t0:.0f} s)")
    print("\nselftest: all checks pass" if ok else f"\nselftest: {sum(1 for c, _ in checks if not c)} FAILED")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "Optional[list[str]]" = None) -> int:
    ap = argparse.ArgumentParser(prog="igvfagent sceptre-igvf", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def base(p, label: str, mudata: bool = True):
        if mudata:
            p.add_argument("--mudata", required=True, help="IGVF CRISPR MuData (.h5mu): mod/gene, mod/guide, uns pairs_to_test")
            p.add_argument("--moi", choices=["low", "high"], help="override guide uns['moi']")
        p.add_argument("--label", default=label)

    def test_opts(p):
        p.add_argument("--side", choices=["both", "left", "right"], default="both")
        p.add_argument("--control-group", choices=["default", "complement", "nt_cells"], default="default")
        p.add_argument("--resampling-mechanism", choices=["default", "crt", "permutations"], default="default")
        p.add_argument("--resampling-approximation", choices=["skew_normal", "no_approximation"], default="skew_normal")
        p.add_argument("--b1", type=int, default=499)
        p.add_argument("--b2", type=int, default=4999)
        p.add_argument("--b3", type=int, default=24999)
        p.add_argument("--seed", type=int, default=4)
        p.add_argument("--alpha", type=float, default=0.1, help="BH level (multiple_testing_alpha)")
        p.add_argument("--formula", default=None, help="formula_object, e.g. '~ log(response_n_umis) + batch' ('default' = the upstream default)")
        p.add_argument("--grna-integration-strategy", choices=["union", "singleton", "bonferroni"], default="union")

    def workflow_opts(p):
        p.add_argument("--method", choices=["default", "mixture", "thresholding", "maximum"], default="default",
                       help="gRNA assignment (default: mixture high MOI, maximum low MOI)")
        p.add_argument("--threshold", type=float, default=None, help="thresholding UMI cutoff (default 5)")
        p.add_argument("--n-nonzero-trt-thresh", type=int, default=7)
        p.add_argument("--n-nonzero-cntrl-thresh", type=int, default=7)
        p.add_argument("--positive-control-pairs", help="TSV grna_target/response_id (or intended_target_name/gene_id)")
        p.add_argument("--discovery-pairs", help="TSV of discovery pairs (default: uns pairs_to_test)")
        p.add_argument("--pair-set", choices=["cis", "trans"], help="construct discovery pairs instead of reading them")
        p.add_argument("--cis-distance", type=int, default=500_000)

    p = sub.add_parser("convert", help="MuData -> sceptre object directory (convert_mudata_to_sceptre_object)")
    base(p, "convert")
    p.add_argument("--remove-collinear-covariates", action="store_true")
    p.add_argument("--collinear-mode", choices=["all-or-nothing", "greedy"], default="all-or-nothing")
    p.add_argument("--target-col", default="intended_target_name")
    p.set_defaults(func=cmd_convert)

    p = sub.add_parser("assign-guides", help="assign_grnas_sceptre: write the guide_assignment layer")
    base(p, "assign_guides")
    p.add_argument("--method", choices=["mixture", "thresholding", "maximum", "default"], default="mixture")
    p.add_argument("--threshold", type=float, default=None)
    p.add_argument("--probability-threshold", type=float, default=0.8)
    p.add_argument("--n-em-rep", type=int, default=5)
    p.add_argument("--collinear-mode", choices=["all-or-nothing", "greedy"], default="all-or-nothing")
    p.add_argument("--out", help="output h5mu (default: run dir)")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_assign_guides)

    p = sub.add_parser("qc", help="cell-wise QC summary (sceptre::run_qc)")
    base(p, "qc")
    p.add_argument("--method", choices=["default", "mixture", "thresholding", "maximum"], default="default")
    p.add_argument("--threshold", type=float, default=None)
    p.add_argument("--p-mito-threshold", type=float, default=0.2)
    p.add_argument("--response-n-umis-range", type=float, nargs=2, default=[0.01, 0.99])
    p.add_argument("--response-n-nonzero-range", type=float, nargs=2, default=[0.01, 0.99])
    p.set_defaults(func=cmd_qc)

    p = sub.add_parser("inference", help="inference_sceptre: uns test_results (or per_element/per_guide results)")
    base(p, "inference")
    test_opts(p)
    p.add_argument("--mode", choices=["sceptreIGVF", "pipeline"], default="sceptreIGVF")
    p.add_argument("--scope", choices=["cis", "trans"], default=None, help="pipeline mode: prefix uns keys cis_/trans_")
    p.add_argument("--out", help="output h5mu (default: run dir)")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_inference)

    p = sub.add_parser("calibration-check", help="negative-control (NT) pairs")
    base(p, "calibration_check")
    test_opts(p)
    workflow_opts(p)
    p.add_argument("--n-calibration-pairs", type=int, default=None)
    p.add_argument("--group-size", type=int, default=None, help="undercover NT group size (default: median targeting group size)")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_calibration_check)

    p = sub.add_parser("power-check", help="positive-control pairs")
    base(p, "power_check")
    test_opts(p)
    workflow_opts(p)
    p.set_defaults(func=cmd_power_check)

    p = sub.add_parser("run", help="whole SCEPTRE workflow -> MuData (sceptre_object_to_mudata)")
    base(p, "sceptre_run")
    test_opts(p)
    workflow_opts(p)
    p.add_argument("--p-mito-threshold", type=float, default=0.2)
    p.add_argument("--n-calibration-pairs", type=int, default=None)
    p.add_argument("--skip-calibration", action="store_true")
    p.add_argument("--out", help="output h5mu (default: run dir)")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("export-mudata", help="sceptre_object_to_mudata_inputs_outputs + save_mudata_list")
    base(p, "export_mudata")
    test_opts(p)
    workflow_opts(p)
    p.add_argument("--num-discovery-pairs", type=int, default=100)
    p.add_argument("--out-dir", help="directory for guide_assignment/ and inference/ (default: run dir)")
    p.add_argument("--prefix", default="")
    p.add_argument("--gene-info", help="TSV of gene rowData (index = gene id)")
    p.add_argument("--guide-capture-method", default=None)
    p.set_defaults(func=cmd_export_mudata)

    p = sub.add_parser("r-runner", help="run the real sceptreIGVF R package (prints the command if R is absent)")
    base(p, "r_runner")
    p.add_argument("--step", choices=["assign", "inference"], required=True)
    p.add_argument("--side", choices=["both", "left", "right"], default=None)
    p.add_argument("--formula", default=None)
    p.add_argument("--out")
    p.set_defaults(func=cmd_r_runner)

    p = sub.add_parser("compare-results", help="compare two test_results tables (e.g. Python vs R)")
    base(p, "compare", mudata=False)
    p.add_argument("--a", required=True)
    p.add_argument("--b", required=True)
    p.add_argument("--keys", nargs="+", default=["intended_target_name", "gene_id"])
    p.add_argument("--alpha", type=float, default=0.1)
    p.set_defaults(func=cmd_compare_results)

    p = sub.add_parser("selftest", help="synthetic screens with planted effects; all subcommands asserted")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
