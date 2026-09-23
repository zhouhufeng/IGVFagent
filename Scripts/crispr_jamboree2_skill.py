#!/usr/bin/env python3
"""2nd IGVF CRISPR Jamboree single-cell Perturb-seq analyses (port of IGVF-CRISPR/CRISPR-JAMBOREE).

Port of https://github.com/IGVF-CRISPR/CRISPR-JAMBOREE (2024 "2nd IGVF CRISPR
Jamboree", No LICENSE file; notebooks + R/Python scripts + Nextflow modules).
Every notebook and script under single-cell/ was read and the analysis it
performs was re-derived in Python (numpy / scipy / pandas / anndata); no code
was copied.  Relationship: port (definitions, thresholds, column names and
MuData schema reproduced; R and Stan models re-implemented).  The bulk/ tree
of the upstream holds only empty placeholder files.  Large example MuDatas and
the guide-counting FASTQ tarball (single-cell/**/*.h5mu, example_data/
data.tar.gz, ~70 MB total) are not vendored; point the subcommands at them.

The shared data format (single-cell/inference/*.pdf, scanpy-demo-python.ipynb):
a MuData with modalities `gene` (.X RNA UMI counts; obs total_gene_umis,
num_expressed_genes; var symbol, gene_chr, gene_start, gene_end) and `guide`
(.X gRNA UMI counts, layers['guide_assignment'] binary; var targeting "TRUE" /
"FALSE", intended_target_name, intended_target_chr/start/end; uns moi
low|high, capture_method CROP-seq|direct capture), global obs (prep_batch ...)
and uns['pairs_to_test'] (intended_target_name, gene_id[, pair_type]).  An
inference module adds uns['test_results'] = pairs_to_test + p_value / log2_fc /
posterior_probability.

Definitions, exactly as upstream computes them
  control cells       low MOI: cells with >= 1 assigned non-targeting gRNA
                      (targeting == FALSE); high MOI: cells without a gRNA for
                      the element (wilcoxon-test.R, negative-binomial.R).
  treatment cells     cells with max(guide_assignment) over the element's
                      gRNAs == 1 ("union" of the element's guides).
  wilcoxon (R)        stats::wilcox.test(treatment counts, control counts),
                      side left|right|both -> less|greater|two.sided; exact
                      for n < 50 without ties, else normal approximation with
                      continuity and tie correction.  Column p_value.
  negbinom (R)        MASS::glm.nb(y ~ element_targeted + offset(log
                      total_gene_umis)), Wald p of element_targeted, low MOI
                      restricted to treatment + NT cells.  Upstream stores
                      `l2fc` = exp(b0 + b1) / exp(b0) = exp(b1), a fold change
                      not a log2 fold change; the port writes log2_fc =
                      b1 / ln 2 and keeps `l2fc` only with --upstream-compat.
  scanpy (Python)     run_scanpy_test.py: log1p of raw counts (no size-factor
                      normalisation), sc.tl.rank_genes_groups(groupby
                      guide_status, groups=[present]) vs all other cells;
                      methods wilcoxon (rank-sum z, no tie correction,
                      p = 2 sf|z|), t-test (Welch), t-test_overestim_var
                      (rest variance divided by n_group); logfoldchanges =
                      log2((expm1(mean_g)+1e-9)/(expm1(mean_rest)+1e-9)).
  sceptre (R pkg)     sceptreIGVF::inference_sceptre defaults: formula
                      ~ log(response_n_nonzero) + log(response_n_umis), side
                      both, grna_integration_strategy union, resampling
                      conditional randomisation (high MOI) / permutations (low
                      MOI), skew-normal approximation of the resampled NB
                      score statistics.  Re-implemented (NB GLM score test +
                      CRT with a logistic propensity model); validated in the
                      self-test, not bit-identical to the R package.
  perturbo            Bayesian NB model (pyro); run when the `perturbo`
                      package is importable (run-perturbo.py settings: batch
                      prep_batch, library size total_gene_umis, covariate
                      total_guide_umis, 20 epochs, lr 0.01, batch 128;
                      q_value -> posterior_probability, loc/ln2 -> log2_fc).
  guide assignment    igvf_guide_assignment.py: --umi-threshold assigns
                      count >= t (default 5); --cleanser fits the CLEANSER
                      zero-truncated mixture per gRNA (CROP-seq: Poisson
                      background + NB signal; direct capture: NB + NB; the Stan
                      priors of cs-/dc-guide-mixture.stan) and stores the
                      posterior PZi, or 1 where PZi >= t.  The port uses the
                      MAP fit instead of the median over 50 NUTS draws.
                      sceptre mixture: two-component Poisson GLM per gRNA,
                      assigned when posterior > 0.8.
  simulation          1_create_sim_mudata.R: pairs_to_test = pairs to perturb
                      with effect_size (default 0.5).  2_fit_negative_binomial.R:
                      DESeq2 estimateSizeFactors(type poscounts) and
                      estimateDispersions(fitType parametric), rowData mean
                      (baseMean), dispersion, disp_outlier_deseq2, colData
                      size_factors.  3_simulate_count_data.R: guide effect ~
                      N(effect_size, guide_var) (guide_var 0.1 in the script;
                      1 for null pairs), one random guide per cell-gene when
                      several apply, mu = mean * size_factor * effect,
                      counts ~ rnbinom(mu, size = 1/dispersion).  Upstream
                      writes the *unsimulated* object (`writeH5MU(mu, ...)`
                      instead of `output`); the port writes the simulation and
                      reproduces the bug with --upstream-compat.
  evaluation          eval_and_clustergram.ipynb: true label effect_size != 1;
                      AUPRC = auc(recall, precision) (trapezoid over the
                      sklearn PR curve), AUROC; continuous: Pearson r and MSE
                      of log2(effect_size) vs estimated log2_fc.  Upstream
                      scores pairs by the raw p_value (larger p = more
                      positive, which inverts the curves); the port scores by
                      -log10(p) and keeps the raw score with --upstream-compat.
                      Clustergram: per gRNA (guide .X != 0) log2(mean gene
                      expression with / without the gRNA), non-finite -> 0,
                      average-linkage Euclidean clustering of both axes.
  volcano / IGV       volcano_and_igv.ipynb: down = log2_fc <= -2 and
                      sig <= 0.01, up = log2_fc >= 2 and sig <= 0.01 (sig =
                      posterior_probability or p_value); IGV bedgraph for
                      promoter pairs (intended_target_name == gene_id) and
                      bedpe element -> gene links, coordinates from gene.var
                      and guide.var.
  network             network.ipynb / evaluation.py: edges element -> gene
                      with |log2_fc| >= min_weight (0.01), per central node,
                      node size -log10(sig); mean_nontargeting_expression per
                      gene over cells with >= 1 non-targeting and no targeting
                      gRNA (guide .X counts).
  guide reference     CountGuides.ipynb / parsing_guide_medata.py: pseudo
                      chromosome per protospacer (seqid chr_<sgRNA_ID>_peseudo_
                      chr_<seq>), BED12 (0, len, thick 1..len, 0,0,255),
                      isoforms.txt (seqid, protospacer), bed2gtf GTF, pseudo
                      genome FASTA; protospacer FASTA with index
                      chr_start_end_strand_protospacer_target (targeting) or
                      type_protospacer; optional 'NNN' + reverse complement.
  guide counting      starSoloGuide.nf: CB = R1[1..16], UMI = R1[17..28],
                      whitelist 1MM correction, protospacer found in R2, UMI
                      de-duplication with 1-mismatch collapsing
                      (soloUMIdedup 1MM_CR); soloOutGuideParsing.py maps
                      feature names (last '_' field = protospacer) back to
                      sgRNA_ID.  STAR itself is optional (command printed).
  seqspec index       seqSpecExtract.ipynb: `seqspec index -t kb|starsolo -m
                      <modality> -r <reads>` positions of barcode / umi / cdna.

Validation against the upstream's own output MuDatas (Gasperini subset,
9704 cells, 110 pairs): scanpy wilcoxon / t-test / t-test_overestim_var
p-values and log2_fc identical (max |diff| 3e-14 and 1e-7); negbinom p and
fold change within 3e-5 relative; R wilcoxon identical (1e-16) with side left,
the setting the upstream demo outputs were produced with.  Simulation: size
factors and baseMean identical (1e-14); dispersions within 0.3% for genes with
mean >= 1, while the upstream values for very low-count genes follow a fitted
trend the parametric port does not reproduce (DESeq2 fitType fallback).

Subcommands
  inspect         MuData summary + schema check against the jamboree input /
                  output specification (eval_setup.ipynb).
  guide-reference guide metadata -> pseudo-genome FASTA, BED12, isoforms, GTF,
                  protospacer FASTA; STAR genomeGenerate when STAR is on PATH.
  count-guides    FASTQ R1/R2 + whitelist + guide metadata -> cells x guides
                  UMI matrix (h5ad + MatrixMarket), or parse a STARsolo
                  Solo.out directory (--solo-out).
  seqspec-index   seqspec YAML -> kb `-x` technology string / STARsolo flags.
  assign-guides   guide_assignment layer by UMI threshold, CLEANSER or the
                  sceptre mixture.
  infer           test pairs_to_test with wilcoxon, negbinom, scanpy-wilcoxon,
                  scanpy-t-test, scanpy-t-test-overestim-var, sceptre, perturbo.
  simulate        DESeq2-style NB fit + simulated Perturb-seq counts with
                  planted effects (optionally extra null pairs).
  evaluate        AUPRC / AUROC / Pearson / MSE of test results against
                  planted effects or control pair types, + clustergram.
  volcano         volcano tables + figure and IGV bedpe / bedgraph tracks.
  network         element -> gene network edges + figure, NT mean expression.
  run             assign (if needed) -> infer -> evaluate -> volcano/IGV ->
                  network, one report.
  selftest        synthetic MuData + FASTQs with planted signal; all
                  subcommands asserted.

Output: Docs/CRISPRJamboree2/<timestamp>_<label>/.

Usage:
    igvfagent crispr-jamboree2 inspect --mudata gasperini_inference_input.h5mu
    igvfagent crispr-jamboree2 guide-reference --guides guide_metadata.tsv --label guides
    igvfagent crispr-jamboree2 count-guides --r1 R1.fastq.gz --r2 R2.fastq.gz --whitelist 737K-august-2016.txt --guides guide_metadata.tsv
    igvfagent crispr-jamboree2 seqspec-index --spec spec.yaml --modality crispr --reads R1,R2 --tool kb
    igvfagent crispr-jamboree2 assign-guides --mudata in.h5mu --method umi-threshold --threshold 5 --out assigned.h5mu
    igvfagent crispr-jamboree2 infer --mudata gasperini_inference_input.h5mu --methods wilcoxon scanpy-wilcoxon sceptre
    igvfagent crispr-jamboree2 simulate --mudata gasperini_inference_input.h5mu --perturb ENSG00000136856:candidate_enh_3:0.5 --n-null-pairs 50
    igvfagent crispr-jamboree2 evaluate --mudata inference_output.h5mu
    igvfagent crispr-jamboree2 run --mudata gasperini_inference_input.h5mu --label gasperini
    igvfagent crispr-jamboree2 selftest --no-plots
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
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1]).resolve()
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
OUT_ROOT = DOCS_DIR / "CRISPRJamboree2"

UPSTREAM_REPO = "IGVF-CRISPR/CRISPR-JAMBOREE"
UPSTREAM_COMMIT = "56d31fd4bc56966500c30b51b1c9b3923cf5e1b0"
INFER_METHODS = ["wilcoxon", "negbinom", "scanpy-wilcoxon", "scanpy-t-test", "scanpy-t-test-overestim-var", "sceptre", "perturbo"]

INK, INK2, AXIS, SURFACE = "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"
BLUE, ORANGE, GREEN, RED, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#e34948", "#96a0b3"

log = logging.getLogger("crispr_jamboree2")


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"crispr_jamboree2_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
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


def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


# ---------------------------------------------------------------------------
# MuData I/O (h5py + anndata; the `mudata` package is optional)
# ---------------------------------------------------------------------------

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
# Statistics
# ---------------------------------------------------------------------------

def bh(p):
    np = _np()
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    n = int(ok.sum())
    if n == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


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


def r_wilcox(x, y, side: str = "both") -> float:
    """stats::wilcox.test(x, y, alternative): exact without ties when n < 50, else normal approx + continuity."""
    np = _np()
    from scipy.stats import mannwhitneyu  # type: ignore
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    alt = {"left": "less", "right": "greater", "both": "two-sided"}[side]
    allv = np.concatenate([x, y])
    ties = len(np.unique(allv)) < len(allv)
    exact = (len(x) < 50 and len(y) < 50 and not ties)
    if not exact and np.all(allv == allv[0]):
        return 1.0
    return float(mannwhitneyu(x, y, alternative=alt, method="exact" if exact else "asymptotic", use_continuity=True).pvalue)


def scanpy_rank_test(Xlog, group, method: str = "wilcoxon"):
    """sc.tl.rank_genes_groups(groups=['present'], reference='rest') on a dense log1p matrix.

    Returns (scores, pvals, logfoldchanges) per column.
    """
    np = _np()
    from scipy import stats  # type: ignore
    g = np.asarray(group, bool)
    n1 = int(g.sum())
    n2 = int((~g).sum())
    m1 = Xlog[g].mean(axis=0)
    m2 = Xlog[~g].mean(axis=0)
    if method == "wilcoxon":
        ranks = stats.rankdata(Xlog, axis=0)
        s = ranks[g].sum(axis=0)
        scores = (s - n1 * (n1 + n2 + 1) / 2.0) / math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
        scores = np.where(np.isfinite(scores), scores, 0.0)
        pvals = 2 * stats.norm.sf(np.abs(scores))
    else:
        v1 = Xlog[g].var(axis=0, ddof=1) if n1 > 1 else np.zeros(Xlog.shape[1])
        v2 = Xlog[~g].var(axis=0, ddof=1) if n2 > 1 else np.zeros(Xlog.shape[1])
        ns_rest = n2 if method == "t-test" else n1
        with np.errstate(divide="ignore", invalid="ignore"):
            a, b = v1 / n1, v2 / ns_rest
            scores = (m1 - m2) / np.sqrt(a + b)
            df = (a + b) ** 2 / (a ** 2 / (n1 - 1) + b ** 2 / (ns_rest - 1))
            pvals = 2 * stats.t.sf(np.abs(scores), df)
        scores = np.where(np.isfinite(scores), scores, 0.0)
        pvals = np.where(np.isfinite(pvals), pvals, 1.0)
    lfc = np.log2((np.expm1(m1) + 1e-9) / (np.expm1(m2) + 1e-9))
    return scores, pvals, lfc


def auc_trapezoid(x, y) -> float:
    np = _np()
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    order = np.lexsort((y, x))
    return float(np.trapz(y[order], x[order])) if hasattr(np, "trapz") else float(np.trapezoid(y[order], x[order]))


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


# ---------------------------------------------------------------------------
# SCEPTRE-style test (NB score statistic + conditional randomisation / permutation)
# ---------------------------------------------------------------------------

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
# Inference methods (wilcoxon-test.R, negative-binomial.R, run_scanpy_test.py, sceptre, PerTurbo)
# ---------------------------------------------------------------------------

def _pairs(md: MuDataLite):
    p = uns_df(md, "pairs_to_test")
    p["gene_id"] = p["gene_id"].astype(str)
    p["intended_target_name"] = p["intended_target_name"].astype(str)
    return p


def infer_wilcoxon(md: MuDataLite, pairs, side: str = "both", moi: Optional[str] = None):
    np = _np()
    moi = moi_of(md, moi)
    ctrl_nt = nt_cells(md.guide)
    res = pairs.copy()
    pv = []
    for _, r in pairs.iterrows():
        trt = element_indicator(md.guide, r["intended_target_name"])
        ctrl = ctrl_nt if moi == "low" else ~trt
        y = gene_counts(md.gene, [r["gene_id"]])[:, 0]
        pv.append(r_wilcox(y[trt], y[ctrl], side))
    res["p_value"] = pv
    return res


def infer_negbinom(md: MuDataLite, pairs, moi: Optional[str] = None, compat: bool = False):
    np = _np()
    moi = moi_of(md, moi)
    ctrl_nt = nt_cells(md.guide)
    umis = total_umis(md.gene)
    res = pairs.copy()
    pv, l2, fc = [], [], []
    for _, r in pairs.iterrows():
        trt = element_indicator(md.guide, r["intended_target_name"])
        keep = (trt | ctrl_nt) if moi == "low" else np.ones(md.gene.n_obs, dtype=bool)
        keep = keep & (umis > 0)
        y = gene_counts(md.gene, [r["gene_id"]])[keep, 0]
        x = trt[keep].astype(float)
        if x.sum() == 0 or x.sum() == len(x) or y.sum() == 0:
            pv.append(np.nan); l2.append(np.nan); fc.append(np.nan)
            continue
        fit = glm_nb(y, np.column_stack([np.ones(len(y)), x]), offset=np.log(umis[keep]))
        pv.append(float(fit["p"][1]))
        l2.append(float(fit["beta"][1] / math.log(2)))
        fc.append(float(math.exp(fit["beta"][1])))
    res["p_value"] = pv
    if compat:
        res["l2fc"] = fc
    else:
        res["log2_fc"] = l2
    return res


def infer_scanpy(md: MuDataLite, pairs, method: str = "wilcoxon"):
    np, pd = _np(), _pd()
    rows = []
    for target in pairs["intended_target_name"].unique():
        sub = pairs[pairs["intended_target_name"] == target]
        genes = list(dict.fromkeys(sub["gene_id"]))
        trt = element_indicator(md.guide, target)
        Xlog = np.log1p(gene_counts(md.gene, genes))
        if trt.sum() == 0 or (~trt).sum() == 0:
            s = p = l = np.full(len(genes), np.nan)
        else:
            s, p, l = scanpy_rank_test(Xlog, trt, method)
        for g, pp, ll in zip(genes, p, l):
            rows.append({"gene_id": g, "intended_target_name": target, "p_value": float(pp), "log2_fc": float(ll)})
    tr = pd.DataFrame(rows)
    return tr.merge(pairs, how="right", on=["intended_target_name", "gene_id"])


def infer_perturbo(md_path: Path, out_h5mu: Path):
    try:
        import mudata  # type: ignore  # noqa: F401
        import perturbo  # type: ignore  # noqa: F401
    except ImportError:
        cmd = f"python run-perturbo.py {md_path} {out_h5mu}  # needs `pip install perturbo mudata` (pyro, scvi-tools)"
        print("perturbo: package not importable; PerTurbo is a pyro variational model and is not re-implemented here.")
        print(f"Command: {cmd}")
        return None
    import mudata as mdmod  # type: ignore
    import numpy as np  # type: ignore
    import pandas as pd  # type: ignore
    import perturbo  # type: ignore
    m = mdmod.read_h5mu(str(md_path))
    m["gene"].obs = (m.obs.join(m["gene"].obs).join(m["guide"].obs)
                     .assign(log1p_total_guide_umis=lambda x: np.log1p(x["total_guide_umis"])))
    m["guide"].X = m["guide"].layers["guide_assignment"]
    ptt = pd.DataFrame(m.uns["pairs_to_test"])
    names = sorted(pd.unique(ptt["intended_target_name"]))
    m.uns["intended_target_names"] = names
    m["gene"].varm["intended_targets"] = (ptt.assign(value=1).pivot(index="gene_id", columns="intended_target_name", values="value")
                                          .reindex(m["gene"].var_names).fillna(0))
    m["guide"].varm["intended_targets"] = pd.get_dummies(m["guide"].var["intended_target_name"]).astype(float)[names]
    perturbo.PERTURBO.setup_mudata(m, batch_key="prep_batch", library_size_key="total_gene_umis",
                                   continuous_covariates_keys=["total_guide_umis"], guide_by_element_key="intended_targets",
                                   gene_by_element_key="intended_targets", modalities={"rna_layer": "gene", "perturbation_layer": "guide"})
    model = perturbo.PERTURBO(m, likelihood="nb")
    model.train(20, lr=0.01, batch_size=128)
    eff = (model.get_element_effects().rename(columns={"element": "intended_target_name", "gene": "gene_id", "q_value": "posterior_probability"})
           .assign(log2_fc=lambda x: x["loc"] / np.log(2)).merge(ptt))
    return eff[["intended_target_name", "gene_id", "posterior_probability", "log2_fc"]]


def run_inference(md: MuDataLite, method: str, side: str = "both", moi: Optional[str] = None, B: int = 499, seed: int = 4,
                  compat: bool = False, md_path: Optional[Path] = None, work: Optional[Path] = None):
    pairs = _pairs(md)
    if method == "wilcoxon":
        return infer_wilcoxon(md, pairs, side, moi)
    if method == "negbinom":
        return infer_negbinom(md, pairs, moi, compat)
    if method.startswith("scanpy-"):
        return infer_scanpy(md, pairs, method[len("scanpy-"):])
    if method == "sceptre":
        res = sceptre_test_pairs(md, pairs, side=side, B=B, moi=moi, seed=seed)
        return pairs.merge(res, how="left", on=["intended_target_name", "gene_id"])
    if method == "perturbo":
        if md_path is None:
            md_path = (work or Path(".")) / "perturbo_input.h5mu"
            write_h5mu(md, md_path)
        return infer_perturbo(md_path, (work or Path(".")) / "perturbo_output.h5mu")
    raise SystemExit(f"unknown method {method}; choose from {INFER_METHODS}")


def sig_column(tr) -> str:
    for c in ("p_value", "posterior_probability", "q_value"):
        if c in tr.columns:
            return c
    raise SystemExit("test_results has no p_value / posterior_probability column")


def lfc_column(tr) -> Optional[str]:
    for c in ("log2_fc", "l2fc"):
        if c in tr.columns:
            return c
    return None


def cmd_infer(args: argparse.Namespace) -> int:
    np = _np()
    md = read_h5mu(Path(args.mudata))
    d = run_dir(args.label)
    methods = args.methods or ["wilcoxon", "negbinom", "scanpy-wilcoxon", "sceptre"]
    rows, results = [], {}
    for m in methods:
        t0 = time.time()
        tr = run_inference(md, m, side=args.side, moi=args.moi, B=args.resamples, seed=args.seed, compat=args.upstream_compat,
                           md_path=Path(args.mudata), work=d)
        if tr is None:
            rows.append([m, "skipped (not installed)", "", ""])
            continue
        results[m] = tr
        write_tsv(tr, d / f"test_results.{m}.tsv")
        sc = sig_column(tr)
        n_sig = int((tr[sc] < 0.05).sum())
        rows.append([m, len(tr), n_sig, f"{time.time() - t0:.1f}s"])
    if results and args.out:
        first = next(iter(results))
        md.uns["test_results"] = results[first]
        write_h5mu(md, Path(args.out))
    write_report(d, "CRISPR Jamboree 2 — Perturb-seq inference",
                 [f"Input: `{args.mudata}` ({md.gene.n_obs} cells, {md.gene.n_vars} genes, {md.guide.n_vars} gRNAs, MOI {moi_of(md, args.moi)}).",
                  md_table(["method", "pairs tested", "p < 0.05", "time"], rows)],
                 {"mudata": args.mudata, "methods": {r[0]: {"n": r[1], "n_p_lt_0.05": r[2]} for r in rows}, "run_dir": d})
    return 0


# ---------------------------------------------------------------------------
# Guide assignment (igvf_guide_assignment.py, CLEANSER Stan models, sceptre mixture)
# ---------------------------------------------------------------------------

def _nb2_logpmf(x, mu, phi):
    np = _np()
    from scipy.special import gammaln  # type: ignore
    return gammaln(x + phi) - gammaln(phi) - gammaln(x + 1) + phi * np.log(phi / (phi + mu)) + x * np.log(mu / (phi + mu))


def _pois_logpmf(x, lam):
    np = _np()
    from scipy.special import gammaln  # type: ignore
    return x * np.log(lam) - lam - gammaln(x + 1)


def _lognormal_lp(v, mu, sigma):
    return -math.log(v) - (math.log(v) - mu) ** 2 / (2 * sigma ** 2)


def cleanser_fit(x, capture: str = "CROP-seq", L=None) -> "tuple[dict, Any]":
    """MAP fit of the CLEANSER zero-truncated mixture (cs-/dc-guide-mixture.stan); returns (params, PZi)."""
    np = _np()
    from scipy.optimize import minimize  # type: ignore
    from scipy.special import logsumexp  # type: ignore
    x = np.asarray(x, float)
    L = np.ones_like(x) if L is None else np.asarray(L, float)
    dc = capture.lower().startswith("direct")
    sig = lambda u, lo, hi: lo + (hi - lo) / (1 + math.exp(-u))  # noqa: E731

    def unpack(th):
        if dc:
            n_mean = 1e-6 + math.exp(th[0]); n_disp = 0.1 + math.exp(th[1]); r = sig(th[2], 1e-6, 0.1)
            mean = 5 + math.exp(th[3]); disp = 0.1 + math.exp(th[4])
            return {"n_nbMean": n_mean, "n_nbDisp": n_disp, "r": r, "nbMean": mean, "nbDisp": disp}
        lam = sig(th[0], 1e-6, 5.0); r = sig(th[1], 1e-6, 0.1); mean = lam + math.exp(th[2]); disp = sig(th[3], 1e-6, 1 - 1e-9)
        return {"lambda": lam, "r": r, "nbMean": mean, "nbDisp": disp}

    def comps(p):
        if dc:
            lp0 = math.log(1 - p["r"]) + _nb2_logpmf(x, p["n_nbMean"] * L, p["n_nbDisp"])
            z0 = math.log(1 - p["r"]) + p["n_nbDisp"] * (np.log(p["n_nbDisp"]) - np.log(p["n_nbDisp"] + p["n_nbMean"] * L))
        else:
            lp0 = math.log(1 - p["r"]) + _pois_logpmf(x, p["lambda"])
            z0 = math.log(1 - p["r"]) - p["lambda"] + 0 * x
        lp1 = math.log(p["r"]) + _nb2_logpmf(x, p["nbMean"] * L, p["nbDisp"])
        z1 = math.log(p["r"]) + p["nbDisp"] * (np.log(p["nbDisp"]) - np.log(p["nbDisp"] + p["nbMean"] * L))
        return lp0, lp1, z0, z1

    def neg_lp(th):
        try:
            p = unpack(th)
        except OverflowError:
            return 1e30
        lp0, lp1, z0, z1 = comps(p)
        pz = np.exp(logsumexp(np.vstack([np.broadcast_to(z0, x.shape), np.broadcast_to(z1, x.shape)]), axis=0))
        ll = np.sum(logsumexp(np.vstack([lp0, lp1]), axis=0) - np.log(np.clip(1 - pz, 1e-300, None)))
        pr = 9 * math.log(1 - p["r"])
        if dc:
            pr += (_lognormal_lp(p["nbDisp"], math.log(3), math.log(2)) + _lognormal_lp(p["n_nbMean"], math.log(1.075), math.log(1.1))
                   + _lognormal_lp(p["n_nbDisp"], math.log(1.1), math.log(1.1)) + _lognormal_lp(p["nbMean"], math.log(100), math.log(3)))
        else:
            pr += (_lognormal_lp(p["lambda"], math.log(1.1), math.log(1.1)) + _lognormal_lp(p["nbMean"], math.log(100), math.log(100))
                   + 9 * math.log(1 - p["nbDisp"]))
        v = -(ll + pr)
        return v if np.isfinite(v) else 1e30

    hi = np.percentile(x, 90) if len(x) else 10.0
    best = None
    for mean0 in (max(hi, 6.0), max(np.median(x), 6.0), 50.0):
        for r0 in (0.05, 0.5):
            if dc:
                th0 = [math.log(0.1), math.log(1.0), 0.0 if r0 == 0.05 else 2.0, math.log(max(mean0 - 5, 1)), math.log(2.9)]
            else:
                th0 = [0.0, 0.0 if r0 == 0.05 else 2.0, math.log(max(mean0 - 1, 1)), -1.0]
            res = minimize(neg_lp, np.array(th0, float), method="Nelder-Mead", options={"maxiter": 4000, "xatol": 1e-6, "fatol": 1e-8})
            if best is None or res.fun < best.fun:
                best = res
    p = unpack(best.x)
    lp0, lp1, _, _ = comps(p)
    pzi = np.exp(lp1 - logsumexp(np.vstack([lp0, lp1]), axis=0))
    p["neg_log_posterior"] = float(best.fun)
    return p, pzi


def assign_cleanser(guide, threshold: Optional[float] = None, capture: Optional[str] = None, pooled: bool = False):
    """Per-gRNA CLEANSER posteriors over the nonzero counts (pooled=True: one fit over all gRNAs, jamboree3 cleanser.py)."""
    np = _np()
    import scipy.sparse as sp  # type: ignore
    cap = capture or _uns_scalar(guide, "capture_method", "CROP-seq")
    C = sp.csc_matrix(guide.X)
    out = sp.lil_matrix(C.shape)
    params = {}
    if pooled:
        coo = sp.coo_matrix(guide.X)
        order = np.lexsort((coo.row, coo.col))
        rows, cols, vals = coo.row[order], coo.col[order], coo.data[order]
        p, pzi = cleanser_fit(vals, cap)
        params["pooled"] = p
        for g in np.unique(cols):
            sel = np.where(cols == g)[0]
            # upstream indexes the pooled PZi by the within-guide position i (cell_info enumerate), reproduced here
            for i, k in enumerate(sel):
                v = pzi[i]
                if threshold is None:
                    out[rows[k], g] = v
                elif v >= threshold:
                    out[rows[k], g] = 1
        return out.tocsr(), params
    for g in range(C.shape[1]):
        col = C.getcol(g)
        idx = col.indices
        vals = col.data
        if len(vals) < 3:
            continue
        p, pzi = cleanser_fit(vals, cap)
        params[str(guide.var_names[g])] = p
        for i, v in zip(idx, pzi):
            if threshold is None:
                if v > 0:
                    out[i, g] = v
            elif v >= threshold:
                out[i, g] = 1
    return out.tocsr(), params


def assign_threshold(guide, threshold: Optional[float] = None):
    import scipy.sparse as sp  # type: ignore
    t = 5 if threshold is None else threshold
    X = sp.csr_matrix(guide.X)
    A = (X >= t).astype(float)
    A.eliminate_zeros()
    return sp.csr_matrix(A)


def assign_sceptre_mixture(guide, gene=None, prob_threshold: float = 0.8, maxit: int = 50):
    """Two-component Poisson GLM mixture per gRNA: log mu = b0 + b*log(response_n_umis) [+ beta if perturbed]."""
    np = _np()
    import scipy.sparse as sp  # type: ignore
    C = sp.csc_matrix(guide.X)
    n = C.shape[0]
    if gene is not None:
        cov = np.log(np.maximum(np.asarray(gene.X.sum(axis=1)).ravel(), 1))
    else:
        cov = np.log(np.maximum(np.asarray(C.sum(axis=1)).ravel(), 1))
    cov = cov - cov.mean()
    out = sp.lil_matrix(C.shape)
    for g in range(C.shape[1]):
        y = np.asarray(C.getcol(g).todense()).ravel().astype(float)
        if (y > 0).sum() < 2:
            continue
        ly = np.log1p(y)
        lo, hi = np.quantile(ly, 0.5), ly.max()
        post = (np.abs(ly - hi) < np.abs(ly - lo)).astype(float)
        if post.sum() == 0:
            continue
        pi = post.mean()
        for _ in range(maxit):
            # M step: weighted Poisson GLM on stacked rows [background | perturbed]
            Xs = np.column_stack([np.ones(2 * n), np.concatenate([cov, cov]), np.concatenate([np.zeros(n), np.ones(n)])])
            ys = np.concatenate([y, y])
            ws = np.concatenate([1 - post, post]) + 1e-8
            b, _, _ = glm_log_irls(ys, Xs, weights=ws)
            pi = min(max(post.mean(), 1e-6), 1 - 1e-6)
            mu0 = np.exp(b[0] + b[1] * cov)
            mu1 = np.exp(b[0] + b[1] * cov + b[2])
            l0 = math.log(1 - pi) + _pois_logpmf(y, mu0)
            l1 = math.log(pi) + _pois_logpmf(y, mu1)
            new = 1 / (1 + np.exp(np.clip(l0 - l1, -700, 700)))
            if np.max(np.abs(new - post)) < 1e-6:
                post = new
                break
            post = new
        if b[2] <= 0:
            post = 1 - post
        for i in np.where(post > prob_threshold)[0]:
            out[i, g] = 1
    return out.tocsr()


def assign_guides(md: MuDataLite, method: str, threshold: Optional[float] = None, capture: Optional[str] = None,
                  pooled: bool = False):
    if method == "umi-threshold":
        return assign_threshold(md.guide, threshold), {}
    if method == "cleanser":
        return assign_cleanser(md.guide, threshold, capture, pooled)
    if method == "sceptre":
        return assign_sceptre_mixture(md.guide, md.gene if "gene" in md else None), {}
    raise SystemExit(f"unknown assignment method {method}")


def cmd_assign(args: argparse.Namespace) -> int:
    np = _np()
    md = read_h5mu(Path(args.mudata))
    d = run_dir(args.label)
    A, params = assign_guides(md, args.method, args.threshold, args.capture_method)
    md.guide.layers["guide_assignment"] = A
    out = Path(args.out) if args.out else d / "guide_assignment_output.h5mu"
    write_h5mu(md, out)
    Ad = dense(A)
    per_cell = (Ad > 0).sum(axis=1)
    per_guide = (Ad > 0).sum(axis=0)
    pd = _pd()
    tab = pd.DataFrame({"guide_id": md.guide.var_names, "intended_target_name": md.guide.var["intended_target_name"].values,
                        "n_cells_assigned": per_guide})
    write_tsv(tab, d / "cells_per_guide.tsv")
    if params:
        write_json(params, d / "cleanser_params.json")
    write_report(d, "CRISPR Jamboree 2 — guide assignment",
                 [f"Method `{args.method}` (threshold {args.threshold}); {int((per_cell > 0).sum())} of {len(per_cell)} cells carry >= 1 gRNA; "
                  f"mean {per_cell.mean():.2f} gRNAs per cell; values are {'binary' if args.threshold is not None or args.method != 'cleanser' else 'posterior PZi'}.",
                  md_table(["guide", "target", "cells"], tab.values.tolist(), 30)],
                 {"method": args.method, "threshold": args.threshold, "cells_with_guide": int((per_cell > 0).sum()),
                  "mean_guides_per_cell": float(per_cell.mean()), "out": out})
    return 0


# ---------------------------------------------------------------------------
# Guide reference, counting, seqspec (counting-guides/)
# ---------------------------------------------------------------------------

def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def normalise_guide_metadata(df):
    """Accept jamboree xlsx (sgRNA_ID, sgRNA_sequences) or IGVF metadata (guide_id / protospacer ...)."""
    df = df.copy()
    if "protospacer" not in df.columns:
        for c in ("sgRNA_sequences", "spacer", "sequence"):
            if c in df.columns:
                df["protospacer"] = df[c]
                break
    if "sgRNA_ID" not in df.columns:
        for c in ("guide_id", "grna_id", "name", "ID"):
            if c in df.columns:
                df["sgRNA_ID"] = df[c]
                break
    if "protospacer" not in df.columns or "sgRNA_ID" not in df.columns:
        raise SystemExit("guide metadata needs a protospacer / sgRNA_sequences column and an sgRNA_ID / guide_id column")
    df["protospacer"] = df["protospacer"].astype(str).str.upper().str.strip()
    df["sgRNA_ID"] = df["sgRNA_ID"].astype(str)
    return df


def guide_index_name(row) -> str:
    """CountGuides.ipynb get_index_str."""
    if "targeting" in row and truthy_str(row["targeting"]) and all(c in row for c in ("guide_chr", "guide_start", "guide_end", "strand", "intended_target_name")):
        return f"{row['guide_chr']}_{int(row['guide_start'])}_{int(row['guide_end'])}_{row['strand']}_{row['protospacer']}_{row['intended_target_name']}"
    return f"{row.get('type', 'guide')}_{row['protospacer']}"


def build_guide_reference(df, outdir: Path, nnn_revcomp: bool = False) -> dict:
    pd = _pd()
    df = normalise_guide_metadata(df)
    seqs = df["protospacer"].map(lambda s: "NNN" + revcomp(s) if nnn_revcomp else s)
    bed = pd.DataFrame({"seqid": [f"chr_{i}_peseudo_chr_{s}" for i, s in zip(df["sgRNA_ID"], seqs)], "star": 0,
                        "end": seqs.map(len).values, "source": seqs.values, "type": seqs.values, "strand": "+", "thickStart": 1})
    bed["thickEnd"] = bed["end"]
    bed["color"] = "0,0,255"
    bed["blockCount"] = 1
    bed["blockSizes"] = bed["end"]
    bed["blockStarts"] = 0
    outdir.mkdir(parents=True, exist_ok=True)
    p_bed = outdir / "guides.bed"
    bed.to_csv(p_bed, sep="\t", header=False, index=False)
    p_iso = outdir / "isoforms.txt"
    bed[["seqid", "source"]].to_csv(p_iso, sep="\t", header=False, index=False)
    p_gtf = outdir / "guides.gtf"
    with open(p_gtf, "w") as fh:  # bed2gtf --isoforms: gene = seqid, transcript = BED name
        for _, r in bed.iterrows():
            attrs = f'gene_id "{r["seqid"]}"; transcript_id "{r["source"]}";'
            fh.write(f"{r['seqid']}\tbed2gtf\tgene\t1\t{r['end']}\t.\t+\t.\tgene_id \"{r['seqid']}\";\n")
            fh.write(f"{r['seqid']}\tbed2gtf\ttranscript\t1\t{r['end']}\t.\t+\t.\t{attrs}\n")
            fh.write(f"{r['seqid']}\tbed2gtf\texon\t1\t{r['end']}\t.\t+\t.\t{attrs} exon_number \"1\";\n")
    p_fa = outdir / "pseudo_genome.fa"
    p_fa.write_text("\n".join(f">{s}\n{q}" for s, q in zip(bed["seqid"], bed["source"])))
    p_ps = outdir / "protospacers.fa"
    with open(p_ps, "w") as fh:
        for _, r in df.iterrows():
            fh.write(f">{guide_index_name(r)}\n{r['protospacer']}\n")
    for p in (p_bed, p_iso, p_gtf, p_fa, p_ps):
        print(f"Wrote: {p}")
    return {"bed": p_bed, "isoforms": p_iso, "gtf": p_gtf, "pseudo_genome": p_fa, "protospacers": p_ps, "n_guides": len(df)}


def cmd_guide_reference(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    df = read_table(Path(args.guides))
    ref = build_guide_reference(df, d, args.nnn_revcomp)
    star_cmd = (f"STAR --runThreadN 4 --runMode genomeGenerate --genomeDir {d / 'new_ref'} --genomeFastaFiles {ref['pseudo_genome']} "
                f"--genomeSAindexNbases {args.sa_index_nbases} --sjdbGTFfile {ref['gtf']}")
    ran = False
    if shutil.which("STAR") and not args.no_star:
        (d / "new_ref").mkdir(exist_ok=True)
        ran = subprocess.run(star_cmd, shell=True).returncode == 0
    else:
        print(f"STAR not on PATH (or --no-star); command: {star_cmd}")
    write_report(d, "CRISPR Jamboree 2 — guide reference",
                 [f"{ref['n_guides']} protospacers -> pseudo-chromosomes; files: guides.bed, isoforms.txt, guides.gtf, pseudo_genome.fa, protospacers.fa.",
                  f"STAR genome: {'built' if ran else 'not built'}.\n\n```\n{star_cmd}\n```"],
                 {"n_guides": ref["n_guides"], "star_ran": ran, "star_command": star_cmd, **{k: v for k, v in ref.items() if k != "n_guides"}})
    return 0


def read_fastq_seqs(path: Path, max_reads: Optional[int] = None):
    with _open_text(path) as fh:
        n = 0
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().rstrip()
            fh.readline()
            fh.readline()
            yield s
            n += 1
            if max_reads and n >= max_reads:
                return


def _hamming1_neighbors(s: str):
    for i, c in enumerate(s):
        for b in "ACGT":
            if b != c:
                yield s[:i] + b + s[i + 1:]


def correct_barcodes(raw_counts: Counter, whitelist: set) -> dict:
    """1MM_multi_Nbase_pseudocounts-style: exact match, else the whitelisted 1-mismatch neighbour with most exact reads."""
    exact = {b: c for b, c in raw_counts.items() if b in whitelist}
    out = {}
    for b in raw_counts:
        if b in whitelist:
            out[b] = b
            continue
        cands = [(exact.get(nb, 0) + 1, nb) for nb in _hamming1_neighbors(b.replace("N", "A")) if nb in whitelist]
        if "N" in b:
            cands += [(exact.get(nb, 0) + 1, nb) for nb in (b.replace("N", x, 1) for x in "ACGT") if nb in whitelist]
        if cands:
            cands.sort(reverse=True)
            out[b] = cands[0][1]
    return out


def dedup_umis_1mm(umi_counts: Counter) -> int:
    """Number of UMIs after collapsing each UMI into a more abundant one at Hamming distance 1 (1MM_CR)."""
    kept: list = []
    for u, _ in sorted(umi_counts.items(), key=lambda kv: (-kv[1], kv[0])):
        if any(sum(a != b for a, b in zip(u, k)) <= 1 for k in kept):
            continue
        kept.append(u)
    return len(kept)


def count_guides_fastq(r1: Path, r2: Path, whitelist: Path, guides, cb_start: int = 1, cb_len: int = 16, umi_start: int = 17,
                       umi_len: int = 12, max_reads: Optional[int] = None):
    """Cells x guides UMI matrix from paired FASTQs (starSoloGuide.nf coordinates, 1-based)."""
    np, pd = _np(), _pd()
    import anndata as ad  # type: ignore
    import scipy.sparse as sp  # type: ignore
    wl = set(x.strip().split()[0] for x in _open_text(whitelist) if x.strip())
    g = normalise_guide_metadata(guides)
    seq2id = {}
    for s, i in zip(g["protospacer"], g["sgRNA_ID"]):
        seq2id[s] = i
        seq2id.setdefault(revcomp(s), i)
    lens = sorted({len(s) for s in seq2id})
    stats_ = Counter()
    raw = []
    for s1, s2 in zip(read_fastq_seqs(r1, max_reads), read_fastq_seqs(r2, max_reads)):
        stats_["reads"] += 1
        cb = s1[cb_start - 1: cb_start - 1 + cb_len]
        umi = s1[umi_start - 1: umi_start - 1 + umi_len]
        hit = None
        for L in lens:
            for k in range(0, len(s2) - L + 1):
                gid = seq2id.get(s2[k:k + L])
                if gid is not None:
                    hit = gid
                    break
            if hit:
                break
        if hit is None:
            continue
        stats_["guide_matched"] += 1
        raw.append((cb, umi, hit))
    cb_map = correct_barcodes(Counter(r[0] for r in raw), wl)
    per = defaultdict(Counter)
    for cb, umi, gid in raw:
        c = cb_map.get(cb)
        if c is None:
            continue
        stats_["valid_barcode"] += 1
        per[(c, gid)][umi] += 1
    cells = sorted({k[0] for k in per})
    gids = list(dict.fromkeys(g["sgRNA_ID"]))
    ci = {c: i for i, c in enumerate(cells)}
    gi = {x: i for i, x in enumerate(gids)}
    M = sp.lil_matrix((len(cells), len(gids)))
    for (c, gid), umis in per.items():
        M[ci[c], gi[gid]] = dedup_umis_1mm(umis)
    var = pd.DataFrame({"protospacer": g.drop_duplicates("sgRNA_ID")["protospacer"].values}, index=gids)
    a = ad.AnnData(X=sp.csr_matrix(M), obs=pd.DataFrame(index=cells), var=var)
    stats_["cells"] = len(cells)
    stats_["umis"] = int(a.X.sum())
    return a, dict(stats_)


def parse_solo_out(solo_dir: Path, guides, compat: bool = False):
    """soloOutGuideParsing.py: Solo.out/Gene/raw matrix -> AnnData with sgRNA_ID names."""
    pd = _pd()
    import anndata as ad  # type: ignore
    import scipy.io  # type: ignore
    import scipy.sparse as sp  # type: ignore

    def find(name):
        hits = [p for p in solo_dir.rglob(name) if "raw" in str(p)] or list(solo_dir.rglob(name))
        if not hits:
            raise SystemExit(f"{name} not found under {solo_dir}")
        return hits[0]
    mtx = sp.csr_matrix(scipy.io.mmread(str(find("matrix.mtx"))))  # features x barcodes
    feats = [ln.rstrip("\n").split("\t") for ln in open(find("features.tsv")) if ln.strip()]
    bars = [ln.strip().split("-")[0] for ln in open(find("barcodes.tsv")) if ln.strip()]
    g = normalise_guide_metadata(guides)
    seq_to_id = dict(zip(g["protospacer"], g["sgRNA_ID"]))
    names = [seq_to_id.get((f[1] if len(f) > 1 else f[0]).split("_")[-1], (f[1] if len(f) > 1 else f[0])) for f in feats]
    if compat:
        return ad.AnnData(X=mtx, obs=pd.DataFrame(index=names), var=pd.DataFrame(index=bars))
    return ad.AnnData(X=sp.csr_matrix(mtx.T), obs=pd.DataFrame(index=bars), var=pd.DataFrame(index=names))


def cmd_count_guides(args: argparse.Namespace) -> int:
    np, pd = _np(), _pd()
    import scipy.io  # type: ignore
    d = run_dir(args.label)
    guides = read_table(Path(args.guides))
    if args.solo_out:
        a = parse_solo_out(Path(args.solo_out), guides, args.upstream_compat)
        stats_ = {"cells": a.n_obs, "umis": int(a.X.sum())}
    else:
        if not (args.r1 and args.r2 and args.whitelist):
            raise SystemExit("--r1, --r2 and --whitelist are required (or --solo-out)")
        star = (f"STAR --soloType CB_UMI_Simple --readFilesIn {args.r2} {args.r1} --genomeDir new_ref --soloCBwhitelist {args.whitelist} "
                f"--soloCBstart {args.cb_start} --soloCBlen {args.cb_len} --soloUMIstart {args.umi_start} --soloUMIlen {args.umi_len} "
                "--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --outFilterScoreMin 15 --outSAMtype BAM SortedByCoordinate "
                "--outSAMattributes CR UR CY UY CB UB --soloBarcodeReadLength 0 --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 "
                "--seedSearchStartLmax 23 --seedMultimapNmax 100 --alignMatesGapMax 1 --alignIntronMax 0 --outSAMunmapped Within --clipAdapterType CellRanger4")
        print(f"Equivalent STARsolo command (not required): {star}")
        a, stats_ = count_guides_fastq(Path(args.r1), Path(args.r2), Path(args.whitelist), guides, args.cb_start, args.cb_len,
                                       args.umi_start, args.umi_len, args.max_reads)
    p = d / "guides_star_solo.h5ad"
    a.write_h5ad(str(p))
    print(f"Wrote: {p}")
    mdir = d / "Solo.out" / "Gene" / "raw"
    mdir.mkdir(parents=True, exist_ok=True)
    scipy.io.mmwrite(str(mdir / "matrix.mtx"), a.X.T if not args.upstream_compat or not args.solo_out else a.X)
    (mdir / "barcodes.tsv").write_text("\n".join(a.obs_names) + "\n")
    (mdir / "features.tsv").write_text("\n".join(f"{v}\t{v}\tGene Expression" for v in a.var_names) + "\n")
    tot = pd.DataFrame({"guides_id": a.var_names, "total_UMIs": np.asarray(a.X.sum(axis=0)).ravel()})
    tot["target"] = tot["guides_id"].str.split("_").str[0]
    write_tsv(tot, d / "guide_totals.tsv")
    write_report(d, "CRISPR Jamboree 2 — guide counting", [md_table(["metric", "value"], stats_.items()),
                                                           md_table(["guide", "UMIs", "target"], tot.values.tolist(), 40)],
                 {"stats": stats_, "h5ad": p})
    return 0


# --- seqspec (minimal YAML reader for seqspec files: tags stripped) ---

def _yaml_scalar(v: str):
    v = v.strip()
    if v.startswith("!"):
        v = v.split(None, 1)[1] if " " in v else ""
    if v in ("null", "~", ""):
        return None
    if (v.startswith("'") and v.endswith("'")) or (v.startswith('"') and v.endswith('"')):
        return v[1:-1]
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return v


def mini_yaml(text: str):
    """Indentation YAML subset used by seqspec (mappings, '- ' lists, !Tags)."""
    try:
        import yaml  # type: ignore

        class _L(yaml.SafeLoader):
            pass
        _L.add_multi_constructor("!", lambda loader, suffix, node: loader.construct_mapping(node, deep=True)
                                 if isinstance(node, yaml.MappingNode) else loader.construct_scalar(node))
        return yaml.load(text, Loader=_L)
    except ImportError:
        pass
    lines = []
    for raw in text.splitlines():
        if not raw.strip() or raw.strip().startswith("#"):
            continue
        s = raw.rstrip()
        if s.strip().startswith("!") and ":" not in s:
            continue
        lines.append((len(s) - len(s.lstrip(" ")), s.strip()))

    def parse(i: int, indent: int):
        if i >= len(lines):
            return None, i
        if lines[i][1].startswith("- "):
            out = []
            while i < len(lines) and lines[i][0] == indent and lines[i][1].startswith("- "):
                item = lines[i][1][2:].strip()
                if item.startswith("!") and ":" not in item:
                    val, i = parse(i + 1, lines[i + 1][0] if i + 1 < len(lines) else indent + 2)
                    out.append(val)
                elif ":" in item:
                    lines[i] = (indent + 2, item)
                    val, i = parse(i, indent + 2)
                    out.append(val)
                else:
                    out.append(_yaml_scalar(item))
                    i += 1
            return out, i
        out = {}
        while i < len(lines) and lines[i][0] == indent and not lines[i][1].startswith("- "):
            k, _, v = lines[i][1].partition(":")
            v = v.strip()
            if v and not (v.startswith("!") and " " not in v):
                out[k.strip()] = _yaml_scalar(v)
                i += 1
            else:
                if i + 1 < len(lines) and (lines[i + 1][0] > indent or (lines[i + 1][0] == indent and lines[i + 1][1].startswith("- "))):
                    val, i = parse(i + 1, lines[i + 1][0])
                    out[k.strip()] = val
                else:
                    out[k.strip()] = None
                    i += 1
        return out, i
    return parse(0, lines[0][0] if lines else 0)[0]


def _leaves(region) -> list:
    subs = region.get("regions") or []
    if not subs:
        return [region]
    out = []
    for r in subs:
        out.extend(_leaves(r))
    return out


def seqspec_index(spec: dict, modality: str, reads: "list[str]") -> dict:
    """Positions of barcode / umi / cdna (region_type) on each read, as `seqspec index`."""
    lib = spec.get("library_spec") or []
    top = next((r for r in lib if r.get("region_id") == modality), None) or (lib[0] if lib else {})
    leaves = _leaves(top)
    rspecs = {r["read_id"]: r for r in (spec.get("sequence_spec") or [])}
    feats = {"barcode": [], "umi": [], "cdna": []}
    per_read = []
    for ri, rid in enumerate(reads):
        rs = rspecs.get(rid) or next((r for r in rspecs.values() if r.get("name") == rid), None)
        if rs is None:
            raise SystemExit(f"read {rid} not in sequence_spec ({list(rspecs)})")
        ids = [l.get("region_id") for l in leaves]
        pi = ids.index(rs["primer_id"]) if rs.get("primer_id") in ids else 0
        order = leaves[pi + 1:] if rs.get("strand", "pos") == "pos" else list(reversed(leaves[:pi]))
        pos, rlen = 0, int(rs.get("max_len") or 0)
        for reg in order:
            if pos >= rlen:
                break
            L = int(reg.get("max_len") or 0)
            end = min(pos + L, rlen)
            rt = reg.get("region_type")
            if rt in feats and L > 0:
                feats[rt].append((ri, pos, end))
                per_read.append({"read": rid, "region_type": rt, "region_id": reg.get("region_id"), "start": pos, "end": end})
            pos += L
    kb = ":".join(",".join(f"{r},{s},{e}" for r, s, e in feats[k]) for k in ("barcode", "umi", "cdna"))
    star = ""
    if feats["barcode"] and feats["umi"]:
        b, u = feats["barcode"][0], feats["umi"][0]
        star = f"--soloType CB_UMI_Simple --soloCBstart {b[1] + 1} --soloCBlen {b[2] - b[1]} --soloUMIstart {u[1] + 1} --soloUMIlen {u[2] - u[1]}"
    onlists = [l["onlist"].get("filename") for l in leaves if isinstance(l.get("onlist"), dict) and l.get("region_type") == "barcode"]
    return {"kb": kb, "starsolo": star, "regions": per_read, "whitelist": onlists[0] if onlists else None}


def cmd_seqspec_index(args: argparse.Namespace) -> int:
    spec = mini_yaml(Path(args.spec).read_text())
    reads = args.reads.split(",")
    idx = seqspec_index(spec, args.modality, reads)
    d = run_dir(args.label)
    print(f"seqspec index -t {args.tool}: {idx['kb'] if args.tool == 'kb' else idx['starsolo']}")
    pd = _pd()
    write_tsv(pd.DataFrame(idx["regions"]), d / "seqspec_regions.tsv")
    write_report(d, "CRISPR Jamboree 2 — seqspec index", [f"kb `-x`: `{idx['kb']}`", f"STARsolo: `{idx['starsolo']}`",
                                                         f"barcode whitelist: {idx['whitelist']}"], idx)
    return 0


# ---------------------------------------------------------------------------
# Simulation (simulation_engine/1-3 R scripts; DESeq2 dispersion fit re-implemented)
# ---------------------------------------------------------------------------

def size_factors_poscounts(counts):
    """DESeq2 estimateSizeFactors(type='poscounts') for a cells x genes matrix."""
    np = _np()
    with np.errstate(divide="ignore"):
        lc = np.where(counts > 0, np.log(np.maximum(counts, 1e-300)), 0.0)
    n = counts.shape[0]
    geo = lc.sum(axis=0) / n  # log geometric mean of positive counts, zeros contribute 0 (poscounts)
    sf = np.empty(n)
    for i in range(n):
        pos = counts[i] > 0
        sf[i] = math.exp(np.median((lc[i] - geo)[pos])) if pos.any() else np.nan
    sf = sf / math.exp(np.mean(np.log(sf[np.isfinite(sf) & (sf > 0)])))
    return sf


def _nb_loglik_alpha(y, mu, log_alpha, cr: bool = True):
    np = _np()
    from scipy.special import gammaln  # type: ignore
    a = math.exp(log_alpha)
    r = 1 / a
    ll = np.sum(gammaln(y + r) - gammaln(r) - gammaln(y + 1) + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu)))
    if cr:
        w = mu / (1 + a * mu)
        ll -= 0.5 * math.log(max(w.sum(), 1e-300))
    return ll


def deseq2_dispersions(counts, sf, fit_type: str = "parametric", disp_type: str = "dispersion") -> dict:
    """DESeq2 estimateDispersions for design ~ 1: gene-wise Cox-Reid MLE, parametric trend a0 + a1/mean, MAP shrinkage."""
    np = _np()
    from scipy.optimize import minimize_scalar  # type: ignore
    from scipy.special import polygamma  # type: ignore
    n, G = counts.shape
    norm = counts / sf[:, None]
    base_mean = norm.mean(axis=0)
    base_var = norm.var(axis=0, ddof=1)
    min_disp = 1e-8
    gw = np.full(G, np.nan)
    for g in range(G):
        y = counts[:, g]
        if y.sum() == 0:
            continue
        mu = sf * base_mean[g]
        mom = max((base_var[g] - base_mean[g] * np.mean(1 / sf)) / max(base_mean[g] ** 2, 1e-12), min_disp)
        res = minimize_scalar(lambda la: -_nb_loglik_alpha(y, mu, la), bounds=(math.log(min_disp), math.log(max(10.0, n))), method="bounded",
                              options={"xatol": 1e-6})
        gw[g] = math.exp(res.x) if res.success else mom
    ok = np.isfinite(gw) & (base_mean > 0)
    # parametric fit: gamma-family GLM (identity link) of gene-wise dispersions on 1/mean, iterated with outlier removal
    a0, a1 = 0.1, 1.0
    if fit_type == "mean":
        a0, a1 = float(np.mean(gw[ok & (gw > 100 * min_disp)])), 0.0
    else:
        use = ok & (gw > 100 * min_disp)
        for _ in range(10):
            fit = a0 + a1 / base_mean
            resid = gw / fit
            sel = use & (resid > 1e-4) & (resid < 15)
            if sel.sum() < 3:
                break
            Xg = np.column_stack([np.ones(sel.sum()), 1 / base_mean[sel]])
            yg = gw[sel]
            b = np.array([a0, a1])
            for _ in range(50):
                m = np.clip(Xg @ b, 1e-8, None)
                w = 1 / m ** 2
                b_new = np.linalg.solve((Xg * w[:, None]).T @ Xg, (Xg * w[:, None]).T @ yg)
                if np.max(np.abs(b_new - b)) < 1e-10:
                    b = b_new
                    break
                b = b_new
            conv = abs(math.log(max(b[0], 1e-12) / max(a0, 1e-12))) < 1e-6 and abs(math.log(max(b[1], 1e-12) / max(a1, 1e-12))) < 1e-6
            a0, a1 = float(max(b[0], 1e-8)), float(max(b[1], 1e-8))
            if conv:
                break
    trend = a0 + a1 / np.maximum(base_mean, 1e-12)
    m, p = n, 1
    resid = np.log(gw[ok]) - np.log(trend[ok])
    mad = 1.4826 * np.median(np.abs(resid - np.median(resid))) if resid.size else 0.0
    prior_var = max(mad ** 2 - float(polygamma(1, (m - p) / 2)), 0.25)
    mapd = np.full(G, np.nan)
    for g in range(G):
        if not ok[g]:
            continue
        y = counts[:, g]
        mu = sf * base_mean[g]
        lt = math.log(trend[g])
        res = minimize_scalar(lambda la: -(_nb_loglik_alpha(y, mu, la) - (la - lt) ** 2 / (2 * prior_var)),
                              bounds=(math.log(min_disp), math.log(max(10.0, n))), method="bounded", options={"xatol": 1e-6})
        mapd[g] = math.exp(res.x)
    outlier = ok & (np.log(np.where(ok, gw, 1)) > np.log(trend) + 2 * math.sqrt(prior_var))
    final = np.where(outlier, gw, mapd)
    chosen = {"dispersion": final, "dispFit": trend, "dispGeneEst": gw, "dispMAP": mapd}[disp_type]
    return {"mean": base_mean, "dispersion": chosen, "disp_outlier_deseq2": outlier, "dispGeneEst": gw, "dispFit": trend,
            "dispMAP": mapd, "size_factors": sf, "trend": (a0, a1), "prior_var": prior_var}


def simulate_counts(md: MuDataLite, pairs_to_perturb, guide_var: float = 0.0, seed: int = 1, compat: bool = False) -> "tuple[MuDataLite, dict]":
    """Fit NB per gene (2_fit_negative_binomial.R) then simulate counts with effects (3_simulate_count_data.R)."""
    np, pd = _np(), _pd()
    import scipy.sparse as sp  # type: ignore
    rng = np.random.default_rng(seed)
    gene = md.gene.copy()
    counts = dense(gene.X).astype(float)
    sf = size_factors_poscounts(counts)
    disp = deseq2_dispersions(counts, sf)
    gene.var["mean"] = disp["mean"]
    gene.var["dispersion"] = disp["dispersion"]
    gene.var["disp_outlier_deseq2"] = disp["disp_outlier_deseq2"]
    gene.obs["size_factors"] = sf
    eff = np.ones(counts.shape)
    A = assignment_matrix(md.guide)
    gv = md.guide.var
    gene_ix = {g: i for i, g in enumerate(gene.var_names)}
    for _, r in pairs_to_perturb.iterrows():
        es = float(r["effect_size"])
        guides = np.where(gv["intended_target_name"].astype(str).values == str(r["intended_target_name"]))[0]
        if r["gene_id"] not in gene_ix or guides.size == 0:
            continue
        gvar = guide_var if es != 1 else 0.0
        geff = {int(k): float(rng.normal(es, gvar)) if gvar > 0 else es for k in guides}
        gi = gene_ix[r["gene_id"]]
        for c in np.where(A[:, guides].any(axis=1))[0]:
            have = [k for k in guides if A[c, k]]
            eff[c, gi] = geff[int(rng.choice(have))]  # several guides for the pair: pick one at random
    mu = disp["mean"][None, :] * sf[:, None] * eff
    size = 1 / np.maximum(disp["dispersion"], 1e-8)
    p = size[None, :] / (size[None, :] + mu)
    sim = rng.negative_binomial(np.broadcast_to(size, mu.shape), np.clip(p, 1e-12, 1.0))
    if not compat:
        gene.X = sp.csr_matrix(sim.astype(np.float64))
        gene.obs["total_gene_umis"] = sim.sum(axis=1).astype(float)
        gene.obs["num_expressed_genes"] = (sim > 0).sum(axis=1)
    mods = dict(md.mods)
    mods["gene"] = gene
    out = MuDataLite(mods, md.obs, dict(md.uns))
    out.uns["pairs_to_test"] = pairs_to_perturb.reset_index(drop=True)
    return out, disp


def cmd_simulate(args: argparse.Namespace) -> int:
    np, pd = _np(), _pd()
    md = read_h5mu(Path(args.mudata))
    d = run_dir(args.label)
    if args.pairs:
        ptp = read_table(Path(args.pairs))
    elif args.perturb:
        ptp = pd.DataFrame([dict(zip(["gene_id", "intended_target_name", "effect_size"], s.split(":"))) for s in args.perturb])
    else:
        raise SystemExit("--pairs TSV (gene_id, intended_target_name, effect_size) or --perturb gene:target:effect is required")
    if "effect_size" not in ptp.columns:
        ptp["effect_size"] = args.effect_size
    ptp["effect_size"] = ptp["effect_size"].astype(float)
    if args.n_null_pairs:
        rng = np.random.default_rng(args.seed + 1)
        targets = sorted(set(md.guide.var["intended_target_name"].astype(str)) - {"non-targeting"})
        taken = set(zip(ptp["gene_id"], ptp["intended_target_name"]))
        extra = []
        for _ in range(args.n_null_pairs * 20):
            if len(extra) >= args.n_null_pairs:
                break
            g, t = str(rng.choice(md.gene.var_names)), str(rng.choice(targets))
            if (g, t) not in taken:
                taken.add((g, t))
                extra.append({"gene_id": g, "intended_target_name": t, "effect_size": 1.0})
        ptp = pd.concat([ptp, pd.DataFrame(extra)], ignore_index=True)
    sim, disp = simulate_counts(md, ptp, args.guide_var, args.seed, args.upstream_compat)
    out = Path(args.out) if args.out else d / "simulation_output.h5mu"
    write_h5mu(sim, out)
    tab = pd.DataFrame({"gene_id": md.gene.var_names, "mean": disp["mean"], "dispersion": disp["dispersion"],
                        "dispGeneEst": disp["dispGeneEst"], "dispFit": disp["dispFit"], "disp_outlier_deseq2": disp["disp_outlier_deseq2"]})
    write_tsv(tab, d / "gene_dispersions.tsv")
    write_tsv(ptp, d / "pairs_to_test.tsv")
    write_report(d, "CRISPR Jamboree 2 — Perturb-seq simulation",
                 [f"Size factors poscounts; parametric trend dispersion = {disp['trend'][0]:.4g} + {disp['trend'][1]:.4g}/mean; "
                  f"prior variance {disp['prior_var']:.3f}; guide_var {args.guide_var}.",
                  md_table(["gene_id", "target", "effect_size"], ptp[["gene_id", "intended_target_name", "effect_size"]].values.tolist(), 30)],
                 {"out": out, "n_pairs": len(ptp), "n_perturbed": int((ptp["effect_size"] != 1).sum()), "trend": disp["trend"],
                  "guide_var": args.guide_var, "upstream_compat": args.upstream_compat})
    return 0


# ---------------------------------------------------------------------------
# Evaluation (eval_and_clustergram, volcano_and_igv, network notebooks)
# ---------------------------------------------------------------------------

POS_TYPES = {"positive_control", "poscontrol", "positive", "direct targeting", "cis"}
NEG_TYPES = {"negative_control", "negcontrol", "negative", "non-targeting", "targeting_negative_control", "safe-targeting"}


def truth_labels(tr):
    """effect_size != 1 when present (simulation), else pair_type positive / negative controls (discovery pairs dropped)."""
    np = _np()
    if "effect_size" in tr.columns:
        return (tr["effect_size"].astype(float) != 1).astype(int).values, np.ones(len(tr), bool), "effect_size != 1"
    if "pair_type" in tr.columns:
        pt = tr["pair_type"].astype(str).str.lower().str.strip()
        lab = pt.isin(POS_TYPES).astype(int).values
        keep = (pt.isin(POS_TYPES) | pt.isin(NEG_TYPES)).values
        return lab, keep, "pair_type positive vs negative controls"
    return None, None, "no truth column"


def evaluate_results(tr, compat: bool = False) -> dict:
    np = _np()
    from scipy import stats  # type: ignore
    lab, keep, how = truth_labels(tr)
    out = {"truth": how}
    if lab is None:
        return out
    sc = sig_column(tr)
    p = tr[sc].astype(float).values
    score = p if compat else -np.log10(np.clip(p, 1e-300, 1.0))
    out.update({"score": f"raw {sc}" if compat else f"-log10({sc})"})
    out.update(binary_metrics(lab[keep], score[keep]))
    lc = lfc_column(tr)
    if "effect_size" in tr.columns and lc:
        te = np.log2(tr["effect_size"].astype(float).values)
        pe = tr[lc].astype(float).values
        ok = np.isfinite(te) & np.isfinite(pe)
        if ok.sum() > 2 and np.std(te[ok]) > 0 and np.std(pe[ok]) > 0:
            out["pearson_r_log2fc"] = float(stats.pearsonr(te[ok], pe[ok])[0])
        out["mse_log2fc"] = float(np.mean((te[ok] - pe[ok]) ** 2)) if ok.any() else float("nan")
    return out


def clustergram_matrix(md: MuDataLite):
    """Per gRNA log2(mean expression in cells with the gRNA / without), guide .X != 0, non-finite -> 0."""
    np, pd = _np(), _pd()
    G = dense(md.guide.X) != 0
    X = dense(md.gene.X).astype(float)
    tot = X.sum(axis=0)
    n = X.shape[0]
    rows = []
    for j in range(G.shape[1]):
        has = np.where(G[:, j])[0]
        wi = X[has].sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            lfc = np.log2((wi / max(len(has), 1)) / ((tot - wi) / max(n - len(has), 1)))
        lfc[~np.isfinite(lfc)] = 0
        if len(has) == 0:
            lfc[:] = 0
        rows.append(lfc)
    return pd.DataFrame(np.array(rows), index=md.guide.var_names, columns=md.gene.var_names)


def cluster_order(M) -> "tuple[list, list]":
    np = _np()
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage  # type: ignore
        r = leaves_list(linkage(M, method="average", metric="euclidean")) if M.shape[0] > 2 else np.arange(M.shape[0])
        c = leaves_list(linkage(M.T, method="average", metric="euclidean")) if M.shape[1] > 2 else np.arange(M.shape[1])
        return list(r), list(c)
    except Exception:
        return list(range(M.shape[0])), list(range(M.shape[1]))


def plot_pr_roc(evals: dict, path: Path):
    plt = _plt()
    if plt is None:
        return None
    ms = [m for m, e in evals.items() if "_curves" in e]
    if not ms:
        return None
    fig, axes = plt.subplots(len(ms), 2, figsize=(10, 4 * len(ms)), squeeze=False)
    for i, m in enumerate(ms):
        c = evals[m]["_curves"]
        axes[i, 0].plot(c["recall"], c["precision"], lw=1.2, color=BLUE, label=f"AUPRC={evals[m]['auprc']:.3f}")
        axes[i, 1].plot(c["fpr"], c["tpr"], lw=1.2, color=BLUE, label=f"AUROC={evals[m]['auroc']:.3f}")
        _style(axes[i, 0], f"{m} - Precision-Recall", "Recall", "Precision")
        _style(axes[i, 1], f"{m} - ROC", "False Positive Rate", "True Positive Rate")
        axes[i, 0].legend(frameon=False, fontsize=8)
        axes[i, 1].legend(frameon=False, fontsize=8)
    return _save(fig, path)


def cmd_evaluate(args: argparse.Namespace) -> int:
    np, pd = _np(), _pd()
    d = run_dir(args.label)
    results = {}
    md = read_h5mu(Path(args.mudata)) if args.mudata else None
    if args.results:
        for p in args.results:
            results[Path(p).name.split(".")[1] if Path(p).name.count(".") >= 2 else Path(p).stem] = read_table(Path(p))
    elif md is not None:
        results["test_results"] = uns_df(md, "test_results")
    evals = {m: evaluate_results(tr, args.upstream_compat) for m, tr in results.items()}
    rows = [[m, e.get("truth"), e.get("n", ""), e.get("n_pos", ""), _fmt(e.get("auprc")), _fmt(e.get("auroc")),
             _fmt(e.get("pearson_r_log2fc")), _fmt(e.get("mse_log2fc"))] for m, e in evals.items()]
    figs = []
    if not args.no_plots:
        f = plot_pr_roc(evals, d / "evaluation_curves.png")
        if f:
            figs.append(f)
    sections = [md_table(["results", "truth", "n", "positives", "AUPRC", "AUROC", "Pearson r (log2FC)", "MSE"], rows)]
    if md is not None and not args.no_clustergram:
        M = clustergram_matrix(md)
        r, c = cluster_order(M.values)
        Mo = M.iloc[r, c]
        write_tsv(Mo.reset_index().rename(columns={"index": "guide_id"}), d / "clustergram_lfc.tsv")
        sections.append(f"Clustergram: {M.shape[0]} gRNAs x {M.shape[1]} genes (log2 mean ratio, average linkage).")
        plt = _plt()
        if plt is not None and not args.no_plots:
            fig, ax = plt.subplots(figsize=(7, 6))
            v = float(np.nanmax(np.abs(Mo.values))) or 1.0
            im = ax.imshow(Mo.values, cmap="coolwarm", vmin=-v, vmax=v, aspect="auto")
            fig.colorbar(im, ax=ax, shrink=0.7, label="log2 FC")
            _style(ax, "Clustergram (gRNA x gene)", "genes", "gRNAs")
            figs.append(_save(fig, d / "clustergram.png"))
    summary = {m: {k: v for k, v in e.items() if k != "_curves"} for m, e in evals.items()}
    write_report(d, "CRISPR Jamboree 2 — evaluation", sections, {"evaluations": summary, "figures": figs})
    return 0


def volcano_tables(tr, lfc_thresh: float = 2.0, sig_thresh: float = 0.01):
    sc = sig_column(tr)
    lc = lfc_column(tr)
    if lc is None:
        raise SystemExit("test_results has no log2_fc column (the R wilcoxon module reports p-values only); use another method's table")
    down = tr[(tr[lc] <= -lfc_thresh) & (tr[sc] <= sig_thresh)]
    up = tr[(tr[lc] >= lfc_thresh) & (tr[sc] <= sig_thresh)]
    return down, up, sc, lc


def igv_tracks(md: MuDataLite, tr):
    """volcano_and_igv.ipynb: bedgraph for promoter pairs, bedpe for element -> gene links."""
    np, pd = _np(), _pd()
    coord = {}
    gv = md.gene.var
    if {"gene_chr", "gene_start", "gene_end"} <= set(gv.columns):
        for idx, r in gv.iterrows():
            try:
                if np.isnan(float(r["gene_start"])) or np.isnan(float(r["gene_end"])):
                    continue
            except (TypeError, ValueError):
                continue
            coord[idx] = [r["gene_chr"], r["gene_start"], r["gene_end"]]
    for _, r in md.guide.var.iterrows():
        t = r["intended_target_name"]
        if t in coord or t == "non-targeting" or "intended_target_chr" not in md.guide.var.columns:
            continue
        coord[t] = [r["intended_target_chr"], r["intended_target_start"], r["intended_target_end"]]
    sc = sig_column(tr)
    lc = lfc_column(tr) or "log2_fc"
    bp, bg = defaultdict(list), defaultdict(list)
    for _, r in tr.iterrows():
        t, g = r["intended_target_name"], r["gene_id"]
        if t not in coord:
            continue
        if t == g:
            for k, v in zip(("chr", "start", "end"), coord[t]):
                bg[k].append(v)
            bg[sc].append(r[sc]); bg["log2_fc"].append(r.get(lc))
        elif g in coord:
            for k, v in zip(("chr1", "start1", "end1"), coord[t]):
                bp[k].append(v)
            for k, v in zip(("chr2", "start2", "end2"), coord[g]):
                bp[k].append(v)
            bp[sc].append(r[sc]); bp["log2_fc"].append(r.get(lc))
    return pd.DataFrame(bp), pd.DataFrame(bg)


def plot_volcano(tr, down, up, sc, lc, lfc_thresh, sig_thresh, path: Path, label_col: str = "gene_id"):
    np = _np()
    plt = _plt()
    if plt is None:
        return None
    fig, ax = plt.subplots(figsize=(6, 4.5))
    y = -np.log10(np.clip(tr[sc].astype(float), 1e-300, 1))
    ax.scatter(tr[lc], y, s=5, color=GREEN, label="Not significant")
    ax.scatter(down[lc], -np.log10(np.clip(down[sc].astype(float), 1e-300, 1)), s=10, color=BLUE, label="Down-regulated")
    ax.scatter(up[lc], -np.log10(np.clip(up[sc].astype(float), 1e-300, 1)), s=10, color=RED, label="Up-regulated")
    for _, r in down.head(15).iterrows():
        ax.text(r[lc], -math.log10(max(float(r[sc]), 1e-300)), str(r[label_col]), fontsize=6, color=INK2)
    for v in (-lfc_thresh, lfc_thresh):
        ax.axvline(v, color=GREY, ls="--", lw=0.8)
    ax.axhline(-math.log10(sig_thresh), color=GREY, ls="--", lw=0.8)
    finite = tr[lc][np.isfinite(tr[lc].astype(float))]
    b = (float(np.max(np.abs(finite))) if len(finite) else 1) + 0.5
    ax.set_xlim(-b, b)
    ax.legend(loc="upper right", frameon=False, fontsize=7)
    _style(ax, "Volcano", "log2 fold change", f"-log10 {sc}")
    return _save(fig, path)


def cmd_volcano(args: argparse.Namespace) -> int:
    md = read_h5mu(Path(args.mudata))
    tr = read_table(Path(args.results)) if args.results else uns_df(md, "test_results")
    d = run_dir(args.label)
    down, up, sc, lc = volcano_tables(tr, args.log2fc_threshold, args.sig_threshold)
    write_tsv(down, d / "volcano_down.tsv")
    write_tsv(up, d / "volcano_up.tsv")
    bp, bg = igv_tracks(md, tr)
    write_tsv(bp, d / "igv_links_bedpe.tsv")
    write_tsv(bg, d / "igv_promoter_bedgraph.tsv")
    figs = []
    if not args.no_plots:
        f = plot_volcano(tr, down, up, sc, lc, args.log2fc_threshold, args.sig_threshold, d / "volcano.png")
        if f:
            figs.append(f)
    write_report(d, "CRISPR Jamboree 2 — volcano and IGV tracks",
                 [f"{len(down)} down-regulated and {len(up)} up-regulated pairs (|{lc}| >= {args.log2fc_threshold}, {sc} <= {args.sig_threshold}).",
                  f"IGV: {len(bg)} promoter bedgraph rows, {len(bp)} bedpe links."],
                 {"n_down": len(down), "n_up": len(up), "n_bedgraph": len(bg), "n_bedpe": len(bp), "figures": figs})
    return 0


def mean_nontargeting_expression(md: MuDataLite):
    """Mean expression over cells with >= 1 non-targeting gRNA count and no targeting gRNA count (guide .X)."""
    np = _np()
    G = dense(md.guide.X) > 0
    nt = ~is_targeting(md.guide.var)
    cells = G[:, nt].any(axis=1) & ~G[:, ~nt].any(axis=1)
    X = md.gene.X[np.where(cells)[0]]
    return np.asarray(X.mean(axis=0)).ravel() if cells.any() else np.full(md.gene.n_vars, np.nan)


def network_edges(tr, central_node: Optional[str], source: str = "intended_target_name", target: str = "gene_id",
                  weight: str = "log2_fc", min_weight: Optional[float] = 0.01):
    df = tr.drop_duplicates()
    if min_weight is not None:
        df = df[df[weight].abs() >= min_weight]
    if central_node is not None:
        df = df[df[source] == central_node]
    return df


def plot_network(edges, central: str, source: str, target: str, weight: str, size_col: Optional[str], ax):
    np = _np()
    nodes = list(dict.fromkeys(list(edges[source]) + list(edges[target])))
    k = len(nodes)
    pos = {n: (math.cos(2 * math.pi * i / max(k, 1)), math.sin(2 * math.pi * i / max(k, 1))) for i, n in enumerate(nodes)}
    import matplotlib.cm as cm  # type: ignore
    wmax = float(edges[weight].abs().max()) or 1.0
    for _, r in edges.iterrows():
        (x0, y0), (x1, y1) = pos[r[source]], pos[r[target]]
        col = cm.coolwarm(0.5 + 0.5 * float(r[weight]) / wmax)
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="->", color=col, lw=1.2))
        ax.text((x0 + x1) / 2, (y0 + y1) / 2, f"{r[weight]:.2f}", fontsize=6, color=INK2)
    sizes = {}
    if size_col and size_col in edges.columns:
        sizes = edges.set_index(target)[size_col].to_dict()
    for n, (x, y) in pos.items():
        ax.scatter([x], [y], s=max(float(sizes.get(n, 1.5)), 0.1) * 100 if sizes else 150, color="#8ec5e8", zorder=3)
        ax.text(x + 0.12, y + 0.05, str(n), fontsize=7, color=INK)
    ax.set_axis_off()
    ax.set_title(f"{central}", loc="left", fontsize=9, fontweight="bold")


def cmd_network(args: argparse.Namespace) -> int:
    np, pd = _np(), _pd()
    md = read_h5mu(Path(args.mudata))
    tr = read_table(Path(args.results)) if args.results else uns_df(md, "test_results")
    d = run_dir(args.label)
    sc = sig_column(tr)
    lc = lfc_column(tr) or "log2_fc"
    tr = tr.copy()
    tr[f"-log10_{sc}"] = -np.log10(np.clip(tr[sc].astype(float), 1e-300, 1))
    centrals = args.central_nodes or list(tr.sort_values(sc)["intended_target_name"].unique()[: args.n_central])
    all_edges = []
    for c in centrals:
        e = network_edges(tr, c, weight=lc, min_weight=args.min_weight)
        e = e.assign(central_node=c)
        all_edges.append(e)
    edges = pd.concat(all_edges, ignore_index=True) if all_edges else pd.DataFrame()
    write_tsv(edges, d / "network_edges.tsv")
    nt = mean_nontargeting_expression(md)
    gv = pd.DataFrame({"gene_id": md.gene.var_names, "mean_nontargeting_expression": nt}).sort_values("mean_nontargeting_expression", ascending=False)
    write_tsv(gv, d / "mean_nontargeting_expression.tsv")
    figs = []
    plt = _plt()
    if plt is not None and not args.no_plots and len(edges):
        cols = 2
        rows_ = (len(centrals) + 1) // 2
        fig, axes = plt.subplots(rows_, cols, figsize=(10, 5 * rows_), squeeze=False)
        for i, c in enumerate(centrals):
            e = edges[edges["central_node"] == c]
            if len(e):
                plot_network(e, c, "intended_target_name", "gene_id", lc, f"-log10_{sc}", axes[i // 2, i % 2])
            else:
                axes[i // 2, i % 2].set_axis_off()
        for j in range(len(centrals), rows_ * cols):
            axes[j // 2, j % 2].set_axis_off()
        figs.append(_save(fig, d / "network_plot.png"))
    write_report(d, "CRISPR Jamboree 2 — network",
                 [f"Central nodes: {', '.join(map(str, centrals))}; {len(edges)} edges with |{lc}| >= {args.min_weight}.",
                  md_table(["gene", "mean NT expression"], gv.head(10).values.tolist())],
                 {"central_nodes": centrals, "n_edges": len(edges), "figures": figs})
    return 0


# ---------------------------------------------------------------------------
# inspect / run
# ---------------------------------------------------------------------------

def schema_check(md: MuDataLite) -> "list[tuple[str, bool]]":
    checks = [("modality gene", "gene" in md), ("modality guide", "guide" in md)]
    if "guide" in md:
        g = md.guide
        checks += [("guide.layers['guide_assignment']", "guide_assignment" in g.layers),
                   ("guide.var.targeting", "targeting" in g.var.columns),
                   ("guide.var.intended_target_name", "intended_target_name" in g.var.columns),
                   ("guide.uns.moi", "moi" in g.uns), ("guide.uns.capture_method", "capture_method" in g.uns)]
    checks.append(("uns.pairs_to_test", "pairs_to_test" in md.uns))
    return checks


def cmd_inspect(args: argparse.Namespace) -> int:
    md = read_h5mu(Path(args.mudata))
    d = run_dir(args.label)
    secs = [md_table(["modality", "cells", "features", "obs columns", "var columns", "layers"],
                     [[k, a.n_obs, a.n_vars, ", ".join(map(str, a.obs.columns)), ", ".join(map(str, a.var.columns)), ", ".join(a.layers.keys())]
                      for k, a in md.mods.items()])]
    checks = schema_check(md)
    secs.append(md_table(["required field", "present"], [[c, "yes" if ok else "NO"] for c, ok in checks]))
    info = {"modalities": {k: [a.n_obs, a.n_vars] for k, a in md.mods.items()}, "schema": dict(checks),
            "global_obs": list(md.obs.columns), "uns": sorted(md.uns)}
    if "guide" in md:
        info["moi"] = moi_of(md)
        info["capture_method"] = _uns_scalar(md.guide, "capture_method", "")
    for key in ("pairs_to_test", "test_results"):
        t = uns_df(md, key, required=False)
        if t is not None:
            info[key] = {"n": len(t), "columns": list(t.columns)}
            if "pair_type" in t.columns:
                info[key]["pair_type"] = t["pair_type"].astype(str).value_counts().to_dict()
            secs.append(f"**{key}** ({len(t)} rows)\n\n" + md_table(list(map(str, t.columns)), t.head(8).values.tolist()))
    write_report(d, f"CRISPR Jamboree 2 — MuData inspection: {Path(args.mudata).name}", secs, info)
    return 0 if all(ok for _, ok in checks) or not args.strict else 1


def cmd_run(args: argparse.Namespace) -> int:
    np, pd = _np(), _pd()
    md = read_h5mu(Path(args.mudata))
    d = run_dir(args.label)
    sections, summary = [], {"mudata": args.mudata, "run_dir": d}
    if "guide_assignment" not in md.guide.layers or args.assign:
        A, _ = assign_guides(md, args.assign or "umi-threshold", args.threshold)
        md.guide.layers["guide_assignment"] = A
        sections.append(f"Guide assignment: `{args.assign or 'umi-threshold'}` (threshold {args.threshold}).")
    evals, results = {}, {}
    rows = []
    for m in args.methods:
        tr = run_inference(md, m, side=args.side, B=args.resamples, seed=args.seed, compat=args.upstream_compat, work=d)
        if tr is None:
            rows.append([m, "skipped", "", "", ""])
            continue
        results[m] = tr
        write_tsv(tr, d / f"test_results.{m}.tsv")
        e = evaluate_results(tr, args.upstream_compat)
        evals[m] = e
        rows.append([m, len(tr), int((tr[sig_column(tr)] < 0.05).sum()), _fmt(e.get("auprc")), _fmt(e.get("auroc"))])
    sections.append(md_table(["method", "pairs", "p < 0.05", "AUPRC", "AUROC"], rows))
    figs = []
    if not args.no_plots:
        f = plot_pr_roc(evals, d / "evaluation_curves.png")
        if f:
            figs.append(f)
    with_lfc = [m for m, t in results.items() if lfc_column(t)]
    if with_lfc:
        first = with_lfc[0]
        tr = results[first]
        md.uns["test_results"] = tr
        down, up, sc, lc = volcano_tables(tr, args.log2fc_threshold, args.sig_threshold)
        bp, bg = igv_tracks(md, tr)
        write_tsv(bp, d / "igv_links_bedpe.tsv")
        write_tsv(bg, d / "igv_promoter_bedgraph.tsv")
        sections.append(f"Volcano ({first}): {len(down)} down, {len(up)} up; IGV {len(bg)} promoter rows, {len(bp)} links.")
        if not args.no_plots:
            f = plot_volcano(tr, down, up, sc, lc, args.log2fc_threshold, args.sig_threshold, d / "volcano.png")
            if f:
                figs.append(f)
        if lc:
            centrals = list(tr.sort_values(sc)["intended_target_name"].unique()[:2])
            edges = pd.concat([network_edges(tr, c, weight=lc).assign(central_node=c) for c in centrals], ignore_index=True)
            write_tsv(edges, d / "network_edges.tsv")
        if args.out:
            write_h5mu(md, Path(args.out))
    summary.update({"methods": {m: {k: v for k, v in e.items() if k != "_curves"} for m, e in evals.items()}, "figures": figs})
    write_report(d, "CRISPR Jamboree 2 — Perturb-seq inference benchmark", sections, summary)
    return 0


# ---------------------------------------------------------------------------
# Synthetic data + self-test
# ---------------------------------------------------------------------------

def synthetic_mudata(dirpath: Path, seed: int = 3, moi: str = "high", n_cells: int = 1200, capture: str = "CROP-seq") -> dict:
    """Cells x 40 genes; 5 elements x 2 gRNAs + 4 NT gRNAs.

    Planted: E1 -> G1 (x0.3), E2 -> G2 (x0.4), promoter of G3 (x0.25, target named after gene G3's id), E4 -> G4 (x2.5 up).
    gRNA UMIs: assigned -> NB(mean 40, size 3); background Poisson(0.4).
    """
    np, pd = _np(), _pd()
    import anndata as ad  # type: ignore
    import scipy.sparse as sp  # type: ignore
    rng = np.random.default_rng(seed)
    n_genes = 150
    genes = [f"ENSG{900000 + i:08d}" for i in range(n_genes)]
    elements = ["E1", "E2", genes[3], "E4", "E5"]
    guide_ids = [f"{e}_g{k}" for e in elements for k in (1, 2)] + [f"NT_{k}" for k in range(4)]
    target_of = {g: (g.rsplit("_g", 1)[0] if not g.startswith("NT_") else "non-targeting") for g in guide_ids}
    ng = len(guide_ids)
    A = np.zeros((n_cells, ng), dtype=bool)
    if moi == "high":
        for c in range(n_cells):
            k = rng.poisson(2.0) + 1
            A[c, rng.choice(ng, size=min(k, ng), replace=False)] = True
    else:
        A[np.arange(n_cells), rng.integers(0, ng, n_cells)] = True
    lib = rng.lognormal(0, 0.3, n_cells)
    base = rng.gamma(2.0, 2.0, n_genes) + 0.5
    mu = lib[:, None] * base[None, :]
    effects = {("E1", genes[1]): 0.3, ("E2", genes[2]): 0.4, (genes[3], genes[3]): 0.25, ("E4", genes[4]): 2.5}
    for (e, g), f in effects.items():
        gi = genes.index(g)
        cells = A[:, [i for i, x in enumerate(guide_ids) if target_of[x] == e]].any(axis=1)
        mu[cells, gi] *= f
    size = 5.0
    X = rng.negative_binomial(size, size / (size + mu)).astype(np.float64)
    guide_umi = np.where(A, rng.negative_binomial(3, 3 / (3 + 40.0), A.shape) + 1, rng.poisson(0.4, A.shape)).astype(np.float64)
    cells = [f"CELL{c:05d}-1" for c in range(n_cells)]
    batch = np.array(["b1" if c % 2 else "b2" for c in range(n_cells)])
    gene = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame({"num_expressed_genes": (X > 0).sum(axis=1), "total_gene_umis": X.sum(axis=1)}, index=cells),
                      var=pd.DataFrame({"symbol": [f"G{i}" for i in range(n_genes)], "gene_chr": "chr1",
                                        "gene_start": [1_000_000 + 50_000 * i for i in range(n_genes)],
                                        "gene_end": [1_000_000 + 50_000 * i + 1 for i in range(n_genes)]}, index=genes))
    coords = {e: (1_000_000 + 50_000 * i - 20_000) for i, e in enumerate(elements)}
    gvar = pd.DataFrame({"targeting": ["FALSE" if target_of[g] == "non-targeting" else "TRUE" for g in guide_ids],
                         "intended_target_name": [target_of[g] for g in guide_ids],
                         "intended_target_chr": ["chr1" if target_of[g] != "non-targeting" else "NA" for g in guide_ids],
                         "intended_target_start": [float(coords.get(target_of[g], np.nan)) for g in guide_ids],
                         "intended_target_end": [float(coords.get(target_of[g], np.nan) + 500) for g in guide_ids],
                         "protospacer": ["".join(rng.choice(list("ACGT"), 20)) for _ in guide_ids]}, index=guide_ids)
    guide = ad.AnnData(X=sp.csr_matrix(guide_umi), obs=pd.DataFrame({"num_expressed_guides": (guide_umi > 0).sum(axis=1),
                                                                    "total_guide_umis": guide_umi.sum(axis=1)}, index=cells), var=gvar)
    guide.layers["guide_assignment"] = sp.csr_matrix(A.astype(np.float64))
    guide.uns["moi"] = np.array([moi], dtype=object)
    guide.uns["capture_method"] = np.array([capture], dtype=object)
    rows = []
    for e in elements:
        for gi in range(0, 12):
            g = genes[gi]
            pt = "positive_control" if (e, g) in effects else ("negative_control" if gi >= 6 else "discovery")
            rows.append({"gene_id": g, "intended_target_name": e, "pair_type": pt})
    pairs = pd.DataFrame(rows)
    md = MuDataLite({"gene": gene, "guide": guide}, pd.DataFrame({"prep_batch": batch}, index=cells), {"pairs_to_test": pairs})
    path = dirpath / f"synthetic_{moi}.h5mu"
    write_h5mu(md, path)
    return {"h5mu": path, "genes": genes, "elements": elements, "effects": effects, "guide_ids": guide_ids, "A": A,
            "guide_umi": guide_umi, "cells": cells, "protospacers": dict(zip(guide_ids, gvar["protospacer"]))}


def synthetic_fastqs(dirpath: Path, W: dict, seed: int = 5, n_cells: int = 30, reads_per: int = 6) -> dict:
    """R1 = CB(16)+UMI(12), R2 = 5 nt + protospacer + scaffold; whitelist; planted cell -> guide counts."""
    np = _np()
    rng = np.random.default_rng(seed)
    bcs = ["".join(rng.choice(list("ACGT"), 16)) for _ in range(n_cells)]
    wl = dirpath / "whitelist.txt"
    wl.write_text("\n".join(bcs + ["".join(rng.choice(list("ACGT"), 16)) for _ in range(50)]) + "\n")
    gids = W["guide_ids"][:6]
    truth = {}
    r1, r2 = [], []
    for ci, bc in enumerate(bcs):
        g = gids[ci % len(gids)]
        n_umi = 3 + ci % 4
        truth[(bc, g)] = n_umi
        for u in range(n_umi):
            umi = "".join(rng.choice(list("ACGT"), 12))
            for k in range(reads_per):
                cb = bc
                if k == 0 and ci % 5 == 0:  # one read per UMI with a sequencing error in the CB -> corrected
                    cb = bc[:7] + ("A" if bc[7] != "A" else "C") + bc[8:]
                umi_k = umi if k != 1 else umi[:11] + ("A" if umi[11] != "A" else "G")  # 1-mismatch UMI -> collapsed
                r1.append(cb + umi_k)
                ps = W["protospacers"][g]
                r2.append("ACGTA" + (ps if k % 2 == 0 else revcomp(ps)) + "GTTTAAGAGCTATGCTGGAAAC")
    for _ in range(20):
        r1.append("".join(rng.choice(list("ACGT"), 28)))
        r2.append("".join(rng.choice(list("ACGT"), 47)))
    p1, p2 = dirpath / "R1.fastq.gz", dirpath / "R2.fastq.gz"
    for p, seqs in ((p1, r1), (p2, r2)):
        with gzip.open(p, "wt") as fh:
            for i, s in enumerate(seqs):
                fh.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")
    return {"r1": p1, "r2": p2, "whitelist": wl, "truth": truth}


SEQSPEC_EXAMPLE = """!Assay
seqspec_version: 0.2.0
assay_id: TEST_10XV3_Guide
modalities:
- crispr
sequence_spec:
- !Read
  read_id: R1
  name: Read 1
  modality: crispr
  primer_id: r1_primer
  min_len: 28
  max_len: 28
  strand: pos
- !Read
  read_id: R2
  name: Read 2
  modality: crispr
  primer_id: r2_primer
  min_len: 83
  max_len: 83
  strand: neg
library_spec:
- !Region
  parent_id: null
  region_id: crispr
  region_type: null
  min_len: 111
  max_len: 111
  regions:
  - !Region
    region_id: r1_primer
    region_type: r1_primer
    min_len: 0
    max_len: 0
    regions: null
  - !Region
    region_id: barcode
    region_type: barcode
    sequence: NNNNNNNNNNNNNNNN
    min_len: 16
    max_len: 16
    onlist: !Onlist
      location: remote
      filename: 737K-august-2016.txt
    regions: null
  - !Region
    region_id: umi
    region_type: umi
    min_len: 12
    max_len: 12
    regions: null
  - !Region
    region_id: cdna
    region_type: cdna
    min_len: 20
    max_len: 20
    regions: null
  - !Region
    region_id: common
    region_type: common
    sequence: CAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCTTAGCTCTTAAAC
    min_len: 63
    max_len: 63
    regions: null
  - !Region
    region_id: r2_primer
    region_type: r2_primer
    min_len: 0
    max_len: 0
    regions: null
"""


def cmd_selftest(args: argparse.Namespace) -> int:
    import tempfile
    np, pd = _np(), _pd()
    checks = []

    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)

    made: list = []

    def last(label):
        r = sorted(OUT_ROOT.glob(f"*_{label}"))[-1]
        made.append(r)
        return r

    try:
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            W = synthetic_mudata(d, moi="high")
            WL = synthetic_mudata(d, moi="low", seed=8)
            md = read_h5mu(W["h5mu"])
            g1, g2, g3, g4 = W["genes"][1:5]

            print("\nstatistics")
            from scipy import stats  # type: ignore
            rng = np.random.default_rng(0)
            x, y = rng.poisson(3, 40).astype(float), rng.poisson(5, 60).astype(float)
            p = r_wilcox(x, y, "both")
            ref = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic", use_continuity=True).pvalue
            check(abs(p - ref) < 1e-12 and p < 0.01 and r_wilcox(x, y, "left") < p, f"wilcox.test: ties -> normal approx + continuity (p {p:.2e}); left-sided smaller")
            check(abs(r_wilcox([1, 2, 3.5], [4, 5, 6, 7], "both") - 2 / 35) < 1e-12, "wilcox.test: exact distribution without ties (p = 2/35)")
            n = 3000
            xx = (rng.random(n) < 0.3).astype(float)
            off = np.log(rng.lognormal(7, 0.3, n))
            mu = np.exp(-6 + off + math.log(0.5) * xx)
            yy = rng.negative_binomial(4, 4 / (4 + mu)).astype(float)
            fit = glm_nb(yy, np.column_stack([np.ones(n), xx]), offset=off)
            check(abs(fit["beta"][1] - math.log(0.5)) < 0.12 and abs(fit["theta"] - 4) < 1.5 and fit["p"][1] < 1e-10,
                  f"glm.nb: recovers log fold change {fit['beta'][1]:.3f} (true {math.log(0.5):.3f}) and theta {fit['theta']:.2f} (true 4)")
            Xl = np.log1p(rng.poisson(4, (200, 3)).astype(float))
            grp = np.arange(200) < 60
            Xl[grp, 0] += 0.8
            s, pw, lf = scanpy_rank_test(Xl, grp, "wilcoxon")
            try:
                import anndata as ad  # type: ignore
                import scanpy as sc  # type: ignore
                a = ad.AnnData(X=Xl.copy())
                a.obs["g"] = pd.Categorical(np.where(grp, "present", "not_present"))
                a.uns["log1p"] = {"base": None}
                for meth in ("wilcoxon", "t-test", "t-test_overestim_var"):
                    sc.tl.rank_genes_groups(a, groupby="g", groups=["present"], method=meth)
                    rg = a.uns["rank_genes_groups"]
                    order = [int(v) for v in rg["names"]["present"]]
                    ref_p = np.empty(3); ref_l = np.empty(3)
                    ref_p[order] = rg["pvals"]["present"]; ref_l[order] = rg["logfoldchanges"]["present"]
                    _, mp, ml = scanpy_rank_test(Xl, grp, meth)
                    check(np.allclose(mp, ref_p, rtol=1e-5, atol=1e-12) and np.allclose(ml, ref_l, rtol=1e-4, atol=1e-5),
                          f"scanpy port matches sc.tl.rank_genes_groups({meth}) p-values and logfoldchanges")
            except ImportError:
                check(pw[0] < 1e-4 and pw[1] > 1e-3, "scanpy wilcoxon: planted shift detected (scanpy not installed for cross-check)")
            bl = binary_metrics([1, 1, 0, 0, 1], [0.9, 0.8, 0.3, 0.2, 0.1])
            check(abs(bl["auroc"] - 4 / 6) < 1e-12, f"binary_metrics: AUROC {bl['auroc']:.3f} = 4/6")

            print("\nguide assignment")
            mdl = read_h5mu(W["h5mu"])
            A_thr = dense(assign_threshold(mdl.guide, 5)) > 0
            truth = W["A"]
            acc_thr = (A_thr == truth).mean()
            check(acc_thr > 0.97, f"umi-threshold >= 5: {acc_thr:.3f} agreement with planted assignment")
            A_cl, params = assign_cleanser(mdl.guide, threshold=None)
            A_cl = dense(A_cl)
            pz_true = A_cl[truth & (W["guide_umi"] > 0)]
            pz_bg = A_cl[~truth & (W["guide_umi"] > 0)]
            check(np.median(pz_true) > 0.9 and np.median(pz_bg) < 0.1 and 0 < A_cl.max() <= 1,
                  f"CLEANSER (CROP-seq, Poisson+NB): posterior PZi median {np.median(pz_true):.3f} in planted cells vs {np.median(pz_bg):.3f} background")
            A_cl8 = dense(assign_cleanser(mdl.guide, threshold=0.8)[0]) > 0
            check(set(np.unique(dense(assign_cleanser(mdl.guide, threshold=0.8)[0]))) <= {0.0, 1.0} and (A_cl8 == truth).mean() > 0.97,
                  f"CLEANSER -t 0.8: binary layer, {(A_cl8 == truth).mean():.3f} agreement")
            p0 = next(iter(params.values()))
            check(0 < p0["lambda"] < 5 and 0 < p0["r"] <= 0.1 and p0["nbMean"] > 10, f"CLEANSER params: lambda {p0['lambda']:.2f}, r {p0['r']:.3f}, nbMean {p0['nbMean']:.1f} within the Stan bounds")
            gdc = mdl.guide.copy()
            rngd = np.random.default_rng(2)
            gdc.X = __import__("scipy.sparse", fromlist=["csr_matrix"]).csr_matrix(
                np.where(truth, rngd.negative_binomial(3, 3 / (3 + 60.0), truth.shape) + 5, rngd.negative_binomial(1.1, 1.1 / (1.1 + 1.0), truth.shape)).astype(float))
            A_dc = dense(assign_cleanser(gdc, threshold=0.5, capture="direct capture")[0]) > 0
            check((A_dc == truth).mean() > 0.95, f"CLEANSER direct capture (NB+NB): {(A_dc == truth).mean():.3f} agreement")
            A_sc = dense(assign_sceptre_mixture(mdl.guide, mdl.gene)) > 0
            check((A_sc == truth).mean() > 0.95, f"sceptre mixture (posterior > 0.8): {(A_sc == truth).mean():.3f} agreement")
            rc = cmd_assign(argparse.Namespace(mudata=str(W["h5mu"]), method="umi-threshold", threshold=5.0, capture_method=None,
                                               out=str(d / "assigned.h5mu"), label="st_j2_assign"))
            last("st_j2_assign")
            check(rc == 0 and "guide_assignment" in read_h5mu(d / "assigned.h5mu").guide.layers, "assign-guides: writes the guide_assignment layer into a readable .h5mu")

            print("\ninference")
            pairs = _pairs(md)
            pos = pairs["pair_type"] == "positive_control"
            neg = pairs["pair_type"] == "negative_control"
            res = {}
            for m in ("wilcoxon", "negbinom", "scanpy-wilcoxon", "scanpy-t-test", "scanpy-t-test-overestim-var", "sceptre"):
                res[m] = run_inference(md, m, B=299)
                tr = res[m]
                check(len(tr) == len(pairs) and tr.loc[pos.values, "p_value"].max() < 1e-4 and (tr.loc[neg.values, "p_value"] < 0.01).mean() < 0.1,
                      f"{m}: planted pairs p <= {tr.loc[pos.values, 'p_value'].max():.1e}; {int((tr.loc[neg.values, 'p_value'] < 0.01).sum())}/{int(neg.sum())} negative controls at p < 0.01")
            nb = res["negbinom"].set_index(["intended_target_name", "gene_id"])
            check(abs(nb.loc[("E1", g1), "log2_fc"] - math.log2(0.3)) < 0.25 and nb.loc[("E4", g4), "log2_fc"] > 1,
                  f"negbinom: log2_fc E1->G1 {nb.loc[('E1', g1), 'log2_fc']:.2f} (true {math.log2(0.3):.2f}); E4->G4 up {nb.loc[('E4', g4), 'log2_fc']:.2f}")
            nbc = infer_negbinom(md, pairs.iloc[:3], compat=True)
            check("l2fc" in nbc.columns and "log2_fc" not in nbc.columns and abs(nbc["l2fc"].iloc[1] - 0.3) < 0.08,
                  "negbinom --upstream-compat: `l2fc` = exp(b1), the fold change upstream stores")
            scp = res["sceptre"].set_index(["intended_target_name", "gene_id"])
            check(abs(scp.loc[(W["genes"][3], g3), "log2_fc"] - math.log2(0.25)) < 0.3 and scp["skew_normal_fit"].mean() > 0.8,
                  f"sceptre: promoter pair log2_fc {scp.loc[(W['genes'][3], g3), 'log2_fc']:.2f} (true -2.00); skew-normal fit used for {scp['skew_normal_fit'].mean():.0%} of pairs")
            null_p = res["sceptre"].loc[neg.values, "p_value"]
            check(0.2 < null_p.median() < 0.8, f"sceptre: negative-control p-values roughly uniform (median {null_p.median():.2f})")
            mdL = read_h5mu(WL["h5mu"])
            wl_res = run_inference(mdL, "wilcoxon")
            ntc = nt_cells(mdL.guide)
            check(wl_res.loc[(wl_res["pair_type"] == "positive_control").values, "p_value"].max() < 1e-3 and 50 < ntc.sum() < 600,
                  f"low MOI: controls are the {int(ntc.sum())} NT cells; planted pairs detected")
            scL = run_inference(mdL, "sceptre", B=299)
            check(scL.loc[(scL["pair_type"] == "positive_control").values, "p_value"].max() < 1e-3, "low MOI sceptre (permutations): planted pairs detected")
            rc = cmd_infer(argparse.Namespace(mudata=str(W["h5mu"]), methods=["wilcoxon", "scanpy-wilcoxon", "perturbo"], side="both", moi=None,
                                              resamples=199, seed=4, upstream_compat=False, out=str(d / "inferred.h5mu"), label="st_j2_infer"))
            r = last("st_j2_infer")
            s = json.loads((r / "summary.json").read_text())
            check(rc == 0 and (r / "test_results.wilcoxon.tsv").is_file() and "test_results" in read_h5mu(d / "inferred.h5mu").uns,
                  f"infer: tables + output MuData with uns['test_results']; perturbo -> {s['methods'].get('perturbo', {}).get('n')}")

            print("\nevaluation")
            e = evaluate_results(res["sceptre"])
            check(e["auroc"] > 0.95 and e["auprc"] > 0.9, f"evaluate: sceptre AUROC {e['auroc']:.3f}, AUPRC {e['auprc']:.3f} (pair_type controls)")
            ec = evaluate_results(res["sceptre"], compat=True)
            check(ec["auroc"] < 0.1, f"evaluate --upstream-compat: raw p_value as the score inverts the ROC (AUROC {ec['auroc']:.3f})")
            M = clustergram_matrix(md)
            check(M.loc["E1_g1", g1] < -0.3 and abs(M.loc["NT_0", g1]) < 0.2 and M.loc["E4_g1", g4] > 0.3,
                  f"clustergram (guide .X != 0, background UMIs included as upstream): E1_g1 x G1 {M.loc['E1_g1', g1]:.2f}, E4_g1 x G4 {M.loc['E4_g1', g4]:.2f}, NT_0 x G1 {M.loc['NT_0', g1]:.2f}")
            rc = cmd_evaluate(argparse.Namespace(mudata=str(d / "inferred.h5mu"), results=None, upstream_compat=False, no_plots=args.no_plots,
                                                 no_clustergram=False, label="st_j2_eval"))
            r = last("st_j2_eval")
            check(rc == 0 and (r / "clustergram_lfc.tsv").is_file(), "evaluate: report + clustergram table")

            print("\nsimulation")
            ptp = pd.DataFrame({"gene_id": [W["genes"][20], W["genes"][21]], "intended_target_name": ["E5", "E4"], "effect_size": [0.5, 0.5]})
            sim, disp = simulate_counts(md, ptp, guide_var=0.1, seed=3)
            base = W["genes"]
            check(np.all(np.isfinite(disp["dispersion"])) and 0.05 < np.median(disp["dispersion"]) < 0.5,
                  f"DESeq2-style dispersions: median {np.median(disp['dispersion']):.3f} (planted NB size 5 -> 0.2)")
            check(abs(np.median(disp["size_factors"]) - 1) < 0.2 and "size_factors" in sim.gene.obs.columns and "disp_outlier_deseq2" in sim.gene.var.columns,
                  "simulate: size_factors in gene.obs; mean / dispersion / disp_outlier_deseq2 in gene.var")
            rc = cmd_simulate(argparse.Namespace(mudata=str(W["h5mu"]), pairs=None, perturb=[f"{base[20]}:E5:0.3", f"{base[21]}:E4:0.4"], effect_size=0.5,
                                                 n_null_pairs=30, guide_var=0.05, seed=5, upstream_compat=False, out=str(d / "sim.h5mu"), label="st_j2_sim"))
            last("st_j2_sim")
            sm = read_h5mu(d / "sim.h5mu")
            ptt = uns_df(sm, "pairs_to_test")
            tr = run_inference(sm, "negbinom")
            ev = evaluate_results(tr)
            check(rc == 0 and len(ptt) == 32 and ev["truth"] == "effect_size != 1" and ev["auroc"] > 0.9,
                  f"simulate + negbinom + evaluate: planted effects recovered (AUROC {ev['auroc']:.3f}, Pearson r {ev.get('pearson_r_log2fc', float('nan')):.2f}, MSE {ev.get('mse_log2fc', float('nan')):.3f})")
            sim_c, _ = simulate_counts(md, ptp, 0.1, 3, compat=True)
            check(np.array_equal(dense(sim_c.gene.X), dense(md.gene.X)), "simulate --upstream-compat: writes the unsimulated counts (upstream writeH5MU(mu) bug)")

            print("\nvolcano / IGV / network")
            rc = cmd_volcano(argparse.Namespace(mudata=str(d / "inferred.h5mu"), results=str(OUT_ROOT / sorted(OUT_ROOT.glob("*_st_j2_infer"))[-1] / "test_results.scanpy-wilcoxon.tsv"),
                                                log2fc_threshold=1.0, sig_threshold=0.01, no_plots=args.no_plots, label="st_j2_volcano"))
            r = last("st_j2_volcano")
            down = pd.read_csv(r / "volcano_down.tsv", sep="\t")
            up = pd.read_csv(r / "volcano_up.tsv", sep="\t")
            bg = pd.read_csv(r / "igv_promoter_bedgraph.tsv", sep="\t")
            bp = pd.read_csv(r / "igv_links_bedpe.tsv", sep="\t")
            check(rc == 0 and set(down["gene_id"]) == {g1, g2, g3} and set(up["gene_id"]) == {g4} and len(bg) == 1 and len(bp) == len(pairs) - 1,
                  f"volcano: down = planted G1/G2/G3, up = G4 (|log2FC| >= 1, p <= 0.01); IGV 1 promoter bedgraph row + {len(bp)} bedpe links")
            rc = cmd_network(argparse.Namespace(mudata=str(d / "inferred.h5mu"), results=str(OUT_ROOT / sorted(OUT_ROOT.glob("*_st_j2_infer"))[-1] / "test_results.scanpy-wilcoxon.tsv"),
                                                central_nodes=["E1"], n_central=2, min_weight=0.01, no_plots=args.no_plots, label="st_j2_net"))
            r = last("st_j2_net")
            edges = pd.read_csv(r / "network_edges.tsv", sep="\t")
            nte = pd.read_csv(r / "mean_nontargeting_expression.tsv", sep="\t")
            check(rc == 0 and set(edges["intended_target_name"]) == {"E1"} and (edges["log2_fc"].abs() >= 0.01).all() and nte["mean_nontargeting_expression"].notna().all(),
                  f"network: {len(edges)} E1 edges with |log2_fc| >= 0.01; mean_nontargeting_expression for all genes")

            print("\nguide reference / counting / seqspec")
            meta = pd.DataFrame({"sgRNA_ID": W["guide_ids"], "sgRNA_sequences": [W["protospacers"][g] for g in W["guide_ids"]]})
            meta_p = d / "guides.tsv"
            meta.to_csv(meta_p, sep="\t", index=False)
            rc = cmd_guide_reference(argparse.Namespace(guides=str(meta_p), nnn_revcomp=False, sa_index_nbases=5, no_star=True, label="st_j2_ref"))
            r = last("st_j2_ref")
            bed = pd.read_csv(r / "guides.bed", sep="\t", header=None)
            fa = (r / "pseudo_genome.fa").read_text().splitlines()
            gtf = (r / "guides.gtf").read_text().splitlines()
            check(rc == 0 and bed.shape == (len(meta), 12) and bed.iloc[0, 0] == f"chr_{meta.sgRNA_ID[0]}_peseudo_chr_{meta.sgRNA_sequences[0]}"
                  and fa[1] == meta.sgRNA_sequences[0] and len(gtf) == 3 * len(meta),
                  "guide-reference: BED12 pseudo-chromosomes, FASTA and bed2gtf GTF (gene/transcript/exon per guide)")
            FQ = synthetic_fastqs(d, W)
            rc = cmd_count_guides(argparse.Namespace(guides=str(meta_p), r1=str(FQ["r1"]), r2=str(FQ["r2"]), whitelist=str(FQ["whitelist"]),
                                                     solo_out=None, cb_start=1, cb_len=16, umi_start=17, umi_len=12, max_reads=None,
                                                     upstream_compat=False, label="st_j2_count"))
            r = last("st_j2_count")
            import anndata as ad  # type: ignore
            a = ad.read_h5ad(r / "guides_star_solo.h5ad")
            got = {(c, g): a[c, g].X.toarray()[0, 0] for (c, g) in FQ["truth"]}
            check(rc == 0 and all(got[k] == v for k, v in FQ["truth"].items()) and a.n_obs == len({k[0] for k in FQ["truth"]}),
                  f"count-guides: {len(FQ['truth'])} planted cell-guide UMI counts recovered exactly (CB 1MM correction, forward + reverse protospacer, UMI 1MM collapsing)")
            rc = cmd_count_guides(argparse.Namespace(guides=str(meta_p), r1=None, r2=None, whitelist=None, solo_out=str(r / "Solo.out"),
                                                     cb_start=1, cb_len=16, umi_start=17, umi_len=12, max_reads=None, upstream_compat=False, label="st_j2_solo"))
            r2 = last("st_j2_solo")
            a2 = ad.read_h5ad(r2 / "guides_star_solo.h5ad")
            check(rc == 0 and a2.shape == a.shape and np.array_equal(a2.X.toarray(), a.X.toarray()), "count-guides --solo-out: Solo.out MatrixMarket parsed back to the same matrix")
            spec_p = d / "spec.yaml"
            spec_p.write_text(SEQSPEC_EXAMPLE)
            idx = seqspec_index(mini_yaml(SEQSPEC_EXAMPLE), "crispr", ["R1", "R2"])
            check(idx["kb"] == "0,0,16:0,16,28:1,63,83" and idx["starsolo"].endswith("--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12")
                  and idx["whitelist"] == "737K-august-2016.txt", f"seqspec-index: kb '{idx['kb']}' (R2 negative strand: cdna after the 63 bp common region)")
            rc = cmd_seqspec_index(argparse.Namespace(spec=str(spec_p), modality="crispr", reads="R1,R2", tool="kb", label="st_j2_seqspec"))
            last("st_j2_seqspec")
            check(rc == 0, "seqspec-index: subcommand")

            print("\ninspect / run")
            rc = cmd_inspect(argparse.Namespace(mudata=str(W["h5mu"]), strict=True, label="st_j2_inspect"))
            r = last("st_j2_inspect")
            s = json.loads((r / "summary.json").read_text())
            check(rc == 0 and all(s["schema"].values()) and s["moi"] == "high" and s["pairs_to_test"]["n"] == len(pairs), "inspect: schema complete, MOI and pairs reported")
            rc = cmd_run(argparse.Namespace(mudata=str(W["h5mu"]), methods=["wilcoxon", "scanpy-wilcoxon", "sceptre"], assign=None, threshold=5.0,
                                            side="both", resamples=199, seed=4, upstream_compat=False, log2fc_threshold=1.0, sig_threshold=0.01,
                                            no_plots=args.no_plots, out=None, label="st_j2_run"))
            r = last("st_j2_run")
            s = json.loads((r / "summary.json").read_text())
            check(rc == 0 and all(s["methods"][m]["auroc"] > 0.9 for m in s["methods"]) and (r / "igv_links_bedpe.tsv").is_file(),
                  "run: infer x3 -> evaluate (AUROC > 0.9 each) -> volcano/IGV -> network")
            if not args.no_plots:
                check(len(s["figures"]) >= 2, f"run: {len(s['figures'])} figures")
    finally:
        for p in made + sorted(OUT_ROOT.glob("*_st_j2_*")):
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
    ap = argparse.ArgumentParser(prog="igvfagent crispr-jamboree2", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("inspect", help="MuData summary + jamboree schema check")
    p.add_argument("--mudata", required=True)
    p.add_argument("--strict", action="store_true", help="exit 1 when a required field is missing")
    p.add_argument("--label", default="inspect")
    p.set_defaults(func=cmd_inspect)

    p = sub.add_parser("guide-reference", help="guide metadata -> pseudo-genome FASTA / BED12 / isoforms / GTF (+ STAR index)")
    p.add_argument("--guides", required=True, help="guide metadata (xlsx/tsv/csv): sgRNA_ID + sgRNA_sequences, or guide_id + protospacer")
    p.add_argument("--nnn-revcomp", action="store_true", help="use 'NNN' + reverse complement of each protospacer (CountGuides.ipynb)")
    p.add_argument("--sa-index-nbases", type=int, default=5)
    p.add_argument("--no-star", action="store_true")
    p.add_argument("--label", default="guide_reference")
    p.set_defaults(func=cmd_guide_reference)

    p = sub.add_parser("count-guides", help="FASTQs -> cells x guides UMI counts (or parse Solo.out)")
    p.add_argument("--guides", required=True)
    p.add_argument("--r1")
    p.add_argument("--r2")
    p.add_argument("--whitelist")
    p.add_argument("--solo-out", help="existing STARsolo Solo.out directory to parse instead")
    p.add_argument("--cb-start", type=int, default=1)
    p.add_argument("--cb-len", type=int, default=16)
    p.add_argument("--umi-start", type=int, default=17)
    p.add_argument("--umi-len", type=int, default=12)
    p.add_argument("--max-reads", type=int)
    p.add_argument("--upstream-compat", action="store_true", help="guides x cells orientation, as soloOutGuideParsing.py writes")
    p.add_argument("--label", default="count_guides")
    p.set_defaults(func=cmd_count_guides)

    p = sub.add_parser("seqspec-index", help="seqspec YAML -> kb -x string / STARsolo flags")
    p.add_argument("--spec", required=True)
    p.add_argument("--modality", required=True)
    p.add_argument("--reads", required=True, help="comma-separated read ids, e.g. R1,R2")
    p.add_argument("--tool", choices=["kb", "starsolo"], default="kb")
    p.add_argument("--label", default="seqspec")
    p.set_defaults(func=cmd_seqspec_index)

    p = sub.add_parser("assign-guides", help="guide_assignment layer by UMI threshold, CLEANSER or sceptre mixture")
    p.add_argument("--mudata", required=True)
    p.add_argument("--method", choices=["umi-threshold", "cleanser", "sceptre"], default="umi-threshold")
    p.add_argument("--threshold", type=float, help="UMI threshold (default 5) or CLEANSER posterior threshold (default: store posteriors)")
    p.add_argument("--capture-method", help="override guide.uns capture_method (CROP-seq | direct capture)")
    p.add_argument("--out")
    p.add_argument("--label", default="assign_guides")
    p.set_defaults(func=cmd_assign)

    p = sub.add_parser("infer", help="test pairs_to_test with the jamboree inference methods")
    p.add_argument("--mudata", required=True)
    p.add_argument("--methods", nargs="+", choices=INFER_METHODS)
    p.add_argument("--side", choices=["left", "right", "both"], default="both")
    p.add_argument("--moi", choices=["low", "high"])
    p.add_argument("--resamples", type=int, default=499, help="sceptre resamples")
    p.add_argument("--seed", type=int, default=4)
    p.add_argument("--upstream-compat", action="store_true", help="negbinom stores `l2fc` = exp(b1) as upstream")
    p.add_argument("--out", help="output .h5mu with uns['test_results'] (first method)")
    p.add_argument("--label", default="infer")
    p.set_defaults(func=cmd_infer)

    p = sub.add_parser("simulate", help="DESeq2-style NB fit + simulated counts with planted effects")
    p.add_argument("--mudata", required=True)
    p.add_argument("--pairs", help="TSV gene_id, intended_target_name[, effect_size]")
    p.add_argument("--perturb", nargs="+", help="gene_id:intended_target_name:effect_size")
    p.add_argument("--effect-size", type=float, default=0.5)
    p.add_argument("--n-null-pairs", type=int, default=0, help="add random pairs with effect_size 1 (for evaluation)")
    p.add_argument("--guide-var", type=float, default=0.1)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--upstream-compat", action="store_true", help="write the unsimulated counts, as the upstream script does")
    p.add_argument("--out")
    p.add_argument("--label", default="simulate")
    p.set_defaults(func=cmd_simulate)

    p = sub.add_parser("evaluate", help="AUPRC / AUROC / Pearson / MSE + clustergram")
    p.add_argument("--mudata")
    p.add_argument("--results", nargs="+", help="test_results TSV/CSV files (default: mudata uns['test_results'])")
    p.add_argument("--upstream-compat", action="store_true", help="score by the raw p-value, as the notebook")
    p.add_argument("--no-clustergram", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="evaluate")
    p.set_defaults(func=cmd_evaluate)

    p = sub.add_parser("volcano", help="volcano tables/figure + IGV bedpe / bedgraph")
    p.add_argument("--mudata", required=True)
    p.add_argument("--results")
    p.add_argument("--log2fc-threshold", type=float, default=2.0)
    p.add_argument("--sig-threshold", type=float, default=0.01)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="volcano")
    p.set_defaults(func=cmd_volcano)

    p = sub.add_parser("network", help="element -> gene network edges + NT mean expression")
    p.add_argument("--mudata", required=True)
    p.add_argument("--results")
    p.add_argument("--central-nodes", nargs="+")
    p.add_argument("--n-central", type=int, default=2)
    p.add_argument("--min-weight", type=float, default=0.01)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="network")
    p.set_defaults(func=cmd_network)

    p = sub.add_parser("run", help="assign -> infer -> evaluate -> volcano/IGV -> network")
    p.add_argument("--mudata", required=True)
    p.add_argument("--methods", nargs="+", choices=INFER_METHODS, default=["wilcoxon", "negbinom", "scanpy-wilcoxon", "sceptre"])
    p.add_argument("--assign", choices=["umi-threshold", "cleanser", "sceptre"], help="(re)assign guides first")
    p.add_argument("--threshold", type=float, default=5.0)
    p.add_argument("--side", choices=["left", "right", "both"], default="both")
    p.add_argument("--resamples", type=int, default=499)
    p.add_argument("--seed", type=int, default=4)
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--log2fc-threshold", type=float, default=2.0)
    p.add_argument("--sig-threshold", type=float, default=0.01)
    p.add_argument("--out")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="jamboree2_run")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("selftest", help="synthetic MuData + FASTQs, all subcommands, assertions")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
