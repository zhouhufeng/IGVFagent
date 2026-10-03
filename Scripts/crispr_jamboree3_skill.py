#!/usr/bin/env python3
"""3rd IGVF CRISPR Jamboree single-cell Perturb-seq pipeline, stage by stage (port of IGVF-CRISPR/CRISPR-jamboree3).

Port of https://github.com/IGVF-CRISPR/CRISPR-jamboree3 (No LICENSE file;
Nextflow DSL2 workflow IGVF_Workflows/ with 45 Python / R helper scripts,
two Stan models, a configuration web form and two tutorial notebooks).  Every
script, process and notebook was read and the logic re-derived in Python
(numpy / scipy / pandas / anndata / h5py); no code was copied.
Relationship: port.  The pipeline is the September-2024 snapshot of the IGVF
single-cell CRISPR pipeline used at the jamboree: seqspec checks -> kallisto |
bustools mapping of RNA, guide and hashing libraries -> AnnData concat -> QC
filtering -> MuData -> Scrublet doublets -> GMM-demux hashing -> guide
assignment (CLEANSER / sceptre) -> pairs to test -> sceptre / PerTurbo
inference -> volcano / network / IGV evaluation -> HTML dashboard.  The binary
steps (kb ref / kb count, GMM-Demux, the sceptre and PerTurbo packages, cmdstan)
are re-implemented or wrapped: `map` writes (and, with kb on PATH, runs) the
exact kb commands.  Not vendored: the 12 MB 10x whitelist, the example
spreadsheets and the compiled cmdstan binaries (point the subcommands at them).

Definitions, exactly as upstream computes them
  configure       parse_interface_configuration.py: configuration.csv rows
                  (tab_name, variable, variable_value, batch_name, read1, read2)
                  -> params { DATASET_HASHING, non-sequence parameters (numeric
                  when int()/float() parse, '.' decides float, else quoted),
                  fastq_files_rna/guide/hashing grouped per batch "r1 r2 ...",
                  test_*_fastq_r1/r2, batch=[...] } + params.covariate_list
                  (batch_name split on ', ' -> batch, cov1, ...);
                  distance_from_center = 1000000 added when absent.
                  prepare_formula.py: covariates with > 1 level joined " + ".
                  On the upstream Workflow_configuration/configuration.csv the
                  port writes a file byte-identical to the shipped
                  IGVF_Workflows/configs/pipeline_input.config.
  seqspec-check   seqSpecCheck.py: first max_reads (100000) reads of R1/R2,
                  per-position A/C/G/T frequency; for R2 every metadata
                  sequence (2nd column) found in a read adds read.index(seq)
                  to position_table.csv ('read name' row per file).
  seqspec-parse   parsing_guide_metadata.py: reads of the modality from
                  `seqspec info`, `seqspec index -m <mod> -t kb -r <reads>`
                  -> representation, whitelist = first `filename:` line;
                  parsed_seqSpec.txt (modality, representation,
                  barcode_whitelist, seqspec_dir).
  map             kb ref (-d human | kite: guide_features.txt = sgRNA_sequences
                  + sgRNA_ID, hashing_table.txt = HTO_sequences + HTO_ID) and
                  kb count --h5ad -x <representation> -o <batch>_ks_<mod>_out
                  (-m 30G RNA, 20G guide / hashing).
  concat          anndata_concat.py: <batch>_ks_*/counts_unfiltered/adata.h5ad
                  sorted by name, batch from the '(.+)_ks_' prefix joined to the
                  parsed covariate table, ad.concat(join='outer',
                  index_unique='_').
  preprocess      preprocess_adata.py: var symbol from
                  cells_x_genes.genes.names.txt, Ensembl version stripped,
                  knee plot, batch_number = factorize(batch) + 1, mt = symbol
                  startswith 'MT-', ribo = RPS/RPL, sc.pp.calculate_qc_metrics
                  (qc_vars mt, ribo; log1p), filter_cells(min_genes 500),
                  filter_genes(min_cells 3), pct_counts_mt < pct_mito (20).
                  Upstream computes mt_prefix ('Mt-' for mouse) but never uses
                  it; the port uses it, --upstream-compat always uses 'MT-'.
  create-mudata   create_mdata(_HASHING).py: guide var merged with the guide
                  metadata on guide_id == sgRNA_ID; intended_target_name,
                  intended_target_chr/start/end, sequence = Target_name, chr,
                  start, end, sgRNA_sequences; var_names = sgRNA_ID|sequence;
                  uns moi, capture_method CROP-seq; obs num_expressed_guides,
                  batch_number, total_guide_umis; gene var gene_chr/start/end
                  from the GTF (version-stripped gene_id); gene obs renamed
                  n_genes_by_counts -> n_counts, pct_counts_mt -> percent_mito,
                  n_genes -> num_expressed_genes, total_counts ->
                  total_gene_umis; barcodes intersected across modalities;
                  global obs = shared obs columns.  Upstream quirks corrected
                  (reproduced with --upstream-compat): metadata columns are
                  assigned positionally after the merge (misaligned when the
                  order differs), and `targeting` is 'TRUE' for every guide,
                  including non-targeting ones.
  doublets        doublets.py: Scrublet defaults (sim_doublet_ratio 2,
                  expected_doublet_rate 0.1, n_neighbors round(0.5 sqrt(n)),
                  min_counts 3, min_cells 3, gene variability pctl 85, z-score,
                  PCA 30), score q*rho/r/(1-rho-q(1-rho-rho/r)) with
                  q = (n_sim_neighbours+1)/(k_adj+2), threshold_minimum of the
                  simulated-doublet score histogram; doublets removed;
                  doublet_scores / predicted_doublets / doublet_info kept.
                  The v-score gene filter is approximated by the
                  mean-adjusted Fano factor (documented deviation).
  demultiplex     GMM-Demux + demultiplex_filter.py + filter_hashing.py:
                  per-HTO two-component Gaussian mixture on log1p counts,
                  cluster = set of positive HTOs ('negative' when none),
                  hto_type = names joined '-', hto_type_split = 'multiplets'
                  when > 1 HTO; negatives dropped; hashing barcodes
                  intersected with the filtered RNA barcodes, per batch.
  assign-guides   cleanser.py: CLEANSER Stan mixtures (cs: Poisson + NB,
                  dc: NB + NB; priors as the .stan files) -> PZi; with -t the
                  layer is 1 where PZi >= t (THRESHOLD 1 in the example
                  config).  Upstream fits ONE model to all nonzero counts
                  pooled over guides and indexes PZi by the within-guide
                  position; the port fits per guide (as CLEANSER itself) and
                  reproduces the pooled indexing with --upstream-compat.
                  assign_grnas_sceptre.R: sceptre mixture (two-component
                  Poisson GLM per gRNA, posterior > 0.8); umi-threshold >= t.
  prepare-inference prepare_inference.py: user pairs (guide_id, gene_name,
                  intended_target_name, pair_type) restricted to genes in the
                  gene modality and guides in guide.var.guide_id, gene_name ->
                  gene_id; create_pairs_to_test.py: every GTF gene on the
                  guide's chromosome with |gene start - guide start| <=
                  distance_from_center (1 Mb), pair_type 'discovery'.
                  Upstream writes gene *symbols* in gene_name, which never
                  match the Ensembl gene ids of the gene modality; the port
                  writes the Ensembl id (symbols with --upstream-compat).
  infer           inference_sceptre.R: sceptre with formula ~ log(
                  response_n_nonzero) + log(response_n_umis) (+ the
                  non-redundant covariates of cov_string when formula_object
                  is not 'default'), assign_grnas(thresholding, 1), QC off,
                  side both, union, skew-normal; results left-joined to
                  pairs_to_test by (intended_target_name, gene_id); log2_fc =
                  log_2_fold_change.  identify_non_redundant_covariates: keep
                  one of each group of covariates with the same level
                  structure, drop single-level ones.  Re-implemented as in
                  crispr-jamboree2 (NB score test + conditional randomisation
                  or permutations, skew-normal fit).  perturbo_inference.py
                  runs when `perturbo` is importable (batch_key 'batch').
  evaluate        volcano_plot.py (log2_fc 1, p 0.05 from evaluation_plot.nf;
                  annotate the 10 most down-regulated by log2_fc and by p);
                  select_nodes.py (first N intended targets by ascending
                  p_value); network_plot.py (edges with |log2_fc| >= 0.1, node
                  size p_value as upstream); igv.py (promoter bedgraph when
                  intended_target_name == GTF gene_name of gene_id, else bedpe;
                  headerless files).
  dashboard       create_dashboard_plots.py thresholds (RNA UMIs 200, 500,
                  1000, 2000, 5000; detected genes 200, 500, 1000, 1500, 2000;
                  guide UMIs 0, 5, 10, 50, 100, 1000; counts are cells with
                  sum > threshold), sgRNA frequencies (cells with assignment >
                  0), guides per cell / cells per guide (sums of guide .X, as
                  upstream); create_dashboard_df.py highlight blocks (barcode
                  intersections, cells after filtering, genes, mean UMIs,
                  significant pairs at p < 0.05 / 0.01, 'Direct targeting' and
                  'Targeting_negative_control' counts; upstream divides the
                  direct-targeting count by all tested pairs, the port reports
                  that and the fraction of direct-targeting pairs);
                  process_json.py kb inspect.json + run_info.json tables;
                  self-contained dashboard.html.

Subcommands
  configure, seqspec-check, seqspec-parse, map, concat, preprocess,
  create-mudata, doublets, demultiplex, assign-guides, prepare-inference,
  infer, evaluate, dashboard  -- one per pipeline stage (above)
  run        preprocess -> [demultiplex] -> create-mudata -> doublets ->
             assign-guides -> prepare-inference -> infer -> evaluate ->
             dashboard, from RNA / guide [/ hashing] count AnnDatas
  selftest   synthetic FASTQs, config table, count matrices, hashing and
             guide metadata with planted signal; all subcommands asserted

Output: Docs/CRISPRJamboree3/<timestamp>_<label>/.

Usage:
    igvfagent crispr-jamboree3 configure --config-table configuration.csv
    igvfagent crispr-jamboree3 seqspec-check --read1 R1.fastq.gz --read2 R2.fastq.gz --metadata guide_metadata.xlsx
    igvfagent crispr-jamboree3 seqspec-parse --yaml multiseq_guide_utsw.yml --modality guide
    igvfagent crispr-jamboree3 map --modality guide --seqspec-yaml multiseq_guide_utsw.yml --metadata guide_metadata.xlsx --fastqs "batch_a:R1.fq.gz R2.fq.gz"
    igvfagent crispr-jamboree3 preprocess --adata rna_concat.h5ad --min-genes 500 --pct-mito 20
    igvfagent crispr-jamboree3 create-mudata --rna filtered_anndata.h5ad --guide guide_concat.h5ad --guide-metadata guide_metadata.xlsx --gtf gencode.v46.annotation.gtf.gz
    igvfagent crispr-jamboree3 run --rna rna.h5ad --guide guide.h5ad --hashing hashing.h5ad --guide-metadata guides.xlsx --gtf genes.gtf.gz --pairs user_pairs_to_test.csv
    igvfagent crispr-jamboree3 selftest --no-plots
"""
from __future__ import annotations

import argparse
import base64
import gzip
import json
import logging
import math
import os
import re
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
OUT_ROOT = DOCS_DIR / "CRISPRJamboree3"

UPSTREAM_REPO = "IGVF-CRISPR/CRISPR-jamboree3"
UPSTREAM_COMMIT = "fea0857ff461889bbf855d61c529977fe8bc71cc"

INK, INK2, AXIS, SURFACE = "#0b0b0b", "#52514e", "#d0cfca", "#fcfcfb"
BLUE, ORANGE, GREEN, RED, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#e34948", "#96a0b3"

log = logging.getLogger("crispr_jamboree3")


# ---------------------------------------------------------------------------
# Plumbing and MuData I/O (h5py + anndata; the `mudata` package is optional)
# ---------------------------------------------------------------------------

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    path = LOG_DIR / f"crispr_jamboree3_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
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


# ---------------------------------------------------------------------------
# Statistics: NB GLM, sceptre-style test, CLEANSER, sceptre mixture (shared with crispr-jamboree2)
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


# ---------------------------------------------------------------------------
# seqspec + FASTQ helpers
# ---------------------------------------------------------------------------

def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


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


# ---------------------------------------------------------------------------
# Result-table helpers
# ---------------------------------------------------------------------------

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


def read_ad(path: Path):
    import anndata as ad  # type: ignore
    return ad.read_h5ad(str(path))


def write_ad(a, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    for df in (a.obs, a.var):
        for c in df.columns:
            if df[c].dtype == object:
                df[c] = df[c].astype(str)
    a.write_h5ad(str(path))
    print(f"Wrote: {path}")
    return path


def strip_version(ids):
    return [str(x).split(".")[0] for x in ids]


# ---------------------------------------------------------------------------
# configure (parse_interface_configuration.py, parse_covariate.py, prepare_formula.py, process_batches.py, process_reads.py)
# ---------------------------------------------------------------------------

def _config_value(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        v = ""
    v = str(v)
    try:
        return float(v) if "." in v else int(v)
    except ValueError:
        return f"'{v}'"


def _format_fastq_files(df, tab_name: str) -> str:
    df = df.dropna(subset=["batch_name"])
    grouped = [" ".join(f"{r1} {r2}" for r1, r2 in zip(x["read1"], x["read2"])) for _, x in df.groupby("batch_name")]
    key = {"scRNA": "fastq_files_rna", "Guides": "fastq_files_guide", "Hash": "fastq_files_hashing"}.get(tab_name, "fastq_files")
    return f"{key} = [\n    " + ",\n    ".join(f'"{s}"' for s in grouped) + "\n    ]"


def _format_test_fastq_files(df, tab_name: str) -> str:
    df = df.dropna(subset=["read1", "read2"])
    r1 = [str(x).strip() for x in df["read1"]]
    r2 = [str(x).strip() for x in df["read2"]]
    keys = {"Guides": ("test_guide_fastq_r1", "test_guide_fastq_r2"), "Hash": ("test_hashing_fastq_r1", "test_hashing_fastq_r2")}.get(tab_name)
    if keys is None:
        return ""
    return f"{keys[0]} = {r1}\n    {keys[1]} = {r2}"


def generate_pipeline_config(table, mapping_tabs=("scRNA", "Guides", "Hash")) -> "tuple[str, dict, dict]":
    """pipeline_input.config text, covariate_list dict and a parameter dict."""
    pd = _pd()
    seqs, tests, hashing = {}, [], "false"
    for tab in mapping_tabs:
        f = table[table["tab_name"] == tab]
        f = f[f["read1"].notna() & f["read2"].notna()] if "read1" in f.columns else f
        if not f.empty:
            seqs[tab] = _format_fastq_files(f, tab)
            if tab in ("Guides", "Hash"):
                tests.append(_format_test_fastq_files(f, tab))
            if tab == "Hash":
                hashing = "true"
    nonseq = table[table["variable"] != "sequence"][["variable", "variable_value"]]
    lines, params, have_dist = [], {}, False
    for _, r in nonseq.iterrows():
        v = _config_value(r["variable_value"])
        if r["variable"] == "distance_from_center":
            have_dist = True
        lines.append(f"{r['variable']} = {v}")
        params[r["variable"]] = v.strip("'") if isinstance(v, str) else v
    if not have_dist:
        lines.append("distance_from_center = 1000000")
        params["distance_from_center"] = 1000000
    nn = table[table["batch_name"].notna()]
    split = nn["batch_name"].apply(lambda x: str(x).split(", "))
    covariate_dict: dict = {}
    batch_only = "batch=[]"
    cov_text = "params.covariate_list = [\n    batch: [],\n".rstrip(",\n") + "\n]"
    if len(split):
        mx = int(split.apply(len).max())
        cols = ["batch"] + [f"cov{i}" for i in range(1, mx)]
        batch = split.apply(lambda x: x[0])
        cov = pd.DataFrame(split.tolist(), columns=cols).drop_duplicates()
        others = cov.drop(columns=["batch"]).to_dict(orient="list")
        cov_text = f"params.covariate_list = [\n    batch: {batch.unique().tolist()},\n"
        batch_only = f"batch={batch.unique().tolist()}"
        for k, v in others.items():
            cov_text += f"    {k}: {v},\n"
        cov_text = cov_text.rstrip(",\n") + "\n]"
        covariate_dict = {"batch": batch.unique().tolist(), **others}
    text = (f"params {{\n    DATASET_HASHING = '{hashing}'\n\n    " + "\n    ".join(lines) + "\n\n    "
            + "".join(f"{c}\n    " for c in seqs.values()) + "\n\n    " + "\n\n    ".join(tests) + f"\n    {batch_only}\n}}\n\n{cov_text}")
    params["DATASET_HASHING"] = hashing
    return text, covariate_dict, params


def prepare_formula(cov_df) -> str:
    return " + ".join(c for c in cov_df.columns if len(set(cov_df[c])) > 1)


def process_batches(fastq: str, directory: str) -> "list[str]":
    out = []
    for b in fastq.split(";"):
        s = " ".join(b.replace("|", " ").replace(",", " ").split())
        out.append(" ".join(os.path.join(directory, f) for f in s.split()))
    return out


def process_reads(directory: str, reads1: str, reads2: str) -> str:
    r1 = re.sub(r"\s", "", reads1).split(";")
    r2 = re.sub(r"\s", "", reads2).split(";")
    return " ".join(f"{os.path.join(directory, a)} {os.path.join(directory, b)}" for a, b in zip(r1, r2))


def cmd_configure(args: argparse.Namespace) -> int:
    pd = _pd()
    d = run_dir(args.label)
    table = pd.read_csv(args.config_table)
    text, cov, params = generate_pipeline_config(table, args.mapping_tabs)
    p = d / "pipeline_input.config"
    p.write_text(text)
    print(f"Wrote: {p}")
    cov_df = pd.DataFrame(cov) if cov else pd.DataFrame()
    write_tsv(cov_df, d / "parse_covariate.csv")
    formula = prepare_formula(cov_df) if len(cov_df) else ""
    (d / "cov_string.txt").write_text(formula)
    print(f"Wrote: {d / 'cov_string.txt'}")
    write_report(d, "CRISPR Jamboree 3 — pipeline configuration",
                 [f"DATASET_HASHING = {params['DATASET_HASHING']}; {len(params)} parameters; covariate formula `{formula}`.",
                  f"```\n{text}\n```"], {"params": params, "covariate_list": cov, "cov_string": formula, "config": p})
    return 0


# ---------------------------------------------------------------------------
# seqspec checks (seqSpecCheck.py, parsing_guide_metadata.py, extract_parsed_seqspec.py)
# ---------------------------------------------------------------------------

def nucleotide_frequencies(seqs):
    np, pd = _np(), _pd()
    L = max((len(s) for s in seqs), default=0)
    arr = np.array([list(s.ljust(L, "-")) for s in seqs]) if seqs else np.empty((0, 0))
    return pd.DataFrame({b: (arr == b).sum(axis=0) / max(len(seqs), 1) for b in "ACGT"}) if L else pd.DataFrame(columns=list("ACGT"))


def find_sequence_positions(seqs, guide_seqs) -> Counter:
    """seqSpecCheck.py: every metadata sequence contained in a read adds read.index(seq)."""
    pos = Counter()
    for s in seqs:
        for n in guide_seqs:
            if n in s:
                pos[s.index(n)] += 1
    return pos


def cmd_seqspec_check(args: argparse.Namespace) -> int:
    pd = _pd()
    d = run_dir(args.label)
    if len(args.read1) != len(args.read2):
        raise SystemExit("The number of read1 and read2 files must be the same.")
    meta = read_table(Path(args.metadata))
    guide_seqs = [str(x) for x in meta.iloc[:, 1].values]
    plt = None if args.no_plots else _plt()
    fig = axs = None
    if plt is not None:
        fig, axs = plt.subplots(len(args.read1), 2, figsize=(12, 3 * len(args.read1)), squeeze=False)
    tables, summary = [], {}
    for i, (r1, r2) in enumerate(zip(args.read1, args.read2)):
        n1 = os.path.basename(r1).replace(".fastq.gz", "")
        n2 = os.path.basename(r2).replace(".fastq.gz", "")
        s1 = list(read_fastq_seqs(Path(r1), args.max_reads))
        s2 = list(read_fastq_seqs(Path(r2), args.max_reads))
        f1, f2 = nucleotide_frequencies(s1), nucleotide_frequencies(s2)
        write_tsv(f1.rename_axis("position").reset_index(), d / f"nucleotide_frequency.{n1}.tsv")
        write_tsv(f2.rename_axis("position").reset_index(), d / f"nucleotide_frequency.{n2}.tsv")
        pos = find_sequence_positions(s2, guide_seqs)
        t = pd.DataFrame(sorted(pos.items(), key=lambda kv: -kv[1]), columns=["Position Index", "Count"])
        tables.append(pd.concat([pd.DataFrame({"Position Index": ["read name"], "Count": [n2]}), t], ignore_index=True))
        summary[n2] = {"reads": len(s2), "reads_with_guide": sum(pos.values()), "modal_position": int(t.iloc[0, 0]) if len(t) else None}
        if axs is not None:
            for j, (f, nm) in enumerate(((f1, n1), (f2, n2))):
                for b, c in zip("ACGT", (RED, GREEN, BLUE, ORANGE)):
                    axs[i, j].plot(f.index, f[b], color=c, label=b, lw=1)
                _style(axs[i, j], f"Nucleotide Frequency Plot - {nm}", "Position in Sequence", "Frequency")
                axs[i, j].legend(frameon=False, fontsize=7)
    tab = pd.concat(tables, ignore_index=True)
    tab.to_csv(d / "position_table.csv", index=False)
    print(f"CSV: {d / 'position_table.csv'}")
    figs = []
    if fig is not None:
        (d / "seqSpec_plots").mkdir(exist_ok=True)
        figs.append(_save(fig, d / "seqSpec_plots" / "seqSpec_check_plots.png"))
    write_report(d, "CRISPR Jamboree 3 — seqSpec check", [md_table(["read 2", "reads", "reads with a guide", "modal position"],
                                                                  [[k, v["reads"], v["reads_with_guide"], v["modal_position"]] for k, v in summary.items()])],
                 {"files": summary, "figures": figs, "position_table": d / "position_table.csv"})
    return 0


def seqspec_reads_for_modality(spec: dict, modality: str) -> "list[str]":
    return [r["read_id"] for r in (spec.get("sequence_spec") or []) if r.get("modality") == modality]


def whitelist_filename(yaml_text: str) -> Optional[str]:
    for line in yaml_text.splitlines():
        if line.strip().startswith("filename"):
            parts = line.strip().split(":", 1)
            if len(parts) > 1:
                return parts[1].strip()
    return None


def parse_seqspec(yaml_path: Path, modalities: "list[str]", directory: str):
    pd = _pd()
    text = yaml_path.read_text()
    spec = mini_yaml(text)
    rows = []
    for m in modalities:
        reads = seqspec_reads_for_modality(spec, m)
        rep = seqspec_index(spec, m, reads)["kb"] if reads else ""
        rows.append({"modality": m, "representation": rep, "barcode_whitelist": whitelist_filename(text), "seqspec_dir": directory})
    return pd.DataFrame(rows)


def cmd_seqspec_parse(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    df = parse_seqspec(Path(args.yaml), args.modality, args.directory or str(Path(args.yaml).parent))
    p = d / f"{'_'.join(args.modality)}_parsed_seqSpec.txt"
    df.to_csv(p, sep="\t", index=False)
    print(f"TSV: {p}")
    rep = "".join(df["representation"].astype(str))
    print(f"chemistry (-x): {rep}")
    write_report(d, "CRISPR Jamboree 3 — seqspec parse", [md_table(list(df.columns), df.values.tolist())],
                 {"parsed": df.to_dict(orient="records"), "chemistry": rep, "file": p})
    return 0


# ---------------------------------------------------------------------------
# map (createGuideRef / createHashingRef / downloadReference / mapping*.nf, guide_table.py, hashing_table.py)
# ---------------------------------------------------------------------------

def feature_table(meta, modality: str, path: Path) -> Path:
    cols = ["sgRNA_sequences", "sgRNA_ID"] if modality == "guide" else ["HTO_sequences", "HTO_ID"]
    missing = [c for c in cols if c not in meta.columns]
    if missing:
        raise SystemExit(f"metadata lacks {missing}")
    meta[cols].to_csv(path, sep="\t", header=False, index=False)
    print(f"Wrote: {path}")
    return path


def kb_commands(modality: str, batches: "list[tuple[str, str]]", chemistry: str, whitelist: str, genome: str = "genome.fa.gz",
                feature_file: Optional[str] = None, species: str = "human", threads: int = 4) -> "list[str]":
    suffix = {"rna": "transcripts", "guide": "guide", "hashing": "hashing"}[modality]
    if modality == "rna":
        idx, t2g = "transcriptome_index.idx", "transcriptome_t2g.txt"
        cmds = [f"kb ref -d {species} -i {idx} -g {t2g} --kallisto $(which kallisto)"]
        mem = "30G"
    else:
        idx = "guide_index.idx" if modality == "guide" else "hashing_index.idx"
        t2g = "t2guide.txt" if modality == "guide" else "t2g_hashing.txt"
        cmds = [f"kb ref -i {idx} -f1 {genome} -g {t2g} --kallisto $(which kallisto) --bustools $(which bustools) --workflow kite {feature_file}"]
        mem = "20G"
    for batch, fq in batches:
        cmds.append(f"kb count -i {idx} -g {t2g} --verbose -w {whitelist} --h5ad --kallisto $(which kallisto) --bustools $(which bustools) "
                    f"-x {chemistry} -o {batch}_ks_{suffix}_out -t {threads} {fq} --overwrite -m {mem}")
    return cmds


def cmd_map(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    chem = args.chemistry
    wl = args.whitelist
    if args.seqspec_yaml:
        text = Path(args.seqspec_yaml).read_text()
        spec = mini_yaml(text)
        mod = {"rna": "rna", "guide": "guide", "hashing": "hashing"}[args.modality]
        reads = seqspec_reads_for_modality(spec, mod)
        chem = chem or seqspec_index(spec, mod, reads)["kb"]
        wl = wl or whitelist_filename(text)
    if not chem:
        raise SystemExit("--chemistry or --seqspec-yaml is required")
    ff = None
    if args.modality != "rna":
        if not args.metadata:
            raise SystemExit("--metadata (guide or hashing spreadsheet) is required for guide / hashing mapping")
        ff = str(feature_table(read_table(Path(args.metadata)), args.modality, d / ("guide_features.txt" if args.modality == "guide" else "hashing_table.txt")))
    batches = []
    for i, s in enumerate(args.fastqs or []):
        b, _, fq = s.partition(":") if ":" in s.split()[0] else (f"batch_{i + 1}", "", s)
        batches.append((b, fq))
    cmds = kb_commands(args.modality, batches, chem, wl or "whitelist.txt", args.genome, ff, args.species, args.threads)
    sh = d / "map_commands.sh"
    sh.write_text("#!/bin/bash\nset -euo pipefail\n" + "\n".join(cmds) + "\n")
    print(f"Wrote: {sh}")
    ran = False
    if args.execute and shutil.which("kb"):
        ran = subprocess.run(["bash", str(sh)], cwd=str(d)).returncode == 0
    elif args.execute:
        print("kb (kallisto | bustools) is not on PATH; commands written, not run.")
    for c in cmds:
        print(f"Command: {c}")
    write_report(d, f"CRISPR Jamboree 3 — kb mapping ({args.modality})", [f"chemistry `-x {chem}`, whitelist `{wl}`; ran: {ran}",
                                                                          "```\n" + "\n".join(cmds) + "\n```"],
                 {"modality": args.modality, "chemistry": chem, "whitelist": wl, "commands": cmds, "ran": ran})
    return 0


# ---------------------------------------------------------------------------
# concat (anndata_concat.py, hashing_concat.py)
# ---------------------------------------------------------------------------

def concat_batches(dirs: "list[Path]", covariates=None):
    import anndata as ad  # type: ignore
    pd = _pd()
    adatas = []
    for p in sorted(dirs, key=lambda x: Path(x).name):
        p = Path(p)
        f = p / "counts_unfiltered" / "adata.h5ad" if p.is_dir() else p
        a = ad.read_h5ad(str(f))
        m = re.search(r"(.+)_ks_", p.name)
        if m and covariates is not None and len(covariates.columns):
            c1 = covariates.columns[0]
            a.obs[c1] = m.group(1)
            cov = covariates.copy()
            cov[c1] = cov[c1].astype(str)
            a.obs = a.obs.join(cov.drop_duplicates(c1).set_index(c1), on=c1)
        elif m:
            a.obs["batch"] = m.group(1)
        adatas.append(a)
    return ad.concat(adatas, join="outer", index_unique="_")


def cmd_concat(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    cov = read_table(Path(args.covariates)) if args.covariates else None
    a = concat_batches([Path(x) for x in args.inputs], cov)
    out = Path(args.out) if args.out else d / "concatenated_adata.h5ad"
    write_ad(a, out)
    write_report(d, "CRISPR Jamboree 3 — AnnData concat", [f"{len(args.inputs)} batches -> {a.n_obs} barcodes x {a.n_vars} features."],
                 {"n_obs": a.n_obs, "n_vars": a.n_vars, "out": out, "obs_columns": list(a.obs.columns)})
    return 0


# ---------------------------------------------------------------------------
# preprocess (preprocess_adata.py)
# ---------------------------------------------------------------------------

def calculate_qc_metrics(a, qc_vars=("mt", "ribo")) -> None:
    """sc.pp.calculate_qc_metrics(qc_vars, log1p=True, inplace=True) columns."""
    np = _np()
    X = a.X
    tot = np.asarray(X.sum(axis=1)).ravel().astype(float)
    ngenes = np.asarray((X > 0).sum(axis=1)).ravel()
    a.obs["n_genes_by_counts"] = ngenes
    a.obs["log1p_n_genes_by_counts"] = np.log1p(ngenes)
    a.obs["total_counts"] = tot
    a.obs["log1p_total_counts"] = np.log1p(tot)
    Xd = dense(X).astype(float) if X.shape[1] <= 5000 else None
    if Xd is not None:
        srt = -np.sort(-Xd, axis=1)
        for n in (50, 100, 200, 500):
            if n <= Xd.shape[1]:
                a.obs[f"pct_counts_in_top_{n}_genes"] = 100 * srt[:, :n].sum(axis=1) / np.maximum(tot, 1e-300)
    for q in qc_vars:
        m = np.asarray(a.var[q], bool)
        t = np.asarray(X[:, np.where(m)[0]].sum(axis=1)).ravel().astype(float)
        a.obs[f"total_counts_{q}"] = t
        a.obs[f"log1p_total_counts_{q}"] = np.log1p(t)
        a.obs[f"pct_counts_{q}"] = 100 * t / np.maximum(tot, 1e-300)
    ncell = np.asarray((X > 0).sum(axis=0)).ravel()
    gtot = np.asarray(X.sum(axis=0)).ravel().astype(float)
    a.var["n_cells_by_counts"] = ncell
    a.var["mean_counts"] = gtot / a.n_obs
    a.var["log1p_mean_counts"] = np.log1p(gtot / a.n_obs)
    a.var["pct_dropout_by_counts"] = 100 * (1 - ncell / a.n_obs)
    a.var["total_counts"] = gtot
    a.var["log1p_total_counts"] = np.log1p(gtot)


def preprocess_rna(a, gene_names=None, min_genes: int = 500, min_cells: int = 3, pct_mito: float = 20.0,
                   reference: str = "human", compat: bool = False):
    np = _np()
    a = a.copy()
    if gene_names is not None:
        if len(gene_names) != a.n_vars:
            raise SystemExit("The number of gene names does not match the number of variables in adata_rna")
        a.var["symbol"] = list(gene_names)
    elif "symbol" not in a.var.columns:
        a.var["symbol"] = list(a.var_names)
    a.var_names = strip_version(a.var_names)
    a.var_names_make_unique()
    if "batch" not in a.obs.columns:
        a.obs["batch"] = "batch_1"
    a.obs["batch_number"] = __import__("pandas").factorize(a.obs["batch"])[0] + 1
    prefix = "MT-" if (compat or reference == "human") else "Mt-"
    a.var["mt"] = a.var["symbol"].astype(str).str.startswith(prefix)
    a.var["ribo"] = a.var["symbol"].astype(str).str.startswith(("RPS", "RPL"))
    calculate_qc_metrics(a, ("mt", "ribo"))
    ngenes = a.obs["n_genes_by_counts"].values
    f = a[ngenes >= min_genes].copy()
    f.obs["n_genes"] = f.obs["n_genes_by_counts"].values
    ncell = np.asarray((f.X > 0).sum(axis=0)).ravel()
    f = f[:, ncell >= min_cells].copy()
    f.var["n_cells"] = np.asarray((f.X > 0).sum(axis=0)).ravel()
    f = f[f.obs["pct_counts_mt"].values < pct_mito].copy()
    return a, f


def knee_plot(a, path: Path, title: str, index_on_y: bool = True):
    np = _np()
    plt = _plt()
    if plt is None:
        return None
    s = np.sort(np.asarray(a.X.sum(axis=1)).ravel())[::-1]
    fig, ax = plt.subplots(figsize=(8, 5))
    if index_on_y:
        ax.plot(np.log1p(s), np.arange(len(s)), marker="o", ms=2, lw=1, color=BLUE)
        _style(ax, title, "Log of UMI Counts", "Barcode number")
    else:
        ax.plot(np.arange(len(s)), np.log1p(s), marker="o", ms=2, lw=1, color=BLUE)
        _style(ax, title, "Barcode Index", "Log of UMI Counts")
    return _save(fig, path)


def cmd_preprocess(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    a = read_ad(Path(args.adata))
    names = None
    if args.gene_names:
        names = [ln.strip() for ln in open(args.gene_names) if ln.strip()]
    full, f = preprocess_rna(a, names, args.min_genes, args.min_cells, args.pct_mito, args.reference, args.upstream_compat)
    write_ad(full, d / "rna_qc_unfiltered.h5ad")
    out = Path(args.out) if args.out else d / "filtered_anndata.h5ad"
    write_ad(f, out)
    figs = []
    if not args.no_plots:
        (d / "figures").mkdir(exist_ok=True)
        figs.append(knee_plot(full, d / "figures" / "knee_plot_scRNA.png", "Knee Plot"))
    write_report(d, "CRISPR Jamboree 3 — scRNA preprocessing",
                 [f"{full.n_obs} barcodes -> {f.n_obs} cells (n_genes >= {args.min_genes}, pct_counts_mt < {args.pct_mito}); "
                  f"{full.n_vars} -> {f.n_vars} genes (>= {args.min_cells} cells)."],
                 {"n_barcodes": full.n_obs, "n_cells": f.n_obs, "n_genes": f.n_vars, "out": out, "figures": [x for x in figs if x]})
    return 0


# ---------------------------------------------------------------------------
# GTF + create-mudata (create_mdata.py, create_mdata_HASHING.py)
# ---------------------------------------------------------------------------

def read_gtf_genes(path: Path):
    """gene records of a GENCODE GTF: chr, start, end, strand, gene_id, gene_name, gene_type."""
    pd = _pd()
    rows = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            rows.append({"chr": f[0], "start": int(f[3]), "end": int(f[4]), "strand": f[6], "gene_id": attrs.get("gene_id", ""),
                         "gene_name": attrs.get("gene_name", attrs.get("gene_id", "")), "gene_type": attrs.get("gene_type", "")})
    return pd.DataFrame(rows)


def _is_nt_row(r) -> bool:
    t = str(r.get("Target_name", "")).lower()
    return t.startswith("non-targeting") or t.startswith("non_targeting") or str(r.get("chr", "")).lower() in {"non_targeting", "non-targeting"} \
        or t.startswith("safe")


def create_mudata(rna, guide, meta, gtf_df, moi: str = "high", hashing=None, compat: bool = False) -> MuDataLite:
    np, pd = _np(), _pd()
    guide = guide.copy()
    rna = rna.copy()
    meta = meta.copy()
    meta["sgRNA_ID"] = meta["sgRNA_ID"].astype(str)
    gid = list(map(str, guide.var["gene_id"])) if "gene_id" in guide.var.columns else list(map(str, guide.var_names))
    if len(meta) != len(gid):
        print(f"The numbers of sgRNA_ID/guide_id are different: There are {len(gid)} in guide anndata, but there are {len(meta)} in guide metadata.")
    keep = [g in set(meta["sgRNA_ID"]) for g in gid]
    guide = guide[:, np.where(keep)[0]].copy()
    gid = [g for g, k in zip(gid, keep) if k]
    m = meta.drop_duplicates("sgRNA_ID").set_index("sgRNA_ID")
    src = meta.reset_index(drop=True).iloc[: len(gid)] if compat else m.loc[gid].reset_index()  # upstream assigns positionally
    var = pd.DataFrame({"guide_id": gid}, index=range(len(gid)))
    var["intended_target_name"] = src["Target_name"].astype(str).values
    var["intended_target_chr"] = src["chr"].astype(str).values
    var["intended_target_start"] = pd.to_numeric(src["start"], errors="coerce").values
    var["intended_target_end"] = pd.to_numeric(src["end"], errors="coerce").values
    var["sequence"] = src["sgRNA_sequences"].astype(str).values
    var["targeting"] = "TRUE" if compat else ["FALSE" if _is_nt_row(r) else "TRUE" for _, r in src.iterrows()]
    var.index = [f"{a}|{b}" for a, b in zip(src["sgRNA_ID"], src["sgRNA_sequences"])]
    guide.var = var
    guide.uns["moi"] = np.array([moi], dtype=object)
    guide.uns["capture_method"] = np.array(["CROP-seq"], dtype=object)
    guide.obs["num_expressed_guides"] = np.asarray((guide.X > 0).sum(axis=1)).ravel()
    if "batch" not in guide.obs.columns:
        guide.obs["batch"] = "batch_1"
    guide.obs["batch_number"] = pd.factorize(guide.obs["batch"])[0] + 1
    guide.obs["total_guide_umis"] = np.asarray(guide.X.sum(axis=1)).ravel().astype(float)
    if gtf_df is not None and len(gtf_df):
        g = gtf_df.copy()
        g["gene_id2"] = strip_version(g["gene_id"])
        g = g.drop_duplicates("gene_id2").set_index("gene_id2")
        for c_out, c_in in (("gene_chr", "chr"), ("gene_start", "start"), ("gene_end", "end")):
            rna.var[c_out] = g[c_in].reindex(strip_version(rna.var_names)).values
    rna.obs = rna.obs.rename(columns={"n_genes_by_counts": "n_counts", "pct_counts_mt": "percent_mito", "n_genes": "num_expressed_genes",
                                      "total_counts": "total_gene_umis"})
    if "total_gene_umis" not in rna.obs.columns:
        rna.obs["total_gene_umis"] = np.asarray(rna.X.sum(axis=1)).ravel().astype(float)
    if "num_expressed_genes" not in rna.obs.columns:
        rna.obs["num_expressed_genes"] = np.asarray((rna.X > 0).sum(axis=1)).ravel()
    common = set(rna.obs_names) & set(guide.obs_names)
    if hashing is not None:
        common &= set(hashing.obs_names)
    bc = sorted(common)
    mods = {"gene": rna[bc].copy(), "guide": guide[bc].copy()}
    if hashing is not None:
        mods["hashing"] = hashing[bc].copy()
    shared = set(mods["guide"].obs.columns) & set(mods["gene"].obs.columns)
    if hashing is not None:
        shared &= set(mods["hashing"].obs.columns)
    obs = mods["guide"].obs.loc[:, sorted(shared)].copy()
    return MuDataLite(mods, obs, {})


def cmd_create_mudata(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    rna = read_ad(Path(args.rna))
    guide = read_ad(Path(args.guide))
    hashing = read_ad(Path(args.hashing)) if args.hashing else None
    meta = read_table(Path(args.guide_metadata))
    gtf = read_gtf_genes(Path(args.gtf)) if args.gtf else None
    md = create_mudata(rna, guide, meta, gtf, args.moi, hashing, args.upstream_compat)
    out = Path(args.out) if args.out else d / "mudata.h5mu"
    write_h5mu(md, out)
    figs = []
    if not args.no_plots:
        (d / "figures").mkdir(exist_ok=True)
        figs.append(knee_plot(guide, d / "figures" / "knee_plot_guide.png", "Knee Plot", index_on_y=False))
    write_report(d, "CRISPR Jamboree 3 — MuData", [md_table(["modality", "cells", "features"], [[k, a.n_obs, a.n_vars] for k, a in md.mods.items()]),
                                                  f"Shared obs columns: {', '.join(md.obs.columns)}."],
                 {"modalities": {k: [a.n_obs, a.n_vars] for k, a in md.mods.items()}, "out": out, "figures": [f for f in figs if f]})
    return 0


# ---------------------------------------------------------------------------
# doublets (doublets.py: Scrublet)
# ---------------------------------------------------------------------------

def threshold_minimum(values, nbins: int = 256, max_num_iter: int = 10000) -> float:
    """skimage.filters.threshold_minimum on a 1-D sample."""
    np = _np()
    from scipy import ndimage  # type: ignore
    hist, edges = np.histogram(values, bins=nbins)
    centers = (edges[:-1] + edges[1:]) / 2

    def maxima(h):
        idx, direction = [], 1
        for i in range(h.shape[0] - 1):
            if direction > 0:
                if h[i + 1] < h[i]:
                    direction = -1
                    idx.append(i)
            elif h[i + 1] > h[i]:
                direction = 1
        return idx
    sm = hist.astype(np.float64)
    mx: list = []
    for counter in range(max_num_iter):
        sm = ndimage.uniform_filter1d(sm, 3)
        mx = maxima(sm)
        if len(mx) < 3:
            break
    if len(mx) != 2:
        raise RuntimeError("Unable to find two maxima in histogram")
    t = int(np.argmin(sm[mx[0]: mx[1] + 1]))
    return float(centers[mx[0] + t])


def scrublet(X, sim_doublet_ratio: float = 2.0, n_neighbors: Optional[int] = None, expected_doublet_rate: float = 0.1,
             min_counts: int = 3, min_cells: int = 3, min_gene_variability_pctl: float = 85, n_prin_comps: int = 30, seed: int = 0) -> dict:
    np = _np()
    from scipy.spatial import cKDTree  # type: ignore
    rng = np.random.default_rng(seed)
    E = dense(X).astype(float)
    n = E.shape[0]
    k = n_neighbors or int(round(0.5 * math.sqrt(n)))
    tot = E.sum(axis=1)
    norm = E / np.maximum(tot, 1)[:, None] * tot.mean()
    keep = (E >= min_counts).sum(axis=0) >= min_cells
    mu = norm.mean(axis=0)
    var = norm.var(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ff = np.where(mu > 0, var / mu, 0)
        pos = keep & (mu > 0)
        c = float(np.median(np.clip((ff[pos] - 1) / mu[pos], 0, None))) if pos.any() else 0.0
        vscore = np.where(mu > 0, ff / (1 + c * mu), 0)
    cut = np.percentile(vscore[keep], min_gene_variability_pctl) if keep.any() else 0
    genes = np.where(keep & (vscore >= cut))[0]
    if genes.size < 2:
        genes = np.where(keep)[0]
    Eg = E[:, genes]
    n_sim = int(round(n * sim_doublet_ratio))
    pairs = rng.integers(0, n, size=(n_sim, 2))
    Es = Eg[pairs[:, 0]] + Eg[pairs[:, 1]]
    tot_s = tot[pairs[:, 0]] + tot[pairs[:, 1]]
    On = Eg / np.maximum(tot, 1)[:, None] * 1e5
    Sn = Es / np.maximum(tot_s, 1)[:, None] * 1e5
    m, s = On.mean(axis=0), On.std(axis=0)
    s[s == 0] = 1
    On, Sn = (On - m) / s, (Sn - m) / s
    npc = min(n_prin_comps, On.shape[1], n - 1)
    U, S, Vt = np.linalg.svd(On - On.mean(axis=0), full_matrices=False)
    P = Vt[:npc].T
    Po, Ps = (On - On.mean(axis=0)) @ P, (Sn - On.mean(axis=0)) @ P
    allp = np.vstack([Po, Ps])
    r = n_sim / n
    k_adj = int(round(k * (1 + r)))
    _, nb = cKDTree(allp).query(allp, k=k_adj + 1)
    nb = nb[:, 1:]
    nd = (nb >= n).sum(axis=1)
    rho = expected_doublet_rate
    q = (nd + 1) / (k_adj + 2)
    Ld = q * rho / r / (1 - rho - q * (1 - rho - rho / r))
    obs, sim = Ld[:n], Ld[n:]
    try:
        thr = threshold_minimum(sim)
    except RuntimeError:
        thr = None
    return {"doublet_scores": obs, "sim_scores": sim, "threshold": thr,
            "predicted_doublets": (obs > thr) if thr is not None else np.zeros(n, bool), "n_genes_used": int(genes.size), "k_adj": k_adj}


def remove_doublets(md: MuDataLite, seed: int = 0) -> "tuple[MuDataLite, dict]":
    sc = scrublet(md.gene.X, seed=seed)
    g = md.gene.copy()
    g.obs["doublet_scores"] = sc["doublet_scores"]
    g.obs["predicted_doublets"] = sc["predicted_doublets"]
    g.obs["doublet_info"] = [str(bool(x)) for x in sc["predicted_doublets"]]
    kept = g.obs_names[~sc["predicted_doublets"]]
    mods = dict(md.mods)
    mods["gene"] = g
    common = set(kept)
    for k in mods:
        common &= set(mods[k].obs_names)
    bc = [b for b in g.obs_names if b in common]
    out = MuDataLite({k: v[bc].copy() for k, v in mods.items()}, None, dict(md.uns))
    shared = set(out.guide.obs.columns) & set(out.gene.obs.columns)
    if "hashing" in out:
        shared &= set(out["hashing"].obs.columns)
    out.obs = out.guide.obs.loc[:, sorted(shared)].copy()
    return out, sc


def cmd_doublets(args: argparse.Namespace) -> int:
    np = _np()
    d = run_dir(args.label)
    md = read_h5mu(Path(args.mudata))
    out_md, sc = remove_doublets(md, args.seed)
    out = Path(args.out) if args.out else d / "mdata_doublets.h5mu"
    write_h5mu(out_md, out)
    figs = []
    plt = None if args.no_plots else _plt()
    if plt is not None:
        fig, axes = plt.subplots(1, 2, figsize=(8, 3))
        for ax, v, t in ((axes[0], sc["doublet_scores"], "Observed transcriptomes"), (axes[1], sc["sim_scores"], "Simulated doublets")):
            ax.hist(v, bins=np.linspace(0, 1, 50), color=BLUE)
            if sc["threshold"] is not None:
                ax.axvline(sc["threshold"], color=RED, lw=1)
            _style(ax, t, "Doublet score", "Prob. density")
        (d / "figures").mkdir(exist_ok=True)
        figs.append(_save(fig, d / "figures" / "doublets_batch.png"))
    n_d = int(sc["predicted_doublets"].sum())
    print("Number of predicted doublets:", n_d)
    write_report(d, "CRISPR Jamboree 3 — Scrublet doublets",
                 [f"Threshold {_fmt(sc['threshold'])}; {n_d} predicted doublets removed; {out_md.gene.n_obs} cells kept."],
                 {"threshold": sc["threshold"], "n_doublets": n_d, "n_cells_kept": out_md.gene.n_obs, "out": out, "figures": figs})
    return 0


# ---------------------------------------------------------------------------
# demultiplex (GMM-demux, demultiplex_prepare.py, demultiplex_filter.py, filter_hashing.py)
# ---------------------------------------------------------------------------

def gmm_positive(x, iters: int = 200):
    """Two-component 1-D Gaussian mixture by EM; returns posterior of the high component."""
    np = _np()
    from scipy.stats import norm  # type: ignore
    lo, hi = np.percentile(x, 25), np.percentile(x, 95)
    mu = np.array([lo, hi if hi > lo else lo + 1])
    sd = np.array([max(x.std(), 1e-3)] * 2)
    w = np.array([0.5, 0.5])
    for _ in range(iters):
        l0 = w[0] * norm.pdf(x, mu[0], sd[0])
        l1 = w[1] * norm.pdf(x, mu[1], sd[1])
        p1 = l1 / np.maximum(l0 + l1, 1e-300)
        p0 = 1 - p1
        w_new = np.array([p0.mean(), p1.mean()])
        mu_new = np.array([np.sum(p0 * x) / max(p0.sum(), 1e-12), np.sum(p1 * x) / max(p1.sum(), 1e-12)])
        sd = np.sqrt(np.array([np.sum(p0 * (x - mu_new[0]) ** 2) / max(p0.sum(), 1e-12), np.sum(p1 * (x - mu_new[1]) ** 2) / max(p1.sum(), 1e-12)]))
        sd = np.maximum(sd, 1e-3)
        done = np.max(np.abs(mu_new - mu)) < 1e-8
        mu, w = mu_new, w_new
        if done:
            break
    if mu[1] < mu[0]:
        p1 = 1 - p1
    return p1


def gmm_demux(H, hto_names: "list[str]", threshold: float = 0.8):
    """Cluster_id / hto_type / Confidence per cell; config = cluster_id -> name."""
    np, pd = _np(), _pd()
    L = np.log1p(dense(H).astype(float))
    post = np.column_stack([gmm_positive(L[:, j]) for j in range(L.shape[1])])
    posb = post > 0.5
    conf = np.prod(np.where(posb, post, 1 - post), axis=1)
    cid = (posb * (2 ** np.arange(L.shape[1]))[None, :]).sum(axis=1)
    names = ["-".join(h for h, b in zip(hto_names, row) if b) or "negative" for row in posb]
    config = sorted(set(zip(cid.tolist(), names)))
    return pd.DataFrame({"Cluster_id": cid, "Confidence": conf, "hto_type": names}), config


def demultiplex_hashing(hashing, rna_barcodes=None):
    np, pd = _np(), _pd()
    h = hashing.copy()
    if rna_barcodes is not None:
        keep = [b for b in h.obs_names if b in set(rna_barcodes)]
        h = h[keep].copy()
    if "batch" not in h.obs.columns:
        h.obs["batch"] = "batch_1"
    parts, reports = [], {}
    for b in h.obs["batch"].astype(str).unique():
        sub = h[h.obs["batch"].astype(str).values == b].copy()
        rep, config = gmm_demux(sub.X, list(map(str, sub.var_names)))
        rep.index = sub.obs_names
        cfg = pd.DataFrame(config, columns=["cluster_id", "hto_type"])
        cfg["hto_type"] = cfg["hto_type"].str.strip()
        cfg["hto_type_split"] = cfg["hto_type"].str.split("-").str.join(",").apply(lambda x: "multiplets" if "," in x else x)
        sub.obs["cluster_id"] = rep["Cluster_id"].values
        sub.obs = sub.obs.merge(cfg, how="left", on="cluster_id").set_index(sub.obs.index)
        reports[b] = (rep, cfg)
        parts.append(sub[sub.obs["hto_type"].values != "negative"].copy())
    import anndata as ad  # type: ignore
    return ad.concat(parts, join="outer") if len(parts) > 1 else parts[0], reports


def cmd_demultiplex(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    h = read_ad(Path(args.hashing))
    rna_bc = list(read_ad(Path(args.rna)).obs_names) if args.rna else None
    demux, reports = demultiplex_hashing(h, rna_bc)
    for b, (rep, cfg) in reports.items():
        rep.to_csv(d / f"{b}_GMM_full.csv")
        cfg[["cluster_id", "hto_type"]].to_csv(d / f"{b}_GMM_full.config", header=False, index=False)
        print(f"CSV: {d / f'{b}_GMM_full.csv'}")
    out = Path(args.out) if args.out else d / "concatenated_hashing_demux.h5ad"
    write_ad(demux, out)
    counts = demux.obs["hto_type_split"].value_counts().to_dict()
    write_report(d, "CRISPR Jamboree 3 — hashing demultiplexing",
                 [f"{h.n_obs} hashing barcodes{' (intersected with RNA)' if rna_bc else ''} -> {demux.n_obs} non-negative cells.",
                  md_table(["hto_type_split", "cells"], sorted(counts.items()))],
                 {"n_input": h.n_obs, "n_kept": demux.n_obs, "hto_type_split": counts, "out": out})
    return 0


# ---------------------------------------------------------------------------
# assign-guides (cleanser.py, assign_grnas_sceptre.R, add_guide_assignment.py)
# ---------------------------------------------------------------------------

def cmd_assign(args: argparse.Namespace) -> int:
    np, pd = _np(), _pd()
    d = run_dir(args.label)
    md = read_h5mu(Path(args.mudata))
    if args.method == "cleanser":
        A, params = assign_cleanser(md.guide, args.threshold, None, pooled=args.upstream_compat)
    elif args.method == "umi-threshold":
        A, params = assign_threshold(md.guide, args.threshold), {}
    else:
        A, params = assign_sceptre_mixture(md.guide, md.gene), {}
    md.guide.layers["guide_assignment"] = A
    out = Path(args.out) if args.out else d / f"{args.method}_assignment_mudata.h5mu"
    write_h5mu(md, out)
    Ad = dense(A)
    per_cell = (Ad > 0).sum(axis=1)
    pd.DataFrame(Ad, index=md.guide.obs_names, columns=md.guide.var_names).T.reset_index(drop=True).to_csv(d / "guide_assignment.csv", index=False) \
        if Ad.size <= 5_000_000 else None
    if params:
        write_json(params, d / "cleanser_params.json")
    write_report(d, "CRISPR Jamboree 3 — guide assignment",
                 [f"`{args.method}` (threshold {args.threshold}{', pooled fit as upstream' if args.upstream_compat and args.method == 'cleanser' else ''}): "
                  f"{int((per_cell > 0).sum())} of {len(per_cell)} cells with >= 1 gRNA, mean {per_cell.mean():.2f}."],
                 {"method": args.method, "cells_with_guide": int((per_cell > 0).sum()), "mean_guides_per_cell": float(per_cell.mean()), "out": out})
    return 0


# ---------------------------------------------------------------------------
# prepare-inference (create_pairs_to_test.py, prepare_inference.py)
# ---------------------------------------------------------------------------

def pairs_from_gtf(guide_var, gtf_df, limit: int = 1_000_000, compat: bool = False):
    np, pd = _np(), _pd()
    g = gtf_df.drop_duplicates("gene_name").copy()
    g["gene_id2"] = strip_version(g["gene_id"])
    rows = []
    for _, r in guide_var.iterrows():
        if limit < 0:
            cand = g
        else:
            cand = g[g["chr"] == str(r["intended_target_chr"])]
            cand = cand[np.abs(cand["start"] - float(r["intended_target_start"])) <= limit] if pd.notna(r["intended_target_start"]) else cand.iloc[:0]
        for _, c in cand.iterrows():
            rows.append({"guide_id": r["guide_id"], "gene_name": c["gene_name"] if compat else c["gene_id2"],
                         "intended_target_name": r["intended_target_name"], "pair_type": "discovery"})
    return pd.DataFrame(rows, columns=["guide_id", "gene_name", "intended_target_name", "pair_type"])


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


def cmd_prepare_inference(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    md = read_h5mu(Path(args.mudata))
    if args.pairs:
        pairs = read_table(Path(args.pairs))
    else:
        if not args.gtf:
            raise SystemExit("--pairs (user pairs CSV) or --gtf (distance-based pairs) is required")
        pairs = pairs_from_gtf(md.guide.var, read_gtf_genes(Path(args.gtf)), args.distance_from_center, args.upstream_compat)
        write_tsv(pairs, d / "pairs_to_test.csv")
    sub = restrict_pairs(pairs, md)
    md.uns["pairs_to_test"] = sub
    out = Path(args.out) if args.out else d / "mudata_inference_input.h5mu"
    write_h5mu(md, out)
    write_report(d, "CRISPR Jamboree 3 — pairs to test", [f"{len(pairs)} candidate pairs -> {len(sub)} in the MuData.",
                                                         md_table(["pair_type", "pairs"], sub["pair_type"].astype(str).value_counts().items()
                                                                  if "pair_type" in sub.columns else [])],
                 {"n_candidates": len(pairs), "n_pairs": len(sub), "out": out})
    return 0


# ---------------------------------------------------------------------------
# infer (inference_sceptre.R, perturbo_inference.py, add_guide_inference.py)
# ---------------------------------------------------------------------------

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


def infer_perturbo_j3(md_path: Path, out_path: Path):
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


def cmd_infer(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    md = read_h5mu(Path(args.mudata))
    cov = args.cov_string or (Path(args.cov_string_file).read_text().strip() if args.cov_string_file else "")
    if args.method == "sceptre":
        tr, used = sceptre_inference(md, args.side, args.formula, cov, args.resamples, args.seed, args.moi)
    else:
        tr, used = infer_perturbo_j3(Path(args.mudata), d / "inference_mudata.h5mu"), ""
        if tr is None:
            write_report(d, "CRISPR Jamboree 3 — inference", ["PerTurbo not installed; nothing run."], {"method": "perturbo", "ran": False})
            return 0
    tr.to_csv(d / "test_results.csv", index=False)
    print(f"CSV: {d / 'test_results.csv'}")
    md.uns["test_results"] = tr
    out = Path(args.out) if args.out else d / "inference_mudata.h5mu"
    write_h5mu(md, out)
    n5 = int((tr["p_value"] < 0.05).sum())
    write_report(d, f"CRISPR Jamboree 3 — {args.method} inference",
                 [f"{len(tr)} pairs; {n5} with p < 0.05; covariates `{used or 'none (formula default)'}`.",
                  md_table(list(tr.columns), tr.sort_values("p_value").head(15).values.tolist())],
                 {"method": args.method, "n_pairs": len(tr), "n_p_lt_0.05": n5, "covariates": used, "out": out})
    return 0


# ---------------------------------------------------------------------------
# evaluate (volcano_plot.py, select_nodes.py, network_plot.py, igv.py)
# ---------------------------------------------------------------------------

def volcano_sets(tr, log2fc: float = 1.0, pval: float = 0.05):
    np = _np()
    down = tr[(tr["log2_fc"] <= -log2fc) & (tr["p_value"] <= pval)]
    up = tr[(tr["log2_fc"] >= log2fc) & (tr["p_value"] <= pval)]
    dc = down.copy()
    dc["log2_fc"] = dc["log2_fc"].replace([np.inf, -np.inf], [1e10, -1e10])
    dc["p_value"] = dc["p_value"].replace([np.inf, -np.inf], [1e10, -1e10])
    annotated = __import__("pandas").concat([dc.sort_values("log2_fc").head(10), dc.sort_values("p_value").head(10)])
    return down, up, annotated


def select_nodes(tr, num_nodes: int) -> "list[str]":
    return list(tr.sort_values("p_value")["intended_target_name"].unique()[:num_nodes])


def igv_files(md: MuDataLite, tr, gtf_df):
    """igv.py: gene coordinates from gene.var, element coordinates from guide.var, promoter when target == GTF gene_name."""
    np, pd = _np(), _pd()
    coord = {}
    gv = md.gene.var
    for idx, r in gv.iterrows():
        try:
            if np.isnan(float(r.get("gene_start", np.nan))) or np.isnan(float(r.get("gene_end", np.nan))):
                continue
        except (TypeError, ValueError):
            continue
        coord[idx] = [r["gene_chr"], r["gene_start"], r["gene_end"]]
    for _, r in md.guide.var.iterrows():
        t = r["intended_target_name"]
        if t in coord or t == "non-targeting":
            continue
        coord[t] = [r["intended_target_chr"], r["intended_target_start"], r["intended_target_end"]]
    name = {}
    if gtf_df is not None:
        name = dict(zip(strip_version(gtf_df["gene_id"]), gtf_df["gene_name"]))
    elif "symbol" in gv.columns:
        name = dict(zip(gv.index, gv["symbol"]))
    bp, bg = defaultdict(list), defaultdict(list)
    for _, r in tr.iterrows():
        t, g = r["intended_target_name"], r["gene_id"]
        if t not in coord:
            continue
        if t == name.get(g):
            for k, v in zip(("chr", "start", "end"), coord[t]):
                bg[k].append(v)
            bg["p_value"].append(r["p_value"]); bg["log2_fc"].append(r["log2_fc"])
        elif g in coord:
            for k, v in zip(("chr1", "start1", "end1"), coord[t]):
                bp[k].append(v)
            for k, v in zip(("chr2", "start2", "end2"), coord[g]):
                bp[k].append(v)
            bp["p_value"].append(r["p_value"]); bp["log2_fc"].append(r["log2_fc"])
    return pd.DataFrame(bp), pd.DataFrame(bg)


def evaluate_inference(md: MuDataLite, tr, gtf_df, outdir: Path, central_nodes=None, num_nodes: int = 2, log2fc: float = 1.0,
                       pval: float = 0.05, min_weight: float = 0.1, no_plots: bool = False) -> dict:
    np, pd = _np(), _pd()
    ev = outdir / "evaluation_output"
    ev.mkdir(parents=True, exist_ok=True)
    tr = tr.copy()
    for c in ("p_value", "log2_fc"):
        tr[c] = pd.to_numeric(tr[c], errors="coerce")
    down, up, ann = volcano_sets(tr, log2fc, pval)
    write_tsv(down, ev / "volcano_down.tsv")
    write_tsv(up, ev / "volcano_up.tsv")
    nodes = central_nodes if central_nodes and central_nodes != ["undefined"] else select_nodes(tr, num_nodes)
    edges = []
    for c in nodes:
        e = network_edges(tr, c, weight="log2_fc", min_weight=min_weight)
        if e.empty:
            print(f"Error: Central node '{c}' not found in the results(filtered by min_weight). Skipping plot.")
            continue
        edges.append(e.assign(central_node=c))
    edges = pd.concat(edges, ignore_index=True) if edges else pd.DataFrame()
    write_tsv(edges, ev / "network_edges.tsv")
    bp, bg = igv_files(md, tr, gtf_df)
    bp.to_csv(ev / "output.bedpe", sep="\t", index=False, header=False)
    bg.to_csv(ev / "output.bedgraph", sep="\t", index=False, header=False)
    print(f"Wrote: {ev / 'output.bedpe'}")
    print(f"Wrote: {ev / 'output.bedgraph'}")
    figs = []
    plt = None if no_plots else _plt()
    if plt is not None:
        fig, ax = plt.subplots(figsize=(6, 4.5))
        ax.scatter(tr["log2_fc"], -np.log10(np.clip(tr["p_value"], 1e-300, 1)), s=5, color=GREEN, label="Not significant")
        ax.scatter(down["log2_fc"], -np.log10(np.clip(down["p_value"], 1e-300, 1)), s=10, color=BLUE, label="Down-regulated")
        ax.scatter(up["log2_fc"], -np.log10(np.clip(up["p_value"], 1e-300, 1)), s=10, color=RED, label="Up-regulated")
        for _, r in ann.drop_duplicates().iterrows():
            ax.text(r["log2_fc"], -math.log10(max(float(r["p_value"]), 1e-300)), f"{r['intended_target_name']} ({r.get('guide_id', '')})", fontsize=6)
        for v in (-log2fc, log2fc):
            ax.axvline(v, color=GREY, ls="--", lw=0.8)
        ax.axhline(-math.log10(pval), color=GREY, ls="--", lw=0.8)
        ax.legend(loc="upper right", frameon=False, fontsize=7)
        _style(ax, "Volcano", "log2 fold change", "p-value (log10)")
        figs.append(_save(fig, ev / "volcano_plot.png"))
        valid = [c for c in nodes if len(edges) and (edges["central_node"] == c).any()]
        if valid:
            rows_ = (len(valid) + 1) // 2
            fig, axes = plt.subplots(rows_, 2, figsize=(15, 7.5 * rows_), squeeze=False)
            for i, c in enumerate(valid):
                plot_network(edges[edges["central_node"] == c], c, "intended_target_name", "gene_id", "log2_fc", "p_value", axes[i // 2, i % 2])
            for j in range(len(valid), rows_ * 2):
                axes[j // 2, j % 2].set_axis_off()
            figs.append(_save(fig, ev / "network_plot.png"))
    print("Number of down-regulated genes based on the thresholds:", len(down))
    print("Number of up-regulated genes based on the thresholds::", len(up))
    return {"n_down": len(down), "n_up": len(up), "central_nodes": nodes, "n_edges": len(edges), "n_bedpe": len(bp), "n_bedgraph": len(bg),
            "annotated": ann.drop_duplicates()[["intended_target_name", "gene_id"]].astype(str).values.tolist(), "figures": figs}


def cmd_evaluate(args: argparse.Namespace) -> int:
    d = run_dir(args.label)
    md = read_h5mu(Path(args.mudata))
    tr = uns_df(md, "test_results")
    gtf = read_gtf_genes(Path(args.gtf)) if args.gtf else None
    ev = evaluate_inference(md, tr, gtf, d, args.central_nodes, args.num_nodes, args.log2_fc, args.p_value, args.min_weight, args.no_plots)
    write_report(d, "CRISPR Jamboree 3 — evaluation",
                 [f"Volcano (|log2_fc| >= {args.log2_fc}, p <= {args.p_value}): {ev['n_down']} down, {ev['n_up']} up.",
                  f"Network central nodes {', '.join(map(str, ev['central_nodes']))}: {ev['n_edges']} edges; IGV {ev['n_bedgraph']} bedgraph, {ev['n_bedpe']} bedpe."],
                 ev)
    return 0


# ---------------------------------------------------------------------------
# dashboard (create_dashboard_plots.py, create_dashboard_df.py, process_json.py, create_dashboard.py)
# ---------------------------------------------------------------------------

def human_format(num) -> str:
    num = float("{:.3g}".format(float(num)))
    mag = 0
    while abs(num) >= 1000:
        mag += 1
        num /= 1000.0
    return "{}{}".format("{:f}".format(num).rstrip("0").rstrip("."), ["", "K", "M", "B", "T"][mag])


def threshold_counts(md: MuDataLite) -> dict:
    np = _np()
    gx = md.gene.X
    tot = np.asarray(gx.sum(axis=1)).ravel()
    ng = np.asarray((gx > 0).sum(axis=1)).ravel()
    gt = np.asarray(md.guide.X.sum(axis=1)).ravel()
    return {"scRNA_barcodes_UMI_thresholds": {t: int((tot > t).sum()) for t in (200, 500, 1000, 2000, 5000)},
            "scRNA_barcodes_detected_genes_thresholds": {t: int((ng > t).sum()) for t in (200, 500, 1000, 1500, 2000)},
            "guides_UMI_thresholds": {t: int((gt > t).sum()) for t in (0, 5, 10, 50, 100, 1000)}}


def kb_json_tables(dirs: "list[Path]") -> list:
    """process_json.py + create_json_df: inspect.json + run_info.json of each <batch>_ks_<modality>_out."""
    pd = _pd()
    out = []
    for p in dirs or []:
        p = Path(p)
        m = re.match(r"(.+)_ks_(transcripts|guide|hashing)_out", p.name)
        if not m or not (p / "inspect.json").is_file() or not (p / "run_info.json").is_file():
            continue
        data = {**json.loads((p / "inspect.json").read_text()), **json.loads((p / "run_info.json").read_text())}
        data = {k: v for k, v in data.items() if k not in ("start_time", "call")}
        disp = dict(data)
        for k in ("n_targets", "n_processed", "n_unique", "n_pseudoaligned", "numRecords", "numReads", "numBarcodes", "numUMIs",
                  "numBarcodeUMIs", "gtRecords", "numBarcodesOnOnlist", "numReadsOnOnlist"):
            if k in disp:
                disp[k] = human_format(disp[k])
        for k in ("meanReadsPerBarcode", "meanUMIsPerBarcode"):
            if k in disp:
                disp[k] = f"{float(disp[k]):.1f}"
        for k in ("p_pseudoaligned", "p_unique", "percentageBarcodesOnOnlist", "percentageReadsOnOnlist"):
            if k in disp:
                disp[k] = f"{float(disp[k]):.1f}%"
        mod = {"transcripts": "scRNA", "guide": "Guide", "hashing": "Hashing"}[m.group(2)]
        value = (f"Total Reads for {mod}: {disp.get('n_processed')}, Paired Reads Mapped: {disp.get('n_pseudoaligned')}, "
                 f"Alignment Percentage: {disp.get('p_pseudoaligned')}, Total Detected {mod} Barcodes (Unfiltered): {disp.get('numBarcodes')}")
        out.append({"modality": mod, "description": m.group(1), "subject": f"Mapping {mod}", "value_display": value,
                    "table": pd.DataFrame({"parameter": list(disp), "value": list(disp.values())})})
    return out


def dashboard_blocks(md: MuDataLite, gene_unfiltered=None, gene_filtered=None, guide_unfiltered=None, position_table=None) -> "list[dict]":
    np, pd = _np(), _pd()
    blocks = []
    inter = len(set(gene_unfiltered.obs_names) & set(guide_unfiltered.obs_names)) if gene_unfiltered is not None and guide_unfiltered is not None else None
    if position_table is not None:
        blocks.append({"modality": "Guide", "subject": "Fastq Overview", "value_display": "", "table": position_table,
                       "table_description": "Summary of Sequence Index: positions where the guide starts on the reads"})
    cn = []
    if inter is not None:
        cn.append(f"Number of guide barcodes (unfiltered) intersecting with scRNA barcodes (unfiltered): {human_format(inter)}")
    if gene_filtered is not None:
        cn.append(f"Number of cells after filtering by the minimal number of genes to consider a cell usable: {human_format(gene_filtered.shape[0])}")
    cn.append(f"Number of cells after filtering doublets: {human_format(md.gene.n_obs)}")
    blocks.append({"modality": "Filtering Summary", "subject": "Cell Number", "value_display": ", ".join(cn)})
    blocks.append({"modality": "Filtering Summary", "subject": "Gene Number",
                   "value_display": f"Number of genes detected after filtering: {human_format(md.gene.n_vars)}, Mean UMI counts per cell after filtering: "
                                    f"{human_format(float(np.asarray(md.gene.X.sum(axis=1)).mean()))}"})
    if inter is not None:
        blocks.append({"modality": "Guide", "subject": "Visualization",
                       "value_display": f"Number of guide barcodes (unfiltered) intersecting with scRNA barcodes (unfiltered): {human_format(inter)}, "
                                        f"% of guides barcodes (unfiltered) intersecting with the scRNA barcode (unfiltered): "
                                        f"{np.round(inter / guide_unfiltered.n_obs * 100, 2)}%"})
    tr = uns_df(md, "test_results", required=False)
    if tr is not None:
        targets = set(md.guide.var["intended_target_name"].astype(str))
        t5 = tr[tr["intended_target_name"].astype(str).isin(targets) & (tr["p_value"] < 0.05)]
        direct = t5[t5["pair_type"] == "Direct targeting"] if "pair_type" in tr.columns else t5.iloc[:0]
        negt = t5[t5["pair_type"] == "Targeting_negative_control"] if "pair_type" in tr.columns else t5.iloc[:0]
        n_direct_all = int((tr["pair_type"] == "Direct targeting").sum()) if "pair_type" in tr.columns else 0
        tab = tr.copy()
        tab["significant"] = tab["p_value"] < 0.05
        blocks.append({"modality": "Inference", "subject": "Guide Inference", "table": tab,
                       "value_display": f"Total tested sgRNA-gene pairs: {len(tr)}, Total tested significant sgRNA-gene pairs(p_value<0.05): {len(t5)}, "
                                        f"Total number of Direct-Targeting pairs presenting significant perturbation effects(p<0.05): {len(direct)}, "
                                        f"Percentage of total tested Direct-Targeting pairs presenting significant perturbation effects(p<0.05): "
                                        f"{np.round(len(direct) / max(len(tr), 1) * 100, 2)}% (upstream denominator: all pairs; of direct-targeting pairs: "
                                        f"{np.round(len(direct) / max(n_direct_all, 1) * 100, 2)}%), "
                                        f"Total number of Negative pairs presenting significant perturbation effect(p<0.05): {len(negt)}",
                       "stats": {"n_pairs": len(tr), "n_sig_005": len(t5), "n_sig_001": int((tr["p_value"] < 0.01).sum()), "n_direct_sig": len(direct),
                                 "n_direct": n_direct_all, "n_negative_sig": len(negt)}})
    if "guide_assignment" in md.guide.layers:
        A = dense(md.guide.layers["guide_assignment"]) > 0
        freq = pd.DataFrame({"sgRNA": md.guide.var_names, "# guide barcodes": A.sum(axis=0)})
        blocks.append({"modality": "Inference", "subject": "Guide Assignment", "table": freq,
                       "value_display": f"Number of Guide barcodes with a positive sgRNA call: {human_format(int(A.any(axis=1).sum()))}, "
                                        f"The median of cells with a positive sgRNA call is: {float(np.median(A.sum(axis=0)))}"})
    gpc = np.asarray(md.guide.X.sum(axis=1)).ravel()
    cpg = np.asarray(md.guide.X.sum(axis=0)).ravel()
    blocks.append({"modality": "Inference", "subject": "Visualization",
                   "value_display": f"Mean guides per cell: {human_format(gpc.mean())}, Mean cells per guide: {human_format(cpg.mean())} "
                                    "(sums of guide UMIs, as upstream)"})
    return blocks


def dashboard_plots(md: MuDataLite, figdir: Path) -> list:
    np = _np()
    plt = _plt()
    if plt is None:
        return []
    figdir.mkdir(parents=True, exist_ok=True)
    figs = []
    for key, vals in threshold_counts(md).items():
        fig, ax = plt.subplots(figsize=(7, 4))
        xs = list(vals)
        ax.bar(range(len(xs)), list(vals.values()), color=BLUE)
        for i, v in enumerate(vals.values()):
            ax.text(i, v, str(v), ha="center", va="bottom", fontsize=8)
        ax.set_xticks(range(len(xs)))
        ax.set_xticklabels(xs)
        _style(ax, key.replace("_", " "), "threshold", "Number of cells")
        figs.append(_save(fig, figdir / f"{key}.png"))
    if "guide_assignment" in md.guide.layers:
        A = dense(md.guide.layers["guide_assignment"]) > 0
        fr = np.sort(A.sum(axis=0))
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.bar(range(len(fr)), fr, color=GREEN)
        _style(ax, "Histogram of the number of sgRNA represented at the guide barcodes", "sgRNAs", "Number of Guide Barcodes with a Given sgRNA")
        figs.append(_save(fig, figdir / "guides_hist_num_sgRNA.png"))
    for vals, name, lab in ((np.asarray(md.guide.X.sum(axis=1)).ravel(), "guides_per_cell_histogram", "Number of Guides per Cell"),
                            (np.asarray(md.guide.X.sum(axis=0)).ravel(), "cells_per_guide_histogram", "Number of Cells per Guide")):
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.hist(vals, bins=50 if "per_cell" in name else 30, color="#8ec5e8", ec="#4682b4")
        _style(ax, f"Histogram of {lab.replace('Number of ', '')}", lab, "Density")
        figs.append(_save(fig, figdir / f"{name}.png"))
    return figs


def dashboard_html(blocks: list, figs: list, title: str) -> str:
    import html as H
    parts = [f"<!doctype html><html><head><meta charset='utf-8'><title>{H.escape(title)}</title><style>"
             ":root{--ink:#0b0b0b;--muted:#52514e;--line:#d0cfca;--bg:#fcfcfb;--hl:#eef4fb}"
             "@media (prefers-color-scheme: dark){:root{--ink:#eceae4;--muted:#b3b0a8;--line:#3b3a36;--bg:#161614;--hl:#1d2733}}"
             "body{font-family:system-ui,sans-serif;background:var(--bg);color:var(--ink);margin:0 auto;max-width:1100px;padding:16px}"
             "h2{font-size:15px;margin:24px 0 6px}.hl{background:var(--hl);padding:8px 10px;border-radius:6px;font-size:13px}"
             "table{border-collapse:collapse;font-size:12px;display:block;overflow-x:auto;max-height:360px}td,th{border:1px solid var(--line);padding:3px 6px}"
             "img{max-width:100%;height:auto;border:1px solid var(--line);margin:6px 0}</style></head><body>",
             f"<h1 style='font-size:20px'>{H.escape(title)}</h1>"]
    for b in blocks:
        parts.append(f"<h2>{H.escape(str(b.get('modality')))} — {H.escape(str(b.get('subject')))}</h2>")
        if b.get("value_display"):
            parts.append(f"<div class='hl'>{H.escape(str(b['value_display']))}</div>")
        t = b.get("table")
        if t is not None and len(t):
            parts.append(t.head(200).to_html(index=False, border=0))
    for f in figs:
        try:
            parts.append(f"<img alt='{H.escape(Path(f).stem)}' src='data:image/png;base64,{base64.b64encode(Path(f).read_bytes()).decode()}'>")
        except OSError:
            pass
    parts.append("</body></html>")
    return "\n".join(parts)


def cmd_dashboard(args: argparse.Namespace) -> int:
    pd = _pd()
    d = run_dir(args.label)
    md = read_h5mu(Path(args.mudata))
    gu = read_ad(Path(args.gene_ann)) if args.gene_ann else None
    gf = read_ad(Path(args.gene_ann_filtered)) if args.gene_ann_filtered else None
    guu = read_ad(Path(args.guide_ann)) if args.guide_ann else None
    pt = pd.read_csv(args.guide_fq_tbl) if args.guide_fq_tbl else None
    blocks = kb_json_tables([Path(x) for x in (args.kb_dirs or [])])
    blocks = dashboard_blocks(md, gu, gf, guu, pt)[:1] + sorted(blocks, key=lambda b: b["description"]) + dashboard_blocks(md, gu, gf, guu, pt)[1:] \
        if pt is not None else sorted(blocks, key=lambda b: b["description"]) + dashboard_blocks(md, gu, gf, guu, pt)
    figs = [] if args.no_plots else dashboard_plots(md, d / "figures")
    extra = [Path(x) for x in (args.extra_figures or []) if Path(x).is_file()]
    html_p = d / "dashboard.html"
    html_p.write_text(dashboard_html(blocks, figs + extra, args.title))
    print(f"Wrote: {html_p}")
    rows = pd.DataFrame([{k: v for k, v in b.items() if k not in ("table", "stats")} for b in blocks])
    write_tsv(rows, d / "dashboard_blocks.tsv")
    th = threshold_counts(md)
    write_report(d, "CRISPR Jamboree 3 — dashboard", [md_table(["modality", "subject", "highlight"], rows[["modality", "subject", "value_display"]].values.tolist()),
                                                      md_table(["plot", "threshold -> cells"], [[k, json.dumps(v)] for k, v in th.items()])],
                 {"thresholds": th, "blocks": len(blocks), "inference": next((b.get("stats") for b in blocks if b.get("stats")), None),
                  "html": html_p, "figures": figs})
    return 0


# ---------------------------------------------------------------------------
# run: the process_mudata / evaluation / dashboard sub-workflows
# ---------------------------------------------------------------------------

def cmd_run(args: argparse.Namespace) -> int:
    np, pd = _np(), _pd()
    d = run_dir(args.label)
    sections, summary = [], {"run_dir": d}
    rna_raw = read_ad(Path(args.rna))
    guide_raw = read_ad(Path(args.guide))
    gtf = read_gtf_genes(Path(args.gtf)) if args.gtf else None
    names = None
    if args.gene_names:
        names = [ln.strip() for ln in open(args.gene_names) if ln.strip()]
    elif gtf is not None:
        sym = dict(zip(strip_version(gtf["gene_id"]), gtf["gene_name"]))
        names = [sym.get(g, g) for g in strip_version(rna_raw.var_names)]
    full, filt = preprocess_rna(rna_raw, names, args.min_genes, args.min_cells, args.pct_mito, args.reference, args.upstream_compat)
    sections.append(f"preprocess: {full.n_obs} barcodes -> {filt.n_obs} cells, {filt.n_vars} genes.")
    hashing = None
    if args.hashing:
        hashing, _ = demultiplex_hashing(read_ad(Path(args.hashing)), list(filt.obs_names))
        sections.append(f"demultiplex: {hashing.n_obs} non-negative hashed cells; {dict(hashing.obs['hto_type_split'].value_counts())}.")
    md = create_mudata(filt, guide_raw, read_table(Path(args.guide_metadata)), gtf, args.moi, hashing, args.upstream_compat)
    sections.append(f"create-mudata: {md.gene.n_obs} cells x {md.gene.n_vars} genes x {md.guide.n_vars} gRNAs.")
    md, sc = remove_doublets(md, args.seed)
    sections.append(f"doublets: threshold {_fmt(sc['threshold'])}, {int(sc['predicted_doublets'].sum())} removed -> {md.gene.n_obs} cells.")
    if args.assignment_method == "cleanser":
        A, _ = assign_cleanser(md.guide, args.threshold, None, pooled=args.upstream_compat)
    elif args.assignment_method == "umi-threshold":
        A = assign_threshold(md.guide, args.threshold)
    else:
        A = assign_sceptre_mixture(md.guide, md.gene)
    md.guide.layers["guide_assignment"] = A
    sections.append(f"assign-guides: {args.assignment_method}; {int((dense(A) > 0).any(axis=1).sum())} cells with a gRNA.")
    if args.pairs:
        pairs = read_table(Path(args.pairs))
    else:
        if gtf is None:
            raise SystemExit("--pairs or --gtf is required")
        pairs = pairs_from_gtf(md.guide.var, gtf, args.distance_from_center, args.upstream_compat)
    md.uns["pairs_to_test"] = restrict_pairs(pairs, md)
    tr, used = sceptre_inference(md, args.side, args.formula, args.cov_string or "", args.resamples, args.seed, args.moi)
    md.uns["test_results"] = tr
    tr.to_csv(d / "test_results.csv", index=False)
    print(f"CSV: {d / 'test_results.csv'}")
    write_h5mu(md, d / "inference_mudata.h5mu")
    sections.append(f"infer (sceptre): {len(tr)} pairs, {int((tr['p_value'] < 0.05).sum())} at p < 0.05.")
    ev = evaluate_inference(md, tr, gtf, d, args.central_nodes, args.num_nodes, 1.0, 0.05, 0.1, args.no_plots)
    sections.append(f"evaluate: {ev['n_down']} down / {ev['n_up']} up; network {', '.join(map(str, ev['central_nodes']))}.")
    blocks = dashboard_blocks(md, full, filt, guide_raw, None)
    figs = [] if args.no_plots else dashboard_plots(md, d / "figures")
    (d / "dashboard.html").write_text(dashboard_html(blocks, figs + ev["figures"], args.label))
    print(f"Wrote: {d / 'dashboard.html'}")
    stats = next((b.get("stats") for b in blocks if b.get("stats")), {})
    summary.update({"n_cells": md.gene.n_obs, "n_doublets": int(sc["predicted_doublets"].sum()), "inference": stats, "evaluation": ev,
                    "covariates": used, "figures": figs + ev["figures"]})
    write_report(d, "CRISPR Jamboree 3 — single-cell Perturb-seq pipeline", sections, summary)
    return 0


# ---------------------------------------------------------------------------
# Synthetic world + self-test
# ---------------------------------------------------------------------------

CONFIG_TABLE_EXAMPLE = '''tab_name,variable,variable_value,batch_name,read1,read2
"GeneralConfig","run_name","Workshop Test","","",""
"GeneralConfig","transcriptome","human","","",""
"scRNA","scRNA_seqspec_yaml","rna.yml","","",""
"scRNA","sequence","","batch_a, lane1","/data/a_rna_R1.fastq.gz ","/data/a_rna_R2.fastq.gz"
"scRNA","sequence","","batch_b, lane2","/data/b_rna_R1.fastq.gz ","/data/b_rna_R2.fastq.gz"
"Guides","Guides_seqspec_yaml","guide.yml","","",""
"Guides","sequence","","batch_a, lane1","/data/a_guide_R1.fastq.gz ","/data/a_guide_R2.fastq.gz"
"Guides","sequence","","batch_b, lane2","/data/b_guide_R1.fastq.gz ","/data/b_guide_R2.fastq.gz"
"Hash","Hash_seqspec_yaml","hash.yml","","",""
"Hash","sequence","","batch_a, lane1","/data/a_hash_R1.fastq.gz ","/data/a_hash_R2.fastq.gz"
"CellGuideAssignment","THRESHOLD","1","","",""
"PerturbationInference","moi","high","","",""
"QualityFilters","pct_mito","20.5","","",""
'''

GUIDE_SEQSPEC = """!Assay
seqspec_version: 0.2.0
assay_id: MULTISEQ_10XV3_Guide
modalities:
- guide
sequence_spec:
- !Read
  read_id: guide_R1.fq
  name: Read 1
  modality: guide
  primer_id: r1_primer
  min_len: 28
  max_len: 28
  strand: pos
- !Read
  read_id: guide_R2.fq
  name: Read 2
  modality: guide
  primer_id: r2_primer
  min_len: 83
  max_len: 83
  strand: neg
library_spec:
- !Region
  parent_id: null
  region_id: guide
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


def synthetic_world(d: Path, seed: int = 11) -> dict:
    """RNA / guide / hashing count AnnDatas (2 batches) with planted QC failures, doublets, hashing, knockdowns."""
    np, pd = _np(), _pd()
    import anndata as ad  # type: ignore
    import scipy.sparse as sp  # type: ignore
    rng = np.random.default_rng(seed)
    n_good, n_empty, n_dying, n_dbl = 900, 80, 40, 60
    n_genes = 300
    targets = [f"TF{i}" for i in range(1, 7)]
    symbols = targets + [f"MT-CO{i}" for i in range(1, 6)] + ["RPS3", "RPL7", "RPS6", "RPL10"] + [f"G{i}" for i in range(n_genes - 15)]
    ens = [f"ENSG{700000 + i:08d}.{1 + i % 3}" for i in range(n_genes)]
    ens0 = strip_version(ens)
    base = rng.gamma(1.5, 1.0, n_genes) * 1.5 + 0.2
    base[6:11] = 1.5
    ctype = rng.integers(0, 3, n_good)
    prog = np.ones((3, n_genes))
    for t in range(3):
        prog[t, 15 + 40 * t: 15 + 40 * t + 40] = 6.0
    lib = rng.lognormal(0, 0.25, n_good) * 3.0
    guide_ids = [f"{t}|{''.join(rng.choice(list('ACGT'), 20))}" for t in targets for _ in range(2)] + \
                [f"non-targeting_{k:05d}|{''.join(rng.choice(list('ACGT'), 20))}" for k in range(4)]
    tname = [g.split("|")[0] for g in guide_ids]
    ng = len(guide_ids)
    A = np.zeros((n_good, ng), bool)
    for c in range(n_good):
        A[c, rng.choice(ng, size=min(ng, 1 + rng.poisson(1.0)), replace=False)] = True
    kd = {"TF1": 0.2, "TF2": 0.3, "TF3": 0.3, "TF4": 0.35, "TF5": 0.35, "TF6": 1.0}
    mu = lib[:, None] * base[None, :] * prog[ctype]
    for t, f in kd.items():
        cells = A[:, [i for i, x in enumerate(tname) if x == t]].any(axis=1)
        mu[cells, targets.index(t)] *= f
    X = rng.negative_binomial(4, 4 / (4 + mu))
    Xe = rng.negative_binomial(4, 4 / (4 + 0.03 * base[None, :] * np.ones((n_empty, 1))))
    mud = 3.0 * base[None, :] * np.ones((n_dying, 1))
    mud[:, 6:11] *= 60
    Xd = rng.negative_binomial(4, 4 / (4 + mud))
    pa = rng.integers(0, n_good, n_dbl)
    pb = np.array([rng.choice(np.where(ctype != ctype[i])[0]) for i in pa])
    Xdb = X[pa] + X[pb]
    Xall = np.vstack([X, Xe, Xd, Xdb]).astype(np.float64)
    N = Xall.shape[0]
    bcs = []
    seen = set()
    while len(bcs) < N:
        b = "".join(rng.choice(list("ACGT"), 16))
        if b not in seen:
            seen.add(b)
            bcs.append(b)
    kind = np.array(["good"] * n_good + ["empty"] * n_empty + ["dying"] * n_dying + ["doublet"] * n_dbl)
    batch = np.array(["batch_a" if i % 2 == 0 else "batch_b" for i in range(N)])
    Aall = np.vstack([A, np.zeros((n_empty + n_dying, ng), bool), A[pa] | A[pb]])
    gumi = np.where(Aall, rng.negative_binomial(3, 3 / (3 + 30.0), Aall.shape) + 2, rng.poisson(0.2, Aall.shape)).astype(np.float64)
    hto = np.zeros((N, 3), bool)
    hid = rng.integers(0, 3, N)
    hto[np.arange(N), hid] = True
    neg = rng.choice(n_good, 40, replace=False)
    hto[neg] = False
    hto[n_good + n_empty + n_dying:] = hto[pa] | hto[pb]
    for i, (a_, b_) in enumerate(zip(pa, pb)):
        r = n_good + n_empty + n_dying + i
        hto[r, hid[a_]] = True
        hto[r, (hid[a_] + 1) % 3] = True
    humi = np.where(hto, rng.negative_binomial(5, 5 / (5 + 150.0), hto.shape) + 20, rng.poisson(3, hto.shape)).astype(np.float64)
    obs = pd.DataFrame({"batch": batch}, index=bcs)
    rna = ad.AnnData(X=sp.csr_matrix(Xall), obs=obs.copy(), var=pd.DataFrame(index=ens))
    guide = ad.AnnData(X=sp.csr_matrix(gumi), obs=obs.copy(), var=pd.DataFrame(index=guide_ids))
    hashing = ad.AnnData(X=sp.csr_matrix(humi), obs=obs.copy(), var=pd.DataFrame(index=["HTO1", "HTO2", "HTO3"]))
    paths = {}
    for nm, a in (("rna", rna), ("guide", guide), ("hashing", hashing)):
        p = d / f"{nm}.h5ad"
        a.write_h5ad(str(p))
        paths[nm] = p
    (d / "gene_names.txt").write_text("\n".join(symbols) + "\n")
    gtf = d / "genes.gtf.gz"
    with gzip.open(gtf, "wt") as fh:
        fh.write("##description: synthetic\n")
        for i, (e, s) in enumerate(zip(ens, symbols)):
            st = 1_000_000 + i * 200_000
            fh.write(f"chr1\tSYN\tgene\t{st}\t{st + 10_000}\t.\t+\t.\tgene_id \"{e}\"; gene_type \"protein_coding\"; gene_name \"{s}\";\n")
            fh.write(f"chr1\tSYN\ttranscript\t{st}\t{st + 10_000}\t.\t+\t.\tgene_id \"{e}\"; transcript_id \"T{i}\"; gene_name \"{s}\";\n")
    rows = []
    for g in guide_ids:
        t = g.split("|")[0]
        if t in targets:
            st = 1_000_000 + targets.index(t) * 200_000 - 100
            rows.append({"sgRNA_ID": g, "sgRNA_sequences": revcomp(g.split("|")[1]), "Target_name": t, "chr": "chr1", "start": st, "end": st + 20, "Set": "Set A"})
        else:
            rows.append({"sgRNA_ID": g, "sgRNA_sequences": revcomp(g.split("|")[1]), "Target_name": t, "chr": "non_targeting", "start": 0, "end": 0, "Set": "Set A"})
    meta = pd.DataFrame(rows).sample(frac=1.0, random_state=1).reset_index(drop=True)  # shuffled: exercises the by-ID mapping
    meta_p = d / "guide_metadata.tsv"
    meta.to_csv(meta_p, sep="\t", index=False)
    prs = []
    for g, t in zip(guide_ids, tname):
        if t in targets:
            prs.append({"guide_id": g, "gene_name": ens0[targets.index(t)], "intended_target_name": t, "pair_type": "Direct targeting"})
            for k in (200, 210, 220):
                prs.append({"guide_id": g, "gene_name": ens0[k], "intended_target_name": t, "pair_type": "Targeting_negative_control"})
        else:
            for k in range(3):
                prs.append({"guide_id": g, "gene_name": ens0[k], "intended_target_name": t, "pair_type": "Non-targeting_negative_control"})
    prs.append({"guide_id": "missing|AAAA", "gene_name": ens0[0], "intended_target_name": "TF1", "pair_type": "Direct targeting"})
    pairs_p = d / "user_pairs_to_test.csv"
    pd.DataFrame(prs).to_csv(pairs_p, index=False)
    # FASTQs for seqspec-check: R2 = 63 nt common region + protospacer (negative-strand read, as the guide seqspec)
    common = "CAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCTTAGCTCTTAAAC"
    fq1, fq2 = d / "guide_R1.fastq.gz", d / "guide_R2.fastq.gz"
    with gzip.open(fq1, "wt") as f1, gzip.open(fq2, "wt") as f2:
        for i in range(200):
            s1 = bcs[i % N] + "".join(rng.choice(list("ACGT"), 12))
            s2 = common + (meta["sgRNA_sequences"].iloc[i % len(meta)] if i < 150 else "".join(rng.choice(list("ACGT"), 20)))
            f1.write(f"@r{i}\n{s1}\n+\n{'I' * len(s1)}\n")
            f2.write(f"@r{i}\n{s2}\n+\n{'I' * len(s2)}\n")
    kb_dirs = []
    for b in ("batch_a", "batch_b"):
        for mod, a in (("transcripts", rna), ("guide", guide)):
            kd_ = d / f"{b}_ks_{mod}_out"
            (kd_ / "counts_unfiltered").mkdir(parents=True)
            sub = a[a.obs["batch"].values == b].copy()
            del sub.obs["batch"]
            sub.write_h5ad(str(kd_ / "counts_unfiltered" / "adata.h5ad"))
            (kd_ / "inspect.json").write_text(json.dumps({"numRecords": 123456, "numReads": 150000, "numBarcodes": int(sub.n_obs),
                                                           "meanReadsPerBarcode": 120.44, "percentageBarcodesOnOnlist": 97.123}))
            (kd_ / "run_info.json").write_text(json.dumps({"n_targets": 300, "n_processed": 150000, "n_pseudoaligned": 120000,
                                                            "n_unique": 110000, "p_pseudoaligned": 80.0, "p_unique": 73.3,
                                                            "call": "kallisto bus ...", "start_time": "now"}))
            kb_dirs.append(kd_)
    cov_p = d / "parse_covariate.csv"
    pd.DataFrame({"batch": ["batch_a", "batch_b"], "cov1": ["lane1", "lane2"]}).to_csv(cov_p, index=False)
    return {**paths, "gtf": gtf, "meta": meta_p, "pairs": pairs_p, "fq1": fq1, "fq2": fq2, "kb_dirs": kb_dirs, "cov": cov_p,
            "kind": kind, "bcs": bcs, "ens0": ens0, "targets": targets, "guide_ids": guide_ids, "hto": hto, "A": Aall, "n_pairs_valid": len(prs) - 1,
            "gene_names": d / "gene_names.txt", "n_good": n_good}


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
            W = synthetic_world(d)
            kind = W["kind"]
            bc_kind = dict(zip(W["bcs"], kind))

            print("\nconfigure")
            cp = d / "configuration.csv"
            cp.write_text(CONFIG_TABLE_EXAMPLE)
            rc = cmd_configure(argparse.Namespace(config_table=str(cp), mapping_tabs=["scRNA", "Guides", "Hash"], label="st_j3_cfg"))
            r = last("st_j3_cfg")
            text = (r / "pipeline_input.config").read_text()
            check(rc == 0 and "DATASET_HASHING = 'true'" in text and "run_name = 'Workshop Test'" in text and "THRESHOLD = 1\n" in text
                  and "pct_mito = 20.5" in text and "distance_from_center = 1000000" in text,
                  "configure: quoted strings, int / float values, DATASET_HASHING from the Hash tab, distance_from_center default")
            check('fastq_files_rna = [\n    "/data/a_rna_R1.fastq.gz  /data/a_rna_R2.fastq.gz",\n    "/data/b_rna_R1.fastq.gz  /data/b_rna_R2.fastq.gz"\n    ]' in text
                  and "test_guide_fastq_r1 = ['/data/a_guide_R1.fastq.gz', '/data/b_guide_R1.fastq.gz']" in text
                  and text.endswith("params.covariate_list = [\n    batch: ['batch_a', 'batch_b'],\n    cov1: ['lane1', 'lane2']\n]")
                  and "batch=['batch_a', 'batch_b']" in text, "configure: per-batch FASTQ groups, test FASTQ lists, covariate_list block")
            check((r / "cov_string.txt").read_text() == "batch + cov1", "prepare_formula: covariates with > 1 level -> 'batch + cov1'")
            check(process_batches("a.fq,b.fq;c.fq|d.fq", "/x") == ["/x/a.fq /x/b.fq", "/x/c.fq /x/d.fq"]
                  and process_reads("/x", "a1.fq; b1.fq", "a2.fq;b2.fq") == "/x/a1.fq /x/a2.fq /x/b1.fq /x/b2.fq", "process_batches / process_reads")

            print("\nseqspec")
            rc = cmd_seqspec_check(argparse.Namespace(read1=[str(W["fq1"])], read2=[str(W["fq2"])], metadata=str(W["meta"]), max_reads=100000,
                                                      no_plots=args.no_plots, label="st_j3_ssc"))
            r = last("st_j3_ssc")
            pt = pd.read_csv(r / "position_table.csv")
            freq = pd.read_csv(r / "nucleotide_frequency.guide_R2.tsv", sep="\t")
            check(rc == 0 and pt.iloc[0, 0] == "read name" and str(pt.iloc[1, 0]) == "63" and int(pt.iloc[1, 1]) == 150,
                  "seqspec-check: 150 planted guide reads found at R2 position 63 (after the 63 nt common region)")
            check(abs(freq.loc[0, "C"] - 1.0) < 1e-12 and abs(freq.iloc[:, 1:].sum(axis=1).iloc[70] - 1) < 1e-12, "seqspec-check: per-position A/C/G/T frequencies")
            yp = d / "guide.yml"
            yp.write_text(GUIDE_SEQSPEC)
            rc = cmd_seqspec_parse(argparse.Namespace(yaml=str(yp), modality=["guide"], directory=None, label="st_j3_ssp"))
            r = last("st_j3_ssp")
            ps = pd.read_csv(r / "guide_parsed_seqSpec.txt", sep="\t")
            check(rc == 0 and ps.loc[0, "representation"] == "0,0,16:0,16,28:1,63,83" and ps.loc[0, "barcode_whitelist"] == "737K-august-2016.txt",
                  "seqspec-parse: kb representation and whitelist filename (parsed_seqSpec.txt)")
            rc = cmd_map(argparse.Namespace(modality="guide", seqspec_yaml=str(yp), chemistry=None, whitelist=None, metadata=str(W["meta"]),
                                            fastqs=["batch_a:a_R1.fq.gz a_R2.fq.gz"], genome="genome.fa.gz", species="human", threads=4,
                                            execute=False, label="st_j3_map"))
            r = last("st_j3_map")
            sh = (r / "map_commands.sh").read_text()
            ft = pd.read_csv(r / "guide_features.txt", sep="\t", header=None)
            check(rc == 0 and "--workflow kite" in sh and "-x 0,0,16:0,16,28:1,63,83 -o batch_a_ks_guide_out" in sh and "-m 20G" in sh
                  and list(ft.iloc[0]) == [pd.read_csv(W["meta"], sep="\t").loc[0, "sgRNA_sequences"], pd.read_csv(W["meta"], sep="\t").loc[0, "sgRNA_ID"]],
                  "map: kb ref (kite) + kb count commands; guide_features.txt = sequence, ID")

            print("\nconcat / preprocess")
            rc = cmd_concat(argparse.Namespace(inputs=[str(p) for p in W["kb_dirs"] if "transcripts" in p.name], covariates=str(W["cov"]),
                                               out=str(d / "rna_concat.h5ad"), label="st_j3_concat"))
            last("st_j3_concat")
            rc_ = read_ad(d / "rna_concat.h5ad")
            check(rc == 0 and rc_.n_obs == len(W["bcs"]) and set(rc_.obs["batch"]) == {"batch_a", "batch_b"} and set(rc_.obs["cov1"]) == {"lane1", "lane2"}
                  and rc_.obs_names[0].endswith("_0"), "concat: batch from the '<batch>_ks_' prefix, covariates joined, index_unique '_'")
            rna = read_ad(W["rna"])
            names = [ln.strip() for ln in open(W["gene_names"]) if ln.strip()]
            full, filt = preprocess_rna(rna, names, 100, 3, 20.0, "human")
            kept = set(filt.obs_names)
            check(all(bc_kind[b] in ("good", "doublet") for b in kept) and sum(bc_kind[b] == "good" for b in kept) == W["n_good"],
                  f"preprocess: {full.n_obs} -> {filt.n_obs} cells; all empty (n_genes < 100) and dying (pct_counts_mt >= 20) barcodes removed, all good cells kept")
            check(list(filt.var_names[:2]) == W["ens0"][:2] and "pct_counts_mt" in filt.obs.columns and "log1p_total_counts" in filt.obs.columns
                  and "n_cells_by_counts" in filt.var.columns and set(full.obs["batch_number"]) == {1, 2},
                  "preprocess: Ensembl version stripped; calculate_qc_metrics columns; batch_number")
            try:
                import scanpy as sc  # type: ignore
                ref = rna.copy()
                ref.var["mt"] = [s.startswith("MT-") for s in names]
                ref.var["ribo"] = [s.startswith(("RPS", "RPL")) for s in names]
                sc.pp.calculate_qc_metrics(ref, qc_vars=["mt", "ribo"], percent_top=(50, 100, 200), inplace=True, log1p=True)
                check(np.allclose(ref.obs["pct_counts_mt"].values, full.obs["pct_counts_mt"].values) and
                      np.allclose(ref.obs["pct_counts_in_top_50_genes"].values, full.obs["pct_counts_in_top_50_genes"].values) and
                      np.allclose(ref.var["pct_dropout_by_counts"].values, full.var["pct_dropout_by_counts"].values),
                      "preprocess: QC metrics identical to sc.pp.calculate_qc_metrics")
            except ImportError:
                pass
            mouse = rna.copy()
            _, fm = preprocess_rna(mouse, names, 100, 3, 20.0, "mouse")
            _, fmc = preprocess_rna(mouse, names, 100, 3, 20.0, "mouse", compat=True)
            check(fm.n_obs > fmc.n_obs, "preprocess: mouse uses the 'Mt-' prefix (MT- genes then not mitochondrial); --upstream-compat keeps 'MT-'")
            rc = cmd_preprocess(argparse.Namespace(adata=str(W["rna"]), gene_names=str(W["gene_names"]), min_genes=100, min_cells=3, pct_mito=20.0,
                                                   reference="human", upstream_compat=False, out=str(d / "filtered_anndata.h5ad"),
                                                   no_plots=args.no_plots, label="st_j3_pre"))
            last("st_j3_pre")
            check(rc == 0 and read_ad(d / "filtered_anndata.h5ad").n_obs == filt.n_obs, "preprocess: subcommand writes filtered_anndata.h5ad")

            print("\ndemultiplex")
            demux, reps = demultiplex_hashing(read_ad(W["hashing"]), list(filt.obs_names))
            obs = demux.obs
            dbl_bcs = [b for b in demux.obs_names if bc_kind[b] == "doublet"]
            good_bcs = [b for b in demux.obs_names if bc_kind[b] == "good"]
            check((obs.loc[dbl_bcs, "hto_type_split"] == "multiplets").mean() > 0.95 and (obs.loc[good_bcs, "hto_type_split"] != "multiplets").mean() > 0.97,
                  f"GMM demux: {(obs.loc[dbl_bcs, 'hto_type_split'] == 'multiplets').mean():.0%} of planted doublets called multiplets")
            n_neg_planted = sum(1 for b in filt.obs_names if bc_kind[b] == "good") - len(good_bcs)
            check(35 <= n_neg_planted <= 45 and "negative" not in set(obs["hto_type"]), f"GMM demux: {n_neg_planted} negatives removed (40 planted)")
            rc = cmd_demultiplex(argparse.Namespace(hashing=str(W["hashing"]), rna=str(d / "filtered_anndata.h5ad"), out=str(d / "demux.h5ad"), label="st_j3_demux"))
            r = last("st_j3_demux")
            check(rc == 0 and (r / "batch_a_GMM_full.csv").is_file() and (r / "batch_a_GMM_full.config").is_file(), "demultiplex: GMM_full.csv / .config per batch")

            print("\ncreate-mudata / doublets")
            meta = read_table(W["meta"])
            gtf = read_gtf_genes(W["gtf"])
            md = create_mudata(filt, read_ad(W["guide"]), meta, gtf, "high", demux)
            gv = md.guide.var
            tf1 = gv[gv["intended_target_name"] == "TF1"].iloc[0]
            check(set(md.mods) == {"gene", "guide", "hashing"} and md.gene.n_obs == demux.n_obs and tf1["intended_target_chr"] == "chr1"
                  and tf1.name == f"{tf1['guide_id']}|{meta.set_index('sgRNA_ID').loc[tf1['guide_id'], 'sgRNA_sequences']}",
                  "create-mudata: 3 modalities on the intersected barcodes; metadata mapped by sgRNA_ID; var_names sgRNA_ID|sequence")
            check((gv["targeting"] == "FALSE").sum() == 4 and md.gene.var.loc[W["ens0"][0], "gene_start"] == 1_000_000
                  and {"total_gene_umis", "num_expressed_genes", "percent_mito"} <= set(md.gene.obs.columns) and "batch" in md.obs.columns
                  and _uns_scalar(md.guide, "moi") == "high",
                  "create-mudata: non-targeting guides targeting=FALSE, GTF coordinates, renamed obs, shared obs, uns moi")
            mdc = create_mudata(filt, read_ad(W["guide"]), meta, gtf, "high", None, compat=True)
            gvc = mdc.guide.var
            mism = sum(1 for g, t in zip(gvc["guide_id"], gvc["intended_target_name"]) if g.split("|")[0] != t)
            check((gvc["targeting"] == "TRUE").all() and mism > 0, f"create-mudata --upstream-compat: all TRUE and positional metadata ({mism} guides get another guide's target)")
            md2, sc_ = remove_doublets(md)
            flagged = [b for b, p in zip(md.gene.obs_names, sc_["predicted_doublets"]) if p]
            tp = sum(bc_kind[b] == "doublet" for b in flagged)
            n_dbl_in = sum(bc_kind[b] == "doublet" for b in md.gene.obs_names)
            check(sc_["threshold"] is not None and tp >= 0.5 * n_dbl_in and tp >= 0.6 * max(len(flagged), 1),
                  f"Scrublet: threshold {sc_['threshold']:.3f}; {tp}/{n_dbl_in} planted doublets flagged, precision {tp / max(len(flagged), 1):.2f}")
            check(md2.gene.n_obs == md.gene.n_obs - len(flagged) and md2.guide.n_obs == md2.gene.n_obs and "doublet_scores" in md2.gene.obs.columns,
                  "doublets: flagged cells removed from every modality")
            p_md = d / "mudata.h5mu"
            write_h5mu(md, p_md)
            rc = cmd_doublets(argparse.Namespace(mudata=str(p_md), seed=0, out=str(d / "mdata_doublets.h5mu"), no_plots=args.no_plots, label="st_j3_dbl"))
            last("st_j3_dbl")
            rc2 = cmd_create_mudata(argparse.Namespace(rna=str(d / "filtered_anndata.h5ad"), guide=str(W["guide"]), hashing=str(d / "demux.h5ad"),
                                                       guide_metadata=str(W["meta"]), gtf=str(W["gtf"]), moi="high", upstream_compat=False,
                                                       out=str(d / "mudata2.h5mu"), no_plots=args.no_plots, label="st_j3_cmd"))
            last("st_j3_cmd")
            check(rc == 0 and rc2 == 0 and read_h5mu(d / "mdata_doublets.h5mu").gene.n_obs == md2.gene.n_obs, "create-mudata / doublets subcommands")

            print("\nguide assignment")
            truth = pd.DataFrame(W["A"], index=W["bcs"], columns=W["guide_ids"]).loc[list(md2.guide.obs_names), list(md2.guide.var["guide_id"])].values
            A_c = dense(assign_cleanser(md2.guide, 0.5)[0]) > 0
            A_p = dense(assign_cleanser(md2.guide, 0.5, pooled=True)[0]) > 0
            A_s = dense(assign_sceptre_mixture(md2.guide, md2.gene)) > 0
            check((A_c == truth).mean() > 0.98 and (A_s == truth).mean() > 0.97,
                  f"assign-guides: per-guide CLEANSER {(A_c == truth).mean():.3f}, sceptre mixture {(A_s == truth).mean():.3f} agreement")
            check((A_p == truth).mean() < (A_c == truth).mean(), f"assign-guides --upstream-compat (pooled fit, within-guide PZi index): {(A_p == truth).mean():.3f} agreement")
            md2.guide.layers["guide_assignment"] = assign_cleanser(md2.guide, 0.5)[0]
            p_md2 = d / "assigned.h5mu"
            write_h5mu(md2, p_md2)
            rc = cmd_assign(argparse.Namespace(mudata=str(d / "mdata_doublets.h5mu"), method="umi-threshold", threshold=2.0, upstream_compat=False,
                                               out=str(d / "assigned_thr.h5mu"), label="st_j3_assign"))
            last("st_j3_assign")
            check(rc == 0 and "guide_assignment" in read_h5mu(d / "assigned_thr.h5mu").guide.layers, "assign-guides: subcommand")

            print("\npairs + inference")
            rc = cmd_prepare_inference(argparse.Namespace(mudata=str(p_md2), pairs=str(W["pairs"]), gtf=None, distance_from_center=1_000_000,
                                                          upstream_compat=False, out=str(d / "mudata_inference_input.h5mu"), label="st_j3_prep"))
            last("st_j3_prep")
            mdi = read_h5mu(d / "mudata_inference_input.h5mu")
            ptt = uns_df(mdi, "pairs_to_test")
            check(rc == 0 and len(ptt) == W["n_pairs_valid"] and "gene_id" in ptt.columns and "gene_name" not in ptt.columns,
                  f"prepare-inference: {len(ptt)} user pairs kept (unknown guide dropped), gene_name -> gene_id")
            gp = pairs_from_gtf(md2.guide.var, gtf, 1_000_000)
            gpc = pairs_from_gtf(md2.guide.var, gtf, 1_000_000, compat=True)
            tf1g = gp[gp["intended_target_name"] == "TF1"]
            check(set(tf1g["gene_name"]) == {W["ens0"][i] for i in range(5)} and len(gp[gp["intended_target_name"].str.startswith("non-targeting")]) == 0
                  and gpc["gene_name"].iloc[0] in {"TF1", "TF2", "TF3", "TF4", "TF5", "TF6"},
                  "create_pairs_to_test: genes within 1 Mb of the guide (Ensembl ids; symbols with --upstream-compat)")
            tr, used = sceptre_inference(mdi, B=299)
            direct = tr[tr["pair_type"] == "Direct targeting"]
            hits = direct.groupby("intended_target_name")["p_value"].max()
            negs = tr[tr["pair_type"] != "Direct targeting"]["p_value"]
            check(all(hits[t] < 1e-4 for t in ("TF1", "TF2", "TF3", "TF4", "TF5")) and hits["TF6"] > 0.01 and (negs < 0.01).mean() < 0.05,
                  f"sceptre: direct targeting TF1-5 p <= {max(hits[t] for t in ('TF1', 'TF2', 'TF3', 'TF4', 'TF5')):.1e}, unperturbed TF6 p {hits['TF6']:.2f}; "
                  f"{(negs < 0.01).mean():.1%} of controls at p < 0.01")
            check(len(tr) == len(ptt) and tr.groupby(["intended_target_name", "gene_id"])["p_value"].nunique().max() == 1
                  and abs(direct[direct["intended_target_name"] == "TF1"]["log2_fc"].iloc[0] - math.log2(0.2)) < 0.4,
                  "infer: element-level results left-joined onto every guide row; TF1 log2_fc ~ log2(0.2)")
            data = pd.DataFrame({"batch": list("aabb"), "lane": list("xxyy"), "chip": list("1122"), "one": list("zzzz"), "other": list("pqpq")})
            check(identify_non_redundant_covariates(data, "batch + lane + chip + one + other") == "batch + other"
                  and identify_non_redundant_covariates(data, "batch") == "batch", "identify_non_redundant_covariates: nested covariates collapsed, single-level dropped")
            tr2, used2 = sceptre_inference(mdi, formula="custom", cov_string="batch", B=199)
            check(used2 == "batch" and tr2["p_value"].notna().all(), "infer: formula with the batch covariate")
            rc = cmd_infer(argparse.Namespace(mudata=str(d / "mudata_inference_input.h5mu"), method="sceptre", side="both", formula="default",
                                              cov_string=None, cov_string_file=None, resamples=199, seed=4, moi=None,
                                              out=str(d / "inference_mudata.h5mu"), label="st_j3_infer"))
            r = last("st_j3_infer")
            check(rc == 0 and (r / "test_results.csv").is_file() and "test_results" in read_h5mu(d / "inference_mudata.h5mu").uns, "infer: test_results.csv + inference_mudata.h5mu")
            rc = cmd_infer(argparse.Namespace(mudata=str(d / "mudata_inference_input.h5mu"), method="perturbo", side="both", formula="default",
                                              cov_string=None, cov_string_file=None, resamples=199, seed=4, moi=None, out=None, label="st_j3_ptb"))
            last("st_j3_ptb")
            check(rc == 0, "infer perturbo: skipped cleanly when the package is absent")

            print("\nevaluate / dashboard")
            rc = cmd_evaluate(argparse.Namespace(mudata=str(d / "inference_mudata.h5mu"), gtf=str(W["gtf"]), central_nodes=None, num_nodes=2,
                                                 log2_fc=1.0, p_value=0.05, min_weight=0.1, no_plots=args.no_plots, label="st_j3_eval"))
            r = last("st_j3_eval")
            s = json.loads((r / "summary.json").read_text())
            bg = pd.read_csv(r / "evaluation_output" / "output.bedgraph", sep="\t", header=None)
            down = pd.read_csv(r / "evaluation_output" / "volcano_down.tsv", sep="\t")
            check(rc == 0 and set(down["intended_target_name"]) >= {"TF1", "TF2"} and set(s["central_nodes"]) <= set(W["targets"]) and len(s["central_nodes"]) == 2
                  and len(bg) == 12 and bg.shape[1] == 5, f"evaluate: volcano down {sorted(set(down['intended_target_name']))}; top-2 central nodes {s['central_nodes']}; "
                  "12 promoter bedgraph rows (target == GTF gene_name), headerless")
            rc = cmd_dashboard(argparse.Namespace(mudata=str(d / "inference_mudata.h5mu"), gene_ann=str(W["rna"]), gene_ann_filtered=str(d / "filtered_anndata.h5ad"),
                                                  guide_ann=str(W["guide"]), guide_fq_tbl=str(last("st_j3_ssc") / "position_table.csv"),
                                                  kb_dirs=[str(p) for p in W["kb_dirs"]], extra_figures=None, title="synthetic", no_plots=args.no_plots,
                                                  label="st_j3_dash"))
            r = last("st_j3_dash")
            s = json.loads((r / "summary.json").read_text())
            blocks = pd.read_csv(r / "dashboard_blocks.tsv", sep="\t")
            mdd = read_h5mu(d / "inference_mudata.h5mu")
            tot = np.asarray(mdd.gene.X.sum(axis=1)).ravel()
            check(rc == 0 and s["thresholds"]["scRNA_barcodes_UMI_thresholds"]["200"] == int((tot > 200).sum()) and s["inference"]["n_direct_sig"] == 10
                  and (blocks["subject"] == "Mapping scRNA").sum() == 2 and (r / "dashboard.html").is_file() and human_format(1234567) == "1.23M",
                  "dashboard: threshold counts, 10 significant direct-targeting pairs, kb mapping blocks, dashboard.html, human_format")

            print("\nrun")
            rc = cmd_run(argparse.Namespace(rna=str(W["rna"]), gene_names=None, guide=str(W["guide"]), hashing=str(W["hashing"]),
                                            guide_metadata=str(W["meta"]), gtf=str(W["gtf"]), pairs=str(W["pairs"]), moi="high", min_genes=100, min_cells=3,
                                            pct_mito=20.0, reference="human", assignment_method="cleanser", threshold=0.5, side="both", formula="default",
                                            cov_string=None, resamples=199, seed=4, distance_from_center=1_000_000, central_nodes=None, num_nodes=2,
                                            upstream_compat=False, no_plots=args.no_plots, label="st_j3_run"))
            r = last("st_j3_run")
            s = json.loads((r / "summary.json").read_text())
            check(rc == 0 and s["inference"]["n_direct_sig"] == 10 and s["n_doublets"] > 0 and (r / "dashboard.html").is_file(),
                  f"run: preprocess -> demux -> MuData -> doublets ({s['n_doublets']}) -> CLEANSER -> pairs -> sceptre -> evaluate -> dashboard")
    finally:
        for p in sorted(OUT_ROOT.glob("*_st_j3_*")):
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
    ap = argparse.ArgumentParser(prog="igvfagent crispr-jamboree3", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("configure", help="configuration.csv -> pipeline_input.config + covariates + formula")
    p.add_argument("--config-table", required=True)
    p.add_argument("--mapping-tabs", nargs="+", default=["scRNA", "Guides", "Hash"])
    p.add_argument("--label", default="configure")
    p.set_defaults(func=cmd_configure)

    p = sub.add_parser("seqspec-check", help="nucleotide frequencies + guide position table from FASTQs")
    p.add_argument("--read1", nargs="+", required=True)
    p.add_argument("--read2", nargs="+", required=True)
    p.add_argument("--metadata", required=True, help="guide (or hashing) metadata; sequences in the 2nd column")
    p.add_argument("--max-reads", type=int, default=100000)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="seqspec_check")
    p.set_defaults(func=cmd_seqspec_check)

    p = sub.add_parser("seqspec-parse", help="seqspec YAML -> parsed_seqSpec.txt (kb representation, whitelist)")
    p.add_argument("--yaml", required=True)
    p.add_argument("--modality", nargs="+", required=True, help="rna | guide | hashing")
    p.add_argument("--directory")
    p.add_argument("--label", default="seqspec_parse")
    p.set_defaults(func=cmd_seqspec_parse)

    p = sub.add_parser("map", help="kb ref / kb count commands (run with --execute when kb is on PATH)")
    p.add_argument("--modality", choices=["rna", "guide", "hashing"], required=True)
    p.add_argument("--seqspec-yaml")
    p.add_argument("--chemistry", help="kb -x string (default: from --seqspec-yaml)")
    p.add_argument("--whitelist")
    p.add_argument("--metadata", help="guide / hashing spreadsheet for the kite feature table")
    p.add_argument("--fastqs", nargs="+", help="'batch:R1 R2 ...' per batch")
    p.add_argument("--genome", default="genome.fa.gz")
    p.add_argument("--species", default="human")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--execute", action="store_true")
    p.add_argument("--label", default="map")
    p.set_defaults(func=cmd_map)

    p = sub.add_parser("concat", help="concatenate per-batch kb AnnDatas with covariates")
    p.add_argument("--inputs", nargs="+", required=True, help="<batch>_ks_*_out directories or .h5ad files")
    p.add_argument("--covariates", help="parse_covariate.csv (first column = batch)")
    p.add_argument("--out")
    p.add_argument("--label", default="concat")
    p.set_defaults(func=cmd_concat)

    p = sub.add_parser("preprocess", help="scRNA QC metrics + filters (min genes, min cells, pct mito)")
    p.add_argument("--adata", required=True)
    p.add_argument("--gene-names", help="cells_x_genes.genes.names.txt")
    p.add_argument("--min-genes", type=int, default=500)
    p.add_argument("--min-cells", type=int, default=3)
    p.add_argument("--pct-mito", type=float, default=20.0)
    p.add_argument("--reference", default="human")
    p.add_argument("--upstream-compat", action="store_true", help="always use the 'MT-' prefix")
    p.add_argument("--out")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="preprocess")
    p.set_defaults(func=cmd_preprocess)

    p = sub.add_parser("create-mudata", help="RNA + guide [+ hashing] AnnDatas + guide metadata + GTF -> MuData")
    p.add_argument("--rna", required=True)
    p.add_argument("--guide", required=True)
    p.add_argument("--hashing")
    p.add_argument("--guide-metadata", required=True)
    p.add_argument("--gtf")
    p.add_argument("--moi", default="high")
    p.add_argument("--upstream-compat", action="store_true", help="positional metadata assignment, targeting TRUE for all")
    p.add_argument("--out")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="create_mudata")
    p.set_defaults(func=cmd_create_mudata)

    p = sub.add_parser("doublets", help="Scrublet doublet scores; doublets removed")
    p.add_argument("--mudata", required=True)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="doublets")
    p.set_defaults(func=cmd_doublets)

    p = sub.add_parser("demultiplex", help="GMM demultiplexing of the hashing library; negatives dropped")
    p.add_argument("--hashing", required=True)
    p.add_argument("--rna", help="filtered RNA AnnData (barcodes to intersect, filter_hashing.py)")
    p.add_argument("--out")
    p.add_argument("--label", default="demultiplex")
    p.set_defaults(func=cmd_demultiplex)

    p = sub.add_parser("assign-guides", help="CLEANSER / sceptre mixture / UMI threshold guide assignment")
    p.add_argument("--mudata", required=True)
    p.add_argument("--method", choices=["cleanser", "sceptre", "umi-threshold"], default="cleanser")
    p.add_argument("--threshold", type=float)
    p.add_argument("--upstream-compat", action="store_true", help="CLEANSER: one pooled fit, as cleanser.py")
    p.add_argument("--out")
    p.add_argument("--label", default="assign_guides")
    p.set_defaults(func=cmd_assign)

    p = sub.add_parser("prepare-inference", help="pairs_to_test from user pairs or GTF distance")
    p.add_argument("--mudata", required=True)
    p.add_argument("--pairs", help="user pairs CSV: guide_id, gene_name, intended_target_name, pair_type")
    p.add_argument("--gtf")
    p.add_argument("--distance-from-center", type=int, default=1_000_000)
    p.add_argument("--upstream-compat", action="store_true", help="GTF pairs carry gene symbols, as create_pairs_to_test.py")
    p.add_argument("--out")
    p.add_argument("--label", default="prepare_inference")
    p.set_defaults(func=cmd_prepare_inference)

    p = sub.add_parser("infer", help="sceptre (re-implemented) or PerTurbo inference")
    p.add_argument("--mudata", required=True)
    p.add_argument("--method", choices=["sceptre", "perturbo"], default="sceptre")
    p.add_argument("--side", choices=["left", "right", "both"], default="both")
    p.add_argument("--formula", default="default", help="'default' or anything else to add the cov_string covariates")
    p.add_argument("--cov-string")
    p.add_argument("--cov-string-file")
    p.add_argument("--resamples", type=int, default=499)
    p.add_argument("--seed", type=int, default=4)
    p.add_argument("--moi", choices=["low", "high"])
    p.add_argument("--out")
    p.add_argument("--label", default="infer")
    p.set_defaults(func=cmd_infer)

    p = sub.add_parser("evaluate", help="volcano + network + IGV tracks")
    p.add_argument("--mudata", required=True)
    p.add_argument("--gtf")
    p.add_argument("--central-nodes", nargs="+")
    p.add_argument("--num-nodes", type=int, default=2)
    p.add_argument("--log2-fc", type=float, default=1.0)
    p.add_argument("--p-value", type=float, default=0.05)
    p.add_argument("--min-weight", type=float, default=0.1)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="evaluate")
    p.set_defaults(func=cmd_evaluate)

    p = sub.add_parser("dashboard", help="dashboard plots, highlight blocks and dashboard.html")
    p.add_argument("--mudata", required=True)
    p.add_argument("--gene-ann")
    p.add_argument("--gene-ann-filtered")
    p.add_argument("--guide-ann")
    p.add_argument("--guide-fq-tbl")
    p.add_argument("--kb-dirs", nargs="+")
    p.add_argument("--extra-figures", nargs="+")
    p.add_argument("--title", default="CRISPR Perturb-seq pipeline dashboard")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="dashboard")
    p.set_defaults(func=cmd_dashboard)

    p = sub.add_parser("run", help="count AnnDatas -> ... -> dashboard")
    p.add_argument("--rna", required=True)
    p.add_argument("--guide", required=True)
    p.add_argument("--hashing")
    p.add_argument("--guide-metadata", required=True)
    p.add_argument("--gene-names", help="cells_x_genes.genes.names.txt (default: symbols from the GTF)")
    p.add_argument("--gtf")
    p.add_argument("--pairs")
    p.add_argument("--moi", default="high")
    p.add_argument("--min-genes", type=int, default=500)
    p.add_argument("--min-cells", type=int, default=3)
    p.add_argument("--pct-mito", type=float, default=20.0)
    p.add_argument("--reference", default="human")
    p.add_argument("--assignment-method", choices=["cleanser", "sceptre", "umi-threshold"], default="sceptre")
    p.add_argument("--threshold", type=float, default=1.0)
    p.add_argument("--side", choices=["left", "right", "both"], default="both")
    p.add_argument("--formula", default="default")
    p.add_argument("--cov-string")
    p.add_argument("--resamples", type=int, default=499)
    p.add_argument("--seed", type=int, default=4)
    p.add_argument("--distance-from-center", type=int, default=1_000_000)
    p.add_argument("--central-nodes", nargs="+")
    p.add_argument("--num-nodes", type=int, default=2)
    p.add_argument("--upstream-compat", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--label", default="jamboree3_run")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("selftest", help="synthetic data, all subcommands, assertions")
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_selftest)

    args = ap.parse_args(argv)
    logp = setup_logging()
    print(f"Log: {logp}")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
