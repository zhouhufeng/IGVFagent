"""10x Multiome analytical methods — clean-room Python implementations.

Implements the analysis steps the 10x Genomics `analysis_guides` repo
covers (and the canonical Stuart-lab Signac vignette extends), translated
to Python with Apache-2.0-compatible dependencies.

References (algorithmic only — no source copied):
  - 10XGenomics/analysis_guides  (no LICENSE — algorithm spec only)
      Epi_multiome_QC_tutorial.ipynb
      Epi_multiome_QC_interactive_visualization.ipynb
  - Stuart-lab Signac (MIT) — TSSEnrichment, NucleosomeSignal,
                                FindMultiModalNeighbors
  - Seurat 5 WNN (Hao et al. Cell 2021)
  - Ma et al. Cell 2020 — SHARE-seq joint QC thresholds (reused here)

Subcommands wired into the `multiome` CLI:
  multiome qc-atac     — Per-barcode ATAC QC: fragments, TSS enrichment,
                         nucleosome signal, FRIP
  multiome joint-qc    — Join ATAC + RNA per-barcode tables and apply
                         Signac-convention thresholds
  multiome lsi         — TF-IDF + truncated SVD on a peak × cell matrix
                         (drops dimension 1, which captures depth)
  multiome wnn         — Joint weighted-nearest-neighbor embedding +
                         UMAP via muon (BSD-3)
  multiome peak2gene   — Per-peak Pearson/Spearman correlation with
                         nearby genes (configurable window)
  multiome showcase    — End-to-end pipeline + 6-panel composite figure
                         + narrative report

License: Apache-2.0. Heavy deps imported lazily (numpy, pandas, scipy,
anndata, muon, pyranges) so the metadata-only `multiome retrieve`
command keeps working on a bare interpreter.
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
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "Multiome10x"
PLOT_DIR = REPORT_DIR / "Plots"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

sys.path.insert(0, str(Path(__file__).resolve().parent))


# ─── Setup ──────────────────────────────────────────────────────────────────

def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"multiome_10x_analyze_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log), logging.StreamHandler(sys.stdout)],
    )
    return log


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in s)


def _require(name: str, hint: str = "") -> Any:
    # importlib.import_module returns the requested submodule (e.g. "scipy.stats"),
    # whereas __import__ returns the top-level package ("scipy") — which silently
    # breaks attribute access like scipy.stats.spearmanr.
    import importlib
    try:
        return importlib.import_module(name)
    except ImportError as exc:
        raise SystemExit(
            f"Missing dependency '{name}'. {hint}\nInstall with: pip install {name}"
        ) from exc


def save_json(label: str, data: Any) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    p = DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}.json"
    p.write_text(json.dumps(data, indent=2, sort_keys=True, default=str))
    return p


# ─── Fragments + bed I/O ────────────────────────────────────────────────────

def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def iter_fragments(path: Path) -> Iterable[tuple[str, int, int, str, int]]:
    """Yield (chrom, start, end, barcode, count) from a Cell Ranger / SHARE-seq
    style fragments BED. Skips comments and short lines.
    """
    with _open_text(path) as fh:
        for ln in fh:
            if not ln or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                yield parts[0], int(parts[1]), int(parts[2]), parts[3], (
                    int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 1
                )
            except ValueError:
                continue


def load_tss(path: Path) -> dict[str, list[int]]:
    """Load TSS positions per chromosome (one BED entry per TSS).

    Accepts (a) 3-column BED with chrom/start/end (uses start as TSS),
    or (b) 6-column BED with strand (uses start if + strand, end-1 if -).
    """
    tss: dict[str, list[int]] = defaultdict(list)
    with _open_text(path) as fh:
        for ln in fh:
            if not ln or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            strand = parts[5] if len(parts) > 5 else "+"
            tss[chrom].append(start if strand != "-" else end - 1)
    for chrom in tss:
        tss[chrom].sort()
    return tss


def load_peaks(path: Path) -> dict[str, list[tuple[int, int]]]:
    peaks: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with _open_text(path) as fh:
        for ln in fh:
            if not ln or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                peaks[parts[0]].append((int(parts[1]), int(parts[2])))
            except ValueError:
                continue
    for c in peaks:
        peaks[c].sort()
    return peaks


def _bisect_hit(intervals: list[tuple[int, int]], pos: int) -> bool:
    import bisect
    if not intervals:
        return False
    starts = [s for s, _ in intervals]
    i = bisect.bisect_right(starts, pos) - 1
    return i >= 0 and intervals[i][0] <= pos < intervals[i][1]


# ─── TSS Enrichment (Signac-equivalent, Python) ────────────────────────────

def compute_tss_enrichment(
    fragments_path: Path,
    tss_bed: Path,
    *,
    half_window: int = 1000,
    flank_inner: int = 900,
    flank_outer: int = 1000,
) -> dict[str, dict[str, Any]]:
    """Per-barcode TSS enrichment, Signac convention.

    Definition:
      TSS_enrichment = mean_count_in_TSS_center_window
                       / max(0.2, mean_count_in_flank_windows)

    The center window is TSS ± `half_window // 10` (= ±100 bp default).
    Flanks are the 100-bp strips at the inner/outer edges of ± 1 kb.
    """
    center_half = max(half_window // 10, 50)  # default ±100 bp
    flank_lo, flank_hi = flank_inner, flank_outer  # default 900..1000 bp

    # Build per-chrom flat intervals (sorted)
    tss = load_tss(tss_bed)
    center_intv: dict[str, list[tuple[int, int]]] = {}
    flank_intv: dict[str, list[tuple[int, int]]] = {}
    for chrom, positions in tss.items():
        ci: list[tuple[int, int]] = []
        fi: list[tuple[int, int]] = []
        for p in positions:
            ci.append((max(0, p - center_half), p + center_half))
            fi.append((max(0, p - flank_hi), max(0, p - flank_lo)))
            fi.append((p + flank_lo, p + flank_hi))
        center_intv[chrom] = sorted(ci)
        flank_intv[chrom] = sorted(fi)

    per_bc: dict[str, dict[str, float]] = defaultdict(
        lambda: {"reads_center": 0, "reads_flank": 0}
    )
    center_width = 2 * center_half
    flank_width = 2 * (flank_hi - flank_lo)

    for chrom, start, end, bc, count in iter_fragments(fragments_path):
        midpoint = (start + end) // 2
        if chrom in center_intv and _bisect_hit(center_intv[chrom], midpoint):
            per_bc[bc]["reads_center"] += count
        if chrom in flank_intv and _bisect_hit(flank_intv[chrom], midpoint):
            per_bc[bc]["reads_flank"] += count

    result: dict[str, dict[str, Any]] = {}
    for bc, d in per_bc.items():
        center_rate = d["reads_center"] / max(center_width, 1)
        flank_rate = max(d["reads_flank"] / max(flank_width, 1),
                          0.2 / max(center_width, 1))
        result[bc] = {
            "reads_tss_center": int(d["reads_center"]),
            "reads_tss_flank": int(d["reads_flank"]),
            "tss_enrichment": float(center_rate / flank_rate) if flank_rate > 0 else 0.0,
        }
    return result


# ─── Nucleosome Signal (Signac-equivalent) ─────────────────────────────────

def compute_nucleosome_signal(fragments_path: Path) -> dict[str, dict[str, Any]]:
    """Per-barcode nucleosome signal.

    Definition: ratio of mono-nucleosome fragments (147 ≤ length ≤ 294 bp)
    to nucleosome-free (length < 147 bp).
    """
    per_bc: dict[str, dict[str, int]] = defaultdict(
        lambda: {"nfr": 0, "mono": 0, "di_plus": 0, "total": 0}
    )
    for _, start, end, bc, count in iter_fragments(fragments_path):
        length = end - start
        rec = per_bc[bc]
        rec["total"] += count
        if length < 147:
            rec["nfr"] += count
        elif length <= 294:
            rec["mono"] += count
        else:
            rec["di_plus"] += count
    result: dict[str, dict[str, Any]] = {}
    for bc, d in per_bc.items():
        denom = max(d["nfr"], 1)
        result[bc] = {
            "n_nfr": d["nfr"], "n_mono": d["mono"], "n_di_plus": d["di_plus"],
            "n_fragments": d["total"],
            "nucleosome_signal": float(d["mono"] / denom),
        }
    return result


# ─── FRIP ───────────────────────────────────────────────────────────────────

def compute_frip(fragments_path: Path, peaks_bed: Path
                 ) -> dict[str, dict[str, Any]]:
    """Per-barcode fraction of reads (fragments) in peaks."""
    peaks = load_peaks(peaks_bed)
    per_bc: dict[str, dict[str, int]] = defaultdict(
        lambda: {"reads_total": 0, "reads_in_peaks": 0}
    )
    for chrom, start, end, bc, count in iter_fragments(fragments_path):
        rec = per_bc[bc]
        rec["reads_total"] += count
        midpoint = (start + end) // 2
        if chrom in peaks and _bisect_hit(peaks[chrom], midpoint):
            rec["reads_in_peaks"] += count
    return {
        bc: {
            "reads_in_peaks": d["reads_in_peaks"],
            "frip": float(d["reads_in_peaks"] / d["reads_total"])
            if d["reads_total"] else 0.0,
        }
        for bc, d in per_bc.items()
    }


# ─── Commands ───────────────────────────────────────────────────────────────

def cmd_qc_atac(args: argparse.Namespace) -> int:
    """Compute per-barcode ATAC QC (fragments, TSS enrichment, nucleosome
    signal, FRIP) from a fragments BED + TSS + peaks."""
    pd = _require("pandas")
    setup_logging()
    frag = Path(args.fragments)
    if not frag.is_file():
        raise SystemExit(f"fragments not found: {frag}")
    tss_metrics = compute_tss_enrichment(frag, Path(args.tss_bed)) if args.tss_bed else {}
    nuc_metrics = compute_nucleosome_signal(frag)
    frip_metrics = compute_frip(frag, Path(args.peaks_bed)) if args.peaks_bed else {}

    all_bc = set(tss_metrics) | set(nuc_metrics) | set(frip_metrics)
    rows = []
    for bc in sorted(all_bc):
        row = {"barcode": bc}
        row.update(tss_metrics.get(bc, {}))
        row.update(nuc_metrics.get(bc, {}))
        row.update(frip_metrics.get(bc, {}))
        rows.append(row)
    df = pd.DataFrame(rows)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_atac_qc.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"ATAC QC table: {out}  ({len(df)} barcodes)")
    if "tss_enrichment" in df.columns and len(df):
        print(f"  median TSS enrichment: {df['tss_enrichment'].median():.2f}")
    if "nucleosome_signal" in df.columns and len(df):
        print(f"  median nucleosome signal: {df['nucleosome_signal'].median():.3f}")
    if "frip" in df.columns and len(df):
        print(f"  median FRIP: {df['frip'].median():.3f}")
    return 0


def cmd_joint_qc(args: argparse.Namespace) -> int:
    """Join ATAC + RNA per-barcode tables and apply Signac-convention
    thresholds.

    Defaults (Signac multiome pbmc vignette, 10x epi multiome QC tutorial):
        ATAC fragments:        1,800  ≤ count   ≤ 100,000
        RNA UMIs:              1,000  ≤ count   ≤ 25,000
        TSS enrichment:        > 1
        Nucleosome signal:     < 2
        FRIP:                  > 0.15
    """
    pd = _require("pandas")
    setup_logging()
    rna = pd.read_csv(args.rna_qc, sep=None, engine="python")
    atac = pd.read_csv(args.atac_qc, sep=None, engine="python")
    merged = rna.merge(atac, on="barcode", how="outer")

    # Apply thresholds, with sensible defaults for missing columns
    def safe_col(col, default):
        return merged[col] if col in merged.columns else pd.Series(default, index=merged.index)

    rna_pass = (
        (safe_col("umis", 0) >= args.min_umis)
        & (safe_col("umis", 0) <= args.max_umis)
        & (safe_col("genes", 0) >= args.min_genes)
        & (safe_col("pct_mt", 0) <= args.max_pct_mt)
    )
    atac_pass = (
        (safe_col("n_fragments", 0) >= args.min_frags)
        & (safe_col("n_fragments", 0) <= args.max_frags)
        & (safe_col("tss_enrichment", 0) >= args.min_tss)
        & (safe_col("nucleosome_signal", 0) <= args.max_nuc)
        & (safe_col("frip", 0) >= args.min_frip)
    )
    merged["rna_pass"] = rna_pass
    merged["atac_pass"] = atac_pass
    merged["qc"] = "neither"
    merged.loc[rna_pass & ~atac_pass, "qc"] = "RNA only"
    merged.loc[~rna_pass & atac_pass, "qc"] = "ATAC only"
    merged.loc[rna_pass & atac_pass, "qc"] = "both"

    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_joint_qc.tsv"
    merged.to_csv(out, sep="\t", index=False)
    counts = merged["qc"].value_counts().to_dict()
    summary = {
        "n_barcodes_total": int(len(merged)),
        "n_rna_pass": int(rna_pass.sum()),
        "n_atac_pass": int(atac_pass.sum()),
        "n_both_pass": int((rna_pass & atac_pass).sum()),
        "qc_counts": counts,
        "thresholds": {
            "rna": {"min_umis": args.min_umis, "max_umis": args.max_umis,
                     "min_genes": args.min_genes, "max_pct_mt": args.max_pct_mt},
            "atac": {"min_frags": args.min_frags, "max_frags": args.max_frags,
                      "min_tss": args.min_tss, "max_nuc": args.max_nuc,
                      "min_frip": args.min_frip},
        },
    }
    save_json(f"multiome_{args.label}_joint_qc", summary)
    print(f"Joint QC: {out}")
    for k, v in counts.items():
        print(f"  {k}: {v}")
    return 0


# ─── TF-IDF + LSI ───────────────────────────────────────────────────────────

def cmd_lsi(args: argparse.Namespace) -> int:
    """TF-IDF normalization + truncated SVD on a peak × cell matrix.

    Drops dim 1 (correlates with sequencing depth), following Signac
    `RunSVD` + `DepthCor` convention.
    """
    np = _require("numpy")
    import importlib
    sp = importlib.import_module("scipy.sparse")
    pd = _require("pandas")
    anndata = _require("anndata")
    from sklearn.decomposition import TruncatedSVD
    setup_logging()

    adata = anndata.read_h5ad(args.input)
    logging.info("Loaded ATAC AnnData: %d cells × %d peaks",
                  adata.n_obs, adata.n_vars)
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.astype("float64")

    # TF-IDF: tf = X / row_sums; idf = log(1 + n_cells / col_sums)
    row_sums = np.asarray(X.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1.0
    tf = sp.diags(1.0 / row_sums) @ X

    col_sums = np.asarray(X.sum(axis=0)).flatten()
    col_sums[col_sums == 0] = 1.0
    idf = np.log1p(X.shape[0] / col_sums)
    tf_idf = tf @ sp.diags(idf)

    # log-transform per Stuart-lab "log_tf" convention (Signac default)
    tf_idf.data = np.log1p(tf_idf.data * 1e4)

    n_comp = int(args.n_components) + 1  # +1 to drop the depth dim
    svd = TruncatedSVD(n_components=n_comp, random_state=int(args.seed),
                        algorithm="arpack")
    U = svd.fit_transform(tf_idf)
    # Drop component 0 (depth) per Signac convention; keep remaining
    keep = U[:, 1:]
    adata.obsm["X_lsi"] = keep
    var_explained = svd.explained_variance_ratio_

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_lsi.h5ad"
    adata.write_h5ad(out)
    print(f"LSI-embedded ATAC AnnData: {out}")
    print(f"  X_lsi shape: {keep.shape}")
    print(f"  variance explained (first 5 retained): "
          f"{[round(v, 3) for v in var_explained[1:6]]}")
    return 0


# ─── WNN joint embedding via muon ──────────────────────────────────────────

def cmd_wnn(args: argparse.Namespace) -> int:
    """Joint weighted-nearest-neighbor embedding + UMAP via muon.

    Equivalent to Seurat 5 `FindMultiModalNeighbors` + `RunUMAP`.
    """
    np = _require("numpy")
    anndata = _require("anndata")
    sc = _require("scanpy")
    mu = _require("muon")
    setup_logging()

    rna = anndata.read_h5ad(args.rna_h5ad)
    atac = anndata.read_h5ad(args.atac_h5ad)
    # Restrict to common barcodes
    common = sorted(set(rna.obs_names) & set(atac.obs_names))
    if len(common) < 50:
        raise SystemExit(
            f"Too few shared barcodes ({len(common)}); RNA and ATAC AnnData "
            f"must share the same per-cell barcodes for WNN."
        )
    rna = rna[common].copy()
    atac = atac[common].copy()

    # Ensure both modalities have a per-modality embedding
    if "X_pca" not in rna.obsm:
        logging.info("Running PCA on RNA modality (no X_pca found)")
        sc.pp.normalize_total(rna, target_sum=1e4)
        sc.pp.log1p(rna)
        sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3"
                                      if rna.X.max() > 50 else "seurat")
        rna = rna[:, rna.var.highly_variable].copy()
        sc.pp.scale(rna, max_value=10)
        sc.tl.pca(rna, n_comps=int(args.n_pca))
    if "X_lsi" not in atac.obsm:
        raise SystemExit(
            "ATAC AnnData has no `obsm['X_lsi']`. Run `multiome lsi` first "
            "to compute it."
        )

    mdata = mu.MuData({"rna": rna, "atac": atac})
    # Per-modality neighbors first
    sc.pp.neighbors(mdata["rna"], use_rep="X_pca", n_neighbors=int(args.n_neighbors))
    sc.pp.neighbors(mdata["atac"], use_rep="X_lsi", n_neighbors=int(args.n_neighbors))
    # Joint WNN
    mu.pp.neighbors(mdata, key_added="wnn")
    mu.tl.umap(mdata, neighbors_key="wnn", min_dist=float(args.min_dist),
                random_state=int(args.seed))
    # Optional joint clustering
    if args.cluster:
        from scanpy.tools._leiden import leiden
        leiden(mdata, resolution=float(args.resolution),
                neighbors_key="wnn", key_added="leiden_wnn")

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_wnn.h5mu"
    mdata.write(out)
    print(f"WNN joint embedding: {out}")
    print(f"  cells: {mdata.n_obs}")
    print(f"  modalities: {list(mdata.mod.keys())}")
    if args.cluster:
        try:
            n_clusters = mdata.obs["leiden_wnn"].nunique()
            print(f"  Leiden clusters: {n_clusters}")
        except Exception:
            pass
    return 0


# ─── Peak-to-gene linkage ───────────────────────────────────────────────────

def cmd_peak_to_gene(args: argparse.Namespace) -> int:
    """Per-peak correlation with nearby genes (Pearson by default).

    For each peak, finds all genes whose TSS is within `--window` bp of
    the peak center, computes the Pearson/Spearman correlation between
    that peak's accessibility (column of ATAC matrix) and the gene's
    expression (column of RNA matrix), and reports correlation + p-value.

    Sign of correlation indicates whether the peak is positively
    (enhancer-like) or negatively (repressor / chromatin closure)
    associated with gene expression.
    """
    np = _require("numpy")
    pd = _require("pandas")
    sp = _require("scipy.stats")
    anndata = _require("anndata")
    setup_logging()

    rna = anndata.read_h5ad(args.rna_h5ad)
    atac = anndata.read_h5ad(args.atac_h5ad)
    # Restrict to common cells
    common = sorted(set(rna.obs_names) & set(atac.obs_names))
    if len(common) < 30:
        raise SystemExit(
            f"Need ≥30 shared barcodes for correlation (got {len(common)})."
        )
    rna = rna[common].copy()
    atac = atac[common].copy()

    # Parse peak coordinates from var_names (chr:start-end) OR from
    # var columns chrom/start/end.
    def parse_peak(name):
        if ":" in name and "-" in name:
            chrom, rest = name.split(":", 1)
            start, end = rest.split("-")
            return chrom, int(start), int(end)
        return None

    peak_coords = []
    for v in atac.var_names:
        c = parse_peak(v)
        if c:
            peak_coords.append(c)
        else:
            peak_coords.append((None, None, None))

    # Parse gene TSS — assume var_names are gene symbols, and we have a
    # supplied TSS BED to resolve them.
    tss_bed = Path(args.tss_bed) if args.tss_bed else None
    gene_tss: dict[str, tuple[str, int]] = {}
    if tss_bed and tss_bed.is_file():
        with _open_text(tss_bed) as fh:
            for ln in fh:
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                try:
                    gene_tss[parts[3]] = (parts[0], int(parts[1]))
                except ValueError:
                    continue

    window = int(args.window)
    use_spearman = args.method == "spearman"

    pd_rna = pd.DataFrame(rna.X.toarray() if hasattr(rna.X, "toarray") else rna.X,
                           index=rna.obs_names, columns=rna.var_names)
    pd_atac = pd.DataFrame(atac.X.toarray() if hasattr(atac.X, "toarray") else atac.X,
                            index=atac.obs_names, columns=atac.var_names)

    rows = []
    n_pairs = 0
    for peak_idx, peak in enumerate(atac.var_names):
        coord = peak_coords[peak_idx]
        if coord[0] is None:
            continue
        chrom, p_start, p_end = coord
        peak_mid = (p_start + p_end) // 2
        for gene, (g_chrom, g_tss) in gene_tss.items():
            if g_chrom != chrom:
                continue
            if abs(g_tss - peak_mid) > window:
                continue
            if gene not in pd_rna.columns:
                continue
            x = pd_atac.iloc[:, peak_idx].to_numpy()
            y = pd_rna[gene].to_numpy()
            if x.std() == 0 or y.std() == 0:
                continue
            if use_spearman:
                rho, p = sp.spearmanr(x, y)
            else:
                rho, p = sp.pearsonr(x, y)
            rows.append({
                "peak": peak, "gene": gene,
                "distance": int(g_tss - peak_mid),
                "correlation": float(rho), "pvalue": float(p),
            })
            n_pairs += 1
            if args.max_pairs and n_pairs >= args.max_pairs:
                break
        if args.max_pairs and n_pairs >= args.max_pairs:
            break

    if not rows:
        raise SystemExit(
            "No peak-gene pairs found. Check that peaks have chr:start-end "
            "names and that --tss-bed maps gene symbols to coordinates."
        )

    df = pd.DataFrame(rows)
    from statsmodels.stats.multitest import multipletests
    _, padj, _, _ = multipletests(df["pvalue"].fillna(1.0),
                                    method="fdr_bh")
    df["padj"] = padj

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_peak2gene.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"Peak-to-gene links: {out}  ({len(df)} pairs)")
    print(f"  significant (padj < 0.05): {(df['padj'] < 0.05).sum()}")
    print(f"  median |correlation|: {df['correlation'].abs().median():.3f}")
    return 0


# ─── End-to-end showcase ────────────────────────────────────────────────────

def cmd_showcase(args: argparse.Namespace) -> int:
    """Run qc-atac + joint-qc + plot composite figure in one tool call."""
    setup_logging()
    ts = time.strftime("%Y%m%d_%H%M%S")
    label = args.label or "multiome_showcase"
    out_dir = REPORT_DIR / f"{ts}_{safe_label(label)}_showcase"
    plots_dir = out_dir / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: ATAC QC
    atac_args = argparse.Namespace(
        fragments=args.fragments,
        tss_bed=args.tss_bed,
        peaks_bed=args.peaks_bed,
        label=label,
    )
    cmd_qc_atac(atac_args)

    # Find the latest atac_qc output
    atac_tsv = sorted(REPORT_DIR.glob(f"*{label}_atac_qc.tsv"))[-1]

    # Step 2: joint QC, if RNA QC is supplied
    if args.rna_qc:
        joint_args = argparse.Namespace(
            rna_qc=args.rna_qc, atac_qc=str(atac_tsv), label=label,
            min_umis=1000, max_umis=25000, min_genes=200, max_pct_mt=0.20,
            min_frags=1800, max_frags=100000, min_tss=1.0, max_nuc=2.0,
            min_frip=0.15,
        )
        cmd_joint_qc(joint_args)

    # Composite figure with TSS / nucleosome / FRIP histograms
    pd = _require("pandas")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    df = pd.read_csv(atac_tsv, sep="\t")
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), facecolor="white")
    fig.suptitle(f"10x Multiome ATAC QC — {label}", fontsize=16,
                  fontweight="bold", y=0.995)
    panels = [
        ("n_fragments", "Fragments per cell (log10)", True),
        ("tss_enrichment", "TSS enrichment", False),
        ("nucleosome_signal", "Nucleosome signal", False),
        ("frip", "FRIP (fraction reads in peaks)", False),
        ("reads_tss_center", "Reads in TSS ±100 bp", True),
        ("reads_in_peaks", "Reads in peaks", True),
    ]
    import numpy as np
    for ax, (col, title, log_x) in zip(axes.flat, panels):
        if col not in df.columns:
            ax.axis("off"); continue
        vals = df[col].dropna().values
        if log_x:
            vals = np.log10(vals[vals > 0] + 1)
            xlabel = title
        else:
            xlabel = title
        ax.hist(vals, bins=40, color="#5C8DAA", alpha=0.85,
                 edgecolor="#1F2933", linewidth=0.4)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("barcodes")
        if col == "tss_enrichment":
            ax.axvline(1.0, color="#C13E3E", ls="--", lw=1.2,
                        label="Signac threshold")
            ax.legend(fontsize=8)
        if col == "nucleosome_signal":
            ax.axvline(2.0, color="#C13E3E", ls="--", lw=1.2,
                        label="Signac threshold")
            ax.legend(fontsize=8)
        if col == "frip":
            ax.axvline(0.15, color="#C13E3E", ls="--", lw=1.2,
                        label="Signac threshold")
            ax.legend(fontsize=8)
    plt.tight_layout()
    comp = plots_dir / "atac_qc_composite.png"
    fig.savefig(comp, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(comp.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Narrative report
    n_bc = len(df)
    med_tss = df["tss_enrichment"].median() if "tss_enrichment" in df.columns else None
    med_nuc = df["nucleosome_signal"].median() if "nucleosome_signal" in df.columns else None
    med_frip = df["frip"].median() if "frip" in df.columns else None
    report = out_dir / "showcase_report.md"
    lines = [
        f"# 10x Multiome QC Showcase — {label}", "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        f"Fragments file: `{args.fragments}`",
        f"Barcodes profiled: **{n_bc:,}**",
        "",
        "## Median per-cell QC",
        "",
        f"- TSS enrichment: **{med_tss:.2f}**" if med_tss is not None else "",
        f"- Nucleosome signal: **{med_nuc:.3f}**" if med_nuc is not None else "",
        f"- FRIP: **{med_frip:.3f}**" if med_frip is not None else "",
        "",
        "## Composite figure",
        "",
        f"![ATAC QC composite](Plots/atac_qc_composite.png)",
        "",
        "## How to read",
        "",
        "1. **TSS enrichment > 1** is the Signac/Seurat convention for "
        "passing cells (red dashed line). Cells below the line are "
        "likely low-quality / ambient.",
        "2. **Nucleosome signal < 2** is the upper bound; very high "
        "values indicate poor Tn5 digestion / over-fragmentation.",
        "3. **FRIP > 0.15** indicates the cell's fragments concentrate "
        "in the peak set (i.e., the cell has interpretable open "
        "chromatin).",
        "4. **Fragments per cell** should be bi/multi-modal — the "
        "highest-fragment peak is real cells; the lower peak is "
        "background.",
        "",
        f"Atac QC table: `{atac_tsv}`",
    ]
    report.write_text("\n".join(l for l in lines if l is not None))
    print(f"Showcase report: {report}")
    print(f"Composite figure: {comp}")
    return 0


# ─── Playbook ───────────────────────────────────────────────────────────────

# ════════════════════════════════════════════════════════════════════════════
# QUADBIO-ABSORBED EXTENSIONS
#
# Clean-room reimplementations of additional analytical methods covered in
# quadbio/scMultiome_analysis_python_vignette (Treutlein lab, ETH Zürich,
# no LICENSE — clean-room only). All math paraphrased from the published
# algorithm descriptions; no source from the vignette is copied.
# ════════════════════════════════════════════════════════════════════════════


def cmd_da_peaks(args: argparse.Namespace) -> int:
    """Differential accessibility on TF-IDF normalized peaks via Wilcoxon.

    Pairs naturally with `multiome peak2gene` — peak2gene gives you
    enhancer→gene candidates; da-peaks tells you WHICH peaks are
    differentially open between groups.
    """
    sc = _require("scanpy")
    mu = _require("muon")
    anndata = _require("anndata")
    pd = _require("pandas")
    setup_logging()
    adata = anndata.read_h5ad(args.input)
    if args.cluster_key not in adata.obs.columns:
        raise SystemExit(
            f"--cluster-key {args.cluster_key!r} not found in adata.obs. "
            f"Available: {list(adata.obs.columns)[:20]}"
        )
    logging.info("Running mu.atac.pp.tfidf …")
    mu.atac.pp.tfidf(adata, scale_factor=1e4)
    logging.info("Wilcoxon rank-genes per %s …", args.cluster_key)
    sc.tl.rank_genes_groups(adata, groupby=args.cluster_key,
                              method="wilcoxon",
                              use_raw=False, layer=None)
    rgg = adata.uns["rank_genes_groups"]
    groups = list(rgg["names"].dtype.names)
    rows = []
    for g in groups:
        for i in range(min(args.top_n, len(rgg["names"][g]))):
            rows.append({
                "group": g,
                "peak": rgg["names"][g][i],
                "logfc": float(rgg["logfoldchanges"][g][i]),
                "score": float(rgg["scores"][g][i]),
                "pvalue": float(rgg["pvals"][g][i]),
                "padj":   float(rgg["pvals_adj"][g][i]),
            })
    df = pd.DataFrame(rows)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_da_peaks.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"Differential-accessibility table: {out}")
    print(f"  groups: {len(groups)} | top peaks per group: {args.top_n} | "
          f"rows: {len(df):,}")
    return 0


def cmd_atac_spectral(args: argparse.Namespace) -> int:
    """Jaccard-similarity Laplacian spectral embedding on an ATAC matrix.

    Alternative to TF-IDF + LSI. The snapATAC2-style approach builds a
    Jaccard similarity graph between cells over binarized peaks,
    computes the symmetric normalized Laplacian, and takes the top
    eigenvectors. Robust to depth differences without an explicit
    "drop dim 1" step.

    Implementation:
        peaks: cell × peak sparse matrix (binarized at >0)
        S    = Jaccard(cells) = (X @ X.T) / (deg + deg.T - X @ X.T)
        L    = I - D^(-1/2) S D^(-1/2)
        emb  = top-k eigenvectors of (I - L) (= D^(-1/2) S D^(-1/2))
    """
    np = _require("numpy")
    import importlib
    sp = importlib.import_module("scipy.sparse")
    anndata = _require("anndata")
    setup_logging()
    adata = anndata.read_h5ad(args.input)
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    # Binarize
    Xb = (X > 0).astype("float64")
    # Sub-sample to keep things tractable on a laptop
    n_cells = Xb.shape[0]
    if args.max_cells and n_cells > args.max_cells:
        rng = np.random.default_rng(int(args.seed))
        idx = rng.choice(n_cells, size=args.max_cells, replace=False)
        Xb = Xb[idx]
        cell_idx = adata.obs_names[idx]
    else:
        cell_idx = adata.obs_names

    # cell × cell intersection
    inter = (Xb @ Xb.T).toarray()
    deg = inter.diagonal()
    union = deg[:, None] + deg[None, :] - inter
    union[union == 0] = 1.0
    jacc = inter / union
    np.fill_diagonal(jacc, 0.0)

    # k-nearest-neighbour sparsification
    k = int(args.n_neighbors)
    if k > 0:
        rows = []
        cols = []
        vals = []
        for i in range(jacc.shape[0]):
            top = np.argpartition(-jacc[i], min(k, jacc.shape[0] - 1))[:k]
            for j in top:
                if jacc[i, j] > 0:
                    rows.append(i); cols.append(j); vals.append(jacc[i, j])
        S = sp.csr_matrix((vals, (rows, cols)), shape=jacc.shape)
        # Symmetrize
        S = 0.5 * (S + S.T)
    else:
        S = sp.csr_matrix(jacc)

    d = np.asarray(S.sum(axis=1)).flatten()
    d[d == 0] = 1.0
    D_inv_sqrt = sp.diags(1.0 / np.sqrt(d))
    L_norm = D_inv_sqrt @ S @ D_inv_sqrt

    from scipy.sparse.linalg import eigsh
    n_comp = int(args.n_components)
    try:
        vals_eig, vecs = eigsh(L_norm, k=n_comp + 1, which="LM")
        # Sort descending, drop the trivial top eigenvector
        order = np.argsort(-vals_eig)
        vecs = vecs[:, order][:, 1:n_comp + 1]
    except Exception as exc:
        raise SystemExit(f"eigsh failed: {exc}")

    # Stitch back into AnnData (handle sub-sample case)
    if args.max_cells and n_cells > args.max_cells:
        full_emb = np.zeros((adata.n_obs, n_comp))
        cell_set = {bc: i for i, bc in enumerate(cell_idx)}
        for j, bc in enumerate(adata.obs_names):
            if bc in cell_set:
                full_emb[j] = vecs[cell_set[bc]]
        adata.obsm["X_spectral"] = full_emb
    else:
        adata.obsm["X_spectral"] = vecs

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_spectral.h5ad"
    adata.write_h5ad(out)
    print(f"Spectral-embedded ATAC AnnData: {out}")
    print(f"  X_spectral shape: {adata.obsm['X_spectral'].shape}")
    return 0


def cmd_chromvar(args: argparse.Namespace) -> int:
    """Clean-room chromVAR-style TF motif activity per cell.

    For each (cell, motif M):
       raw_dev   = (sum accessibility over peaks with motif M) - expected
       expected  = fraction-of-peaks-with-M × cell's total accessibility
       z_score   = (raw_dev - mean(raw_dev over K GC-matched background
                                    motif sets)) / std(raw_dev background)

    Inputs:
      --input        peak × cell h5ad (binarized inside)
      --motif-hits   peak × motif TSV (0/1 binary; columns are motif IDs)
      --gc-content   optional per-peak GC content TSV (gc_content column);
                     if absent, computed as 0.5 (no GC bias correction)

    Output: cell × motif z-score matrix as obsm['X_motif_zscore'] + a
    standalone TSV.
    """
    np = _require("numpy")
    import importlib
    sp = importlib.import_module("scipy.sparse")
    pd = _require("pandas")
    anndata = _require("anndata")
    setup_logging()
    adata = anndata.read_h5ad(args.input)
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    Xb = (X > 0).astype("float64")
    n_cells, n_peaks = Xb.shape

    # Motif hits: peak × motif (binary)
    hits = pd.read_csv(args.motif_hits, sep="\t", index_col=0)
    common = sorted(set(hits.index) & set(adata.var_names))
    if len(common) < 0.5 * len(adata.var_names):
        logging.warning("Only %d / %d peaks have motif annotations",
                          len(common), len(adata.var_names))
    hits = hits.loc[common]
    # Align peaks
    peak_to_idx = {p: i for i, p in enumerate(adata.var_names)}
    peak_idx = [peak_to_idx[p] for p in common]
    Xb_aligned = Xb[:, peak_idx]
    H = sp.csr_matrix(hits.values.astype("float64"))

    # GC content
    if args.gc_content and Path(args.gc_content).is_file():
        gc = pd.read_csv(args.gc_content, sep="\t", index_col=0)
        gc = gc.loc[common, "gc_content"].fillna(0.5).values
    else:
        gc = np.full(len(common), 0.5)
    log_mean_acc = np.log1p(np.asarray(Xb_aligned.sum(axis=0)).flatten())

    # Raw deviations: cell × motif
    cell_total = np.asarray(Xb_aligned.sum(axis=1)).flatten()
    peak_total = np.asarray(Xb_aligned.sum(axis=0)).flatten()
    grand_total = peak_total.sum() or 1.0
    # observed: how often motif M appears in cell c's accessible peaks
    obs = Xb_aligned @ H  # cell × motif
    if sp.issparse(obs):
        obs = obs.toarray()
    # expected: cell_total × (motif total / grand total)
    motif_freq = np.asarray(H.sum(axis=0)).flatten() / grand_total
    expected = cell_total[:, None] * motif_freq[None, :] * len(common)
    raw_dev = obs - expected

    # Background: K=50 GC-matched motif sets per motif
    K = int(args.k_background)
    rng = np.random.default_rng(int(args.seed))
    # Bin peaks by (GC quantile, log_mean_acc quantile) for matched sampling
    gc_bin = pd.qcut(gc, q=10, labels=False, duplicates="drop")
    acc_bin = pd.qcut(log_mean_acc, q=10, labels=False, duplicates="drop")
    bin_key = gc_bin * 100 + acc_bin
    bin_idx = pd.Series(np.arange(len(common))).groupby(bin_key).apply(list).to_dict()
    z_scores = np.zeros_like(raw_dev)
    n_motifs = H.shape[1]
    logging.info("Computing chromVAR deviations: %d cells × %d motifs × K=%d backgrounds",
                  n_cells, n_motifs, K)
    for m in range(n_motifs):
        motif_peaks = np.where(H[:, m].toarray().flatten() > 0)[0]
        if len(motif_peaks) == 0:
            continue
        bg_devs = np.zeros((K, n_cells))
        for k in range(K):
            sampled = []
            for p in motif_peaks:
                b = bin_key[p] if p < len(bin_key) else 0
                pool = bin_idx.get(b, motif_peaks.tolist())
                sampled.append(int(rng.choice(pool)))
            sampled = np.array(sampled)
            H_bg = sp.csr_matrix(
                (np.ones(len(sampled)), (sampled, np.zeros(len(sampled), dtype=int))),
                shape=(len(common), 1),
            )
            obs_bg = Xb_aligned @ H_bg
            obs_bg = obs_bg.toarray().flatten() if sp.issparse(obs_bg) else obs_bg.flatten()
            motif_freq_bg = len(sampled) / grand_total
            expected_bg = cell_total * motif_freq_bg * len(common)
            bg_devs[k] = obs_bg - expected_bg
        bg_mean = bg_devs.mean(axis=0)
        bg_std = bg_devs.std(axis=0)
        bg_std[bg_std == 0] = 1.0
        z_scores[:, m] = (raw_dev[:, m] - bg_mean) / bg_std

    adata.obsm["X_motif_zscore"] = z_scores
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out_h5 = REPORT_DIR / f"{ts}_{safe_label(args.label)}_chromvar.h5ad"
    out_tsv = REPORT_DIR / f"{ts}_{safe_label(args.label)}_motif_zscore.tsv"
    adata.write_h5ad(out_h5)
    z_df = pd.DataFrame(z_scores, index=adata.obs_names, columns=hits.columns)
    z_df.to_csv(out_tsv, sep="\t")
    print(f"chromVAR motif z-scores: {out_h5}")
    print(f"  per-cell TSV: {out_tsv}")
    print(f"  shape: cells={n_cells}, motifs={n_motifs}")
    return 0


def cmd_css(args: argparse.Namespace) -> int:
    """Cluster Similarity Spectrum batch correction (He 2020).

    Per-batch HVG → per-batch Leiden clustering → cluster centroids
    (in PCA space) → represent each cell as its vector of correlations
    to all batch×cluster centroids. The correlation matrix IS the
    corrected embedding.

    Apache-2.0 / BSD-3 alternative to GPL-licensed Harmony.
    """
    np = _require("numpy")
    pd = _require("pandas")
    sc = _require("scanpy")
    anndata = _require("anndata")
    setup_logging()
    adata = anndata.read_h5ad(args.input)
    if args.batch_key not in adata.obs.columns:
        raise SystemExit(
            f"--batch-key {args.batch_key!r} not in adata.obs. "
            f"Available: {list(adata.obs.columns)[:20]}"
        )
    if "X_pca" not in adata.obsm:
        logging.info("Running PCA first …")
        if adata.X.max() > 50:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=int(args.n_hvg),
                                      flavor="seurat")
        adata = adata[:, adata.var.highly_variable].copy()
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=int(args.n_pca))

    batches = adata.obs[args.batch_key].unique()
    logging.info("CSS: %d batches", len(batches))
    centroids = []
    centroid_labels = []
    for b in batches:
        idx = adata.obs[args.batch_key] == b
        sub = adata[idx].copy()
        # Per-batch Leiden
        sc.pp.neighbors(sub, n_neighbors=int(args.n_neighbors),
                          use_rep="X_pca")
        sc.tl.leiden(sub, resolution=float(args.resolution),
                       key_added="leiden_css")
        # Cluster centroids in PCA space
        for c in sub.obs["leiden_css"].unique():
            mask = sub.obs["leiden_css"] == c
            centroid = sub.obsm["X_pca"][mask].mean(axis=0)
            centroids.append(centroid)
            centroid_labels.append(f"{b}_{c}")
    centroids = np.array(centroids)  # (n_centroids, n_pca)

    # For each cell: pearson correlation to all centroids
    pca = adata.obsm["X_pca"]
    pca_centered = pca - pca.mean(axis=1, keepdims=True)
    cent_centered = centroids - centroids.mean(axis=1, keepdims=True)
    num = pca_centered @ cent_centered.T
    pca_norm = np.linalg.norm(pca_centered, axis=1, keepdims=True)
    cent_norm = np.linalg.norm(cent_centered, axis=1, keepdims=True)
    css_emb = num / (pca_norm @ cent_norm.T + 1e-9)
    adata.obsm["X_css"] = css_emb
    adata.uns["css_centroid_labels"] = centroid_labels

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_css.h5ad"
    adata.write_h5ad(out)
    print(f"CSS batch-corrected AnnData: {out}")
    print(f"  X_css shape: {css_emb.shape}  ({len(centroid_labels)} centroids)")
    return 0


def cmd_multivi(args: argparse.Namespace) -> int:
    """MultiVI deep joint VAE via scvi-tools (Ashuach 2023).

    scvi-tools is BSD-3 but a heavy install (~hundreds MB of torch).
    This command requires `pip install scvi-tools`; we don't pin it as
    a hard dependency.
    """
    try:
        import scvi
    except ImportError:
        raise SystemExit(
            "scvi-tools not installed. Install with: "
            "pip install scvi-tools"
        )
    anndata = _require("anndata")
    mu = _require("muon")
    setup_logging()
    rna = anndata.read_h5ad(args.rna_h5ad)
    atac = anndata.read_h5ad(args.atac_h5ad)
    common = sorted(set(rna.obs_names) & set(atac.obs_names))
    if len(common) < 50:
        raise SystemExit(f"Too few shared barcodes ({len(common)}).")
    rna = rna[common].copy()
    atac = atac[common].copy()

    # Concatenate features for MultiVI
    rna.var["modality"] = "Gene Expression"
    atac.var["modality"] = "Peaks"
    mdata = mu.MuData({"rna": rna, "atac": atac})
    scvi.model.MULTIVI.setup_mudata(
        mdata,
        modalities={"rna_layer": "rna", "protein_layer": "atac"},
        batch_key=args.batch_key if args.batch_key in mdata.obs.columns else None,
    )
    mvi = scvi.model.MULTIVI(
        mdata,
        n_genes=(rna.var.modality == "Gene Expression").sum(),
        n_regions=(atac.var.modality == "Peaks").sum(),
    )
    logging.info("Training MultiVI for %d epochs", args.epochs)
    mvi.train(max_epochs=int(args.epochs))

    mdata.obsm["X_multivi"] = mvi.get_latent_representation()
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_multivi.h5mu"
    mdata.write(out)
    print(f"MultiVI joint embedding: {out}")
    print(f"  X_multivi shape: {mdata.obsm['X_multivi'].shape}")
    return 0



def write_playbook() -> Path:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "MULTIOME_10X_ANALYZE.md"
    path.write_text("""# Skill: 10x Multiome Analysis

Clean-room Python implementations of the analytical methods covered in
the 10x Genomics
[`analysis_guides`](https://github.com/10XGenomics/analysis_guides) repo
(Epi multiome QC tutorial + interactive QC visualization), plus the
Stuart-lab Signac vignette extensions (TF-IDF/LSI, WNN, peak-to-gene
linkage).

No source from the 10x repo (which has no LICENSE) is copied; the math
is paraphrased from the published algorithm descriptions.

## Commands

```bash
# 1. Per-barcode ATAC QC from a fragments BED
igvfagent multiome qc-atac \\
    --fragments fragments.tsv.gz \\
    --tss-bed gencode_tss.bed \\
    --peaks-bed atac_peaks.bed \\
    --label run1

# 2. Joint QC: merge ATAC + RNA per-barcode tables, apply thresholds
igvfagent multiome joint-qc \\
    --rna-qc run1_rna_qc.tsv \\
    --atac-qc run1_atac_qc.tsv \\
    --label run1

# 3. TF-IDF + truncated SVD on a peak × cell matrix (drops depth dim 1)
igvfagent multiome lsi --input atac_peaks.h5ad --label run1

# 4. Joint WNN embedding + UMAP via muon
igvfagent multiome wnn \\
    --rna-h5ad rna.h5ad --atac-h5ad atac_with_lsi.h5ad \\
    --cluster --label run1

# 5. Peak-to-gene correlation (Pearson over cells, ±500 kb window)
igvfagent multiome peak2gene \\
    --rna-h5ad rna.h5ad --atac-h5ad atac.h5ad \\
    --tss-bed gencode_tss.bed --window 500000 \\
    --method pearson --label run1

# 6. Single-command showcase
igvfagent multiome showcase --fragments fragments.tsv.gz \\
    --tss-bed tss.bed --peaks-bed peaks.bed \\
    --rna-qc rna_qc.tsv --label run1
```

## QC thresholds (Signac convention)

| Metric | Pass criterion |
|---|---|
| RNA UMIs | 1,000 ≤ count ≤ 25,000 |
| RNA genes | ≥ 200 |
| RNA % mito | ≤ 20% |
| ATAC fragments | 1,800 ≤ count ≤ 100,000 |
| TSS enrichment | > 1 |
| Nucleosome signal | < 2 |
| FRIP | > 0.15 |

All thresholds are overridable via CLI flags.

## References

- 10x Genomics `analysis_guides` (algorithm spec only — no source copied):
  https://github.com/10XGenomics/analysis_guides
- Stuart-lab Signac (MIT) — `TSSEnrichment`, `NucleosomeSignal`,
  `RunSVD`, `FindMultiModalNeighbors`
- Hao et al., Cell 2021 — Seurat 5 WNN method
""")
    return path


# ─── CLI ────────────────────────────────────────────────────────────────────

def add_subparsers(sub) -> None:
    """Add multiome analysis subparsers to an existing argparse `sub`.

    Called from Scripts/multiome_10x_pipeline.py.
    """
    p = sub.add_parser("qc-atac",
        help="Per-barcode ATAC QC (TSS enrichment, nucleosome signal, FRIP).")
    p.add_argument("--fragments", required=True)
    p.add_argument("--tss-bed", default=None)
    p.add_argument("--peaks-bed", default=None)
    p.add_argument("--label", default="multiome_atac_qc")
    p.set_defaults(func=cmd_qc_atac)

    p = sub.add_parser("joint-qc",
        help="Merge ATAC + RNA per-barcode tables and apply joint thresholds.")
    p.add_argument("--rna-qc", required=True)
    p.add_argument("--atac-qc", required=True)
    p.add_argument("--min-umis", type=int, default=1000)
    p.add_argument("--max-umis", type=int, default=25000)
    p.add_argument("--min-genes", type=int, default=200)
    p.add_argument("--max-pct-mt", type=float, default=0.20)
    p.add_argument("--min-frags", type=int, default=1800)
    p.add_argument("--max-frags", type=int, default=100000)
    p.add_argument("--min-tss", type=float, default=1.0)
    p.add_argument("--max-nuc", type=float, default=2.0)
    p.add_argument("--min-frip", type=float, default=0.15)
    p.add_argument("--label", default="multiome_joint_qc")
    p.set_defaults(func=cmd_joint_qc)

    p = sub.add_parser("lsi",
        help="TF-IDF + truncated SVD on a peak × cell matrix.")
    p.add_argument("--input", required=True,
                    help="h5ad with peak × cell counts.")
    p.add_argument("--n-components", type=int, default=50)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--label", default="multiome_lsi")
    p.set_defaults(func=cmd_lsi)

    p = sub.add_parser("wnn",
        help="Joint WNN embedding via muon.")
    p.add_argument("--rna-h5ad", required=True)
    p.add_argument("--atac-h5ad", required=True)
    p.add_argument("--n-pca", type=int, default=50)
    p.add_argument("--n-neighbors", type=int, default=20)
    p.add_argument("--min-dist", type=float, default=0.3)
    p.add_argument("--cluster", action="store_true")
    p.add_argument("--resolution", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--label", default="multiome_wnn")
    p.set_defaults(func=cmd_wnn)

    p = sub.add_parser("peak2gene",
        help="Per-peak correlation with nearby genes.")
    p.add_argument("--rna-h5ad", required=True)
    p.add_argument("--atac-h5ad", required=True)
    p.add_argument("--tss-bed", required=True,
                    help="BED with chrom, start, end, gene_name.")
    p.add_argument("--window", type=int, default=500000,
                    help="bp window around each gene's TSS (default 500 kb).")
    p.add_argument("--method", default="pearson",
                    choices=["pearson", "spearman"])
    p.add_argument("--max-pairs", type=int, default=200000)
    p.add_argument("--label", default="multiome_peak2gene")
    p.set_defaults(func=cmd_peak_to_gene)


    # ── Quadbio extensions ────────────────────────────────────────────
    p = sub.add_parser("da-peaks",
        help="Differential accessibility on TF-IDF peaks via Wilcoxon.")
    p.add_argument("--input", required=True,
                    help="ATAC h5ad with peaks (will be TF-IDF normalized).")
    p.add_argument("--cluster-key", default="leiden_wnn",
                    help="adata.obs column to compare across.")
    p.add_argument("--top-n", type=int, default=50)
    p.add_argument("--label", default="multiome_da_peaks")
    p.set_defaults(func=cmd_da_peaks)

    p = sub.add_parser("atac-spectral",
        help="Jaccard-Laplacian spectral embedding alternative to LSI.")
    p.add_argument("--input", required=True)
    p.add_argument("--n-components", type=int, default=30)
    p.add_argument("--n-neighbors", type=int, default=20)
    p.add_argument("--max-cells", type=int, default=5000)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--label", default="multiome_spectral")
    p.set_defaults(func=cmd_atac_spectral)

    p = sub.add_parser("chromvar",
        help="Clean-room chromVAR-style TF motif activity per cell.")
    p.add_argument("--input", required=True,
                    help="ATAC h5ad with peaks.")
    p.add_argument("--motif-hits", required=True,
                    help="Peak × motif binary TSV.")
    p.add_argument("--gc-content", default=None,
                    help="Optional per-peak GC content TSV.")
    p.add_argument("--k-background", type=int, default=50)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--label", default="multiome_chromvar")
    p.set_defaults(func=cmd_chromvar)

    p = sub.add_parser("css",
        help="Cluster Similarity Spectrum batch correction (He 2020).")
    p.add_argument("--input", required=True)
    p.add_argument("--batch-key", required=True,
                    help="adata.obs column with batch IDs.")
    p.add_argument("--n-hvg", type=int, default=2000)
    p.add_argument("--n-pca", type=int, default=50)
    p.add_argument("--n-neighbors", type=int, default=20)
    p.add_argument("--resolution", type=float, default=1.0)
    p.add_argument("--label", default="multiome_css")
    p.set_defaults(func=cmd_css)

    p = sub.add_parser("multivi",
        help="MultiVI deep joint VAE (scvi-tools, optional dep).")
    p.add_argument("--rna-h5ad", required=True)
    p.add_argument("--atac-h5ad", required=True)
    p.add_argument("--batch-key", default="batch")
    p.add_argument("--epochs", type=int, default=50)
    p.add_argument("--label", default="multiome_multivi")
    p.set_defaults(func=cmd_multivi)

    p = sub.add_parser("showcase",
        help="One-command ATAC QC + joint QC + 6-panel composite figure.")
    p.add_argument("--fragments", required=True)
    p.add_argument("--tss-bed", default=None)
    p.add_argument("--peaks-bed", default=None)
    p.add_argument("--rna-qc", default=None,
                    help="Optional RNA QC TSV to enable joint QC + cell calling.")
    p.add_argument("--label", default="multiome")
    p.set_defaults(func=cmd_showcase)

    p = sub.add_parser("write-analyze-playbook",
        help="Write Docs/Skills/MULTIOME_10X_ANALYZE.md")
    p.set_defaults(func=lambda a: (print(f"Wrote: {write_playbook()}"), 0)[1])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="multiome_10x_analyze",
        description="Clean-room 10x Multiome analytics (QC + LSI + WNN + peak2gene).",
    )
    sub = parser.add_subparsers(dest="cmd")
    add_subparsers(sub)
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help(); return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
