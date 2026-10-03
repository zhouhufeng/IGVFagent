#!/usr/bin/env python3
"""snMCT-seq: the methylation half as well as the RNA half.

snMCT-seq measures methylcytosine AND the transcriptome from the same
nucleus, and the two halves are published as separate matrices. Handing the
dataset to the generic single-cell route analyses whichever matrix happens
to sort first and silently drops the other -- which on IGVFDS4826YNLK means
reporting an RNA clustering and never mentioning the methylation that is the
point of the assay. 933 MeasurementSets on the Portal are snMCT-seq, the
third-largest assay type, so that is a lot of half-answers.

The reason this needs its own script rather than a second call to
`sc-analyze` is that the modalities are not the same kind of number.
Measured on IGVFDS4826YNLK's published matrices:

  RNA          int64 counts, 92.3% zeros, max 7144   -> log-normalise
  methylation  float32, 0.009-2.91, mean 1.01, no zeros, no NaN

The methylation values are already a normalised mC ratio centred on 1, not
counts. Running the count pipeline over them -- library-size scaling then
log1p -- would be normalising a ratio a second time and calling the result
expression. So methylation goes to PCA and clustering directly, which is
what the ALLCools-style snmC workflow does.

With 670 of 704 cells shared between the halves, the clusterings can be
cross-tabulated, which is the question the assay exists to answer: does
methylation state agree with transcriptional state in the same nucleus.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import raw_data_pipeline as rp                                # noqa: E402

ROOT = rp.ROOT
OUT_DIR = ROOT / "Docs" / "MCT"
WORK_DIR = ROOT / "Data" / "MCT"
LOG_DIR = ROOT / "Docs" / "Logs"

# content_type -> modality. These strings are what the Portal actually uses
# on snMCT-seq analysis sets; matching them explicitly is what keeps the
# methylation half from being chosen or dropped by matrix size.
RNA_TYPES = ("cell by gene matrix",)
MC_BIN_TYPES = ("cell by bin methylation count matrix",)
MC_POS_TYPES = ("cell by position methylation count matrix",)
QC_TYPES = ("per-cell quality report",)


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"mct_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(log),
                                  logging.StreamHandler(sys.stdout)])
    return log


def _scanpy():
    import anndata as ad
    import numpy as np
    import scanpy as sc
    return sc, ad, np


def find_modalities(accession: str) -> dict:
    """Locate each modality's file across the set and its analysis sets."""
    out: "dict[str, list[dict]]" = {"rna": [], "mc_bin": [], "mc_pos": [],
                                     "qc": [], "other": []}
    seen: "set[str]" = set()
    candidates = [accession] + [d.get("accession")
                                 for d in rp.derived_analysis_sets(accession)]
    for acc in [c for c in candidates if c]:
        for f in rp.list_files(acc):
            fid = f.get("accession")
            if not fid or fid in seen:
                continue
            seen.add(fid)
            ct = str(f.get("content_type") or "").strip().lower()
            rec = {**f, "from_set": acc}
            if ct in RNA_TYPES:
                out["rna"].append(rec)
            elif ct in MC_BIN_TYPES:
                out["mc_bin"].append(rec)
            elif ct in MC_POS_TYPES:
                out["mc_pos"].append(rec)
            elif ct in QC_TYPES:
                out["qc"].append(rec)
            else:
                out["other"].append(rec)
    out["analysis_sets"] = [c for c in candidates[1:]]        # type: ignore
    return out


def _fetch(rec: dict, dest_dir: Path) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    name = Path(str(rec.get("href") or rec.get("accession"))).name
    dest = dest_dir / name
    if dest.exists():
        print(f"Cached      {rec.get('accession')} ({rp.gb(rec.get('file_size'))} GB)")
    else:
        print(f"Downloading {rec.get('accession')} ({rp.gb(rec.get('file_size'))} GB)…")
        rp.portal_download(rec["href"], dest)
    return dest


def analyse_methylation(path: Path, out: Path, resolution: float,
                         n_pcs: int) -> dict:
    """Cluster cells on the normalised mC ratio matrix.

    No library-size scaling and no log1p: the values are already a ratio
    centred on 1, and treating them as counts would normalise twice. Bins
    are selected by variance across cells -- the methylation equivalent of
    highly-variable genes -- then PCA and Leiden, which is the shape of the
    standard snmC workflow.
    """
    sc, ad, np = _scanpy()
    a = ad.read_h5ad(path)
    summary: "dict[str, Any]" = {"file": path.name,
                                  "n_cells": int(a.n_obs),
                                  "n_bins": int(a.n_vars)}
    X = np.asarray(a.X.todense()) if hasattr(a.X, "todense") else np.asarray(a.X)
    X = np.nan_to_num(X.astype("float32"), nan=1.0)
    summary["value_min"] = float(X.min())
    summary["value_max"] = float(X.max())
    summary["value_mean"] = float(X.mean())

    # Per-cell mean mC ratio: the simplest global summary, and the one that
    # separates cell types in snmC data before any clustering.
    a.obs["mc_ratio_mean"] = X.mean(axis=1)
    var = X.var(axis=0)
    keep = np.argsort(var)[::-1][: min(5000, X.shape[1])]
    summary["bins_used"] = int(len(keep))
    a2 = ad.AnnData(X[:, keep], obs=a.obs.copy())
    sc.pp.scale(a2, max_value=10)
    n_pcs = int(min(n_pcs, min(a2.shape) - 1))
    sc.tl.pca(a2, n_comps=n_pcs)
    sc.pp.neighbors(a2, n_pcs=n_pcs)
    sc.tl.umap(a2)
    sc.tl.leiden(a2, resolution=resolution, key_added="mc_leiden",
                  flavor="igraph", n_iterations=2, directed=False)
    summary["n_pcs"] = n_pcs
    summary["n_clusters"] = int(a2.obs["mc_leiden"].nunique())
    summary["cells_per_cluster"] = {str(k): int(v) for k, v in
                                     a2.obs["mc_leiden"].value_counts().items()}
    out.mkdir(parents=True, exist_ok=True)
    a2.write_h5ad(out / "methylation_processed.h5ad", compression="gzip")
    _plot_methylation(a2, out)
    return summary, a2


def _plot_methylation(a2, out: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return
    um = a2.obsm.get("X_umap")
    if um is None:
        return
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    lab = a2.obs["mc_leiden"].astype(str).values
    for i, c in enumerate(sorted(set(lab))):
        sel = lab == c
        ax[0].scatter(um[sel, 0], um[sel, 1], s=10, label=c, alpha=.8)
    ax[0].set_title("Methylation UMAP, Leiden clusters")
    ax[0].set_xlabel("UMAP1"); ax[0].set_ylabel("UMAP2")
    ax[0].legend(fontsize=7, markerscale=1.5, ncol=2)
    sc2 = ax[1].scatter(um[:, 0], um[:, 1], s=10,
                        c=a2.obs["mc_ratio_mean"].values, cmap="viridis")
    ax[1].set_title("mean mC ratio per cell")
    ax[1].set_xlabel("UMAP1")
    fig.colorbar(sc2, ax=ax[1], label="mC ratio")
    fig.tight_layout()
    fig.savefig(out / "methylation_umap.png", dpi=140)
    plt.close(fig)


def compare_modalities(rna_h5ad: Path, mc_obj, out: Path) -> "Optional[dict]":
    """Cross-tabulate RNA and methylation clusters on the shared cells.

    This is the question snMCT-seq is for: the two halves come from the SAME
    nucleus, so agreement between them is measurable rather than assumed.
    """
    sc, ad, np = _scanpy()
    try:
        r = ad.read_h5ad(rna_h5ad)
    except OSError:
        return None
    key = next((k for k in ("leiden", "clusters", "cluster")
                if k in r.obs.columns), None)
    if key is None:
        return None
    shared = [c for c in mc_obj.obs_names if c in set(r.obs_names)]
    if len(shared) < 20:
        return {"shared_cells": len(shared),
                "note": "too few shared cells to compare"}
    rc = r.obs.loc[shared, key].astype(str)
    mc = mc_obj.obs.loc[shared, "mc_leiden"].astype(str)
    tab: "dict[str, Counter]" = {}
    for a_, b_ in zip(mc.values, rc.values):
        tab.setdefault(a_, Counter())[b_] += 1
    # Adjusted Rand Index: agreement corrected for chance, so "the
    # clusterings look similar" becomes a number.
    try:
        from sklearn.metrics import adjusted_rand_score
        ari = float(adjusted_rand_score(mc.values, rc.values))
    except ImportError:
        ari = float("nan")
    with (out / "modality_crosstab.tsv").open("w", newline="") as fh:
        rna_labels = sorted({x for c in tab.values() for x in c})
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mc_cluster"] + [f"rna_{x}" for x in rna_labels])
        for m_, counts in sorted(tab.items()):
            w.writerow([m_] + [counts.get(x, 0) for x in rna_labels])
    return {"shared_cells": len(shared), "rna_cluster_key": key,
            "adjusted_rand_index": ari,
            "n_mc_clusters": int(mc.nunique()),
            "n_rna_clusters": int(rc.nunique())}


# Covariates whose name marks them as a methylation LEVEL rather than a
# sequencing artefact. Global mCG fraction is the dominant biological axis
# in snmC data, so a clustering that tracks it is working; one that tracks
# read counts or mapping rate is tracking the library prep.
_BIOLOGICAL_HINTS = ("mcg", "mch", "mccc", "frac", "rate", "meth")
_TECHNICAL_HINTS = ("nreads", "readpairs", "mapped", "trim", "q20", "gc_perc",
                     "dup", "unmapped", "len", "index", "plate", "batch",
                     "run", "col", "time")


def confound_check(qc_path: Path, mc_obj, top: int = 6) -> "Optional[dict]":
    """Which per-cell covariates the methylation clustering explains.

    A clustering has to be interrogated, not just reported: mC clusters that
    separate cells by read depth are an artefact of the library, while ones
    that separate by global mCG fraction are the real dominant axis of snmC
    biology. Measured on IGVFDS4826YNLK the clusters explained 91.6% of
    DNA_mCGfrac and 13% of an RNA mapping rate -- which is the good case,
    and worth showing rather than asserting.

    Reports variance explained per covariate, labelled biological or
    technical by name, so a bad case is visible in the output instead of
    needing someone to think of the question.
    """
    _, _, np = _scanpy()
    try:
        text = rp.read_text_maybe_gzip(qc_path)
    except OSError:
        return None
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if not rows:
        return None
    idcol = next((c for c in ("wellprefix", "cell", "barcode", "cell_id")
                   if c in rows[0]), None)
    if idcol is None:
        return None
    qc = {r[idcol]: r for r in rows}
    labels = mc_obj.obs["mc_leiden"].astype(str).values
    names = list(mc_obj.obs_names)
    matched = sum(1 for n in names if n in qc)
    if matched < 50:
        return {"matched_cells": matched,
                "note": "too few cells matched the QC report to test"}
    scored = []
    for col in rows[0].keys():
        if col == idcol:
            continue
        vals, labs = [], []
        for n, lab in zip(names, labels):
            v = qc.get(n, {}).get(col, "")
            if v in ("", "NA", "nan"):
                continue
            try:
                vals.append(float(v))
                labs.append(lab)
            except ValueError:
                continue
        if len(vals) < 50:
            continue
        v = np.asarray(vals, dtype="float64")
        l = np.asarray(labs)
        groups = [v[l == g] for g in sorted(set(labs)) if (l == g).sum() > 3]
        if len(groups) < 2:
            continue
        total = float(((v - v.mean()) ** 2).sum())
        if total <= 0:
            continue
        gm = np.array([g.mean() for g in groups])
        ns = np.array([len(g) for g in groups])
        between = float(((gm - v.mean()) ** 2 * ns).sum())
        low = col.lower()
        kind = ("biological" if any(h in low for h in _BIOLOGICAL_HINTS)
                and not any(h in low for h in _TECHNICAL_HINTS)
                else "technical")
        scored.append({"covariate": col, "variance_explained": between / total,
                        "kind": kind})
    scored.sort(key=lambda r: -r["variance_explained"])
    top_tech = next((r for r in scored if r["kind"] == "technical"), None)
    top_bio = next((r for r in scored if r["kind"] == "biological"), None)
    verdict = "inconclusive"
    if top_bio and (not top_tech or
                    top_bio["variance_explained"] > top_tech["variance_explained"]):
        verdict = (f"clusters track {top_bio['covariate']} "
                   f"({top_bio['variance_explained']:.0%}) — a methylation "
                   f"level, i.e. the expected biological axis")
    elif top_tech:
        verdict = (f"WARNING: clusters track {top_tech['covariate']} "
                   f"({top_tech['variance_explained']:.0%}), a technical "
                   f"covariate — treat them as library artefact until shown "
                   f"otherwise")
    return {"matched_cells": matched, "verdict": verdict,
            "top_covariates": scored[:top]}


def cmd_discover(args: argparse.Namespace) -> int:
    mod = find_modalities(args.accession)
    print(f"Accession:      {args.accession}")
    print(f"AnalysisSets:   {', '.join(mod['analysis_sets']) or 'none'}")
    for key, label in (("rna", "RNA (cell x gene)"),
                        ("mc_bin", "methylation (cell x bin)"),
                        ("mc_pos", "methylation (cell x position)"),
                        ("qc", "per-cell QC report")):
        files = mod[key]
        if files:
            for f in files:
                print(f"  {label:32} {f.get('accession')}  "
                      f"{rp.gb(f.get('file_size'))} GB")
        else:
            print(f"  {label:32} MISSING")
    if mod["other"]:
        print(f"  other files: {len(mod['other'])} "
              f"({', '.join(sorted({str(f.get('content_type')) for f in mod['other']}))})")
    if not mod["rna"] or not mod["mc_bin"]:
        print("\nBoth halves are needed for a joint analysis; one is absent, "
              "so `analyze` will do what it can and say what it skipped.")
    return 0


def cmd_analyze(args: argparse.Namespace) -> int:
    mod = find_modalities(args.accession)
    if not mod["rna"] and not mod["mc_bin"]:
        print(f"Neither an RNA nor a methylation matrix is published for "
              f"{args.accession}. Raw reads would need bisulfite-aware "
              f"alignment, which IGVFagent does not do.")
        return 2
    label = args.label or f"{time.strftime('%Y%m%d_%H%M%S')}_{rp.safe_label(args.accession)}"
    out = OUT_DIR / label
    out.mkdir(parents=True, exist_ok=True)
    work = WORK_DIR / label
    summary: "dict[str, Any]" = {"accession": args.accession,
                                  "analysis_sets": mod["analysis_sets"]}

    # --- RNA half, via the existing single-cell pipeline -----------------
    rna_run = None
    if mod["rna"]:
        rna_path = _fetch(mod["rna"][0], work)
        summary["rna_file"] = mod["rna"][0].get("accession")
        if not args.skip_rna:
            import subprocess
            cmd = [rp._igvfagent(), "sc-analyze", "pipeline",
                   "--input", str(rna_path), "--label", f"{label}_rna"]
            print("RNA: " + " ".join(cmd))
            subprocess.run(cmd, check=False)
            runs = sorted((ROOT / "Docs" / "SingleCell").glob(f"*{label}_rna"),
                           key=lambda p: p.name)
            rna_run = runs[-1] if runs else None
            summary["rna_output"] = str(rna_run) if rna_run else None
    else:
        print("RNA: no cell-by-gene matrix published; skipping that half.")

    # --- methylation half ------------------------------------------------
    mc_obj = None
    if mod["mc_bin"]:
        mc_path = _fetch(mod["mc_bin"][0], work)
        summary["mc_file"] = mod["mc_bin"][0].get("accession")
        print("Methylation: clustering on the mC ratio matrix "
              "(no count normalisation — the values are already ratios)")
        mc_summary, mc_obj = analyse_methylation(
            mc_path, out, args.resolution, args.n_pcs)
        summary["methylation"] = mc_summary
        for k, v in mc_summary.items():
            if k != "cells_per_cluster":
                print(f"  {k}: {v}")
    else:
        print("Methylation: no cell-by-bin matrix published; skipping.")

    # --- is the methylation clustering biology or library prep? ----------
    if mc_obj is not None and mod["qc"]:
        qc_path = _fetch(mod["qc"][0], work)
        conf = confound_check(qc_path, mc_obj)
        if conf:
            summary["confound_check"] = conf
            print("\nCluster interrogation (per-cell QC report):")
            print(f"  {conf.get('verdict')}")
            for r in conf.get("top_covariates", [])[:5]:
                print(f"    {r['variance_explained']:6.1%}  {r['covariate']} "
                      f"[{r['kind']}]")

    # --- joint ------------------------------------------------------------
    if mc_obj is not None and rna_run is not None:
        cmp_ = compare_modalities(rna_run / "processed.h5ad", mc_obj, out)
        if cmp_:
            summary["modality_comparison"] = cmp_
            print("\nJoint (same nuclei, both halves):")
            for k, v in cmp_.items():
                print(f"  {k}: {v}")
    elif mc_obj is not None:
        print("\nNo joint comparison: the RNA half was skipped or produced "
              "no clustered output.")

    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\nOutput: {out}")
    for p in sorted(out.glob("*.png")):
        print(f"Plot: {p}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="mct_analysis",
        description="snMCT-seq: RNA and methylation from the same nuclei.")
    sub = p.add_subparsers(dest="command", required=True)
    d = sub.add_parser("discover", help="Show which modalities are published.")
    d.add_argument("accession")
    a = sub.add_parser("analyze", help="Analyse both halves and compare them.")
    a.add_argument("accession")
    a.add_argument("--resolution", type=float, default=1.0)
    a.add_argument("--n-pcs", type=int, default=30)
    a.add_argument("--skip-rna", action="store_true",
                    help="Methylation only.")
    a.add_argument("--label")
    return p


def main(argv: "Optional[list[str]]" = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging()
    for d in (OUT_DIR, WORK_DIR):
        d.mkdir(parents=True, exist_ok=True)
    return {"discover": cmd_discover, "analyze": cmd_analyze}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
