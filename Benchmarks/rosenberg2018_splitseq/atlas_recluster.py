#!/usr/bin/env python3
"""Rosenberg 2018 Figs. 2, 4, 5 on the full 156,049-nucleus CNS atlas (GSM3017261).

De novo re-clustering: the authors' clustering code is in the Science
supplement, not in a public repository, so this is an independent standard
pipeline (normalize -> log1p -> 3,000 HVGs -> scale -> 50 PCs -> kNN ->
Leiden), with clusters assigned to cell classes by established marker panels
(the paper's "expression of established markers"), never by the authors'
labels. The deposited labels (``cluster_assignment``,
``spinal_cluster_assignment``) are used only afterwards, to score agreement.

Writes atlas_metrics.json into the run directory given as argv[1].
Heavy: run through sbatch (~40-60 GB RAM for the full atlas).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

ROOT = Path(__file__).resolve().parents[2]
H5AD = ROOT / "Data" / "rosenberg2018" / "rosenberg_cns_full.h5ad"
SEED = 0

PANELS = {
    "neuron": ["Snap25", "Syt1", "Rbfox3", "Stmn2", "Tubb3", "Meg3", "Celf4"],
    "astrocyte": ["Aqp4", "Slc1a3", "Gja1", "Aldh1l1", "Slc1a2", "Gfap"],
    "oligo_lineage": ["Olig1", "Olig2", "Sox10", "Pdgfra", "Plp1", "Mbp", "Mog"],
    "microglia_macrophage": ["Cx3cr1", "P2ry12", "C1qa", "Csf1r", "Mrc1"],
    "vascular": ["Cldn5", "Flt1", "Pecam1", "Vtn", "Pdgfrb", "Rgs5"],
    "ependymal_choroid": ["Foxj1", "Ttr", "Kl", "Ccdc153"],
}
CGC = ["Gabra6", "Pax6", "Cbln3", "Barhl1", "Atoh1", "Neurod1", "Zic1"]
CB_INT = ["Pax2", "Tfap2b", "Lhx5", "Lhx1", "Sorcs3", "Kit"]
UNRESOLVED = ("Unresolved", "Unassigned")
MISSING = ("NA", "nan", "")


def resolved_mask(truth):
    t = pd.Series(truth).astype(str).str.strip()
    name = t.str.replace(r"^\d+\s*", "", regex=True)
    return (~name.str.startswith(UNRESOLVED) & ~t.isin(MISSING)).values


def class_by_markers(ad, key):
    scores = {}
    for cls, genes in PANELS.items():
        gs = [g for g in genes if g in ad.var_names]
        sc.tl.score_genes(ad, gs, score_name=f"s_{cls}", random_state=SEED)
        scores[cls] = ad.obs.groupby(key, observed=True)[f"s_{cls}"].mean()
    tab = pd.DataFrame(scores)
    # Panels differ in baseline (the neuron panel holds very abundant genes such
    # as Meg3/Snap25), so each panel's cluster means are standardised across
    # clusters before the argmax; a raw argmax calls microglia "neuron".
    z = (tab - tab.mean()) / tab.std(ddof=0)
    return z.idxmax(axis=1), tab


RES_GRID = (1.0, 2.0, 3.0, 4.0, 5.0)


def cluster(ad, res, target=None):
    """Leiden at ``res``; with ``target``, the smallest grid resolution giving
    >= target clusters (the paper's stated granularity, not its labels)."""
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.raw = ad
    sc.pp.highly_variable_genes(ad, n_top_genes=3000, flavor="seurat")
    sub = ad[:, ad.var.highly_variable].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=50, random_state=SEED)
    sc.pp.neighbors(sub, n_neighbors=30, n_pcs=50, random_state=SEED)
    for r in ((res,) if target is None else RES_GRID):
        sc.tl.leiden(sub, resolution=r, random_state=SEED, key_added="leiden",
                     flavor="igraph", n_iterations=2)
        print(f"leiden resolution {r}: {sub.obs['leiden'].nunique()} clusters", flush=True)
        if target is None or sub.obs["leiden"].nunique() >= target:
            break
    ad.obs["leiden"] = sub.obs["leiden"].values
    ad.uns["leiden_resolution"] = r
    return ad


def agreement(pred, truth):
    """AMI/ARI on nuclei the authors resolved (their Unresolved/Unassigned/NA excluded)."""
    ok = resolved_mask(truth)
    return (float(adjusted_mutual_info_score(truth[ok], pred[ok])),
            float(adjusted_rand_score(truth[ok], pred[ok])), int(ok.sum()))


def main():
    run = Path(sys.argv[1])
    run.mkdir(parents=True, exist_ok=True)
    ad = sc.read_h5ad(H5AD)
    n = ad.n_obs
    authors = ad.obs["cluster_assignment"].astype(str).values
    ad = cluster(ad, 1.0, target=73)
    cls, tab = class_by_markers(ad, "leiden")
    ad.obs["class"] = ad.obs["leiden"].map(cls).astype(str)
    frac = ad.obs["class"].value_counts(normalize=True).to_dict()
    # CGCs: neuronal clusters whose CGC-panel score is the top neuronal-subtype signal
    cg = [g for g in CGC if g in ad.var_names]
    sc.tl.score_genes(ad, cg, score_name="s_cgc", random_state=SEED)
    cgc_by = ad.obs.groupby("leiden", observed=True)["s_cgc"].mean()
    neuro = cls[cls == "neuron"].index
    thr = cgc_by[neuro].mean() + 2 * cgc_by[neuro].std()
    cgc_clusters = [c for c in neuro if cgc_by[c] > thr]
    cgc_frac = float(ad.obs["leiden"].isin(cgc_clusters).mean())
    ci = [g for g in CB_INT if g in ad.var_names]
    sc.tl.score_genes(ad, ci, score_name="s_cbint", random_state=SEED)
    ci_by = ad.obs.groupby("leiden", observed=True)["s_cbint"].mean()
    ci_thr = ci_by[neuro].mean() + 2 * ci_by[neuro].std()
    ci_clusters = [c for c in neuro if ci_by[c] > ci_thr and c not in cgc_clusters]
    cbint_frac = float(ad.obs["leiden"].isin(ci_clusters).mean())
    ami, ari, n_ok = agreement(ad.obs["leiden"].astype(str).values, authors)

    # Fig. 5: spinal cord nuclei re-clustered on their own
    sp = ad.obs["sample_type"].astype(str).str.contains("spine").values
    sad = sc.read_h5ad(H5AD)[sp].copy()
    sad = cluster(sad, 1.0, target=44)
    scls, _ = class_by_markers(sad, "leiden")
    s_truth = sad.obs["spinal_cluster_assignment"].astype(str).values
    s_ami, s_ari, s_n = agreement(sad.obs["leiden"].astype(str).values, s_truth)
    s_named = pd.Series(s_truth).astype(str)[resolved_mask(s_truth)]

    out = {
        "n_nuclei": int(n),
        "n_leiden_clusters": int(ad.obs["leiden"].nunique()),
        "leiden_resolution": float(ad.uns["leiden_resolution"]),
        "spinal_leiden_resolution": float(sad.uns["leiden_resolution"]),
        "n_author_clusters": int(pd.Series(authors).nunique()),
        "neuron_fraction": float(frac.get("neuron", 0.0)),
        "non_neuronal_fraction": float(1.0 - frac.get("neuron", 0.0)),
        "oligo_lineage_fraction": float(frac.get("oligo_lineage", 0.0)),
        "astrocyte_fraction": float(frac.get("astrocyte", 0.0)),
        "cgc_fraction": cgc_frac,
        "cb_interneuron_fraction": cbint_frac,
        "astrocyte_fraction_of_non_neuronal": float(frac.get("astrocyte", 0.0) / max(1e-9, 1.0 - frac.get("neuron", 0.0))),
        "class_fractions": {k: float(v) for k, v in frac.items()},
        "leiden_vs_authors": {"ami": ami, "ari": ari, "n_scored": n_ok},
        "spinal_n_nuclei": int(sp.sum()),
        "spinal_n_leiden_clusters": int(sad.obs["leiden"].nunique()),
        "spinal_n_author_clusters": int(s_named.nunique()),
        "spinal_neuronal_leiden_clusters": int((scls == "neuron").sum()),
        "spinal_vs_authors": {"ami": s_ami, "ari": s_ari, "n_scored": s_n},
        "cluster_marker_scores": tab.round(4).to_dict(orient="index"),
    }
    (run / "atlas_metrics.json").write_text(json.dumps(out, indent=2))
    ad.obs[["sample_type", "cluster_assignment", "leiden", "class"]].to_csv(
        run / "atlas_leiden_calls.tsv.gz", sep="\t")
    print(json.dumps({k: v for k, v in out.items() if k != "cluster_marker_scores"}, indent=2))


if __name__ == "__main__":
    main()
