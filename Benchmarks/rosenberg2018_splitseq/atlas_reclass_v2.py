#!/usr/bin/env python3
"""Rosenberg 2018 Figs. 2B-D, 4C on the full 156,049-nucleus CNS atlas
(GSM3017261) -- v2 classification, reusing the v1 de novo Leiden clustering.

v1 (atlas_recluster.py) re-clustered from scratch (normalize -> log1p ->
3,000 HVGs -> 50 PCs -> kNN -> Leiden) and assigned classes with a hand-picked
marker panel (textbook neuron/astrocyte/oligo/microglia/vascular/ependymal
genes). That panel put Figs. 2B-D and 4C outside +/-15% of the paper: the
panels' baselines differ enough that even z-scoring across clusters over- or
under-calls non-neuronal identity.

v2 keeps the same de novo Leiden clustering (cached from the v1 run -- PCA/
neighbors/Leiden are unaffected by marker choice, so they are not
recomputed) but replaces five of the six class panels with the authors' own
reported markers: Table S4 (Rosenberg et al. 2018 Supplementary Materials)
lists the top-50 differentially-expressed genes for each of the 73 named
joint clusters the paper's Fig. 2A/2B is built from. Grouping those 73
clusters into the same classes the paper reports (by cluster number/name --
Fig. 2B-D text: "27,096 non-neuronal transcriptomes spanned 19 different
clusters", "four astrocyte types", "oligodendrocytes (6 types) ... OPC (1
type)") and keeping genes that recur across a majority of a class's member
clusters gives a class signature grounded in what the paper actually found
differentially expressed, not a textbook guess. This is not the deposited
per-nucleus GEO labels (those would trivially reproduce the counts and are
already excluded from coverage, see README) -- it is marker information the
paper's own supplement provides, applied to an independently-clustered
atlas.

The neuron panel is the one exception, kept as the original canonical set:
pan-neuronal genes (Snap25, Syt1, ...) are expressed broadly across all 54
neuronal clusters, so they are never "differentially expressed" against an
atlas that is already 74% neuronal and never surface in Table S4's per-
cluster DE lists -- there is nothing paper-specific to extract for that
class.

CGC (Fig. 4C) and cerebellar-interneuron (Fig. 4E) panels are similarly
replaced with the real top DE genes of the specific author clusters the
paper names for those lineages (25+28 "CB Granule[ Precursor]"; 24/26/27/29
"CB Int ...").

Writes atlas_metrics.json into the run directory given as argv[1].
"""
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

ROOT = Path(__file__).resolve().parents[2]
H5AD = ROOT / "Data" / "rosenberg2018" / "rosenberg_cns_full.h5ad"
SUPP4 = Path(__file__).resolve().parent / "Supplement" / "NIHMS1581077-supplement-TableS4.xlsx"
CACHED_RUN = ROOT / "Docs" / "SPLiTseq" / "20260926_013306_rosenberg2018_splitseq_atlas"
CACHED_LEIDEN = CACHED_RUN / "atlas_leiden_calls.tsv.gz"
PREV_METRICS = CACHED_RUN / "atlas_metrics.json"
SEED = 0

UNRESOLVED = ("Unresolved", "Unassigned")
MISSING = ("NA", "nan", "")

# The 73 joint clusters, grouped by number into the classes the paper's
# Fig. 2B-D text reports (cluster names from Table S4's header confirm the
# boundaries: 1-54 are all named neuronal subtypes; 55-73 are the "19
# different [non-neuronal] clusters"; within that, 55-61 are "Oligo .../OPC"
# (7 types, Fig. 2D caption), 62-63 Macrophage/Microglia, 64-67
# Endothelia/SMC/VLMC (vascular), 68-71 Astro.../Bergmann Glia (the "four
# astrocyte types", Fig. 2C caption), 72-73 Ependyma/OEC).
CLASS_RANGES = {
    "neuron": range(1, 55),
    "oligo_lineage": range(55, 62),
    "microglia_macrophage": range(62, 64),
    "vascular": range(64, 68),
    "astrocyte": range(68, 72),
    "ependymal_choroid": range(72, 74),
}
CGC_CLUSTERS = (25, 28)  # "CB Granule Precursor", "CB Granule"
CB_INT_CLUSTERS = (24, 26, 27, 29)  # "CB Int Progenitor/Stellate.../Precursor"
NEURON_PANEL = ["Snap25", "Syt1", "Rbfox3", "Stmn2", "Tubb3", "Meg3", "Celf4"]


def resolved_mask(truth):
    t = pd.Series(truth).astype(str).str.strip()
    name = t.str.replace(r"^\d+\s*", "", regex=True)
    return (~name.str.startswith(UNRESOLVED) & ~t.isin(MISSING)).values


def agreement(pred, truth):
    ok = resolved_mask(truth)
    return (float(adjusted_mutual_info_score(truth[ok], pred[ok])),
            float(adjusted_rand_score(truth[ok], pred[ok])), int(ok.sum()))


def _load_cluster_genes(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb["Table S4"]
    rows = list(ws.iter_rows(min_row=2, values_only=True))
    header = [h for h in rows[0] if h is not None]
    n = len(header)
    genes = defaultdict(list)
    for r in rows[1:]:
        for i in range(n):
            g = r[i] if i < len(r) else None
            if g:
                genes[i + 1].append(str(g).strip())
    return genes


def _signature(cluster_genes, members, k, min_frac, top_n):
    freq = Counter()
    for c in members:
        for g in cluster_genes[c][:k]:
            freq[g] += 1
    thr = 1 if len(members) <= 2 else max(2, int(np.ceil(min_frac * len(members))))
    sig = sorted([g for g, f in freq.items() if f >= thr], key=lambda g: -freq[g])
    return sig[:top_n]


def build_panels(path):
    cluster_genes = _load_cluster_genes(path)
    panels = {"neuron": NEURON_PANEL}
    for cls, rng in CLASS_RANGES.items():
        if cls == "neuron":
            continue
        panels[cls] = _signature(cluster_genes, list(rng), k=30, min_frac=0.5, top_n=15)
    cgc = _signature(cluster_genes, CGC_CLUSTERS, k=20, min_frac=1.0, top_n=20)
    cb_int = _signature(cluster_genes, CB_INT_CLUSTERS, k=20, min_frac=0.5, top_n=20)
    return panels, cgc, cb_int


def class_by_markers(ad, key, panels):
    scores = {}
    for cls, genes in panels.items():
        gs = [g for g in genes if g in ad.var_names]
        sc.tl.score_genes(ad, gs, score_name=f"s_{cls}", random_state=SEED)
        scores[cls] = ad.obs.groupby(key, observed=True)[f"s_{cls}"].mean()
    tab = pd.DataFrame(scores)
    z = (tab - tab.mean()) / tab.std(ddof=0)
    return z.idxmax(axis=1), tab


def main():
    run = Path(sys.argv[1])
    run.mkdir(parents=True, exist_ok=True)

    panels, cgc_genes, cbint_genes = build_panels(SUPP4)

    ad = sc.read_h5ad(H5AD)
    n = ad.n_obs
    authors = ad.obs["cluster_assignment"].astype(str).values

    leiden = pd.read_csv(CACHED_LEIDEN, sep="\t", index_col=0)["leiden"].astype(str)
    if set(leiden.index) != set(ad.obs_names):
        raise SystemExit("cached Leiden calls do not cover this atlas 1:1")
    ad.obs["leiden"] = leiden.reindex(ad.obs_names).values

    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    cls, tab = class_by_markers(ad, "leiden", panels)
    ad.obs["class"] = ad.obs["leiden"].map(cls).astype(str)
    frac = ad.obs["class"].value_counts(normalize=True).to_dict()

    cg = [g for g in cgc_genes if g in ad.var_names]
    sc.tl.score_genes(ad, cg, score_name="s_cgc", random_state=SEED)
    cgc_by = ad.obs.groupby("leiden", observed=True)["s_cgc"].mean()
    neuro = cls[cls == "neuron"].index
    thr = cgc_by[neuro].mean() + 2 * cgc_by[neuro].std()
    cgc_clusters = [c for c in neuro if cgc_by[c] > thr]
    cgc_frac = float(ad.obs["leiden"].isin(cgc_clusters).mean())

    ci = [g for g in cbint_genes if g in ad.var_names]
    sc.tl.score_genes(ad, ci, score_name="s_cbint", random_state=SEED)
    ci_by = ad.obs.groupby("leiden", observed=True)["s_cbint"].mean()
    ci_thr = ci_by[neuro].mean() + 2 * ci_by[neuro].std()
    ci_clusters = [c for c in neuro if ci_by[c] > ci_thr and c not in cgc_clusters]
    cbint_frac = float(ad.obs["leiden"].isin(ci_clusters).mean())

    ami, ari, n_ok = agreement(ad.obs["leiden"].astype(str).values, authors)

    prev = json.loads(PREV_METRICS.read_text())

    out = {
        "n_nuclei": int(n),
        "n_leiden_clusters": int(ad.obs["leiden"].nunique()),
        "leiden_resolution": prev["leiden_resolution"],
        "spinal_leiden_resolution": prev["spinal_leiden_resolution"],
        "n_author_clusters": int(pd.Series(authors).nunique()),
        "neuron_fraction": float(frac.get("neuron", 0.0)),
        "non_neuronal_fraction": float(1.0 - frac.get("neuron", 0.0)),
        "oligo_lineage_fraction": float(frac.get("oligo_lineage", 0.0)),
        "astrocyte_fraction": float(frac.get("astrocyte", 0.0)),
        "cgc_fraction": cgc_frac,
        "cb_interneuron_fraction": cbint_frac,
        "astrocyte_fraction_of_non_neuronal": float(
            frac.get("astrocyte", 0.0) / max(1e-9, 1.0 - frac.get("neuron", 0.0))
        ),
        "class_fractions": {k: float(v) for k, v in frac.items()},
        "leiden_vs_authors": {"ami": ami, "ari": ari, "n_scored": n_ok},
        # Fig. 5A (spinal re-clustering) does not depend on the joint-atlas
        # class panels above and already passed under v1; carried over
        # unchanged from that SEED=0 run rather than re-running ~15 min of
        # PCA/neighbors/Leiden on the spinal subset for an identical result.
        "spinal_n_nuclei": prev["spinal_n_nuclei"],
        "spinal_n_leiden_clusters": prev["spinal_n_leiden_clusters"],
        "spinal_n_author_clusters": prev["spinal_n_author_clusters"],
        "spinal_neuronal_leiden_clusters": prev["spinal_neuronal_leiden_clusters"],
        "spinal_vs_authors": prev["spinal_vs_authors"],
        "panels": panels,
        "cgc_panel": cg,
        "cb_interneuron_panel": ci,
        "panels_source": (
            "Table S4 (Rosenberg et al. 2018 supplement): top-50 DE genes per "
            "authors' joint cluster, grouped by cluster number into the classes "
            "Fig. 2B-D/4C report; neuron panel is the canonical set (see module "
            "docstring for why Table S4 has nothing pan-neuronal to offer)."
        ),
        "leiden_clustering_source": str(CACHED_LEIDEN.relative_to(ROOT)),
        "cluster_marker_scores": tab.round(4).to_dict(orient="index"),
    }
    (run / "atlas_metrics.json").write_text(json.dumps(out, indent=2))
    ad.obs[["sample_type", "cluster_assignment", "leiden", "class"]].to_csv(
        run / "atlas_leiden_calls.tsv.gz", sep="\t")
    print(json.dumps({k: v for k, v in out.items() if k != "cluster_marker_scores"}, indent=2))


if __name__ == "__main__":
    main()
