#!/usr/bin/env python3
"""Rosenberg 2018 Fig. 3A-B: regional specificity of neuronal clusters via
Allen Developing Mouse Brain Atlas ISH composites.

Paper: for each cluster, average Allen DMBA ISH over its five most enriched
genes, using the P4 atlas for P2-dominated clusters and P14 for P11-dominated
ones; the composites "confirmed the high regional specificity of most types".

  de     (heavy, sbatch)  top-5 enriched genes per deposited cluster: largest
                          log fold change (cluster vs rest, Wilcoxon adj. p < 0.05,
                          detected in >= 10% of the cluster's nuclei)
  allen  (light, network) per gene: DMBA sagittal ISH at P4 / P14 ->
                          StructureUnionize expression_energy, rolled up to
                          coarse regions; composite = mean over the 5 genes;
                          the cluster's peak region is compared with the
                          region its deposited name states.

Writes fig3_metrics.json into the run directory (argv[2]).
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "Data" / "rosenberg2018"
ALLEN = D / "allen"
API = "https://api.brain-map.org/api/v2/data/query.json?criteria="

# coarse region -> developmental-ontology acronyms whose subtree it covers
# DMBA StructureUnionize rows exist only down to coarse structures (e.g. r1A,
# MPall, SPall; never CbH or OB), so each region is defined on structures that
# carry data: the cerebellum is the alar plate of rhombomere 1 (its anlage), and
# the olfactory bulb has no unionized rows at all and is not scored.
REGIONS = {  # most specific first: a structure counts toward the first region on its ancestry
    "hippocampus": ["MPall"],
    "cortex": ["Pall"],
    "striatum": ["SPall"],
    "hypothalamus_preoptic": ["THy", "PHy", "POTel"],
    "thalamus": ["p2"],
    "diencephalon_other": ["p1", "p3"],
    "midbrain": ["M"],
    "cerebellum": ["r1A"],
    "pons_hindbrain": ["is", "r1", "r2", "PH", "PMH"],
    "medulla": ["MH"],
    "spinal_cord": ["SpC"],
}
# the region prefix of a deposited cluster name (Rosenberg's naming); ambiguous
# prefixes (SC, Migrating Int, SVZ, Cajal-Retzius, glia, Unresolved) are not scored
NAME_RULES = [
    (r"^\d+ (HIPP|SUB) ", "hippocampus"),
    (r"^\d+ (CTX|CLAU) ", "cortex"),
    (r"^\d+ Medium Spiny", "striatum"),
    (r"^\d+ THAL ", "thalamus"),
    (r"^\d+ (MTt |Nigral)", "midbrain"),
    (r"^\d+ (CB |Purkinje)", "cerebellum"),
    (r"^\d+ MD ", "medulla"),
]


def claimed_region(name):
    for pat, reg in NAME_RULES:
        if re.search(pat, name):
            return reg
    return None


def de(run):
    import scanpy as sc
    ad = sc.read_h5ad(D / "rosenberg_cns_full.h5ad")
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.obs["cluster_assignment"] = ad.obs["cluster_assignment"].astype("category")
    sc.tl.rank_genes_groups(ad, "cluster_assignment", method="wilcoxon", pts=True)
    rows = []
    age = ad.obs.groupby("cluster_assignment", observed=True)["sample_type"].agg(
        lambda s: (s.astype(str).str.startswith("p2")).mean())
    for cl in ad.obs["cluster_assignment"].cat.categories:
        # "most highly enriched": largest log fold change among significant
        # genes detected in >= 10% of the cluster's nuclei
        t = sc.get.rank_genes_groups_df(ad, group=cl)
        t["pct_nz_group"] = t["names"].map(ad.uns["rank_genes_groups"]["pts"][cl])
        t = t[(t["pvals_adj"] < 0.05) & (t["pct_nz_group"] >= 0.10)
              & ~t["names"].str.startswith(("mt-", "Rpl", "Rps", "Gm"))]
        genes = t.sort_values("logfoldchanges", ascending=False)["names"].head(5).tolist()
        rows.append({"cluster": cl, "p2_fraction": float(age[cl]), "top5": ",".join(genes),
                     "n": int((ad.obs["cluster_assignment"] == cl).sum())})
    pd.DataFrame(rows).to_csv(run / "cluster_top5.tsv", sep="\t", index=False)


def q(crit):
    for _ in range(4):
        try:
            u = API + urllib.parse.quote(crit, safe=":,[]$='()")
            return json.load(urllib.request.urlopen(u, timeout=60))["msg"]
        except Exception:  # noqa: BLE001 - retry the public API
            time.sleep(3)
    return []


def descendants():
    g = json.load(open(ALLEN / "dev_structure_graph.json"))["msg"][0]
    idx = {}

    def walk(n, anc):
        anc = anc + [n["acronym"]]
        idx[n["id"]] = anc
        for c in n.get("children", []):
            walk(c, anc)
    walk(g, [])
    return idx


def gene_region_energy(gene, age, ancestry, cache):
    key = f"{gene}|{age}"
    if key in cache:
        return cache[key]
    ds = q("model::SectionDataSet,rma::criteria,[failed$eqfalse],products[id$eq3],"
           f"genes[acronym$eq'{gene}'],specimen(donor(age[name$eq'{age}'])),"
           "plane_of_section[name$eq'sagittal']")
    out = None
    if ds:
        su = q(f"model::StructureUnionize,rma::criteria,[section_data_set_id$eq{ds[0]['id']}],"
               "rma::options[num_rows$eqall]")
        e = {}
        for r in su:
            anc = ancestry.get(r["structure_id"], [])
            if any(a.startswith("v_") or a == "tracts" for a in anc):
                continue
            for reg, roots in REGIONS.items():
                if any(a in roots for a in anc):
                    e.setdefault(reg, []).append(r["expression_energy"])
                    break
        out = {reg: float(np.mean(v)) for reg, v in e.items()}
    cache[key] = out
    return out


def allen(run):
    top = pd.read_csv(run / "cluster_top5.tsv", sep="\t")
    ancestry = descendants()
    cpath = ALLEN / "gene_region_energy_cache.json"
    cache = json.loads(cpath.read_text()) if cpath.exists() else {}
    rows = []
    for _, r in top.iterrows():
        claim = claimed_region(r["cluster"])
        if claim is None:
            continue
        age = "P4" if r["p2_fraction"] >= 0.5 else "P14"
        maps = [gene_region_energy(g, age, ancestry, cache) for g in r["top5"].split(",")]
        maps = [m for m in maps if m]
        cpath.write_text(json.dumps(cache))
        if not maps:
            rows.append({"cluster": r["cluster"], "claimed": claim, "age": age, "n_genes_with_ish": 0})
            continue
        comp = pd.DataFrame(maps).mean()
        rows.append({"cluster": r["cluster"], "claimed": claim, "age": age,
                     "n_genes_with_ish": len(maps), "peak": comp.idxmax(),
                     "match": comp.idxmax() == claim,
                     "specificity": float(comp.max() / comp.sum()) if comp.sum() else None})
        print(rows[-1])
    df = pd.DataFrame(rows)
    df.to_csv(run / "fig3_cluster_regions.tsv", sep="\t", index=False)
    scored = df[df["n_genes_with_ish"] > 0].copy()
    scored["match"] = scored["match"].astype(bool)
    out = {"n_clusters_with_region_in_name": int(len(df)),
           "n_scored": int(len(scored)),
           "region_match_rate": float(scored["match"].mean()) if len(scored) else None,
           "chance_rate": 1.0 / len(REGIONS),
           "by_region": scored.groupby("claimed")["match"].mean().round(3).to_dict()}
    (run / "fig3_metrics.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    mode, run = sys.argv[1], Path(sys.argv[2])
    run.mkdir(parents=True, exist_ok=True)
    {"de": de, "allen": allen}[mode](run)
