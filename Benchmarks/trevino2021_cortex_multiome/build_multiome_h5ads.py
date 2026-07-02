#!/usr/bin/env python3
"""Build panel-restricted RNA + ATAC AnnData + TSS BED for `multiome peak2gene`.

Trevino 2021 (GSE162170) ships the paired multiome subset as dense, headerless
count matrices:

  * RNA  : 34,104 genes (Ensembl) x 8,981 cells   (gene-major, cell IDs in header)
  * ATAC : 467,315 peaks x 8,981 cells            (peak-major, positional; row i
                                                    <-> consensus_peaks row i;
                                                    columns share the RNA/metadata
                                                    cell order)

`multiome peak2gene` densifies whatever it is handed, so we restrict both
matrices to a curated cortical-neurogenesis gene panel and the consensus peaks
within +/- WINDOW of those genes' TSS. That keeps the correlation tractable on a
laptop while exercising the *exact* peak->gene cis-correlation method that
produces Trevino's enhancer-gene links.

TSS coordinates (GRCh38) are resolved from the Ensembl REST API and cached.
"""
import gzip
import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import anndata as ad
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
D = ROOT / "Benchmarks/_data/trevino2021_cortex_multiome"
RNA_TSV = D / "GSE162170_multiome_rna_counts.tsv.gz"
ATAC_TSV = D / "GSE162170_multiome_atac_counts.tsv.gz"
PEAKS = D / "GSE162170_multiome_atac_consensus_peaks.txt.gz"
META = D / "GSE162170_multiome_cell_metadata.txt.gz"
CLUSTERS = D / "GSE162170_multiome_cluster_names.txt.gz"
TSS_CACHE = D / "panel_tss_grch38.json"

WINDOW = 500_000  # bp around each TSS (matches the skill default)

# Cortical-neurogenesis marker panel (Trevino 2021 Figs 1-4 lineage genes).
PANEL = [
    # radial glia / progenitors
    "PAX6", "SOX2", "HES1", "VIM", "GLI3", "HOPX", "TNC",
    # cycling progenitors
    "MKI67", "TOP2A",
    # intermediate progenitors
    "EOMES", "NEUROG2", "PPP1R17",
    # excitatory neurons (deep -> upper layer)
    "NEUROD2", "NEUROD6", "BCL11B", "TBR1", "FEZF2", "SATB2", "FOXP2", "STMN2",
    # interneurons
    "DLX2", "GAD1", "GAD2",
    # glia
    "AQP4", "OLIG1", "OLIG2",
]


def resolve_tss() -> dict:
    """symbol -> {ensembl, chrom, tss} from Ensembl REST (GRCh38), cached."""
    if TSS_CACHE.is_file():
        return json.loads(TSS_CACHE.read_text())
    out = {}
    for sym in PANEL:
        url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{sym}?content-type=application/json"
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
            j = json.load(urllib.request.urlopen(req, timeout=30))
            chrom = "chr" + str(j["seq_region_name"])
            tss = int(j["start"]) if int(j["strand"]) == 1 else int(j["end"])
            out[sym] = {"ensembl": j["id"], "chrom": chrom, "tss": tss}
            print(f"  TSS {sym:8s} {j['id']} {chrom}:{tss}")
        except Exception as e:  # noqa
            print(f"  !! {sym}: {e}")
        time.sleep(0.1)
    TSS_CACHE.write_text(json.dumps(out, indent=2))
    return out


def main():
    for f in (RNA_TSV, ATAC_TSV, PEAKS, META):
        if not f.is_file():
            sys.exit(f"missing {f} — run the download steps in run.sh first")

    print("==> resolving panel TSS (Ensembl GRCh38)")
    tss = resolve_tss()
    ens2sym = {v["ensembl"]: s for s, v in tss.items()}
    panel_ens = set(ens2sym)

    # ---- metadata / cell order + cluster names --------------------------------
    meta = pd.read_csv(META, sep="\t")
    cell_ids = meta["Cell.ID"].astype(str).tolist()
    cl = pd.read_csv(CLUSTERS, sep="\t")
    id2name = dict(zip(cl["Cluster.ID"].astype(str), cl["Cluster.Name"].astype(str)))
    # seurat_clusters already carries the "c0".."cN" Cluster.ID strings
    sc_names = [id2name.get(str(c), str(c)) for c in meta["seurat_clusters"]]
    obs = pd.DataFrame({"cluster_id": meta["seurat_clusters"].astype(str).values,
                        "cluster_name": sc_names,
                        "sample_age": meta["Sample.Age"].astype(str).values},
                       index=cell_ids)

    # ---- RNA: keep only panel genes ------------------------------------------
    print("==> reading RNA (panel genes only)")
    rna_rows, rna_genes = [], []
    with gzip.open(RNA_TSV, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")  # cell barcodes
        assert header[0] == cell_ids[0], "RNA header / metadata cell order mismatch"
        for ln in fh:
            gid, _, rest = ln.partition("\t")
            if gid in panel_ens:
                rna_rows.append(np.fromstring(rest, sep="\t", dtype=np.float32))
                rna_genes.append(gid)
    rna_mat = np.vstack(rna_rows).T  # cells x genes
    rna = ad.AnnData(X=sp.csr_matrix(rna_mat), obs=obs.copy())
    rna.var_names = rna_genes
    rna.var["symbol"] = [ens2sym.get(g, g) for g in rna_genes]
    print(f"    RNA: {rna.n_obs} cells x {rna.n_vars} panel genes")

    # ---- peaks within window of any panel TSS --------------------------------
    print("==> selecting consensus peaks near panel TSS")
    peaks = pd.read_csv(PEAKS, sep="\t")
    pmid = ((peaks["start"] + peaks["end"]) // 2).to_numpy()
    pchr = peaks["seqnames"].astype(str).to_numpy()
    keep_idx, keep_names = [], []
    by_chrom = {}
    for v in tss.values():
        by_chrom.setdefault(v["chrom"], []).append(v["tss"])
    for i in range(len(peaks)):
        ts = by_chrom.get(pchr[i])
        if not ts:
            continue
        if min(abs(pmid[i] - t) for t in ts) <= WINDOW:
            keep_idx.append(i)
            keep_names.append(f"{pchr[i]}:{int(peaks['start'][i])}-{int(peaks['end'][i])}")
    keep_set = set(keep_idx)
    print(f"    {len(keep_idx)} of {len(peaks)} peaks within +/-{WINDOW//1000}kb of panel TSS")

    # ---- ATAC: stream, keep only selected peak rows --------------------------
    print("==> streaming ATAC (selected peak rows only)")
    atac_rows = []
    with gzip.open(ATAC_TSV, "rt") as fh:
        for i, ln in enumerate(fh):
            if i in keep_set:
                atac_rows.append(np.fromstring(ln, sep="\t", dtype=np.float32))
    atac_mat = np.vstack(atac_rows).T  # cells x peaks
    atac = ad.AnnData(X=sp.csr_matrix(atac_mat), obs=obs.copy())
    atac.var_names = keep_names
    print(f"    ATAC: {atac.n_obs} cells x {atac.n_vars} peaks")

    # ---- write outputs --------------------------------------------------------
    rna.write_h5ad(D / "multiome_rna_panel.h5ad", compression="gzip")
    atac.write_h5ad(D / "multiome_atac_panel.h5ad", compression="gzip")
    with open(D / "panel_tss.bed", "w") as fh:
        for sym, v in tss.items():
            fh.write(f"{v['chrom']}\t{v['tss']}\t{v['tss']+1}\t{v['ensembl']}\n")
    print(f"==> wrote multiome_rna_panel.h5ad / multiome_atac_panel.h5ad / panel_tss.bed")


if __name__ == "__main__":
    main()
