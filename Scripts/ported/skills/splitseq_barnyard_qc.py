# Ported by `igvfagent port register` on 2026-09-26 from https://github.com/Alex-Rosenberg/split-seq-pipeline @ a711b56dada4
# (split_seq/analysis.py) for rosenberg2018_splitseq. Unreviewed; provenance in Scripts/ported/registry.json.
"""SPLiT-seq species-mixing (barnyard) QC on a combined human+mouse DGE.

Port of Alex-Rosenberg/split-seq-pipeline split_seq/analysis.py:
  * species call per UBC: a species owns a UBC when >90% of its UMIs map to
    that genome, otherwise "multiplet" (generate_single_dge_report, the same
    rule as barnyard()'s counts1 > 9*counts2);
  * per-species "Median UMIs/Cell" and "Median Genes/Cell" over the UBCs
    assigned to that species, and "Number of Cells Detected".
Genes are counted per gene: the GEO DGEs split intronic reads into
``<gene>_INTRONIC_<SPECIES>`` columns, whereas the authors' pipeline assigns
intronic reads to their gene, so those columns are folded into their gene
before counting detected genes (UMI totals are unaffected).
Species purity is the paper's statement "x% of reads in human UBCs aligned to
the human genome": pooled species UMIs over pooled total UMIs of those UBCs.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "Benchmarks").is_dir())


def load_dge(path):
    d = sio.loadmat(str(path))
    genes = pd.Series(d["genes"]).str.strip()
    sample_type = pd.Series(d["sample_type"]).str.strip().values
    barcodes = pd.Series(d["barcodes"]).str.strip().values
    return sp.csr_matrix(d["DGE"]), genes, sample_type, barcodes


def barnyard_qc(X, genes, species=("HUMAN", "MOUSE")):
    base = genes.str.replace("_INTRONIC", "", regex=False).values
    uniq, inv = np.unique(base, return_inverse=True)
    fold = sp.csr_matrix((np.ones(len(base)), (np.arange(len(base)), inv)),
                         shape=(len(base), len(uniq)))
    Xg = X @ fold
    umi, ngenes = {}, {}
    for s in species:
        cols = genes.str.endswith("_" + s).values
        gcols = pd.Series(uniq).str.endswith("_" + s).values
        umi[s] = np.asarray(X[:, cols].sum(1)).ravel()
        ngenes[s] = np.asarray((Xg[:, gcols] > 0).sum(1)).ravel()
    umi = pd.DataFrame(umi)
    total = umi.sum(1).values
    assign = np.array(["multiplet"] * len(umi), dtype=object)
    for s in species:
        assign[np.where((umi[s].values / np.maximum(total, 1)) > 0.9)] = s
    out = {"n_ubcs": int(len(umi)), "n_multiplet": int((assign == "multiplet").sum())}
    out["multiplet_rate"] = out["n_multiplet"] / max(out["n_ubcs"], 1)
    out["single_species_rate"] = 1.0 - out["multiplet_rate"]
    for s in species:
        k = s.lower()
        own = assign == s
        out[f"{k}_n_cells"] = int(own.sum())
        out[f"{k}_median_umis"] = float(np.median(umi[s].values[own])) if own.any() else None
        out[f"{k}_median_genes"] = float(np.median(ngenes[s][own])) if own.any() else None
        out[f"{k}_purity_fraction"] = (float(umi[s].values[own].sum() / total[own].sum())
                                       if own.any() else None)
    return out, assign, umi


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="splitseq-barnyard-qc", description=__doc__.splitlines()[0])
    ap.add_argument("--mat", required=True, action="append",
                    help="GEO SPLiT-seq MATLAB DGE (repeatable)")
    ap.add_argument("--sample-type", default=None,
                    help="Restrict to UBCs of this sample_type (e.g. cell, nucleus)")
    ap.add_argument("--label", default="splitseq_barnyard_qc")
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args(argv)
    outdir = Path(a.outdir) if a.outdir else \
        ROOT / "Docs" / "SPLiTseq" / f"{time.strftime('%Y%m%d_%H%M%S')}_{a.label}"
    outdir.mkdir(parents=True, exist_ok=True)
    summary, rows = {"sample_type": a.sample_type, "libraries": {}}, []
    for m in a.mat:
        X, genes, st, bcs = load_dge(m)
        keep = np.ones(len(st), bool) if a.sample_type is None else (st == a.sample_type)
        res, assign, umi = barnyard_qc(X[keep], genes)
        lib = Path(m).name.split(".")[0]
        summary["libraries"][lib] = res
        rows.append(pd.DataFrame({"library": lib, "barcode": bcs[keep], "sample_type": st[keep],
                                  "human_umis": umi["HUMAN"].values,
                                  "mouse_umis": umi["MOUSE"].values, "species": assign}))
        print(lib, json.dumps(res))
    pd.concat(rows).to_csv(outdir / "species_calls.tsv", sep="\t", index=False)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"Wrote {outdir}/summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
