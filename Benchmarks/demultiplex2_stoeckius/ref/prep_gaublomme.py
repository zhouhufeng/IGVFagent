#!/usr/bin/env python3
"""Gaublomme 2019 nuclei hashing (human_st) inputs for the deMULTIplex2 Table S3
re-benchmark, following deMULTIplex2-benchmark @ ac0fa17 gaublomme_preprocess.R
and gaublomme_benchmarking.R. Inputs come from the regevlab/demuxem docker
image (software/inputs/), as the paper states (Table S2).

  tags : experiment1_human_st_ADT.csv, barcodes suffixed "-1", restricted to
         barcodes with >= 100 RNA UMIs (experiment1_human_st_raw_GRCh38_premrna)
         and > 0 tag UMIs.
  truth: experiment1_demuxlet.best -> AMB = "unknown"; PRB.DBL <= 0.99 and not
         AMB = singlet (SNG.1ST); otherwise "doublet". All demuxlet barcodes are
         scored (confusion_stats indexes calls by the truth's names).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio

src, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
adt = pd.read_csv(src / "experiment1_human_st_ADT.csv", index_col=0).T
adt.index = adt.index + "-1"
mdir = src / "experiment1_human_st_raw_GRCh38_premrna"
m = sio.mmread(str(mdir / "matrix.mtx")).tocsc()
bcs = pd.read_csv(mdir / "barcodes.tsv", header=None)[0].values
umi = np.asarray(m.sum(axis=0)).ravel()
rna_cells = set(bcs[umi >= 100])
tag = adt.loc[[c for c in adt.index if c in rna_cells]]
tag = tag[tag.sum(axis=1) > 0]
tag.to_csv(out / "gaublomme_c1_tags.csv")
d = pd.read_csv(src / "experiment1_demuxlet.best", sep="\t")
amb = d["BEST"].str.contains("AMB")
sng = (d["PRB.DBL"] <= 0.99) & ~amb
truth = np.where(sng, d["SNG.1ST"], np.where(amb, "unknown", "doublet"))
pd.DataFrame({"cell": d["BARCODE"], "truth": truth}).to_csv(out / "gaublomme_truth.csv", index=False)
mp = {t: t for t in tag.columns}
(out / "gaublomme_tag_mapping.json").write_text(json.dumps(mp, indent=1))
man = {"gaublomme": {"captures": ["gaublomme_c1_tags.csv"], "truth": "gaublomme_truth.csv",
                     "tag_mapping": mp, "n_cells": int(len(d))}}
(out / "gaublomme_manifest.json").write_text(json.dumps(man, indent=1))
print("tag matrix", tag.shape, "| truth", pd.Series(truth).value_counts().to_dict())
