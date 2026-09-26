#!/usr/bin/env python3
"""Prepare the Howitt/Maksimovic (Oshlack/hashtag-demux-paper @ 3be94bb) lung
cell-line and BAL batch 1-3 hashtag matrices + genetic-donor truth for the
deMULTIplex2 Table S3 re-benchmark.

Per capture: <out>/<dataset>_c<k>_tags.csv (cells x tags) and one
<out>/<dataset>_truth.csv (cell, truth) for all captures of the dataset, with
cells prefixed "c<k>_" because barcodes repeat across captures. Truth labels
are the genetic donors (vireo) from the Oshlack repo; "Doublet" becomes
"doublet" (confusion_stats' true.multiplet), everything else is kept.
The tag -> donor mapping is the one in Oshlack analysis/*.Rmd
(tag "BAL 0k"/"CL 0k" <-> k-th letter donor).
"""
import argparse
import json
import string
from pathlib import Path

import pandas as pd

DATASETS = {
    "lung_cellline": [("cell_line_data/lmo_counts_capture%d.csv" % k,
                       "cell_line_data/lmo_donors_capture%d.csv" % k) for k in (1, 2, 3)],
    **{f"bal_batch{b}": [(f"BAL_data/batch{b}_c{c}_hto_counts.csv",
                          f"BAL_data/batch{b}_c{c}_donors.csv") for c in (1, 2)]
       for b in (1, 2, 3)},
}


def tag_to_donor(tag: str) -> str:
    """'BAL 09' -> 'BAL I', 'CL 02' -> 'CL B' (Oshlack donor lists)."""
    prefix, num = tag.rsplit(" ", 1)
    return f"{prefix} {string.ascii_uppercase[int(num) - 1]}"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, help="Oshlack/hashtag-demux-paper checkout")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    data, out = Path(a.repo) / "data", Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for ds, files in DATASETS.items():
        truth_rows, captures = [], []
        for k, (cnt, don) in enumerate(files, 1):
            m = pd.read_csv(data / cnt, index_col=0)
            m = m.T  # rows = cells
            m.index = [f"c{k}_{c}" for c in m.index]
            d = pd.read_csv(data / don, index_col=0)
            d.index = [f"c{k}_{c}" for c in d["Barcode"]]
            if not set(d.index) <= set(m.index):
                raise SystemExit(f"{ds} c{k}: donor barcodes missing from counts")
            m = m.loc[d.index]
            p = out / f"{ds}_c{k}_tags.csv"
            m.to_csv(p)
            captures.append(p.name)
            truth_rows.append(pd.DataFrame({"cell": d.index,
                                            "truth": d["genetic_donor"].replace({"Doublet": "doublet"})}))
        truth = pd.concat(truth_rows)
        truth.to_csv(out / f"{ds}_truth.csv", index=False)
        tags = list(m.columns)
        mapping = {t: tag_to_donor(t) for t in tags}
        (out / f"{ds}_tag_mapping.json").write_text(json.dumps(mapping, indent=1))
        manifest[ds] = {"captures": captures, "truth": f"{ds}_truth.csv",
                        "tag_mapping": mapping, "n_cells": int(len(truth))}
        print(ds, len(truth), "cells;", mapping)
    (out / "oshlack_manifest.json").write_text(json.dumps(manifest, indent=1))


if __name__ == "__main__":
    main()
