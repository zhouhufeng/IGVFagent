#!/usr/bin/env python3
"""Rosenberg 2018 Fig. 1B-D from the GEO species-mixing DGEs (GSM3017262-65).

Runs the registered port ``igvfagent splitseq-barnyard-qc`` (port of the
authors' split_seq/analysis.py) on cells and on nuclei, then computes the
Fig. 1D fresh-vs-frozen gene-expression correlation. Writes fig1_metrics.json
into the run directory given as argv[1].

Fig. 1D definition (the paper does not state it; fixed here before comparing):
pseudobulk of human-assigned whole cells, summed UMIs per human feature,
CPM-normalised, Pearson r of log1p(CPM), same-day 3000-UBC library vs frozen
1000-UBC library (the two largest libraries).
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "Data" / "rosenberg2018"
LIBS = {"fresh3000": "GSM3017262_same_day_cells_nuclei_3000_UBCs",
        "fresh300": "GSM3017263_same_day_cells_nuclei_300_UBCs",
        "frozen1000": "GSM3017264_frozen_preserved_cells_nuclei_1000_UBCs",
        "frozen200": "GSM3017265_frozen_preserved_cells_nuclei_200_UBCs"}
IGVF = str(Path(sys.executable).with_name("igvfagent"))


def run_port(sample_type, outdir):
    cmd = [IGVF, "splitseq-barnyard-qc", "--sample-type", sample_type, "--outdir", str(outdir)]
    for lib in LIBS.values():
        cmd += ["--mat", str(D / f"{lib}.mat")]
    subprocess.run(cmd, check=True)
    return json.loads((outdir / "summary.json").read_text())["libraries"]


def pseudobulk(lib):
    d = sio.loadmat(str(D / f"{lib}.mat"))
    g = pd.Series(d["genes"]).str.strip()
    st = pd.Series(d["sample_type"]).str.strip().values
    X = d["DGE"].tocsr()
    hum = g.str.endswith("_HUMAN").values
    h = np.asarray(X[:, hum].sum(1)).ravel()
    m = np.asarray(X[:, g.str.endswith("_MOUSE").values].sum(1)).ravel()
    sel = (st == "cell") & (h > 9 * m)
    return pd.Series(np.asarray(X[sel][:, hum].sum(0)).ravel(), index=g[hum].values)


def main():
    run = Path(sys.argv[1])
    run.mkdir(parents=True, exist_ok=True)
    for lib in LIBS.values():
        if not (D / f"{lib}.mat").exists():
            subprocess.run(["gunzip", "-k", str(D / f"{lib}.mat.gz")], check=True)
    cells = run_port("cell", run / "barnyard_cells")
    nuclei = run_port("nucleus", run / "barnyard_nuclei")
    a, b = pseudobulk(LIBS["fresh3000"]), pseudobulk(LIBS["frozen1000"])
    df = pd.concat([a, b], axis=1).fillna(0)
    cpm = df / df.sum() * 1e6
    r = float(np.corrcoef(np.log1p(cpm.iloc[:, 0]), np.log1p(cpm.iloc[:, 1]))[0, 1])
    L = LIBS
    out = {
        # Fig. 1B: the 1,758-UBC whole-cell barnyard
        "fig1b_single_species_rate": cells[L["fresh3000"]]["single_species_rate"],
        "fig1b_collision_rate": cells[L["fresh3000"]]["multiplet_rate"],
        # saturating-coverage library (text: 15,365 / 5,498 human; 12,243 / 4,497 mouse; 99.6 / 99.0 %)
        "deep_human_median_umis": cells[L["fresh300"]]["human_median_umis"],
        "deep_human_median_genes": cells[L["fresh300"]]["human_median_genes"],
        "deep_mouse_median_umis": cells[L["fresh300"]]["mouse_median_umis"],
        "deep_mouse_median_genes": cells[L["fresh300"]]["mouse_median_genes"],
        "deep_human_purity_fraction": cells[L["fresh300"]]["human_purity_fraction"],
        "deep_mouse_purity_fraction": cells[L["fresh300"]]["mouse_purity_fraction"],
        # Fig. 1C: median human UMIs (fresh cells / frozen cells / nuclei / frozen nuclei)
        "fig1c_fresh_cells_median_umis": cells[L["fresh300"]]["human_median_umis"],
        "fig1c_frozen_cells_median_umis": cells[L["frozen200"]]["human_median_umis"],
        "fig1c_fresh_nuclei_median_umis": nuclei[L["fresh300"]]["human_median_umis"],
        "fig1c_frozen_nuclei_median_umis": nuclei[L["frozen200"]]["human_median_umis"],
        # Fig. 1D
        "fig1d_frozen_vs_fresh_correlation": r,
        "fig1d_n_features": int(len(df)),
    }
    (run / "fig1_metrics.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
