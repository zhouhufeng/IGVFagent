# Ported by `igvfagent port register` on 2026-09-26 from https://github.com/masai1116/SHARE-seq-alignment @ 11f2da021470
# (lib_size_sc_V4_species_mixing.R) for ma2020_shareseq. Unreviewed; provenance in Scripts/ported/registry.json.
#!/usr/bin/env python3
"""SHARE-seq human/mouse species-mixing calls from paired ATAC + RNA counts.

Reimplements the species-mixing readout of Ma et al., Cell 2020 (Fig. 1B-D,
Fig. S1C) on the GEO GSE140203 per-genome matrices. The authors' pipeline
(masai1116/SHARE-seq-alignment, lib_size_sc_V4_species_mixing.R) tabulates
per-barcode hg19 / mm10 unique molecules but leaves the cell call to the
reader, so the rule follows the paper text:

  * ATAC and RNA barcodes of the same cell share R1.R2.R3 and differ by a
    fixed PCR-index (P1) offset per sub-library; the offset is detected as
    the one maximising shared barcodes.
  * Cells are ATAC barcodes above the "steep drop-off" of the barcode-rank
    curve (STAR Methods, library-size comparison), found as the point of the
    log-log rank curve farthest from its chord. The RNA matrices are already
    limited to barcodes with >100 reads (STAR Methods), so no second RNA cut
    is applied unless ``--rna-min`` is given.
  * A cell is human (mouse) when the majority of its molecules are human
    (mouse) in BOTH modalities; a collision is a cell whose human fraction is
    between ``--mixed-lo`` and ``--mixed-hi`` in both modalities.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io


def _atac_totals(mtx, barcodes):
    m = scipy.io.mmread(gzip.open(mtx) if str(mtx).endswith(".gz") else mtx).tocsc()
    with (gzip.open(barcodes, "rt") if str(barcodes).endswith(".gz") else open(barcodes)) as fh:
        bcs = [l.strip() for l in fh if l.strip()]
    return pd.Series(np.asarray(m.sum(axis=0)).ravel(), index=bcs)


def _rna_totals(path):
    d = pd.read_csv(path, sep="\t", index_col=0)
    return d.sum(axis=0)


def _norm(bc):
    return bc.replace(",", ".")


def _shift(bc, off):
    return "%s%02d" % (bc[:-2], int(bc[-2:]) + off)  # barcode ends in P1.NN


def knee(tot):
    tot = np.asarray(tot, dtype=float)
    tot = tot[tot > 0]
    y = np.log10(np.sort(tot)[::-1])
    x = np.log10(np.arange(1, len(y) + 1))
    x0, y0, x1, y1 = x[0], y[0], x[-1], y[-1]
    dist = np.abs((y1 - y0) * x - (x1 - x0) * y + x1 * y0 - y1 * x0) / np.hypot(y1 - y0, x1 - x0)
    return float(10 ** y[int(np.argmax(dist))])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--atac-human-mtx", required=True)
    ap.add_argument("--atac-human-barcodes", required=True)
    ap.add_argument("--atac-mouse-mtx", required=True)
    ap.add_argument("--atac-mouse-barcodes", required=True)
    ap.add_argument("--rna-human", required=True, help="genes x cells TSV (.gz)")
    ap.add_argument("--rna-mouse", required=True)
    ap.add_argument("--atac-min", type=float, help="ATAC cell cut-off (default: knee)")
    ap.add_argument("--rna-min", type=float, default=0.0)
    ap.add_argument("--mixed-lo", type=float, default=0.2)
    ap.add_argument("--mixed-hi", type=float, default=0.8)
    ap.add_argument("--out", required=True, help="output directory")
    a = ap.parse_args(argv)

    ah = _atac_totals(a.atac_human_mtx, a.atac_human_barcodes)
    am = _atac_totals(a.atac_mouse_mtx, a.atac_mouse_barcodes)
    rh = _rna_totals(a.rna_human)
    rm = _rna_totals(a.rna_mouse)
    for s in (ah, am, rh, rm):
        s.index = [_norm(b) for b in s.index]
    atac = pd.DataFrame({"h": ah, "m": am}).fillna(0)
    rna = pd.DataFrame({"h": rh, "m": rm}).fillna(0)
    atac["tot"] = atac.h + atac.m
    rna["tot"] = rna.h + rna.m

    offs = range(-96, 97)
    shared = {o: len(set(_shift(b, o) for b in atac.index) & set(rna.index)) for o in offs}
    off = max(shared, key=shared.get)
    atac.index = [_shift(b, off) for b in atac.index]

    cut = a.atac_min if a.atac_min is not None else knee(atac["tot"])
    j = atac.join(rna, lsuffix="_atac", rsuffix="_rna", how="inner")
    j = j[(j.tot_atac >= cut) & (j.tot_rna >= a.rna_min)].copy()
    j["fh_atac"] = j.h_atac / j.tot_atac
    j["fh_rna"] = j.h_rna / j.tot_rna
    human = (j.fh_atac > 0.5) & (j.fh_rna > 0.5)
    mouse = (j.fh_atac < 0.5) & (j.fh_rna < 0.5)
    mixed = j.fh_atac.between(a.mixed_lo, a.mixed_hi) & j.fh_rna.between(a.mixed_lo, a.mixed_hi)
    j["call"] = np.where(mixed, "collision", np.where(human, "human", np.where(mouse, "mouse", "discordant")))

    own_rna = np.where(j.call == "human", j.h_rna, j.m_rna)
    single = j.call.isin(["human", "mouse"])
    summary = {
        "p1_offset_atac_to_rna": int(off),
        "shared_barcodes_at_offset": int(shared[off]),
        "atac_cell_cutoff": cut,
        "rna_min": a.rna_min,
        "n_cells": int(len(j)),
        "n_human": int((j.call == "human").sum()),
        "n_mouse": int((j.call == "mouse").sum()),
        "n_collision": int((j.call == "collision").sum()),
        "n_discordant": int((j.call == "discordant").sum()),
        "collision_rate": float(100 * (j.call == "collision").mean()),  # percent of cells
        "human_fraction": float((j.call == "human").sum() / max(single.sum(), 1)),
        "mean_rna_umis_single_species": float(own_rna[single.to_numpy()].mean()),
    }
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    j.to_csv(out / "species_calls.tsv", sep="\t")
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
