#!/usr/bin/env python3
"""Convert the deMULTIplex2-bundled Stoeckius PBMC HTO matrix (.RData) to CSV.

deMULTIplex2 ships ``data/stoeckius_pbmc.RData`` — the Stoeckius 2018 8-donor
PBMC cell-hashing tag count matrix (cells x 8 HTO). IGVFagent's ``multiseq``
skill consumes a plain cells x tags CSV, so we convert once (pure Python via
pyreadr; no R needed).
"""
from pathlib import Path
import pyreadr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
D = ROOT / "Benchmarks/_data/demultiplex2_stoeckius"
SRC = D / "stoeckius_pbmc.RData"
DST = D / "stoeckius_tags.csv"

if not SRC.is_file():
    raise SystemExit(f"missing {SRC} — run.sh downloads it from the deMULTIplex2 repo")

obj = pyreadr.read_r(str(SRC))
df = obj["stoeckius_pbmc"]
df.to_csv(DST)
print(f"wrote {DST.name}: {df.shape[0]} cells x {df.shape[1]} HTOs "
      f"({int(df.values.sum()):,} total UMIs)")
