"""Fit the LLCB causal network over the paper's 84-gene panel.

Input:  Data/Weinstock2024/processed/network_input.csv
Output: Data/Weinstock2024/processed/edges.csv
"""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
PROC = ROOT / "Data/Weinstock2024/processed"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from llcb import fit_llcb  # noqa: E402


def main() -> None:
    data = pd.read_csv(PROC / "network_input.csv", dtype={"donor": str})
    targets = sorted(set(data.columns) - {"donor", "intervention"})
    print(f"{len(targets)} target genes, {len(data)} samples")

    edges = fit_llcb(data, targets)
    edges.to_csv(PROC / "edges.csv", index=False)
    print(f"wrote {PROC / 'edges.csv'}  ({len(edges)} directed pairs)")

    print("\nvs. paper Table/Fig 2 (350, 211, 151 edges at |beta|>0.020/0.025/0.030):")
    for threshold, paper_n in ((0.020, 350), (0.025, 211), (0.030, 151)):
        n = (edges["estimate"].abs() > threshold).sum()
        print(f"  |estimate| > {threshold}: {n} edges  (density {n/len(edges):.1%})  [paper: {paper_n}]")

    print("\nvs. paper's own stated LFSR pathology (67% density at LFSR<5e-3,")
    print("which is *why the paper itself rejected LFSR-based edge-calling*):")
    n = (edges["lsfr"] < 5e-3).sum()
    print(f"  lsfr < 0.005: {n} edges  (density {n/len(edges):.1%})  [paper: 67%]")


if __name__ == "__main__":
    main()
