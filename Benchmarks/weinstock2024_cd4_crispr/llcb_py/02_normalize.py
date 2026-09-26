"""Normalize the raw dedup counts with DESeq2-equivalent size factors +
variance-stabilizing transform (VST), via pydeseq2 -- a from-source Python
port of DESeq2's own algorithms (Love, Huber & Anders 2014), used here in
place of the authors' R DESeq2 call since this cluster has no R/Bioconductor
install. Mirrors the paper's own normalization step (their pipeline's
`vst_normalized_counts_transpose.tsv`), just produced with a different but
numerically-equivalent implementation of the same method.

Input:  Data/Weinstock2024/processed/{counts.csv.gz,sample_meta.csv}
Output: Data/Weinstock2024/processed/vst.csv.gz  (genes x samples)
"""
from pathlib import Path

import pandas as pd
from pydeseq2.dds import DeseqDataSet

ROOT = Path(__file__).resolve().parents[3]
PROC = ROOT / "Data/Weinstock2024/processed"


def main() -> None:
    counts = pd.read_csv(PROC / "counts.csv.gz", index_col=0)
    meta = pd.read_csv(PROC / "sample_meta.csv", index_col="sample", dtype={"donor": str})

    # pydeseq2 wants samples x genes, and only non-all-zero genes.
    counts_t = counts.T
    counts_t = counts_t.loc[meta.index]
    keep = counts_t.sum(axis=0) > 0
    print(f"dropping {(~keep).sum()} all-zero genes of {len(keep)}")
    counts_t = counts_t.loc[:, keep]

    dds = DeseqDataSet(
        counts=counts_t,
        metadata=meta,
        design_factors="donor",
        refit_cooks=False,
    )
    dds.vst(use_design=False)  # intercept-only dispersion trend, matches DESeq2's blind VST

    vst = dds.layers["vst_counts"]
    vst_df = pd.DataFrame(vst, index=counts_t.index, columns=counts_t.columns).T
    vst_df.to_csv(PROC / "vst.csv.gz")
    print(f"wrote {PROC / 'vst.csv.gz'}  shape={vst_df.shape}")


if __name__ == "__main__":
    main()
