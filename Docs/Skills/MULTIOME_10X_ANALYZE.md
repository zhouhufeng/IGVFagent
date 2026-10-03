# Skill: 10x Multiome Analysis

Clean-room Python implementations of the analytical methods covered in
the 10x Genomics
[`analysis_guides`](https://github.com/10XGenomics/analysis_guides) repo
(Epi multiome QC tutorial + interactive QC visualization), plus the
Stuart-lab Signac vignette extensions (TF-IDF/LSI, WNN, peak-to-gene
linkage).

No source from the 10x repo (which has no LICENSE) is copied; the math
is paraphrased from the published algorithm descriptions.

## Commands

```bash
# 1. Per-barcode ATAC QC from a fragments BED
igvfagent multiome qc-atac \
    --fragments fragments.tsv.gz \
    --tss-bed gencode_tss.bed \
    --peaks-bed atac_peaks.bed \
    --label run1

# 2. Joint QC: merge ATAC + RNA per-barcode tables, apply thresholds
igvfagent multiome joint-qc \
    --rna-qc run1_rna_qc.tsv \
    --atac-qc run1_atac_qc.tsv \
    --label run1

# 3. TF-IDF + truncated SVD on a peak × cell matrix (drops depth dim 1)
igvfagent multiome lsi --input atac_peaks.h5ad --label run1

# 4. Joint WNN embedding + UMAP via muon
igvfagent multiome wnn \
    --rna-h5ad rna.h5ad --atac-h5ad atac_with_lsi.h5ad \
    --cluster --label run1

# 5. Peak-to-gene correlation (Pearson over cells, ±500 kb window)
igvfagent multiome peak2gene \
    --rna-h5ad rna.h5ad --atac-h5ad atac.h5ad \
    --tss-bed gencode_tss.bed --window 500000 \
    --method pearson --label run1

# 6. Single-command showcase
igvfagent multiome showcase --fragments fragments.tsv.gz \
    --tss-bed tss.bed --peaks-bed peaks.bed \
    --rna-qc rna_qc.tsv --label run1
```

## QC thresholds (Signac convention)

| Metric | Pass criterion |
|---|---|
| RNA UMIs | 1,000 ≤ count ≤ 25,000 |
| RNA genes | ≥ 200 |
| RNA % mito | ≤ 20% |
| ATAC fragments | 1,800 ≤ count ≤ 100,000 |
| TSS enrichment | > 1 |
| Nucleosome signal | < 2 |
| FRIP | > 0.15 |

All thresholds are overridable via CLI flags.

## References

- 10x Genomics `analysis_guides` (algorithm spec only — no source copied):
  https://github.com/10XGenomics/analysis_guides
- Stuart-lab Signac (MIT) — `TSSEnrichment`, `NucleosomeSignal`,
  `RunSVD`, `FindMultiModalNeighbors`
- Hao et al., Cell 2021 — Seurat 5 WNN method
