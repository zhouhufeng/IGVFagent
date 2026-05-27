# Zheng 2024 — in-vivo AAV Perturb-seq, mouse cortex

## Citation

Zheng X, ..., Jin X. *Cell* **187**: 3236–3248 (2024).
DOI: [10.1016/j.cell.2024.04.050](https://doi.org/10.1016/j.cell.2024.04.050) · PMID: 38772369

## Data

* **GEO**: GSE249416 (~30K cells, AAV-delivered Perturb-seq)
* **GitHub**: https://github.com/jinlabneurogenomics

## Headline workflow (paper)

1. AAV-delivered sgRNAs against Foxg1 / Nr2f1 / Tbr1 / Tcf4 in fetal mouse cortex.
2. scRNA-seq + sgRNA identity.
3. Cell-type-resolved TF effects — Foxg1 loss de-represses other TF networks specifically in Layer-6 corticothalamic neurons.

## What IGVFagent reproduces

This is a **local-input-required** benchmark: download GSE249416 h5ad first, then `sc-analyze pipeline --input` runs the QC → normalize → PCA → UMAP → cluster → markers chain.

## Ground truth to spot-check

* Foxg1-KO cells de-repress Tbr1 + Nr2f1 target genes specifically in Layer-6 markers (Fezf2+, Ldb2+)
* Pipeline should produce ≥ 1 markers TSV
