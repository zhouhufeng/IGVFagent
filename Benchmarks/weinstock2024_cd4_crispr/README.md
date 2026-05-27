# Weinstock 2024 — CRISPR-KO causal network in primary CD4+ T cells

## Citation

Weinstock JS, ..., Pritchard JK. *Cell Genomics* **4**: 100693 (2024).
DOI: [10.1016/j.xgen.2024.100693](https://doi.org/10.1016/j.xgen.2024.100693)

## Data

* **GEO**: GSE171737 (84-gene KO + bulk RNA-seq)
* **GitHub**: [weinstockj/RNAseq-perturbation-CD4-pipeline](https://github.com/weinstockj/RNAseq-perturbation-CD4-pipeline) + [weinstockj/LLCB](https://github.com/weinstockj/LLCB)

## Headline workflow (paper)

1. CRISPR-KO of 84 immune-relevant genes in primary CD4+ T cells.
2. Bulk RNA-seq → LLCB causal-network inference.
3. Integrate with autoimmune GWAS — identify KMT2A as a Th17-IL2-JAK-STAT regulator; rs45480496 (upstream of KMT2A) is a Th17-enhancer autoimmune risk variant.

## What IGVFagent reproduces

* `perturb-catalog search-modality --modality crispr-screen --query KMT2A` finds the catalogued perturbation set.
* `network pkn-from-kg` materialises a signed protein-knowledge-network SIF rooted at the 84 KO genes from the proteomics KG; `network steiner` extracts a Prize-Collecting Steiner tree subnetwork connecting the seeds; `network viz` renders the result.

## Ground truth to spot-check

* The Steiner-tree subnetwork should connect KMT2A through STAT5 / JAK signalling toward IL2 targets — recapitulating the paper's Th17-IL2 axis.
* Search-modality result should return ≥ 1 CD4-context CRISPR-KO dataset.
