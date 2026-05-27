# Joung 2025 — TF Perturb-seq in fibroblasts

## Citation

Joung J, ..., Zhang F. *Nat Genet* **57**: 828–838 (2025).
DOI: [10.1038/s41588-025-02283-2](https://doi.org/10.1038/s41588-025-02283-2)

## Data

* **GEO**: GSE273694
* **Single-Cell Portal**: SCP2169
* **Related code**: [josephreplogle/perturb-seq](https://github.com/josephreplogle/perturb-seq)

## Headline workflow (paper)

1. CRISPRa of 1,836 TFs (10,979 sgRNAs) in primary fibroblasts.
2. Perturb-seq readout: scRNA + sgRNA identity.
3. Identify TFs that recover known in-vivo fibroblast states (KLF4, KLF5 drive a universal-state cluster).

## What IGVFagent reproduces

* **Online step**: `perturb-catalog search-modality --modality perturb-seq` queries the Perturbation Catalogue for TF Perturb-seq datasets.
* **Optional local step**: `sc-analyze pipeline` on the downloaded h5ad → per-perturbation pseudobulk DE.

## Ground truth to spot-check

* KLF4 + KLF5 perturbations land in a shared "universal fibroblast state" cluster
* The dataset query should return at least 1 fibroblast-context Perturb-seq dataset
