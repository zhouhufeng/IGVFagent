# Zou 2024 — ChIP-Atlas 3.0 GATA1 hematopoietic case

## Citation

Zou A, Ohta T, Oki S. *Nucleic Acids Res* **52**: W45–W53 (2024).
DOI: [10.1093/nar/gkae358](https://doi.org/10.1093/nar/gkae358)

## Data

* **ChIP-Atlas**: > 70,000 SRA-mined experiments at https://chip-atlas.org
* **Bulk archive**: https://chip-atlas.dbcls.jp/data

## Headline workflow (paper)

1. Aggregate every public ChIP-seq for an antigen × cell-class slice.
2. Assemble an all-peaks BED at the requested -log10(q) threshold.
3. Cross-reference with GWAS / CRISPRi-validated CREs to identify TF-regulated risk loci.

## What IGVFagent reproduces

The `chipatlas` skill is IGVFagent's clean-room HTTP client for the ChIP-Atlas surface.

* `chipatlas list-antigens --genome hg38 --ag-class "TFs and others" --cell-class Blood` → enumerate TFs catalogued for Blood
* `chipatlas search --query "GATA1 K562"` → list reprocessed K562 GATA1 experiments
* `chipatlas assemble-bed --genome hg38 --ag-class "TFs and others" --antigen GATA1 --cell-class Blood --qval 05` → assembled all-peaks BED URL

Online-only.

## Ground truth to spot-check

* The GATA1 Blood assembled BED should resolve a downloadable URL (the actual BED is multi-MB).
* GATA1 occurrence in the antigen browser for cell-class=Blood should be in the top 20 by experiment count.
