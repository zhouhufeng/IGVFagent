# Yao 2024 — ENCODE4 multicenter noncoding CRISPRi screens

## Citation

Yao D, Tycko J, Oh JW, ..., Engreitz JM. *Nat Methods* **21**: 1980–1992 (2024).
DOI: [10.1038/s41592-024-02216-7](https://doi.org/10.1038/s41592-024-02216-7) · PMID: 38504114

## Data

* **ENCODE portal**: 108 CRISPRi screens (K562, HepG2, Jurkat) covering 24.85 Mb of regulatory DNA
* **GitHub**: [EngreitzLab/CRISPRi_noncoding_analysis](https://github.com/EngreitzLab/CRISPRi_noncoding_analysis) — CASA pipeline
* Includes ABC-model validation: GATA1, MYC, FADS1/2 loci in K562

## Headline workflow (paper)

1. CRISPRi tiling screens at scale across three cell lines.
2. Call significant elements using CASA, CRISPR-SURF, MAGeCK, RELICS (benchmarked head-to-head).
3. Cross-reference CRE → gene calls with ABC-model predictions.

## What IGVFagent reproduces

* **Online step** (this benchmark's scope): `encode retrieve` enumerates the ENCODE CRISPRi MeasurementSets / AnalysisSets, building a manifest of each screen's experiment accession + target cell type + readout.
* **Optional local step**: if you've downloaded a single screen's guide-count table, `crispri analyze-local --input <counts.tsv>` runs the IGVFagent log-fold-change / Mann-Whitney pipeline. (Not enabled in this benchmark — see `martyn2025_variant_flowfish` for an end-to-end CRISPR-screen reproducibility test that exercises the analysis step.)

## Ground truth to spot-check

* 108 CRISPRi screens enumerated
* 332 K562 CRE–gene links confirmed across 24.85 Mb perturbed
* GATA1 +24/+58 kb HS-sites should appear in any K562 manifest's preferred-target list
