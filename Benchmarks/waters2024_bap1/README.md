# Waters 2024 — BAP1 saturation genome editing

## Citation

Waters AJ, ..., Adams DJ. *Nat Genet* **56**: 1434–1445 (2024).
DOI: [10.1038/s41588-024-01799-3](https://doi.org/10.1038/s41588-024-01799-3) · PMID: 38969833

## Data

* **MaveDB**: `urn:mavedb:00001210-a-1` (BAP1 SGE, NCI-H226)
* **ENA**: PRJEB72428 (raw reads)
* **GitHub**: [team113sanger/Waters_BAP1_SGE](https://github.com/team113sanger/Waters_BAP1_SGE)

## Headline workflow (paper)

1. Saturation genome editing of 18,108 BAP1 SNVs in NCI-H226 mesothelioma cells.
2. Quantify each variant's effect on cell fitness over time.
3. Train a benign/pathogenic classifier — paper reports >99 % sensitivity, >98 % specificity on a ClinVar truth-set.

## What IGVFagent reproduces

The `mavedb` skill pulls the canonical MaveDB scoreset CSV, parses every HGVSp / per-position-tile / single-letter-WT-pos-ALT row, and maps each variant to **chr / pos / ref / alt** GRCh38 coordinates via the Ensembl REST canonical transcript for BAP1 (ENST00000460680). The output is a TSV + VCF + JSON summary that any downstream consumer can intersect with ClinVar, gnomAD, IGVF Catalog, etc.

This benchmark is **online-only** — no local downloads required; the entire run is `MaveDB → Ensembl REST → on-disk TSV/VCF`.

## Ground truth to spot-check

* 18,108 SNVs perturbed in total
* ~5,665 LOF + 531 GOF variants reported (after fitness filtering)
* The UCH domain (exons 1–9, 15–17) should be the most depleted region
* ~99.8 % of ClinVar-pathogenic truth-set variants are LOF-depleted

`expected.json` declares hard-range checks against the row-count / type-breakdown that `mavedb map-scoreset` emits, plus an artefact-existence check for the VCF.
