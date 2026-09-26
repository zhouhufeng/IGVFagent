# mitra2024_multi_regression

> **Real local run (2026-09-25).** The `multiome_peak2gene` chain has been executed end to end on the paper's own BMMC/GSE194122 subset — see Concordance below. The 12 `[UNCONFIRMED]` checks in `expected.json` remain unconfirmed **by design**, not for lack of effort: they measure a different quantity than our proxy computes (see Honest caveats).

## Paper

**Single-cell multi-ome regression models identify functional and disease-associated enhancers and enable chromatin potential analysis**  
Mitra Sneha; Malik Rohan; Wong Wilfred; Rahman Afsana; Hartemink Alexander J.; Pritykin Yuri; Dey Kushal K.; Leslie Christina S.  
*Nature Genetics* 2024 · doi:[10.1038/s41588-024-01689-8](https://doi.org/10.1038/s41588-024-01689-8) · PMID 38514783 · PMC11018525

Resolver confidence: **1.00** (resolved).

## Data sources found in the paper

| Repository | Accession | In Data Availability | Mentions |
|---|---|---|---|
| ENCODE | `ENCSR233SQS` | ✓ | 3 |
| ENCODE | `ENCSR233SQG` | ✓ | 3 |
| NCBI GEO | `GSE194122` | ✓ | 3 |
| NCBI GEO | `GSE140203` | ✓ | 3 |
| NCBI GEO | `GSE178454` | ✓ | 3 |
| NCBI GEO | `GSE162170` | ✓ | 3 |
| GitHub | `snehamitra/SCARlink` | ✓ | 6 |
| GitHub | `jmacdon/LDblocks_GRCh38` | ✓ | 5 |
| GitHub | `GreenleafLab/brainchromatin` | ✓ | 3 |
| Zenodo | `10.5281/zenodo.10481793` | — | 2 |

## What IGVFagent does

**Route:** `multiome_peak2gene` — 10x Multiome peak→gene cis-regulatory linkage  
**Skill output dir:** `Docs/Multiome10x/`  
**Modelled on:** trevino2021_cortex_multiome / mitra2024_scarlink

`run.sh` (1) downloads `GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad` from GEO (~2.7 GB, no auth), (2) splits it by `var['feature_types']` into RNA/ATAC h5ads restricted to the paper's exact 10 donors and 5 cell types, (3) runs `igvfagent multiome peak2gene` genome-wide (per-peak Pearson correlation within a 500 kb window of each gene's TSS, using `Data/scE2G/resources/CollapsedGeneBounds.hg38.TSS500bp.bed`), and (4) scores IGVFagent's own measured quantities into `mitra2024_multi_regression_concordance_metrics.json`.

```bash
bash Benchmarks/mitra2024_multi_regression/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark mitra2024_multi_regression
```

## Concordance

Run 2026-09-25, 7,642 cells (site1_donor1/2/3, site2_donor1/4/5, site3_donor6/7/10, site4_donor9 × HSC / MK/E prog / Proerythroblast / Erythroblast / Normoblast — exactly the paper's stated BMMC subset):

| Check | IGVFagent measured | Status |
|---|---|---|
| Donors matched to paper's list | 10 / 10 | ✓ |
| Cell types matched to paper's list | 5 / 5 | ✓ |
| Genome-wide scan not truncated | 923,107 candidate pairs (cap 2,000,000) | ✓ |
| peak2gene link table written | `Docs/Multiome10x/*_peak2gene.tsv` | ✓ |
| Genes tested | 11,614 / 13,431 RNA features | — |
| Peaks tested | 114,078 / 116,490 ATAC features | — |
| Significant pairs (padj < 0.05) | 16,650 | — |
| Genes with ≥1 significant link | 4,240 | — |
| Genes with ≥1 significant *positive* link | 3,513 (61% of significant pairs are positive) | — |

**4 / 4 confirmed structural checks pass.** The paper's own numeric claims (1,655 genes for BMMC after SCARlink's filtering; 785→749 genes after its chromatin-potential-specific filter) are listed in `expected.json` but **not** matched against the 4,240/3,513 above — see below.

## Honest caveats

* **This is not SCARlink.** SCARlink fits a per-gene Poisson regression over all nearby peaks jointly (with its own inclusion criteria: minimum expression, minimum peak count, model convergence). `multiome peak2gene` here does a much simpler, more permissive per-peak-per-gene Pearson correlation test. The two pipelines answer related but non-identical questions, and their gene counts are not on the same scale — our 4,240 "genes with a significant link" is not evidence against, or for, the paper's 1,655.
* **Only the BMMC arm is reproduced.** The paper's PBMC multi-ome (from 10X Genomics, no accession given in the text) and the mouse skin / developing cortex / pancreas / pituitary datasets are not fetched or processed by this benchmark.
* **12 unconfirmed checks remain unconfirmed on purpose**, not for lack of effort — see the metric-mismatch point above. Extracted from prose by regex/LLM; each carries its source quote in `expected.json`.
* This proves the chain runs end-to-end on the paper's real, correctly-subset data and produces plausible cis peak→gene links (61% of significant pairs are positive-correlation, consistent with enhancer-activation direction). It does not prove SCARlink's specific reported numbers.

## Provenance

`provenance.json` in this directory holds the full resolve / harvest / route record, including every source URL consulted.
