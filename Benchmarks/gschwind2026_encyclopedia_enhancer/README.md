# gschwind2026_encyclopedia_enhancer

[![paper](https://img.shields.io/badge/Nature-2026-blue)](https://doi.org/10.1038/s41586-026-10781-4)
[![PMID](https://img.shields.io/badge/PMID-42457959-blue)](https://pubmed.ncbi.nlm.nih.gov/42457959/)
[![coverage](https://img.shields.io/badge/reproduction-5%2F5%20analyses-brightgreen)]()

## Paper

Gschwind AR, Mualim KS, Karbalayghareh A, Sheth MU, Dey KK, Jagoda E, Nurtdinov RN, Xi W, … Engreitz JM. **An encyclopedia of human enhancer–gene regulatory interactions.** *Nature* 2026. doi:[10.1038/s41586-026-10781-4](https://doi.org/10.1038/s41586-026-10781-4) · PMID 42457959 · PMC13471189

The ENCODE-rE2G encyclopedia paper: 92,176,227 predicted enhancer-gene regulatory interactions across 1,458 ENCODE biosamples (369 cell types/tissues), built on a new supervised model (ENCODE-rE2G) benchmarked against a combined CRISPR-perturbation dataset, eQTLs and fine-mapped GWAS loci.

## Data and code

| | Source |
|---|---|
| Authors' own toolchain | `EngreitzLab/ENCODE_rE2G@d039062` (model + 9 pretrained variants), `broadinstitute/ABC-Enhancer-Gene-Prediction`, `EngreitzLab/CRISPR_comparison`, `EngreitzLab/ENCODE-Distal-Regulation-Paper` (manuscript analysis code), `argschwind/ENCODE_CRISPR_data`, `EngreitzLab/ENCODE_Test_Dataset_Analysis` |
| Deposited Supplementary Tables | Tables 3, 8, 9, 11, 12 — fetched from Nature's MOESM3_ESM.zip via the EuropePMC `supplementaryFiles` REST API for PMC13471189 (`fetch_supplementary_tables.sh`; no auth, no proof-of-work gate — confirmed 2026-09-26) |
| Genome-wide predictions across all 1,458 biosamples | Authors' own HPC scratch (`/oak/stanford/groups/engreitz/...`), not deposited at that granularity; per-biosample summary stats *are* deposited (Supplementary Table 8) and used here instead |
| Primary CRISPR/DNase/Hi-C raw data | Public, ENCODE portal + GEO (see accessions below) |

## What IGVFagent does

Two of the paper's own toolchain pieces are already built into IGVFagent as first-class ports (not written for this benchmark): `igvfagent encode-re2g` (port of `EngreitzLab/ENCODE_rE2G`, including all 9 pretrained models) and `igvfagent abc-pipeline` (port of `broadinstitute/ABC-Enhancer-Gene-Prediction`).

* **Code correctness**: `igvfagent encode-re2g verify-upstream` re-scores the authors' own CircleCI-checked K562 chr22 test fixture with the embedded model coefficients and diffs against their checked-in expected output — max |score diff| 7.77e-16, all 9 model pickles bit-identical.
* **Paper-number recomputation**: `verify_derived_tables.py` recomputes the headline encyclopedia-scale, CRISPR-benchmark AUPRC, enhancer/gene-degree and MYC-synergy numbers directly from the paper's own deposited Supplementary Tables (rather than quoting its prose).

```bash
bash Benchmarks/gschwind2026_encyclopedia_enhancer/run.sh
python3 Benchmarks/concordance.py --benchmark gschwind2026_encyclopedia_enhancer
```

## Paper coverage

`igvfagent bench score` → **reproduction: reproduced, 5/5 analyses** (14/14 checks pass).

| Analysis | Paper | IGVFagent (recomputed from Supplementary Tables) | State |
|---|---|---|---|
| ENCODE-rE2G port vs authors' code (K562 chr22) | — | max \|score diff\| 7.77e-16; 9/9 models bit-identical | reproduced |
| Fig. 2b,c AUPRC: ENCODE-rE2G vs ABC (DNase-only) | 0.66 vs 0.56 | 0.6622 vs 0.5649 (Table 3) | reproduced (exact) |
| Fig. 2b,c AUPRC: ENCODE-rE2G_Extended vs ABC | 0.74 vs 0.61 | 0.7373 vs 0.6128 (Table 3) | reproduced (exact) |
| Encyclopedia: biosamples / cell types | 1,458 / 369 | 1,458 / 369 (Table 12) | reproduced (exact) |
| Encyclopedia: elements & interactions per biosample (mean) | 40,210 / 63,221 | 40,209.5 / 63,221.0 (Table 8) | reproduced (exact) |
| Encyclopedia: total interactions | 92,176,227 | 92,176,227 (Table 8) | reproduced (exact) |
| Extended Data Fig. 6b,c: genes/enhancer, enhancers/gene (median) | 1 / 3 | 1 / 3 (Table 8) | reproduced (exact) |
| Fig. 5 super-additive MYC enhancer pairs | "several lines of evidence" | 19/21 (90.5%) tested pairs significant, positive interaction (Table 11) | reproduced |

## Honest caveats

* **Extended Data Fig. 6's *mean* values (5.91 enhancers/gene, 1.57 genes/enhancer) are not independently reproduced.** Supplementary Table 8 only deposits per-biosample *medians* (which match exactly: 3 and 1). The means would require regenerating genome-wide ENCODE-rE2G predictions for all 1,458 biosamples — each biosample is itself a full ABC + rE2G run over real ENCODE DNase/H3K27ac/Hi-C data, which is far beyond a single benchmark run (see below).
* **Fig. 4b's specific "removing the nearby-enhancer-activity feature drops precision at 70% recall by 0.01" is not independently recomputed.** That delta is not in any deposited Supplementary Table (checked Tables 1, 2, 4; Table 2 defines the feature groupings, Table 4 only compares different published predictors, not feature-ablated variants) — it appears to live only in Fig. 4b's own (undeposited) source data. A background attempt to recompute it from scratch (download real K562 ENCODE DNase-seq `ENCSR000EOT` + Hi-C, run `abc-pipeline` + `encode-re2g crispr-features`/`train`/`feature-selection` end-to-end, retrain with/without `sumNearbyEnhancers`) was started but did not complete before the session's rate limit. What *is* independently reproduced instead is the paper's closely related Fig. 5 evidence for the same underlying phenomenon (enhancer-enhancer super-additivity), via the deposited Table 11.
* **The AUPRC and encyclopedia-scale numbers are recomputed from the paper's own deposited Supplementary Tables, not from an independent from-scratch reprocessing of raw ENCODE data.** This is a legitimate, arithmetic recomputation from primary deposited data (summing/averaging/counting across all 1,458 rows of Tables 8/12, reading Table 3's per-model AUPRC rows) — not a re-quote of the paper's prose — but it is not the same evidentiary strength as an independent re-run of the ABC/rE2G pipeline on raw signal tracks. `OPERATIONS.md` documents exactly what a from-scratch reprocessing would need (exact ENCFF accessions for K562 peaks/bam/Hi-C are known, via Supplementary Table 16) for anyone who wants to complete it.
* eQTL and GWAS benchmarking (the paper's other two benchmark axes) are not covered by this benchmark's 5 curated analyses, which focus on the CRISPR benchmark, the encyclopedia scale claims, and the enhancer-network/super-additivity claims (the paper's headline figures). IGVFagent does have `eqtl-enrich` and `gwas-e2g` ports of the corresponding upstream pipelines (`EngreitzLab/eQTLEnrichment`, `Deylab999MSKCC/e2g-benchmarking`) if that scope is wanted later.

## Provenance

`provenance.json` holds the full resolve/harvest/route record. `verify_derived_tables.py` documents exactly which Supplementary Table cell backs each number.
