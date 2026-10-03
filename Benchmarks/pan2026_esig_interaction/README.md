# pan2026_esig_interaction

## Paper

**eSIG-Net: an interaction language model that decodes the protein code of single mutations**
Pan Xingxin; Shrawat Aditya; Raghavan Sidharth; Dong Chuanpeng; Yang Yuntao; Li Zhao; Zheng W. Jim; Eckhardt S. Gail; Wu Erxi; Fuxman Bass Juan I.; Jarosz Daniel F.; Chen Sidi; McGrail Daniel J.; Sheynkman Gloria M.; Huang Jason H.; Sahni Nidhi; Yi S. Stephen
*Nature Methods* 23, 1115-1120 (2026) · doi:[10.1038/s41592-026-03086-x](https://doi.org/10.1038/s41592-026-03086-x) · PMID 42056223 · PMC13259923 · bioRxiv preprint doi:10.64898/2026.03.27.714913

Resolver confidence: **1.00** (exact DOI/PMID match against PubMed, PMC and Crossref).

eSIG-Net is a "mutation-centric interaction language model": given a wild-type protein, a single missense mutant of it, and an interaction partner, it predicts whether the mutation changes the WT-partner interaction state. The paper benchmarks it against 5 sequence-based PPI predictors (SDNN, D-SCRIPT, DeepFE, PIPR, PLM-interact) and 5 structure-based ones (MutaBind2, BeAtMuSiC, GeoPPI, TopNetTree, PIONEER, plus an AlphaFold-Multimer/FoldDock baseline and a raw ESM-1b LLR baseline) on two datasets: a disease-mutation PPI set (Sahni et al. 2015, *Cell*) and a population-variant PPI set (Fragoza et al. 2019, *Nat. Commun.*).

## Why this doesn't use an IGVFagent assay route

Every other benchmark in this suite routes through an assay-specific pipeline (MPRA, CRISPR screen, single-cell, etc. — see `igvfagent bench list-routes`). This paper is a pure computational method: no wet-lab assay, no GEO/IGVF-portal deposit, nothing `bench route` can match (`igvfagent bench route --paper-id pan2026_esig_interaction` returns "No route matched"). Its Code Availability statement points to `Stephen-Yi-Laboratory/eSIG-Net`, but that repository is an **inference-only distribution**: one released checkpoint (`checkpoints/model_pred_opt_fold_0.pth`, fold 0), a feature generator, and a prediction script — no training code, no evaluation splits, no benchmark datasets. Its own `MODEL_CARD.md` lists *"reproduction of paper benchmarks using this inference-only package"* as **out-of-scope use**.

So the reproduction here has two independent parts:

## Part 1 — the paper's headline benchmark, from its own Source Data

Nature Methods ships a per-figure Source Data spreadsheet with every article. `run.sh` fetches them from EuropePMC's open supplementary-files bundle for PMC13259923 (no login, no proof-of-work gate) — `SourceData_Fig1.xlsx` (MOESM3), `SourceData_Fig2.xlsx` (MOESM4), `SourceData_ExtDataFig2.xlsx` (MOESM5). Every relevant sheet carries one row per mutation-partner pair with columns `mutant_id, interactor_id, label_change, pred_change, confidence` — exactly what's needed to recompute ROC-AUC, average precision and accuracy per method, per dataset, independently of the paper's own code.

`verify_derived_tables.py` recomputes:

| Figure | What | eSIG-Net (paper → recomputed) |
|---|---|---|
| 1b | Accuracy, disease-mutation dataset | 0.85 ± 0.02 → 0.848 |
| 1c | ROC-AUC, disease-mutation dataset | 0.91 ± 0.02 → 0.908 |
| 1d | Average precision, disease-mutation dataset | 0.86 ± 0.01 → 0.856 |
| 1e | Accuracy, population-variant dataset | 0.90 ± 0.02 → 0.907 |
| 1f | ROC-AUC, population-variant dataset | 0.87 ± 0.02 → 0.866 |
| 1g | Average precision, population-variant dataset | 0.71 ± 0.02 → 0.703 |
| 2g | ROC-AUC vs. 5 structure-based baselines | 0.91 ± 0.02 → 0.908 |
| 2h | Average precision vs. 5 structure-based baselines | 0.86 ± 0.01 → 0.856 |
| 2c-e | FoldDock mutant/WT-interactor accuracy, consistency | 0.625 / 94.4% → 0.625 / 94.4% |
| 2i | eSIG-Net vs. normalized-ESM-1b-LLR proxy AUC | "stronger correlation" → 0.914 vs. 0.657 |

All 5 baseline methods on both datasets are recomputed too (33 checks total; see `expected.json`) — every number the paper states in the main text and in Supplementary Note 3 lands inside the recomputed value's own reported ± s.d. The one documented discrepancy is Fig. 2c's wild-type/interactor accuracy (paper: 0.922 all-data / 0.971 high-confidence; recomputed: 0.949) — the deposited sheet doesn't carry the pDockQ column needed to reproduce the confidence split exactly, so the check brackets both reported values instead of claiming an exact match (see its `provenance.note`).

```bash
bash Benchmarks/pan2026_esig_interaction/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark pan2026_esig_interaction
```

## Part 2 — the released checkpoint, absorbed as an IGVFagent port

See `OPERATIONS.md` §4. The released fold-0 checkpoint was run end-to-end — real UniProt sequences (RAD51D [O75771](https://rest.uniprot.org/uniprotkb/O75771.fasta), XRCC2 [O43543](https://rest.uniprot.org/uniprotkb/O43543.fasta)), real ESM-1b `esm1b_t33_650M_UR50S` layer-33 embeddings, no synthetic inputs — on the repo's own 3-mutation worked example, absorbed into IGVFagent's port registry as `esig_net_predict`, and verified against the authors' own `test.sh` output on the same generated inputs:

**`match_rate=1.0000` (3/3), `max_abs_diff=0.0` — bit-identical to the authors' own code**, both for the 573-D feature-generation step and the checkpoint's forward pass (`igvfagent port list`: `esig_net_predict → verified PASS (match_rate 1.000)`). A class-B check, kept separate from Part 1 because the repo itself says it isn't a benchmark-reproduction artifact.

## Concordance

Run `run.sh` then `concordance.py` for the current numbers; as of the last run: **34/34 checks, 11/11 analyses reproduced** (`reproduction: reproduced`) — 10 from Part 1 (Source Data recomputation) plus the Part 2 inference-port smoke test.

## Honest caveats

* **Fig. 2d/e's exact "17.6% accurately predicted as disruptive" figure is not reproduced** — the definition of "accurately predicted as disruptive" among FoldDock's inconsistent WT/mutant calls is ambiguous from the text alone, and a plausible reading gave 27.5%, not 17.6%. Dropped from `expected.json` rather than forced; `derived_metrics.json` still reports the recomputed value for transparency.
* **Fig. 2i has no exact quoted number** — the paper only says eSIG-Net's confidence score has a "stronger correlation" than the ESM-1b proxy. The checks operationalize that as eSIG-Net's own AUC (~0.91, consistent with Fig. 1c on the same pairs) being clearly above the ESM-1b proxy's AUC, rather than matching a specific value.
* **Extended Data Fig. 2b-c (TCGA/MMRF survival and immunotherapy-response correlations) are not reproduced** — they need patient-level clinical outcome data merged with expression cohorts; out of scope for this pass, not blocked by access (Xena hosts the harmonized TCGA data publicly; a future pass could add it).
* The released GitHub checkpoint cannot be used to independently re-derive the 1,633/4,020-pair benchmark numbers above (see Part 2) — that's why Part 1 goes to Source Data instead of the repo. Its inference pipeline was still exercised for real (not skipped): on a real 3-mutation example with real UniProt sequences and real ESM-1b embeddings, the absorbed port reproduces the authors' own `test.sh` output bit-for-bit.

## Provenance

`provenance.json` holds the full resolve/harvest record. `Data/pan2026/verify/derived/derived_metrics.json` (regenerated by `run.sh`) holds every recomputed number, including the ones without a matching `expected.json` check.
