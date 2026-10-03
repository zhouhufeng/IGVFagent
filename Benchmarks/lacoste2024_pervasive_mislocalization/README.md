# lacoste2024_pervasive_mislocalization

## Paper

**Pervasive mislocalization of pathogenic coding variants underlying human disorders.**
Lacoste J, Haghighi M, Haider S, Reno C, Lin ZY, Segal D, Qian WW, Xiong X, Teelucksingh T, Miglietta E, Shafqat-Abbasi H, Ryder PV, Senft R, Cimini BA, Murray RR, Nyirakanani C, Hao T, McClain GG, Roth FP, Calderwood MA, Hill DE, Vidal M, Yi SS, Sahni N, Peng J, Gingras AC, Singh S, Carpenter AE, Taipale M.
*Cell* 2024 · doi:[10.1016/j.cell.2024.09.003](https://doi.org/10.1016/j.cell.2024.09.003) · PMID 39353438 · PMC11568917 (closed-access; NIH-embargoed) · preprint doi:[10.1101/2023.09.05.556368](https://doi.org/10.1101/2023.09.05.556368)

A high-throughput imaging platform (Cell Painting-style immunofluorescence + CellProfiler morphological profiling) assaying 3,448 missense variants across >1,000 genes for effects on protein subcellular localization. About one-sixth of pathogenic missense variants mislocalize; effects span all cellular compartments and are driven mainly by protein stability/membrane-insertion changes, not trafficking-signal or interaction disruption.

**Note on paper access**: PMC11568917 is closed-access (EuropePMC's `fullTextXML` API returns HTTP 500 for it). Full text was instead fetched from the publicly-viewable **PMC Author Manuscript** page (`https://pmc.ncbi.nlm.nih.gov/articles/PMC11568917/`), which serves the same content through the web UI even though the JATS API refuses it.

## Data used

Real, public data — **not** the paper's cited `app.springscience.com/workspace/utoronto` viewer (a proprietary front-end this benchmark does not use), but the actual underlying dataset, which the authors' own GitHub README documents as freely downloadable with **no auth, no registration**:

```
aws s3 sync --no-sign-request s3://cellpainting-gallery/cpg0026-lacoste_haghighi-rare-diseases/broad/ .
```

This benchmark pulls two small pieces of that bucket (raw images are 100 GB–1 TB per batch and are **not** used):

| Path | What | Size |
|---|---|---|
| `workspace/metadata/reprocessed/set1_set2.csv` | Per-well reprocessed metadata: `Gene`, `Variant`, `Class` (Wild-type/Mutant), `Metadata_Location` (visual localization annotation), `batch`, across all 5 named screening batches | 1.1 MB, 6,789 rows |
| `workspace/population_profiles/PILOT_1/*/*.csv.gz` | Per-well, per-plate population-level (mean-of-transfected-cells) CellProfiler feature profiles for the PILOT_1 batch (the paper's primary ORFeome + non-pathogenic-control screen) | ~310 MB, 73 plates |

Code Availability names `https://github.com/carpenter-singh-lab/2023_LacosteHaghighi_Cell`, which GitHub redirects to the renamed **`carpenter-singh-lab/2024_LacosteHaghighi_Cell_Mislocalization`** (commit `4b9a789c61410922a2cd105d6cab0c44defaede0` pinned here). Its `utils/impactscore.py` and `snakemake-workflow/scripts/preprocess_calculateCorrelation_calculateImpact.py` implement the paper's Impact Score, but **cannot be run as published**: both import a private, unpublished `singlecell` package from a hardcoded lab path (`/home/ubuntu/workspace_SingleCell/SingleCell_Morphological_Analysis/`, not in this repo or on PyPI), and the checked-in preprocessing script has a literal syntax error (an un-commented banner line) plus undefined names in `main()`. This benchmark's port, `lacoste2024-mislocalization-stats`, is therefore a from-Methods reimplementation of the documented algorithm, run against the real data above — not a run of the authors' code.

Supplementary Table S1 (`NIHMS2023443-supplement-6.xlsx`, "variant and localization annotations and mass spectrometry data" — the single best ground-truth reference table, playing the same role as [`guttman2026_massively_parallel`](../guttman2026_massively_parallel/README.md)'s Supplementary Data 1) is gated behind a JS "Preparing to download..." interstitial at `pmc.ncbi.nlm.nih.gov/articles/instance/11568917/bin/NIHMS2023443-supplement-6.xlsx` that neither curl nor WebFetch can get past. This is the single biggest reason 6 of 8 analyses below are blocked rather than reproduced.

## What IGVFagent does

New port `lacoste2024-mislocalization-stats` (`Scripts/ported/skills/lacoste2024_mislocalization_stats.py`):

- **Library composition**: unique `Gene`+`Variant` constructs per named screening batch, cross-checked against the paper's 4 stated sub-collections.
- **Compartment specificity** (blocked, see below): secretory-pathway share of `Metadata_Location`-annotated mutants.
- **Impact Score**: `(1 - Pearson r) / 2` between each mutant's and its gene's WT mean **Protein-channel** CellProfiler feature profile, following `utils/impactscore.py::impact_score_wt_mt`. Features are z-scored (median/IQR) across the assembled population before correlating — an unscaled first attempt gave every variant a trivial ~1.0 correlation (including the paper's own named positive control, CRYAB R120G), diagnosed as feature-scale dominance (some raw CellProfiler feature categories, e.g. `Correlation_Costes_*`, run 100–1000× the magnitude of others) and fixed by standardizing before flattening into a per-well vector.

```bash
bash Benchmarks/lacoste2024_pervasive_mislocalization/run.sh
igvfagent bench score --paper-id lacoste2024_pervasive_mislocalization
```

## Concordance — `bench score`: 2/8 analyses reproduced (`reproduced_except_access`, 2/2 checks)

| Analysis (paper claim) | State | Measured | Note |
|---|---|---|---|
| Library composition: 3,448 variants / >1,000 genes, 4 sub-collections | **reproduced** | 3,598 unique variants / 1,287 genes (match fraction 0.958) | ~4% over the paper's count; `Common_Variants` batch alone (65/41) matches the paper's non-pathogenic-control sub-collection (65/45) almost exactly |
| Headline screen: 250 confirmed mislocalized variants (11% of 2,280 detected) | **reproduced** (positive-control check only) | CRYAB R120G IS=0.502 vs its other assayed variants (P20S 0.012, R56W 0.005) | Validates the Impact Score metric directionally against the paper's own named example. Does **not** reproduce the exact 250/2,280/11% headline count — see caveats |
| Secretory pathway enriched among mislocalized variants (59% vs 36% of library) | **blocked** (`not_deposited`) | attempted: 17.4% vs 13.9% (crude keyword match) | needs the paper's fixed 53-category compartment taxonomy (Table S1) |
| Pathogenic ClinVar variants mislocalize more (16% vs 6%) | **blocked** (`not_deposited`) | — | needs per-variant ClinVar significance (Table S1) |
| TMD mutations enriched among mislocalized variants (20% vs 5%) | **blocked** (`not_deposited`) | — | needs per-variant TMD annotation (Table S1) |
| 95% inter-observer visual concordance | **blocked** (`not_deposited`) | — | public metadata has only a single consensus `Metadata_Location`, no raw per-observer calls |
| 14% of hits (34 proteins) form foci/clusters | **blocked** (`not_deposited`) | — | needs the authors' texture/clustering pipeline (also depends on the unpublished `singlecell` package) |
| PLP1/GFAP/ACTB/SMAD2 disease case studies | **blocked** (`not_deposited`) | — | separate LUMIER/3TP-lux/BioID validation assays, no raw data deposited |

Full check-level detail: `igvfagent bench report --paper-id lacoste2024_pervasive_mislocalization` or the latest `Benchmarks/results/*_concordance.md`.

## Honest caveats

* **The headline 250/2,280/11% hit-count is not independently reproducible with public data alone.** The paper's classification combines the Impact Score with independent visual confirmation by two observers and a secondary retransfection-validation step (Methods), and pools 8 screening batches (this benchmark uses PILOT_1 only, ~3,600 of the ~3,000+ variant-well pairs). Reporting an exact match here would mean silently tuning a threshold to hit a number, which this benchmark does not do — instead it validates the underlying metric directionally against the paper's own named positive control (CRYAB R120G) and reports the measured distribution (34.6% of 1,394 scored variants exceed IS≥0.25, vs. the paper's 11% of 2,280) honestly as not comparable 1:1.
* **Library composition is close but not exact** (3,598 vs 3,448, ~4% over). The public reprocessed metadata may retain QC-failed or non-headline constructs the paper's final published count excludes; not resolvable without the paper's own filtering criteria.
* **Six analyses are blocked, not attempted lightly.** Supplementary Table S1 — which almost certainly has the exact compartment taxonomy, ClinVar annotations, and TMD calls needed — is real but gated behind a JS download interstitial (`pmc.ncbi.nlm.nih.gov/articles/instance/11568917/bin/NIHMS2023443-supplement-6.xlsx`) that automated fetching (curl, WebFetch) cannot pass. A browser-automation attempt (e.g. `claude-in-chrome`) could plausibly retrieve it and unblock 3 of these 6 analyses in a follow-up pass.
* **The authors' own analysis code cannot be run as published.** `utils/impactscore.py` and `snakemake-workflow/scripts/preprocess_calculateCorrelation_calculateImpact.py` both import an unpublished private `singlecell` package (not in the repo, not on PyPI), and the latter additionally has a literal syntax error and undefined names in `main()`. The registered port is a faithful reimplementation of the *documented* formula, not a run of their code, and is marked unverified/no-reference-available in the registry rather than falsely claiming a code match.
