# southard2024_comprehensive_transcription

> This is the corrected, real reproduction for the paper the user originally
> asked for as "Joung 2025 TF Perturb-seq fibroblasts" via
> `https://www.nature.com/articles/s41588-025-02283-2`. That DOI is a 2-page
> Nature Genetics Research Briefing by Thomas M. Norman *summarizing* the
> paper below — it has no author named Joung. See
> [`../joung2025_tf_perturbseq/README.md`](../joung2025_tf_perturbseq/README.md)
> for the full correction (that directory previously carried a fabricated
> citation and fabricated accessions for this paper_id).

## Paper

**Comprehensive transcription factor perturbations recapitulate fibroblast transcriptional states**
Southard Kaden M.; Ardy Rico C.; Tang Anran; O'Sullivan Deirdre D.; Metzner Eli; Guruvayurappan Karthik; Norman Thomas M.
*Nature Genetics* **57**: 2323–2334 (2025) · doi:[10.1038/s41588-025-02284-1](https://doi.org/10.1038/s41588-025-02284-1)
(closed access; preprinted as bioRxiv [10.1101/2024.07.31.606073](https://doi.org/10.1101/2024.07.31.606073), open access, PMC11312553 — harvested from the preprint's full text since the published version isn't open access)

Resolver confidence: **1.00** (resolved).

## Real data and code (verified — see "How this was verified" below)

| Resource | Identifier | Status |
|---|---|---|
| Authors' own analysis code | [norman-lab-msk/TFs_CRISPRa](https://github.com/norman-lab-msk/TFs_CRISPRa) | public, pinned at commit `3637f77` under `Data/PaperCode/norman-lab-msk__TFs_CRISPRa/` by `igvfagent paper-code fetch` |
| Raw sequencing reads (SRA) | [PRJNA1108254](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1108254) | public |
| Processed Hs27 fibroblast dataset | Zenodo [10.5281/zenodo.15200179](https://doi.org/10.5281/zenodo.15200179) | public (4 files, 58MB–9.7GB) |
| Processed RPE-1 dataset | Zenodo [10.5281/zenodo.15213619](https://doi.org/10.5281/zenodo.15213619) | public (4 files, up to 29.7GB) |
| Cellranger raw outputs (Hs27, RPE-1) | Zenodo [15213597](https://doi.org/10.5281/zenodo.15213597), [15211972](https://doi.org/10.5281/zenodo.15211972) | public (21.6GB, 46.8GB) |
| Bulk RNA/ATAC/CUT&RUN characterization | Zenodo [10.5281/zenodo.15215216](https://doi.org/10.5281/zenodo.15215216) | public (bulk RNA-seq TPM/counts + ATAC/CUT&RUN bigwigs, Hs27 + RPE-1) |
| **Raw (pre-regression) final-population objects** | Zenodo, inside the same two records above: `fibroblast_CRISPRa_final_pop.h5ad` (9.7GB) and `RPE1_CRISPRa_final_population.h5ad` (29.7GB) | public — the direct input to the authors' own Step 3/4 regression notebooks (not previously catalogued here) |
| **RPE-1 essentials benchmark (Fig 1)** | Zenodo [10.5281/zenodo.15215414](https://doi.org/10.5281/zenodo.15215414) `RPE1_E150_all_genes.h5ad` (2.6GB) + [15215389](https://doi.org/10.5281/zenodo.15215389) (cellranger outputs) | public — not previously catalogued here; matches the internal `20240227_RPE1_E150_all_genes.hdf` `fig1_overloading_normalization` reads |
| IGVF-format resubmission (in progress) | [norman-lab-msk/igvf-perturbseq](https://github.com/norman-lab-msk/igvf-perturbseq) | code public; no live `IGVFDS…` accession found via Portal search as of this run |

The full [Zenodo `normanlabmsk` community](https://zenodo.org/api/records?communities=normanlabmsk) has 11 records; the 2 above were found via a direct community-listing query (`zenodo.org/api/records?communities=normanlabmsk`) after the GitHub README's own text undersold what's deposited. 2 further community records (`pEM040 sgRNA singlets`, `multiomeperturbseq`) look like a different sub-project from the same lab and are not referenced by any of this paper's own analysis notebooks — not pursued here.

The in-text accessions IGVFagent's harvester found (`GSE186458`, `MiraldiLab/maxATAC_data`, `lmcinnes/enstop`, `moshi4/pyCirclize`, `zenodo.7598955`) are all **third-party tools and reference datasets the paper's methods cite**, not the paper's own data — this preprint's body has no explicit Data/Code Availability section (unusual, but genuine: verified against every section heading in the harvested full text). The real accessions above were found by reading the authors' own GitHub README and the Zenodo API directly, not by the automated harvester.

## What IGVFagent does

**Route:** `perturb_catalog` — Perturbation Catalogue census (CRISPR screen / Perturb-seq), plus a real data-verification step (below)
**Skill output dir:** `Docs/Perturbation/`

```bash
bash Benchmarks/southard2024_comprehensive_transcription/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark southard2024_comprehensive_transcription
```

## How this was verified — real data, not text-only

The paper's own analysis is ~50 Jupyter notebooks operating on raw 10x
cellranger output (up to 46.8GB per cell line) — not realistically
re-runnable end-to-end in an interactive session. Instead,
`verify_guide_library.py` downloads the smallest real processed artefact
(the 58MB guide×droplet UMI matrix from Zenodo 10.5281/zenodo.15200179) and
independently counts its guide-identity axis:

| Quantity | Paper claim (prose) | Measured directly from the authors' deposited matrix | Match |
|---|---:|---:|:---:|
| Guides in final library | 10,979 | **10,979** | ✓ exact |
| Non-targeting control guides | 78 | **78** | ✓ exact |
| Distinct TFs targeted | 1,836 | **1,836** | ✓ exact |

This is a genuine measurement against primary data (not a repeat of the
paper's own text), and it is what two of `expected.json`'s `confirmed: true`
checks are based on.

`verify_activation_counts.py` goes further: it downloads both cell types'
per-guide summary files (`mean_pop.h5ad`, Hs27 1.7GB + RPE-1 1.7GB — the
authors' own per-guide/per-target activation calls, not raw counts) and
pools "target gene activatable in >=1 cell type" across them, using
`obs['expressed'] | obs['active'] | obs['masked_active']` as the per-guide
activation flag — grounded in the authors' *own* activation-count formula
(`is_activated | masked_active | expanded_masked_active`, cell 70 of
`Data/PaperCode/norman-lab-msk__TFs_CRISPRa/src/Code/Analysis of CRISPRa on
target activation/On target activatation and maximal target variation
across cell types.ipynb`) rather than an independently-chosen proxy:

| Quantity | Paper claim | Measured from primary data | Match |
|---|---:|---:|:---:|
| Guides in final library | 10,979 | **10,979** | ✓ exact |
| Non-targeting control guides | 78 | **78** | ✓ exact |
| Distinct TFs targeted | 1,836 | **1,836** | ✓ exact |
| TFs activatable in ≥1 cell type | 1,482 | **1,514** | ✓ within tolerance |
| Genes resistant to activation entirely | 319 | **323** | ✓ within tolerance |

An earlier pass used `obs['expressed']` alone (1,438 / 399 — the 399 fell
outside the paper's tolerance band and was left `[UNCONFIRMED]`), on the
assumption that `active`/`masked_active` were purely guide-clustering QC
flags for the seed-off-target analysis. Reading the authors' own notebook
(above) shows that assumption was wrong: they explicitly OR `masked_active`
into their own "activated" tally. Adding it (`active` stands in for the
`is_activated` regression test and `expanded_masked_active` reclustering,
neither of which is present in this deposited summary file) brings both
figures inside tolerance. This is **still not bit-identical** to the
authors' own per-cell-type gate — see Honest caveats — so both checks are
tied to `fig2c_ontarget_activation` in `expected.json` but that analysis
scores as class C (retrieval/enumeration: a count), not class A/B
(quantitative reproduction), until the actual regression test is ported.

## Raw-data reproduction (Step 1: guide calling from scratch)

Everything above measures the authors' own *already-processed* deposits.
`reproduce_step1_aggregation.py` goes one level deeper: it is a port of the
authors' own `Step 1 - Aggregation of cellranger outputs and guide
thresholding for Hs27 experiment.ipynb`, run against the **raw** cellranger
outputs (Zenodo 10.5281/zenodo.15213597, 21.6GB — the 16 individual
`cellranger count` lanes, downloaded and processed independently of the
`mean_pop.h5ad`/guide-UMI summaries used above).

Two real things fell out of actually doing this:

1. **A bug in the paper's own public code.** `Code/perturbseq/__init__.py`
   imports `util_jmr.py` and `aneuploidy.py`, neither of which exists
   anywhere in the repo at the pinned commit (`3637f77`) — verified directly
   via GitHub's contents API, not a fetch artefact on this end. Patched out
   (both are unused by the Step 1 notebook) rather than silently worked
   around; see the comment left in the patched file.
2. **Deviation, stated plainly:** the deposit has 16 per-lane `cellranger
   count` outputs, not the single merged `cellranger aggr` output the
   notebook's `EXPERIMENT` path expects — no `cellranger` binary is
   available here to redo that merge. The script concatenates the 16 lanes
   itself (barcode-suffixed per lane, the same convention `cellranger aggr`
   uses) with no depth-equalisation, so absolute per-cell UMI counts are not
   expected to be bit-identical to a real `cellranger aggr` run.

Despite that deviation, re-deriving the guide library from scratch —
517,319 raw droplets → guide/GEX feature split → per-cell guide-UMI
thresholding (>5) → dominant-guide assignment — reproduces a guide
library **identical, set-for-set, to the separately-deposited processed
summary**:

| Quantity | From raw cellranger (this run) | From the deposited processed summary | Match |
|---|---:|---:|:---:|
| Guides in final library | **10,979** | 10,979 | ✓ exact |
| Guide identity set | — | — | ✓ Jaccard = 1.0 (bit-for-bit identical) |
| Cells assigned a guide (pre singlet-filter) | 497,004 | — (this file predates singlet filtering) | n/a |

This is the strongest evidence in this benchmark that the pipeline holds
together end-to-end: two *independent* Zenodo releases of the same
underlying experiment (a raw cellranger dump and a processed summary),
reduced by two independently-written pieces of code (the authors' own
notebook, ported here, vs. whatever produced the deposited summary),
land on the exact same 10,979-guide library.

## Concordance

**6 / 6 confirmed checks pass** (guides, non-targeting controls, TF count, pooled-activatable count, resistant count — all measured from the authors' own deposited data, not from prose). **8 further paper-claimed numbers remain `[UNCONFIRMED]`** (perturbation count, cluster counts, AUC values). Of the 10 headline analyses `expected.json` now tracks (`analyses[]`), only `fig2c_ontarget_activation` has a check tied to it, and that check is class C, not A/B — see "Reproduction status" below. See `Docs/Benchmark/*_southard2024_comprehensive_transcription/replication_report.md` for the full table.

## Reproduction status (`analyses[]`)

`expected.json` tracks the paper's 10 headline computational results (one per major figure), each pointing at its upstream notebook(s) under `Data/PaperCode/norman-lab-msk__TFs_CRISPRa/src/Code`. `igvfagent bench plan --paper-id southard2024_comprehensive_transcription` lists them; `bench score`'s "paper coverage" only counts an analysis as reproduced once a **class A/B** check is tied to it via `"analysis": "<id>"`.

| Analysis | Status | What it would take |
|---|---|---|
| `fig2c_ontarget_activation` | **weak** (checks tied, passing, class C) | Promote to class A/B: port the actual Step 3/4 least-squares regression (`is_activated` = FDR-tested coefficient sign) against the now-public `fibroblast_CRISPRa_final_pop.h5ad` (9.7GB) / `RPE1_CRISPRa_final_population.h5ad` (29.7GB) instead of the `expressed\|active\|masked_active` proxy. Large (~11k regressors x ~20k genes x ~1M cells lstsq) — likely an sbatch job, not interactive. |
| `fig4a_tf_hits` | pending | Same regression as above unlocks `masked_active`/TF-hit counts directly (already in the deposited `mean_pop.h5ad`); `expanded_masked_active` needs a further reclustering step on top. Shares the fig2c regression work. |
| `fig4_sparse_pca_programs` | pending | Its only input (`Step 1` notebook) is `mean_pop.hdf` — **already downloaded** (`fibroblast_CRISPRa_mean_pop.h5ad`). No new data needed; port the nonnegative sparse PCA + bootstrap-resample model selection (Steps 1–2). Moderate effort, no new download. |
| `fig1_overloading_normalization` | pending | Needs `RPE1_E150_all_genes.h5ad` (2.6GB, Zenodo 15215414 — newly catalogued above, not yet downloaded) + port of the Step 3 normalization/benchmarking regression. Moderate effort, one new download. |
| `fig6ab_causal_drivers` | pending | Needs the public Human Protein Atlas bulk tissue file (`rna_single_cell_type_tissue.tsv`, external to the paper, downloadable from proteinatlas.org) + the repo's own `permutation-enrichment-main` package. Moderate-large effort. |
| `fig4b_cross_context_conservation` | pending | Needs per-gene ATAC accessibility scores; the paper's bulk ATAC is public only as whole-genome bigwig tracks (Zenodo 15215216) — deriving a per-gene score from those is its own small pipeline (peak/promoter quantification), not a straight port. Also needs bulk RNA-seq (public, same record) for its GO-term/cross-cell-type comparison. Large effort. |
| `fig5de_newly_expressed_markers` | pending | Same ATAC-gene-score gap as fig4b, plus an HPA retina reference file (`retina_norm_df_filtered_90pctl_celltype_scores.csv`) that is itself a derived output of `fig6ab`'s notebook, not raw data. Large effort, depends on fig4b + fig6ab groundwork. |
| `fig3a_seed_clusters` | pending | Needs internal-only pickled clustering intermediates (`*_masked_variables.pickle`) not found in any of the 11 Zenodo records — would need reimplementing the clustering from scratch off the raw guide-embedding data. Large effort. |
| `fig3s5_seed_regression_fitness` | pending | Needs external "published genome-wide CRISPRa/CRISPRi fitness screens" this repo doesn't itself deposit (likely a different paper's screen data, e.g. Replogle/Horlbeck-lineage) — the accession isn't identified yet. Needs its own harvest before any porting. |
| `fig2g_epigenetic_model` | pending | **Blocked at the source**: its cited upstream folder (`Epigenetic determinants of susceptibility to CRISPRa/Models for on target activation from epigentic features`) is genuinely empty in the fetched repo (verified commit `3637f77`) — no notebook to port. Would have to be reimplemented from the Methods text's feature description alone, with no reference implementation to check against. |

Two of these (`fig4_sparse_pca_programs`, needing only already-downloaded data, and `fig1_overloading_normalization`/`fig6ab_causal_drivers`, needing one bounded new download each) are the next-lowest-effort targets. The Step 3/4 regression that would fully close `fig2c_ontarget_activation` (and partially unlock `fig4a_tf_hits`) is high-value but computationally heavy (multi-GB dense least-squares) and warrants an explicit sbatch job rather than an interactive run.

## Honest caveats

* **Run directory is not paper-tagged.** This route's catalogue-census CLI steps accept no `--label`, so `concordance.py` matches the skill's default `summary` directory. All real checks above are pinned to paper-specific artefacts (via `extra_search_dirs`) regardless.
* **GSE237056 does not exist for this paper (or at all).** Verified via NCBI E-utilities (`esearch`): zero hits, and no GEO series is linked to this paper's SRA BioProject (`elink` returns no `gds` linkset). This paper deposited raw reads to SRA only and processed data to Zenodo only — it never used GEO. (This was the original fabricated benchmark's central false claim.)
* **The 319/1,482 figures are matched, not reproduced bit-for-bit.** `active` is used as a stand-in for the authors' `expanded_masked_active` (a further Hs27-only reclustering not in the deposited summary) and their `is_activated` regression test isn't in the summary file at all. The match is real and grounded in the authors' own formula, but promoting `fig2c_ontarget_activation` past class C means porting that regression.
* **8 unconfirmed checks remain** (perturbation count, cluster counts, AUC values) — would need the authors' regression/clustering notebooks (Steps 3-4, sparse PCA, epigenetic model), not yet ported. The raw-data pipeline itself is no longer untouched, though: see "Raw-data reproduction" above — Step 1 (guide calling from 16 real cellranger lanes) is done and cross-validated exactly against the deposited summary; Steps 2-4 (cell-population construction, masked-active regression, on-target/AUC modelling) are the multi-hour/TB-scale remainder.
* **No IGVF Portal accession found.** The Norman lab's `igvf-perturbseq` repo says this data was "formatted for IGVF submission," but Portal search found no matching `MeasurementSet` as of this run — likely still in progress upstream.

## Provenance

`provenance.json` in this directory holds the full resolve / harvest / route record, including every source URL consulted. The repo fetch record (pinned commit, fetch timestamp) is in `Data/PaperCode/norman-lab-msk__TFs_CRISPRa/` (created by `igvfagent paper-code fetch`, outside version control).
