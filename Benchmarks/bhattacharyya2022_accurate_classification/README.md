# bhattacharyya2022_accurate_classification

## Paper

**Accurate Classification of Cardiomyopathy Diagnosis by Chromatin Accessibility.**
Bhattacharyya S, Duan J, Vela RJ, Bhakta M, Bajona P, Mammen PPA, Hon GC, Munshi NV.
*Circulation* 2022 · doi:[10.1161/circulationaha.122.059659](https://doi.org/10.1161/circulationaha.122.059659) · PMID 36095061 · PMC9475804 · IGVF award IGVF0054, grant UM1HG011996

Resolver confidence: **1.00** (exact PMID match).

This is a short Circulation Research Letter (pp. 878–881, one composite Figure with panels A–K). ATAC-seq was performed on cardiomyocyte nuclei from 21 human heart specimens (6 healthy, 5 ischemic cardiomyopathy [ICM], 5 non-ischemic cardiomyopathy [NICM], 5 hypertrophic cardiomyopathy [HCM]), mapped to GRCh38. The paper shows UMAP separation of samples by etiology, calls 1066 differentially accessible regions (healthy vs. cardiomyopathy), derives subtype-specific TF motifs (MEF2/ICM, nuclear hormone receptors/NICM, ETS/HCM) and GREAT GO terms, trains a 1000-tree Random Forest classifier (LOOCV, 100 iterations) that reaches "excellent" AUC at both full and downsampled (20M-read) sequencing depth, and applies the classifier to three simulated real-world diagnostic cases.

## Why every analysis here is blocked, not reproduced

Europe PMC's `fullTextXML` endpoint returns HTTP 500 for this record (PMC9475804 is an NIH-funded Author Manuscript, not the publisher's own open-access XML), so `igvfagent bench harvest`'s automated fetch degraded to abstract-only. The full text was instead fetched directly from `pmc.ncbi.nlm.nih.gov` (confirmed 2026-09-27) and used to hand-build `harvest.json`.

The paper's own words on data: *"Anonymized data have been made publicly available at dbGaP (https://www.ncbi.nlm.nih.gov/gap/)."* That is the only data-availability statement in the text — a generic portal link, no specific `phs` study accession. Targeted NCBI E-utilities searches (`esearch db=gap` and `db=gds`, author names Bhattacharyya/Hon/Munshi crossed with cardiomyopathy/ATAC-seq, run 2026-09-27) found no matching dbGaP study or GEO series; the two dbGaP hits for "Hon AND cardiomyopathy" are unrelated large population cohorts (Framingham `phs000007`, CARDIA `phs000236`), not this study. The publisher page (`ahajournals.org`) returned HTTP 403 (bot-blocked), so a possible exact accession hiding in the paywalled Data Supplement could not be checked either.

There is also no Code Availability statement anywhere in the text — the only stated availability is *"Research materials, experimental procedures, and protocols are available from the corresponding authors upon reasonable request,"* i.e. no deposited code to run as a reference, and no GitHub/Zenodo repository is mentioned.

Every one of this paper's headline results (UMAP clustering, differential-accessibility calling, TF-motif/GREAT enrichment, the Random Forest classifier itself, and the downsampled case-study diagnoses) depends on the raw ATAC-seq reads from these 21 human cardiac specimens. Since that data sits behind dbGaP — an approved Data Access Request requiring institutional certification and IRB sign-off, a legitimate controlled-access gate, not a bot/CAPTCHA block — **none of the 5 planned analyses below can be run at all**, let alone verified against the paper's own numbers.

## Blocked analyses

| Analysis | Figure | Claim | Blocker |
|---|---|---|---|
| `fig_d_differential_accessibility` | Fig. D | 1066 regions differentially accessible, healthy vs. cardiomyopathy | `controlled_access` |
| `fig_c_g_umap_clustering` | Fig. C, G | UMAP resolves samples by etiology; refined peak subset improves it | `controlled_access` |
| `fig_e_f_subtype_motifs_go` | Fig. E, F | Subtype-specific accessible-peak heatmap, TF motifs (MEF2/NHR/ETS), GREAT GO enrichment | `controlled_access` |
| `fig_h_i_random_forest_auc` | Fig. H, I | Random Forest (1000 trees, LOOCV) — AUC robust to downsampling to 20M reads | `controlled_access` |
| `fig_j_case_studies` | Fig. J | 3 simulated diagnostic cases at 20M-read depth (external HCM, NICM pre/post-LVAD, discordant ICM) | `controlled_access` |

See `expected.json`'s `analyses[].blocker.reason` for the full evidence trail behind each block (identical text, since the same dbGaP/no-code situation blocks all five).

## Running it

```bash
bash Benchmarks/bhattacharyya2022_accurate_classification/run.sh
igvfagent bench score --paper-id bhattacharyya2022_accurate_classification
```

There is nothing to compute — `run.sh` only writes the run directory `concordance.py` expects. `igvfagent bench score` reports `reproduction: reproduced_except_access (0/5)`: every planned analysis needs controlled-access data, which is the honest, complete state for this paper, not a partial attempt.

## Coverage

`igvfagent bench report --paper-id bhattacharyya2022_accurate_classification` reproduces the coverage table above. 5 of 5 planned analyses are blocked (`controlled_access`); 0 are reproduced, because none can be attempted.

## Provenance

`provenance.json` holds the resolve record. `harvest.json` (in the timestamped `Docs/Benchmark/*_bhattacharyya2022_accurate_classification/` run directories) was hand-built from the PMC full text fetched directly, not by the automated harvester, because Europe PMC's `fullTextXML` 500'd on this record — see the harvest.json's own `fulltext_source` field for that note.
