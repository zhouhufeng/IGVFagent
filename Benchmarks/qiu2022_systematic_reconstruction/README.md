# qiu2022_systematic_reconstruction

**Systematic reconstruction of cellular trajectories across mouse embryogenesis.**
Qiu C, Cao J, Martin BK, Li T, Welsh IC, Srivatsan S, Huang X, Calderon D, Noble WS, Disteche CM, Murray SA, Spielmann M, Moens CB, Trapnell C, Shendure J.
*Nature Genetics* 2022 · doi:[10.1038/s41588-022-01018-x](https://doi.org/10.1038/s41588-022-01018-x) · PMID 35288709 · PMC8920898 · IGVF0030 · Award UM1HG011966

Also known as **TOME** (Trajectories Of Mammalian Embryogenesis) / MOCA-2. Resolver confidence: **1.00**.

## Data and code the paper's own Data/Code Availability names

| Repository | Accession | What it is |
|---|---|---|
| NCBI GEO | `GSE186069` | New E8.5 sci-RNA-seq3 data (~240k cells, 12 embryos); deposit includes a processed `cell_annotate.csv` + `gene_count.mtx` |
| NCBI GEO | `GSE186068` | Deeper re-sequencing of Cao et al. 2019's E9.5-E13.5 libraries |
| NCBI GEO | `GSE100597`, `GSE109071`, `GSE106587`, `GSE112294` | Earlier mouse scRNA-seq stages folded into TOME |
| ArrayExpress | `E-MTAB-6967` | Pijuan-Sala et al. E6.5-E8.5 gastrulation atlas |
| GitHub | [`ChengxiangQiu/tome_code`](https://github.com/ChengxiangQiu/tome_code) @ `0c7a72b0` | The paper's own analysis scripts (Sections 1-8) |
| Author website | [tome.gs.washington.edu](https://tome.gs.washington.edu) | Per-timepoint Seurat objects (mouse/zebrafish/frog), not on GEO |

`tome_code`'s `help_code/` folder also bundles the paper's own final Supplementary Tables (S6/S14/S18: mouse/zebrafish/frog key-TF calls), which this benchmark uses directly as ground truth.

## Coverage: `reproduced_except_access` — 4/7 analyses reproduced, 3/7 honestly blocked (`not_deposited`/`controlled_access`); 11/11 checks pass

| Analysis (`expected.json` id) | Figure | State | What was done |
|---|---|---|---|
| `fig1_e85_qc_annotation` | Fig. 1 | reproduced | New port `qiu2022_e85_qc_summary`: re-applies the Methods' exact QC filter (UMI>=200, genes>=100, unmatched_rate<0.4) + doublet/low-quality removal to GEO GSE186069's deposited per-cell table. Reproduces **239,533 cells** and **30 cell types** exactly. Class-B: new port `qiu2022_e85_scrublet_doublets` re-runs Scrublet (help_code/run_scrublet.py's exact parameters) on the deposited raw `gene_count.mtx` (49,585 genes x 239,533 cells) and correlates against the authors' own deposited `doublet_score` column: **Pearson r=0.878, Spearman rho=0.930** (all 239,533 cells), correctly recovering a ~1.1% predicted doublet rate at an automatically-set threshold of 0.32. |
| `fig5_keytf_nomination` | Fig. 5 | reproduced | New port `qiu2022_keytf_score`: re-implements Section6_keyTF_Step2's per-edge z-scoring ("determined_score") that combines child-vs-parent and child-vs-sibling differential-expression effect sizes. Verified against the authors' own Table S6 (mouse): **match_rate=1.0000, r=1.0000000000000002** (2783/2784 rows). Also reproduces the paper's own **632 key TFs / 92 cell types** exactly (counted directly from Table S6). |
| `fig6_7_cross_species_keytf` | Fig. 6-7 | reproduced | Same port applied unmodified to the paper's zebrafish table (Table S14): **match_rate=1.0000, r=1.0** (1901/1904 rows). (Table S18, frog, was also scored — r=0.9989 but exact match_rate=0.68; broadly agrees, not wired to a graded check; see below.) |
| `fig2_trajectory_reconstruction` | Fig. 2 | reproduced | New port `qiu2022_trajectory_connection` (+ companion R script `integrate_embed.R`, calling the authors' own `doClusterSeurat()` unmodified except `k.weight` lowered from Seurat's default 100 to 20, needed at this small scale): bootstrap KNN parent-state voting on the authors' own deposited per-timepoint Seurat objects (E3.5→E4.5, E5.25→E5.5). **7/7 child cell types reconstruct the correct parent state**; the 2 comparable Table S6 edge weights (0.925, 0.865) reproduce to 0.96 and 0.93. |
| `fig3_rna_velocity` | Fig. 3 | **blocked** (`not_deposited`) | Needs lineage-subsetted spliced/unspliced (exon/intron-split) count matrices for two specific lineages; GEO deposits only combined gene-level counts. Not released anywhere the paper points to. |
| `fig4_spatial_inference` | Fig. 4 | **blocked** (`controlled_access`) | Needs (1) a third-party paper's GEO-seq spatial data and (2) CIBERSORTx, a registered-account web-only tool whose terms forbid automated batch use. |
| `fig6_7_cross_species_coembedding` | Fig. 6-7 | **blocked** (`not_deposited`) | The full mouse/zebrafish/frog CCA coembedding + NNLS regression needs the zebrafish and frog whole-embryogenesis atlases, published by other groups under their own accessions -- not this paper's own deposits. |

Run with:
```bash
bash Benchmarks/qiu2022_systematic_reconstruction/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark qiu2022_systematic_reconstruction
```
`run.sh` covers Fig. 1's QC/annotation port and Fig. 5/6-7's key-TF port end to end (downloading GEO + the authors' own bundled Supplementary Tables itself). The Scrublet doublet-detection reproduction (2.2 GB raw count matrix; sbatch) and the trajectory-connection reproduction (Seurat CCA integration; sbatch-eligible, run here directly since each stage pair is small) are one-off analyses under `Benchmarks/_data/qiu2022_systematic_reconstruction/` -- see `OPERATIONS.md`.

## Honest caveats / where this differs from the paper

* **Fig. 1 doublet scores**: Scrublet simulates random doublets and uses approximate nearest neighbours (`annoy`); it is not seed-reproducible even across runs of the *same* script version, let alone across the authors' now-superseded `scrublet.compute_doublet_scores()` free-function API (pre-0.2) vs. the currently-installable `scrublet.Scrublet` class API used here (same package, same author, same algorithm -- the free function was removed in later releases). Scored as a correlation against the deposited `doublet_score` column, not an exact match.
* **Fig. 2 trajectory reconstruction** was run on two small, tractable stage pairs (E3.5→E4.5, E5.25→E5.5) using the authors' own deposited per-timepoint Seurat objects and their own `doClusterSeurat()`, rather than the full 19-stage, ~1.66-million-cell graph (`edge_all.rds`, an internal working file never deposited anywhere). `IntegrateData`'s `k.weight` had to be lowered from Seurat's default (100) to 20 because these early, low-cell-count stages don't have 100 anchor cells per batch-pair the way the full atlas does -- a numerical-feasibility change, not an algorithmic one. The KNN search (scikit-learn) and bootstrap RNG (NumPy) also differ from the authors' R (`FNN`/`sample()`), so results agree closely (7/7 correct top-parent calls; the 2 Table-S6-comparable edge weights within ~0.04-0.07 absolute) rather than exactly.
* **Fig. 6-7 frog key-TF table** (Table S18): the same, otherwise-exact port gives r=0.9989 but only 68% of rows match to atol=1e-6/rtol=1e-6 -- some frog-specific edges evidently used a slightly different population for the z-score than the one recoverable from the table alone. Reported honestly, not wired to a graded check.
* **12 prose-derived scaffold checks removed.** The initial `igvfagent bench scaffold` draft extracted 12 `[UNCONFIRMED]` numeric claims by regex from the full text (150k nuclei, 480 samples, 1,658,968 cells, 1000-gene HVG selection, etc.). Of these, the ones this benchmark's actual analyses can confirm (239,533 cells; 30, 92 cell types; 632 key TFs) were promoted into real checks above; the rest (sample counts, downsampling parameters, decrease-direction TF counts with no bundled supplementary table) were dropped as either Methods-only detail or out of this benchmark's chosen scope.

## Provenance

`provenance.json` holds the full resolve/harvest/route record. Ported code: `Scripts/ported/skills/qiu2022_e85_qc_summary.py`, `qiu2022_e85_scrublet_doublets.py`, `qiu2022_keytf_score.py`, `qiu2022_trajectory_connection.py` (registered in `Scripts/ported/registry.json`).
