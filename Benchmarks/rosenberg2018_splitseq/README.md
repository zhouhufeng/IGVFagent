# Rosenberg 2018 — SPLiT-seq, developing mouse brain and spinal cord

[![paper](https://img.shields.io/badge/Science-360:176--182-blue)](https://doi.org/10.1126/science.aam8999)
[![PMID](https://img.shields.io/badge/PMID-29545511-blue)](https://pubmed.ncbi.nlm.nih.gov/29545511/)
[![data](https://img.shields.io/badge/GEO-GSE110823-orange)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110823)
[![coverage](https://img.shields.io/badge/reproduction-7%2F12%20analyses-yellow)]()

## Bottom line

`igvfagent bench score --paper-id rosenberg2018_splitseq` → **reproduction: incomplete, 7/12 analyses** (16/22 checks).

Reproduced: the Fig. 1B species-mixing barnyard (a port of the authors' own code, matching it exactly), per-cell UMI and gene yields (Fig. 1C, cells), the fresh-vs-frozen correlation (Fig. 1D), agreement of a de novo clustering with the authors' 73 clusters (Fig. 2A) and 44 spinal clusters (Fig. 5A), the Allen ISH regional specificity of neuronal clusters (Fig. 3), and the size of the cerebellar interneuron lineage (Fig. 4E).

Not reproduced: per-nucleus UMIs (Fig. 1C, nuclei) and four cell-class proportions (Figs. 2B–D, 4C), which a de novo, marker-panel annotation does not recover within ±15%.

## Citation

Rosenberg AB, Roco CM, Muscat RA, Kuchina A, Sample P, Yao Z, Graybuck LT, Peeler DJ, Mukherjee S, Chen W, Pun SH, Sellers DL, Tasic B, Seelig G. **Single-cell profiling of the developing mouse brain and spinal cord with split-pool barcoding.** *Science* **360**:176–182 (2018). DOI [10.1126/science.aam8999](https://doi.org/10.1126/science.aam8999) · PMID 29545511 · PMC7643870

## Data and code

| Resource | Identifier |
|---|---|
| GEO series | `GSE110823` (6 samples, MATLAB v5 DGEs; the series links the authors' loading gist [Alex-Rosenberg/5ee8b14e…](https://gist.github.com/Alex-Rosenberg/5ee8b14ea580144facad9c2b87cebf10)) |
| CNS atlas | `GSM3017261_150000_CNS_nuclei` — 156,049 nuclei × 26,894 genes, with the authors' `cluster_assignment` (73) and `spinal_cluster_assignment` (44 + NA) |
| Species mixing | `GSM3017262`–`GSM3017265` (same-day / frozen cells + nuclei, human HEK293/HeLa-S3 + mouse NIH/3T3) |
| Authors' code | [Alex-Rosenberg/split-seq-pipeline](https://github.com/Alex-Rosenberg/split-seq-pipeline) @ `a711b56` (a fork of yjzhang/split-seq-pipeline; read processing + QC). The clustering/annotation code is described only in the Science supplement and is not public. |
| Allen DMBA | `api.brain-map.org` — Developing Mouse Brain Atlas sagittal ISH, P4 / P14 |

## Coverage

| Analysis | Paper | IGVFagent | State |
|---|---|---|---|
| Fig. 1B species mixing | 99.9% single-species, 0.1% collisions; purity 99.6% (human) / 99.0% (mouse) | 99.89%, 0.11%; 99.60% / 99.04%. Port `splitseq-barnyard-qc` = authors' `barnyard()` on 20/20 values | reproduced (A + B) |
| Fig. 1C cells | 15,365 UMIs / 5,498 genes (human), 12,243 / 4,497 (mouse); frozen cells 15,078 | 15,295 / 5,458; 12,122 / 4,436; 15,006 (all within 1.4%) | reproduced |
| Fig. 1C nuclei | fresh 12,113; frozen 13,636 median human UMIs | 15,610; 16,754 | **failing** (see below) |
| Fig. 1D fresh vs frozen | Pearson r 0.987 | 0.985 | reproduced |
| Fig. 2A clustering | 73 clusters | 76 Leiden clusters; AMI 0.68 vs the deposited 73 (resolved nuclei) | reproduced (AMI ≥ 0.5 is our threshold) |
| Fig. 2B non-neuronal | 27,096 nuclei (17.4%) | 36.7% | **failing** |
| Fig. 2C astrocytes | 50% of non-neuronal | 40.6% | **failing** |
| Fig. 2D OPC/oligo lineage | 10,087 nuclei (6.5%) | 10.9% | **failing** |
| Fig. 3 regional specificity | "most types" regionally specific (Allen ISH composites) | 28/37 clusters (75.7%) peak in their named region; chance 9% | reproduced |
| Fig. 4C CGC lineage | 15,360 nuclei (9.8%) | 6.5% | **failing** |
| Fig. 4E cerebellar interneurons | 1,890 nuclei (1.21%) | 1.10% | reproduced (fragile, see below) |
| Fig. 5A spinal re-clustering | 30 + 14 clusters | 48 Leiden clusters; AMI 0.80 vs the deposited 44 | reproduced (AMI ≥ 0.5 is our threshold) |

The deposited labels themselves reproduce every count the paper states exactly: 73 clusters; the seven OPC/oligodendrocyte clusters sum to 10,087; clusters 55–73 (non-neuronal) to 27,096; CGC clusters 25 + 28 to 15,360; cerebellar interneuron clusters 24, 26, 27, 29 to 1,890. These are identity checks on the authors' output, so they are not counted as coverage.

## Where the reproduction differs, and why

* **Per-nucleus UMIs (Fig. 1C).** The cell medians match within 0.5%, but the nuclei in the same GEO libraries give 15,610 (fresh) and 16,754 (frozen) median human UMIs, against the paper's 12,113 and 13,636. Using the same >90% species rule for cells and nuclei, no alternative (exonic-only, other thresholds) reaches the paper's values. A likely cause is that Fig. 1C was computed at a matched read depth or on a different processing of the nuclei, but the GEO DGEs don't show which.
* **Cell-class proportions (Figs. 2B–D, 4C).** The authors' clustering code is not public, so the atlas was re-clustered de novo (normalize → log1p → 3,000 HVGs → 50 PCs → kNN → Leiden) and clusters were assigned to classes by marker panels, never by the deposited labels. Two runs, both kept:
  * run 1 (`run1_res1_rawargmax/`): resolution 1.0 (31 clusters), raw panel-score argmax. This called almost everything neuronal (13.1% non-neuronal) because the neuron panel holds very abundant genes (Meg3, Snap25). The 621 microglia and 774 ependymal/OEC nuclei were all called neuron.
  * run 2 (scored): panel scores z-scored across clusters, and the smallest resolution giving ≥ 73 clusters (the paper's stated granularity; 5.0 → 76). This overcorrects (36.7% non-neuronal, 5.9% ependymal/choroid).

  Neither run recovers the proportions within ±15%. No further tuning was done, since that would fit the annotation to the answer. The Fig. 4E pass depends on this choice too: run 1 gave 2.8% (fail), run 2 gives 1.10% (pass).
* **Paper-internal inconsistencies.** "Neurons accounted for 79%" conflicts with "27,096 non-neuronal" out of 156,049, which gives 82.6%; the deposited labels agree with 27,096. The text gives the cerebellar interneuron lineage as 1,517 cells, while the Fig. 4E caption and the deposited labels give 1,890; the check uses 1,890.
* **Fig. 1D definition.** The paper does not state it. The definition was fixed before comparing: log1p CPM of human-cell pseudobulks, same-day 3000-UBC vs frozen 1000-UBC library. The same definition gives r = 0.981 for cells vs nuclei, against the paper's 0.952 (fig. S2, not planned).
* **Fig. 3 method.** For each deposited cluster whose name states a region (37 scorable), the top-5 enriched genes (largest log fold change, Wilcoxon adj. p < 0.05, detected in ≥ 10%) are averaged over Allen DMBA P4 (P2-dominated clusters) or P14 (P11-dominated) ISH expression energy, and the composite's peak region is compared with the name. Developing-atlas unionized rows stop at coarse structures, so the cerebellum is scored on rhombomere-1 alar plate (its anlage). The olfactory bulb has no unionized rows, so its 3 clusters are unscored. An earlier region map that used CbH (absent from the data) scored cerebellum 0/7; it is kept in `run1_regions_cbh/`.

## Ported into IGVFagent

* `igvfagent splitseq-barnyard-qc` — port of `split_seq/analysis.py` (`barnyard()` species calls + `generate_single_dge_report` per-species medians) from Alex-Rosenberg/split-seq-pipeline @ `a711b56`. `bench verify-port` against the authors' unmodified `barnyard()` on the same UBCs: match rate 1.0 over 20 values (`Docs/Benchmark/*_rosenberg2018_splitseq_port_splitseq_barnyard_qc/`). Unreviewed.

## How to reproduce

```bash
bash Benchmarks/rosenberg2018_splitseq/run.sh     # ~20 min; atlas steps need ~16 GB, run inside an allocation
igvfagent bench score --paper-id rosenberg2018_splitseq
```

Scripts: `fig1_metrics.py` (Fig. 1, via the port), `prep_atlas.py` (GEO atlas → AnnData with deposited labels), `atlas_recluster.py` (Figs. 2, 4, 5), `fig3_allen_regions.py` (Fig. 3). `make_figures.py` and the old figure belong to the earlier subsample benchmark and are no longer run.

## Earlier version of this benchmark

The previous README reported a 12,000-nucleus subsample recovering "8/8 CNS lineages". Its checks were class C counts (atlas size, number of Leiden clusters and annotated types), so they covered no paper analysis. It also wrote that the paper has ">100 fine clusters on full atlas". In fact the full-atlas clustering has 73 clusters; the >100 types are 69 brain plus 44 spinal types from two separate clusterings.

## License + provenance

* **Data**: GEO GSE110823 (public), Allen Developing Mouse Brain Atlas API (public); fetched at run time, never redistributed.
* **Code**: IGVFagent Apache-2.0; the authors' pipeline is used only as the pinned reference for `verify-port`.
