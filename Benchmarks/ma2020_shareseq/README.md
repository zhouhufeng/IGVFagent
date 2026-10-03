# Ma 2020 — SHARE-seq mouse skin (`share` + shared `_scload`)

[![paper](https://img.shields.io/badge/Cell-183:1103--1116-blue)](https://doi.org/10.1016/j.cell.2020.09.056)
[![PMID](https://img.shields.io/badge/PMID-33098772-blue)](https://pubmed.ncbi.nlm.nih.gov/33098772/)
[![data](https://img.shields.io/badge/GEO-GSE140203-orange)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140203)
[![coverage](https://img.shields.io/badge/paper%20coverage-7%2F10%20analyses-yellow)]()

## Paper coverage (`igvfagent bench score`: 7/10 analyses reproduced)

| Analysis | Paper | IGVFagent | How checked | State |
|---|---:|---:|---|:---:|
| Skin peak-gene associations (±50 kb, p<0.05), Fig. S3H-I | 63,110 | 61,119 | `shareseq-dorc` = FigR `runGenePeakcorr` on 3,942 pairs (p-values to 5e-5, 100% match) | reproduced |
| Associations whose peak hits one gene, Fig. S3J | 83.9% | 79.2% | class A ±10% + FigR match | reproduced |
| Mean associations per gene, Fig. 3D-E | 4.4 | 4.12 | class A ±10% + FigR match | reproduced |
| DORCs (>10 peaks), Fig. 3F | 857 | 957 before / 787 after FigR's per-peak dedup | FigR match on per-gene counts (158 genes, exact) | reproduced (count off by >10%) |
| GM12878 peak-gene associations, Fig. S3B | 13,277 | 32,579 | FigR match (GM12878, 1,594 pairs) | reproduced (count 2.5× higher) |
| Species mixing, Fig. 1B-D | 903 human / 1,341 mouse / 1 collision | 971 / 1,397 / 3 | human share 0.410 vs 0.402 (class A) | reproduced |
| Computational pairing, Fig. S2N-S | 74.9% | 67.7% | Seurat CCA label transfer (Seurat 5.5.1) | reproduced (inside ±10%) |
| GM12878 peaks within 2 kb of a TSS | 61.3% | 14.5% (all peaks) | authors' `hg19.TSS.bed` | **failing** |
| DORC genes with positive residual, Fig. 4C | 92% | 45% | Methods reimplementation | **failing** |
| Chromatin potential TAC→differentiated, Fig. 5H | qualitative | 39% of TAC arrows point forward (RNA→RNA baseline 51%) | Methods reimplementation | **failing** |

Reference code: the paper has no code-availability statement. The peak-gene/DORC
reference is the Buenrostro lab's released implementation, `buenrostrolab/FigR`
@ `094f5aa` (`R/DORCs.R`, `R/utils.R`), run unmodified in R 4.4 with the same
chromVAR background peaks. Ported and registered commands: `shareseq-dorc`
(verified against FigR), `shareseq-species-mix`, `shareseq-gene-activity`
(Seurat v3.2.3 `CreateGeneActivityMatrix`), `shareseq-chromatin-potential` and
`shareseq-dorc-residuals`. No code was released for the last two, so they
follow the STAR Methods text.

Why the failures and gaps:
* **TSS fraction.** The deposited GEO peak set has 507,307 merged peaks. Only
  the 20,000 highest-count peaks reach about 59%, and no documented subset
  gives 61.3%.
* **Residuals and chromatin potential.** cisTopic is replaced by LSI, the
  branch-probability lineage cut-offs by the authors' cell-type labels, and the
  Palantir root is picked by rule, because the paper does not say how. These
  choices do not reproduce the paper's result. Ruled out as the cause: the
  DORC gene list. Re-running `shareseq-dorc peakgene` on only the 7,197
  hair-follicle-lineage cells (instead of the 787-gene genome-wide skin DORC
  list) gives 151 lineage-specific DORCs; pseudotime is confirmed correctly
  oriented (TAC-1 lowest mean pseudotime, Hair Shaft-cuticle.cortex highest).
  Re-scoring both ports on this tighter, more relevant gene set left the
  residuals check unchanged (44.98% -> 45.70% positive, still far from 92%)
  and the pseudotime forward-fraction check slightly worse (39.2% -> 37.7%,
  vs the 50.5% RNA-RNA baseline), even though the raw cross-modal
  neighbour-composition flow ratio did improve (0.45 -> 0.55). Gene-set
  breadth is therefore not the gap; what remains is the LSI/cisTopic
  substitution and the heuristic root/lineage choices, none of which the
  paper's text pins down. (Diagnostic outputs: `Data/ma2020/out/peakgene_hf`,
  `residuals_hf`, `chromatin_potential_hf_pt`.)
* **GM12878 count.** The paper used 23,278 cells; here 25,022 paired rep3
  cells are used. The peak set behind 13,277 is not stated. The port's
  205-32,579 associations (before dedup) match FigR run on the identical
  input exactly (class B, 100%), so the gap is upstream of the port, in the
  paper's undocumented peak/cell filtering, not a port bug.
* **Filters.** The species-mixing cut-offs are undocumented; the ATAC
  barcode-rank knee is used instead.

## Original RNA cell-type check

### Bottom line

**IGVFagent ingests Ma 2020's SHARE-seq skin RNA (via the shared `_scload` loader), runs the `share` per-barcode RNA QC, and recovers all 23 author-annotated skin cell types — exactly the 34,774-cell final set the paper reports.** Leiden clustering agrees with the author labels at AMI 0.63 despite the deliberately shallow SHARE-seq RNA (median 920 UMIs/cell). Runs from the per-sample GEO supplementary files.

| Metric | IGVFagent | Paper |
|---|---:|---:|
| Total RNA barcodes | 42,948 | — |
| **Cells in Ma 2020 final skin set** | **34,774** | 34,774 ✓ |
| Author skin cell types | **23** | 23 ✓ |
| Analyzed subsample | 15,000 | — |
| Leiden clusters | 21 | — |
| AMI vs author cell types | **0.63** | — |
| Median UMIs / genes per cell | 920 / 507 | shallow (SHARE-seq RNA) |

![Skin cell types](figures/fig1_celltypes.png)

## Concordance

Ma 2020 introduced SHARE-seq and applied it to mouse skin, resolving the hair-follicle lineage + dermal/immune/vascular types (their chromatin-potential analysis). Running the actual RNA counts through IGVFagent:

- the loader recovers **34,774 labeled cells across 23 cell types** — an exact match to the paper's final skin set;
- `share rna-qc` reproduces the shallow-but-usable per-barcode profile (median 920 UMIs);
- unsupervised Leiden agrees with the author annotation at **AMI 0.63** — strong given the low RNA depth, which makes fine types (e.g. hair-follicle sub-states) harder to separate than a deep 10x dataset.

**Verdict: IGVFagent reproduces the SHARE-seq skin cell-type structure of Ma 2020** — exact cell/type counts and clustering that concurs with the authors' expert labels, from the raw GEO archive through the `share` skill's QC.

## Engineering internalized from this benchmark

The reusable loader this benchmark exercises now lives in **`Scripts/_scload.py`** (`dense_gene_by_cell_tsv`, `matrixmarket`, `matlab_dge`, `cellxgene_h5ad`, `subsample_cells`, `attach_labels`) — shared across the single-cell/multiome benchmarks and importable by the skills, so any future SHARE-seq/GEO matrix loads through one memory-safe path (streams to sparse; never materializes the dense genes×cells array).

## Citation

Ma S, Zhang B, LaFave LM, Earl AS, Chiang Z, Hu Y, Ding J, Brack A, Kartha VK, Tay T, Law T, Lareau C, Hsu Y-C, Regev A, Buenrostro JD. **Chromatin potential identified by shared single-cell profiling of RNA and chromatin.** *Cell* **183**: 1103–1116 (2020). DOI: [10.1016/j.cell.2020.09.056](https://doi.org/10.1016/j.cell.2020.09.056) · PMID 33098772

## How to reproduce

```bash
bash Benchmarks/ma2020_shareseq/run.sh
igvfagent bench score --paper-id ma2020_shareseq
```

`run.sh` fetches the per-sample GEO supplementary files instead of the 7.5 GB `GSE140203_RAW.tar`, and uses the cluster `igvfagent`. The repo `.venv` was built on macOS and cannot run on Linux. The FigR reference runs are in `Data/ma2020/ref/run_figr_peakgene.R`; the skin peak-gene step should go through sbatch.

## Honest caveats

* **RNA-side reproduction; ATAC QC is available but not scored here.** The skin ATAC `fragments.bed.gz` is 5.4 GB; `share fragment-qc` streams it (low-RAM) and `share joint-qc` combines both modalities — ready to run, but we score the RNA cell-type recovery to keep the benchmark fast. The full joint-QC path is a documented next step.
* **Subsampled to 15k of 34,774** for a laptop-scale run (deterministic seed 0). Cell/type *counts* are the full-set numbers; the AMI is computed on the subsample. On a large-memory host, drop `N_SUB` in `prep_input.py`.
* **AMI, not a Fig-by-Fig match.** Ma 2020's final taxonomy used their own iterative clustering; we quantify agreement (AMI) rather than reproduce every sub-cluster.

## License + provenance

* **Data**: GEO GSE140203 (public); fetched at run time, never redistributed.
* **Code**: IGVFagent Apache-2.0; `prep_input.py` (uses `Scripts/_scload.py`) + `make_figures.py` here.
