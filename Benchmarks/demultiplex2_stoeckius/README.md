# deMULTIplex2 (Zhu, Conrad & Gartner, Genome Biol 2024) — reproduction

[![paper](https://img.shields.io/badge/Genome%20Biol-25%3A37-blue)](https://doi.org/10.1186/s13059-024-03177-y)
[![PMID](https://img.shields.io/badge/PMID-38291503-blue)](https://pubmed.ncbi.nlm.nih.gov/38291503/)
[![PMC](https://img.shields.io/badge/PMC-PMC10829271-blue)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10829271/)
[![coverage](https://img.shields.io/badge/reproduction-7%2F13%20analyses-yellow)]()

## Bottom line

`bench score`: **reproduction incomplete, 7/13 analyses** (24/24 checks pass).
Reproduced, each with a class-A check against the paper and a class-B check
against the authors' R code on the same input:

* Table S3 deMULTIplex2 F-scores for **lung cell line, BAL1, BAL2, BAL3, Gaublomme**;
* Fig. 2B deMULTIplex2 F-scores on **simulations 1-5**;
* the classifier itself (`multiseq`) against R `demultiplexTags`.

The other six analyses need ground-truth labels that were never deposited
(3) or that would have to be recomputed from raw reads by genotype
demultiplexing (3); see [Not reproduced](#not-reproduced).

## Results

deMULTIplex2 average F-score over samples (the authors' `confusion_stats`).
"R here" is the authors' own code (deMULTIplex2 @ `de48333`, `R/em.R`,
`R/classify.R`, `R/benchmarking.R` sourced unmodified) run on the same inputs
with the paper's per-dataset settings; "IGVFagent" is `igvfagent
demultiplex2-benchmark`. "Calls identical" is the per-cell match between the two
with subsampling switched off (`max.cell.fit` = 1e9), which removes the only
source of randomness that R and numpy cannot share.

| Analysis | Truth | Paper | R here | IGVFagent | Calls identical |
|---|---|---:|---:|---:|---:|
| Table S3 lung cell line (3 captures) | vireo donors, Oshlack repo | 0.880 | 0.882 | 0.880 | 45,957 / 45,957 |
| Table S3 BAL1 (2 captures) | vireo donors | 0.907 | 0.915 | 0.915 | 24,091 / 24,091 |
| Table S3 BAL2 | vireo donors | 0.829 | 0.848 | 0.848 | 48,841 / 48,841 |
| Table S3 BAL3 | vireo donors | 0.774 | 0.783 | 0.785 | 62,306 / 62,306 |
| Table S3 Gaublomme (nuclei hashing) | demuxlet, demuxEM docker image | 0.976 | 0.976 | 0.976 | 2,860 / 2,860 |
| Fig. 2B simulation 1 | simulated | 1.00 | 0.998 | 0.998 | |
| Fig. 2B simulation 2 | simulated | 1.00 | 0.999 | 0.999 | |
| Fig. 2B simulation 3 | simulated | 1.00 | 0.997 | 0.997 | 125,204 / 125,204 |
| Fig. 2B simulation 4 | simulated | 0.99 | 0.994 | 0.994 | (all five) |
| Fig. 2B simulation 5 | simulated | 0.87 | 0.868 | 0.867 | |

Classifier check: on the Stoeckius PBMC matrix bundled with deMULTIplex2
(15,113 cells x 8 HTOs), all 120,904 per-cell, per-tag posteriors from
`multiseq demultiplex` equal R `demultiplexTags` (max |diff| 7.6e-14).

Where the reproduction differs from the paper: on BAL1-3 recall matches Table
S3 to three decimals, but precision comes out 0.017-0.047 higher, and so F is up
to 0.019 higher. This happens with the authors' R code too, so it is a property
of the inputs rather than of the port. The cause is not established. The
lung/BAL scripts were not published with the benchmark repo, so the capture
pooling ("bc_cbn") and the deposited donor-label version may differ from what
the authors used. Class-A tolerance is ±0.025 F for Table S3 and ±0.015 for
Fig. 2B (the figure prints two decimals).

## Ported into IGVFagent

* **`multiseq demultiplex`**: the classifier core was rewritten as a
  line-by-line port of `demultiplexTags` / `fit.em` / `m.step` / `e.step` and
  `MASS::glm.nb` (IRLS + `theta.ml`). The earlier clean-room version differed
  from R in five ways: it used hard rather than soft mixing weights after
  initialisation, maxed rather than replaced the positive-class likelihood
  below its mean, checked `min.cell.fit` after trimming, raised an error on a
  failed refit where R keeps the previous fit, and did not expose
  `max.cell.fit` or the quantile window. New flags: `--max-cell-fit`,
  `--min-cell-fit`, `--converge-threshold`, `--min-quantile-fit`,
  `--max-quantile-fit`.
* **`demultiplex2-benchmark`** (registered port,
  `Scripts/ported/skills/demultiplex2_benchmark.py`): port of
  `R/benchmarking.R` `benchmark_demultiplex2` + `confusion_stats`. Given
  R's own calls it reproduces R's precision, recall and F to machine precision
  on all four Oshlack datasets. `verify-port` match rate 1.000 on every dataset
  above.

## Not reproduced

| Analysis | Paper | Why |
|---|---:|---|
| Table S3 Stoeckius cell line | 0.962 | Truth `hto12$assign_rna` (incl. `Doublet_rna`) came from an undocumented processing ("following the Seurat tutorial", which has no such labels). It is not in GitHub, Zenodo 8429628 (code only) or GEO GSE108313. **not_deposited** |
| Table S3 Winkler PDX | 0.557 | The tag→tumor table `PDX_MULTI-SEQ_Metadata.csv` is not deposited. The GEO GSE211145 metadata pairs only 6 of the 15 batch-3 tags. **not_deposited** |
| Fig. 4 PDX: 63.2% of cells correctly retrieved | 0.632 | Same missing tag→tumor table. **not_deposited** |
| Table S3 Stoeckius PBMC | 0.957 | Truth = the authors' vireo run on the GSE108313 FASTQs; not deposited. Recomputing it needs cellranger + vireo on raw reads, which was not done here. **other** |
| Table S3 McGinnis MULTI-seq / SCMK | 0.970 / 0.823 | Truth = the authors' souporcell run on the GSE161329 RNA reads; not deposited. Same raw-read recomputation needed. **other** |

The three "other" analyses use public raw data, so they keep the verdict at
`incomplete` rather than `reproduced_except_access`. Genotype clustering
(vireo, souporcell) is unsupervised, so rerunning it would give new truth
labels, not the authors' own.

## How to reproduce

```bash
bash Benchmarks/demultiplex2_stoeckius/run.sh           # ~10 min; REF=1 also runs the authors' R code
igvfagent bench score --paper-id demultiplex2_stoeckius
```

`run.sh` fetches the Stoeckius PBMC matrix (deMULTIplex2 `data/`), Oshlack/hashtag-demux-paper
@ `3be94bb` (lung + BAL counts and vireo donors) and, if Singularity is available,
the regevlab/demuxem image (Gaublomme ADT, demuxlet `.best`, RNA matrix). It
regenerates the five simulations with the authors' `simulateTags` in R. Reference
scripts are in `ref/`: `run_demux2_ref.R`, `run_bench_ref.R`,
`make_simulations.R`, `prep_oshlack.py`, `prep_gaublomme.py`, `score_calls.py`.

The original behaviour checks on the bundled PBMC matrix still run (no truth in that
bundle): 15,113 cells, 8/8 donor groups, singlet rate 0.833, pool balance 0.79.

![Classification](figures/fig1_classification.png)

## Citation

* Zhu Q, Conrad DN, Gartner ZJ. **deMULTIplex2: robust sample demultiplexing for scRNA-seq.** *Genome Biology* **25**:37 (2024). DOI [10.1186/s13059-024-03177-y](https://doi.org/10.1186/s13059-024-03177-y) · PMID 38291503 · PMC10829271
* Data: Stoeckius et al. *Genome Biol* 2018 (GSE108313); Howitt et al. *NAR Genom Bioinform* 2023 and Maksimovic et al. 2022 (lung, BAL); Gaublomme et al. *Nat Commun* 2019 (demuxEM image).

## Corrections to the earlier version of this page

* PMID was given as 38217022 (an unrelated sports-medicine paper); the paper is PMID 38291503, article 37, not 20.
* The earlier "reproduction" scored only behaviour ranges (singlet rate "~75-90 %", "consistent with 10x super-loading"). The paper does not state those ranges, and none of those checks compares against the paper. Those claims are removed.
* deMULTIplex2 is licensed CC BY 4.0 (package DESCRIPTION), not MIT.
