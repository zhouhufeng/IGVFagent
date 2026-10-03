# Tabula Sapiens 2.0 skill

Retrieval and reproduction for the Tabula Sapiens 2.0 human cell atlas
(Tabula Sapiens Consortium, *Cell* 2026): **1,136,218 cells, 28 tissues,
24 donors, 182 fine / 38 broad cell types, 701 tissue-cell-type
populations.**

## Which data routes are actually open

| Route | Content | Access |
|---|---|---|
| figshare 27921984 | 28 processed `.h5ad`, **57 GB** | open |
| GEO GSE306755 | count matrices + full 1.1M-cell metadata | open |
| CELLxGENE | 35 per-tissue datasets | open |
| AWS `czb-tabula-sapiens` | raw FASTQ/BAM, **103.16 TB** | **listable, NOT readable** |

The S3 bucket is gated, and this is worth stating plainly because the
open-data registry listing implies otherwise: `ListBucket` succeeds, but
**every `GetObject` returns 403 AccessDenied**, on v1 and v2 alike, with
or without `x-amz-request-payer`. That is the data transfer agreement the
paper describes. `s3-manifest` probes this on each run and reports what
it found rather than assuming.

Measured composition of the 103.16 TB: 43.2 TB BAM, 32.1 TB FASTQ,
27.6 TB STAR per-cell intermediates, and only **0.32 TB** of count
matrices and metrics. Even if it were readable, the interesting 0.3% is
already public on figshare and GEO in a better form.

## Getting the data

```bash
igvfagent tabula pull-figshare --list-only          # see the 28 files
igvfagent tabula pull-figshare --tissues Lung,Heart # or just what you need
igvfagent tabula pull-geo                           # cell metadata (41 MB)
igvfagent tabula pull-cellxgene --list-only
igvfagent tabula s3-manifest                        # inventory, no download
igvfagent tabula status                             # local state vs the paper
```

`status` checks what you have against the paper's own numbers and prints
`OK` / `!!` per quantity -- cells, donors, tissues, cell types,
populations, droplet/FACS split.

## Reproducing the figures

```bash
# Figure 1 -- donors, tissues, composition (needs only the 41 MB metadata)
igvfagent tabula overview

# Figures 2-3 -- TF specificity and what the ubiquitous TFs do.
igvfagent humantfs build-db                          # 1,639 curated TFs
igvfagent tabula tf-matrix --label tf_means          # gene x cell-type means
igvfagent tabula tf-specificity --matrix <run>/mean_expression.npz
igvfagent tabula tf-enrichment --tau-table <run>/tf_tau.tsv

# Regulon ACTIVITY, as opposed to mere expression (SCENIC, clean-room)
igvfagent tabula tf-regulons --tissues Lung,Heart

# Figure 4 -- senescent-cell burden
igvfagent tabula senescence --hallmarks

# Figure 5 -- sex differences, pseudobulked by donor
igvfagent tabula sex-de

# Figure 6 -- donor clinical metadata (the ChatTS core, offline)
igvfagent tabula donors --min-age 60
```

## What it reproduces

Scored by `Benchmarks/quake2026_tabula_sapiens` at **26/26**:

| Quantity | Ours | Paper |
|---|---:|---:|
| Cells (droplet / FACS) | 1,136,218 (1,093,048 / 43,170) | identical |
| Donors / tissues / fine types / populations | 24 / 28 / 182 / 701 | identical |
| Droplet fine / broad types | 175 / 38 | identical |
| Donor ages, sex split, age groups | 22-74, 11M/13F, 7/11/6 | identical |
| TFs with zero expression anywhere | **2: SHOX, ZBED1** | 2: SHOX, ZBED1 |
| TF specific / non-specific (tau > 0.85, 175 cell types) | 882 / 741 | 890 / 745 |
| GO terms for non-specific TFs (padj < 0.02) | 70 | 69 |

Biology checks that are not merely counting: FOXP3 peaks in regulatory
T cells, germ-cell TFs in spermatogenic cells, FOXI1 in ionocytes, and
all eight of the paper's named ubiquitous TFs fall below tau 0.85
(max 0.616).

## Method notes

**tau depends on the number of cell types.** `tau = sum(1 - x/max(x)) /
(N - 1)` over N cell types, computed on the mean **log-normalised**
expression, droplet subset only -- exactly the paper's recipe. But N is
in the denominator, so a run over 46 cell types from 3 tissues gives
systematically different values than the paper's 175. Both `tf-matrix`
and `tf-specificity` print a warning whenever fewer than 28 tissues are
present. Do not compare a partial run's specific/non-specific counts to
the paper's 890 / 745.

**Genes expressed nowhere get NaN, not zero.** A gene with no counts in
any cell type carries no specificity information; forcing it to 0 or 1
would put it at one end of the distribution. The paper's SHOX and ZBED1
are exactly this case, and `tf-specificity` reports them separately.

**Broad cell classes: 40 or 38?** The full metadata has 40; the paper
says 38. Restricting to the droplet subset -- as the paper's analysis
does -- gives 38, and fine types 182 -> 175. Both reconcile exactly.

**Streaming, not loading.** The atlas is 1.1M x 62k. `mean_by_group`
accumulates per-group sums across tissue files, subsetting genes on the
sparse matrix *before* densifying. Restricting to the 1,639 TFs makes
the run roughly 38x cheaper than pulling all genes.

## Provenance

Apache-2.0. Algorithms are reimplemented from published descriptions; no
upstream source is copied or vendored.

| Capability | Reference | License | Approach |
|---|---|---|---|
| Analysis notebooks | [czbiohub-sf/tabula-sapiens](https://github.com/czbiohub-sf/tabula-sapiens) | BSD-3-Clause | clean-room from `paper2/` |
| tau statistic | tspex 0.6.3 | MIT | clean-room, formula from the Methods |
| Ambient RNA | DecontX (celda 1.16.1) | MIT | clean-room variational EM |
| Consensus modules | cNMF 1.5.4 | MIT | clean-room |
| Regulon activity | pySCENIC 0.12.1 | **GPL-3** | clean-room -- not imported |
| Enrichment | GSEApy | MIT | existing `igvfagent enrich` |
| TF list | [Human TFs](https://humantfs.ccbr.utoronto.ca) (Lambert 2018) | see source | `igvfagent humantfs` |
| Alignment | STAR 2.7.11b / CellRanger 7.0.1 | GPL-3 / 10x EULA | **external tools** -- orchestrated, not reimplemented |

Runtime dependencies: numpy, scipy, pandas, matplotlib, scikit-learn,
anndata, scanpy. No GPL runtime dependencies.
