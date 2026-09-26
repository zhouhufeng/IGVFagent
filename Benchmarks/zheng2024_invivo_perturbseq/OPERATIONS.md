# OPERATIONS — Zheng 2024 in-vivo AAV Perturb-seq, mouse cortex

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

> **Local input required.** The online GEO-metadata step always runs.
> The real reproduction (Fig 4F cell-type proportion test) needs the
> paper's own deposited Seurat objects (~2 GB compressed) and a small R
> package bootstrap described below. `run.sh` exits 77 cleanly until
> both are in place.

## Required input

```
Data/Benchmarks/zheng2024_invivo_perturbseq/raw/GSE249416_Perturb_all.qs
Data/Benchmarks/zheng2024_invivo_perturbseq/raw/GSE249416_Perturb_sg.qs   (optional, secondary cohort)
```

### How to obtain it

GEO's own supplementary files *are* the paper's Seurat objects — no
conversion needed, just decompression. The depositor wrapped `qs`'s own
(already-compressed) serialization format in an extra outer gzip, so a
plain `gunzip -k` recovers a file `qs::qread()` can read directly:

```bash
mkdir -p Data/Benchmarks/zheng2024_invivo_perturbseq/raw
cd Data/Benchmarks/zheng2024_invivo_perturbseq/raw

curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249416/suppl/GSE249416_Perturb_all.qs.gz
curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249416/suppl/GSE249416_Perturb_sg.qs.gz
gunzip -k GSE249416_Perturb_all.qs.gz
gunzip -k GSE249416_Perturb_sg.qs.gz

cd ../../../..
```

### R environment bootstrap (one-time)

`qs` and `SeuratObject` are enough to deserialize and read `@meta.data`
from the paper's Seurat objects — the full `Seurat` package is not
needed for this benchmark (no expression-matrix analysis, only per-cell
metadata: cell type calls, gRNA assignment, 10x channel). On a
Harvard-FASRC-style HPC module system:

```bash
module load R/4.4.1-fasrc01
module load cmake/3.25.2-fasrc01      # needed by RcppParallel's bundled TBB build
export R_LIBS_USER=<a writable R library path>

Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org"));
  install.packages(c("SeuratObject","BH","RApiSerialize"), lib=Sys.getenv("R_LIBS_USER"))'

# qs was superseded by qs2 and archived on CRAN; install its last version from
# the archive. It also requires a stringfish DOWNGRADE: CRAN's current
# stringfish (>=0.17) changed its internal sfstring API (dropped c_str() /
# check_if_native_is_ascii()) in a way that breaks qs 0.27.3's C++ against it.
Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org"));
  install.packages("https://cran.r-project.org/src/contrib/Archive/stringfish/stringfish_0.16.0.tar.gz",
                    repos=NULL, type="source", lib=Sys.getenv("R_LIBS_USER"));
  install.packages("https://cran.r-project.org/src/contrib/Archive/qs/qs_0.27.3.tar.gz",
                    repos=NULL, type="source", lib=Sys.getenv("R_LIBS_USER"))'
```

This pulls in Rcpp, RcppEigen, RcppParallel (bundles + compiles Intel
TBB from source — the slowest step, several minutes), stringfish
(bundles PCRE2), sp, spam, future(.apply) and a handful of small pure-R
packages. Expect a first-time bootstrap to take 15–30 minutes total on
a shared HPC login/compute node; nothing here needs GPU or a build
farm, just patience for the TBB/Eigen template compiles.

## Quick run

```bash
export R_LIBS_USER=<the path you installed into above>
bash Benchmarks/zheng2024_invivo_perturbseq/run.sh
python3 Benchmarks/concordance.py --benchmark zheng2024_invivo_perturbseq
```

## What `run.sh` does

```bash
igvfagent geo series --gse GSE249416          # always runs

# then, if the .qs files are on disk:
Rscript extract_metadata.R  <raw_dir> <derived_dir>       # qs::qread -> @meta.data -> CSV
python3 reproduce_fig4f.py  <derived_dir>/all_metadata.csv <run_dir>
```

`extract_metadata.R` loads `GSE249416_Perturb_all.qs` (the paper's own
published Seurat object, no re-processing) and writes its
`@meta.data` data frame straight to CSV — 50,075 cells × 50 columns,
including `CellType` (author cell-type calls), `assignment` (final
demultiplexed gRNA identity), `orig.ident` (10x channel / "replicate"),
`Keep` (the object's own post-QC flag) and `lowQC`.

`reproduce_fig4f.py` re-derives the paper's Fig 4F headline claim from
that metadata: singlet, QC-passed cells only, cell types with ≥200
cells (matching the paper's own exclusion criterion), then an
arcsin-square-root-transformed proportion test with 10x channel as a
fixed effect, gRNA vs. `NonTarget2` control — the paper's own described
recipe (STAR Methods, "Perturbation-associated analysis": *"proportions
were first transformed using arcsin square root transformation, and the
batch (10x channel) was additionally considered as another fixed effect
... using the propeller.ttest function from speckle"*). `speckle` isn't
installed here, so this is a from-scratch Python re-implementation of
that recipe (statsmodels OLS), not a call into the original R package —
expect the same direction and order of magnitude, not the paper's exact
p-values/FDR.

## Where artefacts land

```
Docs/SingleCell/<ts>_zheng2024_invivo_perturbseq/
├── concordance_metrics.json     # the 4 checks scored below
├── proportion_table.csv         # per (gRNA, channel, cell type) proportions
└── cell_type_proportion_tests.csv  # every gRNA x cell-type contrast vs NonTarget2, BH-FDR corrected
```

## Concordance interpretation

| # | Check | What it verifies |
|---|---|---|
| 1 | `n_cells_total` == 50,075 | The downloaded object is the exact cohort the paper describes ("Across five replicates, we obtained a total of 50,075 cells") |
| 2 | `n_channels` == 5 | Matches the paper's "five replicates" |
| 3 | `foxg1_g1_l6it_reduction_rate` ∈ [2, 20]× (paper: 9.9×) | Direction + order of magnitude of the paper's headline Foxg1-gRNA1 → L6-IT depletion |
| 4 | `foxg1_g1_upper_increase_rate` ∈ [1.05, 4]× (paper: 2.0×) | Direction + order of magnitude of the paper's headline Foxg1-gRNA1 → upper-layer enrichment |

Measured on 2026-09-25: `n_cells_total=50075`, `n_channels=5`,
`foxg1_g1_l6it_reduction_rate=8.22` (p=9.1e-4, FDR=0.081 in this
simplified model vs. paper's FDR=9.6e-6), `foxg1_g1_upper_increase_rate=1.33`
(p=0.18, FDR=0.62 vs. paper's FDR=0.014) — same direction and same
order of magnitude as the paper for both, weaker significance than the
paper's own `propeller.ttest`, as expected for a simplified
re-implementation.

## Honest caveats

* The paper's *most* profound Foxg1 effect — "L6-CT sub-cluster 3"
  (FDR<0.028, 1.8–24.4-fold) — is at a finer clustering resolution than
  the object's own `CellType` column exposes. Reproducing that specific
  number would mean re-deriving the paper's exact sub-clustering, which
  this benchmark does not attempt.
* `speckle::propeller.ttest` is not installed; the proportion test here
  is a plain arcsin-sqrt OLS re-implementation, chosen because it is the
  transform + fixed-effect design the paper's Methods actually describe
  in prose, not because it reproduces `speckle`'s variance estimator.
* `GSE249416_Perturb_sg.qs` (11,688 cells, "AAV serotype secondary
  screen and comparing 5' vs 3' scRNA-seq" per the paper's Methods) is
  downloaded and extracted but not currently used in any scored check.

## License + provenance

* **Paper data**: GEO GSE249416 — public.
* **Code**: IGVFagent Apache-2.0; `qs`/`SeuratObject`/`stringfish` are
  GPL-3/MIT R packages installed at runtime, not vendored.
* **Citation**: Zheng X et al. *Cell* **187**: 3236–3248.e23 (2024).
  doi:10.1016/j.cell.2024.04.050 · PMID:38772369
