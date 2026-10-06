# liu2025_kidney_multiome

[![paper](https://img.shields.io/badge/Science-2025-blue)](https://doi.org/10.1126/science.adp4753)
[![PMID](https://img.shields.io/badge/PMID-39913582-blue)](https://pubmed.ncbi.nlm.nih.gov/39913582/)
[![coverage](https://img.shields.io/badge/reproduction-11%2F16%20analyses-yellow)]()

## Paper

Liu H, Abedini A, Ha E, … Susztak K. **Kidney multiome-based genetic scorecard reveals convergent coding and regulatory variants.** *Science* 2025. doi:[10.1126/science.adp4753](https://doi.org/10.1126/science.adp4753) · PMID 39913582 · PMC12013656

## Data and code

| | Source |
|---|---|
| Derived data | Figshare 26299093 (CC BY): multi-ancestry eGFRcrea GWAS, RASQUAL ASE (tubule, glomeruli), bASA, snASA, Open4Gene links + summary statistics, snATAC peak BED, Kidney Disease Genetic Scorecard. `run.sh` checks the md5 of each file. |
| Authors' code | `hbliu/Kidney_Epi_Pri@552c461`, `hbliu/Open4Gene@6e7f36a` |
| Primary data | CMDGA individual-level kidney data and matched WGS: controlled access (AMP consortium sign-in), not used |

## What IGVFagent does

Two of the paper's methods are now built into IGVFagent as ports, and each is checked against the authors' own code:

* `liu2025-gwas-loci` is a port of `eGFR_GWAS/Independent.Loci.eGFR.GWAS.sh`: P < 5e-8, MHC removal, plink clumping against 1000G EUR, 0.1 cM merge.
* `liu2025-open4gene` is a port of `R/Open4Gene.R`, the hurdle negative-binomial peak→gene test.

`verify_derived_tables.py` recomputes the Fig. 2–6 numbers from the deposited tables.

```bash
sbatch -p <partition> -c 4 --mem 32G -t 2:00:00 --wrap "bash Benchmarks/liu2025_kidney_multiome/run.sh"
igvfagent bench score --paper-id liu2025_kidney_multiome
```

## Paper coverage

`igvfagent bench score` → **reproduction: incomplete, 11/16 analyses** (43/46 checks pass; run `Docs/Benchmark/20260926_003843_liu2025_kidney_multiome`).

| Analysis | Paper | IGVFagent | State |
|---|---|---|---|
| Fig. 1A GWAS significant, non-MHC | 115,726 at P<5e-8 → 112,490 after heterogeneity → 99,595 non-MHC | 115,725 at P<5e-8; heterogeneity removal implied by the author's counts 3,235 (author 3,236); 102,520 non-MHC without the per-study rescue test | numbers reproduced; panel not redrawn |
| Fig. 1B independent loci | 1,026 | 1,020; merge = authors' R, cM = plink | reproduced |
| Fig. 2B ASE genes (tubule / glom / union) | 8,573 / 8,725 / 10,398 | 8,573 / 8,725 / 10,398 | reproduced (exact) |
| Fig. 3A bASA peaks / SNPs / local+distal | 19,083 / 523,700 / 8,317 | 19,083 / 523,700 / 43.58 % | reproduced (exact) |
| Fig. 3D bASA vs ASE correlation | 0.51 (1,482 peaks, 1,235 genes) | **0.520** with the authors' selection (1,410 peaks, 1,186 genes); panel redrawn and reviewed | reproduced |
| Fig. 4B snASA SNPs, cell-type specific | 25,766, 86 % (3,615 shared) | 25,766, 86.0 % (3,615) | reproduced (exact) |
| Fig. 4C snASA SNPs in open chromatin | 93.2 % | 93.20 % (24,015/25,766, table S18) | reproduced (exact) |
| Fig. 4F PT snASA vs bASA Spearman | 0.90 | 0.895 | reproduced |
| Fig. 5A Open4Gene method | hurdle NB model | port vs R `pscl::hurdle`: zero β and p identical, 98 % of count z identifiable | reproduced |
| Fig. 5C Open4Gene links | 142,452 links, 125,699 peaks, ~90 % single-gene, 97 % positive, median 8 peaks/gene | 142,452, 125,699, 88.7 %, 97.0 %, **median 4** | 1 check fails; the first author is re-checking the paper's median |
| Fig. 5D GWAS variants in Open4Gene peaks | 7,137 → 1,351 genes, 80.2 % single-gene | 7,344 → 1,376, 80.3 % | reproduced |
| Fig. 6B regulatory scorecard | 24,437 variants, 1,060 genes, top rs35716097 → SLC34A1 (27) | all exact | weak (counts only) |
| Fig. 6D coding variants | 1,363 in 782 genes; 37.0 % multi; ALMS1 16 | 1,363; 782; 37.0 %; ALMS1 16 (table S24) | reproduced (exact) |
| Fig. 6E coding + regulatory convergence | 601 genes, 563 validated (93.8 %) | 601, 563 (93.7 %) | reproduced |
| Fig. 6F GCTA-COJO convergence (161 genes) | — | not attempted | pending |
| Re-derive ASE / bASA / snASA / Open4Gene / CARMA from individual-level data | — | — | blocked (controlled access) |

## What we need from the authors

This reproduction uses only what is public: the Figshare deposit 26299093, Supplementary Tables S1–S31, and the code in `hbliu/Kidney_Epi_Pri` and `hbliu/Open4Gene`. The items below are what still blocks an exact match, ordered by how much each one unlocks. For each item, any one of the listed forms is enough. Derived, non-identifying tables are all we ask for. Where something is controlled access, the last section describes a way to run those steps on the authors' side instead.

### A. Exact numbers that are off by a few percent

| # | Analysis | What is missing | Why it matters | Simplest form to send |
|---|---|---|---|---|
| A1 | **Fig. 1A**: 99,595 significant non-MHC variants. We get 102,520. | The rescue test for heterogeneous variants: 30,865 significant variants have HetISq > 50 or HetPVal < 0.05, and 3,236 of them were removed because none of the 17 original studies reached P < 5×10⁻⁸ on its own. The per-study results are not deposited. | Everything else in the chain already matches: 115,725 significant variants, 13,205 MHC variants, and an implied 3,235 removed. This one test is the whole remaining gap of 2,925 variants. | Either the 99,595-variant list (rsID or chr:pos), or for the 30,865 heterogeneous variants a yes/no column for "P < 5×10⁻⁸ in at least one study". Per-study summary statistics are not needed. |
| A2 | **Fig. 1B**: 1,026 independent loci. We get 1,020. | The variant list from A1, plus confirmation of the clumping reference. We use 1000G phase 3 EUR, as named in `Independent.Loci.eGFR.GWAS.sh`, with plink `--clump-p1 5e-8 --clump-r2 0.1 --clump-kb 10000` and a 0.1 cM merge. | The port already matches the authors' clumping and merge code on the same input. The loci gap follows from the A1 gap and possibly the exact 1000G build and duplicate-ID handling. Table S4 lists the 1,026 loci, so the remainder can be checked locus by locus. | The A1 list, plus the 1000G file actually used (`reduced.ALL.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.EUR.rmdup`: how it was reduced and deduplicated). |
| A3 | **Fig. 3D**: bASA vs ASE. Now 0.520 vs 0.511 with the code the first author sent. | How `bASA.Tub.ASE.TSS` was built: the join between bASA and tubule ASE, which TSS file, and how `DisPeak2TSS` was computed. | Our selection gives 52,375 rows within 10 kb, then 3,592 rows, 1,410 peaks and 1,186 genes. The authors' code prints 55,628, 3,692, 1,482 and 1,235. The correlation already matches; the counts are 3–5% short. | The lines that create `bASA.Tub.ASE.TSS`, or the file `CodingGene_gencode_v35lift37_hg19_ID_sorted.TSS.autosome.bed` plus the bedtools call. |
| A4 | **Fig. 5C**: median peaks per gene, 8 in the paper. We get 4. | The first author is re-checking this number. | Every other Fig. 5C count matches exactly on the deposited significant-link table: 142,452 links, 125,699 peaks, 88.7% single-gene and 97.0% positive. | Confirmation of the number, and which table and filter it was computed on. |
| A5 | **Fig. 6F**: GCTA-COJO convergence, 161 genes. | The LD reference used for COJO. We use 1000G phase 3 EUR (n = 503); the meta-analysis is multi-ancestry. | It decides whether our conditional estimates can match table S28's effect, SE and P exactly, or only approximately. | The name or version of the LD reference panel and the GCTA flags. |
| A6 | **Fig. 6B**: Kidney Disease Genetic Scorecard (24,437 regulatory variants; 1,060 genes). | The code that computes the prioritization score from its evidence columns. | Every count matches, but only because we read the deposited scorecard. Without the scoring code the score can't be recomputed, so this analysis stays "weak". | The scoring script, or a short description of the points each evidence type adds. |

### B. Code for this paper's figures

`hbliu/Kidney_Epi_Pri/Figure_Code/` holds the figure code of the group's earlier kidney paper (meQTL and eQTL Manhattan plots, h2med enrichment, SLC47A1 PheWAS), not this paper's. Apart from the Fig. 3D snippet sent on 2026-10-06, no plotting code for Figs. 1–6 is public. To count as reproduced, every data panel has to be redrawn from our own output and compared side by side with the paper's (see `Benchmarks/figure_rule.py`). We will redraw from the descriptions if needed, but the R scripts behind Figs. 1–6, or the Science "source data" tables for each panel, would let us match axes, filters and groupings exactly rather than by eye.

### C. Data behind panels that the public deposit can't draw

| Panels | What they need | Is it deposited? |
|---|---|---|
| Fig. 2A, 2D, 2G (ASE Manhattan and locus-zooms) | Genome-wide RASQUAL ASE results, all tested SNP–gene pairs with P values, for tubule and glomeruli | No. Figshare has only the significant (FDR 1%) associations. |
| Fig. 2C, 2D, 2G (ASE vs eQTL) | Kidney eQTL results (MatrixeQTL, same samples) | Not in this paper's deposit. Please point us to it if it is the group's published kidney eQTL resource. |
| Fig. 3A (bASA Manhattan, 650,273 bASAs) | The bASA rows plotted | Almost. Figshare has 650,237 rows, 36 fewer than the 650,273 in the panel title; please say which 36 were dropped, or whether the file is a later version. |
| Fig. 3B (enrichment in chromatin states) | The kidney ChromHMM states used | Not deposited. Please name the source. |
| Fig. 3C, 3E, 3F (bASA vs meQTL, pies, USP24 locus) | Kidney meQTL results and the ATAC/H3K27ac/H3K4me3 tracks | Not in this paper's deposit. Please name the source. |
| Fig. 1C, 6F (validation by eGFRcys/BUN) | eGFRcys and BUN GWAS used for validation | Please name the summary statistics used. |
| Fig. 1F (credible-set enrichment) | The non-synonymous, enhancer and promoter annotation files | Please name the annotation sources. |
| Fig. 4I, 5G (cell-type tracks at example loci) | snATAC pseudobulk coverage per cell type | Controlled access (CMDGA); see D. |

### D. Analyses that need individual-level data (controlled access)

Re-deriving ASE (RASQUAL, 613/522 samples), bASA (82 samples), snASA, Open4Gene from cells, and CARMA meta fine-mapping needs genotypes and reads behind CMDGA/AMP access and per-cohort biobank summary statistics. We don't ask for those files. If the authors want these steps independently checked, two options:

1. **Run our ports on their side.** `igvfagent liu2025-open4gene` (the hurdle model, verified against `R/Open4Gene.R`) and `igvfagent liu2025-gwas-loci` run on any cluster with Python. We'd send exact commands; only the summary outputs (counts, effect sizes) need to come back.
2. **A data use agreement** through CMDGA/AMP for the reproduction, if the authors prefer that route.

Either way, the result would be recorded here as "reproduced by the authors with IGVFagent's port", kept separate from what we reproduced ourselves.

## Author review (2026-10-06)

The first author reviewed this reproduction and sent three comments, now applied:

* **Fig. 1A.** Besides P < 5e-8 and the MHC filter, variants with between-study heterogeneity (HetISq > 50 or HetPVal < 0.05) count only if they also reach P < 5e-8 in at least one of the 17 original studies. Their counts: 115,726 → 112,490 (−3,236 heterogeneous) → 99,595 (−12,895 MHC). The deposited meta-analysis carries `HetISq`/`HetPVal`, and every step of that chain reproduces to within one variant:
  * 115,725 variants at P < 5e-8; the author's 115,726 is a line count that includes the header.
  * 13,205 MHC variants before the heterogeneity filter, so 310 of the removed variants were MHC.
  * Implied removal: (102,520 − 99,595) + (13,205 − 12,895) = 3,235 vs 3,236.
  * 30,865 significant variants carry the heterogeneity flag. The author kept all but 3,236 of them because they were also significant in an individual study. That rescue test needs the per-study files, which are not deposited, so 99,595 itself can't be recomputed here. The gap from our 102,520 is entirely this step, not a different input file.
* **Fig. 3D.** The correlation uses only ATAC peaks within 10 kb of the gene's TSS, to keep one peak from targeting several genes, and per peak keeps the top bASA SNP and then the closest TSS. Code: `input/Fig. 3D bASA vs ASE correlation.sh`. Applied in `verify_derived_tables.py` with GENCODE v35lift37 protein-coding TSSs, the authors' annotation, distance measured from the peak edge as `bedtools closest` does:
  * Spearman **0.520** vs 0.511. The earlier unrestricted SNP join gave 0.31.
  * 1,410 peaks and 1,186 genes vs 1,482 and 1,235, and 52,375 rows within 10 kb vs 55,628. The 3–5% shortfall is probably the authors' intermediate table, which wasn't deposited and couldn't be rebuilt exactly.
  * The panel is redrawn by `fig3d_plot.py`, saved to `figures/fig3d_basa_vs_ase.png`, and reviewed against the paper's: same trend, slope and band near ASE = 0.
* **Fig. 5C median peaks per gene (4 vs 8).** The author isn't sure of the paper's 8 and will double-check it. The check stays failing until they confirm.

## Honest caveats

* **The Supplementary Tables workbook (S1-S31) was fetched by hand 2026-09-26** (PMC gates it behind a browser-only proof-of-work challenge; see `OPERATIONS.md` §4) and dropped at `Data/liu2025/supplement/science.adp4753_tables_s1_to_s31.xlsx`. `verify_derived_tables.py` uses it automatically when present. It resolved two of what were three failing checks:
  * **Fig. 4C** was 87.3% against the deposited peak BED; table S18 gives the authors' own per-SNP `Dis2Peak`, and 24,015/25,766 = **93.20%**, an exact match.
  * **Fig. 6D** was an approximation (1,319 variants) read off the 957-locus `Scorecard.xlsx`; table S24 is the actual 1,363-variant coding-variant table behind the paper's number, and it reproduces **exactly**: 1,363 variants, 782 genes, 37.0% multi-variant, ALMS1 top with 16.
  * Table S25 (601 candidate coding+regulatory pairs with a ground-truth `Convergence` flag) and table S28 (the 161 convergent genes with the authors' own GCTA-COJO conditional effect/SE/P and LD r²) further make **Fig. 6F reproducible without any controlled-access data** — see below.
* **Superseded by the author's review above for Fig. 3D.** The investigation below predates it.
* **One check still fails, cause not yet known.** Re-investigated 2026-09-26; these alternatives are ruled out rather than untried:
  * **Fig. 3D (0.31 vs 0.51, need ≥0.49):** the paper does not say how bASA and ASE SNPs are joined, and none of the fetched supplementary tables pair bASA and ASE effect sizes at SNP level (table S16's GWAS-colocalization subset gives 0.73, but that's a different, GWAS-selected population, not the general claim). Tried: restricting bASA to local (in-peak) SNPs only (best: 0.458 all-cause dedup, **0.476** local + `Validated == MatrixeQTL&stratAS`, still short); Pearson instead of Spearman (lower in every variant); matching against glomeruli or tubule+glomeruli-combined ASE instead of tubule alone (all lower than tubule); confirmed Ref/Alt alleles already agree between the two tables (no sign-flip bug). Local-only + both-methods-validated bASA vs tubule ASE is the closest found (0.476) but still below tolerance.
  * **Fig. 5C median peaks per gene (4 vs 8, need ≥7.5):** the surrounding counts (142,452 links, 125,699 peaks) match exactly on the full deposited significant-link table, and the Methods PDF confirms the population definition ("FDR < 0.01 in either the zero model or the count model" — every one of the 171,877 deposited rows already satisfies this) — so the median gap isn't a wrong-population bug. Ruled out: AllCell-only subset (still median 4, and its totals don't match the quoted 142,452/125,699 anyway); genes-per-peak median instead of peaks-per-gene (that's 1, unrelated). The distribution is heavily right-skewed (mean 11, median 4, 90th pct 30), so whatever restricted subset yields a median of 8 wasn't found. Table S20 (top peak per gene) doesn't carry a per-gene peak count either.
* **Fig. 6B is "weak".** Every number matches exactly, but all of them are counts read from the deposited scorecard. None comes from a recomputation.
* **Everything upstream of the deposited tables is out of reach.** RASQUAL, snASA calling and CARMA need individual-level genotypes that are controlled access. What is reproduced is the analyses built on the deposited results, plus two ported methods verified against the authors' code.
* **GWAS counts are near, not exact.** Significant variants and loci (102,520 vs 99,595; 1,020 vs 1,026) differ by a few percent, although the port matches the authors' scripts on the same input. The significant-variant count comes before any LD step, so the deposited GWAS file itself probably differs slightly from the one behind the paper's numbers. The loci count also depends on the LD panel: we use public 1000G EUR (GRCh37).

## Provenance

`provenance.json` holds the resolve / harvest / route record. `expected.json` holds every check with its paper quote.
