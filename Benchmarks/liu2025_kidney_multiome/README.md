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
| Fig. 1A GWAS significant, non-MHC | 99,595 of 13,700,391 | 102,520 of 13,700,391; filter = authors' awk (r=1.0) | reproduced |
| Fig. 1B independent loci | 1,026 | 1,020; merge = authors' R, cM = plink | reproduced |
| Fig. 2B ASE genes (tubule / glom / union) | 8,573 / 8,725 / 10,398 | 8,573 / 8,725 / 10,398 | reproduced (exact) |
| Fig. 3A bASA peaks / SNPs / local+distal | 19,083 / 523,700 / 8,317 | 19,083 / 523,700 / 43.58 % | reproduced (exact) |
| Fig. 3D bASA vs ASE correlation | 0.51 | **0.31** | failing |
| Fig. 4B snASA SNPs, cell-type specific | 25,766, 86 % (3,615 shared) | 25,766, 86.0 % (3,615) | reproduced (exact) |
| Fig. 4C snASA SNPs in open chromatin | 93.2 % | **87.3 %** | failing |
| Fig. 4F PT snASA vs bASA Spearman | 0.90 | 0.895 | reproduced |
| Fig. 5A Open4Gene method | hurdle NB model | port vs R `pscl::hurdle`: zero β and p identical, 98 % of count z identifiable | reproduced |
| Fig. 5C Open4Gene links | 142,452 links, 125,699 peaks, ~90 % single-gene, 97 % positive, median 8 peaks/gene | 142,452, 125,699, 88.7 %, 97.0 %, **median 4** | reproduced (1 check fails) |
| Fig. 5D GWAS variants in Open4Gene peaks | 7,137 → 1,351 genes, 80.2 % single-gene | 7,344 → 1,376, 80.3 % | reproduced |
| Fig. 6B regulatory scorecard | 24,437 variants, 1,060 genes, top rs35716097 → SLC34A1 (27) | all exact | weak (counts only) |
| Fig. 6D coding variants | 1,363 in 782 genes; 37.0 % multi; ALMS1 16 | 1,319; 36.2 %; ALMS1 16 | reproduced |
| Fig. 6E coding + regulatory convergence | 601 genes, 563 validated (93.8 %) | 601, 563 (93.7 %) | reproduced |
| Fig. 6F GCTA-COJO convergence (161 genes) | — | not attempted | pending |
| Re-derive ASE / bASA / snASA / Open4Gene / CARMA from individual-level data | — | — | blocked (controlled access) |

## Honest caveats

* **Three checks fail, and the cause is not yet known.** Re-investigated 2026-09-26; still open, but these alternatives are now ruled out rather than untried:
  * **Fig. 3D (0.31 vs 0.51, need ≥0.49):** the paper does not say how bASA and ASE SNPs are joined. Tried: restricting bASA to local (in-peak) SNPs only (best: 0.458 all-cause dedup, **0.476** local + `Validated == MatrixeQTL&stratAS`, still short); Pearson instead of Spearman (lower in every variant); matching against glomeruli or tubule+glomeruli-combined ASE instead of tubule alone (all lower than tubule); confirmed Ref/Alt alleles already agree between the two tables (no sign-flip bug). Local-only + both-methods-validated bASA vs tubule ASE is the closest found (0.476) but still below tolerance.
  * **Fig. 4C (87.3 % vs 93.2 %, need ≥0.927):** measured against the deposited `Human.Kidney.OpenChromatin.snATAC` BED. Ruled out: coordinate off-by-one (0-based vs 1-based BED shifts the count by 3, not ~1500); a single bad chromosome (the ~87% miss rate is uniform across all 22 autosomes, 82–90%); using the bASA bulk-ATAC peak set instead (far worse, 39%, wrong peak universe — it's the 19,083-peak deposit, not the 524,324-peak snATAC set). Symmetric padding of each peak by 500 bp on both sides lands the fraction at 92.97 %, inside tolerance — but that padding isn't stated anywhere and would be curve-fitting a magic number to the target, so it is not adopted. Likely explanation: the authors' per-analysis peak set differs from the single deposited consensus BED (e.g. wider or per-cell-type calls not deposited).
  * **Fig. 5C median peaks per gene (4 vs 8, need ≥7.5):** the surrounding counts (142,452 links, 125,699 peaks) match exactly on the full deposited significant-link table, confirming that table is the right population — so the median gap isn't a wrong-population bug. Ruled out: AllCell-only subset (still median 4, and its totals don't match the quoted 142,452/125,699 anyway); genes-per-peak median instead of peaks-per-gene (that's 1, unrelated). The distribution is heavily right-skewed (mean 11, median 4, 90th pct 30), so whatever restricted subset yields a median of 8 wasn't found.
* **Fig. 6B is "weak".** Every number matches exactly, but all of them are counts read from the deposited scorecard. None comes from a recomputation.
* **Fig. 6F is not attempted, and this is now a confirmed anti-bot gate, not a dead link.** The supplementary Methods PDF and Tables S1-S31 zip (table S28) exist at fixed PMC `/bin/` URLs (see `OPERATIONS.md` §4) but return a client-side proof-of-work challenge page instead of the file when fetched programmatically — verified 2026-09-26 via direct `curl`/WebFetch. Solving it automatically would be evading that gate, so this port doesn't. A human can open either URL in a browser and drop the files under `Data/liu2025/` to unblock this.
* **Everything upstream of the deposited tables is out of reach.** RASQUAL, snASA calling and CARMA need individual-level genotypes that are controlled access. What is reproduced is the analyses built on the deposited results, plus two ported methods verified against the authors' code.
* **GWAS counts are near, not exact.** Significant variants and loci (102,520 vs 99,595; 1,020 vs 1,026) differ by a few percent, although the port matches the authors' scripts on the same input. The significant-variant count comes before any LD step, so the deposited GWAS file itself probably differs slightly from the one behind the paper's numbers. The loci count also depends on the LD panel: we use public 1000G EUR (GRCh37).

## Provenance

`provenance.json` holds the resolve / harvest / route record. `expected.json` holds every check with its paper quote.
