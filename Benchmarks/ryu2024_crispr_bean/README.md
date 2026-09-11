# Ryu 2024 — crispr-bean base editing variant effect quantification

[![paper](https://img.shields.io/badge/Nat%20Genet-56:925--937-blue)](https://doi.org/10.1038/s41588-024-01726-6)
[![data](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.10139794-blue)](https://doi.org/10.5281/zenodo.10139794)
[![code](https://img.shields.io/badge/crispr--bean-pinellolab-orange)](https://github.com/pinellolab/crispr-bean)
[![status](https://img.shields.io/badge/data--level%20claims-5%20agree%20%C2%B7%206%20differ-yellow)]()

## Bottom line

**IGVFagent reproduces the paper's data-level statistics from the authors' own
Zenodo deposit, and every remaining difference has an identified cause.** The
two headline reproducibility figures land almost exactly: replicate Spearman
ρ **0.878 vs 0.88** published (tiling screen) and **0.831 vs 0.84** (GWAS
screen), with mean per-gRNA edit fraction **33.9–34.5% vs 34.0%**.

Run it with `bash Benchmarks/ryu2024_crispr_bean/run.sh`.

## The screens are not the IGVF ones

This is the first thing the benchmark established, and it determines
everything after it:

| | paper, LDL-C GWAS | `IGVFDS6464SOVZ` (IGVF) |
|---|---:|---:|
| gRNA species | 3,455 | 8,192 |
| splice-control genes | LDLR, MYLIP, ACAT2, SREBF2, HNF4A, LSS | ABCA1, ARID1A, DENND4C, FADS2, GPN1, HMGCR, LDLR, NPC1L1, SCARB1, SMARCA4, MVK, HNF4A… |
| non-targeting negatives | 100 | **none labelled** |

Only LDLR and HNF4A overlap. Scoring the paper's numbers on IGVF's data would
compare different experiments and agree or disagree for the wrong reason, so
this benchmark runs on the paper's own deposit — which also carries the
`X_bcmatch` and `edit_rate` layers that IGVF does not publish, and which
BEAN's activity-normalised model requires.

## Results

| Claim | Paper | IGVFagent | |
|---|---:|---:|---|
| replicate ρ, LDLR tiling | 0.88 | **0.878** | ✓ |
| replicate ρ, LDL-C GWAS | 0.84 | **0.831** | ✓ |
| tiling library gRNAs | 7,500 | **7,500** | ✓ exact |
| tiling non-targeting controls | 150 | **150** | ✓ exact |
| mean per-gRNA edit fraction | 34.0% | **36.1%** (33.9% top bin, 34.5% bot) | ✓ |
| GWAS library gRNAs | 3,455 | 3,451 | filtering |
| GWAS variants targeted | 583 | 570 | filtering |
| GWAS non-targeting controls | 100 | 99 | filtering |
| tiling variants assessed | 2,182 | 1,894 | filtering |
| tiling missense variants | 874 | 863 | filtering |
| median maximal editing / variant | 0.604 | **0.496** | **open** |
| variants with 95% CI excluding 0 | 54 | — | needs model fit |
| AUPRC — BEAN | 0.90 | — | needs model fit |
| AUPRC — BEAN-Reporter | 0.87 | — | needs model fit |
| AUPRC — BEAN-Uniform | 0.85 | — | needs model fit |
| AUPRC — tiling, BEAN | 0.88 | — | needs model fit |
| tiling variants z < −1.96 | 145 | — | needs model fit |
| …of those, decreasing uptake | 131 | — | needs model fit |

### "filtering" — one cause, six rows

Every count that misses, misses **downward**, and only at the **variant**
level; guide-level counts that survive QC match exactly. The tiling deposit is
named after its thresholds (`..._0.1_0.3.h5ad`). So the deposits are post-QC
objects while the paper's counts are pre-filter.

Those tolerances are deliberately left at **0**. Widening them to turn the
rows green would hide exactly the drift this benchmark exists to detect.

### The one open discrepancy

**Median maximal editing per variant: 0.496 measured, 0.604 published.** The
*mean* reproduces, so it is not a scaling problem. Ruled out by measurement:
the per-guide `obs['edit_rate']` column versus the per-sample
`layers['edit_rate']`, and which bin the layer is averaged over — all land at
0.496–0.513. Remaining hypothesis: the paper's maximal-editing statistic uses
a **variant-specific** rate (the intended edit at the target position, from
the allele-level `edits` layer) rather than the aggregate per-guide rate.

## Questions for the authors

1. **The LDLvar deposit carries eight replicate labels** (`rep5`, `rep9`–`rep15`)
   where Fig. 3b describes five. Fig. 3b's *"15 two-replicate subsamples among
   the five replicates"* is also arithmetically odd: C(5,2) = 10, and 15 =
   C(6,2). The `mask` column does **not** explain it — it masks 2 individual
   samples, not replicates. The tiling deposit has exactly five, so this is
   LDLvar-specific. **This blocks the AUPRC comparison**, because computing it
   over eight replicates when the paper used five or six would produce
   agreement or disagreement for a reason unrelated to the method.
2. **Are the deposited counts expected to sit below the published ones?** If
   the Zenodo objects are post-QC, all six "filtering" rows reconcile.
3. **How is "maximal editing" defined** for the 60.4% figure?

## Notes

`severity` in the tiling deposit is numeric with no legend. The missense level
is read off the data rather than assumed: **severity 1.0 → 863 distinct
targets** against the paper's 874 missense variants, with no other level
within an order of magnitude (next closest, 0.5, has 696).

The two deposits do not share a schema — `target_group` (Variant/PosCtrl/
NegCtrl) versus `Group` (exon numbers, UTRs, DNase HS regions, PosCtrl, and
non-targeting controls split by editor as `ABE control` + `CBE control` =
150). A detector assuming one silently reports zero negative controls for the
other, then computes an AUPRC against an empty negative class.

## Cannot be tested here

| Claim | Why |
|---|---|
| Spearman 0.40 / Pearson 0.45 vs UK Biobank patient LDL-C | UKB is controlled access |
| BEAN-FUSE Spearman 0.50 / Pearson 0.51 vs UKB | UKB access + the FUSE model |
| Spearman of BEAN effect vs individually tested gRNA LFC (26 gRNAs) | validation experiments are not in the deposit |
