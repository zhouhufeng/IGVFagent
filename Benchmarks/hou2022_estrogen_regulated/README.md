# hou2022_estrogen_regulated

## Paper

**Analysis of estrogen-regulated enhancer RNAs identifies a functional motif required for enhancer assembly and gene expression.**
Hou TY, Kraus WL.
*Cell Reports* 2022 · doi:[10.1016/j.celrep.2022.110944](https://doi.org/10.1016/j.celrep.2022.110944) · PMID 35705040 · PMC9246336 · IGVF award IGVF0056, grant UM1HG011996

Resolver confidence: **1.00** (exact PMID match).

The paper annotates estrogen (E2)-regulated enhancer RNAs (eRNAs) genome-wide in MCF-7 breast cancer cells (PRO-cap + RNA-seq), discovers a ~40nt sequence motif enriched in them (FERM, via MEME/FIMO), shows by CRISPR/dCas9 targeting that some eRNAs (containing FERM) stimulate their cognate gene's E2 response by recruiting ERα/SRCs and stimulating p300 histone-acetyltransferase activity, identifies BCAS2 as a FERM-interacting protein by mass-spec pulldown, and connects one FERM-containing eRNA (from the PRRX2 locus) to breast-cancer proliferation and patient survival.

**The paper's own words: "This paper does not report original code."** Every analysis below is therefore implemented directly from the STAR★METHODS text (no author repository to clone/run as reference) — the same situation as `nishizaki2020_predicting_effects`'s SEMpl downstream comparisons, which is the structural template this benchmark follows.

## Why this doesn't use an IGVFagent assay route

Automated routing scores this paper highest for `chipatlas`/`spatial_atac_hic`/`sc_analyze` purely from incidental keyword mentions (ChIP-seq, Hi-C, scRNA-seq appear in citations/methods, not as this paper's own primary assay); none actually fits (`igvfagent bench route` finds no confident match). The paper's own assays are PRO-cap, stranded RNA-seq (polyA-depleted/enriched), ChIP-seq, and LC-MS/MS proteomics — reproduced here directly against GEO/METABRIC data rather than through a route template.

## Part 1 — Fig. 1 & Fig. 5: eRNA annotation, PWM, and the FERM motif

*(PRO-cap + RNA-seq -> eRNA calls -> position weight matrix -> MEME/FIMO motif discovery.)*

**Status 2026-10-02: alignment in progress, no Part 1 result yet.** The Sep-27 pipeline run aligned 12 of 14 samples, but a size check against ENA's `fastq_bytes` showed 5 of its raw FASTQs were truncated (SRR14634910, SRR14634913, SRR14634924_2, SRR14634925_1/_2), so the BAMs built from them (E2 PRO-cap rep1 ±TAP and their merges) were set aside in `Data/hou2022_erna/bam/stale_truncated_input/`. Size-verified downloads resume in job `49906815` (`hou2022_dl2`), and the full pipeline (`49906827`, `afterok`) then realigns only what's missing and runs StringTie → MACS2 → TSS → eRNA annotation → FERM motif → Part 3. See OPERATIONS.md §3.

## Part 2 — Fig. 7: Kaplan-Meier survival (PRRX2, ER+ breast cancer)

Claim: *"High expression of PRRX2 leads to poorer relapse-free survival for patients with ER+ breast cancer"* (STAR Methods: "Kaplan-Meier analysis was performed using Kaplan-Meier Plotter (Nagy et al., 2021)... Patients were split by the median expression... options were set to default").

Kaplan-Meier Plotter (kmplot.com) gates its query form behind a CAPTCHA (`captcha.php` in its submission form) — not attempted to bypass. Instead, the `hou2022_prrx2_survival` port (registered into IGVFagent as `igvfagent hou2022-prrx2-survival`, pinned to the cBioPortal software repo since Kaplan-Meier Plotter itself has no scriptable API) reimplements the described analysis (median split, ER-stratified, log-rank test) against **METABRIC** (Curtis et al. 2012 *Nature*; Pereira et al. 2016 *Nat Commun*), a real, public, array-based (Illumina HT-12 v3) breast cancer cohort with matching expression + ER status + relapse-free-survival annotation, served by the cBioPortal REST API (`brca_metabric`, 2509 patients, no CAPTCHA).

**Result: not reproduced.** In 1,506 ER+ METABRIC patients, PRRX2 expression (median split) shows **no association with relapse-free survival**: log-rank p = 0.92. Median RFS was actually marginally *longer* in the high-expression group (218.7 vs. 228.2 months — not a real difference given the p-value, just noting the point estimate isn't even in the claimed direction). BCAS2 (also plotted in the paper's Fig. 7, no specific number quoted) shows a non-significant trend the *other* way in ER+ patients (p = 0.085, higher expression trending toward *better* survival: 285.7 vs. 185.4 months), and no significant association in ER- patients either (p = 0.37).

This is an honest non-reproduction, not a bug: METABRIC is a different patient cohort and, more importantly, Kaplan-Meier Plotter's classic breast-cancer tool draws on its own curated compilation of Affymetrix microarray GEO series (not METABRIC specifically), so an exact match was never guaranteed — but the size of the discrepancy (essentially null vs. a paper-reported significant effect) is worth flagging rather than hiding. `expected.json`'s check for this requires a significant (p<0.05) log-rank result — designed to fail on our real measurement — so the benchmark reports this transparently rather than loosening the tolerance to pass.

```bash
igvfagent hou2022-prrx2-survival --gene PRRX2 --output Data/hou2022_erna/verify/prrx2_survival.json
igvfagent hou2022-prrx2-survival --gene BCAS2 --output Data/hou2022_erna/verify/bcas2_survival.json
```

## Part 3 — Fig. 7: BCAS2 knockdown RNA-seq/ChIP-seq

*(GSE201595 RNA-seq + GSE201596 ChIP-seq, BCAS2 siRNA vs. control, MCF-7 +/- E2.)* Implemented as the last stage of `Data/hou2022_erna/pipeline.sh` (`hou2022_bcas2_knockdown` port, `rnaseq` and `chipseq` subcommands); queued behind Part 1 in job `49906827`, not yet run.

## Blocked analyses

- **Fig. 6 (mass-spec proteomics, 88 FERM-interacting proteins)**: the paper's own deposited accession, MassIVE `MSV000087492`, is a **private dataset** (`https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000087492` returns "MassIVE Private Dataset" as of 2026-09-27) — not a bot-protection gate, a genuine access restriction. Blocked `controlled_access`.
- **Figs. 2-4 (CRISPR/dCas9 eRNA targeting: ChIP-qPCR, RNA FISH, RT-qPCR time-courses, in vitro p300 HAT assay)**: wet-lab bar-chart/microscopy readouts with no raw numeric data deposited anywhere machine-readable. Blocked `not_deposited`.

## Running it

```bash
bash Benchmarks/hou2022_estrogen_regulated/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark hou2022_estrogen_regulated
```

## Coverage

See `igvfagent bench report --paper-id hou2022_estrogen_regulated` for the current coverage table; this README is updated as each part completes. 2 of 9 planned analyses are blocked (access); the rest are real work in progress, not placeholders left unfinished silently — see Parts 1/3 above for current status.

## Provenance

`provenance.json` holds the full resolve/harvest record.
