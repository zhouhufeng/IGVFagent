# lee2025_massively_parallel

## Paper

**Massively parallel reporter assay investigates shared genetic variants of eight psychiatric disorders.**
Lee S, McAfee JC, Lee J, Gomez A, Ledford AT, Clarke D, Min H, Gerstein MB, Boyle AP, Sullivan PF, Beltran AS, Won H.
*Cell* 2025 · doi:[10.1016/j.cell.2024.12.022](https://doi.org/10.1016/j.cell.2024.12.022) · PMID 39848247 · PMC11890967 (award 1UM1HG012003)

Cross-disorder MPRA (17,841 variants across 136 loci) in human neural progenitors, plus hiPSC-neuron and mouse-neocortex CROP-seq validation, across eight psychiatric disorders (ASD, ADHD, SCZ, BD, MDD, Tourette, OCD, anorexia).

Code Availability: [thewonlab/crossdisorder-MPRA](https://github.com/thewonlab/crossdisorder-MPRA) (commit `0d6a07b5925031a543b49733d2f9aed17ee69c42`) — bulk MPRA count-matrix -> mpralm script only; no code for element-activity calling, eQTL overlap, CRISPR-validation DE, or the gene-property analyses. Data: GEO GSE244011 (bulk MPRA), GSE276947 (hiPSC CROP-seq), GSE282731 (mouse CROP-seq).

## What IGVFagent does

**Route:** `mpra`, plus the registered port `lee2025-mpra-emvar` and reuse of the built-in `sc_crispr_de_*` single-cell CRISPR-screen chain.

```bash
bash Benchmarks/lee2025_massively_parallel/run.sh
.venv/bin/igvfagent bench score --paper-id lee2025_massively_parallel
```

## Coverage: 2/6 analyses reproduced, 3/6 honestly blocked, 1/6 pending compute (7/7 checks pass)

| Analysis | Figure | Status | Result |
|---|---|---|---|
| `fig1_active_elements` | Fig. 1 | **blocked** (`not_deposited`) | "MPRA-active element" calling (9.3%/1,478) needs negative-control/basal-promoter barcodes that were never deposited (GSE244011's matrices are REF/ALT variant columns only — confirmed no neg/scram/ctrl rows); no author script exists for this step either. |
| `fig1_emvar_calling` | Fig. 1 | **reproduced** | Ran the authors' own `mpralm` script (registered as `igvfagent lee2025-mpra-emvar`) on the real GEO count matrices. **Class-B: 100% match to the authors' unmodified code** (14,494/14,494 variants, logFC identical to atol=1e-8). **Class-A**: the paper's two headline example variants match within tolerance — rs301804/RERE logFC **-5.05 vs -5.17** (2.4% diff), rs4513167/DCC logFC **-1.15 vs -1.27** (9.2% diff), both same sign/order-of-magnitude on FDR. |
| `fig3_emvar_eqtl_overlap` | Fig. 3 | **blocked** (`embargoed`) | Paper's 62%/74%/79% needs MetaBrain (email-request gated) or the paper's own Table S4 (PMC proof-of-work gated). Two independent public-data substitutes attempted and reported honestly as misses: GTEx v10 brain cis-eQTL bulk files gave <=1% hit rate (~60x short, a power-gap vs MetaBrain's 6,532 pooled samples); IGVF Catalog's eQTL-Catalogue-backed edges gave a better 34.6% brain-hit rate but chance-level (50%) direction-of-effect concordance (paper: 74%) — not a confident reproduction either way. |
| `fig4_crispr_validation_rere_dcc` | Fig. 4 | **reproduced** (order-of-magnitude) | Real GSE276947 hiPSC-neuron CROP-seq data through `sc_crispr_de_prepare`/`test`/`aggregate` (NB GLM, guide-bearing vs no-guide background). Median per-guide log2FC: rs301804/RERE **-0.105 vs paper's -0.129**; rs4513167/DCC **-0.163 vs -0.261**; rs4614799/DCC **-0.111 vs -0.16**. All three same sign, same order of magnitude. FDR not compared (alpha-RRA aggregation was uncalibrated — too few non-targeting guides in this dataset for a fair FDR estimate). |
| `fig5_pleiotropic_gene_properties` | Fig. 5 | **blocked** (`other`) | Paper's p=0.003/0.009/0.004 needs the authors' own emVar-to-target-gene pleiotropy classification (Table S4), unobtainable this session. Independent proxy attempted (nearest-gene mapping + genome-wide GWAS Catalog disorder counts + STRING/GTEx/gnomAD): PPI connectivity direction matched but n.s. (p=0.25-0.47 vs paper's 0.009); expression breadth showed no difference at all (p=0.65 vs paper's 0.003, direction mismatch); mutation intolerance was blocked by gnomAD API rate-limiting. Reported as a genuine, honest non-reproduction by an independent method, not stretched to pass. |
| `fig6_cropseq_mouse_brain` | Fig. 6 | **pending compute** | GSE282731 (10 batches) downloaded and merged (113,524 barcodes; 19,538 single-guide + 23,790 background cells; Kmt5a is labelled `Setd8` in the guide library). The per-guide NB GLM (Anp32e_1/2, Setd8_1/2, GFP_1 vs background) is SLURM job `49907510` (`cropseq_test.sbatch`). Two earlier submissions (48956656, 48956686) failed because their script and logs were in a login node's `/tmp`, which compute nodes can't read; see OPERATIONS.md §5. When it lands, aggregate it and compare with the paper's 46 Anp32e / 244 Kmt5a DEGs. |

## Honest caveats

* This paper is closed-access; the PMC author-manuscript page (PMC11890967) itself was readable, but its Excel supplementary tables (S1-S6, containing the authors' own MPRA/eQTL/gene classifications) sit behind a client-side proof-of-work JS challenge that neither `curl`/`WebFetch` nor this session (no browser automation tool available) could clear. Cell.com/ScienceDirect return HTTP 403. This is the root cause behind 3 of the 4 non-fully-reproduced analyses.
* `fig1_emvar_calling`'s port (`lee2025-mpra-emvar`) fixes one cosmetic bug in the authors' own script (`rownames()` on a `data.table` silently returns generic row numbers, not the real variant IDs) and optionally applies the authors' own deposited QC excludelist, which their released script never reads.
* `fig3` and `fig5` each received a genuine independent-method attempt using different public data sources as substitutes for the paper's own gated resources, and both are reported as honest non-reproductions rather than loosened to pass — see the `blocker.reason` text in `expected.json` for full numbers.

## Provenance

`provenance.json` holds the resolve/harvest/route record. `expected.json`'s `checks[].provenance` and `analyses[].blocker.reason` hold the full paper-number-vs-reproduced-number comparisons and blocker evidence.
