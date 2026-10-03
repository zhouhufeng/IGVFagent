# Schnitzler, Kang et al. 2024 — CAD GWAS convergence onto endothelial cell programs (V2G2P)

[![paper](https://img.shields.io/badge/Nature-628:643--653-blue)](https://doi.org/10.1038/s41586-024-07022-x)
[![PMID](https://img.shields.io/badge/PMID-38326615-blue)](https://pubmed.ncbi.nlm.nih.gov/38326615/)
[![data](https://img.shields.io/badge/GEO-GSE210523-orange)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210523)
[![coverage](https://img.shields.io/badge/reproduction-4%2F11%20analyses-yellow)]()

## Bottom line

`igvfagent bench score --paper-id schnitzler2024_convergence_coronary` → **reproduction: incomplete, 4/11 analyses** (10/10 checks pass).

This is a large, multi-omics Nature paper (CRISPRi Perturb-seq + cNMF gene programs + a bespoke Variant-to-Gene-to-Program pipeline + CRISPRi-FlowFISH validation + protein interaction/structure + zebrafish + bulk RNA-seq under flow). This session reproduced the four analyses whose ground truth is a plain number or fraction, all matched **exactly** to the paper's own stated values, straight from its own deposited Supplementary Tables (not re-derived from raw sequencing):

- **Perturb-seq guide knockdown efficacy** — ported the log2fc/binomial-test formula from Supplementary Table 6's own column legend (no author code exists for this step) and verified it against the authors' own reported log2fc, r=0.99999889, 100% match rate (tolerance 0.1) over 14,758 nonzero-count guides. Confirms the qualitative claim that knockdown is only reliably significant (FDR<0.05) for higher-expressed targets (20.5% of TPM≥100 guides vs ~0% of TPM<10 guides).
- **cNMF gene programs** — the authors' K=60 cNMF run reduces to exactly 50 real programs after excluding 10 batch-associated components (10/60 = 16.7%, matching "we assigned 10 components with Pearson correlation > 0.15 as likely representing batch"); of the 13 programs the authors call endothelial-cell-specific, exactly 5 are significantly CAD-enriched at FDR<0.05 with 2.66–4.63-fold enrichment (paper: "2.6- to 4-fold").
- **V2G2P convergence** — 228 non-lipid-associated CAD GWAS signals (exact denominator match) and 41 unique V2G2P genes (exact), of which 9 have prior mouse-model evidence for atherosclerosis/barrier function (9/41 = 21.95%, exact).
- **V2G2P generalizes to other traits** — reapplying V2G2P to a K562/blood-trait Perturb-seq dataset gives 90 programs, of which 32 are prioritized for at least one GWAS trait (32/90 = 35.6%, exact match to "of the 90 programs... 32 programs were prioritized").

Not yet attempted or not reconciled (7 analyses, still `pending`): ABC enhancer calling from ATAC/ChIP (Table 2), the specific "43 of 228 signals converge on the CCM pathway" figure (my naive gene-signal join gives 58, not reconciled — see below), CRISPRi-FlowFISH validation of the rs1879454 enhancer and TLNRD1 promoter (guide-level MLE data is in hand, Table 28, but element boundaries can't be reliably reconstructed from it — see below), TLNRD1-CCM2 co-IP, AlphaFold2-Multimer structure prediction, the flow/knockdown transcriptional-response overlap (GSE232400), and zebrafish vascular phenotyping.

## Citation

Schnitzler GR, Kang H, Fang S, Angom RS, Lee-Kim VS, Ma XR, Zhou R, Zeng T, Guo K, Taylor MS, Vellarikkal SK, Barry AE, Sias-Garcia O, Bloemendal A, Munson G, Guckelberger P, Nguyen TH, Bergman DT, Hinshaw S, Cheng N, Cleary B, Aragam K, Lander ES, Finucane HK, Mukhopadhyay D, Gupta RM, Engreitz JM. **Convergence of coronary artery disease genes onto endothelial cell programs.** *Nature* **628**:643–653 (2024). DOI [10.1038/s41586-024-07022-x](https://doi.org/10.1038/s41586-024-07022-x) · PMID 38326615 · PMC10921916

## Data and code

| Resource | Identifier |
|---|---|
| GEO superseries | `GSE210523` — Perturb-seq, ATAC-seq, H3K27ac ChIP-seq, RNA-seq in TeloHAEC |
| GEO subseries | `GSE210489` (ATAC-seq), `GSE210491` (ChIP-seq), `GSE210522` (bulk RNA-seq cytokine/CRISPRi), `GSE232400` (bulk RNA-seq flow/MAP3K3), `GSE212396` (pilot scRNA-seq), `GSE210681` (comprehensive Perturb-seq) |
| Authors' code | [EngreitzLab/V2G](https://github.com/EngreitzLab/V2G) (variant-to-gene) and [EngreitzLab/cNMF_pipeline](https://github.com/EngreitzLab/cNMF_pipeline/) (gene-to-program + V2G2P enrichment) — both Snakemake pipelines; neither released was run this session (see below) |
| ABC model | [broadinstitute/ABC-Enhancer-Gene-Prediction](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) |
| Supplementary Tables | 31 tables across 2 Nature ESM files (`MOESM3` = Tables 1–15, `MOESM4` = Tables 16–31), fetched via Nature/Springer's `static-content.springer.com` ESM URLs — not committed to git (copyrighted publisher material; see License below) |
| Full text | PMC10921916, fetched via `eutils efetch db=pmc` (Europe PMC's own fetcher reported this paper closed-access; the NIH manuscript full text was available regardless) |

## Coverage

| Analysis | Paper | IGVFagent | State |
|---|---|---|---|
| Perturb-seq guide knockdown efficacy | guides "effectively knocked down" targets; reliable only for higher-TPM genes | log2fc r=0.99999889 vs authors' own values (100% match, tol 0.1); 20.5% FDR<0.05 at TPM≥100 vs ~0% at TPM<10 | reproduced |
| ABC enhancers in TeloHAEC | ABC score ≥0.015 (Supp. Table 2) | not attempted | pending |
| cNMF gene programs | 50 programs (10/60 batch-excluded); 5 EC-specific programs CAD-enriched 2.6–4x, FDR<0.05 | 50 (exact); 5 (exact), enrichment 2.66–4.63x | reproduced |
| CCM pathway convergence | 43 of 228 non-lipid signals converge on CCM pathway | not reconciled (see below) | pending |
| V2G2P convergence | 228 non-lipid signals; 41 V2G2P genes; 9/41 with prior athero/barrier evidence | 228 (exact); 41 (exact); 9/41=21.95% (exact) | reproduced |
| CRISPRi-FlowFISH, TLNRD1 | rs1879454 enhancer: −21% effect, FDR-corrected P=0.001; promoter: 37 guides, enhancer: 17 guides, 117 negative controls | guide-level MLE data in hand (Table 28); element (guide→enhancer) boundaries not reconstructable from it (see below) | pending |
| TLNRD1–CCM2 co-IP | TLNRD1 interacts with CCM2 | not attempted (imaging/blot readout, no tabular deposit found) | pending |
| TLNRD1–CCM structure (AlphaFold2-Multimer) | predicted interaction interface | not attempted (heavy external tool) | pending |
| Flow transcriptional response | TLNRD1/CCM2 knockdown mimics flow response | not attempted (GSE232400 not yet pulled) | pending |
| Zebrafish vascular phenotypes | tlnrd1/ccm2 knockdown → vascular/cardiac defects | not attempted (Supp. Tables 19–20, 30–31 hold quantified phenotypes, not yet scored) | pending |
| V2G2P generalizes to other traits | of 90 K562 programs, 32 prioritized for ≥1 of 6 blood-related GWAS traits | 90 (exact); 32 (exact) = 35.6% | reproduced |

## Where the reproduction differs, and why

* **"43 of 228 signals converge on CCM/5 programs" not reconciled.** Supplementary Table 1 gives 228 non-lipid signals (exact match) and 41 `CADAssociatedGenes` (exact match), but joining "non-lipid signal × nearby gene where that gene is a CADAssociatedGene" gives 58 signals, not 43. The V2G prioritization in Table 1 has several sub-criteria (`Top2ClosestGeneV2GLink`, `Top2ABCScoreV2GLink`, `CodingVariantGene`) that combine into `CandidateCADGenes` and then `CADAssociatedGenes`; a gene can plausibly sit in more than one signal's ±500kb window without being that signal's own prioritized link. Reconciling this needs the authors' actual V2G pipeline (`EngreitzLab/V2G`), not just its output table — deferred, not fabricated.
* **CRISPRi-FlowFISH element boundaries.** Supplementary Table 28 gives 4-replicate guide-level MLE effect sizes for 684 targeting guides + 117 negative controls tiling a ~157kb window (chr15:81,269,808–81,426,757), but has no guide→element (candidate-enhancer) column. The paper's Methods define elements as MACS2 peaks at a lenient P<0.1 cutoff — those exact peak calls aren't in this table. Gap-based clustering of guide positions was tried and calibrated against the two guide counts the paper states exactly (37 for the TLNRD1 promoter, 17 for the rs1879454 enhancer): no gap threshold recovers both counts simultaneously (best attempt: 29 and 27 respectively, stable across a wide threshold range), meaning simple positional clustering doesn't match the authors' actual peak boundaries closely enough to trust a forced ±21%/P=0.001 check. Reproducing this properly needs re-calling MACS2 peaks on the raw ATAC-seq (`GSE210489`) at the stated lenient threshold — deferred, not forced to a possibly-wrong number.
* **No author code was run this session.** Both `EngreitzLab/V2G` and `EngreitzLab/cNMF_pipeline` were cloned (commits `a3c264f`, `58d066a`) but not executed — they are Snakemake pipelines expecting the full raw Perturb-seq/ATAC-seq/ChIP-seq inputs (tens of GB) and a full external-reference build (1000 Genomes LD, MAGMA, S-LDSC, PoPS). The 4 reproduced analyses instead verify the paper's own deposited *output* tables directly, which is weaker evidence than an independent re-run from raw data but is exact, not approximate.

## Ported into IGVFagent

* `igvfagent perturbseq-guide-knockdown-efficacy` — per-guide target-gene knockdown efficacy (log2fc + binomial-test p-value + BH-FDR) from Perturb-seq pseudobulk counts. No author code exists for this specific step (it's a formula in Supplementary Table 6's own column legend, not in the V2G/cNMF_pipeline repos); `bench verify-port` against the authors' own reported log2fc: r=0.99999889, 100% match rate over 14,758 nonzero-count guides (`Docs/Benchmark/*_schnitzler2024_convergence_coronary_port_perturbseq_guide_knockdown_efficacy/`). Generalizes to any CRISPRi Perturb-seq screen reporting this pseudobulk-count shape. Unreviewed.

## How to reproduce

```bash
# Supplementary tables (not committed — see License) must be fetched first:
mkdir -p Benchmarks/schnitzler2024_convergence_coronary/Supplement
curl -A "Mozilla/5.0" -o Benchmarks/schnitzler2024_convergence_coronary/Supplement/MOESM3.xlsx \
  "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07022-x/MediaObjects/41586_2024_7022_MOESM3_ESM.xlsx"
curl -A "Mozilla/5.0" -o Benchmarks/schnitzler2024_convergence_coronary/Supplement/MOESM4.xlsx \
  "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07022-x/MediaObjects/41586_2024_7022_MOESM4_ESM.xlsx"

python3 Scripts/ported/schnitzler2024_convergence_coronary/perturbseq_guide_knockdown_efficacy.py \
  --input Data/schnitzler2024/tables/table6_guide_knockdown.tsv \
  --out Benchmarks/schnitzler2024_convergence_coronary/work/guide_knockdown
python3 Benchmarks/schnitzler2024_convergence_coronary/table_counts.py Docs/Benchmark/<timestamp>_schnitzler2024_convergence_coronary_table_counts
igvfagent bench score --paper-id schnitzler2024_convergence_coronary
```

Scripts: `table_counts.py` (cNMF/V2G2P counts and fractions, straight off Supplementary Tables 1, 13, 15, 22). The guide-knockdown port lives at `Scripts/ported/schnitzler2024_convergence_coronary/perturbseq_guide_knockdown_efficacy.py` (registered as `igvfagent perturbseq-guide-knockdown-efficacy`).

## Next steps (not done this session)

1. Reconcile the 43/228 CCM-convergence count against the actual `EngreitzLab/V2G` prioritization logic (run the pipeline, or read its source for the exact combination rule).
2. Re-call ATAC-seq peaks from `GSE210489` at MACS2 P<0.1 to reconstruct real CRISPRi-FlowFISH element boundaries, then score the rs1879454 enhancer (-21%, P=0.001) and TLNRD1 promoter checks against Table 28.
3. Pull `GSE232400` (flow/MAP3K3 bulk RNA-seq) and cross-reference against Supplementary Tables 19/29 for the flow-mimicry claim.
4. Score the zebrafish phenotypes (Supplementary Tables 19–20, 30–31 hold quantified data, e.g. ventricular wall thickness, vascular permeability) directly from their own tables, same approach as the cNMF/V2G2P counts above.
5. TLNRD1–CCM2 co-IP and AlphaFold2-Multimer structure prediction are lower priority: the former is a blot/imaging readout with no obvious tabular ground truth, the latter needs a full AlphaFold3/Multimer run (heavy, external).

## License + provenance

* **Data**: GEO GSE210523 and subseries (public), fetched at run time, never redistributed.
* **Supplementary Tables**: Nature/Springer ESM (`MOESM3.xlsx`, `MOESM4.xlsx`), copyrighted publisher material — `Benchmarks/schnitzler2024_convergence_coronary/Supplement/` is gitignored, never redistributed. Re-fetch with the commands above.
* **Full text**: PMC10921916 (`Data/schnitzler2024/supplement/pmc10921916.xml`), gitignored under `Data/*`, fetched for reference only.
* **Code**: IGVFagent Apache-2.0; the authors' V2G/cNMF_pipeline repos are cited as provenance but were not executed this session (see above).
