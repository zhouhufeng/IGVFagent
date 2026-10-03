# rohm2025_activation_imprinted

## Paper

**Activation of the imprinted Prader-Willi syndrome locus by CRISPR-based epigenome editing**
Rohm D, Black JB, McCutcheon SR, Barrera A, Berry SS, Morone DJ, Nuttle X, de Esch CE, Tai DJC, Talkowski ME, Iglesias N, Gersbach CA.
*Cell Genomics* 2025 · doi:[10.1016/j.xgen.2025.100770](https://doi.org/10.1016/j.xgen.2025.100770) · PMID 39947136 · PMC11872474 · IGVF award IGVF0192, grant UM1HG012053

Resolver confidence: **1.00** (exact PMID match).

The paper runs three genome-scale CRISPR tiling screens across the imprinted Prader-Willi syndrome (PWS) locus in an isogenic SNRPN-2A-GFP iPSC reporter system — CRISPRi (dCas9-KRAB) on the paternal (normally active) allele, and two activator screens (VP64-dCas9-VP64, and the DNA-demethylase Tet1c-dCas9) on the maternal (imprinted/silenced) allele — to map the regulatory elements controlling SNRPN, then validates individual hits with isogenic ΔPWS iPSCs (RNA-seq, ATAC-seq, CUT&RUN, bisulfite sequencing) and iPSC-derived neurons.

## Why this doesn't use an IGVFagent assay route

`igvfagent bench route` scores `perturb_catalog` highest (a real signal — this genuinely is a CRISPR-screen paper), but IGVFagent's built-in `crispr-screen`/`gradient-screen` CLI tools are IGVF-**Portal**-specific (they resolve a `FileSet` accession via the Portal API and group its sorted bins by `construct_library_sets`). This paper's screens are deposited on **GEO**, not the IGVF Portal, so those tools return `not found` on GSE285285 etc. — a real "doesn't fit" case, not a bug. The paper's own words: *"Data generated in high-throughput sequencing studies are available through NCBI GEO... Any additional information required to reanalyze the data... is available from the lead contact upon request"* — **no author analysis code exists** for the screen hit-calling, RNA-seq, ATAC-seq, or CUT&RUN steps either.

## Part 1 — Fig. 1: the three CRISPR tiling screens — DONE (3/3 reproduced)

Real data throughout:
- **The paper's own gRNA library** (Table S1/S2 of the paper, `mmc2.xlsx`) — 11,751 gRNAs (531 non-targeting) for the primary screens, 592 for the Tet1c sub-library — with real hg19 genomic coordinates. Unlike the NIHMS-hosted supplementary files elsewhere in this suite, **Cell Genomics' own supplementary files have no bot-gate**; EuropePMC's open bundle serves them directly.
- **Real raw gRNA-seq reads**, all 18 samples (3 replicates × {high-GFP, low-GFP} sorted tail × 3 screens), from GEO GSE285285 (paternal dCas9-KRAB), GSE285293 (maternal VP64dCas9VP64), GSE285289 (Tet1c-dCas9 sub-library).
- **Counted** with the new port `rohm2025_grna_count` (exact 20-mer substring match against the real library — robust to the unpublished amplicon read layout).
- **Hit-called with DESeq2** (`rohm2025_screen_hits` + `rohm2025_deseq2_hits.R`) — the paper's own named tool ("−log10(padj), where padj is the multiple-hypothesis-corrected p value from DESeq2"), comparing the high-GFP vs. low-GFP sorted tails, padj < 0.05.
- **Clustered into regions** (>=3 adjacent significant guides within 1 kb — our own documented choice; no author region-calling code exists) and checked against the paper's own stated region counts via `bench verify-port`.

| Screen | Real hit gRNAs (padj<0.05) | Regions found | Paper's claim | Match |
|---|---|---|---|---|
| Paternal dCas9-KRAB | 393 | **7** | 6 marked regions (pat1–pat6) | close (rtol 0.2) — dominant cluster (up to 308 hits) spans chr15:25.20–25.24 Mb, squarely the SNRPN gene body/PWS-IC |
| Maternal VP64dCas9VP64 | 11 | **2** | 2 distinct clusters (mat1, mat2) | **exact** |
| Maternal Tet1c-dCas9 sub-library | 54 | **1** | mat3 region at PWS-IC | **exact** |

Two independent, coordinate-level confirmations beyond the region *counts* above:
- Our **mat1** cluster (chr15:25,074,435–25,074,809) falls exactly inside the paper's own **4 individually-validated mat1 gRNA coordinates** (Table S1 "gRNA validations" sheet: mat1 g1–g4, chr15:25,074,051–25,074,809) — a real coordinate-level match, not just a count.
- Our **Tet1c mat3** cluster (chr15:25,199,451–25,201,415) is disjoint from VP64's mat1/mat2 (confirming *"no overlap in the gRNA hits between the maternal and paternal screens"* / *"not detected with VP64dCas9VP64"*) and sits **inside** the paternal dCas9-KRAB screen's own dense hit region (chr15:25,197,530–25,225,216) — confirming the paper's separate claim that Tet1c's mat3 hits *"largely overlapped with significant hits in the pat4 region from the dCas9KRAB screen."*

One real data-quality issue, reported honestly: `Tet1c sorted low GFP rep1` (GSM8700018 / SRR31810701) has only **20 reads** deposited — an effective replicate loss, not a bug in our pipeline. The Tet1c screen above uses the 2 remaining low-GFP replicates.

```bash
bash Benchmarks/rohm2025_activation_imprinted/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark rohm2025_activation_imprinted
```

## Parts 2/3 — Figs. 2–3: RNA-seq specificity, ATAC-seq, CUT&RUN — not yet started

Real GEO deposits exist for all three (GSE243185 RNA-seq, GSE285305 ATAC-seq, GSE285284 H3K4me3 CUT&RUN, all raw reads) and are not blocked — just not yet implemented in this pass. See `OPERATIONS.md` for what each needs.

## Not attempted (wet-lab, no sequencing deposit)

Fig. 2/4's qPCR-based transcript-induction percentages (SNRPN ~7%/25% of WT, SNORD116 ~30%, etc.), HCR-FlowFISH single-cell distributions, bisulfite-sequencing methylation medians, and cell-proliferation/differentiation assays are all low-throughput validation experiments reported as summary plots, not deposited as raw re-analyzable data (GEO here only covers the genome-wide screens/RNA-seq/ATAC-seq/CUT&RUN). Not formally blocked in `expected.json` (no analysis was planned for them) since they were never planned as computational targets in the first place.

## Coverage

**3 of 7** planned analyses reproduced (10/10 checks pass) — all three CRISPR screens (Fig. 1). The remaining 4 (Figs. 2–3) need the RNA-seq/ATAC-seq/CUT&RUN pipelines built.

## Provenance

`provenance.json` holds the full resolve/harvest record. Ports: `rohm2025_grna_count`, `rohm2025_screen_hits` (`igvfagent port list`), pinned to DESeq2's own GitHub mirror commit (no author code exists to pin to instead — see OPERATIONS.md sec 4).
