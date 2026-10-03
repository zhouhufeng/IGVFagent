# kedmi2022_instructs_microbiota

## Paper

**A RORγt+ cell instructs gut microbiota-specific Treg cell differentiation.**
Kedmi R, Najar TA, Mesa KR, Grayson A, Kroehling L, Hao Y, Hao S, Pokrovskii M, Xu M, Talbot J, Wang J, Germino J, Lareau CA, Satpathy AT, Anderson MS, Laufer TM, Aifantis I, Bartleson JM, Allen PM, Paidassi H, Gardner JM, Stoeckius M, Littman DR.
*Nature* 2022 · doi:[10.1038/s41586-022-05089-y](https://doi.org/10.1038/s41586-022-05089-y) · PMID 36071167 · PMC9908423 · IGVF award IGVF0069, grant UM1HG012076

Resolver confidence: **1.00** (exact PMID match).

The paper (Littman lab) identifies which antigen-presenting cells (APCs) instruct gut-microbiota-specific regulatory T cell (iTreg) differentiation vs. inflammatory Th17 differentiation in mice colonized with *Helicobacter hepaticus*. Using genetic fate-mapping, conditional knockouts, and CITE-seq, it shows antigen presentation by RORγt+ cells — type 3 innate lymphoid cells (ILC3) and a rare Aire+RORγt+ population called "Janus cells" (JC) — rather than classical dendritic cells, is required and sufficient for iTreg induction, and depends on MHCII, the chemokine receptor CCR7, and the TGF-β activator αv integrin in those RORγt+ cells.

Europe PMC's `fullTextXML` endpoint returns HTTP 500 for this record (PMC9908423 is an NIHMS author manuscript, not the publisher's own open-access XML), so the full text was fetched directly from `pmc.ncbi.nlm.nih.gov` (confirmed 2026-09-27) and used to build `harvest.json` by hand.

## Why this doesn't use an IGVFagent assay route

`igvfagent bench route` finds no match: this is a CITE-seq (RNA+ADT+HTO) WNN clustering paper plus a separate 5' scRNA-seq dataset, not one of IGVFagent's canned assay families (SPLiT-seq, cell hashing alone, ChIP-Atlas, ENCODE FCE, Synapse, figshare/Zenodo, or IGVF Portal discovery). Scaffolded with the `portal_discovery` fallback template, then hand-edited.

## Code and data

- **Code**: the paper's own repository, [`nygctech/Kedmi-CITEseq`](https://github.com/nygctech/Kedmi-CITEseq) (`kedmi.et.al.R`, a Seurat WNN pipeline for Fig. 2a-b's CITE-seq dataset; `Aire scRNAseq.ipynb`, a scanpy-based notebook for Extended Data Fig. 7-8's Aire-reporter scRNA-seq, which imports a private custom package `scrnatools` not available here).
- **Data**: GEO SuperSeries `GSE200148` (`Series_pubmed_id = 36071167`), composed of:
  - `GSE190372` — CITE-seq + cell hashing of MLN cells from Hh-colonized `tdTomato-ON-CD11c` fate-map mice (2 lanes: cDNA/ADT/HTO 1 and 2; processed 10x/kallisto count matrices, not raw reads).
  - `GSE200147` — 5' scRNA-seq of Aire(Adig)-GFP-sorted vs. unsorted MLN cells (2 samples, each a processed cellranger `.h5` filtered matrix).

Both datasets are small, processed, and public — no alignment or raw-read processing needed, unlike most of this suite's from-raw-reads benchmarks.

## Part 1 — Fig. 2a-b / Extended Data Fig. 2: CITE-seq WNN identifies an ILC3/Janus-cell population

**Reproduced.** New port `kedmi2022_citeseq_wnn` (`igvfagent kedmi2022-citeseq-wnn`) is a faithful scanpy re-implementation of the authors' `kedmi.et.al.R` Seurat WNN pipeline: per-lane QC (the paper's own thresholds: `nFeature_RNA` 700-5000, `percent.mito<6`, `nCount_ADT<5000`, `nCount_HTO<4000`), hashtag demultiplexing (median+3·MAD background threshold per tag — the same statistical intent as Seurat's `HTODemux`, not its exact k-means implementation), CLR-normalized ADT, and Leiden clustering on concatenated RNA+ADT PCA space (a documented simplification of Seurat's per-cell modality-weighted WNN, not a silent substitution).

Run independently from the raw deposited `GSE190372` matrices (no author-computed reference table exists to diff against), it finds one cluster (of 21 at Leiden resolution 1.0, 4291 cells total after QC/demux) with:
- mean *Rorc* = **0.81** (log1p-normalized), a **31x enrichment** over the mean of all other clusters (0.026) — a clean ILC3 signature (this cluster also has elevated *Il7r* = 2.81 and *Kit* = 1.11, the canonical ILC3 markers), and
- anti-IA/IE (MHCII) ADT signal in the **top quartile** (75th percentile) of all 21 clusters,

matching the paper's own qualitative claim: *"we identified both type 3 innate lymphoid cells (ILC3) and a recently described Aire+ RORγt+ population, named Janus cells (JC), among the tdTomato+ cells, and these also expressed MHCII."* The paper gives no exact number for this figure (a UMAP + violin plot), so both checks are class-A range checks against data-derived thresholds set well above background, not tuned to the measured value.

## Part 2 — Extended Data Fig. 7-8: Aire-reporter scRNA-seq identifies the Janus cell population directly

**Reproduced.** New port `kedmi2022_etac_janus_scrnaseq` (`igvfagent kedmi2022-etac-janus-scrnaseq`) reruns the intent of the authors' `Aire scRNAseq.ipynb` (standard scanpy QC/HVG/PCA/Leiden, since the notebook's own `scrnatools` package is private and unavailable) on the real deposited `GSE200147` matrices (`GSM6012929` Adig-GFP+-sorted, `GSM6012930` unsorted reference; 21,343 cells total after QC, 26 clusters at resolution 1.0).

It finds one cluster with mean *Aire* = **0.93** (a **15.5x enrichment** over the mean of all other clusters, 0.06) that co-expresses *Rorc* (mean 0.53) — exactly the paper's "Aire+RORγt+" Janus-cell definition — and is **91% composed of GFP+-sorted cells**, vs. an even ~50/50 input mix, i.e. strongly enriched by the Aire-reporter sort exactly as the paper's gating strategy intends. A separate Rorc-high/Aire-low cluster (interpreted as ILC3; only 18% GFP+-sort composition, near the unsorted baseline) expresses *Itgav* (mean 0.57); *Itgb8* reads as comparatively JC-specific (0.449 in the JC cluster vs. 0.015 in the ILC3 cluster) — consistent with the figure's title, "Itgav and Itgb8 expression in ILC3 and JC."

## Blocked analyses (not_deposited)

Four analyses — `fig1_dc_ccr7_dependence`, `fig2ef_mhcii_rorgt_requirement`, `fig3_ilc3_jc_identity_ccr7`, `fig4_integrin_tgfb_activation` — cover the paper's flow-cytometry mouse-genetics figures (Hh-specific T cell proliferation/differentiation percentages, Foxp3/RORγt/T-bet gating, cytokine production, colon histology, across many conditional-knockout mouse lines). These report only aggregate percentages and representative flow plots in the figure panels. No "Source Data" statement, FlowRepository accession, or any other machine-readable per-cell/per-sample numeric deposit exists anywhere in the full text (confirmed by full-text search for "Source Data" / "FlowRepository" / "Supplementary Table", all absent). The paper's only accessions (`GSE200148` SuperSeries) cover the two sequencing-based figures above, not these wet-lab functional assays. Blocked `not_deposited`, not attempted.

A fifth potential item, Fig. 5 (a summary model diagram integrating the paper's genetic and CITE-seq findings), carries no independent numeric or data-backed claim of its own and was dropped from the plan rather than padding the analysis count.

## Running it

```bash
bash Benchmarks/kedmi2022_instructs_microbiota/run.sh
igvfagent bench score --paper-id kedmi2022_instructs_microbiota
```

## Coverage

`igvfagent bench report --paper-id kedmi2022_instructs_microbiota`: **2/6 analyses reproduced** (4/4 confirmed checks pass), **reproduced_except_access** — the remaining 4 analyses are all genuinely `not_deposited` (real, checked evidence above), not partial attempts.

## Provenance

`provenance.json` holds the full resolve/harvest record.
