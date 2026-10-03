# Kaplan et al. 2024 — CRISPRi screen finds the ONECUT1e-664kb enhancer

[![paper](https://img.shields.io/badge/Cell%20Reports-2024-blue)](https://doi.org/10.1016/j.celrep.2024.114640)
[![PMID](https://img.shields.io/badge/PMID-39163202-blue)](https://pubmed.ncbi.nlm.nih.gov/39163202/)
[![status](https://img.shields.io/badge/IGVFagent%20live%20concordance-11%2F11%20checks%2C%204%2F4%20analyses%20reproduced-success)]()

## Bottom line

Kaplan, Wong, Yan, Pulecio et al. (Huangfu lab, MSKCC), **"CRISPR screening uncovers a long-range enhancer for ONECUT1 in pancreatic differentiation and links a diabetes risk variant"**, *Cell Reports* 2024. DOI: [10.1016/j.celrep.2024.114640](https://doi.org/10.1016/j.celrep.2024.114640) (preprint: [10.1101/2024.04.26.591412](https://doi.org/10.1101/2024.04.26.591412)) · PMID: 39163202 · PMCID: PMC11406439.

**This paper states "This paper does not report original code," and its GEO deposit (GSE267330) does not include the primary CRISPRi screen's raw sgRNA counts** — only the downstream validation assays (RNA-seq, H3K27ac ChIP-seq, Hi-C, long-read phasing). Confirmed directly: GSE267330's SuperSeries relations list exactly 4 subseries (GSE261390 RNA-seq [18 samples], GSE261391 long-read-seq [1], GSE261392 H3K27ac ChIP-seq [12], GSE263174 Hi-C [2] = 33 GSM), and all 33 SRA experiments under the paper's own BioProject PRJNA1111076 (via NCBI eutils) map 1:1 onto those same 4 subseries — none of them is the screen.

**What unblocks this paper anyway: the authors' own already-computed screen and DESeq2 results are deposited as Cell Reports supplementary tables**, downloaded directly from ScienceDirect's asset CDN (`ars.els-cdn.com`, no paywall or bot-gate, unlike the PMC/NIHMS mirror some other Cell Press papers use) — Table S1 (region selection + MAGeCK-RRA screen hits), Table S3 (RNA-seq DESeq2), Table S4 (H3K27ac ChIP-seq DESeq2 + allele-specific counts). IGVFagent re-derives the paper's headline numbers directly from these tables and cross-checks two against live public APIs (Ensembl, gnomAD v4). All 4 planned analyses reproduce; 11/11 checks pass.

## What IGVFagent reproduces

| Analysis | Source | IGVFagent (re-derived) | Paper's own claim | Verdict |
|---|---|---|---|:---:|
| CRISPRi screen hits (Fig. 1E) | Table S1 (mmc2.xlsx), Sheet 5 | **38** top hits, **15** linked enhancer-gene pairs, distances **2–664 kb** | "38 top hits...15 enhancer-gene pairs...2 to 664 kb" | ✓ exact |
| ONECUT1 loss, homozygous eKO, PFG stage (Fig. 2/3) | Table S3 (mmc4.xlsx), Sheet 3 — DESeq2 | **96.05%** loss (log2FC=-4.66, padj=3.2e-95) | ">95% loss" (ddPCR) | ✓ |
| ONECUT1 loss, heterozygous eKO, PFG stage | Table S3 (mmc4.xlsx), Sheet 2 — DESeq2 | **53%** loss (log2FC=-1.09, padj=0.0049) | intermediate/allele-dosage effect expected | ✓ |
| H3K27ac allelic imbalance, heterozygous eKO (Fig. 3) | Table S4 (mmc5.xlsx), Sheet 3 — 86 phased SNPs, 21 top-differential peaks | mean deviation from allelic balance: **0.74** (top peaks, Het) vs **0.24** (WT) | "H3K27ac was decreased in cis to the allele with enhancer deletion" | ✓ direction confirmed |
| rs528350911 location + allele frequency (Fig. 4) | Table S1 coordinates + Ensembl REST + gnomAD v4 API (live) | position chr15:53,455,031 falls inside re-derived ONECUT1e-664kb (chr15:53,454,536–53,455,646); MAF **0.0068** (NFE), **0.0039** (all pops) | "directly within ONECUT1e-664kb"; MAF 0.0068 / 0.0039 (gnomAD v4.1.0) | ✓ exact |

**Independent geometric cross-check:** the T2D risk variant rs528350911's own GRCh38 position (chr15:53,455,031, from a live Ensembl REST lookup, entirely independent of the paper) falls inside the ONECUT1e-664kb coordinates this benchmark itself re-derived from Table S1 (chr15:53,454,536–53,455,646) — confirming the two datasets (the screen hit list and the population-genetics variant record) really do refer to the same genomic element.

## Honest caveats

* **The primary CRISPRi screen (Fig. 1) is not reproduced from raw reads.** No raw sgRNA-seq counts are deposited anywhere findable (GEO/SRA/Zenodo/4DN) and the paper reports no original code, so there is no MAGeCK pipeline to re-run. What IS reproduced is the authors' own already-computed MAGeCK-RRA hit list (Table S1), re-tallied and cross-checked against the paper's stated summary numbers — the same "verify against the authors' own already-scored output" pattern IGVFagent uses elsewhere when raw data isn't available (e.g. `martyn2025_variant_flowfish`).
* **The H3K27ac "in cis" claim is checked directionally, not by exactly reproducing the paper's own statistical test.** IGVFagent computes a simple mean-absolute-deviation-from-1.0 allelic-imbalance metric from the authors' own Table S4 Sheet 3 phased-SNP counts; it shows the same qualitative pattern the paper describes (heterozygous cells are far more allele-imbalanced at the top differential peaks than WT cells are anywhere), but this is not the paper's own DESeq2-based peak-level significance test (that's in Table S4 Sheet 2, not used here).
* **The fine-mapping PPA=0.994 itself is cited, not independently re-derived.** It comes from Mahajan et al.'s published T2D GWAS meta-analysis (SuSiE fine-mapping), per the paper's own citation — re-running SuSiE on the raw GWAS summary statistics was out of scope for this pass. What IS independently checked: the variant's genomic containment within the enhancer (geometric) and its two stated gnomAD allele frequencies (both matched a live gnomAD v4 API query almost exactly).
* **Table S2 (mmc3.xlsx: oligo/antibody/ddPCR details) and the Hi-C and long-read-sequencing GEO subseries were not scored** — they support Methods/validation rather than carrying an independent headline number.

## How to reproduce

```bash
bash Benchmarks/kaplan2024_crispr_screening/run.sh
python3 Benchmarks/concordance.py --benchmark kaplan2024_crispr_screening
```

Downloads the paper's 4 Cell Reports supplementary tables (`mmc2.xlsx`-`mmc5.xlsx`, ~5-20 MB each) directly from ScienceDirect's asset CDN, then runs `analyze_kaplan2024.py`, which also makes two small live API calls (Ensembl REST, gnomAD v4 GraphQL) for the variant cross-check. Writes `summary.json` under `Data/Benchmarks/kaplan2024_crispr_screening/runs/<ts>_kaplan2024_crispr_screening/`.

## License + provenance

* **Data**: Cell Reports supplementary materials (Elsevier, CC-BY per journal policy); Ensembl and gnomAD public APIs.
* **Paper code**: none — the paper states "This paper does not report original code."
* **IGVFagent code**: Apache-2.0; `Benchmarks/kaplan2024_crispr_screening/analyze_kaplan2024.py` (new, written for this benchmark).
* **Citation**: Kaplan SJ et al. *Cell Reports* (2024). doi:10.1016/j.celrep.2024.114640 · PMID:39163202 · PMCID:PMC11406439
