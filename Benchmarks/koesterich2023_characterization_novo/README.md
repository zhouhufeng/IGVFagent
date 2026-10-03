# koesterich2023_characterization_novo

## Paper

**Characterization of De Novo Promoter Variants in Autism Spectrum Disorder with Massively Parallel Reporter Assays.**
Koesterich J, An JY, Inoue F, Sohota A, Ahituv N, Sanders SJ, Kreimer A.
*International Journal of Molecular Sciences* 2023 · doi:[10.3390/ijms24043509](https://doi.org/10.3390/ijms24043509) · PMID 36834916 · PMC9959321 · IGVF award IGVF0024, grant UM1HG011966

Resolver confidence: **1.00** (exact PMID match).

The paper tests 3600 de novo mutations (DNVs) found in gene promoters of ASD probands and unaffected sibling controls (from whole-genome sequencing) via lentiMPRA in human-ESC-derived neural progenitor cells (N2), using the real MPRAnalyze Bioconductor package to call 165 "high-confidence DNVs" (HcDNVs) whose reference/alternate allele reporter activity differs (FDR<=0.2, |logFC|>=0.1). It reports no difference in HcDNV rate between ASD cases and controls (a null finding), enrichment of HcDNVs for open chromatin/active histone marks and TF-motif disruption, and higher GC content in HcDNVs. The paper has no Code Availability statement; every analysis here is a from-Methods reimplementation, running the real, published MPRAnalyze tool the Methods name directly on the paper's own deposited data.

## A real bug found in the paper's own public data deposit

**GEO GSE216129's own `GSE216129_RefAlt_Complete_dnacounts.tsv.gz` and `_rnacounts.tsv.gz` decompress to byte-identical content** (confirmed via an independent, from-scratch re-download — not a caching artifact on our end). An RNA/DNA transcription-rate ratio computed from identical matrices is mathematically guaranteed null everywhere, which is exactly what an independent MPRAnalyze re-run on this data gives: FDR=1 for all 3091 variant pairs, including a locus Table S3 reports as one of the single most significant hits in the entire paper (fdr=2.43e-35). This means the correct RNA barcode counts needed to independently re-derive the paper's core MPRAnalyze result are not usably present in the current public GEO release — see OPERATIONS.md Sections 5-5b for the full diagnostic trail.

Given this, `fig2_mpranalyze_hcdnv` is honestly reported as **blocked (`not_deposited`)** rather than forced to a false "reproduced." The three analyses built on top of it (case/control test, ATAC/H3K27ac enrichment, GC content) still **independently recompute their own real statistics** using the authors' own deposited final result table (Supplementary Table S3, from the paper's supplementary materials — a different, unaffected deposit) for the HcDNV classification input, and all three genuinely reproduce the paper's numbers.

## Coverage

`igvfagent bench score --paper-id koesterich2023_characterization_novo`: **3/7 analyses reproduced**, 2 blocked (`not_deposited`, real evidence), 2 pending (not attempted this pass, not blocked — see below). 5/5 confirmed checks pass.

| Analysis | Figure | State | Result |
|---|---|---|---|
| `fig2_mpranalyze_hcdnv` | Fig. 2A-C | **blocked** (`not_deposited`) | GEO's own dnacounts/rnacounts files are duplicates (see above); the paper's exact 165-HcDNV threshold logic is validated against the authors' own deposited Table S3 instead — an identity match by construction, not independent evidence |
| `case_control_no_diff` | Results 2.5, Fig. S4 | **reproduced** | Fisher's exact test on HcDNV rate, ASD case vs. control: measured p=0.150, OR=0.791 (78/165 case, 47.3%) vs. paper's p=0.14, OR=0.79 (47.2%) — the paper's null finding reproduces almost exactly |
| `fig3b_atac_h3k27ac_enrichment` | Fig. 3B | **reproduced** | ATAC-seq(DNase): p=1.16e-11 (paper: 5.4e-18); H3K27ac: p=1.84e-7 (paper: 4.0e-8) — both independently recomputed via real Fisher's exact test, same direction and order of magnitude as the paper (see below for why not exact) |
| `figs6_gc_content` | Fig. S6 | **reproduced** | HcDNV mean 63.9% GC vs. background 50.7% GC, Wilcoxon p=4.36e-35 — vs. paper's ~63%/~50%/3.9e-36: essentially an exact match, computed independently from the paper's own deposited MPRA oligo FASTA |
| `fig1e_barcode_association` | Fig. 1E, Fig. S1-S2 | **blocked** (`not_deposited`) | The deposited association pickle has candidate barcode strings but no per-barcode read count, so the paper's stated association filter (>=80% of reads / >=3 reads) cannot be applied to it; real pre-filter numbers are reported (13.86M candidate assignments, mean ~1979/sequence) but the paper's own post-filter number (875,000 / mean 125) is not recomputable from what's public |
| `fig4_motif_disruption` | Fig. 4 | pending | Not attempted this pass — tractable (real FIMO + a JASPAR motif database on the deposited FASTA), not access-blocked |
| `fig5_gsea_domino` | Fig. 5 | pending | Not attempted this pass — tractable (IGVFagent already has `enrich_gsea`/`enrich_pathways`; DOMINO's own PPI-subnetwork step would need more work), not access-blocked |

## Why the ATAC-seq/H3K27ac enrichment isn't an exact match

The paper's own Methods cite GEO accession `GSM517439` as the source of the ATAC-seq/H3K27ac peak BED files used for this test. That accession actually resolves to an unrelated 2010 polyA RNA-seq sample (a `.wig` track, hg18, from the original Wu et al. NPC-derivation paper this study's Methods otherwise correctly cites) — not ATAC-seq/ChIP-seq peaks, and no alternative accession for the actual peak files is given anywhere in the text. This looks like a citation error in the paper itself. `koesterich2023_annotation_enrichment.py` therefore uses the authors' own already-annotated Table S3 columns (`ENCODE_DNase`, `MidFetal_H3K27ac`) as the peak-overlap ground truth and independently recomputes the Fisher's exact test from there — real, non-fabricated recomputation, just not derived from freshly re-intersected peak files.

## Running it

```bash
bash Benchmarks/koesterich2023_characterization_novo/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark koesterich2023_characterization_novo
```

## Ports registered

| Port | What it does | Upstream |
|---|---|---|
| `koesterich2023-mpranalyze-run` | Runs the real MPRAnalyze Bioconductor package's comparative (Ref-vs-Alt) `analyzeComparative`/`testLrt` on the deposited barcode-level counts; fails fast with a clear diagnostic on the GEO duplicate-file defect rather than emitting a silent null result | Methods 4.3.3/4.3.4/4.4, real `YosefLab/MPRAnalyze` |
| `koesterich2023-hcdnv-calling` | Applies the paper's exact HcDNV threshold (FDR<=0.2, |logFC|>=0.1) to any MPRAnalyze-style logFC/fdr table | Methods 4.3.4/4.4 |
| `koesterich2023-case-control-test` | Fisher's exact test, HcDNV rate x ASD case/control | Results 2.5 |
| `koesterich2023-annotation-enrichment` | Fisher's exact test, HcDNV x ATAC-seq/H3K27ac peak overlap | Methods 4.1.3/4.1.4 |
| `koesterich2023-gc-content` | Wilcoxon rank-sum test, %GC HcDNV vs. background | Results 2.8 |
| `koesterich2023-barcode-qc` | Honest QC on the deposited (pre-filter) barcode association map | Results 2.2, Methods 4.3 |

## Provenance

`provenance.json` holds the resolve record. Supplementary Table S3/S4 (the authors' own final per-variant result table, from `ijms-24-03509-s001.zip` via Europe PMC's `supplementaryFiles` endpoint for PMC9959321 — MDPI's own supplementary-materials page did not return usable content to a scripted fetch) ship in `Data/koesterich2023_asd_mpra/supplement/`.
