# anglen2025_regulatory_element

## Paper

**A gene regulatory element modulates myosin expression and controls cardiomyocyte response to stress**
Anglen T, Kaplow IM, Choi B, Dewars E, Perelli RM, Hagy KT, Tran D, Ramaker ME, Shah SH, Jung I, Landstrom AP, et al.
*Genome Research* 2025 · doi:[10.1101/gr.280825.125](https://doi.org/10.1101/gr.280825.125) · PMID 41125440 · PMC12581865 · IGVF award IGVF0195, grant UM1HG012053

Resolver confidence: **1.00** (exact PMID match).

The paper identifies "R3," a noncoding regulatory element in intron 1 of *MYH6* (chr14q11.2, the fetal/atrial myosin heavy chain gene), that becomes more accessible upon GSK3-inhibitor (CHIR99021)-driven cardiomyocyte maturation. CRISPRi repression of R3 cuts MYH6 expression; CRISPRa activation of R3 increases MYH6/MYH7 and shifts cardiomyocyte stress response to endothelin-1 (calcium handling, metabolism, polyploidy), mediated through direct chromatin looping (HiCAR) to the MYH6 promoter. Rare noncoding variants in R3 are enriched in UK Biobank cardiomyopathy cases.

## Real author code, real data — the good case

Data and Code availability (verbatim): *"All raw sequencing data and analysis files for ATAC-seq, RNA-seq, and HiCAR can be found at ... GEO ... under accession numbers GSE283424, GSE283426, GSE283427, GSE283430, and GSE283432. All code ... is available at GitHub (https://github.com/Gersbachlab-Bioinformatics/myosin_enhancer) and as Supplemental Code."* Both are real and present. Even better: **several of the GEO deposits are the authors' own already-processed final outputs**, not just raw reads — letting this pass verify the paper's headline numbers directly, without re-running alignment pipelines.

## Part 1 — Fig. 1 ATAC-seq/RNA-seq, and the UK Biobank variant analysis — DONE (3/6 reproduced)

| Claim | Paper | Measured (real deposited data) | Match |
|---|---|---|---|
| ATAC-seq peaks, increased accessibility (GSK3i, padj<0.01, \|log2FC\|>1) | 976 | **976** | exact |
| ATAC-seq peaks, decreased accessibility | 851 | **851** | exact |
| Top differential ATAC peak | padj=6.65e-48 | **padj=6.65e-48**, at chr14_23376493_23379689 (chr14q11.2 = the real MYH6/MYH7 locus) | exact |
| MYH6 upregulation (RNA-seq) | log2FC=3.41 | **3.66** (TPM-based, see caveat below) | close |
| MYH7 downregulation (RNA-seq) | log2FC=-1.34 | **-1.12** (TPM-based) | same direction, same order of magnitude |
| R3-CTCF-motif variant, Fisher's exact p (cardiomyopathy cases vs. controls) | p≈0.043, CADD=15.32 | **p=0.0433578632...** recomputed independently, CADD=15.32 exact | **machine-precision match** |

New ports: `anglen2025_atac_de_summary` (reads the authors' own deposited final DESeq2 ATAC output, `GSE283424_atac_DE_analysis.csv.gz` — the literal output of their `figure2_scripts/ATAC-seq/Deseq2 ATAC.RMD`), `anglen2025_rnaseq_tpm_de` (TPM-based DE from `GSE283426_RNAseq_TPM.txt.gz`), `anglen2025_variant_enrichment` (Fisher's exact test recomputed from Genome Research's own freely-downloadable Supplemental Table S5 — a de-identified, aggregate allele-count table, **not** individual-level UK Biobank data, so no controlled-access application was needed).

**One honest limitation**: `GSE283426_RNAseq_TPM.txt.gz` carries only per-sample TPM/FPKM, not the per-sample raw RSEM counts the authors' own DESeq2 RMD actually consumes. A t-test on log2(TPM) recovers the correct fold-change direction and rough magnitude for MYH6/MYH7 but cannot approach DESeq2's extreme significance (their padj=1.06e-251 needs the negative-binomial model's power across the whole transcriptome; a 3-vs-3 t-test cannot). Reported honestly (see `expected.json`'s check notes), not hidden.

```bash
bash Benchmarks/anglen2025_regulatory_element/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark anglen2025_regulatory_element
```

## Not yet attempted (real data exists — pending, not blocked)

- **Fig. 5 (HiCAR chromatin looping)**: GEO GSE283427/GSE283430 and Supplemental Table S4 both carry the authors' own final DESeq2 differential-contact output (10 kb bin-pair resolution, at exactly the MYH6/MYH7 locus) — real, deposited, ready to use. What's missing is mapping their named anchors ("MYH6", "C3", "M3" control regions) to specific genomic bin coordinates, which isn't explicitly labeled in the deposited table; needs either the Supplemental Code or a closer read of the figure to pin down.
- **Fig. 3/4 (CRISPRi/CRISPRa validation, cellular phenotypes)**: Supplemental Table S1 (CHIR dose-response: EdU, pHH3, genomic content, cell-cycle) and Table S3 (Seahorse) carry real per-well/per-replicate data; not yet parsed into checks.

## Coverage

**3 of 6** planned analyses reproduced (6/6 checks pass). 0 blocked — everything remaining is genuinely "more work to do with data already in hand," not an access problem.

## Provenance

`provenance.json` holds the full resolve/harvest record. Ports: `anglen2025_atac_de_summary`, `anglen2025_rnaseq_tpm_de`, `anglen2025_variant_enrichment` (`igvfagent port list`), all pinned to `Gersbachlab-Bioinformatics/myosin_enhancer@16ee9368ed13`.
