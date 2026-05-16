# Skill: SHARE-seq Joint ATAC + RNA QC

Use this skill for end-to-end QC and joint-cell analytics on SHARE-seq
deposits (Ma et al. *Cell* 2020 — simultaneous scATAC + scRNA per cell,
shared round-1/2/3 split-pool barcodes). The math is a clean-room
reimplementation of the QC algorithms in
[broadinstitute/epi-SHARE-seq-pipeline](https://github.com/broadinstitute/epi-SHARE-seq-pipeline)
(MIT, 2021).

This skill is designed to consume the **processed deposits** that IGVF
Portal SHARE-seq AnalysisSets publish — fragments BED + h5ad — and does
not re-run FASTQ → BAM alignment. Use the upstream Broad WDL pipeline
(or the IGVF Portal ingestion service) for raw-FASTQ processing.

## Commands

```bash
# 1. Discover IGVF Portal SHARE-seq AnalysisSets / MeasurementSets
python3 Scripts/share_seq_skill.py pull-portal --limit 50 --label survey

# 2. (Optional) Demultiplex a raw R2 FASTQ with the SHARE-seq round-1/2/3
#    24-mer barcode + 1-Hamming-mismatch + ±1 bp shift correction
python3 Scripts/share_seq_skill.py demultiplex-bcs \
    --fastq sample_R2.fastq.gz --whitelist whitelist_24mer.txt \
    --out sample_R2.tagged.fastq.gz --label demo

# 3. Per-barcode ATAC QC from a fragments BED + (optional) TSS BED + peaks BED
python3 Scripts/share_seq_skill.py fragment-qc \
    --fragments fragments.bed.gz --tss-bed tss.bed --peaks-bed peaks.bed \
    --label demo

# 4. Per-barcode RNA QC from an h5ad sparse gene-count matrix
python3 Scripts/share_seq_skill.py rna-qc \
    --h5ad rna.h5ad --label demo

# 5. Joint per-barcode pass calls (Ma 2020 thresholds)
python3 Scripts/share_seq_skill.py joint-qc \
    --atac-qc Docs/SHAREseq/<ts>_demo_atac_qc.tsv \
    --rna-qc Docs/SHAREseq/<ts>_demo_rna_qc.tsv --label demo

# 6. Multiplet detection by pairwise Jaccard over fragment-coord sets
python3 Scripts/share_seq_skill.py multiplet-detect \
    --fragments fragments.bed.gz --label demo
```

## Input file conventions

| Input | Expected format |
|---|---|
| `fragments.bed[.gz]` | 4-col + count BED (chrom, start, end, barcode, count). Matches the SHARE-seq pipeline `bam_to_fragments.py` schema (`start = read.reference_start + 4; end = start + tlen - 4`). |
| `rna.h5ad` | AnnData with barcode × gene; mitochondrial genes detected by `MT-` / `mt-` prefix. |
| `tss.bed` | BED of TSS positions (any annotation; chrom + start). |
| `peaks.bed` | BED of MACS / iterative-LSI peaks. |

## Default cell-calling thresholds (Ma 2020)

| Modality | Metric | Pass |
|---|---|---|
| RNA | UMIs ≥ `min_umis` (default 100) AND genes ≥ `min_genes` (default 200) | both required |
| ATAC | fragments ≥ `min_frags` (default 100) AND TSS enrichment ≥ `min_tss` (default 4) | both required |
| Joint | RNA_pass AND ATAC_pass | "both" |

## Citation

- **Ma S et al. (2020)** *Cell* 183:1103–1116. "Chromatin potential identified
  by shared single-cell profiling of RNA and chromatin." doi:10.1016/j.cell.2020.09.056
- Pipeline reference: [broadinstitute/epi-SHARE-seq-pipeline](https://github.com/broadinstitute/epi-SHARE-seq-pipeline) (MIT)

## License-clean Python stack

`pandas` (BSD-3), `numpy` (BSD-3), `scipy` (BSD-3), `anndata` (BSD-3),
`statsmodels` (BSD-3). All non-GPL. Only the FASTQ demultiplex is a hand-
rolled scan over gz files using the standard library — no pysam, no
GPL bowtie2 chain at runtime.
