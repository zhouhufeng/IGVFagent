# ENCODE Smoke Analysis

Generated: 2026-05-06 13:16:38 EDT

## Purpose

These public ENCODE checks exercise major input data classes for IGVF interpretation work: accessibility, DNA binding, transcription, signal tracks, and interval/peak files.

## encode_atac_seq_experiments

Chromatin accessibility experiments.

Parameters: `{'type': 'Experiment', 'assay_title': 'ATAC-seq', 'limit': '3'}`
HTTP status: 200
Returned rows: 3
Reported total: 560
Saved response: `<repo>/Data/20260506_131633_encode_atac_seq_experiments.json`

Top assays:
- ATAC-seq: 6

Top biosamples:
- regulatory T cell: 1
- Mus musculus strain C57BL/6NJ regulatory T cell: 1
- T-helper 17 cell: 1
- Mus musculus strain C57BL/6NJ T-helper 17 cell: 1
- T-helper 2 cell: 1
- Mus musculus strain C57BL/6NJ T-helper 2 cell: 1

Top file formats:
- none observed

## encode_chip_seq_experiments

Transcription factor or histone ChIP-seq experiments.

Parameters: `{'type': 'Experiment', 'assay_slims': 'DNA binding', 'limit': '3'}`
HTTP status: 200
Returned rows: 3
Reported total: 13733
Saved response: `<repo>/Data/20260506_131634_encode_chip_seq_experiments.json`

Top assays:
- ChIP-seq: 2
- capture Hi-C: 2
- TF ChIP-seq: 1
- Control ChIP-seq: 1

Top biosamples:
- HepG2: 2
- Homo sapiens HepG2 expressing RNAi targeting H. sapiens FOXA3: 1
- Homo sapiens HepG2: 1
- K562: 1
- Homo sapiens K562: 1

Top file formats:
- none observed

## encode_rna_seq_experiments

Transcriptome profiling experiments.

Parameters: `{'type': 'Experiment', 'assay_slims': 'Transcription', 'limit': '3'}`
HTTP status: 200
Returned rows: 3
Reported total: 6696
Saved response: `<repo>/Data/20260506_131634_encode_rna_seq_experiments.json`

Top assays:
- PRO-cap: 6

Top biosamples:
- HEK293T: 2
- Homo sapiens HEK293T genetically modified (insertion) using CRISPR targeting H. sapiens ZNF143 treated with 500 nM dTAGV-1 for 6 hours: 1
- Homo sapiens HEK293T genetically modified (insertion) using CRISPR targeting H. sapiens ZNF143: 1
- SUDHL5: 1
- Homo sapiens SUDHL5: 1

Top file formats:
- none observed

## encode_dna_accessibility_experiments

DNA accessibility experiments useful for regulatory variant interpretation.

Parameters: `{'type': 'Experiment', 'assay_slims': 'DNA accessibility', 'limit': '3'}`
HTTP status: 200
Returned rows: 3
Reported total: 4597
Saved response: `<repo>/Data/20260506_131635_encode_dna_accessibility_experiments.json`

Top assays:
- DNase-seq: 6

Top biosamples:
- activated naive CD4-positive, alpha-beta T cell: 1
- Homo sapiens activated naive CD4-positive, alpha-beta T cell male adult (43 years) treated with 10 ng/mL Interleukin-2 for 5 days, anti-CD3 and anti-CD28 coated beads for 7 days: 1
- kidney tubule cell: 1
- Homo sapiens kidney tubule cell male adult: 1
- substantia nigra: 1
- Homo sapiens substantia nigra tissue male adult: 1

Top file formats:
- none observed

## encode_bigwig_files

Signal tracks commonly used for regulatory annotation.

Parameters: `{'type': 'File', 'file_format': 'bigWig', 'limit': '3'}`
HTTP status: 200
Returned rows: 3
Reported total: 274659
Saved response: `<repo>/Data/20260506_131636_encode_bigwig_files.json`

Top assays:
- PRO-cap: 3

Top biosamples:
- HEK293T: 3

Top file formats:
- bigWig: 3

## encode_bed_files

Region/peak calls and genomic interval files.

Parameters: `{'type': 'File', 'file_format': 'bed', 'limit': '3'}`
HTTP status: 200
Returned rows: 3
Reported total: 733220
Saved response: `<repo>/Data/20260506_131638_encode_bed_files.json`

Top assays:
- PRO-cap: 3

Top biosamples:
- HEK293T: 3

Top file formats:
- bed: 3
- bed3+: 3
