# [CORRECTED] TF Perturb-seq in primary fibroblasts

> **This benchmark's original citation was fabricated.** It attributed this
> paper to "Joung J, ..., Zhang F., *Nature Genetics* 57:828–838 (2025)" and
> cited data accessions `GSE237056` (embargoed) and `SCP2169`. Neither the
> citation nor the accessions correspond to any real paper matching this
> DOI. The real facts, reverified below, are:
>
> * DOI `10.1038/s41588-025-02283-2` (this benchmark's original identifier,
>   and the URL the user supplied) is a 2-page Nature Genetics **Research
>   Briefing** — a first-person summary written by the senior author,
>   Thomas M. Norman — of the real research article below. It has no author
>   named Joung.
> * The invented title "A transcription factor atlas of directed
>   differentiation" and page range 828–838 belong to an unrelated 2023
>   *Cell* paper (Joung, J. et al., cited only as reference #4 *inside*
>   Norman's briefing) — likely the source of the "Joung 2025" mix-up.
> * `GSE237056` has no resolvable GEO record. `SCP2169` redirects to an
>   unrelated human-tonsil snRNA-seq study, not this paper.
>
> **The real, full reproduction — including a genuine data-backed
> confirmation of the paper's guide-library counts — lives at
> [`Benchmarks/southard2024_comprehensive_transcription/`](../../southard2024_comprehensive_transcription/README.md).**
> This directory is kept (rather than deleted) as a record of the
> correction and to redirect anyone who arrives here via the old paper_id.

[![paper](https://img.shields.io/badge/Nat%20Genet-57:2323--2334-blue)](https://doi.org/10.1038/s41588-025-02284-1)
[![preprint](https://img.shields.io/badge/bioRxiv-10.1101%2F2024.07.31.606073-lightgrey)](https://doi.org/10.1101/2024.07.31.606073)
[![code](https://img.shields.io/badge/code-norman--lab--msk%2FTFs__CRISPRa-green)](https://github.com/norman-lab-msk/TFs_CRISPRa)
[![data](https://img.shields.io/badge/Zenodo-normanlabmsk%20community-orange)](https://zenodo.org/communities/normanlabmsk/)

## Corrected citation

Southard KM, Ardy RC, Tang A, O'Sullivan DD, Metzner E, Guruvayurappan K,
Norman TM. **Comprehensive transcription factor perturbations recapitulate
fibroblast transcriptional states.** *Nature Genetics* **57**: 2323–2334
(2025). DOI: [10.1038/s41588-025-02284-1](https://doi.org/10.1038/s41588-025-02284-1)
(closed access; preprinted as bioRxiv
[10.1101/2024.07.31.606073](https://doi.org/10.1101/2024.07.31.606073),
open access, PMC11312553).

## Real data sources (verified, not fabricated)

| Resource | Identifier | Status |
|---|---|---|
| Raw sequencing reads (SRA) | [PRJNA1108254](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1108254) | public |
| Processed Hs27 fibroblast dataset (Zenodo) | [10.5281/zenodo.15200179](https://doi.org/10.5281/zenodo.15200179) | public — includes the 58MB guide×droplet UMI matrix this benchmark's sibling directory verifies against |
| Processed RPE-1 dataset (Zenodo) | [10.5281/zenodo.15213619](https://doi.org/10.5281/zenodo.15213619) | public |
| Authors' own analysis code | [norman-lab-msk/TFs_CRISPRa](https://github.com/norman-lab-msk/TFs_CRISPRa) | public, pinned commit under `Data/PaperCode/` by `igvfagent paper-code fetch` |
| IGVF-format submission (in progress) | [norman-lab-msk/igvf-perturbseq](https://github.com/norman-lab-msk/igvf-perturbseq) | code public; no live IGVFDS accession found via Portal search as of this benchmark's last run |

See [`southard2024_comprehensive_transcription/expected.json`](../../southard2024_comprehensive_transcription/expected.json)
for the full set of paper-claimed numbers this reproduction checks.

## What was actually confirmed (real, not text-only)

`Benchmarks/southard2024_comprehensive_transcription/verify_guide_library.py`
downloads the authors' own deposited guide×droplet UMI matrix (Zenodo
10.5281/zenodo.15200179) and counts its guide-identity axis directly —
independent of trusting the paper's prose:

| Quantity | Paper claim | Measured from primary data |
|---|---:|---:|
| Guides in final library | 10,979 | **10,979** (exact match) |
| Non-targeting control guides | 78 | **78** (exact match) |
| Distinct TFs targeted | 1,836 | **1,836** (exact match) |

## This directory's remaining role

This directory is not re-scaffolded to duplicate that work. Its `run.sh` /
`expected.json` still perform the harmless, generic Perturbation Catalogue
census check that was never fabricated (only the paper metadata was) — see
`OPERATIONS.md`. For the paper-specific reproduction, use
`southard2024_comprehensive_transcription/`.
