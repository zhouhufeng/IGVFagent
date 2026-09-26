# southard2024_comprehensive_transcription

> This is the corrected, real reproduction for the paper the user originally
> asked for as "Joung 2025 TF Perturb-seq fibroblasts" via
> `https://www.nature.com/articles/s41588-025-02283-2`. That DOI is a 2-page
> Nature Genetics Research Briefing by Thomas M. Norman *summarizing* the
> paper below — it has no author named Joung. See
> [`../joung2025_tf_perturbseq/README.md`](../joung2025_tf_perturbseq/README.md)
> for the full correction (that directory previously carried a fabricated
> citation and fabricated accessions for this paper_id).

## Paper

**Comprehensive transcription factor perturbations recapitulate fibroblast transcriptional states**
Southard Kaden M.; Ardy Rico C.; Tang Anran; O'Sullivan Deirdre D.; Metzner Eli; Guruvayurappan Karthik; Norman Thomas M.
*Nature Genetics* **57**: 2323–2334 (2025) · doi:[10.1038/s41588-025-02284-1](https://doi.org/10.1038/s41588-025-02284-1)
(closed access; preprinted as bioRxiv [10.1101/2024.07.31.606073](https://doi.org/10.1101/2024.07.31.606073), open access, PMC11312553 — harvested from the preprint's full text since the published version isn't open access)

Resolver confidence: **1.00** (resolved).

## Real data and code (verified — see "How this was verified" below)

| Resource | Identifier | Status |
|---|---|---|
| Authors' own analysis code | [norman-lab-msk/TFs_CRISPRa](https://github.com/norman-lab-msk/TFs_CRISPRa) | public, pinned at commit `3637f77` under `Data/PaperCode/norman-lab-msk__TFs_CRISPRa/` by `igvfagent paper-code fetch` |
| Raw sequencing reads (SRA) | [PRJNA1108254](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1108254) | public |
| Processed Hs27 fibroblast dataset | Zenodo [10.5281/zenodo.15200179](https://doi.org/10.5281/zenodo.15200179) | public (4 files, 58MB–9.7GB) |
| Processed RPE-1 dataset | Zenodo [10.5281/zenodo.15213619](https://doi.org/10.5281/zenodo.15213619) | public (4 files, up to 29.7GB) |
| Cellranger raw outputs (Hs27, RPE-1) | Zenodo [15213597](https://doi.org/10.5281/zenodo.15213597), [15211972](https://doi.org/10.5281/zenodo.15211972) | public (21.6GB, 46.8GB) |
| Bulk RNA/ATAC/CUT&RUN characterization | Zenodo [10.5281/zenodo.15215216](https://doi.org/10.5281/zenodo.15215216) | public |
| IGVF-format resubmission (in progress) | [norman-lab-msk/igvf-perturbseq](https://github.com/norman-lab-msk/igvf-perturbseq) | code public; no live `IGVFDS…` accession found via Portal search as of this run |

The in-text accessions IGVFagent's harvester found (`GSE186458`, `MiraldiLab/maxATAC_data`, `lmcinnes/enstop`, `moshi4/pyCirclize`, `zenodo.7598955`) are all **third-party tools and reference datasets the paper's methods cite**, not the paper's own data — this preprint's body has no explicit Data/Code Availability section (unusual, but genuine: verified against every section heading in the harvested full text). The real accessions above were found by reading the authors' own GitHub README and the Zenodo API directly, not by the automated harvester.

## What IGVFagent does

**Route:** `perturb_catalog` — Perturbation Catalogue census (CRISPR screen / Perturb-seq), plus a real data-verification step (below)
**Skill output dir:** `Docs/Perturbation/`

```bash
bash Benchmarks/southard2024_comprehensive_transcription/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark southard2024_comprehensive_transcription
```

## How this was verified — real data, not text-only

The paper's own analysis is ~50 Jupyter notebooks operating on raw 10x
cellranger output (up to 46.8GB per cell line) — not realistically
re-runnable end-to-end in an interactive session. Instead,
`verify_guide_library.py` downloads the smallest real processed artefact
(the 58MB guide×droplet UMI matrix from Zenodo 10.5281/zenodo.15200179) and
independently counts its guide-identity axis:

| Quantity | Paper claim (prose) | Measured directly from the authors' deposited matrix | Match |
|---|---:|---:|:---:|
| Guides in final library | 10,979 | **10,979** | ✓ exact |
| Non-targeting control guides | 78 | **78** | ✓ exact |
| Distinct TFs targeted | 1,836 | **1,836** | ✓ exact |

This is a genuine measurement against primary data (not a repeat of the
paper's own text), and it is what `expected.json`'s two `confirmed: true`
checks are based on.

## Concordance

**3 / 3 confirmed checks pass.** 10 further paper-claimed numbers (perturbation counts, cluster counts, AUC values) remain `[UNCONFIRMED]` — verifying those would require running the authors' notebooks against multi-GB/TB-scale raw or cellranger data, which is outside what this session could execute. They are reported, not scored, per this repo's convention: an unconfirmed check never counts as a pass. See `Docs/Benchmark/*_southard2024_comprehensive_transcription/replication_report.md` for the full table.

## Honest caveats

* **Run directory is not paper-tagged.** This route's catalogue-census CLI steps accept no `--label`, so `concordance.py` matches the skill's default `summary` directory. The two real checks above are pinned to `guide_library_verification.json` (via `extra_search_dirs`), which is paper-specific regardless.
* **10 unconfirmed checks remain.** Promoting them would need the full raw-data pipeline (cellranger + the authors' regression/clustering notebooks), which is a multi-hour/GB-scale undertaking, not something this benchmark attempts to fake.
* **No IGVF Portal accession found.** The Norman lab's `igvf-perturbseq` repo says this data was "formatted for IGVF submission," but Portal search found no matching `MeasurementSet` as of this run — likely still in progress upstream.

## Provenance

`provenance.json` in this directory holds the full resolve / harvest / route record, including every source URL consulted. The repo fetch record (pinned commit, fetch timestamp) is in `Data/PaperCode/norman-lab-msk__TFs_CRISPRa/` (created by `igvfagent paper-code fetch`, outside version control).
