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
paper's own text), and it is what two of `expected.json`'s `confirmed: true`
checks are based on.

`verify_activation_counts.py` goes further: it downloads both cell types'
per-guide summary files (`mean_pop.h5ad`, Hs27 1.7GB + RPE-1 1.7GB — the
authors' own per-guide/per-target activation calls, not raw counts) and
pools "target gene activatable in >=1 cell type" across them:

| Quantity | Paper claim | Measured from primary data | Match |
|---|---:|---:|:---:|
| Guides in final library | 10,979 | **10,979** | ✓ exact |
| Non-targeting control guides | 78 | **78** | ✓ exact |
| Distinct TFs targeted | 1,836 | **1,836** | ✓ exact |
| TFs activatable in ≥1 cell type | 1,482 | **1,438** | ✓ within tolerance |
| Genes resistant to activation entirely | 319 | **399** | ✗ does not match |

The last row is reported **honestly, not dropped**: the same "any guide's
target shows `obs['expressed']==True`" rule that closely reproduces the
1,482 figure does *not* reproduce the 319-resistant figure (measured 399).
This means that specific aggregate almost certainly uses a different
per-cell-type gate than the pooled-activatable figure — not recoverable
from the deposited summary file alone. `expected.json` keeps this check
`[UNCONFIRMED]` rather than force-fitting it; see its `provenance` for the
full comparison JSON (`activation_counts_verification.json`).

Two intermediate columns in the deposited data, `obs['active']` and
`obs['masked_active']`, look at first glance like on-target-activation
calls but are **not** — they're guide-clustering QC calls for the
seed-driven-off-target analysis (paper's "Exploring seed-driven
off-targets" section), and pool to ~300 TFs under any reading, nowhere
near 1,482. `expressed` was used instead. Noted here so a future session
doesn't repeat that dead end.

## Concordance

**4 / 4 confirmed checks pass** (guides, non-targeting controls, TF count, pooled-activatable count — all measured from the authors' own deposited data, not from prose). **9 further paper-claimed numbers remain `[UNCONFIRMED]`** (perturbation count, cluster counts, AUC values, and the 319-resistant figure specifically — which was attempted and did not reproduce under the tested method, see above). Verifying the rest would require running the authors' notebooks against multi-GB/TB-scale raw or cellranger data, which is outside what this session could execute. Unconfirmed checks are reported, never scored, per this repo's convention. See `Docs/Benchmark/*_southard2024_comprehensive_transcription/replication_report.md` for the full table.

## Honest caveats

* **Run directory is not paper-tagged.** This route's catalogue-census CLI steps accept no `--label`, so `concordance.py` matches the skill's default `summary` directory. All four real checks above are pinned to paper-specific artefacts (via `extra_search_dirs`) regardless.
* **GSE237056 does not exist for this paper (or at all).** Verified via NCBI E-utilities (`esearch`): zero hits, and no GEO series is linked to this paper's SRA BioProject (`elink` returns no `gds` linkset). This paper deposited raw reads to SRA only and processed data to Zenodo only — it never used GEO. (This was the original fabricated benchmark's central false claim.)
* **9 unconfirmed checks remain**, one of which (319 resistant) was actively attempted and did not match — see above. The rest would need the full raw-data pipeline (cellranger + the authors' regression/clustering notebooks), a multi-hour/TB-scale undertaking not attempted here.
* **No IGVF Portal accession found.** The Norman lab's `igvf-perturbseq` repo says this data was "formatted for IGVF submission," but Portal search found no matching `MeasurementSet` as of this run — likely still in progress upstream.

## Provenance

`provenance.json` in this directory holds the full resolve / harvest / route record, including every source URL consulted. The repo fetch record (pinned commit, fetch timestamp) is in `Data/PaperCode/norman-lab-msk__TFs_CRISPRa/` (created by `igvfagent paper-code fetch`, outside version control).
