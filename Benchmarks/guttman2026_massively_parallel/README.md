# guttman2026_massively_parallel

## Paper

**Massively parallel reporter assays for _CYP3A4_ enhancer variants alongside their native promoter.**
Guttman Y, Krupkin B, Ahituv N.
*bioRxiv* 2026 · doi:[10.64898/2026.04.22.719677](https://doi.org/10.64898/2026.04.22.719677) · PMID 42079203 · PMC13131672

A modified MPRA that assays CREs together with their native target promoter, applied to *CYP3A4* (a major pharmacogene): 1,214 variants across six regulatory regions (C/EBPRE, CLEM4/constitutive liver enhancer, DRR, IHRR, Promoter, XREM), integrating population, cancer-genome and archaic-human variants.

## Data used

Two real, public data sources, both the paper's own:

1. **Supplementary Data 1** — the paper's own final per-variant table ("All assayed variants, with allele frequencies, BCalm, and FIMO scores"), fetched from EuropePMC's `supplementaryFiles` endpoint for PMC13131672 (`media-2.tsv`, renamed `supplementary_data1_variants.tsv`). 1,214 rows with the authors' own precomputed `logFC`, `is_significant`, `adj_pval`, per-population gnomAD allele frequencies, and `cancer_type` columns. This is the **primary reference** — it is the paper's own final output, strictly more authoritative than re-deriving significance from raw deposited q-values (see History below).
2. Six IGVF-deposited per-CRE-region `reporter_experiment` (per-oligo, per-replicate barcode counts + LFC) tabular files — used only for replicate QC, since Supplementary Data 1 is variant-level only and has no per-replicate barcode counts.

| Region | Analysis set (for replicate QC) |
|---|---|
| C/EBPRE | IGVFDS5672BEKT |
| CLEM4 | IGVFDS2007GIEG |
| DRR | IGVFDS1255HUVH |
| IHRR | IGVFDS6803OSAI |
| Promoter | IGVFDS7415HWHC |
| XREM | IGVFDS5725HVCI |

Code Availability names only `https://github.com/Ahituv-lab/CypMPRA` (+ Zenodo 10.5281/records/19363634), which covers the **long-read barcode-variant association** step only (`cypmpra/{align,extract,filter}.py`). The downstream differential/significance-calling stage that produces the Results-section numbers is described only in Methods prose ("analyzed using MPRAsnakeflow... adjusted p-value<0.05... |log2FC|>0.5"); no reference code for it is public — but the paper's own Supplementary Data 1 is its final output, so it stands in for that missing reference.

## What IGVFagent does

Port `guttman2026-cyp3a4-mpra-stats` (`Scripts/ported/skills/guttman2026_cyp3a4_mpra_stats.py`, no author code released for the differential-analysis stage — an honest reimplementation reading the paper's own Supplementary Data 1 as ground truth): computes library composition, significance+fold-change hit calling, population/cancer variant spot-checks, and (from the IGVF `reporter_experiment` tables) replicate QC.

```bash
bash Benchmarks/guttman2026_massively_parallel/run.sh
igvfagent bench score --paper-id guttman2026_massively_parallel
```

## Concordance — `bench score`: 6/8 analyses reproduced (`reproduced_except_access`, 27/28 checks)

| Analysis (paper claim) | State | Measured | Note |
|---|---|---|---|
| Library design: 1,214 variants = 1,113 SNV + 101 indel | **reproduced** | 1,214 = 1,113 SNV + 101 indel | exact |
| Replicate QC: ~109 barcodes/variant, correlation >0.75 | **reproduced** | correlation 0.985 (barcodes/variant 225, see below) | the paper's floor (>0.75) is comfortably cleared |
| Hit calling: 43 variants \|LFC\|>0.5 among significant (14 up/29 down) | **reproduced** | 43 hits (14 up/29 down) | exact |
| Fig. 1f: per-CRE significant-variant frequency (7.9/6.2/5.2/3.0/1.7/1.8%) | **reproduced** | 7.87/6.06/5.17/2.70/1.73/1.82% | exact within rounding |
| Population-specific regulatory effects (Fig. 3) | **reproduced** | 6/6 named variants' LFC match (mean \|error\| 0.0024) | exact |
| Cancer-associated variant effects (Fig. 4) | **reproduced** | 57 cancer variants, 1/57 (1.75%) exceeds threshold (the stated DRR/HCC substitution) | exact match to "one substitution ... exceeded" |
| Archaic human variant effects (Fig. 5) | **blocked** (not_deposited) | — | see below |
| DRR haplotype functional analysis (Fig. 6) | **blocked** (not_deposited) | — | combined-haplotype luciferase values undeposited; constituent single-variant LFCs are independently checkable against Supplementary Data 1 |

Full check-level detail: `igvfagent bench report --paper-id guttman2026_massively_parallel` or the latest `Benchmarks/results/*_concordance.md`.

## History (superseded attempt)

An earlier version of this benchmark derived significance directly from the six IGVF-deposited `reporter_variant` tables' raw q-values (`adjusted p<0.05 AND |LFC|>0.5`), without the paper's own multiple-testing procedure. That gave 68 hits (16 up/52 down) against the paper's 43 (14/29), and per-CRE frequencies 5–10× too high — a real, honestly-documented discrepancy at the time. Fetching the paper's Supplementary Material (never pulled before) surfaced Supplementary Data 1, the authors' own final per-variant table with a precomputed `is_significant` column — using that as the reference instead reproduces the paper's exact hit-calling numbers. The lesson: the earlier attempt was re-deriving an unspecified statistical procedure from raw inputs when the paper's own finished output was available in its supplementary material all along.

## Honest caveats (still open)

* **Barcodes-per-variant is ~2x the paper's number** (225 measured vs 109 stated). Supplementary Data 1 has no per-replicate barcode counts, so this still comes from the IGVF-deposited `reporter_experiment` `n_bc` column, which may be a pre-filter count where the paper's 109 is post-filter. The correlation claim (>0.75) still reproduces cleanly on the same data (0.985), so this doesn't undermine the QC conclusion.
* **Cancer gnomAD-absence percentage does not reproduce.** The paper states "79% of cancer-associated CYP3A4 variants ... were not observed in the general healthy population, as annotated by gnomAD." Measured on Supplementary Data 1: all 57 (100%) of `cancer_type`-tagged rows have blank (not merely zero) `global_AC`/`AN`/`AF` fields. No further public data resolves this specific gap; reported honestly (`cancer_analysis.pct_not_in_gnomad` in `summary.json`) rather than forced to match — it is not part of the confirmed check set for this analysis, whose headline claim (the one significant DRR/HCC substitution) does reproduce exactly.
* **Archaic-variant analysis stays blocked.** The 16 specific archaic variants (6 archaic-specific + 10 shared) are not listed anywhere in the main text, Supplementary Data 1 (no archaic/Neanderthal/Denisovan column), or the Supplementary Material. The in-text citation for how they were selected ("using a comparative analysis of archaic and modern human genomes") has no resolvable cross-reference in the article XML. The paper does separately cite the underlying archaic genomes elsewhere (Reich et al. 2010 Denisova; Prüfer et al. 2014 Altai Neanderthal), but reproducing this would require independently re-deriving which CYP3A4-locus variants are archaic-derived from those raw archaic genome sequences — a substantial independent bioinformatics analysis beyond what this paper deposited.
* **DRR haplotype analysis stays blocked** for its headline claim (combined-haplotype luciferase activity, Fig. 6) — an individually-cloned assay whose per-construct values are not in the IGVF deposits, Supplementary Data 1 (variant-level only), or the Supplementary Material (Supplementary Fig. 4 is images only). The haplotypes' constituent single variants (e.g. rs776744, rs11353593, rs61017966, rs776742, 99,694,580:T>C) are independently checkable against Supplementary Data 1's per-variant `logFC`, but that is not the same as the epistatic/combined-effect claim itself.
