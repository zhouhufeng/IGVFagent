# OPERATIONS — Martyn 2025 "Rewriting regulatory DNA..." / Variant-EFFECTS

> **Corrected 2026-09-25.** This file previously described a different paper's
> loci (GATA1 +24/+58 kb, MYC chr8:128.74M — those are Yao et al. 2024's
> ENCODE4 CRISPRi loci, see `Benchmarks/yao2024_encode4_crispri/`) and assumed
> this paper's IGVF data was raw guide×bin counts needing the `flowfish`
> MLE pipeline. Neither was true. This paper's own IGVF-hosted files are
> **already-scored per-variant effect tables** (effect_size, p-value, BH-FDR
> columns present) — IGVF's uniform pipeline, not raw reads. See
> `README.md`'s "Corrected" note for the full story.

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

## Two stages

| Stage | Mode | What it tests |
|---|---|---|
| 1. Portal discovery + mechanics demo | always runs | `flowfish pull-portal` (live enumeration) + a synthetic 20-element run through `flowfish simulate/estimate-effects/real-space/score-elements`, which demonstrates the analytical-chain shape but is **not** run on this paper's data |
| 2. Real paper data | always runs, online-only | Downloads the paper's own 3 headline GRCh38 variant-effects tables by their fixed accessions and scores them with `analyze_real_data.py` |

## Run everything

```bash
bash Benchmarks/martyn2025_variant_flowfish/run.sh
python3 Benchmarks/concordance.py --benchmark martyn2025_variant_flowfish
```

Takes ~30 s total, no local input file needed — unlike a raw-counts pipeline, the real-data stage only downloads 3 small (<10 KB gzipped) already-scored TSVs directly from the public IGVF Portal.

## Where artefacts land

```
Docs/FlowFISH/<ts>_martyn2025_variant_flowfish_portal.tsv        # Stage 1: live portal manifest
Docs/FlowFISH/<ts>_martyn2025_variant_flowfish_pipeline_*.tsv    # Stage 1: synthetic mechanics-demo outputs
Data/Benchmarks/martyn2025_variant_flowfish/real_data/<ts>_martyn2025_variant_flowfish/
├── PPIF_promoter_GRCh38.tsv    # IGVFFI4057VSBO, from IGVFDS5056OAGR
├── PPIF_enhancer_GRCh38.tsv    # IGVFFI4333XLOF, from IGVFDS5031MNRR
├── IL2RA_promoter_GRCh38.tsv   # IGVFFI4854DWEG, from IGVFDS1824XDMU
└── summary.json                # per-element n_variants / fraction_significant / effect_size range
```

## Ground-truth ratios found in the real data

| Element | n variants | Significant (FDR<0.05) | Effect-size range |
|---|---:|---:|---:|
| PPIF promoter | 41 | 38 (92.7%) | [-0.221, +0.148] |
| PPIF enhancer | 98 | 50 (51.0%) | [-0.256, +0.206] |
| IL2RA promoter | 87 | 78 (89.7%) | [-0.599, +6.044] |

Independent cross-check: PPIF enhancer sits 60,784 bp from the PPIF promoter TSS — matches the paper's own stated "~60.5 kb upstream" almost exactly.

## Extending this benchmark (not yet done)

Three more real IGVF accessions were located but not yet downloaded/analyzed:

* **PPIF splice site**: IGVFDS8174BPPS, IGVFDS8267JJHO, IGVFDS7090INDM. Fetch the same way as the 3 elements above — `portal_get` the AnalysisSet, find its `variant effects` / `GRCh38` file, download via `https://api.data.igvf.org/tabular-files/<accession>/@@download/<accession>.tsv.gz`.
* **lentiMPRA parallel arm**: IGVFDS1003XTAF, IGVFDS1376WOXJ (PPIF promoter, episomal reporter design rather than endogenous prime-editing). Comparing this to the endogenous PPIF-promoter numbers above would reproduce the paper's stated finding that "effects [are] specific to genomic context" (episomal vs. endogenous).
* IL2RA's own splice-site/enhancer equivalents, if the paper tested more than the promoter — not yet checked for.

## Running through the UI

```
Reproduce the Martyn 2025 Variant-EFFECTS benchmark (paper: "Rewriting
regulatory DNA to dissect and reprogram gene expression", Cell 2025,
PMID 40245860):

1. Call portal_get on IGVFDS5056OAGR, IGVFDS5031MNRR, and IGVFDS1824XDMU
   (PPIF promoter, PPIF enhancer, IL2RA promoter AnalysisSets).
2. For each, find the file with content_type="variant effects" and
   assembly="GRCh38"; download it.
3. For each file, report: number of variants, number significant at
   FDR<0.05 (fdr_nlog10 > -log10(0.05)), and the effect_size range.
4. Confirm the PPIF enhancer sits ~60.5 kb from the PPIF promoter's
   position, cross-validating against the paper's own stated distance.
```

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `curl` returns an HTML error page instead of a gzip | IGVF Portal accession renamed/withdrawn since this was written | Re-run `portal_get` on the parent AnalysisSet (e.g. `IGVFDS5056OAGR`) to find the current file accession before hardcoding a new one |
| `analyze_real_data.py` errors on a KeyError | Portal file schema changed columns | `zcat <file> \| head -1` to see current columns; the script expects `chr, pos, effect_size, fdr_nlog10` |
| Downloaded file is hg19, not GRCh38 | Grabbed the wrong sibling file — every accession publishes both an hg19 and a GRCh38 tabular file | Check `content_type`/`assembly` fields via `portal_get` before picking the `href` to download |

## License + provenance

* **Paper data**: IGVF Portal AnalysisSets/MeasurementSets (Engreitz lab) — public, CC-BY 4.0.
* **Code**: IGVFagent Apache-2.0. `analyze_real_data.py` is new, written for this benchmark; `Scripts/flowfish_pipeline.py` is an unrelated generic clean-room pipeline (see README's caveats).
* **Citation**: Martyn GE et al. *Cell* **188**(12): 3349–3366.e23 (2025). doi:10.1016/j.cell.2025.03.034 · PMID:40245860 · PMCID:PMC12167154
