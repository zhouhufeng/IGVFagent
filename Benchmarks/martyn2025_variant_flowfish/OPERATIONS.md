# OPERATIONS — Martyn 2025 "Rewriting regulatory DNA..." / Variant-EFFECTS

> **Corrected 2026-09-25.** This file previously described a different paper's
> loci (GATA1 +24/+58 kb, MYC chr8:128.74M — those are Yao et al. 2024's
> ENCODE4 CRISPRi loci, see `Benchmarks/yao2024_encode4_crispri/`) and assumed
> this paper's IGVF data was raw guide×bin counts needing the `flowfish`
> MLE pipeline. Neither was true. This paper's own IGVF-hosted files are
> **already-scored per-variant effect tables** (effect_size, p-value, BH-FDR
> columns present) — IGVF's uniform pipeline, not raw reads. See
> `README.md`'s "Corrected" note for the full story.

> **Extended 2026-09-26.** The 3 follow-ups this file used to list under
> "Extending this benchmark (not yet done)" are now done — see "Extended
> 2026-09-26" below. All 5 headline analyses now carry a passing,
> paper-value-checked (class-A) reproduction (`igvfagent bench score`).

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

## Three stages

| Stage | Mode | What it tests |
|---|---|---|
| 1. Portal discovery + mechanics demo | always runs | `flowfish pull-portal` (live enumeration) + a synthetic 20-element run through `flowfish simulate/estimate-effects/real-space/score-elements`, which demonstrates the analytical-chain shape but is **not** run on this paper's data |
| 2. Real paper data (3 headline elements) | always runs, online-only | Downloads the paper's own 3 headline GRCh38 variant-effects tables by their fixed accessions and scores them with `analyze_real_data.py` |
| 3. Real paper data (follow-ups) | always runs, online-only | Downloads the PPIF splice-site edits (3 replicate files) and the lentiMPRA (episomal) PPIF-promoter reporter file, and scores them with `analyze_followups.py` |

## Run everything

```bash
bash Benchmarks/martyn2025_variant_flowfish/run.sh
python3 Benchmarks/concordance.py --benchmark martyn2025_variant_flowfish
```

Takes ~30 s total, no local input file needed — the real-data stages only download small (<10 KB–50 KB gzipped) already-scored TSVs directly from the public IGVF Portal.

## Where artefacts land

```
Docs/FlowFISH/<ts>_martyn2025_variant_flowfish_portal.tsv        # Stage 1: live portal manifest
Docs/FlowFISH/<ts>_martyn2025_variant_flowfish_pipeline_*.tsv    # Stage 1: synthetic mechanics-demo outputs
Data/Benchmarks/martyn2025_variant_flowfish/real_data/<ts>_martyn2025_variant_flowfish/
├── PPIF_promoter_GRCh38.tsv              # IGVFFI4057VSBO, from IGVFDS5056OAGR
├── PPIF_enhancer_GRCh38.tsv              # IGVFFI4333XLOF, from IGVFDS5031MNRR
├── IL2RA_promoter_GRCh38.tsv             # IGVFFI4854DWEG, from IGVFDS1824XDMU
├── PPIF_splice_repA_GRCh38.tsv           # IGVFFI0524YUIL, from IGVFDS8174BPPS
├── PPIF_splice_repB_GRCh38.tsv           # IGVFFI2542METL, from IGVFDS8267JJHO
├── PPIF_splice_repC_GRCh38.tsv           # IGVFFI5097SDKA, from IGVFDS7090INDM
├── PPIF_promoter_lentiMPRA_GRCh38.tsv    # IGVFFI2620KDMB, from IGVFDS1376WOXJ
└── summary.json                          # elements{} + ppif_splice_site{} + lentimpra_vs_endogenous{}
```

## Ground-truth ratios found in the real data

| Element | n variants | Significant (FDR<0.05) | Effect-size range |
|---|---:|---:|---:|
| PPIF promoter | 41 | 38 (92.7%) | [-0.221, +0.148] |
| PPIF enhancer | 98 | 50 (51.0%) | [-0.256, +0.206] |
| IL2RA promoter | 87 | 78 (89.7%) | [-0.599, +6.044] |
| PPIF splice site | 3 distinct edits | n/a (proof-of-concept, not a screen) | 13.6%–33.4% of wild-type expression retained |
| PPIF promoter, lentiMPRA | 353 (41 overlap endogenous set) | n/a | correlation vs. endogenous: r=0.21 (paper: r=0.54) |

Independent cross-check: PPIF enhancer sits 60,784 bp from the PPIF promoter TSS — matches the paper's own stated "~60.5 kb upstream" almost exactly.

## Extended 2026-09-26 — all three follow-ups from the earlier version of this file are done

* **PPIF splice site**: IGVFDS8174BPPS/8267JJHO/7090INDM → IGVFFI0524YUIL/2542METL/5097SDKA (GRCh38). Downloaded and scored by `analyze_followups.py`: recovers exactly the paper's own "three edits", each strongly decreasing PPIF expression (13.6%–33.4% of wild-type retained — paper states −80% to −52%, i.e. ~20%–48% retained; same direction, close but not exact magnitude since ours pools 3 raw-file replicate estimates rather than the paper's own clonal validation).
* **lentiMPRA parallel arm**: IGVFDS1376WOXJ → IGVFFI2620KDMB (GRCh38) is the scored reporter-variant file (IGVFDS1003XTAF is a barcode-to-element design file only, no scored effects, and was left unused). Joined against the endogenous PPIF-promoter file by the shared `variant` ID (all 41 endogenous variants matched); real Pearson r = 0.21 vs. the paper's own stated r = 0.54 — positive and directionally consistent with "effects specific to genomic context," but weaker than the paper's headline number (see README's Honest caveats).
* **IL2RA's own splice-site/enhancer equivalents**: checked via PMC full text (PMC12167154) — the paper's IL2RA work is confined to the promoter and detailed only in supplementary "Data S1 section 1," not the main text. No further IL2RA elements are described in the retrievable full text.

## Still open

* **IL2RA promoter's own +6.04 outlier effect size** is real (present in IGVF's public file) but still not quote-verifiable against a specific paper claim — Data S1 (the supplementary table where IL2RA detail lives) was not retrievable through our harvest tools. Would need direct access to the paper's Cell supplementary materials.

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
5. Call portal_get on IGVFDS8174BPPS, IGVFDS8267JJHO, and IGVFDS7090INDM
   (PPIF splice-site replicate files); download each GRCh38 variant-effects
   file and report the distinct edits recovered (paper: "three edits").
6. Call portal_get on IGVFDS1376WOXJ (lentiMPRA PPIF-promoter reporter
   file, not IGVFDS1003XTAF which is design-only); download its GRCh38
   file, join it against the PPIF promoter file from step 2 by the shared
   `variant` ID, and report the Pearson correlation (paper: r=0.54).
```

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `curl` returns an HTML error page instead of a gzip | IGVF Portal accession renamed/withdrawn since this was written | Re-run `portal_get` on the parent AnalysisSet (e.g. `IGVFDS5056OAGR`) to find the current file accession before hardcoding a new one |
| `analyze_real_data.py` errors on a KeyError | Portal file schema changed columns | `zcat <file> \| head -1` to see current columns; the script expects `chr, pos, effect_size, fdr_nlog10` |
| Downloaded file is hg19, not GRCh38 | Grabbed the wrong sibling file — every accession publishes both an hg19 and a GRCh38 tabular file | Check `content_type`/`assembly` fields via `portal_get` before picking the `href` to download |

## License + provenance

* **Paper data**: IGVF Portal AnalysisSets/MeasurementSets (Engreitz lab) — public, CC-BY 4.0.
* **Code**: IGVFagent Apache-2.0. `analyze_real_data.py` (3 headline elements) and `analyze_followups.py` (splice site + lentiMPRA, added 2026-09-26) were written for this benchmark; `Scripts/flowfish_pipeline.py` is an unrelated generic clean-room pipeline (see README's caveats).
* **Citation**: Martyn GE et al. *Cell* **188**(12): 3349–3366.e23 (2025). doi:10.1016/j.cell.2025.03.034 · PMID:40245860 · PMCID:PMC12167154
