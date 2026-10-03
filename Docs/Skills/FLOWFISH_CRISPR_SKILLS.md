# Skill: CRISPRi Flow-FISH Screen

Use this skill to call regulatory elements from a CRISPRi Flow-FISH screen
(cells sorted by RNA-FISH or protein-FACS readout into N bins per guide,
guide-counts sequenced per bin). The math is a clean-room reimplementation
of [EngreitzLab/CRISPRi-FlowFISH-pipeline](https://github.com/EngreitzLab/CRISPRi-FlowFISH-pipeline)
(MIT, Engreitz Lab 2021), following the published methods in:

- **Fulco CP et al. (2019)** *Nature Genetics* 51:1664–1669 — SI describes
  the per-guide log-normal MLE with EM treatment of an "outside" bin.
- **Nasser J et al. (2021)** *Nature* 593:238–243 — applies the same method
  at scale and defines the element-collapse / `Significant` / `Regulated`
  output convention.

## Commands

```bash
# 1. Discover IGVF Portal CRISPRi Flow-FISH datasets
python3 Scripts/flowfish_crispr_skill.py pull-portal --limit 50 --label survey

# 2. Generate a synthetic screen for smoke testing
python3 Scripts/flowfish_crispr_skill.py simulate \
    --out-dir /tmp/flowfish_smoke --n-elements 40 \
    --guides-per-element 8 --knockdown-frac 0.30 --cells-per-guide 800

# 3. Per-guide log-normal MLE
python3 Scripts/flowfish_crispr_skill.py estimate-effects \
    --counts /tmp/flowfish_smoke/counts.tsv \
    --sortparams /tmp/flowfish_smoke/sortparams.tsv --label demo

# 4. Real-space conversion + null normalization
python3 Scripts/flowfish_crispr_skill.py real-space \
    --input Docs/FlowFISH/<ts>_demo_raw_effects.tsv \
    --target-col target --negative-label negative_control \
    --clamp 5 --label demo

# 5. Per-element collapse + significance
python3 Scripts/flowfish_crispr_skill.py score-elements \
    --effects Docs/FlowFISH/<ts>_demo_real_space.tsv \
    --target-col target --element-col ElementName \
    --negative-label negative_control --min-guides 5 --fdr 0.05 \
    --min-effect 0.10 --label demo
```

## Input schema

`counts.tsv` (TSV, one row per guide):

| Column | Required | Notes |
|---|---|---|
| `GuideID` | yes | Unique guide identifier |
| `target` | yes | Gene symbol / element name / `negative_control` |
| `ElementName` | optional | Mapping to a higher-level regulatory element |
| `Bin1`, `Bin2`, ... | yes | Cell counts per FACS bin |

`sortparams.tsv`: one row per bin, columns `Bin, LowBound, HighBound`
(fluorescence in linear space; the skill converts to log10 internally).

## Output

`*_FullEnhancerScore.tsv` (per regulatory element):

| Column | Meaning |
|---|---|
| `ElementName` | Element identifier |
| `mean_effect` | Mean of `effect_normalized` across guides targeting this element (1.0 = null) |
| `EffectSize` | `1 - mean_effect`, fraction knockdown vs null |
| `log2FC` | `log2(mean_effect)` |
| `n_guides` | Guides used in the test |
| `mean_ctrl, sem_ctrl, n_ctrl` | Negative-control distribution stats |
| `p_mwu, fdr_mwu` | Mann-Whitney U vs negatives, BH-adjusted |
| `p_ttest, fdr_ttest` | Welch t-test vs negatives, BH-adjusted |
| `Significant` | `fdr_mwu < FDR` AND `mean_effect < 1` AND `n_guides ≥ min_guides` |
| `Regulated` | `Significant` AND `EffectSize ≥ min_effect` |

## Notes

- The original pipeline also runs sliding-window smoothing (20-guide /
  ≤750 bp) for tiling screens. Add it as a post-process if needed.
- The qPCR-anchored knockdown calibration (`normalize_flowfish_to_qpcr.py`)
  is not implemented yet — but the `--clamp` and per-element collapse
  steps are sufficient for the common case.
- Replicates: fit the MLE per replicate independently and average the
  resulting `effect_normalized` values per guide before calling
  `score-elements`. (You can do this with `pandas groupby + mean`
  externally; this skill keeps the per-replicate atom and the
  per-element scoring as separate commands.)
