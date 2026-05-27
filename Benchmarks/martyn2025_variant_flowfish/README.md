# Martyn 2025 — Variant-FlowFISH (the canonical FlowFISH reproducibility test)

## Citation

Martyn GE, ..., Engreitz JM. *Cell* **188**: 3349–3366 (2025).
DOI: [10.1016/j.cell.2025.03.034](https://doi.org/10.1016/j.cell.2025.03.034) · PMID: 40245860

## Data

* **IGVF portal**: 15 analysis-sets + 11 construct-library-sets + 3 measurement-sets (Engreitz lab). Examples: `IGVFDS1003XTAF`, `IGVFDS1132MKHT`, `IGVFDS2207IIRU`, `IGVFDS2374NXLW`
* **GitHub**: [EngreitzLab/Variant-FlowFISH](https://github.com/EngreitzLab/Variant-FlowFISH)

## Headline workflow (paper)

1. Variant tiling with a paired guide-and-edit construct library.
2. Sort cells across expression bins; sequence amplicons to recover per-(variant, bin) counts.
3. Per-(variant, allele) log-normal MLE → fold-change vs negative controls → Mann-Whitney + BH-FDR per element.

## Why this matters

IGVFagent's `flowfish` skill is **the clean-room reimplementation of exactly this pipeline**. Reproducing Martyn 2025's published effects through `flowfish` is the most direct test of the absorption's faithfulness.

## What IGVFagent reproduces

* **Online step**: `flowfish pull-portal` enumerates the IGVF FlowFISH MeasurementSets.
* **Local step** (requires the Engreitz lab's published counts table at `Data/Benchmarks/martyn2025_variant_flowfish/flowfish_counts.tsv`):
  * `flowfish estimate-effects` → per-guide log-normal MLE on (guide × bin) counts
  * `flowfish real-space` → fold-change + negative-control rescaling
  * `flowfish score-elements` → Mann-Whitney + Welch + BH-FDR

## Ground truth to spot-check

* GATA1 / MYC enhancer effects should agree within 10-20 % of Martyn's reported values.
* `score-elements` output should call several elements `Significant` and `Regulated`.
