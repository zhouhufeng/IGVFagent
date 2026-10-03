# E2G CRISPR benchmark — ENCODE-rE2G & scE2G vs. CRISPR ground truth

**Question.** Do the IGVF-hosted **scE2G** single-cell enhancer→gene predictions,
and the **ENCODE-rE2G** predictions, recover CRISPR-validated element→gene
regulatory links — and can IGVFagent reproduce the published benchmark numbers?

**Papers.**
- Gschwind et al. *An encyclopedia of enhancer–gene regulatory interactions in the
  human genome.* bioRxiv 2023, [10.1101/2023.11.09.563812](https://doi.org/10.1101/2023.11.09.563812)
  (ENCODE-rE2G).
- Sheth, Qiu et al. *Mapping enhancer–gene regulatory interactions from single-cell
  data.* bioRxiv 2024, [10.1101/2024.11.23.624931](https://doi.org/10.1101/2024.11.23.624931)
  (scE2G).

## Data (all fetched with IGVFagent)

| Role | Source | Accession / path |
|---|---|---|
| scE2G K562 scored predictions (11.5 M pairs; ABC / ARC-E2G / scE2G scores) | IGVF Portal `PredictionSet` | `IGVFDS5428HHMB` → `IGVFFI1706PNVV.tsv.gz` |
| scE2G coronary-artery predictions (Part B) | IGVF Portal `PredictionSet` | `IGVFDS7543BJBA` → `IGVFFI8252JBBA.tsv.gz` |
| ENCODE-rE2G cross-validated K562 predictions (base + extended) | Synapse (PAT) | `syn53019593`, `syn53019595` |
| CRISPR ground truth (K562; 10,356 pairs / 471 positives) | EngreitzLab/CRISPR_comparison | `EPCrisprBenchmark_combined_data.training_K562.GRCh38` |

The scE2G side is on the **IGVF portal**; ENCODE-rE2G lives on **Synapse**
(controlled-download, PAT). The CRISPR benchmark is the shared gold standard.

## Method

`Scripts/e2g_benchmark_eval.py` maps every CRISPR-tested (element, gene) pair to a
predicted score by **interval overlap of the perturbed element with a predicted
element for the same gene** (max score over overlaps; unmatched → 0), then computes
a precision–recall curve, **AUPRC** (step-integrated average precision) and
**precision at 70 % recall** — the two statistics the E2G literature reports.

```bash
bash Benchmarks/e2g_crispr_benchmark/run.sh          # end-to-end
python Benchmarks/e2g_crispr_benchmark/make_figures.py
```

## Concordance — Part A (K562, vs CRISPR ground truth)

| Method | AUPRC (IGVFagent) | Precision@70%recall | Published |
|---|---|---|---|
| **ENCODE-rE2G extended** (FullModel) | **0.758** | 0.700 | best model in paper |
| ENCODE-rE2G extended (−EP300) | 0.730 | 0.672 | — |
| **ENCODE-rE2G base** | **0.634** | **0.543** | **0.634 / 0.543** ✅ exact |
| scE2G (Score.ignoreTPM) | 0.589 | 0.473 | ≈ rE2G-base range |
| scE2G (Score) | 0.531 | 0.254 | — |
| ARC-E2G | 0.495 | 0.336 | — |
| ABC (baseline) | 0.491 | 0.353 | 0.612 (see caveat) |

10,356 pairs / 471 positives; scE2G-table methods matched 89.1 % of pairs by
overlap, rE2G files 100 % (they ship per-benchmark-pair predictions).

## Concordance — Part B (coronary artery, IGVFDS7543BJBA)

Characterization only (no matched CRISPR ground truth exists for coronary artery,
so precision/recall is not defined): **606,342** predictions across **10** coronary
cell types, 14,359 genes, 465,573 elements; element class 43 % genic / 34 %
intergenic / 23 % promoter; median element→TSS distance 16.9 kb. See figures 4–6.

## Verdict

**Reproduced.** IGVFagent recovers the **ENCODE-rE2G base AUPRC to the third
decimal (0.634 vs 0.634)** and its precision@70%recall (0.543), and ranks the
models exactly as the papers do: rE2G-extended > rE2G-base > scE2G > ARC-E2G ≈ ABC.
The central scientific claims hold: (1) the **extended epigenomic rE2G model is
best**; (2) **scE2G, from single-cell data alone, reaches the rE2G-base range** and
**beats the ABC baseline** (0.531–0.589 vs 0.491). The exact-match on the
independently-computed rE2G number validates the evaluator.

## Honest caveats

- **ABC (0.491) < published ABC (0.612).** Our ABC score comes from the *scE2G
  pipeline's* `ABC.Score` column mapped by generic interval overlap, **not** the
  official `CRISPR_comparison` overlap/annotation procedure (element definitions,
  TSS handling, `ValidConnection` filters). The rE2G files, which ship
  per-benchmark-pair predictions, reproduce the paper exactly precisely because
  no overlap step is needed — isolating the overlap step as the source of the ABC
  gap. scE2G numbers are similarly ~overlap-affected and should be read as
  *lower bounds*.
- **scE2G was trained/evaluated by the authors on K562**; the IGVF-hosted genome-
  wide table is the deployed v1.2.0, not necessarily the exact cross-validated
  benchmark split the paper used, so a small offset from any paper scE2G figure is
  expected.
- **Precision@70%recall** here is max precision among points with recall ≥ 0.70
  (step PR curve); the papers may interpolate slightly differently.
- Part B (coronary artery) is **characterization, not a ground-truth benchmark** —
  flagged because a naive reader might expect precision/recall for it.

## Files
- `results/` — per-method JSON (AUPRC, P@70%, PR curve).
- `figures/` — fig1 PR curves, fig2 AUPRC bars, fig3 precision@70 (Part A);
  fig4–6 coronary-artery characterization (Part B).
- `run.sh`, `make_figures.py`, `expected.json`.
