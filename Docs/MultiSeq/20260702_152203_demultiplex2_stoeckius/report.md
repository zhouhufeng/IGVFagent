# MULTI-seq demultiplexing — demultiplex2_stoeckius

- Input: `Benchmarks/_data/demultiplex2_stoeckius/stoeckius_tags.csv`
- Cells: **15,113**, tags: **8** (HTO_A, HTO_B, HTO_C, HTO_D, HTO_E, HTO_F, HTO_G, HTO_H)

## Droplet-type counts

| Type | n | % |
|---|---|---|
| singlet | 12,591 | 83.3% |
| multiplet | 2,471 | 16.4% |
| negative | 51 | 0.3% |

## Per-sample singlet counts

| Tag | Singlets |
|---|---|
| HTO_A | 1,693 |
| HTO_B | 1,705 |
| HTO_C | 1,633 |
| HTO_D | 1,518 |
| HTO_E | 1,350 |
| HTO_F | 1,387 |
| HTO_G | 1,631 |
| HTO_H | 1,674 |

## Plots
- `Plots/tag_histogram.png` — bimodal log10 UMI histogram per tag
- `Plots/tag_heatmap.png` — mean log10(UMI+1) by call group
- `Plots/diagnostics_*.png` — per-tag 4-panel diagnostics
