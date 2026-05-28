# OPERATIONS — Martyn 2025 Variant-FlowFISH (★ direct absorption test)

> **The clearest reproducibility statement in the suite.** IGVFagent's
> `flowfish` skill is a clean-room reimplementation of the **exact**
> Engreitz-lab Variant-FlowFISH pipeline this paper ran. If
> IGVFagent's `flowfish` reproduces Martyn 2025's published effects
> within 10–20 %, we've proved the absorption was faithful.

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

## Two stages

| Stage | Mode | What it tests |
|---|---|---|
| 1. Online discovery | always runs | `flowfish pull-portal` enumerates IGVF Portal FlowFISH MeasurementSets |
| 2. Local analytical | requires `flowfish_counts.tsv` | `flowfish estimate-effects / real-space / score-elements` — the canonical Fulco/Nasser log-normal bin-MLE pipeline |

## Stage 1 — quick run, ~10 s

```bash
bash Benchmarks/martyn2025_variant_flowfish/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark martyn2025_variant_flowfish
```

## What Stage 1 does

```bash
.venv/bin/igvfagent flowfish pull-portal --label martyn2025_variant_flowfish
```

Queries IGVF Portal for FlowFISH MeasurementSets, including the 15
analysis-sets + 11 construct-library-sets + 3 measurement-sets the
Engreitz lab deposited for Martyn 2025 (examples: `IGVFDS1003XTAF`,
`IGVFDS1132MKHT`, `IGVFDS2207IIRU`, `IGVFDS2374NXLW`).

## Where Stage 1 artefacts land

```
Docs/FlowFISH/<ts>_martyn2025_variant_flowfish_*_report.md
Data/Manifests/FlowFISH/<ts>_martyn2025_variant_flowfish_manifest.csv
```

## Stage 2 — local input

```
Data/Benchmarks/martyn2025_variant_flowfish/flowfish_counts.tsv
```

### How to obtain it

**Option A — IGVF Portal batch-download** (recommended)

```bash
mkdir -p Data/Benchmarks/martyn2025_variant_flowfish

# From Stage 1's manifest, pick one of the measurement-sets and
# batch-download its files:
.venv/bin/igvfagent portal batch-download \
    --type MeasurementSet \
    --field-filters "accession=IGVFDS2207IIRU" \
    --fetch --label martyn2025_inputs

# Find the per-(guide, bin) counts table among the downloaded files
# (expect a TSV with columns like guide, bin1, bin2, ..., bin6)
ls Docs/PortalQuery/*_batchdl_martyn2025_inputs/files/
cp Docs/PortalQuery/*_batchdl_martyn2025_inputs/files/*counts*.tsv \
   Data/Benchmarks/martyn2025_variant_flowfish/flowfish_counts.tsv
```

**Option B — Engreitz lab supplementary**

The Variant-FlowFISH paper's supplementary file (Supplementary Table 5
in Martyn 2025) contains the per-variant guide×bin count matrix. Download
the Excel/CSV and convert:

```bash
# After downloading Martyn2025_supplementary_table_5.xlsx:
.venv/bin/python -c "
import openpyxl, csv
wb = openpyxl.load_workbook('Martyn2025_supplementary_table_5.xlsx')
ws = wb['guide_bin_counts']  # confirm sheet name
with open('flowfish_counts.tsv', 'w') as fh:
    w = csv.writer(fh, delimiter='\t')
    for row in ws.iter_rows(values_only=True):
        w.writerow(row)
"
```

**Option C — Synthetic for pipeline-check only**

```bash
.venv/bin/igvfagent flowfish simulate \
    --out-dir Data/Benchmarks/martyn2025_variant_flowfish/synthetic \
    --n-elements 50 --guides-per-element 8 \
    --knockdown-frac 0.30 --cells-per-guide 800 --seed 7

cp Data/Benchmarks/martyn2025_variant_flowfish/synthetic/guide_bin_counts.tsv \
   Data/Benchmarks/martyn2025_variant_flowfish/flowfish_counts.tsv
```

This generates a synthetic 50-element × 8-guide table with a 30 %
knockdown effect — exercises the full pipeline but won't reproduce
paper-specific signals.

## Stage 2 commands (after counts file exists)

```bash
bash Benchmarks/martyn2025_variant_flowfish/run.sh
```

Which fires:

```bash
.venv/bin/igvfagent flowfish estimate-effects --counts <INPUT> --label martyn2025_variant_flowfish
.venv/bin/igvfagent flowfish real-space      --label martyn2025_variant_flowfish
.venv/bin/igvfagent flowfish score-elements  --label martyn2025_variant_flowfish
```

## Where Stage 2 artefacts land

```
Docs/FlowFISH/<ts>_martyn2025_variant_flowfish/
├── per_guide_effects.tsv       log-normal MLE per guide
├── per_guide_realspace.tsv     converted to fold-change vs negative controls
├── per_element_scores.tsv      Mann-Whitney + Welch + BH-FDR per element
├── per_element_calls.tsv       Significant / Regulated / NotSignificant labels
├── volcano.png
└── summary.json
```

## Stage 2 ground-truth spot-checks

| Signal | Expected (paper) |
|---|---|
| Number of guides analyzed | ~5,000+ (paper used a variant-tiling library) |
| Significantly-regulated elements | Tens of elements per locus |
| GATA1 enhancer effects | Paper Fig 3 — known +24/+58 kb HS sites should appear as "Regulated" |
| MYC enhancer effects | Paper Fig 5 — variant edits at the well-characterized MYC enhancer (chr8:128.74M) should show effects within 10–20% of paper |

After Stage 2:

```bash
# Check that the +58 kb GATA1 enhancer region shows up as Regulated
.venv/bin/python -c "
import csv
for r in csv.DictReader(open('Docs/FlowFISH/<ts>_*/per_element_calls.tsv'), delimiter='\t'):
    if 'GATA1' in r.get('element', ''):
        print(r)"
```

## Running through the UI

`flowfish_pull_portal` is NOT registered. `flowfish_estimate_effects`,
`flowfish_real_space`, `flowfish_score_elements` ARE registered.

### Stage 1 via UI

Use the shell for `flowfish pull-portal`. Or as a proxy:

```
List IGVF Portal FlowFISH MeasurementSets via portal_search:
  - type=MeasurementSet
  - field_filters="preferred_assay_titles=Flow-FISH,CRISPRi Flow-FISH"
  - limit=25
  - label="martyn2025_variant_flowfish"

Then for the top 3 results, call portal_get on each to fetch the
@id of associated count tables. Report a table of (accession,
lab, biosample, file count).
```

### Stage 2 via UI

Once `flowfish_counts.tsv` is in place:

```
Run the Martyn 2025 Variant-FlowFISH analytical pipeline. The counts
table is at Data/Benchmarks/martyn2025_variant_flowfish/flowfish_counts.tsv.

Step 1: flowfish_estimate_effects with counts at that path,
        label="martyn2025_variant_flowfish".

Step 2: flowfish_real_space with that label.

Step 3: flowfish_score_elements with that label.

Report:
  - total guides processed
  - elements called Regulated vs NotSignificant (per Fulco/Nasser convention)
  - top 5 elements by effect magnitude
  - whether the GATA1 +24 kb / +58 kb enhancers appear as Regulated
```

UI sidebar: max iterations = 18, temperature = 0.0.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `flowfish estimate-effects` errors with "guide column not found" | Counts file has wrong column names | Expected schema: `guide, bin1, bin2, ..., binN, [target_element]`; remap with `awk` if needed |
| `real-space` step gives all-zero fold-changes | No negative-control guides labelled | Add `--neg-controls non-targeting` (or whatever your label is) |
| `score-elements` reports 0 Regulated | Either too few guides per element, or per-guide effects too noisy | Check `n_guides_per_element` in summary.json; needs ≥ 4 |
| Synthetic results look fine but paper data fails | Real data has a meta-data column structure that the synthetic data lacks | Inspect `flowfish_counts.tsv` columns; consult `igvfagent flowfish simulate --help` for the canonical schema |

## License + provenance

* **Paper data**: IGVF Portal AnalysisSets (Engreitz lab) + GEO — public.
* **Code**: IGVFagent Apache-2.0 over upstream MIT (EngreitzLab/CRISPRi-FlowFISH-pipeline).
* **Citation**: Martyn GE et al. *Cell* **188**: 3349–3366 (2025).
  doi:10.1016/j.cell.2025.03.034 · PMID:40245860
