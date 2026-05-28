#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="martyn2025_variant_flowfish"
SIM_DIR="Data/Benchmarks/$LABEL"
mkdir -p "$SIM_DIR"

# ---- Step 1. Online: pull the IGVF Portal MeasurementSet manifest ----
.venv/bin/igvfagent flowfish pull-portal --limit 500 --label "$LABEL"

# ---- Step 2. End-to-end pipeline demo (synthetic 20-element screen) ----
# Generates a Variant-FlowFISH-style guide×bin counts table + sortparams
.venv/bin/igvfagent flowfish simulate \
    --out-dir "$SIM_DIR" \
    --n-elements 20 --guides-per-element 5 \
    --knockdown-frac 0.5 --cells-per-guide 200 --seed 42

# Step 2a: per-guide MLE on (guide × bin) counts
.venv/bin/igvfagent flowfish estimate-effects \
    --counts "$SIM_DIR/counts.tsv" \
    --sortparams "$SIM_DIR/sortparams.tsv" \
    --label "${LABEL}_pipeline"

# Step 2b: real-space rescaling against negative controls
RAW=$(ls -t Docs/FlowFISH/*${LABEL}_pipeline_raw_effects.tsv | head -1)
.venv/bin/igvfagent flowfish real-space --input "$RAW" --label "${LABEL}_pipeline"

# Step 2c: per-element Mann-Whitney + Welch + BH-FDR
EFFECTS=$(ls -t Docs/FlowFISH/*${LABEL}_pipeline*real_space*.tsv | head -1)
.venv/bin/igvfagent flowfish score-elements --effects "$EFFECTS" --label "${LABEL}_pipeline"

# ---- Step 3. (Optional) Paper-data step ----
INPUT="$SIM_DIR/flowfish_counts.tsv"
if [ -f "$INPUT" ]; then
    .venv/bin/igvfagent flowfish estimate-effects \
        --counts "$INPUT" --sortparams "$SIM_DIR/sortparams.tsv" \
        --label "${LABEL}_paper"
    .venv/bin/igvfagent flowfish real-space --input "$(ls -t Docs/FlowFISH/*${LABEL}_paper_raw_effects.tsv | head -1)" --label "${LABEL}_paper"
    .venv/bin/igvfagent flowfish score-elements --effects "$(ls -t Docs/FlowFISH/*${LABEL}_paper*real_space*.tsv | head -1)" --label "${LABEL}_paper"
    echo "Paper-data step ran on $INPUT"
fi

echo ""
echo "== Martyn 2025 Variant-FlowFISH benchmark complete =="
echo "Generate figures with:"
echo "  .venv/bin/python Benchmarks/martyn2025_variant_flowfish/make_figures.py"
