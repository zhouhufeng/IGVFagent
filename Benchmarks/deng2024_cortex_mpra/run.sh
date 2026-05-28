#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="deng2024_cortex_mpra"

# Step 1 (online): IGVF Portal MPRA-class manifest
.venv/bin/igvfagent mpra portal-manifest --limit 200 --label "$LABEL"

# Step 2 (online): Perturbation Catalogue summary (places MPRA/MAVE
# modality in context). MAVE = catalogue's MPRA-family modality.
.venv/bin/igvfagent perturb-catalog summary

# Step 3 (online, NEW): Synapse retrieval — Deng's primary deposit lives
# at syn21392931 (PsychENCODE NeuREs study). Pulls anonymous metadata +
# enumerates the MPRA_CapstoneII sub-folder which contains the paper's
# DNA+RNA count fastqs. The download step needs SYNAPSE_AUTH_TOKEN
# once the PsychENCODE Data-Use Agreement is accepted.
.venv/bin/igvfagent synapse entity --syn syn21392931 --annotations
.venv/bin/igvfagent synapse walk --syn syn21392931 --max-depth 3 \
    --max-children 50 --label "${LABEL}_neures"
.venv/bin/igvfagent synapse children --syn syn51090452 \
    --label "${LABEL}_mpra_capstone2"

# Step 3 (local): activity + volcano if the published count table is on disk.
# Deng 2024 deposits via PsychENCODE Synapse (syn21392931); the per-oligo
# DNA + RNA count table needs to be fetched manually + placed at INPUT.
INPUT="Data/Benchmarks/$LABEL/cortex_mpra_counts.tsv"
TARGETS="Data/Benchmarks/$LABEL/top_targets.txt"
if [ -f "$INPUT" ]; then
    .venv/bin/igvfagent mpra activity --counts "$INPUT" --label "$LABEL"
    .venv/bin/igvfagent mpra volcano  --label "$LABEL"
    if [ -f "$TARGETS" ]; then
        .venv/bin/igvfagent enrich ora --genes "$TARGETS" --label "${LABEL}_pathways"
    fi
    echo "Local analytical step ran on $INPUT"
fi

echo ""
echo "== Deng 2024 cortex lentiMPRA benchmark — online steps complete =="
echo "Generate figures with:"
echo "  .venv/bin/python Benchmarks/deng2024_cortex_mpra/make_figures.py"
