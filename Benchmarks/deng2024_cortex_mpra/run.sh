#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
# The repo's .venv when it runs on this machine, else the igvfagent on PATH
# (a .venv copied from another platform cannot run here).
BIN="$ROOT/.venv/bin"
"$BIN/igvfagent" --help >/dev/null 2>&1 || BIN="$(dirname "$(readlink -f "$(command -v igvfagent)")")"
LABEL="deng2024_cortex_mpra"

# Step 1 (online): IGVF Portal MPRA-class manifest
$BIN/igvfagent mpra portal-manifest --limit 200 --label "$LABEL"

# Step 2 (online): Perturbation Catalogue summary (places MPRA/MAVE
# modality in context). MAVE = catalogue's MPRA-family modality.
$BIN/igvfagent perturb-catalog summary

# Step 3 (online, NEW): Synapse retrieval — Deng's primary deposit lives
# at syn21392931 (PsychENCODE NeuREs study). Pulls anonymous metadata +
# enumerates the MPRA_CapstoneII sub-folder which contains the paper's
# DNA+RNA count fastqs. The download step needs SYNAPSE_AUTH_TOKEN
# once the PsychENCODE Data-Use Agreement is accepted.
$BIN/igvfagent synapse entity --syn syn21392931 --annotations
$BIN/igvfagent synapse walk --syn syn21392931 --max-depth 3 \
    --max-children 50 --label "${LABEL}_neures"
$BIN/igvfagent synapse children --syn syn51090452 \
    --label "${LABEL}_mpra_capstone2"

# Step 3 (local): activity + volcano if the published count table is on disk.
# Deng 2024 deposits via PsychENCODE Synapse (syn21392931); the per-oligo
# DNA + RNA count table needs to be fetched manually + placed at INPUT.
INPUT="Data/Benchmarks/$LABEL/cortex_mpra_counts.tsv"
TARGETS="Data/Benchmarks/$LABEL/top_targets.txt"
if [ -f "$INPUT" ]; then
    $BIN/igvfagent mpra activity --counts "$INPUT" --label "$LABEL"
    $BIN/igvfagent mpra volcano  --label "$LABEL"
    if [ -f "$TARGETS" ]; then
        $BIN/igvfagent enrich ora --genes "$TARGETS" --label "${LABEL}_pathways"
    fi
    echo "Local analytical step ran on $INPUT"
fi

echo ""
echo "== Deng 2024 cortex lentiMPRA benchmark — online steps complete =="
echo "Generate figures with:"
echo "  $BIN/python Benchmarks/deng2024_cortex_mpra/make_figures.py"
