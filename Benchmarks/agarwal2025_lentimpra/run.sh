#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
# The repo's .venv when it runs on this machine, else the igvfagent on PATH
# (a .venv copied from another platform cannot run here).
BIN="$ROOT/.venv/bin"
"$BIN/igvfagent" --help >/dev/null 2>&1 || BIN="$(dirname "$(readlink -f "$(command -v igvfagent)")")"
LABEL="agarwal2025_lentimpra"
INPUT="Data/Benchmarks/$LABEL/oligo_counts.tsv"

# Online step: discover IGVF MPRA datasets
$BIN/igvfagent mpra portal-manifest --limit 50 --label "$LABEL" || true

# Local step is gated on the user having downloaded the count table.
if [ ! -f "$INPUT" ]; then
    echo ""
    echo "[$LABEL] Local input not found at $INPUT"
    echo "  Download per-oligo DNA+RNA counts from GEO GSE142696 (or an IGVF MPRA"
    echo "  AnalysisSet from the manifest above) and place it at $INPUT, then re-run."
    exit 77
fi

$BIN/igvfagent mpra activity --counts "$INPUT" --label "$LABEL"
$BIN/igvfagent mpra qc       --counts "$INPUT" --label "$LABEL"
$BIN/igvfagent mpra volcano  --label "$LABEL"

echo ""
echo "== Agarwal 2025 lentiMPRA benchmark complete =="
echo "Score with: $BIN/python Benchmarks/concordance.py --benchmark $LABEL"
