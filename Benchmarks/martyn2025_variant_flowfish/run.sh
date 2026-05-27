#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="martyn2025_variant_flowfish"
INPUT="Data/Benchmarks/$LABEL/flowfish_counts.tsv"

# Online step: discover IGVF FlowFISH MeasurementSets
.venv/bin/igvfagent flowfish pull-portal --label "$LABEL"

# Local step: per-paper counts (Engreitz lab supplementary or IGVF AnalysisSet)
if [ ! -f "$INPUT" ]; then
    echo ""
    echo "[$LABEL] Local guide×bin counts table not found at $INPUT"
    echo "  Download Engreitz Variant-FlowFISH counts and place at $INPUT, then re-run."
    exit 77
fi

.venv/bin/igvfagent flowfish estimate-effects --counts "$INPUT" --label "$LABEL"
.venv/bin/igvfagent flowfish real-space --label "$LABEL"
.venv/bin/igvfagent flowfish score-elements --label "$LABEL"

echo ""
echo "== Martyn 2025 Variant-FlowFISH benchmark complete =="
echo "Score with: .venv/bin/python Benchmarks/concordance.py --benchmark $LABEL"
