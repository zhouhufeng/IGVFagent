#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="deng2024_cortex_mpra"
INPUT="Data/Benchmarks/$LABEL/cortex_mpra_counts.tsv"

if [ ! -f "$INPUT" ]; then
    echo "[$LABEL] Local input not found at $INPUT"
    echo "  Download GSE236018 from GEO and place per-oligo counts at $INPUT, then re-run."
    exit 77
fi

.venv/bin/igvfagent mpra activity --counts "$INPUT" --label "$LABEL"
.venv/bin/igvfagent mpra volcano  --label "$LABEL"

# Optional: ORA on the top-N targets
TARGETS="Data/Benchmarks/$LABEL/top_targets.txt"
if [ -f "$TARGETS" ]; then
    .venv/bin/igvfagent enrich ora --genes "$TARGETS" --label "${LABEL}_pathways"
fi

echo ""
echo "== Deng 2024 cortex MPRA benchmark complete =="
echo "Score with: .venv/bin/python Benchmarks/concordance.py --benchmark $LABEL"
