#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="yao2024_encode4"

# Online step: enumerate ENCODE CRISPRi screens
.venv/bin/igvfagent encode retrieve \
    --assay-type "CRISPR screen" \
    --label "$LABEL" || true

echo ""
echo "== Yao 2024 ENCODE4 CRISPRi manifest pull complete =="
echo "Optional next step (requires local counts table):"
echo "   .venv/bin/igvfagent crispri analyze-local --input <counts.tsv> --label $LABEL"
echo "Score with: .venv/bin/python Benchmarks/concordance.py --benchmark $LABEL"
