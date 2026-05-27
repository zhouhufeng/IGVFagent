#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="mitra2024_scarlink"

# Online step: enumerate IGVF multiome AnalysisSets
.venv/bin/igvfagent multiome retrieve --count 5 --label "$LABEL"

echo ""
echo "== Mitra 2024 SCARlink benchmark — online step complete =="
echo "Optional next step (requires a downloaded 10x Multiome mtx bundle):"
echo "   .venv/bin/igvfagent multiome process-local --input <path-to-bundle>"
echo "   .venv/bin/igvfagent multiome peak2gene --label $LABEL"
echo "Score with: .venv/bin/python Benchmarks/concordance.py --benchmark $LABEL"
