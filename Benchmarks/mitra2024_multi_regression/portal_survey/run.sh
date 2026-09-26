#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
cd "$ROOT"
# The repo's .venv when it runs on this machine, else the igvfagent on PATH
# (a .venv copied from another platform cannot run here).
BIN="$ROOT/.venv/bin"
"$BIN/igvfagent" --help >/dev/null 2>&1 || BIN="$(dirname "$(readlink -f "$(command -v igvfagent)")")"
LABEL="mitra2024_scarlink"

# Online step: enumerate IGVF multiome AnalysisSets
$BIN/igvfagent multiome retrieve --count 5 --label "$LABEL"

echo ""
echo "== Mitra 2024 SCARlink benchmark — online step complete =="
echo "Optional next step (requires a downloaded 10x Multiome mtx bundle):"
echo "   $BIN/igvfagent multiome process-local --input <path-to-bundle>"
echo "   $BIN/igvfagent multiome peak2gene --label $LABEL"
echo "Score with: $BIN/python Benchmarks/concordance.py --benchmark mitra2024_multi_regression"
