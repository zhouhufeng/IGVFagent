#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
# The repo's .venv when it runs on this machine, else the igvfagent on PATH
# (a .venv copied from another platform cannot run here).
BIN="$ROOT/.venv/bin"
"$BIN/igvfagent" --help >/dev/null 2>&1 || BIN="$(dirname "$(readlink -f "$(command -v igvfagent)")")"
LABEL="zou2024_chipatlas_gata1"

# Browse: how many GATA1 experiments are in ChIP-Atlas for Blood
$BIN/igvfagent chipatlas list-antigens \
    --genome hg38 --ag-class "TFs and others" --cell-class Blood --limit 25

# Free-text search — ChIP-Atlas uses exact-token matching, so we try
# the broader "GATA1" rather than "GATA1 K562".
$BIN/igvfagent chipatlas search --query "GATA1" --genome hg38 --limit 10

# Assemble the all-peaks BED URL for GATA1 / Blood / q=05.
# Note: chipatlas assemble-bed does not currently accept --label;
# outputs land at the default Docs/ChIPAtlas/<ts>_assemble/ path.
$BIN/igvfagent chipatlas assemble-bed \
    --genome hg38 --ag-class "TFs and others" --antigen GATA1 \
    --cell-class Blood --qval 05

echo ""
echo "== Zou 2024 ChIP-Atlas GATA1 hematopoietic benchmark complete =="
echo "Score with: $BIN/python Benchmarks/concordance.py --benchmark $LABEL"
