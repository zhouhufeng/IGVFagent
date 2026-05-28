#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="zou2024_chipatlas_gata1"

# Browse: how many GATA1 experiments are in ChIP-Atlas for Blood
.venv/bin/igvfagent chipatlas list-antigens \
    --genome hg38 --ag-class "TFs and others" --cell-class Blood --limit 25

# Free-text search — ChIP-Atlas uses exact-token matching, so we try
# the broader "GATA1" rather than "GATA1 K562".
.venv/bin/igvfagent chipatlas search --query "GATA1" --genome hg38 --limit 10

# Assemble the all-peaks BED URL for GATA1 / Blood / q=05.
# Note: chipatlas assemble-bed does not currently accept --label;
# outputs land at the default Docs/ChIPAtlas/<ts>_assemble/ path.
.venv/bin/igvfagent chipatlas assemble-bed \
    --genome hg38 --ag-class "TFs and others" --antigen GATA1 \
    --cell-class Blood --qval 05

echo ""
echo "== Zou 2024 ChIP-Atlas GATA1 hematopoietic benchmark complete =="
echo "Score with: .venv/bin/python Benchmarks/concordance.py --benchmark $LABEL"
