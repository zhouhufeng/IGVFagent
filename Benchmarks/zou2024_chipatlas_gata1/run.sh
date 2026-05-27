#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="zou2024_chipatlas_gata1"

# Browse: how many GATA1 experiments are in ChIP-Atlas for Blood
.venv/bin/igvfagent chipatlas list-antigens \
    --genome hg38 --ag-class "TFs and others" --cell-class Blood --limit 25

# Search by name
.venv/bin/igvfagent chipatlas search --query "GATA1 K562" --genome hg38 --limit 10

# Assemble the all-peaks BED URL for GATA1 / Blood / q=05
.venv/bin/igvfagent chipatlas assemble-bed \
    --genome hg38 --ag-class "TFs and others" --antigen GATA1 \
    --cell-class Blood --qval 05 --label "$LABEL"

echo ""
echo "== Zou 2024 ChIP-Atlas GATA1 hematopoietic benchmark complete =="
echo "Score with: .venv/bin/python Benchmarks/concordance.py --benchmark $LABEL"
