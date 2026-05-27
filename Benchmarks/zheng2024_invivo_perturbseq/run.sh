#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="zheng2024_invivo_perturbseq"
INPUT="Data/Benchmarks/$LABEL/GSE249416.h5ad"

if [ ! -f "$INPUT" ]; then
    echo "[$LABEL] Local h5ad not found at $INPUT"
    echo "  Download GSE249416 from GEO and convert to h5ad (anndata), then re-run."
    exit 77
fi

.venv/bin/igvfagent sc-analyze pipeline --input "$INPUT" --label "$LABEL"
.venv/bin/igvfagent sc-analyze markers  --input "$INPUT" --label "$LABEL"

echo ""
echo "== Zheng 2024 in-vivo Perturb-seq benchmark complete =="
echo "Score with: .venv/bin/python Benchmarks/concordance.py --benchmark $LABEL"
