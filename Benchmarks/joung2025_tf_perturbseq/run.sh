#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="joung2025_tf_perturbseq"

.venv/bin/igvfagent perturb-catalog search-modality \
    --modality perturb-seq \
    --query KLF4 \
    --label "$LABEL"

echo ""
echo "== Joung 2025 TF Perturb-seq benchmark — online step complete =="
echo "Optional next step (requires GSE273694 h5ad):"
echo "   .venv/bin/igvfagent sc-analyze pipeline --input <GSE273694.h5ad> --label $LABEL"
echo "Score with: .venv/bin/python Benchmarks/concordance.py --benchmark $LABEL"
