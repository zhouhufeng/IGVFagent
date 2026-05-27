#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="weinstock2024_cd4_crispr"

.venv/bin/igvfagent perturb-catalog search-modality \
    --modality crispr-screen \
    --query KMT2A \
    --label "$LABEL"

# Build a context-specific PKN + Steiner subnetwork — exercises the network skill
.venv/bin/igvfagent network pkn-from-kg --label "${LABEL}_pkn" || true

echo ""
echo "== Weinstock 2024 CD4 CRISPR network benchmark complete =="
echo "Score with: .venv/bin/python Benchmarks/concordance.py --benchmark $LABEL"
