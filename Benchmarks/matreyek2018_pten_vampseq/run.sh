#!/usr/bin/env bash
# Verified MaveDB smoke-test for the IGVFagent reproducibility suite.
# This benchmark is the canary: if it passes, the MaveDB + Ensembl REST
# + codon-mapping pipeline all work end-to-end on your machine.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="matreyek2018_pten_vampseq"

.venv/bin/igvfagent mavedb map-scoreset \
    --urn urn:mavedb:00000013-a-1 \
    --gene PTEN \
    --label "$LABEL"

echo ""
echo "== Matreyek 2018 PTEN VAMP-seq smoke-test complete =="
echo "Score with: .venv/bin/python Benchmarks/concordance.py --benchmark $LABEL"
