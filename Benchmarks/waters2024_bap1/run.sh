#!/usr/bin/env bash
# Reproduce the Waters 2024 BAP1 SGE mapping via igvfagent.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"

LABEL="waters2024_bap1"

# Waters 2024 BAP1 SGE scoreset URN — VERIFIED 2026-05-28 via the
# MaveDB POST /api/v1/score-sets/search endpoint. The unified scoreset
# (urn:mavedb:00000662-0-1) covers all 17 per-exon sub-experiments and
# has the full 18,108 SNVs the paper reports. PMID 38969833.
URN="${BAP1_URN:-urn:mavedb:00000662-0-1}"

if ! curl -sf "https://api.mavedb.org/api/v1/score-sets/${URN}/scores" > /dev/null; then
    echo "[$LABEL] URN ${URN} not yet verified — see TODO_VERIFY note above."
    echo "  Run the suite's verified smoke-test instead:"
    echo "    bash Benchmarks/matreyek2018_pten_vampseq/run.sh"
    exit 77
fi

.venv/bin/igvfagent mavedb map-scoreset \
    --urn "$URN" \
    --gene BAP1 \
    --label "$LABEL"

echo ""
echo "== Waters 2024 BAP1 SGE benchmark complete =="
echo "Score with: .venv/bin/python Benchmarks/concordance.py --benchmark $LABEL"
