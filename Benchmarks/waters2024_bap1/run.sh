#!/usr/bin/env bash
# Reproduce the Waters 2024 BAP1 SGE mapping via igvfagent.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"

LABEL="waters2024_bap1"

# TODO_VERIFY: the per-exon BAP1 SGE scoreset URNs (urn:mavedb:00000662-a-*)
# need to be looked up individually on MaveDB — the live `score-sets/<urn>/scores`
# endpoint requires the exact scoreset URN, not the experiment URN. Inspect
# `https://www.mavedb.org/#/experiments/urn:mavedb:00000662-0` to find the
# 17 per-exon scoreset URNs and re-run.
URN="${BAP1_URN:-urn:mavedb:00000662-a-1}"

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
