#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
# The repo's .venv when it runs on this machine, else the igvfagent on PATH
# (a .venv copied from another platform cannot run here).
BIN="$ROOT/.venv/bin"
"$BIN/igvfagent" --help >/dev/null 2>&1 || BIN="$(dirname "$(readlink -f "$(command -v igvfagent)")")"
LABEL="buckley2024_vhl"

# Buckley 2024 VHL SGE scoreset URN — VERIFIED 2026-05-28 via the
# MaveDB POST /api/v1/score-sets/search endpoint. PMID 38969834.
URN="${VHL_URN:-urn:mavedb:00000675-a-1}"

if ! curl -sf "https://api.mavedb.org/api/v1/score-sets/${URN}/scores" > /dev/null; then
    echo "[$LABEL] URN ${URN} not yet verified — see TODO_VERIFY note above."
    echo "  Run the suite's verified smoke-test instead:"
    echo "    bash Benchmarks/matreyek2018_pten_vampseq/run.sh"
    exit 77
fi

$BIN/igvfagent mavedb map-scoreset \
    --urn "$URN" \
    --gene VHL \
    --label "$LABEL"

# Cross-reference (catalog calls — anonymous)
$BIN/igvfagent catalog get-entity VHL || true
$BIN/igvfagent catalog find-associations VHL --relationship genetic --limit 10 || true

echo ""
echo "== Buckley 2024 VHL SGE benchmark complete =="
echo "Score with: $BIN/python Benchmarks/concordance.py --benchmark $LABEL"
