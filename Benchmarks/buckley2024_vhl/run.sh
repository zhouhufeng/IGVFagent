#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="buckley2024_vhl"

# TODO_VERIFY: look up the canonical VHL SGE scoreset URN on
# https://www.mavedb.org and override via VHL_URN env var, or edit below.
URN="${VHL_URN:-urn:mavedb:00001183-a-1}"

if ! curl -sf "https://api.mavedb.org/api/v1/score-sets/${URN}/scores" > /dev/null; then
    echo "[$LABEL] URN ${URN} not yet verified — see TODO_VERIFY note above."
    echo "  Run the suite's verified smoke-test instead:"
    echo "    bash Benchmarks/matreyek2018_pten_vampseq/run.sh"
    exit 77
fi

.venv/bin/igvfagent mavedb map-scoreset \
    --urn "$URN" \
    --gene VHL \
    --label "$LABEL"

# Cross-reference (catalog calls — anonymous)
.venv/bin/igvfagent catalog get-entity VHL || true
.venv/bin/igvfagent catalog find-associations VHL --relationship genetic --limit 10 || true

echo ""
echo "== Buckley 2024 VHL SGE benchmark complete =="
echo "Score with: .venv/bin/python Benchmarks/concordance.py --benchmark $LABEL"
