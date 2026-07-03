#!/usr/bin/env bash
# deMULTIplex2 / Stoeckius 2018 cell-hashing — full local reproduction via multiseq.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="demultiplex2_stoeckius"
D="Benchmarks/_data/$LABEL"
PY=".venv/bin/python"
IGVF=".venv/bin/igvfagent"
mkdir -p "$D"

# pyreadr reads the bundled .RData without needing R.
$PY -c "import pyreadr" 2>/dev/null || $PY -m pip install -q pyreadr

# 1) Fetch the Stoeckius 8-donor PBMC HTO matrix bundled with deMULTIplex2 (public, MIT).
if [ ! -f "$D/stoeckius_pbmc.RData" ]; then
  curl -sL -o "$D/stoeckius_pbmc.RData" \
    "https://raw.githubusercontent.com/Gartner-Lab/deMULTIplex2/main/data/stoeckius_pbmc.RData"
fi

# 2) .RData -> cells x 8-HTO CSV
$PY "$HERE/prep_input.py"

# 3) The reproduction: deMULTIplex2's NB-GLM + EM classifier (IGVFagent port).
$IGVF multiseq demultiplex --input "$D/stoeckius_tags.csv" --label "$LABEL"

# 4) Score classification concordance + figure.
$PY "$HERE/make_figures.py"

echo ""
echo "== deMULTIplex2 / Stoeckius benchmark complete =="
echo "Score with: $PY Benchmarks/concordance.py --benchmark $LABEL"
