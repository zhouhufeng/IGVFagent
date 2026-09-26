#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="weinstock2024_cd4_crispr"

BIN="$ROOT/.venv/bin/igvfagent"
"$BIN" --version >/dev/null 2>&1 || BIN="$(command -v igvfagent)"

# 1. Pull the Perturbation Catalogue summary + KMT2A modality search into a
#    single labelled run dir (Docs/Perturbation/<ts>_<label>/), which is what
#    Benchmarks/concordance.py expects to find.
"$BIN" perturb-catalog pipeline --gene KMT2A --label "$LABEL" --dataset-limit 50

# 2. Pull the Weinstock-specific GEO sub-series metadata
"$BIN" geo series --gse GSE171674

echo ""
echo "== Weinstock 2024 CD4 CRISPR benchmark — online steps complete =="
echo "Optional next step (requires the IGVF KG mirror to be warm):"
echo "  $BIN network pkn-from-kg --label ${LABEL}_pkn"
echo "  $BIN network steiner --seeds KMT2A,STAT5A,IL2 --label ${LABEL}_steiner"
echo "Generate figures with:"
echo "  python3 Benchmarks/weinstock2024_cd4_crispr/make_figures.py"
