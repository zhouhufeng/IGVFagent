#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="weinstock2024_cd4_crispr"

# 1. Pull the Perturbation Catalogue summary — landing-page stats
.venv/bin/igvfagent perturb-catalog summary

# 2. Modality-scoped search for the paper's KMT2A focus locus
.venv/bin/igvfagent perturb-catalog search-modality \
    --modality crispr-screen \
    --query KMT2A \
    --dataset-limit 50

# 3. Pull the Weinstock-specific GEO sub-series metadata
.venv/bin/igvfagent geo series --gse GSE171674

echo ""
echo "== Weinstock 2024 CD4 CRISPR benchmark — online steps complete =="
echo "Optional next step (requires the IGVF KG mirror to be warm):"
echo "  .venv/bin/igvfagent network pkn-from-kg --label ${LABEL}_pkn"
echo "  .venv/bin/igvfagent network steiner --seeds KMT2A,STAT5A,IL2 --label ${LABEL}_steiner"
echo "Generate figures with:"
echo "  .venv/bin/python Benchmarks/weinstock2024_cd4_crispr/make_figures.py"
