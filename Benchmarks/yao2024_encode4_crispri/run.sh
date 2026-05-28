#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="yao2024_encode4"

# Online step: enumerate every ENCODE FunctionalCharacterizationExperiment
# CRISPR-screen FCE. Yao 2024 paper headline: 108 noncoding CRISPRi screens
# (K562 + HepG2 + Jurkat). The live ENCODE total is now larger because the
# database has grown post-paper-cutoff — see Concordance section in README.md.
.venv/bin/igvfagent encode retrieve \
    --assay "CRISPR screen" --limit 500 --label "$LABEL"

echo ""
echo "== Yao 2024 ENCODE4 CRISPRi manifest pull complete =="
echo "Manifest under Data/Manifests/ENCODE/*_${LABEL}_experiments.csv"
echo "Generate figures with:"
echo "    .venv/bin/python Benchmarks/yao2024_encode4_crispri/make_figures.py"
echo "Score with:"
echo "    .venv/bin/python Benchmarks/concordance.py --benchmark yao2024_encode4_crispri"
