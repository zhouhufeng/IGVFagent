#!/usr/bin/env bash
# Generic Perturbation Catalogue census only. This benchmark's original paper
# citation was fabricated — the real reproduction lives in
# Benchmarks/southard2024_comprehensive_transcription/ (see README.md here).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
cd "$ROOT"
LABEL="joung2025_tf_perturbseq"

# The repo-local .venv may have been built on another machine; fall back to
# whatever igvfagent is on $PATH.
IGVF=".venv/bin/igvfagent"
"$IGVF" --help >/dev/null 2>&1 || IGVF="$(command -v igvfagent)"

# Step 1. Perturbation Catalogue landing-page summary (always reachable).
"$IGVF" perturb-catalog summary

# Step 2. Modality-scoped Perturb-seq query. The upstream Cloud Run backend
# sometimes times out on this endpoint, so a failure doesn't stop the run.
"$IGVF" perturb-catalog search-modality \
    --modality perturb-seq --query KLF4 --dataset-limit 20 || \
    echo "[$LABEL] perturb-seq search timed out (upstream issue)."

echo ""
echo "== $LABEL — catalogue census complete =="
echo "For the real paper reproduction run:"
echo "  bash Benchmarks/southard2024_comprehensive_transcription/run.sh"
