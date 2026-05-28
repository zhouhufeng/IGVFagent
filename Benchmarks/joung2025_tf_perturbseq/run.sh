#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="joung2025_tf_perturbseq"

# Step 1. Perturbation Catalogue landing-page summary (always reachable).
# Reports the per-modality dataset counts including the Perturb-seq universe
# that Joung 2025 contributes to once their embargo lifts.
.venv/bin/igvfagent perturb-catalog summary

# Step 2. Modality-scoped Perturb-seq query.
# NOTE: the upstream Google-Cloud-Run backend frequently times out on the
# perturb-seq endpoint — we wrap it in `|| true` so the benchmark continues
# even on a transient 504. The summary above gives us the catalogue total.
.venv/bin/igvfagent perturb-catalog search-modality \
    --modality perturb-seq --query KLF4 --dataset-limit 20 || \
    echo "[$LABEL] perturb-seq search timed out (upstream issue); see Honest caveats."

# Step 3. (Optional) Run sc-analyze on a local Joung h5ad once SCP2169 is
# downloaded.
INPUT="Data/Benchmarks/$LABEL/joung_tf_perturbseq.h5ad"
if [ -f "$INPUT" ]; then
    .venv/bin/igvfagent sc-analyze pipeline --input "$INPUT" --label "$LABEL"
fi

echo ""
echo "== Joung 2025 TF Perturb-seq benchmark — online steps complete =="
echo "Generate figures with:"
echo "  .venv/bin/python Benchmarks/joung2025_tf_perturbseq/make_figures.py"
