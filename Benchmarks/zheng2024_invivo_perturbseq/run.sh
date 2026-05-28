#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="zheng2024_invivo_perturbseq"

# ---- Online step (always runs) ----
# Pull the GEO Series metadata + file inventory for Zheng 2024
.venv/bin/igvfagent geo series --gse GSE249416

# ---- Local step (if h5ad is on disk) ----
INPUT="Data/Benchmarks/$LABEL/GSE249416.h5ad"
if [ -f "$INPUT" ]; then
    .venv/bin/igvfagent sc-analyze pipeline --input "$INPUT" --label "$LABEL"
    .venv/bin/igvfagent sc-analyze markers  --input "$INPUT" --label "$LABEL"
    echo ""
    echo "== Zheng 2024 in-vivo Perturb-seq benchmark complete (with local h5ad) =="
else
    echo ""
    echo "== Zheng 2024 — online GEO step complete =="
    echo "Local sc-analyze step skipped — no $INPUT on disk."
    echo "To download + convert R Seurat -> h5ad:"
    echo "  .venv/bin/igvfagent geo download --gse GSE249416 --pattern '*Perturb*qs.gz'"
    echo "  (then in R: qs::qread, SeuratObject::as.SingleCellExperiment, zellkonverter::writeH5AD)"
fi
echo "Generate figures with:"
echo "  .venv/bin/python Benchmarks/zheng2024_invivo_perturbseq/make_figures.py"
