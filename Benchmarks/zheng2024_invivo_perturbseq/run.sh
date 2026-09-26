#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="zheng2024_invivo_perturbseq"

# ---- Online step (always runs) ----
# Pull the GEO Series metadata + file inventory for Zheng 2024
igvfagent geo series --gse GSE249416

# ---- Local step (real reproduction, if the paper's own deposit is on disk) ----
RAW="Data/Benchmarks/$LABEL/raw"
ALL_QS="$RAW/GSE249416_Perturb_all.qs"
DERIVED="Data/Benchmarks/$LABEL/derived"
OUT="Docs/SingleCell/$(date +%Y%m%d_%H%M%S)_${LABEL}"

if [ -f "$ALL_QS" ]; then
    mkdir -p "$DERIVED" "$OUT"
    Rscript "$HERE/extract_metadata.R" "$RAW" "$DERIVED"
    python3 "$HERE/reproduce_fig4f.py" "$DERIVED/all_metadata.csv" "$OUT"
    echo ""
    echo "== Zheng 2024 in-vivo Perturb-seq: Fig 4F reproduction complete =="
    echo "   Wrote $OUT/concordance_metrics.json"
else
    echo ""
    echo "== Zheng 2024 — online GEO step complete =="
    echo "Local reproduction skipped — no $ALL_QS on disk."
    echo "See OPERATIONS.md for the download + R environment bootstrap:"
    echo "  1. Download GSE249416_Perturb_all.qs.gz and GSE249416_Perturb_sg.qs.gz"
    echo "     from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249416/suppl/"
    echo "     into $RAW/"
    echo "  2. gunzip -k them (the .qs stream is wrapped in an extra outer gzip)"
    echo "  3. module load R + install qs/SeuratObject (recipe in OPERATIONS.md)"
    echo "  4. re-run this script"
    exit 77
fi
