#!/usr/bin/env bash
# Wang 2026 Spatial-ATAC-Hi-C — full local reproduction via the
# `spatial-hic` skill against planted ground truth.
#
# Runs in seconds and needs no network: `prep_input.py` synthesises a
# 36-pixel dataset whose contact statistics are calibrated to the
# paper's reported bands and which carries a known CN gain, a known
# clone-specific loop, known A/B blocks and a known TSS enrichment.
# The chain then has to find all of them.
#
# The optional GEO leg (--with-geo) additionally verifies that the real
# GSE307620 series is reachable and enumerates its deposits.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"; cd "$ROOT"

LABEL="wang2026_spatial_atac_hic"
D="Benchmarks/_data/$LABEL"
PY="${PY:-python3}"
IGVF="${IGVF:-$PY -m igvfagent.cli}"
WITH_GEO=0
[ "${1:-}" = "--with-geo" ] && WITH_GEO=1

# 1) Build the input with exact ground truth.
$PY "$HERE/prep_input.py"

# 2) Demultiplex the barcoded pairs into the spatial pixel grid.
$IGVF spatial-hic pixel-demux \
    --pairs "$D/sample.pairs.gz" \
    --barcode-a "$D/barcodes_A.txt" --barcode-b "$D/barcodes_B.txt" \
    --layout BA --label "${LABEL}_demux"

PIX="$(ls -dt Docs/SpatialATACHiC/2*_${LABEL}_demux/pixels | head -1)"
echo "pixels: $PIX"

# 3) Per-pixel contact QC + TSS enrichment.
$IGVF spatial-hic qc --pairs-dir "$PIX" \
    --chrom-sizes "$D/chrom.sizes" \
    --fragments "$D/fragments.tsv.gz" --gene-model "$D/genes.gtf" \
    --min-contacts 100 --label "${LABEL}_qc"

# 4) Gene-level scores in both modalities, on a shared pixel key.
$IGVF spatial-hic gas --fragments "$D/fragments.tsv.gz" \
    --gene-model "$D/genes.gtf" \
    --barcode-a "$D/barcodes_A.txt" --barcode-b "$D/barcodes_B.txt" \
    --min-cells 1 --label "${LABEL}_gas"
$IGVF spatial-hic gad --pairs-dir "$PIX" \
    --gene-model "$D/genes.gtf" --min-cells 1 --label "${LABEL}_gad"

# 5) A/B compartments on chr2, oriented by gene density.
$IGVF spatial-hic compartment --pairs-dir "$PIX" \
    --chrom-sizes "$D/chrom.sizes" --chroms chr2 --resolution 200000 \
    --gene-model "$D/genes.gtf" --label "${LABEL}_compartment"

# 6) Copy number: pseudobulk + per-pixel clone structure.
$IGVF spatial-hic cnv --pairs-dir "$PIX" \
    --chrom-sizes "$D/chrom.sizes" --resolution 500000 \
    --per-pixel --segment --smooth --min-pixel-contacts 50 \
    --label "${LABEL}_cnv"

# 7) Loops: quantify, ANOVA across clones, APA pileup.
$IGVF spatial-hic loops --pairs-dir "$PIX" \
    --chrom-sizes "$D/chrom.sizes" --bedpe "$D/loops.bedpe" \
    --clusters "$D/clusters.tsv" --resolution 10000 --apa \
    --label "${LABEL}_loops"

# 8) Spatial rendering of a marker gene.
GAS="$(ls -t Docs/SpatialATACHiC/2*_${LABEL}_gas/gene_activity_score.tsv | head -1)"
$IGVF spatial-hic viz --table "$GAS" --column GENEX --label "${LABEL}_viz"

# 9) Optional: confirm the real GEO series is reachable.
if [ "$WITH_GEO" = 1 ]; then
    $IGVF spatial-hic pull-geo --gse GSE307620 --label "${LABEL}_geo"
fi

# 10) Score recovery against the planted truth.
$PY "$HERE/make_figures.py"

echo
echo "Score with: $PY Benchmarks/concordance.py --benchmark $LABEL"
