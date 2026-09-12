#!/usr/bin/env bash
# Wang 2026 Spatial-ATAC-Hi-C — reproduction on the paper's OWN GEO deposit.
#
# Companion to Benchmarks/wang2026_spatial_atac_hic, which proves the
# spatial-hic math recovers planted truth on synthetic data in ~20 s with no
# network. This one is the other half: the same chain on GSE307620, scored
# against the numbers the paper published.
#
# The metadata leg needs only network. The analysis legs need one sample's
# pairs + fragments staged locally, because the series ships a 9.9 GB tar.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
LABEL=wang2026_spatial_atac_hic_geo
DATA="Benchmarks/_data/$LABEL"
PAIRS="${PAIRS:-$DATA/sample.pairs.gz}"
FRAGMENTS="${FRAGMENTS:-$DATA/fragments.tsv.gz}"

echo "== 1. inventory the GEO deposit (no download) =="
igvfagent spatial-hic pull-geo --gse GSE307620 --label "$LABEL"

if [ ! -f "$PAIRS" ]; then
  cat <<MSG

Metadata leg complete. The analysis legs need real data:

  mkdir -p $DATA
  # from GSE307620, stage one sample's pairs + matching ATAC fragments as
  #   $PAIRS
  #   $FRAGMENTS
  # then re-run this script.

Scoring now reports the GEO check and leaves the rest failed-for-want-of-data,
which is the honest outcome — it will NOT report a pass.
MSG
  exit 0
fi

echo "== 2. 50x50 pixel grid (Fig. 1a,b) =="
igvfagent spatial-hic pixel-demux  --pairs "$PAIRS" --label "$LABEL"
echo "== 3. per-pixel QC (Fig. 1e-g) =="
igvfagent spatial-hic qc           --pairs "$PAIRS" --fragments "$FRAGMENTS" --label "$LABEL"
echo "== 4. gene activity score, ATAC (Fig. 2b) =="
igvfagent spatial-hic gas          --fragments "$FRAGMENTS" --label "$LABEL"
echo "== 5. gene-associated domain score, Hi-C (Fig. 2a) =="
igvfagent spatial-hic gad          --pairs "$PAIRS" --label "$LABEL"
echo "== 6. A/B compartments at 100 kb (Fig. 3a) =="
igvfagent spatial-hic compartment  --pairs "$PAIRS" --resolution 100000  --label "$LABEL"
echo "== 7. per-pixel CNV at 5 Mb, diploid baseline (Figs. 4,5) =="
igvfagent spatial-hic cnv          --pairs "$PAIRS" --resolution 5000000 --per-pixel --label "$LABEL"

echo
echo "Score with: python3 Benchmarks/concordance.py --benchmark $LABEL"
