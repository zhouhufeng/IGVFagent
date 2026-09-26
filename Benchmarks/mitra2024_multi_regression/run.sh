#!/usr/bin/env bash
# mitra2024_multi_regression — Mitra 2024 (SCARlink), Nat Genet
# Full local reproduction of the BMMC peak->gene cis-regulatory linkage
# on the paper's own donor/cell-type subset, via `igvfagent multiome peak2gene`.
#
# Paper (Data availability): "BMMC data were part of the NeurIPS 2021 open
# problem, and the dataset was downloaded from GEO (GSE194122). We used
# BMMC samples labeled as site1_donor1, site1_donor2, site1_donor3,
# site2_donor1, site2_donor4, site2_donor5, site3_donor10, site3_donor6,
# site3_donor7 and site4_donor9 and the cell types HSC, MK/E progenitor,
# proerythroblast, erythroblast and normoblast."
#
# NOTE: the paper's PBMC multi-ome (downloaded "from 10X Genomics", no
# accession given) and its mouse skin / cortex / pancreas / pituitary
# datasets are NOT reproduced here — only the BMMC/GSE194122 arm, which is
# the one dataset the paper names with a concrete, fetchable accession.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="mitra2024_multi_regression"
GSE="GSE194122"
RAW="Data/IGVF/GEO/Downloads/$GSE/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad"
D="Benchmarks/_data/$LABEL"
TSS_BED="Data/scE2G/resources/CollapsedGeneBounds.hg38.TSS500bp.bed"
PY=".venv/bin/python"
IGVF=".venv/bin/igvfagent"

# 1) Fetch the BMMC multiome processed h5ad (public, no auth; ~2.7 GB).
if [ ! -f "$RAW" ]; then
  echo "== downloading $GSE multiome BMMC h5ad =="
  "$IGVF" geo download --gse "$GSE" --only suppl --pattern "multiome_BMMC_processed" --max-download-gb 5
fi
if [ ! -f "$RAW" ]; then
  echo "ERROR: expected $RAW after download — check the GEO download log." >&2
  exit 1
fi

# 2) Split into RNA (GEX) + ATAC h5ads restricted to the paper's exact
#    10 donors / 5 cell types (see Data availability statement above).
if [ ! -f "$D/mitra2024_bmmc_rna.h5ad" ] || [ ! -f "$D/mitra2024_bmmc_atac.h5ad" ]; then
  "$PY" "$HERE/scripts/split_bmmc_rna_atac.py"
fi

# 3) The reproduction: per-peak Pearson correlation with each gene's RNA,
#    genome-wide (max-pairs set well above the ~920k pairs a genome-wide
#    500kb-window scan on this subset actually produces).
"$IGVF" multiome peak2gene \
  --rna-h5ad "$D/mitra2024_bmmc_rna.h5ad" \
  --atac-h5ad "$D/mitra2024_bmmc_atac.h5ad" \
  --tss-bed "$TSS_BED" \
  --window 500000 --method pearson --max-pairs 2000000 \
  --label "$LABEL"

# 4) Score concordance metrics (IGVFagent's own measured quantities —
#    see the README's Honest caveats for why these aren't a 1:1 check
#    against the paper's SCARlink-model-specific gene counts).
"$PY" "$HERE/scripts/score_metrics.py"

echo ""
echo "== mitra2024_multi_regression benchmark complete =="
echo "Score with: $PY Benchmarks/concordance.py --benchmark $LABEL"
