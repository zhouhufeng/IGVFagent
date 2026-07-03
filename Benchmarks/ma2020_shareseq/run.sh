#!/usr/bin/env bash
# Ma 2020 SHARE-seq skin — full local reproduction via the share skill + _scload.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"; cd "$ROOT"
LABEL="ma2020_shareseq"; D="Benchmarks/_data/$LABEL"; PY=".venv/bin/python"; IGVF=".venv/bin/igvfagent"
TAR="$D/GSE140203_RAW.tar"
mkdir -p "$D"
# 1) Fetch the GSE140203 RAW tar (7.5 GB; resilient resume).
if [ ! -f "$TAR" ] || ! tar -tf "$TAR" >/dev/null 2>&1; then
  for i in $(seq 1 100); do
    curl -sL -C - --retry 5 -o "$TAR" \
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140203/suppl/GSE140203_RAW.tar"
    tar -tf "$TAR" >/dev/null 2>&1 && break; sleep 4
  done
fi
# 2) Extract the skin RNA counts + author cell-type labels.
tar -xf "$TAR" -C "$D" GSM4156608_skin.late.anagen.rna.counts.txt.gz GSM4156597_skin_celltype.txt.gz
# 3) Build the RNA AnnData via the shared _scload loader (dense TSV -> sparse).
$PY "$HERE/prep_input.py"
# 4) The reproduction: SHARE-seq per-barcode RNA QC (Ma 2020 / Broad epi-SHARE-seq).
$IGVF share rna-qc --h5ad "$D/shareseq_skin_rna.h5ad" --label "$LABEL"
# 5) Score cell-type concordance (Leiden vs 23 author skin types) + figure.
$PY "$HERE/make_figures.py"
echo "Score with: $PY Benchmarks/concordance.py --benchmark $LABEL"
