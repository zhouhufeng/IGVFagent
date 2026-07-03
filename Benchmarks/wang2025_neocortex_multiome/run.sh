#!/usr/bin/env bash
# Wang 2025 developing human neocortex multiome — full local reproduction (sc-analyze).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"; cd "$ROOT"
LABEL="wang2025_neocortex_multiome"; D="Benchmarks/_data/$LABEL"; PY=".venv/bin/python"; IGVF=".venv/bin/igvfagent"
H="$D/wang_multiome.h5ad"; mkdir -p "$D"
# 1) Fetch the CELLxGENE multiome h5ad (collection ad2149fc; ~2.6 GB, resumable).
URL="https://datasets.cellxgene.cziscience.com/a4310202-4dc8-4e1b-a96d-d9675f5b14d1.h5ad"
if [ ! -f "$H" ] || [ "$(stat -c%s "$H" 2>/dev/null)" -lt 2780000000 ]; then
  for i in $(seq 1 60); do curl -sL -C - --retry 5 -o "$H" "$URL"
    [ "$(stat -c%s "$H" 2>/dev/null)" -ge 2780000000 ] && break; sleep 4; done
fi
# 2) Memory-safe subsample -> raw counts + gene symbols + author labels.
$PY "$HERE/prep_input.py"
# 3) The reproduction: Scanpy QC -> HVG -> PCA -> UMAP -> Leiden -> markers.
$IGVF sc-analyze pipeline --input "$D/wang_rna_20k.h5ad" --label "$LABEL" \
  --resolution 1.5 --n-hvg 2000 --n-pcs 50 --skip-tsne --sample-col cell_type \
  --highlight-genes PAX6,EOMES,NEUROD2,SATB2,BCL11B,GAD1,AQP4,OLIG1
# 4) Score cluster concordance vs the 29 author cell types + figure.
$PY "$HERE/make_figures.py"
echo "Score with: $PY Benchmarks/concordance.py --benchmark $LABEL"
