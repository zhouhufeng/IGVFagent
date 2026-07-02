#!/usr/bin/env bash
# Travaglini 2020 Human Lung Cell Atlas — full local reproduction via sc-analyze.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="travaglini2020_lung"
DATA="Benchmarks/_data/travaglini2020_lung"
H5AD="$DATA/lung_atlas.h5ad"
PY=".venv/bin/python"
IGVF=".venv/bin/igvfagent"

mkdir -p "$DATA"

# 1) Fetch the CELLxGENE-distributed h5ad (Krasnow Lab Human Lung Cell Atlas, 10X;
#    collection DOI 10.1038/s41586-020-2922-4). ~596 MB, public, no auth.
if [ ! -f "$H5AD" ]; then
  echo "== downloading lung atlas h5ad (~596 MB) =="
  URL="https://datasets.cellxgene.cziscience.com/f5568ea3-c249-4e4e-91f8-46abc30a5612.h5ad"
  curl -L -o "$H5AD" "$URL"
fi

# 2) CELLxGENE -> Scanpy prep (raw counts -> X, Ensembl -> gene symbols).
$PY "$HERE/prep_input.py"

# 3) The reproduction: IGVFagent's Scanpy pipeline (QC -> HVG -> PCA -> UMAP ->
#    Leiden -> Wilcoxon markers). Author cell_type labels are carried through for
#    scoring; canonical lung markers are overlaid on the UMAP.
$IGVF sc-analyze pipeline \
  --input "$DATA/lung_for_scanalyze.h5ad" \
  --label "$LABEL" \
  --resolution 1.5 --n-hvg 2000 --n-pcs 50 --skip-tsne \
  --sample-col cell_type \
  --highlight-genes SFTPC,AGER,SCGB1A1,FOXJ1,PECAM1,MARCO,CD3E,DCN

# 4) Score concordance (ARI / AMI / homogeneity vs 46 author cell types) + figures.
$PY "$HERE/make_figures.py"

echo ""
echo "== Travaglini 2020 lung benchmark complete =="
echo "Score with: $PY Benchmarks/concordance.py --benchmark $LABEL"
