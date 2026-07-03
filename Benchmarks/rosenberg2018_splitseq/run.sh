#!/usr/bin/env bash
# Rosenberg 2018 SPLiT-seq CNS atlas — full local reproduction via the splitseq skill.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="rosenberg2018_splitseq"
D="Benchmarks/_data/$LABEL"
PY=".venv/bin/python"
IGVF=".venv/bin/igvfagent"
mkdir -p "$D"

# 1) Fetch the GSE110823 RAW tar (public; contains per-sample MATLAB DGE files).
if [ ! -f "$D/GSE110823_RAW.tar" ]; then
  curl -sL -C - -o "$D/GSE110823_RAW.tar" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110823/suppl/GSE110823_RAW.tar"
fi

# scipy.io reads MATLAB v5; scikit-misc powers the seurat_v3 HVG flavor.
$PY -c "import skmisc" 2>/dev/null || $PY -m pip install -q scikit-misc

# 2) MATLAB DGE -> AnnData; subsample the 156,049-nucleus atlas to a memory-safe 12k.
$PY "$HERE/prep_input.py"

# 3) The reproduction: SPLiT-seq analyze — QC -> normalize -> Harmony -> UMAP ->
#    Leiden -> auto-annotate against the bundled mouse-brain marker panel.
$IGVF splitseq analyze \
  --input "$D/rosenberg_cns_12k.h5ad" \
  --tissue brain --resolution 1.0 --min-cells 10 --batch-key sample_id \
  --label "$LABEL"

# 4) Score CNS cell-type recovery + figure.
$PY "$HERE/make_figures.py"

echo ""
echo "== Rosenberg 2018 SPLiT-seq benchmark complete =="
echo "Score with: $PY Benchmarks/concordance.py --benchmark $LABEL"
