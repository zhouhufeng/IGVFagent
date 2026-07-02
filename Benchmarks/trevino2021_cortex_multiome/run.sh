#!/usr/bin/env bash
# Trevino 2021 developing-cortex multiome — full local reproduction of the
# peak->gene cis-regulatory linkage via `igvfagent multiome peak2gene`.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="trevino2021_cortex_multiome"
D="Benchmarks/_data/$LABEL"
B="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162170/suppl"
PY=".venv/bin/python"
IGVF=".venv/bin/igvfagent"
mkdir -p "$D"

# 1) Fetch the paired-multiome subset of GSE162170 (public, no auth).
for f in multiome_rna_counts.tsv.gz multiome_atac_counts.tsv.gz \
         multiome_atac_consensus_peaks.txt.gz multiome_cell_metadata.txt.gz \
         multiome_cluster_names.txt.gz; do
  if [ ! -f "$D/GSE162170_$f" ] || ! gzip -t "$D/GSE162170_$f" 2>/dev/null; then
    echo "== downloading $f =="
    curl -sL -o "$D/GSE162170_$f" "$B/GSE162170_$f"
  fi
done

# 2) Build panel-restricted RNA + ATAC AnnData + TSS BED (Ensembl GRCh38 TSS).
#    Restricts the 467k-peak / 34k-gene matrices to a 26-gene cortical-lineage
#    panel + peaks within 500 kb of their TSS so peak2gene runs on a laptop.
$PY "$HERE/build_multiome_h5ads.py"

# 3) The reproduction: per-peak Spearman correlation with each panel gene's RNA.
$IGVF multiome peak2gene \
  --rna-h5ad "$D/multiome_rna_panel.h5ad" \
  --atac-h5ad "$D/multiome_atac_panel.h5ad" \
  --tss-bed "$D/panel_tss.bed" \
  --window 500000 --method spearman \
  --label "$LABEL"

# 4) Score cis-linkage concordance + figures.
$PY "$HERE/make_figures.py"

echo ""
echo "== Trevino 2021 multiome benchmark complete =="
echo "Score with: $PY Benchmarks/concordance.py --benchmark $LABEL"
