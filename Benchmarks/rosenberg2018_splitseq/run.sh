#!/usr/bin/env bash
# Rosenberg 2018 SPLiT-seq (Science 360:176) — claim-by-claim reproduction from
# the authors' GEO deposit (GSE110823). See README.md for the coverage table.
#
#   Fig. 1B-D  species-mixing libraries GSM3017262-65 -> splitseq-barnyard-qc
#              (port of Alex-Rosenberg/split-seq-pipeline split_seq/analysis.py)
#   Figs. 2,4,5 full 156,049-nucleus atlas GSM3017261 -> de novo re-clustering
#   Fig. 3     Allen Developing Mouse Brain Atlas ISH composites per cluster
#
# The atlas steps need ~16 GB RAM and ~10 min; run them inside an allocation
# (sbatch/salloc), not on a login node.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="rosenberg2018_splitseq"
D="Data/rosenberg2018"
mkdir -p "$D"

# The repo-local .venv may have been built on another machine; fall back to
# whatever igvfagent is on $PATH, and use the python next to it.
IGVF=".venv/bin/igvfagent"
"$IGVF" --help >/dev/null 2>&1 || IGVF="$(command -v igvfagent)"
PY="$(dirname "$IGVF")/python"
TS="$(date +%Y%m%d_%H%M%S)"

# 1) GEO GSE110823 per-sample MATLAB DGEs (public).
GEO="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3017nnn"
for s in GSM3017261/suppl/GSM3017261_150000_CNS_nuclei \
         GSM3017262/suppl/GSM3017262_same_day_cells_nuclei_3000_UBCs \
         GSM3017263/suppl/GSM3017263_same_day_cells_nuclei_300_UBCs \
         GSM3017264/suppl/GSM3017264_frozen_preserved_cells_nuclei_1000_UBCs \
         GSM3017265/suppl/GSM3017265_frozen_preserved_cells_nuclei_200_UBCs; do
  f="$D/$(basename "$s").mat.gz"
  [ -f "$f" ] || curl -sfL -o "$f" "$GEO/$s.mat.gz"
done

# 2) Fig. 1B-D: runs the registered port on cells and nuclei + Fig. 1D correlation.
"$PY" "$HERE/fig1_metrics.py" "Docs/SPLiTseq/${TS}_${LABEL}_fig1"

# 3) Figs. 2, 4, 5: atlas -> AnnData with the deposited labels, then re-cluster.
[ -f "$D/rosenberg_cns_full.h5ad" ] || "$PY" "$HERE/prep_atlas.py"
"$PY" "$HERE/atlas_recluster.py" "Docs/SPLiTseq/${TS}_${LABEL}_atlas"

# 4) Fig. 3: top-5 genes per deposited cluster -> Allen DMBA composites (network).
"$PY" "$HERE/fig3_allen_regions.py" de    "Docs/SPLiTseq/${TS}_${LABEL}_fig3"
"$PY" "$HERE/fig3_allen_regions.py" allen "Docs/SPLiTseq/${TS}_${LABEL}_fig3"

echo ""
echo "== Rosenberg 2018 SPLiT-seq benchmark complete =="
echo "Score with: igvfagent bench score --paper-id $LABEL"
