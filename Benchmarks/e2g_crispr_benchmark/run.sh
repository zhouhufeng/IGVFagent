#!/usr/bin/env bash
# E2G CRISPR benchmark — ENCODE-rE2G & scE2G vs CRISPR ground truth.
# Fetches data with IGVFagent, scores predictions, regenerates figures.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
DATA="Data/E2G_benchmark"
GT="$DATA/crispr_ground_truth/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz"
mkdir -p "$DATA/scE2G_K562" "$DATA/rE2G_K562" "$DATA/crispr_ground_truth" "$HERE/results"

# 1) scE2G K562 scored table (IGVF portal) — 325 MB
[ -f "$DATA/scE2G_K562/IGVFFI1706PNVV.tsv.gz" ] || \
  curl -sL "https://api.data.igvf.org/tabular-files/IGVFFI1706PNVV/@@download/IGVFFI1706PNVV.tsv.gz" \
       -o "$DATA/scE2G_K562/IGVFFI1706PNVV.tsv.gz"

# 2) CRISPR ground truth (EngreitzLab/CRISPR_comparison)
[ -f "$GT" ] || curl -sL \
  "https://raw.githubusercontent.com/EngreitzLab/CRISPR_comparison/main/resources/crispr_data/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz" \
  -o "$GT"

# 3) ENCODE-rE2G cross-validated K562 predictions (Synapse; needs SYNAPSE_AUTH_TOKEN)
#    igvfagent synapse download --syn syn53019593 --out-dir "$DATA/rE2G_K562"
#    igvfagent synapse download --syn syn53019595 --out-dir "$DATA/rE2G_K562"

# 4) Score scE2G table (scE2G / ABC / ARC-E2G)
python3 Scripts/e2g_benchmark_eval.py --ground-truth "$GT" \
  --predictions "$DATA/scE2G_K562/IGVFFI1706PNVV.tsv.gz" \
  --score-cols "Score,ABC.Score,ARC.E2G.Score,Score.ignoreTPM" \
  --out "$HERE/results/k562_scE2G_vs_crispr.json"

# 5) Score ENCODE-rE2G base + extended
python3 Scripts/e2g_benchmark_eval.py --ground-truth "$GT" \
  --predictions "$DATA/rE2G_K562/ENCODE-rE2G_CrossValidated_K562_EPCrisprBenchmark_ensemble_data_GRCh38.predictions.tsv.gz.tsv.gz" \
  --gene-col TargetGene --chr-col chr --start-col start --end-col end \
  --score-cols "Full.Score" --out "$HERE/results/k562_rE2G_base_vs_crispr.json"
python3 Scripts/e2g_benchmark_eval.py --ground-truth "$GT" \
  --predictions "$DATA/rE2G_K562/ENCODE-rE2G_extended_CrossValidated_K562_EPCrisprBenchmark_ensemble_data_GRCh38.predictions.tsv.gz" \
  --gene-col TargetGene --chr-col chr --start-col start --end-col end \
  --score-cols "FullModel.Score,FullModel_minus_EP300.Score" \
  --out "$HERE/results/k562_rE2G_extended_vs_crispr.json"

# 6) Figures
python3 "$HERE/make_figures.py"
echo "E2G CRISPR benchmark complete."
