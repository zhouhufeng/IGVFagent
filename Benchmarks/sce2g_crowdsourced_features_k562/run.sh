#!/usr/bin/env bash
# scE2G crowdsourced-feature benchmark — every feature column of the 16 IGVF
# K562 E2G feature tables on Synapse, scored one at a time against the K562
# CRISPR element-gene ground truth, next to the released scE2G K562 model.
#
#   bash Benchmarks/sce2g_crowdsourced_features_k562/run.sh
#   .venv/bin/python Benchmarks/concordance.py --benchmark sce2g_crowdsourced_features_k562
#
# Needs SYNAPSE_AUTH_TOKEN (or Docs/Secret/SYNAPSE_AUTH_TOKEN.txt) for the
# Synapse folder: it is access-controlled. About 6 GB is downloaded once; the
# two EPCOT tables are 2 GB each. Total runtime after download ~25 min on a
# laptop (the EPCOT batch dominates).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
PY="${PYTHON:-.venv/bin/python}"
IGVF="${IGVFAGENT:-.venv/bin/igvfagent}"
export PYTHONPATH="$ROOT/Scripts${PYTHONPATH:+:$PYTHONPATH}"

FEAT="Data/scE2G/Synapse/syn73717888"
RES="Data/scE2G/resources"
CRISPR="$RES/EPCrisprBenchmark_ensemble_data_GRCh38.intGENCODEv43.tsv.gz"
SCE2G_DIR="Data/scE2G/scE2G_K562"
SCE2G="$SCE2G_DIR/IGVFFI1706PNVV.tsv.gz"
mkdir -p "$FEAT" "$RES" "$SCE2G_DIR" "$HERE/results"

# 0) Synapse token: env var wins; else the gitignored secret file.
if [ -z "${SYNAPSE_AUTH_TOKEN:-}" ] && [ -f Docs/Secret/SYNAPSE_AUTH_TOKEN.txt ]; then
  SYNAPSE_AUTH_TOKEN="$(tr -d '[:space:]' < Docs/Secret/SYNAPSE_AUTH_TOKEN.txt)"
  export SYNAPSE_AUTH_TOKEN
fi

# 1) The 16 crowdsourced K562 feature tables (Synapse folder syn73717888).
if [ "$(ls "$FEAT"/*.gz 2>/dev/null | wc -l | tr -d ' ')" -lt 16 ]; then
  [ -n "${SYNAPSE_AUTH_TOKEN:-}" ] || { echo "SYNAPSE_AUTH_TOKEN is required to download syn73717888"; exit 77; }
  "$IGVF" synapse walk --syn syn73717888 --max-depth 2 --max-children 200 >/dev/null
  MANIFEST="$(ls -t Data/Manifests/Synapse/*syn73717888_walk.csv | head -1)"
  for id in $(awk -F, 'NR>1 && $4=="file"{print $2}' "$MANIFEST"); do
    "$IGVF" synapse download --syn "$id" --out-dir "$FEAT" | grep -i "downloaded" || true
  done
fi

# 2) CRISPR ground truth, as shipped inside the scE2G repository.
[ -f "$CRISPR" ] || curl -sL -o "$CRISPR" \
  "https://raw.githubusercontent.com/EngreitzLab/scE2G/main/resources/EPCrisprBenchmark_ensemble_data_GRCh38.intGENCODEv43.tsv.gz"

# 3) Released scE2G v1.2 K562 predictions (IGVF PredictionSet IGVFDS5428HHMB), the reference.
[ -f "$SCE2G" ] || curl -sL -o "$SCE2G" \
  "https://api.data.igvf.org/tabular-files/IGVFFI1706PNVV/@@download/IGVFFI1706PNVV.tsv.gz"

# 4) What is in the folder.
"$PY" Scripts/sce2g_workbench_skill.py inventory --dir "$FEAT" --max-rows 3000000 --label k562_features_inventory

# 5) Every feature column of every table, in four batches (EPCOT tables are 2 GB each).
B="$PY Scripts/sce2g_workbench_skill.py benchmark --all-features --crispr $CRISPR --bootstrap 200"
$B --label k562_features_A \
  --predictions ChromHMM="$FEAT/K562_E2G_table.with.ChromHMM.state.overlap.features_CrhmmBool.txt.gz" \
  --predictions ArchR="$FEAT/K562_ERR9847049_Multiome_K562_ArchR.tsv.gz" \
  --predictions Cicero="$FEAT/K562_ERR9847049_Multiome_K562_Cicero.tsv.gz" \
  --predictions SCENT="$FEAT/K562_ERR9847049_Multiome_K562_SCENT.tsv.gz" \
  --predictions Signac="$FEAT/K562_ERR9847049_Multiome_K562_Signac.tsv.gz"
$B --label k562_features_B \
  --predictions ENCODEstats="$FEAT/uniformly_processed_K562_candidate_e2g_pairs_gene_ENCODE_stats.tsv.gz" \
  --predictions phastCons="$FEAT/uniformly_processed_K562_candidate_e2g_pairs_phastcons.tsv.gz" \
  --predictions phyloP="$FEAT/uniformly_processed_K562_candidate_e2g_pairs_phylop.tsv.gz" \
  --predictions Pinloop="$FEAT/uniformly_processed_K562_candidate_e2g_pairs_pinloop.tsv.gz" \
  --predictions pLI_LOEUF="$FEAT/uniformly_processed_K562_candidate_e2g_pairs_pli_loeuf.tsv.gz"
$B --label k562_features_C \
  --predictions sHet="$FEAT/uniformly_processed_K562_candidate_e2g_pairs_shet.tsv.gz" \
  --predictions TFgene="$FEAT/uniformly_processed_K562_candidate_e2g_pairs_tf_genes_marked.tsv.gz" \
  --predictions Alu="$FEAT/uniformly_processed_K562.Alu_features.tsv.gz" \
  --predictions Motif="$FEAT/uniformly_processed_K562.MotifDensityFeatures.tsv.gz"
$B --label k562_features_D_epcot \
  --predictions EPCOT="$FEAT/K562_ERR9847049_Multiome_K562_EPCOT.tsv.gz" \
  --predictions EPCOT_woK562="$FEAT/K562_ERR9847049_Multiome_K562_EPCOT_wok562.tsv.gz"

# 6) The released scE2G model's four scores as the reference predictors.
"$PY" Scripts/sce2g_workbench_skill.py benchmark --crispr "$CRISPR" --bootstrap 200 \
  --label k562_sce2g_reference --pred-config "$HERE/pred_config_reference.txt" \
  --predictions scE2G="$SCE2G" --predictions scE2G_ignoreTPM="$SCE2G" \
  --predictions ABC="$SCE2G" --predictions ARC_E2G="$SCE2G"

# 6b) Published anchor: ENCODE-rE2G cross-validated K562 predictions (Synapse), which ship one
#     prediction per benchmark pair, so no overlap step is involved and the published AUPRC
#     (0.634 base / 0.543 precision at 70% recall) must come back exactly.
RE2G="Data/E2G_benchmark/rE2G_K562"
mkdir -p "$RE2G"
if [ -n "${SYNAPSE_AUTH_TOKEN:-}" ]; then
  for s in syn53019593 syn53019595; do
    ls "$RE2G"/*.gz >/dev/null 2>&1 && [ "$(ls "$RE2G"/*.gz | wc -l | tr -d ' ')" -ge 2 ] || \
      "$IGVF" synapse download --syn "$s" --out-dir "$RE2G" | grep -i downloaded || true
  done
fi
RE2G_BASE="$RE2G/ENCODE-rE2G_CrossValidated_K562_EPCrisprBenchmark_ensemble_data_GRCh38.predictions.tsv.gz.tsv.gz"
RE2G_EXT="$RE2G/ENCODE-rE2G_extended_CrossValidated_K562_EPCrisprBenchmark_ensemble_data_GRCh38.predictions.tsv.gz"
ANCHOR=()
if [ -f "$RE2G_BASE" ] && [ -f "$RE2G_EXT" ]; then
  "$PY" Scripts/sce2g_workbench_skill.py benchmark --crispr "$CRISPR" --bootstrap 200 \
    --label k562_re2g_published_anchor --pred-config "$HERE/pred_config_re2g_anchor.txt" \
    --predictions rE2G_base="$RE2G_BASE" --predictions rE2G_ext="$RE2G_EXT" --predictions rE2G_ext_noEP300="$RE2G_EXT"
  ANCHOR=(k562_re2g_published_anchor)
else
  echo "rE2G anchor skipped (needs SYNAPSE_AUTH_TOKEN for syn53019593 / syn53019595)"
fi

# 7) One ranked table, figure, report and summary.json (what concordance.py scores).
runs=()
for lbl in k562_features_A k562_features_B k562_features_C k562_features_D_epcot k562_sce2g_reference "${ANCHOR[@]}"; do
  runs+=(--runs "$(ls -d Docs/scE2G/*_"$lbl" | tail -1)")
done
"$PY" Scripts/sce2g_workbench_skill.py merge "${runs[@]}" --top 30 --label k562_crowdsourced_feature_benchmark

# 8) Committed figures, and the side-by-side with the upstream pipeline's own run
#    (results/upstream/, produced with run_upstream.sh).
"$PY" "$HERE/make_figures.py"
"$PY" "$HERE/compare_upstream.py" || true
echo "scE2G crowdsourced-feature benchmark complete."
