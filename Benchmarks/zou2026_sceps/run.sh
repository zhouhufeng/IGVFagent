#!/usr/bin/env bash
# scEPS benchmark — reproduce the SEA-AD microglia AD-association (Cognitive Status).
# Requires the scEPS conda env (scanpy/anndata/statsmodels) + MAGMA scores
# (see Data/MAGMA/run_magma.sh).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
D=Data/scEPS_benchmark/SEAAD
SK=Scripts/sceps_skill.py

# 0) inputs: SEA-AD Microglia h5ad (CELLxGENE) + AD MAGMA scores (Data/MAGMA/out/AD.magma.txt)
#    prep: subset to Dementia/No-dementia, add numeric CS, rebuild kNN on X_scVI (see benchmark notes)

# 1) estimate — per-neighborhood d-statistics (2000 anchors)
python3 $SK estimate --adata $D/SEAAD_Microglia_prepped.h5ad \
  --donor-id-col donor_id --pheno CS --gene-id-col feature_name \
  --gene-list Data/MAGMA/out/AD.magma.txt --auto-gene-selection \
  --scale-pheno --scale-pheno-neighborhood \
  --num-bootstrap-disattenuation 20 --num-bootstrap-regression 100 \
  --stop-idx 2000 --out $D/sceps_AD_microglia_2k.txt.gz

# 2) cluster — neighborhood blocks for the block bootstrap
python3 $SK cluster --adata $D/SEAAD_Microglia_prepped.h5ad \
  --donor-id-col donor_id --neighbors-use-rep X_scVI --num-kmeans-cluster 50 \
  --out $D/sceps_clusters.txt.gz

# 3) aggregate — cell-type-level mean + block-bootstrap significance
python3 $SK aggregate --sceps-result "$D/sceps_AD_microglia_2k.txt.gz" \
  --neighborhood-clusters $D/sceps_clusters.txt.gz --num-bootstrap 1000 \
  --out $HERE/results/sceps_AD_microglia_2k.aggregated.txt

echo "scEPS benchmark complete — see results/ + figures/"
