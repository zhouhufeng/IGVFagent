#!/usr/bin/env bash
# Run the UPSTREAM pipeline (EngreitzLab/CRISPR_comparison) on the same inputs, as the gold
# standard for `igvfagent sce2g benchmark`. This is how results/upstream/ was produced on
# 2026-09-22 (Ubuntu VM, Docker, condaforge/miniforge3, CRISPR_comparison @ 5058742).
#
#   bash Benchmarks/sce2g_crowdsourced_features_k562/run_upstream.sh /path/to/workdir
#
# Needs: docker, ~8 GB free, the inputs listed below in <workdir>/CRISPR_comparison/inputs/
# (the feature tables and rE2G predictions come from Synapse; see run.sh for the downloads).
# Runtime: ~10 min to build the R environment, then ~25 min for the five comparisons.
set -euo pipefail
W="${1:?workdir}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$W" && cd "$W"
[ -d CRISPR_comparison ] || git clone -q https://github.com/EngreitzLab/CRISPR_comparison.git
cd CRISPR_comparison
mkdir -p inputs config
cp "$HERE"/upstream_run/config.yml config/config.yml
cp "$HERE"/upstream_run/pred_config_*.txt "$HERE"/upstream_run/celltype_map_sce2g.txt inputs/
cp "$HERE"/upstream_run/run_in_docker.sh .
# Inputs (same files run.sh uses):
#   inputs/EPCrisprBenchmark_ensemble_data_GRCh38.intGENCODEv43.tsv.gz   (scE2G repo)
#   inputs/IGVFFI1706PNVV.tsv.gz                                          (IGVF portal, scE2G K562)
#   inputs/ENCODE-rE2G_*CrossValidated_K562*.predictions.tsv.gz           (Synapse syn53019593/5)
#   inputs/uniformly_processed_K562_candidate_e2g_pairs_pinloop.tsv.gz    (Synapse syn73717888)
#   inputs/K562_ERR9847049_Multiome_K562_Signac.tsv.gz, ..._SCENT.tsv.gz  (Synapse syn73717888)
for f in EPCrisprBenchmark_ensemble_data_GRCh38.intGENCODEv43.tsv.gz IGVFFI1706PNVV.tsv.gz; do
  [ -f "inputs/$f" ] || echo "missing inputs/$f"
done
# The R environment is created once from workflow/envs/r_crispr_comparison.yml with strict
# channel priority; snakemake then runs with it on PATH. Snakemake's own --use-conda failed
# twice here (env creation under non-strict priority, then `conda env export` on a path env),
# so the conda: directives are dropped from the rules for this run. Nothing else is changed.
docker run --rm -v "$PWD:/work" -w /work condaforge/miniforge3:latest bash -lc \
  "conda config --set channel_priority strict && mamba env create -y --file workflow/envs/r_crispr_comparison.yml --prefix /work/.snakemake/conda/manual_r_env"
sed -i.bak '/conda: "\.\.\/envs\/r_crispr_comparison.yml"/d' workflow/rules/crispr_comparison.smk
docker run --rm --name crispr_comparison -v "$PWD:/work" -w /work condaforge/miniforge3:latest bash /work/run_in_docker.sh
ls results/*/performance_summary.txt
