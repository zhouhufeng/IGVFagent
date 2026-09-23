#!/usr/bin/env bash
set -uo pipefail
cd /work
source /opt/conda/etc/profile.d/conda.sh
conda config --set channel_priority strict
conda env list | grep -q "^smk " || mamba create -y -q -n smk -c conda-forge -c bioconda "snakemake-minimal>=7.32,<8" >/tmp/smk_install.log 2>&1 || { tail -20 /tmp/smk_install.log; exit 1; }
conda activate smk
# The R environment was built once from workflow/envs/r_crispr_comparison.yml (strict channel
# priority); snakemake runs the R scripts with it on PATH instead of managing conda itself.
export PATH=/work/.snakemake/conda/manual_r_env/bin:$PATH
Rscript --version
T=""
for c in re2g_anchor sce2g_igvf sce2g_igvf_notssfilter features_igvf features_igvf_notssfilter; do T="$T results/$c/${c}_crispr_comparison.html"; done
snakemake -j 2 --rerun-incomplete --keep-going $T
echo "SNAKEMAKE_EXIT=$?"
