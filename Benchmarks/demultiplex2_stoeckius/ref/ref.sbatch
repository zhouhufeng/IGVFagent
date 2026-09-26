#!/bin/bash
#SBATCH -J demux2_ref
#SBATCH -p hsph,shared
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o Data/zhu2024/work/logs/%x_%j.out
# usage: sbatch ref.sbatch <tags.csv> <out_prefix> [max.cell.fit] [max.iter]
export PATH=/n/sw/helmod-rocky8/apps/Core/R/4.4.1-fasrc01/bin:$PATH
Rscript Benchmarks/demultiplex2_stoeckius/ref/run_demux2_ref.R Data/zhu2024/pkg_repo "$@"
