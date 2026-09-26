#!/bin/bash
#SBATCH -J demux2_bench_ref
#SBATCH -p hsph,shared
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o Data/zhu2024/work/logs/%x_%j.out
# usage: sbatch bench_ref.sbatch <dataset> [max.cell.fit]
export PATH=/n/sw/helmod-rocky8/apps/Core/R/4.4.1-fasrc01/bin:$PATH
W=Data/zhu2024/work/oshlack
Rscript Benchmarks/demultiplex2_stoeckius/ref/run_bench_ref.R Data/zhu2024/pkg_repo $W $W/oshlack_manifest.json "$1" $W/ref_"$1"_${2:-1000} ${2:-1000}
