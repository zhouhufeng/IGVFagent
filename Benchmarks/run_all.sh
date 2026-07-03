#!/usr/bin/env bash
# Drive the IGVFagent reproducibility benchmark suite.
#
# Usage:
#   bash Benchmarks/run_all.sh                 # online + local (will skip local if inputs absent)
#   bash Benchmarks/run_all.sh --online-only   # only papers needing no local downloads
#   bash Benchmarks/run_all.sh --quick         # the 3 cleanest only
#
# Each per-paper run.sh writes outputs under Docs/<skill>/<ts>_<label>/.
# Final concordance is scored via concordance.py after this script finishes.

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-all}"

# Three modes: full / online-only / quick
case "$MODE" in
    --online-only|-o)
        PAPERS=(
            matreyek2018_pten_vampseq    # verified smoke-test
            waters2024_bap1
            buckley2024_vhl
            joung2025_tf_perturbseq
            weinstock2024_cd4_crispr
            zou2024_chipatlas_gata1
            yao2024_encode4_crispri      # metadata-only step is online
            mitra2024_scarlink           # multiome retrieve is online
            martyn2025_variant_flowfish  # pull-portal step is online
        )
        ;;
    --quick|-q)
        PAPERS=(
            matreyek2018_pten_vampseq    # always works
            martyn2025_variant_flowfish
            mitra2024_scarlink
        )
        ;;
    *)
        PAPERS=(
            matreyek2018_pten_vampseq
            waters2024_bap1
            buckley2024_vhl
            yao2024_encode4_crispri
            mitra2024_scarlink
            agarwal2025_lentimpra
            joung2025_tf_perturbseq
            zheng2024_invivo_perturbseq
            deng2024_cortex_mpra
            weinstock2024_cd4_crispr
            martyn2025_variant_flowfish
            zou2024_chipatlas_gata1
            travaglini2020_lung            # full local repro: sc-analyze (CELLxGENE h5ad)
            trevino2021_cortex_multiome    # full local repro: multiome peak2gene (GSE162170)
            demultiplex2_stoeckius         # full local repro: multiseq (deMULTIplex2 data)
            rosenberg2018_splitseq         # full local repro: splitseq (GSE110823 MATLAB DGE)
        )
        ;;
esac

echo "=== IGVFagent reproducibility benchmark suite ==="
echo "Mode:     $MODE"
echo "Papers:   ${#PAPERS[@]}"
echo "Started:  $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo ""

PASS=0
FAIL=0
SKIP=0
FAILED_PAPERS=()
SKIPPED_PAPERS=()

for paper in "${PAPERS[@]}"; do
    run_script="$HERE/$paper/run.sh"
    if [ ! -x "$run_script" ] && [ ! -r "$run_script" ]; then
        echo ">> [$paper] missing run.sh — SKIP"
        SKIP=$((SKIP + 1))
        SKIPPED_PAPERS+=("$paper")
        continue
    fi
    echo ">> [$paper] running …"
    if bash "$run_script"; then
        echo ">> [$paper] OK"
        PASS=$((PASS + 1))
    else
        rc=$?
        if [ "$rc" -eq 77 ]; then
            # Conventional exit code for "skipped — missing local input"
            echo ">> [$paper] SKIPPED (missing local input)"
            SKIP=$((SKIP + 1))
            SKIPPED_PAPERS+=("$paper")
        else
            echo ">> [$paper] FAILED rc=$rc"
            FAIL=$((FAIL + 1))
            FAILED_PAPERS+=("$paper")
        fi
    fi
    echo ""
done

echo "=== Done ==="
echo "  passed:  $PASS"
echo "  failed:  $FAIL  ${FAILED_PAPERS[*]:-}"
echo "  skipped: $SKIP  ${SKIPPED_PAPERS[*]:-}"
echo ""
echo "Now score concordance:"
echo "  .venv/bin/python Benchmarks/concordance.py --all"
exit 0
