#!/usr/bin/env bash
# Reproduce Ryu et al. 2024 (crispr-bean) on the paper's own deposited data.
#
#   bash Benchmarks/ryu2024_crispr_bean/run.sh
#
# Downloads ~0.9 GB from Zenodo on first run and caches it. The measurement
# steps need `anndata`, which lives in the BEAN venv rather than the app venv
# -- see Deploy/install-bean.sh. BEAN_PY overrides the interpreter.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"

SKILL="${SKILL:-igvfagent bean-benchmark}"
# anndata is not in the app venv; the BEAN venv has it.
BEAN_PY="${BEAN_PY:-/workspace/opt/bean-venv/bin/python}"
MOD="${MOD:-Scripts/bean_benchmark_skill.py}"

run() {                       # prefer the BEAN venv, fall back to the skill
  if [[ -x "$BEAN_PY" ]]; then "$BEAN_PY" "$MOD" "$@"; else $SKILL "$@"; fi
}

echo "== 1. fetch the paper's deposited screen objects (Zenodo 10139794)"
run fetch all

echo "== 2. describe them BEFORE scoring anything"
# A replication that starts by computing a metric cannot tell "the method
# disagrees" from "this is not the data the paper described".
run describe ldlvar
run describe ldlrcds

echo "== 3. measure what the deposits alone support (no model fit)"
run measure all

echo "== 4. report, claim by claim, against the published values"
run report
