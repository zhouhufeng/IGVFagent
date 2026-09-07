#!/usr/bin/env bash
# Tabula Sapiens 2.0 (Cell 2026) — reproduction via the `tabula` skill.
#
# Two tiers, because the inputs differ by three orders of magnitude:
#
#   metadata tier (default)  41 MB from GEO GSE306755. Reproduces every
#                            dataset-level number in the paper — cells,
#                            donors, tissues, cell types, populations,
#                            droplet/FACS split, demographics — in under
#                            a minute, with no atlas download.
#   atlas tier (--full)      + 57 GB of per-tissue h5ads from figshare.
#                            Adds TF specificity (Figure 2), TF GO
#                            enrichment (Figure 3) and senescence
#                            (Figure 4).
#
# Raw reads are NOT part of either tier. The AWS bucket holds 103 TB but
# every GetObject returns 403: a signed data transfer agreement is
# required for donor genetic privacy. `tabula s3-manifest` documents it.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"; cd "$ROOT"

LABEL="quake2026_tabula_sapiens"
PY="${PY:-$ROOT/.venv/bin/python}"
IGVF="${IGVF:-$PY -m igvfagent.cli}"
FULL=0
[ "${1:-}" = "--full" ] && FULL=1

[ -x "$PY" ] || { echo "No interpreter at $PY. Create one: python -m venv .venv"; exit 2; }

# 1) The Human TF database — 1,639 curated TFs, the paper's own list.
$IGVF humantfs build-db --force >/dev/null 2>&1 || true
$IGVF humantfs info | head -12

# 2) Cell metadata (41 MB) — everything Figure 1 needs.
$IGVF tabula pull-geo

# 3) Local state against the paper's numbers.
$IGVF tabula status

# 4) Figure 1 — donors, tissues, composition.
$IGVF tabula overview --label "${LABEL}_overview"

if [ "$FULL" = 1 ]; then
    # 5) The 57 GB atlas. Resumable; skips complete files.
    $IGVF tabula pull-figshare

    # 6) Figure 2 — TF specificity. Two steps so the expensive
    #    mean-expression matrix is cached and reusable.
    $IGVF tabula tf-matrix --label "${LABEL}_tf_matrix"
    M="$(ls -t Docs/TabulaSapiens/2*_${LABEL}_tf_matrix/mean_expression.npz | head -1)"
    $IGVF tabula tf-specificity --matrix "$M" --label "${LABEL}_tf_specificity"

    # 7) Figure 3 — GO enrichment of the non-specific TFs.
    T="$(ls -t Docs/TabulaSapiens/2*_${LABEL}_tf_specificity/tf_tau.tsv | head -1)"
    $IGVF tabula tf-enrichment --tau-table "$T" --label "${LABEL}_tf_enrichment"

    # 8) Figure 4 — senescent-cell burden.
    $IGVF tabula senescence --label "${LABEL}_senescence"
fi

# 9) Score everything that ran against the paper.
$PY "$HERE/make_figures.py"

echo
echo "Score with: $PY Benchmarks/concordance.py --benchmark $LABEL"
[ "$FULL" = 0 ] && echo "(metadata tier only — re-run with --full for Figures 2-4)"
exit 0
