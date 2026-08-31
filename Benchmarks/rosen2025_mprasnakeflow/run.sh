#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="rosen2025_mprasnakeflow"
DATA="Data/Benchmarks/$LABEL"
API="https://api.data.igvf.org"

# Inputs, straight from the IGVF portal (public, no auth):
#   IGVFFI4856GPLO  reporter experiment barcode  -> per-barcode DNA/RNA counts
#   IGVFFI3491STGH  reporter experiment          -> the AUTHORS' MPRAsnakeflow
#                                                   output for the same data
# Both belong to AnalysisSet IGVFDS1933XFAF (lentiMPRA, Ahituv lab,
# processed by the Kircher lab, i.e. the paper's own pipeline run).
BARCODES="$DATA/IGVFFI4856GPLO.tsv.gz"
PUBLISHED="$DATA/IGVFFI3491STGH.tsv.gz"

mkdir -p "$DATA"
# Complexity claims are checked on two further published barcode files:
#   IGVFFI0032GUOI  8K-neurons  (ASD library)
#   IGVFFI8345QIJJ  80K-neurons (NGN2)
for ACC in IGVFFI4856GPLO IGVFFI3491STGH IGVFFI0032GUOI IGVFFI8345QIJJ; do
    if [ ! -f "$DATA/$ACC.tsv.gz" ]; then
        echo "[$LABEL] downloading $ACC ..."
        curl -sfL -o "$DATA/$ACC.tsv.gz" \
            "$API/tabular-files/$ACC/@@download/$ACC.tsv.gz" || {
            echo "[$LABEL] download failed for $ACC — the portal may be"
            echo "  unreachable or the accession withdrawn."
            exit 77
        }
    fi
done

# Wide IGVF file -> one narrow count file per replicate + assignment
python3 "$HERE/prepare.py" "$BARCODES" "$DATA/inputs"

# The published file was produced with minDNACounts=1 / minRNACounts=1
# (no DNA pseudocount); using the tool defaults inflates dna_counts by
# exactly n_bc per oligo.
.venv/bin/igvfagent mpraflow pipeline \
    --assignment "$DATA/inputs/assignment.tsv" \
    --replicate 1="$DATA/inputs/rep1.tsv" \
    --replicate 2="$DATA/inputs/rep2.tsv" \
    --replicate 3="$DATA/inputs/rep3.tsv" \
    --min-dna-counts 1 --min-rna-counts 1 \
    --threshold 10 --label "$LABEL"

RUN=$(ls -dt Docs/MPRASnakeflow/*_"$LABEL" | head -1)
python3 "$HERE/compare.py" \
    "$RUN/master_table.all.tsv.gz" "$PUBLISHED" "$RUN/summary.json"

# Reproduce the paper's stated library-complexity numbers
echo ""
echo "[$LABEL] checking published complexity claims ..."
python3 "$HERE/claims.py" "$DATA" "$RUN/claims.json"

# Merge the claim results into the scored summary
python3 - "$RUN/summary.json" "$RUN/claims.json" <<'PYMERGE'
import json, sys
s = json.load(open(sys.argv[1])); c = json.load(open(sys.argv[2]))
s.update({k: c[k] for k in ("claims_checked", "claims_exact")})
s["claims"] = c["claims"]
json.dump(s, open(sys.argv[1], "w"), indent=2)
PYMERGE

echo ""
echo "[$LABEL] run dir: $RUN"
