#!/usr/bin/env bash
# Mirror the IGVF Knowledge Graph onto the deployment, within a disk budget.
#
#   bash Deploy/mirror-kg.sh --check          # plan only, pull nothing
#   bash Deploy/mirror-kg.sh                  # mirror what fits
#   BUDGET_GB=500 bash Deploy/mirror-kg.sh    # allow more
#
# WHY A BUDGET. The full graph is 1,928 GB across 62 collections and 11.6
# billion documents. The data volume is 2.0 TB with ~1.74 TB free, and now
# also holds container images and the FASTQ cache, so a full mirror would
# leave no headroom and fail somewhere in the middle -- after hours of
# transfer, with no way to tell a partial mirror from a complete one.
#
# Four collections are 87% of the total:
#   variants                    940.9 GB   1.87 B docs
#   variants_variants           531.2 GB   5.93 B docs
#   coding_variants_phenotypes  127.8 GB   1.10 B docs
#   coding_variants              88.1 GB   0.94 B docs
#
# Excluding those, the remaining 58 collections are ~240 GB and still cover
# 1.8 billion documents -- the genes, elements, phenotype and linkage edges
# that exploration actually touches. The four giants stay queryable live
# through the Catalog API, which already works for per-variant lookups.
set -euo pipefail

BUDGET_GB="${BUDGET_GB:-300}"
FLOOR_GB="${FLOOR_GB:-200}"        # never fill the volume below this
CONTAINER="${CONTAINER:-igvfagent-app}"
CREDS="Docs/Secret/ArangoDB-logins.txt"
CHECK_ONLY=0
[[ "${1:-}" == "--check" ]] && CHECK_ONLY=1

die()  { printf '\nFAILED: %s\n' "$1" >&2; exit 1; }
step() { printf '\n=== %s ===\n' "$1"; }

# Credentials come from the file when present, else the environment. The
# file is gitignored and so never reaches the VM by `git pull`, and copying
# it there just to read it twice would leave a third copy of a secret on
# disk for no gain.
if [[ -f "$CREDS" ]]; then
  USER_=$(grep -i '^username:' "$CREDS" | cut -d: -f2  | tr -d ' \r')
  PASS_=$(grep -i '^password:' "$CREDS" | cut -d: -f2- | tr -d ' \r')
else
  USER_="${IGVF_ARANGO_USER:-}"
  PASS_="${IGVF_ARANGO_PASSWORD:-}"
fi
[[ -n "$USER_" && -n "$PASS_" ]] || die \
  "no ArangoDB credentials: $CREDS is absent and IGVF_ARANGO_USER / \
IGVF_ARANGO_PASSWORD are unset"
DEX=(docker exec -e "IGVF_ARANGO_USER=$USER_" -e "IGVF_ARANGO_PASSWORD=$PASS_" "$CONTAINER")

step "1. Inventory the graph"
"${DEX[@]}" igvfagent kg-mirror inventory >/dev/null 2>&1 \
  || die "inventory failed — check ArangoDB credentials and egress"
INV=/mnt/igvf-data/Data/Warehouse/KG/_inventory.csv
[[ -f "$INV" ]] || die "inventory csv not written at $INV"
python3 - "$INV" "$BUDGET_GB" <<'PY' > /tmp/kg_plan.txt
import csv, sys
inv, budget = sys.argv[1], float(sys.argv[2])
rows = []
for r in csv.DictReader(open(inv)):
    r["_gb"] = float(r["bytes"]) / 1e9
    r["_n"] = int(r["documents"])
    rows.append(r)
# Smallest first: maximises how many collections land inside the budget, and
# gets useful breadth on disk early in case the run is interrupted.
rows.sort(key=lambda r: r["_gb"])
cum = 0.0
for r in rows:
    if cum + r["_gb"] > budget:
        continue
    cum += r["_gb"]
    print(f"{r['collection']}\t{r['_gb']:.3f}\t{r['_n']}")
print(f"#TOTAL\t{cum:.1f}", file=sys.stderr)
PY
PLANNED=$(wc -l < /tmp/kg_plan.txt | tr -d ' ')
PLAN_GB=$(awk -F'\t' '{s+=$2} END{printf "%.1f", s}' /tmp/kg_plan.txt)
TOTAL_COLS=$(( $(wc -l < "$INV" | tr -d ' ') - 1 ))
echo "  $TOTAL_COLS collections in the graph"
echo "  $PLANNED fit within the ${BUDGET_GB} GB budget  (~${PLAN_GB} GB)"
echo "  excluded (too large): $(awk -F'\t' 'NR>1{print $1}' "$INV" | comm -23 <(awk -F'\t' 'NR>1{print $1}' "$INV" | sort) <(cut -f1 /tmp/kg_plan.txt | sort) | tr '\n' ' ')" 2>/dev/null || true
df -h /mnt/igvf-data | awk 'NR>1{printf "  volume: %s used of %s (%s), %s free\n", $3,$2,$5,$4}'

if [[ "$CHECK_ONLY" == 1 ]]; then
  echo
  echo "Plan (smallest first):"
  head -12 /tmp/kg_plan.txt | awk -F'\t' '{printf "  %8.3f GB  %12d docs  %s\n", $2,$3,$1}'
  [[ "$PLANNED" -gt 12 ]] && echo "  … and $((PLANNED-12)) more"
  echo
  echo "--check complete. Nothing was pulled."
  exit 0
fi

step "2. Mirror, newest state resumed automatically"
LOG=/tmp/kg_mirror_progress.log
: > "$LOG"
i=0
while IFS=$'\t' read -r col gb docs; do
  i=$((i+1))
  FREE=$(df -BG --output=avail /mnt/igvf-data | tail -1 | tr -dc '0-9')
  if [[ "$FREE" -lt "$FLOOR_GB" ]]; then
    echo "STOPPING: only ${FREE} GB free, floor is ${FLOOR_GB} GB. "\
"$((PLANNED-i+1)) collection(s) not mirrored." | tee -a "$LOG"
    break
  fi
  printf '[%d/%d] %-42s %8.3f GB  %s docs  (%s GB free)\n' \
    "$i" "$PLANNED" "$col" "$gb" "$docs" "$FREE" | tee -a "$LOG"
  # kg-mirror keeps per-collection state on disk, so a re-run resumes
  # rather than restarting -- which is what makes a multi-hour job over a
  # flaky link survivable.
  if ! "${DEX[@]}" igvfagent kg-mirror pull --collection "$col" >>"$LOG" 2>&1; then
    echo "  WARNING: $col failed; continuing with the rest" | tee -a "$LOG"
  fi
done < /tmp/kg_plan.txt

step "3. Verify"
"${DEX[@]}" igvfagent kg-mirror verify 2>&1 | tail -20 | tee -a "$LOG"
df -h /mnt/igvf-data | awk 'NR>1{printf "  volume: %s used of %s (%s)\n", $3,$2,$5}'
echo
echo "Progress log: $LOG"
