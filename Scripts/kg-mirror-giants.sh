#!/usr/bin/env bash
# Hand the remaining billion-row collections off from offset paging to
# planned byte-range streaming.
#
# WHY. `kg-mirror pull` pages by offset, so cost per batch grows with depth:
# genomic_elements_genes decayed from 3,400 rows/s to 530 over 147 M rows.
# The collections still to mirror hold 11.1 B rows between them, which at
# that rate is not a wait, it is a no. `pull-planned` streams each byte-wise
# key range with a single cursor, so throughput stays flat, and the ranges
# are independent enough to run several at once.
#
# It waits for the in-flight offset pull to finish its current collection
# first. That collection was 79% done when this was written, and its
# offset-written shards do not carry over to the planned layout, so
# interrupting it would throw that away.
#
# Parallelism is deliberately modest. The Catalog is a shared consortium
# service and the hosted VM is mirroring from it too; JOBS cursors is
# already several times the load one offset pull was making.
set -uo pipefail

ROOT=/media/hzhou/HSA/research/projects/IGVFagent
cd "$ROOT" || exit 1
set -a; . ./.env; set +a

WAIT_FOR="${WAIT_FOR:-genomic_elements_genes}"
JOBS="${JOBS:-4}"
TARGET_ROWS="${TARGET_ROWS:-25000000}"
POLL="${POLL:-300}"
STATE="$ROOT/Data/Warehouse/KG/_state"
PIDFILE="$ROOT/Docs/Logs/.kg-mirror.pid"
BIN="$ROOT/.venv/bin/igvfagent"
LOG="$ROOT/Docs/Logs/kg-giants-$(date +%Y%m%d_%H%M%S).log"

# Smallest first, so the cheap ones bank progress before variants_variants.
COLLECTIONS=(variants_proteins coding_variants variants_coding_variants
             coding_variants_phenotypes variants variants_variants)

mkdir -p "$(dirname "$LOG")"
log() { printf '%s  %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG"; }

# One driver at a time. Two would plan the same collection twice and race on
# the same range indexes, and there are two ways to start one: this script
# directly, or kg-mirror-run.sh once the offset phase is done (which is also
# what @reboot and kg-mirror-watch.sh call). The lock is held for the life of
# the process via fd 9 and released when it exits, however it exits.
exec 9>"$ROOT/Docs/Logs/.kg-giants.lock"
if ! flock -n 9; then
  echo "another giants driver holds the lock; nothing to do"
  exit 0
fi

done_p() {  # done_p <collection> -> 0 if its state says done
  python3 - "$1" <<'PY'
import json, sys, pathlib
p = pathlib.Path("Data/Warehouse/KG/_state") / f"{sys.argv[1]}.json"
try:
    sys.exit(0 if json.loads(p.read_text()).get("status") == "done" else 1)
except Exception:
    sys.exit(1)
PY
}

n_ranges() {  # n_ranges <collection> -> range count of its plan
  python3 - "$1" <<'PY'
import json, sys, pathlib
p = pathlib.Path("Data/Warehouse/KG/_state") / f"{sys.argv[1]}__ranges.json"
try:
    print(len(json.loads(p.read_text())["ranges"]))
except Exception:
    print(0)
PY
}

# planned_done_p <collection> -> 0 if every range in the plan is done.
#
# The collection-level state file is NOT the answer here: pull-planned writes
# only per-range state (<collection>__<tag>.json), and nothing in the planned
# path ever stamps <collection>.json as done -- that is the offset path's
# bookkeeping. Asking done_p after a planned sweep therefore reports failure
# even when every range succeeded. Completeness is the conjunction over the
# plan's own ranges, and a range whose streamed count disagreed with the plan
# is deliberately not done (status=count_mismatch), so it is caught here too.
planned_done_p() {
  python3 - "$1" <<'PY'
import json, sys, pathlib
st = pathlib.Path("Data/Warehouse/KG/_state")
name = sys.argv[1]
plan = st / f"{name}__ranges.json"
try:
    ranges = json.loads(plan.read_text())["ranges"]
except Exception:
    sys.exit(1)
if not ranges:
    sys.exit(1)
for r in ranges:
    f = st / f"{name}__{r['tag']}.json"
    try:
        if json.loads(f.read_text()).get("status") != "done":
            sys.exit(1)
    except Exception:
        sys.exit(1)
sys.exit(0)
PY
}

# n_incomplete <collection> -> how many ranges are not done (for the log)
n_incomplete() {
  python3 - "$1" <<'PY'
import json, sys, pathlib
st = pathlib.Path("Data/Warehouse/KG/_state")
name = sys.argv[1]
try:
    ranges = json.loads((st / f"{name}__ranges.json").read_text())["ranges"]
except Exception:
    print(-1); sys.exit(0)
bad = 0
for r in ranges:
    try:
        if json.loads((st / f"{name}__{r['tag']}.json").read_text()).get("status") != "done":
            bad += 1
    except Exception:
        bad += 1
print(bad)
PY
}

log "giants driver starting (JOBS=$JOBS TARGET_ROWS=$TARGET_ROWS); waiting on $WAIT_FOR"

# ---- 1. wait for the offset pull to finish the collection it is inside ----
# The ONLY exit from this wait is $WAIT_FOR reaching done. Taking over
# because the offset pull looks dead would be a trap: kg-mirror-watch.sh
# restarts it within the hour and writes a new pid, so a crash here would
# have this driver kill the *replacement* and walk away from a collection
# sitting at 79%. A dead offset pull is the watch's job to fix; this just
# says so in the log and keeps waiting.
absent=0
while ! done_p "$WAIT_FOR"; do
  oldpid=$(cat "$PIDFILE" 2>/dev/null || echo)
  if [ -z "$oldpid" ] || ! kill -0 "$oldpid" 2>/dev/null; then
    absent=$((absent + 1))
    [ $((absent % 12)) -eq 1 ] && log "offset pull not running and $WAIT_FOR unfinished; waiting for kg-mirror-watch.sh to restart it"
  else
    absent=0
  fi
  sleep "$POLL"
done
log "$WAIT_FOR is done -- taking over"

# ---- 2. stop the offset pull, so it cannot start a giant the slow way ----
oldpid=$(cat "$PIDFILE" 2>/dev/null || echo)
if [ -n "$oldpid" ] && kill -0 "$oldpid" 2>/dev/null; then
  log "stopping offset pull-all (pid $oldpid)"
  kill "$oldpid" 2>/dev/null
  for _ in $(seq 1 30); do kill -0 "$oldpid" 2>/dev/null || break; sleep 2; done
  kill -0 "$oldpid" 2>/dev/null && { log "pid $oldpid ignored SIGTERM; SIGKILL"; kill -9 "$oldpid" 2>/dev/null; }
fi
# Claim the pidfile: kg-mirror-watch.sh restarts whatever it does not find
# alive, and what it would restart is the offset pull we just stopped.
echo $$ > "$PIDFILE"
log "claimed $PIDFILE as pid $$"

# ---- 3. plan + pull each remaining collection ----
failed=0
for coll in "${COLLECTIONS[@]}"; do
  if done_p "$coll" || planned_done_p "$coll"; then
    log "$coll: already complete, skipping"; continue
  fi

  log "$coll: planning ranges (target $TARGET_ROWS rows)"
  if ! "$BIN" kg-mirror plan-ranges --collection "$coll" --target-rows "$TARGET_ROWS" >>"$LOG" 2>&1; then
    log "$coll: PLANNING FAILED -- skipping this collection"
    failed=$((failed + 1)); continue
  fi
  n=$(n_ranges "$coll")
  if [ "${n:-0}" -eq 0 ]; then
    log "$coll: plan has no ranges -- skipping"
    failed=$((failed + 1)); continue
  fi
  log "$coll: $n ranges; pulling $JOBS at a time"

  # Already-finished ranges return immediately, so a blind sweep is also the
  # retry: re-running costs one state-file read per done range.
  for i in $(seq 0 $((n - 1))); do
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n 2>/dev/null || sleep 5; done
    ( "$BIN" kg-mirror pull-planned --collection "$coll" --index "$i" >>"$LOG" 2>&1 \
        || log "$coll[index $i]: range FAILED (exit $?)" ) &
  done
  wait
  if planned_done_p "$coll"; then
    log "$coll: COLLECTION DONE ($n ranges)"
  else
    log "$coll: sweep finished with $(n_incomplete "$coll") of $n ranges incomplete -- rerun to retry them"
    failed=$((failed + 1))
  fi
done

log "giants driver finished; $failed collection(s) left incomplete"
exit $([ "$failed" -eq 0 ] && echo 0 || echo 1)
