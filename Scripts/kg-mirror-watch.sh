#!/usr/bin/env bash
# Hourly watch for the LOCAL KG mirror (Scripts/kg-mirror-watch.sh).
#
# Liveness comes from the PID file written by kg-mirror-run.sh, not from
# `pgrep -f`: the pattern matched the checking shell's own command line and
# cheerfully reported a mirror that the OOM killer had taken hours earlier.
#
# Progress comes from the mtime of the mirror's own log, so the check stays
# O(1) as the warehouse grows toward ~1.8 TB and keeps ticking through a
# collection whose batches commit minutes apart. Warehouse bytes are
# reported, never used to decide a stall.
#
# A dead mirror is restarted rather than merely reported: the run is weeks
# long, pull-all resumes from per-collection state, and the common cause is
# an OOM that a retry survives. RESTART_CAP bounds a crash loop -- past it,
# the watch stops restarting and just complains.
set -uo pipefail

ROOT=/media/hzhou/HSA/research/projects/IGVFagent
PIDFILE=$ROOT/Docs/Logs/.kg-mirror.pid
LOG=$ROOT/Docs/Logs/watch-local.log
RESTARTS=$ROOT/Docs/Logs/.kg-mirror.restarts
STALL_HOURS=6
DISK_PCT=90
RESTART_CAP=3          # per calendar day

mkdir -p "$(dirname "$LOG")"
TS=$(date '+%Y-%m-%d %H:%M:%S'); NOW=$(date +%s); TODAY=$(date +%F)

PID=$(cat "$PIDFILE" 2>/dev/null || echo "")
if [ -n "$PID" ] && kill -0 "$PID" 2>/dev/null; then ALIVE=1; else ALIVE=0; fi

MLOG=$(ls -t "$ROOT"/Docs/Logs/kg-mirror-local-*.log 2>/dev/null | head -1)
if [ -n "$MLOG" ]; then
  AGE_H=$(( (NOW - $(stat -c %Y "$MLOG")) / 3600 ))
  COLL=$(grep -oE 'Pulling [a-z_]+' "$MLOG" 2>/dev/null | tail -1 | cut -d' ' -f2)
else
  AGE_H=0; COLL=none
fi

# A quiet log is not a stalled mirror. In the planned-range phase a cursor
# logs once per finished range -- 25 M rows, often hours -- so log mtime
# alone would cry stall through perfectly healthy work. Shards, though, land
# continuously. `-quit` stops at the first hit, so this costs a directory
# walk and not a scan of the whole warehouse.
FRESH_SHARD=$(find "$ROOT/Data/Warehouse/KG" -name '*.parquet' \
                -newermt "-${STALL_HOURS} hours" -print -quit 2>/dev/null)
[ -n "$FRESH_SHARD" ] && AGE_H=0
SIZE=$(du -sh "$ROOT/Data/Warehouse/KG" 2>/dev/null | cut -f1)
read -r _ _ _ _ pct _ < <(df -h "$ROOT" | tail -1); PCT=${pct%\%}

printf '%s  warehouse=%-8s coll=%-34s alive=%s disk=%s%% idle=%sh\n' \
  "$TS" "${SIZE:-0}" "${COLL:-?}" "$ALIVE" "$PCT" "$AGE_H" >> "$LOG"

alert() {
  printf '%s  ALERT: %s\n' "$TS" "$1" >> "$LOG"
  local uid; uid=$(id -u)
  DISPLAY=:0 DBUS_SESSION_BUS_ADDRESS="unix:path=/run/user/$uid/bus" \
    notify-send -u critical "IGVFagent KG mirror" "$1" 2>/dev/null
}

if [ "$ALIVE" -eq 0 ]; then
  read -r RDAY RCOUNT < "$RESTARTS" 2>/dev/null || { RDAY=$TODAY; RCOUNT=0; }
  [ "${RDAY:-}" = "$TODAY" ] || { RDAY=$TODAY; RCOUNT=0; }
  if [ "${RCOUNT:-0}" -lt "$RESTART_CAP" ]; then
    echo "$TODAY $((RCOUNT + 1))" > "$RESTARTS"
    OUT=$("$ROOT/Scripts/kg-mirror-run.sh" 2>&1 | tail -1)
    alert "Mirror was dead on '${COLL}' - restarted ($((RCOUNT + 1))/$RESTART_CAP today). $OUT"
  else
    alert "Mirror dead on '${COLL}' and already restarted $RCOUNT times today - NOT restarting. Check $MLOG and dmesg for OOM."
  fi
elif [ "$AGE_H" -ge "$STALL_HOURS" ]; then
  alert "Mirror alive (pid $PID) but its log is ${AGE_H}h cold on '${COLL}' - likely stalled."
fi

[ "${PCT:-0}" -ge "$DISK_PCT" ] && alert "Local disk ${PCT}% full."

# Right phase? Once the last ordinary collection is done, only billion-row
# collections remain and offset paging cannot finish them. If pull-all is
# still the live process at that point, the handoff to kg-mirror-giants.sh
# did not happen (its waiter died, say) and the mirror is grinding a giant
# the one way that never ends. That is a silent failure: alive=1, shards
# landing, nothing else complains. Alert rather than act -- restarting into
# the other phase means killing a running pull, which is the user's call.
if [ "$ALIVE" -eq 1 ] && python3 - <<'PY'
import json, pathlib, sys
p = pathlib.Path("/media/hzhou/HSA/research/projects/IGVFagent/Data/Warehouse/KG/_state/genomic_elements_genes.json")
try:
    sys.exit(0 if json.loads(p.read_text()).get("status") == "done" else 1)
except Exception:
    sys.exit(1)
PY
then
  if tr '\0' ' ' < "/proc/$PID/cmdline" 2>/dev/null | grep -q 'pull-all'; then
    alert "Offset phase is done but pid $PID is still pull-all -- the planned-range handoff did not happen. Start it: Scripts/kg-mirror-giants.sh"
  fi
fi
exit 0
