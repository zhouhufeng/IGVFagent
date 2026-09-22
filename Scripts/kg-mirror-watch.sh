#!/usr/bin/env bash
# Hourly watch for the LOCAL KG mirror (Scripts/kg-mirror-watch.sh).
#
# The mirror is a ~117-day run, so the failure that matters is not a crash --
# it is a silent stall: the process alive but no longer writing.
#
# Liveness is read from the mtime of the mirror's own log rather than from
# the size of the warehouse. The log is a single file, so the check stays
# O(1) as the warehouse grows toward ~1.9 TB, and it keeps ticking during a
# giant collection whose per-batch commits are minutes apart. Warehouse
# bytes are reported for progress, not used for the stall decision.
set -uo pipefail

ROOT=/media/hzhou/HSA/research/projects/IGVFagent
WAREHOUSE=$ROOT/Data/Warehouse/KG
LOG=$ROOT/Docs/Logs/watch-local.log
STALL_HOURS=6
DISK_PCT=90

mkdir -p "$(dirname "$LOG")"
TS=$(date '+%Y-%m-%d %H:%M:%S')
NOW=$(date +%s)

MLOG=$(ls -t "$ROOT"/Docs/Logs/kg-mirror-local-*.log 2>/dev/null | head -1)
if [ -n "$MLOG" ]; then
  AGE_H=$(( (NOW - $(stat -c %Y "$MLOG")) / 3600 ))
  COLL=$(grep -oE 'Pulling [a-z_]+' "$MLOG" 2>/dev/null | tail -1 | cut -d' ' -f2)
  DONE=$(grep -c 'Pulling ' "$MLOG" 2>/dev/null)
else
  AGE_H=0; COLL=none; DONE=0
fi

SIZE=$(du -sh "$WAREHOUSE" 2>/dev/null | cut -f1)
pgrep -f 'igvfagent kg-mirror' >/dev/null && ALIVE=1 || ALIVE=0
read -r _ _ _ _ pct _ < <(df -h "$ROOT" | tail -1); PCT=${pct%\%}

printf '%s  warehouse=%-8s coll=%-34s n=%-3s alive=%s disk=%s%% idle=%sh\n' \
  "$TS" "${SIZE:-0}" "${COLL:-?}" "$DONE" "$ALIVE" "$PCT" "$AGE_H" >> "$LOG"

alert() {
  printf '%s  ALERT: %s\n' "$TS" "$1" >> "$LOG"
  local uid; uid=$(id -u)
  DISPLAY=:0 DBUS_SESSION_BUS_ADDRESS="unix:path=/run/user/$uid/bus" \
    notify-send -u critical "IGVFagent KG mirror" "$1" 2>/dev/null
}

if [ "$ALIVE" -eq 0 ]; then
  alert "Mirror NOT running (warehouse ${SIZE:-0}, ${DONE} collections started). Resume: cd $ROOT && set -a; . ./.env; set +a; .venv/bin/igvfagent kg-mirror pull-all --include-giants"
elif [ "$AGE_H" -ge "$STALL_HOURS" ]; then
  alert "Mirror alive but its log has not been written for ${AGE_H}h (on '${COLL}') - likely stalled."
fi

[ "${PCT:-0}" -ge "$DISK_PCT" ] && alert "Local disk ${PCT}% full."
exit 0
