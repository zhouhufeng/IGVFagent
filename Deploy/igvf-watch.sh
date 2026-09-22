#!/usr/bin/env bash
# Daily health watch for igvfagent-prod.
#   - warns once per crossing when /mnt/igvf-data passes THRESHOLD%
#   - warns if the KG mirror stopped before it finished
# The log is the reliable channel; mail from a cloud VM may be filtered.
set -uo pipefail

MOUNT=/mnt/igvf-data
THRESHOLD=85
TO=hufengzhou@g.harvard.edu
LOG=$MOUNT/Docs/Logs/watch.log
STATE=/var/lib/igvf-watch.state

mkdir -p "$(dirname "$LOG")" 2>/dev/null
touch "$STATE" 2>/dev/null

read -r _ size used avail pct _ < <(df -h "$MOUNT" | tail -1)
PCT=${pct%\%}
KG=$(du -sh "$MOUNT/Data/KG" 2>/dev/null | cut -f1)
if pgrep -f 'mirror-kg.sh' >/dev/null; then MIRROR=running; else MIRROR=stopped; fi

TS=$(date '+%Y-%m-%d %H:%M:%S')
printf '%s  disk=%s/%s (%s%%)  free=%s  KG=%s  mirror=%s\n' \
  "$TS" "$used" "$size" "$PCT" "$avail" "${KG:-n/a}" "$MIRROR" >> "$LOG"

notify() {
  grep -qx "$1" "$STATE" 2>/dev/null && return 0
  echo "$1" >> "$STATE"
  printf '%s\n' "$3" | mail -s "$2" "$TO" 2>/dev/null
  printf '%s  ALERT SENT: %s\n' "$TS" "$2" >> "$LOG"
}

if [ "$PCT" -ge "$THRESHOLD" ]; then
  notify "disk$PCT" "[igvfagent-prod] disk ${PCT}% - time to extend the volume" \
"$MOUNT is ${PCT}% full (${used} of ${size}, ${avail} free).

Plan: extend the Cinder volume igvfagent-data in Horizon, then:
  sudo xfs_growfs $MOUNT

Stop the stack and detach the volume before extending."
else
  sed -i '/^disk/d' "$STATE" 2>/dev/null
fi

if [ "$MIRROR" = stopped ] && ! grep -q 'MIRROR COMPLETE' "$MOUNT"/Docs/Logs/kg-mirror-*.log 2>/dev/null; then
  notify "mirror-$(date +%F)" "[igvfagent-prod] KG mirror is not running" \
"No mirror-kg.sh process, and no completion marker in the mirror logs.
Latest log: $(ls -t "$MOUNT"/Docs/Logs/kg-mirror-*.log 2>/dev/null | head -1)

Re-run to resume - per-collection state is on disk:
  cd /srv/igvfagent && bash Deploy/mirror-kg.sh"
fi
