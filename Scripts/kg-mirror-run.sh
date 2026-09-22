#!/usr/bin/env bash
# Launch the local KG mirror detached, and record its PID.
#
# No --include-giants. That flag does not add `variants`/`variants_variants`
# -- pull-all already mirrors those -- it only re-adds the two collections
# the inventory marks skip_default, genes_coding_variants_scores{,_grp}.
# Those hold ~68k documents in 36-50 GB, so roughly half a megabyte each,
# and a single batch of them OOM-killed this box on 2026-09-22. They are
# skipped by design; mirror them alone, with a small --batch-size, if they
# are ever actually needed.
#
# The PID goes to a file because `pgrep -f` matched the checking shell's own
# command line and reported a dead mirror as alive.
set -uo pipefail

ROOT=/media/hzhou/HSA/research/projects/IGVFagent
PIDFILE=$ROOT/Docs/Logs/.kg-mirror.pid
cd "$ROOT" || exit 1

if [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE" 2>/dev/null)" 2>/dev/null; then
  echo "already running (pid $(cat "$PIDFILE"))"; exit 0
fi

mkdir -p Docs/Logs
LOG=$ROOT/Docs/Logs/kg-mirror-local-$(date +%Y%m%d_%H%M%S).log

setsid nohup bash -c '
  cd '"$ROOT"' || exit 1
  set -a; . ./.env; set +a
  exec .venv/bin/igvfagent kg-mirror pull-all
' > "$LOG" 2>&1 < /dev/null &

echo $! > "$PIDFILE"
echo "launched pid $(cat "$PIDFILE")  log=$LOG"
