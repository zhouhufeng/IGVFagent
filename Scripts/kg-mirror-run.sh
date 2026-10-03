#!/usr/bin/env bash
# Launch the local KG mirror detached, and record its PID.
#
# Two phases, and this picks the right one. While any ordinary collection is
# outstanding, `pull-all` mirrors them offset-paged, which is fine at their
# size. Once those are done, only the billion-row collections remain, and
# offset paging cannot finish them (cost per batch grows with depth:
# genomic_elements_genes decayed 3,400 -> 530 rows/s over 147 M rows, and
# 11.1 B rows remain). Those go to kg-mirror-giants.sh, which plans byte-wise
# key ranges and streams each with one cursor at flat throughput.
#
# This matters for kg-mirror-watch.sh: what it restarts after a crash has to
# be the phase that is actually current, not whichever one ran first.
#
# No --include-giants on pull-all. That flag does not add
# variants/variants_variants -- pull-all already mirrors those -- it only
# re-adds the two collections the inventory marks skip_default,
# genes_coding_variants_scores{,_grp}: ~68k documents in 36-50 GB, about
# half a megabyte each, which OOM-killed this box on 2026-09-22.
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

# The offset phase is over once the last ordinary collection is done; all
# that is left then is the planned-range work.
offset_phase_done() {
  python3 - <<'PY'
import json, pathlib, sys
p = pathlib.Path("Data/Warehouse/KG/_state/genomic_elements_genes.json")
try:
    sys.exit(0 if json.loads(p.read_text()).get("status") == "done" else 1)
except Exception:
    sys.exit(1)
PY
}

if offset_phase_done; then
  LOG=$ROOT/Docs/Logs/kg-mirror-local-$(date +%Y%m%d_%H%M%S).log
  setsid nohup "$ROOT/Scripts/kg-mirror-giants.sh" > "$LOG" 2>&1 < /dev/null &
  echo $! > "$PIDFILE"
  echo "launched giants driver pid $(cat "$PIDFILE")  log=$LOG"
else
  LOG=$ROOT/Docs/Logs/kg-mirror-local-$(date +%Y%m%d_%H%M%S).log
  setsid nohup bash -c '
    cd '"$ROOT"' || exit 1
    set -a; . ./.env; set +a
    exec .venv/bin/igvfagent kg-mirror pull-all
  ' > "$LOG" 2>&1 < /dev/null &
  echo $! > "$PIDFILE"
  echo "launched offset pull-all pid $(cat "$PIDFILE")  log=$LOG"
fi
