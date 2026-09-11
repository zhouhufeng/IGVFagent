#!/usr/bin/env bash
# Wait for the KG mirror to finish, then build the BEAN-0.2.9 benchmark image.
#
# Installed and launched on the VM by:
#     bash Deploy/build-bean029-when-idle.sh --install
# which copies it there and starts it detached. Run directly it does the
# waiting and building itself.
#
# WHY A WATCHER RATHER THAN JUST BUILDING. The mirror is I/O-bound and has
# hours left; a docker build pulling a torch wheel and compiling a Cython
# extension alongside it slows both and has already halved the mirror's
# throughput once (1490 -> 513 rows/s) when work overlapped.
#
# IT MUST OUTLIVE THE SESSION THAT STARTS IT. setsid detaches it from the
# controlling terminal, so an ssh disconnect does not kill it -- that is
# exactly how the mirror itself was lost earlier: it was started from a shell
# whose environment held the only copy of the credentials, and killing the
# shell killed the job.
set -uo pipefail

LOG=/tmp/bean029_build.log
STATE=/tmp/bean029_build.state
POLL="${POLL:-300}"                 # seconds between checks
MAX_WAIT="${MAX_WAIT:-172800}"      # give up after 48 h rather than spin forever
REPO=/srv/igvfagent
IMAGE=igvfagent-bean029

say() { printf '%s  %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

if [[ "${1:-}" == "--install" ]]; then
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1
  # shellcheck source=Deploy/vm.sh
  source Deploy/vm.sh
  KEY="$(vm_key)" || { echo "no VM key — bash Deploy/vm.sh doctor" >&2; exit 1; }
  SSH=(ssh -i "$KEY" -o IdentitiesOnly=yes -o StrictHostKeyChecking=accept-new "$(vm_host)")
  echo "== copying the watcher and the Dockerfile to the VM"
  base64 < Deploy/build-bean029-when-idle.sh | tr -d '\n' \
    | "${SSH[@]}" "base64 -d > /tmp/build-bean029-when-idle.sh"
  base64 < Deploy/Dockerfile.bean029 | tr -d '\n' \
    | "${SSH[@]}" "base64 -d > $REPO/Deploy/Dockerfile.bean029"
  echo "== launching it detached (setsid: survives this ssh session ending)"
  # No pkill here. The remote shell's OWN command line contains this
  # script's name, so `pkill -f build-bean029-when-idle` matches the shell
  # that is about to launch the watcher and kills it before it starts --
  # which is exactly what happened the first time. Stale watchers are
  # stopped by PID instead, resolved before the new one is launched.
  "${SSH[@]}" "for pid in \$(pgrep -f 'bash /tmp/build-bean029-when-idle.sh'); do \
       kill \$pid 2>/dev/null; done; \
     setsid nohup bash /tmp/build-bean029-when-idle.sh > /dev/null 2>&1 < /dev/null & \
     sleep 3; \
     pgrep -f 'bash /tmp/build-bean029-when-idle.sh' >/dev/null \
       && echo '  watcher running' || echo '  WATCHER DID NOT START'"
  echo "== check on it later with:"
  echo "     bash Deploy/vm.sh ssh 'cat $STATE; tail -20 $LOG'"
  exit 0
fi

echo "waiting" > "$STATE"
say "watcher started; polling every ${POLL}s for the KG mirror to finish"
waited=0
while pgrep -f "[m]irror-kg.sh" >/dev/null; do
  if (( waited >= MAX_WAIT )); then
    say "GIVING UP: mirror still running after ${MAX_WAIT}s. Not building."
    echo "gave_up_waiting" > "$STATE"
    exit 1
  fi
  sleep "$POLL"
  waited=$(( waited + POLL ))
  # A heartbeat every half hour, so a stuck watcher is distinguishable from a
  # patient one when someone reads the log days later.
  if (( waited % 1800 == 0 )); then
    say "still waiting (${waited}s); mirror at: $(tail -1 /tmp/kg_mirror_progress.log 2>/dev/null | cut -c1-80)"
  fi
done

say "mirror finished after ${waited}s of waiting; starting the image build"
echo "building" > "$STATE"
cd "$REPO" || { say "FAILED: no $REPO"; echo "failed_no_repo" > "$STATE"; exit 1; }

if ! docker build -f Deploy/Dockerfile.bean029 -t "$IMAGE" . >> "$LOG" 2>&1; then
  say "BUILD FAILED — see $LOG"
  echo "build_failed" > "$STATE"
  exit 1
fi
say "build ok: $IMAGE"

# Prove the stack is what the paper ran, not merely that the image exists.
# The Dockerfile asserts this too, but re-checking here means the state file
# is trustworthy on its own.
if docker run --rm --entrypoint python "$IMAGE" -c "
import torch, pyro, numpy, bean
assert torch.__version__.startswith('1.'), torch.__version__
assert numpy.__version__.startswith('1.'), numpy.__version__
print('torch', torch.__version__, 'pyro', pyro.__version__, 'numpy', numpy.__version__)
" >> "$LOG" 2>&1; then
  say "stack verified: torch 1.x + numpy 1.x + bean 0.2.9"
else
  say "IMAGE BUILT BUT STACK CHECK FAILED — see $LOG"
  echo "stack_check_failed" > "$STATE"
  exit 1
fi

# The run that motivated all of this: MixtureNormal, which NaNs on torch 2
# under both bean 1.2.2 and bean 0.2.9. If it fits here, torch was the cause.
say "running BEAN-Reporter (MixtureNormal) on the paper's LDLvar deposit"
echo "running_reporter" > "$STATE"
OUT=/mnt/igvf-data/Docs/BEANbenchmark/v029_torch1_reporter
rm -rf "$OUT"
if docker run --rm -v /mnt/igvf-data:/workspace "$IMAGE" \
      variant /workspace/Data/Benchmarks/BEANpaper/bean_count_LDLvar_annotated.h5ad \
      -o "${OUT/\/mnt\/igvf-data//workspace}" --ignore-bcmatch >> "$LOG" 2>&1; then
  say "MixtureNormal FIT SUCCEEDED under torch 1.x — the NaN was torch 2"
  echo "reporter_ok" > "$STATE"
else
  # A failure here is a real result too: it would mean torch was not the
  # cause, and the next suspect is the data itself.
  say "MixtureNormal still failed under torch 1.x — torch was NOT the cause"
  echo "reporter_failed" > "$STATE"
fi
say "done. results under $OUT"
