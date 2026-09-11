#!/usr/bin/env bash
# Copy a module into the running container, and restart it if -- and only if
# -- the change cannot take effect without one.
#
#   bash Deploy/ship.sh Scripts/kg_visualizer.py
#   bash Deploy/ship.sh Scripts/base_editing_screen.py Scripts/_stats.py
#   bash Deploy/ship.sh --check Scripts/kg_visualizer.py   # say, change nothing
#
# WHY THIS EXISTS. The README used to state the rule as a list of three
# exceptions -- "NOT _agent.py, _llm.py or _tools.py" -- and that list is
# wrong the moment a fourth module is imported by the app. It cost a live
# fix: kg_visualizer.py was copied in, announced as live because it was not
# on the list, and was not live at all. The Streamlit process had imported it
# at start-up and cached it in sys.modules; only `streamlit_app.py` itself is
# re-read per rerun.
#
# So the question is not "is this module on a list" but "does the long-lived
# process import it", which is computable. This script computes it from the
# app's own import graph rather than from memory.
set -euo pipefail

cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." || exit 1
# shellcheck source=Deploy/vm.sh
source Deploy/vm.sh

CHECK_ONLY=0
if [[ "${1:-}" == "--check" ]]; then CHECK_ONLY=1; shift; fi
[[ $# -gt 0 ]] || { echo "usage: bash Deploy/ship.sh [--check] Scripts/<module>.py ..." >&2; exit 2; }

DEST=/opt/venv/lib/python3.11/site-packages/igvfagent
CONTAINER="${CONTAINER:-igvfagent-app}"

# Which modules does the Streamlit app pull into its own process? Computed by
# walking the import graph from streamlit_app.py through the local package,
# so a module that becomes imported later is classified correctly without
# anyone remembering to update a list.
IN_PROCESS=$(python3 - <<'PY'
import ast, pathlib, sys
root = pathlib.Path("Scripts")
seen, queue = set(), ["streamlit_app"]
local = {p.stem for p in root.glob("*.py")}
while queue:
    name = queue.pop()
    if name in seen or name not in local:
        continue
    seen.add(name)
    try:
        tree = ast.parse((root / f"{name}.py").read_text())
    except (OSError, SyntaxError):
        continue
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                queue.append(a.name.split(".")[-1])
        elif isinstance(node, ast.ImportFrom):
            # `from igvfagent import kg_visualizer` and `import kg_visualizer`
            # both matter; so does `from . import x`.
            if node.module in (None, "igvfagent") or (
                    node.module or "").startswith("igvfagent"):
                for a in node.names:
                    queue.append(a.name.split(".")[-1])
            base = (node.module or "").split(".")[-1]
            if base in local:
                queue.append(base)
print(" ".join(sorted(seen)))
PY
)

needs_restart=0
declare -a COPIED=()
for path in "$@"; do
  [[ -f "$path" ]] || { echo "no such file: $path" >&2; exit 2; }
  stem=$(basename "$path" .py)
  if [[ "$stem" == "streamlit_app" ]]; then
    # The one genuine exception: Streamlit re-executes its MAIN SCRIPT on
    # every rerun. It is the root of the import graph, so the walk above
    # necessarily contains it -- but it is the only member that is re-read.
    verdict="main script — re-read on the next interaction"
    printf '  %-34s %s\n' "$stem" "$verdict"
    COPIED+=("$path")
    continue
  fi
  if grep -qw -- "$stem" <<<"$IN_PROCESS"; then
    verdict="IN-PROCESS — needs a restart to take effect"
    needs_restart=1
  else
    verdict="fresh subprocess per call — live on copy"
  fi
  printf '  %-34s %s\n' "$stem" "$verdict"
  COPIED+=("$path")
done

if [[ "$CHECK_ONLY" == 1 ]]; then
  echo
  [[ "$needs_restart" == 1 ]] \
    && echo "--check: at least one module needs a restart. Nothing copied." \
    || echo "--check: all of these go live on copy. Nothing copied."
  exit 0
fi

KEY_FILE="$(vm_key)" || { echo "no VM key — bash Deploy/vm.sh doctor" >&2; exit 1; }
SSH=(ssh -i "$KEY_FILE" -o IdentitiesOnly=yes -o StrictHostKeyChecking=accept-new "$(vm_host)")

for path in "${COPIED[@]}"; do
  base=$(basename "$path")
  base64 < "$path" | tr -d '\n' | "${SSH[@]}" "base64 -d > /tmp/$base"
  "${SSH[@]}" "docker cp /tmp/$base $CONTAINER:$DEST/$base >/dev/null && rm -f /tmp/$base"
done
"${SSH[@]}" "docker exec $CONTAINER sh -c 'rm -rf $DEST/__pycache__'"
echo "  copied ${#COPIED[@]} module(s)"

if [[ "$needs_restart" == 1 ]]; then
  echo "  restarting (one or more modules are imported in-process)…"
  # The KG mirror runs as a docker exec and dies with the container. It
  # resumes from per-collection state, but it has to be relaunched, and
  # silently losing it is how a multi-hour job goes missing.
  MIRROR_WAS_RUNNING=$("${SSH[@]}" 'pgrep -f "[m]irror-kg.sh" >/dev/null && echo yes || echo no')
  [[ "$MIRROR_WAS_RUNNING" == yes ]] && "${SSH[@]}" 'pkill -f "[m]irror-kg.sh" || true'
  "${SSH[@]}" "docker restart $CONTAINER >/dev/null"
  "${SSH[@]}" "sleep 15; docker exec $CONTAINER python3 -c \"
import urllib.request
print('  health:', urllib.request.urlopen('http://127.0.0.1:8501/_stcore/health', timeout=20).status)\""
  if [[ "$MIRROR_WAS_RUNNING" == yes ]]; then
    "${SSH[@]}" 'cd /srv/igvfagent && nohup bash Deploy/mirror-kg.sh > /tmp/kg_mirror_run.log 2>&1 & sleep 6; pgrep -f "[m]irror-kg.sh" >/dev/null && echo "  KG mirror relaunched" || echo "  WARNING: mirror did NOT relaunch"'
  fi
else
  echo "  no restart needed — these take effect on the next call"
fi

echo
echo "Reminder: a hot copy lives in the container, not the image. Follow it"
echo "with Deploy/redeploy.sh at the next window or a rebuild reverts it."
