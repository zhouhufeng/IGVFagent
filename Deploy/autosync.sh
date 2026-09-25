#!/usr/bin/env bash
# Keep a working checkout on a Mac up to date with GitHub, in the background.
#
#   bash Deploy/autosync.sh install     # every 5 min via launchd (this checkout)
#   bash Deploy/autosync.sh uninstall
#   bash Deploy/autosync.sh status      # agent loaded? last log lines
#   bash Deploy/autosync.sh once        # one sync now (what the agent runs)
#
# Run `install` once on each machine. A sync only ever fast-forwards the
# current branch from its upstream, carrying uncommitted edits along
# (autostash). It never merges, rebases, pushes or discards anything: when
# the machine has local commits AND GitHub has new ones, or the stash would
# not re-apply cleanly, it stops and posts a notification so you resolve it
# by hand. Pushing stays a deliberate act.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LABEL="org.igvfagent.autosync.$(printf '%s' "$ROOT" | shasum | cut -c1-8)"
PLIST="$HOME/Library/LaunchAgents/$LABEL.plist"
LOG="$HOME/Library/Logs/igvfagent-autosync.log"
INTERVAL="${AUTOSYNC_INTERVAL:-300}"

notify() {
    osascript -e "display notification \"$1\" with title \"IGVFagent git sync\"" >/dev/null 2>&1 || true
}

sync_once() {
    cd "$ROOT"
    local ts; ts=$(date '+%F %T')
    # Mid-merge/rebase/cherry-pick: someone is working, leave it alone.
    local gd; gd=$(git rev-parse --git-dir)
    if [[ -e "$gd/MERGE_HEAD" || -d "$gd/rebase-merge" || -d "$gd/rebase-apply" || -e "$gd/CHERRY_PICK_HEAD" ]]; then
        echo "$ts skip: merge/rebase in progress"; return 0
    fi
    local up; up=$(git rev-parse --abbrev-ref --symbolic-full-name '@{u}' 2>/dev/null) || {
        echo "$ts skip: $(git branch --show-current) has no upstream"; return 0; }
    git fetch --quiet --prune origin || { echo "$ts fetch failed (offline?)"; return 0; }
    local behind ahead
    behind=$(git rev-list --count "HEAD..$up")
    ahead=$(git rev-list --count "$up..HEAD")
    if [[ "$behind" == 0 ]]; then
        [[ "$ahead" != 0 ]] && echo "$ts up to date; $ahead local commit(s) not pushed"
        return 0
    fi
    if [[ "$ahead" != 0 ]]; then
        echo "$ts DIVERGED: $ahead local, $behind on $up; not touching it"
        notify "Diverged from $up ($ahead local, $behind remote): pull --rebase by hand"
        return 0
    fi
    # An incoming commit that touches a file edited here would leave conflict
    # markers in it; wait for the edit to be committed (then it is DIVERGED).
    local overlap
    overlap=$(comm -12 <(git diff --name-only HEAD | sort) <(git diff --name-only HEAD "$up" | sort) | head -5)
    if [[ -n "$overlap" ]]; then
        echo "$ts WAITING: $behind commit(s) on $up change files edited here: $(echo $overlap)"
        notify "$up changes files you are editing ($(echo $overlap | cut -c1-60)): commit, then pull --rebase"
        return 0
    fi
    if git -c merge.autoStash=true merge --ff-only --quiet "$up" >/dev/null 2>&1; then
        echo "$ts fast-forwarded $behind commit(s) from $up -> $(git rev-parse --short HEAD)"
    else
        echo "$ts FAILED to fast-forward (uncommitted edits conflict?); see git status / git stash list"
        notify "Could not fast-forward from $up: check git status and git stash list"
    fi
}

case "${1:-}" in
    once) mkdir -p "$(dirname "$LOG")"; sync_once >>"$LOG" 2>&1; tail -1 "$LOG" ;;
    install)
        mkdir -p "$(dirname "$PLIST")" "$(dirname "$LOG")"
        cat >"$PLIST" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
  <key>Label</key><string>$LABEL</string>
  <key>ProgramArguments</key><array>
    <string>/bin/bash</string><string>$ROOT/Deploy/autosync.sh</string><string>once</string>
  </array>
  <key>StartInterval</key><integer>$INTERVAL</integer>
  <key>RunAtLoad</key><true/>
  <key>EnvironmentVariables</key><dict>
    <key>PATH</key><string>/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin</string>
  </dict>
  <key>StandardOutPath</key><string>/dev/null</string>
  <key>StandardErrorPath</key><string>$LOG</string>
</dict></plist>
EOF
        launchctl bootout "gui/$(id -u)/$LABEL" 2>/dev/null || true
        launchctl bootstrap "gui/$(id -u)" "$PLIST"
        echo "installed $LABEL: syncs $ROOT every ${INTERVAL}s; log $LOG" ;;
    uninstall)
        launchctl bootout "gui/$(id -u)/$LABEL" 2>/dev/null || true
        rm -f "$PLIST"; echo "removed $LABEL" ;;
    status)
        launchctl print "gui/$(id -u)/$LABEL" >/dev/null 2>&1 && echo "loaded: $LABEL" || echo "not loaded"
        [[ -f "$LOG" ]] && tail -5 "$LOG" || true ;;
    *) sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; exit 2 ;;
esac
