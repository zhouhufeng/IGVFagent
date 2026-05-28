#!/usr/bin/env bash
# Interactive terminal-side dialog with IGVFagent.
#
# Usage:
#   bash Scripts/igvfagent_repl.sh                         # zsh/bash, OpenAI gpt-5 default
#   IGVF_LLM_MODEL=gpt-4o-mini bash Scripts/igvfagent_repl.sh   # faster + cheaper
#   IGVF_LLM_BACKEND=ollama IGVF_LLM_MODEL=qwen3.6:35b-a3b-coding-bf16 \
#       bash Scripts/igvfagent_repl.sh                     # local, no network
#
# Commands inside the REPL:
#   :q  /  :quit  /  :exit       leave the REPL
#   :model <name>                switch model on the fly
#   :iter <n>                    change max-iterations cap (default 12)
#   :quiet                       toggle per-step trace printing
#   :history                     list this session's transcripts
#   :help                        print this menu

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT"

# Source .env files (repo, home, or both) so the API key reaches the agent.
[ -f "$ROOT/.env" ] && { set -a; source "$ROOT/.env"; set +a; }
[ -f "$HOME/.env" ] && { set -a; source "$HOME/.env"; set +a; }

# Find a working `igvfagent` executable. Tries, in order:
#   1. $ROOT/.venv/bin/igvfagent          (script's own repo .venv)
#   2. .venv/bin/igvfagent in any worktree under $ROOT/.claude/worktrees/
#   3. system `igvfagent` on PATH
#   4. .venv/bin/python -m igvfagent      (as a last resort)
find_igvfagent() {
    if [ -x "$ROOT/.venv/bin/igvfagent" ]; then
        echo "$ROOT/.venv/bin/igvfagent"; return 0
    fi
    if [ -d "$ROOT/.claude/worktrees" ]; then
        for d in "$ROOT/.claude/worktrees"/*/; do
            if [ -x "$d.venv/bin/igvfagent" ]; then
                echo "$d.venv/bin/igvfagent"; return 0
            fi
        done
    fi
    if command -v igvfagent >/dev/null 2>&1; then
        echo "$(command -v igvfagent)"; return 0
    fi
    if [ -x "$ROOT/.venv/bin/python" ]; then
        echo "$ROOT/.venv/bin/python -m igvfagent"; return 0
    fi
    return 1
}

IGVFAGENT_BIN="$(find_igvfagent)" || {
    printf 'error: cannot find an `igvfagent` executable.\n' >&2
    printf '  searched:\n' >&2
    printf '    %s/.venv/bin/igvfagent\n' "$ROOT" >&2
    printf '    %s/.claude/worktrees/*/.venv/bin/igvfagent\n' "$ROOT" >&2
    printf '    $(command -v igvfagent)\n' >&2
    printf '\nTo fix:\n' >&2
    printf '  cd %s/.claude/worktrees/festive-volhard-60dea7/\n' "$ROOT" >&2
    printf '  bash Scripts/igvfagent_repl.sh\n' >&2
    printf '\n(or `pip install -e .` into a .venv at the repo root)\n' >&2
    exit 2
}

BACKEND="${IGVF_LLM_BACKEND:-openai}"
MODEL="${IGVF_LLM_MODEL:-gpt-5}"
MAX_ITER="${IGVF_AGENT_MAX_ITER:-12}"
MAX_TOK="${IGVF_AGENT_MAX_TOK:-4096}"
TEMP="${IGVF_AGENT_TEMP:-0.0}"
QUIET_FLAG=""

printf '\n┌──────────────────────────────────────────────────────────────┐\n'
printf '│ IGVFagent terminal dialog                                    │\n'
printf '│ backend: %-10s · model: %-30s │\n' "$BACKEND" "$MODEL"
printf '│ max_iter: %-3d · max_tok: %-5d · temp: %-3s · :help for menu │\n' "$MAX_ITER" "$MAX_TOK" "$TEMP"
printf '└──────────────────────────────────────────────────────────────┘\n'
printf '  exec: %s\n\n' "$IGVFAGENT_BIN"

while true; do
    # Read a multi-line-friendly prompt (single line, ends on Enter)
    printf '\033[1;36m›\033[0m '
    if ! IFS= read -r line; then
        printf '\n(EOF — bye)\n'
        break
    fi
    case "$line" in
        ""|:|':help')
            cat <<EOF
  :q | :quit | :exit       leave the REPL
  :model <name>            switch LLM (e.g. :model gpt-4o-mini)
  :iter <n>                set max-iterations (e.g. :iter 25)
  :quiet                   toggle per-step trace
  :history                 list this session's transcripts
EOF
            continue
            ;;
        ':q'|':quit'|':exit')
            printf 'bye.\n'
            break
            ;;
        ':model '*)
            MODEL="${line#:model }"
            printf '  ✓ model -> %s\n' "$MODEL"
            continue
            ;;
        ':iter '*)
            MAX_ITER="${line#:iter }"
            printf '  ✓ max_iter -> %s\n' "$MAX_ITER"
            continue
            ;;
        ':quiet')
            if [ -z "$QUIET_FLAG" ]; then
                QUIET_FLAG="--quiet"
                printf '  ✓ quiet mode ON\n'
            else
                QUIET_FLAG=""
                printf '  ✓ quiet mode OFF\n'
            fi
            continue
            ;;
        ':history')
            ls -t Docs/Agent/ 2>/dev/null | head -10 | sed 's/^/  /'
            continue
            ;;
    esac
    # Fire the agent — every turn is an independent run; if you want
    # cross-turn memory, use the Streamlit UI which keeps a session.
    $IGVFAGENT_BIN ask \
        --backend "$BACKEND" \
        --model "$MODEL" \
        --max-iterations "$MAX_ITER" \
        --max-tokens "$MAX_TOK" \
        --temperature "$TEMP" \
        $QUIET_FLAG \
        "$line"
    echo
done
