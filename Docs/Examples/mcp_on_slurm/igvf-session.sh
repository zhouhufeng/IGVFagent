#!/usr/bin/env bash
# Claude Code + the IGVFagent MCP server on a SLURM COMPUTE node.
#
# Run from a login node, inside tmux so the session survives a dropped
# connection:
#
#   tmux new -s igvf                          # later: tmux attach -t igvf
#   bash igvf-session.sh                      # the defaults below
#   CPUS=32 MEM=128G HOURS=48 PARTITION=shared bash igvf-session.sh
#
# An MCP server over stdio runs where its client runs, so everything the
# agent's tools do (downloads, alignment, R/Python pipelines) happens on the
# allocated node, not on the login node.
#
# Edit the four paths/defaults marked EDIT for your cluster.
set -euo pipefail

PARTITION="${PARTITION:-test}"        # EDIT: a partition you may use interactively
CPUS="${CPUS:-6}"
MEM="${MEM:-64G}"
HOURS="${HOURS:-8}"                   # must fit the partition's time limit
REPO="${IGVFAGENT_REPO:-$HOME/IGVFagent}"                       # EDIT: your checkout
ENV_BIN="${IGVFAGENT_ENV_BIN:-$HOME/envs/igvfagent/bin}"        # EDIT: env with `igvfagent`
CLAUDE="${CLAUDE_BIN:-$HOME/.local/bin/claude}"                 # EDIT: Claude Code binary

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  # Inside the allocation: environment first on PATH (the repo's .mcp.json
  # starts `igvfagent mcp serve` by name), then Claude Code.
  cd "$REPO"
  export PATH="$ENV_BIN:$PATH"
  # Many clusters set OMP_NUM_THREADS=1 for every job; let tools use the CPUs
  # this job was given.
  export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
  echo "IGVFagent session on $(hostname): ${SLURM_CPUS_PER_TASK:-?} CPUs, job $SLURM_JOB_ID"
  exec "$CLAUDE"
fi

echo "Requesting $CPUS CPUs, $MEM, ${HOURS} h on '$PARTITION' ..."
exec srun --pty -p "$PARTITION" -c "$CPUS" --mem="$MEM" \
  -t "$((HOURS / 24))-$((HOURS % 24)):00:00" -J igvf-session bash "$0"
