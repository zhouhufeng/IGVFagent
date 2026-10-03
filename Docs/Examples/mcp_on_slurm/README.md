# Running the IGVFagent MCP server on a SLURM cluster

How to use IGVFagent's tools from Claude Code on an HPC cluster, so that paper
reproductions and raw-data pipelines run on a **compute node** rather than a
login node. Tested on Harvard FASRC (Cannon) on 2026-09-25; the same pattern
works on any SLURM cluster.

**The one thing to know:** an MCP server over stdio is started by its client
and runs on the same machine. If Claude Code runs on a login node, every tool
the agent calls (downloads, alignment, R and Python pipelines) runs on the login
node too. So start Claude Code **inside a SLURM allocation**, and the MCP server
and all its tools run on that compute node.

```
login node                         compute node (SLURM allocation)
──────────                         ───────────────────────────────
tmux ─► igvf-session.sh ─► srun ─► claude ─► igvfagent mcp serve ─► tools
                                   (Claude Code)   (MCP server)      (run here)
```

## Contents

- [`igvf-session.sh`](igvf-session.sh): the launcher. It asks SLURM for a node and starts Claude Code there, with the environment set up.

## One-time setup

1. **An IGVFagent environment built on the cluster.** Don't use one copied from a laptop: a macOS `.venv` fails on Linux with "Exec format error". For example:
   ```bash
   conda create -y -p ~/envs/igvfagent python=3.11 pip
   ~/envs/igvfagent/bin/pip install -e '/path/to/IGVFagent[analysis,ui,llm,warehouse]'
   ```
2. **Claude Code installed and logged in** on the cluster, e.g. `~/.local/bin/claude`. Your home directory is shared by login and compute nodes, so logging in once is enough.
3. **The MCP server registered.** The repository's [`.mcp.json`](../../../.mcp.json) does this for everyone who opens the repo in Claude Code. It starts `igvfagent mcp serve` by name, so the environment must be first on `PATH`; the launcher arranges that.
   - Keep **one** registration. If you also ran `claude mcp add igvfagent …` yourself, Claude Code warns that the server is defined in two scopes. Remove one with `claude mcp remove igvfagent -s local`.
4. **Edit the four `EDIT` lines** in `igvf-session.sh`: the repository path, the environment's `bin/`, the Claude Code binary, and a partition you may use interactively.

## Starting a session

```bash
tmux new -s igvf                                   # survives a dropped connection
bash Docs/Examples/mcp_on_slurm/igvf-session.sh    # defaults: 6 CPUs, 64 GB, 8 h on "test"
CPUS=32 MEM=128G HOURS=48 PARTITION=shared bash Docs/Examples/mcp_on_slurm/igvf-session.sh
```

When the node is allocated, Claude Code opens in the repository. **The first time**, it asks whether to approve the `igvfagent` MCP server from `.mcp.json`; answer yes. Then:

- `/mcp` lists the server and its tools (128 by default).
- Ask as usual, e.g. *"Reproduce the key results of <paper> with the IGVFagent tools."* Tool calls appear as `mcp__igvfagent__…`, and they run on the compute node.
- **Detach:** press Ctrl-b, then d. **Return:** `tmux attach -t igvf`, on the **same login node**, because tmux lives there.

### Checking it headlessly

To confirm a node can start the server without an interactive session:

```bash
srun -p test -c 2 --mem=8G -t 00:15:00 bash -c \
  'cd /path/to/IGVFagent && export PATH=/path/to/envs/igvfagent/bin:$PATH && claude mcp list'
# igvfagent: igvfagent mcp serve - ✔ Connected
#   (or "Pending approval" until you approve it once in an interactive session)
```

## Choosing resources

| Work | Suggested request |
|---|---|
| Catalog / Portal queries, dataset explanation, small analyses | 2–4 CPUs, 16 GB, a few hours |
| Single-cell or multiome analysis from processed matrices | 8–16 CPUs, 64–128 GB |
| Paper reproduction from raw reads (FASTQ → alignment → analysis) | 16–32 CPUs, 128 GB+, long enough for the pipeline, on a partition with local scratch |

- **Match the time to the partition's limit.** FASRC's `test` partition allows 12 h, `shared` 3 days, lab partitions often more.
- **Put scratch on fast storage.** A single raw-data reproduction can write hundreds of GB: one spatial Hi-C run wrote 600 GB of uncompressed FASTQ. Keep IGVFagent's `Data/` and `Docs/` on project storage with room to spare.
- **One session is one node.** For work bigger than one node, ask the agent to submit the heavy steps as their own SLURM jobs (`sbatch`), and keep the session for planning and checking results.

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `Exec format error` when a tool runs | the environment was copied from another OS (e.g. a laptop's `.venv`). Build one on the cluster |
| `igvfagent: … ✘ Failed to connect` | `igvfagent` is not on `PATH` in that shell. Start through the launcher, or put the environment's `bin/` first on `PATH` |
| "defined in multiple scopes with different endpoints" | two registrations (project `.mcp.json` and a personal `claude mcp add`). Keep one |
| "Pending approval" | approve the project server once in an interactive `claude` session |
| Tools run slowly, one core busy | the cluster set `OMP_NUM_THREADS=1`. The launcher raises it to the job's CPU count. `nproc` also follows that variable, so it can report 1 even when more CPUs are allocated |
| The process is killed on the login node | work ran outside an allocation. Start through the launcher |
| Downloads fail on the compute node | some clusters block outbound network on compute nodes; ask your admins which partitions allow it (FASRC does) |

## See also

- [Ways to use IGVFagent](../../Guide/interfaces.md): every interface, hosted and local, including the MCP server and Claude Code.
- [LLM backends](../../Guide/llm-backends.md): running IGVFagent's own agent on Claude Code (`--backend claude_cli`) instead of an API key.
