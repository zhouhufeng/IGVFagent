# Long-running jobs

Reproducing a paper, running a pipeline from raw reads or downloading and
analysing a large dataset cannot finish in one chat reply. IGVFagent runs such
tasks as **jobs**. A job keeps working in the background, plans on disk, has
each stage checked, and is reviewed by an independent verifier before it
reports done. The design follows Paper2Agent's workflow (stage markers that
depend on evidence, bounded retries, fresh verifiers, resumable state).

**On this page**

- [What changes compared with a chat reply](#what-changes-compared-with-a-chat-reply)
- [Using jobs in the browser](#using-jobs-in-the-browser)
- [Using jobs from the command line](#using-jobs-from-the-command-line)
- [Two orchestrators](#two-orchestrators)
- [How a job decides it is done](#how-a-job-decides-it-is-done)
- [Budgets and safety](#budgets-and-safety)

## What changes compared with a chat reply

| | Chat reply | Job |
|---|---|---|
| Where it runs | inside the web request | its own background process; closing the page does not stop it |
| Length | one agent run (25 iterations on the hosted site) | rounds of 40 iterations each, up to a wall-clock budget (default 4 h) |
| Plan | in the reply text | `plan.json` on disk: stages, success criteria, evidence, checks |
| "Done" | the model stops calling tools | every stage's evidence exists and its check passes, and an independent verifier approves the report |
| Failures | listed in the answer | the stage reverts to failed with the reason; up to 5 repair attempts, then it must be marked blocked with a reason |
| Continue later | starts cold | "continue" resumes from the first unfinished stage |
| After a restart | lost | the job resumes itself after a redeploy |

## Using jobs in the browser

- One **Orchestrator** choice in the sidebar's Model section applies to every message, short or long. *IGVF Agent* answers short questions directly and runs long tasks as jobs. *Claude Code agent* runs every message through Claude Code: a short question becomes a small job (3 rounds, 20 minutes) whose answer appears in the chat.
- Under **🕒 Long tasks**, *Auto* (the default) turns messages that ask to reproduce, re-run, process from raw reads or run end to end into jobs. *Always* sends every message to a job; *Never* keeps single replies.
- **🗂️ Projects & jobs** in the sidebar is where work is organised. New chats and new jobs are filed into the active project automatically.
  - The project tab lists that project's jobs (status, stages done, ⏹ stop, ▶ resume; click one to open it above the chat) and its chats (click to reopen the answer).
  - *Unfiled* lists your jobs and chats that are in no project; ➕ files one into the active project.
  - *Manage* renames, shares or creates projects.
- The jobs panel above the chat shows each of your jobs: its stages with harness-verified state, the latest events, the verifier's verdict and, when finished, the report and files. **Stop** and **Resume** buttons are there too.
- Saying **continue** (or "keep going", "resume the reproduction") resumes your latest unfinished job instead of starting over.
- The agent can also start a job itself (`job_start`) when a request turns out to need one.

## Using jobs from the command line

```bash
igvfagent job start "Reproduce Matreyek 2018 PTEN VAMP-seq from MaveDB and score it" \
    --backend anthropic --model claude-sonnet-5 --budget-minutes 60
igvfagent job list
igvfagent job status J20260924...          # plan, verdict, report
igvfagent job logs J20260924... --last 40  # every tool call and stage change
igvfagent job stop J20260924...            # after the current round (--now to kill it)
igvfagent job resume J20260924... --note "use the GRCh38 coordinates"
igvfagent job resume-interrupted           # after a restart (run by Deploy/redeploy.sh)
```

A job's files are in `Data/Jobs/<id>/`: `job.json` (state), `plan.json`,
`events.jsonl`, `rounds/NN.json` (each round's prompt and answer),
`answer.md` (the final report) and `worker.log`.

## Two orchestrators

- **IGVFagent orchestrator** (default): IGVFagent's own loop, with the job protocol and every tool the backend accepts (the whole registry on Anthropic). Inside a job it also gets:
  - `plan_set`, `plan_update` and `plan_show` for the plan;
  - `delegate_tasks`, which runs up to four sub-agents in parallel, each in a fresh context;
  - `job_wait`, which waits on a detached pipeline or an output file.
- **Claude Code agent on the Anthropic API**: Claude Code in headless mode drives the job with its own session, sub-agents and to-do list. It calls every IGVFagent tool over MCP (`igvfagent mcp serve`), and its session is resumed between rounds. The same plan checks and verifier apply, and the API cost is tracked per job. It needs the `claude` CLI and `ANTHROPIC_API_KEY`. On a shared deployment it gets IGVFagent's tools and read-only file tools but no shell; `--allow-shell` is for local installs.

## How a job decides it is done

1. The agent records a plan of 3–10 stages, each with success criteria, the evidence files it will produce and, where possible, a check.
2. When it marks a stage done, the **harness** runs the check. Available checks:
   - `files`: the files exist;
   - `json`: a value in a JSON file meets a condition (for example `r >= 0.8`);
   - `rows`: a table has at least N rows;
   - `concordance`: the benchmark scorer's checks pass.

   Missing or empty evidence, or a failed check, reverts the stage to failed, and the reason starts the next round. Accepted evidence is hashed; if it changes later, the stage must be verified again.
3. When every stage is done or blocked, a **fresh model call** that shares none of the job's context reviews the task, the plan, excerpts of the evidence and the report. If it finds unsupported claims, its issues start a repair round (at most two cycles).
4. The job ends *done*, or *done with blocked stages* when something needed data or credentials it could not get; those stages are named in the report.

## Budgets and safety

- Defaults: 240 minutes, 12 rounds, 40 iterations per round, 5 attempts per stage (`IGVF_JOB_BUDGET_MIN`, `IGVF_JOB_MAX_ROUNDS`, `IGVF_JOB_ROUND_ITERATIONS`, `IGVF_JOB_STAGE_ATTEMPTS`). A job that runs out of budget stops as *budget exhausted* and can be resumed.
- `--max-usd` caps the API cost of a Claude Code job.
- At most `IGVF_MAX_AGENT_JOBS` (default 4) jobs run at once; others wait in the queue.
- Jobs belong to their owner: on the hosted site each user sees and controls only their own jobs.
- Stage evidence must be inside the workspace; secret paths never count as evidence.
- `Deploy/redeploy.sh` waits for running jobs, and resumes any job a restart interrupted.

---

[← Documentation index](README.md) · [Project README](../../README.md)
