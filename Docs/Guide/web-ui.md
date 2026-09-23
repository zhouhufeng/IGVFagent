# The browser UI, accounts and history

What each tab of the web UI does, how sign-in works on the hosted site, and how past results and projects are kept.

**On this page**

- [The browser UI](#the-browser-ui)
- [Accounts, approval and sign-in](#accounts-approval-and-sign-in)
- [Projects and permanent history](#projects-and-permanent-history)

## The browser UI

The Streamlit UI exposes the **same** ReAct agent as `igvfagent ask`,
with an interactive sidebar for backend / model / tool-subset selection,
a streaming progress trace as the agent plans and calls tools, and inline
rendering of any artefacts the tools produce.

The page has four tabs, and **🐛 Report a bug** sits in the sidebar so it
is visible from all of them:

| Tab | Sub-tabs | What it is for |
|---|---|---|
| 💬 Chat | | Ask the agent; answers, figures and live jobs. |
| 🕸 Knowledge & networks | Explore · Your networks · Sources | Search the integrated KG (`Data/KG/local_kg.sqlite`) and walk a node's neighbourhood; view subnetworks inferred from your data (`network viz` runs); **Sources** shows every Catalog collection side by side: mirrored rows vs Catalog documents, and whether a merge step has brought it into the graph (operators and admins can start `kg-integrate` steps from there). |
| 🔬 Data viewers | Runs and files · Single-cell · Spatial-ATAC-Hi-C | One browser over every run directory (`Docs/<Skill>/<timestamp>_<label>/`) and every Portal fetch (`Data/Processed/<accession>/`), routing `.h5ad` files and Spatial-ATAC-Hi-C / network runs to their viewers. Signed-in users see only runs from their own or shared sessions. |
| 📊 Validation | | Reproducibility benchmark figures, suite dashboard, and an About panel (build id, upstream pins). |

Supported artefact viewers (rendered directly in the chat, no terminal
round-trip):

- **PNG / JPG / SVG / GIF** — inline gallery (UMAP, dot-plot, volcano, IGV-style snapshots).
- **CSV / TSV** — first 400 rows as a sortable `st.dataframe`, plus a download button.
- **PDF** — embedded base64 iframe (FAVOR / cCRE reports, advanced-variant exports).
- **Markdown reports** — rendered in place; backtick-quoted file paths inside the body are auto-followed so the underlying CSV / JSONL / PDF / PNG that the report references each get their own inline viewer.
- **JSONL** — first 50 records normalized into a DataFrame, falls back to raw JSON.
- **JSON** — pretty-printed up to 200KB.
- **HTML** — sandboxed `st.components.v1.html` embed.
- **TXT / LOG** — first 50KB in a code block.
- Anything else — download button.

When the **Claude Code CLI** backend is selected, the model picker is
restricted to the current Claude tiers (`claude-opus-5-5`, the default,
`claude-sonnet-5`, `claude-fable-5-1`) plus a `(custom...)`
escape hatch — picking a retired model id is no longer possible from
the dropdown.

## Accounts, approval and sign-in

The hosted deployment has real accounts. Identity comes from the **Genohub
Discourse community** — signup, email verification and moderation already live
there — and access to the agent requires membership of a Discourse group that
administrators control. Signing up and being approved are two separate steps on
purpose: anyone may join the forum; only approved members may drive an agent
that runs analysis pipelines on a shared machine.

Approving someone is adding them to the group. Revoking is removing them.
There is no second password to distribute and no credential store of our own —
authentication sits in the gateway (`Deploy/auth/gate.py`, stdlib only), never
in the app, because the app is the thing being protected.

Setup, cutover and day-to-day administration: [`Docs/AUTH.md`](../AUTH.md).

A local install has no gateway and therefore no login — one person at the
keyboard, nothing hidden, exactly as before.

## Projects and permanent history

The knowledge graph remembers *what is true*. A separate store remembers *what
was done*: `Data/History/history.sqlite` holds one durable row per agent run —
the question, the answer, the artefact paths, the accessions — plus every skill
invocation and download, all in one FTS5 search index.

Two things follow from that.

**A dataset is analysed once.** Ask about an accession that has been worked on
before and the agent looks it up instead of recomputing it:

```bash
igvfagent project recall IGVFDS5414UFNC   # every result ever produced about it
igvfagent project search "spatial ATAC Hi-C"   # full text over questions AND answers
igvfagent project show Docs/Agent/20260914_150552_…   # replay one session
```

The agent has the same three as tools (`history_recall`, `history_search`,
`history_show`) and its first instruction is to use them before re-running work.

**Results can be grouped into a project that outlives the session.**

```bash
igvfagent project create "Wang 2026 reproduction" --use
igvfagent project add --kind dataset --ref IGVFDS5414UFNC --title "primary set"
igvfagent project items
igvfagent project rename "Wang 2026 reproduction" --to "Spatial ATAC-Hi-C repro"
```

While a project is active, every answer is filed into it automatically — in the
CLI and in the web UI, which has a **🗂️ Project** panel in the sidebar.
Renaming is safe: items reference the project's immutable id, and every former
name stays resolvable, so a reference written down months ago still works.

**Results are shared, organisation is private.** Any past answer is readable
by anyone signed in and is recalled automatically when someone asks about the
same accession — paying twice for the same analysis to hide it from a colleague
helps nobody, and each run still records who produced it. Projects are the
private part: visible to their owner and whoever they share them with
(`igvfagent project share <username>`), where members can add but only the
owner can rename, archive or change membership. `IGVF_HISTORY_SHARED=0` flips
answers back to strict per-user privacy. Without authentication in front,
nothing is filtered at all.

**Nothing here is ever deleted.** That is enforced by `BEFORE DELETE` triggers
on every history table, not by convention — removing an item from a project
marks it removed and keeps the row searchable, and archiving a project hides it
from the listing and is reversible. The store is a separate SQLite file from
the knowledge graph on purpose: the KG is bulk-rebuilt by `kg-integrate` and is
disposable scratch, while history must outlive all of it.

Index work that predates the store (safe to re-run; also rebuilds the index):

```bash
igvfagent project backfill
igvfagent project stats
```

---

[← Documentation index](README.md) · [Project README](../../README.md)
