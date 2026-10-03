---
name: igvf-portal-facet-filter
description: Iteratively narrow IGVF Portal items by their facet histograms, never dumping the full facet payload at once. Use when a user asks "what's in IGVF for ...", "show me the breakdown of X", or wants to drill from a coarse type down to a concrete record set.
argument-hint: <ItemType>
---

# IGVF Portal facet-driven filter loop

A progressive-disclosure workflow that walks the user from "I want to find some X" down to a concrete, faceted record set, without ever flooding the chat with every facet's full term distribution. Re-derived under Apache-2.0 from the workflow described in [IGVF-DACC/igvf-portal-mcp](https://github.com/IGVF-DACC/igvf-portal-mcp) (MIT) and retargeted at IGVFagent's `portal` CLI.

## The four-step loop

### Step 1 — Top-level summary (count + facet menu only)

Start by counting the total items of the requested type and listing **only the names of the available facets** — not their term values. This keeps the first response short.

Run:

```bash
igvfagent portal facets --type "$ARGUMENTS" --label step1
```

Report:
- the total item count;
- the **list of facet titles + field names** that have more than one value (skip facets with 0 or 1 terms — they cannot discriminate);
- one sentence reminding the user that **any documented field** is filterable, not just the facet fields surfaced here.

Ask the user: which facet do they want to drill into, or do they want the full filterable-field list?

### Step 2 — Expand only what the user asks for

When the user picks a facet name, re-run `portal facets` and show the **top terms + counts** *only for the chosen facets*. Do not expand the others.

If the user instead asks for "what else can I filter on?", run

```bash
igvfagent portal endpoint-params <collection>
```

(where `<collection>` is the snake-case form of the item type — e.g. `measurement-sets` for `MeasurementSet`) and present the full `agent_param ↔ search_field` table. Mention that the dotted `search_field` form is what goes into `--field-filters`.

### Step 3 — Apply a filter and re-count

When the user picks a filter value, re-run `portal facets` with the chosen `--field-filters` accumulated into a `;`-separated DSL clause:

```bash
igvfagent portal facets --type "$ARGUMENTS" \
                          --field-filters "<accumulated_clauses>" \
                          --label step3
```

Show:
- the new total count;
- the **updated facet menu** (titles + field names only — same shape as Step 1);
- a brief note of which clauses are currently active.

Loop back to Step 2 until the user is satisfied.

### Step 4 — Materialise the result set

Once filters look right, fetch the actual records (or a TSV report, depending on volume):

```bash
# JSON for a small set (< 1k items)
igvfagent portal search --type "$ARGUMENTS" \
                          --field-filters "<accumulated_clauses>" \
                          --limit all

# TSV for a large set (everything via /report.tsv)
igvfagent portal report --type "$ARGUMENTS" \
                          --field-filters "<accumulated_clauses>" \
                          --limit all \
                          --label final
```

For FileSet types (`MeasurementSet` / `AnalysisSet` / etc.), offer to fetch the pre-signed download manifest:

```bash
igvfagent portal batch-download --type "$ARGUMENTS" \
                                  --field-filters "<accumulated_clauses>" \
                                  --label final
# add --fetch to download every referenced file
```

## Filter-DSL refresher

`--field-filters` clauses are `;`-joined; commas separate list values; `!=` is negation; `gte:`/`lte:`/`gt:`/`lt:` are range ops:

```
preferred_assay_titles=10x multiome,SHARE-seq;status=released;file_size=gte:1000000;lab.@id!=/labs/x
```

## Output etiquette

- Never dump every facet's term distribution at once — that's the whole point of this loop.
- Skip single-value facets in every iteration; they convey no choice.
- Round counts > 10k to thousands ("12.4 K results") for readability.
- If the user provides no `$ARGUMENTS`, ask which item type they want and remind them they can run `igvfagent portal list-types` to see the canonical 50+ choices.
