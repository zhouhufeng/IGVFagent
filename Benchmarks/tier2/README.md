# Tier-2 cases — planning and tool selection

Each case is a question plus a **gold plan**: the tools a domain expert
considers necessary, which arguments are load-bearing, and which orderings are
dependencies. Skills are held fixed and assumed correct (that is Tier 1); what
is scored here is the decision layer.

Run with `igvfagent eval-tiers tier2 --all`.

## Schema

| Field | Meaning |
|---|---|
| `id` | case identifier |
| `question` | the natural-language request, verbatim |
| `rationale` | why this case exists — usually a real observed failure |
| `required_tools` | every one must be called |
| `any_of` | list of groups; at least one member of each group must be called |
| `required_args` | `{tool: [arg, …]}` — arguments that must be bound |
| `ordering` | `[[before, after], …]` dependency pairs |
| `forbidden_tools` | calling any of these is a failure |
| `inject_failures` | tools to force-fail when testing recovery |
| `expected_in_answer` | strings a correct answer should contain (Tier 3 uses these) |

## Scoring

- **Recall** — fraction of `required_tools` plus satisfied `any_of` groups.
  Weighted above precision: a missing essential tool invalidates the
  conclusion, an extra call mostly costs tokens.
- **Precision** — fraction of calls that were in the gold plan.
- **Argument score** — fraction of `required_args` actually bound.
- **Ordering violations** — a dependency pair observed in the wrong order.
  A pair where either tool was never called is *not* counted as a violation;
  that is a recall miss and double-penalising conflates two failures.
- **Recovery** — with `--inject-failures`, whether the run still reached a
  conclusion after a tool was forced to fail.

`any_of` exists so that a defensible alternative route is not scored as a
miss: several tools legitimately reach the same evidence.

## Why these four cases

Each encodes a failure actually observed in this system, so the suite is a
regression test rather than a hypothetical:

- `encode_experiment_identity` — the agent reported RNA-seq in K562 for a
  DNase-seq experiment in GM12878, because the tool returned only file paths.
- `read_back_report` — the agent said "18 diseases" and could not name them,
  having no tool to read the file it had just written.
- `variant_list_annotation` — a pasted variant list could not reach any tool,
  so the agent tried to author its own.
- `gene_context_basic` — the control: a single-entity lookup that should be
  one call.
