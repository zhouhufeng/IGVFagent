# Extending IGVFagent

Add your own tools, skills, prompt skills and playbooks without editing the core code.

**On this page**

- [Add a custom tool (YAML manifest — no code)](#add-a-custom-tool-yaml-manifest--no-code)
- [Add a custom skill (Python subcommand)](#add-a-custom-skill-python-subcommand)
- [Add a prompt skill or playbook](#add-a-prompt-skill-or-playbook)

IGVFagent is **expandable by design**: the 118 built-in skills and 603
registered tools are a starting point, not a ceiling. A user-extension
framework absorbs **your own skills and tools** at startup — no core-code
edits, no re-install, no registration step. Anything you drop into an
extension directory becomes a first-class citizen: it shows up in
`igvfagent --help` and `igvfagent tools`, the `igvfagent ask` ReAct agent
can plan with and call it, the Streamlit UI lists it in the tool picker,
and its outputs are harvested into the local knowledge graph exactly like
built-in results.

There are four extension surfaces, from zero-code to full-code:

| Surface | You write | Where it plugs in |
|---|---|---|
| **Custom tool** | one YAML/JSON manifest (no code) | LLM-callable tool in `ask` + UI, wrapping *any* executable or an existing `igvfagent` subcommand |
| **Custom skill** | one Python file with `main()` | first-class `igvfagent <name>` subcommand with the post-run KG harvest |
| **Prompt skill** | one `SKILL.md` file | Claude Code slash-skill orchestrating existing CLI steps |
| **Playbook** | one YAML file | deterministic multi-step tool chain, run via `igvfagent playbook` |

Extensions are discovered from these locations (scanned in order; first
definition of a name wins, and built-in names can never be shadowed):

1. every directory in `$IGVF_USER_EXT_DIR` (`:`-separated) — for testing
   or shared lab locations
2. `~/.igvfagent/` — per-user, works from any checkout
3. `<repo>/UserExtensions/` — per-checkout, committable so a whole lab
   shares one extension set through git

Each location uses the same two subfolders: `tools/` (manifests) and
`skills/` (Python modules). Copy-paste-ready templates live in
[`Docs/Examples/user_extensions/`](../Examples/user_extensions), and

```bash
igvfagent extensions          # what was discovered, from where, + any
igvfagent extensions --json   # skipped/malformed definitions explained
```

shows exactly what the framework picked up. Discovery is defensive: a
malformed manifest or a skill that fails to import is skipped with a
diagnostic in `igvfagent extensions` and can never break the core CLI or
agent runtime.

## Add a custom tool (YAML manifest — no code)

A *tool* is what the LLM agent calls during `igvfagent ask`. To teach the
agent a new capability, describe your command in a manifest — here a lab
script that computes GC content:

```bash
mkdir -p ~/.igvfagent/tools
cat > ~/.igvfagent/tools/gc_content.yaml <<'EOF'
name: gc_content
description: >
  Compute per-sequence GC content for a FASTA file and write a TSV
  summary. Use when the user asks for GC% or base-composition QC.

# argv of YOUR program — any language, any location
command: ["python3", "/home/me/bin/gc_content.py"]

parameters:               # JSON Schema for what the LLM may pass
  type: object
  properties:
    fasta:  {type: string, description: Path to the input FASTA.}
    window: {type: integer, default: 0}
  required: [fasta]

positional: [fasta]        # emitted as a positional argument
flag_map:
  window: "--window"       # emitted as `--window <value>`
EOF

igvfagent extensions               # confirm it was absorbed
igvfagent tools | grep -A3 gc_content
igvfagent ask "run GC content QC on Data/my_seqs.fa"
```

The only contract your program must honour: read arguments from argv,
print results to stdout, exit 0 on success, and announce output files as
`Report: <path>` / `Manifest: <path>` lines — the same convention every
built-in follows, which is how the agent chains your artefacts into its
next step and how the localstore harvester grows the local KG from them.

Instead of `command:`, use `cli:` to re-surface an existing `igvfagent`
subcommand under your own name, defaults, and description (see the
`kg_gene_quick.yaml` template) — useful for lab-specific shortcuts.
Optional mapping fields mirror the built-in registry: `positional`,
`flag_map` (parameter → flag), `flag_repeat` (list values repeat the
flag), `bool_flags` (bare flags). One manifest may also carry several
definitions under a top-level `tools:` list.

## Add a custom skill (Python subcommand)

A *skill* is a human-facing subcommand. Drop a Python file with a
`main()` into `skills/` and it becomes `igvfagent <name>` (filename stem,
`_` → `-`):

```bash
mkdir -p ~/.igvfagent/skills
cp Docs/Examples/user_extensions/skills/variant_bed_export.py ~/.igvfagent/skills/

igvfagent --help                    # now lists it under "User skills"
igvfagent variant-bed-export my_variants.tsv --out variants.bed
```

Conventions worth copying from the template:

- the module docstring's first line is the description shown in
  `igvfagent --help` and `igvfagent extensions`;
- `main()` owns its own argparse parser and returns an int exit code —
  the dispatcher hands over `sys.argv` exactly as for built-ins;
- print every artefact as `Report: <path>` so downstream steps and the
  local-KG harvest (which runs after user skills too) pick it up.

To expose the same logic to the LLM agent as well, add a small `cli:`- or
`command:`-style manifest next to it — the skill/tool split is the same
one the built-ins use (human CLI surface vs. curated LLM tool surface).

## Add a prompt skill or playbook

The remaining two surfaces need no Python at all:

- **Prompt skills** — drop a `.claude/skills/<name>/SKILL.md` (YAML
  frontmatter + a markdown workflow whose steps are `igvfagent …`
  commands) and Claude Code picks it up with zero registration. The
  eight shipped `igvf-*` skills under `.claude/skills/` are the
  reference implementations; see
  `Docs/Skills/PROMPT_SKILLS_INDEX.md`.
- **Playbooks** — a YAML file listing `steps: [{tool, args}]` over any
  registered tools (including your user tools), with `${param}`
  interpolation and `--param key=value` overrides. Put it in
  `Docs/Playbooks/` and run `igvfagent playbook <name>`; the schema is
  documented at the top of `Scripts/_playbook.py`.

Together the four surfaces mean IGVFagent can absorb a lab's entire
private toolbox — QC scripts, internal APIs, bespoke pipelines — while
keeping the audited, reproducible execution model of the built-ins.

---

[← Documentation index](README.md) · [Project README](../../README.md)
