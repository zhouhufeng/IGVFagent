# User-extension templates

Copy-paste-ready templates for plugging your own **tools** and **skills**
into IGVFagent — no core-code edits, no re-install. Full tutorial:
[README.md → Extending IGVFagent](../../../README.md#extending-igvfagent).

## Install

Copy any template into one of the auto-discovered locations:

```bash
# per-user (works from any checkout)
mkdir -p ~/.igvfagent/tools ~/.igvfagent/skills
cp tools/gc_content.yaml        ~/.igvfagent/tools/
cp skills/variant_bed_export.py ~/.igvfagent/skills/

# or per-checkout (committable alongside the repo)
mkdir -p UserExtensions/tools UserExtensions/skills
```

Then verify — no restart or registration step needed:

```bash
igvfagent extensions          # what was discovered, from where, any problems
igvfagent tools | grep -A3 gc_content
igvfagent variant-bed-export --help
```

## Files

| Template | Kind | Shows |
|---|---|---|
| `tools/gc_content.yaml` | tool (`command:`) | wrapping **any executable on your machine** as an agent tool |
| `tools/kg_gene_quick.yaml` | tool (`cli:`) | re-surfacing an `igvfagent` subcommand with your own defaults / description |
| `skills/variant_bed_export.py` | skill | a full `igvfagent <name>` subcommand with argparse + `Report:` artefact contract |
