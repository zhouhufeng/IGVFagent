# Claude Code prompt-skill suite

In addition to the per-skill playbooks in this directory (one per CLI
subsystem), IGVFagent ships **seven workflow prompt skills** under
`.claude/skills/<name>/SKILL.md`. These are designed to be auto-loaded
by Claude Code / any MCP-aware client when the user asks a matching
question, and they orchestrate IGVFagent's `portal` + `catalog`
subcommands behind the scenes.

Clean-room reimplementations of the SKILL.md prompts published in
[IGVF-DACC/igvf-portal-mcp](https://github.com/IGVF-DACC/igvf-portal-mcp)
and [IGVF-DACC/igvf-catalog-mcp](https://github.com/IGVF-DACC/igvf-catalog-mcp)
(both MIT). The upstream prompts orchestrate the DACC MCP servers'
tool names (`igvf_portal_facets`, `igvf_catalog_find_associations`,
etc.); ours orchestrate IGVFagent's CLI surface (`igvfagent portal
facets`, `igvfagent catalog find-associations`, etc.).

## The seven skills

| Skill | Trigger question | Underlying CLI |
|---|---|---|
| `igvf-portal-facet-filter` | "what's in IGVF for X" / "show me the breakdown" | `portal facets / endpoint-params / search / report / batch-download` |
| `igvf-catalog-variant-report` | "interpret rs…" / "what does this variant do?" | `catalog resolve-id / find-associations / find-ld` |
| `igvf-catalog-gene-dossier` | "tell me about gene X" / "build a dossier for BRCA1" | `catalog get-entity / find-associations / search-region` |
| `igvf-catalog-dissect-locus` | "which gene does this GWAS hit affect?" / "fine-map this locus" | `catalog resolve-id / find-associations / find-ld / search-region` |
| `igvf-catalog-regulatory-landscape` | "how is X regulated?" / "regulatory architecture of …" | `catalog get-entity / search-region / find-associations` |
| `igvf-catalog-disease-genes` | "what genes cause X?" / "genetic basis of …" | `catalog get-entity / find-associations` |
| `igvf-catalog-ld-compare` | "compare LD across EUR/AFR/EAS/SAS/AMR" / "PRS portability for this SNP" | `catalog resolve-id / find-ld` (×5 ancestries) |

## How they fit together

The skills cross-reference each other through their final "Pairs well
with" sections. A typical multi-step interaction:

```
       (user) "fine-map rs429358 for Alzheimer's"
          ↓
   igvf-catalog-variant-report      (deep dive on the lead)
          ↓
   igvf-catalog-dissect-locus       (rank candidate causal genes)
          ↓
   igvf-catalog-gene-dossier APOE   (build dossier for the nominee)
          ↓
   igvf-catalog-regulatory-landscape APOE microglia
                                    (mechanism in the disease-relevant tissue)
          ↓
   igvf-catalog-disease-genes "Alzheimer's disease"
                                    (place APOE in the wider gene network)
```

Plus the cross-ancestry option:

```
   igvf-catalog-ld-compare rs429358     (if EUR fine-mapping is ambiguous)
```

## How Claude Code picks them up

Each `.claude/skills/<name>/SKILL.md` carries a YAML frontmatter with
`name`, `description`, and `argument-hint`. Claude Code surfaces these
as slash commands and auto-triggers them when the user's prompt
matches the `description`. No additional registration is needed once
the repo is opened in Claude Code.

The corresponding bare CLI commands work everywhere — you can invoke
the same workflows from a plain shell without Claude Code involved at
all. For example the `gene-dossier` skill is essentially a structured
sequence of `igvfagent catalog ...` calls; the SKILL.md just makes the
sequence prompt-driven rather than scripted.

## License posture

Apache-2.0 IGVFagent ⊃ MIT upstream prompts. The MCP server source code
and prompt text are licensed permissively, so verbatim copy would be
allowed with attribution — but we paraphrased each prompt and
retargeted at IGVFagent's CLI surface (different tool names, different
defaults, additional cross-pointers to IGVFagent-only skills like
`network steiner`, `enrich pathways`, `ccre`, and `enhancer`).
