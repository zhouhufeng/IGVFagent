# Test Skills

A tiered smoke-test sweep of every IGVFagent skill. Assumes you have
already followed [`PREPARATION_SKILLS.md`](PREPARATION_SKILLS.md). Run
every command from the project root; outputs land in `Docs/<Skill>/…`
and `Data/Manifests/<Skill>/…` (always **relative** to the project
root).

## How to use this guide

- Tiers are cheapest-first. Stop after the tier that matches how heavily
  you plan to use the agent.
- Every command below is independent and safe to re-run; outputs are
  timestamped, never overwriting.
- Replace `.venv/bin/igvfagent` with plain `igvfagent` if you've
  activated the venv with `source .venv/bin/activate`.
- A run is "passing" if the command exits 0 and prints output paths
  under `Docs/…` and/or `Data/…`. Most skills also write a
  `report.md` you can open to inspect results.

## Tier 0 — Sanity

```bash
.venv/bin/igvfagent --version           # 0.1.0
.venv/bin/igvfagent --list              # full skill list
.venv/bin/igvfagent backends            # registered LLM providers
.venv/bin/igvfagent tools               # tool registry the agent runtime sees
.venv/bin/igvfagent client check        # Catalog/Portal/ENCODE ping
```

Expected: Catalog endpoints 200, Portal 403 if `IGVF_PORTAL_COOKIE`
unset. ENCODE 200.

## Tier 1 — `write-playbook` (no network, ~10 s total)

Every skill ships a `write-playbook` subcommand that regenerates its
playbook under `Docs/Skills/`. These are the cheapest possible smoke
tests for the CLI dispatcher and the per-skill module import path.

```bash
.venv/bin/igvfagent advanced-variant write-playbook
.venv/bin/igvfagent singlecell        write-playbook
.venv/bin/igvfagent multiome          write-playbook
.venv/bin/igvfagent splitseq          write-playbook
.venv/bin/igvfagent ccre              write-playbook
.venv/bin/igvfagent enhancer          write-playbook
.venv/bin/igvfagent mpra              write-playbook
.venv/bin/igvfagent crispri           write-playbook
.venv/bin/igvfagent kg                write-playbook
.venv/bin/igvfagent portal-kg         write-playbook
.venv/bin/igvfagent ref               write-playbook
.venv/bin/igvfagent specialized       write-playbook
.venv/bin/igvfagent explain           write-playbook
.venv/bin/igvfagent encode            write-playbook
.venv/bin/igvfagent se-targets        write-playbook
.venv/bin/igvfagent rnaseq            write-playbook
.venv/bin/igvfagent geo               write-playbook
.venv/bin/igvfagent proteomics        write-playbook
```

Expected: 18 lines of "Wrote / Playbook: Docs/Skills/…" output, exit 0.

## Tier 2 — Public network smokes

Hits Catalog, ENCODE, GEO, SCREEN, and literature APIs. No auth
required. ~3 min total.

```bash
# IGVF Catalog client
.venv/bin/igvfagent client gene TP53 --limit 3
.venv/bin/igvfagent client variant rs429358 --limit 3
.venv/bin/igvfagent client encode-search --type Experiment \
   --param assay_title=ATAC-seq --param limit=3

# Catalog + ENCODE smoke
.venv/bin/igvfagent data catalog-smoke --limit 3
.venv/bin/igvfagent data encode-smoke  --limit 3

# Knowledge Graph (HTTP API; no Arango password needed)
.venv/bin/igvfagent kg gene    APOE        --depth 1 --limit 5 --label smoke_apoe
.venv/bin/igvfagent kg variant rs429358    --limit 5            --label smoke_var
.venv/bin/igvfagent kg region  chr19:44903000-44912000 --limit 5 --label smoke_region

# Enhancer–gene linkage
.venv/bin/igvfagent enhancer overview  --source catalog --limit 5
.venv/bin/igvfagent enhancer pull-sets --region chr19:44903000-44912000 --gene APOE --limit 5

# MPRA / CRISPRi catalog overviews
.venv/bin/igvfagent mpra    pull --source catalog --limit 5
.venv/bin/igvfagent crispri pull --source catalog --limit 5

# cCRE linkage manifest (SCREEN page is rate-limited; this falls back gracefully)
.venv/bin/igvfagent ccre linkage-manifest --source all --limit 10 --hydrate-limit 5

# Specialized assay catalog
.venv/bin/igvfagent specialized smoke --skill all --limit 3

# Literature search (Crossref/S2 may 400/429 transiently; OpenAlex + PubMed are reliable)
.venv/bin/igvfagent ref search --query "enhancer-gene linkage rE2G" --top 5

# GEO retrieval
.venv/bin/igvfagent geo search --query "GM12878 RNA-seq" --limit 5

# ENCODE pipeline discovery
.venv/bin/igvfagent encode retrieve --assay 'Histone ChIP-seq' \
   --target H3K27ac --biosample K562 --assembly GRCh38 --limit 3

# Dataset explainer
.venv/bin/igvfagent explain explain IGVFDS2544COZH
```

Expected: each command writes a manifest under `Data/Manifests/<Skill>/`
and/or a report under `Docs/<Skill>/`.

## Tier 3 — Pipelines on the example variant list

Uses the bundled `Data/Input/VariantList/example_variants.csv`. ~2 min.

```bash
# Plain variant annotation against Catalog evidence
.venv/bin/igvfagent variant \
  --input Data/Input/VariantList/example_variants.csv --max-rows 10

# Variant annotation with cCRE / FAVOR overlay
.venv/bin/igvfagent ccre annotate-variants \
  --input Data/Input/VariantList/example_variants.csv --max-rows 10

# Integrated variant pipeline (logistic + plots + markdown report)
# Requires matplotlib/pandas/numpy/scipy from Preparation §3a.
.venv/bin/igvfagent advanced-variant run \
  --input Data/Input/VariantList/example_variants.csv --label smoke_example

# MPRA / CRISPRi joins against the same variant list
.venv/bin/igvfagent mpra    analyze-local \
  --input Data/Input/VariantList/example_variants.csv --label smoke_example
.venv/bin/igvfagent crispri analyze-local \
  --input Data/Input/VariantList/example_variants.csv --label smoke_example
```

Expected outputs:

- `Docs/VariantAnnotation/<ts>_example_variants_annotation_report.md`
- `Docs/cCRELinkage/<ts>_example_variants_favor_igvf_variant_annotation_report.md`
- `Docs/AdvancedVariantAnalysis/<ts>_smoke_example_analysis_report.md` + SVG plots
- `Docs/MPRA/<ts>_smoke_example_analysis_report.md`
- `Docs/CRISPRi/<ts>_smoke_example_analysis_report.md` + 4 SVG plots

## Tier 4 — Medium-heavy network skills

Real-but-bounded network work. ~3 min total when nothing is cached.

```bash
# Proteomics — version probe only (no downloads yet)
.venv/bin/igvfagent proteomics versions

# SPLiT-seq retrieval (Portal AnalysisSets via api.data.igvf.org)
.venv/bin/igvfagent splitseq retrieve --limit 5 --label smoke_split

# Reference skill — design + learn (pulls literature across 5 sources)
.venv/bin/igvfagent ref design --data-type parse_split_seq --assay-title "Parse SPLiT-seq"
.venv/bin/igvfagent ref learn  --topic "10x multiome human putamen" --limit 5 --top 3

# Super-enhancer pipeline — discovery half only
.venv/bin/igvfagent se-targets discover --biosample K562 --assembly GRCh38 --target H3K27ac

# Single-cell smoke against ENCODE (Portal smoke needs auth)
.venv/bin/igvfagent singlecell smoke --skill all --source encode --limit 3

# Enhancer linkage comparison (demo-if-empty so it works without cached corpora)
.venv/bin/igvfagent enhancer compare-sets --include-local-catalog --demo-if-empty

# 10x Multiome retrieval (uses api.data.igvf.org; works without Portal cookie)
.venv/bin/igvfagent multiome retrieve --count 5
```

## Tier 5 — Agent runtime (LLM-driven)

Verifies the plan→act loop end-to-end. Pick the backend that matches
your setup. The example below uses `claude_cli` because it requires
zero extra credentials if you already use Claude Code.

```bash
# No tool call (pure LLM round-trip — fast)
.venv/bin/igvfagent ask --backend claude_cli \
  "List exactly three IGVF Catalog entity types and return only the list, nothing else."

# One tool call (KG traversal end-to-end)
.venv/bin/igvfagent ask --backend claude_cli --max-iterations 3 --tool kg_gene \
  "Pull a one-paragraph summary of gene APOE from the IGVF Catalog. Use kg_gene with depth=1, limit=10."
```

Expected: each `ask` run writes `Docs/Agent/<ts>_<query-slug>/transcript.json`
+ `report.md`. The second command also creates a
`Docs/KGTraversal/<ts>_gene_APOE/` evidence pack via the tool.

For the other backends:

```bash
# Cloud Anthropic (best function-calling quality)
export ANTHROPIC_API_KEY=...
.venv/bin/igvfagent ask --backend anthropic --model claude-sonnet-4-5 "…"

# Local Ollama (free, private)
ollama serve & && ollama pull qwen3:8b
.venv/bin/igvfagent ask --backend ollama --model qwen3:8b "…"
```

## Tier 6 — Skipped by default

Not exercised by this guide because they require either a large
download, a long pipeline, or a credential the public smoke can't
provide. Run individually when you have the inputs:

- `proteomics download/build-kg/pipeline` — 100s of MB across BioGRID /
  IntAct / Reactome / KEGG. Use `proteomics download --source biogrid`
  to test one source at a time.
- `proteomics vampseq-pull / vampseq-analyze` — pulls MaveDB scoresets
  (~10 MB) and runs the deep VAMP-seq analysis.
- `encode download / browser / bigwig-* / hic-* / loops-analyze` —
  needs actual BED/bigWig/.mcool files downloaded first.
- `multiome process-local / splitseq process-local` — needs the
  matched download manifest from a prior `retrieve`/`download`.
- `rnaseq qc/pca/deg/pipeline` — needs a real counts matrix +
  sample sheet.
- `portal-kg pull/annotate/enrich` — Portal-cookie required for the
  full annotate step; `pull` alone works against
  `api.data.igvf.org`.
- `igvfagent ui` — interactive Streamlit, not a one-shot CLI test.

## Known issues / mismatches to be aware of

- `igvfagent client encode-search` does **not** accept `--limit`. Pass
  the limit via `--param limit=N` instead. (The README's quick-start
  block currently shows `--limit 3`, which fails.)
- `igvfagent geo search` uses `--limit`, not `--max-results`. The README
  shows `--max-results` in places.
- `igvfagent advanced-variant run` hard-requires matplotlib + pandas.
  If you installed only `[llm]`, you'll see
  `matplotlib is required for plots` and an exit code 1 before any
  output is written. Install Tier-3a's extra packages first.
- The IGVF Portal `data.igvf.org` site returns 403 without
  `IGVF_PORTAL_COOKIE`. The read-only `api.data.igvf.org` endpoint is
  unauthenticated and is what most skills actually hit, so most public
  workflows work without a cookie.
- Crossref / Semantic Scholar can return transient 400/429 on
  literature queries — the Reference skill tolerates this and falls
  back to OpenAlex + PubMed + arXiv.

## Run-the-whole-thing one-liner

If you want to fire every tier of this guide in order from a clean
checkout:

```bash
# Tier 0
.venv/bin/igvfagent --version && .venv/bin/igvfagent client check

# Tier 1 (write-playbook sweep)
for s in advanced-variant singlecell multiome splitseq ccre enhancer \
         mpra crispri kg portal-kg ref specialized explain encode \
         se-targets rnaseq geo proteomics; do
   .venv/bin/igvfagent "$s" write-playbook >/dev/null && echo "OK $s"
done
```

Add Tier 2/3/4 commands as needed.
