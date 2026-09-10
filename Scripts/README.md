# Scripts

All IGVF agent skills and source code live here. Every script is a CLI tool
that can be invoked from the repository root.

## Core client

```bash
python3 Scripts/igvf_client.py check
python3 Scripts/igvf_client.py catalog-api /
python3 Scripts/igvf_client.py catalog-files --limit 25
python3 Scripts/igvf_client.py gene TP53 --limit 10
python3 Scripts/igvf_client.py variant rs58658771 --limit 10
python3 Scripts/igvf_client.py encode-search --type Experiment --param assay_title=ATAC-seq
python3 Scripts/igvf_client.py aql "FOR doc IN genes LIMIT 5 RETURN doc"
```

## Catalog / Portal / ENCODE overviews

```bash
python3 Scripts/igvf_data_skills.py catalog-smoke --limit 10
python3 Scripts/igvf_data_skills.py overview --limit 25
python3 Scripts/igvf_data_skills.py encode-overview --limit 25
python3 Scripts/igvf_data_skills.py encode-smoke --limit 10
python3 Scripts/igvf_data_skills.py encode-export-csv --type Experiment --param assay_title=ATAC-seq
python3 Scripts/igvf_frontpage_summary.py refresh --update-readme
```

## Variant annotation

Provide your own variant CSV with `--input`. A small illustrative example is
included at `Data/Input/VariantList/example_variants.csv`.

```bash
python3 Scripts/annotate_variant_list.py --input Data/Input/VariantList/example_variants.csv --max-rows 10
python3 Scripts/ccre_linkage_annotation_skills.py annotate-variants --input Data/Input/VariantList/example_variants.csv --max-rows 10
```

## Advanced variant analysis (integrated functional + experimental modeling)

End-to-end pipeline that combines IGVF Catalog evidence, ENCODE cCRE class,
predicted-functional composites, optional user experimental data, logistic
models, and per-gene Miami / volcano / overlap plots into a research report.

```bash
# Annotation + composite + plots only (no experimental data)
python3 Scripts/advanced_variant_analysis.py run \
  --input Data/Input/VariantList/example_variants.csv \
  --label example_locus

# Joining a user CRISPRi/MPRA/GWAS table and modeling an outcome
python3 Scripts/advanced_variant_analysis.py run \
  --input Data/Input/VariantList/my_variants.csv \
  --experimental Data/Input/Experimental/my_crispri.csv \
  --outcome BEAN_pval_lt_.05 \
  --gene-list LDLR,PCSK9,APOE \
  --label my_crispri_v1

python3 Scripts/advanced_variant_analysis.py write-playbook
```

## Single-cell, multiome

```bash
python3 Scripts/single_cell_data_skills.py smoke --skill all --source encode --limit 5
python3 Scripts/single_cell_data_skills.py manifest --skill scrna --source both --limit 25
python3 Scripts/single_cell_data_skills.py download-examples
python3 Scripts/single_cell_data_skills.py analyze-examples --max-cells 12000
python3 Scripts/single_cell_data_skills.py write-playbook
python3 Scripts/multiome_10x_pipeline.py retrieve --count 20 --fetch-file-details
python3 Scripts/multiome_10x_pipeline.py process-local --file-manifest Data/Manifests/Multiome10x/<files.csv> --download-manifest Data/Manifests/Multiome10x/<download_manifest.csv>
python3 Scripts/multiome_10x_pipeline.py write-playbook
python3 Scripts/multiome_research_demo.py
python3 Scripts/igvf_specialized_data_skills.py smoke --skill all --limit 5
python3 Scripts/igvf_specialized_data_skills.py manifest --skill all --limit 25
python3 Scripts/igvf_specialized_data_skills.py download-plan --skill all --limit 25
python3 Scripts/igvf_specialized_data_skills.py write-playbook
```

## cCRE / linkage / FAVOR

```bash
python3 Scripts/ccre_linkage_annotation_skills.py screen-manifest
python3 Scripts/ccre_linkage_annotation_skills.py screen-download --only PLS --download --max-download-gb 1
python3 Scripts/ccre_linkage_annotation_skills.py linkage-manifest --source all --limit 100 --hydrate-limit 50
python3 Scripts/ccre_linkage_annotation_skills.py cosmic-from-favor --region chr19:44851820-44908922
python3 Scripts/ccre_linkage_annotation_skills.py browser-demo --region chr19:44850000-44910000
python3 Scripts/ccre_linkage_annotation_skills.py write-playbook
```

## Data illustration and interpretation

```bash
python3 Scripts/data_illustration_interpretation.py explain '<igvf-portal-url>/curated-sets/IGVFDS2544COZH/'
python3 Scripts/data_illustration_interpretation.py explain '<encode-portal-url>/search/?type=Annotation&searchTerm=encode-re2g&status!=archived'
python3 Scripts/data_illustration_interpretation.py explain IGVFDS2544COZH --download --max-download-gb 2
python3 Scripts/data_illustration_interpretation.py write-playbook
```

## Pathway databases (KEGG + Reactome + WikiPathways, integrated)

Pulls current releases, normalises identifiers through NCBI `gene_info`,
unifies pathways across databases and loads them into the local KG.
Full documentation: [`Docs/PATHWAYDB.md`](../Docs/PATHWAYDB.md).

```bash
igvfagent pathwaydb sources                     # databases, licences, coverage
igvfagent pathwaydb pull --kgml                 # download current releases
igvfagent pathwaydb build --relations           # normalise, unify, ingest
igvfagent pathwaydb query --genes CLU,BIN1,PICALM
igvfagent pathwaydb evaluate --agreement        # score the unification criteria
igvfagent pathwaydb status
igvfagent intpath status                        # the 2012 IntPath release
```

Licensed databases that cannot be fetched anonymously (BioCyc/HumanCyc) are
integrated from a local export:
`igvfagent pathwaydb build --extra-gmt HumanCyc=~/humancyc.gmt`

## Enhancer-gene linkage

```bash
python3 Scripts/enhancer_gene_linkage_skills.py overview --source catalog --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py overview --source encode --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py pull-sets --region chr1:903900-904900 --gene SAMD11 --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py compare-sets --include-local-catalog --demo-if-empty
python3 Scripts/enhancer_gene_linkage_skills.py write-playbook
```

## MPRA / STARR

```bash
python3 Scripts/mpra_data_skills.py pull --source catalog --limit 25
python3 Scripts/mpra_data_skills.py portal-manifest --limit 100 --label igvf_portal_mpra_many
python3 Scripts/mpra_data_skills.py analyze-local --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra
python3 Scripts/mpra_data_skills.py literature-demo --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra_literature_demo
python3 Scripts/mpra_data_skills.py write-playbook
```

## CRISPRi / CRISPR-FACS / Perturb-seq

```bash
python3 Scripts/crispri_data_skills.py pull --source catalog --limit 25
python3 Scripts/crispri_data_skills.py analyze-local --input Data/Input/VariantList/example_variants.csv --label my_locus_crispri
python3 Scripts/crispri_data_skills.py write-playbook
```

## Raw reads → count matrix → analysis

The routing front door for any dataset's raw data. `plan` reports the route,
the download size and whether an aligner is present, and downloads nothing;
`run` executes it. It REFUSES assays whose reads are not a transcript
library, naming the right tool instead — only 59.5% of the Portal's 11,070
MeasurementSets are transcript assays, so a default of "quantify it as RNA"
is wrong for 4,486 of them.

```bash
igvfagent raw-pipeline plan          IGVFDS7013XXYV      # route + cost, no download
igvfagent raw-pipeline run           IGVFDS7013XXYV --label multiome_run
igvfagent raw-pipeline run           IGVFDS7013XXYV --detach     # long jobs
igvfagent raw-pipeline status                                    # detached runs
igvfagent raw-pipeline assay-coverage                            # which assays are classified
igvfagent raw-pipeline guide-library IGVFDS6464SOVZ              # is the sgRNA library reachable?
igvfagent raw-pipeline guide-count   IGVFDS6464SOVZ              # count guides in one library
```

## IGVF Portal queries

```bash
igvfagent portal search --type MeasurementSet --limit 25
igvfagent portal get IGVFDS6464SOVZ
igvfagent portal schema MeasurementSet
igvfagent portal list-types
igvfagent portal facets --type MeasurementSet
igvfagent portal report --type MeasurementSet
```

## IGVF Catalog queries

```bash
igvfagent catalog get-entity ENSG00000141510
igvfagent catalog resolve-id rs7412
igvfagent catalog variant-evidence --variant rs7412
igvfagent catalog variant-enhancers --variant rs1250566
igvfagent catalog search-region --region chr19:44900000-45000000
igvfagent catalog find-associations --from-id ENSG00000130164
igvfagent catalog find-ld --variant rs7412
```

`variant-evidence` ranks variants by how many DISTINCT assay types support
them, which is the question "is this variant well evidenced" actually asks.

## Variant lists

```bash
igvfagent variant-list annotate --input Data/Input/VariantList/example_variants.csv
igvfagent variant-list parse    --input my_variants.txt
```

## Enrichment (ORA / GSEA)

```bash
igvfagent enrich ora      --genes TP53,MDM2,CDKN1A
igvfagent enrich gsea     --ranked my_scores.tsv
igvfagent enrich go       --genes TP53,MDM2,CDKN1A
igvfagent enrich pathways --genes TP53,MDM2,CDKN1A
igvfagent enrich showcase                              # cell-cycle positive control
```

## Single-cell effect-size / power estimation (scEPS)

```bash
igvfagent sceps estimate  --input processed.h5ad
igvfagent sceps cluster   --input processed.h5ad
igvfagent sceps aggregate
```

## Open4Gene enhancer–gene links

```bash
igvfagent open4gene link --gene PCSK9
```

## Documents (PDF / DOCX / tables)

```bash
igvfagent document read    --path Data/Uploads/<session>/paper.pdf
igvfagent document analyze --path Data/Uploads/<session>/manifest.csv
igvfagent document plan    --path Data/Uploads/<session>/manifest.csv
```

## figshare retrieval

```bash
igvfagent figshare search   --query "MPRA"
igvfagent figshare article  --id 12345678
igvfagent figshare files    --id 12345678
igvfagent figshare download --id 12345678
```

## Paper replication benchmarks

End-to-end: resolve a paper, harvest its accessions and claims, route them
to an IGVFagent chain, scaffold and run it, then score the reproduction
against the paper's stated numbers.

```bash
igvfagent bench resolve  --query "Rosen 2025 MPRA"
igvfagent bench harvest  --paper-id rosen2025
igvfagent bench route    --paper-id rosen2025
igvfagent bench scaffold --paper-id rosen2025
igvfagent bench run      --paper-id rosen2025
igvfagent bench score    --paper-id rosen2025
igvfagent bench report   --paper-id rosen2025
```

## MCP server (expose IGVFagent to other agents)

```bash
igvfagent mcp list          # tools that would be exposed
igvfagent mcp manifest      # MCP manifest JSON
igvfagent mcp serve         # run the server
```

## CRISPR screens from raw reads (three shapes, three tools)

A sorted CRISPR screen is not analysable one bin at a time: the measurement
IS the comparison between bins, so each tool finds the screen's siblings and
analyses them together. Which tool applies depends on the screen's shape and
on its guide library, and the wrong one refuses rather than undercounting.

```bash
# 1. TAIL SORT (bottom20% vs top20%). Handles both bin-naming conventions
#    (bottom20/top20 and Bot20/Top20) and the unsorted Bulk bins some add.
igvfagent crispr-screen discover IGVFDS5542IBUS
igvfagent crispr-screen analyze  IGVFDS5542IBUS --tail 20 --label ldlr

# 2. LETTERED-BIN GRADIENT (BinA..BinF, no tails). 554 Portal MeasurementSets
#    have this shape. Scores each construct's frequency-weighted mean bin.
igvfagent gradient-screen discover IGVFDS8710ZSOZ
igvfagent gradient-screen analyze  IGVFDS8710ZSOZ --label kitlg

# 3. BASE EDITING (ABE/CBE). Follows crispr-bean's masked matching, because
#    a base editor edits the guide's own locus and exact matching discards
#    those reads: 36.7% assigned vs 62.5% on IGVFDS6464SOVZ.
igvfagent bean discover IGVFDS6464SOVZ        # editor, and which BEAN inputs IGVF publishes
igvfagent bean count    IGVFDS6464SOVZ        # measure which masking recovers reads
igvfagent bean analyze  IGVFDS6464SOVZ --label ldl_abe
igvfagent bean analyze  IGVFDS5542IBUS --run-bean   # hand the counts to real BEAN
```

`--run-bean` writes BEAN's own four input tables, runs `bean create-screen`
and `bean run sorting variant`, and reports BEAN's per-target posteriors
(`mu`, `mu_sd`, `mu_z`) next to the per-guide scores. `bean discover` says up
front whether it can complete for a given screen, because two things are
measured rather than assumed:

- **An unsorted bin is required** as `--control-condition`, and it is found
  by shared `construct_library_sets` rather than by alias prefix. The
  depositor does not use one naming scheme per screen: `IGVFDS6464SOVZ`'s
  tails are `18loci_uptake_Rep*_bottom20` and its unsorted bins are
  `B*_18loci_Rep*_bulk`, so grouping on the parsed series name split one
  screen into two and lost all four. With them, `bean run` returns 1,659
  target posteriors. Where a screen truly has none, the counts and
  `bean_screen.h5ad` are still written and `bean run` is declined rather than
  handed a tail bin as its baseline.
- **BEAN's activity normalisation is unavailable on IGVF data.** Both it and
  `--scale-by-acc` route through MixtureNormal, which reads
  `screen.layers['X_bcmatch']` — counts assigned by guide barcode, which IGVF
  does not publish. The models that avoid that read reject
  `--guide-activity-col`, so BEAN supplies the posterior under its Normal
  model and `log2_per_edit` stays this tool's own activity-aware estimate.
  The output labels which model produced which number.

The counting key is chosen by measurement in all three: a prime-editing
library shares one spacer across every variant it installs, so keying on
spacer collapses 1,741 pegRNAs onto 52. Counts are cached per library, so a
re-analysis at a different `--tail` returns in seconds.

## Submitting to the IGVF Portal

The submission loop is monthly: submit, wait for DACC audits, learn at the
next meeting what was missing. The same properties are missing every time, so
those findings are encoded as checks that run in seconds. Every rule cites
the meeting that produced it.

```bash
igvfagent submit doctor                             # creds, connectivity, which lab
igvfagent submit audit IGVFDS2581EDPS --recurse     # the recurring-findings checklist
igvfagent submit inputs IGVFFI0856UJHP              # flag revoked/archived inputs
igvfagent submit validate mydata.bed.gz --bed3      # gzip, ragged rows, 0-based half-open
igvfagent submit template prediction_set            # required properties as a TSV header
igvfagent submit fileset-files IGVFDS5361GWAK       # paste-ready JSON accession array
igvfagent portal-qc audits --lab "Lin lab"          # what my lab still needs to fix
igvfagent portal-qc provenance --dead-status revoked,archived
```

`audit` exits non-zero when release is blocked, so it drops into CI.

**The revoked-input question.** When an input is corrected upstream, the only
thing you need to know is whether the correction touched your data:

```bash
igvfagent submit crosscheck --old IGVFFI7160EKDK                             --new IGVFFI1678CDBR --mine IGVFFI0856UJHP
```

It keys on `(chrom, pos) -> {(ref, alt)}`, because a reference-allele mismatch
is corrected **in place** — the coordinate does not move while ref/alt change,
so comparing positions alone reports exactly this case as "no change". If none
of the revised variants reach your file, repoint `derived_from` and release;
otherwise the file needs regenerating. It refuses to give a verdict unless
*unrevised* input variants also appear in your file, since otherwise the two
files' conventions differ and a clean result would be false.

**Rehearse on staging, not sandbox.** `api.sandbox.igvf.org` answers HTTP 410
now; `-m sandbox` warns and redirects. Add `-m staging` to any command.

**Schema requirements hide in disjunctions.** `prediction_set` has no
top-level `required` at all — `lab`, `award` and `file_set_type` are always
needed, then `samples` **or** `donors`, expressed as a `oneOf`. `template`
reads the live schema and reports both parts, so it cannot drift from what the
Portal accepts.

## Saturation genome editing (SGE)

```bash
igvfagent sge design  IGVFDS4629JYPY
igvfagent sge count   IGVFDS4629JYPY
igvfagent sge score   IGVFDS4629JYPY --label palb2
igvfagent sge analyze IGVFDS4629JYPY
```

## snMCT-seq (RNA + methylation from the same nuclei)

```bash
igvfagent mct discover IGVFDS4826YNLK
igvfagent mct analyze  IGVFDS4826YNLK --label mct_run
```

Both halves are analysed and cross-compared. Methylation ratios are NOT
log-normalised, which is what makes a methylome look like an expression
matrix.

## Catalog GRN and protein-variant effects

```bash
igvfagent grn network --regulator ENSG00000141510
igvfagent grn protein-variants --protein P04637
```

## eQTL / MPRA concordance

```bash
igvfagent eqtl-mpra tissues --genes PCSK9,LDLR
igvfagent eqtl-mpra eqtl --genes PCSK9 --tissue Liver
igvfagent eqtl-mpra concordance --genes PCSK9,LDLR --tissue Liver
```

## Pathway / interaction network for a gene list

```bash
igvfagent pathway-viz network --genes TP53,MDM2,CDKN1A --sources reactome,kegg
```

## scE2G bulk ingestion into the local KG

```bash
igvfagent sce2g-kg pull --region chr19:44900000-45000000
```

## Cross-source variant verification

```bash
igvfagent variant-verify clinvar --annotated Data/Annotated/VariantList/my_run.csv
```

Re-checks annotations against current ClinVar rather than trusting a cached
call.

## Reading back artefacts

The agent's own read path — how it inspects what a previous step wrote,
instead of re-running the analysis to see the answer.

```bash
igvfagent artifact ls   --path Docs/CrisprScreen
igvfagent artifact read --path Docs/CrisprScreen/ldl_abe/summary.json
igvfagent artifact grep --pattern "significant_fdr" --path Docs
```

## Skill cards (machine-readable skill specifications)

```bash
igvfagent skillcard list
igvfagent skillcard show --skill crispr-screen
igvfagent skillcard export
igvfagent skillcard coverage        # which skills have validation
```

## Agent-authored extensions

IGVFagent writing its own tools and skills. **An authored skill is Python
this host later executes**, so on a shared deployment this turns "knows the
password" into "can run code" — gated by `IGVF_ALLOW_AGENT_AUTHORING`, and
visitor uploads are separately gated by `IGVF_ALLOW_UPLOAD_EXTENSIONS`.

```bash
igvfagent extauthor list
igvfagent extauthor write-tool      # manifest for an existing CLI
igvfagent extauthor write-skill     # new Python skill module
igvfagent extauthor validate
igvfagent extauthor remove
```

## Evaluation tiers and playbook freezing

Three-tier evaluation — skill correctness, planning quality, conclusion
validity — and freezing a recorded session into a deterministic playbook so
the same question re-runs the same way.

```bash
igvfagent eval-tiers list
igvfagent eval-tiers tier2
igvfagent eval-tiers tier3

igvfagent playbook-freeze list
igvfagent playbook-freeze freeze
igvfagent playbook-freeze check
```

## Functional-assay calibration → ACMG/AMP evidence

```bash
python3 Scripts/excalibr_skill.py thresholds --prior 0.1
python3 Scripts/excalibr_skill.py prepare --pillar MSH2_Jia_2021.csv --name MSH2_Jia_2021
python3 Scripts/excalibr_skill.py run --pillar MSH2_Jia_2021.csv --name MSH2_Jia_2021 --components 2 3 --n-bootstraps 1000 --fits-per-bootstrap 100
python3 Scripts/excalibr_skill.py assign --calibration MSH2_Jia_2021_2c_calibration.json --scores my_variants.csv
python3 Scripts/excalibr_skill.py selftest
python3 Scripts/excalibr_skill.py write-playbook
```
