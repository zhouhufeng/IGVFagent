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
```

The counting key is chosen by measurement in all three: a prime-editing
library shares one spacer across every variant it installs, so keying on
spacer collapses 1,741 pegRNAs onto 52. Counts are cached per library, so a
re-analysis at a different `--tail` returns in seconds.

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

## Functional-assay calibration → ACMG/AMP evidence

```bash
python3 Scripts/excalibr_skill.py thresholds --prior 0.1
python3 Scripts/excalibr_skill.py prepare --pillar MSH2_Jia_2021.csv --name MSH2_Jia_2021
python3 Scripts/excalibr_skill.py run --pillar MSH2_Jia_2021.csv --name MSH2_Jia_2021 --components 2 3 --n-bootstraps 1000 --fits-per-bootstrap 100
python3 Scripts/excalibr_skill.py assign --calibration MSH2_Jia_2021_2c_calibration.json --scores my_variants.csv
python3 Scripts/excalibr_skill.py selftest
python3 Scripts/excalibr_skill.py write-playbook
```
