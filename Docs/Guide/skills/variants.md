# Variants and variant effects

Variant annotation, advanced variant analysis, cCRE and FAVOR, regulatory-region annotation, MaveDB, saturation genome editing and assay calibration to ACMG/AMP evidence.

**On this page**

- [Variant annotation](#variant-annotation)
- [Advanced variant analysis](#advanced-variant-analysis)
- [Enhancer / regulatory-region cCRE annotation](#enhancer--regulatory-region-ccre-annotation)
- [Saturation genome editing (SGE)](#saturation-genome-editing-sge)
- [cCRE, FAVOR, IGV-style browser views](#ccre-favor-igv-style-browser-views)
- [MaveDB mapping (incl. SGE cDNA path)](#mavedb-mapping-incl-sge-cdna-path)
- [Assay calibration → ACMG/AMP evidence (exCALIBR)](#assay-calibration--acmgamp-evidence-excalibr)
- [Variant lists and your own data](#variant-lists-and-your-own-data)

## Variant annotation

Annotate any user variant CSV against IGVF Catalog evidence (CADD, QTL,
phenotypes, regulatory elements, predictions). Provide your own variant list;
a tiny illustrative example is shipped at
`Data/Input/VariantList/example_variants.csv`.

```bash
python3 Scripts/annotate_variant_list.py \
  --input Data/Input/VariantList/example_variants.csv \
  --max-rows 10
python3 Scripts/ccre_linkage_annotation_skills.py annotate-variants \
  --input Data/Input/VariantList/example_variants.csv \
  --max-rows 10
```

## Advanced variant analysis

End-to-end pipeline that combines IGVF Catalog evidence, ENCODE cCRE class,
predicted-functional composites, optional user experimental data, logistic
models, and per-gene Miami / volcano / overlap plots into a research report.
See `Docs/Skills/ADVANCED_VARIANT_ANALYSIS_SKILLS.md` for full details.

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

Outputs land in `Docs/AdvancedVariantAnalysis/` (annotated CSV, summary stats,
logistic model JSON, markdown report) and `Docs/AdvancedVariantAnalysis/Plots/`
(SVG plots).

## Enhancer / regulatory-region cCRE annotation

Upload a list of **regions** and annotate every interval against the ENCODE
SCREEN candidate cis-regulatory element registry. This is interval overlap —
distinct from `ccre annotate-variants`, which is point-level.

```bash
# Once: download + index the registry (V4 default, 2,348,854 elements)
igvfagent enhancer-annot build-db --registry V4
igvfagent enhancer-annot db-stats

# Which sheet of the workbook actually holds coordinates?
igvfagent enhancer-annot inspect --input Data/Input/EnhancerList/library.xlsx

# Annotate
igvfagent enhancer-annot annotate \
    --input Data/Input/EnhancerList/library.xlsx \
    --sheet High.sgRNA0820 --label my_library
```

Accepts `.xlsx` / `.xlsm` (any sheet), BED, TSV, CSV, gzipped. Coordinate
columns are auto-detected and a headerless BED is recognised by shape; `.xlsx`
is read with `zipfile` + `xml.etree` so no spreadsheet dependency is needed and
million-row sheets stream. A CRISPR library repeats each enhancer once per
sgRNA, so rows are deduplicated on coordinates by default.

Per region you get: number of overlapping cCREs, a ranked representative class
(PLS > pELS > dELS > CA-H3K4me3 > CA-CTCF > CA-TF > CA > TF), every class seen,
bases covered as a **union** (overlapping elements are not double-counted),
coverage fraction, and — when nothing overlaps — the nearest element and its
distance, which is what separates a cCRE desert from a near miss.

Results are written to `Data/Output/EnhancerAnnotation/<timestamp>_<label>/`
as `annotated_regions.csv`, `regions_without_ccre.csv`, `annotated_regions.bed`
and `summary.json`.

> V3 and V4 use **different class vocabularies** (V3: `DNase-H3K4me3` /
> `CTCF-only` / `CTCF-bound`; V4: `CA-H3K4me3` / `CA-CTCF` / `CA-TF` / `CA` /
> `TF`) and V4 has more than twice as many elements, so annotations are not
> comparable across versions. The registry version is recorded in every
> `summary.json`.

## Saturation genome editing (SGE)

```bash
igvfagent sge design  IGVFDS4629JYPY
igvfagent sge count   IGVFDS4629JYPY
igvfagent sge score   IGVFDS4629JYPY --label palb2
igvfagent sge analyze IGVFDS4629JYPY
```

SGE reads are a fixed amplicon, so `raw-pipeline` refuses them for a stated
reason: quantifying them against a transcriptome would restate the amplicon
design, not measure anything.

## cCRE, FAVOR, IGV-style browser views

```bash
python3 Scripts/ccre_linkage_annotation_skills.py screen-manifest
python3 Scripts/ccre_linkage_annotation_skills.py screen-download \
  --only PLS --download --max-download-gb 1
python3 Scripts/ccre_linkage_annotation_skills.py linkage-manifest \
  --source all --limit 100 --hydrate-limit 50
python3 Scripts/ccre_linkage_annotation_skills.py linkage-download \
  --manifest Data/Manifests/cCRELinkage/<manifest.csv> --only rE2G --download
python3 Scripts/ccre_linkage_annotation_skills.py cosmic-from-favor \
  --region chr19:44851820-44908922
python3 Scripts/ccre_linkage_annotation_skills.py browser-demo \
  --region chr19:44850000-44910000
python3 Scripts/ccre_linkage_annotation_skills.py write-playbook
```

Always run a `*-manifest` command before `*-download`. Full SCREEN cCRE,
rE2G, and single-cell linkage corpora can be many gigabytes.

## MaveDB mapping (incl. SGE cDNA path)

Map [MaveDB](https://www.mavedb.org) multiplexed-assay scoresets to genomic
coordinates via the Ensembl REST API (no UTA / SeqRepo / BLAT dependency). In
addition to the protein-coordinate VAMP-seq path, a dedicated **SGE (Saturation
Genome Editing) path** parses the full HGVS-c grammar used by SGE scoresets
(CDS, intronic, 5′UTR, 3′UTR) and emits VCF-4.2 with score-bearing INFO fields.
This is the path the **Waters 2024 BAP1** and **Buckley 2024 VHL** benchmarks
exercise. Playbook: [`Docs/Skills/MAVEDB_MAPPING_SKILL.md`](../../Skills/MAVEDB_MAPPING_SKILL.md).

```bash
igvfagent mavedb map-scoreset --urn urn:mavedb:00000097-0-1   # PTEN VAMP-seq
igvfagent mavedb map-scoreset --urn <BAP1-SGE-urn> --sge      # SGE cDNA path
```

## Assay calibration → ACMG/AMP evidence (exCALIBR)

Turn a raw multiplexed-assay score into a **clinically usable evidence
strength**. This is the last mile of the MAVE chain: `mavedb map-scoreset`
gives a variant genomic coordinates, and `calibrate` says how much a given
assay score is actually worth as PS3 / BS3 evidence — supporting, moderate,
strong or very strong — instead of leaving the reader with a bare number.

Clean-room reimplementation of [exCALIBR](https://github.com/rosstewart/exCALIBR)
(MIT), which implements the gene-based calibration method of Zeiberg et al.
(*bioRxiv* 2025.04.29.651326) on top of Tavtigian's Bayesian reading of the
ACMG/AMP guidelines. The chain: label variants into P/LP, B/LB, gnomAD-population
and synonymous samples → fit a **constrained skew-normal mixture** by EM (shared
components, per-sample weights, monotone density-ratio constraint enforced by
binary search inside every M-step) → **bootstrap** it → EM-estimate the prior
P(pathogenic | population) → build LR⁺(score) → solve for **Tavtigian's C** →
emit the score window that earns each evidence strength. Playbook:
[`Docs/Skills/ASSAY_CALIBRATION_SKILL.md`](../../Skills/ASSAY_CALIBRATION_SKILL.md).

```bash
# 0) Evidence thresholds alone — instant, no data needed
igvfagent calibrate thresholds --prior 0.1
#    -> C = 348; LR+ 2.08 (supporting) / 4.32 (moderate) / 18.7 (strong) / 348 (very strong)

# 1) Label a scoreset into the four calibration samples
igvfagent calibrate prepare --pillar MSH2_Jia_2021.csv --name MSH2_Jia_2021
#    or from IGVFagent's own MAVE chain, joined to a ClinVar release:
igvfagent calibrate prepare --mapped Docs/MaveDB/<run>/mapped.tsv \
    --clinvar-tsv variant_summary.txt.gz --gnomad-tsv gnomad_sites.tsv

# 2) Calibrate (long job: progress heartbeat + resumable ledger,
#    2c-vs-3c model selection, calibration JSON + figure)
igvfagent calibrate run --pillar MSH2_Jia_2021.csv --name MSH2_Jia_2021 \
    --components 2 3 --n-bootstraps 1000 --fits-per-bootstrap 100
igvfagent calibrate run --table scores.csv --name MSH2 --resume   # continue

# 3) Interpret new variants with the calibration
igvfagent calibrate assign --calibration MSH2_Jia_2021_2c_calibration.json \
    --scores my_variants.csv
#    -> score, evidence_points, acmg_evidence ("PS3 moderate", "BS3 strong", …)

# 4) Validate the port on this machine (~1 min)
igvfagent calibrate selftest
```

Every numeric component of the default path was checked against the upstream
implementation:
Tavtigian's C matches over a 12-prior grid (including the `original` and
`strict` rule variants), the constrained-EM iterates are **bit-identical** for
60 consecutive steps at 2 and 3 components, and the prior EM, LR⁺ → point-range
conversion, and paired model-selection test all agree to machine precision. The
Pillar/IGVF-format loader reproduces upstream's variant labelling exactly
(identical variant-ID sets per sample across ClinVar releases and star
thresholds). On the MSH2 (Jia 2021) example both implementations select the same
model and the same set of evidence strengths, with thresholds inside ~1–2 % of
the score range at equal bootstrap budgets. Unlike upstream the run is fully
seeded, so a rerun reproduces the calibration exactly.

## Variant lists and your own data

The variant-annotation skills accept any user-provided CSV via `--input`.
Recognized identifier columns (case-insensitive):

- `rsid` / `rsID` / `dbSNP` — e.g. `rs58658771`
- `chrom`, `pos`, `ref`, `alt` — GRCh38 coordinates and alleles
- `spdi` — NCBI SPDI string
- `hgvs` — e.g. `NC_000019.10:g.44908822C>T`

Additional user columns (locus name, phenotype, prior, notes, experimental
effect, p-value) are preserved through annotation and used by
`advanced_variant_analysis.py` when you pass `--experimental`/`--outcome`.

A minimal `Data/Input/VariantList/example_variants.csv` is included for
smoke testing. Replace it with your own list — the `.gitignore` excludes
everything in this folder except the README and the example CSV.

**Do not commit confidential or pre-publication variant lists.**

## AlphaGenome predictions and Atlas scores (`alphagenome`)

Access to [AlphaGenome](https://github.com/google-deepmind/alphagenome), Google DeepMind's sequence-to-function model, from inside IGVFagent. You can predict expression, accessibility, histone and TF binding, splicing and 3D contacts over a region or gene. You can compare REF and ALT for a variant, score variant effects with the recommended scorers, run in silico mutagenesis, or look up the pre-computed **AlphaGenome Atlas** scores (including the AlphaGenome Variant Impact score) without running the model. Variants can be rsIDs, SPDI, `chr-pos-ref-alt` or `chr:pos:ref>alt`; rsIDs and gene symbols are resolved through the IGVF Catalog, and tissues are chosen by ontology CURIE (`alphagenome metadata --search heart` lists them).

**Access.** An API key is required ([get one](https://alphagenome.google/api); free for non-commercial use). Set `ALPHAGENOME_API_KEY`, or save it to `Docs/Secret/ALPHAGENOME_API_KEY.txt`; it is never printed. Outputs are for non-commercial research only, must not be used to train other models, and are not for clinical use ([terms](https://alphagenome.google/terms)). The SDK needs Python 3.10 or newer: on a 3.9 install, `igvfagent alphagenome setup --install` creates a separate environment once, and every AlphaGenome command then runs there automatically. The hosted container (Python 3.11) has the SDK built in.

Coordinates follow the SDK: intervals are 0-based and half-open (like BED); variant positions are 1-based. In a shell, quote the `chr:pos:ref>alt` form (`'chr22:36201698:A>C'`) or use `chr22-36201698-A-C`: an unquoted `>` is a redirect. A prediction window is 16 KB, 100 KB, 500 KB or 1 MB (default), centred on the region or variant.

```bash
igvfagent alphagenome setup --install --ping           # SDK environment + key check
igvfagent alphagenome metadata --search heart          # tracks and ontology CURIEs
igvfagent alphagenome score-variants --variants rs429358 chr22-36201698-A-C --ontology UBERON:0000948
igvfagent alphagenome predict-variant --variant rs429358 --outputs RNA_SEQ DNASE --ontology UBERON:0000948
igvfagent alphagenome predict-interval --gene APOE --outputs RNA_SEQ --ontology UBERON:0002107
igvfagent alphagenome ism --ism-interval chr19:44908670-44908700 --scorer DNASE --ontology UBERON:0000948
igvfagent alphagenome atlas-variants --variants rs429358 --genes APOE
```

Agent tools: `alphagenome_setup`, `alphagenome_metadata`, `alphagenome_predict_interval`, `alphagenome_predict_variant`, `alphagenome_score_variants`, `alphagenome_score_interval`, `alphagenome_ism`, `alphagenome_atlas_scorers`, `alphagenome_atlas_variants`, `alphagenome_atlas_interval`, `alphagenome_selftest`.

---

[← Documentation index](../README.md) · [Project README](../../../README.md)
