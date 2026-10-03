# Skill: Proteomics & PPI knowledge graph

End-to-end proteomics skill for IGVFagent. Aggregates **BioGRID**,
**IntAct**, **HuRI**, **Reactome**, **KEGG**, **UniProt id-mapping**,
and **IGVF Portal protein assays** into a local SQLite knowledge graph
that is incrementally update-aware. Provides summary statistics,
network visualizations (degree distribution, hub bar plot, ego graph),
literature surveys for VAMP-seq family / MAVE / semi-qY2H / DUAL-IPA via
the Reference skill, and per-assay example figures from real Portal
files.

## Subcommands

### download
```
igvfagent proteomics download --source all
igvfagent proteomics download --source biogrid
igvfagent proteomics download --source intact
igvfagent proteomics download --source huri
igvfagent proteomics download --source reactome
igvfagent proteomics download --source kegg --kegg-max-pathways 400
igvfagent proteomics download --source igvf
igvfagent proteomics download --source uniprot   # idmap, optional
```
The skill stores `Data/Proteomics/_versions.json` with version, URL, sha256,
and record count per source. BioGRID is auto-resolved against the latest
`Release-Archive` listing. IntAct uses `Last-Modified` on
`psimitab/intact-micluster.txt`. HuRI is mostly static.

### versions / update
```
igvfagent proteomics versions
igvfagent proteomics update
```
`update` only re-fetches sources where the upstream version differs from
the local one (delta refresh).

### igvf-protein
```
igvfagent proteomics igvf-protein
```
Pulls all 214 protein-slim MeasurementSets, the 14 PPI-score files, the
DUAL-IPA stability file, the UniProt protein references, and the protein
language model files from the IGVF Portal. Saves metadata JSON + actual
TSV/CSV file content under `Data/Proteomics/Sources/IGVF/`.

### build-kg
```
igvfagent proteomics build-kg --sources all
igvfagent proteomics build-kg --sources biogrid,reactome,igvf
```
Ingests downloaded sources into `Data/Proteomics/KG/proteomics.sqlite`.
Schema: `interactions`, `proteins`, `pathways`, `pathway_membership`,
`igvf_evidence`, `id_map`, `versions`. Edges deduped on
`(id_a, id_b, source, source_id)`.

### kg-stats
```
igvfagent proteomics kg-stats --label initial
```
Summary report under `Docs/Proteomics/<ts>_initial/`:
total interactions, distinct proteins, pathways, IGVF evidence,
breakdowns by source / evidence type / detection method, top-30 hubs.

### kg-visualize
```
igvfagent proteomics kg-visualize --gene TP53 --label tp53
igvfagent proteomics kg-visualize --label snapshot
```
Produces `Plots/degree_distribution.png`, `Plots/top_hubs.png`,
`Plots/per_source.png`, and (when `--gene` given) an ego graph
`Plots/ego_<GENE>.png` (requires networkx).

### assay-survey
```
igvfagent proteomics assay-survey --label may2026
```
Calls the Reference skill (PubMed + Semantic Scholar + OpenAlex) to
retrieve current studies on:
- VAMP-seq (MultiSTEP)
- VAMP-seq
- MAVE
- semi-qY2H v1 / v2 / v3 (semi-quantitative yeast two-hybrid)
- DUAL-IPA
Filters to Nature / Cell / Science family journals (Nat Genet, Nat Methods,
Nat Biotechnol, Mol Cell, Mol Syst Biol, Nat Commun) and writes
`literature_survey.md` and `.json`.

### assay-figures
```
igvfagent proteomics assay-figures --label demos
```
Generates per-assay example histograms from the actual IGVF Portal files
pulled by `igvf-protein`. One PNG per assay under
`Docs/Proteomics/<ts>_demos_assay_figures/Plots/`.

### vampseq-pull / vampseq-analyze / vampseq-inventory
```
# 1) Pull canonical published scoresets from MaveDB (PTEN, TPMT, VKOR,
#    PRKN, CYP2C9, NUDT15) — full per-replicate score CSVs
igvfagent proteomics vampseq-pull
igvfagent proteomics vampseq-pull --gene PTEN

# 2) Deep analysis — produces the Fowler-lab-style plot suite per gene
igvfagent proteomics vampseq-analyze --gene PTEN --label pten_deep     --pdb-id 1d5r                          # optional PyMOL .pml export
igvfagent proteomics vampseq-analyze       # all 6 catalogued targets

# 3) Inventory the IGVF Portal raw VAMP-seq experiments by decoding the
#    alias scheme (<lab>:<GENE>-DMS-<antibody>-Tile<i>-Replicate<j>-Bin<k>).
#    Produces tile×bin and antibody×tile coverage matrices per gene.
igvfagent proteomics vampseq-inventory --label igvf_f9
```
The deep analysis follows the canonical VAMP-seq pipeline distilled from
Matreyek et al. *Nat Genet* 2018, Suiter *eLife* 2020, Clausen *Nat Commun*
2024, and Coyote-Maestas *Nat Commun* 2024 (MultiSTEP), augmented to
match the public **FowlerLab/VAMPseq** analysis Rmd plot suite:

  1. **Distribution** — overlay missense / synonymous / nonsense densities
     anchored at WT=1, nonsense=0.
  2. **Residue × AA heatmap** — the iconic VAMP-seq view; one cell per
     (position, substituted AA), color = abundance score.
  3. **Per-position mean ± IQR with domain track** — annotates which
     domain each residue belongs to (e.g. PTEN phosphatase / C2 / C-tail,
     PRKN Ubl / RING0 / RING1 / IBR / RING2).
  4. **Per-position median + 3-residue moving average** — Fowler-lab
     `ma(score, n=3)` style smoothing alongside the per-residue median.
  5. **Replicate concordance** — Pearson r between rep-1 and rep-2 per
     variant.
  6. **N × N replicate scatter matrix** — emitted when ≥ 3 replicates
     are present (e.g. PTEN with 8 reps). Lower triangle = scatter,
     diagonal = per-rep histogram, upper triangle = Pearson r.
  7. **Abundance class breakdown** — categorical bar (low / low-int /
     intermediate / WT-like / hyper-abundant), matching the Matreyek
     2018 `abundance_class` convention.
  8. **Cumulative ranked variants with 95 % CI band** — sorted score
     curve, CI envelope from `lower_ci` / `upper_ci` (or normal-approx
     from `se`), low-abundance fraction (score < 0.5) shaded.
  9. **Nonsense-by-position scatter** — the canonical "do truncations
     crash the score?" QC plot from Matreyek 2018 Fig 1.
  10. **Biophysical-feature Spearman ρ panel** — emitted only when the
      scoreset carries RSA / B-factor / hydrophobicity / Grantham /
      PSIC conservation / ΔΔG / Tm columns (most public MaveDB
      scoresets don't; the FowlerLab supplementary PTEN / TPMT tables
      do).
  11. **PyMOL .pml export** — opt-in via `--pdb-id`. Writes a ready-to-
      source `.pml` that loads the structure, sets per-residue
      B-factor to the median missense abundance score, and applies the
      blue → white → red colorscale. Default PDB IDs are pinned to
      the canonical Fowler-lab structures (PTEN → `1d5r`,
      TPMT → `2bzg`).

Catalogued MaveDB targets (URN, paper, length, domains) are in
`MAVEDB_VAMPSEQ_CATALOG` in `proteomics_skill.py` — extend this dict to
analyze additional published scoresets.

The `vampseq-inventory` command decodes the IGVF Portal `aliases` field
to build a coverage matrix across the 144 MultiSTEP MeasurementSets
(currently all targeting **F9 / Coagulation Factor IX** across 3 tiles ×
4 bins × 4 replicates × 5 antibody readouts) and the 36 plain VAMP-seq
sets (CYP2C19, G6PD).

### pipeline
```
igvfagent proteomics pipeline --label may2026 --gene TP53 \
  --sources biogrid,intact,huri,reactome,kegg,igvf
```
End-to-end: download → build-kg → kg-stats → kg-visualize → assay-figures
→ assay-survey. Use `--skip-download` to reuse cached files,
`--skip-literature` to skip the network-dependent literature scan,
`--max-rows N` to cap parser ingestion (useful for smoke tests).

### write-playbook
```
igvfagent proteomics write-playbook
```
Writes this document to `Docs/Skills/PROTEOMICS_SKILLS.md`.

## Storage

```
Data/Proteomics/
  _versions.json
  Sources/
    BioGRID/BIOGRID-ALL-X.Y.Z.tab3.zip
    IntAct/intact-micluster.txt
    HuRI/{HuRI.tsv, HI-union.tsv, Lit-BM.tsv}
    Reactome/{ReactomePathways.txt, UniProt2Reactome.txt,
              reactome.homo_sapiens.interactions.tab-delimited.txt}
    KEGG/{pathways_hsa.json, pathway_membership_hsa.json}
    UniProt/HUMAN_9606_idmapping_selected.tab.gz
    IGVF/{measurement_sets.json, ppi_files_meta.json,
          stability_files_meta.json, protein_reference_meta.json,
          protein_language_models.json, ppi_files/IGVFFI*}
  KG/proteomics.sqlite

Docs/Proteomics/<ts>_<label>/
  report.md
  Plots/*.png
  literature_survey.md (if assay-survey)
  per_assay_figures/Plots/*.png (if assay-figures)
```

## Cross-skill chaining

- `proteomics build-kg` → `kg gene <SYM>` to combine the proteomics PPI
  graph with the IGVF Catalog Knowledge Graph (variants, regulatory
  elements, pathways).
- `proteomics assay-survey` → `ref design --query "<assay>"` to draft
  a study design that mirrors the canonical analysis flow used by
  recent Nature/Cell/Science studies.
- `rnaseq deg` significant DEGs → query
  `proteomics kg-visualize --gene <SYM>` to see what proteins the
  upregulated genes interact with.
- `se-targets pipeline` → for each SE-target gene, look up upstream
  PPI partners with `proteomics kg-visualize` to suggest functional
  cofactors.

## Notes

- KEGG REST is throttled (default 0.4s/req) to respect their TOS;
  bulk distribution requires a license.
- IntAct full PSI-MITAB (`intact.zip`) is ~700MB — the skill defaults
  to the smaller `intact-micluster.txt`. Set `INTACT_BASE` env to
  override.
- HuRI uses Ensembl gene IDs. Run `proteomics download --source uniprot`
  to populate `id_map` for cross-resource harmonization.
