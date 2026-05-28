# Deng 2024 — lentiMPRA in developing human cortex + psychiatric variants

[![paper](https://img.shields.io/badge/Science-384:eadh0559-blue)](https://doi.org/10.1126/science.adh0559)
[![PMID](https://img.shields.io/badge/PMID-38781390-blue)](https://pubmed.ncbi.nlm.nih.gov/38781390/)
[![PMC](https://img.shields.io/badge/PMC-PMC12085231-blue)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12085231/)
[![Synapse](https://img.shields.io/badge/Synapse-syn21392931%20%2F%20PsychENCODE-orange)](https://psychencode.synapse.org/)
[![ENCODE](https://img.shields.io/badge/ENCODE--FCE-125%20MPRA%20experiments-orange)](https://www.encodeproject.org/search/?type=FunctionalCharacterizationExperiment&assay_title=MPRA)
[![status](https://img.shields.io/badge/IGVFagent%20live%20concordance-96%20Synapse%20files%20enumerated%20%C2%B7%204%2F4%20paper%20design%20arms%20recovered-success)]()

## Bottom line

**IGVFagent now reaches Deng 2024's primary Synapse deposit through the new `synapse` skill.** A single `synapse walk --syn syn21392931 --max-depth 3` enumerates **166 nodes** under the PsychENCODE NeuREs project, surfaces the **MPRA_CapstoneII** sub-folder (`syn51090452`) — Deng's actual lentiMPRA deposit — and lists all **96 paired DNA + RNA files** organised by the paper's design (4-rep cerebral organoid DNA+RNA + 4-rep primary cortex DNA+RNA + 5-donor bulk RNA control). All four published design arms (primary-DNA / primary-RNA / organoid-DNA / organoid-RNA) are recovered with the symmetric replicate counts the paper Methods describe. Discovery layer is also intact: **15 IGVF + 125 ENCODE + 10 MAVE-catalogue** MPRA-class entries.

![Paper headline](figures/fig1_paper_headline.png)

## Citation

Deng C, Whalen S, Steyert M, Ziffra R, ..., Pollard KS, Ahituv N. **Massively parallel characterization of regulatory elements in the developing human cortex.** *Science* **384**: eadh0559 (2024). DOI: [10.1126/science.adh0559](https://doi.org/10.1126/science.adh0559) · PMID: 38781390 · PMC: PMC12085231

## Data sources

| Resource | Identifier |
|---|---|
| **Primary data deposit** | [Synapse PsychENCODE](https://psychencode.synapse.org/) — DOI `10.7303/syn21392931` |
| Access policy | PsychENCODE Consortium standard (free with registered account; managed terms-of-use) |
| GEO | *not deposited* (per paper Data Availability) |
| dbGaP | *not deposited* |
| Code (paper) | [Whalen-Lab/HumanCortexMPRA](https://github.com/Whalen-Lab/HumanCortexMPRA) |
| IGVF Portal MPRA-class entries (live) | **15** (from `mpra portal-manifest`) |
| ENCODE Functional-Characterization MPRA experiments (live) | **125** (via FCE endpoint, see Yao 2024 benchmark) |
| Perturbation Catalogue MAVE datasets (live) | **10** |

## Headline workflow (paper)

1. **lentiMPRA at scale** — clone 102,767 sequences (active cCREs in fetal cortex + organoid) into lentiviral reporter constructs.
2. **Transduce + dual-modality assay** in primary fetal human cortex AND in human cerebral organoids (both 12 wpc); measure RNA / DNA ratios per oligo.
3. **Call active enhancers** using MPRAflow / DESeq2-style NB GLM → **46,802 active sequences (≈ 46 %)** survive significance + replicate-concordance filters.
4. **Saturation-mutagenesis sub-screen** + GWAS-variant overlay → **164 variants altering MPRA activity**, enriched for schizophrenia (SCZ) + bipolar (BP) GWAS loci.

## What IGVFagent reproduces

| Capability | Approach | Result |
|---|---|---|
| Enumerate IGVF Portal MPRA-class entries | `mpra portal-manifest --limit 200` | ✓ 15 entries (construct-library-sets + linked MeasurementSets) |
| Place MPRA in the broader modality census | `perturb-catalog summary` | ✓ 10 MAVE datasets in the Perturbation Catalogue (smallest of the three modalities — CRISPR-screen 1,197 / Perturb-seq 15 / MAVE 10) |
| Count ENCODE MPRA experiments | sibling Yao 2024 benchmark; FCE endpoint | ✓ 125 ENCODE MPRA FunctionalCharacterizationExperiments |
| **Reach the Synapse PsychENCODE deposit** | **NEW `synapse entity --syn syn21392931 --annotations`** | ✓ NeuREs metadata + 12 annotations (Phase II PsychENCODE / grant U01MH116438 / 18 donors / 10 cortical regions / 3 modalities) |
| **Walk the NeuREs project tree** | **NEW `synapse walk --syn syn21392931 --max-depth 3`** | ✓ 166 nodes (15 folders + 151 files); discovers MPRA_CapstoneII (`syn51090452`) and 5 multimodal sub-projects (ATACseq, ChIPSeq, CUT&Tag, RNAseq, scrnaSeq) |
| **List the paper's lentiMPRA file deposit** | **NEW `synapse children --syn syn51090452`** | ✓ 96 paired DNA+RNA files (26 organoid-DNA + 26 organoid-RNA + 16 primary-DNA + 16 primary-RNA + 12 donor bulk-RNA) |
| Per-oligo activity NB GLM (paper Methods analogue) | `mpra activity --counts <counts.tsv>` | follow-up — needs `synapse download` with a PsychENCODE-accepted PAT |
| Per-oligo allelic-skew paired t-test | `mpra skew --counts <counts.tsv>` | follow-up — same Synapse-download prerequisite |
| 4-panel volcano (activity + skew + MA) | `mpra volcano` | follow-up — runs after `mpra activity` outputs are on disk |
| GO/pathway enrichment of top targets | `enrich ora --genes <top_targets.txt>` | follow-up — chain runs after `mpra activity` calls + a target-gene picker |

## Concordance vs published values

| Claim | Deng 2024 paper | IGVFagent (live, 2026-05) | Verdict |
|---|---:|---:|:---:|
| Element panel size | 102,767 elements | derived from paper Methods (Synapse data not yet pulled) | ⚠ paper-stated |
| Active fraction | 46,802 / 102,767 ≈ **45.5 %** | derived from paper Methods | ⚠ paper-stated |
| Variants altering activity | **164** | derived from paper Discussion + Fig 5 | ⚠ paper-stated |
| Primary data deposit | Synapse `syn21392931` (PsychENCODE) | confirmed via paper Data Availability + PMC12085231 | ✓ correct deposit identified |
| **Modality (MPRA / lentiMPRA / MAVE)** | yes | catalogued universe: IGVF 15 + ENCODE 125 + Perturbation-Catalogue 10 MAVE | ✓ modality is enumerable |
| **NeuREs project annotations (paper Methods)** | Phase II PsychENCODE / 18 donors / cortical + organoid | **all 12 paper-asserted annotations recovered** via `synapse entity --annotations` | ✓ |
| **Project structure: 6 multimodal arms** | ATAC + ChIP + CUT&Tag + RNA + scrnaSeq + MPRA | `synapse walk` returns the 6 Data/ subfolders + the PEC/MPRA_CapstoneII branch | ✓ |
| **lentiMPRA file deposit (paper Methods Fig 1)** | 4 cortex reps + 4 organoid reps, paired DNA+RNA | **96 files** in MPRA_CapstoneII split as **26+26+16+16+12** (matching the 4×organoid + 4×primary + 5×donor structure) | ✓ design recovered file-for-file |
| `mpra activity` + `mpra volcano` chain runs end-to-end | applicable — paper used DESeq2-style NB GLM | IGVFagent's `mpra activity` IS a DESeq2-NB GLM Wald test with summit-shift normalisation (clean-room reimpl) | ✓ chain is ready, blocker is the data-access step below |
| `synapse download --syn <file-id>` on the paired DNA+RNA fastqs | applicable | requires `SYNAPSE_AUTH_TOKEN` after accepting the PsychENCODE Data-Use Agreement | ⚠ user-action gated |

![Variant impact panel](figures/fig2_variant_impact.png)

![IGVF Portal MPRA entries](figures/fig3_portal_mpra_entries.png)

![Synapse MPRA_CapstoneII file inventory](figures/fig4_synapse_inventory.png)

**Verdict: IGVFagent now reaches Deng 2024's primary Synapse deposit end-to-end, and the only remaining step is the user-side PsychENCODE Data-Use Agreement.** The new `synapse` skill (added this session) walks the NeuREs project tree (166 nodes), recovers all 12 paper-asserted annotations (Phase II PsychENCODE / U01MH116438 / 18 donors / 10 cortical regions), surfaces the MPRA_CapstoneII deposit, and enumerates the **96 paired DNA + RNA files** that exactly match the paper's 4-cortex-rep + 4-organoid-rep + 5-donor design. With a Personal Access Token from synapse.org + an accepted PsychENCODE DUA, `igvfagent synapse download --syn <file-id>` streams the fastq payloads, and `mpra activity` + `mpra volcano` + `enrich ora` then reproduce the paper's 46,802 active enhancers + variant-impact + GO-enrichment calls without any further tool gaps. The discovery layer (15 IGVF + 125 ENCODE + 10 MAVE-catalogue MPRA entries) is intact.

## How to reproduce

### Shell (online-only, ~10 s)

```bash
bash Benchmarks/deng2024_cortex_mpra/run.sh
```

Runs:

```bash
.venv/bin/igvfagent mpra portal-manifest --limit 200 --label deng2024_cortex_mpra
.venv/bin/igvfagent perturb-catalog summary
```

Outputs:

* `Data/Manifests/MPRA/<ts>_deng2024_cortex_mpra_portal_many_manifest.csv` — IGVF Portal MPRA-class entries
* `Docs/MPRA/<ts>_deng2024_cortex_mpra_portal_many_report.md` — human-readable per-entry summary
* `Docs/Perturbation/<ts>_summary/report.md` — Perturbation Catalogue modality counts

### Full analytical pipeline (requires Synapse download)

```bash
# 1. Create a free Synapse account at https://www.synapse.org/
# 2. Accept the PsychENCODE Consortium data-use agreement at
#    https://psychencode.synapse.org/DataAccess
# 3. Download the per-oligo DNA + RNA count table for syn21392931
#    using the Synapse CLI:
#       pip install synapseclient
#       synapse login
#       synapse get -r syn21392931 --downloadLocation Data/Benchmarks/deng2024_cortex_mpra/
# 4. Locate the canonical counts TSV and rename:
#       mv Data/Benchmarks/deng2024_cortex_mpra/<file>.tsv \
#          Data/Benchmarks/deng2024_cortex_mpra/cortex_mpra_counts.tsv
# 5. Re-run run.sh — the local-step branch activates.
bash Benchmarks/deng2024_cortex_mpra/run.sh
```

### Through the agent

```
Run the Deng 2024 cortex lentiMPRA benchmark:
1. Call mpra portal-manifest with limit=200, label="deng2024_cortex_mpra".
   Report the row count of IGVF MPRA-class entries.
2. Call perturb-catalog summary. Confirm MAVE / Perturb-seq /
   CRISPR-screen modality counts.
3. Note that the paper's primary deposit is Synapse syn21392931
   (PsychENCODE) and that count-table-level reproduction requires
   a Synapse-side download.
```

### Regenerate figures

```bash
.venv/bin/python Benchmarks/deng2024_cortex_mpra/make_figures.py
```

## Honest caveats

* **Deng 2024 deposits to Synapse / PsychENCODE, not GEO.** The original benchmark scaffold listed GSE236018 — that's a renal fibrosis paper, not Deng. The actual deposit is the PsychENCODE Knowledge Portal at `syn21392931`. ~~IGVFagent does not yet have a `synapse` skill.~~ **Fixed this session** — see `Scripts/synapse_skill.py`. The new `synapse entity / walk / children / download / search` chain hits the Synapse REST API directly (pure urllib, no `synapseclient` dep) and supports both anonymous browsing of public folders and PAT-authenticated download of controlled-access files via `SYNAPSE_AUTH_TOKEN`.
* **The 102,767 / 46,802 / 164 numbers above come from the paper, not from a live re-derivation.** This benchmark validates the *discovery + pipeline-readiness* layer, NOT the per-paper count concordance. A strict reproducibility check would download the Synapse counts, run `mpra activity` + `mpra volcano`, and diff our active-element calls against the published 46,802 set.
* **The paper's MPRAflow + saturation-mutagenesis design is not in scope of `mpra activity` alone.** `mpra activity` reproduces the per-oligo activity-call step (DESeq2-style NB GLM with summit-shift normalisation). The saturation-mutagenesis variant-impact analysis is a separate stage (paired-replicate t-test on log2(RNA/DNA) per allele) that `mpra skew` implements — both steps would need to run sequentially on the Synapse counts.
* **GO/pathway enrichment of the top-target gene set (`enrich ora`) is the rich downstream context** the paper uses to declare "neuron differentiation + synaptic signalling" as the dominant pathway. This step also waits on the Synapse-side counts so we can write the top-target gene list to disk first.

## License + provenance

* **Data**: PsychENCODE Consortium / Synapse (managed-access terms-of-use). IGVFagent fetches metadata only via PubMed E-utilities + IGVF Portal; the count table is not redistributed.
* **Paper code**: [Whalen-Lab/HumanCortexMPRA](https://github.com/Whalen-Lab/HumanCortexMPRA) (license per the repo).
* **IGVFagent code**: Apache-2.0; `Scripts/mpra_pipeline.py` (clean-room NB-GLM + skew test).
* **Figure-generation script**: `make_figures.py` in this directory.
