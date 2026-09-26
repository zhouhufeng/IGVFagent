# IGVF Papers Benchmark: Publications with Raw Data on the IGVF Portal

**Source:** `IGVF Manuscript Tracking List-Sept24.xlsx` (the IGVF consortium’s own manuscript tracker, 7 sheets, 288 tracked manuscript entries as of the Sept 24 snapshot).

**Scope of this refined list:** a manuscript qualifies only if it is genuinely public (published in a journal, or posted as a preprint) **and** its own award group has self-reported that raw data was submitted to the IGVF Portal. Both conditions are checked against the tracker’s own fields — nothing here is inferred from outside sources.

## Method, and an important distinction the tracker itself makes

The tracker carries two *different* signals that are easy to conflate:

1. **`Data Added on Portal`** (self-reported, free text, by the awardee) — the only field that actually speaks to *raw experimental data* being deposited under that publication.
2. **A cross-reference against the `Paper on portal` sheet** (PMID/title matched to a live IGVF Portal `/publications/` record with a UUID) — this only confirms a **citation/reference record** exists on the Portal, not that a dataset is linked to it.

The two disagree far more often than they agree. Of the 32 published-or-preprint papers with a confirmed Portal citation record, **30 self-report their data field as blank or "Not Provided."** A citation existing on the Portal is *not* evidence that raw data was submitted — so signal (2) alone cannot be used to answer "does this paper have raw data on the Portal," and this list uses signal (1) as the actual criterion, with signal (2) shown only as a supporting cross-check where it happens to agree.

**Filter applied:**
- `Status` normalized to `Published` or `Pre-Print Only` (typo/case variants folded in) — this is what "published/arxived" means here. `Accepted` (17 papers, awaiting formal publication) and everything earlier in the pipeline (`Submitted`, `In Revision`, `In Preparation`, 122 papers combined) are excluded by this wording; they are not lost, just out of scope for this specific list.
- `Data Added on Portal` classified as an affirmative signal: literal `Yes`/`Y`, or free text unambiguously describing data present or actively landing (e.g. *"planning for inclusion with catalog release 1"*, *"the bulk of these data are on the production server"*). Free text is preserved verbatim in the table below rather than collapsed to a bare "yes" — the difference between "Yes" and "in progress" matters if you plan to actually pull the data.

**Result:** of 288 tracked manuscripts, 137 are published or preprint-only, and **9** of those self-report raw data on the IGVF Portal.

| Confidence | Meaning |
|---|---|
| ⭐⭐⭐ confirmed | Self-report is affirmative **and** independently corroborated by a matching Portal citation record |
| ⭐⭐ self-reported | Self-report is a plain "Yes," no independent citation match (often because the citation simply isn’t in the `Paper on portal` sheet yet) |
| ⭐ in progress | Self-report describes data as partially landed / still being prepared, not a finished "Yes" |

## Refined list — published/preprint papers with raw data on the IGVF Portal

| # | Title | Journal | Status | Award(s) | Lead Author(s) | PMID (PMCID) | Preprint | Publication link | Data on IGVF Portal (verbatim) | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Mapping the convergence of genes for coronary artery disease onto endothelial cell programs | Nature | Published | UM1HG011972 | Gavin Schnitzler, Helen Kang | 38326615 (PMC10921916) | Mapping the convergence of genes for coronary artery disease onto endothelial cell programs \| bioRxiv ⚠️*(not a URL in the source sheet)* (2022-11-01) | [https://www.nature.com/articles/s41586-024-07022-x](https://www.nature.com/articles/s41586-024-07022-x) | CRISPRi-Perturb-seq and CRISPRi-FlowFISH in teloHAEC | ⭐⭐⭐ confirmed |
| 2 | CRISPR screening uncovers a long-range enhancer for ONECUT1 in pancreatic differentiation and links a diabetes risk variant | Cell Reports | Published | U01HG012051 | Samual Kaplan | 39163202 (PMC11406439) | [https://doi.org/10.1101/2024.04.26.591412](https://doi.org/10.1101/2024.04.26.591412) (2024-04-29) | [https://www.sciencedirect.com/science/article/pii/S2211124724009902?via%3Dihub](https://www.sciencedirect.com/science/article/pii/S2211124724009902?via%3Dihub) | Yes ⚠️ *not-on-portal audit list also matches this PMID/title — worth a manual check* | ⭐⭐ self-reported |
| 3 | Massively parallel characterization of transcriptional regulatory elements in three diverse human cell types | Nature | Published | 1UM1HG011966 | Vikram Agarwal, Fumitaka Inoue | 39814889 (PMC11903340) | [https://www.biorxiv.org/content/10.1101/2023.03.05.531189v1](https://www.biorxiv.org/content/10.1101/2023.03.05.531189v1) (2023-03-05) | [https://www.nature.com/articles/s41586-024-08430-9](https://www.nature.com/articles/s41586-024-08430-9) | Data in ENCODE ⚠️ *not-on-portal audit list also matches this PMID/title — worth a manual check* | ⭐⭐ self-reported |
| 4 | An encyclopedia of enhancer-gene regulatory interactions in the human genome | Nature | Published | UM1HG011972, U01HG012069 | Andreas Gschwind, Kristy Mualim, Alireza Karbalayghareh, Maya Sheth, Kushal Dey, Evelyn Jagoda, Ramil Nurtdinov, Wang Xi | 42457959 (Not Available Yet) | An encyclopedia of enhancer-gene regulatory interactions in the human genome \| bioRxiv ⚠️*(not a URL in the source sheet)* (2023-09-11) | [https://www.nature.com/articles/s41586-026-10781-4](https://www.nature.com/articles/s41586-026-10781-4) | Yes | ⭐⭐ self-reported |
| 5 | Massively parallel regulatory assays for CYP3A4 enhancer variants alongside their native promoter | Nature Communications | Pre-Print Only | 1UM1HG011966 | Yelena Guttman | Not Available Yet | [https://www.biorxiv.org/content/10.64898/2026.04.22.719677v1](https://www.biorxiv.org/content/10.64898/2026.04.22.719677v1) (2026-04-24) | Not Available Yet ⚠️*(not a URL in the source sheet)* | Yes | ⭐⭐ self-reported |
| 6 | eSIG-Net: an interaction language model that decodes the protein code of single mutations | Nature Methods | Published |  | Sidharth Raghavan, Xingxin Pan | 42056223 (PMC13259923) | [https://www.biorxiv.org/content/10.64898/2026.03.27.714913v1](https://www.biorxiv.org/content/10.64898/2026.03.27.714913v1) | [https://www.nature.com/articles/s41592-026-03086-x](https://www.nature.com/articles/s41592-026-03086-x) | Yes ⚠️ *not-on-portal audit list also matches this PMID/title — worth a manual check* | ⭐⭐ self-reported |
| 7 | Predicting the effects of SNPs on transcription factor binding affinity | Bioinformatics | Published | U01HG011952 | Sierra Nishizaki |  |  | [https://doi.org/10.1093/bioinformatics/btz612](https://doi.org/10.1093/bioinformatics/btz612) | yes; planning for inclusion with catalog release 1 | ⭐ in progress |
| 8 | Pervasive mislocalization of pathogenic coding variants underlying human disorders | Cell | Published | UM1HG011989 | Jessica Lacoste | 39353438 (PMC11669423) | [https://www.biorxiv.org/content/10.1101/2023.09.05.556368v1](https://www.biorxiv.org/content/10.1101/2023.09.05.556368v1) (2023-09-05) | [https://www.cell.com/cell/fulltext/S0092-8674(24)01021-3](https://www.cell.com/cell/fulltext/S0092-8674(24)01021-3) | In progress ⚠️ *not-on-portal audit list also matches this PMID/title — worth a manual check* | ⭐ in progress |
| 9 | Multiplex profiling of developmental cis-regulatory elements with quantitative, single-cell expression reporters | Nature Methods | Published | UM1HG011966 | J-B Lalanne, S. G. Regalado | 38724692 (PMC11166576) | Multiplex profiling of developmental enhancers with quantitative, single-cell expression reporters \| bioRxiv ⚠️*(not a URL in the source sheet)* (2022-12-10) | [https://www.nature.com/articles/s41592-024-02260-3](https://www.nature.com/articles/s41592-024-02260-3) | The bulk of these data are on the production server; seqspec files being worked on | ⭐ in progress |

**Notes on specific rows:**
- **IGVF0004** (*Data in ENCODE*): the data lives on the **ENCODE** portal, not literally the IGVF Portal — same underlying platform family, but a different deposit. Included here because that's the source's own answer to the question; flagged so it isn't mistaken for an IGVF Portal deposit.
- **IGVF0001, IGVF0017, IGVF0003**: the `PrePrint Server URL` field in the source sheet contains the paper's *title* rather than an actual link for these three rows — a data-entry artifact in the tracker itself, not something introduced here.
- **IGVF0016, IGVF00100, IGVF0004, IGVF0203**: flagged ⚠️ above — their PMID or title also matches an entry in the tracker's independent "Accepted IGVF papers not on portal" audit list. This is most likely because that audit list was matched against the **preprint's** PMID rather than the final published one (three of these four have both), but it's a genuine discrepancy in the source data worth a manual look before relying on it.

## Appendix — published/preprint papers with a confirmed Portal citation, but no data self-report

These 30 papers are independently confirmed to have a citation record on the IGVF Portal (PMID or title match against the `Paper on portal` sheet), but the tracker’s own `Data Added on Portal` field is blank or "Not Provided." A Portal citation commonly exists before, or independent of, any dataset being linked to it (see Method above) — so these are **candidates worth checking directly against the live Portal**, not confirmed raw-data submissions. Not merged into the refined list above for that reason.

| Title | Journal | PMID | Tracking # | Award(s) |
|---|---|---|---|---|
| Scalable functional assays for the interpretation of human genetic variation | Annual Review of Genetics | 36055970 | IGVF0012 | UM1HG011989 |
| Computational identification of clonal cells in single-cell CRISPR screens. | BCM Genomics | 35168568 | IGVF0059 | UM1HG011996 |
| Computational and experimental methods for classifying variants of unknown clinical significance | CSH: Molecular Case Studies | 35483875 | IGVF0028 | UM1HG011966 |
| Analysis of estrogen-regulated enhancer RNAs identifies a functional motif required for enhancer assembly and gene expression. | Cell | 35705040 | IGVF0056 | UM1HG011996 |
| Lymph node colonization induces tumor-immune tolerance to promote distant metastasis. | Cell | 35525247 | IGVF0074 | UM1HG012076 |
| Massively parallel reporter assay investigates shared genetic variants of eight psychiatric disorders | Cell | 39848247 | IGVF0136 | 1UM1HG012003 |
| Activation of the imprinted Prader-Willi syndrome locus by CRISPR-based epigenome editing | Cell Genomics | 39947136 | IGVF0192 | UM1HG012053 |
| High-throughput characterization of the role of non-B DNA motifs on promoter function. | Cell Genomics | 35573091 | IGVF0027 | UM1HG011966 |
| Accurate Classification of Cardiomyopathy Diagnosis by Chromatin Accessibility | Circulation | 36095061 | IGVF0054 | UM1HG011996 |
| ZEB2 Shapes the Epigenetic Landscape of Atherosclerosis | Circulation | 34990206 | IGVF0043 | UM1HG011972 |
| High-throughput techniques enable advances in the roles of DNA and RNA secondary structures in transcriptional and post-transcriptional gene regulation. | Genome Biology | 35851062 | IGVF0026 | UM1HG011966 |
| A gene regulatory element modulates myosin expression and controls cardiomyocyte response to stress | Genome Research | 41125440 | IGVF0195 | UM1HG012053 |
| Thymic and extrathymic Aire-expressing cells in maternal-fetal tolerance. | Immunological Reviews | 35535447 | IGVF0073 | UM1HG012076 |
| Characterization of De Novo Promoter Variants in Autism Spectrum Disorder with Massively Parallel Reporter Assays | Int J Mol Sci | 36834916 | IGVF0024 | UM1HG011966 |
| Focus on your locus with a massively parallel reporter assay. | J Neurodev Disord. | 36085003 | IGVF0060 | UM1HG012003 |
| Convergence and Divergence of Rare Genetic Disorders on Brain Phenotypes: A Review. | JAMA Psychiatry | 35767289 | IGVF0061 | UM1HG012003 |
| A RORγt+ cell instructs gut microbiota-specific Treg cell differentiation. | Nature | 36071167 | IGVF0069 | UM1HG012076 |
| Transition to a mesenchymal state in neuroblastoma confers resistance to anti-GD2 antibody via reduced expression of ST8SIA1. | Nature Cancer | 35817829 | IGVF0071 | UM1HG012076 |
| Massively parallel reporter perturbation assays uncover temporal regulatory architecture during neural differentiation | Nature Communications | 35315433 | IGVF0029 | UM1HG011966 |
| Systematic reconstruction of cellular trajectories across mouse embryogenesis. | Nature Genetics | 35288709 | IGVF0030 | UM1HG011966 |
| Transcriptional and epigenetic regulators of human CD8+ T cell function identified through orthogonal CRISPR screens. | Nature Genetics | 37945901 | IGVF0062 | UM1HG012053 |
| Valid inference for machine learning-assisted genome-wide association studies | Nature Genetics |  | IGVF0097 | U01HG012039 |
| Epigenetic regulation of T cell exhaustion. | Nature Immunology | 35624210 | IGVF0072 | UM1HG012076 |
| Runx3 drives a CD8+ T cell tissue residency program that is absent in CD4+ T cells. | Nature Immunology | 35882933 | IGVF0070 | UM1HG012076 |
| NEAT-seq: simultaneous profiling of intra-nuclear proteins, chromatin accessibility and gene expression in single cells | Nature Methods | 35501385 | IGVF0042 | UM1HG011972 |
| Large-scale discovery of neural enhancers for cis-regulation therapies | Not Provided |  | IGVF0223 | 1UM1HG011966 |
| Multi-scale dissection, compaction and derivatization of mammalian developmental enhancers | Not Provided |  | IGVF0225 | 1UM1HG011966 |
| A single factor elicits multilineage reprogramming of astrocytes in the adult mouse striatum. | Proc Natl Acad Sci U S A. | 35254903 | IGVF0058 | UM1HG011996 |
| Mechanosensitive genomic enhancers potentiate the cellular response to matrix stiffness | Science | 40997217 | IGVF0196 | UM1HG012053 |
| Long-read sequencing transcriptome quantification with lr-kallisto | plos computational biology | 41325434 | IGVF0140 | HG012077 |

## Everything excluded, and why

| Reason excluded | Count |
|---|---|
| Not yet published or preprinted (`In Preparation`, `Submitted`, `In Revision`, `Accepted`, `Requested Information`, blank/duplicate) | 151 |
| Published/preprint, but data field is a clear negative (`No`, `Not Provided`, `N/A`) or genuinely blank with no independent Portal citation either | 98 |
| Published/preprint, blank data field, but *has* an independent Portal citation — see Appendix | 30 |

---
*Generated from `Benchmarks/IGVFpapers/IGVF Manuscript Tracking List-Sept24.xlsx` by IGVFagent. Re-run the classification if the tracker is updated — the source spreadsheet, not this file, is the record of truth.*
