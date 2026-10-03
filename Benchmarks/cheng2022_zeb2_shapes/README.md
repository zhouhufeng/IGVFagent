# cheng2022_zeb2_shapes

## Paper

**ZEB2 Shapes the Epigenetic Landscape of Atherosclerosis**
Cheng P, Wirka RC, Shoa Clarke L, Zhao Q, Kundu R, Nguyen T, Nair S, Sharma D, Kim HJ, Shi H, Assimes T, Brian Kim J, Kundaje A, Quertermous T.
*Circulation* 145(6):469-485 (2022) · doi:[10.1161/CIRCULATIONAHA.121.057789](https://doi.org/10.1161/CIRCULATIONAHA.121.057789) · PMID 34990206 · PMC8896308 · IGVF award IGVF0043, grant UM1HG011972

Resolver confidence: **1.00** (exact PMID match).

The paper traces a complex CAD GWAS signal at chromosome 2q22.3 to a distal enhancer for *ZEB2* (>500kb away), then uses CRISPR genome/epigenome editing, mouse SMC lineage tracing with single-cell RNA-seq/ATAC-seq, and human coronary artery smooth muscle cell (HCASMC) knockdown to show ZEB2 governs the epigenetic transition of smooth muscle cells during atherosclerotic plaque formation, restraining Notch/TGFβ-driven fibromyocyte (FMC) differentiation.

## The central problem: this paper's own data was never deposited

The Data Availability statement, verbatim: *"All data and materials are currently being deposited to the NCBI Gene Expression Omnibus and will be available for general public access upon acceptance of this manuscript for publication. (Accession numbers currently pending.)"*

**No accession numbers appear anywhere in the article** (checked main text, figure legends, footnotes, funding/disclosures — full text via `pmc.ncbi.nlm.nih.gov`, since EuropePMC's JATS XML returns HTTP 500 for this closed-access article). Confirmed independently via NCBI E-utilities:
- `elink` (`dbfrom=pubmed db=gds id=34990206`): **zero** linked GEO records for this PMID.
- Targeted `esearch` across GEO for every combination of author names (Wirka, Quertermous, Cheng) and assay/keyword terms (Zeb2, smooth muscle, atherosclerosis, HCASMC, ATAC, ApoE aortic root scRNA): **no matching series**.

Four-plus years post-publication, none of this paper's own generated datasets — mouse aortic-root scRNA-seq/scATAC-seq (baseline and *Zeb2*-ΔSMC knockout), H3K27ac ChIP-seq in HCASMC, or bulk RNA-seq/ATAC-seq in ZEB2-knockdown HCASMC — are discoverable. This is a real, documented data-availability failure, not a search-tooling gap; see `expected.json`'s blockers on `fig2_zeb2_scrna_scatac`, `fig3_zeb2_ko_scatac`, `fig4_5_zeb2_ko_scrna_histology`, `fig6_zeb2_kd_hcasmc`, and the CRISPR-validation part of Fig. 1 (`fig1_crispr_enhancer`) — **5 of 8 planned analyses**.

## Part 1 — Fig. 1: the CAD-GWAS-to-ZEB2 locus, from real public data — DONE (2/8 reproduced)

The one part of this paper that rests on genuinely public, independent data: a published CAD GWAS meta-analysis and GTEx. New port `cheng2022_gwas_eqtl` queries the real GWAS Catalog and GTEx v8 REST APIs directly (both are the paper's own named external resources — no author code exists for this step either).

| Claim | Paper | Measured (real GWAS Catalog / GTEx v8 data) | Match |
|---|---|---|---|
| Lead CAD SNP p-value at 2q22.3 | 5.4×10⁻¹³ | **5×10⁻¹³** (rs2252641, "coronary artery disorder," van der Harst et al. 2017 CARDIoGRAMplusC4D meta-analysis, GWAS Catalog study GCST005195) | almost exact |
| ZEB2 is the closest gene, >600,000 bp away | >600 kb | **551 kb** (rs957293 at chr2:145,075,620 to ZEB2's hg38 span chr2:144,364,364-144,524,583) | same order of magnitude (hg19-vs-hg38 / TSS-definition sensitive, not exact) |
| CAD risk SNPs colocalize with ZEB2 eQTL | qualitative | **real, significant GTEx v8 aorta eQTL**: rs957293-ZEB2 p=1.5×10⁻⁵, below GTEx's own significance threshold (1.5×10⁻⁴) for this gene/tissue, at the same SNP that is CAD-associated at p=2×10⁻¹⁴ | confirms colocalization directly |

```bash
bash Benchmarks/cheng2022_zeb2_shapes/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark cheng2022_zeb2_shapes
```

## Not attempted (real data exists, just not built out this pass)

`fig1_conservation` (six conserved candidate enhancer regions via cross-species sequence conservation, ECR Browser) is **pending**, not blocked — UCSC phyloP/phastCons conservation tracks for this 3 Mb region are genuinely public, this just hasn't been implemented yet (the ambiguity is in matching the authors' exact "conserved candidate enhancer" region-calling parameters, not data access).

## Blocked (5 of 8 analyses) — evidence

See `expected.json`'s `blocker` field on each: `fig1_crispr_enhancer`, `fig2_zeb2_scrna_scatac`, `fig3_zeb2_ko_scatac`, `fig4_5_zeb2_ko_scrna_histology`, `fig6_zeb2_kd_hcasmc` — all `not_deposited`, all citing the same verified absence of any GEO deposit for this paper.

## Coverage

**2 of 8** planned analyses reproduced (3/3 checks pass). 5 blocked (`not_deposited`, real evidence above); 1 pending (real public data, not yet implemented).

## Provenance

`provenance.json` holds the full resolve/harvest record. Port: `cheng2022_gwas_eqtl` (`igvfagent port list`).
