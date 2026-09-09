# Evaluation of the Hosted IGVFagent Using Public MultiFlow-Relevant Data

**System:** Hosted IGVFagent, with the interface displaying version 0.2.9  
**Scope:** End-to-end testing through the hosted web interface, using public metadata only

## Executive summary

We evaluated the hosted IGVFagent with four focused tasks designed to resemble data-discovery and reasoning steps relevant to MultiFlow, while maintaining a strict boundary around unpublished project data. The tests covered public multi-omic metadata retrieval, RNA/ATAC sample pairing, uploaded-table validation, multi-turn experimental-design reasoning, and kidney-specific regulatory-evidence retrieval.

The first three tests produced useful and largely accurate results. IGVFagent correctly reconstructed the sample structure of the public GEO series GSE213151, audited a sanitized public manifest, retained the attachment across conversational turns, and generated the requested holdout designs without separating paired RNA and ATAC measurements.

The fourth test revealed substantive scientific-correctness problems. The gene-centric workflow used a region-level endpoint whose rows can target genes other than the queried gene. The final answer did not consistently enforce exact target-gene matching. It also treated a truncated set of 40 records as sufficient to make negative claims. This produced false-negative conclusions for WT1 and HNF4A and misassigned non-target-gene edges to GATA3 and SOX9. Prediction scores were also discussed using significance language despite null adjusted P values and null significance fields.

The most important recommendations are therefore to enforce exact target-gene filtering, prevent absence claims on truncated retrievals, add record-level answer verification, separate predictive scores from inferential significance, and isolate generated artefacts by run. These should be addressed before the gene-regulatory workflow is used for scientific prioritization.

## 1. Data-safety boundary

No unpublished MultiFlow result was uploaded to the hosted service. Specifically, the evaluation did not upload or expose:

- MultiFlow predictions or reconstructed profiles;
- unpublished benchmark metrics or model rankings;
- candidate regulatory relationships or prioritized genes;
- donor-level, replicate-level, or other sensitive metadata;
- project H5AD, RDS, count-matrix, or result-table files.

The only uploaded file was `GSE213151_public_manifest.csv`, a 2,043-byte table derived entirely from the public GEO series GSE213151. It contains ten public sample labels, twenty public GSM accessions, and public GEO filenames. Its SHA-256 checksum was:

`d6ce821ac90c09456b665c3d13512579ccb4004c93fa121c5199ee7c4dd123a0`

The hosted interface identifies its workspace as shared. Unpublished MultiFlow files should therefore continue to be tested only in a local or otherwise isolated private deployment.

## 2. Terminology and evaluation criteria

| Canonical term | Meaning in this report |
|---|---|
| Biological sample | One cell-line-by-day condition represented by a paired RNA and ATAC measurement in the public manifest |
| Pairing | Assignment of exactly one RNA GSM and one ATAC GSM to the same biological sample |
| Target gene | The gene referenced in the `gene` field of a specific regulatory-edge record, not merely a gene whose genomic locus was used as the region query |
| Kidney-relevant evidence | A record whose biological context explicitly identifies kidney or renal cortex; adrenal records were not treated as kidney evidence |
| Observed evidence | A record classified by the source as observed data, such as a perturbation result |
| Prediction | A model-derived linkage record, such as ENCODE-rE2G, which must not be described as statistically significant unless the source provides an explicit significance field and threshold |
| Absence claim | A statement that no qualifying evidence exists, which requires an exhaustive query rather than a truncated result set |

Each round was assessed for factual correctness, adherence to the prompt, appropriate treatment of uncertainty, data provenance, efficiency, and potential scientific impact of any error.

## 3. Overview of the focused tests

| Round | Capability tested | Execution | Verdict |
|---|---|---:|---|
| 1 | Live GEO discovery and RNA/ATAC pairing | 7 iterations, 11 tool calls, approximately 95 s | Pass |
| 2 | Audit of an uploaded public manifest | 2 iterations, 1 tool call, approximately 45 s | Mostly pass |
| 3 | Multi-turn memory and leakage-aware holdout design | 2 iterations, 1 tool call, approximately 50 s | Mostly pass |
| 4 | Kidney-specific regulatory-evidence retrieval for six genes | 21 iterations, 101 tool calls, approximately 405 s | Fail |

The reported times were measured approximately from prompt submission to completion in the hosted interface. They should be interpreted as user-observed end-to-end latency rather than controlled backend benchmarks.

## 4. Round 1: GEO metadata discovery and modality pairing

### 4.1 Prompt

> Use live GEO metadata for GSE213151. Metadata only; do not download data files. Build a manifest paired by biological sample with these columns: sample ID, cell line, differentiation day, RNA GSM, ATAC GSM, RNA file type, and ATAC file type. Report the total number of biological samples, whether every sample has both modalities, and whether GEO directly supplies (a) raw_feature_bc_matrix.h5 and (b) filtered_peak_bc_matrix.h5. Do not infer file availability: support every conclusion with exact GEO filenames/accessions and mark uncertainty explicitly.

### 4.2 Expected result

- Twenty GSM records should form ten paired biological samples.
- AN1-1 should be represented at d7, d12, d19, and d26.
- H9 should be represented at d7, d12, d16, d19, and d26.
- BJFF6 should be represented only at d26.
- Every biological sample should contain one RNA GSM and one ATAC GSM.
- RNA files should be named `*_filtered_feature_bc_matrix.h5`.
- ATAC files should be named `*_atac_fragments.tsv.gz`, with tabix index files.
- GEO should not be described as directly supplying either `raw_feature_bc_matrix.h5` or `filtered_peak_bc_matrix.h5`.
- The approximately 19-GB series archive should not be downloaded.

### 4.3 Observed result

All substantive checks passed. IGVFagent reconstructed all ten biological samples, assigned the correct cell lines and differentiation days, and paired every RNA GSM with the corresponding ATAC GSM. It correctly identified the available file types and explicitly stated that raw RNA matrices and filtered ATAC peak matrices were not directly listed by GEO.

The agent did not download the large series archive. It retrieved only the approximately 2.6-KB GEO file listing required for metadata inspection.

### 4.4 Assessment

**Verdict: Pass.** The answer was source-grounded, respected the metadata-only constraint, and did not invent unavailable files.

One interface issue remained: generated artefact paths were duplicated in the Artefacts panel.

## 5. Round 2: Audit of the uploaded public manifest

### 5.1 Prompt

> Audit the attached GSE213151 manifest without using outside knowledge. Verify row count, uniqueness of sample IDs and GSMs, one-to-one RNA/ATAC pairing, cell-line × day coverage, and missing-file flags. Return a compact pass/fail table. Explicitly identify which cell lines lack d16 and which time points cannot support cross-cell-line comparisons. Do not download anything and do not alter the file.

### 5.2 Expected result

- Ten rows and ten unique sample IDs.
- Ten unique RNA GSMs and ten unique ATAC GSMs.
- Complete one-to-one RNA/ATAC pairing.
- Only H9 has a d16 sample.
- d7, d12, d19, and d26 support an AN1-1 versus H9 comparison.
- BJFF6 appears only at d26 and cannot support a within-line longitudinal trajectory.
- `raw_rna_matrix_listed=false` and `atac_peak_matrix_listed=false` for every row.

### 5.3 Observed result

IGVFagent correctly reported the row count, all uniqueness checks, complete modality pairing, and the cell-line-by-day coverage. It also correctly identified the limitations of BJFF6 and the absence of raw RNA and filtered ATAC peak matrices in the manifest.

The uploaded file was read without modification, and the system did not download external data.

### 5.4 Issue: validation failure versus known data limitation

The answer displayed the two expected `false` file-availability fields as failed checks. This presentation conflates two different concepts:

1. a malformed or internally inconsistent manifest; and
2. a valid manifest that truthfully records an unavailable input.

The latter is a study-design or data-availability limitation, not a validation failure. A more appropriate status would be “pass with limitation” or “expected missing input.”

### 5.5 Assessment

**Verdict: Mostly pass.** The underlying reasoning was correct, but the pass/fail semantics could cause users to misinterpret a valid manifest as erroneous.

## 6. Round 3: Multi-turn memory and holdout design

### 6.1 Prompt

> Using only the attached manifest, design two evaluation folds: (1) leave H9 d16 out for temporal interpolation and (2) hold AN1-1 out at d19 for cross-cell-line generalization. List the exact training and test sample IDs for each fold. Explain which comparisons are valid and flag any leakage if the held-out sample or its paired modality is included in training. Do not invent biological replicates.

### 6.2 Observed result

IGVFagent retained the attachment across turns and did not request that the file be uploaded again.

For the temporal-interpolation fold, it correctly assigned `H9_d16` to the test set and excluded both associated modalities from training:

- RNA: `GSM6573634`
- ATAC: `GSM6573635`

It correctly identified H9 d12 and H9 d19 as the nearest observed temporal anchors.

For the held-out line-by-day combination, it correctly assigned `AN1-1_d19` to the test set:

- RNA: `GSM6573630`
- ATAC: `GSM6573631`

It retained `H9_d19` as a same-day cross-line reference and correctly described the task as generalization to a held-out line-by-day combination, not generalization to an entirely unseen cell line.

In both folds, RNA and ATAC from the held-out biological sample remained in the same split. The answer also recognized that no biological-replicate column was present and did not fabricate replicates.

### 6.3 Issue: overstatement of leakage exclusion

The answer stated that replicate leakage was impossible because each condition occupied a single row. That conclusion is stronger than the manifest supports. A single row per condition rules out detectable duplicate-row or sample-level leakage, but hidden relatedness cannot be assessed without donor, replicate, batch, or derivation identifiers.

A defensible statement would be:

> No duplicate-sample or paired-modality leakage is detectable from the supplied manifest. Biological-replicate, donor, and batch leakage cannot be assessed because the required identifiers are absent.

### 6.4 Assessment

**Verdict: Mostly pass.** Conversational memory and split construction were effective. The principal weakness was uncertainty calibration rather than the split assignments themselves.

## 7. Round 4: Kidney-relevant regulatory evidence

### 7.1 Prompt

> For the kidney-development genes PAX2, LHX1, WT1, HNF4A, GATA3 and SOX9, use IGVF Catalog and ENCODE only. For each gene, return the Catalog gene ID and up to three kidney-relevant regulatory-element or perturbation links. Separate cis from trans links and significant from nonsignificant evidence. Include biosample, effect direction, adjusted p-value when available, and source accession/URL. If kidney-specific evidence is absent, say so rather than substituting evidence from another tissue.

### 7.2 Correct components of the answer

All six Catalog gene identifiers were correct:

| Gene | Catalog gene ID |
|---|---|
| PAX2 | `ENSG00000075891` |
| LHX1 | `ENSG00000273706` |
| WT1 | `ENSG00000184937` |
| HNF4A | `ENSG00000101076` |
| GATA3 | `ENSG00000107485` |
| SOX9 | `ENSG00000125398` |

The answer also recognized that the kidney-relevant rows it reported were primarily ENCODE-rE2G predictions and that these records generally did not contain adjusted P values. It did not directly relabel cardiac or T-cell perturbation data as kidney-specific experimental evidence.

### 7.3 Failure 1: false absence claims caused by truncated retrieval

The underlying gene workflow queried linkage records associated with the genomic region of each gene and requested only 40 rows. Region-based linkage results are not restricted to the queried target gene. Consequently, the first 40 rows can contain other target genes, other tissues, or both.

IGVFagent nevertheless treated the limited result set as sufficient to report that WT1 and HNF4A lacked kidney-specific regulatory records.

For an independent check, the same official Catalog endpoint was queried with a limit of 500. Each query returned exactly 500 rows, so even this expanded retrieval may remain truncated. We then required an exact match to the queried Ensembl gene ID and conservatively retained biological contexts explicitly containing “kidney” or “renal cortex.” Adrenal contexts were excluded.

| Gene | Target-gene rows among the first 500 | Explicit kidney/renal-cortex target-gene rows among the first 500 |
|---|---:|---:|
| PAX2 | 475 | 3 |
| LHX1 | 419 | 44 |
| WT1 | 390 | 28 |
| HNF4A | 346 | 1 |
| GATA3 | 235 | 14 |
| SOX9 | 181 | 6 |

Because all six responses reached the 500-row limit, these are lower-bound audit counts rather than complete database totals. They are nevertheless sufficient to disprove the reported absence of WT1 and HNF4A kidney records.

Examples include:

| Gene | Biological context | Source | Model score | Record class | Adjusted P value | Significant field |
|---|---|---|---:|---|---|---|
| WT1 | left kidney, ENCDO863JOG | `ENCFF899XZX` | 0.7911 | prediction | null | null |
| HNF4A | kidney, ENCDO528BHB | `ENCFF179DXS` | 0.2675 | prediction | null | null |
| GATA3 | right kidney, ENCDO245RIZ | `ENCFF902EDT` | 0.9909 | prediction | null | null |
| SOX9 | kidney, ENCDO163CLX | `ENCFF448PFC` | 0.9302 | prediction | null | null |

The WT1 and HNF4A absence statements were therefore false negatives.

### 7.4 Failure 2: non-target genes were assigned to the queried gene

The relevant endpoint returns element-to-gene edges associated with regulatory elements overlapping the submitted genomic region. A returned row is evidence for the gene identified in that row's `gene` field. It is not automatically evidence for the gene whose locus was used to construct the region query.

In the GATA3 and SOX9 sections, IGVFagent included records whose target gene ID did not match GATA3 or SOX9. These rows were retained “for completeness” and counted among the requested links. This changes the biological meaning of the evidence and can incorrectly imply support for an element-to-query-gene relationship.

The required operation order is:

1. resolve the query symbol to a stable gene ID;
2. retrieve candidate records with complete pagination;
3. require `record.gene == queried_gene_id`;
4. apply tissue or cell-type filtering;
5. separate observed evidence from predictions;
6. rank, deduplicate, and apply the requested top-k limit.

### 7.5 Failure 3: prediction scores were discussed using significance language

The ENCODE-rE2G records in the audit had `class=prediction`, numerical model scores, `p_value_adj=null`, and `significant=null`. The response acknowledged the missing adjusted P values but also suggested that “significant” meant passing the model's own threshold. No threshold field or versioned threshold definition was supplied.

Prediction scores and statistical significance should be represented separately:

- Observed records may be categorized as significant or nonsignificant only when the source provides a significance field or a documented inferential rule.
- Prediction records should be reported with their model name and score.
- Prediction records should not be described as statistically significant, experimentally validated, or nonsignificant when those fields are null.
- Any model threshold must be traceable to an explicit source field or versioned method document.

### 7.6 Failure 4: current-run artefacts were not isolated

The Artefacts panel displayed approximately 70 entries during the six-gene run. It included reports for APOE, TP53, BRCA1, BRCA2, and other entities unrelated to the six requested genes. Multiple paths were also displayed more than once.

This observation demonstrates at least a provenance and attribution defect: the interface presented files not produced for, or requested by, the current run as current-run artefacts. In a shared workspace, the same behavior creates a potential data-isolation concern even when the observed filenames themselves refer only to public analyses.

Artefacts should be associated with an explicit run ID and populated from a run-specific allowlist. The interface should not infer ownership by scanning a shared directory or collecting every path mentioned in prior outputs. Paths should also be normalized before deduplication.

### 7.7 Efficiency

The six-gene metadata request required 21 agent iterations and 101 tool calls, taking approximately 6 min 45 s. Many calls reread or searched generated artefacts rather than performing a single structured retrieval and filter operation.

For this task, the efficient workflow is deterministic after symbol resolution: batch retrieval, pagination, exact target-gene filtering, ontology-based tissue filtering, evidence classification, ranking, and table generation. Caching endpoint schemas and batching gene requests should substantially reduce latency and cost.

### 7.8 Assessment

**Verdict: Fail.** Identifier resolution was accurate, but the central gene–tissue–regulatory-edge mapping was not reliable. The errors are scientifically consequential and could alter candidate prioritization.

## 8. Prioritized recommendations

### P0: Correctness and provenance

1. **Enforce exact target-gene matching.** Every element-to-gene edge must be checked against the resolved query gene ID before it can enter the answer. Add the GATA3 and SOX9 cases to the regression suite.
2. **Prohibit absence claims on truncated results.** Tool outputs should expose `returned_count`, `total_count` or a continuation token, and `truncated`. When exhaustive retrieval is unavailable, the answer must say, “No qualifying record was found among the N retrieved records,” rather than “No evidence exists.”
3. **Add a structured answer verifier.** Each reported gene, element, biosample, score, effect, adjusted P value, significance status, and source URL should map to the same source record. The verifier should reject cross-record field assembly.
4. **Isolate artefacts by run and session.** Assign an immutable run ID, maintain a run-specific artefact allowlist, normalize paths, and deduplicate before display. Unrelated files in a shared workspace must never be attributed to the current run.

### P1: Scientific interpretation

5. **Separate prediction from inferential significance.** An rE2G score is not an adjusted P value. Null significance fields should remain “not provided,” not be converted to significant or nonsignificant labels.
6. **Use ontology-aware tissue filtering.** Kidney and renal cortex should be defined through the appropriate UBERON/CL hierarchy or a documented ontology closure. Substring matching is unsafe because, for example, “adrenal” contains the character sequence “renal.”
7. **Calibrate uncertainty to available metadata.** If donor, replicate, or batch identifiers are absent, state that leakage cannot be assessed at those levels.
8. **Separate validation status from scientific availability.** A valid `false` availability flag is not a schema failure. Use distinct statuses for integrity failure, expected missingness, and design limitation.

### P2: Efficiency and usability

9. **Batch and cache structured queries.** Resolve multiple genes in one plan, retrieve all required pages, cache schema information, and avoid repeatedly reading generated reports.
10. **Apply a declared ranking rule.** A suitable default is: exact target-gene match, target-tissue ontology match, observed evidence before predictions, significant observed evidence before nonsignificant evidence, prediction score ordering, then stable deduplication.
11. **Expose a compact reproducibility record.** Report the endpoint and parameters, retrieval time, database or schema version when available, returned count, truncation status, ranking rule, and source accession for each selected row.
12. **Run a final instruction-compliance check.** Numeric limits, requested columns, tissue restrictions, and evidence-class definitions should be validated before the answer is shown.

## 9. Suggested reproducible issue reports

The findings can be separated into three focused GitHub issues:

1. **P0: Gene-centric KG answer mixes non-target genes and gives false-negative kidney evidence**
   - Reproduce with the six-gene prompt in Section 7.1.
   - Assert exact target-gene equality for every output row.
   - Assert that a truncated query cannot support an unconditional absence statement.

2. **P0: Artefacts panel includes unrelated cross-run files and duplicate paths**
   - Run a six-gene query in a workspace containing prior gene reports.
   - Assert that only artefacts associated with the current run ID are displayed.
   - Assert canonical-path uniqueness.

3. **P1: Uploaded-manifest audit conflates expected missingness with validation failure**
   - Upload the public manifest used in Round 2.
   - Assert that valid `false` availability flags are reported as data limitations rather than integrity failures.

The gene-centric findings provide concrete regression cases for the existing structured-grounding and answer-verification discussion in GitHub issue #18. The successful reuse of the uploaded manifest across turns provides a positive regression result relevant to the conversational-memory discussion in issue #19.

## 10. Supplementary baseline acceptance tests

Before the four MultiFlow-focused rounds, four smaller public-data acceptance tests were conducted on 8 September 2026:

| Test | Main observation | Verdict |
|---|---|---|
| APOE Catalog lookup | Correct Ensembl ID, coordinates, provenance, and returned regulatory-edge fields; edge selection was not biologically ranked | Pass with ranking caveat |
| Deliberately nonexistent identifiers | Correctly returned no gene and no IGVF object rather than fabricating a match | Pass |
| ENCODE K562 H3K27ac lookup | Correctly identified `ENCSR000AKP`; organism was omitted despite being available in nested metadata | Mostly pass |
| Human brain 10x multiome discovery | Correctly separated MeasurementSet and AnalysisSet objects, but returned ten records when the global maximum requested was five | Partial pass |

Additional interface observations from these baseline tests were:

- the chat input appeared approximately 55 s after the model status first showed “Loaded”;
- artefact paths were duplicated across completed runs;
- four concurrent fresh sessions did not expose the chat input within 150 s. This last observation should be reproduced under controlled conditions before it is treated as a confirmed concurrency defect.

## 11. Overall conclusion

IGVFagent showed promising performance for public metadata retrieval, accession-level grounding, small-table auditing, conversational reuse of an uploaded file, and basic paired-modality experimental design. These capabilities could be useful for preparing public datasets and checking analysis inputs.

The current gene-centric regulatory workflow, however, is not yet reliable enough for candidate selection or biological interpretation. Region-level retrieval, incomplete pagination, missing target-gene constraints, and unverified synthesis can produce conclusions that are fluent but biologically incorrect. Run-level artefact attribution also requires attention because the hosted workspace is shared.

The most valuable next step is a deterministic, structured verification layer that sits between retrieval and natural-language synthesis. With exact entity constraints, pagination-aware absence logic, ontology-based context filtering, and run-scoped provenance, the system would retain its useful exploratory interface while substantially improving scientific trustworthiness.

## 12. Public reproducibility resources

- GEO GSE213151: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213151>
- IGVF Catalog gene API example: <https://api.catalogkg.igvf.org/api/genes?name=GATA3&limit=1>
- Existing IGVFagent issue #18: <https://github.com/zhouhufeng/IGVFagent/issues/18>
- Existing IGVFagent issue #19: <https://github.com/zhouhufeng/IGVFagent/issues/19>

No GitHub issue was submitted as part of this evaluation, and no destructive action was performed in the hosted workspace.
