# IGVFagent test with released SC-islet multiome data

Test date: September 15, 2026. Hosted service: https://igvfagent.genohub.org/.

## Summary

We ran three consecutive tests using a released IGVF SC-islet dataset: metadata retrieval, quantitative RNA-matrix QC, and interpretation of what the available inputs support.

The basic data-reading and per-barcode calculations worked. We downloaded the agent's QC table and independently checked all 211,901 rows against the checksum-verified public original. Every barcode, total count, and detected-gene count matched.

The requested workflow nevertheless remained incomplete: IGVFagent could not apply a simple two-column threshold and report the resulting count and medians. Its inspection tool also returned an incorrect whole-matrix maximum, and the follow-up answer contradicted earlier results about control metadata and whether filtering had been completed.

These observations separate a working numerical component from problems in tool selection, tool output, and cross-turn state tracking. This was not a load test, and there is no evidence from this run that server overload caused the failures.

## 1. Dataset and privacy boundary

We selected a public dataset already used in our project, not an unpublished project output. The experimental system is H1-hESC-derived islet organoids, with a three-cytokine treatment for six hours. It is treatment-response 10x multiome, not a CRISPR Perturb-seq experiment.

- AnalysisSet: [IGVFDS2239CRKD](https://data.igvf.org/analysis-sets/IGVFDS2239CRKD/).
- RNA MeasurementSet: [IGVFDS9331KIOE](https://data.igvf.org/measurement-sets/IGVFDS9331KIOE/).
- ATAC MeasurementSet: [IGVFDS6979KCTQ](https://data.igvf.org/measurement-sets/IGVFDS6979KCTQ/).
- Shared biosample: IGVFSM1994SNJP.
- RNA matrix: [IGVFFI6543TDVW](https://data.igvf.org/matrix-files/IGVFFI6543TDVW/), released, unfiltered h5ad.
- File size: 173,829,200 bytes.
- Public MD5: `057a6b74cd30c966bea968058cd0306f`.

The live metadata describes interferon gamma at 10 ng/mL, interleukin-1 beta at 0.5 ng/mL, and tumor necrosis factor at 1 ng/mL, each for six hours. RNA and ATAC inputs refer to the same biosample. Actual barcode-level pairing was not tested because no ATAC data were loaded.

Only public accessions and test instructions were submitted. No local data files were uploaded. The server fetched the released RNA matrix itself. No unpublished predictions, rankings, performance results, internal processed matrices, or private donor-level files were supplied. Our independent reference used an existing local copy of the original public file, after verifying its size and MD5 against live IGVF metadata.

The quantitative test allowed only this h5ad, with a 200 MB data-download limit. BAM, FASTQ, fragment, archive, and additional dataset downloads were prohibited, as were software installation and unrelated downstream analyses. The trace records one failed HTTP 400 download attempt followed by one successful 173,829,200-byte fetch. It does not show a successful download of any other data file. We did not independently meter all server network traffic.

## 2. Test execution

All three prompts were submitted sequentially in one authenticated browser conversation. The UI reported build code `bbf9c7f2126e`, backend `anthropic`, and model `claude-sonnet-5`. The build code is recorded as displayed; it was not established to be a public Git commit hash.

| Test | Browser-observed elapsed time | Iterations / tool calls | Outcome |
|---|---:|---:|---|
| Metadata and experimental design | 58.7 seconds | 4 / 8 | Main requested facts correct; control-selection wording needs care |
| Bounded RNA-matrix QC | 217.3 seconds | 13 / 16 | File read and per-barcode QC correct; threshold summary not completed; incorrect maximum reported |
| Follow-up interpretation | 28.7 seconds | 1 / 0 | Main feasibility decisions correct; two factual contradictions |

Elapsed times include browser polling and are not isolated server-compute benchmarks. A `completed: true` field in the local test harness means the agent turn ended and was captured, not that the scientific task passed. The QC turn explicitly ended with `max_iterations_wrapped` and six failed tool calls.

Full prompts are preserved in [prompts.json](prompts.json). Unedited agent answers and visible execution traces are saved in the `metadata`, `matrix_qc`, and `interpretation` subdirectories.

## 3. Independent quantitative check

The reference calculation operated on the entire sparse matrix, without making a dense barcode-by-gene array. Per-barcode total counts were row sums of X; detected genes were counts of entries greater than zero. The requested rule was:

`total_counts >= 1000 AND n_genes_by_counts >= 200`

The thresholds are an explicit test specification, not a recommended cell-calling method. Passing barcodes must not be described as validated cells.

| Quantity | Independent reference | IGVFagent result |
|---|---:|---|
| Original barcode rows | 211,901 | 211,901 — correct |
| Gene columns | 62,757 | 62,757 — correct |
| Median total counts before filtering | 1 | 1 — correct |
| Median detected genes before filtering | 1 | 1 — correct |
| Barcodes passing the specified rule | 1,858 | Not computed in the agent answer |
| Median total counts after filtering | 3,443 | Not computed in the agent answer |
| Median detected genes after filtering | 1,857 | Not computed in the agent answer |
| Maximum value in the full X matrix | 3,724 | 14 — incorrect as a whole-matrix maximum |

Additional structure checks matched: X uses float64 storage but contains nonnegative integer-valued counts; the available layers are `ambiguous`, `mature`, and `nascent`; no separate `.raw` object or ATAC embedding is present. `obs` and `var` have no annotation columns. Their index names are `barcode` and `gene_id`, respectively—index fields should not be mistaken for annotation columns.

The downloaded hosted QC table contains 211,901 rows in exactly the same barcode order as the original matrix. There were **zero total-count mismatches and zero detected-gene-count mismatches**. Independently applying the requested filter to that hosted table also gives 1,858 passing barcodes and post-filter medians of 3,443 and 1,857. These final aggregations were performed by the tester, not by IGVFagent.

Evidence: [independent reference](independent_qc_reference.json), [row-by-row comparison](hosted_artifact_comparison.json), and [hosted QC table](hosted_artifacts/20260915_200617_IGVFFI6543TDVW_qc_rna_qc.tsv).

## 4. Findings

### A. Correct per-barcode QC could not be turned into a simple threshold summary

`share_rna_qc` successfully produced the correct per-barcode table. The remaining operation required no additional input data: select rows satisfying two numeric conditions, count them, and calculate medians.

The execution trace shows:

1. An attempted custom sparse-QC tool was rejected as overlapping with existing tools.
2. The agent then successfully ran `share_rna_qc`.
3. A `warehouse_query` aggregation failed because the warehouse was not initialized.
4. A small TSV aggregation tool was rejected twice as overlapping with `enrich_gsea`, based on the shared words `stats` and `tsv`.
5. The detailed UI trace recorded repeated all-failed iterations at iteration 13, followed by wrap-up. Although the final stop label was `max_iterations_wrapped`, the trace does not show exhaustion of all 25 advertised iterations.

This is a workflow/capability-routing failure in the tested deployment, not evidence of invalid input counts. The final QC answer appropriately admitted that the threshold summary was incomplete and did not invent the missing values.

Suggested change: provide a bounded table-filtering and aggregation operation that does not require a separately initialized warehouse or dynamic skill creation. Tool-overlap checks should compare actual operations and input/output contracts, not just generic words. Existing safety checks should remain in place; this recommendation is not to disable them.

### B. The inspection tool returned 14 as X_max, but the full-matrix maximum is 3,724

The raw `h5ad_inspect` tool output itself contains `X_max: 14.0`, and the final answer repeats it without a sampling qualifier. The independently verified public matrix contains a value of 3,724 at zero-based row 21,043, column 34,645:

- Barcode: `ACGCCTTTCACCGGTA_IGVFSM1994SNJP`.
- Gene ID: `ENSG00000251562.10`.

The maximum in the first 100 rows happens to be 14. This supports a sampling-related explanation, but the exact implementation of the deployed inspector was not inspected, so that mechanism remains an inference. What is established is that the tool output is not the maximum of the complete matrix. This is not simply an arithmetic mistake introduced by the language model.

Suggested change: return explicit scope fields, such as `sampled`, `rows_examined`, and `sample_X_max`, or compute the full sparse-matrix reduction. The answer should not present a sampled statistic as a whole-dataset statistic.

### C. The follow-up answer lost the distinction between completed and incomplete operations

The QC answer clearly said it had not completed the threshold step. In the next turn, however, the agent stated that the barcode QC thresholds had been computed.

The underlying per-barcode metrics did exist, but the requested pass/fail aggregation did not. The follow-up therefore overstated what had been completed.

Suggested change: persist operation-level status—such as `metrics_completed`, `filter_not_run`, and `summary_missing`—and use it when answering follow-up questions. A failed or partial step should never become completed through conversational summarization alone.

### D. The follow-up incorrectly denied the presence of control links in the metadata

The first answer identified `control_file_sets`. The follow-up then said that neither input MeasurementSet named an untreated counterpart. That is inconsistent with both the earlier answer and the live metadata.

For example, the RNA record lists:

- `IGVFDS0270LKNE`: alias includes `control_rep1_10x-multiome_6hr`.
- `IGVFDS3434LMTM`: alias includes `control_rep2_10x-multiome_6hr`.
- `IGVFDS2052PCVB`: alias includes `control_rep3_10x-multiome_0hr`.

The correct distinction is: **control accessions are linked in the metadata, but their matrices were not downloaded and are not included in the current single-AnalysisSet input.** The conclusion that the loaded data alone cannot support treatment-versus-control DE remains correct.

There is also a smaller issue in the initial answer: it characterized the linked controls collectively as 0-hour controls and suggested a 0-hour accession as a next step, although the links include 6-hour controls. No incorrect DE comparison was actually run. Future control selection should verify treatment, time point, biological replicate, and batch; a control label alone is insufficient.

Suggested change: separately track `linked`, `retrieved_metadata`, and `loaded_matrix` status for controls. Prefer an explicitly justified time-matched comparison when appropriate, rather than silently mixing treatment and time effects.

## 5. Other failures and what can be attributed

- **File-accession routing:** `raw_pipeline_plan` said IGVFFI6543TDVW was not returned by the Portal, even though `explain_dataset` resolved the same released file and a subsequent fetch succeeded. The router had automatically passed a MatrixFile accession into the planning path. This suggests a resource-type/routing mismatch; the trace does not support labeling the accession invalid or requesting new private-data credentials. The precise implementation cause was not established.
- **Initial download error:** the first download path recorded HTTP 400, while `igvf_fetch_file` subsequently succeeded. Recovery worked. A client/redirect issue, upstream response, or another transport detail cannot be distinguished from the captured record alone.
- **Submission-validation failure:** `igvf_submit_validate` was called on the h5ad and returned a traceback. The captured traceback ends before the terminal exception, so the exact cause cannot be determined. Submission validation was not part of the requested QC task.
- **Failure summaries:** the condensed `failed_calls` entries for tool authoring foregrounded a separate tool-shadowing warning. The full tool outputs show that the actual refusals were overlap decisions. Diagnostic summaries should preserve the decisive error, not merely the first stderr line.
- **Checksum:** the agent explicitly admitted that it had not recomputed the downloaded file's MD5. Our local reference was checksum-verified, and all hosted per-barcode counts match, but those checks do not establish a byte-for-byte MD5 match for the server's copy. That integrity check remains unverified on the hosted side.

The numerical/reference disagreement and conversation-state contradictions are established findings. The exact causes of the HTTP 400 and truncated validation traceback are not. No memory-exhaustion error, server-crash evidence, or load-dependent failure was observed.

## 6. What worked

- Correct classification as cytokine-treatment multiome, not CRISPR Perturb-seq.
- Correct main metadata: experimental system, six-hour treatment, RNA and ATAC accessions, shared biosample, and released output inventory.
- Successful recovery to download and read the public h5ad.
- Correct sparse per-barcode total counts and detected-gene counts across every row checked.
- Honest disclosure of the incomplete quantitative step in its immediate answer.
- Correct distinction between a paired-multiome dataset and an RNA-only loaded h5ad.
- No additional tool calls in the final interpretation turn, as requested.
- No invented guide assignments, DE results, clustering, or regulatory networks.

## 7. Suggested regression checks

1. Run this exact matrix through the specified QC rule and assert 1,858 passing barcodes, with post-filter medians 3,443 and 1,857. Label them as QC-passing barcodes, not called cells.
2. Assert full-matrix X_max = 3,724, or explicitly identify a sampled statistic and its scope.
3. After a deliberately incomplete filtering step, ask what was completed; the answer must retain the incomplete status.
4. Present metadata containing linked controls without loading control matrices; the answer must preserve both facts.
5. Test the same valid accession through discovery, planning, and download. A MatrixFile should resolve through its parent AnalysisSet when required, or receive a clear unsupported-resource message.
6. Make simple, read-only aggregation usable without requiring warehouse initialization. A TSV threshold calculation must not be classified as GSEA merely because both mention statistics and TSV files.

These checks would test execution and evidence consistency more directly than adding more domain-specific skills alone.

## 8. Evidence and limitations

Primary local evidence:

- [Live public Portal records](public_portal_records.json).
- [Exact prompts](prompts.json).
- [Metadata answer](metadata/report.md).
- [QC answer](matrix_qc/report.md) and [visible execution trace](matrix_qc/page.txt).
- [Interpretation answer](interpretation/report.md).
- [Raw QC tool outputs](hosted_artifacts/tool_outputs.json).
- [Independent QC reference](independent_qc_reference.json).
- [Hosted-table comparison](hosted_artifact_comparison.json).

The source records were fetched from the official API at `https://api.data.igvf.org`, using the four accession-specific endpoints listed above. Browser screenshots, run timing records, and the downloaded QC transcript are also retained locally.

This was one dataset and one three-turn test sequence, not a repeated performance benchmark or a broad estimate of error frequency. No ATAC matrix, controls, or additional biological replicates were analyzed. Therefore, the test does not validate multimodal integration, DE, cell-type annotation, CRISPR guide assignment, or biological perturbation-effect estimation. It also does not establish user-isolation guarantees. No issue, email, or other message was posted to the maintainers during this test.
