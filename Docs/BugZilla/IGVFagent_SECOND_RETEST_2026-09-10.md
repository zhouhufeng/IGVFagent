# IGVFagent Hosted Regression Retest After the Second Update

**Retest date:** 10 September 2026  
**Hosted interface:** v0.2.9  
**Hosted code fingerprint:** `b677a58977ca`  
**Backend/model:** Anthropic / Claude Sonnet 5  

## Executive summary

The second update fixes the most important scientific errors from the original test. The hosted agent now finds kidney-related WT1 and HNF4A rE2G predictions, filters element-to-gene rows by the requested target gene, preserves the fields of each reported link within one source record, and distinguishes computational predictions from experimental perturbation evidence.

Two substantive problems remain:

1. The agent repeatedly describes selected rows as the highest-scoring or “top 3” records even though it obtained them from bounded `grep`/file excerpts rather than sorting the complete manifest. The reported rows are real and correctly attributed, but the ranking claims are false for GATA3, SOX9 and WT1.
2. The exact original six-gene query did not complete within approximately 25 minutes. Several single-gene queries also required 3–6 minutes, and two reached the 25-iteration limit.

The cross-run Artefacts leakage appears fixed. Duplicate rendering of current-run files remains.

## Data-safety boundary

No unpublished results, internal project files, candidate rankings, performance measurements, donor-level information, or private data were uploaded. The hosted tests used gene symbols and public IGVF Catalog/ENCODE records only. No file was uploaded in this retest.

## 1. Deployment verification

The hosted UI now displays a content fingerprint in the sidebar and run footer. The tested server reported:

```text
v0.2.9
code b677a58977ca
```

The same fingerprint is reproduced by repository states at commits `00a3275` and `408b748`. Therefore, the hosted code includes the relevant fixes introduced by:

- `a948ca5`: pagination, exact target-gene filtering, prediction/observation separation and safe absence wording;
- `8b4381d`: deployed-code fingerprint;
- `2a67207`: run-scoped Artefacts;
- `5569e9f`: record-level verification against cross-record field assembly.

The current `main` branch has a different top-level Python content hash (`04cf4dd9fcb0` at `31aff70`), but the subsequent differences are documentation, attribution text and comments rather than changes to the regulatory-evidence logic tested here.

## 2. Local regression tests

The latest public source was cloned independently and the relevant regression suites were run.

| Test suite | Result |
|---|---:|
| `Benchmarks/test_kg_linkage.py` | 46 passed, 0 failed |
| `Benchmarks/test_artefact_isolation.py` | 28 passed, 0 failed |

The linkage suite covers Catalog pagination, exact target-gene filtering, observation/prediction separation, safe handling of truncated retrieval and rejection of cross-record field assembly. The Artefacts suite covers run scoping, path normalization, cross-run isolation, upload isolation and disabling executable extension uploads on the shared deployment.

## 3. Hosted tests

### 3.1 Exact six-gene regression prompt

The original public-data prompt for PAX2, LHX1, WT1, HNF4A, GATA3 and SOX9 was submitted without modification.

- Chat ready: approximately 7.5 seconds
- Model ready: approximately 11.7 seconds from navigation
- Status after 24 minutes 41 seconds: still running
- Final outcome: the browser test was terminated at approximately 25 minutes without a completed answer
- Saved completed run: none

This run cannot be used to score the six-gene scientific answer, because no final answer was produced. It is, however, a performance and robustness failure for the same query that previously completed in approximately 230 seconds.

### 3.2 HNF4A positive-control query

**Hosted run:** approximately 197 seconds, 25 iterations, 35 tool calls, `stop=max_iterations_wrapped`.

Independent Catalog audit:

- 4,367 region rows retrieved over 9 pages;
- retrieval exhausted at the Catalog layer;
- 2,398 rows retained for target `ENSG00000101076`;
- 1,969 rows for other target genes removed;
- 5 explicit kidney rows after excluding adrenal matches;
- all 2,398 retained rows were predictions, with null adjusted-P-value and significance fields.

The hosted answer correctly reported three real kidney HNF4A records from `ENCFF179DXS`, including the promoter score `0.9999957863` and two intragenic scores `0.2674612296` and `0.6511874483`. Every reported row targeted `ENSG00000101076`; no fields were borrowed from another record. The answer also correctly stated that these are ENCODE-rE2G predictions rather than statistically significant perturbation results.

Remaining issue: the final answer called retrieval “TRUNCATED” because the model only inspected part of the generated file. The underlying Catalog retrieval was exhaustive. These are two different statuses and should be reported separately.

### 3.3 GATA3 target-attribution query

**Hosted run:** approximately 206 seconds, 16 iterations, 30 tool calls, `stop=complete`.

Independent Catalog audit:

- 14,783 region rows retrieved over 30 pages;
- 7,934 rows retained for target `ENSG00000107485`;
- 6,849 rows for other genes removed;
- 490 kidney/renal rows after excluding adrenal matches;
- all 7,934 retained linkage rows were predictions.

The three hosted rows were all genuine single records with the correct GATA3 target, biosample, score and ENCODE file accession. The previous wrong-gene attribution problem is therefore fixed for this test.

However, the answer stated that it ranked the exhaustive kidney/renal set by score and returned the top three. It returned scores `0.9909405`, `0.9270458` and `0.9175458`. These are not the top three under its declared rule. The complete target-filtered manifest contains higher-scoring kidney records, including:

| True score order (examples) | Biosample | Source |
|---:|---|---|
| 0.9999999981 | kidney glomerular epithelial cell | ENCFF929DQM |
| 0.9999999974 | renal cortical epithelial cell | ENCFF829IVH |
| 0.9999988818 | kidney | ENCFF517HGF |

The likely cause is that the agent used bounded `grep_artifacts` matches and file excerpts rather than a deterministic sort of the full CSV.

### 3.4 SOX9 target-attribution and perturbation query

**Hosted run:** approximately 86 seconds, 7 iterations, 11 tool calls, `stop=complete`.

Independent Catalog audit:

- 6,396 region rows retrieved over 13 pages;
- 2,934 rows retained for target `ENSG00000125398`;
- 3,462 rows for other genes removed;
- 2,933 predictions and 1 observed-data record;
- 158 kidney/renal prediction rows after excluding adrenal matches.

The hosted answer correctly kept all three predicted rows on the SOX9 target and preserved their biosample, score and accession within individual source records. It also correctly reported the single observed Perturb-seq record: cardiac muscle context, `log2FC=-1.235099346`, adjusted P value `2.8663e-07`, `significant=true`, source `IGVFFI0830FXFI`. It did not mislabel that non-kidney perturbation result as kidney evidence.

The ranking claim was again wrong. The answer said the three predictions were chosen as the highest-scoring records per distinct source accession, but it returned approximately `0.9999893`, `0.9302299` and `0.8048358`. Higher-scoring kidney records from distinct accessions include:

| True score order (examples) | Biosample | Source |
|---:|---|---|
| 0.9999993987 | right kidney | ENCFF229GDZ |
| 0.9999990862 | right kidney | ENCFF947VVS |
| 0.9999988238 | right kidney | ENCFF050KZD |

Thus, record integrity is fixed, but evidence selection/ranking is not yet reliable.

### 3.5 WT1 positive-control and ranking query

**Hosted run:** approximately 5 minutes 48 seconds, 25 iterations, 34 tool calls, `stop=max_iterations_wrapped`.

Independent Catalog audit:

- 4,741 region rows retrieved over 10 pages;
- retrieval exhausted at the Catalog layer;
- 4,458 rows retained for target `ENSG00000184937`;
- 283 rows for other genes removed;
- 4,457 predictions and 1 observed-data record;
- 402 records matched the answer's kidney/renal/nephron/glomerular/podocyte terms after adrenal exclusion.

The answer correctly returned only WT1-targeted rows and correctly described the single non-kidney observed Perturb-seq record (`IGVFFI0830FXFI`). The top-ranked prediction (`ENCFF083XAN`, score `0.9999999974`) was correct.

Ranks 2 and 3 were not the next highest records under the stated score-descending rule. Higher records existed for `ENCFF955XUL` (`0.9999999969`) and `ENCFF829IVH` (`0.9999999954`). The answer partly acknowledged that its file inspection was incomplete, but still labeled the table “Top 3,” which is internally inconsistent.

The answer also incorrectly described Catalog retrieval as truncated. The API traversal itself was exhaustive; only the model's subsequent `read_artifact` view was truncated. Finally, its suggested follow-up command used unsupported arguments (`--filter` and `--export_full`) and does not match the actual `kg_traversal_skill.py` CLI.

## 4. Artefacts panel

No unrelated APOE, TP53, BRCA1 or BRCA2 files appeared in the tested runs. Cross-run artefact leakage therefore appears fixed.

Current-run duplicates remain. For example, the GATA3 panel rendered `gene_GATA3_report.md` and `evidence_pack.json` twice, and several manifest files appeared both as direct artefacts and as children of an expanded run directory. The SOX9 panel showed the same pattern. The displayed count (for example, “Artefacts (8)”) also did not match the much longer visible list after expansion.

The source currently deduplicates before directory expansion. Expanded children are converted to absolute paths while already-reported files may remain relative, allowing the same resolved file to reappear. A final resolved-path deduplication after `_expand_artefact_dirs()` should fix this. A regression case should include both a run directory and one of its explicitly reported child files.

## 5. What is now fixed

- The deployed build can be identified from saved runs.
- HNF4A and WT1 kidney rE2G predictions are no longer falsely reported as absent.
- Nearby-element rows are filtered by the row's actual target gene.
- GATA3 and SOX9 test rows no longer point to other genes.
- Prediction scores are not described as adjusted P values or statistical significance.
- Experimental Perturb-seq rows are separated from computational rE2G predictions.
- Cross-run Artefacts leakage was not reproduced.
- Relevant local regression suites pass.

## 6. Recommended changes

### High priority

1. **Replace grep-based ranking with a deterministic table operation.** Filter the complete target-specific manifest by an explicit tissue ontology rule, sort numerically by score and then apply the requested top-k/deduplication rule. The selected rows should be written to a small `selected_evidence.csv` file.
2. **Verify the final selected table, not only the source manifest.** Run `verify_table()` on every final row and reject the answer if any gene, element, biosample, score, significance field or accession does not map to one source record.
3. **Separate retrieval completeness from inspection completeness.** Report, for example, `Catalog retrieval: exhaustive (10 pages)` and `Agent inspection: partial due to display/read limit`. Do not collapse these into one “truncated” label.
4. **Optimize multi-gene queries.** Avoid repeatedly reading multi-megabyte manifests through the LLM. Add a structured `select_regulatory_evidence` tool that performs target filtering, tissue filtering, sorting and top-k selection in one call and returns a compact result.

### Medium priority

5. Deduplicate Artefacts after directory expansion using normalized resolved paths.
6. Validate suggested CLI commands against the live argument parser before displaying them.
7. When the iteration limit is reached, return a clearly labeled partial answer with per-gene completion status rather than a generic wrapped completion message.
8. Set `IGVF_GIT_SHA` during deployment so saved reports show both the content fingerprint and the originating Git commit.

## Overall assessment

The update resolves the original high-severity attribution and false-absence errors. The system is now substantially safer for regulatory-evidence lookup because the returned records themselves are correctly tied to the requested genes and are labeled as predictions versus observations.

It is not yet reliable for ranked evidence selection or larger multi-gene requests. A user can trust that the reported GATA3/SOX9/WT1/HNF4A rows exist and belong to those genes in the tested cases, but should not trust a claim that they are the highest-scoring records unless the ranking is computed programmatically from the full manifest.
