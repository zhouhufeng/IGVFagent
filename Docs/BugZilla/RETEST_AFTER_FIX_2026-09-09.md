# IGVFagent Retest After the 9 September Fix

**Retest date:** 9 September 2026  
**Hosted interface version shown:** v0.2.9  
**Public-data prompt:** The same six-gene kidney regulatory-evidence prompt used in the original evaluation

## Result

The source-code fix is technically sound, but the hosted site did not appear to be running that fix during this retest.

## Audit of the retest procedure

The main scientific conclusion is supported, but three methodological details should be stated explicitly:

1. Two browser-automation attempts ended before a final answer was produced. They were excluded from the evaluation. The counted run is the saved 18:22 UTC run whose transcript reports `stop_reason=complete`, 12 iterations, and 76 tool calls.
2. The saved transcript contains the exact original six-gene prompt. This was therefore a like-for-like regression test rather than a modified prompt that might change tool routing.
3. The independent Catalog check concerns **predicted regulatory links**, not experimentally validated perturbation effects. The relevant WT1 and HNF4A rows have `class=prediction`, a model score, and null adjusted-P-value and significance fields. The hosted answer remains incorrect because it claimed that no kidney regulatory-element links existed at all. A correct answer could say that no kidney perturbation evidence was found while separately reporting the kidney rE2G predictions.

An initial exploratory tissue filter used the substring `renal`, which can also match `adrenal`. That issue was identified and corrected before the reported audit counts were finalized. The final checks require an explicit `kidney` or `renal cortex` context and exclude adrenal records.

The Artefacts observation is based on one completed retest. No unrelated APOE, TP53, or BRCA files appeared in that run, so cross-run attribution appears improved, but one run is not enough to call the issue permanently resolved.

## Source-code verification

Commit `a948ca5` directly addresses the reported linkage problems. It adds:

- pagination using the Catalog API's working `page=` parameter;
- explicit truncation metadata and safe absence wording;
- exact target-gene filtering;
- rejection of API error bodies as data rows;
- separate counts and outputs for observations and predictions;
- 28 regression cases covering these behaviors.

All 28 regression cases passed locally.

A live Catalog query run through the latest local code also behaved correctly for HNF4A:

- 4,367 region-level rows retrieved over nine pages;
- 2,398 rows retained for target gene `ENSG00000101076`;
- 1,969 rows for other target genes dropped;
- five retained rows had an explicit kidney or renal-cortex context;
- zero wrong-target rows remained;
- retrieval was reported as exhausted, not truncated;
- all 2,398 rows were correctly labeled as predictions with null adjusted P values and null significance fields.

## Hosted regression run

The exact original six-gene prompt was submitted again through the hosted interface.

- Backend/model: Anthropic, Claude Sonnet 5
- Iterations: 12
- Tool calls: 76
- Approximate completion time: 230 seconds
- Stop reason: complete

The final hosted answer still stated:

> none of the six genes has kidney-tissue regulatory-element or perturbation evidence in the IGVF Catalog

It again reported no kidney-specific regulatory-element evidence for WT1 and HNF4A, so the main scientific correctness test still failed.

## Evidence that the hosted container is using the earlier implementation

The hosted transcript showed that its generated gene reports still contained the earlier heading:

> Enhancer-gene linkage in this region

The updated code instead produces a heading that explicitly states:

> edges whose TARGET GENE is `<gene_id>`

The hosted report also lacked the new retrieval summary containing target-gene row counts, dropped non-target rows, page count, and truncation status. Its tool sequence subsequently read `regulatory_elements.csv`, but did not read the new `linkage_predictions.csv` output needed to answer the regulatory-linkage question.

These differences indicate that the current hosted container was not built from commit `a948ca5`, even though the repository contains the fix.

## Improvements observed

- The chat input appeared after approximately 7–8 seconds, compared with roughly 55 seconds in the earlier test.
- The model loaded in approximately four seconds.
- No unrelated APOE, TP53, BRCA1, or BRCA2 artefacts were observed in the current run's produced-file list. Run-level artefact attribution therefore appears improved in this single retest, although the report still linked 83 current-run files.
- The six-gene run decreased from 21 iterations and 101 tool calls to 12 iterations and 76 tool calls, with approximate runtime decreasing from 405 to 230 seconds.

## Recommended next step

1. Rebuild and recreate the hosted container from commit `a948ca5` or later rather than only pulling the repository.
2. Verify the running container's code hash after deployment.
3. Display the deployed Git commit SHA in the UI alongside v0.2.9.
4. Rerun the exact six-gene prompt and confirm that the final answer reads the target-filtered `linkage_predictions.csv` files.

One small source-level issue remains in the latest local code: the generated report's `Ensembl ID` line is blank because the renderer does not fall back to the metadata `_id` field, although the correct `_id` is used internally for target-gene filtering. This does not invalidate the linkage fix but should be corrected for report completeness.
