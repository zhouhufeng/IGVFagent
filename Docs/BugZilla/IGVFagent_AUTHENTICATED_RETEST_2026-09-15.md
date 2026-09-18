# IGVFagent Authenticated Hosted Retest

**Retest date:** 15 September 2026  
**Hosted interface:** <https://igvfagent.genohub.org/>  
**Hosted build code:** `bbf9c7f2126e`  
**Backend/model:** Anthropic / Claude Sonnet 5  
**Account:** newly registered and approved Genohub Community user

## Executive summary

The new registration, approval and SSO flow worked. Anonymous visitors were
sent to the access page, the former shared HTTP Basic credentials no longer
opened the application, and the approved account reached the agent normally.

The most important workflow improvement was also confirmed. The exact
six-gene kidney query that previously ran for more than 25 minutes without a
final answer now completed in about 168 seconds. It used the new
`kg_genes_batch` and `rank_artifact` tools, and the ranked GATA3, SOX9, WT1 and
HNF4A results agreed with the independently audited score ordering from the
previous retest.

Three important issues remain:

1. Pre-account history is deliberately visible to every approved user. The
   new account could open and reuse reports and manifests created under the
   old shared login. This is consistent with the current migration policy,
   but it is unsafe unless every legacy run and upload is known to be public.
2. History reuse has no quality gate. The GATA3 query reused an old answer
   whose claimed top-three ranking was already known to be wrong, rather than
   running the new deterministic ranking workflow.
3. The six-gene answer confused `neg_log10_pvalue` with adjusted P value and
   described a value of 240 as `p_adj`. An adjusted P value cannot exceed 1.

No unpublished or internal files were uploaded. All prompts used public gene
symbols and public IGVF Catalog/ENCODE records.

## 1. Access-control checks

| Check | Result |
|---|---|
| Anonymous visit | Redirected to the IGVF Agent access page |
| Old shared HTTP Basic login | No longer grants access |
| Genohub Community SSO | Worked |
| Approved-group membership | Worked; the agent UI loaded as the new user |
| Model initialization | Worked; Claude Sonnet 5 loaded normally |

This confirms admission control, not server-capacity scaling. Registration
limits who can submit work, but it does not by itself bound simultaneous jobs
from approved users. The source includes a one-slot lock for the main
single-cell analysis pipeline, but a general per-user/global scheduler was not
demonstrated by this test.

## 2. HNF4A single-gene control

**Elapsed time:** 182.1 s  
**Agent work:** 22 iterations, 38 tool calls  
**Stop status:** `complete_with_failures`  
**Failed calls:** 1 (`catalog_find_associations`, exit 2 after an unsupported
argument combination)

### What worked

- The answer used the correct target gene, `ENSG00000101076`.
- It returned real kidney ENCODE-rE2G rows from `ENCFF179DXS`.
- It correctly described the rows as predictions with model scores, not as
  statistically significant perturbation results.
- It reported the Catalog traversal as exhaustive.

### Remaining problems

- Most verification was performed by opening the legacy 10 September evidence
  pack and repeatedly slicing or grepping it. The new account was able to read
  that pre-account run.
- The answer first said the kidney-filtered set contained exactly three unique
  pairs, then stated that a fourth kidney element existed. Those statements
  cannot both be true under the stated uniqueness rule.
- One failed tool invocation caused the otherwise useful answer to be marked
  incomplete.
- Several artefacts had identical display names (`gene_HNF4A_report.md` and
  `evidence_pack.json`) from different runs. Even if their full paths differ,
  the panel does not make their provenance clear.

## 3. GATA3 single-gene ranking control

**Elapsed time:** 88.9 s  
**Agent work:** 17 iterations, 20 tool calls  
**Stop status:** `complete`

The scientific row attribution remained correct: all returned rows target
`ENSG00000107485`, and predictions were distinguished from observations.

However, this was not a fresh analysis. All substantive evidence came from the
legacy 10 September run. The agent said it had re-verified the result by
grepping the old manifest and again reported the following as the top three:

| Reported score | Source |
|---:|---|
| 0.9909 | `ENCFF902EDT` |
| 0.9270 | `ENCFF101ISC` |
| 0.9175 | `ENCFF433JOJ` |

These are not the top three rows in the exhaustive target-filtered manifest.
The previous independent full-table audit found:

| True score order | Biosample | Source |
|---:|---|---|
| 0.9999999981 | kidney glomerular epithelial cell | `ENCFF929DQM` |
| 0.9999999974 | renal cortical epithelial cell | `ENCFF829IVH` |
| 0.9999988818 | kidney | `ENCFF517HGF` |

Thus the new ranking tool exists, but history recall can bypass it and revive a
superseded result.

## 4. Exact six-gene kidney query

**Genes:** PAX2, LHX1, WT1, HNF4A, GATA3 and SOX9  
**Elapsed time:** 167.8 s  
**Agent work:** 9 iterations, 26 tool calls  
**Stop status:** `complete`

### Major improvement

The agent used one purpose-built `kg_genes_batch` call followed by
deterministic `rank_artifact` calls for each gene. It no longer drove six
independent `kg_gene` workflows through the LLM loop.

The previous exact prompt was still running after approximately 24 minutes 41
seconds and was terminated without a final answer. The new run completed in
under three minutes, making it more than 8.8 times faster than that lower bound.

The top-three results for the four genes independently audited in the previous
round were now correct:

- GATA3: `ENCFF929DQM`, `ENCFF829IVH`, `ENCFF517HGF`
- SOX9: `ENCFF229GDZ`, `ENCFF947VVS`, `ENCFF050KZD`
- WT1: `ENCFF083XAN`, `ENCFF955XUL`, `ENCFF829IVH`
- HNF4A: the three highest kidney rows from `ENCFF179DXS`

The answer also correctly separated cis ENCODE-rE2G predictions from
non-kidney Perturb-seq/CRISPR-screen edges and did not substitute cardiac or
T-cell observations as kidney evidence.

### Scientific field-label error

The GATA3 trans summary described one result as having “`p_adj` capped at
240.” The underlying GRN tool explicitly returns `neg_log10_pvalue`, not
`p_adj`. A value of 240 is plausible for a capped negative-log P-value score
but impossible for an adjusted P value. This should be treated as a semantic
field-mapping error in the final answer.

## 5. Source-level regression checks

The current public source was cloned independently and its standalone tests
were run:

| Test | Result |
|---|---:|
| Regulatory linkage regression | 46 passed, 0 failed |
| Artefact isolation/deduplication regression | 36 passed, 0 failed |
| History visibility and project sharing checks | all checks passed |

The history tests confirm that post-account work is private by owner unless a
project is shared. They also confirm that all pre-account records are assigned
the empty legacy owner and intentionally remain visible to every signed-in
user. The hosted observation therefore matches the implemented policy rather
than showing a failure of the owner filter.

## 6. Recommendations

### Highest priority

1. **Quarantine legacy history instead of treating it as automatically
   public.** Inventory pre-account sessions and uploads, mark only reviewed
   public benchmark runs as shared, and make all unreviewed legacy material
   admin-only. A shared-password run is not evidence that its inputs were
   public.
2. **Add quality state to cached runs.** History records should carry
   `unreviewed`, `validated`, `failed` or `superseded`. The agent should reuse
   only validated results produced by a compatible workflow/data version.
3. **Validate scientific field semantics before finalization.** Enforce ranges
   and units: `0 <= p_adj <= 1`; `neg_log10_pvalue >= 0`; model scores must be
   labeled separately. Reject or rewrite a final table when a value is placed
   in an incompatible field.

### Workflow reliability

4. Route both single-gene and multi-gene regulatory-evidence questions through
   the same deterministic selection pipeline. HNF4A and GATA3 should not fall
   back to dozens of `read_artifact`, `doc_slice` and `grep_artifacts` calls.
5. Validate tool arguments against the live schema before execution so the
   agent does not consume an iteration on unsupported flags.
6. Show artefacts grouped by run with their relative parent directory, rather
   than displaying several unrelated files under the same basename.

### Capacity control

7. Keep registration/approval as admission control, but add a central job
   scheduler for heavy tools: per-user concurrency limits, a global capacity
   limit, queue position, cancellation, timeouts and stale-job recovery. Apply
   it through the tool executor rather than separately inside individual
   analysis skills.

## Overall assessment

The update is a substantial improvement for the exact multi-gene workflow:
the former non-completing query now finishes in under three minutes and uses a
deterministic ranking operation that returns the correct audited top rows. The
new account system also works as an access gate.

The main risks are now history governance and final-answer semantics. The
legacy visibility policy can expose any pre-account material to every approved
user, cached failed analyses can override newer corrected workflows, and the
GRN answer still mislabeled a negative-log P-value score as adjusted P value.
These should be addressed before inviting a broader user group.
