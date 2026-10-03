# belk2022_epigenetic_regulation

## Paper

**Epigenetic regulation of T cell exhaustion.**
Belk JA, Daniel B, Satpathy AT.
*Nature Immunology* 23:848–860 (2022) · doi:[10.1038/s41590-022-01224-z](https://doi.org/10.1038/s41590-022-01224-z) · PMID 35624210 · PMC10439681 · IGVF0072, UM1HG012076

Resolver confidence: **1.00** (exact PMID lookup).

## Why this benchmark has no analyses to reproduce

This paper is a **review article**, not a primary research report — confirmed two independent ways:

1. **Europe PMC `pubTypeList`** for PMID 35624210 includes `Review` (alongside `Journal Article`).
2. **The paper's own abstract** states outright: *"Here, we review foundational technologies to profile the epigenome at multiple scales... We discuss how these technologies have elucidated the development and epigenetic regulation of exhausted T cells..."*

Consequences that follow directly from that:

- **No Data Availability or Code Availability statement**, no deposited accessions, no authors' repository — there is nothing for `igvfagent bench harvest` to find, and it correctly found 0 accessions / 0 numeric-claim candidates.
- **No IGVF assay route matches** (`igvfagent bench route` returned "no route matched"); the scaffold below was force-generated with the generic `portal_discovery` fallback purely so a benchmark directory could exist, but that route is not actually applicable — `run.sh` exits 77 unconditionally rather than issuing a meaningless, paper-unrelated Portal query.
- **Every quantitative figure the text cites belongs to a different paper**, e.g.:
  - ~4,500 differentially accessible regions in basal-cell-carcinoma T_EX cells — Yost et al./Satpathy lab primary data, not this review's own.
  - 555 chromatin regions partially reversed by anti-PD-L1 — Philip et al. 2017 *Nature*.
  - ~1,200 DNA-methylation events accompanying the T_EX transition — Sen et al. 2016 *Science* / Pauken et al. 2016 *Science*, roughly.
  - 182 regulatory elements changing after antigen removal — Yates et al. 2021 *Sci Immunol*.

  Porting or "verifying" any of these numbers under `belk2022`'s name would misattribute other groups' primary experiments as this review's own reproduction target, which the skill's own "no invented science" rule forbids.

## What this means for `bench score`

```
igvfagent bench score --paper-id belk2022_epigenetic_regulation
  reproduction: incomplete 0/1
```

This is a known, deliberate limit of the scoring framework, not unfinished work: `Benchmarks/concordance.py`'s `judge_coverage()` only lets `controlled_access` / `embargoed` / `not_deposited` blockers count a paper as done without a passing check (`reproduced_except_access`); a `blocker.kind: "other"` analysis is treated as "still work to do" by design (see the comment at `concordance.py:279-282`). There is no "review article, nothing to reproduce" terminal state in the framework yet. The single planned analysis (`review_no_original_analysis`) is blocked with `kind: "other"` and a full documented reason in `expected.json`, which is the most honest representation available — an empty `analyses: []` plan would score `reproduction: unplanned`, which reads as "not yet looked at" and is *less* accurate than a documented, explicit block.

## Recommendation

If IGVF wants primary-data coverage for the T cell exhaustion epigenetics literature this review surveys, the actionable targets are the **primary papers it cites** (each has its own accessions/code): Philip et al. 2017 *Nature* (PMID 28723893), Pauken et al. 2016 *Science*, Sen et al. 2016 *Science*, Ghoneim et al. 2017 *Cell*, and Yates et al. 2021 *Sci Immunol*. Those are candidates for their own `igvf-replicate-paper` runs; this review itself is not.

## Provenance

`provenance.json` holds the full resolve/harvest/route record and every source URL consulted, including the Europe PMC `pubTypeList` lookup that established the Review classification.
