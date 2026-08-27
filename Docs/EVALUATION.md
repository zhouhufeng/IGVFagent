# Evaluating IGVFagent — three layers, measured separately

An agent system can fail in three independent ways, and a single aggregate
score hides which one occurred:

1. **Skill correctness** — does `open4gene` compute what the paper computed?
2. **Planning correctness** — given a question, does the agent choose the right
   skills, in a workable order, with correct arguments, and recover when one
   fails?
3. **Conclusion validity** — is the biological claim in the final answer
   actually supported by the artefacts produced?

A system can pass (1) completely and still fail (2) and (3). These are
therefore reported separately and never averaged together.

---

## Tier 1 — Skill correctness *(implemented)*

`Benchmarks/<paper_id>/` — twenty-one papers, each with `run.sh`,
`expected.json` machine-readable checks, and `OPERATIONS.md`. Scored by
`Benchmarks/concordance.py`.

**What it establishes.** The implementations reproduce published results from
public data. For the five full single-cell / multiome reproductions this
includes downloading the data and running the real analytical chain.

**What it does not establish — and this is a scope statement, not a caveat.**
Every `run.sh` invokes skills directly as shell commands:

```bash
python3 Scripts/figshare_skill.py download --id 26299093 ...
```

There is **no model in the loop**. Tier 1 measures library correctness. It is
silent on whether the agent could have selected those skills unaided.

**Known selection bias.** Papers were included partly because their analysis is
expressible in the existing command vocabulary, which makes the suite
in-distribution by construction. It is evidence that the covered methods are
right, not that coverage is representative.

---

## Tier 2 — Planning and tool selection *(implemented)*

Measures the decision layer with the skills held fixed and known-good.

**Unit.** A question plus a *gold plan*: the set of tool calls a domain expert
considers necessary and sufficient, with which arguments matter and which are
free.

**Metrics.**

| Metric | Definition |
|---|---|
| Tool-selection precision / recall | called vs. gold-required tools |
| Argument correctness | fraction of gold-required arguments bound correctly |
| Order validity | were dependencies respected (data pulled before analysed) |
| Recovery rate | after an induced tool failure, does the run still reach a valid answer |
| Redundancy | calls beyond the gold set that add no evidence |

Recall matters more than precision: a missed essential tool invalidates the
conclusion, while an extra call mostly costs tokens.

**Failure injection.** Tier 2 must include deliberately failing tools —
network error, empty result, malformed argument — because *recovery* is a
planning behaviour the current suite never exercises. A run that cannot
proceed when an archive returns 429 has a planning defect, not a skill defect.

**Held-out protocol.** To address the in-distribution problem, Tier 2 questions
must be written by people who did **not** build the system and who are not
shown the tool list. Questions authored by builders are, unavoidably, questions
the builders know are answerable. Held-out and builder-authored questions are
reported separately.

---

## Tier 3 — Conclusion validity *(implemented)*

Tier 3 asks the only question a biologist cares about: **is the answer right?**

**Why the current metrics cannot answer it.** `Scripts/_eval.py` computes
`jaccard_artefacts` (do two runs produce the same files) and `tfidf_cosine`
(does the prose look similar). Both are *agreement* measures. Two runs can
agree perfectly and both be wrong; two runs can produce different file sets and
reach the same correct conclusion. Neither number is evidence of correctness.

**Verifier agent.** An independent model instance receives the question, the
final answer, and the produced artefacts — but *not* the reasoning trace — and
judges:

| Judgement | Question |
|---|---|
| Supported | Is each factual claim traceable to a specific artefact? |
| Unsupported | Which claims have no artefact backing? |
| Contradicted | Does any artefact contradict a claim? |
| Equivalence | Do two runs reach biologically equivalent conclusions, independent of wording and file names? |

Equivalence is the metric that replaces prose similarity: "BRCA1 is associated
with Fanconi anemia" and "BRCA1 variants cause FA complementation group S" are
textually dissimilar and biologically concordant. Jaccard and cosine both score
that pair poorly; a verifier should score it as agreement.

**Independence requirement.** The verifier must not be the same instance that
produced the answer, must not see the reasoning trace, and its own agreement
with human expert judgement should be reported — an unvalidated judge only
moves the problem.

---

## Reporting

Results are reported per tier and never collapsed into one number. Free-form
planning and pinned workflows are reported as **separate conditions**, because
the difference between them is a finding rather than noise: pinned workflows
reproduce artefacts exactly, free-form planning does not, and that gap
quantifies the cost of flexibility. Claiming a single "reproducibility" figure
across both would hide the very thing worth measuring.

See [`Docs/THREAT_MODEL.md`](THREAT_MODEL.md) for the parallel distinction
between auditability, repeatability, and validity.

---

## Running the tiers

```bash
igvfagent eval-tiers list                       # the Tier-2 cases
igvfagent eval-tiers tier2 --all                # planning only
igvfagent eval-tiers tier2 --all --inject-failures
igvfagent eval-tiers tier3 --case <case.json> --verifier-model claude-opus-5
```

Cases live in `Benchmarks/tier2/*.json`; the schema is in that directory's
README. Each existing case encodes a failure actually observed in this system,
so the suite is a regression test rather than a hypothetical.

Use a **different** model for the verifier than the one under test. A model
asked to check its own output is a weak judge, and the verifier never sees the
reasoning trace precisely so it cannot be led by the chain that produced an
error.

## First results

Four cases, Claude Sonnet 5 planning, Opus 5 verifying.

**Planning (Tier 2).** All four reached the required tools — recall 1.0
throughout. Precision varies sharply: 1.0 on the focused cases, 0.33 and 0.10
on the broader ones, where the agent made ten calls to satisfy a one-call
plan. Precision is reported but does not gate the verdict, since redundant
calls cost tokens rather than correctness.

**Recovery is the weak spot, and it is a real finding.** With a load-bearing
tool forced to fail:

| induced failure | recovered |
|---|---|
| `explain_dataset` | yes |
| `catalog_get_entity` | yes |
| `read_artifact` | **no** — retried 5×, ended `max_iterations_wrapped` |
| `annotate_variant_list` | **no** — retried 3×, ended `max_iterations_wrapped` |

When an *alternative route exists*, the agent finds it. When the failing tool
is the only route, it retries until the budget is exhausted instead of
reporting the failure and concluding. Nothing in the Tier-1 suite exercises
this path, because Tier 1 never lets a tool fail.

**Conclusion validity (Tier 3) caught an auditability gap.** On the ENCODE
identity case the verifier returned `PARTIALLY_SUPPORTED`, flagging *assembly*
and *status* as unsupported. Both claims were **true**, and the model had seen
them — but they were printed to stdout and never written into any artefact, so
they were unevidenced in the record. Persisting the identity block into the
report moved the verdict to `VALID` (7 supported, 0 contradicted).

That is the distinction the threat model draws, caught mechanically: an answer
can be correct and still not auditable. It is also the clearest argument for
Tier 3 — no agreement metric would have found it, because a second run would
have produced the same unevidenced-but-correct claim and scored perfect
agreement.

The one remaining unsupported claim, *"GM12878 is a lymphoblastoid cell
line"*, is correct background knowledge absent from the artefacts. Flagging it
is the desired behaviour: the verifier separates evidence from model prior.

## Limitations of the harness itself

- **Four cases is a seed, not a suite.** They are regression tests for known
  failures, so they are in-distribution by construction — the same criticism
  that applies to Tier 1.
- **Gold plans are authored by the system's developers.** Held-out cases
  should come from people who did not build it and who are not shown the tool
  list; that protocol is specified above and not yet executed.
- **The verifier is unvalidated.** Its agreement with human expert judgement
  has not been measured, so it moves the problem rather than closing it. It
  should be calibrated against expert labels before any headline number rests
  on it.
