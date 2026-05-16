# Skill: network-integration — context-specific subnetwork MILP

The **integration layer** of the IGVF data warehouse. Implements
**CARNIVAL** (signed signaling subnetwork) and **prize-collecting
Steiner tree** as pure-cvxpy mixed-integer programs.

The math is re-implemented from scratch in this module — no upstream
library import. Algorithm reference:
``Docs/Architecture/INTEGRATION_LAYER_REFERENCE.md`` (which documents
the CORNETO formulation we re-implement, with attribution).

License posture:
- This module: **Apache-2.0** (matches IGVFagent).
- CORNETO (the reference framework): GPL-3.0. **Not imported, not
  vendored.** The math is reused under standard practice — algorithms
  aren't copyrightable; source code is.

## Subcommands

### demo
```
igvfagent network demo --beta 0.05 --solver SCIP
```
Self-test on a synthetic 6-edge cascade
(EGFR → SOS1 → RAS → RAF → MEK → ERK → MYC). CARNIVAL must recover
all six edges. The inferred subnetwork is also appended to the
warehouse's ``edges`` table with ``upstream='network:demo'``.

### pkn-from-kg
```
igvfagent network pkn-from-kg --label reactome_pkn --taxon 9606
```
Materialise a signed PKN SIF from
``Data/Proteomics/KG/proteomics.sqlite``. The PPI sources currently
record physical interaction without functional sign — default sign is
**+1**. Reactome direction can re-sign in a follow-up.

### carnival
```
igvfagent network carnival \
    --perturbations Data/Network/perts.csv \
    --measurements  Data/Network/degs.csv \
    --beta 0.2 --solver SCIP --label perturb_seq_demo
```
**Inputs (CSV, headers required)**:

  * ``perts.csv`` — ``gene,sign`` with ``sign ∈ {-1, +1}`` (−1 for
    knock-down / knockout).
  * ``degs.csv`` — ``gene,score`` with signed log2 fold-change.

**Outputs (under ``Docs/Network/<ts>_<label>/``)**:

  * ``selected_subnetwork.sif`` — the minimum-cost upstream subnetwork
    that explains perturbations → measurements.
  * ``summary.json`` — objective value, status, sizes.
  * ``report.md`` — first-50-edges table.

Selected edges are also appended to the central DuckDB warehouse's
``edges`` table with ``upstream='network:carnival:<label>'`` so they're
queryable like any other source.

### steiner
```
igvfagent network steiner \
    --terminals Data/Network/vamp_prizes.csv \
    --pkn-limit 5000 --edge-cost 1.0 --solver SCIP
```
Prize-collecting Steiner tree (undirected). Good fit for VAMP-seq
abundance prizes on a PPI: each gene's prize is its summed |abundance
change| across all variants, and the optimal tree is the most
informative connecting subnetwork.

### pkn-from-sif
```
igvfagent network pkn-from-sif --input pkn.sif
```
Sanity-check a SIF file: edge count + sign distribution.

### write-playbook
```
igvfagent network write-playbook
```

## Algorithm summary (full spec in the reference doc)

### CARNIVAL MILP

Decision vars (per vertex / per edge, binary):
- ``N_act, N_inh`` — vertex activated / inhibited (one sign at most).
- ``R_act, R_inh`` — edge sends activation / inhibition downstream.
- ``Fi`` — edge selected (== ``R_act`` OR ``R_inh``).

Key constraints:
- One sign per vertex / per edge.
- Edge active ⇔ at least one role assigned.
- **Sign consistency**: a +edge propagates the upstream sign; a −edge
  swaps it.
- **Cascade**: non-boundary vertex activation requires an incoming
  edge whose role matches. Boundary (perturbation) vertices are
  pinned by hard constraints.
- (Optional MTZ acyclicity not enforced in v1 — sensible default for
  small dense PKNs; enable for large sparse ones to avoid feedback
  loops.)

Objective:
```
minimize  Σ_v |score_v| · |(N_act_v − N_inh_v) − sign(score_v)|
        + β · Σ_e Fi_e
        + λ_V · Σ_v (N_act_v + N_inh_v)
```

### Steiner MILP

Decision vars:
- ``y_v`` — vertex selected (binary).
- ``x_e`` — edge selected (binary, undirected pair).

Constraints:
- ``x_e ≤ y_u, y_v`` (edge requires both endpoints).
- (Connectivity enforced via simpler tree relaxation in v1 — full
  single-commodity flow encoding is documented in the reference doc
  as a follow-up.)

Objective:
```
maximize  Σ_v prize_v · y_v  −  cost · Σ_e x_e
```

## Cross-skill chaining

- ``perturb-catalog pipeline`` → ``network carnival`` — Perturb-seq
  DEGs become the CARNIVAL gap-to-data target.
- ``proteomics build-kg`` → ``network pkn-from-kg`` →
  ``network carnival`` — the PPI graph from BioGRID + Reactome is
  the prior.
- ``rnaseq deg`` → ``network carnival`` — bulk DEGs drive the
  loss term.
- ``proteomics vampseq-analyze`` → ``network steiner`` —
  per-gene |abundance change| sums become prizes.
- ``network *`` → ``warehouse query
  "SELECT * FROM edges WHERE upstream LIKE 'network:%'"`` — every
  inferred subnetwork lives in the central warehouse alongside the
  upstream sources.

## Dependencies

  * Required: ``cvxpy>=1.5``, ``numpy``, ``pyscipopt>=5.0`` (free
    MILP).
  * Optional: ``gurobipy`` (academic-free, faster).
  * Optional warehouse write-back: ``duckdb>=0.10``.

All Apache-2-compatible. No GPL dependencies.

## References

- **CORNETO** — Rodriguez-Mier P. et al., *Nat Mach Intell* 2025.
  DOI: [10.1038/s42256-025-01069-9](https://www.nature.com/articles/s42256-025-01069-9).
  Source: [saezlab/corneto](https://github.com/saezlab/corneto)
  (GPL-3.0; reference only, not imported).
- **CARNIVAL** — Liu, A. et al. *npj Syst Biol Appl* 2019.
- Full algorithmic reference + extended bibliography:
  [`Docs/Architecture/INTEGRATION_LAYER_REFERENCE.md`](../Architecture/INTEGRATION_LAYER_REFERENCE.md).
