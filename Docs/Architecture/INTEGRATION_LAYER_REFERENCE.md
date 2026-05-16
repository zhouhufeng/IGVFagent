# Integration-layer reference — CORNETO et al.

This is the **algorithmic reference** for IGVFagent's network-integration
layer (`Scripts/network_integration_skill.py`). It distills the math
behind the **CORNETO** framework (Rodriguez-Mier et al., *Nat Mach
Intell* 2025; https://github.com/saezlab/corneto) and the published
methods CORNETO itself unifies — CARNIVAL, COSMOS, PCSF /
prize-collecting Steiner, OmniPath subnetwork extraction, FBA / iMAT,
shortest-path. IGVFagent's own implementation is a **clean-room
re-implementation in cvxpy**: same math, original Apache-2 source, no
runtime dependency on CORNETO (which is GPL-3.0). Math is not
copyrightable; this document records the reference algorithms and our
implementation decisions so future authors can extend the layer safely.

This doc is the result of a careful study of the CORNETO source tree
(`saezlab/corneto`) and the companion manuscript-experiments repo.

---

## 1. Problem statement

CORNETO frames context-specific biological network inference as **one
mathematical primitive: a mixed-integer linear program (MILP) over a
directed signed prior-knowledge graph (PKN)**. The decision variables
select an interpretable subnetwork constrained by:

- **Network-flow conservation** at interior vertices.
- **Sign consistency** — an edge can propagate activation only if its
  upstream vertex's sign matches the edge's sign.
- **Structured-sparsity penalties** — L0 on edges, L0 on vertices,
  L1 on flow.
- **Gap-to-data loss** — signed Hamming or L1 distance between the
  inferred vertex signs at "measurement" nodes and the observed
  scores.

Previously-published methods that map to the same template:

| Method | Inputs | Output | Notes |
|---|---|---|---|
| **CARNIVAL** | TF activities (signed); signed PKN (e.g. OmniPath) | Upstream subnetwork that explains TF activities | Headline method. |
| **COSMOS** | Multi-omic: metabolic + transcriptomic + signaling | Multi-layer subnetwork | Heterogeneous PKN. |
| **PCSF / Steiner** | Vertex prizes (e.g. abundance change); edge costs | Min-cost connecting subnetwork | Undirected variant common. |
| **CellNOpt-ILP** | Boolean phosphoproteomic; signed PKN | Boolean signaling model | Drops continuous flow. |
| **FBA / iMAT** | Genome-scale stoichiometry; gene expression bounds | Flux distribution | Pure LP. |
| **Shortest path / OmniPath** | Source / sink sets; weighted PKN | k shortest paths | Polynomial special case. |

The value-add over the patchwork of standalone Bioconductor / R
packages: **(a)** exact MILP optimality guarantees vs. heuristics,
**(b)** multi-sample joint inference (one subnetwork per condition
with shared sparsity), letting variant / perturbation / tissue
conditions borrow strength.

---

## 2. CARNIVAL — the headline algorithm

### Decision variables (per condition `c`)

| Symbol | Type | Per | Meaning |
|---|---|---|---|
| `F_e`         | continuous, [lb, ub] | edge `e` | network flow magnitude |
| `Fi_e`        | binary | edge | flow indicator (`F_e > 0`) |
| `N_act_v`     | binary | vertex | vertex activated (sign +1) |
| `N_inh_v`     | binary | vertex | vertex inhibited (sign −1) |
| `R_act_e`     | binary | edge   | edge sends activation downstream |
| `R_inh_e`     | binary | edge   | edge sends inhibition downstream |
| `L_v`         | continuous | vertex | DAG layer (MTZ) |

For a signed PKN edge `e = (u, sign_e, v)`:

- `sign_e = +1` → edge is activating
- `sign_e = -1` → edge is inhibiting

### Constraints

```
# (1) one sign per vertex / per edge
N_act_v + N_inh_v   ≤ 1                            ∀v
R_act_e + R_inh_e   ≤ 1                            ∀e

# (2) signal propagation ⇒ flow
R_act_e + R_inh_e   ≤ Fi_e                          ∀e

# (3) sign consistency at each edge e = (u, sign_e, v)
  if sign_e = +1:
    R_act_e ≤ N_act_u
    R_inh_e ≤ N_inh_u
  if sign_e = −1:
    R_act_e ≤ N_inh_u
    R_inh_e ≤ N_act_u

# (4) downstream activation needs ≥1 incoming activating edge
N_act_v ≤ Σ_{e ∈ in(v)} R_act_e                    ∀ non-boundary v
N_inh_v ≤ Σ_{e ∈ in(v)} R_inh_e                    ∀ non-boundary v

# (5) boundary (perturbations): force sign
if pert_v = +1:  N_act_v = 1
if pert_v = −1:  N_inh_v = 1

# (6) acyclicity — Miller–Tucker–Zemlin (optional, makes the
#     selected subnetwork a DAG; rooted at the boundary `_s`)
L_v − L_u ≥ Fi_e − (|V| − 1) · (1 − Fi_e)         ∀ edge e=(u, _, v)
1 ≤ L_v ≤ |V|                                      ∀v

# (7) network-flow conservation at interior vertices
A · F = 0     (A = vertex-incidence matrix; rows sum to 0 in interior)
F_e ≥ ε · (R_act_e + R_inh_e)                       (signal forces flow)
F_e ≤ M · Fi_e                                       (no flow when no signal)
```

Constraint (7) is the "Flow" primitive in CORNETO. The simpler "no
flow conservation, just cascade rule (4)" relaxation produces sensible
solutions on small graphs and is what IGVFagent's first
implementation uses; the full flow form is a follow-up.

### Objective

```
minimize  Σ_v∈measurements |score_v| · | (N_act_v − N_inh_v) − sign(score_v) |
        + β · Σ_e Fi_e                  # L0 on selected edges
        + λ_V · Σ_v (N_act_v + N_inh_v)  # L0 on selected vertices
        + λ_F · Σ_e F_e                  # L1 on flow
```

The gap-to-data is **signed**: a measurement with positive log2FC
contributes loss when the inferred vertex isn't activated; a
measurement with negative log2FC contributes loss when the vertex
isn't inhibited.

### Multi-sample variant

For conditions `c = 1 … C`, every per-vertex / per-edge variable is
duplicated per condition (`N_act_{v,c}`, `R_act_{e,c}`, …) and the
loss is summed over conditions. The L0 penalties on edges and
vertices share a single L0 term Σ_e Fi_e where `Fi_e = max_c
(R_act_{e,c} + R_inh_{e,c})`, encouraging a sparse *backbone* that
covers all conditions. This is the "shared sparsity" that lets a
Perturb-seq run with N guides recover one upstream signaling core.

### Solver

Any MILP-capable CVXPY backend works. CORNETO recommends Gurobi
(academic-free); IGVFagent defaults to **SCIP** (BSD-style, free) for
reproducibility on a laptop. HiGHS is also free but had value-readback
issues with the CORNETO 1.0.0a backend in our testing.

---

## 3. Prize-collecting Steiner tree

For when sign / direction isn't available and you only have **per-gene
prizes** (e.g. summed |abundance change| over all VAMP-seq variants in
a gene; GWAS hit p-value).

### Variables

| Symbol | Type | Per | Meaning |
|---|---|---|---|
| `y_v` | binary | vertex | vertex selected |
| `x_e` | binary | edge   | edge selected |

### Constraints

```
# (1) edges select only when both endpoints are selected
x_e ≤ y_u                                  e = (u, v)
x_e ≤ y_v

# (2) connectivity (Steiner tree): there must be a tree spanning the
#     selected vertices. Encode via a single-commodity flow from a
#     pseudo-root r:
A · F = b      where b_r = (#selected − 1), b_v = −1 if y_v=1, 0 else
F_e ≥ 0
F_e ≤ M · x_e
```

### Objective

```
maximize  Σ_v prize_v · y_v   −   Σ_e cost_e · x_e
```

Cost defaults to 1 per edge (so the objective is `Σ y_v · prize_v − |E_selected|`).

---

## 4. CORNETO's source layout — what we drew from

| Path in `saezlab/corneto` | What we learned |
|---|---|
| `corneto/methods/carnival.py` | `runVanillaCarnival`, `heuristic_carnival`, `get_selected_edges` |
| `corneto/methods/signaling.py` | `create_flow_graph`, `signflow_constraints`, `signflow`, `default_sign_loss` |
| `corneto/methods/steiner.py` | `exact_steiner_tree`, `create_exact_multi_steiner_tree` |
| `corneto/methods/future/carnival.py` | newer `CarnivalFlow`, `CarnivalILP`, `multi_carnival` |
| `corneto/methods/metabolism/fba.py` | `FBAProblem` (LP for flux balance) |
| `corneto/backend/_base.py` | `Backend.Flow`, `AcyclicFlow`, `Indicator`, `NonZeroIndicator`, `Problem`, `Variable` |
| `corneto/graph/_graph.py` | `Graph` data structure |
| `corneto/_constants.py` | `Solver`, `VarType`, `Direction`, `VAR_FLOW` |

CORNETO supports two backends (CVXPY and PICOS); IGVFagent uses only
**CVXPY** because it's the simpler integration and the larger ecosystem
of MILP solvers exposes it.

---

## 5. Manuscript experiments — case studies CORNETO ships

Six end-to-end notebooks in `saezlab/corneto-manuscript-experiments`:

1. **FBA** — exact LP flux balance on COBRA models.
2. **iMAT** — gene-expression-constrained metabolic-flux MILP.
3. **CARNIVAL (signed signaling)** — synthetic-PKN benchmark on
   CollecTRI; sweeps over λ and graph size; single- vs multi-sample.
4. **CARNIVAL on CPTAC LUAD** — real Omnipath PKN + decoupler-inferred
   TF activities (CollecTRI) + kinase activities (phosformer) across
   many patients via `multi_carnival(G, patient_data, lambd=int_lam)`.
5. **SteinerTrees** — `exact_steiner_tree` vs NetworkX and PCSF-fast.
6. **Showcase** — drug-response PCST on CPTAC + Fragpipe; cytokine
   single-cell Steiner.

All use the same `load → preprocess → assemble flow graph → solve →
visualise` pipeline. This is the template IGVFagent's `network`
skill follows.

---

## 6. Mapping to IGVFagent data sources

| CARNIVAL input | IGVFagent source |
|---|---|
| Signed directed PKN | `Data/Proteomics/KG/proteomics.sqlite` (BioGRID + IntAct + HuRI + Reactome interactions; signs default to +1 because most PPI sources don't annotate functional sign — Reactome direction can be re-signed in a follow-up pass). |
| Perturbations (signed) | CRISPRi / Perturb-seq guide-target gene with sign = −1 (knock-down). Source: `perturb-catalog dataset-rows`. |
| Measurements (signed) | DEGs (signed log2FC) from `rnaseq deg`, or per-target effect rows from the Perturbation Catalogue. |
| Multi-sample sharing | One condition per perturbed gene in a Perturb-seq run. |

| Steiner input | IGVFagent source |
|---|---|
| Per-gene prizes | Summed |abundance change| from `proteomics vampseq-analyze` outputs, GWAS catalog p-values, or DEG significance scores. |
| Edge costs | Default 1; future variant uses |1 − confidence_score| from the interactions table. |
| PPI graph | Same proteomics KG as CARNIVAL. |

The **subnetwork outputs flow back into the central DuckDB warehouse**
as `edges` rows tagged with `upstream='network:carnival:<label>'` /
`upstream='network:steiner:<label>'`. The Phase-4 embedding alignment
treats these inferred edges as positive contrastive pairs in the
condition-specific latent space — the same edge between two proteins
that would be unsupervised noise becomes a supervised positive once
the integration layer has selected it.

---

## 7. Why we re-implemented rather than depended on CORNETO

- **License**. CORNETO is GPL-3.0; IGVFagent is Apache-2.0. Static
  linking / source copy of GPL code into Apache code is incompatible.
  Runtime dependency works but adds a license surface to anyone who
  pip-installs IGVFagent.
- **Stability**. CORNETO 1.0.0a0 (current as of 2025) has a default
  PICOS backend with value-readback bugs on small MILPs. We work
  around with CVXPY-only paths.
- **Self-contained**. The MILP for CARNIVAL is ~150 lines of cvxpy.
  Re-implementing keeps every line of the integration layer under our
  own audit and Apache license.
- **Educational**. The reformulation lives in one file, side-by-side
  with the application code that feeds and consumes it.

**Math is not copyrightable** — algorithms can be re-implemented
freely. **Source code is.** We do not copy CORNETO source. We cite
CORNETO as the reference for the formulation and use it as a
reference implementation for behaviour comparison.

---

## 8. References

- Rodriguez-Mier P. et al. *CORNETO: a unified framework for the
  inference of context-specific networks.* **Nat Mach Intell** 2025.
  DOI: [10.1038/s42256-025-01069-9](https://www.nature.com/articles/s42256-025-01069-9).
- Liu, A. et al. *CARNIVAL: from inputs (TFs, perturbations) to
  pathways: inferring upstream signaling networks.* **npj Syst Biol
  Appl** 2019.
- Dugourd, A. et al. *COSMOS: causal mechanism of signal transduction
  and metabolism.* **Mol Syst Biol** 2021.
- Akhmedov, M. et al. *PCSF: prize-collecting Steiner forest from
  interactome data.* **Bioinformatics** 2017.
- Türei, D. et al. *OmniPath: guidelines and gateway for
  literature-curated signaling pathway resources.* **Nat Methods** 2016.

CORNETO source: [saezlab/corneto](https://github.com/saezlab/corneto)
· manuscript experiments: [saezlab/corneto-manuscript-experiments](https://github.com/saezlab/corneto-manuscript-experiments)
