# Skill: CORNETO network integration

The **integration layer** of the IGVF integrated data warehouse. Wraps
the **CORNETO** framework (Saez-Rodriguez lab,
[*Nat Mach Intell* 2025](https://www.nature.com/articles/s42256-025-01069-9),
GPL-3.0, https://github.com/saezlab/corneto) to infer context-specific
subnetworks by mixed-integer constrained optimization over a signed
directed prior-knowledge graph.

CORNETO reformulates many separate biological-network-inference tools
(CARNIVAL, COSMOS, PCSF/Steiner-tree, OmniPath subnetwork extraction,
FBA/iMAT, shortest-path) as instances of one MILP template:

  * Binary edge / vertex activation indicators + continuous flow
  * Flow conservation + sign consistency constraints
  * Objective = data-loss + L0 sparsity (edges, vertices) + L1 flow

## License

CORNETO is **GPL-3.0**. This skill installs CORNETO as a runtime
dependency and calls it via its public API — no source is copied into
the (Apache-2-licensed) IGVFagent repo. Install:

```
pip install corneto cvxpy 'pyscipopt>=5.0'
# optional faster solver (academic-free):
pip install gurobipy
```

## Subcommands

### demo
```
igvfagent corneto demo --beta 0.05 --solver SCIP
```
Self-test on a synthetic signed cascade (EGFR → SOS1 → RAS → RAF →
MEK → ERK → MYC) — proves CARNIVAL is wired end-to-end and writes
the selected subnetwork back into the warehouse with
`upstream='corneto:demo'`.

### pkn-from-kg
```
igvfagent corneto pkn-from-kg --label reactome --taxon 9606
```
Materialise a signed PKN SIF from `Data/Proteomics/KG/proteomics.sqlite`.
The PPI sources we have today (BioGRID, IntAct, HuRI) record physical
interaction without functional sign, so the default sign is **+1**
(positive regulation). Reactome edges are direction-aware and could
be re-signed via a follow-up; for now they're +1 too.

### carnival
```
# Build a PKN on the fly from the proteomics KG, run CARNIVAL with
# perturbed genes (signed) and DEGs (signed log2FC):
igvfagent corneto carnival \
    --perturbations Data/Corneto/perts.csv \
    --measurements  Data/Corneto/degs.csv \
    --beta 0.2 --solver SCIP --label perturb_seq_demo

# Or pass a pre-built PKN:
igvfagent corneto carnival --pkn pkn.sif \
    --perturbations perts.csv --measurements degs.csv
```
**Inputs** (CSV, headers required):

  * `perts.csv`: `gene,sign` with `sign ∈ {-1, +1}` — +1 for activation,
    −1 for knock-down / knockout.
  * `degs.csv`: `gene,score` with signed log2 fold-change.

Output: `selected_subnetwork.sif` (the minimum-cost upstream subnetwork
that explains the perturbations → measurements) + a `summary.json` +
a markdown report. The selected edges are also appended to the
warehouse's `edges` table with `upstream='corneto:carnival:<label>'`,
so they're queryable like any other source.

### steiner
```
igvfagent corneto steiner \
    --pkn ppi.sif --terminals prizes.csv \
    --root TP53 --solver SCIP
```
Prize-collecting Steiner tree. Good fit for **VAMP-seq variant prizes
on a PPI**: each gene's "prize" is its summed absolute abundance change
across all variants, and the optimal tree is the most informative
sub-network connecting the high-prize genes.

### multi-carnival (planned)
Joint inference across many conditions with shared sparsity (one
subnetwork per Perturb-seq guide, all sharing edge / vertex L0
penalties). Available via `corneto.methods.future.carnival.multi_carnival`
upstream; not yet wired here.

### write-playbook
```
igvfagent corneto write-playbook
```

## Cross-skill chaining

- `perturb-catalog pipeline` → `corneto carnival` — turn a Perturb-seq
  dataset's DEGs into an upstream signaling subnetwork.
- `proteomics build-kg` → `corneto pkn-from-kg` → `corneto carnival` —
  the PPI graph from BioGRID + Reactome becomes the prior for context
  inference.
- `rnaseq deg` → `corneto carnival` — DEGs from a bulk comparison
  drive the gap-to-data loss in CARNIVAL.
- `corneto carnival` → `warehouse query "SELECT * FROM edges WHERE
  upstream LIKE 'corneto:%'"` — every inferred subnetwork is a
  first-class citizen of the central warehouse, ready for downstream
  embedding extraction + foundation-model training.

## Why this design

- CORNETO is **the right primitive** for the integration layer because
  it unifies the math behind dozens of separate biological-network
  tools, with MILP optimality guarantees and multi-sample joint
  inference.
- We don't reinvent it — we install it as a runtime dependency,
  honour its GPL by not copying source, and use the public API.
- Selected subnetworks flow back into the **DuckDB warehouse**
  alongside every other producer, so the downstream embedding /
  foundation-model layer treats CORNETO-inferred edges as just
  another evidence stream.

## References

- Rodriguez-Mier P. et al. *CORNETO: a unified framework for the
  inference of context-specific networks.* **Nat Mach Intell** 2025.
  DOI: [10.1038/s42256-025-01069-9](https://www.nature.com/articles/s42256-025-01069-9).
- Source: [saezlab/corneto](https://github.com/saezlab/corneto)
- Manuscript experiments:
  [saezlab/corneto-manuscript-experiments](https://github.com/saezlab/corneto-manuscript-experiments)
