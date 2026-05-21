# Skill: Network Integration Visualization

Publication-grade visualization for `network` skill outputs (CARNIVAL /
Steiner / external SIFs).

## Inputs

- `--sif <path>` — signed-SIF triples `src \t sign \t dst`. Required.
- `--prizes <csv>` (optional) — per-node `prize / sign / role` table.
  Columns recognized (case-insensitive): `node|gene|symbol`, `prize|weight|score`,
  `sign|direction`, `role` (perturbation / measurement / inferred).
- `--pathways <csv>` (optional) — `node, pathway` long-form table used
  for the pathway-enrichment bar panel.
- `--perturbation NAME` / `--measurement NAME` — tag specific nodes
  (repeatable) without needing a CSV.

## Outputs

All under `Docs/Network/<ts>_<label>_viz/`:

- `Plots/graph_force.png/.svg` — force-directed signed graph
- `Plots/pathway_enrichment.png/.svg` — top-K enriched pathways
- `Plots/degree_histogram.png/.svg` — degree distribution
- `Plots/edge_breakdown.png/.svg` — edge sign counts
- `Plots/composite_publication_figure.png/.svg` — 2×2 publication mosaic
- `graph_interactive.html` — interactive vis.js HTML (with `--html`)
- `node_summary.csv` — per-node degree, sign, prize, role, pathways
- `viz_report.md` — narrative report linking to all artefacts

## Example

```bash
# Visualize the EGFR -> MYC demo cascade
igvfagent network viz \
    --sif Docs/Network/<ts>_demo/selected_subnetwork.sif \
    --perturbation EGFR --measurement MYC \
    --layout spring --html --label demo
```

License: Apache-2.0. Deps: networkx (BSD), matplotlib (PSF), pyvis (MIT), pandas (BSD-3).
