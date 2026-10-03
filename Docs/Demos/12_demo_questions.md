# IGVFagent — 12-question live demo

**Total run time:** ~22 minutes if you talk through each step.
**Setup:** Streamlit UI at <http://127.0.0.1:8501> with Ollama + Gemma loaded.

The 12 questions are organized as 5 acts with a visual-climax cadence —
each act has at least one "big plot" moment. Question text is verbatim:
copy-paste it straight into the chat box.

Before starting, pre-warm one h5ad so Q6 is instant:
```bash
python3 -c "import scanpy as sc; sc.datasets.pbmc3k().write_h5ad('/tmp/pbmc3k.h5ad')"
```

---

## Act 1 · Orient (~3 min)

### Q1. Front-page IGVF stats
```
What's currently in the IGVF Portal? Give me a front-page summary —
total datasets, top assays, top biosamples, recent releases.
```
- **Tool the agent should pick:** `igvf_frontpage_summary` (or live Portal queries)
- **Watch for:** counts that match the live Portal homepage; agent narrates the
  modality mix (10x multiome, SHARE-seq, MPRA, Perturb-seq, etc.).
- **~5 s.**

### Q2. Knowledge Graph evidence pack for APOE
```
Tell me everything in the IGVF Knowledge Graph about APOE — variants,
transcripts, proteins, regulatory elements, diseases, pathways.
```
- **Tool:** `kg_gene` with `symbol="APOE"`.
- **Watch for:** the agent makes ONE tool call and gets a comprehensive
  multi-hop evidence pack. Open the generated `report.md` from the artefact
  pane — it has all five categories filled.
- **~10 s.**

---

## Act 2 · Visual climax: regulatory genomics (~6 min)

### Q3. rE2G enhancer-gene linkage browser for BRCA1
```
Show me the rE2G enhancer-gene linkage arcs around BRCA1 and render the
browser visualization with SCREEN cCREs underneath.
```
- **Tool:** `encode_browser` with `region="chr17:43000000-43200000"`,
  `with_ccre=True`, `re2g="auto"`, `re2g_gene_filter="BRCA1"`.
- **Watch for:** 150+ score-coloured Bézier arcs all converging on the
  BRCA1 TSS, SCREEN cCRE track below, dashed vertical TSS marker. The
  companion `.rE2G.csv` is the audit trail.
- **★★★ Big visual moment** — this capability shipped this week.
- **~10 s.**

### Q4. GM12878 super-enhancers → target genes
```
Discover GM12878 H3K27ac super-enhancers, rank them with a hockey-stick
plot, and link the top 10 to their target genes via Catalog rE2G,
loops, and proximity.
```
- **Tool:** `se_targets_pipeline` with `biosample="GM12878"`,
  `target="H3K27ac"`, `assembly="GRCh38"`, `top_n=10`, `label="demo"`.
- **Watch for:** ROSE-style call → ~1 SE at the MALAT1/NEAT1 locus →
  4-stream linkage (loops + Catalog rE2G + proximity + cCRE composition)
  → hockey-stick + top-SE browser PNGs.
- **~60 s** (downloads ENCODE peak BED on first run; subsequent runs are fast).

---

## Act 3 · Single-cell (~5 min)

### Q5. Discover macrophage SPLiT-seq datasets
```
Find IGVF Parse SPLiT-seq datasets profiling macrophages in mouse and
write the per-pool manifest.
```
- **Tool:** `splitseq_retrieve` (or `singlecell` discovery).
- **Watch for:** the agent narrows to mouse macrophage, returns a clean
  manifest. Good intro for Q6 — "we discovered the data, now let's analyse some."
- **~5 s.**

### Q6. PBMC3k UMAP pipeline (Scanpy)
```
Run a full single-cell analysis on /tmp/pbmc3k.h5ad — QC, normalize,
PCA, UMAP, t-SNE, Leiden clustering, and marker genes. Overlay CD3D,
CD8A, MS4A1, and LYZ on the UMAP.
```
- **Tool:** `sc_pipeline` with `input="/tmp/pbmc3k.h5ad"`,
  `highlight_genes="CD3D,CD8A,MS4A1,LYZ"`, `resolution=0.6`.
- **Watch for:** ~20-second pipeline → 6 Leiden clusters → 4-panel UMAP
  with CD3D (T-cells), CD8A (CD8+ subset), MS4A1 (B-cells), LYZ
  (monocytes) lighting up their canonical clusters.
- **★★★ Big visual moment.** After the answer renders, click to the
  🔬 **Single-cell** tab in the UI — the same `processed.h5ad` is
  auto-discovered; flip embedding to t-SNE, add `APOE`, `TREM2` to the
  gene overlay live.
- **~20 s.**

---

## Act 4 · Functional & perturbation (~5 min)

### Q7. PTEN VAMP-seq deep analysis
```
Pull the published PTEN VAMP-seq scoreset from MaveDB and produce the
canonical 6-panel deep analysis: distribution, residue × amino-acid
heatmap, per-residue mean with domain track, replicate concordance,
abundance classes, and the ranked variant curve.
```
- **Tool:** `proteomics_vampseq_analyze` with `gene="PTEN"`.
- **Watch for:** 4,112 missense variants → iconic residue × AA heatmap
  (phosphatase domain shows the classic vulnerability pattern; C-tail
  is tolerant) + per-residue mean with PIP4-bind / Phosphatase / C2 /
  C-tail track. Replicate Pearson r ≈ 0.65.
- **★★★ Big visual moment.**
- **~15 s.**

### Q8. BRCA1 in the Perturbation Catalogue
```
What perturbation data exists in the Perturbation Catalogue for BRCA1
— across MAVE, CRISPR screens, and Perturb-seq? Include the top GSEA
hallmarks.
```
- **Tool:** `perturb_catalog_pipeline` with `gene="BRCA1"`.
- **Watch for:** 1,191 CRISPR-screen + 7 Perturb-seq datasets; canonical
  HALLMARK_E2F_TARGETS / G2M_CHECKPOINT / MYC_TARGETS signature (textbook
  for BRCA1).
- **~5 s.**

### Q9. Variants → cCRE → logistic → Miami
```
Take the variants in
Data/Input/VariantList/example_variants.csv, score them against IGVF
Catalog evidence and ENCODE cCRE classes, fit a logistic model, and
produce the volcano and Miami plots with the full markdown report.
```
- **Tool:** `advanced_variant_pipeline` with
  `variants="Data/Input/VariantList/example_variants.csv"`.
- **Watch for:** volcano + Miami + evidence-overlap PNGs + research-grade
  markdown report. This is the integrative variant workflow that closes
  the loop from raw variant → mechanistic interpretation.
- **~30 s.**

---

## Act 5 · Cross-evidence integration & literature (~3 min)

### Q10. Validate three genes against the literature
```
Validate APOE, TREM2, and LDLR against published Alzheimer and
cardiovascular literature. Flag any genes where the IGVF Catalog
evidence disagrees with the prior literature.
```
- **Tool:** `ref_validate` with `terms=["APOE","TREM2","LDLR"]`,
  `context=["Alzheimer","cardiovascular"]`.
- **Watch for:** multi-source pull (PubMed + Semantic Scholar +
  OpenAlex), per-gene literature corroboration count, agent narrates
  agreement vs disagreement.
- **~30 s.**

### Q11. Explain a real IGVF dataset
```
Explain IGVFDS7013XXYV in plain language — what assay, what biosample,
what files are available, and what the report says.
```
- **Tool:** `explain` (data illustration / interpretation).
- **Watch for:** the agent reads the dataset record + linked report and
  narrates it in plain English. Proves it can do comprehension, not
  just retrieval.
- **~5 s.**

### Q12. Capstone: integrative proteomics workflow
```
Build a proteomics protein-protein interaction knowledge graph from
Reactome and IGVF, then show me TP53's interaction neighborhood and
survey recent Nature, Cell, and Science papers on VAMP-seq.
```
- **Tool chain:** `proteomics_download` (reactome+igvf) → `proteomics_build_kg`
  → `proteomics_kg_visualize` with `gene="TP53"` → `proteomics_assay_survey`.
- **Watch for:** the agent chains 3–4 tool calls itself (ReAct in action),
  produces a network plot of TP53's PPI neighbourhood, then a literature
  table for VAMP-seq family papers. This is the orchestration headliner.
- **~90 s** (Reactome is cached after the first download).

---

## Notes for the live audience

- **Every question** the agent answers is also a CLI command you could
  run from the terminal. Use the "view CLI" toggle in the artefact pane
  to show the underlying `igvfagent <skill>` invocation.
- **The Plan → Action → Results → Evaluation loop** is visible in the
  chat: each tool call shows up as a `🛠️` line with arguments + result
  summary before the agent's narrative continues.
- **Three tabs** in the UI to flip between:
  - 💬 **Chat** — all 12 questions land here
  - 🕸 **Knowledge Graph** — after Q12, the local proteomics KG is alive
  - 🔬 **Single-cell** — after Q6, the PBMC3k h5ad is auto-discovered

## Recovery scripts (if anything breaks mid-demo)

```bash
# Re-warm PBMC3k
python3 -c "import scanpy as sc; sc.datasets.pbmc3k().write_h5ad('/tmp/pbmc3k.h5ad')"

# Re-build proteomics KG (Reactome is small + fast)
igvfagent proteomics download --source reactome
igvfagent proteomics build-kg --sources reactome,igvf

# Reset UI cleanly
lsof -nP -iTCP:8501 -sTCP:LISTEN -t | xargs kill
igvfagent ui --host 127.0.0.1 --port 8501
```
