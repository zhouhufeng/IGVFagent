#!/usr/bin/env python3
"""Build a SHARE-seq skin RNA AnnData via the shared _scload loaders.

Dogfoods Scripts/_scload.py (internalized benchmark loaders): streams the dense
gene-by-cell RNA counts, attaches Ma 2020's author cell-type labels, and
subsamples to a memory-safe size for `share rna-qc` + clustering concordance.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Scripts"))
import _scload  # noqa: E402

D = ROOT / "Benchmarks/_data/ma2020_shareseq"
RNA = D / "GSM4156608_skin.late.anagen.rna.counts.txt.gz"
CT = D / "GSM4156597_skin_celltype.txt.gz"
N_SUB = 15000
SEED = 0

if not RNA.is_file():
    sys.exit(f"missing {RNA} — extract it from GSE140203_RAW.tar first (see run.sh)")

print("==> loading dense gene×cell RNA counts via _scload")
a = _scload.dense_gene_by_cell_tsv(RNA)
# RNA header barcodes use commas (R1.01,R2.01,…); celltype table uses dots.
a.obs_names = [b.replace(",", ".") for b in a.obs_names]
print(f"    RNA: {a.n_obs:,} cells × {a.n_vars:,} genes")

# attach Ma 2020 author cell types (join on rna.bc)
ct = pd.read_csv(CT, sep="\t")
ct["rna.bc"] = ct["rna.bc"].astype(str)
ct = ct.drop_duplicates("rna.bc").set_index("rna.bc")
a.obs["cell_type"] = ct["celltype"].reindex(a.obs_names).values
labeled = a[a.obs["cell_type"].notna()].copy()
print(f"    labeled cells (in Ma 2020 final set): {labeled.n_obs:,} "
      f"({labeled.obs['cell_type'].nunique()} cell types)")

sub = _scload.subsample_cells(labeled, N_SUB, seed=SEED)
sub.obs["sample_id"] = "skin_late_anagen"
out = D / "shareseq_skin_rna.h5ad"
sub.write_h5ad(out, compression="gzip")
print(f"==> wrote {out.name}: {sub.n_obs:,} cells × {sub.n_vars:,} genes")
