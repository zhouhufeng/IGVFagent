"""Atlas I/O for the Tabula Sapiens skill.

Loading and aggregating the Tabula Sapiens 2.0 deposits without ever
holding 1.1 million cells x 62,000 genes in memory at once.

The published atlas ships as 28 per-tissue ``.h5ad`` files (57 GB
total), each carrying:

  layers  ``raw_counts``, ``decontXcounts``, ``log_normalized``,
          ``scale_data``
  obsm    ``X_pca``, ``X_scvi``, ``X_umap``, per-tissue variants
  obs     donor, tissue, method (10X / smartseq), assay,
          ``cell_ontology_class`` (182 fine types), ``broad_cell_class``,
          compartment, free_annotation, age, sex, ethnicity

That the *upstream* DecontX and scVI outputs are stored alongside the
raw counts is what lets this skill validate its own clean-room
reimplementations against the published ones rather than merely running.

The central routine here is :func:`mean_by_group`, which streams every
tissue file and accumulates a gene x cell-type mean of the
log-normalised layer. That matrix is the single input to the tau
statistic, and therefore to Figures 2 and 3.

License: Apache-2.0.
"""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Any, Iterator, Optional

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
ATLAS_DIR = ROOT / "Data" / "SingleCell" / "tabular-sapiens"
FIGSHARE_DIR = ATLAS_DIR / "figshare_v2"
GEO_DIR = ATLAS_DIR / "geo"

# The paper's own reference numbers, used to sanity-check a local copy.
PAPER = {
    "cells_total": 1_136_218,
    "cells_droplet": 1_093_048,
    "cells_facs": 43_170,
    "donors": 24,
    "tissues": 28,
    "fine_cell_types": 182,
    "fine_cell_types_droplet": 175,
    "broad_cell_types_droplet": 38,
    "populations": 701,
    "tf_total": 1639,
    "tf_in_gtf": 1637,
    "tau_threshold": 0.85,
    "tf_specific": 890,
    "tf_nonspecific": 745,
    "senescent_cells": 48_114,
    "senescent_frac": 0.044,
    "sag_genes": 3_792,
    "cnmf_k": 34,
    "senescence_pathways": 17,
}

_TISSUE_RE = re.compile(r"^([A-Za-z_]+?)_TSP")


def tissue_of(path: "str | Path") -> str:
    """Tissue name from a figshare filename, e.g. ``Large_Intestine``."""
    m = _TISSUE_RE.match(Path(path).name)
    return m.group(1) if m else Path(path).stem


def atlas_files(directory: "Optional[str | Path]" = None,
                *, tissues: "Optional[list[str]]" = None) -> "list[Path]":
    """Per-tissue ``.h5ad`` files, optionally filtered by tissue name."""
    d = Path(directory) if directory else FIGSHARE_DIR
    if d.is_file():
        return [d]
    if not d.is_dir():
        raise SystemExit(
            f"No atlas directory at {d}.\n"
            f"Fetch it first:  igvfagent tabula pull-figshare")
    files = sorted(d.glob("*.h5ad"))
    if not files:
        raise SystemExit(f"No .h5ad files under {d}")

    # Skip partially-downloaded files. A 57 GB fetch is usually still
    # running while analysis starts, and a truncated h5ad fails deep
    # inside HDF5 with "truncated file: eof = ..." rather than anywhere
    # informative. Size-check against the download manifest when present.
    manifest = d.parent / "_manifests" / "figshare_v2_files.tsv"
    if manifest.is_file():
        expected: "dict[str, int]" = {}
        for line in manifest.read_text().splitlines():
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2].isdigit():
                expected[parts[1]] = int(parts[2])
        complete, partial = [], []
        for f in files:
            exp = expected.get(f.name)
            (partial if (exp is not None and f.stat().st_size != exp)
             else complete).append(f)
        if partial:
            logging.warning(
                "skipping %d incomplete download(s): %s",
                len(partial), ", ".join(f.name.split("_TSP")[0] for f in partial))
        files = complete
        if not files:
            raise SystemExit(
                f"Every .h5ad under {d} is incomplete. Let "
                f"`igvfagent tabula pull-figshare` finish first.")
    if tissues:
        want = {t.strip().lower() for t in tissues if t.strip()}
        files = [f for f in files if tissue_of(f).lower() in want]
        if not files:
            raise SystemExit(
                f"No tissue matched {sorted(want)}. Available: "
                f"{sorted({tissue_of(f) for f in sorted(d.glob('*.h5ad'))})}")
    return files


def cell_metadata(path: "Optional[str | Path]" = None) -> Any:
    """The full 1.1M-cell metadata table from GEO GSE306755.

    Far cheaper than opening 28 h5ads when only ``obs`` is needed --
    which is the case for the donor/tissue composition of Figure 1.
    """
    import pandas as pd
    if path:
        p = Path(path)
    else:
        cands = sorted(GEO_DIR.glob("*metadata*.csv.gz")) + \
                sorted(GEO_DIR.glob("*metadata*.csv"))
        if not cands:
            raise SystemExit(
                f"No cell metadata under {GEO_DIR}.\n"
                f"Fetch it:  igvfagent tabula pull-geo")
        p = cands[0]
    df = pd.read_csv(p, low_memory=False)
    if "Unnamed: 0" in df.columns:
        df = df.rename(columns={"Unnamed: 0": "cell_id"})
    return df


def open_atlas(path: "str | Path", *, backed: bool = True) -> Any:
    """Open one tissue h5ad, backed by default so obs is cheap."""
    import anndata as ad
    return ad.read_h5ad(str(path), backed="r" if backed else None)


def read_layer_rows(path: "str | Path", layer: str, start: int, stop: int,
                    *, gene_idx: Any = None) -> Any:
    """Dense row block of an on-disk CSR layer, without loading the layer.

    anndata's ``backed='r'`` only keeps ``.X`` on disk -- touching
    ``.layers[...]`` materialises the WHOLE layer in RAM. For a 45k-cell
    tissue that is gigabytes, and streaming 28 of them OOMs a 31 GB box
    (it did). The h5ad layer is a plain CSR triple, so slicing rows
    directly off HDF5 costs only the nonzeros in the requested block.

    Parameters
    ----------
    path, layer
        h5ad file and layer name.
    start, stop
        Half-open row range.
    gene_idx
        Optional column subset, applied while still sparse.

    Returns
    -------
    Dense ``(stop-start, n_selected_genes)`` float64 array.
    """
    import h5py
    import numpy as np

    with h5py.File(str(path), "r") as f:
        g = f["layers"][layer]
        n_genes = int(g.attrs["shape"][1])
        indptr = g["indptr"][start:stop + 1]
        lo, hi = int(indptr[0]), int(indptr[-1])
        if hi <= lo:
            cols = n_genes if gene_idx is None else len(gene_idx)
            return np.zeros((stop - start, cols), dtype=np.float64)
        data = g["data"][lo:hi].astype(np.float64)
        indices = g["indices"][lo:hi].astype(np.int64)

    rel = indptr - indptr[0]
    n_rows = stop - start
    if gene_idx is None:
        out = np.zeros((n_rows, n_genes))
        for r in range(n_rows):
            a, b = int(rel[r]), int(rel[r + 1])
            out[r, indices[a:b]] = data[a:b]
        return out

    # Column subset: map selected gene columns to 0..k-1 and drop the rest
    # before densifying, which is the whole point of doing this by hand.
    gene_idx = np.asarray(gene_idx, dtype=np.int64)
    lookup = np.full(n_genes, -1, dtype=np.int64)
    lookup[gene_idx] = np.arange(gene_idx.size)
    out = np.zeros((n_rows, gene_idx.size))
    for r in range(n_rows):
        a, b = int(rel[r]), int(rel[r + 1])
        if b <= a:
            continue
        cols = lookup[indices[a:b]]
        keep = cols >= 0
        if keep.any():
            out[r, cols[keep]] = data[a:b][keep]
    return out


def iter_tissues(files: "list[Path]", *, backed: bool = True
                 ) -> "Iterator[tuple[str, Any]]":
    """Yield ``(tissue, AnnData)`` for each file, closing as it goes."""
    for f in files:
        a = open_atlas(f, backed=backed)
        try:
            yield tissue_of(f), a
        finally:
            if getattr(a, "isbacked", False) and a.file is not None:
                a.file.close()


def mean_by_group(
    files: "list[Path]",
    *,
    group_key: str = "cell_ontology_class",
    layer: str = "log_normalized",
    method: "Optional[str]" = "10X",
    genes: "Optional[list[str]]" = None,
    gene_key: str = "gene_symbol",
    chunk: int = 20_000,
) -> "tuple[Any, list[str], list[str]]":
    """Mean expression per gene per group, streamed across tissue files.

    This is the matrix Figures 2 and 3 are computed from. The paper takes
    "the mean of the log-normalized expression of each TF across the 175
    cell types detected in the droplet microfluidic subset", so the
    defaults are ``log_normalized``, ``cell_ontology_class`` and
    ``method='10X'``.

    Accumulates a running sum and count per group rather than
    concatenating: the full atlas is 1.1M x 62k, which does not fit in
    memory, but a genes x cell-types matrix is trivially small.

    Parameters
    ----------
    files
        Per-tissue h5ad paths.
    group_key
        ``obs`` column to average within (fine or broad cell type).
    layer
        Layer to average. ``None`` uses ``.X``.
    method
        Restrict to one assay technology; ``None`` keeps all.
    genes
        Restrict to these gene symbols. Hugely faster for a TF-only run.
    gene_key
        ``var`` column holding the symbol.
    chunk
        Cells read per block from the backed file.

    Returns
    -------
    ``(means, gene_names, group_names)`` with ``means`` shaped
    ``(n_genes, n_groups)``; groups absent everywhere are dropped.
    """
    import numpy as np
    import pandas as pd

    want = {g.upper() for g in genes} if genes else None
    sums: "dict[str, Any]" = {}
    counts: "dict[str, int]" = {}
    gene_names: "Optional[list[str]]" = None
    gene_idx: "Optional[Any]" = None

    for path in files:
        a = open_atlas(path, backed=True)
        try:
            obs = a.obs
            if method and "method" in obs.columns:
                keep_mask = (obs["method"].astype(str) == method).to_numpy()
            else:
                keep_mask = np.ones(a.n_obs, dtype=bool)
            if group_key not in obs.columns:
                logging.warning("%s has no %s; skipped", path.name, group_key)
                continue
            if not keep_mask.any():
                continue

            if gene_names is None:
                var = a.var
                symbols = (var[gene_key].astype(str).to_numpy()
                           if gene_key in var.columns
                           else var.index.astype(str).to_numpy())
                if want is not None:
                    gene_idx = np.flatnonzero(
                        np.isin(np.char.upper(symbols.astype(str)),
                                list(want)))
                    if gene_idx.size == 0:
                        raise SystemExit(
                            "none of the requested genes are in the atlas")
                else:
                    gene_idx = np.arange(symbols.size)
                gene_names = [str(symbols[i]) for i in gene_idx]

            groups = obs[group_key].astype(str).to_numpy()
            # Contiguous slices, not fancy row indexing: a backed h5ad
            # reads a slice in one go but gathers arbitrary rows one at a
            # time. Subset genes on the SPARSE matrix before densifying --
            # densifying all ~62k genes to keep 1.6k of them was costing
            # a ~38x memory and time blow-up.
            for start in range(0, a.n_obs, chunk):
                stop = min(start + chunk, a.n_obs)
                sel = keep_mask[start:stop]
                if not sel.any():
                    continue
                if layer:
                    M = read_layer_rows(path, layer, start, stop,
                                        gene_idx=gene_idx)
                else:
                    M = a[start:stop].X
                    if hasattr(M, "toarray"):
                        M = M[:, gene_idx].toarray()
                    else:
                        M = np.asarray(M)[:, gene_idx]
                M = np.asarray(M, dtype=np.float64)[sel]
                g = groups[start:stop][sel]
                # One pass with integer codes beats a boolean mask per
                # group when a chunk spans many cell types.
                labs, codes = np.unique(g, return_inverse=True)
                acc = np.zeros((labs.size, M.shape[1]))
                np.add.at(acc, codes, M)
                n_per = np.bincount(codes, minlength=labs.size)
                for j, lab in enumerate(labs):
                    if lab in sums:
                        sums[lab] += acc[j]
                        counts[lab] += int(n_per[j])
                    else:
                        sums[lab] = acc[j].copy()
                        counts[lab] = int(n_per[j])
            logging.info("  %-22s %6d cells", tissue_of(path), int(keep_mask.sum()))
        finally:
            if getattr(a, "isbacked", False) and a.file is not None:
                a.file.close()

    if not sums:
        raise SystemExit("no cells matched; check --method and --group-key")
    group_names = sorted(sums)
    means = np.zeros((len(gene_names or []), len(group_names)))
    for j, lab in enumerate(group_names):
        means[:, j] = sums[lab] / max(counts[lab], 1)
    return means, (gene_names or []), group_names


def pseudobulk(
    files: "list[Path]",
    *,
    group_keys: "list[str]",
    layer: str = "raw_counts",
    method: "Optional[str]" = None,
    chunk: int = 20_000,
) -> "tuple[Any, list[str], Any]":
    """Summed counts per combination of ``group_keys``.

    Cells within a donor are not independent replicates, so any
    differential test across donors or sexes must be run on pseudobulks
    rather than on cells -- treating cells as replicates is what
    produces implausibly small p-values.

    Returns ``(counts, gene_names, obs)`` where ``counts`` is
    ``(n_groups, n_genes)`` and ``obs`` is a DataFrame of the group keys.
    """
    import numpy as np
    import pandas as pd

    acc: "dict[tuple, Any]" = {}
    ncell: "dict[tuple, int]" = {}
    gene_names: "Optional[list[str]]" = None

    for path in files:
        a = open_atlas(path, backed=True)
        try:
            obs = a.obs
            missing = [k for k in group_keys if k not in obs.columns]
            if missing:
                logging.warning("%s lacks %s; skipped", path.name, missing)
                continue
            mask = (np.ones(a.n_obs, dtype=bool) if not method
                    or "method" not in obs.columns
                    else (obs["method"].astype(str) == method).to_numpy())
            if not mask.any():
                continue
            if gene_names is None:
                var = a.var
                gene_names = list(
                    var["gene_symbol"].astype(str) if "gene_symbol" in var.columns
                    else var.index.astype(str))
            keys = list(zip(*[obs[k].astype(str).to_numpy() for k in group_keys]))
            idx = np.flatnonzero(mask)
            for start in range(0, idx.size, chunk):
                rows = idx[start:start + chunk]
                block = a[rows]
                M = block.layers[layer] if layer else block.X
                if hasattr(M, "toarray"):
                    M = M.toarray()
                M = np.asarray(M, dtype=np.float64)
                sub = [keys[i] for i in rows]
                for kk in set(sub):
                    m = np.array([s == kk for s in sub])
                    s = M[m].sum(axis=0)
                    if kk in acc:
                        acc[kk] += s
                        ncell[kk] += int(m.sum())
                    else:
                        acc[kk] = s
                        ncell[kk] = int(m.sum())
        finally:
            if getattr(a, "isbacked", False) and a.file is not None:
                a.file.close()

    if not acc:
        raise SystemExit("no cells matched for pseudobulk")
    labels = sorted(acc)
    counts = np.vstack([acc[k] for k in labels])
    obs = pd.DataFrame(labels, columns=group_keys)
    obs["n_cells"] = [ncell[k] for k in labels]
    return counts, (gene_names or []), obs


def gene_flags(
    files: "list[Path]",
    genes: "list[str]",
    *,
    layer: str = "raw_counts",
    gene_key: str = "gene_symbol",
    chunk: int = 20_000,
) -> Any:
    """Per-cell boolean expression of named genes, plus obs annotation.

    The senescence analysis needs exactly this: which cells are
    ``CDKN2A+`` and ``MKI67-``. Streaming a handful of gene columns over
    1.1M cells is cheap; loading the matrix is not.

    Returns a DataFrame indexed by cell with one boolean column per gene
    plus donor / tissue / cell type / age / sex.
    """
    import numpy as np
    import pandas as pd

    want = [g.upper() for g in genes]
    keep_obs = ["donor", "tissue", "cell_ontology_class", "broad_cell_class",
                "compartment", "method", "age", "sex"]
    frames = []
    for path in files:
        a = open_atlas(path, backed=True)
        try:
            var = a.var
            symbols = (var[gene_key].astype(str).to_numpy()
                       if gene_key in var.columns
                       else var.index.astype(str).to_numpy())
            up = np.char.upper(symbols.astype(str))
            cols = {g: np.flatnonzero(up == g) for g in want}
            missing = [g for g, ix in cols.items() if ix.size == 0]
            if missing:
                logging.warning("%s: %s not found", path.name, missing)
            idx_all = np.concatenate([ix for ix in cols.values() if ix.size]) \
                if any(ix.size for ix in cols.values()) else np.array([], int)
            if idx_all.size == 0:
                continue

            out = {g: np.zeros(a.n_obs, dtype=bool) for g in want}
            # Same reason as mean_by_group: anndata materialises a whole
            # layer even in backed mode, so read the CSR rows off HDF5
            # and keep only the handful of gene columns we asked for.
            order = np.argsort(idx_all)
            sorted_idx = idx_all[order]
            pos = {g: np.searchsorted(sorted_idx, ix) for g, ix in cols.items()
                   if ix.size}
            for start in range(0, a.n_obs, chunk):
                stop = min(start + chunk, a.n_obs)
                if layer:
                    M = read_layer_rows(path, layer, start, stop,
                                        gene_idx=sorted_idx)
                else:
                    M = a[start:stop].X
                    if hasattr(M, "toarray"):
                        M = M[:, sorted_idx].toarray()
                    else:
                        M = np.asarray(M)[:, sorted_idx]
                M = np.asarray(M)
                for g, p_ in pos.items():
                    out[g][start:stop] = (M[:, p_] > 0).any(axis=1)

            df = pd.DataFrame(out, index=a.obs_names.astype(str))
            for c in keep_obs:
                if c in a.obs.columns:
                    # Categorical, not object: these columns hold a handful
                    # of repeated labels across 1.1M cells, and Python str
                    # objects for all of them cost gigabytes once 28 tissue
                    # frames are held at once.
                    df[c] = pd.Categorical(a.obs[c].astype(str).to_numpy())
            df["tissue_file"] = pd.Categorical([tissue_of(path)] * len(df))
            frames.append(df)
            logging.info("  %-22s %7d cells", tissue_of(path), a.n_obs)
        finally:
            if getattr(a, "isbacked", False) and a.file is not None:
                a.file.close()

    if not frames:
        raise SystemExit("no tissue yielded any of the requested genes")
    return pd.concat(frames, axis=0)
