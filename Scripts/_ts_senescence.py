"""Senescence analysis for the Tabula Sapiens skill (Figure 4).

Reproduces the CDKN2A+ MKI67- workflow of Tabula Sapiens 2.0
(Tabula Sapiens Consortium, *Cell* 2026), clean-room from the published
Methods and the ``paper2/Figure4`` notebooks (BSD-3-Clause upstream; no
source copied).

The paper's logic, and why each step is shaped the way it is:

1. **Define the population.** Cells carrying ``CDKN2A`` transcripts and
   lacking ``MKI67`` (a proliferation marker) -- 48,114 cells, ~4.4% of
   the atlas, over 145 cell types in 25 tissues. This is an enrichment,
   not a pure sort: the paper is explicit that the population "may
   include some non-senescent cells and exclude senescent cells that
   exhibit alternative senescence programs".

2. **Balance the comparison.** Tissues where one sex has only a single
   donor are dropped first, so a sex effect cannot masquerade as a
   senescence effect. That removes 52,149 cells across six tissues,
   leaving 1.08M from 25 tissues and 21 donors.

3. **Find senescence-associated genes per stratum.** Differential
   expression of CDKN2A+ MKI67- against the rest, *within* each
   donor x tissue x broad-cell-type stratum. Stratifying is the point:
   a gene that merely marks a cell type or a donor would otherwise
   dominate.

4. **Keep what replicates.** A gene is a SAG only if it clears the
   thresholds in at least half the donors for some tissue-and-cell-type
   -- 3,792 genes.

5. **Ask whether the program is universal.** Count in how many cell
   types each SAG is up. The answer in the paper is that it is not:
   after CDKN2A, the most universal gene (CDKN2B) reaches only ~60% of
   cell types, and canonical SASP factors such as IL6 and IL1B do not
   even clear the bar.

License: Apache-2.0.
"""

from __future__ import annotations

import logging
from typing import Any, Optional

# Thresholds are the paper's, named rather than inlined so a run that
# deviates is visible in the output rather than buried in a call site.
MIN_LOG2FC = 0.5
MAX_PADJ = 0.001
MIN_PCT_IN_GROUP = 0.5
MIN_DONOR_FRACTION = 0.5
PREVALENCE_LOG2FC = 0.5

SENESCENCE_MARKER = "CDKN2A"
PROLIFERATION_MARKER = "MKI67"

# Canonical hallmark genes the paper checks the population against.
HALLMARKS = {
    "cell_cycle_arrest": ["CDKN2A", "CDKN1A", "CDKN2B", "TP53", "RB1"],
    "dna_damage": ["H2AX", "H2AFX", "TP53BP1", "ATM"],
    "sasp_regulator": ["NFKB1", "TGFB1", "CEBPB", "GATA4"],
    "sasp_factor": ["IL6", "IL1B", "CXCL8", "SERPINE1", "MIF",
                    "TIMP2", "HMGB1", "IGFBP7"],
    "anti_apoptotic": ["BCL2L1", "BIRC3", "BCL2"],
    "lysosomal": ["GLB1", "ATP6V0D1", "ATP6V1D", "ATP6V1A"],
}


def _np():
    import numpy as np
    return np


# ---------------------------------------------------------------------------
# Population definition
# ---------------------------------------------------------------------------

def balanced_strata(obs: Any, *, tissue_key: str = "tissue",
                    sex_key: str = "sex", donor_key: str = "donor") -> Any:
    """Boolean mask keeping only tissue x sex strata with >1 donor.

    The paper drops these before anything else: with a single donor per
    tissue-sex, donor identity and sex are perfectly confounded, and any
    apparent senescence difference could be either.
    """
    import pandas as pd
    df = obs
    counts = (df.drop_duplicates(subset=[donor_key, tissue_key, sex_key])
                .groupby([tissue_key, sex_key]).size())
    keep_pairs = {k for k, v in counts.items() if v > 1}
    pairs = list(zip(df[tissue_key].astype(str), df[sex_key].astype(str)))
    return pd.Series([p in keep_pairs for p in pairs], index=df.index)


def call_senescent(flags: Any, *,
                   marker: str = SENESCENCE_MARKER,
                   proliferation: str = PROLIFERATION_MARKER) -> Any:
    """``CDKN2A+ MKI67-`` boolean, from a per-cell gene-flag table."""
    for c in (marker, proliferation):
        if c not in flags.columns:
            raise SystemExit(f"gene flag table has no {c!r} column")
    return flags[marker].astype(bool) & (~flags[proliferation].astype(bool))


# ---------------------------------------------------------------------------
# Differential expression
# ---------------------------------------------------------------------------

def rank_genes(
    X: Any,
    group: Any,
    *,
    gene_names: "list[str]",
    min_log2fc: float = MIN_LOG2FC,
    max_padj: float = MAX_PADJ,
    min_pct_in_group: float = MIN_PCT_IN_GROUP,
) -> "list[dict]":
    """Wilcoxon rank-sum of ``group`` vs the rest, one stratum.

    Mirrors ``scanpy.tl.rank_genes_groups(method='wilcoxon')`` followed
    by the paper's filters: log2FC > 0.5, adjusted p < 0.001, and the
    gene detected in at least 50% of the senescent cells.

    The ``min_pct_in_group`` filter matters more than it looks. Without
    it, a gene detected in three senescent cells and zero others passes
    the fold-change and p-value tests on noise alone.

    Parameters
    ----------
    X
        ``(n_cells, n_genes)`` log-normalised expression for one stratum.
    group
        Boolean mask of the test group (the senescent cells).
    gene_names
        Column names of ``X``.

    Returns
    -------
    List of ``{gene, log2fc, pvalue, padj, pct_in, pct_out}`` for genes
    that pass every filter, sorted by descending fold change.
    """
    np = _np()
    from scipy import stats

    try:
        from igvfagent._ts_stats import benjamini_hochberg
    except Exception:
        from _ts_stats import benjamini_hochberg

    A = np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=float)
    g = np.asarray(group, dtype=bool)
    n_in, n_out = int(g.sum()), int((~g).sum())
    if n_in < 3 or n_out < 3:
        return []

    inn, out = A[g], A[~g]
    pct_in = (inn > 0).mean(axis=0)
    pct_out = (out > 0).mean(axis=0)

    # Only test genes that could pass the detection filter; the Wilcoxon
    # is the expensive part and most genes are hopeless.
    cand = np.flatnonzero(pct_in >= min_pct_in_group)
    if cand.size == 0:
        return []

    mean_in = inn[:, cand].mean(axis=0)
    mean_out = out[:, cand].mean(axis=0)
    log2fc = np.log2((np.expm1(mean_in) + 1e-9) / (np.expm1(mean_out) + 1e-9))

    pvals = np.ones(cand.size)
    for j, c in enumerate(cand):
        a, b = inn[:, c], out[:, c]
        if a.std() == 0 and b.std() == 0:
            continue
        try:
            pvals[j] = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
        except ValueError:
            pvals[j] = 1.0
    padj = benjamini_hochberg(pvals)

    hits = []
    for j, c in enumerate(cand):
        if not np.isfinite(padj[j]) or padj[j] >= max_padj:
            continue
        if log2fc[j] <= min_log2fc:
            continue
        hits.append({
            "gene": gene_names[c],
            "log2fc": float(log2fc[j]),
            "pvalue": float(pvals[j]),
            "padj": float(padj[j]),
            "pct_in": float(pct_in[c]),
            "pct_out": float(pct_out[c]),
        })
    hits.sort(key=lambda h: -h["log2fc"])
    return hits


def replicated_sags(
    per_stratum: "dict[tuple, list[dict]]",
    *,
    min_donor_fraction: float = MIN_DONOR_FRACTION,
) -> "tuple[list[str], dict[str, dict]]":
    """Genes that clear the filters in >= half the donors of some stratum.

    ``per_stratum`` is keyed ``(donor, tissue, broad_cell_type)``. The
    paper requires replication across donors, not merely a good p-value
    in one: a stratum is one donor's cells, and single-donor findings are
    exactly what a 24-donor atlas is meant to guard against.

    Returns ``(sag_list, detail)`` where detail maps gene -> per-stratum
    support.
    """
    from collections import defaultdict

    donors_of: "dict[tuple, set]" = defaultdict(set)
    hits_of: "dict[tuple, dict[str, set]]" = defaultdict(lambda: defaultdict(set))
    for (donor, tissue, ct), hits in per_stratum.items():
        donors_of[(tissue, ct)].add(donor)
        for h in hits:
            hits_of[(tissue, ct)][h["gene"]].add(donor)

    detail: "dict[str, dict]" = {}
    sags: "set[str]" = set()
    for key, gene_map in hits_of.items():
        n_donor = len(donors_of[key])
        if n_donor == 0:
            continue
        need = max(1, int(round(min_donor_fraction * n_donor)))
        for gene, donors in gene_map.items():
            if len(donors) >= need:
                sags.add(gene)
                d = detail.setdefault(gene, {"strata": [], "n_strata": 0})
                d["strata"].append({"tissue": key[0], "cell_type": key[1],
                                    "donors": len(donors),
                                    "donors_total": n_donor})
                d["n_strata"] += 1
    return sorted(sags), detail


def cell_type_prevalence(
    per_cell_type_lfc: "dict[str, dict[str, float]]",
    *,
    min_log2fc: float = PREVALENCE_LOG2FC,
) -> "dict[str, int]":
    """How many cell types each gene is up in.

    The paper's Figure 4C axis. A gene "universal" to senescence would
    score near the number of cell types; the finding is that almost
    nothing does -- CDKN2B, the runner-up to CDKN2A, reaches only ~60%.
    """
    from collections import Counter
    out: "Counter[str]" = Counter()
    for _ct, lfcs in per_cell_type_lfc.items():
        for gene, lfc in lfcs.items():
            if lfc is not None and lfc > min_log2fc:
                out[gene] += 1
    return dict(out)


def hallmark_coverage(
    per_cell_type_lfc: "dict[str, dict[str, float]]",
    *,
    hallmarks: "Optional[dict[str, list[str]]]" = None,
    min_log2fc: float = PREVALENCE_LOG2FC,
) -> "dict[str, dict[str, Any]]":
    """Per cell type, which senescence hallmarks show any upregulated gene.

    The SenNet guidance the paper follows asks for a multi-pronged
    definition rather than a single marker. The paper's claim is that
    "in most cell types at least two additional hallmarks of senescence
    manifested themselves"; this quantifies that directly.
    """
    hm = hallmarks or HALLMARKS
    out: "dict[str, dict[str, Any]]" = {}
    for ct, lfcs in per_cell_type_lfc.items():
        present = {}
        for name, genes in hm.items():
            hit = [g for g in genes
                   if lfcs.get(g) is not None and lfcs[g] > min_log2fc]
            present[name] = hit
        out[ct] = {"hallmarks_present": sum(1 for v in present.values() if v),
                   "detail": present}
    return out
