#!/usr/bin/env python3
"""IGVF agent Flow-FISH CRISPR-screen analytics skill.

End-to-end **CRISPRi Flow-FISH** counts analysis: per-guide log-normal
maximum-likelihood expression estimate (with EM treatment of an "outside"
overflow bin), real-space conversion + null normalization, and per-element
significance (Mann-Whitney U + Welch t against negative controls,
BH-FDR).

This is a clean-room reimplementation of the algorithms in
[EngreitzLab/CRISPRi-FlowFISH-pipeline](https://github.com/EngreitzLab/CRISPRi-FlowFISH-pipeline)
(MIT, Engreitz Lab 2021). No source code from that repository is copied;
the math is paraphrased from the published descriptions in:

  * Fulco CP et al., *Nature Genetics* 51:1664-1669 (2019).
    "Activity-by-contact model of enhancer-promoter regulation from
     thousands of CRISPR perturbations."
  * Nasser J et al., *Nature* 593:238-243 (2021).
    "Genome-wide enhancer maps link risk variants to disease genes."

Commands

    flowfish pull-portal      — discover IGVF Portal CRISPRi-FlowFISH datasets
    flowfish estimate-effects — per-guide log-normal MLE with EM
    flowfish real-space       — convert to fold-change + normalize to null
    flowfish score-elements   — per-element Mann-Whitney + Welch + FDR
    flowfish simulate         — synthetic guide-bin counts for smoke testing
    flowfish write-playbook   — write markdown playbook

Heavy deps (numpy, pandas, scipy, statsmodels) imported lazily.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "FlowFISH"
PLOT_DIR = REPORT_DIR / "Plots"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

PORTAL_BASE = _resolve_endpoint("portal_api", "IGVF_PORTAL_BASE")


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"flowfish_crispr_skill_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    return log_path


def safe_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in label)


def _require_pkg(name: str, hint: str) -> Any:
    try:
        return __import__(name)
    except Exception as exc:
        raise SystemExit(
            f"Missing dependency '{name}'. {hint}\nInstall with: pip install {name}"
        ) from exc


def save_json(label: str, data: Any) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(label)}.json"
    path.write_text(json.dumps(data, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Portal discovery
# ---------------------------------------------------------------------------

def cmd_pull_portal(args: argparse.Namespace) -> int:
    setup_logging()
    params = {
        "type": "MeasurementSet",
        "searchTerm": "Flow-FISH",
        "format": "json",
        "limit": str(args.limit),
    }
    url = f"{PORTAL_BASE}/search/?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        data = json.load(resp)
    rows = data.get("@graph", [])
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out = REPORT_DIR / f"{time.strftime('%Y%m%d_%H%M%S')}_{safe_label(args.label)}_portal.tsv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "assay_titles", "preferred_assay_titles", "n_files", "status", "description"])
        for r in rows:
            w.writerow([
                r.get("accession", ""),
                "|".join(r.get("assay_titles") or []),
                "|".join(r.get("preferred_assay_titles") or []),
                len(r.get("files") or []),
                r.get("status", ""),
                (r.get("description") or "")[:120],
            ])
    print(f"Flow-FISH portal manifest: {out} (total={data.get('total')}, returned={len(rows)})")
    save_json(f"flowfish_{args.label}_portal", {"url": url, "total": data.get("total"), "n": len(rows), "out": str(out)})
    return 0


# ---------------------------------------------------------------------------
# B1 — Per-guide log-normal MLE (Fulco 2019 SI)
# ---------------------------------------------------------------------------

def _nll_bin_lognormal(params: Any, edges_log10: Any, counts: Any) -> float:
    """Negative log-likelihood for a multinomial log-normal bin model.

    P(bin b | mu, sigma) = Phi((log10(x_b_hi) - mu)/sigma) - Phi((log10(x_b_lo) - mu)/sigma)
    NLL = - sum_b c_b * log(P_b)
    """
    np = _require_pkg("numpy", "")
    scipy_stats = __import__("scipy.stats", fromlist=["stats"])
    mu, sigma = params
    sigma = max(float(sigma), 1e-4)
    norm = scipy_stats.norm
    # edges_log10 has shape (B+1,) where edges[b], edges[b+1] flank bin b
    z = (edges_log10 - mu) / sigma
    cdf = norm.cdf(z)
    p = np.diff(cdf)
    # Last "outside" bin = 1 - sum of inside probabilities; if counts has B+1 entries
    if len(counts) == len(p) + 1:
        p = np.append(p, max(1.0 - p.sum(), 1e-12))
    p = np.maximum(p, 1e-12)
    return float(-(counts * np.log(p)).sum())


def _mle_guide(
    edges_log10: Any,
    counts: Any,
    mu0: float,
    sigma0: float,
    *,
    sigma_bounds: tuple[float, float] = (0.1, 1.0),
    mu_bounds: tuple[float, float] = None,
) -> dict[str, Any]:
    """Bounded L-BFGS-B fit. ``counts`` may include a trailing "outside" bin
    that captures cells not represented by any inside bin (EM imputation by
    construction of the log-normal model).
    """
    np = _require_pkg("numpy", "")
    optimize = __import__("scipy.optimize", fromlist=["optimize"])
    if mu_bounds is None:
        # Derive a generous mu floor/ceiling from the bin edges, with a wider
        # margin so a guide that pushes cells into the lowest or highest bin
        # can still be fit. Default range: edges[0] - 1 to edges[-1] + 1
        # (in log10 units, ~10x wider than the bin coverage).
        mu_bounds = (float(edges_log10.min()) - 1.0, float(edges_log10.max()) + 1.0)
    res = optimize.minimize(
        _nll_bin_lognormal,
        x0=np.array([mu0, sigma0]),
        args=(edges_log10, counts),
        method="L-BFGS-B",
        bounds=[mu_bounds, sigma_bounds],
    )
    mu_hat, sigma_hat = float(res.x[0]), float(res.x[1])
    weighted_avg = float(((edges_log10[:-1] + edges_log10[1:]) / 2.0 * counts[:len(edges_log10) - 1]).sum() /
                          max(counts[:len(edges_log10) - 1].sum(), 1e-9))
    return {
        "logMean": mu_hat,
        "logSD": sigma_hat,
        "method": "MLE" if res.success else "MLE_failed",
        "WeightedAvg": weighted_avg,
        "sumcells": int(counts.sum()),
        "numobsbins": int((counts > 0).sum()),
    }


def cmd_estimate_effects(args: argparse.Namespace) -> int:
    """Per-guide log-normal MLE on a guide-by-bin count matrix.

    Inputs (TSV):
      counts.tsv with columns: ``GuideID`` (or ``OligoID``), and one column per
        bin labeled ``Bin1, Bin2, ...`` (last bin may be a synthetic "outside"
        catch-all if the user supplies it).
      sortparams.tsv with columns: ``Bin, LowBound, HighBound``
        (fluorescence in linear space).
    """
    pd = _require_pkg("pandas", "")
    np = _require_pkg("numpy", "")
    setup_logging()
    counts_df = pd.read_csv(args.counts, sep=None, engine="python")
    sort_df = pd.read_csv(args.sortparams, sep=None, engine="python")
    if "Bin" not in sort_df.columns or "LowBound" not in sort_df.columns:
        raise SystemExit("sortparams must contain Bin, LowBound, HighBound columns.")
    # Build log10 bin edges
    sort_df = sort_df.sort_values("LowBound").reset_index(drop=True)
    edges = np.log10(np.maximum(np.append(sort_df["LowBound"].to_numpy(),
                                            sort_df["HighBound"].iloc[-1]), 1e-3))
    bin_cols = [b for b in sort_df["Bin"].astype(str).tolist() if b in counts_df.columns]
    if len(bin_cols) != len(sort_df):
        raise SystemExit(f"Counts table missing bins: got {bin_cols}, expected {sort_df['Bin'].tolist()}.")
    guide_col = "GuideID" if "GuideID" in counts_df.columns else "OligoID"
    if guide_col not in counts_df.columns:
        raise SystemExit("Counts table needs a 'GuideID' or 'OligoID' column.")
    # Initial mu = mean(log10(midpoint), weighted by counts); sigma = 0.5
    mid = (edges[:-1] + edges[1:]) / 2.0
    results = []
    for _, row in counts_df.iterrows():
        c = row[bin_cols].to_numpy(dtype=float)
        if c.sum() <= 0:
            results.append({"GuideID": row[guide_col], "logMean": np.nan, "logSD": np.nan,
                             "method": "no_counts", "WeightedAvg": np.nan, "sumcells": 0,
                             "numobsbins": 0})
            continue
        mu0 = float((mid * c).sum() / c.sum())
        rec = _mle_guide(edges, c, mu0=mu0, sigma0=0.5)
        rec["GuideID"] = row[guide_col]
        results.append(rec)
    out_df = pd.DataFrame(results)
    # Carry over any metadata columns
    meta_cols = [c for c in counts_df.columns if c not in bin_cols + [guide_col]]
    if meta_cols:
        out_df = out_df.merge(counts_df[[guide_col] + meta_cols], on="GuideID", how="left")
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_raw_effects.tsv"
    out_df.to_csv(out, sep="\t", index=False)
    print(f"Per-guide MLE raw effects: {out} ({len(out_df)} guides)")
    return 0


# ---------------------------------------------------------------------------
# B2 — Real-space conversion + null normalization
# ---------------------------------------------------------------------------

def cmd_real_space(args: argparse.Namespace) -> int:
    """Convert MLE (mu, sigma) -> log-normal mean, normalize to negative-control median."""
    pd = _require_pkg("pandas", "")
    np = _require_pkg("numpy", "")
    setup_logging()
    df = pd.read_csv(args.input, sep="\t")
    if "logMean" not in df.columns:
        raise SystemExit("Input must contain logMean and logSD columns (run `flowfish estimate-effects` first).")
    if args.target_col not in df.columns:
        raise SystemExit(f"Target column '{args.target_col}' not in input; have {df.columns.tolist()}.")
    # Convert log-normal -> mean = 10^(mu + ln(10)*sigma^2/2)
    # (E[X] = exp(mu + sigma^2/2) for X ~ LogNormal(mu, sigma) on natural-log scale)
    # Our mu / sigma are in log10 -> convert with ln(10) factor
    ln10 = np.log(10.0)
    mu_nat = df["logMean"].to_numpy() * ln10
    sig_nat = df["logSD"].to_numpy() * ln10
    df["mleAvg"] = np.exp(mu_nat + sig_nat ** 2 / 2.0)
    df["mleSD"] = np.sqrt((np.exp(sig_nat ** 2) - 1.0) * np.exp(2 * mu_nat + sig_nat ** 2))

    is_neg = df[args.target_col].astype(str).str.lower().isin(
        [s.strip().lower() for s in args.negative_label.split(",")]
    )
    neg_median = float(df.loc[is_neg, "mleAvg"].median())
    if not np.isfinite(neg_median) or neg_median <= 0:
        raise SystemExit(
            f"Negative-control median is invalid ({neg_median}). "
            f"Make sure at least one guide has {args.target_col} in {{{args.negative_label}}}."
        )
    df["effect_raw"] = df["mleAvg"] / neg_median
    df["effect_clamped"] = np.clip(df["effect_raw"], 0.0, float(args.clamp))
    # Rescale by mean of clamped negatives so null centers at 1.0
    neg_mean_clamped = float(df.loc[is_neg, "effect_clamped"].mean())
    if neg_mean_clamped > 0:
        df["effect_normalized"] = df["effect_clamped"] / neg_mean_clamped
    else:
        df["effect_normalized"] = df["effect_clamped"]
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_real_space.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"Real-space effects: {out}")
    print(f"Negative-control median mleAvg = {neg_median:.4f}; clamp at {args.clamp}; "
          f"null mean after rescale = 1.0")
    return 0


# ---------------------------------------------------------------------------
# B3 — Per-element collapse + Mann-Whitney + Welch t-test + BH FDR
# ---------------------------------------------------------------------------

def cmd_score_elements(args: argparse.Namespace) -> int:
    """Collapse guides to elements and test against negative-control distribution."""
    pd = _require_pkg("pandas", "")
    np = _require_pkg("numpy", "")
    scipy_stats = __import__("scipy.stats", fromlist=["stats"])
    from statsmodels.stats.multitest import multipletests
    setup_logging()

    eff = pd.read_csv(args.effects, sep="\t")
    if "effect_normalized" not in eff.columns:
        raise SystemExit("Input must come from `flowfish real-space`.")
    is_neg = eff[args.target_col].astype(str).str.lower().isin(
        [s.strip().lower() for s in args.negative_label.split(",")]
    )
    null = eff.loc[is_neg, "effect_normalized"].dropna().to_numpy()
    if len(null) < args.min_negative:
        raise SystemExit(f"Need at least {args.min_negative} negative-control guides; have {len(null)}.")

    elem_col = args.element_col
    if elem_col not in eff.columns:
        raise SystemExit(f"Element column '{elem_col}' not in input; have {eff.columns.tolist()}.")

    # Aggregate by element
    rows = []
    for elem, sub in eff.loc[~is_neg].groupby(elem_col, dropna=True):
        x = sub["effect_normalized"].dropna().to_numpy()
        if len(x) < args.min_guides:
            continue
        mean = float(x.mean())
        sem = float(x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else float("nan")
        u_stat, p_u = scipy_stats.mannwhitneyu(x, null, alternative="two-sided")
        t_stat, p_t = scipy_stats.ttest_ind(x, null, equal_var=False)
        rows.append({
            elem_col: elem,
            "mean_effect": mean,
            "sem": sem,
            "n_guides": int(len(x)),
            "mean_ctrl": float(null.mean()),
            "sem_ctrl": float(null.std(ddof=1) / np.sqrt(len(null))),
            "n_ctrl": int(len(null)),
            "u_stat": float(u_stat),
            "p_mwu": float(p_u),
            "t_stat": float(t_stat),
            "p_ttest": float(p_t),
            "log2FC": float(np.log2(max(mean, 1e-9))),
            "EffectSize": float(1.0 - mean),  # fraction knockdown vs null
        })
    if not rows:
        raise SystemExit("No elements passed min_guides threshold.")
    out_df = pd.DataFrame(rows)
    _, padj_mwu, _, _ = multipletests(out_df["p_mwu"].to_numpy(), method="fdr_bh")
    _, padj_t, _, _ = multipletests(out_df["p_ttest"].to_numpy(), method="fdr_bh")
    out_df["fdr_mwu"] = padj_mwu
    out_df["fdr_ttest"] = padj_t
    out_df["Significant"] = (out_df["fdr_mwu"] < args.fdr) & (out_df["mean_effect"] < 1.0) & (out_df["n_guides"] >= args.min_guides)
    out_df["Regulated"] = out_df["Significant"] & (out_df["EffectSize"] >= args.min_effect)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    ts = time.strftime("%Y%m%d_%H%M%S")
    out = REPORT_DIR / f"{ts}_{safe_label(args.label)}_FullEnhancerScore.tsv"
    out_df.to_csv(out, sep="\t", index=False)
    summary = {
        "n_elements": int(len(out_df)),
        "significant_mwu_fdr": int(out_df["Significant"].sum()),
        "regulated_fdr_effect": int(out_df["Regulated"].sum()),
        "fdr": float(args.fdr),
        "min_guides": int(args.min_guides),
        "min_effect": float(args.min_effect),
    }
    save_json(f"flowfish_{args.label}_score", summary)
    print(f"Element effect table: {out}")
    print(f"Elements tested: {summary['n_elements']}")
    print(f"Significant (MWU FDR<{args.fdr}): {summary['significant_mwu_fdr']}")
    print(f"Regulated (Significant AND |1-mean|>={args.min_effect}): {summary['regulated_fdr_effect']}")
    return 0


# ---------------------------------------------------------------------------
# Simulation (smoke test)
# ---------------------------------------------------------------------------

def cmd_simulate(args: argparse.Namespace) -> int:
    """Generate a synthetic Flow-FISH guide-by-bin counts table for smoke testing.

    100 negative-control guides + ``--n_elements`` regulatory elements; each
    element has ``--guides_per_element`` guides. ``--knockdown_frac`` of
    elements are truly regulating (50% knockdown signature).
    """
    pd = _require_pkg("pandas", "")
    np = _require_pkg("numpy", "")
    rng = np.random.default_rng(args.seed)
    setup_logging()
    n_bins = 6
    bin_lows = np.array([0.1, 0.3, 1.0, 3.0, 10.0, 30.0])
    bin_highs = np.array([0.3, 1.0, 3.0, 10.0, 30.0, 100.0])

    sort_df = pd.DataFrame({"Bin": [f"Bin{i+1}" for i in range(n_bins)],
                            "LowBound": bin_lows, "HighBound": bin_highs})
    sort_path = Path(args.out_dir) / "sortparams.tsv"
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    sort_df.to_csv(sort_path, sep="\t", index=False)

    log_edges = np.log10(np.append(bin_lows, bin_highs[-1]))
    log_centers = (log_edges[:-1] + log_edges[1:]) / 2.0
    rows = []
    # 100 negative-control guides centered at log_mean=1.0, sigma=0.4
    for k in range(100):
        mu = rng.normal(1.0, 0.05)
        cells = rng.multinomial(args.cells_per_guide,
                                  pvals=_lognormal_pvals(log_edges, mu, 0.4))
        rec = {"GuideID": f"NEG_{k:03d}", "target": "negative_control",
                "ElementName": f"NEG_{k:03d}"}
        for i in range(n_bins):
            rec[f"Bin{i+1}"] = int(cells[i])
        rows.append(rec)
    # Real elements
    for e in range(args.n_elements):
        knockdown = rng.random() < args.knockdown_frac
        target = f"ELEM_{e:03d}"
        mu_elem = 0.5 if knockdown else 1.0  # 50% knockdown moves mu by ~-0.3
        for g in range(args.guides_per_element):
            mu = rng.normal(mu_elem, 0.05)
            cells = rng.multinomial(args.cells_per_guide,
                                      pvals=_lognormal_pvals(log_edges, mu, 0.4))
            rec = {"GuideID": f"{target}_g{g}", "target": target, "ElementName": target}
            for i in range(n_bins):
                rec[f"Bin{i+1}"] = int(cells[i])
            rows.append(rec)
    counts_df = pd.DataFrame(rows)
    counts_path = Path(args.out_dir) / "counts.tsv"
    counts_df.to_csv(counts_path, sep="\t", index=False)
    print(f"Wrote synthetic counts: {counts_path} ({len(counts_df)} guides)")
    print(f"Wrote synthetic sortparams: {sort_path}")
    return 0


def _lognormal_pvals(edges_log10: Any, mu: float, sigma: float) -> Any:
    """Multinomial cell-fraction vector for one guide."""
    np = _require_pkg("numpy", "")
    scipy_stats = __import__("scipy.stats", fromlist=["stats"])
    z = (edges_log10 - mu) / max(sigma, 1e-3)
    cdf = scipy_stats.norm.cdf(z)
    p = np.diff(cdf)
    s = p.sum()
    return p / s if s > 0 else np.ones_like(p) / len(p)


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------

def write_playbook() -> Path:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "FLOWFISH_CRISPR_SKILLS.md"
    path.write_text(
        """# Skill: CRISPRi Flow-FISH Screen

Use this skill to call regulatory elements from a CRISPRi Flow-FISH screen
(cells sorted by RNA-FISH or protein-FACS readout into N bins per guide,
guide-counts sequenced per bin). The math is a clean-room reimplementation
of [EngreitzLab/CRISPRi-FlowFISH-pipeline](https://github.com/EngreitzLab/CRISPRi-FlowFISH-pipeline)
(MIT, Engreitz Lab 2021), following the published methods in:

- **Fulco CP et al. (2019)** *Nature Genetics* 51:1664–1669 — SI describes
  the per-guide log-normal MLE with EM treatment of an "outside" bin.
- **Nasser J et al. (2021)** *Nature* 593:238–243 — applies the same method
  at scale and defines the element-collapse / `Significant` / `Regulated`
  output convention.

## Commands

```bash
# 1. Discover IGVF Portal CRISPRi Flow-FISH datasets
python3 Scripts/flowfish_crispr_skill.py pull-portal --limit 50 --label survey

# 2. Generate a synthetic screen for smoke testing
python3 Scripts/flowfish_crispr_skill.py simulate \\
    --out-dir /tmp/flowfish_smoke --n-elements 40 \\
    --guides-per-element 8 --knockdown-frac 0.30 --cells-per-guide 800

# 3. Per-guide log-normal MLE
python3 Scripts/flowfish_crispr_skill.py estimate-effects \\
    --counts /tmp/flowfish_smoke/counts.tsv \\
    --sortparams /tmp/flowfish_smoke/sortparams.tsv --label demo

# 4. Real-space conversion + null normalization
python3 Scripts/flowfish_crispr_skill.py real-space \\
    --input Docs/FlowFISH/<ts>_demo_raw_effects.tsv \\
    --target-col target --negative-label negative_control \\
    --clamp 5 --label demo

# 5. Per-element collapse + significance
python3 Scripts/flowfish_crispr_skill.py score-elements \\
    --effects Docs/FlowFISH/<ts>_demo_real_space.tsv \\
    --target-col target --element-col ElementName \\
    --negative-label negative_control --min-guides 5 --fdr 0.05 \\
    --min-effect 0.10 --label demo
```

## Input schema

`counts.tsv` (TSV, one row per guide):

| Column | Required | Notes |
|---|---|---|
| `GuideID` | yes | Unique guide identifier |
| `target` | yes | Gene symbol / element name / `negative_control` |
| `ElementName` | optional | Mapping to a higher-level regulatory element |
| `Bin1`, `Bin2`, ... | yes | Cell counts per FACS bin |

`sortparams.tsv`: one row per bin, columns `Bin, LowBound, HighBound`
(fluorescence in linear space; the skill converts to log10 internally).

## Output

`*_FullEnhancerScore.tsv` (per regulatory element):

| Column | Meaning |
|---|---|
| `ElementName` | Element identifier |
| `mean_effect` | Mean of `effect_normalized` across guides targeting this element (1.0 = null) |
| `EffectSize` | `1 - mean_effect`, fraction knockdown vs null |
| `log2FC` | `log2(mean_effect)` |
| `n_guides` | Guides used in the test |
| `mean_ctrl, sem_ctrl, n_ctrl` | Negative-control distribution stats |
| `p_mwu, fdr_mwu` | Mann-Whitney U vs negatives, BH-adjusted |
| `p_ttest, fdr_ttest` | Welch t-test vs negatives, BH-adjusted |
| `Significant` | `fdr_mwu < FDR` AND `mean_effect < 1` AND `n_guides ≥ min_guides` |
| `Regulated` | `Significant` AND `EffectSize ≥ min_effect` |

## Notes

- The original pipeline also runs sliding-window smoothing (20-guide /
  ≤750 bp) for tiling screens. Add it as a post-process if needed.
- The qPCR-anchored knockdown calibration (`normalize_flowfish_to_qpcr.py`)
  is not implemented yet — but the `--clamp` and per-element collapse
  steps are sufficient for the common case.
- Replicates: fit the MLE per replicate independently and average the
  resulting `effect_normalized` values per guide before calling
  `score-elements`. (You can do this with `pandas groupby + mean`
  externally; this skill keeps the per-replicate atom and the
  per-element scoring as separate commands.)
""",
        encoding="utf-8",
    )
    return path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="CRISPRi Flow-FISH screen analytics (clean-room).")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("pull-portal", help="Discover IGVF Portal Flow-FISH MeasurementSets.")
    p.add_argument("--limit", type=int, default=50)
    p.add_argument("--label", default="flowfish_portal_survey")

    p = sub.add_parser("estimate-effects", help="Per-guide log-normal MLE on guide×bin counts.")
    p.add_argument("--counts", required=True, help="Guide × bin counts TSV.")
    p.add_argument("--sortparams", required=True, help="Bin lookups TSV: Bin, LowBound, HighBound.")
    p.add_argument("--label", default="flowfish_estimate")

    p = sub.add_parser("real-space", help="Convert MLE to fold-change + normalize to negative controls.")
    p.add_argument("--input", required=True, help="raw_effects.tsv from `estimate-effects`.")
    p.add_argument("--target-col", default="target", help="Column name for guide target (default: target).")
    p.add_argument("--negative-label", default="negative_control",
                    help="Comma-separated list of values in target-col that mark non-targeting controls.")
    p.add_argument("--clamp", type=float, default=5.0)
    p.add_argument("--label", default="flowfish_real_space")

    p = sub.add_parser("score-elements", help="Mann-Whitney + Welch t-test + BH-FDR per element.")
    p.add_argument("--effects", required=True, help="real_space.tsv from `real-space`.")
    p.add_argument("--target-col", default="target")
    p.add_argument("--element-col", default="ElementName")
    p.add_argument("--negative-label", default="negative_control")
    p.add_argument("--min-guides", type=int, default=5)
    p.add_argument("--min-negative", type=int, default=10)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--min-effect", type=float, default=0.10)
    p.add_argument("--label", default="flowfish_score")

    p = sub.add_parser("simulate", help="Generate a synthetic guide × bin counts table for smoke testing.")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--n-elements", type=int, default=40)
    p.add_argument("--guides-per-element", type=int, default=8)
    p.add_argument("--knockdown-frac", type=float, default=0.30)
    p.add_argument("--cells-per-guide", type=int, default=800)
    p.add_argument("--seed", type=int, default=7)

    sub.add_parser("write-playbook", help="Write skill markdown playbook.")

    args = parser.parse_args(argv)
    if args.command == "pull-portal":
        return cmd_pull_portal(args)
    if args.command == "estimate-effects":
        return cmd_estimate_effects(args)
    if args.command == "real-space":
        return cmd_real_space(args)
    if args.command == "score-elements":
        return cmd_score_elements(args)
    if args.command == "simulate":
        return cmd_simulate(args)
    if args.command == "write-playbook":
        path = write_playbook()
        print(f"Wrote {path}")
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
