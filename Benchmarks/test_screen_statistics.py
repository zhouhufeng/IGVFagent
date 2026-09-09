#!/usr/bin/env python3
"""Are the CRISPR screen p-values defensible at small replicate counts?

The first version of score_screen used a normal reference and no variance
floor. On the LDLR screen it reported its top hit -- a variant with TWO
replicate observations that agreed to within 0.016 log2 units -- at
p = 0.00e+00, and 41 targets under FDR 0.05. Two numbers agreeing by luck
were being read as infinite confidence.

Nothing caught it, because every test was about routing and parsing. These
cases are about the arithmetic, and they need no network.
"""
import math
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import Scripts                                     # noqa: E402
sys.modules.setdefault("igvfagent", Scripts)
from Scripts import crispr_screen_analysis as cs   # noqa: E402

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:52} {detail}")
    if not ok:
        FAILURES.append(name)


# ── the reference distribution ──────────────────────────────────────────────

# The exact case that produced p = 0. df = 1 is Cauchy; its tails must stay
# fat enough that a two-point agreement cannot manufacture certainty.
p1 = cs._t_sf(215.4, 1)
check("|t|=215, df=1 is not p=0", p1 > 1e-4, f"p={p1:.3e}")
check("|t|=215, df=1 beats the normal by orders of magnitude",
      p1 > cs._norm_sf(215.4) * 1e9, f"t={p1:.2e} vs normal={cs._norm_sf(215.4):.2e}")

# Heavier tails than the normal at every small df, converging upward.
check("t tails exceed normal tails at df=1..8",
      all(cs._t_sf(3.0, df) > cs._norm_sf(3.0) for df in range(1, 9)))
check("t p-value decreases as df grows (tails thin out)",
      all(cs._t_sf(3.0, df) >= cs._t_sf(3.0, df + 1) - 1e-12
          for df in range(1, 20)))
check("t p-value decreases as |t| grows",
      all(cs._t_sf(t, 3) > cs._t_sf(t + 0.5, 3) for t in (0.5, 1, 2, 4, 8)))
check("df<1 yields p=1 rather than a crash", cs._t_sf(5.0, 0) == 1.0)
check("non-finite t yields p=1", cs._t_sf(float("inf"), 3) == 1.0)
check("p stays a probability", all(0.0 <= cs._t_sf(t, df) <= 1.0
                                   for t in (0, .1, 1, 10, 1e6)
                                   for df in (1, 2, 3, 10)))

# Known values, so a future rewrite cannot quietly change the reference.
check("df=1 matches Cauchy closed form",
      abs(cs._t_sf(1.0, 1) - 0.5) < 1e-9, f"{cs._t_sf(1.0, 1):.6f}")
check("df=2 matches its closed form",
      abs(cs._t_sf(2.0, 2) - (1 - 2 / math.sqrt(6))) < 1e-9)


# ── percentile helper ──────────────────────────────────────────────────────

check("percentile of a singleton", cs._percentile([4.0], 0.1) == 4.0)
check("percentile of an empty list is 0", cs._percentile([], 0.5) == 0.0)
check("median of 1..5", cs._percentile([1, 2, 3, 4, 5], 0.5) == 3)
check("10th percentile interpolates",
      abs(cs._percentile([0, 10], 0.10) - 1.0) < 1e-9)


# ── score_screen end to end, on synthetic counts ───────────────────────────

def synth(per_rep, reps=(1, 2, 3, 4), total=100_000):
    """Build (per_bin, bins) where per_rep[gid] = [(lo, hi), ...] per rep."""
    bins, per_bin = [], {}
    for i, r in enumerate(reps):
        for side in ("bottom", "top"):
            acc = f"ACC{r}{side[:3]}"
            bins.append({"accession": acc, "rep": r, "side": side, "pct": 20})
            per_bin[acc] = Counter()
    for gid, obs in per_rep.items():
        for r, (lo, hi) in zip(reps, obs):
            per_bin[f"ACC{r}bot"][gid] = lo
            per_bin[f"ACC{r}top"][gid] = hi
    # a filler construct so the bin totals are realistic
    for acc in per_bin:
        per_bin[acc]["filler"] = total
    return per_bin, bins


# One variant seen in only TWO replicates, with near-identical ratios: the
# shape of the bogus top hit. It must not come out as overwhelming evidence.
per_bin, bins = synth({"g_quiet": [(2000, 400), (2010, 402), (0, 0), (0, 0)]})
tgt = {"g_quiet": "V_quiet", "filler": "V_filler"}
rows, mod = cs.score_screen(per_bin, bins, tgt, {}, 20, 10)
quiet = next(r for r in rows if r["target"] == "V_quiet")
check("two agreeing replicates: p is not ~0",
      quiet["p_value"] > 1e-3,
      f"n={quiet['n_guide_obs']} p={quiet['p_value']:.3e}")
check("two agreeing replicates: sd was floored",
      quiet["sd_used"] is not None and quiet["sd"] is not None
      and quiet["sd_used"] >= quiet["sd"],
      f"sd={quiet['sd']} -> {quiet['sd_used']}")
check("df is reported as k-1", quiet["df"] == quiet["n_guide_obs"] - 1)

# A real, consistent effect across four replicates should still be found.
per_bin, bins = synth({"g_real": [(4000, 500), (3800, 520), (4200, 480),
                                  (3900, 510)]})
rows, mod = cs.score_screen(per_bin, bins,
                            {"g_real": "V_real", "filler": "V_filler"}, {},
                            20, 10)
real = next(r for r in rows if r["target"] == "V_real")
check("a strong 4-replicate effect is still significant",
      real["p_value"] < 0.05 and real["mean_log2_low_over_high"] > 1,
      f"log2={real['mean_log2_low_over_high']:.2f} p={real['p_value']:.3e}")

# A construct with no real difference must not be called.
per_bin, bins = synth({"g_null": [(1000, 1000), (900, 1100), (1100, 950),
                                  (1000, 1050)]})
rows, mod = cs.score_screen(per_bin, bins,
                            {"g_null": "V_null", "filler": "V_filler"}, {},
                            20, 10)
null = next(r for r in rows if r["target"] == "V_null")
check("a null variant is not significant", null["p_value"] > 0.05,
      f"log2={null['mean_log2_low_over_high']:.2f} p={null['p_value']:.3e}")

# Replicates missing a tail must be excluded, not silently half-counted.
per_bin, bins = synth({"g": [(1000, 500)] * 4})
bins = [b for b in bins if not (b["rep"] == 4 and b["side"] == "top")]
rows, mod = cs.score_screen(per_bin, bins, {"g": "V", "filler": "F"}, {},
                            20, 10)
v = next(r for r in rows if r["target"] == "V")
check("a replicate lacking one tail is excluded", v["n_guide_obs"] == 3,
      f"n={v['n_guide_obs']}")

# BH must stay monotone and bounded.
qs = cs.benjamini_hochberg([0.001, 0.01, 0.02, 0.5, 0.9])
check("BH is monotone non-decreasing",
      all(qs[i] <= qs[i + 1] + 1e-12 for i in range(len(qs) - 1)), f"{qs}")
check("BH q-values stay <= 1", all(q <= 1.0 + 1e-12 for q in qs))

print(f"\n{22} cases, {len(FAILURES)} failure(s)")
sys.exit(1 if FAILURES else 0)
