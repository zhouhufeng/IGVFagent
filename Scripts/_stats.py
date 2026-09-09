#!/usr/bin/env python3
"""Statistics shared by the screen analyses.

These live here rather than in one screen module because two analyses need
them -- the tail-comparison screens (crispr_screen_analysis) and the
expression-gradient screens (allelic_screen_analysis) -- and a second copy
of a test statistic is a second copy that can drift. The first version of
this arithmetic reported p = 0 from two replicates agreeing by luck; it
should exist exactly once.

Benchmarks/test_screen_statistics.py covers all of it.
"""
from __future__ import annotations

import math

def _norm_sf(z: float) -> float:
    """Two-sided normal tail probability, via erfc. No scipy needed."""
    return math.erfc(abs(z) / math.sqrt(2.0))


def _t_sf(t: float, df: int) -> float:
    """Two-sided Student-t tail probability.

    The normal is the wrong reference here and it is not a subtle error. A
    variant with two replicate observations that happen to agree to within
    0.016 log2 units gets se = 0.011 and z = -215, which the normal reports
    as p = 0 exactly -- infinite confidence from two numbers. The t
    distribution with df = k-1 is the textbook correction: at df = 1 it is
    Cauchy, whose tails are heavy enough that |t| = 215 is p ~ 3e-3 rather
    than 0.

    scipy is present in both the local environment and the deployed image,
    but the fallback matters: without it this would silently revert to the
    normal and to implausible p-values, so it degrades to a documented
    df-scaled approximation instead of pretending nothing changed.
    """
    if df < 1 or not math.isfinite(t):
        return 1.0
    try:
        from scipy import stats
        return float(2.0 * stats.t.sf(abs(t), df))
    except ImportError:
        pass
    if df == 1:                       # Cauchy, exactly
        return 2.0 * (0.5 - math.atan(abs(t)) / math.pi)
    if df == 2:                       # closed form
        return 1.0 - abs(t) / math.sqrt(2.0 + t * t)
    # Otherwise the normal on a variance-inflated statistic. Approximate,
    # and conservative in the direction that matters: it does not turn a
    # small sample into certainty.
    return _norm_sf(abs(t) / math.sqrt(df / max(df - 2.0, 1.0)))


def _percentile(xs: "list[float]", q: float) -> float:
    """Linear-interpolated percentile. Small helper, avoids a numpy import
    in a function that otherwise only needs the standard library."""
    if not xs:
        return 0.0
    ys = sorted(xs)
    if len(ys) == 1:
        return ys[0]
    i = (len(ys) - 1) * q
    lo, hi = int(math.floor(i)), int(math.ceil(i))
    return ys[lo] + (ys[hi] - ys[lo]) * (i - lo)


def benjamini_hochberg(pvals: "list[float]") -> "list[float]":
    """BH-adjusted p-values, order preserved."""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [0.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = n - rank + 1
        val = min(prev, pvals[idx] * n / i)
        adj[idx] = val
        prev = val
    return adj


def moderated_t(vals: "list[float]", sd_floor: float,
                 sd_counting: float = 0.0) -> dict:
    """Mean, moderated sd, t and p for one target's replicate observations.

    The moderation is the whole point. An sd computed from two or four
    observations is a weak estimate, and when those observations happen to
    agree it collapses toward zero and hands back an arbitrarily large t.
    Three lower bounds are applied, whichever binds hardest:

      the target's own empirical sd,
      `sd_counting`, from the read depths behind the observations, which is
        a per-target bound needing no reference distribution at all, and
      `sd_floor`, a low percentile of the sds seen across the whole screen.

    The reference is Student t with df = k-1, not the normal: with k = 2 the
    normal turns a coincidence into certainty.
    """
    k = len(vals)
    if k == 0:
        return {"n": 0, "mean": float("nan"), "sd": None, "sd_used": None,
                "t": 0.0, "df": 0, "p_value": 1.0}
    mean = sum(vals) / k
    if k > 1:
        var = sum((v - mean) ** 2 for v in vals) / (k - 1)
        sd = math.sqrt(var)
        sd_used = max(sd, sd_floor, sd_counting)
        se = sd_used / math.sqrt(k)
        tstat = mean / se if se > 0 else 0.0
        p = _t_sf(tstat, k - 1) if se > 0 else 1.0
    else:
        sd = sd_used = float("nan")
        tstat, p = 0.0, 1.0
    return {"n": k, "mean": mean,
            "sd": None if math.isnan(sd) else sd,
            "sd_used": None if math.isnan(sd_used) else sd_used,
            "t": tstat, "df": max(k - 1, 0), "p_value": p}
