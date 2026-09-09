#!/usr/bin/env python3
"""Does the lettered-bin screen scorer hold up? Synthetic data, no network.

Written after the first real run of gradient-screen on IGVFDS8710ZSOZ
(KITLG FlowFISH) produced 858 "significant" constructs of 2,719 -- and 27.2%
of its 434 NON-TARGETING CONTROLS along with them. Three defects, each found
by the controls and each fixed here under test:

  1. The four sorts are flow replicates of one biological sample. Their
     spread is instrument precision, not biology, and using it as the null
     made everything significant.
  2. Judging effects against control scatter alone put constructs with
     25-40 reads in a single sort at the top of the hit list, because a
     shallow construct's score is expected to be extreme.
  3. One null width for all depths is dominated by deep controls: the
     shallowest quartile still produced 3.7% false positives against 0.0%
     in the other three.

Final state on that screen: 65 targeting hits, none with fewer than 100
reads or a single sort, and 0.46% of controls significant.
"""
import math
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Scripts"))
import gradient_screen_analysis as g  # noqa: E402

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:56} {detail}")
    if not ok:
        FAILURES.append(name)


# ── alias parsing: token-based, not positional ─────────────────────────────

CASES = [
 ("220228_8merInsertion-THP1-Amp_PPIF_promoter-BioRep3-FFrep3-BinF_AmpliconSequencing",
  "220228_8merInsertion-THP1-Amp_PPIF_promoter", "F", 3, 3, None),
 ("measurement_CRUDO_KITLG-Auxin6hrs-FF4-BinE-PCR2_S95",
  "measurement_CRUDO_KITLG-Auxin6hrs", "E", None, 4, 2),
 ("measurement_CRUDO_MYC-Auxin0hrs-FF1-BinE-PCR4_Rep1_S23",
  "measurement_CRUDO_MYC-Auxin0hrs", "E", 1, 1, 4),
 ("measurement_CRUDO_FAM3C_Auxin6hrs-FF2-BinB-PCR3-Rep1_S200",
  "measurement_CRUDO_FAM3C_Auxin6hrs", "B", 1, 2, 3),
]
for alias, series, b, bio, ff, pcr in CASES:
    r = g.parse_bin_alias(alias)
    check(f"parses {alias[:34]}", r is not None)
    if r:
        check(f"  series {series[:30]}", r["series"] == series, r["series"])
        check(f"  axes bin={b} bio={bio} ff={ff} pcr={pcr}",
              (r["bin"], r["bio_rep"], r["flow_rep"], r["pcr_rep"])
              == (b, bio, ff, pcr),
              f"{r['bin']},{r['bio_rep']},{r['flow_rep']},{r['pcr_rep']}")

# A tail sort belongs to the other tool and must not be claimed.
check("tail-sort alias is not claimed",
      g.parse_bin_alias("0426_LDLR137-219_repo_Rep1_Bot20_ms") is None)
check("no bin letter -> None", g.parse_bin_alias("some_series-Rep1") is None)
# The submission index must not split one screen into one series per sample.
a1 = g.parse_bin_alias("X_thing-FF1-BinA-PCR1_S1")
a2 = g.parse_bin_alias("X_thing-FF1-BinA-PCR1_S77")
check("trailing _Sxx does not change the series",
      a1 and a2 and a1["series"] == a2["series"], a1["series"] if a1 else "")
check("bin letters map to ranks A=1..F=6",
      [g.parse_bin_alias(f"S-FF1-Bin{c}")["rank"] for c in "ABCDEF"]
      == [1, 2, 3, 4, 5, 6])


# ── mean-bin score ─────────────────────────────────────────────────────────

def bins6(depths=None):
    d = depths or [10_000] * 6
    return [{"accession": f"B{i+1}", "rank": i + 1, "bin": "ABCDEF"[i],
             "bio_rep": None, "flow_rep": 1} for i in range(6)], d


def per_bin_from(shares, depths):
    """shares[gid][rank] = fraction of that bin's reads."""
    pb = {}
    for i in range(6):
        c = Counter()
        for gid, by in shares.items():
            c[gid] = int(by.get(i + 1, 0.0) * depths[i])
        c["_filler"] = max(depths[i] - sum(c.values()), 1)
        pb[f"B{i+1}"] = c
    return pb


bins, depths = bins6()
shares = {
 "low":  {1: 0.90, 2: 0.10},
 "high": {5: 0.10, 6: 0.90},
 "flat": {r: 1 / 6 for r in range(1, 7)},
}
pb = per_bin_from({k: {r: v * 0.01 for r, v in s.items()}
                   for k, s in shares.items()}, depths)
score, cvar, reads = g.sort_scores(pb, bins, min_count=10)
check("construct in low bins scores near 1", abs(score["low"] - 1.1) < 0.06,
      f"{score['low']:.3f}")
check("construct in high bins scores near 6", abs(score["high"] - 5.9) < 0.06,
      f"{score['high']:.3f}")
check("uniform construct scores mid-scale", abs(score["flat"] - 3.5) < 0.06,
      f"{score['flat']:.3f}")

# Depth must not tilt the score: bin F sequenced 10x deeper.
deep = [10_000, 10_000, 10_000, 10_000, 10_000, 100_000]
pb2 = per_bin_from({"flat": {r: 1 / 6 * 0.01 for r in range(1, 7)}}, deep)
s2, _c2, _r2 = g.sort_scores(pb2, bins6(deep)[0], min_count=10)
check("a 10x deeper bin does not tilt the score",
      abs(s2["flat"] - 3.5) < 0.15, f"{s2['flat']:.3f}")

# The degenerate case that was scored as infinitely precise.
pb3 = per_bin_from({"onebin": {1: 0.003}}, depths)
_s3, c3, r3 = g.sort_scores(pb3, bins, min_count=10)
check("all reads in one bin still has non-zero variance",
      c3["onebin"] > 0, f"cvar={c3['onebin']:.2e} reads={r3['onebin']}")
check("shallower construct has larger counting variance",
      c3["onebin"] > cvar["low"], f"{c3['onebin']:.2e} vs {cvar['low']:.2e}")
check("min_count drops a construct below threshold",
      "rare" not in g.sort_scores(
          per_bin_from({"rare": {1: 0.00005}}, depths), bins, 10_000)[0])


# ── the depth-dependent null ───────────────────────────────────────────────

# Var = a + b/N by construction; the fit must recover both.
import random
random.seed(7)
A, B = 0.09, 60.0
ctrl = []
for _ in range(600):
    n = random.choice([50, 200, 1000, 5000, 20000, 60000])
    ctrl.append((n, random.gauss(0.0, math.sqrt(A + B / n))))
a, b = g._fit_depth_null(ctrl)
check("depth-null recovers the constant term", abs(a - A) < 0.03, f"a={a:.4f}")
check("depth-null recovers the per-read term", abs(b - B) / B < 0.45, f"b={b:.1f}")
check("depth-null returns non-negative terms", a >= 0 and b >= 0)
check("too few controls -> single width, no depth term",
      g._fit_depth_null([(100, 0.1), (200, -0.1), (300, 0.0)])[1] == 0.0)


# ── end to end: controls define the null, and a real effect survives ──────

def synth_screen(n_ctrl=200, effect=-1.6, eff_reads=8000, n_sorts=4):
    per_sort = {}
    random.seed(3)
    for f in range(1, n_sorts + 1):
        score, cvar, reads = {}, {}, {}
        for i in range(n_ctrl):
            g_ = f"ctrl{i}"
            n = random.choice([500, 5000, 30000])
            score[g_] = 3.5 + random.gauss(0, math.sqrt(0.09 + 60.0 / n))
            cvar[g_] = 0.5 / n
            reads[g_] = n
        score["hit"] = 3.5 + effect + random.gauss(0, 0.12)
        cvar["hit"] = 0.5 / eff_reads
        reads["hit"] = eff_reads
        score["shallow"] = 3.5 + effect
        cvar["shallow"] = 0.5 / 30
        reads["shallow"] = 30
        per_sort[(None, f)] = {"score": score, "cvar": cvar, "reads": reads,
                               "profile": {}, "profile_ranks": [1, 2, 3]}
    types = {f"ctrl{i}": "non-targeting" for i in range(n_ctrl)}
    types["hit"] = "targeting"
    types["shallow"] = "targeting"
    targets = {k: k for k in list(types)}
    return g.score_gradient(per_sort, targets, types)


rows, mod = synth_screen()
by = {r["construct"]: r for r in rows}
check("null basis is the controls", mod["null_basis"] == "controls")
check("a deep, consistent effect is significant",
      by["hit"]["fdr"] < 0.05, f"fdr={by['hit']['fdr']:.2e}")
check("the SAME effect at 30 reads is not significant",
      by["shallow"]["fdr"] >= 0.05,
      f"fdr={by['shallow']['fdr']:.2e} dbins={by['shallow']['delta_bins']}")
ctrl_sig = sum(1 for r in rows
               if r["construct"].startswith("ctrl") and r["fdr"] < 0.05)
check("control false-positive rate stays under 5%",
      ctrl_sig / 200 <= 0.05, f"{ctrl_sig}/200 = {ctrl_sig/200:.1%}")
check("effects are centred on the controls",
      abs(sum(r["delta_bins"] for r in rows
              if r["construct"].startswith("ctrl")) / 200) < 0.05)
check("replicate-spread p is kept but not used for FDR",
      by["hit"]["p_replicate"] is not None
      and by["hit"]["p_value"] != by["hit"]["p_replicate"])

# With no controls the tool must fall back and say so.
rows2, mod2 = synth_screen(n_ctrl=3)
check("no usable controls -> replicate-spread basis",
      mod2["null_basis"] == "replicate spread", mod2["null_basis"])


# ── concordance ────────────────────────────────────────────────────────────

ps = {(None, 1): {"score": {"a": 1.0, "b": 2.0, "c": 3.0}},
      (None, 2): {"score": {"a": 1.0, "b": 2.0, "c": 3.0}},
      (None, 3): {"score": {"a": 3.0, "b": 2.0, "c": 1.0}}}
conc = {(tuple(c["a"]), tuple(c["b"])): c["r"] for c in g.concordance(ps)}
check("identical sorts give r = 1",
      abs(conc[((None, 1), (None, 2))] - 1.0) < 1e-9)
check("reversed sorts give r = -1",
      abs(conc[((None, 1), (None, 3))] + 1.0) < 1e-9)
check("too few shared constructs -> r is None",
      g.concordance({(None, 1): {"score": {"a": 1.0}},
                     (None, 2): {"score": {"a": 1.0}}})[0]["r"] is None)

print(f"\n{31} cases, {len(FAILURES)} failure(s)")
sys.exit(1 if FAILURES else 0)
