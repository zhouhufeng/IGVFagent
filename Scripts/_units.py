"""Range checks for statistical fields, so an impossible number cannot ship.

A hosted retest reported the agent describing a value as "`p_adj` capped at
240". An adjusted P value cannot exceed 1. The underlying tool was right — it
returned `neg_log10_pvalue`, for which 240 is perfectly ordinary — and the
error appeared only in the final prose, where the field was renamed to
something it was not.

That failure mode is specific and worth catching mechanically: nothing about
"240" looks wrong until you know which field it belongs to, so a reader has no
chance, while the *name plus the value* is contradictory on its face.

Two entry points:

* :func:`check_value` — one field, one number.
* :func:`scan_text` — a finished answer; finds `field = value` claims in prose
  or markdown tables and reports the impossible ones.

This deliberately does NOT rewrite anything. It reports, and the caller
decides: the agent is told to correct itself, and a reviewer sees the same
list. Silent correction would hide exactly the confusion that needs fixing.
"""
from __future__ import annotations

import math
import re

# field name -> (low, high, what it is, what the value probably was)
#
# Only fields whose range is a matter of definition, not convention. A p-value
# cannot be 240 in any pipeline; an expression count can be anything, so
# counts are absent here on purpose.
RANGES: "dict[str, tuple[float, float, str, str]]" = {
    "p_adj":        (0.0, 1.0, "adjusted P value", "neg_log10_pvalue"),
    "padj":         (0.0, 1.0, "adjusted P value", "neg_log10_pvalue"),
    "p_value":      (0.0, 1.0, "P value", "neg_log10_pvalue"),
    "pvalue":       (0.0, 1.0, "P value", "neg_log10_pvalue"),
    "p_val":        (0.0, 1.0, "P value", "neg_log10_pvalue"),
    "fdr":          (0.0, 1.0, "false discovery rate", "neg_log10_pvalue"),
    "q_value":      (0.0, 1.0, "q value", "neg_log10_pvalue"),
    "qvalue":       (0.0, 1.0, "q value", "neg_log10_pvalue"),
    "neg_log10_pvalue": (0.0, math.inf, "-log10 P", "a raw P value"),
    "abc_score":    (0.0, 1.0, "ABC score (a share)", "a raw signal"),
    "probability":  (0.0, 1.0, "probability", "a percentage"),
    "fraction":     (0.0, 1.0, "fraction", "a percentage"),
    "percent":      (0.0, 100.0, "percentage", "a fraction"),
    "correlation":  (-1.0, 1.0, "correlation", "a covariance"),
    "pearson_r":    (-1.0, 1.0, "Pearson r", "a covariance"),
    "spearman_rho": (-1.0, 1.0, "Spearman rho", "a covariance"),
    "auc":          (0.0, 1.0, "AUC", "a percentage"),
    "auroc":        (0.0, 1.0, "AUROC", "a percentage"),
    "mapping_rate": (0.0, 1.0, "mapping rate", "a percentage"),
    "duplication_rate": (0.0, 1.0, "duplication rate", "a percentage"),
}

# Scanning is driven FROM the known field names, not by a generic
# "word = number" pattern. A generic pattern loses: in "the summary shows
# p_adj capped at 240" it matches `summary` as the field, finds 240, discards
# the pair because `summary` has no declared range, and moves past the number
# -- so the one claim that mattered is never examined. Anchoring on the
# handful of names that HAVE a range removes that competition entirely.
_NUMBER = r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"

# field, then at most a short run of non-numeric text, then the number. The
# filler is non-greedy and the sign lives in the number, so "= -3" yields -3
# rather than 3 with the minus quietly eaten by the filler.
_AFTER_FIELD = re.compile(
    r"[^0-9\n]{0,24}?(?P<value>" + _NUMBER + r")")


def _table_claims(text: str) -> "list[tuple[str, str]]":
    """(field, value) pairs from markdown tables, read via the header row.

    A table row is `| GATA3 | 240 |`; which column 240 belongs to is stated
    only in the header, so a row read on its own attributes the number to the
    gene name. The header is the thing that names the field.
    """
    out: "list[tuple[str, str]]" = []
    header: "list[str]" = []
    for line in (text or "").splitlines():
        line = line.strip()
        if not line.startswith("|"):
            header = []
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if set("".join(cells)) <= set("-: "):
            continue                       # the |---|---| separator row
        if not header:
            header = [_norm(c) for c in cells]
            continue
        for name, cell in zip(header, cells):
            if name in RANGES:
                m = re.fullmatch(r"\**\s*(" + _NUMBER + r")\s*\**", cell)
                if m:
                    out.append((name, m.group(1)))
    return out


def _norm(field: str) -> str:
    return field.strip().strip("*`").lower().replace("-", "_")


def check_value(field: str, value) -> "str | None":
    """Complaint string if this value is impossible for this field, else None."""
    spec = RANGES.get(_norm(field))
    if spec is None:
        return None
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v):
        return None
    low, high, what, likely = spec
    if low <= v <= high:
        return None
    bound = ("at most " + _fmt(high)) if v > high else ("at least " + _fmt(low))
    return (f"{field} = {_fmt(v)} is impossible: {what} must be {bound}. "
            f"This is usually {likely} under the wrong name — check the tool "
            f"output and use the field the tool actually returned.")


def scan_text(text: str, limit: int = 12) -> "list[str]":
    """Impossible `field = value` claims in a finished answer."""
    text = text or ""
    claims: "list[tuple[str, str]]" = list(_table_claims(text))
    # Table lines are handled above, by header. Running the prose scan over
    # them too reads across cell boundaries: in a header
    # `| gene | p_adj | neg_log10_pvalue |` it finds p_adj, runs on into the
    # next column, and reports `p_adj = 10` out of the digits in "log10".
    prose = "\n".join("" if ln.lstrip().startswith("|") else ln
                      for ln in text.splitlines())
    for field in RANGES:
        for m in re.finditer(r"(?:^|[\s`*|(\[,;])" + re.escape(field) +
                             r"\b(?:\*\*|`)?", prose,
                             re.IGNORECASE | re.MULTILINE):
            tail = _AFTER_FIELD.match(prose, m.end())
            if tail:
                claims.append((field, tail.group("value")))
    out, seen = [], set()
    for field, value in claims:
        key = (_norm(field), value)
        if key in seen:
            continue
        seen.add(key)
        complaint = check_value(field, value)
        if complaint:
            out.append(complaint)
            if len(out) >= limit:
                break
    return out


def _fmt(v: float) -> str:
    if v in (math.inf, -math.inf):
        return "infinite"
    return f"{v:g}"
