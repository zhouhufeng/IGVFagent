#!/usr/bin/env python3
"""Score an existing call table with the ported confusion_stats (used to check
the port of R/benchmarking.R::confusion_stats on the R reference's own calls).
Usage: score_calls.py <calls.csv cell,call> <truth.csv> <tag_mapping.json> <out.json>"""
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd

spec = importlib.util.spec_from_file_location(
    "d2b", Path(__file__).resolve().parents[3] / "Scripts/ported/skills/demultiplex2_benchmark.py")
if not Path(spec.origin).is_file():
    spec = importlib.util.spec_from_file_location(
        "d2b", Path(__file__).resolve().parents[3] / "Data/zhu2024/port/demultiplex2_benchmark.py")
d2b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d2b)

calls_p, truth_p, map_p, out_p = sys.argv[1:5]
c = pd.read_csv(calls_p, keep_default_na=False)
calls = pd.Series(c["call"].replace({"NA": None}).values, index=c["cell"].astype(str))
t = pd.read_csv(truth_p)
m = json.loads(Path(map_p).read_text())
cs = d2b.confusion_stats(calls, pd.Series(t["truth"].values, index=t["cell"].astype(str)),
                         pd.DataFrame({"tag": list(m), "true_label": list(m.values())}))
Path(out_p).write_text(json.dumps({**cs["singlet_avg_stats"], "doublet": cs["doublet_avg_stats"],
                                   "per_tag": cs["tag_stats"]}, indent=2))
print(out_p, cs["singlet_avg_stats"])
