#!/usr/bin/env python3
"""IGVFagent's scores next to the upstream CRISPR_comparison pipeline's, predictor by predictor.

Reads every ``results/upstream/<comparison>/performance_summary.txt`` (written
by EngreitzLab/CRISPR_comparison, run on the same inputs; see README) and the
latest merged IGVFagent run under Docs/scE2G/*_k562_crowdsourced_feature_benchmark,
and writes ``results/upstream_vs_igvfagent.tsv``.

The upstream AUPRC is caTools::trapz over yardstick's tie-aware precision-recall
points with the first (recall 0) and last rows dropped; IGVFagent reports that as
``auprc_crispr_comparison`` (column igvfagent_auprc_cc here) next to the step-rule
``auprc``. Precision at 70% recall is defined identically on both sides.
"""
import csv
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
UP = os.path.join(HERE, "results", "upstream")

# upstream comparison -> which IGVFagent predictor ids it corresponds to
MAP = {
    "re2g_anchor": {"rE2G_base.Full.Score": "rE2G_base", "rE2G_ext.FullModel.Score": "rE2G_ext",
                    "rE2G_ext_noEP300.FullModel_minus_EP300.Score": "rE2G_ext_noEP300"},
    "sce2g_igvf_notssfilter": {"scE2G.Score": "scE2G", "scE2G_ignoreTPM.Score.ignoreTPM": "scE2G_ignoreTPM",
                               "ABC.ABC.Score": "ABC", "ARC_E2G.ARC.E2G.Score": "ARC_E2G"},
    "features_igvf_notssfilter": {"Pinloop.Pinloop": "Pinloop:Pinloop", "Signac.Signac_Score": "Signac:Signac_Score",
                                  "SCENT.SCENT_beta": "SCENT:SCENT_beta"},
    "sce2g_igvf": {"scE2G.Score": "scE2G", "scE2G_ignoreTPM.Score.ignoreTPM": "scE2G_ignoreTPM",
                   "ABC.ABC.Score": "ABC", "ARC_E2G.ARC.E2G.Score": "ARC_E2G"},
    "features_igvf": {"Pinloop.Pinloop": "Pinloop:Pinloop", "Signac.Signac_Score": "Signac:Signac_Score",
                      "SCENT.SCENT_beta": "SCENT:SCENT_beta"},
}
TSS_FILTERED = {"sce2g_igvf", "features_igvf", "re2g_anchor"}   # upstream default filter_pred_tss = True


def main() -> int:
    merged = sorted(glob.glob(os.path.join(ROOT, "Docs", "scE2G", "*_k562_crowdsourced_feature_benchmark")))
    if not merged:
        sys.exit("no merged IGVFagent run found")
    summ = json.load(open(os.path.join(merged[-1], "summary.json")))
    preds = summ["predictors"]

    def mine(pid):
        return preds.get(pid.replace(".", "_"))

    rows = []
    for comp, ids in MAP.items():
        f = os.path.join(UP, comp, "performance_summary.txt")
        if not os.path.isfile(f):
            continue
        for r in csv.DictReader(open(f), delimiter="\t"):
            pid = ids.get(r["pred_uid"])
            if pid is None:
                if r["pred_uid"].startswith("baseline."):
                    rows.append({"upstream_comparison": comp, "predictor": r["pred_uid"], "tss_filter": comp in TSS_FILTERED,
                                 "upstream_auprc": round(float(r["AUPRC"]), 4), "upstream_prec70": round(float(r["PrecMinSens"]), 4),
                                 "igvfagent_auprc_cc": "", "igvfagent_auprc_step": "", "igvfagent_prec70": "",
                                 "delta_auprc_cc": "", "note": "upstream-only baseline (computed per pair, no overlap)"})
                continue
            m = mine(pid)
            up_a, up_p = float(r["AUPRC"]), float(r["PrecMinSens"])
            rows.append({"upstream_comparison": comp, "predictor": pid, "tss_filter": comp in TSS_FILTERED,
                         "upstream_auprc": round(up_a, 4), "upstream_prec70": round(up_p, 4),
                         "igvfagent_auprc_cc": round(m["auprc_crispr_comparison"], 4) if m else "",
                         "igvfagent_auprc_step": round(m["auprc"], 4) if m else "",
                         "igvfagent_prec70": round(m["precision_at_70_recall"], 4) if m else "",
                         "delta_auprc_cc": round(m["auprc_crispr_comparison"] - up_a, 4) if m else "",
                         "note": "" if m else "not in IGVFagent merged run"})
    out = os.path.join(HERE, "results", "upstream_vs_igvfagent.tsv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['upstream_comparison']:26s} {r['predictor']:22s} tssfilter={str(r['tss_filter']):5s} "
              f"upstream {r['upstream_auprc']:.4f}  igvfagent(cc)   {r['igvfagent_auprc_cc'] or 'nan':>7} "
              f"delta {r['delta_auprc_cc'] if r['delta_auprc_cc'] != '' else 'nan':>8}  "
              f"P@70 {r['upstream_prec70']:.4f} vs {r['igvfagent_prec70'] or 'nan'}")
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
