"""Offline checks for absorbing a paper's code (port_skill.py) and for judging a
reproduction by paper coverage rather than by a check tally: registration and
its refusals, `port find`, `bench verify-port`'s comparison, and
concordance.judge_coverage. Writes only under a temporary directory."""
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "Benchmarks"))
import benchmark_skill as bs  # noqa: E402
import concordance  # noqa: E402
import port_skill  # noqa: E402

PORT_SRC = '''"""Mean of a numeric column."""
import argparse


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--values", nargs="+", type=float, required=True)
    a = ap.parse_args(argv)
    print(sum(a.values) / len(a.values))
    return 0
'''


def main() -> int:
    ok = True

    def check(label, cond):
        nonlocal ok
        print(f"  {'ok  ' if cond else 'FAIL'}  {label}")
        ok &= bool(cond)

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        port_skill.PORTED = td / "ported"
        port_skill.REGISTRY = port_skill.PORTED / "registry.json"
        src = td / "mean_port.py"
        src.write_text(PORT_SRC)
        args = dict(name="column_mean_demo", paper_id="demo2026", repo="https://github.com/x/y",
                    commit="abc1234", upstream=["R/mean.R"], source_file=str(src),
                    description="Mean of a column, ported from R/mean.R",
                    keywords="average,summary")

        os.environ.pop("IGVF_ALLOW_AGENT_AUTHORING", None)
        try:
            port_skill.register(**args)
            check("register refuses without IGVF_ALLOW_AGENT_AUTHORING", False)
        except port_skill._ea.AuthoringDisabled:
            check("register refuses without IGVF_ALLOW_AGENT_AUTHORING", True)

        os.environ["IGVF_ALLOW_AGENT_AUTHORING"] = "1"
        for bad, why in ((dict(upstream=[]), "no upstream source"),
                         (dict(commit="main"), "unpinned commit"),
                         (dict(name="bench"), "a built-in name")):
            try:
                port_skill.register(**{**args, **bad})
                check(f"register refuses {why}", False)
            except port_skill._ea.InvalidExtension:
                check(f"register refuses {why}", True)

        out = port_skill.register(**args)
        skill = port_skill.PORTED / "skills" / "column_mean_demo.py"
        check("register writes the module and manifest",
              skill.is_file() and (port_skill.PORTED / "tools" / "column_mean_demo.json").is_file())
        reg = port_skill._load_registry()["column_mean_demo"]
        check("provenance recorded, unreviewed", reg["commit"] == "abc1234"
              and reg["upstream"] == ["R/mean.R"] and reg["reviewed"] is False)
        check("module header names the source", "R/mean.R" in skill.read_text())
        check("command name", out["command"] == "igvfagent column-mean-demo")
        check("find matches by keyword", [h["name"] for h in port_skill.find("column average")]
              == ["column_mean_demo"])
        check("find ignores unrelated queries", port_skill.find("chromatin loops") == [])
        port_skill.remove("column_mean_demo")
        check("remove deletes files and registry entry",
              not skill.exists() and "column_mean_demo" not in port_skill._load_registry())
        os.environ.pop("IGVF_ALLOW_AGENT_AUTHORING", None)

    ref = {"a": 1.0, "b": 2.0, "c": 3.0, "d": 4.0}
    same = bs.compare_values(ref, dict(ref), rtol=1e-6, atol=1e-8)
    check("identical outputs match fully", same["match_rate"] == 1.0 and same["max_abs_diff"] == 0)
    half = bs.compare_values(ref, {"a": 1.0, "b": 2.0}, rtol=1e-6, atol=1e-8)
    check("missing port values count as mismatches", half["match_rate"] == 0.5
          and half["n_missing_in_port"] == 2)
    off = bs.compare_values(ref, {**ref, "d": 5.0}, rtol=1e-6, atol=1e-8)
    check("a wrong value is caught", off["match_rate"] == 0.75 and off["worst"][0]["key"] == "d")

    spec = {"analyses": [{"id": "fig2"}, {"id": "fig3"},
                         {"id": "fig4", "blocker": {"kind": "controlled_access", "reason": "dbGaP"}}]}
    strong = {"passed": True, "class": "B", "analysis": "fig2"}
    count = {"passed": True, "class": "C", "analysis": "fig3"}
    v = concordance.judge_coverage(spec, [strong, count])
    check("a count check does not cover an analysis",
          v["reproduction"] == "incomplete" and v["coverage"] == "1/3")
    v = concordance.judge_coverage(spec, [strong, {**strong, "analysis": "fig3"}])
    check("access-blocked analyses allow reproduced_except_access",
          v["reproduction"] == "reproduced_except_access")
    other = {"analyses": [{"id": "fig2"}, {"id": "fig5", "blocker": {"kind": "other", "reason": "slow"}}]}
    check("a non-access blocker keeps the paper incomplete",
          concordance.judge_coverage(other, [strong])["reproduction"] == "incomplete")
    check("no plan means unplanned, never reproduced",
          concordance.judge_coverage({}, [strong])["reproduction"] == "unplanned")

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
