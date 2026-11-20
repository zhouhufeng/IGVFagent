"""Offline checks for changing existing Portal objects and reading staging ones.

    python3 Scripts/test_portal_patch.py

* `portal get` resolves @ids, accessions, UUIDs and Portal URLs to a Portal and a path (a TST
  accession pasted onto the prod host was a DNS error);
* the agent's `igvf_submit_sheet_patch` and `portal_get` tools build the commands their CLIs take;
* a clean `submit sheet-patch` dry run saves a plan under Data/Submitter, and only a plan that holds
  its own checked sheet, for the same Portal, turns into an `--apply` command.

No network: the profile the dry run validates against is a stub.
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

FAIL = []


def check(name: str, cond: bool) -> None:
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}")
    if not cond:
        FAIL.append(name)


def raises(fn) -> bool:
    try:
        fn()
    except ValueError:
        return True
    return False


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="portalpatch_")).resolve()
    os.environ["IGVF_PROJECT_ROOT"] = str(tmp)
    for k in ("IGVF_STAGING_API_KEY", "IGVF_STAGING_SECRET_KEY", "IGVF_API_KEY", "IGVF_SECRET_KEY"):
        os.environ.pop(k, None)
    os.environ["IGVF_NO_KEY_FILES"] = "1"

    print("portal get: which Portal, which path")
    import portal_query_skill as pq
    r = pq.resolve_ref
    check("a TST accession is a staging object", r("TSTFI48600603") == ("staging", "/TSTFI48600603/"))
    check("an IGVF accession is prod", r("IGVFFI1234ABCD") == ("prod", "/IGVFFI1234ABCD/"))
    check("an @id keeps its path; TST still means staging",
          r("/tabular-files/TSTFI48600603/") == ("staging", "/tabular-files/TSTFI48600603/"))
    check("a staging URL names staging", r("https://api.staging.igvf.org/tabular-files/TSTFI48600603/")
          == ("staging", "/tabular-files/TSTFI48600603/"))
    check("a UI URL names its Portal too", r("https://data.igvf.org/files/IGVFFI1234ABCD/") == ("prod", "/files/IGVFFI1234ABCD/"))
    check("mode given explicitly wins over the accession", r("IGVFDS3909HJKS", "staging") == ("staging", "/IGVFDS3909HJKS/"))
    check("a UUID is a path of its own", r("6a9e7f1c-3b2d-4c5e-8f9a-0b1c2d3e4f5a")[1] == "/6a9e7f1c-3b2d-4c5e-8f9a-0b1c2d3e4f5a/")
    check("every path starts and ends with /", all(r(x)[1].startswith("/") and r(x)[1].endswith("/")
                                                  for x in ("TSTFI1", "files/IGVFFI1", "/x")))
    check("a URL on another host is refused", raises(lambda: r("https://example.org/files/IGVFFI1/")))
    check("a staging URL with mode=prod is refused", raises(lambda: r("https://api.staging.igvf.org/x/", "prod")))
    check("an empty id is refused", raises(lambda: r(" / ")))

    print("the agent's tools")
    import _tools as T
    tools = {t.name: t for t in T.list_tools()}
    check("igvf_submit_sheet_patch is registered", "igvf_submit_sheet_patch" in tools)
    check("its name gives it the user's own Portal keys", "igvf_submit_sheet_patch".startswith("igvf_submit_"))
    check("it cannot send: no apply parameter",
          "apply" not in tools["igvf_submit_sheet_patch"].parameters.get("properties", {}))
    check("sheet-patch argv: sheet first, --columns, mode last",
          T._build_argv(tools["igvf_submit_sheet_patch"], {"sheet": "tabular_file.tsv", "columns": "description",
                                                           "all": False, "mode": "staging"})
          == ["igvfagent", "submit", "sheet-patch", "tabular_file.tsv", "--columns", "description", "-m", "staging"])
    check("sheet-patch argv: all is a bare flag",
          T._build_argv(tools["igvf_submit_sheet_patch"], {"sheet": "x.tsv", "all": True})[-1] == "--all")
    check("portal_get passes mode",
          T._build_argv(tools["portal_get"], {"resource_id": "TSTFI48600603", "mode": "staging"})
          == ["igvfagent", "portal", "get", "TSTFI48600603", "--mode", "staging"])

    print("sheet-patch dry run saves a plan the Submit button can send")
    import igvf_metadata_sheet as ms
    import igvf_submission_skill as ss
    sheet = tmp / "tabular_file.tsv"
    sheet.write_text("accession\tdescription\nTSTFI48600603\tA better description\n")
    ms._load_profile = lambda portal, name: {
        "title": "Tabular file", "identifyingProperties": ["uuid", "accession", "aliases"],
        "properties": {"accession": {"type": "string"}, "description": {"type": "string"},
                       "uuid": {"type": "string"}, "aliases": {"type": "array", "items": {"type": "string"}}}}
    args = argparse.Namespace(sheet=str(sheet), sheet_name=None, profile=None, columns="description", all=False,
                              remove=None, apply=False, yes=False, out=None, admin=False, mode="staging")
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = ms._submit(ss.Portal(mode="staging"), args, "PATCH")
    text = out.getvalue()
    check("the dry run passes and sends nothing", rc == 0 and "would send 1 property" in text and "nothing sent" in text)
    line = next((ln for ln in text.splitlines() if ln.startswith("Patch plan: ")), "")
    check("it prints the line the app looks for", line.endswith("(1 object(s), staging)"))
    plan_path = Path(line.split()[2]) if line else tmp / "missing"
    check("the plan is under Data/Submitter", tmp / "Data" / "Submitter" in plan_path.parents)
    plan = json.loads(plan_path.read_text()) if plan_path.is_file() else {}
    check("the plan holds a copy of the sheet as checked", Path(plan.get("sheet", "")).parent == plan_path.parent
          and Path(plan.get("sheet", "")).read_text() == sheet.read_text())
    cmd = ms.patch_plan_command(plan_path, "staging")
    check("a good plan becomes the --apply command", isinstance(cmd, list) and cmd[:3] == ["submit", "sheet-patch", plan["sheet"]]
          and "--apply" in cmd and cmd[cmd.index("--columns") + 1] == "description" and "--yes" not in cmd)
    check("a plan for staging is not sent to prod", isinstance(ms.patch_plan_command(plan_path, "prod"), str))
    bad = dict(plan, sheet=str(sheet))                          # the user's sheet, which may have changed since
    plan_path.write_text(json.dumps(bad))
    check("a plan whose sheet is outside its folder is refused", isinstance(ms.patch_plan_command(plan_path, "staging"), str))
    check("prod plans carry --yes", "--yes" in ms.patch_plan_argv(dict(plan, mode="prod")))
    sheet.write_text("accession\tdescription\nTSTFI48600603\t\n")
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        ms._submit(ss.Portal(mode="staging"), args, "PATCH")
    check("nothing to send: no plan", "Patch plan:" not in out.getvalue())

    # ── saved plans are offered to the person who made them, however the dry run ran ──
    meta = json.loads((plan_path.parent / ss.PLAN_META).read_text()) if (plan_path.parent / ss.PLAN_META).is_file() else {}
    check("a patch dry run saves plan.meta.json (kind, count, mode)", meta.get("kind") == "patch" and meta.get("n") == 1
          and meta.get("mode") == "staging" and meta.get("plan") == "patch.json")
    root = tmp / "offer"
    sub = root / "Data" / "Submitter"

    def plan(name: str, owner: str, age_h: float = 0.0, kind: str = "create") -> Path:
        d = sub / name
        d.mkdir(parents=True, exist_ok=True)
        f = d / ("drafts.json" if kind == "create" else "patch.json")
        f.write_text("[]")
        old = os.environ.get("IGVF_ACTING_USER")
        os.environ["IGVF_ACTING_USER"] = owner
        try:
            mp = ss.save_plan_meta(f, kind, 2, "staging")
        finally:
            os.environ.pop("IGVF_ACTING_USER") if old is None else os.environ.__setitem__("IGVF_ACTING_USER", old)
        if age_h:
            ss.mark_plan(mp, created=time.time() - age_h * 3600)
        return mp

    (sub / "drafts").mkdir(parents=True)                       # folders without a meta file are ignored
    plan("20261120_170436_plan_staging", "HZHOU")              # a deep job's plan
    m2 = plan("20261120_174439_plan_staging", "alice")
    plan("20261120_180000_plan_staging", "bob")
    got = ss.latest_unsent_plan(root, "alice")
    check("a person is offered their own latest plan", got is not None and got["plan"].endswith("174439_plan_staging/drafts.json")
          and got["n"] == 2 and got["mode"] == "staging" and got["kind"] == "create")
    check("a job owner is offered the plan their job made", (ss.latest_unsent_plan(root, "HZHOU") or {}).get("plan", "")
          .endswith("170436_plan_staging/drafts.json"))
    check("nobody is offered someone else's plan", ss.latest_unsent_plan(root, "carol") is None
          and ss.latest_unsent_plan(root, "") is None)
    ss.mark_plan(m2, last_exit=1)
    check("a plan whose sending stopped is offered again", ss.latest_unsent_plan(root, "alice") is not None)
    ss.mark_plan(m2, applied_at="2026-11-20T18:00:00")
    check("a sent plan is not offered again", ss.latest_unsent_plan(root, "alice") is None)
    plan("20261120_181000_plan_staging", "alice")
    m4 = plan("20261120_182000_patch_staging", "alice", kind="patch")
    check("the newest plan wins, a patch plan included", (ss.latest_unsent_plan(root, "alice") or {}).get("kind") == "patch")
    ss.mark_plan(m4, discarded_at="2026-11-20T18:30:00")
    check("a discarded newest plan also sets the older ones aside", ss.latest_unsent_plan(root, "alice") is None)
    plan("20261118_090000_plan_staging", "dave", age_h=60)
    check("a plan older than two days is not offered", ss.latest_unsent_plan(root, "dave") is None)

    # `submit send`: the command-line way to send a saved plan (what the web app's button does)
    def stops(fn) -> bool:
        try:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                fn()
        except SystemExit:
            return True
        return False

    def send_args(plan_file, mode=None, yes=False):
        return argparse.Namespace(plan=str(plan_file), mode=mode, yes=yes)

    calls: "list[list[str]]" = []
    real_main = ss.main
    ss.main = lambda argv=None: (calls.append(list(argv)), 0)[1]
    try:
        cp = sub / "20261120_183000_plan_staging" / "drafts.json"
        cp.parent.mkdir(parents=True)
        cp.write_text("[]")
        cm = ss.save_plan_meta(cp, "create", 2, "staging")
        check("send: a missing plan is refused", stops(lambda: ss.cmd_send(None, send_args(sub / "nope.json"))))
        check("send: a plan is sent only to the Portal it was checked on",
              stops(lambda: ss.cmd_send(None, send_args(cp, mode="prod", yes=True))))
        with contextlib.redirect_stdout(io.StringIO()):
            rc = ss.cmd_send(None, send_args(cp))
        check("send: a create plan becomes create --apply --upload on its own Portal",
              rc == 0 and calls[-1] == ["-m", "staging", "create", "--objects", str(cp), "--apply", "--upload"])
        check("send: a sent plan is marked, so the web app stops offering it",
              bool(json.loads(cm.read_text()).get("applied_at")) and ss.latest_unsent_plan(root, "") is None)
        pp = sub / "20261120_184000_plan_prod" / "drafts.json"
        pp.parent.mkdir(parents=True)
        pp.write_text("[]")
        ss.save_plan_meta(pp, "create", 1, "prod")
        check("send: prod needs --yes", stops(lambda: ss.cmd_send(None, send_args(pp))))
        with contextlib.redirect_stdout(io.StringIO()):
            ss.cmd_send(None, send_args(pp, yes=True))
        check("send: and passes it on", calls[-1][-1] == "--yes" and calls[-1][:2] == ["-m", "prod"])
        with contextlib.redirect_stdout(io.StringIO()):
            ms_plan = json.loads(plan_path.read_text()) if plan_path.is_file() else {}
            good = dict(ms_plan, sheet=str(plan_path.parent / Path(ms_plan.get("sheet", "x")).name))
            plan_path.write_text(json.dumps(good))
            rc = ss.cmd_send(None, send_args(plan_path))
        check("send: a patch plan becomes its saved sheet-patch command",
              rc == 0 and calls[-1][:2] == ["sheet-patch", good.get("sheet")] and "--apply" in calls[-1])
    finally:
        ss.main = real_main

    print(f"\n{'FAILED: ' + ', '.join(FAIL) if FAIL else 'all checks passed'}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
