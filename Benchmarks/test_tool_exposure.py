#!/usr/bin/env python3
"""A capability the CLI has but the agent's schema does not is invisible.

`igvfagent bean analyze` grew a `--run-bean` flag that drives the real
crispr-bean end to end. The agent-facing declaration was never updated, so
the model had no parameter to set: asked for BEAN, it ran the frequentist
path and then correctly reported that no BEAN output existed and that "no
crispr_bean tool appears to be available". It was right, and the cause was
this gap, not the installation.

The general failure is worth guarding: a flag that exists in the CLI and
matters to a user request must be reachable from the tool schema, and a bool
flag must actually reach argv.
"""
import inspect
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Scripts"))
import _tools  # noqa: E402

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:64} {detail}")
    if not ok:
        FAILURES.append(name)


def tool(name):
    hits = [t for t in _tools._TOOLS if t.name == name]
    return hits[0] if hits else None


# ── the specific regression ───────────────────────────────────────────────
t = tool("base_editing_screen_analyze")
check("the analyze tool is registered", t is not None)
props = t.parameters["properties"]
check("run_bean is a declared parameter", "run_bean" in props)
check("it is a boolean", props.get("run_bean", {}).get("type") == "boolean")
check("it is declared as a bare flag", "run_bean" in t.bool_flags)
argv = _tools._build_argv(t, {"accession": "IGVFDS6464SOVZ", "run_bean": True})
check("run_bean=true emits --run-bean", "--run-bean" in argv, " ".join(argv))
argv_off = _tools._build_argv(t, {"accession": "IGVFDS6464SOVZ",
                                   "run_bean": False})
check("run_bean=false emits nothing", "--run-bean" not in argv_off)
check("omitting it emits nothing",
      "--run-bean" not in _tools._build_argv(t, {"accession": "X"}))
argv_full = _tools._build_argv(t, {"accession": "X", "run_bean": True,
                                    "bean_iter": 300, "bean_mode": "variant"})
check("bean_iter maps to --bean-iter", "--bean-iter" in argv_full)
check("bean_mode maps to --bean-mode", "--bean-mode" in argv_full)
check("the accession stays positional, not a flag",
      argv_full[3] == "X", " ".join(argv_full))

# The description is the only thing that makes the model SET the flag.
d = t.description
check("the description tells the model to set run_bean for BEAN",
      "run_bean=true" in d)
check("and says what happens without it",
      "frequentist" in d.lower() or "without run_bean" in d.lower())
check("and does not still claim it merely prints the command",
      "prints the real BEAN command" not in d)

# ── the general rule: declared flags must exist in the CLI ────────────────
# A flag_map entry pointing at an option the CLI does not accept fails with
# argparse exit 2 at call time, which reads to the model as a broken tool.
def cli_options(module, subcommand):
    """Options the CLI really accepts, read from the skill's own parser.

    cli.py cannot be run as a script (it uses a package-relative import), so
    an earlier version of this helper always failed and printed "cross-check
    skipped" -- a guard that never runs is worse than no guard, because it
    reads as a pass. The skill module is imported and its parser inspected
    instead, which is also faster than spawning a process.
    """
    import importlib
    mod = importlib.import_module(module)
    parser = mod.build_parser()
    subs = [a for a in parser._actions
            if hasattr(a, "choices") and isinstance(a.choices, dict)]
    if not subs or subcommand not in subs[0].choices:
        return None
    sp = subs[0].choices[subcommand]
    out = set()
    for a in sp._actions:
        out.update(a.option_strings)
    return out or None


opts = cli_options("base_editing_screen", "analyze")
if opts:
    for param, flag in t.flag_map.items():
        check(f"--{param}'s flag {flag} exists in the CLI", flag in opts)
    for b in t.bool_flags:
        flag = t.flag_map.get(b, "--" + b.replace("_", "-"))
        check(f"bool flag {flag} exists in the CLI", flag in opts)
else:
    print("note: could not read `bean analyze --help`; CLI cross-check skipped")

# ── no tool may declare a parameter it cannot pass ────────────────────────
bad = []
for tt in _tools._TOOLS:
    declared = set(tt.parameters.get("properties", {}))
    passable = (set(tt.positional) | set(tt.flag_map)
                | set(tt.bool_flags) | set(tt.flag_repeat))
    # Anything declared but not positional and not in a flag map still gets a
    # --dashed-name by convention, so this only catches genuine orphans in
    # tools that use an explicit flag_map for everything else.
    orphan = declared - passable
    if orphan and tt.flag_map and len(orphan) > len(declared) / 2:
        bad.append((tt.name, sorted(orphan)))
check("no tool declares mostly-unpassable parameters", not bad, str(bad[:2]))

# ── a capability needs a NAME, not just a flag ───────────────────────────
# Exposing run_bean was not enough. Asked for crispr-bean, the agent still
# called base_editing_screen_analyze without it and reported "no
# crispr_bean_sorting or similar tool appears to be available" -- it was
# searching for a tool NAMED after BEAN. A flag on a differently-named tool
# is not discoverable by a model that searches the registry by name, which
# is how anyone would look.
bt = tool("crispr_bean_analyze")
check("a tool named after BEAN exists", bt is not None)
check("its name contains 'bean' so a name search finds it",
      "bean" in (bt.name if bt else ""))
bargv = _tools._build_argv(bt, {"accession": "IGVFDS6464SOVZ"})
check("it routes to the bayesian subcommand",
      bargv[1:3] == ["bean", "bayesian"], " ".join(bargv))
check("the caller does not have to remember a flag",
      "--run-bean" not in bargv and "run_bean" not in bt.parameters["properties"])
d2 = bt.description
for phrase in ("crispr-bean", "posterior", "credible interval"):
    check(f"the description names {phrase!r}", phrase in d2.lower())
check("it steers away from the frequentist tool",
      "Do NOT use base_editing_screen_analyze" in d2)
check("it pre-empts the NOT activity-normalised surprise",
      "NOT activity-normalised" in d2)
# The subcommand must exist and force the flag on, or the tool is a no-op.
import base_editing_screen as _bes  # noqa: E402
_a = _bes.build_parser().parse_args(["bayesian", "IGVFDS6464SOVZ"])
check("`bean bayesian` forces run_bean on", _a.run_bean is True)
check("`bean analyze` still does not", 
      _bes.build_parser().parse_args(["analyze", "X"]).run_bean is False)
check("both dispatch to the same implementation",
      "bayesian" in inspect.getsource(_bes.main))
_opts = cli_options("base_editing_screen", "bayesian")
if _opts:
    for _p, _f in bt.flag_map.items():
        check(f"{_f} exists on the bayesian subcommand", _f in _opts)

print(f"\n{len(FAILURES)} failure(s)")
sys.exit(1 if FAILURES else 0)
