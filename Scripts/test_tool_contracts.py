"""The tool catalogue's own contracts, asserted.

These are the failures that cost a whole agent iteration and surface to the
user as "3 tool call(s) did not succeed" with no clue why. A retest on
2026-09-15 lost an iteration to exactly one of them: the schema declared
`verbose` as a boolean, the tool did not list it among the bare flags, so the
wrapper emitted `--verbose true` and argparse exited 2. Eleven other tools had
the same defect waiting.

    python3 Scripts/test_tool_contracts.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _tools  # noqa: E402

FAILED: "list[str]" = []


def fail(msg: str) -> None:
    print(f"  FAIL  {msg}")
    FAILED.append(msg)


def main() -> int:
    tools = _tools._TOOLS
    print(f"\n{len(tools)} tools in the catalogue")

    print("\nevery boolean parameter is a BARE flag")
    n = 0
    for t in tools:
        for name, spec in ((t.parameters or {}).get("properties") or {}).items():
            if not isinstance(spec, dict) or spec.get("type") != "boolean":
                continue
            n += 1
            if name not in t.bool_flags:
                fail(f"{t.name}.{name} is boolean but not in bool_flags — "
                     f"the wrapper would emit `--{name.replace('_','-')} true` "
                     f"and argparse exits 2")
    print(f"  checked {n} boolean parameter(s)")

    print("\nbare flags emit no value, and are omitted when false")
    for t in tools:
        for name in t.bool_flags:
            argv_true = _tools._build_argv(t, {name: True})
            flag = t.flag_map.get(name, "--" + name.replace("_", "-"))
            if flag not in argv_true:
                fail(f"{t.name}.{name}: true did not emit {flag}")
            elif argv_true[argv_true.index(flag) + 1:][:1] in (["true"], ["True"]):
                fail(f"{t.name}.{name}: emitted a value after {flag}")
            if flag in _tools._build_argv(t, {name: False}):
                fail(f"{t.name}.{name}: false still emitted {flag}")

    print("\nevery required parameter is declared in properties")
    for t in tools:
        props = (t.parameters or {}).get("properties") or {}
        for req in (t.parameters or {}).get("required") or []:
            if req not in props:
                fail(f"{t.name}: required parameter {req!r} has no schema")

    print("\nevery parameter is routed — positional, mapped, or defaulted")
    for t in tools:
        props = (t.parameters or {}).get("properties") or {}
        for name in props:
            if name in t.positional or name in t.flag_map:
                continue
            # falls back to --name-with-dashes, which is legitimate
            argv = _tools._build_argv(t, {name: "x"})
            if "--" + name.replace("_", "-") not in argv and name not in t.bool_flags:
                fail(f"{t.name}.{name}: not positional, not mapped, and no "
                     f"default flag emitted")

    print("\ntool names are unique")
    seen = set()
    for t in tools:
        if t.name in seen:
            fail(f"duplicate tool name {t.name!r}")
        seen.add(t.name)

    print()
    if FAILED:
        print(f"{len(FAILED)} contract violation(s)")
        return 1
    print("all contracts hold")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
