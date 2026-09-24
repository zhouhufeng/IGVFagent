"""Every extension promoted into Scripts/promoted/ stays honest: its manifest
loads as a tool, its module imports and answers --help, the tool is registered
as a built-in, and no high-risk pattern appeared that the reviewer did not
accept (extension_review.py). Passes trivially while nothing is promoted."""
import contextlib
import importlib.util
import io
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _userext  # noqa: E402
import extension_review as er  # noqa: E402


def main() -> int:
    ok, n = True, 0
    reg_p = HERE / "promoted" / "registry.json"
    reg = json.loads(reg_p.read_text()) if reg_p.exists() else {}
    import _tools
    for path in sorted((HERE / "promoted" / "tools").glob("*.json")):
        n += 1
        spec = _userext._normalize_tool(_userext._load_manifest(path), path)
        good = bool(spec) and spec["name"] in _tools._BY_NAME and spec["name"] not in _tools._USER_TOOL_NAMES
        print(f"  {'ok  ' if good else 'FAIL'}  tool {path.stem}: loads and is a built-in")
        ok &= good
        prov = json.loads(path.read_text()).get("x-promoted") or {}
        good = bool(prov.get("reviewer")) and path.stem in reg
        print(f"  {'ok  ' if good else 'FAIL'}  tool {path.stem}: provenance recorded")
        ok &= good
    for path in sorted((HERE / "promoted" / "skills").glob("*.py")):
        n += 1
        spec = importlib.util.spec_from_file_location(f"promoted_{path.stem}", path)
        mod = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(mod)
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                try:
                    mod.main(["--help"])
                    code = 0
                except SystemExit as e:
                    code = e.code or 0
            good = code == 0
        except Exception as e:  # noqa: BLE001
            print(f"        {type(e).__name__}: {e}")
            good = False
        print(f"  {'ok  ' if good else 'FAIL'}  skill {path.stem}: imports and answers --help")
        ok &= good
        entry = next((v for v in reg.values() if any(f.endswith(path.name) for f in v.get("files", []))), {})
        new_high = {f["flag"] for f in er.scan(path.read_text()) if f["level"] == "high"} - set(
            entry.get("flags_at_review") or [])
        good = not new_high and (not entry.get("flags_at_review") or entry.get("accepted_risk") or not {
            f["flag"] for f in er.scan(path.read_text()) if f["level"] == "high"})
        print(f"  {'ok  ' if good else 'FAIL'}  skill {path.stem}: no unreviewed high-risk pattern"
              + (f" ({', '.join(new_high)})" if new_high else ""))
        ok &= bool(good)
    print(f"{n} promoted item(s)")
    print("all checks pass" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
