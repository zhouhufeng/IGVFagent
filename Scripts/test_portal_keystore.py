"""Each web user's Portal keys: one slot per host, both kept, never mixed up.

Offline: the Portal hosts are replaced by a fake that knows one staging pair
and one prod pair.

    python3 Scripts/test_portal_keystore.py
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

_TMP = Path(tempfile.mkdtemp(prefix="igvf-keys-"))
os.environ["IGVF_PROJECT_ROOT"] = str(_TMP)
os.environ["IGVF_KEYSTORE_SECRET"] = "test-master-secret"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import portal_keystore as K  # noqa: E402

FAILED: "list[str]" = []
STAGING = ("ABCD1234", "abcd1234efgh5678")
PROD = ("WXYZ9876", "wxyz9876stuv5432")


def check(name: str, got, want) -> None:
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + ("" if ok else f"   got {got!r}, want {want!r}"))
    if not ok:
        FAILED.append(name)


def fake_verify(key: str, secret: str, mode: str, timeout: float = 20):
    if (mode, (key, secret)) in (("staging", STAGING), ("prod", PROD)):
        return {"account": "Ada Lovelace", "submits_for": ["xihong-lin"]}
    return None


def main() -> int:
    K.verify = fake_verify
    user = "ada"

    print("typed into the keys panel")
    ok, msg = K.register(user, "staging", "abc", "def")
    check("a malformed pair is refused before any request", (ok, "does not look like" in msg), (False, True))
    ok, msg = K.register(user, "prod", *STAGING)
    check("a staging pair offered as prod is refused", ok, False)
    check("...and the message says where it works", "works on staging" in msg, True)
    check("staging pair saved", K.register(user, "staging", *STAGING)[0], True)
    check("prod pair saved", K.register(user, "prod", *PROD)[0], True)
    check("both slots kept side by side", sorted(h["mode"] for h in K.listing(user)), ["prod", "staging"])
    check("each host gets its own pair", (K.load(user, "staging"), K.load(user, "prod")), (STAGING, PROD))
    env = K.env_for(user)
    check("submission tools get both, per host",
          (env["IGVF_STAGING_API_KEY"], env["IGVF_PROD_API_KEY"]), (STAGING[0], PROD[0]))
    check("the listing never shows a secret", any(STAGING[1] in str(h) or PROD[1] in str(h)
                                                  for h in K.listing(user)), False)
    check("secrets are encrypted at rest", STAGING[1] in K._path(user).read_text(), False)

    print("\npasted into the chat")
    clean, notes = K.intercept(f"here are my igvf keys {PROD[0]} {PROD[1]}", "grace")
    check("the pair is cut out of the message", PROD[1] in clean, False)
    check("...and saved under the host that accepted it", [h["mode"] for h in K.listing("grace")], ["prod"])

    print("\nreplacing and deleting one host leaves the other alone")
    K.delete(user, "staging")
    check("staging deleted, prod kept", [h["mode"] for h in K.listing(user)], ["prod"])
    check("staging saved again", K.register(user, "staging", *STAGING)[0], True)
    check("...prod unchanged", K.load(user, "prod"), PROD)

    print()
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED: {', '.join(FAILED)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
