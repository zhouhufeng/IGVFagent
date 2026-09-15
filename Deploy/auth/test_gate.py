"""End-to-end test of the authentication gate, over real HTTP.

Runs the gate on a loopback port and plays the part of both the browser and
Discourse. Worth having as a file rather than a one-off check: this is the
component that decides who may drive an agent that executes analysis CLIs, and
every case below is one where getting it wrong means either locking everyone
out or letting everyone in.

    python3 Deploy/auth/test_gate.py
"""
from __future__ import annotations

import base64
import http.client
import os
import sys
import threading
import time
import urllib.parse
from pathlib import Path

os.environ.setdefault("IGVF_SSO_SECRET", "sso-secret-for-tests")
os.environ.setdefault("IGVF_SESSION_SECRET", "session-secret-for-tests")
os.environ.setdefault("IGVF_AUTH_PORT", "9911")
os.environ.setdefault("IGVF_PUBLIC_URL", "https://igvfagent.example.org")
os.environ.setdefault("IGVF_DISCOURSE_URL", "https://discussion.example.org")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gate  # noqa: E402

from http.server import ThreadingHTTPServer  # noqa: E402

PORT = int(os.environ["IGVF_AUTH_PORT"])
FAILED: "list[str]" = []


def check(name: str, got, want) -> None:
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {name}"
          + ("" if ok else f"   got {got!r}, want {want!r}"))
    if not ok:
        FAILED.append(name)


def request(path: str, cookie: str = "") -> "tuple[int, dict, bytes]":
    conn = http.client.HTTPConnection("127.0.0.1", PORT, timeout=5)
    headers = {"Cookie": f"{gate.COOKIE}={cookie}"} if cookie else {}
    conn.request("GET", path, headers=headers)
    resp = conn.getresponse()
    body = resp.read()
    out = (resp.status, dict(resp.getheaders()), body)
    conn.close()
    return out


def discourse_callback(nonce: str, *, groups: str = "igvfagent",
                       username: str = "tester", admin: str = "false",
                       secret: "bytes | None" = None) -> str:
    """Build the URL Discourse would send the browser back to."""
    payload = urllib.parse.urlencode({
        "nonce": nonce, "email": f"{username}@example.org",
        "username": username, "name": username.title(),
        "external_id": "42", "groups": groups, "admin": admin,
    }).encode()
    sso = base64.b64encode(payload).decode()
    sig = gate._sig(secret if secret is not None else gate.SSO_SECRET,
                    sso.encode())
    return "/_auth/callback?" + urllib.parse.urlencode({"sso": sso, "sig": sig})


def main() -> int:
    server = ThreadingHTTPServer(("127.0.0.1", PORT), gate.Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    time.sleep(0.3)

    print("\nunauthenticated")
    status, _, _ = request("/_auth/verify")
    check("verify with no cookie is 401", status, 401)
    status, _, _ = request("/_auth/verify", cookie="not-a-real-token")
    check("verify with junk cookie is 401", status, 401)

    print("\nlogin redirect")
    status, headers, _ = request("/_auth/login")
    check("login redirects", status, 302)
    loc = headers.get("Location", "")
    check("redirects to the Discourse SSO provider", loc.startswith(
        "https://discussion.example.org/session/sso_provider?"), True)
    args = urllib.parse.parse_qs(urllib.parse.urlparse(loc).query)
    sso_out = args["sso"][0]
    check("outbound payload is signed with the SSO secret",
          gate._sig(gate.SSO_SECRET, sso_out.encode()), args["sig"][0])
    sent = urllib.parse.parse_qs(base64.b64decode(sso_out).decode())
    check("return url points back at our callback",
          sent["return_sso_url"][0],
          "https://igvfagent.example.org/_auth/callback")
    nonce = sent["nonce"][0]

    print("\nforged and replayed callbacks")
    status, _, _ = request(discourse_callback(nonce, secret=b"wrong-secret"))
    check("a callback signed with the wrong secret is refused", status, 403)
    status, _, _ = request("/_auth/callback?sso=x")
    check("a callback with no signature is refused", status, 400)
    status, headers, _ = request(discourse_callback("never-issued-nonce"))
    check("a callback with an unknown nonce restarts login", status, 302)
    check("...by bouncing to /_auth/login",
          headers.get("Location"), "/_auth/login")

    print("\napproval")
    status, headers, body = request(
        discourse_callback(nonce, groups="trust_level_0,staff"))
    check("a real user outside the group is refused", status, 403)
    check("...and is told which group to ask for",
          b"igvfagent" in body, True)
    check("...and gets no session cookie",
          any("Set-Cookie" in k for k in headers), False)

    print("\nhappy path")
    status, headers, _ = request("/_auth/login")
    nonce = urllib.parse.parse_qs(base64.b64decode(urllib.parse.parse_qs(
        urllib.parse.urlparse(headers["Location"]).query)["sso"][0]).decode()
    )["nonce"][0]
    status, headers, _ = request(discourse_callback(nonce))
    check("an approved user is let in", status, 302)
    check("...and lands on the app", headers.get("Location"), "/")
    setc = headers.get("Set-Cookie", "")
    token = setc.split("=", 1)[1].split(";")[0]
    check("cookie is HttpOnly", "HttpOnly" in setc, True)
    check("cookie is Secure", "Secure" in setc, True)
    check("cookie is SameSite=Lax", "SameSite=Lax" in setc, True)

    print("\nthe issued session")
    status, headers, _ = request("/_auth/verify", cookie=token)
    check("verify now passes", status, 200)
    check("identity is handed to the app", headers.get("X-IGVF-User"), "tester")
    check("email is handed to the app",
          headers.get("X-IGVF-Email"), "tester@example.org")
    check("admin flag is off for a normal member",
          headers.get("X-IGVF-Admin"), "0")

    print("\nreplay of a spent nonce")
    status, headers, _ = request(discourse_callback(nonce))
    check("the same callback cannot be used twice", status, 302)
    check("...it restarts login instead",
          headers.get("Location"), "/_auth/login")

    print("\ntampering with a valid cookie")
    encoded, _, sig = token.rpartition(".")
    forged = gate._b64u(b'{"a":true,"e":"","i":"","n":"","u":"root","x":9999999999}')
    status, _, _ = request("/_auth/verify", cookie=f"{forged}.{sig}")
    check("swapping the body invalidates the signature", status, 401)
    status, _, _ = request("/_auth/verify", cookie=f"{encoded}.{'0' * 64}")
    check("swapping the signature is refused", status, 401)

    print("\nbreak-glass")
    check("bootstrap is absent unless a token is configured",
          request("/_auth/bootstrap?token=anything")[0], 404)
    gate.BOOTSTRAP_TOKEN = "break-glass-token"
    check("a wrong bootstrap token is refused",
          request("/_auth/bootstrap?token=wrong")[0], 403)
    status, headers, _ = request("/_auth/bootstrap?token=break-glass-token")
    check("the right one issues a session", status, 302)
    bt = headers.get("Set-Cookie", "").split("=", 1)[1].split(";")[0]
    status, vheaders, _ = request("/_auth/verify", cookie=bt)
    check("...that verifies", status, 200)
    check("...as a non-admin", vheaders.get("X-IGVF-Admin"), "0")
    gate.BOOTSTRAP_TOKEN = ""

    print("\nsign out")
    status, headers, _ = request("/_auth/logout")
    check("logout clears the cookie",
          "Max-Age=0" in headers.get("Set-Cookie", ""), True)

    server.shutdown()
    print()
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED: {', '.join(FAILED)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
