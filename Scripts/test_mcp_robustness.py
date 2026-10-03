"""The MCP server survives a tool that reads stdin, and handles calls
concurrently.

Regression for a hosted job that burned 4.5 hours and $26: an agent-authored
extension read its content from stdin when an argument was missing; run from
the MCP server it inherited the server's JSON-RPC stream, swallowed it, and
every later call hung until the client's timeout. Tools now get stdin=DEVNULL
and the server answers calls in parallel with a per-call timeout.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main() -> int:
    ok = True

    def check(name, cond):
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'ok  ' if cond else 'FAIL'}  {name}")

    with tempfile.TemporaryDirectory() as td:
        tools = Path(td) / "tools"
        tools.mkdir()
        py = sys.executable
        (tools / "stdin_eater.json").write_text(json.dumps({
            "name": "stdin_eater", "description": "reads all of stdin (bad extension)",
            "parameters": {"type": "object", "properties": {}},
            "command": [py, "-c", "import sys; d=sys.stdin.read(); print('read', len(d))"]}))
        (tools / "slow_tool.json").write_text(json.dumps({
            "name": "slow_tool", "description": "sleeps 3 s",
            "parameters": {"type": "object", "properties": {}},
            "command": [py, "-c", "import time; time.sleep(3); print('slow done')"]}))
        (tools / "fast_tool.json").write_text(json.dumps({
            "name": "fast_tool", "description": "answers at once",
            "parameters": {"type": "object", "properties": {}},
            "command": [py, "-c", "print('fast done')"]}))
        env = {**os.environ, "IGVF_USER_EXT_DIR": td, "PYTHONPATH": str(HERE)}
        srv = subprocess.Popen([py, str(HERE / "mcp_server_skill.py"), "serve", "--all"],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                               text=True, env=env, bufsize=1)

        def send(obj):
            srv.stdin.write(json.dumps(obj) + "\n")
            srv.stdin.flush()

        send({"jsonrpc": "2.0", "id": 0, "method": "initialize", "params": {}})
        json.loads(srv.stdout.readline())
        send({"jsonrpc": "2.0", "id": 1, "method": "tools/call", "params": {"name": "stdin_eater", "arguments": {}}})
        send({"jsonrpc": "2.0", "id": 2, "method": "tools/call", "params": {"name": "slow_tool", "arguments": {}}})
        send({"jsonrpc": "2.0", "id": 3, "method": "tools/call", "params": {"name": "fast_tool", "arguments": {}}})
        send({"jsonrpc": "2.0", "id": 4, "method": "ping"})
        t0, got, order = time.time(), {}, []
        import select
        while len(got) < 4 and time.time() - t0 < 60:
            r, _, _ = select.select([srv.stdout], [], [], 1)
            if r:
                msg = json.loads(srv.stdout.readline())
                got[msg["id"]] = msg
                order.append(msg["id"])
        srv.stdin.close()
        srv.wait(timeout=30)
        text = lambda i: ((got.get(i) or {}).get("result") or {}).get("content", [{}])[0].get("text", "")
        check("a tool that reads stdin gets EOF instead of the protocol stream", text(1).startswith("read 0"))
        check("the server keeps answering after it", 3 in got and 4 in got and "fast done" in text(3))
        check("calls run concurrently: the fast call is not queued behind the slow one",
              2 in got and order.index(3) < order.index(2))
    print("all checks pass" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
