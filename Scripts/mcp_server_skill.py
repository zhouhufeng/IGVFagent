"""Serve the IGVFagent skill registry over MCP.

Makes the skill library reusable outside this agent. A third-party client —
Claude Desktop, an IDE, another agent framework — connects over the Model
Context Protocol and calls IGVF skills directly, without importing IGVFagent's
runtime or reimplementing its archive handling.

This is the reuse path that turns a curated tool registry into shared
infrastructure rather than a feature of one program. The registry already
carries exactly what MCP needs — name, description, JSON-Schema parameters —
so the mapping is direct and no schema is restated here. A tool added to
IGVFagent appears over MCP automatically.

Run it::

    igvfagent mcp serve                 # stdio, the usual MCP transport
    igvfagent mcp list                  # what would be exposed
    igvfagent mcp manifest              # client config snippet

Client configuration (Claude Desktop's claude_desktop_config.json)::

    {"mcpServers": {"igvfagent": {"command": "igvfagent",
                                   "args": ["mcp", "serve"]}}}

Implemented directly against the JSON-RPC 2.0 wire format rather than an SDK:
the protocol surface needed here is small, and a stdlib implementation keeps
the base package dependency-free — the same reason the skills themselves avoid
non-essential dependencies.

**Exposure is the same code path the agent uses.** Tool calls run through
``_tools.execute``, so an external client is bound by the same declared
schemas and the same subprocess execution as the in-app agent. It is not a
wider door.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

__all__ = ["main", "serve_stdio"]

PROTOCOL_VERSION = "2024-11-05"
SERVER_NAME = "igvfagent"


def _load():
    try:
        from igvfagent import _tools, _llm, __version__
    except Exception:  # pragma: no cover
        import _tools, _llm  # type: ignore
        try:
            from __init__ import __version__  # type: ignore
        except Exception:
            __version__ = "0.0.0"
    return _tools, _llm, __version__


def _exposed_tools(*, all_tools: bool = False):
    """Tools to advertise.

    Defaults to the canonical set the agent itself sees, so an external client
    and the in-app agent are looking at the same registry. ``--all`` lifts the
    cap, which is safe here because MCP clients do not share OpenAI's
    128-function limit — that cap exists for cross-backend parity, not for MCP.
    """
    _tools, _llm, _v = _load()
    tools = _tools.list_tools()
    dicts = [t.to_dict() for t in tools]
    if all_tools:
        return dicts
    return _llm.canonical_tools(dicts) or dicts


def _mcp_tool(entry: dict) -> dict:
    params = entry.get("parameters") or {"type": "object", "properties": {}}
    desc = (entry.get("description") or "").lstrip("★ ").strip()
    return {"name": entry["name"], "description": desc, "inputSchema": params}


# ---------------------------------------------------------------------------
# JSON-RPC
# ---------------------------------------------------------------------------

def _result(req_id, payload):
    return {"jsonrpc": "2.0", "id": req_id, "result": payload}


def _error(req_id, code, message, data=None):
    err = {"code": code, "message": message}
    if data:
        err["data"] = data
    return {"jsonrpc": "2.0", "id": req_id, "error": err}


def handle(request: dict, *, all_tools: bool = False) -> "dict | None":
    """Handle one JSON-RPC request. Returns None for notifications."""
    method = request.get("method")
    req_id = request.get("id")
    params = request.get("params") or {}
    _tools, _llm, version = _load()

    # Notifications carry no id and must not be answered.
    if req_id is None and method and method.startswith("notifications/"):
        return None

    if method == "initialize":
        return _result(req_id, {
            "protocolVersion": PROTOCOL_VERSION,
            "capabilities": {"tools": {"listChanged": False}},
            "serverInfo": {"name": SERVER_NAME, "version": version},
            "instructions": (
                "IGVFagent skills for functional genomics: IGVF Portal and "
                "Catalog, ENCODE, GEO, FAVOR, MaveDB, single-cell and "
                "multiome analysis, variant annotation, knowledge-graph "
                "traversal. Each tool runs a real analysis and writes "
                "artefacts to disk; paths are returned in the result."),
        })

    if method == "tools/list":
        return _result(req_id, {
            "tools": [_mcp_tool(t) for t in _exposed_tools(all_tools=all_tools)]})

    if method == "tools/call":
        name = params.get("name")
        args = params.get("arguments") or {}
        if not name:
            return _error(req_id, -32602, "missing tool name")
        try:
            res = _tools.execute(name, args, timeout=_call_timeout(name, args))
        except KeyError:
            return _error(req_id, -32602, f"unknown tool: {name}")
        except Exception as e:
            return _error(req_id, -32603, f"{type(e).__name__}: {e}")

        stdout = (res.get("stdout") or "").strip()
        stderr = (res.get("stderr") or "").strip()
        exit_code = int(res.get("exit_code") or 0)
        arts = res.get("artifacts") or {}
        flat = [p for v in arts.values() for p in v] if isinstance(arts, dict) else []

        text = stdout or "(no output)"
        if flat:
            text += "\n\nArtefacts:\n" + "\n".join(f"- {p}" for p in flat)
        if exit_code != 0 and stderr:
            text += f"\n\nstderr:\n{stderr[:2000]}"

        # isError reflects the tool's exit status so a client can distinguish
        # "ran and failed" from "ran and returned nothing".
        return _result(req_id, {"content": [{"type": "text", "text": text}],
                                 "isError": exit_code != 0})

    if method == "ping":
        return _result(req_id, {})

    if req_id is None:
        return None
    return _error(req_id, -32601, f"method not found: {method}")


def _call_timeout(name: str, args: dict) -> float:
    """Every call is bounded, so a stuck tool cannot hold a client forever.
    Waiting tools get their own wait plus a margin."""
    base = float(os.environ.get("IGVF_MCP_TOOL_TIMEOUT", "3600"))
    try:
        wait_min = float(args.get("timeout_min") or 0)
    except (TypeError, ValueError):
        wait_min = 0.0
    if name == "job_wait":
        wait_min = wait_min or float(os.environ.get("IGVF_JOB_WAIT_MAX_MIN", "240"))
    if name == "delegate_tasks":
        wait_min = max(wait_min, float(os.environ.get("IGVF_DELEGATE_TIMEOUT_MIN", "90")))
    return max(base, wait_min * 60 + 300)


def serve_stdio(*, all_tools: bool = False) -> int:
    """Read JSON-RPC from stdin, write responses to stdout.

    Anything the server needs to say about itself goes to stderr: stdout is
    the protocol channel, and one stray print corrupts the stream.
    """
    print(f"{SERVER_NAME} MCP server ready on stdio "
          f"({len(_exposed_tools(all_tools=all_tools))} tools)", file=sys.stderr)
    # Requests are handled concurrently: a client issues parallel tool calls,
    # and one long call (delegate_tasks, a download) must not stall the rest.
    # Responses carry their request id, so their order does not matter; the
    # lock keeps each JSON line whole on stdout.
    import concurrent.futures
    import threading
    out_lock = threading.Lock()
    pool = concurrent.futures.ThreadPoolExecutor(
        max_workers=int(os.environ.get("IGVF_MCP_WORKERS", "8")))

    def emit(obj: dict) -> None:
        with out_lock:
            sys.stdout.write(json.dumps(obj) + "\n")
            sys.stdout.flush()

    def work(request: dict) -> None:
        try:
            response = handle(request, all_tools=all_tools)
        except Exception as e:  # never let one bad call kill the server
            traceback.print_exc(file=sys.stderr)
            response = _error(request.get("id"), -32603,
                              f"internal error: {type(e).__name__}: {e}")
        if response is not None:
            emit(response)

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            request = json.loads(line)
        except json.JSONDecodeError:
            emit(_error(None, -32700, "parse error"))
            continue
        if request.get("method") in ("initialize", "tools/list") or request.get("id") is None:
            work(request)  # handshake and notifications stay in order
        else:
            pool.submit(work, request)
    pool.shutdown(wait=True)
    return 0


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def cmd_list(args) -> int:
    tools = _exposed_tools(all_tools=args.all)
    for t in tools:
        desc = (t.get("description") or "").lstrip("★ ")
        print(f"  {t['name']:34s} {desc[:82]}")
    print(f"\n{len(tools)} tool(s) would be exposed over MCP")
    return 0


def cmd_manifest(args) -> int:
    exe = "igvfagent"
    cfg = {"mcpServers": {SERVER_NAME: {"command": exe,
                                         "args": ["mcp", "serve"]}}}
    if args.all:
        cfg["mcpServers"][SERVER_NAME]["args"].append("--all")
    print(json.dumps(cfg, indent=2))
    print("\n# Claude Desktop: merge into claude_desktop_config.json",
          file=sys.stderr)
    print("# Other clients: any MCP client supporting stdio transport",
          file=sys.stderr)
    return 0


def cmd_serve(args) -> int:
    return serve_stdio(all_tools=args.all)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="igvfagent mcp",
        description="Serve the IGVFagent skill registry over the Model "
                    "Context Protocol so other agents can call IGVF skills.")
    sub = p.add_subparsers(dest="cmd", required=True)
    for name, helptext in (("serve", "Run the MCP server on stdio"),
                            ("list", "List the tools that would be exposed"),
                            ("manifest", "Print a client config snippet")):
        s = sub.add_parser(name, help=helptext)
        s.add_argument("--all", action="store_true",
                       help="Expose every registered tool, not just the "
                            "canonical capped set (MCP clients have no "
                            "128-function limit)")
    args = p.parse_args(argv)
    return {"serve": cmd_serve, "list": cmd_list,
            "manifest": cmd_manifest}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
