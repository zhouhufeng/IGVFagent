"""A Claude Code job whose session cannot be resumed (it lived in a container
a redeploy replaced) continues in a fresh session from its plan instead of
failing round after round with error_during_execution."""
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        os.environ["IGVF_JOBS_DIR"] = td
        import importlib
        import agent_jobs as aj
        importlib.reload(aj)
        job = aj.create_job("Reproduce X", owner="", orchestrator="claude_code")
        aj.update_job(job["id"], session_id="OLD-SESSION")
        job = aj.load_job(job["id"])
        calls = []
        r = aj.ClaudeCodeRunner(None)

        def fake_once(j, prompt):
            calls.append(j.get("session_id"))
            if j.get("session_id"):
                return aj.RoundResult(stop_reason="error_during_execution", error="error_during_execution")
            return aj.RoundResult(answer="continued", stop_reason="success", session_id="NEW")
        r._run_once = fake_once
        res = r.run_round(job, "JOB CONTINUES", [])
        ok1 = calls == ["OLD-SESSION", ""] and res.answer == "continued"
        ok2 = aj.load_job(job["id"]).get("session_id") == ""
        ok3 = any(e["kind"] == "session" for e in aj.read_events(job["id"], 20))
        for name, ok in (("a lost session is retried once in a fresh session", ok1),
                         ("the stale session id is dropped", ok2), ("the switch is recorded as an event", ok3)):
            print(f"  {'ok  ' if ok else 'FAIL'}  {name}")
        allok = ok1 and ok2 and ok3
    print("all checks pass" if allok else "SOME CHECKS FAILED")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
