#!/usr/bin/env python3
"""One heavy analysis at a time, so a second one cannot OOM the container.

Measured on the deployment: a full single-cell pipeline on a 1.46 GB h5ad
peaked at 18.2 GB of the container's 22 GB limit. Two of those overlapping
exceeds the cap, and Docker's response is to OOM-kill the container -- which
ends EVERY user's session, not just the second job's. On a shared demo that
turns one person's click into everyone's outage.

So heavy work takes a slot. The second caller waits (or is told plainly that
it is waiting) instead of racing for memory. A wait affects one user; an OOM
affects all of them.

The lock is a file holding the holder's pid and what it is doing, because a
lock nobody can explain is worse than none: an operator who finds work
blocked needs to see what holds it. A holder whose process is gone has its
lock broken automatically -- a crashed analysis must not wedge the site
until someone notices.
"""

from __future__ import annotations

import contextlib
import json
import os
import time
from pathlib import Path
from typing import Iterator, Optional

ROOT = Path(os.environ.get("IGVF_PROJECT_ROOT")
            or Path(__file__).resolve().parents[1]).resolve()
LOCK_DIR = ROOT / "Data" / ".locks"

# How many heavy analyses may run at once. 1 by default: the peak footprint
# measured is 18 GB against a 22 GB cap, so there is no room for a second.
SLOTS = max(1, int(os.environ.get("IGVF_HEAVY_SLOTS", "1")))
# How long to wait for a slot before giving up, in seconds.
WAIT = float(os.environ.get("IGVF_HEAVY_WAIT", "1800"))


def _slot_path(i: int) -> Path:
    return LOCK_DIR / f"heavy_{i}.json"


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(int(pid), 0)
    except (OSError, TypeError, ValueError):
        return False
    return True


def _read(path: Path) -> "Optional[dict]":
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def holders() -> "list[dict]":
    """Who currently holds a slot. Stale entries are reported as such."""
    out = []
    for i in range(SLOTS):
        rec = _read(_slot_path(i))
        if not rec:
            continue
        rec["stale"] = not _pid_alive(rec.get("pid", -1))
        out.append(rec)
    return out


def _try_claim(i: int, what: str, label: str) -> bool:
    """Claim slot i atomically, breaking it only if its holder is gone."""
    path = _slot_path(i)
    LOCK_DIR.mkdir(parents=True, exist_ok=True)
    rec = _read(path)
    if rec is not None:
        if _pid_alive(rec.get("pid", -1)):
            return False
        # The holder died. Breaking the lock here is safe and necessary:
        # a crashed analysis would otherwise block the slot indefinitely.
        try:
            path.unlink()
        except OSError:
            return False
    payload = json.dumps({"pid": os.getpid(), "what": what, "label": label,
                          "started": time.strftime("%Y-%m-%d %H:%M:%S")})
    try:
        # O_EXCL is the atomic part: two processes racing here, exactly one
        # creates the file and the other sees FileExistsError.
        fd = os.open(str(path), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        return False
    except OSError:
        return False
    with os.fdopen(fd, "w") as fh:
        fh.write(payload)
    return True


@contextlib.contextmanager
def heavy_slot(what: str, label: str = "", wait: Optional[float] = None,
                on_wait=None) -> Iterator[bool]:
    """Hold a heavy-analysis slot for the duration of the block.

    Yields True when a slot was taken, False when the wait expired -- the
    caller decides whether that is fatal. It never raises on contention,
    because a queued job failing with a traceback reads like a bug rather
    than like a queue.
    """
    deadline = time.time() + (WAIT if wait is None else wait)
    mine = -1
    announced = False
    while True:
        for i in range(SLOTS):
            if _try_claim(i, what, label):
                mine = i
                break
        if mine >= 0 or time.time() >= deadline:
            break
        if not announced and on_wait:
            on_wait(holders())
            announced = True
        time.sleep(3)
    try:
        yield mine >= 0
    finally:
        if mine >= 0:
            try:
                rec = _read(_slot_path(mine))
                if rec and rec.get("pid") == os.getpid():
                    _slot_path(mine).unlink()
            except OSError:
                pass


def describe() -> str:
    """One line for a status display."""
    hs = holders()
    live = [h for h in hs if not h.get("stale")]
    if not live:
        return f"0/{SLOTS} heavy slots in use"
    parts = ", ".join(f"{h.get('what')} ({h.get('label') or 'unlabelled'}) "
                      f"since {h.get('started')}" for h in live)
    return f"{len(live)}/{SLOTS} heavy slots in use: {parts}"


__all__ = ["heavy_slot", "holders", "describe", "SLOTS", "WAIT"]
