"""Compute-resource snapshot for the UI: CPU, memory, disk.

WHY NOT psutil ALONE. Inside a container psutil reports the HOST's cpu_count
and total memory, not the limits the container actually runs under. On this
deployment that is 8 vCPU / 30 GB host against a 6 CPU / 22 GB cgroup, so a
psutil-only panel would show the app using a third of "available" memory when
it is really using half. cgroup v2 is read first and psutil is the fallback
for a local (non-container) install.

WHY memory.current IS NOT THE ANSWER. It counts page cache. Measured here it
read 13.8 GB against a 22 GB limit -- 63%, which looks like the container is
nearly full -- while `docker stats` said 2.25 GiB / 22 GiB. The difference is
inactive_file: 11.4 GB of reclaimable cache from reading FASTQ and parquet.
Subtracting it is what Docker itself does, and it is the number a reader
should see, because cache is returned under pressure rather than being a
shortfall.
"""
from __future__ import annotations

import os
import shutil
import time
from pathlib import Path
from typing import Any, Optional

_CG = Path("/sys/fs/cgroup")
# Previous CPU sample, so a percentage can be computed from the delta.
_prev: "Optional[tuple[float, int]]" = None


def _read(path: Path) -> Optional[str]:
    try:
        return path.read_text().strip()
    except OSError:
        return None


def _cgroup_mem() -> "Optional[tuple[float, Optional[float]]]":
    """(used_gb, limit_gb) from cgroup v2, cache-corrected."""
    cur = _read(_CG / "memory.current")
    if cur is None:
        return None
    try:
        used = int(cur)
    except ValueError:
        return None
    stat = _read(_CG / "memory.stat") or ""
    for line in stat.splitlines():
        if line.startswith("inactive_file "):
            try:
                used -= int(line.split()[1])
            except (IndexError, ValueError):
                pass
            break
    raw_max = _read(_CG / "memory.max")
    limit = None
    if raw_max and raw_max != "max":
        try:
            limit = int(raw_max) / 1e9
        except ValueError:
            limit = None
    return max(used, 0) / 1e9, limit


def _cgroup_cpu_limit() -> Optional[float]:
    raw = _read(_CG / "cpu.max")
    if not raw:
        return None
    parts = raw.split()
    if len(parts) != 2 or parts[0] == "max":
        return None
    try:
        return int(parts[0]) / int(parts[1])
    except (ValueError, ZeroDivisionError):
        return None


def _cgroup_cpu_percent() -> Optional[float]:
    """CPU use as a percentage of the container's OWN quota.

    Reported against the quota, not against host cores: "80% of the 6 CPUs
    this container may use" is actionable; "60% of 8 host cores" is not, since
    the container can never reach the other two.
    """
    global _prev
    stat = _read(_CG / "cpu.stat")
    if not stat:
        return None
    usec = None
    for line in stat.splitlines():
        if line.startswith("usage_usec"):
            try:
                usec = int(line.split()[1])
            except (IndexError, ValueError):
                return None
            break
    if usec is None:
        return None
    now = time.monotonic()
    prev = _prev
    _prev = (now, usec)
    if prev is None:
        return None                      # first call has no interval yet
    dt = now - prev[0]
    if dt <= 0.05:                       # too short to be meaningful
        return None
    quota = _cgroup_cpu_limit() or (os.cpu_count() or 1)
    busy = (usec - prev[1]) / 1e6        # CPU-seconds used in the interval
    return max(0.0, min(100.0, 100.0 * busy / (dt * quota)))


def snapshot(data_path: str = "/workspace/Data") -> "dict[str, Any]":
    """One reading. Every field may be None when it cannot be determined."""
    out: "dict[str, Any]" = {"source": "cgroup"}

    mem = _cgroup_mem()
    if mem is None:
        out["source"] = "psutil"
        try:
            import psutil
            vm = psutil.virtual_memory()
            out["mem_used_gb"] = (vm.total - vm.available) / 1e9
            out["mem_limit_gb"] = vm.total / 1e9
        except Exception:                                   # noqa: BLE001
            out["mem_used_gb"] = out["mem_limit_gb"] = None
    else:
        out["mem_used_gb"], out["mem_limit_gb"] = mem

    out["cpu_limit"] = _cgroup_cpu_limit() or os.cpu_count()
    pct = _cgroup_cpu_percent()
    if pct is None and out["source"] == "psutil":
        try:
            import psutil
            pct = psutil.cpu_percent(interval=None) or None
        except Exception:                                   # noqa: BLE001
            pct = None
    out["cpu_percent"] = pct

    if out.get("mem_used_gb") is not None and out.get("mem_limit_gb"):
        out["mem_percent"] = 100.0 * out["mem_used_gb"] / out["mem_limit_gb"]
    else:
        out["mem_percent"] = None

    p = Path(data_path)
    while not p.exists() and p != p.parent:
        p = p.parent
    try:
        du = shutil.disk_usage(str(p))
        out["disk_free_gb"] = du.free / 1e9
        out["disk_total_gb"] = du.total / 1e9
        out["disk_percent"] = 100.0 * (du.total - du.free) / du.total
    except OSError:
        out["disk_free_gb"] = out["disk_total_gb"] = out["disk_percent"] = None
    return out


__all__ = ["snapshot"]
