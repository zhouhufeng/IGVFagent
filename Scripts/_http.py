"""Shared HTTP helpers for IGVFagent skills.

The single feature here is :func:`prefer_ipv4_dns` — a process-wide
monkeypatch of :func:`socket.getaddrinfo` that biases name resolution
toward IPv4 addresses. Why this exists:

Several IGVF services (notably ``api.catalogkg.igvf.org``,
``data.igvf.org``, ``api.data.igvf.org``) publish both ``A`` and
``AAAA`` DNS records. On networks where IPv6 is not actually routable
(common on macOS laptops, university networks behind IPv4-only NAT,
cellular tethering, …) Python's stdlib :mod:`urllib` silently waits
for the IPv6 socket-connect to time out (10–40 s) before falling
back to IPv4. The ``timeout=`` parameter passed to
:func:`urllib.request.urlopen` does NOT cover that connect-phase
wait, so the call appears to hang even when callers think they've
asked for a 5 s ceiling.

`curl` implements RFC 8305 Happy-Eyeballs (parallel IPv4+IPv6 with
sub-second fallback). Python's stdlib does not. Until Python ships
HE, the cheapest workaround is to skip IPv6 entirely whenever the
calling network is IPv4-only. We do that here by filtering
:func:`socket.getaddrinfo` results to ``AF_INET`` (IPv4).

Opt out by setting ``IGVF_PREFER_IPV4=0`` in the environment if you
have a legitimately dual-stack network and want IPv6 attempted.

This module is import-once-safe: applying the monkeypatch twice is a
no-op.
"""

from __future__ import annotations

import os
import socket
from typing import Any


_ORIG_GETADDRINFO = None


def prefer_ipv4_dns() -> bool:
    """Install the IPv4-preferred resolver. Returns True if installed.

    Set ``IGVF_PREFER_IPV4=0`` to disable.
    """
    global _ORIG_GETADDRINFO
    if os.environ.get("IGVF_PREFER_IPV4", "1") == "0":
        return False
    if _ORIG_GETADDRINFO is not None:
        # Already installed
        return True
    _ORIG_GETADDRINFO = socket.getaddrinfo

    def _ipv4_only(host: Any, port: Any, family: int = 0,
                    *args: Any, **kwargs: Any) -> Any:
        # If caller explicitly asked for a non-zero family, honor it.
        # Otherwise force AF_INET (IPv4) so the resolver returns only A
        # records and we never attempt an unroutable IPv6 connect.
        if family == 0:
            family = socket.AF_INET
        assert _ORIG_GETADDRINFO is not None  # for mypy
        return _ORIG_GETADDRINFO(host, port, family, *args, **kwargs)

    socket.getaddrinfo = _ipv4_only  # type: ignore[assignment]
    return True


# Apply at import time so any IGVFagent skill that imports this module
# (and downstream modules that import a skill that imports this) gets
# the IPv4-preferred resolver for the rest of the process.
prefer_ipv4_dns()
