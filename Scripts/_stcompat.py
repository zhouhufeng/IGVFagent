"""Streamlit version-compatibility helpers.

Streamlit is replacing the boolean ``use_container_width=`` argument with
a string ``width=`` ("stretch" / "content"). The migration is staggered
in two dimensions at once, which is what makes a blanket find-and-replace
actively harmful:

  * **By version.** Old releases accept only ``use_container_width``.
    Current releases accept both and warn on the old one. Post-removal
    releases accept only ``width`` and raise ``TypeError`` on the old one.
  * **By widget.** Even within one release the coverage is uneven. On
    Streamlit 1.50, ``st.button`` / ``st.dataframe`` / ``st.image`` /
    ``st.pyplot`` accept both, but ``st.plotly_chart`` and
    ``st.altair_chart`` accept **only** ``use_container_width``.

So "just switch everything to width=" breaks the chart widgets on the
very version that prompted the change, and "leave it on
use_container_width" breaks every widget once the removal lands. Pinning
Streamlit avoids the problem only until the pin is raised.

``fit()`` resolves the argument per call, against the signature of the
actual function being invoked, so one source tree runs unmodified on
old, current, and post-removal Streamlit:

    from _stcompat import fit
    st.dataframe(df, **fit(st.dataframe), height=300)
    st.plotly_chart(figure, **fit(st.plotly_chart))
    st.image(path, **fit(st.image, stretch=False))

If a future release drops *both* spellings, ``fit()`` returns an empty
dict and the call still succeeds with that widget's default sizing —
degraded layout, never a crash.
"""

from __future__ import annotations

import inspect
from typing import Any, Callable, Dict

__all__ = ["fit", "supports"]


def _params(fn: Callable[..., Any]):
    """Best-effort parameter mapping. Streamlit decorates several widgets,
    so treat an unreadable signature as 'unknown' rather than failing."""
    try:
        return inspect.signature(fn).parameters
    except (TypeError, ValueError):
        return None


# Cache is keyed on the function object. Streamlit's widget callables are
# module-level singletons, so this is a handful of entries at most, and
# it keeps the per-render cost off the hot path.
_MODE_CACHE: "Dict[int, str]" = {}


def _mode(fn: Callable[..., Any]) -> str:
    key = id(fn)
    cached = _MODE_CACHE.get(key)
    if cached is not None:
        return cached

    params = _params(fn)
    if params is None:
        mode = "none"
    elif "width" in params:
        # Prefer the modern spelling wherever it exists: on releases that
        # accept both, this is also what silences the deprecation warning.
        mode = "width"
    elif "use_container_width" in params:
        mode = "ucw"
    else:
        mode = "none"

    _MODE_CACHE[key] = mode
    return mode


def supports(fn: Callable[..., Any]) -> str:
    """Which spelling ``fn`` takes: ``'width'``, ``'ucw'`` or ``'none'``."""
    return _mode(fn)


def fit(fn: Callable[..., Any], stretch: bool = True) -> "Dict[str, Any]":
    """Keyword arguments that make ``fn`` fill its container (or not).

    ``stretch=True`` fills the available width; ``stretch=False`` sizes to
    the content. Splat the result into the call:
    ``st.dataframe(df, **fit(st.dataframe))``.
    """
    mode = _mode(fn)
    if mode == "width":
        return {"width": "stretch" if stretch else "content"}
    if mode == "ucw":
        return {"use_container_width": stretch}
    return {}
