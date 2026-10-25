"""Numba availability shim for the built-in coalescent backend."""

from __future__ import annotations

__all__ = ["HAVE_NUMBA", "_jit"]

try:
    from numba import njit

    HAVE_NUMBA = True

    def _jit(func):
        return njit(cache=True)(func)

except ImportError:  # pure-Python fallback: correct, just slower
    HAVE_NUMBA = False

    def _jit(func):
        return func
