import importlib.util

import pytest

from mixmogam import _fast


@pytest.fixture
def numpy_fallback(monkeypatch):
    """Load an isolated copy of a module with JIT disabled.

    The process's production module and its compiled dispatchers stay
    untouched for other tests.
    """
    def load(module):
        with monkeypatch.context() as context:
            context.setattr(_fast, "HAS_NUMBA", False)
            spec = importlib.util.spec_from_file_location(f"{module.__name__}_numpy_fallback",
                                                          module.__file__)
            copy = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(copy)
        return copy
    return load


def require_numba_threads(n=2):
    """Skip when the process was started with fewer than ``n`` Numba threads."""
    numba = pytest.importorskip("numba")
    if numba.config.NUMBA_NUM_THREADS < n:
        pytest.skip(f"requires NUMBA_NUM_THREADS >= {n} before Python starts")
    return n
