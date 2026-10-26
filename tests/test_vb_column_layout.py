"""The column-major (P, n) VB workspace matches the sample-major one.

Fits of 2 to ``_COLUMN_MAJOR_MAX`` model columns keep the residual, its
storage-precision image and the updates column-major. The arithmetic is the
sample-major layout's with each GEMM transposed: Accelerate (macOS) rounds
these as the sample-major products from a few thousand samples, so the fits
are identical there; other BLAS builds may round the transposed GEMMs apart.
The dense oracle of tests/test_vb_workspace.py pins the sample-major layout.
"""

import importlib.util
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

from mixmogam import _fast, _vb
from mixmogam._loco import LocoGenotypes
from mixmogam.genotypes import Genotypes

EXACT = sys.platform == "darwin"


def _fallback(monkeypatch, name):
    """An isolated copy of the module with JIT disabled."""
    with monkeypatch.context() as context:
        context.setattr(_fast, "HAS_NUMBA", False)
        spec = importlib.util.spec_from_file_location(name, _vb.__file__)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


def _sample_major(monkeypatch, *modules):
    for module in modules or (_vb,):
        monkeypatch.setattr(module, "_COLUMN_MAJOR_MAX", 1)


def _assert_same(got, want, single):
    if EXACT:
        np.testing.assert_array_equal(got, want)
    else:  # transposed GEMMs may round apart: float32 products, then the fit
        tolerance = 1e-4 if single else 1e-10
        np.testing.assert_allclose(got, want, rtol=tolerance,
                                   atol=tolerance * np.max(np.abs(want), initial=0.0))


def _assert_same_fit(actual, expected, single):
    assert actual.keys() == expected.keys()
    if EXACT:
        assert actual["iterations"] == expected["iterations"]
        assert actual["converged"] is expected["converged"]
    for key in ("beta", "resid", "prediction", "rel_change"):
        if key in expected:
            _assert_same(actual[key], expected[key], single)
    if expected["mask"] is None:
        assert actual["mask"] is None
    else:
        np.testing.assert_array_equal(actual["mask"], expected["mask"])


@pytest.mark.parametrize("numba", [True, False])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("folded", [False, True])
@pytest.mark.parametrize("p", [1, 2, 3, 4, 5, 8, 11])
def test_column_refresh_matches_sample_major_refresh(monkeypatch, numba, dtype, folded, p):
    if numba and not _fast.HAS_NUMBA:
        pytest.skip("requires Numba")
    module = _vb if numba else _fallback(monkeypatch, "mixmogam._vb_column_refresh_fallback")
    rng = np.random.default_rng(1301 + p)
    n = 301
    R = rng.standard_normal((n, p))
    work = R.astype(dtype)
    delta = rng.standard_normal((n, p)).astype(dtype)
    folds = rng.integers(0, 3, n) if folded else None
    col_fold = rng.integers(-1, 3, p)
    change = rng.random(p)
    initial = R.copy()
    columns = [np.ascontiguousarray(a.T) for a in (R, work, delta)] + [change.copy()]
    module._refresh_residual(R, work, delta, folds, col_fold, change)
    module._refresh_columns(columns[0], columns[1], columns[2], folds, col_fold, columns[3])
    # Four columns per pass, each summed in sample order: the same bits.
    np.testing.assert_array_equal(columns[0], R.T)
    np.testing.assert_array_equal(columns[1], work.T)
    np.testing.assert_array_equal(columns[3], change)
    for c, f in enumerate(col_fold):
        held = folds == f if folded else np.zeros(n, dtype=bool)
        np.testing.assert_array_equal(columns[0][c, held], initial[held, c])
        assert np.all(columns[0][c, ~held] != initial[~held, c])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("p", [3, 10, 25])
def test_column_products_multiply_the_whole_decode_buffer(dtype, p):
    rng = np.random.default_rng(1307 + p)
    n, b = 5000, 128
    rows = rng.standard_normal((b, n)).astype(dtype)
    work = rng.standard_normal((p, n)).astype(dtype)
    flat = np.empty(b * p, dtype=dtype)
    full = _vb._column_products(work, rows, b, flat).copy()
    single = dtype == np.float32
    _assert_same(full, rows @ np.ascontiguousarray(work.T), single)
    for k in (1, 2, 45, 127):
        got = _vb._column_products(work, rows, k, flat)
        assert got.shape == (k, p) and np.shares_memory(got, flat)
        # The sample-major products of the block alone; one row as its GEMV.
        _assert_same(got, rows[:k] @ np.ascontiguousarray(work.T), single)
        if k > 1:
            # Rows past the block only fill discarded columns of one GEMM.
            np.testing.assert_array_equal(got, full[:k])
            stale = rows.copy()
            stale[k:] = rng.standard_normal((b - k, n))
            np.testing.assert_array_equal(_vb._column_products(work, stale, k, flat), full[:k])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_column_backward_splits_samples_like_the_serial_product(dtype):
    rng = np.random.default_rng(1311)
    n, p = 1003, 5
    Z = rng.standard_normal((128, n)).astype(dtype)
    with ThreadPoolExecutor(max_workers=4) as pool:
        for workers in (2, 3, 4):
            workspace = _vb._GemmWorkspace(n, p, 128, dtype, pool, workers, True)
            assert not hasattr(workspace, "work")  # no sample-major scratch
            for k in (128, 45):
                changes = rng.standard_normal((p, k)).astype(dtype)
                expected = np.empty((p, n), dtype=dtype)
                np.matmul(changes, Z[:k], out=expected)
                actual = np.full((p, n), np.nan, dtype=dtype)
                workspace.backward_columns(changes, Z[:k], actual)
                # Each product sums the same k terms whichever worker owns it.
                _assert_same(actual, expected, dtype == np.float32)
    assert Z.flags.writeable


def test_raw_slices_decode_into_the_leading_rows_of_a_given_buffer():
    rng = np.random.default_rng(1319)
    n, m = 57, 300
    G = rng.binomial(2, rng.uniform(0.1, 0.5, m), size=(n, m)).astype(np.int8)
    groups = np.repeat([0, 1, 2], [130, 129, 41])
    lg = LocoGenotypes(Genotypes(G, chromosome=groups), groups, np.ones((n, 1)) / np.sqrt(n))
    expected = [(idx.copy(), g, z.copy()) for idx, g, z in lg.raw_slices(64)]
    out = np.full((64, n), 7.0, dtype=lg.dtype)
    before = out.copy()
    sizes = []
    for (idx, g, z), (want_idx, want_g, want_z) in zip(lg.raw_slices(64, out=out), expected, strict=True):
        k = idx.size
        assert z.base is out and z.shape == (k, n)
        np.testing.assert_array_equal(idx, want_idx)
        assert g == want_g
        np.testing.assert_array_equal(z, want_z)
        np.testing.assert_array_equal(out[k:], before[k:])  # the other rows are left alone
        before = out.copy()
        sizes.append(k)
    assert sizes == [64, 64, 2, 64, 64, 1, 41]


@pytest.fixture(scope="module")
def panel():
    """5,000 samples, three groups whose last sub-blocks hold 1, 45 and
    102 variants (a one-row tail goes through GEMV in the sample-major
    layout), missing calls, three folds, a covariate basis of rank 3."""
    rng = np.random.default_rng(1321)
    n = 5000
    sizes = (257, 173, 230)
    m = sum(sizes)
    G = rng.binomial(2, rng.uniform(0.05, 0.5, m), size=(n, m)).astype(np.int8)
    G[rng.random(G.shape) < 0.01] = -1
    groups = np.repeat(np.arange(3), sizes)
    Q = np.linalg.qr(np.column_stack([np.ones(n), rng.standard_normal((n, 2))]))[0]
    causal = rng.choice(m, 40, replace=False)
    g = np.where(G[:, causal] < 0, 0, G[:, causal]) @ rng.standard_normal(40)
    y = 0.5 ** 0.5 * (g - g.mean()) / g.std() + 0.5 ** 0.5 * rng.standard_normal(n)
    y -= Q @ (Q.T @ y)
    folds = rng.permutation(n) % 3
    scale = np.exp(rng.normal(0, 0.3, n))
    return G, groups, Q, y, folds, scale


def _engine(panel, dtype=np.float32, packed=False, weighted=False, intercept=False, cache=10**8,
            n_threads=1):
    G, groups, Q, _, folds, scale = panel
    Q = Q[:, :1] if intercept else Q
    gt = Genotypes(G, chromosome=groups, packed=packed)
    lg = LocoGenotypes(gt, groups, Q, dtype=dtype, row_scale=scale if weighted else None)
    return _vb.VBEngine(lg, folds=folds, gram_cache_bytes=cache, n_threads=n_threads)


def _layout_fit(engine, y, col_fold, col_group, max_iter=200, **kwargs):
    p = col_fold.size
    prior = np.column_stack([np.linspace(0.05, 0.9, p), np.full(p, 30.0),
                             np.linspace(0.2, 1.0, p) / engine.lg.m])
    return engine.fit(np.repeat(y[:, None], p, axis=1), col_fold, col_group, _vb.PRIOR_ENET,
                      prior, np.linspace(0.6, 0.9, p), snp_scale=np.linspace(0.5, 1.5, engine.lg.m),
                      max_iter=max_iter, tol=1e-9, **kwargs)


# Cross-validation columns (held-out folds), LOCO columns, the two mixed,
# and HRATT's cross-fitted layout (every fold held out once).
COLUMNS = {"cv": (np.array([0, 0, 0, 0]), np.full(4, -1)),
           "loco": (np.full(3, -1), np.arange(3)),
           "mixed": (np.array([2, 0, -1, 1, -1, 2, 0]), np.array([-1, 1, 2, 0, 0, -1, 2])),
           "crossfit": (np.arange(3), np.full(3, -1))}


# With sampling weights and no covariates, Q'R~ from the transposed residual
# (a GEMV with one column) would round apart, and the fit with it.
@pytest.mark.parametrize("variant", ["plain", "uncached", "float64", "packed", "weighted",
                                     "weighted_intercept", "first_sweep", "effects_only"])
@pytest.mark.parametrize("columns", list(COLUMNS))
def test_column_major_fit_matches_sample_major(panel, monkeypatch, columns, variant):
    y = panel[3]
    col_fold, col_group = COLUMNS[columns]
    options = {"plain": {}, "uncached": {"cache": 0}, "float64": {"dtype": np.float64},
               "packed": {"packed": True}, "weighted": {"weighted": True},
               "weighted_intercept": {"weighted": True, "intercept": True},
               "first_sweep": {}, "effects_only": {}}[variant]
    fit = {"first_sweep": {"max_iter": 1}, "effects_only": {"_compute_prediction": False}}.get(variant, {})
    single = options.get("dtype", np.float32) == np.float32
    calls = []
    original = _vb._refresh_columns
    monkeypatch.setattr(_vb, "_refresh_columns", lambda *args: calls.append(1) or original(*args))
    actual = _layout_fit(_engine(panel, **options), y, col_fold, col_group, **fit)
    assert calls
    _sample_major(monkeypatch)
    calls.clear()
    expected = _layout_fit(_engine(panel, **options), y, col_fold, col_group, **fit)
    assert not calls
    assert expected["converged"] is (variant != "first_sweep")
    _assert_same_fit(actual, expected, single)


def test_numpy_fallback_column_major_fit_matches_sample_major(panel, monkeypatch):
    fallback = _fallback(monkeypatch, "mixmogam._vb_column_fit_fallback")
    G, groups, Q, y, folds, _ = panel
    lg = LocoGenotypes(Genotypes(G, chromosome=groups), groups, Q)
    col_fold, col_group = COLUMNS["mixed"]
    actual = _layout_fit(fallback.VBEngine(lg, folds=folds), y, col_fold, col_group, max_iter=8)
    _sample_major(monkeypatch, fallback)
    expected = _layout_fit(fallback.VBEngine(lg, folds=folds), y, col_fold, col_group, max_iter=8)
    _assert_same_fit(actual, expected, True)


@pytest.mark.parametrize("p, columns", [(1, False), (2, True), (32, True), (33, False)])
def test_layout_follows_the_number_of_model_columns(panel, monkeypatch, p, columns):
    seen = []
    for name in ("_refresh_columns", "_refresh_residual"):
        original = getattr(_vb, name)
        monkeypatch.setattr(_vb, name, lambda *args, _name=name, _f=original: seen.append(_name) or _f(*args))
    engine = _engine(panel)
    col_fold = np.arange(p) % 3
    _layout_fit(engine, panel[3], col_fold, np.full(p, -1), max_iter=1, _compute_prediction=False)
    assert set(seen) == {"_refresh_columns" if columns else "_refresh_residual"}


@pytest.mark.parametrize("columns", ["mixed", "crossfit"])
def test_threaded_column_major_fit_matches_serial(panel, monkeypatch, columns):
    numba = pytest.importorskip("numba")
    if numba.config.NUMBA_NUM_THREADS < 4:
        pytest.skip("requires NUMBA_NUM_THREADS >= 4 before Python starts")
    pytest.importorskip("threadpoolctl")
    y = panel[3]
    col_fold, col_group = COLUMNS[columns]
    expected = _layout_fit(_engine(panel), y, col_fold, col_group)
    # Below the production threshold: the worker pool splits each backward
    # product's samples, and every product keeps the serial reduction.
    monkeypatch.setattr(_vb, "_GEMM_MIN_SAMPLES", 1)
    pools = []
    original = _vb.ThreadPoolExecutor

    def track(**kwargs):
        pools.append(original(**kwargs))
        return pools[-1]

    monkeypatch.setattr(_vb, "ThreadPoolExecutor", track)
    actual = _layout_fit(_engine(panel, n_threads=4), y, col_fold, col_group)
    assert len(pools) == 1
    _assert_same_fit(actual, expected, True)
