"""Compiled covariate corrections match NumPy's, column set by column set.

VBEngine.fit corrects each sub-block's products, and the changes of each
update, by rank-q covariate terms. When NumPy's calls per column set cost
more than their arithmetic (``_CORRECTION_LOOP_MAX``), two Numba loops
correct every column in one call each. Their sums run over the covariates
and SNPs in order instead of BLAS's, so fits move at rounding level; with
only an intercept and no sampling weights the corrections are near no-ops
and these fits are unchanged. Without Numba, NumPy's products remain.
"""


import numpy as np
import pytest

from mixmogam import _fast, _vb
from mixmogam._loco import LocoGenotypes
from mixmogam.genotypes import Genotypes

needs_numba = pytest.mark.skipif(not _fast.HAS_NUMBA, reason="requires Numba")


def _column_sets(gidx):
    sets = []
    for g in np.unique(gidx):
        cols = np.flatnonzero(gidx == g)
        sets.append((int(g), slice(None) if cols.size == gidx.size else cols))
    return sets


def _numpy_products(U, coefs, coef, A, s_proj, train, gidx):
    """VBEngine.fit's NumPy products, one column set at a time."""
    for g, cols in _column_sets(gidx):
        T = coefs[0] - coefs[g] if g else coefs[0]
        U[:, cols] += T @ A[:, cols] - coef @ (s_proj[:, cols] + train[g] @ A[:, cols])


def _numpy_changes(coefs, coef, d, s_proj, A, train, gidx, change):
    e = coef.T @ d
    for g, cols in _column_sets(gidx):
        T = coefs[0] - coefs[g] if g else coefs[0]
        W = T.T @ d[:, cols]
        s_proj[:, cols] -= W
        ec = e[:, cols]
        change[cols] += (np.einsum("kp,kp->p", ec, train[g] @ ec)
                         - 2.0 * np.einsum("kp,kp->p", W, ec))
    A += e


# Gram slots per column: HRATT's hold-out fit (one slot), LOCO columns (no
# held-out fold, slot 0), a mix that leaves slot 3 unused, HRATT's
# cross-fit (a slot per column) and its strong-structure refit.
SLOTS = {"holdout": np.full(4, 1), "loco": np.zeros(3, dtype=np.int64),
         "mixed": np.array([4, 1, 0, 2, 0, 4, 1]), "crossfit": np.arange(1, 6),
         "refit": np.repeat(np.arange(1, 6), 3)}


@needs_numba
@pytest.mark.parametrize("k", [128, 45, 1])
@pytest.mark.parametrize("q", [1, 2, 5, 11])
@pytest.mark.parametrize("layout", list(SLOTS))
def test_compiled_corrections_match_numpy_column_sets(layout, q, k):
    rng = np.random.default_rng(1409 + 7 * q + k)
    gidx = SLOTS[layout].astype(np.int64)
    P, slots = gidx.size, int(gidx.max()) + 1
    coefs = rng.standard_normal((slots, k, q))
    train = np.eye(q) - 0.05 * rng.standard_normal((slots, q, q))
    train[0] = np.eye(q)
    unused = np.setdiff1d(np.arange(1, slots), gidx)
    coefs[unused] = train[unused] = np.nan  # never read
    coef = rng.standard_normal((k, q))
    A = 1e-2 * rng.standard_normal((q, P))
    s_proj = rng.standard_normal((q, P))
    U = np.asfortranarray(rng.standard_normal((128, P)))
    d = rng.standard_normal((k, P)).astype(np.float32).astype(np.float64)
    change = rng.random(P)
    expected = [U.copy(order="F"), A.copy(), s_proj.copy(), change.copy()]
    _numpy_products(expected[0][:k], coefs, coef, expected[1], expected[2], train, gidx)
    _numpy_changes(coefs, coef, d, expected[2], expected[1], train, gidx, expected[3])
    tmp, e = np.full(q, np.nan), np.empty((q, P))
    _vb._correct_products(U[:k], coefs, coef, A, s_proj, train, gidx, tmp)
    np.matmul(coef.T, d, out=e)
    _vb._correct_changes(coefs, d, e, s_proj, A, train, gidx, change, tmp)
    np.testing.assert_array_equal(A, expected[1])  # the same e, added
    np.testing.assert_array_equal(U[k:], expected[0][k:])  # rows past the block untouched
    for got, want in zip((U, s_proj, change), (expected[0], expected[2], expected[3])):
        np.testing.assert_allclose(got, want, rtol=0, atol=1e-13 * np.abs(want).max())
    if q == 1:  # single products, no sums: the products are NumPy's
        np.testing.assert_array_equal(U, expected[0])


@pytest.fixture(scope="module")
def panel():
    """600 samples, three groups (sub-blocks of 128, then 42, 12 and 110
    variants), missing calls, five folds, an intercept and four covariates,
    and sampling weights."""
    rng = np.random.default_rng(1423)
    n, sizes = 600, (170, 140, 110)
    m = sum(sizes)
    G = rng.binomial(2, rng.uniform(0.05, 0.5, m), size=(n, m)).astype(np.int8)
    G[rng.random(G.shape) < 0.01] = -1
    groups = np.repeat(np.arange(3), sizes)
    X = np.column_stack([np.ones(n), rng.standard_normal((n, 4))])
    causal = rng.choice(m, 30, replace=False)
    g = np.where(G[:, causal] < 0, 0, G[:, causal]) @ rng.standard_normal(30)
    y = 0.5 ** 0.5 * (g - g.mean()) / g.std() + 0.5 ** 0.5 * rng.standard_normal(n)
    folds = rng.permutation(n) % 5
    scale = np.exp(rng.normal(0, 0.4, n))
    return G, groups, X, y, folds, scale


def _genotypes(panel, q=5, weighted=False, packed=False, dtype=np.float32):
    G, groups, X, y, _, scale = panel
    row_scale = scale if weighted else None
    Q = np.linalg.qr(X[:, :q] * (1.0 if row_scale is None else row_scale[:, None]))[0] if q else None
    lg = LocoGenotypes(Genotypes(G, chromosome=groups, packed=packed), groups, Q, dtype=dtype,
                       row_scale=row_scale)
    return lg, y if Q is None else y - Q @ (Q.T @ y)


def _fit(engine, y, col_fold, col_group, **kwargs):
    p = col_fold.size
    prior = np.column_stack([np.linspace(0.05, 0.9, p), np.full(p, 30.0),
                             np.linspace(0.2, 1.0, p) / engine.lg.m])
    return engine.fit(np.repeat(y[:, None], p, axis=1), col_fold, col_group, _vb.PRIOR_ENET,
                      prior, np.linspace(0.6, 0.9, p), snp_scale=np.linspace(0.5, 1.5, engine.lg.m),
                      max_iter=200, tol=1e-9, **kwargs)


def _assert_close_fit(actual, expected, single):
    assert actual["converged"] and expected["converged"]
    assert abs(actual["iterations"] - expected["iterations"]) <= 1
    tolerance = 1e-7 if single else 1e-12
    for key in ("beta", "resid", "prediction"):
        if key in expected:
            np.testing.assert_allclose(actual[key], expected[key], rtol=0,
                                       atol=tolerance * np.abs(expected[key]).max())


def _spy(monkeypatch, module=_vb):
    calls = []
    original = module._correct_products
    monkeypatch.setattr(module, "_correct_products", lambda *args: calls.append(1) or original(*args))
    return calls


# (col_fold, col_group): the layouts of tests/test_vb_column_layout.py, with
# five folds, and HRATT's strong-structure refit (folds x groups).
COLUMNS = {"cv": (np.zeros(4, dtype=np.int64), np.full(4, -1)),
           "loco": (np.full(3, -1), np.arange(3)),
           "mixed": (np.array([2, 0, -1, 1, -1, 2, 0]), np.array([-1, 1, 2, 0, 0, -1, 2])),
           "crossfit": (np.arange(5), np.full(5, -1)),
           "refit": (np.repeat(np.arange(5), 3), np.tile(np.arange(3), 5))}


@needs_numba
@pytest.mark.parametrize("variant", ["covariates", "intercept", "weighted", "weighted_intercept",
                                     "float64", "packed", "uncached", "partial"])
@pytest.mark.parametrize("columns", list(COLUMNS))
def test_compiled_fit_matches_numpy_corrections(panel, monkeypatch, columns, variant):
    options = {"intercept": {"q": 1}, "weighted": {"weighted": True},
               "weighted_intercept": {"weighted": True, "q": 1}, "float64": {"dtype": np.float64},
               "packed": {"packed": True}}.get(variant, {})
    lg, y = _genotypes(panel, **options)
    col_fold, col_group = COLUMNS[columns]
    # Partial: room for the first sub-block's Grams (every slot the fit reads).
    first = (np.unique(col_fold[col_fold >= 0]).size + 1) * (128 * 129 // 2 * lg.dtype.itemsize + 128 * 5 * 8)
    cache = {"uncached": 0, "partial": 1.5 * first}.get(variant, 1e9)
    calls = _spy(monkeypatch)
    actual = _fit(_vb.VBEngine(lg, folds=panel[4], gram_cache_bytes=cache), y, col_fold, col_group)
    assert calls
    monkeypatch.setattr(_vb, "_CORRECTION_LOOP_MAX", 0)
    calls.clear()
    engine = _vb.VBEngine(lg, folds=panel[4], gram_cache_bytes=cache)
    expected = _fit(engine, y, col_fold, col_group)
    assert not calls
    if variant == "partial":
        assert 0 < len(engine._grams) < sum(1 for _ in lg.slices(128))
    if variant == "intercept":
        # No sampling weights: the corrections vanish to rounding, and
        # these fits are unchanged.
        for key in ("beta", "resid", "prediction", "rel_change"):
            np.testing.assert_array_equal(actual[key], expected[key])
    _assert_close_fit(actual, expected, options.get("dtype", np.float32) == np.float32)


@needs_numba
@pytest.mark.parametrize("refit", [False, True])
@pytest.mark.parametrize("budget", [1e9, 3e5])
def test_compiled_corrections_after_joint_preparation(panel, monkeypatch, refit, budget):
    # HRATT's engines store no all-sample Gram slot; the hold-out fit reads
    # slot 0 as its Gram row and coefs[1] as its held-out coefficients.
    lg, y = _genotypes(panel, weighted=True)
    cv_folds = np.ones(lg.n, dtype=np.int64)
    cv_folds[np.random.default_rng(29).choice(lg.n, lg.n // 10, replace=False)] = 0

    def fits():
        cv = _vb.VBEngine(lg, folds=cv_folds, gram_cache_bytes=budget)
        crossfit = _vb.VBEngine(lg, folds=panel[4], gram_cache_bytes=budget)
        _vb._prepare_crossfit_grams(cv, crossfit)
        G = lg.n_groups
        layout = (np.repeat(np.arange(5), G), np.tile(np.arange(G), 5)) if refit else COLUMNS["crossfit"]
        return (_fit(cv, y, np.zeros(10, dtype=np.int64), np.full(10, -1)),
                _fit(crossfit, y, *layout))

    calls = _spy(monkeypatch)
    actual = fits()
    assert calls
    monkeypatch.setattr(_vb, "_CORRECTION_LOOP_MAX", 0)
    for got, want in zip(actual, fits()):
        _assert_close_fit(got, want, True)


@needs_numba
@pytest.mark.parametrize("columns", ["mixed", "crossfit", "refit"])
def test_threaded_fit_runs_the_same_compiled_corrections(panel, monkeypatch, columns):
    numba = pytest.importorskip("numba")
    if numba.config.NUMBA_NUM_THREADS < 4:
        pytest.skip("requires NUMBA_NUM_THREADS >= 4 before Python starts")
    lg, y = _genotypes(panel, weighted=True)
    col_fold, col_group = COLUMNS[columns]
    calls = _spy(monkeypatch)
    expected = _fit(_vb.VBEngine(lg, folds=panel[4]), y, col_fold, col_group)
    actual = _fit(_vb.VBEngine(lg, folds=panel[4], n_threads=4), y, col_fold, col_group)
    assert calls
    # Serial corrections around parallel sweeps: the serial fit's bits.
    assert actual["iterations"] == expected["iterations"]
    for key in ("beta", "resid", "prediction", "rel_change"):
        np.testing.assert_array_equal(actual[key], expected[key])
    monkeypatch.setattr(_vb, "_CORRECTION_LOOP_MAX", 0)
    _assert_close_fit(actual, _fit(_vb.VBEngine(lg, folds=panel[4], n_threads=4), y, col_fold, col_group),
                      True)


@needs_numba
@pytest.mark.parametrize("limit", [5, 128])
@pytest.mark.parametrize("columns, q", [("cv", 5), ("cv", 1), ("crossfit", 5), ("refit", 5),
                                        ("mixed", 5), ("loco", 0)])
def test_corrections_compile_when_numpy_calls_per_set_cost_more(panel, monkeypatch, columns, q, limit):
    lg, y = _genotypes(panel, q=q)
    col_fold, col_group = COLUMNS[columns]
    monkeypatch.setattr(_vb, "_CORRECTION_LOOP_MAX", limit)
    calls = _spy(monkeypatch)
    _fit(_vb.VBEngine(lg, folds=panel[4]), y, col_fold, col_group, _compute_prediction=False)
    sets = np.unique(col_fold).size
    # One column set per held-out fold (or none): compiled while q x columns
    # <= limit x sets, so always for a column per fold.
    assert bool(calls) == (0 < q * col_fold.size <= limit * sets)


def test_numpy_fallback_keeps_the_column_set_products(panel, monkeypatch, numpy_fallback):
    fallback = numpy_fallback(_vb)
    lg, y = _genotypes(panel, weighted=True)
    monkeypatch.setattr(fallback, "_CORRECTION_LOOP_MAX", 10**9)
    calls = _spy(monkeypatch, fallback)
    col_fold, col_group = COLUMNS["crossfit"]
    actual = _fit(fallback.VBEngine(lg, folds=panel[4]), y, col_fold, col_group, _compute_prediction=False)
    assert not calls and actual["converged"]
