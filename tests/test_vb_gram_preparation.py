"""Gram preparation: packed storage, compiled gathers and HRATT's joint pass.

The engine stores each Gram as a packed upper triangle and widens its rows
into reused buffers; neither changes a stored value. HRATT's hold-out and
cross-fitted engines are prepared together (one pass, one budget); the
cross-fitted Grams are those of a separate engine, the hold-out Grams differ
from a separate engine's only by summation order, and neither fit depends
on the budget.
"""

import importlib.util
import warnings

import numpy as np
import pytest

from mixmogam import _fast, _vb
from mixmogam._loco import LocoGenotypes
from mixmogam._vb import PRIOR_ENET, VBEngine
from mixmogam.genotypes import Genotypes


def _packed(grams):
    rows, cols = np.triu_indices(grams.shape[-1])
    return np.ascontiguousarray(grams[:, rows, cols])


def _numpy_gram(eng, Zs, idx, folds, used=None):
    """The Gram arithmetic before compiled gathers and packing: NumPy's
    fancy gathers and fresh float64 copies, dense (folds + 1, k, k)."""
    c, Q = eng.lg._projection[idx], eng.lg.Q
    q, k = Q.shape[1], Zs.shape[0]
    G = np.empty((folds.size + 1, k, k), dtype=eng.lg.dtype)
    coefs = np.empty((folds.size + 1, k, q))
    held = [None] * folds.size
    if folds.size > 1 and folds.size == eng.n_folds:
        order, bounds, Qp, QtQ = eng._fold_layout()
        Z64 = Zs[:, order].astype(np.float64)
        full = np.zeros((k, k))
        for i, f in enumerate(folds):
            rows = slice(bounds[f], bounds[f + 1])
            Zt = Z64[:, rows]
            held[i] = Zt @ Zt.T
            full += held[i]
            if q:
                coefs[i + 1] = Zt @ Qp[rows]
        if q:
            coefs[0] = coefs[1:].sum(axis=0)
    else:
        Z64 = Zs.astype(np.float64)
        full = Z64 @ Z64.T
        if q:
            coefs[0] = Z64 @ Q
    if q:
        cross = coefs[0] @ c.T
        full -= cross + cross.T - c @ c.T
    G[0] = full
    for i, f in enumerate(folds):
        if used is not None and i + 1 not in used:
            continue
        H = held[i]
        if H is None:
            mask = eng.folds == f
            Zt = Z64[:, mask]
            H = Zt @ Zt.T
            Qt = Q[mask]
            if q:
                coefs[i + 1] = Zt @ Qt
            QtQ_f = Qt.T @ Qt
        else:
            QtQ_f = QtQ[f]
        if q:
            cross = coefs[i + 1] @ c.T
            H -= cross + cross.T - c @ QtQ_f @ c.T
        G[i + 1] = full - H
    return G, coefs


def _problem(dtype=np.float32, q=2, weighted=False, packed=False, n=91, m=411, seed=883):
    rng = np.random.default_rng(seed)
    G = rng.binomial(2, rng.uniform(0.05, 0.9, m), size=(n, m)).astype(np.int8)
    G[rng.random(G.shape) < 0.02] = -1
    G[:, 5] = 1
    groups = np.arange(m) % 3  # noncontiguous variants per group
    gt = Genotypes(G, chromosome=groups, packed=packed)
    scale = np.exp(rng.normal(0, 0.4, n)) if weighted else None
    X = np.column_stack([np.ones(n), rng.standard_normal((n, max(q, 1) - 1))])[:, :q]
    Q = np.linalg.qr(X * (1.0 if scale is None else scale[:, None]))[0] if q else None
    return LocoGenotypes(gt, groups, Q, block=97, dtype=dtype, row_scale=scale)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("q, weighted, packed", [(0, False, False), (1, False, True),
                                                 (2, True, False), (7, False, False)])
def test_gram_matches_numpy_gathers_exactly(dtype, q, weighted, packed):
    lg = _problem(dtype, q, weighted, packed)
    folds = np.random.default_rng(5).permutation(lg.n) % 4
    eng = VBEngine(lg, folds=folds, sub_block=40)
    sizes = set()
    for idx, _, Z in eng._subblocks():
        for k in sorted({1, 2, 17, idx.size}):
            if k > idx.size:
                continue
            sizes.add(k)
            Zk, ik = Z[:k], idx[:k]
            for request in (np.arange(4), np.array([1]), np.array([0, 3])):
                expected, expected_coefs = _numpy_gram(eng, Zk, ik, request)
                grams, coefs = eng._gram(Zk, ik, request)
                assert grams.dtype == lg.dtype and grams.shape == (request.size + 1, k * (k + 1) // 2)
                np.testing.assert_array_equal(grams, _packed(expected))
                np.testing.assert_array_equal(coefs, expected_coefs)
                used = {request.size}  # only the last fold's slot
                grams, coefs = eng._gram(Zk, ik, request, used)
                expected, expected_coefs = _numpy_gram(eng, Zk, ik, request, used)
                np.testing.assert_array_equal(grams[request.size], _packed(expected)[request.size])
                np.testing.assert_array_equal(coefs[[0, request.size]], expected_coefs[[0, request.size]])
    assert sizes == {1, 2, 17, 40}  # 17: the tail of a 97-variant block
    # Exactness rests on NumPy's gathers being F-ordered, which the reused
    # buffers copy; pin that layout.
    Z = next(eng._subblocks())[2]
    order = eng._fold_layout()[0]
    assert Z[:, order].astype(np.float64).flags.f_contiguous
    assert Z.astype(np.float64)[:, folds == 1].flags.f_contiguous
    eng._prepare_grams(np.array([0, 2]))
    assert not eng._scratch  # gather buffers end with the preparation


def _pair(lg, budget, n_threads=1, warn=False, module=_vb):
    rng = np.random.default_rng(29)
    cv_folds = np.ones(lg.n, dtype=np.int64)
    cv_folds[rng.choice(lg.n, lg.n // 10, replace=False)] = 0
    folds = rng.permutation(lg.n) % 5
    cv = module.VBEngine(lg, folds=cv_folds, gram_cache_bytes=budget, n_threads=n_threads,
                         _warn_on_partial=warn)
    crossfit = module.VBEngine(lg, folds=folds, gram_cache_bytes=budget, n_threads=n_threads,
                               _warn_on_partial=warn)
    module._prepare_crossfit_grams(cv, crossfit)
    return cv, crossfit


def _sizes(lg, item):
    """Per sub-block bytes of the hold-out and cross-fitted caches."""
    q = lg.Q.shape[1]
    k = np.array([idx.size for idx, _ in lg.slices(128)])
    tri = k * (k + 1) // 2 * item
    return tri + 2 * k * q * 8, 5 * tri + 6 * k * q * 8


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("q, weighted", [(0, False), (1, False), (3, True)])
def test_joint_preparation_matches_separate_engines(dtype, q, weighted):
    lg = _problem(dtype, q, weighted, n=137, m=611)
    cv, crossfit = _pair(lg, 1e9)
    assert not cv._scratch and not crossfit._scratch  # gather buffers end with the pass
    separate_cv = VBEngine(lg, folds=cv.folds)
    separate = VBEngine(lg, folds=crossfit.folds)
    eps = np.finfo(dtype).eps
    assert len(cv._grams) == len(crossfit._grams) == len(list(lg.slices(128)))
    for b, (idx, _, Z) in enumerate(crossfit._subblocks()):
        k = idx.size
        # No all-sample slot: one hold-out Gram, five cross-fitted ones.
        assert cv._grams[b].shape == (1, k * (k + 1) // 2)
        assert crossfit._grams[b].shape == (5, k * (k + 1) // 2)
        grams, coefs = separate._gram(Z, idx, np.arange(5))
        np.testing.assert_array_equal(crossfit._grams[b], grams[1:])
        np.testing.assert_array_equal(crossfit._coefs[b], coefs)
        # The hold-out Gram sums the products in another order.
        grams, coefs = separate_cv._gram(Z, idx, np.array([0]))
        np.testing.assert_allclose(cv._grams[b], grams[1:], rtol=0,
                                   atol=64 * eps * np.abs(grams).max())
        np.testing.assert_array_equal(cv._coefs[b][1], coefs[1])
        np.testing.assert_allclose(cv._coefs[b][0], coefs[0], rtol=0,
                                   atol=1e-12 * max(np.abs(coefs).max(initial=0), 1))
        # A sub-block rebuilt during a fit equals the cached one.
        for eng, request in ((cv, np.array([0])), (crossfit, np.arange(5))):
            grams, coefs = eng._gram(Z, idx, request)
            np.testing.assert_array_equal(grams, eng._grams[b])
            np.testing.assert_array_equal(coefs, eng._coefs[b])


def test_shared_budget_caches_every_hold_out_gram_first():
    lg = _problem(q=1, n=67, m=1100)
    cv_size, fit_size = _sizes(lg, 4)
    total = cv_size.size
    for budget, n_cv, n_fit in ((0, 0, 0),
                                (cv_size[:3].sum(), 3, 0),
                                (cv_size.sum() - 1, total - 1, 0),
                                (cv_size.sum() + fit_size[:2].sum() + 1, total, 2),
                                (cv_size.sum() + fit_size.sum(), total, total)):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cv, crossfit = _pair(lg, budget, warn=True)
        assert len(cv._grams or []) == n_cv and len(crossfit._grams or []) == n_fit
        partial = [w for w in caught if "Gram cache holds" in str(w.message)]
        if 0 < budget and n_fit < total:
            assert len(partial) == 1
            assert (f"holds {n_cv} of {total} sub-blocks for the hold-out fit and {n_fit} "
                    "for the cross-fitted fits") in str(partial[0].message)
            assert "shared by both" in str(partial[0].message)
        else:
            assert not partial


@pytest.mark.parametrize("packed, q", [(False, 0), (False, 2), (True, 2)])
def test_crossfit_extends_partial_cache_on_first_sweep(packed, q, monkeypatch):
    lg = _problem(q=q, packed=packed, n=67, m=1100)
    cv_size, fit_size = _sizes(lg, lg.dtype.itemsize)
    budget = int(cv_size.sum() + fit_size[:2].sum() + 1)
    expected_blocks = int(np.searchsorted(np.cumsum(fit_size), budget, side="right"))
    assert 2 < expected_blocks < fit_size.size
    rng = np.random.default_rng(711)
    y = rng.normal(size=(lg.n, 1))
    fits, builds = [], []
    for extend in (False, True):
        cv, eng = _pair(lg, budget)
        assert len(eng._grams) == 2
        del cv  # its share of the budget is now available
        counts = []
        original = eng._gram

        def counted(*args, **kwargs):
            counts.append(1)
            return original(*args, **kwargs)

        monkeypatch.setattr(eng, "_gram", counted)
        # Only two folds are used. Newly cached blocks must nevertheless
        # contain every slot for a subsequent fit using the other folds.
        folds = np.array([0, 1])
        fit = eng.fit(np.repeat(y, folds.size, axis=1), folds, np.full(folds.size, -1),
                      PRIOR_ENET, np.tile([0.2, 1.0, 0.3], (folds.size, 1)),
                      np.full(folds.size, 0.7), snp_scale=np.full(lg.m, 0.5 / lg.m),
                      max_iter=3, tol=0, _extend_grams=extend)
        fits.append(fit)
        builds.append(len(counts))
        assert len(eng._grams) == (expected_blocks if extend else 2)
        assert sum(g.nbytes + c.nbytes for g, c in zip(eng._grams, eng._coefs)) <= budget
        if extend:
            for b, (idx, _, z) in enumerate(eng._subblocks()):
                if b >= len(eng._grams):
                    break
                grams, coefs = original(z, idx, np.arange(5))
                np.testing.assert_array_equal(eng._grams[b], grams)
                np.testing.assert_array_equal(eng._coefs[b], coefs)
    _assert_same(fits[0], fits[1])
    assert builds[0] - builds[1] == 2 * (expected_blocks - 2)


def _fits(cv, crossfit, Y, refit=False):
    P = 4
    prior = np.column_stack([np.linspace(0.1, 0.9, P), np.full(P, 30.0), np.full(P, 0.4 / cv.lg.m)])
    scale = np.linspace(0.2, 1.4, cv.lg.m)
    first = cv.fit(np.repeat(Y[:, :1], P, axis=1), np.zeros(P, dtype=np.int64), np.full(P, -1),
                   PRIOR_ENET, prior, np.full(P, 0.8), snp_scale=scale, max_iter=60, tol=1e-9)
    groups = cv.lg.n_groups
    col_fold = np.repeat(np.arange(5), groups) if refit else np.arange(5)
    col_group = np.tile(np.arange(groups), 5) if refit else np.full(5, -1)
    second = crossfit.fit(np.repeat(Y[:, :1], col_fold.size, axis=1), col_fold, col_group,
                          PRIOR_ENET, np.tile(prior[1], (col_fold.size, 1)),
                          np.full(col_fold.size, 0.8), snp_scale=scale, max_iter=60, tol=1e-9)
    return first, second


def _assert_same(actual, expected):
    assert actual["iterations"] == expected["iterations"]
    for key in ("beta", "rel_change", "prediction"):
        np.testing.assert_array_equal(actual[key], expected[key])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("refit", [False, True])
def test_fits_do_not_depend_on_the_budget(dtype, refit):
    lg = _problem(dtype, q=2, n=97, m=700)
    Y = np.random.default_rng(3).standard_normal((lg.n, 1))
    cv_size, fit_size = _sizes(lg, np.dtype(dtype).itemsize)
    expected = _fits(*_pair(lg, 1e9), Y, refit)
    assert expected[0]["converged"] and expected[1]["converged"]
    # Neither cache complete, the hold-out cache complete with part of the
    # cross-fit's, and none: rebuilt sub-blocks use the cached arithmetic.
    for budget in (cv_size[:2].sum(), cv_size.sum() + fit_size[:3].sum(), 0):
        for actual, want in zip(_fits(*_pair(lg, budget), Y, refit), expected):
            _assert_same(actual, want)
    # The cross-fitted fits are those of a separately prepared engine.
    folds = _pair(lg, 0)[1].folds
    _assert_same(_fits(_pair(lg, 1e9)[0], VBEngine(lg, folds=folds), Y, refit)[1], expected[1])


@pytest.mark.parametrize("n_threads", [2, 4])
def test_parallel_sweeps_read_the_packed_slots(n_threads):
    numba = pytest.importorskip("numba")
    if numba.config.NUMBA_NUM_THREADS < n_threads:
        pytest.skip(f"requires NUMBA_NUM_THREADS >= {n_threads} before Python starts")
    lg = _problem(q=2, n=97, m=700)
    Y = np.random.default_rng(3).standard_normal((lg.n, 1))
    cv_size, _ = _sizes(lg, 4)
    for budget in (1e9, cv_size[:2].sum()):
        expected = _fits(*_pair(lg, budget), Y, refit=True)
        for actual, want in zip(_fits(*_pair(lg, budget, n_threads=n_threads), Y, refit=True), expected):
            _assert_same(actual, want)


def test_numpy_fallback_prepares_the_same_fits(monkeypatch):
    with monkeypatch.context() as context:
        context.setattr(_fast, "HAS_NUMBA", False)
        spec = importlib.util.spec_from_file_location("mixmogam._vb_gram_fallback", _vb.__file__)
        fallback = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(fallback)
    lg = _problem(q=2, n=61, m=300)
    Y = np.random.default_rng(3).standard_normal((lg.n, 1))
    cv_size, _ = _sizes(lg, 4)
    cv, crossfit = _pair(lg, 1e9, module=fallback)
    reference = _pair(lg, 1e9)
    for eng, compiled in zip((cv, crossfit), reference):
        for grams, want in zip(eng._grams, compiled._grams):
            np.testing.assert_array_equal(grams, want)
    expected = _fits(cv, crossfit, Y)
    for actual, want in zip(_fits(*_pair(lg, cv_size[:2].sum(), module=fallback), Y), expected):
        _assert_same(actual, want)


def test_full_data_fit_after_joint_preparation_stores_the_full_slot():
    lg = _problem(q=2, n=97, m=500)
    Y = np.random.default_rng(4).standard_normal((lg.n, 3))
    cv, crossfit = _pair(lg, 1e9)
    args = (Y, np.array([-1, 2, -1]), np.array([0, -1, 1]), PRIOR_ENET,
            np.tile([0.4, 30.0, 0.4 / lg.m], (3, 1)), np.full(3, 0.8))
    expected = VBEngine(lg, folds=crossfit.folds).fit(*args, max_iter=60, tol=1e-9)
    actual = crossfit.fit(*args, max_iter=60, tol=1e-9)
    assert crossfit._full_slot and crossfit._grams[0].shape[0] == 2
    _assert_same(actual, expected)
