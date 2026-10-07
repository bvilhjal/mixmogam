"""Model-column parallelism preserves the serial arithmetic and fit decisions."""

import importlib.util

import numpy as np
import pytest

from mixmogam import _fast, _vb, twostep
from mixmogam._loco import LocoGenotypes
from mixmogam.genotypes import Genotypes


def _require_threads(n_threads):
    numba = pytest.importorskip("numba")
    if numba.config.NUMBA_NUM_THREADS < n_threads:
        pytest.skip(f"requires NUMBA_NUM_THREADS >= {n_threads} before Python starts")
    return numba


@pytest.mark.parametrize("n_threads", [1, 2, 4])
@pytest.mark.parametrize("prior_type", [_vb.PRIOR_MIXTURE, _vb.PRIOR_ENET])
@pytest.mark.parametrize("layout", ["C", "F", "partial_F"])
def test_parallel_sweep_is_exact_with_skips_zero_scales_and_warm_effects(n_threads, prior_type, layout):
    _require_threads(n_threads)
    rng = np.random.default_rng(737)
    b, p = 19, 7
    Z = rng.standard_normal((3, b, 31))
    grams = Z @ Z.transpose(0, 2, 1)
    grams[:, 5, :] = 0
    grams[:, :, 5] = 0  # nonpositive diagonal leaves the existing effect alone
    gidx = np.array([0, 1, 2, 0, 2, 1, 0])
    skip = np.array([False, True, False, False, False, True, False])
    prior = np.column_stack([np.linspace(0, 1, p),
                             np.full(p, 0.04 if prior_type == 0 else 7.0),
                             np.linspace(0, 0.015, p)])
    scale = np.linspace(0, 1.5, b)
    scale[11] = 0
    s2e = np.linspace(0.4, 1.2, p)
    initial_u = rng.standard_normal((b, p))
    initial_beta = rng.standard_normal((b, p)) * 0.01
    expected = (initial_u.copy(), initial_beta.copy(), np.full((b, p), np.nan))
    # A tail block slices the first axis of the full F-order workspace,
    # producing an arbitrary-stride view instead of an F-contiguous array.
    def workspace(source):
        buffer = np.empty((b + 5 if layout == "partial_F" else b, p),
                          order="C" if layout == "C" else "F")
        view = buffer[:b]
        view[:] = source
        if layout == "partial_F":
            assert not view.flags.c_contiguous and not view.flags.f_contiguous
        return view

    actual = (workspace(expected[0]), expected[1].copy(), workspace(expected[2]))
    _vb._sweep_block(expected[0], expected[1], grams, gidx, skip, prior_type,
                     prior, scale, s2e, expected[2])
    with _vb._numba_thread_limit(n_threads):
        _vb._sweep_block_parallel(actual[0], actual[1], grams, gidx, skip, prior_type,
                                  prior, scale, s2e, actual[2])
    for got, want in zip(actual, expected):
        np.testing.assert_array_equal(got, want)
    assert np.all(actual[2][:, skip] == 0)
    assert np.all(actual[1][[0, 11]][:, ~skip] == 0)


@pytest.fixture
def phensim_problem():
    """Small structured panel: three populations, four traits, plain NumPy."""
    rng = np.random.default_rng(741)
    n, m = 81, 397
    labels = np.arange(n) % 3
    base = rng.uniform(0.1, 0.5, m)
    freq = np.clip(base + rng.normal(0, 0.06, (3, m)), 0.02, 0.98)
    G = rng.binomial(2, freq[labels]).astype(np.int8)
    Z = (G - G.mean(0)) / np.maximum(G.std(0), 1e-8)
    Y = np.empty((n, 4))
    for c in range(4):
        causal = rng.choice(m, 20, replace=False)
        g = Z[:, causal] @ rng.standard_normal(20)
        Y[:, c] = 0.3 ** 0.5 * g / g.std() + 0.7 ** 0.5 * rng.standard_normal(n)
    Q = np.linalg.qr(np.column_stack([np.ones(n), labels == 1, labels == 2]))[0]
    Y -= Q @ (Q.T @ Y)
    # Missing and constant markers exercise the same imputation/projection in
    # both paths; the complete data above generated every phenotype.
    G = G.copy()
    G[::7, ::13] = -1
    G[:, 7] = 0
    groups = np.arange(G.shape[1]) % 3
    gt = Genotypes(G, chromosome=groups)
    folds = np.arange(G.shape[0]) % 3
    return gt, groups, Q, folds, Y


def _fit_engine(engine, Y, col_fold, col_group, prior_type, max_iter=150, zero_scale=False):
    p = Y.shape[1]
    prior = np.tile([0.3, 0.003 if prior_type == 0 else 25.0, 0.0003], (p, 1))
    prior[:, 0] = np.linspace(0, 1, p)
    scale = np.zeros(engine.lg.m) if zero_scale else np.linspace(0, 1.3, engine.lg.m)
    return engine.fit(Y, col_fold, col_group, prior_type, prior, np.linspace(0.7, 1.0, p),
                      snp_scale=scale, max_iter=max_iter, tol=1e-9)


def _assert_fit_equal(actual, expected):
    assert actual["iterations"] == expected["iterations"]
    assert actual["converged"] is expected["converged"]
    for key in ("beta", "resid", "prediction", "rel_change", "mask"):
        if expected[key] is None:
            assert actual[key] is None
        else:
            np.testing.assert_array_equal(actual[key], expected[key])


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("prior_type", [_vb.PRIOR_MIXTURE, _vb.PRIOR_ENET])
@pytest.mark.parametrize("cache_bytes", [0, 10**8])
def test_parallel_fits_match_serial_folds_loco_and_cached_reuse(phensim_problem, dtype, prior_type, cache_bytes):
    numba = _require_threads(4)
    gt, groups, Q, folds, Y = phensim_problem
    lg = LocoGenotypes(gt, groups, Q, block=173, dtype=dtype)
    serial = _vb.VBEngine(lg, folds=folds, gram_cache_bytes=cache_bytes)
    col_fold, col_group = np.array([2, 0, -1, 2]), np.array([-1, 1, 2, 0])
    expected_cv = _fit_engine(serial, Y, col_fold, col_group, prior_type)
    assert expected_cv["converged"]
    expected_loco = _fit_engine(serial, Y, np.full(4, -1), col_group, prior_type)
    original_threads = numba.get_num_threads()
    for n_threads in (2, 4):
        engine = _vb.VBEngine(lg, folds=folds, gram_cache_bytes=cache_bytes, n_threads=n_threads)
        _assert_fit_equal(_fit_engine(engine, Y, col_fold, col_group, prior_type), expected_cv)
        _assert_fit_equal(_fit_engine(engine, Y, np.full(4, -1), col_group, prior_type), expected_loco)
        assert numba.get_num_threads() == original_threads
        assert any(idx.size == 128 for idx, _, _ in engine._subblocks())
        assert any(idx.size < 128 for idx, _, _ in engine._subblocks())
    for c, (fold, group) in enumerate(zip(col_fold, col_group)):
        assert np.all(expected_cv["beta"][groups == group, c] == 0)
        assert np.all(expected_cv["resid"][folds == fold, c] == 0)


def test_parallel_first_sweep_retains_nonconvergence(phensim_problem):
    _require_threads(2)
    gt, groups, Q, folds, Y = phensim_problem
    lg = LocoGenotypes(gt, groups, Q, block=173)
    col_fold, col_group = np.array([2, 0, -1, 2]), np.array([-1, 1, 2, 0])
    expected = _fit_engine(_vb.VBEngine(lg, folds=folds), Y, col_fold, col_group,
                           _vb.PRIOR_ENET, max_iter=1)
    actual = _fit_engine(_vb.VBEngine(lg, folds=folds, n_threads=2), Y, col_fold, col_group,
                         _vb.PRIOR_ENET, max_iter=1)
    assert not expected["converged"] and expected["iterations"] == 1
    _assert_fit_equal(actual, expected)


@pytest.mark.parametrize("prior_type", [_vb.PRIOR_MIXTURE, _vb.PRIOR_ENET])
def test_parallel_zero_genetic_scale_is_point_mass(phensim_problem, prior_type):
    _require_threads(2)
    gt, groups, Q, _, Y = phensim_problem
    lg = LocoGenotypes(gt, groups, Q, block=173)
    for n_threads in (1, 2):
        fit = _fit_engine(_vb.VBEngine(lg, n_threads=n_threads), Y, np.full(4, -1),
                          np.full(4, -1), prior_type, zero_scale=True)
        assert fit["converged"] and fit["iterations"] == 1
        assert np.all(fit["beta"] == 0) and np.all(fit["rel_change"] == 0)
        np.testing.assert_array_equal(fit["resid"], Y)


def test_parallel_failure_restores_callers_thread_mask(phensim_problem, monkeypatch):
    numba = _require_threads(4)
    gt, groups, Q, _, Y = phensim_problem
    engine = _vb.VBEngine(LocoGenotypes(gt, groups, Q, block=173), n_threads=4)

    def interrupted(*args):
        assert numba.get_num_threads() == 4
        raise RuntimeError("interrupted residual update")

    monkeypatch.setattr(_vb, "_refresh_residual", interrupted)
    with _vb._numba_thread_limit(2):
        with pytest.raises(RuntimeError, match="interrupted residual update"):
            _fit_engine(engine, Y, np.full(4, -1), np.full(4, -1), _vb.PRIOR_ENET)
        assert numba.get_num_threads() == 2


def test_single_column_retains_serial_kernels(phensim_problem, monkeypatch):
    numba = _require_threads(2)
    gt, groups, Q, _, Y = phensim_problem
    lg = LocoGenotypes(gt, groups, Q, block=173)

    def unwanted_parallel(*args):
        raise AssertionError("a single column used a parallel kernel")

    monkeypatch.setattr(_vb, "_sweep_block_parallel", unwanted_parallel)
    before = numba.get_num_threads()
    expected = _fit_engine(_vb.VBEngine(lg), Y[:, :1], np.array([-1]), np.array([-1]), _vb.PRIOR_ENET)
    actual = _fit_engine(_vb.VBEngine(lg, n_threads=2), Y[:, :1], np.array([-1]), np.array([-1]), _vb.PRIOR_ENET)
    _assert_fit_equal(actual, expected)
    assert numba.get_num_threads() == before


@pytest.mark.parametrize("n_threads", [0, -1, True, np.bool_(False), 1.0, 2.5, "2", None])
def test_invalid_threads_fail_before_genotype_access(n_threads):
    with pytest.raises(ValueError, match="positive integer"):
        _vb.VBEngine(object(), n_threads=n_threads)
    with pytest.raises(ValueError, match="positive integer"):
        twostep.hratt(None, None, n_threads=n_threads)


def test_numba_thread_limit_is_validated_before_genotype_access():
    numba = _require_threads(1)
    requested = numba.config.NUMBA_NUM_THREADS + 1
    with pytest.raises(ValueError, match="Numba thread limit"):
        _vb.VBEngine(object(), n_threads=requested)
    with pytest.raises(ValueError, match="Numba thread limit"):
        twostep.hratt(None, None, n_threads=requested)
    assert _vb._validate_n_threads(np.int64(1)) == 1


def test_numpy_fallback_only_requires_numba_for_multiple_threads(monkeypatch):
    with monkeypatch.context() as context:
        context.setattr(_fast, "HAS_NUMBA", False)
        spec = importlib.util.spec_from_file_location("mixmogam._vb_parallel_fallback", _vb.__file__)
        fallback = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(fallback)
    assert fallback._validate_n_threads(1) == 1
    with pytest.raises(ImportError, match="requires Numba"):
        fallback.VBEngine(object(), n_threads=2)


def test_hratt_parallel_matches_serial_with_phensim_data(phensim_problem):
    _require_threads(4)
    gt, _, Q, _, Y = phensim_problem
    options = dict(X=Q[:, 1:], heritability_method="he", alphas=(-1.0,),
                   grid=[(0.0, 1.0), (0.3, 0.3), (0.5, 0.5)],
                   n_calibration=8, vb_max_iter=200, random_state=761, block=173)
    expected = twostep.hratt(Y[:, 0], gt, n_threads=1, **options)
    assert expected.extra["cv_converged"] and expected.extra["loco_converged"]
    for n_threads in (2, 4):
        actual = twostep.hratt(Y[:, 0], gt, n_threads=n_threads, **options)
        for key in ("beta", "se", "f_stat", "p"):
            # The missing/constant-variant fixture also changes a few
            # near-zero projected genotypes by float64 reduction rounding.
            np.testing.assert_allclose(getattr(actual, key), getattr(expected, key),
                                       rtol=32 * np.finfo(float).eps,
                                       atol=4 * np.finfo(float).eps)
        assert actual.extra.keys() == expected.extra.keys()
        for key, value in expected.extra.items():
            if isinstance(value, np.ndarray):
                if key == "cv_mse":
                    np.testing.assert_allclose(actual.extra[key], value,
                                               rtol=32 * np.finfo(float).eps,
                                               atol=4 * np.finfo(float).eps)
                else:
                    np.testing.assert_array_equal(actual.extra[key], value)
            else:
                assert actual.extra[key] == value


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("b", [128, 32, 7])
@pytest.mark.parametrize("p", [2, 6])
def test_gemm_workspace_matches_serial_and_float64_reference(dtype, b, p):
    from concurrent.futures import ThreadPoolExecutor
    pytest.importorskip("scipy.linalg.blas")  # load before discovering BLAS pools
    control = pytest.importorskip("threadpoolctl")
    rng = np.random.default_rng(877)
    n = 503
    Z = rng.standard_normal((128, n)).astype(dtype)[:b]
    work = rng.standard_normal((n, p)).astype(dtype)
    changes = rng.standard_normal((128, p)).astype(dtype)[:b]
    expected_forward = Z.astype(float) @ work.astype(float)
    expected_backward = Z.astype(float).T @ changes.astype(float)
    tolerance = 5e-5 if dtype == np.float32 else 2e-12
    with control.threadpool_limits(limits=1, user_api="blas"):
        serial_forward, serial_backward = Z @ work, Z.T @ changes
        for workers in (2, 4):
            with ThreadPoolExecutor(max_workers=workers) as pool:
                workspace = _vb._GemmWorkspace(n, p, 128, dtype, pool, workers)
                output = np.empty((128, p), dtype=dtype)[:b]
                actual_forward = workspace.forward(Z, work, output)
                actual_backward = np.empty((n, p), dtype=dtype)
                workspace.backward(Z, changes, actual_backward)
                if b >= 32:
                    assert actual_forward.flags.f_contiguous
                    assert np.shares_memory(actual_forward, workspace.products)
                else:
                    assert actual_forward is output
                for actual, serial, reference in (
                        (actual_forward, serial_forward, expected_forward),
                        (actual_backward, serial_backward, expected_backward)):
                    # Different BLAS builds may select different arithmetic
                    # for these geometries. Check against higher precision,
                    # without treating the existing float32 result as truth.
                    np.testing.assert_allclose(actual, reference, rtol=tolerance, atol=tolerance)
                    np.testing.assert_allclose(actual, serial, rtol=tolerance, atol=tolerance)
    assert Z.flags.writeable and changes.flags.writeable


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_large_gemm_dispatch_preserves_fit_and_reuses_one_pool(phensim_problem, monkeypatch, dtype):
    _require_threads(4)
    pytest.importorskip("threadpoolctl")
    gt, groups, Q, folds, Y = phensim_problem
    lg = LocoGenotypes(gt, groups, Q, block=173, dtype=dtype)
    col_fold, col_group = np.array([2, 0, -1, 2]), np.array([-1, 1, 2, 0])
    expected = _fit_engine(_vb.VBEngine(lg, folds=folds), Y, col_fold, col_group, _vb.PRIOR_ENET)
    monkeypatch.setattr(_vb, "_GEMM_MIN_SAMPLES", 1)
    original = _vb.ThreadPoolExecutor
    pools = []

    def track_pool(**kwargs):
        pool = original(**kwargs)
        pools.append(pool)
        return pool

    monkeypatch.setattr(_vb, "ThreadPoolExecutor", track_pool)
    actual = _fit_engine(_vb.VBEngine(lg, folds=folds, n_threads=4), Y, col_fold, col_group, _vb.PRIOR_ENET)
    assert len(pools) == 1
    assert expected["converged"] and actual["converged"]
    assert expected["iterations"] == actual["iterations"]
    tolerance = 2e-6 if dtype == np.float32 else 5e-12
    for field in ("beta", "resid", "prediction"):
        np.testing.assert_allclose(actual[field], expected[field], rtol=tolerance, atol=tolerance * 0.01)
    np.testing.assert_allclose(actual["rel_change"], expected["rel_change"], rtol=1e-3, atol=1e-14)
    np.testing.assert_array_equal(actual["mask"], expected["mask"])
    with pytest.raises(RuntimeError, match="shutdown"):
        pools[0].submit(lambda: None)


def test_large_gemm_failure_joins_workers_and_restores_thread_controls(phensim_problem, monkeypatch):
    import threading
    numba = _require_threads(4)
    control = pytest.importorskip("threadpoolctl")
    gt, groups, Q, _, Y = phensim_problem
    engine = _vb.VBEngine(LocoGenotypes(gt, groups, Q, block=173), n_threads=4)
    monkeypatch.setattr(_vb, "_GEMM_MIN_SAMPLES", 1)
    workers = []

    def fail_worker(left, right, out):
        workers.append(threading.current_thread())
        assert not left.flags.writeable and not right.flags.writeable
        raise RuntimeError("injected GEMM worker failure")

    monkeypatch.setattr(_vb, "_matmul_rows", fail_worker)
    with control.threadpool_limits(limits=2, user_api="blas"), _vb._numba_thread_limit(2):
        before = {p["filepath"]: p["num_threads"] for p in control.threadpool_info()}
        with pytest.raises(RuntimeError, match="injected GEMM worker failure"):
            _fit_engine(engine, Y, np.full(4, -1), np.full(4, -1), _vb.PRIOR_ENET)
        assert numba.get_num_threads() == 2
        assert {p["filepath"]: p["num_threads"] for p in control.threadpool_info()} == before
    assert workers and all(not worker.is_alive() for worker in workers)


def test_large_gemm_loads_blas_before_setting_its_limit(phensim_problem, monkeypatch):
    from contextlib import contextmanager
    _require_threads(2)
    control = pytest.importorskip("threadpoolctl")
    gt, groups, Q, _, Y = phensim_problem
    lg = LocoGenotypes(gt, groups, Q, block=173)
    monkeypatch.setattr(_vb, "_GEMM_MIN_SAMPLES", 1)
    initialized = []
    constructor = _vb._GemmWorkspace
    original_limits = control.threadpool_limits

    def workspace(*args):
        result = constructor(*args)
        initialized.append(result)
        return result

    @contextmanager
    def limits(*args, **kwargs):
        assert initialized  # SciPy BLAS resolution precedes pool discovery.
        with original_limits(*args, **kwargs):
            yield

    monkeypatch.setattr(_vb, "_GemmWorkspace", workspace)
    monkeypatch.setattr(control, "threadpool_limits", limits)
    _fit_engine(_vb.VBEngine(lg, n_threads=2), Y, np.full(4, -1), np.full(4, -1),
                _vb.PRIOR_ENET, max_iter=1)


def test_missing_threadpoolctl_retains_serial_gemms(phensim_problem, monkeypatch):
    import builtins
    _require_threads(2)
    gt, groups, Q, _, Y = phensim_problem
    lg = LocoGenotypes(gt, groups, Q, block=173)
    expected = _fit_engine(_vb.VBEngine(lg), Y, np.full(4, -1), np.full(4, -1), _vb.PRIOR_ENET)
    monkeypatch.setattr(_vb, "_GEMM_MIN_SAMPLES", 1)
    original_import = builtins.__import__

    def unavailable(name, *args, **kwargs):
        if name == "threadpoolctl":
            raise ImportError("optional BLAS controller unavailable")
        return original_import(name, *args, **kwargs)

    def no_pool(**kwargs):
        raise AssertionError("pool started without BLAS control")

    monkeypatch.setattr(builtins, "__import__", unavailable)
    monkeypatch.setattr(_vb, "ThreadPoolExecutor", no_pool)
    actual = _fit_engine(_vb.VBEngine(lg, n_threads=2), Y, np.full(4, -1), np.full(4, -1), _vb.PRIOR_ENET)
    _assert_fit_equal(actual, expected)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("parent_size", [67, 173])
def test_parallel_decode_retains_parent_boundaries(phensim_problem, monkeypatch, dtype, parent_size):
    _require_threads(2)
    gt, groups, Q, folds, Y = phensim_problem
    lg = LocoGenotypes(gt, groups, Q, block=parent_size, dtype=dtype, n_threads=2)
    engine = _vb.VBEngine(lg, folds=folds, n_threads=2)
    original = lg._decode
    decoded = []

    def track_decode(idx, **kwargs):
        assert idx.size <= 128
        decoded.append(idx.copy())
        return original(idx, **kwargs)

    with monkeypatch.context() as context:
        context.setattr(lg, "_decode", track_decode)
        # Streamed subblocks share one buffer: copy what outlives an iteration.
        blocks = [(i, g, b.copy()) for i, g, b in engine._subblocks()]
    assert len(blocks) == len(decoded)
    np.testing.assert_array_equal(np.sort(np.concatenate([i for i, _, _ in blocks])), np.arange(lg.m))
    calls = np.asarray(gt.G, dtype=np.float64)
    for (idx, group, block), observed in zip(blocks, decoded):
        np.testing.assert_array_equal(observed, idx)
        # Each subblock lies inside one parent block of its own group.
        assert any(g == group and np.isin(idx, parent).all() for parent, g in lg._blocks)
        # Unprojected standardized calls; no-calls decode to zero.
        g = calls[:, idx]
        sd = np.where(lg.sd[idx] > 0, lg.sd[idx], 1.0)
        expected = np.where(g != -1, (g - lg.mean[idx]) / sd, 0.0).T
        np.testing.assert_allclose(block, expected, rtol=0, atol=8 * np.finfo(dtype).eps * np.abs(expected).max())
