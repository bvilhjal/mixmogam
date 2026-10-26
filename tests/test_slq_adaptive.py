"""Adaptive Lanczos stopping in stochastic REML (``slq_tol``): a stopped run
is exactly a shorter run, it never stops on a grid end, and ``slq_tol=0``
keeps the fixed step count."""

from types import SimpleNamespace

import numpy as np
import pytest

from mixmogam import LMM, lmm, twostep
from mixmogam._slq import lanczos_quadrature_batch
from mixmogam.genotypes import Genotypes
from mixmogam.simulate import simulate_genotypes, simulate_kinship, simulate_traits


def _rules_equal(a, b):
    assert len(a) == len(b)
    for ra, rb in zip(a, b):
        np.testing.assert_array_equal(ra.theta, rb.theta)
        np.testing.assert_array_equal(ra.weights, rb.weights)
        assert ra.norm2 == rb.norm2


def test_stopped_batch_is_the_shorter_run():
    rng = np.random.default_rng(21)
    A = rng.standard_normal((120, 120))
    A = A @ A.T / 120 + np.diag(np.r_[np.full(3, 20.0), np.zeros(117)])
    V = rng.standard_normal((120, 4))
    V[:, 3] = np.linalg.eigh(A)[1][:, -2:].sum(axis=1)  # breaks down after 2 steps
    seen, offered = [], {}

    def stop(j, rules):
        seen.append(j)
        if j == 7:
            offered[j] = rules()
        return j == 10

    stopped = lanczos_quadrature_batch(lambda X: A @ X, V, 24, stop)
    assert seen == list(range(1, 11))
    assert stopped[3].theta.size == 2
    _rules_equal(stopped, lanczos_quadrature_batch(lambda X: A @ X, V, 10))
    _rules_equal(offered[7], lanczos_quadrature_batch(lambda X: A @ X, V, 7))
    # never asked at the cap: the full run returns its rules
    seen.clear()
    full = lanczos_quadrature_batch(lambda X: A @ X, V, 24, lambda j, rules: seen.append(j))
    assert seen == list(range(1, 24))
    _rules_equal(full, lanczos_quadrature_batch(lambda X: A @ X, V, 24))


@pytest.fixture(scope="module")
def tall():
    """More samples than variants and a spread spectrum (LD-like column
    scales): the short rules put h2 on the grid end for a dozen steps."""
    rng = np.random.default_rng(1)
    n, m = 1000, 250
    scale = np.arange(1, m + 1, dtype=float) ** -2.0
    Z = rng.standard_normal((n, m)) * np.sqrt(scale / scale.mean())
    Z -= Z.mean(axis=0)
    K = Z @ Z.T / m
    g = Z @ rng.standard_normal(m)
    y = np.sqrt(0.5) * g / g.std() + np.sqrt(0.5) * rng.standard_normal(n)
    return K, y


def _counted(K):
    """A streaming kinship that counts its passes (one per Lanczos step)."""
    op = SimpleNamespace(n=K.shape[0], trace=float(np.trace(K)), passes=0)

    def matmul(x):
        op.passes += 1
        return K @ x

    op.matmul = matmul
    return op


def _fit(K, y, **options):
    op = _counted(K)
    fit = LMM(y, K=op, random_state=0)._fit_variance_components(solver="slq", slq_deflate=0,
                                                                **options)
    return fit, op.passes


def test_guarded_stop_waits_for_an_interior_root(tall):
    K, y = tall
    full, passes = _fit(K, y, slq_tol=0)
    assert passes == 24 and full.newton_used
    # An unguarded rule (two changes below the tolerance) stopped here: the
    # 6-, 7- and 8-step fits sit on the grid end, so h2 does not change.
    ends = [_fit(K, y, slq_steps=j, slq_tol=0)[0] for j in (6, 7, 8)]
    assert all(f.pseudo_heritability > 0.9999 and not f.newton_used for f in ends)
    assert len({f.pseudo_heritability for f in ends}) == 1
    assert abs(full.pseudo_heritability - 0.55) < 0.05

    adaptive, stopped = _fit(K, y)
    assert 8 < stopped < 24
    assert adaptive.newton_used
    assert abs(adaptive.pseudo_heritability - full.pseudo_heritability) < 1e-6
    # the stopped run is exactly the shorter fixed run
    shorter, _ = _fit(K, y, slq_steps=stopped, slq_tol=0)
    assert adaptive == shorter


@pytest.fixture(scope="module")
def structured():
    G = simulate_genotypes(n=500, m=4000, n_pop=6, pop_fst=0.35, seed=21)
    return simulate_kinship(G), simulate_traits(G, h2=0.5, n_causal=10, seed=22)["y"]


@pytest.mark.parametrize("seed", [0, 2])
def test_deflated_stop_matches_all_steps(structured, monkeypatch, seed):
    """With deflated outliers the vectorized check fits the same h2 as the
    final scalar evaluation of the same rules."""
    K, y = structured
    kept, checks, passes = [], [], []
    extremes, refine, batch = lmm._separated_extremes, LMM._refine_delta, lmm.lanczos_quadrature_batch

    def spy_extremes(values, vectors, ratio=4.0):
        out = extremes(values, vectors, ratio)
        kept.append(out[0].size)
        return out

    def spy_refine(self, deltas, lls, dlls, ll_at, dll_at, tol):
        out = refine(self, deltas, lls, dlls, ll_at, dll_at, tol)
        checks.append((lls, out))
        return out

    def spy_batch(matvec, V0, steps, stop=None):
        count = []
        rules = batch(lambda X: count.append(1) or matvec(X), V0, steps, stop)
        passes.append(len(count))
        return rules

    full = LMM(y, K=K, random_state=seed).fit(solver="slq", slq_tol=0)
    monkeypatch.setattr(lmm, "_separated_extremes", spy_extremes)
    monkeypatch.setattr(LMM, "_refine_delta", spy_refine)
    monkeypatch.setattr(lmm, "lanczos_quadrature_batch", spy_batch)
    adaptive = LMM(y, K=K, random_state=seed).fit(solver="slq")
    assert kept[-1] > 0  # outliers deflated
    stopped = passes[-1]
    assert stopped < 24 and len(checks) == stopped - 6 + 2  # checks from step 6, then the fit
    assert abs(adaptive.pseudo_heritability - full.pseudo_heritability) < 1e-6
    (last_lls, last), (final_lls, final) = checks[-2:]
    np.testing.assert_allclose(last_lls, final_lls, rtol=1e-12)
    assert last[0] == pytest.approx(final[0], rel=1e-12) and last[2] and final[2]
    monkeypatch.undo()
    shorter = LMM(y, K=K, random_state=seed).fit(solver="slq", slq_tol=0, slq_steps=stopped)
    assert adaptive.pseudo_heritability == shorter.pseudo_heritability


def test_ml_and_zero_tolerance_run_all_steps(tall):
    K, y = tall
    ml, passes = _fit(K, y, method="ml")
    assert passes == 48  # the phenotype run and ML's separate trace run
    assert ml == _fit(K, y, method="ml", slq_tol=0)[0]
    assert _fit(K, y, slq_steps=12, slq_tol=0)[1] == 12


@pytest.mark.parametrize("bad", [-1e-6, np.nan, np.inf])
def test_slq_tol_must_be_nonnegative(tall, bad):
    K, y = tall
    with pytest.raises(ValueError, match="slq_tol"):
        LMM(y, K=K).fit(solver="slq", slq_tol=bad)


def test_public_fit_caches_by_slq_tol(tall):
    K, y = tall
    model = LMM(y, K=K, random_state=0)
    adaptive = model.fit(solver="slq", slq_deflate=0)
    fixed = model.fit(solver="slq", slq_deflate=0, slq_tol=0)
    assert fixed.delta != adaptive.delta and model.fit_result.delta == fixed.delta
    assert abs(adaptive.pseudo_heritability - fixed.pseudo_heritability) < 1e-6


def test_two_step_methods_pass_slq_tol(monkeypatch):
    rng = np.random.default_rng(71)
    n, m = 120, 160
    dosage = rng.binomial(2, rng.uniform(0.1, 0.9, m), size=(n, m)).astype(np.int8)
    gt = Genotypes(dosage, chromosome=np.repeat(np.arange(4), m // 4))
    y = 0.3 * dosage[:, 0] + rng.standard_normal(n)
    seen = []
    original = twostep.fit_variance_components

    def spy(*args, **kwargs):
        seen.append(kwargs.get("slq_tol"))
        return original(*args, **kwargs)

    monkeypatch.setattr(twostep, "fit_variance_components", spy)
    kwargs = dict(alphas=(-1.0, 0.0), grid=[(0.0, 1.0)], he_probes=3, n_calibration=8,
                  loco_folds=1, vb_max_iter=300, random_state=7, block=64)
    twostep.hratt(y, gt, alpha_method="reml", slq_tol=0, **kwargs)
    twostep.hratt(y, gt, **kwargs)
    twostep.bolt_inf(y, gt, n_calibration=8, random_state=2, block=64, slq_tol=0)
    assert seen == [0, 0, 1e-6, 0]
