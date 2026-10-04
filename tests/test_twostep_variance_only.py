"""Two-step variance estimation does not complete an unused GLS null fit."""

import numpy as np
import pytest

from mixmogam.genotypes import Genotypes
from mixmogam.lmm import LMM
from mixmogam.twostep import _KOp, _setup, fit_variance_components


def setup(cache_bytes):
    rng = np.random.default_rng(852)
    dosage = rng.binomial(2, rng.uniform(.1, .9, 72), size=(48, 72)).astype(np.int8)
    gt = Genotypes(dosage, chromosome=np.repeat(np.arange(3), 24))
    X = rng.normal(size=(48, 1))
    y = .4 * dosage[:, 5] + .5 * X[:, 0] + rng.normal(size=48)
    return _setup(y, gt, X, 25, 19, cache_bytes)


@pytest.mark.parametrize("cache_bytes", [0, 1e9])
@pytest.mark.parametrize("weighted", [False, True])
def test_operator_trace_is_lazy_exact_and_cached(monkeypatch, cache_bytes, weighted):
    st = setup(cache_bytes)
    rng = np.random.default_rng(946)
    weights = rng.uniform(.2, 1.5, st.lg.m) if weighted else None
    blocks = list(st.lg.blocks())
    Z = np.empty((st.lg.m, st.lg.n), dtype=st.lg.dtype)
    for idx, _, block in blocks:
        Z[idx] = block
    if weights is None:
        expected_trace = st.lg.trace
        dense_weights = np.ones(st.lg.m)
    else:
        # Preserve the original eager operator's storage-precision weights,
        # float64 accumulation, and per-block summation order exactly.
        # An explicit loop: Python's sum() compensates float rounding (3.12+).
        expected_trace = 0.0
        for idx, _, block in blocks:
            expected_trace += float(np.einsum("ij,ij,i->", block, block,
                                              weights[idx].astype(block.dtype), dtype=np.float64))
        expected_trace /= float(weights.sum())
        dense_weights = weights.astype(st.lg.dtype).astype(np.float64)
    divisor = st.lg.m if weights is None else weights.sum()
    dense_K = (Z.astype(np.float64).T * dense_weights) @ Z.astype(np.float64) / divisor
    assert expected_trace == pytest.approx(np.trace(dense_K), rel=3e-15)
    original_blocks = st.lg.blocks
    original_products = st.lg._streamed_products
    passes = []

    def observed_blocks(reuse=False):
        passes.append(True)
        yield from original_blocks(reuse=reuse)

    def observed_products(*args):
        passes.append(True)
        return original_products(*args)

    monkeypatch.setattr(st.lg, "blocks", observed_blocks)
    monkeypatch.setattr(st.lg, "_streamed_products", observed_products)
    op = _KOp(st.lg, weights)
    assert not passes  # Construction does not scan the weighted genotypes.
    P = rng.normal(size=(st.lg.n, 2))
    # Streamed products cast the projected operand to float32: one more rounding.
    np.testing.assert_allclose(op.matmul(P), dense_K @ P, rtol=2e-5,
                               atol=1e-6 if cache_bytes == 0 else 5e-7)
    assert len(passes) == 1
    assert op._trace is None  # Products do not implicitly ask for a trace.
    assert op.trace == expected_trace
    assert len(passes) == (2 if weighted else 1)

    def forbid_pass(reuse=False):
        raise AssertionError("a cached trace must not scan genotypes again")

    monkeypatch.setattr(st.lg, "blocks", forbid_pass)
    assert op.trace == expected_trace


@pytest.mark.parametrize("cache_bytes", [0, 1e9])
@pytest.mark.parametrize("weighted", [False, True])
def test_two_step_variances_match_complete_fit_without_gls_or_trace(monkeypatch, cache_bytes, weighted):
    st = setup(cache_bytes)
    weights = np.linspace(.2, 1.5, st.lg.m) if weighted else None
    options = dict(ngrids=12, slq_probes=4, slq_steps=24, slq_deflate=4)
    complete = LMM(st.y, X=st.X, K=_KOp(st.lg, weights), add_intercept=False,
                   random_state=19).fit(solver="slq", **options)

    def forbid_gls(*args, **kwargs):
        raise AssertionError("two-step variance estimation must not compute unused GLS coefficients")

    def forbid_trace(*args):
        raise AssertionError("SLQ variance estimation must not scan a weighted kinship trace")

    monkeypatch.setattr(LMM, "_gls_beta_cg", forbid_gls)
    monkeypatch.setattr(_KOp, "trace", property(forbid_trace))
    variance = fit_variance_components(st, weights=weights, random_state=19, **options)
    for name in ("method", "delta", "vg", "ve", "ll", "pseudo_heritability",
                 "n_grid", "newton_used", "solver"):
        assert getattr(variance, name) == getattr(complete, name)
    assert not hasattr(variance, "beta")
