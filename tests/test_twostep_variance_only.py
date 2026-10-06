"""Two-step variance estimation does not complete an unused GLS null fit."""

import numpy as np
import pytest

from mixmogam.genotypes import Genotypes
from mixmogam.lmm import LMM
from mixmogam.twostep import _KOp, _setup, fit_variance_components


def setup(packed):
    rng = np.random.default_rng(852)
    dosage = rng.binomial(2, rng.uniform(.1, .9, 72), size=(48, 72)).astype(np.int8)
    gt = Genotypes(dosage, chromosome=np.repeat(np.arange(3), 24), packed=packed)
    X = rng.normal(size=(48, 1))
    y = .4 * dosage[:, 5] + .5 * X[:, 0] + rng.normal(size=48)
    return _setup(y, gt, X, 25, 19)


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("weighted", [False, True])
def test_operator_trace_is_lazy_exact_and_cached(monkeypatch, packed, weighted):
    st = setup(packed)
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
        # The prepared float64 projected norms, weighted exactly.
        expected_trace = float(weights @ st.lg.zz) / float(weights.sum())
        dense_weights = weights
    divisor = st.lg.m if weights is None else weights.sum()
    dense_K = (Z.astype(np.float64).T * dense_weights) @ Z.astype(np.float64) / divisor
    # The dense matrix uses float32 storage rows; the trace exact norms.
    assert expected_trace == pytest.approx(np.trace(dense_K), rel=1e-6)
    # Every genotype pass of the operator decodes through raw_slices.
    original_slices = st.lg.raw_slices
    passes = []

    def observed_slices(*args, **kwargs):
        passes.append(True)
        yield from original_slices(*args, **kwargs)

    monkeypatch.setattr(st.lg, "raw_slices", observed_slices)
    op = _KOp(st.lg, weights)
    assert not passes  # Construction does not scan the weighted genotypes.
    P = rng.normal(size=(st.lg.n, 2))
    # Both routes cast the projected operand to float32: one more rounding.
    np.testing.assert_allclose(op.matmul(P), dense_K @ P, rtol=2e-5, atol=1e-6)
    assert len(passes) == 1
    assert op._trace is None  # Products do not implicitly ask for a trace.
    assert op.trace == expected_trace
    assert len(passes) == 1  # Prepared norms: the trace needs no genotype pass.

    def forbid_pass(*args, **kwargs):
        raise AssertionError("a cached trace must not scan genotypes again")

    monkeypatch.setattr(st.lg, "raw_slices", forbid_pass)
    assert op.trace == expected_trace


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("weighted", [False, True])
def test_two_step_variances_match_complete_fit_without_gls_or_trace(monkeypatch, packed, weighted):
    st = setup(packed)
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
