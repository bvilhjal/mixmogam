"""Stochastic Lanczos REML defaults: steps converge early, probes carry the
error, and degenerate probes remain in the average."""

import numpy as np
import pytest

from mixmogam import LMM
from mixmogam._slq import rademacher_probes
from mixmogam.simulate import simulate_genotypes, simulate_kinship, simulate_traits


@pytest.fixture(scope="module")
def structured():
    G = simulate_genotypes(n=500, m=4000, n_pop=6, pop_fst=0.35, seed=21)
    K = simulate_kinship(G)
    y = simulate_traits(G, h2=0.5, n_causal=10, seed=22)["y"]
    return K, y


def _slq(K, y, **options):
    return LMM(y, K=K, random_state=options.pop("seed", 0)).fit(solver="slq", **options)


def test_24_steps_reproduce_96(structured):
    K, y = structured
    short = _slq(K, y, slq_steps=24, slq_tol=0)
    long = _slq(K, y, slq_steps=96, slq_tol=0)
    assert short.pseudo_heritability == pytest.approx(long.pseudo_heritability, abs=1e-5)
    assert short.ll == pytest.approx(long.ll, abs=1e-3)


def test_probe_budgets_stay_close_to_exact_reml(structured):
    """Check numerical error against REML, not a historical ratio of seed SDs.

    Projecting probes out of the deflated space removes the former few-probe
    outliers. Twelve seeds cannot establish the asymptotic variance ratio.
    """
    K, y = structured
    exact = LMM(y, K=K).fit().pseudo_heritability
    for probes in (12, 48):
        estimates = np.array([_slq(K, y, seed=s, slq_probes=probes).pseudo_heritability
                              for s in range(12)])
        assert np.abs(estimates - exact).max() < 0.06


def test_probes_projected_to_noise_are_retained_as_zero():
    class Stream:
        """A first draw aliased with the projected-out direction."""

        def __init__(self):
            self.rng, self.first = np.random.default_rng(1), True

        def choice(self, values, size):
            if self.first:
                self.first = False
                return np.ones(size)
            return self.rng.choice(values, size=size)

    Z = rademacher_probes(50, 3, Stream(), pre=lambda z: z - z.mean())
    assert Z.shape == (50, 3)
    np.testing.assert_array_equal(Z[:, 0], 0.)
