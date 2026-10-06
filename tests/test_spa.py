"""Saddlepoint p-values of retrospective score statistics against exact
distributions."""

import itertools

import numpy as np
import pytest
from scipy import stats

from mixmogam._spa import _cgf, spa_pvalue

CALLS = np.array([0.0, 1.0, 2.0, 0.0])


def hwe(maf):
    return np.array([(1 - maf) ** 2, 2 * maf * (1 - maf), maf**2, 0.0])


def test_cgf_value_slope_and_curvature():
    rng = np.random.default_rng(31)
    a = rng.normal(size=50)
    freqs = np.array([0.81, 0.18, 0.01, 0.0])
    values = CALLS - CALLS @ freqs
    k0, k1, k2 = _cgf(0.0, a, values, freqs)
    assert k0 == pytest.approx(0.0, abs=1e-13) and k1 == pytest.approx(0.0, abs=1e-12)
    assert k2 == pytest.approx((freqs @ values**2) * np.sum(a * a), rel=1e-12)
    h = 1e-5
    for t in (-0.7, 0.3, 1.2):
        _, k1, k2 = _cgf(t, a, values, freqs)
        assert (_cgf(t + h, a, values, freqs)[0] - _cgf(t - h, a, values, freqs)[0]) / (2 * h) == pytest.approx(k1, rel=1e-7)
        assert (_cgf(t + h, a, values, freqs)[1] - _cgf(t - h, a, values, freqs)[1]) / (2 * h) == pytest.approx(k2, rel=1e-7)


def lattice(a, freqs):
    """Exact distribution of sum a_i g_i for integer a_i (any sign) and
    genotypes g_i in {0, 1, 2} with frequencies ``freqs``: the lowest
    attainable value and the probabilities from there on."""
    low = int(np.sum(np.minimum(2 * a, 0)))
    dist = np.zeros(int(np.sum(np.abs(2 * a))) + 1)
    dist[0] = 1.0
    span = 0
    for ai in a.astype(np.int64):
        new = np.zeros_like(dist)
        for g, pg in enumerate(freqs[:3]):
            shift = ai * g - min(2 * ai, 0)
            new[shift : shift + span + 1] += pg * dist[: span + 1]
        dist, span = new, span + abs(2 * ai)
    return low, dist


def lattice_errors(a, freqs):
    """log10 errors of SPA and normal two-sided tails at continuity-corrected
    lattice points with exact tails in [1e-7, 1e-2]."""
    low, dist = lattice(a, freqs)
    upper, lower = np.cumsum(dist[::-1])[::-1], np.cumsum(dist)
    mean = (CALLS @ freqs) * a.sum()
    sd = np.sqrt(np.sum(a * a) * (freqs @ (CALLS - CALLS @ freqs) ** 2))
    spa, normal = [], []
    for k in range(dist.size):
        x = low + k - 0.5 - mean  # continuity-corrected point below lattice value low + k
        if x <= 0.5 * sd:
            continue
        below = int(np.floor(mean - x)) - low
        exact = upper[k] + (lower[below] if below >= 0 else 0.0)
        if not 1e-7 <= exact <= 1e-2:
            continue
        approx = spa_pvalue(np.array([x]), a, CALLS[None, :], freqs[None, :])[0][0]
        spa.append(abs(np.log10(approx / exact)))
        normal.append(abs(np.log10(2 * stats.norm.sf(x / sd) / exact)))
    return np.array(spa), np.array(normal)


def heavy_tailed(n=400, seed=37):
    """Weighted residuals of outcome-dependent selection: a few samples carry
    the score."""
    rng = np.random.default_rng(seed)
    return np.round(np.exp(rng.normal(0.0, 1.0, n)) * rng.choice([-1, 1], n))


def case_control(n, cases):
    """Residuals y - mu of a case-control sample: cases (n - cases) / cases,
    controls -1."""
    return np.r_[np.full(cases, (n - cases) / cases), -np.ones(n - cases)]


# Fewer expected case alleles make the distribution clumpy and every smooth
# approximation coarse; the normal one is off by orders of magnitude.
@pytest.mark.parametrize("a, maf, spa_bound, normal_floor", [
    (heavy_tailed(), 0.02, 0.1, 1.0),
    (case_control(2000, 100), 0.02, 0.02, 3.0),  # 5% cases, 4 expected case alleles
    (case_control(1000, 20), 0.05, 0.2, 3.0),  # 2% cases, 2 expected case alleles
], ids=["heavy-tailed", "cases-5pct", "cases-2pct"])
def test_lattice_tails_match_the_exact_distribution(a, maf, spa_bound, normal_floor):
    spa, normal = lattice_errors(a, hwe(maf))
    assert spa.size >= 10 and spa.max() < spa_bound
    assert normal.max() > normal_floor


def test_two_sided_tails_match_enumeration_with_real_coefficients_and_no_calls():
    rng = np.random.default_rng(38)
    n = 10
    a = rng.normal(size=n) * np.exp(rng.normal(0.0, 0.7, n))
    raw = np.array([0.0, 1.0, 2.0, 0.6])  # three calls and the no-call value
    freqs = np.array([0.7, 0.2, 0.05, 0.05])
    mean = raw @ freqs
    G = np.array(list(itertools.product(range(4), repeat=n)))
    prob = np.prod(freqs[G], axis=1)
    S = (raw[G] - mean) @ a
    sd = np.sqrt(np.sum(a * a) * (freqs @ (raw - mean) ** 2))
    errors = []
    for z in np.linspace(1.5, 4.5, 31):
        exact = prob[S >= z * sd].sum() + prob[S <= -z * sd].sum()
        if not 1e-6 <= exact <= 1e-1:
            continue
        approx = spa_pvalue(np.array([z * sd]), a, raw[None, :], freqs[None, :])[0][0]
        errors.append(abs(np.log10(approx / exact)))
    assert len(errors) >= 15
    assert max(errors) < 0.25 and np.median(errors) < 0.1


def test_symmetry_and_score_scaling():
    rng = np.random.default_rng(33)
    k, n = 6, 200
    a = rng.normal(size=n)
    values = np.tile(CALLS, (k, 1))
    freqs = np.tile(hwe(0.05), (k, 1))
    sd = np.sqrt(np.sum(a * a) * (freqs[0] @ (values[0] - values[0] @ freqs[0]) ** 2))
    u = sd * np.array([2.5, -3.0, 4.0, -5.0, 6.0, 7.5])
    p, log_p = spa_pvalue(u, a, values, freqs)
    assert np.all((p > 0) & (p < 1)) and np.all(np.isfinite(log_p))
    np.testing.assert_allclose(np.exp(log_p), p, rtol=1e-14)
    # Two-sided: negating the score, the coefficients or the genotype
    # values leaves the p-value unchanged.
    np.testing.assert_allclose(spa_pvalue(-u, a, values, freqs)[0], p, rtol=1e-12)
    np.testing.assert_allclose(spa_pvalue(u, -a, values, freqs)[0], p, rtol=1e-9)
    np.testing.assert_allclose(spa_pvalue(u, a, -values, freqs)[0], p, rtol=1e-9)
    # Frequencies are normalized, and var_ratio scales the scores first.
    np.testing.assert_allclose(spa_pvalue(u, a, values, 3 * freqs, var_ratio=0.8)[0],
                               spa_pvalue(u * np.sqrt(0.8), a, values, freqs)[0], rtol=1e-12)
    # The low-frequency allele skews the score: heavier tail on one side.
    assert not np.allclose(p, 2 * stats.norm.sf(np.abs(u) / sd), rtol=0.05)


def test_extreme_scores_and_scores_beyond_the_support():
    rng = np.random.default_rng(35)
    a = rng.uniform(0.5, 1.5, 400) * rng.choice([-1, 1], 400)
    freqs = np.array([0.99, 0.01, 0.0, 0.0])
    sd = np.sqrt(np.sum(a * a) * 0.0099)
    p, log_p = spa_pvalue(np.array([3.0, 8.0]) * sd, a, np.tile(CALLS, (2, 1)), np.tile(freqs, (2, 1)))
    assert np.all(np.isfinite(log_p)) and np.all((p > 0) & (p < 1))
    p, log_p = spa_pvalue(np.array([2.0 * np.sum(np.abs(a))]), a, CALLS[None, :], freqs[None, :])
    assert p[0] == 0.0 and log_p[0] == -np.inf


def test_input_shapes_are_checked():
    with pytest.raises(ValueError, match="equal"):
        spa_pvalue(np.ones(1), np.ones(5), np.zeros((1, 4)), np.ones((1, 3)))
    with pytest.raises(ValueError, match="frequencies"):
        spa_pvalue(np.ones(1), np.ones(5), np.zeros((1, 2)), np.array([[0.5, -0.5]]))


@pytest.mark.numba
def test_parallel_rows_are_exact():
    numba = pytest.importorskip("numba")
    if numba.config.NUMBA_NUM_THREADS < 2:
        pytest.skip("requires NUMBA_NUM_THREADS >= 2 before Python starts")
    rng = np.random.default_rng(36)
    a = rng.normal(size=300)
    values = np.tile(CALLS, (40, 1))
    freqs = np.column_stack([rng.dirichlet([20, 4, 1], 40), np.zeros(40)])
    u = rng.normal(0.0, 4.0, 40) * np.sqrt(np.sum(a * a) * 0.2)
    serial = spa_pvalue(u, a, values, freqs)
    parallel = spa_pvalue(u, a, values, freqs, n_threads=2)
    np.testing.assert_array_equal(parallel[0], serial[0])
    np.testing.assert_array_equal(parallel[1], serial[1])
