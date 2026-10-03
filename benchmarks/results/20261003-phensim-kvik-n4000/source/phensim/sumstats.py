"""Summary-statistic simulation from genotype-level truth.

The RSS (regression-on-summary-statistics) oracle shared across the
family's benchmark suites: given per-variant effects ``beta`` and
block-diagonal population LD, marginal GWAS effects are
``bhat = R beta + N(0, R / n)``. Plus the effect-size architectures the
draws start from, a marginal GWAS scan over individual-level genotypes,
and a reference-panel LD-noise generator.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple, Union

import numpy as np

from phensim._numba import _jit

__all__ = [
    "simulate_effects",
    "simulate_effects_pair",
    "genetic_correlation",
    "simulate_sumstats",
    "simulate_sumstats_pair",
    "gwas_scan",
    "shake_ld",
    "prepare_blocks",
]

#: Blocks are ``(R, ix)`` pairs: a correlation matrix and the variant
#: indices it covers. Functions consume them in sequence, drawing one
#: block of RNG noise per entry in block order.
Block = Tuple[np.ndarray, np.ndarray]


def _input_eps(R) -> float:
    """Unit roundoff of ``R``'s own floating dtype (float64 otherwise)."""
    dtype = np.asarray(R).dtype
    return float(np.finfo(dtype if np.issubdtype(dtype, np.floating) else np.float64).eps)


def _as_blocks(blocks: Sequence[Block]) -> list:
    """LDpred3-style coverage and tiled dense-correlation checks.

    Preserve input order (and hence RNG order); never assemble genome-wide
    dense LD. Dense floating correlations are required, not encoded D8/LR8.
    PSD is checked by the factorization used by each consumer. Tolerances
    follow the input precision: float32 LD is judged at float32 roundoff.
    """
    normalized = []
    for item in blocks:
        if not isinstance(item, (tuple, list)) or len(item) != 2:
            raise ValueError("each block must be an (R, ix) pair")
        R, ix = map(np.asarray, item)
        if ix.ndim != 1 or ix.size == 0 or not np.issubdtype(ix.dtype, np.integer):
            raise ValueError("LD block indices must be non-empty integer vectors")
        if np.any(ix < 0) or np.unique(ix).size != ix.size:
            raise ValueError("LD block indices must be unique and nonnegative")
        if R.shape != (ix.size, ix.size):
            raise ValueError("each block must be (R, ix) with R len(ix) x len(ix)")
        if not np.issubdtype(R.dtype, np.number) or np.iscomplexobj(R):
            raise ValueError("LD must contain real numeric correlations")
        rtol = max(1e-7, 4 * _input_eps(R))
        exact_symmetry = True
        for start in range(0, ix.size, 256):
            band, transpose = R[start:start + 256], R[:, start:start + 256].T
            if not np.isfinite(band).all():
                raise ValueError("LD correlations must be finite")
            if not np.allclose(band, transpose, rtol=rtol, atol=1e-10):
                raise ValueError("LD correlations must be symmetric")
            if np.any(band < -1 - rtol) or np.any(band > 1 + rtol):
                raise ValueError("LD correlations must lie in [-1, 1]; decode encoded LD first")
            exact_symmetry &= np.array_equal(band, transpose)
        if not np.allclose(np.diag(R), 1.0, rtol=rtol, atol=1e-10):
            raise ValueError("LD correlations must have a unit diagonal")
        if not exact_symmetry:
            # Canonicalize roundoff only, so signal and noise use the same R;
            # floating LD keeps its precision for the PSD tolerance.
            if not np.issubdtype(R.dtype, np.floating):
                R = R.astype(np.float64)
            R = (R + R.T) * 0.5
        normalized.append((R, ix))
    m = sum(ix.size for _, ix in normalized)
    if m == 0:
        raise ValueError("LD blocks must cover at least one variant")
    # Allocate by the number of supplied indices, not their maximum: a bad
    # index must not cause a huge allocation before it can be rejected.
    seen = np.zeros(m, dtype=bool)
    for _, ix in normalized:
        if np.any(ix >= m):
            raise ValueError("LD blocks must cover every index in 0..m-1 exactly once")
        if np.any(seen[ix]):
            raise ValueError("LD block indices must not overlap or repeat")
        seen[ix] = True
    if not seen.all():
        raise ValueError("LD blocks must cover every index in 0..m-1 exactly once")
    return normalized


def _chol(R: np.ndarray) -> np.ndarray:
    """A covariance factor: Cholesky for PD, eigenfactor for singular PSD.

    Negative eigenvalues beyond floating-point roundoff are errors: float64
    algorithmic error, or the input's own rounding (``k * eps``) when it is
    stored at lower precision. Never add diagonal noise: a singular LD
    block has a genuine null space.
    """
    eps_in = _input_eps(R)
    R = np.asarray(R, dtype=np.float64)
    try:
        return np.linalg.cholesky(R)
    except np.linalg.LinAlgError:
        values, vectors = np.linalg.eigh(R)
        tolerance = R.shape[0] * max(
            64 * np.finfo(float).eps * max(1.0, values[-1]), eps_in)
        if values[0] < -tolerance:
            raise ValueError("LD correlations must be positive semidefinite") from None
        return vectors * np.sqrt(np.maximum(values, 0.0))


@dataclass(frozen=True)
class _PreparedBlocks:
    """A validated, factored snapshot of LD blocks, ready for reuse.

    ``entries`` holds ``(R, ix, factor)`` triples of read-only copies --
    later edits to the caller's originals cannot reach them. Treat the
    object as immutable. ``m`` is the covered variant count.
    """
    entries: tuple
    m: int

    def __iter__(self):
        return iter(self.entries)

    def __len__(self):
        return len(self.entries)

    def __getitem__(self, item):
        return self.entries[item]


def _consumer_blocks(blocks: Sequence[Block]):
    """``(entries, m)`` for a block consumer.

    A prepared object yields its factored triples; raw ``(R, ix)`` pairs
    are validated as usual and yielded ``(R, ix, None)`` without copying
    so the caller factorizes one block at a time.
    """
    if isinstance(blocks, _PreparedBlocks):
        return blocks.entries, blocks.m
    normalized = _as_blocks(blocks)
    entries = [(R, ix, None) for R, ix in normalized]
    return entries, sum(ix.size for _, ix in normalized)


def prepare_blocks(blocks: Sequence[Block]) -> _PreparedBlocks:
    """Validate LD blocks once and snapshot them with their factors.

    The same ``_as_blocks`` coverage/shape checks and ``_chol`` PSD
    factorization every block consumer performs, done up front. All four
    block consumers accept the result and skip their own validation and
    factorization; raw ``(R, ix)`` lists keep working on every call.
    Passing an already-prepared object returns it unchanged.
    """
    if isinstance(blocks, _PreparedBlocks):
        return blocks
    entries = []
    for R, ix in _as_blocks(blocks):
        factor = _chol(R)                          # at the input precision
        R = np.array(R, dtype=np.float64)          # snapshot, not a view
        ix = np.array(ix)                          # same integer dtype
        R.flags.writeable = False
        ix.flags.writeable = False
        factor.flags.writeable = False
        entries.append((R, ix, factor))
    m = sum(ix.size for _, ix, _ in entries)
    return _PreparedBlocks(tuple(entries), m)


def _effects_vector(beta, m):
    beta = np.asarray(beta, dtype=float)
    if beta.shape != (m,) or not np.isfinite(beta).all():
        raise ValueError("beta must be a finite vector matching the LD blocks")
    return beta


def _sample_size(n, m):
    n = np.asarray(n, dtype=float)
    if n.shape not in ((), (m,)) or not np.isfinite(n).all() or np.any(n <= 0):
        raise ValueError("n must be positive and finite, scalar or one value per variant")
    return n


def simulate_effects(
    blocks: Sequence[Block],
    h2: float = 0.5,
    n_causal: Optional[int] = None,
    architecture: str = "sparse",
    maf: Optional[np.ndarray] = None,
    alpha: float = -0.3,
    seed: Union[int, np.random.Generator, None] = 0,
) -> np.ndarray:
    """Effect sizes with ``beta' R beta = h2`` over the block-diagonal LD.

    ``architecture`` draws the effect *shape*; the result is then rescaled
    so the population genetic variance under ``blocks`` hits ``h2``
    exactly:

    - ``'sparse'``: ``n_causal`` random normal effects (required);
    - ``'polygenic'``: every variant;
    - ``'maf'``: every variant, with per-allele effect variance
      proportional to ``[2 f (1-f)]^alpha`` -- standardized effects scaled
      by ``[2 f (1-f)]^((1+alpha)/2)`` (needs ``maf``). ``alpha = -1`` is
      flat on the standardized scale; this is ldpred3's ``alpha`` (SBayesS
      ``S``) convention;
    - ``'equal'``: ``n_causal`` random-sign, equal-magnitude effects.

    ``maf`` is a per-variant array indexed like the blocks.
    """
    if architecture not in ("sparse", "polygenic", "maf", "equal"):
        raise ValueError(f"unknown architecture {architecture!r}")
    if not 0 <= h2 <= 1:
        raise ValueError("h2 must be in [0, 1]")
    entries, m = _consumer_blocks(blocks)
    if architecture in ("sparse", "equal") and n_causal is None:
        raise ValueError(f"architecture {architecture!r} needs n_causal")
    if architecture == "maf" and maf is None:
        raise ValueError("architecture 'maf' needs per-variant maf")
    # This consumer does not draw LD noise, but its variance claim still
    # requires valid PSD blocks. Use the same check as the noise samplers;
    # prepared factors already prove them and are not retained here.
    for R, _ix, factor in entries:
        if factor is None:
            _chol(R)
    rng = np.random.default_rng(seed)
    beta = np.zeros(m)
    if architecture in ("sparse", "equal"):
        if (isinstance(n_causal, (bool, np.bool_)) or not isinstance(n_causal, (int, np.integer))
                or n_causal < 0):
            raise ValueError("n_causal must be a nonnegative integer")
        causal = rng.choice(m, size=min(int(n_causal), m), replace=False)
        if architecture == "equal":
            beta[causal] = np.sign(rng.standard_normal(causal.size))
        else:
            beta[causal] = rng.standard_normal(causal.size)
    elif architecture == "polygenic":
        beta = rng.standard_normal(m)
    else:
        f = np.asarray(maf, dtype=float)
        if f.shape != (m,) or not np.isfinite(f).all() or np.any((f <= 0) | (f >= 1)) or not np.isfinite(alpha):
            raise ValueError("maf must be a finite length-m vector in (0, 1), and alpha finite")
        beta = rng.standard_normal(m) * (2.0 * f * (1.0 - f)) ** ((1.0 + alpha) / 2.0)
    if not np.isfinite(beta).all():
        raise ValueError("effect architecture produced non-finite effects")
    if h2 == 0:
        return np.zeros(m)
    var = sum(beta[ix] @ (np.asarray(R, np.float64) @ beta[ix])
              for R, ix, _ in entries)
    if var <= 0:
        raise ValueError("zero genetic variance; check n_causal / architecture")
    return beta * np.sqrt(h2 / var)


def _quadratic(a, b, entries) -> float:
    """``a' R b`` summed over the block-diagonal LD."""
    return sum(a[ix] @ (np.asarray(R, np.float64) @ b[ix]) for R, ix, _ in entries)


def _count(name, value) -> int:
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer))
            or value < 0):
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def simulate_effects_pair(
    blocks: Sequence[Block],
    h2_a: float = 0.5,
    h2_b: float = 0.5,
    rho: float = 0.5,
    *,
    p: Optional[float] = None,
    n_causal=None,
    n_shared: Optional[int] = None,
    seed: Union[int, np.random.Generator, None] = 0,
):
    """Two traits' standardized effects with correlated shared effects.

    Shared causal effects are drawn ``N(0, [[1, rho], [rho, 1]])``; each
    trait is then scaled so ``beta' R beta`` hits its ``h2`` under the
    block LD. Pass exactly one causal layout:

    - ``p``: every variant is causal for both traits with probability
      ``p`` (at least one is forced); ``p=1`` makes all variants causal
      without a draw. The genetic correlation target is ``rho``.
    - ``n_causal`` (an int, or ``(n_a, n_b)``) with ``n_shared``: exact
      counts, ``n_shared`` variants causal for both and the rest for one
      trait only (the MiXeR four-state truth). The target is
      ``rho * n_shared / sqrt(n_a n_b)`` under equal per-variant variance.

    Finite draws under LD realize a different genetic correlation; score
    estimates against :func:`genetic_correlation`. The draw orders are
    those of bipred's benchmarks (``rg_architectures.sim_effects`` and
    ``mixer_overlap._sim_mixture``), which these reproduce bit for bit.
    Returns ``(beta_a, beta_b)``.
    """
    for name, h2 in (("h2_a", h2_a), ("h2_b", h2_b)):
        if not 0 <= h2 <= 1:
            raise ValueError(f"{name} must be in [0, 1]")
    rho = float(rho)
    if not -1.0 <= rho <= 1.0:
        raise ValueError("rho must be in [-1, 1]")
    if (p is None) == (n_causal is None):
        raise ValueError("pass exactly one of p or n_causal")
    entries, m = _consumer_blocks(blocks)
    for R, _ix, factor in entries:
        if factor is None:
            _chol(R)  # the h2 scaling needs valid PSD blocks
    if abs(rho) < 1:
        L = np.linalg.cholesky([[1.0, rho], [rho, 1.0]])
    else:
        L = np.array([[1.0, 0.0], [rho, 0.0]])
    rng = np.random.default_rng(seed)
    beta_a, beta_b = np.zeros(m), np.zeros(m)
    if p is not None:
        p = float(p)
        if not 0.0 < p <= 1.0:
            raise ValueError("p must be in (0, 1]")
        causal = np.ones(m, bool) if p == 1.0 else rng.random(m) < p
        if not causal.any():
            causal[rng.integers(m)] = True
        raw = L @ rng.standard_normal((2, int(causal.sum())))
        beta_a[causal], beta_b[causal] = raw[0], raw[1]
    else:
        n_a, n_b = ((n_causal, n_causal) if np.ndim(n_causal) == 0
                    else tuple(n_causal))
        n_a, n_b = _count("n_causal", n_a), _count("n_causal", n_b)
        n_shared = _count("n_shared", 0 if n_shared is None else n_shared)
        if n_shared > min(n_a, n_b) or n_a + n_b - n_shared > m:
            raise ValueError("n_shared must not exceed either count, nor the union m")
        picks = rng.choice(m, n_a + n_b - n_shared, replace=False)
        shared = picks[:n_shared]
        only_a = picks[n_shared:n_a]
        only_b = picks[n_a:]
        beta_a[only_a] = rng.standard_normal(n_a - n_shared)
        beta_b[only_b] = rng.standard_normal(n_b - n_shared)
        if n_shared:
            raw = L @ rng.standard_normal((2, n_shared))
            beta_a[shared], beta_b[shared] = raw[0], raw[1]
    out = []
    for beta, h2 in ((beta_a, h2_a), (beta_b, h2_b)):
        if h2 == 0:
            out.append(np.zeros(m))
            continue
        var = _quadratic(beta, beta, entries)
        if not var > 0:
            raise ValueError("zero genetic variance; check the causal layout")
        out.append(beta * np.sqrt(h2 / var))
    return out[0], out[1]


def genetic_correlation(beta_a, beta_b, blocks: Sequence[Block]) -> float:
    """Realized genetic correlation ``a'Rb / sqrt(a'Ra b'Rb)`` under the LD.

    NaN when either genetic variance is not positive.
    """
    entries, m = _consumer_blocks(blocks)
    beta_a = _effects_vector(beta_a, m)
    beta_b = _effects_vector(beta_b, m)
    v_a = float(_quadratic(beta_a, beta_a, entries))
    v_b = float(_quadratic(beta_b, beta_b, entries))
    if not (v_a > 0 and v_b > 0):
        return float("nan")
    return float(_quadratic(beta_a, beta_b, entries)) / np.sqrt(v_a * v_b)


def _factor_options(entries, jitter, factors) -> tuple:
    """``(factors, jitter)`` validated; a ``None`` factor derives from the LD."""
    try:
        jitter = float(jitter)
    except (TypeError, ValueError):
        raise ValueError("jitter must be a finite nonnegative scalar") from None
    if not np.isfinite(jitter) or jitter < 0:
        raise ValueError("jitter must be a finite nonnegative scalar")
    if factors is None:
        return [None] * len(entries), jitter
    if jitter:
        raise ValueError("pass jitter or factors, not both")
    factors = [np.asarray(F, dtype=np.float64) for F in factors]
    if len(factors) != len(entries):
        raise ValueError("factors must hold one matrix per LD block")
    for F, (_R, ix, _f) in zip(factors, entries):
        if F.ndim != 2 or F.shape[0] != len(ix) or F.shape[1] < 1 or not np.isfinite(F).all():
            raise ValueError("each factor must be a finite (len(ix), r) matrix")
    return factors, jitter


def _noise_factor(R, factor, jitter, supplied):
    """The block's noise factor: supplied, ``chol(R + jitter I)``, or the
    validated PSD factor of ``R`` (prepared or computed here)."""
    if supplied is not None:
        return supplied
    if jitter:
        try:
            return np.linalg.cholesky(np.asarray(R, np.float64) + jitter * np.eye(R.shape[0]))
        except np.linalg.LinAlgError:
            raise ValueError("LD + jitter * I must be positive definite") from None
    return _chol(R) if factor is None else factor


def simulate_sumstats(
    beta: np.ndarray,
    blocks: Sequence[Block],
    n,
    seed: Union[int, np.random.Generator, None] = 0,
    *,
    jitter: float = 0.0,
    factors: Optional[Sequence[np.ndarray]] = None,
) -> np.ndarray:
    """Marginal effects from the LDpred model: ``R beta + N(0, R / n)``.

    ``n`` is the GWAS sample size -- a scalar, or a per-variant vector
    indexed by the same ``ix`` as the blocks (heterogeneous N). One RNG
    draw per block, in block order. Blocks must tile 0..m-1 exactly once
    and contain finite, symmetric, unit-diagonal PSD correlations. For
    heterogeneous N the noise covariance is D R D, D_jj = 1/sqrt(n_j);
    this is an oracle model, not a model of arbitrary sample missingness.

    The noise factor is the exact PSD factor of ``R`` by default.
    ``jitter`` uses ``chol(R + jitter I)`` instead, and ``factors`` takes
    one caller-made ``(k, r)`` factor per block (drawing ``r`` normals),
    e.g. an eigenvalue-clipped root of thresholded LD; then ``R`` need
    not be PSD. The signal always uses ``R`` itself. These reproduce the
    family's benchmark draws (``chol(R + 1e-4 I)``, ...) bit for bit.
    """
    entries, m = _consumer_blocks(blocks)
    supplied, jitter = _factor_options(entries, jitter, factors)
    beta = _effects_vector(beta, m)
    n = _sample_size(n, m)
    rng = np.random.default_rng(seed)
    bhat = np.empty(m)
    per_variant = np.ndim(n) > 0
    for (R, ix, factor), F in zip(entries, supplied):
        factor = _noise_factor(R, factor, jitter, F)
        R = np.asarray(R, np.float64)
        noise = factor @ rng.standard_normal(factor.shape[1])
        bhat[ix] = (
            R @ beta[ix]
            + noise / np.sqrt(n[ix] if per_variant else n)
        )
    return bhat


def simulate_sumstats_pair(
    beta_a: np.ndarray,
    beta_b: np.ndarray,
    blocks: Sequence[Block],
    n,
    noise_correlation: Optional[float] = None,
    seed: Union[int, np.random.Generator, None] = 0,
    *,
    overlap: Optional[float] = None,
    n_b=None,
    jitter: float = 0.0,
    factors: Optional[Sequence[np.ndarray]] = None,
):
    """Two GWAS marginal-effect vectors with correlated sampling noise.

    ``noise_correlation`` is rho in ``Cov(noise_a, noise_b) = rho R/n``
    (default 0), or rho D R D for per-variant N. It is not the fraction
    of shared participants: for equal-size studies with independent
    standardized residuals, even complete overlap gives rho=0. Under
    the conditional RSS model, rho is overlap fraction times residual
    correlation. ``overlap`` is a deprecated spelling for the historical
    noise correlation; it warns rather than silently reinterpreting old
    calls. ``n`` is trait A's sample size and ``n_b`` trait B's (default
    ``n``); with unequal sizes the noise covariance is
    ``rho R / sqrt(n_a n_b)``. ``jitter`` and ``factors`` choose the noise
    factor as in :func:`simulate_sumstats`. Per block, ``z_a`` then
    ``z_b`` are drawn. Returns ``(bhat_a, bhat_b)``.
    """
    if overlap is not None:
        if noise_correlation is not None:
            raise ValueError("pass noise_correlation, not both noise_correlation and overlap")
        warnings.warn("overlap means noise correlation, not participant overlap; use "
                      "noise_correlation explicitly", FutureWarning, stacklevel=2)
        noise_correlation = overlap
    rho = 0.0 if noise_correlation is None else float(noise_correlation)
    if not -1.0 <= rho <= 1.0:
        raise ValueError("noise_correlation must be in [-1, 1]")
    entries, m = _consumer_blocks(blocks)
    supplied, jitter = _factor_options(entries, jitter, factors)
    rng = np.random.default_rng(seed)
    beta_a = _effects_vector(beta_a, m)
    beta_b = _effects_vector(beta_b, m)
    n = _sample_size(n, m)
    n_b = n if n_b is None else _sample_size(n_b, m)
    bhat_a = np.empty(m)
    bhat_b = np.empty(m)
    scale = np.sqrt(1.0 - rho**2)
    for (R, ix, factor), F in zip(entries, supplied):
        chol = _noise_factor(R, factor, jitter, F)
        R = np.asarray(R, np.float64)
        z1 = rng.standard_normal(chol.shape[1])
        z2 = rng.standard_normal(chol.shape[1])
        rootn_a = np.sqrt(n[ix] if np.ndim(n) else n)
        rootn_b = np.sqrt(n_b[ix] if np.ndim(n_b) else n_b)
        bhat_a[ix] = R @ beta_a[ix] + (chol @ z1) / rootn_a
        bhat_b[ix] = R @ beta_b[ix] + (chol @ (rho * z1 + scale * z2)) / rootn_b
    return bhat_a, bhat_b


@_jit
def _normal_pvalues(z: np.ndarray) -> np.ndarray:
    """Per-variant two-sided normal p-values ``erfc(|z| / sqrt(2))``.

    JIT-compiled when Numba is installed; identical values to the plain
    ``math.erfc`` loop on every input, including infinities and NaN.
    """
    out = np.empty(z.shape[0], dtype=np.float64)
    for i in range(z.shape[0]):
        out[i] = math.erfc(abs(z[i]) / math.sqrt(2.0))
    return out


def _called_varies(values: np.ndarray, ok: np.ndarray) -> np.ndarray:
    """Per column: do the called entries take at least two distinct values?"""
    return (np.where(ok, values, np.inf).min(axis=0)
            < np.where(ok, values, -np.inf).max(axis=0))


def gwas_scan(
    G: np.ndarray, y: np.ndarray
) -> dict:
    """Marginal GWAS scan of a phenotype over genotype columns.

    Each variant uses its called samples (missing = NaN or negative).
    Both genotype and phenotype are standardized within that subset, so
    beta is Pearson r and ``se = sqrt((1-r^2)/(n_called-2))``. ``z`` is
    the OLS t statistic, retained under its historical name; ``p`` uses
    the large-sample normal approximation ``erfc(|z|/sqrt(2))``, not an
    exact finite-sample t test. Perfect associations give signed infinity
    and p=0. Untestable variants (<3 calls or a constant genotype/called
    phenotype) have NaN outputs. The phenotype must be finite. Returns
    ``{"beta", "se", "z", "p"}`` on the standardized scale.
    """
    Gd = np.asarray(G, dtype=np.float64)
    y = np.asarray(y, dtype=float)
    if Gd.ndim != 2 or y.ndim != 1 or Gd.shape[0] != y.size:
        raise ValueError("G and y must have the same number of samples")
    if not np.isfinite(y).all() or y.size < 3 or y.std() == 0 or np.isinf(Gd).any():
        raise ValueError("need a finite nonconstant phenotype, at least 3 samples, and no infinite genotypes")
    n, m = Gd.shape
    miss = (Gd < 0) | np.isnan(Gd)
    ok = ~miss
    cen = np.where(miss, 0.0, Gd)
    cnt = ok.sum(axis=0)
    cen -= cen.sum(axis=0) / np.maximum(cnt, 1)
    cen[miss] = 0.0
    ss_g = np.einsum("ij,ij->j", cen, cen)
    # A global shift improves stability without changing any subset's OLS.
    yc = y - y.mean()
    sum_y = np.einsum("ij,i->j", ok, yc)
    ss_y = np.einsum("ij,i->j", ok, yc * yc) - sum_y**2 / np.maximum(cnt, 1)
    num = cen.T @ yc
    # Constancy is tested exactly on the called values: the subtracted sums
    # above leave roundoff (pseudo-random r) where the data have none.
    valid = (cnt >= 3) & (ss_g > 0) & (ss_y > 0)
    for start in range(0, m, 4096):
        cols = slice(start, start + 4096)
        valid[cols] &= _called_varies(Gd[:, cols], ok[:, cols])
        valid[cols] &= _called_varies(y[:, None], ok[:, cols])
    r, se, z = (np.full(m, np.nan) for _ in range(3))
    r[valid] = np.clip(num[valid] / np.sqrt(ss_g[valid] * ss_y[valid]), -1.0, 1.0)
    se[valid] = np.sqrt(np.maximum(1 - r[valid]**2, 0) / (cnt[valid] - 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        z[valid] = r[valid] / se[valid]
    p = _normal_pvalues(z)
    return {"beta": r, "se": se, "z": z, "p": p}


def _shaken_correlation(factor, n_ref: int, chunk_size: int, rng) -> np.ndarray:
    """Sample correlation of ``X = Z factor'`` accumulated in row chunks.

    Stable one-pass Chan update of the running count, column mean and
    centred second-moment matrix ``M2`` -- ``O(chunk * k + k^2)`` work
    storage instead of the ``n_ref`` rows. Chunk boundaries change only
    floating-point reduction grouping, so results agree with the
    full-panel path to roundoff while the RNG stream is identical.
    """
    k = factor.shape[0]
    count = 0
    mean = np.zeros(k)
    M2 = np.zeros((k, k))
    for start in range(0, n_ref, chunk_size):
        r = min(chunk_size, n_ref - start)
        X = rng.standard_normal((r, k)) @ factor.T
        local_mean = X.mean(0)
        X -= local_mean
        delta = local_mean - mean
        total = count + r
        M2 += X.T @ X + np.outer(delta, delta) * (count * r / total)
        mean += delta * (r / total)
        count = total
    sd = np.sqrt(np.maximum(np.diag(M2), 0.0))
    sd[sd == 0] = 1.0
    return M2 / sd[:, None] / sd[None, :]


def shake_ld(
    blocks: Sequence[Block],
    n_ref: Optional[int],
    seed: Union[int, np.random.Generator, None] = 0,
    *,
    chunk_size: Optional[int] = None,
    shrink: float = 0.0,
    jitter: float = 0.0,
):
    """Reference-panel LD: the truth, or a finite noisy panel of it.

    ``n_ref=None`` returns the blocks symmetrised (the exact population
    LD). Otherwise each block draws a Wishart panel ``X = Z F'`` with
    ``n_ref`` rows, standardizes its columns and returns ``X'X / n_ref``,
    the sample correlation (unit diagonal to roundoff) -- exactly the
    mismatch a finite reference panel hands an LD-based method. ``F`` is
    the PSD factor of ``R``, or ``chol(R + jitter I)``. ``shrink`` returns
    ``(1 - shrink) LD + shrink I``. With ``jitter=1e-4`` and a cast to
    float32 this is the family's benchmark panel (ldpred3/gwfm
    ``panel_genome``, bipred ``ref_panel``) bit for bit.

    ``chunk_size`` opts into accumulating the panel in row chunks of at
    most ``chunk_size`` samples (a centered one-pass update,
    ``O(chunk * k + k^2)`` storage); ``None`` or a value ``>= n_ref`` is
    the full-panel path, and chunked draws consume the same RNG stream
    but can differ at roundoff. Returns a new ``(R, ix)`` list.
    """
    entries, _m = _consumer_blocks(blocks)
    if n_ref is not None and (isinstance(n_ref, (bool, np.bool_))
                              or not isinstance(n_ref, (int, np.integer)) or n_ref < 2):
        raise ValueError("n_ref must be an integer at least 2")
    if chunk_size is not None and (isinstance(chunk_size, (bool, np.bool_))
                                   or not isinstance(chunk_size, (int, np.integer))
                                   or chunk_size < 1):
        raise ValueError("chunk_size must be a positive integer")
    _supplied, jitter = _factor_options(entries, jitter, None)
    try:
        shrink = float(shrink)
    except (TypeError, ValueError):
        raise ValueError("shrink must be in [0, 1]") from None
    if not 0.0 <= shrink <= 1.0:
        raise ValueError("shrink must be in [0, 1]")
    rng = np.random.default_rng(seed)
    out = []
    for R, ix, factor in entries:
        factor = _noise_factor(R, factor, jitter, None)  # also validates R
        R = np.asarray(R, dtype=np.float64)
        if n_ref is None:
            pass  # R is the symmetrised population LD
        elif chunk_size is None or int(chunk_size) >= n_ref:
            X = rng.standard_normal((int(n_ref), len(ix))) @ factor.T
            sd = X.std(0)
            X = (X - X.mean(0)) / np.where(sd > 0, sd, 1.0)
            R = (X.T @ X) / n_ref
        else:
            R = _shaken_correlation(factor, int(n_ref), int(chunk_size), rng)
        if shrink:
            R = (1.0 - shrink) * R + shrink * np.eye(len(ix))
        R = (R + R.T) / 2.0
        out.append((R, ix.copy()))
    return out
