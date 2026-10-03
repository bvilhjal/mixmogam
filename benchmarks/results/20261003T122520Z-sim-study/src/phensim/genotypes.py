"""Genotype simulators: independent, structured, haplotype-block and
coalescent backends."""

from __future__ import annotations

from typing import Union

import numpy as np

from phensim._common import norm_isf

__all__ = [
    "simulate_independent",
    "simulate_population_structure",
    "simulate_haplotype_blocks",
    "simulate_ar1_blocks",
    "realistic_block_sizes",
    "simulate_coalescent",
    "simulate_by_mutation_rate",
    "resolve_backend",
]


def _positive_int(name: str, value) -> int:
    """A strictly positive Python/NumPy integer (no bools or floats)."""
    if (isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer)) or value < 1):
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _finite_scalar(name: str, value) -> float:
    """A finite scalar float."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a finite scalar") from None
    if not np.isfinite(v):
        raise ValueError(f"{name} must be a finite scalar")
    return v


def _maf_scalar(value) -> float:
    """A minor allele frequency in [0, 0.5] (endpoints valid)."""
    maf = _finite_scalar("maf", value)
    if not 0.0 <= maf <= 0.5:
        raise ValueError("maf must be a finite value in [0, 0.5]")
    return maf


def _coalescent_params(Ne, recomb_rate, mut_rate, min_maf):
    """Shared scalar validation for the coalescent-family simulators."""
    Ne = _finite_scalar("Ne", Ne)
    recomb_rate = _finite_scalar("recomb_rate", recomb_rate)
    mut_rate = _finite_scalar("mut_rate", mut_rate)
    min_maf = _finite_scalar("min_maf", min_maf)
    if Ne <= 0:
        raise ValueError("Ne must be positive and finite")
    if recomb_rate < 0 or mut_rate < 0:
        raise ValueError("recomb_rate and mut_rate must be finite and nonnegative")
    if not 0 <= min_maf < 0.5:
        raise ValueError("min_maf must be a finite value in [0, 0.5)")
    return Ne, recomb_rate, mut_rate, min_maf


def resolve_backend(backend: str = "auto") -> str:
    """Pick a coalescent backend: 'numba', 'msprime' or 'auto'."""
    if backend not in ("auto", "numba", "msprime"):
        raise ValueError("backend must be 'auto', 'numba' or 'msprime'")
    if backend != "auto":
        return backend
    from phensim._numba import HAVE_NUMBA

    if HAVE_NUMBA:
        return "numba"
    try:
        import msprime  # noqa: F401

        return "msprime"
    except ImportError:
        return "numba"  # pure-Python fallback (correct, just slower)


def _draw_freqs(
    m: int, freq_dist: str, rng: np.random.Generator
) -> np.ndarray:
    """Per-site minor allele frequencies under a named SFS shape."""
    if freq_dist == "beta":
        f = rng.beta(0.35, 1.0, m) * 0.5
    elif freq_dist == "uniform":
        f = rng.uniform(0.01, 0.5, m)
    elif freq_dist == "rare":
        f = rng.exponential(0.02, m).clip(1e-4, 0.5)
    elif freq_dist == "common":
        f = rng.uniform(0.1, 0.5, m)
    else:
        raise ValueError(f"unknown freq_dist {freq_dist!r}")
    return f


def simulate_independent(
    n: int,
    m: int,
    maf: float = 0.3,
    freq_dist: str = "fixed",
    seed: Union[int, np.random.Generator, None] = 0,
) -> np.ndarray:
    """(n, m) int8 dosages with independent SNPs.

    ``freq_dist='fixed'`` uses a single ``maf`` for every site; 'beta',
    'uniform', 'rare' and 'common' draw per-site minor allele frequencies
    from the corresponding allele-frequency spectrum shape. ``maf`` is
    validated only for 'fixed'; it is ignored by the other shapes.
    """
    n = _positive_int("n", n)
    m = _positive_int("m", m)
    rng = np.random.default_rng(seed)
    if freq_dist == "fixed":
        f = np.full(m, _maf_scalar(maf))
    else:
        f = _draw_freqs(m, freq_dist, rng)
    flip = rng.random(m) < 0.5
    p = np.where(flip, 1 - f, f)
    return rng.binomial(2, p[None, :], size=(n, m)).astype(np.int8)


def simulate_population_structure(
    n: int,
    m: int,
    n_pops: int = 3,
    fst: float = 0.1,
    maf: float = 0.3,
    model: str = "normal",
    seed: Union[int, np.random.Generator, None] = 0,
):
    """Independent SNPs with diverged per-population allele frequencies.

    Returns ``(G, pop_labels)``. ``model`` picks the drift model for the
    per-population frequencies around the shared ``maf``-centred base:

    - ``'normal'``: ``p_k ~ N(p, fst * p (1 - p))`` -- the cheap
      approximation, adequate for LMM-confounding studies;
    - ``'balding-nichols'``: ``p_k ~ Beta(p(1-fst)/fst, (1-p)(1-fst)/fst)``,
      the exact Balding--Nichols drift distribution (clipped polymorphic).

    Use :func:`simulate_coalescent` when LD realism matters.
    """
    if model not in ("normal", "balding-nichols"):
        raise ValueError("model must be 'normal' or 'balding-nichols'")
    if not 0.0 < fst < 1.0:
        raise ValueError("fst must be in (0, 1)")
    n = _positive_int("n", n)
    m = _positive_int("m", m)
    n_pops = _positive_int("n_pops", n_pops)
    maf = _maf_scalar(maf)
    rng = np.random.default_rng(seed)
    base = np.clip(maf + rng.normal(0, 0.05, m), 0.05, 0.95)
    if model == "normal":
        spread = np.sqrt(fst * base * (1 - base))
        freqs = np.clip(
            base[:, None] + rng.normal(0, 1, (m, n_pops)) * spread[:, None],
            0.01, 0.99,
        )
    else:
        a = base * (1 - fst) / fst
        b = (1 - base) * (1 - fst) / fst
        freqs = np.clip(
            rng.beta(a[:, None], b[:, None], size=(m, n_pops)), 0.01, 0.99
        )
    labels = rng.integers(0, n_pops, n)
    p = freqs[:, labels].T
    return rng.binomial(2, p).astype(np.int8), labels


def simulate_ar1_blocks(
    n: int,
    block_sizes,
    maf: Union[float, np.ndarray] = 0.3,
    rho: float = 0.9,
    seed: Union[int, np.random.Generator, None] = 0,
    *,
    method: str = "cholesky",
):
    """Dosages with within-block AR(1) LD via a latent Gaussian model.

    Each block of ``k`` SNPs gets two latent Gaussian haplotypes per
    person, ``z ~ N(0, C)`` with ``C_ij = rho**|i-j|``, thresholded at
    the frequency-implied quantile and summed to 0/1/2 dosages. Smooth
    geometric LD decay within blocks, sharp decay between them. ``maf``
    is the counted allele's frequency, a scalar or per-site array in
    ``[0, 1]`` (values above 0.5 count the major allele); ``block_sizes`` a sequence of
    block lengths (see :func:`realistic_block_sizes` for right-skewed
    geometry). Returns ``(G, blocks)`` with ``G`` int8 ``(n, m)`` and
    ``blocks`` column-index arrays.

    ``method='cholesky'`` (default) draws ``z`` through the factorized
    ``C + 1e-8 I``; ``method='scan'`` is an opt-in O(nk) forward
    recursion, ``z_j = rho z_{j-1} + sqrt(1-rho^2) eps_j``, sampling the
    *unjittered* AR(1) covariance ``C`` -- the distributions differ only
    by the default's 1e-8 diagonal jitter, and at ``rho = +/-1`` the
    recursion is exactly deterministic. The RNG call order (two ``(n, k)``
    standard-normal draws per block, blocks in sequence) is fixed under
    both methods; given the same generator and ``maf`` array the Cholesky
    path reproduces the ldpred3 benchmark simulator it was extracted from
    bit for bit. The two methods consume the same draws but produce
    different seeded genotypes.
    """
    if method not in ("cholesky", "scan"):
        raise ValueError("method must be 'cholesky' or 'scan'")
    n = _positive_int("n", n)
    sizes = np.asarray(block_sizes)
    if (sizes.ndim != 1 or sizes.size == 0
            or not np.issubdtype(sizes.dtype, np.integer) or np.any(sizes < 1)):
        raise ValueError(
            "block_sizes must be a non-empty vector of positive integer lengths")
    # Check in the original precision and with Python-int summation: an
    # overflowing length or total must fail, not wrap into a negative
    # allocation or a truncated column count.
    if np.any(sizes > np.iinfo(np.intp).max):
        raise ValueError("block_sizes lengths must fit in memory")
    m = sum(int(k) for k in sizes)
    if m > np.iinfo(np.intp).max:
        raise ValueError("block_sizes total length must fit in memory")
    block_sizes = sizes.astype(np.int64)
    rho = _finite_scalar("rho", rho)
    if not -1.0 <= rho <= 1.0:
        raise ValueError("rho must be in [-1, 1]")
    maf = np.asarray(maf, dtype=float)
    if maf.ndim == 0:
        maf = np.full(m, float(maf))
    if (maf.shape != (m,) or not np.isfinite(maf).all()
            or np.any((maf < 0) | (maf > 1))):
        raise ValueError("maf must be scalar or a finite length-m vector in [0, 1]")
    rng = np.random.default_rng(seed)
    innovation_sd = np.sqrt(1.0 - rho * rho)
    G = np.empty((n, m), dtype=np.int8)
    blocks = []
    col = 0
    for k in block_sizes:
        k = int(k)
        if method == "cholesky":
            idx = np.arange(k)
            corr = rho ** np.abs(idx[:, None] - idx[None, :])
            chol = np.linalg.cholesky(corr + 1e-8 * np.eye(k))
        thr = norm_isf(maf[col:col + k])
        hap_sum = np.zeros((n, k))
        for _ in range(2):  # two haplotypes -> dosage 0/1/2
            eps = rng.standard_normal((n, k))
            if method == "cholesky":
                z = eps @ chol.T
            else:
                z = np.empty((n, k))
                z[:, 0] = eps[:, 0]
                for j in range(1, k):
                    z[:, j] = rho * z[:, j - 1] + innovation_sd * eps[:, j]
            hap_sum += (z > thr)
        G[:, col:col + k] = hap_sum.astype(np.int8)
        blocks.append(np.arange(col, col + k))
        col += k
    return G, blocks


def realistic_block_sizes(m: int, n_blocks: int, *, cv: float = 0.9,
                          seed: Union[int, np.random.Generator, None] = 0):
    """Partition ``m`` SNPs into ``n_blocks`` right-skewed LD blocks.

    Block *lengths* are log-normal with coefficient of variation ``cv`` --
    a tunable synthetic approximation to right-skewed
    recombination-delimited blocks, stressing the few large blocks that
    dominate quadratic LD work. Returns an int array summing to exactly
    ``m`` (rounding drift is repaired at the largest/smallest blocks).
    """
    for name, value in (("m", m), ("n_blocks", n_blocks)):
        _positive_int(name, value)
    if not np.isfinite(cv) or cv < 0:
        raise ValueError("cv must be finite and nonnegative")
    n_blocks = min(int(n_blocks), int(m))
    rng = np.random.default_rng(seed)
    sigma = float(np.sqrt(np.log(1.0 + cv * cv)))  # log-normal CV -> sigma
    w = rng.lognormal(mean=-0.5 * sigma * sigma, sigma=sigma, size=n_blocks)
    sizes = np.maximum(1, np.round(w / w.sum() * m)).astype(np.int64)
    # fix rounding drift so the sizes sum to exactly m
    drift = int(sizes.sum() - m)
    order = np.argsort(sizes)  # adjust largest/smallest
    i = 0
    while drift != 0:
        j = order[-1 - (i % n_blocks)] if drift > 0 else order[i % n_blocks]
        if drift > 0 and sizes[j] > 1:
            sizes[j] -= 1
            drift -= 1
        elif drift < 0:
            sizes[j] += 1
            drift += 1
        i += 1
    return sizes


def simulate_haplotype_blocks(
    n: int,
    m: int,
    block_size: int = 100,
    n_founders: int = 20,
    mutation_rate: float = 0.001,
    seed: Union[int, np.random.Generator, None] = 0,
) -> np.ndarray:
    """LD-structured genotypes from founder haplotype copying.

    Each block of ``block_size`` SNPs is founded by ``n_founders``
    haplotypes; descendants inherit a founder haplotype with per-site
    mutation flips. Fast (one pass, no coalescent) with genuine
    haplotypic LD within blocks and sharp decay between blocks.
    Returns exactly ``m`` columns, including a shorter final block.
    """
    for name, value in (("n", n), ("m", m), ("block_size", block_size), ("n_founders", n_founders)):
        _positive_int(name, value)
    if not 0 <= mutation_rate <= 1:
        raise ValueError("mutation_rate must be in [0, 1]")
    rng = np.random.default_rng(seed)
    G = np.empty((n, m), dtype=np.int8)
    for start in range(0, m, block_size):
        stop = min(start + block_size, m)
        founders = rng.binomial(1, 0.3, size=(n_founders, stop - start))
        parents = rng.integers(0, n_founders, size=(n, 2))
        hap = founders[parents]  # (n, 2, block_size)
        flip = rng.random(hap.shape) < mutation_rate
        hap = np.where(flip, 1 - hap, hap)
        G[:, start:stop] = (
            hap[:, 0, :] + hap[:, 1, :]
        ).astype(np.int8)
    return G


def _coalescent_dosages(n, seq_len, *, recomb_rate, mut_rate, Ne, seed, backend):
    """One coalescent replicate -> ``(dos, af)`` via the chosen backend.

    ``seed`` is ``None`` or an integer in ``[1, 2**31)``, the range both
    backends honour: the built-in kernel masks seeds to 31 bits (so larger
    or negative seeds would alias silently) and msprime rejects 0.
    """
    if seed is not None:
        if (isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer))
                or not 1 <= seed < 2**31):
            raise ValueError("seed must be None or an integer in [1, 2**31)")
        seed = int(seed)
    if backend == "numba":
        from phensim._coalescent import simulate_dosages

        if seed is None:
            seed = int(np.random.default_rng().integers(1, 2**31 - 1))
        dos, _pos, af = simulate_dosages(
            n, seq_len, recomb_rate=recomb_rate, mut_rate=mut_rate, Ne=Ne, seed=seed
        )
        return dos, af
    try:
        import msprime
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "the msprime backend needs msprime (pip install phensim[msprime])"
        ) from e
    ms_seed = None if seed is None else int(seed)
    ts = msprime.sim_ancestry(
        samples=n,
        ploidy=2,
        population_size=Ne,
        recombination_rate=recomb_rate,
        sequence_length=int(seq_len),
        discrete_genome=False,
        random_seed=ms_seed,
    )
    mts = msprime.sim_mutations(
        ts, rate=mut_rate, random_seed=ms_seed, discrete_genome=False,
        model=msprime.BinaryMutationModel()
    )
    # int8 dosages decoded site by site: tskit's genotype_matrix() is int32
    # (sites, 2n) over every site, rare ones included, which with the pair
    # sum peaked near 8 GB at n = 4,000 for m = 50,000 common SNPs
    dos = np.empty((mts.num_sites, n), dtype=np.int8)
    for v in mts.variants(copy=False):
        g = v.genotypes  # 0/1 per haplotype; an individual's two are adjacent
        np.add(g[0::2], g[1::2], out=dos[v.site.id], casting="unsafe")
    dos = dos.T  # (n, sites), 0/1/2
    af = dos.mean(axis=0) / 2.0
    return dos, af


def simulate_coalescent(
    n: int,
    m: int,
    block_size: int = 200,
    *,
    Ne: int = 10_000,
    recomb_rate: float = 1e-8,
    mut_rate: float = 1e-8,
    min_maf: float = 0.01,
    seed: Union[int, None] = None,
    backend: str = "auto",
):
    """Coalescent genotypes with recombination-driven LD.

    Human-like defaults (Ne = 10,000; recombination and mutation rates
    1e-8/bp/generation). The sequence length grows until at least ``m``
    common SNPs (MAF > ``min_maf``) exist; the first ``m`` are kept and
    cut into contiguous blocks of ``block_size``.

    ``backend``: ``'numba'`` (built-in JIT coalescent,
    :mod:`phensim._coalescent`), ``'msprime'`` (the msprime C library) or
    ``'auto'`` (built-in when Numba is available, else msprime, else the
    pure-Python built-in). Returns ``(G, blocks)`` with ``G`` int8
    ``(n, m')`` sample-major dosages and ``blocks`` contiguous index
    arrays; ``m'`` is ``m`` rounded down to a multiple of ``block_size``.
    """
    n = _positive_int("n", n)
    m = _positive_int("m", m)
    block_size = _positive_int("block_size", block_size)
    if block_size > m:
        raise ValueError("block_size must not exceed m")
    Ne, recomb_rate, mut_rate, min_maf = _coalescent_params(
        Ne, recomb_rate, mut_rate, min_maf)
    if mut_rate == 0:
        raise ValueError("mut_rate must be positive to reach a SNP-count target")
    backend = resolve_backend(backend)
    rng = np.random.default_rng(seed)
    seq_len = max(1e6, m / 1200 * 1e6)  # ~1200 common SNPs per Mb to start
    G = None
    for _ in range(7):
        rep_seed = int(rng.integers(1, 2**31 - 1))
        dos, af = _coalescent_dosages(
            n,
            seq_len,
            recomb_rate=recomb_rate,
            mut_rate=mut_rate,
            Ne=Ne,
            seed=rep_seed,
            backend=backend,
        )
        dos = dos[:, (af > min_maf) & (af < 1 - min_maf)]
        if dos.shape[1] >= m:
            G = dos
            break
        seq_len *= 1.8
    if G is None or G.shape[1] < m:
        raise RuntimeError(
            "coalescent simulation produced too few common SNPs; "
            "increase sequence length / Ne"
        )
    n_blocks = m // block_size
    m2 = n_blocks * block_size
    G = np.ascontiguousarray(G[:, :m2].astype(np.int8))
    blocks = [
        np.arange(i * block_size, (i + 1) * block_size) for i in range(n_blocks)
    ]
    return G, blocks


def simulate_by_mutation_rate(
    n: int,
    seq_len: float,
    *,
    recomb_rate: float = 1e-8,
    mut_rate: float = 1e-8,
    Ne: int = 10_000,
    min_maf: float = 0.01,
    seed: Union[int, None] = None,
    backend: str = "auto",
) -> np.ndarray:
    """Coalescent genotypes on a *fixed* segment; density set by mutation.

    Unlike :func:`simulate_coalescent` (which grows the segment to hit a
    SNP-count target, coupling LD extent to SNP count), the segment
    length and recombination rate here fix the LD structure and the
    mutation rate controls how many variants sit on it. With a fixed
    seed the genealogy is identical across mutation rates, so raising
    the rate is the same chromosome with more discovered variants.
    Returns ``G`` int8 ``(n, k)``; ``k`` emerges from the rate. Columns
    are in physical order, so contiguous slices are contiguous LD.
    ``seed`` is ``None`` or an integer in ``[1, 2**31)`` on either backend.
    """
    n = _positive_int("n", n)
    seq_len = _finite_scalar("seq_len", seq_len)
    if seq_len < 1:
        raise ValueError("seq_len must be finite and at least 1")
    Ne, recomb_rate, mut_rate, min_maf = _coalescent_params(
        Ne, recomb_rate, mut_rate, min_maf)
    backend = resolve_backend(backend)
    dos, af = _coalescent_dosages(
        n,
        seq_len,
        recomb_rate=recomb_rate,
        mut_rate=mut_rate,
        Ne=Ne,
        seed=seed,
        backend=backend,
    )
    dos = dos[:, (af > min_maf) & (af < 1 - min_maf)]
    return np.ascontiguousarray(dos.astype(np.int8))
