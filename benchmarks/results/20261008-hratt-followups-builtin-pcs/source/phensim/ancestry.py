"""Several populations and admixed samples, with ancestry truth.

``drift_frequencies`` draws per-population allele frequencies from a shared
ancestral spectrum; ``simulate_populations`` samples discrete populations
with exact sizes and population-specific AR(1) LD; ``simulate_admixed``
builds local-ancestry mosaics over the same model; and
``simulate_split_coalescent`` runs a population split, optionally with an
admixture pulse, through msprime. See docs/technical.md for the models.
"""
from __future__ import annotations

from typing import Union

import numpy as np

from phensim._common import norm_isf
# Sibling caches hash genotypes.py, so this module only reads from it.
from phensim.genotypes import (
    _coalescent_params,
    _finite_scalar,
    _positive_int,
    simulate_ar1_blocks,
)

__all__ = [
    "drift_frequencies",
    "simulate_populations",
    "simulate_admixed",
    "simulate_split_coalescent",
]


def _per_pop(name, value, k, lo, hi, *, hi_open=False):
    """A scalar or length-``k`` vector of finite values in ``[lo, hi]``."""
    try:
        arr = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a scalar or one value per population") from None
    if arr.ndim == 0:
        arr = np.full(k, float(arr))
    upper = arr >= hi if hi_open else arr > hi
    if arr.shape != (k,) or not np.isfinite(arr).all() or np.any(arr < lo) or np.any(upper):
        bound = f"[{lo:g}, {hi:g}{')' if hi_open else ']'}"
        raise ValueError(f"{name} must be a scalar or one value per population in {bound}")
    return arr


def _drift_freqs(base, fst, model, rng):
    """Unclipped ``(m, K)`` frequencies drifted from ``base``, one F per column.

    The same models as ``simulate_population_structure``; with a constant
    F this is its stream. Columns with ``F = 0`` return ``base``.
    """
    zero = fst == 0
    if zero.all():
        return np.repeat(base[:, None], fst.size, axis=1)
    F = np.where(zero, 0.5, fst)[None, :]  # placeholder keeps Beta valid
    p = base[:, None]
    if model == "normal":
        freqs = p + rng.normal(0, 1, (base.size, fst.size)) * np.sqrt(F * p * (1 - p))
    else:
        freqs = rng.beta(p * (1 - F) / F, (1 - p) * (1 - F) / F)
    freqs[:, zero] = base[:, None]
    return freqs


def _sample_sizes(n):
    """One positive integer per population; at least one population."""
    if isinstance(n, (int, np.integer)) or np.ndim(n) != 1 or len(n) == 0:
        raise ValueError("n must be a sequence with one sample size per population")
    return [_positive_int("n", k) for k in n]


def _freq_matrix(freqs, k=None):
    f = np.asarray(freqs, dtype=float)
    if f.ndim != 2 or not f.size or (k is not None and f.shape[0] != k):
        raise ValueError("freqs must have shape (populations, variants)")
    if not np.isfinite(f).all() or np.any((f < 0) | (f > 1)):
        raise ValueError("freqs must be finite frequencies in [0, 1]")
    return f


def _block_vector(sizes, m):
    s = np.asarray(sizes)
    if (s.ndim != 1 or not s.size or not np.issubdtype(s.dtype, np.integer)
            or np.any(s < 1) or sum(int(v) for v in s) != m):
        raise ValueError("block_sizes must be positive integer lengths summing to the variant count")
    return s.astype(np.int64)


def drift_frequencies(
    m: int,
    n_pops: int = 2,
    fst=0.1,
    *,
    ancestral=None,
    model: str = "balding-nichols",
    min_freq: float = 0.01,
    seed: Union[int, np.random.Generator, None] = 0,
) -> np.ndarray:
    """Per-population allele frequencies drifted from one ancestral spectrum.

    Returns ``(n_pops, m)`` float64 frequencies of the counted allele.
    ``ancestral`` gives the shared frequency of each site, in ``(0, 1)``;
    by default it is drawn ``U(0.05, 0.95)``. ``fst`` is a scalar or one
    value per population in ``[0, 1)``, so populations can drift by
    different amounts. ``model='balding-nichols'`` draws the exact
    ``Beta(p(1-F)/F, (1-p)(1-F)/F)``; ``'normal'`` the approximation
    ``N(p, F p (1 - p))``. ``F = 0`` copies the ancestral frequency.
    Every result is clipped to ``[min_freq, 1 - min_freq]``.
    """
    m = _positive_int("m", m)
    n_pops = _positive_int("n_pops", n_pops)
    if model not in ("normal", "balding-nichols"):
        raise ValueError("model must be 'normal' or 'balding-nichols'")
    fst = _per_pop("fst", fst, n_pops, 0.0, 1.0, hi_open=True)
    min_freq = _finite_scalar("min_freq", min_freq)
    if not 0 <= min_freq < 0.5:
        raise ValueError("min_freq must be in [0, 0.5)")
    rng = np.random.default_rng(seed)
    if ancestral is None:
        base = rng.uniform(0.05, 0.95, m)
    else:
        base = np.asarray(ancestral, dtype=float)
        if base.shape != (m,) or not np.isfinite(base).all() or np.any((base <= 0) | (base >= 1)):
            raise ValueError("ancestral must be a length-m vector of frequencies in (0, 1)")
    freqs = _drift_freqs(base, fst, model, rng)
    return np.clip(freqs, min_freq, 1 - min_freq).T


def simulate_populations(
    n,
    freqs,
    block_sizes,
    *,
    rho=0.9,
    method: str = "scan",
    phased: bool = False,
    seed: Union[int, np.random.Generator, None] = 0,
):
    """Discrete populations with exact sizes and their own LD.

    ``n`` gives one sample size per population and ``freqs`` the
    ``(populations, m)`` counted-allele frequencies, e.g. from
    :func:`drift_frequencies`. ``block_sizes`` is one vector of block
    lengths summing to ``m``, shared by all populations, or a list with one
    such vector per population. ``rho`` is the latent AR(1) correlation, a
    scalar or one value per population. Each population is drawn by
    :func:`simulate_ar1_blocks` (``method`` as there) in order, continuing
    one random stream, and the rows are concatenated.

    Returns ``(G, labels)``: ``G`` int8 ``(sum(n), m)`` dosages, or
    ``(sum(n), 2, m)`` haplotypes with ``phased=True``, and ``labels`` the
    population index of each row (sorted, ``n[k]`` rows each).
    """
    sizes_n = _sample_sizes(n)
    k = len(sizes_n)
    f = _freq_matrix(freqs, k)
    m = f.shape[1]
    per_pop = (not np.isscalar(block_sizes) and len(block_sizes) == k
               and all(np.ndim(b) == 1 for b in block_sizes))
    blocks =([_block_vector(b, m) for b in block_sizes] if per_pop
              else [_block_vector(block_sizes, m)] * k)
    rho = _per_pop("rho", rho, k, -1.0, 1.0)
    rng = np.random.default_rng(seed)
    G = np.empty((sum(sizes_n), 2, m) if phased else (sum(sizes_n), m), dtype=np.int8)
    row = 0
    for pop, size in enumerate(sizes_n):
        G[row:row + size], _ = simulate_ar1_blocks(
            size, blocks[pop], maf=f[pop], rho=rho[pop], seed=rng,
            method=method, phased=phased)
        row += size
    return G, np.repeat(np.arange(k), sizes_n)


def _genetic_map(cm, chromosome, m):
    """Junction probabilities per gap from a map: 1 at each chromosome start."""
    pos = np.asarray(cm, dtype=float)
    if pos.ndim == 0:
        if not np.isfinite(pos) or pos < 0:
            raise ValueError("cm spacing must be finite and nonnegative")
        pos = np.arange(m) * float(pos)
    if pos.shape != (m,) or not np.isfinite(pos).all() or np.any(pos < 0):
        raise ValueError("cm must be a nonnegative spacing or a length-m map of finite positions")
    chrom = np.zeros(m, dtype=np.int64) if chromosome is None else np.asarray(chromosome)
    if chrom.shape != (m,):
        raise ValueError("chromosome must have one entry per variant")
    start = np.r_[True, chrom[1:] != chrom[:-1]]
    if len(np.unique(chrom)) != start.sum():
        raise ValueError("each chromosome must occur in one contiguous block")
    gap = np.diff(pos, prepend=pos[0])
    if np.any(gap[~start] < 0):
        raise ValueError("cm must be nondecreasing within each chromosome")
    return gap, start


def simulate_admixed(
    n: int,
    freqs,
    block_sizes,
    proportions,
    *,
    generations: float = 10,
    cm=0.01,
    chromosome=None,
    rho=0.9,
    phased: bool = False,
    seed: Union[int, np.random.Generator, None] = 0,
):
    """Admixed genotypes with local-ancestry truth.

    Each haplotype is a mosaic of tracts. Junctions between tracts form a
    Poisson process along the genetic map at rate ``generations`` per
    Morgan, the recombination of a single admixture pulse that many
    generations ago. Each tract's ancestry is drawn from the individual's
    ``proportions``, and each chromosome starts afresh. Inside a tract the
    haplotype follows its ancestry's model from :func:`simulate_populations`
    (``method='scan'``): frequencies ``freqs[k]``, latent AR(1) correlation
    ``rho[k]``, restarted at block edges. A junction starts a new founder
    haplotype, so local LD is ancestry-specific and zero across junctions,
    while shared ancestry tracts create admixture LD between nearby blocks.

    ``proportions`` is one ``(K,)`` vector for everyone or an ``(n, K)``
    array of per-individual rows; rows are normalized to sum to one.
    ``cm`` is a constant spacing in centimorgans between adjacent variants
    or a length-m map, nondecreasing within each contiguous
    ``chromosome`` label. ``block_sizes`` are shared LD block lengths
    summing to m; ``rho`` is a scalar or one value per ancestry.

    Returns ``(G, local_ancestry)``: ``G`` int8 ``(n, m)`` dosages, or
    ``(n, 2, m)`` haplotypes with ``phased=True`` (the same draws: phased
    output sums to the dosages), and ``local_ancestry`` int8 ``(n, 2, m)``,
    the ancestry index of each haplotype at each variant.
    """
    n = _positive_int("n", n)
    f = _freq_matrix(freqs)
    k, m = f.shape
    sizes = _block_vector(block_sizes, m)
    alpha = np.asarray(proportions, dtype=float)
    if alpha.ndim == 1:
        alpha = np.broadcast_to(alpha, (n, alpha.size))
    if (alpha.shape != (n, k) or not np.isfinite(alpha).all() or np.any(alpha < 0)
            or np.any(alpha.sum(axis=1) <= 0)):
        raise ValueError("proportions must be nonnegative, one per ancestry, globally or per individual")
    alpha = alpha / alpha.sum(axis=1, keepdims=True)
    generations = _finite_scalar("generations", generations)
    if generations <= 0:
        raise ValueError("generations must be positive")
    rho = _per_pop("rho", rho, k, -1.0, 1.0)
    gap, start = _genetic_map(cm, chromosome, m)
    p_junction = np.where(start, 1.0, -np.expm1(-generations * gap / 100.0))
    thr = norm_isf(f, clip=False)
    innovation = np.sqrt(1.0 - rho * rho)
    # Two haplotypes per person, adjacent rows; ancestry cut points per row.
    cuts = np.repeat(np.cumsum(alpha, axis=1)[:, :-1], 2, axis=0)
    rng = np.random.default_rng(seed)
    h = 2 * n
    hap = np.empty((h, m), dtype=np.int8)
    anc = np.empty((h, m), dtype=np.int8)
    current = np.zeros(h, dtype=np.int64)
    col = 0
    for size in sizes:
        size = int(size)
        jump = rng.random((h, size)) < p_junction[col:col + size]
        pick = (rng.random((h, size))[:, :, None] >= cuts[:, None, :]).sum(axis=2)
        eps = rng.standard_normal((h, size))
        z = eps[:, 0].copy()  # the latent chain restarts at each block edge
        for j in range(size):
            current = np.where(jump[:, j], pick[:, j], current)
            if j:
                z = np.where(jump[:, j], eps[:, j],
                             rho[current] * z + innovation[current] * eps[:, j])
            anc[:, col + j] = current
            hap[:, col + j] = z > thr[current, col + j]
        col += size
    local =anc.reshape(n, 2, m)
    G = hap.reshape(n, 2, m)
    return (G if phased else G.sum(axis=1, dtype=np.int8)), local


def simulate_split_coalescent(
    n,
    m: int,
    block_size: int = 200,
    *,
    fst: float = 0.1,
    admixed: int = 0,
    proportions=None,
    generations: float = 10,
    Ne: float = 10_000,
    recomb_rate: float = 1e-8,
    mut_rate: float = 1e-8,
    min_maf: float = 0.01,
    seed: Union[int, None] = None,
):
    """Coalescent genotypes for K populations split from one ancestor (msprime).

    The populations (``n[k]`` diploids each, K >= 2, all of size ``Ne``)
    split together from an ancestral population of size ``Ne`` at
    ``t = 2 Ne fst / (1 - fst)`` generations, so the expected Hudson F_ST
    between any two is ``fst``. Population-specific LD and frequencies
    follow from the genealogy. With ``admixed > 0`` a further population of
    that many diploids forms ``generations`` ago (before the split, i.e.
    ``generations < t``) as a pulse from the K populations in
    ``proportions``. Its local ancestry is the source population of each
    haplotype's ancestor between the pulse and the split.

    As in :func:`simulate_coalescent`, the segment grows until at least
    ``m`` SNPs have pooled minor-allele frequency above ``min_maf``; the
    first ``m'`` (``m`` rounded down to a multiple of ``block_size``) are
    kept as contiguous blocks. Dosages count the derived allele.

    Returns a dict with ``G`` int8 ``(N, m')``, ``population`` (row labels
    ``0..K-1``, and ``K`` for admixed rows, which come last), ``positions``
    (bp), ``blocks`` and ``local_ancestry`` int8 ``(admixed, 2, m')`` or
    ``None``. Requires ``pip install phensim[msprime]``.
    """
    sizes = _sample_sizes(n)
    k = len(sizes)
    if k < 2:
        raise ValueError("n must give at least two populations")
    m = _positive_int("m", m)
    block_size = _positive_int("block_size", block_size)
    if block_size > m:
        raise ValueError("block_size must not exceed m")
    Ne, recomb_rate, mut_rate, min_maf = _coalescent_params(Ne, recomb_rate, mut_rate, min_maf)
    if mut_rate == 0:
        raise ValueError("mut_rate must be positive to reach a SNP-count target")
    fst = _finite_scalar("fst", fst)
    if not 0 < fst < 1:
        raise ValueError("fst must be in (0, 1)")
    split = 2 * Ne * fst / (1 - fst)
    if isinstance(admixed, (bool, np.bool_)) or not isinstance(admixed, (int, np.integer)) or admixed < 0:
        raise ValueError("admixed must be a nonnegative integer")
    admixed = int(admixed)
    if admixed:
        alpha = np.asarray(proportions, dtype=float) if proportions is not None else None
        if (alpha is None or alpha.shape != (k,) or not np.isfinite(alpha).all()
                or np.any(alpha < 0) or alpha.sum() <= 0):
            raise ValueError("admixed samples need nonnegative proportions, one per population")
        alpha = alpha / alpha.sum()
        generations = _finite_scalar("generations", generations)
        if not 0 < generations < split:
            raise ValueError(f"generations must be in (0, split time {split:.4g})")
    try:
        import msprime
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "simulate_split_coalescent needs msprime (pip install phensim[msprime])") from e

    names = [f"pop{i}" for i in range(k)]
    demography = msprime.Demography()
    for name in names + (["admixed"] if admixed else []) + ["ancestral"]:
        demography.add_population(name=name, initial_size=Ne)
    if admixed:
        demography.add_admixture(time=generations, derived="admixed", ancestral=names,
                                 proportions=alpha.tolist())
        demography.add_census(time=0.5 * (generations + split))
    demography.add_population_split(time=split, derived=names, ancestral="ancestral")
    demography.sort_events()
    samples = dict(zip(names, sizes))
    if admixed:
        samples["admixed"] = admixed
    total = sum(sizes) + admixed
    m2 = (m // block_size) * block_size
    rng = np.random.default_rng(seed)
    seq_len = max(1e6, m / 1200 * 1e6)
    for _ in range(7):
        rep_seed = int(rng.integers(1, 2**31 - 1))
        ts = msprime.sim_ancestry(
            samples=samples, demography=demography, sequence_length=int(seq_len),
            recombination_rate=recomb_rate, discrete_genome=False,
            random_seed=rep_seed, record_provenance=False)
        mts = msprime.sim_mutations(
            ts, rate=mut_rate, random_seed=rep_seed, discrete_genome=False,
            model=msprime.BinaryMutationModel())
        dos = np.empty((m2, total), dtype=np.int8)
        positions = np.empty(m2)
        kept = 0
        lo, hi = 2 * total * min_maf, 2 * total * (1 - min_maf)
        for v in mts.variants(copy=False):
            g = v.genotypes  # haplotypes; an individual's two are adjacent
            if lo < g.sum() < hi:
                np.add(g[0::2], g[1::2], out=dos[kept], casting="unsafe")
                positions[kept] = v.site.position
                kept += 1
                if kept == m2:
                    break
        if kept == m2:
            break
        seq_len *= 1.8
    else:
        raise RuntimeError("split coalescent produced too few common SNPs; increase Ne or mut_rate")
    tables = mts.tables
    node_pop = tables.nodes.population
    sample_nodes = mts.samples()
    population = node_pop[sample_nodes[0::2]].astype(np.int64)
    local = None
    if admixed:
        target = sample_nodes[2 * sum(sizes):]  # admixed samples come last
        census = np.flatnonzero(tables.nodes.flags & msprime.NODE_IS_CEN_EVENT)
        links = tables.link_ancestors(target, census)
        row = np.full(mts.num_nodes, -1)
        row[target] = np.arange(target.size)
        local = np.full((target.size, m2), -1, dtype=np.int8)
        left = np.searchsorted(positions, links.left)
        right = np.searchsorted(positions, links.right)
        for child, parent, a, b in zip(links.child, links.parent, left, right):
            local[row[child], a:b] = node_pop[parent]
        if np.any(local < 0):  # pragma: no cover - census covers every lineage
            raise RuntimeError("local ancestry left a site unassigned")
        local = local.reshape(admixed, 2, m2)
    blocks = [np.arange(i, i + block_size) for i in range(0, m2, block_size)]
    return {
        "G": np.ascontiguousarray(dos.T),
        "population": population,
        "positions": positions,
        "blocks": blocks,
        "local_ancestry": local,
    }
