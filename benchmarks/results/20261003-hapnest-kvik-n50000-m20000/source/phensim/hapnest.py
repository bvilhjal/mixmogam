"""HAPNEST's age-dependent haplotype-copying model, in bounded memory.

This independent implementation follows Wharrie et al. (2023),
doi:10.1093/bioinformatics/btad535. It is not the upstream Julia program
or its ABC parameter-fitting pipeline. See docs/hapnest.md for conventions.
"""
from __future__ import annotations

import numpy as np

from phensim._numba import HAVE_NUMBA, _jit
from phensim.genotypes import _positive_int

__all__ = ["simulate_hapnest", "iter_hapnest"]


def _fill_numpy(out, H, cm, ages, boundaries, donors, offsets, groups, ne, rho, rng):
    """Simple segment-wise reference implementation; no compiled code."""
    for i, group in enumerate(groups):
        count = offsets[group + 1] - offsets[group]
        out[i] = 0
        for phase in range(2):
            for left, right in zip(boundaries[:-1], boundaries[1:]):
                start = left
                while start < right:
                    age = rng.gamma(2.0, ne[group] / count)
                    length = rng.exponential(1 / (2 * age * rho[group])) if rho[group] else np.inf
                    donor = donors[offsets[group] + rng.integers(0, count)]
                    # Upstream includes the first marker beyond the breakpoint.
                    stop = min(right, start + np.searchsorted(
                        cm[start:right], cm[start] + length, side="right") + 1)
                    # Binary alleles fit int8; convert only this segment so
                    # uint64 references do not promote accumulation to float64.
                    alleles = H[donor, phase, start:stop].astype(np.int8, copy=False)
                    out[i, start:stop] += alleles * (age < ages[start:stop])
                    start = stop


@_jit
def _fill_numba(out, H, cm, ages, boundaries, donors, offsets, groups, ne, rho, rng):
    """Fuse random segment generation, age filtering and diploid accumulation."""
    for i in range(len(groups)):
        group = groups[i]
        count = offsets[group + 1] - offsets[group]
        for j in range(out.shape[1]):
            out[i, j] = 0
        for phase in range(2):
            for chromosome in range(len(boundaries) - 1):
                start, right = boundaries[chromosome], boundaries[chromosome + 1]
                while start < right:
                    age = rng.gamma(2.0, ne[group] / count)
                    length = rng.exponential(1 / (2 * age * rho[group])) if rho[group] else np.inf
                    donor = donors[offsets[group] + rng.integers(0, count)]
                    # Bounded binary search: no scan through the genetic map.
                    lo, hi = start, right
                    target = cm[start] + length
                    while lo < hi:
                        mid = (lo + hi) // 2
                        if cm[mid] <= target:
                            lo = mid + 1
                        else:
                            hi = mid
                    stop = min(right, lo + 1)
                    for j in range(start, stop):
                        if age < ages[j]:
                            out[i, j] += H[donor, phase, j]
                    start = stop


def _labels(value, size, name):
    if value is None:
        return np.zeros(size, dtype=np.int64)
    a = np.asarray(value)
    if a.shape != (size,) or a.dtype.kind not in "iu" or np.any(a > np.iinfo(np.int64).max) or np.any(a < 0):
        raise ValueError(f"{name} must be a length-{size} vector of nonnegative integers")
    return a.astype(np.int64, copy=False)


def _prepare(reference, n, genetic_map, mutation_age, reference_populations,
             sample_populations, chromosome, ne, rho, backend):
    n = _positive_int("n", n)
    H = np.asarray(reference)
    if H.ndim != 3 or H.shape[1] != 2 or not H.shape[0] or not H.shape[2] or H.dtype.kind not in "biu":
        raise ValueError("reference must have shape (reference individuals, 2, variants) with binary integer haplotypes")
    # Validate even disk-backed references without allocating a full boolean copy.
    for start in range(0, H.shape[0], max(1, (1 << 20) // (2 * H.shape[2]))):
        tile = H[start:start + max(1, (1 << 20) // (2 * H.shape[2]))]
        if tile.min() < 0 or tile.max() > 1:
            raise ValueError("reference haplotypes must contain only 0 and 1")
    m = H.shape[2]
    cm, ages = np.asarray(genetic_map, dtype=float), np.asarray(mutation_age, dtype=float)
    if cm.shape != (m,) or not np.isfinite(cm).all() or np.any(cm < 0):
        raise ValueError("genetic_map must contain finite nonnegative centimorgan positions")
    if ages.shape != (m,) or np.isnan(ages).any() or np.any(ages < 0):
        raise ValueError("mutation_age must contain nonnegative generations (infinity disables age filtering)")
    chrom = np.zeros(m, dtype=int) if chromosome is None else np.asarray(chromosome)
    if chrom.shape != (m,):
        raise ValueError("chromosome must have one entry per variant")
    boundaries = np.r_[0, np.flatnonzero(chrom[1:] != chrom[:-1]) + 1, m].astype(np.int64)
    if len(np.unique(chrom)) != len(boundaries) - 1:
        raise ValueError("each chromosome must occur in one contiguous block")
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        if np.any(np.diff(cm[left:right]) < 0):
            raise ValueError("genetic_map must be nondecreasing within each chromosome")
    refs = _labels(reference_populations, H.shape[0], "reference_populations")
    pops = np.unique(refs)
    if not np.array_equal(pops, np.arange(len(pops))):
        raise ValueError("reference population codes must be consecutive from zero")
    groups = _labels(sample_populations, n, "sample_populations")
    if np.any(groups >= len(pops)):
        raise ValueError("every sample population must have reference donors")
    params = []
    for name, value in (("ne", ne), ("rho", rho)):
        arr = np.asarray(value, dtype=float)
        if arr.ndim == 0:
            arr = np.full(len(pops), float(arr))
        if arr.shape != (len(pops),) or not np.isfinite(arr).all() or np.any(arr < 0) or (name == "ne" and np.any(arr == 0)):
            raise ValueError(f"{name} must be finite {'positive' if name == 'ne' else 'nonnegative'}, scalar or one value per reference population")
        params.append(arr)
    donors = np.argsort(refs, kind="stable")
    offsets = np.r_[0, np.cumsum(np.bincount(refs))]
    if backend not in ("auto", "numpy", "numba"):
        raise ValueError("backend must be 'auto', 'numpy' or 'numba'")
    if backend == "numba" and not H.dtype.isnative:
        raise ValueError("backend='numba' requires reference haplotypes in native byte order; use backend='numpy' or convert the reference")
    if backend == "numba" and not HAVE_NUMBA:
        raise ImportError("backend='numba' requires phensim[fast]")
    fill = _fill_numba if HAVE_NUMBA and backend != "numpy" and H.dtype.isnative else _fill_numpy
    return fill, (H, cm, ages, boundaries, donors, offsets, groups, *params)


def iter_hapnest(reference, n, genetic_map, mutation_age, *, reference_populations=None,
                 sample_populations=None, chromosome=None, ne=10000, rho=0.7,
                 seed=0, batch_size=256, backend="auto"):
    """Yield independent int8 dosage batches under the HAPNEST model.

    ``reference`` is a binary array (reference diploids, 2, variants),
    optionally memory-mapped. Genetic-map positions are in cM; mutation
    ages and segment ages are in generations. Population codes are
    consecutive integers from zero, with one code per reference/output
    diploid. Each output samples donors only from its specified group.
    ``ne`` and ``rho`` are scalar or one value per reference population.
    This supports discrete populations, not a demographic admixture model.

    Each segment draws T ~ Gamma(shape=2, scale=ne/N_reference_diploids),
    then length ~ Exp(rate=2*T*rho) in cM. Only alleles older than T are
    copied. Like upstream, haplotype 1 copies reference phase 1 and
    haplotype 2 copies phase 2; segments include the first marker beyond
    the breakpoint. Chromosome boundaries always reset the process.

    The default uses optional Numba for native-byte-order references and
    otherwise falls back to NumPy without copying the full reference.
    ``backend='numpy'`` is the simple reference path; explicit Numba
    requires native byte order. Both use the same NumPy Generator stream, consume
    draws in the same order, and are independent of batch size. Working
    storage is O(batch_size * markers + samples + reference individuals);
    no segment table or full synthetic haplotype matrix is constructed.
    Returned batches have independent storage and may be retained.
    """
    batch_size = _positive_int("batch_size", batch_size)
    fill, params = _prepare(reference, n, genetic_map, mutation_age,
                            reference_populations, sample_populations, chromosome, ne, rho, backend)
    H, cm, ages, boundaries, donors, offsets, groups, ne, rho = params
    rng = np.random.default_rng(seed)
    for start in range(0, n, batch_size):
        out = np.empty((min(batch_size, n - start), H.shape[2]), dtype=np.int8)
        fill(out, H, cm, ages, boundaries, donors, offsets, groups[start:start + batch_size], ne, rho, rng)
        yield out


def simulate_hapnest(reference, n, genetic_map, mutation_age, *, out=None, **kwargs):
    """Materialize :func:`iter_hapnest` into an int8 (n, m) array.

    Supply a writable int8 array or memmap as ``out`` to avoid allocating
    the complete result in RAM. Validation precedes output modification.
    Use the iterator when downstream work can consume sample batches.
    """
    batches = iter_hapnest(reference, n, genetic_map, mutation_age, **kwargs)
    first = next(batches)  # validate before allocating/modifying the full result
    shape = (n, first.shape[1])
    if out is None:
        out = np.empty(shape, dtype=np.int8)
    elif not isinstance(out, np.ndarray) or out.shape != shape or out.dtype != np.int8 or not out.flags.writeable:
        raise ValueError("out must be a writable int8 (n, variants) array")
    if np.may_share_memory(out, np.asarray(reference)):
        raise ValueError("out must not overlap the reference haplotypes")
    out[:len(first)] = first
    start = len(first)
    for batch in batches:
        out[start:start + len(batch)] = batch
        start += len(batch)
    return out
