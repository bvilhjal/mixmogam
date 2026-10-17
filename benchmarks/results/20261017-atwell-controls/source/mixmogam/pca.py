"""Ancestry principal components of a thinned marker panel."""

from __future__ import annotations

from typing import Optional, Union

import numpy as np

from mixmogam.genotypes import Genotypes

__all__ = ["principal_components"]


def principal_components(
    G: Union[Genotypes, np.ndarray],
    k: int = 10,
    *,
    min_maf: float = 0.05,
    n_markers: Optional[int] = 2000,
    random_state=0,
    block: int = 2048,
    max_memory_gib: float = 4.0,
    return_info: bool = False,
):
    """Leading ancestry principal components for use as covariates.

    ``G`` is a :class:`~mixmogam.genotypes.Genotypes` container or an
    (n_samples, n_variants) hard-call matrix. The panel is thinned to at
    most ``n_markers`` variants with sample MAF at least ``min_maf``, drawn
    at random without replacement (``random_state``), since PCs are usually
    computed from a thinned panel rather than every variant; ``n_markers=None``
    keeps every variant that passes the MAF filter (the full standardized
    panel then has to fit in memory as float64). Each retained variant is
    standardized to mean zero and standard deviation one (no-calls are
    mean-imputed first), and the result is the leading left singular
    vectors of that standardized thinned panel: the columns are
    ``numpy.linalg.svd``'s U factors (unit norm), so the output matches the
    SVD of the panel selected in ``info["variants"]``. Column signs are
    fixed (the largest-magnitude entry of each is positive, the first such
    entry on ties) so the values are reproducible across platforms and
    independent of allele coding. As covariates their scale is irrelevant
    (regression refits the coefficients); multiply by ``sqrt(n)`` for
    unit-variance scores.

    With ``return_info=True`` the result is ``(pcs, info)``, where ``info``
    holds ``variants`` (the thinned variant indices), ``singular_values``
    and ``explained_variance_ratio``.

    Before scanning allele frequencies or allocating the panel, a conservative
    working-memory estimate must fit ``max_memory_gib`` (default 4 GiB).
    It includes the requested panel, decoding, SVD factors and workspace,
    but excludes resident input calls and library/allocator overhead. This
    is an allocation guard, not a measured peak-RSS guarantee. Lower
    ``n_markers`` or supply external PCs if the panel exceeds the budget;
    increasing the budget leaves the float64 SVD unchanged.

    ``random_state=None`` draws fresh thinning randomness; every other
    value seeds ``numpy.random.default_rng``. Before this was built in,
    each benchmark computed its own PCs (randomized subspace iteration in
    `hratt_followups.py`, a seeded Lanczos in `ldak_kvik_comparison.py`),
    and callers had to assemble the thinning and standardization by hand.
    """
    if not isinstance(G, Genotypes):
        G = Genotypes(G)
    if isinstance(k, (bool, np.bool_)) or not isinstance(k, (int, np.integer)) or k < 1:
        raise ValueError("k must be a positive integer")
    if not (isinstance(min_maf, (int, float, np.floating)) and not isinstance(min_maf, (bool, np.bool_))
            and np.isfinite(min_maf) and 0.0 <= min_maf <= 0.5):
        raise ValueError("min_maf must lie in [0, 0.5]")
    if n_markers is not None and (isinstance(n_markers, (bool, np.bool_))
                                  or not isinstance(n_markers, (int, np.integer)) or n_markers < 1):
        raise ValueError("n_markers must be a positive integer or None")
    if not isinstance(block, (int, np.integer)) or isinstance(block, (bool, np.bool_)) or block < 1:
        raise ValueError("block must be a positive integer")
    if (isinstance(max_memory_gib, (bool, np.bool_))
            or not isinstance(max_memory_gib, (int, float, np.integer, np.floating))
            or not np.isfinite(max_memory_gib) or max_memory_gib <= 0):
        raise ValueError("max_memory_gib must be positive and finite")

    n = G.n_samples
    m = G.n_variants if n_markers is None else min(n_markers, G.n_variants)
    rank = min(n, m)
    # Panel + its LAPACK copy + U, Vt, divide-and-conquer workspace,
    # and the last decoded block, which remains live at the SVD call.
    svd_bytes = 8 * (2 * n * m + n * rank + rank * m + 4 * rank**2 + n * min(block, m))
    # allele_freqs uses blocks of 8192 int8 calls, a mask and a masked copy.
    frequency_bytes = 3 * n * min(8192, G.n_variants) + 24 * G.n_variants
    estimated_bytes = max(svd_bytes, frequency_bytes)
    if estimated_bytes > max_memory_gib * 1024**3:
        raise ValueError(f"PCA working-memory estimate {estimated_bytes / 1024**3:.2f} GiB "
                         f"exceeds max_memory_gib={max_memory_gib:g} for {n} samples and "
                         f"up to {m} markers; lower n_markers, use external PCs, or "
                         "explicitly increase max_memory_gib")

    maf = G.allele_freqs(minor=True)
    common = np.flatnonzero(maf >= min_maf)  # NaN MAF (all-missing) never passes
    if common.size == 0:
        raise ValueError(f"no variants with sample MAF >= {min_maf:g}; lower min_maf")
    if n_markers is None:
        take = common
    else:
        rng = np.random.default_rng(random_state)
        # The draw is the one hratt_followups.top_pcs made, so its thinned
        # panels (and therefore its S4 covariates) are unchanged.
        take = np.sort(rng.choice(common, min(n_markers, common.size), replace=False))
    if k > min(n, take.size):
        raise ValueError(f"k={k} exceeds the {min(n, take.size)} components available "
                         f"from {take.size} thinned variants and {n} samples")

    Z = np.empty((n, take.size), dtype=np.float64)
    at = 0
    for blockZ in G.iter_snp_blocks(block, dtype=np.float64, impute="mean", variant_indices=take):
        Z[:, at:at + blockZ.shape[0]] = blockZ.T
        at += blockZ.shape[0]
    mean = Z.mean(axis=0)
    sd = Z.std(axis=0)
    if np.any(sd == 0):
        raise ValueError("the thinned panel contains monomorphic variants; filter them or raise min_maf")
    Z -= mean
    Z /= sd
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    pcs = U[:, :k]
    pivot = np.argmax(np.abs(pcs), axis=0)
    signs = np.sign(pcs[pivot, np.arange(pcs.shape[1])])
    signs[signs == 0] = 1.0
    pcs = pcs * signs
    if not return_info:
        return pcs
    total = float(s @ s)
    info = {"variants": take,
            "singular_values": s[:k],
            "explained_variance_ratio": (s[:k] ** 2) / total if total > 0 else np.zeros(k)}
    return pcs, info
