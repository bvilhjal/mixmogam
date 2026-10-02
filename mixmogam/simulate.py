"""Genotype and trait simulators for calibration and power testing."""

from __future__ import annotations

from typing import Union

import numpy as np

__all__ = ["simulate_genotypes", "simulate_traits", "simulate_kinship"]


def simulate_genotypes(
    n: int,
    m: int,
    maf: float = 0.3,
    n_pop: int = 1,
    pop_fst: float = 0.0,
    seed: Union[int, np.random.Generator, None] = 0,
) -> np.ndarray:
    """(m, n) int8 genotype matrix in {0, 1, 2}.

    With ``n_pop > 1`` allele frequencies are drawn per population
    (``pop_fst`` scaling the between-population frequency spread, giving
    structure-rich kinships) and individuals inherit their population's
    frequencies.
    """
    rng = np.random.default_rng(seed)
    if n_pop <= 1 or pop_fst == 0.0:
        freqs = np.full(m, maf) + rng.normal(0, 0.05, m)
    else:
        base = np.clip(maf + rng.normal(0, 0.05, m), 0.05, 0.95)
        spread = pop_fst * base * (1 - base)
        freqs = np.clip(
            base[:, None] + rng.normal(0, 1, (m, n_pop)) * np.sqrt(spread)[:, None],
            0.01,
            0.99,
        )
    pops = rng.integers(0, n_pop, n)
    if n_pop <= 1 or pop_fst == 0.0:
        p = np.broadcast_to(freqs[:, None], (m, n))
    else:
        p = freqs[:, pops]
    G = rng.binomial(2, p).astype(np.int8)
    return G


def simulate_traits(
    G: np.ndarray,
    h2: float = 0.5,
    n_causal: int = 10,
    seed: Union[int, np.random.Generator, None] = 1,
    effect_dist: str = "normal",
) -> dict:
    """Simulate a standardized trait with polygenic background and QTLs.

    Returns ``{"y", "u", "causal", "effects"}``: ``u`` is the full
    infinitesimal genetic value (u ~ N(0, h2 K) in distribution, drawn
    through the genotype projection so it matches the GRM the model will
    fit), on top of which ``n_causal`` QTL effects of total variance
    ``h2`` replace a like share of the background. ``effect_dist='equal'``
    gives every causal locus the same (absolute) effect so per-locus power
    is deterministic; ``'normal'`` draws N(0,1) effects (v1 style).
    """
    rng = np.random.default_rng(seed)
    Gd = np.asarray(G, dtype=np.float64)
    m, n = Gd.shape
    Z = Gd - Gd.mean(axis=1, keepdims=True)
    sd = Z.std(axis=1, keepdims=True)
    Z = Z / np.where(sd > 0, sd, 1.0)

    u = Z.T @ rng.standard_normal(m)
    u *= np.sqrt(h2) / u.std()

    causal = rng.choice(m, size=min(n_causal, m), replace=False)
    if effect_dist == "equal":
        effects = np.ones(causal.size) * rng.choice([-1.0, 1.0])
    else:
        effects = rng.normal(0, 1, causal.size)
    q = Z[causal].T @ effects
    q *= np.sqrt(h2) / q.std()
    # combine: half infinitesimal, half QTL (both var h2/2)
    u = (u + q) / np.sqrt(2)
    e = rng.standard_normal(n) * np.sqrt(1 - h2)
    y = u + e
    y = (y - y.mean()) / y.std()
    return {"y": y, "u": u, "causal": causal, "effects": effects}


def simulate_kinship(G: np.ndarray) -> np.ndarray:
    """Additive GRM from a SNP-major genotype matrix, mean diagonal 1."""
    Gd = np.asarray(G, dtype=np.float64)
    Z = Gd - Gd.mean(axis=1, keepdims=True)
    sd = Z.std(axis=1, keepdims=True)
    Z = Z / np.where(sd > 0, sd, 1.0)
    from mixmogam.kinship import scale_k

    K = (Z.T @ Z) / Gd.shape[0]
    return scale_k(K)
