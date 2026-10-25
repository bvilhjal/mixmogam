"""Phenotype simulators on top of genotype data.

All quantitative-trait simulators draw the infinitesimal component as
``u ~ N(0, sigma2 K)`` on the empirical GRM, so the data-generating
covariance matches what a mixed model will fit. ``simulate_trait`` (and
the binary/GxE wrappers delegating to it) and
``simulate_correlated_traits`` draw matrix-free with no supplied
kinship -- innovations in marker space through an exact factor of the
scaled GRM -- while a caller-supplied ``K`` keeps the classic
eigendecomposition draw. ``simulate_confounded_trait`` only materializes
the GRM when no explicit environmental axis is supplied, sharing one
eigendecomposition between its default axis and the background.
"""

from __future__ import annotations

import warnings
from statistics import NormalDist
from typing import Optional, Union

import numpy as np

from phensim._common import norm_isf
from phensim.kinship import _called_standardized

__all__ = [
    "simulate_trait",
    "simulate_binary_trait",
    "simulate_confounded_trait",
    "simulate_gxe_trait",
    "simulate_correlated_traits",
    "ascertain_case_control",
    "n_eff_case_control",
    "h2_liability",
]


def _grm(G: np.ndarray) -> np.ndarray:
    """Additive GRM, standardized columns, mean diagonal 1."""
    from phensim.kinship import grm

    return grm(G)


def _standardized(x: np.ndarray) -> np.ndarray:
    sd = x.std()
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("cannot standardize a constant or non-finite component")
    return (x - x.mean()) / sd


def _trait_genotypes(G):
    """Validated complete dosages. Integer and float arrays keep their
    dtype: consumers convert what they read, and a float64 copy of an
    int8 matrix costs eight times its size."""
    G = np.asarray(G)
    if G.dtype.kind not in "iuf":
        G = np.asarray(G, dtype=np.float64)
    if (G.ndim != 2 or G.shape[0] < 2 or G.shape[1] == 0
            or (G.dtype.kind == "f" and not np.isfinite(G).all())):
        raise ValueError("G must be a finite matrix with at least two samples and one variant")
    if G.min() < 0 or G.max() > 2:
        raise ValueError("trait genotypes must be complete dosages in [0, 2]")
    return G


def _causal_indices(G, n_causal, causal, rng):
    """Choose observed-variable columns; never silently drop explicit indices."""
    m = G.shape[1]
    if causal is None:
        if (isinstance(n_causal, (bool, np.bool_)) or not isinstance(n_causal, (int, np.integer))
                or n_causal < 0):
            raise ValueError("n_causal must be a nonnegative integer")
        if n_causal == 0:
            return np.empty(0, dtype=np.intp)
        # Exact constancy, with O(m) workspace even for mapped/blocked G.
        eligible = np.flatnonzero(G.min(axis=0) != G.max(axis=0))
        return rng.choice(eligible, size=min(n_causal, eligible.size), replace=False)
    causal = np.asarray(causal)
    if (causal.ndim != 1 or not np.issubdtype(causal.dtype, np.integer)
            or np.any(causal < 0) or np.any(causal >= m) or np.unique(causal).size != causal.size):
        raise ValueError("causal must contain distinct in-range integer indices")
    for j in causal:
        if G[:, j].min() == G[:, j].max():
            raise ValueError(f"causal variant {j} is constant in the supplied genotypes")
    return causal


def _scaled_score(Z, effects, variance):
    """Scale the coefficients and their score together on the liability scale."""
    score = Z @ effects
    sd = score.std()
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("positive QTL/interaction variance requires a nonconstant causal score")
    factor = np.sqrt(variance) / sd
    return score * factor, effects * factor


def _trait_kinship(K, n: int) -> np.ndarray:
    """A caller-supplied kinship: finite, symmetric, ``(n, n)``.

    The input dtype is preserved -- an eigendecomposition of a float32
    kinship must not silently promote to float64, since callers may rely
    on the older seeded draw's arithmetic.
    """
    K = np.asarray(K)
    if np.iscomplexobj(K) or not np.issubdtype(K.dtype, np.number):
        raise ValueError("K must be a real numeric kinship matrix")
    if (K.shape != (n, n) or not np.isfinite(K).all()
            or not np.allclose(K, K.T, rtol=1e-7, atol=1e-10)):
        raise ValueError("K must be a finite symmetric (n, n) kinship matrix")
    return K


def _background_factor(G: np.ndarray):
    """Matrix-free factor of the EMMAX-scaled GRM: ``(Z, a, b)`` with
    ``F F' = K`` for ``F = [a Z, b 1]``.

    ``Z`` is the Yang-2010 called-only standardized genotype matrix.
    Every column of ``Z`` sums to zero, so the raw Gram ``ZZ'/m`` has a
    zero grand mean and EMMAX mean off-diagonal ``-S/(m n (n-1))`` with
    ``S = sum(Z**2)``; the scaled GRM is therefore exactly
    ``((n-1)/S) ZZ' + J/n``.
    """
    Z, _cnt = _called_standardized(G)  # G passed _trait_genotypes
    n = Z.shape[0]
    if n < 2:
        raise ValueError("a genetic background needs at least two samples")
    S = np.einsum("ij,ij->", Z, Z)
    if not np.isfinite(S) or S <= 0:
        raise ValueError(
            "no polymorphic genotype variance for a genetic background")
    marker_scale = np.sqrt((n - 1) / S)
    common_scale = 1.0 / np.sqrt(n)
    return Z, marker_scale, common_scale


def _draw_background(factor, variance: float, rng) -> np.ndarray:
    """``u ~ N(0, variance * K)`` from ``m + 1`` innovations, no ``n x n``
    matrix. The shared scalar innovation carries the ``J/n`` part of the
    scaled GRM, so the draw is not centred afterwards."""
    Z, marker_scale, common_scale = factor
    w = rng.standard_normal(Z.shape[1])
    common = rng.standard_normal()
    return np.sqrt(variance) * (marker_scale * (Z @ w) + common_scale * common)


def _draw_background_blocked(G, variance, rng, block_size):
    """Same factor and innovations as the dense draw, accumulated in tiles."""
    n, m = G.shape
    w, common = rng.standard_normal(m), rng.standard_normal()
    score, ss = np.zeros(n), 0.0
    for start in range(0, m, block_size):
        Z, _ = _called_standardized(G[:, start:start + block_size])
        score += Z @ w[start:start + block_size]
        ss += np.einsum("ij,ij->", Z, Z)
    if not np.isfinite(ss) or ss <= 0:
        raise ValueError("no polymorphic genotype variance for a genetic background")
    return np.sqrt(variance) * (np.sqrt((n - 1) / ss) * score + common / np.sqrt(n))


def simulate_trait(
    G: np.ndarray,
    h2: float = 0.5,
    n_causal: int = 20,
    architecture: str = "mixed",
    effect_dist: str = "normal",
    causal: Optional[np.ndarray] = None,
    K: Optional[np.ndarray] = None,
    seed: Union[int, np.random.Generator, None] = 1,
    *,
    genotype_block_size: Optional[int] = None,
) -> dict:
    """Quantitative trait with a model-consistent genetic component.

    ``architecture`` splits the target heritability ``h2``:

    - ``'mixed'``: half infinitesimal (u ~ N(0, (h2/2) K)), half from
      ``n_causal`` QTLs;
    - ``'infinitesimal'``: all h2 through the kinship;
    - ``'qtl'``: all h2 at the causal variants, no background.

    ``effect_dist`` is ``'normal'`` (random effect sizes) or ``'equal'``
    (same absolute effect per causal variant -- deterministic per-locus
    power). Automatic causal selection samples only columns that vary
    across the supplied samples, capped at the number available. Explicit
    ``causal`` indices must all vary; constant columns raise ``ValueError``.
    Returns standardized ``y`` and raw ``liability = u + q + e``.
    ``effects`` act on the centred, unit-SD causal genotype columns and
    reconstruct ``q`` on the liability scale. Divide them by
    ``liability.std()`` for effects on the standardized-y scale. With no
    QTL variance (``'infinitesimal'`` or ``h2=0``) ``causal`` still lists
    the drawn variants, all with zero effect -- select on ``effects != 0``.
    The background draw is exact either way: with ``K=None`` it is
    matrix-free -- ``m + 1`` innovations through an exact factor of the
    scaled GRM, so no ``n x n`` matrix or eigendecomposition is formed --
    while a supplied ``K`` must be finite and symmetric and is drawn
    through its eigendecomposition without rescaling. For positive-semidefinite
    K, a mean diagonal of one makes the average marginal background variance
    equal to its h2 allocation; scaling K scales that variance. Component variances target h2;
    finite-sample variances and covariances need not give an exactly
    realized heritability. G must be complete. With ``K=None``, setting
    ``genotype_block_size`` bounds the float64 working matrix to that many
    variants. This keeps the same innovations and covariance, with only
    summation roundoff differences from the default dense factor.
    """
    return _simulate_trait(
        G, h2=h2, n_causal=n_causal, architecture=architecture,
        effect_dist=effect_dist, causal=causal, K=K, seed=seed,
        genotype_block_size=genotype_block_size)


def _simulate_trait(
    G: np.ndarray,
    h2: float = 0.5,
    n_causal: int = 20,
    architecture: str = "mixed",
    effect_dist: str = "normal",
    causal: Optional[np.ndarray] = None,
    K: Optional[np.ndarray] = None,
    seed: Union[int, np.random.Generator, None] = 1,
    *,
    eigendecomposition=None,
    genotype_block_size=None,
) -> dict:
    """Body of :func:`simulate_trait`.

    ``eigendecomposition`` is a private escape hatch for wrappers that
    already eigendecomposed the same kinship; it is not part of the
    public API.
    """
    if architecture not in ("mixed", "infinitesimal", "qtl"):
        raise ValueError(f"unknown architecture {architecture!r}")
    if not 0.0 <= h2 <= 1.0:
        raise ValueError("h2 must be in [0, 1]")
    if effect_dist not in ("normal", "equal"):
        raise ValueError("effect_dist must be 'normal' or 'equal'")
    if genotype_block_size is not None:
        from phensim.genotypes import _positive_int
        genotype_block_size = _positive_int("genotype_block_size", genotype_block_size)
    rng = np.random.default_rng(seed)
    Gd = _trait_genotypes(G)
    n = Gd.shape[0]

    h2_bg = {"mixed": h2 / 2, "infinitesimal": h2, "qtl": 0.0}[architecture]
    h2_qtl = h2 - h2_bg
    if h2_qtl > 0:
        if causal is None and n_causal == 0:
            raise ValueError(
                "a positive-QTL architecture requires at least one causal variant")
        if causal is not None and np.asarray(causal).size == 0:
            raise ValueError(
                "a positive-QTL architecture requires at least one causal variant")

    u = np.zeros(n)
    if h2_bg > 0:
        if K is None:
            if genotype_block_size is None:
                u = _draw_background(_background_factor(Gd), h2_bg, rng)
            else:
                u = _draw_background_blocked(Gd, h2_bg, rng, genotype_block_size)
        else:
            K = _trait_kinship(K, n)
            if eigendecomposition is None:
                eigendecomposition = np.linalg.eigh(K)
            lam, U = eigendecomposition
            lam = np.maximum(lam, 0.0)
            u = U @ (np.sqrt(lam * h2_bg) * rng.standard_normal(n))

    causal = _causal_indices(Gd, n_causal, causal, rng)
    if h2_qtl > 0 and causal.size == 0:
        raise ValueError(
            "positive QTL variance requires at least one polymorphic causal variant")
    if h2_qtl > 0:
        if effect_dist == "equal":
            effects = np.sign(rng.standard_normal(causal.size))
        else:
            effects = rng.normal(0, 1, causal.size)
        q, effects = _scaled_score(_standardize_cols(Gd, causal), effects, h2_qtl)
    else:
        effects = np.zeros(causal.size)
        q = np.zeros(n)

    e = rng.standard_normal(n) * np.sqrt(max(1.0 - h2, 0.0))
    liability = u + q + e
    return {
        "y": _standardized(liability),
        "liability": liability,
        "u": u,
        "q": q,
        "e": e,
        "causal": causal,
        "effects": effects,
    }


def simulate_binary_trait(
    G: np.ndarray,
    prevalence: float = 0.05,
    **trait_kwargs,
) -> dict:
    """Liability-threshold case/control trait.

    Simulates a quantitative liability via :func:`simulate_trait`, then
    thresholds at the standard-normal ``1-prevalence`` quantile. This
    targets prevalence for a unit-normal liability; sparse QTL or
    structured liabilities can give a different case fraction. Returns
    the trait dict with ``y`` binary and ``case_control`` added.
    """
    if not 0 < prevalence < 1:
        raise ValueError("prevalence must be in (0, 1)")
    tr = simulate_trait(G, **trait_kwargs)
    thresh = norm_isf(prevalence)
    cases = tr["liability"] > thresh
    tr["case_control"] = cases.astype(np.int8)
    tr["y"] = cases.astype(np.float64)
    return tr


def simulate_confounded_trait(
    G: np.ndarray,
    confounding_strength: float = 0.6,
    h2: float = 0.5,
    n_causal: int = 10,
    K: Optional[np.ndarray] = None,
    seed: Union[int, np.random.Generator, None] = 2,
    *,
    environment: Optional[np.ndarray] = None,
    architecture: str = "mixed",
    effect_dist: str = "normal",
    causal: Optional[np.ndarray] = None,
    genotype_block_size: Optional[int] = None,
) -> dict:
    """Trait driven partly by population structure.

    With ``s=confounding_strength``, the raw liability is
    ``sqrt(s) * standardized(environment) + sqrt(1-s) * base_liability``.
    These are component targets, not realized variance fractions: the
    environment and genetic score can covary. ``h2`` describes the base
    trait, whose architecture, effects and causal indices follow
    :func:`simulate_trait`.

    An explicit finite, nonconstant ``environment`` has one value per
    person (for example a population-specific exposure). It uses the
    matrix-free background draw unless K is supplied. If omitted, the
    leading kinship eigenvector supplies the axis, with its largest-magnitude
    loading made positive, sharing the eigendecomposition with the background.
    This fixes the axis's sign only: tied eigenvalues and the background
    eigenbasis can still prevent identical seeded draws across LAPACK builds.
    This default axis also exists in unstructured data; its presence
    alone does not establish a population-stratified genotype model.
    """
    if not 0 <= confounding_strength <= 1:
        raise ValueError("confounding_strength must be in [0, 1]")
    Gd = _trait_genotypes(G)
    rng = np.random.default_rng(seed)
    eigendecomposition = None
    if environment is None:
        K = _grm(Gd) if K is None else _trait_kinship(K, Gd.shape[0])
        lam, U = np.linalg.eigh(K)
        axis = U[:, -1]
        lead = _standardized(axis if axis[np.argmax(np.abs(axis))] >= 0 else -axis)
        eigendecomposition = (lam, U)
    else:
        lead = np.asarray(environment, dtype=float)
        if lead.shape != (Gd.shape[0],):
            raise ValueError("environment must have one value per sample")
        lead = _standardized(lead)
    tr = _simulate_trait(
        Gd,
        h2=h2,
        n_causal=n_causal,
        architecture=architecture,
        effect_dist=effect_dist,
        causal=causal,
        K=K,
        eigendecomposition=eigendecomposition,
        genotype_block_size=genotype_block_size,
        seed=int(rng.integers(1, 2**31 - 1)),
    )
    s = confounding_strength
    structure = np.sqrt(s) * lead
    factor = np.sqrt(1 - s)
    liability = structure + factor * tr["liability"]
    return {
        "y": _standardized(liability),
        "liability": liability,
        "u": factor * tr["u"],
        "q": factor * tr["q"],
        "e": factor * tr["e"],
        "structure": structure,
        "causal": tr["causal"],
        "effects": factor * tr["effects"],
    }


def simulate_gxe_trait(
    G: np.ndarray,
    E: Optional[np.ndarray] = None,
    h2: float = 0.5,
    interaction_h2: float = 0.2,
    n_causal: int = 5,
    seed: Union[int, np.random.Generator, None] = 3,
) -> dict:
    """Trait with genotype-environment interaction effects.

    The component targets are ``h2-interaction_h2`` additive,
    ``interaction_h2`` interaction and ``1-h2`` residual variance.
    ``E`` is standardized (drawn normal if omitted). Correlated components
    and finite samples can change realized variance fractions.
    Returns raw ``liability = u + q + interaction + e`` and standardized
    ``y``. ``effects`` and ``interaction_effects`` multiply standardized
    causal genotypes and their products with E, respectively.
    """
    if not 0 <= interaction_h2 <= h2 <= 1:
        raise ValueError("require 0 <= interaction_h2 <= h2 <= 1")
    if n_causal == 0 and h2 > 0:
        raise ValueError("positive additive or interaction variance requires at least one causal variant")
    rng = np.random.default_rng(seed)
    Gd = _trait_genotypes(G)
    n = Gd.shape[0]
    if E is None:
        E = rng.standard_normal(n)
    E = np.asarray(E, dtype=float)
    if E.shape != (n,):
        raise ValueError("E must have one value per sample")
    E = _standardized(E)

    # Continue the same stream: resetting an integer seed would reuse the
    # environment's innovations in the genetic/residual draw.
    base = simulate_trait(Gd, h2=h2 - interaction_h2, n_causal=n_causal, seed=rng)
    causal = base["causal"]
    inter, inter_effects = np.zeros(n), np.zeros(causal.size)
    if interaction_h2 > 0:
        if causal.size == 0:
            raise ValueError("positive interaction variance requires at least one polymorphic causal variant")
        inter, inter_effects = _scaled_score(
            _standardize_cols(Gd, causal) * E[:, None],
            np.sign(rng.standard_normal(causal.size)), interaction_h2)
    base_residual = 1 - h2 + interaction_h2
    e = base["e"] * np.sqrt((1 - h2) / base_residual) if base_residual > 0 else base["e"]
    liability = base["u"] + base["q"] + inter + e
    base.update(y=_standardized(liability), liability=liability, e=e,
                environment=E, interaction=inter, interaction_effects=inter_effects)
    return base


def simulate_correlated_traits(
    G: np.ndarray,
    h2_a: float = 0.5,
    h2_b: float = 0.5,
    rg: float = 0.6,
    n_causal: int = 20,
    seed: Union[int, np.random.Generator, None] = 4,
) -> dict:
    """Two traits with a target genetic correlation ``rg``.

    Both the infinitesimal and QTL components carry ``rg``: combine each
    trait-A draw with an independent draw using ``rg`` and
    ``sqrt(1-rg**2)``. Each component carries half the genetic variance.
    This targets the total genetic correlation, with finite-sample
    variation; at rg = +/-1 the genetic values are exactly proportional.
    ``g_a``/``g_b`` and ``liability_a``/``liability_b`` expose the raw
    genetic values and liabilities. ``u_a``/``u_b`` retain the unscaled
    background draws; ``y_a``/``y_b`` are standardized phenotypes.
    Causal variants are sampled only from columns that vary across samples,
    capped at the number available, as in :func:`simulate_trait`.
    """
    if not -1.0 <= rg <= 1.0:
        raise ValueError("rg must be in [-1, 1]")
    if not 0 <= h2_a <= 1 or not 0 <= h2_b <= 1:
        raise ValueError("h2_a and h2_b must be in [0, 1]")
    G = _trait_genotypes(G)
    if n_causal == 0:
        raise ValueError(
            "correlated-trait architectures require at least one causal variant")
    rng = np.random.default_rng(seed)
    factor = _background_factor(G)
    shared = _draw_background(factor, 1.0, rng)
    idio = _draw_background(factor, 1.0, rng)
    for v in (shared, idio):
        sd = v.std()
        if not np.isfinite(sd) or sd <= 0:
            raise ValueError(
                "cannot standardize a constant or non-finite component")
        v /= sd

    n = G.shape[0]
    causal = _causal_indices(G, n_causal, None, rng)
    Za = _standardize_cols(G, causal)
    qa, _ = _scaled_score(Za, rng.normal(size=causal.size), 1.0)
    qb, _ = _scaled_score(Za, rng.normal(size=causal.size), 1.0)
    correlated_bg = rg * shared + np.sqrt(1 - rg**2) * idio
    correlated_qtl = rg * qa + np.sqrt(1 - rg**2) * qb

    g_a = np.sqrt(h2_a * 0.5) * (shared + qa)
    g_b = np.sqrt(h2_b * 0.5) * (correlated_bg + correlated_qtl)
    e_a = rng.standard_normal(n) * np.sqrt(1 - h2_a)
    e_b = rng.standard_normal(n) * np.sqrt(1 - h2_b)
    ya = g_a + e_a
    yb = g_b + e_b
    return {
        "y_a": _standardized(ya),
        "y_b": _standardized(yb),
        "u_a": shared,
        "u_b": correlated_bg,
        "g_a": g_a,
        "g_b": g_b,
        "liability_a": ya,
        "liability_b": yb,
        "causal": causal,
    }


def _standardize_cols(G: np.ndarray, idx: np.ndarray) -> np.ndarray:
    Z = np.asarray(G[:, idx], dtype=np.float64)
    Z = Z - Z.mean(axis=0, keepdims=True)
    sd = Z.std(axis=0, keepdims=True)
    return Z / np.where(sd > 0, sd, 1.0)


# --------------------------------------------------------------------------- #
# Ascertainment and case/control bookkeeping
# --------------------------------------------------------------------------- #
def ascertain_case_control(
    trait: Union[dict, np.ndarray],
    n_cases: int,
    n_controls: int,
    seed: Union[int, np.random.Generator, None] = 5,
) -> dict:
    """Sample exact case/control counts from a liability-threshold trait.

    ``trait`` is the dict from :func:`simulate_binary_trait` (or any dict
    with ``case_control`` and ``liability``). Cases and controls are drawn
    without replacement to the *exact* requested counts -- the
    all-cases-plus-k-controls register design and the balanced cohort, the
    two ascertainment schemes every case/control method faces. Returns
    ``{"index", "case_control", "liability"}`` aligned to the sampled
    participants; raises ``ValueError`` when the population holds fewer
    cases or controls than requested.
    """
    if isinstance(trait, dict):
        if "case_control" not in trait or "liability" not in trait:
            raise ValueError(
                "a trait dict must contain 'case_control' and 'liability'")
        cc = np.asarray(trait["case_control"])
        liab = np.asarray(trait["liability"], dtype=float)
    else:
        cc = np.asarray(trait)
        liab = None
    if cc.ndim != 1 or not np.isin(cc, (0, 1)).all():
        raise ValueError("case_control must be a 1-D vector of exact 0/1 values")
    if liab is not None and (liab.shape != cc.shape or not np.isfinite(liab).all()):
        raise ValueError("liability must be finite with one entry per participant")
    cc = cc.astype(np.intp)
    for name, count in (("n_cases", n_cases), ("n_controls", n_controls)):
        if (isinstance(count, (bool, np.bool_))
                or not isinstance(count, (int, np.integer)) or count < 0):
            raise ValueError(f"{name} must be a nonnegative integer")
    cases = np.flatnonzero(cc == 1)
    controls = np.flatnonzero(cc == 0)
    if cases.size < n_cases:
        raise ValueError(
            f"population has {cases.size} cases, {n_cases} requested; "
            "raise the prevalence or the population size"
        )
    if controls.size < n_controls:
        raise ValueError(
            f"population has {controls.size} controls, {n_controls} requested"
        )
    rng = np.random.default_rng(seed)
    index = np.concatenate(
        [
            rng.choice(cases, size=int(n_cases), replace=False),
            rng.choice(controls, size=int(n_controls), replace=False),
        ]
    )
    out = {
        "index": index,
        "case_control": cc[index],
    }
    if liab is not None:
        out["liability"] = liab[index]
    return out


def n_eff_case_control(n_case, n_control):
    """Effective sample size of a case/control GWAS: ``4/(1/N_case + 1/N_control)``.

    Equals total N for a balanced study and tends to ``4*N_case`` with
    many more controls than cases. Accepts scalars or arrays.
    """
    n_case = np.asarray(n_case, dtype=float)
    n_control = np.asarray(n_control, dtype=float)
    if np.any(n_case <= 0) or np.any(n_control <= 0):
        raise ValueError("n_case and n_control must be positive")
    n = 4.0 / (1.0 / n_case + 1.0 / n_control)
    return float(n) if n.ndim == 0 else n


def h2_liability(h2_observed, prevalence, *, prop_cases=None):
    """Convert observed-scale SNP h² to the liability scale (Lee et al. 2011).

    For population prevalence ``K``, study case fraction ``P``, threshold
    ``t = Phi^-1(1-K)``, and ``z = phi(t)``::

        h²_liab = h²_obs * [K(1-K)]² / (z² * P(1-P))

    ``prop_cases=None`` defaults ``P`` to ``K`` **with a warning**: that
    default is correct only when the study sample mirrors the population
    case fraction. A balanced case/control GWAS of a 1% trait must pass
    ``prop_cases=0.5``; the silent default overstates h² there by
    ``P(1-P)/(K(1-K))`` (about 25x).
    """
    K = float(prevalence)
    if not 0.0 < K < 1.0:
        raise ValueError("prevalence must be in (0, 1)")
    if prop_cases is None:
        warnings.warn(
            "h2_liability(prop_cases=None) assumes the study sample mirrors the "
            "population case fraction (P=K). For an ascertained case/control "
            "study pass the GWAS case fraction explicitly (balanced: "
            "prop_cases=0.5); the silent default overstates liability h2 by "
            "P(1-P)/(K(1-K)) there.", UserWarning, stacklevel=2)
    P = K if prop_cases is None else float(prop_cases)
    if not 0.0 < P < 1.0:
        raise ValueError("prop_cases must be in (0, 1)")
    nd = NormalDist()
    # Avoid forming 1-K: Phi^-1(1-K) = -Phi^-1(K), including tiny K.
    t = -nd.inv_cdf(K)
    z = nd.pdf(t)
    factor = (K * (1.0 - K)) ** 2 / (z * z * P * (1.0 - P))
    h2 = np.asarray(h2_observed, dtype=float) * factor
    return float(h2) if h2.ndim == 0 else h2
