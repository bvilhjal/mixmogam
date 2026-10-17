"""Two-step mixed-model association: BOLT-LMM and HRATT.

Both methods first build, for every LOCO group (chromosome), a residual
phenotype from which the polygenic effects of the *other* chromosomes have
been conditioned out, and then test each SNP against the residual of its
own LOCO group with a retrospective score statistic,

    chi2_j = n_eff (z_j' w)^2 / (z_j' z_j * |w|^2),

calibrated by a single genome-wide constant. They differ in the residual
and in the calibration:

- BOLT-LMM-inf (Loh et al. 2015, Nat Genet): w = V_{-g}^{-1} y by
  conjugate gradients; the constant matches the exact prospective
  statistic (z' V^{-1} y)^2 / (z' V^{-1} z) at 30 random SNPs.
- BOLT-LMM: w = y minus the variational-Bayes posterior-mean prediction
  under a two-Gaussian mixture prior (hyperparameters by 5-fold CV over
  BOLT's 18-point grid); calibrated so its LD Score regression intercept
  matches BOLT-LMM-inf's. Falls back to BOLT-LMM-inf when the mixture
  does not predict better in CV.
- HRATT, the Heritability-weighted Residual Association Two-step Test:
  this package's method inspired by LDAK-KVIK (Hof & Speed 2025, Nat
  Genet) and following its design: w = y minus an elastic-net
  LOCO polygenic score under the LDAK-Thin heritability model, cross-fitted
  as in REGENIE with each fold's predictor trained on the other folds;
  lambda = 1 unless a test of inter-chromosome correlation finds strong
  structure, in which case lambda is matched to the GRAMMAR-Gamma-calibrated
  ridge statistic on weakly associated SNPs.

The single calibration constant rests on the assumption that the
prospective denominator z_j' V^{-1} z_j is proportional to z_j' z_j across
SNPs. BOLT-LMM's authors flag it as possibly invalid outside human data.
Every result here therefore reports the spread of the ratio over the
calibration SNPs (``calibration_cv``): with strong structure it is large,
and a single constant cannot calibrate every SNP.

HRATT also tests case-control outcomes and sampling-weighted samples: the
score z_j'a of a weighted residual a against its variance when genotypes
are drawn given the covariates, zvar_j |a|^2 rho_j, with saddlepoint tails
(see :func:`hratt`).

Deviations from the reference implementations are documented on each
function. Variance components default to this package's stochastic Lanczos
REML (HE with sampling weights); HRATT also offers projected
single-component randomized HE, distinct from LDAK-KVIK's partitioned HE and
optional Monte Carlo REML.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import warnings

import numpy as np
from scipy import linalg, stats

from mixmogam._binary import null_logistic
from mixmogam._cg import SpectralPreconditioner, batched_pcg
from mixmogam._loco import LocoGenotypes, loco_groups
from mixmogam._slq import randomized_eigh_op
from mixmogam._vb import PRIOR_ENET, PRIOR_MIXTURE, VBEngine, _validate_n_threads
from mixmogam.ldscore import ld_scores, ldsc_intercept
from mixmogam.lmm import LMM, _design_matrix
from mixmogam.results import GwasResult

__all__ = ["bolt_inf", "bolt", "hratt", "structure_test", "BOLT_GRID", "HRATT_GRID",
           "HRATT_ALPHAS"]

# BOLT-LMM: f2 = spike share of the prior variance, p = slab probability
BOLT_GRID = [(f2, p) for f2 in (0.5, 0.3, 0.1) for p in (0.5, 0.2, 0.1, 0.05, 0.02, 0.01)]
# LDAK-KVIK's elastic-net grid: (p, F) = (lasso share, ridge variance share)
HRATT_GRID = [(0.5, 0.5), (0.5, 0.3), (0.5, 0.1), (0.1, 0.5), (0.1, 0.3), (0.1, 0.1),
             (0.01, 0.5), (0.01, 0.3), (0.01, 0.1), (0.0, 1.0)]
HRATT_ALPHAS = (-1.0, -0.75, -0.5, -0.25, 0.0)


# ----------------------------------------------------------------------
# Shared machinery
# ----------------------------------------------------------------------


@dataclass
class _Setup:
    y: np.ndarray  # phenotype (working response of the row-scaled model)
    y_p: np.ndarray  # covariate-projected phenotype
    X: np.ndarray  # design incl. intercept (row-scaled with the working weights)
    lg: LocoGenotypes
    labels: list
    n_eff: int  # n - q
    w: Optional[np.ndarray] = None  # sampling weights, mean 1 (None: unweighted)
    v: Optional[np.ndarray] = None  # working weights, mean 1 (None: constant)
    s: Optional[np.ndarray] = None  # row scale sqrt(v) (None: unscaled)
    y_raw: Optional[np.ndarray] = None  # the phenotype as given
    X_raw: Optional[np.ndarray] = None  # the unscaled design
    n_struct: Optional[float] = None  # Kish size of v, less q: the structure test's sample size
    mu0_clipped: int = 0  # binary: working probabilities clipped in step 1
    Q_raw: Optional[np.ndarray] = None  # orthonormal basis of the unscaled design


def _normalize_weights(sample_weights, n: int) -> Optional[np.ndarray]:
    """Sampling weights scaled to mean one (None stays None)."""
    if sample_weights is None:
        return None
    w = np.asarray(sample_weights, dtype=np.float64)
    if w.shape != (n,) or not np.isfinite(w).all() or np.any(w <= 0):
        raise ValueError("sample_weights must hold one finite positive weight per sample; "
                         "drop samples with zero weight first")
    return w / w.mean()


def _setup(y, gt, X, max_loco_groups, block, n_threads=1, *, sample_weights=None,
           trait: str = "quantitative", mu_clip: float = 1e-4) -> _Setup:
    """The (row-scaled) working model shared by the two-step steps.

    Without weights and for quantitative traits this is the phenotype, the
    design and the genotypes as given. Sampling weights w and, for binary
    traits, the logistic working weights mu0 (1 - mu0) of the covariate-only
    null fit define working weights v (mean one); the model is then fitted
    in row-scaled form, y~ = s y, X~ = s X and z~ = s z with s = sqrt(v), so
    every solver keeps its K + delta I form. The binary working response is
    (y - mu0) / (mu0 (1 - mu0)), whose scaled form is the Pearson residual.
    """
    y = np.asarray(y, dtype=np.float64).ravel()
    if not np.isfinite(y).all():
        raise ValueError("y contains NaN or infinite values")
    if y.size != gt.n_samples:
        raise ValueError("y and genotypes disagree on the sample count")
    Xd = _design_matrix(X, y.size, add_intercept=True)
    w = _normalize_weights(sample_weights, y.size)
    clipped = 0
    if trait == "quantitative":
        v, y_work = w, y
    elif trait == "binary":
        if not np.all((y == 0.0) | (y == 1.0)):
            raise ValueError("a binary trait needs outcomes coded 0 (control) and 1 (case)")
        null = null_logistic(y, Xd, weights=w)
        mu = np.clip(null["mu"], mu_clip, 1.0 - mu_clip)
        clipped = int(np.count_nonzero(mu != null["mu"]))
        if clipped:
            warnings.warn(f"{clipped} working probabilities of the null logistic fit were clipped to "
                          f"[{mu_clip:g}, {1 - mu_clip:g}] for step 1", RuntimeWarning, stacklevel=3)
        variance = mu * (1.0 - mu)
        v = variance if w is None else w * variance
        y_work = (y - mu) / variance
    else:
        raise ValueError("trait must be 'quantitative' or 'binary'")
    s = None
    if v is not None:
        v = v / v.mean()
        if np.ptp(v) > 0:
            s = np.sqrt(v)
        else:  # constant working weights: the unscaled model, exactly
            v = None
    Xs = Xd if s is None else Xd * s[:, None]
    Q = linalg.qr(Xs, mode="economic")[0]
    groups, labels = loco_groups(gt.chromosome, max_loco_groups)
    if len(labels) < 2:
        raise ValueError("LOCO needs variants on at least two chromosomes/groups")
    lg = LocoGenotypes(gt, groups, Q, block=block, n_threads=n_threads, row_scale=s)
    ys = y_work if s is None else s * y_work
    y_p = ys - Q @ (Q.T @ ys)
    q = Xd.shape[1]
    n_struct = None if v is None else float(np.sum(v) ** 2 / np.sum(v * v)) - q
    return _Setup(y=ys, y_p=y_p, X=Xs, lg=lg, labels=labels, n_eff=y.size - q,
                  w=w, v=v, s=s, y_raw=y, X_raw=Xd, n_struct=n_struct,
                  mu0_clipped=clipped,
                  Q_raw=Q if s is None else linalg.qr(Xd, mode="economic")[0])


class _KOp:
    """Kinship operator adapter for :class:`LMM` (K = Z' W Z / sum W)."""

    def __init__(self, lg: LocoGenotypes, weights=None):
        self.lg = lg
        self.n = lg.n
        self.weights = weights
        self._trace = None

    @property
    def trace(self):
        # The variants' projected squared norms are prepared: no genotype pass.
        if self._trace is None and self.weights is None:
            self._trace = self.lg.trace
        elif self._trace is None:
            self._trace = float(self.weights @ self.lg.zz) / float(np.sum(self.weights))
        return self._trace

    def matmul(self, P):
        return self.lg.matmul(P, self.weights)


def fit_variance_components(st: _Setup, weights=None, random_state=0, basis=None,
                            **slq) -> object:
    """REML variance components with the stochastic Lanczos solver.

    The kinship is Z' W Z / sum(W) over the covariate-projected
    standardized genotypes, applied as a streaming operator. BOLT-LMM uses
    Monte Carlo REML; LDAK-KVIK starts with partitioned randomized HE and
    optionally refines its heritability by Monte Carlo REML. ``basis``
    (:func:`_spectral_basis`, same weights) supplies the eigenpairs the
    solver deflates, which a later CG preconditioner then reuses.
    """
    op = _KOp(st.lg, weights) if basis is None else basis["op"]
    lmm = LMM(st.y, X=st.X, K=op, add_intercept=False, random_state=random_state)
    # Two-step methods use variance components, never the null-model GLS
    # coefficients or prediction methods of a complete LMFit.
    deflation = None if basis is None else (basis["values"], basis["vectors"])
    return lmm._fit_variance_components(solver="slq", deflation_basis=deflation, **slq)


SPECTRAL_BASIS_K = 128  # the REML deflation width; preconditioners use the top 64


def _spectral_basis(st: _Setup, weights=None, random_state=0) -> dict:
    """Top eigenpairs of the (weighted) kinship, computed once.

    The genotype blocks are covariate-projected, so the kinship already
    equals the REML operator S K S: the Lanczos REML deflation and the
    conjugate-gradient preconditioner can share one randomized subspace
    iteration (the deflation's own: width 128, seed ``random_state``)
    instead of each running its own passes over the genotypes.
    """
    op = _KOp(st.lg, weights)
    values, vectors = randomized_eigh_op(op.matmul, st.lg.n,
                                         min(SPECTRAL_BASIS_K, st.lg.n - 1),
                                         random_state=random_state)
    return {"op": op, "values": values, "vectors": vectors}


def _preconditioner(st: _Setup, basis: dict) -> SpectralPreconditioner:
    return SpectralPreconditioner.from_eigenpairs(basis["values"], basis["vectors"], st.lg.n,
                                                  basis["op"].trace, k=min(64, st.lg.n - 2))


def _loco_solve(st: _Setup, delta: float, rhs: np.ndarray, col_group: np.ndarray,
                weights=None, tol: float = 1e-6, pre=None) -> tuple[np.ndarray, dict]:
    """(K_{-g} + delta I)^{-1} rhs[:, r] for g = col_group[r], batched."""
    if pre is None:
        op = _KOp(st.lg, weights)
        pre = SpectralPreconditioner(op.matmul, st.lg.n, op.trace, k=min(64, st.lg.n - 2))

    def apply(P, cols):
        return st.lg.matmul_loco(P, col_group[cols], weights) + delta * P

    solution, info = batched_pcg(apply, rhs, lambda R: pre(R, delta), tol=tol, max_iter=2000)
    if not info["converged"]:
        raise RuntimeError("LOCO conjugate-gradient solve did not converge; association statistics unavailable")
    return solution, info


def _allele_probabilities(cu, zz, Q, mean, sd):
    """F-shrunk covariate predictions of allele frequencies, in (0, 1)."""
    n, q = Q.shape
    fit2 = np.einsum("ij,ij->i", cu, cu)
    with np.errstate(divide="ignore", invalid="ignore"):
        F = fit2 * (n - q) / ((q - 1) * np.maximum(zz - fit2, 0.0))
        shrink = np.where(np.isfinite(F), np.clip(1.0 - 1.0 / F, 0.0, 1.0), 1.0)
    shrink[~(fit2 > 0)] = 0.0
    p = (cu @ Q.T) * (0.5 * sd * shrink)[:, None] + (0.5 * mean)[:, None]
    np.clip(p, 0.5 / n, 1.0 - 0.5 / n, out=p)
    return p


def _ancestry_ratio(cu: np.ndarray, zz: np.ndarray, Q: np.ndarray, mean: np.ndarray,
                    sd: np.ndarray, a2: np.ndarray) -> np.ndarray:
    """Covariate-specific genotype variance of retrospective scores, as a
    ratio to the pooled one (rows: variants).

    Sample i's allele frequency p_ij is half of variant j's genotype count
    fitted on the covariates, the fit's deviations from the mean shrunk by
    (F - 1) / F (zero if F < 1) with F its regression F statistic, so that
    covariates unrelated to the genotype leave the pooled frequency; it is
    kept in [1/(2n), 1 - 1/(2n)]. ``cu`` = z_j'Q are the standardized
    genotypes' projections on the orthonormal design basis ``Q`` (the
    intercept first), ``zz`` their sums of squares, and ``mean`` and ``sd``
    undo the standardization. pi_ij = p_ij (1 - p_ij) is sample i's binomial
    genotype variance up to a constant, as in SPAmix (Ma et al. 2025).
    Returns rho_j = (sum_i a_i^2 pi_ij / sum_i a_i^2) / mean_i pi_ij, the
    factor by which the variance of sum_i a_i z_ij departs from
    zvar_j |a|^2 when allele frequencies follow the covariates (ancestry
    principal components among them) and ``a2`` (a_i^2) varies with
    ancestry.
    """
    if Q.shape[1] < 2:
        return np.ones(cu.shape[0])
    p = _allele_probabilities(cu, zz, Q, mean, sd)
    p *= 1.0 - p
    return (p @ a2) / (float(np.sum(a2)) * np.mean(p, axis=1))


def _retro_stats(st: _Setup, W: np.ndarray, retrospective: bool = False) -> dict:
    """Uncalibrated retrospective chi2 of every variant against its group's
    column of W: n_eff (z'w)^2 / (z'z |w|^2).

    ``retrospective`` adds the test of a row-scaled (sampling-weighted)
    model. With a_g = s (I - QQ') w_g, the weighted residual (orthogonal to
    the covariates), the score z_j'a_g has variance zvar_j |a_g|^2 when the
    genotypes are exchangeable given the covariates, whatever the weights;
    zvar_j = (sum z^2 - |Q_X' z|^2) / n_eff is the unweighted genotype
    variance after covariates. When allele frequencies vary with the
    covariates the variance is zvar_j |a_g|^2 rho_j instead
    (:func:`_ancestry_ratio`). ``chi2_retro`` = (z_j'a_g)^2 / ``V`` with
    ``V`` = zvar_j |a_g|^2 rho_j, and ``A`` holds the a_g.
    """
    m = st.lg.m
    num = np.zeros(m)
    zz = np.array(st.lg.zz)
    wn = np.einsum("ij,ij->j", W, W)
    Wp = st.lg.project(W)  # Z_j w = z_j (I - QQ') w
    if retrospective:
        A = Wp if st.s is None else Wp * st.s[:, None]
        A2 = A * A
        inv_s = None if st.s is None else 1.0 / st.s
        zvar, rho = np.zeros(m), np.zeros(m)
    # Widen at most 16 MiB of variant rows, retaining each full sample
    # reduction. The minimum one-row tile and native BLAS scratch are separate.
    tile = max(1, (16 * 1024**2) // max(8 * st.lg.n, 1))
    for idx, g, Z in st.lg.raw_slices(tile):
        Z64 = Z.astype(np.float64, copy=False)
        num[idx] = Z64 @ Wp[:, g]
        if retrospective:
            z = Z64 if inv_s is None else Z64 * inv_s
            cu = z @ st.Q_raw
            zsq = np.einsum("ij,ij->i", z, z)
            zvar[idx] = (zsq - np.einsum("ij,ij->i", cu, cu)) / st.n_eff
            rho[idx] = _ancestry_ratio(cu, zsq, st.Q_raw, st.lg.mean[idx], st.lg.sd[idx], A2[:, g])
    wg = wn[st.lg.groups]
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = st.n_eff * num * num / (zz * wg)
    chi2[~(zz > 0)] = np.nan
    out = {"chi2": chi2, "num": num, "zz": zz, "wnorm2": wg}
    if retrospective:
        V = zvar * np.einsum("ij,ij->j", A, A)[st.lg.groups] * rho
        with np.errstate(divide="ignore", invalid="ignore"):
            chi2_retro = num * num / V
        chi2_retro[~((zz > 0) & (V > 0))] = np.nan
        out.update({"V": V, "chi2_retro": chi2_retro, "A": A, "rho": rho})
    return out


def _genotype_spa(st: _Setup, A: np.ndarray, idx: np.ndarray, U: np.ndarray, V: np.ndarray,
                  lam: float, n_threads: int, two_sided: str = "distance") -> tuple:
    """Retrospective saddlepoint p-values (and logs) of scores
    ``U`` = z_j'a_g for variants ``idx``.

    With covariates, genotypes are independent Binomial(2, p_ij) draws at
    the covariate-predicted allele frequencies, conditional on the observed
    missing-call mask. This assumes Hardy-Weinberg equilibrium conditional
    on those covariates. An intercept-only design retains the empirical
    genotype distribution. Scores are first centred and scaled by
    sqrt(lam K''(0) / V), so that the bulk agrees with the calibrated
    normal statistic; ``two_sided`` as in :func:`mixmogam._spa.spa_pvalue`.
    """
    from mixmogam._spa import spa_pvalue

    n = st.lg.n
    p = np.full(idx.size, np.nan)
    log_p = np.full(idx.size, np.nan)
    groups = st.lg.groups[idx]
    for g in np.unique(groups):
        a = np.ascontiguousarray(A[:, g])
        a2 = float(a @ a)
        sel = np.flatnonzero(groups == g)
        # Bound the (variants, samples, 4) probability tile to 8 MiB.
        tile = min(256, max(1, (8 * 1024**2) // (32 * n)))
        for start in range(0, sel.size, tile):
            rows = sel[start : start + tile]
            z = st.lg._decode(idx[rows], scaled=False)
            if st.Q_raw.shape[1] > 1:
                z64 = z.astype(np.float64)
                js = idx[rows]
                af = _allele_probabilities(z64 @ st.Q_raw, np.einsum("ij,ij->i", z64, z64),
                                            st.Q_raw, st.lg.mean[js], st.lg.sd[js])
                called = np.asarray(st.lg.gt.G[:, js]).T != -1
                freqs = np.zeros((rows.size, n, 4))
                freqs[:, :, 0] = (1 - af) ** 2 * called
                freqs[:, :, 1] = 2 * af * (1 - af) * called
                freqs[:, :, 2] = af**2 * called
                freqs[:, :, 3] = ~called
                values = np.column_stack([np.zeros(rows.size), np.ones(rows.size),
                                           np.full(rows.size, 2), st.lg.mean[js]])
                variance = (2 * af * (1 - af) * called) @ (a * a)
                raw_score = U[rows] * st.lg.sd[js] + st.lg.mean[js] * a.sum()
                ratio = lam * variance / (V[rows] * st.lg.sd[js] ** 2)
                p[rows], log_p[rows] = spa_pvalue(raw_score, a, values, freqs, var_ratio=ratio,
                                                  two_sided=two_sided, n_threads=n_threads)
                continue
            values, freqs = np.zeros((rows.size, 4)), np.zeros((rows.size, 4))
            for r in range(rows.size):
                vals, counts = np.unique(z[r], return_counts=True)
                values[r, : vals.size], freqs[r, : vals.size] = vals, counts / n
            mean = np.sum(values * freqs, axis=1, keepdims=True)
            var = np.sum((values - mean) ** 2 * freqs, axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = lam * a2 * var / V[rows]
            p[rows], log_p[rows] = spa_pvalue(U[rows], a, values, freqs, var_ratio=ratio,
                                              two_sided=two_sided, n_threads=n_threads)
    return p, log_p


def _genotype_tails(st: _Setup, A: np.ndarray, U: np.ndarray, V: np.ndarray,
                    chi2_normal: np.ndarray, lam: float, threshold: float, n_threads: int,
                    two_sided: str = "distance") -> tuple:
    """:func:`_saddlepoint_tails` with the genotype saddlepoint of scores
    ``U`` against the coefficients ``A`` (one column per group)."""

    def spa(idx):
        return _genotype_spa(st, A, idx, U[idx], V[idx], lam, n_threads, two_sided)

    return _saddlepoint_tails(chi2_normal, threshold, spa)


def _saddlepoint_tails(chi2_normal: np.ndarray, threshold: float, spa) -> tuple:
    """p-values and chi2 statistics with saddlepoint tails above |z| =
    ``threshold``: ``spa(idx)`` returns (p, log p) of those variants."""
    from scipy.special import ndtri_exp

    p = stats.chi2.sf(chi2_normal, 1)
    chi2 = chi2_normal.copy()
    tail = np.flatnonzero(np.isfinite(chi2_normal) & (chi2_normal > threshold * threshold))
    n_spa = 0
    if tail.size:
        p_spa, log_p = spa(tail)
        ok = np.isfinite(log_p)
        p[tail[ok]] = p_spa[ok]
        chi2[tail[ok]] = ndtri_exp(log_p[ok] - np.log(2.0)) ** 2
        n_spa = int(ok.sum())
    return p, chi2, n_spa


CALIBRATION_OVERSAMPLE = 1.25


def _calibration_draw(st: _Setup, n_cal: int, rng) -> np.ndarray:
    """Calibration candidates, drawn before the LOCO residuals exist.

    BOLT-LMM-inf calibrates on random SNPs whose GRAMMAR chi2, a function of
    the residuals, is below 5. Polymorphic SNPs drawn in random order and
    filtered afterwards (rejection sampling: the same distribution) let the
    candidates' prospective solves share the residuals' CG passes. About a
    quarter more are drawn than needed; a short second solve covers a
    shortfall.
    """
    if not isinstance(n_cal, (int, np.integer)) or n_cal < 1:
        raise ValueError("n_calibration must be a positive integer")
    pool = np.nonzero(st.lg.sd > 0)[0]
    if pool.size == 0:
        raise ValueError("calibration needs polymorphic SNPs outside the covariate span")
    size = min(pool.size, int(np.ceil(CALIBRATION_OVERSAMPLE * n_cal)) + 4)
    return rng.choice(pool, size=size, replace=False)


def _solve_with_calibration(st: _Setup, delta: float, rhs: np.ndarray, col_group: np.ndarray,
                            cand: np.ndarray, weights=None, tol: float = 1e-6, pre=None):
    """LOCO solves for ``rhs`` and for the candidates' genotypes, batched.

    Returns the ``rhs`` solutions, the candidates' rows and solutions, and
    the CG information of the shared run.
    """
    Zc = st.lg.rows(cand)
    groups = np.concatenate([np.asarray(col_group, dtype=np.int64), st.lg.groups[cand]])
    X, info = _loco_solve(st, delta, np.column_stack([rhs, Zc.T]), groups,
                          weights=weights, tol=tol, pre=pre)
    r = rhs.shape[1]
    return X[:, :r], {"cand": cand, "Z": Zc, "V": X[:, r:]}, info


def _calibrate_inf(st: _Setup, vg: float, delta: float, rs: dict, n_cal: int, rng,
                   drawn: dict, weights=None, pre=None, tol: float = 1e-6) -> dict:
    """BOLT-LMM-inf calibration constant from exact prospective statistics.

    For ``n_cal`` random SNPs with GRAMMAR chi2 < 5, c_j = D_prosp / D_retro
    with D_prosp = vg z'(K_{-g} + delta I)^{-1} z and D_retro = (z'z / n_eff)
    |w|^2; the constant is mean(c_j). ``cv`` = sd(c_j) / mean(c_j) measures
    how far the proportional-denominator assumption is from holding.
    ``drawn`` holds the candidates solved together with the residuals
    (:func:`_solve_with_calibration`); the first ``n_cal`` that qualify are
    used, as from uniform sampling among all qualifying SNPs.
    """
    qualifies = (rs["chi2"] < 5.0) & (rs["zz"] > 0)
    if not qualifies.any():
        qualifies = rs["zz"] > 0  # as before: no SNP below 5, use any
    if not qualifies.any():
        raise ValueError("calibration needs polymorphic SNPs outside the covariate span")
    cand, Zc, V = drawn["cand"], drawn["Z"], drawn["V"]
    keep = np.nonzero(qualifies[cand])[0][:n_cal]
    sel, Zsel, Vsel = cand[keep], Zc[keep], V[:, keep]
    short = min(n_cal, int(qualifies.sum())) - sel.size
    if short > 0:
        rest = np.setdiff1d(np.nonzero(qualifies)[0], cand)
        extra = rng.choice(rest, size=min(short, rest.size), replace=False)
        Ze = st.lg.rows(extra)
        Ve, _ = _loco_solve(st, delta, Ze.T, st.lg.groups[extra], weights=weights, tol=tol, pre=pre)
        sel, Zsel, Vsel = np.r_[sel, extra], np.vstack([Zsel, Ze]), np.hstack([Vsel, Ve])
    order = np.argsort(sel)
    sel, Zsel, Vsel = sel[order], Zsel[order], Vsel[:, order]
    d_prosp = vg * np.einsum("ij,ji->i", Zsel, Vsel)
    d_retro = rs["zz"][sel] / st.n_eff * rs["wnorm2"][sel]
    ratio = d_prosp / d_retro
    c = float(np.mean(ratio))
    return {"c": c, "cv": float(np.std(ratio) / c), "ratios": ratio, "snps": sel,
            "d_prosp": d_prosp}


def _loco_eigh(st: _Setup, k: int, n_iter: int = 4, oversampling: int = 12,
               random_state=0, weights=None) -> list:
    """Top-k eigenpairs of every LOCO kinship K_{-g}, batched.

    Randomized subspace iteration (Halko et al. 2011) run for all groups at
    once: each pass over the genotypes applies every K_{-g} to its own
    block of columns (:meth:`LocoGenotypes.matmul_loco`). Returns one
    ``(values, vectors, lam_bar)`` per group; ``lam_bar`` is the mean of the
    remaining eigenvalues, (trace K_{-g} - sum values) / (n - k).
    """
    lg = st.lg
    n, G = lg.n, lg.n_groups
    k = max(1, min(k, n - 2))
    L = min(k + oversampling, n)
    col_group = np.repeat(np.arange(G), L)
    rng = np.random.default_rng(random_state)

    def orth(Y):
        return np.concatenate(
            [np.linalg.qr(Y[:, g * L : (g + 1) * L], mode="reduced")[0] for g in range(G)],
            axis=1)

    Q = orth(lg.matmul_loco(rng.standard_normal((n, G * L)), col_group, weights))
    for _ in range(n_iter - 1):
        Q = orth(lg.matmul_loco(Q, col_group, weights))
    KQ = lg.matmul_loco(Q, col_group, weights)
    w = np.ones(lg.m) if weights is None else np.asarray(weights, dtype=np.float64)
    sq = np.bincount(lg.groups, weights=w * lg.zz, minlength=G)
    wsum = np.bincount(lg.groups, weights=w, minlength=G)
    out = []
    for g in range(G):
        Qg, KQg = Q[:, g * L : (g + 1) * L], KQ[:, g * L : (g + 1) * L]
        B = Qg.T @ KQg
        vals, vecs = np.linalg.eigh(0.5 * (B + B.T))
        order = np.argsort(vals)[::-1][:k]
        vals = np.maximum(vals[order], 0.0)
        trace_g = (sq.sum() - sq[g]) / (wsum.sum() - wsum[g])
        out.append((vals, Qg @ vecs[:, order], max((trace_g - vals.sum()) / (n - k), 0.0)))
    return out


def _spectral_quadform(st: _Setup, bases: list, delta: float) -> np.ndarray:
    """z_j' M_g z_j for every variant, with its group's
    M_g = U (L + d)^{-1} U' + (I - U U') / (lam_bar + d).

    The low-rank-plus-flat approximation of z'(K_{-g} + delta I)^{-1} z:
    exact on the leading (structure) eigenvectors of the LOCO kinship, a
    constant on its bulk. ``bases`` holds one ``(values, vectors, lam_bar)``
    per LOCO group (or a single entry used for every group).
    """
    out = np.zeros(st.lg.m)
    # Z_j u = z_j (I - QQ') u. Slices arrive in group order, so one projected
    # basis is held at a time; rows are widened in tiles of at most 32 MiB.
    tile = max(1, (32 * 1024**2) // max(8 * st.lg.n, 1))
    current, basis = None, None
    for idx, g, Z in st.lg.raw_slices():
        key = g if len(bases) > 1 else 0
        if key != current:
            current, basis = key, None  # release the previous projection first
            vals, U, lam_bar = bases[key]
            basis = (vals, st.lg.project(U), lam_bar)
        vals, U, lam_bar = basis
        for start in range(0, idx.size, tile):
            rows = idx[start : start + tile]
            P = Z[start : start + tile].astype(np.float64) @ U
            top = np.einsum("ij,ij->i", P, P)
            out[rows] = (P * P) @ (1.0 / (vals + delta)) + (st.lg.zz[rows] - top) / (lam_bar + delta)
    return out


def _result(st: _Setup, gt, chi2, beta_z, se_z, extra: dict, p=None) -> GwasResult:
    if any(extra.get(key) is False for key in ("cv_converged", "loco_converged")):
        warnings.warn("variational fit did not converge; inspect cv_converged/loco_converged before using results",
                      RuntimeWarning, stacklevel=3)
    sd = st.lg.sd
    with np.errstate(divide="ignore", invalid="ignore"):
        beta = np.where(sd > 0, beta_z / sd, np.nan)
        se = np.where(sd > 0, se_z / sd, np.nan)
    if p is None:
        p = stats.chi2.sf(chi2, 1)
    res = GwasResult(
        chromosome=np.asarray(gt.chromosome),
        position=np.asarray(gt.position),
        p=p,
        variant_ids=np.asarray(gt.variant_ids),
        f_stat=chi2,
        beta=beta,
        se=se,
        af=st.lg.mean / 2.0,
        effect_allele=gt.allele1, other_allele=gt.allele2,
        n=st.lg.n,
    )
    res.extra.update({"statistic": "chi2", "n": st.lg.n,
                      "n_loco_groups": st.lg.n_groups, "loco_groups": st.labels})
    res.extra.update(extra)
    return res


SPECTRAL_WIDTHS = (64, 128, 256, 512)
SPECTRAL_CV_TARGET = 0.03


def _spectral_denominator(st: _Setup, vg: float, delta: float, cal: dict,
                          n_spectral="auto", weights=None) -> dict:
    """Structure-aware denominators D_j = c' vg z' M_g z (mixmogam extension).

    M_g is the low-rank-plus-flat approximation of (K_{-g} + delta I)^{-1}
    on the top-k eigenvectors of the SNP's own LOCO kinship
    (:func:`_loco_eigh`, :func:`_spectral_quadform`); c' is calibrated on the
    same prospective SNPs as BOLT-LMM-inf's constant. The LOCO basis
    matters: a basis from the full kinship contains each SNP's own LD
    block, which reintroduces proximal contamination into the denominator.

    ``n_spectral="auto"`` widens k through 64, 128, 256, 512 until the
    spread of the ratio over the calibration SNPs falls below 3%: the bulk
    of an LD-rich kinship is not flat, and the calibration SNPs measure how
    much of it the basis misses. An integer fixes k.
    """
    widths = (SPECTRAL_WIDTHS if n_spectral == "auto" else (int(n_spectral),))
    if min(widths) < 1:
        raise ValueError("n_spectral must be a positive integer or 'auto'")
    widths = tuple(dict.fromkeys(min(k, st.lg.n - 2) for k in widths))
    best = None
    for k in widths:
        bases = _loco_eigh(st, k, weights=weights)
        A = _spectral_quadform(st, bases, delta)
        ratio = cal["d_prosp"] / (vg * A[cal["snps"]])
        c_s = float(np.mean(ratio))
        best = {"D": c_s * vg * A, "c": c_s, "cv": float(np.std(ratio) / c_s),
                "ratios": ratio, "k": k}
        if best["cv"] <= SPECTRAL_CV_TARGET:
            break
    return best


def _inf_core(st: _Setup, fit, n_cal: int, rng, tol: float,
              denominator: str = "constant", n_spectral="auto", basis=None) -> dict:
    """BOLT-LMM steps 1a-1b on a prepared setup.

    ``denominator="constant"`` is BOLT-LMM-inf: D_j = c (z'z / n_eff) |w|^2.
    ``"spectral"`` uses :func:`_spectral_denominator`. ``kappa`` = D_constant
    / D is the per-SNP factor that turns the constant-calibrated statistic
    into the structure-aware one.
    """
    if denominator not in ("constant", "spectral"):
        raise ValueError(f"unknown denominator {denominator!r}")
    G = st.lg.n_groups
    if basis is not None:
        pre = _preconditioner(st, basis)
    else:
        op = _KOp(st.lg)
        pre = SpectralPreconditioner(op.matmul, st.lg.n, op.trace, k=min(64, st.lg.n - 2))
    cand = _calibration_draw(st, n_cal, rng)
    U, drawn, info = _solve_with_calibration(st, fit.delta, np.repeat(st.y_p[:, None], G, axis=1),
                                             np.arange(G), cand, tol=tol, pre=pre)
    rs = _retro_stats(st, U)
    cal = _calibrate_inf(st, fit.vg, fit.delta, rs, n_cal, rng, drawn, pre=pre, tol=tol)
    D_const = cal["c"] * rs["zz"] / st.n_eff * rs["wnorm2"]
    if denominator == "constant":
        D = D_const
    else:
        sp = _spectral_denominator(st, fit.vg, fit.delta, cal, n_spectral)
        cal = {**cal, "c_constant": cal["c"], "cv_constant": cal["cv"],
               "c": sp["c"], "cv": sp["cv"], "ratios": sp["ratios"], "k": sp["k"]}
        D = sp["D"]
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = rs["num"] ** 2 / D
        # prospective-equivalent effect sizes: D approximates vg z'V^-1 z
        beta_z = rs["num"] * fit.vg / D
        se_z = fit.vg / np.sqrt(D)
        kappa = D_const / D
    chi2[~(rs["zz"] > 0)] = np.nan
    return {"U": U, "rs": rs, "cal": cal, "chi2": chi2, "beta_z": beta_z,
            "se_z": se_z, "kappa": kappa, "cg_iterations": info["iterations"]}


def _fit_extra(fit) -> dict:
    return {"delta": fit.delta, "vg": fit.vg, "ve": fit.ve,
            "pseudo_heritability": fit.pseudo_heritability}


# ----------------------------------------------------------------------
# BOLT-LMM
# ----------------------------------------------------------------------


def bolt_inf(y, gt, X=None, *, max_loco_groups: int = 25, n_calibration: int = 30,
             denominator: str = "constant", n_spectral="auto",
             cg_tol: float = 1e-6, random_state: int = 0,
             block: int = 4096) -> GwasResult:
    """BOLT-LMM-inf: infinitesimal mixed-model association with LOCO.

    Steps 1a-1b of Loh et al. (2015): variance components on all SNPs (not
    refitted per LOCO group, BOLT's default), V_{-g}^{-1} y for every LOCO
    group by batched conjugate gradients, the retrospective statistic for
    every SNP, and the calibration constant from exact prospective
    statistics at ``n_calibration`` random SNPs. No n x n matrix is formed.

    ``denominator="spectral"`` (mixmogam extension, not BOLT-LMM) replaces
    the proportional denominator by a structure-aware one; see
    :func:`_inf_core`.
    """
    rng = np.random.default_rng(random_state)
    st = _setup(y, gt, X, max_loco_groups, block)
    basis = _spectral_basis(st, random_state=random_state)
    fit = fit_variance_components(st, random_state=random_state, basis=basis)
    core = _inf_core(st, fit, n_calibration, rng, cg_tol, denominator, n_spectral, basis=basis)
    extra = {"method": "bolt-inf" if denominator == "constant" else "bolt-inf-spectral",
             "denominator": denominator, **_fit_extra(fit),
             "calibration": core["cal"]["c"], "calibration_cv": core["cal"]["cv"],
             "calibration_ratios": core["cal"]["ratios"],
             "spectral_k": core["cal"].get("k"),
             "cg_iterations": core["cg_iterations"]}
    return _result(st, gt, core["chi2"], core["beta_z"], core["se_z"], extra)


def bolt(y, gt, X=None, *, max_loco_groups: int = 25, n_calibration: int = 30,
         denominator: str = "constant", n_spectral="auto",
         grid=BOLT_GRID, n_folds: int = 5, min_cv_gain: float = 0.01,
         ld_window_bp: int = 1_000_000, vb_max_iter: int = 100, vb_tol: float = 1e-5,
         cg_tol: float = 1e-6, random_state: int = 0,
         block: int = 4096) -> GwasResult:
    """BOLT-LMM: Gaussian-mixture (non-infinitesimal) association with LOCO.

    Steps 2a-2b of Loh et al. (2015) on top of :func:`bolt_inf`. The prior
    is p N(0, (1 - f2) s2 / p) + (1 - p) N(0, f2 s2 / (1 - p)) with
    s2 = vg / M (total prior variance fixed by the REML fit); (f2, p) is
    chosen by ``n_folds``-fold cross-validated prediction over ``grid``. The
    mixture is used only when its CV prediction R^2 beats the
    infinitesimal model's by the relative margin ``min_cv_gain`` (BOLT's
    exact threshold is not given in the paper; this is mixmogam's choice).
    LOCO residuals come from the same variational iteration, and the
    statistic is calibrated by matching LD Score regression intercepts to
    BOLT-LMM-inf. Effects and standard errors come from that infinitesimal
    fit even when the mixture supplies the p-values, as in BOLT-LMM.

    ``denominator="spectral"`` (mixmogam extension, heuristic for the
    mixture statistic) multiplies each SNP's mixture statistic by the
    structure-aware factor kappa_j of the infinitesimal model before the
    calibration, and calibrates against the structure-aware BOLT-LMM-inf.
    """
    rng = np.random.default_rng(random_state)
    st = _setup(y, gt, X, max_loco_groups, block)
    basis = _spectral_basis(st, random_state=random_state)
    fit = fit_variance_components(st, random_state=random_state, basis=basis)
    inf = _inf_core(st, fit, n_calibration, rng, cg_tol, denominator, n_spectral, basis=basis)
    n, m, G = st.lg.n, st.lg.m, st.lg.n_groups
    s2 = fit.vg / m

    def prior_row(f2, p):
        return (p, (1.0 - f2) * s2 / p, f2 * s2 / (1.0 - p))

    # Step 2a: cross-validated choice of (f2, p)
    folds = rng.permutation(n) % n_folds
    eng = VBEngine(st.lg, folds=folds)
    P = len(grid) * n_folds
    col_fold = np.tile(np.arange(n_folds), len(grid))
    prior = np.array([prior_row(*grid[c // n_folds]) for c in range(P)])
    cv = eng.fit(np.repeat(st.y_p[:, None], P, axis=1), col_fold, np.full(P, -1),
                 PRIOR_MIXTURE, prior, np.full(P, fit.ve), max_iter=vb_max_iter, tol=vb_tol)
    pred = cv["prediction"]
    sse = np.zeros(len(grid))
    for c in range(P):
        test = folds == col_fold[c]
        sse[c // n_folds] += float(np.sum((st.y_p[test] - pred[test, c]) ** 2))
    mse = sse / n
    r2 = 1.0 - mse / float(np.mean(st.y_p**2))
    best = int(np.argmin(mse))
    i_inf = grid.index((0.5, 0.5)) if (0.5, 0.5) in grid else None
    r2_inf = r2[i_inf] if i_inf is not None else 0.0
    use_mixture = r2_inf > 0 and (r2[best] - r2_inf) > min_cv_gain * r2_inf
    extra = {"method": "bolt", "effect_method": "bolt-inf", "denominator": denominator, **_fit_extra(fit),
             "calibration_inf": inf["cal"]["c"], "calibration_cv": inf["cal"]["cv"],
             "calibration_ratios": inf["cal"]["ratios"], "spectral_k": inf["cal"].get("k"),
             "cv_grid": list(grid), "cv_r2": r2, "cv_best": grid[best],
             "cv_iterations": cv["iterations"], "cv_converged": cv["converged"],
             "use_mixture": bool(use_mixture)}
    # The CV fit's (m, folds x grid) effects and (n, folds x grid) residuals,
    # predictions and masks are not needed once its scores are taken.
    del cv, pred
    if not use_mixture:
        extra["note"] = "mixture did not beat the infinitesimal model in CV; BOLT-LMM-inf statistics"
        return _result(st, gt, inf["chi2"], inf["beta_z"], inf["se_z"], extra)

    # Step 2b: LOCO posterior-mean residuals and LD Score calibration
    loco = eng.fit(np.repeat(st.y_p[:, None], G, axis=1), np.full(G, -1), np.arange(G),
                    PRIOR_MIXTURE, np.tile(prior_row(*grid[best]), (G, 1)),
                    np.full(G, fit.ve), max_iter=vb_max_iter, tol=vb_tol)
    # Only the residuals of the LOCO fit are used; the engine's Gram cache
    # and the (m, groups) effects go before the LD scores need memory.
    loco_iterations, loco_converged = loco["iterations"], loco["converged"]
    resid = loco["resid"]
    del loco, eng
    rs = _retro_stats(st, resid)
    del resid
    if denominator == "spectral":
        rs["chi2"] = rs["chi2"] * inf["kappa"]
    ell = ld_scores(st.lg, window_bp=ld_window_bp)
    ok = np.isfinite(ell)
    ell_cv = float(np.std(ell[ok]) / max(np.mean(ell[ok]), 1e-12))
    if ell_cv >= 0.2:
        a_inf = ldsc_intercept(inf["chi2"], ell, n)
        a_raw = ldsc_intercept(rs["chi2"], ell, n)
        c_mix = a_raw["intercept"] / a_inf["intercept"]
        how = "ldsc"
        extra.update({"ldsc_intercept_inf": a_inf["intercept"],
                      "ldsc_intercept_mixture_raw": a_raw["intercept"]})
    else:
        # without LD-score variation the LDSC intercept is not identified
        # (intercept and slope are collinear); match the bulk instead
        fin = np.isfinite(rs["chi2"]) & np.isfinite(inf["chi2"])
        c_mix = float(np.median(rs["chi2"][fin]) / np.median(inf["chi2"][fin]))
        how = "median (LD scores uninformative)"
    chi2 = rs["chi2"] / c_mix
    # As in BOLT-LMM, report infinitesimal effect estimates alongside the
    # mixture test. In-sample mixture residuals give attenuated slopes;
    # calibrating their significance does not calibrate those slopes.
    extra.update({"calibration": c_mix, "calibration_method": how,
                  "ld_score_cv": ell_cv, "loco_iterations": loco_iterations,
                  "loco_converged": loco_converged})
    return _result(st, gt, chi2, inf["beta_z"], inf["se_z"], extra)


# ----------------------------------------------------------------------
# HRATT, after LDAK-KVIK
# ----------------------------------------------------------------------


def structure_test(st: _Setup, n_snps: int = 512, rng=None) -> dict:
    """LDAK-KVIK's test for strong structure (operation 1a).

    Mean squared correlation between ``n_snps`` evenly spread SNPs on
    different chromosomes. Without structure n r^2 averages ~1, so the
    statistic is ``excess = n_eff * mean(r^2) - 1`` (LDAK-KVIK: n rho^2 > 0.1);
    structure is strong when excess > 0.1 and its z-score (standard error
    from the spread over SNPs, conservatively treating each SNP as one
    unit) exceeds 3. In a row-scaled (weighted) model the correlations of
    scaled genotypes average sum v^2 / (sum v)^2 without structure, so the
    sample size is the Kish size of the working weights, less q.
    """
    rng = np.random.default_rng(0) if rng is None else rng
    m = st.lg.m
    k = min(n_snps, m)
    base = np.linspace(0, m - 1, k)
    jitter = rng.uniform(-0.5, 0.5, k) * (m / k)
    sel = np.unique(np.clip(np.round(base + jitter), 0, m - 1).astype(np.int64))
    Z = st.lg.rows(sel)
    norms = np.sqrt(np.einsum("ij,ij->i", Z, Z))
    norms[norms == 0] = 1.0
    R = (Z @ Z.T) / np.outer(norms, norms)
    chrom = np.asarray(st.lg.gt.chromosome)[sel]
    diff = chrom[:, None] != chrom[None, :]
    r2 = np.where(diff, R * R, np.nan)
    n_ref = st.n_eff if st.n_struct is None else st.n_struct
    per_snp = n_ref * np.nanmean(r2, axis=1)
    per_snp = per_snp[np.isfinite(per_snp)]
    excess = float(np.mean(per_snp) - 1.0)
    se = float(np.std(per_snp) / np.sqrt(max(per_snp.size, 1)))
    z = excess / se if se > 0 else np.inf
    return {"excess": excess, "z": float(z), "strong": bool(excess > 0.1 and z > 3.0),
            "n_snps": int(sel.size)}


def _he_alpha(st: _Setup, f: np.ndarray, alphas, n_probes: int, rng, *, fit_h2=False,
              design_weights=None) -> dict:
    """Pick the LDAK-Thin power by randomized Haseman-Elston regression.

    For each alpha, K_alpha = Z' W_alpha Z / sum(W_alpha) with
    W_alpha = [f (1 - f)]^(1 + alpha). HE regresses the off-diagonal
    products y_i y_j on K_ij, and at the least-squares h2 its residual sum
    of squares falls by (sum_{i != j} y_i y_j K_ij)^2 / sum_{i != j} K_ij^2,
    which is the fit score maximized over alpha. tr(K^2) comes from
    ``n_probes`` Rademacher probes shared by every alpha (common random
    numbers); one pass over the genotypes serves all alphas. This is the
    single-component version of LDAK-KVIK's partitioned randomized HE.

    ``design_weights`` (sampling weights D, mean one) make the variance fit
    design-consistent in the row-scaled model: the noise covariance is
    ve S D S and each diagonal moment is weighted by 1/D (off-diagonal
    moments already carry the pair weights w_i w_k). Alpha selection, an
    off-diagonal score, is unchanged.
    """
    if not isinstance(n_probes, (int, np.integer)) or n_probes < 2:
        raise ValueError("he_probes must be an integer of at least two")
    alphas = np.asarray(alphas, dtype=np.float64)
    if alphas.ndim != 1 or not alphas.size or not np.isfinite(alphas).all():
        raise ValueError("alphas must be a non-empty finite sequence")
    lg = st.lg
    n, A = lg.n, len(alphas)
    Wt = np.column_stack([(f * (1.0 - f)) ** (1.0 + a) for a in alphas])  # (m, A)
    ys = st.y_p / np.sqrt(np.sum(st.y_p**2) / st.n_eff)
    P = np.column_stack([ys, rng.choice(np.array([-1.0, 1.0]), size=(n, n_probes))])
    # Z P = z (I - QQ') P and sum Z' t = (I - QQ') sum z' t: the probes and
    # the accumulated products are projected once, never the genotypes.
    P32 = lg.project(P).astype(lg.dtype)
    q, c = lg.Q.shape[1], P.shape[1]
    KP = np.zeros((A, n, c))
    # The diagonal of Z' W Z from unprojected rows: with Z_j = z_j - c_j Q',
    # sum_j w_j Z_ij^2 = sum_j w_j z_ij^2 - 2 Q_i . X_i + Q_i S Q_i', where
    # X = z' W C rides along the product GEMM and S = C' W C is q x q.
    cross = np.zeros((A, n, q))
    S = np.zeros((A, q, q))
    diag = np.zeros((n, A))
    # Per-slice buffers, reused (grown to the largest slice): the products
    # with appended coefficients, their weighted copy, the n-row product and
    # the squared sample tiles.
    prod = np.empty((n, c + q), dtype=lg.dtype)
    T = None
    for idx, _, Z in lg.raw_slices():
        k = idx.size
        if T is None or T.shape[0] < k:
            T = np.empty((k, c + q), dtype=lg.dtype)
            Tw = np.empty_like(T)
            tile = max(1, (16 * 1024**2) // max(k * T.itemsize, 1))
            square = np.empty((k, min(tile, n)), dtype=lg.dtype)
        coef = lg._projection[idx]
        np.matmul(Z, P32, out=T[:k, :c])
        T[:k, c:] = coef
        wb = Wt[idx].astype(lg.dtype)
        for a in range(A):
            np.multiply(T[:k], wb[:, a, None], out=Tw[:k])
            np.matmul(Z.T, Tw[:k], out=prod)
            KP[a] += prod[:, :c]
            cross[a] += prod[:, c:]
            S[a] += (coef.T * Wt[idx, a]) @ coef
        # Square bounded sample tiles; each product reduces the whole slice.
        for start in range(0, n, tile):
            width = min(tile, n - start)
            Zs = np.multiply(Z[:, start : start + width], Z[:, start : start + width],
                             out=square[:k, :width])
            diag[start : start + width] += Zs.T @ wb
    tot = Wt.sum(axis=0)
    for a in range(A):
        KP[a] = lg.project(KP[a])
        if q:
            diag[:, a] += np.einsum("ik,ik->i", lg.Q @ S[a] - 2.0 * cross[a], lg.Q)
    KP /= tot[:, None, None]
    diag /= tot
    scores = np.empty(A)
    h2 = np.empty(A)
    for a in range(A):
        yky = float(ys @ KP[a][:, 0] - np.sum(ys * ys * diag[:, a]))
        k2 = float(np.mean(np.sum(KP[a][:, 1:] ** 2, axis=0)) - np.sum(diag[:, a] ** 2))
        if not np.isfinite([yky, k2]).all() or k2 <= 0:
            scores[a], h2[a] = -np.inf, np.nan
            continue
        scores[a] = yky * yky / k2
        h2[a] = yky / k2
    if not np.isfinite(scores).any():
        raise ValueError("HE alpha selection has no identifiable finite candidate; increase probes or use alpha_method='reml'")
    scores[~np.isfinite(scores)] = -np.inf
    best = int(np.argmax(scores))
    result = {"alpha": alphas[best], "weights": Wt[:, best], "scores": scores, "h2_he": h2}
    if fit_h2:
        # Reuse the products already needed to select alpha. Residual noise
        # has covariance I - QQ' (S D S with design weights D), so fit both K
        # and the noise term; reusing the off-diagonal slope above would
        # ignore the noise term's off-diagonals.
        from mixmogam._he import fit_he_moments, fit_projected_he

        probe_norm2 = np.sum(KP[best, :, 1:] ** 2, axis=0)
        yky_full = float(ys @ KP[best, :, 0])
        if design_weights is None:
            result["variance_fit"] = fit_projected_he(
                yky=yky_full, y2=float(ys @ ys), trace_k=float(diag[:, best].sum()),
                probe_norm2=probe_norm2, df=st.n_eff)
        else:
            D = np.asarray(design_weights, dtype=np.float64)
            d = diag[:, best]
            delta = 1.0 / D - 1.0
            Qb = lg.Q
            qn = np.einsum("ik,ik->i", Qb, Qb)
            B = (Qb * D[:, None]).T @ Qb  # Q' D Q
            N = D * (1.0 - 2.0 * qn) + np.einsum("ik,ik->i", Qb @ B, Qb)  # diag(S D S)
            y2 = ys * ys
            result["variance_fit"] = fit_he_moments(
                t_kk=float(np.mean(probe_norm2)) + float(np.sum(delta * d * d)),
                t_kn=float(np.sum(D * d) + np.sum(delta * d * N)),
                t_nn=float(np.sum(D * D) - 2.0 * np.sum(D * D * qn) + np.sum(B * B)
                           + np.sum(delta * N * N)),
                t_ky=float(yky_full + np.sum(delta * d * y2)),
                t_ny=float(np.sum(D * y2) + np.sum(delta * N * y2)),
                probe_norm2=probe_norm2)
    return result


def hratt(y, gt, X=None, *, trait: str = "quantitative", sample_weights=None,
          max_loco_groups: int = 25, alphas=HRATT_ALPHAS,
          alpha_method: str = "he", he_probes: int = 32,
          heritability_method: str = "auto",
          grid=HRATT_GRID, cv_fraction: float = 0.1, loco_folds: int = 5,
          n_calibration: int = 30, denominator: str = "constant", n_spectral="auto",
          structure_snps: int = 512, vb_max_iter: int = 100, vb_tol: float = 1e-5,
          cg_tol: float = 1e-6, spa_threshold: float = 2.0, spa_two_sided: str = "distance",
          random_state: int = 0,
          block: int = 4096, n_threads: int = 1) -> GwasResult:
    """HRATT, the Heritability-weighted Residual Association Two-step Test.

    Elastic-net LOCO polygenic scores as offsets, OLS score tests. Inspired
    by LDAK-KVIK (Hof & Speed 2025) and following its technical
    documentation, but not LDAK-KVIK: variance fitting, the variational
    sweep and the calibration solves are this package's own (below), and
    results differ from the reference program's. (1a) structure test;
    (1b) the LDAK-Thin power alpha is chosen among ``alphas`` by randomized
    Haseman-Elston regression (``alpha_method="he"``, all alphas in one
    pass; :func:`_he_alpha`) or by the REML likelihood of each weighted
    kinship (``"reml"``, one REML fit per alpha), and h2 is then estimated
    by REML at the chosen alpha, giving
    per-SNP heritabilities
    h2_j = w_j h2 / W with w_j = [f_j (1 - f_j)]^(1 + alpha) (no LD
    thinning: every SNP keeps weight w_j); (1d) the elastic-net prior
    p Laplace + (1 - p) N(0, F h2_j / (1 - p)) with (p, F) chosen on a
    ``cv_fraction`` hold-out over ``grid``; (1e) cross-fitted LOCO scores by
    variational Bayes: ``loco_folds`` (k, default 5) genome-wide fits, each
    on all but one fold of the samples, and a sample's score for group c is
    its fold's fit without group c's variants (REGENIE's level-0 design).
    Variance components and prior selection still use the whole trait.
    With strong structure each
    group's model is refitted without its variants instead (k x groups fits,
    ``extra["loco_refit"]``), since dropping one chromosome's share of a
    genome-wide fit leaves ancestry in the residual; (2) U_j = (z'(y - b_c
    P_c))^2 / (z'z s^2) scaled by lambda, with b_c the score's
    least-squares coefficient kept in [0, 1] (its out-of-fold calibration
    slope: an over-dispersed score is shrunk, none inflated). lambda = 1
    without strong structure; otherwise lambda = lambda' sum T / sum U over
    SNPs with lambda U < 6 and lambda' T < 6, where T are the ridge-score
    statistics and lambda' their GRAMMAR-Gamma calibration (the BOLT-LMM-inf
    constant). The phenotype is analysed on the unit-variance scale; effect
    sizes are returned in phenotype units. In-sample LOCO fits
    (``loco_folds=1``: LDAK-KVIK's design and this package's earlier
    default) absorb part of each effect into the score, which attenuates
    effect estimates (by about a third in simulations, with unchanged
    p-values there), and can give rare binary outcomes offsets that nearly
    separate the cases.

    ``heritability_method="he"`` reuses the selected alpha's existing HE
    products to fit the covariance ``vg K + ve (I - QQ')``, with nonnegative
    variance components. It reports the unconstrained estimates, boundary
    status and Monte Carlo precision in ``extra["he_variance"]``. This is
    single-component HE, not LDAK's partitioned estimator with large-effect
    exclusions. Alpha selection remains the same off-diagonal HE score.
    Nonidentified or numerically invalid HE fits raise an error; increase
    ``he_probes`` or, without sampling weights, use REML in that case. The
    default, ``"auto"``, is REML, or HE with sampling weights. The existing
    variational noise floor of 0.001 on
    the unit-variance scale still applies, and is reported for HE fits.
    LDAK optionally revises its HE estimate by Monte
    Carlo REML (by default for smaller samples or strong structure).

    ``n_threads`` controls optional parallel genotype preparation and Numba
    coordinate sweeps across independent candidate models and LOCO fits.
    Model sweeps retain their SNP order, and genotype preparation gives the
    same values at every thread count. The default is one; larger values require the ``fast`` extra and cannot exceed Numba's
    configured thread limit. Large fits also parallelize residual matrix
    products over sample rows, temporarily limiting detected BLAS pools to
    one thread. Apple Accelerate is not detected by threadpoolctl; set
    VECLIB_MAXIMUM_THREADS=1 before Python starts to avoid nested threading.

    ``denominator="spectral"`` (mixmogam extension, heuristic for the
    elastic-net statistic): with strong structure, U and T are multiplied by
    the structure-aware factors kappa_j of the weighted ridge model, T is
    then already calibrated (lambda' = 1), and lambda matches the two.

    ``sample_weights`` (finite and positive, e.g. inverse probabilities of
    participation; scaled to mean one) fit step 1 by weighted least squares
    in row-scaled form, with ``heritability_method="auto"`` choosing HE whose
    noise term and diagonal moments are made design-consistent (REML has no
    likelihood under sampling weights). Each variant is tested by the score
    of its weighted least-squares slope against the LOCO residual, z'a with
    a = w r (the weighted residual), and the variance that score has when
    genotypes are drawn given the covariates: zvar |a|^2 rho_j, with zvar
    the unweighted genotype variance after covariates and rho_j the factor
    from allele frequencies that vary with the covariates (as in SPAmix;
    :func:`_ancestry_ratio`; one when the covariates do not predict the
    genotype). Weights that vary with ancestry therefore need ancestry
    principal components among the covariates; HRATT warns when strong
    structure remains after them. At unit weights and rho_j = 1 this is the
    unweighted statistic. Above |z| = ``spa_threshold`` the saddlepoint
    approximation uses covariate-specific Binomial(2, p_ij) genotypes
    (conditional Hardy-Weinberg equilibrium); an intercept-only design
    uses empirical genotype frequencies. It gives both tails
    beyond +-|u| (``spa_two_sided="distance"``, as in SAIGE and REGENIE) or
    twice the tail beyond u (``"doubled"``, LDAK's default: half the level
    in each tail of a skewed null). lambda multiplies the statistic.
    Effects are weighted least-squares slopes. The Huber-White sandwich is
    not used: with heavy-tailed weighted residuals a few samples carry each
    score, and it is anti-conservative for low-frequency variants (design
    notes).

    ``trait="binary"`` (outcomes 0/1) fits step 1 as weighted linear
    regression of the working response (y - mu0) / (mu0 (1 - mu0)), with
    working weights mu0 (1 - mu0) of the covariate-only null logistic fit
    (times any sampling weights): one IRLS step, as LDAK-KVIK's binary step 1
    is "based on weighted linear regression". Step 2 refits the null
    logistic model per LOCO group with the logit-scale LOCO polygenic score
    as a covariate, its coefficient kept in [0, 1] as above
    (``extra["offset_slope"]``; scores that anti-predict, or whose fit
    fails, are dropped and counted in ``offsets_dropped``; in-sample scores
    are fixed offsets) and tests the logistic score
    z'a, a = w (y - mu) (w = 1 without weights), retrospectively as above,
    with the genotype saddlepoint approximation above |z| =
    ``spa_threshold`` (default 2; ``np.inf`` disables it). This variance
    does not trust the fitted probabilities, as the model-based variance
    sum mu (1 - mu) g~^2 does (design notes). Effects are one-step log odds
    ratios per counted allele, conditional on the score, not logistic MLEs.
    Their standard errors describe the null-score linearization and are
    independent of the tail rule. CSV/DataFrame exports label these
    ``beta_one_step`` and ``se_null_score``; posterior probabilities are
    unavailable for these approximations.
    ``denominator="spectral"`` is unavailable with weights or binary traits.
    """
    if trait not in ("quantitative", "binary"):
        raise ValueError("trait must be 'quantitative' or 'binary'")
    if heritability_method not in ("auto", "he", "reml"):
        raise ValueError("heritability_method must be 'auto', 'he' or 'reml'")
    weighted = sample_weights is not None
    if weighted and heritability_method == "reml":
        raise ValueError("REML has no likelihood under sampling weights; use heritability_method='he' or 'auto'")
    if weighted and alpha_method == "reml":
        raise ValueError("alpha_method='reml' is not available with sampling weights")
    if (weighted or trait == "binary") and denominator == "spectral":
        raise ValueError("denominator='spectral' is not available with sampling weights or binary traits")
    if heritability_method == "auto":
        heritability_method = "he" if weighted else "reml"
    if heritability_method == "he" and alpha_method != "he":
        raise ValueError("heritability_method='he' requires alpha_method='he'")
    if not (isinstance(spa_threshold, (int, float, np.floating)) and spa_threshold >= 0):
        raise ValueError("spa_threshold must be a nonnegative number (np.inf disables SPA)")
    if spa_two_sided not in ("distance", "doubled"):
        raise ValueError("spa_two_sided must be 'distance' or 'doubled'")
    if denominator not in ("constant", "spectral"):
        raise ValueError(f"unknown denominator {denominator!r}")
    if alpha_method not in ("he", "reml"):
        raise ValueError(f"unknown alpha_method {alpha_method!r}; use 'he' or 'reml'")
    if not isinstance(loco_folds, (int, np.integer)) or isinstance(loco_folds, bool) or loco_folds < 1:
        raise ValueError("loco_folds must be a positive integer (1: in-sample LOCO fits)")
    n_threads = _validate_n_threads(n_threads)
    if loco_folds > gt.n_samples:
        raise ValueError("loco_folds exceeds the number of samples")
    if weighted:
        sample_weights = _normalize_weights(sample_weights, gt.n_samples)
    rng = np.random.default_rng(random_state)
    st = _setup(y, gt, X, max_loco_groups, block, n_threads=n_threads,
                sample_weights=sample_weights, trait=trait)
    n = st.lg.n
    if linalg.norm(st.y_p) <= 10 * np.finfo(float).eps * linalg.norm(st.y):
        raise ValueError("phenotype has no residual variation after covariate adjustment")
    sy = float(np.sqrt(np.sum(st.y_p**2) / st.n_eff))
    if not np.isfinite(sy) or sy <= np.sqrt(np.finfo(float).tiny):
        raise ValueError("phenotype has no residual variation after covariate adjustment")
    ys = st.y_p / sy

    struct = structure_test(st, structure_snps, rng)
    if struct["strong"] and (weighted or trait == "binary"):
        warnings.warn("strong population structure remains after the covariates; the tests of "
                      "weighted and binary traits model allele frequencies by the covariates, "
                      "so include ancestry principal components among them",
                      RuntimeWarning, stacklevel=2)
    design = st.w if (st.w is not None and st.v is not None) else None
    vc = _hratt_variance(st, struct, alphas, alpha_method, heritability_method, he_probes, rng,
                         random_state, design_weights=design)
    h2, w = vc["h2"], vc["w"]
    h2j = w * h2 / float(np.sum(w))
    s2e = max(1.0 - h2, 1e-3)

    binary = trait == "binary"
    held = _holdout(n, cv_fraction, rng, strata=st.y_raw if binary else None)
    folds = None
    if loco_folds > 1:
        folds = _crossfit_folds(n, int(loco_folds), rng, strata=st.y_raw if binary else None)
    scores = _hratt_fit_scores(st, ys, h2j, s2e, grid, held, n_threads, vb_max_iter, vb_tol,
                               keep_prediction=binary, folds=folds, refit=struct["strong"])
    if st.v is not None:  # the held-out error per unit of working weight
        scores["mse"] = scores["mse"] * held.size / float(np.sum(st.v[held]))
    if binary and not (struct["strong"] and h2 > 0):
        rsU, U = None, None  # step 2 is logistic; no working-model statistic needed
        scores.pop("resid")
    else:
        rsU = _retro_stats(st, scores.pop("resid"), retrospective=weighted and not binary)
        U = rsU["chi2"]

    extra = {"method": "hratt", "denominator": denominator, "alpha": vc["alpha"], "h2": h2,
             "alpha_method": alpha_method, "alpha_scores": vc["alpha_scores"],
             "structure": struct,
             "cv_grid": list(grid), "cv_mse": scores["mse"], "cv_best": grid[scores["best"]],
             "cv_converged": scores["cv_converged"],
             "loco_folds": int(loco_folds), "loco_refit": bool(folds is not None and struct["strong"]),
             "loco_iterations": scores["loco_iterations"],
             "loco_converged": scores["loco_converged"]}
    if "offset_slope" in scores and not binary:
        extra["offset_slope"] = scores["offset_slope"]
    if heritability_method == "he":
        extra.update({"heritability_method": "he", "he_variance": vc["he"]["variance_fit"],
                      "vb_residual_variance": s2e,
                      "vb_noise_floor_applied": bool(1.0 - h2 < 1e-3)})
    if trait != "quantitative" or weighted:
        extra.update({"trait": trait, "heritability_method": heritability_method})
    if weighted:
        kish = float(np.sum(st.w) ** 2 / np.sum(st.w * st.w))
        extra.update({"kish_n": kish, "design_effect": n / kish})
    lam = 1.0
    if struct["strong"] and h2 > 0:
        # At the HE residual-zero boundary, use the same declared noise
        # floor for ridge calibration as for VB. A zero ridge would leave
        # the projected LOCO operator singular in the covariate directions.
        residual_variance = s2e if heritability_method == "he" else 1.0 - h2
        lam, U, lam_extra = _hratt_lambda(st, ys, U, h2, residual_variance, w, vc["basis"],
                                          n_calibration, denominator, n_spectral, cg_tol, rng)
        extra.update(lam_extra)
    extra["lambda"] = lam
    if binary:
        return _hratt_binary_step2(st, gt, scores.pop("prediction"), sy, lam, spa_threshold,
                                   n_threads, extra, crossfit=folds is not None,
                                   two_sided=spa_two_sided)
    p = None
    if weighted:
        # lambda, a structure correction of the model-based statistic,
        # multiplies the retrospective one; saddlepoint tails above the
        # threshold.
        p, chi2, n_spa = _genotype_tails(st, rsU["A"], rsU["num"], rsU["V"],
                                         lam * rsU["chi2_retro"], lam, spa_threshold, n_threads,
                                         spa_two_sided)
        extra.update({"spa_threshold": float(spa_threshold), "spa_two_sided": spa_two_sided,
                      "n_spa": n_spa})
    else:
        chi2 = lam * U
    with np.errstate(divide="ignore", invalid="ignore"):
        beta_z = lam * rsU["num"] / rsU["zz"] * sy
        se_z = (np.sqrt(lam * rsU["V"]) / rsU["zz"] * sy if weighted else
                np.abs(beta_z) / np.sqrt(chi2))
    return _result(st, gt, chi2, beta_z, se_z, extra, p=p)


def _hratt_variance(st: _Setup, struct: dict, alphas, alpha_method: str, heritability_method: str,
                    he_probes: int, rng, random_state, design_weights=None) -> dict:
    """(1b) The LDAK-Thin power alpha and the variance components at it.

    Returns the chosen ``alpha``, its variant weights ``w``, ``h2``, the
    spectral ``basis`` shared with the ridge solves under strong structure
    (or None), the per-alpha ``alpha_scores`` and the HE products ``he``
    (or None).
    """
    f = np.clip(st.lg.mean / 2.0, 1e-6, 1 - 1e-6)
    # Under strong structure the ridge solves need a preconditioner: REML
    # computes its deflation basis explicitly so that they can share it.
    basis, fit, he = None, None, None
    if alpha_method == "he":
        he = _he_alpha(st, f, list(alphas), he_probes, rng,
                       fit_h2=heritability_method == "he", design_weights=design_weights)
        alpha, w = he["alpha"], he["weights"]
        alpha_scores = he["scores"]
        if heritability_method == "he":
            h2 = he["variance_fit"]["h2"]
        else:
            basis = _spectral_basis(st, w, random_state) if struct["strong"] else None
            fit = fit_variance_components(st, weights=w, random_state=random_state, basis=basis)
    elif alpha_method == "reml":
        fits = []
        for a in alphas:
            w_a = (f * (1.0 - f)) ** (1.0 + a)
            basis_a = _spectral_basis(st, w_a, random_state) if struct["strong"] else None
            fit_a = fit_variance_components(st, weights=w_a, random_state=random_state,
                                            basis=basis_a)
            fits.append((fit_a.ll, a, w_a, fit_a, basis_a))
        _, alpha, w, fit, basis = max(fits, key=lambda t: t[0])
        alpha_scores = np.array([t[0] for t in fits])
    else:
        raise ValueError(f"unknown alpha_method {alpha_method!r}; use 'he' or 'reml'")
    if heritability_method == "reml":
        h2 = float(fit.pseudo_heritability)
    return {"alpha": alpha, "w": w, "h2": h2, "basis": basis,
            "alpha_scores": alpha_scores, "he": he}


def _holdout(n: int, fraction: float, rng, strata=None) -> np.ndarray:
    """(1d) Indices of the cross-validation hold-out.

    Without ``strata`` this is one draw of ``round(fraction * n)`` samples.
    With ``strata`` (one label per sample), each stratum contributes its
    rounded share, at least one sample, so rare classes are represented.
    """
    if strata is None:
        return rng.choice(n, size=max(1, int(round(fraction * n))), replace=False)
    strata = np.asarray(strata)
    held = [rng.choice(members, size=max(1, int(round(fraction * members.size))), replace=False)
            for members in (np.flatnonzero(strata == label) for label in np.unique(strata))]
    return np.sort(np.concatenate(held))


def _crossfit_folds(n: int, k: int, rng, strata=None) -> np.ndarray:
    """(1e) Fold labels 0..k-1 of the cross-fitted LOCO fits: a random split
    into near-equal folds, made within each stratum (one label per sample)
    when ``strata`` is given, so that rare cases spread over the folds.
    Strata are taken smallest first (then by first member), not by label, so
    that swapping case and control labels leaves the folds unchanged."""
    folds = np.empty(n, dtype=np.int64)
    members = [np.arange(n)] if strata is None else [np.flatnonzero(strata == v) for v in np.unique(strata)]
    members.sort(key=lambda idx: (idx.size, idx[0]))
    start = 0
    for idx in members:
        folds[rng.permutation(idx)] = (start + np.arange(idx.size)) % k
        start = (start + idx.size) % k
    return folds


def _hratt_fit_scores(st: _Setup, ys: np.ndarray, h2j: np.ndarray, s2e: float, grid, held,
                      n_threads: int, vb_max_iter: int, vb_tol: float,
                      keep_prediction: bool = False, folds=None, refit: bool = False) -> dict:
    """(1d)-(1e) Hold-out choice of the elastic-net prior, then the LOCO fits.

    Returns the LOCO residuals (``resid``), the chosen grid index ``best``
    with the held-out ``mse``, and convergence records; with
    ``keep_prediction`` also the LOCO predictions. The cross-validation
    effects are released before the LOCO fit, and the LOCO effects and the
    engine's Gram cache before returning.

    With ``folds`` (labels 0..k-1, :func:`_crossfit_folds`) the LOCO scores
    are cross-fitted: as in REGENIE's level 0, k genome-wide fits,
    fold f's on the other folds (conditional on the selected model),
    group g's score of a sample in fold f being fold f's prediction without
    group g's own variants (:func:`_out_of_fold_loco`); or, with ``refit``
    (strong structure), k x groups fits, each group's without its variants.
    A genome-wide fit spreads the ancestry signal over every chromosome, so
    under strong structure dropping one group's share leaves ancestry in the
    residual (null lambda_GC 1.5-4.4 at Fst 0.15 without principal
    components, against 0.9-1.1 refitted). Each group's residual is then
    the projected phenotype less its score times the score's least-squares
    coefficient kept in [0, 1] (``offset_slope``: its out-of-fold
    calibration slope, which shrinks an over-dispersed score but never
    undoes the prior's shrinkage). Without polygenic variance (all ``h2j``
    zero) there are no scores.
    """
    n, G = st.lg.n, st.lg.n_groups

    def prior_row(p, F):
        lam = np.sqrt(2.0 * p / (1.0 - F)) if p > 0 else 0.0
        v = F / (1.0 - p) if p < 1 else 0.0
        return (p, lam, v)

    cv_folds = np.ones(n, dtype=np.int64)
    cv_folds[held] = 0  # fold 0 is held out; fold 1 is never held out
    eng = VBEngine(st.lg, folds=cv_folds, n_threads=n_threads)
    Pn = len(grid)
    prior = np.array([prior_row(*g_) for g_ in grid])
    cv = eng.fit(np.repeat(ys[:, None], Pn, axis=1), np.zeros(Pn, dtype=np.int64),
                 np.full(Pn, -1), PRIOR_ENET, prior, np.full(Pn, s2e), snp_scale=h2j,
                 max_iter=vb_max_iter, tol=vb_tol)
    mse = np.mean((ys[held, None] - cv["prediction"][held]) ** 2, axis=0)
    best = int(np.argmin(mse))
    cv_converged = cv["converged"]
    del cv  # the (m, grid) effects and (n, grid) arrays are not needed again

    # (1e) LOCO elastic-net scores
    if folds is None:
        loco = eng.fit(np.repeat(ys[:, None], G, axis=1), np.full(G, -1), np.arange(G),
                       PRIOR_ENET, np.tile(prior_row(*grid[best]), (G, 1)), np.full(G, s2e),
                       snp_scale=h2j, max_iter=vb_max_iter, tol=vb_tol)
        resid, prediction, slope = loco["resid"], loco["prediction"], None
    elif not np.any(h2j > 0):  # no polygenic variance, no scores
        eng, loco = None, {"iterations": 0, "converged": True}
        prediction, slope = np.zeros((n, G)), np.zeros(G)
        resid = np.repeat(ys[:, None], G, axis=1)
    else:
        del eng  # the cross-fit has its own folds and Gram cache
        k = int(folds.max()) + 1
        eng = VBEngine(st.lg, folds=folds, n_threads=n_threads)
        col_fold = np.repeat(np.arange(k), G) if refit else np.arange(k)
        col_group = np.tile(np.arange(G), k) if refit else np.full(k, -1)
        P = col_fold.size
        loco = eng.fit(np.repeat(ys[:, None], P, axis=1), col_fold, col_group,
                       PRIOR_ENET, np.tile(prior_row(*grid[best]), (P, 1)), np.full(P, s2e),
                       snp_scale=h2j, max_iter=vb_max_iter, tol=vb_tol)
        if refit:
            prediction = np.empty((n, G))
            for c in range(P):
                rows = folds == col_fold[c]
                prediction[rows, col_group[c]] = loco["prediction"][rows, c]
        else:
            prediction = _out_of_fold_loco(st.lg, loco["beta"], folds)
        Pp = st.lg.project(prediction)  # ys is projected already
        norm2 = np.einsum("ij,ij->j", Pp, Pp)
        slope = np.zeros(G)
        np.divide(Pp.T @ ys, norm2, out=slope, where=norm2 > 0)
        # Shrink an over-dispersed score, never inflate one: a slope above
        # one would undo the prior's shrinkage, and one below zero is noise.
        np.clip(slope, 0.0, 1.0, out=slope)
        resid = ys[:, None] - Pp * slope
    # Only the LOCO residuals (and predictions) are used; the engine's Gram
    # cache and the (m, groups) effects go before the calibration solves.
    out = {"resid": resid, "best": best, "mse": mse, "cv_converged": cv_converged,
           "loco_iterations": loco["iterations"], "loco_converged": loco["converged"]}
    if slope is not None:
        out["offset_slope"] = slope
    if keep_prediction:
        out["prediction"] = prediction
    del loco, eng
    return out


def _out_of_fold_loco(lg: LocoGenotypes, beta: np.ndarray, folds: np.ndarray) -> np.ndarray:
    """(1e) Out-of-fold LOCO scores (n, groups) from k genome-wide fits.

    ``beta`` (m, k) holds fold f's effects in column f. A sample in fold f
    gets, for group g, fold f's fitted value Z'beta_f without group g's
    variants: its total less the group's own share, with the covariate
    projection Q (c'beta) split the same way. One pass over the genotypes.
    """
    n, G, k = lg.n, lg.n_groups, beta.shape[1]
    q = lg.Q.shape[1]
    own = np.zeros((n, G))  # each sample's group shares of z'beta under its fold's fit
    coef = np.zeros((k, q, G))  # c'beta per fold and group
    pick = np.arange(n)
    for idx, g, Z in lg.raw_slices():
        bb = beta[idx]
        rows = max(1, (16 * 1024**2) // (8 * idx.size))
        for start in range(0, n, rows):
            stop = min(start + rows, n)
            part = Z[:, start:stop].T.astype(np.float64) @ bb
            own[start:stop, g] += part[pick[: stop - start], folds[start:stop]]
        if q:
            coef[:, :, g] += (lg._projection[idx].T @ bb).T
    out = own.sum(axis=1, keepdims=True) - own
    if q:
        for f in range(k):
            rows = folds == f
            out[rows] -= lg.Q[rows] @ (coef[f].sum(axis=1, keepdims=True) - coef[f])
    return out


def _hratt_lambda(st: _Setup, ys: np.ndarray, U: np.ndarray, h2: float, residual_variance: float,
                  w: np.ndarray, basis, n_calibration: int, denominator: str, n_spectral, cg_tol,
                  rng) -> tuple:
    """(2) The calibration lambda under strong structure.

    Ridge-prior LOCO scores: their OLS statistic is the BOLT-LMM-inf
    retrospective statistic, lambda' its prospective calibration. Returns
    ``(lambda, U, extra)``; with the spectral denominator ``U`` comes back
    multiplied by the structure-aware factors kappa.
    """
    n, G = st.lg.n, st.lg.n_groups
    delta = residual_variance / h2
    # The ridge residuals and the calibration SNPs share one batched CG
    # run, preconditioned by the REML deflation basis when there is one.
    if basis is not None:
        pre = _preconditioner(st, basis)
    else:
        op = _KOp(st.lg, w)
        pre = SpectralPreconditioner(op.matmul, n, op.trace, k=min(64, n - 2))
    cand = _calibration_draw(st, n_calibration, rng)
    Wr, drawn, _ = _solve_with_calibration(st, delta, np.repeat(ys[:, None], G, axis=1),
                                           np.arange(G), cand, weights=w, tol=cg_tol, pre=pre)
    rsT = _retro_stats(st, Wr)
    cal = _calibrate_inf(st, h2, delta, rsT, n_calibration, rng, drawn, weights=w,
                         pre=pre, tol=cg_tol)
    lam_p = 1.0 / cal["c"]
    T = rsT["chi2"]
    extra = {}
    if denominator == "spectral":
        sp = _spectral_denominator(st, h2, delta, cal, n_spectral, weights=w)
        kappa = (cal["c"] * rsT["zz"] / st.n_eff * rsT["wnorm2"]) / sp["D"]
        T = T * kappa / cal["c"]  # = num_T^2 / D_spectral: calibrated
        U = U * kappa
        lam_p = 1.0
        cal = {**cal, "cv": sp["cv"], "ratios": sp["ratios"]}
        extra["spectral_k"] = sp["k"]
    lam = lam_p
    for _ in range(20):
        S = (lam * U < 6) & (lam_p * T < 6) & np.isfinite(U) & np.isfinite(T)
        new = lam_p * float(np.sum(T[S])) / float(np.sum(U[S]))
        if abs(new - lam) < 1e-10:
            break
        lam = new
    extra.update({"lambda_prime": lam_p, "calibration_cv": cal["cv"],
                  "calibration_ratios": cal["ratios"]})
    return lam, U, extra


def _hratt_binary_step2(st: _Setup, gt, prediction, sy, lam, spa_threshold, n_threads,
                        extra, crossfit: bool = False, two_sided: str = "distance") -> GwasResult:
    """Logistic score tests of a binary trait with the LOCO polygenic scores
    as offsets (step 2 of binary HRATT).

    Per LOCO group the null logistic model is refitted with the group's
    score (sampling-weighted when there are weights): as a covariate with a
    fitted coefficient when the scores are cross-fitted (``crossfit``), as
    a fixed offset otherwise. The score z_j'a_g of the residual a_g = w (y -
    mu_g) is tested against its variance when genotypes are drawn given the
    covariates, zvar_j |a_g|^2 rho_j, as in the weighted quantitative test;
    lambda multiplies the statistic, and above |z| = ``spa_threshold`` the
    genotype saddlepoint approximation (``two_sided`` tails) replaces the
    normal tail. Unlike the model-based variance sum mu (1 - mu) g~^2, this
    does not assume the fitted probabilities: an in-sample LOCO offset that
    absorbs part of a rare outcome shrinks score and variance together (a
    pilot with 50 cases in 5,000 gave lambda_GC 0.55 with the model-based
    variance, 1.04 with this one). Effects are one-step log odds ratios
    U / J per counted allele, with null-score uncertainty sqrt(V/lambda)/J.
    """
    from mixmogam._binary import group_null_fits, loco_offsets, score_pass

    offsets = loco_offsets(prediction, st.s, sy)
    del prediction
    y = st.y_raw
    try:
        fits = group_null_fits(y, st.X_raw, offsets, weights=st.w, free_offset=crossfit)
    except ValueError as exc:
        if crossfit:
            raise
        raise ValueError(f"{exc}; with in-sample LOCO scores (loco_folds=1) the offset can nearly "
                         "separate rare cases; use loco_folds >= 2") from None
    sp = score_pass(st.lg, y, st.X_raw, fits, weights=st.w, s=st.s, Q=st.Q_raw)
    U, J, A = sp["U"], sp["J"], sp["A"]
    w = np.ones(y.size) if st.w is None else st.w
    V = sp["zvar"] * np.einsum("ij,ij->j", A, A)[st.lg.groups] * sp["rho"]
    valid = (J > 0) & (V > 0) & (st.lg.zz > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2_normal = np.where(valid, lam * U * U / V, np.nan)
    p, chi2, n_spa = _genotype_tails(st, A, U, V, chi2_normal, lam, spa_threshold, n_threads,
                                     two_sided)
    with np.errstate(divide="ignore", invalid="ignore"):
        beta_z = np.where(valid, U / J, np.nan)
        # Uncertainty of the null-score linearization, independent of the
        # choice of tail approximation. This is not a fitted logistic MLE.
        se_z = np.where(valid, np.sqrt(V / lam) / J, np.nan)
    off = offsets - offsets.mean(axis=0)
    extra.update({"n_cases": int(y.sum()),
                  "n_controls": int(y.size - y.sum()),
                  "prevalence": float(np.sum(w * y) / np.sum(w)),
                  "null_converged": bool(all(f["converged"] for f in fits)),
                  "mu0_clipped": st.mu0_clipped, "offset_sd": float(off.std()),
                  "spa_threshold": float(spa_threshold), "spa_two_sided": two_sided,
                  "n_spa": n_spa, "effect_method": "one-step-logistic"})
    if crossfit:
        extra.update({"offset_slope": np.array([f["slope"] for f in fits]),
                      "offsets_dropped": int(sum(f["dropped"] for f in fits))})
    return _result(st, gt, chi2, beta_z, se_z, extra, p=p)
