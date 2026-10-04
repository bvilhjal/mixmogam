"""Two-step mixed-model association: BOLT-LMM and LDAK-KVIK.

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
- LDAK-KVIK (Hof & Speed 2025, Nat Genet): w = y minus an elastic-net
  LOCO polygenic score under the LDAK-Thin heritability model; lambda = 1
  unless a test of inter-chromosome correlation finds strong structure, in
  which case lambda is matched to the GRAMMAR-Gamma-calibrated ridge
  statistic on weakly associated SNPs.

The single calibration constant rests on the assumption that the
prospective denominator z_j' V^{-1} z_j is proportional to z_j' z_j across
SNPs. BOLT-LMM's authors flag it as possibly invalid outside human data.
Every result here therefore reports the spread of the ratio over the
calibration SNPs (``calibration_cv``): with strong structure it is large,
and a single constant cannot calibrate every SNP.

Deviations from the reference implementations are documented on each
function. Variance components default to this package's stochastic Lanczos
REML; KVIK also offers projected single-component randomized HE, distinct
from the reference program's partitioned HE and optional Monte Carlo REML.
"""

from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np
from scipy import linalg, stats

from mixmogam._cg import SpectralPreconditioner, batched_pcg
from mixmogam._loco import LocoGenotypes, loco_groups
from mixmogam._slq import randomized_eigh_op
from mixmogam._vb import PRIOR_ENET, PRIOR_MIXTURE, VBEngine, _validate_n_threads
from mixmogam.ldscore import ld_scores, ldsc_intercept
from mixmogam.lmm import LMM, _design_matrix
from mixmogam.results import GwasResult

__all__ = ["bolt_inf", "bolt", "kvik", "structure_test", "BOLT_GRID", "KVIK_GRID"]

# BOLT-LMM: f2 = spike share of the prior variance, p = slab probability
BOLT_GRID = [(f2, p) for f2 in (0.5, 0.3, 0.1) for p in (0.5, 0.2, 0.1, 0.05, 0.02, 0.01)]
# LDAK-KVIK elastic net: (p, F) = (lasso share, ridge variance share)
KVIK_GRID = [(0.5, 0.5), (0.5, 0.3), (0.5, 0.1), (0.1, 0.5), (0.1, 0.3), (0.1, 0.1),
             (0.01, 0.5), (0.01, 0.3), (0.01, 0.1), (0.0, 1.0)]
KVIK_ALPHAS = (-1.0, -0.75, -0.5, -0.25, 0.0)


# ----------------------------------------------------------------------
# Shared machinery
# ----------------------------------------------------------------------


@dataclass
class _Setup:
    y: np.ndarray  # phenotype
    y_p: np.ndarray  # covariate-projected phenotype
    X: np.ndarray  # design incl. intercept
    lg: LocoGenotypes
    labels: list
    n_eff: int  # n - q


def _setup(y, gt, X, max_loco_groups, block, cache_bytes, n_threads=1) -> _Setup:
    y = np.asarray(y, dtype=np.float64).ravel()
    if not np.isfinite(y).all():
        raise ValueError("y contains NaN or infinite values")
    if y.size != gt.n_samples:
        raise ValueError("y and genotypes disagree on the sample count")
    Xd = _design_matrix(X, y.size, add_intercept=True)
    Q = linalg.qr(Xd, mode="economic")[0]
    groups, labels = loco_groups(gt.chromosome, max_loco_groups)
    if len(labels) < 2:
        raise ValueError("LOCO needs variants on at least two chromosomes/groups")
    lg = LocoGenotypes(gt, groups, Q, block=block, cache_bytes=cache_bytes,
                       n_threads=n_threads)
    y_p = y - Q @ (Q.T @ y)
    return _Setup(y=y, y_p=y_p, X=Xd, lg=lg, labels=labels, n_eff=y.size - Xd.shape[1])


class _KOp:
    """Kinship operator adapter for :class:`LMM` (K = Z' W Z / sum W)."""

    def __init__(self, lg: LocoGenotypes, weights=None):
        self.lg = lg
        self.n = lg.n
        self.weights = weights
        self._trace = None

    @property
    def trace(self):
        # SLQ only needs products. Compute this separate genotype pass when
        # a spectral preconditioner actually requests the trace.
        if self._trace is None and self.weights is None:
            self._trace = self.lg.trace
        elif self._trace is None:
            tot = 0.0
            for idx, _, Z in self.lg.blocks(reuse=True):
                tot += float(np.einsum("ij,ij,i->", Z, Z, self.weights[idx].astype(Z.dtype),
                                       dtype=np.float64))
            self._trace = tot / float(np.sum(self.weights))
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


def _retro_stats(st: _Setup, W: np.ndarray) -> dict:
    """Uncalibrated retrospective chi2 of every variant against its group's
    column of W: n_eff (z'w)^2 / (z'z |w|^2)."""
    m = st.lg.m
    num = np.zeros(m)
    zz = np.zeros(m)
    wn = np.einsum("ij,ij->j", W, W)
    # Widen at most 16 MiB of variant rows, retaining each full sample
    # reduction. The minimum one-row tile and native BLAS scratch are separate.
    tile = max(1, (16 * 1024**2) // max(8 * st.lg.n, 1))
    for idx, g, Z in st.lg.blocks(reuse=True):
        for start in range(0, idx.size, tile):
            take = idx[start : start + tile]
            Zs = Z[start : start + tile]
            num[take] = Zs.astype(np.float64, copy=False) @ W[:, g]
            zz[take] = np.einsum("ij,ij->i", Zs, Zs, dtype=np.float64)
    wg = wn[st.lg.groups]
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = st.n_eff * num * num / (zz * wg)
    chi2[~(zz > 0)] = np.nan
    return {"chi2": chi2, "num": num, "zz": zz, "wnorm2": wg}


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
    sq = np.zeros(G)
    wsum = np.zeros(G)
    for idx, g, Z in lg.blocks(reuse=True):
        w = np.ones(idx.size) if weights is None else weights[idx]
        sq[g] += float(np.einsum("ij,ij,i->", Z.astype(np.float64), Z.astype(np.float64), w))
        wsum[g] += float(w.sum())
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
    for idx, g, Z in st.lg.blocks(reuse=True):
        vals, U, lam_bar = bases[g] if len(bases) > 1 else bases[0]
        Z64 = Z.astype(np.float64)
        P = Z64 @ U
        top = np.einsum("ij,ij->i", P, P)
        zz = np.einsum("ij,ij->i", Z64, Z64)
        out[idx] = (P * P) @ (1.0 / (vals + delta)) + (zz - top) / (lam_bar + delta)
    return out


def _result(st: _Setup, gt, chi2, beta_z, se_z, extra: dict) -> GwasResult:
    if any(extra.get(key) is False for key in ("cv_converged", "loco_converged")):
        warnings.warn("variational fit did not converge; inspect cv_converged/loco_converged before using results",
                      RuntimeWarning, stacklevel=3)
    sd = st.lg.sd
    with np.errstate(divide="ignore", invalid="ignore"):
        beta = np.where(sd > 0, beta_z / sd, np.nan)
        se = np.where(sd > 0, se_z / sd, np.nan)
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
             cg_tol: float = 1e-6, random_state: int = 0, block: int = 4096,
             cache_bytes: float = 4e9) -> GwasResult:
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
    st = _setup(y, gt, X, max_loco_groups, block, cache_bytes)
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
         cg_tol: float = 1e-6, random_state: int = 0, block: int = 4096,
         cache_bytes: float = 4e9) -> GwasResult:
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
    BOLT-LMM-inf.

    ``denominator="spectral"`` (mixmogam extension, heuristic for the
    mixture statistic) multiplies each SNP's mixture statistic by the
    structure-aware factor kappa_j of the infinitesimal model before the
    calibration, and calibrates against the structure-aware BOLT-LMM-inf.
    """
    rng = np.random.default_rng(random_state)
    st = _setup(y, gt, X, max_loco_groups, block, cache_bytes)
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
    extra = {"method": "bolt", "denominator": denominator, **_fit_extra(fit),
             "calibration_inf": inf["cal"]["c"], "calibration_cv": inf["cal"]["cv"],
             "calibration_ratios": inf["cal"]["ratios"], "spectral_k": inf["cal"].get("k"),
             "cv_grid": list(grid), "cv_r2": r2, "cv_best": grid[best],
             "cv_iterations": cv["iterations"], "cv_converged": cv["converged"],
             "use_mixture": bool(use_mixture)}
    if not use_mixture:
        extra["note"] = "mixture did not beat the infinitesimal model in CV; BOLT-LMM-inf statistics"
        return _result(st, gt, inf["chi2"], inf["beta_z"], inf["se_z"], extra)

    # Step 2b: LOCO posterior-mean residuals and LD Score calibration
    loco = eng.fit(np.repeat(st.y_p[:, None], G, axis=1), np.full(G, -1), np.arange(G),
                    PRIOR_MIXTURE, np.tile(prior_row(*grid[best]), (G, 1)),
                    np.full(G, fit.ve), max_iter=vb_max_iter, tol=vb_tol)
    del loco["prediction"]  # Only CV needs the fitted values after residual reconstruction.
    rs = _retro_stats(st, loco["resid"])
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
    with np.errstate(divide="ignore", invalid="ignore"):
        beta_z = rs["num"] / rs["zz"]
        se_z = np.abs(beta_z) / np.sqrt(chi2)
    extra.update({"calibration": c_mix, "calibration_method": how,
                  "ld_score_cv": ell_cv, "loco_iterations": loco["iterations"],
                  "loco_converged": loco["converged"]})
    return _result(st, gt, chi2, beta_z, se_z, extra)


# ----------------------------------------------------------------------
# LDAK-KVIK
# ----------------------------------------------------------------------


def structure_test(st: _Setup, n_snps: int = 512, rng=None) -> dict:
    """LDAK-KVIK's test for strong structure (operation 1a).

    Mean squared correlation between ``n_snps`` evenly spread SNPs on
    different chromosomes. Without structure n r^2 averages ~1, so the
    statistic is ``excess = n_eff * mean(r^2) - 1`` (KVIK: n rho^2 > 0.1);
    structure is strong when excess > 0.1 and its z-score (standard error
    from the spread over SNPs, conservatively treating each SNP as one
    unit) exceeds 3.
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
    per_snp = st.n_eff * np.nanmean(r2, axis=1)
    per_snp = per_snp[np.isfinite(per_snp)]
    excess = float(np.mean(per_snp) - 1.0)
    se = float(np.std(per_snp) / np.sqrt(max(per_snp.size, 1)))
    z = excess / se if se > 0 else np.inf
    return {"excess": excess, "z": float(z), "strong": bool(excess > 0.1 and z > 3.0),
            "n_snps": int(sel.size)}


def _he_alpha(st: _Setup, f: np.ndarray, alphas, n_probes: int, rng, *, fit_h2=False) -> dict:
    """Pick the LDAK-Thin power by randomized Haseman-Elston regression.

    For each alpha, K_alpha = Z' W_alpha Z / sum(W_alpha) with
    W_alpha = [f (1 - f)]^(1 + alpha). HE regresses the off-diagonal
    products y_i y_j on K_ij, and at the least-squares h2 its residual sum
    of squares falls by (sum_{i != j} y_i y_j K_ij)^2 / sum_{i != j} K_ij^2,
    which is the fit score maximized over alpha. tr(K^2) comes from
    ``n_probes`` Rademacher probes shared by every alpha (common random
    numbers); one pass over the genotypes serves all alphas. This is the
    single-component version of LDAK-KVIK's partitioned randomized HE.
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
    P32 = P.astype(lg.dtype)
    KP = np.zeros((A, n, P.shape[1]))
    diag = np.zeros((n, A))
    for idx, _, Z in lg.blocks(reuse=True):
        T = Z @ P32
        wb = Wt[idx].astype(lg.dtype)
        for a in range(A):
            KP[a] += (Z.T @ (T * wb[:, a, None])).astype(np.float64)
        # Square sample tiles instead of the whole genotype block. Each
        # product still reduces over all variants in the block, in the same
        # storage precision. Budget the square plus its narrow product/cast.
        # One sample is the minimum tile; native BLAS scratch is separate.
        per_sample = Z.shape[0] * Z.dtype.itemsize + A * (Z.dtype.itemsize + 8)
        tile = max(1, (16 * 1024**2) // max(per_sample, 1))
        for start in range(0, n, tile):
            Zs = Z[:, start : start + tile]
            diag[start : start + tile] += ((Zs * Zs).T @ wb).astype(np.float64, copy=False)
    tot = Wt.sum(axis=0)
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
        from mixmogam._he import fit_projected_he

        # Reuse the products already needed to select alpha. Residual noise
        # has covariance I - QQ', so fit both K and that projector; simply
        # reusing the off-diagonal slope above would ignore its off-diagonals.
        result["variance_fit"] = fit_projected_he(
            yky=float(ys @ KP[best, :, 0]), y2=float(ys @ ys),
            trace_k=float(diag[:, best].sum()),
            probe_norm2=np.sum(KP[best, :, 1:] ** 2, axis=0), df=st.n_eff)
    return result


def kvik(y, gt, X=None, *, max_loco_groups: int = 25, alphas=KVIK_ALPHAS,
         alpha_method: str = "he", he_probes: int = 32,
         heritability_method: str = "reml",
         grid=KVIK_GRID, cv_fraction: float = 0.1, n_calibration: int = 30,
         denominator: str = "constant", n_spectral="auto",
         structure_snps: int = 512, vb_max_iter: int = 100, vb_tol: float = 1e-5,
         cg_tol: float = 1e-6, random_state: int = 0, block: int = 4096,
         cache_bytes: float = 4e9, n_threads: int = 1) -> GwasResult:
    """LDAK-KVIK: elastic-net LOCO polygenic scores as offsets, OLS score tests.

    Following the LDAK-KVIK technical documentation: (1a) structure test;
    (1b) the LDAK-Thin power alpha is chosen among ``alphas`` by randomized
    Haseman-Elston regression (``alpha_method="he"``, all alphas in one
    pass; :func:`_he_alpha`) or by the REML likelihood of each weighted
    kinship (``"reml"``, one REML fit per alpha), and h2 is then estimated
    by REML at the chosen alpha, giving
    per-SNP heritabilities
    h2_j = w_j h2 / W with w_j = [f_j (1 - f_j)]^(1 + alpha) (no LD
    thinning: every SNP keeps weight w_j); (1d) the elastic-net prior
    p Laplace + (1 - p) N(0, F h2_j / (1 - p)) with (p, F) chosen on a
    ``cv_fraction`` hold-out over ``grid``; (1e) LOCO scores by variational
    Bayes; (2) U_j = (z'(y - P_c))^2 / (z'z s^2) scaled by lambda. lambda = 1
    without strong structure; otherwise lambda = lambda' sum T / sum U over
    SNPs with lambda U < 6 and lambda' T < 6, where T are the ridge-score
    statistics and lambda' their GRAMMAR-Gamma calibration (the BOLT-LMM-inf
    constant). The phenotype is analysed on the unit-variance scale; effect
    sizes are returned in phenotype units.

    ``heritability_method="he"`` reuses the selected alpha's existing HE
    products to fit the covariance ``vg K + ve (I - QQ')``, with nonnegative
    variance components. It reports the unconstrained estimates, boundary
    status and Monte Carlo precision in ``extra["he_variance"]``. This is
    single-component HE, not LDAK's partitioned estimator with large-effect
    exclusions. Alpha selection remains the same off-diagonal HE score.
    Nonidentified or numerically invalid HE fits raise an error; increase
    ``he_probes`` or use ``heritability_method="reml"`` in that case. REML
    remains the default. The existing variational noise floor of 0.001 on
    the unit-variance scale still applies, and is reported for HE fits.
    LDAK optionally revises its HE estimate by Monte
    Carlo REML (by default for smaller samples or strong structure).

    ``n_threads`` controls optional parallel genotype preparation and Numba
    coordinate sweeps across independent candidate models and LOCO fits.
    Model sweeps retain their SNP order; float64 genotype preparation may
    differ from NumPy by rounding at the final storage precision. The default
    is one; larger values require the ``fast`` extra and cannot exceed Numba's
    configured thread limit. Large fits also parallelize residual matrix
    products over sample rows, temporarily limiting detected BLAS pools to
    one thread. Apple Accelerate is not detected by threadpoolctl; set
    VECLIB_MAXIMUM_THREADS=1 before Python starts to avoid nested threading.

    ``denominator="spectral"`` (mixmogam extension, heuristic for the
    elastic-net statistic): with strong structure, U and T are multiplied by
    the structure-aware factors kappa_j of the weighted ridge model, T is
    then already calibrated (lambda' = 1), and lambda matches the two.
    """
    if heritability_method not in ("he", "reml"):
        raise ValueError("heritability_method must be 'he' or 'reml'")
    if heritability_method == "he" and alpha_method != "he":
        raise ValueError("heritability_method='he' requires alpha_method='he'")
    n_threads = _validate_n_threads(n_threads)
    rng = np.random.default_rng(random_state)
    st = _setup(y, gt, X, max_loco_groups, block, cache_bytes, n_threads=n_threads)
    n, G = st.lg.n, st.lg.n_groups
    if denominator not in ("constant", "spectral"):
        raise ValueError(f"unknown denominator {denominator!r}")
    if linalg.norm(st.y_p) <= 10 * np.finfo(float).eps * linalg.norm(st.y):
        raise ValueError("phenotype has no residual variation after covariate adjustment")
    sy = float(np.sqrt(np.sum(st.y_p**2) / st.n_eff))
    if not np.isfinite(sy) or sy <= np.sqrt(np.finfo(float).tiny):
        raise ValueError("phenotype has no residual variation after covariate adjustment")
    ys = st.y_p / sy

    struct = structure_test(st, structure_snps, rng)

    # (1b) alpha, then variance components at that alpha.
    f = np.clip(st.lg.mean / 2.0, 1e-6, 1 - 1e-6)
    # Under strong structure the ridge solves below need a preconditioner:
    # REML computes its deflation basis explicitly so that they can share it.
    basis = None
    if alpha_method == "he":
        he = _he_alpha(st, f, list(alphas), he_probes, rng,
                       fit_h2=heritability_method == "he")
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
    h2j = w * h2 / float(np.sum(w))
    s2e = max(1.0 - h2, 1e-3)

    def prior_row(p, F):
        lam = np.sqrt(2.0 * p / (1.0 - F)) if p > 0 else 0.0
        v = F / (1.0 - p) if p < 1 else 0.0
        return (p, lam, v)

    # (1d) hold-out choice of (p, F)
    held = rng.choice(n, size=max(1, int(round(cv_fraction * n))), replace=False)
    folds = np.ones(n, dtype=np.int64)
    folds[held] = 0  # fold 0 is held out; fold 1 is never held out
    eng = VBEngine(st.lg, folds=folds, n_threads=n_threads)
    Pn = len(grid)
    prior = np.array([prior_row(*g_) for g_ in grid])
    cv = eng.fit(np.repeat(ys[:, None], Pn, axis=1), np.zeros(Pn, dtype=np.int64),
                 np.full(Pn, -1), PRIOR_ENET, prior, np.full(Pn, s2e), snp_scale=h2j,
                 max_iter=vb_max_iter, tol=vb_tol)
    pred = cv["prediction"]
    mse = np.mean((ys[held, None] - pred[held]) ** 2, axis=0)
    best = int(np.argmin(mse))

    # (1e) LOCO elastic-net scores and the OLS-on-offset statistics U
    loco = eng.fit(np.repeat(ys[:, None], G, axis=1), np.full(G, -1), np.arange(G),
                    PRIOR_ENET, np.tile(prior_row(*grid[best]), (G, 1)), np.full(G, s2e),
                    snp_scale=h2j, max_iter=vb_max_iter, tol=vb_tol)
    del loco["prediction"]
    rsU = _retro_stats(st, loco["resid"])
    U = rsU["chi2"]

    extra = {"method": "kvik", "denominator": denominator, "alpha": alpha, "h2": h2,
             "alpha_method": alpha_method, "alpha_scores": alpha_scores,
             "structure": struct,
             "cv_grid": list(grid), "cv_mse": mse, "cv_best": grid[best],
             "cv_converged": cv["converged"],
             "loco_iterations": loco["iterations"], "loco_converged": loco["converged"]}
    if heritability_method == "he":
        extra.update({"heritability_method": "he", "he_variance": he["variance_fit"],
                      "vb_residual_variance": s2e,
                      "vb_noise_floor_applied": bool(1.0 - h2 < 1e-3)})
    lam = 1.0
    if struct["strong"] and h2 > 0:
        # ridge-prior LOCO scores: their OLS statistic is the BOLT-LMM-inf
        # retrospective statistic, lambda' its prospective calibration
        # At the HE residual-zero boundary, use the same declared noise
        # floor for ridge calibration as for VB. A zero ridge would leave
        # the projected LOCO operator singular in the covariate directions.
        residual_variance = s2e if heritability_method == "he" else 1.0 - h2
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
    extra["lambda"] = lam
    chi2 = lam * U
    with np.errstate(divide="ignore", invalid="ignore"):
        beta_z = lam * rsU["num"] / rsU["zz"] * sy
        se_z = np.abs(beta_z) / np.sqrt(chi2)
    return _result(st, gt, chi2, beta_z, se_z, extra)
