"""Stepwise multi-locus mixed model (MLMM; Segura et al. 2012, Nat Genet).

Forward inclusion of the most significant SNP as a fixed-effect cofactor,
re-estimating the variance components after every step, then backward
elimination from the last forward model by dropping the least significant
cofactor. Every visited model is scored and the final model is chosen by

- ``"ebic"``: extended BIC (Chen & Chen 2008) -2 ll + k' log n + 2 log C(m, k),
- ``"mbonf"``: the largest model whose cofactors all have conditional
  p-values below the Bonferroni threshold alpha / m (ties between a forward
  and a backward model of equal size go to the smaller maximum p-value),
- ``"bic"``: plain BIC, reported for completeness. Segura et al. found it
  too tolerant for GWAS and recommend EBIC or mBonf.

Log-likelihoods for the criteria are maximum likelihood (ML): REML
likelihoods are not comparable between models with different fixed
effects. The scans and conditional cofactor tests use the REML variance
components, as in v1. This is a port of the v1 ``linear_models.mlmm``
semantics (preserved under the ``v1.0-legacy`` tag) on the batched engine.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
from scipy import special, stats

from mixmogam.genotypes import MISSING
from mixmogam.lmm import LMM

__all__ = ["mlmm", "CRITERIA"]

CRITERIA = ("ebic", "mbonf", "bic")


def _dosage(gt, j: int) -> np.ndarray:
    """Mean-imputed float64 dosage column of variant ``j``."""
    g = np.asarray(gt.G[:, j]).astype(np.float64)
    miss = g == MISSING
    if miss.any():
        g[miss] = g[~miss].mean() if (~miss).any() else 0.0
    return g


def _log_choose(m: int, k: int) -> float:
    return float(special.gammaln(m + 1) - special.gammaln(k + 1) - special.gammaln(m - k + 1))


class _Model:
    """One visited model: REML fit for scanning, ML log-likelihood for scoring."""

    def __init__(self, y, X0, cof_cols, K, base: Optional[LMM]):
        X = X0 if not cof_cols else np.column_stack([X0] + cof_cols)
        self.lmm = LMM(y, X=X, K=K, add_intercept=False)
        if base is not None and K is not None:
            self.lmm._eig = base._eig  # spectrum of K is X-independent
            self.lmm._basis_cache = base._basis_cache
        self.ll = self.lmm.fit(method="ml").ll
        self.fit = self.lmm.fit(method="reml", recompute=True)

    def cofactor_pvalues(self, n_cof: int) -> np.ndarray:
        """Conditional F-test p-values of the last ``n_cof`` columns of X.

        One GLS regression at the REML variance components: the drop-one
        F statistic of column i is beta_i^2 / (s^2 [(X'V^-1 X)^-1]_ii).
        """
        if n_cof == 0:
            return np.empty(0)
        lmm, fit = self.lmm, self.fit
        Xw = lmm._apply_inv_sqrt(lmm.X, fit.delta, np.float64)
        yw = lmm._apply_inv_sqrt(lmm.y, fit.delta, np.float64)
        XtX = Xw.T @ Xw
        beta = np.linalg.solve(XtX, Xw.T @ yw)
        resid = yw - Xw @ beta
        df = lmm.n - lmm.q
        s2 = float(resid @ resid) / df
        diag = np.diag(np.linalg.inv(XtX))[-n_cof:]
        f = beta[-n_cof:] ** 2 / (s2 * diag)
        return stats.f.sf(f, 1, df)


def mlmm(
    y,
    gt,
    K: Optional[np.ndarray] = None,
    X: Optional[np.ndarray] = None,
    max_steps: int = 10,
    alpha: float = 0.05,
    h2_stop: float = 1e-3,
    backward: bool = True,
    dtype=np.float64,
    verbose: bool = False,
) -> dict:
    """Forward-backward multi-locus mixed model.

    Parameters
    ----------
    y : (n,) phenotype
    gt : Genotypes (sample-major container); missing calls are mean-imputed
    K : kinship; None gives the stepwise linear model (SWLM)
    X : extra covariates (an intercept is always included)
    max_steps : maximum number of forward steps
    alpha : family-wise level of the mBonf criterion (threshold alpha / m)
    h2_stop : forward inclusion stops once the pseudo-heritability has
        been below this value for two consecutive steps (v1 semantics)
    backward : run backward elimination from the last forward model

    Returns
    -------
    dict with ``steps`` (one record per visited model: ``action`` "start",
    "+" or "-", ``cofactors`` as variant indices, ``ll`` (ML), ``bic``,
    ``ebic``, ``max_cof_p``, ``pseudo_heritability`` and, for forward
    models, ``min_p`` of the next scan), ``selected`` (criterion ->
    cofactor list), ``selected_step`` (criterion -> index into ``steps``)
    and ``cofactors`` (the EBIC selection).
    """
    y = np.asarray(y, dtype=np.float64).ravel()
    n = y.size
    m = gt.n_variants
    X0 = np.ones((n, 1)) if X is None else np.column_stack(
        [np.ones(n), np.asarray(X, dtype=np.float64).reshape(n, -1)]
    )
    q0 = X0.shape[1]
    threshold = alpha / m
    log_n = np.log(n)

    base = LMM(y, X=X0, K=K, add_intercept=False) if K is not None else None
    if base is not None:
        base.eigen()  # one eigendecomposition shared by every visited model

    steps: list[dict] = []

    def record(action, cofs, model, max_p, extra=None):
        k = len(cofs)
        n_par = q0 + 1 + k  # covariates, variance scale, cofactors
        bic = -2.0 * model.ll + n_par * log_n
        rec = {
            "step": len(steps),
            "action": action,
            "cofactors": list(cofs),
            "ll": model.ll,
            "bic": bic,
            "ebic": bic + 2.0 * _log_choose(m, k),
            "max_cof_p": max_p,
            "pseudo_heritability": model.fit.pseudo_heritability,
        }
        if extra:
            rec.update(extra)
        steps.append(rec)
        if verbose:
            print(f"step {rec['step']} {action} k={k} ll={model.ll:.2f} "
                  f"ebic={rec['ebic']:.2f} max_cof_p={max_p:.2e} "
                  f"h2={rec['pseudo_heritability']:.3f}")
        return rec

    cofs: list[int] = []
    cols: list[np.ndarray] = []
    model = _Model(y, X0, cols, K, base)
    rec = record("start", cofs, model, 0.0)
    low_h2 = 0
    for _ in range(max_steps):
        if n - model.lmm.q <= 1:
            break
        scan = model.lmm.scan(gt, dtype=dtype)
        ps = np.asarray(scan["ps"], dtype=np.float64).copy()
        ps[cofs] = np.nan
        if not np.isfinite(ps).any():
            break
        j = int(np.nanargmin(ps))
        rec["min_p"] = float(ps[j])
        cofs.append(j)
        cols.append(_dosage(gt, j))
        model = _Model(y, X0, cols, K, base)
        cof_p = model.cofactor_pvalues(len(cofs))
        rec = record("+", cofs, model, float(cof_p.max()))
        if K is not None and model.fit.pseudo_heritability < h2_stop:
            low_h2 += 1
            if low_h2 >= 2:
                break
        else:
            low_h2 = 0

    if backward:
        while len(cofs) > 1:
            cof_p = model.cofactor_pvalues(len(cofs))
            drop = int(np.argmax(cof_p))  # least significant = smallest F
            del cofs[drop]
            del cols[drop]
            model = _Model(y, X0, cols, K, base)
            cof_p = model.cofactor_pvalues(len(cofs))
            record("-", cofs, model, float(cof_p.max()))

    selected_step = {}
    for crit in ("ebic", "bic"):
        selected_step[crit] = int(np.argmin([s[crit] for s in steps]))
    best = None
    for s in steps:
        if s["max_cof_p"] >= threshold:
            continue
        key = (len(s["cofactors"]), -s["max_cof_p"])
        if best is None or key > best[0]:
            best = (key, s["step"])
    selected_step["mbonf"] = best[1]
    selected = {c: list(steps[i]["cofactors"]) for c, i in selected_step.items()}
    return {
        "steps": steps,
        "selected": selected,
        "selected_step": selected_step,
        "cofactors": selected["ebic"],
        "threshold": threshold,
    }
