"""Stepwise multi-locus mixed model (MLMM), v1's signature feature.

Forward selection with BIC-based model-size choice and backward
elimination, each step rescanning with the current cofactors as
covariates on the batched engine. The stepwise criterion follows v1:
BIC = log-likelihood - (k/2) log(n) with the model REML likelihood.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from mixmogam.lmm import LMM

__all__ = ["mlmm", "mlmm_bic"]


def _bic(ll: float, n_par: int, n: int) -> float:
    return ll - 0.5 * n_par * np.log(n)


def _fit_ll(y, X, K) -> float:
    return LMM(y, X=X, K=K).fit().ll


def mlmm(
    y,
    gt,
    K: Optional[np.ndarray] = None,
    max_steps: int = 10,
    candidate_fraction: float = 0.02,
    min_maf: float = 0.05,
    backward: bool = True,
    verbose: bool = False,
) -> dict:
    """Forward-backward multi-locus mixed model.

    Parameters
    ----------
    y : (n,) phenotype
    gt : Genotypes (sample-major container)
    K : kinship; None gives the plain linear-model stepwise
    max_steps : maximum forward steps
    candidate_fraction : fraction of the genome scanned per forward step
        (v1 scanned everything; the fraction caps cost and mirrors its
        cofactor-candidate heuristics)

    Returns
    -------
    dict with ``cofactors`` (variant indices), ``steps`` (per-step log),
    ``best_step`` (BIC-selected model size), ``bics``.
    """
    y = np.asarray(y, dtype=np.float64).ravel()
    n = y.size
    G_all = gt.G.astype(np.float64)
    cofactors: list[int] = []
    steps = []
    lls = [_fit_ll(y, None, K)]
    bics = [_bic(lls[0], 1, n)]
    best_step = 0
    p_threshold = 1.0 / gt.n_variants

    for step in range(max_steps):
        X = G_all[:, cofactors] if cofactors else None
        lmm = LMM(y, X=X, K=K)
        lmm.fit()
        scan = lmm.scan(gt, dtype=np.float64)
        order = np.argsort(scan["ps"])
        n_cand = max(int(candidate_fraction * gt.n_variants), 50)
        chosen = None
        for j in order[:n_cand]:
            if j in cofactors:
                continue
            chosen = int(j)
            break
        if chosen is None or scan["ps"][chosen] > p_threshold:
            break  # nothing even nominally significant left
        cofactors.append(chosen)
        ll = _fit_ll(y, G_all[:, cofactors], K)
        lls.append(ll)
        bic = _bic(ll, len(cofactors) + 1, n)
        bics.append(bic)
        steps.append(
            {"step": step + 1, "added": chosen, "ll": ll, "bic": bic,
             "p": float(scan["ps"][chosen])}
        )
        if verbose:
            print(f"step {step+1}: +{chosen} p={scan['ps'][chosen]:.2e} bic={bic:.2f}")
        if bic > bics[best_step]:
            best_step = len(cofactors)

    if backward and len(cofactors) > 1:
        # drop cofactors that no longer help
        improved = True
        while improved and len(cofactors) > 1:
            improved = False
            for c in list(cofactors):
                trial = [x for x in cofactors if x != c]
                ll = _fit_ll(y, G_all[:, trial], K)
                if _bic(ll, len(trial) + 1, n) > _bic(
                    lls[len(cofactors)], len(cofactors) + 1, n
                ):
                    cofactors = trial
                    improved = True
                    break

    return {
        "cofactors": cofactors,
        "steps": steps,
        "bics": bics,
        "best_step": best_step,
    }


def mlmm_bic(y, gt, K=None, max_steps: int = 10) -> dict:
    """Convenience wrapper returning the BIC-optimal cofactor set."""
    res = mlmm(y, gt, K=K, max_steps=max_steps)
    keep = res["cofactors"][: res["best_step"]] if res["best_step"] else []
    res["cofactors_optimal"] = keep
    return res
