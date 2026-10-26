"""Stage-1 whole-genome predictors for the two-step engines.

The two-step statistic of :mod:`mixmogam.twostep` tests each SNP against a
LOCO residual from which the polygenic effects of the other chromosomes have
been conditioned out; the predictor that supplies those effects is the
stage-1 model. This module holds the contract both implementations share and
an optional LDpred3 implementation of it.

The contract. A stage-1 fit returns a dictionary with the keys
:func:`mixmogam.twostep._hratt_fit_scores` returns: ``resid`` (n, groups)
LOCO residuals, ``mse`` (candidates,) hold-out mean squared errors, ``best``
the chosen candidate, and the convergence records; ``prediction`` (n,
groups) and ``offset_slope`` (groups,) accompany cross-fitted fits. The
residual construction is shared (:func:`residualize_score`): each group's
LOCO score enters with its least-squares slope on the projected phenotype
kept in [0, 1].

Effect units. Effects are per unit of the covariate-projected standardized
rows of :class:`mixmogam._loco.LocoGenotypes` -- the rows
:func:`mixmogam.twostep._out_of_fold_loco` multiplies them with -- counting
the allele whose dosage ``Genotypes.G`` encodes. LDpred3 fits effects per
*unit-variance* column (:func:`ld_blocks` standardizes the rows to build its
correlation blocks), so the fitted effects divide by the per-variant column
scale :func:`column_scale` before leaving this module.

LDpred3 stage-1 (:func:`ldpred3_stage1_scores`). LDpred3 fits posterior-mean
effects from standardized marginal effects and per-block LD (its
``beta_hat = R beta + N(0, R/N)`` working model on standardized effects).
The summary statistics here are covariate-adjusted marginal regressions of
the working phenotype on the projected standardized rows
(:func:`fold_marginals`), one per cross-fit training set so that no sample's
phenotype enters its own score, and the LD blocks are the correlations of
those same rows (:func:`ld_blocks`; LDpred3's own panel builder expects
non-negative dosages with negatives as missing calls and cannot carry them).
The model is block-independent given its hyperparameters, so the LOCO score
of a sample is exactly its genome-wide fitted value less its own group's
share: the subtraction in
:func:`~mixmogam.twostep._out_of_fold_loco` *is* the refit, unlike the
variational engine whose block updates interact through the residual.
Three consequences of that model are accepted and documented here:

1. Cross-block LD is ignored; long-range-LD regions are fitted blockwise.
2. Shared hyperparameters (the genome-wide ``h2`` of the infinitesimal fit,
   the ``(s, lambda)`` lassosum2 selects across blocks) leak the tested
   chromosome's variants into the other chromosomes' effects second-order,
   the same order BOLT-LMM and REGENIE accept. LOCO groups are kept
   block-pure so nothing else crosses.
3. The LD reference is the analysis cohort itself and is shared across
   folds: only genotypes, never phenotypes, of a held-out fold enter its
   prediction through the LD, which is the usual unsupervised-preprocessing
   allowance.

The fitters are deterministic (``ldpred3_inf``'s closed form, lassosum2's
coordinate descent with pseudo-validation on the fold's own rows), so the
stage-1 statistics carry no Monte Carlo noise. LDpred3 is an optional
dependency: it is imported only when a ``stage1='ldpred3-*'`` fit is
requested, and remains outside this package's install requirements.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["STAGE1_METHODS", "MarginalFits", "fold_marginals", "column_scale",
           "residualize_score", "ld_blocks", "ldpred3_stage1_scores"]

STAGE1_METHODS = ("enet", "ldpred3-inf", "ldpred3-lassosum2", "auto")


@dataclass
class MarginalFits:
    """Covariate-adjusted marginal regressions of the working phenotype on
    the projected standardized rows, one per training set.

    ``beta`` and ``se`` are (m, sets) OLS slopes and standard errors of the
    regression through the projected space (no intercept: the covariates,
    intercept included, are projected out), ``n_eff`` (sets,) the rows less
    the covariate count, ``labels`` (sets,) names the sets, and ``scale``
    (m,) is the per-variant scale of the unit-variance columns LDpred3 fits
    in (:func:`column_scale`).
    """

    beta: np.ndarray
    se: np.ndarray
    n_eff: np.ndarray
    labels: tuple
    scale: np.ndarray


def _held_mask(n: int, held) -> np.ndarray:
    """Boolean membership mask of the hold-out ``held`` (indices or mask)."""
    held = np.asarray(held)
    mask = np.zeros(n, dtype=bool)
    if held.ndim != 1:
        raise ValueError("held must be one-dimensional")
    if held.size and held.dtype.kind == "b":
        if held.size != n:
            raise ValueError("a boolean held mask needs one entry per sample")
        mask[:] = held
    else:
        mask[np.asarray(held, dtype=np.int64)] = True
    return mask


def fold_marginals(lg, ys: np.ndarray, folds, held) -> MarginalFits:
    """Marginal regressions per cross-fit training set and their complements.

    The sets are the training complements of each fold (``folds`` labels
    0..k-1 per sample) and of the model-selection hold-out ``held`` (indices),
    then each set's complement -- the validation rows lassosum2 selects its
    penalty on. One genotype pass accumulates every set's cross-products, so
    k additional fits cost no additional genotype traffic beyond this.

    ``held`` rows take part in the fold training complements: the hold-out
    selects the prior, the fold fits see every phenotype except their own
    fold's, as in the variational path. The covariate basis is the whole
    cohort's (genotypes and design only enter it), so a training subset's
    regression is the projected slope on its own rows.
    """
    ys = np.asarray(ys, dtype=np.float64).ravel()
    if ys.shape != (lg.n,):
        raise ValueError("ys must hold one projected phenotype value per sample")
    held_mask = _held_mask(lg.n, held)
    folds = None if folds is None else np.asarray(folds, dtype=np.int64)
    if folds is not None and folds.shape != (lg.n,):
        raise ValueError("folds must hold one label per sample")
    k = 0 if folds is None else int(folds.max()) + 1

    masks, labels = [], []
    for f in range(k):
        masks.append(folds != f)
        labels.append(f"train: complement of fold {f}")
    masks.append(~held_mask)
    labels.append("train: complement of the hold-out")
    for f in range(k):
        masks.append(folds == f)
        labels.append(f"validation: fold {f}")
    masks.append(held_mask)
    labels.append("validation: hold-out")
    M = np.column_stack([mask.astype(np.float64) for mask in masks])
    counts = M.sum(axis=0)
    q = lg.Q.shape[1]
    n_eff = counts - q
    if np.any(n_eff <= 0):
        bad = [labels[i] for i in np.flatnonzero(n_eff <= 0)]
        raise ValueError(f"training sets too small for {q} covariates: {', '.join(bad)}")

    m, t = lg.m, M.shape[1]
    W = ys[:, None] * M
    s1 = np.zeros((m, t))
    s2 = np.zeros((m, t))
    mean = np.zeros(m)
    tile = max(1, int(2**21 // max(lg.n, 1)))  # ~16 MiB of float64 rows
    for start in range(0, m, tile):
        idx = np.arange(start, min(start + tile, m))
        Z = lg.rows(idx)  # (b, n) projected standardized rows
        s1[idx] = Z @ W
        s2[idx] = (Z * Z) @ M
        mean[idx] = Z.mean(axis=1)
    beta = np.where(s2 > 0, s1 / np.where(s2 > 0, s2, 1.0), 0.0)
    y2 = (ys * ys) @ M
    rss = np.maximum(y2[None, :] - beta * s1, 0.0)
    df = counts - q - 1
    if np.any(df <= 0):
        raise ValueError("training sets leave no residual degrees of freedom")
    se2 = (rss / df[None, :]) / np.where(s2 > 0, s2, 1.0)
    se = np.sqrt(np.where(s2 > 0, se2, 0.0))
    return MarginalFits(beta=beta, se=se, n_eff=n_eff, labels=tuple(labels),
                        scale=column_scale(lg, mean))


def column_scale(lg, mean: np.ndarray) -> np.ndarray:
    """Per-variant scale of the unit-variance columns LDpred3 fits in: the
    centered root-mean-square of the projected standardized rows,
    ``sqrt((sum z^2 - n mean^2) / n)``, from the prepared norms ``lg.zz`` and
    the column means. Effects fit per unit-variance column divide by this to
    land in the :class:`mixmogam._loco.LocoGenotypes` row units; the mean
    subtraction the standardization also applies shifts predictions by one
    constant per fit, which the covariate projection of
    :func:`residualize_score` removes. Monomorphic projected columns get
    scale 0 and carry no effects. The prepared norms and a pass's own sums
    accumulate in different orders; the two agree to float summation
    precision.
    """
    mean = np.asarray(mean, dtype=np.float64)
    n = lg.n
    ss = lg.zz - n * mean * mean
    return np.sqrt(np.maximum(ss, 0.0) / max(n, 1))


def residualize_score(lg, ys: np.ndarray, prediction: np.ndarray) -> dict:
    """LOCO scores enter the residual with their least-squares slope in [0, 1].

    ``prediction`` (n, groups) holds each sample's LOCO scores; each group's
    slope of the projected phenotype on its projected score shrinks an
    over-dispersed score and never undoes the prior's shrinkage (a slope
    above one) or flips an anti-predicting one (below zero). Returns
    ``{"resid", "offset_slope"}``.
    """
    ys = np.asarray(ys, dtype=np.float64).ravel()
    Pp = lg.project(np.asarray(prediction, dtype=np.float64))
    norm2 = np.einsum("ij,ij->j", Pp, Pp)
    slope = np.zeros(Pp.shape[1])
    np.divide(Pp.T @ ys, norm2, out=slope, where=norm2 > 0)
    np.clip(slope, 0.0, 1.0, out=slope)
    return {"resid": ys[:, None] - Pp * slope, "offset_slope": slope}


def ld_blocks(lg, order: np.ndarray, scale: np.ndarray, block_size: int = 500) -> list:
    """Per-block correlation matrices of the projected standardized rows.

    Returns ``[(R, idx)]`` in LDpred3's block format: ``R`` is the float32
    correlation matrix of the unit-variance columns of the variants ``idx``,
    where ``idx`` indexes the group-sorted permutation ``order`` (so blocks
    are group-pure and tile the permutation). The rows are standardized with
    the ``scale`` of :func:`column_scale`, the same columns the effect
    conversion undoes; monomorphic columns get an identity row and no
    effects. LDpred3's own ``compute_ld_blocks`` is not used: its panel
    input is a dosage convention (non-negative calls, negatives as missing)
    that signed standardized rows cannot pass, and these correlations of the
    projected rows are the covariate-adjusted LD the marginal regressions
    imply.
    """
    if not isinstance(block_size, (int, np.integer)) or block_size <= 0:
        raise ValueError("block_size must be a positive integer")
    order = np.asarray(order, dtype=np.int64)
    scale = np.asarray(scale, dtype=np.float64)[order]
    groups = lg.groups[order]  # non-decreasing: blocks stop at group changes
    blocks = []
    start = 0
    while start < order.size:
        stop = min(start + int(block_size), order.size)
        stop = min(stop, int(np.searchsorted(groups, groups[start], side="right")))
        idx = np.arange(start, stop)
        Z = lg.rows(order[idx])  # (b, n) projected standardized rows
        # rows() owns its float64 buffer: center and scale it in place.
        # Selecting all polymorphic rows twice for the GEMM otherwise
        # keeps several additional (block_size, n) arrays alive.
        Z -= Z.mean(axis=1, keepdims=True)
        good = scale[idx] > 0
        unit = Z if np.all(good) else Z[good]
        unit /= scale[idx][good, None]
        R = np.zeros((idx.size, idx.size))
        R[np.ix_(good, good)] = (unit @ unit.T) / lg.n
        np.clip(R, -1.0, 1.0, out=R)
        R = (R + R.T) / 2.0
        np.fill_diagonal(R, 1.0)
        blocks.append((R.astype(np.float32), idx))
        start = stop
    return blocks


def _import_ldpred3():
    try:
        import ldpred3
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError(
            "stage1 'ldpred3-*' needs the ldpred3 package "
            "(github.com/bvilhjal/ldpred3); install it separately "
            "(mixmogam keeps it out of its install requirements)") from exc
    return ldpred3


def _fit_effects(ldpred3, blocks, beta_std, n_eff, method, h2, fit_kwargs) -> np.ndarray:
    """One LDpred3 fit's standardized effects, in the block cover's order."""
    if method == "inf":
        if fit_kwargs:
            raise ValueError(f"stage1 'ldpred3-inf' takes no fit keywords, got {sorted(fit_kwargs)}")
        return np.asarray(ldpred3.ldpred3_by_blocks(
            blocks, beta_std, n_eff, method="inf", h2=h2), dtype=np.float64)
    return np.asarray(ldpred3.lassosum2(
        blocks, beta_std, **fit_kwargs).beta_est, dtype=np.float64)


def ldpred3_stage1_scores(st, ys: np.ndarray, h2: float, folds, held,
                          *, method: str, keep_prediction: bool = False,
                          block_size: int = 500, **fit_kwargs) -> dict:
    """Contract dictionary from cross-fitted LDpred3 LOCO scores.

    ``method`` is ``"ldpred3-inf"`` (the infinitesimal closed form at the
    genome-wide ``h2`` of the variance fit) or ``"ldpred3-lassosum2"``
    (penalty selected per fold by pseudo-validation on that fold's own rows).
    ``fit_kwargs`` (with ``block_size``) pass to lassosum2. The hold-out MSE
    of a fit excluding ``held`` is the candidate's model-selection score,
    on the same samples as the elastic-net candidates'.
    """
    lg = st.lg
    k = 0 if folds is None else int(np.max(np.asarray(folds))) + 1
    if not k:
        raise ValueError("stage1 'ldpred3-*' needs cross-fitted scores (loco_folds >= 2)")
    folds = np.asarray(folds, dtype=np.int64)
    if h2 <= 0:  # no polygenic variance: no scores, as in the elastic-net path
        held_mask = _held_mask(lg.n, held)
        ys = np.asarray(ys, dtype=np.float64).ravel()
        G = lg.n_groups
        out = {"resid": np.repeat(ys[:, None], G, axis=1),
               "offset_slope": np.zeros(G),
               "mse": np.array([float(np.mean(ys[held_mask] ** 2))]),
               "best": 0, "cv_converged": True, "loco_iterations": 0,
               "loco_converged": True, "cv_grid": [method], "cv_best": method}
        if keep_prediction:
            out["prediction"] = np.zeros((lg.n, G))
        return out
    ldpred3 = _import_ldpred3()
    marg = fold_marginals(lg, ys, folds, held)
    order = np.argsort(lg.groups, kind="stable")
    # Column indices of MarginalFits: k fold training complements, the
    # hold-out complement, then the matching validation sets.
    n_sets = len(marg.labels)
    train_cols = list(range(k + 1))
    val_cols = list(range(k + 1, n_sets))
    scale = marg.scale[order]
    blocks = ld_blocks(lg, order, marg.scale, block_size)

    effects = np.zeros((lg.m, len(train_cols)))
    for j, c in enumerate(train_cols):
        beta_std, _ = ldpred3.standardize_betas(
            marg.beta[:, c], marg.se[:, c], float(marg.n_eff[c]))
        kwargs = dict(fit_kwargs)
        if method == "ldpred3-lassosum2":
            val_std, _ = ldpred3.standardize_betas(
                marg.beta[:, val_cols[j]], marg.se[:, val_cols[j]],
                float(marg.n_eff[val_cols[j]]))
            kwargs["validation_beta"] = val_std[order]
        e = _fit_effects(ldpred3, blocks, beta_std[order], float(marg.n_eff[c]),
                         "inf" if method == "ldpred3-inf" else "lassosum2",
                         h2, kwargs)
        with np.errstate(divide="ignore", invalid="ignore"):
            e = np.where(scale > 0, e / np.where(scale > 0, scale, 1.0), 0.0)
        effects[order, j] = e  # back to variant index order

    # Held-out candidate score: the fit excluding the hold-out, genome-wide.
    from mixmogam.twostep import _out_of_fold_loco  # late: the patchable seam
    held_mask = _held_mask(lg.n, held)
    fitted = lg.project(_fitted_unprojected(lg, effects[:, k]))
    mse = np.array([float(np.mean((ys[held_mask] - fitted[held_mask]) ** 2))])

    prediction = _out_of_fold_loco(lg, effects[:, :k], folds)
    out = {"mse": mse, "best": 0, "cv_converged": True,
           "loco_iterations": 0, "loco_converged": True,
           "cv_grid": [method], "cv_best": method}
    out.update(residualize_score(lg, ys, prediction))
    if keep_prediction:
        out["prediction"] = prediction
    return out


def _fitted_unprojected(lg, beta: np.ndarray) -> np.ndarray:
    """z' beta over one pass of the unprojected rows (the projection of the
    result is the projected fitted value)."""
    out = np.zeros(lg.n)
    for take, g, z in lg.raw_slices():
        out += z.T @ beta[take]
    return out
