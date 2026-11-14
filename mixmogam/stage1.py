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
The model is block-independent given its hyperparameters, so for
``ldpred3-inf`` (shared genome-wide h2) the LOCO score of a sample is
exactly its genome-wide fitted value less its own group's share: the
subtraction in :func:`~mixmogam.twostep._out_of_fold_loco` *is* the
leave-chromosome-out refit. ``ldpred3-lassosum2`` and ``ldpred3-fast``
choose hyperparameters on all chromosomes, so that subtraction is not
a full refit; with evidence of relatives, :func:`ldpred3_stage1_scores`
re-fits leave-group-out (``refit=True``), the same k x groups pattern
as the elastic net. The variational engine always needs that refit,
because its block updates interact through the residual.
Three consequences of that model are accepted and documented here:

1. Cross-block LD is ignored; long-range-LD regions are fitted blockwise.
2. Shared hyperparameters (the genome-wide ``h2`` of the infinitesimal fit,
   the ``(s, lambda)`` lassosum2 selects across blocks) leak the tested
   chromosome's variants into the other chromosomes' effects. For inf
   that leak is second-order, the same order BOLT-LMM and REGENIE
   accept, and the relatives auto-refit is skipped. For lassosum2 and
   fast the leak is larger (in-model null mean chi-square 1.122 without
   refit in review simulations); the relatives auto-refit is applied.
   LOCO groups are kept block-pure so nothing else crosses.
3. The LD reference is the analysis cohort itself and is shared across
   folds: only genotypes, never phenotypes, of a held-out fold enter its
   prediction through the LD, which is the usual unsupervised-preprocessing
   allowance.

Per-fold LD (``stage1_kwargs["ld_split"] = "train"``). By default (``"all"``)
every fold's fit reads one cover built from all n rows, while its marginals
use only that fold's training rows. ``"train"`` gives each training set its
own cover: the sample partition (fold x hold-out membership) has 2k parts,
each part's Gram (and, for hard calls, its call and covariate products) is
accumulated in the same single pass, and a training set's Gram is the total
less its excluded parts, re-centred and renormalised on its own rows. Build
flops are unchanged (the parts partition the samples); storage and cache are
k + 1 covers, keyed additionally by the partition, source calls and covariate
projection. Only exact matches can reuse training covers; whole-cohort
``"all"`` covers retain their approximate reuse with overlap diagnostics.

D8 spectrum (``stage1_kwargs["ld_psd"]``). D8 rounding adds a symmetric error
E with ||E||_2 about 2 sigma sqrt(k), sigma = 1 / (127 sqrt 12), and pushes
the smallest eigenvalue of a (near-)singular block below zero. lassosum2
fits (1 - s) R + s I, PSD iff s / (1 - s) >= -lambda_min; the default
smallest s = 0.2 tolerates lambda_min >= -0.25. ``"check"`` (the D8 default)
records each block's lambda_min (LAPACK up to 2,000 variants, Lanczos
above) and warns when that tolerance is exceeded; ``"floor"`` raises
lassosum2's smallest s to (c / (1 + c)), c = -min lambda_min + 1e-3 (exact,
no requantising), and for the other fits replaces each negative block by
(R + c I) / (1 + c), requantised and checked; ``"ridge"`` does the latter for
every fit; ``"off"`` skips the check (the f32 default).

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

STAGE1_METHODS = ("enet", "ldpred3-inf", "ldpred3-lassosum2", "ldpred3-fast", "auto")

# stage1_kwargs keys this module consumes (the rest go to the fitter):
# ``ld`` (see :mod:`mixmogam._ld8`), ``ld_cache`` (None, an LDCache or a
# directory), ``ld_tile`` (samples per hard-call tile), ``ld_split`` and
# ``ld_psd`` (module notes), and for
# ``ldpred3-fast`` the layout and chain controls below.
_LD_KEYS = ("ld", "ld_cache", "ld_tile", "ld_split", "ld_psd")
_FAST_KEYS = ("width", "landmarks", "reach", "burn_in", "num_iter", "seed", "fast_block_size",
              "p_init", "layout_seed")
_FAST_DEFAULTS = {"width": 2500, "landmarks": 64, "fast_block_size": 6000, "burn_in": 300,
                  "num_iter": 100, "seed": 0, "p_init": 0.02, "reach": None, "layout_seed": 0}


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
    from mixmogam._ld8 import _validate_block_size
    block_size = _validate_block_size(block_size)
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


def _chain_mixing(path) -> dict:
    """Single-chain diagnostics of a sampled scalar path: its two halves
    treated as two chains (ldpred3.diagnostics.scalar_mixing: rank-normalised
    split R-hat and bulk ESS)."""
    x = np.asarray(path, dtype=np.float64).ravel()
    x = x[np.isfinite(x)]
    n = (x.size // 4) * 4
    if n < 16:
        return {"rhat": float("nan"), "ess_bulk": float("nan")}
    try:
        from ldpred3.diagnostics import scalar_mixing
        mix = scalar_mixing({"v": x[-n:].reshape(2, n // 2)})
        return {"rhat": float(mix["rhat"]["v"]), "ess_bulk": float(mix["ess_bulk"]["v"])}
    except Exception:  # noqa: BLE001  (diagnostics only)
        a, b = x[-n:].reshape(2, n // 2)
        w = (a.var(ddof=1) + b.var(ddof=1)) / 2
        bvar = (n // 2) * np.var([a.mean(), b.mean()], ddof=1)
        vp = (n // 2 - 1) / (n // 2) * w + bvar / (n // 2)
        return {"rhat": float(np.sqrt(vp / w)) if w > 0 else float("nan"),
                "ess_bulk": float("nan")}


def _fit_effects(ldpred3, blocks, beta_std, n_eff, method, h2, fit_kwargs, diag=None):
    """One LDpred3 fit's standardized effects, in the block cover's order;
    ``diag`` (a list) collects the fast chain's h2 / p estimates and mixing."""
    if method == "fast":
        import warnings
        kw = dict(fit_kwargs)
        # ldpred3 0.10.4 emits this once per fit; HRATT already warned once
        # per stage-1 call (below), so swallow the per-fold copies.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="LDpred3-fast is experimental", category=UserWarning)
            res = ldpred3.ldpred3_fast(
                blocks, beta_std, n_eff, h2_init=float(min(max(h2, 1e-3), 0.99)),
                p_init=kw["p_init"], burn_in=kw["burn_in"], num_iter=kw["num_iter"],
                seed=kw["seed"], width=kw["width"], landmarks=kw["landmarks"],
                reach=kw["reach"], layout_seed=kw["layout_seed"])
        if diag is not None:
            h2m = _chain_mixing(res.h2_path) if res.h2_path is not None else {}
            p_path = (res.p_rao_blackwell_path if res.p_rao_blackwell_path is not None
                      else res.p_path)
            pm = _chain_mixing(p_path) if p_path is not None else {}
            diag.append({"h2_est": float(res.h2_est), "p_est": float(res.p_est),
                         "n_iter": int(res.n_iter), "rhat_h2": h2m.get("rhat"),
                         "ess_h2": h2m.get("ess_bulk"), "rhat_p": pm.get("rhat"),
                         "ess_p": pm.get("ess_bulk")})
        return np.asarray(res.beta_est, dtype=np.float64)
    if method == "inf":
        if fit_kwargs:
            raise ValueError(f"stage1 'ldpred3-inf' takes no fit keywords, got {sorted(fit_kwargs)}")
        return np.asarray(ldpred3.ldpred3_by_blocks(
            blocks, beta_std, n_eff, method="inf", h2=h2), dtype=np.float64)
    return np.asarray(ldpred3.lassosum2(
        blocks, beta_std, **fit_kwargs).beta_est, dtype=np.float64)


def _stage1_blocks(lg, order, scale, method, block_size, ld, ld_cache, ld_tile, fast_opt,
                   ldpred3, info, split=None):
    """The LD covers a fit reads, from ``ld_cache`` when it holds them: one
    cover (``split`` None, the whole cohort's rows) or one per training set
    (``split = (labels, n_parts, excluded)``, :func:`mixmogam._ld8.split_parts`)."""
    import time

    from mixmogam import _ld8

    cache = _ld8.resolve_cache(ld_cache)
    size = (_ld8._validate_block_size(fast_opt["fast_block_size"], "fast_block_size")
            if method == "ldpred3-fast" else _ld8._validate_block_size(block_size))
    extra = None
    if method == "ldpred3-fast":
        extra = ("fast", fast_opt["width"], fast_opt["landmarks"], fast_opt["reach"],
                 fast_opt["layout_seed"], fast_opt.get("_psd"))
    n_covers = 1 if split is None else len(split[2])
    if split is not None and cache is not None:
        # Exact training LD needs all calls and the projection, not just
        # the sampled fingerprints used to estimate cohort overlap.
        from mixmogam._stream import panel_fingerprint
        extra = (extra, "split", _ld8.split_cache_extra(split[0]), panel_fingerprint(lg.gt),
                 _ld8._covariate_digest(lg))
    keys = None
    if cache is not None:
        base = _ld8.ld_cache_key(lg, order, size, ld, extra)
        keys = [base] if split is None else [f"{base}-s{j}" for j in range(n_covers)]
        t0 = time.perf_counter()
        got = [cache.get(key) for key in keys]
        if all(c is not None for c in got):
            chk = _ld8.reuse_check(lg, got[0][1])
            info.update(ld_cache="hit", ld_key=keys[0], ld_load_s=time.perf_counter() - t0,
                        ld_bytes=sum(_ld8.cover_nbytes(c[0]) for c in got),
                        ld_cache_overlap=chk["overlap"], ld_cache_coverage=chk["coverage"],
                        ld_cache_fst=chk["fst"], ld_cache_exact=chk["exact"],
                        ld_cache_covariates_match=chk["covariates_match"])
            reusable = chk["exact"] if split is not None else chk["overlap"] > 0.0
            if reusable:
                for w in chk["warnings"]:
                    import warnings
                    warnings.warn(f"LD cache reuse: {w}", RuntimeWarning, stacklevel=4)
                return [c[0] for c in got]
            info["ld_cache"] = "stale"  # train covers need exact reuse; all needs overlap
    c0, t0, cells0 = time.process_time(), time.perf_counter(), _ld8.COUNTER["cells"]
    if split is None:
        covers = [_ld8.build_blocks(lg, order, scale, size, ld, ld_tile)]
    else:
        covers = _ld8.build_split_blocks(lg, order, size, ld, split[0], split[1], split[2],
                                         ld_tile)
    info.update(ld_build_s=time.perf_counter() - t0, ld_build_cpu_s=time.process_time() - c0,
                ld_hardcall_cells=_ld8.COUNTER["cells"] - cells0)
    if method == "ldpred3-fast":
        c0, t0 = time.process_time(), time.perf_counter()
        covers = [_psd_cover(c, fast_opt.get("_psd", "off"), info, method, None) for c in covers]
        covers = [ldpred3.fast_ld_blocks(c, width=fast_opt["width"],
                                         landmarks=fast_opt["landmarks"],
                                         reach=fast_opt["reach"], seed=fast_opt["layout_seed"])
                  for c in covers]
        info["ld_build_s"] += time.perf_counter() - t0
        info["ld_build_cpu_s"] += time.process_time() - c0
    info.update(ld_cache=("rebuilt" if info.get("ld_cache") == "stale" else "miss")
                if cache is not None else "off",
                ld_bytes=sum(_ld8.cover_nbytes(c) for c in covers))
    if cache is not None:
        t0 = time.perf_counter()
        meta = _ld8.cover_meta(lg)
        for key, c in zip(keys, covers):
            cache.put(key, c, meta)
        info.update(ld_key=keys[0], ld_store_s=time.perf_counter() - t0)
    return covers


_PSD_MODES = ("off", "check", "floor", "ridge")
_PSD_EPS = 1e-3


def _psd_cover(blocks, mode, info, method, ctx):
    """Spectrum check / floor of one dense cover (``stage1_kwargs['ld_psd']``).

    ``"check"`` records each block's lambda_min (LAPACK below 2,000
    variants, Lanczos above) in ``info['ld_spectrum']``; ``"ridge"`` also
    floors every block with lambda_min < 0 to (R + c I) / (1 + c),
    requantised for D8, with c doubled from -lambda_min + 1e-3 until the
    stored block is PSD; ``"floor"`` is ``"ridge"`` except for lassosum2,
    whose own LD shrinkage (1 - s) R + s I is PSD for s / (1 - s) >=
    -lambda_min, so it raises the smallest s instead (exact, no requantising;
    :func:`_lassosum2_s_floor`)."""
    import time
    import warnings

    from mixmogam import _ld8

    if mode == "off" or not blocks or not all(isinstance(R, np.ndarray) for R, _ in blocks):
        return blocks
    t0 = time.perf_counter()
    spec = _ld8.cover_spectrum(blocks)
    lam = spec["lambda_min"]
    rec = info.setdefault("ld_spectrum", [])
    entry = {k: v for k, v in spec.items() if k != "lambda_min"}
    out = blocks
    if mode == "ridge" or (mode == "floor" and method != "ldpred3-lassosum2"):
        out, cs = [], []
        for (R, idx), l in zip(blocks, lam):
            R2, c, _ = _ld8.ridge_block_d8(R, float(l), margin=_PSD_EPS) if l < 0 else (R, 0.0, l)
            out.append((R2, idx))
            cs.append(c)
        cs = np.asarray(cs)
        entry.update(ridge_blocks=int(np.sum(cs > 0)), ridge_max=float(cs.max()),
                     ridge_mean_nonzero=float(cs[cs > 0].mean()) if np.any(cs > 0) else 0.0)
    entry["spectrum_s"] = time.perf_counter() - t0
    rec.append(entry)
    if mode == "check" and method == "ldpred3-lassosum2" and ctx is not None:
        s_min = min(ctx)
        if spec["min"] < -s_min / (1.0 - s_min):
            warnings.warn(f"an LD block has lambda_min {spec['min']:.3f} below -s/(1-s) = "
                          f"{-s_min / (1 - s_min):.3f} of lassosum2's smallest s = {s_min}; "
                          "its (1-s)R + sI is indefinite (stage1_kwargs ld_psd='floor' raises s)",
                          RuntimeWarning, stacklevel=4)
    return out


def _lassosum2_s_floor(s_seq, lam_min: float):
    """The lassosum2 s grid raised so (1 - s) R + s I is PSD for every block:
    s >= c / (1 + c), c = -lam_min + 1e-3 (unchanged when lam_min >= 0)."""
    if lam_min >= 0:
        return tuple(s_seq), 0.0
    c = -float(lam_min) + _PSD_EPS
    s_min = c / (1.0 + c)
    return tuple(sorted({float(max(s, s_min)) for s in s_seq})), s_min


def ldpred3_stage1_scores(st, ys: np.ndarray, h2: float, folds, held,
                          *, method: str, keep_prediction: bool = False,
                          block_size: int = 500, refit: bool = False,
                          **fit_kwargs) -> dict:
    """Contract dictionary from cross-fitted LDpred3 LOCO scores.

    ``method`` is ``"ldpred3-inf"`` (the infinitesimal closed form at the
    genome-wide ``h2`` of the variance fit), ``"ldpred3-lassosum2"``
    (penalty selected per fold by pseudo-validation on that fold's own rows)
    or ``"ldpred3-fast"`` (experimental LDpred3-auto on the fast layout:
    int8 dense sub-blocks plus an int8 Nystrom landmark factor, within
    group-pure blocks of ``fast_block_size`` variants, so the LOCO
    subtraction stays exact given the shared h2 and p; a seeded Gibbs
    chain per fold; warns once and records ``ldpred3_version`` and
    ``experimental=True`` in ``stage1_info``).
    ``fit_kwargs`` (with ``block_size``) pass to lassosum2, except the LD
    options: ``ld`` (``"f32"``, the default float32 blocks, or an int8 /
    hard-call representation of :mod:`mixmogam._ld8`), ``ld_cache``
    (reuse a cover across traits) and ``ld_tile``; ``ldpred3-fast`` takes
    ``width``, ``landmarks``, ``reach``, ``burn_in``, ``num_iter``,
    ``seed``, ``p_init``, ``layout_seed`` and ``fast_block_size``. The
    ``ld_split`` (``"all"`` | ``"train"``) and ``ld_psd`` (``"off"`` |
    ``"check"`` | ``"floor"`` | ``"ridge"``) are described in the module notes. The
    hold-out MSE of a fit excluding ``held`` is the candidate's
    model-selection score, on the same samples as the elastic-net
    candidates'. With ``refit`` (the relatives / structure auto-refit
    for lassosum2 and fast), each fold is re-fit k x groups times,
    leaving out the tested group's blocks and re-selecting that method's
    hyperparameters, instead of subtracting the group's share of a
    genome-wide fit. ``ldpred3-inf`` does not need this when h2 is
    shared. ``stage1_info`` in the result records the LD build or
    cache hit, its bytes and the fit times. The returned dict includes
    ``refit``, the flag that actually ran.
    """
    import time

    lg = st.lg
    k = 0 if folds is None else int(np.max(np.asarray(folds))) + 1
    if not k:
        raise ValueError("stage1 'ldpred3-*' needs cross-fitted scores (loco_folds >= 2)")
    folds = np.asarray(folds, dtype=np.int64)
    fit_kwargs = dict(fit_kwargs)
    ld = fit_kwargs.pop("ld", "f32")
    ld_cache = fit_kwargs.pop("ld_cache", None)
    ld_tile = fit_kwargs.pop("ld_tile", None)
    ld_split = fit_kwargs.pop("ld_split", "all")
    ld_psd = fit_kwargs.pop("ld_psd", None)
    from mixmogam._ld8 import LD_KINDS, _validate_block_size
    if ld not in LD_KINDS:
        raise ValueError(f"stage1_kwargs['ld'] must be one of {list(LD_KINDS)}, got {ld!r}")
    if ld_split not in ("all", "train"):
        raise ValueError(f"stage1_kwargs['ld_split'] must be 'all' or 'train', got {ld_split!r}")
    if ld_psd is None:
        ld_psd = "check" if ld.startswith("d8") else "off"
    if ld_psd not in _PSD_MODES:
        raise ValueError(f"stage1_kwargs['ld_psd'] must be one of {list(_PSD_MODES)}, "
                         f"got {ld_psd!r}")
    fast_opt = dict(_FAST_DEFAULTS)
    if method == "ldpred3-fast":
        for key_ in _FAST_KEYS:
            if key_ in fit_kwargs:
                fast_opt[key_] = fit_kwargs.pop(key_)
        if fit_kwargs:
            raise ValueError(f"stage1 'ldpred3-fast' got unknown keywords {sorted(fit_kwargs)}")
        _validate_block_size(fast_opt["fast_block_size"], "fast_block_size")
    else:
        _validate_block_size(block_size)
    if h2 <= 0:  # no polygenic variance: no scores, as in the elastic-net path
        held_mask = _held_mask(lg.n, held)
        ys = np.asarray(ys, dtype=np.float64).ravel()
        G = lg.n_groups
        out = {"resid": np.repeat(ys[:, None], G, axis=1),
               "offset_slope": np.zeros(G),
               "mse": np.array([float(np.mean(ys[held_mask] ** 2))]),
               "best": 0, "cv_converged": True, "loco_iterations": 0,
               "loco_converged": True, "cv_grid": [method], "cv_best": method,
               "refit": False}
        if keep_prediction:
            out["prediction"] = np.zeros((lg.n, G))
        return out
    ldpred3 = _import_ldpred3()
    # FastLD stores int8 (dense sub-blocks + int8 Nystrom factor) whatever
    # the dense source it is built from.
    info = {"ld": "fast-int8" if method == "ldpred3-fast" else ld, "ld_source": ld,
            "method": method, "ld_split": ld_split, "ld_psd": ld_psd,
            "ldpred3_version": getattr(ldpred3, "__version__", None)}
    if method == "ldpred3-fast":
        import warnings
        info["experimental"] = True
        try:
            from ldpred3.methods import FAST_EXPERIMENTAL_WARNING as _fast_msg
        except ImportError:  # pragma: no cover - ldpred3 < 0.10.4
            _fast_msg = ("LDpred3-fast is experimental: int8 rounding can make the "
                         "working LD indefinite even when the source is positive "
                         "semidefinite. The independent-Gibbs update does not repair "
                         "this. Prediction and architecture estimates require "
                         "validation against native-LD fits.")
        warnings.warn(
            f"stage1 'ldpred3-fast' is experimental (ldpred3 {ldpred3.__version__}): {_fast_msg}",
            UserWarning, stacklevel=3)
    t0 = time.perf_counter()
    marg = fold_marginals(lg, ys, folds, held)
    info["marginals_s"] = time.perf_counter() - t0
    order = np.argsort(lg.groups, kind="stable")
    # Column indices of MarginalFits: k fold training complements, the
    # hold-out complement, then the matching validation sets.
    n_sets = len(marg.labels)
    train_cols = list(range(k + 1))
    val_cols = list(range(k + 1, n_sets))
    scale = marg.scale[order]
    split = None
    if ld_split == "train":
        from mixmogam._ld8 import split_parts
        split = split_parts(lg.n, folds, _held_mask(lg.n, held))
    if method == "ldpred3-fast":
        fast_opt["_psd"] = ld_psd if ld_psd in ("floor", "ridge") else "off"
    covers = _stage1_blocks(lg, order, marg.scale, method, block_size, ld, ld_cache, ld_tile,
                            fast_opt, ldpred3, info, split)
    fast_opt.pop("_psd", None)
    s_default = None
    if method == "ldpred3-lassosum2":
        import inspect
        s_default = tuple(fit_kwargs.get(
            "s_seq", inspect.signature(ldpred3.lassosum2).parameters["s_seq"].default))
    if method != "ldpred3-fast":
        covers = [_psd_cover(c, ld_psd, info, method, s_default) for c in covers]
        if method == "ldpred3-lassosum2" and ld_psd == "floor":
            lam_min = min(e["min"] for e in info["ld_spectrum"])
            fit_kwargs["s_seq"], info["s_floor"] = _lassosum2_s_floor(s_default, lam_min)
            info["s_seq"] = list(fit_kwargs["s_seq"])

    c0, t0 = time.process_time(), time.perf_counter()
    fast_diag = [] if method == "ldpred3-fast" else None
    effects = np.zeros((lg.m, len(train_cols)))
    for j, c in enumerate(train_cols):
        blocks = covers[0] if len(covers) == 1 else covers[j]
        beta_std, _ = ldpred3.standardize_betas(
            marg.beta[:, c], marg.se[:, c], float(marg.n_eff[c]))
        kwargs = dict(fit_kwargs)
        if method == "ldpred3-lassosum2":
            val_std, _ = ldpred3.standardize_betas(
                marg.beta[:, val_cols[j]], marg.se[:, val_cols[j]],
                float(marg.n_eff[val_cols[j]]))
            kwargs["validation_beta"] = val_std[order]
        if method == "ldpred3-fast":
            kwargs = {**fast_opt, "seed": int(fast_opt["seed"]) + j}
        e = _fit_effects(ldpred3, blocks, beta_std[order], float(marg.n_eff[c]),
                         {"ldpred3-inf": "inf", "ldpred3-lassosum2": "lassosum2",
                          "ldpred3-fast": "fast"}[method], h2, kwargs, fast_diag)
        with np.errstate(divide="ignore", invalid="ignore"):
            e = np.where(scale > 0, e / np.where(scale > 0, scale, 1.0), 0.0)
        effects[order, j] = e  # back to variant index order
    info.update(solve_s=time.perf_counter() - t0, solve_cpu_s=time.process_time() - c0,
                n_fits=len(train_cols))
    if fast_diag is not None:
        info["fast_chains"] = fast_diag
        rh = [d["rhat_h2"] for d in fast_diag if d.get("rhat_h2") is not None]
        info["fast_max_rhat_h2"] = float(np.nanmax(rh)) if rh and np.any(np.isfinite(rh)) else None
        if info["fast_max_rhat_h2"] is not None and info["fast_max_rhat_h2"] > 1.1:
            import warnings
            warnings.warn(f"ldpred3-fast: h2 chain split R-hat {info['fast_max_rhat_h2']:.2f} > 1.1; "
                          "raise burn_in / num_iter", RuntimeWarning, stacklevel=3)

    # Held-out candidate score: the fit excluding the hold-out, genome-wide.
    from mixmogam.twostep import _out_of_fold_loco  # late: the patchable seam
    t0 = time.perf_counter()
    held_mask = _held_mask(lg.n, held)
    fitted = lg.project(_fitted_unprojected(lg, effects[:, k]))
    mse = np.array([float(np.mean((ys[held_mask] - fitted[held_mask]) ** 2))])

    used_refit = bool(refit and h2 > 0)
    if used_refit:
        prediction, n_leave = _ldpred3_leave_group_scores(
            ldpred3, lg, covers, order, scale, marg, train_cols, folds, method, h2,
            fit_kwargs, fast_opt, fast_diag)
        info["n_fits"] = int(info.get("n_fits", 0)) + n_leave
    else:
        prediction = _out_of_fold_loco(lg, effects[:, :k], folds)
    info["predict_s"] = time.perf_counter() - t0
    out = {"mse": mse, "best": 0, "cv_converged": True,
           "loco_iterations": 0, "loco_converged": True,
           "cv_grid": [method], "cv_best": method, "stage1_info": info,
           "refit": used_refit}
    out.update(residualize_score(lg, ys, prediction))
    if keep_prediction:
        out["prediction"] = prediction
    return out



def _ldpred3_leave_group_scores(ldpred3, lg, covers, order, scale, marg, train_cols, folds,
                                method, h2, fit_kwargs, fast_opt, fast_diag):
    """k x groups leave-group-out LDpred3 scores (``refit=True``).

    Each block is group-pure. Dropping group ``g``'s blocks and re-fitting
    (re-selecting lassosum2 / fast hyperparameters on the rest) is the
    analogue of the elastic-net k x groups refit. Block indices are
    recompacted: LDpred3 requires a cover that tiles 0..m_keep-1.
    """
    short = {"ldpred3-inf": "inf", "ldpred3-lassosum2": "lassosum2",
             "ldpred3-fast": "fast"}[method]
    k = int(np.max(np.asarray(folds))) + 1
    G = lg.n_groups
    order = np.asarray(order, dtype=np.int64)
    layout = covers[0]
    block_group = []
    for _, idx in layout:
        idx = np.asarray(idx).ravel()
        block_group.append(int(lg.groups[order[idx[0]]]))
    block_group = np.asarray(block_group, dtype=np.int64)
    prediction = np.empty((lg.n, G))
    n_leave = 0
    n_sets = len(marg.labels)
    val_cols = list(range(k + 1, n_sets))
    for j in range(k):
        blocks = covers[0] if len(covers) == 1 else covers[j]
        beta_std, _ = ldpred3.standardize_betas(
            marg.beta[:, train_cols[j]], marg.se[:, train_cols[j]],
            float(marg.n_eff[train_cols[j]]))
        beta_ord = beta_std[order]
        if short == "lassosum2":
            val_std, _ = ldpred3.standardize_betas(
                marg.beta[:, val_cols[j]], marg.se[:, val_cols[j]],
                float(marg.n_eff[val_cols[j]]))
            val_ord = val_std[order]
        rows = np.flatnonzero(np.asarray(folds) == j)
        for g in range(G):
            keep = [i for i, bg in enumerate(block_group) if bg != g]
            if not keep:
                prediction[rows, g] = 0.0
                continue
            sub, cat = [], []
            offset = 0
            for i in keep:
                R, idx = blocks[i]
                idx = np.asarray(idx).ravel()
                sub.append((R, np.arange(offset, offset + idx.size)))
                cat.append(idx)
                offset += idx.size
            cat = np.concatenate(cat)
            kwargs = dict(fit_kwargs)
            if short == "lassosum2":
                kwargs["validation_beta"] = val_ord[cat]
            if short == "fast":
                kwargs = {**fast_opt, "seed": int(fast_opt["seed"]) + j * (G + 1) + 1 + g}
            e = _fit_effects(ldpred3, sub, beta_ord[cat], float(marg.n_eff[train_cols[j]]),
                             short, h2, kwargs, fast_diag)
            n_leave += 1
            with np.errstate(divide="ignore", invalid="ignore"):
                e = np.where(scale[cat] > 0, e / np.where(scale[cat] > 0, scale[cat], 1.0), 0.0)
            beta = np.zeros(lg.m)
            beta[order[cat]] = e
            # Unprojected z'beta, same units as _out_of_fold_loco.
            prediction[rows, g] = _fitted_unprojected(lg, beta)[rows]
    return prediction, n_leave


def _fitted_unprojected(lg, beta: np.ndarray) -> np.ndarray:
    """z' beta over one pass of the unprojected rows (the projection of the
    result is the projected fitted value)."""
    out = np.zeros(lg.n)
    for take, g, z in lg.raw_slices():
        out += z.T @ beta[take]
    return out
