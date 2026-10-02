"""Reference oracle: float64 port of the mixmogam v1 EMMA/EMMAX math.

This is a deliberately naive, per-SNP implementation mirroring
``linear_models.py`` from the v1.0-legacy tag (pinv-based hat matrix,
S(K+I)S eigenspace, log-delta grid + Newton, per-SNP least squares with a
residualized transformed SNP). It exists to certify that the v2 batched
engine computes the same statistics.
"""

import numpy as np
from scipy import linalg, optimize, stats


def _eigen_L(K):
    evals, evecs = linalg.eigh(np.asarray(K, dtype=np.float64))
    return evals, evecs


def _eigen_R(X, K, n):
    X = np.asarray(X, dtype=np.float64)
    q = X.shape[1]
    X_pinv = np.linalg.pinv(X.T @ X)
    hat = X @ X_pinv @ X.T
    S = np.eye(n) - hat
    M = S @ (np.asarray(K, dtype=np.float64) + np.eye(n)) @ S
    evals, evecs = linalg.eigh(M)
    return evals[q:] - 1.0, evecs[:, q:].T  # vectors as rows, like v1


def _rell(delta, lam, sq_etas):
    p = lam.size
    v = lam + delta
    return 0.5 * (
        p * (np.log(p / (2.0 * np.pi)) - 1.0 - np.log(np.sum(sq_etas / v)))
        - np.sum(np.log(v))
    )


def _redll(delta, lam, sq_etas):
    p = lam.size
    v = lam + delta
    v1 = sq_etas / v
    return p * np.sum(v1 / v) / np.sum(v1) - np.sum(1.0 / v)


def emma_fit(y, K, X=None, ngrids=100, llim=-10.0, ulim=10.0):
    """REML fit returning (delta, vg, ve, ll, pseudo_h2)."""
    y = np.asarray(y, dtype=np.float64).ravel()
    n = y.size
    if X is None:
        X = np.ones((n, 1))
    lam_R, U_R = _eigen_R(X, K, n)
    etas = U_R @ y
    sq_etas = etas * etas

    log_deltas = np.linspace(llim, ulim, ngrids + 1)
    deltas = np.exp(log_deltas)
    lamm = lam_R[:, None] + deltas[None, :]
    s1 = np.sum(sq_etas[:, None] / lamm, axis=0)
    s2 = np.sum(np.log(lamm), axis=0)
    p = lam_R.size
    lls = 0.5 * (p * (np.log(p / (2.0 * np.pi)) - 1.0 - np.log(s1)) - s2)
    s3 = np.sum(sq_etas[:, None] / (lamm * lamm), axis=0)
    s4 = np.sum(1.0 / lamm, axis=0)
    dlls = 0.5 * (p * s3 / s1 - s4)

    max_ll_i = int(np.argmax(lls))
    crosses = np.nonzero((dlls[1:] < 0.0) & (dlls[:-1] > 0.0))[0]
    if crosses.size == 0:
        opt_delta = float(deltas[max_ll_i])
        opt_ll = float(lls[max_ll_i])
    else:
        i = int(crosses[np.argmax((lls[1:] + lls[:-1])[crosses] * 0.5)])
        opt_delta = 0.5 * (deltas[i] + deltas[i + 1])
        try:
            new_delta = optimize.newton(
                _redll, opt_delta, args=(lam_R, sq_etas), tol=1e-6, maxiter=100
            )
        except RuntimeError:
            new_delta = opt_delta
        if (
            deltas[i] - 1e-6 <= new_delta <= deltas[i + 1] + 1e-6
            or (i == 0 and 0.0 < new_delta <= deltas[1] + 1e-6)
            or (i == len(deltas) - 2 and new_delta >= deltas[i] - 1e-6)
        ):
            opt_delta = new_delta
        opt_ll = _rell(opt_delta, lam_R, sq_etas)
        if opt_ll < lls[max_ll_i]:
            opt_delta = float(deltas[max_ll_i])
            opt_ll = float(lls[max_ll_i])
    vg = float(np.sum(sq_etas / (lam_R + opt_delta)) / p)
    return opt_delta, vg, vg * opt_delta, opt_ll, 1.0 / (1.0 + opt_delta)


def emmax_scan(snps, y, K, X=None, delta=None):
    """Per-SNP EMMAX F-test exactly as v1's ``_emmax_f_test_`` did.

    Returns (ps, f_stats, delta).
    """
    snps = np.asarray(snps, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64).ravel()
    n = y.size
    if X is None:
        X = np.ones((n, 1))
    q = X.shape[1]
    if delta is None:
        delta, _, _, _, _ = emma_fit(y, K, X)

    evals_L, evecs_L = _eigen_L(K)
    H_sqrt_inv = (np.diag(1.0 / np.sqrt(evals_L + delta))) @ evecs_L.T

    h0_X = H_sqrt_inv @ X
    Y = H_sqrt_inv @ y
    h0_betas, h0_rss, _, _ = linalg.lstsq(h0_X, Y)
    Y = Y - h0_X @ h0_betas

    Q, _ = linalg.qr(h0_X, mode="economic")
    M = H_sqrt_inv.T @ (np.eye(n) - Q @ Q.T)

    num_snps = snps.shape[0]
    rss_list = np.full(num_snps, h0_rss)
    for j in range(num_snps):
        xj = snps[j] @ M  # transformed, residualized SNP (row)
        betas, rss, rank, _ = linalg.lstsq(xj[:, None], Y)
        if rss:
            rss_list[j] = float(np.squeeze(rss))

    n_p = n - (q + 1)
    rss_ratio = h0_rss / rss_list
    f_stats = (rss_ratio - 1.0) * n_p
    ps = stats.f.sf(f_stats, 1, n_p)
    return ps, f_stats, delta
