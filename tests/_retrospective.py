"""Dense oracle shared by the retrospective-test suites."""

import numpy as np


def dense_probabilities(g, Xd):
    """Allele-frequency regression using an independent dense least-squares fit."""
    n, q = Xd.shape
    mean = g.mean(axis=0)
    if q == 1:
        return np.broadcast_to(mean / 2, g.shape)
    fitted = Xd @ np.linalg.lstsq(Xd, g, rcond=None)[0] - mean
    r2 = np.sum(fitted**2, axis=0) / np.sum((g - mean) ** 2, axis=0)
    F = (r2 / (q - 1)) / ((1 - r2) / (n - q))
    return np.clip((mean + np.clip(1 - 1 / F, 0, 1) * fitted) / 2, 0.5 / n, 1 - 0.5 / n)


def dense_rho(g, Xd, a):
    """The covariate-specific genotype variance factor of columns ``g``
    (genotype counts) for residuals ``a``: allele frequencies fitted on
    ``Xd``, their deviations shrunk by (F - 1) / F, binomial variances
    weighted by a^2 against their plain mean."""
    p = dense_probabilities(g, Xd)
    pi = p * (1 - p)
    return (a * a) @ pi / (np.sum(a * a) * pi.mean(axis=0))
