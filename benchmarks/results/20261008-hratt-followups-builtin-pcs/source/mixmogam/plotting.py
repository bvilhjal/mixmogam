"""Manhattan and QQ plots (matplotlib optional)."""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

__all__ = ["plot_manhattan", "plot_qq", "qq_quantiles"]


def _plt():
    try:
        import matplotlib

        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt

        return plt
    except ImportError as err:  # pragma: no cover
        raise ImportError(
            "plotting requires matplotlib; install with "
            'pip install "mixmogam[plot]"'
        ) from err


def qq_quantiles(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Expected vs observed -log10(p) quantiles (v1 analyze_gwas_results)."""
    p = np.asarray(p, dtype=np.float64)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    p = np.sort(p)
    n = p.size
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(p)
    return expected, observed


def plot_qq(
    p: np.ndarray,
    ax=None,
    title: Optional[str] = None,
    savepath: Optional[str] = None,
):
    """Log-scale QQ plot of p-values with the uniform reference line."""
    plt = _plt()
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    exp, obs = qq_quantiles(p)
    ax.scatter(exp, obs, s=4, color="#2166ac", alpha=0.7)
    lim = max(exp.max(), obs.max()) * 1.05
    ax.plot([0, lim], [0, lim], color="#b2182b", lw=1)
    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    if title:
        ax.set_title(title)
    if savepath:
        ax.figure.savefig(savepath, dpi=150, bbox_inches="tight")
    return ax


def plot_manhattan(
    result,
    ax=None,
    bonferroni: bool = True,
    highlight: Optional[Sequence] = None,
    title: Optional[str] = None,
    savepath: Optional[str] = None,
):
    """Manhattan plot from a :class:`mixmogam.results.GwasResult`.

    Chromosomes alternate colors; ``highlight`` indexes variants to mark.
    """
    plt = _plt()
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))
    chrom = np.asarray(result.chromosome)
    pos = np.asarray(result.position, dtype=np.int64)
    nlp = -np.log10(np.clip(result.p, 1e-300, 1.0))

    chroms = np.unique(chrom)
    offsets = {}
    x = np.empty(nlp.size)
    prev_end = 0.0
    for c in chroms:
        sel = chrom == c
        span = pos[sel].max() - pos[sel].min() if sel.sum() > 1 else 1
        offsets[c] = prev_end - pos[sel].min()
        x[sel] = pos[sel] + offsets[c]
        prev_end += span * 1.02
    colors = ["#4393c3", "#f4a582", "#92c5de", "#d6604d", "#2166ac", "#b2182b"]
    for i, c in enumerate(chroms):
        sel = chrom == c
        ax.scatter(
            x[sel], nlp[sel], s=2, color=colors[i % len(colors)], alpha=0.7, lw=0
        )
    if bonferroni and nlp.size:
        ax.axhline(-np.log10(0.05 / nlp.size), color="#444444", ls="--", lw=0.8)
    if highlight is not None:
        idx = np.asarray(highlight, dtype=int)
        ax.scatter(x[idx], nlp[idx], s=14, facecolors="none", edgecolors="#b2182b")
    ticks = [offsets[c] + (pos[chrom == c].max() + pos[chrom == c].min()) / 2 for c in chroms]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(c) for c in chroms])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log10(p)")
    if title:
        ax.set_title(title)
    if savepath:
        ax.figure.savefig(savepath, dpi=150, bbox_inches="tight")
    return ax
