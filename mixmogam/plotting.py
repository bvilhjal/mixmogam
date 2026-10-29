"""Manhattan and QQ plots (matplotlib optional)."""

from __future__ import annotations

from typing import TYPE_CHECKING, Sequence

import numpy as np
from scipy import stats

if TYPE_CHECKING:
    from matplotlib.axes import Axes

    from mixmogam.results import GwasResult

__all__ = ["plot_manhattan", "plot_qq", "qq_quantiles"]

_POINT_COLORS = ("#08519c", "#6baed6")
_ACCENT = "#cb181d"
_LINE = "#252525"
_CHI2_1DF_MEDIAN = float(stats.chi2.ppf(0.5, 1))  # 0.4549


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


def _style(ax: Axes) -> None:
    """Despine, light y-grid behind the points, smaller tick labels."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", color="#e6e6e6", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=9)


def _position_scale(pos: np.ndarray) -> tuple[float, str]:
    """Scale positions to Mb/kb when the panel is large enough to need it."""
    top = float(pos.max()) if pos.size else 0.0
    if top >= 1e6:
        return 1e6, "Mb"
    if top >= 1e3:
        return 1e3, "kb"
    return 1.0, "bp"


def _genomic_control(p: np.ndarray) -> float:
    """lambda_GC from p-values (same formula as GwasResult.genomic_control)."""
    p = np.asarray(p, dtype=np.float64)
    p = p[np.isfinite(p)]
    if p.size == 0:
        return float("nan")
    chi2 = stats.chi2.isf(np.clip(p, 1e-300, 1.0), 1)
    return float(np.median(chi2) / _CHI2_1DF_MEDIAN)


def _line_label(ax: Axes, y: float, text: str) -> None:
    """Small right-edge label for a horizontal threshold line."""
    from matplotlib.transforms import blended_transform_factory

    trans = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(0.995, y, f" {text}", transform=trans, ha="right", va="bottom",
            fontsize=7.5, color="#555555")


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
    ax: Axes | None = None,
    title: str | None = None,
    savepath: str | None = None,
    band: bool = True,
    show_lambda: bool = True,
) -> Axes:
    """Log-scale QQ plot of p-values with the uniform reference line.

    With ``band`` (default) a shaded 95% pointwise band from the Beta
    order-statistic distribution is drawn behind the points; with
    ``show_lambda`` (default) the genomic inflation factor is annotated.
    """
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5), layout="constrained")
    exp, obs = qq_quantiles(p)
    n = exp.size
    if band and n:
        i = np.arange(1, n + 1)
        lo = -np.log10(stats.beta.ppf(0.975, i, n - i + 1))
        hi = -np.log10(stats.beta.ppf(0.025, i, n - i + 1))
        ax.fill_between(exp, lo, hi, color=_POINT_COLORS[1], alpha=0.25, lw=0, zorder=0)
    ax.scatter(exp, obs, s=6, color=_POINT_COLORS[0], alpha=0.75, lw=0, zorder=2)
    top = max(float(exp.max()) if n else 0.0, float(obs.max()) if n else 0.0)
    lim = top * 1.05
    ax.plot([0, lim], [0, lim], color=_ACCENT, lw=1.1, zorder=1)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel(r"Expected $-\log_{10}(p)$", fontsize=11)
    ax.set_ylabel(r"Observed $-\log_{10}(p)$", fontsize=11)
    if show_lambda:
        lam = _genomic_control(np.asarray(p, dtype=np.float64))
        ax.text(0.04, 0.96, rf"$\lambda_{{GC}}$ = {lam:.2f}", transform=ax.transAxes,
                va="top", ha="left", fontsize=10, color="#333333")
    if title:
        ax.set_title(title, fontsize=12)
    _style(ax)
    if savepath:
        ax.figure.savefig(savepath, dpi=150, bbox_inches="tight")
    return ax


def plot_manhattan(
    result: GwasResult,
    ax: Axes | None = None,
    bonferroni: bool = True,
    alpha: float = 0.05,
    suggestive: float | None = 1e-5,
    label_top: int = 5,
    highlight: Sequence[int] | None = None,
    title: str | None = None,
    savepath: str | None = None,
) -> Axes:
    """Manhattan plot from a :class:`mixmogam.results.GwasResult`.

    Chromosomes alternate two colors on a gap-separated axis, with
    positions scaled to Mb/kb for large panels. With ``bonferroni``
    (default) a dashed line marks ``alpha`` / n variants; ``suggestive``
    adds a dotted line (default 1e-5, ``None`` disables). ``label_top``
    annotates the smallest-p variants (variant ID, or chr:pos) with
    staggered leader lines; ``highlight`` indexes variants to circle.
    """
    if not 0 < alpha < 1:
        raise ValueError("alpha must be in (0, 1)")
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(11, 4.5), layout="constrained")
    chrom = np.asarray(result.chromosome)
    pos = np.asarray(result.position, dtype=np.float64)
    scale, unit = _position_scale(pos)
    pos = pos / scale
    nlp = -np.log10(np.clip(result.p, 1e-300, 1.0))

    # Vectorized chromosome layout: group positions by chromosome, then lay
    # the groups out left to right with a 2% span gap.
    chroms, inv = np.unique(chrom, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    pos_sorted = pos[order]
    counts = np.bincount(inv, minlength=chroms.size)
    starts = np.concatenate(([0], np.cumsum(counts)[:-1]))
    cmin = pos_sorted[starts]
    cmax = np.maximum.reduceat(pos_sorted, starts)
    span = np.maximum(cmax - cmin, 1.0)
    offsets = np.concatenate(([0.0], np.cumsum(span * 1.02)))[:-1]
    x = pos + offsets[inv] - cmin[inv]

    for i in range(chroms.size):
        sel = inv == i
        ax.scatter(x[sel], nlp[sel], s=3, color=_POINT_COLORS[i % 2], alpha=0.8, lw=0, zorder=2)

    ymax = float(nlp.max()) if nlp.size else 0.0
    if bonferroni and nlp.size:
        y_b = -np.log10(alpha / nlp.size)
        ax.axhline(y_b, color=_ACCENT, ls="--", lw=1.0, zorder=1)
        _line_label(ax, y_b, "genome-wide")
        ymax = max(ymax, y_b)
    if suggestive is not None and suggestive > 0:
        y_s = -np.log10(suggestive)
        ax.axhline(y_s, color=_LINE, ls=":", lw=0.9, zorder=1)
        _line_label(ax, y_s, "suggestive")
    if highlight is not None:
        idx = np.asarray(highlight, dtype=int)
        ax.scatter(x[idx], nlp[idx], s=28, facecolors="none", edgecolors=_ACCENT, lw=1.2, zorder=3)
    if label_top and nlp.size:
        top = np.argsort(np.asarray(result.p, dtype=np.float64))[:label_top]
        vids = result.variant_ids
        raw_pos = np.asarray(result.position)
        for rank, j in enumerate(top):
            label = str(vids[j]) if vids is not None else f"chr{chrom[j]}:{int(raw_pos[j])}"
            if len(label) > 24:
                label = label[:21] + "..."
            ax.annotate(label, (x[j], nlp[j]), xytext=(0, 7 + (rank % 3) * 10),
                        textcoords="offset points", ha="center", fontsize=7, color="#222222",
                        arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.6), zorder=4)

    pad = (x.max() - x.min()) * 0.01 if x.size else 1.0
    ax.set_xlim(x.min() - pad, x.max() + pad)
    ax.set_ylim(0, ymax * 1.15 + 0.5)
    ax.set_xticks(offsets + (cmin + cmax) / 2)
    ax.set_xticklabels([str(c) for c in chroms])
    ax.set_xlabel(f"Chromosome (position in {unit})" if scale != 1.0 else "Chromosome",
                  fontsize=11)
    ax.set_ylabel(r"$-\log_{10}(p)$", fontsize=11)
    if title:
        ax.set_title(title, fontsize=12)
    _style(ax)
    if savepath:
        ax.figure.savefig(savepath, dpi=150, bbox_inches="tight")
    return ax
