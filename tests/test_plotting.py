"""Plotting API tests: top-level exports, GwasResult methods, the plot CLI."""

import numpy as np
import pytest

import mixmogam
from mixmogam import Genotypes, gwas
from mixmogam.cli import main
from mixmogam.results import GwasResult
from mixmogam.simulate import simulate_genotypes, simulate_traits

mpl = pytest.importorskip("matplotlib")
mpl.use("Agg")


@pytest.fixture(scope="module")
def result():
    G = simulate_genotypes(n=300, m=1200, n_pop=3, pop_fst=0.3, seed=41)
    gt = Genotypes(G.T, chromosome=np.repeat([1, 2, 3], 400),
                   position=np.tile(np.arange(400) * 1000, 3))
    y = simulate_traits(G, h2=0.5, n_causal=4, seed=42)["y"]
    return gwas(y, gt)


def test_top_level_exports():
    assert mixmogam.plot_manhattan is not None
    assert mixmogam.plot_qq is not None
    from mixmogam.plotting import plot_manhattan, plot_qq

    assert mixmogam.plot_manhattan is plot_manhattan
    assert mixmogam.plot_qq is plot_qq
    assert "plot_manhattan" in dir(mixmogam)


def test_qq_quantiles():
    from mixmogam.plotting import qq_quantiles

    rng = np.random.default_rng(0)
    exp, obs = qq_quantiles(rng.uniform(1e-8, 1.0, 500))
    assert exp.size == obs.size == 500
    assert np.all(np.diff(exp) < 0)  # expected -log10 quantiles decrease with i
    assert np.isfinite(exp).all() and np.isfinite(obs).all()
    # non-finite and out-of-range p-values are dropped
    exp2, obs2 = qq_quantiles(np.array([0.5, np.nan, 0.0, -1.0, 2.0, 0.25]))
    assert exp2.size == 2


def test_result_plot_methods(result, tmp_path):
    man = tmp_path / "manhattan.png"
    qq = tmp_path / "qq.png"
    ax = result.plot_manhattan(savepath=str(man), title="sim")
    assert ax is not None and man.stat().st_size > 1000
    ax = result.plot_qq(savepath=str(qq), title="sim")
    assert ax is not None and qq.stat().st_size > 1000


def test_plot_manhattan_alpha_and_highlight(result, tmp_path):
    from mixmogam.plotting import plot_manhattan

    out = tmp_path / "m.png"
    plot_manhattan(result, alpha=0.01, highlight=[0, 10], savepath=str(out))
    assert out.stat().st_size > 1000
    with pytest.raises(ValueError, match="alpha"):
        plot_manhattan(result, alpha=1.5)


def test_qq_confidence_band_toggle(result):
    from matplotlib.collections import PolyCollection

    from mixmogam.plotting import plot_qq

    ax = plot_qq(result.p)
    assert any(isinstance(c, PolyCollection) for c in ax.collections)
    ax = plot_qq(result.p, band=False)
    assert not any(isinstance(c, PolyCollection) for c in ax.collections)


def test_manhattan_labels_and_threshold_lines(result):
    from mixmogam.plotting import plot_manhattan

    ax = plot_manhattan(result, label_top=3, suggestive=1e-5)
    # 3 SNP annotations + "genome-wide" and "suggestive" line labels
    assert len(ax.texts) == 5
    ax = plot_manhattan(result, label_top=0, suggestive=None, bonferroni=False)
    assert len(ax.texts) == 0


def test_qq_lambda_annotation(result):
    from mixmogam.plotting import plot_qq

    ax = plot_qq(result.p)
    assert any("lambda" in t.get_text() for t in ax.texts)
    ax = plot_qq(result.p, show_lambda=False)
    assert not ax.texts


def test_manhattan_position_scaling(result):
    """Large positions are scaled to Mb and the axis says so."""
    from mixmogam.plotting import plot_manhattan

    big = result.take(np.arange(0, len(result), 10))
    big.position = big.position * 1000  # now in the Mb range
    ax = plot_manhattan(big, label_top=0, bonferroni=False, suggestive=None)
    assert "Mb" in ax.get_xlabel()


def test_manhattan_chromosome_blocks_do_not_overlap(result):
    """The vectorized layout keeps each chromosome's x-range disjoint."""
    from mixmogam.plotting import plot_manhattan

    ax = plot_manhattan(result, bonferroni=False)
    offsets = np.concatenate([c.get_offsets()[:, 0] for c in ax.collections])
    chrom = np.asarray(result.chromosome)
    for c in np.unique(chrom):
        x = offsets[chrom == c]
        assert x.min() >= 0
    # blocks are ordered by chromosome and separated by gaps
    edges = [(offsets[chrom == c].min(), offsets[chrom == c].max()) for c in np.unique(chrom)]
    assert all(lo_next > hi_prev for (_, hi_prev), (lo_next, _) in zip(edges, edges[1:]))


def test_plot_cli_writes_both_pngs(result, tmp_path, capsys):
    csv = tmp_path / "gwas.csv"
    result.write_csv(str(csv))
    prefix = str(tmp_path / "fig")
    rc = main(["plot", "--in", str(csv), "--out-prefix", prefix, "--title", "cli"])
    assert rc == 0
    assert (tmp_path / "fig_manhattan.png").stat().st_size > 1000
    assert (tmp_path / "fig_qq.png").stat().st_size > 1000
    assert "lambda_GC" in capsys.readouterr().out


def test_plot_cli_no_threshold_and_gzip(result, tmp_path):
    csv = tmp_path / "gwas.tsv.gz"
    result.write_csv(str(csv))
    prefix = str(tmp_path / "g")
    rc = main(["plot", "-i", str(csv), "-o", prefix, "--no-threshold", "--alpha", "0.01"])
    assert rc == 0
    assert (tmp_path / "g_manhattan.png").exists()


def test_plot_cli_missing_file_is_actionable(tmp_path, capsys):
    rc = main(["plot", "--in", str(tmp_path / "nope.csv"), "-o", str(tmp_path / "x")])
    assert rc == 2
    assert "error" in capsys.readouterr().err


def test_plot_cli_roundtrip_matches_api(result, tmp_path):
    """The CLI plots the same values the API would read back from the CSV."""
    csv = tmp_path / "gwas.csv"
    result.write_csv(str(csv))
    back = GwasResult.read_csv(str(csv))
    assert np.array_equal(back.p, result.p)
    assert np.array_equal(back.chromosome, result.chromosome)
    assert np.array_equal(back.position, result.position)
