"""int8 (D8) LD backends, the LD-block cache and the LDpred3-fast stage 1
(branch stage1-ldpred3-int8; :mod:`mixmogam._ld8`).

Algebra of the hard-call build against explicit projected rows; D8 against
the float32 baseline; model_snps / relatedness with every LDpred3 stage 1;
LOCO scores across representations and against the elastic net; the cache;
and the memory bound of the hard-call build (no n x k float block).
"""

import tracemalloc
import warnings

import numpy as np
import pytest

from mixmogam import _ld8, stage1, twostep
from mixmogam._packed import PackedCalls
from mixmogam.genotypes import Genotypes
from tests._families import family_panel
from tests.test_stage1_hm3 import _hm3_trait_panel, hm3map  # noqa: F401  (fixture)
from tests.test_stage1_subset import ld_panel

ldpred3 = pytest.importorskip("ldpred3")


def _missing(G, rate, seed):
    G = G.copy()
    G[np.random.default_rng(seed).random(G.shape) < rate] = -1
    return G


def _lg(G, chrom, X=None, packed=False, y=None):
    calls = PackedCalls.pack(np.ascontiguousarray(G)) if packed else np.ascontiguousarray(G)
    gt = Genotypes(calls, chromosome=chrom, position=np.arange(G.shape[1]))
    y = np.random.default_rng(0).normal(size=G.shape[0]) if y is None else y
    return twostep._setup(y, gt, X, 25, 4096).lg


def _projected_corr(lg, idx):
    Z = lg.rows(idx)
    Z -= Z.mean(axis=1, keepdims=True)
    d = np.einsum("ij,ij->i", Z, Z)
    return (Z @ Z.T) / np.sqrt(np.outer(d, d))


@pytest.fixture(scope="module")
def ldp():
    G, gt = ld_panel(700, 1600, seed=11)
    return _missing(G, 0.01, 1), gt


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("q", [0, 4])
@pytest.mark.parametrize("tile", [None, 97])
def test_hardcall_matches_projected_float_ld(ldp, packed, q, tile):
    G, gt = ldp
    X = None if q == 0 else np.random.default_rng(2).normal(size=(G.shape[0], q))
    lg = _lg(G, gt.chromosome, X, packed)
    idx = np.arange(37, 437)
    D8, R = _ld8.hardcall_block(lg, idx, tile=tile, return_float=True)
    Rf = _projected_corr(lg, idx)
    # Before rounding: the baseline's rows are float32, (b) is exact.
    assert np.max(np.abs(R - Rf)) < 2e-6
    assert np.allclose(np.diag(R), 1.0) and np.array_equal(R, R.T)
    # D8: within half a code of exact, at most one code from the float32 path.
    assert D8.dtype == np.int8 and np.array_equal(D8, D8.T) and np.all(np.diag(D8) == 127)
    assert np.max(np.abs(D8 / 127.0 - R)) <= 0.5 / 127 + 1e-6
    assert np.max(np.abs(D8.astype(int) - _ld8.quantize_d8(Rf).astype(int))) <= 1
    # float32 hard-call blocks equal the exact matrix to float32 precision.
    R32 = _ld8.hardcall_block(lg, idx, tile=tile, quantize=False)
    assert R32.dtype == np.float32 and np.max(np.abs(R32 - Rf)) < 2e-6


def test_hardcall_counts_are_exact_integers(ldp):
    G, gt = ldp
    lg = _lg(G, gt.chromosome, packed=True)
    idx = np.arange(0, 60)
    XX, XO, OO, called, sums = _ld8.hardcall_counts(lg, idx, tile=50)
    g = G[:, idx].astype(np.int64)
    o = (g >= 0).astype(np.int64)
    x = np.where(g >= 0, g, 0)
    assert np.array_equal(XX, x.T @ x) and np.array_equal(XO, x.T @ o)
    assert np.array_equal(OO, o.T @ o)
    assert np.array_equal(called, o.sum(0)) and np.array_equal(sums, x.sum(0))


def test_monomorphic_and_quantize_convention():
    rng = np.random.default_rng(5)
    G = rng.binomial(2, 0.3, size=(200, 40)).astype(np.int8)
    G[:, 7] = 1  # monomorphic
    lg = _lg(G, np.repeat([1, 2], 20))
    D8 = _ld8.hardcall_block(lg, np.arange(40))
    assert D8[7, 7] == 127 and not np.any(np.delete(D8[7], 7))
    R = np.array([[1.0, 0.5039, -1.2], [0.5039, 1.0, 0.0], [-1.2, 0.0, 1.0]])
    Q = _ld8.quantize_d8(R)
    assert Q[0, 1] == 64 and Q[0, 2] == -127 and np.all(np.diag(Q) == 127)
    from ldpred3._common import _quantize_dense
    Rs = np.clip(R, -1, 1)
    assert np.array_equal(_quantize_dense(Rs.astype(np.float32)), _ld8.quantize_d8(Rs))


def test_hardcall_build_forms_no_float_sample_block(ldp):
    """(b) keeps O(k T + k^2) bytes: far below one n x k float32 block."""
    G, gt = ld_panel(20000, 400, seed=4)
    G = _missing(G, 0.01, 3)
    lg = _lg(G, gt.chromosome, packed=True)
    idx = np.arange(0, 100)
    n, k = lg.n, idx.size
    _ld8.hardcall_block(lg, idx, tile=256)  # compile first
    tracemalloc.start()
    _ld8.hardcall_block(lg, idx, tile=256)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    # tiles of 256 samples plus a dozen k x k float64 work matrices
    assert peak < 3 * k * 256 * 4 + 16 * k * k * 8 < 0.25 * n * k * 4, (peak, n * k * 4)
    # The float baseline forms (k, n) float64 rows: 8 n k bytes.
    tracemalloc.start()
    stage1.ld_blocks(lg, idx, np.ones(lg.m), block_size=100)
    _, peak_f = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert peak_f > n * k * 4 > 4 * peak


def test_build_blocks_kinds_share_one_layout(ldp):
    G, gt = ldp
    lg = _lg(G, gt.chromosome, np.random.default_rng(1).normal(size=(G.shape[0], 2)))
    order = np.argsort(lg.groups, kind="stable")
    scale = np.sqrt(lg.zz / lg.n)
    f32 = _ld8.build_blocks(lg, order, scale, 300, "f32")
    proj = _ld8.build_blocks(lg, order, scale, 300, "d8-proj")
    hard = _ld8.build_blocks(lg, order, scale, 300, "d8-hardcall")
    for (a, ia), (b, ib), (c, ic) in zip(f32, proj, hard):
        assert np.array_equal(ia, ib) and np.array_equal(ia, ic)
        assert a.dtype == np.float32 and b.dtype == c.dtype == np.int8
        assert np.array_equal(b, _ld8.quantize_d8(a))
        assert np.max(np.abs(b.astype(int) - c.astype(int))) <= 1
    assert _ld8.cover_nbytes(proj) * 4 == _ld8.cover_nbytes(f32)
    with pytest.raises(ValueError, match="ld must be one of"):
        _ld8.build_blocks(lg, order, scale, 300, "bf16")


# --------------------------------------------------------------- hratt end to end
@pytest.fixture(scope="module")
def trait_panel():
    G, gt = ld_panel(800, 1200, seed=21)
    rng = np.random.default_rng(22)
    caus = rng.choice(1200, 60, replace=False)
    Z = G[:, caus] - G[:, caus].mean(0)
    g = Z @ rng.normal(size=60)
    y = 0.75 * g / g.std() + 0.66 * rng.standard_normal(800)
    X = rng.normal(size=(800, 2))
    return gt, y, X


def _hratt(y, gt, X=None, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return twostep.hratt(y, gt, X, loco_folds=2, random_state=0, **kw)


def _chi2(r):
    return np.asarray(r.f_stat, dtype=float)


def _corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)  # monomorphic variants test as NaN
    return np.corrcoef(a[ok], b[ok])[0, 1]


def test_default_ld_is_the_unchanged_float32_path(trait_panel):
    gt, y, X = trait_panel
    a = _hratt(y, gt, X, stage1="ldpred3-inf", stage1_kwargs={"block_size": 200})
    b = _hratt(y, gt, X, stage1="ldpred3-inf", stage1_kwargs={"block_size": 200, "ld": "f32"})
    assert np.array_equal(_chi2(a), _chi2(b), equal_nan=True)
    assert a.extra["stage1_info"]["ld"] == "f32"
    e = _hratt(y, gt, X)
    assert "stage1_info" not in e.extra and "stage1" not in e.extra


@pytest.mark.parametrize("method", ["ldpred3-inf", "ldpred3-lassosum2"])
def test_d8_backends_agree_with_float32_and_enet(trait_panel, method):
    gt, y, X = trait_panel
    out = {}
    for ld in ("f32", "d8-proj", "d8-hardcall", "f32-hardcall"):
        out[ld] = _hratt(y, gt, X, stage1=method, stage1_kwargs={"block_size": 200, "ld": ld})
    ref = _chi2(out["f32"])
    for ld in ("d8-proj", "d8-hardcall"):
        c = _chi2(out[ld])
        assert _corr(ref, c) > 0.999
        assert abs(np.nanmedian(c) / np.nanmedian(ref) - 1) < 0.01
        assert out[ld].extra["stage1_info"]["ld_bytes"] * 4 == out["f32"].extra["stage1_info"]["ld_bytes"]
    np.testing.assert_allclose(_chi2(out["f32-hardcall"]), ref, rtol=2e-3, atol=2e-3)
    enet = _chi2(_hratt(y, gt, X))
    assert _corr(enet, _chi2(out["d8-hardcall"])) > 0.95


@pytest.mark.parametrize("method", ["ldpred3-inf", "ldpred3-lassosum2"])
@pytest.mark.parametrize("ld", ["f32", "d8-hardcall"])
def test_model_snps_and_relatedness_with_ldpred3(method, ld):
    G, gt, y, pairs, causal = family_panel(n_fam=60, kids=3, m=3000, n_chrom=4, seed=3)
    kw = dict(stage1=method, stage1_kwargs={"block_size": 200, "ld": ld})
    sub = np.arange(0, 3000, 2)
    r = _hratt(y, gt, model_snps=sub, relatedness="auto", **kw)
    assert r.extra["stage1"] == method and r.extra["model_snps"]["n"] == sub.size
    rel = r.extra["relatedness"]
    assert rel["correction"] == "kappa" and np.all(np.isfinite(rel["kappa"]))
    assert rel["n_pairs"] > 0
    assert np.all(np.isfinite(_chi2(r)))
    r2 = _hratt(y, gt, model_snps="auto", relatedness=pairs, **kw)
    assert r2.extra["model_snps"]["source"].startswith("auto")
    assert np.all(np.isfinite(_chi2(r2)))


def test_model_snps_hm3_with_ldpred3(hm3map):  # noqa: F811
    from mixmogam._hm3 import load_hm3
    hm = load_hm3(hm3map)
    _, gt, y = _hm3_trait_panel(hm)
    for method in ("ldpred3-inf", "ldpred3-lassosum2"):
        r = _hratt(y, gt, stage1=method, model_snps="hm3",
                   model_snps_options={"map": hm3map, "min_snps": 100},
                   stage1_kwargs={"block_size": 100, "ld": "d8-hardcall"})
        assert r.extra["model_snps"]["source"] == "hm3" and r.extra["stage1"] == method


def test_ld_cache_reuses_one_build_across_traits(trait_panel, tmp_path):
    gt, y, X = trait_panel
    cache = _ld8.LDCache()
    kw = dict(stage1="ldpred3-inf", stage1_kwargs={"block_size": 200, "ld": "d8-hardcall",
                                                  "ld_cache": cache})
    a = _hratt(y, gt, X, **kw)
    y2 = np.random.default_rng(9).permutation(y)
    b = _hratt(y2, gt, X, **kw)
    assert a.extra["stage1_info"]["ld_cache"] == "miss"
    assert b.extra["stage1_info"]["ld_cache"] == "hit" and cache.hits == 1
    fresh = _hratt(y2, gt, X, stage1="ldpred3-inf",
                   stage1_kwargs={"block_size": 200, "ld": "d8-hardcall"})
    assert np.array_equal(_chi2(b), _chi2(fresh), equal_nan=True)
    # other covariates: another key
    c = _hratt(y, gt, X[:, :1], **kw)
    assert c.extra["stage1_info"]["ld_cache"] == "miss"
    # on disk: a new process (new LDCache) reads the memory-mapped payload
    for ld in ("d8-hardcall", "f32"):
        d = tmp_path / ld
        kw_disk = dict(stage1="ldpred3-lassosum2",
                       stage1_kwargs={"block_size": 200, "ld": ld, "ld_cache": str(d)})
        first = _hratt(y, gt, X, **kw_disk)
        second = _hratt(y, gt, X, **kw_disk)
        assert first.extra["stage1_info"]["ld_cache"] == "miss"
        assert second.extra["stage1_info"]["ld_cache"] == "hit"
        assert np.array_equal(_chi2(first), _chi2(second), equal_nan=True)
        assert _ld8.LDCache(str(d)).disk_bytes() >= first.extra["stage1_info"]["ld_bytes"]


def test_ldpred3_fast_loco_stage1(trait_panel, tmp_path):
    gt, y, X = trait_panel
    opts = {"ld": "d8-hardcall", "fast_block_size": 300, "width": 120, "landmarks": 8,
            "burn_in": 40, "num_iter": 40, "seed": 1}
    a = _hratt(y, gt, X, stage1="ldpred3-fast", stage1_kwargs=opts)
    b = _hratt(y, gt, X, stage1="ldpred3-fast", stage1_kwargs=opts)
    assert a.extra["stage1"] == "ldpred3-fast"
    assert np.array_equal(_chi2(a), _chi2(b), equal_nan=True)  # seeded chains
    inf = _hratt(y, gt, X, stage1="ldpred3-inf", stage1_kwargs={"block_size": 300})
    assert _corr(_chi2(a), _chi2(inf)) > 0.95
    # FastLD covers cache (pickled) too
    kw = dict(stage1="ldpred3-fast", stage1_kwargs={**opts, "ld_cache": str(tmp_path / "f")})
    c1, c2 = _hratt(y, gt, X, **kw), _hratt(y, gt, X, **kw)
    assert c2.extra["stage1_info"]["ld_cache"] == "hit"
    assert np.array_equal(_chi2(c1), _chi2(c2), equal_nan=True) and np.array_equal(_chi2(c1), _chi2(a), equal_nan=True)
    with pytest.raises(ValueError, match="unknown keywords"):
        _hratt(y, gt, X, stage1="ldpred3-fast", stage1_kwargs={"bogus": 1})


def test_fast_blocks_stay_group_pure(trait_panel):
    """LOCO stays an exact subtraction: every fast block lies in one group."""
    gt, y, X = trait_panel
    lg = _lg(np.asarray(gt.G), gt.chromosome, X)
    order = np.argsort(lg.groups, kind="stable")
    blocks = _ld8.build_blocks(lg, order, np.ones(lg.m), 500, "d8-hardcall")
    fast = ldpred3.fast_ld_blocks(blocks, width=150, landmarks=8)
    for (F, idx) in fast:
        assert np.unique(lg.groups[order[idx]]).size == 1
        assert F.m == idx.size
