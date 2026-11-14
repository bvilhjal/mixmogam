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

pytestmark = pytest.mark.external

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
    assert b.extra["stage1_info"]["ld_cache_exact"]
    # other covariates: the cohort cover is reused (SG Q3), flagged inexact
    c = _hratt(y, gt, X[:, :1], **kw)
    ci = c.extra["stage1_info"]
    assert ci["ld_cache"] == "hit" and not ci["ld_cache_exact"]
    assert not ci["ld_cache_covariates_match"] and ci["ld_cache_overlap"] == 1.0
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
    assert a.extra["stage1_info"]["ld"] == "fast-int8"
    assert a.extra["stage1_info"]["ld_source"] == "d8-hardcall"
    assert a.extra["stage1_info"]["ldpred3_version"] == ldpred3.__version__
    assert a.extra["stage1_info"]["experimental"] is True
    ch = a.extra["stage1_info"]["fast_chains"]
    assert len(ch) == 3 and all(np.isfinite(d["h2_est"]) and d["rhat_h2"] is not None for d in ch)
    # FastLD covers cache on disk as npz + json (no pickle)
    kw = dict(stage1="ldpred3-fast", stage1_kwargs={**opts, "ld_cache": str(tmp_path / "f")})
    c1, c2 = _hratt(y, gt, X, **kw), _hratt(y, gt, X, **kw)
    assert c2.extra["stage1_info"]["ld_cache"] == "hit"
    assert np.array_equal(_chi2(c1), _chi2(c2), equal_nan=True) and np.array_equal(_chi2(c1), _chi2(a), equal_nan=True)
    import os
    files = os.listdir(tmp_path / "f")
    assert not any(f.endswith(".pkl") for f in files) and any(f.endswith(".npz") for f in files)
    # a fresh LDCache instance (as a new process) loads the npz cover
    c3 = _hratt(y, gt, X, stage1="ldpred3-fast",
                stage1_kwargs={**opts, "ld_cache": _ld8.LDCache(str(tmp_path / "f"))})
    assert c3.extra["stage1_info"]["ld_cache"] == "hit"
    assert np.array_equal(_chi2(c3), _chi2(a), equal_nan=True)
    with pytest.raises(ValueError, match="unknown keywords"):
        _hratt(y, gt, X, stage1="ldpred3-fast", stage1_kwargs={"bogus": 1})


def test_ldpred3_fast_warns_once_and_records_experimental(trait_panel):
    """ldpred3 0.10.4 warns that fast is experimental; HRATT surfaces it
    once per call (not once per fold) and records it in stage1_info."""
    gt, y, X = trait_panel
    opts = {"ld": "d8-hardcall", "fast_block_size": 300, "width": 120, "landmarks": 8,
            "burn_in": 20, "num_iter": 20, "seed": 1}
    with pytest.warns(UserWarning, match=r"ldpred3-fast.*experimental") as rec:
        twostep.hratt(y, gt, X, loco_folds=2, random_state=0, stage1="ldpred3-fast",
                      stage1_kwargs=opts)
    # one HRATT warning; the per-fold ldpred3 warnings are swallowed
    msgs = [str(w.message) for w in rec.list if "experimental" in str(w.message).lower()]
    assert sum(m.startswith("stage1 'ldpred3-fast'") for m in msgs) == 1


def test_fast_ldpred3_imports_on_pinned_package():
    from ldpred3.fast import FastLD
    import ldpred3.fast as fast_mod
    from ldpred3.diagnostics import scalar_mixing
    assert FastLD is fast_mod.FastLD and callable(scalar_mixing)
    assert ldpred3.__version__  # recorded in stage1_info


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


# --------------------------------------------------------------------------- Phase 2
def _train_corr(lg, idx, rows):
    Z = lg.rows(idx)[:, rows]
    Z = Z - Z.mean(axis=1, keepdims=True)
    d = np.einsum("ij,ij->i", Z, Z)
    return (Z @ Z.T) / np.sqrt(np.outer(d, d))


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("q", [0, 4])
def test_split_train_ld_matches_explicit_rows(ldp, packed, q):
    """ld_split='train': each training set's cover is the centred correlation
    of the projected rows restricted to that set, from one pass, for the
    float path and the hard-call algebra (Q_P'Q_P, u = Q_P'z terms)."""
    G, gt = ldp
    n = G.shape[0]
    X = None if q == 0 else np.random.default_rng(2).normal(size=(n, q))
    lg = _lg(G, gt.chromosome, X, packed)
    rng = np.random.default_rng(5)
    folds = rng.integers(0, 3, n)
    held = rng.random(n) < 0.2
    labels, n_parts, excluded = _ld8.split_parts(n, folds, held)
    order = np.arange(lg.m)
    sets = [folds != f for f in range(3)] + [~held]
    for kind in ("f32", "f32-hardcall"):
        covers = _ld8.build_split_blocks(lg, order, 300, kind, labels, n_parts, excluded, tile=97)
        assert len(covers) == 4
        for cover, rows in zip(covers, sets):
            for R, idx in cover[:3]:
                ref = _train_corr(lg, order[idx], rows)
                np.fill_diagonal(ref, 1.0)
                assert np.max(np.abs(R - ref)) < 3e-6, kind
    d8 = _ld8.build_split_blocks(lg, order, 300, "d8-hardcall", labels, n_parts, excluded)
    f32 = _ld8.build_split_blocks(lg, order, 300, "f32-hardcall", labels, n_parts, excluded)
    for (Rd, _), (Rf, _) in zip(d8[1], f32[1]):
        assert Rd.dtype == np.int8 and np.max(np.abs(Rd / 127.0 - Rf)) <= 0.5 / 127 + 1e-6


def test_split_train_stage1_cache_and_numerics():
    G, gt, y = _family_free_panel()
    caches = _ld8.LDCache()
    base = twostep.hratt(y, gt, None, loco_folds=2, stage1="ldpred3-inf", random_state=1,
                        stage1_kwargs={"ld": "d8-hardcall"})
    r1 = twostep.hratt(y, gt, None, loco_folds=2, stage1="ldpred3-inf", random_state=1,
                      stage1_kwargs={"ld": "d8-hardcall", "ld_split": "train", "ld_cache": caches})
    info = r1.extra["stage1_info"]
    assert info["ld_split"] == "train" and info["ld_cache"] == "miss"
    # k = 2 folds + hold-out complement: three covers.
    assert info["ld_bytes"] == 3 * base.extra["stage1_info"]["ld_bytes"]
    r2 = twostep.hratt(y + 0.0, gt, None, loco_folds=2, stage1="ldpred3-inf", random_state=1,
                      stage1_kwargs={"ld": "d8-hardcall", "ld_split": "train", "ld_cache": caches})
    assert r2.extra["stage1_info"]["ld_cache"] == "hit"
    assert np.allclose(r1.f_stat, r2.f_stat, equal_nan=True)
    # A different fold split must not hit the cache.
    r3 = twostep.hratt(y, gt, None, loco_folds=2, stage1="ldpred3-inf", random_state=2,
                      stage1_kwargs={"ld": "d8-hardcall", "ld_split": "train", "ld_cache": caches})
    assert r3.extra["stage1_info"]["ld_cache"] == "miss"
    assert _corr(_chi2(r1), _chi2(base)) > 0.999


def _family_free_panel(n=900, m=1200, seed=3):
    G, gt = ld_panel(n, m, seed=seed)
    rng = np.random.default_rng(seed)
    beta = np.zeros(m)
    c = rng.choice(m, 40, replace=False)
    beta[c] = rng.normal(size=40)
    g = (G - G.mean(0)) @ beta
    y = g / g.std() * np.sqrt(0.5) + rng.normal(size=n) * np.sqrt(0.5)
    return G, gt, y


def test_spectrum_guard_and_ridge(ldp):
    G, gt = ldp
    # A rank-deficient block: more variants than samples makes R singular,
    # so the D8 rounding pushes lambda_min below zero.
    lg = _lg(G[:150], gt.chromosome)
    idx = np.arange(0, 400)
    D8, R = _ld8.hardcall_block(lg, idx, return_float=True)
    lam_exact = _ld8.block_lambda_min(R)
    lam8 = _ld8.block_lambda_min(D8)
    assert lam_exact > -1e-8 and lam8 < -1e-3
    # Lanczos above the dense threshold agrees with LAPACK.
    assert abs(_ld8.block_lambda_min(D8, dense_max=10) - lam8) < 1e-3
    R2, c, after = _ld8.ridge_block_d8(D8, lam8, margin=1e-3)
    assert R2.dtype == np.int8 and np.all(np.diag(R2) == 127) and after >= 0 and c > 0
    s_seq, s_min = stage1._lassosum2_s_floor((0.2, 0.5, 0.9), -0.4)
    assert abs(s_min - 0.401 / 1.401) < 1e-12 and min(s_seq) == s_min
    assert stage1._lassosum2_s_floor((0.2, 0.5), 0.0) == ((0.2, 0.5), 0.0)


@pytest.mark.parametrize("mode", ["check", "floor", "ridge"])
def test_psd_modes_run_in_stage1(mode):
    G, gt, y = _family_free_panel()
    for meth in ("ldpred3-inf", "ldpred3-lassosum2"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = twostep.hratt(y, gt, None, loco_folds=2, stage1=meth, random_state=1,
                             stage1_kwargs={"ld": "d8-proj", "ld_psd": mode})
        info = r.extra["stage1_info"]
        assert info["ld_psd"] == mode and info["ld_spectrum"][0]["n_blocks"] >= 1
        if mode == "floor" and meth == "ldpred3-lassosum2":
            assert "s_floor" in info
    with pytest.raises(ValueError, match="ld_psd"):
        twostep.hratt(y, gt, None, loco_folds=2, stage1="ldpred3-inf",
                     stage1_kwargs={"ld": "d8-proj", "ld_psd": "bogus"})


def _two_pop(n=800, m=1200, fst=0.1, seed=0):
    rng = np.random.default_rng(seed)
    p = rng.uniform(0.1, 0.5, m)
    a, b = p * (1 - fst) / fst, (1 - p) * (1 - fst) / fst
    f = np.clip(np.stack([rng.beta(a, b), rng.beta(a, b)]), 0.02, 0.98)
    pop = rng.random(n) < 0.5
    fi = f[pop.astype(int)]
    G = ((rng.random((n, m)) < fi).astype(np.int8) + (rng.random((n, m)) < fi)).astype(np.int8)
    gt = Genotypes(G, chromosome=np.repeat([1, 2, 3], m // 3), position=np.tile(np.arange(m // 3), 3))
    c = rng.choice(m, 30, replace=False)
    g = (G[:, c] - G[:, c].mean(0)) @ rng.normal(size=30)
    y = g / g.std() * 0.6 + rng.normal(size=n) * 0.8 + 0.8 * pop
    return gt, y, pop


def test_strong_structure_gate_and_auto_fallback():
    """Strong structure after the covariates: every LDpred3 stage 1 refuses
    (block-diagonal LD cannot absorb ancestry) with advice to add PCs, the
    override runs with a warning, the elastic net and auto warn and auto
    keeps the refitted elastic net; ancestry covariates lift the gate."""
    gt, y, pop = _two_pop()
    for s1 in ("ldpred3-inf", "ldpred3-lassosum2", "ldpred3-fast"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(ValueError, match="(?s)strong population structure.*principal components"):
                twostep.hratt(y, gt, None, loco_folds=2, random_state=0, stage1=s1)
    with pytest.warns(RuntimeWarning, match="elastic-net stage 1 with the lambda calibration"):
        r = twostep.hratt(y, gt, None, loco_folds=2, random_state=0, stage1="auto")
    sel = r.extra["stage1_selection"]
    assert sel["chosen"] == "enet" and sel["reason"].startswith("strong structure")
    assert r.extra["loco_refit"]
    with pytest.warns(RuntimeWarning, match="null lambda_GC still reached 1.01-1.13"):
        twostep.hratt(y, gt, None, loco_folds=2, random_state=0, stage1="enet")
    with pytest.warns(RuntimeWarning, match="allow_strong_structure=True"):
        r = twostep.hratt(y, gt, None, loco_folds=2, random_state=0, stage1="ldpred3-inf",
                          stage1_kwargs={"allow_strong_structure": True})
    assert r.extra["stage1"] == "ldpred3-inf"
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        warnings.filterwarnings("ignore", message="(?!.*population structure)")
        r = twostep.hratt(y, gt, pop[:, None].astype(float), loco_folds=2, random_state=0,
                          stage1="ldpred3-inf")
    assert not r.extra["structure"]["strong"] and r.extra["stage1"] == "ldpred3-inf"


def test_strong_structure_warning_with_relatedness():
    """With relatedness kappa_g replaces lambda, so the structure warning
    says nothing removes the ancestry; LDpred3 still refuses."""
    gt, y, pop = _two_pop()
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        r = twostep.hratt(y, gt, None, loco_folds=2, random_state=0, stage1="enet",
                          relatedness="auto")
    msgs = [str(x.message) for x in rec if "kappa_g replaces the lambda calibration" in str(x.message)]
    assert len(msgs) == 1
    assert "LDpred3" in msgs[0] and "in simulations" in msgs[0]
    assert r.extra["structure"]["strong"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError, match="principal components"):
            twostep.hratt(y, gt, None, loco_folds=2, random_state=0, stage1="ldpred3-inf",
                          relatedness="auto")
        # A supplied pedigree / pair file (not KING-auto) also mentions
        # missing or unmodelled relatives as a possible cause.
        pairs = (np.array([0]), np.array([1]), np.array([0.5]))
        with pytest.raises(ValueError, match="relatives missing from it"):
            twostep.hratt(y, gt, None, loco_folds=2, random_state=0, stage1="ldpred3-inf",
                          relatedness=pairs)


def test_weights_gate_and_auto_fallback():
    G, gt, y = _family_free_panel()
    w = np.random.default_rng(0).uniform(0.5, 2.0, y.size)
    with pytest.raises(ValueError, match="unweighted samples"):
        twostep.hratt(y, gt, None, loco_folds=2, stage1="ldpred3-inf", sample_weights=w)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = twostep.hratt(y, gt, None, loco_folds=2, random_state=1, stage1="auto",
                          sample_weights=w)
    assert r.extra["stage1_selection"]["chosen"] == "enet"
    assert r.extra["stage1_selection"]["reason"].startswith("sample_weights")


def test_binary_values_as_quantitative_warn_with_ldpred3():
    """Count the 0/1 warning on the ldpred3 path: a duplicate must fail."""
    G, gt, y = _family_free_panel()
    yb = (y > np.median(y)).astype(float)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always", UserWarning)
        twostep.hratt(yb, gt, None, loco_folds=2, random_state=1, stage1="ldpred3-inf")
    n_two = sum("two values" in str(w.message) for w in rec)
    assert n_two == 1, n_two
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        twostep.hratt(y, gt, None, loco_folds=2, random_state=1, stage1="ldpred3-inf")


def test_cache_cohort_reuse_across_sample_subsets(trait_panel):
    """Key = variants (model-SNP set) + layout + versions, not samples or
    covariates: a trait on 80 % of the samples reuses the cohort cover and
    records overlap / coverage / F_ST; a disjoint panel with the same
    variants is rebuilt."""
    gt, y, X = trait_panel
    cache = _ld8.LDCache()
    kw = dict(stage1="ldpred3-inf", stage1_kwargs={"block_size": 200, "ld": "d8-hardcall",
                                                  "ld_cache": cache})
    _hratt(y, gt, X, **kw)
    keep = np.random.default_rng(3).random(y.size) < 0.8
    G = np.asarray(gt.G)
    sub = Genotypes(np.ascontiguousarray(G[keep]), chromosome=gt.chromosome, position=gt.position,
                    variant_ids=gt.variant_ids)
    r = _hratt(y[keep], sub, X[keep], **kw)
    info = r.extra["stage1_info"]
    assert info["ld_cache"] == "hit" and not info["ld_cache_exact"]
    assert info["ld_cache_overlap"] > 0.99 and 0.7 < info["ld_cache_coverage"] < 0.9
    assert info["ld_cache_fst"] < 0.01
    fresh = _hratt(y[keep], sub, X[keep], stage1="ldpred3-inf",
                   stage1_kwargs={"block_size": 200, "ld": "d8-hardcall"})
    assert _corr(_chi2(r), _chi2(fresh)) > 0.999
    # a different panel (no shared samples) with the same variant layout
    other, _ = ld_panel(800, 1200, seed=99)
    og = Genotypes(np.ascontiguousarray(other), chromosome=gt.chromosome, position=gt.position,
                   variant_ids=gt.variant_ids)
    r2 = _hratt(y, og, X, **kw)
    assert r2.extra["stage1_info"]["ld_cache"] == "rebuilt"


def test_cache_layout_validated_and_versioned(trait_panel, tmp_path):
    import json
    import os

    gt, y, X = trait_panel
    d = tmp_path / "c"
    kw = dict(stage1="ldpred3-inf", stage1_kwargs={"block_size": 200, "ld": "d8-hardcall",
                                                  "ld_cache": str(d)})
    _hratt(y, gt, X, **kw)
    js = [f for f in os.listdir(d) if f.endswith(".json")][0]
    with open(d / js) as fh:
        lay = json.load(fh)
    assert lay["versions"]["ldpred3"] == ldpred3.__version__ and "quantizer" in lay["versions"]
    bad = dict(lay, blocks=[[0, 7, 0]] + lay["blocks"][1:])
    with open(d / js, "w") as fh:
        json.dump(bad, fh)
    with pytest.raises(ValueError, match="layout"):
        _ld8.LDCache(str(d)).get(js[:-5])
    # another ldpred3 version: the stored cover is ignored (a miss)
    with open(d / js, "w") as fh:
        json.dump(dict(lay, versions={**lay["versions"], "ldpred3": "0.0"}), fh)
    assert _ld8.LDCache(str(d)).get(js[:-5]) is None


def test_hardcall_needs_no_private_projection(ldp, monkeypatch):
    G, gt = ldp
    lg = _lg(G, gt.chromosome, np.random.default_rng(2).normal(size=(G.shape[0], 4)))
    idx = np.arange(0, 200)
    ref = _projected_corr(lg, idx)
    lg._projection = None  # the hard-call build must not read it
    _, R = _ld8.hardcall_block(lg, idx, return_float=True)
    assert np.max(np.abs(R - ref)) < 2e-6


def test_ldpred3_relatives_skip_refit_and_apply_kappa_once():
    """Explicit ldpred3-inf + relatives: no LOCO refit (block-exact LOCO
    given shared h2), and λ_g = 1/(κ_g c_X) is applied a single time.
    Family folds are skipped."""
    # Same construction as test_model_snps_and_relatedness_with_ldpred3:
    # relatives, but not strong structure (so the LDpred3 gate does not fire).
    G, gt, y, pairs, causal = family_panel(n_fam=60, kids=3, m=3000, n_chrom=4, seed=3)
    kw = dict(stage1="ldpred3-inf", stage1_kwargs={"block_size": 200})
    r = _hratt(y, gt, relatedness=pairs, relatedness_options={"crossfit_folds": "family"},
               **kw)
    assert r.extra["stage1"] == "ldpred3-inf"
    assert not r.extra["structure"]["strong"]
    assert r.extra["loco_refit"] is False
    assert r.extra["loco_refit_reason"] == "n/a (block-exact LOCO)"
    assert r.extra["crossfit_folds"] == "random"
    ff = r.extra["relatedness"].get("family_folds") or {}
    assert ff.get("applied") is False and ff.get("reason") == "n/a (block-exact LOCO)"
    base = _hratt(y, gt, relatedness=pairs, relatedness_correction="none", **kw)
    k, cx = r.extra["relatedness"]["kappa"], r.extra["relatedness"]["c_X"]
    np.testing.assert_allclose(r.f_stat, base.f_stat / (k[gt.chromosome - 1] * cx),
                               rtol=1e-10)


def test_auto_picking_inf_reports_actual_loco_refit_not_the_decision():
    """When stage1='auto' picks the ldpred3-inf candidate, relatives
    would have triggered a refit for enet, but inf is block-exact and
    is fit with refit=False. extra['loco_refit'] must record what
    actually ran (False), not the pre-selection decision."""
    from tests._families import family_panel

    G, gt, y, pairs, causal = family_panel(n_fam=60, kids=3, m=3000, n_chrom=4, seed=3)
    orig = twostep._hratt_fit_scores

    def high_mse(*a, **k):
        out = orig(*a, **k)
        mse = np.asarray(out["mse"], dtype=np.float64)
        out["mse"] = mse + 1e6
        return out

    twostep._hratt_fit_scores = high_mse
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = twostep.hratt(y, gt, None, loco_folds=2, random_state=1,
                              stage1="auto", relatedness=pairs,
                              model_snps_h2_check=False,
                              stage1_kwargs={"block_size": 200})
    finally:
        twostep._hratt_fit_scores = orig
    assert r.extra["stage1_selection"]["chosen"] == "ldpred3-inf"
    assert r.extra.get("relatedness", {}).get("relatives_evidence") is True
    assert r.extra["loco_refit"] is False
    assert r.extra["stage1"] == "ldpred3-inf"


def test_ldpred3_lassosum2_refit_with_relatives():
    """lassosum2 is not block-exact: relatives auto-refit must run, and
    extra['loco_refit'] must report the refit that actually ran."""
    G, gt, y, pairs, causal = family_panel(n_fam=60, kids=3, m=3000, n_chrom=4, seed=3)
    kw = dict(stage1="ldpred3-lassosum2", stage1_kwargs={"block_size": 200},
              model_snps_h2_check=False)
    r = _hratt(y, gt, relatedness=pairs, **kw)
    assert r.extra["stage1"] == "ldpred3-lassosum2"
    assert not r.extra["structure"]["strong"]
    assert r.extra["loco_refit"] is True
    assert r.extra["loco_refit_reason"] == "relatives"
    off = _hratt(y, gt, relatedness=pairs, relatedness_options={"loco_refit": False}, **kw)
    assert off.extra["loco_refit"] is False
    assert not np.allclose(r.f_stat, off.f_stat, rtol=1e-5, atol=1e-5)
    inf = _hratt(y, gt, relatedness=pairs, stage1="ldpred3-inf",
                 stage1_kwargs={"block_size": 200}, model_snps_h2_check=False)
    assert inf.extra["loco_refit"] is False
    assert inf.extra["loco_refit_reason"] == "n/a (block-exact LOCO)"
