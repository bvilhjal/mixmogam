"""Spec A (SG specs-for-wizard.md): relatedness= kappa_g correction of HRATT."""
import warnings

import numpy as np
import pytest

from mixmogam._relatedness import (Relatedness, align_pairs, covariate_factor, king_pairs,
                                   read_grm_sparse, read_kin0, unrelated_set)
from mixmogam.twostep import hratt
from tests._families import family_panel


def _null_mean(res, causal):
    c = np.delete(res.f_stat, causal)
    return float(np.mean(c[np.isfinite(c)]))


@pytest.fixture(scope="module")
def fam():
    G, gt, y, pairs, causal = family_panel(seed=0)
    sub = np.sort(np.random.default_rng(9).choice(gt.n_variants, gt.n_variants // 20,
                                                  replace=False))
    return G, gt, y, pairs, causal, sub


def _run(y, gt, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return hratt(y, gt, random_state=1, model_snps_h2_check=False, **kw)


# 1. formula
def test_kappa_and_cx_match_brute_force():
    rng = np.random.default_rng(0)
    n = 6
    rel = Relatedness(n, [0, 2], [1, 5], [0.5, 0.25])
    Phi = np.eye(n)
    Phi[0, 1] = Phi[1, 0] = 0.5
    Phi[2, 5] = Phi[5, 2] = 0.25
    A = rng.standard_normal((n, 3))
    brute = np.einsum("ig,ij,jg->g", A, Phi, A) / np.einsum("ig,ig->g", A, A)
    np.testing.assert_allclose(rel.kappa(A), brute, rtol=0, atol=1e-12)
    Q, _ = np.linalg.qr(np.c_[np.ones(n), rng.standard_normal(n)])
    cx, tr = covariate_factor(rel, Q, n - 2)
    assert abs(tr - np.trace(Q.T @ Phi @ Q)) < 1e-12
    assert abs(cx - (n - 2) / (n - np.trace(Q.T @ Phi @ Q))) < 1e-12


# 2. no pairs: kappa = 1 and p identical to correction 'none'
def test_no_pairs_is_byte_identical_to_none(fam):
    G, gt, y, pairs, causal, sub = fam
    empty = (np.array([], int), np.array([], int), np.array([]))
    a = _run(y, gt, model_snps=sub, relatedness=empty, relatedness_correction="kappa")
    b = _run(y, gt, model_snps=sub, relatedness=empty, relatedness_correction="none")
    np.testing.assert_array_equal(a.extra["relatedness"]["kappa"], 1.0)
    assert a.extra["relatedness"]["c_X"] == 1.0
    assert a.p.tobytes() == b.p.tobytes()
    assert a.extra["lambda"] == b.extra["lambda"] == 1.0


# 3. calibration (quantitative), 5% and full stage 1
@pytest.mark.parametrize("subset", [True, False])
def test_nuclear_families_calibrated(fam, subset):
    G, gt, y, pairs, causal, sub = fam
    ms = sub if subset else None
    base = _run(y, gt, model_snps=ms, relatedness=pairs, relatedness_correction="none")
    kap = _run(y, gt, model_snps=ms, relatedness=pairs)
    assert _null_mean(base, causal) > 1.1  # relatives inflate the uncorrected statistic
    assert 0.97 <= _null_mean(kap, causal) <= 1.03
    assert kap.extra["lambda"] == 1.0
    assert np.all(kap.extra["relatedness"]["kappa"] > 1.05)
    s = kap.extra["structure"]
    assert s["on"] == "unrelated" and s["n_unrelated"] < gt.n_samples
    assert s["excess"] < s["all_samples"]["excess"]  # relatives trip the all-sample test only


# 4. binary + SPA: per-group lambda_g = 1/kappa_g
def test_binary_kappa_and_spa():
    G, gt, y, pairs, causal = family_panel(seed=2, prevalence=0.2)
    sub = np.sort(np.random.default_rng(3).choice(gt.n_variants, gt.n_variants // 20,
                                                  replace=False))
    base = _run(y, gt, trait="binary", model_snps=sub, relatedness=pairs,
                relatedness_correction="none")
    kap = _run(y, gt, trait="binary", model_snps=sub, relatedness=pairs, spa_threshold=1.0)
    k = kap.extra["relatedness"]["kappa"]
    assert k.shape == (4,) and np.all(k > 1.02)
    assert _null_mean(base, causal) > 1.05
    assert 0.95 <= _null_mean(kap, causal) <= 1.05
    assert kap.extra["n_spa"] > 0
    # the normal-tail part is exactly the uncorrected statistic over kappa_g
    plain = _run(y, gt, trait="binary", model_snps=sub, relatedness=pairs,
                 spa_threshold=1e6)
    base_n = _run(y, gt, trait="binary", model_snps=sub, relatedness=pairs,
                  relatedness_correction="none", spa_threshold=1e6)
    g = gt.chromosome - 1
    cx = plain.extra["relatedness"]["c_X"]
    assert 1.0 < cx < 1.01  # intercept only: c_X = (n - 1) / (n - tr(Q'PhiQ)) slightly above 1
    np.testing.assert_allclose(plain.f_stat, base_n.f_stat / (k[g] * cx), rtol=1e-12)


# P2: c_X on the binary and weighted paths (family-level covariates make it large)
@pytest.mark.parametrize("path", ["binary", "weighted"])
def test_cx_on_binary_and_weighted_paths(path):
    G, gt, y, pairs, causal = family_panel(seed=2, prevalence=0.2 if path == "binary" else None)
    fam_id = np.repeat(np.arange(200), 6)
    X = np.c_[np.random.default_rng(1).standard_normal(200)[fam_id],
              np.random.default_rng(2).standard_normal(200)[fam_id]]
    kw = dict(trait="binary") if path == "binary" else dict(
        sample_weights=np.random.default_rng(3).uniform(0.5, 2.0, y.size))
    kap = _run(y, gt, X=X, relatedness=pairs, spa_threshold=1e6, **kw)
    base = _run(y, gt, X=X, relatedness=pairs, relatedness_correction="none", spa_threshold=1e6,
                **kw)
    rel = kap.extra["relatedness"]
    cx, k = rel["c_X"], rel["kappa"]
    Q = np.linalg.qr(np.c_[np.ones(y.size), X])[0]
    ref, _ = covariate_factor(Relatedness(y.size, *pairs), Q, y.size - 3)
    assert cx > 1.002 and abs(cx - ref) < 1e-12
    np.testing.assert_allclose(kap.f_stat, base.f_stat / (k[gt.chromosome - 1] * cx), rtol=1e-10)
    assert 0.95 <= _null_mean(kap, causal) <= 1.05


# 5. household environment (unascertained)
def test_household_environment_calibrated():
    G, gt, y, pairs, causal = family_panel(seed=4, c2=0.2, h2=0.4)
    sub = np.sort(np.random.default_rng(5).choice(gt.n_variants, gt.n_variants // 20,
                                                  replace=False))
    kap = _run(y, gt, model_snps=sub, relatedness=pairs)
    assert 0.97 <= _null_mean(kap, causal) <= 1.03


# 5b. covariates: c_X > 1 and calibrated
def test_covariates_cx(fam):
    G, gt, y, pairs, causal, sub = fam
    fam_id = np.repeat(np.arange(200), 6)
    # family-level covariates make Q'PhiQ large: c_X clearly above 1
    X = np.c_[np.random.default_rng(1).standard_normal(200)[fam_id],
              np.random.default_rng(2).standard_normal(200)[fam_id]]
    kap = _run(y, gt, X=X, model_snps=sub, relatedness=pairs)
    assert kap.extra["relatedness"]["c_X"] > 1.002
    assert 0.97 <= _null_mean(kap, causal) <= 1.03


# 6. unrelated samples
def test_unrelated_samples_unchanged():
    G, gt, y, pairs, causal = family_panel(seed=6, kids=0, n_fam=600)
    a = _run(y, gt, relatedness="auto")
    b = _run(y, gt)
    k = a.extra["relatedness"]["kappa"]
    assert np.all(np.abs(k - 1) < 0.01)
    assert np.nanmax(np.abs(a.f_stat - b.f_stat) / np.maximum(b.f_stat, 1)) < 1e-2


# 6c. noise-floor guard
def test_king_noise_floor_raises_threshold(fam):
    G, gt, y, pairs, causal, sub = fam
    res = king_pairs(gt, np.arange(1000), 0.025)
    assert res["threshold"] > 0.025 and res["noise_sd"] > 0.025 / 4
    with pytest.warns(RuntimeWarning, match="noise SD"):
        a = hratt(y, gt, model_snps=np.setdiff1d(np.arange(gt.n_variants), np.arange(1000)),
                  relatedness="auto", random_state=1, model_snps_h2_check=False,
                  relatedness_options={"king_max_variants": 1000})
    assert a.extra["relatedness"]["king"]["threshold"] > 0.025


# 7. inputs
def test_readers_alignment_and_duplicates(tmp_path, fam):
    G, gt, y, pairs, causal, sub = fam
    ids = np.asarray(gt.sample_ids)
    i, j, r = pairs
    with open(tmp_path / "x.kin0", "w") as fh:
        fh.write("#FID1\tIID1\tFID2\tIID2\tNSNP\tHETHET\tIBS0\tKINSHIP\n")
        for a_, b_, v in zip(i, j, r):
            fh.write(f"f\t{ids[b_]}\tf\t{ids[a_]}\t100\t0.1\t0\t{v / 2}\n")
        fh.write("f\tnobody\tf\ts0\t100\t0.1\t0\t0.25\n")
    perm = np.random.default_rng(0).permutation(ids.size)
    rel, dropped = read_kin0(str(tmp_path / "x.kin0"), ids[perm])
    assert dropped == 1 and rel.n_pairs == i.size
    inv = np.argsort(perm)
    u = rel.upper.tocoo()
    got = set(zip(*sorted(zip(*[(min(perm[a_], perm[b_]), max(perm[a_], perm[b_]))
                                for a_, b_ in zip(u.row, u.col)]))))
    assert got == set(zip(*sorted(zip(*[(min(a_, b_), max(a_, b_)) for a_, b_ in zip(i, j)]))))
    assert np.all(inv[perm] == np.arange(ids.size))
    with open(tmp_path / "g.grm.id", "w") as fh:
        for s_ in ids:
            fh.write(f"f\t{s_}\n")
    with open(tmp_path / "g.grm.sp", "w") as fh:
        for k_ in range(3):
            fh.write(f"{k_}\t{k_}\t1.0\n")
        for a_, b_, v in zip(i, j, r):
            fh.write(f"{b_}\t{a_}\t{v}\n")
    rel2, d2 = read_grm_sparse(str(tmp_path / "g"), ids)
    assert d2 == 0 and rel2.n_pairs == i.size
    res = _run(y, gt, model_snps=sub, relatedness=str(tmp_path / "g.grm.sp"))
    ref = _run(y, gt, model_snps=sub, relatedness=pairs)
    np.testing.assert_allclose(res.f_stat, ref.f_stat, rtol=1e-12)
    with pytest.raises(ValueError, match="more than once"):
        align_pairs(["s0", "s1"], ["s1", "s0"], [0.5, 0.5], ids)


# 8. KING on non-model variants vs all variants
def test_king_non_model_matches_all(fam):
    G, gt, y, pairs, causal, sub = fam
    a = _run(y, gt, model_snps=sub, relatedness="auto")
    b = _run(y, gt, model_snps=sub, relatedness=pairs)
    assert "non-model" in a.extra["relatedness"]["king_on"]
    np.testing.assert_allclose(a.extra["relatedness"]["kappa"], b.extra["relatedness"]["kappa"],
                               atol=0.005)


# 9. guard
def test_kappa_floor_clips_and_warns(fam):
    G, gt, y, pairs, causal, sub = fam
    i, j, r = pairs
    with pytest.warns(RuntimeWarning, match="clipped"):
        res = hratt(y, gt, model_snps=sub, relatedness=(i, j, -1.5 * np.ones_like(r)),
                    random_state=1, model_snps_h2_check=False)
    assert np.all(res.extra["relatedness"]["kappa"] >= 0.2)


def test_unrelated_set_is_maximal_and_unrelated():
    rel = Relatedness(7, [0, 0, 0, 3, 4], [1, 2, 3, 4, 5], np.full(5, 0.5))
    keep = unrelated_set(rel)
    K = rel.Koff.toarray() != 0
    assert not K[np.ix_(keep, keep)].any()
    for s_ in np.setdiff1d(np.arange(7), keep):
        assert K[s_, keep].any()  # maximal: every dropped sample has a kept relative
    assert 6 in keep and keep.size == 5  # drops 0 (3 pairs), then 4


def test_gls_is_opt_in_and_api_errors(fam):
    G, gt, y, pairs, causal, sub = fam
    g = _run(y, gt, model_snps=sub, relatedness=pairs, relatedness_correction="gls",
             relatedness_options={"denominator": "gamma"})
    assert g.extra["relatedness"]["correction"] == "gls" and g.extra["lambda"] == 1.0
    with pytest.raises(ValueError, match="needs relatedness"):
        hratt(y, gt, relatedness_correction="kappa")
    with pytest.raises(ValueError, match="external pairs"):
        hratt(y, gt, relatedness="auto", relatedness_options={"king_max_n": 100})
    with pytest.raises(ValueError, match="gls"):
        hratt((y > 0).astype(float), gt, trait="binary", relatedness=pairs,
              relatedness_correction="gls")


# P1 (SG specA-benchmark §2b): relatives keep the per-group LOCO refit; without
# it the joint cross-fit leaks into in-model null SNPs (household c2).
def test_relatives_keep_loco_refit_in_model_nulls_calibrated():
    G, gt, y, pairs, causal = family_panel(seed=0, n_fam=200, kids=4, m=8000, c2=0.3, h2=0.5)
    on = _run(y, gt, relatedness=pairs)
    off = _run(y, gt, relatedness=pairs, relatedness_options={"loco_refit": False})
    assert on.extra["loco_refit"] and on.extra["loco_refit_reason"] == "relatives"
    assert not off.extra["loco_refit"] and off.extra["loco_refit_reason"] is None
    assert not on.extra["structure"]["strong"]  # the unrelated-set test does not fire
    assert on.extra["relatedness"]["n_pairs_above_threshold"] == pairs[0].size
    assert on.extra["relatedness"]["s2_above_threshold"] > 1.0  # 14 pairs of r=0.5 per 6
    assert _null_mean(off, causal) > 1.03  # every null here is a model SNP
    assert _null_mean(on, causal) < 1.02


def test_loco_refit_option_and_no_close_pairs():
    G, gt, y, pairs, causal = family_panel(seed=6, kids=0, n_fam=300, m=4000)
    far = (np.array([0, 2]), np.array([1, 3]), np.array([0.04, 0.03]))  # phi <= 0.025
    a = _run(y, gt, relatedness=far)
    assert not a.extra["loco_refit"] and a.extra["relatedness"]["n_pairs_above_threshold"] == 0
    b = _run(y, gt, relatedness=far, relatedness_options={"loco_refit": True})
    assert b.extra["loco_refit"] and b.extra["loco_refit_reason"] == "forced"
    few = (np.array([0]), np.array([1]), np.array([0.5]))  # one pair: s2 = 2 * 0.25 / 600 < 0.005
    c = _run(y, gt, relatedness=few)
    assert not c.extra["loco_refit"] and c.extra["relatedness"]["n_pairs_above_threshold"] == 1
    d = _run(y, gt, relatedness=few, relatedness_options={"loco_refit_min_s2": 0.0})
    assert d.extra["loco_refit"] and d.extra["loco_refit_reason"] == "relatives"
    with pytest.raises(ValueError, match="loco_refit"):
        hratt(y, gt, relatedness=far, relatedness_options={"loco_refit": "yes"})


# P3: the unrelated set splits close pairs only, with a minimum size; KING
# 'auto' refuses haploid / inbred coding.
def test_unrelated_set_threshold_and_min_size():
    rel = Relatedness(5, [0, 1, 2], [1, 2, 3], [0.5, 0.05, 0.02])  # phi 0.25, 0.025, 0.01
    assert unrelated_set(rel).size == 3  # every pair split
    keep = unrelated_set(rel, min_phi=0.0442)
    assert keep.size == 4 and not ({0, 1} <= set(keep))
    G, gt, y, pairs, causal = family_panel(seed=3, n_fam=30, kids=4, m=2000, n_causal=100)
    with pytest.warns(RuntimeWarning, match="unrelated set"):
        r = hratt(y, gt, relatedness=pairs, random_state=1, model_snps_h2_check=False)
    s = r.extra["structure"]
    assert s["on"].startswith("all") and s["n_unrelated"] < 100
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        r2 = hratt(y, gt, relatedness=pairs, random_state=1, model_snps_h2_check=False,
                   relatedness_options={"unrelated_min_n": 10})
    assert r2.extra["structure"]["on"] == "unrelated"


@pytest.mark.parametrize("coding", ["haploid01", "inbred02"])
def test_king_auto_refuses_haploid_or_inbred(coding):
    rng = np.random.default_rng(0)
    n, m = 300, 3000
    f = rng.uniform(0.1, 0.5, m)
    hap = (rng.random((n, m)) < f).astype(np.int8)
    from mixmogam import Genotypes
    G = hap if coding == "haploid01" else 2 * hap
    gt = Genotypes(G, chromosome=np.repeat([1, 2], m // 2), position=np.arange(m) * 10_000 + 1)
    y = rng.standard_normal(n)
    with pytest.raises(ValueError, match="KING-robust assumes outbred diploid"):
        hratt(y, gt, relatedness="auto", random_state=1, model_snps_h2_check=False)


def test_king_auto_records_inbreeding_on_outbred(fam):
    G, gt, y, pairs, causal, sub = fam
    a = _run(y, gt, model_snps=sub, relatedness="auto")
    k = a.extra["relatedness"]["king"]
    assert abs(k["median_F"]) < 0.05 and k["hom_alt_fraction"] > 0.05
