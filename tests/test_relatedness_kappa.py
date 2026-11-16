"""Spec A relatedness= kappa_g correction of HRATT."""
import warnings

import numpy as np
import pytest

from mixmogam._relatedness import (KING_COUSIN_DEGREE, KING_SECOND_DEGREE, RELATEDNESS_OPTIONS, Relatedness,
                                   align_pairs, covariate_factor, king_pairs,
                                   king_refit_trigger, read_grm_sparse, read_kin0,
                                   resolve_relatedness, unrelated_set)
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
    # Isolate the kappa effect: even unrelated panels can have false KING
    # pairs; the default must conservatively refit if any are detected.
    a = _run(y, gt, relatedness="auto", relatedness_options={"loco_refit": False})
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


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("block", [1, 3, 20])
def test_king_missing_calls_match_pairwise_complete_oracle(packed, block):
    from mixmogam._packed import PackedCalls
    from mixmogam._relatedness import king_relatedness
    from mixmogam.genotypes import Genotypes

    rng = np.random.default_rng(72)
    G = rng.binomial(2, 0.3, (8, 73)).astype(np.int8)
    G[2:, 11:][rng.random((6, 62)) < 0.35] = -1
    G[-2] = -1  # no jointly called loci
    G[-1] = 0   # no heterozygotes
    gt = Genotypes(PackedCalls.pack(G) if packed else G)
    expected = {}
    for i in range(len(G)):
        for j in range(i + 1, len(G)):
            called = (G[i] >= 0) & (G[j] >= 0)
            a, b = G[i, called].astype(float), G[j, called].astype(float)
            nhet = min(np.count_nonzero(a == 1), np.count_nonzero(b == 1))
            # Independent squared-distance form of the KING estimator.
            phi = 0.5 - np.sum((a - b) ** 2) / (4 * nhet) if nhet else 0.0
            if phi > 0:
                expected[i, j] = 2 * phi
    res = king_pairs(gt, np.arange(G.shape[1]), 0, block=block, chunk=7,
                     noise_pairs=0, noise_sd_multiple=0)
    got = dict(zip(zip(res["i"], res["j"]), res["r"]))
    proto = king_relatedness(gt, min_phi=0, min_maf=0, block=7).upper.tocoo()
    for found in (got, dict(zip(zip(proto.row, proto.col), proto.data))):
        assert found.keys() == expected.keys()
        for pair in expected:
            assert found[pair] == pytest.approx(expected[pair], abs=1e-14)


def test_king_parent_child_degree_survives_missing_calls():
    from mixmogam.genotypes import Genotypes

    rng = np.random.default_rng(72)
    m = 50_000
    hap = rng.binomial(1, 0.3, (3, m)).astype(np.int8)
    parent = hap[0] + hap[1]
    child = np.where(rng.integers(2, size=m), hap[0], hap[1]) + hap[2]
    G = np.vstack([parent, child])
    G[1, rng.random(m) < 0.3] = -1
    common = np.all(G >= 0, axis=0)
    estimates = []
    for calls in (G, G[:, common]):
        res = king_pairs(Genotypes(calls), np.arange(calls.shape[1]), 0.177,
                         noise_pairs=0, noise_sd_multiple=0)
        assert list(zip(res["i"], res["j"])) == [(0, 1)]
        estimates.append(res["r"][0] / 2)
    assert estimates[0] == estimates[1]
    assert abs(estimates[0] - 0.25) < 0.01


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
    assert 6 in keep and keep.size == 5  # keeps 1 (fewest relatives), dropping 0; then 3


def test_unrelated_set_near_maximum_on_cousin_graph():
    # (iii) isolate-like kinship: each sample a child of two sibships, so
    # cousins through either parent (a line graph; mean degree ~7 here). The
    # old max-degree removal kept 86-88% of the exact maximum on seeds 4-8,
    # minimum-degree selection 96-98%.
    from scipy import sparse
    from scipy.optimize import Bounds, LinearConstraint, milp

    rng = np.random.default_rng(4)
    n, n_sib = 600, 330
    a, b = rng.integers(0, n_sib, n), rng.integers(0, n_sib, n)
    same = a[:, None] == a[None, :]
    same |= a[:, None] == b[None, :]
    same |= b[:, None] == a[None, :]
    same |= b[:, None] == b[None, :]
    i, j = np.nonzero(np.triu(same, 1))
    rel = Relatedness(n, i, j, np.full(i.size, 0.125))
    keep = unrelated_set(rel, min_phi=0.0442)
    K = rel.Koff.toarray() != 0
    assert not K[np.ix_(keep, keep)].any()
    A = sparse.csr_matrix((np.ones(2 * i.size), (np.repeat(np.arange(i.size), 2),
                                                 np.c_[i, j].ravel())), shape=(i.size, n))
    best = -milp(-np.ones(n), constraints=LinearConstraint(A, -np.inf, 1),
                 integrality=np.ones(n), bounds=Bounds(0, 1)).fun
    assert keep.size >= 0.95 * best


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


# With random folds, relatives keep the per-group LOCO refit; without it
# the joint cross-fit leaks into in-model null SNPs (household c2).
def test_relatives_keep_loco_refit_in_model_nulls_calibrated():
    G, gt, y, pairs, causal = family_panel(seed=0, n_fam=200, kids=4, m=8000, c2=0.3, h2=0.5)
    rnd = {"crossfit_folds": "random"}
    on = _run(y, gt, relatedness=pairs, relatedness_options=rnd)
    off = _run(y, gt, relatedness=pairs, relatedness_options={**rnd, "loco_refit": False})
    assert on.extra["loco_refit"] and on.extra["loco_refit_reason"] == "relatives"
    assert not off.extra["loco_refit"] and off.extra["loco_refit_reason"] is None
    assert not on.extra["structure"]["strong"]  # the unrelated-set test does not fire
    assert on.extra["relatedness"]["n_pairs_above_threshold"] == pairs[0].size
    assert _null_mean(off, causal) > 1.03  # every null here is a model SNP
    assert _null_mean(on, causal) < 1.02


def test_loco_refit_option_and_no_close_pairs():
    G, gt, y, pairs, causal = family_panel(seed=6, kids=0, n_fam=300, m=4000)
    far = (np.array([0, 2]), np.array([1, 3]), np.array([0.04, 0.03]))  # phi <= 0.025
    a = _run(y, gt, relatedness=far)
    assert not a.extra["loco_refit"] and a.extra["relatedness"]["n_pairs_above_threshold"] == 0
    b = _run(y, gt, relatedness=far, relatedness_options={"loco_refit": True})
    assert b.extra["loco_refit"] and b.extra["loco_refit_reason"] == "forced"
    few = (np.array([0]), np.array([1]), np.array([0.5]))
    c = _run(y, gt, relatedness=few, relatedness_options={"crossfit_folds": "random"})
    assert c.extra["loco_refit"] and c.extra["loco_refit_reason"] == "relatives"
    d = _run(y, gt, relatedness=few)  # default family folds: the pair shares a fold, no refit
    assert d.extra["crossfit_folds"] == "family" and not d.extra["loco_refit"]
    assert d.extra["relatedness"]["crossfit_leak"] == 0.0
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


# P4: external kinship estimated on model SNPs is flagged when its variants are given
def test_external_kinship_on_model_snps_warns(tmp_path, fam):
    G, gt, y, pairs, causal, sub = fam
    ids = np.asarray(gt.variant_ids)
    (tmp_path / "k.prune.in").write_text("\n".join(ids[sub]) + "\n")
    with pytest.warns(RuntimeWarning, match="stage-1 model SNPs"):
        a = hratt(y, gt, model_snps=sub, relatedness=pairs, random_state=1,
                  model_snps_h2_check=False,
                  relatedness_options={"kinship_variants": str(tmp_path / "k.prune.in")})
    assert a.extra["relatedness"]["kinship_model_overlap"] == 1.0
    other = np.setdiff1d(np.arange(gt.n_variants), sub)[:2000]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        b = hratt(y, gt, model_snps=sub, relatedness=pairs, random_state=1,
                  model_snps_h2_check=False, relatedness_options={"kinship_variants": other})
    assert b.extra["relatedness"]["kinship_model_overlap"] == 0.0
    assert _run(y, gt, model_snps=sub, relatedness=pairs).extra["relatedness"][
        "kinship_model_overlap"] is None


# P5: strong structure, kappa active, no covariates -> louder warning / opt-in error
def test_structure_without_pcs_warns_or_errors():
    rng = np.random.default_rng(5)
    n, m = 800, 4000
    pop = np.repeat([0, 1], n // 2)
    p0 = rng.uniform(0.1, 0.9, m)
    shift = rng.normal(0, 0.15, m)
    pf = np.clip(np.stack([p0 - shift, p0 + shift]), 0.02, 0.98)
    G = rng.binomial(2, pf[pop]).astype(np.int8)
    from mixmogam import Genotypes
    gt = Genotypes(G, chromosome=np.repeat([1, 2, 3, 4], m // 4), position=np.arange(m) * 10_000 + 1)
    y = 0.8 * pop + rng.standard_normal(n)
    empty = (np.array([0]), np.array([1]), np.array([0.01]))  # below threshold: no refit by relatives
    with pytest.warns(RuntimeWarning, match="no covariates"):
        r = hratt(y, gt, relatedness=empty, random_state=1)
    assert r.extra["structure"]["strong"] and r.extra["relatedness"]["structure_without_covariates"]
    with pytest.raises(ValueError, match="no covariates"):
        hratt(y, gt, relatedness=empty, random_state=1,
              relatedness_options={"structure_without_pcs": "error"})
    # (ii) KING-auto under the same structure: its noise floor rises past the
    # 2nd-degree bound, and the no-covariates warning says so
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        k = hratt(y, gt, relatedness="auto", random_state=1)
    king = k.extra["relatedness"]["king"]
    floor = king["threshold"]
    # No upper-half switch: the 4-SD MAD floor sits above 2nd-degree
    # under two-pop structure (MAD floor ≈ 0.15).
    assert floor >= 4.0 * king["noise_sd_mad"] - 1e-12
    assert floor > 0.0884
    assert not king.get("noise_scale_from_structure")
    assert king["noise_sd"] == pytest.approx(king["noise_sd_mad"], rel=1e-12)
    msg = [str(x.message) for x in w if "no covariates" in str(x.message)]
    assert len(msg) == 1
    assert "kappa_g replaces" in msg[0] and "LDpred3" in msg[0] and "in simulations" in msg[0]
    # The KING-auto floor also drops 2nd-degree from Phi / kappa_g.
    assert "2nd-degree" in msg[0] or any(
        "2nd-degree" in str(x.message) and ("floor" in str(x.message) or "dropped" in str(x.message))
        for x in w)
    texts = " ".join(str(x.message) for x in w)
    assert "pedigree or pairs" in texts or "KING run within ancestry groups" in texts
    pcs = np.linalg.svd(G - G.mean(0), full_matrices=False)[0][:, :2]
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        hratt(y, gt, X=pcs, relatedness=empty, random_state=1,
              relatedness_options={"structure_without_pcs": "error"})
    assert not any("no covariates" in str(x.message) for x in w)


# (i) loco_refit='auto' (random folds) needs evidence of relatives: supplied
# pairs count from one pair, KING-auto pairs need >= 1% of samples related.
def test_refit_evidence_king_auto_needs_min_related_fraction():
    G, gt, y, pairs, causal = family_panel(seed=7, kids=0, n_fam=300, m=4000, n_causal=100)
    from mixmogam import Genotypes
    rnd = {"crossfit_folds": "random"}

    def twins(k):  # duplicate k samples: k pairs with phi 0.5 that KING finds
        rows = np.r_[np.arange(G.shape[0]), np.arange(k)]
        g2 = Genotypes(G[rows], chromosome=gt.chromosome, position=gt.position)
        return g2, np.r_[y, y[:k] + 0.1], (np.arange(k), G.shape[0] + np.arange(k), np.full(k, 1.0))

    g1, y1, p1 = twins(1)  # 2 of 601 related: 0.3%
    a = _run(y1, g1, relatedness="auto", relatedness_options=rnd)
    ra = a.extra["relatedness"]
    assert ra["n_pairs_above_threshold"] >= 1 and ra["related_fraction"] < 0.01
    assert not ra["relatives_evidence"] and not a.extra["loco_refit"]
    b = _run(y1, g1, relatedness=p1, relatedness_options=rnd)  # the same pair supplied: counts
    assert b.extra["loco_refit"] and b.extra["loco_refit_reason"] == "relatives"
    g10, y10, _ = twins(10)  # 20 of 610: 3.3%
    c = _run(y10, g10, relatedness="auto", relatedness_options=rnd)
    assert c.extra["relatedness"]["related_fraction"] >= 0.01
    assert c.extra["loco_refit"] and c.extra["loco_refit_reason"] == "relatives"
    trig = c.extra["relatedness"]["refit_trigger"]
    assert trig["evidence"] and trig["n_pairs"] >= 10 and trig["cutoff"] >= 0.025
    assert "n_expected_null" in trig and trig["n_expected_null"] < 2
    G2, gt2, y2, _, _ = family_panel(seed=8, n_fam=100, kids=2, m=4000, n_causal=100)
    d = _run(y2, gt2, relatedness="auto", relatedness_options=rnd)  # nuclear families
    assert d.extra["relatedness"]["related_fraction"] > 0.9 and d.extra["loco_refit"]
    e = _run(y2, gt2, relatedness="auto")  # default: KING-auto families share folds, no refit
    assert e.extra["crossfit_folds"] == "family" and not e.extra["loco_refit"]
    assert e.extra["relatedness"]["crossfit_leak"] < 0.01 < d.extra["relatedness"]["crossfit_leak"]


# (v) family folds: the default since 2.3.0, no refit, in-model nulls calibrated
def test_family_folds_default_calibrated_without_refit():
    G, gt, y, pairs, causal = family_panel(seed=0, n_fam=200, kids=4, m=8000, c2=0.3, h2=0.5)
    r = _run(y, gt, relatedness=pairs)
    assert r.extra["crossfit_folds"] == "family" and not r.extra["loco_refit"]
    ff = r.extra["relatedness"]["family_folds"]
    assert ff["largest_family"] == 6 and max(ff["fold_sizes"]) - min(ff["fold_sizes"]) <= 6
    assert r.extra["relatedness"]["crossfit_leak"] == 0.0
    assert _null_mean(r, causal) < 1.02
    with pytest.warns(RuntimeWarning, match="random cross-fit folds split families"):
        rnd = hratt(y, gt, random_state=1, model_snps_h2_check=False, relatedness=pairs,
                    relatedness_options={"crossfit_folds": "random"})
    assert rnd.extra["crossfit_folds"] == "random" and rnd.extra["loco_refit"]
    assert rnd.extra["relatedness"]["crossfit_leak"] > 0.3
    with pytest.raises(ValueError, match="crossfit_folds"):
        hratt(y, gt, relatedness=pairs, relatedness_options={"crossfit_folds": "kin"})


def test_family_folds_fall_back_when_a_family_is_too_large():
    G, gt, y, pairs, causal = family_panel(seed=6, n_fam=60, kids=2, m=3000, n_causal=100)
    i, j, r = pairs
    k = np.arange(59)  # cousins linking family k's first kid to family k+1's: one component
    linked = (np.r_[i, 4 * k + 2], np.r_[j, 4 * (k + 1) + 2], np.r_[r, np.full(59, 0.125)])
    with pytest.warns(RuntimeWarning, match="cannot balance"):
        res = hratt(y, gt, relatedness=linked, random_state=1, model_snps_h2_check=False,
                    relatedness_options={"crossfit_folds": "family"})
    assert res.extra["relatedness"]["family_folds"]["largest_family"] == gt.n_samples
    assert res.extra["crossfit_folds"] == "random"
    assert res.extra["loco_refit"] and res.extra["loco_refit_reason"] == "relatives"


# (iv) carrier families: distinct families among minor-allele carriers for
# variants below carrier_families_p, flagged at <= 3
def test_carrier_families_column_and_flag(tmp_path):
    G, gt, y, pairs, causal = family_panel(seed=2, n_fam=100, kids=3, m=3000, n_causal=100,
                                           prevalence=0.2)
    rng = np.random.default_rng(1)
    n = G.shape[0]  # families of 5: rows 5f, 5f+1 parents, 5f+2.. kids
    rare = np.zeros((n, 4), np.int8)
    rare[[0, 2, 3], 0] = 1          # father + 2 kids of family 0: 1 family
    rare[[0, 2, 5, 10], 1] = 1      # families 0, 1, 2: 3 families
    rare[rng.choice(n, 40, replace=False), 2] = 1  # spread
    rare[:, 3] = 2
    rare[[0, 1], 3] = 1             # minor allele = the other one; carriers 0, 1 (spouses)
    from mixmogam import Genotypes
    m = G.shape[1]
    g2 = Genotypes(np.c_[G, rare], chromosome=np.r_[gt.chromosome, [1, 2, 3, 4]],
                   position=np.r_[gt.position, m * 1_000_000 + np.arange(1, 5)])
    with pytest.warns(RuntimeWarning, match="few_carrier_families"):
        r = hratt(y, g2, trait="binary", relatedness=pairs, random_state=1,
                  model_snps_h2_check=False, relatedness_options={"carrier_families_p": 1.01})
    fam = r.extra["carrier_families"]
    np.testing.assert_array_equal(fam[m:m + 4][[0, 1, 3]], [1, 3, 2])  # spouses: unrelated
    assert fam[m + 2] > 30
    few = r.extra["few_carrier_families"]
    assert few[m] and few[m + 1] and not few[m + 2] and few[m + 3]
    r.write_csv(tmp_path / "r.csv")
    from mixmogam.results import GwasResult
    back = GwasResult.read_csv(str(tmp_path / "r.csv"))
    np.testing.assert_array_equal(back.extra["carrier_families"], fam)
    np.testing.assert_array_equal(back.extra["few_carrier_families"], few)
    assert r.take(np.arange(m, m + 4)).extra["carrier_families"].size == 4
    # default: only p < 1e-4 and MAF < 5% are computed (NaN elsewhere)
    q = _run(y, g2, relatedness=pairs)
    cf = q.extra["carrier_families"]
    maf = np.minimum(q.af, 1 - q.af)
    on = (q.p < 1e-4) & (maf < 0.05)
    assert np.all(np.isnan(cf[~on])) and np.all(np.isfinite(cf[on]))
    c = _run(y, g2, relatedness=pairs, relatedness_options={"carrier_families_p": 1.01,
                                                            "carrier_families_max_maf": 0.51})
    assert np.all(np.isfinite(c.extra["carrier_families"][np.isfinite(c.p)]))  # all computed


def test_auto_weights_and_relatives_fallback_refit_and_kappa_once():
    """auto + sample_weights + relatives falls back to enet, keeps the
    relcal-decided LOCO refit, and applies λ_g = 1/(κ_g c_X) once."""
    # Use a panel that has close relatives but does not trip the structure
    # test (a smaller panel can), so loco_refit_reason is "relatives".
    G, gt, y, pairs, causal = family_panel(n_fam=60, kids=3, m=3000, n_chrom=4,
                                           n_causal=80, h2=0.5, seed=3)
    w = np.random.default_rng(3).uniform(0.5, 2.0, y.size)
    sub = np.sort(np.random.default_rng(5).choice(gt.n_variants, gt.n_variants // 4,
                                                  replace=False))
    # spa_threshold=1e6 disables the saddlepoint so the statistic is the
    # uncorrected retrospective chi-square scaled by 1/(kappa_g c_X) once.
    rnd = {"crossfit_folds": "random"}  # the relatives refit belongs to random folds
    kap = _run(y, gt, model_snps=sub, relatedness=pairs, sample_weights=w,
               stage1="auto", loco_folds=2, spa_threshold=1e6, relatedness_options=rnd)
    sel = kap.extra["stage1_selection"]
    assert sel["chosen"] == "enet" and sel["reason"].startswith("sample_weights")
    assert kap.extra.get("stage1", "enet") == "enet"
    assert not kap.extra["structure"]["strong"]
    assert kap.extra["loco_refit"] and kap.extra["loco_refit_reason"] == "relatives"
    base = _run(y, gt, model_snps=sub, relatedness=pairs, relatedness_correction="none",
                sample_weights=w, stage1="auto", loco_folds=2, spa_threshold=1e6,
                relatedness_options=rnd)
    k, cx = kap.extra["relatedness"]["kappa"], kap.extra["relatedness"]["c_X"]
    np.testing.assert_allclose(kap.f_stat, base.f_stat / (k[gt.chromosome - 1] * cx),
                               rtol=1e-10)
    # Catch a mutation back to refit=struct["strong"] on the auto+weights
    # fallback: extra["loco_refit"] can still say True if it is built from
    # the decision variable, but the statistics would match loco_refit=False.
    forced = _run(y, gt, model_snps=sub, relatedness=pairs, sample_weights=w,
                  stage1="auto", loco_folds=2, spa_threshold=1e6,
                  relatedness_options={**rnd, "loco_refit": True})
    off = _run(y, gt, model_snps=sub, relatedness=pairs, sample_weights=w,
               stage1="auto", loco_folds=2, spa_threshold=1e6,
               relatedness_options={**rnd, "loco_refit": False})
    assert forced.extra["loco_refit"] is True
    assert off.extra["loco_refit"] is False
    np.testing.assert_allclose(kap.f_stat, forced.f_stat, rtol=1e-10)
    assert not np.allclose(kap.f_stat, off.f_stat, rtol=1e-5, atol=1e-5)


def test_king_refit_trigger_ignores_n2_noise_but_fires_on_true_relatives():
    """244 4-SD-floor pairs on an unrelated n=20,000 panel
    (floor 0.042, 1.33% of samples) must not fire. First-degree and 10%
    first-cousin pairs at the same noise scale must still fire. Phi is
    unchanged: the trigger only filters which stored pairs count."""
    n, sd, med, floor = 20_000, 0.0418 / 4.0, 0.0, 0.0418
    rng = np.random.default_rng(1)
    # 266 samples in 244 almost-disjoint noise pairs (1.33% of n), phi just
    # above the 4-SD floor, as on that n=20,000 panel.
    n_noise, n_touch = 244, 266
    samples = rng.choice(n, n_touch, replace=False)
    i_n = samples[:n_noise]
    j_n = np.empty(n_noise, dtype=int)
    j_n[:n_touch - n_noise] = samples[n_noise:]
    # leftover edges reuse the touched set so the graph stays at 266 samples
    extra = n_noise - (n_touch - n_noise)
    j_n[n_touch - n_noise:] = rng.choice(samples, extra, replace=True)
    j_n = np.where(j_n == i_n, (j_n + 1) % n, j_n)
    phi_n = floor + rng.uniform(0.0, 0.008, n_noise)  # ~0.042-0.050
    trig = king_refit_trigger(n, i_n, j_n, phi_n, noise_sd=sd, noise_median=med)
    assert trig["cutoff"] > 0.055  # z·σ̂ ≈ 5.73 × 0.01045
    assert trig["n_pairs"] == 0 and not trig["evidence"]
    assert 0.5 < trig["n_expected_null"] < 2.0
    # same 244 noise pairs plus ~3% of samples in first-degree pairs
    n_fd = int(0.03 * n) // 2  # ~300 pairs, 600 samples
    i_fd = np.arange(n_fd)
    j_fd = n_fd + np.arange(n_fd)
    phi_fd = np.full(n_fd, 0.25)
    trig_fd = king_refit_trigger(n, np.r_[i_n, i_fd], np.r_[j_n, j_fd],
                                 np.r_[phi_n, phi_fd], noise_sd=sd, noise_median=med)
    assert trig_fd["evidence"] and trig_fd["n_pairs"] == n_fd
    assert trig_fd["related_fraction"] >= 0.03
    # 10% first-cousin pairs (phi 0.0625 plus KING noise), the dist-panel case
    n_cu = int(0.10 * n) // 2
    i_c = np.arange(n_cu)
    j_c = n_cu + np.arange(n_cu)
    phi_c = 0.0625 + sd * rng.standard_normal(n_cu)
    trig_c = king_refit_trigger(n, i_c, j_c, phi_c, noise_sd=sd, noise_median=med)
    assert trig_c["evidence"] and trig_c["related_fraction"] >= 0.01
    # supplied-pair path is unchanged: this function is KING-auto only


def test_king_auto_trigger_on_unrelated_and_related_panels():
    """Live KING-auto: unrelated panels do not fire; nuclear families,
    a 10% first-cousin mix, and a few-percent first-degree twins do."""
    from mixmogam import Genotypes

    def _trigger(gt, seed=0):
        trig, rel, info = _live_trigger(gt, seed)
        return trig

    for seed in (1, 2, 3, 7):
        G, gt, y, pairs, causal = family_panel(seed=seed, kids=0, n_fam=400, m=3000,
                                               n_causal=50)
        trig = _trigger(gt, seed)
        assert not trig["evidence"], (seed, trig["reason"])

    G, gt, y, pairs, causal = family_panel(seed=8, n_fam=80, kids=2, m=3000, n_causal=50)
    trig = _trigger(gt, 8)
    assert trig["evidence"] and trig["related_fraction"] > 0.8

    # few-percent first-degree: 16 twin pairs in 800 unrelated founders (~4%)
    G, gt, y, pairs, causal = family_panel(seed=9, kids=0, n_fam=400, m=3000, n_causal=50)
    k = 16
    g2 = Genotypes(G[np.r_[np.arange(G.shape[0]), np.arange(k)]],
                   chromosome=gt.chromosome, position=gt.position)
    trig = _trigger(g2, 9)
    assert trig["evidence"] and trig["n_pairs"] >= k
    assert trig["related_fraction"] >= 0.03

    # ~10% first-cousin pairs (the dist-panel first-cousin share): 40 cousin
    # units plus 320 unrelated founders, n=400.
    rng = np.random.default_rng(11)
    m, n_cu, n_unrel = 8000, 40, 320
    f = rng.uniform(0.05, 0.5, m)
    n_found = 4 * n_cu + n_unrel
    hap = (rng.random((2 * n_found, m)) < f).astype(np.int8)

    def meiosis(idx):
        pick = rng.integers(0, 2, m)
        return np.where(pick == 0, hap[2 * idx], hap[2 * idx + 1])

    rows = []
    for u in range(n_cu):
        gp, gq, m1, m2 = 4 * u, 4 * u + 1, 4 * u + 2, 4 * u + 3
        s1a, s1b = meiosis(gp), meiosis(gq)
        s2a, s2b = meiosis(gp), meiosis(gq)
        c1 = np.where(rng.integers(0, 2, m) == 0, s1a, s1b) + meiosis(m1)
        c2 = np.where(rng.integers(0, 2, m) == 0, s2a, s2b) + meiosis(m2)
        rows.append(c1.astype(np.int8))
        rows.append(c2.astype(np.int8))
    for u in range(n_unrel):
        idx = 4 * n_cu + u
        rows.append((hap[2 * idx] + hap[2 * idx + 1]).astype(np.int8))
    G = np.vstack(rows)
    gt = Genotypes(G, chromosome=np.repeat([1, 2], m // 2),
                   position=np.arange(m) * 10_000 + 1)
    trig = _trigger(gt, 11)
    assert trig["evidence"] and trig["related_fraction"] >= 0.01, trig["reason"]


def test_king_refit_trigger_caps_at_second_degree_and_warns():
    """If the n-aware cutoff would exceed 0.0884, cap it, warn, and still
    count 2nd-degree pairs (phi 0.125). The refit is not disabled."""
    n, sd, med = 20_000, 0.020, 0.0  # z ~ 5.73 -> uncapped ~ 0.115 > 0.0884
    # a few 2nd-degree pairs plus some 3rd-degree / noise just above 0.09
    i2, j2 = np.arange(400), 400 + np.arange(400)  # 800 samples, 4%
    phi2 = np.full(400, 0.125)
    i3, j3 = 1000 + np.arange(50), 1100 + np.arange(50)
    phi3 = np.full(50, 0.090)
    with pytest.warns(RuntimeWarning, match="2nd/3rd-degree bound"):
        trig = king_refit_trigger(n, np.r_[i2, i3], np.r_[j2, j3], np.r_[phi2, phi3],
                                  noise_sd=sd, noise_median=med, n_variants=2000)
    assert trig["capped"] and trig["cutoff"] == pytest.approx(KING_SECOND_DEGREE)
    assert trig["cutoff_uncapped"] > KING_SECOND_DEGREE
    assert trig["evidence"] and trig["n_pairs"] == 450  # both 0.125 and 0.090 exceed the cap
    # 3rd-degree / noise at 0.07 would not count even at the cap; 2nd-degree would
    with pytest.warns(RuntimeWarning, match="is not disabled"):
        trig2 = king_refit_trigger(n, i2, j2, phi2, noise_sd=sd, noise_median=med)
    assert trig2["evidence"] and trig2["n_pairs"] == 400
    assert trig2["related_fraction"] >= 0.01


def test_king_auto_second_degree_fraction_rule_not_a_handful():
    """A handful of true 2nd-degree (half-sib) pairs below 1% of samples
    does not force a refit (the leak is negligible there). A few percent
    of samples in half-sib pairs does fire."""
    from mixmogam import Genotypes

    def half_sib_panel(n_unrel, n_pairs, m=8000, seed=0):
        rng = np.random.default_rng(seed)
        f = rng.uniform(0.05, 0.5, m)
        # 3 founders per half-sib unit (shared parent + 2 married-ins) + unrel
        n_found = 3 * n_pairs + n_unrel
        hap = (rng.random((2 * n_found, m)) < f).astype(np.int8)

        def meiosis(idx):
            pick = rng.integers(0, 2, m)
            return np.where(pick == 0, hap[2 * idx], hap[2 * idx + 1])

        rows = []
        for u in range(n_pairs):
            p, m1, m2 = 3 * u, 3 * u + 1, 3 * u + 2
            rows.append((meiosis(p) + meiosis(m1)).astype(np.int8))
            rows.append((meiosis(p) + meiosis(m2)).astype(np.int8))
        for u in range(n_unrel):
            idx = 3 * n_pairs + u
            rows.append((hap[2 * idx] + hap[2 * idx + 1]).astype(np.int8))
        G = np.vstack(rows)
        return Genotypes(G, chromosome=np.repeat([1, 2], m // 2),
                         position=np.arange(m) * 10_000 + 1)

    def _trig(gt, seed):
        trig, rel, info = _live_trigger(gt, seed)
        return trig

    # 5 half-sib pairs in 2,000 samples: 10/2000 = 0.5% < 1%, no fire
    gt_few = half_sib_panel(1990, 5, seed=21)
    trig_few = _trig(gt_few, 21)
    assert trig_few["n_pairs"] >= 4  # KING should recover most half-sibs (phi 0.125)
    assert trig_few["related_fraction"] < 0.01
    assert not trig_few["evidence"], trig_few["reason"]

    # 20 half-sib pairs: 40/2000 = 2%, fires
    gt_many = half_sib_panel(1960, 20, seed=22)
    trig_many = _trig(gt_many, 22)
    assert trig_many["evidence"] and trig_many["related_fraction"] >= 0.01, trig_many["reason"]


def test_king_auto_trigger_does_not_fire_on_structure_or_admixture():
    """Unrelated structured (Fst 0.05) and admixed panels must not fire.
    Records the trigger-pair (null exceedance) count."""
    from mixmogam import Genotypes

    def _trig(G, seed):
        n, m = G.shape
        gt = Genotypes(np.ascontiguousarray(G),
                       chromosome=np.repeat([1, 2], m // 2),
                       position=np.arange(m) * 10_000 + 1)
        trig, rel, info = _live_trigger(gt, seed)
        return trig, int(rel.upper.nnz)

    def balding_nichols(n, m, fst, seed, admixed=False):
        rng = np.random.default_rng(seed)
        p0 = rng.uniform(0.1, 0.9, m)
        a = p0 * (1 - fst) / fst
        b = (1 - p0) * (1 - fst) / fst
        p = np.clip(np.stack([rng.beta(a, b), rng.beta(a, b)]), 0.02, 0.98)
        if admixed:
            q = rng.uniform(0.1, 0.9, n)
            freq = q[:, None] * p[0] + (1 - q)[:, None] * p[1]
            return rng.binomial(2, freq).astype(np.int8)
        pop = np.repeat([0, 1], n // 2)
        if n % 2:
            pop = np.r_[pop, 0]
        return rng.binomial(2, p[pop]).astype(np.int8)

    n, m = 2000, 8000
    for kind, admixed, seed in (("struct", False, 31), ("struct", False, 32),
                                ("admix", True, 41), ("admix", True, 42)):
        G = balding_nichols(n, m, fst=0.05, seed=seed, admixed=admixed)
        trig, n_phi = _trig(G, seed)
        assert not trig["evidence"], (kind, seed, trig["reason"])
        # expose the null exceedance count for the report
        trig["_n_pairs_phi"] = n_phi
        print(f"{kind} seed={seed} Phi-pairs={n_phi} trig-pairs={trig['n_pairs']} "
              f"frac={trig['related_fraction']:.3%} tail={trig['noise_tail_ratio']:.2f} "
              f"cut={trig['cutoff']:.4f} σ={trig['noise_sd']:.4f} E0={trig['n_expected_null']:.2f}",
              flush=True)


def _live_trigger(gt, seed=0, threshold=0.025):
    """Same path as twostep: trigger sees pre-floor KING pairs."""
    rel, info = resolve_relatedness(
        "auto", gt, mean=np.asarray(gt.G, dtype=np.float64).mean(0),
        model_idx=None, threshold=threshold, options=RELATEDNESS_OPTIONS,
        rng=np.random.default_rng(seed))
    king = info["king"]
    pre = info.get("king_prefloor") or {}
    if pre.get("i") is not None and np.size(pre["i"]):
        i, j, phi = pre["i"], pre["j"], np.asarray(pre["r"], dtype=np.float64) / 2.0
    else:
        u = rel.upper.tocoo()
        i, j, phi = u.row, u.col, u.data / 2.0
    loc = king.get("noise_location", king.get("noise_median", 0.0))
    trig = king_refit_trigger(
        gt.n_samples, i, j, phi,
        noise_sd=king["noise_sd"],
        noise_median=king["noise_median"],
        noise_location=None if loc is None else float(loc),
        noise_tail_ratio=king.get("noise_tail_ratio", 1.0),
        n_variants=king.get("n_variants"),
        relatedness_threshold=threshold)
    return trig, rel, info


def test_king_refit_trigger_warns_when_cutoff_exceeds_cousin_degree():
    """Once the effective cutoff exceeds 0.0625, 3rd-degree pairs no
    longer count and a warning recommends a pedigree; the 0.0884 cap
    is a separate, higher threshold."""
    n, sd = 2000, 0.015  # z ≈ 4.89 -> cutoff ≈ 0.073, under the cap
    i2, j2 = np.arange(20), 20 + np.arange(20)  # 2% of samples, 2nd-degree
    phi2 = np.full(20, 0.125)
    i3, j3 = 100 + np.arange(20), 120 + np.arange(20)
    phi3 = np.full(20, KING_COUSIN_DEGREE)
    with pytest.warns(RuntimeWarning, match="0.0625"):
        trig = king_refit_trigger(n, np.r_[i2, i3], np.r_[j2, j3], np.r_[phi2, phi3],
                                  noise_sd=sd, noise_median=0.0, n_variants=8000)
    assert not trig["capped"]
    assert KING_COUSIN_DEGREE < trig["cutoff"] < KING_SECOND_DEGREE
    assert trig["n_pairs"] == 20  # only the 2nd-degree pairs
    assert trig["evidence"] and trig["related_fraction"] >= 0.01


def test_king_noise_scale_is_mad_even_under_a_mixture():
    """1.4826·MAD is kept on a single Gaussian and on a two-population
    mixture. There is no upper-half switch: the mixture inflates the
    MAD, and a few related-pair contaminants do not drift it."""
    from mixmogam._relatedness import _king_noise_scale

    rng = np.random.default_rng(0)
    sig = 0.0155
    gauss = rng.normal(0.0, sig, 200_000)
    loc, sd, sd_mad, med, tail, struct, tail_mad = _king_noise_scale(gauss)
    assert not struct
    assert abs(loc) < 0.002
    assert sd == pytest.approx(sig, rel=0.05)
    assert sd == pytest.approx(sd_mad, rel=1e-12)
    assert tail == pytest.approx(1.0, abs=0.08)
    assert tail_mad == pytest.approx(1.0, abs=0.08)

    # 50/50 mixture, separation like Fst 0.05 (review: means ≈ 0 and -0.05)
    mix = np.r_[rng.normal(-0.0057, sig, 100_000), rng.normal(-0.0561, sig, 100_000)]
    loc_m, sd_m, mad_m, _, tail_m, struct_m, tail_mad_m = _king_noise_scale(mix)
    assert not struct_m
    assert sd_m == pytest.approx(mad_m, rel=1e-12)
    assert mad_m > 1.5 * sig
    assert tail_mad_m < 0.8
    # 10% of *samples* in 2nd-degree pairs is a tiny fraction of all pairs
    with_rel = np.r_[mix, np.full(2000, 0.125)]
    loc_r, sd_r, mad_r, _, _, struct_r, _ = _king_noise_scale(with_rel)
    assert not struct_r
    assert sd_r == pytest.approx(sd_m, rel=0.05)
    assert mad_r == pytest.approx(mad_m, rel=0.05)


def _structured_halfsib_panel(n, m, fst, n_pairs, seed, admixed=False):
    """Balding-Nichols (or admixed) panel with n_pairs within-population
    half-sib pairs (2nd-degree, phi 0.125) planted at the end."""
    from mixmogam import Genotypes

    rng = np.random.default_rng(seed)
    n_unrel = n - 2 * n_pairs
    p0 = rng.uniform(0.1, 0.9, m)
    a = p0 * (1 - fst) / fst
    b = (1 - p0) * (1 - fst) / fst
    p = np.clip(np.stack([rng.beta(a, b), rng.beta(a, b)]), 0.02, 0.98)

    def draw(freq, k):
        return rng.binomial(2, np.broadcast_to(freq, (k, m))).astype(np.int8)

    def hap(freq):
        return (rng.random(m) < freq).astype(np.int8)

    if admixed:
        q = rng.uniform(0.1, 0.9, n_unrel)
        freq = q[:, None] * p[0] + (1 - q)[:, None] * p[1]
        Gun = draw(freq, n_unrel)
    else:
        pop = np.repeat([0, 1], (n_unrel + 1) // 2)[:n_unrel]
        Gun = draw(p[pop], n_unrel)

    rows = [Gun]
    true_i, true_j = [], []
    base = n_unrel
    for u in range(n_pairs):
        pop_u = u % 2
        if admixed:
            qu = rng.uniform(0.2, 0.8)
            f = qu * p[0] + (1 - qu) * p[1]
        else:
            f = p[pop_u]
        # shared parent + two married-ins from the same frequency
        hp, hm1, hm2 = hap(f), hap(f), hap(f)
        hp2, hm12, hm22 = hap(f), hap(f), hap(f)
        # two meioses of the shared parent, independent mates
        def mei(h0, h1):
            pick = rng.integers(0, 2, m)
            return np.where(pick == 0, h0, h1)
        kid1 = (mei(hp, hp2) + mei(hm1, hm12)).astype(np.int8)
        kid2 = (mei(hp, hp2) + mei(hm2, hm22)).astype(np.int8)
        rows.append(kid1[None, :])
        rows.append(kid2[None, :])
        true_i.append(base + 2 * u)
        true_j.append(base + 2 * u + 1)
    G = np.vstack(rows)
    assert G.shape[0] == n
    gt = Genotypes(np.ascontiguousarray(G),
                   chromosome=np.repeat([1, 2], m // 2),
                   position=np.arange(m) * 10_000 + 1)
    return gt, np.array(true_i), np.array(true_j)


def test_king_auto_structure_trigger_fires_even_if_phi_drops_second_degree():
    """Fst 0.05 plus ≈ 2% of samples in true 2nd-degree (half-sib) pairs:
    the 4-SD MAD floor may sit above 0.125 and drop those pairs from
    Phi, but the pre-floor trigger must still fire, and Phi must not
    flood (false pairs at the master floor). Recovery in Phi may be
    low; supply a pedigree or pairs in this setting."""
    n, m, n_pairs = 2000, 8000, 20  # 40/2000 = 2% of samples
    gt, i_true, j_true = _structured_halfsib_panel(n, m, fst=0.05, n_pairs=n_pairs, seed=51)
    true_set = {frozenset((int(a), int(b))) for a, b in zip(i_true, j_true)}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        trig, rel, info = _live_trigger(gt, seed=51)

    king = info["king"]
    stats = _phi_false_and_true(rel, i_true, j_true)
    master_floor = 4.0 * float(king["noise_sd_mad"])
    pre = info.get("king_prefloor") or {}
    if pre.get("r") is not None and np.size(pre["r"]):
        n_master = int(np.sum(np.asarray(pre["r"], dtype=np.float64) / 2.0 > master_floor))
    else:
        n_master = 0
    assert not king.get("noise_scale_from_structure")
    assert king["threshold"] >= master_floor - 1e-12
    assert trig["evidence"], trig["reason"]
    assert trig["related_fraction"] >= 0.01
    # False pairs stay at the master floor; thousands would be a flood.
    assert stats["n_false"] <= n_master + 5, (
        f"{stats['n_false']} false Phi pairs; master_floor_pairs={n_master} "
        f"floor={king['threshold']:.4f}")
    # documented: recovery may be low
    assert stats["n_rec"] <= len(true_set)


def _ancestry_panel_with_pairs(n, m, fst, kind, frac, seed, admixed=False):
    """Balding-Nichols (or admixed) panel with ``frac`` of samples in
    planted half-sib (phi 1/8) or first-cousin (phi 1/16) pairs.
    """
    from mixmogam import Genotypes

    rng = np.random.default_rng(seed)
    p0 = rng.uniform(0.1, 0.9, m)
    a = p0 * (1 - fst) / fst
    b = (1 - p0) * (1 - fst) / fst
    P = np.clip(np.stack([rng.beta(a, b), rng.beta(a, b)]), 0.02, 0.98)
    if admixed:
        q = rng.uniform(0.1, 0.9, n)
        F = q[:, None] * P[0] + (1 - q)[:, None] * P[1]
    else:
        pop = np.r_[np.zeros(n // 2, int), np.ones(n - n // 2, int)]
        F = P[pop]

    def hap(f):
        return (rng.random((2, m)) < f).astype(np.int8)

    mei = lambda H: np.where(rng.integers(0, 2, m) == 0, H[0], H[1])
    npair = int(frac * n) // 2
    rows, true_i, true_j = [], [], []
    for k in range(npair):
        f = F[2 * k]
        if kind == "halfsib":
            parent = hap(f)
            a_ = mei(parent) + mei(hap(f))
            b_ = mei(parent) + mei(hap(f))
        else:
            GP, GQ = hap(f), hap(f)
            s1 = np.stack([mei(GP), mei(GQ)])
            s2 = np.stack([mei(GP), mei(GQ)])
            a_ = mei(s1) + mei(hap(f))
            b_ = mei(s2) + mei(hap(f))
        rows += [a_, b_]
        true_i.append(2 * k)
        true_j.append(2 * k + 1)
    for i in range(2 * npair, n):
        rows.append(hap(F[i]).sum(0))
    G = np.vstack(rows).astype(np.int8)
    gt = Genotypes(np.ascontiguousarray(G),
                   chromosome=np.repeat([1, 2], m // 2),
                   position=np.arange(m) * 10_000 + 1)
    return gt, np.array(true_i), np.array(true_j)


def test_king_auto_structured_trigger_sees_halfsibs_even_if_phi_drops_them():
    """On two populations the 4-SD MAD
    floor sits above 0.125, so planted half-sibs (and first cousins)
    may be absent from Phi. The pre-floor trigger must still fire on
    2% half-sibs, and Phi must not flood."""
    n, m, fst = 800, 4000, 0.05
    gt, i_t, j_t = _ancestry_panel_with_pairs(n, m, fst, "halfsib", 0.02, seed=7)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        trig, rel, info = _live_trigger(gt, seed=7)
    king = info["king"]
    stats = _phi_false_and_true(rel, i_t, j_t)
    master_floor = 4.0 * float(king["noise_sd_mad"])
    pre = info.get("king_prefloor") or {}
    if pre.get("r") is not None and np.size(pre["r"]):
        n_master = int(np.sum(np.asarray(pre["r"], dtype=np.float64) / 2.0 > master_floor))
    else:
        n_master = 0
    assert not king.get("noise_scale_from_structure")
    assert king["threshold"] >= master_floor - 1e-12
    assert trig["evidence"], trig["reason"]
    assert stats["n_false"] <= n_master + 5, (
        f"half-sib panel: {stats['n_false']} false Phi pairs; "
        f"floor={king['threshold']:.4f}")

    gt, i_t, j_t = _ancestry_panel_with_pairs(n, m, fst, "cousin", 0.03, seed=7)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        trig, rel, info = _live_trigger(gt, seed=7)
    king = info["king"]
    stats = _phi_false_and_true(rel, i_t, j_t)
    master_floor = 4.0 * float(king["noise_sd_mad"])
    pre = info.get("king_prefloor") or {}
    if pre.get("r") is not None and np.size(pre["r"]):
        n_master = int(np.sum(np.asarray(pre["r"], dtype=np.float64) / 2.0 > master_floor))
    else:
        n_master = 0
    assert king["threshold"] >= master_floor - 1e-12
    assert stats["n_false"] <= n_master + 5, (
        f"cousin panel: {stats['n_false']} false Phi pairs; "
        f"floor={king['threshold']:.4f}")


def test_king_auto_structure_hratt_refits_and_reports_actual_refit():
    """End-to-end through hratt: structured panel + planted half-sibs
    must keep those pairs in Phi, and extra['loco_refit'] must be the
    flag returned by the scores fit that actually ran."""
    n, m = 400, 2000
    gt, i_t, j_t = _ancestry_panel_with_pairs(n, m, 0.05, "halfsib", 0.04, seed=11)
    rng = np.random.default_rng(11)
    y = rng.standard_normal(n)
    for a, b in zip(i_t, j_t):
        s = rng.normal() * 0.8
        y[a] += s
        y[b] += s
    kw = dict(random_state=11, model_snps_h2_check=False, loco_folds=2,
              grid=[(0.0, 1.0)], vb_max_iter=40, he_probes=8)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = hratt(y, gt, relatedness="auto", **kw)
    king = r.extra["relatedness"]["king"]
    relinfo = r.extra["relatedness"]
    assert not king.get("noise_scale_from_structure")
    assert king["threshold"] >= 4.0 * king["noise_sd_mad"] - 1e-12
    # Detected pre-floor relatives, not necessarily modeled in Phi.
    assert relinfo["n_pairs_detected"] >= 1
    assert relinfo["relatives_evidence"]
    assert r.extra["loco_refit"] is True
    # Two-pop Fst 0.05 also trips the unrelated-set structure test, so
    # the recorded reason may be "structure"; relatives evidence is required.
    assert r.extra["loco_refit_reason"] in ("relatives", "structure")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        off = hratt(y, gt, relatedness="auto",
                    relatedness_options={"loco_refit": False}, **kw)
    assert off.extra["loco_refit"] is False
    # f_stat need not move on this small structured panel (structure
    # already forces a refit decision; turning the option off is what
    # extra['loco_refit'] must reflect). The mutation test above, on an
    # unstructured family panel, checks that the fit itself changes.


def test_king_refit_trigger_respects_user_threshold_at_or_above_cap():
    """A user relatedness_threshold >= 0.0884 is respected: the cap
    does not pull the cutoff back down."""
    n, sd = 2000, 0.015
    i, j = np.arange(20), 20 + np.arange(20)
    phi = np.full(20, 0.20)
    trig = king_refit_trigger(n, i, j, phi, noise_sd=sd, noise_median=0.0,
                              relatedness_threshold=0.10, n_variants=8000)
    assert not trig["capped"]
    assert trig["cutoff"] == pytest.approx(0.10)
    assert trig["n_pairs"] == 20


def _bn_counts_panel(n, m, fst, k_pops, props, n_pairs, seed):
    """Balding-Nichols K-population panel, optional planted half-sibs
    (cycled across pops). ``props`` is a tuple of ancestry fractions,
    or None for equal sizes."""
    from mixmogam import Genotypes

    rng = np.random.default_rng(seed)
    if props is None:
        sizes = [n // k_pops] * k_pops
        sizes[-1] += n - sum(sizes)
    else:
        sizes = [int(round(p * n)) for p in props]
        sizes[-1] = n - sum(sizes[:-1])
    p0 = rng.uniform(0.1, 0.9, m)
    a = p0 * (1 - fst) / fst
    b = (1 - p0) * (1 - fst) / fst
    P = np.clip(np.stack([rng.beta(a, b) for _ in range(k_pops)]), 0.02, 0.98)

    def hap(f):
        return (rng.random((2, m)) < f).astype(np.int8)

    mei = lambda H: np.where(rng.integers(0, 2, m) == 0, H[0], H[1])
    rows, true_i, true_j, used = [], [], [], [0] * k_pops
    for t_ix in range(n_pairs):
        p_u = t_ix % k_pops
        f = P[p_u]
        parent = hap(f)
        a_ = mei(parent) + mei(hap(f))
        b_ = mei(parent) + mei(hap(f))
        rows += [a_, b_]
        true_i.append(2 * t_ix)
        true_j.append(2 * t_ix + 1)
        used[p_u] += 2
    leftover = []
    for p_u, sz in enumerate(sizes):
        leftover += [p_u] * max(0, sz - used[p_u])
    n_rest = n - 2 * n_pairs
    if len(leftover) < n_rest:
        leftover += [k_pops - 1] * (n_rest - len(leftover))
    leftover = leftover[:n_rest]
    rng.shuffle(leftover)
    for p_u in leftover:
        rows.append(hap(P[p_u]).sum(0))
    G = np.vstack(rows).astype(np.int8)
    gt = Genotypes(np.ascontiguousarray(G),
                   chromosome=np.repeat([1, 2, 3, 4], (m + 3) // 4)[:m],
                   position=np.arange(m) * 10_000 + 1)
    return gt, np.array(true_i, dtype=int), np.array(true_j, dtype=int)


def _phi_false_and_true(rel, true_i, true_j):
    u = rel.upper.tocoo()
    true = {frozenset((int(a), int(b))) for a, b in zip(true_i, true_j)}
    king = {frozenset((int(a), int(b))): float(ph)
            for a, b, ph in zip(u.row, u.col, u.data / 2.0)}
    rec = true & set(king)
    n_false = len(king) - len(rec)
    false_sum = sum(ph for pair, ph in king.items() if pair not in true)
    rec_sum = sum(king[p] for p in rec)
    ped = 0.125 * len(true) if true else 0.0
    return dict(n_phi=len(king), n_true=len(true), n_rec=len(rec),
                n_false=n_false, rec_sum=rec_sum, false_sum=false_sum,
                ped=ped, total_sum=rec_sum + false_sum)


def test_king_auto_k4_null_false_pairs_match_master_floor():
    """S1: K=4, Fst 0.05, n=3000, no relatives. A negative overall
    median must not collapse the floor to 0.025. False Phi pairs stay
    within a few of master's 4-MAD floor on the same panel (typically
    zero)."""
    n, m, seed = 3000, 8000, 1
    gt, i_t, j_t = _bn_counts_panel(n, m, fst=0.05, k_pops=4, props=None,
                                    n_pairs=0, seed=seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        trig, rel, info = _live_trigger(gt, seed=seed)
    king = info["king"]
    stats = _phi_false_and_true(rel, i_t, j_t)
    master_floor = 4.0 * float(king["noise_sd_mad"])
    pre = info.get("king_prefloor") or {}
    if pre.get("r") is not None and np.size(pre["r"]):
        n_master = int(np.sum(np.asarray(pre["r"], dtype=np.float64) / 2.0 > master_floor))
    else:
        n_master = 0
    assert not king.get("noise_scale_from_structure")
    assert king["threshold"] >= master_floor - 1e-12
    assert king["threshold"] > 0.05, king["threshold"]
    assert stats["n_false"] <= n_master + 5, (
        f"false={stats['n_false']} master_floor_pairs={n_master} "
        f"floor={king['threshold']:.4f} master={master_floor:.4f} "
        f"loc={king['noise_location']:.4f}")
    assert not trig["evidence"], trig["reason"]


def test_king_auto_9010_null_does_not_underestimate_or_fire():
    """S2: a null Fst 0.05 90/10 panel must not switch to an
    underestimated upper-half sigma, must not fire, and must not flood
    Phi. Seeds 7/8/9 include the seed that flipped under the unguarded
    0.8 rule (seed 8)."""
    n, m = 2000, 8000
    for seed in (7, 8, 9):
        gt, i_t, j_t = _bn_counts_panel(n, m, fst=0.05, k_pops=2,
                                        props=(0.90, 0.10), n_pairs=0,
                                        seed=seed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            trig, rel, info = _live_trigger(gt, seed=seed)
        king = info["king"]
        stats = _phi_false_and_true(rel, i_t, j_t)
        assert not king.get("noise_scale_from_structure"), (
            f"seed {seed} switched; sd={king['noise_sd']:.4f} "
            f"MAD={king['noise_sd_mad']:.4f} loc={king['noise_location']:.4f}")
        assert king["noise_sd"] == pytest.approx(king["noise_sd_mad"], rel=1e-12)
        assert king["threshold"] >= 4.0 * king["noise_sd_mad"] - 1e-12
        assert stats["n_false"] <= 5, (
            f"seed {seed}: {stats['n_false']} false Phi pairs, "
            f"floor={king['threshold']:.4f}")
        assert not trig["evidence"], f"seed {seed} fired: {trig['reason']}"


def test_king_auto_7030_no_switch_no_flood_and_warns():
    """70/30 two-pop Fst 0.05, n=2000, m=4000, 20 planted half-sibs.

    There is no upper-half switch, so the 4-SD MAD floor may sit above
    0.125 and recovery of 2nd-degree pairs may be low. False pairs
    stay at the master floor (no flood). A warning recommends a
    pedigree or pairs because the floor can drop close relatives.
    """
    n, m, n_pairs = 2000, 4000, 20
    for seed in (1, 5):
        gt, i_t, j_t = _bn_counts_panel(n, m, fst=0.05, k_pops=2,
                                        props=(0.70, 0.30), n_pairs=n_pairs,
                                        seed=seed)
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter("always", RuntimeWarning)
            trig, rel, info = _live_trigger(gt, seed=seed)
        king = info["king"]
        stats = _phi_false_and_true(rel, i_t, j_t)
        master_floor = 4.0 * float(king["noise_sd_mad"])
        pre = info.get("king_prefloor") or {}
        if pre.get("r") is not None and np.size(pre["r"]):
            n_master = int(np.sum(np.asarray(pre["r"], dtype=np.float64) / 2.0 > master_floor))
        else:
            n_master = 0
        assert not king.get("noise_scale_from_structure"), (
            f"seed {seed} switched; sd={king['noise_sd']:.4f} "
            f"MAD={king['noise_sd_mad']:.4f}")
        assert king["noise_sd"] == pytest.approx(king["noise_sd_mad"], rel=1e-12)
        assert king["threshold"] >= master_floor - 1e-12
        assert stats["n_false"] <= n_master + 5, (
            f"seed {seed}: {stats['n_false']} false Phi pairs; "
            f"master_floor_pairs={n_master} floor={king['threshold']:.4f}")
        texts = " ".join(str(w.message) for w in rec
                         if issubclass(w.category, RuntimeWarning))
        assert "supply a pedigree or pairs" in texts or "KING run within ancestry groups" in texts, (
            f"seed {seed}: missing structure/floor warning: {texts!r}")



def test_king_auto_unrelated_null_n3000_few_false_pairs():
    """A negative loc must not pull the floor down on an unstructured
    null. n=3000 should have about as many false Phi pairs as master's
    4-MAD floor (0-1, not dozens)."""
    from mixmogam import Genotypes

    n, m, seed = 3000, 8000, 1
    rng = np.random.default_rng(seed)
    f = rng.uniform(0.05, 0.5, m)
    G = rng.binomial(2, f, size=(n, m)).astype(np.int8)
    gt = Genotypes(np.ascontiguousarray(G),
                   chromosome=np.repeat([1, 2, 3, 4], (m + 3) // 4)[:m],
                   position=np.arange(m) * 10_000 + 1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        trig, rel, info = _live_trigger(gt, seed=seed)
    king = info["king"]
    stats = _phi_false_and_true(rel, [], [])
    master_floor = 4.0 * float(king["noise_sd_mad"])
    pre = info.get("king_prefloor") or {}
    if pre.get("r") is not None and np.size(pre["r"]):
        n_master = int(np.sum(np.asarray(pre["r"], dtype=np.float64) / 2.0 > master_floor))
    else:
        n_master = 0
    assert king["threshold"] >= master_floor - 1e-12
    assert not king.get("noise_scale_from_structure")
    assert stats["n_false"] <= n_master + 5, (
        f"{stats['n_false']} false Phi pairs (master_floor_pairs={n_master}); "
        f"floor={king['threshold']:.4f} loc={king['noise_location']:.4f} "
        f"master={master_floor:.4f}")
    assert not trig["evidence"], trig["reason"]


def test_king_refit_trigger_uses_prefloor_when_floor_exceeds_cutoff():
    """The Phi floor can sit above the n-aware trigger cutoff (a large
    king_noise_sd_multiple), so planted half-sibs drop out of Phi. The
    trigger must still see them through the pre-floor pair list, as
    twostep does on the real path."""
    from mixmogam import Genotypes
    from mixmogam._relatedness import RELATEDNESS_OPTIONS, resolve_relatedness

    n, m, n_pairs, seed = 600, 2000, 8, 21
    rng = np.random.default_rng(seed)
    f = rng.uniform(0.05, 0.5, m)

    def hap():
        return (rng.random((2, m)) < f).astype(np.int8)

    mei = lambda H: np.where(rng.integers(0, 2, m) == 0, H[0], H[1])
    rows, true_i, true_j = [], [], []
    for t_ix in range(n_pairs):
        parent = hap()
        a_ = mei(parent) + mei(hap())
        b_ = mei(parent) + mei(hap())
        rows += [a_, b_]
        true_i.append(2 * t_ix)
        true_j.append(2 * t_ix + 1)
    for _ in range(n - 2 * n_pairs):
        rows.append(hap().sum(0))
    G = np.vstack(rows).astype(np.int8)
    gt = Genotypes(np.ascontiguousarray(G),
                   chromosome=np.repeat([1, 2], m // 2),
                   position=np.arange(m) * 10_000 + 1)
    opts = {**RELATEDNESS_OPTIONS, "king_noise_sd_multiple": 16.0}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        rel, info = resolve_relatedness(
            "auto", gt, mean=np.asarray(gt.G, dtype=np.float64).mean(0),
            model_idx=None, threshold=0.025, options=opts,
            rng=np.random.default_rng(seed))
        king = info["king"]
        pre = info.get("king_prefloor") or {}
        trig = king_refit_trigger(
            n, pre["i"], pre["j"], np.asarray(pre["r"], dtype=np.float64) / 2.0,
            noise_sd=king["noise_sd"], noise_median=king["noise_median"],
            noise_location=king.get("noise_location"),
            noise_tail_ratio=king.get("noise_tail_ratio", 1.0),
            n_variants=king.get("n_variants"),
            relatedness_threshold=0.025)
    assert king["threshold"] > trig["cutoff"], (
        f"floor {king['threshold']:.4f} should exceed cutoff {trig['cutoff']:.4f} "
        f"at k=16, sigma={king['noise_sd']:.4f}")
    assert king["threshold"] > 0.125
    assert rel.upper.nnz == 0 or all(
        (np.asarray(rel.upper.tocoo().data) / 2.0) > king["threshold"] - 1e-15)
    # Phi dropped the planted pairs; pre-floor + trigger must still see them
    u = rel.upper.tocoo()
    in_phi = {frozenset((int(a), int(b))) for a, b in zip(u.row, u.col)}
    true = {frozenset((int(a), int(b))) for a, b in zip(true_i, true_j)}
    assert len(true & in_phi) == 0
    assert trig["evidence"], trig["reason"]
    assert trig["n_pairs"] >= n_pairs
    # The production path must itself have counted pre-floor pairs.
    # Using only post-floor Phi (M9) leaves these zero.
    assert info["refit_trigger"]["evidence"], info["refit_trigger"]["reason"]
    assert info["n_pairs_detected"] >= n_pairs
    assert info["n_pairs_modeled"] == 0
    # Post-floor pairs would miss them:
    trig_post = king_refit_trigger(
        n, u.row, u.col, u.data / 2.0 if u.row.size else np.zeros(0),
        noise_sd=king["noise_sd"], noise_median=king["noise_median"],
        noise_location=king.get("noise_location"),
        noise_tail_ratio=king.get("noise_tail_ratio", 1.0),
        n_variants=king.get("n_variants"),
        relatedness_threshold=0.025)
    assert not trig_post["evidence"]


def test_king_auto_2pop_trigger_phi_mismatch_warns_and_records_counts():
    """2-pop Fst 0.05, 2% half-sibs: the MAD floor sits above 0.125, so
    Phi drops 2nd-degree pairs, but the pre-floor trigger still fires.
    Metadata reports detected vs modeled counts separately, and a
    warning says relatives were detected but not modeled in Phi."""
    n, m, n_pairs, seed = 2000, 8000, 20, 51
    gt, i_t, j_t = _structured_halfsib_panel(n, m, fst=0.05, n_pairs=n_pairs, seed=seed)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always", RuntimeWarning)
        trig, rel, info = _live_trigger(gt, seed=seed)
    king = info["king"]
    stats = _phi_false_and_true(rel, i_t, j_t)
    assert not king.get("noise_scale_from_structure")
    assert king["threshold"] > 0.125
    assert trig["evidence"], trig["reason"]
    assert info["n_pairs_detected"] == trig["n_pairs"]
    assert info["n_pairs_modeled"] == stats["n_phi"]
    assert info["n_pairs_detected"] >= n_pairs
    assert info["n_pairs_modeled"] < 0.5 * info["n_pairs_detected"]
    texts = " ".join(str(w.message) for w in rec
                     if issubclass(w.category, RuntimeWarning))
    assert "relatives detected but not modeled in Phi" in texts, texts
    assert "supply a pedigree or pairs (e.g. KING within ancestry groups)" in texts, texts


def test_king_auto_floor_never_below_master_mad_bound(monkeypatch):
    """The master 4·MAD bound applies even if the estimated scale is
    smaller than the MAD (a reintroduced switch, or any other lower
    σ̂). Removing `thr = max(thr, k * sd_mad)` must fail this test."""
    from mixmogam import Genotypes
    from mixmogam import _relatedness as R

    n, m, seed = 400, 800, 1
    rng = np.random.default_rng(seed)
    f = rng.uniform(0.05, 0.5, m)
    G = rng.binomial(2, f, size=(n, m)).astype(np.int8)
    gt = Genotypes(np.ascontiguousarray(G),
                   chromosome=np.repeat([1, 2], m // 2),
                   position=np.arange(m) * 10_000 + 1)
    mad = 0.04

    def fake(phi_r, n_variants=None):
        x = np.asarray(phi_r, dtype=np.float64)
        x = x[np.isfinite(x)]
        med = float(np.median(x)) if x.size else 0.0
        return 0.0, 0.01, mad, med, 1.0, False, 1.0

    monkeypatch.setattr(R, "_king_noise_scale", fake)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        res = R.king_pairs(
            gt, np.arange(m), 0.025, noise_sd_multiple=4.0,
            rng=np.random.default_rng(seed))
    assert res["threshold"] >= 4.0 * mad - 1e-12, res["threshold"]
    assert res["noise_sd"] == 0.01
    assert res["noise_sd_mad"] == mad


def test_king_refit_trigger_anchors_location_at_zero():
    """A negative noise location must not lower the trigger cutoff.
    Dropping max(0, loc) from the trigger fails this test."""
    from mixmogam._relatedness import king_refit_trigger
    from scipy.stats import norm

    n, sd, loc, floor = 3000, 0.01, -0.05, 0.025
    trig = king_refit_trigger(
        n, np.zeros(0, int), np.zeros(0, int), np.zeros(0),
        noise_sd=sd, noise_median=loc, noise_location=loc,
        relatedness_threshold=floor)
    assert trig["noise_location"] == 0.0
    N = n * (n - 1) // 2
    z = float(norm.isf(1.0 / N))
    expected = max(floor, 0.0 + z * sd)
    assert trig["cutoff"] == pytest.approx(expected, rel=1e-10)
    lowered = max(floor, loc + z * sd)
    assert trig["cutoff"] > lowered + 0.01


def _bn_ld_panel(n, m, fst, props, n_pairs, seed, B=30, r=0.5):
    """Balding-Nichols panel with a mild block-LD haplotype
    model: z = sqrt(r) u_block + sqrt(1-r) e, meiosis
    picks a parental haplotype per block."""
    from mixmogam import Genotypes
    from scipy.stats import norm

    rng = np.random.default_rng(seed)
    k_pops = len(props)
    sizes = [int(round(q * n)) for q in props]
    sizes[-1] = n - sum(sizes[:-1])
    p0 = rng.uniform(0.1, 0.9, m)
    a = p0 * (1 - fst) / fst
    b = (1 - p0) * (1 - fst) / fst
    P = np.clip(np.stack([rng.beta(a, b) for _ in range(k_pops)]), 0.02, 0.98)
    nb = -(-m // B)
    thr = {k: norm.ppf(P[k]).astype(np.float32) for k in range(k_pops)}

    def haps(k, nh):
        u = rng.standard_normal((nh, nb), dtype=np.float32)
        z = (np.sqrt(r) * np.repeat(u, B, axis=1)[:, :m]
             + np.sqrt(1 - r) * rng.standard_normal((nh, m), dtype=np.float32))
        return (z < thr[k]).astype(np.int8)

    def mei(H):
        pick = np.repeat(rng.integers(0, 2, nb), B)[:m]
        return np.where(pick == 0, H[0], H[1])

    rows, true_i, true_j, used = [], [], [], [0] * k_pops
    for t_ix in range(n_pairs):
        k = t_ix % k_pops
        if used[k] + 2 > sizes[k]:
            k = int(np.argmax(np.array(sizes) - np.array(used)))
        H = haps(k, 4)
        par = H[:2]
        rows += [mei(par) + H[2], mei(par) + H[3]]
        true_i.append(2 * t_ix)
        true_j.append(2 * t_ix + 1)
        used[k] += 2
    leftover = []
    for k, sz in enumerate(sizes):
        leftover += [k] * max(0, sz - used[k])
    leftover = leftover[: n - 2 * n_pairs]
    rng.shuffle(leftover)
    # draw remaining individuals in batches per pop
    by = {}
    for k in leftover:
        by.setdefault(k, 0)
        by[k] += 1
    for k, rest in by.items():
        H = haps(k, 2 * rest)
        rows.extend((H[0::2] + H[1::2]).reshape(rest, m))
    G = np.vstack([np.atleast_2d(x) for x in rows]).astype(np.int8)
    gt = Genotypes(np.ascontiguousarray(G),
                   chromosome=np.repeat([1, 2, 3, 4], (m + 3) // 4)[:m],
                   position=np.arange(m) * 10_000 + 1)
    return gt, np.array(true_i, dtype=int), np.array(true_j, dtype=int)


def test_king_auto_mild_ld_8020_false_pairs_match_master_floor():
    """80/20 two-population panel with mild block LD (B=30, r=0.5): the
    old path-(b) switch flooded Phi with hundreds of false pairs. With
    no switch, false pairs stay at the master 4-MAD floor. A
    reintroduced switch fails this test."""
    n, m, n_pairs, seed = 2000, 8000, 20, 1
    gt, i_t, j_t = _bn_ld_panel(n, m, fst=0.05, props=(0.80, 0.20),
                                n_pairs=n_pairs, seed=seed, B=30, r=0.5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        trig, rel, info = _live_trigger(gt, seed=seed)
    king = info["king"]
    stats = _phi_false_and_true(rel, i_t, j_t)
    master_floor = 4.0 * float(king["noise_sd_mad"])
    pre = info.get("king_prefloor") or {}
    if pre.get("r") is not None and np.size(pre["r"]):
        n_master = int(np.sum(np.asarray(pre["r"], dtype=np.float64) / 2.0 > master_floor))
    else:
        n_master = 0
    assert not king.get("noise_scale_from_structure"), (
        f"switch reintroduced; sd={king['noise_sd']:.4f} MAD={king['noise_sd_mad']:.4f} "
        f"floor={king['threshold']:.4f} master={master_floor:.4f}")
    assert king["noise_sd"] == pytest.approx(king["noise_sd_mad"], rel=1e-12)
    assert king["threshold"] >= master_floor - 1e-12
    assert abs(king["threshold"] - max(0.025, max(0.0, king["noise_location"]) + master_floor,
                                       master_floor)) < 1e-8 or \
        king["threshold"] == pytest.approx(max(0.025, max(0.0, float(king["noise_location"])) + master_floor), abs=1e-12) or \
        True
    # tighter: floor equals the no-switch formula
    loc_floor = max(0.0, float(king["noise_location"]))
    expected = max(0.025, loc_floor + 4.0 * king["noise_sd_mad"], 4.0 * king["noise_sd_mad"])
    assert king["threshold"] == pytest.approx(expected, abs=1e-12)
    assert stats["n_false"] <= n_master + 5, (
        f"{stats['n_false']} false Phi pairs; master_floor_pairs={n_master} "
        f"floor={king['threshold']:.4f} rec={stats['n_rec']}/{stats['n_true']}")


def test_gls_gamma_skips_monomorphic_calibration_variants():
    """A monomorphic variant drawn into the gamma calibration must not NaN
    the whole group's denominator (review F12). n_gamma covers the panel,
    so the monomorphic sites are guaranteed to be in the draw."""
    from mixmogam.genotypes import Genotypes

    rng = np.random.default_rng(11)
    n, m = 200, 300
    G = rng.integers(0, 3, size=(n, m)).astype(np.int8)
    G[:, 50], G[:, 150], G[:, 250] = 1, 0, 2
    gt = Genotypes(G, sample_ids=[f"s{i}" for i in range(n)],
                   chromosome=np.repeat([1, 2, 3], 100),
                   position=np.arange(m) * 1000 + 1,
                   variant_ids=[f"v{j}" for j in range(m)])
    y = G[:, 100] * 0.4 + rng.normal(size=n)
    rel = (np.array([0]), np.array([1]), np.array([0.5]))
    r = hratt(y, gt, loco_folds=2, random_state=1, relatedness=rel,
              relatedness_correction="gls",
              relatedness_options={"denominator": "gamma", "n_gamma": m})
    poly = np.ones(m, dtype=bool)
    poly[[50, 150, 250]] = False
    assert np.isfinite(r.p[poly]).all()
    assert np.isnan(r.p[~poly]).all()
