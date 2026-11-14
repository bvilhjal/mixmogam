"""Spec A (SG specs-for-wizard.md): relatedness= kappa_g correction of HRATT."""
import warnings

import numpy as np
import pytest

from mixmogam._relatedness import (KING_SECOND_DEGREE, RELATEDNESS_OPTIONS, Relatedness,
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
    c = _run(y, gt, relatedness=few)
    assert c.extra["loco_refit"] and c.extra["loco_refit_reason"] == "relatives"
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
    floor = k.extra["relatedness"]["king"]["threshold"]
    assert floor > 0.0884
    msg = [str(x.message) for x in w if "no covariates" in str(x.message)]
    assert len(msg) == 1 and "noise floor" in msg[0] and "degree relatives" in msg[0]
    assert "kappa_g replaces" in msg[0] and "LDpred3" in msg[0] and "in simulations" in msg[0]
    assert any("2nd-degree bound" in str(x.message) for x in w)
    pcs = np.linalg.svd(G - G.mean(0), full_matrices=False)[0][:, :2]
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        hratt(y, gt, X=pcs, relatedness=empty, random_state=1,
              relatedness_options={"structure_without_pcs": "error"})
    assert not any("no covariates" in str(x.message) for x in w)


# (i) loco_refit='auto' needs evidence of relatives: supplied pairs count from
# one pair, KING-auto pairs need >= 1% of samples related.
def test_refit_evidence_king_auto_needs_min_related_fraction():
    G, gt, y, pairs, causal = family_panel(seed=7, kids=0, n_fam=300, m=4000, n_causal=100)
    from mixmogam import Genotypes

    def twins(k):  # duplicate k samples: k pairs with phi 0.5 that KING finds
        rows = np.r_[np.arange(G.shape[0]), np.arange(k)]
        g2 = Genotypes(G[rows], chromosome=gt.chromosome, position=gt.position)
        return g2, np.r_[y, y[:k] + 0.1], (np.arange(k), G.shape[0] + np.arange(k), np.full(k, 1.0))

    g1, y1, p1 = twins(1)  # 2 of 601 related: 0.3%
    a = _run(y1, g1, relatedness="auto")
    ra = a.extra["relatedness"]
    assert ra["n_pairs_above_threshold"] >= 1 and ra["related_fraction"] < 0.01
    assert not ra["relatives_evidence"] and not a.extra["loco_refit"]
    b = _run(y1, g1, relatedness=p1)  # the same pair supplied: counts
    assert b.extra["loco_refit"] and b.extra["loco_refit_reason"] == "relatives"
    g10, y10, _ = twins(10)  # 20 of 610: 3.3%
    c = _run(y10, g10, relatedness="auto")
    assert c.extra["relatedness"]["related_fraction"] >= 0.01
    assert c.extra["loco_refit"] and c.extra["loco_refit_reason"] == "relatives"
    trig = c.extra["relatedness"]["refit_trigger"]
    assert trig["evidence"] and trig["n_pairs"] >= 10 and trig["cutoff"] >= 0.025
    assert "n_expected_null" in trig and trig["n_expected_null"] < 2
    G2, gt2, y2, _, _ = family_panel(seed=8, n_fam=100, kids=2, m=4000, n_causal=100)
    d = _run(y2, gt2, relatedness="auto")  # nuclear families: KING-auto still refits
    assert d.extra["relatedness"]["related_fraction"] > 0.9 and d.extra["loco_refit"]


# (v) family folds: opt-in, no refit, in-model nulls calibrated
def test_family_folds_opt_in_calibrated_without_refit():
    G, gt, y, pairs, causal = family_panel(seed=0, n_fam=200, kids=4, m=8000, c2=0.3, h2=0.5)
    r = _run(y, gt, relatedness=pairs, relatedness_options={"crossfit_folds": "family"})
    assert r.extra["crossfit_folds"] == "family" and not r.extra["loco_refit"]
    ff = r.extra["relatedness"]["family_folds"]
    assert ff["largest_family"] == 6 and max(ff["fold_sizes"]) - min(ff["fold_sizes"]) <= 6
    assert _null_mean(r, causal) < 1.02
    assert _run(y, gt, relatedness=pairs).extra["crossfit_folds"] == "random"
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
    kap = _run(y, gt, model_snps=sub, relatedness=pairs, sample_weights=w,
               stage1="auto", loco_folds=2, spa_threshold=1e6)
    sel = kap.extra["stage1_selection"]
    assert sel["chosen"] == "enet" and sel["reason"].startswith("sample_weights")
    assert kap.extra.get("stage1", "enet") == "enet"
    assert not kap.extra["structure"]["strong"]
    assert kap.extra["loco_refit"] and kap.extra["loco_refit_reason"] == "relatives"
    base = _run(y, gt, model_snps=sub, relatedness=pairs, relatedness_correction="none",
                sample_weights=w, stage1="auto", loco_folds=2, spa_threshold=1e6)
    k, cx = kap.extra["relatedness"]["kappa"], kap.extra["relatedness"]["c_X"]
    np.testing.assert_allclose(kap.f_stat, base.f_stat / (k[gt.chromosome - 1] * cx),
                               rtol=1e-10)


def test_king_refit_trigger_ignores_n2_noise_but_fires_on_true_relatives():
    """SG p20k §4: 244 4-SD-floor pairs on an unrelated n=20,000 panel
    (floor 0.042, 1.33% of samples) must not fire. First-degree and 10%
    first-cousin pairs at the same noise scale must still fire. Phi is
    unchanged: the trigger only filters which stored pairs count."""
    n, sd, med, floor = 20_000, 0.0418 / 4.0, 0.0, 0.0418
    rng = np.random.default_rng(1)
    # 266 samples in 244 almost-disjoint noise pairs (1.33% of n), phi just
    # above the 4-SD floor, as in the SG table.
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
        rel, info = resolve_relatedness(
            "auto", gt, mean=np.asarray(gt.G, dtype=np.float64).mean(0),
            model_idx=None, threshold=0.025, options=RELATEDNESS_OPTIONS,
            rng=np.random.default_rng(seed))
        king = info["king"]
        u = rel.upper.tocoo()
        return king_refit_trigger(gt.n_samples, u.row, u.col, u.data / 2.0,
                                  noise_sd=king["noise_sd"],
                                  noise_median=king["noise_median"],
                                  noise_tail_ratio=king.get("noise_tail_ratio", 1.0),
                                  n_variants=king.get("n_variants"))

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
    rng = np.random.default_rng(0)
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
        rel, info = resolve_relatedness(
            "auto", gt, mean=np.asarray(gt.G, dtype=np.float64).mean(0),
            model_idx=None, threshold=0.025, options=RELATEDNESS_OPTIONS,
            rng=np.random.default_rng(seed))
        king = info["king"]
        u = rel.upper.tocoo()
        return king_refit_trigger(gt.n_samples, u.row, u.col, u.data / 2.0,
                                  noise_sd=king["noise_sd"],
                                  noise_median=king["noise_median"],
                                  noise_tail_ratio=king.get("noise_tail_ratio", 1.0),
                                  n_variants=king.get("n_variants"))

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
        rel, info = resolve_relatedness(
            "auto", gt, mean=G.mean(0).astype(np.float64), model_idx=None,
            threshold=0.025, options=RELATEDNESS_OPTIONS,
            rng=np.random.default_rng(seed))
        king = info["king"]
        u = rel.upper.tocoo()
        trig = king_refit_trigger(n, u.row, u.col, u.data / 2.0,
                                  noise_sd=king["noise_sd"],
                                  noise_median=king["noise_median"],
                                  noise_tail_ratio=king.get("noise_tail_ratio", 1.0),
                                  n_variants=king.get("n_variants"))
        return trig, int(u.row.size)

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

