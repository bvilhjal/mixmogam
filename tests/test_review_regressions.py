"""Independent file-format, linear-algebra and identity checks from the review."""

import csv

import numpy as np
import pytest
from scipy import linalg, stats

from mixmogam import Genotypes, LMM, gwas
from mixmogam.io.plink import read_plink, write_plink
from mixmogam.kinship import GenotypeKinship, realized_relationship, windowed_kinships
from mixmogam.results import GwasResult


@pytest.fixture
def problem():
    rng = np.random.default_rng(195)
    calls = rng.binomial(2, 0.4, (40, 60))
    gt = Genotypes(calls, chromosome=np.repeat([1, 2], 30))
    y = 4 * rng.normal(size=40) + 5 * calls[:, 1]
    return gt, realized_relationship(gt, dtype=np.float64), y


def test_plink_matches_documented_bits(tmp_path):
    """0xe4 contains 00, 01, 10, 11 in PLINK sample order."""
    p = tmp_path / "known"
    p.with_suffix(".fam").write_text("".join(f"0 s{i} 0 0 0 -9\n" for i in range(4)))
    p.with_suffix(".bim").write_text("X rs1 0 10 T C\n")
    p.with_suffix(".bed").write_bytes(bytes([0x6C, 0x1B, 0x01, 0xE4]))
    gt = read_plink(str(p))
    np.testing.assert_array_equal(gt.G[:, 0], [2, -1, 1, 0])
    assert str(gt.chromosome[0]) == "X"
    assert gt.allele1[0] == "T" and gt.allele2[0] == "C"
    write_plink(gt, str(tmp_path / "written"))
    assert (tmp_path / "written.bed").read_bytes() == p.with_suffix(".bed").read_bytes()
    assert (tmp_path / "written.bim").read_text().split()[-2:] == ["T", "C"]


def test_mac_counts_alleles_and_af_has_a_direction():
    gt = Genotypes([[0, 2, -1], [0, 2, -1], [1, 1, -1], [1, 1, 2]])
    np.testing.assert_array_equal(gt.filter_variants(min_mac=2).variant_ids, ["V0", "V1"])
    np.testing.assert_allclose(gt.allele_freqs(), [0.25, 0.75, 1.0])
    np.testing.assert_allclose(gt.allele_freqs(minor=True), [0.25, 0.25, 0.0])


@pytest.mark.parametrize("calls", [[[0.7]], [[128]], [[np.nan]], [[np.inf]], [[-2]]])
def test_invalid_hard_calls_are_not_silently_cast(calls):
    with pytest.raises(ValueError, match="hard calls"):
        Genotypes(calls)


def test_alignment_preserves_missing_identity():
    gt = Genotypes([[2, 0], [1, 2]], sample_ids=["a", "b"])
    aligned, found = gt.align_samples(["missing", "b"], strict=False)
    np.testing.assert_array_equal(found, [False, True])
    np.testing.assert_array_equal(aligned.G, [[-1, -1], [1, 2]])
    with pytest.raises(ValueError, match="unique"):
        Genotypes([[0], [1]], sample_ids=["a", "a"])
    with pytest.raises(ValueError, match="position"):
        Genotypes([[0, 1]], position=[10])


def test_blup_is_conditional_gaussian_mean_and_scales_linearly(problem):
    gt, K, y = problem
    fit = LMM(y, K=K).fit()
    V = fit.vg * K + fit.ve * np.eye(y.size)
    expected = fit.vg * K @ np.linalg.solve(V, y - fit.model.X @ fit.beta)
    np.testing.assert_allclose(fit.blup(), expected, rtol=1e-8, atol=1e-9)
    np.testing.assert_allclose(LMM(3 * y, K=K).fit().blup(), 3 * expected, rtol=1e-6)
    np.testing.assert_allclose(fit.predict(), fit.model.X @ fit.beta + expected)


def test_ml_variance_and_fit_method_cache(problem):
    _, K, y = problem
    model = LMM(y, K=K)
    model.fit()
    fit = model.fit(method="ml")
    assert fit.method == "ml"
    resid = y - model.X @ fit.beta
    expected = resid @ np.linalg.solve(K + fit.delta * np.eye(y.size), resid) / y.size
    assert fit.vg == pytest.approx(expected, rel=1e-9)


def test_superseded_fit_cannot_silently_use_new_parameters(problem):
    from mixmogam.scan import permutation_min_p

    gt, K, y = problem
    model = LMM(y, K=K)
    old = model.fit(method="ml")
    current = model.fit(method="reml")
    for operation in (old.blup, old.predict, lambda: old.scan(gt),
                      lambda: permutation_min_p(old, gt, n_perm=2)):
        with pytest.raises(ValueError, match="superseded"):
            operation()
    assert np.isfinite(current.blup()).all()


@pytest.mark.parametrize("mixed", [False, True])
@pytest.mark.parametrize("constant", [False, True])
def test_fully_explained_phenotype_is_rejected(problem, mixed, constant):
    gt, K, _ = problem
    x = np.arange(gt.n_samples, dtype=float)
    y = np.ones_like(x) if constant else 7 + 3 * x
    model = LMM(y, X=x, K=K if mixed else None)
    with pytest.raises(ValueError, match="no residual variation"):
        model.fit()
    if not mixed:
        with pytest.raises(ValueError, match="no residual variation"):
            model.scan(gt)


def test_exact_fit_is_invariant_to_large_intercept(problem):
    _, K, y = problem
    plain = LMM(y, K=K).fit()
    shifted = LMM(y + 1e9, K=K).fit()
    assert np.isfinite(shifted.ll)
    assert shifted.delta == pytest.approx(plain.delta, rel=1e-5)
    assert shifted.vg == pytest.approx(plain.vg, rel=1e-6)


@pytest.mark.parametrize("method", ["ml", "reml"])
def test_no_kinship_fit_is_ordinary_linear_model(problem, method):
    _, _, y = problem
    fit = LMM(y).fit(method=method)
    df = y.size if method == "ml" else y.size - 1
    assert fit.vg == 0 and fit.pseudo_heritability == 0
    assert fit.ve == pytest.approx(np.sum((y - y.mean()) ** 2) / df)
    np.testing.assert_allclose(fit.predict(), np.full(y.size, y.mean()))


def test_weighted_operator_matches_direct_weighted_grm(problem):
    gt, _, y = problem
    weights = np.geomspace(0.01, 10, gt.n_variants)
    expected = realized_relationship(gt, weights=weights, dtype=np.float64)
    for cache_bytes in (0, 10**8):
        op = GenotypeKinship(gt, weights=weights, dtype=np.float64, cache_bytes=cache_bytes)
        np.testing.assert_allclose(op @ y, expected @ y, atol=1e-10)
        np.testing.assert_allclose(op.diagonal(), np.diag(expected), atol=1e-12)


@pytest.mark.parametrize("jump", [10, 25])
def test_window_complements_with_overlap_or_gaps(problem, jump):
    gt, _, _ = problem
    for wi, local, rest in windowed_kinships(gt, 20, jump, scale=False, dtype=np.float64):
        start = wi * jump
        on = (np.arange(gt.n_variants) >= start) & (np.arange(gt.n_variants) < start + 20)
        expected = realized_relationship(gt.variant_mask(np.flatnonzero(~on)), scale=False, dtype=np.float64)
        np.testing.assert_allclose(rest, expected, atol=1e-12)


def test_posterior_matches_integrated_normal_and_is_subset_invariant():
    beta, se = np.array([0.05, 0.5]), np.array([0.1, 0.1])
    prior, W = np.array([0.1, 0.1]), 0.04
    res = GwasResult(np.ones(2), np.arange(2), np.ones(2), beta=beta, se=se)
    res.posterior_probabilities(prior, prior_variance=W)
    bf = stats.norm.pdf(beta, scale=np.sqrt(se**2 + W)) / stats.norm.pdf(beta, scale=se)
    expected = bf * prior / (1 - prior + bf * prior)
    np.testing.assert_allclose(res.extra["ppa"], expected)
    single = res.take([1]).posterior_probabilities(prior[1:], prior_variance=W)
    assert single.extra["ppa"][0] == pytest.approx(expected[1])
    with pytest.raises(ValueError, match="prior_variance"):
        res.posterior_probabilities(prior)


def test_csv_preserves_identifiers_effects_and_precision(tmp_path):
    path = tmp_path / "result.csv"
    res = GwasResult(np.array(["X"]), np.array([10]), np.array([1.234567890123456e-90]),
                     variant_ids=np.array(['rs,"quoted"']), beta=np.array([-0.3212345678901]),
                     se=np.array([0.012345678901]), rss=np.array([2.5]),
                     effect_allele=np.array(["T"]), other_allele=np.array(["C"]))
    res.write_csv(path)
    assert list(csv.DictReader(path.open()))[0]["variant_id"] == 'rs,"quoted"'
    back = GwasResult.read_csv(path)
    assert len(back) == 1
    for name in ("p", "variant_ids", "beta", "se", "rss", "effect_allele", "other_allele"):
        np.testing.assert_array_equal(getattr(back, name), getattr(res, name))


def test_gwas_rejects_impossible_loco_and_unused_options(problem):
    gt, _, y = problem
    with pytest.raises(ValueError, match="two"):
        gwas(y, Genotypes(gt.G))
    with pytest.raises(TypeError, match="unexpected"):
        gwas(y, gt, denominator="spectral")
    with pytest.raises(ValueError, match="sample count"):
        gwas(y[:-1], gt)


def test_invalid_covariance_and_design_rejected(problem):
    _, K, y = problem
    with pytest.raises(ValueError, match="rank"):
        LMM(y, K=K, X=np.ones(y.size))
    with pytest.raises(ValueError, match="finite"):
        LMM(y, K=K * np.nan)
    with pytest.raises(ValueError, match="positive semidefinite"):
        LMM(y, K=-np.eye(y.size)).fit()


def test_collinear_snp_has_no_estimable_effect(problem):
    gt, K, y = problem
    model = LMM(y, K=K, X=gt.G[:, 0])
    res = model.fit().scan(gt, dtype=np.float64, with_betas=True)
    assert np.isnan(res["ps"][0]) and np.isnan(res["betas"][0])


def test_residual_permutations_preserve_null_geometry(problem):
    from mixmogam.scan import _permuted_residuals

    _, K, y = problem
    X = np.column_stack([np.arange(y.size), np.arange(y.size) == 0])
    model = LMM(y, K=K, X=X)
    model.fit()
    fac = model._scan_factors(np.float64)
    Rp = _permuted_residuals(fac["Q"], fac["r"], 7, np.random.default_rng(11))
    np.testing.assert_allclose(fac["Q"].T @ Rp, 0, atol=1e-11)
    np.testing.assert_allclose(np.sum(Rp**2, axis=0), fac["rss0"], rtol=1e-12)
    # Independent explicit complement: small test only, never in production.
    U = linalg.qr(fac["Q"], mode="full")[0][:, model.q:]
    xi = U.T @ fac["r"]
    perms = np.argsort(np.random.default_rng(11).random((xi.size, 7)), axis=0)
    np.testing.assert_allclose(Rp, U @ xi[perms], atol=1e-11)


def test_standard_tped_allele_pairs(tmp_path):
    from mixmogam.io.plink import read_tped

    p = tmp_path / "pairs"
    p.with_suffix(".tfam").write_text("".join(f"0 s{i} 0 0 0 -9\n" for i in range(4)))
    p.with_suffix(".tped").write_text("Y rs1 0 10 A A A G G G 0 0\n")
    gt = read_tped(str(p))
    np.testing.assert_array_equal(gt.G[:, 0], [2, 1, 0, -1])
    assert gt.allele1[0] == "A" and gt.allele2[0] == "G"


def test_regmap_merges_variants_and_aligns_file_headers(tmp_path):
    from mixmogam.io.regmap import read_regmap

    a, b = tmp_path / "one.csv", tmp_path / "two.csv"
    a.write_text("Chromosome,Position,a,b\n1,1,0,2\n1,2,1,2\n")
    b.write_text("Chromosome,Position,b,a\n2,1,0,1\n")
    gt = read_regmap([a, b])
    np.testing.assert_array_equal(gt.G, [[0, 1, 1], [2, 2, 0]])
    np.testing.assert_array_equal(gt.chromosome, [1, 1, 2])


def test_hdf5_keeps_alleles_and_string_chromosomes(tmp_path):
    pytest.importorskip("h5py")
    from mixmogam.io.hdf5 import read_hdf5, write_hdf5

    gt = Genotypes([[0, 2], [1, -1]], chromosome=["X", "Y"],
                   allele1=["A", "C"], allele2=["G", "T"])
    path = tmp_path / "alleles.h5"
    write_hdf5(gt, path)
    back = read_hdf5(path)
    for name in ("G", "chromosome", "allele1", "allele2"):
        np.testing.assert_array_equal(getattr(back, name), getattr(gt, name))


def test_replicates_cannot_overwrite_or_corrupt_shared_axis(tmp_path):
    from mixmogam.io.phenofile import read_phenotypes
    from mixmogam.phenotypes import Phenotypes

    path = tmp_path / "repeat.csv"
    path.write_text("phenotype_id,sample_id,value,replicate_id\nt,a,1,1\nt,a,2,2\n")
    with pytest.raises(ValueError, match="duplicate observation"):
        read_phenotypes(path)
    ph = Phenotypes(["a", "a", "b", "b"])
    ph.replicates = np.array([1, 1, 2, 2])
    ph.add("x", [1, 2, 3, 4])
    ph.add("y", [4, 3, 2, 1])
    with pytest.raises(ValueError, match="single untransformed"):
        ph.convert_to_averages("x")
    assert len(ph.sample_ids) == len(ph.values("x")) == len(ph.values("y")) == 4


def test_small_spectral_scan_and_forwarded_block(problem, monkeypatch):
    from mixmogam import twostep

    gt, _, y = problem
    res = twostep.bolt_inf(y, gt, denominator="spectral")
    assert np.isfinite(res.p).all()
    assert 0 < res.extra["spectral_k"] <= y.size - 2
    called = {}

    def fake(*args, **kwargs):
        called.update(kwargs)
        return res

    monkeypatch.setattr(twostep, "bolt_inf", fake)
    assert gwas(y, gt, method="bolt-inf", block=7) is res
    assert called["block"] == 7


def test_unconverged_loco_is_rejected(problem, monkeypatch):
    from mixmogam import twostep

    gt, _, y = problem
    st = twostep._setup(y, gt, None, 25, 16)
    monkeypatch.setattr(twostep, "batched_pcg", lambda *a, **k: (np.zeros((40, 2)), {"converged": False}))
    with pytest.raises(RuntimeError, match="did not converge"):
        twostep._loco_solve(st, 1.0, np.tile(y[:, None], (1, 2)), np.arange(2), pre=lambda x: x)


def test_operator_prediction_includes_genetic_value(problem):
    gt, _, y = problem
    op = GenotypeKinship(gt, dtype=np.float64)
    fit = LMM(y, K=op, random_state=5).fit(slq_deflate=0)
    K = op @ np.eye(y.size)
    u = K @ np.linalg.solve(K + fit.delta * np.eye(y.size), y - fit.model.X @ fit.beta)
    np.testing.assert_allclose(fit.blup(), u, rtol=1e-6, atol=1e-8)
    np.testing.assert_allclose(fit.predict(), fit.model.X @ fit.beta + u, rtol=1e-6)


def test_linear_stepwise_does_not_stop_for_zero_genetic_variance(problem):
    from mixmogam.stepwise import mlmm

    gt, _, y = problem
    out = mlmm(y, gt, max_steps=4, backward=False)
    assert len([s for s in out["steps"] if s["action"] == "+"]) == 4
