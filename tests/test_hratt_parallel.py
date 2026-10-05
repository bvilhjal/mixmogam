"""Public HRATT threading preserves the complete association calculation."""

import numpy as np
import pytest

from mixmogam import gwas, twostep
from mixmogam.genotypes import Genotypes


@pytest.mark.parametrize("threads", [0, -1, 1.5, True, None, "2"])
def test_invalid_thread_count_fails_before_genotype_preparation(threads, monkeypatch):
    def unexpected(*args, **kwargs):
        raise AssertionError("invalid thread count reached genotype preparation")

    monkeypatch.setattr(twostep, "_setup", unexpected)
    with pytest.raises(ValueError, match="n_threads"):
        twostep.hratt(None, None, n_threads=threads)


@pytest.mark.numba
@pytest.mark.parametrize("threads", [2, 4])
@pytest.mark.parametrize("fst", [0.0, 0.08])
def test_parallel_gwas_preserves_inference_on_phensim_data(threads, fst):
    numba = pytest.importorskip("numba")
    phensim = pytest.importorskip("phensim")
    if threads > numba.config.NUMBA_NUM_THREADS:
        pytest.skip("run with NUMBA_NUM_THREADS >= 4 to exercise all thread counts")
    G, _ = phensim.simulate_population_structure(
        192, 180, n_pops=2, fst=fst, model="balding-nichols", seed=962,
        block_sizes=[30] * 6, rho=0.4)
    Z = (G - G.mean(axis=0)) / G.std(axis=0)
    K = Z @ Z.T / G.shape[1]
    _, U = np.linalg.eigh(K)
    X = U[:, -2:] if fst else None
    phenotype = phensim.simulate_confounded_trait(
        G, h2=0.3, confounding_strength=0.2 if fst else 0.0,
        environment=U[:, -1], architecture="infinitesimal", n_causal=0,
        K=K, seed=965)
    gt = Genotypes(G, chromosome=np.repeat(np.arange(6), 30))
    options = dict(method="hratt", X=X, heritability_method="he",
                   vb_max_iter=300, random_state=967)
    before = numba.get_num_threads()
    serial = gwas(phenotype["liability"], gt, n_threads=1, **options)
    parallel = gwas(phenotype["liability"], gt, n_threads=threads, **options)
    assert numba.get_num_threads() == before
    for name, value in vars(serial).items():
        if isinstance(value, np.ndarray):
            if name in ("beta", "se"):
                # The per-variant SD reduction changes the final effect-unit
                # conversion by a few ulps; standardized inference is exact.
                np.testing.assert_allclose(getattr(parallel, name), value,
                                           rtol=32 * np.finfo(float).eps, atol=0)
            else:
                np.testing.assert_array_equal(getattr(parallel, name), value)
    for key in ("h2", "alpha", "cv_mse", "cv_best", "lambda", "loco_iterations",
                "cv_converged", "loco_converged"):
        np.testing.assert_array_equal(parallel.extra[key], serial.extra[key])
    assert serial.extra["cv_converged"] and serial.extra["loco_converged"]
