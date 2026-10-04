"""Check a wheel installed outside the checkout; no benchmark timings."""

import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys

import numpy as np
import mixmogam
from mixmogam import Genotypes, gwas


def main():
    source = Path(sys.argv[1]).resolve()
    installed = Path(mixmogam.__file__).resolve().parent
    assert not installed.is_relative_to(source)
    assert "site-packages" in str(installed)
    assert mixmogam.__version__ == importlib.metadata.version("mixmogam") == "2.0.0.dev5"
    hashes = {}
    for path in sorted((source / "mixmogam").rglob("*.py")):
        relative = path.relative_to(source / "mixmogam")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert hashlib.sha256((installed / relative).read_bytes()).hexdigest() == digest
        hashes[str(relative)] = digest
    rng = np.random.default_rng(10473)
    G = rng.binomial(2, rng.uniform(0.1, 0.5, 120), size=(144, 120)).astype(np.int8)
    gt = Genotypes(G, chromosome=np.repeat(np.arange(6), 20))
    y = (G - G.mean(0)) @ rng.normal(0, 0.09, 120) + rng.normal(size=144)
    fits = {}
    for method in ("exact", "bolt-inf", "kvik"):
        options = {} if method == "exact" else {"random_state": 37}
        if method == "kvik":
            options.update(vb_max_iter=300)
        fits[method] = gwas(y, gt, method=method, **options)
    options = dict(method="kvik", heritability_method="he", n_threads=4,
                   random_state=37, vb_max_iter=300)
    fits["kvik-he-cached"] = gwas(y, gt, cache_bytes=4e9, **options)
    fits["kvik-he-uncached"] = gwas(y, gt, cache_bytes=0, **options)
    fits["exact-4-threads"] = gwas(y, gt, method="exact", n_threads=4)
    for result in fits.values():
        assert np.isfinite(result.p).all() and ((result.p >= 0) & (result.p <= 1)).all()
        assert np.isfinite(result.beta).all() and np.isfinite(result.se).all()
        if "cv_converged" in result.extra:
            assert result.extra["cv_converged"] and result.extra["loco_converged"]
    for name in ("p", "beta", "se"):
        np.testing.assert_array_equal(getattr(fits["kvik-he-cached"], name),
                                      getattr(fits["kvik-he-uncached"], name))
        np.testing.assert_array_equal(getattr(fits["exact"], name),
                                      getattr(fits["exact-4-threads"], name))
    print(json.dumps({"version": mixmogam.__version__, "installed_path": str(installed),
                      "package_source_sha256": hashes,
                      "valid_variant_pvalues_by_method": {key: int(np.isfinite(fit.p).sum())
                                                           for key, fit in fits.items()},
                      "he_cache_arrays_equal": True,
                      "exact_thread_arrays_equal": True,
                      "all_variational_fits_converged": True}, indent=2))


if __name__ == "__main__":
    main()
