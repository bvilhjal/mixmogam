"""What the Gram cache buys: BOLT CV and KVIK-like fits with and without it.
python gram_value.py SOURCE PANEL_DIR"""
import json
import sys
import time

sys.path.insert(0, sys.argv[1])
import numpy as np
from mixmogam import Genotypes, twostep
from mixmogam._vb import VBEngine, PRIOR_MIXTURE, PRIOR_ENET

data = sys.argv[2]
spec = json.load(open(f"{data}/panel.json"))["spec"]
G = np.load(f"{data}/G.npy")
y, X = np.load(f"{data}/y.npy"), np.load(f"{data}/X.npy")
per = G.shape[1] // spec["chromosomes"]
gt = Genotypes(G, chromosome=np.repeat(np.arange(1, spec["chromosomes"] + 1), per))
del G
st = twostep._setup(y, gt, X, 25, 4096, 0) if "cache_bytes" in twostep._setup.__code__.co_varnames else twostep._setup(y, gt, X, 25, 4096)
n, m = st.lg.n, st.lg.m
rng = np.random.default_rng(0)
s2 = 0.4 / m
# BOLT CV: 18 priors x 5 folds
folds = rng.permutation(n) % 5
grid = twostep.BOLT_GRID
P = len(grid) * 5
col_fold = np.tile(np.arange(5), len(grid))
prior = np.array([(p, (1 - f2) * s2 / p, f2 * s2 / (1 - p)) for f2, p in grid for _ in range(5)])
for budget in (1e9, 0):
    eng = VBEngine(st.lg, folds=folds, gram_cache_bytes=budget)
    t = time.perf_counter()
    out = eng.fit(np.repeat(st.y_p[:, None], P, axis=1), col_fold, np.full(P, -1), PRIOR_MIXTURE,
                  prior, np.full(P, 0.6), max_iter=100, tol=1e-5)
    print(f"bolt CV budget {budget:.0e}: {time.perf_counter() - t:.1f} s, sweeps {out['iterations']}, cached blocks "
          f"{0 if eng._grams is None else len(eng._grams)}", flush=True)
# KVIK CV: 10 priors, one 10% fold
held = rng.choice(n, n // 10, replace=False)
folds = np.ones(n, dtype=np.int64)
folds[held] = 0
Pk = len(twostep.KVIK_GRID)
prior = np.array([(p, np.sqrt(2 * p / (1 - F)) if p > 0 else 0.0, F / (1 - p) if p < 1 else 0.0)
                  for p, F in twostep.KVIK_GRID])
for budget in (1e9, 0):
    eng = VBEngine(st.lg, folds=folds, gram_cache_bytes=budget)
    t = time.perf_counter()
    out = eng.fit(np.repeat(st.y_p[:, None], Pk, axis=1), np.zeros(Pk, dtype=np.int64), np.full(Pk, -1),
                  PRIOR_ENET, prior, np.full(Pk, 0.6), snp_scale=np.full(m, 0.4 / m), max_iter=100, tol=1e-5)
    print(f"kvik CV budget {budget:.0e}: {time.perf_counter() - t:.1f} s, sweeps {out['iterations']}", flush=True)
