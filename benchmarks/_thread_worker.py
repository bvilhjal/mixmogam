"""One (fit, kinship, scan) measurement under externally-set BLAS threads."""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from mixmogam import LMM
from mixmogam.kinship import realized_relationship
from sim_data import make_dataset

ds = make_dataset(n=2000, m=20000, h2=0.5, n_causal=20, seed=301)
gt, y = ds["gt"], ds["y"]

t0 = time.perf_counter()
K = realized_relationship(gt)
t_kin = time.perf_counter() - t0
t0 = time.perf_counter()
fit = LMM(y, K=K).fit()
t_fit = time.perf_counter() - t0
t0 = time.perf_counter()
fit.scan(gt, dtype=np.float32)
t_scan = time.perf_counter() - t0
print(f"{t_kin:.4f},{t_fit:.4f},{t_scan:.4f}")
