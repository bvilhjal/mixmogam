"""Sparse close-relative relatedness for HRATT's stage 2 (``relatedness=``).

A :class:`Relatedness` holds the off-diagonal relatedness r_ij = 2 phi_ij
(phi: kinship coefficient) of the sample pairs that are related, e.g. from a
pedigree (:func:`pedigree_relatedness`) or from KING-robust kinship
(:func:`king_relatedness`; Manichaikul et al. 2010). It must not be estimated
from the stage-1 model SNPs' fit; see :func:`hratt`.

Two corrections use it (``relatedness_correction``):

``"kappa"``: for a null variant unlinked to the trait loci, with genotype
covariance sigma^2 R (R = I + offdiag(r)) given the pedigree, and the
projected LOCO residual w_g of its group,
E[(z'w_g)^2 | w_g] = sigma^2 w_g' R w_g, so the HRATT statistic
U = n_eff (z'w_g)^2 / (z'z |w_g|^2) has expectation
kappa_g = 1 + w_g' offdiag(r) w_g / |w_g|^2. Dividing U by kappa_g is the
retrospective variance of MASTOR (Jakobsdottir & McPeek 2013) reduced to one
factor per LOCO group; it replaces lambda.

``"vc"``: a second variance component. Cov(w_g) = s_e I + s_f R with s_f, s_e
moment-estimated from the pairs (Haseman-Elston on w_g), and the GLS score
test chi2 = (z'P w_g)^2 / (z'P z), P = V^-1 - V^-1 Q (Q'V^-1 Q)^-1 Q'V^-1
(fastGWA-like; Jiang et al. 2019). ``denominator="exact"`` computes z'P z per
variant with the block-diagonal V^-1 (one dense inverse per related
cluster); ``"gamma"`` uses gamma_g z'(I - QQ')z with gamma_g the mean ratio
over random calibration variants (GRAMMAR-Gamma style; Svishcheva et al.
2012), which costs nothing per variant.
"""

from __future__ import annotations


import numpy as np
from scipy import linalg, sparse
from scipy.sparse.csgraph import connected_components

__all__ = ["Relatedness", "pedigree_relatedness", "king_relatedness", "king_pairs",
           "king_refit_trigger",
           "as_relatedness", "read_kin0", "read_grm_sparse", "align_pairs", "unrelated_set",
           "covariate_factor", "resolve_relatedness"]

import os
import warnings

# KING kinship boundaries (Manichaikul et al. 2010): 3rd degree from 2^-4.5.
KING_THIRD_DEGREE = 0.0442
KING_SECOND_DEGREE = 0.0884  # KING's 2nd-degree lower bound (phi 1/8 nominal)
KING_FIRST_DEGREE = 0.177  # KING's 1st-degree lower bound (phi 1/4 nominal)


def dropped_degree(floor: float) -> str:
    """Closest relatives a KING floor of ``floor`` can drop: '' when it is
    below the 2nd-degree bound."""
    if floor > KING_FIRST_DEGREE:
        return "1st-degree relatives (phi 0.25: parents, offspring, full sibs)"
    if floor > KING_SECOND_DEGREE:
        return "2nd-degree relatives (phi 0.125)"
    return ""


class Relatedness:
    """Off-diagonal relatedness r_ij = 2 phi_ij of related sample pairs.

    ``i``, ``j``: sample indices (0-based, gt's sample order), ``r``: values.
    Pairs may be given once in either order; duplicates are summed."""

    def __init__(self, n: int, i, j, r, *, source: str = ""):
        i = np.asarray(i, dtype=np.int64).ravel()
        j = np.asarray(j, dtype=np.int64).ravel()
        r = np.asarray(r, dtype=np.float64).ravel()
        if not (i.size == j.size == r.size):
            raise ValueError("relatedness pairs: i, j and r differ in length")
        if i.size and (min(i.min(), j.min()) < 0 or max(i.max(), j.max()) >= n):
            raise ValueError("relatedness pairs index samples outside 0..n-1")
        if not np.all(np.isfinite(r)):
            raise ValueError("relatedness values must be finite")
        keep = (i != j) & (r != 0)
        a, b, v = np.minimum(i, j)[keep], np.maximum(i, j)[keep], r[keep]
        up = sparse.coo_matrix((v, (a, b)), shape=(n, n)).tocsr()
        up.sum_duplicates()
        up.eliminate_zeros()
        self.n = int(n)
        self.upper = up
        self.Koff = (up + up.T).tocsr()
        self.source = source
        self.n_clusters, self.labels = connected_components(self.Koff, directed=False)
        sizes = np.bincount(self.labels)
        self.cluster_sizes = sizes

    @property
    def n_pairs(self) -> int:
        return int(self.upper.nnz)

    def summary(self) -> dict:
        v = self.upper.data
        s = self.cluster_sizes
        return {"source": self.source, "n_pairs": self.n_pairs,
                "mean_sq_per_sample": float(2 * np.sum(v * v) / self.n),
                "related_samples": int(np.count_nonzero(s[self.labels] > 1)),
                "largest_cluster": int(s.max()) if s.size else 0,
                "sum_cluster_sq": int(np.sum(s.astype(np.int64) ** 2))}

    def threshold(self, min_r: float) -> "Relatedness":
        u = self.upper.tocoo()
        k = u.data >= min_r
        return Relatedness(self.n, u.row[k], u.col[k], u.data[k], source=f"{self.source}>={min_r:g}")

    def threshold_phi(self, phi: float) -> "Relatedness":
        """Pairs with phi = r / 2 > ``phi``."""
        u = self.upper.tocoo()
        k = u.data / 2.0 > phi
        if k.all():
            return self
        return Relatedness(self.n, u.row[k], u.col[k], u.data[k], source=self.source)

    def quad(self, W: np.ndarray) -> np.ndarray:
        """w' offdiag(r) w for each column of W (n, k)."""
        W = np.asarray(W, dtype=np.float64)
        if W.ndim == 1:
            W = W[:, None]
        return np.einsum("ij,ij->j", W, self.Koff @ W)

    def kappa(self, Wp: np.ndarray) -> np.ndarray:
        """kappa_g = 1 + w' offdiag(r) w / |w|^2 per column of the projected residuals."""
        ww = np.einsum("ij,ij->j", Wp, Wp)
        with np.errstate(divide="ignore", invalid="ignore"):
            k = 1.0 + self.quad(Wp) / ww
        return np.where(ww > 0, k, 1.0)

    def family_vc(self, Wp: np.ndarray) -> tuple:
        """Moment fit of Cov(w) = s_e I + s_f R per column: s_f from the pairs
        (sum r_ij w_i w_j / sum r_ij^2), s_e = mean w^2 - s_f. s_f is kept in
        [-0.5, 0.95] x mean w^2 (negative: the fit over-explains relatives)."""
        den = float(np.sum(self.upper.data ** 2))
        tot = np.einsum("ij,ij->j", Wp, Wp) / Wp.shape[0]
        sf = self.quad(Wp) / 2.0 / den if den > 0 else np.zeros(Wp.shape[1])
        sf = np.clip(sf, -0.5 * tot, 0.95 * tot)
        return sf, tot - sf

    def vinv(self, sf: float, se: float, max_cluster: int = 5000) -> sparse.csr_matrix:
        """Block-diagonal (s_e I + s_f R)^-1, one dense inverse per cluster."""
        if self.cluster_sizes.max(initial=1) > max_cluster:
            raise ValueError(f"a related cluster has {int(self.cluster_sizes.max()):,} samples "
                             f"(> {max_cluster:,}); use denominator='gamma' or a stricter threshold")
        order = np.argsort(self.labels, kind="stable")
        bounds = np.r_[0, np.flatnonzero(np.diff(self.labels[order])) + 1, self.n]
        rows, cols, vals = [], [], []
        single = []
        for a, b in zip(bounds[:-1], bounds[1:]):
            mem = order[a:b]
            if mem.size == 1:
                single.append(mem[0])
                continue
            V = (se + sf) * np.eye(mem.size) + sf * self.Koff[mem][:, mem].toarray()
            ev = linalg.eigvalsh(V)
            floor = 1e-3 * (se + sf if se + sf > 0 else 1.0)
            if ev.min() < floor:
                V += (floor - ev.min()) * np.eye(mem.size)
            Vi = linalg.inv(V)
            rr, cc = np.meshgrid(mem, mem, indexing="ij")
            rows.append(rr.ravel()); cols.append(cc.ravel()); vals.append(Vi.ravel())
        if single:
            s = np.asarray(single)
            rows.append(s); cols.append(s); vals.append(np.full(s.size, 1.0 / (se + sf)))
        return sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                                 shape=(self.n, self.n))


def as_relatedness(obj, n: int) -> Relatedness:
    """A :class:`Relatedness` from one, from (i, j, r) arrays, or from a
    scipy sparse (n, n) matrix of relatedness (its diagonal is ignored)."""
    if isinstance(obj, Relatedness):
        if obj.n != n:
            raise ValueError(f"relatedness is for {obj.n} samples, the panel has {n}")
        return obj
    if sparse.issparse(obj):
        if obj.shape != (n, n):
            raise ValueError(f"relatedness matrix shape {obj.shape} != ({n}, {n})")
        u = sparse.triu(obj, k=1).tocoo()
        lo = sparse.tril(obj, k=-1).tocoo()
        if lo.nnz and not u.nnz:
            u = lo
        return Relatedness(n, u.row, u.col, u.data, source="matrix")
    if isinstance(obj, (tuple, list)) and len(obj) == 3:
        return Relatedness(n, *obj, source="pairs")
    raise TypeError("relatedness must be a Relatedness, (i, j, r) arrays or a sparse matrix")


def pedigree_relatedness(father, mother, genotyped, *, min_r: float = 0.0,
                         inbreeding: bool = False) -> Relatedness:
    """Pedigree relatedness (the numerator relationship A = 2 x kinship) among
    the genotyped individuals, by Henderson's decomposition A = T D T'.

    ``father``, ``mother``: parent row index per individual (-1 unknown),
    parents listed before their children. ``genotyped``: row indices of the
    genotyped individuals, in the panel's sample order. Row i of T holds the
    gene-flow from each ancestor (1/2 per generation), so A_ij among the
    genotyped is a sparse product that is non-zero only for pairs with a
    common ancestor: cost O(sum over ancestors of their genotyped
    descendants squared), no genotype data. ``inbreeding=False`` uses
    D = 1, 3/4, 1/2 (no / one / two known parents); with ``True`` it uses
    the Meuwissen-Luo D from the recursively computed F (slower)."""
    father = np.asarray(father, dtype=np.int64)
    mother = np.asarray(mother, dtype=np.int64)
    N = father.size
    if mother.size != N:
        raise ValueError("father and mother differ in length")
    if np.any(father >= np.arange(N)) or np.any(mother >= np.arange(N)):
        raise ValueError("parents must be listed before their children")
    # rows of T as sparse dicts would be slow in Python for large N; build T row by row in CSR
    indptr = [0]
    indices: list = []
    data: list = []
    rows_cache = []
    F = np.zeros(N)
    d = np.empty(N)
    for i in range(N):
        acc: dict = {i: 1.0}
        for p in (father[i], mother[i]):
            if p >= 0:
                ri, rv = rows_cache[p]
                for a, v in zip(ri, rv):
                    acc[a] = acc.get(a, 0.0) + 0.5 * v
        nk = int(father[i] >= 0) + int(mother[i] >= 0)
        if inbreeding:
            fs = F[father[i]] if father[i] >= 0 else -1.0
            fm = F[mother[i]] if mother[i] >= 0 else -1.0
            d[i] = (1.0 if nk == 0 else 0.75 - 0.25 * (fs if father[i] >= 0 else fm) if nk == 1
                    else 0.5 - 0.25 * (fs + fm))
        else:
            d[i] = (1.0, 0.75, 0.5)[nk]
        ks = np.fromiter(acc.keys(), dtype=np.int64, count=len(acc))
        vs = np.fromiter(acc.values(), dtype=np.float64, count=len(acc))
        if inbreeding:
            F[i] = float(np.sum(vs * vs * d[ks])) - 1.0
        rows_cache.append((ks, vs))
        indices.append(ks); data.append(vs); indptr.append(indptr[-1] + ks.size)
    T = sparse.csr_matrix((np.concatenate(data), np.concatenate(indices), np.asarray(indptr)),
                          shape=(N, N))
    gen = np.asarray(genotyped, dtype=np.int64)
    Tg = T[gen]
    A = (Tg @ sparse.diags(d) @ Tg.T).tocoo()
    k = (A.row < A.col) & (A.data > min_r)
    return Relatedness(gen.size, A.row[k], A.col[k], A.data[k], source="pedigree")


def king_relatedness(gt, variants=None, *, min_phi: float = KING_THIRD_DEGREE,
                     min_maf: float = 0.05, block: int = 2048) -> Relatedness:
    """KING-robust kinship (Manichaikul et al. 2010, between-family
    estimator) of all sample pairs from the calls of ``variants`` (indices;
    default all), keeping pairs with phi > ``min_phi`` (default: 3rd degree
    and closer) as relatedness 2 phi. Dense O(n^2) accumulators: a
    prototype for n up to a few 10^4; biobank-scale KING runs elsewhere
    (KING, PLINK 2 --make-king-table) and is passed in as pairs."""
    n = gt.n_samples
    idx = np.arange(gt.n_variants) if variants is None else np.asarray(variants)
    hh = np.zeros((n, n))
    opp = np.zeros((n, n))
    nhet = np.zeros(n)
    used = 0
    for i in range(0, idx.size, block):
        take = idx[i:i + block]
        G = np.asarray(gt.G[:, take]).astype(np.int8)  # (n, k) calls, -1 missing
        called = G >= 0
        f = np.where(called, G, 0).sum(0) / np.maximum(2 * called.sum(0), 1)
        ok = np.minimum(f, 1 - f) >= min_maf
        G = G[:, ok]
        used += int(ok.sum())
        H = (G == 1).astype(np.float32)
        A = (G == 0).astype(np.float32)
        B = (G == 2).astype(np.float32)
        hh += H @ H.T
        ab = A @ B.T
        opp += ab + ab.T
        nhet += H.sum(1)
    nmin = np.maximum(np.minimum(nhet[:, None], nhet[None, :]), 1.0)
    phi = (hh - 2 * opp) / (2 * nmin) + 0.5 - (nhet[:, None] + nhet[None, :]) / (4 * nmin)
    iu = np.triu_indices(n, 1)
    v = phi[iu]
    s = v > min_phi
    rel = Relatedness(n, iu[0][s], iu[1][s], 2 * v[s], source=f"king({used} variants)")
    return rel


RELATEDNESS_OPTIONS = {"denominator": "exact", "n_gamma": 1000, "max_cluster": 5000,
                       "king_max_variants": 50_000, "king_min_maf": 0.05, "king_max_n": 50_000,
                       "king_noise_pairs": 1_000_000, "king_noise_sd_multiple": 4.0,
                       "kappa_floor": 0.2, "loco_refit": "auto",
                       "unrelated_min_phi": KING_THIRD_DEGREE, "unrelated_min_n": 100,
                       "king_max_abs_f": 0.5, "kinship_variants": None,
                       "kinship_max_model_overlap": 0.1, "structure_without_pcs": "warn",
                       "loco_refit_min_related": 0.01, "crossfit_folds": "random",
                       "carrier_families_p": 1e-4, "carrier_families_max": 3,
                       "carrier_families_max_maf": 0.05}

# KING-robust from a sample block against all samples, keeping pairs above a
# threshold: O(block x n) memory, never the n x n matrix.


def _king_indicators(G: np.ndarray):
    """float32 het / hom-ref / hom-alt indicators of int8 calls (-1 missing)."""
    return ((G == 1).astype(np.float32), (G == 0).astype(np.float32),
            (G == 2).astype(np.float32))


def _king_variants(gt, mean, exclude, min_maf: float, cap: int) -> np.ndarray:
    """Common variants (MAF >= min_maf from the prepared means) not in
    ``exclude``, at most ``cap`` evenly spread."""
    half = np.asarray(mean, dtype=np.float64) / 2.0
    ok = np.minimum(half, 1.0 - half) >= min_maf
    if exclude is not None:
        ok[np.asarray(exclude, dtype=np.int64)] = False
    idx = np.flatnonzero(ok)
    if idx.size > cap:
        idx = idx[np.unique(np.linspace(0, idx.size - 1, cap).round().astype(np.int64))]
    return idx


def king_pairs(gt, variants, threshold: float, *, block: int = 2048, chunk: int = 4096,
               noise_pairs: int = 1_000_000, noise_sd_multiple: float = 4.0, rng=None) -> dict:
    """KING-robust kinship (Manichaikul et al. 2010, between-family
    estimator) of every sample pair on ``variants``, blockwise: sample rows
    in blocks of ``block`` against all samples, variants in chunks of
    ``chunk`` (memory O(block x n + n x chunk)). The noise floor is the
    robust SD (1.4826 MAD) of phi over a seeded random subsample of
    ``noise_pairs`` pairs (default 1e6, or all pairs if fewer exist),
    taken from the same block products (no extra genotype pass, O(K)
    memory). The kept pairs have phi > max(threshold, noise_sd_multiple x
    SD). An empirical upper-quantile / Gaussian ratio on that subsample
    (``noise_tail_ratio``) is also returned, for a heavier-than-Gaussian
    tail check on the LOCO-refit trigger. Each row block decodes the KING
    variants once: ceil(n / block) reads of them, 3 n^2 k flops (float32
    GEMM; the symmetric half is not exploited).

    Returns {"i", "j", "r" (= 2 phi), "threshold", "threshold_requested",
    "noise_sd", "noise_median", "noise_tail_ratio", "noise_n_pairs",
    "n_variants", "seconds"}. Missing calls count as neither
    het nor hom (KING-robust ignores them; heterozygote counts per sample
    use the called variants only)."""
    import time

    t0 = time.perf_counter()
    rng = np.random.default_rng(0) if rng is None else rng
    n = gt.n_samples
    idx = np.asarray(variants, dtype=np.int64)
    if idx.size == 0:
        raise ValueError("KING: no common variants to compute kinship on")
    chunks = [idx[c:c + chunk] for c in range(0, idx.size, chunk)]

    def calls(take):
        return np.asarray(gt.G[:, take]).astype(np.int8, copy=False)

    # Random pairs for the noise floor, read off the block products (no
    # extra pass). Heterozygote counts complete during the first row block,
    # whose chunk sweep decodes every sample. Pairs above the requested
    # threshold are kept; the floor-raised threshold filters them at the end.
    k = min(int(noise_pairs), n * (n - 1) if n > 1 else 0)
    ra = rng.integers(0, n, k)
    rb = rng.integers(0, n, k)
    keep = ra != rb
    ra, rb = ra[keep], rb[keep]
    phi_r = np.empty(ra.size)
    nhet = np.zeros(n)
    hexp = np.zeros(n)  # expected heterozygotes per sample under HWE
    n_hom_alt = 0
    n_called = 0
    I, J, V = [], [], []
    for s in range(0, n, block):
        e = min(n, s + block)
        hh = np.zeros((e - s, n), dtype=np.float64)
        opp = np.zeros((e - s, n), dtype=np.float64)
        for take in chunks:
            Gc = calls(take)
            H, A, B = _king_indicators(Gc)
            if s == 0:
                nhet += H.sum(1)
                called = (Gc >= 0)
                f = np.where(called, Gc, 0).sum(0) / np.maximum(2 * called.sum(0), 1)
                hexp += called.astype(np.float32) @ (2 * f * (1 - f)).astype(np.float32)
                n_hom_alt += int(B.sum())
                n_called += int(called.sum())
            hh += H[s:e] @ H.T
            opp += A[s:e] @ B.T + B[s:e] @ A.T
        if s == 0:
            # Inbreeding F = 1 - observed / expected heterozygosity per
            # sample; KING-robust assumes outbred diploid calls (0/1/2).
            with np.errstate(divide="ignore", invalid="ignore"):
                F = 1.0 - nhet / hexp
            F_med = float(np.nanmedian(F)) if np.isfinite(F).any() else 0.0
            nhet = np.maximum(nhet, 1.0)
        nm = np.minimum(nhet[s:e, None], nhet[None, :])
        phi = (hh - 2 * opp) / (2 * nm) + 0.5 - (nhet[s:e, None] + nhet[None, :]) / (4 * nm)
        inb = (ra >= s) & (ra < e)
        phi_r[inb] = phi[ra[inb] - s, rb[inb]]
        r, c = np.nonzero(phi > threshold)
        r = r + s
        up = c > r
        I.append(r[up]); J.append(c[up]); V.append(phi[r[up] - s, c[up]])
    med = float(np.median(phi_r)) if phi_r.size else 0.0
    sd = float(1.4826 * np.median(np.abs(phi_r - med))) if phi_r.size else 0.0
    tail_ratio = _king_tail_ratio(phi_r, med, sd)
    thr = float(max(threshold, noise_sd_multiple * sd))
    I, J, V = np.concatenate(I), np.concatenate(J), np.concatenate(V)
    ok = V > thr
    I, J, V = I[ok], J[ok], V[ok]
    return {"i": I, "j": J, "r": 2.0 * V,
            "threshold": thr, "threshold_requested": float(threshold), "noise_sd": sd,
            "noise_median": med, "noise_tail_ratio": tail_ratio,
            "noise_n_pairs": int(phi_r.size),
            "n_variants": int(idx.size), "seconds": time.perf_counter() - t0,
            "median_F": F_med, "hom_alt_fraction": n_hom_alt / max(n_called, 1)}


def _king_tail_ratio(phi_r, med: float, sd: float) -> float:
    """Empirical upper-quantile / Gaussian ratio on a pair subsample.

    KING-robust noise is not exactly Gaussian (it depends on MAF, LD and
    ancestry). A ratio > 1 means a heavier-than-Gaussian upper tail; the
    LOCO-refit trigger then inflates σ̂ by this factor. The quantile is
    taken at a moderately central alpha (2%, ~2.05 σ under a Gaussian)
    so a few percent of true related pairs in the subsample cannot inflate
    it, while a heavier-than-Gaussian noise tail still can. Returns 1 when
    the subsample is too small or σ̂ is 0.
    """
    from scipy.stats import norm

    phi_r = np.asarray(phi_r, dtype=np.float64)
    if phi_r.size < 1_000 or sd <= 0.0:
        return 1.0
    alpha = 0.02
    q_emp = float(np.quantile(phi_r, 1.0 - alpha))
    z_a = float(norm.isf(alpha))
    den = z_a * sd
    if den <= 0.0:
        return 1.0
    return float((q_emp - med) / den)


def king_refit_trigger(n: int, i, j, phi, *, noise_sd: float, noise_median: float = 0.0,
                       relatedness_threshold: float = 0.025, min_related: float = 0.01,
                       expected_null: float = 1.0, noise_tail_ratio: float = 1.0,
                       n_variants: int | None = None) -> dict:
    """Whether KING-auto pairs are evidence of relatives for the LOCO refit.

    The pairs that enter Phi are not changed. For the trigger only, a pair
    counts as related if ``phi > max(relatedness_threshold, median + z σ_eff)``,
    with σ_eff = σ̂ · max(1, tail_ratio) (the empirical-to-Gaussian upper
    quantile ratio, so heavier-than-Gaussian tails raise the cutoff) and
    z the standard-normal upper quantile at which the expected number of
    null pairs among n(n-1)/2 is ``expected_null`` (default 1; z is about
    √(2 ln(n²/2))). The cutoff is capped at ``KING_SECOND_DEGREE`` (0.0884,
    the 2nd/3rd-degree boundary) so 2nd-degree relatives still count; if
    the uncapped noise cutoff would exceed the cap, a warning is issued
    that KING-auto has too few markers to distinguish 3rd-degree relatives
    from noise, and a pedigree or external pairs should be supplied. The
    refit itself is not disabled.

    Evidence then uses the same ``min_related`` sample-fraction rule as
    before (default 1%), applied to this quieter set. A handful of true
    close pairs below that fraction does not fire: the in-model leak is
    negligible there (≤ ~0.0014 at 1% related, ≤ 0.006 at 10%, in the
    CW mixed-family panels), so the 1.3-3.8x stage-1 cost of a refit is
    not worth it. n² growth of 4-SD-floor noise pairs is what this
    cutoff stops from tripping the 1% rule (SG p20k §4: 244 pairs, 1.33%
    of samples, on an unrelated n=20,000 panel).
    """
    from scipy.stats import norm

    n = int(n)
    i = np.asarray(i)
    j = np.asarray(j)
    phi = np.asarray(phi, dtype=np.float64)
    n_pairs_all = n * (n - 1) // 2
    sd = float(noise_sd)
    med = float(noise_median)
    tail_ratio = float(noise_tail_ratio) if np.isfinite(noise_tail_ratio) else 1.0
    sd_eff = sd * max(1.0, tail_ratio)
    cap = float(KING_SECOND_DEGREE)
    floor = float(relatedness_threshold)
    cutoff_uncapped = floor
    capped = False
    if n < 2 or n_pairs_all < 1:
        z = 0.0
        cutoff = floor
        n_exp = 0.0
    else:
        # P(Z > z) = expected_null / N, so E[null exceedances] = expected_null
        # under the (possibly tail-inflated) Gaussian.
        p = min(max(float(expected_null) / n_pairs_all, 0.0), 1.0)
        z = 0.0 if p >= 1.0 else (float("inf") if p <= 0.0 else float(norm.isf(p)))
        if not np.isfinite(z) or sd_eff <= 0.0:
            cutoff_uncapped = floor
        else:
            cutoff_uncapped = float(max(floor, med + z * sd_eff))
        capped = cutoff_uncapped > cap and floor < cap
        cutoff = cutoff_uncapped if not capped else cap
        z_eff = (cutoff - med) / sd_eff if sd_eff > 0.0 else float("inf")
        n_exp = float(n_pairs_all * norm.sf(z_eff)) if np.isfinite(z_eff) else 0.0
        if capped:
            markers = f" on {n_variants:,} variants" if n_variants else ""
            warnings.warn(
                f"relatedness='auto': the n-aware KING-auto refit-trigger cutoff would be "
                f"phi {cutoff_uncapped:.4f} (median + z·σ̂{markers}, σ̂={sd:.4f}"
                + (f", tail ratio {tail_ratio:.2f}" if tail_ratio > 1.0 else "")
                + f"), above KING's 2nd/3rd-degree bound ({cap:g}). The trigger is "
                f"capped at {cap:g} and may not distinguish 3rd-degree relatives from "
                "noise; pass a pedigree or external KING / GCTA pairs. The LOCO refit "
                "is not disabled",
                RuntimeWarning, stacklevel=2)
    ok = phi > cutoff if phi.size else np.zeros(0, dtype=bool)
    n_pairs = int(np.count_nonzero(ok))
    if n_pairs:
        n_rel = int(np.unique(np.r_[i[ok], j[ok]]).size)
    else:
        n_rel = 0
    frac = n_rel / n if n else 0.0
    evidence = n_rel >= float(min_related) * n if n else False
    reason = (f"{n_pairs} pair(s) above trigger cutoff {cutoff:.4f}"
              + (f" (capped from {cutoff_uncapped:.4f})" if capped else "")
              + f" (median {med:.4f} + z={z:.2f}·σ_eff, σ̂={sd:.4f}"
              + (f", tail ratio {tail_ratio:.2f}" if tail_ratio != 1.0 else "")
              + f"; expected null {n_exp:.2f} of {n_pairs_all:,} pairs; "
              f"{n_rel} samples, {frac:.2%} of n={n:,}; "
              f"{'evidence' if evidence else 'no evidence'} at min_related={min_related:g})")
    return {"noise_sd": sd, "noise_median": med, "noise_tail_ratio": tail_ratio,
            "sd_eff": sd_eff, "z": z, "cutoff": cutoff,
            "cutoff_uncapped": float(cutoff_uncapped), "capped": bool(capped),
            "n_pairs": n_pairs, "n_expected_null": n_exp, "n_related": n_rel,
            "related_fraction": frac, "evidence": bool(evidence), "reason": reason}


def align_pairs(id1, id2, r, sample_ids, *, source: str = "pairs") -> tuple:
    """Map ID pairs to row indices of ``sample_ids``. Pairs with a sample
    outside the analysis set are dropped and counted; a pair listed twice
    (either order) is an error. Returns (Relatedness, n_dropped)."""
    pos = {str(s): i for i, s in enumerate(np.asarray(sample_ids).astype(str))}
    id1 = np.asarray(id1).astype(str)
    id2 = np.asarray(id2).astype(str)
    a = np.array([pos.get(x, -1) for x in id1], dtype=np.int64)
    b = np.array([pos.get(x, -1) for x in id2], dtype=np.int64)
    ok = (a >= 0) & (b >= 0) & (a != b)
    dropped = int(np.count_nonzero(~ok))
    a, b, r = a[ok], b[ok], np.asarray(r, dtype=np.float64)[ok]
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    key = lo * len(pos) + hi
    if np.unique(key).size != key.size:
        raise ValueError(f"{source}: a sample pair is listed more than once")
    return Relatedness(len(pos), lo, hi, r, source=source), dropped


def read_kin0(path: str, sample_ids) -> tuple:
    """A KING ``.kin0`` / ``.kin`` or PLINK 2 ``.kin0`` table (header with
    ID1/IID1, ID2/IID2 and Kinship/KINSHIP); relatedness 2 x Kinship.
    Returns (Relatedness, n_dropped)."""
    with open(path) as fh:
        head = fh.readline().lstrip("#").split()
        rows = [ln.split() for ln in fh if ln.strip()]
    up = [h.upper() for h in head]

    def col(*names):
        for nm in names:
            if nm in up:
                return up.index(nm)
        raise ValueError(f"{path}: no column {' / '.join(names)} in the header")

    c1, c2, ck = col("IID1", "ID1"), col("IID2", "ID2"), col("KINSHIP")
    if not rows:
        return Relatedness(len(sample_ids), [], [], [], source=os.path.basename(path)), 0
    arr = np.array(rows, dtype=object)
    return align_pairs(arr[:, c1], arr[:, c2], 2.0 * arr[:, ck].astype(float), sample_ids,
                       source=os.path.basename(path))


def read_grm_sparse(prefix: str, sample_ids) -> tuple:
    """A GCTA sparse GRM (``prefix.grm.sp``: 0-based index pairs into
    ``prefix.grm.id`` and the relatedness; diagonal entries are ignored).
    Returns (Relatedness, n_dropped)."""
    if prefix.endswith(".grm.sp"):
        prefix = prefix[:-len(".grm.sp")]
    ids = np.array([ln.split()[1] for ln in open(prefix + ".grm.id") if ln.strip()])
    sp = np.loadtxt(prefix + ".grm.sp", ndmin=2)
    if sp.size == 0:
        return Relatedness(len(sample_ids), [], [], [], source=os.path.basename(prefix)), 0
    a, b = sp[:, 0].astype(np.int64), sp[:, 1].astype(np.int64)
    off = a != b
    return align_pairs(ids[a[off]], ids[b[off]], sp[off, 2], sample_ids,
                       source=os.path.basename(prefix) + ".grm.sp")


def unrelated_set(rel: Relatedness, min_phi: float = 0.0) -> np.ndarray:
    """A maximal unrelated set (no kept pair with phi = r / 2 > ``min_phi``;
    0: every listed pair): repeatedly keep the open sample with the fewest
    open relatives and drop those relatives (ties: the lower index). Returns
    the kept row indices, sorted.

    Minimum-degree selection, not maximum-degree removal: on the SG's
    isolate panel (2,000 sampled first cousins of each other, 3.6 cousins
    each, phi 1/16) removal kept 775 of an exact maximum of 822-824 (MILP),
    selection 795-828. A maximum independent set is NP-hard in general;
    this greedy is within a few % of it on the simulated pedigrees."""
    import heapq

    if min_phi > 0:
        rel = rel.threshold_phi(min_phi)
    K = rel.Koff.tocsr()
    ind, ptr = K.indices, K.indptr
    deg = np.diff(ptr).astype(np.int64)
    state = np.zeros(rel.n, dtype=np.int8)  # 0 open, 1 kept, 2 dropped
    state[deg == 0] = 1
    heap = [(int(d), i) for i, d in enumerate(deg) if d > 0]
    heapq.heapify(heap)
    while heap:
        d, i = heapq.heappop(heap)
        if state[i] or d != deg[i]:
            continue
        state[i] = 1
        for j in ind[ptr[i]:ptr[i + 1]]:
            if state[j] == 0:
                state[j] = 2
                for k in ind[ptr[j]:ptr[j + 1]]:
                    if state[k] == 0:
                        deg[k] -= 1
                        heapq.heappush(heap, (int(deg[k]), int(k)))
    return np.flatnonzero(state == 1)


def carrier_families(lg, rel: Relatedness, idx, threshold: float) -> np.ndarray:
    """(iv) Distinct families among the minor-allele carriers of variants
    ``idx`` (``lg``: their LocoGenotypes): connected components of the
    carriers' own pairs with phi > ``threshold``. Unrelated carriers count
    one each; carriers linked only through non-carriers count apart (in a
    cousin-rich isolate the full pair graph is one giant component, which
    would call every variant one family). One decode of ``idx``."""
    from scipy.sparse.csgraph import connected_components

    from mixmogam._standardize import decode_table

    idx = np.asarray(idx, dtype=np.int64)
    out = np.zeros(idx.size, dtype=np.int64)
    if not idx.size:
        return out
    K = rel.threshold_phi(threshold).Koff.tocsr()
    minor_is_counted = lg.mean[idx] / 2.0 <= 0.5
    lut = np.where(minor_is_counted[:, None], np.array([0.0, 1.0, 1.0, 0.0]),
                   np.array([1.0, 1.0, 0.0, 0.0])).astype(np.float32)
    step = 256  # decode in row tiles: memory O(256 n)
    for a in range(0, idx.size, step):
        b = min(idx.size, a + step)
        C = np.empty((b - a, lg.n), dtype=np.float32)
        decode_table(lg.gt.G, idx[a:b], lut[a:b], C, 1, compiled=lg._compiled)
        for r in range(b - a):
            car = np.flatnonzero(C[r] > 0)
            out[a + r] = (car.size if car.size <= 1 else
                          connected_components(K[car][:, car], directed=False)[0])
    return out


def covariate_factor(rel: Relatedness, Q: np.ndarray, n_eff: float) -> tuple:
    """c_X = n_eff / (n - tr(Q' Phi Q)), Phi = I + offdiag(r), Q the
    orthonormal covariate basis: tr(Q'PhiQ) = q + sum_{i != k} r_ik (QQ')_ik.
    Returns (c_X, tr(Q' Phi Q))."""
    q = Q.shape[1]
    u = rel.upper.tocoo()
    tr = float(q + 2.0 * np.sum(u.data * np.einsum("ij,ij->i", Q[u.row], Q[u.col]))) if q else 0.0
    den = rel.n - tr
    return (float(n_eff / den) if den > 0 else 1.0), tr


def resolve_relatedness(spec, gt, *, mean, model_idx, threshold: float, options: dict,
                        rng=None) -> tuple:
    """The analysis' Relatedness from ``relatedness=`` (see :func:`hratt`),
    with its record. Precedence is the caller's: a pedigree / register
    Relatedness or pairs as given; a KING ``.kin0`` or GCTA ``.grm.sp`` path;
    ``"auto"``: KING-robust on common non-model variants (all variants
    without a subset), capped and noise-floor guarded."""
    n = gt.n_samples
    info = {"dropped_pairs": 0}
    if isinstance(spec, str) and spec == "auto":
        if n > options["king_max_n"]:
            raise ValueError(
                f"relatedness='auto' computes KING-robust kinship over all {n:,}^2/2 pairs; "
                f"above {options['king_max_n']:,} samples pass external pairs instead (KING "
                "--related, PLINK 2 --make-king-table, or GCTA --make-bK-sparse) as a "
                ".kin0 / .grm.sp path or (i, j, r) arrays")
        idx = _king_variants(gt, mean, model_idx, options["king_min_maf"],
                             options["king_max_variants"])
        on = "non-model variants" if model_idx is not None else "all variants"
        if idx.size < 1000 and model_idx is not None:
            warnings.warn(f"relatedness='auto': only {idx.size} common non-model variants; "
                          "KING on few variants has a high noise floor", RuntimeWarning,
                          stacklevel=3)
        res = king_pairs(gt, idx, threshold, noise_pairs=options["king_noise_pairs"],
                         noise_sd_multiple=options["king_noise_sd_multiple"], rng=rng)
        F, homalt = res["median_F"], res["hom_alt_fraction"]
        coding = ("haploid or 0/1-coded calls (no homozygous-alternative calls)" if homalt == 0
                  else f"median inbreeding F = {F:+.2f}")
        if abs(F) > options["king_max_abs_f"] or homalt == 0:
            raise ValueError(
                f"relatedness='auto': {coding}; KING-robust assumes outbred diploid 0/1/2 "
                "genotypes and is invalid here (it calls most pairs related). Pass a pedigree "
                "Relatedness or external pairs, or raise relatedness_options['king_max_abs_f']")
        if abs(F) > 0.05:
            warnings.warn(f"relatedness='auto': {coding}; KING-robust kinship is biased for "
                          "inbred or excess-heterozygosity samples", RuntimeWarning, stacklevel=3)
        if res["i"].size > 5 * n:
            warnings.warn(f"relatedness='auto': {res['i'].size:,} pairs for {n:,} samples (> 5 "
                          "per sample) is implausible for close relatives; check the genotype "
                          "coding and population structure, or pass a pedigree or external "
                          "pairs", RuntimeWarning,
                          stacklevel=3)
        if res["threshold"] > threshold:
            # Population structure inflates the random-pair spread (the SG's
            # Fst 0.05 panel without PCs: SD 0.039, floor 0.158), and above
            # 0.0884 the floor drops 2nd-degree relatives from kappa_g.
            lost = dropped_degree(res["threshold"])
            warnings.warn(
                f"relatedness='auto': the KING noise SD on {res['n_variants']:,} variants is "
                f"{res['noise_sd']:.4f}, so pairs with phi <= {res['threshold']:.4f} "
                f"({options['king_noise_sd_multiple']:g} SD) are dropped instead of <= "
                f"{threshold:g}"
                + (f"; this floor is above KING's 2nd-degree bound (0.0884), so {lost} "
                   "can be dropped too, usually because population structure inflates the "
                   "random-pair spread. Pass pedigree or external (PC-aware) kinship pairs"
                   if lost else
                   "; pass pedigree or external KING pairs for distant relatives"),
                RuntimeWarning, stacklevel=3)
        rel = Relatedness(n, res["i"], res["j"], res["r"], source=f"king-auto ({on})")
        info.update({"source": "king-auto", "king": {k: res[k] for k in (
            "threshold", "threshold_requested", "noise_sd", "noise_median",
            "noise_tail_ratio", "noise_n_pairs", "n_variants",
            "seconds", "median_F", "hom_alt_fraction")},
            "king_on": on})
        return rel, info
    if isinstance(spec, (str, os.PathLike)):
        path = os.fspath(spec)
        if path.endswith(".grm.sp") or os.path.exists(path + ".grm.sp"):
            rel, dropped = read_grm_sparse(path, gt.sample_ids)
        else:
            rel, dropped = read_kin0(path, gt.sample_ids)
        rel = rel.threshold_phi(threshold)
        info.update({"source": rel.source, "dropped_pairs": dropped})
    else:
        rel = as_relatedness(spec, n)
        info["source"] = rel.source or "pairs"
    info.update(_kinship_overlap(gt, options, model_idx))
    return rel, info


def _kinship_overlap(gt, options: dict, model_idx) -> dict:
    """P4: external kinship estimated on the stage-1 model SNPs under-corrects
    (KING on model SNPs left the null mean chi-square at up to 1.09 in the
    SG's grid: the fit anti-correlates the residuals of pairs that look
    related on its own SNPs). With ``options["kinship_variants"]`` (the
    variants the pairs were estimated on: indices, a mask, IDs, or a path
    to a file of IDs, e.g. PLINK ``.prune.in``) the overlap with the model
    SNPs is recorded and a warning issued above ``kinship_max_model_overlap``
    of them. Without it nothing can be checked (recorded as None)."""
    kv = options.get("kinship_variants")
    if kv is None:
        return {"kinship_model_overlap": None}
    from mixmogam._thin import _resolve_indices

    if isinstance(kv, (str, os.PathLike)):
        with open(os.fspath(kv)) as fh:
            kv = np.array([ln.split()[0] for ln in fh if ln.strip()])
    kidx = _resolve_indices(gt, kv, "kinship_variants", strict=False)
    model = (np.arange(gt.n_variants) if model_idx is None
             else np.asarray(model_idx, dtype=np.int64))
    if kidx.size == 0:
        warnings.warn("relatedness_options['kinship_variants'] matched no panel variant; the "
                      "kinship-vs-model-SNP overlap was not checked", RuntimeWarning, stacklevel=4)
        return {"kinship_model_overlap": None, "kinship_variants_matched": 0}
    frac = float(np.isin(kidx, model).mean())
    # Full stage 1 (no subset): every variant is a model SNP and KING on the
    # whole panel calibrated in the SG grid; only subsets are guarded.
    if model_idx is not None and frac > options["kinship_max_model_overlap"]:
        warnings.warn(
            f"relatedness: {frac:.0%} of the {kidx.size:,} variants the kinship was estimated "
            "on are stage-1 model SNPs; kinship from model SNPs biases kappa_g down (null mean "
            "chi-square up to 1.09 in simulation, worst with small subsets). Estimate it on "
            "common variants that are not model SNPs, or use a pedigree", RuntimeWarning,
            stacklevel=4)
    return {"kinship_model_overlap": frac, "kinship_variants_matched": int(kidx.size)}


def gls_stats(lg, Wp: np.ndarray, rel: Relatedness, *, denominator: str = "exact",
              n_gamma: int = 1000, max_cluster: int = 5000, rng=None) -> dict:
    """Fix (b): GLS score tests of every variant of ``lg`` against its group's
    column of the projected LOCO residuals ``Wp`` under Cov(w_g) = s_e I + s_f R
    (see the module docstring). Returns chi2, num (z'P w), den (z'P z) and the
    per-group fits."""
    Q = lg.Q if lg.Q.shape[1] else np.ones((lg.n, 1)) / np.sqrt(lg.n)
    G = Wp.shape[1]
    sf, se = rel.family_vc(Wp)
    vinv, b, QVQi = [], np.empty_like(Wp), []
    for g in range(G):
        Vi = rel.vinv(float(sf[g]), float(se[g]), max_cluster=max_cluster)
        ViQ = Vi @ Q
        M = linalg.inv(Q.T @ ViQ)
        w = Wp[:, g]
        Viw = Vi @ w
        b[:, g] = Viw - ViQ @ (M @ (Q.T @ Viw))   # P w
        vinv.append((Vi, ViQ, M))
        QVQi.append(M)
    m = lg.m
    num = np.zeros(m)
    den = np.full(m, np.nan)
    tile = max(1, (16 * 1024**2) // max(8 * lg.n, 1))
    exact = denominator == "exact"
    for idx, g, Z64 in lg.wide_slices(tile):
        num[idx] = Z64 @ b[:, g]
        if exact:
            Vi, ViQ, M = vinv[g]
            VZ = (Vi @ Z64.T)                      # (n, k)
            c = ViQ.T @ Z64.T                      # (q, k) = Q'V^-1 z
            den[idx] = np.einsum("ij,ji->i", Z64, VZ) - np.einsum("ik,ij,jk->k", c, M, c)
    gamma = None
    if not exact:
        rng = rng or np.random.default_rng(0)
        cal = np.sort(rng.choice(m, min(int(n_gamma), m), replace=False))
        Zc = lg.rows(cal)                          # projected rows
        gamma = np.empty(G)
        for g in range(G):
            Vi, ViQ, M = vinv[g]
            VZ = Vi @ Zc.T
            c = ViQ.T @ Zc.T
            zpz = np.einsum("ij,ji->i", Zc, VZ) - np.einsum("ik,ij,jk->k", c, M, c)
            gamma[g] = float(np.mean(zpz / lg.zz[cal]))
        den = gamma[lg.groups] * lg.zz
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = num * num / den
    chi2[~(den > 0)] = np.nan
    return {"chi2": chi2, "num": num, "den": den, "s_f": sf, "s_e": se, "gamma": gamma}
