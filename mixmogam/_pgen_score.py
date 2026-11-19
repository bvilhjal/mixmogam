"""HRATT score products on PGEN's lists of exceptions to a homozygous call."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np

from mixmogam._fast import HAS_NUMBA
from mixmogam._packed import PackedCalls
from mixmogam.genotypes import Genotypes
from mixmogam.io._hardcalls import block_variants
from mixmogam.io.pgen import _calls, _metadata, _pgen_reader
from mixmogam.results import GwasResult


def _moments(indices, values, flip, remap, Q, q_sum, W, w_sum, scores):
    """Called-sample moments, projected norm and all traits' score numerators.

    The modal contribution uses a sum over samples prepared once. Each
    exception replaces that sample's contribution. Missing calls have a
    standardized value of zero. Only scores use the float32 decoding table.
    """
    if indices.size != values.size:
        raise ValueError("PGEN sparse hard-call lengths disagree")
    modal = 2 if flip else 0
    counts = np.zeros(4, dtype=np.int64)
    counts[modal] = Q.shape[0]
    for j in range(indices.size):
        value = values[j]
        if indices[j] >= remap.size or (value != -9 and not 0 <= value <= 2):
            raise ValueError("invalid PGEN sparse hard call")
        if remap[indices[j]] >= 0:
            code = 3 if value == -9 else (2 - value if flip else value)
            counts[modal] -= 1
            counts[code] += 1
    called = max(1, Q.shape[0] - counts[3])
    mean = (counts[1] + 2 * counts[2]) / called
    raw = np.array([-mean, 1.0 - mean, 2.0 - mean, 0.0])
    sumsq = 0.0
    for code in range(3):
        sumsq += counts[code] * raw[code] * raw[code]
    sd = np.sqrt(sumsq / called)
    divisor = sd if sd else 1.0
    table = (raw / divisor).astype(np.float32).astype(np.float64)
    projection = raw[modal] * q_sum
    for k in range(scores.size):
        scores[k] = table[modal] * w_sum[k]
    for j in range(indices.size):
        row = remap[indices[j]]
        if row < 0:
            continue
        value = values[j]
        code = 3 if value == -9 else (2 - value if flip else value)
        dq = raw[code] - raw[modal]
        dw = table[code] - table[modal]
        for k in range(Q.shape[1]):
            projection[k] += dq * Q[row, k]
        for k in range(scores.size):
            scores[k] += dw * W[row, k]
    norm = sumsq / (divisor * divisor)
    for k in range(projection.size):
        norm -= (projection[k] / divisor) ** 2
    return mean, sd, norm


if HAS_NUMBA:
    from numba import njit
    _moments = njit(cache=True)(_moments)


class _Scores:
    def __init__(self, runner):
        head = runner.fits[0]
        self.Q = head._panel["Q"]
        self.q_sum = self.Q.sum(axis=0)
        n, groups = head._state["resid"].shape
        self.W = np.empty((groups, n, len(runner.fits)))
        for k, (fit, cache) in enumerate(zip(runner.fits, runner.caches)):
            residual = fit._state["resid"]
            wp = self.W[:, :, k].T
            if "Wp" in cache:
                wp[:] = cache["Wp"]
            else:
                correction = (self.Q @ (self.Q.T @ residual) if self.Q.shape[1]
                              else residual.mean(axis=0))
                np.subtract(residual, correction, out=wp)
            cache["Wp"] = wp
            cache["wn"] = np.einsum("ij,ij->j", residual, residual)
        self.w_sum = self.W.sum(axis=1)

    def scan(self, runner, reader, start, meta, alleles, remap):
        from mixmogam.twostep import _hratt_test

        indices = np.flatnonzero(alleles >= 0)
        groups = runner._groups(0, meta["chromosome"][indices])
        mean, sd, zz = np.empty((3, indices.size))
        num = np.empty((indices.size, len(runner.fits)))
        for i, j in enumerate(indices):
            allele, group = int(alleles[j]), groups[i]
            rows, calls = reader.read_sparse8(start + j, allele_idx=allele)
            mean[i], sd[i], zz[i] = _moments(rows, calls, allele == 0, remap,
                self.Q, self.q_sum, self.W[group], self.w_sum[group], num[i])
        # A nearly exact projection can change sign with summation order.
        # Let the same dense kernel decide its rank and finite-value mask.
        tol = np.sqrt(np.finfo(float).eps)
        keep = zz > tol * self.Q.shape[0]
        # The historical standard error is |beta|/sqrt(chi2). Exactly
        # zero and a dot product rounded just above zero have different masks.
        for k, cache in enumerate(runner.caches):
            bound = np.sqrt(self.Q.shape[0] * cache["wn"][groups])
            keep &= np.abs(num[:, k]) > tol * bound
        keep |= sd == 0
        if not keep.all():
            indices, groups = indices[keep], groups[keep]
            mean, sd, zz, num = mean[keep], sd[keep], zz[keep], num[keep]
        if not indices.size:
            return indices, None
        data = SimpleNamespace(**{name: values[indices] for name, values in meta.items()})
        lg = SimpleNamespace(mean=mean, sd=sd, n=self.Q.shape[0],
                             n_groups=self.W.shape[0], groups=groups)
        results = []
        for k, (fit, cache) in enumerate(zip(runner.fits, runner.caches)):
            with np.errstate(divide="ignore", invalid="ignore"):
                chi2 = fit._setup.n_eff * num[:, k] ** 2 / (zz * cache["wn"][groups])
            chi2[~(zz > 0)] = np.nan
            stats = {"chi2": chi2, "num": num[:, k], "zz": zz}
            results.append(_hratt_test(fit._state, replace(fit._setup, lg=lg), data, rsU=stats))
        return indices, results


def pgen_blocks(runner, prefix, packed):
    eligible = (HAS_NUMBA and len(runner.shared) == 1
                and runner.fits[0]._panel["dtype"] == "float32"
                and all(not f._state["binary"] and not f._state["weighted"]
                        and f._state["rel"] is None and f._state["correction"] != "gls"
                        for f in runner.fits))
    with _pgen_reader(prefix, None) as (reader, pvar, ids, m, autosomes):
        eligible = eligible and hasattr(reader, "has_sparse") and hasattr(reader, "read_sparse8")
        step = block_variants(len(ids), runner.chunk_bytes)
        remap = None
        for start in range(0, m, step):
            stop = min(start + step, m)
            meta = _metadata(pvar, start, stop, autosomes)
            sparse = None
            dense_idx = slice(None)
            if eligible:
                alleles = np.array([1 if reader.has_sparse(j, allele_idx=1) else
                                    0 if reader.has_sparse(j, allele_idx=0) else -1
                                    for j in range(start, stop)], dtype=np.int8)
                if np.any(alleles >= 0):
                    if remap is None:
                        rows = runner._sample_index(ids)[0]
                        if rows is None:
                            remap = np.arange(len(ids))
                        else:
                            remap = np.full(len(ids), -1, dtype=np.int64)
                            remap[rows] = np.arange(rows.size)
                    if runner.pgen_scores is None:
                        runner.pgen_scores = _Scores(runner)
                    used, sparse = runner.pgen_scores.scan(runner, reader, start, meta, alleles, remap)
                    if used.size == stop - start:
                        yield sparse
                        continue
                    if used.size:
                        dense_idx = np.setdiff1d(np.arange(stop - start), used, assume_unique=True)
            calls = _calls(reader, start, stop, len(ids),
                           indices=None if sparse is None else start + dense_idx)
            if packed:
                calls = PackedCalls.pack(calls)
            gt = Genotypes._from_validated(calls, sample_ids=ids,
                **{key: value[dense_idx] for key, value in meta.items()})
            for dense in runner._call_blocks(gt):
                if sparse is None:
                    yield dense
                else:
                    order = np.argsort(np.r_[used, dense_idx])
                    yield [GwasResult.concat([s, d]).take(order) for s, d in zip(sparse, dense)]
