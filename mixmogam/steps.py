"""HRATT in two steps: fit step 1 once, test any genotypes in step 2.

:func:`mixmogam.twostep.hratt` fits the LOCO polygenic scores (step 1) and
tests every variant of the same panel (step 2) in one call. At biobank scale
the two belong apart, as in REGENIE and LDAK-KVIK:

- :func:`hratt_step1` fits step 1 on a panel (typically the model SNPs'
  hard calls) and returns a :class:`HrattStep1`, which :meth:`HrattStep1.save`
  writes to one ``.npz`` file. What step 2 needs is small: the LOCO
  residuals (n x groups), lambda or kappa_g, and the working model.
- :func:`hratt_step2` tests any genotypes against it: other files (one per
  chromosome, say), samples in any order (step 1's are selected), hard calls
  or imputed dosages (:class:`~mixmogam.genotypes.Dosages`), read in blocks
  of bounded memory, for several traits at once, writing results as it goes.

Every step-2 statistic is per variant given step 1, so hard calls get the
values :func:`~mixmogam.twostep.hratt` gives the same variants, whatever
the blocks (up to the summation order of one product, ~1e-16 relative,
when the storage differs). Dosages are standardized and tested as values,
in float64; integer dosages give the hard-call statistics up to float32
rounding of the stored calls (~1e-7 relative).
"""

from __future__ import annotations

from dataclasses import fields, replace
import json
import warnings
from types import SimpleNamespace
from typing import Optional

import numpy as np

from mixmogam.genotypes import Dosages, Genotypes
from mixmogam.results import GwasResult

__all__ = ["HrattStep1", "hratt_step1", "hratt_step2"]

FORMAT = "mixmogam-hratt-step1"
FORMAT_VERSION = 1
# Bytes of genotype storage (dosages: float64 values) per step-2 block.
CHUNK_BYTES = 256 * 1024**2


class HrattStep1:
    """A fitted HRATT step 1: what step 2 tests variants against.

    Made by :func:`hratt_step1`; :meth:`save` and :meth:`load` keep it in a
    file. ``sample_ids`` are the analysed samples in model order,
    ``loco_groups`` the chromosome groups (step 2 tests a variant against
    its chromosome's group), ``extra`` the step-1 records that results carry
    (h2, alpha, structure, lambda, relatedness, ...).

    In the process that fitted it, it also keeps the step-1 panel's prepared
    genotypes, so step 2 on that panel (the same ``Genotypes`` object) skips
    a pass; a loaded one does not. ``del`` it to release them.
    """

    def __init__(self, state: dict, setup, panel: dict, cache: Optional[dict] = None):
        self._state, self._setup, self._panel, self._cache = state, setup, panel, cache

    @property
    def sample_ids(self) -> np.ndarray:
        return self._panel["sample_ids"]

    @property
    def loco_groups(self) -> list:
        return list(self._panel["labels"])

    @property
    def trait(self) -> str:
        return self._state["trait"]

    @property
    def extra(self) -> dict:
        return dict(self._state["extra"])

    @property
    def n(self) -> int:
        return int(self._panel["sample_ids"].size)

    def __repr__(self):
        return (f"HrattStep1(trait={self.trait!r}, n={self.n}, "
                f"groups={len(self._panel['labels'])})")

    # ------------------------------------------------------------------
    def save(self, path) -> None:
        """Write step 1 to ``path`` (NumPy ``.npz``; arrays and a JSON
        record, no pickles), to be read by :meth:`load`."""
        from mixmogam import __version__

        arrays: dict = {}
        setup = {f.name: getattr(self._setup, f.name) for f in fields(self._setup)
                 if f.name != "lg"}
        record = {"format": FORMAT, "format_version": FORMAT_VERSION, "mixmogam": __version__,
                  "state": _encode(self._state, arrays, "state"),
                  "setup": _encode(setup, arrays, "setup"),
                  "panel": _encode(self._panel, arrays, "panel")}
        arrays["__record__"] = np.frombuffer(json.dumps(record).encode(), dtype=np.uint8)
        with open(path, "wb") as fh:
            np.savez(fh, **arrays)

    @classmethod
    def load(cls, path) -> "HrattStep1":
        """Read a step 1 written by :meth:`save`."""
        from mixmogam.twostep import _Setup

        with np.load(path, allow_pickle=False) as data:
            arrays = {key: data[key] for key in data.files}
        record = json.loads(arrays.pop("__record__").tobytes().decode())
        if record.get("format") != FORMAT or record.get("format_version") != FORMAT_VERSION:
            raise ValueError(f"{path} is not a mixmogam HRATT step-1 file "
                             f"(format {FORMAT!r} version {FORMAT_VERSION})")
        state = _decode(record["state"], arrays)
        setup = _Setup(lg=None, **_decode(record["setup"], arrays))
        panel = _decode(record["panel"], arrays)
        return cls(state, setup, panel)


def hratt_step1(y, gt, X=None, **options) -> HrattStep1:
    """Fit HRATT's step 1 (:func:`mixmogam.twostep.hratt`, same arguments)
    and return it for :func:`hratt_step2`.

    The structure test, relatedness, variance components, prior selection,
    LOCO scores and their calibration (lambda, kappa_g) are those of
    :func:`~mixmogam.twostep.hratt` on ``gt``. ``model_snps`` works as there:
    step 1 then fits on the subset of ``gt``. ``denominator="spectral"`` has
    per-variant factors only for step-1 variants and is not available.
    ``relatedness_correction="gls"`` with ``relatedness_options["denominator"]
    = "gamma"`` fixes gamma here, on ``gt``'s variants.
    """
    from mixmogam import twostep

    if options.get("denominator", "constant") == "spectral":
        raise ValueError("denominator='spectral' has per-variant factors only for the step-1 "
                         "panel; use mixmogam.twostep.hratt for it")
    state, cache = twostep._hratt_stage1(y, gt, X, **options)
    st = cache["st_test"]
    test_cache: dict = {}
    if state["correction"] == "gls" and state["rel_opt"]["denominator"] == "gamma":
        from mixmogam._relatedness import gls_gamma, gls_prepare

        rel, opt = state["rel"], state["rel_opt"]
        Wp = st.lg.project(state["resid"])
        Q = st.lg.Q if st.lg.Q.shape[1] else np.ones((st.lg.n, 1)) / np.sqrt(st.lg.n)
        vinv, b = gls_prepare(Q, Wp, rel, state["extra"]["relatedness"]["s_f"],
                              state["extra"]["relatedness"]["s_e"], opt["max_cluster"])
        seed = 0 if state["random_state"] is None else state["random_state"]
        state["gls_gamma"] = gls_gamma(st.lg, vinv, opt["n_gamma"],
                                       np.random.default_rng([seed, 0x7663]))
        test_cache.update({"Wp": Wp, "gls": (vinv, b)})
    extra = state["extra"]
    if any(extra.get(key) is False for key in ("cv_converged", "loco_converged")):
        warnings.warn("variational fit did not converge; inspect cv_converged/loco_converged "
                      "before using results", RuntimeWarning, stacklevel=2)
    panel = {"sample_ids": np.asarray(gt.sample_ids), "labels": list(st.labels),
             "Q": st.lg.Q, "n_groups": int(st.lg.n_groups), "block": int(st.lg.block),
             "dtype": st.lg.dtype.name}
    return HrattStep1(state, replace(st, lg=None), panel,
                      cache={"gt": gt, "test": test_cache, **cache})


def hratt_step2(step1, genotypes, *, n_threads: int = 1, chunk_bytes: float = CHUNK_BYTES,
                out=None, genome_build: Optional[str] = None):
    """Test variants against a fitted step 1 (HRATT's step 2).

    ``step1``: a :class:`HrattStep1` (or a list of them: several traits,
    whose blocks are read once). ``genotypes``: a
    :class:`~mixmogam.genotypes.Genotypes`, a
    :class:`~mixmogam.genotypes.Dosages`, or an iterable of them (files per
    chromosome; :func:`mixmogam.io.bgen.iter_bgen_dosages`). Each must hold
    every step-1 sample (others are ignored; any order) and only variants on
    chromosomes of step 1's LOCO groups. Variants are read in blocks of about
    ``chunk_bytes`` of storage (a memory-mapped bed file is read once, with
    drop-behind), so memory is bounded by one block and the per-sample
    arrays, whatever the number of variants.

    Hard calls get the statistics of :func:`~mixmogam.twostep.hratt`.
    Dosages are tested for unweighted quantitative traits with the lambda or
    kappa_g correction (binary, weighted and GLS tests need hard calls); the
    carrier-family columns of relatedness runs are NaN for them.

    Returns a :class:`~mixmogam.results.GwasResult` (a list, for a list of
    step 1s). With ``out`` (a path, or one per step 1) the blocks are
    written as they are tested and a summary ``{"variants": count, "out":
    paths}`` is returned instead. ``genome_build`` declares the results'
    coordinate assembly (no liftover).
    """
    single = isinstance(step1, HrattStep1)
    fits = [step1] if single else list(step1)
    if not fits or not all(isinstance(f, HrattStep1) for f in fits):
        raise TypeError("step1 must be a HrattStep1 or a list of them")
    def finish(res):
        res.genome_build = genome_build
        return res

    paths = None
    if out is not None:
        paths = [out] if isinstance(out, (str, bytes)) or hasattr(out, "__fspath__") else list(out)
        if len(paths) != len(fits):
            raise ValueError("out needs one path per step 1")
    # The step-1 panel itself, in the fitting process: hratt()'s exact pass.
    if (isinstance(genotypes, Genotypes) and all(f._cache is not None and f._cache["gt"] is genotypes
                                                for f in fits)):
        from mixmogam.twostep import _hratt_test

        results = [finish(_hratt_test(f._state, f._cache["st_test"], genotypes, n_threads,
                                      rsU=f._cache["rsU"], U=f._cache["U"],
                                      cache=dict(f._cache["test"]))) for f in fits]
        return _finish(results, paths, single)
    sources = [genotypes] if isinstance(genotypes, (Genotypes, Dosages)) else genotypes
    runner = _Runner(fits, n_threads, float(chunk_bytes))
    written = 0
    parts = [[] for _ in fits]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        for source in sources:
            for block in runner.blocks(source):
                for k, res in enumerate(block):
                    finish(res)
                    if paths is None:
                        parts[k].append(res)
                    else:
                        res.write_csv(paths[k], append=written > 0)
                written += len(block[0])
    seen = set()
    for w in caught:  # each block repeats the same warnings: show each once
        key = (w.category, str(w.message))
        if key not in seen:
            seen.add(key)
            warnings.warn_explicit(w.message, w.category, w.filename, w.lineno)
    if written == 0:
        raise ValueError("no variants to test")
    if paths is not None:
        return {"variants": written, "out": paths}
    results = [GwasResult.concat(p) for p in parts]
    return results[0] if single else results


def _finish(results, paths, single):
    if paths is None:
        return results[0] if single else results
    for res, path in zip(results, paths):
        res.write_csv(path)
    return {"variants": len(results[0]), "out": paths}


class _Runner:
    """Step 2's block loop. Each block of a source is read once; step 1s
    fitted on the same samples, covariates, row scale and LOCO groups share
    its sample selection and prepared genotypes (the moments are the costly
    pass), others get their own."""

    def __init__(self, fits, n_threads: int, chunk_bytes: float):
        from mixmogam._standardize import _validate_n_threads

        self.fits = fits
        self.n_threads = _validate_n_threads(n_threads)
        self.chunk_bytes = max(chunk_bytes, 1.0)
        self.caches = [{} for _ in fits]
        self.shared = []  # lists of fit indices
        for k, f in enumerate(fits):
            for members in self.shared:
                g = fits[members[0]]
                if (np.array_equal(f._panel["sample_ids"], g._panel["sample_ids"])
                        and _labels(f) == _labels(g)
                        and f._panel["dtype"] == g._panel["dtype"]
                        and f._panel["block"] == g._panel["block"]
                        and np.array_equal(f._panel["Q"], g._panel["Q"])
                        and _same(f._setup.s, g._setup.s)):
                    members.append(k)
                    break
            else:
                self.shared.append([k])
        self.ids = [np.asarray(fits[m[0]]._panel["sample_ids"]).astype(str) for m in self.shared]
        self.group_of = [{_chrom_key(c): g for g, label in enumerate(fits[m[0]]._panel["labels"])
                          for c in label} for m in self.shared]
        self._ids_ref, self._index = None, None

    def _sample_index(self, ids) -> list:
        """Per sharing group, the rows of its samples in ``ids`` (None: the
        same samples in the same order). Kept while ``ids`` is unchanged."""
        if ids is self._ids_ref:
            return self._index
        ids_s = np.asarray(ids).astype(str)
        where = None
        index = []
        for want in self.ids:
            if want.size == ids_s.size and np.array_equal(want, ids_s):
                index.append(None)
                continue
            if where is None:
                where = {s: i for i, s in enumerate(ids_s)}
            missing = [s for s in want if s not in where]
            if missing:
                raise KeyError(f"{len(missing)} step-1 samples are absent from the genotypes "
                               f"(e.g. {missing[:3]})")
            index.append(np.fromiter((where[s] for s in want), dtype=np.int64, count=want.size))
        self._ids_ref, self._index = ids, index
        return index

    def _groups(self, h: int, chromosome) -> np.ndarray:
        chroms, inverse = np.unique(np.asarray(chromosome), return_inverse=True)
        keys = [_chrom_key(c) for c in chroms]
        unknown = [str(c) for c, k in zip(chroms, keys) if k not in self.group_of[h]]
        if unknown:
            raise ValueError(f"variants on chromosomes {unknown} have no LOCO group in step 1; "
                             "test only the chromosomes step 1 was fitted on")
        return np.array([self.group_of[h][k] for k in keys], dtype=np.int64)[inverse.ravel()]

    def blocks(self, source):
        if isinstance(source, Dosages):
            yield from self._dosage_blocks(source)
        elif isinstance(source, Genotypes):
            yield from self._call_blocks(source)
        else:
            raise TypeError("genotypes must be Genotypes, Dosages or an iterable of them")

    def _call_blocks(self, gt):
        from mixmogam._loco import LocoGenotypes
        from mixmogam._packed import PackedCalls
        from mixmogam._stream import Readahead
        from mixmogam.genotypes import _copy_samples, _variant_major
        from mixmogam.twostep import _hratt_test

        index = self._sample_index(gt.sample_ids)
        groups_all = [self._groups(h, gt.chromosome) for h in range(len(self.shared))]
        packed = gt.packed
        row_bytes = (gt.n_samples + 3) // 4 if packed else gt.n_samples
        step = max(1, int(self.chunk_bytes // max(row_bytes, 1)))
        advise = Readahead.for_array(gt.G.data, step * row_bytes) if packed else None
        if advise is not None:
            advise.begin()
        for a in range(0, gt.n_variants, step):
            b = min(a + step, gt.n_variants)
            if advise is not None:
                advise.at(np.arange(a, b))
            block = gt.G.take_variants(slice(a, b)) if packed else gt.G[:, a:b]
            out = [None] * len(self.fits)
            for h, members in enumerate(self.shared):
                calls, rows = block, index[h]
                if rows is not None:
                    if isinstance(calls, PackedCalls):
                        calls = calls.take_samples(rows, n_threads=self.n_threads, in_memory=True)
                    else:
                        order = "F" if _variant_major(calls) else "C"
                        calls = _copy_samples(calls, rows, np.empty((rows.size, b - a),
                                                                    dtype=np.int8, order=order))
                head = self.fits[members[0]]
                sub = Genotypes._from_validated(
                    calls, sample_ids=head._panel["sample_ids"],
                    chromosome=gt.chromosome[a:b], position=gt.position[a:b],
                    variant_ids=gt.variant_ids[a:b],
                    allele1=None if gt.allele1 is None else gt.allele1[a:b],
                    allele2=None if gt.allele2 is None else gt.allele2[a:b])
                lg = LocoGenotypes(sub, groups_all[h][a:b], head._panel["Q"],
                                   block=head._panel["block"], dtype=head._panel["dtype"],
                                   n_threads=self.n_threads, row_scale=head._setup.s,
                                   n_groups=head._panel["n_groups"])
                for k in members:
                    f = self.fits[k]
                    out[k] = _hratt_test(f._state, replace(f._setup, lg=lg), sub, self.n_threads,
                                         cache=self.caches[k])
                del lg, sub, calls
            yield out

    def _dosage_blocks(self, dos):
        for f in self.fits:
            st = f._state
            if st["binary"] or st["weighted"] or st["correction"] == "gls":
                raise ValueError("dosages are tested for unweighted quantitative traits with the "
                                 "lambda or kappa_g correction; use hard calls for binary, "
                                 "weighted or GLS tests")
        index = self._sample_index(dos.sample_ids)
        groups_all = [self._groups(h, dos.chromosome) for h in range(len(self.shared))]
        n = max(ids.size for ids in self.ids)
        step = max(1, int(self.chunk_bytes // max(8 * n, 1)))
        for a in range(0, dos.n_variants, step):
            b = min(a + step, dos.n_variants)
            meta = SimpleNamespace(chromosome=dos.chromosome[a:b], position=dos.position[a:b],
                                   variant_ids=dos.variant_ids[a:b],
                                   allele1=None if dos.allele1 is None else dos.allele1[a:b],
                                   allele2=None if dos.allele2 is None else dos.allele2[a:b])
            out = [None] * len(self.fits)
            for h, members in enumerate(self.shared):
                rows = index[h]
                if rows is None:
                    D = np.array(dos.values[:, a:b], dtype=np.float64)
                else:
                    # Row selection already owns its buffer.
                    D = dos.values[rows, a:b].astype(np.float64, copy=False)
                mean, sd = _standardize_dosages(D)
                groups = groups_all[h][a:b]
                zz = _dosage_norms(D, self.fits[members[0]]._panel["Q"])
                bounds = np.r_[0, np.flatnonzero(groups[1:] != groups[:-1]) + 1, groups.size]
                for k in members:
                    out[k] = _dosage_test(self.fits[k], D, mean, sd, zz, groups, bounds, meta,
                                          self.caches[k])
                del D
            yield out


def _standardize_dosages(D: np.ndarray):
    """Standardize the (n, k) block in place as calls are: mean and SD over
    the observed values (population SD), missing values at the mean (0)."""
    missing = np.isnan(D)
    count = np.maximum((~missing).sum(axis=0), 1)
    D[missing] = 0.0
    mean = D.sum(axis=0) / count
    D -= mean
    D[missing] = 0.0
    sd = np.sqrt(np.einsum("ij,ij->j", D, D) / count)
    D /= np.where(sd > 0, sd, 1.0)
    return mean, sd


def _dosage_norms(Z: np.ndarray, Q: np.ndarray) -> np.ndarray:
    C = Z.T @ Q
    return np.einsum("ij,ij->j", Z, Z) - np.einsum("ij,ij->i", C, C)


def _dosage_scores(Z: np.ndarray, Wp: np.ndarray, groups, bounds) -> np.ndarray:
    """Dot each dosage column with its LOCO residual."""
    num = np.empty(groups.size)
    if bounds.size - 1 <= 2 * Wp.shape[1]:
        for lo, hi in zip(bounds[:-1], bounds[1:]):
            w = np.ascontiguousarray(Wp[:, groups[lo]])
            np.einsum("ij,i->j", Z[:, lo:hi], w, out=num[lo:hi])
    else:
        # Interleaved chromosomes can leave hundreds of tiny runs. Gather
        # their residuals in tiles of at most 8 MiB, or one column.
        width = max(1, (8 * 1024**2) // Wp[:, 0].nbytes)
        for lo in range(0, groups.size, width):
            hi = min(lo + width, groups.size)
            np.einsum("ij,ij->j", Z[:, lo:hi], Wp[:, groups[lo:hi]], out=num[lo:hi])
    return num


def _dosage_test(fit: HrattStep1, Z: np.ndarray, mean, sd, zz, groups, bounds, meta, cache) -> GwasResult:
    """Step 2 of an unweighted quantitative trait for standardized dosage
    rows ``Z`` (n, k): the statistics of :func:`mixmogam.twostep._hratt_test`
    with the values in place of decoded calls."""
    from mixmogam.twostep import _result

    state, setup = fit._state, fit._setup
    Q = fit._panel["Q"]
    if "Wp" not in cache:
        W = state["resid"]
        cache["Wp"] = W - Q @ (Q.T @ W) if Q.shape[1] else W - W.mean(axis=0)
    if "wn" not in cache:
        cache["wn"] = np.einsum("ij,ij->j", state["resid"], state["resid"])
    Wp, wn = cache["Wp"], cache["wn"]
    num = _dosage_scores(Z, Wp, groups, bounds)
    with np.errstate(divide="ignore", invalid="ignore"):
        U = setup.n_eff * num * num / (zz * wn[groups])
    U[~(zz > 0)] = np.nan
    lam, sy, kap = state["lam"], state["sy"], state["kap_g"]
    with np.errstate(divide="ignore", invalid="ignore"):
        if kap is not None:
            chi2 = U / (kap[groups] * state["c_X"])
            beta_z = num / zz * sy
        else:
            chi2 = lam * U
            beta_z = lam * num / zz * sy
        se_z = np.abs(beta_z) / np.sqrt(chi2)
    extra = dict(state["extra"])
    extra["dosage"] = True
    if state["rel"] is not None:  # carriers are not defined for dosages: the columns stay
        k = Z.shape[1]
        extra.update({"carrier_families": np.full(k, np.nan),
                      "few_carrier_families": np.zeros(k, dtype=bool),
                      "carrier_families_p": float(state["rel_opt"]["carrier_families_p"]),
                      "n_few_carrier_families": 0})
    lg = SimpleNamespace(sd=sd, mean=mean, n=Z.shape[0], n_groups=fit._panel["n_groups"])
    return _result(SimpleNamespace(lg=lg, labels=fit._panel["labels"]), meta, chi2, beta_z,
                   se_z, extra)


def _labels(fit) -> list:
    return [tuple(_chrom_key(c) for c in label) for label in fit._panel["labels"]]


def _chrom_key(c) -> str:
    return str(c).removeprefix("chr").removeprefix("CHR")


def _same(a, b) -> bool:
    if a is None or b is None:
        return a is None and b is None
    return np.array_equal(a, b)


# ----------------------------------------------------------------------
# Step-1 files: arrays in the npz, everything else as JSON.
def _encode(obj, arrays: dict, path: str):
    from mixmogam._relatedness import Relatedness

    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if isinstance(obj, (int, float)):
        return obj
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        key = f"a{len(arrays)}"
        if obj.dtype == object:
            if not all(isinstance(v, str) for v in obj.ravel()):
                raise TypeError(f"{path}: object arrays must hold strings")
            arrays[key] = obj.astype(str)
            return {"__array__": key, "object": True}
        arrays[key] = obj
        return {"__array__": key}
    if isinstance(obj, Relatedness):
        u = obj.upper.tocoo()
        return {"__relatedness__": {"n": obj.n, "source": obj.source,
                                    "i": _encode(u.row, arrays, path),
                                    "j": _encode(u.col, arrays, path),
                                    "r": _encode(u.data, arrays, path)}}
    if isinstance(obj, tuple):
        return {"__tuple__": [_encode(v, arrays, f"{path}[{i}]") for i, v in enumerate(obj)]}
    if isinstance(obj, list):
        return [_encode(v, arrays, f"{path}[{i}]") for i, v in enumerate(obj)]
    if isinstance(obj, dict):
        if all(isinstance(k, str) and not k.startswith("__") for k in obj):
            return {k: _encode(v, arrays, f"{path}.{k}") for k, v in obj.items()}
        return {"__items__": [[_encode(k, arrays, path), _encode(v, arrays, f"{path}[{k!r}]")]
                              for k, v in obj.items()]}
    raise TypeError(f"{path}: cannot store a {type(obj).__name__} in a step-1 file")


def _decode(obj, arrays: dict):
    from mixmogam._relatedness import Relatedness

    if isinstance(obj, list):
        return [_decode(v, arrays) for v in obj]
    if not isinstance(obj, dict):
        return obj
    if "__array__" in obj:
        a = arrays[obj["__array__"]]
        return a.astype(object) if obj.get("object") else a
    if "__tuple__" in obj:
        return tuple(_decode(v, arrays) for v in obj["__tuple__"])
    if "__items__" in obj:
        return {_hashable(_decode(k, arrays)): _decode(v, arrays) for k, v in obj["__items__"]}
    if "__relatedness__" in obj:
        r = obj["__relatedness__"]
        return Relatedness(r["n"], _decode(r["i"], arrays), _decode(r["j"], arrays),
                           _decode(r["r"], arrays), source=r["source"])
    return {k: _decode(v, arrays) for k, v in obj.items()}


def _hashable(key):
    return tuple(key) if isinstance(key, list) else key
