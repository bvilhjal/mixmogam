"""Stage-1 model-SNP subsets for HRATT (BOLT-LMM ``--modelSnps``,
LDAK-KVIK's ``--thin-common`` output, REGENIE step 1).

A subset is either supplied (variant indices, a boolean mask or variant
IDs) or chosen by :func:`thin_variants`: variants with MAF above
``min_maf``, then greedy pairwise LD pruning in genome order, keeping a
variant unless its squared correlation with an already kept variant within
``window_bp`` on the same chromosome reaches ``max_r2``. Correlations are
estimated on ``ld_samples`` samples drawn at random (seed ``ld_seed``,
default 0) so that a file sorted by cohort, family or ancestry does not
bias them; ``ld_seed=None`` takes the first ``ld_samples`` samples
instead (one contiguous block of each packed row: a small fraction of a
memory-mapped BED file is read, where random samples touch every page of
each row); ``ld_samples=None`` uses every sample.
"""

from __future__ import annotations

from typing import Optional

import re
import warnings

import numpy as np

from mixmogam._fast import HAS_NUMBA
from mixmogam._packed import CODE_CALLS, PackedCalls

if HAS_NUMBA:
    from numba import njit

__all__ = ["THIN_DEFAULTS", "thin_variants", "resolve_model_snps", "resolve_include"]

THIN_DEFAULTS = {"min_maf": 0.01, "max_r2": 0.5, "window_bp": 100_000,
                 "ld_samples": 10_000, "ld_seed": 0, "max_snps": None, "chunk": 256}


def _greedy(pos, ok, prev_pos, prev_r2, cur_r2, max_r2, window):
    """Sequential keep decisions for one chunk (positions sorted).

    prev_r2 (k, p): r2 of the chunk's variants with the p variants kept
    before the chunk (positions prev_pos); cur_r2 (k, k) within the chunk."""
    k = pos.size
    keep = np.zeros(k, dtype=np.bool_)
    for i in range(k):
        if not ok[i]:
            continue
        bad = False
        for j in range(prev_pos.size):
            if pos[i] - prev_pos[j] <= window and prev_r2[i, j] >= max_r2:
                bad = True
                break
        if not bad:
            for j in range(i):
                if keep[j] and pos[i] - pos[j] <= window and cur_r2[i, j] >= max_r2:
                    bad = True
                    break
        keep[i] = not bad
    return keep


if HAS_NUMBA:
    _greedy = njit(cache=True)(_greedy)


def _ld_rows_packed(data, idx, s, out):
    """Compiled two-bit route of :func:`_ld_rows`: decode, centre and scale
    the first ``s`` samples of each row in one sweep over its bytes."""
    calls = np.empty(s, dtype=np.int8)
    for r in range(idx.size):
        row = data[idx[r]]
        total = 0
        cnt = 0
        for i in range(s):
            g = CODE_CALLS[(row[i >> 2] >> (2 * (i & 3))) & 3]
            calls[i] = g
            if g >= 0:
                total += g
                cnt += 1
        mean = total / max(cnt, 1)
        ss = 0.0
        for i in range(s):
            if calls[i] >= 0:
                d = calls[i] - mean
                ss += d * d
        sd = np.sqrt(ss / max(cnt, 1))
        inv = 1.0 / sd if sd > 0 else 1.0
        for i in range(s):
            out[r, i] = np.float32((calls[i] - mean) * inv) if calls[i] >= 0 else np.float32(0.0)


def _ld_rows_packed_at(data, idx, samples, out):
    """:func:`_ld_rows_packed` on the samples ``samples`` (sorted)."""
    s = samples.size
    calls = np.empty(s, dtype=np.int8)
    for r in range(idx.size):
        row = data[idx[r]]
        total = 0
        cnt = 0
        for t in range(s):
            i = samples[t]
            g = CODE_CALLS[(row[i >> 2] >> (2 * (i & 3))) & 3]
            calls[t] = g
            if g >= 0:
                total += g
                cnt += 1
        mean = total / max(cnt, 1)
        ss = 0.0
        for t in range(s):
            if calls[t] >= 0:
                d = calls[t] - mean
                ss += d * d
        sd = np.sqrt(ss / max(cnt, 1))
        inv = 1.0 / sd if sd > 0 else 1.0
        for t in range(s):
            out[r, t] = np.float32((calls[t] - mean) * inv) if calls[t] >= 0 else np.float32(0.0)


if HAS_NUMBA:
    _ld_rows_packed = njit(cache=True)(_ld_rows_packed)
    _ld_rows_packed_at = njit(cache=True)(_ld_rows_packed_at)


def _ld_sample_set(n, ld_samples, ld_seed):
    """Sorted sample indices for r2, or None for the first ``s`` samples."""
    if ld_samples is None or ld_samples >= n or ld_seed is None:
        return None
    rng = np.random.default_rng(int(ld_seed))
    return np.sort(rng.choice(n, int(ld_samples), replace=False)).astype(np.int64)


def _ld_rows(gt, idx, ld_samples, samples=None):
    """Standardized float32 calls (k, s) of variants ``idx`` on the first
    ``s`` samples, or on ``samples`` (sorted indices) when given; no-calls
    at the mean (0)."""
    n = gt.n_samples
    if samples is not None:
        s = samples.size
        if HAS_NUMBA and isinstance(gt.G, PackedCalls):
            z = np.empty((idx.size, s), dtype=np.float32)
            _ld_rows_packed_at(gt.G.data, np.ascontiguousarray(idx, dtype=np.int64), samples, z)
            return z, s
        if isinstance(gt.G, PackedCalls):
            raw = np.asarray(gt.G.data[idx])
            byte, shift = samples >> 2, (2 * (samples & 3)).astype(np.uint8)
            calls = CODE_CALLS[(raw[:, byte] >> shift) & 3]
        else:
            calls = np.asarray(gt.G[np.ix_(samples, idx)]).T
        return _standardize_rows(calls), s
    s = n if ld_samples is None else min(int(ld_samples), n)
    if HAS_NUMBA and isinstance(gt.G, PackedCalls):
        z = np.empty((idx.size, s), dtype=np.float32)
        _ld_rows_packed(gt.G.data, np.ascontiguousarray(idx, dtype=np.int64), s, z)
        return z, s
    if isinstance(gt.G, PackedCalls):
        nb = (s + 3) // 4
        raw = np.asarray(gt.G.data[idx, :nb])  # only the first nb bytes of each row
        codes = (raw[:, :, None] >> np.array([0, 2, 4, 6], dtype=np.uint8)) & 3
        calls = CODE_CALLS[codes.reshape(idx.size, -1)[:, :s]]
    else:
        calls = np.asarray(gt.G[:s, idx]).T
    return _standardize_rows(calls), s


def _standardize_rows(calls):
    z = calls.astype(np.float32)
    called = calls >= 0
    cnt = np.maximum(called.sum(axis=1), 1)
    z[~called] = 0.0
    mean = z.sum(axis=1) / cnt
    z -= mean[:, None]
    z[~called] = 0.0
    sd = np.sqrt((z * z).sum(axis=1) / cnt)
    z /= np.where(sd > 0, sd, 1.0)[:, None]
    return z


def thin_variants(gt, maf: Optional[np.ndarray] = None, *, min_maf: float = 0.01,
                  max_r2: float = 0.5, window_bp: int = 100_000,
                  ld_samples: Optional[int] = 10_000, ld_seed: Optional[int] = 0,
                  max_snps: Optional[int] = None, chunk: int = 256) -> np.ndarray:
    """Sorted indices of an LD-thinned set of common variants.

    ``maf`` (m,) is the minor allele frequency over all samples (computed
    when None, one pass); variants at or below ``min_maf`` are never kept.
    Within each chromosome, in position order, a variant is kept unless its
    r2 (on ``ld_samples`` random samples, seed ``ld_seed``; the first
    ``ld_samples`` when ``ld_seed`` is None) with a kept variant at most
    ``window_bp`` before it reaches ``max_r2``: PLINK ``--indep-pairwise``
    with a base-pair window, pruning toward earlier variants. ``max_snps``
    caps the result by keeping an evenly strided subset of the thinned set.
    """
    if not (0.0 <= min_maf < 0.5):
        raise ValueError("min_maf must lie in [0, 0.5)")
    if not (0.0 < max_r2 <= 1.0):
        raise ValueError("max_r2 must lie in (0, 1]")
    if window_bp < 0:
        raise ValueError("window_bp must be non-negative")
    if ld_samples is not None and ld_samples < 2:
        raise ValueError("ld_samples must be at least 2 (or None for all samples)")
    if ld_seed is not None and (not isinstance(ld_seed, (int, np.integer)) or isinstance(ld_seed, bool)
                                or ld_seed < 0):
        raise ValueError("ld_seed must be a nonnegative integer or None (the first ld_samples samples)")
    if max_snps is not None and max_snps < 1:
        raise ValueError("max_snps must be a positive integer or None")
    if maf is None:
        maf = gt.allele_freqs(minor=True)
    maf = np.asarray(maf, dtype=np.float64)
    if maf.shape != (gt.n_variants,):
        raise ValueError("maf must hold one frequency per variant")
    chrom = np.asarray(gt.chromosome)
    pos = np.asarray(gt.position, dtype=np.int64)
    samples = _ld_sample_set(gt.n_samples, ld_samples, ld_seed)
    kept_all = []
    for c in np.unique(chrom):
        members = np.flatnonzero(chrom == c)
        members = members[np.argsort(pos[members], kind="stable")]
        prev_idx = np.empty(0, dtype=np.int64)
        prev_pos = np.empty(0, dtype=np.int64)
        prev_z = None
        for start in range(0, members.size, chunk):
            idx = members[start : start + chunk]
            p = pos[idx]
            ok = maf[idx] > min_maf
            if not ok.any():
                continue
            z, s = _ld_rows(gt, idx, ld_samples, samples)
            # kept variants still inside the window of this chunk's first variant
            live = prev_pos >= p[0] - window_bp
            prev_idx, prev_pos = prev_idx[live], prev_pos[live]
            prev_z = None if prev_z is None else prev_z[live]
            if prev_z is not None and prev_z.shape[0]:
                prev_r2 = ((z @ prev_z.T) / s) ** 2
            else:
                prev_r2 = np.zeros((idx.size, 0), dtype=np.float32)
            cur_r2 = ((z @ z.T) / s) ** 2
            keep = _greedy(p, ok, prev_pos, np.ascontiguousarray(prev_r2),
                           np.ascontiguousarray(cur_r2), float(max_r2), int(window_bp))
            kept = idx[keep]
            kept_all.append(kept)
            prev_idx = np.r_[prev_idx, kept]
            prev_pos = np.r_[prev_pos, p[keep]]
            prev_z = z[keep] if prev_z is None else np.vstack([prev_z, z[keep]])
    out = np.sort(np.concatenate(kept_all)) if kept_all else np.empty(0, dtype=np.int64)
    if max_snps is not None and out.size > max_snps:
        out = out[np.unique(np.linspace(0, out.size - 1, int(max_snps)).round().astype(np.int64))]
    return out.astype(np.int64)


def resolve_model_snps(gt, model_snps, maf=None, options: Optional[dict] = None) -> tuple:
    """(indices, source, settings) of the stage-1 variants.

    ``model_snps``: ``"auto"`` (:func:`thin_variants` with ``options``),
    ``"hm3"`` (the panel's HapMap3 variants, :func:`mixmogam._hm3.match_hm3`;
    options ``build``, ``map``, ``min_frac``, ``min_snps``, ``min_maf``,
    ``drop_ambiguous``, plus thinning options for the fallback), an integer
    index array, a boolean mask of length m, or variant IDs (str).

    With ``"hm3"``, fewer than ``min_frac`` (0.5) of the map's variants on
    the panel's chromosomes or
    fewer than ``min_snps`` (300,000) matched means the panel is not a
    HapMap3-covering human panel (an unimputed array, another build,
    another species): HRATT warns and falls back to ``"auto"``; the source
    is then ``"auto (hm3 fallback)"`` and ``settings["hm3"]`` keeps the
    match report.
    """
    options = dict(options or {})
    if isinstance(model_snps, str):
        if model_snps == "auto":
            unknown = set(options) - set(THIN_DEFAULTS)
            if unknown:
                raise ValueError(f"unknown model_snps_options: {sorted(unknown)}")
            settings = {**THIN_DEFAULTS, **options}
            idx = thin_variants(gt, maf, **settings)
            return idx, "auto", settings
        if model_snps == "hm3":
            from mixmogam._hm3 import HM3_OPTIONS, format_report, match_hm3

            unknown = set(options) - set(THIN_DEFAULTS) - set(HM3_OPTIONS)
            if unknown:
                raise ValueError(f"unknown model_snps_options: {sorted(unknown)}")
            hopt = {**HM3_OPTIONS, **{k: v for k, v in options.items() if k in HM3_OPTIONS}}
            idx, rep = match_hm3(gt, build=hopt["build"], map=hopt["map"], maf=maf,
                                 min_maf=hopt["min_maf"], drop_ambiguous=hopt["drop_ambiguous"])
            text = format_report(rep)
            settings = {"hm3": rep, "hm3_options": hopt, "log": text}
            low = (idx.size < hopt["min_frac"] * rep["n_hm3_on_panel_chromosomes"]
                   or idx.size < hopt["min_snps"])
            if not low:
                return idx, "hm3", settings
            thin_opt = {k: v for k, v in options.items() if k in THIN_DEFAULTS}
            thin_settings = {**THIN_DEFAULTS, **thin_opt}
            warnings.warn(f"model_snps='hm3' matched only {idx.size:,} of the "
                          f"{rep['n_hm3_on_panel_chromosomes']:,} HapMap3 variants on the panel's "
                          f"chromosomes (needs {100 * hopt['min_frac']:.0f}% and "
                          f"{hopt['min_snps']:,}); "
                          "falling back to model_snps='auto' (MAF/LD thinning). Is this an imputed "
                          "human panel, and is the build right?\n" + text,
                          RuntimeWarning, stacklevel=3)
            idx = thin_variants(gt, maf, **thin_settings)
            settings.update(thin_settings)
            return idx, "auto (hm3 fallback)", settings
        raise ValueError("model_snps must be 'auto', 'hm3', variant indices, a boolean mask or variant IDs")
    if options:
        raise ValueError("model_snps_options apply to model_snps='auto' or 'hm3' only")
    idx = _resolve_indices(gt, model_snps, "model_snps", strict=True)
    if idx.size < 2:
        raise ValueError("model_snps must select at least two variants")
    return idx, "user", None


def _resolve_indices(gt, spec, name, strict):
    arr = np.asarray(spec)
    m = gt.n_variants
    if arr.dtype == bool:
        if arr.shape != (m,):
            raise ValueError(f"a boolean {name} mask needs one entry per variant")
        return np.flatnonzero(arr)
    if arr.dtype.kind in "iu":
        idx = np.unique(arr.astype(np.int64))
        if idx.size and (idx[0] < 0 or idx[-1] >= m):
            raise ValueError(f"{name} index out of range")
        return idx
    if arr.dtype.kind in "USO":
        ids = np.asarray(gt.variant_ids).astype(str)
        want = np.unique(arr.astype(str))
        hit = np.isin(ids, want)
        missing = np.setdiff1d(want, ids[hit])
        if missing.size and strict:
            raise ValueError(f"{missing.size} {name} IDs are not among the variants, "
                             f"e.g. {list(missing[:3])}")
        return np.flatnonzero(hit)
    raise ValueError(f"{name} must be 'auto', 'hm3', variant indices, a boolean mask or variant IDs")


_RANGE = re.compile(r"^(?:chr)?([0-9A-Za-z]+):(\d+)(?:-(\d+))?$", re.IGNORECASE)


def resolve_include(gt, include) -> tuple:
    """(sorted indices, record) of the variants forced into the model SNPs.

    ``include``: variant indices, a boolean mask, or strings, each a
    variant ID or (when no variant has that ID) a region ``chr:pos`` or
    ``chr:start-end`` (inclusive; 'chr1', '1', 'X' and '23' name the same
    chromosome). IDs and regions that match nothing are counted and
    reported (with a warning), not an error: lists of known loci rarely
    match a panel entirely.
    """
    arr = np.asarray(include)
    if arr.dtype.kind not in "USO":
        idx = _resolve_indices(gt, arr, "model_snps_include", strict=True)
        return idx, {"n": int(idx.size), "unmatched": []}
    from mixmogam._hm3 import normalize_chromosome

    items = [str(x).strip() for x in arr.reshape(-1) if str(x).strip()]
    ids = np.asarray(gt.variant_ids).astype(str)
    is_id = np.isin(np.asarray(items, dtype=str), ids) if items else np.zeros(0, bool)
    sel = np.isin(ids, np.asarray(items, dtype=str)[is_id]) if items else np.zeros(gt.n_variants, bool)
    chrom = normalize_chromosome(gt.chromosome)
    pos = np.asarray(gt.position, dtype=np.int64)
    unmatched = []
    for item, known in zip(items, is_id):
        if known:
            continue
        m = _RANGE.match(item)
        if not m:
            unmatched.append(item)
            continue
        c = int(normalize_chromosome(np.array([m.group(1)]))[0])
        lo = int(m.group(2))
        hi = int(m.group(3)) if m.group(3) else lo
        if hi < lo:
            raise ValueError(f"model_snps_include region {item!r} ends before it starts")
        hit = (chrom == c) & (pos >= lo) & (pos <= hi)
        if not hit.any():
            unmatched.append(item)
        sel |= hit
    if unmatched:
        warnings.warn(f"model_snps_include: {len(unmatched)} of {len(items)} entries match no "
                      f"variant, e.g. {unmatched[:3]}", RuntimeWarning, stacklevel=3)
    idx = np.flatnonzero(sel)
    return idx, {"n": int(idx.size), "unmatched": unmatched}


def read_snp_list(path) -> np.ndarray:
    """Variant IDs, one per line (first whitespace field; '#' comments and blanks skipped)."""
    ids = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                ids.append(line.split()[0])
    if not ids:
        raise ValueError(f"no variant IDs in {path}")
    return np.array(ids)
