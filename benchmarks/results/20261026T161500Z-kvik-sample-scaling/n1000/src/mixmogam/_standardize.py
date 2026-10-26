"""Genotype preparation and decoding over independent variants.

Preparation computes, per variant, the called-sample mean and standard
deviation from exact call counts, the covariate coefficients c = z Q of the
standardized values and the squared norm of the projected values z - c Q'.
The Numba kernels take tiles of 64 variants, from int8 calls or two-bit
codes: one contiguous pass counts each variant's calls, then eight
variants at a time are centered into 256-sample float64 rows (16 KB of
scratch) whose covariate sums stay in registers. Sample-major int8 is
read a sample at a time instead. Each variant sums its samples in order,
so every thread count, storage layout and format gives the same values.
No standardized or projected block is written. With per-sample squared
weights (the row scaling of weighted HRATT), the projected norm is that
of the scaled values, sum_i v_i z_i^2 - |c|^2, accumulated over the same
rows.

Every later pass decodes through a per-variant table: hard calls 0/1/2 can
only standardize to three values, computed once in float64, so decoding
selects them (no-calls are 0) without arithmetic. One M2 Pro thread decodes
about 0.16 ns per genotype from variant-major int8 (PLINK order), 0.2 ns
from two-bit codes and 1 ns from sample-major int8. Covariates are
removed by the consumers, never per variant.
"""

from __future__ import annotations

import math

import numpy as np

from mixmogam._fast import HAS_NUMBA
from mixmogam._packed import PackedCalls
from mixmogam._vb import _numba_thread_limit, _validate_n_threads

if HAS_NUMBA:
    from numba import njit, prange
else:  # the NumPy fallbacks below never invoke these kernels
    prange = range

# Decoded genotypes per NumPy fallback chunk: the table gather needs an
# intp index per element.
_NUMPY_DECODE_CELLS = 1 << 22
# Variants per preparation tile; samples and variants per tile of
# sample-major decoding.
_TILE = 64
# Samples per centered row of preparation (eight rows per pass).
_SPAN = 256
# Two-bit codes 0-3 of each byte value, counted in 16-bit fields: summed
# over at most _CODE_BYTES bytes, no field carries into the next or into
# the sign bit.
_BYTE_CODES = sum(np.int64(1) << (16 * ((np.arange(256, dtype=np.int64) >> shift) & 3))
                  for shift in (0, 2, 4, 6))
_CODE_BYTES = 8191


def _count_moments(n0, n1, n2):
    """Called-sample mean, sum of squared deviations and divisor of a
    variant from its call counts: exact integers, no cancellation."""
    called = n0 + n1 + n2
    divisor = max(called, 1)
    mu = (n1 + 2 * n2) / divisor
    d0, d1, d2 = 0.0 - mu, 1.0 - mu, 2.0 - mu
    return mu, n0 * d0 * d0 + n1 * d1 * d1 + n2 * d2 * d2, divisor


def _count_calls(col, counts, t):
    """Call counts of one variant's int8 calls into counts[t], from integer
    sums that vectorize where per-call comparisons do not:
    sum g = n1 + 2 n2 - n_miss and sum g^2 = n1 + 4 n2 + n_miss."""
    total = 0
    squares = 0
    missing = 0
    for i in range(col.size):
        value = col[i]
        total += value
        squares += value * value
        missing += value < 0
    twos = (squares - total) // 2 - missing
    ones = total + missing - 2 * twos
    counts[t, 0] = col.size - missing - ones - twos
    counts[t, 1] = ones
    counts[t, 2] = twos


def _count_codes(row, n, counts, t):
    """Call counts of one variant's n two-bit codes into counts[t]: whole
    bytes through the per-byte code counts, the last byte's samples one
    at a time (its unused bits are ignored)."""
    full = n // 4
    codes = np.zeros(4, dtype=np.int64)
    for b0 in range(0, full, _CODE_BYTES):
        fields = 0
        for b in range(b0, min(b0 + _CODE_BYTES, full)):
            fields += _BYTE_CODES[row[b]]
        for code in range(4):
            codes[code] += (fields >> (16 * code)) & 0xFFFF
    for i in range(4 * full, n):
        codes[(row[i >> 2] >> (2 * (i & 3))) & 3] += 1
    # code 0 = call 2, 1 = no-call, 2 = call 1, 3 = call 0
    counts[t, 0], counts[t, 1], counts[t, 2] = codes[3], codes[2], codes[0]


def _sums_one(rows, Q, k, i0, i1, projection, t0):
    """projection[k, t0:t0 + 8] += sum_i rows[u, i - i0] Q[i, k] over
    samples i0:i1 in order: eight independent sums held in registers."""
    a0, a1, a2, a3 = projection[k, t0], projection[k, t0 + 1], projection[k, t0 + 2], projection[k, t0 + 3]
    a4, a5, a6, a7 = projection[k, t0 + 4], projection[k, t0 + 5], projection[k, t0 + 6], projection[k, t0 + 7]
    r0, r1, r2, r3, r4, r5, r6, r7 = rows[0], rows[1], rows[2], rows[3], rows[4], rows[5], rows[6], rows[7]
    for i in range(i1 - i0):
        qik = Q[i0 + i, k]
        a0 += r0[i] * qik
        a1 += r1[i] * qik
        a2 += r2[i] * qik
        a3 += r3[i] * qik
        a4 += r4[i] * qik
        a5 += r5[i] * qik
        a6 += r6[i] * qik
        a7 += r7[i] * qik
    projection[k, t0], projection[k, t0 + 1], projection[k, t0 + 2], projection[k, t0 + 3] = a0, a1, a2, a3
    projection[k, t0 + 4], projection[k, t0 + 5], projection[k, t0 + 6], projection[k, t0 + 7] = a4, a5, a6, a7


def _sums_two(rows, Q, k, i0, i1, projection, t0):
    """:func:`_sums_one` for covariates k and k + 1 at once: each loaded
    value serves two sums (sixteen in registers)."""
    a0, a1, a2, a3 = projection[k, t0], projection[k, t0 + 1], projection[k, t0 + 2], projection[k, t0 + 3]
    a4, a5, a6, a7 = projection[k, t0 + 4], projection[k, t0 + 5], projection[k, t0 + 6], projection[k, t0 + 7]
    b0, b1, b2, b3 = projection[k + 1, t0], projection[k + 1, t0 + 1], projection[k + 1, t0 + 2], projection[k + 1, t0 + 3]
    b4, b5, b6, b7 = projection[k + 1, t0 + 4], projection[k + 1, t0 + 5], projection[k + 1, t0 + 6], projection[k + 1, t0 + 7]
    r0, r1, r2, r3, r4, r5, r6, r7 = rows[0], rows[1], rows[2], rows[3], rows[4], rows[5], rows[6], rows[7]
    for i in range(i1 - i0):
        qa = Q[i0 + i, k]
        qb = Q[i0 + i, k + 1]
        c = r0[i]
        a0 += c * qa
        b0 += c * qb
        c = r1[i]
        a1 += c * qa
        b1 += c * qb
        c = r2[i]
        a2 += c * qa
        b2 += c * qb
        c = r3[i]
        a3 += c * qa
        b3 += c * qb
        c = r4[i]
        a4 += c * qa
        b4 += c * qb
        c = r5[i]
        a5 += c * qa
        b5 += c * qb
        c = r6[i]
        a6 += c * qa
        b6 += c * qb
        c = r7[i]
        a7 += c * qa
        b7 += c * qb
    projection[k, t0], projection[k, t0 + 1], projection[k, t0 + 2], projection[k, t0 + 3] = a0, a1, a2, a3
    projection[k, t0 + 4], projection[k, t0 + 5], projection[k, t0 + 6], projection[k, t0 + 7] = a4, a5, a6, a7
    projection[k + 1, t0], projection[k + 1, t0 + 1], projection[k + 1, t0 + 2], projection[k + 1, t0 + 3] = b0, b1, b2, b3
    projection[k + 1, t0 + 4], projection[k + 1, t0 + 5], projection[k + 1, t0 + 6], projection[k + 1, t0 + 7] = b4, b5, b6, b7


def _sums_weighted(rows, sqw, i0, i1, wsq, t0):
    """wsq[t0:t0 + 8] += sum_i sqw_i rows[u, i - i0]^2 over samples i0:i1
    in order."""
    a0, a1, a2, a3 = wsq[t0], wsq[t0 + 1], wsq[t0 + 2], wsq[t0 + 3]
    a4, a5, a6, a7 = wsq[t0 + 4], wsq[t0 + 5], wsq[t0 + 6], wsq[t0 + 7]
    r0, r1, r2, r3, r4, r5, r6, r7 = rows[0], rows[1], rows[2], rows[3], rows[4], rows[5], rows[6], rows[7]
    for i in range(i1 - i0):
        vi = sqw[i0 + i]
        a0 += vi * r0[i] * r0[i]
        a1 += vi * r1[i] * r1[i]
        a2 += vi * r2[i] * r2[i]
        a3 += vi * r3[i] * r3[i]
        a4 += vi * r4[i] * r4[i]
        a5 += vi * r5[i] * r5[i]
        a6 += vi * r6[i] * r6[i]
        a7 += vi * r7[i] * r7[i]
    wsq[t0], wsq[t0 + 1], wsq[t0 + 2], wsq[t0 + 3] = a0, a1, a2, a3
    wsq[t0 + 4], wsq[t0 + 5], wsq[t0 + 6], wsq[t0 + 7] = a4, a5, a6, a7


def _accumulate_rows(rows, Q, i0, i1, projection, sqw, wsq, t0):
    """Add samples i0:i1 of eight centered variants t0:t0 + 8 (``rows``)
    to their covariate sums (and weighted squares)."""
    q = Q.shape[1]
    for k in range(0, q - 1, 2):
        _sums_two(rows, Q, k, i0, i1, projection, t0)
    if q % 2:
        _sums_one(rows, Q, q - 1, i0, i1, projection, t0)
    if sqw.size > 0:
        _sums_weighted(rows, sqw, i0, i1, wsq, t0)


def _finish_tile(counts, mu, mean, sd, coefficients, zz, r0, width, projection, q,
                 wsq, weighted):
    """SD, covariate coefficients and projected squared norm of a tile;
    ``weighted`` takes the norm of the row-scaled values from ``wsq``."""
    for t in range(width):
        j = r0 + t
        _, sumsq, divisor = _count_moments(counts[t, 0], counts[t, 1], counts[t, 2])
        mean[j] = mu[t]
        sigma = math.sqrt(sumsq / divisor)
        sd[j] = sigma
        scale = sigma if sigma > 0 else 1.0
        norm = (wsq[t] if weighted else sumsq) / (scale * scale)
        for k in range(q):
            c = projection[k, t] / scale
            coefficients[j, k] = c
            norm -= c * c
        zz[j] = norm


def _moments_tile(G, idx, Q, mean, sd, coefficients, zz, r0, r1, projection, rows,
                  sqw, wsq):
    """Moments of variants r0:r1 of int8 calls.

    The mean and variance follow from each variant's call counts; the
    covariate coefficients of the standardized values are
    c = (sum_i (g_i - mu) Q_i) / sd and the projected squared norm is
    sum z^2 - |c|^2, so no standardized vector is formed. Eight variants
    at a time are centered into ``rows`` (no-calls are 0), a span of
    samples per pass, and every sum runs over the samples in order, so a
    variant's values do not depend on the tile, storage layout or format.
    A non-empty ``sqw`` (one squared row scale per sample) also accumulates
    sum_i sqw_i (g_i - mu)^2. ``projection`` holds the tile's (q, 64)
    sums; padding rows of a partial group are zero and their sums unused.
    """
    n, q = Q.shape
    width = r1 - r0
    weighted = sqw.size > 0
    counts = np.empty((width, 3), dtype=np.int64)
    for t in range(width):
        _count_calls(G[:, idx[r0 + t]], counts, t)
    mu = np.empty(width)
    for t in range(width):
        mu[t] = _count_moments(counts[t, 0], counts[t, 1], counts[t, 2])[0]
    projection[:] = 0.0
    wsq[:] = 0.0
    if q or weighted:
        for t0 in range(0, width, 8):
            for i0 in range(0, n, _SPAN):
                i1 = min(i0 + _SPAN, n)
                for u in range(8):
                    out = rows[u]
                    if t0 + u < width:
                        col = G[i0:i1, idx[r0 + t0 + u]]
                        m = mu[t0 + u]
                        for i in range(i1 - i0):
                            value = col[i]
                            out[i] = value - m if value >= 0 else 0.0
                    else:
                        out[:i1 - i0] = 0.0
                _accumulate_rows(rows, Q, i0, i1, projection, sqw, wsq, t0)
    _finish_tile(counts, mu, mean, sd, coefficients, zz, r0, width, projection, q,
                 wsq, weighted)


def _moments_samples_tile(G, idx, Q, mean, sd, coefficients, zz, r0, r1, projection,
                          centered, sqw, wsq):
    """:func:`_moments_tile` for sample-major int8 calls, a sample at a
    time: column reads would stride across the samples' rows, so each
    sample's calls are read along its row. The same counts and the same
    per-sample arithmetic, so identical values."""
    n, q = Q.shape
    width = r1 - r0
    weighted = sqw.size > 0
    counts = np.zeros((width, 3), dtype=np.int64)
    for i in range(n):
        for t in range(width):
            value = G[i, idx[r0 + t]]
            if value >= 0:
                counts[t, value] += 1
    mu = np.empty(width)
    for t in range(width):
        mu[t] = _count_moments(counts[t, 0], counts[t, 1], counts[t, 2])[0]
    projection[:, :width] = 0.0
    wsq[:width] = 0.0
    for i in range(n if q or weighted else 0):  # nothing to sum without either
        for t in range(width):
            value = G[i, idx[r0 + t]]
            centered[t] = value - mu[t] if value >= 0 else 0.0
        if weighted:
            vi = sqw[i]
            for t in range(width):
                wsq[t] += vi * centered[t] * centered[t]
        for k in range(q):
            qik = Q[i, k]
            for t in range(width):
                projection[k, t] += centered[t] * qik
    _finish_tile(counts, mu, mean, sd, coefficients, zz, r0, width, projection, q,
                 wsq, weighted)


def _moments_packed_tile(data, n, idx, Q, mean, sd, coefficients, zz, r0, r1,
                         projection, rows, deviation, sqw, wsq):
    """:func:`_moments_tile` for two-bit PLINK codes (variant-major rows):
    the same counts and the same per-sample arithmetic, so identical values."""
    q = Q.shape[1]
    width = r1 - r0
    weighted = sqw.size > 0
    counts = np.empty((width, 3), dtype=np.int64)
    for t in range(width):
        _count_codes(data[idx[r0 + t]], n, counts, t)
    mu = np.empty(width)
    for t in range(width):
        mu[t] = _count_moments(counts[t, 0], counts[t, 1], counts[t, 2])[0]
        # code 0 = call 2, 1 = no-call, 2 = call 1, 3 = call 0
        deviation[t, 0] = 2 - mu[t]
        deviation[t, 1] = 0.0
        deviation[t, 2] = 1 - mu[t]
        deviation[t, 3] = 0 - mu[t]
    projection[:] = 0.0
    wsq[:] = 0.0
    if q or weighted:
        for t0 in range(0, width, 8):
            for i0 in range(0, n, _SPAN):  # _SPAN is whole bytes: i0 starts a byte
                i1 = min(i0 + _SPAN, n)
                for u in range(8):
                    out = rows[u]
                    t = t0 + u
                    if t < width:
                        row = data[idx[r0 + t]]
                        d0, d1, d2, d3 = deviation[t, 0], deviation[t, 1], deviation[t, 2], deviation[t, 3]
                        full = (i1 - i0) // 4
                        for b in range(full):
                            byte = row[(i0 >> 2) + b]
                            for s in range(4):
                                code = (byte >> (2 * s)) & 3
                                out[4 * b + s] = (d3 if code == 3 else d2) if code >= 2 else (d1 if code == 1 else d0)
                        for i in range(i0 + 4 * full, i1):
                            out[i - i0] = deviation[t, (row[i >> 2] >> (2 * (i & 3))) & 3]
                    else:
                        out[:i1 - i0] = 0.0
                _accumulate_rows(rows, Q, i0, i1, projection, sqw, wsq, t0)
    _finish_tile(counts, mu, mean, sd, coefficients, zz, r0, width, projection, q,
                 wsq, weighted)


def _moments_any_tile(G, idx, Q, mean, sd, coefficients, zz, r0, r1, projection, rows,
                      centered, sqw, wsq):
    """The int8 tile kernel for the storage order (identical values)."""
    if G.strides[0] > G.strides[1]:  # sample-major
        _moments_samples_tile(G, idx, Q, mean, sd, coefficients, zz, r0, r1, projection,
                              centered, sqw, wsq)
    else:
        _moments_tile(G, idx, Q, mean, sd, coefficients, zz, r0, r1, projection, rows,
                      sqw, wsq)


def _moments_columns_serial(G, idx, Q, mean, sd, coefficients, zz, sqw):
    projection = np.empty((Q.shape[1], _TILE))
    rows = np.empty((8, _SPAN))
    centered = np.empty(_TILE)
    wsq = np.empty(_TILE)
    for r0 in range(0, idx.size, _TILE):
        _moments_any_tile(G, idx, Q, mean, sd, coefficients, zz, r0,
                          min(r0 + _TILE, idx.size), projection, rows, centered, sqw, wsq)


# Separate from the serial kernel: Numba caches per function, and a parallel
# kernel would start the threading layer even for one worker.
def _moments_columns_parallel(G, idx, Q, mean, sd, coefficients, zz, sqw):
    for t in prange((idx.size + _TILE - 1) // _TILE):
        projection = np.empty((Q.shape[1], _TILE))
        rows = np.empty((8, _SPAN))
        centered = np.empty(_TILE)
        wsq = np.empty(_TILE)
        _moments_any_tile(G, idx, Q, mean, sd, coefficients, zz, t * _TILE,
                          min(t * _TILE + _TILE, idx.size), projection, rows, centered, sqw, wsq)


def _moments_packed_serial(data, n, idx, Q, mean, sd, coefficients, zz, sqw):
    projection = np.empty((Q.shape[1], _TILE))
    rows = np.empty((8, _SPAN))
    deviation = np.empty((_TILE, 4))
    wsq = np.empty(_TILE)
    for r0 in range(0, idx.size, _TILE):
        _moments_packed_tile(data, n, idx, Q, mean, sd, coefficients, zz, r0,
                             min(r0 + _TILE, idx.size), projection, rows, deviation,
                             sqw, wsq)


def _moments_packed_parallel(data, n, idx, Q, mean, sd, coefficients, zz, sqw):
    for t in prange((idx.size + _TILE - 1) // _TILE):
        projection = np.empty((Q.shape[1], _TILE))
        rows = np.empty((8, _SPAN))
        deviation = np.empty((_TILE, 4))
        wsq = np.empty(_TILE)
        _moments_packed_tile(data, n, idx, Q, mean, sd, coefficients, zz, t * _TILE,
                             min(t * _TILE + _TILE, idx.size), projection, rows, deviation,
                             sqw, wsq)


def _select_call(v, l0, l1, l2, missing):
    """The value of validated call ``v`` (-1, 0, 1 or 2): the lookup row's
    entries chosen by comparisons, which vectorize where a table gather
    does not (0.16 against 0.26 ns per genotype on one M2 Pro thread)."""
    return (l2 if v == 2 else l1) if v >= 1 else (l0 if v == 0 else missing)


def _decode_column(col, l0, l1, l2, missing, out):
    for i in range(col.size):
        out[i] = _select_call(col[i], l0, l1, l2, missing)


def _decode_variants(G, idx, lut, out, start):
    """Variant-major storage: each variant's calls are contiguous. Rows of
    ``out`` hold samples start:start + out.shape[1]."""
    stop = start + out.shape[1]
    for r in range(idx.size):
        _decode_column(G[start:stop, idx[r]], lut[r, 0], lut[r, 1], lut[r, 2], lut[r, 3], out[r])


def _decode_variants_parallel(G, idx, lut, out, start):
    stop = start + out.shape[1]
    for r in prange(idx.size):
        _decode_column(G[start:stop, idx[r]], lut[r, 0], lut[r, 1], lut[r, 2], lut[r, 3], out[r])


def _decode_tile(G, idx, lut, out, i0, i1, start):
    for r0 in range(0, idx.size, _TILE):
        for r in range(r0, min(r0 + _TILE, idx.size)):
            j = idx[r]
            for i in range(i0, i1):
                out[r, i - start] = lut[r, G[i, j] & 3]


def _decode_samples(G, idx, lut, out, start):
    """Sample-major storage: decode square tiles, so the strided reads of a
    tile's sample rows and its row writes both stay in cache."""
    stop = start + out.shape[1]
    for i0 in range(start, stop, _TILE):
        _decode_tile(G, idx, lut, out, i0, min(i0 + _TILE, stop), start)


def _decode_samples_parallel(G, idx, lut, out, start):
    stop = start + out.shape[1]
    for t in prange((out.shape[1] + _TILE - 1) // _TILE):
        i0 = start + t * _TILE
        _decode_tile(G, idx, lut, out, i0, min(i0 + _TILE, stop), start)


def _decode_packed_rows(data, idx, lut, out, r0, r1, code_lut):
    """Rows r0:r1 from two-bit codes; ``lut`` is in call order (call & 3)."""
    n = out.shape[1]
    full = n // 4
    for r in range(r0, r1):
        code_lut[0], code_lut[1], code_lut[2], code_lut[3] = lut[r, 2], lut[r, 3], lut[r, 1], lut[r, 0]
        row = data[idx[r]]
        for b in range(full):
            byte = row[b]
            i = 4 * b
            out[r, i] = code_lut[byte & 3]
            out[r, i + 1] = code_lut[(byte >> 2) & 3]
            out[r, i + 2] = code_lut[(byte >> 4) & 3]
            out[r, i + 3] = code_lut[byte >> 6]
        for i in range(4 * full, n):
            out[r, i] = code_lut[(row[i >> 2] >> (2 * (i & 3))) & 3]


def _decode_packed_span(data, idx, lut, out, r0, r1, code_lut, start):
    """Rows r0:r1 of samples start:start + out.shape[1] from two-bit codes.
    Whole bytes serve four samples; a range that starts or ends inside a
    byte decodes those samples one at a time. (Full-width decoding uses
    :func:`_decode_packed_rows`, which Numba compiles much tighter.)"""
    stop = start + out.shape[1]
    first = min(stop, (start + 3) // 4 * 4)  # the first sample of a whole byte
    last = max(first, stop // 4 * 4)  # and the end of the whole bytes
    for r in range(r0, r1):
        code_lut[0], code_lut[1], code_lut[2], code_lut[3] = lut[r, 2], lut[r, 3], lut[r, 1], lut[r, 0]
        row = data[idx[r]]
        for i in range(start, first):
            out[r, i - start] = code_lut[(row[i >> 2] >> (2 * (i & 3))) & 3]
        for b in range(first // 4, last // 4):
            byte = row[b]
            i = 4 * b - start
            out[r, i] = code_lut[byte & 3]
            out[r, i + 1] = code_lut[(byte >> 2) & 3]
            out[r, i + 2] = code_lut[(byte >> 4) & 3]
            out[r, i + 3] = code_lut[byte >> 6]
        for i in range(last, stop):
            out[r, i - start] = code_lut[(row[i >> 2] >> (2 * (i & 3))) & 3]


def _decode_packed(data, idx, lut, out):
    """Variant-major two-bit codes: one byte serves four samples."""
    _decode_packed_rows(data, idx, lut, out, 0, idx.size, np.empty(4, dtype=out.dtype))


def _decode_packed_parallel(data, idx, lut, out):
    for t in prange((idx.size + _TILE - 1) // _TILE):
        _decode_packed_rows(data, idx, lut, out, t * _TILE, min(t * _TILE + _TILE, idx.size),
                            np.empty(4, dtype=out.dtype))


def _decode_packed_span_all(data, idx, lut, out, start):
    _decode_packed_span(data, idx, lut, out, 0, idx.size, np.empty(4, dtype=out.dtype), start)


def _decode_packed_span_parallel(data, idx, lut, out, start):
    for t in prange((idx.size + _TILE - 1) // _TILE):
        _decode_packed_span(data, idx, lut, out, t * _TILE, min(t * _TILE + _TILE, idx.size),
                            np.empty(4, dtype=out.dtype), start)


def _gather_variants(G, idx, lut, rows, scale, out):
    """Variant-major storage: samples ``rows`` (in that order) of each
    variant's contiguous calls into the columns of ``out``. A non-empty
    ``scale`` (one value per requested sample) multiplies each value in
    the lookup table's precision."""
    scaled = scale.size > 0
    for r in range(idx.size):
        j = idx[r]
        l0, l1, l2, missing = lut[r, 0], lut[r, 1], lut[r, 2], lut[r, 3]
        for h in range(rows.size):
            value = _select_call(G[rows[h], j], l0, l1, l2, missing)
            out[r, h] = value * scale[h] if scaled else value


def _gather_samples(G, idx, lut, rows, scale, out):
    """Sample-major storage: each requested sample's calls are read along
    its own row."""
    scaled = scale.size > 0
    for h in range(rows.size):
        i = rows[h]
        for r in range(idx.size):
            value = lut[r, G[i, idx[r]] & 3]
            out[r, h] = value * scale[h] if scaled else value


def _gather_packed(data, idx, lut, rows, scale, out):
    """Two-bit codes (variant-major rows): each requested sample's code is
    read from its byte."""
    scaled = scale.size > 0
    code_lut = np.empty(4, dtype=lut.dtype)
    for r in range(idx.size):
        code_lut[0], code_lut[1], code_lut[2], code_lut[3] = lut[r, 2], lut[r, 3], lut[r, 1], lut[r, 0]
        row = data[idx[r]]
        for h in range(rows.size):
            i = rows[h]
            value = code_lut[(row[i >> 2] >> (2 * (i & 3))) & 3]
            out[r, h] = value * scale[h] if scaled else value


def _project_rows(z, coefficients, Qt, out):
    """out = z - c Q' row by row, the q-term correction summed in a fixed
    order (c_0 Q_0, then + c_k Q_k) along the samples, so a row's values do
    not depend on which other rows are projected with it."""
    k_rows, n = z.shape
    q = Qt.shape[0]
    correction = np.empty(n)
    for r in range(k_rows):
        for i in range(n):
            correction[i] = coefficients[r, 0] * Qt[0, i]
        for k in range(1, q):
            ck = coefficients[r, k]
            for i in range(n):
                correction[i] += ck * Qt[k, i]
        for i in range(n):
            out[r, i] = z[r, i] - correction[i]


if HAS_NUMBA:
    _count_moments = njit(cache=True)(_count_moments)
    _count_calls = njit(cache=True)(_count_calls)
    _count_codes = njit(cache=True)(_count_codes)
    _sums_one = njit(cache=True)(_sums_one)
    _sums_two = njit(cache=True)(_sums_two)
    _sums_weighted = njit(cache=True)(_sums_weighted)
    _accumulate_rows = njit(cache=True)(_accumulate_rows)
    _finish_tile = njit(cache=True)(_finish_tile)
    _moments_tile = njit(cache=True)(_moments_tile)
    _moments_samples_tile = njit(cache=True)(_moments_samples_tile)
    _moments_any_tile = njit(cache=True)(_moments_any_tile)
    _moments_packed_tile = njit(cache=True)(_moments_packed_tile)
    _moments_packed_serial = njit(cache=True)(_moments_packed_serial)
    _moments_packed_parallel = njit(parallel=True, cache=True)(_moments_packed_parallel)
    _decode_packed_rows = njit(cache=True)(_decode_packed_rows)
    _decode_packed_span = njit(cache=True)(_decode_packed_span)
    _decode_packed = njit(cache=True)(_decode_packed)
    _decode_packed_parallel = njit(parallel=True, cache=True)(_decode_packed_parallel)
    _decode_packed_span_all = njit(cache=True)(_decode_packed_span_all)
    _decode_packed_span_parallel = njit(parallel=True, cache=True)(_decode_packed_span_parallel)
    _project_rows = njit(cache=True)(_project_rows)
    _moments_columns_serial = njit(cache=True)(_moments_columns_serial)
    _moments_columns_parallel = njit(parallel=True, cache=True)(_moments_columns_parallel)
    _select_call = njit(inline="always", cache=True)(_select_call)
    _decode_column = njit(cache=True)(_decode_column)
    _decode_variants = njit(cache=True)(_decode_variants)
    _decode_variants_parallel = njit(parallel=True, cache=True)(_decode_variants_parallel)
    _decode_tile = njit(cache=True)(_decode_tile)
    _decode_samples = njit(cache=True)(_decode_samples)
    _decode_samples_parallel = njit(parallel=True, cache=True)(_decode_samples_parallel)
    _gather_variants = njit(cache=True)(_gather_variants)
    _gather_samples = njit(cache=True)(_gather_samples)
    _gather_packed = njit(cache=True)(_gather_packed)


def prepare_moments(G, idx, Q, n_threads=1, sq_weights=None):
    """Float64 means, SDs, covariate coefficients and projected squared norms.

    Inputs are the validated hard calls (an int8 array or
    :class:`~mixmogam._packed.PackedCalls`) and the copied covariate basis
    held by LocoGenotypes. The mean and variance come from exact call
    counts and one pass accumulates the covariate products of the centered
    calls. Each variant is summed sequentially by the same arithmetic, so
    every thread count, layout and storage format gives the same values.
    Without Numba this raises ImportError; ``LocoGenotypes._prepare`` then
    prepares with NumPy, which can differ by reduction rounding.

    ``sq_weights`` (one value per sample) are the squares of a row scale s:
    with ``Q`` already scaled (s Q~), the coefficients are those of the
    scaled values and the norms sum_i s_i^2 z_i^2 - |c|^2. Means and SDs stay
    unweighted.
    """
    n_threads = _validate_n_threads(n_threads)
    if not HAS_NUMBA:
        raise ImportError("compiled preparation requires Numba")
    mean, sd, zz = np.empty(idx.size), np.empty(idx.size), np.empty(idx.size)
    coefficients = np.empty((idx.size, Q.shape[1]))
    packed = isinstance(G, PackedCalls)
    sqw = (np.empty(0) if sq_weights is None
           else np.ascontiguousarray(sq_weights, dtype=np.float64))
    args = (((G.data, G.shape[0]) if packed else (G,))
            + (idx, Q, mean, sd, coefficients, zz, sqw))
    if n_threads == 1:
        (_moments_packed_serial if packed else _moments_columns_serial)(*args)
    else:
        with _numba_thread_limit(min(n_threads, max(idx.size, 1))):
            (_moments_packed_parallel if packed else _moments_columns_parallel)(*args)
    return mean, sd, coefficients, zz


def decode_table(G, idx, lut, out, n_threads=1, compiled=True, start=0):
    """Decode variants ``idx`` into rows of ``out`` from their lookup rows.

    ``lut`` (k, 4) holds each requested variant's standardized values of
    calls 0, 1 and 2, then 0 for a no-call (index ``code & 3``), already in
    the storage precision; ``out`` may be wider (float64 rows hold the same
    values). Its columns are samples ``start`` to ``start + out.shape[1]``.
    The compiled kernels follow the storage layout; without Numba (or with
    ``compiled=False``, for storage that is not an in-memory array) the
    gather is done in bounded NumPy chunks.
    """
    stop = start + out.shape[1]
    if not (HAS_NUMBA and compiled):
        step = max(1, _NUMPY_DECODE_CELLS // max(out.shape[1], 1))
        for first in range(0, idx.size, step):
            codes = np.asarray(G[start:stop, idx[first : first + step]]).T.astype(np.intp) & 3
            out[first : first + codes.shape[0]] = np.take_along_axis(
                lut[first : first + codes.shape[0]], codes, axis=1)
        return out
    if isinstance(G, PackedCalls):
        whole = start == 0 and out.shape[1] == G.shape[0]  # the tight full-width kernels
        if n_threads > 1:
            with _numba_thread_limit(min(_validate_n_threads(n_threads), max(idx.size, 1))):
                if whole:
                    _decode_packed_parallel(G.data, idx, lut, out)
                else:
                    _decode_packed_span_parallel(G.data, idx, lut, out, start)
        else:
            if whole:
                _decode_packed(G.data, idx, lut, out)
            else:
                _decode_packed_span_all(G.data, idx, lut, out, start)
        return out
    sample_major = G.strides[0] > G.strides[1]
    if n_threads > 1:
        n_threads = _validate_n_threads(n_threads)
        kernel = _decode_samples_parallel if sample_major else _decode_variants_parallel
        work = (out.shape[1] + _TILE - 1) // _TILE if sample_major else idx.size
        with _numba_thread_limit(min(n_threads, max(work, 1))):
            kernel(G, idx, lut, out, start)
    else:
        (_decode_samples if sample_major else _decode_variants)(G, idx, lut, out, start)
    return out


def gather_table(G, idx, lut, rows, out, scale=None, compiled=True):
    """:func:`decode_table` for the samples ``rows`` (int64, any order):
    column h of ``out`` receives the value decode_table gives sample
    rows[h] of each variant ``idx``.

    A ``scale`` (one value per requested sample) multiplies each value in
    the storage precision before ``out`` receives it, as LocoGenotypes
    scales decoded rows. The compiled kernels read only the requested
    samples' calls, serially: a gather of a tenth of the samples costs
    less than a parallel decode of all of them. Without Numba (or with
    ``compiled=False``) NumPy gathers the calls of bounded variant chunks.
    """
    if scale is not None:
        scale = np.ascontiguousarray(scale, dtype=lut.dtype)
    if not (HAS_NUMBA and compiled):
        low, high = (int(rows.min()), int(rows.max()) + 1) if rows.size else (0, 0)
        step = max(1, _NUMPY_DECODE_CELLS // max(high - low, rows.size, 1))
        array = isinstance(G, np.ndarray)
        for first in range(0, idx.size, step):
            take = idx[first : first + step]
            if not array:  # other storage reads the samples' span, then picks
                codes = np.asarray(G[low:high, take])[rows - low].T
            elif G.strides[0] > G.strides[1]:  # sample-major: along each sample
                codes = G[np.ix_(rows, take)].T
            else:  # variant-major: along each variant
                codes = G.T[np.ix_(take, rows)]
            values = np.take_along_axis(lut[first : first + take.size],
                                        codes.astype(np.intp) & 3, axis=1)
            if scale is not None:
                values *= scale
            out[first : first + take.size] = values
        return out
    factors = np.empty(0, dtype=lut.dtype) if scale is None else scale  # empty: unscaled
    if isinstance(G, PackedCalls):
        _gather_packed(G.data, idx, lut, rows, factors, out)
    elif G.strides[0] > G.strides[1]:
        _gather_samples(G, idx, lut, rows, factors, out)
    else:
        _gather_variants(G, idx, lut, rows, factors, out)
    return out
