"""Compare packed BLAS with a four-state Mailman Gram product.

Grouping consecutive variants encodes each sample as a base-four pattern.
The fourth digit is a missing hard call with standardized value zero. For
Z P, sum the right-hand sides by pattern, then successively marginalize one
digit. For Z' T, expand a table of contributions and gather by pattern.

The experiment includes covariate projection, row scaling and SNP weights.
Factoring the sample scale out of float32 genotype decoding changes the
rounding order. The recorded relative error is the Frobenius norm of the
product difference divided by the norm of the shipping product.

Liberty and Zucker, The Mailman algorithm, Yale technical report 1402.
https://cs.yale.edu/homes/el327/papers/mailmanAlgorithm.pdf
"""

import argparse
import json
from pathlib import Path
import sys
from time import perf_counter
from statistics import median
import numpy as np
from numba import njit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mixmogam import Genotypes
from mixmogam._loco import LocoGenotypes
from _harness import power_state


@njit
def product(codes, lut, snp_weights, P, digits):
    blocks, n = codes.shape
    width = P.shape[1]
    hist = np.empty((4**digits, width))
    lookup = np.empty_like(hist)
    out = np.zeros((n, width))
    for group in range(blocks):
        begin = group * digits
        b = min(digits, lut.shape[0] - begin)
        count = 4**b
        hist[:count, :] = 0
        for i in range(n):
            for c in range(width):
                hist[codes[group, i], c] += P[i, c]
        # Summing out the low digit leaves a smaller histogram for the next row.
        scores = np.zeros((b, width))
        for j in range(b):
            for prefix in range(count // 4):
                for c in range(width):
                    subtotal = 0.
                    for h in range(4):
                        value = hist[4 * prefix + h, c]
                        subtotal += value
                        scores[j, c] += lut[begin + j, h] * value
                    hist[prefix, c] = subtotal
            scores[j] *= snp_weights[begin + j]
            count //= 4
        # Reverse the recursion to get the score for each original pattern.
        lookup[0] = 0
        for j in range(b - 1, -1, -1):
            for prefix in range(count - 1, -1, -1):
                for c in range(width):
                    before = lookup[prefix, c]
                    for h in range(4):
                        lookup[4 * prefix + h, c] = before + lut[begin + j, h] * scores[j, c]
            count *= 4
        for i in range(n):
            for c in range(width):
                out[i, c] += lookup[codes[group, i], c]
    return out


def encode(gt, digits):
    codes = np.zeros(((gt.n_variants + digits - 1) // digits, gt.n_samples), np.uint16)
    for block, call in enumerate(codes):
        begin = block * digits
        G = np.asarray(gt.G[:, begin:begin + digits])
        for j in range(G.shape[1]):
            call |= ((G[:, j] & 3).astype(np.uint16) << (2*j))
    return codes


def elapsed(fn, *args):
    t = perf_counter()
    value = fn(*args)
    return perf_counter() - t, value


def main():
    power_state()
    rng = np.random.default_rng(839)
    parser = argparse.ArgumentParser()
    parser.add_argument("--tall", action="store_true", help="test 65,537 samples and 257 variants")
    parser.add_argument("--out", type=Path, help="write all paired repetitions as JSON")
    args = parser.parse_args()
    n, m = (65537, 257) if args.tall else (8193, 1025)
    G = rng.binomial(2, rng.uniform(.05, .5, size=m), size=(n, m)).astype('i1')
    G[rng.random(G.shape) < .015] = -1
    gt = Genotypes(G, chromosome=1 + (np.arange(m) >= m//2), packed=True)
    Q, _ = np.linalg.qr(np.column_stack([np.ones(n), rng.normal(size=(n, 3))]))
    scale = rng.uniform(.7, 1.3, n)
    lg = LocoGenotypes(gt, (np.arange(m) >= m // 2).astype(int), Q, row_scale=scale)
    weights = rng.uniform(.5, 1.5, m)
    rows = []
    for digits in (4, 5, 6):
        encode_sec, codes = elapsed(encode, gt, digits)
        for width in (1, 8, 48, 140):
            P = rng.normal(size=(n, width))
            def mailman():
                projected = lg.project(P).astype(np.float32) * lg._scale_store[:, None]
                raw = product(codes, lg._lut, weights.astype(np.float32), projected, digits)
                return lg.project(raw * lg._scale_store[:, None]) / weights.sum()
            def packed():
                return lg.matmul(P, weights)
            dst, ref = mailman(), packed()
            error = np.linalg.norm(dst - ref) / np.linalg.norm(ref)
            assert np.isfinite(error), error
            timings = [[], []]
            for rep in range(3):
                for k in ((0, 1) if rep % 2 == 0 else (1, 0)):
                    sec, _ = elapsed((packed, mailman)[k])
                    timings[k].append(sec)
            rows.append(dict(n=n, m=m, digits=digits, width=width, encode=encode_sec,
                bytes=codes.nbytes, packed=timings[0], mailman=timings[1], relative_error=error,
                meets_tolerance=bool(error < 4e-6)))
            print(digits, width, median(timings[0]), median(timings[1]), encode_sec, error, flush=True)
    if args.out is not None:
        args.out.write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
