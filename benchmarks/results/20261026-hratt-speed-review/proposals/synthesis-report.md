# Speeding up mixmogam HRATT: final recommendation

## 1. Where the time goes and what is realistic

At one thread, the two VB fits take most of an HE run: 23.1 of 31.0 s at 5k x 400k ("p3") and 18.0 of 22.7 s at 50k x 20k ("hapnest"). These are the CV grid fit and the 5-fold cross-fitted LOCO fit.

- **Inside the fits:** the fit's own GEMMs, Numba sweeps and refresh take about 10 s at both sizes. Gram work takes 8.0 s at p3 and 4.2 s at hapnest. At p3 only 2,518 of 3,140 sub-blocks fit the cache, so the rest are rebuilt every LOCO sweep. Decoding takes 5.2 s and 4.2 s and already runs at about its floor (about 0.16 ns per call).
- **Setup:** `prepare_moments` takes 3.4 s and 1.7 s; `_he_alpha` takes 2.8 s and 2.2 s.
- **Tail:** CV predict, out-of-fold LOCO and `_retro_stats` take 2.1 s and 1.1 s.
- **REML:** adds 24 Lanczos passes over the genotypes, 26.9 s at p3 and 19.1 s at hapnest, close to half the run.

The changes below that need no statistical decision save about 10.6 s at p3 and 8.6–9.6 s at hapnest. That is about a third of an HE run and a fifth of a REML run. It brings HRATT-HE level with LDAK-KVIK at 5k x 400k and ahead of it at 50k x 20k. REML stays 1.7–1.9x slower than KVIK unless you approve guarded adaptive Lanczos stopping (C1). With C1, p3-REML falls to about 28 s.

## 2. Recommended changes, ranked by default-path seconds saved (one thread)

Columns: "p3 HE / REML" and "hapnest HE / REML" are projected seconds saved. Equivalence gives the max |Δlog10 p|.

| # | id | change | location | measured ratio | p3 HE / REML | hapnest HE / REML | equivalence | memory | effort |
|---|---|---|---|---|---|---|---|---|---|
| 1 | C1 (**needs approval**) | Guarded adaptive Lanczos stopping in REML; 24 becomes a cap; `slq_tol` (0 = fixed 24) | `_slq.py:lanczos_quadrature_batch` gets a stop callback from `lmm.py:LMM._slq_pack` | SLQ 3.03 wide, 1.55 tall; whole REML 1.36 wide, 1.28 tall | 0 / 17 | 0 / 4.4 | rounding: 2.9e-5 (deflated strong structure), ≤3.4e-6 otherwise; 0 decisions changed | lower (−124 MB peak on tall) | M |
| 2 | C2 (includes C3, C4) | Joint Gram build for the CV and cross-fit engines: one decode, fold syrks summed, one hold-out syrk; packed triangles without slot 0; CV-first budget; packed prange sweep | `twostep.py:_hratt_fit_scores`; `_vb.py:_prepare_grams/_gram/_fold_layout` and the three sweep kernels | Gram work 2.65 (prep 1.97 at n=5k, 3.09 at 50k); whole wide HE 1.13 (1.35 with a p3-like cache), tall HE 1.15 | 5.8 / 6.1 | 3.1 / 3.1 (could be 2–3) | rounding: bit-identical in 18 unweighted runs, 1.4e-13 with sampling weights | global peak lower (−57 / −72 MB measured; p3 about 1.0 → 0.64 GB); CV-phase peak +0.22 GB at p3 | L |
| 3 | C5 | (P, n) layout for R, work and delta when P ≤ ~32; four-column refresh; full-buffer tail multiply; old layout kept above 32 | `_vb.py:VBEngine.fit`, `_refresh_residual`, `_GemmWorkspace`; `raw_slices` must provide a zeroed 128-row buffer | fit 1.07–1.12 (5k), 1.30–1.40 (20k), 1.80 (50k); whole tall HE 1.05 | 0.7 / 0.63 | 4.0 / 4.9 (floor 3.5) | bit-identical on Accelerate; the threaded n ≥ 50k path differs from today's threaded path by 4.4e-6 but now equals serial | neutral | M |
| 4 | C7 | `prepare_moments`: exact counts, centered decode, one BLAS product with Q (Option A). Option B is a bit-identical tiled kernel | `_standardize.py:prepare_moments`; `_loco.py:LocoGenotypes._prepare` | A 2.1 (in situ 1.9–2.2); B 1.27–1.56 | 1.72 / 1.72 (B: 0.8–0.9) | 0.86 / 0.89 (B: 0.4–0.5) | A rounding 2.0e-13 (1.4e-6 on the strong-structure refit path, where any last-digit change does this); B bit-identical | peak unchanged; a 64·n float64 tile per thread, capped by n | M–L |
| 5 | C8 | Two Numba kernels for the rank-q covariate corrections, with a dispatch rule: always for the cross-fit layout, small q only for the CV and refit fits | `_vb.py:VBEngine.fit` (the `column_sets` blocks) | corrections 5.3x (CV 1.6–2.0, cross-fit 10–12.8); whole wide HE 1.047 | 1.11 / 1.08 | 0.12 / 0.12 | rounding: identical for unweighted intercept-only fits, ≤9e-14 with weights or covariates, ≤1.5e-7 on the refit with covariates | nil | S–M |
| 6 | C11 (**precision choice**) | Accumulate out-of-fold LOCO shares inside the cross-fit sweeps; `_out_of_fold_loco` keeps its genotype pass for stage1 | `_vb.py:_refresh_residual` (plus the NumPy fallback), `fit(_own=)`; `twostep.py:_out_of_fold_loco(own=None)` | component 1.048 wide, 1.087 tall | 0.69 / 0.72 | 0.40 / 0.40 | rounding 2.6e-6; 0 decisions changed | at most n·G·8 bytes (measured −7.5 / +1.5 MiB) | S–M |
| 7 | C9 | CV hold-out MSE from a held-rows-only predict (Numba gather-decode; NumPy `np.ix_` fallback) | `twostep.py:_hratt_fit_scores`; `VBEngine.predict(rows=)`; kernels in `_standardize.py` | 5.03 (3.80 at 50k) | 0.67 / 0.66 | 0.33 / 0.33 | bit-identical | −13 MiB transient | M |
| 8 | C12 | `_he_alpha`: one stacked GEMM over all alphas per sample tile, tile size `-(-n // A)` | `twostep.py:_he_alpha` | 1.24 wide, 1.17 tall, 1.11–1.14 at 50k | 0.55 / 0.55 | 0.21 / 0.21 | bit-identical (with `n // A` it is 1.2e-10 relative when n % A == 1) | +0.5 MB | S |
| 9 | C10 (n_threads > 1 only) | Scoped `kmp_set_blocktime(max(prev, 1))` during parallel VB fits; skipped if `KMP_BLOCKTIME` or `OMP_WAIT_POLICY` is set | `_vb.py:VBEngine.fit` | whole wide HE 1.15 at nt=4 (VB 1.17); 1.07 at nt=2 | 0 (2.1 / 2.2 at nt=4) | 0 (0.2 at nt=4) | bit-identical | none; about 2x process CPU time | S |

If C2 is judged too invasive, take C4 instead (the full-column gather kernel: 1.39x on `_gram`, bit-identical, 2.2 / 2.3 s at p3 and 1.4 s at hapnest) plus C3 (packed storage without slot 0). The gram proposer put this bit-identical package at 4.2 / 4.4 s at p3 and 1.5 s at hapnest. Do not add C4 to C2: C2 already contains it, and its compiled gather must use C4's full-column loop order.

## 3. Projected wall time (one thread)

Savings applied, in seconds:
- **p3-HE:** 5.8 + 0.7 + 1.72 + 1.11 + 0.67 + 0.55 = 10.55; with C11 (+0.69) = 11.24.
- **p3-REML:** 6.1 + 0.63 + 1.72 + 1.08 + 0.66 + 0.55 = 10.74; with C11 (+0.72) = 11.46; with C1 (+17) = 28.46.
- **hapnest-HE:** 3.1 + 4.0 + 0.86 + 0.12 + 0.33 + 0.21 = 8.62; with C11 (+0.40) = 9.02.
- **hapnest-REML:** 3.1 + 4.9 + 0.89 + 0.12 + 0.33 + 0.21 = 9.55; with C11 (+0.40) = 9.95; with C1 (+4.4) = 14.35.

| | p3 HE | p3 REML | hapnest HE | hapnest REML |
|---|---|---|---|---|
| profile total now | 31.0 | 56.6 | 22.7 | 43.8 |
| after C2, C5, C7A, C8, C9, C12 | 20.5 | 45.9 | 14.1 | 34.3 |
| + C11 | 19.8 | 45.1 | 13.7 | 33.9 |
| + C1 (guarded) | 19.8 | 28.1 | 13.7 | 29.5 |
| benchmark wall now | 33.8 | ~56 | 27.2 | – |
| benchmark wall minus all savings | 22.6 | ~27.5 | 18.2 | – |
| LDAK-KVIK | 23.7 | 23.7 | 20.1 | 20.1 |

How the overlaps are handled:
- C2 already contains C3 and C4, so they are not added.
- C9 and C11 remove different decode passes (the CV predict pass and the OOF pass), so both count. C13, which removes the Gram-prep decode, is unverified and not counted.
- C5 and C8 change different parts of each sub-block visit (GEMMs and refresh versus the q-corrections), so both count.
- Not counted: C6 (its default does not change, and its saving shrinks as C2, C5 and C8 make CV sweeps cheaper), C10 (does nothing at one thread), and the SLQ-side items C14, C17, C20 and C24 (unverified, and they shrink after C1).

Uncertainty:
- At hapnest, C2 may realize only 2–3 s and C5's floor is 3.5 s. Taking 2 and 3.5 makes hapnest-HE 22.7 − (9.02 − 1.1 − 0.5) = 15.3 s on the profile basis, or 19.8 s on the benchmark basis. That is still at or under KVIK's 20.1 s.
- The benchmark walls are higher than the profile totals at HE, so subtracting seconds from them is the conservative reading. Scaling proportionally instead gives 33.8 × 19.76/31.0 = 21.5 s and 27.2 × 13.68/22.7 = 16.4 s.
- The first run's verdicts (C2, C4, C7) were measured at load average 8–15, and this run's at about 2–13. Whole-run ratios on the reduced panels are noisy, so the projections scale component timings onto the idle full-size profiles. Re-profile at full size on an idle machine after landing C2 and after C5 + C8.
- Multi-core: C10 adds about 2.1 s at p3 with n_threads=4, but nt=4 is still only about 1.29x faster than nt=1 on wide. The underlying cause is two small parallel regions per 128-SNP sub-block. A real fix needs fewer, larger parallel regions.

## 4. Changes that alter the estimator or a numerical default

- **C1, adaptive Lanczos stopping (REML).** It replaces the settled 24 fixed steps with "stop once h2 has converged, cap 24".
  - The proposer's rule must not ship: with n > m it stops while h2 sits on the grid boundary. Example: HAPNEST n/m 4, h2 0.5, where it returned h2 0.99995 against 0.503, max |Δlog10 p| 19.2, λGC 2.351 → 2.187, and 70 of 446 decisions changed at 0.05/m.
  - The guarded rule fixed every case tried: h2 error ≤1e-6 (probe Monte Carlo SD is 2–6e-3), max |Δlog10 p| ≤2.9e-5, 0 decisions changed. Guard: at least 8 passes, two consecutive |Δh2| < 1e-6, and each of the last three checks found a bracketed interior root.
  - It saves nothing where 24 steps are not themselves converged (h2 0.8–0.9 with n > m).
  - **Your decision:** whether `slq_tol` is on by default (17 s at p3-REML, 4.4 s at hapnest-REML) or opt-in. Either way it needs a CHANGELOG entry, and the per-check grid evaluation should be vectorized (29–58 ms per check).
- **C11, OOF scores.** The out-of-fold scores become sums of float32 sweep increments instead of an exact float64 recompute from the final effects. Max |Δlog10 p| is 2.6e-6, λGC moves ≤6e-7, and no call changes. Binary offsets inherit the change.
  - **Your decision:** accept this precision trade. The code comment "Exact residual from the final effects" and the docstrings would change.
- **C7 Option A.** Covariate coefficients change at ≤8e-13 absolute and come to depend on tile composition and BLAS. The docstring promise "each variant sums its samples in order whatever the tile" has to go.
  - **Your decision:** A (about 1.7 s at p3), or the bit-identical B (about 0.85 s).
- **C6, a separate CV tolerance (cv_tol = 1e-4).** The speed reproduces (CV component 1.21–1.50) but the equivalence does not.
  - The grid choice flipped in 4 of 32 wide cases (n = 5,000), always away from the converged choice. Max |Δlog10 p| was 1.45 and at most one decision changed in ordinary runs.
  - On the unconverged strong-structure path, 37–48 decisions changed.
  - Recommendation: add `cv_tol` defaulting to `vb_tol` (no behaviour change) and document 1e-4 as a speed option for large n (tall: 0 flips in 32). If enabled it would save 1.5 / 1.2 s at p3 and 3.1 / 3.1 s at hapnest. A margin-based rule (stop at 1e-4 only when the hold-out winner leads clearly) is untested.
- **C2, not an estimator change.** Record it as rounding-level in the CHANGELOG (1.4e-13 with weights; CV hold-out MSE up to 5e-8 relative, with no grid flip in 19 runs). `gram_cache_bytes` becomes a budget shared by both engines, and the CV phase holds both caches (+0.22 GB at p3, while the global peak falls).
- **Bit-identity depends on Accelerate.** For C5, C9 and C12 it rests on Accelerate's GEMM behaviour. On OpenBLAS or MKL (Linux CI), pin these by tolerance.

## 5. Refuted or unverified ideas

**Refuted (measured or argued):**
- 64-row sub-blocks: net about 0 on wide (0.872 vs 0.865), rounding 2.6e-6, and worse once C2 lands.
- Exact integer Grams in float32: 0.69x at n=5k; breaks with missing calls and non-intercept covariates.
- Float32 syrk with an exact float64 diagonal: 0.59–0.89x, and about 100x less accurate off the diagonal.
- Fewer transcendentals in the posterior mean: already settled in the roadmap (1.8e-7 log10 p for about 0.3 s).
- Sweep micro-optimizations (hoisted log p and similar): about 0.1 s at p3.
- Folding the q-corrections into the sweep and refresh kernels: about 0.12 s, unmeasured; revisit after C5 and C8.
- Warm-starting the cross-fit from the CV effects: leaks held-out-fold information (0.05 log10 p, one flip at 0.05/m).
- Looser LOCO tolerance: adds about 50% to the existing convergence error (0.0077 log10 p).
- C4's proposed 16-row gather kernel: half the saving of the full-column kernel, because NumPy's gather is cheap on the real fold-sorted order.
- Fused decode straight to float64 Grams: ignores `row_scale` (a bug: 0.005 log10 p with weights), fails on packed calls, and C-order output changes rounding.
- C1's unguarded stopping rule: see section 4.
- C5's sample-major `s_proj` copy: not needed for bit-identity.
- Process-wide blocktime: slows the BLAS-threaded SLQ passes (1.69 → 2.07 s).
- Helper-thread prefetch without nogil sweeps: 1.02–1.19x slower.
- Float32 BLAS accumulation of the SLQ backward GEMM: 1.005, no gain.
- Fewer SLQ probes or columns: GEMM cost is nearly flat in column count, so passes are the lever.
- Computing the first REML Lanczos product inside the HE pass: net about 0.
- Float32 cache of decoded rows across REML passes: about 0.9 s at hapnest for about 1 GB.
- Over-relaxed sweeps: LOCO needs more sweeps, and p moves by 0.006–0.02 in log10.
- Fold 0 as the CV hold-out: an estimator change, and C2 dominates it.
- Float32 retro pass: puts float32 error into the final statistic.
- `_GEMM_MIN_SAMPLES = 0` at n=20k: 1.00–1.05x, and betas change at 1e-7.
- Parallel decode of 128-row sub-blocks: 110 → 101 µs; the decode is memory-bound.
- Sample-major decode: the forward GEMM is 4x faster but the backward is 4–6x slower.
- CV hold-out accumulated in the sweeps instead of C9: rounding-level, and no gain over C9.

**Unverified (proposer evidence only):**
- C13, Grams built on first use in sweep 1: about 0.64 / 0.31 s, halved after C2; never prototyped inside fit.
- C14, letting Accelerate thread the SLQ GEMMs (a docs change): about 3–4 s for multi-core REML; the comparison crossed processes, and it is about 60% smaller after C1.
- C15, overlapping CV Gram prep with SLQ: about 3 / 2 s at n_threads > 1; conflicts with C2 on memory.
- C16, parallel Gram prep: 1.2–2.0 s versus nt=4; largely overlaps C2; +307 MB at 50k with 4 workers.
- C17, Numba mixed-precision accumulation in `_products`: not resolved in a whole run; shrinks after C1; its HE part is superseded by C12.
- C18, prefetching the next VB sub-block: 0.85 s at p3 with nt=4; the two proposers' measurements conflict.
- C19, a thread pool over `_he_alpha`'s alphas: superseded by C12.
- C20, a contiguous Krylov basis for reorthogonalization: 0.7 s at hapnest-REML (0.3 s after C1), rounding-level.
- C21, double-buffered `raw_slices`: 0 at p3; hapnest is estimated, not measured.
- C22, dense float32 kinship for small-n REML: an alternative to C1 (about 0 after C1); keep in mind for strong structure.
- C23, freezing and pruning CV grid columns: a heuristic estimator change on top of C6; about 0.4 s.
- C24, dropping the redundant covariate projections in the REML operator: about 0.29 s at hapnest-REML, with a coupling risk.
- C25, the chi2(1) tail via `erfc`: 0.19 s at p3, 1e-13; trivial.
- C26, contiguous variant ranges as views: about 0.13 s, within noise; bundle with C5.
- C27, fused decode-and-dot in `_retro_stats`: about 0.1 s at one thread, load-dependent.
- C28, the refit path predicting only each fold's rows: off the default path; about 2 s at p3 under strong structure; cheap once C9 exists.

## 6. Suggested order

1. **C12:** small, bit-identical, isolated. Add an `np.array_equal` test that includes n % 5 == 1.
2. **C9:** bit-identical. Its gather kernels are reused later by C28 and by an exact fallback for C11. Update the predict-count test in `tests/test_hratt_crossfit.py`.
3. **C2,** with C4's full-column gather, the CV-first budget and the packed prange kernel. It rewrites Gram storage and the sweep kernels, so it should land before the other fit changes. The fallback is C4 + C3.
4. **C5** (dispatch at P ≤ 32), then **C8** (with its dispatch rule). Both edit `VBEngine.fit`. Then re-profile at full size on an idle machine.
5. **C7,** Option A or B according to your decision.
6. **C10:** independent and small.
7. **After your decisions:**
   - C1 is independent of the VB work. If approved, it can go first, since it is the largest REML gain.
   - C11 goes after C5, because both edit `_refresh_residual`.
   - C6 lands as an opt-in `cv_tol` keyword.

Prototypes for each verdict:
- C1: `/tmp/hratt-speed/work/verify2-C1/adapt.py`
- C2: `/tmp/hratt-speed/work/verify-C2/c2.py`
- C4: `/tmp/hratt-speed/work/verify-C4/c4.py`
- C5: `/tmp/hratt-speed/work/verify2-C5/c5v2.py`
- C6: `/tmp/hratt-speed/work/verify2-C6/proto.py`
- C7: `/tmp/hratt-speed/work/verify-C7/` (`proto.py`, `proto2.py`)
- C8: `/tmp/hratt-speed/work/verify2-C8/proto.py`
- C9: `/tmp/hratt-speed/work/verify2-C9/proto.py`
- C10: `/tmp/hratt-speed/work/verify2-C10/bt.py`
- C11: `/tmp/hratt-speed/work/verify2-C11/proto.py`
- C12: `/tmp/hratt-speed/work/verify2-C12/proto.py`

Profiles: `/tmp/hratt-speed/prof/`.