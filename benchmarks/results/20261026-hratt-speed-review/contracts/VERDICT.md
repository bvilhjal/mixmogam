# Review: contracts (branch head `63c2afb6`)

Scope: the public surface and its documentation across the nine commits
`f2b993b0..63c2afb6`, on top of master `dba1ce8d`. The full diff is
`vb.diff` (26 files, +3188/-321).

## Checks

1. **Full test suite.** `full_suite.log`: **1358 passed, 6 skipped**
   (optional reader dependencies absent: bgen_reader, zstandard, pysam),
   0 failed. Includes slow and integration tests. Ruff was clean in each
   change's own verification.

2. **CHANGELOG.** Every one of the nine changes has an entry under
   `[Unreleased]`, each with its measured effect and its numerical
   character. The two entries that need care are honest:
   - `slq_tol` (default 1e-6) **changes a numerical default** — the
     entry says so, gives the p-value movement band (8e-6, 3.7e-5 on the
     strong-structure refit) against the probe Monte Carlo SD of 2-6e-3,
     and documents `slq_tol=0` as the bit-exact opt-out. The unguarded
     rule's failure mode (h2 stuck at the grid end when n > m) is
     recorded as the reason for the root guard.
   - `cv_tol` is opt-in, defaults to `vb_tol`, and the entry states
     plainly that `cv_tol=1e-4` changes the estimator (4 of 32 grid
     choices flipped on 5,000 x 40,000), so no one can mistake it for a
     free speed-up.

3. **Public API.** New keywords: `cv_tol` (`hratt`), `slq_tol`
   (`LMM.fit`, `hratt`, `bolt_inf`, `bolt`), `VBEngine.predict(beta,
   rows=...)`. All default to the previous behaviour (`cv_tol=None`,
   `slq_tol` restores 24 steps only at 0 — the one deliberate default
   change, documented above). `gram_cache_bytes` changes meaning
   (one shared budget, packed triangles, hold-out Grams cached first,
   ~620,000 variants per GiB instead of ~320,000); README and
   `docs/design.md` were updated to match, and the old 400,000-variant
   archive text was rewritten in the past tense so it stays true.

4. **Docs.** README (Gram-cache numbers), `docs/design.md` (Gram
   preparation, packed triangles, adaptive stopping), `docs/quickstart.md`
   (the spinning-workers CPU cost and how to opt out) all match the code.
   `benchmarks/slq_defaults.py` and `benchmarks/scale_study.py` were
   adjusted for the shared cache budget.

5. **Tests.** Six new test files (`test_slq_adaptive.py`,
   `test_hratt_cv_tol.py`, `test_vb_gram_preparation.py`,
   `test_vb_covariate_corrections.py`, `test_vb_column_layout.py`,
   `test_predict_rows.py`) plus updates to seven existing ones. The new
   tests pin the guards (the root-bracket condition, the hold-out-only
   scope of `cv_tol`, the packed-triangle reads) rather than just the
   happy path.

## Verdict: pass

No contract break: every new knob preserves the previous result at its
default except `slq_tol`, whose default change is documented, bounded,
reversible and guarded. No finding requires a code change.

One note for the future, not a defect: `test_decode_and_whitening.py`
and `test_spa.py` ask for 2 threads and fail when the whole environment
is pinned to `NUMBA_NUM_THREADS=1`; they fail the same way on master.
