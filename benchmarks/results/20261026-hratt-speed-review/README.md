# HRATT speed-up review, 26 October 2026

Three independent reviews of the nine-change speed branch
(`f2b993b0..63c2afb6`, on top of `dba1ce8d`) before it landed on
master. The changes themselves, with their component measurements, are
in the `[Unreleased]` section of the CHANGELOG; this archive records the
review evidence.

The reviews were started by the implementation workflow and interrupted
by a session limit; they were completed on resume from their own
artifacts (`review-numerics/`, `review-contracts/`,
`review-performance/` under `/tmp/hratt-speed`, copied here).

## Numerics (`numerics/`)

Whole branch against master, one thread, fresh process per fit, sixteen
configurations covering every HRATT path plus BOLT-LMM-inf, covariates,
two-bit storage, missing calls, sampling weights and binary traits.

**Pass.** No call at 0.05/m changed in any configuration; the genomic
control lambda agrees to six decimals everywhere; the largest p-value
movement is 3.7e-5 in log10, on the strong-structure REML refit path,
exactly the band the adaptive Lanczos stopping documented. Verdict and
raw rows: `numerics/VERDICT.md`, `numerics/matrix_t1.jsonl`.

## Contracts (`contracts/`)

**Pass.** The full suite passed **1358 tests, none failed** (6 skips are
absent optional readers). Every change has a CHANGELOG entry with its
measured effect; the one numerical default that moved (`slq_tol`,
default 1e-6) is documented, bounded against its Monte Carlo error, and
reversible with `slq_tol=0`; `cv_tol` is opt-in and its entry states
that 1e-4 changes the estimator. README, design notes and quickstart
match the code, including the new `gram_cache_bytes` meaning (~620,000
variants per GiB). Verdict and suite log: `contracts/`.

## Performance (`performance/`)

Paired A/B, seven pairs per case, fresh processes, load-gated.

| Case | Threads | Base s | New s | Ratio |
|---|---|---:|---:|---:|
| wide (5,000 x 40,000) HE | 1 | 5.20 | 3.88 | 1.33 |
| wide REML | 1 | 8.26 | 5.04 | 1.64 |
| tall (20,000 x 10,000) HE | 1 | 3.59 | 2.31 | 1.56 |
| tall REML | 1 | 8.45 | 5.28 | 1.62 |
| wide HE | 4 | 5.10 | 3.18 | 1.61 |
| tall HE | 4 | 3.73 | 2.24 | 1.61 |

**Pass.** Every configuration is at least 1.33x faster; the four-thread
ratios now match the one-thread ones, so the worker-launch latency is
gone. Verdict, raw timings and driver log: `performance/`.

## Outcome

The branch landed on master unchanged as a fast-forward to `63c2afb6`.
No review finding required a code change. The two documented caveats
remain open by design: `slq_tol` is a numerical default change with an
exact opt-out, and the spinning-workers setting costs 1.5-2x process
CPU time during parallel fits (documented in the quickstart, honoured by
the OpenMP runtime where it already spins).
