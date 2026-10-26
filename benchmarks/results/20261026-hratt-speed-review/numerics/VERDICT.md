# Review: numerics (base `dba1ce8d` vs branch head `63c2afb6`)

Verification of the whole nine-change branch against master, one thread,
fresh process per fit, same panels and seeds. Sixteen configurations
covering every HRATT path (HE and REML, cross-fit, strong-structure
refit, in-sample), BOLT-LMM-inf, covariates, two-bit storage, missing
calls, sampling weights, binary traits and `alpha_method="reml"`.
Rows 1-10 were produced before the session limit; rows 11-16 were
completed on resume. Raw rows: `matrix_t1.jsonl`.

| # | case | identical p | max abs Δlog10 p | decisions changed at 0.05/m |
|---|---|---|---|---:|
| 1 | wide he | yes | 0 | 0 |
| 2 | wide reml | no | 1.5e-6 | 0 |
| 3 | tall he | yes | 0 | 0 |
| 4 | tall reml | no | 5.1e-6 | 0 |
| 5 | wide-cov he | no | 5.7e-14 | 0 |
| 6 | wide-cov reml | no | 5.8e-7 | 0 |
| 7 | wide-packed he | yes | 0 | 0 |
| 8 | wide-missing reml | no | 6.8e-7 | 0 |
| 9 | tall-weights he | no | 1.2e-13 | 0 |
| 10 | tall-binary he | yes | 0 | 0 |
| 11 | tall-binary reml | no | 8.2e-6 | 0 |
| 12 | struct he | yes | 0 | 0 |
| 13 | struct reml | no | 3.7e-5 | 0 |
| 14 | struct-pc3 reml | no | 8.5e-7 | 0 |
| 15 | wide reml-alpha | no | 6.3e-7 | 0 |
| 16 | tall bolt-inf | no | 7.2e-6 | 0 |

## Verdict: pass

- **No call at 0.05/m changed in any of the 16 configurations**, and the
  genomic-control lambda agrees to six decimals everywhere.
- Every non-zero difference sits inside the band the per-change
  verifications already documented: the REML rows move at 1e-6 to 8e-6
  from the adaptive Lanczos stopping (`slq_tol`, whose own verification
  bounded the strong-structure refit at 3.7e-5 — exactly what case 13
  reproduces), the covariate and weight rows at 1e-13 to 1e-7 from
  reassociated covariate corrections, and the HE-only and packed rows
  are bit-identical.
- `max_rel_beta` reaches 0.069 on case 11, but that is relative error on
  effects whose p-values agree to 8e-6; the decisions and lambda are
  unchanged. The same character appears in the per-change reviews.
- The fits also run faster and leaner in these correctness runs (for
  example wide `reml-alpha` 20.4 s to 10.8 s, peak RSS 866 to 529 MB),
  which corroborates the performance review on the same panels.
- `slq_tol=0` restores the fixed 24 steps bit for bit (per the change's
  own test), so the one numerical default that moved has an exact
  opt-out.

No finding requires a code change.
