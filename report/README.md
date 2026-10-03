# Historical report status

`mixmogam_report.tex`, its PDF, figures and tables describe the archived
pre-review implementation. They were not regenerated during the
`2.0.0.dev1` cleanup. The external LDAK comparison used incorrect PLINK
encoding and must be rerun before its claims are cited. Genotype filtering,
missingness handling and default permutations have also changed.

The current [critical review](../docs/review-2026-10-03.md) takes precedence
over the historical report and `docs/methods.pdf`. Rebuilding figures from
the same old CSVs would not validate the corrected implementation.
