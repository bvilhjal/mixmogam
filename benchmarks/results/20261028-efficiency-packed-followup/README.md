# Packed sample selection: paired follow-up

Four alternating fresh-process pairs compare the original source with the
final direct two-bit selection and in-place alignment mask. The archived
protocol, source trees, hashes, power state and resource measurements are
retained here. The synthetic panel has 12,000 samples and 6,000 variants,
including 1% absent IDs in the requested order.

The median paired time ratio was 1.42 (median 0.240 versus 0.167 s); the
four pair ratios were 1.32, 1.52, 1.42 and 1.41. Median process peak RSS was
283.8 versus 220.7 MiB. All packed output bytes and found-sample counts
matched. These are alignment measurements on one synthetic panel, not a
claim about full GWAS speed.

`gwas-parity/` contains five additional complete public-HRATT comparisons
against these frozen source trees (before the later curvature safeguard), including weighted, binary, packed and
strong-structure refit paths. All variational fits converged, their chosen
models and iteration counts agreed, and the maximum absolute log10 p
change was 1.44e-13. No decision at 0.05 / m changed; masks for NaN and
infinite values matched. This is software parity, not a calibration study.
