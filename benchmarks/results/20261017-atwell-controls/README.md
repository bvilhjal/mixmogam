# Atwell positive controls, 17 October 2026

Prespecified in [plan.md](plan.md), which was written before the first fit.
Run with `benchmarks/atwell_controls.py` on mixmogam 2.0.0.dev6 at
59a4d1bd, with the uncommitted streaming `read_regmap` (git status in
`environment.json`). One thread. The run took 774 s with 1.8 GiB peak
resident memory.

**A1 failed: `exact` recovered six of seven controls. A2 passed: `bolt-inf`
and `hratt` recovered the same six.** OLS with 10 principal components also
recovered those six. No method recovered MOT1 for molybdenum.

| Trait | Gene | n | exact | bolt-inf | hratt | OLS+PCs |
|---|---|---|---|---|---|---|
| avrPphB | RPS5 | 90 | 2e-22 (1) | 3e-13 (1) | 5e-11 (1) | 6e-18 (1) |
| avrRpm1 | RPM1 | 84 | 8e-16 (1) | 2e-10 (2) | 2e-12 (1) | 3e-15 (1) |
| avrRpt2 | RPS2 | 89 | 2e-16 (1) | 7e-12 (1) | 1e-11 (1) | 3e-15 (1) |
| avrB | RPM1 | 87 | 7e-13 (2) | 3e-09 (2) | 1e-09 (1) | 8e-13 (1) |
| FRI expression | FRI | 163 | 3e-10 (1) | 8e-09 (1) | 2e-08 (1) | 1e-11 (1) |
| Na | HKT1 | 93 | 4e-10 (1) | 6e-08 (2) | 7e-09 (1) | 5e-10 (1) |
| **Mo** | **MOT1** | 93 | **3e-06 (6)** | **1e-05 (12)** | **1e-06 (10)** | **2e-06 (5)** |
| FLC expression\* | FLC | 166 | 2e-03 (684) | 8e-03 (2674) | 9e-03 (1163) | 4e-03 (2645) |
| FLC expression\* | FRI | 166 | 5e-03 (1513) | 1e-02 (3835) | 6e-03 (763) | 3e-03 (2153) |

Each cell gives the smallest p within 50 kb of the gene body and, in
parentheses, that variant's genome-wide rank. The thresholds 0.05/m are
2.4e-7 to 2.5e-7, with m = 201,844 to 211,693 tested variants. Bold marks
a miss. \*No criterion: two genes regulate FLC expression.

- **Where the signals sit.** The best window variant lies inside the gene
  body for RPM1 and RPS2 under every method, and for RPS5 under `exact`.
  The other three methods' best RPS5 variant is 1.4 kb outside it. For FRI
  it is 36 bp upstream. For HKT1 it is 7.0 kb from the gene body, or 0.8 kb
  under OLS+PCs. Every phenotyped accession is in the panel. So the
  sample join and the positions hold on the real files. Allele coding
  cannot be checked this way: swapping it leaves p unchanged.
- **avrB.** `exact` and `bolt-inf` rank chr3:2,337,838 first, 109 kb past
  RPM1. `bolt-inf` does the same for avrRpm1. The RPM1 variant
  (chr3:2,227,817) is second in each case.
- **Mo.** The best MOT1 window variant is chr2:10,942,360, 7.8 kb from
  the gene body, at p = 2.8e-6 under `exact`. The strongest
  `exact` signal, chr2:9,604,507 at p = 1.2e-7, clears the threshold
  1.3 Mb from MOT1. See
  the exploratory section.
- **FLC expression.** Neither window is near the threshold. `exact` and
  `bolt-inf` rank chr1:6,148,689 first, at p = 3.9e-9 and 1.5e-15.
  `bolt-inf`'s lambda_GC is 0.87 on that trait.
- **lambda_GC** (`scans.csv`, not a null calibration): `exact` 0.99–1.04,
  `bolt-inf` 0.87–1.05, `hratt` 0.93–1.20, OLS+PCs 1.08–1.22.
- **Heritability.** The median LOCO REML pseudo-heritability under `exact`
  is at its upper bound (0.99995) for all four 0/1 `avr` traits. It is 0.45
  for FRI, 0.76 for Na, 0.93 for Mo and 0.92 for FLC.
- **Time.** `hratt` took 11–144 s per trait, against under 2 s for every
  other method. The slow fits did not converge (exploratory section).
  Python showed that warning once per call site, so `run.log` has only the
  first.

## Exploratory follow-ups (not prespecified)

`exploratory/explore.py`, run after the results above, with the run's
definitions otherwise unchanged. A1 stays failed.

- **MOT1 is found on the log scale.** With log10 Mo values, the best MOT1
  window variant is chr2:10,934,335, inside the gene body. It ranks first
  genome-wide (p = 6.4e-10). On raw values the best window p was 2.8e-6,
  and 2.0e-6 with a minor-allele count of two accessions instead of five
  (`mo98.csv`). So the miss comes from the prespecified raw scale, not
  from the frequency filter. Baxter et al. (2008) trace low shoot Mo to a
  MOT1 promoter deletion in seven of the accessions with the lowest Mo.
  Here seven accessions carry the minor allele at chr2:10,934,335, with
  median Mo 0.47 (0.26–1.05) against 2.53 (0.34–10.31) for the rest. On
  the raw scale the high tail of the common allele (5–10) dominates the
  variance; the log scale does not let it. A future real-data plan should
  fix each trait's transformation in advance rather than treat raw values
  as neutral.
- **HRATT fails to converge when h2 reaches its upper bound.** That
  happened for the four `avr` traits and FLC expression, where HRATT's h2
  was 0.992–0.99995: the cross-validation and LOCO variational fits ran to
  the 100-iteration cap. Those fits took 101–143 s. With h2 0.47–0.88 (FRI,
  Na, Mo) both fits converged in 6–19 iterations and took 12–30 s
  (`hratt_convergence.csv`). The unconverged fits still found every
  control `exact` found. At these sample sizes (84–166) `auto` uses
  `exact`.

## Files

- `windows.csv`: one row per trait, gene window and method.
- `scans.csv`: one row per trait and method (n, m, lambda_GC, top variant,
  seconds, `exact` heritability).
- `criteria.json`: A1, A2 and the recovered sets.
- `top_hits.npz`: the 1,000 smallest p per scan, with variant IDs.
- `run.log`: the console output.
- `exploratory/`: the follow-up scripts (`explore.py`, `mot1_alleles.py`),
  the log, `mo98.csv` and `hratt_convergence.csv`.
- `environment.json`: versions, threads, input hashes and the snapshot.
- `source/`: the package and the two benchmark scripts that ran. The
  in-tree driver later gained two things: it copies itself into the
  archive, and it records HRATT's `cv_converged` and `loco_converged`.
  Neither changes a p-value.

The panel is not copied. It is `at_data/at_genotypes.zip` (sha256 in
`environment.json`).
