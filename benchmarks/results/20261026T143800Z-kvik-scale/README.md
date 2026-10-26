# Matched comparison with LDAK-KVIK after the speed branch, 26 October 2026

Same protocol, same panel (5,000 x 400,000 founder blocks, seed 20261028),
same traits and seeds and one thread as
[`20261025T153917Z-kvik-scale`](../20261025T153917Z-kvik-scale/README.md),
rerun on master `6a63f215` (the nine-change speed branch, reviewed in
[`20261026-hratt-speed-review`](../20261026-hratt-speed-review/README.md)).
The reference program is the same executable (`~/bin/ldak63`,
SHA256 in the manifest), so its rows measure host noise; the mixmogam
rows measure the branch.

## Resources (wall seconds, peak RSS from `os.wait4`)

| Case | Method | 25 Oct | 26 Oct | Change | RSS 25 Oct | RSS 26 Oct |
|---|---|---:|---:|---:|---:|---:|
| mixed | LDAK-KVIK | 23.7 | 25.7 | (same program) | 1.19 | 1.11 |
| mixed | mixmogam HRATT-HE | 33.8 | **28.9** | 1.17x | 3.10 | **2.61** |
| mixed | mixmogam BOLT-LMM-inf | 51.0 | **34.6** | 1.47x | 2.38 | **2.14** |
| null | LDAK-KVIK | 27.7 | 29.4 | (same program) | 1.18 | 1.06 |
| null | mixmogam HRATT-HE | 29.8 | **22.5** | 1.32x | 3.28 | **2.80** |
| null | mixmogam BOLT-LMM-inf | 50.4 | **31.9** | 1.58x | 2.38 | **2.19** |

HRATT-HE is now level with the reference program on the mixed trait and
faster than it on the null trait. BOLT-LMM-inf gains most (its REML
stage stops early and its Grams are packed). The memory drop is the
packed-triangle Gram cache: the default 1 GB budget now holds the whole
400,000-variant panel, so the partial-cache warning the 25 October run
recorded is gone.

## Statistics are unchanged

Every agreement number reproduces to the digits shown on 25 October:
Spearman 0.9982 (mixed) and 0.9991 (null) for HRATT-HE against the
reference program, the same max abs dlog10 p (16.7 mixed, 0.119 null),
the same 95/100 and 96/100 top-hit overlaps, lambda_GC identical to
four decimals for all three programs, and the same 8 of 10 causal
variants detected on the mixed trait. The speed branch moved no
association result, as its reviews verified; this confirms it end to
end against the external program.

## Caveats

- One host, one day per side: the reference program's own 8% day-to-day
  movement (23.7 to 25.7 s) is the noise floor for cross-program ratios
  here; the mixmogam gains (17-58%) are outside it.
- The panel is a controlled LD-block model, not human demography.
