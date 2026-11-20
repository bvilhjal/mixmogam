# Atwell positive controls: prespecified check

Frozen before the first fit, 17 October 2026. Runner:
`benchmarks/atwell_controls.py`.

## Purpose

Simulations cannot contain the plumbing errors of real files: sample IDs
that do not join, positions on another assembly, swapped allele coding,
phenotypes parsed into the wrong rows. Atwell et al. (2010, Nature 465:627)
measured these phenotypes on this panel. For seven of them a single gene
with a large effect is known from mutant and cloning work: RPM1 recognises
AvrRpm1 and AvrB, RPS2 recognises AvrRpt2, and RPS5 recognises AvrPphB.
HKT1 sets leaf sodium, MOT1 molybdenum, and FRI expression is tested at
FRI itself. Finding them is a sanity check. Passing it says nothing about
tail calibration.

## Data

- Genotypes: the bundled 250K RegMap panel (`at_data/at_genotypes.zip`,
  `all_chromosomes_binary.csv`), 1,307 accessions and 214,051 SNPs. Inbred
  0/1 calls are recoded 0/2.
- Phenotypes: `at_data/199_phenotypes.csv`, raw values, no transformation.
  The four `avr` traits are 0/1 hypersensitive responses, analysed as
  quantitative traits.
- Per trait: the accessions with both a genotype and a value. Variants with
  minor allele count >= 10 (at least five accessions under 0/2 coding) and
  at most 10% missing calls in that sample.
- Gene coordinates: TAIR10, from `at_data/gene_annotations.zip`. The panel's
  assembly is not recorded. All its positions lie inside the TAIR10
  chromosome lengths, which fits TAIR8 too.

## Controls

| Id | Trait | Gene | AGI | Chr | Gene body (TAIR10) |
|---|---|---|---|---|---|
| 32 | avrPphB | RPS5 | AT1G12220 | 1 | 4,144,572–4,147,911 |
| 33 | avrRpm1 | RPM1 | AT3G07040 | 3 | 2,225,952–2,229,024 |
| 34 | avrRpt2 | RPS2 | AT4G26090 | 4 | 13,224,338–13,227,872 |
| 35 | avrB | RPM1 | AT3G07040 | 3 | 2,225,952–2,229,024 |
| 44 | FRI (expression) | FRI | AT4G00650 | 4 | 269,026–271,503 |
| 16 | Na23 | HKT1 | AT4G10310 | 4 | 6,391,984–6,395,877 |
| 30 | Mo98 | MOT1 | AT2G25680 | 2 | 10,933,061–10,934,551 |

Reported without a criterion: trait 43, FLC expression, at FLC (AT5G10140,
chr 5, 3,173,497–3,179,448) and at FRI. Both genes regulate FLC
expression, so neither one is the single expected locus.

**Window.** The gene body plus 50 kb on each side. That covers local LD and
a possible offset between assemblies. The distance of the best variant from
the gene body is reported, so a tighter window can be read off the table.

**Recovered.** The smallest p in the window is at most 0.05 / m, where m is
the number of variants tested for that trait.

## Methods

Each with default options and one thread:

- `exact`: EMMAX with per-chromosome LOCO REML. `auto` routes n <= 5,000
  here.
- `bolt-inf`
- `hratt`
- OLS with 10 built-in principal components (`principal_components`,
  default thinning) of the trait's sample.

## Criteria

- **A1.** `exact` recovers all seven controls.
- **A2.** `bolt-inf` and `hratt` each recover every control that `exact`
  recovers. An error counts as a miss.

Reported without a criterion:

- the OLS+PCs recoveries;
- per trait and method, genome-wide lambda_GC;
- per window, the best variant, its distance from the gene body, its p and
  its genome-wide rank;
- the genome-wide top variant;
- the `exact` pseudo-heritability (median over LOCO groups).

Lambda_GC here is not a null calibration. These traits have major genes
and structure.

A miss is investigated and reported. It is not tuned away. Filters,
windows and the transformation stay as written. Any later analysis is
labelled exploratory.
