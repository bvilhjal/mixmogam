"""Thin ``mixmogam`` entry point: files in, CSV out, the same values as the API.

``mixmogam gwas`` reads standard genotype files plus a phenotype file (and
optional covariate, weight and PC columns) and writes a summary-statistic
CSV/TSV (optionally gzip-compressed); ``mixmogam pca``
writes built-in ancestry components (:mod:`mixmogam.pca`) as a covariate
file. Every value is the library's own: the CSV is byte-identical to the
one the same call in Python writes.
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional, Sequence

import numpy as np

from mixmogam import __version__

__all__ = ["main"]

_METHODS = ("auto", "exact", "bolt-inf", "bolt", "hratt")


def _genotype_args(ap):
    ap.add_argument("--geno", required=True, help="genotype file or PLINK/EIGENSTRAT prefix")
    ap.add_argument("--geno-format", default="auto",
                    choices=("auto", "plink", "pgen", "vcf", "bgen", "eigenstrat", "tped", "hdf5"))
    ap.add_argument("--packed", action="store_true", help="keep calls in two-bit storage")
    ap.add_argument("--mmap", action="store_true", help="map PLINK BED calls read-only (implies packed)")
    ap.add_argument("--sample", help="BGEN sample file when sample IDs are not embedded")
    ap.add_argument("--hard-call-threshold", type=float,
                    help="explicit BGEN posterior threshold for converting probabilities to hard calls")


def _read_genotypes(args):
    from mixmogam.io import read_genotypes

    return read_genotypes(args.geno, format=args.geno_format, packed=args.packed,
                          mmap=args.mmap, hard_call_threshold=args.hard_call_threshold,
                          sample_path=args.sample)


def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="mixmogam", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=f"mixmogam {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    gwas = sub.add_parser("gwas",
                          help="association scan of one trait (LOCO mixed models)",
                          description="Run mixmogam.gwas on files and write its result CSV.")
    _genotype_args(gwas)
    gwas.add_argument("--pheno", required=True,
                      help="phenotype file: wide (sample_id + named traits) or long (see read_phenotypes)")
    gwas.add_argument("--pheno-name", help="trait to analyse; required unless the file has exactly one")
    gwas.add_argument("--pheno-format", choices=("auto", "wide", "long", "plink"), default="auto",
                      help="phenotype layout; plink explicitly treats -9 as missing")
    gwas.add_argument("--pheno-id-columns", type=int, choices=(1, 2),
                      help="ID-column count for headerless phenotypes (2 for FID IID)")
    gwas.add_argument("--covar", help="covariate file (see read_covariates)")
    gwas.add_argument("--covar-format", choices=("auto", "wide", "plink"), default="auto",
                      help="covariate missing-value policy; wide retains genuine numeric -9")
    gwas.add_argument("--weights", help="one-column file of sampling weights (inverse probabilities)")
    gwas.add_argument("--pcs", type=int, default=0, metavar="K",
                      help="append K built-in ancestry PCs to the covariates (min_maf 0.05, "
                           "2000 markers, seed 0; use the pca command for other settings)")
    gwas.add_argument("--method", default="auto", choices=_METHODS,
                      help="auto: exact up to 5000 samples, bolt-inf above (hratt when "
                           "--trait or --weights is given); default auto")
    gwas.add_argument("--trait", choices=("quantitative", "binary"),
                      help="binary needs --method hratt and 0/1 outcomes")
    gwas.add_argument("--binary-coding", choices=("0/1", "1/2"), default="0/1",
                      help="binary outcome coding; 1/2 means PLINK control/case, with 0/-9 missing")
    gwas.add_argument("--stage1", default="enet",
                      choices=("enet", "ldpred3-inf", "ldpred3-lassosum2", "auto"),
                      help="hratt stage-1 whole-genome predictor (default enet); "
                           "ldpred3-* and auto need the separate ldpred3 package, "
                           "quantitative traits and --stage1 implies --method hratt")
    gwas.add_argument("--random-state", type=int, default=None,
                      help="two-step method seed (default: the method's own, 0)")
    gwas.add_argument("--n-threads", type=int, default=1,
                      help="parallel genotype preparation (values are identical at every count)")
    gwas.add_argument("--max-loco-groups", type=int, default=25)
    gwas.add_argument("--out", "-o", required=True, help="summary-statistic .csv or .tsv (optional .gz)")
    gwas.add_argument("--genome-build", help="declared coordinate assembly (e.g. GRCh38 or TAIR10); no liftover")
    gwas.add_argument("--pca-max-memory-gib", type=float, default=4.0,
                      help="working-memory budget for --pcs (default 4 GiB)")

    pca = sub.add_parser("pca", help="ancestry PCs of the genotypes, as a covariate file",
                         description="Run mixmogam.pca.principal_components and write a covariate file.")
    _genotype_args(pca)
    pca.add_argument("--out", "-o", required=True, help="covariate file (columns PC1..PCk)")
    pca.add_argument("--k", type=int, default=10, help="number of components (default 10)")
    pca.add_argument("--min-maf", type=float, default=0.05,
                     help="sample MAF the PC panel is thinned to (default 0.05)")
    pca.add_argument("--n-markers", type=int, default=2000,
                     help="markers in the thinned panel (default 2000)")
    pca.add_argument("--seed", type=int, default=0, help="thinning seed (default 0)")
    pca.add_argument("--max-memory-gib", type=float, default=4.0,
                     help="PCA working-memory budget (default 4 GiB)")
    return ap


def _cmd_pca(args) -> int:
    from mixmogam.io.covariates import write_covariates
    from mixmogam.pca import principal_components

    gt = _read_genotypes(args)
    pcs, info = principal_components(gt, k=args.k, min_maf=args.min_maf,
                                     n_markers=args.n_markers, random_state=args.seed,
                                     max_memory_gib=args.max_memory_gib,
                                     return_info=True)
    write_covariates(args.out, gt.sample_ids, pcs,
                     names=[f"PC{j}" for j in range(1, args.k + 1)])
    print(f"wrote {args.k} principal components for {gt.n_samples} samples to {args.out} "
          f"(from {info['variants'].size} variants)")
    return 0


def _cmd_gwas(args) -> int:
    from mixmogam.association import gwas
    from mixmogam.io.covariates import read_covariates
    from mixmogam.io.phenofile import read_phenotypes
    from mixmogam.pca import principal_components

    if args.binary_coding == "1/2" and args.trait != "binary":
        raise ValueError("--binary-coding 1/2 requires --trait binary")
    gt = _read_genotypes(args)
    ph = read_phenotypes(args.pheno, format=args.pheno_format, id_columns=args.pheno_id_columns)
    names = ph.pids()
    if args.pheno_name is not None:
        if args.pheno_name not in ph:
            raise ValueError(f"unknown trait {args.pheno_name!r}; the file has {names}")
        trait_name = args.pheno_name
    elif len(names) == 1:
        trait_name = names[0]
    else:
        raise ValueError(f"the phenotype file has {len(names)} traits {names}; choose one with --pheno-name")
    y = ph.align(gt.sample_ids, trait_name)
    if args.binary_coding == "1/2":
        if np.any(np.isfinite(y) & ~np.isin(y, [-9, 0, 1, 2])):
            raise ValueError("PLINK binary phenotypes must be 1/2, with 0/-9 or NA for missing")
        y = np.where((y == 1) | (y == 2), y - 1, np.nan)

    columns = []
    if args.covar:
        columns.append(read_covariates(args.covar, sample_ids=gt.sample_ids, format=args.covar_format).values)
    w = None
    if args.weights:
        weights = read_covariates(args.weights, sample_ids=gt.sample_ids)
        if len(weights.names) != 1:
            raise ValueError(f"--weights needs exactly one column, found {list(weights.names)}")
        w = weights.values[:, 0]

    # Complete-case samples: the API fits finite phenotypes and positive
    # weights only (covariates are checked finite by the reader).
    complete = np.isfinite(y)
    if w is not None:
        complete &= np.isfinite(w) & (w > 0)
    dropped = int(np.count_nonzero(~complete))
    if not complete.any():
        raise ValueError("no samples have a phenotype (and a positive weight)")
    if dropped:
        print(f"using {int(complete.sum())} of {gt.n_samples} samples "
              f"({dropped} without a phenotype value or with a non-positive weight)")
    if np.count_nonzero(complete) == gt.n_samples:
        subset = gt
    else:
        subset = gt.filter_samples(np.flatnonzero(complete))
    y = y[complete]
    columns = [values[complete] for values in columns]
    if w is not None:
        w = w[complete]
    if args.pcs:
        columns.append(principal_components(subset, k=args.pcs, max_memory_gib=args.pca_max_memory_gib))
    X = np.column_stack(columns) if columns else None

    method = args.method
    if method == "auto":
        # Resolve exactly as mixmogam.gwas does (trait=, sample_weights= and
        # stage1= are HRATT options, so they move auto to hratt).
        from mixmogam.association import EXACT_N_AUTO
        if args.trait is not None or w is not None or args.stage1 != "enet":
            method = "hratt"
        else:
            method = "exact" if y.size <= EXACT_N_AUTO else "bolt-inf"
    kwargs = {"n_threads": args.n_threads}
    if args.random_state is not None and method != "exact":
        kwargs["random_state"] = args.random_state
    if args.trait is not None:
        kwargs["trait"] = args.trait
    if args.stage1 != "enet":
        kwargs["stage1"] = args.stage1
    if w is not None:
        kwargs["sample_weights"] = w
    result = gwas(y, subset, X=X, method=method, max_loco_groups=args.max_loco_groups,
                  genome_build=args.genome_build, **kwargs)
    result.write_csv(args.out)
    print(f"wrote {len(result)} variants for {result.n} samples to {args.out}")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    try:
        return _cmd_pca(args) if args.command == "pca" else _cmd_gwas(args)
    except (ValueError, KeyError, OSError, TypeError, ImportError) as exc:
        print(f"mixmogam: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
