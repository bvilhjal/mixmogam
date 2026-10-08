"""Thin ``mixmogam`` entry point: files in, CSV out, the same values as the API.

``mixmogam gwas`` reads PLINK 1 genotypes plus a phenotype file (and
optional covariate, weight and PC columns) and writes the summary-statistic
CSV of :meth:`~mixmogam.results.GwasResult.write_csv`; ``mixmogam pca``
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


def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="mixmogam", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=f"mixmogam {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    gwas = sub.add_parser("gwas",
                          help="association scan of one trait (LOCO mixed models)",
                          description="Run mixmogam.gwas on files and write its result CSV.")
    gwas.add_argument("--geno", required=True, help="PLINK 1 prefix (.bed/.bim/.fam)")
    gwas.add_argument("--pheno", required=True,
                      help="phenotype file: wide (sample_id + named traits) or long (see read_phenotypes)")
    gwas.add_argument("--pheno-name", help="trait to analyse; required unless the file has exactly one")
    gwas.add_argument("--covar", help="covariate file (see read_covariates)")
    gwas.add_argument("--weights", help="one-column file of sampling weights (inverse probabilities)")
    gwas.add_argument("--pcs", type=int, default=0, metavar="K",
                      help="append K built-in ancestry PCs to the covariates (min_maf 0.05, "
                           "2000 markers, seed 0; use the pca command for other settings)")
    gwas.add_argument("--method", default="auto", choices=_METHODS,
                      help="auto: exact up to 5000 samples, bolt-inf above (hratt when "
                           "--trait or --weights is given); default auto")
    gwas.add_argument("--trait", choices=("quantitative", "binary"),
                      help="binary needs --method hratt and 0/1 outcomes")
    gwas.add_argument("--random-state", type=int, default=None,
                      help="two-step method seed (default: the method's own, 0)")
    gwas.add_argument("--n-threads", type=int, default=1,
                      help="parallel genotype preparation (values are identical at every count)")
    gwas.add_argument("--max-loco-groups", type=int, default=25)
    gwas.add_argument("--out", "-o", required=True, help="summary-statistic CSV")

    pca = sub.add_parser("pca", help="ancestry PCs of the genotypes, as a covariate file",
                         description="Run mixmogam.pca.principal_components and write a covariate file.")
    pca.add_argument("--geno", required=True, help="PLINK 1 prefix (.bed/.bim/.fam)")
    pca.add_argument("--out", "-o", required=True, help="covariate file (columns PC1..PCk)")
    pca.add_argument("--k", type=int, default=10, help="number of components (default 10)")
    pca.add_argument("--min-maf", type=float, default=0.05,
                     help="sample MAF the PC panel is thinned to (default 0.05)")
    pca.add_argument("--n-markers", type=int, default=2000,
                     help="markers in the thinned panel (default 2000)")
    pca.add_argument("--seed", type=int, default=0, help="thinning seed (default 0)")
    return ap


def _cmd_pca(args) -> int:
    from mixmogam.io.covariates import write_covariates
    from mixmogam.io.plink import read_plink
    from mixmogam.pca import principal_components

    gt = read_plink(args.geno)
    pcs, info = principal_components(gt, k=args.k, min_maf=args.min_maf,
                                     n_markers=args.n_markers, random_state=args.seed,
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
    from mixmogam.io.plink import read_plink
    from mixmogam.pca import principal_components

    gt = read_plink(args.geno)
    ph = read_phenotypes(args.pheno)
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

    columns = []
    if args.covar:
        columns.append(read_covariates(args.covar, sample_ids=gt.sample_ids).values)
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
    if w is not None:
        w = w[complete]
    if args.pcs:
        columns.append(principal_components(subset, k=args.pcs))
    X = np.column_stack(columns) if columns else None

    method = args.method
    if method == "auto":
        # Resolve exactly as mixmogam.gwas does (trait= and sample_weights=
        # are HRATT options, so they move auto to hratt).
        from mixmogam.association import EXACT_N_AUTO
        if args.trait is not None or w is not None:
            method = "hratt"
        else:
            method = "exact" if y.size <= EXACT_N_AUTO else "bolt-inf"
    kwargs = {"n_threads": args.n_threads}
    if args.random_state is not None and method != "exact":
        kwargs["random_state"] = args.random_state
    if args.trait is not None:
        kwargs["trait"] = args.trait
    if w is not None:
        kwargs["sample_weights"] = w
    result = gwas(y, subset, X=X, method=method, max_loco_groups=args.max_loco_groups, **kwargs)
    result.write_csv(args.out)
    print(f"wrote {len(result)} variants for {result.n} samples to {args.out}")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    try:
        return _cmd_pca(args) if args.command == "pca" else _cmd_gwas(args)
    except (ValueError, KeyError, OSError, TypeError) as exc:
        print(f"mixmogam: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
