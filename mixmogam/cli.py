"""Thin ``mixmogam`` entry point: files in, CSV out, the same values as the API.

``mixmogam gwas`` reads standard genotype files plus a phenotype file (and
optional covariate, weight and PC columns) and writes a summary-statistic
CSV/TSV (optionally gzip-compressed); ``mixmogam pca``
writes built-in ancestry components (:mod:`mixmogam.pca`) as a covariate
file; ``mixmogam plot`` turns a result CSV into Manhattan and QQ PNGs
(needs the ``plot`` extra); ``mixmogam hm3`` shows, matches or rebuilds the
bundled HapMap3 model-SNP map. Every value is the library's own: the CSV is
byte-identical to the one the same call in Python writes.
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
    gwas.add_argument("--weights", help="sample ID (or FID IID) then one column of sampling weights (inverse probabilities), "
                              "read like a covariate file")
    gwas.add_argument("--pcs", type=int, default=0, metavar="K",
                      help="append K built-in ancestry PCs to the covariates (min_maf 0.05, "
                           "2000 markers, seed 0; use the pca command for other settings)")
    gwas.add_argument("--method", default="auto", choices=_METHODS,
                      help="auto: exact up to 5000 samples, bolt-inf above (hratt when "
                           "--trait, --weights, --stage1, --model-snps or a "
                           "--relatedness* option is given); default auto")
    gwas.add_argument("--trait", choices=("quantitative", "binary"),
                      help="binary needs --method hratt and 0/1 outcomes")
    gwas.add_argument("--binary-coding", choices=("0/1", "1/2"), default="0/1",
                      help="binary outcome coding; 1/2 means PLINK control/case, with 0/-9 missing")
    gwas.add_argument("--stage1", default="enet",
                      choices=("enet", "ldpred3-inf", "ldpred3-lassosum2", "auto"),
                      help="hratt stage-1 whole-genome predictor (default enet); "
                           "ldpred3-* and auto need the separate ldpred3 package, "
                           "quantitative traits and --stage1 implies --method hratt")
    gwas.add_argument("--model-snps", metavar="FILE|hm3|auto",
                      help="hratt stage-1 model SNPs; stage 2 tests every variant; implies --method "
                           "hratt. 'hm3' (recommended for human data): the panel's HapMap3 variants "
                           "(ldpred3's 1,444,196 HapMap3+, bundled; --hm3-map hm3-ldpred2 for LDpred2's "
                           "1,054,330), matched by rsID and/or position, falling "
                           "back to 'auto' when fewer than half of those on the panel's chromosomes "
                           "(or 300,000) match; 'auto': MAF/LD "
                           "thinning (non-human or unspecified data); or a file of variant IDs, one "
                           "per line. Default: every variant")
    gwas.add_argument("--hm3-build", choices=("auto", "GRCh37", "GRCh38", "rsid"),
                      help="--model-snps hm3: coordinate build for position matching (default: "
                           "--genome-build when it names GRCh37/38, else detected from the overlap; "
                           "'rsid' matches rsIDs only)")
    gwas.add_argument("--hm3-map", metavar="NAME|NPZ",
                      help="--model-snps hm3: a bundled map, 'hm3-ldpred3' (default; ldpred3's "
                           "HapMap3+) or 'hm3-ldpred2' (LDpred2's map_hm3), or a map written by "
                           "'mixmogam hm3 build'")
    gwas.add_argument("--hm3-min-frac", type=float,
                      help="--model-snps hm3: fall back to auto below this matched fraction of the map's "
                           "variants on the panel's chromosomes (default 0.5)")
    gwas.add_argument("--hm3-min-snps", type=int,
                      help="--model-snps hm3: fall back to auto below this many matches (default 300000)")
    gwas.add_argument("--model-snps-include", metavar="FILE",
                      help="variants forced into the model SNPs (with --model-snps): one per line, a "
                           "variant ID or a region chr:pos or chr:start-end; kept out of their own "
                           "chromosome's LOCO scores like every model SNP")
    gwas.add_argument("--no-model-snps-h2-check", action="store_true",
                      help="skip the subset-vs-all-variants HE h2 check (one extra panel pass)")
    gwas.add_argument("--model-snps-stream", choices=("off", "on", "auto"),
                      help="keep the model-SNP subset in an on-disk two-bit cache read sequentially "
                           "by every stage-1 pass instead of in RAM; identical results (default off; "
                           "auto: above --model-snps-stream-threshold)")
    gwas.add_argument("--model-snps-cache", metavar="DIR",
                      help="--model-snps-stream: keep the cache file in DIR and reuse it while the panel, "
                           "the subset and the samples are unchanged (default: temporary, deleted)")
    gwas.add_argument("--model-snps-stream-threshold", type=float, metavar="GB",
                      help="--model-snps-stream auto: stream when the subset needs more than GB "
                           "gigabytes (n k / 4 bytes; default 8 GiB)")
    gwas.add_argument("--relatedness", metavar="FILE|auto",
                      help="hratt: close relatives for the kappa_g correction: a KING / PLINK 2 .kin0 "
                           "table, a GCTA sparse GRM prefix (.grm.sp + .grm.id), or 'auto' (KING-robust "
                           "on common non-model variants, n <= 50,000)")
    gwas.add_argument("--relatedness-threshold", type=float, metavar="PHI",
                      help="--relatedness: keep pairs with KING kinship > PHI (default 0.025; 'auto' "
                           "raises it to 4 x the KING noise SD)")
    gwas.add_argument("--relatedness-correction", choices=("auto", "kappa", "lambda", "none", "gls"),
                      help="hratt: kappa (default with --relatedness), lambda (default without), none, "
                           "or gls (sparse family covariance, opt-in, unweighted quantitative)")
    gwas.add_argument("--kinship-variants", metavar="FILE",
                      help="--relatedness FILE: IDs of the variants the kinship was estimated on "
                           "(one per line, e.g. PLINK .prune.in); warns when more than 10%% are "
                           "stage-1 model SNPs")
    gwas.add_argument("--thin-min-maf", type=float, help="--model-snps auto: MAF threshold (default 0.01)")
    gwas.add_argument("--thin-max-r2", type=float, help="--model-snps auto: r2 threshold (default 0.5)")
    gwas.add_argument("--thin-window-kb", type=float, help="--model-snps auto: window in kb (default 100)")
    gwas.add_argument("--thin-ld-samples", type=int,
                      help="--model-snps auto: samples used for r2, 0 for all (default 10000)")
    gwas.add_argument("--thin-ld-seed", type=int,
                      help="--model-snps auto: seed of the random r2 samples; -1 takes the first "
                           "samples (reads less of a memory-mapped file; default 0)")
    gwas.add_argument("--thin-max-snps", type=int, help="--model-snps auto: cap on model SNPs (default none)")
    gwas.add_argument("--random-state", type=int, default=None,
                      help="two-step method seed (default: the method's own, 0)")
    gwas.add_argument("--n-threads", type=int, default=1,
                      help="parallel genotype preparation and sweeps for the exact and hratt methods "
                           "(values are identical at every count, up to rounding in parallel products)")
    gwas.add_argument("--max-loco-groups", type=int, default=25)
    gwas.add_argument("--out", "-o", required=True, help="summary-statistic .csv or .tsv (optional .gz)")
    gwas.add_argument("--genome-build", help="declared coordinate assembly (e.g. GRCh38 or TAIR10); no liftover")
    gwas.add_argument("--pca-max-memory-gib", type=float, default=4.0,
                      help="working-memory budget for --pcs (default 4 GiB)")
    gwas.add_argument("--pca-solver", choices=("exact", "randomized"), default="exact",
                      help="PCA eigensolver for --pcs (default exact SVD)")

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
    pca.add_argument("--solver", choices=("exact", "randomized"), default="exact",
                     help="PCA eigensolver (default exact SVD)")

    hm3 = sub.add_parser("hm3", help="the bundled HapMap3 model-SNP map: show, match a panel, rebuild",
                         description="info: provenance and size of the map. match: match a "
                                     "genotype panel and print the report (optionally write the "
                                     "matched variant IDs). build: write a map from an .rds "
                                     "(needs the rdata package) or a delimited table with columns "
                                     "chr, pos (GRCh37), pos_hg38, a0, a1, rsid.")
    hsub = hm3.add_subparsers(dest="hm3_command", required=True)
    h_info = hsub.add_parser("info", help="provenance and size of the map")
    h_info.add_argument("--map", help="hm3-ldpred3 (default), hm3-ldpred2, or a map written by "
                                      "'mixmogam hm3 build'")
    h_match = hsub.add_parser("match", help="match a genotype panel and print the report")
    _genotype_args(h_match)
    h_match.add_argument("--map", help="hm3-ldpred3 (default), hm3-ldpred2, or a map written by "
                                       "'mixmogam hm3 build'")
    h_match.add_argument("--build", default="auto", choices=("auto", "GRCh37", "GRCh38", "rsid"))
    h_match.add_argument("--min-maf", type=float, default=0.01,
                         help="drop matches at or below this sample MAF (default 0.01; 0 keeps all)")
    h_match.add_argument("--out", "-o", help="write the matched variant IDs, one per line")
    h_build = hsub.add_parser("build", help="write a map from an LDpred2-style table")
    h_build.add_argument("--input", "-i", required=True, help=".rds or delimited table with a header")
    h_build.add_argument("--out", "-o", required=True, help="output .npz")
    h_build.add_argument("--name", default="custom", help="map name recorded in the file")
    h_build.add_argument("--source", default="", help="source URL/description recorded in the file")
    h_build.add_argument("--licence", default="", help="licence recorded in the file")
    h_build.add_argument("--require-rsid", action="store_true",
                         help="skip rows whose ID is not rs<digits> (default: keep them, matched by "
                              "position only)")

    plot = sub.add_parser("plot", help="Manhattan and QQ plots of a result CSV",
                          description="Read a result CSV/TSV (optional .gz) written by "
                                      "mixmogam gwas or GwasResult.write_csv and write "
                                      "<out-prefix>_manhattan.png and <out-prefix>_qq.png. "
                                      "Requires the plot extra (matplotlib).")
    plot.add_argument("--in", "-i", dest="infile", required=True,
                      help="result CSV/TSV (optional .gz) to plot")
    plot.add_argument("--out-prefix", "-o", required=True,
                      help="output prefix; writes <prefix>_manhattan.png and <prefix>_qq.png")
    plot.add_argument("--alpha", type=float, default=0.05,
                      help="genome-wide significance level of the Manhattan threshold line (default 0.05)")
    plot.add_argument("--no-threshold", action="store_true",
                      help="omit the significance line from the Manhattan plot")
    plot.add_argument("--suggestive", type=float, default=1e-5,
                      help="suggestive significance line on the Manhattan plot (default 1e-5; 0 omits it)")
    plot.add_argument("--label-top", type=int, default=5, metavar="N",
                      help="annotate the N smallest-p variants on the Manhattan plot (default 5; 0 omits)")
    plot.add_argument("--title", help="title for both figures")
    return ap


def _cmd_pca(args) -> int:
    from mixmogam.io.covariates import write_covariates
    from mixmogam.pca import principal_components

    gt = _read_genotypes(args)
    pcs, info = principal_components(gt, k=args.k, min_maf=args.min_maf,
                                     n_markers=args.n_markers, random_state=args.seed,
                                     max_memory_gib=args.max_memory_gib,
                                     solver=args.solver, return_info=True)
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
        columns.append(principal_components(subset, k=args.pcs,
                                            max_memory_gib=args.pca_max_memory_gib,
                                            solver=args.pca_solver))
    X = np.column_stack(columns) if columns else None

    method = args.method
    if method == "auto":
        # --trait, --weights, --stage1, --model-snps and the --relatedness*
        # flags are HRATT-only options, so they move auto to hratt (the
        # Python API raises TypeError for them with method="auto" instead).
        from mixmogam.association import EXACT_N_AUTO
        if (args.trait is not None or w is not None or args.stage1 != "enet"
                or args.model_snps is not None or args.relatedness is not None
                or args.relatedness_threshold is not None
                or args.relatedness_correction is not None
                or args.kinship_variants is not None):
            method = "hratt"
        else:
            method = "exact" if y.size <= EXACT_N_AUTO else "bolt-inf"
    kwargs = {}
    if method in ("exact", "hratt"):
        kwargs["n_threads"] = args.n_threads
    elif args.n_threads != 1:
        raise ValueError(f"--n-threads applies to the exact and hratt methods, not {method}")
    if args.random_state is not None and method != "exact":
        kwargs["random_state"] = args.random_state
    if args.trait is not None:
        kwargs["trait"] = args.trait
    if args.stage1 != "enet":
        kwargs["stage1"] = args.stage1
    if w is not None:
        kwargs["sample_weights"] = w
    thin = {"min_maf": args.thin_min_maf, "max_r2": args.thin_max_r2,
            "window_bp": None if args.thin_window_kb is None else int(round(args.thin_window_kb * 1000)),
            "ld_samples": args.thin_ld_samples, "ld_seed": args.thin_ld_seed,
            "max_snps": args.thin_max_snps}
    thin = {k: v for k, v in thin.items() if v is not None}
    if "ld_samples" in thin and thin["ld_samples"] == 0:
        thin["ld_samples"] = None
    if "ld_seed" in thin and thin["ld_seed"] < 0:
        thin["ld_seed"] = None
    hm3 = {"map": args.hm3_map, "min_frac": args.hm3_min_frac, "min_snps": args.hm3_min_snps,
           "build": args.hm3_build}
    hm3 = {k: v for k, v in hm3.items() if v is not None}
    if args.model_snps == "hm3" and "build" not in hm3 and args.genome_build:
        from mixmogam._hm3 import _BUILD_ALIAS
        b = _BUILD_ALIAS.get(args.genome_build.strip().lower().replace("_", ""))
        if b is not None:
            hm3["build"] = b
    if args.model_snps is not None:
        if args.model_snps in ("auto", "hm3"):
            kwargs["model_snps"] = args.model_snps
            if hm3 and args.model_snps != "hm3":
                raise ValueError("--hm3-* options apply to --model-snps hm3 only")
            opts = {**thin, **hm3}
            if opts:
                kwargs["model_snps_options"] = opts
        else:
            from mixmogam._thin import read_snp_list
            if thin or hm3:
                raise ValueError("--thin-* options apply to --model-snps auto (or the hm3 "
                                 "fallback) and --hm3-* to --model-snps hm3")
            kwargs["model_snps"] = read_snp_list(args.model_snps)
        if args.model_snps_include:
            from mixmogam._thin import read_snp_list
            kwargs["model_snps_include"] = read_snp_list(args.model_snps_include)
        if args.no_model_snps_h2_check:
            kwargs["model_snps_h2_check"] = False
        stream = {"mode": args.model_snps_stream or "off"}
        if args.model_snps_cache is not None:
            stream["cache_dir"] = args.model_snps_cache
        if args.model_snps_stream_threshold is not None:
            stream["threshold_bytes"] = args.model_snps_stream_threshold * 1e9
        if stream["mode"] != "off":
            kwargs["model_snps_stream"] = stream
        elif len(stream) > 1:
            raise ValueError("--model-snps-cache and --model-snps-stream-threshold need "
                             "--model-snps-stream on or auto")
    elif (thin or hm3 or args.model_snps_include or args.no_model_snps_h2_check
          or args.model_snps_stream not in (None, "off") or args.model_snps_cache is not None
          or args.model_snps_stream_threshold is not None):
        raise ValueError("--thin-*, --hm3-*, --model-snps-include, --no-model-snps-h2-check and "
                         "--model-snps-stream/-cache/-stream-threshold need --model-snps")
    if args.relatedness is not None:
        kwargs["relatedness"] = args.relatedness
    if args.relatedness_threshold is not None:
        kwargs["relatedness_threshold"] = args.relatedness_threshold
    if args.relatedness_correction is not None:
        kwargs["relatedness_correction"] = args.relatedness_correction
    if args.kinship_variants is not None:
        if args.relatedness is None:
            raise ValueError("--kinship-variants needs --relatedness")
        kwargs["relatedness_options"] = {"kinship_variants": args.kinship_variants}
    result = gwas(y, subset, X=X, method=method, max_loco_groups=args.max_loco_groups,
                  genome_build=args.genome_build, **kwargs)
    info = result.extra.get("model_snps")
    if info is not None:
        opts = info.get("options") or {}
        if "log" in opts:
            print(opts["log"])
        print(f"stage 1: {info['n']:,} model SNPs ({info['source']})"
              + (f", {info['include']['added']:,} added by --model-snps-include"
                 if "include" in info else ""))
        if "h2_check" in info:
            c = info["h2_check"]
            print(f"stage 1: HE h2 {c['h2_he_subset']:.3f} on the model SNPs, "
                  f"{c['h2_he_full']:.3f} on every variant (ratio {c['ratio']:.2f})")
    result.write_csv(args.out)
    print(f"wrote {len(result)} variants for {result.n} samples to {args.out}")
    return 0


def _cmd_hm3(args) -> int:
    from mixmogam import _hm3

    if args.hm3_command == "info":
        hm = _hm3.load_hm3(args.map)
        print(f"{hm['chr'].size:,} variants on chromosomes "
              f"{', '.join(str(c) for c in np.unique(hm['chr']))}; "
              f"{int(np.count_nonzero(hm['pos38'] < 0)):,} without a GRCh38 position, "
              f"{int(np.count_nonzero(hm['rsid'] < 0)):,} without an rsID, "
              f"{len(hm['siblings']):,} at multi-variant sites")
        print(f"file: {_hm3.resolve_hm3_map(args.map)}")
        print(f"bundled maps: {', '.join(_hm3.HM3_MAPS)} (default {_hm3.DEFAULT_HM3_MAP})")
        for k, v in hm["meta"].items():
            print(f"{k}: {v}")
        return 0
    if args.hm3_command == "build":
        out = _hm3.build_hm3_map(args.input, args.out, name=args.name, source=args.source,
                                 licence=args.licence, require_rsid=args.require_rsid)
        print(f"wrote {out['n']:,} variants to {args.out} ({out['skipped']:,} rows skipped, "
              f"{out['no_grch38']:,} without a GRCh38 position, {out['no_rsid']:,} without an "
              f"rsID, {out['duplicated_rsids']:,} duplicated rsIDs)")
        return 0
    gt = _read_genotypes(args)
    idx, rep = _hm3.match_hm3(gt, build=args.build, map=args.map, min_maf=args.min_maf)
    print(_hm3.format_report(rep))
    if args.out:
        with open(args.out, "w") as fh:
            for v in np.asarray(gt.variant_ids)[idx]:
                fh.write(f"{v}\n")
        print(f"wrote {idx.size:,} variant IDs to {args.out}")
    return 0


def _cmd_plot(args) -> int:
    from mixmogam.plotting import plot_manhattan, plot_qq
    from mixmogam.results import GwasResult

    result = GwasResult.read_csv(args.infile)
    manhattan = f"{args.out_prefix}_manhattan.png"
    qq = f"{args.out_prefix}_qq.png"
    plot_manhattan(result, bonferroni=not args.no_threshold, alpha=args.alpha,
                   suggestive=args.suggestive if args.suggestive > 0 else None,
                   label_top=args.label_top, title=args.title, savepath=manhattan)
    plot_qq(result.p, title=args.title, savepath=qq)
    print(f"wrote {manhattan} and {qq} "
          f"({len(result)} variants, lambda_GC={result.genomic_control():.3f})")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "pca":
            return _cmd_pca(args)
        if args.command == "plot":
            return _cmd_plot(args)
        if args.command == "hm3":
            return _cmd_hm3(args)
        return _cmd_gwas(args)
    except (ValueError, KeyError, OSError, TypeError, ImportError) as exc:
        print(f"mixmogam: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
