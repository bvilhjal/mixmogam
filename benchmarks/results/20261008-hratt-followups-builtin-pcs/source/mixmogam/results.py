"""GWAS results container: tables, priors/PPA, genomic control, power."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy import special, stats

__all__ = ["GwasResult"]

_CHI2_1DF_MEDIAN = float(stats.chi2.ppf(0.5, 1))  # 0.4549


@dataclass
class GwasResult:
    """Association scan results with variant metadata."""

    chromosome: np.ndarray
    position: np.ndarray
    p: np.ndarray
    variant_ids: Optional[np.ndarray] = None
    f_stat: Optional[np.ndarray] = None
    beta: Optional[np.ndarray] = None
    se: Optional[np.ndarray] = None
    rss: Optional[np.ndarray] = None
    af: Optional[np.ndarray] = None
    var_perc: Optional[np.ndarray] = None
    extra: dict = field(default_factory=dict)
    effect_allele: Optional[np.ndarray] = None
    other_allele: Optional[np.ndarray] = None
    n: Optional[int] = None  # samples analysed; exported as the CSV ``n`` column

    def __post_init__(self):
        if self.n is not None:
            if isinstance(self.n, (bool, np.bool_)) or not isinstance(self.n, (int, np.integer)) or self.n <= 0:
                raise ValueError("n must be a positive sample count")
            self.n = int(self.n)
        self.p = np.asarray(self.p, dtype=float)
        if self.p.ndim != 1:
            raise ValueError("p must be one-dimensional")
        for name in ("chromosome", "position", "variant_ids", "f_stat", "beta", "se",
                     "rss", "af", "var_perc", "effect_allele", "other_allele"):
            value = getattr(self, name)
            if value is not None:
                value = np.asarray(value)
                if value.shape != self.p.shape:
                    raise ValueError(f"{name} must have one entry per p-value")
                setattr(self, name, value)

    @classmethod
    def from_scan(cls, scan: dict, gt, fit=None) -> "GwasResult":
        """Build from ``LMM.scan`` output and a Genotypes container."""
        res = cls(
            chromosome=np.asarray(gt.chromosome),
            position=np.asarray(gt.position),
            p=scan["ps"],
            f_stat=scan.get("f_stats"),
            rss=scan.get("rss"),
            var_perc=scan.get("var_perc"),
            beta=scan.get("betas"),
            se=scan.get("ses"),
            variant_ids=np.asarray(gt.variant_ids),
            af=gt.allele_freqs(),
            effect_allele=gt.allele1, other_allele=gt.allele2,
        )
        for name in ("n", "h0_rss"):
            if name in scan:
                res.extra[name] = scan[name]
        if "n" in scan:
            res.n = scan["n"]
        if fit is not None:
            res.extra["pseudo_heritability"] = fit.pseudo_heritability
            res.extra["delta"] = fit.delta
            res.extra["ll"] = fit.ll
        return res

    def __len__(self):
        return self.p.size

    # ------------------------------------------------------------------
    # Statistics
    # ------------------------------------------------------------------

    def genomic_control(self) -> float:
        """Genomic inflation factor lambda_GC (Devlin & Roeder 1999).

        The median 1-df chi-square statistic over its null median (0.455),
        with statistics recovered from the p-values so F and chi-square
        scans are treated alike. Above 1 means inflation. (The 2026-10-02
        revision returned median(p)/0.5 here, which runs the other way:
        values below 1 meant inflation.)
        """
        p = np.asarray(self.p, dtype=np.float64)
        p = p[np.isfinite(p)]
        if p.size == 0:
            return float("nan")
        chi2 = stats.chi2.isf(np.clip(p, 1e-300, 1.0), 1)
        return float(np.median(chi2) / _CHI2_1DF_MEDIAN)

    def neg_log10_p(self) -> np.ndarray:
        """-log10 p per variant, p capped at 1 (p = 0 gives inf)."""
        return -np.log10(np.minimum(self.p, 1.0))

    def bonferroni_threshold(self, alpha: float = 0.05) -> float:
        """``alpha`` divided by the number of variants in the result."""
        return alpha / self.p.size

    def top_snps(self, k: int = 10) -> "GwasResult":
        """The ``k`` variants with the smallest p-values, smallest first."""
        idx = np.argsort(self.p)[:k]
        return self.take(idx)

    def take(self, idx) -> "GwasResult":
        """Results for the variants at ``idx`` (indices or a boolean mask);
        per-variant posterior extras are subset too, other extras copied."""
        idx = np.atleast_1d(idx)
        def sel(a):
            return None if a is None else np.asarray(a)[idx]

        extra = dict(self.extra)
        for key in ("ppa", "log_bf", "priors", "prior_variance"):
            if key in extra:
                extra[key] = np.asarray(extra[key])[idx]
        return GwasResult(
            chromosome=sel(self.chromosome),
            position=sel(self.position),
            p=sel(self.p),
            variant_ids=sel(self.variant_ids),
            f_stat=sel(self.f_stat),
            beta=sel(self.beta),
            se=sel(self.se),
            rss=sel(self.rss),
            af=sel(self.af),
            var_perc=sel(self.var_perc),
            extra=extra,
            effect_allele=sel(self.effect_allele), other_allele=sel(self.other_allele),
            n=self.n,
        )

    def power_analysis(self, causal, window: int = 0, alpha=None) -> dict:
        """SNP-level power summary against known causal variants.

        ``discovered`` counts causal variants that are significant (or have
        a significant SNP within ``window`` bp). ``n_significant_noncausal``
        counts significant SNPs that are not themselves causal -- LD tags of
        a causal variant included, so it is NOT a false-positive count; use
        :meth:`locus_summary` for locus-level power and FDR.
        """
        if alpha is None:
            alpha = self.bonferroni_threshold()
        sig = self.p < alpha
        causal = np.atleast_1d(causal)
        discovered = 0
        for c in causal:
            hit = sig[c]
            if not hit and window > 0:
                lo = self.position[c] - window
                hi = self.position[c] + window
                near = (self.chromosome == self.chromosome[c]) & (
                    (self.position >= lo) & (self.position <= hi)
                )
                hit = bool(sig[near].any())
            discovered += int(hit)
        noncausal = sig.copy()
        noncausal[causal] = False
        return {
            "n_causal": int(causal.size),
            "discovered": discovered,
            "power": discovered / max(causal.size, 1),
            "n_significant": int(sig.sum()),
            "n_significant_noncausal": int(noncausal.sum()),
        }

    def locus_summary(self, causal, window: int, alpha=None) -> dict:
        """Locus-level power and false discoveries.

        Significant SNPs are merged into loci by single linkage: same
        chromosome, consecutive significant SNPs at most ``window`` bp
        apart. A locus is true when a causal variant lies within ``window``
        bp of its span; a causal variant is discovered when a significant
        SNP lies within ``window`` bp of it. ``fdr`` is false loci over all
        loci (0 when nothing is significant).
        """
        if alpha is None:
            alpha = self.bonferroni_threshold()
        causal = np.atleast_1d(np.asarray(causal, dtype=np.int64))
        sig = np.nonzero(self.p < alpha)[0]
        chrom = np.asarray(self.chromosome)
        pos = np.asarray(self.position, dtype=np.int64)
        c_chrom, c_pos = chrom[causal], pos[causal]

        loci = []  # (chromosome, start, stop)
        for c in np.unique(chrom[sig]):
            p_c = np.sort(pos[sig[chrom[sig] == c]])
            breaks = np.nonzero(np.diff(p_c) > window)[0]
            starts = np.r_[0, breaks + 1]
            stops = np.r_[breaks, p_c.size - 1]
            loci.extend((c, p_c[a], p_c[b]) for a, b in zip(starts, stops))
        n_true = 0
        for c, lo, hi in loci:
            on = c_chrom == c
            n_true += int(np.any((c_pos[on] >= lo - window) & (c_pos[on] <= hi + window)))
        discovered = 0
        for c, p in zip(c_chrom, c_pos):
            near = sig[(chrom[sig] == c) & (np.abs(pos[sig] - p) <= window)]
            discovered += int(near.size > 0)
        n_loci = len(loci)
        return {
            "n_causal": int(causal.size),
            "discovered": discovered,
            "power": discovered / max(causal.size, 1),
            "n_significant": int(sig.size),
            "n_loci": n_loci,
            "n_true_loci": n_true,
            "n_false_loci": n_loci - n_true,
            "fdr": (n_loci - n_true) / n_loci if n_loci else 0.0,
        }

    # ------------------------------------------------------------------
    # Bayesian priors (v1 snp_priors / PPA machinery)
    # ------------------------------------------------------------------

    def posterior_probabilities(self, priors, use_f=None, *, prior_variance=None) -> "GwasResult":
        """Marginal association probabilities under a normal effect prior.

        The alternative is beta_j ~ N(0, W_j), with W supplied explicitly
        in squared phenotype-units per counted-allele. Using the normal
        approximation beta_hat | beta ~ N(beta, se^2), the Bayes factor is
        the ratio N(beta_hat; 0, se^2 + W) / N(beta_hat; 0, se^2).
        These are single-variant association probabilities, not joint
        fine-mapping probabilities. Untestable variants retain NaN.

        The former RSS likelihood ratio was not a Bayes factor and its
        cross-variant rescaling changed posterior odds. It is no longer
        used; calls must supply ``prior_variance`` and beta/se estimates.
        """
        if self.extra.get("effect_method") == "one-step-logistic":
            raise ValueError("one-step logistic score estimates do not provide the fitted normal "
                             "likelihood required for posterior probabilities")
        if prior_variance is None:
            raise ValueError("prior_variance is required: specify the normal effect-prior variance")
        if use_f is not None:
            raise ValueError("use_f is no longer supported; supply beta, se and prior_variance")
        if self.beta is None or self.se is None:
            raise ValueError("posterior probabilities require beta and se")
        prior = np.broadcast_to(np.asarray(priors, dtype=float), self.p.shape)
        W = np.broadcast_to(np.asarray(prior_variance, dtype=float), self.p.shape)
        if not np.isfinite(prior).all() or np.any((prior < 0) | (prior > 1)):
            raise ValueError("priors must be finite probabilities in [0, 1]")
        if not np.isfinite(W).all() or np.any(W < 0):
            raise ValueError("prior_variance must be finite and non-negative")
        ok = np.isfinite(self.beta) & np.isfinite(self.se) & (self.se > 0)
        log_bf = np.full(self.p.shape, np.nan)
        V = self.se[ok] ** 2
        ratio = W[ok] / V
        log_bf[ok] = 0.5 * (-np.log1p(ratio) + (self.beta[ok] / self.se[ok]) ** 2 * ratio / (1 + ratio))
        with np.errstate(divide="ignore", invalid="ignore"):
            ppa = special.expit(log_bf + np.log(prior) - np.log1p(-prior))
        ppa[ok & (prior == 0)] = 0
        ppa[ok & (prior == 1)] = 1
        self.extra.update(log_bf=log_bf, ppa=ppa, prior_variance=W.copy(),
                          priors=prior.copy(), posterior_method="normal-approximation")
        return self

    # ------------------------------------------------------------------
    # IO
    # ------------------------------------------------------------------

    def to_dataframe(self):
        """The variant columns :meth:`write_csv` writes, as a pandas DataFrame."""
        import pandas as pd

        data = {
            "chromosome": self.chromosome,
            "position": self.position,
            "p": self.p,
        }
        if self.n is not None:
            data["n"] = np.full(self.p.size, self.n)
        for key, arr in [
            ("variant_id", self.variant_ids),
            ("f_stat", self.f_stat),
            ("beta", self.beta),
            ("se", self.se),
            ("af", self.af),
            ("var_perc", self.var_perc),
            ("rss", self.rss),
            ("effect_allele", self.effect_allele),
            ("other_allele", self.other_allele),
        ]:
            if arr is not None:
                if self.extra.get("effect_method") == "one-step-logistic":
                    key = _SCORE_COLUMNS.get(key, key)
                data[key] = arr
        return pd.DataFrame(data)

    def write_csv(self, path: str) -> None:
        """Write variant columns with round-trip precision.

        The sample count ``n`` is written as a constant column (the headers
        spell out what downstream sumstats readers map to their sample-size
        field). Binary score approximations use beta_one_step/se_null_score
        headers so that their interpretation survives export. Other extra
        metadata is not serialized.
        """
        columns = {"chromosome": self.chromosome, "position": self.position, "p": self.p}
        if self.n is not None:
            columns["n"] = np.full(len(self), self.n)
        for header, name in _CSV_OPTIONAL.items():
            arr = getattr(self, name)
            if arr is not None:
                if self.extra.get("effect_method") == "one-step-logistic":
                    header = _SCORE_COLUMNS.get(header, header)
                columns[header] = arr
        with open(path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(columns)
            writer.writerows(zip(*columns.values()))

    @classmethod
    def read_csv(cls, path: str) -> "GwasResult":
        """Read current columns and the legacy chromosome/position/p aliases."""
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            names = reader.fieldnames or []
            rows = list(reader)

        def col(*candidates, dtype=float, required=False):
            name = next((c for c in candidates if c in names), None)
            if name is None:
                if required:
                    raise ValueError(f"CSV is missing {candidates[0]}")
                return None
            return np.array([row[name] for row in rows], dtype=dtype)

        chrom = col("chromosome", "chromosomes", dtype=str, required=True)
        if chrom.size and all(c.lstrip("-").isdigit() for c in chrom):
            chrom = chrom.astype(np.int64)
        n_vals = col("n", "n_eff")
        n = None
        if n_vals is not None and n_vals.size:
            if not np.isfinite(n_vals).all() or not np.all(n_vals == n_vals[0]):
                raise ValueError("n must be a constant finite sample count")
            if float(n_vals[0]) != int(n_vals[0]):
                raise ValueError("n must be a whole sample count")
            n = int(n_vals[0])
        kwargs = {name: col(header, dtype=str if name in
                  ("variant_ids", "effect_allele", "other_allele") else float)
                  for header, name in _CSV_OPTIONAL.items()}
        kwargs["f_stat"] = col("f_stat", "f_stats")
        kwargs["af"] = col("af", "mafs")
        if any(name in names for name in _SCORE_COLUMNS.values()):
            if "beta" in names or "se" in names:
                raise ValueError("CSV mixes fitted effects with one-step score estimates")
            kwargs.update(beta=col("beta_one_step"), se=col("se_null_score"),
                          extra={"effect_method": "one-step-logistic"})
        return cls(chromosome=chrom,
                   position=col("position", "positions", dtype=np.int64, required=True),
                   p=col("p", "ps", "scores", required=True), n=n, **kwargs)


_CSV_OPTIONAL = {
    "variant_id": "variant_ids", "f_stat": "f_stat", "beta": "beta",
    "se": "se", "rss": "rss", "af": "af", "var_perc": "var_perc",
    "effect_allele": "effect_allele", "other_allele": "other_allele",
}

_SCORE_COLUMNS = {"beta": "beta_one_step", "se": "se_null_score"}
