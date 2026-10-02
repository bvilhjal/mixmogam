"""GWAS results container: tables, priors/PPA, genomic control, power."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

__all__ = ["GwasResult"]


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
        )
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
        """Median-based inflation factor lambda_GC."""
        return float(np.median(self.p) / 0.5)

    def neg_log10_p(self) -> np.ndarray:
        return -np.log10(np.minimum(self.p, 1.0))

    def bonferroni_threshold(self, alpha: float = 0.05) -> float:
        return alpha / self.p.size

    def top_snps(self, k: int = 10) -> "GwasResult":
        idx = np.argsort(self.p)[:k]
        return self.take(idx)

    def take(self, idx) -> "GwasResult":
        def sel(a):
            return None if a is None else np.asarray(a)[idx]

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
            extra=dict(self.extra),
        )

    def power_analysis(self, causal, window: int = 0, alpha=None) -> dict:
        """Power/FDR-style summary against known causal variants.

        Returns counts: discovered causal variants (p below threshold, or
        within ``window`` bp of such a SNP) and false positives.
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
        return {
            "n_causal": int(causal.size),
            "discovered": discovered,
            "power": discovered / max(causal.size, 1),
            "n_significant": int(sig.sum()),
            "n_false_positive": int(sig.sum() - discovered),
        }

    # ------------------------------------------------------------------
    # Bayesian priors (v1 snp_priors / PPA machinery)
    # ------------------------------------------------------------------

    def posterior_probabilities(
        self, priors: np.ndarray, use_f: bool = True
    ) -> "GwasResult":
        """Attach Bayes factors / posterior probabilities from SNP priors.

        Follows v1: log BF_j = (n/2) log(RSS0 / RSS_j); posterior odds =
        BF * prior_odds. Requires the scan's ``rss`` and the null RSS in
        ``extra['h0_rss']``.
        """
        if self.rss is None or "h0_rss" not in self.extra:
            raise ValueError("posterior probabilities need rss and extra['h0_rss']")
        h0 = self.extra["h0_rss"]
        log_bf = 0.5 * self.n_samples_log() * np.log(h0 / self.rss)
        bf = np.exp(log_bf - log_bf.max())  # stabilized
        prior_odds = np.asarray(priors) / np.maximum(1 - np.asarray(priors), 1e-12)
        odds = bf * prior_odds
        ppa = odds / (1.0 + odds)
        out = self  # store in extra to keep the dataclass lean
        out.extra["log_bf"] = log_bf
        out.extra["ppa"] = ppa
        return out

    def n_samples_log(self) -> float:
        return float(self.extra.get("n", 0.0))

    # ------------------------------------------------------------------
    # IO
    # ------------------------------------------------------------------

    def to_dataframe(self):
        import pandas as pd

        data = {
            "chromosome": self.chromosome,
            "position": self.position,
            "p": self.p,
        }
        for key, arr in [
            ("variant_id", self.variant_ids),
            ("f_stat", self.f_stat),
            ("beta", self.beta),
            ("se", self.se),
            ("af", self.af),
            ("var_perc", self.var_perc),
        ]:
            if arr is not None:
                data[key] = arr
        return pd.DataFrame(data)

    def write_csv(self, path: str) -> None:
        cols = ["chromosome", "position", "p"]
        arrays = [self.chromosome, self.position, self.p]
        for key, arr in [
            ("f_stat", self.f_stat),
            ("beta", self.beta),
            ("se", self.se),
            ("af", self.af),
            ("var_perc", self.var_perc),
        ]:
            if arr is not None:
                cols.append(key)
                arrays.append(arr)
        with open(path, "w") as fh:
            fh.write(",".join(cols) + "\n")
            for row in zip(*arrays):
                fh.write(",".join(_fmt(v) for v in row) + "\n")

    @classmethod
    def read_csv(cls, path: str) -> "GwasResult":
        data = np.genfromtxt(
            path, delimiter=",", names=True, dtype=None, encoding="utf-8"
        )
        names = data.dtype.names
        # v1 header normalization (chromosome->chromosomes, maf->macs)
        def col(*cands):
            for c in cands:
                if c in names:
                    return data[c]
            return None

        return cls(
            chromosome=col("chromosomes", "chromosome"),
            position=col("position", "positions"),
            p=col("p", "ps", "scores"),
            f_stat=col("f_stat", "f_stats"),
            af=col("af", "mafs"),
            var_perc=col("var_perc"),
        )


def _fmt(v) -> str:
    if isinstance(v, (str, np.str_)):
        return str(v)
    if isinstance(v, (float, np.floating)):
        return f"{v:.6g}"
    return str(v)
