"""Phenotype container with the v1 transformation suite."""

from __future__ import annotations

from typing import Callable, Dict, Optional, Sequence

import numpy as np
from scipy import stats

__all__ = ["Phenotypes"]


def _anscombe(x):
    return 2.0 * np.sqrt(x + 3.0 / 8.0)


def _arcsin_sqrt(x):
    return np.arcsin(np.sqrt(np.clip(x, 0, 1)))


_TRANSFORMS: Dict[str, tuple[Callable, Callable]] = {
    "identity": (lambda x: x, lambda x: x),
    "log": (np.log, np.exp),
    "sqrt": (np.sqrt, np.square),
    "sqr": (np.square, np.sqrt),
    "exp": (np.exp, np.log),
    "anscombe": (_anscombe, lambda x: np.square(x / 2.0) - 3.0 / 8.0),
    "arcsin_sqrt": (_arcsin_sqrt, lambda x: np.square(np.sin(x))),
}


class Phenotypes:
    """Named traits over a shared sample axis with tracked transformations.

    ``data[pid]`` holds ``{"values": float array, "transformation": str,
    "raw_values": array or None}``. Values may contain NaN (missing);
    replicate-aware averaging is available when replicate ids were given
    at parse time.
    """

    def __init__(self, sample_ids: Sequence[str]):
        self.sample_ids = np.asarray(sample_ids)
        self.data: Dict[str, dict] = {}
        self.replicates: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def add(self, pid: str, values) -> None:
        values = np.asarray(values, dtype=np.float64).ravel()
        if values.size != self.sample_ids.size:
            raise ValueError(
                f"trait {pid!r}: {values.size} values for "
                f"{self.sample_ids.size} samples"
            )
        self.data[pid] = {"values": values, "transformation": "identity",
                          "raw_values": None}

    def __contains__(self, pid: str) -> bool:
        return pid in self.data

    def __len__(self) -> int:
        return len(self.data)

    def pids(self):
        return list(self.data.keys())

    # ------------------------------------------------------------------
    # Access and alignment
    # ------------------------------------------------------------------

    def values(self, pid: str) -> np.ndarray:
        return self.data[pid]["values"]

    def transformation(self, pid: str) -> str:
        return self.data[pid]["transformation"]

    def complete(self, pid: str):
        """(sample_ids, values) restricted to non-missing entries."""
        v = self.values(pid)
        ok = np.isfinite(v)
        return self.sample_ids[ok], v[ok]

    def align(self, sample_ids: Sequence[str], pid: str):
        """Values for ``sample_ids`` (NaN where absent or missing)."""
        order = {s: i for i, s in enumerate(self.sample_ids)}
        out = np.full(len(sample_ids), np.nan)
        for j, s in enumerate(sample_ids):
            if s in order:
                out[j] = self.values(pid)[order[s]]
        return out

    # ------------------------------------------------------------------
    # Transformations
    # ------------------------------------------------------------------

    def transform(self, pid: str, name: str, revert: bool = False) -> None:
        if name not in _TRANSFORMS:
            raise ValueError(f"unknown transformation {name!r}")
        rec = self.data[pid]
        fwd, inv = _TRANSFORMS[name]
        if revert:
            if rec["raw_values"] is None:
                raise ValueError(f"trait {pid!r} has no transformation to revert")
            rec["values"] = rec["raw_values"]
            rec["raw_values"] = None
            rec["transformation"] = "identity"
        else:
            rec["raw_values"] = rec["values"].copy()
            v = rec["values"]
            ok = np.isfinite(v)
            out = v.copy()
            out[ok] = fwd(v[ok])
            rec["values"] = out
            rec["transformation"] = name

    def box_cox(self, pid: str, lam: Optional[float] = None) -> float:
        """Box-Cox with lambda chosen by Shapiro-Wilk normality (v1 style)."""
        rec = self.data[pid]
        v = rec["values"]
        ok = np.isfinite(v) & (v > 0)
        if not ok.all():
            raise ValueError("box_cox requires positive, complete values")
        if lam is None:
            lams = np.linspace(-2, 2, 41)
            best, best_p = 0.0, -np.inf
            for l in lams:
                t = self._bc(v, l)
                p = stats.shapiro(t).pvalue
                if p > best_p:
                    best, best_p = l, p
            lam = float(best)
        rec["raw_values"] = rec["values"].copy()
        rec["values"] = self._bc(v, lam)
        rec["transformation"] = f"box_cox({lam:.3f})"
        return lam

    @staticmethod
    def _bc(x, lam):
        if abs(lam) < 1e-8:
            return np.log(x)
        return (np.power(x, lam) - 1.0) / lam

    def most_normal(self, pid: str) -> str:
        """Apply the most-normalizing transform from the standard set.

        Candidates: identity, log, sqrt, anscombe, arcsin_sqrt (the latter
        two only for non-negative data). Returns the winning name.
        """
        rec = self.data[pid]
        v = rec["values"]
        ok = np.isfinite(v)
        pos = v[ok]
        candidates = ["identity", "log", "sqrt"]
        if (pos >= 0).all():
            candidates += ["anscombe", "arcsin_sqrt"]
        best, best_p = "identity", -np.inf
        for name in candidates:
            fwd, _ = _TRANSFORMS[name]
            try:
                t = fwd(pos)
            except ValueError:
                continue
            if not np.isfinite(t).all():
                continue
            p = stats.shapiro(t).pvalue
            if p > best_p:
                best, best_p = name, p
        if best != "identity":
            self.transform(pid, best)
        return best

    # ------------------------------------------------------------------
    # Replicates / incidence
    # ------------------------------------------------------------------

    def convert_to_averages(self, pid: str) -> None:
        """Average values over replicate groups in place."""
        if self.replicates is None:
            raise ValueError("no replicate ids recorded")
        v = self.values(pid)
        reps = self.replicates
        uniq, inv = np.unique(reps, return_inverse=True)
        sums = np.bincount(inv, weights=np.nan_to_num(v))
        cnts = np.bincount(inv, weights=np.isfinite(v).astype(float))
        with np.errstate(invalid="ignore", divide="ignore"):
            means = np.where(cnts > 0, sums / np.maximum(cnts, 1), np.nan)
        keep = np.ones(self.sample_ids.size, dtype=bool)
        first = {}
        for i, u in enumerate(inv):
            if u not in first:
                first[u] = i
            else:
                keep[i] = False
        rec = self.data[pid]
        rec["values"] = means[inv[keep]]
        self.sample_ids = self.sample_ids[keep]
        self.replicates = None

    def incidence_matrix(self, sample_ids: Sequence[str]) -> np.ndarray:
        """Z (n_obs, n_unique) mapping breeding values to observations."""
        order = {s: i for i, s in enumerate(self.sample_ids)}
        cols = [order[s] for s in sample_ids if s in order]
        Z = np.zeros((len(cols), self.sample_ids.size))
        Z[np.arange(len(cols)), cols] = 1.0
        return Z
