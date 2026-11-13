"""HapMap3 model SNPs for HRATT's stage 1 (``model_snps="hm3"``).

The bundled map (``mixmogam/data/hm3_ldpred2.npz``) is a compact copy of
the 1,054,330 HapMap3 variants of LDpred2 (Privé et al. 2020; the
``map_hm3.rds`` of figshare item 19213299, "European LD reference (with
blocks)", CC BY 4.0, file https://figshare.com/ndownloader/files/36360900):
chromosome (1-22), GRCh37 and GRCh38 positions, the rsID and the two
alleles (a0, a1). LDpred2 already removed strand-ambiguous (A/T, C/G)
SNPs from it. The UK Biobank frequencies, LD scores and block IDs of the
original are not bundled. ``mixmogam hm3 build`` writes the same format
from another table (e.g. LDpred2's HapMap3+ ``map_hm3_plus.rds``) for
``model_snps_options={"map": path}``.

Matching a panel (:func:`match_hm3`) uses rsIDs (variant IDs of the form
``rs<digits>``) and chromosome:position in a build (GRCh37 or GRCh38;
detected from the overlap when not given). Stage 1 only needs to know
*which* variants to model, not their orientation, so allele flips and
strand flips are counted, never corrected; alleles (when the panel has
them) only reject position matches that are a different variant at the
same site (an indel, another allele of a multi-allelic site).
"""

from __future__ import annotations

import re
import warnings
from functools import lru_cache
from typing import Optional

import numpy as np

__all__ = ["HM3_OPTIONS", "load_hm3", "match_hm3", "normalize_chromosome", "build_hm3_map"]

HM3_OPTIONS = {"build": "auto", "map": None, "min_frac": 0.5, "min_snps": 300_000,
               "min_maf": 0.01, "drop_ambiguous": False}

_BASES = {"A": 0, "C": 1, "G": 2, "T": 3}
_NA38 = -1
_CHR_ALIAS = {"X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26}
_BUILD_ALIAS = {"grch37": "GRCh37", "hg19": "GRCh37", "b37": "GRCh37", "37": "GRCh37",
                "grch38": "GRCh38", "hg38": "GRCh38", "b38": "GRCh38", "38": "GRCh38"}
_RS = re.compile(r"^rs(\d+)$", re.IGNORECASE)


def _default_map_path():
    from importlib.resources import files
    return files("mixmogam") / "data" / "hm3_ldpred2.npz"


@lru_cache(maxsize=4)
def _load(path: str) -> dict:
    with np.load(path, allow_pickle=False) as z:
        chrom = z["chr"].astype(np.int16)
        pos37 = np.cumsum(z["dpos37"].astype(np.int64))
        d38 = z["d38"].astype(np.int64)
        pos38 = np.where(d38 == np.iinfo(np.int32).min, _NA38, pos37 + d38)
        al = z["alleles"]
        meta = {k[5:]: str(z[k]) for k in z.files if k.startswith("meta_")}
        return {"chr": chrom, "pos37": pos37, "pos38": pos38,
                "rsid": z["rsid"].astype(np.int64), "a0": (al >> 2).astype(np.int8),
                "a1": (al & 3).astype(np.int8), "meta": meta}


def load_hm3(path=None) -> dict:
    """The HapMap3 map as arrays: ``chr`` (1-22), ``pos37``, ``pos38``
    (-1 where GRCh38 has no position), ``rsid`` (the number after 'rs'),
    allele codes ``a0``/``a1`` (A, C, G, T = 0-3) and its provenance ``meta``."""
    return _load(str(_default_map_path() if path is None else path))


def normalize_chromosome(chrom) -> np.ndarray:
    """Integer chromosome codes: '1'/'chr1'/1 -> 1, X/23 -> 23, Y/24 -> 24,
    XY/25 -> 25, MT/M/26 -> 26; anything else -> -1."""
    arr = np.asarray(chrom)
    if arr.dtype.kind in "iu":
        return arr.astype(np.int64)
    uniq, inv = np.unique(arr.astype(str), return_inverse=True)
    codes = np.empty(uniq.size, dtype=np.int64)
    for i, c in enumerate(uniq):
        s = c.strip().upper()
        if s.startswith("CHR"):
            s = s[3:]
        if s.isdigit():
            codes[i] = int(s)
        else:
            codes[i] = _CHR_ALIAS.get(s, -1)
    return codes[inv.reshape(-1)]


def _allele_codes(a) -> np.ndarray:
    """A/C/G/T -> 0-3 (case-insensitive), anything else (indels, 0, .) -> -1."""
    arr = np.asarray(a).astype(str)
    uniq, inv = np.unique(arr, return_inverse=True)
    codes = np.array([_BASES.get(u.strip().upper(), -1) for u in uniq], dtype=np.int8)
    return codes[inv.reshape(-1)]


def _rs_numbers(ids) -> np.ndarray:
    arr = np.asarray(ids).astype(str)
    out = np.full(arr.size, -1, dtype=np.int64)
    for j, s in enumerate(arr):
        m = _RS.match(s.strip())
        if m:
            out[j] = int(m.group(1))
    return out


def _lookup(keys_sorted, order, query):
    """Index into the unsorted table of each query key, or -1."""
    pos = np.searchsorted(keys_sorted, query)
    pos = np.minimum(pos, max(keys_sorted.size - 1, 0))
    hit = (keys_sorted.size > 0) & (keys_sorted[pos] == query) if keys_sorted.size else np.zeros(query.size, bool)
    out = np.full(query.size, -1, dtype=np.int64)
    out[hit] = order[pos[hit]]
    return out


def _site_keys(chrom, pos):
    return chrom.astype(np.int64) * (1 << 32) + pos.astype(np.int64)


def _position_hits(hm, build, chrom, pos, valid):
    col = hm["pos37"] if build == "GRCh37" else hm["pos38"]
    ok = col >= 0
    keys = _site_keys(hm["chr"][ok], col[ok])
    order = np.flatnonzero(ok)
    srt = np.argsort(keys, kind="stable")
    hits = np.full(chrom.size, -1, dtype=np.int64)
    hits[valid] = _lookup(keys[srt], order[srt], _site_keys(chrom[valid], pos[valid]))
    return hits


def match_hm3(gt, *, build: str = "auto", map=None, maf=None, min_maf: float = 0.01,
              drop_ambiguous: bool = False) -> tuple:
    """(indices, report) of the panel's HapMap3 variants.

    A panel variant matches an HM3 variant through its rsID (ID
    ``rs<digits>``, on the same chromosome when the panel's chromosome is
    known), through chromosome:position in ``build`` ("GRCh37", "GRCh38",
    aliases hg19/hg38/37/38; "auto" picks the build whose positions overlap
    the panel more), or both (they must name the same HM3 variant, else the
    variant is dropped). Chromosomes may be '1', 'chr1' or 1; X is X, 23 or
    chrX (the bundled map is autosomal). With panel alleles, a match needs
    the same two bases, in either order (a flip) or complemented (a strand
    flip); a strand-ambiguous pair (A/T, C/G) cannot reveal its strand and
    is kept unless ``drop_ambiguous``. Several panel variants on one HM3
    variant (duplicates, split multi-allelic sites with the same bases)
    keep the first. Variants with MAF at or below ``min_maf`` (``maf``
    over all samples; computed when None and ``min_maf`` > 0) are dropped.
    """
    hm = load_hm3(map)
    m = gt.n_variants
    n_hm = hm["chr"].size
    chrom = normalize_chromosome(gt.chromosome)
    pos = np.asarray(gt.position, dtype=np.int64)
    valid = (chrom > 0) & (pos > 0)
    report = {"n_panel": int(m), "n_hm3": int(n_hm), "map": hm["meta"].get("name", "custom"),
              "panel_chromosomes_unrecognized": int(np.count_nonzero(chrom < 0)),
              "panel_non_autosomal": int(np.count_nonzero((chrom > 22))),}

    # build
    req = str(build)
    overlap = {}
    for b in ("GRCh37", "GRCh38"):
        overlap[b] = int(np.count_nonzero(_position_hits(hm, b, chrom, pos, valid) >= 0))
    report["position_overlap"] = overlap
    if req.lower() == "auto":
        best = max(overlap, key=overlap.get)
        used = best if overlap[best] > 0 else None
        if used is not None and overlap["GRCh37"] == overlap["GRCh38"]:
            used = None  # cannot tell; rsIDs only
    elif req.lower() in ("none", "rsid"):
        used = None
    else:
        used = _BUILD_ALIAS.get(req.lower().replace("_", ""))
        if used is None:
            raise ValueError(f"unknown build {build!r}; use 'auto', 'GRCh37' or 'GRCh38'")
        other = "GRCh38" if used == "GRCh37" else "GRCh37"
        if overlap[other] > 2 * max(overlap[used], 1) and overlap[other] >= 50:
            warnings.warn(f"HM3 matching: the panel's positions overlap HapMap3 in {other} "
                          f"({overlap[other]}) far more than in the requested {used} "
                          f"({overlap[used]}); check the build", RuntimeWarning, stacklevel=3)
    report["build_requested"], report["build"] = req, used

    pos_hit = (_position_hits(hm, used, chrom, pos, valid) if used is not None
               else np.full(m, -1, dtype=np.int64))
    rs = _rs_numbers(gt.variant_ids)
    srt = np.argsort(hm["rsid"], kind="stable")
    rs_hit = np.full(m, -1, dtype=np.int64)
    has_rs = rs >= 0
    rs_hit[has_rs] = _lookup(hm["rsid"][srt], srt, rs[has_rs])
    dropped = {}

    def drop(mask, reason):
        c = int(np.count_nonzero(mask))
        if c:
            dropped[reason] = dropped.get(reason, 0) + c

    # rsID hits on another chromosome than the panel says
    wrong_chr = (rs_hit >= 0) & (chrom > 0) & (hm["chr"][np.maximum(rs_hit, 0)] != chrom)
    drop(wrong_chr & (pos_hit < 0), "rsID is on another chromosome in HM3")
    rs_hit[wrong_chr] = -1
    both = (rs_hit >= 0) & (pos_hit >= 0)
    conflict = both & (rs_hit != pos_hit)
    drop(conflict, "rsID and position name different HM3 variants")
    hit = np.where(rs_hit >= 0, rs_hit, pos_hit)
    hit[conflict] = -1
    by = np.full(m, "", dtype=object)
    report["matched_by"] = {
        "rsid+position": int(np.count_nonzero(both & ~conflict)),
        "rsid only": int(np.count_nonzero((rs_hit >= 0) & (pos_hit < 0))),
        "position only": int(np.count_nonzero((rs_hit < 0) & (pos_hit >= 0))),
    }
    del by
    report["panel_with_rsid"] = int(np.count_nonzero(has_rs))

    # alleles
    cand = np.flatnonzero(hit >= 0)
    h = hit[cand]
    a0, a1 = hm["a0"][h], hm["a1"][h]
    ambiguous = a0 == 3 - a1
    allele_report = {}
    if gt.allele1 is not None and gt.allele2 is not None:
        p1 = _allele_codes(np.asarray(gt.allele1)[cand])
        p2 = _allele_codes(np.asarray(gt.allele2)[cand])
        snp = (p1 >= 0) & (p2 >= 0)
        same = snp & (p1 == a1) & (p2 == a0)
        swap = snp & (p1 == a0) & (p2 == a1)
        c1, c2 = np.where(p1 >= 0, 3 - p1, -1), np.where(p2 >= 0, 3 - p2, -1)
        strand = snp & ~ambiguous & (c1 == a1) & (c2 == a0)
        strand_swap = snp & ~ambiguous & (c1 == a0) & (c2 == a1)
        ok = same | swap | strand | strand_swap
        allele_report = {"same": int(np.count_nonzero(same)), "swapped": int(np.count_nonzero(swap)),
                         "strand flip": int(np.count_nonzero(strand)),
                         "strand flip, swapped": int(np.count_nonzero(strand_swap))}
        drop(~snp, "panel alleles are not two single bases (indel or missing allele)")
        drop(snp & ~ok, "alleles differ from HM3 (another variant at the site)")
    else:
        ok = np.ones(cand.size, dtype=bool)
        allele_report = {"unchecked (panel has no alleles)": int(cand.size)}
    amb_ok = ok & ambiguous
    allele_report["strand-ambiguous (A/T, C/G)"] = int(np.count_nonzero(amb_ok))
    if drop_ambiguous:
        drop(amb_ok, "strand-ambiguous (drop_ambiguous)")
        ok &= ~ambiguous
    report["alleles"] = allele_report
    cand, h = cand[ok], h[ok]

    # duplicates: one panel variant per HM3 variant (the first)
    first = np.unique(h, return_index=True)[1]
    keep = np.zeros(cand.size, dtype=bool)
    keep[first] = True
    drop(~keep, "duplicate (another panel variant already matched this HM3 variant)")
    cand = cand[keep]

    if min_maf and min_maf > 0:
        if maf is None:
            maf = gt.allele_freqs(minor=True)
        maf = np.asarray(maf, dtype=np.float64)
        low = ~(maf[cand] > min_maf)
        drop(low, f"MAF <= {min_maf} in this sample")
        cand = cand[~low]
    report["dropped"] = dropped
    on_panel = np.isin(hm["chr"], np.unique(chrom[chrom > 0]))
    report["n_hm3_on_panel_chromosomes"] = int(np.count_nonzero(on_panel))
    report["n_matched"] = int(cand.size)
    report["frac_hm3"] = float(cand.size / n_hm) if n_hm else 0.0
    denom = report["n_hm3_on_panel_chromosomes"]
    report["frac_hm3_panel_chromosomes"] = float(cand.size / denom) if denom else 0.0
    return np.sort(cand).astype(np.int64), report


def format_report(rep: dict) -> str:
    """Multi-line text of a :func:`match_hm3` report."""
    lines = [f"HM3 model SNPs ({rep['map']}): {rep['n_matched']:,} of {rep['n_hm3']:,} HM3 variants "
             f"({100 * rep['frac_hm3']:.1f}%; {100 * rep['frac_hm3_panel_chromosomes']:.1f}% of the "
             f"{rep['n_hm3_on_panel_chromosomes']:,} on the panel's chromosomes) matched in a panel "
             f"of {rep['n_panel']:,}",
             f"  build: {rep['build'] or 'none (rsIDs only)'} (requested {rep['build_requested']}; "
             f"position overlap GRCh37 {rep['position_overlap']['GRCh37']:,}, "
             f"GRCh38 {rep['position_overlap']['GRCh38']:,})",
             "  matched by: " + ", ".join(f"{k} {v:,}" for k, v in rep["matched_by"].items()),
             "  alleles: " + ", ".join(f"{k} {v:,}" for k, v in rep["alleles"].items())]
    if rep["dropped"]:
        lines.append("  dropped: " + "; ".join(f"{k}: {v:,}" for k, v in rep["dropped"].items()))
    else:
        lines.append("  dropped: none")
    if rep["panel_non_autosomal"]:
        lines.append(f"  panel variants on X/Y/MT (not in the map): {rep['panel_non_autosomal']:,}")
    if rep["panel_chromosomes_unrecognized"]:
        lines.append(f"  panel variants on unrecognized chromosomes: {rep['panel_chromosomes_unrecognized']:,}")
    return "\n".join(lines)


def _read_table(path):
    """chr, pos, pos_hg38, a0, a1, rsid columns of an .rds (needs the
    ``rdata`` package) or a delimited text table with a header."""
    p = str(path)
    if p.endswith(".rds"):
        try:
            import rdata
        except ImportError as exc:  # pragma: no cover - optional
            raise ImportError("reading .rds needs the rdata package (pip install rdata); "
                              "or export the table to TSV in R") from exc
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            df = rdata.read_rds(p)
        return {c: np.asarray(df[c]) for c in df.columns}
    import csv
    import gzip

    opener = gzip.open if p.endswith(".gz") else open
    with opener(p, "rt") as fh:
        lines = [line for line in fh if line.strip() and not line.startswith("#")]
    if not lines:
        raise ValueError(f"no rows in {p}")
    delim = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else None)
    if delim is None:
        rows = [line.split() for line in lines]
    else:
        rows = list(csv.reader(lines, delimiter=delim))
    header, rows = rows[0], rows[1:]
    cols = list(zip(*rows)) if rows else [[] for _ in header]
    return {h: np.asarray(c) for h, c in zip(header, cols)}


def build_hm3_map(table, out, *, name: str = "custom", source: str = "", licence: str = "",
                  extra_meta: Optional[dict] = None) -> dict:
    """Write a map in the bundled format from a table with columns ``chr``,
    ``pos`` (GRCh37), ``pos_hg38`` (optional, NA allowed), ``a0``, ``a1``
    and ``rsid`` (``rs<digits>``); ``table`` is a dict of columns or a path
    (.rds, or delimited text with a header). Rows that are not single-base
    SNPs on chromosomes 1-26 with an rsID are skipped. Returns a summary."""
    t = _read_table(table) if not isinstance(table, dict) else table
    chrom = normalize_chromosome(t["chr"])
    pos37 = np.asarray(t["pos"]).astype(np.float64)
    p38 = t.get("pos_hg38")
    pos38 = (np.full(chrom.size, np.nan) if p38 is None else
             np.array([np.nan if (x is None or str(x) in ("", "NA", "nan", "<NA>")) else float(x)
                       for x in np.asarray(p38, dtype=object)]))
    a0, a1 = _allele_codes(t["a0"]), _allele_codes(t["a1"])
    rs = _rs_numbers(t["rsid"])
    ok = (chrom > 0) & np.isfinite(pos37) & (a0 >= 0) & (a1 >= 0) & (a0 != a1) & (rs >= 0)
    chrom, pos37, pos38, a0, a1, rs = chrom[ok], pos37[ok], pos38[ok], a0[ok], a1[ok], rs[ok]
    order = np.lexsort((pos37, chrom))
    chrom, pos37, pos38, a0, a1, rs = (x[order] for x in (chrom, pos37, pos38, a0, a1, rs))
    p37 = pos37.astype(np.int64)
    d38 = np.where(np.isfinite(pos38), pos38 - pos37, np.iinfo(np.int32).min)
    meta = {"name": name, "source": source, "licence": licence, "n": str(chrom.size),
            **{k: str(v) for k, v in (extra_meta or {}).items()}}
    np.savez_compressed(out, chr=chrom.astype(np.int8),
                        dpos37=np.diff(p37, prepend=0).astype(np.int32),
                        d38=d38.astype(np.int32), rsid=rs.astype(np.uint32),
                        alleles=(a0.astype(np.uint8) * 4 + a1.astype(np.uint8)).astype(np.uint8),
                        **{f"meta_{k}": np.array(v) for k, v in meta.items()})
    _load.cache_clear()
    return {"n": int(chrom.size), "skipped": int(np.count_nonzero(~ok)),
            "no_grch38": int(np.count_nonzero(~np.isfinite(pos38)))}
