"""Streamed stage-1 subsets: an on-disk two-bit cache of the model SNPs.

HRATT's stage 1 on a variant subset (``model_snps``) normally copies the
subset's calls into RAM, n k / 4 bytes: about 180 GB for 1.44M HapMap3+
variants at n = 5e5. With ``model_snps_stream`` the rows are instead written
once (one extraction pass that reads only the subset's rows of the panel) to
a page-aligned file in PLINK BED row layout and memory-mapped read-only as
:class:`~mixmogam._packed.PackedCalls`. Every compiled kernel then reads the
same bytes it would read from the in-RAM copy, so all results are
bit-identical; each stage-1 pass is one forward sequential scan of the file.

:class:`Readahead` keeps those scans I/O-overlapped and resident memory
bounded: before a slice is decoded the kernel is asked to read the next
window (``MADV_WILLNEED``, asynchronous), and the mapping of slices already
passed is dropped (``MADV_DONTNEED``; clean file pages, re-read from the page
cache or disk if needed again).

Cache files are keyed by a fingerprint of the panel, the subset's indices
and the samples, and reused only when the header matches and the stored
per-variant means and SDs equal the panel's freshly prepared moments (which
the panel's setup computes anyway), so an edited panel is rebuilt.
"""

from __future__ import annotations

import hashlib
import json
import mmap
import os
import tempfile
import time
from typing import Optional

import numpy as np

from mixmogam._packed import PackedCalls

__all__ = ["Readahead", "stream_options", "open_subset_cache",
           "panel_fingerprint"]

VERSION = 1
MAGIC = b"MIXMOGAM-STAGE1-SUBSET\n"
HEADER_BYTES = 4096
_CHUNK_BYTES = 64 * 1024**2
DEFAULT_WINDOW = 64 * 1024**2
DEFAULT_THRESHOLD = 8 * 1024**3
_PAGE = mmap.PAGESIZE


class Readahead:
    """madvise read-ahead / drop-behind over the rows of a mapped 2-D array.

    ``begin()`` starts a sequential pass; ``at(rows)`` is called with the
    ascending row indices about to be decoded."""

    def __init__(self, mm: mmap.mmap, base: int, row_bytes: int, n_rows: int, window: int):
        self.mm, self.base, self.row_bytes, self.n_rows = mm, base, row_bytes, n_rows
        self.window = max(int(window), _PAGE)
        self.end = base + row_bytes * n_rows
        self.calls = 0
        self._behind = self._ahead = base
        self._last = -1

    @classmethod
    def for_array(cls, data, window: int) -> Optional["Readahead"]:
        """For a C-contiguous 2-D ``numpy.memmap`` (or a view of one), else None."""
        if not isinstance(data, np.memmap) or data.ndim != 2 or not data.flags.c_contiguous:
            return None
        owner = data
        while getattr(owner, "_mmap", None) is None and isinstance(owner.base, np.memmap):
            owner = owner.base
        mm = getattr(owner, "_mmap", None)
        if mm is None or not hasattr(mm, "madvise"):
            return None
        start = np.frombuffer(mm, dtype=np.uint8).ctypes.data
        base = data.ctypes.data - start
        return cls(mm, base, int(data.strides[0]), int(data.shape[0]), window)

    def _advise(self, flag, lo: int, hi: int) -> None:
        lo = max(lo - lo % _PAGE, 0)
        hi = min(hi, len(self.mm))
        if hi > lo:
            try:
                self.mm.madvise(flag, lo, hi - lo)
                self.calls += 1
            except (OSError, ValueError):  # advice only; never fatal
                pass

    def begin(self) -> None:
        self._behind = self._ahead = self.base
        self._last = -1
        self._advise(mmap.MADV_SEQUENTIAL, self.base, self.end)

    def at(self, rows) -> None:
        if len(rows) == 0:
            return
        lo = self.base + int(rows[0]) * self.row_bytes
        hi = self.base + (int(rows[-1]) + 1) * self.row_bytes
        if int(rows[0]) < self._last or self._last < 0:  # a new pass or a jump back
            self._behind = self._ahead = lo - lo % _PAGE
        self._last = int(rows[-1])
        # drop the mapping of full pages strictly behind this slice
        cut = lo - lo % _PAGE
        if cut > self._behind:
            self._advise(mmap.MADV_DONTNEED, self._behind, cut)
            self._behind = cut
        # keep about one window requested beyond the slice, re-advised in
        # half-window steps (the kernel reads it while this slice is decoded)
        if self._ahead < hi:
            self._ahead = lo
        target = min(hi + self.window, self.end)
        if target > self._ahead and (target - self._ahead >= self.window // 2 or self._ahead < hi):
            self._advise(mmap.MADV_WILLNEED, self._ahead, target)
            self._ahead = target


def stream_options(value) -> Optional[dict]:
    """Normalize ``model_snps_stream``: None/False (off), True (on),
    ``"auto"``, or a dict with ``mode`` (on/off/auto), ``cache_dir``,
    ``threshold_bytes``, ``window_bytes``, ``verify``."""
    if value is None or value is False:
        return None
    if value is True:
        value = {"mode": "on"}
    elif isinstance(value, str):
        value = {"mode": value}
    if not isinstance(value, dict):
        raise ValueError("model_snps_stream must be True, False, 'on', 'off', 'auto' or a dict")
    known = {"mode", "cache_dir", "threshold_bytes", "window_bytes", "verify"}
    unknown = set(value) - known
    if unknown:
        raise ValueError(f"unknown model_snps_stream keys: {sorted(unknown)}")
    out = {"mode": value.get("mode", "on"), "cache_dir": value.get("cache_dir"),
           "threshold_bytes": float(value.get("threshold_bytes", DEFAULT_THRESHOLD)),
           "window_bytes": int(value.get("window_bytes", DEFAULT_WINDOW)),
           "verify": bool(value.get("verify", False))}
    if out["mode"] not in ("on", "off", "auto"):
        raise ValueError("model_snps_stream mode must be 'on', 'off' or 'auto'")
    if out["mode"] == "off":
        return None
    if not np.isfinite(out["threshold_bytes"]) or out["threshold_bytes"] < 0:
        raise ValueError("model_snps_stream threshold_bytes must be finite and non-negative")
    if out["window_bytes"] < 0:
        raise ValueError("model_snps_stream window_bytes must be non-negative")
    if out["cache_dir"] is not None:
        out["cache_dir"] = os.fspath(out["cache_dir"])
        if not os.path.isdir(out["cache_dir"]):
            raise ValueError(f"model_snps_stream cache_dir {out['cache_dir']!r} is not a directory")
    return out


def _sha(*parts) -> str:
    h = hashlib.sha256()
    for p in parts:
        h.update(p if isinstance(p, (bytes, bytearray, memoryview)) else str(p).encode())
        h.update(b"\0")
    return h.hexdigest()


def panel_fingerprint(gt) -> dict:
    """Identity of the panel's calls: file, size, mtime and offset for a
    memory-mapped BED; otherwise the shape and a hash of 64 spread variants."""
    G = gt.G
    data = G.data if isinstance(G, PackedCalls) else None
    if isinstance(data, np.memmap) and getattr(data, "filename", None):
        st = os.stat(data.filename)
        return {"kind": "bed-mmap", "path": os.path.realpath(data.filename), "size": st.st_size,
                "mtime_ns": st.st_mtime_ns, "offset": int(data.offset), "shape": list(data.shape)}
    n, m = G.shape
    take = np.unique(np.linspace(0, m - 1, num=min(64, m)).astype(np.int64)) if m else np.empty(0, np.int64)
    rows = (np.ascontiguousarray(G.data[take]) if data is not None
            else np.ascontiguousarray(np.asarray(G[:, take])))
    return {"kind": "array", "packed": data is not None, "shape": [int(n), int(m)],
            "sample_sha256": _sha(rows.tobytes())}


def _moments_sha(mean, sd) -> str:
    return _sha(np.ascontiguousarray(mean, dtype=np.float64).tobytes(),
                np.ascontiguousarray(sd, dtype=np.float64).tobytes())


def open_subset_cache(gt, idx: np.ndarray, mean: np.ndarray, sd: np.ndarray,
                      options: dict) -> tuple:
    """(PackedCalls over a read-only mapping of rows ``idx``, info dict).

    ``mean``/``sd`` are the panel's prepared moments of ``idx`` (validation
    of a reused cache). With ``options["cache_dir"]`` the file is kept and
    reused; without it, it lives in the temp directory and is unlinked once
    mapped."""
    idx = np.ascontiguousarray(idx, dtype=np.int64)
    n = int(gt.n_samples)
    nb = (n + 3) // 4
    k = int(idx.size)
    fp = panel_fingerprint(gt)
    ids = np.asarray(gt.sample_ids).astype(str)
    meta = {"version": VERSION, "n": n, "k": k, "row_bytes": nb, "header_bytes": HEADER_BYTES,
            "panel": fp, "idx_sha256": _sha(idx.tobytes()),
            "samples_sha256": _sha("\n".join(ids.tolist())),
            "moments_sha256": _moments_sha(mean, sd)}
    key = _sha(json.dumps({k_: meta[k_] for k_ in ("version", "n", "k", "panel", "idx_sha256",
                                                    "samples_sha256")}, sort_keys=True))[:32]
    keep = options.get("cache_dir") is not None
    directory = options["cache_dir"] if keep else tempfile.gettempdir()
    path = os.path.join(directory, f"mixmogam-stage1-{key}.bed")
    info = {"path": path, "bytes": k * nb, "kept": keep, "key": key}
    hit = False
    if keep and os.path.exists(path) and os.path.exists(path + ".json"):
        try:
            with open(path + ".json") as fh:
                old = json.load(fh)
            same = all(old.get(f) == meta[f] for f in meta)
            hit = same and os.path.getsize(path) == HEADER_BYTES + k * nb
            if hit and options.get("verify"):
                hit = _file_sha(path) == old.get("data_sha256")
            info["reason"] = "hit" if hit else "stale"
        except (OSError, ValueError):
            info["reason"] = "unreadable"
    else:
        info["reason"] = "new"
    t0 = time.perf_counter()
    if not hit:
        if not keep:
            fd, path = tempfile.mkstemp(prefix="mixmogam-stage1-", suffix=".bed", dir=directory)
            os.close(fd)
            info["path"] = path
        meta["data_sha256"] = _write_rows(gt, idx, nb, path + ".tmp")
        os.replace(path + ".tmp", path)
        if keep:
            with open(path + ".json.tmp", "w") as fh:
                json.dump(meta, fh, indent=1, sort_keys=True)
            os.replace(path + ".json.tmp", path + ".json")
        info["data_sha256"] = meta["data_sha256"]
    info["extract_seconds"] = time.perf_counter() - t0
    info["cache_hit"] = hit
    mm = np.memmap(path, dtype=np.uint8, mode="r", offset=HEADER_BYTES, shape=(k, nb))
    if not keep:
        os.unlink(path)  # the mapping keeps the data until it goes
    return PackedCalls(mm, n), info


def _file_sha(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        fh.seek(HEADER_BYTES)
        for chunk in iter(lambda: fh.read(_CHUNK_BYTES), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_rows(gt, idx: np.ndarray, nb: int, path: str) -> str:
    """One pass over the subset's rows of the panel: header, then rows."""
    h = hashlib.sha256()
    header = bytearray(HEADER_BYTES)
    head = MAGIC + json.dumps({"version": VERSION, "n": int(gt.n_samples), "k": int(idx.size),
                               "row_bytes": nb}).encode()
    header[: len(head)] = head
    G = gt.G
    packed = isinstance(G, PackedCalls)
    source = Readahead.for_array(G.data, _CHUNK_BYTES) if packed else None
    if source is not None:
        source.begin()
    step = max(1, _CHUNK_BYTES // max(nb, 1))
    with open(path, "wb") as fh:
        fh.write(bytes(header))
        for start in range(0, idx.size, step):
            take = idx[start : start + step]
            if source is not None and np.all(np.diff(take) > 0):
                source.at(take)
            rows = (np.ascontiguousarray(G.data[take]) if packed
                    else PackedCalls.pack(np.ascontiguousarray(np.asarray(G[:, take]))).data)
            rows = np.ascontiguousarray(rows, dtype=np.uint8)
            h.update(rows.data)
            fh.write(rows.data)
    if source is not None:
        source.begin()  # restart markers (no state kept)
    return h.hexdigest()
