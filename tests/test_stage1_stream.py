"""Streamed stage-1 subsets (model_snps_stream=): equivalence, cache, edges."""

import json
import mmap
import os

import numpy as np
import pytest

def _require_numba_threads(n=2):
    numba = pytest.importorskip("numba")
    if numba.config.NUMBA_NUM_THREADS < n:
        pytest.skip(f"requires NUMBA_NUM_THREADS >= {n} before Python starts")
    return n


from mixmogam import _stream
from mixmogam._loco import LocoGenotypes
from mixmogam._packed import PackedCalls
from mixmogam.association import gwas
from mixmogam.cli import main
from mixmogam.genotypes import Genotypes
from mixmogam.io.plink import read_plink, write_plink

from tests.test_stage1_subset import ld_panel

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def _panel(n=401, m=903, seed=11, missing=0.0, packed=True, structure=False):
    G, gt = ld_panel(n, m, n_chrom=3, seed=seed)
    rng = np.random.default_rng(seed)
    if structure:  # two populations with shifted allele frequencies
        G2, _ = ld_panel(n, m, n_chrom=3, seed=seed + 100)
        half = n // 2
        G = np.vstack([G[:half], G2[half:]])
    if missing:
        G = G.copy()
        G[rng.random(G.shape) < missing] = -1
    calls = PackedCalls.pack(np.ascontiguousarray(G)) if packed else np.ascontiguousarray(G)
    gt = Genotypes(calls, chromosome=gt.chromosome, position=gt.position, variant_ids=gt.variant_ids)
    Gf = np.where(G < 0, 0, G).astype(float)
    causal = rng.choice(m, 25, replace=False)
    g = (Gf[:, causal] - Gf[:, causal].mean(0)) @ rng.normal(0, 1, 25)
    y = g / g.std() + rng.standard_normal(n)
    if structure:
        y += 0.8 * (np.arange(n) >= n // 2)
    idx = np.sort(rng.choice(m, m // 4, replace=False))
    return gt, y, idx


def _same(a, b):
    assert np.array_equal(a.p, b.p, equal_nan=True)
    assert np.array_equal(a.beta, b.beta, equal_nan=True)
    assert np.array_equal(a.se, b.se, equal_nan=True)
    for key in ("h2", "lambda", "alpha"):
        assert a.extra.get(key) == b.extra.get(key)


@pytest.mark.parametrize("kw", [{}, {"block": 1}, {"block": 7}, {"block": 100000},
                                {"trait": "binary"}, {"weights": True}, {"n_threads": 2}])
def test_stream_bitwise_equals_in_memory(kw, tmp_path):
    if kw.get("n_threads", 1) and kw.get("n_threads", 1) > 1:
        _require_numba_threads(int(kw["n_threads"]))
    gt, y, idx = _panel(missing=0.02)
    kw = dict(kw)
    if kw.pop("weights", False):
        kw["sample_weights"] = np.random.default_rng(1).uniform(0.5, 2.0, gt.n_samples)
    if kw.get("trait") == "binary":
        y = (y > np.median(y)).astype(float)
    a = gwas(y, gt, method="hratt", model_snps=idx, **kw)
    b = gwas(y, gt, method="hratt", model_snps=idx, **kw,
             model_snps_stream={"mode": "on", "cache_dir": str(tmp_path), "window_bytes": 4096})
    _same(a, b)
    s = b.extra["model_snps"]["stream"]
    assert s["streamed"] and s["bytes"] == idx.size * ((gt.n_samples + 3) // 4)


def test_stream_strong_structure_loco_and_lambda(tmp_path):
    gt, y, idx = _panel(n=400, structure=True)
    a = gwas(y, gt, method="hratt", model_snps=idx, structure_snps=300)
    assert a.extra["structure"]["strong"]
    b = gwas(y, gt, method="hratt", model_snps=idx, structure_snps=300,
             model_snps_stream={"mode": "on", "cache_dir": str(tmp_path)})
    _same(a, b)


def test_int8_panel_and_temp_cache_unlinked(tmp_path, monkeypatch):
    gt, y, idx = _panel(packed=False, missing=0.01)
    monkeypatch.setattr("tempfile.tempdir", str(tmp_path))
    a = gwas(y, gt, method="hratt", model_snps=idx)
    b = gwas(y, gt, method="hratt", model_snps=idx, model_snps_stream=True)
    _same(a, b)
    s = b.extra["model_snps"]["stream"]
    assert not s["kept"] and not os.path.exists(s["path"]) and os.listdir(tmp_path) == []


def test_cache_roundtrip_reuse_and_invalidation(tmp_path):
    gt, y, idx = _panel()
    opt = {"mode": "on", "cache_dir": str(tmp_path), "verify": True}
    a = gwas(y, gt, method="hratt", model_snps=idx, model_snps_stream=opt)
    s = a.extra["model_snps"]["stream"]
    assert s["reason"] == "new" and not s["cache_hit"]
    # bytes on disk are the subset's packed rows after a page-sized header
    raw = np.fromfile(s["path"], dtype=np.uint8)
    nb = (gt.n_samples + 3) // 4
    assert raw.size == _stream.HEADER_BYTES + idx.size * nb
    assert raw[: len(_stream.MAGIC)].tobytes() == _stream.MAGIC
    assert np.array_equal(raw[_stream.HEADER_BYTES:].reshape(idx.size, nb), gt.G.data[idx])
    meta = json.load(open(s["path"] + ".json"))
    assert meta["k"] == idx.size and meta["n"] == gt.n_samples
    b = gwas(y, gt, method="hratt", model_snps=idx, model_snps_stream=opt)
    assert b.extra["model_snps"]["stream"]["cache_hit"]
    _same(a, b)
    # another subset: another key
    c = gwas(y, gt, method="hratt", model_snps=idx[1:], model_snps_stream=opt)
    assert c.extra["model_snps"]["stream"]["path"] != s["path"]
    # A corrupted payload is caught even when the source identity matches.
    with open(s["path"], "r+b") as fh:
        fh.seek(_stream.HEADER_BYTES + 5)
        byte = fh.read(1)
        fh.seek(-1, 1)
        fh.write(bytes([byte[0] ^ 255]))
    e = gwas(y, gt, method="hratt", model_snps=idx, model_snps_stream=opt)
    assert e.extra["model_snps"]["stream"]["reason"] == "stale"
    _same(a, e)
    # an edited panel with the same fingerprint: moments differ -> rebuilt
    G = gt.G.data.copy()
    G[idx[0]] = G[idx[1]]
    gt2 = Genotypes(PackedCalls(G, gt.n_samples), chromosome=gt.chromosome,
                    position=gt.position, variant_ids=gt.variant_ids)
    old = _stream.panel_fingerprint
    try:
        _stream.panel_fingerprint = lambda g, take=None: old(gt, take)
        d = gwas(y, gt2, method="hratt", model_snps=idx, model_snps_stream=opt)
    finally:
        _stream.panel_fingerprint = old
    sd = d.extra["model_snps"]["stream"]
    assert sd["path"] == s["path"] and sd["reason"] == "stale" and not sd["cache_hit"]
    _same(d, gwas(y, gt2, method="hratt", model_snps=idx))


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("verify", [False, True])
def test_cache_detects_permuted_calls_with_identical_moments(tmp_path, monkeypatch, packed, verify):
    rng = np.random.default_rng(430)
    n, m = 120, 200
    G = rng.binomial(2, 0.3, (n, m)).astype(np.int8)
    probes = np.unique(np.linspace(0, m - 1, 64).astype(np.int64))
    idx = np.setdiff1d(np.arange(m), probes)
    changed = G.copy()
    changed[:, idx] = changed[rng.permutation(n)[:, None], idx]
    panels = [Genotypes(PackedCalls.pack(g) if packed else g) for g in (G, changed)]
    Q = np.ones((n, 1)) / np.sqrt(n)
    lgs = [LocoGenotypes(gt, np.zeros(m, int), Q) for gt in panels]
    np.testing.assert_array_equal(lgs[0].mean, lgs[1].mean)
    np.testing.assert_array_equal(lgs[0].sd, lgs[1].sd)
    fp = _stream.panel_fingerprint(panels[1], idx)
    monkeypatch.setattr(_stream, "_CHUNK_BYTES", 1024)  # hash across many chunks
    assert _stream.panel_fingerprint(panels[1], idx) == fp
    opt = {"cache_dir": str(tmp_path), "verify": verify}

    def cached(j):
        return _stream.open_subset_cache(panels[j], idx, lgs[j].mean[idx], lgs[j].sd[idx], opt)

    _, old = cached(0)
    calls, new = cached(1)
    assert new["key"] != old["key"] and not new["cache_hit"]
    np.testing.assert_array_equal(np.asarray(calls[:, :]), changed[:, idx])
    assert cached(1)[1]["cache_hit"]
    # Calls outside the requested subset cannot change its payload.
    if packed:
        panels[1].G.data[probes[0]] ^= 255
    else:
        panels[1].G[:, probes[0]] = 2 - panels[1].G[:, probes[0]]
    assert cached(1)[1]["cache_hit"]


def test_memmapped_bed_panel(tmp_path):
    gt, y, idx = _panel(n=203, missing=0.01)
    write_plink(gt, str(tmp_path / "p"))
    gm = read_plink(str(tmp_path / "p"), mmap=True)
    assert isinstance(gm.G.data, np.memmap)
    fp = _stream.panel_fingerprint(gm)
    assert fp["kind"] == "bed-mmap" and fp["offset"] == 3
    a = gwas(y, gm, method="hratt", model_snps=idx)
    b = gwas(y, gm, method="hratt", model_snps=idx,
             model_snps_stream={"mode": "on", "cache_dir": str(tmp_path), "window_bytes": 1})
    _same(a, b)


def test_fingerprint_hashes_mapping_views_and_private_changes(tmp_path):
    gt, _, _ = _panel(n=40, m=90)
    write_plink(gt, str(tmp_path / "p"))
    gm = read_plink(str(tmp_path / "p"), mmap=True)
    views = [Genotypes(PackedCalls(gm.G.data[s:s + 40], 40)) for s in (0, 40)]
    assert _stream.panel_fingerprint(views[0]) != _stream.panel_fingerprint(views[1])
    private = np.memmap(tmp_path / "p.bed", dtype=np.uint8, mode="c", offset=3, shape=(90, 10))
    gc = Genotypes(PackedCalls(private, 40))
    before = _stream.panel_fingerprint(gc)
    private[1] ^= 255
    assert _stream.panel_fingerprint(gc) != before
    np.testing.assert_array_equal(gm.G.data, gt.G.data)


def test_auto_threshold(tmp_path):
    gt, y, idx = _panel()
    nbytes = gt.n_samples * idx.size / 4
    lo = gwas(y, gt, method="hratt", model_snps=idx,
              model_snps_stream={"mode": "auto", "threshold_bytes": nbytes + 1})
    hi = gwas(y, gt, method="hratt", model_snps=idx,
              model_snps_stream={"mode": "auto", "threshold_bytes": nbytes - 1,
                                 "cache_dir": str(tmp_path)})
    assert not lo.extra["model_snps"]["stream"]["streamed"]
    assert hi.extra["model_snps"]["stream"]["streamed"]
    _same(lo, hi)


def test_options_validation():
    gt, y, idx = _panel(n=101, m=303)
    with pytest.raises(ValueError, match="set model_snps"):
        gwas(y, gt, method="hratt", model_snps_stream=True)
    with pytest.raises(ValueError, match="unknown model_snps_stream keys"):
        gwas(y, gt, method="hratt", model_snps=idx, model_snps_stream={"dir": "x"})
    with pytest.raises(ValueError, match="mode"):
        gwas(y, gt, method="hratt", model_snps=idx, model_snps_stream="sometimes")
    with pytest.raises(ValueError, match="not a directory"):
        gwas(y, gt, method="hratt", model_snps=idx, model_snps_stream={"cache_dir": "/nonexistent/x"})
    with pytest.raises(TypeError, match="method='hratt' only"):
        gwas(y, gt, method="exact", model_snps_stream=True)
    assert _stream.stream_options(False) is None and _stream.stream_options("off") is None


def test_moments_reuse_matches_preparation():
    gt, y, idx = _panel(missing=0.03)
    groups = np.repeat(np.arange(3), gt.n_variants // 3)
    rng = np.random.default_rng(0)
    Q = np.linalg.qr(np.column_stack([np.ones(gt.n_samples), rng.standard_normal(gt.n_samples)]))[0]
    full = LocoGenotypes(gt, groups, Q, block=64)
    sub = gt.variant_mask(idx)
    a = LocoGenotypes(sub, groups[idx], Q, block=64, n_groups=3)
    b = LocoGenotypes(sub, groups[idx], Q, block=64, n_groups=3,
                      moments=(full.mean[idx], full.sd[idx], full._projection[idx], full.zz[idx]))
    for name in ("mean", "sd", "_projection", "zz", "_lut"):
        assert np.array_equal(getattr(a, name), getattr(b, name)), name
    assert a.trace == b.trace
    with pytest.raises(ValueError, match="moments"):
        LocoGenotypes(sub, groups[idx], Q, n_groups=3, moments=(full.mean, full.sd, full._projection, full.zz))


def test_readahead_bounds_and_passes(tmp_path):
    """Drop-behind advice: the mapping behind the current slice is released,
    and a jump back (a new pass) restarts the markers; values unchanged."""
    gt, y, idx = _panel(n=4001, m=303, seed=2)
    calls, info = _stream.open_subset_cache(gt, idx, np.zeros(idx.size), np.zeros(idx.size),
                                            {"cache_dir": str(tmp_path)})
    ra = _stream.Readahead.for_array(calls.data, 8192)
    assert ra is not None and ra.row_bytes == 1001
    assert ra.base == _stream.HEADER_BYTES % mmap.ALLOCATIONGRANULARITY
    ra.begin()
    for s in range(0, idx.size, 9):
        ra.at(np.arange(s, min(idx.size, s + 9)))
    assert ra._behind > 0 and ra.calls > 2
    ra.at(np.arange(0, 9))
    assert ra._behind == ra.base - ra.base % mmap.PAGESIZE
    assert np.array_equal(np.asarray(calls.data), gt.G.data[idx])
    assert _stream.Readahead.for_array(np.zeros((3, 3), np.uint8), 10) is None


def test_cli_stream_flags(tmp_path):
    gt, y, idx = _panel(n=151, m=303)
    write_plink(gt, str(tmp_path / "p"))
    (tmp_path / "y.txt").write_text(
        "sample_id\ty\n" + "".join(f"{s}\t{v:.17g}\n" for s, v in zip(gt.sample_ids, y)))
    (tmp_path / "snps.txt").write_text("\n".join(np.asarray(gt.variant_ids)[idx]) + "\n")
    cache = tmp_path / "cache"
    cache.mkdir()
    base = ["gwas", "--geno", str(tmp_path / "p"), "--pheno", str(tmp_path / "y.txt")]
    sel = ["--model-snps", str(tmp_path / "snps.txt")]
    assert main([*base, *sel, "--out", str(tmp_path / "a.csv")]) == 0
    assert main([*base, *sel, "--model-snps-stream", "on", "--model-snps-cache", str(cache),
                 "--out", str(tmp_path / "b.csv")]) == 0
    assert (tmp_path / "a.csv").read_bytes() == (tmp_path / "b.csv").read_bytes()
    assert any(f.endswith(".bed") for f in os.listdir(cache))
    assert main([*base, *sel, "--model-snps-stream", "auto", "--model-snps-stream-threshold", "1e-9",
                 "--out", str(tmp_path / "c.csv")]) == 0
    assert (tmp_path / "a.csv").read_bytes() == (tmp_path / "c.csv").read_bytes()
    assert main([*base, "--model-snps-stream", "on", "--out", str(tmp_path / "d.csv")]) != 0
    assert main([*base, *sel, "--model-snps-cache", str(cache), "--out", str(tmp_path / "e.csv")]) != 0


@pytest.mark.parametrize("n_threads", [1, 2])
def test_prefetch_thread_and_advice_bitwise(tmp_path, n_threads):
    if n_threads > 1:
        _require_numba_threads(n_threads)
    gt, y, idx = _panel(missing=0.02)
    a = gwas(y, gt, method="hratt", model_snps=idx, n_threads=n_threads)
    for mode in ("thread", "advise"):
        b = gwas(y, gt, method="hratt", model_snps=idx, n_threads=n_threads,
                 model_snps_stream={"mode": "on", "cache_dir": str(tmp_path), "prefetch": mode,
                                    "window_bytes": 8192})
        _same(a, b)
    with pytest.raises(ValueError, match="prefetch"):
        _stream.stream_options({"prefetch": "magic"})


def test_prefetcher_reads_ahead_restarts_and_stops(tmp_path):
    import gc
    import threading
    import time

    path = tmp_path / "f.bin"
    path.write_bytes(bytes(range(256)) * 4096)  # 1 MiB
    fd = os.open(path, os.O_RDONLY)
    p = _stream.Prefetcher(fd, chunk=65536, close=True)
    p.request(0, 4096, 1 << 20)
    for _ in range(200):
        if p.stats["bytes"] >= (1 << 20) - 4096:
            break
        time.sleep(0.01)
    assert p.stats["bytes"] >= (1 << 20) - 4096
    p.request(0, 4096, 65536)          # a new pass: restarts at the consumer
    time.sleep(0.05)
    assert p.stats["restarts"] == 0
    n0 = threading.active_count()
    del p
    gc.collect()
    time.sleep(0.1)
    assert threading.active_count() <= n0 - 1
    with pytest.raises(OSError):
        os.fstat(fd)                   # closed with the prefetcher


def test_readahead_thread_on_temp_cache(tmp_path, monkeypatch):
    gt, y, idx = _panel(n=4001, m=303, seed=2)
    monkeypatch.setattr("tempfile.tempdir", str(tmp_path))
    calls, info = _stream.open_subset_cache(gt, idx, np.zeros(idx.size), np.zeros(idx.size),
                                            {"cache_dir": None})
    assert not os.path.exists(info["path"])            # unlinked, still readable by fd
    ra = _stream.Readahead.for_array(calls.data, 1 << 16, prefetch="thread")
    assert ra.prefetcher is not None
    assert ra.file_offset + ra.base == _stream.HEADER_BYTES
    ra.begin()
    # Check the thread on the first request, which asks for a window beyond
    # the slice. Driving the whole pass first is racy: if the thread is not
    # scheduled before the last slice, the final request has target == pos
    # (end of file), nothing is left to prefetch and bytes stays 0 forever.
    import time
    ra.at(np.arange(0, 9))
    deadline = time.monotonic() + 5
    while ra.prefetcher.stats["bytes"] == 0 and time.monotonic() < deadline:
        time.sleep(0.01)
    assert ra.prefetcher.stats["bytes"] > 0
    for s in range(9, idx.size, 9):
        ra.at(np.arange(s, min(idx.size, s + 9)))
    assert np.array_equal(np.asarray(calls.data), gt.G.data[idx])


def test_subset_cache_concurrent_same_key_writers(tmp_path, monkeypatch):
    """Parallel traits sharing a cache_dir write the same key; neither
    writer may crash on a shared temp file, and the published cache must
    be a valid hit (review F13)."""
    import threading

    gt, y, idx = _panel(n=60, m=42, seed=13, packed=False)
    mean, sd = np.asarray(gt.G).mean(axis=0), np.asarray(gt.G).std(axis=0) + 0.1
    opts = _stream.stream_options({"mode": "on", "cache_dir": str(tmp_path)})
    barrier = threading.Barrier(2, timeout=30)
    real_write = _stream._write_rows

    def synced_write(gt_, idx_, nb, path):
        sha = real_write(gt_, idx_, nb, path)
        barrier.wait()  # both writers finish writing before either publishes
        return sha

    monkeypatch.setattr(_stream, "_write_rows", synced_write)
    errors, results = [], []

    def writer():
        try:
            calls, _ = _stream.open_subset_cache(gt, idx, mean[idx], sd[idx], opts)
            results.append(np.asarray(calls.data))
        except Exception as e:
            errors.append(e)

    threads = [threading.Thread(target=writer) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(30)
    assert not errors, errors
    assert len(results) == 2
    np.testing.assert_array_equal(results[0], results[1])
    calls, info = _stream.open_subset_cache(gt, idx, mean[idx], sd[idx], opts)
    assert info["cache_hit"]
    np.testing.assert_array_equal(np.asarray(calls.data), results[0])
    assert not [f for f in os.listdir(tmp_path) if f.endswith(".tmp")]
