"""HRATT in two steps (mixmogam.steps) and the scalability pieces under it:
fused score numerators, mapped sample filtering, block results."""

import gzip
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from mixmogam.association import gwas
from mixmogam.cli import main
from mixmogam.genotypes import Dosages, Genotypes
from mixmogam.io.plink import read_plink, write_plink
from mixmogam._loco import LocoGenotypes, loco_groups
from mixmogam._packed import PackedCalls
from mixmogam.results import GwasResult
from mixmogam.simulate import simulate_genotypes, simulate_traits
from mixmogam.steps import HrattStep1, hratt_step1, hratt_step2
from mixmogam.twostep import hratt

from .conftest import require_numba_threads

N, M = 360, 1200


@pytest.fixture(scope="module")
def data():
    G = simulate_genotypes(n=N, m=M, seed=13).T.copy(order="F")  # (n, m)
    ids = np.array([f"s{i}" for i in range(N)])
    gt = Genotypes(G, sample_ids=ids, chromosome=np.repeat([1, 2, 3, 4], M // 4),
                   position=np.tile(np.arange(M // 4) * 100, 4),
                   variant_ids=np.array([f"v{j}" for j in range(M)], dtype=object),
                   allele1=np.array(["A"] * M), allele2=np.array(["C"] * M))
    y = simulate_traits(G.T, h2=0.5, n_causal=15, seed=14)["y"]
    X = np.random.default_rng(15).standard_normal((N, 2))
    return G, gt, y, X


def _fields(res):
    return {f: np.asarray(getattr(res, f), dtype=float) for f in ("p", "f_stat", "beta", "se")}


def _assert_same(a, b, rtol=0.0):
    for name, x in _fields(a).items():
        y = _fields(b)[name]
        if rtol == 0.0:
            np.testing.assert_array_equal(x, y, err_msg=name)
        else:
            np.testing.assert_allclose(x, y, rtol=rtol, atol=0, err_msg=name)


def _per_chromosome(G, gt, rng, extra=5):
    """The panel as one packed file per chromosome, with the samples
    shuffled and ``extra`` samples step 1 never saw."""
    allG = np.vstack([G, rng.integers(0, 3, size=(extra, G.shape[1])).astype(np.int8)])
    ids = np.r_[gt.sample_ids, [f"x{i}" for i in range(extra)]]
    perm = rng.permutation(ids.size)
    out = []
    for c in np.unique(gt.chromosome):
        cols = np.flatnonzero(gt.chromosome == c)
        out.append(Genotypes(PackedCalls.pack(np.asfortranarray(allG[perm][:, cols])),
                             sample_ids=ids[perm], chromosome=gt.chromosome[cols],
                             position=gt.position[cols], variant_ids=gt.variant_ids[cols],
                             allele1=gt.allele1[cols], allele2=gt.allele2[cols]))
    return out


QUANT = [{}, {"heritability_method": "he", "model_snps": np.arange(0, M, 3)},
         {"relatedness": (np.arange(0, 60, 2), np.arange(1, 61, 2), np.full(30, 0.5))}]


@pytest.mark.parametrize("options", QUANT, ids=["reml", "he-subset", "kappa"])
def test_split_steps_equal_hratt_exactly(data, options, tmp_path):
    G, gt, y, X = data
    ref = hratt(y, gt, X, **options)
    step1 = hratt_step1(y, gt, X, **options)
    _assert_same(hratt_step2(step1, gt), ref)  # the step-1 panel itself
    step1.save(tmp_path / "s1.npz")
    loaded = HrattStep1.load(tmp_path / "s1.npz")
    assert np.array_equal(loaded.sample_ids, step1.sample_ids)
    parts = _per_chromosome(G, gt, np.random.default_rng(1))
    # blocks of ~37 variants, samples reordered and selected per block
    streamed = hratt_step2(loaded, parts, chunk_bytes=37 * ((N + 5 + 3) // 4))
    _assert_same(streamed, ref)
    assert np.array_equal(streamed.variant_ids, gt.variant_ids)
    assert streamed.extra["lambda"] == ref.extra["lambda"]
    assert streamed.extra["h2"] == ref.extra["h2"]


@pytest.mark.parametrize("options", [{"trait": "binary"}, {"sample_weights": "w"},
                                     {"relatedness_correction": "gls"}],
                         ids=["binary", "weighted", "gls"])
def test_split_steps_other_tests_match_hratt(data, options):
    G, gt, y, X = data
    options = dict(options)
    if options.get("trait") == "binary":
        y = (y > np.quantile(y, 0.7)).astype(float)
    if "sample_weights" in options:
        options["sample_weights"] = np.random.default_rng(2).uniform(0.5, 2.0, N)
    if "relatedness_correction" in options:
        options["relatedness"] = (np.arange(0, 60, 2), np.arange(1, 61, 2), np.full(30, 0.5))
    ref = hratt(y, gt, X, **options)
    step1 = hratt_step1(y, gt, X, **options)
    _assert_same(hratt_step2(step1, gt), ref)
    # Products over other block shapes round apart in the last bits: the
    # measured worst is 9.4e-13 relative (gls f_stat, M2 Pro, 2026-11-14);
    # 3e-12 is the CHANGELOG 2.1.0 bound, defended here (review F8).
    streamed = hratt_step2(step1, _per_chromosome(G, gt, np.random.default_rng(3)),
                           chunk_bytes=50 * ((N + 5 + 3) // 4))
    _assert_same(streamed, ref, rtol=3e-12)


def test_step2_reads_a_mapped_bed_file_without_copying_the_panel(data, tmp_path):
    G, gt, y, X = data
    write_plink(gt, str(tmp_path / "panel"))
    mapped = read_plink(str(tmp_path / "panel"), mmap=True)
    step1 = hratt_step1(y, gt, X)
    _assert_same(hratt_step2(step1, mapped, chunk_bytes=100 * ((N + 3) // 4)), hratt(y, gt, X))


@pytest.mark.parametrize("format", ["vcf", "pgen"])
def test_cli_step2_streams_input_before_testing(data, tmp_path, monkeypatch, format):
    import mixmogam.io as io

    pg = pytest.importorskip("pgenlib") if format == "pgen" else None
    G, gt, y, X = data
    step1 = hratt_step1(y, gt, X, heritability_method="he")
    ref = hratt_step2(step1, gt)
    fit = tmp_path / "fit.npz"
    step1.save(fit)
    path = tmp_path / f"calls.{format}"
    if pg is None:
        with path.open("w") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                     + "\t".join(gt.sample_ids) + "\n")
            codes = np.array(["0/0", "0/1", "1/1", "./."])
            for j in range(M):
                fh.write(f"{gt.chromosome[j]}\t{gt.position[j]}\t{gt.variant_ids[j]}\tC\tA\t.\t.\t.\tGT\t"
                         + "\t".join(codes[G[:, j]]) + "\n")
    else:
        path.with_suffix(".psam").write_text("#IID\n" + "\n".join(gt.sample_ids) + "\n")
        path.with_suffix(".pvar").write_text("#CHROM\tPOS\tID\tREF\tALT\n" + "".join(
            f"{gt.chromosome[j]}\t{gt.position[j]}\t{gt.variant_ids[j]}\tC\tA\n" for j in range(M)))
        with pg.PgenWriter(bytes(path), N, variant_ct=M, nonref_flags=False) as writer:
            for j in range(M):
                calls = G[:, j].copy()
                calls[calls < 0] = -9
                writer.append_biallelic(calls)
    sizes = []
    if format == "pgen":
        import mixmogam._pgen_score as scorer
        source, name = scorer, "pgen_blocks"
    else:
        source, name = io, "iter_genotypes"
    read = getattr(source, name)

    def bounded(*args, **kwargs):
        for block in read(*args, **kwargs):
            count = len(block[0]) if format == "pgen" else block.n_variants
            sizes.append(count)
            assert count <= 37
            yield block

    def forbidden(*args, **kwargs):
        raise AssertionError("step 2 loaded a whole file")

    monkeypatch.setattr(source, name, bounded)
    monkeypatch.setattr(io, "read_genotypes", forbidden)
    out = tmp_path / "result.tsv"
    assert main(["step2", "--step1", str(fit), "--geno", str(path), "--packed",
                 "--chunk-mb", str(37 * N / 1024**2), "--out", str(out)]) == 0
    assert sum(sizes) == M and len(sizes) > 1
    got = GwasResult.read_csv(out)
    _assert_same(got, ref, rtol=1e-12)
    np.testing.assert_array_equal(got.variant_ids, ref.variant_ids)
    np.testing.assert_array_equal(got.effect_allele, ref.effect_allele)
    np.testing.assert_array_equal(got.other_allele, ref.other_allele)
    assert "gram_cache" in HrattStep1.load(fit).extra


def test_cli_step2_dosage_wiring(data, tmp_path):
    """--dosage streams BGEN dosages through the CLI and its validation is
    CLI-only wiring: the underlying paths are tested elsewhere (review F6)."""
    from .bgen_fixture import write_bgen

    G, gt, y, X = data
    k = 80
    step1 = hratt_step1(y, gt, X, heritability_method="he")
    fit = tmp_path / "fit.npz"
    step1.save(fit)
    D = G[:, :k].astype(np.float64)
    D[D < 0] = np.nan
    dos = Dosages(D, sample_ids=gt.sample_ids, chromosome=gt.chromosome[:k],
                  position=gt.position[:k], variant_ids=gt.variant_ids[:k],
                  allele1=gt.allele1[:k], allele2=gt.allele2[:k])
    ref = hratt_step2(step1, dos)

    records = []
    for j in range(k):
        calls = G[:, j]
        probs = 255 * np.column_stack([calls == 0, calls == 1])
        flags = np.where(calls < 0, 130, 2).astype(np.uint8)
        records.append((gt.variant_ids[j], gt.chromosome[j], gt.position[j], gt.allele2[j],
                        gt.allele1[j], 8, probs, flags))
    path = tmp_path / "d.bgen"
    write_bgen(path, gt.sample_ids, records)
    out = tmp_path / "dos.tsv"
    assert main(["step2", "--step1", str(fit), "--geno", str(path), "--dosage",
                 "--chunk-mb", str(20 * N * 8 / 1024**2), "--out", str(out)]) == 0
    got = GwasResult.read_csv(out)
    np.testing.assert_allclose(got.p, ref.p, rtol=1e-9)
    np.testing.assert_array_equal(got.variant_ids, ref.variant_ids)
    assert main(["step2", "--step1", str(fit), "--geno", str(path), "--dosage",
                 "--geno-format", "vcf", "--out", str(out)]) == 2
    assert main(["step2", "--step1", str(fit), "--geno", str(path), "--dosage",
                 "--hard-call-threshold", "0.9", "--out", str(out)]) == 2


def test_several_traits_with_their_own_samples_share_the_blocks(data, tmp_path):
    G, gt, y, X = data
    keep = np.flatnonzero(np.arange(N) % 9 != 0)  # trait 2 misses every ninth sample
    y2 = y[keep] + np.random.default_rng(4).standard_normal(keep.size)
    gt2 = gt.filter_samples(keep)
    fits = [hratt_step1(y, gt, X), hratt_step1(y2, gt2, X[keep]),
            hratt_step1(y * 2.0 + 1.0, gt, X)]
    parts = _per_chromosome(G, gt, np.random.default_rng(5))
    results = hratt_step2(fits, parts, chunk_bytes=80 * ((N + 8) // 4))
    _assert_same(results[0], hratt(y, gt, X))
    _assert_same(results[1], hratt(y2, gt2, X[keep]))
    _assert_same(results[2], hratt(y * 2.0 + 1.0, gt, X))
    paths = [tmp_path / f"t{k}.tsv.gz" for k in range(3)]
    summary = hratt_step2(fits, parts, chunk_bytes=80 * ((N + 8) // 4), out=paths)
    assert summary["variants"] == M
    for path, res in zip(paths, results):
        back = GwasResult.read_csv(path)
        np.testing.assert_array_equal(back.p, res.p)
        np.testing.assert_array_equal(back.variant_ids, res.variant_ids)


def test_integer_dosages_give_the_hard_call_statistics(data):
    G, gt, y, X = data
    step1 = hratt_step1(y, gt, X, relatedness=(np.arange(0, 60, 2), np.arange(1, 61, 2),
                                               np.full(30, 0.5)))
    D = G.astype(np.float64)
    D[D < 0] = np.nan
    dos = Dosages(D, sample_ids=gt.sample_ids, chromosome=gt.chromosome, position=gt.position,
                  variant_ids=gt.variant_ids, allele1=gt.allele1, allele2=gt.allele2)
    got = hratt_step2(step1, dos, chunk_bytes=8 * N * 70)
    ref = hratt_step2(step1, gt)
    # The calls are tested at float32 precision, dosages in float64: the
    # measured worst is 4.5e-6 relative (M2 Pro, 2026-11-14); 1e-5 is the
    # CHANGELOG 2.1.0 "about 1e-5" bound, defended here (review F8).
    np.testing.assert_allclose(np.log10(got.p), np.log10(ref.p), rtol=1e-5, atol=1e-6)
    np.testing.assert_allclose(got.beta, ref.beta, rtol=1e-5, atol=1e-7)
    np.testing.assert_allclose(got.af, ref.af, rtol=1e-12)
    assert got.extra["dosage"] and np.isnan(got.extra["carrier_families"]).all()
    binary = hratt_step1((y > 0).astype(float), gt, X, trait="binary")
    with pytest.raises(ValueError, match="unweighted quantitative"):
        hratt_step2(binary, dos)


def test_interleaved_dosage_tests_share_norms_and_keep_trait_specific_scores(data, monkeypatch):
    from mixmogam import steps

    _, gt, y, X = data
    rng = np.random.default_rng(41)
    pairs = (np.arange(0, 60, 2), np.arange(1, 61, 2), np.full(30, 0.5))
    fits = [hratt_step1(y, gt, X), hratt_step1(y + rng.normal(size=N), gt, X, relatedness=pairs)]
    picks = (np.arange(16)[:, None] + np.arange(4) * (M // 4)).ravel()
    rows = rng.permutation(N)
    D = rng.uniform(0.0, 2.0, size=(N, picks.size))
    D[::5, ::3] = np.nan
    D[:, 2], D[:, 5] = np.nan, 1.0
    dos = Dosages(D[rows], gt.sample_ids[rows], gt.chromosome[picks], variant_ids=gt.variant_ids[picks])

    def dense_scores(Z, Wp, groups, bounds):
        return np.einsum("ij,ij->j", Z, Wp[:, groups])

    with monkeypatch.context() as mp:
        mp.setattr(steps, "_dosage_scores", dense_scores)
        reference = [hratt_step2(fit, dos) for fit in fits]

    norms = steps._dosage_norms
    blocks = []

    def count_norms(Z, Q):
        blocks.append(Z.shape)
        return norms(Z, Q)

    monkeypatch.setattr(steps, "_dosage_norms", count_norms)
    got = hratt_step2(fits, dos)
    assert blocks == [D.shape]
    np.testing.assert_array_equal(dos.values, D[rows])
    for observed, expected in zip(got, reference):
        for field in ("p", "f_stat", "beta", "se", "af"):
            np.testing.assert_allclose(getattr(observed, field), getattr(expected, field),
                                       rtol=1e-12, atol=1e-14, err_msg=field)
        np.testing.assert_array_equal(observed.variant_ids, dos.variant_ids)
        assert np.isnan(observed.p[[2, 5]]).all()


def test_step2_rejects_unknown_chromosomes_missing_samples_and_spectral(data):
    G, gt, y, X = data
    step1 = hratt_step1(y, gt, X)
    other = Genotypes(G[:, :10], sample_ids=gt.sample_ids, chromosome=np.full(10, 9))
    with pytest.raises(ValueError, match=r"chromosomes \['9'\]"):
        hratt_step2(step1, other)
    fewer = gt.filter_samples(np.arange(N - 1))
    with pytest.raises(KeyError, match="1 step-1 samples are absent"):
        hratt_step2(step1, fewer)
    with pytest.raises(ValueError, match="spectral"):
        hratt_step1(y, gt, X, denominator="spectral")


def test_step1_file_holds_no_pickles(data, tmp_path):
    G, gt, y, X = data
    step1 = hratt_step1(y, gt, X, model_snps=np.arange(0, M, 2))
    step1.save(tmp_path / "s1.npz")
    with np.load(tmp_path / "s1.npz", allow_pickle=False) as f:
        assert all(f[k].dtype != object for k in f.files)
    loaded = HrattStep1.load(tmp_path / "s1.npz")
    assert loaded.extra["model_snps"]["n"] == step1.extra["model_snps"]["n"]
    assert loaded.loco_groups == step1.loco_groups


def test_cli_save_step1_and_step2_match_gwas(data, tmp_path, monkeypatch):
    G, gt, y, X = data
    write_plink(gt, str(tmp_path / "panel"))
    for c in (1, 2, 3, 4):
        write_plink(gt.variant_mask(np.flatnonzero(gt.chromosome == c)), str(tmp_path / f"chr{c}"))
    with open(tmp_path / "pheno.txt", "w") as fh:
        fh.write("IID t1 t2\n")
        for i in range(N):
            fh.write(f"{gt.sample_ids[i]} {'NA' if i % 11 == 0 else y[i]} {y[i] * 0.5 + i % 3}\n")
    common = ["--geno", str(tmp_path / "panel"), "--mmap", "--pheno", str(tmp_path / "pheno.txt")]
    assert main(["gwas", *common, "--pheno-name", "t1", "--method", "hratt",
                 "--out", str(tmp_path / "ref.tsv")]) == 0
    writes = []
    write = GwasResult.write_csv
    def record(self, *args, **kwargs):
        writes.append(len(self))
        return write(self, *args, **kwargs)
    monkeypatch.setattr(GwasResult, "write_csv", record)
    assert main(["gwas", *common, "--pheno-name", "t1", "--save-step1", str(tmp_path / "t1.npz"),
                 "--out", str(tmp_path / "both.tsv"), "--chunk-mb", "0.02"]) == 0
    assert max(writes) < M and sum(writes) == M
    assert (tmp_path / "ref.tsv").read_bytes() == (tmp_path / "both.tsv").read_bytes()
    assert main(["gwas", *common, "--pheno-name", "t2", "--save-step1",
                 str(tmp_path / "t2.npz")]) == 0
    assert main(["step2", "--step1", str(tmp_path / "t1.npz"), str(tmp_path / "t2.npz"),
                 "--geno", *[str(tmp_path / f"chr{c}") for c in (1, 2, 3, 4)], "--mmap",
                 "--out", str(tmp_path / "s1.tsv.gz"), str(tmp_path / "s2.tsv.gz"),
                 "--chunk-mb", "0.02"]) == 0
    ref, got = GwasResult.read_csv(tmp_path / "ref.tsv"), GwasResult.read_csv(tmp_path / "s1.tsv.gz")
    np.testing.assert_array_equal(got.p, ref.p)
    np.testing.assert_array_equal(got.variant_ids, ref.variant_ids)
    assert len(GwasResult.read_csv(tmp_path / "s2.tsv.gz")) == M
    assert main(["step2", "--step1", str(tmp_path / "t1.npz"), "--geno", str(tmp_path / "panel"),
                 "--out", str(tmp_path / "a.tsv"), str(tmp_path / "b.tsv")]) == 2


# ----------------------------------------------------------------------
# Building blocks


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("scaled", [False, True])
@pytest.mark.parametrize("n", [5, 203])
def test_score_numerators_equal_the_decoded_products(packed, scaled, n):
    rng = np.random.default_rng(n)
    m = 300
    G = rng.binomial(2, rng.uniform(0.02, 0.5, m), size=(n, m)).astype(np.int8)
    G[rng.random((n, m)) < 0.05] = -1
    chrom = np.repeat([1, 2, 3], m // 3)
    gt = Genotypes(PackedCalls.pack(G) if packed else np.asfortranarray(G), chromosome=chrom)
    Q = np.linalg.qr(np.column_stack([np.ones(n), rng.standard_normal(n)]))[0]
    lg = LocoGenotypes(gt, loco_groups(chrom)[0], Q, block=37,
                       row_scale=rng.uniform(0.3, 2.0, n) if scaled else None)
    W = rng.standard_normal((n, lg.n_groups))
    want = np.zeros(m)
    for idx, g, Z in lg.wide_slices(16):
        want[idx] = Z @ W[:, g]
    np.testing.assert_allclose(lg.score_numerators(W), want, rtol=1e-13,
                               atol=1e-13 * np.abs(want).max())


def test_score_numerators_are_thread_invariant():
    require_numba_threads(4)
    rng = np.random.default_rng(8)
    G = rng.binomial(2, 0.3, size=(501, 400)).astype(np.int8)
    gt = Genotypes(PackedCalls.pack(G), chromosome=np.repeat([1, 2], 200))
    W = rng.standard_normal((501, 2))
    one = LocoGenotypes(gt, loco_groups(gt.chromosome)[0], n_threads=1).score_numerators(W)
    four = LocoGenotypes(gt, loco_groups(gt.chromosome)[0], n_threads=4).score_numerators(W)
    np.testing.assert_array_equal(one, four)


def test_sample_major_calls_decode_instead():
    G = np.random.default_rng(9).binomial(2, 0.3, size=(40, 30)).astype(np.int8)
    gt = Genotypes(np.ascontiguousarray(G), chromosome=np.repeat([1, 2], 15))
    assert LocoGenotypes(gt, loco_groups(gt.chromosome)[0]).score_numerators(
        np.zeros((40, 2))) is None


@pytest.mark.parametrize("n_threads", [1, 2])
def test_filtering_a_mapped_panel_stays_mapped(tmp_path, n_threads):
    if n_threads > 1:
        require_numba_threads(n_threads)
    rng = np.random.default_rng(10)
    G = rng.binomial(2, 0.3, size=(101, 57)).astype(np.int8)
    G[rng.random(G.shape) < 0.1] = -1
    write_plink(Genotypes(G, sample_ids=np.array([f"s{i}" for i in range(101)])),
                str(tmp_path / "p"))
    mapped = read_plink(str(tmp_path / "p"), mmap=True)
    keep = rng.permutation(101)[:77]
    spill = tmp_path / "spill"
    spill.mkdir()
    sub = mapped.filter_samples(keep, n_threads=n_threads, directory=str(spill))
    assert isinstance(sub.G.data, np.memmap)
    assert not os.listdir(spill)  # unlinked once mapped
    np.testing.assert_array_equal(np.asarray(sub.G), G[keep])
    np.testing.assert_array_equal(sub.sample_ids, mapped.sample_ids[keep])
    in_memory = read_plink(str(tmp_path / "p"), packed=True).filter_samples(keep)
    assert not isinstance(in_memory.G.data, np.memmap)
    np.testing.assert_array_equal(in_memory.G.data, sub.G.data)


def test_concat_and_appended_csv_round_trip(tmp_path):
    def part(k, lo):
        return GwasResult(chromosome=np.full(k, 1), position=np.arange(lo, lo + k),
                          p=np.linspace(0.1, 0.9, k), beta=np.arange(k) * 0.1,
                          variant_ids=np.array([f"v{lo + i}" for i in range(k)]), n=10,
                          extra={"n_spa": 2, "carrier_families": np.arange(k, dtype=float)})
    parts = [part(3, 0), part(4, 3)]
    parts[1].extra.pop("carrier_families")
    whole = GwasResult.concat(parts)
    assert len(whole) == 7 and whole.extra["n_spa"] == 4
    assert np.isnan(whole.extra["carrier_families"][3:]).all()
    path = tmp_path / "r.tsv.gz"
    parts[0].write_csv(path)
    part(4, 3).write_csv(path, append=True)
    with gzip.open(path, "rt") as fh:
        assert sum(line.startswith("chromosome") for line in fh) == 1
    back = GwasResult.read_csv(path)
    np.testing.assert_array_equal(back.position, np.arange(7))
    with pytest.raises(ValueError, match="share n"):
        GwasResult.concat([part(2, 0), GwasResult(chromosome=[1], position=[1], p=[0.5], n=3)])


def test_bgen_dosages_are_expected_counts(tmp_path, monkeypatch):
    import struct

    import mixmogam.io.bgen as module

    path = tmp_path / "d.bgen"
    path.write_bytes(struct.pack("<IIII4sI", 20, 20, 3, 4, b"bgen", 9 | (1 << 31)))
    probs = np.array([[[1, 0, 0], [0, 0.5, 0.5], [0.2, 0.3, 0.5]],
                      [[0, 1, 0], [0, 0, 1], [0.5, 0.5, 0]],
                      [[np.nan] * 3, [1, 0, 0], [0, 0, 1]],
                      [[0, 0, 1], [0.1, 0.8, 0.1], [1, 0, 0]]])
    missing = np.zeros((4, 3), dtype=bool)
    missing[2, 0] = True
    reads = []

    class Reader:
        nsamples, nvariants = 4, 3
        samples = np.array(["a", "b", "c", "d"])
        nalleles = np.full(3, 2)
        ncombinations = np.full(3, 3)
        phased = np.zeros(3, dtype=bool)
        allele_ids = np.array(["A,G"] * 3)
        ids = np.array(["x", "y", "z"])
        chromosomes = np.array(["1", "1", "2"])
        positions = np.array([10, 20, 30])

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def read(self, index, **kwargs):
            reads.append(index[1])
            return (probs[index].copy(), missing[index].copy(),
                    np.full(probs[index].shape[:2], 2, dtype=np.int32))

    monkeypatch.setitem(sys.modules, "bgen_reader",
                        SimpleNamespace(open_bgen=lambda *a, **k: Reader()))
    blocks = list(module.iter_bgen_dosages(path, chunk_variants=2, backend="reader"))
    assert [b.n_variants for b in blocks] == [2, 1] and len(reads) == 2
    values = np.hstack([b.values for b in blocks])
    want = probs[:, :, 1] + 2 * probs[:, :, 2]
    want[2, 0] = np.nan
    np.testing.assert_allclose(values, want)
    assert list(blocks[1].chromosome) == ["2"] and list(blocks[0].allele1) == ["G", "G"]


def test_dosages_validate_their_range():
    with pytest.raises(ValueError, match=r"\[0, 2\]"):
        Dosages(np.array([[0.0, 2.5]]), sample_ids=["a"], chromosome=[1, 1])
    with pytest.raises(ValueError, match="float"):
        Dosages(np.array([[0, 1]]), sample_ids=["a"], chromosome=[1, 1])
    assert Dosages(np.array([[np.nan, 1.0]]), sample_ids=["a"], chromosome=[1, 1]).n_variants == 2


def test_gwas_and_split_api_agree_through_gwas(data):
    """gwas(method='hratt') stays hratt(): the split is internal."""
    G, gt, y, X = data
    _assert_same(gwas(y, gt, X, method="hratt", block=4096), hratt(y, gt, X))


def test_pgen_sparse_and_dense_scores_share_order_and_rank_decisions(data, tmp_path, monkeypatch):
    pg = pytest.importorskip("pgenlib")
    pytest.importorskip("numba")
    from mixmogam.io.pgen import iter_pgen
    import mixmogam._pgen_score as scorer

    _, gt, y, X = data
    rng = np.random.default_rng(784)
    n = gt.n_samples
    rare = np.zeros(n, dtype=np.int8)
    rare[:3] = 1
    fitted = rng.permutation(n)[:-15]
    train = gt.filter_samples(fitted)
    design = np.column_stack([X, rare])[fitted]
    opts = dict(heritability_method="he", random_state=17)
    fits = [hratt_step1(t[fitted], train, design, **opts) for t in (y, y + 0.3*X[:, 0])]
    calls = np.column_stack([rare, rng.binomial(2, .4, n), 2-rare, rare, rare, 0*rare,
                             rng.binomial(2, .25, n), rare]).astype(np.int8)
    calls[11:16, [2]] = -9
    calls[20:23, 3] = 2
    calls[15:18, 4] = -9
    calls[:, -1] = -9
    chrom = 1 + np.arange(calls.shape[1]) % 3
    path = tmp_path / "target"
    order = rng.permutation(n)
    path.with_suffix(".psam").write_text("#IID\n" + "\n".join(gt.sample_ids[order]) + "\n")
    path.with_suffix(".pvar").write_text("#CHROM POS ID REF ALT\n" + "".join(
        f"{c} {j+1} v{j} A C\n" for j,c in enumerate(chrom)))
    with pg.PgenWriter(bytes(path.with_suffix(".pgen")), n, variant_ct=calls.shape[1]) as writer:
        for col in calls.T:
            writer.append_biallelic(col[order])
    count = []
    moments = scorer._moments
    def recorded(rows, *args):
        count.append(len(rows))
        return moments(rows, *args)
    monkeypatch.setattr(scorer, "_moments", recorded)
    budget = 5 * n
    for packed in (False, True):
        expected = hratt_step2(fits, iter_pgen(path, packed=packed, chunk_bytes=budget),
                               chunk_bytes=budget)
        observed = hratt_step2(fits, path, packed=packed, chunk_bytes=budget)
        for actual, ref in zip(observed, expected):
            for key in ("p", "beta", "se", "f_stat", "af"):
                np.testing.assert_allclose(getattr(actual, key), getattr(ref, key),
                                           rtol=3e-9, atol=1e-12, equal_nan=True)
            for key in ("variant_ids", "chromosome", "position", "effect_allele", "other_allele"):
                np.testing.assert_array_equal(getattr(actual, key), getattr(ref, key))
    assert count and max(count) < n // 4
    # A weighted test must retain its variance, SPA and missing-call semantics.
    weighted = hratt_step1(y[fitted], train, design, sample_weights=1+np.abs(y[fitted]), **opts)
    before = len(count)
    expected = hratt_step2(weighted, iter_pgen(path, chunk_bytes=budget), chunk_bytes=budget)
    _assert_same(hratt_step2(weighted, path, chunk_bytes=budget), expected)
    assert len(count) == before


def test_streamed_output_bounds_the_original_panel(data, tmp_path, monkeypatch):
    _, gt, y, X = data
    fit = hratt_step1(y, gt, X, heritability_method="he")
    want = hratt_step2(fit, gt)
    writes = []
    write = GwasResult.write_csv
    def record(self, *args, **kwargs):
        writes.append(len(self))
        return write(self, *args, **kwargs)
    monkeypatch.setattr(GwasResult, "write_csv", record)
    out = tmp_path / "out.tsv"
    row = (gt.n_samples + 3)//4 if gt.packed else gt.n_samples
    assert hratt_step2(fit, gt, chunk_bytes=17*row, out=out)["variants"] == len(want)
    assert max(writes) == 17 and sum(writes) == len(want)
    _assert_same(GwasResult.read_csv(out), want, rtol=1e-11)
