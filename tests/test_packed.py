"""Two-bit genotype storage: a drop-in for int8 calls with identical results."""

import numpy as np
import pytest

from mixmogam import gwas
from mixmogam._loco import LocoGenotypes
from mixmogam._packed import PackedCalls
from mixmogam.genotypes import Genotypes, MISSING
from mixmogam.io.plink import read_plink, write_plink


def _calls(rng, n, m, missing=0.05):
    G = rng.binomial(2, rng.uniform(0.05, 0.6, m), size=(n, m)).astype(np.int8)
    G[rng.random((n, m)) < missing] = MISSING
    return G


@pytest.mark.parametrize("n, m", [(0, 3), (1, 1), (3, 5), (4, 7), (7, 3), (1001, 37)])
def test_round_trips_indexing_and_selection(n, m):
    rng = np.random.default_rng(n + m)
    G = rng.choice(np.array([-1, 0, 1, 2], dtype=np.int8), size=(n, m))
    P = PackedCalls.pack(G)
    assert P.shape == G.shape and P.dtype == np.int8 and P.data.shape == (m, (n + 3) // 4)
    np.testing.assert_array_equal(np.asarray(P), G)
    if n % 4:  # PLINK's unused trailing bits are zero
        assert not np.any(P.data[:, -1] >> (2 * (n % 4)))
    if not n:
        return
    for cols in (1 % m, -1, slice(1, None), [m - 1, 0, m - 1], rng.random(m) < 0.5):
        np.testing.assert_array_equal(P[:, cols], G[:, cols])
    rows = rng.permutation(n)[: max(1, n // 2)]
    np.testing.assert_array_equal(P[rows], G[rows])
    np.testing.assert_array_equal(np.asarray(P.take_samples(rows)), G[rows])
    np.testing.assert_array_equal(np.asarray(P.take_variants([m - 1, 0])), G[:, [m - 1, 0]])
    drop = rng.random(n) < 0.3
    expected = G.copy()
    expected[drop] = MISSING
    np.testing.assert_array_equal(np.asarray(P.with_missing(drop)), expected)
    with pytest.raises(IndexError):
        P[:, m]


def test_genotypes_methods_match_int8_storage():
    rng = np.random.default_rng(3)
    G = _calls(rng, 37, 29)
    G[:, 3] = MISSING
    dense, packed = Genotypes(G), Genotypes(G, packed=True)
    assert packed.packed and not dense.packed
    assert packed.G.nbytes == 29 * 10
    np.testing.assert_array_equal(packed.allele_freqs(), dense.allele_freqs())
    np.testing.assert_array_equal(packed.missing_rates(), dense.missing_rates())
    for got, want in ((packed.filter_variants(min_mac=2, max_missing=0.2),
                       dense.filter_variants(min_mac=2, max_missing=0.2)),
                      (packed.filter_samples([5, 0, 9]), dense.filter_samples([5, 0, 9])),
                      (packed.variant_mask(np.arange(29) % 2 == 0), dense.variant_mask(np.arange(29) % 2 == 0))):
        assert got.packed and not want.packed
        np.testing.assert_array_equal(np.asarray(got.G), want.G)
        np.testing.assert_array_equal(got.variant_ids, want.variant_ids)
    wanted = ["S4", "absent", "S0"]
    got, found = packed.align_samples(wanted, strict=False)
    want, want_found = dense.align_samples(wanted, strict=False)
    assert got.packed
    np.testing.assert_array_equal(found, want_found)
    np.testing.assert_array_equal(np.asarray(got.G), want.G)
    for a, b in zip(packed.iter_snp_blocks(block=7), dense.iter_snp_blocks(block=7)):
        np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize("n", [9, 12])
def test_plink_packed_and_mapped_reads_match_and_write_back(tmp_path, n):
    rng = np.random.default_rng(n)
    G = _calls(rng, n, 23)
    source = Genotypes(G, chromosome=np.repeat([1, 2], [11, 12]), position=np.arange(23) * 10)
    write_plink(source, str(tmp_path / "dense"))
    dense = read_plink(str(tmp_path / "dense"))
    for options in ({"packed": True}, {"mmap": True}):
        packed = read_plink(str(tmp_path / "dense"), **options)
        assert packed.packed
        assert isinstance(packed.G.data, np.memmap) == bool(options.get("mmap"))
        np.testing.assert_array_equal(np.asarray(packed.G), dense.G)
        np.testing.assert_array_equal(packed.position, dense.position)
        np.testing.assert_array_equal(np.asarray(read_plink(str(tmp_path / "dense"), max_variants=5,
                                                           **options).G), dense.G[:, :5])
        write_plink(packed, str(tmp_path / "again"))
        for suffix in ("bed", "bim", "fam"):
            assert (tmp_path / f"again.{suffix}").read_bytes() == (tmp_path / f"dense.{suffix}").read_bytes()
    assert read_plink(str(tmp_path / "dense"), max_variants=0, mmap=True).n_variants == 0
    with open(tmp_path / "dense.bed", "r+b") as fh:
        fh.truncate(3 + 2 * ((n + 3) // 4))
    with pytest.raises(ValueError, match="truncated"):
        read_plink(str(tmp_path / "dense"), mmap=True)


def test_preparation_ignores_unused_padding_bits():
    rng = np.random.default_rng(5)
    G = _calls(rng, 13, 11)
    gt = Genotypes(G, packed=True)
    groups = np.arange(11) % 2
    clean = LocoGenotypes(gt, groups, dtype=np.float64)
    noisy = gt.G.data.copy()
    noisy[:, -1] |= 0b11111100  # 13 samples: one used slot in the last byte
    noisy_gt = Genotypes(PackedCalls(noisy, 13))
    np.testing.assert_array_equal(np.asarray(noisy_gt.G), G)
    lg = LocoGenotypes(noisy_gt, groups, dtype=np.float64)
    for name in ("mean", "sd", "_projection", "zz"):
        np.testing.assert_array_equal(getattr(lg, name), getattr(clean, name))
    for (_, _, a), (_, _, b) in zip(lg.blocks(), clean.blocks()):
        np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize("method, options", [
    ("kvik", {"heritability_method": "he"}),
    ("bolt-inf", {}),
    ("bolt", {"min_cv_gain": -1.0}),
])
def test_two_step_associations_are_identical_on_packed_calls(method, options):
    rng = np.random.default_rng(11)
    n, m = 301, 480  # n is not a multiple of four
    G = _calls(rng, n, m, missing=0.02)
    chromosome = np.repeat(np.arange(1, 5), m // 4)
    position = np.tile(np.arange(m // 4) * 1000, 4)
    X = rng.normal(size=(n, 2))
    y = (G[:, :20].clip(0) @ rng.normal(size=20)) * 0.2 + X @ [0.3, -0.2] + rng.normal(size=n)
    results = [gwas(y, Genotypes(G, chromosome=chromosome, position=position, packed=packed),
                    X=X, method=method, random_state=3, **options) for packed in (False, True)]
    for name in ("p", "beta", "se"):
        np.testing.assert_array_equal(getattr(results[1], name), getattr(results[0], name))
