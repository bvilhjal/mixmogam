"""QC agrees with independent call counts while bounding decoded input."""

import numpy as np
import pytest

import mixmogam.genotypes as genotypes
from mixmogam._packed import PackedCalls


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("budget", [1, 17, 10_000])
def test_qc_counts_and_filtered_metadata(monkeypatch, packed, budget):
    calls = np.array([[0, 2, -1, 1, 0], [1, 2, -1, -1, 0],
                      [2, 2, -1, 2, -1], [-1, 2, -1, 1, 0]], dtype=np.int8)
    gt = genotypes.Genotypes(calls, packed=packed,
                             allele1=["A"] * 5, allele2=["T"] * 5,
                             variant_ids=["a", "b", "c", "d", "e"])
    monkeypatch.setattr(genotypes, "_QC_BLOCK_CELLS", budget)
    decode = PackedCalls._decode_columns
    def bounded_decode(self, indices):
        assert indices.size * calls.shape[0] <= max(budget, calls.shape[0])
        return decode(self, indices)
    monkeypatch.setattr(PackedCalls, "_decode_columns", bounded_decode)
    np.testing.assert_array_equal(gt.allele_freqs(), [0.5, 1, np.nan, 2 / 3, 0])
    np.testing.assert_allclose(gt.allele_freqs(minor=True), [0.5, 0, np.nan, 1 / 3, 0])
    np.testing.assert_array_equal(gt.missing_rates(), [0.25, 0, 1, 0.25, 0.25])
    filtered = gt.filter_variants(min_mac=2, max_missing=0.25)
    assert filtered.packed == packed
    assert list(filtered.variant_ids) == ["a", "d"]
    assert list(filtered.allele1) == ["A", "A"]
    # Inspect one column at a time to retain the test's tiny decode budget.
    for j, original in enumerate([0, 3]):
        np.testing.assert_array_equal(filtered.G[:, j], calls[:, original])


@pytest.mark.parametrize("packed", [False, True])
def test_empty_variant_qc(packed):
    gt = genotypes.Genotypes(np.empty((5, 0), dtype=np.int8), packed=packed)
    assert gt.allele_freqs().size == gt.missing_rates().size == 0
    assert gt.filter_variants().G.shape == (5, 0)


@pytest.mark.parametrize("key", [slice(2, 17, 3), slice(None, None, -1),
                               slice(-3, None), slice(15, 2, -2), slice(3, 3)])
def test_packed_slice_indices_follow_numpy(key):
    from mixmogam._packed import _column_indices

    np.testing.assert_array_equal(_column_indices(key, 19), np.arange(19)[key])


def test_packed_small_slice_does_not_allocate_genome_wide_indices():
    import tracemalloc
    from mixmogam._packed import _column_indices

    tracemalloc.start()
    try:
        selected = _column_indices(slice(100, 104), 1_000_000)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    np.testing.assert_array_equal(selected, [100, 101, 102, 103])
    assert peak < 100_000  # A whole-genome index costs 8 MB for this slice.
