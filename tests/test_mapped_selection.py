"""Selections must not turn mapped cohorts into resident genotype copies."""

import numpy as np
import pytest

from mixmogam import Genotypes
from mixmogam._packed import PackedCalls
from mixmogam.io.plink import read_plink, write_plink


@pytest.fixture
def panel(tmp_path):
    calls = (np.arange(7 * 19).reshape(7, 19) % 4 - 1).astype(np.int8)
    prefix = tmp_path / "source"
    write_plink(Genotypes(calls), prefix)
    return calls, read_plink(prefix, mmap=True)


@pytest.mark.parametrize("indices", [np.arange(19), np.arange(2, 8), [16, 2, 16, 0],
                                     np.arange(19) % 3 == 0, [], [-1, 0]])
def test_mapped_variants_preserve_values_order_and_storage(panel, indices, tmp_path, monkeypatch):
    import mixmogam._packed as packed
    monkeypatch.setattr(packed, "_SPILL_BYTES", 4)  # at most two packed rows per write
    calls, gt = panel
    if len(indices) == 0:
        indices = np.array([], dtype=int)
    sub = gt.variant_mask(indices, directory=tmp_path)
    np.testing.assert_array_equal(np.asarray(sub.G), calls[:, indices])
    np.testing.assert_array_equal(sub.variant_ids, gt.variant_ids[indices])
    if sub.n_variants:
        assert isinstance(sub.G.data, np.memmap)
        assert not sub.G.data.flags.writeable
    assert not list(tmp_path.glob("mixmogam-variants-*"))
    np.testing.assert_array_equal(np.asarray(gt.G), calls)


def test_identity_and_contiguous_readonly_selections_do_not_copy(panel, monkeypatch):
    _, gt = panel
    def forbidden(*args, **kwargs):
        pytest.fail("identity or contiguous selection copied mapped calls")
    monkeypatch.setattr(PackedCalls, "_copy_variants_to_file", forbidden)
    monkeypatch.setattr(PackedCalls, "_take_samples_to_file", forbidden)
    for sub in (gt.filter_variants(), gt.filter_samples(np.arange(gt.n_samples)),
                gt.align_samples(gt.sample_ids)[0], gt.variant_mask(np.arange(2, 8))):
        assert np.shares_memory(sub.G.data, gt.G.data)
        assert isinstance(sub.G.data, np.memmap)


@pytest.mark.parametrize("compiled", [True, False])
def test_mapped_sample_selection_and_missing_alignment(panel, tmp_path, monkeypatch, compiled):
    import mixmogam._packed as packed
    if compiled and not packed.HAS_NUMBA:
        pytest.skip("requires Numba")
    monkeypatch.setattr(packed, "HAS_NUMBA", compiled)
    monkeypatch.setattr(packed, "_SPILL_BYTES", 64)
    calls, gt = panel
    sub = gt.filter_samples([5, 1, 3], directory=tmp_path)
    assert isinstance(sub.G.data, np.memmap)
    np.testing.assert_array_equal(np.asarray(sub.G), calls[[5, 1, 3]])
    aligned, found = gt.align_samples([gt.sample_ids[5], "absent", gt.sample_ids[1]], strict=False)
    assert isinstance(aligned.G.data, np.memmap)
    np.testing.assert_array_equal(found, [True, False, True])
    np.testing.assert_array_equal(np.asarray(aligned.G), np.vstack([calls[5], -np.ones(19), calls[1]]))
    np.testing.assert_array_equal(np.asarray(gt.G), calls)


def test_writable_selection_does_not_alias():
    gt = Genotypes(np.array([[0, 1], [2, 0]], dtype=np.int8), packed=True)
    for sub in (gt.variant_mask(np.arange(2)), gt.filter_samples([0, 1]),
                gt.align_samples(gt.sample_ids)[0]):
        assert not np.shares_memory(gt.G.data, sub.G.data)
        sub.G.data[:] = 0
        np.testing.assert_array_equal(np.asarray(gt.G), [[0, 1], [2, 0]])
