"""Small allocation helpers for readers of diploid hard calls."""

from __future__ import annotations

import numpy as np

from mixmogam._packed import PackedCalls


def variant_limit(value):
    if value is not None and (isinstance(value, (bool, np.bool_))
                              or not isinstance(value, (int, np.integer)) or value < 0):
        raise ValueError("max_variants must be a non-negative integer")
    return None if value is None else int(value)


def allocate_calls(n, m, packed):
    if packed:
        return PackedCalls(np.empty((m, (n + 3) // 4), dtype=np.uint8), n)
    return np.empty((n, m), dtype=np.int8, order="F")


def store_calls(out, start, calls):
    """Store a samples-by-variants block without an unpacked panel copy."""
    stop = start + calls.shape[1]
    if isinstance(out, PackedCalls):
        out.data[start:stop] = PackedCalls.pack(calls).data
    else:
        out[:, start:stop] = calls


def allocate_metadata(m):
    return {name: np.empty(m, dtype=np.int64 if name == "position" else object)
            for name in ("chromosome", "position", "variant_ids", "allele1", "allele2")}


def set_metadata(meta, j, chromosome, position, variant_id, ref, alt):
    if not ref or ref == "." or not alt or alt in (".", "<NON_REF>", "<*>") or "," in alt:
        raise ValueError(f"variant {variant_id}: only biallelic records with a defined ALT are supported")
    if ref == alt:
        raise ValueError(f"variant {variant_id}: REF and ALT must differ")
    try:
        position = int(position)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"variant {variant_id}: invalid position") from exc
    if position < 0:
        raise ValueError(f"variant {variant_id}: position must be non-negative")
    meta["chromosome"][j] = chromosome
    meta["position"][j] = position
    meta["variant_ids"][j] = variant_id
    meta["allele1"][j] = alt
    meta["allele2"][j] = ref


def validate_samples(samples):
    samples = np.asarray(samples, dtype=str)
    if samples.ndim != 1 or np.any(samples == "") or np.unique(samples).size != samples.size:
        raise ValueError("sample IDs must be nonempty and unique")
    return samples
