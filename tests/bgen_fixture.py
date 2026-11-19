"""Write the diploid layout-2 fields independently of either BGEN reader."""

import struct
import zlib

import numpy as np


def write_bgen(path, samples, records, compression="zlib", embedded=True):
    def string(value, size):
        data = str(value).encode()
        return len(data).to_bytes(size, "little") + data

    n = len(samples)
    ids = b"".join(string(s, 2) for s in samples)
    sample_block = struct.pack("<II", 8 + len(ids), n) + ids if embedded else b""
    flag = {None: 0, "zlib": 1, "zstd": 2}[compression]
    if flag == 2:
        try:
            from compression.zstd import compress
        except ImportError:
            import pytest
            compress = pytest.importorskip("zstandard").ZstdCompressor().compress
    else:
        compress = zlib.compress

    with path.open("wb") as fh:
        fh.write(struct.pack("<IIII4sI", 20 + len(sample_block), 20, len(records), n, b"bgen",
                             8 | (int(embedded) << 31) | flag))
        fh.write(sample_block)
        for variant, chrom, pos, ref, alt, bits, pairs, flags in records:
            fh.write(string(variant, 2) + string("rs1", 2) + string(chrom, 2))
            fh.write(struct.pack("<IH", pos, 2) + string(ref, 4) + string(alt, 4))
            if bits in (8, 16, 32):
                prob = np.asarray(pairs, dtype=f"<u{bits // 8}").tobytes()
            else:
                encoded = sum(int(v) << (i * bits) for i, v in enumerate(np.ravel(pairs)))
                prob = encoded.to_bytes((2 * n * bits + 7) // 8, "little")
            raw = struct.pack("<IHBB", n, 2, 2, 2) + bytes(flags) + bytes([0, bits]) + prob
            data = struct.pack("<I", len(raw)) + compress(raw) if flag else raw
            fh.write(struct.pack("<I", len(data)) + data)
