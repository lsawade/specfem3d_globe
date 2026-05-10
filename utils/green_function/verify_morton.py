#!/usr/bin/env python3
"""
Verify Morton encoding in centroids.bin against a Python reimplementation.

Usage:
    python verify_morton.py [path/to/centroids.bin]

Default path: DB_TESTING/db_gen/OUTPUT_FILES/gf_database/elements/centroids.bin
"""

import struct
import sys
import numpy as np


def spread_bits(v):
    """Spread 21 low bits of v into every-third-bit positions."""
    v = v & 0x1FFFFF  # 21 bits
    v = (v | (v << 32)) & 0x001F00000000FFFF
    v = (v | (v << 16)) & 0x001F0000FF0000FF
    v = (v | (v <<  8)) & 0x100F00F00F00F00F
    v = (v | (v <<  4)) & 0x10C30C30C30C30C3
    v = (v | (v <<  2)) & 0x1249249249249249
    return v


def morton_encode_3d(cx, cy, cz):
    """Quantize float32 coordinates from [-1,1] to 21-bit integers and
    bit-interleave into a 63-bit Morton code.

    Must match the Fortran implementation in green_function_morton.F90.
    """
    cx32 = np.float32(cx)
    cy32 = np.float32(cy)
    cz32 = np.float32(cz)
    MORTON_BINS = 2**21
    ix = int(np.float32((cx32 + np.float32(1.0)) * np.float32(0.5) * np.float32(MORTON_BINS - 1)))
    iy = int(np.float32((cy32 + np.float32(1.0)) * np.float32(0.5) * np.float32(MORTON_BINS - 1)))
    iz = int(np.float32((cz32 + np.float32(1.0)) * np.float32(0.5) * np.float32(MORTON_BINS - 1)))
    ix = max(0, min(ix, MORTON_BINS - 1))
    iy = max(0, min(iy, MORTON_BINS - 1))
    iz = max(0, min(iz, MORTON_BINS - 1))
    return spread_bits(ix) | (spread_bits(iy) << 1) | (spread_bits(iz) << 2)


def test_corner_cases():
    """Test spread_bits and morton_encode_3d with known values."""
    print("Corner case tests:")
    assert spread_bits(0) == 0, f"spread_bits(0) = {spread_bits(0)}"
    print(f"  spread_bits(0) = {spread_bits(0)} [OK]")

    assert spread_bits(1) == 1, f"spread_bits(1) = {spread_bits(1)}"
    print(f"  spread_bits(1) = {spread_bits(1)} [OK]")

    assert spread_bits(7) == 73, f"spread_bits(7) = {spread_bits(7)}"
    print(f"  spread_bits(7) = {spread_bits(7)} [OK]")

    m = morton_encode_3d(-1, -1, -1)
    assert m == 0, f"morton(-1,-1,-1) = {m:016X}"
    print(f"  morton(-1,-1,-1) = {m:016X} [OK]")

    m = morton_encode_3d(1, 1, 1)
    assert m == (1 << 63) - 1, f"morton(1,1,1) = {m:016X}"
    print(f"  morton( 1, 1, 1) = {m:016X} [OK]")

    print()


def verify_centroids(filepath):
    """Read centroids.bin and verify each Morton code against Python encoding."""
    with open(filepath, "rb") as f:
        data = f.read()

    n_records = len(data) // 32
    print(f"Records in centroids.bin: {n_records}")

    mismatches = 0
    for i in range(n_records):
        offset = i * 32
        morton_stored, cx, cy, cz = struct.unpack("<q3d", data[offset : offset + 32])
        morton_computed = morton_encode_3d(cx, cy, cz)
        if morton_stored != morton_computed:
            mismatches += 1
            print(
                f"  MISMATCH elem {i}: stored={morton_stored:016X} "
                f"computed={morton_computed:016X}"
            )
            print(f"    coords: {cx:.15e} {cy:.15e} {cz:.15e}")

    if mismatches == 0:
        print(f"All {n_records} Morton codes match Python implementation!")
    else:
        print(f"{mismatches}/{n_records} mismatches found")
        sys.exit(1)


if __name__ == "__main__":
    test_corner_cases()

    if len(sys.argv) > 1:
        filepath = sys.argv[1]
    else:
        filepath = "DB_TESTING/db_gen/OUTPUT_FILES/gf_database/elements/centroids.bin"

    verify_centroids(filepath)
