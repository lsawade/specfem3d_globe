#!/usr/bin/env python3
"""
Decode Morton codes back to approximate coordinates.

Usage:
    python decode_morton.py <morton_hex> [morton_hex ...]
    python decode_morton.py --file path/to/centroids.bin

Examples:
    python decode_morton.py 113B5BC4F83F7641
    python decode_morton.py --file DB_TESTING/db_gen/OUTPUT_FILES/gf_database/elements/centroids.bin
"""

import struct
import sys
import numpy as np


MORTON_BITS = 21
MORTON_BINS = 2**MORTON_BITS  # 2,097,152


def compact_bits(v):
    """Reverse of spread_bits: extract every-third bit into contiguous low bits."""
    v = v & 0x1249249249249249
    v = (v | (v >>  2)) & 0x10C30C30C30C30C3
    v = (v | (v >>  4)) & 0x100F00F00F00F00F
    v = (v | (v >>  8)) & 0x001F0000FF0000FF
    v = (v | (v >> 16)) & 0x001F00000000FFFF
    v = (v | (v >> 32)) & 0x00000000001FFFFF
    return v


def morton_decode_3d(morton):
    """Decode a 63-bit Morton code back to approximate (x, y, z) in [-1, 1]."""
    ix = compact_bits(morton)
    iy = compact_bits(morton >> 1)
    iz = compact_bits(morton >> 2)

    # map from [0, MORTON_BINS - 1] back to [-1, 1]
    x = ix / (MORTON_BINS - 1) * 2.0 - 1.0
    y = iy / (MORTON_BINS - 1) * 2.0 - 1.0
    z = iz / (MORTON_BINS - 1) * 2.0 - 1.0

    return x, y, z


def xyz_to_latlon(x, y, z):
    """Convert Cartesian unit-sphere coordinates to latitude/longitude (degrees).

    The coordinates are normalized to the unit sphere (r ~ 1 for surface elements,
    smaller for deeper elements).
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    lat = np.degrees(np.arcsin(z / r)) if r > 0 else 0.0
    lon = np.degrees(np.arctan2(y, x))
    return lat, lon, r


def decode_hex(hex_str):
    """Decode a Morton hex string and print coordinates."""
    morton = int(hex_str, 16)
    x, y, z = morton_decode_3d(morton)
    lat, lon, r = xyz_to_latlon(x, y, z)
    print(f"  Morton: {hex_str}")
    print(f"  Cartesian: x={x:.6f}  y={y:.6f}  z={z:.6f}")
    print(f"  Spherical: lat={lat:.4f}  lon={lon:.4f}  r={r:.6f}")
    print()
    return x, y, z


def decode_centroids_file(filepath):
    """Decode all entries in a centroids.bin file and compare with stored coords."""
    with open(filepath, "rb") as f:
        data = f.read()

    n_records = len(data) // 32
    print(f"Decoding {n_records} records from {filepath}\n")
    print(f"{'Morton hex':>20s}  {'dx':>10s} {'dy':>10s} {'dz':>10s}"
          f"  {'lat':>8s} {'lon':>9s} {'r':>8s}")
    print("-" * 85)

    max_err = 0.0
    for i in range(n_records):
        offset = i * 32
        morton_stored, cx, cy, cz = struct.unpack("<q3d", data[offset:offset + 32])
        hex_str = f"{morton_stored:016X}"

        xd, yd, zd = morton_decode_3d(morton_stored)
        lat, lon, r = xyz_to_latlon(xd, yd, zd)

        dx = xd - cx
        dy = yd - cy
        dz = zd - cz
        err = np.sqrt(dx**2 + dy**2 + dz**2)
        max_err = max(max_err, err)

        print(f"{hex_str}  {dx:+10.2e} {dy:+10.2e} {dz:+10.2e}"
              f"  {lat:8.3f} {lon:9.3f} {r:8.5f}")

    print("-" * 85)
    print(f"Max coordinate error from quantization: {max_err:.2e}")
    print(f"  (21-bit bin width: {2.0 / (MORTON_BINS - 1):.2e})")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    if sys.argv[1] == "--file":
        if len(sys.argv) < 3:
            print("Error: --file requires a path argument")
            sys.exit(1)
        decode_centroids_file(sys.argv[2])
    else:
        for hex_str in sys.argv[1:]:
            decode_hex(hex_str.strip())
