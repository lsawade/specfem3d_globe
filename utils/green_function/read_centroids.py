#!/usr/bin/env python3
"""
Example: read centroids.bin and print the Morton codes + centroid coordinates.

Usage:
    python read_centroids.py /path/to/gf_database/centroids.bin
"""

import argparse
import struct
import sys
from pathlib import Path


def read_centroids(filepath):
    """Read centroids.bin and return list of (morton_hex, cx, cy, cz)."""
    record_size = 32  # uint64 + 3x float64
    data = filepath.read_bytes()

    if len(data) % record_size != 0:
        print(f"Warning: file size {len(data)} not a multiple of {record_size}")

    n = len(data) // record_size
    records = []
    for i in range(n):
        morton, cx, cy, cz = struct.unpack_from("<Qddd", data, i * record_size)
        morton_hex = f"{morton:016X}"
        records.append((morton_hex, cx, cy, cz))

    return records


def main():
    parser = argparse.ArgumentParser(description="Read centroids.bin")
    parser.add_argument("filepath", type=Path, help="Path to centroids.bin")
    args = parser.parse_args()

    if not args.filepath.is_file():
        print(f"Error: {args.filepath} not found")
        sys.exit(1)

    records = read_centroids(args.filepath)
    print(f"{'Morton hex':<18s} {'cx':>20s} {'cy':>20s} {'cz':>20s}")
    print("-" * 80)
    for morton_hex, cx, cy, cz in records:
        print(f"{morton_hex:<18s} {cx:>20.12e} {cy:>20.12e} {cz:>20.12e}")

    print(f"\n{len(records)} records")


if __name__ == "__main__":
    main()
