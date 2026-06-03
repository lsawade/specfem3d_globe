#!/usr/bin/env python3
"""
Cross-validate a Green function database against forward simulations.

Runs both force and CMT validation from a single base directory,
auto-deriving all parameters from the standard directory layout.

Expected directory structure:
    basedir/
    ├── db_base/DATA/STATIONS              (station list for GF database)
    ├── db_base/OUTPUT_FILES/values_from_mesher.h  (min resolved period)
    ├── GFDB/                              (GF database)
    ├── forward/OUTPUT_FILES/              (force validation forward sim)
    ├── forward_cmt/OUTPUT_FILES/          (CMT validation forward sim)
    └── validation_data/CMTSOLUTION        (CMT source for validation)

Usage:
    python xvalidate.py EXAMPLES/green_function_database/regional
"""

import argparse
import re
import sys
from pathlib import Path

# Import cross_validate from sibling module
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gf_cross_validate import cross_validate


def parse_station_from_file(stations_path):
    """Read the first station from a SPECFEM STATIONS file.

    Returns NET.STA string (e.g. 'IU.SJG').
    """
    with open(stations_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                return f"{parts[1]}.{parts[0]}"
    raise ValueError(f"No station found in {stations_path}")


def parse_min_period(values_from_mesher_path):
    """Extract minimum resolved period from values_from_mesher.h.

    Parses the line:
      ! the (approximate) minimum period resolved will be =  69.4968262  (s)
    """
    text = values_from_mesher_path.read_text()
    m = re.search(r"minimum period resolved will be\s*=\s*([\d.Ee+-]+)", text)
    if m:
        return float(m.group(1))
    raise ValueError(
        f"Could not parse minimum period from {values_from_mesher_path}"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Cross-validate GF database (force + CMT) from base directory"
    )
    parser.add_argument(
        "basedir", type=Path,
        help="Base directory (e.g., EXAMPLES/green_function_database/regional)",
    )
    parser.add_argument(
        "--force-only", action="store_true",
        help="Run only the force source validation",
    )
    parser.add_argument(
        "--cmt-only", action="store_true",
        help="Run only the CMT validation",
    )
    args = parser.parse_args()

    basedir = args.basedir.resolve()

    # Derive all paths from directory structure
    gf_db_path = basedir / "GFDB"
    forward_force_output = basedir / "forward" / "OUTPUT_FILES"
    forward_cmt_output = basedir / "forward_cmt" / "OUTPUT_FILES"
    cmtsolution = basedir / "validation_data" / "CMTSOLUTION"
    stations_file = basedir / "db_base" / "DATA" / "STATIONS"
    values_from_mesher = basedir / "db_base" / "OUTPUT_FILES" / "values_from_mesher.h"
    output_dir = basedir / "validation_output"

    # Validate required paths
    missing = []
    for label, path in [
        ("GF database", gf_db_path),
        ("STATIONS file", stations_file),
        ("values_from_mesher.h", values_from_mesher),
    ]:
        if not path.exists():
            missing.append(f"  {label}: {path}")

    if not args.force_only:
        if not forward_cmt_output.exists():
            missing.append(f"  CMT forward output: {forward_cmt_output}")
        if not cmtsolution.exists():
            missing.append(f"  CMTSOLUTION: {cmtsolution}")

    if not args.cmt_only:
        if not forward_force_output.exists():
            missing.append(f"  Force forward output: {forward_force_output}")

    if missing:
        print("ERROR: Missing required paths:")
        print("\n".join(missing))
        sys.exit(1)

    # Parse derived parameters
    station = parse_station_from_file(stations_file)
    highpass_period = parse_min_period(values_from_mesher)

    print(f"Base directory: {basedir}")
    print(f"Station: {station}")
    print(f"Highpass period: {highpass_period:.1f} s (from mesh)")

    output_dir.mkdir(exist_ok=True)

    # Run force validation
    if not args.cmt_only:
        print("\n" + "=" * 60)
        print("FORCE SOURCE VALIDATION")
        print("=" * 60)
        cross_validate(
            gf_db_path=gf_db_path,
            forward_output=forward_force_output,
            station=station,
            force=True,
            highpass_period=highpass_period,
            output=output_dir / "xvalidate_force.svg",
        )

    # Run CMT validation
    if not args.force_only:
        print("\n" + "=" * 60)
        print("CMT SOURCE VALIDATION")
        print("=" * 60)
        cross_validate(
            gf_db_path=gf_db_path,
            forward_output=forward_cmt_output,
            station=station,
            force=False,
            cmtsolution=cmtsolution,
            highpass_period=highpass_period,
            output=output_dir / "xvalidate_cmt.svg",
        )

    print(f"\nValidation complete. Plots saved to {output_dir}/")


if __name__ == "__main__":
    main()
