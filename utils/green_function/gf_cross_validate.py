#!/usr/bin/env python3
"""
Cross-validate Green function database against a forward simulation.

Reads the GF displacement from the HDF5 database, interpolates to the
source location using Lagrange interpolation on GLL nodes, applies
reciprocity rotation, and compares against the forward seismogram.

The forward simulation must use USE_FORCE_POINT_SOURCE = .true. with the
source placed at a location covered by the GF database (i.e., one of
the GF_LOCATIONS entries).

Usage:
    python gf_cross_validate.py <gf_database_path> <forward_output_path> [options]

Example:
    python gf_cross_validate.py \
        DB_TESTING/db_gen/OUTPUT_FILES/gf_database \
        DB_TESTING/forward/OUTPUT_FILES \
        --forward-solver-output DB_TESTING/forward/OUTPUT_FILES/output_solver.txt \
        --station IU.SJG
"""

import argparse
import os
import re
import sys
from pathlib import Path
from math import sin, cos, radians, sqrt, pi

import numpy as np

try:
    import h5py
except ImportError:
    print("Error: h5py is required. Install with: pip install h5py")
    sys.exit(1)

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Error: matplotlib is required. Install with: pip install matplotlib")
    sys.exit(1)

try:
    from obspy import read as obspy_read
except ImportError:
    print("Error: obspy is required. Install with: pip install obspy")
    sys.exit(1)


# GLL points for NGLL=5 (Gauss-Lobatto-Legendre quadrature nodes on [-1,1])
GLL5 = np.array([-1.0, -0.6546536707079771, 0.0, 0.6546536707079771, 1.0])


def lagrange_basis(xi, nodes):
    """Compute Lagrange basis function values at point xi for given nodes.

    Returns array of length len(nodes) with L_i(xi) for each i.
    """
    n = len(nodes)
    L = np.ones(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                L[i] *= (xi - nodes[j]) / (nodes[i] - nodes[j])
    return L


def interpolate_gll(data, xi, eta, gamma):
    """Interpolate data on a 5x5x5 GLL element to point (xi, eta, gamma).

    Parameters
    ----------
    data : ndarray, shape (nt, 5, 5, 5, 3, 3)
        Displacement on GLL nodes. Dims are (time, k, j, i, disp, force)
        where k,j,i are GLL node indices (reversed from Fortran i,j,k).
    xi, eta, gamma : float
        Reference coordinates in [-1, 1].

    Returns
    -------
    result : ndarray, shape (nt, 3, 3)
        Interpolated values.
    """
    Li = lagrange_basis(xi, GLL5)      # basis for i-direction (dim 3)
    Lj = lagrange_basis(eta, GLL5)     # basis for j-direction (dim 2)
    Lk = lagrange_basis(gamma, GLL5)   # basis for k-direction (dim 1)

    # Contract over GLL dimensions (dims 1, 2, 3)
    # data[t, k, j, i, d, f] -> sum over i,j,k with Li, Lj, Lk
    result = np.einsum("tkjidf,i,j,k->tdf", data, Li, Lj, Lk)
    return result


def parse_solver_output(filepath):
    """Parse source info from specfem3D output_solver.txt.

    Returns dict with: xi, eta, gamma, x, y, z, nu1, nu2, nu3, lat, lon, depth
    """
    text = Path(filepath).read_text()

    info = {}

    # Parse xi, eta, gamma
    m = re.search(
        r"at xi,eta,gamma coordinates\s*=\s*([-\dE.+]+)\s+([-\dE.+]+)\s+([-\dE.+]+)",
        text,
    )
    if m:
        info["xi"] = float(m.group(1))
        info["eta"] = float(m.group(2))
        info["gamma"] = float(m.group(3))

    # Parse (x, y, z)
    m = re.search(
        r"at \(x,y,z\)\s*=\s*([-\dE.+]+)\s+([-\dE.+]+)\s+([-\dE.+]+)", text
    )
    if m:
        info["x"] = float(m.group(1))
        info["y"] = float(m.group(2))
        info["z"] = float(m.group(3))

    # Parse nu1, nu2, nu3 (rotation matrix rows: North, East, Vertical)
    for name in ["nu1", "nu2", "nu3"]:
        m = re.search(
            rf"{name}\s*=\s*([-\dE.+]+)\s+([-\dE.+]+)\s+([-\dE.+]+)", text
        )
        if m:
            info[name] = np.array(
                [float(m.group(1)), float(m.group(2)), float(m.group(3))]
            )

    # Parse original source position
    m = re.search(r"latitude:\s*([-\dE.+]+)", text)
    if m:
        info["lat"] = float(m.group(1))
    m = re.search(r"longitude:\s*([-\dE.+]+)", text)
    if m:
        info["lon"] = float(m.group(1))

    # Parse force direction
    m = re.search(r"East direction:\s*([-\dE.+]+)", text)
    if m:
        info["force_E"] = float(m.group(1))
    m = re.search(r"North direction:\s*([-\dE.+]+)", text)
    if m:
        info["force_N"] = float(m.group(1))
    m = re.search(r"Vertical direction:\s*([-\dE.+]+)", text)
    if m:
        info["force_Z"] = float(m.group(1))

    return info


def lagrange_deriv(xi, nodes):
    """Compute derivative of Lagrange basis functions at point xi."""
    n = len(nodes)
    dL = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                prod = 1.0 / (nodes[i] - nodes[j])
                for k in range(n):
                    if k != i and k != j:
                        prod *= (xi - nodes[k]) / (nodes[i] - nodes[k])
                dL[i] += prod
    return dL


def locate_point_in_element(xyz_elem, target, maxiter=20, tol=1e-10):
    """Find reference coordinates (xi, eta, gamma) of a point in a spectral element.

    Uses Newton-Raphson iteration to invert the GLL coordinate mapping.

    Parameters
    ----------
    xyz_elem : ndarray, shape (5, 5, 5, 3)
        GLL node coordinates of the element.
    target : ndarray, shape (3,)
        Target point in Cartesian coordinates.
    maxiter : int
        Maximum Newton iterations.
    tol : float
        Convergence tolerance on residual norm.

    Returns
    -------
    xi, eta, gamma : float
        Reference coordinates. If all are in [-1, 1], the point is inside.
    converged : bool
        Whether Newton iteration converged.
    """
    xi, eta, gamma = 0.0, 0.0, 0.0

    for _ in range(maxiter):
        Li = lagrange_basis(xi, GLL5)
        Lj = lagrange_basis(eta, GLL5)
        Lk = lagrange_basis(gamma, GLL5)

        # Current mapped position
        current = np.array([
            np.einsum("kji,i,j,k->", xyz_elem[:, :, :, d], Li, Lj, Lk)
            for d in range(3)
        ])

        residual = target - current
        if np.linalg.norm(residual) < tol:
            return xi, eta, gamma, True

        # Jacobian via derivatives of Lagrange basis
        dLi = lagrange_deriv(xi, GLL5)
        dLj = lagrange_deriv(eta, GLL5)
        dLk = lagrange_deriv(gamma, GLL5)

        J = np.zeros((3, 3))
        for d in range(3):
            J[d, 0] = np.einsum("kji,i,j,k->", xyz_elem[:, :, :, d], dLi, Lj, Lk)
            J[d, 1] = np.einsum("kji,i,j,k->", xyz_elem[:, :, :, d], Li, dLj, Lk)
            J[d, 2] = np.einsum("kji,i,j,k->", xyz_elem[:, :, :, d], Li, Lj, dLk)

        dxi = np.linalg.solve(J, residual)
        xi += dxi[0]
        eta += dxi[1]
        gamma += dxi[2]

        # Clamp to avoid divergence
        xi = np.clip(xi, -1.5, 1.5)
        eta = np.clip(eta, -1.5, 1.5)
        gamma = np.clip(gamma, -1.5, 1.5)

    return xi, eta, gamma, False


def find_containing_element(gf_db_path, source_xyz):
    """Find the element that contains the source point.

    Uses Newton-Raphson point location on candidate elements (sorted by
    centroid distance). Returns the first element where the source maps to
    reference coordinates within [-1, 1].

    Returns (morton_hex, xi, eta, gamma, centroid_distance).
    """
    elements_dir = gf_db_path / "elements"
    candidates = []

    for elem_dir in sorted(elements_dir.iterdir()):
        if not elem_dir.is_dir():
            continue
        coord_file = elem_dir / "coordinates.h5"
        if not coord_file.exists():
            continue

        with h5py.File(coord_file, "r") as f:
            cx = float(np.asarray(f.attrs["cx"]).flat[0])
            cy = float(np.asarray(f.attrs["cy"]).flat[0])
            cz = float(np.asarray(f.attrs["cz"]).flat[0])

        dist = sqrt(
            (cx - source_xyz[0]) ** 2
            + (cy - source_xyz[1]) ** 2
            + (cz - source_xyz[2]) ** 2
        )
        candidates.append((dist, elem_dir.name))

    # Sort by centroid distance and check nearest candidates
    candidates.sort()

    for dist, morton in candidates[:20]:  # check top 20 nearest
        coord_file = elements_dir / morton / "coordinates.h5"
        with h5py.File(coord_file, "r") as f:
            xyz_elem = f["xyz"][:]

        xi, eta, gamma, converged = locate_point_in_element(xyz_elem, source_xyz)

        if converged and all(abs(v) <= 1.01 for v in (xi, eta, gamma)):
            return morton, xi, eta, gamma, dist

    # Fallback: return closest centroid (may not contain the point)
    morton = candidates[0][1]
    dist = candidates[0][0]
    print(f"  WARNING: no element contains the source, using closest centroid {morton}")
    return morton, 0.0, 0.0, 0.0, dist


def butterworth_highpass_sos(order, fc, fs):
    """Compute Butterworth highpass SOS coefficients.

    Returns ndarray of shape (order//2, 6) with [b0, b1, b2, a0=1, a1, a2] per section.
    """
    nsections = order // 2
    sos = np.zeros((nsections, 6))

    # Pre-warp cutoff frequency
    wc = 2.0 * fs * np.tan(pi * fc / fs)

    for k in range(1, nsections + 1):
        # Analog Butterworth pole angle
        angle = pi * (2 * k + order - 1) / (2 * order)
        pole_real = wc * cos(angle)

        # Bilinear transform: highpass has s^2 numerator in analog domain
        # H_hp(s) = s^2 / (s^2 - 2*Re(p)*s + |p|^2) for each conjugate pair
        # After bilinear transform z = (1 + s/(2fs)) / (1 - s/(2fs)):
        b0 = 4.0 * fs * fs
        b1 = -8.0 * fs * fs
        b2 = 4.0 * fs * fs

        a0 = 4.0 * fs * fs - 4.0 * fs * pole_real + wc * wc
        a1 = 2.0 * wc * wc - 8.0 * fs * fs
        a2 = 4.0 * fs * fs + 4.0 * fs * pole_real + wc * wc

        # Normalize
        sos[k - 1] = [b0 / a0, b1 / a0, b2 / a0, 1.0, a1 / a0, a2 / a0]

    return sos


def butterworth_sos(order, fc, fs):
    """Compute Butterworth lowpass SOS coefficients (matching Fortran implementation).

    Returns ndarray of shape (order//2, 6) with [b0, b1, b2, a0=1, a1, a2] per section.
    """
    nsections = order // 2
    sos = np.zeros((nsections, 6))

    # Pre-warp cutoff frequency
    wc = 2.0 * fs * np.tan(pi * fc / fs)
    wc2 = wc * wc

    for k in range(1, nsections + 1):
        # Analog Butterworth pole angle
        angle = pi * (2 * k + order - 1) / (2 * order)
        pole_real = wc * cos(angle)
        pole_imag = wc * sin(angle)

        # Bilinear transform coefficients
        b0 = wc2
        b1 = 2.0 * wc2
        b2 = wc2

        a0 = 4.0 * fs * fs - 4.0 * fs * pole_real + wc2
        a1 = 2.0 * wc2 - 8.0 * fs * fs
        a2 = 4.0 * fs * fs + 4.0 * fs * pole_real + wc2

        # Normalize
        sos[k - 1] = [b0 / a0, b1 / a0, b2 / a0, 1.0, a1 / a0, a2 / a0]

    return sos


def sos_filter_forward(x, b0, b1, b2, a1, a2):
    """Apply single SOS section forward (direct form II transposed)."""
    y = np.zeros_like(x)
    w1 = 0.0
    w2 = 0.0
    for i in range(len(x)):
        yi = b0 * x[i] + w1
        w1 = b1 * x[i] - a1 * yi + w2
        w2 = b2 * x[i] - a2 * yi
        y[i] = yi
    return y


def sosfiltfilt(x, sos):
    """Zero-phase SOS filter matching Fortran implementation.

    For each section: forward pass, then reverse + forward + reverse (backward pass).
    This matches the Fortran gf_sosfiltfilt which processes one section at a time.
    """
    y = x.copy()
    for isec in range(sos.shape[0]):
        b0, b1, b2, _, a1, a2 = sos[isec]
        # Forward pass
        y = sos_filter_forward(y, b0, b1, b2, a1, a2)
        # Backward pass (reverse, filter, reverse)
        y = y[::-1].copy()
        y = sos_filter_forward(y, b0, b1, b2, a1, a2)
        y = y[::-1].copy()
    return y


def filter_with_taper(x, sos, taper_fraction=0.05, pad_factor=1.0):
    """Taper signal edges to zero, zero-pad, filter, and trim.

    The taper brings the signal smoothly to zero at both ends so the
    transition to zero-padding is seamless. The padding then absorbs
    IIR filter edge transients. After filtering, the result is trimmed
    back to the original length.
    """
    n = len(x)
    n_taper = max(1, int(taper_fraction * n))
    n_pad = int(pad_factor * n)

    # Cosine taper to bring edges to zero
    taper = np.ones(n)
    taper[:n_taper] = 0.5 * (1.0 - np.cos(np.linspace(0, pi, n_taper)))
    taper[-n_taper:] = 0.5 * (1.0 - np.cos(np.linspace(pi, 0, n_taper)))
    tapered = x * taper

    # Zero-pad on both sides
    padded = np.concatenate([np.zeros(n_pad), tapered, np.zeros(n_pad)])

    # Filter and trim back to original length
    filtered = sosfiltfilt(padded, sos)
    return filtered[n_pad:n_pad + n]


def main():
    parser = argparse.ArgumentParser(
        description="Cross-validate GF database against forward simulation"
    )
    parser.add_argument(
        "gf_db_path", type=Path, help="Path to GF database directory"
    )
    parser.add_argument(
        "forward_output", type=Path, help="Path to forward OUTPUT_FILES directory"
    )
    parser.add_argument(
        "--forward-solver-output",
        type=Path,
        default=None,
        help="Path to forward output_solver.txt (default: <forward_output>/output_solver.txt)",
    )
    parser.add_argument(
        "--station",
        type=str,
        default="IU.SJG",
        help="Station name as NET.STA (default: IU.SJG)",
    )
    parser.add_argument(
        "--gf-solver-output",
        type=Path,
        default=None,
        help="Path to GF output_solver.txt (to get t0_gf; default: auto-detect)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output plot path (default: gf_cross_validation.png next to forward output)",
    )
    parser.add_argument(
        "--highpass-period",
        type=float,
        default=None,
        help="Apply zero-phase Butterworth highpass filter at this period (in seconds). "
             "Removes energy at periods longer than this value.",
    )
    parser.add_argument(
        "--lowpass-period",
        type=float,
        default=None,
        help="Override the lowpass filter period (in seconds). "
             "Default uses the GF database anti-alias cutoff frequency.",
    )
    args = parser.parse_args()

    gf_db = args.gf_db_path
    fwd_dir = args.forward_output
    station = args.station

    if args.forward_solver_output:
        solver_output = args.forward_solver_output
    else:
        solver_output = fwd_dir / "output_solver.txt"

    if args.output:
        output_plot = Path(args.output)
    else:
        output_plot = fwd_dir.parent / "gf_cross_validation.png"

    # ---------------------------------------------------------------
    # 1. Parse forward simulation output
    # ---------------------------------------------------------------
    print("Parsing forward solver output...")
    src = parse_solver_output(solver_output)
    print(f"  Source (x, y, z): ({src['x']:.8f}, {src['y']:.8f}, {src['z']:.8f})")
    print(f"  Source (xi, eta, gamma): ({src['xi']:.8f}, {src['eta']:.8f}, {src['gamma']:.8f})")
    print(f"  nu1 (N): {src['nu1']}")
    print(f"  nu2 (E): {src['nu2']}")
    print(f"  nu3 (Z): {src['nu3']}")

    # Determine which force component was used in forward simulation
    force_dir = np.array([src.get("force_E", 0), src.get("force_N", 0), src.get("force_Z", 0)])
    print(f"  Force direction (E, N, Z): {force_dir}")

    # Build the force direction in Cartesian from nu vectors
    # force_cart = force_E * nu_E + force_N * nu_N + force_Z * nu_Z
    # where nu_N = nu1, nu_E = nu2, nu_Z = nu3
    nu_force = (
        src.get("force_N", 0) * src["nu1"]
        + src.get("force_E", 0) * src["nu2"]
        + src.get("force_Z", 0) * src["nu3"]
    )
    print(f"  Force direction (Cartesian): {nu_force}")

    # ---------------------------------------------------------------
    # 2. Read mesh info
    # ---------------------------------------------------------------
    print("\nReading mesh info...")
    with h5py.File(gf_db / "mesh_info.h5", "r") as f:
        dt = float(np.asarray(f.attrs["dt"]).flat[0])
        nstep = int(np.asarray(f.attrs["nstep"]).flat[0])
        subsample_step = int(np.asarray(f.attrs["subsample_step"]).flat[0])
        nt_sub = int(np.asarray(f.attrs["nt_subsampled"]).flat[0])
    print(f"  dt={dt}, nstep={nstep}, subsample_step={subsample_step}, nt_sub={nt_sub}")

    # Read station metadata for STF filter params
    sta_file = gf_db / "stations" / f"{station}.h5"
    with h5py.File(sta_file, "r") as f:
        f_cutoff = float(np.asarray(f.attrs["f_cutoff"]).flat[0])
        hdur = float(np.asarray(f.attrs["hdur"]).flat[0])
    print(f"  f_cutoff={f_cutoff:.4f} Hz, hdur={hdur:.6f} s")

    # ---------------------------------------------------------------
    # 3. Find containing element
    # ---------------------------------------------------------------
    print("\nFinding containing element...")
    source_xyz = np.array([src["x"], src["y"], src["z"]])
    morton_hex, xi, eta, gamma, dist = find_containing_element(gf_db, source_xyz)
    print(f"  Containing element: {morton_hex} (centroid dist: {dist:.6e})")
    print(f"  Located at (xi={xi:.6f}, eta={eta:.6f}, gamma={gamma:.6f})")

    # ---------------------------------------------------------------
    # 4. Read GF displacement and interpolate
    # ---------------------------------------------------------------
    print("\nReading GF displacement...")
    gf_file = gf_db / "elements" / morton_hex / f"{station}.h5"
    with h5py.File(gf_file, "r") as f:
        # Python shape: (nt_sub, NGLLZ, NGLLY, NGLLX, 3_disp, 3_force)
        displ = f["displacement"][:]
    print(f"  Displacement shape: {displ.shape}, dtype: {displ.dtype}")

    # Interpolate to the exact source point using reference coordinates
    print(f"  Interpolating at (xi={xi:.6f}, eta={eta:.6f}, gamma={gamma:.6f})...")

    # displ shape: (nt, k, j, i, disp_comp, force_comp)
    # We want to contract over k, j, i dimensions (indices 1, 2, 3)
    # Result: G[nt, disp_comp, force_comp]
    G = interpolate_gll(displ, xi, eta, gamma)
    print(f"  Green tensor shape: {G.shape}")  # should be (nt_sub, 3, 3)

    # ---------------------------------------------------------------
    # 5. Apply reciprocity
    # ---------------------------------------------------------------
    # G[t, disp_xyz, force_NEZ] is the Cartesian displacement at the source
    # from a geographic force at the station.
    #
    # By reciprocity, the geographic displacement at the station from a
    # Cartesian force F_cart at the source is:
    #   u_geo_b(station) = sum_j F_cart_j * G[t, j, b]
    #
    # Our forward force is in geographic direction -> Cartesian via nu:
    #   F_cart = force_E * nu_E + force_N * nu_N + force_Z * nu_Z = nu_force
    #
    # So: u_b(station) = nu_force . G[:, :, b]

    print("\nApplying reciprocity...")
    gf_N = G[:, :, 0] @ nu_force  # geographic N at station
    gf_E = G[:, :, 1] @ nu_force  # geographic E at station
    gf_Z = G[:, :, 2] @ nu_force  # geographic Z at station
    print(f"  GF trace ranges: N=[{gf_N.min():.6e}, {gf_N.max():.6e}], "
          f"E=[{gf_E.min():.6e}, {gf_E.max():.6e}], "
          f"Z=[{gf_Z.min():.6e}, {gf_Z.max():.6e}]")

    # ---------------------------------------------------------------
    # 6. Convolve GF reconstruction with corrected STF
    # ---------------------------------------------------------------
    # The GF database has hdur_gf baked into the stored displacement.
    # The forward simulation uses hdur_fwd (from f0 / SOURCE_DECAY_MIMIC_TRIANGLE).
    # By Tarantola (2005, eq. 6.21), the convolution of two Gaussians gives a
    # Gaussian with sigma_total = sqrt(sigma_1^2 + sigma_2^2).
    # So we convolve the GF reconstruction with a corrected Gaussian:
    #   hdur_corrected = sqrt(hdur_fwd^2 - hdur_gf^2)
    # to produce the equivalent of the forward trace.

    SOURCE_DECAY_MIMIC_TRIANGLE = 1.628

    # Parse forward hdur from solver output
    fwd_text = Path(solver_output).read_text()
    m = re.search(r"Gaussian half duration:\s*([\d.Ee+-]+)", fwd_text)
    if m:
        hdur_fwd = float(m.group(1))
    else:
        # Fallback: read f0 from FORCESOLUTION
        m2 = re.search(r"half duration:\s*([\d.Ee+-]+)\s*seconds", fwd_text)
        hdur_fwd = float(m2.group(1)) / SOURCE_DECAY_MIMIC_TRIANGLE if m2 else hdur

    hdur_gf = hdur  # from station metadata

    print(f"\nSTF correction (Tarantola 2005, eq. 6.21):")
    print(f"  hdur_fwd (forward Gaussian):  {hdur_fwd:.6f} s")
    print(f"  hdur_gf  (GF database):       {hdur_gf:.6f} s")

    if hdur_fwd > hdur_gf:
        hdur_corrected = np.sqrt(hdur_fwd**2 - hdur_gf**2)
    else:
        print("  WARNING: hdur_fwd <= hdur_gf, no correction needed")
        hdur_corrected = 0.0

    print(f"  hdur_corrected:               {hdur_corrected:.6f} s")

    # Build corrected Gaussian STF and convolve with GF traces
    # Time axis for GF (subsampled dt)
    dt_sub = dt * subsample_step
    # STF kernel: symmetric Gaussian centered at t=0
    # Extent: +-4*hdur_corrected should be sufficient
    if hdur_corrected > 0:
        stf_half_len = int(4.0 * hdur_corrected / dt_sub)
        stf_t = np.arange(-stf_half_len, stf_half_len + 1) * dt_sub
        stf = np.exp(-(stf_t / hdur_corrected)**2) / (np.sqrt(pi) * hdur_corrected)
        stf *= dt_sub  # normalize for discrete convolution

        print(f"  STF kernel: {len(stf)} samples, peak={stf.max():.6e}")

        # Convolve GF traces with corrected STF using FFT for proper handling.
        # The STF is causal-shifted: we want the peak of the Gaussian at t=0
        # to correspond to no time shift in the output.
        from scipy.signal import fftconvolve
        gf_N = fftconvolve(gf_N, stf, mode="same")
        gf_E = fftconvolve(gf_E, stf, mode="same")
        gf_Z = fftconvolve(gf_Z, stf, mode="same")
        print(f"  After convolution - GF trace ranges:")
        print(f"    N=[{gf_N.min():.6e}, {gf_N.max():.6e}]")
        print(f"    E=[{gf_E.min():.6e}, {gf_E.max():.6e}]")
        print(f"    Z=[{gf_Z.min():.6e}, {gf_Z.max():.6e}]")

    # ---------------------------------------------------------------
    # 7. Read and filter forward seismograms
    # ---------------------------------------------------------------
    print("\nReading forward seismograms...")
    fwd_traces = {}
    for comp in ["BXN", "BXE", "BXZ"]:
        sac_file = fwd_dir / f"{station}.{comp}.sem.sac"
        tr = obspy_read(str(sac_file))[0]
        time = tr.times() + tr.stats.sac.b  # seconds from origin time
        fwd_traces[comp] = {"time": time, "amp": tr.data.astype(np.float64)}
        print(f"  {comp}: {len(tr.data)} samples, "
              f"amp range [{tr.data.min():.6e}, {tr.data.max():.6e}]")

    # Apply Butterworth filter to match GF anti-alias filtering, then subsample
    lp_fc = 1.0 / args.lowpass_period if args.lowpass_period else f_cutoff
    lp_filter = args.lowpass_period is not None
    print(f"\nApplying Butterworth lowpass filter (order=4, fc={lp_fc:.4f} Hz"
          f" [T={1/lp_fc:.1f}s], fs={1/dt:.4f} Hz)...")
    sos = butterworth_sos(4, lp_fc, 1.0 / dt)

    fwd_filtered = {}
    for comp in ["BXN", "BXE", "BXZ"]:
        amp = fwd_traces[comp]["amp"].astype(np.float64)
        filtered = filter_with_taper(amp, sos)
        fwd_filtered[comp] = filtered

    # Subsample forward traces
    fwd_sub = {}
    for comp in ["BXN", "BXE", "BXZ"]:
        fwd_sub[comp] = fwd_filtered[comp][::subsample_step]
    time_fwd = fwd_traces["BXN"]["time"][::subsample_step]

    # Construct GF time axis
    # GF stores at it = subsample_step, 2*subsample_step, ...
    # Time for snapshot isnap (1-based): t = (isnap * subsample_step - 1) * dt - t0_gf
    # Parse t0_gf from GF solver output
    t0_gf = None
    gf_solver_candidates = [
        args.gf_solver_output,
        fwd_dir.parent / "simulations" / station / "OUTPUT_FILES" / "output_solver.txt",
        fwd_dir.parent / "simulations" / station / "OUTPUT_FILES" / "output_solver_N.txt",
    ]
    for gf_solver_path in gf_solver_candidates:
        if gf_solver_path is not None and gf_solver_path.exists():
            gf_solver_text = gf_solver_path.read_text()
            m = re.search(r"start time\s*:\s*([-\dE.+]+)", gf_solver_text)
            if m:
                t0_gf = -float(m.group(1))  # start time is -t0
                print(f"  Found GF t0 from: {gf_solver_path}")
                break
    if t0_gf is None:
        # Fallback: t0_gf = T_min_period / 2
        t0_gf = hdur_gf * 5.0  # approximate
        print(f"  WARNING: could not find GF t0, estimating t0_gf={t0_gf:.4f}")

    dt_sub = dt * subsample_step
    time_gf = np.array([(isnap * subsample_step - 1) * dt - t0_gf
                         for isnap in range(1, nt_sub + 1)])

    print(f"  Forward subsampled: {len(time_fwd)} samples, t=[{time_fwd[0]:.2f}, {time_fwd[-1]:.2f}]")
    print(f"  GF: {nt_sub} samples, t=[{time_gf[0]:.2f}, {time_gf[-1]:.2f}]")
    print(f"  t0_gf={t0_gf:.4f}, t0_fwd={-fwd_traces['BXN']['time'][0]:.4f}")

    # Align traces using overlapping time window
    # Find common time range
    t_start = max(time_fwd[0], time_gf[0])
    t_end = min(time_fwd[-1], time_gf[-1])
    print(f"  Common time window: [{t_start:.2f}, {t_end:.2f}] s")

    # Use GF time axis as reference (it has fewer samples)
    # Interpolate forward trace onto GF time points within the common window
    mask_gf = (time_gf >= t_start) & (time_gf <= t_end)
    time_common = time_gf[mask_gf]

    gf_traces = {
        "BXN": gf_N[mask_gf],
        "BXE": gf_E[mask_gf],
        "BXZ": gf_Z[mask_gf],
    }

    fwd_interp = {}
    for comp in ["BXN", "BXE", "BXZ"]:
        fwd_interp[comp] = np.interp(time_common, time_fwd, fwd_sub[comp])

    n_compare = len(time_common)
    print(f"  Comparison samples: {n_compare}")

    # ---------------------------------------------------------------
    # 7b. Optional lowpass filter on GF traces
    # ---------------------------------------------------------------
    # The forward traces are always lowpass-filtered before subsampling (section 7).
    # The GF traces already have the database anti-alias filter baked in.
    # When --lowpass-period is explicitly set, apply the same lowpass to GF traces
    # so both are filtered identically.
    if lp_filter:
        lp_fc_sub = 1.0 / args.lowpass_period
        lp_fs_sub = 1.0 / dt_sub
        lp_order = 4
        print(f"\nApplying Butterworth lowpass filter to GF traces (order={lp_order}, "
              f"period={args.lowpass_period:.1f} s, fc={lp_fc_sub:.6f} Hz, fs={lp_fs_sub:.4f} Hz)...")
        lp_sos_sub = butterworth_sos(lp_order, lp_fc_sub, lp_fs_sub)
        for comp in ["BXN", "BXE", "BXZ"]:
            gf_traces[comp] = filter_with_taper(gf_traces[comp], lp_sos_sub)

    # ---------------------------------------------------------------
    # 7c. Optional highpass filter
    # ---------------------------------------------------------------
    if args.highpass_period is not None:
        hp_fc = 1.0 / args.highpass_period  # convert period to frequency
        hp_fs = 1.0 / dt_sub  # sampling rate of subsampled traces
        hp_order = 4
        print(f"\nApplying Butterworth highpass filter (order={hp_order}, "
              f"period={args.highpass_period:.1f} s, fc={hp_fc:.6f} Hz, fs={hp_fs:.4f} Hz)...")
        hp_sos = butterworth_highpass_sos(hp_order, hp_fc, hp_fs)
        for comp in ["BXN", "BXE", "BXZ"]:
            gf_traces[comp] = filter_with_taper(gf_traces[comp], hp_sos)
            fwd_interp[comp] = filter_with_taper(fwd_interp[comp], hp_sos)

    # ---------------------------------------------------------------
    # 8. Compute residuals
    # ---------------------------------------------------------------
    print("\nResiduals:")
    for comp in ["BXN", "BXE", "BXZ"]:
        fwd = fwd_interp[comp]
        gf = gf_traces[comp]
        residual = np.sqrt(np.sum((fwd - gf) ** 2))
        ref = np.sqrt(np.sum(fwd ** 2))
        rel_err = residual / ref if ref > 0 else float("inf")
        print(f"  {comp}: ||fwd - gf|| / ||fwd|| = {rel_err:.6e} "
              f"(abs: {residual:.6e}, ref: {ref:.6e})")

    # ---------------------------------------------------------------
    # 9. Plot comparison
    # ---------------------------------------------------------------
    print(f"\nPlotting to {output_plot}...")
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    comp_labels = {"BXN": "North", "BXE": "East", "BXZ": "Vertical"}

    for ax, comp in zip(axes, ["BXN", "BXE", "BXZ"]):
        fwd = fwd_interp[comp]
        gf = gf_traces[comp]

        residual = np.sqrt(np.sum((fwd - gf) ** 2))
        ref = np.sqrt(np.sum(fwd ** 2))
        rel_err = residual / ref if ref > 0 else float("inf")

        ax.plot(time_common, fwd, "b-", linewidth=0.8, label="Forward (filtered)")
        ax.plot(time_common, gf, "r--", linewidth=0.8, label="GF reconstructed")
        ax.set_ylabel(f"{comp} ({comp_labels[comp]})")
        ax.set_title(f"{comp} — relative error: {rel_err:.4e}")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time (s)")
    filter_info = ""
    if args.lowpass_period:
        filter_info += f", lowpass T<{args.lowpass_period:.0f}s"
    if args.highpass_period:
        filter_info += f", highpass T>{args.highpass_period:.0f}s"
    fig.suptitle(
        f"GF Cross-Validation: {station}{filter_info}\n"
        f"Source at element {morton_hex}\n"
        f"(xi={xi:.4f}, eta={eta:.4f}, gamma={gamma:.4f})",
        fontsize=11,
    )
    plt.tight_layout()
    plt.savefig(output_plot, dpi=150)
    print(f"Saved {output_plot}")
    plt.close()


if __name__ == "__main__":
    main()
