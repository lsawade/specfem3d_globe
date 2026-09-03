#!/usr/bin/env python3
"""
Cross-validate Green function database against a forward simulation.

Supports two source types:
  - Force source (--force): Reads GF displacement, interpolates to source,
    applies reciprocity with force direction. The forward simulation must
    use USE_FORCE_POINT_SOURCE = .true.
  - Moment tensor / CMT (default): Computes strain from GF displacement
    gradients, rotates the moment tensor to Cartesian coordinates, contracts
    strain with MT, and time-integrates to convert from Gaussian to Heaviside
    STF. The forward simulation must use USE_FORCE_POINT_SOURCE = .false.

In both cases, the source must be placed at a location covered by the GF
database (i.e., one of the GF_LOCATIONS entries).

Usage:
    python gf_cross_validate.py <gf_database_path> <forward_output_path> [options]

Examples:
    # CMT (default):
    python gf_cross_validate.py GFDB/ forward_cmt/OUTPUT_FILES --station IU.SJG

    # Force source:
    python gf_cross_validate.py GFDB/ forward/OUTPUT_FILES --force --station IU.SJG
"""

import argparse
import re
import sys
from pathlib import Path
from math import sin, cos, radians, sqrt, pi

import numpy as np
from scipy.signal import fftconvolve
from scipy.integrate import cumulative_trapezoid

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


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# GLL points for NGLL=5 (Gauss-Lobatto-Legendre quadrature nodes on [-1,1])
GLL5 = np.array([-1.0, -0.6546536707079771, 0.0, 0.6546536707079771, 1.0])

# SPECFEM3D non-dimensionalization constants (matching setup/constants.h.in)
R_PLANET = 6371000.0        # Earth radius in meters
RHOAV = 5514.3              # Average density in kg/m^3
GRAV = 6.67384e-11          # Gravitational constant in m^3 kg^-1 s^-2

# Moment tensor scaling: M_nondim = M_dyne_cm / SCALE_M
# From get_cmt.f90: scaleM = 1.d7 * RHOAV * (R_PLANET**5) * PI * GRAV * RHOAV
SCALE_M = 1.0e7 * RHOAV * (R_PLANET ** 5) * pi * GRAV * RHOAV

SOURCE_DECAY_MIMIC_TRIANGLE = 1.628


# ---------------------------------------------------------------------------
# GLL interpolation and derivatives
# ---------------------------------------------------------------------------

def lagrange_basis(xi, nodes):
    """Compute Lagrange basis function values at point xi for given nodes."""
    n = len(nodes)
    L = np.ones(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                L[i] *= (xi - nodes[j]) / (nodes[i] - nodes[j])
    return L


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
    Li = lagrange_basis(xi, GLL5)
    Lj = lagrange_basis(eta, GLL5)
    Lk = lagrange_basis(gamma, GLL5)
    return np.einsum("tkjidf,i,j,k->tdf", data, Li, Lj, Lk)


# ---------------------------------------------------------------------------
# Element location (Newton-Raphson)
# ---------------------------------------------------------------------------

def locate_point_in_element(xyz_elem, target, maxiter=20, tol=1e-10):
    """Find reference coordinates (xi, eta, gamma) of a point in a spectral element.

    Uses Newton-Raphson iteration to invert the GLL coordinate mapping.

    Parameters
    ----------
    xyz_elem : ndarray, shape (5, 5, 5, 3)
        GLL node coordinates of the element.
    target : ndarray, shape (3,)
        Target point in Cartesian coordinates.

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

        current = np.array([
            np.einsum("kji,i,j,k->", xyz_elem[:, :, :, d], Li, Lj, Lk)
            for d in range(3)
        ])

        residual = target - current
        if np.linalg.norm(residual) < tol:
            return xi, eta, gamma, True

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

    candidates.sort()

    for dist, morton in candidates[:20]:
        coord_file = elements_dir / morton / "coordinates.h5"
        with h5py.File(coord_file, "r") as f:
            xyz_elem = f["xyz"][:]

        xi, eta, gamma, converged = locate_point_in_element(xyz_elem, source_xyz)

        if converged and all(abs(v) <= 1.01 for v in (xi, eta, gamma)):
            return morton, xi, eta, gamma, dist

    morton = candidates[0][1]
    dist = candidates[0][0]
    print(f"  WARNING: no element contains the source, using closest centroid {morton}")
    return morton, 0.0, 0.0, 0.0, dist


# ---------------------------------------------------------------------------
# Butterworth filtering (matching Fortran implementation)
# ---------------------------------------------------------------------------

def butterworth_highpass_sos(order, fc, fs):
    """Compute Butterworth highpass SOS coefficients."""
    nsections = order // 2
    sos = np.zeros((nsections, 6))
    wc = 2.0 * fs * np.tan(pi * fc / fs)

    for k in range(1, nsections + 1):
        angle = pi * (2 * k + order - 1) / (2 * order)
        pole_real = wc * cos(angle)

        b0 = 4.0 * fs * fs
        b1 = -8.0 * fs * fs
        b2 = 4.0 * fs * fs

        a0 = 4.0 * fs * fs - 4.0 * fs * pole_real + wc * wc
        a1 = 2.0 * wc * wc - 8.0 * fs * fs
        a2 = 4.0 * fs * fs + 4.0 * fs * pole_real + wc * wc

        sos[k - 1] = [b0 / a0, b1 / a0, b2 / a0, 1.0, a1 / a0, a2 / a0]

    return sos


def butterworth_sos(order, fc, fs):
    """Compute Butterworth lowpass SOS coefficients (matching Fortran implementation)."""
    nsections = order // 2
    sos = np.zeros((nsections, 6))
    wc = 2.0 * fs * np.tan(pi * fc / fs)
    wc2 = wc * wc

    for k in range(1, nsections + 1):
        angle = pi * (2 * k + order - 1) / (2 * order)
        pole_real = wc * cos(angle)

        b0 = wc2
        b1 = 2.0 * wc2
        b2 = wc2

        a0 = 4.0 * fs * fs - 4.0 * fs * pole_real + wc2
        a1 = 2.0 * wc2 - 8.0 * fs * fs
        a2 = 4.0 * fs * fs + 4.0 * fs * pole_real + wc2

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
    """Zero-phase SOS filter matching Fortran implementation."""
    y = x.copy()
    for isec in range(sos.shape[0]):
        b0, b1, b2, _, a1, a2 = sos[isec]
        y = sos_filter_forward(y, b0, b1, b2, a1, a2)
        y = y[::-1].copy()
        y = sos_filter_forward(y, b0, b1, b2, a1, a2)
        y = y[::-1].copy()
    return y


def filter_with_taper(x, sos, taper_fraction=0.05, pad_factor=1.0):
    """Taper signal edges to zero, zero-pad, filter, and trim."""
    n = len(x)
    n_taper = max(1, int(taper_fraction * n))
    n_pad = int(pad_factor * n)

    taper = np.ones(n)
    taper[:n_taper] = 0.5 * (1.0 - np.cos(np.linspace(0, pi, n_taper)))
    taper[-n_taper:] = 0.5 * (1.0 - np.cos(np.linspace(pi, 0, n_taper)))
    tapered = x * taper

    padded = np.concatenate([np.zeros(n_pad), tapered, np.zeros(n_pad)])
    filtered = sosfiltfilt(padded, sos)
    return filtered[n_pad:n_pad + n]


# ---------------------------------------------------------------------------
# Source file parsers
# ---------------------------------------------------------------------------

def parse_solver_output(filepath):
    """Parse source info from specfem3D output_solver.txt.

    Returns dict with: xi, eta, gamma, x, y, z, nu1, nu2, nu3, lat, lon, depth
    """
    text = Path(filepath).read_text()
    info = {}

    m = re.search(
        r"at xi,eta,gamma coordinates\s*=\s*([-\dE.+]+)\s+([-\dE.+]+)\s+([-\dE.+]+)",
        text,
    )
    if m:
        info["xi"] = float(m.group(1))
        info["eta"] = float(m.group(2))
        info["gamma"] = float(m.group(3))

    m = re.search(
        r"at \(x,y,z\)\s*=\s*([-\dE.+]+)\s+([-\dE.+]+)\s+([-\dE.+]+)", text
    )
    if m:
        info["x"] = float(m.group(1))
        info["y"] = float(m.group(2))
        info["z"] = float(m.group(3))

    for name in ["nu1", "nu2", "nu3"]:
        m = re.search(
            rf"{name}\s*=\s*([-\dE.+]+)\s+([-\dE.+]+)\s+([-\dE.+]+)", text
        )
        if m:
            info[name] = np.array(
                [float(m.group(1)), float(m.group(2)), float(m.group(3))]
            )

    m = re.search(r"latitude:\s*([-\dE.+]+)", text)
    if m:
        info["lat"] = float(m.group(1))
    m = re.search(r"longitude:\s*([-\dE.+]+)", text)
    if m:
        info["lon"] = float(m.group(1))

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


def parse_cmtsolution(filepath):
    """Parse CMTSOLUTION file for moment tensor components and source parameters.

    Returns dict with: Mrr, Mtt, Mpp, Mrt, Mrp, Mtp (dyne-cm),
    half_duration, time_shift (seconds), latitude, longitude, depth (km).
    """
    info = {}
    with open(filepath) as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()
        if line.startswith("time shift:"):
            info["time_shift"] = float(line.split(":")[1])
        elif line.startswith("half duration:"):
            info["half_duration"] = float(line.split(":")[1])
        elif line.startswith("latitude:"):
            info["latitude"] = float(line.split(":")[1])
        elif line.startswith("longitude:"):
            info["longitude"] = float(line.split(":")[1])
        elif line.startswith("depth:"):
            info["depth"] = float(line.split(":")[1])
        elif line.startswith("Mrr:"):
            info["Mrr"] = float(line.split(":")[1])
        elif line.startswith("Mtt:"):
            info["Mtt"] = float(line.split(":")[1])
        elif line.startswith("Mpp:"):
            info["Mpp"] = float(line.split(":")[1])
        elif line.startswith("Mrt:"):
            info["Mrt"] = float(line.split(":")[1])
        elif line.startswith("Mrp:"):
            info["Mrp"] = float(line.split(":")[1])
        elif line.startswith("Mtp:"):
            info["Mtp"] = float(line.split(":")[1])

    return info


# ---------------------------------------------------------------------------
# Moment tensor rotation
# ---------------------------------------------------------------------------

def rotate_mt_to_xyz(Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, theta, phi):
    """Rotate moment tensor from spherical (r, theta, phi) to Cartesian (x, y, z).

    Matches SPECFEM3D locate_sources.f90 lines 258-273 exactly.
    """
    sint = sin(theta)
    cost = cos(theta)
    sinp = sin(phi)
    cosp = cos(phi)

    Mxx = (sint*sint*cosp*cosp*Mrr + cost*cost*cosp*cosp*Mtt + sinp*sinp*Mpp
           + 2.0*sint*cost*cosp*cosp*Mrt - 2.0*sint*sinp*cosp*Mrp
           - 2.0*cost*sinp*cosp*Mtp)

    Myy = (sint*sint*sinp*sinp*Mrr + cost*cost*sinp*sinp*Mtt + cosp*cosp*Mpp
           + 2.0*sint*cost*sinp*sinp*Mrt + 2.0*sint*sinp*cosp*Mrp
           + 2.0*cost*sinp*cosp*Mtp)

    Mzz = cost*cost*Mrr + sint*sint*Mtt - 2.0*sint*cost*Mrt

    Mxy = (sint*sint*sinp*cosp*Mrr + cost*cost*sinp*cosp*Mtt - sinp*cosp*Mpp
           + 2.0*sint*cost*sinp*cosp*Mrt
           + sint*(cosp*cosp - sinp*sinp)*Mrp
           + cost*(cosp*cosp - sinp*sinp)*Mtp)

    Mxz = (sint*cost*cosp*Mrr - sint*cost*cosp*Mtt
           + (cost*cost - sint*sint)*cosp*Mrt
           - cost*sinp*Mrp + sint*sinp*Mtp)

    Myz = (sint*cost*sinp*Mrr - sint*cost*sinp*Mtt
           + (cost*cost - sint*sint)*sinp*Mrt
           + cost*cosp*Mrp - sint*cosp*Mtp)

    return Mxx, Myy, Mzz, Mxy, Mxz, Myz


# ---------------------------------------------------------------------------
# Strain computation
# ---------------------------------------------------------------------------

def compute_inverse_jacobian(xyz_elem, xi, eta, gamma):
    """Compute the inverse Jacobian at a point within a spectral element."""
    Li = lagrange_basis(xi, GLL5)
    Lj = lagrange_basis(eta, GLL5)
    Lk = lagrange_basis(gamma, GLL5)
    dLi = lagrange_deriv(xi, GLL5)
    dLj = lagrange_deriv(eta, GLL5)
    dLk = lagrange_deriv(gamma, GLL5)

    J = np.zeros((3, 3))
    for d in range(3):
        J[d, 0] = np.einsum("kji,i,j,k->", xyz_elem[:, :, :, d], dLi, Lj, Lk)
        J[d, 1] = np.einsum("kji,i,j,k->", xyz_elem[:, :, :, d], Li, dLj, Lk)
        J[d, 2] = np.einsum("kji,i,j,k->", xyz_elem[:, :, :, d], Li, Lj, dLk)

    return np.linalg.inv(J)


def compute_strain(displ, xi, eta, gamma, Jinv):
    """Compute strain tensor from GF displacement at a point within an element.

    Parameters
    ----------
    displ : ndarray, shape (nt, 5, 5, 5, 3, 3)
        GF displacement on GLL nodes: (time, k, j, i, disp_comp, force_comp).
    xi, eta, gamma : float
        Reference coordinates of the target point.
    Jinv : ndarray, shape (3, 3)
        Inverse Jacobian: Jinv[ref_coord, phys_coord].

    Returns
    -------
    strain : ndarray, shape (nt, 6, 3)
        Strain in Voigt notation [xx, yy, zz, xy, xz, yz] for each
        force component (N, E, Z).
    """
    Li = lagrange_basis(xi, GLL5)
    Lj = lagrange_basis(eta, GLL5)
    Lk = lagrange_basis(gamma, GLL5)
    dLi = lagrange_deriv(xi, GLL5)
    dLj = lagrange_deriv(eta, GLL5)
    dLk = lagrange_deriv(gamma, GLL5)

    dG_dxi = np.einsum("tkjidb,i,j,k->tdb", displ, dLi, Lj, Lk)
    dG_deta = np.einsum("tkjidb,i,j,k->tdb", displ, Li, dLj, Lk)
    dG_dgamma = np.einsum("tkjidb,i,j,k->tdb", displ, Li, Lj, dLk)

    nt = displ.shape[0]
    dGdx = np.zeros((nt, 3, 3, 3))  # (time, disp_comp, phys_deriv, force_comp)
    for j in range(3):
        dGdx[:, :, j, :] = (dG_dxi * Jinv[0, j]
                             + dG_deta * Jinv[1, j]
                             + dG_dgamma * Jinv[2, j])

    # Symmetrize to get strain in Voigt notation: [xx, yy, zz, xy, xz, yz]
    strain = np.zeros((nt, 6, 3))
    strain[:, 0, :] = dGdx[:, 0, 0, :]                                  # eps_xx
    strain[:, 1, :] = dGdx[:, 1, 1, :]                                  # eps_yy
    strain[:, 2, :] = dGdx[:, 2, 2, :]                                  # eps_zz
    strain[:, 3, :] = 0.5 * (dGdx[:, 0, 1, :] + dGdx[:, 1, 0, :])     # eps_xy
    strain[:, 4, :] = 0.5 * (dGdx[:, 0, 2, :] + dGdx[:, 2, 0, :])     # eps_xz
    strain[:, 5, :] = 0.5 * (dGdx[:, 1, 2, :] + dGdx[:, 2, 1, :])     # eps_yz

    return strain


# ---------------------------------------------------------------------------
# High-level reconstruction functions
# ---------------------------------------------------------------------------

def read_mesh_info(gf_db):
    """Read mesh timing parameters from the GF database.

    Returns dict with keys: dt, nstep, subsample_step, nt_sub, and t0 (the
    time of the first time step relative to the source origin; None for older
    databases that did not store it).
    """
    with h5py.File(gf_db / "mesh_info.h5", "r") as f:
        info = {
            "dt": float(np.asarray(f.attrs["dt"]).flat[0]),
            "nstep": int(np.asarray(f.attrs["nstep"]).flat[0]),
            "subsample_step": int(np.asarray(f.attrs["subsample_step"]).flat[0]),
            "nt_sub": int(np.asarray(f.attrs["nt_subsampled"]).flat[0]),
            "t0": (float(np.asarray(f.attrs["t0"]).flat[0])
                   if "t0" in f.attrs else None),
        }
    return info


def read_station_meta(gf_db, station):
    """Read station metadata from the GF database.

    Returns dict with keys: f_cutoff, hdur, factor_force_nondim.
    """
    with h5py.File(gf_db / "stations" / f"{station}.h5", "r") as f:
        meta = {
            "f_cutoff": float(np.asarray(f.attrs["f_cutoff"]).flat[0]),
            "hdur": float(np.asarray(f.attrs["hdur"]).flat[0]),
            "factor_force_nondim": float(
                np.asarray(f.attrs["factor_force_source"]).flat[0]
            ),
        }
    return meta


def reconstruct_force(displ, xi, eta, gamma, src):
    """Reconstruct seismograms from GF database using force source reciprocity.

    Returns ndarray of shape (3, nt) with components [N, E, Z].
    """
    nu_force = (
        src.get("force_N", 0) * src["nu1"]
        + src.get("force_E", 0) * src["nu2"]
        + src.get("force_Z", 0) * src["nu3"]
    )

    G = interpolate_gll(displ, xi, eta, gamma)  # (nt, 3_disp, 3_force)

    gf = np.zeros((3, G.shape[0]))
    for i in range(3):
        gf[i] = G[:, :, i] @ nu_force
    return gf


def reconstruct_cmt(displ, xi, eta, gamma, gf_db, morton_hex,
                    src, cmt_path, factor_force_nondim, dt_sub):
    """Reconstruct seismograms from GF database using moment tensor source.

    Computes strain, rotates the MT, contracts, and time-integrates
    to convert from Gaussian to Heaviside STF.

    Returns ndarray of shape (3, nt) with components [N, E, Z].
    """
    # Inverse Jacobian
    coord_file = gf_db / "elements" / morton_hex / "coordinates.h5"
    with h5py.File(coord_file, "r") as f:
        xyz_elem = f["xyz"][:]
    Jinv = compute_inverse_jacobian(xyz_elem, xi, eta, gamma)

    # Strain from displacement gradients
    strain = compute_strain(displ, xi, eta, gamma, Jinv)  # (nt, 6, 3)

    # Parse CMTSOLUTION
    cmt = parse_cmtsolution(cmt_path)

    # Rotate moment tensor from spherical to Cartesian
    r = sqrt(src["x"]**2 + src["y"]**2 + src["z"]**2)
    theta = np.arccos(src["z"] / r)
    phi = np.arctan2(src["y"], src["x"])

    Mxx, Myy, Mzz, Mxy, Mxz, Myz = rotate_mt_to_xyz(
        cmt["Mrr"], cmt["Mtt"], cmt["Mpp"],
        cmt["Mrt"], cmt["Mrp"], cmt["Mtp"],
        theta, phi,
    )

    # Contract strain with moment tensor (Voigt notation)
    voigt_factor = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])
    Mx = np.array([Mxx, Myy, Mzz, Mxy, Mxz, Myz])
    mt_scale = 1.0 / (SCALE_M * factor_force_nondim)
    weights = voigt_factor * Mx * mt_scale  # shape (6,)

    gf = np.zeros((3, strain.shape[0]))
    for i in range(3):
        gf[i] = np.sum(weights[:, None] * strain[:, :, i].T, axis=0)

    # Time-integrate: Gaussian -> Heaviside STF conversion.
    # Trapezoidal rule (second-order, phase-correct). A left-endpoint Riemann
    # sum (np.cumsum * dt) is first-order and lags the true integral by half a
    # sample (dt_sub/2); on the coarse subsampled grid that is ~0.2 s, a
    # visible sub-sample phase error. The integrand is band-limited well below
    # the subsampled Nyquist, so the trapezoid is essentially exact here.
    for i in range(3):
        gf[i] = cumulative_trapezoid(gf[i], dx=dt_sub, initial=0.0)

    return gf, cmt


def apply_stf_correction(gf, dt_sub, hdur_fwd, hdur_gf):
    """Convolve GF reconstruction with corrected source time function.

    Applies Tarantola (2005) eq. 6.21 STF correction in-place.

    Parameters
    ----------
    gf : ndarray, shape (3, nt)
        GF traces to correct.
    dt_sub : float
        Subsampled time step.
    hdur_fwd : float
        Forward simulation Gaussian half-duration.
    hdur_gf : float
        GF database half-duration.
    """
    if hdur_fwd <= hdur_gf:
        print("  WARNING: hdur_fwd <= hdur_gf, no STF correction applied")
        return

    hdur_corrected = np.sqrt(hdur_fwd**2 - hdur_gf**2)
    stf_half_len = int(4.0 * hdur_corrected / dt_sub)
    stf_t = np.arange(-stf_half_len, stf_half_len + 1) * dt_sub
    stf = np.exp(-(stf_t / hdur_corrected)**2) / (np.sqrt(pi) * hdur_corrected)
    stf *= dt_sub

    print(f"  STF correction: hdur_corrected={hdur_corrected:.4f} s, "
          f"kernel={len(stf)} samples")

    for i in range(3):
        gf[i] = fftconvolve(gf[i], stf, mode="same")


def read_forward_seismograms(fwd_dir, station):
    """Read forward simulation seismograms (auto-detect channel codes).

    Returns
    -------
    traces : dict
        Keys are channel codes (e.g. 'BXN'), values are dicts with
        'time' (ndarray) and 'amp' (ndarray, float64).
    channel_prefix : str
        Detected prefix (e.g. 'BX', 'BH', 'MX').
    """
    for prefix in ["BX", "BH", "MX"]:
        test_file = fwd_dir / f"{station}.{prefix}N.sem.sac"
        if test_file.exists():
            channel_prefix = prefix
            break
    else:
        print("ERROR: Could not find forward seismogram files.")
        sys.exit(1)

    traces = {}
    for suffix in ["N", "E", "Z"]:
        comp = f"{channel_prefix}{suffix}"
        sac_file = fwd_dir / f"{station}.{comp}.sem.sac"
        tr = obspy_read(str(sac_file))[0]
        time = tr.times() + tr.stats.sac.b
        traces[comp] = {"time": time, "amp": tr.data.astype(np.float64)}

    return traces, channel_prefix


def sinc_interp(y_coarse, time_coarse, time_fine):
    """Whittaker-Shannon (finite-support sinc) interpolation onto a fine axis.

    Reconstructs the band-limited continuous signal from its uniform samples via

        recon(t) = sum_n y[n] * sinc((t - t_n) / Ts)

    where ``Ts`` is the coarse sampling interval. This is the ideal band-limited
    interpolant. Because the sum runs only over the finite set of stored samples
    (no circular wrap), it makes no periodicity assumption and is therefore safe
    at the edges of a finite, non-periodic transient -- unlike FFT-based
    ``resample``, whose periodic sinc wraps the endpoints and rings. It also has
    a flat passband all the way to the band edge, unlike a cubic spline, which
    droops over the upper part of the band.

    A non-zero endpoint level or tilt (e.g. the DC pedestal left by the
    Gaussian->Heaviside cumsum in the CMT reconstruction) is a low-frequency
    component that the truncated, slowly-decaying sinc kernel reproduces poorly
    near the trace boundaries, producing a large spurious swing at the edges. To
    avoid this, the straight line through the first and last samples is removed
    before the sinc sum and added back (exactly, since a line is reproduced
    exactly on the fine grid) afterward. This leaves an endpoint-zeroed,
    well-behaved signal for the band-limited interpolation.

    Parameters
    ----------
    y_coarse : ndarray, shape (n_coarse,)
        Samples on the coarse, uniformly spaced grid.
    time_coarse : ndarray, shape (n_coarse,)
        Coarse time axis (uniform spacing ``Ts``).
    time_fine : ndarray, shape (n_fine,)
        Target (fine) time axis.

    Returns
    -------
    ndarray, shape (n_fine,)
        ``y_coarse`` reconstructed on ``time_fine``.
    """
    # Remove the endpoint-connecting linear trend (kills the CMT DC pedestal and
    # any tilt) before interpolating; restore it on the fine grid afterward.
    slope = (y_coarse[-1] - y_coarse[0]) / (time_coarse[-1] - time_coarse[0])
    trend_coarse = y_coarse[0] + slope * (time_coarse - time_coarse[0])
    trend_fine = y_coarse[0] + slope * (time_fine - time_coarse[0])
    y_detrended = y_coarse - trend_coarse

    Ts = time_coarse[1] - time_coarse[0]
    # M[i, n] = sinc((t_fine[i] - t_coarse[n]) / Ts); shape (n_fine, n_coarse).
    # For the regional example this is ~18000 x 106, a small matmul. If a much
    # finer/longer case makes M too large, chunk time_fine into blocks.
    M = np.sinc((time_fine[:, None] - time_coarse[None, :]) / Ts)
    return M @ y_detrended + trend_fine


def align_and_filter(gf, time_gf, fwd_traces, dt, subsample_step, dt_sub,
                     f_cutoff, lowpass_period, highpass_period):
    """Lowpass-filter forward traces, upsample GF to forward rate, align, and filter.

    The comparison is carried out on the forward simulation's full-resolution
    time axis. The forward traces are kept at their native sampling rate (only
    anti-alias lowpassed to match the band-limited GF), and the subsampled GF
    reconstruction is upsampled onto the forward axis with finite-support sinc
    (Whittaker-Shannon) interpolation -- the ideal band-limited interpolant for
    the finite, non-periodic GF transient. This gives an honest full-resolution
    reconstruction error and avoids the edge ringing that FFT-based resampling
    produces on a non-periodic signal.

    Parameters
    ----------
    gf : ndarray, shape (3, nt_gf)
        GF reconstruction [N, E, Z].
    time_gf : ndarray
        GF time axis.
    fwd_traces : dict
        Forward seismogram traces from read_forward_seismograms().
    dt, subsample_step, dt_sub : float
        Timing parameters.
    f_cutoff : float
        GF database anti-alias cutoff frequency (Hz).
    lowpass_period, highpass_period : float or None
        Optional filter periods.

    Returns
    -------
    time_common : ndarray
        Common time axis after alignment.
    gf_aligned : dict
        GF traces keyed by channel code.
    fwd_aligned : dict
        Forward traces keyed by channel code.
    comps : list
        Channel code list ['BXN', 'BXE', 'BXZ'].
    comp_labels : dict
        Mapping from channel code to label name.
    """
    comps = list(fwd_traces.keys())
    comp_labels = {comps[0]: "North", comps[1]: "East", comps[2]: "Vertical"}

    # Lowpass filter forward traces to the GF band (kept at full resolution).
    # This band-limits the forward trace to the same cutoff the GF database
    # stored, so the comparison measures reconstruction error, not a band
    # mismatch. lp_fc defaults to the database anti-alias cutoff.
    lp_fc = 1.0 / lowpass_period if lowpass_period else f_cutoff
    sos = butterworth_sos(4, lp_fc, 1.0 / dt)

    fwd_filtered = {}
    for comp in comps:
        fwd_filtered[comp] = filter_with_taper(
            fwd_traces[comp]["amp"].astype(np.float64), sos
        )

    time_fwd = fwd_traces[comps[0]]["time"]

    # Align traces using the overlapping time window. The common axis is the
    # forward (full-resolution) time axis.
    t_start = max(time_fwd[0], time_gf[0])
    t_end = min(time_fwd[-1], time_gf[-1])

    mask_fwd = (time_fwd >= t_start) & (time_fwd <= t_end)
    time_common = time_fwd[mask_fwd]

    fwd_aligned = {comp: fwd_filtered[comp][mask_fwd] for comp in comps}

    # Upsample the GF reconstruction onto the forward axis with finite-support
    # sinc interpolation (band-limited, non-periodic -> no edge ringing).
    gf_aligned = {
        comps[i]: sinc_interp(gf[i], time_gf, time_common)
        for i in range(3)
    }

    # Lowpass the reconstructed GF to the same band as the forward trace. The GF
    # is already band-limited by construction, but this cleans any residual
    # band-edge wiggle and keeps both traces in an identical band.
    for comp in comps:
        gf_aligned[comp] = filter_with_taper(gf_aligned[comp], sos)

    # Optional highpass filter on both, at the forward sampling rate.
    if highpass_period is not None:
        hp_sos = butterworth_highpass_sos(4, 1.0 / highpass_period, 1.0 / dt)
        for comp in comps:
            gf_aligned[comp] = filter_with_taper(gf_aligned[comp], hp_sos)
            fwd_aligned[comp] = filter_with_taper(fwd_aligned[comp], hp_sos)

    return time_common, gf_aligned, fwd_aligned, comps, comp_labels


def compute_residuals(gf_traces, fwd_traces, comps):
    """Compute relative RMS error per component.

    Returns dict of {comp: rel_err}.
    """
    residuals = {}
    for comp in comps:
        fwd = fwd_traces[comp]
        gf = gf_traces[comp]
        rms_err = np.sqrt(np.sum((fwd - gf) ** 2))
        rms_ref = np.sqrt(np.sum(fwd ** 2))
        residuals[comp] = rms_err / rms_ref if rms_ref > 0 else float("inf")
    return residuals


def plot_comparison(time_common, gf_traces, fwd_traces, comps, comp_labels,
                    residuals, station, source_type, morton_hex, xi, eta, gamma,
                    highpass_period, lowpass_period, output_plot):
    """Plot GF vs forward comparison and save."""
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 7.5), sharex=True)

    for ax, comp in zip(axes, comps):
        ax.plot(time_common, fwd_traces[comp], "b-", linewidth=0.8,
                label="Forward (filtered)")
        ax.plot(time_common, gf_traces[comp], "r--", linewidth=0.8,
                label="GF reconstructed")
        ax.set_ylabel(f"{comp} ({comp_labels[comp]})")
        ax.set_title(f"{comp} — relative error: {residuals[comp]:.4e}")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time (s)")

    filter_info = ""
    if lowpass_period:
        filter_info += f", lowpass T<{lowpass_period:.0f}s"
    if highpass_period:
        filter_info += f", highpass T>{highpass_period:.0f}s"

    fig.suptitle(
        f"GF Cross-Validation ({source_type}): {station}{filter_info}\n"
        f"Source at element {morton_hex}\n"
        f"(xi={xi:.4f}, eta={eta:.4f}, gamma={gamma:.4f})",
        fontsize=11,
    )
    plt.tight_layout()
    plt.savefig(output_plot, dpi=150)
    plt.close()


# ---------------------------------------------------------------------------
# Main cross-validation logic
# ---------------------------------------------------------------------------

def cross_validate(
    gf_db_path,
    forward_output,
    station,
    force=False,
    cmtsolution=None,
    highpass_period=None,
    lowpass_period=None,
    output=None,
    forward_solver_output=None,
    gf_solver_output=None,
):
    """Cross-validate a GF database against a forward simulation.

    Parameters
    ----------
    gf_db_path : Path
        Path to the GF database directory.
    forward_output : Path
        Path to the forward simulation OUTPUT_FILES directory.
    station : str
        Station name as NET.STA.
    force : bool
        If True, use force source workflow; otherwise use CMT.
    cmtsolution : Path or None
        Path to CMTSOLUTION file. Auto-detected if None and force=False.
    highpass_period : float or None
        Highpass filter period in seconds.
    lowpass_period : float or None
        Override lowpass filter period in seconds.
    output : Path or None
        Output plot path. Auto-generated if None.
    forward_solver_output : Path or None
        Path to forward output_solver.txt. Auto-detected if None.
    gf_solver_output : Path or None
        Path to GF output_solver.txt for t0_gf. Auto-detected if None.

    Returns
    -------
    output_plot : Path
        Path to the saved plot.
    """
    gf_db = Path(gf_db_path)
    fwd_dir = Path(forward_output)

    solver_output = forward_solver_output or fwd_dir / "output_solver.txt"
    output_plot = Path(output) if output else fwd_dir.parent / "gf_cross_validation.png"

    source_type = "Force" if force else "CMT"
    print(f"\nSource type: {source_type}")

    # 1. Parse forward simulation output
    print("Parsing forward solver output...")
    src = parse_solver_output(solver_output)
    print(f"  Source (x,y,z): ({src['x']:.6f}, {src['y']:.6f}, {src['z']:.6f})")

    # 2. Read mesh info and station metadata
    mesh = read_mesh_info(gf_db)
    dt = mesh["dt"]
    subsample_step = mesh["subsample_step"]
    nt_sub = mesh["nt_sub"]
    dt_sub = dt * subsample_step

    sta_meta = read_station_meta(gf_db, station)
    print(f"  f_cutoff={sta_meta['f_cutoff']:.4f} Hz, hdur={sta_meta['hdur']:.6f} s")

    # 3. Find containing element
    print("Finding containing element...")
    source_xyz = np.array([src["x"], src["y"], src["z"]])
    morton_hex, xi, eta, gamma, dist = find_containing_element(gf_db, source_xyz)
    print(f"  Element: {morton_hex} (xi={xi:.4f}, eta={eta:.4f}, gamma={gamma:.4f})")

    # 4. Read GF displacement
    gf_file = gf_db / "elements" / morton_hex / f"{station}.h5"
    with h5py.File(gf_file, "r") as f:
        displ = f["displacement"][:]

    # 5. Reconstruct seismograms
    if force:
        gf = reconstruct_force(displ, xi, eta, gamma, src)
    else:
        # Auto-detect CMTSOLUTION if not provided
        if cmtsolution is None:
            candidates = [
                fwd_dir.parent / "DATA" / "CMTSOLUTION",
                fwd_dir / ".." / "DATA" / "CMTSOLUTION",
            ]
            for c in candidates:
                if c.resolve().exists():
                    cmtsolution = c.resolve()
                    break
            if cmtsolution is None:
                print("ERROR: Could not find CMTSOLUTION file. "
                      "Use --cmtsolution to specify.")
                sys.exit(1)

        gf, cmt = reconstruct_cmt(
            displ, xi, eta, gamma, gf_db, morton_hex,
            src, cmtsolution, sta_meta["factor_force_nondim"], dt_sub,
        )

    # 6. STF correction
    if force:
        fwd_text = Path(solver_output).read_text()
        m = re.search(r"Gaussian half duration:\s*([\d.Ee+-]+)", fwd_text)
        if m:
            hdur_fwd = float(m.group(1))
        else:
            m2 = re.search(r"half duration:\s*([\d.Ee+-]+)\s*seconds", fwd_text)
            hdur_fwd = (float(m2.group(1)) / SOURCE_DECAY_MIMIC_TRIANGLE
                        if m2 else sta_meta["hdur"])
    else:
        hdur_fwd = cmt["half_duration"] / SOURCE_DECAY_MIMIC_TRIANGLE

    apply_stf_correction(gf, dt_sub, hdur_fwd, sta_meta["hdur"])

    # 7. Read forward seismograms
    fwd_traces, channel_prefix = read_forward_seismograms(fwd_dir, station)

    # 8. Construct GF time axis.
    # t0_gf is the time of the first time step relative to the source origin
    # (simulation time is (it-1)*dt - t0). It sets the absolute time origin of
    # the reconstructed trace, so an error here is a rigid time shift between
    # the GF and forward traces.
    t0_gf = mesh.get("t0")

    time_gf = np.array([(isnap * subsample_step - 1) * dt - t0_gf
                         for isnap in range(1, nt_sub + 1)])

    # 9. Align, filter, compare
    time_common, gf_aligned, fwd_aligned, comps, comp_labels = align_and_filter(
        gf, time_gf, fwd_traces, dt, subsample_step, dt_sub,
        sta_meta["f_cutoff"], lowpass_period, highpass_period,
    )

    # 10. Compute and print residuals
    residuals = compute_residuals(gf_aligned, fwd_aligned, comps)
    print("Residuals:")
    for comp in comps:
        print(f"  {comp}: {residuals[comp]:.6e}")

    # 11. Plot
    output_plot.parent.mkdir(parents=True, exist_ok=True)
    plot_comparison(
        time_common, gf_aligned, fwd_aligned, comps, comp_labels,
        residuals, station, source_type, morton_hex, xi, eta, gamma,
        highpass_period, lowpass_period, output_plot,
    )
    print(f"Saved {output_plot}")

    return output_plot


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

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
        "--forward-solver-output", type=Path, default=None,
        help="Path to forward output_solver.txt "
             "(default: <forward_output>/output_solver.txt)",
    )
    parser.add_argument(
        "--station", type=str, default="IU.SJG",
        help="Station name as NET.STA (default: IU.SJG)",
    )
    parser.add_argument(
        "--gf-solver-output", type=Path, default=None,
        help="Path to GF output_solver.txt (to get t0_gf; default: auto-detect)",
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Output plot path (default: gf_cross_validation.png next to forward output)",
    )
    parser.add_argument(
        "--force", action="store_true", default=False,
        help="Use force source workflow (default: moment tensor / CMT workflow).",
    )
    parser.add_argument(
        "--cmtsolution", type=Path, default=None,
        help="Path to CMTSOLUTION file "
             "(default: auto-detect from forward output directory).",
    )
    parser.add_argument(
        "--highpass-period", type=float, default=None,
        help="Apply zero-phase Butterworth highpass filter at this period (seconds).",
    )
    parser.add_argument(
        "--lowpass-period", type=float, default=None,
        help="Override the lowpass filter period (seconds). "
             "Default uses the GF database anti-alias cutoff frequency.",
    )
    args = parser.parse_args()

    cross_validate(
        gf_db_path=args.gf_db_path,
        forward_output=args.forward_output,
        station=args.station,
        force=args.force,
        cmtsolution=args.cmtsolution,
        highpass_period=args.highpass_period,
        lowpass_period=args.lowpass_period,
        output=Path(args.output) if args.output else None,
        forward_solver_output=args.forward_solver_output,
        gf_solver_output=args.gf_solver_output,
    )


if __name__ == "__main__":
    main()
