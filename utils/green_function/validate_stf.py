"""
Cross-validation of the Fortran Green function STF against Python (scipy).

Reads the gf_stf.txt output from the solver and compares the filtered STF
against an independently computed Butterworth-filtered Gaussian using scipy.
"""

import numpy as np
from scipy.signal import butter, sosfiltfilt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os


def parse_gf_stf(filename):
    """Parse the gf_stf.txt file and extract parameters + data."""
    params = {}
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                line = line[1:].strip()
                if '=' in line:
                    key, val = line.split('=', 1)
                    key = key.strip()
                    val = val.strip()
                    try:
                        if '.' in val or 'E' in val.upper():
                            params[key] = float(val)
                        else:
                            params[key] = int(val)
                    except ValueError:
                        params[key] = val

    # Fortran ES format can produce values like '1.234-260' (missing 'E')
    # when the exponent is 3+ digits. Fix these before parsing.
    import re
    lines = []
    with open(filename) as f2:
        for line in f2:
            if not line.startswith('#'):
                # fix missing 'E' in exponents: match digit followed by +/- and digits
                line = re.sub(r'(\d)([-+])(\d{2,})', r'\1E\2\3', line)
                lines.append(line)
    from io import StringIO
    data = np.loadtxt(StringIO(''.join(lines)))
    return params, data


def gauss_stf(t, hdur):
    """Specfem-style Gaussian STF."""
    return np.exp(-(t ** 2) / (hdur ** 2)) / (np.sqrt(np.pi) * hdur)


def main():
    stf_file = 'db_gen/OUTPUT_FILES/gf_stf.txt'
    if not os.path.exists(stf_file):
        print(f'Error: {stf_file} not found')
        print('Run the solver in DB_TESTING/db_gen first.')
        sys.exit(1)

    # parse Fortran output
    params, data = parse_gf_stf(stf_file)
    t = data[:, 0]
    stf_unfiltered_f90 = data[:, 1]
    stf_filtered_f90 = data[:, 2]

    hdur = params['hdur']
    f_cutoff = params['f_cutoff']
    dt = params['dt']
    nstep = params['nstep']
    t0 = params['t0']
    order = params['filter_order']

    print(f'Parameters from Fortran:')
    print(f'  hdur = {hdur:.6f} s')
    print(f'  f_cutoff = {f_cutoff:.6f} Hz')
    print(f'  dt = {dt:.6f} s')
    print(f'  nstep = {nstep}')
    print(f'  t0 = {t0:.6f} s')
    print(f'  filter_order = {order}')
    print()

    # compute Python reference
    fs = 1.0 / dt
    stf_unfiltered_py = gauss_stf(t, hdur)

    # Butterworth lowpass via scipy
    sos_py = butter(order, f_cutoff, btype='low', fs=fs, output='sos')
    stf_filtered_py = sosfiltfilt(sos_py, stf_unfiltered_py)

    # compare
    diff_unfiltered = np.max(np.abs(stf_unfiltered_f90 - stf_unfiltered_py))
    diff_filtered = np.max(np.abs(stf_filtered_f90 - stf_filtered_py))
    peak_f90 = np.max(np.abs(stf_filtered_f90))
    peak_py = np.max(np.abs(stf_filtered_py))
    rel_diff = diff_filtered / peak_py if peak_py > 0 else float('inf')

    print(f'Comparison:')
    print(f'  max |diff| unfiltered: {diff_unfiltered:.2e}')
    print(f'  max |diff| filtered:   {diff_filtered:.2e}')
    print(f'  peak Fortran:          {peak_f90:.6f}')
    print(f'  peak Python:           {peak_py:.6f}')
    print(f'  relative diff:         {rel_diff:.2e}')
    print()

    if rel_diff < 1e-6:
        print('PASS: Fortran and Python STFs match to within 1e-6 relative error.')
    elif rel_diff < 1e-3:
        print('CLOSE: Fortran and Python STFs differ slightly (check plot).')
    else:
        print('FAIL: Significant difference between Fortran and Python STFs.')

    # plot
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    axes[0].plot(t, stf_unfiltered_f90, 'b-', label='Fortran', linewidth=1.5)
    axes[0].plot(t, stf_unfiltered_py, 'r--', label='Python', linewidth=1.5)
    axes[0].set_ylabel('Amplitude')
    axes[0].set_title('Unfiltered Gaussian STF')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, stf_filtered_f90, 'b-', label='Fortran', linewidth=1.5)
    axes[1].plot(t, stf_filtered_py, 'r--', label='Python', linewidth=1.5)
    axes[1].set_ylabel('Amplitude')
    axes[1].set_title(f'Butterworth-filtered STF (order={order}, fc={f_cutoff:.3f} Hz)')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t, stf_filtered_f90 - stf_filtered_py, 'k-', linewidth=1.0)
    axes[2].set_ylabel('Difference')
    axes[2].set_xlabel('Time (s)')
    axes[2].set_title('Fortran - Python')
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('stf_cross_validation.png', dpi=150)
    print(f'\nPlot saved to stf_cross_validation.png')


if __name__ == '__main__':
    main()
