#!/usr/bin/env python
"""Bundle all M12 binary matrices into a single file for fast loading.

Creates two files:
  M12_bundle_N{N}.bin   - all regular M12 matrices concatenated
  M12_bundle_N{N}_ortho.bin - all orthogonal M12 matrices concatenated

The order matches the m12_pairs array in nonlinear_pt.c.
Each matrix has size (N+1)*(N+2) doubles.
"""
import numpy as np
import os
import sys

CLASSDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Must match m12_pairs[] in nonlinear_pt.c exactly
M12_PAIRS = [
    ("", "matter"),
    ("/matter_multipoles", "matter_multipoles-vv0_f2.txt"),
    ("/matter_multipoles", "matter_multipoles-vv0_f3.txt"),
    ("/matter_multipoles", "matter_multipoles-vd0_f1.txt"),
    ("/matter_multipoles", "matter_multipoles-vd0_f2.txt"),
    ("/matter_multipoles", "matter_multipoles-dd0_f0.txt"),
    ("/matter_multipoles", "matter_multipoles-dd0_f1.txt"),
    ("/matter_multipoles", "matter_multipoles-vv2_f3.txt"),
    ("/matter_multipoles", "matter_multipoles-vd2_f2.txt"),
    ("/matter_multipoles", "matter_multipoles-vv4_f3.txt"),
    ("/matter_multipoles", "matter_multipoles-vd4_f2.txt"),
    ("/matter_mu_powers", "matter_mu_powers-vd2_f1.txt"),
    ("/matter_mu_powers", "matter_mu_powers-vd2_f2.txt"),
    ("/matter_mu_powers", "matter_mu_powers-dd2_f1.txt"),
    ("/matter_mu_powers", "matter_mu_powers-vv4_f2.txt"),
    ("/matter_mu_powers", "matter_mu_powers-vd4_f2.txt"),
    ("/matter_mu_powers", "matter_mu_powers-vv6_f3.txt"),
    ("/bias_real_space", "bias_real_space-b2.txt"),
    ("/bias_real_space", "bias_real_space-G2.txt"),
    ("/bias_multipoles", "bias_multipoles-b2_vv0_f1.txt"),
    ("/bias_multipoles", "bias_multipoles-bG2_vv0_f1.txt"),
]

def bundle(N):
    Nstr = f"N{N}"
    m12_size = (N + 1) * (N + 2)

    reg_data = []
    ortho_data = []

    for subdir, filename in M12_PAIRS:
        # Regular
        reg_path = f"{CLASSDIR}/pt_matrices/compute_matrices_python{subdir}/M12oneline_{Nstr}-{filename}.dat.bin"
        if not os.path.exists(reg_path):
            print(f"Missing: {reg_path}")
            sys.exit(1)
        d = np.fromfile(reg_path, dtype=np.float64)
        assert len(d) == m12_size, f"{reg_path}: expected {m12_size}, got {len(d)}"
        reg_data.append(d)

        # Orthogonal
        ortho_path = f"{CLASSDIR}/ORTHOGONAL_NG_matrices/computation_of_matrices{subdir}/M12oneline_{Nstr}-{filename}-orthogonal.dat.bin"
        if not os.path.exists(ortho_path):
            print(f"Missing: {ortho_path}")
            sys.exit(1)
        d = np.fromfile(ortho_path, dtype=np.float64)
        assert len(d) == m12_size, f"{ortho_path}: expected {m12_size}, got {len(d)}"
        ortho_data.append(d)

    # Write bundled files
    reg_bundle = np.concatenate(reg_data)
    ortho_bundle = np.concatenate(ortho_data)

    reg_out = f"{CLASSDIR}/pt_matrices/M12_bundle_{Nstr}.bin"
    ortho_out = f"{CLASSDIR}/pt_matrices/M12_bundle_{Nstr}_ortho.bin"

    reg_bundle.tofile(reg_out)
    ortho_bundle.tofile(ortho_out)

    n = len(M12_PAIRS)
    print(f"Bundled {n} matrices (each {m12_size} doubles = {m12_size*8/1024:.1f}KB)")
    print(f"  Regular: {reg_out} ({os.path.getsize(reg_out)/1024:.1f}KB)")
    print(f"  Ortho:   {ortho_out} ({os.path.getsize(ortho_out)/1024:.1f}KB)")

if __name__ == "__main__":
    for N in [128, 256, 512]:
        try:
            print(f"\n=== N={N} ===")
            bundle(N)
        except (FileNotFoundError, SystemExit) as e:
            print(f"  Skipped N={N}: {e}")
