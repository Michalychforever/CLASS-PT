#!/usr/bin/env python
"""
Generate reference data for CLASS-PT regression tests.

Run from the CLASS-PT root directory:
    python tests/generate_reference_data.py

This creates tests/reference_data.npz containing saved output values
for all configurations used in test_classpt.py.
"""
import sys, os

BASEDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASEDIR, 'python'))
os.chdir(BASEDIR)

from classy import Class
import numpy as np

# Identical to test_classpt.py constants
COSMO_PARAMS = {
    'A_s': 2.089e-9, 'n_s': 0.9649, 'tau_reio': 0.052,
    'omega_b': 0.02237, 'omega_cdm': 0.12, 'h': 0.6736,
    'YHe': 0.2425, 'N_ur': 2.0328, 'N_ncdm': 1, 'm_ncdm': 0.06,
}
Z_PK = 0.61
B1, B2, BG2, BGAMMA3 = 2.0, -1.0, 0.1, -0.1
CS, CS0, CS2, CS4, PSHOT, B4 = 1.0, 5.0, 15.0, -5.0, 5e3, 100.0
PSHOT_NBAR, A0_NBAR, A2_NBAR = 5e3, 0., 0.

def make_kvecs(h):
    kvec = np.logspace(-3, np.log10(3), 300)
    return kvec, kvec * h

# Reference k-indices: 9 points spanning the k-grid
ref_idx = np.array([10, 50, 100, 130, 150, 180, 200, 230, 260])

# ============================================================================
# Create all instances in EXACTLY the same order as _create_shared_instances()
# in test_classpt.py.  This is critical because CLASS-PT has C-level global
# state that affects fNL values depending on instance creation order.
# ============================================================================

print("Creating instances (same order as test_classpt.py)...")

# 1. Full config (RSD + IR + Bias + AP + fNL)
M_full = Class()
M_full.set(COSMO_PARAMS)
M_full.set({
    'output': 'mPk', 'non linear': 'PT',
    'IR resummation': 'Yes', 'Bias tracers': 'Yes',
    'cb': 'Yes', 'RSD': 'Yes', 'AP': 'Yes', 'Omfid': '0.31',
    'z_pk': Z_PK, 'PNG': 'Yes',
})
M_full.compute()
print("  [1/6] full config")

# 2. No-IR, no-bias, no-RSD
M_noIR = Class()
M_noIR.set(COSMO_PARAMS)
M_noIR.set({
    'output': 'mPk', 'non linear': 'PT',
    'IR resummation': 'No', 'Bias tracers': 'No',
    'cb': 'Yes', 'z_pk': Z_PK,
})
M_noIR.compute()
print("  [2/6] noIR config")

# 3. Bias, no-RSD (real-space biased tracers)
M_bias = Class()
M_bias.set(COSMO_PARAMS)
M_bias.set({
    'output': 'mPk', 'non linear': 'PT',
    'IR resummation': 'No', 'Bias tracers': 'Yes',
    'RSD': 'No', 'cb': 'Yes', 'z_pk': Z_PK,
})
M_bias.compute()
print("  [3/6] bias_real config")

# 4. RSD + IR + Bias, no AP
M_noAP = Class()
M_noAP.set(COSMO_PARAMS)
M_noAP.set({
    'output': 'mPk', 'non linear': 'PT',
    'IR resummation': 'Yes', 'Bias tracers': 'Yes',
    'cb': 'Yes', 'RSD': 'Yes', 'AP': 'No',
    'z_pk': Z_PK,
})
M_noAP.compute()
print("  [4/6] noAP config")

# 5. CMB lensing with PT
M_cmb = Class()
M_cmb.set(COSMO_PARAMS)
M_cmb.set({
    'output': 'mPk,lCl,tCl,pCl', 'lensing': 'Yes', 'l_switch_limber': 9,
    'non linear': 'PT', 'cb': 'No',
    'IR resummation': 'No', 'Bias tracers': 'No', 'RSD': 'No',
    'z_pk': Z_PK,
})
M_cmb.compute()
print("  [5/6] cmb_pt config")

# 6. CMB lensing linear (no non-linear)
M_cmb_lin = Class()
M_cmb_lin.set(COSMO_PARAMS)
M_cmb_lin.set({
    'output': 'mPk,lCl,tCl,pCl', 'lensing': 'Yes', 'l_switch_limber': 9,
    'non linear': 'none', 'z_pk': Z_PK,
})
M_cmb_lin.compute()
print("  [6/7] cmb_lin config")

# 7. External P(k) with RSD + IR + Bias
ext_file = os.path.join(BASEDIR, 'tests', 'external_Pk_test.dat')
if os.path.exists(ext_file):
    COSMO_EXT = {k: v for k, v in COSMO_PARAMS.items() if k not in ('A_s', 'n_s')}
    M_ext = Class()
    M_ext.set(COSMO_EXT)
    M_ext.set({
        'output': 'mPk', 'non linear': 'PT',
        'IR resummation': 'Yes', 'Bias tracers': 'Yes',
        'cb': 'Yes', 'RSD': 'Yes', 'AP': 'No', 'z_pk': Z_PK,
        'P_k_ini type': 'external_Pk',
        'command': 'cat ' + os.path.abspath(ext_file),
    })
    M_ext.compute()
    print("  [7/7] ext_pk config")
else:
    M_ext = None
    print("  [7/7] ext_pk SKIPPED (file not found)")
print("All instances created.\n")

ref = {}
ref['ref_idx'] = ref_idx

# ============================================================================
# Full config reference data
# ============================================================================
h = M_full.h()
kvec, khvec = make_kvecs(h)
M_full.initialize_output(khvec, Z_PK, len(khvec))
pk_mult = M_full.get_pk_mult(khvec, Z_PK, len(khvec))
fz = M_full.scale_independent_growth_factor_f(Z_PK)

ref['ref_kvec'] = kvec[ref_idx]
ref['h'] = np.float64(h)
ref['sigma8'] = np.float64(M_full.sigma8())
ref['Omega_m'] = np.float64(M_full.Omega_m())
ref['fz'] = np.float64(fz)
ref['full_pk_mult'] = pk_mult[:, ref_idx]

# Real-space convenience methods
ref['pk_mm_real'] = M_full.pk_mm_real(CS)[ref_idx]
ref['pk_gm_real'] = M_full.pk_gm_real(B1, B2, BG2, BGAMMA3, CS, CS0)[ref_idx]
ref['pk_gg_real'] = M_full.pk_gg_real(B1, B2, BG2, BGAMMA3, CS, CS0, PSHOT)[ref_idx]

# Matter multipoles
ref['pk_mm_l0'] = M_full.pk_mm_l0(CS0)[ref_idx]
ref['pk_mm_l2'] = M_full.pk_mm_l2(CS2)[ref_idx]
ref['pk_mm_l4'] = M_full.pk_mm_l4(CS4)[ref_idx]

# Galaxy multipoles
ref['pk_gg_l0'] = M_full.pk_gg_l0(B1, B2, BG2, BGAMMA3, CS0, PSHOT_NBAR, A0_NBAR, A2_NBAR, B4)[ref_idx]
ref['pk_gg_l2'] = M_full.pk_gg_l2(B1, B2, BG2, BGAMMA3, CS2, A2_NBAR, B4)[ref_idx]
ref['pk_gg_l4'] = M_full.pk_gg_l4(B1, B2, BG2, BGAMMA3, CS4, B4)[ref_idx]

# fNL equilateral
ref['pk_mm_fNL_real'] = M_full.pk_mm_fNL_real()[ref_idx]
ref['pk_gm_fNL_real'] = M_full.pk_gm_fNL_real(B1, B2, BG2)[ref_idx]
ref['pk_gg_fNL_real'] = M_full.pk_gg_fNL_real(B1, B2, BG2)[ref_idx]
ref['pk_mm_fNL_l0'] = M_full.pk_mm_fNL_l0()[ref_idx]
ref['pk_mm_fNL_l2'] = M_full.pk_mm_fNL_l2()[ref_idx]
ref['pk_mm_fNL_l4'] = M_full.pk_mm_fNL_l4()[ref_idx]
ref['pk_gg_fNL_l0'] = M_full.pk_gg_fNL_l0(B1, B2, BG2)[ref_idx]
ref['pk_gg_fNL_l2'] = M_full.pk_gg_fNL_l2(B1, B2, BG2)[ref_idx]
ref['pk_gg_fNL_l4'] = M_full.pk_gg_fNL_l4(B1, B2, BG2)[ref_idx]

# fNL orthogonal
ref['pk_mm_fNL_real_ortho'] = M_full.pk_mm_fNL_real_ortho()[ref_idx]
ref['pk_gm_fNL_real_ortho'] = M_full.pk_gm_fNL_real_ortho(B1, B2, BG2)[ref_idx]
ref['pk_gg_fNL_real_ortho'] = M_full.pk_gg_fNL_real_ortho(B1, B2, BG2)[ref_idx]
ref['pk_mm_fNL_l0_ortho'] = M_full.pk_mm_fNL_l0_ortho()[ref_idx]
ref['pk_mm_fNL_l2_ortho'] = M_full.pk_mm_fNL_l2_ortho()[ref_idx]
ref['pk_mm_fNL_l4_ortho'] = M_full.pk_mm_fNL_l4_ortho()[ref_idx]
ref['pk_gg_fNL_l0_ortho'] = M_full.pk_gg_fNL_l0_ortho(B1, B2, BG2)[ref_idx]
ref['pk_gg_fNL_l2_ortho'] = M_full.pk_gg_fNL_l2_ortho(B1, B2, BG2)[ref_idx]
ref['pk_gg_fNL_l4_ortho'] = M_full.pk_gg_fNL_l4_ortho(B1, B2, BG2)[ref_idx]
print(f"Full config: {sum(1 for k in ref if not k.startswith(('noIR','noAP','bias','cmb','nw_','alpha')))} keys")

# ============================================================================
# noIR config reference data
# ============================================================================
h2 = M_noIR.h()
kvec2, khvec2 = make_kvecs(h2)
M_noIR.initialize_output(khvec2, Z_PK, len(khvec2))
pk_mult2 = M_noIR.get_pk_mult(khvec2, Z_PK, len(khvec2))

ref['noIR_pk_mm_real'] = M_noIR.pk_mm_real(CS)[ref_idx]
# Save valid pk_mult components (skip NaN/sentinel-value components)
known_nan = {49, 73, 85}
valid_indices = []
for i in range(96):
    if i in known_nan:
        continue
    vals = pk_mult2[i, ref_idx]
    if np.all(np.isfinite(vals)) and not np.all(np.abs(vals) > 999000):
        valid_indices.append(i)
        ref[f'noIR_pk_mult_{i}'] = vals
ref['noIR_valid_indices'] = np.array(valid_indices)
print(f"noIR config: pk_mm_real + {len(valid_indices)} pk_mult components")

# ============================================================================
# noAP config reference data
# ============================================================================
h4 = M_noAP.h()
kvec4, khvec4 = make_kvecs(h4)
M_noAP.initialize_output(khvec4, Z_PK, len(khvec4))
pk_mult4 = M_noAP.get_pk_mult(khvec4, Z_PK, len(khvec4))

ref['noAP_pk_mult'] = pk_mult4[:, ref_idx]
ref['noAP_pk_mm_real'] = M_noAP.pk_mm_real(CS)[ref_idx]
ref['noAP_pk_gm_real'] = M_noAP.pk_gm_real(B1, B2, BG2, BGAMMA3, CS, CS0)[ref_idx]
ref['noAP_pk_gg_real'] = M_noAP.pk_gg_real(B1, B2, BG2, BGAMMA3, CS, CS0, PSHOT)[ref_idx]
ref['noAP_pk_mm_l0'] = M_noAP.pk_mm_l0(CS0)[ref_idx]
ref['noAP_pk_mm_l2'] = M_noAP.pk_mm_l2(CS2)[ref_idx]
ref['noAP_pk_mm_l4'] = M_noAP.pk_mm_l4(CS4)[ref_idx]
ref['noAP_pk_gg_l0'] = M_noAP.pk_gg_l0(B1, B2, BG2, BGAMMA3, CS0, PSHOT_NBAR, A0_NBAR, A2_NBAR, B4)[ref_idx]
ref['noAP_pk_gg_l2'] = M_noAP.pk_gg_l2(B1, B2, BG2, BGAMMA3, CS2, A2_NBAR, B4)[ref_idx]
ref['noAP_pk_gg_l4'] = M_noAP.pk_gg_l4(B1, B2, BG2, BGAMMA3, CS4, B4)[ref_idx]
print(f"noAP config: 10 convenience methods + pk_mult")

# ============================================================================
# bias_real config reference data
# ============================================================================
h3 = M_bias.h()
kvec3, khvec3 = make_kvecs(h3)
M_bias.initialize_output(khvec3, Z_PK, len(khvec3))

ref['bias_real_pk_mm_real'] = M_bias.pk_mm_real(CS)[ref_idx]
ref['bias_real_pk_gm_real'] = M_bias.pk_gm_real(B1, B2, BG2, BGAMMA3, CS, CS0)[ref_idx]
ref['bias_real_pk_gg_real'] = M_bias.pk_gg_real(B1, B2, BG2, BGAMMA3, CS, CS0, PSHOT)[ref_idx]
print(f"bias_real config: 3 convenience methods")

# ============================================================================
# CMB lensing reference data
# ============================================================================
cl = M_cmb.lensed_cl(2500)
ell_ref = np.array([10, 50, 100, 200, 500, 1000, 1500, 2000, 2500])
ref['cmb_ell_ref'] = ell_ref
ref['cmb_lensed_tt'] = cl['tt'][ell_ref]
ref['cmb_lensed_ee'] = cl['ee'][ell_ref]
ref['cmb_lensed_pp'] = cl['pp'][ell_ref]

h_cmb = M_cmb.h()
k_test = np.array([0.01, 0.05, 0.1, 0.2, 0.5]) * h_cmb
pk_pt = []
for k in k_test:
    r = M_cmb.pk(k, Z_PK)
    pk_pt.append(r[0] if isinstance(r, (list, tuple)) else r)
ref['cmb_k_test'] = k_test
ref['cmb_pk_pt'] = np.array(pk_pt)
print(f"CMB config: 3 Cl spectra + {len(k_test)} P(k) values")

# ============================================================================
# No-wiggle and alpha_rs reference data (using noAP instance)
# ============================================================================
k_nw = np.array([0.05, 0.1, 0.15, 0.2, 0.25]) * h4
pk_w, pk_nw = [], []
for k in k_nw:
    rw = M_noAP.pk(k, Z_PK, no_wiggle=False)
    rnw = M_noAP.pk(k, Z_PK, no_wiggle=True)
    pk_w.append(rw[0] if isinstance(rw, (list, tuple)) else rw)
    pk_nw.append(rnw[0] if isinstance(rnw, (list, tuple)) else rnw)
ref['nw_k_test'] = k_nw
ref['nw_pk_wiggle'] = np.array(pk_w)
ref['nw_pk_nowiggle'] = np.array(pk_nw)
alphas = [0.95, 1.0, 1.05]
k_alpha = 0.1 * h4
pk_alpha = []
for a in alphas:
    r = M_noAP.pk(k_alpha, Z_PK, alpha_rs=a)
    pk_alpha.append(r[0] if isinstance(r, (list, tuple)) else r)
ref['alpha_k_test'] = np.float64(k_alpha)
ref['alpha_values'] = np.array(alphas)
ref['alpha_pk'] = np.array(pk_alpha)
print(f"No-wiggle/alpha_rs: {len(k_nw)} k-points, {len(alphas)} alpha values")

# ============================================================================
# Save
# ============================================================================
outfile = os.path.join(BASEDIR, 'tests', 'reference_data.npz')
np.savez(outfile, **ref)
print(f"\nSaved {len(ref)} keys to {outfile}")

# Cleanup
instances = [M_full, M_noIR, M_bias, M_noAP, M_cmb, M_cmb_lin]
if M_ext is not None:
    instances.append(M_ext)
for M in instances:
    M.struct_cleanup()
    M.empty()
print("Done.")
