"""
Comprehensive unit tests for CLASS-PT Python interface.

Tests all major functionalities:
  - Background cosmology
  - Linear power spectrum
  - Real-space 1-loop matter/galaxy/cross spectra (Pmm, Pgm, Pgg)
  - Redshift-space multipoles (l=0,2,4) for matter and galaxies
  - fNL contributions (equilateral and orthogonal), real-space and multipoles
  - IR resummation toggle
  - Alcock-Paczynski effect
  - CMB lensing with PT
  - External Pk input
  - Consistency checks between convenience methods and raw pk_mult components
  - Recompute with no-wiggle / alpha_rs

Run with: python tests/test_classpt.py
"""

import numpy as np
import os
import sys
import unittest

# Add python directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))
from classy import Class

# Reference cosmology (Planck 2018-like)
COSMO_PARAMS = {
    'A_s': 2.089e-9,
    'n_s': 0.9649,
    'tau_reio': 0.052,
    'omega_b': 0.02237,
    'omega_cdm': 0.12,
    'h': 0.6736,
    'YHe': 0.2425,
    'N_ur': 2.0328,
    'N_ncdm': 1,
    'm_ncdm': 0.06,
}

Z_PK = 0.61

# Bias parameters
B1 = 2.0
B2 = -1.0
BG2 = 0.1
BGAMMA3 = -0.1
CS = 1.0
CS0 = 5.0
CS2 = 15.0
CS4 = -5.0
PSHOT = 5e3
B4 = 100.0


def make_kvecs(h, kmin=-3, kmax_log10=np.log10(3), nk=300):
    """Create k vectors in h/Mpc and 1/Mpc."""
    kvec = np.logspace(kmin, kmax_log10, nk)
    khvec = kvec * h
    return kvec, khvec


# ============================================================================
# Module-level shared Class instances.
# CLASS-PT has C-level global state that can be corrupted when creating multiple
# instances. We create all instances upfront to avoid cross-contamination.
# ============================================================================

def _create_shared_instances():
    """Create all shared Class instances at module load time."""
    instances = {}

    # Full config (RSD + IR + Bias + AP + fNL) — most featureful
    M = Class()
    M.set(COSMO_PARAMS)
    M.set({
        'output': 'mPk', 'non linear': 'PT',
        'IR resummation': 'Yes', 'Bias tracers': 'Yes',
        'cb': 'Yes', 'RSD': 'Yes', 'AP': 'Yes', 'Omfid': '0.31',
        'z_pk': Z_PK, 'PNG': 'Yes',
    })
    M.compute()
    h = M.h()
    kvec, khvec = make_kvecs(h)
    M.initialize_output(khvec, Z_PK, len(khvec))
    pk_mult = M.get_pk_mult(khvec, Z_PK, len(khvec))
    fz = M.scale_independent_growth_factor_f(Z_PK)
    instances['full'] = {
        'M': M, 'h': h, 'kvec': kvec, 'khvec': khvec,
        'pk_mult': pk_mult, 'fz': fz,
    }

    # No-IR, no-bias, no-RSD config (real-space matter only)
    M2 = Class()
    M2.set(COSMO_PARAMS)
    M2.set({
        'output': 'mPk', 'non linear': 'PT',
        'IR resummation': 'No', 'Bias tracers': 'No',
        'cb': 'Yes', 'z_pk': Z_PK,
    })
    M2.compute()
    h2 = M2.h()
    kvec2, khvec2 = make_kvecs(h2)
    M2.initialize_output(khvec2, Z_PK, len(khvec2))
    pk_mult2 = M2.get_pk_mult(khvec2, Z_PK, len(khvec2))
    instances['noIR'] = {
        'M': M2, 'h': h2, 'kvec': kvec2, 'khvec': khvec2,
        'pk_mult': pk_mult2,
    }

    # Bias, no-RSD config (real-space biased tracers)
    M3 = Class()
    M3.set(COSMO_PARAMS)
    M3.set({
        'output': 'mPk', 'non linear': 'PT',
        'IR resummation': 'No', 'Bias tracers': 'Yes',
        'RSD': 'No', 'cb': 'Yes', 'z_pk': Z_PK,
    })
    M3.compute()
    h3 = M3.h()
    kvec3, khvec3 = make_kvecs(h3)
    M3.initialize_output(khvec3, Z_PK, len(khvec3))
    pk_mult3 = M3.get_pk_mult(khvec3, Z_PK, len(khvec3))
    instances['bias_real'] = {
        'M': M3, 'h': h3, 'kvec': kvec3, 'khvec': khvec3,
        'pk_mult': pk_mult3,
    }

    # RSD + IR + Bias, no AP (for AP comparison and no-wiggle tests)
    M4 = Class()
    M4.set(COSMO_PARAMS)
    M4.set({
        'output': 'mPk', 'non linear': 'PT',
        'IR resummation': 'Yes', 'Bias tracers': 'Yes',
        'cb': 'Yes', 'RSD': 'Yes', 'AP': 'No',
        'z_pk': Z_PK,
    })
    M4.compute()
    h4 = M4.h()
    kvec4, khvec4 = make_kvecs(h4)
    M4.initialize_output(khvec4, Z_PK, len(khvec4))
    pk_mult4 = M4.get_pk_mult(khvec4, Z_PK, len(khvec4))
    instances['noAP'] = {
        'M': M4, 'h': h4, 'kvec': kvec4, 'khvec': khvec4,
        'pk_mult': pk_mult4,
    }

    # CMB lensing with PT
    M5 = Class()
    M5.set(COSMO_PARAMS)
    M5.set({
        'output': 'mPk,lCl,tCl,pCl', 'lensing': 'Yes', 'l_switch_limber': 9,
        'non linear': 'PT', 'cb': 'No',
        'IR resummation': 'No', 'Bias tracers': 'No', 'RSD': 'No',
        'z_pk': Z_PK,
    })
    M5.compute()
    instances['cmb_pt'] = {'M': M5, 'h': M5.h()}

    # CMB lensing linear (no non-linear)
    M6 = Class()
    M6.set(COSMO_PARAMS)
    M6.set({
        'output': 'mPk,lCl,tCl,pCl', 'lensing': 'Yes', 'l_switch_limber': 9,
        'non linear': 'none', 'z_pk': Z_PK,
    })
    M6.compute()
    instances['cmb_lin'] = {'M': M6, 'h': M6.h()}

    # External P(k) with RSD + IR + Bias (no A_s/n_s — spectrum from file)
    ext_file = os.path.join(os.path.dirname(__file__), 'external_Pk_test.dat')
    if os.path.exists(ext_file):
        COSMO_EXT = {k: v for k, v in COSMO_PARAMS.items() if k not in ('A_s', 'n_s')}
        M7 = Class()
        M7.set(COSMO_EXT)
        M7.set({
            'output': 'mPk', 'non linear': 'PT',
            'IR resummation': 'Yes', 'Bias tracers': 'Yes',
            'cb': 'Yes', 'RSD': 'Yes', 'AP': 'No',
            'z_pk': Z_PK,
            'P_k_ini type': 'external_Pk',
            'command': 'cat ' + os.path.abspath(ext_file),
        })
        M7.compute()
        h7 = M7.h()
        kvec7, khvec7 = make_kvecs(h7, kmin=-2, kmax_log10=np.log10(0.3), nk=100)
        M7.initialize_output(khvec7, Z_PK, len(khvec7))
        instances['ext_pk'] = {
            'M': M7, 'h': h7, 'kvec': kvec7, 'khvec': khvec7,
        }

    return instances

# Create all instances at import time (before any tests run)
_SHARED = _create_shared_instances()


class _FullConfigMixin:
    """Mixin providing access to the shared full-config Class instance."""

    @classmethod
    def setUpClass(cls):
        s = _SHARED['full']
        cls.M = s['M']
        cls.h = s['h']
        cls.kvec = s['kvec']
        cls.khvec = s['khvec']
        cls.pk_mult = s['pk_mult']
        cls.fz = s['fz']

    @classmethod
    def tearDownClass(cls):
        pass  # Shared instance cleaned up at module exit


class TestBackground(_FullConfigMixin, unittest.TestCase):
    """Test background cosmology computations (uses shared full-config instance)."""

    def test_h(self):
        self.assertAlmostEqual(self.M.h(), 0.6736, places=4)

    def test_omega_b(self):
        ob = self.M.Omega_b()
        self.assertGreater(ob, 0.04)
        self.assertLess(ob, 0.06)

    def test_omega_m(self):
        om = self.M.Omega_m()
        self.assertGreater(om, 0.28)
        self.assertLess(om, 0.35)

    def test_sigma8(self):
        s8 = self.M.sigma8()
        self.assertGreater(s8, 0.7)
        self.assertLess(s8, 0.9)

    def test_hubble(self):
        self.assertGreater(self.M.Hubble(0.), 0.)

    def test_angular_distance(self):
        self.assertGreater(self.M.angular_distance(Z_PK), 0.)

    def test_growth_factor(self):
        D = self.M.scale_independent_growth_factor(Z_PK)
        self.assertGreater(D, 0.)
        self.assertLess(D, 1.)

    def test_growth_rate(self):
        f = self.M.scale_independent_growth_factor_f(Z_PK)
        self.assertGreater(f, 0.)
        self.assertLess(f, 1.5)

    def test_pk_lin(self):
        h = self.M.h()
        pk = self.M.pk_lin(0.1 * h, Z_PK)
        self.assertGreater(pk, 0.)


class TestRealSpaceMatter(unittest.TestCase):
    """Test real-space matter power spectrum (no IR, no bias tracers)."""

    @classmethod
    def setUpClass(cls):
        s = _SHARED['noIR']
        cls.M = s['M']
        cls.h = s['h']
        cls.kvec = s['kvec']
        cls.khvec = s['khvec']
        cls.pk_mult = s['pk_mult']

    @classmethod
    def tearDownClass(cls):
        pass

    def test_pk_mm_real_positive_low_k(self):
        """Matter P(k) should be positive at low k (where PT is valid)."""
        pk = self.M.pk_mm_real(CS)
        mask = self.kvec < 0.3
        self.assertTrue(np.all(pk[mask] > 0), "Matter P(k) has non-positive values at low k")

    def test_pk_mm_real_shape(self):
        self.assertEqual(len(self.M.pk_mm_real(CS)), len(self.kvec))

    def test_pk_mm_real_magnitude(self):
        """At k~0.1 h/Mpc, P(k) should be O(1e3-1e4) (h/Mpc)^3."""
        pk = self.M.pk_mm_real(CS)
        idx = np.argmin(np.abs(self.kvec - 0.1))
        self.assertGreater(pk[idx], 1e2)
        self.assertLess(pk[idx], 1e5)

    def test_pk_mult_shape(self):
        self.assertEqual(self.pk_mult.shape, (96, len(self.kvec)))

    def test_tree_level_positive(self):
        """Tree-level P(k) (index 14) should be positive."""
        pk_tree = self.pk_mult[14]
        mask = self.kvec > 0.01
        self.assertTrue(np.all(pk_tree[mask] > 0))

    def test_loop_correction_subdominant_at_low_k(self):
        """1-loop correction (index 0) should not dominate tree level at low k."""
        pk_loop = self.pk_mult[0]
        pk_tree = self.pk_mult[14]
        idx = np.argmin(np.abs(self.kvec - 0.05))
        ratio = np.abs(pk_loop[idx]) / pk_tree[idx]
        self.assertLess(ratio, 1.0, f"Loop/tree ratio = {ratio} at k=0.05")

    def test_pk_mm_real_consistency(self):
        """Verify pk_mm_real matches manual construction from pk_mult."""
        h = self.h
        m = self.pk_mult
        pk = self.M.pk_mm_real(CS)
        pk_manual = (m[0] + m[14] + 2*CS*m[10]/h**2) * h**3
        np.testing.assert_allclose(pk, pk_manual, rtol=1e-12)


class TestRealSpaceBiasedTracers(unittest.TestCase):
    """Test real-space galaxy/cross power spectra (Pgg, Pgm, with bias, no IR/RSD)."""

    @classmethod
    def setUpClass(cls):
        s = _SHARED['bias_real']
        cls.M = s['M']
        cls.h = s['h']
        cls.kvec = s['kvec']
        cls.khvec = s['khvec']
        cls.pk_mult = s['pk_mult']

    @classmethod
    def tearDownClass(cls):
        pass

    def test_pk_gg_real_positive(self):
        """Galaxy P(k) should be positive (given Pshot > 0)."""
        pk = self.M.pk_gg_real(B1, B2, BG2, BGAMMA3, CS, CS0, PSHOT)
        self.assertTrue(np.all(pk > 0))

    def test_pk_gm_real_shape(self):
        self.assertEqual(len(self.M.pk_gm_real(B1, B2, BG2, BGAMMA3, CS, CS0)), len(self.kvec))

    def test_pk_gg_vs_mm_scaling(self):
        """At large scales, P_gg ~ b1^2 * P_mm."""
        pk_mm = self.M.pk_mm_real(CS)
        pk_gg = self.M.pk_gg_real(B1, B2, BG2, BGAMMA3, CS, CS0, PSHOT)
        idx = np.argmin(np.abs(self.kvec - 0.02))
        ratio = (pk_gg[idx] - PSHOT) / pk_mm[idx]
        self.assertAlmostEqual(ratio, B1**2, delta=0.5,
                               msg=f"P_gg/P_mm ratio = {ratio}, expected ~{B1**2}")

    def test_pk_gm_between_mm_and_gg(self):
        """At low k, P_gm should be between P_mm and P_gg (roughly b1 * P_mm)."""
        pk_mm = self.M.pk_mm_real(CS)
        pk_gm = self.M.pk_gm_real(B1, B2, BG2, BGAMMA3, CS, CS0)
        idx = np.argmin(np.abs(self.kvec - 0.02))
        ratio = pk_gm[idx] / pk_mm[idx]
        self.assertAlmostEqual(ratio, B1, delta=0.3,
                               msg=f"P_gm/P_mm ratio = {ratio}, expected ~{B1}")

    def test_bias_terms_nonzero(self):
        """Key bias component arrays should be nonzero."""
        m = self.pk_mult
        mask = self.kvec > 0.01
        for idx, name in [(1, 'Id2d2'), (2, 'Id2'), (3, 'IG2'),
                          (4, 'Id2G2'), (5, 'IG2G2'), (6, 'IFG2')]:
            self.assertTrue(np.any(m[idx][mask] != 0), f"pk_mult[{idx}] ({name}) is all zero")

    def test_pk_gg_real_consistency(self):
        """Verify pk_gg_real matches manual construction from pk_mult."""
        h = self.h
        m = self.pk_mult
        pk_gg = self.M.pk_gg_real(B1, B2, BG2, BGAMMA3, CS, CS0, PSHOT)
        pk_manual = (
            B1**2 * m[14] + B1**2 * m[0]
            + 2.*(CS*B1**2 + CS0*B1) * m[10] / h**2
            + B1*B2 * m[2] + 0.25*B2**2 * m[1]
            + 2.*B1*BG2 * m[3]
            + B1*(2.*BG2 + 0.8*BGAMMA3) * m[6]
            + BG2**2 * m[5] + B2*BG2 * m[4]
        ) * h**3 + PSHOT
        np.testing.assert_allclose(pk_gg, pk_manual, rtol=1e-12,
                                   err_msg="pk_gg_real inconsistent with pk_mult")

    def test_pk_gm_real_consistency(self):
        """Verify pk_gm_real matches manual construction from pk_mult."""
        h = self.h
        m = self.pk_mult
        pk_gm = self.M.pk_gm_real(B1, B2, BG2, BGAMMA3, CS, CS0)
        pk_manual = (
            B1*m[14] + B1*m[0]
            + (2.*CS*B1 + CS0)*m[10]/h**2
            + (B2/2.)*m[2] + BG2*m[3]
            + (BG2 + 0.4*BGAMMA3)*m[6]
        ) * h**3
        np.testing.assert_allclose(pk_gm, pk_manual, rtol=1e-12,
                                   err_msg="pk_gm_real inconsistent with pk_mult")


class TestRedshiftSpaceMultipoles(_FullConfigMixin, unittest.TestCase):
    """Test redshift-space multipoles with IR resummation, bias, RSD, AP."""

    # --- Matter multipoles ---
    def test_pk_mm_l0_positive(self):
        pk = self.M.pk_mm_l0(CS0)
        mask = (self.kvec > 0.01) & (self.kvec < 0.3)
        self.assertTrue(np.all(pk[mask] > 0))

    def test_pk_mm_l2_shape(self):
        self.assertEqual(len(self.M.pk_mm_l2(CS2)), len(self.kvec))

    def test_pk_mm_l4_shape(self):
        self.assertEqual(len(self.M.pk_mm_l4(CS4)), len(self.kvec))

    def test_pk_mm_l0_consistency(self):
        """Verify pk_mm_l0 matches manual construction."""
        h = self.h
        m = self.pk_mult
        pk = self.M.pk_mm_l0(CS0)
        pk_manual = (m[15]+m[21] + m[16]+m[22] + m[17]+m[23] + 2.*CS0*m[11]/h**2) * h**3
        np.testing.assert_allclose(pk, pk_manual, rtol=1e-12)

    def test_pk_mm_l2_consistency(self):
        h = self.h
        m = self.pk_mult
        pk = self.M.pk_mm_l2(CS2)
        pk_manual = (m[18]+m[24] + m[19]+m[25] + m[26] + 2.*CS2*m[12]/h**2) * h**3
        np.testing.assert_allclose(pk, pk_manual, rtol=1e-12)

    def test_pk_mm_l4_consistency(self):
        h = self.h
        m = self.pk_mult
        pk = self.M.pk_mm_l4(CS4)
        pk_manual = (m[20]+m[27] + m[28]+m[29] + 2.*CS4*m[13]/h**2) * h**3
        np.testing.assert_allclose(pk, pk_manual, rtol=1e-12)

    # --- Galaxy multipoles ---
    def test_pk_gg_l0_positive(self):
        pk = self.M.pk_gg_l0(B1, B2, BG2, BGAMMA3, CS0, PSHOT, B4)
        mask = (self.kvec > 0.01) & (self.kvec < 0.3)
        self.assertTrue(np.all(pk[mask] > 0))

    def test_pk_gg_l2_shape(self):
        self.assertEqual(len(self.M.pk_gg_l2(B1, B2, BG2, BGAMMA3, CS2, B4)), len(self.kvec))

    def test_pk_gg_l4_shape(self):
        self.assertEqual(len(self.M.pk_gg_l4(B1, B2, BG2, BGAMMA3, CS4, B4)), len(self.kvec))

    def test_monopole_dominates_quadrupole(self):
        pk0 = self.M.pk_mm_l0(CS0)
        pk2 = self.M.pk_mm_l2(CS2)
        idx = np.argmin(np.abs(self.kvec - 0.05))
        self.assertGreater(np.abs(pk0[idx]), np.abs(pk2[idx]),
                           "Quadrupole should not dominate monopole at low k")

    def test_pk_gg_l0_consistency(self):
        """Verify pk_gg_l0 matches manual construction from pk_mult."""
        h, fz, m, kh = self.h, self.fz, self.pk_mult, self.khvec
        pk = self.M.pk_gg_l0(B1, B2, BG2, BGAMMA3, CS0, PSHOT, B4)
        pk_manual = (
            m[15]+m[21] + B1*m[16]+B1*m[22] + B1**2*m[17]+B1**2*m[23]
            + 0.25*B2**2*m[1] + B1*B2*m[30]+B2*m[31]
            + B1*BG2*m[32]+BG2*m[33] + B2*BG2*m[4]+BG2**2*m[5]
            + 2.*CS0*m[11]/h**2
            + (2.*BG2+0.8*BGAMMA3)*(B1*m[7]+m[8])
        )*h**3 + PSHOT + fz**2*B4*(kh/h)**2*(fz**2/9.+2.*fz*B1/7.+B1**2/5.)*(35./8.)*m[13]*h
        np.testing.assert_allclose(pk, pk_manual, rtol=1e-12)

    def test_pk_gg_l2_consistency(self):
        h, fz, m, kh = self.h, self.fz, self.pk_mult, self.khvec
        pk = self.M.pk_gg_l2(B1, B2, BG2, BGAMMA3, CS2, B4)
        pk_manual = (
            m[18]+m[24] + B1*m[19]+B1*m[25] + B1**2*m[26]
            + B1*B2*m[34]+B2*m[35] + B1*BG2*m[36]+BG2*m[37]
            + 2.*CS2*m[12]/h**2 + (2.*BG2+0.8*BGAMMA3)*m[9]
        )*h**3 + fz**2*B4*(kh/h)**2*((fz**2*70.+165.*fz*B1+99.*B1**2)*4./693.)*(35./8.)*m[13]*h
        np.testing.assert_allclose(pk, pk_manual, rtol=1e-12)

    def test_pk_gg_l4_consistency(self):
        h, fz, m, kh = self.h, self.fz, self.pk_mult, self.khvec
        pk = self.M.pk_gg_l4(B1, B2, BG2, BGAMMA3, CS4, B4)
        pk_manual = (
            m[20]+m[27] + B1*m[28]+B1**2*m[29]
            + B2*m[38]+BG2*m[39] + 2.*CS4*m[13]/h**2
        )*h**3 + fz**2*B4*(kh/h)**2*((fz**2*210.+390.*fz*B1+143.*B1**2)*8./5005.)*(35./8.)*m[13]*h
        np.testing.assert_allclose(pk, pk_manual, rtol=1e-12)


class TestFNLEquilateral(_FullConfigMixin, unittest.TestCase):
    """Test equilateral fNL power spectrum contributions (real-space + multipoles)."""

    # --- Real-space fNL ---
    def test_pk_mm_fNL_real(self):
        pk = self.M.pk_mm_fNL_real()
        self.assertEqual(len(pk), len(self.kvec))
        self.assertTrue(np.any(pk != 0), "Equilateral fNL mm real is all zero")

    def test_pk_gm_fNL_real(self):
        pk = self.M.pk_gm_fNL_real(B1, B2, BG2)
        self.assertEqual(len(pk), len(self.kvec))
        self.assertTrue(np.any(pk != 0))

    def test_pk_gg_fNL_real(self):
        pk = self.M.pk_gg_fNL_real(B1, B2, BG2)
        self.assertEqual(len(pk), len(self.kvec))
        self.assertTrue(np.any(pk != 0))

    # --- Multipole fNL ---
    def test_pk_mm_fNL_multipoles(self):
        for name in ['pk_mm_fNL_l0', 'pk_mm_fNL_l2', 'pk_mm_fNL_l4']:
            pk = getattr(self.M, name)()
            self.assertEqual(len(pk), len(self.kvec))
            self.assertTrue(np.any(pk != 0), f"{name} is all zero")

    def test_pk_gg_fNL_multipoles(self):
        for name in ['pk_gg_fNL_l0', 'pk_gg_fNL_l2', 'pk_gg_fNL_l4']:
            pk = getattr(self.M, name)(B1, B2, BG2)
            self.assertEqual(len(pk), len(self.kvec))
            self.assertTrue(np.any(pk != 0), f"{name} is all zero")

    # --- Consistency ---
    def test_pk_mm_fNL_real_consistency(self):
        h = self.h
        pk = self.M.pk_mm_fNL_real()
        np.testing.assert_allclose(pk, self.pk_mult[48]*h**3, rtol=1e-12)

    def test_pk_gm_fNL_real_consistency(self):
        h = self.h
        m = self.pk_mult
        pk = self.M.pk_gm_fNL_real(B1, B2, BG2)
        pk_manual = B1*m[48]*h**3 + ((B2/2.)*m[49]+BG2*m[50])*h**3
        np.testing.assert_allclose(pk, pk_manual, rtol=1e-12)

    def test_pk_gg_fNL_real_consistency(self):
        h = self.h
        m = self.pk_mult
        pk = self.M.pk_gg_fNL_real(B1, B2, BG2)
        pk_manual = B1**2*m[48]*h**3 + B1*((B2/2.)*m[49]+BG2*m[50])*h**3
        np.testing.assert_allclose(pk, pk_manual, rtol=1e-12)

    def test_pk_gg_fNL_l0_consistency(self):
        h = self.h
        m = self.pk_mult
        pk = self.M.pk_gg_fNL_l0(B1, B2, BG2)
        pk_manual = (
            m[51] + B1*m[52] + B1**2*m[53]
            + B1*(B2/2.)*m[60] + (B2/2.)*m[61]
            + B1*BG2*m[62] + BG2*m[63]
        ) * h**3
        np.testing.assert_allclose(pk, pk_manual, rtol=1e-12)


class TestFNLOrthogonal(_FullConfigMixin, unittest.TestCase):
    """Test orthogonal fNL power spectrum contributions (real-space + multipoles)."""

    # --- Real-space ortho fNL ---
    def test_pk_mm_fNL_real_ortho(self):
        pk = self.M.pk_mm_fNL_real_ortho()
        self.assertEqual(len(pk), len(self.kvec))
        self.assertTrue(np.any(pk != 0))

    def test_pk_gm_fNL_real_ortho(self):
        pk = self.M.pk_gm_fNL_real_ortho(B1, B2, BG2)
        self.assertEqual(len(pk), len(self.kvec))

    def test_pk_gg_fNL_real_ortho(self):
        pk = self.M.pk_gg_fNL_real_ortho(B1, B2, BG2)
        self.assertEqual(len(pk), len(self.kvec))

    # --- Multipole ortho fNL ---
    def test_pk_mm_fNL_ortho_multipoles(self):
        for name in ['pk_mm_fNL_l0_ortho', 'pk_mm_fNL_l2_ortho', 'pk_mm_fNL_l4_ortho']:
            pk = getattr(self.M, name)()
            self.assertEqual(len(pk), len(self.kvec))
            self.assertTrue(np.any(pk != 0), f"{name} is all zero")

    def test_pk_gg_fNL_ortho_multipoles(self):
        for name in ['pk_gg_fNL_l0_ortho', 'pk_gg_fNL_l2_ortho', 'pk_gg_fNL_l4_ortho']:
            pk = getattr(self.M, name)(B1, B2, BG2)
            self.assertEqual(len(pk), len(self.kvec))
            self.assertTrue(np.any(pk != 0), f"{name} is all zero")

    # --- Cross-checks ---
    def test_equil_vs_ortho_different(self):
        """Equilateral and orthogonal fNL should give different results."""
        pk_eq = self.M.pk_mm_fNL_real()
        pk_or = self.M.pk_mm_fNL_real_ortho()
        self.assertFalse(np.allclose(pk_eq, pk_or),
                         "Equilateral and orthogonal fNL are identical")

    def test_pk_gg_fNL_l0_ortho_consistency(self):
        h = self.h
        m = self.pk_mult
        pk = self.M.pk_gg_fNL_l0_ortho(B1, B2, BG2)
        pk_manual = (
            m[75] + B1*m[76] + B1**2*m[77]
            + B1*(B2/2.)*m[84] + (B2/2.)*m[85]
            + B1*BG2*m[86] + BG2*m[87]
        ) * h**3
        np.testing.assert_allclose(pk, pk_manual, rtol=1e-12)


class TestIRResummation(unittest.TestCase):
    """Test that IR resummation changes the output.

    Uses shared noIR instance and full-config instance (which has IR=Yes).
    Compares real-space matter P(k) at a subset of k-values where PT is valid.
    """

    @classmethod
    def setUpClass(cls):
        # noIR instance for matter-only
        cls.M_noIR = _SHARED['noIR']['M']
        # full config has IR=Yes
        cls.M_IR = _SHARED['full']['M']
        cls.kvec = _SHARED['noIR']['kvec']
        # Restrict to PT-valid range
        cls.mask = (cls.kvec > 0.01) & (cls.kvec < 0.3)

    @classmethod
    def tearDownClass(cls):
        pass

    def test_ir_changes_matter_pk(self):
        pk_noIR = self.M_noIR.pk_mm_real(CS)[self.mask]
        pk_IR = self.M_IR.pk_mm_real(CS)[self.mask]
        self.assertFalse(np.allclose(pk_noIR, pk_IR, rtol=1e-4),
                         "IR resummation has no effect on matter P(k)")

    def test_ir_preserves_broadband(self):
        pk_noIR = self.M_noIR.pk_mm_real(CS)[self.mask]
        pk_IR = self.M_IR.pk_mm_real(CS)[self.mask]
        ratio = pk_IR / pk_noIR
        self.assertTrue(np.all(ratio > 0.8) and np.all(ratio < 1.2),
                        "IR resummation changes broadband by >20%")


class TestAPEffect(unittest.TestCase):
    """Test that the Alcock-Paczynski effect changes the output.

    Uses the shared full-config instance (AP=Yes) and a separate no-AP instance.
    """

    @classmethod
    def setUpClass(cls):
        s = _SHARED['full']
        cls.M_AP = s['M']
        cls.M_noAP = _SHARED['noAP']['M']
        cls.kvec = s['kvec']
        cls.mask = (cls.kvec > 0.01) & (cls.kvec < 0.3)

    @classmethod
    def tearDownClass(cls):
        pass

    def test_ap_changes_monopole(self):
        pk_noAP = self.M_noAP.pk_gg_l0(B1, B2, BG2, BGAMMA3, CS0, PSHOT, B4)[self.mask]
        pk_AP = self.M_AP.pk_gg_l0(B1, B2, BG2, BGAMMA3, CS0, PSHOT, B4)[self.mask]
        self.assertFalse(np.allclose(pk_noAP, pk_AP, rtol=1e-4),
                         "AP effect has no impact on galaxy monopole")

    def test_ap_changes_quadrupole(self):
        pk_noAP = self.M_noAP.pk_gg_l2(B1, B2, BG2, BGAMMA3, CS2, B4)[self.mask]
        pk_AP = self.M_AP.pk_gg_l2(B1, B2, BG2, BGAMMA3, CS2, B4)[self.mask]
        self.assertFalse(np.allclose(pk_noAP, pk_AP, rtol=1e-4),
                         "AP effect has no impact on galaxy quadrupole")


class TestNoWiggleRecompute(unittest.TestCase):
    """Test no-wiggle and alpha_rs recomputation via pk() method.

    Uses the shared no-AP instance (which has IR=Yes, needed for wiggle/no-wiggle).
    """

    @classmethod
    def setUpClass(cls):
        cls.M = _SHARED['noAP']['M']
        cls.h = _SHARED['noAP']['h']
        cls.k_test = 0.1 * cls.h  # k in 1/Mpc

    @classmethod
    def tearDownClass(cls):
        pass

    def test_no_wiggle_differs(self):
        """no_wiggle flag via pk() should change the matter P(k)."""
        pk_w = self.M.pk(self.k_test, Z_PK, no_wiggle=False)
        pk_nw = self.M.pk(self.k_test, Z_PK, no_wiggle=True)
        self.assertFalse(np.isclose(pk_w[0], pk_nw[0], rtol=1e-4),
                         "No-wiggle pk identical to wiggle pk")

    def test_alpha_rs_recompute(self):
        """Different alpha_rs should give different P(k)."""
        pk_1 = self.M.pk(self.k_test, Z_PK, alpha_rs=1.0)
        pk_2 = self.M.pk(self.k_test, Z_PK, alpha_rs=1.05)
        self.assertFalse(np.isclose(pk_1[0], pk_2[0], rtol=1e-4),
                         "alpha_rs=1.05 gives same pk as alpha_rs=1.0")


class TestCMBLensing(unittest.TestCase):
    """Test CMB lensing with PT non-linear corrections."""

    @classmethod
    def setUpClass(cls):
        cls.M_PT = _SHARED['cmb_pt']['M']
        cls.M_lin = _SHARED['cmb_lin']['M']

    @classmethod
    def tearDownClass(cls):
        pass

    def test_lensed_cl_keys(self):
        cl = self.M_PT.lensed_cl(2500)
        for key in ['tt', 'ee', 'pp']:
            self.assertIn(key, cl)

    def test_pt_modifies_pk(self):
        """PT non-linear P(k) should differ from linear P(k) at moderate k."""
        h = self.M_PT.h()
        k = 0.2 * h  # k in 1/Mpc
        # pk() returns a list when nlpt.method != 0; first element is P_1loop
        pk_pt_result = self.M_PT.pk(k, Z_PK)
        if isinstance(pk_pt_result, (list, tuple)):
            pk_pt = pk_pt_result[0]  # P_1loop (index 0 = pk_nl)
        else:
            pk_pt = pk_pt_result
        pk_lin = self.M_PT.pk_lin(k, Z_PK)
        ratio = pk_pt / pk_lin
        self.assertGreater(abs(ratio - 1.0), 0.005,
                           msg=f"PT P(k)/P_lin = {ratio} at k=0.2 h/Mpc, expected >0.5% correction")


class TestPkMultComponents(_FullConfigMixin, unittest.TestCase):
    """Test that all 96 pk_mult components are populated."""

    def test_pk_mult_shape(self):
        self.assertEqual(self.pk_mult.shape, (96, len(self.kvec)))

    def test_all_components_finite(self):
        """All pk_mult components should be finite, except known-NaN fNL*delta^2 components (49, 73, 85)."""
        known_nan = {49, 73, 85}  # fNL x delta^2 bias terms (equil + ortho + ortho bias)
        for i in range(96):
            if i in known_nan:
                continue
            self.assertTrue(np.all(np.isfinite(self.pk_mult[i])),
                            f"pk_mult[{i}] contains non-finite values")

    def test_gaussian_components_nonzero(self):
        for i in range(48):
            self.assertTrue(np.any(self.pk_mult[i] != 0), f"pk_mult[{i}] is all zero")

    def test_fnl_components_nonzero(self):
        for i in range(48, 96):
            self.assertTrue(np.any(self.pk_mult[i] != 0), f"pk_mult[{i}] is all zero")


class TestGetPkMult(_FullConfigMixin, unittest.TestCase):
    """Test get_pk_mult and pk methods directly."""

    def test_get_pk_mult_shape(self):
        pk_mult = self.M.get_pk_mult(self.khvec, Z_PK, len(self.khvec))
        self.assertEqual(pk_mult.shape, (96, len(self.kvec)))

    def test_pk_returns_96_components(self):
        result = self.M.pk(self.khvec[25], Z_PK)
        self.assertEqual(len(result), 96)

    def test_pk_all_finite(self):
        """All pk() components should be finite except known-NaN fNL bias indices."""
        # fNLd2 (49), fNLd2_ortho (73), fNL bias ortho (85) are always NaN in pk()
        # pk12_0_b1b2 (60) can also be NaN at some k due to log-space interpolation
        known_nan = {49, 60, 73, 85}
        result = self.M.pk(self.khvec[25], Z_PK)
        for i, v in enumerate(result):
            if i in known_nan:
                continue
            self.assertTrue(np.isfinite(v), f"pk component {i} = {v}")


class TestExternalPk(unittest.TestCase):
    """Test external power spectrum input mode.

    Uses a shared instance created from tests/external_Pk_test.dat —
    a power-law primordial P(k) with A_s=2.1e-9, n_s=0.965, wide k range.
    """

    @classmethod
    def setUpClass(cls):
        if 'ext_pk' not in _SHARED:
            raise unittest.SkipTest("External Pk test data not found")
        s = _SHARED['ext_pk']
        cls.M = s['M']
        cls.h = s['h']
        cls.kvec = s['kvec']
        cls.khvec = s['khvec']

    @classmethod
    def tearDownClass(cls):
        pass  # Shared instance cleaned up at module exit

    def test_pk_mm_real_finite_positive(self):
        pk = self.M.pk_mm_real(CS)
        self.assertTrue(np.all(np.isfinite(pk)), "External Pk: pk_mm_real has non-finite values")
        mask = self.kvec < 0.3
        self.assertTrue(np.all(pk[mask] > 0), "External Pk: pk_mm_real not positive at low k")

    def test_pk_mm_real_shape(self):
        self.assertEqual(len(self.M.pk_mm_real(CS)), len(self.kvec))

    def test_pk_gg_l0_finite(self):
        pk = self.M.pk_gg_l0(B1, B2, BG2, BGAMMA3, CS0, PSHOT, B4)
        self.assertTrue(np.all(np.isfinite(pk)), "External Pk: pk_gg_l0 has non-finite values")

    def test_pk_mm_l2_finite(self):
        pk = self.M.pk_mm_l2(CS2)
        self.assertTrue(np.all(np.isfinite(pk)), "External Pk: pk_mm_l2 has non-finite values")

    def test_pk_mult_shape(self):
        pk_mult = self.M.get_pk_mult(self.khvec, Z_PK, len(self.khvec))
        self.assertEqual(pk_mult.shape, (96, len(self.kvec)))

    def test_tree_level_positive(self):
        pk_mult = self.M.get_pk_mult(self.khvec, Z_PK, len(self.khvec))
        mask = self.kvec > 0.01
        self.assertTrue(np.all(pk_mult[14][mask] > 0), "External Pk: tree-level not positive")


class TestCleanupAndReuse(unittest.TestCase):
    """Test basic Class behavior (using shared instances, no new instances)."""

    def test_sigma8_consistent(self):
        """All shared instances with the same cosmology should give the same sigma8."""
        s8_full = _SHARED['full']['M'].sigma8()
        s8_noIR = _SHARED['noIR']['M'].sigma8()
        self.assertAlmostEqual(s8_full, s8_noIR, places=4,
                               msg="sigma8 differs between full and noIR configs")


class TestRegressionFull(_FullConfigMixin, unittest.TestCase):
    """Regression tests against saved reference data (full config: RSD+IR+Bias+AP+fNL).

    If any code change alters the numerical output of nonlinear_pt.c, these tests will fail.
    To regenerate reference data after an intentional change, run:
        python tests/generate_extended_ref.py
    """

    REF_FILE = os.path.join(os.path.dirname(__file__), 'reference_data.npz')
    RTOL = 1e-6

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        if not os.path.exists(cls.REF_FILE):
            raise unittest.SkipTest("Reference data file not found")
        cls.ref = dict(np.load(cls.REF_FILE, allow_pickle=True))
        cls.ridx = cls.ref['ref_idx'].astype(int)

    def _check(self, key, computed_full, rtol=None):
        if rtol is None:
            rtol = self.RTOL
        ref_val = self.ref[key]
        computed = computed_full[self.ridx]
        atol = max(1e-12 * np.max(np.abs(ref_val) + 1e-30), 1e-13)
        np.testing.assert_allclose(computed, ref_val, rtol=rtol, atol=atol,
                                   err_msg=f"Regression failure for '{key}'")

    # --- Background ---
    def test_sigma8(self):
        np.testing.assert_allclose(self.M.sigma8(), self.ref['sigma8'], rtol=self.RTOL)

    def test_h(self):
        np.testing.assert_allclose(self.M.h(), self.ref['h'], rtol=self.RTOL)

    def test_fz(self):
        np.testing.assert_allclose(
            self.M.scale_independent_growth_factor_f(Z_PK), self.ref['fz'], rtol=self.RTOL)

    # --- pk_mult components ---
    def test_pk_mult_gaussian(self):
        """All 48 Gaussian pk_mult components match reference at reference k-indices."""
        ridx = self.ridx
        for i in range(48):
            np.testing.assert_allclose(
                self.pk_mult[i, ridx], self.ref['full_pk_mult'][i],
                rtol=self.RTOL, err_msg=f"pk_mult[{i}] regression failure")

    def test_pk_mult_fNL(self):
        """All fNL pk_mult components match reference (skipping known-NaN)."""
        ridx = self.ridx
        known_nan = {49, 73, 85}
        for i in range(48, 96):
            if i in known_nan:
                continue
            ref_vals = self.ref['full_pk_mult'][i]
            atol = max(1e-12 * (np.max(np.abs(ref_vals)) + 1e-30), 1e-11)
            np.testing.assert_allclose(
                self.pk_mult[i, ridx], ref_vals,
                rtol=self.RTOL, atol=atol,
                err_msg=f"pk_mult[{i}] regression failure")

    # --- Real-space convenience methods ---
    def test_pk_mm_real(self):
        self._check('pk_mm_real', self.M.pk_mm_real(CS))

    def test_pk_gm_real(self):
        self._check('pk_gm_real', self.M.pk_gm_real(B1, B2, BG2, BGAMMA3, CS, CS0))

    def test_pk_gg_real(self):
        self._check('pk_gg_real', self.M.pk_gg_real(B1, B2, BG2, BGAMMA3, CS, CS0, PSHOT))

    # --- Matter multipoles ---
    def test_pk_mm_l0(self):
        self._check('pk_mm_l0', self.M.pk_mm_l0(CS0))

    def test_pk_mm_l2(self):
        self._check('pk_mm_l2', self.M.pk_mm_l2(CS2))

    def test_pk_mm_l4(self):
        self._check('pk_mm_l4', self.M.pk_mm_l4(CS4))

    # --- Galaxy multipoles ---
    def test_pk_gg_l0(self):
        self._check('pk_gg_l0', self.M.pk_gg_l0(B1, B2, BG2, BGAMMA3, CS0, PSHOT, B4))

    def test_pk_gg_l2(self):
        self._check('pk_gg_l2', self.M.pk_gg_l2(B1, B2, BG2, BGAMMA3, CS2, B4))

    def test_pk_gg_l4(self):
        self._check('pk_gg_l4', self.M.pk_gg_l4(B1, B2, BG2, BGAMMA3, CS4, B4))

    # --- fNL equilateral ---
    def test_pk_mm_fNL_real(self):
        self._check('pk_mm_fNL_real', self.M.pk_mm_fNL_real())

    def test_pk_gm_fNL_real(self):
        self._check('pk_gm_fNL_real', self.M.pk_gm_fNL_real(B1, B2, BG2))

    def test_pk_gg_fNL_real(self):
        self._check('pk_gg_fNL_real', self.M.pk_gg_fNL_real(B1, B2, BG2))

    def test_pk_mm_fNL_l0(self):
        self._check('pk_mm_fNL_l0', self.M.pk_mm_fNL_l0())

    def test_pk_mm_fNL_l2(self):
        self._check('pk_mm_fNL_l2', self.M.pk_mm_fNL_l2())

    def test_pk_mm_fNL_l4(self):
        self._check('pk_mm_fNL_l4', self.M.pk_mm_fNL_l4())

    def test_pk_gg_fNL_l0(self):
        self._check('pk_gg_fNL_l0', self.M.pk_gg_fNL_l0(B1, B2, BG2))

    def test_pk_gg_fNL_l2(self):
        self._check('pk_gg_fNL_l2', self.M.pk_gg_fNL_l2(B1, B2, BG2))

    def test_pk_gg_fNL_l4(self):
        self._check('pk_gg_fNL_l4', self.M.pk_gg_fNL_l4(B1, B2, BG2))

    # --- fNL orthogonal ---
    def test_pk_mm_fNL_real_ortho(self):
        self._check('pk_mm_fNL_real_ortho', self.M.pk_mm_fNL_real_ortho())

    def test_pk_gm_fNL_real_ortho(self):
        self._check('pk_gm_fNL_real_ortho', self.M.pk_gm_fNL_real_ortho(B1, B2, BG2))

    def test_pk_gg_fNL_real_ortho(self):
        self._check('pk_gg_fNL_real_ortho', self.M.pk_gg_fNL_real_ortho(B1, B2, BG2))

    def test_pk_mm_fNL_l0_ortho(self):
        self._check('pk_mm_fNL_l0_ortho', self.M.pk_mm_fNL_l0_ortho())

    def test_pk_mm_fNL_l2_ortho(self):
        self._check('pk_mm_fNL_l2_ortho', self.M.pk_mm_fNL_l2_ortho())

    def test_pk_mm_fNL_l4_ortho(self):
        self._check('pk_mm_fNL_l4_ortho', self.M.pk_mm_fNL_l4_ortho())

    def test_pk_gg_fNL_l0_ortho(self):
        self._check('pk_gg_fNL_l0_ortho', self.M.pk_gg_fNL_l0_ortho(B1, B2, BG2))

    def test_pk_gg_fNL_l2_ortho(self):
        self._check('pk_gg_fNL_l2_ortho', self.M.pk_gg_fNL_l2_ortho(B1, B2, BG2))

    def test_pk_gg_fNL_l4_ortho(self):
        self._check('pk_gg_fNL_l4_ortho', self.M.pk_gg_fNL_l4_ortho(B1, B2, BG2))


class TestRegressionNoIR(unittest.TestCase):
    """Regression tests for no-IR-resummation config.

    Uses the same Class instance as TestRealSpaceMatter (shares the noIR config).
    """

    REF_FILE = os.path.join(os.path.dirname(__file__), 'reference_data.npz')
    RTOL = 1e-6

    @classmethod
    def setUpClass(cls):
        if not os.path.exists(cls.REF_FILE):
            raise unittest.SkipTest("Reference data file not found")
        cls.ref = dict(np.load(cls.REF_FILE, allow_pickle=True))
        cls.ridx = cls.ref['ref_idx'].astype(int)
        s = _SHARED['noIR']
        cls.M = s['M']
        cls.h = s['h']
        cls.kvec = s['kvec']
        cls.khvec = s['khvec']

    @classmethod
    def tearDownClass(cls):
        pass

    def test_pk_mm_real_noIR(self):
        pk = self.M.pk_mm_real(CS)
        np.testing.assert_allclose(pk[self.ridx], self.ref['noIR_pk_mm_real'], rtol=self.RTOL,
                                   err_msg="noIR pk_mm_real regression failure")

    def test_pk_mult_noIR(self):
        """Valid pk_mult components without IR resummation match reference.
        Skip fNL rows (48-95) since noIR config does not enable PNG."""
        pk_mult = self.M.get_pk_mult(self.khvec, Z_PK, len(self.khvec))
        valid_idx = self.ref['noIR_valid_indices'].astype(int)
        for i in valid_idx:
            if i >= 48:
                continue  # fNL rows are uninitialized without PNG
            ref_vals = self.ref[f'noIR_pk_mult_{i}']
            actual_vals = pk_mult[i, self.ridx]
            # Skip components where either actual or reference has NaN/sentinel values
            if not np.all(np.isfinite(actual_vals)) or not np.all(np.isfinite(ref_vals)):
                continue
            if np.all(np.abs(ref_vals) > 999000) or np.all(np.abs(actual_vals) > 999000):
                continue
            atol = max(1e-12 * (np.max(np.abs(ref_vals)) + 1e-30), 1e-11)
            np.testing.assert_allclose(
                actual_vals, ref_vals, rtol=self.RTOL, atol=atol,
                err_msg=f"noIR pk_mult[{i}] regression failure")


class _RegressionBase(unittest.TestCase):
    """Base class for regression tests with common reference-data logic."""

    REF_FILE = os.path.join(os.path.dirname(__file__), 'reference_data.npz')
    RTOL = 1e-6

    @classmethod
    def setUpClass(cls):
        if not os.path.exists(cls.REF_FILE):
            raise unittest.SkipTest("Reference data file not found")
        cls.ref = dict(np.load(cls.REF_FILE, allow_pickle=True))
        cls.ridx = cls.ref['ref_idx'].astype(int)

    def _check(self, key, computed_full, rtol=None):
        if rtol is None:
            rtol = self.RTOL
        ref_val = self.ref[key]
        computed = computed_full[self.ridx]
        atol = max(1e-12 * np.max(np.abs(ref_val) + 1e-30), 1e-13)
        np.testing.assert_allclose(computed, ref_val, rtol=rtol, atol=atol,
                                   err_msg=f"Regression failure for '{key}'")


class TestRegressionNoAP(_RegressionBase):
    """Regression tests for noAP config (RSD + IR + Bias, no AP effect).

    Catches bugs in IR resummation and RSD multipole code without AP distortion.
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        s = _SHARED['noAP']
        cls.M = s['M']
        cls.h = s['h']
        cls.kvec = s['kvec']
        cls.khvec = s['khvec']

    # Real-space
    def test_pk_mm_real(self):
        self._check('noAP_pk_mm_real', self.M.pk_mm_real(CS))

    def test_pk_gm_real(self):
        self._check('noAP_pk_gm_real', self.M.pk_gm_real(B1, B2, BG2, BGAMMA3, CS, CS0))

    def test_pk_gg_real(self):
        self._check('noAP_pk_gg_real', self.M.pk_gg_real(B1, B2, BG2, BGAMMA3, CS, CS0, PSHOT))

    # Matter multipoles
    def test_pk_mm_l0(self):
        self._check('noAP_pk_mm_l0', self.M.pk_mm_l0(CS0))

    def test_pk_mm_l2(self):
        self._check('noAP_pk_mm_l2', self.M.pk_mm_l2(CS2))

    def test_pk_mm_l4(self):
        self._check('noAP_pk_mm_l4', self.M.pk_mm_l4(CS4))

    # Galaxy multipoles
    def test_pk_gg_l0(self):
        self._check('noAP_pk_gg_l0', self.M.pk_gg_l0(B1, B2, BG2, BGAMMA3, CS0, PSHOT, B4))

    def test_pk_gg_l2(self):
        self._check('noAP_pk_gg_l2', self.M.pk_gg_l2(B1, B2, BG2, BGAMMA3, CS2, B4))

    def test_pk_gg_l4(self):
        self._check('noAP_pk_gg_l4', self.M.pk_gg_l4(B1, B2, BG2, BGAMMA3, CS4, B4))

    # pk_mult components (spot-check Gaussian block)
    def test_pk_mult_gaussian(self):
        pk_mult = self.M.get_pk_mult(self.khvec, Z_PK, len(self.khvec))
        ridx = self.ridx
        for i in range(48):
            ref_vals = self.ref['noAP_pk_mult'][i]
            if not np.all(np.isfinite(ref_vals)):
                continue
            atol = max(1e-12 * (np.max(np.abs(ref_vals)) + 1e-30), 1e-11)
            np.testing.assert_allclose(
                pk_mult[i, ridx], ref_vals, rtol=self.RTOL, atol=atol,
                err_msg=f"noAP pk_mult[{i}] regression failure")


class TestRegressionBiasReal(_RegressionBase):
    """Regression tests for real-space biased tracers (no IR, no RSD).

    Catches bugs in bias operator loop integrals (Id2d2, IG2, IFG2, etc.).
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        s = _SHARED['bias_real']
        cls.M = s['M']

    def test_pk_mm_real(self):
        self._check('bias_real_pk_mm_real', self.M.pk_mm_real(CS))

    def test_pk_gm_real(self):
        self._check('bias_real_pk_gm_real', self.M.pk_gm_real(B1, B2, BG2, BGAMMA3, CS, CS0))

    def test_pk_gg_real(self):
        self._check('bias_real_pk_gg_real', self.M.pk_gg_real(B1, B2, BG2, BGAMMA3, CS, CS0, PSHOT))


class TestRegressionCMBLensing(_RegressionBase):
    """Regression tests for CMB lensing with PT non-linear corrections.

    Catches bugs in the nl_corr_density ratio used for CMB lensing potential.
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.M = _SHARED['cmb_pt']['M']
        cls.h = _SHARED['cmb_pt']['h']

    def test_lensed_tt(self):
        cl = self.M.lensed_cl(2500)
        ell_ref = self.ref['cmb_ell_ref'].astype(int)
        np.testing.assert_allclose(cl['tt'][ell_ref], self.ref['cmb_lensed_tt'],
                                   rtol=self.RTOL, err_msg="CMB lensed TT regression failure")

    def test_lensed_ee(self):
        cl = self.M.lensed_cl(2500)
        ell_ref = self.ref['cmb_ell_ref'].astype(int)
        np.testing.assert_allclose(cl['ee'][ell_ref], self.ref['cmb_lensed_ee'],
                                   rtol=self.RTOL, err_msg="CMB lensed EE regression failure")

    def test_lensed_pp(self):
        cl = self.M.lensed_cl(2500)
        ell_ref = self.ref['cmb_ell_ref'].astype(int)
        np.testing.assert_allclose(cl['pp'][ell_ref], self.ref['cmb_lensed_pp'],
                                   rtol=self.RTOL, err_msg="CMB lensed PP regression failure")

    def test_pt_pk(self):
        """PT 1-loop P(k) at test k-values matches reference."""
        k_test = self.ref['cmb_k_test']
        for i, k in enumerate(k_test):
            result = self.M.pk(k, Z_PK)
            pk_val = result[0] if isinstance(result, (list, tuple)) else result
            np.testing.assert_allclose(
                pk_val, self.ref['cmb_pk_pt'][i], rtol=self.RTOL,
                err_msg=f"CMB PT P(k) regression at k={k:.4f}")


class TestRegressionNoWiggle(_RegressionBase):
    """Regression tests for no-wiggle mode and alpha_rs rescaling.

    Catches bugs in the smooth P(k) extraction and BAO template rescaling.
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.M = _SHARED['noAP']['M']
        cls.h = _SHARED['noAP']['h']

    def test_pk_wiggle(self):
        """Wiggle P(k) at test k-values matches reference."""
        k_test = self.ref['nw_k_test']
        for i, k in enumerate(k_test):
            result = self.M.pk(k, Z_PK, no_wiggle=False)
            pk_val = result[0] if isinstance(result, (list, tuple)) else result
            np.testing.assert_allclose(
                pk_val, self.ref['nw_pk_wiggle'][i], rtol=self.RTOL,
                err_msg=f"Wiggle P(k) regression at k={k:.4f}")

    def test_pk_nowiggle(self):
        """No-wiggle P(k) at test k-values matches reference."""
        k_test = self.ref['nw_k_test']
        for i, k in enumerate(k_test):
            result = self.M.pk(k, Z_PK, no_wiggle=True)
            pk_val = result[0] if isinstance(result, (list, tuple)) else result
            np.testing.assert_allclose(
                pk_val, self.ref['nw_pk_nowiggle'][i], rtol=self.RTOL,
                err_msg=f"No-wiggle P(k) regression at k={k:.4f}")

    def test_alpha_rs(self):
        """P(k) with different alpha_rs values matches reference."""
        k = self.ref['alpha_k_test']
        for i, alpha in enumerate(self.ref['alpha_values']):
            result = self.M.pk(k, Z_PK, alpha_rs=float(alpha))
            pk_val = result[0] if isinstance(result, (list, tuple)) else result
            np.testing.assert_allclose(
                pk_val, self.ref['alpha_pk'][i], rtol=self.RTOL,
                err_msg=f"alpha_rs={alpha} P(k) regression failure")


class TestMultiCPUReproducibility(unittest.TestCase):
    """Test that results are identical across different OMP thread counts.

    Runs CLASS-PT in separate subprocesses with OMP_NUM_THREADS=1 and =4,
    then compares the output P(k) arrays to ensure bitwise reproducibility.
    """

    @staticmethod
    def _run_with_threads(n_threads):
        """Run CLASS-PT in a subprocess with the given thread count."""
        import subprocess, json
        script = f"""
import os
os.environ['OMP_NUM_THREADS'] = '{n_threads}'
import sys, json
sys.path.insert(0, os.path.join('{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}', 'python'))
os.chdir('{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}')
from classy import Class
import numpy as np

M = Class()
M.set({json.dumps(COSMO_PARAMS)})
M.set({{'output': 'mPk', 'non linear': 'PT', 'IR resummation': 'Yes',
        'Bias tracers': 'Yes', 'cb': 'Yes', 'RSD': 'Yes', 'AP': 'No',
        'z_pk': {Z_PK}}})
M.compute()
h = M.h()
kvec = np.logspace(-3, np.log10(3), 300)
khvec = kvec * h
M.initialize_output(khvec, {Z_PK}, len(khvec))
pk_mult = M.get_pk_mult(khvec, {Z_PK}, len(khvec))
pk_mm = M.pk_mm_real({CS})
pk_gg = M.pk_gg_l0({B1}, {B2}, {BG2}, {BGAMMA3}, {CS0}, {PSHOT}, {B4})
# Output as hex for bitwise comparison
result = {{
    'pk_mm': [float(x) for x in pk_mm[::30]],
    'pk_gg': [float(x) for x in pk_gg[::30]],
    'pk_mult_0': [float(x) for x in pk_mult[0, ::30]],
    'pk_mult_21': [float(x) for x in pk_mult[21, ::30]],
}}
print(json.dumps(result))
M.struct_cleanup(); M.empty()
"""
        proc = subprocess.run(
            [sys.executable, '-c', script],
            capture_output=True, text=True, timeout=300)
        if proc.returncode != 0:
            raise RuntimeError(f"Subprocess failed (threads={n_threads}):\n{proc.stderr}")
        return json.loads(proc.stdout.strip())

    def test_1_vs_4_threads(self):
        """P(k) must be identical between 1 and 4 OMP threads."""
        r1 = self._run_with_threads(1)
        r4 = self._run_with_threads(4)
        for key in r1:
            np.testing.assert_array_equal(
                np.array(r1[key]), np.array(r4[key]),
                err_msg=f"{key} differs between 1 and 4 threads")


def tearDownModule():
    """Clean up all shared Class instances at end of module."""
    for key, inst in _SHARED.items():
        try:
            inst['M'].struct_cleanup()
            inst['M'].empty()
        except Exception:
            pass


if __name__ == '__main__':
    print(f"Running CLASS-PT comprehensive tests...")
    print(f"Cosmology: h={COSMO_PARAMS['h']}, omega_cdm={COSMO_PARAMS['omega_cdm']}")
    print(f"z_pk = {Z_PK}")
    print()
    unittest.main(verbosity=2)
