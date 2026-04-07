/** @file nonlinear_pt_output.c
 *
 * Output/interpolation routines for the PT power spectra.
 *
 * These functions provide:
 * - nonlinear_pt_bias_at_z_i(): Copy PT spectra from pnlpt arrays at a given z index
 * - nonlinear_pt_pk_at_k_and_z(): Interpolate PT spectra to arbitrary k at a given z
 * - nonlinear_pt_pk_mult_at_kvec_and_z(): Batch interpolation of all 96 PT components
 */

#include "nonlinear_pt.h"

#define _NLPT_NUM_COMPONENTS_ 96

/**
 * Return pointer to the ln_pk_* array for component index c (0..95).
 * The returned pointer has layout [i_z * k_size + index_k].
 *
 * The ordering matches the parameter order of nonlinear_pt_bias_at_z_i()
 * and nonlinear_pt_pk_at_k_and_z(), i.e. component 0 = P_nl, 1 = Id2d2, etc.
 */
static double * _nlpt_component_array(struct nonlinear_pt *pnlpt, int c) {
    switch (c) {
    /* Gaussian real-space (0-14) */
    case  0: return pnlpt->ln_pk_nl;
    case  1: return pnlpt->ln_pk_Id2d2;
    case  2: return pnlpt->ln_pk_Id2;
    case  3: return pnlpt->ln_pk_IG2;
    case  4: return pnlpt->ln_pk_Id2G2;
    case  5: return pnlpt->ln_pk_IG2G2;
    case  6: return pnlpt->ln_pk_IFG2;
    case  7: return pnlpt->ln_pk_IFG2_0b1;
    case  8: return pnlpt->ln_pk_IFG2_0;
    case  9: return pnlpt->ln_pk_IFG2_2;
    case 10: return pnlpt->ln_pk_CTR;
    case 11: return pnlpt->ln_pk_CTR_0;
    case 12: return pnlpt->ln_pk_CTR_2;
    case 13: return pnlpt->ln_pk_CTR_4;
    case 14: return pnlpt->ln_pk_Tree;
    /* Tree-level RSD multipoles (15-20) */
    case 15: return pnlpt->ln_pk_Tree_0_vv;
    case 16: return pnlpt->ln_pk_Tree_0_vd;
    case 17: return pnlpt->ln_pk_Tree_0_dd;
    case 18: return pnlpt->ln_pk_Tree_2_vv;
    case 19: return pnlpt->ln_pk_Tree_2_vd;
    case 20: return pnlpt->ln_pk_Tree_4_vv;
    /* One-loop RSD multipoles (21-29) */
    case 21: return pnlpt->ln_pk_0_vv;
    case 22: return pnlpt->ln_pk_0_vd;
    case 23: return pnlpt->ln_pk_0_dd;
    case 24: return pnlpt->ln_pk_2_vv;
    case 25: return pnlpt->ln_pk_2_vd;
    case 26: return pnlpt->ln_pk_2_dd;
    case 27: return pnlpt->ln_pk_4_vv;
    case 28: return pnlpt->ln_pk_4_vd;
    case 29: return pnlpt->ln_pk_4_dd;
    /* Bias multipoles (30-47) */
    case 30: return pnlpt->ln_pk_0_b1b2;
    case 31: return pnlpt->ln_pk_0_b2;
    case 32: return pnlpt->ln_pk_0_b1bG2;
    case 33: return pnlpt->ln_pk_0_bG2;
    case 34: return pnlpt->ln_pk_2_b1b2;
    case 35: return pnlpt->ln_pk_2_b2;
    case 36: return pnlpt->ln_pk_2_b1bG2;
    case 37: return pnlpt->ln_pk_2_bG2;
    case 38: return pnlpt->ln_pk_4_b2;
    case 39: return pnlpt->ln_pk_4_bG2;
    case 40: return pnlpt->ln_pk_4_b1b2;
    case 41: return pnlpt->ln_pk_4_b1bG2;
    case 42: return pnlpt->ln_pk_Id2d2_2;
    case 43: return pnlpt->ln_pk_Id2G2_2;
    case 44: return pnlpt->ln_pk_IG2G2_2;
    case 45: return pnlpt->ln_pk_Id2d2_4;
    case 46: return pnlpt->ln_pk_Id2G2_4;
    case 47: return pnlpt->ln_pk_IG2G2_4;
    /* fNL equilateral (48-71) */
    case 48: return pnlpt->ln_pk_nl_fNL;
    case 49: return pnlpt->ln_pk_fNLd2;
    case 50: return pnlpt->ln_pk_fNLG2;
    case 51: return pnlpt->ln_pk_fNL_0_vv;
    case 52: return pnlpt->ln_pk_fNL_0_vd;
    case 53: return pnlpt->ln_pk_fNL_0_dd;
    case 54: return pnlpt->ln_pk_fNL_2_vv;
    case 55: return pnlpt->ln_pk_fNL_2_vd;
    case 56: return pnlpt->ln_pk_fNL_2_dd;
    case 57: return pnlpt->ln_pk_fNL_4_vv;
    case 58: return pnlpt->ln_pk_fNL_4_vd;
    case 59: return pnlpt->ln_pk_fNL_4_dd;
    case 60: return pnlpt->ln_pk12_0_b1b2;
    case 61: return pnlpt->ln_pk12_0_b2;
    case 62: return pnlpt->ln_pk12_0_b1bG2;
    case 63: return pnlpt->ln_pk12_0_bG2;
    case 64: return pnlpt->ln_pk12_2_b1b2;
    case 65: return pnlpt->ln_pk12_2_b2;
    case 66: return pnlpt->ln_pk12_2_b1bG2;
    case 67: return pnlpt->ln_pk12_2_bG2;
    case 68: return pnlpt->ln_pk12_4_b1b2;
    case 69: return pnlpt->ln_pk12_4_b2;
    case 70: return pnlpt->ln_pk12_4_b1bG2;
    case 71: return pnlpt->ln_pk12_4_bG2;
    /* fNL orthogonal (72-95) */
    case 72: return pnlpt->ln_pk_nl_fNL_ortho;
    case 73: return pnlpt->ln_pk_fNLd2_ortho;
    case 74: return pnlpt->ln_pk_fNLG2_ortho;
    case 75: return pnlpt->ln_pk_fNL_0_vv_ortho;
    case 76: return pnlpt->ln_pk_fNL_0_vd_ortho;
    case 77: return pnlpt->ln_pk_fNL_0_dd_ortho;
    case 78: return pnlpt->ln_pk_fNL_2_vv_ortho;
    case 79: return pnlpt->ln_pk_fNL_2_vd_ortho;
    case 80: return pnlpt->ln_pk_fNL_2_dd_ortho;
    case 81: return pnlpt->ln_pk_fNL_4_vv_ortho;
    case 82: return pnlpt->ln_pk_fNL_4_vd_ortho;
    case 83: return pnlpt->ln_pk_fNL_4_dd_ortho;
    case 84: return pnlpt->ln_pk12_0_b1b2_ortho;
    case 85: return pnlpt->ln_pk12_0_b2_ortho;
    case 86: return pnlpt->ln_pk12_0_b1bG2_ortho;
    case 87: return pnlpt->ln_pk12_0_bG2_ortho;
    case 88: return pnlpt->ln_pk12_2_b1b2_ortho;
    case 89: return pnlpt->ln_pk12_2_b2_ortho;
    case 90: return pnlpt->ln_pk12_2_b1bG2_ortho;
    case 91: return pnlpt->ln_pk12_2_bG2_ortho;
    case 92: return pnlpt->ln_pk12_4_b1b2_ortho;
    case 93: return pnlpt->ln_pk12_4_b2_ortho;
    case 94: return pnlpt->ln_pk12_4_b1bG2_ortho;
    case 95: return pnlpt->ln_pk12_4_bG2_ortho;
    default: return NULL;
    }
}

/**
 * Pack 96 output pointers into an array matching the component ordering
 * of _nlpt_component_array(). This macro is used in bias_at_z_i and
 * pk_at_k_and_z to avoid repeating the same 96-element initializer.
 *
 * Arguments: P = prefix of variable names (e.g. output_tot or pk_tot).
 * The caller must have local variables named P, P_Id2d2, P_Id2, etc.
 */
#define _NLPT_PACK_OUTPUTS_(P) {                                            \
    P, P##_Id2d2, P##_Id2, P##_IG2, P##_Id2G2, P##_IG2G2,                 \
    P##_IFG2, P##_IFG2_0b1, P##_IFG2_0, P##_IFG2_2,                      \
    P##_CTR, P##_CTR_0, P##_CTR_2, P##_CTR_4, P##_Tree,                   \
    P##_Tree_0_vv, P##_Tree_0_vd, P##_Tree_0_dd,                          \
    P##_Tree_2_vv, P##_Tree_2_vd, P##_Tree_4_vv,                          \
    P##_0_vv, P##_0_vd, P##_0_dd,                                         \
    P##_2_vv, P##_2_vd, P##_2_dd,                                         \
    P##_4_vv, P##_4_vd, P##_4_dd,                                         \
    P##_0_b1b2, P##_0_b2, P##_0_b1bG2, P##_0_bG2,                        \
    P##_2_b1b2, P##_2_b2, P##_2_b1bG2, P##_2_bG2,                        \
    P##_4_b2, P##_4_bG2, P##_4_b1b2, P##_4_b1bG2,                        \
    P##_2_b2b2, P##_2_b2bG2, P##_2_bG2bG2,                               \
    P##_4_b2b2, P##_4_b2bG2, P##_4_bG2bG2,                               \
    P##_pk_nl_fNL, P##_pk_fNLd2, P##_pk_fNLG2,                           \
    P##_pk_l_fNL_0_vv, P##_pk_l_fNL_0_vd, P##_pk_l_fNL_0_dd,            \
    P##_pk_l_fNL_2_vv, P##_pk_l_fNL_2_vd, P##_pk_l_fNL_2_dd,            \
    P##_pk_l_fNL_4_vv, P##_pk_l_fNL_4_vd, P##_pk_l_fNL_4_dd,            \
    P##_pk12_l_0_b1b2, P##_pk12_l_0_b2,                                   \
    P##_pk12_l_0_b1bG2, P##_pk12_l_0_bG2,                                 \
    P##_pk12_l_2_b1b2, P##_pk12_l_2_b2,                                   \
    P##_pk12_l_2_b1bG2, P##_pk12_l_2_bG2,                                 \
    P##_pk12_l_4_b1b2, P##_pk12_l_4_b2,                                   \
    P##_pk12_l_4_b1bG2, P##_pk12_l_4_bG2,                                 \
    P##_pk_nl_fNL_ortho, P##_pk_fNLd2_ortho, P##_pk_fNLG2_ortho,         \
    P##_pk_l_fNL_0_vv_ortho, P##_pk_l_fNL_0_vd_ortho,                    \
    P##_pk_l_fNL_0_dd_ortho,                                               \
    P##_pk_l_fNL_2_vv_ortho, P##_pk_l_fNL_2_vd_ortho,                    \
    P##_pk_l_fNL_2_dd_ortho,                                               \
    P##_pk_l_fNL_4_vv_ortho, P##_pk_l_fNL_4_vd_ortho,                    \
    P##_pk_l_fNL_4_dd_ortho,                                               \
    P##_pk12_l_0_b1b2_ortho, P##_pk12_l_0_b2_ortho,                       \
    P##_pk12_l_0_b1bG2_ortho, P##_pk12_l_0_bG2_ortho,                     \
    P##_pk12_l_2_b1b2_ortho, P##_pk12_l_2_b2_ortho,                       \
    P##_pk12_l_2_b1bG2_ortho, P##_pk12_l_2_bG2_ortho,                     \
    P##_pk12_l_4_b1b2_ortho, P##_pk12_l_4_b2_ortho,                       \
    P##_pk12_l_4_b1bG2_ortho, P##_pk12_l_4_bG2_ortho                      \
}


int nonlinear_pt_bias_at_z_i(
    struct nonlinear_pt *pnlpt,
    enum linear_or_logarithmic mode,
    int i_z,
    double * output_tot,
    double * output_tot_Id2d2,
    double * output_tot_Id2,
    double * output_tot_IG2,
    double * output_tot_Id2G2,
    double * output_tot_IG2G2,
    double * output_tot_IFG2,
    double * output_tot_IFG2_0b1,
    double * output_tot_IFG2_0,
    double * output_tot_IFG2_2,
    double * output_tot_CTR,
    double * output_tot_CTR_0,
    double * output_tot_CTR_2,
    double * output_tot_CTR_4,
    double * output_tot_Tree,
    double * output_tot_Tree_0_vv,
    double * output_tot_Tree_0_vd,
    double * output_tot_Tree_0_dd,
    double * output_tot_Tree_2_vv,
    double * output_tot_Tree_2_vd,
    double * output_tot_Tree_4_vv,
    double * output_tot_0_vv,
    double * output_tot_0_vd,
    double * output_tot_0_dd,
    double * output_tot_2_vv,
    double * output_tot_2_vd,
    double * output_tot_2_dd,
    double * output_tot_4_vv,
    double * output_tot_4_vd,
    double * output_tot_4_dd,
    double * output_tot_0_b1b2,
    double * output_tot_0_b2,
    double * output_tot_0_b1bG2,
    double * output_tot_0_bG2,
    double * output_tot_2_b1b2,
    double * output_tot_2_b2,
    double * output_tot_2_b1bG2,
    double * output_tot_2_bG2,
    double * output_tot_4_b2,
    double * output_tot_4_bG2,
    double * output_tot_4_b1b2,
    double * output_tot_4_b1bG2,
    double * output_tot_2_b2b2,
    double * output_tot_2_b2bG2,
    double * output_tot_2_bG2bG2,
    double * output_tot_4_b2b2,
    double * output_tot_4_b2bG2,
    double * output_tot_4_bG2bG2,
    double * output_tot_pk_nl_fNL,
    double * output_tot_pk_fNLd2,
    double * output_tot_pk_fNLG2,
    double * output_tot_pk_l_fNL_0_vv,
    double * output_tot_pk_l_fNL_0_vd,
    double * output_tot_pk_l_fNL_0_dd,
    double * output_tot_pk_l_fNL_2_vv,
    double * output_tot_pk_l_fNL_2_vd,
    double * output_tot_pk_l_fNL_2_dd,
    double * output_tot_pk_l_fNL_4_vv,
    double * output_tot_pk_l_fNL_4_vd,
    double * output_tot_pk_l_fNL_4_dd,
    double * output_tot_pk12_l_0_b1b2,
    double * output_tot_pk12_l_0_b2,
    double * output_tot_pk12_l_0_b1bG2,
    double * output_tot_pk12_l_0_bG2,
    double * output_tot_pk12_l_2_b1b2,
    double * output_tot_pk12_l_2_b2,
    double * output_tot_pk12_l_2_b1bG2,
    double * output_tot_pk12_l_2_bG2,
    double * output_tot_pk12_l_4_b1b2,
    double * output_tot_pk12_l_4_b2,
    double * output_tot_pk12_l_4_b1bG2,
    double * output_tot_pk12_l_4_bG2,
    double * output_tot_pk_nl_fNL_ortho,
    double * output_tot_pk_fNLd2_ortho,
    double * output_tot_pk_fNLG2_ortho,
    double * output_tot_pk_l_fNL_0_vv_ortho,
    double * output_tot_pk_l_fNL_0_vd_ortho,
    double * output_tot_pk_l_fNL_0_dd_ortho,
    double * output_tot_pk_l_fNL_2_vv_ortho,
    double * output_tot_pk_l_fNL_2_vd_ortho,
    double * output_tot_pk_l_fNL_2_dd_ortho,
    double * output_tot_pk_l_fNL_4_vv_ortho,
    double * output_tot_pk_l_fNL_4_vd_ortho,
    double * output_tot_pk_l_fNL_4_dd_ortho,
    double * output_tot_pk12_l_0_b1b2_ortho,
    double * output_tot_pk12_l_0_b2_ortho,
    double * output_tot_pk12_l_0_b1bG2_ortho,
    double * output_tot_pk12_l_0_bG2_ortho,
    double * output_tot_pk12_l_2_b1b2_ortho,
    double * output_tot_pk12_l_2_b2_ortho,
    double * output_tot_pk12_l_2_b1bG2_ortho,
    double * output_tot_pk12_l_2_bG2_ortho,
    double * output_tot_pk12_l_4_b1b2_ortho,
    double * output_tot_pk12_l_4_b2_ortho,
    double * output_tot_pk12_l_4_b1bG2_ortho,
    double * output_tot_pk12_l_4_bG2_ortho
) {

    int index_k, c;
    int nk = pnlpt->ln_k_size;
    int offset = i_z * pnlpt->k_size;
    double * outputs[_NLPT_NUM_COMPONENTS_] = _NLPT_PACK_OUTPUTS_(output_tot);

    for (c = 0; c < _NLPT_NUM_COMPONENTS_; c++) {
        double * src = _nlpt_component_array(pnlpt, c);
        for (index_k = 0; index_k < nk; index_k++)
            outputs[c][index_k] = src[offset + index_k];
    }

    if (mode == linear) {
        for (c = 0; c < _NLPT_NUM_COMPONENTS_; c++)
            for (index_k = 0; index_k < nk; index_k++)
                outputs[c][index_k] = exp(outputs[c][index_k]);
    }

    return _SUCCESS_;
}


int nonlinear_pt_pk_at_k_and_z(
    struct nonlinear_pt * pnlpt,
    double k,
    double z,
    double * pk_tot,
    double * pk_tot_Id2d2,
    double * pk_tot_Id2,
    double * pk_tot_IG2,
    double * pk_tot_Id2G2,
    double * pk_tot_IG2G2,
    double * pk_tot_IFG2,
    double * pk_tot_IFG2_0b1,
    double * pk_tot_IFG2_0,
    double * pk_tot_IFG2_2,
    double * pk_tot_CTR,
    double * pk_tot_CTR_0,
    double * pk_tot_CTR_2,
    double * pk_tot_CTR_4,
    double * pk_tot_Tree,
    double * pk_tot_Tree_0_vv,
    double * pk_tot_Tree_0_vd,
    double * pk_tot_Tree_0_dd,
    double * pk_tot_Tree_2_vv,
    double * pk_tot_Tree_2_vd,
    double * pk_tot_Tree_4_vv,
    double * pk_tot_0_vv,
    double * pk_tot_0_vd,
    double * pk_tot_0_dd,
    double * pk_tot_2_vv,
    double * pk_tot_2_vd,
    double * pk_tot_2_dd,
    double * pk_tot_4_vv,
    double * pk_tot_4_vd,
    double * pk_tot_4_dd,
    double * pk_tot_0_b1b2,
    double * pk_tot_0_b2,
    double * pk_tot_0_b1bG2,
    double * pk_tot_0_bG2,
    double * pk_tot_2_b1b2,
    double * pk_tot_2_b2,
    double * pk_tot_2_b1bG2,
    double * pk_tot_2_bG2,
    double * pk_tot_4_b2,
    double * pk_tot_4_bG2,
    double * pk_tot_4_b1b2,
    double * pk_tot_4_b1bG2,
    double * pk_tot_2_b2b2,
    double * pk_tot_2_b2bG2,
    double * pk_tot_2_bG2bG2,
    double * pk_tot_4_b2b2,
    double * pk_tot_4_b2bG2,
    double * pk_tot_4_bG2bG2,
    double * pk_tot_fNL,
    double * pk_tot_fNLd2,
    double * pk_tot_fNLG2,
    double * pk_tot_fNL_0_vv,
    double * pk_tot_fNL_0_vd,
    double * pk_tot_fNL_0_dd,
    double * pk_tot_fNL_2_vv,
    double * pk_tot_fNL_2_vd,
    double * pk_tot_fNL_2_dd,
    double * pk_tot_fNL_4_vv,
    double * pk_tot_fNL_4_vd,
    double * pk_tot_fNL_4_dd,
    double * pk_tot_fNL_0_b1b2,
    double * pk_tot_fNL_0_b2,
    double * pk_tot_fNL_0_b1bG2,
    double * pk_tot_fNL_0_bG2,
    double * pk_tot_fNL_2_b1b2,
    double * pk_tot_fNL_2_b2,
    double * pk_tot_fNL_2_b1bG2,
    double * pk_tot_fNL_2_bG2,
    double * pk_tot_fNL_4_b1b2,
    double * pk_tot_fNL_4_b2,
    double * pk_tot_fNL_4_b1bG2,
    double * pk_tot_fNL_4_bG2,
    double * pk_tot_fNL_ortho,
    double * pk_tot_fNLd2_ortho,
    double * pk_tot_fNLG2_ortho,
    double * pk_tot_fNL_0_vv_ortho,
    double * pk_tot_fNL_0_vd_ortho,
    double * pk_tot_fNL_0_dd_ortho,
    double * pk_tot_fNL_2_vv_ortho,
    double * pk_tot_fNL_2_vd_ortho,
    double * pk_tot_fNL_2_dd_ortho,
    double * pk_tot_fNL_4_vv_ortho,
    double * pk_tot_fNL_4_vd_ortho,
    double * pk_tot_fNL_4_dd_ortho,
    double * pk_tot_fNL_0_b1b2_ortho,
    double * pk_tot_fNL_0_b2_ortho,
    double * pk_tot_fNL_0_b1bG2_ortho,
    double * pk_tot_fNL_0_bG2_ortho,
    double * pk_tot_fNL_2_b1b2_ortho,
    double * pk_tot_fNL_2_b2_ortho,
    double * pk_tot_fNL_2_b1bG2_ortho,
    double * pk_tot_fNL_2_bG2_ortho,
    double * pk_tot_fNL_4_b1b2_ortho,
    double * pk_tot_fNL_4_b2_ortho,
    double * pk_tot_fNL_4_b1bG2_ortho,
    double * pk_tot_fNL_4_bG2_ortho
) {

    int c;
    double pk_buf[_NLPT_NUM_COMPONENTS_];

    class_test((k < exp(pnlpt->ln_k[0])) || (k > exp(pnlpt->ln_k[pnlpt->ln_k_size-1])),
               pnlpt->error_message,
               "k=%e out of bounds [%e:%e]",
               k, 0., exp(pnlpt->ln_k[pnlpt->ln_k_size-1]));

    class_call(nonlinear_pt_pk_mult_at_kvec_and_z(pnlpt, &k, 1, z, pk_buf),
               pnlpt->error_message, pnlpt->error_message);

    /* The macro cannot be used here because pk_at_k_and_z uses a different
     * naming convention for fNL parameters (pk_tot_fNL_* vs pk_tot_pk_*).
     * We must list explicitly. */
    double * outputs[_NLPT_NUM_COMPONENTS_] = {
        pk_tot, pk_tot_Id2d2, pk_tot_Id2, pk_tot_IG2,
        pk_tot_Id2G2, pk_tot_IG2G2, pk_tot_IFG2,
        pk_tot_IFG2_0b1, pk_tot_IFG2_0, pk_tot_IFG2_2,
        pk_tot_CTR, pk_tot_CTR_0, pk_tot_CTR_2, pk_tot_CTR_4,
        pk_tot_Tree,
        pk_tot_Tree_0_vv, pk_tot_Tree_0_vd, pk_tot_Tree_0_dd,
        pk_tot_Tree_2_vv, pk_tot_Tree_2_vd, pk_tot_Tree_4_vv,
        pk_tot_0_vv, pk_tot_0_vd, pk_tot_0_dd,
        pk_tot_2_vv, pk_tot_2_vd, pk_tot_2_dd,
        pk_tot_4_vv, pk_tot_4_vd, pk_tot_4_dd,
        pk_tot_0_b1b2, pk_tot_0_b2, pk_tot_0_b1bG2, pk_tot_0_bG2,
        pk_tot_2_b1b2, pk_tot_2_b2, pk_tot_2_b1bG2, pk_tot_2_bG2,
        pk_tot_4_b2, pk_tot_4_bG2, pk_tot_4_b1b2, pk_tot_4_b1bG2,
        pk_tot_2_b2b2, pk_tot_2_b2bG2, pk_tot_2_bG2bG2,
        pk_tot_4_b2b2, pk_tot_4_b2bG2, pk_tot_4_bG2bG2,
        pk_tot_fNL, pk_tot_fNLd2, pk_tot_fNLG2,
        pk_tot_fNL_0_vv, pk_tot_fNL_0_vd, pk_tot_fNL_0_dd,
        pk_tot_fNL_2_vv, pk_tot_fNL_2_vd, pk_tot_fNL_2_dd,
        pk_tot_fNL_4_vv, pk_tot_fNL_4_vd, pk_tot_fNL_4_dd,
        pk_tot_fNL_0_b1b2, pk_tot_fNL_0_b2,
        pk_tot_fNL_0_b1bG2, pk_tot_fNL_0_bG2,
        pk_tot_fNL_2_b1b2, pk_tot_fNL_2_b2,
        pk_tot_fNL_2_b1bG2, pk_tot_fNL_2_bG2,
        pk_tot_fNL_4_b1b2, pk_tot_fNL_4_b2,
        pk_tot_fNL_4_b1bG2, pk_tot_fNL_4_bG2,
        pk_tot_fNL_ortho, pk_tot_fNLd2_ortho, pk_tot_fNLG2_ortho,
        pk_tot_fNL_0_vv_ortho, pk_tot_fNL_0_vd_ortho,
        pk_tot_fNL_0_dd_ortho,
        pk_tot_fNL_2_vv_ortho, pk_tot_fNL_2_vd_ortho,
        pk_tot_fNL_2_dd_ortho,
        pk_tot_fNL_4_vv_ortho, pk_tot_fNL_4_vd_ortho,
        pk_tot_fNL_4_dd_ortho,
        pk_tot_fNL_0_b1b2_ortho, pk_tot_fNL_0_b2_ortho,
        pk_tot_fNL_0_b1bG2_ortho, pk_tot_fNL_0_bG2_ortho,
        pk_tot_fNL_2_b1b2_ortho, pk_tot_fNL_2_b2_ortho,
        pk_tot_fNL_2_b1bG2_ortho, pk_tot_fNL_2_bG2_ortho,
        pk_tot_fNL_4_b1b2_ortho, pk_tot_fNL_4_b2_ortho,
        pk_tot_fNL_4_b1bG2_ortho, pk_tot_fNL_4_bG2_ortho
    };

    for (c = 0; c < _NLPT_NUM_COMPONENTS_; c++)
        *outputs[c] = pk_buf[c];

    return _SUCCESS_;
}


int nonlinear_pt_pk_mult_at_kvec_and_z(
    struct nonlinear_pt * pnlpt,
    double * kvec,
    int n_k,
    double z,
    double * pk_mult
) {
    int i_z, c, ik, nk;
    int last_index;
    int offset;
    double val;
    double * data = NULL;
    double * spline_c = NULL;
    double * ln_kvec = NULL;

    nk = pnlpt->ln_k_size;

    for (ik = 0; ik < n_k; ik++) {
        class_test((kvec[ik] < exp(pnlpt->ln_k[0])) || (kvec[ik] > exp(pnlpt->ln_k[nk-1])),
                   pnlpt->error_message,
                   "nonlinear_pt_pk_mult_at_kvec_and_z: k=%e out of bounds [%e:%e]",
                   kvec[ik], 0., exp(pnlpt->ln_k[nk-1]));
    }

    for (i_z = 0; i_z < pnlpt->z_pk_num; i_z++)
        if (pnlpt->z_pk[i_z] == z) break;
    if (i_z == pnlpt->z_pk_num) {
        sprintf(pnlpt->error_message,
                "nonlinear_pt_pk_mult_at_kvec_and_z: z=%e not in z_pk list", z);
        return _FAILURE_;
    }

    offset = i_z * pnlpt->k_size;

    class_alloc(ln_kvec, n_k * sizeof(double), pnlpt->error_message);
    for (ik = 0; ik < n_k; ik++)
        ln_kvec[ik] = log(kvec[ik]);

    class_alloc(data, nk * sizeof(double), pnlpt->error_message);
    class_alloc(spline_c, nk * sizeof(double), pnlpt->error_message);

    for (c = 0; c < _NLPT_NUM_COMPONENTS_; c++) {
        double * src = _nlpt_component_array(pnlpt, c);
        for (ik = 0; ik < nk; ik++)
            data[ik] = src[offset + ik];

        class_call(array_spline_table_lines(
                       pnlpt->ln_k, nk, data, 1, spline_c,
                       _SPLINE_NATURAL_, pnlpt->error_message),
                   pnlpt->error_message, pnlpt->error_message);

        last_index = 0;
        for (ik = 0; ik < n_k; ik++) {
            class_call(array_interpolate_spline(
                           pnlpt->ln_k, nk, data, spline_c, 1,
                           ln_kvec[ik], &last_index, &val, 1,
                           pnlpt->error_message),
                       pnlpt->error_message, pnlpt->error_message);
            pk_mult[c * n_k + ik] = exp(val);
        }
    }

    free(data);
    free(spline_c);
    free(ln_kvec);

    return _SUCCESS_;
}
