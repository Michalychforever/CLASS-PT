/** @file nonlinear_pt_output.c
 *
 * Output/interpolation routines for the PT power spectra.
 * Ported from the v2 spectra.c functions spectra_pk_nl_bias_at_z_i()
 * and spectra_pk_nl_at_k_and_z() to work with the v3 CLASS codebase.
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
 */
static double * _nlpt_component_array(struct nonlinear_pt *pnlpt, int c) {
    switch (c) {
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
    case 15: return pnlpt->ln_pk_Tree_0_vv;
    case 16: return pnlpt->ln_pk_Tree_0_vd;
    case 17: return pnlpt->ln_pk_Tree_0_dd;
    case 18: return pnlpt->ln_pk_Tree_2_vv;
    case 19: return pnlpt->ln_pk_Tree_2_vd;
    case 20: return pnlpt->ln_pk_Tree_4_vv;
    case 21: return pnlpt->ln_pk_0_vv;
    case 22: return pnlpt->ln_pk_0_vd;
    case 23: return pnlpt->ln_pk_0_dd;
    case 24: return pnlpt->ln_pk_2_vv;
    case 25: return pnlpt->ln_pk_2_vd;
    case 26: return pnlpt->ln_pk_2_dd;
    case 27: return pnlpt->ln_pk_4_vv;
    case 28: return pnlpt->ln_pk_4_vd;
    case 29: return pnlpt->ln_pk_4_dd;
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
    /* fNL equilateral */
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
    /* fNL orthogonal */
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

int nonlinear_pt_bias_at_z_i(
                       struct nonlinear_pt *pnlpt,
                       enum linear_or_logarithmic mode,
                       int i_z, /* number(s) of redshift at which P(k,z) and T_i(k,z) should be written*/
                       double * output_tot, /* array with argument output_tot[index_k] (must be already allocated) */
                       double * output_tot_Id2d2, //1
                       double * output_tot_Id2, //2
                            double * output_tot_IG2, //3
                            double * output_tot_Id2G2, //4
                            double * output_tot_IG2G2, //5
                            double * output_tot_IFG2, //6
                              double * output_tot_IFG2_0b1, //7
                              double * output_tot_IFG2_0, //8
                              double * output_tot_IFG2_2, //9
                            double * output_tot_CTR, //10
                              double * output_tot_CTR_0, //11
                              double * output_tot_CTR_2, //12
                              double * output_tot_CTR_4, //13
				double * output_tot_Tree, //14
                              double * output_tot_Tree_0_vv, //15
                              double * output_tot_Tree_0_vd, //16
                              double * output_tot_Tree_0_dd, //17
                              double * output_tot_Tree_2_vv, //18
                              double * output_tot_Tree_2_vd, //19
                              double * output_tot_Tree_4_vv, //20
                  double * output_tot_0_vv, //21
                              double * output_tot_0_vd, //22
                              double * output_tot_0_dd, //23
                              double * output_tot_2_vv, //24
                              double * output_tot_2_vd, //25
                              double * output_tot_2_dd, //26
                              double * output_tot_4_vv, //27
                              double * output_tot_4_vd, //28
                              double * output_tot_4_dd, //29
                              double * output_tot_0_b1b2, //30
                              double * output_tot_0_b2, //31
                              double * output_tot_0_b1bG2, //32
                              double * output_tot_0_bG2, //33
                              double * output_tot_2_b1b2, //34
                              double * output_tot_2_b2, //35
                              double * output_tot_2_b1bG2, //36
                              double * output_tot_2_bG2, //37
                              double * output_tot_4_b2, //38
                              double * output_tot_4_bG2, //39
                              double * output_tot_4_b1b2, //40
                              double * output_tot_4_b1bG2, //41
                              double * output_tot_2_b2b2, //42
                              double * output_tot_2_b2bG2, //43
                              double * output_tot_2_bG2bG2, //44
                              double * output_tot_4_b2b2, //45
                              double * output_tot_4_b2bG2, //46
                              double * output_tot_4_bG2bG2, //47
                              double * output_tot_pk_nl_fNL, //48
                              double * output_tot_pk_fNLd2, //49
                              double * output_tot_pk_fNLG2, //50
                              double * output_tot_pk_l_fNL_0_vv, //51
                              double * output_tot_pk_l_fNL_0_vd, //52
                              double * output_tot_pk_l_fNL_0_dd, //53
                              double * output_tot_pk_l_fNL_2_vv, //54
                              double * output_tot_pk_l_fNL_2_vd, //55
                              double * output_tot_pk_l_fNL_2_dd, //56
                              double * output_tot_pk_l_fNL_4_vv, //57
                              double * output_tot_pk_l_fNL_4_vd, //58
                              double * output_tot_pk_l_fNL_4_dd, //59
                              double * output_tot_pk12_l_0_b1b2, //60
                              double * output_tot_pk12_l_0_b2, //61
                              double * output_tot_pk12_l_0_b1bG2, //62
                              double * output_tot_pk12_l_0_bG2, //63
                              double * output_tot_pk12_l_2_b1b2, //64
                              double * output_tot_pk12_l_2_b2, //65
                              double * output_tot_pk12_l_2_b1bG2, //66
                              double * output_tot_pk12_l_2_bG2, //67
                              double * output_tot_pk12_l_4_b1b2, //68
                              double * output_tot_pk12_l_4_b2, //69
                              double * output_tot_pk12_l_4_b1bG2, //70
                              double * output_tot_pk12_l_4_bG2, //71
                              double * output_tot_pk_nl_fNL_ortho, //72
                              double * output_tot_pk_fNLd2_ortho, //73
                              double * output_tot_pk_fNLG2_ortho, //74
                              double * output_tot_pk_l_fNL_0_vv_ortho, //75
                              double * output_tot_pk_l_fNL_0_vd_ortho, //76
                              double * output_tot_pk_l_fNL_0_dd_ortho, //77
                              double * output_tot_pk_l_fNL_2_vv_ortho, //78
                              double * output_tot_pk_l_fNL_2_vd_ortho, //79
                              double * output_tot_pk_l_fNL_2_dd_ortho, //80
                              double * output_tot_pk_l_fNL_4_vv_ortho, //81
                              double * output_tot_pk_l_fNL_4_vd_ortho, //82
                              double * output_tot_pk_l_fNL_4_dd_ortho, //83
                              double * output_tot_pk12_l_0_b1b2_ortho, //84
                              double * output_tot_pk12_l_0_b2_ortho, //85
                              double * output_tot_pk12_l_0_b1bG2_ortho, //86
                              double * output_tot_pk12_l_0_bG2_ortho, //87
                              double * output_tot_pk12_l_2_b1b2_ortho, //88
                              double * output_tot_pk12_l_2_b2_ortho, //89
                              double * output_tot_pk12_l_2_b1bG2_ortho, //90
                              double * output_tot_pk12_l_2_bG2_ortho, //91
                              double * output_tot_pk12_l_4_b1b2_ortho, //92
                              double * output_tot_pk12_l_4_b2_ortho, //93
                              double * output_tot_pk12_l_4_b1bG2_ortho, //94
                              double * output_tot_pk12_l_4_bG2_ortho //95
                       ) {

    int index_k, c;
    int nk = pnlpt->ln_k_size;
    int offset = i_z * pnlpt->k_size;

    /* Pack all 96 output pointers into an array for looping */
    double * outputs[_NLPT_NUM_COMPONENTS_] = {
        output_tot, output_tot_Id2d2, output_tot_Id2, output_tot_IG2,
        output_tot_Id2G2, output_tot_IG2G2, output_tot_IFG2,
        output_tot_IFG2_0b1, output_tot_IFG2_0, output_tot_IFG2_2,
        output_tot_CTR, output_tot_CTR_0, output_tot_CTR_2, output_tot_CTR_4,
        output_tot_Tree,
        output_tot_Tree_0_vv, output_tot_Tree_0_vd, output_tot_Tree_0_dd,
        output_tot_Tree_2_vv, output_tot_Tree_2_vd, output_tot_Tree_4_vv,
        output_tot_0_vv, output_tot_0_vd, output_tot_0_dd,
        output_tot_2_vv, output_tot_2_vd, output_tot_2_dd,
        output_tot_4_vv, output_tot_4_vd, output_tot_4_dd,
        output_tot_0_b1b2, output_tot_0_b2, output_tot_0_b1bG2, output_tot_0_bG2,
        output_tot_2_b1b2, output_tot_2_b2, output_tot_2_b1bG2, output_tot_2_bG2,
        output_tot_4_b2, output_tot_4_bG2, output_tot_4_b1b2, output_tot_4_b1bG2,
        output_tot_2_b2b2, output_tot_2_b2bG2, output_tot_2_bG2bG2,
        output_tot_4_b2b2, output_tot_4_b2bG2, output_tot_4_bG2bG2,
        /* fNL equilateral */
        output_tot_pk_nl_fNL, output_tot_pk_fNLd2, output_tot_pk_fNLG2,
        output_tot_pk_l_fNL_0_vv, output_tot_pk_l_fNL_0_vd, output_tot_pk_l_fNL_0_dd,
        output_tot_pk_l_fNL_2_vv, output_tot_pk_l_fNL_2_vd, output_tot_pk_l_fNL_2_dd,
        output_tot_pk_l_fNL_4_vv, output_tot_pk_l_fNL_4_vd, output_tot_pk_l_fNL_4_dd,
        output_tot_pk12_l_0_b1b2, output_tot_pk12_l_0_b2,
        output_tot_pk12_l_0_b1bG2, output_tot_pk12_l_0_bG2,
        output_tot_pk12_l_2_b1b2, output_tot_pk12_l_2_b2,
        output_tot_pk12_l_2_b1bG2, output_tot_pk12_l_2_bG2,
        output_tot_pk12_l_4_b1b2, output_tot_pk12_l_4_b2,
        output_tot_pk12_l_4_b1bG2, output_tot_pk12_l_4_bG2,
        /* fNL orthogonal */
        output_tot_pk_nl_fNL_ortho, output_tot_pk_fNLd2_ortho, output_tot_pk_fNLG2_ortho,
        output_tot_pk_l_fNL_0_vv_ortho, output_tot_pk_l_fNL_0_vd_ortho,
        output_tot_pk_l_fNL_0_dd_ortho,
        output_tot_pk_l_fNL_2_vv_ortho, output_tot_pk_l_fNL_2_vd_ortho,
        output_tot_pk_l_fNL_2_dd_ortho,
        output_tot_pk_l_fNL_4_vv_ortho, output_tot_pk_l_fNL_4_vd_ortho,
        output_tot_pk_l_fNL_4_dd_ortho,
        output_tot_pk12_l_0_b1b2_ortho, output_tot_pk12_l_0_b2_ortho,
        output_tot_pk12_l_0_b1bG2_ortho, output_tot_pk12_l_0_bG2_ortho,
        output_tot_pk12_l_2_b1b2_ortho, output_tot_pk12_l_2_b2_ortho,
        output_tot_pk12_l_2_b1bG2_ortho, output_tot_pk12_l_2_bG2_ortho,
        output_tot_pk12_l_4_b1b2_ortho, output_tot_pk12_l_4_b2_ortho,
        output_tot_pk12_l_4_b1bG2_ortho, output_tot_pk12_l_4_bG2_ortho
    };

    /* Copy data from struct arrays (logarithmic format) */
    for (c = 0; c < _NLPT_NUM_COMPONENTS_; c++) {
        double * src = _nlpt_component_array(pnlpt, c);
        for (index_k = 0; index_k < nk; index_k++)
            outputs[c][index_k] = src[offset + index_k];
    }

    /* Convert to linear format if requested */
    if (mode == linear) {
        for (c = 0; c < _NLPT_NUM_COMPONENTS_; c++) {
            for (index_k = 0; index_k < nk; index_k++)
                outputs[c][index_k] = exp(outputs[c][index_k]);
        }
    }

    return _SUCCESS_;
}


int nonlinear_pt_pk_at_k_and_z(
                             struct nonlinear_pt * pnlpt,
                             double k,
                             double z,
                             double * pk_tot, /* pointer to a single number (must be already allocated) */
                             double * pk_tot_Id2d2, //1
                             double * pk_tot_Id2, //2
                             double * pk_tot_IG2, //3
                             double * pk_tot_Id2G2, //4
                             double * pk_tot_IG2G2, //5
                             double * pk_tot_IFG2, //6
                             double * pk_tot_IFG2_0b1, //7
                             double * pk_tot_IFG2_0, //8
                             double * pk_tot_IFG2_2, //9
                             double * pk_tot_CTR, //10
                             double * pk_tot_CTR_0, //11
                             double * pk_tot_CTR_2, //12
                             double * pk_tot_CTR_4, //13
                             double * pk_tot_Tree, //14
                             double * pk_tot_Tree_0_vv, //15
                             double * pk_tot_Tree_0_vd, //16
                             double * pk_tot_Tree_0_dd, //17
                             double * pk_tot_Tree_2_vv, //18
                             double * pk_tot_Tree_2_vd, //19
                             double * pk_tot_Tree_4_vv, //20
                             double * pk_tot_0_vv, //21
                             double * pk_tot_0_vd, //22
                             double * pk_tot_0_dd, //23
                             double * pk_tot_2_vv, //24
                             double * pk_tot_2_vd, //25
                             double * pk_tot_2_dd, //26
                             double * pk_tot_4_vv, //27
                             double * pk_tot_4_vd, //28
                             double * pk_tot_4_dd, //29
                             double * pk_tot_0_b1b2, //30
                             double * pk_tot_0_b2, //31
                             double * pk_tot_0_b1bG2, //32
                             double * pk_tot_0_bG2, //33
                             double * pk_tot_2_b1b2, //34
                             double * pk_tot_2_b2, //35
                             double * pk_tot_2_b1bG2, //36
                             double * pk_tot_2_bG2, //37
                             double * pk_tot_4_b2, //38
                             double * pk_tot_4_bG2, //39
                             double * pk_tot_4_b1b2, //40
                             double * pk_tot_4_b1bG2, //41
                             double * pk_tot_2_b2b2, //42
                             double * pk_tot_2_b2bG2, //43
                             double * pk_tot_2_bG2bG2, //44
                             double * pk_tot_4_b2b2, //45
                             double * pk_tot_4_b2bG2, //46
                             double * pk_tot_4_bG2bG2, //47
                             double * pk_tot_fNL, //48
                             double * pk_tot_fNLd2, //49
                             double * pk_tot_fNLG2, //50
                             double * pk_tot_fNL_0_vv, //51
                             double * pk_tot_fNL_0_vd, //52
                             double * pk_tot_fNL_0_dd, //53
                             double * pk_tot_fNL_2_vv, //54
                             double * pk_tot_fNL_2_vd, //55
                             double * pk_tot_fNL_2_dd, //56
                             double * pk_tot_fNL_4_vv, //57
                             double * pk_tot_fNL_4_vd, //58
                             double * pk_tot_fNL_4_dd, //59
                             double * pk_tot_fNL_0_b1b2, //60
                             double * pk_tot_fNL_0_b2, //61
                             double * pk_tot_fNL_0_b1bG2, //62
                             double * pk_tot_fNL_0_bG2, //63
                             double * pk_tot_fNL_2_b1b2, //64
                             double * pk_tot_fNL_2_b2, //65
                             double * pk_tot_fNL_2_b1bG2, //66
                             double * pk_tot_fNL_2_bG2, //67
                             double * pk_tot_fNL_4_b1b2, //68
                             double * pk_tot_fNL_4_b2, //69
                             double * pk_tot_fNL_4_b1bG2, //70
                             double * pk_tot_fNL_4_bG2, //71
                             double * pk_tot_fNL_ortho, //72
                             double * pk_tot_fNLd2_ortho, //73
                             double * pk_tot_fNLG2_ortho, //74
                             double * pk_tot_fNL_0_vv_ortho, //75
                             double * pk_tot_fNL_0_vd_ortho, //76
                             double * pk_tot_fNL_0_dd_ortho, //77
                             double * pk_tot_fNL_2_vv_ortho, //78
                             double * pk_tot_fNL_2_vd_ortho, //79
                             double * pk_tot_fNL_2_dd_ortho, //80
                             double * pk_tot_fNL_4_vv_ortho, //81
                             double * pk_tot_fNL_4_vd_ortho, //82
                             double * pk_tot_fNL_4_dd_ortho, //83
                             double * pk_tot_fNL_0_b1b2_ortho, //84
                             double * pk_tot_fNL_0_b2_ortho, //85
                             double * pk_tot_fNL_0_b1bG2_ortho, //86
                             double * pk_tot_fNL_0_bG2_ortho, //87
                             double * pk_tot_fNL_2_b1b2_ortho, //88
                             double * pk_tot_fNL_2_b2_ortho, //89
                             double * pk_tot_fNL_2_b1bG2_ortho, //90
                             double * pk_tot_fNL_2_bG2_ortho, //91
                             double * pk_tot_fNL_4_b1b2_ortho, //92
                             double * pk_tot_fNL_4_b2_ortho, //93
                             double * pk_tot_fNL_4_b1bG2_ortho, //94
                             double * pk_tot_fNL_4_bG2_ortho //95
                             ) {

    double pk_buf[_NLPT_NUM_COMPONENTS_];

    /* Check that k is in valid range */
    class_test((k < exp(pnlpt->ln_k[0])) || (k > exp(pnlpt->ln_k[pnlpt->ln_k_size-1])),
               pnlpt->error_message,
               "k=%e out of bounds [%e:%e]",k,0.,exp(pnlpt->ln_k[pnlpt->ln_k_size-1]));

    /* Delegate to the batch function for a single k */
    class_call(nonlinear_pt_pk_mult_at_kvec_and_z(pnlpt, &k, 1, z, pk_buf),
               pnlpt->error_message, pnlpt->error_message);

    /* Distribute results to individual output pointers */
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
        /* fNL equilateral */
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
        /* fNL orthogonal */
        pk_tot_fNL_ortho, pk_tot_fNLd2_ortho, pk_tot_fNLG2_ortho,
        pk_tot_fNL_0_vv_ortho, pk_tot_fNL_0_vd_ortho, pk_tot_fNL_0_dd_ortho,
        pk_tot_fNL_2_vv_ortho, pk_tot_fNL_2_vd_ortho, pk_tot_fNL_2_dd_ortho,
        pk_tot_fNL_4_vv_ortho, pk_tot_fNL_4_vd_ortho, pk_tot_fNL_4_dd_ortho,
        pk_tot_fNL_0_b1b2_ortho, pk_tot_fNL_0_b2_ortho,
        pk_tot_fNL_0_b1bG2_ortho, pk_tot_fNL_0_bG2_ortho,
        pk_tot_fNL_2_b1b2_ortho, pk_tot_fNL_2_b2_ortho,
        pk_tot_fNL_2_b1bG2_ortho, pk_tot_fNL_2_bG2_ortho,
        pk_tot_fNL_4_b1b2_ortho, pk_tot_fNL_4_b2_ortho,
        pk_tot_fNL_4_b1bG2_ortho, pk_tot_fNL_4_bG2_ortho
    };

    int c;
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
    double val;
    double * data = NULL;
    double * spline_c = NULL;
    double * ln_kvec = NULL;

    nk = pnlpt->ln_k_size;

    /* Check that all k values are in valid range */
    for (ik = 0; ik < n_k; ik++) {
        class_test((kvec[ik] < exp(pnlpt->ln_k[0])) || (kvec[ik] > exp(pnlpt->ln_k[nk-1])),
                   pnlpt->error_message,
                   "nonlinear_pt_pk_mult_at_kvec_and_z: k=%e out of bounds [%e:%e]",
                   kvec[ik], 0., exp(pnlpt->ln_k[nk-1]));
    }

    /* Find z-index */
    for (i_z = 0; i_z < pnlpt->z_pk_num; i_z++) {
        if (pnlpt->z_pk[i_z] == z) break;
    }
    if (i_z == pnlpt->z_pk_num) {
        sprintf(pnlpt->error_message,
                "nonlinear_pt_pk_mult_at_kvec_and_z: z=%e not in z_pk list", z);
        return _FAILURE_;
    }

    int offset = i_z * pnlpt->k_size;

    /* Pre-compute ln(k) array */
    class_alloc(ln_kvec, n_k * sizeof(double), pnlpt->error_message);
    for (ik = 0; ik < n_k; ik++)
        ln_kvec[ik] = log(kvec[ik]);

    /* Allocate working buffers: one for data copy, one for spline coefficients */
    class_alloc(data, nk * sizeof(double), pnlpt->error_message);
    class_alloc(spline_c, nk * sizeof(double), pnlpt->error_message);

    /* For each component: copy from struct, compute spline, evaluate at all k */
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
