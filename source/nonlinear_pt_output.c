/** @file nonlinear_pt_output.c
 *
 * Output/interpolation routines for the PT power spectra.
 * Ported from the v2 spectra.c functions spectra_pk_nl_bias_at_z_i()
 * and spectra_pk_nl_at_k_and_z() to work with the v3 CLASS codebase.
 *
 * These functions provide:
 * - nonlinear_pt_bias_at_z_i(): Copy PT spectra from pnlpt arrays at a given z index
 * - nonlinear_pt_pk_at_k_and_z(): Interpolate PT spectra to arbitrary k at a given z
 */

#include "nonlinear_pt.h"

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
                              double * output_tot_4_bG2bG2, //47 //GC!
                              //GC -> account for the +1... 48+3+9+12=72...
                              double * output_tot_pk_nl_fNL, //49
                              double * output_tot_pk_fNLd2, //50
                              double * output_tot_pk_fNLG2, //51
                              //GC -> in real space...
                              double * output_tot_pk_l_fNL_0_vv, //52
                              double * output_tot_pk_l_fNL_0_vd, //53
                              double * output_tot_pk_l_fNL_0_dd, //54
                              double * output_tot_pk_l_fNL_2_vv, //55
                              double * output_tot_pk_l_fNL_2_vd, //56
                              double * output_tot_pk_l_fNL_2_dd, //57
                              double * output_tot_pk_l_fNL_4_vv, //58
                              double * output_tot_pk_l_fNL_4_vd, //59
                              double * output_tot_pk_l_fNL_4_dd, //60
                              //GC -> multipoles with matter...
                              double * output_tot_pk12_l_0_b1b2, //61
                              double * output_tot_pk12_l_0_b2, //62
                              double * output_tot_pk12_l_0_b1bG2, //63
                              double * output_tot_pk12_l_0_bG2, //64
                              double * output_tot_pk12_l_2_b1b2, //65
                              double * output_tot_pk12_l_2_b2, //66
                              double * output_tot_pk12_l_2_b1bG2, //67
                              double * output_tot_pk12_l_2_bG2, //68
                              double * output_tot_pk12_l_4_b1b2, //69
                              double * output_tot_pk12_l_4_b2, //70
                              double * output_tot_pk12_l_4_b1bG2, //71
                              double * output_tot_pk12_l_4_bG2, //72 //GC: ORTHOGONAL...
                              //GC -> multipoles with biased tracers...
                                //GC: ORTHOGONAL -- start
                              double * output_tot_pk_nl_fNL_ortho, //73
                              double * output_tot_pk_fNLd2_ortho, //74
                              double * output_tot_pk_fNLG2_ortho, //75
                              //GC -> in real space...
                              double * output_tot_pk_l_fNL_0_vv_ortho, //76
                              double * output_tot_pk_l_fNL_0_vd_ortho, //77
                              double * output_tot_pk_l_fNL_0_dd_ortho, //78
                              double * output_tot_pk_l_fNL_2_vv_ortho, //79
                              double * output_tot_pk_l_fNL_2_vd_ortho, //80
                              double * output_tot_pk_l_fNL_2_dd_ortho, //81
                              double * output_tot_pk_l_fNL_4_vv_ortho, //82
                              double * output_tot_pk_l_fNL_4_vd_ortho, //83
                              double * output_tot_pk_l_fNL_4_dd_ortho, //84
                              //GC -> multipoles with matter...
                              double * output_tot_pk12_l_0_b1b2_ortho, //85
                              double * output_tot_pk12_l_0_b2_ortho, //86
                              double * output_tot_pk12_l_0_b1bG2_ortho, //87
                              double * output_tot_pk12_l_0_bG2_ortho, //88
                              double * output_tot_pk12_l_2_b1b2_ortho, //89
                              double * output_tot_pk12_l_2_b2_ortho, //90
                              double * output_tot_pk12_l_2_b1bG2_ortho, //91
                              double * output_tot_pk12_l_2_bG2_ortho, //92
                              double * output_tot_pk12_l_4_b1b2_ortho, //93
                              double * output_tot_pk12_l_4_b2_ortho, //94
                              double * output_tot_pk12_l_4_b1bG2_ortho, //95
                              double * output_tot_pk12_l_4_bG2_ortho // 96: which is right, 72 + 24 (which is 3+9+12)...
                                //GC: ORTHOGONAL -- finish
                       ) {
    
    /** Summary: */
    
    /** - define local variables */
    
    int last_index;
    int last_index2;
    int last_index3;
    int last_index4;
    int last_index5;
    int last_index6;
    int last_index7;
    int last_index8;
    int last_index9;
    int index_k;
    double tau,ln_tau;
    
    /** - second step: for both modes (linear or logarithmic), store the spectrum in logarithmic format in the output array(s) */
    
    /** - --> (a) if only values at tau=tau_today are stored and we want P(k,z=0), no need to interpolate */
      
        for (index_k=0; index_k<pnlpt->ln_k_size; index_k++) {
            output_tot[index_k] = pnlpt->ln_pk_nl[i_z*pnlpt->k_size+index_k];
            output_tot_Id2d2[index_k] = pnlpt->ln_pk_Id2d2[i_z*pnlpt->k_size+index_k];
            output_tot_Id2[index_k] = pnlpt->ln_pk_Id2[i_z*pnlpt->k_size+index_k];
            output_tot_IG2[index_k] = pnlpt->ln_pk_IG2[i_z*pnlpt->k_size+index_k];
            output_tot_Id2G2[index_k] = pnlpt->ln_pk_Id2G2[i_z*pnlpt->k_size+index_k];
            output_tot_IG2G2[index_k] = pnlpt->ln_pk_IG2G2[i_z*pnlpt->k_size+index_k];
            output_tot_IFG2[index_k] = pnlpt->ln_pk_IFG2[i_z*pnlpt->k_size+index_k];
            output_tot_IFG2_0b1[index_k] = pnlpt->ln_pk_IFG2_0b1[i_z*pnlpt->k_size+index_k];
            output_tot_IFG2_0[index_k] = pnlpt->ln_pk_IFG2_0[i_z*pnlpt->k_size+index_k];
            output_tot_IFG2_2[index_k] = pnlpt->ln_pk_IFG2_2[i_z*pnlpt->k_size+index_k];
            output_tot_CTR[index_k] = pnlpt->ln_pk_CTR[i_z*pnlpt->k_size+index_k];
            output_tot_CTR_0[index_k] = pnlpt->ln_pk_CTR_0[i_z*pnlpt->k_size+index_k];
            output_tot_CTR_2[index_k] = pnlpt->ln_pk_CTR_2[i_z*pnlpt->k_size+index_k];
            output_tot_CTR_4[index_k] = pnlpt->ln_pk_CTR_4[i_z*pnlpt->k_size+index_k];
            output_tot_Tree[index_k] = pnlpt->ln_pk_Tree[i_z*pnlpt->k_size+index_k];
            
            output_tot_Tree_0_vv[index_k] = pnlpt->ln_pk_Tree_0_vv[i_z*pnlpt->k_size+index_k];
            output_tot_Tree_0_vd[index_k] = pnlpt->ln_pk_Tree_0_vd[i_z*pnlpt->k_size+index_k];
            output_tot_Tree_0_dd[index_k] = pnlpt->ln_pk_Tree_0_dd[i_z*pnlpt->k_size+index_k];
            output_tot_Tree_2_vv[index_k] = pnlpt->ln_pk_Tree_2_vv[i_z*pnlpt->k_size+index_k];
            output_tot_Tree_2_vd[index_k] = pnlpt->ln_pk_Tree_2_vd[i_z*pnlpt->k_size+index_k];
            output_tot_Tree_4_vv[index_k] = pnlpt->ln_pk_Tree_4_vv[i_z*pnlpt->k_size+index_k];
            
            output_tot_0_vv[index_k] = pnlpt->ln_pk_0_vv[i_z*pnlpt->k_size+index_k];
            output_tot_0_vd[index_k] = pnlpt->ln_pk_0_vd[i_z*pnlpt->k_size+index_k];
            output_tot_0_dd[index_k] = pnlpt->ln_pk_0_dd[i_z*pnlpt->k_size+index_k];
            output_tot_2_vv[index_k] = pnlpt->ln_pk_2_vv[i_z*pnlpt->k_size+index_k];
            output_tot_2_vd[index_k] = pnlpt->ln_pk_2_vd[i_z*pnlpt->k_size+index_k];
            output_tot_2_dd[index_k] = pnlpt->ln_pk_2_dd[i_z*pnlpt->k_size+index_k];
            output_tot_4_vv[index_k] = pnlpt->ln_pk_4_vv[i_z*pnlpt->k_size+index_k];
            output_tot_4_vd[index_k] = pnlpt->ln_pk_4_vd[i_z*pnlpt->k_size+index_k];
            output_tot_4_dd[index_k] = pnlpt->ln_pk_4_dd[i_z*pnlpt->k_size+index_k];
            output_tot_0_b1b2[index_k] = pnlpt->ln_pk_0_b1b2[i_z*pnlpt->k_size+index_k];
            output_tot_0_b1bG2[index_k] = pnlpt->ln_pk_0_b1bG2[i_z*pnlpt->k_size+index_k];
            output_tot_0_b2[index_k] = pnlpt->ln_pk_0_b2[i_z*pnlpt->k_size+index_k];
            output_tot_0_bG2[index_k] = pnlpt->ln_pk_0_bG2[i_z*pnlpt->k_size+index_k];
            
            output_tot_2_b1b2[index_k] = pnlpt->ln_pk_2_b1b2[i_z*pnlpt->k_size+index_k];
            output_tot_2_b1bG2[index_k] = pnlpt->ln_pk_2_b1bG2[i_z*pnlpt->k_size+index_k];
            output_tot_2_b2[index_k] = pnlpt->ln_pk_2_b2[i_z*pnlpt->k_size+index_k];
            output_tot_2_bG2[index_k] = pnlpt->ln_pk_2_bG2[i_z*pnlpt->k_size+index_k];
            
            output_tot_4_b2[index_k] = pnlpt->ln_pk_4_b2[i_z*pnlpt->k_size+index_k];
            output_tot_4_bG2[index_k] = pnlpt->ln_pk_4_bG2[i_z*pnlpt->k_size+index_k];
            
            output_tot_4_b1b2[index_k] = pnlpt->ln_pk_4_b1b2[i_z*pnlpt->k_size+index_k];
            output_tot_4_b1bG2[index_k] = pnlpt->ln_pk_4_b1bG2[i_z*pnlpt->k_size+index_k];
            
            output_tot_2_b2b2[index_k] = pnlpt->ln_pk_Id2d2_2[i_z*pnlpt->k_size+index_k];
            output_tot_2_b2bG2[index_k] = pnlpt->ln_pk_Id2G2_2[i_z*pnlpt->k_size+index_k];
            output_tot_2_bG2bG2[index_k] = pnlpt->ln_pk_IG2G2_2[i_z*pnlpt->k_size+index_k];
            
            output_tot_4_b2b2[index_k] = pnlpt->ln_pk_Id2d2_4[i_z*pnlpt->k_size+index_k];
            output_tot_4_b2bG2[index_k] = pnlpt->ln_pk_Id2G2_4[i_z*pnlpt->k_size+index_k];
            output_tot_4_bG2bG2[index_k] = pnlpt->ln_pk_IG2G2_4[i_z*pnlpt->k_size+index_k];
            
            
           //GC!!!
            
            output_tot_pk_nl_fNL[index_k] = pnlpt->ln_pk_nl_fNL[i_z*pnlpt->k_size+index_k];
            output_tot_pk_fNLd2[index_k] = pnlpt->ln_pk_fNLd2[i_z*pnlpt->k_size+index_k];
            output_tot_pk_fNLG2[index_k] = pnlpt->ln_pk_fNLG2[i_z*pnlpt->k_size+index_k];
            
            
            //GC: ORTHOGONAL -- start


            output_tot_pk_nl_fNL_ortho[index_k] = pnlpt->ln_pk_nl_fNL_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk_fNLd2_ortho[index_k] = pnlpt->ln_pk_fNLd2_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk_fNLG2_ortho[index_k] = pnlpt->ln_pk_fNLG2_ortho[i_z*pnlpt->k_size+index_k];


            //GC: ORTHOGONAL -- finish


            //GC...
            
            output_tot_pk_l_fNL_0_vv[index_k] = pnlpt->ln_pk_fNL_0_vv[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_0_vd[index_k] = pnlpt->ln_pk_fNL_0_vd[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_0_dd[index_k] = pnlpt->ln_pk_fNL_0_dd[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_2_vv[index_k] = pnlpt->ln_pk_fNL_2_vv[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_2_vd[index_k] = pnlpt->ln_pk_fNL_2_vd[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_2_dd[index_k] = pnlpt->ln_pk_fNL_2_dd[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_4_vv[index_k] = pnlpt->ln_pk_fNL_4_vv[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_4_vd[index_k] = pnlpt->ln_pk_fNL_4_vd[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_4_dd[index_k] = pnlpt->ln_pk_fNL_4_dd[i_z*pnlpt->k_size+index_k];
            
            
            //GC: ORTHOGONAL -- start


            output_tot_pk_l_fNL_0_vv_ortho[index_k] = pnlpt->ln_pk_fNL_0_vv_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_0_vd_ortho[index_k] = pnlpt->ln_pk_fNL_0_vd_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_0_dd_ortho[index_k] = pnlpt->ln_pk_fNL_0_dd_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_2_vv_ortho[index_k] = pnlpt->ln_pk_fNL_2_vv_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_2_vd_ortho[index_k] = pnlpt->ln_pk_fNL_2_vd_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_2_dd_ortho[index_k] = pnlpt->ln_pk_fNL_2_dd_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_4_vv_ortho[index_k] = pnlpt->ln_pk_fNL_4_vv_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_4_vd_ortho[index_k] = pnlpt->ln_pk_fNL_4_vd_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk_l_fNL_4_dd_ortho[index_k] = pnlpt->ln_pk_fNL_4_dd_ortho[i_z*pnlpt->k_size+index_k];
            


            //GC: ORTHOGONAL -- finish

            
            //GC...
            
            output_tot_pk12_l_0_b1b2[index_k] = pnlpt->ln_pk12_0_b1b2[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_0_b2[index_k] = pnlpt->ln_pk12_0_b2[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_0_b1bG2[index_k] = pnlpt->ln_pk12_0_b1bG2[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_0_bG2[index_k] = pnlpt->ln_pk12_0_bG2[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_2_b1b2[index_k] = pnlpt->ln_pk12_2_b1b2[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_2_b2[index_k] = pnlpt->ln_pk12_2_b2[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_2_b1bG2[index_k] = pnlpt->ln_pk12_2_b1bG2[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_2_bG2[index_k] = pnlpt->ln_pk12_2_bG2[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_4_b1b2[index_k] = pnlpt->ln_pk12_4_b1b2[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_4_b2[index_k] = pnlpt->ln_pk12_4_b2[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_4_b1bG2[index_k] = pnlpt->ln_pk12_4_b1bG2[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_4_bG2[index_k] = pnlpt->ln_pk12_4_bG2[i_z*pnlpt->k_size+index_k];
            
            
            //GC: ORTHOGONAL -- start


            output_tot_pk12_l_0_b1b2_ortho[index_k] = pnlpt->ln_pk12_0_b1b2_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_0_b2_ortho[index_k] = pnlpt->ln_pk12_0_b2_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_0_b1bG2_ortho[index_k] = pnlpt->ln_pk12_0_b1bG2_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_0_bG2_ortho[index_k] = pnlpt->ln_pk12_0_bG2_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_2_b1b2_ortho[index_k] = pnlpt->ln_pk12_2_b1b2_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_2_b2_ortho[index_k] = pnlpt->ln_pk12_2_b2_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_2_b1bG2_ortho[index_k] = pnlpt->ln_pk12_2_b1bG2_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_2_bG2_ortho[index_k] = pnlpt->ln_pk12_2_bG2_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_4_b1b2_ortho[index_k] = pnlpt->ln_pk12_4_b1b2_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_4_b2_ortho[index_k] = pnlpt->ln_pk12_4_b2_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_4_b1bG2_ortho[index_k] = pnlpt->ln_pk12_4_b1bG2_ortho[i_z*pnlpt->k_size+index_k];
            output_tot_pk12_l_4_bG2_ortho[index_k] = pnlpt->ln_pk12_4_bG2_ortho[i_z*pnlpt->k_size+index_k];


            //GC: ORTHOGONAL -- finish


            
        }

    /** - fourth step: eventually convert to linear format */
    
    if (mode == linear) {
        for (index_k=0; index_k<pnlpt->ln_k_size; index_k++) {
            output_tot[index_k] = exp(output_tot[index_k]);
            output_tot_Id2d2[index_k] = exp(output_tot_Id2d2[index_k]);
            output_tot_Id2[index_k] = exp(output_tot_Id2[index_k]);
            output_tot_IG2[index_k] = exp(output_tot_IG2[index_k]);
            output_tot_Id2G2[index_k] = exp(output_tot_Id2G2[index_k]);
            output_tot_IG2G2[index_k] = exp(output_tot_IG2G2[index_k]);
            output_tot_IFG2[index_k] = exp(output_tot_IFG2[index_k]);
            output_tot_IFG2_0b1[index_k] = exp(output_tot_IFG2_0b1[index_k]);
            output_tot_IFG2_0[index_k] = exp(output_tot_IFG2_0[index_k]);
            output_tot_IFG2_2[index_k] = exp(output_tot_IFG2_2[index_k]);
            output_tot_CTR[index_k] = exp(output_tot_CTR[index_k]);
            output_tot_CTR_0[index_k] = exp(output_tot_CTR_0[index_k]);
            output_tot_CTR_2[index_k] = exp(output_tot_CTR_2[index_k]);
            output_tot_CTR_4[index_k] = exp(output_tot_CTR_4[index_k]);
            output_tot_Tree[index_k] = exp(output_tot_Tree[index_k]);
            output_tot_Tree_0_vv[index_k] = exp(output_tot_Tree_0_vv[index_k]);
            output_tot_Tree_0_vd[index_k] = exp(output_tot_Tree_0_vd[index_k]);
            output_tot_Tree_0_dd[index_k] = exp(output_tot_Tree_0_dd[index_k]);
            output_tot_Tree_2_vv[index_k] = exp(output_tot_Tree_2_vv[index_k]);
            output_tot_Tree_2_vd[index_k] = exp(output_tot_Tree_2_vd[index_k]);
            output_tot_Tree_4_vv[index_k] = exp(output_tot_Tree_4_vv[index_k]);
            output_tot_0_vv[index_k] = exp(output_tot_0_vv[index_k]);
            output_tot_0_vd[index_k] = exp(output_tot_0_vd[index_k]);
            output_tot_0_dd[index_k] = exp(output_tot_0_dd[index_k]);
            output_tot_2_vv[index_k] = exp(output_tot_2_vv[index_k]);
            output_tot_2_vd[index_k] = exp(output_tot_2_vd[index_k]);
            output_tot_2_dd[index_k] = exp(output_tot_2_dd[index_k]);
            output_tot_4_vv[index_k] = exp(output_tot_4_vv[index_k]);
            output_tot_4_vd[index_k] = exp(output_tot_4_vd[index_k]);
            output_tot_4_dd[index_k] = exp(output_tot_4_dd[index_k]);
            output_tot_0_b1b2[index_k] = exp(output_tot_0_b1b2[index_k]);
            output_tot_0_b1bG2[index_k] = exp(output_tot_0_b1bG2[index_k]);
            output_tot_0_b2[index_k] = exp(output_tot_0_b2[index_k]);
            output_tot_0_bG2[index_k] = exp(output_tot_0_bG2[index_k]);
            
            output_tot_2_b1b2[index_k] = exp(output_tot_2_b1b2[index_k]);
            output_tot_2_b1bG2[index_k] = exp(output_tot_2_b1bG2[index_k]);
            output_tot_2_b2[index_k] = exp(output_tot_2_b2[index_k]);
            output_tot_2_bG2[index_k] = exp(output_tot_2_bG2[index_k]);
            
            output_tot_4_b2[index_k] = exp(output_tot_4_b2[index_k]);
            output_tot_4_bG2[index_k] = exp(output_tot_4_bG2[index_k]);
            output_tot_4_b1b2[index_k] = exp(output_tot_4_b1b2[index_k]);
            output_tot_4_b1bG2[index_k] = exp(output_tot_4_b1bG2[index_k]);
            
            output_tot_2_b2b2[index_k] = exp(output_tot_2_b2b2[index_k]);
            output_tot_2_b2bG2[index_k] = exp(output_tot_2_b2bG2[index_k]);
            output_tot_2_bG2bG2[index_k] = exp(output_tot_2_bG2bG2[index_k]);
            
            output_tot_4_b2b2[index_k] = exp(output_tot_4_b2b2[index_k]);
            output_tot_4_b2bG2[index_k] = exp(output_tot_4_b2bG2[index_k]);
            output_tot_4_bG2bG2[index_k] = exp(output_tot_4_bG2bG2[index_k]);
            
            //GC!!!
            
            output_tot_pk_nl_fNL[index_k] = exp(output_tot_pk_nl_fNL[index_k]);
            output_tot_pk_fNLd2[index_k] = exp(output_tot_pk_fNLd2[index_k]);
            output_tot_pk_fNLG2[index_k] = exp(output_tot_pk_fNLG2[index_k]);
            
            
            //GC: ORTHOGONAL -- start


            output_tot_pk_nl_fNL_ortho[index_k] = exp(output_tot_pk_nl_fNL_ortho[index_k]);
            output_tot_pk_fNLd2_ortho[index_k] = exp(output_tot_pk_fNLd2_ortho[index_k]);
            output_tot_pk_fNLG2_ortho[index_k] = exp(output_tot_pk_fNLG2_ortho[index_k]);
            


            //GC: ORTHOGONAL -- finish


            
            //GC!
            
            output_tot_pk_l_fNL_0_vv[index_k] = exp(output_tot_pk_l_fNL_0_vv[index_k]);
            output_tot_pk_l_fNL_0_vd[index_k] = exp(output_tot_pk_l_fNL_0_vd[index_k]);
            output_tot_pk_l_fNL_0_dd[index_k] = exp(output_tot_pk_l_fNL_0_dd[index_k]);
            output_tot_pk_l_fNL_2_vv[index_k] = exp(output_tot_pk_l_fNL_2_vv[index_k]);
            output_tot_pk_l_fNL_2_vd[index_k] = exp(output_tot_pk_l_fNL_2_vd[index_k]);
            output_tot_pk_l_fNL_2_dd[index_k] = exp(output_tot_pk_l_fNL_2_dd[index_k]);
            output_tot_pk_l_fNL_4_vv[index_k] = exp(output_tot_pk_l_fNL_4_vv[index_k]);
            output_tot_pk_l_fNL_4_vd[index_k] = exp(output_tot_pk_l_fNL_4_vd[index_k]);
            output_tot_pk_l_fNL_4_dd[index_k] = exp(output_tot_pk_l_fNL_4_dd[index_k]);
            
            
            
            //GC: ORTHOGONAL -- start


            output_tot_pk_l_fNL_0_vv_ortho[index_k] = exp(output_tot_pk_l_fNL_0_vv_ortho[index_k]);
            output_tot_pk_l_fNL_0_vd_ortho[index_k] = exp(output_tot_pk_l_fNL_0_vd_ortho[index_k]);
            output_tot_pk_l_fNL_0_dd_ortho[index_k] = exp(output_tot_pk_l_fNL_0_dd_ortho[index_k]);
            output_tot_pk_l_fNL_2_vv_ortho[index_k] = exp(output_tot_pk_l_fNL_2_vv_ortho[index_k]);
            output_tot_pk_l_fNL_2_vd_ortho[index_k] = exp(output_tot_pk_l_fNL_2_vd_ortho[index_k]);
            output_tot_pk_l_fNL_2_dd_ortho[index_k] = exp(output_tot_pk_l_fNL_2_dd_ortho[index_k]);
            output_tot_pk_l_fNL_4_vv_ortho[index_k] = exp(output_tot_pk_l_fNL_4_vv_ortho[index_k]);
            output_tot_pk_l_fNL_4_vd_ortho[index_k] = exp(output_tot_pk_l_fNL_4_vd_ortho[index_k]);
            output_tot_pk_l_fNL_4_dd_ortho[index_k] = exp(output_tot_pk_l_fNL_4_dd_ortho[index_k]);


            //GC: ORTHOGONAL -- finish


            
            
            //GC!
            
            output_tot_pk12_l_0_b1b2[index_k] = exp(output_tot_pk12_l_0_b1b2[index_k]);
            output_tot_pk12_l_0_b2[index_k] = exp(output_tot_pk12_l_0_b2[index_k]);
            output_tot_pk12_l_0_b1bG2[index_k] = exp(output_tot_pk12_l_0_b1bG2[index_k]);
            output_tot_pk12_l_0_bG2[index_k] = exp(output_tot_pk12_l_0_bG2[index_k]);
            output_tot_pk12_l_2_b1b2[index_k] = exp(output_tot_pk12_l_2_b1b2[index_k]);
            output_tot_pk12_l_2_b2[index_k] = exp(output_tot_pk12_l_2_b2[index_k]);
            output_tot_pk12_l_2_b1bG2[index_k] = exp(output_tot_pk12_l_2_b1bG2[index_k]);
            output_tot_pk12_l_2_bG2[index_k] = exp(output_tot_pk12_l_2_bG2[index_k]);
            output_tot_pk12_l_4_b1b2[index_k] = exp(output_tot_pk12_l_4_b1b2[index_k]);
            output_tot_pk12_l_4_b2[index_k] = exp(output_tot_pk12_l_4_b2[index_k]);
            output_tot_pk12_l_4_b1bG2[index_k] = exp(output_tot_pk12_l_4_b1bG2[index_k]);
            output_tot_pk12_l_4_bG2[index_k] = exp(output_tot_pk12_l_4_bG2[index_k]);
            
            
            //GC: ORTHOGONAL -- start


            output_tot_pk12_l_0_b1b2_ortho[index_k] = exp(output_tot_pk12_l_0_b1b2_ortho[index_k]);
            output_tot_pk12_l_0_b2_ortho[index_k] = exp(output_tot_pk12_l_0_b2_ortho[index_k]);
            output_tot_pk12_l_0_b1bG2_ortho[index_k] = exp(output_tot_pk12_l_0_b1bG2_ortho[index_k]);
            output_tot_pk12_l_0_bG2_ortho[index_k] = exp(output_tot_pk12_l_0_bG2_ortho[index_k]);
            output_tot_pk12_l_2_b1b2_ortho[index_k] = exp(output_tot_pk12_l_2_b1b2_ortho[index_k]);
            output_tot_pk12_l_2_b2_ortho[index_k] = exp(output_tot_pk12_l_2_b2_ortho[index_k]);
            output_tot_pk12_l_2_b1bG2_ortho[index_k] = exp(output_tot_pk12_l_2_b1bG2_ortho[index_k]);
            output_tot_pk12_l_2_bG2_ortho[index_k] = exp(output_tot_pk12_l_2_bG2_ortho[index_k]);
            output_tot_pk12_l_4_b1b2_ortho[index_k] = exp(output_tot_pk12_l_4_b1b2_ortho[index_k]);
            output_tot_pk12_l_4_b2_ortho[index_k] = exp(output_tot_pk12_l_4_b2_ortho[index_k]);
            output_tot_pk12_l_4_b1bG2_ortho[index_k] = exp(output_tot_pk12_l_4_b1bG2_ortho[index_k]);
            output_tot_pk12_l_4_bG2_ortho[index_k] = exp(output_tot_pk12_l_4_bG2_ortho[index_k]);


            //GC: ORTHOGONAL -- finish

            

        }
    }
    
    return _SUCCESS_;
    
}


int nonlinear_pt_pk_at_k_and_z(
                             struct nonlinear_pt * pnlpt,
                             double k,
                             double z,
                             double * pk_tot, /* pointer to a single number (must be already allocated) */
                             //GC -> the first counts as well... The issue now is that I do not understand if I really only need to modify this function... The key is the function spectra_pk_nl_at_z -> it calls the spectra that... No, I need also spectra_pk_nl_bias_at_z_i since this is used below...
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
                             //GC!
                             double * pk_tot_fNL,
                             double * pk_tot_fNLd2,
                             double * pk_tot_fNLG2,
                             //GC!!!
                             double * pk_tot_fNL_0_vv,
                             double * pk_tot_fNL_0_vd,
                             double * pk_tot_fNL_0_dd,
                             double * pk_tot_fNL_2_vv,
                             double * pk_tot_fNL_2_vd,
                             double * pk_tot_fNL_2_dd,
                             double * pk_tot_fNL_4_vv,
                             double * pk_tot_fNL_4_vd,
                             double * pk_tot_fNL_4_dd,
                             //GC!
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
                             double * pk_tot_fNL_4_bG2, //GC: ORTHOGONAL...
                                //GC: ORTHOGONAL -- start
                             double * pk_tot_fNL_ortho,
                             double * pk_tot_fNLd2_ortho,
                             double * pk_tot_fNLG2_ortho,
                             //GC!!!
                             double * pk_tot_fNL_0_vv_ortho,
                             double * pk_tot_fNL_0_vd_ortho,
                             double * pk_tot_fNL_0_dd_ortho,
                             double * pk_tot_fNL_2_vv_ortho,
                             double * pk_tot_fNL_2_vd_ortho,
                             double * pk_tot_fNL_2_dd_ortho,
                             double * pk_tot_fNL_4_vv_ortho,
                             double * pk_tot_fNL_4_vd_ortho,
                             double * pk_tot_fNL_4_dd_ortho,
                             //GC!
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
                                //GC: ORTHOGONAL -- finish
                             ) {

  /** Summary: */

  /** - define local variables */
    

  int last_index;
    int last_index2;
    int last_index3;
    int last_index4;
    int last_index5;
    int last_index6;
    int last_index7;
    int last_index8;
    int last_index9;
    
    int last_index10;
    int last_index11;
    int last_index12;
    
    //GC -> I do not understand why you need all of these... For example, last_index12 is not used... For now, I keep last_index_fNL_, and then I will see what to do... It seems that only one index is needed... This makes sense with how it is coded in nonlinear_pt.c... But here the index is not even initialized... I will follow this...
    
    //int last_index_fNL_;
    int last_index_fNL;


  double * spectrum_at_z = NULL;
  double * spectrum_Id2d2_at_z = NULL;
    double * spectrum_Id2_at_z = NULL;
    double * spectrum_IG2_at_z = NULL;
      double * spectrum_Id2G2_at_z = NULL;
    double * spectrum_IG2G2_at_z = NULL;
    double * spectrum_IFG2_at_z = NULL;
    double * spectrum_IFG2_0b1_at_z = NULL;
    double * spectrum_IFG2_0_at_z = NULL;
    double * spectrum_IFG2_2_at_z = NULL;
     double * spectrum_CTR_at_z = NULL;
    double * spectrum_CTR_0_at_z = NULL;
    double * spectrum_CTR_2_at_z = NULL;
    double * spectrum_CTR_4_at_z = NULL;
     double * spectrum_Tree_at_z = NULL;
    double * spectrum_Tree_0_vv_at_z = NULL;
    double * spectrum_Tree_0_vd_at_z = NULL;
    double * spectrum_Tree_0_dd_at_z = NULL;
    double * spectrum_Tree_2_vv_at_z = NULL;
    double * spectrum_Tree_2_vd_at_z = NULL;
    double * spectrum_Tree_4_vv_at_z = NULL;
  double * spline;
  double * spline_Id2d2;
    double * spline_Id2;
    double * spline_IG2;
      double * spline_Id2G2;
     double * spline_IG2G2;
    double * spline_IFG2;
    double * spline_IFG2_0b1;
    double * spline_IFG2_0;
    double * spline_IFG2_2;
    double * spline_CTR;
    double * spline_CTR_2;
    double * spline_CTR_4;
    double * spline_CTR_0;
    double * spline_Tree;
    
    double * spline_Tree_0_vv;
    double * spline_Tree_0_vd;
    double * spline_Tree_0_dd;
    double * spline_Tree_2_vv;
    double * spline_Tree_2_vd;
    double * spline_Tree_4_vv;
    
    double * spectrum_0_vv_at_z = NULL;
    double * spline_0_vv;
    double * spectrum_0_vd_at_z = NULL;
    double * spline_0_vd;
    double * spectrum_0_dd_at_z = NULL;
    double * spline_0_dd;
    
    double * spectrum_2_vv_at_z = NULL;
    double * spline_2_vv;
    double * spectrum_2_vd_at_z = NULL;
    double * spline_2_vd;
    double * spectrum_2_dd_at_z = NULL;
    double * spline_2_dd;
    
    double * spectrum_4_vv_at_z = NULL;
    double * spline_4_vv;
    double * spectrum_4_vd_at_z = NULL;
    double * spline_4_vd;
    double * spectrum_4_dd_at_z = NULL;
    double * spline_4_dd;
    
    double * spectrum_0_b1b2_at_z = NULL;
    double * spline_0_b1b2;
    double * spectrum_0_b2_at_z = NULL;
    double * spline_0_b2;
    double * spectrum_0_b1bG2_at_z = NULL;
    double * spline_0_b1bG2;
    double * spectrum_0_bG2_at_z = NULL;
    double * spline_0_bG2;
    
    double * spectrum_2_b1b2_at_z = NULL;
    double * spline_2_b1b2;
    double * spectrum_2_b2_at_z = NULL;
    double * spline_2_b2;
    double * spectrum_2_b1bG2_at_z = NULL;
    double * spline_2_b1bG2;
    double * spectrum_2_bG2_at_z = NULL;
    double * spline_2_bG2;
    
    double * spectrum_4_b2_at_z = NULL;
    double * spline_4_b2;
    double * spectrum_4_bG2_at_z = NULL;
    double * spline_4_bG2;
    double * spectrum_4_b1b2_at_z = NULL;
    double * spline_4_b1b2;
    double * spectrum_4_b1bG2_at_z = NULL;
    double * spline_4_b1bG2;
    
    double * spectrum_2_b2b2_at_z = NULL;
    double * spline_2_b2b2;
    double * spectrum_2_b2bG2_at_z = NULL;
    double * spline_2_b2bG2;
    double * spectrum_2_bG2bG2_at_z = NULL;
    double * spline_2_bG2bG2;
    
    double * spectrum_4_b2b2_at_z = NULL;
    double * spline_4_b2b2;
    double * spectrum_4_b2bG2_at_z = NULL;
    double * spline_4_b2bG2;
    double * spectrum_4_bG2bG2_at_z = NULL;
    double * spline_4_bG2bG2;
    
    //GC!
    
    double * spectrum_fNL_at_z = NULL;
    double * spline_fNL;

    double * spectrum_fNLd2_at_z = NULL;
    double * spline_fNLd2;

    double * spectrum_fNLG2_at_z = NULL;
    double * spline_fNLG2;

                                 //GC: ORTHOGONAL -- start

                                 double * spectrum_fNL_at_z_ortho = NULL;
                                 double * spline_fNL_ortho;

                                 double * spectrum_fNLd2_at_z_ortho = NULL;
                                 double * spline_fNLd2_ortho;

                                 double * spectrum_fNLG2_at_z_ortho = NULL;
                                 double * spline_fNLG2_ortho;

                                 //GC: ORTHOGONAL -- finish

    
    //GC!
    
    double * spectrum_fNL_0_vv_at_z = NULL;
    double * spline_fNL_0_vv;

    double * spectrum_fNL_0_vd_at_z = NULL;
    double * spline_fNL_0_vd;

    double * spectrum_fNL_0_dd_at_z = NULL;
    double * spline_fNL_0_dd;

    double * spectrum_fNL_2_vv_at_z = NULL;
    double * spline_fNL_2_vv;

    double * spectrum_fNL_2_vd_at_z = NULL;
    double * spline_fNL_2_vd;

    double * spectrum_fNL_2_dd_at_z = NULL;
    double * spline_fNL_2_dd;

    double * spectrum_fNL_4_vv_at_z = NULL;
    double * spline_fNL_4_vv;

    double * spectrum_fNL_4_vd_at_z = NULL;
    double * spline_fNL_4_vd;

    double * spectrum_fNL_4_dd_at_z = NULL;
    double * spline_fNL_4_dd;
                                 
                                 
                                 //GC: ORTHOGONAL -- start


                                 double * spectrum_fNL_0_vv_at_z_ortho = NULL;
                                 double * spline_fNL_0_vv_ortho;

                                 double * spectrum_fNL_0_vd_at_z_ortho = NULL;
                                 double * spline_fNL_0_vd_ortho;

                                 double * spectrum_fNL_0_dd_at_z_ortho = NULL;
                                 double * spline_fNL_0_dd_ortho;

                                 double * spectrum_fNL_2_vv_at_z_ortho = NULL;
                                 double * spline_fNL_2_vv_ortho;

                                 double * spectrum_fNL_2_vd_at_z_ortho = NULL;
                                 double * spline_fNL_2_vd_ortho;

                                 double * spectrum_fNL_2_dd_at_z_ortho = NULL;
                                 double * spline_fNL_2_dd_ortho;

                                 double * spectrum_fNL_4_vv_at_z_ortho = NULL;
                                 double * spline_fNL_4_vv_ortho;

                                 double * spectrum_fNL_4_vd_at_z_ortho = NULL;
                                 double * spline_fNL_4_vd_ortho;

                                 double * spectrum_fNL_4_dd_at_z_ortho = NULL;
                                 double * spline_fNL_4_dd_ortho;


                                 //GC: ORTHOGONAL -- finish


        
    //GC!
    
    double * spectrum_fNL_0_b1b2_at_z = NULL;
    double * spline_fNL_0_b1b2;

    double * spectrum_fNL_0_b2_at_z = NULL;
    double * spline_fNL_0_b2;

    double * spectrum_fNL_0_b1bG2_at_z = NULL;
    double * spline_fNL_0_b1bG2;

    double * spectrum_fNL_0_bG2_at_z = NULL;
    double * spline_fNL_0_bG2;


    
    
    
    double * spectrum_fNL_2_b1b2_at_z = NULL;
    double * spline_fNL_2_b1b2;

    double * spectrum_fNL_2_b2_at_z = NULL;
    double * spline_fNL_2_b2;

    double * spectrum_fNL_2_b1bG2_at_z = NULL;
    double * spline_fNL_2_b1bG2;

    double * spectrum_fNL_2_bG2_at_z = NULL;
    double * spline_fNL_2_bG2;


    
    
    double * spectrum_fNL_4_b1b2_at_z = NULL;
    double * spline_fNL_4_b1b2;

    double * spectrum_fNL_4_b2_at_z = NULL;
    double * spline_fNL_4_b2;

    double * spectrum_fNL_4_b1bG2_at_z = NULL;
    double * spline_fNL_4_b1bG2;

    double * spectrum_fNL_4_bG2_at_z = NULL;
    double * spline_fNL_4_bG2;
                                 
                                 
                                 
                                 
                                 //GC: ORTHOGONAL -- start

                                 
                                 
                                 
                                 double * spectrum_fNL_0_b1b2_at_z_ortho = NULL;
                                 double * spline_fNL_0_b1b2_ortho;

                                 double * spectrum_fNL_0_b2_at_z_ortho = NULL;
                                 double * spline_fNL_0_b2_ortho;

                                 double * spectrum_fNL_0_b1bG2_at_z_ortho = NULL;
                                 double * spline_fNL_0_b1bG2_ortho;

                                 double * spectrum_fNL_0_bG2_at_z_ortho = NULL;
                                 double * spline_fNL_0_bG2_ortho;


                                 
                                 
                                 
                                 double * spectrum_fNL_2_b1b2_at_z_ortho = NULL;
                                 double * spline_fNL_2_b1b2_ortho;

                                 double * spectrum_fNL_2_b2_at_z_ortho = NULL;
                                 double * spline_fNL_2_b2_ortho;

                                 double * spectrum_fNL_2_b1bG2_at_z_ortho = NULL;
                                 double * spline_fNL_2_b1bG2_ortho;

                                 double * spectrum_fNL_2_bG2_at_z_ortho = NULL;
                                 double * spline_fNL_2_bG2_ortho;


                                 
                                 
                                 double * spectrum_fNL_4_b1b2_at_z_ortho = NULL;
                                 double * spline_fNL_4_b1b2_ortho;

                                 double * spectrum_fNL_4_b2_at_z_ortho = NULL;
                                 double * spline_fNL_4_b2_ortho;

                                 double * spectrum_fNL_4_b1bG2_at_z_ortho = NULL;
                                 double * spline_fNL_4_b1bG2_ortho;

                                 double * spectrum_fNL_4_bG2_at_z_ortho = NULL;
                                 double * spline_fNL_4_bG2_ortho;




                                 //GC: ORTHOGONAL -- finish


    
    


  /** - check that k is in valid range [0:kmax] (the test for z will be done when calling spectra_pk_at_z()) */

  //printf("k=%e\n",k);
  class_test((k < exp(pnlpt->ln_k[0])) || (k > exp(pnlpt->ln_k[pnlpt->ln_k_size-1])),
             pnlpt->error_message,
             "k=%e out of bounds [%e:%e]",k,0.,exp(pnlpt->ln_k[pnlpt->ln_k_size-1]));

  /** - compute P(k,z) (in logarithmic format for more accurate interpolation) */
  class_alloc(spectrum_at_z,
              pnlpt->ln_k_size*sizeof(double),
              pnlpt->error_message);
    
    class_alloc(spectrum_Id2d2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_Id2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_IG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_Id2G2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_IG2G2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_IFG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_IFG2_0b1_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_IFG2_0_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_IFG2_2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_CTR_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_CTR_0_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_CTR_2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_CTR_4_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_Tree_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_Tree_0_vv_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_Tree_0_vd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_Tree_0_dd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_Tree_2_vv_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_Tree_2_vd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_Tree_4_vv_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_0_vv_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_0_vd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_0_dd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_2_vv_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_2_vd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_2_dd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_4_vv_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_4_vd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_4_dd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_0_b1b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_0_b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_0_b1bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_0_bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_2_b1b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_2_b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_2_b1bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_2_bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_4_b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_4_bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_4_b1b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_4_b1bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_2_b2b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_2_b2bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_2_bG2bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    class_alloc(spectrum_4_b2b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_4_b2bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    class_alloc(spectrum_4_bG2bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);
    
    
    //GC!
    
    class_alloc(spectrum_fNL_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNLd2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNLG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

                                 
                                 
                                 //GC: ORTHOGONAL -- start


                                 class_alloc(spectrum_fNL_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNLd2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNLG2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);


                                 //GC: ORTHOGONAL -- finish

    
    
    //GC!
    
    class_alloc(spectrum_fNL_0_vv_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_0_vd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_0_dd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    
    
    class_alloc(spectrum_fNL_2_vv_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_2_vd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_2_dd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    
    
    class_alloc(spectrum_fNL_4_vv_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_4_vd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_4_dd_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    
                                 //GC: ORTHOGONAL -- start


                                 class_alloc(spectrum_fNL_0_vv_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_0_vd_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_0_dd_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 
                                 
                                 class_alloc(spectrum_fNL_2_vv_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_2_vd_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_2_dd_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 
                                 
                                 class_alloc(spectrum_fNL_4_vv_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_4_vd_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_4_dd_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);


                                 //GC: ORTHOGONAL -- finish

                                 
    
    
    //GC!
    
    
    class_alloc(spectrum_fNL_0_b1b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_0_b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_0_b1bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_0_bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    
    
    
    class_alloc(spectrum_fNL_2_b1b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_2_b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_2_b1bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_2_bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    
    
    
    class_alloc(spectrum_fNL_4_b1b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_4_b2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_4_b1bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    class_alloc(spectrum_fNL_4_bG2_at_z,
                pnlpt->ln_k_size*sizeof(double),
                pnlpt->error_message);

    
                                 //GC: ORTHOGONAL -- start


                                 class_alloc(spectrum_fNL_0_b1b2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_0_b2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_0_b1bG2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_0_bG2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 
                                 
                                 
                                 class_alloc(spectrum_fNL_2_b1b2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_2_b2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_2_b1bG2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_2_bG2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 
                                 
                                 
                                 class_alloc(spectrum_fNL_4_b1b2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_4_b2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_4_b1bG2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);

                                 class_alloc(spectrum_fNL_4_bG2_at_z_ortho,
                                             pnlpt->ln_k_size*sizeof(double),
                                             pnlpt->error_message);


                                 //GC: ORTHOGONAL -- finish

    
    

/*Correspondance between requested redshift z and number of value in z_pk input*/
int i_z;
for (i_z=0; i_z<pnlpt->z_pk_num; i_z++) {  
  if (pnlpt->z_pk[i_z] == z)
    break;
  if (i_z == pnlpt->z_pk_num-1) 
    printf("\nRequested redshift z=%e is not calculated by class. Please check input class parameter z_pk.\n",z);
}
/*printf("i_z=%d\n",i_z);*/
    
       class_call(nonlinear_pt_bias_at_z_i(
                                  pnlpt,
                                  logarithmic,
                                  i_z,
                                  spectrum_at_z,
                                  spectrum_Id2d2_at_z,
                                  spectrum_Id2_at_z,
                                       spectrum_IG2_at_z,
                                       spectrum_Id2G2_at_z,
                                       spectrum_IG2G2_at_z,
                                       spectrum_IFG2_at_z,
                                            spectrum_IFG2_0b1_at_z,
                                            spectrum_IFG2_0_at_z,
                                            spectrum_IFG2_2_at_z,
                                       spectrum_CTR_at_z,
                                            spectrum_CTR_0_at_z,
                                            spectrum_CTR_2_at_z,
                                            spectrum_CTR_4_at_z,
				                       spectrum_Tree_at_z,
                                            spectrum_Tree_0_vv_at_z,
                                            spectrum_Tree_0_vd_at_z,
                                            spectrum_Tree_0_dd_at_z,
                                            spectrum_Tree_2_vv_at_z,
                                            spectrum_Tree_2_vd_at_z,
                                            spectrum_Tree_4_vv_at_z,
                                            spectrum_0_vv_at_z,
                                            spectrum_0_vd_at_z,
                                            spectrum_0_dd_at_z,
                                            spectrum_2_vv_at_z,
                                            spectrum_2_vd_at_z,
                                            spectrum_2_dd_at_z,
                                            spectrum_4_vv_at_z,
                                            spectrum_4_vd_at_z,
                                            spectrum_4_dd_at_z,
                                            spectrum_0_b1b2_at_z,
                                            spectrum_0_b2_at_z,
                                            spectrum_0_b1bG2_at_z,
                                            spectrum_0_bG2_at_z,
                                            spectrum_2_b1b2_at_z,
                                            spectrum_2_b2_at_z,
                                            spectrum_2_b1bG2_at_z,
                                            spectrum_2_bG2_at_z,
                                            spectrum_4_b2_at_z,
                                            spectrum_4_bG2_at_z,
                                            spectrum_4_b1b2_at_z,
                                            spectrum_4_b1bG2_at_z,
                                            spectrum_2_b2b2_at_z,
                                            spectrum_2_b2bG2_at_z,
                                            spectrum_2_bG2bG2_at_z,
                                            spectrum_4_b2b2_at_z,
                                            spectrum_4_b2bG2_at_z,
                                            spectrum_4_bG2bG2_at_z, //GC!
                                            //GC!
                                            spectrum_fNL_at_z,
                                            spectrum_fNLd2_at_z,
                                            spectrum_fNLG2_at_z,
                                            //GC!
                                            spectrum_fNL_0_vv_at_z,
                                            spectrum_fNL_0_vd_at_z,
                                            spectrum_fNL_0_dd_at_z,
                                            spectrum_fNL_2_vv_at_z,
                                            spectrum_fNL_2_vd_at_z,
                                            spectrum_fNL_2_dd_at_z,
                                            spectrum_fNL_4_vv_at_z,
                                            spectrum_fNL_4_vd_at_z,
                                            spectrum_fNL_4_dd_at_z,
                                            //GC!
                                            spectrum_fNL_0_b1b2_at_z,
                                            spectrum_fNL_0_b2_at_z,
                                            spectrum_fNL_0_b1bG2_at_z,
                                            spectrum_fNL_0_bG2_at_z,
                                            spectrum_fNL_2_b1b2_at_z,
                                            spectrum_fNL_2_b2_at_z,
                                            spectrum_fNL_2_b1bG2_at_z,
                                            spectrum_fNL_2_bG2_at_z,
                                            spectrum_fNL_4_b1b2_at_z,
                                            spectrum_fNL_4_b2_at_z,
                                            spectrum_fNL_4_b1bG2_at_z,
                                            spectrum_fNL_4_bG2_at_z, //GC: ORTHOGONAL...
                                            //GC: ORTHOGONAL -- start
                                            spectrum_fNL_at_z_ortho,
                                            spectrum_fNLd2_at_z_ortho,
                                            spectrum_fNLG2_at_z_ortho,
                                            //GC!
                                            spectrum_fNL_0_vv_at_z_ortho,
                                            spectrum_fNL_0_vd_at_z_ortho,
                                            spectrum_fNL_0_dd_at_z_ortho,
                                            spectrum_fNL_2_vv_at_z_ortho,
                                            spectrum_fNL_2_vd_at_z_ortho,
                                            spectrum_fNL_2_dd_at_z_ortho,
                                            spectrum_fNL_4_vv_at_z_ortho,
                                            spectrum_fNL_4_vd_at_z_ortho,
                                            spectrum_fNL_4_dd_at_z_ortho,
                                            //GC!
                                            spectrum_fNL_0_b1b2_at_z_ortho,
                                            spectrum_fNL_0_b2_at_z_ortho,
                                            spectrum_fNL_0_b1bG2_at_z_ortho,
                                            spectrum_fNL_0_bG2_at_z_ortho,
                                            spectrum_fNL_2_b1b2_at_z_ortho,
                                            spectrum_fNL_2_b2_at_z_ortho,
                                            spectrum_fNL_2_b1bG2_at_z_ortho,
                                            spectrum_fNL_2_bG2_at_z_ortho,
                                            spectrum_fNL_4_b1b2_at_z_ortho,
                                            spectrum_fNL_4_b2_at_z_ortho,
                                            spectrum_fNL_4_b1bG2_at_z_ortho,
                                            spectrum_fNL_4_bG2_at_z_ortho
                                            //GC: ORTHOGONAL -- finish
                                            ),
               pnlpt->error_message,
               pnlpt->error_message);

    

  /** - get its second derivatives with spline, then interpolate, then convert to linear format */

  class_alloc(spline,
              sizeof(double)*pnlpt->ln_k_size,
              pnlpt->error_message);
    
    class_alloc(spline_Id2d2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_Id2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_IG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_Id2G2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_IG2G2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_IFG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_IFG2_0b1,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_IFG2_0,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_IFG2_2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_CTR,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_CTR_0,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_CTR_2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_CTR_4,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_Tree,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_Tree_0_vv,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    class_alloc(spline_Tree_0_vd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    class_alloc(spline_Tree_0_dd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    class_alloc(spline_Tree_2_vv,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    class_alloc(spline_Tree_2_vd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    class_alloc(spline_Tree_4_vv,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_0_vv,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_0_vd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_0_dd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_2_vv,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_2_vd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_2_dd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_4_vv,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_4_vd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_4_dd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    

    class_alloc(spline_0_b1b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_0_b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_0_b1bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_0_bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    
    class_alloc(spline_2_b1b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_2_b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_2_b1bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_2_bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_4_b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_4_bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_4_b1b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_4_b1bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_2_b2b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    class_alloc(spline_2_b2bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    class_alloc(spline_2_bG2bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    class_alloc(spline_4_b2b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    class_alloc(spline_4_b2bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    class_alloc(spline_4_bG2bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
    
    
    //GC!
    
    
    class_alloc(spline_fNL,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNLd2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    
    class_alloc(spline_fNLG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);
                                 
                                 
                                 //GC: ORTHOGONAL -- start



                                 class_alloc(spline_fNL_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNLd2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 
                                 class_alloc(spline_fNLG2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 
                                 //GC: ORTHOGONAL -- finish



    
    //GC!
    
    
    class_alloc(spline_fNL_0_vv,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_0_vd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_0_dd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);


    
    
    class_alloc(spline_fNL_2_vv,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_2_vd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_2_dd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);


    
    class_alloc(spline_fNL_4_vv,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_4_vd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_4_dd,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);


                                 //GC: ORTHOGONAL -- start


                                 class_alloc(spline_fNL_0_vv_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_0_vd_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_0_dd_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);


                                 
                                 
                                 class_alloc(spline_fNL_2_vv_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_2_vd_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_2_dd_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);


                                 
                                 class_alloc(spline_fNL_4_vv_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_4_vd_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_4_dd_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);


                                 //GC: ORTHOGONAL -- finish

    
    
    
    //GC!
    
    
    class_alloc(spline_fNL_0_b1b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_0_b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_0_b1bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_0_bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);


    
    
    class_alloc(spline_fNL_2_b1b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_2_b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_2_b1bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_2_bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);


    
    
    class_alloc(spline_fNL_4_b1b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_4_b2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_4_b1bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);

    class_alloc(spline_fNL_4_bG2,
                sizeof(double)*pnlpt->ln_k_size,
                pnlpt->error_message);


    
                                 //GC: ORTHOGONAL -- start


                                 class_alloc(spline_fNL_0_b1b2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_0_b2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_0_b1bG2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_0_bG2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);


                                 
                                 
                                 class_alloc(spline_fNL_2_b1b2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_2_b2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_2_b1bG2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_2_bG2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);


                                 
                                 
                                 class_alloc(spline_fNL_4_b1b2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_4_b2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_4_b1bG2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);

                                 class_alloc(spline_fNL_4_bG2_ortho,
                                             sizeof(double)*pnlpt->ln_k_size,
                                             pnlpt->error_message);


                                 //GC: ORTHOGONAL -- finish

    
    
    
    /*[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]*/
    
    
    
    
    
    
    
    
  class_call(array_spline_table_lines(pnlpt->ln_k,
                                      pnlpt->ln_k_size,
                                      spectrum_at_z,
                                      1,
                                      spline,
                                      _SPLINE_NATURAL_,
                                      pnlpt->error_message),
             pnlpt->error_message,
             pnlpt->error_message);

  class_call(array_interpolate_spline(pnlpt->ln_k,
                                      pnlpt->ln_k_size,
                                      spectrum_at_z,
                                      spline,
                                      1,
                                      log(k),
                                      &last_index,
                                      pk_tot,
                                      1,
                                      pnlpt->error_message),
             pnlpt->error_message,
             pnlpt->error_message);

  *pk_tot = exp(*pk_tot);
     free(spectrum_at_z);
      free(spline);
    
    // Id2d2
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Id2d2_at_z,
                                        1,
                                        spline_Id2d2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Id2d2_at_z,
                                        spline_Id2d2,
                                        1,
                                        log(k),
                                        &last_index2,
                                        pk_tot_Id2d2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_Id2d2 = exp(*pk_tot_Id2d2);
    
  free(spectrum_Id2d2_at_z);
  free(spline_Id2d2);
    
    // Id2
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Id2_at_z,
                                        1,
                                        spline_Id2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Id2_at_z,
                                        spline_Id2,
                                        1,
                                        log(k),
                                        &last_index3,
                                        pk_tot_Id2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_Id2 = exp(*pk_tot_Id2);
    
    free(spectrum_Id2_at_z);
    free(spline_Id2);
    
    
    // IG2
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_IG2_at_z,
                                        1,
                                        spline_IG2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_IG2_at_z,
                                        spline_IG2,
                                        1,
                                        log(k),
                                        &last_index4,
                                        pk_tot_IG2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_IG2 = exp(*pk_tot_IG2);
    
    free(spectrum_IG2_at_z);
    free(spline_IG2);
    
    
    
    // Id2G2
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Id2G2_at_z,
                                        1,
                                        spline_Id2G2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Id2G2_at_z,
                                        spline_Id2G2,
                                        1,
                                        log(k),
                                        &last_index5,
                                        pk_tot_Id2G2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_Id2G2 = exp(*pk_tot_Id2G2);
    
    free(spectrum_Id2G2_at_z);
    free(spline_Id2G2);
    
    
    // IG2G2
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_IG2G2_at_z,
                                        1,
                                        spline_IG2G2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_IG2G2_at_z,
                                        spline_IG2G2,
                                        1,
                                        log(k),
                                        &last_index6,
                                        pk_tot_IG2G2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_IG2G2 = exp(*pk_tot_IG2G2);
    
    free(spectrum_IG2G2_at_z);
    free(spline_IG2G2);
    
    
    // IFG2
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_IFG2_at_z,
                                        1,
                                        spline_IFG2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_IFG2_at_z,
                                        spline_IFG2,
                                        1,
                                        log(k),
                                        &last_index7,
                                        pk_tot_IFG2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_IFG2 = exp(*pk_tot_IFG2);
    
    free(spectrum_IFG2_at_z);
    free(spline_IFG2);
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_IFG2_0b1_at_z,
                                        1,
                                        spline_IFG2_0b1,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_IFG2_0b1_at_z,
                                        spline_IFG2_0b1,
                                        1,
                                        log(k),
                                        &last_index7,
                                        pk_tot_IFG2_0b1,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_IFG2_0b1 = exp(*pk_tot_IFG2_0b1);
    
    free(spectrum_IFG2_0b1_at_z);
    free(spline_IFG2_0b1);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_IFG2_0_at_z,
                                        1,
                                        spline_IFG2_0,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_IFG2_0_at_z,
                                        spline_IFG2_0,
                                        1,
                                        log(k),
                                        &last_index7,
                                        pk_tot_IFG2_0,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_IFG2_0 = exp(*pk_tot_IFG2_0);
    
    free(spectrum_IFG2_0_at_z);
    free(spline_IFG2_0);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_IFG2_2_at_z,
                                        1,
                                        spline_IFG2_2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_IFG2_2_at_z,
                                        spline_IFG2_2,
                                        1,
                                        log(k),
                                        &last_index7,
                                        pk_tot_IFG2_2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_IFG2_2 = exp(*pk_tot_IFG2_2);
    
    free(spectrum_IFG2_2_at_z);
    free(spline_IFG2_2);
    
    
    // CTR
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_CTR_at_z,
                                        1,
                                        spline_CTR,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_CTR_at_z,
                                        spline_CTR,
                                        1,
                                        log(k),
                                        &last_index8,
                                        pk_tot_CTR,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_CTR = exp(*pk_tot_CTR);
    
    free(spectrum_CTR_at_z);
    free(spline_CTR);
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_CTR_0_at_z,
                                        1,
                                        spline_CTR_0,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_CTR_0_at_z,
                                        spline_CTR_0,
                                        1,
                                        log(k),
                                        &last_index8,
                                        pk_tot_CTR_0,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_CTR_0 = exp(*pk_tot_CTR_0);
    
    free(spectrum_CTR_0_at_z);
    free(spline_CTR_0);
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_CTR_2_at_z,
                                        1,
                                        spline_CTR_2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_CTR_2_at_z,
                                        spline_CTR_2,
                                        1,
                                        log(k),
                                        &last_index8,
                                        pk_tot_CTR_2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_CTR_2 = exp(*pk_tot_CTR_2);
    
    free(spectrum_CTR_2_at_z);
    free(spline_CTR_2);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_CTR_4_at_z,
                                        1,
                                        spline_CTR_4,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_CTR_4_at_z,
                                        spline_CTR_4,
                                        1,
                                        log(k),
                                        &last_index8,
                                        pk_tot_CTR_4,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_CTR_4 = exp(*pk_tot_CTR_4);
    
    free(spectrum_CTR_4_at_z);
    free(spline_CTR_4);
    

   // Tree
   //
   //
   //

    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_at_z,
                                        1,
                                        spline_Tree,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_at_z,
                                        spline_Tree,
                                        1,
                                        log(k),
                                        &last_index9,
                                        pk_tot_Tree,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);


    *pk_tot_Tree = exp(*pk_tot_Tree);

    free(spectrum_Tree_at_z);
    free(spline_Tree);
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_0_vv_at_z,
                                        1,
                                        spline_Tree_0_vv,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_0_vv_at_z,
                                        spline_Tree_0_vv,
                                        1,
                                        log(k),
                                        &last_index9,
                                        pk_tot_Tree_0_vv,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_Tree_0_vv = exp(*pk_tot_Tree_0_vv);
    free(spectrum_Tree_0_vv_at_z);
    free(spline_Tree_0_vv);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_0_vd_at_z,
                                        1,
                                        spline_Tree_0_vd,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_0_vd_at_z,
                                        spline_Tree_0_vd,
                                        1,
                                        log(k),
                                        &last_index9,
                                        pk_tot_Tree_0_vd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_Tree_0_vd = exp(*pk_tot_Tree_0_vd);
    free(spectrum_Tree_0_vd_at_z);
    free(spline_Tree_0_vd);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_0_dd_at_z,
                                        1,
                                        spline_Tree_0_dd,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_0_dd_at_z,
                                        spline_Tree_0_dd,
                                        1,
                                        log(k),
                                        &last_index9,
                                        pk_tot_Tree_0_dd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    *pk_tot_Tree_0_dd = exp(*pk_tot_Tree_0_dd);
    free(spectrum_Tree_0_dd_at_z);
    free(spline_Tree_0_dd);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_2_vv_at_z,
                                        1,
                                        spline_Tree_2_vv,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_2_vv_at_z,
                                        spline_Tree_2_vv,
                                        1,
                                        log(k),
                                        &last_index9,
                                        pk_tot_Tree_2_vv,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    *pk_tot_Tree_2_vv = exp(*pk_tot_Tree_2_vv);
    free(spectrum_Tree_2_vv_at_z);
    free(spline_Tree_2_vv);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_2_vd_at_z,
                                        1,
                                        spline_Tree_2_vd,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_2_vd_at_z,
                                        spline_Tree_2_vd,
                                        1,
                                        log(k),
                                        &last_index9,
                                        pk_tot_Tree_2_vd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    *pk_tot_Tree_2_vd = exp(*pk_tot_Tree_2_vd);
    free(spectrum_Tree_2_vd_at_z);
    free(spline_Tree_2_vd);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_4_vv_at_z,
                                        1,
                                        spline_Tree_4_vv,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_Tree_4_vv_at_z,
                                        spline_Tree_4_vv,
                                        1,
                                        log(k),
                                        &last_index9,
                                        pk_tot_Tree_4_vv,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    *pk_tot_Tree_4_vv = exp(*pk_tot_Tree_4_vv);
    free(spectrum_Tree_4_vv_at_z);
    free(spline_Tree_4_vv);
    
    
    // RSD
    //
    //
    //
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_vv_at_z,
                                        1,
                                        spline_0_vv,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_vv_at_z,
                                        spline_0_vv,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_0_vv,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_0_vv = exp(*pk_tot_0_vv);
    free(spectrum_0_vv_at_z);
    free(spline_0_vv);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_vd_at_z,
                                        1,
                                        spline_0_vd,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_vd_at_z,
                                        spline_0_vd,
                                        1,
                                        log(k),
                                        &last_index11,
                                        pk_tot_0_vd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_0_vd = exp(*pk_tot_0_vd);
    free(spectrum_0_vd_at_z);
    free(spline_0_vd);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_dd_at_z,
                                        1,
                                        spline_0_dd,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_dd_at_z,
                                        spline_0_dd,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_0_dd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_0_dd = exp(*pk_tot_0_dd);
    free(spectrum_0_dd_at_z);
    free(spline_0_dd);
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_vv_at_z,
                                        1,
                                        spline_2_vv,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_vv_at_z,
                                        spline_2_vv,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_2_vv,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_2_vv = exp(*pk_tot_2_vv);
    free(spectrum_2_vv_at_z);
    free(spline_2_vv);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_vd_at_z,
                                        1,
                                        spline_2_vd,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_vd_at_z,
                                        spline_2_vd,
                                        1,
                                        log(k),
                                        &last_index11,
                                        pk_tot_2_vd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_2_vd = exp(*pk_tot_2_vd);
    free(spectrum_2_vd_at_z);
    free(spline_2_vd);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_dd_at_z,
                                        1,
                                        spline_2_dd,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_dd_at_z,
                                        spline_2_dd,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_2_dd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_2_dd = exp(*pk_tot_2_dd);
    free(spectrum_2_dd_at_z);
    free(spline_2_dd);

    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_4_vv_at_z,
                                        1,
                                        spline_4_vv,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_4_vv_at_z,
                                        spline_4_vv,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_4_vv,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_4_vv = exp(*pk_tot_4_vv);
    free(spectrum_4_vv_at_z);
    free(spline_4_vv);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_4_vd_at_z,
                                        1,
                                        spline_4_vd,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_4_vd_at_z,
                                        spline_4_vd,
                                        1,
                                        log(k),
                                        &last_index11,
                                        pk_tot_4_vd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_4_vd = exp(*pk_tot_4_vd);
    free(spectrum_4_vd_at_z);
    free(spline_4_vd);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_4_dd_at_z,
                                        1,
                                        spline_4_dd,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_4_dd_at_z,
                                        spline_4_dd,
                                        1,
                                        log(k),
                                        &last_index11,
                                        pk_tot_4_dd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_4_dd = exp(*pk_tot_4_dd);
    free(spectrum_4_dd_at_z);
    free(spline_4_dd);

    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_b1b2_at_z,
                                        1,
                                        spline_0_b1b2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_b1b2_at_z,
                                        spline_0_b1b2,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_0_b1b2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_0_b1b2 = exp(*pk_tot_0_b1b2);
    free(spectrum_0_b1b2_at_z);
    free(spline_0_b1b2);
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_b2_at_z,
                                        1,
                                        spline_0_b2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_b2_at_z,
                                        spline_0_b2,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_0_b2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_0_b2 = exp(*pk_tot_0_b2);
    free(spectrum_0_b2_at_z);
    free(spline_0_b2);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_b1bG2_at_z,
                                        1,
                                        spline_0_b1bG2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_b1bG2_at_z,
                                        spline_0_b1bG2,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_0_b1bG2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_0_b1bG2 = exp(*pk_tot_0_b1bG2);
    free(spectrum_0_b1bG2_at_z);
    free(spline_0_b1bG2);
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_bG2_at_z,
                                        1,
                                        spline_0_bG2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_0_bG2_at_z,
                                        spline_0_bG2,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_0_bG2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_0_bG2 = exp(*pk_tot_0_bG2);
    free(spectrum_0_bG2_at_z);
    free(spline_0_bG2);
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_b1b2_at_z,
                                        1,
                                        spline_2_b1b2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_b1b2_at_z,
                                        spline_2_b1b2,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_2_b1b2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_2_b1b2 = exp(*pk_tot_2_b1b2);
    free(spectrum_2_b1b2_at_z);
    free(spline_2_b1b2);
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_b2_at_z,
                                        1,
                                        spline_2_b2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_b2_at_z,
                                        spline_2_b2,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_2_b2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_2_b2 = exp(*pk_tot_2_b2);
    free(spectrum_2_b2_at_z);
    free(spline_2_b2);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_b1bG2_at_z,
                                        1,
                                        spline_2_b1bG2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_b1bG2_at_z,
                                        spline_2_b1bG2,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_2_b1bG2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_2_b1bG2 = exp(*pk_tot_2_b1bG2);
    free(spectrum_2_b1bG2_at_z);
    free(spline_2_b1bG2);
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_bG2_at_z,
                                        1,
                                        spline_2_bG2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_2_bG2_at_z,
                                        spline_2_bG2,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_2_bG2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_2_bG2 = exp(*pk_tot_2_bG2);
    free(spectrum_2_bG2_at_z);
    free(spline_2_bG2);
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_4_b2_at_z,
                                        1,
                                        spline_4_b2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_4_b2_at_z,
                                        spline_4_b2,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_4_b2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_4_b2 = exp(*pk_tot_4_b2);
    free(spectrum_4_b2_at_z);
    free(spline_4_b2);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_4_bG2_at_z,
                                        1,
                                        spline_4_bG2,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_4_bG2_at_z,
                                        spline_4_bG2,
                                        1,
                                        log(k),
                                        &last_index10,
                                        pk_tot_4_bG2,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    *pk_tot_4_bG2 = exp(*pk_tot_4_bG2);
    free(spectrum_4_bG2_at_z);
    free(spline_4_bG2);

    
    class_call(array_spline_table_lines(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_4_b1b2_at_z,1,spline_4_b1b2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_4_b1b2_at_z,spline_4_b1b2,1,log(k),&last_index10,pk_tot_4_b1b2,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    *pk_tot_4_b1b2 = exp(*pk_tot_4_b1b2);
    free(spectrum_4_b1b2_at_z);
    free(spline_4_b1b2);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_4_b1bG2_at_z,1,spline_4_b1bG2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_4_b1bG2_at_z,spline_4_b1bG2,1,log(k),&last_index10,pk_tot_4_b1bG2,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    *pk_tot_4_b1bG2 = exp(*pk_tot_4_b1bG2);
    free(spectrum_4_b1bG2_at_z);
    free(spline_4_b1bG2);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_2_b2b2_at_z,1,spline_2_b2b2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_2_b2b2_at_z,spline_2_b2b2,1,log(k),&last_index10,pk_tot_2_b2b2,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    *pk_tot_2_b2b2 = exp(*pk_tot_2_b2b2);
    free(spectrum_2_b2b2_at_z);
    free(spline_2_b2b2);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_2_b2bG2_at_z,1,spline_2_b2bG2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_2_b2bG2_at_z,spline_2_b2bG2,1,log(k),&last_index10,pk_tot_2_b2bG2,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    *pk_tot_2_b2bG2 = exp(*pk_tot_2_b2bG2);
    free(spectrum_2_b2bG2_at_z);
    free(spline_2_b2bG2);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_2_bG2bG2_at_z,1,spline_2_bG2bG2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_2_bG2bG2_at_z,spline_2_bG2bG2,1,log(k),&last_index10,pk_tot_2_bG2bG2,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    *pk_tot_2_bG2bG2 = exp(*pk_tot_2_bG2bG2);
    free(spectrum_2_bG2bG2_at_z);
    free(spline_2_bG2bG2);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_4_b2b2_at_z,1,spline_4_b2b2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_4_b2b2_at_z,spline_4_b2b2,1,log(k),&last_index10,pk_tot_4_b2b2,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    *pk_tot_4_b2b2 = exp(*pk_tot_4_b2b2);
    free(spectrum_4_b2b2_at_z);
    free(spline_4_b2b2);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_4_b2bG2_at_z,1,spline_4_b2bG2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_4_b2bG2_at_z,spline_4_b2bG2,1,log(k),&last_index10,pk_tot_4_b2bG2,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    *pk_tot_4_b2bG2 = exp(*pk_tot_4_b2bG2);
    free(spectrum_4_b2bG2_at_z);
    free(spline_4_b2bG2);
    
    class_call(array_spline_table_lines(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_4_bG2bG2_at_z,1,spline_4_bG2bG2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    class_call(array_interpolate_spline(pnlpt->ln_k,pnlpt->ln_k_size,spectrum_4_bG2bG2_at_z,spline_4_bG2bG2,1,log(k),&last_index10,pk_tot_4_bG2bG2,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    *pk_tot_4_bG2bG2 = exp(*pk_tot_4_bG2bG2);
    free(spectrum_4_bG2bG2_at_z);
    free(spline_4_bG2bG2);
    
    
    /*()()()()()()()()()()()()()()()()()()()()()()()*/
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_at_z, //spectrum_at_z_...,
                                        spline_fNL, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL = exp(*pk_tot_fNL);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL);  //free(spline_...);



    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNLd2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNLd2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNLd2_at_z, //spectrum_at_z_...,
                                        spline_fNLd2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNLd2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNLd2 = exp(*pk_tot_fNLd2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNLd2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNLd2);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNLG2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNLG2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNLG2_at_z, //spectrum_at_z_...,
                                        spline_fNLG2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNLG2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNLG2 = exp(*pk_tot_fNLG2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNLG2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNLG2);  //free(spline_...);

                                 
    
                                 //GC: ORTHOGONAL -- start


                                 class_call(array_spline_table_lines(pnlpt->ln_k,
                                                                     pnlpt->ln_k_size,
                                                                     spectrum_fNL_at_z_ortho, //spectrum_at_z_...,
                                                                     1,
                                                                     spline_fNL_ortho, //spline_...,
                                                                     _SPLINE_NATURAL_,
                                                                     pnlpt->error_message),
                                            pnlpt->error_message,
                                            pnlpt->error_message);

                                 class_call(array_interpolate_spline(pnlpt->ln_k,
                                                                     pnlpt->ln_k_size,
                                                                     spectrum_fNL_at_z_ortho, //spectrum_at_z_...,
                                                                     spline_fNL_ortho, //spline_...,
                                                                     1,
                                                                     log(k),
                                                                     &last_index_fNL,
                                                                     pk_tot_fNL_ortho, //pk_tot_...,
                                                                     1,
                                                                     pnlpt->error_message),
                                            pnlpt->error_message,
                                            pnlpt->error_message);

                                 *pk_tot_fNL_ortho = exp(*pk_tot_fNL_ortho);//*pk_tot_... = exp(*pk_tot_...);
                                 free(spectrum_fNL_at_z_ortho);  //free(spectrum_at_z_...);
                                 free(spline_fNL_ortho);  //free(spline_...);

                                 
                                 
                                 class_call(array_spline_table_lines(pnlpt->ln_k,
                                                                     pnlpt->ln_k_size,
                                                                     spectrum_fNLd2_at_z_ortho, //spectrum_at_z_...,
                                                                     1,
                                                                     spline_fNLd2_ortho, //spline_...,
                                                                     _SPLINE_NATURAL_,
                                                                     pnlpt->error_message),
                                            pnlpt->error_message,
                                            pnlpt->error_message);

                                 class_call(array_interpolate_spline(pnlpt->ln_k,
                                                                     pnlpt->ln_k_size,
                                                                     spectrum_fNLd2_at_z_ortho, //spectrum_at_z_...,
                                                                     spline_fNLd2_ortho, //spline_...,
                                                                     1,
                                                                     log(k),
                                                                     &last_index_fNL,
                                                                     pk_tot_fNLd2_ortho, //pk_tot_...,
                                                                     1,
                                                                     pnlpt->error_message),
                                            pnlpt->error_message,
                                            pnlpt->error_message);

                                 *pk_tot_fNLd2_ortho = exp(*pk_tot_fNLd2_ortho);//*pk_tot_... = exp(*pk_tot_...);
                                 free(spectrum_fNLd2_at_z_ortho);  //free(spectrum_at_z_...);
                                 free(spline_fNLd2_ortho);  //free(spline_...);

                                 
                                 
                                 class_call(array_spline_table_lines(pnlpt->ln_k,
                                                                     pnlpt->ln_k_size,
                                                                     spectrum_fNLG2_at_z_ortho, //spectrum_at_z_...,
                                                                     1,
                                                                     spline_fNLG2_ortho, //spline_...,
                                                                     _SPLINE_NATURAL_,
                                                                     pnlpt->error_message),
                                            pnlpt->error_message,
                                            pnlpt->error_message);

                                 class_call(array_interpolate_spline(pnlpt->ln_k,
                                                                     pnlpt->ln_k_size,
                                                                     spectrum_fNLG2_at_z_ortho, //spectrum_at_z_...,
                                                                     spline_fNLG2_ortho, //spline_...,
                                                                     1,
                                                                     log(k),
                                                                     &last_index_fNL,
                                                                     pk_tot_fNLG2_ortho, //pk_tot_...,
                                                                     1,
                                                                     pnlpt->error_message),
                                            pnlpt->error_message,
                                            pnlpt->error_message);

                                 *pk_tot_fNLG2_ortho = exp(*pk_tot_fNLG2_ortho);//*pk_tot_... = exp(*pk_tot_...);
                                 free(spectrum_fNLG2_at_z_ortho);  //free(spectrum_at_z_...);
                                 free(spline_fNLG2_ortho);  //free(spline_...);


                                 //GC: ORTHOGONAL -- finish

                                 
                                 
    
            /*==============================================================*/
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_vv_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_vv, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_vv_at_z, //spectrum_at_z_...,
                                        spline_fNL_0_vv, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_vv, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_vv = exp(*pk_tot_fNL_0_vv);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_vv_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_0_vv);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_vd_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_vd, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_vd_at_z, //spectrum_at_z_...,
                                        spline_fNL_0_vd, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_vd, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_vd = exp(*pk_tot_fNL_0_vd);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_vd_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_0_vd);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_dd_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_dd, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_dd_at_z, //spectrum_at_z_...,
                                        spline_fNL_0_dd, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_dd, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_dd = exp(*pk_tot_fNL_0_dd);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_dd_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_0_dd);  //free(spline_...);

    
    

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_vv_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_vv, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_vv_at_z, //spectrum_at_z_...,
                                        spline_fNL_2_vv, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_vv, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_vv = exp(*pk_tot_fNL_2_vv);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_vv_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_2_vv);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_vd_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_vd, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_vd_at_z, //spectrum_at_z_...,
                                        spline_fNL_2_vd, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_vd, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_vd = exp(*pk_tot_fNL_2_vd);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_vd_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_2_vd);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_dd_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_dd, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_dd_at_z, //spectrum_at_z_...,
                                        spline_fNL_2_dd, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_dd, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_dd = exp(*pk_tot_fNL_2_dd);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_dd_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_2_dd);  //free(spline_...);

    
    

    
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_vv_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_vv, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_vv_at_z, //spectrum_at_z_...,
                                        spline_fNL_4_vv, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_vv, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_vv = exp(*pk_tot_fNL_4_vv);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_vv_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_4_vv);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_vd_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_vd, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_vd_at_z, //spectrum_at_z_...,
                                        spline_fNL_4_vd, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_vd, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_vd = exp(*pk_tot_fNL_4_vd);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_vd_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_4_vd);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_dd_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_dd, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_dd_at_z, //spectrum_at_z_...,
                                        spline_fNL_4_dd, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_dd, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_dd = exp(*pk_tot_fNL_4_dd);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_dd_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_4_dd);  //free(spline_...);

                                 
                                 
                                 
                                 
                                 
      //GC: ORTHOGONAL -- start


                                 
                                 
       class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_vv_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_vv_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_vv_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_0_vv_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_vv_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_vv_ortho = exp(*pk_tot_fNL_0_vv_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_vv_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_0_vv_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_vd_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_vd_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_vd_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_0_vd_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_vd_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_vd_ortho = exp(*pk_tot_fNL_0_vd_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_vd_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_0_vd_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_dd_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_dd_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_dd_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_0_dd_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_dd_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_dd_ortho = exp(*pk_tot_fNL_0_dd_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_dd_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_0_dd_ortho);  //free(spline_...);

    
    

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_vv_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_vv_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_vv_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_2_vv_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_vv_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_vv_ortho = exp(*pk_tot_fNL_2_vv_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_vv_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_2_vv_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_vd_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_vd_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_vd_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_2_vd_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_vd_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_vd_ortho = exp(*pk_tot_fNL_2_vd_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_vd_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_2_vd_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_dd_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_dd_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_dd_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_2_dd_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_dd_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_dd_ortho = exp(*pk_tot_fNL_2_dd_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_dd_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_2_dd_ortho);  //free(spline_...);

    
    

    
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_vv_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_vv_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_vv_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_4_vv_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_vv_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_vv_ortho = exp(*pk_tot_fNL_4_vv_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_vv_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_4_vv_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_vd_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_vd_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_vd_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_4_vd_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_vd_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_vd_ortho = exp(*pk_tot_fNL_4_vd_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_vd_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_4_vd_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_dd_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_dd_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_dd_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_4_dd_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_dd_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_dd_ortho = exp(*pk_tot_fNL_4_dd_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_dd_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_4_dd_ortho);  //free(spline_...);
                              
                                 
                                 
                                 
                                 
                                 
                                 


    //GC: ORTHOGONAL -- finish
                               
                                 
                                 
                                 
                                 
                                 
                                 
                                 
    
    
                /*===========================================================*/
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_b1b2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_b1b2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_b1b2_at_z, //spectrum_at_z_...,
                                        spline_fNL_0_b1b2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_b1b2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_b1b2 = exp(*pk_tot_fNL_0_b1b2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_b1b2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_0_b1b2);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_b2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_b2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_b2_at_z, //spectrum_at_z_...,
                                        spline_fNL_0_b2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_b2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_b2 = exp(*pk_tot_fNL_0_b2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_b2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_0_b2);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_b1bG2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_b1bG2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_b1bG2_at_z, //spectrum_at_z_...,
                                        spline_fNL_0_b1bG2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_b1bG2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_b1bG2 = exp(*pk_tot_fNL_0_b1bG2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_b1bG2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_0_b1bG2);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_bG2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_bG2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_bG2_at_z, //spectrum_at_z_...,
                                        spline_fNL_0_bG2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_bG2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_bG2 = exp(*pk_tot_fNL_0_bG2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_bG2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_0_bG2);  //free(spline_...);

    
    

    
    
    
    
    
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_b1b2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_b1b2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_b1b2_at_z, //spectrum_at_z_...,
                                        spline_fNL_2_b1b2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_b1b2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_b1b2 = exp(*pk_tot_fNL_2_b1b2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_b1b2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_2_b1b2);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_b2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_b2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_b2_at_z, //spectrum_at_z_...,
                                        spline_fNL_2_b2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_b2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_b2 = exp(*pk_tot_fNL_2_b2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_b2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_2_b2);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_b1bG2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_b1bG2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_b1bG2_at_z, //spectrum_at_z_...,
                                        spline_fNL_2_b1bG2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_b1bG2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_b1bG2 = exp(*pk_tot_fNL_2_b1bG2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_b1bG2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_2_b1bG2);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_bG2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_bG2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_bG2_at_z, //spectrum_at_z_...,
                                        spline_fNL_2_bG2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_bG2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_bG2 = exp(*pk_tot_fNL_2_bG2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_bG2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_2_bG2);  //free(spline_...);

    
    

    
    
    
    
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_b1b2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_b1b2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_b1b2_at_z, //spectrum_at_z_...,
                                        spline_fNL_4_b1b2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_b1b2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_b1b2 = exp(*pk_tot_fNL_4_b1b2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_b1b2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_4_b1b2);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_b2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_b2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_b2_at_z, //spectrum_at_z_...,
                                        spline_fNL_4_b2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_b2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_b2 = exp(*pk_tot_fNL_4_b2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_b2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_4_b2);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_b1bG2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_b1bG2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_b1bG2_at_z, //spectrum_at_z_...,
                                        spline_fNL_4_b1bG2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_b1bG2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_b1bG2 = exp(*pk_tot_fNL_4_b1bG2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_b1bG2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_4_b1bG2);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_bG2_at_z, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_bG2, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_bG2_at_z, //spectrum_at_z_...,
                                        spline_fNL_4_bG2, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_bG2, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_bG2 = exp(*pk_tot_fNL_4_bG2);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_bG2_at_z);  //free(spectrum_at_z_...);
    free(spline_fNL_4_bG2);  //free(spline_...);

                                 
                                 
                                 
                                 
                                 
                                 //GC: ORTHOGONAL -- start

                                 
                                 
                                 
     class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_b1b2_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_b1b2_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_b1b2_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_0_b1b2_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_b1b2_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_b1b2_ortho = exp(*pk_tot_fNL_0_b1b2_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_b1b2_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_0_b1b2_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_b2_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_b2_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_b2_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_0_b2_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_b2_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_b2_ortho = exp(*pk_tot_fNL_0_b2_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_b2_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_0_b2_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_b1bG2_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_b1bG2_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_b1bG2_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_0_b1bG2_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_b1bG2_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_b1bG2_ortho = exp(*pk_tot_fNL_0_b1bG2_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_b1bG2_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_0_b1bG2_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_bG2_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_0_bG2_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_0_bG2_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_0_bG2_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_0_bG2_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_0_bG2_ortho = exp(*pk_tot_fNL_0_bG2_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_0_bG2_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_0_bG2_ortho);  //free(spline_...);

    
    

    
    
    
    
    
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_b1b2_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_b1b2_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_b1b2_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_2_b1b2_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_b1b2_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_b1b2_ortho = exp(*pk_tot_fNL_2_b1b2_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_b1b2_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_2_b1b2_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_b2_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_b2_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_b2_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_2_b2_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_b2_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_b2_ortho = exp(*pk_tot_fNL_2_b2_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_b2_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_2_b2_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_b1bG2_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_b1bG2_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_b1bG2_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_2_b1bG2_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_b1bG2_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_b1bG2_ortho = exp(*pk_tot_fNL_2_b1bG2_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_b1bG2_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_2_b1bG2_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_bG2_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_2_bG2_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_2_bG2_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_2_bG2_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_2_bG2_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_2_bG2_ortho = exp(*pk_tot_fNL_2_bG2_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_2_bG2_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_2_bG2_ortho);  //free(spline_...);

    
    

    
    
    
    
    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_b1b2_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_b1b2_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_b1b2_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_4_b1b2_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_b1b2_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_b1b2_ortho = exp(*pk_tot_fNL_4_b1b2_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_b1b2_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_4_b1b2_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_b2_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_b2_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_b2_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_4_b2_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_b2_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_b2_ortho = exp(*pk_tot_fNL_4_b2_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_b2_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_4_b2_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_b1bG2_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_b1bG2_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_b1bG2_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_4_b1bG2_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_b1bG2_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_b1bG2_ortho = exp(*pk_tot_fNL_4_b1bG2_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_b1bG2_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_4_b1bG2_ortho);  //free(spline_...);

    
    
    class_call(array_spline_table_lines(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_bG2_at_z_ortho, //spectrum_at_z_...,
                                        1,
                                        spline_fNL_4_bG2_ortho, //spline_...,
                                        _SPLINE_NATURAL_,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_call(array_interpolate_spline(pnlpt->ln_k,
                                        pnlpt->ln_k_size,
                                        spectrum_fNL_4_bG2_at_z_ortho, //spectrum_at_z_...,
                                        spline_fNL_4_bG2_ortho, //spline_...,
                                        1,
                                        log(k),
                                        &last_index_fNL,
                                        pk_tot_fNL_4_bG2_ortho, //pk_tot_...,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    *pk_tot_fNL_4_bG2_ortho = exp(*pk_tot_fNL_4_bG2_ortho);//*pk_tot_... = exp(*pk_tot_...);
    free(spectrum_fNL_4_bG2_at_z_ortho);  //free(spectrum_at_z_...);
    free(spline_fNL_4_bG2_ortho);  //free(spline_...);
                                
                                 
                                 


                                 //GC: ORTHOGONAL -- finish

        
    

  return _SUCCESS_;

}


/**
 * Vectorized interpolation of PT power spectra at an array of k values.
 *
 * This function is equivalent to calling nonlinear_pt_pk_at_k_and_z() for
 * each k value, but avoids the per-k overhead of allocating/freeing ~192
 * arrays and recomputing 96 splines. Instead, it:
 * - Calls nonlinear_pt_bias_at_z_i() once to get all 96 spectra
 * - Computes spline coefficients once per component
 * - Evaluates the spline at all k values
 *
 * Output layout: pk_mult[c * n_k + ik] = exp(spline(ln(kvec[ik]))) for
 * component c and k-index ik. The 96 components are in the same order as
 * the output pointers of nonlinear_pt_pk_at_k_and_z().
 */
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

    /* Find z-index */
    for (i_z = 0; i_z < pnlpt->z_pk_num; i_z++) {
        if (pnlpt->z_pk[i_z] == z) break;
    }
    if (i_z == pnlpt->z_pk_num) {
        sprintf(pnlpt->error_message,
                "nonlinear_pt_pk_mult_at_kvec_and_z: z=%e not in z_pk list", z);
        return _FAILURE_;
    }

    /* Allocate 96 spectrum arrays for bias_at_z_i */
    double * spectra[96];
    for (c = 0; c < 96; c++) {
        class_alloc(spectra[c], nk * sizeof(double), pnlpt->error_message);
    }

    /* Copy data from struct arrays */
    class_call(nonlinear_pt_bias_at_z_i(
        pnlpt, logarithmic, i_z,
        spectra[0], spectra[1], spectra[2], spectra[3], spectra[4],
        spectra[5], spectra[6], spectra[7], spectra[8], spectra[9],
        spectra[10], spectra[11], spectra[12], spectra[13], spectra[14],
        spectra[15], spectra[16], spectra[17], spectra[18], spectra[19],
        spectra[20], spectra[21], spectra[22], spectra[23], spectra[24],
        spectra[25], spectra[26], spectra[27], spectra[28], spectra[29],
        spectra[30], spectra[31], spectra[32], spectra[33], spectra[34],
        spectra[35], spectra[36], spectra[37], spectra[38], spectra[39],
        spectra[40], spectra[41], spectra[42], spectra[43], spectra[44],
        spectra[45], spectra[46], spectra[47],
        spectra[48], spectra[49], spectra[50],
        spectra[51], spectra[52], spectra[53], spectra[54], spectra[55],
        spectra[56], spectra[57], spectra[58], spectra[59],
        spectra[60], spectra[61], spectra[62], spectra[63],
        spectra[64], spectra[65], spectra[66], spectra[67],
        spectra[68], spectra[69], spectra[70], spectra[71],
        spectra[72], spectra[73], spectra[74],
        spectra[75], spectra[76], spectra[77], spectra[78], spectra[79],
        spectra[80], spectra[81], spectra[82], spectra[83],
        spectra[84], spectra[85], spectra[86], spectra[87],
        spectra[88], spectra[89], spectra[90], spectra[91],
        spectra[92], spectra[93], spectra[94], spectra[95]
    ), pnlpt->error_message, pnlpt->error_message);

    /* Pre-compute ln(k) array */
    class_alloc(ln_kvec, n_k * sizeof(double), pnlpt->error_message);
    for (ik = 0; ik < n_k; ik++)
        ln_kvec[ik] = log(kvec[ik]);

    /* Allocate single spline buffer (reused for each component) */
    class_alloc(spline_c, nk * sizeof(double), pnlpt->error_message);

    /* For each component: compute spline once, evaluate at all k */
    for (c = 0; c < 96; c++) {
        class_call(array_spline_table_lines(
            pnlpt->ln_k, nk, spectra[c], 1, spline_c,
            _SPLINE_NATURAL_, pnlpt->error_message),
            pnlpt->error_message, pnlpt->error_message);

        last_index = 0;
        for (ik = 0; ik < n_k; ik++) {
            class_call(array_interpolate_spline(
                pnlpt->ln_k, nk, spectra[c], spline_c, 1,
                ln_kvec[ik], &last_index, &val, 1,
                pnlpt->error_message),
                pnlpt->error_message, pnlpt->error_message);
            pk_mult[c * n_k + ik] = exp(val);
        }

        free(spectra[c]);
    }

    free(spline_c);
    free(ln_kvec);

    return _SUCCESS_;
}

