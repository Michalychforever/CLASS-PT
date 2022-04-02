/** @file nonlinear.h Documented includes for trg module */
#include <unistd.h>
#include <errno.h>
#include "primordial.h"

#ifndef __NONLINEAR_PT__
#define __NONLINEAR_PT__

/**
 * Maximum number of values of redshift at which the spectra will be
 * written in output files
 */

#define _Z_PK_NUM_MAX_ 100


enum non_linear_method_pt {nlpt_none,nlpt_spt};
enum non_linear_irresumm_pt {irres_yes,irres_no};
enum non_linear_bias_pt {bias_yes,bias_no};
enum non_linear_rsd_pt {rsd_yes,rsd_no};
// enum non_linear_cb_pt {cb_yes,cb_no};
// enum non_linear_fftlogn {fftlog_fast,fftlog_norm};
enum non_linear_AP_effect_pt {AP_effect_yes,AP_effect_no};
enum non_linear_fNL_equil_ortho {fNL_equil_ortho_yes,fNL_equil_ortho_no};

//enum halofit_integral_type {halofit_integral_one, halofit_integral_two, halofit_integral_three};
//enum halofit_statement {ok, too_small};

/**
 * Structure containing all information on non-linear spectra.
 *
 * Once initialized by nonlinear_init(), contains a table for all two points correlation functions
 * and for all the ai,bj functions (containing the three points correlation functions), for each
 * time and wave-number.
 */

struct nonlinear_pt {

  /** @name - input parameters initialized by user in input module
      (all other quantities are computed in this module, given these
      parameters and the content of the 'precision', 'background',
      'thermo', 'primordial' and 'spectra' structures) */

  //@{

enum non_linear_method_pt method; /**< method for computing non-linear corrections (none, Halogit, etc.) */
enum non_linear_irresumm_pt irres;
enum non_linear_bias_pt bias;
enum non_linear_rsd_pt rsd;
// enum non_linear_cb_pt cb;
// enum non_linear_fftlogn norm;
enum non_linear_AP_effect_pt AP_effect;
enum non_linear_fNL_equil_ortho fNL_equil_ortho_switch;

  //@}

  /** @name - table non-linear corrections for matter density, sqrt(P_NL(k,z)/P_NL(k,z)) */

  //@{

    // M.I. Global collable variables have to be defined here (they are called with prefix pnlpt->)

    int z_pk_num;
    double z_pk[_Z_PK_NUM_MAX_];

  int k_size;      /**< k_size = total number of k values */
  int ln_k_size;   /**< k_size = total number of k values */
  double * k;      /**< k[index_k] = list of k values */
  double * ln_k;      /**< k[index_k] = list of k values */
  int tau_size;    /**< tau_size = number of values */
  double * tau;    /**< tau[index_tau] = list of time values */
  double * ln_tau;    /**< ln_tau is absolutely different from psp->ln_tau!! It is a simple log of tau with size 609*/

  int index_md_scalars,ic_size,tp_size;
  double * dd_sources_tp_delta_m;
  double * dd_sources_tp_delta_cb;
  double * sources_tp_delta_m;
  double * sources_tp_delta_cb;

short fast_output;
short cb;

//  double * pk_nl_out;

  //@{

  int k_size_cmb;  /**< k_size_cmb[index_md] number of k values used
                        for CMB calculations, requiring a fine
                        sampling in k-space */

  int k_size_cl;  /**< k_size_cl[index_md] number of k values used
                       for non-CMB \f$ C_l \f$ calculations, requiring a coarse
                       sampling in k-space. */

  double k_min;     /**< minimum value (over all modes) */
  double k_max;     /**< maximum value (over all modes) */

  //@}

double OmfidAP;

double * M13_oneline;
double * M22_oneline;
double * M22basic_oneline;
double * IFG2_oneline;
//double * M12_oneline;    //GC!


double * M12_oneline;

    //GC: ORTHOGONAL -- start

double * M12_oneline_ortho;
    
    //GC: ORTHOGONAL -- finish
    
/**/
double * M12_oneline_matter_multipoles_vv0_f2;
double * M12_oneline_matter_multipoles_vv0_f3;
double * M12_oneline_matter_multipoles_vd0_f1;
double * M12_oneline_matter_multipoles_vd0_f2;
double * M12_oneline_matter_multipoles_dd0_f0;
double * M12_oneline_matter_multipoles_dd0_f1;
double * M12_oneline_matter_multipoles_vv2_f3;
double * M12_oneline_matter_multipoles_vd2_f2;
double * M12_oneline_matter_multipoles_vv4_f3;
double * M12_oneline_matter_multipoles_vd4_f2;
/**/
    
    //GC: ORTHOGONAL -- start

    double * M12_oneline_matter_multipoles_vv0_f2_ortho;
    double * M12_oneline_matter_multipoles_vv0_f3_ortho;
    double * M12_oneline_matter_multipoles_vd0_f1_ortho;
    double * M12_oneline_matter_multipoles_vd0_f2_ortho;
    double * M12_oneline_matter_multipoles_dd0_f0_ortho;
    double * M12_oneline_matter_multipoles_dd0_f1_ortho;
    double * M12_oneline_matter_multipoles_vv2_f3_ortho;
    double * M12_oneline_matter_multipoles_vd2_f2_ortho;
    double * M12_oneline_matter_multipoles_vv4_f3_ortho;
    double * M12_oneline_matter_multipoles_vd4_f2_ortho;
    
    //GC: ORTHOGONAL -- finish

 
double * M12_oneline_matter_mu_powers_vd2_f1;
double * M12_oneline_matter_mu_powers_vd2_f2;
double * M12_oneline_matter_mu_powers_dd2_f1;
double * M12_oneline_matter_mu_powers_vv4_f2;
double * M12_oneline_matter_mu_powers_vd4_f2;
double * M12_oneline_matter_mu_powers_vv6_f3;

    
    //GC: ORTHOGONAL -- start

double * M12_oneline_matter_mu_powers_vd2_f1_ortho;
double * M12_oneline_matter_mu_powers_vd2_f2_ortho;
double * M12_oneline_matter_mu_powers_dd2_f1_ortho;
double * M12_oneline_matter_mu_powers_vv4_f2_ortho;
double * M12_oneline_matter_mu_powers_vd4_f2_ortho;
double * M12_oneline_matter_mu_powers_vv6_f3_ortho;

    
    //GC: ORTHOGONAL -- finish

    
/**/
double * M12_oneline_bias_real_space_b2;
double * M12_oneline_bias_real_space_bG2;
/**/

    
    //GC: ORTHOGONAL -- start

    
    double * M12_oneline_bias_real_space_b2_ortho;
    double * M12_oneline_bias_real_space_bG2_ortho;

    
    //GC: ORTHOGONAL -- finish

 
double * M12_oneline_bias_multipoles_b2_vv0_f1;
double * M12_oneline_bias_multipoles_bG2_vv0_f1;
    
    //GC: ORTHOGONAL -- finish

    double * M12_oneline_bias_multipoles_b2_vv0_f1_ortho;
    double * M12_oneline_bias_multipoles_bG2_vv0_f1_ortho;
        
    //GC: ORTHOGONAL -- finish

    
    //GC --> go now to the complexified ones...
    
    //GC --> THESE ALL REMAIN...

double complex * M13_oneline_complex;
double complex * M22_oneline_complex;
double complex * M22basic_oneline_complex;
double complex * IFG2_oneline_complex;

double complex * M22_oneline_0_vv_complex;
double complex * M13_0_vv_oneline_complex;

double complex * M22_oneline_0_vd_complex;
double complex * M13_0_vd_oneline_complex;

double complex * M22_oneline_0_dd_complex;
double complex * M13_0_dd_oneline_complex;

    double complex * M22_oneline_2_vv_complex;
    double complex * M13_2_vv_oneline_complex;

    double complex * M22_oneline_4_vv_complex;
    double complex * M13_4_vv_oneline_complex;

    double complex * M22_oneline_2_vd_complex;
    double complex * M13_2_vd_oneline_complex;

    double complex * M22_oneline_4_vd_complex;
    double complex * M13_4_vd_oneline_complex;
    double complex * M22_oneline_4_dd_complex;

    double complex * M22_oneline_2_dd_complex;
    double complex * M13_2_dd_oneline_complex;

double complex * M22_0_b1b2_oneline_complex;
double complex * M22_0_b2_oneline_complex;
double complex * M22_0_b1bG2_oneline_complex;
    double complex * M22_0_bG2_oneline_complex;


    double complex * M22_2_b1b2_oneline_complex;
    double complex * M22_2_b2_oneline_complex;
    double complex * M22_2_b1bG2_oneline_complex;
    double complex * M22_2_bG2_oneline_complex;

    double complex * M22_4_b2_oneline_complex;
    double complex * M22_4_bG2_oneline_complex;

double complex * M_Id2;
double complex * M_IG2;
double complex * M_Id2G2;
double complex * M_IG2G2;

double complex * M22_oneline_mu2_vd_complex;
double complex * M22_oneline_mu2_dd_complex;
    double complex * M22_oneline_mu4_vv_complex;
    double complex * M22_oneline_mu4_vd_complex;

    double complex * M22_oneline_mu4_dd_complex;
    double complex * M22_oneline_mu6_vv_complex;
    double complex * M22_oneline_mu6_vd_complex;
    double complex * M22_oneline_mu8_complex;

double complex * M13_mu2_dd_oneline_complex;
    double complex * M13_mu2_vd_oneline_complex;
    double complex * M13_mu4_vv_oneline_complex;
    double complex * M13_mu4_vd_oneline_complex;
    double complex * M13_mu6_oneline_complex;
    
    
    //GC!
    
    
    double complex * M12_oneline_complex;
    
    //GC: ORTHOGONAL -- start

    double complex * M12_oneline_complex_ortho;

    //GC: ORTHOGONAL -- finish



    double complex * M12_oneline_complex_matter_multipoles_vv0_f2;
    double complex * M12_oneline_complex_matter_multipoles_vv0_f3;
    double complex * M12_oneline_complex_matter_multipoles_vd0_f1;
    double complex * M12_oneline_complex_matter_multipoles_vd0_f2;
    double complex * M12_oneline_complex_matter_multipoles_dd0_f0;
    double complex * M12_oneline_complex_matter_multipoles_dd0_f1;
    double complex * M12_oneline_complex_matter_multipoles_vv2_f3;
    double complex * M12_oneline_complex_matter_multipoles_vd2_f2;
    double complex * M12_oneline_complex_matter_multipoles_vv4_f3;
    double complex * M12_oneline_complex_matter_multipoles_vd4_f2;
    
    
    //GC: ORTHOGONAL -- start
    
    
    double complex * M12_oneline_complex_matter_multipoles_vv0_f2_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vv0_f3_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vd0_f1_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vd0_f2_ortho;
    double complex * M12_oneline_complex_matter_multipoles_dd0_f0_ortho;
    double complex * M12_oneline_complex_matter_multipoles_dd0_f1_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vv2_f3_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vd2_f2_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vv4_f3_ortho;
    double complex * M12_oneline_complex_matter_multipoles_vd4_f2_ortho;


    //GC: ORTHOGONAL -- finish



    double complex * M12_oneline_complex_matter_mu_powers_vd2_f1;
    double complex * M12_oneline_complex_matter_mu_powers_vd2_f2;
    double complex * M12_oneline_complex_matter_mu_powers_dd2_f1;
    double complex * M12_oneline_complex_matter_mu_powers_vv4_f2;
    double complex * M12_oneline_complex_matter_mu_powers_vd4_f2;
    double complex * M12_oneline_complex_matter_mu_powers_vv6_f3;
    
    
    //GC: ORTHOGONAL -- start
    
    
    double complex * M12_oneline_complex_matter_mu_powers_vd2_f1_ortho;
    double complex * M12_oneline_complex_matter_mu_powers_vd2_f2_ortho;
    double complex * M12_oneline_complex_matter_mu_powers_dd2_f1_ortho;
    double complex * M12_oneline_complex_matter_mu_powers_vv4_f2_ortho;
    double complex * M12_oneline_complex_matter_mu_powers_vd4_f2_ortho;
    double complex * M12_oneline_complex_matter_mu_powers_vv6_f3_ortho;


    //GC: ORTHOGONAL -- finish



    double complex * M12_oneline_complex_bias_real_space_b2;
    double complex * M12_oneline_complex_bias_real_space_bG2;
    
    
    //GC: ORTHOGONAL -- start

    double complex * M12_oneline_complex_bias_real_space_b2_ortho;
    double complex * M12_oneline_complex_bias_real_space_bG2_ortho;

    //GC: ORTHOGONAL -- finish



    double complex * M12_oneline_complex_bias_multipoles_b2_vv0_f1;
    double complex * M12_oneline_complex_bias_multipoles_bG2_vv0_f1;
    
    
    //GC: ORTHOGONAL -- start

    double complex * M12_oneline_complex_bias_multipoles_b2_vv0_f1_ortho;
    double complex * M12_oneline_complex_bias_multipoles_bG2_vv0_f1_ortho;
    

    //GC: ORTHOGONAL -- finish



    double complex * M12_oneline_0_vv_complex;
    double complex * M12_oneline_0_vd_complex;
    double complex * M12_oneline_0_dd_complex;
    double complex * M12_oneline_2_vv_complex;
    double complex * M12_oneline_2_vd_complex;
    double complex * M12_oneline_2_dd_complex;
    double complex * M12_oneline_4_vv_complex;
    double complex * M12_oneline_4_vd_complex;
    double complex * M12_oneline_mu2_vd_complex;
    double complex * M12_oneline_mu2_dd_complex;
    double complex * M12_oneline_mu4_vv_complex;
    double complex * M12_oneline_mu4_vd_complex;
    double complex * M12_oneline_mu6_vv_complex;
    
    
    //GC: ORTHOGONAL -- start

    double complex * M12_oneline_0_vv_complex_ortho;
    double complex * M12_oneline_0_vd_complex_ortho;
    double complex * M12_oneline_0_dd_complex_ortho;
    double complex * M12_oneline_2_vv_complex_ortho;
    double complex * M12_oneline_2_vd_complex_ortho;
    double complex * M12_oneline_2_dd_complex_ortho;
    double complex * M12_oneline_4_vv_complex_ortho;
    double complex * M12_oneline_4_vd_complex_ortho;
    double complex * M12_oneline_mu2_vd_complex_ortho;
    double complex * M12_oneline_mu2_dd_complex_ortho;
    double complex * M12_oneline_mu4_vv_complex_ortho;
    double complex * M12_oneline_mu4_vd_complex_ortho;
    double complex * M12_oneline_mu6_vv_complex_ortho;

    //GC: ORTHOGONAL -- finish

    
    
    double complex * M12_2_bG2_oneline_complex;
    double complex * M12_2_b2_oneline_complex;
    double complex * M12_0_bG2_oneline_complex;
    double complex * M12_0_b1bG2_oneline_complex;
    double complex * M12_0_b2_oneline_complex;
    double complex * M12_0_b1b2_oneline_complex;
    
    
    //GC: ORTHOGONAL -- start

    double complex * M12_2_bG2_oneline_complex_ortho;
    double complex * M12_2_b2_oneline_complex_ortho;
    double complex * M12_0_bG2_oneline_complex_ortho;
    double complex * M12_0_b1bG2_oneline_complex_ortho;
    double complex * M12_0_b2_oneline_complex_ortho;
    double complex * M12_0_b1b2_oneline_complex_ortho;
    

    //GC: ORTHOGONAL -- finish


    double complex * M_fNLd2;
    double complex * M_fNLG2;
    
    
    //GC: ORTHOGONAL -- start

    double complex * M_fNLd2_ortho;
    double complex * M_fNLG2_ortho;

    //GC: ORTHOGONAL -- finish
    
    
    //GC!


double * ln_pk_nl;     /*For classy[i_z*pnlpt->k_size+index_k]*/
double * ln_pk_Id2d2;
    double * ln_pk_Id2d2_2;
    double * ln_pk_Id2d2_4;
double * ln_pk_Id2;
double * ln_pk_IG2;
double * ln_pk_Id2G2;
    double * ln_pk_Id2G2_2;
    double * ln_pk_Id2G2_4;
double * ln_pk_IG2G2;
    double * ln_pk_IG2G2_2;
    double * ln_pk_IG2G2_4;
double * ln_pk_IFG2;
    double * ln_pk_IFG2_0b1;
    double * ln_pk_IFG2_0;
    double * ln_pk_IFG2_2;
double * ln_pk_CTR;
    double * ln_pk_CTR_0;
    double * ln_pk_CTR_2;
    double * ln_pk_CTR_4;
double * ln_pk_Tree;
    double * ln_pk_Tree_0_vv;
    double * ln_pk_Tree_0_vd;
    double * ln_pk_Tree_0_dd;
    double * ln_pk_Tree_2_vv;
    double * ln_pk_Tree_2_vd;
    double * ln_pk_Tree_4_vv;

double * ln_pk_0_vv;
double * ln_pk_0_vd;
double * ln_pk_0_dd;

    double * ln_pk_2_vv;
    double * ln_pk_2_vd;
    double * ln_pk_2_dd;

    double * ln_pk_4_vv;
    double * ln_pk_4_vd;

    double * ln_pk_4_dd;

double * ln_pk_0_b1b2;
double * ln_pk_0_b1bG2;
double * ln_pk_0_b2;
    double * ln_pk_0_bG2;

    double * ln_pk_2_b1b2;
    double * ln_pk_2_b1bG2;
    double * ln_pk_2_b2;
    double * ln_pk_2_bG2;

    double * ln_pk_4_b2;
    double * ln_pk_4_bG2;
    double * ln_pk_4_b1b2;
    double * ln_pk_4_b1bG2;
    
    //GC!
        
    double * ln_pk_fNL_0_vv;
    double * ln_pk_fNL_0_vd;
    double * ln_pk_fNL_0_dd;

    double * ln_pk_fNL_2_vv;
    double * ln_pk_fNL_2_vd;
    double * ln_pk_fNL_2_dd;

    double * ln_pk_fNL_4_vv;
    double * ln_pk_fNL_4_vd;
    double * ln_pk_fNL_4_dd;
    
    
    //GC: ORTHOGONAL -- start


    double * ln_pk_fNL_0_vv_ortho;
    double * ln_pk_fNL_0_vd_ortho;
    double * ln_pk_fNL_0_dd_ortho;

    double * ln_pk_fNL_2_vv_ortho;
    double * ln_pk_fNL_2_vd_ortho;
    double * ln_pk_fNL_2_dd_ortho;

    double * ln_pk_fNL_4_vv_ortho;
    double * ln_pk_fNL_4_vd_ortho;
    double * ln_pk_fNL_4_dd_ortho;

    
    //GC: ORTHOGONAL -- finish


    
    double * ln_pk12_0_b1b2;
    double * ln_pk12_0_b2;
    double * ln_pk12_0_b1bG2;
    double * ln_pk12_0_bG2;
    double * ln_pk12_2_b1b2;
    double * ln_pk12_2_b2;
    double * ln_pk12_2_b1bG2;
    double * ln_pk12_2_bG2;
    double * ln_pk12_4_b1b2;
    double * ln_pk12_4_b2;
    double * ln_pk12_4_b1bG2;
    double * ln_pk12_4_bG2;
    
    
    
    //GC: ORTHOGONAL -- start
    
    
    double * ln_pk12_0_b1b2_ortho;
    double * ln_pk12_0_b2_ortho;
    double * ln_pk12_0_b1bG2_ortho;
    double * ln_pk12_0_bG2_ortho;
    double * ln_pk12_2_b1b2_ortho;
    double * ln_pk12_2_b2_ortho;
    double * ln_pk12_2_b1bG2_ortho;
    double * ln_pk12_2_bG2_ortho;
    double * ln_pk12_4_b1b2_ortho;
    double * ln_pk12_4_b2_ortho;
    double * ln_pk12_4_b1bG2_ortho;
    double * ln_pk12_4_bG2_ortho;


    //GC: ORTHOGONAL -- finish


    
    double * ln_pk_nl_fNL;
    double * ln_pk_fNLd2;
    double * ln_pk_fNLG2;
    
    
    //GC: ORTHOGONAL -- start
    
    
    double * ln_pk_nl_fNL_ortho;
    double * ln_pk_fNLd2_ortho;
    double * ln_pk_fNLG2_ortho;


    //GC: ORTHOGONAL -- finish


    //GC!

    double * growthf;   /*For RSD effect[i_z]*/
    double * hratio_array;   /*For RSD effect[i_z]*/
    double * Dratio_array;   /*For RSD effect[i_z]*/

    double * gauss;
    double * gauss_x;
    double * gauss_w;

    // double *koff;
    // double *Poff;
    // double *offtab;

/*< nl_corr_density[index_tau * ppt->k_size + index_k] */
    double * nl_corr_density;
    /*
    double * nl_corr_density;
    double * nl_corr_Id2d2;
    double * nl_corr_IG2;
    double * nl_corr_Id2;
    double * nl_corr_Id2G2;
    double * nl_corr_IG2G2;
    double * nl_corr_IFG2;
    double * nl_corr_CTR;
    double * nl_corr_Tree;
    */


    char input_pk[500];
    int replace_pk;
    int replace_background;
    int no_wiggle;
    int wiggle_only;
    double alpha_rs;
    double replace_Hz_value;
    double replace_DAz_value;
    double replace_Dz_value;
    double replace_fz_value;

//  double * k_nl;  /**< wavenumber at which non-linear corrections become important, defined differently by different non_linear_method's */
  int index_tau_min_nl; /**< index of smallest value of tau at which nonlinear corrections have been computed (so, for tau<tau_min_nl, the array nl_corr_density only contains some factors 1 */

  //@}

  /** @name - technical parameters */

  //@{

  short nonlinear_pt_verbose;  	/**< amount of information written in standard output */

  ErrorMsg error_message; 	/**< zone for writing error messages */

  //@}
};

/********************************************************************************/


extern void zspmv_(char*, int*, double complex*, double complex*, double complex*, int*, double complex*, double complex*, int*);
extern double complex zdotu_(int*, double complex*, int*, double complex*, int*);

/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif


  int nonlinear_pt_init(
                     struct precision *ppr,
                     struct background *pba,
                     struct thermo *pth,
                     struct perturbs *ppt,
                     struct primordial *ppm,
                     struct nonlinear_pt *pnlpt
                     );

  int perturb_get_k_list_nl(
                         struct precision * ppr,
                         struct background * pba,
                         struct thermo * pth,
                         struct perturbs * ppt,
                         struct nonlinear_pt *pnlpt
                         );

  int nonlinear_pt_free(
                     struct nonlinear_pt *pnlpt
                     );

  int nonlinear_pt_pk_l(
                     struct background *pba,
                     struct perturbs *ppt,
                     struct primordial *ppm,
                     struct nonlinear_pt *pnlpt,
                     int index_tau,
                     double *pk_l,
                     double *lnk,
                     double *lnpk,
                     double *ddlnpk); //,
                     //double *tk_l,
                     //double *lntk,
                     //double *ddlntk); //GC

//GC!

//GC: the function that extracts the primordial PS to get the transfer function {\cal M}(k)...

    int nonlinear_pt_pPRIMk_l(
                              //struct background *pba,
                              struct perturbs *ppt,
                              struct primordial *ppm,
                              struct nonlinear_pt *pnlpt,
                              double *lnk,
                              double *pPRIMk_l,
                              double *lnpPRIMk,
                              double *ddlnpPRIMk);

//GC!


    int nonlinear_pt_loop(
                          struct precision *ppr,
                          struct background *pba,
                          struct primordial *ppm,
                          struct thermo *pth,
                          struct nonlinear_pt *pnlpt,
                          double tau,
                          double f,
                          double hratio,
                          double Dratio,
                          double *pk_l,
                          double *pPRIMk_l, //GC! Never used...
                          double *pk_l_0_vv,
                          double *pk_l_0_vd,
                          double *pk_l_0_dd,
                          double *pk_l_2_vv,
                          double *pk_l_2_vd,
                          double *pk_l_2_dd,
                          double *pk_l_4_vv,
                          double *pk_l_4_vd,
                          double *pk_l_4_dd,
                          double *pk_nl,
                          double *pk_nl_fNL, //GC!
                          double *pk_fNLd2, //GC!
                          double *pk_fNLG2, //GC!
                          //GC: ORTHOGONAL -- start
                          double *pk_nl_fNL_ortho, //GC!
                          double *pk_fNLd2_ortho, //GC!
                          double *pk_fNLG2_ortho, //GC!
                          //GC: ORTHOGONAL -- finish
                          double *pk_Id2d2,
                          double *pk_Id2d2_2,
                          double *pk_Id2d2_4,
                          double *pk_l_0_b1b2,
                          double *pk_l_0_b2,
                          double *pk_l_0_b1bG2,
                          double *pk_l_0_bG2,
                          double *pk_l_2_b1b2,
                          double *pk_l_2_b2,
                          double *pk_l_2_b1bG2,
                          double *pk_l_2_bG2,
                          double *pk_l_4_b2,
                          double *pk_l_4_bG2,
                          double *pk_l_4_b1b2,
                          double *pk_l_4_b1bG2,
                          double *pk_Id2,
                          double *pk_IG2,
                          double *pk_Id2G2,
                          double *pk_Id2G2_2,
                          double *pk_Id2G2_4,
                          double *pk_IG2G2,
                          double *pk_IG2G2_2,
                          double *pk_IG2G2_4,
                          double *pk_IFG2,
                          double *pk_IFG2_0b1,
                          double *pk_IFG2_0,
                          double *pk_IFG2_2,
                          double *pk_CTR,
                          double *pk_CTR_0,
                          double *pk_CTR_2,
                          double *pk_CTR_4,
                          double *pk_Tree,
                          double *pk_Tree_0_vv,
                          double *pk_Tree_0_vd,
                          double *pk_Tree_0_dd,
                          double *pk_Tree_2_vv,
                          double *pk_Tree_2_vd,
                          double *pk_Tree_4_vv,
                          //GC!
                          double *pk_l_fNL_0_vv,
                          double *pk_l_fNL_0_vd,
                          double *pk_l_fNL_0_dd,
                          double *pk_l_fNL_2_vv,
                          double *pk_l_fNL_2_vd,
                          double *pk_l_fNL_2_dd,
                          double *pk_l_fNL_4_vv,
                          double *pk_l_fNL_4_vd,
                          double *pk_l_fNL_4_dd,
                          double *pk12_l_0_b1b2,
                          double *pk12_l_0_b2,
                          double *pk12_l_0_b1bG2,
                          double *pk12_l_0_bG2,
                          double *pk12_l_2_b1b2,
                          double *pk12_l_2_b2,
                          double *pk12_l_2_b1bG2,
                          double *pk12_l_2_bG2,
                          double *pk12_l_4_b1b2,
                          double *pk12_l_4_b2,
                          double *pk12_l_4_b1bG2,
                          double *pk12_l_4_bG2,
                          //GC: ORTHOGONAL -- start
                          double *pk_l_fNL_0_vv_ortho,
                          double *pk_l_fNL_0_vd_ortho,
                          double *pk_l_fNL_0_dd_ortho,
                          double *pk_l_fNL_2_vv_ortho,
                          double *pk_l_fNL_2_vd_ortho,
                          double *pk_l_fNL_2_dd_ortho,
                          double *pk_l_fNL_4_vv_ortho,
                          double *pk_l_fNL_4_vd_ortho,
                          double *pk_l_fNL_4_dd_ortho,
                          double *pk12_l_0_b1b2_ortho,
                          double *pk12_l_0_b2_ortho,
                          double *pk12_l_0_b1bG2_ortho,
                          double *pk12_l_0_bG2_ortho,
                          double *pk12_l_2_b1b2_ortho,
                          double *pk12_l_2_b2_ortho,
                          double *pk12_l_2_b1bG2_ortho,
                          double *pk12_l_2_bG2_ortho,
                          double *pk12_l_4_b1b2_ortho,
                          double *pk12_l_4_b2_ortho,
                          double *pk12_l_4_b1bG2_ortho,
                          double *pk12_l_4_bG2_ortho,
                          //GC: ORTHOGONAL -- finish
                          //GC!
                          double *lnk_l,
                          double *lnpk_l,
                          double *lnpPRIMk_l //GC!
                          );


    /**
     * Function definitions of the FFT tool, used by transfer.c and spectra.c
     * For more information, see fft.c
     */

#ifndef FFT_DEFINED
#define FFT_DEFINED
    /**
     * Compute FFT of two real inputs and store into two complex outputs of sizes N
     * Assumes arrays are allocated and of size N
     */
    void FFT_real(double* input_real_1, double* input_real_2, double* output_real_1, double* output_imag_1, double* output_real_2, double* output_imag_2, int N);
    /**
     * Compute FFT of two real inputs and store into two complex outputs of sizes N
     * Assumes arrays are allocated and of size N
     * Only gives up nonredundant part for real FFT which have c_(N-n)=c_n
     */
    void FFT_real_short(double* input_real_1, double* input_real_2, double* output_real_1, double* output_imag_1, double* output_real_2, double* output_imag_2, int N);
    /**
     * Compute FFT of single complex input and stores into single complex output of size N
     * Assumes arrays are allocated and of size N
     */
    void FFT(double* input_real, double* input_imag, double* output_real, double* output_imag, int N, int stepsize);
    void DCT_real(double* input_real,double* input_imag,double* output_real,double* output_imag,int N);
#endif


#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
/* @endcond */
