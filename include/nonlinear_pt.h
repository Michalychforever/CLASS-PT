/** @file nonlinear.h Documented includes for trg module */

#include "primordial.h"

#ifndef __NONLINEAR_PT__
#define __NONLINEAR_PT__

#define _M_EV_TOO_BIG_FOR_HALOFIT_ 10. /**< above which value of non-CDM mass (in eV) do we stop trusting halofit? */

/**
 * Maximum number of values of redshift at which the spectra will be
 * written in output files
 */

#define _Z_PK_NUM_MAX_ 100

enum non_linear_method_pt {nlpt_none,nlpt_spt};
enum non_linear_irresumm_pt {irres_yes,irres_no};
enum non_linear_bias_pt {bias_yes,bias_no};

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
  //@}

  /** @name - table non-linear corrections for matter density, sqrt(P_NL(k,z)/P_NL(k,z)) */

  //@{

    // M.I. Global collable variables have to be defined here (they are called with prefix pnlpt->)
    
    int z_pk_num;
    double z_pk[_Z_PK_NUM_MAX_];
    
  int k_size;      /**< k_size = total number of k values */
  double * k;      /**< k[index_k] = list of k values */
  int tau_size;    /**< tau_size = number of values */
  double * tau;    /**< tau[index_tau] = list of time values */

//  double * pk_nl_out;
    
    
      double * M13_oneline;
      double * M22_oneline;
    double * M22basic_oneline;
    double * IFG2_oneline;
    
    
  double * nl_corr_density;/**< nl_corr_density[index_tau * ppt->k_size + index_k] */
  double * nl_corr_Id2d2;
    double * nl_corr_IG2;
    double * nl_corr_Id2;
    double * nl_corr_Id2G2;
    double * nl_corr_IG2G2;
        double * nl_corr_IFG2;
    double * nl_corr_CTR;
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

  int nonlinear_pt_free(
                     struct nonlinear_pt *pnlpt
                     );

  int nonlinear_pt_pk_l(struct perturbs *ppt,
                     struct primordial *ppm,
                     struct nonlinear_pt *pnlpt,
                     int index_tau,
                     double *pk_l,
                     double *lnk,
                     double *lnpk,
                     double *ddlnpk);
    
    int nonlinear_pt_loop(
                          struct precision *ppr,
                          struct background *pba,
                          struct primordial *ppm,
                          struct nonlinear_pt *pnlpt,
                          double tau,
                          double *pk_l,
                          double *pk_nl,
                          double *pk_Id2d2,
                          double *pk_Id2,
                          double *pk_IG2,
                          double *pk_Id2G2,
                          double *pk_IG2G2,
                          double *pk_IFG2,
                          double *pk_CTR,
                          double *lnk_l,
                          double *lnpk_l,
                          double *ddlnpk_l
                          );
    

    
    int nonlinear_pt_p13UV_integrate(
                                    struct nonlinear_pt *pnlpt,
                                    struct background *pba,
                                    double * integrand_array,
                                    int integrand_size,
                                    int ia_size,
                                    int index_ia_k,
                                    int index_ia_pk,
                                    int index_ia_sum,
                                    int index_ia_ddsum,
                                    double * sum
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
