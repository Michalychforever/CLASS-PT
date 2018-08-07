/** @file nonlinear.c Documented nonlinear module
 *
 * Julien Lesgourgues, 6.03.2014
 *
 * New module replacing an older one present up to version 2.0 The new
 * module is located in a better place in the main, allowing it to
 * compute non-linear correction to \f$ C_l\f$'s and not just \f$ P(k)\f$. It will
 * also be easier to generalize to new methods.  The old implementation
 * of one-loop calculations and TRG calculations has been dropped from
 * this version, they can still be found in older versions.
 *
 */
#include "nonlinear.h"

/*#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#undef I
#include <fftw3.h>
 */

int nonlinear_k_nl_at_z(
                        struct background *pba,
                        struct nonlinear * pnl,
                        double z,
                        double * k_nl
                        ) {

  double tau;

  class_call(background_tau_of_z(pba,
                                 z,
                                 &tau),
             pba->error_message,
             pnl->error_message);

  if (pnl->tau_size == 1) {
    *k_nl = pnl->k_nl[0];
  }
  else {
    class_call(array_interpolate_two(pnl->tau,
                                     1,
                                     0,
                                     pnl->k_nl,
                                     1,
                                     pnl->tau_size,
                                     tau,
                                     k_nl,
                                     1,
                                     pnl->error_message),
               pnl->error_message,
               pnl->error_message);
  }

  return _SUCCESS_;
}

int nonlinear_init(
                   struct precision *ppr,
                   struct background *pba,
                   struct thermo *pth,
                   struct perturbs *ppt,
                   struct primordial *ppm,
                   struct nonlinear *pnl
                   ) {

  int index_ncdm;
  int index_k;
  int index_tau;
  double *pk_l;
  double *pk_nl;
  double *lnk_l;
  double *lnpk_l;
  double *ddlnpk_l;
  short print_warning=_FALSE_;
  double * pvecback;
  int last_index;
  double a,z;
  enum halofit_statement halofit_found_k_max;

  /** Summary
   *
   * (a) First deal with the case where non non-linear corrections requested */

  if (pnl->method == nl_none) {
    if (pnl->nonlinear_verbose > 0)
      printf("No non-linear spectra requested. Nonlinear module skipped.\n");
  }

  /** (b) Compute for HALOFIT non-linear spectrum */

  else if (pnl->method == nl_halofit) {
    if (pnl->nonlinear_verbose > 0)
      printf("Computing non-linear matter power spectrum with Halofit (including update Takahashi et al. 2012 and Bird 2014)\n");

    if (pba->has_ncdm) {
      for (index_ncdm=0;index_ncdm < pba->N_ncdm; index_ncdm++){
        if (pba->m_ncdm_in_eV[index_ncdm] >  _M_EV_TOO_BIG_FOR_HALOFIT_)
          fprintf(stdout,"Warning: Halofit is proved to work for CDM, and also with a small HDM component thanks to Bird et al.'s update. But it sounds like you are running with a WDM component of mass %f eV, which makes the use of Halofit suspicious.\n",pba->m_ncdm_in_eV[index_ncdm]);
      }
    }
      
    /** - copy list of (k,tau) from perturbation module */

    pnl->k_size = ppt->k_size[ppt->index_md_scalars];
    class_alloc(pnl->k,pnl->k_size*sizeof(double),pnl->error_message);
    for (index_k=0; index_k<pnl->k_size; index_k++)
      pnl->k[index_k] = ppt->k[ppt->index_md_scalars][index_k];

    pnl->tau_size = ppt->tau_size;
    class_alloc(pnl->tau,pnl->tau_size*sizeof(double),pnl->error_message);
    for (index_tau=0; index_tau<pnl->tau_size; index_tau++)
      pnl->tau[index_tau] = ppt->tau_sampling[index_tau];

    class_alloc(pnl->nl_corr_density,pnl->tau_size*pnl->k_size*sizeof(double),pnl->error_message);
    class_alloc(pnl->k_nl,pnl->tau_size*sizeof(double),pnl->error_message);

    class_alloc(pk_l,pnl->k_size*sizeof(double),pnl->error_message);
    class_alloc(pk_nl,pnl->k_size*sizeof(double),pnl->error_message);

    class_alloc(lnk_l,pnl->k_size*sizeof(double),pnl->error_message);
    class_alloc(lnpk_l,pnl->k_size*sizeof(double),pnl->error_message);
    class_alloc(ddlnpk_l,pnl->k_size*sizeof(double),pnl->error_message);

    /** - loop over time */

      //M: I'm switching this loop off!
      //M: uncomment:
      /*for (index_tau = pnl->tau_size-1; index_tau>=0; index_tau--) {*/
      
      for (index_tau = pnl->tau_size-1; index_tau > pnl->tau_size-2; index_tau--) {

      /* get P_L(k) at this time */
      class_call(nonlinear_pk_l(ppt,ppm,pnl,index_tau,pk_l,lnk_l,lnpk_l,ddlnpk_l),
                 pnl->error_message,
                 pnl->error_message);

       /* get P_NL(k) at this time */
      if (print_warning == _FALSE_) {
          class_call(nonlinear_pt(ppr,
                                     pba,
                                     ppm,
                                     pnl,
                                     pnl->tau[index_tau],
                                     pk_l,
                                     pk_nl,
                                     lnk_l,
                                     lnpk_l,
                                     ddlnpk_l
//                                     &halofit_found_k_max
                                  ),
                   pnl->error_message,
                   pnl->error_message);
            
            for (index_k=0; index_k<pnl->k_size; index_k++) {
                pnl->nl_corr_density[index_tau * pnl->k_size + index_k] = sqrt(pk_nl[index_k]/pk_l[index_k]);
            }


      }
      
      else {
        for (index_k=0; index_k<pnl->k_size; index_k++) {
          pnl->nl_corr_density[index_tau * pnl->k_size + index_k] = 1.;
        }
      }
    }
      
    
    free(pk_l);
    free(pk_nl);
    free(lnk_l);
    free(lnpk_l);
    free(ddlnpk_l);
  }

  else {
    class_stop(pnl->error_message,
               "Your non-linear method variable is set to %d, out of the range defined in nonlinear.h",pnl->method);
  }

  return _SUCCESS_;
    
} //end of non-linear init

int nonlinear_free(
                   struct nonlinear *pnl
                   ) {

  if (pnl->method > nl_none) {

    if (pnl->method == nl_halofit) {
      /* free here */
      free(pnl->k);
      free(pnl->tau);
      free(pnl->nl_corr_density);
      free(pnl->k_nl);
    }
  }

  return _SUCCESS_;

}// end of nonlinear_free

// This function executes the linear power spectrum
int nonlinear_pk_l(
                   struct perturbs *ppt,
                   struct primordial *ppm,
                   struct nonlinear *pnl,
                   int index_tau,
                   double *pk_l,
                   double *lnk,
                   double *lnpk,
                   double *ddlnpk) {

  int index_md;
  int index_k;
  int index_ic1,index_ic2,index_ic1_ic2;
  double * primordial_pk;
  double source_ic1,source_ic2;

  index_md = ppt->index_md_scalars;

  class_alloc(primordial_pk,ppm->ic_ic_size[index_md]*sizeof(double),pnl->error_message);

  for (index_k=0; index_k<pnl->k_size; index_k++) {

    class_call(primordial_spectrum_at_k(ppm,
                                        index_md,
                                        linear,
                                        pnl->k[index_k],
                                        primordial_pk),
               ppm->error_message,
               pnl->error_message);

    pk_l[index_k] = 0;

    /* part diagonal in initial conditions */
    for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {

      index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_md]);

      source_ic1 = ppt->sources[index_md]
        [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
        [index_tau * ppt->k_size[index_md] + index_k];
        
       // source_ic1 are transfer functions
        
      pk_l[index_k] += 2.*_PI_*_PI_/pow(pnl->k[index_k],3)
        *source_ic1*source_ic1
        *primordial_pk[index_ic1_ic2];
    }

      
    /* part non-diagonal in initial conditions */
    for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
      for (index_ic2 = index_ic1+1; index_ic2 < ppm->ic_size[index_md]; index_ic2++) {

        index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic2,ppm->ic_size[index_md]);

        if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {

          source_ic1 = ppt->sources[index_md]
            [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
            [index_tau * ppt->k_size[index_md] + index_k];

          source_ic2 = ppt->sources[index_md]
            [index_ic2 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
            [index_tau * ppt->k_size[index_md] + index_k];

          pk_l[index_k] += 2.*2.*_PI_*_PI_/pow(pnl->k[index_k],3)
            *source_ic1*source_ic2
            *primordial_pk[index_ic1_ic2]; // extra 2 factor (to include the symmetric term ic2,ic1)

        }
      }
    }

    lnk[index_k] = log(pnl->k[index_k]);
    lnpk[index_k] = log(pk_l[index_k]);
    //  printf("%e  %e\n",pnl->k[index_k],pk_l[index_k]);
  }

  class_call(array_spline_table_columns(lnk,
                                        pnl->k_size,
                                        lnpk,
                                        1,
                                        ddlnpk,
                                        _SPLINE_NATURAL_,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  free(primordial_pk);

  return _SUCCESS_;

} // end of non-linear_Pk_l. This function executes the linear power spectrum

//This function executes the halofit
int nonlinear_halofit(
                      struct precision *ppr,
                      struct background *pba,
                      struct primordial *ppm,
                      struct nonlinear *pnl,
                      double tau,
                      double *pk_l,
                      double *pk_nl,
                      double *lnk_l,
                      double *lnpk_l,
                      double *ddlnpk_l,
                      double *k_nl,
                      enum halofit_statement * halofit_found_k_max
                      ) {

  double Omega_m,Omega_v,fnu,Omega0_m, w0, dw_over_da_fld, integral_fld;

  /** Determine non linear ratios (from pk) **/

  int index_k;
  double pk_lin,pk_quasi,pk_halo,rk;
  double sigma,rknl,rneff,rncur,d1,d2,sigmav;
  double diff,xlogr1,xlogr2,rmid;

  double gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3;
  double pk_linaa;
  double y;
  double f1a,f2a,f3a,f1b,f2b,f3b,frac;

  double * pvecback;

  int last_index;
  int counter;
  double sum1,sum2,sum3,sum0,sum01;
  double anorm;

  double *integrand_array;
  int integrand_size;
  int index_ia_k;
  int index_ia_pk;
  int index_ia_sum;
  int index_ia_ddsum;
  /*
  int index_ia_sum2;
  int index_ia_ddsum2;
  int index_ia_sum3;
  int index_ia_ddsum3;
  */
  int ia_size;
  int index_ia;

  double k_integrand;
  double lnpk_integrand;

  double R;

  //M: this array has the size of the background vector
  class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);

  Omega0_m = (pba->Omega0_cdm + pba->Omega0_b + pba->Omega0_ncdm_tot + pba->Omega0_dcdm);

  /* Halofit needs w0 = w_fld today */
  class_call(background_w_fld(pba,pba->a_today,&w0,&dw_over_da_fld,&integral_fld), pba->error_message, pnl->error_message);

  fnu      = pba->Omega0_ncdm_tot/Omega0_m;
  anorm    = 1./(2*pow(_PI_,2));
    
     //fprintf(stderr, "Misha is here.\n");
    //M: here there's a small loop

  /*      Until the 17.02.2015 the values of k used for integrating sigma(R) quantities needed by Halofit where the same as in the perturbation module.
          Since then, we sample these integrals on more values, in order to get more precise integrals (thanks Matteo Zennaro for noticing the need for this).

     We create a temporary integrand_array which columns will be:
     - k in 1/Mpc
     - just linear P(k) in Mpc**3
     - 1/(2(pi**2)) P(k) k**2 exp(-(kR)**2) or 1/(2(pi**2)) P(k) k**2 2 (kR) exp(-(kR)**2) or 1/(2(pi**2)) P(k) k**2 4 (kR)(1-kR) exp(-(kR)**2)
     - second derivative of previous line with spline
  */

  index_ia=0;
  class_define_index(index_ia_k,     _TRUE_,index_ia,1);
  class_define_index(index_ia_pk,    _TRUE_,index_ia,1);
  class_define_index(index_ia_sum,   _TRUE_,index_ia,1);
  class_define_index(index_ia_ddsum, _TRUE_,index_ia,1);
  ia_size = index_ia;

   //M: here the integrand_array is tabulated. First one computes the size of the integrand (number of points I presume)
    
  integrand_size=(int)(log(pnl->k[pnl->k_size-1]/pnl->k[0])/log(10.)*ppr->halofit_k_per_decade)+1;

    //printf("Misha is here.\n", integrand_size);
   // printf("%d\n", integrand_size);
    
  class_alloc(integrand_array,integrand_size*ia_size*sizeof(double),pnl->error_message);

  //fprintf(stderr,"Omega_m=%e,  fnu=%e\n",Omega0_m,fnu);

  /* we fill integrand_array with values of k and P(k) using interpolation */

  last_index=0;

  for (index_k=0; index_k < integrand_size; index_k++) {

    k_integrand=pnl->k[0]*pow(10.,index_k/ppr->halofit_k_per_decade);

    class_call(array_interpolate_spline(lnk_l,
                                        pnl->k_size,
                                        lnpk_l,           /*what is being interpolated*/
                                        ddlnpk_l,
                                        1,
                                        log(k_integrand),
                                        &last_index,
                                        &lnpk_integrand,
                                        1,
                                        pnl->error_message),
               pnl->error_message,
               pnl->error_message);

    integrand_array[index_k*ia_size + index_ia_k] = k_integrand;
   //   printf("%d\n", k_integrand);
    integrand_array[index_k*ia_size + index_ia_pk] = exp(lnpk_integrand);

  }
    //printf("k_max=%e\n", integrand_array[(integrand_size-1)*ia_size + index_ia_k]);
    //k_max=5.717018
 /*
  ia_size=4
  index_ia_k=0
  index_ia_pk=1
  integrand_size=473
  */
/*
printf("ia_size=%i\n", ia_size);
printf("index_ia_k=%i\n", index_ia_k);
printf("index_ia_pk=%i\n", index_ia_pk);
printf("integrand_size=%i\n", integrand_size);
 */
    
// class_call only calls a function!
    
  class_call(background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback),
             pba->error_message,
             pnl->error_message);

  Omega_m = pvecback[pba->index_bg_Omega_m];
  Omega_v = 1.-pvecback[pba->index_bg_Omega_m]-pvecback[pba->index_bg_Omega_r];

  // for debugging:
  printf("Call Halofit at z=%e\n",pba->a_today/pvecback[pba->index_bg_a]-1.);
  // M: Here's a loop over redshifts !
    
  /* minimum value of R such that the integral giving sigma_R is
     converged.  The parameter halofit_sigma_precision should be
     understood as follows: we trust our calculation of sigma(R) as
     long as the integral reaches a value k_max such that the factor
     exp(-(Rk_max)**2) is already as low as halofit_sigma_precisio,
     shoing that the integreal is converged.  In practise this
     condition is tested only for R_max, the highest value of R in our
     bisection algorithm. Hence a smaller value of
     halofit_sigma_precision will lead to a more precise halofit
     result at the *highest* redshift at which halofit can make
     computations, at the expense of requiring a larger k_max; but
     this parameter is not relevant for the precision on P_nl(k,z) at
     other redshifts, so there is normally no need to change i
   */

  R=sqrt(-log(ppr->halofit_sigma_precision))/integrand_array[(integrand_size-1)*ia_size + index_ia_k];

  class_call(nonlinear_halofit_integrate(
                                         pnl,
                                         integrand_array,
                                         integrand_size,
                                         ia_size,
                                         index_ia_k,
                                         index_ia_pk,
                                         index_ia_sum,
                                         index_ia_ddsum,
                                         R,
                                         halofit_integral_one,
                                         &sum1
                                         ),
             pnl->error_message,
             pnl->error_message);

  sigma  = sqrt(sum1);

  /* the following error should not stop the code: it will arrive
     inevitably at some large redshift, and then the code should not
     stop, but just give up computing P_NL(k,z). This is why we have a
     special error handling here (using class_test_except and free()
     commands to avoid memory leaks, and calling this whole function
     not through a class_call) */

  /*
  class_test_except(sigma < 1.,
                    pnl->error_message,
                    free(pvecback);free(integrand_array),
                    "Your k_max=%g 1/Mpc is too small for Halofit to find the non-linearity scale z_nl at z=%g. Increase input parameter P_k_max_h/Mpc or P_k_max_1/Mpc",
                    pnl->k[pnl->k_size-1],
                    pba->a_today/pvecback[pba->index_bg_a]-1.);
  */

    
  if (sigma < 1.) {
    * halofit_found_k_max = too_small;
    free(pvecback);
    free(integrand_array);
    return _SUCCESS_;
  }
  else {
    * halofit_found_k_max = ok;
  }
     

  xlogr1 = log(R)/log(10.);

  /* maximum value of R in the bisection algorithm leading to the
     determination of R_nl.  For this value we can make a
     conservaitive guess: 1/halofit_min_k_nonlinear, where
     halofit_min_k_nonlinear is the minimum value of k at which we ask
     halofit to give us an estimate of P_nl(k,z). By assumption we
     treat all smaller k's as linear, so we know that
     sigma(1/halofit_min_k_nonlinear) must be <<1 (and if it is not
     the test below will alert us) */

  R=1./ppr->halofit_min_k_nonlinear;

  /* corresponding value of sigma_R */
  class_call(nonlinear_halofit_integrate(
                                         pnl,
                                         integrand_array,
                                         integrand_size,
                                         ia_size,
                                         index_ia_k,
                                         index_ia_pk,
                                         index_ia_sum,
                                         index_ia_ddsum,
                                         R,
                                         halofit_integral_one,
                                         &sum1
                                         ),
             pnl->error_message,
             pnl->error_message);

  sigma  = sqrt(sum1);
    
  class_test(sigma > 1.,
             pnl->error_message,
             "Your input value for the precision parameter halofit_min_k_nonlinear=%e is too large, such that sigma(R=1/halofit_min_k_nonlinear)=% > 1. For self-consistency, it should have been <1. Decrease halofit_min_k_nonlinear",
             ppr->halofit_min_k_nonlinear,sigma);

  xlogr2 = log(R)/log(10.);

  counter = 0;
  do {
    rmid = pow(10,(xlogr2+xlogr1)/2.0);
    counter ++;

    class_call(nonlinear_halofit_integrate(
                                           pnl,
                                           integrand_array,
                                           integrand_size,
                                           ia_size,
                                           index_ia_k,
                                           index_ia_pk,
                                           index_ia_sum,
                                           index_ia_ddsum,
                                           rmid,
                                           halofit_integral_one,
                                           &sum1
                                           ),
             pnl->error_message,
             pnl->error_message);

    sigma  = sqrt(sum1);

    diff = sigma - 1.0;

    if (diff > ppr->halofit_tol_sigma){
      xlogr1=log10(rmid);
    }
    else if (diff < -ppr->halofit_tol_sigma) {
      xlogr2 = log10(rmid);
    }

    /* The first version of this test woukld let the code continue: */
    /*
    class_test_except(counter > _MAX_IT_,
                      pnl->error_message,
                      free(pvecback);free(integrand_array),
                      "could not converge within maximum allowed number of iterations");
    */
    /* ... but in this situation it sounds better to make it stop and return an error! */
    class_test(counter > _MAX_IT_,
               pnl->error_message,
               "could not converge within maximum allowed number of iterations");

  } while (fabs(diff) > ppr->halofit_tol_sigma);

  /* evaluate all the other integrals at R=rmid */

  class_call(nonlinear_halofit_integrate(
                                         pnl,
                                         integrand_array,
                                         integrand_size,
                                         ia_size,
                                         index_ia_k,
                                         index_ia_pk,
                                         index_ia_sum,
                                         index_ia_ddsum,
                                         rmid,
                                         halofit_integral_two,
                                         &sum2
                                         ),
             pnl->error_message,
             pnl->error_message);

  class_call(nonlinear_halofit_integrate(
                                         pnl,
                                         integrand_array,
                                         integrand_size,
                                         ia_size,
                                         index_ia_k,
                                         index_ia_pk,
                                         index_ia_sum,
                                         index_ia_ddsum,
                                         rmid,
                                         halofit_integral_three,
                                         &sum3
                                         ),
             pnl->error_message,
             pnl->error_message);
    
    
  sigma  = sqrt(sum1);
  d1 = -sum2/sum1;
  d2 = -sum2*sum2/sum1/sum1 - sum3/sum1;

  rknl  = 1./rmid;
  rneff = -3.-d1;
  rncur = -d2;

  *k_nl = rknl;

  //fprintf(stderr,"Here\n");
    
// here I'm still in the loop
    
    
// loop over k here - Multiply matrices here then
    
  for (index_k = 0; index_k < pnl->k_size; index_k++){

    rk = pnl->k[index_k];

    if (rk > ppr->halofit_min_k_nonlinear) {

      pk_lin = pk_l[index_k]*pow(pnl->k[index_k],3)*anorm;

      /* in original halofit, this is the beginning of the function halofit() */

      /*SPB11: Standard halofit underestimates the power on the smallest
       * scales by a factor of two. Add an extra correction from the
       * simulations in Bird, Viel,Haehnelt 2011 which partially accounts for
       * this.*/
      /*SPB14: This version of halofit is an updated version of the fit to the massive neutrinos
       * based on the results of Takahashi 2012, (arXiv:1208.2701).
       */
      gam=0.1971-0.0843*rneff+0.8460*rncur;
      a=1.5222+2.8553*rneff+2.3706*rneff*rneff+0.9903*rneff*rneff*rneff+ 0.2250*rneff*rneff*rneff*rneff-0.6038*rncur+0.1749*Omega_v*(1.+w0);
      a=pow(10,a);
      b=pow(10, (-0.5642+0.5864*rneff+0.5716*rneff*rneff-1.5474*rncur+0.2279*Omega_v*(1.+w0)));
      c=pow(10, 0.3698+2.0404*rneff+0.8161*rneff*rneff+0.5869*rncur);
      xmu=0.;
      xnu=pow(10,5.2105+3.6902*rneff);
      alpha=fabs(6.0835+1.3373*rneff-0.1959*rneff*rneff-5.5274*rncur);
      beta=2.0379-0.7354*rneff+0.3157*pow(rneff,2)+1.2490*pow(rneff,3)+0.3980*pow(rneff,4)-0.1682*rncur + fnu*(1.081 + 0.395*pow(rneff,2));

      if(fabs(1-Omega_m)>0.01) { /*then omega evolution */
        f1a=pow(Omega_m,(-0.0732));
        f2a=pow(Omega_m,(-0.1423));
        f3a=pow(Omega_m,(0.0725));
        f1b=pow(Omega_m,(-0.0307));
        f2b=pow(Omega_m,(-0.0585));
        f3b=pow(Omega_m,(0.0743));
        frac=Omega_v/(1.-Omega_m);
        f1=frac*f1b + (1-frac)*f1a;
        f2=frac*f2b + (1-frac)*f2a;
        f3=frac*f3b + (1-frac)*f3a;
      }
      else {
        f1=1.;
        f2=1.;
        f3=1.;
      }

      y=(rk/rknl);
      pk_halo = a*pow(y,f1*3.)/(1.+b*pow(y,f2)+pow(f3*c*y,3.-gam));
      pk_halo=pk_halo/(1+xmu*pow(y,-1)+xnu*pow(y,-2))*(1+fnu*(0.977-18.015*(Omega0_m-0.3)));
      // rk is in 1/Mpc, 47.48and 1.5 in Mpc**-2, so we need an h**2 here (Credits Antonio J. Cuesta)
      pk_linaa=pk_lin*(1+fnu*47.48*pow(rk/pba->h,2)/(1+1.5*pow(rk/pba->h,2)));
      pk_quasi=pk_lin*pow((1+pk_linaa),beta)/(1+pk_linaa*alpha)*exp(-y/4.0-pow(y,2)/8.0);


        pk_nl[index_k] = (pk_halo+pk_quasi)/pow(pnl->k[index_k],3)/anorm;
  
//        printf("Misha is here.\n");

      /* in original halofit, this is the end of the function halofit() */
    }
    else {
      pk_nl[index_k] = pk_l[index_k];
    }
  }

  free(pvecback);
  free(integrand_array);
  return _SUCCESS_;
}

/* beginning of my own function */
 
 int nonlinear_pt(
                       struct precision *ppr,
                       struct background *pba,
                       struct primordial *ppm,
                       struct nonlinear *pnl,
                       double tau,
                       double *pk_l,
                       double *pk_nl,
                       double *lnk_l,
                       double *lnpk_l,
                       double *ddlnpk_l
                       ) {

int index_k;
double pk_lin,pk_quasi,pk_halo,rk;
double sigmav;
     
double * pvecback;

int last_index;
double sum0;
double anorm;

double *integrand_array;
int integrand_size;
int index_ia_k;
int index_ia_pk;
int index_ia_sum;
int index_ia_ddsum;

int ia_size;
int index_ia;

double k_integrand;
double lnpk_integrand;
     
//* halofit_found_k_max = ok;
     
class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);
     
anorm    = 1./(2*pow(_PI_,2));
 

index_ia=0;
class_define_index(index_ia_k,     _TRUE_,index_ia,1);
class_define_index(index_ia_pk,    _TRUE_,index_ia,1);
class_define_index(index_ia_sum,   _TRUE_,index_ia,1);
class_define_index(index_ia_ddsum, _TRUE_,index_ia,1);
ia_size = index_ia;

//M: here the integrand_array is tabulated. First one computes the size of the integrand (number of points I presume)

integrand_size=(int)(log(pnl->k[pnl->k_size-1]/pnl->k[0])/log(10.)*ppr->halofit_k_per_decade)+1;
class_alloc(integrand_array,integrand_size*ia_size*sizeof(double),pnl->error_message);

/* we fill integrand_array with values of k and P(k) using interpolation */

last_index=0;

for (index_k=0; index_k < integrand_size; index_k++) {
    
    k_integrand=pnl->k[0]*pow(10.,index_k/ppr->halofit_k_per_decade);
    
    class_call(array_interpolate_spline(lnk_l,
                                        pnl->k_size,
                                        lnpk_l,           /*what is being interpolated*/
                                        ddlnpk_l,
                                        1,
                                        log(k_integrand),
                                        &last_index,
                                        &lnpk_integrand,
                                        1,
                                        pnl->error_message),
               pnl->error_message,
               pnl->error_message);
    
    integrand_array[index_k*ia_size + index_ia_k] = k_integrand;
    integrand_array[index_k*ia_size + index_ia_pk] = exp(lnpk_integrand);
    
}

// class_call only calls a function!

class_call(background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback),
                pba->error_message,
                pnl->error_message);
printf("Call Halofit at z=%e\n",pba->a_today/pvecback[pba->index_bg_a]-1.);

class_call(nonlinear_p13UV_integrate(
                                     pnl,
                                     pba,
                                     integrand_array,
                                     integrand_size,
                                     ia_size,
                                     index_ia_k,
                                     index_ia_pk,
                                     index_ia_sum,
                                     index_ia_ddsum,
                                     &sum0
                                     ),
           pnl->error_message,
           pnl->error_message);
printf("Sigma_v=%f\n",sum0);
sigmav = sum0;
     
     int Nmax = 240;
     double Nmaxd = Nmax * 1.;
     double b = -0.3;
     double kcut = 68.;
     double kmin = 0.0001 * pba->h;
     double kmax = 50. * pba->h;
     
     //     double kmax = pnl->k[pnl->k_size-1];
     
/*
     if ( pnl->k[pnl->k_size-1] *1. < kcut ) {
         kmax = pnl->k[pnl->k_size-1];
     }
     else {
         kmax = kcut;
     }
 */
     
     // double kmin = pnl->k[0];
     double Delta = log(kmax / kmin) / (Nmax - 1);
     double js[Nmax+1], kdisc[Nmax], Pdisc[Nmax], Pbin[Nmax];
     double complex etam[Nmax+1];
     for (size_t i=0; i< Nmax +1 ; i++){
         js[i] = i - Nmaxd/2;
         etam[i] = b + 2 * M_PI * _Complex_I * js[i]/Nmaxd / Delta ;
     }
     
     for (size_t i=0; i< Nmax +1 ; i++){
         js[i] = i - Nmaxd/2;
         etam[i] = b + 2 * M_PI * _Complex_I * js[i]/Nmaxd / Delta ;
     //    printf("%e %e\n",creal(etam[i]),cimag(etam[i]));
     }
     
     for (size_t i=0; i < Nmax; i++) {
         
         kdisc[i] = kmin * exp(i * Delta);
         k_integrand = kdisc[i];
         
         class_call(array_interpolate_spline(lnk_l,
                                             pnl->k_size,
                                             lnpk_l,           
                                             ddlnpk_l,
//                                             lnpk_l,
                                             1,
                                             log(k_integrand),
                                             &last_index,
                                             &lnpk_integrand,
                                             1,
                                             pnl->error_message),
                    pnl->error_message,
                    pnl->error_message);
         
         Pdisc[i] = exp(lnpk_integrand) * exp(-1.* i * b* Delta);
         Pbin[i] = exp(lnpk_integrand);
         
//         printf("%f %f\n",kdisc[i],Pdisc[i]);
         
     }
     
     // here we compute the FFT coefficients
     
     /*
     int Nx = 240;
     double kx[Nx], Px[Nx];
     
     FILE *myFile2;
     myFile2 = fopen("/Users/michalychforever/Dropbox/Marko/P_test.dat", "r");
     
     if (myFile2 == NULL){
         printf("Error Reading File test...\n");
         exit (0);
     }
     
     for (size_t i=0; i < Nx; i++){
         fscanf(myFile2, "%lf%lf", &kx[i], &Px[i]);
     }
     fclose(myFile2);
     
    // fftw_complex in[Nx], out[Nx];
     
     fftw_complex *in, *out;
     fftw_plan p;
     
     in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
     out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
     
     for (size_t i=0; i< Nx ; i++){
         in[i][0] = Px[i];
         in[i][1] = 0;
     }
     
     p = fftw_plan_dft_1d(Nmax, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
     fftw_execute(p);
     
     for (size_t i=0; i< Nmax+1; i++){
         printf("%f %f\n",out[i][0],out[i][1]);
        }

     printf("End of cn's\n");
      */
     
     
      fftw_complex *in, *out;
      fftw_plan p;
     
     in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
     out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nmax);
     
//     double input_real1[Nmax], input_real2[Nmax];
     for (size_t i=0; i< Nmax ; i++){
         in[i][0] = Pdisc[i];
//         input_real1[i] = in[i][0];
         in[i][1] = 0;
//         input_real2[i] = 0;
     }
     
      p = fftw_plan_dft_1d(Nmax, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
      fftw_execute(p);
//      fftw_destroy_plan(p);
    
     
     /*
     double fft_coeff_real1[Nmax+1], fft_coeff_imag1[Nmax+1], fft_coeff_real2[Nmax+1], fft_coeff_imag2[Nmax+1];
     
     FFT_real_short(input_real1,
                    input_real2,
                    fft_coeff_real1,
                    fft_coeff_imag1,
                    fft_coeff_real2,
                    fft_coeff_imag2,
                    Nmax);
     */
     
     /*
     printf("Start of cn's\n");
     
     for (size_t i=0; i< Nmax+1; i++){
         printf("%f %f\n",out[i][0],out[i][1]);
     }
     
     printf("End of cn's\n");
     */
     
     
     // bad stuff happens here
     
     // decompose untill I get the same results as in my small c code
     
     double complex cmsym[Nmax+1];
     for (size_t i=0; i< Nmax+1; i++){
         if (i < Nmax/2) {
            cmsym[i]= cpow(kmin,-etam[i]) * (out[Nmax/2 - i][0] - _Complex_I * out[Nmax/2 - i][1])/Nmaxd;
         }
         else {
            cmsym[i]= cpow(kmin,-etam[i]) * (out[i - Nmax/2][0] + _Complex_I * out[i - Nmax/2][1])/Nmaxd;
         }
     }
     
     
     cmsym[0] = cmsym[0]/2.;
     cmsym[Nmax] = cmsym[Nmax]/2.;
     
     /*
     
     printf("Start of cnsym's\n");
     
     for (size_t i=0; i< Nmax+1; i++){
         printf("%f %f\n",creal(cmsym[i]),cimag(cmsym[i]));
     }
     
     printf("End of cnsym's\n");
    
     */
     
     
     /*
     double Pout[Nmax];
     for (size_t j=0; j < Nmax; j++){
         Pout[j] = 0.;
         for (size_t i=0; i < Nmax+1; i++){
             Pout[j] += creal(cmsym[i]* cpow(kdisc[j], etam[i]));
         }
         printf("%f %f\n",kdisc[j],Pout[j]);
     }
      */
     
     
     
     fftw_free(in);
     fftw_free(out);
     fftw_destroy_plan(p);
     
     
// here we input the PT matrices
     
     double M22_red_in[(Nmax+1)*(Nmax+2)/2][2];
     FILE *myFile22;
 //    myFile = fopen("pt_matrices/M22_red.dat", "r");
     
     myFile22 = fopen("/Users/michalychforever/Dropbox/Marko/M22_red.dat", "r");
     
     if (myFile22 == NULL){
         printf("Error Reading File M22_red\n");
         exit (0);
     }
     
     for (size_t i=0; i < (Nmax+1)*(Nmax+2)/2; i++){
         fscanf(myFile22, "%lf%lf", &M22_red_in[i][0], &M22_red_in[i][1]);
     //    printf("%e\n",M22_red_in[i][0]);
     }
     
     fclose(myFile22);
     
     double complex M22[Nmax+1][Nmax+1];
     
     int count=0;
     for (size_t i=0; i < Nmax+1; i++){
         for (size_t j=0; j <= i; j++){
             M22[i][j] = M22_red_in[count][0] + _Complex_I * M22_red_in[count][1];
             count++;
         }
     }
     
     for (size_t i=0; i < Nmax+1; i++){
         for (size_t j=i+1; j < Nmax+1; j++){
             M22[i][j] = M22[j][i];
         }
     }
     
     double M13in[Nmax+1][2];
     FILE *myFile13;
     
     myFile13 = fopen("/Users/michalychforever/Dropbox/Marko/M13.dat", "r");
     
     if (myFile13 == NULL){
         printf("Error Reading File M13\n");
         exit (0);
     }
     
     for (size_t i=0; i < Nmax+1; i++){
         fscanf(myFile13, "%lf%lf", &M13in[i][0], &M13in[i][1]);
     }
     
     fclose(myFile13);

     
     double complex M13[Nmax+1];
     for (size_t i=0; i < Nmax+1; i++){
         M13[i] = M13in[i][0]+ _Complex_I * M13in[i][1];
     }
     
     //End of inputting PT matrices
     
     /* Multiplying the matrices */
     
     double cutoff = 10. * pba->h;
     
     double complex b13[Nmax+1], f13[Nmax];
     double P13[Nmax], P13UV[Nmax];
     for (size_t j=0; j < Nmax; j++){
         for (size_t i=0; i < Nmax+1; i++){
             b13[i] = cmsym[i] * cpow(kdisc[j], etam[i]);
//             printf("%e\n", creal(b13[i]));
             f13[j] += b13[i] * M13[i];
         }
         P13UV[j] = - 61. * Pbin[j] * pow(kdisc[j], 2.) * sigmav / 105.;
         P13[j] = (creal(cpow(kdisc[j], 3) * Pbin[j] * f13[j]) + P13UV[j]) * exp(-pow(kdisc[j]/cutoff, 6.));
//                  printf("%f  %f\n",kdisc[j],P13[j]);
     }
     
     
     double complex b22[Nmax+1], f22[Nmax];
     double P22[Nmax], P1loop[Nmax];
     for (size_t j=0; j < Nmax; j++){
         for (size_t i=0; i < Nmax+1; i++){
             b22[i] = cmsym[i] * cpow(kdisc[j], etam[i]);
             f22[j] += b22[i] * M22[i][i] * b22[i];
                          for (size_t l=0; l<i; l++){
                              f22[j] += b22[i] * M22[i][l] * b22[l] * 2.;
                          }
         }
         //         printf("%e\n", creal(f22[j]));
         P22[j] = creal(pow(kdisc[j], 3) * f22[j] * exp(-pow(kdisc[j]/cutoff, 6)) );
         P1loop[j] = 0*Pbin[j] + P13[j] + P22[j];
        printf("%e  %e\n",kdisc[j],P1loop[j]);
     }
     
     

     
      //        printf("%e  %e\n",kdisc[i],Pbin[i]);

          // loop over k here - Multiply matrices here then
    
     
     double kinterp1loop;
     double pk_1loop;
     
     last_index=0;
     for (index_k=0; index_k < pnl->k_size; index_k++) {
         
         kinterp1loop = pnl->k[index_k];
         
         class_call(array_interpolate_spline(kdisc,
                                             Nmax,
                                             P1loop,
                                             P1loop,
                                             1,
                                             kinterp1loop,
                                             &last_index,
                                             &pk_1loop,
                                             1,
                                             pnl->error_message),
                    pnl->error_message,
                    pnl->error_message);
         
         pk_lin = pk_l[index_k]*pow(pnl->k[index_k],3)*anorm;
         pk_nl[index_k] = pk_1loop;
     }
     
/*
     
     for (index_k = 0; index_k < pnl->k_size; index_k++){
         pk_lin = pk_l[index_k]*pow(pnl->k[index_k],3)*anorm;
         pk_nl[index_k] = (pk_lin - pow(pnl->k[index_k],2) * pk_lin * sigmav)/pow(pnl->k[index_k],3)/anorm;
     }
 */
     
     

free(pvecback);
free(integrand_array);
return _SUCCESS_;
}

/* end of my own function */



/* in original halofit, this is equivalent to the function wint() */

// This function integrates
int nonlinear_halofit_integrate(
                                struct nonlinear *pnl,
                                double * integrand_array,
                                int integrand_size,
                                int ia_size,
                                int index_ia_k,
                                int index_ia_pk,
                                int index_ia_sum,
                                int index_ia_ddsum,
                                double R,
                                enum halofit_integral_type type,
                                double * sum
                                ) {

  double k,pk,x2,integrand;
  int index_k;
  double anorm = 1./(2*pow(_PI_,2));

    for (index_k=0; index_k < integrand_size; index_k++) {
      k = integrand_array[index_k*ia_size + index_ia_k];
      pk = integrand_array[index_k*ia_size + index_ia_pk];
      x2 = k*k*R*R;

      integrand = pk*k*k*anorm*exp(-x2);
      if (type == halofit_integral_two) integrand *= 2.*x2;
      if (type == halofit_integral_three) integrand *= 4.*x2*(1.-x2);

      integrand_array[index_k*ia_size + index_ia_sum] = integrand;
    }

    /* fill in second derivatives */
    class_call(array_spline(integrand_array,
                            ia_size,
                            integrand_size,
                            index_ia_k,
                            index_ia_sum,
                            index_ia_ddsum,
                            _SPLINE_NATURAL_,
                            pnl->error_message),
               pnl->error_message,
               pnl->error_message);

    /* integrate */
    class_call(array_integrate_all_spline(integrand_array,
                                          ia_size,
                                          integrand_size,
                                          index_ia_k,
                                          index_ia_sum,
                                          index_ia_ddsum,
                                          sum,
                                          pnl->error_message),
               pnl->error_message,
               pnl->error_message);

  return _SUCCESS_;
}

// when you define the arguments of the function, don't forget to also change $nonlinear.h$ file

int nonlinear_p13UV_integrate(
                                struct nonlinear *pnl,
                                struct background *pba,
                                double * integrand_array,
                                int integrand_size,
                                int ia_size,
                                int index_ia_k,
                                int index_ia_pk,
                                int index_ia_sum,
                                int index_ia_ddsum,
                                double * sum
                                ) {
    
    double k,pk,integrand;
    int index_k;
    double bnorm = 1./(6. * pow(_PI_,2)) * pow(1 * pba->h,0);
    
    for (index_k=0; index_k < integrand_size; index_k++) {
        k = integrand_array[index_k*ia_size + index_ia_k];
        pk = integrand_array[index_k*ia_size + index_ia_pk];
        
        integrand = pk * bnorm;
      //  printf("%e  %e\n",k,pk);
        integrand_array[index_k*ia_size + index_ia_sum] = integrand;
    }
    
    /* fill in second derivatives */
    class_call(array_spline(integrand_array,
                            ia_size,
                            integrand_size,
                            index_ia_k,
                            index_ia_sum,
                            index_ia_ddsum,
                            _SPLINE_NATURAL_,
                            pnl->error_message),
               pnl->error_message,
               pnl->error_message);
    
    /* integrate */
    class_call(array_integrate_all_spline(integrand_array,
                                          ia_size,
                                          integrand_size,
                                          index_ia_k,
                                          index_ia_sum,
                                          index_ia_ddsum,
                                          sum,
                                          pnl->error_message),
               pnl->error_message,
               pnl->error_message);
    
    return _SUCCESS_;
    
    
    
    /** @file fft.c Documented fast fourier transform module
     *
     * Nils Schoenberg, 12.10.2017
     *
     * This module computes the fast fourier transform (FFT) of any function
     */
#include "nonlinear.h"
#include <math.h>
    //Mathematical constants used in this file
#define MATH_PI_2 6.2831853071795864769252867665590057683943387987502
#define INV_SQRT_2 1/sqrt(2)
    //Function implementations
    /**
     * Computes the fast fourier transform of some arbitrary input array input_real and input_double of sizes N
     * Returns the output by writing into output_real and output_imag
     *
     * It is assumed that all the arrays are allocated and of size N
     *
     * For recursion there is a stepsize parameter
     * If the full FFT should be calculated, set
     *
     *          ** stepsize = 1 **
     *
     * */
    void FFT(double* input_real, double* input_imag, double* output_real, double* output_imag, int N, int stepsize){
        //Larger base cases up to N == 8 have proven advantageous
        if (N == 8){
            //FFT N==4 EVEN
            //FFT(input_real, input_imag, output_real, output_imag, 4, 2 * stepsize);
            // i = 0
            double temp_even_real = input_real[0] + input_real[4 * stepsize];
            double temp_even_imag = input_imag[0] + input_imag[4 * stepsize];
            double temp_odd_real = input_real[2*stepsize] + input_real[6 * stepsize];
            double temp_odd_imag = input_imag[2*stepsize] + input_imag[6 * stepsize];
            output_real[0] = temp_even_real + temp_odd_real;
            output_imag[0] = temp_even_imag + temp_odd_imag;
            output_real[2] = temp_even_real - temp_odd_real;
            output_imag[2] = temp_even_imag - temp_odd_imag;
            // i = 1
            temp_even_real = input_real[0] - input_real[4 * stepsize];
            temp_even_imag = input_imag[0] - input_imag[4 * stepsize];
            temp_odd_real = input_real[2*stepsize] - input_real[6 * stepsize];
            temp_odd_imag = input_imag[2*stepsize] - input_imag[6 * stepsize];
            output_real[1] = temp_even_real + temp_odd_imag;
            output_imag[1] = temp_even_imag - temp_odd_real;
            output_real[3] = temp_even_real - temp_odd_imag;
            output_imag[3] = temp_even_imag + temp_odd_real;
            //FFT N==4 ODD
            //FFT(input_real + stepsize, input_imag + stepsize, output_real + 4, output_imag + 4, 4, 2 * stepsize);
            // i = 0
            temp_even_real = input_real[stepsize] + input_real[5 * stepsize];
            temp_even_imag = input_imag[stepsize] + input_imag[5 * stepsize];
            temp_odd_real = input_real[3*stepsize] + input_real[7 * stepsize];
            temp_odd_imag = input_imag[3*stepsize] + input_imag[7 * stepsize];
            output_real[4] = temp_even_real + temp_odd_real;
            output_imag[4] = temp_even_imag + temp_odd_imag;
            output_real[6] = temp_even_real - temp_odd_real;
            output_imag[6] = temp_even_imag - temp_odd_imag;
            // i = 1
            temp_even_real = input_real[stepsize] - input_real[5 * stepsize];
            temp_even_imag = input_imag[stepsize] - input_imag[5 * stepsize];
            temp_odd_real = input_real[3*stepsize] - input_real[7 * stepsize];
            temp_odd_imag = input_imag[3*stepsize] - input_imag[7 * stepsize];
            output_real[5] = temp_even_real + temp_odd_imag;
            output_imag[5] = temp_even_imag - temp_odd_real;
            output_real[7] = temp_even_real - temp_odd_imag;
            output_imag[7] = temp_even_imag + temp_odd_real;
            //FINAL FFT N==8
            //i==0
            temp_even_real = output_real[0];
            temp_even_imag = output_imag[0];
            temp_odd_real = output_real[4];
            temp_odd_imag = output_imag[4];
            output_real[0] = temp_even_real + temp_odd_real;
            output_imag[0] = temp_even_imag + temp_odd_imag;
            output_real[4] = temp_even_real - temp_odd_real;
            output_imag[4] = temp_even_imag - temp_odd_imag;
            //i==1
            temp_even_real = output_real[1];
            temp_even_imag = output_imag[1];
            temp_odd_real = output_real[5];
            temp_odd_imag = output_imag[5];
            output_real[1] = temp_even_real + INV_SQRT_2*temp_odd_real + INV_SQRT_2*temp_odd_imag;
            output_imag[1] = temp_even_imag + INV_SQRT_2*temp_odd_imag - INV_SQRT_2*temp_odd_real;
            output_real[5] = temp_even_real - INV_SQRT_2*temp_odd_real - INV_SQRT_2*temp_odd_imag;
            output_imag[5] = temp_even_imag + INV_SQRT_2*temp_odd_real - INV_SQRT_2*temp_odd_imag;
            //i==2
            temp_even_real = output_real[2];
            temp_even_imag = output_imag[2];
            temp_odd_real = output_real[6];
            temp_odd_imag = output_imag[6];
            output_real[2] = temp_even_real + temp_odd_imag;
            output_imag[2] = temp_even_imag - temp_odd_real;
            output_real[6] = temp_even_real - temp_odd_imag;
            output_imag[6] = temp_even_imag + temp_odd_real;
            //i==3
            temp_even_real = output_real[3];
            temp_even_imag = output_imag[3];
            temp_odd_real = output_real[7];
            temp_odd_imag = output_imag[7];
            output_real[3] = temp_even_real - INV_SQRT_2*temp_odd_real + INV_SQRT_2*temp_odd_imag;
            output_imag[3] = temp_even_imag - INV_SQRT_2*temp_odd_imag - INV_SQRT_2*temp_odd_real;
            output_real[7] = temp_even_real + INV_SQRT_2*temp_odd_real - INV_SQRT_2*temp_odd_imag;
            output_imag[7] = temp_even_imag + INV_SQRT_2*temp_odd_real + INV_SQRT_2*temp_odd_imag;
            return;
        }
        else if (N < 8){
            if (N == 4){
                // i = 0
                double temp_even_real = input_real[0] + input_real[2 * stepsize];
                double temp_even_imag = input_imag[0] + input_imag[2 * stepsize];
                double temp_odd_real = input_real[stepsize] + input_real[3 * stepsize];
                double temp_odd_imag = input_imag[stepsize] + input_imag[3 * stepsize];
                output_real[0] = temp_even_real + temp_odd_real;
                output_imag[0] = temp_even_imag + temp_odd_imag;
                output_real[2] = temp_even_real - temp_odd_real;
                output_imag[2] = temp_even_imag - temp_odd_imag;
                // i = 1
                temp_even_real = input_real[0] - input_real[2 * stepsize];
                temp_even_imag = input_imag[0] - input_imag[2 * stepsize];
                temp_odd_real = input_real[stepsize] - input_real[3 * stepsize];
                temp_odd_imag = input_imag[stepsize] - input_imag[3 * stepsize];
                output_real[1] = temp_even_real + temp_odd_imag;
                output_imag[1] = temp_even_imag - temp_odd_real;
                output_real[3] = temp_even_real - temp_odd_imag;
                output_imag[3] = temp_even_imag + temp_odd_real;
                return;
            }
            else
                if (N == 2){
                    output_real[0] = input_real[0] + input_real[stepsize];
                    output_real[1] = input_real[0] - input_real[stepsize];
                    output_imag[0] = input_imag[0] + input_imag[stepsize];
                    output_imag[1] = input_imag[0] - input_imag[stepsize];
                    return;
                }
                else
                    if (N == 1){
                        output_real[0] = input_real[0];
                        output_imag[0] = input_imag[0];
                        return;
                    }
        }
        else{
            //Use the butterfly algorithm to compute the fft of even and odd sets
            FFT(input_real, input_imag, output_real, output_imag, N / 2, 2 * stepsize);
            FFT(input_real + stepsize, input_imag + stepsize, output_real + N / 2, output_imag + N / 2, N / 2, 2 * stepsize);
            //Reunite even and odd sets
            for (int i = 0 , j = N/2; i < N / 2; ++i,++j){
                double temp_even_real = output_real[i];
                double temp_even_imag = output_imag[i];
                double temp_odd_real = output_real[j];
                double temp_odd_imag = output_imag[j];
                //These twiddle factors cos_val and sin_val could be read from file instead
                //It will depend on many things whether or not that is advantagous
                // (e.g. on the number of frequqencies N)
                double cos_val = cos(i*MATH_PI_2 / ((double)(N)));
                double sin_val = sin(i*MATH_PI_2 / ((double)(N)));
                output_real[i] = temp_even_real + cos_val*temp_odd_real + sin_val*temp_odd_imag;
                output_imag[i] = temp_even_imag + cos_val*temp_odd_imag - sin_val*temp_odd_real;
                output_real[j] = temp_even_real - cos_val*temp_odd_real - sin_val*temp_odd_imag;
                output_imag[j] = temp_even_imag + sin_val*temp_odd_real - cos_val*temp_odd_imag;
            }
            return;
        }
    }
    /**
     * Computes the fast fourier transform of some purely real inputs input_real_1 and input_real_2 of size N
     * Returns the output by writing into output_real_i and output_imag_i for input_real_i with i=1,2
     *
     * It is assumed that all the arrays are allocated and of size N
     *
     * Returns full output, arrays of size N
     * */
    void FFT_real(double* input_real_1, double* input_real_2, double* output_real_1,double* output_imag_1, double* output_real_2,double* output_imag_2, int N){
        FFT(input_real_1, input_real_2, output_real_1, output_real_2, N,1);
        //output_real_1[0] remains the same
        double temp1,temp2;
        output_imag_1[0] = 0.0;
        output_imag_2[0] = 0.0;
        output_imag_1[N/2] = 0.0;
        output_imag_2[N/2] = 0.0;
        for (int i = 1; i < N/2; ++i){
            temp1 = output_real_1[i];
            temp2 = output_real_2[i];
            output_real_1[i] = 0.5*(temp1 + output_real_1[N - i]);
            output_real_2[i] = 0.5*(temp2 + output_real_2[N - i]);
            output_imag_1[i] = 0.5*(temp2 - output_real_2[N - i]);
            output_imag_2[i] = 0.5*(output_real_1[N - i] - temp1);
        }
        for (int i = 0; i < N / 2; ++i){
            output_real_1[i + N / 2] = output_real_1[N / 2 - i];
            output_real_2[i + N / 2] = output_real_2[N / 2 - i];
            output_imag_1[i + N / 2] = - output_imag_1[N / 2 - i];
            output_imag_2[i + N / 2] = - output_imag_2[N / 2 - i];
        }
    }
    /**
     * Computes the fast fourier transform of some purely real inputs input_real_1 and input_real_2 of size N
     * Returns the output by writing into output_real_i and output_imag_i for input_real_i with i=1,2
     *
     * It is assumed that all the arrays are allocated and of size N
     * 
     * Only returns N/2 output arrays (still have to be allocated at size N)
     * 
     * For any real fourier transformation c_(-n) = c_n and thus c_(N-n) = c_n for finite fourier transformation of size N
     * */
    void FFT_real_short(double* input_real_1, double* input_real_2, double* output_real_1, double* output_imag_1, double* output_real_2, double* output_imag_2, int N){
        //Only computes first N/2 elements, since others are related by symmetry
        FFT(input_real_1, input_real_2, output_real_1, output_real_2, N, 1);
        double temp1,temp2;
        output_imag_1[0] = 0.0;
        output_imag_2[0] = 0.0;
        output_imag_1[N/2] = 0.0;
        output_imag_2[N/2] = 0.0;
        for (int i = 1; i < N/2; ++i){
            temp1 = output_real_1[i];
            temp2 = output_real_2[i];
            output_real_1[i] = 0.5*(temp1 + output_real_1[N - i]);
            output_real_2[i] = 0.5*(temp2 + output_real_2[N - i]);
            output_imag_1[i] = 0.5*(temp2 - output_real_2[N - i]);
            output_imag_2[i] = 0.5*(output_real_1[N - i] - temp1);
        }
    }
    /**
     * Computes the discrete cosine transform of some purely real input input of size N
     * Returns the output by writing into output_real and output_imag
     * 
     * CAUTION :: It is assumed that all the arrays are allocated and of size 2*N  !!
     * 
     * The DCT can be found in the output_real array in the first N elements
     * */
#include <stdio.h>
    void DCT_real(double* input_real,double* input_imag,double* output_real,double* output_imag,int N){
        int i;double c,s;double temp;
        for(i=0;i<N;++i){
            input_real[i+N]=0;
            input_imag[i]=0;
            input_imag[i+N]=0;
        }
        FFT(input_real,input_imag,output_real,output_imag,2*N,1);
        for(i=0;i<N;++i){
            c = cos(-0.25*i*MATH_PI_2/N);
            s = sin(-0.25*i*MATH_PI_2/N);
            temp = output_real[i];
            output_real[i] = c*temp-s*output_imag[i];
            output_imag[i] = c*output_imag[i]+s*temp;
        }
    }

    
}

