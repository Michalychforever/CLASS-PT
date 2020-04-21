/** @file nonlinear_pt.c Documented nonlinear perturbation theory module
 *
 * Mikhail M. Ivanov, 17.06.2018
 *
 * Here is a module performing one-loop perturbation theory calculations for dark matter and biased tracers
 * with IR resummation in real space and redshift spaces. We use the FFTLog for a fast evaluation of PT loop integrals in the EdS approximation.
 * TODO: exact time dependence, one-loop bispecrum and two-loop power spectrum.
 *
 */
#include "nonlinear_pt.h"
#include "fft.h"

int perturb_get_k_list_nl(
                        struct precision * ppr,
                        struct background * pba,
                        struct thermo * pth,
                        struct perturbs * ppt,
                        struct nonlinear_pt *pnlpt
                        ) {
  int index_k, index_k_output, index_mode;
  double k,k_min=0.,k_rec,step,tau1;
  double k_max_cmb;
  double k_max_cl;
  double k_max=0.;
  double scale2;
  double * tmp_k_list;
  int newk_size, index_newk, add_k_output_value;

  /** Summary: */

  class_test(ppr->k_step_transition == 0.,
             pnlpt->error_message,
             "stop to avoid division by zero");

  class_test(pth->rs_rec == 0.,
             pnlpt->error_message,
             "stop to avoid division by zero");

  /** - scalar modes */

  if (ppt->has_scalars == _TRUE_) {


    k_min=ppr->k_min_tau0/pba->conformal_age;

    /* first value */
    if (pba->sgnK == 0) {
      /* K<0 (flat)  : start close to zero */
      k_min=ppr->k_min_tau0/pba->conformal_age;
    }

    else if (pba->sgnK == -1) {
       // K<0 (open)  : start close to sqrt(-K)
         // (in transfer modules, for scalars, this will correspond to q close to zero;
         // for vectors and tensors, this value is even smaller than the minimum necessary value) 
      k_min=sqrt(-pba->K+pow(ppr->k_min_tau0/pba->conformal_age/pth->angular_rescaling,2));

    }
    else if (pba->sgnK == 1) {
      /* K>0 (closed): start from q=sqrt(k2+(1+m)K) equal to 3sqrt(K), i.e. k=sqrt((8-m)K) */
      k_min = sqrt((8.-1.e-4)*pba->K);
    }

    // printf("k_min=%le\n",k_min);
    /** - --> find k_max (as well as k_max_cmb[ppt->index_md_scalars], k_max_cl[ppt->index_md_scalars]) */

    k_rec = 2. * _PI_ / pth->rs_rec; /* comoving scale corresponding to sound horizon at recombination */

    k_max_cmb = k_min;
    k_max_cl = k_min;
    k_max = k_min;

    /* find k_max: */

    if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_))
      k_max = MAX(k_max,ppt->k_max_for_pk);

    if (ppt->has_nl_corrections_based_on_delta_m == _TRUE_)
      k_max = MAX(k_max,ppr->halofit_min_k_max);

    // printf("k_max=%le\n",k_max);

    /** - --> test that result for k_min, k_max make sense */

    class_test(k_min<0.,
               pnlpt->error_message,
               "buggy definition of k_min");

    class_test(k_max<0.,
               pnlpt->error_message,
               "buggy definition of k_max");

    class_test(k_max<k_min,
               pnlpt->error_message,
               "buggy definition of k_min and/or k_max");

    /* if K>0, the transfer function will be calculated for discrete
       integer values of nu=3,4,5,... where nu=sqrt(k2+(1+m)K) and
       m=0,1,2 for scalars/vectors/tensors. However we are free to
       define in the perturbation module some arbitrary values of k:
       later on, the transfer module will interpolate at values of k
       corresponding exactly to integer values of nu. Hence, apart
       from the value of k_min and the step size in the vicinity of
       k_min, we define exactly the same sampling in the three cases
       K=0, K<0, K>0 */

    /* allocate array with, for the moment, the largest possible size */
    class_alloc(pnlpt->k,
                ((int)((k_max_cmb-k_min)/k_rec/MIN(ppr->k_step_super,ppr->k_step_sub))+
                 (int)(MAX(ppr->k_per_decade_for_pk,ppr->k_per_decade_for_bao)*log(k_max/k_min)/log(10.))+3)
                *sizeof(double),pnlpt->error_message);

    /* first value */

    index_k=0;
    k = k_min;
    pnlpt->k[index_k] = k;
    index_k++;

    /* values until k_max_cmb[ppt->index_md_scalars] */

    while (k < k_max_cmb) {

      /* the linear step is not constant, it has a step-like shape,
         centered around the characteristic scale set by the sound
         horizon at recombination (associated to the comoving wavenumber
         k_rec) */

      step = (ppr->k_step_super
              + 0.5 * (tanh((k-k_rec)/k_rec/ppr->k_step_transition)+1.)
              * (ppr->k_step_sub-ppr->k_step_super)) * k_rec;

      /* there is one other thing to take into account in the step
         size. There are two other characteristic scales that matter for
         the sampling: the Hubble scale today, k0=a0H0, and eventually
         curvature scale sqrt(|K|). We define "scale2" as the sum of the
         squared Hubble radius and squared curvature radius. We need to
         increase the sampling for k<sqrt(scale2), in order to get the
         first mutipoles accurate enough. The formula below reduces it
         gradually in the k-->0 limit, by up to a factor 10. The actual
         stepsize is still fixed by k_step_super, this is just a
         reduction factor. */

      scale2 = pow(pba->a_today*pba->H0,2)+fabs(pba->K);

      step *= (k*k/scale2+1.)/(k*k/scale2+1./ppr->k_step_super_reduction);

      class_test(step / k < ppr->smallest_allowed_variation,
                 pnlpt->error_message,
                 "k step =%e < machine precision : leads either to numerical error or infinite loop",
                 step * k_rec);

      k += step;

      class_test(k <= pnlpt->k[index_k-1],
                 pnlpt->error_message,
                 "consecutive values of k should differ and should be in growing order");

      pnlpt->k[index_k] = k;

      index_k++;
    }

    pnlpt->k_size_cmb = index_k;

    /* values until k_max_cl[ppt->index_md_scalars] */

    while (k < k_max_cl) {

      k *= pow(10.,1./(ppr->k_per_decade_for_pk
                       +(ppr->k_per_decade_for_bao-ppr->k_per_decade_for_pk)
                       *(1.-tanh(pow((log(k)-log(ppr->k_bao_center*k_rec))/log(ppr->k_bao_width),4)))));

      pnlpt->k[index_k] = k;
      index_k++;
    }

    pnlpt->k_size_cl = index_k;

    /* values until k_max */

    while (k < k_max) {

      k *= pow(10.,1./(ppr->k_per_decade_for_pk
                       +(ppr->k_per_decade_for_bao-ppr->k_per_decade_for_pk)
                       *(1.-tanh(pow((log(k)-log(ppr->k_bao_center*k_rec))/log(ppr->k_bao_width),4)))));

      pnlpt->k[index_k] = k;
      index_k++;
    }

    pnlpt->k_size = index_k;

    class_realloc(pnlpt->k,
                  pnlpt->k,
                  pnlpt->k_size*sizeof(double),
                  pnlpt->error_message);
  }

  /** - vector modes skipped */ 

  /** - tensor modes skipped */

  /** - If user asked for k_output_values, add those to all k lists: */
  if (ppt->k_output_values_num>0){
    /* Allocate storage */
    class_alloc(ppt->index_k_output_values,sizeof(double)*1*ppt->k_output_values_num,pnlpt->error_message);

    /** - --> Find indices in ppt->k[index_md] corresponding to 'k_output_values'.
        We are assuming that ppt->k is sorted and growing, and we have made sure
        that ppt->k_output_values is also sorted and growing.*/
    //for (index_mode=0; index_mode<ppt->md_size; index_mode++){

      newk_size = pnlpt->k_size+ppt->k_output_values_num;

      class_alloc(tmp_k_list,sizeof(double)*newk_size,pnlpt->error_message);

      index_k=0;
      index_k_output=0;
      for (index_newk=0; index_newk<newk_size; index_newk++){
        /** - --> Decide if we should add k_output_value now. This has to be this complicated, since we
            can only compare the k-values when both indices are in range.*/
        if (index_k >= pnlpt->k_size)
          add_k_output_value = _TRUE_;
        else if (index_k_output >= ppt->k_output_values_num)
          add_k_output_value = _FALSE_;
        else if (ppt->k_output_values[index_k_output] < pnlpt->k[index_k])
          add_k_output_value = _TRUE_;
        else
          add_k_output_value = _FALSE_;

        if (add_k_output_value == _TRUE_){
          tmp_k_list[index_newk] = ppt->k_output_values[index_k_output];
          ppt->index_k_output_values[0*ppt->k_output_values_num+index_k_output]=index_newk;
          index_k_output++;
        }
        else{
          tmp_k_list[index_newk] = pnlpt->k[index_k];
          index_k++;
        }
      }

      free(pnlpt->k);
      pnlpt->k = tmp_k_list;
      pnlpt->k_size = newk_size;

      index_k = newk_size-1;
      while (pnlpt->k[index_k] > k_max_cl)
        index_k--;
      pnlpt->k_size_cl = MIN(index_k+2,pnlpt->k_size);

      index_k = newk_size-1;
      while (pnlpt->k[index_k] > k_max_cmb)
        index_k--;
      pnlpt->k_size_cmb = MIN(index_k+2,pnlpt->k_size);

      /** - --> The two MIN statements are here because in a normal run, the cl and cmb
          arrays contain a single k value larger than their respective k_max.
          We are mimicking this behavior. */
    //}
  }

  /* For testing, can be useful to print the k list in a file:

  FILE * out=fopen("output/k","w");

  for (index_k=0; index_k < ppt->k_size[0]; index_k++) {

    fprintf(out,"%e\n",ppt->k[0][index_k],pba->K);

  }
     fclose(out);
  */

  /** - finally, find the global k_min and k_max for the ensemble of all modes 9scalars, vectors, tensors) */

  pnlpt->k_min = _HUGE_;
  pnlpt->k_max = 0.;
  if (ppt->has_scalars == _TRUE_) {
    pnlpt->k_min = MIN(pnlpt->k_min,pnlpt->k[0]); /* first value, inferred from perturbations structure */
    pnlpt->k_max = MAX(pnlpt->k_max,pnlpt->k[pnlpt->k_size-1]); /* last value, inferred from perturbations structure */
  }

  //free(k_max_cmb);
  //free(k_max_cl);

  return _SUCCESS_;

}

int nonlinear_pt_init(
                   struct precision *ppr,
                   struct background *pba,
                   struct thermo *pth,
                   struct perturbs *ppt,
                   struct primordial *ppm,
                   struct nonlinear_pt *pnlpt
                   ) {

  int index_ncdm;
  int index_k;
  int index_tau;
  double *pk_l;
  double *pk_l_0_vv;
  double *pk_l_0_vd;
  double *pk_l_0_dd;
    double *pk_l_2_vv;
    double *pk_l_2_vd;
    double *pk_l_2_dd;
    double *pk_l_4_vv;
    double *pk_l_4_vd;
    double *pk_l_4_dd;
  double *pk_nl;
  double *pk_Id2d2;
  double *pk_Id2d2_2;
  double *pk_Id2d2_4;
  double *pk_l_0_b1b2;
  double *pk_l_0_b2;
  double *pk_l_0_b1bG2;
    double *pk_l_0_bG2;
    double *pk_l_2_b1b2;
    double *pk_l_2_b2;
    double *pk_l_2_b1bG2;
    double *pk_l_2_bG2;
    double *pk_l_4_b2;
    double *pk_l_4_bG2;
    double *pk_l_4_b1b2;
    double *pk_l_4_b1bG2;
  double *pk_Id2;
  double *pk_IG2;
  double *pk_Id2G2;
    double *pk_Id2G2_2;
    double *pk_Id2G2_4;
  double *pk_IG2G2;
    double *pk_IG2G2_2;
    double *pk_IG2G2_4;
  double *pk_IFG2;
    double *pk_IFG2_0b1;
    double *pk_IFG2_0;
    double *pk_IFG2_2;
  double *pk_CTR;
  double *pk_CTR_0;
  double *pk_CTR_2;
  double *pk_CTR_4;
  double *pk_Tree;
    double *pk_Tree_0_vv;
    double *pk_Tree_0_vd;
    double *pk_Tree_0_dd;
    double *pk_Tree_2_vv;
    double *pk_Tree_2_vd;
    double *pk_Tree_4_vv;
  double *lnk_l;
  double *lnpk_l;
  double *ddlnpk_l;
  double *tau_req;
  short print_warning=_FALSE_;
  //double * pvecback;
  int last_index;
  double a,z;

last_index=0;

  /** Summary
   *
   * (a) First deal with the case where non non-linear corrections requested */

  if (pnlpt->method == nlpt_none) {
    if (pnlpt->nonlinear_pt_verbose > 0)
      printf("No non-linear spectra requested. Nonlinear PT module skipped.\n");
  }

  /** (b) Compute the non-linear spectrum */

  else if (pnlpt->method == nlpt_spt) {
    if (pnlpt->nonlinear_pt_verbose > 0)
      printf("Computing non-linear matter and bias power spectra in perturbation theory.\n");

      
    /** - copy list of (k,tau) from perturbation module */

    pnlpt->index_md_scalars=ppt->index_md_scalars;
    pnlpt->ic_size=ppt->ic_size[ppt->index_md_scalars];
    pnlpt->tp_size=ppt->tp_size[ppt->index_md_scalars];

    if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {

      class_call(perturb_get_k_list_nl(ppr,
                                pba,
                                pth,
                                ppt,
                                pnlpt),
             pnlpt->error_message,
             pnlpt->error_message);

    } else {

      pnlpt->k_size = ppt->k_size[ppt->index_md_scalars];
      class_alloc(pnlpt->k,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      for (index_k=0; index_k<pnlpt->k_size; index_k++){
        pnlpt->k[index_k] = ppt->k[ppt->index_md_scalars][index_k];
      }

    }

    pnlpt->ln_k_size = pnlpt->k_size;
    class_alloc(pnlpt->ln_k,pnlpt->ln_k_size*sizeof(double),pnlpt->error_message);
    for (index_k=0; index_k<pnlpt->k_size; index_k++){
      pnlpt->ln_k[index_k] = log(pnlpt->k[index_k]);
    }


       

      pnlpt->tau_size = ppt->tau_size;
      class_alloc(pnlpt->tau,pnlpt->tau_size*sizeof(double),pnlpt->error_message);
      class_alloc(pnlpt->ln_tau,pnlpt->tau_size*sizeof(double),pnlpt->error_message);
      for (index_tau=0; index_tau<pnlpt->tau_size; index_tau++){
          pnlpt->tau[index_tau] = ppt->tau_sampling[index_tau];
          pnlpt->ln_tau[index_tau] = log(ppt->tau_sampling[index_tau]);
      }

      
      /*
      pnlpt->tau_size = pnlpt->z_pk_num;
      
    class_alloc(pnlpt->tau,pnlpt->tau_size*sizeof(double),pnlpt->error_message);
      
      for (index_tau=0; index_tau<pnlpt->tau_size; index_tau++){
          
          class_call(background_tau_of_z(pba,pnlpt->z_pk[pnlpt->tau_size - index_tau - 1],&pnlpt->tau[index_tau]),
                     pba->error_message,
                     pnlpt->error_message);
          printf("%lf %lf\n",pnlpt->tau[index_tau],pnlpt->z_pk[pnlpt->tau_size - index_tau - 1]);
          
      }
       */
    
//    class_alloc(pnlpt->pk_nl_out,pnlpt->tau_size*pnlpt->k_size*sizeof(double),pnlpt->error_message);
  
/*Needed for transfers*/
class_alloc(pnlpt->nl_corr_density,pnlpt->tau_size*ppt->k_size[pnlpt->index_md_scalars]*sizeof(double),pnlpt->error_message);

    class_alloc(pk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_l_0_vv,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_l_0_vd,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_l_0_dd,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      
      class_alloc(pk_l_2_vv,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_l_2_vd,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_l_2_dd,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      
      class_alloc(pk_l_4_vv,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_l_4_vd,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_l_4_dd,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      
      class_alloc(pk_l_0_b1b2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_l_0_b2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_l_0_b1bG2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_l_0_bG2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      
      class_alloc(pk_l_2_b1b2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_l_2_b2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_l_2_b1bG2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_l_2_bG2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      
    class_alloc(pk_l_4_b2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_l_4_bG2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_l_4_b1b2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_l_4_b1bG2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      
    class_alloc(pk_nl,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_Id2d2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_Id2d2_2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_Id2d2_4,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      
    class_alloc(pk_Id2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_IG2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_Id2G2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
       class_alloc(pk_Id2G2_2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
       class_alloc(pk_Id2G2_4,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_IG2G2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_IG2G2_2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_IG2G2_4,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_IFG2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_IFG2_0b1,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_IFG2_0,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_IFG2_2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_CTR,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_CTR_0,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_CTR_2,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_CTR_4,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(pk_Tree,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_Tree_0_vv,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_Tree_0_vd,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_Tree_0_dd,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_Tree_2_vv,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_Tree_2_vd,pnlpt->k_size*sizeof(double),pnlpt->error_message);
      class_alloc(pk_Tree_4_vv,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(lnk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(lnpk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message);
    class_alloc(ddlnpk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message);

      
    //  double tau_req[pnlpt->z_pk_num];
    //  double deltatau[pnlpt->z_pk_num];
      

      class_alloc(tau_req,pnlpt->z_pk_num * sizeof(double),pnlpt->error_message);
      //class_alloc(deltatau,pnlpt->z_pk_num * sizeof(double),pnlpt->error_message);

      int i_z=0;
      for (i_z=0; i_z<pnlpt->z_pk_num; i_z++) {
          //                        printf("redshift requested %lf\n", pnlpt->z_pk[i]);
          class_call(background_tau_of_z(pba,pnlpt->z_pk[i_z],&tau_req[i_z]),
                     pba->error_message,
                     pnlpt->error_message);
      }
      
      /** Inputing the PT matrices */
          
          char file2openM22[256];
          char file2openM13[256];
          char file2openM22basic[256];
          char file2openM13basic[256];

          if (ppr->nmax_nlpt==128){
          sprintf(file2openM22,"%s/pt_matrices/M22oneline_N128.dat",__CLASSDIR__);
          sprintf(file2openM13,"%s/pt_matrices/M13oneline_N128.dat",__CLASSDIR__);
          sprintf(file2openM22basic,"%s/pt_matrices/M22basiconeline_N128.dat",__CLASSDIR__);
          sprintf(file2openM13basic,"%s/pt_matrices/IFG2oneline_N128.dat",__CLASSDIR__);
        }
        else {
          sprintf(file2openM22,"%s/pt_matrices/M22oneline_N256_packed.dat",__CLASSDIR__);
          sprintf(file2openM13,"%s/pt_matrices/M13oneline_N256.dat",__CLASSDIR__);
          sprintf(file2openM22basic,"%s/pt_matrices/M22basiconeline_N256_packed.dat",__CLASSDIR__);
          sprintf(file2openM13basic,"%s/pt_matrices/IFG2oneline_N256.dat",__CLASSDIR__);
          }
      
      int index_M22=0;
      
      class_alloc(pnlpt->M22_oneline,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)) * sizeof(double),pnlpt->error_message);

      FILE *myFile22;
      
          myFile22 = fopen(file2openM22,"r");
          
      if (myFile22 == NULL){
          printf("Error Reading File M22oneline_....dat\n");
          exit (0);
      }
      
      for (index_M22=0; index_M22 < (ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2); index_M22++){
          fscanf(myFile22, "%lf", &pnlpt->M22_oneline[index_M22]);
      }
      
      fclose(myFile22);
      
      int index_M13=0;
      
      class_alloc(pnlpt->M13_oneline,((ppr->nmax_nlpt+1)*2) * sizeof(double),pnlpt->error_message);
      
      FILE *myFile13;
          
      myFile13 = fopen(file2openM13,"r");
          
      if (myFile13 == NULL){
          printf("Error Reading File M13oneline_....dat\n");
          exit (0);
      }
      
      for (index_M13=0; index_M13 < 2*(ppr->nmax_nlpt+1); index_M13++){
          fscanf(myFile13, "%lf", &pnlpt->M13_oneline[index_M13]);
      }
      
      fclose(myFile13);
      
      class_alloc(pnlpt->M22basic_oneline,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)) * sizeof(double),pnlpt->error_message);
      

      FILE *myFile22_basic;
      
      int index_M22_basic=0;
      
        myFile22_basic = fopen(file2openM22basic,"r");
          
      if (myFile22_basic == NULL){
          printf("Error Reading File M22basiconeline_....dat\n");
          exit (0);
      }
      
      for (index_M22_basic=0; index_M22_basic < (ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2); index_M22_basic++){
          fscanf(myFile22_basic, "%lf", &pnlpt->M22basic_oneline[index_M22_basic]);
      }
      
      fclose(myFile22_basic);
      
      int index_MIFG2=0;
      
      class_alloc(pnlpt->IFG2_oneline,((ppr->nmax_nlpt+1)*2) * sizeof(double),pnlpt->error_message);
      
      FILE *myFile_IFG2;
          
      myFile_IFG2 = fopen(file2openM13basic,"r");
          
      if (myFile_IFG2 == NULL){
          printf("Error Reading File IFG2oneline_....dat\n");
          exit (0);
      }
      
      for (index_MIFG2=0; index_MIFG2 < 2*(ppr->nmax_nlpt+1); index_MIFG2++){
          fscanf(myFile_IFG2, "%lf", &pnlpt->IFG2_oneline[index_MIFG2]);
      }
      
      fclose(myFile_IFG2);
      
      
      
      char file2openGauss[256];
      sprintf(file2openGauss,"%s/pt_matrices/gauss_tab.dat",__CLASSDIR__);
      int index_gauss = 0;
      
      class_alloc(pnlpt->gauss_w,40* sizeof(double),pnlpt->error_message);
      class_alloc(pnlpt->gauss_x,40* sizeof(double),pnlpt->error_message);
      class_alloc(pnlpt->gauss,80* sizeof(double),pnlpt->error_message);
      
      FILE *Gauss_file;
      
      Gauss_file = fopen(file2openGauss,"r");
      
      if (Gauss_file == NULL){
          printf("Error Reading File gauss_tab.dat\n");
          exit (0);
      }
      
      for (index_gauss=0; index_gauss < 80; index_gauss++){
          fscanf(Gauss_file, "%lf", &pnlpt->gauss[index_gauss]);
      }
      
      fclose(Gauss_file);
      
      for (index_gauss=0;index_gauss<40;index_gauss++){
          pnlpt->gauss_x[index_gauss] = pnlpt->gauss[index_gauss];
          pnlpt->gauss_w[index_gauss] = pnlpt->gauss[40+index_gauss];
          //  printf("%lf %lf\n",gauss_x[index_gauss],gauss_w[index_gauss]);
      }
      
      
      
      // This is a place for future optimization !
      
      class_alloc(pnlpt->M22_oneline_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M13_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22basic_oneline_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_0_b1b2_oneline_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_0_b2_oneline_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_0_b1bG2_oneline_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_0_bG2_oneline_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      
      class_alloc(pnlpt->M22_2_b1b2_oneline_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_2_b2_oneline_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_2_b1bG2_oneline_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_2_bG2_oneline_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      
            class_alloc(pnlpt->M22_4_b2_oneline_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
            class_alloc(pnlpt->M22_4_bG2_oneline_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      
      class_alloc(pnlpt->IFG2_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      
      class_alloc(pnlpt->M_IG2G2,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M_Id2,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M_IG2,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M_Id2G2,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(double complex),pnlpt->error_message);
      class_alloc(pnlpt->M13_0_vv_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_0_vv_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M13_0_vd_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_0_vd_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M13_0_dd_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_0_dd_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      
      class_alloc(pnlpt->M13_2_vv_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_2_vv_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      
      class_alloc(pnlpt->M13_2_vd_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_2_vd_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      
      class_alloc(pnlpt->M13_2_dd_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_2_dd_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      
      class_alloc(pnlpt->M13_4_vv_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_4_vv_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      
      class_alloc(pnlpt->M13_4_vd_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_4_vd_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_4_dd_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      
      class_alloc(pnlpt->M13_mu2_dd_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M13_mu2_vd_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M13_mu4_vv_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
       class_alloc(pnlpt->M13_mu4_vd_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M13_mu6_oneline_complex,(ppr->nmax_nlpt+1)*sizeof(complex double),pnlpt->error_message);
      
      class_alloc(pnlpt->M22_oneline_mu2_vd_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_mu2_dd_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_mu4_vd_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_mu4_vv_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_mu4_dd_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_mu6_vv_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_mu6_vd_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      class_alloc(pnlpt->M22_oneline_mu8_complex,((ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2)*sizeof(complex double),pnlpt->error_message);
      
      int count = 0;
      
      for (count=0; count < (ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2; count++){
          pnlpt->M22_oneline_complex[count] = pnlpt->M22_oneline[count] + _Complex_I * pnlpt->M22_oneline[count + (ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2];
      }

    // printf("%i\n",ppr->nmax_nlpt); 

      for (count=0; count < ppr->nmax_nlpt+1; count++){
          pnlpt->M13_oneline_complex[count] = pnlpt->M13_oneline[count] + _Complex_I * pnlpt->M13_oneline[count+ppr->nmax_nlpt+1];
     //     printf("%le\n",creal( pnlpt->M13_oneline_complex[count]));
      }

      
      for (count=0; count < ppr->nmax_nlpt+1; count++){
          pnlpt->IFG2_oneline_complex[count] = pnlpt->IFG2_oneline[count] + _Complex_I * pnlpt->IFG2_oneline[count+ppr->nmax_nlpt+1];
      }
      
      for (count=0; count < (ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2; count++){
          pnlpt->M22basic_oneline_complex[count] = pnlpt->M22basic_oneline[count] + _Complex_I * pnlpt->M22basic_oneline[count + (ppr->nmax_nlpt+1)*(ppr->nmax_nlpt+2)/2];
      }
      
      /* Pt matrices uploaded */
      
      /*It is used by spectra_pk_nl_bias_at_z (for classy) */
class_alloc(pnlpt->ln_pk_nl,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_Id2d2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_Id2d2_2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_Id2d2_4,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      
class_alloc(pnlpt->ln_pk_Id2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_IG2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_Id2G2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_Id2G2_2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_Id2G2_4,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);

      
class_alloc(pnlpt->ln_pk_IG2G2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_IG2G2_2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_IG2G2_4,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);

      
class_alloc(pnlpt->ln_pk_IFG2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_IFG2_0b1,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_IFG2_0,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_IFG2_2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_CTR,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_CTR_0,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_CTR_2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_CTR_4,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_Tree,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_Tree_0_vv,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_Tree_0_vd,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_Tree_0_dd,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_Tree_2_vv,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_Tree_2_vd,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_Tree_4_vv,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_0_vv,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_0_vd,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_0_dd,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);

      class_alloc(pnlpt->ln_pk_2_vv,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_2_vd,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_2_dd,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      
      class_alloc(pnlpt->ln_pk_4_vv,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_4_vd,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_4_dd,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      
class_alloc(pnlpt->ln_pk_0_b1b2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_0_b2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_0_b1bG2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
class_alloc(pnlpt->ln_pk_0_bG2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      
      class_alloc(pnlpt->ln_pk_2_b1b2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_2_b2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_2_b1bG2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_2_bG2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      
      class_alloc(pnlpt->ln_pk_4_b2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_4_bG2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_4_b1b2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);
      class_alloc(pnlpt->ln_pk_4_b1bG2,sizeof(double)*pnlpt->z_pk_num*pnlpt->k_size,pnlpt->error_message);


int index_md;
index_md=pnlpt->index_md_scalars;
int index_ic=0;

// printf("ppt->k_size[index_md]=%i\n",ppt->k_size[index_md]);


       // printf("pnlpt->k[pnlpt->k_size-1]=%le\n",pnlpt->k[pnlpt->k_size-1]);
       // printf("ppt->k[index_md][ppt->k_size[index_md]-1]=%le\n",ppt->k[index_md][ppt->k_size[index_md]-1]);

//int start2=clock();
/*Begin of new part for nonlinear_pt_pk_l*/
if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {
  if (pnlpt->cb == _TRUE_) {
    class_alloc(pnlpt->dd_sources_tp_delta_cdm,
                    ppt->k_size[index_md]*ppt->tau_size*sizeof(double),
                    pnlpt->error_message);
    class_alloc(pnlpt->dd_sources_tp_delta_b,
                    ppt->k_size[index_md]*ppt->tau_size*sizeof(double),
                    pnlpt->error_message);
    class_alloc(pnlpt->sources_tp_delta_cdm,
                    ppt->k_size[index_md]*ppt->tau_size*sizeof(double),
                    pnlpt->error_message);
    class_alloc(pnlpt->sources_tp_delta_b,
                    ppt->k_size[index_md]*ppt->tau_size*sizeof(double),
                    pnlpt->error_message);

    class_call(array_spline_table_columns(ppt->k[index_md],
                                               ppt->k_size[index_md],
                                               ppt->sources[index_md]
                                               [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_cdm],
                                               ppt->tau_size,
                                               pnlpt->dd_sources_tp_delta_cdm,
                                               _SPLINE_EST_DERIV_,
                                               pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

    class_call(array_spline_table_columns(ppt->k[index_md],
                                               ppt->k_size[index_md],
                                               ppt->sources[index_md]
                                               [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_b],
                                               ppt->tau_size,
                                               pnlpt->dd_sources_tp_delta_b,
                                               _SPLINE_EST_DERIV_,
                                               pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

    for (index_tau = pnlpt->tau_size-1; index_tau>=0; index_tau--) {

      for (index_k=0; index_k<pnlpt->k_size; index_k++) {

         if (pnlpt->k[index_k]<=ppt->k[index_md][ppt->k_size[index_md]-1]){ 
        class_call(array_interpolate_spline_one_column(ppt->k[index_md],
                                                ppt->k_size[index_md],
                                                ppt->sources[index_md]
                                                [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_cdm],
                                                pnlpt->tau_size,
                                                index_tau,
                                                pnlpt->dd_sources_tp_delta_cdm,
                                                pnlpt->k[index_k],
                                                &pnlpt->sources_tp_delta_cdm[index_tau*pnlpt->k_size+index_k],
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
  class_call(array_interpolate_spline_one_column(ppt->k[index_md],
                                                ppt->k_size[index_md],
                                                ppt->sources[index_md]
                                                [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_b],
                                                pnlpt->tau_size,
                                                index_tau,
                                                pnlpt->dd_sources_tp_delta_b,
                                                pnlpt->k[index_k],
                                                &pnlpt->sources_tp_delta_b[index_tau*pnlpt->k_size+index_k],
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
}// condition for the new "if" 
/* Very crude padding, but the physics should not depend on these scales anlyway */
/* This is done only for the non-zero spatial curvature, for which the kmax from perturbations does not coincude with the kmax from the spectra */
else{
  pnlpt->sources_tp_delta_b[index_tau*pnlpt->k_size+index_k] = pnlpt->sources_tp_delta_b[index_tau*pnlpt->k_size+index_k-1];
  pnlpt->sources_tp_delta_cdm[index_tau*pnlpt->k_size+index_k] = pnlpt->sources_tp_delta_cdm[index_tau*pnlpt->k_size+index_k-1];
 }
  // printf(" pnlpt->sources_tp_delta_b[index_tau*pnlpt->k_size+index_k] = %le\n", pnlpt->sources_tp_delta_b[index_tau*pnlpt->k_size+index_k]);
      }
    }
  } else {
    class_alloc(pnlpt->dd_sources_tp_delta_m,
                    ppt->k_size[index_md]*ppt->tau_size*sizeof(double),
                    pnlpt->error_message);
    class_alloc(pnlpt->sources_tp_delta_m,
                    ppt->k_size[index_md]*ppt->tau_size*sizeof(double),
                    pnlpt->error_message);

    class_call(array_spline_table_columns(ppt->k[index_md],
                                               ppt->k_size[index_md],
                                               ppt->sources[index_md]
                                               [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_m],
                                               ppt->tau_size,
                                               pnlpt->dd_sources_tp_delta_m,
                                               _SPLINE_EST_DERIV_,
                                               pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
    for (index_tau = pnlpt->tau_size-1; index_tau>=0; index_tau--) {
      for (index_k=0; index_k<pnlpt->k_size; index_k++) {

                 if (pnlpt->k[index_k]<=ppt->k[index_md][ppt->k_size[index_md]-1]){ 
      class_call(array_interpolate_spline_one_column(ppt->k[index_md],
                                                ppt->k_size[index_md],
                                                ppt->sources[index_md]
                                                [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_m],
                                                pnlpt->tau_size,
                                                index_tau,
                                                pnlpt->dd_sources_tp_delta_m,
                                                pnlpt->k[index_k],
                                                &pnlpt->sources_tp_delta_m[index_tau*pnlpt->k_size+index_k],
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
      }// condition for the new "if" 
/* Very crude padding, but the physics should not depend on these scales anlyway */
/* This is done only for curvature, for which the kmax from perturbations does not coincude with the kmax from the spectra */
else{
  pnlpt->sources_tp_delta_m[index_tau*pnlpt->k_size+index_k] = pnlpt->sources_tp_delta_m[index_tau*pnlpt->k_size+index_k-1];
 }

      }
    }
  }
/*End cb*/
}

/*End of new part for nonlinear_pt_pk_l*/
//int end2=clock();
//printf("new part is run over %d ms\n",end2-start2);

      /*determination of full tables of P_L(k,tau) lnP_L(k,tau) and ddP_L(k,tau)*/
double *lnpk_l_full;
double *pk_l_full;
double *ddlnpk_l_full;
/*double *ddddlnpk_l_full;*/
   class_alloc(lnpk_l_full,sizeof(double)*pnlpt->tau_size*pnlpt->k_size,pnlpt->error_message); 
   class_alloc(pk_l_full,sizeof(double)*pnlpt->tau_size*pnlpt->k_size,pnlpt->error_message);
class_alloc(ddlnpk_l_full,sizeof(double)*pnlpt->tau_size*pnlpt->k_size,pnlpt->error_message);
/*class_alloc(ddddlnpk_l_full,sizeof(double)*pnlpt->tau_size*pnlpt->k_size,pnlpt->error_message);*/

/*begin cycle over tau*/
for (index_tau = pnlpt->tau_size-1; index_tau>=0; index_tau--) {

      /* get P_L(k) at this time */
          
          class_call(nonlinear_pt_pk_l(pba,ppt,ppm,pnlpt,index_tau,pk_l,lnk_l,lnpk_l,ddlnpk_l),
                     pnlpt->error_message,
                     pnlpt->error_message);

      /* get P_L(k,tau) lnP_L(k,tau) and ddP_L(k,tau) */

    for (index_k=0; index_k<pnlpt->k_size; index_k++) {
    lnpk_l_full[index_tau * pnlpt->k_size + index_k]=lnpk_l[index_k];
    pk_l_full[index_tau * pnlpt->k_size + index_k]=pk_l[index_k];
    ddlnpk_l_full[index_tau * pnlpt->k_size + index_k]=ddlnpk_l[index_k];
  //  pnlpt->nl_corr_density[index_tau * pnlpt->k_size + index_k] = 1.;
    }
}
/*end cycle over tau*/
//int end1=clock();
//printf("tau cycle is run over %d ms\n",end1-start1);

//int end1=clock();
//printf("tau cycle is run over %d ms\n",end1-start1);

/* interpolation of P_L(k,tau) lnP_L(k,tau) and ddP_L(k,tau) at tau_req[i_z]*/
double *ln_pk_l_at_z_req;
double *pk_l_at_z_req;

double Dref = 0;

class_alloc(ln_pk_l_at_z_req,sizeof(double)*pnlpt->k_size,pnlpt->error_message);
class_alloc(pk_l_at_z_req,sizeof(double)*pnlpt->k_size,pnlpt->error_message);

/* Retriveing the logarithmic growth factor f for RSD */

class_alloc(pnlpt->growthf,sizeof(double)*pnlpt->z_pk_num,pnlpt->error_message);
      class_alloc(pnlpt->hratio_array,sizeof(double)*pnlpt->z_pk_num,pnlpt->error_message);
      class_alloc(pnlpt->Dratio_array,sizeof(double)*pnlpt->z_pk_num,pnlpt->error_message);
      

if (pnlpt->rsd == rsd_yes) {

double *pvecbackf;
int last_indexf;
class_alloc(pvecbackf,pba->bg_size*sizeof(double),pnlpt->error_message);        

// these things are needed for Alcock-Pazynsky stuff
int j;
int Nz = 2000;
//double Omfid = 0.307115;
//double Omfid = 0.3;
double Omfid = pnlpt->OmfidAP;
double dz = 0;
double kmsMpc = 3.33564095198145e-6;
//double Omtrue = pba->Omega0_ncdm_tot + pba->Omega0_cdm + pba->Omega0_b;
    
    if (pnlpt->nonlinear_pt_verbose > 0)
    if (pnlpt->AP_effect == AP_effect_yes){
        printf("Computing the Alcock-Paczynski effect for fiducial cosmology with Om=%lf\n",Omfid);
    }
    else {
        printf("No Alcock-Paczynski effect.\n");
    }
    
for (i_z=0; i_z<pnlpt->z_pk_num; i_z++) {
    
    double Da = 0;
    double Dfid = 0;
    double hnew = 0;
    double hfid = 0;
    
   
class_call(background_at_tau(pba,tau_req[i_z],pba->long_info,pba->inter_normal,&last_indexf,pvecbackf),
    pnlpt->error_message,
    pnlpt->error_message);

 pnlpt->growthf[i_z] =  pvecbackf[pba->index_bg_f];
 Dref = pvecbackf[pba->index_bg_D];

 //printf("Dref=%lf\n",Dref);
    
if (pnlpt->AP_effect == AP_effect_yes){
    
   // printf("z_pk[i_z]=%lf\n",pnlpt->z_pk[i_z]);
    
if (pnlpt->z_pk[i_z]== 0.)
    {
        pnlpt->hratio_array[i_z] = 1.;
        pnlpt->Dratio_array[i_z] = 1.;
        hnew = 1.;
    }
    else {
   //     hnew = pow((Omtrue*pow((1.+pnlpt->z_pk[i_z]),3.) + (1. - Omtrue)+(pba->Omega0_g)*pow((1.+pnlpt->z_pk[i_z]),4.)),0.5);
        hfid = pow((Omfid*pow((1.+pnlpt->z_pk[i_z]),3.) + (1. - Omfid)+(pba->Omega0_g)*pow((1.+pnlpt->z_pk[i_z]),4.)),0.5);
   //     printf("Hclass=%le\n",pvecbackf[pba->index_bg_H]/kmsMpc/100/pba->h);
  //      printf("Hmy = %le\n",hnew);
        hnew = pvecbackf[pba->index_bg_H]/kmsMpc/100/pba->h;
        pnlpt->hratio_array[i_z] = hnew/hfid;
        
       // pnlpt->hratio_array[i_z] = pow(((pba->Omega0_cdm+pba->Omega0_b)*pow((1.+pnlpt->z_pk[i_z]),3.) + (1. - pba->Omega0_cdm - pba->Omega0_b)+(pba->Omega0_g)*pow((1.+pnlpt->z_pk[i_z]),4.)),0.5)/pow((Omfid*pow((1.+pnlpt->z_pk[i_z]),3.) + (1. - Omfid)+(pba->Omega0_g)*pow((1.+pnlpt->z_pk[i_z]),4.)),0.5);
    dz = pnlpt->z_pk[i_z]/(1.*Nz-1.);
    
for (j=1; j<Nz; j ++) {
 //   Da = Da + dz*(1./pow((Omtrue*pow((1.+ dz*j),3.) + (1. -Omtrue)+(pba->Omega0_g)*pow((1.+pnlpt->z_pk[i_z]),4.)),0.5)+1./pow((Omtrue*pow((1.+ dz*(j-1)),3.) + (1. - Omtrue)+(pba->Omega0_g)*pow((1.+pnlpt->z_pk[i_z]),4.)),0.5))/2.;
    Dfid = Dfid + dz*(1./pow((Omfid*pow((1.+ dz*j),3.) + (1. - Omfid)+(pba->Omega0_g)*pow((1.+pnlpt->z_pk[i_z]),4.)),0.5)+1./pow((Omfid*pow((1.+dz*(j-1)),3.) + (1. - Omfid)+(pba->Omega0_g)*pow((1.+pnlpt->z_pk[i_z]),4.)),0.5))/2.;
                    }
 //   printf("Dclass=%le\n",pvecbackf[pba->index_bg_ang_distance]*kmsMpc*100*pba->h*(1.+pnlpt->z_pk[i_z]));
//    printf("Dmy = %le\n",Da);
    Da = pvecbackf[pba->index_bg_ang_distance]*kmsMpc*100*pba->h*(1.+pnlpt->z_pk[i_z]);
    pnlpt->Dratio_array[i_z] = Da/Dfid;
    //    printf(" h/hfid=%le\n", hnew/hfid);
    //    printf(" Da/Dfid=%le\n", Da/Dfid);
    }
}
else {
    pnlpt->hratio_array[i_z] = 1.;
    pnlpt->Dratio_array[i_z] = 1.;
}
 //   printf("pba->Omega0_ncdm_tot = %lf\n",pba->Omega0_ncdm_tot);
 //   printf("pba->M_ncdm[0] = %lf\n",pba->M_ncdm[0]);
  //  printf("pba->m_ncdm_in_eV[0] = %lf\n",pba->m_ncdm_in_eV[0]);
    
//    printf("hratio = %lf\n",hratio);
//    printf("Dfid = %lf\n",Dfid);

    // Check that Da and H(z) are computed correctly
    
//printf("Om=%lf\n",pba->Omega0_cdm+pba->Omega0_b);
//printf("pvecbackf[pba->index_bg_a]=%lf\n",pvecbackf[pba->index_bg_a]);
//printf("z_pk[i_z]=%lf\n",pnlpt->z_pk[i_z]);
//printf("pvecbackf[pba->index_bg_H]=%le\n",pvecbackf[pba->index_bg_H]/kmsMpc);
//printf("pvecbackf[pba->index_bg_ang_distance]=%le\n",pvecbackf[pba->index_bg_ang_distance]*kmsMpc);
//printf("D = %le\n",Da/100/pba->h/(1.+pnlpt->z_pk[i_z]));
//printf("H = %le\n",hnew*100*pba->h);
    
  //  printf("z_pk[i_z]=%lf\n",pnlpt->z_pk[i_z]);
   // printf("f[i_z]=%lf\n",pvecbackf[pba->index_bg_f]);

          }

free(pvecbackf);
}
      
else {

  double *pvecbackf;
  int last_indexf = 0;
  class_alloc(pvecbackf,pba->bg_size*sizeof(double),pnlpt->error_message);

    for (i_z=0; i_z<pnlpt->z_pk_num; i_z++) {
    pnlpt->growthf[i_z] = 1.;
         pnlpt->hratio_array[i_z] = 1.;
         pnlpt->Dratio_array[i_z] = 1.;
       //  printf("No RSD computed\n");


for (i_z=0; i_z<pnlpt->z_pk_num; i_z++) {
class_call(background_at_tau(pba,tau_req[i_z],pba->long_info,pba->inter_normal,&last_indexf,pvecbackf),
    pnlpt->error_message,
    pnlpt->error_message);
    Dref = pvecbackf[pba->index_bg_D];
  }


    }
    free(pvecbackf);
}

      
/* end of RSD specification */
      
      //ln_pk_l_at_z_req is the array of the linear PS

      last_index = 0;
for (i_z=0; i_z<pnlpt->z_pk_num; i_z++) {  

    class_call(array_interpolate_spline(pnlpt->ln_tau,
                                                pnlpt->tau_size,
                                                lnpk_l_full,
                                                ddlnpk_l_full,
                                                pnlpt->k_size,
                                                log(tau_req[i_z]),
                                                &last_index,
                                                ln_pk_l_at_z_req,
                                                pnlpt->k_size,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);

    for (index_k=0; index_k<pnlpt->k_size; index_k++) {        
    ln_pk_l_at_z_req[index_k] = ln_pk_l_at_z_req[index_k];
    pk_l_at_z_req[index_k]=exp(ln_pk_l_at_z_req[index_k]);
  //      printf("%le %le\n",khere,exp(ln_pk_l_at_z_req[index_k]));
    }

       /* get P_NL(k) at tau_req */
       
        if (print_warning == _FALSE_) {
int start=clock();
          class_call(nonlinear_pt_loop(ppr,
                                     pba,
                                     ppm,
                                     pth,
                                     pnlpt,
                                     tau_req[i_z],
                                     pnlpt->growthf[i_z],
                                       pnlpt->hratio_array[i_z],
                                       pnlpt->Dratio_array[i_z],
                                     pk_l_at_z_req,
                                     pk_l_0_vv,
                                     pk_l_0_vd,
                                     pk_l_0_dd,
                                     pk_l_2_vv,
                                     pk_l_2_vd,
                                     pk_l_2_dd,
                                     pk_l_4_vv,
                                     pk_l_4_vd,
                                     pk_l_4_dd,
                                     pk_nl,
                                     pk_Id2d2,
                                       pk_Id2d2_2,
                                       pk_Id2d2_4,
                                     pk_l_0_b1b2,
                                     pk_l_0_b2,
                                     pk_l_0_b1bG2,
                                     pk_l_0_bG2,
                                       pk_l_2_b1b2,
                                       pk_l_2_b2,
                                       pk_l_2_b1bG2,
                                       pk_l_2_bG2,
                                       pk_l_4_b2,
                                       pk_l_4_bG2,
                                       pk_l_4_b1b2,
                                       pk_l_4_b1bG2,
                                     pk_Id2,
                                     pk_IG2,
                                     pk_Id2G2,
                                       pk_Id2G2_2,
                                       pk_Id2G2_4,
                                     pk_IG2G2,
                                       pk_IG2G2_2,
                                       pk_IG2G2_4,
                                     pk_IFG2,
                                       pk_IFG2_0b1,
                                       pk_IFG2_0,
                                       pk_IFG2_2,
                                     pk_CTR,
                                       pk_CTR_0,
                                       pk_CTR_2,
                                       pk_CTR_4,
				                     pk_Tree,
                                       pk_Tree_0_vv,
                                       pk_Tree_0_vd,
                                       pk_Tree_0_dd,
                                       pk_Tree_2_vv,
                                       pk_Tree_2_vd,
                                       pk_Tree_4_vv,
                                     lnk_l,
                                     ln_pk_l_at_z_req
                                  ),
                   pnlpt->error_message,
                   pnlpt->error_message);
int end=clock();
            

if (pnlpt->nonlinear_pt_verbose > 0)
printf("Module nonlinear_pt_loop takes %d musec.\n",end-start);

           for (index_k=0; index_k<pnlpt->k_size; index_k++) {
                pnlpt->ln_pk_nl[i_z*pnlpt->k_size+index_k] = log(pk_nl[index_k]);
                pnlpt->ln_pk_Id2d2[i_z*pnlpt->k_size+index_k] = log(pk_Id2d2[index_k]);
               pnlpt->ln_pk_Id2d2_2[i_z*pnlpt->k_size+index_k] = log(pk_Id2d2_2[index_k]);
               pnlpt->ln_pk_Id2d2_4[i_z*pnlpt->k_size+index_k] = log(pk_Id2d2_4[index_k]);
               
                pnlpt->ln_pk_Id2[i_z*pnlpt->k_size+index_k] = log(pk_Id2[index_k]);
                pnlpt->ln_pk_IG2[i_z*pnlpt->k_size+index_k] = log(pk_IG2[index_k]);
                pnlpt->ln_pk_Id2G2[i_z*pnlpt->k_size+index_k] = log(pk_Id2G2[index_k]);
               pnlpt->ln_pk_Id2G2_2[i_z*pnlpt->k_size+index_k] = log(pk_Id2G2_2[index_k]);
               pnlpt->ln_pk_Id2G2_4[i_z*pnlpt->k_size+index_k] = log(pk_Id2G2_4[index_k]);
               
                pnlpt->ln_pk_IG2G2[i_z*pnlpt->k_size+index_k] = log(pk_IG2G2[index_k]);
               pnlpt->ln_pk_IG2G2_2[i_z*pnlpt->k_size+index_k] = log(pk_IG2G2_2[index_k]);
               pnlpt->ln_pk_IG2G2_4[i_z*pnlpt->k_size+index_k] = log(pk_IG2G2_4[index_k]);
               
                pnlpt->ln_pk_IFG2[i_z*pnlpt->k_size+index_k] = log(pk_IFG2[index_k]);
               
               pnlpt->ln_pk_IFG2_0[i_z*pnlpt->k_size+index_k] = log(pk_IFG2_0[index_k]);
               pnlpt->ln_pk_IFG2_0b1[i_z*pnlpt->k_size+index_k] = log(pk_IFG2_0b1[index_k]);
               pnlpt->ln_pk_IFG2_2[i_z*pnlpt->k_size+index_k] = log(pk_IFG2_2[index_k]);
               
                pnlpt->ln_pk_CTR[i_z*pnlpt->k_size+index_k] = log(pk_CTR[index_k]);
               pnlpt->ln_pk_CTR_0[i_z*pnlpt->k_size+index_k] = log(pk_CTR_0[index_k]);
               
               pnlpt->ln_pk_CTR_2[i_z*pnlpt->k_size+index_k] = log(pk_CTR_2[index_k]);
               pnlpt->ln_pk_CTR_4[i_z*pnlpt->k_size+index_k] = log(pk_CTR_4[index_k]);
               
		        pnlpt->ln_pk_Tree[i_z*pnlpt->k_size+index_k] = log(pk_Tree[index_k]);
               pnlpt->ln_pk_Tree_0_vv[i_z*pnlpt->k_size+index_k] = log(pk_Tree_0_vv[index_k]);
               pnlpt->ln_pk_Tree_0_vd[i_z*pnlpt->k_size+index_k] = log(pk_Tree_0_vd[index_k]);
               pnlpt->ln_pk_Tree_0_dd[i_z*pnlpt->k_size+index_k] = log(pk_Tree_0_dd[index_k]);
               pnlpt->ln_pk_Tree_2_vv[i_z*pnlpt->k_size+index_k] = log(pk_Tree_2_vv[index_k]);
               pnlpt->ln_pk_Tree_2_vd[i_z*pnlpt->k_size+index_k] = log(pk_Tree_2_vd[index_k]);
               pnlpt->ln_pk_Tree_4_vv[i_z*pnlpt->k_size+index_k] = log(pk_Tree_4_vv[index_k]);
                pnlpt->ln_pk_0_vv[i_z*pnlpt->k_size+index_k] = log(pk_l_0_vv[index_k]);
                pnlpt->ln_pk_0_vd[i_z*pnlpt->k_size+index_k] = log(pk_l_0_vd[index_k]);
                pnlpt->ln_pk_0_dd[i_z*pnlpt->k_size+index_k] = log(pk_l_0_dd[index_k]);
                pnlpt->ln_pk_0_b1b2[i_z*pnlpt->k_size+index_k] = log(pk_l_0_b1b2[index_k]);
               pnlpt->ln_pk_0_b2[i_z*pnlpt->k_size+index_k] = log(pk_l_0_b2[index_k]);
               pnlpt->ln_pk_0_b1bG2[i_z*pnlpt->k_size+index_k] = log(pk_l_0_b1bG2[index_k]);
               pnlpt->ln_pk_0_bG2[i_z*pnlpt->k_size+index_k] = log(pk_l_0_bG2[index_k]);
               
               pnlpt->ln_pk_2_b1b2[i_z*pnlpt->k_size+index_k] = log(pk_l_2_b1b2[index_k]);
               pnlpt->ln_pk_2_b2[i_z*pnlpt->k_size+index_k] = log(pk_l_2_b2[index_k]);
               pnlpt->ln_pk_2_b1bG2[i_z*pnlpt->k_size+index_k] = log(pk_l_2_b1bG2[index_k]);
               pnlpt->ln_pk_2_bG2[i_z*pnlpt->k_size+index_k] = log(pk_l_2_bG2[index_k]);
               
               pnlpt->ln_pk_4_b2[i_z*pnlpt->k_size+index_k] = log(pk_l_4_b2[index_k]);
               pnlpt->ln_pk_4_bG2[i_z*pnlpt->k_size+index_k] = log(pk_l_4_bG2[index_k]);
               pnlpt->ln_pk_4_b1b2[i_z*pnlpt->k_size+index_k] = log(pk_l_4_b1b2[index_k]);
               pnlpt->ln_pk_4_b1bG2[i_z*pnlpt->k_size+index_k] = log(pk_l_4_b1bG2[index_k]);
               
               pnlpt->ln_pk_2_vv[i_z*pnlpt->k_size+index_k] = log(pk_l_2_vv[index_k]);
               pnlpt->ln_pk_2_vd[i_z*pnlpt->k_size+index_k] = log(pk_l_2_vd[index_k]);
               pnlpt->ln_pk_2_dd[i_z*pnlpt->k_size+index_k] = log(pk_l_2_dd[index_k]);
               
               pnlpt->ln_pk_4_vv[i_z*pnlpt->k_size+index_k] = log(pk_l_4_vv[index_k]);
               pnlpt->ln_pk_4_vd[i_z*pnlpt->k_size+index_k] = log(pk_l_4_vd[index_k]);
               pnlpt->ln_pk_4_dd[i_z*pnlpt->k_size+index_k] = log(pk_l_4_dd[index_k]);
            
            }
      }
}      

 //printf("Dlast=%f \n",Dref);

  /* This is required for lensing with one-loop ! */

/*Begin lensing module*/
if (ppt->has_cl_cmb_lensing_potential == _TRUE_){
/*Lensing turn on*/
double Dplus = 0;
double *pvecbackD;
int last_indexD = 0;
class_alloc(pvecbackD,pba->bg_size*sizeof(double),pnlpt->error_message);

double * dd_pk_Tree;
double * dd_pk_nl;
double * dd_pk_l_at_z_req;
double * pk_Tree_int;
double * pk_nl_int;
double * pk_l_at_z_req_int;
// double * pk_ctr_int;

last_index = 0;

if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {
class_alloc(dd_pk_Tree, pnlpt->k_size*sizeof(double), pnlpt->error_message);
class_alloc(dd_pk_nl, pnlpt->k_size*sizeof(double), pnlpt->error_message);
class_alloc(dd_pk_l_at_z_req, pnlpt->k_size*sizeof(double), pnlpt->error_message);
class_alloc(pk_Tree_int, ppt->k_size[pnlpt->index_md_scalars]*sizeof(double), pnlpt->error_message);
class_alloc(pk_nl_int, ppt->k_size[pnlpt->index_md_scalars]*sizeof(double), pnlpt->error_message);
class_alloc(pk_l_at_z_req_int, ppt->k_size[pnlpt->index_md_scalars]*sizeof(double), pnlpt->error_message);
// class_alloc(pk_ctr_int, ppt->k_size[pnlpt->index_md_scalars]*sizeof(double), pnlpt->error_message);

    class_call(array_spline_table_columns(pnlpt->k,
                                              pnlpt->k_size,
                                              pk_Tree,
                                              1,
                                              dd_pk_Tree,
                                              _SPLINE_EST_DERIV_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
    class_call(array_spline_table_columns(pnlpt->k,
                                              pnlpt->k_size,
                                              pk_nl,
                                              1,
                                              dd_pk_nl,
                                              _SPLINE_EST_DERIV_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
    class_call(array_spline_table_columns(pnlpt->k,
                                              pnlpt->k_size,
                                              pk_l_at_z_req,
                                              1,
                                              dd_pk_l_at_z_req,
                                              _SPLINE_EST_DERIV_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
//printf("ppt->k_max=%g  pnlpt->k_max=%g\n",ppt->k[pnlpt->index_md_scalars][ppt->k_size[pnlpt->index_md_scalars]-1],pnlpt->k[pnlpt->k_size-1]);
    for (index_k=0; index_k < ppt->k_size[pnlpt->index_md_scalars]; index_k++){
      if (ppt->k[pnlpt->index_md_scalars][index_k]<pnlpt->k[pnlpt->k_size-1]) {
        class_call(array_interpolate_spline(pnlpt->k,
                                                    pnlpt->k_size,
                                                    pk_Tree,
                                                    dd_pk_Tree,
                                                    1,
                                                    ppt->k[pnlpt->index_md_scalars][index_k],
                                                    &last_index,
                                                    &pk_Tree_int[index_k],
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);
        class_call(array_interpolate_spline(pnlpt->k,
                                                    pnlpt->k_size,
                                                    pk_nl,
                                                    dd_pk_nl,
                                                    1,
                                                    ppt->k[pnlpt->index_md_scalars][index_k],
                                                    &last_index,
                                                    &pk_nl_int[index_k],
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);
        class_call(array_interpolate_spline(pnlpt->k,
                                                    pnlpt->k_size,
                                                    pk_l_at_z_req,
                                                    dd_pk_l_at_z_req,
                                                    1,
                                                    ppt->k[pnlpt->index_md_scalars][index_k],
                                                    &last_index,
                                                    &pk_l_at_z_req_int[index_k],
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);
        ppt->k[pnlpt->index_md_scalars][index_k];
        // pk_ctr_int[index_k] = 0.;
        // pk_ctr_int[index_k] = 2*pk_Tree_int[index_k]*ppt->k[pnlpt->index_md_scalars][index_k]*ppt->k[pnlpt->index_md_scalars][index_k]/(1.+ppt->k[pnlpt->index_md_scalars][index_k]*ppt->k[pnlpt->index_md_scalars][index_k]);

      } else {
        pk_Tree_int[index_k] = pk_Tree_int[index_k-1];
        pk_nl_int[index_k] = pk_nl_int[index_k-1];
        pk_l_at_z_req_int[index_k] = pk_l_at_z_req_int[index_k-1];
        // pk_ctr_int[index_k] = 0.;
        // pk_ctr_int[index_k] = 2*pk_Tree_int[index_k]*ppt->k[pnlpt->index_md_scalars][index_k-1]*ppt->k[pnlpt->index_md_scalars][index_k-1]/(1.+ppt->k[pnlpt->index_md_scalars][index_k-1]*ppt->k[pnlpt->index_md_scalars][index_k-1]);
      }
    }
}
        

for (index_tau = pnlpt->tau_size-1; index_tau>=0; index_tau--) {

/*
class_call(nonlinear_pt_pk_l(pba,ppt,ppm,pnlpt,index_tau,pk_l,lnk_l,lnpk_l,ddlnpk_l),
                     pnlpt->error_message,
                     pnlpt->error_message);
                     */

class_call(background_at_tau(pba,pnlpt->tau[index_tau],pba->long_info,pba->inter_normal,&last_indexD,pvecbackD),pnlpt->error_message,pnlpt->error_message);

    Dplus =  pvecbackD[pba->index_bg_D];
    //printf("D+=%f \n",Dplus);
    if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {
    for (index_k=0; index_k<ppt->k_size[pnlpt->index_md_scalars]; index_k++) 
        //pnlpt->nl_corr_density[index_tau * ppt->k_size[pnlpt->index_md_scalars] + index_k] = 1.;
        pnlpt->nl_corr_density[index_tau * ppt->k_size[pnlpt->index_md_scalars] + index_k] = sqrt((pk_Tree_int[index_k]+Dplus*Dplus*(pk_nl_int[index_k]-5000. 
          // -pk_ctr_int[index_k]
          // -2.*pow(ppt->k[pnlpt->index_md_scalars][index_k],2.)*pk_Tree_int[index_k]/(1.+pow(ppt->k[pnlpt->index_md_scalars][index_k],2.)) 
          )/Dref/Dref)/pk_l_at_z_req_int[index_k]);
    } else {
    for (index_k=0; index_k<ppt->k_size[pnlpt->index_md_scalars]; index_k++) 
        pnlpt->nl_corr_density[index_tau * ppt->k_size[pnlpt->index_md_scalars] + index_k] = sqrt((pk_Tree[index_k]+Dplus*Dplus*(pk_nl[index_k]-5000. 
          // -pk_ctr_int[index_k]
       // -2.*pow(ppt->k[pnlpt->index_md_scalars][index_k],2.)*pk_Tree_int[index_k]/(1.+pow(ppt->k[pnlpt->index_md_scalars][index_k],2.)) 
          )/Dref/Dref)/pk_l_at_z_req[index_k]);
    }
    //printf("index_tau=%i \n",index_tau);
    //printf("index_k=%i \n",index_k);
    //printf("pnlpt->nl_corr_density[index_tau * pnlpt->k_size + index_k]=%f \n",pnlpt->nl_corr_density[index_tau * pnlpt->k_size + index_k]);
    
}
free(pvecbackD);

if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {
free(dd_pk_Tree);
free(dd_pk_nl);
free(dd_pk_l_at_z_req);

free(pk_Tree_int);
// free(pk_ctr_int);
free(pk_nl_int);
free(pk_l_at_z_req_int);
}

} else {
/*Lensing turn off*/
  for (index_tau = pnlpt->tau_size-1; index_tau>=0; index_tau--) {
    for (index_k=0; index_k<ppt->k_size[pnlpt->index_md_scalars]; index_k++) {
      pnlpt->nl_corr_density[index_tau * ppt->k_size[pnlpt->index_md_scalars] + index_k] = 1.;
    }
  }
}
/*End lensing module*/

    free(pk_l);
    free(pk_l_0_vv);
    free(pk_l_0_vd);
    free(pk_l_0_dd);
    free(pk_l_2_vv);
    free(pk_l_2_vd);
    free(pk_l_2_dd);
    free(pk_l_4_vv);
    free(pk_l_4_vd);
    free(pk_l_4_dd);
      
    free(pk_l_0_b1b2);
    free(pk_l_0_b2);
    free(pk_l_0_b1bG2);
    free(pk_l_0_bG2);
      
      free(pk_l_2_b1b2);
      free(pk_l_2_b2);
      free(pk_l_2_b1bG2);
      free(pk_l_2_bG2);
      
      free(pk_l_4_b2);
      free(pk_l_4_bG2);
      free(pk_l_4_b1b2);
      free(pk_l_4_b1bG2);
      
    free(pk_nl);
    free(pk_Id2d2);
      free(pk_Id2d2_2);
      free(pk_Id2d2_4);
    free(pk_IG2);
    free(pk_Id2G2);
    free(pk_IG2G2);
      
      free(pk_Id2G2_2);
      free(pk_IG2G2_2);
      
      free(pk_Id2G2_4);
      free(pk_IG2G2_4);
      
    free(pk_IFG2);
      free(pk_IFG2_0);
      free(pk_IFG2_0b1);
      free(pk_IFG2_2);
      
    free(pk_Id2);
    free(pk_CTR);
      free(pk_CTR_0);
      free(pk_CTR_2);
      free(pk_CTR_4);
    free(pk_Tree);
      free(pk_Tree_0_vv);
      free(pk_Tree_0_vd);
      free(pk_Tree_0_dd);
      free(pk_Tree_2_vv);
      free(pk_Tree_2_vd);
      free(pk_Tree_4_vv);
    free(lnk_l);
    free(lnpk_l);
    free(ddlnpk_l);
    free(tau_req);
      
    free(lnpk_l_full);
    free(pk_l_full);
    free(ddlnpk_l_full);
    free(ln_pk_l_at_z_req);
    free(pk_l_at_z_req);
     

if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {
  if (pnlpt->cb == _TRUE_) {
    free(pnlpt->dd_sources_tp_delta_cdm);
    free(pnlpt->sources_tp_delta_cdm);
    free(pnlpt->dd_sources_tp_delta_b);
    free(pnlpt->sources_tp_delta_b);
  } else {
    free(pnlpt->dd_sources_tp_delta_m);
    free(pnlpt->sources_tp_delta_m);
  }
}

     if (pnlpt->nonlinear_pt_verbose > 0)
     printf(" 'nonlinear_pt_init' module executed successfully\n");
  }
  
  else {
    class_stop(pnlpt->error_message,
               "Your non-linear method variable is set to %d, out of the range defined in nonlinear_pt.h",pnlpt->method);
  }
    
    //end of the loop over time before 'else' above

  return _SUCCESS_;
    
}


//end of non-linear init

int nonlinear_pt_free(
                   struct nonlinear_pt *pnlpt
                   ) {

  if (pnlpt->method > nlpt_none) {

    if (pnlpt->method == nlpt_spt) {
      /* free here */
        
        free(pnlpt->M13_oneline);
        free(pnlpt->M22_oneline);

        free(pnlpt->M22basic_oneline);
        free(pnlpt->IFG2_oneline);
        free(pnlpt->k);
        free(pnlpt->ln_k);
        free(pnlpt->tau);
        free(pnlpt->ln_tau);

        free(pnlpt->M22_oneline_0_vv_complex);
        free(pnlpt->M13_0_vv_oneline_complex);
        free(pnlpt->M22_oneline_0_vd_complex);
        free(pnlpt->M13_0_vd_oneline_complex);
        free(pnlpt->M22_oneline_0_dd_complex);
        free(pnlpt->M13_0_dd_oneline_complex);
        
        free(pnlpt->M22_oneline_2_vv_complex);
        free(pnlpt->M13_2_vv_oneline_complex);
        free(pnlpt->M22_oneline_2_vd_complex);
        free(pnlpt->M13_2_vd_oneline_complex);
        free(pnlpt->M22_oneline_2_dd_complex);
        free(pnlpt->M13_2_dd_oneline_complex);
        
        free(pnlpt->M22_oneline_4_vv_complex);
        free(pnlpt->M13_4_vv_oneline_complex);
        free(pnlpt->M22_oneline_4_vd_complex);
        free(pnlpt->M13_4_vd_oneline_complex);
        free(pnlpt->M22_oneline_4_dd_complex);
        
        free(pnlpt->M22_0_b1b2_oneline_complex);
        free(pnlpt->M22_0_b2_oneline_complex);
        free(pnlpt->M22_0_b1bG2_oneline_complex);
        free(pnlpt->M22_0_bG2_oneline_complex);
        
        free(pnlpt->M22_2_b1b2_oneline_complex);
        free(pnlpt->M22_2_b2_oneline_complex);
        free(pnlpt->M22_2_b1bG2_oneline_complex);
        free(pnlpt->M22_2_bG2_oneline_complex);
        
        free(pnlpt->M22_4_b2_oneline_complex);
        free(pnlpt->M22_4_bG2_oneline_complex);
        
        free(pnlpt->M13_oneline_complex);
        free(pnlpt->M22_oneline_complex);
        free(pnlpt->M22basic_oneline_complex);
        free(pnlpt->IFG2_oneline_complex);

        
        free(pnlpt->M13_mu2_vd_oneline_complex);
        free(pnlpt->M13_mu2_dd_oneline_complex);
        free(pnlpt->M13_mu4_vv_oneline_complex);
        free(pnlpt->M13_mu4_vd_oneline_complex);
        free(pnlpt->M13_mu6_oneline_complex);
        
        free(pnlpt->M22_oneline_mu2_vd_complex);
        free(pnlpt->M22_oneline_mu2_dd_complex);
        free(pnlpt->M22_oneline_mu4_vv_complex);
        free(pnlpt->M22_oneline_mu4_vd_complex);
        free(pnlpt->M22_oneline_mu4_dd_complex);
        free(pnlpt->M22_oneline_mu6_vv_complex);
        free(pnlpt->M22_oneline_mu6_vd_complex);
        free(pnlpt->M22_oneline_mu8_complex);
        
        free(pnlpt->M_Id2);
        free(pnlpt->M_IG2);
        free(pnlpt->M_Id2G2);
        free(pnlpt->M_IG2G2);
        
        free(pnlpt->growthf);
        free(pnlpt->hratio_array);
        free(pnlpt->Dratio_array);

        free(pnlpt->nl_corr_density);
        
        free(pnlpt->ln_pk_nl);
        free(pnlpt->ln_pk_Id2d2);
        free(pnlpt->ln_pk_Id2d2_2);
        free(pnlpt->ln_pk_Id2d2_4);
        
        free(pnlpt->ln_pk_Id2);
        free(pnlpt->ln_pk_IG2);
        free(pnlpt->ln_pk_Id2G2);
        free(pnlpt->ln_pk_IG2G2);
        
        free(pnlpt->ln_pk_Id2G2_2);
        free(pnlpt->ln_pk_IG2G2_2);
        
        free(pnlpt->ln_pk_Id2G2_4);
        free(pnlpt->ln_pk_IG2G2_4);
        
        free(pnlpt->ln_pk_IFG2);
        free(pnlpt->ln_pk_IFG2_0);
        free(pnlpt->ln_pk_IFG2_0b1);
        free(pnlpt->ln_pk_IFG2_2);
        free(pnlpt->ln_pk_CTR);
        free(pnlpt->ln_pk_CTR_0);
        free(pnlpt->ln_pk_CTR_2);
        free(pnlpt->ln_pk_CTR_4);
        free(pnlpt->ln_pk_Tree);
        free(pnlpt->ln_pk_Tree_0_vv);
        free(pnlpt->ln_pk_Tree_0_vd);
        free(pnlpt->ln_pk_Tree_0_dd);
        free(pnlpt->ln_pk_Tree_2_vv);
        free(pnlpt->ln_pk_Tree_2_vd);
        free(pnlpt->ln_pk_Tree_4_vv);
        
        free(pnlpt->ln_pk_0_vv);
        free(pnlpt->ln_pk_0_vd);
        free(pnlpt->ln_pk_0_dd);
        
        free(pnlpt->ln_pk_2_vv);
        free(pnlpt->ln_pk_2_vd);
        free(pnlpt->ln_pk_2_dd);
        
        free(pnlpt->ln_pk_4_vv);
        free(pnlpt->ln_pk_4_vd);
        free(pnlpt->ln_pk_4_dd);
        
        free(pnlpt->ln_pk_0_b1b2);
        free(pnlpt->ln_pk_0_b1bG2);
        free(pnlpt->ln_pk_0_b2);
        free(pnlpt->ln_pk_0_bG2);
        
        free(pnlpt->ln_pk_2_b1b2);
        free(pnlpt->ln_pk_2_b1bG2);
        free(pnlpt->ln_pk_2_b2);
        free(pnlpt->ln_pk_2_bG2);
        
        free(pnlpt->ln_pk_4_b2);
        free(pnlpt->ln_pk_4_bG2);
        free(pnlpt->ln_pk_4_b1b2);
        free(pnlpt->ln_pk_4_b1bG2);
        
        free(pnlpt->gauss_x);
        free(pnlpt->gauss);
        free(pnlpt->gauss_w);
        
 //       printf(" 'nonlinear_pt_free' module executed successfully\n");

    }
  }
  return _SUCCESS_;

}

// This function gets the linear power spectrum
int nonlinear_pt_pk_l(
                   struct background *pba,
                   struct perturbs *ppt,
                   struct primordial *ppm,
                   struct nonlinear_pt *pnlpt,
                   int index_tau,
                   double *pk_l,
                   double *lnk,
                   double *lnpk,
                   double *ddlnpk) {

  int index_md;
  int index_ic=0;
  int index_type;
  int index_k;
  int index_ic1,index_ic2,index_ic1_ic2;
  double * primordial_pk;
  double source_ic1,source_ic2;

  index_md = ppt->index_md_scalars;

  class_alloc(primordial_pk,ppm->ic_ic_size[index_md]*sizeof(double),pnlpt->error_message);

  for (index_k=0; index_k<pnlpt->k_size; index_k++) {

    class_call(primordial_spectrum_at_k(ppm,
                                        index_md,
                                        linear,
                                        pnlpt->k[index_k],
                                        primordial_pk),
               ppm->error_message,
               pnlpt->error_message);

    pk_l[index_k] = 0;

    /* part diagonal in initial conditions */
    for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
      index_ic1_ic2 = index_symmetric_matrix(index_ic1,index_ic1,ppm->ic_size[index_md]);

    if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {
      if (pnlpt->cb == _TRUE_) {
        source_ic1 = (pba->Omega0_cdm*pnlpt->sources_tp_delta_cdm[index_tau*pnlpt->k_size+index_k]+pba->Omega0_b*pnlpt->sources_tp_delta_b[index_tau*pnlpt->k_size+index_k])/(pba->Omega0_cdm+pba->Omega0_b);
      } else {
        source_ic1 = pnlpt->sources_tp_delta_m[index_tau*pnlpt->k_size+index_k];
      }
    } else {
      if (pnlpt->cb == _TRUE_) {
        source_ic1 = (pba->Omega0_cdm*ppt->sources[index_md][index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_cdm][index_tau * ppt->k_size[index_md] + index_k]+pba->Omega0_b*ppt->sources[index_md][index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_b][index_tau * ppt->k_size[index_md] + index_k])/(pba->Omega0_cdm+pba->Omega0_b);
      } else {
        source_ic1 = ppt->sources[index_md][index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m][index_tau * ppt->k_size[index_md] + index_k];
      }
    }
//printf("Omegacdm=%f   Omegab=%f",pba->Omega0_cdm,pba->Omega0_b);
//printf("k=%f source_m=%f   source_cdmb=%f\n",pnlpt->k[index_k],ppt->sources[index_md][index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m][index_tau * ppt->k_size[index_md] + index_k],source_ic1);
       
       // source_ic1 are transfer functions
        
      pk_l[index_k] += 2.*_PI_*_PI_/pow(pnlpt->k[index_k],3)
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
          pk_l[index_k] += 2.*2.*_PI_*_PI_/pow(pnlpt->k[index_k],3)
            *source_ic1*source_ic2
            *primordial_pk[index_ic1_ic2]; // extra 2 factor (to include the symmetric term ic2,ic1)
        }
      }
    }

    lnk[index_k] = log(pnlpt->k[index_k]);
    lnpk[index_k] = log(pk_l[index_k]);
 //   printf("%e  %e\n",pnlpt->k[index_k],pk_l[index_k]);
  }

    class_call(array_spline_table_columns(lnk,
                                          pnlpt->k_size,
                                          lnpk,
                                          1,
                                          ddlnpk,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

  free(primordial_pk);

  return _SUCCESS_;

} // end of non-linear_Pk_l. This function executes the linear power spectrum


/* beginning of the main function */
 
 int nonlinear_pt_loop(
                       struct precision *ppr,
                       struct background *pba,
                       struct primordial *ppm,
                       struct thermo * pth,
                       struct nonlinear_pt *pnlpt,
                       double tau,
                       double f,
                       double hratio,
                       double Dratio,
                       double *pk_l,
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
                       double *lnk_l,
                       double *lnpk_l
                       ) {

int index_k = 0;
int index_j = 0;
int index_i = 0;
int index_l = 0;
     
double sigmav=0.;
double * pvecback;

int last_index = 0;
int last_index2 = 0;

     
class_alloc(pvecback,pba->bg_size*sizeof(double),pnlpt->error_message);

// class_call only calls a function!

class_call(background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback),
                pba->error_message,
                pnlpt->error_message);

// printf("pnlpt->cb=%i\n",pnlpt->cb);

if (pnlpt->nonlinear_pt_verbose > 0)
printf("Computing one-loop power spectra at z=%e\n",pba->a_today/pvecback[pba->index_bg_a]-1.);
//printf("Alcock-Paczynski effect!\n");

free(pvecback);
     
int Nmax = ppr->nmax_nlpt;
int Nmaxf = ppr->nmax_nlpt+1;
double Nmaxd = Nmax * 1.;
double kmin = 0.00005 * pba->h;
double kmax = 100. * pba->h;
     
     /* If you generate new PT matrices, don't forget to choose kmin and kmax appropriately ! */
     
//double kmin = pnlpt->k[0];
//double kmax = pnlpt->k[pnlpt->k_size-1];
     
double *js;
class_alloc(js,(Nmax+1)*sizeof(double),pnlpt->error_message);
double *kdisc;
class_alloc(kdisc,Nmax * sizeof(double),pnlpt->error_message);
double *Pbin;
class_alloc(Pbin,Nmax * sizeof(double),pnlpt->error_message);
double *Pdisc;
class_alloc(Pdisc,Nmax * sizeof(double),pnlpt->error_message);
double *Ptree;
class_alloc(Ptree,Nmax * sizeof(double),pnlpt->error_message);

     /*
double *pk_10;
class_alloc(pk_10,pnlpt->k_size * sizeof(double),pnlpt->error_message);
double *pk_10b1;
class_alloc(pk_10b1,pnlpt->k_size * sizeof(double),pnlpt->error_message);
double *pk_12;
class_alloc(pk_12,pnlpt->k_size * sizeof(double),pnlpt->error_message);*/

double *P10;
class_alloc(P10,Nmax*sizeof(double),pnlpt->error_message);
double *P10b1;
class_alloc(P10b1,Nmax*sizeof(double),pnlpt->error_message);
double *P12;
class_alloc(P12,Nmax*sizeof(double),pnlpt->error_message);
     
     /*
     int i_kdisc1;
     for (i_kdisc1=0; i_kdisc1< Nmax; i_kdisc1++){
         P_IFG2_0b1_x[i_kdisc1] = 0.;
         P_IFG2_0[i_kdisc1] = 0.;
         P_IFG2_2[i_kdisc1] = 0.;
     }
     */
     
double *Pnw;
class_alloc(Pnw,Nmax * sizeof(double),pnlpt->error_message);
     
double *Pw;
class_alloc(Pw,Nmax * sizeof(double),pnlpt->error_message);
     
double *dd_Pnw;
class_alloc(dd_Pnw,sizeof(double)*Nmax,pnlpt->error_message);
double *dd_Pw;
class_alloc(dd_Pw,sizeof(double)*Nmax,pnlpt->error_message);
     
double SigmaBAO = 0.;
double deltaSigmaBAO = 0.;
     // double kmin = pnl->k[0];
double Delta = log(kmax / kmin) / (Nmaxd - 1);
     
     // Interpolating the linear power spectrum
     
double lnpk_out = 0.;
double *myddlnpk;
class_alloc(myddlnpk,sizeof(double)*pnlpt->k_size,pnlpt->error_message);
class_call(array_spline_table_columns(lnk_l,
                                      pnlpt->k_size,
                                      lnpk_l,
                                      1,
                                      myddlnpk,
                                      _SPLINE_NATURAL_,
                                      pnlpt->error_message),
                pnlpt->error_message,
                pnlpt->error_message);
	    
// printf("kmin * exp(0 * Delta)=%f\n",kmin);
// printf("kmin * exp((Nmax-1) * Delta)=%f\n",kmin * exp((Nmax-1) * Delta));
// printf("kmax =%f\n",kmax);

last_index=0;
double delta_array = 1.e-9;
int i_kdisc = 0;

for (i_kdisc=0; i_kdisc< Nmax; i_kdisc++){

         kdisc[i_kdisc] = kmin * exp(i_kdisc * Delta);

         if (kdisc[i_kdisc] >= exp(lnk_l[0]) ){

         class_call(array_interpolate_spline(lnk_l,
                                             pnlpt->k_size,
                                             lnpk_l,
                                             myddlnpk,
                                             1,
                                             log(kdisc[i_kdisc]),
                                             &last_index,
                                             &lnpk_out,
                                             1,
                                             pnlpt->error_message),
                    pnlpt->error_message,
                    pnlpt->error_message);
     

         Pdisc[i_kdisc] = exp(lnpk_out);
       } // end of the condition that our k's are within the range computed by class
 //        printf("%le %le\n",kdisc[i],Pdisc[i]);
       else{
        Pdisc[i_kdisc] = exp(lnpk_l[0])*pow(kdisc[i_kdisc]/exp(lnk_l[0]), ppm->n_s);
       }
         
     }

     int Nside; 

     /* 
     To be safe we're shifting the range of kmaxes and kmins used in the analysis in order to be able to compute the AP, 
     this should be increased if needed
     */

     if (pnlpt->AP_effect == AP_effect_yes){
     Nside = 10;
      }
     else{
     Nside = 1;
     }
     double kmaxnew = kdisc[Nmax-1-Nside];
     double kminnew = kdisc[Nside];
     
		/*
		double kdiff;
		double kdiff2;
		kdiff = kmin * exp((Nmax-1) * Delta) - kmax;
		kdiff2 = kmin * exp((Nmax-1) * Delta) - exp(lnk_l[pnlpt->k_size-1]);
		printf("%le\n",kdiff);
		*/

     /*
		lnpk_out = 0;
                kdisc[Nmax - 1] = kmax;
         class_call(array_interpolate_spline(lnk_l,
                                             pnlpt->k_size,
                                             lnpk_l,
                                             myddlnpk,
                                             1,
                                             log(kmax),
                                             &last_index,
                                             &lnpk_out,
                                             1,
                                             pnlpt->error_message),
                    pnlpt->error_message,
                    pnlpt->error_message);


         Pdisc[Nmax - 1] = exp(lnpk_out);
      */

     
     
//     kdisc[0] = kdisc[0] + delta_array;
//     kdisc[Nmax - 1] = kdisc[Nmax - 1] -delta_array;


     
     // Computing sigmav
    sigmav = 0.;
    for (index_k=0; index_k< pnlpt->k_size-1; index_k++){
         sigmav += (lnk_l[index_k+1]-lnk_l[index_k]) *(exp(lnk_l[index_k+1])*exp(lnpk_l[index_k+1]) + exp(lnpk_l[index_k])*exp(lnk_l[index_k]))/2. / (6. * pow(M_PI,2.));
     }
    // printf("%le\n",sigmav);

     int irindex = 0;
     
     if (pnlpt->irres == irres_yes) {
         
     irindex = 1;
     
         
     if (pnlpt->nonlinear_pt_verbose > 0)
     printf("Performing IR resummation...\n");
         
     //   IR-1) Computing the DFST-II of log(kP)
     
     int Nirby4 = 262144;
     int Nir = 65536;
     int Nirover2 = 32768;
     double Nird = Nir * 1.;
     // double kmin2 = 0.0001 * pba->h;
     // double kmax2 = 10.* pba->h;
     double kmin2 = 0.00007;
     double kmax2 = 7;
     double *logkPdiscr;
     class_alloc(logkPdiscr,Nir * sizeof(double),pnlpt->error_message);
     double logPbin2;
     double kbin2;
         
         double *input_realv2;
         class_alloc(input_realv2,Nirby4 * sizeof(double),pnlpt->error_message);
         double *input_imagv2;
         class_alloc(input_imagv2,Nirby4 * sizeof(double),pnlpt->error_message);
         double *output_realv2;
         class_alloc(output_realv2,Nirby4 * sizeof(double),pnlpt->error_message);
         double *output_imagv2;
         class_alloc(output_imagv2,Nirby4 * sizeof(double),pnlpt->error_message);

        
     last_index=0;
     int index_ir = 0;
     for (index_ir=0; index_ir< Nir; index_ir++){

         kbin2 = kmin2 + index_ir * (kmax2 - kmin2)/(Nird-1.);
         
         if (kbin2>=exp(lnk_l[0])){
         class_call(array_interpolate_spline(lnk_l,
                                             pnlpt->k_size,
                                             lnpk_l,
                                             myddlnpk,
                                             1,
                                             log(kbin2),
                                             &last_index,
                                             &logPbin2,
                                             1,
                                             pnlpt->error_message),
                    pnlpt->error_message,
                    pnlpt->error_message);

         logkPdiscr[index_ir] = log(kbin2) + logPbin2;
       }

       else{
        logkPdiscr[index_ir] = log(kbin2) + lnpk_l[0]+ (ppm->n_s)*(log(kbin2)-lnk_l[0]);
       }


         input_realv2[2*index_ir+1] = pow(-1.,1.*index_ir)*logkPdiscr[index_ir];
         input_realv2[Nirby4-2*index_ir-1] = pow(-1.,1.*index_ir)*logkPdiscr[index_ir];

         input_realv2[2*index_ir] = 0.;
         input_realv2[Nirby4-2*index_ir-2] = 0.;
         
         input_imagv2[2*index_ir+1] = 0.;
         input_imagv2[Nirby4-2*index_ir-1] = 0.;
         input_imagv2[2*index_ir] = 0.;
         input_imagv2[Nirby4-2*index_ir-2] = 0.;
         
     }

         int stepsize = 1;
         
         FFT(input_realv2,input_imagv2,output_realv2,output_imagv2,Nirby4,stepsize);
         
 //        FFT(input_real_test,input_imag_test,output_real_test,output_imag_test,Nirby4,stepsize);
         
         double *out_ir;
         class_alloc(out_ir,Nirby4 * sizeof(double),pnlpt->error_message);
         
	int index_ir2 = 0;
         for (index_ir2=0; index_ir2< Nir; index_ir2++){
             out_ir[index_ir2] = output_realv2[Nir - index_ir2 - 1];
         }
         
         free(input_realv2);
         free(input_imagv2);
         free(output_realv2);
         free(output_imagv2);
         
         double *cmodd;
         class_alloc(cmodd,Nirover2 * sizeof(double),pnlpt->error_message);
         double *cmeven;
         class_alloc(cmeven,Nirover2 * sizeof(double),pnlpt->error_message);
         int *ivar;
         class_alloc(ivar,Nirover2 * sizeof(int),pnlpt->error_message);
     
     index_ir = 0;	    
     for (index_ir=0; index_ir<Nirover2; index_ir++){
         ivar[index_ir] = index_ir + 1;
         cmodd[index_ir] = out_ir[2*index_ir];
         cmeven[index_ir] = out_ir[2*index_ir+1];
     }
     
     
     //   IR-2) Erasing the BAO bump from the odd and even Fourier harmonics and smoothy interpolating the remaining Fourier coefficients
         
     int Nleft = 120;
     int Nright = 240;
     // int Nright = 220;
     int Nthrow = Nright - Nleft;
     int Nnew = Nirover2 - Nthrow;
         double *cmoddnw;
         class_alloc(cmoddnw,Nnew * sizeof(double),pnlpt->error_message);
         double *cmevennw;
         class_alloc(cmevennw,Nnew * sizeof(double),pnlpt->error_message);
         double *inew;
         class_alloc(inew,Nnew * sizeof(double),pnlpt->error_message);
     
     index_ir = 0;
     for (index_ir=0; index_ir<Nleft; index_ir++){
         cmoddnw[index_ir] = cmodd[index_ir];
         cmevennw[index_ir] = cmeven[index_ir];
         inew[index_ir] = index_ir+1.;
     }
     
     index_ir = 0;
     for (index_ir = Nleft; index_ir<Nnew; index_ir++){
         cmoddnw[index_ir] = cmodd[index_ir + Nthrow];
         cmevennw[index_ir] = cmeven[index_ir + Nthrow];
         inew[index_ir] = index_ir+1. + Nthrow;
     }
         
         
         last_index=0;
         double *dd_cmoddnw;
         double cmodd_newval;
         class_alloc(dd_cmoddnw,sizeof(double)*Nnew,pnlpt->error_message);
         class_call(array_spline_table_columns(inew,
                                               Nnew,
                                               cmoddnw,
                                               1,
                                               dd_cmoddnw,
                                               _SPLINE_NATURAL_,
                                               pnlpt->error_message),
                    pnlpt->error_message,
                    pnlpt->error_message);
         double *dd_cmevennw;
         double cmeven_newval;
         class_alloc(dd_cmevennw,sizeof(double)*Nnew,pnlpt->error_message);
         class_call(array_spline_table_columns(inew,
                                               Nnew,
                                               cmevennw,
                                               1,
                                               dd_cmevennw,
                                               _SPLINE_NATURAL_,
                                               pnlpt->error_message),
                    pnlpt->error_message,
                    pnlpt->error_message);
     
         double *cmnew;
         class_alloc(cmnew,Nir * sizeof(double),pnlpt->error_message);
         
     last_index=0;
     last_index2=0;
     index_ir = 0;
     for (index_ir=0; index_ir<Nirover2; index_ir++ ){
         
         class_call(array_interpolate_spline(inew,
                                             Nnew,
                                             cmoddnw,
                                             dd_cmoddnw,
                                             1,
                                             1.*index_ir + 1.,
                                             &last_index,
                                             &cmodd_newval,
                                             1,
                                             pnlpt->error_message),
                    pnlpt->error_message,
                    pnlpt->error_message);
         
         class_call(array_interpolate_spline(inew,
                                             Nnew,
                                             cmevennw,
                                             dd_cmevennw,
                                             1,
                                             1.*index_ir + 1.,
                                             &last_index2,
                                             &cmeven_newval,
                                             1,
                                             pnlpt->error_message),
                    pnlpt->error_message,
                    pnlpt->error_message);
         
         cmnew[index_ir*2] = cmodd_newval;
         cmnew[index_ir*2+1] = cmeven_newval;
     }
     
     //   IR-3) Inverse DST-II (= DST-III/(2*N)) and interpolating P_nw
         
         double *out_2;
         class_alloc(out_2,Nir * sizeof(double),pnlpt->error_message);
         double *input_realv3;
         class_alloc(input_realv3,Nirby4 * sizeof(double),pnlpt->error_message);
         double *input_imagv3;
         class_alloc(input_imagv3,Nirby4 * sizeof(double),pnlpt->error_message);
         double *output_realv3;
         class_alloc(output_realv3,Nirby4 * sizeof(double),pnlpt->error_message);
         double *output_imagv3;
         class_alloc(output_imagv3,Nirby4 * sizeof(double),pnlpt->error_message);
         
         input_realv3[0] = cmnew[Nir-1]*0.5;
         input_realv3[Nir] = 0.;
         input_realv3[2*Nir] = -1.*cmnew[Nir-1]*0.5;
         input_realv3[3*Nir] = 0.;
         
         input_imagv3[0] = 0.;
         input_imagv3[Nir] = 0.;
         input_imagv3[2*Nir] = 0.;
         input_imagv3[3*Nir] = 0.;
         
         index_ir = 0;
         for (index_ir=1; index_ir<Nir; index_ir++){
             input_realv3[index_ir] = 0.5*cmnew[Nir-1-index_ir];
             input_realv3[4*Nir - index_ir] =0.5*cmnew[Nir-1-index_ir];
             input_realv3[2*Nir - index_ir] = -0.5*cmnew[Nir-1-index_ir];
             input_realv3[2*Nir + index_ir] = -0.5*cmnew[Nir-1-index_ir];
             
             input_imagv3[index_ir] = 0.;
             input_imagv3[4*Nir - index_ir] = 0.;
             input_imagv3[2*Nir - index_ir] = 0.;
             input_imagv3[2*Nir + index_ir] = 0.;
         }
         
         FFT(input_realv3,input_imagv3,output_realv3,output_imagv3,Nirby4,stepsize);
     
       //  double out_3[Nir];
         
         index_ir = 0;
         for (index_ir=0; index_ir< Nir; index_ir++){
             out_2[index_ir] = pow(-1.,index_ir)*output_realv3[2*index_ir+1];
         //    out_3[i] = output_realv3[2*i+1]/logkPdiscr[i];
        //     printf("%e\n",out_3[i]);
             
         }
         
         free(input_realv3);
         free(input_imagv3);
         free(output_realv3);
         free(output_imagv3);
         
         
         double *Pnw_ir;
         class_alloc(Pnw_ir,Nir * sizeof(double),pnlpt->error_message);
         double *knw_ir;
         class_alloc(knw_ir,Nir * sizeof(double),pnlpt->error_message);
         
         index_ir = 0;
     for (index_ir=0; index_ir<Nir; index_ir++){
         knw_ir[index_ir] = kmin2 + index_ir * (kmax2 - kmin2)/(Nird-1.);
         Pnw_ir[index_ir] = exp(out_2[index_ir]/(2*Nird))/(knw_ir[index_ir]);
         //       printf("%le %le\n",knw[i],Pnw[i]);
     }
         
     free(dd_cmoddnw);
     free(dd_cmevennw);
         
         double *ddPnw;
         double Pnwval;
         class_alloc(ddPnw,sizeof(double)*Nir,pnlpt->error_message);
         class_call(array_spline_table_columns(knw_ir,
                                               Nir,
                                               Pnw_ir,
                                               1,
                                               ddPnw,
                                               _SPLINE_NATURAL_,
                                               pnlpt->error_message),
                    pnlpt->error_message,
                    pnlpt->error_message);
         
     //   IR-4) Computing the BAO damping factor
     
     double rbao  = pth->rs_d;  
  // double rbao = 110./pba->h; 

     // rbao = pth->rs_d;
    // printf("pth->rs_d=%lf\n",pth->rs_d);
    // printf("rbao=%lf\n",rbao);


     int Nint2 = 500;
         
         double *qint2;
         class_alloc(qint2,(Nint2+1) * sizeof(double),pnlpt->error_message);
         
         double *IntegrandBAO;
         class_alloc(IntegrandBAO,(Nint2+1) * sizeof(double),pnlpt->error_message);
         
         double *IntegrandBAO2;
         class_alloc(IntegrandBAO2,(Nint2+1) * sizeof(double),pnlpt->error_message);
         
         last_index=0;
         
         double ks = 0.2 * pba->h;
         //double ks = 0.5 * pba->h;
         index_ir = 0;
         for (index_ir=0; index_ir<Nint2+1; index_ir++){
             
             qint2[index_ir] = kmin2 * exp(index_ir * log(ks/kmin2) / (Nint2));
             
             
             class_call(array_interpolate_spline(knw_ir,
                                                 Nir,
                                                 Pnw_ir,
                                                 ddPnw,
                                                 1,
                                                 qint2[index_ir],
                                                 &last_index,
                                                 &Pnwval,
                                                 1,
                                                 pnlpt->error_message),
                        pnlpt->error_message,
                        pnlpt->error_message);
             
             
             IntegrandBAO[index_ir] = Pnwval * (1. - 3.*sin(qint2[index_ir] * rbao)/(qint2[index_ir] * rbao) + 6.*(sin(qint2[index_ir] * rbao)/pow((qint2[index_ir] * rbao),3.) - cos(qint2[index_ir] * rbao)/pow((qint2[index_ir] * rbao),2.)));
             
             IntegrandBAO2[index_ir] = -1.*Pnwval * (3.*cos(qint2[index_ir] * rbao)*rbao*qint2[index_ir]+(-3.+pow((qint2[index_ir] * rbao),2.))*sin(qint2[index_ir] * rbao))/pow((qint2[index_ir] * rbao),3.);
         }
     
     SigmaBAO = 0.;
     deltaSigmaBAO = 0.;
         for (index_ir=0; index_ir<Nint2; index_ir++){
             SigmaBAO += (log( qint2[index_ir+1])-log( qint2[index_ir])) *(qint2[index_ir+1]*IntegrandBAO[index_ir+1] + qint2[index_ir]*IntegrandBAO[index_ir])/2. / (6. * pow(M_PI,2.));
             deltaSigmaBAO += (log( qint2[index_ir+1])-log( qint2[index_ir])) *(qint2[index_ir+1]*IntegrandBAO2[index_ir+1] + qint2[index_ir]*IntegrandBAO2[index_ir])/2. / (2. * pow(M_PI,2.));
         }
//     double SigmaBAOh = SigmaBAO * pow(pba->h,2.);
 //    printf("Sigma_BAO(ks=0.2 h/Mpc)=%lf (Mpc/h)^2\n",SigmaBAO * pow(pba->h,2.));
   //  printf("deltaSigma_BAO(ks=0.2 h/Mpc)=%lf (Mpc/h)^2\n",deltaSigmaBAO * pow(pba->h,2.));

         // done for excersise purposes. REMOVE! 
       // SigmaBAO = SigmaBAO/4.;
       // deltaSigmaBAO = deltaSigmaBAO/4.;




     //   IR-5) Computing the LO IR resummed power spectrum
         
         double Pnwval2;
         last_index=0;
         int index_kdisc = 0;
         for (index_kdisc=0; index_kdisc<Nmax; index_kdisc++){
         
             if (kdisc[index_kdisc]<= kmax2 && kdisc[index_kdisc] >= kmin2) {
                 
                 class_call(array_interpolate_spline(knw_ir,
                                                     Nir,
                                                     Pnw_ir,
                                                     ddPnw,
                                                     1,
                                                     kdisc[index_kdisc],
                                                     &last_index,
                                                     &Pnwval2,
                                                     1,
                                                     pnlpt->error_message),
                            pnlpt->error_message,
                            pnlpt->error_message);
                 
                 Pnw[index_kdisc] = Pnwval2;
                 Pw[index_kdisc] = Pdisc[index_kdisc] - Pnw[index_kdisc];
                 // (uncomment here if you want to test the fake NW part only)
            //     Pw[index_kdisc] = 0.;
                 Pbin[index_kdisc] = Pnw[index_kdisc] + Pw[index_kdisc] * exp(-SigmaBAO * pow(kdisc[index_kdisc],2.));
             }
         
         else {
                Pnw[index_kdisc] = Pdisc[index_kdisc];
                Pw[index_kdisc] = 0.;
                Pbin[index_kdisc] = Pdisc[index_kdisc];
         }
         Ptree[index_kdisc] = Pnw[index_kdisc] + Pw[index_kdisc] * exp(-SigmaBAO * pow(kdisc[index_kdisc],2.))*(1. + SigmaBAO * pow(kdisc[index_kdisc],2.));
     }
     
         free(ddPnw);
         free(out_ir);
         free(out_2);
         free(cmnew);
         free(logkPdiscr);
         free(cmeven);
         free(cmodd);
         free(ivar);
         free(inew);
         free(cmoddnw);
         free(cmevennw);
         free(knw_ir);
         free(Pnw_ir);
         free(qint2);
         free(IntegrandBAO);
         free(IntegrandBAO2);
       //  free(Pw);
       //  free(Pnw);
} /* End of IR resummation conditional expression */
     
     else{
         
        if (pnlpt->nonlinear_pt_verbose > 0)
         printf("IR resummation skipped.\n");
         int index_kd = 0;
         for (index_kd=0; index_kd<Nmax; index_kd++){
          Pbin[index_kd] = Pdisc[index_kd];
          Ptree[index_kd] = Pbin[index_kd];
             Pnw[index_kd] = Pdisc[index_kd];
             Pw[index_kd] =0.;
        }
       }
     
     
     class_call(array_spline_table_columns(kdisc,Nmax,Pnw, 1, dd_Pnw,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
     last_index=0;
     double Pnw_ap_out=0;
     

     class_call(array_spline_table_columns(kdisc,Nmax,Pw, 1, dd_Pw,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
     double Pw_ap_out=0;
     
     // here we compute the FFT coefficients
     
     double complex *etam;
     class_alloc(etam,(Nmax+1)*sizeof(complex double),pnlpt->error_message);
     
     int index_c = 0;
     double b = -0.3;
     for (index_c=0; index_c< Nmax +1 ; index_c++){
         js[index_c] = index_c - Nmaxd/2;
         etam[index_c] = b + 2. * M_PI * _Complex_I * js[index_c]/Nmaxd / Delta ;
     }
     
     double *input_real;
     class_alloc(input_real,(Nmax)*sizeof(double),pnlpt->error_message);
     double *input_imag;
     class_alloc(input_imag,(Nmax)*sizeof(double),pnlpt->error_message);
     double *output_real;
     class_alloc(output_real,(Nmax)*sizeof(double),pnlpt->error_message);
     double *output_imag;
     class_alloc(output_imag,(Nmax)*sizeof(double),pnlpt->error_message);

     int stepsize = 1;
     
     int index_kd = 0;
     for (index_kd=0; index_kd< Nmax ; index_kd++){
         input_real[index_kd] = Pbin[index_kd]* exp(-1.* index_kd * b* Delta);
         input_imag[index_kd] = 0.;
     }

     
     /*
     for (size_t i=0; i< N; i++){
         input_real_1[i] =in[i][0];
         input_real_2[i] =in[i][1];
     }*/
     
 //    FFT_real(input_real_1,input_real_2,output_real_1,output_imag_1,output_real_2,output_imag_2,N);
 //    FFT(input_real,input_imag,output_real,output_imag,N,stepsize);
     
     
     FFT(input_real,input_imag,output_real,output_imag,Nmax,stepsize);
     
 /*    output_real_1[Nmax] = output_real_1[0]; */
     
/*
     output_real_1[Nmax] = 0.;
     for (size_t i=0; i< Nmax; i++){
         output_real_1[Nmax] += in[i][0];
         //((double)(Nmax))
     }
 */
     
   /* printf("Start of cn's\n");
      for (size_t i=0; i< Nmax+1; i++){
          printf("%e %e\n",output_real[i],output_imag[i]);
     // printf("%f %f\n",out[i][0],out[i][1]);
      //    printf("%e %e\n",output_real_1[i],output_imag_1[i]);
     //    printf("%e %e\n",out[i][0],out[i][1]);
      }
      printf("End of cn's\n");
    */
     
     double complex *cmsym;
     class_alloc(cmsym,(Nmax+1)*sizeof(complex double),pnlpt->error_message);
     
     index_c = 0;
     for (index_c=0; index_c< Nmax+1; index_c++){
         if (index_c < Nmax/2) {
             cmsym[index_c]= cpow(kmin,-etam[index_c]) * (output_real[Nmax/2 - index_c] - _Complex_I * output_imag[Nmax/2 - index_c])/Nmaxd;
         }
         else {
             cmsym[index_c]= cpow(kmin,-etam[index_c]) * (output_real[index_c - Nmax/2] + _Complex_I * output_imag[index_c - Nmax/2])/Nmaxd;
         }
     }
     
     cmsym[0] = cmsym[0]/2.;
     cmsym[Nmax] = cmsym[Nmax]/2.;

     free(input_real);
     free(input_imag);
     free(output_real);
     free(output_imag);
     
     /*
     printf("Start of cnsym's\n");
     
     for (size_t i=0; i< Nmax+1; i++){
         printf("%le  %le\n",creal(cmsym[i]),cimag(cmsym[i]));
     }
     
     printf("End of cnsym's\n");
      */
     
// here we input the precomputed PT matrices
     
     double complex nu1,nu2,nu12;
     
     double cutoff = 3. * pba->h;
     double *P13;
     double *P13UV;
     double *P1loop;
     double complex *f13;
     class_alloc(P13,Nmax*sizeof(double),pnlpt->error_message);
     class_alloc(P13UV,Nmax*sizeof(double),pnlpt->error_message);
     class_alloc(P1loop,Nmax*sizeof(double),pnlpt->error_message);
     class_alloc(f13,Nmax*sizeof(complex double),pnlpt->error_message);
     
     double complex *f22;
     double *P22;
     double *P_CTR;
     class_alloc(f22,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P22,Nmax*sizeof(double),pnlpt->error_message);
     class_alloc(P_CTR,Nmax*sizeof(double),pnlpt->error_message);

     int count=0;
     int start1=clock();
     char uplo = 'L';
     int inc = 1;
     double complex alpha = 1.;
     double complex beta = 0.;

     double complex *x;
     double complex *x_w;
     double complex *y;
     class_alloc(x,(Nmax+1)*sizeof(complex double),pnlpt->error_message);
     class_alloc(x_w,(Nmax+1)*sizeof(complex double),pnlpt->error_message);
     class_alloc(y,(Nmax+1)*sizeof(complex double),pnlpt->error_message);

//      int Ntest = 2;    
//      double complex *test1;
//      double complex *test2;
// class_alloc(test1,Ntest*sizeof(complex double),pnlpt->error_message);
// class_alloc(test2,Ntest*sizeof(complex double),pnlpt->error_message); 


// for (count=0; count<Ntest; count++){
// test2[count] = 2.+ 5.*count;
// test1[count] = 1.+ 5.*count;
// printf("%le %le\n",creal(test1[count]),creal(test2[count]));
// }

// double complex tt;
// tt = zdotu_(&Ntest, test1, &inc, test2, &inc);

// printf("%le\n",tt);

// free(test1);
// free(test2);

   //  if (pnlpt->rsd == rsd_yes && pnlpt->rsd_only == rsd_only_no || pnlpt->rsd == rsd_no) {
     
     for (index_j=0; index_j < Nmax; index_j++){
         f13[index_j]=0.;
         for (count=0; count < Nmax+1; count++){
             x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
     //        printf(" %le \n", x[count]);
         }
        f13[index_j]=zdotu_(&Nmaxf, x, &inc, pnlpt->M13_oneline_complex, &inc);
        P13UV[index_j] = -61. * Pbin[index_j] * pow(kdisc[index_j], 2.) * sigmav / 105.;
        P13[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13[index_j] * Pbin[index_j]) + P13UV[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
         }
    index_j = 0;
    count = 0;
     for (index_j=0; index_j < Nmax; index_j++){
         for (count=0; count < Nmax+1; count++){
             x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
         }
         
         zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_complex, x, &inc, &beta, y, &inc);
         f22[index_j]=zdotu_(&Nmaxf, x, &inc, y, &inc);
         //printf("f22_real=%.18le f22_imag=%.18le\n",creal(f22[index_j]),cimag(f22[index_j]));
         P22[index_j] = creal(cpow(kdisc[index_j], 3.) * f22[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
         P1loop[index_j] = 0.*Ptree[index_j] +(P13[index_j] + P22[index_j]);
 //        P_CTR[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j];
//         printf("%le %le\n",kdisc[j],P22[j]);
     }
     //} //End of RSD only condition
         
     for (index_j=0; index_j < Nmax; index_j++){
         P_CTR[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j];
     }
     
     double *ddpk_nl;
     class_alloc(ddpk_nl,sizeof(double)*Nmax,pnlpt->error_message);
     double *ddpk_CTR;
     class_alloc(ddpk_CTR,sizeof(double)*Nmax,pnlpt->error_message);
     double *ddpk_Tree;
     class_alloc(ddpk_Tree,sizeof(double)*Nmax,pnlpt->error_message);
     
     
     class_call(array_spline_table_columns(kdisc,
                                           Nmax,
                                           P1loop,
                                           1,
                                           ddpk_nl,
                                           _SPLINE_NATURAL_,
                                           pnlpt->error_message),
                pnlpt->error_message,
                pnlpt->error_message);


     class_call(array_spline_table_columns(kdisc,
                                           Nmax,
                                           P_CTR,
                                           1,
                                           ddpk_CTR,
                                           _SPLINE_NATURAL_,
                                           pnlpt->error_message),
                pnlpt->error_message,
                pnlpt->error_message); 
 
     class_call(array_spline_table_columns(kdisc,
                                           Nmax,
                                           Ptree,
                                           1,
                                           ddpk_Tree,
                                           _SPLINE_NATURAL_,
                                           pnlpt->error_message),
                pnlpt->error_message,
                pnlpt->error_message);

     
     double pk_nl_out;
     double pk_CTR_out;
     double pk_Tree_out;
 
     last_index=0;
     last_index2=0;
     int last_index3 = 0;
     
     for (index_k=0; index_k < pnlpt->k_size; index_k++){
         
         if (pnlpt->k[index_k]>=kmin && pnlpt->k[index_k]<=kmax){
         
         class_call(array_interpolate_spline(kdisc,
                                             Nmax,
                                             P1loop,
                                             ddpk_nl,
                                             1,
                                             pnlpt->k[index_k],
                                             &last_index2,
                                             &pk_nl_out,
                                             1,
                                             pnlpt->error_message),
                    pnlpt->error_message,
                    pnlpt->error_message);
         
         class_call(array_interpolate_spline(kdisc,
                                             Nmax,
                                             P_CTR,
                                             ddpk_CTR,
                                             1,
                                             pnlpt->k[index_k],
                                             &last_index,
                                             &pk_CTR_out,
                                             1,
                                             pnlpt->error_message),
                    pnlpt->error_message,
                    pnlpt->error_message);
         

         class_call(array_interpolate_spline(kdisc,
                                             Nmax,
                                             Ptree,
                                             ddpk_Tree,
                                             1,
                                             pnlpt->k[index_k],
                                             &last_index3,
                                             &pk_Tree_out,
                                             1,
                                             pnlpt->error_message),
                    pnlpt->error_message,
                    pnlpt->error_message);

             pk_nl[index_k] = pk_nl_out + 5000.;
             pk_CTR[index_k] = pk_CTR_out;
	         pk_Tree[index_k] = pk_Tree_out;
         }
         
         else {
             //pk_nl[index_k] = exp(lnpk_l[index_k]);
          if (pnlpt->k[index_k]<kmin){
             pk_nl[index_k] = 5000.  - 61. * exp(lnpk_l[index_k]+2.*lnk_l[index_k]) * sigmav/ 105.;
           }
          if (pnlpt->k[index_k]>kmax){
             pk_nl[index_k] = 5000.;
           }
             pk_CTR[index_k] = exp(lnpk_l[index_k]+2.*lnk_l[index_k]) ;
             pk_Tree[index_k] = exp(lnpk_l[index_k]);
         }
   //  printf("%i %f %f \n",index_k, pk_Tree[index_k],pk_nl[index_k]-5000.);
     }
     
if (pnlpt->rsd == rsd_yes){
    
    if (pnlpt->nonlinear_pt_verbose > 0)
    printf("Computing RSD...\n");
    if (pnlpt->nonlinear_pt_verbose > 0)
    printf("Logarithmic growth factor f=%f\n",f);
    
    double *Ptree_0_vv;
    class_alloc(Ptree_0_vv,Nmax * sizeof(double),pnlpt->error_message);
    double *Ptree_0_vd;
    class_alloc(Ptree_0_vd,Nmax * sizeof(double),pnlpt->error_message);
    double *Ptree_0_dd;
    class_alloc(Ptree_0_dd,Nmax * sizeof(double),pnlpt->error_message);
    double *Ptree_2_vv;
    class_alloc(Ptree_2_vv,Nmax * sizeof(double),pnlpt->error_message);
    double *Ptree_2_vd;
    class_alloc(Ptree_2_vd,Nmax * sizeof(double),pnlpt->error_message);
    double *Ptree_4_vv;
    class_alloc(Ptree_4_vv,Nmax * sizeof(double),pnlpt->error_message);
    
    
    double *P13_0_vv;
    double *P13UV_0_vv;
    double *P1loop_0_vv;
    class_alloc(P13_0_vv,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_0_vv,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P1loop_0_vv,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_0_vv;
    class_alloc(P22_0_vv,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P13_0_vd;
    double *P13UV_0_vd;
    double *P1loop_0_vd;
    class_alloc(P13_0_vd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_0_vd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P1loop_0_vd,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_0_vd;
    class_alloc(P22_0_vd,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P13_0_dd;
    double *P13UV_0_dd;
    double *P1loop_0_dd;
    class_alloc(P13_0_dd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_0_dd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P1loop_0_dd,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_0_dd;
    class_alloc(P22_0_dd,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P13_2_vv;
    double *P13UV_2_vv;
    double *P1loop_2_vv;
    class_alloc(P13_2_vv,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_2_vv,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P1loop_2_vv,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_2_vv;
    class_alloc(P22_2_vv,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P13_2_vd;
    double *P13UV_2_vd;
    double *P1loop_2_vd;
    class_alloc(P13_2_vd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_2_vd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P1loop_2_vd,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_2_vd;
    class_alloc(P22_2_vd,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P13_2_dd;
    double *P13UV_2_dd;
    double *P1loop_2_dd;
    class_alloc(P13_2_dd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_2_dd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P1loop_2_dd,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_2_dd;
    class_alloc(P22_2_dd,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P13_4_vv;
    double *P13UV_4_vv;
    double *P1loop_4_vv;
    class_alloc(P13_4_vv,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_4_vv,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P1loop_4_vv,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_4_vv;
    class_alloc(P22_4_vv,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P13_4_vd;
    double *P13UV_4_vd;
    double *P1loop_4_vd;
    class_alloc(P13_4_vd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_4_vd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P1loop_4_vd,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_4_vd;
    class_alloc(P22_4_vd,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P22_4_dd;
    class_alloc(P22_4_dd,Nmax*sizeof(double),pnlpt->error_message);
    double *P1loop_4_dd;
    class_alloc(P1loop_4_dd,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P_CTR_0;
    class_alloc(P_CTR_0,Nmax*sizeof(double),pnlpt->error_message);
    double *P_CTR_2;
    class_alloc(P_CTR_2,Nmax*sizeof(double),pnlpt->error_message);
    double *P_CTR_4;
    class_alloc(P_CTR_4,Nmax*sizeof(double),pnlpt->error_message);
    
    
  //  if (pnlpt->irres == irres_no){
    
  //  if (pnlpt->irres == irres_yes){
    if (irindex == 0){
    
  //  printf("Computing RSD without IR resummation...\n");
        
 // Computing P_{vv} contribution
    
    count = 0;
    index_l = 0;
    index_i = 0;
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu1 = -0.5*etam[index_i];
            nu2 = -0.5*etam[index_l];
            nu12 = nu1+nu2;
            
            pnlpt->M22_oneline_0_vv_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*(f*f*(14.*f*f*(24.-8.*nu2-15.*nu2*nu2+5.*nu2*nu2*nu2+5.*nu1*nu1*nu1*(1.+7.*nu2)+5.*nu1*nu1*(-3.-10.*nu2+14.*nu2*nu2)+nu1*(-8.-24.*nu2-50.*nu2*nu2+35.*nu2*nu2*nu2))+18.*f*(36.-8.*nu2+70.*nu1*nu1*nu1*nu2-23.*nu2*nu2+nu1*nu1*(-23.-94.*nu2+140.*nu2*nu2)+nu1*(-8.-42.*nu2-94.*nu2*nu2+70.*nu2*nu2*nu2))+9.*(50.-9.*nu2+98.*nu1*nu1*nu1*nu2-35.*nu2*nu2+7.*nu1*nu1*(-5.-18.*nu2+28.*nu2*nu2)+nu1*(-9.-66.*nu2-126.*nu2*nu2+98.*nu2*nu2*nu2))))/8820.;
            
            
            count++;
            }
        }
    
    double complex *f13_0_vv;
    class_alloc(f13_0_vv,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_0_vv;
    class_alloc(f22_0_vv,Nmax*sizeof(complex double),pnlpt->error_message);
    index_i = 0;
    for (index_i=0; index_i<Nmax+1; index_i++){
        nu1=-0.5*etam[index_i];
        pnlpt->M13_0_vv_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i]*112./(1.+9.*nu1)*(3.*(f*f)*(7.*(-5. + 3.*nu1) + 6.*f*(-7. + 5.*nu1)))/3920.;
        }
    
    index_j = 0;
    count= 0;
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
            }
        f13_0_vv[index_j]=zdotu_(&Nmaxf, x, &inc, pnlpt->M13_0_vv_oneline_complex, &inc);
        P13UV_0_vv[index_j] = -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*(f*f*(441.+566.*f+175.*f*f)/1225.);
        P13_0_vv[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_0_vv[index_j] * Pbin[index_j]) + P13UV_0_vv[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
        }
        
    index_j = 0;
    count= 0;
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_0_vv_complex, x, &inc, &beta, y, &inc);
        f22_0_vv[index_j]=zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_0_vv[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_0_vv[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        P1loop_0_vv[index_j] = Pbin[index_j]*(f*f/5.)*0.+(P13_0_vv[index_j] + P22_0_vv[index_j]);
        P_CTR_0[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j];
        Ptree_0_vv[index_j] = Pbin[index_j]*(f*f/5.);
    //    printf("%le %le\n",kdisc[j],P22[j]);
        }

    
 // Computing P_{vd} contribution
    
    count = 0;
    index_l= 0;
    index_i= 0;
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu12 =-0.5*etam[index_i]-0.5*etam[index_l];
            nu1=-0.5*etam[index_i];
            nu2=-0.5*etam[index_l];
            pnlpt->M22_oneline_0_vd_complex[count] = pnlpt->M22_oneline_complex[count]*196./((nu1*nu2*(98.*nu12*nu12-14.*nu12+36.)-91.*nu12*nu12+3.*nu12+58.))*(f*(21.*f*f*(6.+3.*nu2-10.*nu2*nu2+2.*nu2*nu2*nu2+2.*nu1*nu1*nu1*(1.+5.*nu2)+2*nu1*nu1*(-5.-2.*nu2+10.*nu2*nu2)+nu1*(3.-24.*nu2-4.*nu2*nu2+10.*nu2*nu2*nu2))+14.*f*(18.+11.*nu2+42.*nu1*nu1*nu1*nu2-31.*nu2*nu2+nu1*nu1*(-31.-22.*nu2+84.*nu2*nu2)+nu1*(11.-74.*nu2-22.*nu2*nu2+42.*nu2*nu2*nu2))+5.*(46.+13.*nu2+98.*nu1*nu1*nu1*nu2-63.*nu2*nu2+7.*nu1*nu1*(-9.-10.*nu2+28.*nu2*nu2)+nu1*(13.-138.*nu2-70.*nu2*nu2+98.*nu2*nu2*nu2))))/1470.;
            count++;
            }
        }

    index_i=0;

    for (index_i=0; index_i<Nmax+1; index_i++){
        nu1=-0.5*etam[index_i];
        pnlpt->M13_0_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i]*112./(1.+9.*nu1)*(f*(-35. - 18.*f + 45.*nu1 + 54.*f*nu1))/840.;
        }
    
    double complex *f13_0_vd;
    class_alloc(f13_0_vd,Nmax*sizeof(complex double),pnlpt->error_message);
    
    double complex *f22_0_vd;
    class_alloc(f22_0_vd,Nmax*sizeof(complex double),pnlpt->error_message);
    index_j=0;
    count = 0;
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        f13_0_vd[index_j]=zdotu_(&Nmaxf, x, &inc, pnlpt->M13_0_vd_oneline_complex, &inc);
        P13UV_0_vd[index_j] = -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*(2.*f*(625. + 558.*f + 315.*f*f)/1575.);
        P13_0_vd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_0_vd[index_j] * Pbin[index_j]) + P13UV_0_vd[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
        }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_0_vd_complex, x, &inc, &beta, y, &inc);
        f22_0_vd[index_j]=zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_0_vd[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_0_vd[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        P1loop_0_vd[index_j] = Pbin[index_j]*(f*2./3.)*0.+(P13_0_vd[index_j] + P22_0_vd[index_j]);
        Ptree_0_vd[index_j] = Pbin[index_j]*(f*2./3.);
        }

     // Computing P_{dd} contribution
    
    count = 0;
    index_l=0;
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu12 =-0.5*etam[index_i]-0.5*etam[index_l];
            nu1=-0.5*etam[index_i];
            nu2=-0.5*etam[index_l];
            pnlpt->M22_oneline_0_dd_complex[count] = pnlpt->M22_oneline_complex[count]*196./((nu1*nu2*(98.*nu12*nu12-14.*nu12+36.)-91.*nu12*nu12+3.*nu12+58.))*(98.*f*f*(4.-2.*nu2-5.*nu2*nu2+nu2*nu2*nu2+nu1*nu1*nu1*(1.+3.*nu2)+nu1*nu1*(-5.+2.*nu2+6.*nu2*nu2)+nu1*(-2.-4.*nu2+2.*nu2*nu2+3.*nu2*nu2*nu2))+70.*f*(10.-nu2+14.*nu1*nu1*nu1*nu2-17.*nu2*nu2+nu1*nu1*(-17.+6.*nu2+28.*nu2*nu2)+nu1*(-1.-22.*nu2+6.*nu2*nu2+14.*nu2*nu2*nu2))+15.*(58.+3.*nu2+98.*nu1*nu1*nu1*nu2-91.*nu2*nu2+7.*nu1*nu1*(-13.-2.*nu2+28.*nu2*nu2)+nu1*(3.-146.*nu2-14.*nu2*nu2+98.*nu2*nu2*nu2)))/2940.;
            count++;
        }
        }
    index_i=0;
    for (index_i=0; index_i<Nmax+1; index_i++){
        nu1=-0.5*etam[index_i];
        pnlpt->M13_0_dd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i]/(1.+9.*nu1)*(1.+9.*nu1+6.*f*(1.+nu1));
        }
    

    double complex *f13_0_dd;
    class_alloc(f13_0_dd,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_0_dd;
    class_alloc(f22_0_dd,Nmax*sizeof(complex double),pnlpt->error_message);
    index_j=0;
    count=0;
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        f13_0_dd[index_j]=zdotu_(&Nmaxf, x, &inc, pnlpt->M13_0_dd_oneline_complex, &inc);
        P13UV_0_dd[index_j] = -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*((61. -2.*f + 35.*f*f)/105.);
        P13_0_dd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_0_dd[index_j] * Pbin[index_j]) + P13UV_0_dd[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
        }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_0_dd_complex, x, &inc, &beta, y, &inc);
        f22_0_dd[index_j]=zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_0_dd[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_0_dd[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        P1loop_0_dd[index_j] = Pbin[index_j]*0.+(P13_0_dd[index_j] + P22_0_dd[index_j]);
        Ptree_0_dd[index_j] = Pbin[index_j];
    }
    
    // Computing P_{vv} contribution - Quadrupole
    
    count = 0;
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            //  nu12 =-0.5*etam[index_i]-0.5*etam[index_l];
            nu1 = -0.5*etam[index_i];
            nu2 = -0.5*etam[index_l];
            
            nu12 = nu1+nu2;
            
            pnlpt->M22_oneline_2_vv_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*(f*f*(396.*(50.-9.*nu2+98.*nu1*nu1*nu1*nu2-35.*nu2*nu2+7.*nu1*nu1*(-5.-18.*nu2+28.*nu2*nu2)+nu1*(-9.-66.*nu2-126.*nu2*nu2+98.*nu2*nu2*nu2))+231.*f*(142.-21.*nu2+280.*nu1*nu1*nu1*nu2-106.*nu2*nu2+2.*nu1*nu1*(-53.-174.*nu2+280.*nu2*nu2)+nu1*(-21.-204.*nu2-348.*nu2*nu2+280.*nu2*nu2*nu2))+49.*f*f*(336.-62.*nu2-255.*nu2*nu2+50.*nu2*nu2*nu2+10.*nu1*nu1*nu1*(5.+56.*nu2)+5.*nu1*nu1*(-51.-142.*nu2+224.*nu2*nu2)+nu1*(-62.-486.*nu2-710.*nu2*nu2+560.*nu2*nu2*nu2))))/135828.;
            
            count++;
        }
        }
    
    double complex *f13_2_vv;
    class_alloc(f13_2_vv,Nmax*sizeof(complex double),pnlpt->error_message);
    //    double *P_CTR0;
    //    class_alloc(P_CTR0,Nmax*sizeof(double),pnlpt->error_message);
    
    double complex *f22_2_vv;
    class_alloc(f22_2_vv,Nmax*sizeof(complex double),pnlpt->error_message);
    
    for (index_i=0; index_i<Nmax+1; index_i++){
        nu1=-0.5*etam[index_i];
        pnlpt->M13_2_vv_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i]*112./(1.+9.*nu1)*(3.*f*f*(-5.+3.*nu1+f*(-6.+5.*nu1)))/196.;
    }
    
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        f13_2_vv[index_j]=zdotu_(&Nmaxf, x, &inc, pnlpt->M13_2_vv_oneline_complex, &inc);
        P13UV_2_vv[index_j] = -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*(2.*f*f*(54.+74.*f+25.*f*f)/105.);
        P13_2_vv[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_2_vv[index_j] * Pbin[index_j]) + P13UV_2_vv[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
    }
    
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_2_vv_complex, x, &inc, &beta, y, &inc);
        f22_2_vv[index_j]=zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_2_vv[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_2_vv[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        P1loop_2_vv[index_j] = Pbin[index_j]*(f*f*4./7.)*0.+(P13_2_vv[index_j] + P22_2_vv[index_j]);
        P_CTR_2[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j]*f*2./3.;
        Ptree_2_vv[index_j] = Pbin[index_j]*(f*f*4./7.);
        //    printf("%le %le\n",kdisc[j],P22[j]);
    }

    // Computing P_{vd} contribution - Quadrupole
    
    count = 0;
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            //  nu12 =-0.5*etam[index_i]-0.5*etam[index_l];
            nu1 = -0.5*etam[index_i];
            nu2 = -0.5*etam[index_l];
            
            nu12 = nu1+nu2;
            
            pnlpt->M22_oneline_2_vd_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*(f*(7.*f*f*(22.+11.*nu2-40.*nu2*nu2+4.*nu2*nu2*nu2+nu1*nu1*nu1*(4.+40.*nu2)+8.*nu1*nu1*(-5.-nu2+10.*nu2*nu2)+nu1*(11.-88.*nu2-8.*nu2*nu2+40.*nu2*nu2*nu2))+4.*(46.+13.*nu2+98.*nu1*nu1*nu1*nu2-63.*nu2*nu2+7.*nu1*nu1*(-9.-10.*nu2+28.*nu2*nu2)+nu1*(13.-138.*nu2-70.*nu2*nu2 + 98.*nu2*nu2*nu2))+f*(306.+161.*nu2+672.*nu1*nu1*nu1*nu2-538.*nu2*nu2+2.*nu1*nu1*(-269.-134.*nu2+672.*nu2*nu2)+nu1*(161.-1196.*nu2-268.*nu2*nu2+672.*nu2*nu2*nu2))))/588.;
            
            count++;
        }
    }
    
    double complex *f13_2_vd;
    class_alloc(f13_2_vd,Nmax*sizeof(complex double),pnlpt->error_message);
    //    double *P_CTR0;
    //    class_alloc(P_CTR0,Nmax*sizeof(double),pnlpt->error_message);
    
    double complex *f22_2_vd;
    class_alloc(f22_2_vd,Nmax*sizeof(complex double),pnlpt->error_message);
    
    for (index_i=0; index_i<Nmax+1; index_i++){
        nu1=-0.5*etam[index_i];
        pnlpt->M13_2_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i]*112./(1.+9.*nu1)*(f*(-49.-9.*f+63.*nu1+108.*f*nu1))/588.;
    }
    
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        f13_2_vd[index_j]= zdotu_(&Nmaxf, x, &inc, pnlpt->M13_2_vd_oneline_complex, &inc);
        P13UV_2_vd[index_j] = -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*4.*f*(175.+180.*f+126.*f*f)/441.;
        P13_2_vd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_2_vd[index_j] * Pbin[index_j]) + P13UV_2_vd[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
    }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_2_vd_complex, x, &inc, &beta, y, &inc);
        f22_2_vd[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_2_vd[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_2_vd[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        P1loop_2_vd[index_j] = Pbin[index_j]*(f*4./3.)*0.+(P13_2_vd[index_j] + P22_2_vd[index_j]);
        Ptree_2_vd[index_j] = Pbin[index_j]*(f*4./3.);
        //    P_CTR[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j];
        //printf("%le %le\n",kdisc[j],P22[j]);
    }
    
    // Computing P_{dd} contribution - Quadrupole
    
    count = 0;
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            //  nu12 =-0.5*etam[index_i]-0.5*etam[index_l];
            nu1 = -0.5*etam[index_i];
            nu2 = -0.5*etam[index_l];
            
            nu12 = nu1+nu2;
            
            pnlpt->M22_oneline_2_dd_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*f*(4.*(10.-nu2+14.*nu1*nu1*nu1*nu2-17.*nu2*nu2+nu1*nu1*(-17.+6.*nu2+28.*nu2*nu2)+nu1*(-1.-22.*nu2+6.*nu2*nu2+14.*nu2*nu2*nu2))+f*(26.-13.*nu2-37.*nu2*nu2+2.*nu2*nu2*nu2+nu1*nu1*nu1*(2.+24.*nu2)+nu1*nu1*(-37.+22.*nu2+48.*nu2*nu2)+nu1*(-13.-26.*nu2+22.*nu2*nu2+24.*nu2*nu2*nu2)))/84.;
            
            count++;
        }
    }
        
    double complex *f13_2_dd;
    class_alloc(f13_2_dd,Nmax*sizeof(complex double),pnlpt->error_message);
    //    double *P_CTR0;
    //    class_alloc(P_CTR0,Nmax*sizeof(double),pnlpt->error_message);
    
    double complex *f22_2_dd;
    class_alloc(f22_2_dd,Nmax*sizeof(complex double),pnlpt->error_message);
    
    for (index_i=0; index_i<Nmax+1; index_i++){
        nu1=-0.5*etam[index_i];
        pnlpt->M13_2_dd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i]*112./(1.+9.*nu1)*(3.*f*(1.+nu1))/28.;
    }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        f13_2_dd[index_j]=zdotu_(&Nmaxf, x, &inc, pnlpt->M13_2_dd_oneline_complex, &inc);
        P13UV_2_dd[index_j] = -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*(2.*f*(35.*f-2.)/105.);
        P13_2_dd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_2_dd[index_j] * Pbin[index_j]) + P13UV_2_dd[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
    }
    
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_2_dd_complex, x, &inc, &beta, y, &inc);
        f22_2_dd[index_j]=zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_2_dd[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_2_dd[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        P1loop_2_dd[index_j] = P13_2_dd[index_j] + P22_2_dd[index_j];
        //    P_CTR[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j];
        //    printf("%le %le\n",kdisc[j],P22[j]);
    }

    // Computing P_{vv} contribution - Hexadecapole
    
    count = 0;
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu1 = -0.5*etam[index_i];
            nu2 = -0.5*etam[index_l];
            nu12 = nu1+nu2;
            
            pnlpt->M22_oneline_4_vv_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*(f*f*(1144.*(50.+98.*nu1*nu1*nu1*nu2-nu2*(9.+35.*nu2)+7.*nu1*nu1*(-5.+2.*nu2*(-9.+14.*nu2))+nu1*(-9.+2.*nu2*(-33.+7.*nu2*(-9.+7.*nu2))))+147.*f*f*(483.+40.*nu1*nu1*nu1*(-1.+28.*nu2)-2.*nu2*(-57.+10.*nu2*(29.+2.*nu2))+20.*nu1*nu1*(-29.+2.*nu2*(-25.+56.*nu2))+2.*nu1*(57.+2.*nu2*(-327.+10.*nu2*(-25.+28.*nu2))))+728.*f*(206.+420.*nu1*nu1*nu1*nu2+(7.-208.*nu2)*nu2+8.*nu1*nu1*(-26.+nu2*(-53.+105.*nu2))+nu1*(7.+4.*nu2*(-108.+nu2*(-106.+105.*nu2))))))/980980.;
            
            count++;
        }
    }
    
    double complex *f13_4_vv;
    class_alloc(f13_4_vv,Nmax*sizeof(complex double),pnlpt->error_message);
    //    double *P_CTR0;
    //    class_alloc(P_CTR0,Nmax*sizeof(double),pnlpt->error_message);
    
    double complex *f22_4_vv;
    class_alloc(f22_4_vv,Nmax*sizeof(complex double),pnlpt->error_message);
    
    for (index_i=0; index_i<Nmax+1; index_i++){
        nu1=-0.5*etam[index_i];
        pnlpt->M13_4_vv_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i]*112./(1.+9.*nu1)*(3.*f*f*(-55.+33.*nu1+f*(-66.+90.*nu1)))/5390.;
    }
    
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        f13_4_vv[index_j]=zdotu_(&Nmaxf, x, &inc, pnlpt->M13_4_vv_oneline_complex, &inc);
        P13UV_4_vv[index_j] = -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*(24.*f*f*(33.+58.*f+25.*f*f)/1925.);
        P13_4_vv[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_4_vv[index_j] * Pbin[index_j]) + P13UV_4_vv[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
    }
    
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_4_vv_complex, x, &inc, &beta, y, &inc);
        f22_4_vv[index_j]=zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_4_vv[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_4_vv[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        P1loop_4_vv[index_j] = Pbin[index_j]*(f*f*8./35.)*0.+(P13_4_vv[index_j] + P22_4_vv[index_j]);
        P_CTR_4[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j]*(f*f*8./35.);
        Ptree_4_vv[index_j] = Pbin[index_j]*(f*f*8./35.);
        //    printf("%le %le\n",kdisc[j],P22[j]);
    }
    
    // Computing P_{vd} contribution - Hexadecapole
    
    count = 0;
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            //  nu12 =-0.5*etam[index_i]-0.5*etam[index_l];
            nu1 = -0.5*etam[index_i];
            nu2 = -0.5*etam[index_l];
            
            nu12 = nu1+nu2;
            
            pnlpt->M22_oneline_4_vd_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*f*f*(14.*f*(26.+13.*nu2-60.*nu2*nu2-8.*nu2*nu2*nu2+nu1*nu1*nu1*(-8.+60.*nu2)+4.*nu1*nu1*(-15.+4.*nu2+30.*nu2*nu2)+nu1*(13.-104.*nu2+16.*nu2*nu2+60.*nu2*nu2*nu2))+11.*(58.+21.*nu2+112.*nu1*nu1*nu1*nu2-106.*nu2*nu2+2.*nu1*nu1*(-53.-6.*nu2+112.*nu2*nu2)+nu1*(21.-204.*nu2-12.*nu2*nu2+112.*nu2*nu2*nu2)))/2695.;
            
            count++;
        }
    }
    
    double complex *f13_4_vd;
    class_alloc(f13_4_vd,Nmax*sizeof(complex double),pnlpt->error_message);
    //    double *P_CTR0;
    //    class_alloc(P_CTR0,Nmax*sizeof(double),pnlpt->error_message);
    
    double complex *f22_4_vd;
    class_alloc(f22_4_vd,Nmax*sizeof(complex double),pnlpt->error_message);
        
    for (index_i=0; index_i<Nmax+1; index_i++){
        nu1=-0.5*etam[index_i];
        pnlpt->M13_4_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i]*112./(1.+9.*nu1)*9.*(f*f*(1.+2.*nu1))/245.;
    }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        f13_4_vd[index_j]= zdotu_(&Nmaxf, x, &inc, pnlpt->M13_4_vd_oneline_complex, &inc);
        P13UV_4_vd[index_j] = -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*16.*f*f*(22.+35.*f)/1225.;
        P13_4_vd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_4_vd[index_j] * Pbin[index_j]) + P13UV_4_vd[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
    }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_4_vd_complex, x, &inc, &beta, y, &inc);
        f22_4_vd[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_4_vd[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_4_vd[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        P1loop_4_vd[index_j] = (P13_4_vd[index_j] + P22_4_vd[index_j]);
        //    P_CTR[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j];
        //printf("%le %le\n",kdisc[j],P22[j]);
    }
        
        
        // Computing P_{dd} contribution - Hexadecapole
        
        count = 0;
        for (index_l=0; index_l < Nmax+1; index_l++){
            for (index_i=index_l; index_i < Nmax+1; index_i++){
                nu1 = -0.5*etam[index_i];
                nu2 = -0.5*etam[index_l];
                
                nu12 = nu1+nu2;
                
                pnlpt->M22_oneline_4_dd_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*f*f*(2.*nu1-1.)*(2.*nu2-1.)*(1.+nu12)*(2.+nu12)/35.;
                
                count++;
            }
        }
        
        double complex *f22_4_dd;
        class_alloc(f22_4_dd,Nmax*sizeof(complex double),pnlpt->error_message);
        
        for (index_j=0; index_j < Nmax; index_j++){
            for (count=0; count < Nmax+1; count++){
                x[count]= cmsym[count]* cpow(kdisc[index_j], etam[count]);
            }
            zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_4_dd_complex, x, &inc, &beta, y, &inc);
            f22_4_dd[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
            P22_4_dd[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_4_dd[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
            P1loop_4_dd[index_j] = P22_4_dd[index_j];
        }
        
        free(f22_4_dd);
        free(f13_0_vd);
        free(f13_0_vv);
        free(f22_0_vv);
        free(f22_0_vd);
        free(f13_0_dd);
        free(f22_0_dd);
        free(f22_2_vv);
        free(f22_2_vd);
        free(f22_2_dd);
        free(f13_2_vd);
        free(f13_2_dd);
        free(f13_2_vv);
        free(f13_4_vv);
        free(f13_4_vd);
        free(f22_4_vd);
        free(f22_4_vv);
        
    }// end of IR resummation expression

     
    
if (irindex == 1){
   // printf("Computing IR resummed spectra in redshift space...\n");
    
    //Computing FFT for wiggly and non-wiggly parts
    
    double *input_real_nw;
    class_alloc(input_real_nw,(Nmax)*sizeof(double),pnlpt->error_message);
    double *input_imag_nw;
    class_alloc(input_imag_nw,(Nmax)*sizeof(double),pnlpt->error_message);
    double *output_real_nw;
    class_alloc(output_real_nw,(Nmax)*sizeof(double),pnlpt->error_message);
    double *output_imag_nw;
    class_alloc(output_imag_nw,(Nmax)*sizeof(double),pnlpt->error_message);
    
    double *input_real_w;
    class_alloc(input_real_w,(Nmax)*sizeof(double),pnlpt->error_message);
    double *input_imag_w;
    class_alloc(input_imag_w,(Nmax)*sizeof(double),pnlpt->error_message);
    double *output_real_w;
    class_alloc(output_real_w,(Nmax)*sizeof(double),pnlpt->error_message);
    double *output_imag_w;
    class_alloc(output_imag_w,(Nmax)*sizeof(double),pnlpt->error_message);
    
    index_kd = 0;
    for (index_kd=0; index_kd< Nmax ; index_kd++){
        input_real_nw[index_kd] = Pnw[index_kd]* exp(-1.* index_kd * b* Delta);
        input_imag_nw[index_kd] = 0.;
        input_real_w[index_kd] = Pw[index_kd]* exp(-1.* index_kd * b* Delta);
        input_imag_w[index_kd] = 0.;
    }
    
    
    FFT(input_real_nw,input_imag_nw,output_real_nw,output_imag_nw,Nmax,stepsize);
    
    double complex *cmsym_nw;
    class_alloc(cmsym_nw,(Nmax+1)*sizeof(complex double),pnlpt->error_message);
    
    index_c = 0;
    for (index_c=0; index_c< Nmax+1; index_c++){
        if (index_c < Nmax/2) {
            cmsym_nw[index_c]= cpow(kmin,-etam[index_c]) * (output_real_nw[Nmax/2 - index_c] - _Complex_I * output_imag_nw[Nmax/2 - index_c])/Nmaxd;
        }
        else {
            cmsym_nw[index_c]= cpow(kmin,-etam[index_c]) * (output_real_nw[index_c - Nmax/2] + _Complex_I * output_imag_nw[index_c - Nmax/2])/Nmaxd;
        }
    }
    
    cmsym_nw[0] = cmsym_nw[0]/2.;
    cmsym_nw[Nmax] = cmsym_nw[Nmax]/2.;
    
    FFT(input_real_w,input_imag_w,output_real_w,output_imag_w,Nmax,stepsize);
    
    double complex *cmsym_w;
    class_alloc(cmsym_w,(Nmax+1)*sizeof(complex double),pnlpt->error_message);
    
    index_c = 0;
    for (index_c=0; index_c< Nmax+1; index_c++){
        if (index_c < Nmax/2) {
            cmsym_w[index_c] = cpow(kmin,-etam[index_c]) * (output_real_w[Nmax/2 - index_c] - _Complex_I * output_imag_w[Nmax/2 - index_c])/Nmaxd;
        }
        else {
            cmsym_w[index_c] = cpow(kmin,-etam[index_c]) * (output_real_w[index_c - Nmax/2] + _Complex_I * output_imag_w[index_c - Nmax/2])/Nmaxd;
        }
    }
    
    cmsym_w[0] = cmsym_w[0]/2.;
    cmsym_w[Nmax] = cmsym_w[Nmax]/2.;
    
    
    free(input_real_nw);
    free(input_imag_nw);
    free(output_real_nw);
    free(output_imag_nw);
    
    free(input_real_w);
    free(input_imag_w);
    free(output_real_w);
    free(output_imag_w);
    
    // Matrix multiplication
    
    
    double *P13_mu0_dd;
    double *P13UV_mu0_dd;
    class_alloc(P13_mu0_dd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_mu0_dd,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu0_dd;
    class_alloc(P22_mu0_dd,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P13_mu2_dd;
    double *P13UV_mu2_dd;
    class_alloc(P13_mu2_dd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_mu2_dd,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P13_mu2_vd;
    double *P13UV_mu2_vd;
    class_alloc(P13_mu2_vd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_mu2_vd,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P13_mu4_vd;
    double *P13UV_mu4_vd;
    class_alloc(P13_mu4_vd,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_mu4_vd,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P13_mu4_vv;
    double *P13UV_mu4_vv;
    class_alloc(P13_mu4_vv,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_mu4_vv,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P13_mu6;
    double *P13UV_mu6;
    class_alloc(P13_mu6,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P13UV_mu6,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P22_mu2_vd;
    class_alloc(P22_mu2_vd,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu2_dd;
    class_alloc(P22_mu2_dd,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu4_vv;
    class_alloc(P22_mu4_vv,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu4_vd;
    class_alloc(P22_mu4_vd,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu4_dd;
    class_alloc(P22_mu4_dd,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu6_vv;
    class_alloc(P22_mu6_vv,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu6_vd;
    class_alloc(P22_mu6_vd,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu8;
    class_alloc(P22_mu8,Nmax*sizeof(double),pnlpt->error_message);
    
    
    count = 0;
    index_l=0;
    index_i=0;
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu1 = -0.5*etam[index_i];
            nu2 = -0.5*etam[index_l];
            nu12 = nu1+nu2;
            
            pnlpt->M22_oneline_mu2_vd_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*(-1.)*f*(7.*f*(-1.+2.*nu1)*(-1.+2.*nu2)*(6.+7.*nu12)-4.*(46.+13.*nu2+98.*nu1*nu1*nu1*nu2-63.*nu2*nu2+7.*nu1*nu1*(-9.-10.*nu2+28.*nu2*nu2)+nu1*(13.-138.*nu2-70.*nu2*nu2+98.*nu2*nu2*nu2)))/392.;
            
            pnlpt->M22_oneline_mu2_dd_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*f*(7.*f*(2.+2.*nu1*nu1*nu1-nu2-nu2*nu2+2.*nu2*nu2*nu2-nu1*nu1*(1.+2.*nu2)-nu1*(1.+2.*nu2+2.*nu2*nu2))+4.*(10.-nu2+14.*nu1*nu1*nu1*nu2-17.*nu2*nu2+nu1*nu1*(-17.+6.*nu2+28.*nu2*nu2)+nu1*(-1.-22.*nu2+6.*nu2*nu2+14.*nu2*nu2*nu2)))/56.;
            
            pnlpt->M22_oneline_mu4_vv_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*f*f*(147.*f*f*(-1.+2.*nu1)*(-1.+2.*nu2)-28.*f*(-1.+2.*nu1)*(-1.+2.*nu2)*(-2.+7.*nu12)+8.*(50.-9.*nu2+98*nu1*nu1*nu1*nu2-35.*nu2*nu2+7.*nu1*nu1*(-5.-18.*nu2+28*nu2*nu2)+nu1*(-9.-66.*nu2-126*nu2*nu2+98.*nu2*nu2*nu2)))/1568.;
            
            pnlpt->M22_oneline_mu4_vd_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*f*f*(58.+21.*nu2+112.*nu1*nu1*nu1*nu2-106.*nu2*nu2+2.*nu1*nu1*(-53.-6.*nu2+112.*nu2*nu2)+7.*f*(2.+nu1+4.*nu1*nu1*nu1+nu2-8.*nu1*nu2-8.*nu1*nu1*nu2-8.*nu1*nu2*nu2+4.*nu2*nu2*nu2)+nu1*(21.-204.*nu2-12.*nu2*nu2+112.*nu2*nu2*nu2))/56.;
            
            pnlpt->M22_oneline_mu4_dd_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*f*f*(2.*nu1-1.)*(2.*nu2-1.)*(2.+nu1*nu1+3.*nu2+nu2*nu2+nu1*(3.+2.*nu2))/8.;
            
            pnlpt->M22_oneline_mu6_vv_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*f*f*f*(7.*f*(1.+4.* nu1*nu1*nu1+nu1*nu1*(2.-12.*nu2)+2.*nu2+2.*nu2*nu2+4.*nu2*nu2*nu2-2.*nu1*(-1.+4.*nu2+6.*nu2*nu2))+2.*(26.+9.*nu2+56.*nu1*nu1*nu1*nu2-38.*nu2*nu2+2.*nu1*nu1*(-19.-18.*nu2+56.*nu2*nu2)+nu1*(9.-84.*nu2-36.*nu2*nu2+56.*nu2*nu2*nu2)))/112.;
            
            pnlpt->M22_oneline_mu6_vd_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*f*f*f*(2.*nu1-1.)*(2.*nu2-1.)*(2.+2.*nu1*nu1+5.*nu2+2.*nu2*nu2+nu1*(5.+4.*nu2))/8.;
            
            pnlpt->M22_oneline_mu8_complex[count] = (pnlpt->M22_oneline_complex[count]*196./(98.*nu1*nu2*nu12*nu12-91.*nu12*nu12+36.*nu1*nu2-14.*nu1*nu2*nu12+3.*nu12+58.))*f*f*f*f*(2.*nu1-1.)*(2.*nu2-1.)*(3.+4.*nu1*nu1+8.*nu2+4.*nu2*nu2+8.*nu1*(1.+nu2))/32.;
            count++;
        }
    }

    for (index_i=0; index_i<Nmax+1; index_i++){
        nu1=-0.5*etam[index_i];
        pnlpt->M13_mu2_dd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i]*2./(1.+9.*nu1)*9.*f*(1.+nu1);
        pnlpt->M13_mu2_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i]*2./(1.+9.*nu1)*(-1.)*(f*(7.+9.*f-9.*nu1));
        pnlpt->M13_mu4_vv_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i]/(1.+9.*nu1)*(-1.)*3.*f*f*(5.+6.*f-3.*nu1);
        pnlpt->M13_mu4_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i]*2./(1.+9.*nu1)*9.*f*f*(1.+2.*nu1);
        pnlpt->M13_mu6_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i]*2./(1.+9.*nu1)*9.*f*f*f*nu1;
    }
    
    double complex *f13_mu0_dd;
    class_alloc(f13_mu0_dd,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu0_dd;
    class_alloc(f22_mu0_dd,Nmax*sizeof(complex double),pnlpt->error_message);
    
    double complex *f13_mu2_vd;
    class_alloc(f13_mu2_vd,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f13_mu2_dd;
    class_alloc(f13_mu2_dd,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f13_mu4_vv;
    class_alloc(f13_mu4_vv,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f13_mu4_vd;
    class_alloc(f13_mu4_vd,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f13_mu6;
    class_alloc(f13_mu6,Nmax*sizeof(complex double),pnlpt->error_message);
    
    double complex *f22_mu2_vd;
    class_alloc(f22_mu2_vd,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu2_dd;
    class_alloc(f22_mu2_dd,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu4_vd;
    class_alloc(f22_mu4_vd,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu4_dd;
    class_alloc(f22_mu4_dd,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu4_vv;
    class_alloc(f22_mu4_vv,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu6_vv;
    class_alloc(f22_mu6_vv,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu6_vd;
    class_alloc(f22_mu6_vd,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu8;
    class_alloc(f22_mu8,Nmax*sizeof(complex double),pnlpt->error_message);
    
    count=0;
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym_nw[count]* cpow(kdisc[index_j], etam[count]);
        }
        
        f13_mu0_dd[index_j]=zdotu_(&Nmaxf, x, &inc, pnlpt->M13_oneline_complex, &inc);
        P13UV_mu0_dd[index_j] = -61. * Pnw[index_j] * pow(kdisc[index_j], 2.) * sigmav / 105.;
        P13_mu0_dd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu0_dd[index_j] * Pnw[index_j]) + P13UV_mu0_dd[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
    
        P13UV_mu2_dd[index_j] = -1.*Pnw[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*f*(105.*f-6.)/105.;
        P13UV_mu2_vd[index_j] = -1.*Pnw[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*f*(250.+144.*f)/105.;
        P13UV_mu4_vv[index_j] = -1.*Pnw[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*f*f*(63.+48.*f)/35.;
        P13UV_mu4_vd[index_j] = -1.*Pnw[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*f*f*(44.+70.*f)/35.;
        P13UV_mu6[index_j] = -1.*Pnw[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*f*f*f*(46.+35.*f)/35.;
        
        f13_mu2_dd[index_j]= zdotu_(&Nmaxf, x, &inc, pnlpt->M13_mu2_dd_oneline_complex, &inc);
        f13_mu2_vd[index_j]= zdotu_(&Nmaxf, x, &inc, pnlpt->M13_mu2_vd_oneline_complex, &inc);
        f13_mu4_vv[index_j]= zdotu_(&Nmaxf, x, &inc, pnlpt->M13_mu4_vv_oneline_complex, &inc);
        f13_mu4_vd[index_j]= zdotu_(&Nmaxf, x, &inc, pnlpt->M13_mu4_vd_oneline_complex, &inc);
        f13_mu6[index_j]= zdotu_(&Nmaxf, x, &inc, pnlpt->M13_mu6_oneline_complex, &inc);
        
        P13_mu2_dd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu2_dd[index_j] * Pnw[index_j]) + P13UV_mu2_dd[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
        P13_mu2_vd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu2_vd[index_j] * Pnw[index_j]) + P13UV_mu2_vd[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
        P13_mu4_vv[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu4_vv[index_j] * Pnw[index_j]) + P13UV_mu4_vv[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
        P13_mu4_vd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu4_vd[index_j] * Pnw[index_j]) + P13UV_mu4_vd[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
        P13_mu6[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu6[index_j] * Pnw[index_j]) + P13UV_mu6[index_j]) * exp(-pow(kdisc[index_j]/cutoff, 6.));
    }
    
    count=0;
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym_nw[count]* cpow(kdisc[index_j], etam[count]);
        }
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_complex, x, &inc, &beta, y, &inc);
        f22_mu0_dd[index_j]=zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu0_dd[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu0_dd[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu2_vd_complex, x, &inc, &beta, y, &inc);
        f22_mu2_vd[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu2_vd[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu2_vd[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu2_dd_complex, x, &inc, &beta, y, &inc);
        f22_mu2_dd[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu2_dd[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu2_dd[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu4_vv_complex, x, &inc, &beta, y, &inc);
        f22_mu4_vv[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu4_vv[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu4_vv[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu4_vd_complex, x, &inc, &beta, y, &inc);
        f22_mu4_vd[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu4_vd[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu4_vd[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu4_dd_complex, x, &inc, &beta, y, &inc);
        f22_mu4_dd[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu4_dd[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu4_dd[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu6_vv_complex, x, &inc, &beta, y, &inc);
        f22_mu6_vv[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu6_vv[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu6_vv[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu6_vd_complex, x, &inc, &beta, y, &inc);
        f22_mu6_vd[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu6_vd[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu6_vd[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu8_complex, x, &inc, &beta, y, &inc);
        f22_mu8[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu8[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu8[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        }
    
    // Matrix multiplication in case of Pw
    
    
    double *P13_mu0_dd_w;
    double *P22_mu0_dd_w;
    class_alloc(P13_mu0_dd_w,Nmax*sizeof(double),pnlpt->error_message);
    class_alloc(P22_mu0_dd_w,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P13_mu2_dd_w;
    class_alloc(P13_mu2_dd_w,Nmax*sizeof(double),pnlpt->error_message);
    double *P13_mu2_vd_w;
    class_alloc(P13_mu2_vd_w,Nmax*sizeof(double),pnlpt->error_message);
    double *P13_mu4_vd_w;
    class_alloc(P13_mu4_vd_w,Nmax*sizeof(double),pnlpt->error_message);
    double *P13_mu4_vv_w;
    class_alloc(P13_mu4_vv_w,Nmax*sizeof(double),pnlpt->error_message);
    double *P13_mu6_w;
    class_alloc(P13_mu6_w,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P22_mu2_vd_w;
    class_alloc(P22_mu2_vd_w,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu2_dd_w;
    class_alloc(P22_mu2_dd_w,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu4_vv_w;
    class_alloc(P22_mu4_vv_w,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu4_vd_w;
    class_alloc(P22_mu4_vd_w,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu4_dd_w;
    class_alloc(P22_mu4_dd_w,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu6_vv_w;
    class_alloc(P22_mu6_vv_w,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu6_vd_w;
    class_alloc(P22_mu6_vd_w,Nmax*sizeof(double),pnlpt->error_message);
    double *P22_mu8_w;
    class_alloc(P22_mu8_w,Nmax*sizeof(double),pnlpt->error_message);
    
    double complex *f13_mu0_dd_w;
    class_alloc(f13_mu0_dd_w,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu0_dd_w;
    class_alloc(f22_mu0_dd_w,Nmax*sizeof(complex double),pnlpt->error_message);
    
    double complex *f13_mu2_vd_w;
    class_alloc(f13_mu2_vd_w,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f13_mu2_dd_w;
    class_alloc(f13_mu2_dd_w,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f13_mu4_vv_w;
    class_alloc(f13_mu4_vv_w,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f13_mu4_vd_w;
    class_alloc(f13_mu4_vd_w,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f13_mu6_w;
    class_alloc(f13_mu6_w,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu2_vd_w;
    class_alloc(f22_mu2_vd_w,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu2_dd_w;
    class_alloc(f22_mu2_dd_w,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu4_vd_w;
    class_alloc(f22_mu4_vd_w,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu4_dd_w;
    class_alloc(f22_mu4_dd_w,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu4_vv_w;
    class_alloc(f22_mu4_vv_w,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu6_vv_w;
    class_alloc(f22_mu6_vv_w,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu6_vd_w;
    class_alloc(f22_mu6_vd_w,Nmax*sizeof(complex double),pnlpt->error_message);
    double complex *f22_mu8_w;
    class_alloc(f22_mu8_w,Nmax*sizeof(complex double),pnlpt->error_message);
    
    count=0;
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
         //   x[count]= cmsym_nw[count]* cpow(kdisc[index_j], etam[count]);
            x_w[count]= cmsym_w[count]* cpow(kdisc[index_j], etam[count]);
        }
        
        f13_mu2_dd_w[index_j]= zdotu_(&Nmaxf, x_w, &inc, pnlpt->M13_mu2_dd_oneline_complex, &inc);
        f13_mu2_vd_w[index_j]= zdotu_(&Nmaxf, x_w, &inc, pnlpt->M13_mu2_vd_oneline_complex, &inc);
        f13_mu4_vv_w[index_j]= zdotu_(&Nmaxf, x_w, &inc, pnlpt->M13_mu4_vv_oneline_complex, &inc);
        f13_mu4_vd_w[index_j]= zdotu_(&Nmaxf, x_w, &inc, pnlpt->M13_mu4_vd_oneline_complex, &inc);
        f13_mu6_w[index_j]= zdotu_(&Nmaxf, x_w, &inc, pnlpt->M13_mu6_oneline_complex, &inc);
        f13_mu0_dd_w[index_j]=zdotu_(&Nmaxf, x_w, &inc, pnlpt->M13_oneline_complex, &inc);

        P13_mu0_dd_w[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu0_dd_w[index_j] * Pnw[index_j])) * exp(-pow(kdisc[index_j]/cutoff, 6.));
        P13_mu2_dd_w[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu2_dd_w[index_j] * Pnw[index_j])) * exp(-pow(kdisc[index_j]/cutoff, 6.));
        P13_mu2_vd_w[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu2_vd_w[index_j] * Pnw[index_j])) * exp(-pow(kdisc[index_j]/cutoff, 6.));
        P13_mu4_vv_w[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu4_vv_w[index_j] * Pnw[index_j])) * exp(-pow(kdisc[index_j]/cutoff, 6.));
        P13_mu4_vd_w[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu4_vd_w[index_j] * Pnw[index_j])) * exp(-pow(kdisc[index_j]/cutoff, 6.));
        P13_mu6_w[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu6_w[index_j] * Pnw[index_j])) * exp(-pow(kdisc[index_j]/cutoff, 6.));
    }
    
    count=0;
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x[count]= cmsym_nw[count]* cpow(kdisc[index_j], etam[count]);
            x_w[count]= 2.*cmsym_w[count]* cpow(kdisc[index_j], etam[count]);
        }
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_complex, x_w, &inc, &beta, y, &inc);
        f22_mu0_dd_w[index_j]=zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu0_dd_w[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu0_dd_w[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu2_vd_complex, x_w, &inc, &beta, y, &inc);
        f22_mu2_vd_w[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu2_vd_w[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu2_vd_w[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu2_dd_complex, x_w, &inc, &beta, y, &inc);
        f22_mu2_dd_w[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu2_dd_w[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu2_dd_w[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu4_vv_complex, x_w, &inc, &beta, y, &inc);
        f22_mu4_vv_w[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu4_vv_w[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu4_vv_w[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu4_vd_complex, x_w, &inc, &beta, y, &inc);
        f22_mu4_vd_w[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu4_vd_w[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu4_vd_w[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu4_dd_complex, x_w, &inc, &beta, y, &inc);
        f22_mu4_dd_w[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu4_dd_w[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu4_dd_w[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu6_vv_complex, x_w, &inc, &beta, y, &inc);
        f22_mu6_vv_w[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu6_vv_w[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu6_vv_w[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu6_vd_complex, x_w, &inc, &beta, y, &inc);
        f22_mu6_vd_w[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu6_vd_w[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu6_vd_w[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
        
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu8_complex, x_w, &inc, &beta, y, &inc);
        f22_mu8_w[index_j]= zdotu_(&Nmaxf, x, &inc, y, &inc);
        P22_mu8_w[index_j]= creal(cpow(kdisc[index_j], 3.) * f22_mu8_w[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );
    }
    
    double *dd_P13_mu4_vv;
    class_alloc(dd_P13_mu4_vv,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P13_mu4_vv, 1, dd_P13_mu4_vv,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P13_mu4_vv_ap_out=0;
    
    double *dd_P22_mu4_vv;
    class_alloc(dd_P22_mu4_vv,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu4_vv, 1, dd_P22_mu4_vv,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu4_vv_ap_out=0;
    
    double *dd_P22_mu4_vv_w;
    class_alloc(dd_P22_mu4_vv_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu4_vv_w, 1, dd_P22_mu4_vv_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu4_vv_w_ap_out=0;
    
    double *dd_P13_mu4_vv_w;
    class_alloc(dd_P13_mu4_vv_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P13_mu4_vv_w, 1, dd_P13_mu4_vv_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P13_mu4_vv_w_ap_out=0;
    
    double *dd_P13_mu6;
    class_alloc(dd_P13_mu6,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P13_mu6, 1, dd_P13_mu6,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P13_mu6_ap_out=0;

    double *dd_P22_mu6_vv;
    class_alloc(dd_P22_mu6_vv,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu6_vv, 1, dd_P22_mu6_vv,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu6_vv_ap_out=0;
    
    double *dd_P22_mu6_vv_w;
    class_alloc(dd_P22_mu6_vv_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu6_vv_w, 1, dd_P22_mu6_vv_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu6_vv_w_ap_out=0;
    
    double *dd_P13_mu6_w;
    class_alloc(dd_P13_mu6_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P13_mu6_w, 1, dd_P13_mu6_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P13_mu6_w_ap_out=0;
    
    double *dd_P22_mu8;
    class_alloc(dd_P22_mu8,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu8, 1, dd_P22_mu8,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu8_ap_out=0;
    
    double *dd_P22_mu8_w;
    class_alloc(dd_P22_mu8_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu8_w, 1, dd_P22_mu8_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu8_w_ap_out=0;
    
    double *dd_P13_mu0_dd;
    class_alloc(dd_P13_mu0_dd,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P13_mu0_dd, 1, dd_P13_mu0_dd,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P13_mu0_dd_ap_out=0;

    double *dd_P22_mu0_dd;
    class_alloc(dd_P22_mu0_dd,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu0_dd, 1, dd_P22_mu0_dd,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu0_dd_ap_out=0;
    
    double *dd_P13_mu0_dd_w;
    class_alloc(dd_P13_mu0_dd_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P13_mu0_dd_w, 1, dd_P13_mu0_dd_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P13_mu0_dd_w_ap_out=0;
    
    double *dd_P22_mu0_dd_w;
    class_alloc(dd_P22_mu0_dd_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu0_dd_w, 1, dd_P22_mu0_dd_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu0_dd_w_ap_out=0;
    
    double *dd_P22_mu2_dd;
    class_alloc(dd_P22_mu2_dd,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu2_dd, 1, dd_P22_mu2_dd,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu2_dd_ap_out=0;
    
    double *dd_P13_mu2_dd;
    class_alloc(dd_P13_mu2_dd,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P13_mu2_dd, 1, dd_P13_mu2_dd,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P13_mu2_dd_ap_out=0;
    
    double *dd_P22_mu2_dd_w;
    class_alloc(dd_P22_mu2_dd_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu2_dd_w, 1, dd_P22_mu2_dd_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu2_dd_w_ap_out=0;
    
    double *dd_P13_mu2_dd_w;
    class_alloc(dd_P13_mu2_dd_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P13_mu2_dd_w, 1, dd_P13_mu2_dd_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P13_mu2_dd_w_ap_out=0;
    
    double *dd_P22_mu4_dd;
    class_alloc(dd_P22_mu4_dd,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu4_dd, 1, dd_P22_mu4_dd,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu4_dd_ap_out=0;
    
    double *dd_P22_mu4_dd_w;
    class_alloc(dd_P22_mu4_dd_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu4_dd_w, 1, dd_P22_mu4_dd_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu4_dd_w_ap_out=0;
    
    double *dd_P13_mu2_vd;
    class_alloc(dd_P13_mu2_vd,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P13_mu2_vd, 1, dd_P13_mu2_vd,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P13_mu2_vd_ap_out=0;
    
    double *dd_P22_mu2_vd;
    class_alloc(dd_P22_mu2_vd,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu2_vd, 1, dd_P22_mu2_vd,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu2_vd_ap_out=0;
    
    double *dd_P22_mu2_vd_w;
    class_alloc(dd_P22_mu2_vd_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu2_vd_w, 1, dd_P22_mu2_vd_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu2_vd_w_ap_out=0;
    
    double *dd_P13_mu2_vd_w;
    class_alloc(dd_P13_mu2_vd_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P13_mu2_vd_w, 1, dd_P13_mu2_vd_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P13_mu2_vd_w_ap_out=0;
    
    double *dd_P13_mu4_vd;
    class_alloc(dd_P13_mu4_vd,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P13_mu4_vd, 1, dd_P13_mu4_vd,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P13_mu4_vd_ap_out=0;
    
    double *dd_P22_mu4_vd;
    class_alloc(dd_P22_mu4_vd,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu4_vd, 1, dd_P22_mu4_vd,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu4_vd_ap_out=0;
    
    double *dd_P22_mu4_vd_w;
    class_alloc(dd_P22_mu4_vd_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu4_vd_w, 1, dd_P22_mu4_vd_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu4_vd_w_ap_out=0;
    
    double *dd_P13_mu4_vd_w;
    class_alloc(dd_P13_mu4_vd_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P13_mu4_vd_w, 1, dd_P13_mu4_vd_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P13_mu4_vd_w_ap_out=0;
    
    double *dd_P22_mu6_vd;
    class_alloc(dd_P22_mu6_vd,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu6_vd, 1, dd_P22_mu6_vd,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu6_vd_ap_out=0;
    
    double *dd_P22_mu6_vd_w;
    class_alloc(dd_P22_mu6_vd_w,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P22_mu6_vd_w, 1, dd_P22_mu6_vd_w,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P22_mu6_vd_w_ap_out=0;
    
   //     printf("Until now all good\n");

        // Numerical integration over mu

    double Sigmatot=0.;
    double LegendreP0=1.;
    double LegendreP2=0.;
    double LegendreP4 =0.;
    double P1loopvv=0.;
    double P1loopvd=0.;
    double P1loopdd=0.;
    double Exp = 0.;
    double p_tree=0.;
    double p_tree_vv=0.;
    double p_tree_vd=0.;
    double p_tree_dd=0.;
    double P13ratio=0.;
    double Pctr0 = 0.;
    double Pctr2 = 0.;
    double Pctr4 = 0.;
    double P1=0.;
    double P1b1=0.;
    double mu = 0.;
    double ktrue = 0.;
    double mutrue = 0.;
    last_index=0;
    
    int index_gauss2 = 0;

    double P1loopvv_ap_ir=0.;
    double P1loopvd_ap_ir=0.;
    double P1loopdd_ap_ir=0.;
    
    
//    printf("Dratio=%lf\n",Dratio);
//    printf("hratio=%lf\n",hratio);
    for (index_j=0; index_j < Nside; index_j++){
    
    P1loop_0_vv[index_j] = 0.;
    P1loop_0_dd[index_j] = 0.;
    P1loop_0_vd[index_j] = 0.;
    P1loop_2_vv[index_j] = 0.;
    P1loop_2_dd[index_j] = 0.;
    P1loop_2_vd[index_j] = 0.;
    P1loop_4_vv[index_j] = 0.;
    P1loop_4_dd[index_j] = 0.;
    P1loop_4_vd[index_j] = 0.;
    P_CTR_0[index_j] = 0.;
    P_CTR_2[index_j] = 0.;
    P_CTR_4[index_j] = 0.;
    Ptree_0_vv[index_j] = 0.;
    Ptree_0_vd[index_j] = 0.;
    Ptree_0_dd[index_j] = 0.;
    Ptree_2_vv[index_j] = 0.;
    Ptree_2_vd[index_j] = 0.;
    Ptree_4_vv[index_j] = 0.;
    
    P1loop_0_vv[Nmax-1-index_j] = 0.;
    P1loop_0_dd[Nmax-1-index_j] = 0.;
    P1loop_0_vd[Nmax-1-index_j] = 0.;
    P1loop_2_vv[Nmax-1-index_j] = 0.;
    P1loop_2_dd[Nmax-1-index_j] = 0.;
    P1loop_2_vd[Nmax-1-index_j] = 0.;
    P1loop_4_vv[Nmax-1-index_j] = 0.;
    P1loop_4_dd[Nmax-1-index_j] = 0.;
    P1loop_4_vd[Nmax-1-index_j] = 0.;
    P_CTR_0[Nmax-1-index_j] = 0.;
    P_CTR_2[Nmax-1-index_j] = 0.;
    P_CTR_4[Nmax-1-index_j] = 0.;
    Ptree_0_vv[Nmax-1-index_j] = 0.;
    Ptree_0_vd[Nmax-1-index_j] = 0.;
    Ptree_0_dd[Nmax-1-index_j] = 0.;
    Ptree_2_vv[Nmax-1-index_j] = 0.;
    Ptree_2_vd[Nmax-1-index_j] = 0.;
    Ptree_4_vv[Nmax-1-index_j] = 0.;
    }

//    for (index_j=0; index_j < Nmax; index_j++){
        for (index_j=Nside; index_j < Nmax-Nside; index_j++){

        P1loop_0_vv[index_j] = 0.;
        P1loop_0_dd[index_j] = 0.;
        P1loop_0_vd[index_j] = 0.;
        
        P1loop_2_vv[index_j] = 0.;
        P1loop_2_dd[index_j] = 0.;
        P1loop_2_vd[index_j] = 0.;
        
        P1loop_4_vv[index_j] = 0.;
        P1loop_4_dd[index_j] = 0.;
        P1loop_4_vd[index_j] = 0.;
        
        P_CTR_0[index_j] = 0.;
        P_CTR_2[index_j] = 0.;
        P_CTR_4[index_j] = 0.;
        
        P10b1[index_j] = 0.;
        P10[index_j] = 0.;
        P12[index_j] = 0.;
        
        Ptree_0_vv[index_j] = 0.;
        Ptree_0_vd[index_j] = 0.;
        Ptree_0_dd[index_j] = 0.;
        Ptree_2_vv[index_j] = 0.;
        Ptree_2_vd[index_j] = 0.;
        Ptree_4_vv[index_j] = 0.;
        
        for (index_gauss2=0; index_gauss2 < 40; index_gauss2++){
        
        mu = pnlpt->gauss_x[index_gauss2];
            
        if (pnlpt->AP_effect == AP_effect_yes){
        mutrue = mu*hratio/pow((1./Dratio/Dratio + (hratio*hratio - 1./Dratio/Dratio)*mu*mu),0.5);
        ktrue = kdisc[index_j]*pow((1./Dratio/Dratio + (hratio*hratio - 1./Dratio/Dratio)*mu*mu),0.5);
        }
            
        else {
            mutrue = mu;
            ktrue =kdisc[index_j];
        }

            
        class_call(array_interpolate_spline(kdisc,Nmax,Pnw,dd_Pnw,1,ktrue,&last_index,&Pnw_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
        class_call(array_interpolate_spline(kdisc,Nmax,Pw,dd_Pw,1,ktrue,&last_index,&Pw_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu4_vv,dd_P22_mu4_vv,1,ktrue,&last_index,&P22_mu4_vv_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P13_mu4_vv,dd_P13_mu4_vv,1,ktrue,&last_index,&P13_mu4_vv_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu4_vv_w,dd_P22_mu4_vv_w,1,ktrue,&last_index,&P22_mu4_vv_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P13_mu4_vv_w,dd_P13_mu4_vv_w,1,ktrue,&last_index,&P13_mu4_vv_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P13_mu6,dd_P13_mu6,1,ktrue,&last_index,&P13_mu6_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu6_vv,dd_P22_mu6_vv,1,ktrue,&last_index,&P22_mu6_vv_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu6_vv_w,dd_P22_mu6_vv_w,1,ktrue,&last_index,&P22_mu6_vv_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P13_mu6_w,dd_P13_mu6_w,1,ktrue,&last_index,&P13_mu6_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu8,dd_P22_mu8,1,ktrue,&last_index,&P22_mu8_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu8_w,dd_P22_mu8_w,1,ktrue,&last_index,&P22_mu8_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu0_dd,dd_P22_mu0_dd,1,ktrue,&last_index,&P22_mu0_dd_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P13_mu0_dd,dd_P13_mu0_dd,1,ktrue,&last_index,&P13_mu0_dd_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
            class_call(array_interpolate_spline(kdisc,Nmax,P13_mu0_dd_w,dd_P13_mu0_dd_w,1,ktrue,&last_index,&P13_mu0_dd_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu0_dd_w,dd_P22_mu0_dd_w,1,ktrue,&last_index,&P22_mu0_dd_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu2_dd,dd_P22_mu2_dd,1,ktrue,&last_index,&P22_mu2_dd_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
            class_call(array_interpolate_spline(kdisc,Nmax,P13_mu2_dd,dd_P13_mu2_dd,1,ktrue,&last_index,&P13_mu2_dd_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu2_dd_w,dd_P22_mu2_dd_w,1,ktrue,&last_index,&P22_mu2_dd_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
            class_call(array_interpolate_spline(kdisc,Nmax,P13_mu2_dd_w,dd_P13_mu2_dd_w,1,ktrue,&last_index,&P13_mu2_dd_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu4_dd,dd_P22_mu4_dd,1,ktrue,&last_index,&P22_mu4_dd_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu4_dd_w,dd_P22_mu4_dd_w,1,ktrue,&last_index,&P22_mu4_dd_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);

            class_call(array_interpolate_spline(kdisc,Nmax,P13_mu2_vd,dd_P13_mu2_vd,1,ktrue,&last_index,&P13_mu2_vd_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu2_vd,dd_P22_mu2_vd,1,ktrue,&last_index,&P22_mu2_vd_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu2_vd_w,dd_P22_mu2_vd_w,1,ktrue,&last_index,&P22_mu2_vd_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P13_mu2_vd_w,dd_P13_mu2_vd_w,1,ktrue,&last_index,&P13_mu2_vd_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P13_mu4_vd,dd_P13_mu4_vd,1,ktrue,&last_index,&P13_mu4_vd_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu4_vd,dd_P22_mu4_vd,1,ktrue,&last_index,&P22_mu4_vd_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu4_vd_w,dd_P22_mu4_vd_w,1,ktrue,&last_index,&P22_mu4_vd_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P13_mu4_vd_w,dd_P13_mu4_vd_w,1,ktrue,&last_index,&P13_mu4_vd_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu6_vd,dd_P22_mu6_vd,1,ktrue,&last_index,&P22_mu6_vd_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu6_vd_w,dd_P22_mu6_vd_w,1,ktrue,&last_index,&P22_mu6_vd_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P22_mu6_vd_w,dd_P22_mu6_vd_w,1,ktrue,&last_index,&P22_mu6_vd_w_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
        LegendreP2 = (3.*pow(mu,2.)-1.)/2.;
        LegendreP4 = (35.*pow(mu,4.)-30.*pow(mu,2.)+3.)/8.;
            
            
            Sigmatot = SigmaBAO*(1.+f*mutrue*mutrue*(2.+f))+f*f*mutrue*mutrue*(mutrue*mutrue-1.)*deltaSigmaBAO;
            Exp = exp(-Sigmatot * pow(ktrue,2.));
                        
            
            p_tree = (Pnw_ap_out + (1. + Sigmatot*pow(ktrue,2.))* Pw_ap_out * Exp);
            
            P13ratio = 1.+(Pw_ap_out/Pnw_ap_out)*Exp;
            P1b1 = (Pnw[index_j]+Pw[index_j]*Exp)*pnlpt->gauss_w[index_gauss2];
            P1 = (Pnw[index_j]+Pw[index_j]*Exp)*f*pow(pnlpt->gauss_x[index_gauss2],2.)*pnlpt->gauss_w[index_gauss2];
            
            p_tree_vv = (Pnw_ap_out + (1. + Sigmatot*pow(ktrue,2.))* Pw_ap_out * Exp)*pow(f*pow(mutrue,2.),2.)*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            p_tree_vd = (Pnw_ap_out + (1. + Sigmatot*pow(ktrue,2.))* Pw_ap_out * Exp)*2.*f*pow(mutrue,2.)*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            p_tree_dd = (Pnw_ap_out + (1. + Sigmatot*pow(ktrue,2.))* Pw_ap_out * Exp)*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            
            Pctr0 = pow(ktrue,2.)*(Pnw_ap_out+Pw_ap_out*Exp)*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            Pctr2 = (Pnw_ap_out+Pw_ap_out*Exp)*pnlpt->gauss_w[index_gauss2]*f*pow(mutrue*ktrue,2.)*hratio/Dratio/Dratio;
            Pctr4 = pow(ktrue,2.)*(Pnw_ap_out+Pw_ap_out*Exp)*pnlpt->gauss_w[index_gauss2]*f*f*pow(mutrue,4.)*hratio/Dratio/Dratio;
            
            P1loopvv = (p_tree*0. +(P13_mu4_vv_ap_out*P13ratio  + P22_mu4_vv_ap_out + (P22_mu4_vv_w_ap_out+P13_mu4_vv_w_ap_out)*Exp)*pow(mutrue,4.) + (P13_mu6_ap_out*P13ratio + P22_mu6_vv_ap_out+ (P22_mu6_vv_w_ap_out+P13_mu6_w_ap_out)*Exp)*pow(mutrue,6.)+(P22_mu8_ap_out+P22_mu8_w_ap_out*Exp)*pow(mutrue,8.))*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            
            P1loopdd = ((p_tree*0. + P22_mu0_dd_ap_out + P13_mu0_dd_ap_out*P13ratio + (P13_mu0_dd_w_ap_out+P22_mu0_dd_w_ap_out)*Exp)+(P22_mu2_dd_ap_out + P13_mu2_dd_ap_out*P13ratio + (P22_mu2_dd_w_ap_out+P13_mu2_dd_w_ap_out)*Exp)*pow(mutrue,2.)+(P22_mu4_dd_ap_out+P22_mu4_dd_w_ap_out*Exp)*pow(mutrue,4.))*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            
            P1loopvd = ((p_tree*2.*f*0. + P13_mu2_vd_ap_out*P13ratio + P22_mu2_vd_ap_out+ (P22_mu2_vd_w_ap_out+P13_mu2_vd_w_ap_out)*Exp)*pow(mutrue,2.)+(P13_mu4_vd_ap_out*P13ratio + P22_mu4_vd_ap_out + (P22_mu4_vd_w_ap_out+P13_mu4_vd_w_ap_out)*Exp)*pow(mutrue,4.)+(P22_mu6_vd_ap_out+P22_mu6_vd_w_ap_out*Exp)*pow(mutrue,6.))*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
                   

                   
        P1loopdd_ap_ir = ((p_tree+ P22_mu0_dd_ap_out + P13_mu0_dd_ap_out*P13ratio + (P13_mu0_dd_w_ap_out+P22_mu0_dd_w_ap_out)*Exp)+(P22_mu2_dd_ap_out + P13_mu2_dd_ap_out*P13ratio + (P22_mu2_dd_w_ap_out+P13_mu2_dd_w_ap_out)*Exp)*pow(mutrue,2.)+(P22_mu4_dd_ap_out+P22_mu4_dd_w_ap_out*Exp)*pow(mutrue,4.))*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
        P1loopvd_ap_ir = ((p_tree*2.*f + P13_mu2_vd_ap_out*P13ratio + P22_mu2_vd_ap_out+ (P22_mu2_vd_w_ap_out+P13_mu2_vd_w_ap_out)*Exp)*pow(mutrue,2.)+(P13_mu4_vd_ap_out*P13ratio + P22_mu4_vd_ap_out + (P22_mu4_vd_w_ap_out+P13_mu4_vd_w_ap_out)*Exp)*pow(mutrue,4.)+(P22_mu6_vd_ap_out+P22_mu6_vd_w_ap_out*Exp)*pow(mutrue,6.))*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
        P1loop_0_vv[index_j] +=  P1loopvv*LegendreP0/2.;
        P1loop_2_vv[index_j] +=  P1loopvv*LegendreP2*2.5;
        P1loop_4_vv[index_j] +=  P1loopvv*LegendreP4*4.5;
            
        P1loop_0_dd[index_j] +=  P1loopdd*LegendreP0/2.;
        P1loop_2_dd[index_j] +=  P1loopdd_ap_ir*LegendreP2*2.5;
        P1loop_4_dd[index_j] +=  P1loopdd_ap_ir*LegendreP4*4.5;

        P1loop_0_vd[index_j] +=  P1loopvd*LegendreP0/2.;
        P1loop_2_vd[index_j] +=  P1loopvd*LegendreP2*2.5;
        P1loop_4_vd[index_j] +=  P1loopvd_ap_ir*LegendreP4*4.5;
                                
            /* this was used before the breakdown into tree and one-loop
        P1loop_0_vv[index_j] +=  P1loopvv*LegendreP0/2.;
        P1loop_2_vv[index_j] +=  P1loopvv*LegendreP2*2.5;
        P1loop_4_vv[index_j] +=  P1loopvv*LegendreP4*4.5;
            
        P1loop_0_dd[index_j] +=  P1loopdd*LegendreP0/2.;
        P1loop_2_dd[index_j] +=  P1loopdd*LegendreP2*2.5;
        P1loop_4_dd[index_j] +=  P1loopdd*LegendreP4*4.5;

        P1loop_0_vd[index_j] +=  P1loopvd*LegendreP0/2.;
        P1loop_2_vd[index_j] +=  P1loopvd*LegendreP2*2.5;
        P1loop_4_vd[index_j] +=  P1loopvd*LegendreP4*4.5;
        */
            
        P_CTR_0[index_j] += Pctr0*LegendreP0/2.;
        P_CTR_2[index_j] += Pctr2*LegendreP2*2.5;
        P_CTR_4[index_j] += Pctr4*LegendreP4*4.5;
            
            Ptree_0_vv[index_j] += p_tree_vv*LegendreP0/2.;
            Ptree_0_vd[index_j] += p_tree_vd*LegendreP0/2.;
            Ptree_0_dd[index_j] += p_tree_dd*LegendreP0/2.;
            Ptree_2_vv[index_j] += p_tree_vv*LegendreP2*2.5;
            Ptree_2_vd[index_j] += p_tree_vd*LegendreP2*2.5;
            Ptree_4_vv[index_j] += p_tree_vv*LegendreP4*4.5;
            
        P10b1[index_j] +=  P1b1*LegendreP0/2.;
        P10[index_j] +=  P1*LegendreP0/2.;
        P12[index_j] +=  P1*LegendreP2*2.5;
            
           /* P1loopvv = ((Pbin[index_j]*f*f + P13_mu4_vv[index_j]+P22_mu4_vv[index_j])*pow(pnlpt->gauss_x[index_gauss2],4.) + (P13_mu6[index_j]+P22_mu6_vv[index_j])*pow(pnlpt->gauss_x[index_gauss2],6.)+P22_mu8[index_j]*pow(pnlpt->gauss_x[index_gauss2],8.))*pnlpt->gauss_w[index_gauss2];
            P1loopdd = (Pbin[index_j] + P13[index_j] + P22[index_j] + (P13_mu2_dd[index_j] + P22_mu2_dd[index_j])*pow(pnlpt->gauss_x[index_gauss2],2.) + P22_mu4_dd[index_j]*pow(pnlpt->gauss_x[index_gauss2],4.))*pnlpt->gauss_w[index_gauss2];
            P1loopvd = (((Pbin[index_j])*2.*f + P13_mu2_vd[index_j] + P22_mu2_vd[index_j])*pow(pnlpt->gauss_x[index_gauss2],2.)+(P13_mu4_vd[index_j]+P22_mu4_vd[index_j])*pow(pnlpt->gauss_x[index_gauss2],4.)+P22_mu6_vd[index_j]*pow(pnlpt->gauss_x[index_gauss2],6.))*pnlpt->gauss_w[index_gauss2];*/
           //  P1loopvv = (((Pnw[index_j] + (1. + Sigmatot*kdisc[index_j]*kdisc[index_j])* Pw[index_j] * exp(-Sigmatot * pow(kdisc[index_j],2.)))*f*f + P13_mu4_vv[index_j]+P22_mu4_vv[index_j])*pow(gauss_x[index_gauss],4.) + (P13_mu6[index_j]+P22_mu6_vv[index_j])*pow(gauss_x[index_gauss],6.)+P22_mu8[index_j]*pow(gauss_x[index_gauss],8.))*gauss_w[index_gauss];
          //   P1loopdd = ((Pnw[index_j] + (1. + Sigmatot*kdisc[index_j]*kdisc[index_j])* Pw[index_j] * exp(-Sigmatot * pow(kdisc[index_j],2.))) + P13[index_j] + P22[index_j] + (P13_mu2_dd[index_j] + P22_mu2_dd[index_j])*pow(gauss_x[index_gauss],2.) + P22_mu4_dd[index_j]*pow(gauss_x[index_gauss],4.))*gauss_w[index_gauss];
          //   P1loopvd = (((Pnw[index_j] + (1. + Sigmatot*kdisc[index_j]*kdisc[index_j])* Pw[index_j] * exp(-Sigmatot * pow(kdisc[index_j],2.)))*2.*f + P13_mu2_vd[index_j] + P22_mu2_vd[index_j])*pow(gauss_x[index_gauss],2.)+(P13_mu4_vd[index_j]+P22_mu4_vd[index_j])*pow(gauss_x[index_gauss],4.)+P22_mu6_vd[index_j]*pow(gauss_x[index_gauss],6.))*gauss_w[index_gauss];
        }
        
  //      printf("%lf %lf %lf %lf\n",kdisc[index_j],P1loop_4_vv[index_j],P1loop_4_vd[index_j],P1loop_4_dd[index_j]);
        
//        printf("%lf %lf %lf %lf\n",kdisc[index_j],P_CTR_0[index_j],P_CTR_2[index_j],P_CTR_4[index_j]);
        
        P1loop_4_dd[index_j] =  P1loop_4_dd[index_j];
        P1loop_4_vd[index_j] =  P1loop_4_vd[index_j];
        P1loop_2_dd[index_j] =  P1loop_2_dd[index_j];
    }
    

    free(dd_P13_mu4_vv);
    free(dd_P22_mu4_vv);
    free(dd_P22_mu4_vv_w);
    free(dd_P13_mu4_vv_w);
    free(dd_P13_mu6);
    free(dd_P22_mu6_vv);
    free(dd_P22_mu6_vv_w);
    free(dd_P13_mu6_w);
    free(dd_P22_mu8);
    free(dd_P22_mu8_w);
    free(dd_P13_mu0_dd);
    free(dd_P22_mu0_dd);
    free(dd_P13_mu0_dd_w);
    free(dd_P22_mu0_dd_w);
    free(dd_P22_mu2_dd);
    free(dd_P13_mu2_dd);
    free(dd_P22_mu2_dd_w);
    free(dd_P13_mu2_dd_w);
    free(dd_P22_mu4_dd);
    free(dd_P22_mu4_dd_w);
    free(dd_P13_mu2_vd);
    free(dd_P22_mu2_vd);
    free(dd_P22_mu2_vd_w);
    free(dd_P13_mu2_vd_w);
    free(dd_P13_mu4_vd);
    free(dd_P22_mu4_vd);
    free(dd_P22_mu4_vd_w);
    free(dd_P13_mu4_vd_w);
    free(dd_P22_mu6_vd);
    free(dd_P22_mu6_vd_w);
    

    free(f22_mu0_dd);
    free(f22_mu0_dd_w);
    free(f13_mu0_dd);
    free(f13_mu0_dd_w);
    free(P13_mu0_dd_w);
    free(P13_mu0_dd);
    free(P13UV_mu0_dd);
    free(P22_mu0_dd);
    free(P22_mu0_dd_w);
    
    free(f22_mu8);
    free(f22_mu6_vd);
    free(f22_mu6_vv);
    free(f22_mu4_vv);
    free(f22_mu4_vd);
    free(f22_mu4_dd);
    free(f22_mu2_vd);
    free(f22_mu2_dd);
    
    free(P22_mu8);
    free(P22_mu6_vd);
    free(P22_mu6_vv);
    free(P22_mu4_vv);
    free(P22_mu4_vd);
    free(P22_mu4_dd);
    free(P22_mu2_vd);
    free(P22_mu2_dd);
    
    free(f13_mu4_vd);
    free(f13_mu4_vv);
    free(f13_mu2_vd);
    free(f13_mu2_dd);
    free(f13_mu6);
    
    free(P13_mu4_vd);
    free(P13_mu4_vv);
    free(P13_mu2_vd);
    free(P13_mu2_dd);
    free(P13_mu6);
    
    free(P13UV_mu4_vd);
    free(P13UV_mu4_vv);
    free(P13UV_mu2_vd);
    free(P13UV_mu2_dd);
    free(P13UV_mu6);
    
    free(f22_mu8_w);
    free(f22_mu6_vd_w);
    free(f22_mu6_vv_w);
    free(f22_mu4_vv_w);
    free(f22_mu4_vd_w);
    free(f22_mu4_dd_w);
    free(f22_mu2_vd_w);
    free(f22_mu2_dd_w);
    
    free(P22_mu8_w);
    free(P22_mu6_vd_w);
    free(P22_mu6_vv_w);
    free(P22_mu4_vv_w);
    free(P22_mu4_vd_w);
    free(P22_mu4_dd_w);
    free(P22_mu2_vd_w);
    free(P22_mu2_dd_w);
    
    free(f13_mu4_vd_w);
    free(f13_mu4_vv_w);
    free(f13_mu2_vd_w);
    free(f13_mu2_dd_w);
    free(f13_mu6_w);
    
    free(P13_mu4_vd_w);
    free(P13_mu4_vv_w);
    free(P13_mu2_vd_w);
    free(P13_mu2_dd_w);
    free(P13_mu6_w);
    
    free(cmsym_w);
    free(cmsym_nw);
    
    
//    printf("%lf\n",SigmaBAO);
    
    }// end of second IR resummation condition
    
    
// Constructing the final output spectra
    
    double *ddpk_nl_0_vv;
    class_alloc(ddpk_nl_0_vv,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P1loop_0_vv,
                                          1,
                                          ddpk_nl_0_vv,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    last_index=0;
    pk_nl_out=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P1loop_0_vv,
                                                ddpk_nl_0_vv,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_nl_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            pk_l_0_vv[index_k] = pk_nl_out + 1.e7;
        }
        else{
            
            pk_l_0_vv[index_k] =  -1.*exp(lnpk_l[index_k]+2.*lnk_l[index_k])*sigmav*(f*f*(441.+566.*f+175.*f*f)/1225.) + 1.e7;
          //  pk_l_0_vv[index_k] = exp(lnpk_l[index_k])*(f*f/5.);
            
         //   -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*(f*f*(441.+566.*f+175.*f*f)/1225.)
            
        }
        //     printf("%le %le\n",pnlpt->k[index_k]/pba->h,pk_l_0_vv[index_k]*pow(pba->h,3));
    }
    
    double *ddpk_nl_0_vd;
    class_alloc(ddpk_nl_0_vd,sizeof(double)*Nmax,pnlpt->error_message);
    
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P1loop_0_vd,
                                          1,
                                          ddpk_nl_0_vd,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    last_index=0;
    pk_nl_out=0;
    
    
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P1loop_0_vd,
                                                ddpk_nl_0_vd,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_nl_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            
            
            pk_l_0_vd[index_k] = pk_nl_out + 1.e7;
        }
        else {
      //      pk_l_0_vd[index_k] = (f*2./3.)*exp(lnpk_l[index_k]);
            pk_l_0_vd[index_k] = -1.*exp(lnpk_l[index_k]+2.*lnk_l[index_k])*sigmav*(2.*f*(625. + 558.*f + 315.*f*f)/1575.) + 1.e7;
                    }
        //         printf("%le %le\n",pnlpt->k[index_k]/pba->h,pk_l_0_vd[index_k]*pow(pba->h,3));
    }
    
    
    double *ddpk_nl_0_dd;
    class_alloc(ddpk_nl_0_dd,sizeof(double)*Nmax,pnlpt->error_message);
    
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P1loop_0_dd,
                                          1,
                                          ddpk_nl_0_dd,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    last_index=0;
    pk_nl_out=0;
    
    
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P1loop_0_dd,
                                                ddpk_nl_0_dd,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_nl_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            
            pk_l_0_dd[index_k] = pk_nl_out + 1.e7;
        }
        else {
          //  pk_l_0_dd[index_k] = exp(lnpk_l[index_k]);
            pk_l_0_dd[index_k] = -1.*exp(lnpk_l[index_k]+2.*lnk_l[index_k])*sigmav*((61. -2.*f + 35.*f*f)/105.) + 1.e7;
        }
      //  -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*((61. -2.*f + 35.*f*f)/105.);
     //     printf("%le %le\n",pnlpt->k[index_k]/pba->h,pk_l_0_dd[index_k]*pow(pba->h,3));
    }
    
    
    double *ddpk_nl_2_vv;
    class_alloc(ddpk_nl_2_vv,sizeof(double)*Nmax,pnlpt->error_message);
    
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P1loop_2_vv,
                                          1,
                                          ddpk_nl_2_vv,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    last_index=0;
    pk_nl_out=0;
    
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P1loop_2_vv,
                                                ddpk_nl_2_vv,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_nl_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            
            pk_l_2_vv[index_k] = pk_nl_out + 1.e7;
        }
        
        else{
            pk_l_2_vv[index_k] = -1.*exp(lnpk_l[index_k]+2.*lnk_l[index_k])*sigmav*(2.*f*f*(54.+74.*f+25.*f*f)/105.) + 1.e7;
            //pk_l_2_vv[index_k] = exp(lnpk_l[index_k])*(f*f*4./7.);
           // -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*(2.*f*f*(54.+74.*f+25.*f*f)/105.);
        }
    }
    
    
    double *ddpk_nl_2_vd;
    class_alloc(ddpk_nl_2_vd,sizeof(double)*Nmax,pnlpt->error_message);
    
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P1loop_2_vd,
                                          1,
                                          ddpk_nl_2_vd,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    last_index=0;
    pk_nl_out=0;
    
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P1loop_2_vd,
                                                ddpk_nl_2_vd,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_nl_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            
            pk_l_2_vd[index_k] = pk_nl_out + 1.e7;
        }
        
        else{
            pk_l_2_vd[index_k] = -1.*exp(lnpk_l[index_k]+2.*lnk_l[index_k])*sigmav*4.*f*(175.+180.*f+126.*f*f)/441. + 1.e7;
        }
        //    -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*4.*f*(175.+180.*f+126.*f*f)/441.;
    }
    
    double *ddpk_nl_2_dd;
    class_alloc(ddpk_nl_2_dd,sizeof(double)*Nmax,pnlpt->error_message);
    
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P1loop_2_dd,
                                          1,
                                          ddpk_nl_2_dd,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    last_index=0;
    pk_nl_out=0;
    
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P1loop_2_dd,
                                                ddpk_nl_2_dd,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_nl_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            
            pk_l_2_dd[index_k] = pk_nl_out + 1.e7;
        }
        
        else{
            pk_l_2_dd[index_k] =  -1.*exp(lnpk_l[index_k]+2.*lnk_l[index_k])*sigmav*(2.*f*(35.*f-2.)/105.) + 1.e7;
        }
        
        //    -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*(2.*f*(35.*f-2.)/105.);
    }
    
    
    double *ddpk_nl_4_vv;
    class_alloc(ddpk_nl_4_vv,sizeof(double)*Nmax,pnlpt->error_message);
    
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P1loop_4_vv,
                                          1,
                                          ddpk_nl_4_vv,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    last_index=0;
    pk_nl_out=0;
    
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P1loop_4_vv,
                                                ddpk_nl_4_vv,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_nl_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            
            pk_l_4_vv[index_k] = pk_nl_out + 1.e7;
        }
        
        else{
            pk_l_4_vv[index_k] =  -1.*exp(lnpk_l[index_k]+2.*lnk_l[index_k])*sigmav*(24.*f*f*(33.+58.*f+25.*f*f)/1925.) + 1.e7;
            //-1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*(24.*f*f*(33.+58.*f+25.*f*f)/1925.);
        }
    }
    
    
    double *ddpk_nl_4_vd;
    class_alloc(ddpk_nl_4_vd,sizeof(double)*Nmax,pnlpt->error_message);
    
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P1loop_4_vd,
                                          1,
                                          ddpk_nl_4_vd,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    last_index=0;
    pk_nl_out=0;
    
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P1loop_4_vd,
                                                ddpk_nl_4_vd,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_nl_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            
            pk_l_4_vd[index_k] = pk_nl_out + 1.e7;
        }
        
        else{
            pk_l_4_vd[index_k] = -1.*exp(lnpk_l[index_k]+2.*lnk_l[index_k])*sigmav*16.*f*f*(22.+35.*f)/1225.+ 1.e7;
        }
        //    -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*16.*f*f*(22.+35.*f)/1225.;
    }
    
    double *ddpk_nl_4_dd;
    class_alloc(ddpk_nl_4_dd,sizeof(double)*Nmax,pnlpt->error_message);
    
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P1loop_4_dd,
                                          1,
                                          ddpk_nl_4_dd,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    
    last_index=0;
    pk_nl_out=0;
    double ir4dd;
    class_call(array_interpolate_spline(kdisc,
                                        Nmax,
                                        P1loop_4_dd,
                                        ddpk_nl_4_dd,
                                        1,
                                        kminnew,
                                        &last_index,
                                        &ir4dd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    last_index=0;
    pk_nl_out=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P1loop_4_dd,
                                                ddpk_nl_4_dd,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_nl_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            
            pk_l_4_dd[index_k] = pk_nl_out*exp(-pow(pnlpt->k[index_k]/3.,4.))+1.e7;
        }
        
        else{
            pk_l_4_dd[index_k] = 1.e7 + ir4dd*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.));
        }
    }
    
    double *ddpk_CTR_0;
    class_alloc(ddpk_CTR_0,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P_CTR_0,
                                          1,
                                          ddpk_CTR_0,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    last_index=0;
    pk_nl_out=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_CTR_0,
                                                ddpk_CTR_0,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_nl_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            
            pk_CTR_0[index_k] = pk_nl_out;
        }
        else{
            pk_CTR_0[index_k] = exp(lnpk_l[index_k]+2.*lnk_l[index_k]);
        }
    }
    
    double *ddpk_CTR_2;
    class_alloc(ddpk_CTR_2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P_CTR_2,
                                          1,
                                          ddpk_CTR_2,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    last_index=0;
    pk_nl_out=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_CTR_2,
                                                ddpk_CTR_2,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_nl_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            pk_CTR_2[index_k] = pk_nl_out;
        }
        else{
            pk_CTR_2[index_k] = exp(lnpk_l[index_k]+2.*lnk_l[index_k])*f*2./3.;
        //    pk_CTR_2[index_k] = exp(lnpk_l[index_k]+2.*lnk_l[index_k])*(f*f*2./3.*2.2*2.2+8.*2.2*pow(f,3.)/7.+10.*pow(f,4.)/21.);
        }
    }
    
    double *ddpk_CTR_4;
    class_alloc(ddpk_CTR_4,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_CTR_4,1,ddpk_CTR_4,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    pk_nl_out=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_CTR_4,
                                                ddpk_CTR_4,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_nl_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            pk_CTR_4[index_k] = pk_nl_out;
        }
        
        else{
            pk_CTR_4[index_k] = exp(lnpk_l[index_k]+2.*lnk_l[index_k])*f*f*8./35.;
        }
    }
    
    double *ddpk_tree_0_vv;
    class_alloc(ddpk_tree_0_vv,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,Ptree_0_vv,1,ddpk_tree_0_vv,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    pk_nl_out=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,Ptree_0_vv,ddpk_tree_0_vv,1,pnlpt->k[index_k],&last_index,&pk_nl_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_Tree_0_vv[index_k] = pk_nl_out + 1.e7;
        }
        else{
            pk_Tree_0_vv[index_k] = exp(lnpk_l[index_k])*f*f/5. + 1.e7;
        }
    }
    free(ddpk_tree_0_vv);
    free(Ptree_0_vv);
    
    double *ddpk_tree_0_vd;
    class_alloc(ddpk_tree_0_vd,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,Ptree_0_vd,1,ddpk_tree_0_vd,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    pk_nl_out=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,Ptree_0_vd,ddpk_tree_0_vd,1,pnlpt->k[index_k],&last_index,&pk_nl_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_Tree_0_vd[index_k] = pk_nl_out + 1.e7;
        }
        else{
            pk_Tree_0_vd[index_k] = exp(lnpk_l[index_k])*2.*f/3. + 1.e7;
        }
    }
    free(ddpk_tree_0_vd);
    free(Ptree_0_vd);
    
    double *ddpk_tree_0_dd;
    class_alloc(ddpk_tree_0_dd,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,Ptree_0_dd,1,ddpk_tree_0_dd,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    pk_nl_out=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,Ptree_0_dd,ddpk_tree_0_dd,1,pnlpt->k[index_k],&last_index,&pk_nl_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_Tree_0_dd[index_k] = pk_nl_out + 1.e7;
        }
        else{
            pk_Tree_0_dd[index_k] = exp(lnpk_l[index_k]) + 1.e7;
        }
    }
    free(ddpk_tree_0_dd);
    free(Ptree_0_dd);
    
    double *ddpk_tree_2_vv;
    class_alloc(ddpk_tree_2_vv,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,Ptree_2_vv,1,ddpk_tree_2_vv,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    pk_nl_out=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,Ptree_2_vv,ddpk_tree_2_vv,1,pnlpt->k[index_k],&last_index,&pk_nl_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_Tree_2_vv[index_k] = pk_nl_out + 1.e7;
        }
        else{
            pk_Tree_2_vv[index_k] = exp(lnpk_l[index_k])*4.*f*f/7. + 1.e7;
        }
    }
    free(ddpk_tree_2_vv);
    free(Ptree_2_vv);
    
    double *ddpk_tree_2_vd;
    class_alloc(ddpk_tree_2_vd,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,Ptree_2_vd,1,ddpk_tree_2_vd,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    pk_nl_out=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,Ptree_2_vd,ddpk_tree_2_vd,1,pnlpt->k[index_k],&last_index,&pk_nl_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_Tree_2_vd[index_k] = pk_nl_out + 1.e7;
        }
        else{
            pk_Tree_2_vd[index_k] = exp(lnpk_l[index_k])*4.*f/3. + 1.e7;
        }
    }
    free(ddpk_tree_2_vd);
    free(Ptree_2_vd);
    
    double *ddpk_tree_4_vv;
    class_alloc(ddpk_tree_4_vv,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,Ptree_4_vv,1,ddpk_tree_4_vv,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    pk_nl_out=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,Ptree_4_vv,ddpk_tree_4_vv,1,pnlpt->k[index_k],&last_index,&pk_nl_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_Tree_4_vv[index_k] = pk_nl_out + 1.e7;
        }
        else{
            pk_Tree_4_vv[index_k] = exp(lnpk_l[index_k])*8.*f*f/35. + 1.e7;
        }
    }
    free(ddpk_tree_4_vv);
    free(Ptree_4_vv);

    
    free(ddpk_nl_2_vv);
    free(ddpk_nl_2_vd);
    free(ddpk_nl_2_dd);
    
    free(ddpk_nl_4_vv);
    free(ddpk_nl_4_vd);
    free(ddpk_nl_4_dd);
    
    free(ddpk_nl_0_vv);
    free(ddpk_nl_0_vd);
    free(ddpk_nl_0_dd);
    
        free(P13UV_0_vv);
        free(P13_0_vv);
        free(P1loop_0_vv);
        free(P22_0_vv);
        free(P13UV_0_vd);
        free(P13_0_vd);
        free(P1loop_0_vd);
        free(P22_0_vd);
        free(P13UV_0_dd);
        free(P13_0_dd);
        free(P1loop_0_dd);
        free(P22_0_dd);
        free(P13UV_2_vv);
        free(P13_2_vv);
        free(P1loop_2_vv);
        free(P22_2_vv);
        free(P13UV_2_vd);
        free(P13_2_vd);
        free(P1loop_2_vd);
        free(P22_2_vd);
        free(P13UV_2_dd);
        free(P13_2_dd);
        free(P1loop_2_dd);
        free(P22_2_dd);
        free(P13UV_4_vv);
        free(P13_4_vv);
        free(P1loop_4_vv);
        free(P22_4_vv);
        free(P13UV_4_vd);
        free(P13_4_vd);
        free(P1loop_4_vd);
        free(P22_4_vd);
        free(P1loop_4_dd);
        free(P22_4_dd);
    
        free(P_CTR_0);
        free(P_CTR_2);
        free(P_CTR_4);
    
    free(ddpk_CTR_0);
    free(ddpk_CTR_2);
    free(ddpk_CTR_4);
    
    /*
     for (index_k=0; index_k < pnlpt->k_size; index_k++){
     pk_l_0_vv[index_k] = (f*f/5.)*exp(lnpk_l[index_k]);
     pk_l_0_vd[index_k] = (f*2./3.)*exp(lnpk_l[index_k]);
     pk_l_0_dd[index_k] = exp(lnpk_l[index_k]);
     
     pk_l_2_vv[index_k] = (f*f*4./7.)*exp(lnpk_l[index_k]);
     pk_l_2_vd[index_k] = (f*4./3.)*exp(lnpk_l[index_k]);
     pk_l_2_dd[index_k] = (f*4./3.)*exp(lnpk_l[index_k]);
     
     pk_l_4_vv[index_k] = (f*f*8./35.)*exp(lnpk_l[index_k]);
     }
     */
    
    
}// end of RSD conditional expression
    
     free(f13);
     free(x);
     free(x_w);
     free(f22);
     free(Pdisc);
     free(ddpk_nl);
     free(ddpk_CTR);
     free(ddpk_Tree);
     free(etam);
     free(cmsym);
     free(P22);
     free(P13);
     free(P13UV);
     free(P_CTR);
     free(P1loop);
     free(Ptree);
     free(myddlnpk);

     /* Computing the power spectra for biased tracers. For this reason we have to compute the FFTLog coefficients for a new 'bias' exponent b2 */
     
    if (pnlpt->bias == bias_yes) {
     
    if (pnlpt->nonlinear_pt_verbose > 0)
     printf("Computing the spectra for biased tracers...\n");

     double complex *etam2;
     class_alloc(etam2,(Nmax+1)*sizeof(complex double),pnlpt->error_message);
     
     double b2 = -1.6000001;
     int index_c2 = 0;
     for (index_c2=0; index_c2< Nmax +1 ; index_c2++){
         etam2[index_c2] = b2 + 2. * M_PI * _Complex_I * js[index_c2]/Nmaxd / Delta ;
     }

        double *input_real_bias;
        class_alloc(input_real_bias,(Nmax)*sizeof(double),pnlpt->error_message);
        double *input_imag_bias;
        class_alloc(input_imag_bias,(Nmax)*sizeof(double),pnlpt->error_message);
        double *output_real_bias;
        class_alloc(output_real_bias,(Nmax)*sizeof(double),pnlpt->error_message);
        double *output_imag_bias;
        class_alloc(output_imag_bias,(Nmax)*sizeof(double),pnlpt->error_message);
        
        index_c2=0;
        for (index_c2=0; index_c2< Nmax ; index_c2++){
            input_real_bias[index_c2] = Pbin[index_c2]* exp(-1.* index_c2 * b2* Delta);
            input_imag_bias[index_c2] = 0.;
        }
        
        FFT(input_real_bias,input_imag_bias,output_real_bias,output_imag_bias,Nmax,stepsize);

        double complex *cmsym2;
        class_alloc(cmsym2,(Nmax+1)*sizeof(complex double),pnlpt->error_message);
        
        //     double complex cmsym[Nmax+1];
        index_c2 = 0;
        for (index_c2=0; index_c2< Nmax+1; index_c2++){
            if (index_c2 < Nmax/2) {
                cmsym2[index_c2]= cpow(kmin,-etam2[index_c2]) * (output_real_bias[Nmax/2 - index_c2] - _Complex_I * output_imag_bias[Nmax/2 - index_c2])/Nmaxd;
            }
            else {
                cmsym2[index_c2]= cpow(kmin,-etam2[index_c2]) * (output_real_bias[index_c2 - Nmax/2] + _Complex_I * output_imag_bias[index_c2 - Nmax/2])/Nmaxd;
            }
        }
        
        cmsym2[0] = cmsym2[0]/2.;
        cmsym2[Nmax] = cmsym2[Nmax]/2.;

        free(input_real_bias);
        free(input_imag_bias);
        free(output_real_bias);
        free(output_imag_bias);
        
        
        
          /* Computing Id2d2 */
        
     double complex *f22_Id2d2;
     double *P_Id2d2;
     class_alloc(f22_Id2d2,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P_Id2d2,Nmax*sizeof(double),pnlpt->error_message);
     
     double epsilon_for_logs = 1.*pow(10.,-6.);
     
     int count2=0;
     index_i = 0;
     index_l = 0;
     index_j = 0; 
     count = 0;

     double complex *x2; 
     class_alloc(x2,(Nmax+1)*sizeof(complex double),pnlpt->error_message);

     for (index_j=0; index_j < Nmax; index_j++){
 
     for (count=0; count < Nmax+1; count++){
         x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
     }
     zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22basic_oneline_complex, x2, &inc, &beta, y, &inc);
     f22_Id2d2[index_j] = 2. * zdotu_(&Nmaxf, x2, &inc, y, &inc);
     
     P_Id2d2[index_j] = fabs(creal(cpow(kdisc[index_j], 3.) * f22_Id2d2[index_j]) - creal(cpow(kdisc[0], 3.) * f22_Id2d2[0]) + epsilon_for_logs);
     }

        
        /*  Having eta2 we compute the rest of the PT matrices */
        
        count2 = 0;
        for (index_l=0; index_l < Nmax+1; index_l++){
            for (index_i=index_l; index_i < Nmax+1; index_i++){
                pnlpt->M_IG2G2[count2] = pnlpt->M22basic_oneline_complex[count2] * ((3.+etam2[index_i]+etam2[index_l])*(1.+etam2[index_i]+etam2[index_l])/((-0.5*etam2[index_i])*(-0.5*etam2[index_l])*(1.-0.5*etam2[index_i])*(1.-0.5*etam2[index_l])));
                
                pnlpt->M_Id2[count2] = pnlpt->M22basic_oneline_complex[count2] * ((3. + etam2[index_i] + etam2[index_l])*(4. + 3.5 *(etam2[index_i] + etam2[index_l]))/(14.*(-0.5*etam2[index_l])*(-0.5 *etam2[index_i])));

                pnlpt->M_Id2G2[count2] = pnlpt->M22basic_oneline_complex[count2] * ((3.+etam2[index_i]+etam2[index_l])/((-0.5*etam2[index_i])*(-0.5*etam2[index_l])));
                
                pnlpt->M_IG2[count2] = pnlpt->M22basic_oneline_complex[count2] * (-1.*(3.+etam2[index_i]+etam2[index_l])*(1.+etam2[index_i]+etam2[index_l])*(6. - 3.5 *(etam2[index_i]+etam2[index_l]))/(28.*(1.-0.5*etam2[index_i])*(1.-0.5*etam2[index_l])*(-0.5*etam2[index_i])*(-0.5*etam2[index_l])));
                
                
                count2++;
            }
        }

        
    //if (pnlpt->rsd == rsd_yes && pnlpt->rsd_only == rsd_only_no || pnlpt->rsd == rsd_no) {
     /* Computing Id2 */
     
     double complex *f22_Id2;
     double *P_Id2;
     class_alloc(f22_Id2,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P_Id2,Nmax*sizeof(double),pnlpt->error_message);

     index_j = 0;
     count = 0;
        /*
        count2 = 0;
        for (index_l=0; index_l < Nmax+1; index_l++){
            for (index_i=index_l; index_i < Nmax+1; index_i++){
                pnlpt->M_Id2[count2] = pnlpt->M22basic_oneline_complex[count2] * ((3. + etam2[index_i] + etam2[index_l])*(4. + 3.5 *(etam2[index_i] + etam2[index_l]))/(14.*(-0.5*etam2[index_l])*(-0.5 *etam2[index_i])));
                count2++;
            }
        }*/

     for (index_j=0; index_j < Nmax; index_j++){

         for (count=0; count < Nmax+1; count++){
             x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
         }
         
         zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M_Id2, x2, &inc, &beta, y, &inc);
         f22_Id2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
         P_Id2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_Id2[index_j]);
     }

     /* Computing IG2 */
     
     double complex *f22_IG2;
     double *P_IG2;
     class_alloc(f22_IG2,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P_IG2,Nmax*sizeof(double),pnlpt->error_message);
   
        
     index_j = 0;
     count = 0;
        /*
        count2 = 0;
        for (index_l=0; index_l < Nmax+1; index_l++){
            for (index_i=index_l; index_i < Nmax+1; index_i++){
                pnlpt->M_IG2[count2] = pnlpt->M22basic_oneline_complex[count2] * (-1.*(3.+etam2[index_i]+etam2[index_l])*(1.+etam2[index_i]+etam2[index_l])*(6. - 3.5 *(etam2[index_i]+etam2[index_l]))/(28.*(1.-0.5*etam2[index_i])*(1.-0.5*etam2[index_l])*(-0.5*etam2[index_i])*(-0.5*etam2[index_l])));
                count2++;
            }
        }*/
        

     for (index_j=0; index_j < Nmax; index_j++){
         for (count=0; count < Nmax+1; count++){
             x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
         };
         zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M_IG2, x2, &inc, &beta, y, &inc);
         f22_IG2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
         //printf("f22_IG2_real=%.18le f22_IG2_imag=%.18le\n",creal(f22_IG2[index_j]),cimag(f22_IG2[index_j]));

         P_IG2[index_j] = fabs(creal(cpow(kdisc[index_j], 3) * f22_IG2[index_j]));
         //printf("%le %le\n",kdisc[j],P_Id2[j]);
     }
        
        
        double *ddpk_PId2;
        double pk_Id2_out = 0.; 
        class_alloc(ddpk_PId2,sizeof(double)*Nmax,pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_Id2,
                                              1,
                                              ddpk_PId2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index=0;
        for (index_k=0; index_k < pnlpt->k_size; index_k++){
            
            if (pnlpt->k[index_k]<=kmax && pnlpt->k[index_k]>=kmin){
                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P_Id2,
                                                    ddpk_PId2,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_Id2_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);
                pk_Id2[index_k] = 10. + pk_Id2_out;
            }
            
            else {
                pk_Id2[index_k] = 10.;
            }
        }
        free(ddpk_PId2);
        
        free(f22_Id2);
        free(P_Id2);
        
        double *ddpk_IG2;
        double pk_IG2_out = 0.;
        class_alloc(ddpk_IG2,sizeof(double)*Nmax,pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_IG2,
                                              1,
                                              ddpk_IG2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index=0;
        for (index_k=0; index_k < pnlpt->k_size; index_k++){
            
            if (pnlpt->k[index_k]<=kmax && pnlpt->k[index_k]>=kmin){
                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P_IG2,
                                                    ddpk_IG2,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_IG2_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);
                pk_IG2[index_k] = pk_IG2_out;
            }
            else {
                pk_IG2[index_k] = epsilon_for_logs;
            }
        }
        free(ddpk_IG2);
        
        free(f22_IG2);
        free(P_IG2);
     
      //  }// end of RSD_only condition
     
     /* Computing Id2G2 */
     
     double complex *f22_Id2G2;
     double *P_Id2G2;
     class_alloc(f22_Id2G2,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P_Id2G2,Nmax*sizeof(double),pnlpt->error_message);
     

     index_j = 0;
        /*
        count2 = 0;
        for (index_l=0; index_l < Nmax+1; index_l++){
            for (index_i=index_l; index_i < Nmax+1; index_i++){
                pnlpt->M_Id2G2[count2] = pnlpt->M22basic_oneline_complex[count2] * ((3.+etam2[index_i]+etam2[index_l])/((-0.5*etam2[index_i])*(-0.5*etam2[index_l])));
                count2++;
            }
        }*/

     for (index_j=0; index_j < Nmax; index_j++){


         for (count=0; count < Nmax+1; count++){
             x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
         }
         zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M_Id2G2, x2, &inc, &beta, y, &inc);
         f22_Id2G2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);

         P_Id2G2[index_j] = fabs(creal(cpow(kdisc[index_j], 3.) * f22_Id2G2[index_j]));
     }
        
     /* Computing IG2G2 */
     
     double complex *f22_IG2G2;
     double *P_IG2G2;
     class_alloc(f22_IG2G2,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P_IG2G2,Nmax*sizeof(double),pnlpt->error_message);


     index_j = 0;
     count = 0;
        /*
        count2 = 0;
        for (index_l=0; index_l < Nmax+1; index_l++){
            for (index_i=index_l; index_i < Nmax+1; index_i++){
                pnlpt->M_IG2G2[count2] = pnlpt->M22basic_oneline_complex[count2] * ((3.+etam2[index_i]+etam2[index_l])*(1.+etam2[index_i]+etam2[index_l])/((-0.5*etam2[index_i])*(-0.5*etam2[index_l])*(1.-0.5*etam2[index_i])*(1.-0.5*etam2[index_l])));
                count2++;
            }
        }*/

     for (index_j=0; index_j < Nmax; index_j++){

         for (count=0; count < Nmax+1; count++){
             x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
         }
         zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M_IG2G2, x2, &inc, &beta, y, &inc);
         f22_IG2G2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);

         P_IG2G2[index_j] = fabs(creal(cpow(kdisc[index_j], 3.) * f22_IG2G2[index_j]));
     }
     
     /* Computing IFG2 */
     
     double complex *f13_IFG2;
     double *P_IFG2;
     class_alloc(f13_IFG2,Nmax*sizeof(complex double),pnlpt->error_message);
     class_alloc(P_IFG2,Nmax*sizeof(double),pnlpt->error_message);

     for (index_j=0; index_j < Nmax; index_j++){
         /*f13_IFG2[index_j] = 0.;
         for (index_i=0; index_i < Nmax+1; index_i++){
             f13_IFG2[index_j] += cmsym2[index_i]* cpow(kdisc[index_j], etam2[index_i]) * (pnlpt->IFG2_oneline[index_i]+ _Complex_I * pnlpt->IFG2_oneline[index_i + Nmax+1]);
         }*/

         for (count=0; count < Nmax+1; count++){
             x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
         };
         
         f13_IFG2[index_j]=zdotu_(&Nmaxf, x2, &inc, pnlpt->IFG2_oneline_complex, &inc);
//         printf("f13_IFG2_real=%le f13_IFG2_imag=%le\n",creal(f13_IFG2[index_j]),cimag(f13_IFG2[index_j]));
         P_IFG2[index_j] = fabs(creal(cpow(kdisc[index_j], 3.) * f13_IFG2[index_j] * Pbin[index_j]));
     }
        
        double *ddpk_IFG2;
        double pk_IFG2_out = 0.;
        class_alloc(ddpk_IFG2,sizeof(double)*Nmax,pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_IFG2,
                                              1,
                                              ddpk_IFG2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index=0;
        for (index_k=0; index_k < pnlpt->k_size; index_k++){
              if (pnlpt->k[index_k]<=kmax && pnlpt->k[index_k]>=kmin){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_IFG2,
                                                ddpk_IFG2,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_IFG2_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            pk_IFG2[index_k] = pk_IFG2_out;
              }
              else {
               pk_IFG2[index_k] = epsilon_for_logs;
              }
        }
        
        
if (pnlpt->rsd == rsd_yes) {
    
    /* Computing Pb1b2 correction in RSD */
    
    double complex *f22_0_b1b2;
    double *P_0_b1b2;
    class_alloc(f22_0_b1b2,Nmax*sizeof(complex double),pnlpt->error_message);
    class_alloc(P_0_b1b2,Nmax*sizeof(double),pnlpt->error_message);
    
    count2 = 0;
    index_i = 0;
    index_l = 0;
    index_j = 0;
    
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu1 = -0.5*etam2[index_i];
            nu2 = -0.5*etam2[index_l];
            
            pnlpt->M22_0_b1b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*(-12.+7.*(3.+f)*nu1+7.*(3.+f)*nu2)/(42.*nu1*nu2);
            
            /*
            
            pnlpt->M22_0_b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (7.*f*f*(12. + 6.*nu1*nu1 - 17.*nu2 + 6.*nu2*nu2 + nu1*(-17. + 12.*nu2))+5.*f*(24. + 14.*nu1*nu1 - 37.*nu2 + 14.*nu2*nu2 + nu1*(-37. + 28.*nu2)))/(210.*nu1*nu2);
            
            pnlpt->M22_0_b1bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*(-1.+2.*nu1+2.*nu2)*(7.*f*(2.+nu1+nu2)+3.*(6.+7.*nu1+7.*nu2))/(42.*nu1*(1. + nu1)*nu2 *(1. + nu2));
            
            pnlpt->M22_0_bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*(-1.+2.*nu1+2.*nu2)*(-10.*f + 7.*f*(5.*(nu1 + nu2) + f*(-2. + 3.*nu1 + 3.*nu2)))/(210.*nu1*(1. + nu1)*nu2 *(1. + nu2));
            
            pnlpt->M22_2_b1b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*f*(nu1+nu2)/(3.*nu1*nu2);
            
            pnlpt->M22_2_b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*f*(-16.+14.*(nu1+nu2)+f*(-13. + 12.*(nu1 + nu2)))/(42.*nu1*nu2);
            
            pnlpt->M22_2_b1bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*(-1.+2.*nu1+2.*nu2)*f*(2.+nu1+nu2)/(3.*nu1*(1. + nu1)*nu2 *(1. + nu2));

            pnlpt->M22_2_bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*(-1.+2.*nu1+2.*nu2)*f*(-2.-f+7.*(nu1 + nu2)+ 6.*f*(nu1 + nu2))/(21.*nu1*(1. + nu1)*nu2 *(1. + nu2));
            
            pnlpt->M22_4_b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*(-1.+2.*nu1+2.*nu2)*2.*f*f/(35.*nu1*nu2);
            
                pnlpt->M22_4_bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*(-1.+2.*nu1+2.*nu2)*4.*f*f*(1.+nu1+nu2)/(35.*nu1*(1. + nu1)*nu2 *(1. + nu2));
            */
            
            count2++;
        }
    }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_0_b1b2_oneline_complex, x2, &inc, &beta, y, &inc);
        f22_0_b1b2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
        P_0_b1b2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_0_b1b2[index_j]);
    }
    
    
    /* Computing b2 correction in RSD */
    
    
    double complex *f22_0_b2;
    double *P_0_b2;
    class_alloc(f22_0_b2,Nmax*sizeof(complex double),pnlpt->error_message);
    class_alloc(P_0_b2,Nmax*sizeof(double),pnlpt->error_message);
    
    count2=0;
    index_i = 0;
    index_l = 0;
    index_j = 0;
    
    
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu1 = -0.5*etam2[index_i];
            nu2 = -0.5*etam2[index_l];
            pnlpt->M22_0_b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (7.*f*f*(12. + 6.*nu1*nu1 - 17.*nu2 + 6.*nu2*nu2 + nu1*(-17. + 12.*nu2))+5.*f*(24. + 14.*nu1*nu1 - 37.*nu2 + 14.*nu2*nu2 + nu1*(-37. + 28.*nu2)))/(210.*nu1*nu2);
            count2++;
        }
    }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_0_b2_oneline_complex, x2, &inc, &beta, y, &inc);
        f22_0_b2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
        P_0_b2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_0_b2[index_j]);
    }
    
    /* Computing b1bG2 correction in RSD */
    
    
    double complex *f22_0_b1bG2;
    double *P_0_b1bG2;
    class_alloc(f22_0_b1bG2,Nmax*sizeof(complex double),pnlpt->error_message);
    class_alloc(P_0_b1bG2,Nmax*sizeof(double),pnlpt->error_message);
    
    count2=0;
    index_i = 0;
    index_l = 0;
    index_j = 0;
    
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu1 = -0.5*etam2[index_i];
            nu2 = -0.5*etam2[index_l];
            pnlpt->M22_0_b1bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*(-1.+2.*nu1+2.*nu2)*(7.*f*(2.+nu1+nu2)+3.*(6.+7.*nu1+7.*nu2))/(42.*nu1*(1. + nu1)*nu2 *(1. + nu2));
            count2++;
        }
    }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_0_b1bG2_oneline_complex, x2, &inc, &beta, y, &inc);
        f22_0_b1bG2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
        P_0_b1bG2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_0_b1bG2[index_j]);
    }
    
    /* Computing bG2 correction in RSD */
    
    
    double complex *f22_0_bG2;
    double *P_0_bG2;
    class_alloc(f22_0_bG2,Nmax*sizeof(complex double),pnlpt->error_message);
    class_alloc(P_0_bG2,Nmax*sizeof(double),pnlpt->error_message);
    
    count2=0;
    index_i = 0;
    index_l = 0;
    index_j = 0;
    
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu1 = -0.5*etam2[index_i];
            nu2 = -0.5*etam2[index_l];
            pnlpt->M22_0_bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*(-1.+2.*nu1+2.*nu2)*(-10.*f + 7.*f*(5.*(nu1 + nu2) + f*(-2. + 3.*nu1 + 3.*nu2)))/(210.*nu1*(1. + nu1)*nu2 *(1. + nu2));
            count2++;
        }
    }
    
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_0_bG2_oneline_complex, x2, &inc, &beta, y, &inc);
        P_0_bG2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
        P_0_bG2[index_j] = creal(cpow(kdisc[index_j], 3) * P_0_bG2[index_j]);
    }
    
    
     /* Computing Pb1b2 correction for the Quadrupole */
    
    double complex *f22_2_b1b2;
    double *P_2_b1b2;
    class_alloc(f22_2_b1b2,Nmax*sizeof(complex double),pnlpt->error_message);
    class_alloc(P_2_b1b2,Nmax*sizeof(double),pnlpt->error_message);
    
    count2=0;
    index_i = 0;
    index_l = 0;
    index_j = 0;
    
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu1 = -0.5*etam2[index_i];
            nu2 = -0.5*etam2[index_l];
            pnlpt->M22_2_b1b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*f*(nu1+nu2)/(3.*nu1*nu2);
            count2++;
        }
    }
    
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_2_b1b2_oneline_complex, x2, &inc, &beta, y, &inc);
        f22_2_b1b2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
        P_2_b1b2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_2_b1b2[index_j]);
    }
    
    /* Computing Pb2 correction for the Quadrupole */
    
    double complex *f22_2_b2;
    double *P_2_b2;
    class_alloc(f22_2_b2,Nmax*sizeof(complex double),pnlpt->error_message);
    class_alloc(P_2_b2,Nmax*sizeof(double),pnlpt->error_message);
    
    count2=0;
    index_i = 0;
    index_l = 0;
    index_j = 0;
    
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu1 = -0.5*etam2[index_i];
            nu2 = -0.5*etam2[index_l];
            pnlpt->M22_2_b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*f*(-16.+14.*(nu1+nu2)+f*(-13. + 12.*(nu1 + nu2)))/(42.*nu1*nu2);
            count2++;
        }
    }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_2_b2_oneline_complex, x2, &inc, &beta, y, &inc);
        f22_2_b2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
        P_2_b2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_2_b2[index_j]);
    }
    
    
    /* Computing b1bG2 correction for Quadrupole */
    
    
    double complex *f22_2_b1bG2;
    double *P_2_b1bG2;
    class_alloc(f22_2_b1bG2,Nmax*sizeof(complex double),pnlpt->error_message);
    class_alloc(P_2_b1bG2,Nmax*sizeof(double),pnlpt->error_message);
    
    count2=0;
    index_i = 0;
    index_l = 0;
    index_j = 0;

    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu1 = -0.5*etam2[index_i];
            nu2 = -0.5*etam2[index_l];
            pnlpt->M22_2_b1bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*(-1.+2.*nu1+2.*nu2)*f*(2.+nu1+nu2)/(3.*nu1*(1. + nu1)*nu2 *(1. + nu2));
            count2++;
        }
    }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_2_b1bG2_oneline_complex, x2, &inc, &beta, y, &inc);
        f22_2_b1bG2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
        P_2_b1bG2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_2_b1bG2[index_j]);
    }

    /* Computing bG2 correction for Quadrupole */
    
    
    double complex *f22_2_bG2;
    double *P_2_bG2;
    class_alloc(f22_2_bG2,Nmax*sizeof(complex double),pnlpt->error_message);
    class_alloc(P_2_bG2,Nmax*sizeof(double),pnlpt->error_message);
    
    
    count2 = 0;
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu1 = -0.5*etam2[index_i];
            nu2 = -0.5*etam2[index_l];
            pnlpt->M22_2_bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*(-1.+2.*nu1+2.*nu2)*f*(-2.-f+7.*(nu1 + nu2)+ 6.*f*(nu1 + nu2))/(21.*nu1*(1. + nu1)*nu2 *(1. + nu2));
            count2++;
        }
    }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_2_bG2_oneline_complex, x2, &inc, &beta, y, &inc);
        P_2_bG2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
        P_2_bG2[index_j] = creal(cpow(kdisc[index_j], 3) * P_2_bG2[index_j]);
    }
    
    
    /* Computing Pb2 correction for Hexadecapole */
    
    double complex *f22_4_b2;
    double *P_4_b2;
    class_alloc(f22_4_b2,Nmax*sizeof(complex double),pnlpt->error_message);
    class_alloc(P_4_b2,Nmax*sizeof(double),pnlpt->error_message);
    
    
    count2=0;
    index_i = 0;
    index_l = 0;
    index_j = 0;
    
    
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu1 = -0.5*etam2[index_i];
            nu2 = -0.5*etam2[index_l];
            pnlpt->M22_4_b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*(-1.+2.*nu1+2.*nu2)*2.*f*f/(35.*nu1*nu2);
            count2++;
        }
    }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_4_b2_oneline_complex, x2, &inc, &beta, y, &inc);
        f22_4_b2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
        P_4_b2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_4_b2[index_j]);
    }
    
    /* Computing bG2 correction for Hexadecapole */
    
    double complex *f22_4_bG2;
    double *P_4_bG2;
    class_alloc(f22_4_bG2,Nmax*sizeof(complex double),pnlpt->error_message);
    class_alloc(P_4_bG2,Nmax*sizeof(double),pnlpt->error_message);
    
    
    count2 = 0;
    index_i = 0;
    index_l = 0;
    index_j = 0;
    
    for (index_l=0; index_l < Nmax+1; index_l++){
        for (index_i=index_l; index_i < Nmax+1; index_i++){
            nu1 = -0.5*etam2[index_i];
            nu2 = -0.5*etam2[index_l];
            pnlpt->M22_4_bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3.+2.*nu1+2.*nu2)*(-1.+2.*nu1+2.*nu2)*4.*f*f*(1.+nu1+nu2)/(35.*nu1*(1. + nu1)*nu2 *(1. + nu2));
            count2++;
        }
    }
    
    for (index_j=0; index_j < Nmax; index_j++){
        for (count=0; count < Nmax+1; count++){
            x2[count]= cmsym2[count]* cpow(kdisc[index_j], etam2[count]);
        }
        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_4_bG2_oneline_complex, x2, &inc, &beta, y, &inc);
        P_4_bG2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
        P_4_bG2[index_j] = creal(cpow(kdisc[index_j], 3) * P_4_bG2[index_j]);
    }
    
    /* Computing IFG2 for monopole and quadrupole */
    
    /*
    if (irindex == 0) {
       // printf("Computing IFG2 without IR resummation! \n");
        for (index_j=0; index_j < Nmax; index_j++){
        P_IFG2_0b1_x[index_j] = P_IFG2[index_j];
        P_IFG2_0[index_j] = P_IFG2[index_j]*f/3.;
        P_IFG2_2[index_j] = P_IFG2[index_j]*2.*f/3.;
        }
    }
    
    if (irindex == 1){
       // printf("Computing IFG2 with IR resummation! \n");
    for (index_j=0; index_j < Nmax; index_j++){
        P_IFG2_0b1_x[index_j] = P_IFG2[index_j]*(P10b1[index_j]/Pbin[index_j]);
        P_IFG2_0[index_j] = P_IFG2[index_j]*(P10[index_j]/Pbin[index_j]);
        P_IFG2_2[index_j] = P_IFG2[index_j]*(P12[index_j]/Pbin[index_j]);
     //   printf("%le %le %le %le\n",kdisc[index_j],P_IFG2_0b1[index_j],P_IFG2_0[index_j],P_IFG2_2[index_j]);
    }
    }
     */
    
    
   // Numerical integration for bias tracers
    
    double *P_Id2d2_2;
    class_alloc(P_Id2d2_2,Nmax*sizeof(double),pnlpt->error_message);
    double *P_Id2d2_4;
    class_alloc(P_Id2d2_4,Nmax*sizeof(double),pnlpt->error_message);
    double *P_Id2G2_2;
    class_alloc(P_Id2G2_2,Nmax*sizeof(double),pnlpt->error_message);
    double *P_Id2G2_4;
    class_alloc(P_Id2G2_4,Nmax*sizeof(double),pnlpt->error_message);
    double *P_IG2G2_2;
    class_alloc(P_IG2G2_2,Nmax*sizeof(double),pnlpt->error_message);
    double *P_IG2G2_4;
    class_alloc(P_IG2G2_4,Nmax*sizeof(double),pnlpt->error_message);
    double *P_4_b1b2;
    class_alloc(P_4_b1b2,Nmax*sizeof(double),pnlpt->error_message);
    double *P_4_b1bG2;
    class_alloc(P_4_b1bG2,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P_IFG2_0b1_x;
    class_alloc(P_IFG2_0b1_x,Nmax*sizeof(double),pnlpt->error_message);
    double *P_IFG2_0;
    class_alloc(P_IFG2_0,Nmax*sizeof(double),pnlpt->error_message);
    double *P_IFG2_2;
    class_alloc(P_IFG2_2,Nmax*sizeof(double),pnlpt->error_message);
    
    double *P_Id2d2_new;
    class_alloc(P_Id2d2_new,sizeof(double)*Nmax,pnlpt->error_message);
    double *P_Id2G2_new;
    class_alloc(P_Id2G2_new,sizeof(double)*Nmax,pnlpt->error_message);
    double *P_IG2G2_new;
    class_alloc(P_IG2G2_new,sizeof(double)*Nmax,pnlpt->error_message);
    
    double *P_0_bG2_new;
    class_alloc(P_0_bG2_new,sizeof(double)*Nmax,pnlpt->error_message);
    double *P_2_bG2_new;
    class_alloc(P_2_bG2_new,sizeof(double)*Nmax,pnlpt->error_message);
    double *P_4_bG2_new;
    class_alloc(P_4_bG2_new,sizeof(double)*Nmax,pnlpt->error_message);
    
    double *P_IFG2_new;
    class_alloc(P_IFG2_new,sizeof(double)*Nmax,pnlpt->error_message);
    
    
    double *P_0_b1b2_new;
    class_alloc(P_0_b1b2_new,sizeof(double)*Nmax,pnlpt->error_message);
    double *P_2_b1b2_new;
    class_alloc(P_2_b1b2_new,sizeof(double)*Nmax,pnlpt->error_message);
    double *P_0_b1bG2_new;
    class_alloc(P_0_b1bG2_new,sizeof(double)*Nmax,pnlpt->error_message);
    double *P_2_b1bG2_new;
    class_alloc(P_2_b1bG2_new,sizeof(double)*Nmax,pnlpt->error_message);
    double *P_0_b2_new;
    class_alloc(P_0_b2_new,sizeof(double)*Nmax,pnlpt->error_message);
    double *P_2_b2_new;
    class_alloc(P_2_b2_new,sizeof(double)*Nmax,pnlpt->error_message);
    double *P_4_b2_new;
    class_alloc(P_4_b2_new,sizeof(double)*Nmax,pnlpt->error_message);
    
    
    for (index_j=0; index_j < Nmax; index_j++){
        P_Id2d2_new[index_j] = P_Id2d2[index_j];
        P_Id2G2_new[index_j] = P_Id2G2[index_j];
        P_IG2G2_new[index_j] = P_IG2G2[index_j];
        P_0_b1b2_new[index_j] = P_0_b1b2[index_j];
        P_2_b1b2_new[index_j] = P_2_b1b2[index_j];
        P_0_b1bG2_new[index_j] = P_0_b1bG2[index_j];
        P_2_b1bG2_new[index_j] = P_2_b1bG2[index_j];
        P_4_b2_new[index_j] = P_4_b2[index_j];
        P_2_b2_new[index_j] = P_2_b2[index_j];
        P_0_b2_new[index_j] = P_0_b2[index_j];
        P_4_bG2_new[index_j] = P_4_bG2[index_j];
        P_2_bG2_new[index_j] = P_2_bG2[index_j];
        P_0_bG2_new[index_j] = P_0_bG2[index_j];
        P_IFG2_new[index_j] = P_IFG2[index_j];
    }
    
    double *dd_P_Id2d2;
    class_alloc(dd_P_Id2d2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_Id2d2, 1, dd_P_Id2d2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_Id2d2_ap_out=0;
    
    double *dd_P_Id2G2;
    class_alloc(dd_P_Id2G2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_Id2G2, 1, dd_P_Id2G2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_Id2G2_ap_out=0;
    
    double *dd_P_IG2G2;
    class_alloc(dd_P_IG2G2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_IG2G2, 1, dd_P_IG2G2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_IG2G2_ap_out=0;
    
    double *dd_P_0_b1b2;
    class_alloc(dd_P_0_b1b2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_0_b1b2_new, 1, dd_P_0_b1b2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_0_b1b2_ap_out=0;
    
    double *dd_P_2_b1b2;
    class_alloc(dd_P_2_b1b2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_2_b1b2_new, 1, dd_P_2_b1b2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_2_b1b2_ap_out=0;
    
    double *dd_P_0_b1bG2;
    class_alloc(dd_P_0_b1bG2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_0_b1bG2, 1, dd_P_0_b1bG2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_0_b1bG2_ap_out=0;
    
    double *dd_P_2_b1bG2;
    class_alloc(dd_P_2_b1bG2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_2_b1bG2, 1, dd_P_2_b1bG2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_2_b1bG2_ap_out=0;
    
    double *dd_P_0_b2;
    class_alloc(dd_P_0_b2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_0_b2_new, 1, dd_P_0_b2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_0_b2_ap_out=0;
    double *dd_P_2_b2;
    class_alloc(dd_P_2_b2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_2_b2_new, 1, dd_P_2_b2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_2_b2_ap_out=0;
    double *dd_P_4_b2;
    class_alloc(dd_P_4_b2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_4_b2_new, 1, dd_P_4_b2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_4_b2_ap_out=0;
    
    double *dd_P_0_bG2;
    class_alloc(dd_P_0_bG2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_0_bG2_new, 1, dd_P_0_bG2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_0_bG2_ap_out=0;
    double *dd_P_2_bG2;
    class_alloc(dd_P_2_bG2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_2_bG2_new, 1, dd_P_2_bG2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_2_bG2_ap_out=0;
    double *dd_P_4_bG2;
    class_alloc(dd_P_4_bG2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_4_bG2_new, 1, dd_P_4_bG2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_4_bG2_ap_out=0;
    
    double *dd_P_IFG2;
    class_alloc(dd_P_IFG2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_IFG2_new, 1, dd_P_IFG2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double P_IFG2_ap_out=0;
    
    double *dd_Pbin;
    class_alloc(dd_Pbin,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,Pbin, 1, dd_Pbin,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    double Pbin_ap_out=0;

    double Pd2d2_in=0.;
    double Pd2G2_in=0.;
    double PG2G2_in=0.;
    double Pb1b2_in=0.;
    double Pb1bG2_in=0.;
    double Pb2_in=0.;
    double PbG2_in=0.;
    double IFG2_in=0.;

    double LegendreP0 = 1.;
    double LegendreP2=0.;
    double LegendreP4=0.;
    
    double LegendreP2true=0.;
    double LegendreP4true=0.;
    
    int index_gauss2 = 0;
    double mu=0.;
    double mutrue = 0.;
    double ktrue = 0.;
    last_index=0;
    double Sigmatot = 0.;
    double p_lo = 0.;
    double Exp = 0.;
    
    /*
    hratio = 1.;
    Dratio = 1.;
    */
    
    for (index_j=0; index_j < Nside; index_j++){
    P_IFG2_0b1_x[index_j] = 0.;
    P_IFG2_0[index_j] = 0.;
    P_IFG2_2[index_j] = 0.;
    
    P_IFG2_0b1_x[Nmax-1-index_j] = 0.;
    P_IFG2_0[Nmax-1-index_j] = 0.;
    P_IFG2_2[Nmax-1-index_j] = 0.;
    
    P_Id2d2_2[index_j] = 0.;
    P_Id2d2_4[index_j] = 0.;
    P_Id2d2_2[Nmax-1-index_j] = 0.;
    P_Id2d2_4[Nmax-1-index_j] = 0.;
    
    P_Id2G2_2[index_j] = 0.;
    P_Id2G2_4[index_j] = 0.;
    P_Id2G2_2[Nmax-1-index_j] = 0.;
    P_Id2G2_4[Nmax-1-index_j] = 0.;
    
    P_IG2G2_2[index_j] = 0.;
    P_IG2G2_4[index_j] = 0.;
    P_IG2G2_2[Nmax-1-index_j] = 0.;
    P_IG2G2_4[Nmax-1-index_j] = 0.;
    
    P_4_b1b2[index_j] = 0.;
    P_4_b1bG2[index_j] = 0.;
    P_4_b1b2[Nmax-1-index_j] = 0.;
    P_4_b1bG2[Nmax-1-index_j] = 0.;
    }
    
    for (index_j=Nside; index_j < Nmax-Nside; index_j++){
        
        P_IFG2_0b1_x[index_j] = 0.;
        P_IFG2_0[index_j] = 0.;
        P_IFG2_2[index_j] = 0.;
        
        P_Id2d2_2[index_j] = 0.;
        P_Id2d2_4[index_j] = 0.;
        P_Id2d2[index_j] = 0.;
        
        P_Id2G2[index_j] = 0.;
        P_Id2G2_2[index_j] = 0.;
        P_Id2G2_4[index_j] = 0.;
        
        P_IG2G2[index_j] = 0.;
        P_IG2G2_2[index_j] = 0.;
        P_IG2G2_4[index_j] = 0.;
        
        P_0_b1b2[index_j] = 0.;
        P_2_b1b2[index_j] = 0.;
        
        P_0_b1bG2[index_j] = 0.;
        P_2_b1bG2[index_j] = 0.;
        
        P_0_b2[index_j] = 0.;
        P_2_b2[index_j] = 0.;
        P_4_b2[index_j] = 0.;
        
        P_0_bG2[index_j] = 0.;
        P_2_bG2[index_j] = 0.;
        P_4_bG2[index_j] = 0.;
        
        P_4_b1b2[index_j] = 0.;
        P_4_b1bG2[index_j] = 0.;
        
        Pnw_ap_out = 0.;
        Pw_ap_out = 0.;
        
        for (index_gauss2=0; index_gauss2 < 40; index_gauss2++){

            /*
            mu = pnlpt->gauss_x[index_gauss2];
            mutrue = mu*hratio/pow((1./Dratio/Dratio + (hratio*hratio - 1./Dratio/Dratio)*mu*mu),0.5);
            ktrue = kdisc[index_j]*pow((1./Dratio/Dratio + (hratio*hratio - 1./Dratio/Dratio)*mu*mu),0.5);
            
            
            mutrue = mu;
            ktrue = kdisc[index_j];
             */
            
            mu = pnlpt->gauss_x[index_gauss2];
            
            if (pnlpt->AP_effect == AP_effect_yes){
                mutrue = mu*hratio/pow((1./Dratio/Dratio + (hratio*hratio - 1./Dratio/Dratio)*mu*mu),0.5);
                ktrue = kdisc[index_j]*pow((1./Dratio/Dratio + (hratio*hratio - 1./Dratio/Dratio)*mu*mu),0.5);
            }
            
            else {
                mutrue = mu;
                ktrue =kdisc[index_j];
            }
                        
            LegendreP2 = (3.*pow(mu,2.)-1.)/2.;
            LegendreP4 = (35.*pow(mu,4.)-30.*pow(mu,2.)+3.)/8.;
            LegendreP2true = (3.*pow(mutrue,2.)-1.)/2.;
            LegendreP4true = (35.*pow(mutrue,4.)-30.*pow(mutrue,2.)+3.)/8.;
            
            class_call(array_interpolate_spline(kdisc,Nmax,Pnw,dd_Pnw,1,ktrue,&last_index,&Pnw_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,Pw,dd_Pw,1,ktrue,&last_index,&Pw_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,Pbin,dd_Pbin,1,ktrue,&last_index,&Pbin_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
            class_call(array_interpolate_spline(kdisc,Nmax,P_IFG2_new,dd_P_IFG2,1,ktrue,&last_index,&P_IFG2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
            class_call(array_interpolate_spline(kdisc,Nmax,P_Id2d2_new,dd_P_Id2d2,1,ktrue,&last_index,&P_Id2d2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P_Id2G2_new,dd_P_Id2G2,1,ktrue,&last_index,&P_Id2G2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P_IG2G2_new,dd_P_IG2G2,1,ktrue,&last_index,&P_IG2G2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
            class_call(array_interpolate_spline(kdisc,Nmax,P_0_b1b2_new,dd_P_0_b1b2,1,ktrue,&last_index,&P_0_b1b2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P_2_b1b2_new,dd_P_2_b1b2,1,ktrue,&last_index,&P_2_b1b2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
            class_call(array_interpolate_spline(kdisc,Nmax,P_0_b1bG2_new,dd_P_0_b1bG2,1,ktrue,&last_index,&P_0_b1bG2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P_2_b1bG2_new,dd_P_2_b1bG2,1,ktrue,&last_index,&P_2_b1bG2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
             class_call(array_interpolate_spline(kdisc,Nmax,P_0_b2_new,dd_P_0_b2,1,ktrue,&last_index,&P_0_b2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
             class_call(array_interpolate_spline(kdisc,Nmax,P_2_b2_new,dd_P_2_b2,1,ktrue,&last_index,&P_2_b2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
             class_call(array_interpolate_spline(kdisc,Nmax,P_4_b2_new,dd_P_4_b2,1,ktrue,&last_index,&P_4_b2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
            class_call(array_interpolate_spline(kdisc,Nmax,P_0_bG2_new,dd_P_0_bG2,1,ktrue,&last_index,&P_0_bG2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P_2_bG2_new,dd_P_2_bG2,1,ktrue,&last_index,&P_2_bG2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            class_call(array_interpolate_spline(kdisc,Nmax,P_4_bG2_new,dd_P_4_bG2,1,ktrue,&last_index,&P_4_bG2_ap_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            
  //          Sigmatot = SigmaBAO*(1.+f*mutrue*mutrue*(2.+f))+f*f*mutrue*mutrue*(mutrue*mutrue-1.)*deltaSigmaBAO;
  //          Exp = exp(-Sigmatot * pow(ktrue,2.));
 //           P13ratio = 1.+(Pw_ap_out/Pnw_ap_out)*Exp;
 //           P1b1 = (Pnw[index_j]+Pw[index_j]*Exp)*pnlpt->gauss_w[index_gauss2];
 //           P1 = (Pnw[index_j]+Pw[index_j]*Exp)*f*pow(pnlpt->gauss_x[index_gauss2],2.)*pnlpt->gauss_w[index_gauss2];
            
            Sigmatot = SigmaBAO*(1.+f*mutrue*mutrue*(2.+f))+f*f*mutrue*mutrue*(mutrue*mutrue-1.)*deltaSigmaBAO;
            Exp = exp(-Sigmatot * pow(ktrue,2.));
            p_lo = (Pnw_ap_out + Pw_ap_out * Exp)/Pbin_ap_out;
            
            IFG2_in = p_lo*P_IFG2_ap_out*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            
            Pd2d2_in = P_Id2d2_ap_out*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio ;
            Pd2G2_in = P_Id2G2_ap_out*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio ;
            PG2G2_in = P_IG2G2_ap_out*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio ;
            Pb1b2_in = (P_0_b1b2_ap_out + LegendreP2true*P_2_b1b2_ap_out)*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            Pb1bG2_in = (P_0_b1bG2_ap_out + LegendreP2true*P_2_b1bG2_ap_out)*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            Pb2_in = (P_0_b2_ap_out + LegendreP2true*P_2_b2_ap_out + LegendreP4true*P_4_b2_ap_out)*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            PbG2_in = (P_0_bG2_ap_out + LegendreP2true*P_2_bG2_ap_out + LegendreP4true*P_4_bG2_ap_out)*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            /*
            P1loopdd = ((p_tree + P22_mu0_dd_ap_out + P13_mu0_dd_ap_out*P13ratio + (P13_mu0_dd_w_ap_out+P22_mu0_dd_w_ap_out)*Exp)+(P22_mu2_dd_ap_out + P13_mu2_dd_ap_out*P13ratio + (P22_mu2_dd_w_ap_out+P13_mu2_dd_w_ap_out)*Exp)*pow(mutrue,2.)+(P22_mu4_dd_ap_out+P22_mu4_dd_w_ap_out*Exp)*pow(mutrue,4.))*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            
            P1loopvd = ((p_tree*2.*f + P13_mu2_vd_ap_out*P13ratio + P22_mu2_vd_ap_out+ (P22_mu2_vd_w_ap_out+P13_mu2_vd_w_ap_out)*Exp)*pow(mutrue,2.)+(P13_mu4_vd_ap_out*P13ratio + P22_mu4_vd_ap_out + (P22_mu4_vd_w_ap_out+P13_mu4_vd_w_ap_out)*Exp)*pow(mutrue,4.)+(P22_mu6_vd_ap_out+P22_mu6_vd_w_ap_out*Exp)*pow(mutrue,6.))*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            */
            P_IFG2_0b1_x[index_j] += IFG2_in*LegendreP0/2.;
            P_IFG2_0[index_j] += IFG2_in*f*pow(mutrue,2.)*LegendreP0/2.;
            P_IFG2_2[index_j] += IFG2_in*f*pow(mutrue,2.)*LegendreP2*2.5;
            P_Id2d2[index_j] += Pd2d2_in*LegendreP0/2.;
            P_Id2d2_2[index_j] += Pd2d2_in*LegendreP2*2.5;
            P_Id2d2_4[index_j] += Pd2d2_in*LegendreP4*4.5;
            P_Id2G2[index_j] += Pd2G2_in*LegendreP0/2.;
            P_Id2G2_2[index_j] += Pd2G2_in*LegendreP2*2.5;
            P_Id2G2_4[index_j] += Pd2G2_in*LegendreP4*4.5;
            P_IG2G2[index_j] += PG2G2_in*LegendreP0/2.;
            P_IG2G2_2[index_j] += PG2G2_in*LegendreP2*2.5;
            P_IG2G2_4[index_j] += PG2G2_in*LegendreP4*4.5;
            P_0_b1b2[index_j] += Pb1b2_in*LegendreP0/2.;
            P_2_b1b2[index_j] += Pb1b2_in*LegendreP2*2.5;
            P_0_b1bG2[index_j] += Pb1bG2_in*LegendreP0/2.;
            P_2_b1bG2[index_j] += Pb1bG2_in*LegendreP2*2.5;
            P_0_b2[index_j] +=  Pb2_in*LegendreP0/2.;
            P_2_b2[index_j] +=  Pb2_in*LegendreP2*2.5;
            P_4_b2[index_j] +=  Pb2_in*LegendreP4*4.5;
            P_0_bG2[index_j] +=  PbG2_in*LegendreP0/2.;
            P_2_bG2[index_j] +=  PbG2_in*LegendreP2*2.5;
            P_4_bG2[index_j] +=  PbG2_in*LegendreP4*4.5;
            P_4_b1b2[index_j] += Pb1b2_in*LegendreP4*4.5;
            P_4_b1bG2[index_j] += Pb1bG2_in*LegendreP4*4.5;;
            
     }
   //     printf("%lf %lf %lf\n",P_IFG2_0b1_x[index_j],P_IFG2_0[index_j],P_IFG2_2[index_j]);
     }
    
    free(P_Id2d2_new);
    free(P_Id2G2_new);
    free(P_IG2G2_new);
    
    free(P_0_b1b2_new);
    free(P_2_b1b2_new);
    free(P_0_b1bG2_new);
    free(P_2_b1bG2_new);
    free(P_0_b2_new);
    free(P_2_b2_new);
    free(P_4_b2_new);
    free(P_0_bG2_new);
    free(P_2_bG2_new);
    free(P_4_bG2_new);
    
    free(dd_P_Id2d2);
    free(dd_P_Id2G2);
    free(dd_P_IG2G2);
    free(dd_P_0_b1b2);
    free(dd_P_2_b1b2);
    free(dd_P_0_b1bG2);
    free(dd_P_2_b1bG2);
    free(dd_P_0_b2);
    free(dd_P_2_b2);
    free(dd_P_4_b2);
    free(dd_P_0_bG2);
    free(dd_P_2_bG2);
    free(dd_P_4_bG2);
    
    free(dd_P_IFG2);
    free(dd_Pbin);
    
    double pk_Id2d2_2_out=0.;
    double *ddpk_PId2d2_2;
    class_alloc(ddpk_PId2d2_2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_Id2d2_2,1,ddpk_PId2d2_2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,P_Id2d2_2,ddpk_PId2d2_2,1,pnlpt->k[index_k],&last_index,&pk_Id2d2_2_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_Id2d2_2[index_k] = pk_Id2d2_2_out + 1.e7;
           // printf("%lf\n",pk_Id2d2_2[index_k]);
        }
        else {
            pk_Id2d2_2[index_k] = 1.e7;
        }
    }
    free(ddpk_PId2d2_2);

    double pk_Id2d2_4_out=0;
    double *ddpk_PId2d2_4;
    class_alloc(ddpk_PId2d2_4,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_Id2d2_4,1,ddpk_PId2d2_4,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,P_Id2d2_4,ddpk_PId2d2_4,1,pnlpt->k[index_k],&last_index,&pk_Id2d2_4_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_Id2d2_4[index_k] = pk_Id2d2_4_out + 1.e7;
        }
        else {
            pk_Id2d2_4[index_k] = 1.e7;
        }
    }
    free(ddpk_PId2d2_4);
    
    double pk_Id2G2_2_out=0.;
    double *ddpk_PId2G2_2;
    class_alloc(ddpk_PId2G2_2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_Id2G2_2,1,ddpk_PId2G2_2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,P_Id2G2_2,ddpk_PId2G2_2,1,pnlpt->k[index_k],&last_index,&pk_Id2G2_2_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_Id2G2_2[index_k] = pk_Id2G2_2_out + 1.e7;
        }
        else {
            pk_Id2G2_2[index_k] = 0.*epsilon_for_logs + 1.e7;
        }
    }
    free(ddpk_PId2G2_2);
    
    
    double pk_Id2G2_4_out=0.;
    double *ddpk_PId2G2_4;
    class_alloc(ddpk_PId2G2_4,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_Id2G2_4,1,ddpk_PId2G2_4,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,P_Id2G2_4,ddpk_PId2G2_4,1,pnlpt->k[index_k],&last_index,&pk_Id2G2_4_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_Id2G2_4[index_k] = pk_Id2G2_4_out + 1.e7;
        }
        else {
            pk_Id2G2_4[index_k] = 0.*epsilon_for_logs + 1.e7;
        }
    }
    free(ddpk_PId2G2_4);
    
    
    double pk_IG2G2_2_out=0.;
    double *ddpk_PIG2G2_2;
    class_alloc(ddpk_PIG2G2_2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_IG2G2_2,1,ddpk_PIG2G2_2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,P_IG2G2_2,ddpk_PIG2G2_2,1,pnlpt->k[index_k],&last_index,&pk_IG2G2_2_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_IG2G2_2[index_k] = pk_IG2G2_2_out + 1.e7;
        }
        else {
            pk_IG2G2_2[index_k] = 0.*epsilon_for_logs + 1.e7;
        }
//        printf("%le %le %le\n", pk_Id2d2_2[index_k], pk_Id2G2_2[index_k],pk_IG2G2_2[index_k]);
    }
    free(ddpk_PIG2G2_2);
    
    double pk_IG2G2_4_out=0.;
    double *ddpk_PIG2G2_4;
    class_alloc(ddpk_PIG2G2_4,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_IG2G2_4,1,ddpk_PIG2G2_4,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,P_IG2G2_4,ddpk_PIG2G2_4,1,pnlpt->k[index_k],&last_index,&pk_IG2G2_4_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_IG2G2_4[index_k] = pk_IG2G2_4_out + 1.e7;
        }
        else {
            pk_IG2G2_4[index_k] = 0.*epsilon_for_logs + 1.e7;
        }
 //       printf("%lf\n",pk_IG2G2_4[index_k]);
 //       printf("%le %le %le\n", pk_Id2d2_4[index_k], pk_Id2G2_4[index_k],pk_IG2G2_4[index_k]);

    }
    free(ddpk_PIG2G2_4);
    
    
    
    double pk_0_b1b2_out=0;
    double *ddpk_0_b1b2;
    class_alloc(ddpk_0_b1b2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_0_b1b2,1,ddpk_0_b1b2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,P_0_b1b2,ddpk_0_b1b2,1,pnlpt->k[index_k],&last_index,&pk_0_b1b2_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_l_0_b1b2[index_k] = pk_0_b1b2_out + 1.e7;
        }
        else {
            pk_l_0_b1b2[index_k] = epsilon_for_logs + 1.e7;
        }
  //      printf("%lf\n", pk_l_0_b1b2[index_k]);
    }
    free(ddpk_0_b1b2);
    free(f22_0_b1b2);
    free(P_0_b1b2);
    
    double pk_0_b2_out=0;
    double *ddpk_0_b2;
    class_alloc(ddpk_0_b2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P_0_b2,
                                          1,
                                          ddpk_0_b2,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_0_b2,
                                                ddpk_0_b2,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_0_b2_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            pk_l_0_b2[index_k] = pk_0_b2_out + 1.e7;
        }
        else {
            pk_l_0_b2[index_k] = epsilon_for_logs + 1.e7;
        }
    }
    free(ddpk_0_b2);
    free(f22_0_b2);
    free(P_0_b2);
    
    double pk_0_b1bG2_out;
    double *ddpk_0_b1bG2;
    class_alloc(ddpk_0_b1bG2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P_0_b1bG2,
                                          1,
                                          ddpk_0_b1bG2,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_0_b1bG2,
                                                ddpk_0_b1bG2,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_0_b1bG2_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            pk_l_0_b1bG2[index_k] = pk_0_b1bG2_out + 1.e7;
        }
        else {
            pk_l_0_b1bG2[index_k] = epsilon_for_logs + 1.e7;
        }
    }
    free(ddpk_0_b1bG2);
    free(f22_0_b1bG2);
    free(P_0_b1bG2);
    
    double pk_2_b1b2_out;
    double *ddpk_2_b1b2;
    class_alloc(ddpk_2_b1b2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P_2_b1b2,
                                          1,
                                          ddpk_2_b1b2,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_2_b1b2,
                                                ddpk_2_b1b2,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_2_b1b2_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            pk_l_2_b1b2[index_k] = pk_2_b1b2_out + 1.e7;
        }
        else {
            pk_l_2_b1b2[index_k] = epsilon_for_logs + 1.e7;
        }
    }
    free(ddpk_2_b1b2);
    free(f22_2_b1b2);
    free(P_2_b1b2);
    
    
    double pk_4_b1b2_out;
    double *ddpk_4_b1b2;
    class_alloc(ddpk_4_b1b2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_4_b1b2,1,ddpk_4_b1b2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,P_4_b1b2,ddpk_4_b1b2,1,pnlpt->k[index_k],&last_index,&pk_4_b1b2_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_l_4_b1b2[index_k] = pk_4_b1b2_out + 1.e7;
        }
        else {
            pk_l_4_b1b2[index_k] = epsilon_for_logs + 1.e7;
        }
    }
    free(ddpk_4_b1b2);
    
    
    double pk_0_bG2_out;
    double *ddpk_0_bG2;
    class_alloc(ddpk_0_bG2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P_0_bG2,
                                          1,
                                          ddpk_0_bG2,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_0_bG2,
                                                ddpk_0_bG2,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_0_bG2_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            pk_l_0_bG2[index_k] = pk_0_bG2_out + 1.e7;
        }
        else {
            pk_l_0_bG2[index_k] = epsilon_for_logs + 1.e7;
        }
    }
    free(ddpk_0_bG2);
    free(f22_0_bG2);
    free(P_0_bG2);
    
    double pk_2_b2_out;
    double *ddpk_2_b2;
    class_alloc(ddpk_2_b2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P_2_b2,
                                          1,
                                          ddpk_2_b2,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_2_b2,
                                                ddpk_2_b2,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_2_b2_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            pk_l_2_b2[index_k] = pk_2_b2_out + 1.e7;
        }
        else {
            pk_l_2_b2[index_k] = epsilon_for_logs + 1.e7;
        }
    }
    free(ddpk_2_b2);
    free(f22_2_b2);
    free(P_2_b2);

    
    double *ddpk_2_b1bG2;
    double pk_2_b1bG2_out;
    class_alloc(ddpk_2_b1bG2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P_2_b1bG2,
                                          1,
                                          ddpk_2_b1bG2,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_2_b1bG2,
                                                ddpk_2_b1bG2,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_2_b1bG2_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            pk_l_2_b1bG2[index_k] = pk_2_b1bG2_out + 1.e7;
        }
        else {
            pk_l_2_b1bG2[index_k] = epsilon_for_logs + 1.e7;
        }
    }
    free(ddpk_2_b1bG2);
    free(f22_2_b1bG2);
    free(P_2_b1bG2);
    
    double pk_4_b1bG2_out;
    double *ddpk_4_b1bG2;
    class_alloc(ddpk_4_b1bG2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_4_b1bG2,1,ddpk_4_b1bG2,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,P_4_b1bG2,ddpk_4_b1bG2,1,pnlpt->k[index_k],&last_index,&pk_4_b1bG2_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_l_4_b1bG2[index_k] = pk_4_b1bG2_out + 1.e7;
        }
        else {
            pk_l_4_b1bG2[index_k] = epsilon_for_logs + 1.e7;
        }
  //      printf("%le %le\n", pk_l_4_b1b2[index_k],pk_l_4_b1bG2[index_k]);
    }
    free(ddpk_4_b1bG2);


    double *ddpk_2_bG2;
    double pk_2_bG2_out;
    class_alloc(ddpk_2_bG2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P_2_bG2,
                                          1,
                                          ddpk_2_bG2,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_2_bG2,
                                                ddpk_2_bG2,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_2_bG2_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            pk_l_2_bG2[index_k] = pk_2_bG2_out + 1.e7;
        }
        else {
            pk_l_2_bG2[index_k] = epsilon_for_logs + 1.e7;
        }
    }
    free(ddpk_2_bG2);
    free(f22_2_bG2);
    free(P_2_bG2);

    double pk_4_b2_out;
    double *ddpk_4_b2;
    class_alloc(ddpk_4_b2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P_4_b2,
                                          1,
                                          ddpk_4_b2,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_4_b2,
                                                ddpk_4_b2,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_4_b2_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            pk_l_4_b2[index_k] = pk_4_b2_out + 1.e7;
        }
        else {
            pk_l_4_b2[index_k] = epsilon_for_logs + 1.e7;
        }
    }
    free(ddpk_4_b2);
    free(f22_4_b2);
    free(P_4_b2);

    double pk_4_bG2_out;
    double *ddpk_4_bG2;
    class_alloc(ddpk_4_bG2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P_4_bG2,
                                          1,
                                          ddpk_4_bG2,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>=kminnew){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_4_bG2,
                                                ddpk_4_bG2,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &pk_4_bG2_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            pk_l_4_bG2[index_k] = pk_4_bG2_out + 1.e7;
        }
        else {
            pk_l_4_bG2[index_k] = epsilon_for_logs + 1.e7;
        }
    }
    free(ddpk_4_bG2);
    free(f22_4_bG2);
    free(P_4_bG2);
    
    
    double fg2_out=0;
    double *ddpk_IFG2_0b1_mmm;
    class_alloc(ddpk_IFG2_0b1_mmm,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_IFG2_0b1_x,1,ddpk_IFG2_0b1_mmm,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    int last_index_2=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,P_IFG2_0b1_x,ddpk_IFG2_0b1_mmm,1,pnlpt->k[index_k],&last_index_2,&fg2_out,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_IFG2_0b1[index_k] = fg2_out + 1.e7;
            //    printf("%lf\n", fg2_out);
        }
        else{
            pk_IFG2_0b1[index_k] = exp(lnpk_l[index_k])*0. + 1.e7;
        }
     //   printf("%lf\n",pk_Id2d2_2[index_k]);
      //        printf("%le %le %le\n", pk_IFG2_0b1[index_k], pk_Id2d2_2[index_k],pk_Id2d2_4[index_k]);
    }
    free(ddpk_IFG2_0b1_mmm);
    free(P_IFG2_new);
    
    
    double fg2_out_2=0.;
    double *ddpk_IFG2_0;
    class_alloc(ddpk_IFG2_0,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,Nmax,P_IFG2_0,1,ddpk_IFG2_0,_SPLINE_NATURAL_,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,Nmax,P_IFG2_0,ddpk_IFG2_0,1,pnlpt->k[index_k],&last_index,&fg2_out_2,1,pnlpt->error_message),pnlpt->error_message,pnlpt->error_message);
            pk_IFG2_0[index_k] = fg2_out_2 + 1.e7;
        }
        else{
            pk_IFG2_0[index_k] = exp(lnpk_l[index_k])*f/3.*0. + 1.e7;
        }
        //      printf("%lf\n", pk_IFG2_0[index_k]);
    }
    free(ddpk_IFG2_0);
    
    double fg2_out_3=0.;
    double *ddpk_IFG2_2;
    class_alloc(ddpk_IFG2_2,sizeof(double)*Nmax,pnlpt->error_message);
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P_IFG2_2,
                                          1,
                                          ddpk_IFG2_2,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    last_index=0;
    for (index_k=0; index_k < pnlpt->k_size; index_k++){
        if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P_IFG2_2,
                                                ddpk_IFG2_2,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index,
                                                &fg2_out_3,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
  //          printf("%lf\n", fg2_out_3);
            pk_IFG2_2[index_k] = fg2_out_3 + 1.e7;
        }
        else{
            pk_IFG2_2[index_k] = exp(lnpk_l[index_k])*2*f/3.*0. + 1.e7;
        }
   //     printf("%le %le %le\n", pk_IFG2_0b1[index_k], pk_IFG2_0[index_k],pk_IFG2_2[index_k]);

        
    }
    
    free(ddpk_IFG2_2);
    
    
    free(P_IFG2_0b1_x);
    free(P_IFG2_0);
    free(P_IFG2_2);
    
    free(P_Id2G2_2);
    free(P_Id2d2_2);
    free(P_IG2G2_2);
    free(P_Id2d2_4);
    free(P_Id2G2_4);
    free(P_IG2G2_4);
    free(P_4_b1b2);
    free(P_4_b1bG2);
    
}//end of RSD conditional expression
        
       /* double kmaxnew = kdisc[Nmax-2];
        double kminnew = kdisc[1];*/
        
        double pk_Id2d2_out;
        double *ddpk_PId2d2;
        class_alloc(ddpk_PId2d2,sizeof(double)*Nmax,pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_Id2d2,
                                              1,
                                              ddpk_PId2d2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        
        last_index=0;
        for (index_k=0; index_k < pnlpt->k_size; index_k++){
            
            if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P_Id2d2,
                                                    ddpk_PId2d2,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_Id2d2_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);
                pk_Id2d2[index_k] = pk_Id2d2_out + 1.e7;
            }
            else {
                pk_Id2d2[index_k] = epsilon_for_logs + 1.e7;
            }
        }
        free(ddpk_PId2d2);
        free(f22_Id2d2);
        free(P_Id2d2);

        double pk_Id2G2_out;
        double *ddpk_Id2G2;
        class_alloc(ddpk_Id2G2,sizeof(double)*Nmax,pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_Id2G2,
                                              1,
                                              ddpk_Id2G2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index=0;
        for (index_k=0; index_k < pnlpt->k_size; index_k++){
            
            if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
                
                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P_Id2G2,
                                                    ddpk_Id2G2,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_Id2G2_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);
                pk_Id2G2[index_k] = pk_Id2G2_out + 1.e7;
            }
            
            else {
                pk_Id2G2[index_k] = epsilon_for_logs + 1.e7;
            }
        }
        free(ddpk_Id2G2);
        
        
        free(f22_Id2G2);
        free(P_Id2G2);
        
        double pk_IG2G2_out;
        double *ddpk_IG2G2;
        class_alloc(ddpk_IG2G2,sizeof(double)*Nmax,pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_IG2G2,
                                              1,
                                              ddpk_IG2G2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index=0;
        for (index_k=0; index_k < pnlpt->k_size; index_k++){
            
            if (pnlpt->k[index_k]<=kmaxnew && pnlpt->k[index_k]>= kminnew){
                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P_IG2G2,
                                                    ddpk_IG2G2,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_IG2G2_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);
                pk_IG2G2[index_k] = pk_IG2G2_out + 1.e7;
            }
            else {
                pk_IG2G2[index_k] = epsilon_for_logs + 1.e7;
            }
        }
        free(ddpk_IG2G2);
        free(f22_IG2G2);
        free(P_IG2G2);

        
        free(ddpk_IFG2);
        free(f13_IFG2);
        free(P_IFG2);
     
        //old here
     
free(etam2);
free(cmsym2);
free(x2);
        
        /*
        int* addressOfX = &index_k;
        printf("%p\n", &index_k);
         */ 

} // end of bias conditional expression
     
     else {
         
         if (pnlpt->nonlinear_pt_verbose > 0)
         printf("No bias tracers requested.\n");
         
        double epsilon_for_logs = 1.*pow(10.,-6.);
         
         for (index_k=0; index_k < pnlpt->k_size; index_k++){
             pk_Id2d2[index_k] = epsilon_for_logs ;
             pk_Id2[index_k] = epsilon_for_logs;
             pk_IG2[index_k] = epsilon_for_logs ;
             pk_Id2G2[index_k] = epsilon_for_logs ;
             pk_IG2G2[index_k] = epsilon_for_logs;
             pk_IFG2[index_k] = epsilon_for_logs ;
         }
         
     }
int end1=clock();
     
if (pnlpt->nonlinear_pt_verbose > 0)     
printf("All matrices are calculated in %d mus\n",end1-start1);

free(js);
free(kdisc);
free(Pbin);
free(y);
     
free(Pw);
free(Pnw);
free(dd_Pnw);
free(dd_Pw);
     
free(P10b1);
free(P10);
free(P12);
     
return _SUCCESS_;
}

