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
    struct precision *ppr,
    struct background *pba,
    struct thermo *pth,
    //struct transfer * ptr //GC!
    struct perturbs *ppt,
    struct nonlinear_pt *pnlpt)
{
    int index_k, index_k_output, index_mode;
    double k, k_min = 0., k_rec, step, tau1;
    double k_max_cmb;
    double k_max_cl;
    double k_max = 0.;
    double scale2;
    double *tmp_k_list;
    int newk_size, index_newk, add_k_output_value;

    /** Summary: */

    class_test(ppr->k_step_transition == 0.,
               pnlpt->error_message,
               "stop to avoid division by zero");

    class_test(pth->rs_rec == 0.,
               pnlpt->error_message,
               "stop to avoid division by zero");

    /** - scalar modes */

    if (ppt->has_scalars == _TRUE_)
    {

        k_min = ppr->k_min_tau0 / pba->conformal_age;

        /* first value */
        if (pba->sgnK == 0)
        {
            /* K<0 (flat)  : start close to zero */
            k_min = ppr->k_min_tau0 / pba->conformal_age;
        }

        else if (pba->sgnK == -1)
        {
            // K<0 (open)  : start close to sqrt(-K)
            // (in transfer modules, for scalars, this will correspond to q close to zero;
            // for vectors and tensors, this value is even smaller than the minimum necessary value)
            k_min = sqrt(-pba->K + pow(ppr->k_min_tau0 / pba->conformal_age / pth->angular_rescaling, 2));
        }
        else if (pba->sgnK == 1)
        {
            /* K>0 (closed): start from q=sqrt(k2+(1+m)K) equal to 3sqrt(K), i.e. k=sqrt((8-m)K) */
            k_min = sqrt((8. - 1.e-4) * pba->K);
        }

        // printf("k_min=%le\n",k_min);
        /** - --> find k_max (as well as k_max_cmb[ppt->index_md_scalars], k_max_cl[ppt->index_md_scalars]) */

        k_rec = 2. * _PI_ / pth->rs_rec; /* comoving scale corresponding to sound horizon at recombination */

        k_max_cmb = k_min;
        k_max_cl = k_min;
        k_max = k_min;

        /* find k_max: */

        if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_))
            k_max = MAX(k_max, ppt->k_max_for_pk);

        if (ppt->has_nl_corrections_based_on_delta_m == _TRUE_)
            k_max = MAX(k_max, ppr->halofit_min_k_max);

        // printf("k_max=%le\n",k_max);

        /** - --> test that result for k_min, k_max make sense */

        class_test(k_min < 0.,
                   pnlpt->error_message,
                   "buggy definition of k_min");

        class_test(k_max < 0.,
                   pnlpt->error_message,
                   "buggy definition of k_max");

        class_test(k_max < k_min,
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
                    ((int)((k_max_cmb - k_min) / k_rec / MIN(ppr->k_step_super, ppr->k_step_sub)) +
                     (int)(MAX(ppr->k_per_decade_for_pk, ppr->k_per_decade_for_bao) * log(k_max / k_min) / log(10.)) + 3) *
                        sizeof(double),
                    pnlpt->error_message);

        /* first value */

        index_k = 0;
        k = k_min;
        pnlpt->k[index_k] = k;
        index_k++;

        /* values until k_max_cmb[ppt->index_md_scalars] */

        while (k < k_max_cmb)
        {

            /* the linear step is not constant, it has a step-like shape,
         centered around the characteristic scale set by the sound
         horizon at recombination (associated to the comoving wavenumber
         k_rec) */

            step = (ppr->k_step_super + 0.5 * (tanh((k - k_rec) / k_rec / ppr->k_step_transition) + 1.) * (ppr->k_step_sub - ppr->k_step_super)) * k_rec;

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

            scale2 = pow(pba->a_today * pba->H0, 2) + fabs(pba->K);

            step *= (k * k / scale2 + 1.) / (k * k / scale2 + 1. / ppr->k_step_super_reduction);

            class_test(step / k < ppr->smallest_allowed_variation,
                       pnlpt->error_message,
                       "k step =%e < machine precision : leads either to numerical error or infinite loop",
                       step * k_rec);

            k += step;

            class_test(k <= pnlpt->k[index_k - 1],
                       pnlpt->error_message,
                       "consecutive values of k should differ and should be in growing order");

            pnlpt->k[index_k] = k;

            index_k++;
        }

        pnlpt->k_size_cmb = index_k;

        /* values until k_max_cl[ppt->index_md_scalars] */

        while (k < k_max_cl)
        {

            k *= pow(10., 1. / (ppr->k_per_decade_for_pk + (ppr->k_per_decade_for_bao - ppr->k_per_decade_for_pk) * (1. - tanh(pow((log(k) - log(ppr->k_bao_center * k_rec)) / log(ppr->k_bao_width), 4)))));

            pnlpt->k[index_k] = k;
            index_k++;
        }

        pnlpt->k_size_cl = index_k;

        /* values until k_max */

        while (k < k_max)
        {

            k *= pow(10., 1. / (ppr->k_per_decade_for_pk + (ppr->k_per_decade_for_bao - ppr->k_per_decade_for_pk) * (1. - tanh(pow((log(k) - log(ppr->k_bao_center * k_rec)) / log(ppr->k_bao_width), 4)))));

            pnlpt->k[index_k] = k;
            index_k++;
        }

        pnlpt->k_size = index_k;

        class_realloc(pnlpt->k,
                      pnlpt->k,
                      pnlpt->k_size * sizeof(double),
                      pnlpt->error_message);
    }

    /** - vector modes skipped */

    /** - tensor modes skipped */

    /** - If user asked for k_output_values, add those to all k lists: */
    if (ppt->k_output_values_num > 0)
    {
        /* Allocate storage */
        class_alloc(ppt->index_k_output_values, sizeof(double) * 1 * ppt->k_output_values_num, pnlpt->error_message);

        /** - --> Find indices in ppt->k[index_md] corresponding to 'k_output_values'.
        We are assuming that ppt->k is sorted and growing, and we have made sure
        that ppt->k_output_values is also sorted and growing.*/
        //for (index_mode=0; index_mode<ppt->md_size; index_mode++){

        newk_size = pnlpt->k_size + ppt->k_output_values_num;

        class_alloc(tmp_k_list, sizeof(double) * newk_size, pnlpt->error_message);

        index_k = 0;
        index_k_output = 0;
        for (index_newk = 0; index_newk < newk_size; index_newk++)
        {
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

            if (add_k_output_value == _TRUE_)
            {
                tmp_k_list[index_newk] = ppt->k_output_values[index_k_output];
                ppt->index_k_output_values[0 * ppt->k_output_values_num + index_k_output] = index_newk;
                index_k_output++;
            }
            else
            {
                tmp_k_list[index_newk] = pnlpt->k[index_k];
                index_k++;
            }
        }

        free(pnlpt->k);
        pnlpt->k = tmp_k_list;
        pnlpt->k_size = newk_size;

        index_k = newk_size - 1;
        while (pnlpt->k[index_k] > k_max_cl)
            index_k--;
        pnlpt->k_size_cl = MIN(index_k + 2, pnlpt->k_size);

        index_k = newk_size - 1;
        while (pnlpt->k[index_k] > k_max_cmb)
            index_k--;
        pnlpt->k_size_cmb = MIN(index_k + 2, pnlpt->k_size);

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

    /** - finally, find the global k_min and k_max for the ensemble of all modes scalars, vectors, tensors) */

    pnlpt->k_min = _HUGE_;
    pnlpt->k_max = 0.;
    if (ppt->has_scalars == _TRUE_)
    {
        pnlpt->k_min = MIN(pnlpt->k_min, pnlpt->k[0]);                 /* first value, inferred from perturbations structure */
        pnlpt->k_max = MAX(pnlpt->k_max, pnlpt->k[pnlpt->k_size - 1]); /* last value, inferred from perturbations structure */
    }

    //free(k_max_cmb);
    //free(k_max_cl);

    return _SUCCESS_;
}

//GC -> the function above can remain totally unchanged... This below instead must be changed. It allocates stuff, AND ALSO CALLS nonlinear_pt_loop, which is the main function...

int nonlinear_pt_init(
    struct precision *ppr,
    struct background *pba,
    struct thermo *pth,
    struct perturbs *ppt,
    struct primordial *ppm,
    struct nonlinear_pt *pnlpt)
{

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
    //double *tk_l; //GC!
    //double *lntk_l; //GC!
    //double *ddlntk_l; //GC!

    /* -> this we need: indeed, this function CALLS nonlinear_pt_pk_l, which uses these variables to get the linear PS. I will likely create a new function that extracts only the primordial PS. In any case, I would still need these... */

    double *pPRIMk_l;     //GC!
    double *lnpPRIMk_l;   //GC!
    double *ddlnpPRIMk_l; //GC!

    /**/

    //GC -> then I need to allocate the rest, the outputs...

    //double *pk_nl_fNL_out; //GC: this shouldn't be here...

    double *pk_l_fNL_0_vv;
    double *pk_l_fNL_0_vd;
    double *pk_l_fNL_0_dd;
    double *pk_l_fNL_2_vv;
    double *pk_l_fNL_2_vd;
    double *pk_l_fNL_2_dd;
    double *pk_l_fNL_4_vv;
    double *pk_l_fNL_4_vd;
    double *pk_l_fNL_4_dd;

    //GC: ORTHOGONAL -- start

    double *pk_l_fNL_0_vv_ortho;
    double *pk_l_fNL_0_vd_ortho;
    double *pk_l_fNL_0_dd_ortho;
    double *pk_l_fNL_2_vv_ortho;
    double *pk_l_fNL_2_vd_ortho;
    double *pk_l_fNL_2_dd_ortho;
    double *pk_l_fNL_4_vv_ortho;
    double *pk_l_fNL_4_vd_ortho;
    double *pk_l_fNL_4_dd_ortho;

    //GC: ORTHOGONAL -- finish

    double *pk12_l_0_b1b2;
    double *pk12_l_0_b2;
    double *pk12_l_0_b1bG2;
    double *pk12_l_0_bG2;
    double *pk12_l_2_b1b2;
    double *pk12_l_2_b2;
    double *pk12_l_2_b1bG2;
    double *pk12_l_2_bG2;
    double *pk12_l_4_b1b2;
    double *pk12_l_4_b2;
    double *pk12_l_4_b1bG2;
    double *pk12_l_4_bG2;

    //GC: ORTHOGONAL -- start

    double *pk12_l_0_b1b2_ortho;
    double *pk12_l_0_b2_ortho;
    double *pk12_l_0_b1bG2_ortho;
    double *pk12_l_0_bG2_ortho;
    double *pk12_l_2_b1b2_ortho;
    double *pk12_l_2_b2_ortho;
    double *pk12_l_2_b1bG2_ortho;
    double *pk12_l_2_bG2_ortho;
    double *pk12_l_4_b1b2_ortho;
    double *pk12_l_4_b2_ortho;
    double *pk12_l_4_b1bG2_ortho;
    double *pk12_l_4_bG2_ortho;

    //GC: ORTHOGONAL -- finish

    //GC -> shouldn't I allocate also the ones for matter only? Yes, it is pk_nl_fNL, which takes the place of pk_nl...

    double *pk_nl_fNL;

    //GC -> as it is now I also need the following ->

    double *pk_fNLd2;
    double *pk_fNLG2;

    //GC: ORTHOGONAL -- start

    double *pk_nl_fNL_ortho;

    double *pk_fNLd2_ortho;
    double *pk_fNLG2_ortho;

    //GC: ORTHOGONAL -- finish

    short print_warning = _FALSE_;
    //double * pvecback;
    int last_index;
    double a, z;

    last_index = 0;

    /** Summary
   *
   * (a) First deal with the case where non non-linear corrections requested */

    if (pnlpt->method == nlpt_none)
    {
        if (pnlpt->nonlinear_pt_verbose > 0)
            printf("No non-linear spectra requested. Nonlinear PT module skipped.\n");
    }

    /** (b) Compute the non-linear spectrum */

    else if (pnlpt->method == nlpt_spt)
    {
        if (pnlpt->nonlinear_pt_verbose > 0)
            printf("Computing non-linear matter and bias power spectra in perturbation theory.\n");

        /** - copy list of (k,tau) from perturbation module */

        pnlpt->index_md_scalars = ppt->index_md_scalars;
        pnlpt->ic_size = ppt->ic_size[ppt->index_md_scalars];
        pnlpt->tp_size = ppt->tp_size[ppt->index_md_scalars];

        if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_)
        {

            class_call(perturb_get_k_list_nl(ppr,
                                             pba,
                                             pth,
                                             ppt,
                                             pnlpt),
                       pnlpt->error_message,
                       pnlpt->error_message);
        }
        else
        {

            pnlpt->k_size = ppt->k_size[ppt->index_md_scalars];
            class_alloc(pnlpt->k, pnlpt->k_size * sizeof(double), pnlpt->error_message);
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                pnlpt->k[index_k] = ppt->k[ppt->index_md_scalars][index_k];
            }
        }

        pnlpt->ln_k_size = pnlpt->k_size;
        class_alloc(pnlpt->ln_k, pnlpt->ln_k_size * sizeof(double), pnlpt->error_message);
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {
            pnlpt->ln_k[index_k] = log(pnlpt->k[index_k]);
        }

        pnlpt->tau_size = ppt->tau_size;
        class_alloc(pnlpt->tau, pnlpt->tau_size * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->ln_tau, pnlpt->tau_size * sizeof(double), pnlpt->error_message);
        for (index_tau = 0; index_tau < pnlpt->tau_size; index_tau++)
        {
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

        //GC -> this is some place where actually I need to allocate new stuff... [[[The stuff above is unchanged...]]]

        /*Needed for transfers*/
        class_alloc(pnlpt->nl_corr_density, pnlpt->tau_size * ppt->k_size[pnlpt->index_md_scalars] * sizeof(double), pnlpt->error_message);

        class_alloc(pk_l, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_0_vv, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_0_vd, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_0_dd, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        class_alloc(pk_l_2_vv, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_2_vd, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_2_dd, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        class_alloc(pk_l_4_vv, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_4_vd, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_4_dd, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        class_alloc(pk_l_0_b1b2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_0_b2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_0_b1bG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_0_bG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        class_alloc(pk_l_2_b1b2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_2_b2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_2_b1bG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_2_bG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        class_alloc(pk_l_4_b2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_4_bG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_4_b1b2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_4_b1bG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        class_alloc(pk_nl, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_Id2d2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_Id2d2_2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_Id2d2_4, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        class_alloc(pk_Id2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_IG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_Id2G2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_Id2G2_2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_Id2G2_4, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_IG2G2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_IG2G2_2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_IG2G2_4, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_IFG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_IFG2_0b1, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_IFG2_0, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_IFG2_2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_CTR, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_CTR_0, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_CTR_2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_CTR_4, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_Tree, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_Tree_0_vv, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_Tree_0_vd, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_Tree_0_dd, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_Tree_2_vv, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_Tree_2_vd, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_Tree_4_vv, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(lnk_l, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(lnpk_l, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(ddlnpk_l, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        //class_alloc(tk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message); //GC!
        //class_alloc(lntk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message); //GC!
        //class_alloc(ddlntk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message); //GC!

        //GC!

        class_alloc(pPRIMk_l, pnlpt->k_size * sizeof(double), pnlpt->error_message);     //GC!
        class_alloc(lnpPRIMk_l, pnlpt->k_size * sizeof(double), pnlpt->error_message);   //GC!
        class_alloc(ddlnpPRIMk_l, pnlpt->k_size * sizeof(double), pnlpt->error_message); //GC!
        /* -> see line ~ 370... Notice that pk_l is allocated as first, for some reason... */

        class_alloc(pk_l_fNL_0_vv, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_fNL_0_vd, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_fNL_0_dd, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        class_alloc(pk_l_fNL_2_vv, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_fNL_2_vd, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_fNL_2_dd, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        class_alloc(pk_l_fNL_4_vv, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_fNL_4_vd, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_fNL_4_dd, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pk_l_fNL_0_vv_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_fNL_0_vd_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_fNL_0_dd_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        class_alloc(pk_l_fNL_2_vv_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_fNL_2_vd_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_fNL_2_dd_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        class_alloc(pk_l_fNL_4_vv_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_fNL_4_vd_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_l_fNL_4_dd_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        class_alloc(pk12_l_0_b1b2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_0_b2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_0_b1bG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_0_bG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_2_b1b2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_2_b2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_2_b1bG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_2_bG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_4_b1b2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_4_b2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_4_b1bG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_4_bG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pk12_l_0_b1b2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_0_b2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_0_b1bG2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_0_bG2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_2_b1b2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_2_b2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_2_b1bG2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_2_bG2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_4_b1b2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_4_b2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_4_b1bG2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk12_l_4_bG2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        class_alloc(pk_nl_fNL, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_fNLd2, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_fNLG2, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pk_nl_fNL_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_fNLd2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);
        class_alloc(pk_fNLG2_ortho, pnlpt->k_size * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        //GC!

        //  double tau_req[pnlpt->z_pk_num];
        //  double deltatau[pnlpt->z_pk_num];

        class_alloc(tau_req, pnlpt->z_pk_num * sizeof(double), pnlpt->error_message);
        //class_alloc(deltatau,pnlpt->z_pk_num * sizeof(double),pnlpt->error_message);

        int i_z = 0;
        for (i_z = 0; i_z < pnlpt->z_pk_num; i_z++)
        {
            //                        printf("redshift requested %lf\n", pnlpt->z_pk[i]);
            class_call(background_tau_of_z(pba, pnlpt->z_pk[i_z], &tau_req[i_z]),
                       pba->error_message,
                       pnlpt->error_message);
        }

        //GC: this above remains the same...

        /** Inputing the PT matrices */

        //GC: as I thought, the basic matrix is only I(\nu1,nu2), while the not-basic is the one that includes the kernels for matter only. I can check in the python script that the factor with 196 is the one for P22, and this is for the M22oneline, NOT BASIC... So I only need M12oneline and M12basiconeline... No, I should actually only need M12oneline to extract from python... THERE IS ONE DIFFERENCE/PROBLEM -> for some reason, right now it seems that Misha uses M22oneline, NOT the basic (i.e. I(n1,n2)), as a baseline ("pnlpt->M22_oneline_complex[count]") for the matrices below. While I have "FullSimplified" with respect to I(n1,n2)... So I have to be careful -> most likely I want some equivalent of "pnlpt->M22basic_oneline_complex[count]". I am almost 100% sure of this. The fact that the not-basic matrix is the one for P22 can be checked on python, and on python I can also check better if the basic is only I(nu1,nu2) -> in this case, I would only need python to compute the non-basic P12... Ok -> I checked that M22 in python is exactly the matrix for 2*F2^2. Then I use my new matrix there for 2*F2*SNG and I call it M12oneline_N128.dat (and higher number of points). There is one final confusion -> why do they have a packed version ONLY for the 256 case? Does it mean that they actually only run the 256 case and didn't bother to run the others? There is another problem -> and maybe this is the crux of the matter => the bias! First, I checked that indeed the bias must enter also in the matrix itself, so it is correct that I need to supply it in the python script, along with the number of points (ppr->nmax_nlpt). Then, the point is that for \delta^2\delta^2, that I would have used as the basic matrix, I need a different bias (all have -0.8 apart from this, that has -1.25). Something similar happens for the biased tracers -> so, to avoid having to compute TWO "basic" matrices, one just defines the baseline as the one for P22 (or P12 in my case). If so, I need to compute again the function of (n1,n2) that multiplies them below, and I reserve a separate matrix for the \delta^2 contribution to P12. For this reason, I wait to write down the hexadecapole that I am currently missing...

        char file2openM22[256];
        char file2openM13[256];
        char file2openM22basic[256];
        char file2openM13basic[256];

        //GC!
        char file2openM12_matter[256];

        //GC: ORTHOGONAL -- start

        char file2openM12_matter_ortho[256];

        //GC: ORTHOGONAL -- finish

        /*-----------------------------------------------*/
        /**/
        //GC: CAREFUL, I do not add the ortho subscript to make it quicker!!! Since I only need to check the orthogonal ones...
        char file2openM12_matter_multipoles_vv0_f2[256];
        char file2openM12_matter_multipoles_vv0_f3[256];
        char file2openM12_matter_multipoles_vd0_f1[256];
        char file2openM12_matter_multipoles_vd0_f2[256];
        char file2openM12_matter_multipoles_dd0_f0[256];
        char file2openM12_matter_multipoles_dd0_f1[256];
        char file2openM12_matter_multipoles_vv2_f3[256];
        char file2openM12_matter_multipoles_vd2_f2[256];
        char file2openM12_matter_multipoles_vv4_f3[256];
        char file2openM12_matter_multipoles_vd4_f2[256];
        /**/
        /*-----------------------------------------------*/
        /**/
        //GC: CAREFUL, I do not add the ortho subscript to make it quicker!!! Since I only need to check the orthogonal ones...
        char file2openM12_matter_multipoles_vv0_f2_ortho[256];
        char file2openM12_matter_multipoles_vv0_f3_ortho[256];
        char file2openM12_matter_multipoles_vd0_f1_ortho[256];
        char file2openM12_matter_multipoles_vd0_f2_ortho[256];
        char file2openM12_matter_multipoles_dd0_f0_ortho[256];
        char file2openM12_matter_multipoles_dd0_f1_ortho[256];
        char file2openM12_matter_multipoles_vv2_f3_ortho[256];
        char file2openM12_matter_multipoles_vd2_f2_ortho[256];
        char file2openM12_matter_multipoles_vv4_f3_ortho[256];
        char file2openM12_matter_multipoles_vd4_f2_ortho[256];
        /**/
        /*-----------------------------------------------*/
        char file2openM12_matter_mu_powers_vd2_f1[256];
        char file2openM12_matter_mu_powers_vd2_f2[256];
        char file2openM12_matter_mu_powers_dd2_f1[256];
        char file2openM12_matter_mu_powers_vv4_f2[256];
        char file2openM12_matter_mu_powers_vd4_f2[256];
        char file2openM12_matter_mu_powers_vv6_f3[256];
        //GC: ORTHOGONAL -- start
        char file2openM12_matter_mu_powers_vd2_f1_ortho[256];
        char file2openM12_matter_mu_powers_vd2_f2_ortho[256];
        char file2openM12_matter_mu_powers_dd2_f1_ortho[256];
        char file2openM12_matter_mu_powers_vv4_f2_ortho[256];
        char file2openM12_matter_mu_powers_vd4_f2_ortho[256];
        char file2openM12_matter_mu_powers_vv6_f3_ortho[256];
        //GC: ORTHOGONAL -- finish
        /*-----------------------------------------------*/
        /**/
        //GC: CAREFUL, I do not add the ortho subscript to make it quicker!!! Since I only need to check the orthogonal ones...
        char file2openM12_bias_real_space_b2[256];
        char file2openM12_bias_real_space_bG2[256];
        /**/
        /*-----------------------------------------------*/
        /**/
        //GC: CAREFUL, I do not add the ortho subscript to make it quicker!!! Since I only need to check the orthogonal ones...
        char file2openM12_bias_real_space_b2_ortho[256];
        char file2openM12_bias_real_space_bG2_ortho[256];
        /**/
        /*-----------------------------------------------*/
        char file2openM12_bias_multipoles_b2_vv0_f1[256];
        char file2openM12_bias_multipoles_bG2_vv0_f1[256];
        //GC: ORTHOGONAL -- start
        char file2openM12_bias_multipoles_b2_vv0_f1_ortho[256];
        char file2openM12_bias_multipoles_bG2_vv0_f1_ortho[256];
        //GC: ORTHOGONAL -- finish

        //GC!

        //int startmatrix=clock();

        //GC -> no need to store in RAM. They are just a bunch of arrays. Write them in a structure and have this code read those... Essentially, instead of e.g. "fscanf(myFile22, "%lf", &pnlpt->M22_oneline[index_M22]);", there would be a command that says pnlpt->M22_oneline[index_M22] is equal to the other array... Actually, this is stupid. It does not make any sense to have a separate file, except for readability... One could literally write the matrices below here... However, I do not know how to do it so this point is moot so far... You can have struct matrices *mtrx or something, and then write pnlpt->M22_oneline[index_M22] = mtrx->something[index_M22] + _Complex_I * whatever... One could look at input.c for how it is done... But for now, I think I will just import the matrices -> since we fix cosmology (but not biases) the loss of time is likely not too tragic, since we need to explore a smaller parameter space... Also, we only need the AP with IR-resummation in redshift space, so not all matrices need to be imported... I need to import the non-IR-resummed multipoles just to check the numerics...
        
        //GC - SWITCH HERE!
        //GC: the quickest way is to remove the opening of the files first... I think also to just keep allocating. Only now they will have junk inside...
        
        //GC -> define a master variable "switch" here, as a string, to then run tests...
        
        int SWITCH_index = 0;
        
        if (pnlpt->fNL_equil_ortho_switch == fNL_equil_ortho_yes)
            
        {
            
            SWITCH_index = 1;
            
        }
        
        //printf("%d\n",SWITCH_index); //GC - SWITCH!

        if (ppr->nmax_nlpt == 128)
        {
            sprintf(file2openM22, "%s/pt_matrices/M22oneline_N128.dat", __CLASSDIR__);
            sprintf(file2openM13, "%s/pt_matrices/M13oneline_N128.dat", __CLASSDIR__);
            sprintf(file2openM22basic, "%s/pt_matrices/M22basiconeline_N128.dat", __CLASSDIR__);
            sprintf(file2openM13basic, "%s/pt_matrices/IFG2oneline_N128.dat", __CLASSDIR__);
            //sprintf(file2openM12,"%s/pt_matrices/M12oneline_N128.dat",__CLASSDIR__);    //GC!
            
            //GC: SWITCH HERE!!!
            
            if (SWITCH_index == 1)
                
            {
                
                //printf("%d\n",SWITCH_index); //GC - SWITCH!
            
            /*-----------------------------------------------*/
            sprintf(file2openM12_matter, "%s/pt_matrices/compute_matrices_python/M12oneline_N128-matter.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- start
            sprintf(file2openM12_matter_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/M12oneline_N128-matter-orthogonal.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
            /**/
            sprintf(file2openM12_matter_multipoles_vv0_f2, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N128-matter_multipoles-vv0_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv0_f3, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N128-matter_multipoles-vv0_f3.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd0_f1, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N128-matter_multipoles-vd0_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd0_f2, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N128-matter_multipoles-vd0_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_dd0_f0, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N128-matter_multipoles-dd0_f0.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_dd0_f1, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N128-matter_multipoles-dd0_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv2_f3, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N128-matter_multipoles-vv2_f3.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd2_f2, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N128-matter_multipoles-vd2_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv4_f3, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N128-matter_multipoles-vv4_f3.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd4_f2, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N128-matter_multipoles-vd4_f2.txt.dat", __CLASSDIR__);
            /**/
            /*-----------------------------------------------*/
            /**/
            //GC: ORTHOGONAL -- start
            sprintf(file2openM12_matter_multipoles_vv0_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N128-matter_multipoles-vv0_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv0_f3_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N128-matter_multipoles-vv0_f3.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd0_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N128-matter_multipoles-vd0_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd0_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N128-matter_multipoles-vd0_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_dd0_f0_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N128-matter_multipoles-dd0_f0.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_dd0_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N128-matter_multipoles-dd0_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv2_f3_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N128-matter_multipoles-vv2_f3.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd2_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N128-matter_multipoles-vd2_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv4_f3_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N128-matter_multipoles-vv4_f3.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd4_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N128-matter_multipoles-vd4_f2.txt-orthogonal.dat", __CLASSDIR__);
            /**/
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
            sprintf(file2openM12_matter_mu_powers_vd2_f1, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N128-matter_mu_powers-vd2_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vd2_f2, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N128-matter_mu_powers-vd2_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_dd2_f1, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N128-matter_mu_powers-dd2_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vv4_f2, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N128-matter_mu_powers-vv4_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vd4_f2, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N128-matter_mu_powers-vd4_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vv6_f3, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N128-matter_mu_powers-vv6_f3.txt.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- start
            sprintf(file2openM12_matter_mu_powers_vd2_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N128-matter_mu_powers-vd2_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vd2_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N128-matter_mu_powers-vd2_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_dd2_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N128-matter_mu_powers-dd2_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vv4_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N128-matter_mu_powers-vv4_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vd4_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N128-matter_mu_powers-vd4_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vv6_f3_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N128-matter_mu_powers-vv6_f3.txt-orthogonal.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
            /**/
            sprintf(file2openM12_bias_real_space_b2, "%s/pt_matrices/compute_matrices_python/bias_real_space/M12oneline_N128-bias_real_space-b2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_bias_real_space_bG2, "%s/pt_matrices/compute_matrices_python/bias_real_space/M12oneline_N128-bias_real_space-G2.txt.dat", __CLASSDIR__);
            /**/
            /*-----------------------------------------------*/
            //GC: ORTHOGONAL -- start
            /**/
            sprintf(file2openM12_bias_real_space_b2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/bias_real_space/M12oneline_N128-bias_real_space-b2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_bias_real_space_bG2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/bias_real_space/M12oneline_N128-bias_real_space-G2.txt-orthogonal.dat", __CLASSDIR__);
            /**/
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
            sprintf(file2openM12_bias_multipoles_b2_vv0_f1, "%s/pt_matrices/compute_matrices_python/bias_multipoles/M12oneline_N128-bias_multipoles-b2_vv0_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_bias_multipoles_bG2_vv0_f1, "%s/pt_matrices/compute_matrices_python/bias_multipoles/M12oneline_N128-bias_multipoles-bG2_vv0_f1.txt.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- start
            sprintf(file2openM12_bias_multipoles_b2_vv0_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/bias_multipoles/M12oneline_N128-bias_multipoles-b2_vv0_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_bias_multipoles_bG2_vv0_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/bias_multipoles/M12oneline_N128-bias_multipoles-bG2_vv0_f1.txt-orthogonal.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
                
            }
        }
        else if (ppr->nmax_nlpt == 512)
        {
            sprintf(file2openM22, "%s/pt_matrices/M22oneline_N512.dat", __CLASSDIR__);
            sprintf(file2openM13, "%s/pt_matrices/M13oneline_N512.dat", __CLASSDIR__);
            sprintf(file2openM22basic, "%s/pt_matrices/M22basiconeline_N512.dat", __CLASSDIR__);
            sprintf(file2openM13basic, "%s/pt_matrices/IFG2oneline_N512.dat", __CLASSDIR__);
            //sprintf(file2openM12,"%s/pt_matrices/M12oneline_N512.dat",__CLASSDIR__);    //GC!
            
            if (SWITCH_index == 1)
                
            {
                
                //printf("%d\n",SWITCH_index); //GC - SWITCH!
            
            /*-----------------------------------------------*/
            sprintf(file2openM12_matter, "%s/pt_matrices/compute_matrices_python/M12oneline_N512-matter.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- start
            sprintf(file2openM12_matter_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/M12oneline_N512-matter-orthogonal.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
            /**/
            sprintf(file2openM12_matter_multipoles_vv0_f2, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N512-matter_multipoles-vv0_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv0_f3, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N512-matter_multipoles-vv0_f3.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd0_f1, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N512-matter_multipoles-vd0_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd0_f2, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N512-matter_multipoles-vd0_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_dd0_f0, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N512-matter_multipoles-dd0_f0.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_dd0_f1, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N512-matter_multipoles-dd0_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv2_f3, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N512-matter_multipoles-vv2_f3.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd2_f2, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N512-matter_multipoles-vd2_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv4_f3, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N512-matter_multipoles-vv4_f3.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd4_f2, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N512-matter_multipoles-vd4_f2.txt.dat", __CLASSDIR__);
            /**/
            /*-----------------------------------------------*/
            //GC: ORTHOGONAL -- start
            /**/
            sprintf(file2openM12_matter_multipoles_vv0_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N512-matter_multipoles-vv0_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv0_f3_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N512-matter_multipoles-vv0_f3.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd0_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N512-matter_multipoles-vd0_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd0_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N512-matter_multipoles-vd0_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_dd0_f0_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N512-matter_multipoles-dd0_f0.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_dd0_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N512-matter_multipoles-dd0_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv2_f3_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N512-matter_multipoles-vv2_f3.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd2_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N512-matter_multipoles-vd2_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv4_f3_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N512-matter_multipoles-vv4_f3.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd4_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N512-matter_multipoles-vd4_f2.txt-orthogonal.dat", __CLASSDIR__);
            /**/
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
            sprintf(file2openM12_matter_mu_powers_vd2_f1, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N512-matter_mu_powers-vd2_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vd2_f2, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N512-matter_mu_powers-vd2_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_dd2_f1, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N512-matter_mu_powers-dd2_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vv4_f2, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N512-matter_mu_powers-vv4_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vd4_f2, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N512-matter_mu_powers-vd4_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vv6_f3, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N512-matter_mu_powers-vv6_f3.txt.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- start
            sprintf(file2openM12_matter_mu_powers_vd2_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N512-matter_mu_powers-vd2_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vd2_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N512-matter_mu_powers-vd2_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_dd2_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N512-matter_mu_powers-dd2_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vv4_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N512-matter_mu_powers-vv4_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vd4_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N512-matter_mu_powers-vd4_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vv6_f3_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N512-matter_mu_powers-vv6_f3.txt-orthogonal.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
            /**/
            sprintf(file2openM12_bias_real_space_b2, "%s/pt_matrices/compute_matrices_python/bias_real_space/M12oneline_N512-bias_real_space-b2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_bias_real_space_bG2, "%s/pt_matrices/compute_matrices_python/bias_real_space/M12oneline_N512-bias_real_space-G2.txt.dat", __CLASSDIR__);
            /**/
            /*-----------------------------------------------*/
            //GC: ORTHOGONAL -- start
            /**/
            sprintf(file2openM12_bias_real_space_b2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/bias_real_space/M12oneline_N512-bias_real_space-b2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_bias_real_space_bG2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/bias_real_space/M12oneline_N512-bias_real_space-G2.txt-orthogonal.dat", __CLASSDIR__);
            /**/
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
            sprintf(file2openM12_bias_multipoles_b2_vv0_f1, "%s/pt_matrices/compute_matrices_python/bias_multipoles/M12oneline_N512-bias_multipoles-b2_vv0_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_bias_multipoles_bG2_vv0_f1, "%s/pt_matrices/compute_matrices_python/bias_multipoles/M12oneline_N512-bias_multipoles-bG2_vv0_f1.txt.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- start
            sprintf(file2openM12_bias_multipoles_b2_vv0_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/bias_multipoles/M12oneline_N512-bias_multipoles-b2_vv0_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_bias_multipoles_bG2_vv0_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/bias_multipoles/M12oneline_N512-bias_multipoles-bG2_vv0_f1.txt-orthogonal.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
            
            }
                
        }
        else
        {
            sprintf(file2openM22, "%s/pt_matrices/M22oneline_N256_packed.dat", __CLASSDIR__);
            sprintf(file2openM13, "%s/pt_matrices/M13oneline_N256.dat", __CLASSDIR__);
            sprintf(file2openM22basic, "%s/pt_matrices/M22basiconeline_N256_packed.dat", __CLASSDIR__);
            sprintf(file2openM13basic, "%s/pt_matrices/IFG2oneline_N256.dat", __CLASSDIR__);
            //sprintf(file2openM12,"%s/pt_matrices/M12oneline_N256.dat",__CLASSDIR__);    //GC!
            
            if (SWITCH_index == 1)
                
            {
                
                //printf("%d\n",SWITCH_index); //GC - SWITCH!
            
            /*-----------------------------------------------*/
            sprintf(file2openM12_matter, "%s/pt_matrices/compute_matrices_python/M12oneline_N256-matter.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- start
            sprintf(file2openM12_matter_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/M12oneline_N256-matter-orthogonal.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
            /**/
            sprintf(file2openM12_matter_multipoles_vv0_f2, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N256-matter_multipoles-vv0_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv0_f3, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N256-matter_multipoles-vv0_f3.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd0_f1, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N256-matter_multipoles-vd0_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd0_f2, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N256-matter_multipoles-vd0_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_dd0_f0, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N256-matter_multipoles-dd0_f0.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_dd0_f1, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N256-matter_multipoles-dd0_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv2_f3, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N256-matter_multipoles-vv2_f3.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd2_f2, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N256-matter_multipoles-vd2_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv4_f3, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N256-matter_multipoles-vv4_f3.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd4_f2, "%s/pt_matrices/compute_matrices_python/matter_multipoles/M12oneline_N256-matter_multipoles-vd4_f2.txt.dat", __CLASSDIR__);
            /**/
            /*-----------------------------------------------*/
            //GC: ORTHOGONAL -- start
            /**/
            sprintf(file2openM12_matter_multipoles_vv0_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N256-matter_multipoles-vv0_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv0_f3_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N256-matter_multipoles-vv0_f3.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd0_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N256-matter_multipoles-vd0_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd0_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N256-matter_multipoles-vd0_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_dd0_f0_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N256-matter_multipoles-dd0_f0.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_dd0_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N256-matter_multipoles-dd0_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv2_f3_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N256-matter_multipoles-vv2_f3.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd2_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N256-matter_multipoles-vd2_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vv4_f3_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N256-matter_multipoles-vv4_f3.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_multipoles_vd4_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_multipoles/M12oneline_N256-matter_multipoles-vd4_f2.txt-orthogonal.dat", __CLASSDIR__);
            /**/
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
            sprintf(file2openM12_matter_mu_powers_vd2_f1, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N256-matter_mu_powers-vd2_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vd2_f2, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N256-matter_mu_powers-vd2_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_dd2_f1, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N256-matter_mu_powers-dd2_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vv4_f2, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N256-matter_mu_powers-vv4_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vd4_f2, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N256-matter_mu_powers-vd4_f2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vv6_f3, "%s/pt_matrices/compute_matrices_python/matter_mu_powers/M12oneline_N256-matter_mu_powers-vv6_f3.txt.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- start
            sprintf(file2openM12_matter_mu_powers_vd2_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N256-matter_mu_powers-vd2_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vd2_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N256-matter_mu_powers-vd2_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_dd2_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N256-matter_mu_powers-dd2_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vv4_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N256-matter_mu_powers-vv4_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vd4_f2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N256-matter_mu_powers-vd4_f2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_matter_mu_powers_vv6_f3_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/matter_mu_powers/M12oneline_N256-matter_mu_powers-vv6_f3.txt-orthogonal.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
            /**/
            sprintf(file2openM12_bias_real_space_b2, "%s/pt_matrices/compute_matrices_python/bias_real_space/M12oneline_N256-bias_real_space-b2.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_bias_real_space_bG2, "%s/pt_matrices/compute_matrices_python/bias_real_space/M12oneline_N256-bias_real_space-G2.txt.dat", __CLASSDIR__);
            /**/
            /*-----------------------------------------------*/
            //GC: ORTHOGONAL -- start
            /**/
            sprintf(file2openM12_bias_real_space_b2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/bias_real_space/M12oneline_N256-bias_real_space-b2.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_bias_real_space_bG2_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/bias_real_space/M12oneline_N256-bias_real_space-G2.txt-orthogonal.dat", __CLASSDIR__);
            /**/
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
            sprintf(file2openM12_bias_multipoles_b2_vv0_f1, "%s/pt_matrices/compute_matrices_python/bias_multipoles/M12oneline_N256-bias_multipoles-b2_vv0_f1.txt.dat", __CLASSDIR__);
            sprintf(file2openM12_bias_multipoles_bG2_vv0_f1, "%s/pt_matrices/compute_matrices_python/bias_multipoles/M12oneline_N256-bias_multipoles-bG2_vv0_f1.txt.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- start
            sprintf(file2openM12_bias_multipoles_b2_vv0_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/bias_multipoles/M12oneline_N256-bias_multipoles-b2_vv0_f1.txt-orthogonal.dat", __CLASSDIR__);
            sprintf(file2openM12_bias_multipoles_bG2_vv0_f1_ortho, "%s/ORTHOGONAL_NG_matrices/computation_of_matrices/bias_multipoles/M12oneline_N256-bias_multipoles-bG2_vv0_f1.txt-orthogonal.dat", __CLASSDIR__);
            //GC: ORTHOGONAL -- finish
            /*-----------------------------------------------*/
                
            }
        }

        int index_M22 = 0;

        class_alloc(pnlpt->M22_oneline, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile22;

        myFile22 = fopen(file2openM22, "r");

        if (myFile22 == NULL)
        {
            printf("Error Reading File M22oneline_....dat\n");
            exit(0);
        }

        for (index_M22 = 0; index_M22 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M22++)
        {
            fscanf(myFile22, "%lf", &pnlpt->M22_oneline[index_M22]);
        }

        fclose(myFile22);

        int index_M13 = 0;

        class_alloc(pnlpt->M13_oneline, ((ppr->nmax_nlpt + 1) * 2) * sizeof(double), pnlpt->error_message);

        FILE *myFile13;

        myFile13 = fopen(file2openM13, "r");

        if (myFile13 == NULL)
        {
            printf("Error Reading File M13oneline_....dat\n");
            exit(0);
        }

        for (index_M13 = 0; index_M13 < 2 * (ppr->nmax_nlpt + 1); index_M13++)
        {
            fscanf(myFile13, "%lf", &pnlpt->M13_oneline[index_M13]);
        }

        fclose(myFile13);

        class_alloc(pnlpt->M22basic_oneline, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile22_basic;

        int index_M22_basic = 0;

        myFile22_basic = fopen(file2openM22basic, "r");

        if (myFile22_basic == NULL)
        {
            printf("Error Reading File M22basiconeline_....dat\n");
            exit(0);
        }

        for (index_M22_basic = 0; index_M22_basic < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M22_basic++)
        {
            fscanf(myFile22_basic, "%lf", &pnlpt->M22basic_oneline[index_M22_basic]);
        }

        fclose(myFile22_basic);

        int index_MIFG2 = 0;

        class_alloc(pnlpt->IFG2_oneline, ((ppr->nmax_nlpt + 1) * 2) * sizeof(double), pnlpt->error_message);

        FILE *myFile_IFG2;

        myFile_IFG2 = fopen(file2openM13basic, "r");

        if (myFile_IFG2 == NULL)
        {
            printf("Error Reading File IFG2oneline_....dat\n");
            exit(0);
        }

        for (index_MIFG2 = 0; index_MIFG2 < 2 * (ppr->nmax_nlpt + 1); index_MIFG2++)
        {
            fscanf(myFile_IFG2, "%lf", &pnlpt->IFG2_oneline[index_MIFG2]);
        }

        fclose(myFile_IFG2);
        
        
        if (SWITCH_index == 1)
            
            //GC - SWITCH -> now the issue is that if I do like this I also do not allocate these matrices... E.g. I do not allocate M12_oneline and so on. This is bad because later I deallocate it... And these objects are in the .h -> the quickest way is to do an else if, and in the else just call the class_alloc... The not opening and scanning the files is the bulk of the issues...
            
        {

        //GC!

        /*-----------------------------------------------*/

        int index_M12 = 0;

        class_alloc(pnlpt->M12_oneline, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12;

        myFile12 = fopen(file2openM12_matter, "r");

        if (myFile12 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12 = 0; index_M12 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12++)
        {
            fscanf(myFile12, "%lf", &pnlpt->M12_oneline[index_M12]);
        }

        fclose(myFile12);

        //GC: ORTHOGONAL -- start

        int index_M12_ortho = 0;

        class_alloc(pnlpt->M12_oneline_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_ortho;

        myFile12_ortho = fopen(file2openM12_matter_ortho, "r");

        if (myFile12_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_ortho = 0; index_M12_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_ortho++)
        {
            fscanf(myFile12_ortho, "%lf", &pnlpt->M12_oneline_ortho[index_M12_ortho]);
        }

        fclose(myFile12_ortho);

        //GC: ORTHOGONAL -- finish

        /*-----------------------------------------------*/

        /**/

        int index_M12_matter_multipoles_vv0_f2 = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv0_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vv0_f2;

        myFile12_matter_multipoles_vv0_f2 = fopen(file2openM12_matter_multipoles_vv0_f2, "r");

        if (myFile12_matter_multipoles_vv0_f2 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vv0_f2 = 0; index_M12_matter_multipoles_vv0_f2 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vv0_f2++)
        {
            fscanf(myFile12_matter_multipoles_vv0_f2, "%lf", &pnlpt->M12_oneline_matter_multipoles_vv0_f2[index_M12_matter_multipoles_vv0_f2]);
        }

        fclose(myFile12_matter_multipoles_vv0_f2);

        int index_M12_matter_multipoles_vv0_f3 = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv0_f3, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vv0_f3;

        myFile12_matter_multipoles_vv0_f3 = fopen(file2openM12_matter_multipoles_vv0_f3, "r");

        if (myFile12_matter_multipoles_vv0_f3 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vv0_f3 = 0; index_M12_matter_multipoles_vv0_f3 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vv0_f3++)
        {
            fscanf(myFile12_matter_multipoles_vv0_f3, "%lf", &pnlpt->M12_oneline_matter_multipoles_vv0_f3[index_M12_matter_multipoles_vv0_f3]);
        }

        fclose(myFile12_matter_multipoles_vv0_f3);

        int index_M12_matter_multipoles_vd0_f1 = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd0_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vd0_f1;

        myFile12_matter_multipoles_vd0_f1 = fopen(file2openM12_matter_multipoles_vd0_f1, "r");

        if (myFile12_matter_multipoles_vd0_f1 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vd0_f1 = 0; index_M12_matter_multipoles_vd0_f1 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vd0_f1++)
        {
            fscanf(myFile12_matter_multipoles_vd0_f1, "%lf", &pnlpt->M12_oneline_matter_multipoles_vd0_f1[index_M12_matter_multipoles_vd0_f1]);
        }

        fclose(myFile12_matter_multipoles_vd0_f1);

        int index_M12_matter_multipoles_vd0_f2 = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd0_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vd0_f2;

        myFile12_matter_multipoles_vd0_f2 = fopen(file2openM12_matter_multipoles_vd0_f2, "r");

        if (myFile12_matter_multipoles_vd0_f2 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vd0_f2 = 0; index_M12_matter_multipoles_vd0_f2 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vd0_f2++)
        {
            fscanf(myFile12_matter_multipoles_vd0_f2, "%lf", &pnlpt->M12_oneline_matter_multipoles_vd0_f2[index_M12_matter_multipoles_vd0_f2]);
        }

        fclose(myFile12_matter_multipoles_vd0_f2);

        int index_M12_matter_multipoles_dd0_f0 = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_dd0_f0, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_dd0_f0;

        myFile12_matter_multipoles_dd0_f0 = fopen(file2openM12_matter_multipoles_dd0_f0, "r");

        if (myFile12_matter_multipoles_dd0_f0 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_dd0_f0 = 0; index_M12_matter_multipoles_dd0_f0 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_dd0_f0++)
        {
            fscanf(myFile12_matter_multipoles_dd0_f0, "%lf", &pnlpt->M12_oneline_matter_multipoles_dd0_f0[index_M12_matter_multipoles_dd0_f0]);
        }

        fclose(myFile12_matter_multipoles_dd0_f0);

        int index_M12_matter_multipoles_dd0_f1 = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_dd0_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_dd0_f1;

        myFile12_matter_multipoles_dd0_f1 = fopen(file2openM12_matter_multipoles_dd0_f1, "r");

        if (myFile12_matter_multipoles_dd0_f1 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_dd0_f1 = 0; index_M12_matter_multipoles_dd0_f1 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_dd0_f1++)
        {
            fscanf(myFile12_matter_multipoles_dd0_f1, "%lf", &pnlpt->M12_oneline_matter_multipoles_dd0_f1[index_M12_matter_multipoles_dd0_f1]);
        }

        fclose(myFile12_matter_multipoles_dd0_f1);

        int index_M12_matter_multipoles_vv2_f3 = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv2_f3, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vv2_f3;

        myFile12_matter_multipoles_vv2_f3 = fopen(file2openM12_matter_multipoles_vv2_f3, "r");

        if (myFile12_matter_multipoles_vv2_f3 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vv2_f3 = 0; index_M12_matter_multipoles_vv2_f3 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vv2_f3++)
        {
            fscanf(myFile12_matter_multipoles_vv2_f3, "%lf", &pnlpt->M12_oneline_matter_multipoles_vv2_f3[index_M12_matter_multipoles_vv2_f3]);
        }

        fclose(myFile12_matter_multipoles_vv2_f3);

        int index_M12_matter_multipoles_vd2_f2 = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd2_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vd2_f2;

        myFile12_matter_multipoles_vd2_f2 = fopen(file2openM12_matter_multipoles_vd2_f2, "r");

        if (myFile12_matter_multipoles_vd2_f2 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vd2_f2 = 0; index_M12_matter_multipoles_vd2_f2 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vd2_f2++)
        {
            fscanf(myFile12_matter_multipoles_vd2_f2, "%lf", &pnlpt->M12_oneline_matter_multipoles_vd2_f2[index_M12_matter_multipoles_vd2_f2]);
        }

        fclose(myFile12_matter_multipoles_vd2_f2);

        int index_M12_matter_multipoles_vv4_f3 = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv4_f3, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vv4_f3;

        myFile12_matter_multipoles_vv4_f3 = fopen(file2openM12_matter_multipoles_vv4_f3, "r");

        if (myFile12_matter_multipoles_vv4_f3 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vv4_f3 = 0; index_M12_matter_multipoles_vv4_f3 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vv4_f3++)
        {
            fscanf(myFile12_matter_multipoles_vv4_f3, "%lf", &pnlpt->M12_oneline_matter_multipoles_vv4_f3[index_M12_matter_multipoles_vv4_f3]);
        }

        fclose(myFile12_matter_multipoles_vv4_f3);

        int index_M12_matter_multipoles_vd4_f2 = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd4_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vd4_f2;

        myFile12_matter_multipoles_vd4_f2 = fopen(file2openM12_matter_multipoles_vd4_f2, "r");

        if (myFile12_matter_multipoles_vd4_f2 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vd4_f2 = 0; index_M12_matter_multipoles_vd4_f2 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vd4_f2++)
        {
            fscanf(myFile12_matter_multipoles_vd4_f2, "%lf", &pnlpt->M12_oneline_matter_multipoles_vd4_f2[index_M12_matter_multipoles_vd4_f2]);
        }

        fclose(myFile12_matter_multipoles_vd4_f2);

        /**/

        //GC: ORTHOGONAL -- start

        int index_M12_matter_multipoles_vv0_f2_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv0_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vv0_f2_ortho;

        myFile12_matter_multipoles_vv0_f2_ortho = fopen(file2openM12_matter_multipoles_vv0_f2_ortho, "r");

        if (myFile12_matter_multipoles_vv0_f2_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vv0_f2_ortho = 0; index_M12_matter_multipoles_vv0_f2_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vv0_f2_ortho++)
        {
            fscanf(myFile12_matter_multipoles_vv0_f2_ortho, "%lf", &pnlpt->M12_oneline_matter_multipoles_vv0_f2_ortho[index_M12_matter_multipoles_vv0_f2_ortho]);
        }

        fclose(myFile12_matter_multipoles_vv0_f2_ortho);

        int index_M12_matter_multipoles_vv0_f3_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv0_f3_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vv0_f3_ortho;

        myFile12_matter_multipoles_vv0_f3_ortho = fopen(file2openM12_matter_multipoles_vv0_f3_ortho, "r");

        if (myFile12_matter_multipoles_vv0_f3_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vv0_f3_ortho = 0; index_M12_matter_multipoles_vv0_f3_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vv0_f3_ortho++)
        {
            fscanf(myFile12_matter_multipoles_vv0_f3_ortho, "%lf", &pnlpt->M12_oneline_matter_multipoles_vv0_f3_ortho[index_M12_matter_multipoles_vv0_f3_ortho]);
        }

        fclose(myFile12_matter_multipoles_vv0_f3_ortho);

        int index_M12_matter_multipoles_vd0_f1_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd0_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vd0_f1_ortho;

        myFile12_matter_multipoles_vd0_f1_ortho = fopen(file2openM12_matter_multipoles_vd0_f1_ortho, "r");

        if (myFile12_matter_multipoles_vd0_f1_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vd0_f1_ortho = 0; index_M12_matter_multipoles_vd0_f1_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vd0_f1_ortho++)
        {
            fscanf(myFile12_matter_multipoles_vd0_f1_ortho, "%lf", &pnlpt->M12_oneline_matter_multipoles_vd0_f1_ortho[index_M12_matter_multipoles_vd0_f1_ortho]);
        }

        fclose(myFile12_matter_multipoles_vd0_f1_ortho);

        int index_M12_matter_multipoles_vd0_f2_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd0_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vd0_f2_ortho;

        myFile12_matter_multipoles_vd0_f2_ortho = fopen(file2openM12_matter_multipoles_vd0_f2_ortho, "r");

        if (myFile12_matter_multipoles_vd0_f2_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vd0_f2_ortho = 0; index_M12_matter_multipoles_vd0_f2_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vd0_f2_ortho++)
        {
            fscanf(myFile12_matter_multipoles_vd0_f2_ortho, "%lf", &pnlpt->M12_oneline_matter_multipoles_vd0_f2_ortho[index_M12_matter_multipoles_vd0_f2_ortho]);
        }

        fclose(myFile12_matter_multipoles_vd0_f2_ortho);

        int index_M12_matter_multipoles_dd0_f0_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_dd0_f0_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_dd0_f0_ortho;

        myFile12_matter_multipoles_dd0_f0_ortho = fopen(file2openM12_matter_multipoles_dd0_f0_ortho, "r");

        if (myFile12_matter_multipoles_dd0_f0_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_dd0_f0_ortho = 0; index_M12_matter_multipoles_dd0_f0_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_dd0_f0_ortho++)
        {
            fscanf(myFile12_matter_multipoles_dd0_f0_ortho, "%lf", &pnlpt->M12_oneline_matter_multipoles_dd0_f0_ortho[index_M12_matter_multipoles_dd0_f0_ortho]);
        }

        fclose(myFile12_matter_multipoles_dd0_f0_ortho);

        int index_M12_matter_multipoles_dd0_f1_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_dd0_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_dd0_f1_ortho;

        myFile12_matter_multipoles_dd0_f1_ortho = fopen(file2openM12_matter_multipoles_dd0_f1_ortho, "r");

        if (myFile12_matter_multipoles_dd0_f1_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_dd0_f1_ortho = 0; index_M12_matter_multipoles_dd0_f1_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_dd0_f1_ortho++)
        {
            fscanf(myFile12_matter_multipoles_dd0_f1_ortho, "%lf", &pnlpt->M12_oneline_matter_multipoles_dd0_f1_ortho[index_M12_matter_multipoles_dd0_f1_ortho]);
        }

        fclose(myFile12_matter_multipoles_dd0_f1_ortho);

        int index_M12_matter_multipoles_vv2_f3_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv2_f3_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vv2_f3_ortho;

        myFile12_matter_multipoles_vv2_f3_ortho = fopen(file2openM12_matter_multipoles_vv2_f3_ortho, "r");

        if (myFile12_matter_multipoles_vv2_f3_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vv2_f3_ortho = 0; index_M12_matter_multipoles_vv2_f3_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vv2_f3_ortho++)
        {
            fscanf(myFile12_matter_multipoles_vv2_f3_ortho, "%lf", &pnlpt->M12_oneline_matter_multipoles_vv2_f3_ortho[index_M12_matter_multipoles_vv2_f3_ortho]);
        }

        fclose(myFile12_matter_multipoles_vv2_f3_ortho);

        int index_M12_matter_multipoles_vd2_f2_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd2_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vd2_f2_ortho;

        myFile12_matter_multipoles_vd2_f2_ortho = fopen(file2openM12_matter_multipoles_vd2_f2_ortho, "r");

        if (myFile12_matter_multipoles_vd2_f2_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vd2_f2_ortho = 0; index_M12_matter_multipoles_vd2_f2_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vd2_f2_ortho++)
        {
            fscanf(myFile12_matter_multipoles_vd2_f2_ortho, "%lf", &pnlpt->M12_oneline_matter_multipoles_vd2_f2_ortho[index_M12_matter_multipoles_vd2_f2_ortho]);
        }

        fclose(myFile12_matter_multipoles_vd2_f2_ortho);

        int index_M12_matter_multipoles_vv4_f3_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv4_f3_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vv4_f3_ortho;

        myFile12_matter_multipoles_vv4_f3_ortho = fopen(file2openM12_matter_multipoles_vv4_f3_ortho, "r");

        if (myFile12_matter_multipoles_vv4_f3_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vv4_f3_ortho = 0; index_M12_matter_multipoles_vv4_f3_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vv4_f3_ortho++)
        {
            fscanf(myFile12_matter_multipoles_vv4_f3_ortho, "%lf", &pnlpt->M12_oneline_matter_multipoles_vv4_f3_ortho[index_M12_matter_multipoles_vv4_f3_ortho]);
        }

        fclose(myFile12_matter_multipoles_vv4_f3_ortho);

        int index_M12_matter_multipoles_vd4_f2_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd4_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_multipoles_vd4_f2_ortho;

        myFile12_matter_multipoles_vd4_f2_ortho = fopen(file2openM12_matter_multipoles_vd4_f2_ortho, "r");

        if (myFile12_matter_multipoles_vd4_f2_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_multipoles_vd4_f2_ortho = 0; index_M12_matter_multipoles_vd4_f2_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_multipoles_vd4_f2_ortho++)
        {
            fscanf(myFile12_matter_multipoles_vd4_f2_ortho, "%lf", &pnlpt->M12_oneline_matter_multipoles_vd4_f2_ortho[index_M12_matter_multipoles_vd4_f2_ortho]);
        }

        fclose(myFile12_matter_multipoles_vd4_f2_ortho);

        //GC: ORTHOGONAL -- finish

        /*-----------------------------------------------*/

        int index_M12_matter_mu_powers_vd2_f1 = 0;

        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vd2_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_mu_powers_vd2_f1;

        myFile12_matter_mu_powers_vd2_f1 = fopen(file2openM12_matter_mu_powers_vd2_f1, "r");

        if (myFile12_matter_mu_powers_vd2_f1 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_mu_powers_vd2_f1 = 0; index_M12_matter_mu_powers_vd2_f1 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_mu_powers_vd2_f1++)
        {
            fscanf(myFile12_matter_mu_powers_vd2_f1, "%lf", &pnlpt->M12_oneline_matter_mu_powers_vd2_f1[index_M12_matter_mu_powers_vd2_f1]);
        }

        fclose(myFile12_matter_mu_powers_vd2_f1);

        //GC: ORTHOGONAL -- start

        int index_M12_matter_mu_powers_vd2_f1_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vd2_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_mu_powers_vd2_f1_ortho;

        myFile12_matter_mu_powers_vd2_f1_ortho = fopen(file2openM12_matter_mu_powers_vd2_f1_ortho, "r");

        if (myFile12_matter_mu_powers_vd2_f1_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_mu_powers_vd2_f1_ortho = 0; index_M12_matter_mu_powers_vd2_f1_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_mu_powers_vd2_f1_ortho++)
        {
            fscanf(myFile12_matter_mu_powers_vd2_f1_ortho, "%lf", &pnlpt->M12_oneline_matter_mu_powers_vd2_f1_ortho[index_M12_matter_mu_powers_vd2_f1_ortho]);
        }

        fclose(myFile12_matter_mu_powers_vd2_f1_ortho);

        //GC: ORTHOGONAL -- finish

        int index_M12_matter_mu_powers_vd2_f2 = 0;

        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vd2_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_mu_powers_vd2_f2;

        myFile12_matter_mu_powers_vd2_f2 = fopen(file2openM12_matter_mu_powers_vd2_f2, "r");

        if (myFile12_matter_mu_powers_vd2_f2 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_mu_powers_vd2_f2 = 0; index_M12_matter_mu_powers_vd2_f2 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_mu_powers_vd2_f2++)
        {
            fscanf(myFile12_matter_mu_powers_vd2_f2, "%lf", &pnlpt->M12_oneline_matter_mu_powers_vd2_f2[index_M12_matter_mu_powers_vd2_f2]);
        }

        fclose(myFile12_matter_mu_powers_vd2_f2);

        //GC: ORTHOGONAL -- start

        int index_M12_matter_mu_powers_vd2_f2_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vd2_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_mu_powers_vd2_f2_ortho;

        myFile12_matter_mu_powers_vd2_f2_ortho = fopen(file2openM12_matter_mu_powers_vd2_f2_ortho, "r");

        if (myFile12_matter_mu_powers_vd2_f2_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_mu_powers_vd2_f2_ortho = 0; index_M12_matter_mu_powers_vd2_f2_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_mu_powers_vd2_f2_ortho++)
        {
            fscanf(myFile12_matter_mu_powers_vd2_f2_ortho, "%lf", &pnlpt->M12_oneline_matter_mu_powers_vd2_f2_ortho[index_M12_matter_mu_powers_vd2_f2_ortho]);
        }

        fclose(myFile12_matter_mu_powers_vd2_f2_ortho);

        //GC: ORTHOGONAL -- finish

        int index_M12_matter_mu_powers_dd2_f1 = 0;

        class_alloc(pnlpt->M12_oneline_matter_mu_powers_dd2_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_mu_powers_dd2_f1;

        myFile12_matter_mu_powers_dd2_f1 = fopen(file2openM12_matter_mu_powers_dd2_f1, "r");

        if (myFile12_matter_mu_powers_dd2_f1 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_mu_powers_dd2_f1 = 0; index_M12_matter_mu_powers_dd2_f1 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_mu_powers_dd2_f1++)
        {
            fscanf(myFile12_matter_mu_powers_dd2_f1, "%lf", &pnlpt->M12_oneline_matter_mu_powers_dd2_f1[index_M12_matter_mu_powers_dd2_f1]);
        }

        fclose(myFile12_matter_mu_powers_dd2_f1);

        //GC: ORTHOGONAL -- start

        int index_M12_matter_mu_powers_dd2_f1_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_mu_powers_dd2_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_mu_powers_dd2_f1_ortho;

        myFile12_matter_mu_powers_dd2_f1_ortho = fopen(file2openM12_matter_mu_powers_dd2_f1_ortho, "r");

        if (myFile12_matter_mu_powers_dd2_f1_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_mu_powers_dd2_f1_ortho = 0; index_M12_matter_mu_powers_dd2_f1_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_mu_powers_dd2_f1_ortho++)
        {
            fscanf(myFile12_matter_mu_powers_dd2_f1_ortho, "%lf", &pnlpt->M12_oneline_matter_mu_powers_dd2_f1_ortho[index_M12_matter_mu_powers_dd2_f1_ortho]);
        }

        fclose(myFile12_matter_mu_powers_dd2_f1_ortho);

        //GC: ORTHOGONAL -- finish

        int index_M12_matter_mu_powers_vv4_f2 = 0;

        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vv4_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_mu_powers_vv4_f2;

        myFile12_matter_mu_powers_vv4_f2 = fopen(file2openM12_matter_mu_powers_vv4_f2, "r");

        if (myFile12_matter_mu_powers_vv4_f2 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_mu_powers_vv4_f2 = 0; index_M12_matter_mu_powers_vv4_f2 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_mu_powers_vv4_f2++)
        {
            fscanf(myFile12_matter_mu_powers_vv4_f2, "%lf", &pnlpt->M12_oneline_matter_mu_powers_vv4_f2[index_M12_matter_mu_powers_vv4_f2]);
        }

        fclose(myFile12_matter_mu_powers_vv4_f2);

        //GC: ORTHOGONAL -- start

        int index_M12_matter_mu_powers_vv4_f2_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vv4_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_mu_powers_vv4_f2_ortho;

        myFile12_matter_mu_powers_vv4_f2_ortho = fopen(file2openM12_matter_mu_powers_vv4_f2_ortho, "r");

        if (myFile12_matter_mu_powers_vv4_f2_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_mu_powers_vv4_f2_ortho = 0; index_M12_matter_mu_powers_vv4_f2_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_mu_powers_vv4_f2_ortho++)
        {
            fscanf(myFile12_matter_mu_powers_vv4_f2_ortho, "%lf", &pnlpt->M12_oneline_matter_mu_powers_vv4_f2_ortho[index_M12_matter_mu_powers_vv4_f2_ortho]);
        }

        fclose(myFile12_matter_mu_powers_vv4_f2_ortho);

        //GC: ORTHOGONAL -- finish

        int index_M12_matter_mu_powers_vd4_f2 = 0;

        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vd4_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_mu_powers_vd4_f2;

        myFile12_matter_mu_powers_vd4_f2 = fopen(file2openM12_matter_mu_powers_vd4_f2, "r");

        if (myFile12_matter_mu_powers_vd4_f2 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_mu_powers_vd4_f2 = 0; index_M12_matter_mu_powers_vd4_f2 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_mu_powers_vd4_f2++)
        {
            fscanf(myFile12_matter_mu_powers_vd4_f2, "%lf", &pnlpt->M12_oneline_matter_mu_powers_vd4_f2[index_M12_matter_mu_powers_vd4_f2]);
        }

        fclose(myFile12_matter_mu_powers_vd4_f2);

        //GC: ORTHOGONAL -- start

        int index_M12_matter_mu_powers_vd4_f2_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vd4_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_mu_powers_vd4_f2_ortho;

        myFile12_matter_mu_powers_vd4_f2_ortho = fopen(file2openM12_matter_mu_powers_vd4_f2_ortho, "r");

        if (myFile12_matter_mu_powers_vd4_f2_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_mu_powers_vd4_f2_ortho = 0; index_M12_matter_mu_powers_vd4_f2_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_mu_powers_vd4_f2_ortho++)
        {
            fscanf(myFile12_matter_mu_powers_vd4_f2_ortho, "%lf", &pnlpt->M12_oneline_matter_mu_powers_vd4_f2_ortho[index_M12_matter_mu_powers_vd4_f2_ortho]);
        }

        fclose(myFile12_matter_mu_powers_vd4_f2_ortho);

        //GC: ORTHOGONAL -- finish

        int index_M12_matter_mu_powers_vv6_f3 = 0;

        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vv6_f3, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_mu_powers_vv6_f3;

        myFile12_matter_mu_powers_vv6_f3 = fopen(file2openM12_matter_mu_powers_vv6_f3, "r");

        if (myFile12_matter_mu_powers_vv6_f3 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_mu_powers_vv6_f3 = 0; index_M12_matter_mu_powers_vv6_f3 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_mu_powers_vv6_f3++)
        {
            fscanf(myFile12_matter_mu_powers_vv6_f3, "%lf", &pnlpt->M12_oneline_matter_mu_powers_vv6_f3[index_M12_matter_mu_powers_vv6_f3]);
        }

        fclose(myFile12_matter_mu_powers_vv6_f3);

        //GC: ORTHOGONAL -- start

        int index_M12_matter_mu_powers_vv6_f3_ortho = 0;

        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vv6_f3_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_matter_mu_powers_vv6_f3_ortho;

        myFile12_matter_mu_powers_vv6_f3_ortho = fopen(file2openM12_matter_mu_powers_vv6_f3_ortho, "r");

        if (myFile12_matter_mu_powers_vv6_f3_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_matter_mu_powers_vv6_f3_ortho = 0; index_M12_matter_mu_powers_vv6_f3_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_matter_mu_powers_vv6_f3_ortho++)
        {
            fscanf(myFile12_matter_mu_powers_vv6_f3_ortho, "%lf", &pnlpt->M12_oneline_matter_mu_powers_vv6_f3_ortho[index_M12_matter_mu_powers_vv6_f3_ortho]);
        }

        fclose(myFile12_matter_mu_powers_vv6_f3_ortho);

        //GC: ORTHOGONAL -- finish

        /*-----------------------------------------------*/

        /**/

        int index_M12_bias_real_space_b2 = 0;

        class_alloc(pnlpt->M12_oneline_bias_real_space_b2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_bias_real_space_b2;

        myFile12_bias_real_space_b2 = fopen(file2openM12_bias_real_space_b2, "r");

        if (myFile12_bias_real_space_b2 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_bias_real_space_b2 = 0; index_M12_bias_real_space_b2 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_bias_real_space_b2++)
        {
            fscanf(myFile12_bias_real_space_b2, "%lf", &pnlpt->M12_oneline_bias_real_space_b2[index_M12_bias_real_space_b2]);
        }

        fclose(myFile12_bias_real_space_b2);

        int index_M12_bias_real_space_bG2 = 0;

        class_alloc(pnlpt->M12_oneline_bias_real_space_bG2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_bias_real_space_bG2;

        myFile12_bias_real_space_bG2 = fopen(file2openM12_bias_real_space_bG2, "r");

        if (myFile12_bias_real_space_bG2 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_bias_real_space_bG2 = 0; index_M12_bias_real_space_bG2 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_bias_real_space_bG2++)
        {
            fscanf(myFile12_bias_real_space_bG2, "%lf", &pnlpt->M12_oneline_bias_real_space_bG2[index_M12_bias_real_space_bG2]);
        }

        fclose(myFile12_bias_real_space_bG2);

        /**/

        //GC: ORTHOGONAL -- start

        int index_M12_bias_real_space_b2_ortho = 0;

        class_alloc(pnlpt->M12_oneline_bias_real_space_b2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_bias_real_space_b2_ortho;

        myFile12_bias_real_space_b2_ortho = fopen(file2openM12_bias_real_space_b2_ortho, "r");

        if (myFile12_bias_real_space_b2_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_bias_real_space_b2_ortho = 0; index_M12_bias_real_space_b2_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_bias_real_space_b2_ortho++)
        {
            fscanf(myFile12_bias_real_space_b2_ortho, "%lf", &pnlpt->M12_oneline_bias_real_space_b2_ortho[index_M12_bias_real_space_b2_ortho]);
        }

        fclose(myFile12_bias_real_space_b2_ortho);

        int index_M12_bias_real_space_bG2_ortho = 0;

        class_alloc(pnlpt->M12_oneline_bias_real_space_bG2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_bias_real_space_bG2_ortho;

        myFile12_bias_real_space_bG2_ortho = fopen(file2openM12_bias_real_space_bG2_ortho, "r");

        if (myFile12_bias_real_space_bG2_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_bias_real_space_bG2_ortho = 0; index_M12_bias_real_space_bG2_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_bias_real_space_bG2_ortho++)
        {
            fscanf(myFile12_bias_real_space_bG2_ortho, "%lf", &pnlpt->M12_oneline_bias_real_space_bG2_ortho[index_M12_bias_real_space_bG2_ortho]);
        }

        fclose(myFile12_bias_real_space_bG2_ortho);

        //GC: ORTHOGONAL -- finish

        /*-----------------------------------------------*/

        int index_M12_bias_multipoles_b2_vv0_f1 = 0;

        class_alloc(pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_bias_multipoles_b2_vv0_f1;

        myFile12_bias_multipoles_b2_vv0_f1 = fopen(file2openM12_bias_multipoles_b2_vv0_f1, "r");

        if (myFile12_bias_multipoles_b2_vv0_f1 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_bias_multipoles_b2_vv0_f1 = 0; index_M12_bias_multipoles_b2_vv0_f1 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_bias_multipoles_b2_vv0_f1++)
        {
            fscanf(myFile12_bias_multipoles_b2_vv0_f1, "%lf", &pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1[index_M12_bias_multipoles_b2_vv0_f1]);
        }

        fclose(myFile12_bias_multipoles_b2_vv0_f1);

        //GC: ORTHOGONAL -- start

        int index_M12_bias_multipoles_b2_vv0_f1_ortho = 0;

        class_alloc(pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_bias_multipoles_b2_vv0_f1_ortho;

        myFile12_bias_multipoles_b2_vv0_f1_ortho = fopen(file2openM12_bias_multipoles_b2_vv0_f1_ortho, "r");

        if (myFile12_bias_multipoles_b2_vv0_f1_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_bias_multipoles_b2_vv0_f1_ortho = 0; index_M12_bias_multipoles_b2_vv0_f1_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_bias_multipoles_b2_vv0_f1_ortho++)
        {
            fscanf(myFile12_bias_multipoles_b2_vv0_f1_ortho, "%lf", &pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1_ortho[index_M12_bias_multipoles_b2_vv0_f1_ortho]);
        }

        fclose(myFile12_bias_multipoles_b2_vv0_f1_ortho);

        //GC: ORTHOGONAL -- finish

        int index_M12_bias_multipoles_bG2_vv0_f1 = 0;

        class_alloc(pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_bias_multipoles_bG2_vv0_f1;

        myFile12_bias_multipoles_bG2_vv0_f1 = fopen(file2openM12_bias_multipoles_bG2_vv0_f1, "r");

        if (myFile12_bias_multipoles_bG2_vv0_f1 == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_bias_multipoles_bG2_vv0_f1 = 0; index_M12_bias_multipoles_bG2_vv0_f1 < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_bias_multipoles_bG2_vv0_f1++)
        {
            fscanf(myFile12_bias_multipoles_bG2_vv0_f1, "%lf", &pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1[index_M12_bias_multipoles_bG2_vv0_f1]);
        }

        fclose(myFile12_bias_multipoles_bG2_vv0_f1);

        //GC: ORTHOGONAL -- start

        int index_M12_bias_multipoles_bG2_vv0_f1_ortho = 0;

        class_alloc(pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        FILE *myFile12_bias_multipoles_bG2_vv0_f1_ortho;

        myFile12_bias_multipoles_bG2_vv0_f1_ortho = fopen(file2openM12_bias_multipoles_bG2_vv0_f1_ortho, "r");

        if (myFile12_bias_multipoles_bG2_vv0_f1_ortho == NULL)
        {
            printf("Error Reading File M12oneline_....dat\n");
            exit(0);
        }

        for (index_M12_bias_multipoles_bG2_vv0_f1_ortho = 0; index_M12_bias_multipoles_bG2_vv0_f1_ortho < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2); index_M12_bias_multipoles_bG2_vv0_f1_ortho++)
        {
            fscanf(myFile12_bias_multipoles_bG2_vv0_f1_ortho, "%lf", &pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1_ortho[index_M12_bias_multipoles_bG2_vv0_f1_ortho]);
        }

        fclose(myFile12_bias_multipoles_bG2_vv0_f1_ortho);

        //GC: ORTHOGONAL -- finish

        /*-----------------------------------------------*/

        //int endmatrix=clock(); //GC!

        //printf("Reading the matrices takes %d musec.\n",endmatrix-startmatrix); //GC!

        //GC!

        //GC -> before I had it in the middle of the Gauss stuff... Now it is in a position that makes sense...
            
        }
        
        else //GC - SWITCH -> CHECK SYNTAX...
            
            
        {
            
          
        
        class_alloc(pnlpt->M12_oneline, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv0_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv0_f3, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd0_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd0_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_dd0_f0, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_dd0_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv2_f3, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd2_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv4_f3, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd4_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv0_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv0_f3_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd0_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd0_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_dd0_f0_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_dd0_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv2_f3_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd2_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vv4_f3_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_multipoles_vd4_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vd2_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vd2_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vd2_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vd2_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_mu_powers_dd2_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_mu_powers_dd2_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vv4_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vv4_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vd4_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vd4_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vv6_f3, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_matter_mu_powers_vv6_f3_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);

        class_alloc(pnlpt->M12_oneline_bias_real_space_b2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_bias_real_space_bG2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_bias_real_space_b2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_bias_real_space_bG2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)) * sizeof(double), pnlpt->error_message);
            
            
        }
            
            

        char file2openGauss[256];
        sprintf(file2openGauss, "%s/pt_matrices/gauss_tab.dat", __CLASSDIR__);
        int index_gauss = 0;

        class_alloc(pnlpt->gauss_w, 40 * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->gauss_x, 40 * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->gauss, 80 * sizeof(double), pnlpt->error_message);

        FILE *Gauss_file;

        Gauss_file = fopen(file2openGauss, "r");

        if (Gauss_file == NULL)
        {
            printf("Error Reading File gauss_tab.dat\n");
            exit(0);
        }

        for (index_gauss = 0; index_gauss < 80; index_gauss++)
        {
            fscanf(Gauss_file, "%lf", &pnlpt->gauss[index_gauss]);
        }

        fclose(Gauss_file);

        for (index_gauss = 0; index_gauss < 40; index_gauss++)
        {
            pnlpt->gauss_x[index_gauss] = pnlpt->gauss[index_gauss];
            pnlpt->gauss_w[index_gauss] = pnlpt->gauss[40 + index_gauss];
            //  printf("%lf %lf\n",gauss_x[index_gauss],gauss_w[index_gauss]);
        }

        //GC -> this is a place for future optimization! Also Misha wrote this indeed, before my modifications... I keep allocating, since the bulk is importing the files and actually computing stuff...

        class_alloc(pnlpt->M22_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M13_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22basic_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_0_b1b2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_0_b2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_0_b1bG2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_0_bG2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        class_alloc(pnlpt->M22_2_b1b2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_2_b2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_2_b1bG2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_2_bG2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        class_alloc(pnlpt->M22_4_b2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_4_bG2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        class_alloc(pnlpt->IFG2_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);

        class_alloc(pnlpt->M_IG2G2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M_Id2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M_IG2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M_Id2G2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(double complex), pnlpt->error_message);
        class_alloc(pnlpt->M13_0_vv_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_0_vv_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M13_0_vd_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_0_vd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M13_0_dd_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_0_dd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        class_alloc(pnlpt->M13_2_vv_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_2_vv_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        class_alloc(pnlpt->M13_2_vd_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_2_vd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        class_alloc(pnlpt->M13_2_dd_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_2_dd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        class_alloc(pnlpt->M13_4_vv_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_4_vv_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        class_alloc(pnlpt->M13_4_vd_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_4_vd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_4_dd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        class_alloc(pnlpt->M13_mu2_dd_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M13_mu2_vd_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M13_mu4_vv_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M13_mu4_vd_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M13_mu6_oneline_complex, (ppr->nmax_nlpt + 1) * sizeof(complex double), pnlpt->error_message);

        class_alloc(pnlpt->M22_oneline_mu2_vd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_mu2_dd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_mu4_vd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_mu4_vv_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_mu4_dd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_mu6_vv_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_mu6_vd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M22_oneline_mu8_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC -> my allocations start here...

        //GC -> what do I need? M12_oneline_complex, M12_oneline_0_vv_complex, M12_oneline_0_vd_complex, M12_oneline_0_dd_complex, M12_oneline_2_vv_complex, M12_oneline_2_vd_complex, M12_oneline_2_dd_complex, M12_oneline_4_vv_complex, M12_oneline_4_vd_complex, M12_oneline_mu2_vd_complex, M12_oneline_mu2_dd_complex, M12_oneline_mu4_vv_complex, M12_oneline_mu4_vd_complex, M12_oneline_mu6_vv_complex (the mu0_dd is just the simple standard one I load...). Then, for the bias stuff I need -> M12_2_bG2_oneline_complex, M12_2_b2_oneline_complex, M12_0_bG2_oneline_complex, M12_0_b1bG2_oneline_complex, M12_0_b2_oneline_complex, M12_0_b1b2_oneline_complex... AND ALSO M_fNLd2, M_fNLG2...

        class_alloc(pnlpt->M12_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pnlpt->M12_oneline_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        //GC -> this I keep... I.e., I could remove also these operations, but it becomes a [legit nightmare] of what to remove in the code besides this, and makes everything prone to new mistakes. In this way, I just say that these matrices are all equal to the matter-matter in real space... Actually, this must be done above... This remains as is... No, it can be donw below... Here I say that the complex matrices are constructed from the "non-complex" matter one...

        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vv0_f3, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vd0_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_dd0_f0, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vv2_f3, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vd2_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vv4_f3, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vd4_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vv0_f3_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vd0_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_dd0_f0_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vv2_f3_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vd2_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vv4_f3_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_multipoles_vd4_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        class_alloc(pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_mu_powers_dd2_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_mu_powers_vv4_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_mu_powers_vd4_f2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_mu_powers_vv6_f3, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_mu_powers_dd2_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_mu_powers_vv4_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_mu_powers_vd4_f2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_matter_mu_powers_vv6_f3_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        class_alloc(pnlpt->M12_oneline_complex_bias_real_space_b2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_bias_real_space_bG2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pnlpt->M12_oneline_complex_bias_real_space_b2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_bias_real_space_bG2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        class_alloc(pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        //GC!

        class_alloc(pnlpt->M12_oneline_0_vv_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_0_vd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_0_dd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_2_vv_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_2_vd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_2_dd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_4_vv_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_4_vd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_mu2_vd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_mu2_dd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_mu4_vv_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_mu4_vd_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_mu6_vv_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pnlpt->M12_oneline_0_vv_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_0_vd_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_0_dd_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_2_vv_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_2_vd_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_2_dd_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_4_vv_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_4_vd_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_mu2_vd_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_mu2_dd_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_mu4_vv_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_mu4_vd_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_oneline_mu6_vv_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        class_alloc(pnlpt->M12_2_bG2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_2_b2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_0_bG2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_0_b1bG2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_0_b2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_0_b1b2_oneline_complex, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pnlpt->M12_2_bG2_oneline_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_2_b2_oneline_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_0_bG2_oneline_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_0_b1bG2_oneline_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_0_b2_oneline_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M12_0_b1b2_oneline_complex_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        class_alloc(pnlpt->M_fNLd2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M_fNLG2, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pnlpt->M_fNLd2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);
        class_alloc(pnlpt->M_fNLG2_ortho, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        //GC -> my allocations end here...
        
        int count = 0;

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M22_oneline_complex[count] = pnlpt->M22_oneline[count] + _Complex_I * pnlpt->M22_oneline[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        // printf("%i\n",ppr->nmax_nlpt);

        for (count = 0; count < ppr->nmax_nlpt + 1; count++)
        {
            pnlpt->M13_oneline_complex[count] = pnlpt->M13_oneline[count] + _Complex_I * pnlpt->M13_oneline[count + ppr->nmax_nlpt + 1];
            //     printf("%le\n",creal( pnlpt->M13_oneline_complex[count]));
        }

        for (count = 0; count < ppr->nmax_nlpt + 1; count++)
        {
            pnlpt->IFG2_oneline_complex[count] = pnlpt->IFG2_oneline[count] + _Complex_I * pnlpt->IFG2_oneline[count + ppr->nmax_nlpt + 1];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M22basic_oneline_complex[count] = pnlpt->M22basic_oneline[count] + _Complex_I * pnlpt->M22basic_oneline[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        //GC!
        
        //GC - SWITCH HERE!
        //GC -> here my idea is to switch this completely off... Indeed, I do **not** care if these matrices are not filled with garbage... Whatever they are filled with, they should never generate segmentation faults? The blas routines should be able to handle without completely crashing, if I feed them garbage matrices...
        
        if (SWITCH_index == 1)
            
        {


        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex[count] = pnlpt->M12_oneline[count] + _Complex_I * pnlpt->M12_oneline[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        //GC: ORTHOGONAL -- start

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_ortho[count] = pnlpt->M12_oneline_ortho[count] + _Complex_I * pnlpt->M12_oneline_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        //GC: ORTHOGONAL -- finish

        /*-----------------------------------------------*/

        //GC -> here is where I make the modifications? No, see above? Yes, they are here... NOTICE THAT THIS IS NOT SUPER EFFICIENT... I am still making a lot of useless operations... But these should NOT AT ALL BE THE MAIN SOURCE OF SLOWING DOWN!!! So I will keep them for now... This can be optimized in the future...

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2[count] = pnlpt->M12_oneline_matter_multipoles_vv0_f2[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vv0_f2[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vv0_f3[count] = pnlpt->M12_oneline_matter_multipoles_vv0_f3[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vv0_f3[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1[count] = pnlpt->M12_oneline_matter_multipoles_vd0_f1[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vd0_f1[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vd0_f2[count] = pnlpt->M12_oneline_matter_multipoles_vd0_f2[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vd0_f2[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_dd0_f0[count] = pnlpt->M12_oneline_matter_multipoles_dd0_f0[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_dd0_f0[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1[count] = pnlpt->M12_oneline_matter_multipoles_dd0_f1[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_dd0_f1[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vv2_f3[count] = pnlpt->M12_oneline_matter_multipoles_vv2_f3[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vv2_f3[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vd2_f2[count] = pnlpt->M12_oneline_matter_multipoles_vd2_f2[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vd2_f2[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vv4_f3[count] = pnlpt->M12_oneline_matter_multipoles_vv4_f3[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vv4_f3[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vd4_f2[count] = pnlpt->M12_oneline_matter_multipoles_vd4_f2[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vd4_f2[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        //GC: ORTHOGONAL -- start

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2_ortho[count] = pnlpt->M12_oneline_matter_multipoles_vv0_f2_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vv0_f2_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vv0_f3_ortho[count] = pnlpt->M12_oneline_matter_multipoles_vv0_f3_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vv0_f3_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1_ortho[count] = pnlpt->M12_oneline_matter_multipoles_vd0_f1_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vd0_f1_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vd0_f2_ortho[count] = pnlpt->M12_oneline_matter_multipoles_vd0_f2_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vd0_f2_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_dd0_f0_ortho[count] = pnlpt->M12_oneline_matter_multipoles_dd0_f0_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_dd0_f0_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1_ortho[count] = pnlpt->M12_oneline_matter_multipoles_dd0_f1_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_dd0_f1_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vv2_f3_ortho[count] = pnlpt->M12_oneline_matter_multipoles_vv2_f3_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vv2_f3_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vd2_f2_ortho[count] = pnlpt->M12_oneline_matter_multipoles_vd2_f2_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vd2_f2_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vv4_f3_ortho[count] = pnlpt->M12_oneline_matter_multipoles_vv4_f3_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vv4_f3_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_multipoles_vd4_f2_ortho[count] = pnlpt->M12_oneline_matter_multipoles_vd4_f2_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_multipoles_vd4_f2_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        //GC: ORTHOGONAL -- finish

        /*-----------------------------------------------*/

        //GC -> here it remains the same...

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f1[count] = pnlpt->M12_oneline_matter_mu_powers_vd2_f1[count] + _Complex_I * pnlpt->M12_oneline_matter_mu_powers_vd2_f1[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2[count] = pnlpt->M12_oneline_matter_mu_powers_vd2_f2[count] + _Complex_I * pnlpt->M12_oneline_matter_mu_powers_vd2_f2[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_mu_powers_dd2_f1[count] = pnlpt->M12_oneline_matter_mu_powers_dd2_f1[count] + _Complex_I * pnlpt->M12_oneline_matter_mu_powers_dd2_f1[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        //GC -> every time that Xcode spazzes out you have likely lost a bracket...

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_mu_powers_vv4_f2[count] = pnlpt->M12_oneline_matter_mu_powers_vv4_f2[count] + _Complex_I * pnlpt->M12_oneline_matter_mu_powers_vv4_f2[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_mu_powers_vd4_f2[count] = pnlpt->M12_oneline_matter_mu_powers_vd4_f2[count] + _Complex_I * pnlpt->M12_oneline_matter_mu_powers_vd4_f2[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_mu_powers_vv6_f3[count] = pnlpt->M12_oneline_matter_mu_powers_vv6_f3[count] + _Complex_I * pnlpt->M12_oneline_matter_mu_powers_vv6_f3[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        //GC: ORTHOGONAL -- start

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f1_ortho[count] = pnlpt->M12_oneline_matter_mu_powers_vd2_f1_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_mu_powers_vd2_f1_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2_ortho[count] = pnlpt->M12_oneline_matter_mu_powers_vd2_f2_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_mu_powers_vd2_f2_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_mu_powers_dd2_f1_ortho[count] = pnlpt->M12_oneline_matter_mu_powers_dd2_f1_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_mu_powers_dd2_f1_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_mu_powers_vv4_f2_ortho[count] = pnlpt->M12_oneline_matter_mu_powers_vv4_f2_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_mu_powers_vv4_f2_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_mu_powers_vd4_f2_ortho[count] = pnlpt->M12_oneline_matter_mu_powers_vd4_f2_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_mu_powers_vd4_f2_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_matter_mu_powers_vv6_f3_ortho[count] = pnlpt->M12_oneline_matter_mu_powers_vv6_f3_ortho[count] + _Complex_I * pnlpt->M12_oneline_matter_mu_powers_vv6_f3_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        //GC: ORTHOGONAL -- finish

        /*-----------------------------------------------*/

        //GC -> here it changes...

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_bias_real_space_b2[count] = pnlpt->M12_oneline_bias_real_space_b2[count] + _Complex_I * pnlpt->M12_oneline_bias_real_space_b2[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_bias_real_space_bG2[count] = pnlpt->M12_oneline_bias_real_space_bG2[count] + _Complex_I * pnlpt->M12_oneline_bias_real_space_bG2[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        //GC: ORTHOGONAL -- start

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_bias_real_space_b2_ortho[count] = pnlpt->M12_oneline_bias_real_space_b2_ortho[count] + _Complex_I * pnlpt->M12_oneline_bias_real_space_b2_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_bias_real_space_bG2_ortho[count] = pnlpt->M12_oneline_bias_real_space_bG2_ortho[count] + _Complex_I * pnlpt->M12_oneline_bias_real_space_bG2_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        //GC: ORTHOGONAL -- finish

        /*-----------------------------------------------*/

        //GC -> here it remains the same...

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1[count] = pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1[count] + _Complex_I * pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1[count] = pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1[count] + _Complex_I * pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        //GC: ORTHOGONAL -- start

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1_ortho[count] = pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1_ortho[count] + _Complex_I * pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        for (count = 0; count < (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2; count++)
        {
            pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1_ortho[count] = pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1_ortho[count] + _Complex_I * pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1_ortho[count + (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2];
        }

        //GC: ORTHOGONAL -- finish
            
            
        }

        /*-----------------------------------------------*/

        //GC!

        //GC -> no more changes needed... There is nothing to deallocate or things like that... The only thing to deallocate was the matrices that are the "non-complex" ones...

        /* Pt matrices uploaded */

        /*It is used by spectra_pk_nl_bias_at_z (for classy) */
        //GC -> will need to play around here as well...
        class_alloc(pnlpt->ln_pk_nl, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_Id2d2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_Id2d2_2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_Id2d2_4, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        class_alloc(pnlpt->ln_pk_Id2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_IG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_Id2G2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_Id2G2_2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_Id2G2_4, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        class_alloc(pnlpt->ln_pk_IG2G2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_IG2G2_2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_IG2G2_4, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        class_alloc(pnlpt->ln_pk_IFG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_IFG2_0b1, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_IFG2_0, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_IFG2_2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_CTR, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_CTR_0, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_CTR_2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_CTR_4, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_Tree, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_Tree_0_vv, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_Tree_0_vd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_Tree_0_dd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_Tree_2_vv, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_Tree_2_vd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_Tree_4_vv, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_0_vv, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_0_vd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_0_dd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        class_alloc(pnlpt->ln_pk_2_vv, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_2_vd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_2_dd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        class_alloc(pnlpt->ln_pk_4_vv, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_4_vd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_4_dd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        class_alloc(pnlpt->ln_pk_0_b1b2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_0_b2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_0_b1bG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_0_bG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        class_alloc(pnlpt->ln_pk_2_b1b2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_2_b2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_2_b1bG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_2_bG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        class_alloc(pnlpt->ln_pk_4_b2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_4_bG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_4_b1b2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_4_b1bG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        //GC!

        //GC -> it does NOT require new coding, in the sense that below these are literally defined as the logs of things I compute... But still, I need to allocate them... They "eat" the ones with "_l_"...

        class_alloc(pnlpt->ln_pk_fNL_0_vv, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNL_0_vd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNL_0_dd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        class_alloc(pnlpt->ln_pk_fNL_2_vv, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNL_2_vd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNL_2_dd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        class_alloc(pnlpt->ln_pk_fNL_4_vv, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNL_4_vd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNL_4_dd, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pnlpt->ln_pk_fNL_0_vv_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNL_0_vd_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNL_0_dd_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        class_alloc(pnlpt->ln_pk_fNL_2_vv_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNL_2_vd_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNL_2_dd_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        class_alloc(pnlpt->ln_pk_fNL_4_vv_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNL_4_vd_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNL_4_dd_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        class_alloc(pnlpt->ln_pk12_0_b1b2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_0_b2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_0_b1bG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_0_bG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_2_b1b2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_2_b2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_2_b1bG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_2_bG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_4_b1b2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_4_b2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_4_b1bG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_4_bG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pnlpt->ln_pk12_0_b1b2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_0_b2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_0_b1bG2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_0_bG2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_2_b1b2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_2_b2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_2_b1bG2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_2_bG2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_4_b1b2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_4_b2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_4_b1bG2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk12_4_bG2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        class_alloc(pnlpt->ln_pk_nl_fNL, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNLd2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNLG2, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        class_alloc(pnlpt->ln_pk_nl_fNL_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNLd2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pnlpt->ln_pk_fNLG2_ortho, sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size, pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        //GC!
        
        //GC - SWITCH -> here I keep... Everything that allocates stuff I keep...

        int index_md;
        index_md = pnlpt->index_md_scalars;
        int index_ic = 0;

        // printf("ppt->k_size[index_md]=%i\n",ppt->k_size[index_md]);

        // printf("pnlpt->k[pnlpt->k_size-1]=%le\n",pnlpt->k[pnlpt->k_size-1]);
        // printf("ppt->k[index_md][ppt->k_size[index_md]-1]=%le\n",ppt->k[index_md][ppt->k_size[index_md]-1]);

        //int start2=clock();
        /*Begin of new part for nonlinear_pt_pk_l*/
        if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_)
        {
            if (pnlpt->cb == _TRUE_)
            {
                class_alloc(pnlpt->dd_sources_tp_delta_cb,
                            ppt->k_size[index_md] * ppt->tau_size * sizeof(double),
                            pnlpt->error_message);
                class_alloc(pnlpt->sources_tp_delta_cb,
                            ppt->k_size[index_md] * ppt->tau_size * sizeof(double),
                            pnlpt->error_message);

                class_call(array_spline_table_columns(ppt->k[index_md],
                                                      ppt->k_size[index_md],
                                                      ppt->sources[index_md]
                                                                  [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_cb],
                                                      ppt->tau_size,
                                                      pnlpt->dd_sources_tp_delta_cb,
                                                      _SPLINE_EST_DERIV_,
                                                      pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                for (index_tau = pnlpt->tau_size - 1; index_tau >= 0; index_tau--)
                {

                    for (index_k = 0; index_k < pnlpt->k_size; index_k++)
                    {

                        if (pnlpt->k[index_k] <= ppt->k[index_md][ppt->k_size[index_md] - 1])
                        {
                            class_call(array_interpolate_spline_one_column(ppt->k[index_md],
                                                                           ppt->k_size[index_md],
                                                                           ppt->sources[index_md]
                                                                                       [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_cb],
                                                                           pnlpt->tau_size,
                                                                           index_tau,
                                                                           pnlpt->dd_sources_tp_delta_cb,
                                                                           pnlpt->k[index_k],
                                                                           &pnlpt->sources_tp_delta_cb[index_tau * pnlpt->k_size + index_k],
                                                                           pnlpt->error_message),
                                       pnlpt->error_message,
                                       pnlpt->error_message);

                        } // condition for the new "if"
                        /* Very crude padding, but the physics should not depend on these scales anlyway */
                        /* This is done only for the non-zero spatial curvature, for which the kmax from perturbations does not coincude with the kmax from the spectra */
                        else
                        {
                            pnlpt->sources_tp_delta_cb[index_tau * pnlpt->k_size + index_k] = pnlpt->sources_tp_delta_cb[index_tau * pnlpt->k_size + index_k - 1];
                        }
                        // printf(" pnlpt->sources_tp_delta_b[index_tau*pnlpt->k_size+index_k] = %le\n", pnlpt->sources_tp_delta_b[index_tau*pnlpt->k_size+index_k]);
                    }
                }
            }
            else
            {
                class_alloc(pnlpt->dd_sources_tp_delta_m,
                            ppt->k_size[index_md] * ppt->tau_size * sizeof(double),
                            pnlpt->error_message);
                class_alloc(pnlpt->sources_tp_delta_m,
                            ppt->k_size[index_md] * ppt->tau_size * sizeof(double),
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
                for (index_tau = pnlpt->tau_size - 1; index_tau >= 0; index_tau--)
                {
                    for (index_k = 0; index_k < pnlpt->k_size; index_k++)
                    {

                        if (pnlpt->k[index_k] <= ppt->k[index_md][ppt->k_size[index_md] - 1])
                        {
                            class_call(array_interpolate_spline_one_column(ppt->k[index_md],
                                                                           ppt->k_size[index_md],
                                                                           ppt->sources[index_md]
                                                                                       [index_ic * ppt->tp_size[index_md] + ppt->index_tp_delta_m],
                                                                           pnlpt->tau_size,
                                                                           index_tau,
                                                                           pnlpt->dd_sources_tp_delta_m,
                                                                           pnlpt->k[index_k],
                                                                           &pnlpt->sources_tp_delta_m[index_tau * pnlpt->k_size + index_k],
                                                                           pnlpt->error_message),
                                       pnlpt->error_message,
                                       pnlpt->error_message);
                        } // condition for the new "if"
                        /* Very crude padding, but the physics should not depend on these scales anlyway */
                        /* This is done only for curvature, for which the kmax from perturbations does not coincide with the kmax from the spectra */
                        else
                        {
                            pnlpt->sources_tp_delta_m[index_tau * pnlpt->k_size + index_k] = pnlpt->sources_tp_delta_m[index_tau * pnlpt->k_size + index_k - 1];
                        }
                    }
                }
            }
            /*End cb*/
        }

        //GC -> this part above remains the same!!! Now the real thing begins...

        /*End of new part for nonlinear_pt_pk_l*/
        //int end2=clock();
        //printf("new part is run over %d ms\n",end2-start2);

        /*determination of full tables of P_L(k,tau) lnP_L(k,tau) and ddP_L(k,tau)*/
        double *lnpk_l_full;
        double *pk_l_full;
        double *ddlnpk_l_full;
        //double *lntk_l_full; //GC!
        //double *tk_l_full; //GC!
        //double *ddlntk_l_full; //GC!

        /*double *ddddlnpk_l_full;*/
        class_alloc(lnpk_l_full, sizeof(double) * pnlpt->tau_size * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pk_l_full, sizeof(double) * pnlpt->tau_size * pnlpt->k_size, pnlpt->error_message);
        class_alloc(ddlnpk_l_full, sizeof(double) * pnlpt->tau_size * pnlpt->k_size, pnlpt->error_message);
        /*class_alloc(ddddlnpk_l_full,sizeof(double)*pnlpt->tau_size*pnlpt->k_size,pnlpt->error_message);*/

        //class_alloc(lntk_l_full,sizeof(double)*pnlpt->tau_size*pnlpt->k_size,pnlpt->error_message); //GC!
        //class_alloc(tk_l_full,sizeof(double)*pnlpt->tau_size*pnlpt->k_size,pnlpt->error_message); //GC!
        //class_alloc(ddlntk_l_full,sizeof(double)*pnlpt->tau_size*pnlpt->k_size,pnlpt->error_message); //GC!

        //GC: the one BELOW is the cycle over \tau, and the function nonlinear_pt_pk_l that extracts the linear power spectrum is called repeatedly. The crucial point is that it is idiotic to call the PRIMORDIAL power spectrum as a function of \tau. So I need to create a new function. Notice that the allocation of the "_full" variables is an allocation of variables of size pnlpt->tau_size*pnlpt->k_size, so we do not need to allocate them for the primordial power...

        /*
      
      double *pPRIMk_l; //GC!
      double *lnpPRIMk_l; //GC!
      double *ddlnpPRIMk_l; //GC!

      
      class_alloc(pPRIMk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message); //GC!
      class_alloc(lnpPRIMk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message); //GC!
      class_alloc(ddlnpPRIMk_l,pnlpt->k_size*sizeof(double),pnlpt->error_message); //GC!

      */

        /*begin cycle over tau*/
        for (index_tau = pnlpt->tau_size - 1; index_tau >= 0; index_tau--)
        {

            /* get P_L(k) at this time */

            class_call(nonlinear_pt_pk_l(pba, ppt, ppm, pnlpt, index_tau, pk_l, lnk_l, lnpk_l, ddlnpk_l), //GC!
                                                                                                          //class_call(nonlinear_pt_pk_l(pba,ppt,ppm,pnlpt,index_tau,pk_l,lnk_l,lnpk_l,ddlnpk_l,tk_l,lntk_l,ddlntk_l),
                                                                                                          //class_call(nonlinear_pt_pk_l(pba,ppt,ppm,pnlpt,index_tau,pk_l,lnk_l,lnpk_l,ddlnpk_l,pPRIMk_l,lnpPRIMk_l,ddlnpPRIMk_l),
                       pnlpt->error_message,
                       pnlpt->error_message);

            /* get P_L(k,tau) lnP_L(k,tau) and ddP_L(k,tau) */

            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                lnpk_l_full[index_tau * pnlpt->k_size + index_k] = lnpk_l[index_k];
                pk_l_full[index_tau * pnlpt->k_size + index_k] = pk_l[index_k];
                ddlnpk_l_full[index_tau * pnlpt->k_size + index_k] = ddlnpk_l[index_k];
                //  pnlpt->nl_corr_density[index_tau * pnlpt->k_size + index_k] = 1.;

                //GC -> see that the _full ones are filled...

                //lntk_l_full[index_tau * pnlpt->k_size + index_k]=lntk_l[index_k]; //GC!
                //tk_l_full[index_tau * pnlpt->k_size + index_k]=tk_l[index_k]; //GC!
                //ddlntk_l_full[index_tau * pnlpt->k_size + index_k]=ddlntk_l[index_k]; //GC!
            }
        }
        /*end cycle over tau*/
        //int end1=clock();
        //printf("tau cycle is run over %d ms\n",end1-start1);

        //int end1=clock();
        //printf("tau cycle is run over %d ms\n",end1-start1);

        /*
      
       //GC -> need to call my new function now!!!
      
       */

        class_call(nonlinear_pt_pPRIMk_l(      //pba,
                       ppt, ppm, pnlpt, lnk_l, //lnk_l, //GC -> SHOULD ALREADY HAVE BEEN FILLED BEFORE? Boh, let's keep it anyway...
                       pPRIMk_l, lnpPRIMk_l, ddlnpPRIMk_l),
                   pnlpt->error_message,
                   pnlpt->error_message);

        //GC -> DO NOT NEED A LOOP HERE, since I just fill in pPRIMk_l,lnpPRIMk_l,ddlnpPRIMk_l with this function, and I do not have _full ones...

        /* interpolation of P_L(k,tau) lnP_L(k,tau) and ddP_L(k,tau) at tau_req[i_z]*/
        double *ln_pk_l_at_z_req;
        double *pk_l_at_z_req;

        //double *ln_tk_l_at_z_req; //GC!
        //double *tk_l_at_z_req; //GC!

        //GC -> here it allocates a pk_l_at_z_req quantity and its log that will be OVERWRITTEN in the cycle below on i_z... I will allocate a primordial PS that will not be overwritten, then...

        /*
      
      double *ln_pPRIMk_l_req;
      double *pPRIMk_l_req;
      
      
      */

        double *ln_pPRIMk_l_req;
        double *pPRIMk_l_req;

        class_alloc(ln_pPRIMk_l_req, sizeof(double) * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pPRIMk_l_req, sizeof(double) * pnlpt->k_size, pnlpt->error_message);

        //GC!

        double Dref = 0;

        class_alloc(ln_pk_l_at_z_req, sizeof(double) * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pk_l_at_z_req, sizeof(double) * pnlpt->k_size, pnlpt->error_message);

        //class_alloc(ln_tk_l_at_z_req,sizeof(double)*pnlpt->k_size,pnlpt->error_message); //GC!
        //class_alloc(tk_l_at_z_req,sizeof(double)*pnlpt->k_size,pnlpt->error_message); //GC!

        /* Retriveing the logarithmic growth factor f for RSD */

        class_alloc(pnlpt->growthf, sizeof(double) * pnlpt->z_pk_num, pnlpt->error_message);
        class_alloc(pnlpt->hratio_array, sizeof(double) * pnlpt->z_pk_num, pnlpt->error_message);
        class_alloc(pnlpt->Dratio_array, sizeof(double) * pnlpt->z_pk_num, pnlpt->error_message);

        //GC -> no need to do anything here...

        if (pnlpt->rsd == rsd_yes)
        {

            double *pvecbackf;
            int last_indexf;
            class_alloc(pvecbackf, pba->bg_size * sizeof(double), pnlpt->error_message);

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
                if (pnlpt->AP_effect == AP_effect_yes)
                {
                    printf("Computing the Alcock-Paczynski effect for fiducial cosmology with Om=%lf\n", Omfid);
                }
                else
                {
                    printf("No Alcock-Paczynski effect.\n");
                }

            for (i_z = 0; i_z < pnlpt->z_pk_num; i_z++)
            {

                double Da = 0;
                double Dfid = 0;
                double hnew = 0;
                double hfid = 0;

                class_call(background_at_tau(pba, tau_req[i_z], pba->long_info, pba->inter_normal, &last_indexf, pvecbackf),
                           pnlpt->error_message,
                           pnlpt->error_message);

                // Add new fz and Dref values from input file, if necessary.
                if (pnlpt->replace_background)
                {
                    pnlpt->growthf[i_z] = pnlpt->replace_fz_value;
                    Dref = pnlpt->replace_Dz_value;
                }
                else
                {
                    pnlpt->growthf[i_z] = pvecbackf[pba->index_bg_f];
                    Dref = pvecbackf[pba->index_bg_D];
                }

                //printf("Dref=%lf\n",Dref);

                if (pnlpt->AP_effect == AP_effect_yes)
                {

                    //printf("z_pk[i_z]=%lf\n",pnlpt->z_pk[i_z]);

                    if (pnlpt->z_pk[i_z] == 0.)
                    {
                        pnlpt->hratio_array[i_z] = 1.;
                        pnlpt->Dratio_array[i_z] = 1.;
                        hnew = 1.;
                    }
                    else
                    {
                        //hnew = pow((Omtrue*pow((1.+pnlpt->z_pk[i_z]),3.) + (1. - Omtrue)+(pba->Omega0_g)*pow((1.+pnlpt->z_pk[i_z]),4.)),0.5);
                        hfid = pow((Omfid * pow((1. + pnlpt->z_pk[i_z]), 3.) + (1. - Omfid) + (pba->Omega0_g) * pow((1. + pnlpt->z_pk[i_z]), 4.)), 0.5);
                        //printf("Hclass=%le\n",pvecbackf[pba->index_bg_H]/kmsMpc/100/pba->h);
                        //printf("Hmy = %le\n",hnew);
                        if (pnlpt->replace_background)
                        {
                            hnew = pnlpt->replace_Hz_value / kmsMpc / 100.0 / pba->h;
                        }
                        else
                        {
                            hnew = pvecbackf[pba->index_bg_H] / kmsMpc / 100.0 / pba->h;
                        }
                        pnlpt->hratio_array[i_z] = hnew / hfid;

                        // pnlpt->hratio_array[i_z] = pow(((pba->Omega0_cdm+pba->Omega0_b)*pow((1.+pnlpt->z_pk[i_z]),3.) + (1. - pba->Omega0_cdm - pba->Omega0_b)+(pba->Omega0_g)*pow((1.+pnlpt->z_pk[i_z]),4.)),0.5)/pow((Omfid*pow((1.+pnlpt->z_pk[i_z]),3.) + (1. - Omfid)+(pba->Omega0_g)*pow((1.+pnlpt->z_pk[i_z]),4.)),0.5);
                        dz = pnlpt->z_pk[i_z] / (1. * Nz - 1.);

                        for (j = 1; j < Nz; j++)
                        {
                            //Da = Da + dz*(1./pow((Omtrue*pow((1.+ dz*j),3.) + (1. -Omtrue)+(pba->Omega0_g)*pow((1.+pnlpt->z_pk[i_z]),4.)),0.5)+1./pow((Omtrue*pow((1.+ dz*(j-1)),3.) + (1. - Omtrue)+(pba->Omega0_g)*pow((1.+pnlpt->z_pk[i_z]),4.)),0.5))/2.;
                            Dfid = Dfid + dz * (1. / pow((Omfid * pow((1. + dz * j), 3.) + (1. - Omfid) + (pba->Omega0_g) * pow((1. + pnlpt->z_pk[i_z]), 4.)), 0.5) + 1. / pow((Omfid * pow((1. + dz * (j - 1)), 3.) + (1. - Omfid) + (pba->Omega0_g) * pow((1. + pnlpt->z_pk[i_z]), 4.)), 0.5)) / 2.;
                        }
                        //printf("Dclass=%le\n",pvecbackf[pba->index_bg_ang_distance]*kmsMpc*100*pba->h*(1.+pnlpt->z_pk[i_z]));
                        //printf("Dmy = %le\n",Da);
                        if (pnlpt->replace_background)
                        {
                            Da = pnlpt->replace_DAz_value * kmsMpc * 100 * pba->h * (1. + pnlpt->z_pk[i_z]);
                        }
                        else
                        {
                            Da = pvecbackf[pba->index_bg_ang_distance] * kmsMpc * 100 * pba->h * (1. + pnlpt->z_pk[i_z]);
                        }

                        pnlpt->Dratio_array[i_z] = Da / Dfid;
                        //    printf(" h/hfid=%le\n", hnew/hfid);
                        //    printf(" Da/Dfid=%le\n", Da/Dfid);
                    }
                }
                else
                {
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

                //GC!

                /*
    
    printf("Om=%.16e\n",pba->Omega0_cdm+pba->Omega0_b);
    printf("Omfid=%.16e\n",Omfid);
    printf("pvecbackf[pba->index_bg_a]=%.16e\n",pvecbackf[pba->index_bg_a]);
    printf("z_pk[i_z]=%.16e\n",pnlpt->z_pk[i_z]);
    printf("pvecbackf[pba->index_bg_H]=%.16e\n",pvecbackf[pba->index_bg_H]/kmsMpc);
    printf("pvecbackf[pba->index_bg_ang_distance]=%.16e\n",pvecbackf[pba->index_bg_ang_distance]*kmsMpc);
    printf("D = %.16e\n",Da/100/pba->h/(1.+pnlpt->z_pk[i_z])); //GC -> this is zero if you do not require the AP...
    printf("H = %.16e\n",hnew*100*pba->h); //GC -> this is zero if you do not require the AP...

    */

                //GC!
            }

            free(pvecbackf);
        }

        else
        {

            double *pvecbackf;
            int last_indexf = 0;
            class_alloc(pvecbackf, pba->bg_size * sizeof(double), pnlpt->error_message);

            for (i_z = 0; i_z < pnlpt->z_pk_num; i_z++)
            {
                pnlpt->growthf[i_z] = 1.;
                pnlpt->hratio_array[i_z] = 1.;
                pnlpt->Dratio_array[i_z] = 1.;
                //  printf("No RSD computed\n");

                for (i_z = 0; i_z < pnlpt->z_pk_num; i_z++)
                {
                    class_call(background_at_tau(pba, tau_req[i_z], pba->long_info, pba->inter_normal, &last_indexf, pvecbackf),
                               pnlpt->error_message,
                               pnlpt->error_message);
                    Dref = pvecbackf[pba->index_bg_D];
                }
            }
            free(pvecbackf);
        }

        //GC -> no need to do anything above!!!

        /* end of RSD specification */

        //ln_pk_l_at_z_req is the array of the linear PS

        /*
       
       
        //GC -> I think that also the loop below is useless, since it is a loop over a redshift index... And it fills the "_req" variables that I decided to not have for now...
       
       
       */

        last_index = 0;
        for (i_z = 0; i_z < pnlpt->z_pk_num; i_z++)
        {

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

            //    class_call(array_interpolate_spline(pnlpt->ln_tau,
            //                                                pnlpt->tau_size,
            //                                                lntk_l_full,
            //                                                ddlntk_l_full,
            //                                                pnlpt->k_size,
            //                                                log(tau_req[i_z]),
            //                                                &last_index,
            //                                                ln_tk_l_at_z_req,
            //                                                pnlpt->k_size,
            //                                                pnlpt->error_message),
            //                       pnlpt->error_message,
            //                       pnlpt->error_message); //GC!

            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                ln_pk_l_at_z_req[index_k] = ln_pk_l_at_z_req[index_k];
                pk_l_at_z_req[index_k] = exp(ln_pk_l_at_z_req[index_k]);
                //ln_tk_l_at_z_req[index_k] = ln_tk_l_at_z_req[index_k]; //GC!
                //tk_l_at_z_req[index_k]=exp(ln_tk_l_at_z_req[index_k]); //GC!
                //      printf("%le %le\n",khere,exp(ln_pk_l_at_z_req[index_k]));

                //GC!

                ln_pPRIMk_l_req[index_k] = lnpPRIMk_l[index_k];
                pPRIMk_l_req[index_k] = exp(ln_pPRIMk_l_req[index_k]);

                //GC -> I feel this is very redundant...

                //GC!
            }

            /* get P_NL(k) at tau_req */

            if (print_warning == _FALSE_)
            {
                //int start=clock();

                //GC: ORTHOGONAL -- start
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
                                             pPRIMk_l_req, //GC!
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
                                             pk_nl_fNL, //GC!
                                             pk_fNLd2,  //GC!
                                             pk_fNLG2,  //GC!
                                             //GC: ORTHOGONAL -- start
                                             pk_nl_fNL_ortho, //GC!
                                             pk_fNLd2_ortho,  //GC!
                                             pk_fNLG2_ortho,  //GC!
                                             //GC: ORTHOGONAL -- finish
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
                                             //GC!
                                             pk_l_fNL_0_vv,
                                             pk_l_fNL_0_vd,
                                             pk_l_fNL_0_dd,
                                             pk_l_fNL_2_vv,
                                             pk_l_fNL_2_vd,
                                             pk_l_fNL_2_dd,
                                             pk_l_fNL_4_vv,
                                             pk_l_fNL_4_vd,
                                             pk_l_fNL_4_dd,
                                             pk12_l_0_b1b2,
                                             pk12_l_0_b2,
                                             pk12_l_0_b1bG2,
                                             pk12_l_0_bG2,
                                             pk12_l_2_b1b2,
                                             pk12_l_2_b2,
                                             pk12_l_2_b1bG2,
                                             pk12_l_2_bG2,
                                             pk12_l_4_b1b2,
                                             pk12_l_4_b2,
                                             pk12_l_4_b1bG2,
                                             pk12_l_4_bG2,
                                             //GC!
                                             //GC: ORTHOGONAL -- start
                                             pk_l_fNL_0_vv_ortho,
                                             pk_l_fNL_0_vd_ortho,
                                             pk_l_fNL_0_dd_ortho,
                                             pk_l_fNL_2_vv_ortho,
                                             pk_l_fNL_2_vd_ortho,
                                             pk_l_fNL_2_dd_ortho,
                                             pk_l_fNL_4_vv_ortho,
                                             pk_l_fNL_4_vd_ortho,
                                             pk_l_fNL_4_dd_ortho,
                                             pk12_l_0_b1b2_ortho,
                                             pk12_l_0_b2_ortho,
                                             pk12_l_0_b1bG2_ortho,
                                             pk12_l_0_bG2_ortho,
                                             pk12_l_2_b1b2_ortho,
                                             pk12_l_2_b2_ortho,
                                             pk12_l_2_b1bG2_ortho,
                                             pk12_l_2_bG2_ortho,
                                             pk12_l_4_b1b2_ortho,
                                             pk12_l_4_b2_ortho,
                                             pk12_l_4_b1bG2_ortho,
                                             pk12_l_4_bG2_ortho,
                                             //GC: ORTHOGONAL -- finish
                                             lnk_l,
                                             ln_pk_l_at_z_req, //GC!!!
                                             ln_pPRIMk_l_req   //lnpPRIMk_l //GC!!!
                                             ),
                           pnlpt->error_message,
                           pnlpt->error_message);
                //int end=clock();

                //GC: this below needs to be accounted for...

                //if (pnlpt->nonlinear_pt_verbose > 0)
                //printf("Module nonlinear_pt_loop takes %d musec.\n",end-start);

                for (index_k = 0; index_k < pnlpt->k_size; index_k++)
                {
                    pnlpt->ln_pk_nl[i_z * pnlpt->k_size + index_k] = log(pk_nl[index_k]);
                    pnlpt->ln_pk_Id2d2[i_z * pnlpt->k_size + index_k] = log(pk_Id2d2[index_k]);
                    pnlpt->ln_pk_Id2d2_2[i_z * pnlpt->k_size + index_k] = log(pk_Id2d2_2[index_k]);
                    pnlpt->ln_pk_Id2d2_4[i_z * pnlpt->k_size + index_k] = log(pk_Id2d2_4[index_k]);

                    pnlpt->ln_pk_Id2[i_z * pnlpt->k_size + index_k] = log(pk_Id2[index_k]);
                    pnlpt->ln_pk_IG2[i_z * pnlpt->k_size + index_k] = log(pk_IG2[index_k]);
                    pnlpt->ln_pk_Id2G2[i_z * pnlpt->k_size + index_k] = log(pk_Id2G2[index_k]);
                    pnlpt->ln_pk_Id2G2_2[i_z * pnlpt->k_size + index_k] = log(pk_Id2G2_2[index_k]);
                    pnlpt->ln_pk_Id2G2_4[i_z * pnlpt->k_size + index_k] = log(pk_Id2G2_4[index_k]);

                    pnlpt->ln_pk_IG2G2[i_z * pnlpt->k_size + index_k] = log(pk_IG2G2[index_k]);
                    pnlpt->ln_pk_IG2G2_2[i_z * pnlpt->k_size + index_k] = log(pk_IG2G2_2[index_k]);
                    pnlpt->ln_pk_IG2G2_4[i_z * pnlpt->k_size + index_k] = log(pk_IG2G2_4[index_k]);

                    pnlpt->ln_pk_IFG2[i_z * pnlpt->k_size + index_k] = log(pk_IFG2[index_k]);

                    pnlpt->ln_pk_IFG2_0[i_z * pnlpt->k_size + index_k] = log(pk_IFG2_0[index_k]);
                    pnlpt->ln_pk_IFG2_0b1[i_z * pnlpt->k_size + index_k] = log(pk_IFG2_0b1[index_k]);
                    pnlpt->ln_pk_IFG2_2[i_z * pnlpt->k_size + index_k] = log(pk_IFG2_2[index_k]);

                    pnlpt->ln_pk_CTR[i_z * pnlpt->k_size + index_k] = log(pk_CTR[index_k]);
                    pnlpt->ln_pk_CTR_0[i_z * pnlpt->k_size + index_k] = log(pk_CTR_0[index_k]);

                    pnlpt->ln_pk_CTR_2[i_z * pnlpt->k_size + index_k] = log(pk_CTR_2[index_k]);
                    pnlpt->ln_pk_CTR_4[i_z * pnlpt->k_size + index_k] = log(pk_CTR_4[index_k]);

                    pnlpt->ln_pk_Tree[i_z * pnlpt->k_size + index_k] = log(pk_Tree[index_k]);
                    pnlpt->ln_pk_Tree_0_vv[i_z * pnlpt->k_size + index_k] = log(pk_Tree_0_vv[index_k]);
                    pnlpt->ln_pk_Tree_0_vd[i_z * pnlpt->k_size + index_k] = log(pk_Tree_0_vd[index_k]);
                    pnlpt->ln_pk_Tree_0_dd[i_z * pnlpt->k_size + index_k] = log(pk_Tree_0_dd[index_k]);
                    pnlpt->ln_pk_Tree_2_vv[i_z * pnlpt->k_size + index_k] = log(pk_Tree_2_vv[index_k]);
                    pnlpt->ln_pk_Tree_2_vd[i_z * pnlpt->k_size + index_k] = log(pk_Tree_2_vd[index_k]);
                    pnlpt->ln_pk_Tree_4_vv[i_z * pnlpt->k_size + index_k] = log(pk_Tree_4_vv[index_k]);
                    pnlpt->ln_pk_0_vv[i_z * pnlpt->k_size + index_k] = log(pk_l_0_vv[index_k]);
                    pnlpt->ln_pk_0_vd[i_z * pnlpt->k_size + index_k] = log(pk_l_0_vd[index_k]);
                    pnlpt->ln_pk_0_dd[i_z * pnlpt->k_size + index_k] = log(pk_l_0_dd[index_k]);
                    pnlpt->ln_pk_0_b1b2[i_z * pnlpt->k_size + index_k] = log(pk_l_0_b1b2[index_k]);
                    pnlpt->ln_pk_0_b2[i_z * pnlpt->k_size + index_k] = log(pk_l_0_b2[index_k]);
                    pnlpt->ln_pk_0_b1bG2[i_z * pnlpt->k_size + index_k] = log(pk_l_0_b1bG2[index_k]);
                    pnlpt->ln_pk_0_bG2[i_z * pnlpt->k_size + index_k] = log(pk_l_0_bG2[index_k]);

                    pnlpt->ln_pk_2_b1b2[i_z * pnlpt->k_size + index_k] = log(pk_l_2_b1b2[index_k]);
                    pnlpt->ln_pk_2_b2[i_z * pnlpt->k_size + index_k] = log(pk_l_2_b2[index_k]);
                    pnlpt->ln_pk_2_b1bG2[i_z * pnlpt->k_size + index_k] = log(pk_l_2_b1bG2[index_k]);
                    pnlpt->ln_pk_2_bG2[i_z * pnlpt->k_size + index_k] = log(pk_l_2_bG2[index_k]);

                    pnlpt->ln_pk_4_b2[i_z * pnlpt->k_size + index_k] = log(pk_l_4_b2[index_k]);
                    pnlpt->ln_pk_4_bG2[i_z * pnlpt->k_size + index_k] = log(pk_l_4_bG2[index_k]);
                    pnlpt->ln_pk_4_b1b2[i_z * pnlpt->k_size + index_k] = log(pk_l_4_b1b2[index_k]);
                    pnlpt->ln_pk_4_b1bG2[i_z * pnlpt->k_size + index_k] = log(pk_l_4_b1bG2[index_k]);

                    pnlpt->ln_pk_2_vv[i_z * pnlpt->k_size + index_k] = log(pk_l_2_vv[index_k]);
                    pnlpt->ln_pk_2_vd[i_z * pnlpt->k_size + index_k] = log(pk_l_2_vd[index_k]);
                    pnlpt->ln_pk_2_dd[i_z * pnlpt->k_size + index_k] = log(pk_l_2_dd[index_k]);

                    pnlpt->ln_pk_4_vv[i_z * pnlpt->k_size + index_k] = log(pk_l_4_vv[index_k]);
                    pnlpt->ln_pk_4_vd[i_z * pnlpt->k_size + index_k] = log(pk_l_4_vd[index_k]);
                    pnlpt->ln_pk_4_dd[i_z * pnlpt->k_size + index_k] = log(pk_l_4_dd[index_k]);

                    //GC!

                    pnlpt->ln_pk_fNL_0_vv[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_0_vv[index_k]);
                    pnlpt->ln_pk_fNL_0_vd[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_0_vd[index_k]);
                    pnlpt->ln_pk_fNL_0_dd[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_0_dd[index_k]);

                    pnlpt->ln_pk_fNL_2_vv[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_2_vv[index_k]);
                    pnlpt->ln_pk_fNL_2_vd[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_2_vd[index_k]);
                    pnlpt->ln_pk_fNL_2_dd[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_2_dd[index_k]);

                    pnlpt->ln_pk_fNL_4_vv[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_4_vv[index_k]);
                    pnlpt->ln_pk_fNL_4_vd[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_4_vd[index_k]);
                    pnlpt->ln_pk_fNL_4_dd[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_4_dd[index_k]);

                    pnlpt->ln_pk12_0_b1b2[i_z * pnlpt->k_size + index_k] = log(pk12_l_0_b1b2[index_k]);
                    pnlpt->ln_pk12_0_b2[i_z * pnlpt->k_size + index_k] = log(pk12_l_0_b2[index_k]);
                    pnlpt->ln_pk12_0_b1bG2[i_z * pnlpt->k_size + index_k] = log(pk12_l_0_b1bG2[index_k]);
                    pnlpt->ln_pk12_0_bG2[i_z * pnlpt->k_size + index_k] = log(pk12_l_0_bG2[index_k]);
                    pnlpt->ln_pk12_2_b1b2[i_z * pnlpt->k_size + index_k] = log(pk12_l_2_b1b2[index_k]);
                    pnlpt->ln_pk12_2_b2[i_z * pnlpt->k_size + index_k] = log(pk12_l_2_b2[index_k]);
                    pnlpt->ln_pk12_2_b1bG2[i_z * pnlpt->k_size + index_k] = log(pk12_l_2_b1bG2[index_k]);
                    pnlpt->ln_pk12_2_bG2[i_z * pnlpt->k_size + index_k] = log(pk12_l_2_bG2[index_k]);
                    pnlpt->ln_pk12_4_b1b2[i_z * pnlpt->k_size + index_k] = log(pk12_l_4_b1b2[index_k]);
                    pnlpt->ln_pk12_4_b2[i_z * pnlpt->k_size + index_k] = log(pk12_l_4_b2[index_k]);
                    pnlpt->ln_pk12_4_b1bG2[i_z * pnlpt->k_size + index_k] = log(pk12_l_4_b1bG2[index_k]);
                    pnlpt->ln_pk12_4_bG2[i_z * pnlpt->k_size + index_k] = log(pk12_l_4_bG2[index_k]);

                    pnlpt->ln_pk_nl_fNL[i_z * pnlpt->k_size + index_k] = log(pk_nl_fNL[index_k]);
                    pnlpt->ln_pk_fNLd2[i_z * pnlpt->k_size + index_k] = log(pk_fNLd2[index_k]);
                    pnlpt->ln_pk_fNLG2[i_z * pnlpt->k_size + index_k] = log(pk_fNLG2[index_k]);

                    //GC: ORTHOGONAL -- start

                    pnlpt->ln_pk_fNL_0_vv_ortho[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_0_vv_ortho[index_k]);
                    pnlpt->ln_pk_fNL_0_vd_ortho[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_0_vd_ortho[index_k]);
                    pnlpt->ln_pk_fNL_0_dd_ortho[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_0_dd_ortho[index_k]);

                    pnlpt->ln_pk_fNL_2_vv_ortho[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_2_vv_ortho[index_k]);
                    pnlpt->ln_pk_fNL_2_vd_ortho[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_2_vd_ortho[index_k]);
                    pnlpt->ln_pk_fNL_2_dd_ortho[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_2_dd_ortho[index_k]);

                    pnlpt->ln_pk_fNL_4_vv_ortho[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_4_vv_ortho[index_k]);
                    pnlpt->ln_pk_fNL_4_vd_ortho[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_4_vd_ortho[index_k]);
                    pnlpt->ln_pk_fNL_4_dd_ortho[i_z * pnlpt->k_size + index_k] = log(pk_l_fNL_4_dd_ortho[index_k]);

                    pnlpt->ln_pk12_0_b1b2_ortho[i_z * pnlpt->k_size + index_k] = log(pk12_l_0_b1b2_ortho[index_k]);
                    pnlpt->ln_pk12_0_b2_ortho[i_z * pnlpt->k_size + index_k] = log(pk12_l_0_b2_ortho[index_k]);
                    pnlpt->ln_pk12_0_b1bG2_ortho[i_z * pnlpt->k_size + index_k] = log(pk12_l_0_b1bG2_ortho[index_k]);
                    pnlpt->ln_pk12_0_bG2_ortho[i_z * pnlpt->k_size + index_k] = log(pk12_l_0_bG2_ortho[index_k]);
                    pnlpt->ln_pk12_2_b1b2_ortho[i_z * pnlpt->k_size + index_k] = log(pk12_l_2_b1b2_ortho[index_k]);
                    pnlpt->ln_pk12_2_b2_ortho[i_z * pnlpt->k_size + index_k] = log(pk12_l_2_b2_ortho[index_k]);
                    pnlpt->ln_pk12_2_b1bG2_ortho[i_z * pnlpt->k_size + index_k] = log(pk12_l_2_b1bG2_ortho[index_k]);
                    pnlpt->ln_pk12_2_bG2_ortho[i_z * pnlpt->k_size + index_k] = log(pk12_l_2_bG2_ortho[index_k]);
                    pnlpt->ln_pk12_4_b1b2_ortho[i_z * pnlpt->k_size + index_k] = log(pk12_l_4_b1b2_ortho[index_k]);
                    pnlpt->ln_pk12_4_b2_ortho[i_z * pnlpt->k_size + index_k] = log(pk12_l_4_b2_ortho[index_k]);
                    pnlpt->ln_pk12_4_b1bG2_ortho[i_z * pnlpt->k_size + index_k] = log(pk12_l_4_b1bG2_ortho[index_k]);
                    pnlpt->ln_pk12_4_bG2_ortho[i_z * pnlpt->k_size + index_k] = log(pk12_l_4_bG2_ortho[index_k]);

                    pnlpt->ln_pk_nl_fNL_ortho[i_z * pnlpt->k_size + index_k] = log(pk_nl_fNL_ortho[index_k]);
                    pnlpt->ln_pk_fNLd2_ortho[i_z * pnlpt->k_size + index_k] = log(pk_fNLd2_ortho[index_k]);
                    pnlpt->ln_pk_fNLG2_ortho[i_z * pnlpt->k_size + index_k] = log(pk_fNLG2_ortho[index_k]);

                    //GC: ORTHOGONAL -- finish

                    //GC!
                }
            }
        }

        //printf("Dlast=%f \n",Dref);

        /* This is required for lensing with one-loop ! */

        /*Begin lensing module*/
        if (ppt->has_cl_cmb_lensing_potential == _TRUE_)
        {
            /*Lensing turn on*/
            double Dplus = 0;
            double *pvecbackD;
            int last_indexD = 0;
            class_alloc(pvecbackD, pba->bg_size * sizeof(double), pnlpt->error_message);

            double *dd_pk_Tree;
            double *dd_pk_nl;
            double *dd_pk_l_at_z_req;
            double *pk_Tree_int;
            double *pk_nl_int;
            double *pk_l_at_z_req_int;
            // double * pk_ctr_int;

            last_index = 0;
            
            double large_for_logs_matter = 5. * pow(10., 3.); //GC!

            if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_)
            {
                class_alloc(dd_pk_Tree, pnlpt->k_size * sizeof(double), pnlpt->error_message);
                class_alloc(dd_pk_nl, pnlpt->k_size * sizeof(double), pnlpt->error_message);
                class_alloc(dd_pk_l_at_z_req, pnlpt->k_size * sizeof(double), pnlpt->error_message);
                class_alloc(pk_Tree_int, ppt->k_size[pnlpt->index_md_scalars] * sizeof(double), pnlpt->error_message);
                class_alloc(pk_nl_int, ppt->k_size[pnlpt->index_md_scalars] * sizeof(double), pnlpt->error_message);
                class_alloc(pk_l_at_z_req_int, ppt->k_size[pnlpt->index_md_scalars] * sizeof(double), pnlpt->error_message);
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
                for (index_k = 0; index_k < ppt->k_size[pnlpt->index_md_scalars]; index_k++)
                {
                    if (ppt->k[pnlpt->index_md_scalars][index_k] < pnlpt->k[pnlpt->k_size - 1])
                    {
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
                    }
                    else
                    {
                        pk_Tree_int[index_k] = pk_Tree_int[index_k - 1];
                        pk_nl_int[index_k] = pk_nl_int[index_k - 1];
                        pk_l_at_z_req_int[index_k] = pk_l_at_z_req_int[index_k - 1];
                        // pk_ctr_int[index_k] = 0.;
                        // pk_ctr_int[index_k] = 2*pk_Tree_int[index_k]*ppt->k[pnlpt->index_md_scalars][index_k-1]*ppt->k[pnlpt->index_md_scalars][index_k-1]/(1.+ppt->k[pnlpt->index_md_scalars][index_k-1]*ppt->k[pnlpt->index_md_scalars][index_k-1]);
                    }
                }
            }

            for (index_tau = pnlpt->tau_size - 1; index_tau >= 0; index_tau--)
            {

                /*
class_call(nonlinear_pt_pk_l(pba,ppt,ppm,pnlpt,index_tau,pk_l,lnk_l,lnpk_l,ddlnpk_l),
                     pnlpt->error_message,
                     pnlpt->error_message);
                     */

                class_call(background_at_tau(pba, pnlpt->tau[index_tau], pba->long_info, pba->inter_normal, &last_indexD, pvecbackD), pnlpt->error_message, pnlpt->error_message);

                Dplus = pvecbackD[pba->index_bg_D];
                //printf("D+=%f \n",Dplus);
                if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_)
                {
                    for (index_k = 0; index_k < ppt->k_size[pnlpt->index_md_scalars]; index_k++)
                        //pnlpt->nl_corr_density[index_tau * ppt->k_size[pnlpt->index_md_scalars] + index_k] = 1.;
                        pnlpt->nl_corr_density[index_tau * ppt->k_size[pnlpt->index_md_scalars] + index_k] = sqrt(fabs((pk_Tree_int[index_k] + Dplus * Dplus * (pk_nl_int[index_k] - large_for_logs_matter
                                                                                                                                                                // -pk_ctr_int[index_k]
                                                                                                                                                                // -2.*pow(ppt->k[pnlpt->index_md_scalars][index_k],2.)*pk_Tree_int[index_k]/(1.+pow(ppt->k[pnlpt->index_md_scalars][index_k],2.))
                                                                                                                                                                ) /
                                                                                                                                                   Dref / Dref)) /
                                                                                                                  pk_l_at_z_req_int[index_k]);
                }
                else
                {
                    for (index_k = 0; index_k < ppt->k_size[pnlpt->index_md_scalars]; index_k++)
                        pnlpt->nl_corr_density[index_tau * ppt->k_size[pnlpt->index_md_scalars] + index_k] = sqrt(fabs((pk_Tree[index_k] + Dplus * Dplus * (pk_nl[index_k] - large_for_logs_matter
                                                                                                                                                            // -pk_ctr_int[index_k]
                                                                                                                                                            // -2.*pow(ppt->k[pnlpt->index_md_scalars][index_k],2.)*pk_Tree_int[index_k]/(1.+pow(ppt->k[pnlpt->index_md_scalars][index_k],2.))
                                                                                                                                                            ) /
                                                                                                                                               Dref / Dref)) /
                                                                                                                  pk_l_at_z_req[index_k]);
                }
                //printf("index_tau=%i \n",index_tau);
                //printf("index_k=%i \n",index_k);
                //printf("pnlpt->nl_corr_density[index_tau * pnlpt->k_size + index_k]=%f \n",pnlpt->nl_corr_density[index_tau * pnlpt->k_size + index_k]);
            }
            free(pvecbackD);

            if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_)
            {
                free(dd_pk_Tree);
                free(dd_pk_nl);
                free(dd_pk_l_at_z_req);

                free(pk_Tree_int);
                // free(pk_ctr_int);
                free(pk_nl_int);
                free(pk_l_at_z_req_int);
            }
        }
        else
        {
            /*Lensing turn off*/
            for (index_tau = pnlpt->tau_size - 1; index_tau >= 0; index_tau--)
            {
                for (index_k = 0; index_k < ppt->k_size[pnlpt->index_md_scalars]; index_k++)
                {
                    pnlpt->nl_corr_density[index_tau * ppt->k_size[pnlpt->index_md_scalars] + index_k] = 1.;
                }
            }
        }

        //GC -> THE ABOVE ONE I DO NOT TOUCH!!! In principle of course it should be modified, since at one-loop I change the matter PS...
        //GC: THIS WILL BE DONE IN FUTURE WORK!!!

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

        //free(tk_l); //GC!
        //free(lntk_l); //GC!
        //free(lntk_l_full); //GC!
        //free(tk_l_full); //GC!
        //free(ddlntk_l_full); //GC!
        //free(ln_tk_l_at_z_req); //GC!
        //free(tk_l_at_z_req); //GC!

        //GC!

        free(ln_pPRIMk_l_req); //GC: CAREFUL!!!
        free(pPRIMk_l_req);    //GC: CAREFUL!!!

        free(pk_l_fNL_0_vv);
        free(pk_l_fNL_0_vd);
        free(pk_l_fNL_0_dd);

        free(pk_l_fNL_2_vv);
        free(pk_l_fNL_2_vd);
        free(pk_l_fNL_2_dd);

        free(pk_l_fNL_4_vv);
        free(pk_l_fNL_4_vd);
        free(pk_l_fNL_4_dd);

        free(pk12_l_0_b1b2);
        free(pk12_l_0_b2);
        free(pk12_l_0_b1bG2);
        free(pk12_l_0_bG2);
        free(pk12_l_2_b1b2);
        free(pk12_l_2_b2);
        free(pk12_l_2_b1bG2);
        free(pk12_l_2_bG2);
        free(pk12_l_4_b1b2);
        free(pk12_l_4_b2);
        free(pk12_l_4_b1bG2);
        free(pk12_l_4_bG2);

        free(pk_nl_fNL);
        free(pk_fNLd2);
        free(pk_fNLG2);

        //GC: ORTHOGONAL -- start

        free(pk_l_fNL_0_vv_ortho);
        free(pk_l_fNL_0_vd_ortho);
        free(pk_l_fNL_0_dd_ortho);

        free(pk_l_fNL_2_vv_ortho);
        free(pk_l_fNL_2_vd_ortho);
        free(pk_l_fNL_2_dd_ortho);

        free(pk_l_fNL_4_vv_ortho);
        free(pk_l_fNL_4_vd_ortho);
        free(pk_l_fNL_4_dd_ortho);

        free(pk12_l_0_b1b2_ortho);
        free(pk12_l_0_b2_ortho);
        free(pk12_l_0_b1bG2_ortho);
        free(pk12_l_0_bG2_ortho);
        free(pk12_l_2_b1b2_ortho);
        free(pk12_l_2_b2_ortho);
        free(pk12_l_2_b1bG2_ortho);
        free(pk12_l_2_bG2_ortho);
        free(pk12_l_4_b1b2_ortho);
        free(pk12_l_4_b2_ortho);
        free(pk12_l_4_b1bG2_ortho);
        free(pk12_l_4_bG2_ortho);

        free(pk_nl_fNL_ortho);
        free(pk_fNLd2_ortho);
        free(pk_fNLG2_ortho);

        //GC: ORTHOGONAL -- finish

        free(pPRIMk_l);
        free(lnpPRIMk_l);
        //GC -> not these if I use them as input later? No, I already by now called it, so free everything...
        free(ddlnpPRIMk_l);

        //GC!

        //GC -> now everything should be freed...

        if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_)
        {
            if (pnlpt->cb == _TRUE_)
            {
                free(pnlpt->dd_sources_tp_delta_cb);
                free(pnlpt->sources_tp_delta_cb);
            }
            else
            {
                free(pnlpt->dd_sources_tp_delta_m);
                free(pnlpt->sources_tp_delta_m);
            }
        }

        if (pnlpt->nonlinear_pt_verbose > 0)
            printf(" 'nonlinear_pt_init' module executed successfully\n");
    }

    else
    {
        class_stop(pnlpt->error_message,
                   "Your non-linear method variable is set to %d, out of the range defined in nonlinear_pt.h", pnlpt->method);
    }

    //end of the loop over time before 'else' above

    return _SUCCESS_;
}

//GC -> nothing to fix above...

//end of non-linear init

int nonlinear_pt_free(
    struct nonlinear_pt *pnlpt)
{

    if (pnlpt->method > nlpt_none)
    {

        if (pnlpt->method == nlpt_spt)
        {
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

            //GC!

            //GC -> what do I need? Here I only had the new matrices to free...

            free(pnlpt->M12_oneline); //GC: CAREFUL!!!

            //GC: ORTHOGONAL -- start

            free(pnlpt->M12_oneline_ortho); //GC: CAREFUL!!!

            //GC: ORTHOGONAL -- finish

            //GC!

            /**/
            free(pnlpt->M12_oneline_matter_multipoles_vv0_f2);
            free(pnlpt->M12_oneline_matter_multipoles_vv0_f3);
            free(pnlpt->M12_oneline_matter_multipoles_vd0_f1);
            free(pnlpt->M12_oneline_matter_multipoles_vd0_f2);
            free(pnlpt->M12_oneline_matter_multipoles_dd0_f0);
            free(pnlpt->M12_oneline_matter_multipoles_dd0_f1);
            free(pnlpt->M12_oneline_matter_multipoles_vv2_f3);
            free(pnlpt->M12_oneline_matter_multipoles_vd2_f2);
            free(pnlpt->M12_oneline_matter_multipoles_vv4_f3);
            free(pnlpt->M12_oneline_matter_multipoles_vd4_f2);
            /**/

            //GC: ORTHOGONAL -- start

            free(pnlpt->M12_oneline_matter_multipoles_vv0_f2_ortho);
            free(pnlpt->M12_oneline_matter_multipoles_vv0_f3_ortho);
            free(pnlpt->M12_oneline_matter_multipoles_vd0_f1_ortho);
            free(pnlpt->M12_oneline_matter_multipoles_vd0_f2_ortho);
            free(pnlpt->M12_oneline_matter_multipoles_dd0_f0_ortho);
            free(pnlpt->M12_oneline_matter_multipoles_dd0_f1_ortho);
            free(pnlpt->M12_oneline_matter_multipoles_vv2_f3_ortho);
            free(pnlpt->M12_oneline_matter_multipoles_vd2_f2_ortho);
            free(pnlpt->M12_oneline_matter_multipoles_vv4_f3_ortho);
            free(pnlpt->M12_oneline_matter_multipoles_vd4_f2_ortho);

            //GC: ORTHOGONAL -- finish

            free(pnlpt->M12_oneline_matter_mu_powers_vd2_f1);
            free(pnlpt->M12_oneline_matter_mu_powers_vd2_f2);
            free(pnlpt->M12_oneline_matter_mu_powers_dd2_f1);
            free(pnlpt->M12_oneline_matter_mu_powers_vv4_f2);
            free(pnlpt->M12_oneline_matter_mu_powers_vd4_f2);
            free(pnlpt->M12_oneline_matter_mu_powers_vv6_f3);

            //GC: ORTHOGONAL -- start

            free(pnlpt->M12_oneline_matter_mu_powers_vd2_f1_ortho);
            free(pnlpt->M12_oneline_matter_mu_powers_vd2_f2_ortho);
            free(pnlpt->M12_oneline_matter_mu_powers_dd2_f1_ortho);
            free(pnlpt->M12_oneline_matter_mu_powers_vv4_f2_ortho);
            free(pnlpt->M12_oneline_matter_mu_powers_vd4_f2_ortho);
            free(pnlpt->M12_oneline_matter_mu_powers_vv6_f3_ortho);

            //GC: ORTHOGONAL -- finish

            /**/
            free(pnlpt->M12_oneline_bias_real_space_b2);
            free(pnlpt->M12_oneline_bias_real_space_bG2);
            /**/

            //GC: ORTHOGONAL -- start

            free(pnlpt->M12_oneline_bias_real_space_b2_ortho);
            free(pnlpt->M12_oneline_bias_real_space_bG2_ortho);

            //GC: ORTHOGONAL -- finish

            free(pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1);
            free(pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1);

            //GC: ORTHOGONAL -- start

            free(pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1_ortho);
            free(pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1_ortho);

            //GC: ORTHOGONAL -- finish

            //GC!

            //GC -> EVERYTHING BELOW REMAINS... Because the "complex" matrices still exist!!! It is clear that this was the simplest way to modify the code in such a way that I do not spend time reading the disk... It could be optimized more, but it is something for another time...

            free(pnlpt->M12_oneline_complex);

            //GC: ORTHOGONAL -- start

            free(pnlpt->M12_oneline_complex_ortho);

            //GC: ORTHOGONAL -- finish

            //GC!

            free(pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vv0_f3);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vd0_f2);
            free(pnlpt->M12_oneline_complex_matter_multipoles_dd0_f0);
            free(pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vv2_f3);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vd2_f2);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vv4_f3);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vd4_f2);

            free(pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f1);
            free(pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2);
            free(pnlpt->M12_oneline_complex_matter_mu_powers_dd2_f1);
            free(pnlpt->M12_oneline_complex_matter_mu_powers_vv4_f2);
            free(pnlpt->M12_oneline_complex_matter_mu_powers_vd4_f2);
            free(pnlpt->M12_oneline_complex_matter_mu_powers_vv6_f3);

            free(pnlpt->M12_oneline_complex_bias_real_space_b2);
            free(pnlpt->M12_oneline_complex_bias_real_space_bG2);

            free(pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1);
            free(pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1);

            //GC: ORTHOGONAL -- start

            free(pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2_ortho);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vv0_f3_ortho);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1_ortho);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vd0_f2_ortho);
            free(pnlpt->M12_oneline_complex_matter_multipoles_dd0_f0_ortho);
            free(pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1_ortho);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vv2_f3_ortho);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vd2_f2_ortho);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vv4_f3_ortho);
            free(pnlpt->M12_oneline_complex_matter_multipoles_vd4_f2_ortho);

            free(pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f1_ortho);
            free(pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2_ortho);
            free(pnlpt->M12_oneline_complex_matter_mu_powers_dd2_f1_ortho);
            free(pnlpt->M12_oneline_complex_matter_mu_powers_vv4_f2_ortho);
            free(pnlpt->M12_oneline_complex_matter_mu_powers_vd4_f2_ortho);
            free(pnlpt->M12_oneline_complex_matter_mu_powers_vv6_f3_ortho);

            free(pnlpt->M12_oneline_complex_bias_real_space_b2_ortho);
            free(pnlpt->M12_oneline_complex_bias_real_space_bG2_ortho);

            free(pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1_ortho);
            free(pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1_ortho);

            //GC: ORTHOGONAL -- finish

            //GC!

            free(pnlpt->M12_oneline_0_vv_complex);
            free(pnlpt->M12_oneline_0_vd_complex);
            free(pnlpt->M12_oneline_0_dd_complex);
            free(pnlpt->M12_oneline_2_vv_complex);
            free(pnlpt->M12_oneline_2_vd_complex);
            free(pnlpt->M12_oneline_2_dd_complex);
            free(pnlpt->M12_oneline_4_vv_complex);
            free(pnlpt->M12_oneline_4_vd_complex);
            free(pnlpt->M12_oneline_mu2_vd_complex);
            free(pnlpt->M12_oneline_mu2_dd_complex);
            free(pnlpt->M12_oneline_mu4_vv_complex);
            free(pnlpt->M12_oneline_mu4_vd_complex);
            free(pnlpt->M12_oneline_mu6_vv_complex);

            free(pnlpt->M12_2_bG2_oneline_complex);
            free(pnlpt->M12_2_b2_oneline_complex);
            free(pnlpt->M12_0_bG2_oneline_complex);
            free(pnlpt->M12_0_b1bG2_oneline_complex);
            free(pnlpt->M12_0_b2_oneline_complex);
            free(pnlpt->M12_0_b1b2_oneline_complex);

            free(pnlpt->M_fNLd2);
            free(pnlpt->M_fNLG2);

            //GC: ORTHOGONAL -- start

            free(pnlpt->M12_oneline_0_vv_complex_ortho);
            free(pnlpt->M12_oneline_0_vd_complex_ortho);
            free(pnlpt->M12_oneline_0_dd_complex_ortho);
            free(pnlpt->M12_oneline_2_vv_complex_ortho);
            free(pnlpt->M12_oneline_2_vd_complex_ortho);
            free(pnlpt->M12_oneline_2_dd_complex_ortho);
            free(pnlpt->M12_oneline_4_vv_complex_ortho);
            free(pnlpt->M12_oneline_4_vd_complex_ortho);
            free(pnlpt->M12_oneline_mu2_vd_complex_ortho);
            free(pnlpt->M12_oneline_mu2_dd_complex_ortho);
            free(pnlpt->M12_oneline_mu4_vv_complex_ortho);
            free(pnlpt->M12_oneline_mu4_vd_complex_ortho);
            free(pnlpt->M12_oneline_mu6_vv_complex_ortho);

            free(pnlpt->M12_2_bG2_oneline_complex_ortho);
            free(pnlpt->M12_2_b2_oneline_complex_ortho);
            free(pnlpt->M12_0_bG2_oneline_complex_ortho);
            free(pnlpt->M12_0_b1bG2_oneline_complex_ortho);
            free(pnlpt->M12_0_b2_oneline_complex_ortho);
            free(pnlpt->M12_0_b1b2_oneline_complex_ortho);

            free(pnlpt->M_fNLd2_ortho);
            free(pnlpt->M_fNLG2_ortho);

            //GC: ORTHOGONAL -- finish

            //GC!

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

            //GC!

            //GC -> here only to free the logs for classy...

            free(pnlpt->ln_pk_fNL_0_vv);
            free(pnlpt->ln_pk_fNL_0_vd);
            free(pnlpt->ln_pk_fNL_0_dd);

            free(pnlpt->ln_pk_fNL_2_vv);
            free(pnlpt->ln_pk_fNL_2_vd);
            free(pnlpt->ln_pk_fNL_2_dd);

            free(pnlpt->ln_pk_fNL_4_vv);
            free(pnlpt->ln_pk_fNL_4_vd);
            free(pnlpt->ln_pk_fNL_4_dd);

            free(pnlpt->ln_pk12_0_b1b2);
            free(pnlpt->ln_pk12_0_b2);
            free(pnlpt->ln_pk12_0_b1bG2);
            free(pnlpt->ln_pk12_0_bG2);
            free(pnlpt->ln_pk12_2_b1b2);
            free(pnlpt->ln_pk12_2_b2);
            free(pnlpt->ln_pk12_2_b1bG2);
            free(pnlpt->ln_pk12_2_bG2);
            free(pnlpt->ln_pk12_4_b1b2);
            free(pnlpt->ln_pk12_4_b2);
            free(pnlpt->ln_pk12_4_b1bG2);
            free(pnlpt->ln_pk12_4_bG2);

            free(pnlpt->ln_pk_nl_fNL);
            free(pnlpt->ln_pk_fNLd2);
            free(pnlpt->ln_pk_fNLG2);

            //GC: ORTHOGONAL -- start

            free(pnlpt->ln_pk_fNL_0_vv_ortho);
            free(pnlpt->ln_pk_fNL_0_vd_ortho);
            free(pnlpt->ln_pk_fNL_0_dd_ortho);

            free(pnlpt->ln_pk_fNL_2_vv_ortho);
            free(pnlpt->ln_pk_fNL_2_vd_ortho);
            free(pnlpt->ln_pk_fNL_2_dd_ortho);

            free(pnlpt->ln_pk_fNL_4_vv_ortho);
            free(pnlpt->ln_pk_fNL_4_vd_ortho);
            free(pnlpt->ln_pk_fNL_4_dd_ortho);

            free(pnlpt->ln_pk12_0_b1b2_ortho);
            free(pnlpt->ln_pk12_0_b2_ortho);
            free(pnlpt->ln_pk12_0_b1bG2_ortho);
            free(pnlpt->ln_pk12_0_bG2_ortho);
            free(pnlpt->ln_pk12_2_b1b2_ortho);
            free(pnlpt->ln_pk12_2_b2_ortho);
            free(pnlpt->ln_pk12_2_b1bG2_ortho);
            free(pnlpt->ln_pk12_2_bG2_ortho);
            free(pnlpt->ln_pk12_4_b1b2_ortho);
            free(pnlpt->ln_pk12_4_b2_ortho);
            free(pnlpt->ln_pk12_4_b1bG2_ortho);
            free(pnlpt->ln_pk12_4_bG2_ortho);

            free(pnlpt->ln_pk_nl_fNL_ortho);
            free(pnlpt->ln_pk_fNLd2_ortho);
            free(pnlpt->ln_pk_fNLG2_ortho);

            //GC: ORTHOGONAL -- finish

            //GC!

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
    double *ddlnpk)
{
    //, //GC!
    //double *pPRIMk_l, //GC!
    //double *lnpPRIMk, //GC!
    //double *ddlnpPRIMk) //GC!
    //double *tk_l, //GC! -> this maybe must remain even if we have the square root??? Very likely...
    //double *lntk, //GC!
    //double *ddlntk) //GC!
    //{ //GC!

    int index_md;
    int index_ic = 0;
    int index_type;
    int index_k;
    int index_ic1, index_ic2, index_ic1_ic2;
    double *primordial_pk;
    double source_ic1, source_ic2;

    // Replace linear spectrum with a given input file
    if (pnlpt->replace_pk)
    {
        // open input interpolation file
        char line[100000];
        FILE *fp;
        errno = 0;
        fp = fopen(pnlpt->input_pk, "r");

        if (fp == NULL)
        {
            fprintf(stderr, "Interpolation file %s not found\n", pnlpt->input_pk);
            fprintf(stderr, "Trying again in 3sec\n");
            sleep(3);
            fp = fopen(pnlpt->input_pk, "r");
            if (fp == NULL)
            {
                fprintf(stderr, "Interpolation file %s still not found\n Aborting!\n", pnlpt->input_pk);
                fprintf(stderr, "The error code is %d.\n", errno);
                fprintf(stderr, "Error opening file: %s\n", strerror(errno));
                exit(0);
            }
        }
        // Count lines to construct the correct size
        int nline = 0;
        while (fgets(line, 10000, fp) != NULL)
        {
            if (line[0] == '#')
                continue; // comment line
            if (line[0] == '\n')
                continue;
            nline++;
        }
        rewind(fp); // restart file

        // Now allocate memory to the weights array
        double kint[nline], pkint[nline], kmin, kmax;

        int line_count = 0; // line counter
        int counter = 0;    // counts which element in line

        // Read in values to file
        while (fgets(line, 100000, fp) != NULL)
        {
            // Select required lines in file
            if (line[0] == '#')
                continue;
            if (line[0] == '\n')
                continue;

            // Split into variables
            char *split_string;
            split_string = strtok(line, "\t");
            counter = 0;

            // Iterate over line
            while (split_string != NULL)
            {
                if (counter == 0)
                {
                    kint[line_count] = atof(split_string);
                    //fprintf(stderr,"%.4e\n",kint[line_count]);
                }
                if (counter == 1)
                {
                    pkint[line_count] = atof(split_string);
                }
                if (counter > 1)
                {
                    fprintf(stderr, "Incorrect file format");
                    abort();
                }
                split_string = strtok(NULL, "\t");
                counter++;
            }
            line_count++;
        }
        kmin = kint[0];
        kmax = kint[nline - 1];
        double this_pk[1], this_k;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            lnk[index_k] = log(pnlpt->k[index_k]);

            this_k = pnlpt->k[index_k];
            this_pk[0] = -1;
            if (this_k < kmin)
                this_pk[0] = pkint[0];
            else if (this_k > kmax)
                this_pk[0] = pkint[nline - 1];
            else
            {
                int segment;
                for (segment = 0; segment < nline - 1; segment++)
                {
                    //printf("%d %d\n",segment,nline);
                    //printf("%.4e %.4e %.4e\n",kint[segment],kint[segment+1],this_k);
                    if ((kint[segment + 1] >= this_k) && (kint[segment] <= this_k))
                    {
                        this_pk[0] = pkint[segment] + (pkint[segment + 1] - pkint[segment]) * (this_k - kint[segment]) / (kint[segment + 1] - kint[segment]);
                        break;
                    }
                }
                if (this_pk[0] == -1)
                    abort();
            }
            pk_l[index_k] = this_pk[0];
            lnpk[index_k] = log(this_pk[0]);
        }
        fclose(fp);
    }
    // Compute linear power from CLASS as usual
    else
    {

        index_md = ppt->index_md_scalars;

        class_alloc(primordial_pk, ppm->ic_ic_size[index_md] * sizeof(double), pnlpt->error_message);

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            class_call(primordial_spectrum_at_k(ppm,
                                                index_md,
                                                linear,
                                                pnlpt->k[index_k],
                                                primordial_pk),
                       ppm->error_message,
                       pnlpt->error_message);

            pk_l[index_k] = 0;
            //pPRIMk_l[index_k] = 0; //GC! Also, notice that this was copied roughly straightforwardly from nonlinear.c...

            /* part diagonal in initial conditions */
            for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++)
            {
                index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic1, ppm->ic_size[index_md]);

                if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_)
                {
                    if (pnlpt->cb == _TRUE_)
                    {
                        source_ic1 = pnlpt->sources_tp_delta_cb[index_tau * pnlpt->k_size + index_k];
                    }
                    else
                    {
                        source_ic1 = pnlpt->sources_tp_delta_m[index_tau * pnlpt->k_size + index_k];
                    }
                }
                else
                {
                    if (pnlpt->cb == _TRUE_)
                    {
                        source_ic1 = ppt->sources[index_md][index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_cb][index_tau * ppt->k_size[index_md] + index_k];
                    }
                    else
                    {
                        source_ic1 = ppt->sources[index_md][index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m][index_tau * ppt->k_size[index_md] + index_k];
                    }
                }
                //printf("Omegacdm=%f   Omegab=%f",pba->Omega0_cdm,pba->Omega0_b);
                //printf("k=%f source_m=%f   source_cdmb=%f\n",pnlpt->k[index_k],ppt->sources[index_md][index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m][index_tau * ppt->k_size[index_md] + index_k],source_ic1);

                // source_ic1 are transfer functions

                pk_l[index_k] += 2. * _PI_ * _PI_ / pow(pnlpt->k[index_k], 3) * source_ic1 * source_ic1 * primordial_pk[index_ic1_ic2];

                //tk_l[index_k] += source_ic1; //GC! [[This must be checked!!!]]
                //pPRIMk_l[index_k] += 2.*_PI_*_PI_/pow(pnlpt->k[index_k],3)*primordial_pk[index_ic1_ic2]; //GC!
            }

            /* part non-diagonal in initial conditions */
            for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++)
            {
                for (index_ic2 = index_ic1 + 1; index_ic2 < ppm->ic_size[index_md]; index_ic2++)
                {
                    index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic2, ppm->ic_size[index_md]);
                    if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_)
                    {
                        source_ic1 = ppt->sources[index_md]
                                                 [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
                                                 [index_tau * ppt->k_size[index_md] + index_k];
                        source_ic2 = ppt->sources[index_md]
                                                 [index_ic2 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
                                                 [index_tau * ppt->k_size[index_md] + index_k];
                        pk_l[index_k] += 2. * 2. * _PI_ * _PI_ / pow(pnlpt->k[index_k], 3) * source_ic1 * source_ic2 * primordial_pk[index_ic1_ic2]; // extra 2 factor (to include the symmetric term ic2,ic1)
                        //tk_l[index_k] += 2.*source_ic1; //GC! [[This must be checked!!!]]
                        //pPRIMk_l[index_k] += 2.*2.*_PI_*_PI_/pow(pnlpt->k[index_k],3)*primordial_pk[index_ic1_ic2]; //GC!
                    }
                }
            }

            lnk[index_k] = log(pnlpt->k[index_k]);
            lnpk[index_k] = log(pk_l[index_k]);
            //lntk[index_k] = log(tk_l[index_k]); //GC!
            //tk_l[index_k] = pow(pk_l[index_k]/2./_PI_/_PI_*pow(pnlpt->k[index_k],3)/...,0.5); //GC -> check that square root is same as tricky method with source_ic1. Recall that primordial P is P_\phi, not P_\zeta... Need to check the 5/3 however... Notice that they will NOT be the same because sum of squares different from square of sums... Use the square root and keep track of 5/3... NOTICE THAT THIS IS TRICKY BECAUSE I NEED TO KEEP TRACK OF THE PRIMORDIAL PS ALSO... If I could extract the total matter transfer function directly, come on... Maybe I can email Tram? Notice that this sum is NOT an else, it is summing the diagonal terms plus the off-diagonal ones... The problem is that the sum runs also over the primordial_pk, so "extracting it out" is difficult... Maybe one can bypass this completely and NOT modify this function at all??? Boh... No, because I still need to remove the primordial PS... So maybe here the output should be the primordial PS only... In this way I can obtain the transfer function by dividing the power spectrum by the primordial PS. The primordial PS is obtained by doing the same loop, but without source_ic1 and source_ic2 stuff... And then I can obtain everything... Maybe...
            //GC -> this way of computing pPRIM and then using it to get the transfer function is WRONG if I have isocurvature. That is, if in the diagonal piece I have anything but \zeta, AND in the non-diagonal terms I have anything but zero. Not having any isocurvature should mean that all those are zero... Notice that the "cb" thing, which has only CDM and baryons, is an IF statement... So essentially I would have <(source_ic1*ic_1 + source_ic2*ic_2 + ...)^2> = source_ic1^2<ic_1^2> + source_ic2^2<ic_2^2> + 2*source_ic1*source_ic2<ic_1ic_2> + ... But now all ic past 1 are zero, so this is automatically what I want once I take the square root. So the only thing I need to understand is if primordial_pk has P_\zeta. This seems to be the case? If so, to get P_\phi I need some other factors, but at that point it is better to rewrite the overall factor to match to \zeta... How do I check whether this is P_\zeta? I look at primordial.c, and the slides, and this is FOR SURE the power spectrum of {\cal R}... Great...
            //lntk[index_k] = log(tk_l[index_k]); //GC!

            //lnpPRIMk[index_k] = log(pPRIMk_l[index_k]); //GC!

            //printf("%lf",pk_l[index_k]);
            //printf("\n");
        }
        free(primordial_pk);
    }

    class_call(array_spline_table_columns(lnk,
                                          pnlpt->k_size,
                                          lnpk,
                                          1, //GC -> means that only 1D interpolation...
                                          ddlnpk,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    /*
    
    class_call(array_spline_table_columns(lnk,
                                          pnlpt->k_size,
                                          lnpPRIMk,
                                          1,
                                          ddlnpPRIMk,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message); //GC!
     
     */

    return _SUCCESS_;

} // end of non-linear_Pk_l. This function executes the linear power spectrum

//GC -> AFTER ALL THE MILLION COMMENTS, THE FUNCTION ABOVE IS THE SAME AS THE ORIGINAL!!!

//GC!

int nonlinear_pt_pPRIMk_l(
    //struct background *pba,
    struct perturbs *ppt,
    struct primordial *ppm,
    struct nonlinear_pt *pnlpt,
    //int index_tau,
    //double *pk_l,
    double *lnk,
    //double *lnpk,
    //double *ddlnpk,
    double *pPRIMk_l,   //GC!
    double *lnpPRIMk,   //GC!
    double *ddlnpPRIMk) //GC!
                        //double *tk_l, //GC! -> this maybe must remain even if we have the square root??? Very likely...
                        //double *lntk, //GC!
                        //double *ddlntk) //GC!
{

    int index_md;
    int index_ic = 0;
    int index_type;
    int index_k;
    int index_ic1, index_ic2, index_ic1_ic2;
    double *primordial_pk;
    //double source_ic1,source_ic2; //GC!

    //GC: this part is not needed...

    /*
                                      
  // Replace linear spectrum with a given input file //
  if (pnlpt->replace_pk){
        // open input interpolation file
        char line[100000];
        FILE *fp;
        fp = fopen(pnlpt->input_pk,"r");

        if (fp==NULL){
            fprintf(stderr,"Interpolation file %s not found\n",pnlpt->input_pk);
            abort();
        }
        // Count lines to construct the correct size
        int nline = 0;
        while (fgets(line,10000,fp)!=NULL){
            if (line[0]=='#') continue; // comment line
            if (line[0]=='\n') continue;
                nline++;
            }
        rewind(fp); // restart file

        // Now allocate memory to the weights array
        double kint[nline], pkint[nline], kmin, kmax;

        int line_count=0; // line counter
        int counter=0; // counts which element in line

        // Read in values to file
        while (fgets(line,100000,fp)!=NULL) {
            // Select required lines in file
            if (line[0]=='#') continue;
            if (line[0]=='\n') continue;

            // Split into variables
            char * split_string;
            split_string = strtok(line, "\t");
            counter=0;

            // Iterate over line
            while (split_string!=NULL){
                if(counter==0){
                    kint[line_count]=atof(split_string);
                    //fprintf(stderr,"%.4e\n",kint[line_count]);
                    }
                if(counter==1){
                    pkint[line_count]=atof(split_string);
                    }
                if(counter>1){
                    fprintf(stderr,"Incorrect file format");
                    abort();
                }
                split_string = strtok(NULL,"\t");
                counter++;
            }
            line_count++;
        }
        kmin = kint[0];
        kmax = kint[nline-1];
        double this_pk[1], this_k;

        for (index_k=0; index_k<pnlpt->k_size; index_k++) {

          lnk[index_k] = log(pnlpt->k[index_k]);

          this_k = pnlpt->k[index_k];
          this_pk[0] = -1;
          if (this_k<kmin) this_pk[0] = pkint[0];
          else if (this_k>kmax) this_pk[0] = pkint[nline-1];
          else{
            int segment;
            for (segment = 0; segment < nline-1;segment++){
              //printf("%d %d\n",segment,nline);
              //printf("%.4e %.4e %.4e\n",kint[segment],kint[segment+1],this_k);
              if ((kint[segment+1]>=this_k)&&(kint[segment]<=this_k)){
                  this_pk[0] = pkint[segment]+(pkint[segment+1]-pkint[segment])*(this_k-kint[segment])/(kint[segment+1]-kint[segment]);
                  break;
              }
            }
            if (this_pk[0]==-1) abort();
          }
          pk_l[index_k] = this_pk[0];
          lnpk[index_k] = log(this_pk[0]);
            
          pPRIMk_l[index_k] = this_pk[0];
          lnpPRIMk[index_k] = log(this_pk[0]); //GC! -> hackfraud way to avoid issues if code enters here... Maybe it is worth to make this better? I can actually make it pass some input transfer or something. Let's see with Misha... This is NOT a priority at all... We agreed to leave this as is, when sampling this is NEVER USED...
            
        }
  }
  // Compute linear power from CLASS as usual
  
    */

    //else{

    index_md = ppt->index_md_scalars;

    class_alloc(primordial_pk, ppm->ic_ic_size[index_md] * sizeof(double), pnlpt->error_message);

    for (index_k = 0; index_k < pnlpt->k_size; index_k++)
    {

        class_call(primordial_spectrum_at_k(ppm,
                                            index_md,
                                            linear,
                                            pnlpt->k[index_k],
                                            primordial_pk),
                   ppm->error_message,
                   pnlpt->error_message);

        //pk_l[index_k] = 0; //GC!
        pPRIMk_l[index_k] = 0; //GC! Also, notice that this was copied roughly straightforwardly from nonlinear.c...

        /* part diagonal in initial conditions */
        for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++)
        {
            index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic1, ppm->ic_size[index_md]);

            /*
            
        if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {
          if (pnlpt->cb == _TRUE_) {
            source_ic1 = pnlpt->sources_tp_delta_cb[index_tau*pnlpt->k_size+index_k];
          } else {
            source_ic1 = pnlpt->sources_tp_delta_m[index_tau*pnlpt->k_size+index_k];
          }
        } else {
          if (pnlpt->cb == _TRUE_) {
            source_ic1 = ppt->sources[index_md][index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_cb][index_tau * ppt->k_size[index_md] + index_k];
          } else {
            source_ic1 = ppt->sources[index_md][index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m][index_tau * ppt->k_size[index_md] + index_k];
          }
        }
             
             */

            //printf("Omegacdm=%f   Omegab=%f",pba->Omega0_cdm,pba->Omega0_b);
            //printf("k=%f source_m=%f   source_cdmb=%f\n",pnlpt->k[index_k],ppt->sources[index_md][index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m][index_tau * ppt->k_size[index_md] + index_k],source_ic1);

            // source_ic1 are transfer functions

            //pk_l[index_k] += 2.*_PI_*_PI_/pow(pnlpt->k[index_k],3)
            //  *source_ic1*source_ic1
            //  *primordial_pk[index_ic1_ic2];

            //GC!

            //tk_l[index_k] += source_ic1; //GC! [[This must be checked!!!]]
            //pPRIMk_l[index_k] += 2.*_PI_*_PI_/pow(pnlpt->k[index_k],3)*primordial_pk[index_ic1_ic2]/ppm->A_s; //GC!
            pPRIMk_l[index_k] += 2. * _PI_ * _PI_ * (pnlpt->k[index_k]) * primordial_pk[index_ic1_ic2] / ppm->A_s; //GC!
        }

        /* part non-diagonal in initial conditions */
        for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++)
        {
            for (index_ic2 = index_ic1 + 1; index_ic2 < ppm->ic_size[index_md]; index_ic2++)
            {
                index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic2, ppm->ic_size[index_md]);
                if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_)
                {

                    /*
                source_ic1 = ppt->sources[index_md]
                [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
                [index_tau * ppt->k_size[index_md] + index_k];
              source_ic2 = ppt->sources[index_md]
                [index_ic2 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
                [index_tau * ppt->k_size[index_md] + index_k];
                 
                 */

                    //pk_l[index_k] += 2.*2.*_PI_*_PI_/pow(pnlpt->k[index_k],3)
                    //  *source_ic1*source_ic2
                    //  *primordial_pk[index_ic1_ic2]; // extra 2 factor (to include the symmetric term ic2,ic1)
                    //tk_l[index_k] += 2.*source_ic1; //GC! [[This must be checked!!!]]
                    //pPRIMk_l[index_k] += 2.*2.*_PI_*_PI_/pow(pnlpt->k[index_k],3)*primordial_pk[index_ic1_ic2]/ppm->A_s; //GC!
                    pPRIMk_l[index_k] += 2. * 2. * _PI_ * _PI_ * (pnlpt->k[index_k]) * primordial_pk[index_ic1_ic2] / ppm->A_s; //GC!
                }
            }
        }

        lnk[index_k] = log(pnlpt->k[index_k]);
        //lnpk[index_k] = log(pk_l[index_k]);
        //lntk[index_k] = log(tk_l[index_k]); //GC!
        //tk_l[index_k] = pow(pk_l[index_k]/2./_PI_/_PI_*pow(pnlpt->k[index_k],3)/...,0.5); //GC -> check that square root is same as tricky method with source_ic1. Recall that primordial P is P_\phi, not P_\zeta... Need to check the 5/3 however... Notice that they will NOT be the same because sum of squares different from square of sums... Use the square root and keep track of 5/3... NOTICE THAT THIS IS TRICKY BECAUSE I NEED TO KEEP TRACK OF THE PRIMORDIAL PS ALSO... If I could extract the total matter transfer function directly, come on... Maybe I can email Tram? Notice that this sum is NOT an else, it is summing the diagonal terms plus the off-diagonal ones... The problem is that the sum runs also over the primordial_pk, so "extracting it out" is difficult... Maybe one can bypass this completely and NOT modify this function at all??? Boh... No, because I still need to remove the primordial PS... So maybe here the output should be the primordial PS only... In this way I can obtain the transfer function by dividing the power spectrum by the primordial PS. The primordial PS is obtained by doing the same loop, but without source_ic1 and source_ic2 stuff... And then I can obtain everything... Maybe...
        //GC -> this way of computing pPRIM and then using it to get the transfer function is WRONG if I have isocurvature. That is, if in the diagonal piece I have anything but \zeta, AND in the non-diagonal terms I have anything but zero. Not having any isocurvature should mean that all those are zero... Notice that the "cb" thing, which has only CDM and baryons, is an IF statement... So essentially I would have <(source_ic1*ic_1 + source_ic2*ic_2 + ...)^2> = source_ic1^2<ic_1^2> + source_ic2^2<ic_2^2> + 2*source_ic1*source_ic2<ic_1ic_2> + ... But now all ic past 1 are zero, so this is automatically what I want once I take the square root. So the only thing I need to understand is if primordial_pk has P_\zeta. This seems to be the case? If so, to get P_\phi I need some other factors, but at that point it is better to rewrite the overall factor to match to \zeta... How do I check whether this is P_\zeta? I look at primordial.c, and the slides, and this is FOR SURE the power spectrum of {\cal R}... Great...
        //lntk[index_k] = log(tk_l[index_k]); //GC!

        lnpPRIMk[index_k] = log(pPRIMk_l[index_k]); //GC!

        //printf("%lf",pPRIMk_l[index_k]); //GC -> plotting this does NOT give me the linear PS, cavolo!!! A questo punto estraggo semplicemente A_s * stuff... In this way I also have under control h powers... Here pnlpt->k[index_k] is in units of h/Mpc: I need to check this... I then need: ppm->ns, ppm->k_pivot and ppm->A_s, I would say... Notice that Misha calls the primordial PS at pnlpt->k[index_k]. Whatever it is, it is correct... I need to think more about my case with the pivot... I don't really get how he takes into account the powers of h, since CLASS works with Mpc...
        //printf("\n");

        //printf("%lf %lf",pnlpt->k[index_k],pPRIMk_l[index_k]);
        //printf("%.16f %.16f",pnlpt->k[index_k],pPRIMk_l[index_k]);
        //printf("\n");
    }

    free(primordial_pk);

    //}

    /*
                                      
    class_call(array_spline_table_columns(lnk,
                                          pnlpt->k_size,
                                          lnpk,
                                          1, //GC -> means that only 1D interpolation...
                                          ddlnpk,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
     
     */

    class_call(array_spline_table_columns(lnk,
                                          pnlpt->k_size,
                                          lnpPRIMk,
                                          1,
                                          ddlnpPRIMk,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message); //GC!

    return _SUCCESS_;

} // end of non-linear_Pk_l. This function executes the linear power spectrum

/* beginning of the main function */

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
    double *pk_l,     //GC -> means this linear... Maybe this is not used -> it is never used below here... For IR resummation take square roots and stop at linear order in (1+x)^1/2...
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
    double *pk_nl,     //GC -> this means 1-loop matter PS...
    double *pk_nl_fNL, //GC!
    double *pk_fNLd2,  //GC!
    double *pk_fNLG2,  //GC!
    //GC: ORTHOGONAL -- start
    double *pk_nl_fNL_ortho, //GC!
    double *pk_fNLd2_ortho,  //GC!
    double *pk_fNLG2_ortho,  //GC!
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
    double *pk12_l_0_b2, //GC: there was a space between * and pk12_l_0_b2: irrelevant of course!!!
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
    double *lnpPRIMk_l)
{
    //, //GC -> this is log of linear PS...
    /* double *lnpPRIMk_l //GC: maybe this is not needed since I do not use it as input??? No, I need it... */
    //) {

    int index_k = 0;
    int index_j = 0;
    int index_i = 0;
    int index_l = 0;

    double sigmav = 0.;
    double *pvecback;

    int last_index = 0;
    int last_index2 = 0;

    class_alloc(pvecback, pba->bg_size * sizeof(double), pnlpt->error_message);

    // class_call only calls a function!

    class_call(background_at_tau(pba, tau, pba->long_info, pba->inter_normal, &last_index, pvecback),
               pba->error_message,
               pnlpt->error_message);

    // printf("pnlpt->cb=%i\n",pnlpt->cb);

    if (pnlpt->nonlinear_pt_verbose > 0)
        printf("Computing one-loop power spectra at z=%e\n", pba->a_today / pvecback[pba->index_bg_a] - 1.);
    //printf("Alcock-Paczynski effect!\n");

    free(pvecback);

    int Nmax = ppr->nmax_nlpt;
    int Nmaxf = ppr->nmax_nlpt + 1;
    double Nmaxd = Nmax * 1.;
    double kmin = 0.00005 * pba->h;
    double kmax = 100. * pba->h;

    /* If you generate new PT matrices, don't forget to choose kmin and kmax appropriately ! */

    //double kmin = pnlpt->k[0];
    //double kmax = pnlpt->k[pnlpt->k_size-1];

    double *js;
    class_alloc(js, (Nmax + 1) * sizeof(double), pnlpt->error_message);
    double *kdisc;
    class_alloc(kdisc, Nmax * sizeof(double), pnlpt->error_message);
    double *Pbin;
    class_alloc(Pbin, Nmax * sizeof(double), pnlpt->error_message);
    double *Pdisc;
    class_alloc(Pdisc, Nmax * sizeof(double), pnlpt->error_message);
    double *Ptree;
    class_alloc(Ptree, Nmax * sizeof(double), pnlpt->error_message);

    //GC!
    //double *PPRIMbin;
    //class_alloc(PPRIMbin,Nmax * sizeof(double),pnlpt->error_message); //GC: should be unused...
    double *PPRIMdisc;
    class_alloc(PPRIMdisc, Nmax * sizeof(double), pnlpt->error_message); //GC: CAREFUL!!!
    //double *PPRIMtree;
    //class_alloc(PPRIMtree,Nmax * sizeof(double),pnlpt->error_message); //GC: should be unused...

    //GC!
    double *Tbin;
    class_alloc(Tbin, Nmax * sizeof(double), pnlpt->error_message); //GC: CAREFUL!!!
    //double *Tdisc;
    //class_alloc(Tdisc,Nmax * sizeof(double),pnlpt->error_message); //GC -> useless...
    //double *Ttree;
    //class_alloc(Ttree,Nmax * sizeof(double),pnlpt->error_message);
    //GC: NOT NEEDED!!!

    /*
double *pk_10;
class_alloc(pk_10,pnlpt->k_size * sizeof(double),pnlpt->error_message);
double *pk_10b1;
class_alloc(pk_10b1,pnlpt->k_size * sizeof(double),pnlpt->error_message);
double *pk_12;
class_alloc(pk_12,pnlpt->k_size * sizeof(double),pnlpt->error_message);*/

    double *P10;
    class_alloc(P10, Nmax * sizeof(double), pnlpt->error_message);
    double *P10b1;
    class_alloc(P10b1, Nmax * sizeof(double), pnlpt->error_message);
    double *P12;
    class_alloc(P12, Nmax * sizeof(double), pnlpt->error_message);
    //GC: this shouldn't need to change...

    /*
     int i_kdisc1;
     for (i_kdisc1=0; i_kdisc1< Nmax; i_kdisc1++){
         P_IFG2_0b1_x[i_kdisc1] = 0.;
         P_IFG2_0[i_kdisc1] = 0.;
         P_IFG2_2[i_kdisc1] = 0.;
     }
     */

    double *Pnw;
    class_alloc(Pnw, Nmax * sizeof(double), pnlpt->error_message);

    double *Pw;
    class_alloc(Pw, Nmax * sizeof(double), pnlpt->error_message);

    double *dd_Pnw;
    class_alloc(dd_Pnw, sizeof(double) * Nmax, pnlpt->error_message);
    double *dd_Pw;
    class_alloc(dd_Pw, sizeof(double) * Nmax, pnlpt->error_message);

    //GC!

    double *Tnw;
    class_alloc(Tnw, Nmax * sizeof(double), pnlpt->error_message); //GC: CAREFUL!!!

    double *Tw;
    class_alloc(Tw, Nmax * sizeof(double), pnlpt->error_message); //GC: CAREFUL!!!

    double *dd_Tnw;
    class_alloc(dd_Tnw, sizeof(double) * Nmax, pnlpt->error_message); //GC: CAREFUL!!!
    double *dd_Tw;
    class_alloc(dd_Tw, sizeof(double) * Nmax, pnlpt->error_message); //GC: CAREFUL!!!
    //GC: should be unused? No, needed: notice that it is with a "_"...

    double SigmaBAO = 0.;
    double deltaSigmaBAO = 0.;
    // double kmin = pnl->k[0];
    double Delta = log(kmax / kmin) / (Nmaxd - 1);

    // Interpolating the linear power spectrum
    //GC: and primordial one, to get the transfer function...

    //GC -> something wrong happens here, before computing the FFT coefficients...

    double lnpk_out = 0.;
    double *myddlnpk;
    double lnpPRIMk_out = 0.; //GC!
    double *myddlnpPRIMk;     //GC!
    class_alloc(myddlnpk, sizeof(double) * pnlpt->k_size, pnlpt->error_message);
    class_call(array_spline_table_columns(lnk_l,
                                          pnlpt->k_size,
                                          lnpk_l,
                                          1,
                                          myddlnpk,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    class_alloc(myddlnpPRIMk, sizeof(double) * pnlpt->k_size, pnlpt->error_message); //GC!
    class_call(array_spline_table_columns(lnk_l,
                                          pnlpt->k_size,
                                          lnpPRIMk_l,
                                          1,
                                          myddlnpPRIMk,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
    // printf("kmin * exp(0 * Delta)=%f\n",kmin);
    // printf("kmin * exp((Nmax-1) * Delta)=%f\n",kmin * exp((Nmax-1) * Delta));
    // printf("kmax =%f\n",kmax);

    //GC -> here I interpolate the primordial PS, to get the transfer function. Should be ok...

    last_index = 0;
    double delta_array = 1.e-9;
    int i_kdisc = 0;

    for (i_kdisc = 0; i_kdisc < Nmax; i_kdisc++)
    {

        kdisc[i_kdisc] = kmin * exp(i_kdisc * Delta);

        if (kdisc[i_kdisc] >= exp(lnk_l[0]))
        {

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

            class_call(array_interpolate_spline(lnk_l, //GC!
                                                pnlpt->k_size,
                                                lnpPRIMk_l,
                                                myddlnpPRIMk,
                                                1,
                                                log(kdisc[i_kdisc]),
                                                &last_index,
                                                &lnpPRIMk_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);

            Pdisc[i_kdisc] = exp(lnpk_out);
            PPRIMdisc[i_kdisc] = exp(lnpPRIMk_out); //GC!

            //GC: essentially for me the trick will be to actually export the PRIMORDIAL PS from the function above...
        } // end of the condition that our k's are within the range computed by class
          //        printf("%le %le\n",kdisc[i],Pdisc[i]);
        else
        {
            Pdisc[i_kdisc] = exp(lnpk_l[0]) * pow(kdisc[i_kdisc] / exp(lnk_l[0]), ppm->n_s);         //GC ~> THIS MUST BE CHECKED!!! Ah, but this is an else... So there is no issue here... At low k I only have the primordial stuff (notice that this most likely screws up if you have running)...
            PPRIMdisc[i_kdisc] = exp(lnpPRIMk_l[0]) * pow(kdisc[i_kdisc] / exp(lnk_l[0]), ppm->n_s); //GC!
        }

        //printf("%lf %lf",kdisc[i_kdisc],PPRIMdisc[i_kdisc]); //GC ~> this doesn't give me correctly the primordial power spectrum? No, it actually does: the only issue is that at low k I get some values that are zero. This is a problem... Now I divided by A_s. But still this seems an idiotic way to do it...
        //printf("%lf %lf",kdisc[i_kdisc],Pdisc[i_kdisc]); //GC -> this correctly gives me the linear power spectrum...
        //printf("\n");
    }

    int Nside;

    /*
     To be safe we're shifting the range of kmaxes and kmins used in the analysis in order to be able to compute the AP,
     this should be increased if needed
     */

    if (pnlpt->AP_effect == AP_effect_yes)
    {
        Nside = 10;
    }
    else
    {
        Nside = 1;
    }
    double kmaxnew = kdisc[Nmax - 1 - Nside];
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
    for (index_k = 0; index_k < pnlpt->k_size - 1; index_k++)
    {
        sigmav += (lnk_l[index_k + 1] - lnk_l[index_k]) * (exp(lnk_l[index_k + 1]) * exp(lnpk_l[index_k + 1]) + exp(lnpk_l[index_k]) * exp(lnk_l[index_k])) / 2. / (6. * pow(M_PI, 2.));
    }
    // printf("%le\n",sigmav);

    int irindex = 0;

    if (pnlpt->irres == irres_yes)
    {

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
        class_alloc(logkPdiscr, Nir * sizeof(double), pnlpt->error_message);
        double logPbin2;
        double kbin2;

        double *input_realv2;
        class_alloc(input_realv2, Nirby4 * sizeof(double), pnlpt->error_message);
        double *input_imagv2;
        class_alloc(input_imagv2, Nirby4 * sizeof(double), pnlpt->error_message);
        double *output_realv2;
        class_alloc(output_realv2, Nirby4 * sizeof(double), pnlpt->error_message);
        double *output_imagv2;
        class_alloc(output_imagv2, Nirby4 * sizeof(double), pnlpt->error_message);

        last_index = 0;
        int index_ir = 0;
        for (index_ir = 0; index_ir < Nir; index_ir++)
        {

            kbin2 = kmin2 + index_ir * (kmax2 - kmin2) / (Nird - 1.);

            if (kbin2 >= exp(lnk_l[0]))
            {
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

            else
            {
                logkPdiscr[index_ir] = log(kbin2) + lnpk_l[0] + (ppm->n_s) * (log(kbin2) - lnk_l[0]);
            }

            input_realv2[2 * index_ir + 1] = pow(-1., 1. * index_ir) * logkPdiscr[index_ir];
            input_realv2[Nirby4 - 2 * index_ir - 1] = pow(-1., 1. * index_ir) * logkPdiscr[index_ir];

            input_realv2[2 * index_ir] = 0.;
            input_realv2[Nirby4 - 2 * index_ir - 2] = 0.;

            input_imagv2[2 * index_ir + 1] = 0.;
            input_imagv2[Nirby4 - 2 * index_ir - 1] = 0.;
            input_imagv2[2 * index_ir] = 0.;
            input_imagv2[Nirby4 - 2 * index_ir - 2] = 0.;
        }

        int stepsize = 1;

        FFT(input_realv2, input_imagv2, output_realv2, output_imagv2, Nirby4, stepsize);

        //        FFT(input_real_test,input_imag_test,output_real_test,output_imag_test,Nirby4,stepsize);

        double *out_ir;
        class_alloc(out_ir, Nirby4 * sizeof(double), pnlpt->error_message);

        int index_ir2 = 0;
        for (index_ir2 = 0; index_ir2 < Nir; index_ir2++)
        {
            out_ir[index_ir2] = output_realv2[Nir - index_ir2 - 1];
        }

        free(input_realv2);
        free(input_imagv2);
        free(output_realv2);
        free(output_imagv2);

        double *cmodd;
        class_alloc(cmodd, Nirover2 * sizeof(double), pnlpt->error_message);
        double *cmeven;
        class_alloc(cmeven, Nirover2 * sizeof(double), pnlpt->error_message);
        int *ivar;
        class_alloc(ivar, Nirover2 * sizeof(int), pnlpt->error_message);

        index_ir = 0;
        for (index_ir = 0; index_ir < Nirover2; index_ir++)
        {
            ivar[index_ir] = index_ir + 1;
            cmodd[index_ir] = out_ir[2 * index_ir];
            cmeven[index_ir] = out_ir[2 * index_ir + 1];
        }

        //   IR-2) Erasing the BAO bump from the odd and even Fourier harmonics and smoothy interpolating the remaining Fourier coefficients

        int Nleft = 120;
        int Nright = 240;
        // int Nright = 220;
        int Nthrow = Nright - Nleft;
        int Nnew = Nirover2 - Nthrow;
        double *cmoddnw;
        class_alloc(cmoddnw, Nnew * sizeof(double), pnlpt->error_message);
        double *cmevennw;
        class_alloc(cmevennw, Nnew * sizeof(double), pnlpt->error_message);
        double *inew;
        class_alloc(inew, Nnew * sizeof(double), pnlpt->error_message);

        index_ir = 0;
        for (index_ir = 0; index_ir < Nleft; index_ir++)
        {
            cmoddnw[index_ir] = cmodd[index_ir];
            cmevennw[index_ir] = cmeven[index_ir];
            inew[index_ir] = index_ir + 1.;
        }

        index_ir = 0;
        for (index_ir = Nleft; index_ir < Nnew; index_ir++)
        {
            cmoddnw[index_ir] = cmodd[index_ir + Nthrow];
            cmevennw[index_ir] = cmeven[index_ir + Nthrow];
            inew[index_ir] = index_ir + 1. + Nthrow;
        }

        last_index = 0;
        double *dd_cmoddnw;
        double cmodd_newval;
        class_alloc(dd_cmoddnw, sizeof(double) * Nnew, pnlpt->error_message);
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
        class_alloc(dd_cmevennw, sizeof(double) * Nnew, pnlpt->error_message);
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
        class_alloc(cmnew, Nir * sizeof(double), pnlpt->error_message);

        last_index = 0;
        last_index2 = 0;
        index_ir = 0;
        for (index_ir = 0; index_ir < Nirover2; index_ir++)
        {

            class_call(array_interpolate_spline(inew,
                                                Nnew,
                                                cmoddnw,
                                                dd_cmoddnw,
                                                1,
                                                1. * index_ir + 1.,
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
                                                1. * index_ir + 1.,
                                                &last_index2,
                                                &cmeven_newval,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);

            cmnew[index_ir * 2] = cmodd_newval;
            cmnew[index_ir * 2 + 1] = cmeven_newval;
        }

        //   IR-3) Inverse DST-II (= DST-III/(2*N)) and interpolating P_nw

        double *out_2;
        class_alloc(out_2, Nir * sizeof(double), pnlpt->error_message);
        double *input_realv3;
        class_alloc(input_realv3, Nirby4 * sizeof(double), pnlpt->error_message);
        double *input_imagv3;
        class_alloc(input_imagv3, Nirby4 * sizeof(double), pnlpt->error_message);
        double *output_realv3;
        class_alloc(output_realv3, Nirby4 * sizeof(double), pnlpt->error_message);
        double *output_imagv3;
        class_alloc(output_imagv3, Nirby4 * sizeof(double), pnlpt->error_message);

        input_realv3[0] = cmnew[Nir - 1] * 0.5;
        input_realv3[Nir] = 0.;
        input_realv3[2 * Nir] = -1. * cmnew[Nir - 1] * 0.5;
        input_realv3[3 * Nir] = 0.;

        input_imagv3[0] = 0.;
        input_imagv3[Nir] = 0.;
        input_imagv3[2 * Nir] = 0.;
        input_imagv3[3 * Nir] = 0.;

        index_ir = 0;
        for (index_ir = 1; index_ir < Nir; index_ir++)
        {
            input_realv3[index_ir] = 0.5 * cmnew[Nir - 1 - index_ir];
            input_realv3[4 * Nir - index_ir] = 0.5 * cmnew[Nir - 1 - index_ir];
            input_realv3[2 * Nir - index_ir] = -0.5 * cmnew[Nir - 1 - index_ir];
            input_realv3[2 * Nir + index_ir] = -0.5 * cmnew[Nir - 1 - index_ir];

            input_imagv3[index_ir] = 0.;
            input_imagv3[4 * Nir - index_ir] = 0.;
            input_imagv3[2 * Nir - index_ir] = 0.;
            input_imagv3[2 * Nir + index_ir] = 0.;
        }

        FFT(input_realv3, input_imagv3, output_realv3, output_imagv3, Nirby4, stepsize);

        //  double out_3[Nir];

        index_ir = 0;
        for (index_ir = 0; index_ir < Nir; index_ir++)
        {
            out_2[index_ir] = pow(-1., index_ir) * output_realv3[2 * index_ir + 1];
            //    out_3[i] = output_realv3[2*i+1]/logkPdiscr[i];
            //     printf("%e\n",out_3[i]);
        }

        free(input_realv3);
        free(input_imagv3);
        free(output_realv3);
        free(output_imagv3);

        double *Pnw_ir;
        class_alloc(Pnw_ir, Nir * sizeof(double), pnlpt->error_message);
        double *knw_ir;
        class_alloc(knw_ir, Nir * sizeof(double), pnlpt->error_message);

        index_ir = 0;
        for (index_ir = 0; index_ir < Nir; index_ir++)
        {
            knw_ir[index_ir] = kmin2 + index_ir * (kmax2 - kmin2) / (Nird - 1.);
            Pnw_ir[index_ir] = exp(out_2[index_ir] / (2 * Nird)) / (knw_ir[index_ir]);
            //       printf("%le %le\n",knw[i],Pnw[i]);
        }

        free(dd_cmoddnw);
        free(dd_cmevennw);

        double *ddPnw;
        double Pnwval;
        class_alloc(ddPnw, sizeof(double) * Nir, pnlpt->error_message);
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

        double rbao = pth->rs_d;
        // double rbao = 110./pba->h;

        // rbao = pth->rs_d;
        // printf("pth->rs_d=%lf\n",pth->rs_d);
        // printf("rbao=%lf\n",rbao);

        int Nint2 = 500;

        double *qint2;
        class_alloc(qint2, (Nint2 + 1) * sizeof(double), pnlpt->error_message);

        double *IntegrandBAO;
        class_alloc(IntegrandBAO, (Nint2 + 1) * sizeof(double), pnlpt->error_message);

        double *IntegrandBAO2;
        class_alloc(IntegrandBAO2, (Nint2 + 1) * sizeof(double), pnlpt->error_message);

        last_index = 0;

        double ks = 0.2 * pba->h;
        //double ks = 0.5 * pba->h;
        index_ir = 0;
        for (index_ir = 0; index_ir < Nint2 + 1; index_ir++)
        {

            qint2[index_ir] = kmin2 * exp(index_ir * log(ks / kmin2) / (Nint2));

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

            IntegrandBAO[index_ir] = Pnwval * (1. - 3. * sin(qint2[index_ir] * rbao) / (qint2[index_ir] * rbao) + 6. * (sin(qint2[index_ir] * rbao) / pow((qint2[index_ir] * rbao), 3.) - cos(qint2[index_ir] * rbao) / pow((qint2[index_ir] * rbao), 2.)));

            IntegrandBAO2[index_ir] = -1. * Pnwval * (3. * cos(qint2[index_ir] * rbao) * rbao * qint2[index_ir] + (-3. + pow((qint2[index_ir] * rbao), 2.)) * sin(qint2[index_ir] * rbao)) / pow((qint2[index_ir] * rbao), 3.);
        }

        SigmaBAO = 0.;
        deltaSigmaBAO = 0.;
        for (index_ir = 0; index_ir < Nint2; index_ir++)
        {
            SigmaBAO += (log(qint2[index_ir + 1]) - log(qint2[index_ir])) * (qint2[index_ir + 1] * IntegrandBAO[index_ir + 1] + qint2[index_ir] * IntegrandBAO[index_ir]) / 2. / (6. * pow(M_PI, 2.));
            deltaSigmaBAO += (log(qint2[index_ir + 1]) - log(qint2[index_ir])) * (qint2[index_ir + 1] * IntegrandBAO2[index_ir + 1] + qint2[index_ir] * IntegrandBAO2[index_ir]) / 2. / (2. * pow(M_PI, 2.));
        }
        //     double SigmaBAOh = SigmaBAO * pow(pba->h,2.);
        //    printf("Sigma_BAO(ks=0.2 h/Mpc)=%lf (Mpc/h)^2\n",SigmaBAO * pow(pba->h,2.));
        //  printf("deltaSigma_BAO(ks=0.2 h/Mpc)=%lf (Mpc/h)^2\n",deltaSigmaBAO * pow(pba->h,2.));

        // done for excersise purposes. REMOVE!
        // SigmaBAO = SigmaBAO/4.;
        // deltaSigmaBAO = deltaSigmaBAO/4.;

        // OLIVER: old hack for creating no-wiggle spectra
        //SigmaBAO *= 1000;
        //deltaSigmaBAO *= 1000;

        //   IR-5) Computing the LO IR resummed power spectrum

        double Pnwval2;
        double Pnwval_rescaled;
        double Pwval_rescaled;
        last_index = 0;
        int index_kdisc = 0;

        for (index_kdisc = 0; index_kdisc < Nmax; index_kdisc++)
        {

            if (kdisc[index_kdisc] <= kmax2 && kdisc[index_kdisc] >= kmin2 && pnlpt->alpha_rs * kdisc[index_kdisc] <= kmax2 && pnlpt->alpha_rs * kdisc[index_kdisc] >= kmin2)
            {

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
                class_call(array_interpolate_spline(knw_ir,
                                                    Nir,
                                                    Pnw_ir,
                                                    ddPnw,
                                                    1,
                                                    pnlpt->alpha_rs * kdisc[index_kdisc],
                                                    &last_index,
                                                    &Pnwval_rescaled,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                if (kdisc[index_kdisc] >= exp(lnk_l[0]))
                {

                    class_call(array_interpolate_spline(lnk_l,
                                                        pnlpt->k_size,
                                                        lnpk_l,
                                                        myddlnpk,
                                                        1,
                                                        log(pnlpt->alpha_rs * kdisc[index_kdisc]),
                                                        &last_index,
                                                        &lnpk_out,
                                                        1,
                                                        pnlpt->error_message),
                               pnlpt->error_message,
                               pnlpt->error_message);

                    Pwval_rescaled = exp(lnpk_out);
                } // end of the condition that our k's are within the range computed by class
                  //        printf("%le %le\n",kdisc[i],Pdisc[i]);
                else
                {
                    Pwval_rescaled = exp(lnpk_l[0]) * pow(pnlpt->alpha_rs * kdisc[index_kdisc] / exp(lnk_l[0]), ppm->n_s);
                }

                //nb: artificial no-wiggle not implemented for fNL with alpha_rs not equal to 1!
                Pnw[index_kdisc] = Pnwval2;
                Tnw[index_kdisc] = pow(Pnwval2 / PPRIMdisc[index_kdisc], 0.5) * (5. / 3.);
                Pw[index_kdisc] = Pwval_rescaled - Pnwval_rescaled;
                Tw[index_kdisc] = Tnw[index_kdisc] * (Pw[index_kdisc] / Pnwval2 / 2. - Pw[index_kdisc] * Pw[index_kdisc] / Pnwval2 / Pnwval2 / 8. + Pw[index_kdisc] * Pw[index_kdisc] * Pw[index_kdisc] / Pnwval2 / Pnwval2 / Pnwval2 / 16.); //GC: MAKE PLOTS WITH EH TRANSFER FUNCTION TO BE SURE OF SMOOTHNESS... This risks to not work good... No, I checked that for \Omega_b crazy large, everything is fine. But actually it would be good to check for small \Omega_b as well. Indeed I care about the relative difference... But no, small \Omega_b will make the expansion of the square root converge faster, there is no doubt about it -> I actually checked this "by eye". Things become tricky for \Omega_b = 100 what it is now, but that is crazy... Notice that here I go twice beyond leading order... See notebook!!!

                // (uncomment here if you want to test the fake NW part only)
                if (pnlpt->no_wiggle)
                {
                    Pw[index_kdisc] = 0.;
                }
                if (pnlpt->wiggle_only)
                {
                    Pnw[index_kdisc] = 0.;
                }
                Pbin[index_kdisc] = Pnw[index_kdisc] + Pw[index_kdisc] * exp(-SigmaBAO * pow(kdisc[index_kdisc], 2.));
                Tbin[index_kdisc] = Tnw[index_kdisc] + Tw[index_kdisc] * exp(-SigmaBAO * pow(kdisc[index_kdisc], 2.)); //GC!
            }

            else
            {

                Pnw[index_kdisc] = Pdisc[index_kdisc];
                Pw[index_kdisc] = 0.;
                Pbin[index_kdisc] = Pdisc[index_kdisc];

                //Tnw[index_kdisc] = pow(Pdisc[index_kdisc]/PPRIMdisc[index_kdisc],0.5)/kdisc[index_kdisc]/kdisc[index_kdisc]; //GC!
                Tnw[index_kdisc] = pow(Pdisc[index_kdisc] / PPRIMdisc[index_kdisc], 0.5) * (5. / 3.); //GC!
                Tw[index_kdisc] = 0.;                                                                 //GC!
                Tbin[index_kdisc] = Tnw[index_kdisc];                                                 //pow(Pdisc[index_kdisc]/PPRIMdisc[index_kdisc],0.5); //GC -> why compute twice, also above???
            }
            Ptree[index_kdisc] = Pnw[index_kdisc] + Pw[index_kdisc] * exp(-SigmaBAO * pow(kdisc[index_kdisc], 2.)) * (1. + SigmaBAO * pow(kdisc[index_kdisc], 2.));
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

    else
    {

        if (pnlpt->nonlinear_pt_verbose > 0)
            printf("IR resummation skipped.\n");
        int index_kd = 0;
        for (index_kd = 0; index_kd < Nmax; index_kd++)
        {
            Pbin[index_kd] = Pdisc[index_kd];
            Ptree[index_kd] = Pbin[index_kd];
            Pnw[index_kd] = Pdisc[index_kd];
            Pw[index_kd] = 0.;
            //Tbin[index_kd] = pow(Pdisc[index_kd]/PPRIMdisc[index_kd],0.5)/kdisc[index_kd]/kdisc[index_kd]; //GC!
            //GC -> I was thinking of putting the 1/k^2 here... But it becomes too complicated... I think the best thing is adding like a -2 in the bias, and it should essentially correspond to multiply by 1/k^2 the transfers without too much hassle... This seems the BEST way to do it, with a MINIMAL number of changes to the code, and a deeper use of my understanding of bias... No, cannot do this...
            Tbin[index_kd] = pow(Pdisc[index_kd] / PPRIMdisc[index_kd], 0.5) * (5. / 3.); //GC!
            Tnw[index_kd] = Tbin[index_kd];
            Tw[index_kd] = 0.;
            //GC! Here I do my thing with square roots...

            //printf("%.16f %.16f",kdisc[index_kd],Tbin[index_kd]);
            //printf("\n");
        }
    }

    class_call(array_spline_table_columns(kdisc, Nmax, Pnw, 1, dd_Pnw, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
    last_index = 0;
    double Pnw_ap_out = 0;

    class_call(array_spline_table_columns(kdisc, Nmax, Pw, 1, dd_Pw, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
    double Pw_ap_out = 0;

    //GC!

    class_call(array_spline_table_columns(kdisc, Nmax, Tnw, 1, dd_Tnw, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
    last_index = 0;
    double Tnw_ap_out = 0; //GC ~> maybe not needed, check later!!!

    class_call(array_spline_table_columns(kdisc, Nmax, Tw, 1, dd_Tw, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
    double Tw_ap_out = 0; //GC ~> maybe not needed, check later!!!

    // here we compute the FFT coefficients

    double complex *etam;
    double complex *etam_transfer; //GC!
    class_alloc(etam, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);
    class_alloc(etam_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message); //GC!

    int index_c = 0;
    double b = -0.3;
    double b_transfer = -0.8;
    //double b_transfer = -0.8000000001;
    //- 2.0; //GC! See above... I need to put a minus 2 since my transfer is just M, not M/k^2... Here -0.8 is whatever I had in the Mathematica code... I saw that essentially only for \delta^2 I need something different, -1.25...
    for (index_c = 0; index_c < Nmax + 1; index_c++)
    {
        js[index_c] = index_c - Nmaxd / 2;
        etam[index_c] = b + 2. * M_PI * _Complex_I * js[index_c] / Nmaxd / Delta;
        etam_transfer[index_c] = b_transfer + 2. * M_PI * _Complex_I * js[index_c] / Nmaxd / Delta; //GC!
        //GC -> I need etam_transfer because I have another bias!!! This is now saved and used below...
    }

    double *input_real;
    class_alloc(input_real, (Nmax) * sizeof(double), pnlpt->error_message);
    double *input_imag;
    class_alloc(input_imag, (Nmax) * sizeof(double), pnlpt->error_message);
    double *output_real;
    class_alloc(output_real, (Nmax) * sizeof(double), pnlpt->error_message);
    double *output_imag;
    class_alloc(output_imag, (Nmax) * sizeof(double), pnlpt->error_message);

    double *input_real_transfer;
    class_alloc(input_real_transfer, (Nmax) * sizeof(double), pnlpt->error_message);
    double *input_imag_transfer;
    class_alloc(input_imag_transfer, (Nmax) * sizeof(double), pnlpt->error_message);
    double *output_real_transfer;
    class_alloc(output_real_transfer, (Nmax) * sizeof(double), pnlpt->error_message);
    double *output_imag_transfer;
    class_alloc(output_imag_transfer, (Nmax) * sizeof(double), pnlpt->error_message);

    int stepsize = 1;

    int index_kd = 0;
    for (index_kd = 0; index_kd < Nmax; index_kd++)
    {
        input_real[index_kd] = Pbin[index_kd] * exp(-1. * index_kd * b * Delta);
        input_imag[index_kd] = 0.;
        //GC!
        input_real_transfer[index_kd] = Tbin[index_kd] * exp(-1. * index_kd * b_transfer * Delta);
        input_imag_transfer[index_kd] = 0.;
        //printf("%.16f",Tbin[index_kd]); //GC!
        //printf("%e",Tbin[index_kd]); //GC!
        //printf("%.20f",Tbin[index_kd]); //GC!
        //printf("%.16f",Tbin[index_kd] * exp(-1.* index_kd * b_transfer* Delta));
        //printf("%.16e",Tnw[index_kd]); //GC!
        //printf("\n");
    }

    //printf("{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}");
    //printf("\n");
    //printf("{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}{--}");

    //printf("\n");
    //printf("\n");
    //printf("\n");
    //GC!!!

    /*
     for (size_t i=0; i< N; i++){
         input_real_1[i] =in[i][0];
         input_real_2[i] =in[i][1];
     }*/

    //    FFT_real(input_real_1,input_real_2,output_real_1,output_imag_1,output_real_2,output_imag_2,N);
    //    FFT(input_real,input_imag,output_real,output_imag,N,stepsize);

    FFT(input_real, input_imag, output_real, output_imag, Nmax, stepsize);
    FFT(input_real_transfer, input_imag_transfer, output_real_transfer, output_imag_transfer, Nmax, stepsize); //GC!

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
    class_alloc(cmsym, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);
    double complex *cmsym_transfer; //GC!
    class_alloc(cmsym_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

    index_c = 0;
    for (index_c = 0; index_c < Nmax + 1; index_c++)
    {
        if (index_c < Nmax / 2)
        {
            cmsym[index_c] = cpow(kmin, -etam[index_c]) * (output_real[Nmax / 2 - index_c] - _Complex_I * output_imag[Nmax / 2 - index_c]) / Nmaxd;
            cmsym_transfer[index_c] = cpow(kmin, -etam_transfer[index_c]) * (output_real_transfer[Nmax / 2 - index_c] - _Complex_I * output_imag_transfer[Nmax / 2 - index_c]) / Nmaxd; //GC!
        }
        else
        {
            cmsym[index_c] = cpow(kmin, -etam[index_c]) * (output_real[index_c - Nmax / 2] + _Complex_I * output_imag[index_c - Nmax / 2]) / Nmaxd;
            cmsym_transfer[index_c] = cpow(kmin, -etam_transfer[index_c]) * (output_real_transfer[index_c - Nmax / 2] + _Complex_I * output_imag_transfer[index_c - Nmax / 2]) / Nmaxd; //GC!
        }
    }

    cmsym[0] = cmsym[0] / 2.;
    cmsym[Nmax] = cmsym[Nmax] / 2.;
    cmsym_transfer[0] = cmsym_transfer[0] / 2.;
    cmsym_transfer[Nmax] = cmsym_transfer[Nmax] / 2.; //GC!

    free(input_real);
    free(input_imag);
    free(output_real);
    free(output_imag);

    free(input_real_transfer);
    free(input_imag_transfer);
    free(output_real_transfer);
    free(output_imag_transfer); //GC! NO NEED TO DEALLOCATE THEM SINCE I USE THEM BELOW? Ah no, I can -> they are used to get the cmsym, so once I do NOT deallocate THEM I am fine...

    /*
     printf("Start of cnsym's\n");

     for (size_t i=0; i< Nmax+1; i++){
         printf("%le  %le\n",creal(cmsym[i]),cimag(cmsym[i]));
     }

     printf("End of cnsym's\n");
      */

    // here we input the precomputed PT matrices

    double complex nu1, nu2, nu12;

    double epsilon_for_logs_fNL = 1. * pow(10., -6.); //GC!!!
    //double large_for_logs_fNL = 1.*pow(10.,3.); //GC!!!

    //GC: ORTHOGONAL -- start
    double large_for_logs_fNL = 5. * pow(10., 4.); //GC!!!
    //GC: ORTHOGONAL -- finish
    
    double large_for_logs_matter = 5. * pow(10., 3.); //GC!
    
    double large_for_logs_big = 1000000.; //GC!
    double large_for_logs_small = 10.; //GC!

    double cutoff = 3. * pba->h; //GC -> what is this???
    double *P13;
    double *P13UV;
    double *P1loop;
    double complex *f13;
    class_alloc(P13, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(P13UV, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(P1loop, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(f13, Nmax * sizeof(complex double), pnlpt->error_message);

    double complex *f22;
    double *P22;
    double *P_CTR;
    class_alloc(f22, Nmax * sizeof(complex double), pnlpt->error_message);
    class_alloc(P22, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(P_CTR, Nmax * sizeof(double), pnlpt->error_message);

    double *P12_fNL;     //GC!
    double complex *f12; //GC!
    class_alloc(P12_fNL, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(f12, Nmax * sizeof(complex double), pnlpt->error_message);

    //GC: what is the best way to do for orthogonal? Maybe there is no need to create new variables... Well, there is no need to create new f12 stuff since y_transfer gets always overwritten, but I feel that I need a new P12_fNL_ortho...

    //GC: ORTHOGONAL -- start

    double *P12_fNL_ortho; //GC!
    class_alloc(P12_fNL_ortho, Nmax * sizeof(double), pnlpt->error_message);

    //GC: ORTHOGONAL -- finish

    int count = 0;
    //int start1=clock();
    char uplo = 'L';
    int inc = 1;
    double complex alpha = 1.;
    double complex beta = 0.;
    
    //GC - SWITCH -> later the switch index will come from outside so this will be commented...
    
    int SWITCH_index = 0;
    
    if (pnlpt->fNL_equil_ortho_switch == fNL_equil_ortho_yes)
        
    {
        
        SWITCH_index = 1;
        
    }


    double complex *x;
    double complex *x_w;
    double complex *y;
    class_alloc(x, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);
    class_alloc(x_w, (Nmax + 1) * sizeof(complex double), pnlpt->error_message); //GC: used for wiggly stuff later. Will need to allocate mine... It is freed MUCH BELOW...
    class_alloc(y, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);
    //GC: this is used to compute the AP stuff??? Unclear. Even if it was used to compute the UV piece of P13, it is not needed now... No, it seems used in the BLAS stuff below... NOTICE THAT beta is ZERO, so the use of y is useless. I assume Misha allocates it because it is necessary to use BLAS, which is fast! Anyway, I can keep this as is... HOWEVER, I need to use a x for myself, see below. SO, for now I try to allocate only x_transfer... Since the sizes are the same, I can use the same allocated y (CHECK WITH MISHA!!!)... Notice that also Misha reuses it... No, wait. y is used by zsmpv, which overwrites it (instructions -> "Before entry, the incremented array Y must contain the n element vector y. On exit, Y is overwritten by the updated vector y.")... I need one for me as well... Will fix it when I work out x_w as well...
    double complex *x_transfer; //GC!
    class_alloc(x_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

    double complex *x_w_transfer; //GC!

    class_alloc(x_w_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message); //GC!

    double complex *y_transfer; //GC!
    class_alloc(y_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

    //GC -> will add my stuff here...

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

    for (index_j = 0; index_j < Nmax; index_j++)
    {
        f13[index_j] = 0.;
        for (count = 0; count < Nmax + 1; count++)
        {
            x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
            //        printf(" %le \n", x[count]);
        }
        f13[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_oneline_complex, &inc);
        P13UV[index_j] = -61. * Pbin[index_j] * pow(kdisc[index_j], 2.) * sigmav / 105.;
        P13[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13[index_j] * Pbin[index_j]) + P13UV[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
    }
    index_j = 0;
    count = 0;
    for (index_j = 0; index_j < Nmax; index_j++)
    {
        for (count = 0; count < Nmax + 1; count++)
        {
            x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
            
            //GC - SWITCH ~> put switch here...
            
            if (SWITCH_index == 1) {
            
            x_transfer[count] = cmsym_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]); //GC!
                
                //printf("%d\n",SWITCH_index); //GC - SWITCH!
                
            }
        }

        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_complex, x, &inc, &beta, y, &inc);
        f22[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
        //printf("f22_real=%.18le f22_imag=%.18le\n",creal(f22[index_j]),cimag(f22[index_j]));
        P22[index_j] = creal(cpow(kdisc[index_j], 3.) * f22[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
        P1loop[index_j] = 0. * Ptree[index_j] + (P13[index_j] + P22[index_j]);
        //        P_CTR[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j];
        //         printf("%le %le\n",kdisc[j],P22[j]);

        if (SWITCH_index == 1) {

            //GC - SWITCH ~> put switch also here...

        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_complex, x_transfer, &inc, &beta, y_transfer, &inc);

        f12[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);
        P12_fNL[index_j] = Tbin[index_j] * //[FACTOR] *
                           creal(cpow(kdisc[index_j], 3.) * f12[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
        //kdisc[index_j] * f12[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) ); //GC! HERE CHECK CUTOFF AND OVERALL POWER! IT SHOULD BE THREE ALSO NOW, BUT CHECK! ADD HERE THE OVERALL SCALING WITH THE TRANSFER!!! ALSO, recall that I need to divide the transfer by k^2 to match my formulas... This must be done above... Cutoff is something to remove some artifact, that at least for my choices of bias seemed to never be there... The overall factors in fron ot P12 must also be added here... FINISH TO CODE THIS UP, AND THEN EXTRACT THE MATTER TRANSFER FUNCTION AS DISCUSSED WITH MISHA!!! No, I cannot play that trick if I want to extract the wiggly part using the decomposition of Misha... I should be able to extract it by taking Tbin[index_kd]? Notice that index_j goes to Nmax, from 0. And Tbin[index_kdisc] -> has also the same range. If I multiply by Tbin[...], however, I am damping its wiggles... Let's do it for now... It is correct if we are computing the 12 diagram with the IR-resummed transfer... Ahah, actually I can see it above for the P13, he needs to multiply by the PS outside, and he uses the Pbin, called as Pbin of index_j. GREAT!!! NOTICE THAT I STILL HAVE TO DIVIDE BY k^2... I DID IT ALREADY HERE (removing the k^3 and substituting it with a simple k), but I must do it also above where I compute the coefficients. Here k should be already in units of h/Mpc... Nice, now my division is done at the bias level... Notice that this works also for local NG... No, I need to use the biases I have since I am sure they work...

        //GC -> TEST!!!
        //P12_fNL[index_j] = creal(cpow(kdisc[index_j], 3.) * f12[index_j] * exp(-pow(kdisc[index_j]/cutoff, 6.)) );

        //printf("%.16f %.16f %.16f",kdisc[index_j],Tbin[index_j],P12_fNL[index_j]);
        //printf("%e %e %e",kdisc[index_j],Tbin[index_j],P12_fNL[index_j]);
        //printf("%.20f %.20f %.20f",kdisc[index_j],Tbin[index_j],P12_fNL[index_j]);
        //GC -> so now I have the correct points to compute the integral, and the correct sampling points. NO INTERPOLATION...
        //printf("%.16e %.16e %.16e",kdisc[index_j],Tbin[index_j],P12_fNL[index_j]);
        //printf("\n");

        //GC!

        //GC: ORTHOGONAL -- start

        zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);
        f12[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);
        P12_fNL_ortho[index_j] = Tbin[index_j] * //[FACTOR] *
                                 creal(cpow(kdisc[index_j], 3.) * f12[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));

        //printf("%.16e %.16e %.16e",kdisc[index_j],Tbin[index_j],P12_fNL_ortho[index_j]);
        //printf("\n");
        
        //GC - SWITCH ~> put switch up to here...

        }
            
        //GC: ORTHOGONAL -- finish
    }
    //} //End of RSD only condition

    for (index_j = 0; index_j < Nmax; index_j++)
    {
        P_CTR[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j];
    }

    double *ddpk_nl;
    class_alloc(ddpk_nl, sizeof(double) * Nmax, pnlpt->error_message);
    double *ddpk_nl_fNL;
    class_alloc(ddpk_nl_fNL, sizeof(double) * Nmax, pnlpt->error_message); //GC!

    //GC: ORTHOGONAL -- start

    double *ddpk_nl_fNL_ortho;
    class_alloc(ddpk_nl_fNL_ortho, sizeof(double) * Nmax, pnlpt->error_message); //GC!

    //GC: ORTHOGONAL -- finish

    double *ddpk_CTR;
    class_alloc(ddpk_CTR, sizeof(double) * Nmax, pnlpt->error_message);
    double *ddpk_Tree;
    class_alloc(ddpk_Tree, sizeof(double) * Nmax, pnlpt->error_message);
    
    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P1loop,
                                          1,
                                          ddpk_nl,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    
    //GC - SWITCH ~> put switch here... This indeed could be a problem. Now this function wants to eat P12_fNL. But this is utter garbage in general... How am I sure that it does not give segmentation faults? Will first do a run with matter only...
    
    //GC: seems to work (no segmentation fault) -> now I am feeding junk, and I print as a check below...


    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P12_fNL,
                                          1,
                                          ddpk_nl_fNL,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    //GC!

    //GC: ORTHOGONAL -- start

    class_call(array_spline_table_columns(kdisc,
                                          Nmax,
                                          P12_fNL_ortho,
                                          1,
                                          ddpk_nl_fNL_ortho,
                                          _SPLINE_NATURAL_,
                                          pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    //GC: ORTHOGONAL -- finish

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
    double pk_nl_fNL_out; //GC!

    //GC: ORTHOGONAL -- start

    double pk_nl_fNL_out_ortho; //GC!

    //GC: ORTHOGONAL -- finish

    double pk_CTR_out;
    double pk_Tree_out;

    last_index = 0;
    last_index2 = 0;
    int last_index3 = 0;
    int last_index4 = 0; //GC!

    //GC: ORTHOGONAL -- start

    int last_index4_ortho = 0; //GC!

    //GC: ORTHOGONAL -- finish

    //printf("\n"); //GC!
    //printf("[*][*][*][*][*][*][*][*][*][*][*][*][*]\n"); //GC!
    //printf("[*][*][*][*][*][*][*][*][*][*][*][*][*]\n"); //GC!
    //printf("\n"); //GC!

    for (index_k = 0; index_k < pnlpt->k_size; index_k++)
    {

        if (pnlpt->k[index_k] >= kmin && pnlpt->k[index_k] <= kmax)
        {

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
                                                P12_fNL,
                                                ddpk_nl_fNL,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index4,
                                                &pk_nl_fNL_out,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            //GC!

            //GC: ORTHOGONAL -- start

            class_call(array_interpolate_spline(kdisc,
                                                Nmax,
                                                P12_fNL_ortho,
                                                ddpk_nl_fNL_ortho,
                                                1,
                                                pnlpt->k[index_k],
                                                &last_index4_ortho,
                                                &pk_nl_fNL_out_ortho,
                                                1,
                                                pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);

            //GC: ORTHOGONAL -- finish

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

            pk_nl[index_k] = pk_nl_out + large_for_logs_matter; //GC -> because log can become negative...
            //pk_nl_fNL[index_k] = pk_nl_fNL_out + 5000.; //GC -> because log can become negative... WILL NEED TO CHECK NUMBERS, if 5000 is enough!!! Would something change if I had chosen 50000? Would it have caused errors?
            pk_nl_fNL[index_k] = pk_nl_fNL_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;

            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_nl_fNL[index_k]-1.*large_for_logs_fNL); //GC! TEST!!!
            
            //GC - SWITCH -> SEEMS TO WORK FINE!
            
            //GC: ORTHOGONAL -- start

            pk_nl_fNL_ortho[index_k] = pk_nl_fNL_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;

            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_nl_fNL_ortho[index_k]-1.*large_for_logs_fNL); //GC! TEST!!!
            
            //GC - SWITCH -> SEEMS TO WORK FINE -> indeed it prints random numbers, that before I removed the large_for_logs_fNL were overshadowed by it -> -1.4036503372991333e-02, 6.6630490593878808e-02, -1.2951487905045988e-08, 5.1819733362425955e-20, etc. Notice that even if I get a number larger than large_for_logs_fNL, I KNOW that it does not give segmentation fault. It just gives some NaNs from spline stuff, as we know...
            
            //GC: ORTHOGONAL -- finish
	    if(pk_CTR_out<=0) pk_CTR_out=1e-16;
	    if(pk_Tree_out<=0) pk_Tree_out=1e-16;
            pk_CTR[index_k] = pk_CTR_out;
            pk_Tree[index_k] = pk_Tree_out;
        }

        else
        {
            //pk_nl[index_k] = exp(lnpk_l[index_k]);
            if (pnlpt->k[index_k] < kmin)
            {
                pk_nl[index_k] = large_for_logs_matter - 61. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav / 105.;
                //pk_nl_fNL[index_k] = 5000.;
                pk_nl_fNL[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!

                //GC: ORTHOGONAL -- start

                pk_nl_fNL_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!

                //GC: ORTHOGONAL -- finish
            }
            if (pnlpt->k[index_k] > kmax)
            {
                pk_nl[index_k] = large_for_logs_matter;
                //pk_nl_fNL[index_k] = 5000.;
                pk_nl_fNL[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!

                //GC: ORTHOGONAL -- start

                pk_nl_fNL_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!

                //GC: ORTHOGONAL -- finish
            }
            pk_CTR[index_k] = exp(lnpk_l[index_k] + 2. * lnk_l[index_k]);
            pk_Tree[index_k] = exp(lnpk_l[index_k]);
        }
        //  printf("%i %f %f\n",index_k, pk_Tree[index_k],pk_nl[index_k]-5000.);
        //printf("%.16f %.16f\n",pnlpt->k[index_k],pk_nl_fNL[index_k]-0.*large_for_logs_fNL); //GC!
        //printf("%e %e\n",pnlpt->k[index_k],pk_nl_fNL[index_k]-0.*large_for_logs_fNL); //GC!
        //printf("%.20f %.20f\n",pnlpt->k[index_k],pk_nl_fNL[index_k]-0.*large_for_logs_fNL); //GC!

        //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_nl_fNL[index_k]-0.*large_for_logs_fNL); //GC!
        //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_nl_fNL_ortho[index_k]-1.*large_for_logs_fNL); //GC!
    }

    if (pnlpt->rsd == rsd_yes)
    { //GC: screwed up a bracket. Should end at the line "// end of RSD conditional expression"...
        //GC -> need to go step-by-step: it seems that here the calculation is done for RSDs on matter only... I think that 0, 2, 4 means the power of \mu, while vv, vd and dd means the power of b1... So actually this is even different -> I add here all the contributions that come from a LINEARLY BIASED TRACER IN REDSHIFT SPACE!!! Time to get these matrices! Notice that I should account for FoG in Z1? Notice also that I was missing the factors of f when combining in the notebook? No... Anyway, no FoG here yet... But the FoG are STOCHASTIC, so why do we change Z1??? Boh... Good idea to not add them yet... For the counterterms see below...

        //GC -> I have, before any change to counterterms, vv_0, vd_0, dd_0; vv_0, vd_2, dd_2; vv_4, vd_4, and that's it...

        if (pnlpt->nonlinear_pt_verbose > 0)
            printf("Computing RSD...\n");
        if (pnlpt->nonlinear_pt_verbose > 0)
            //printf("Logarithmic growth factor f=%f\n",f); //GC!
            printf("Logarithmic growth factor f=%.16e\n", f);

        double *Ptree_0_vv;
        class_alloc(Ptree_0_vv, Nmax * sizeof(double), pnlpt->error_message);
        double *Ptree_0_vd;
        class_alloc(Ptree_0_vd, Nmax * sizeof(double), pnlpt->error_message);
        double *Ptree_0_dd;
        class_alloc(Ptree_0_dd, Nmax * sizeof(double), pnlpt->error_message);
        double *Ptree_2_vv;
        class_alloc(Ptree_2_vv, Nmax * sizeof(double), pnlpt->error_message);
        double *Ptree_2_vd;
        class_alloc(Ptree_2_vd, Nmax * sizeof(double), pnlpt->error_message);
        double *Ptree_4_vv;
        class_alloc(Ptree_4_vv, Nmax * sizeof(double), pnlpt->error_message);

        double *P13_0_vv;
        double *P13UV_0_vv;
        double *P1loop_0_vv;
        class_alloc(P13_0_vv, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P13UV_0_vv, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P1loop_0_vv, Nmax * sizeof(double), pnlpt->error_message);
        double *P22_0_vv;
        class_alloc(P22_0_vv, Nmax * sizeof(double), pnlpt->error_message);

        double *P12_0_vv;
        class_alloc(P12_0_vv, Nmax * sizeof(double), pnlpt->error_message);
        //GC -> this contains nothing. There isn't anything that has no b1, and no \mu... Same for 0_vd... No, there are multipoles... Some comments below are wrong/useless...

        //GC: ORTHOGONAL -- start

        double *P12_0_vv_ortho;
        class_alloc(P12_0_vv_ortho, Nmax * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        double *P13_0_vd;
        double *P13UV_0_vd;
        double *P1loop_0_vd;
        class_alloc(P13_0_vd, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P13UV_0_vd, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P1loop_0_vd, Nmax * sizeof(double), pnlpt->error_message);
        double *P22_0_vd;
        class_alloc(P22_0_vd, Nmax * sizeof(double), pnlpt->error_message);

        double *P12_0_vd;
        class_alloc(P12_0_vd, Nmax * sizeof(double), pnlpt->error_message);
        //GC!

        //GC: ORTHOGONAL -- start

        double *P12_0_vd_ortho;
        class_alloc(P12_0_vd_ortho, Nmax * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        double *P13_0_dd;
        double *P13UV_0_dd;
        double *P1loop_0_dd;
        class_alloc(P13_0_dd, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P13UV_0_dd, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P1loop_0_dd, Nmax * sizeof(double), pnlpt->error_message);
        double *P22_0_dd;
        class_alloc(P22_0_dd, Nmax * sizeof(double), pnlpt->error_message);

        double *P12_0_dd;
        class_alloc(P12_0_dd, Nmax * sizeof(double), pnlpt->error_message);
        //GC -> this contains the term from b1^2 F2...

        //GC: ORTHOGONAL -- start

        double *P12_0_dd_ortho;
        class_alloc(P12_0_dd_ortho, Nmax * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        double *P13_2_vv;
        double *P13UV_2_vv;
        double *P1loop_2_vv;
        class_alloc(P13_2_vv, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P13UV_2_vv, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P1loop_2_vv, Nmax * sizeof(double), pnlpt->error_message);
        double *P22_2_vv;
        class_alloc(P22_2_vv, Nmax * sizeof(double), pnlpt->error_message);

        double *P12_2_vv;
        class_alloc(P12_2_vv, Nmax * sizeof(double), pnlpt->error_message);
        //GC -> this is also zero...

        //GC: ORTHOGONAL -- start

        double *P12_2_vv_ortho;
        class_alloc(P12_2_vv_ortho, Nmax * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        double *P13_2_vd;
        double *P13UV_2_vd;
        double *P1loop_2_vd;
        class_alloc(P13_2_vd, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P13UV_2_vd, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P1loop_2_vd, Nmax * sizeof(double), pnlpt->error_message);
        double *P22_2_vd;
        class_alloc(P22_2_vd, Nmax * sizeof(double), pnlpt->error_message);

        double *P12_2_vd;
        class_alloc(P12_2_vd, Nmax * sizeof(double), pnlpt->error_message);
        //GC -> this contains a piece b1*f*\mu^2*F2, a piece b1*f*\mu^2*G2, and a piece that comes from the shift part in Z2 with f^2 that generates a \mu^2, and then multiplied by the b1 in Z1...

        //GC: ORTHOGONAL -- start

        double *P12_2_vd_ortho;
        class_alloc(P12_2_vd_ortho, Nmax * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        double *P13_2_dd;
        double *P13UV_2_dd;
        double *P1loop_2_dd;
        class_alloc(P13_2_dd, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P13UV_2_dd, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P1loop_2_dd, Nmax * sizeof(double), pnlpt->error_message);
        double *P22_2_dd;
        class_alloc(P22_2_dd, Nmax * sizeof(double), pnlpt->error_message);

        double *P12_2_dd;
        class_alloc(P12_2_dd, Nmax * sizeof(double), pnlpt->error_message);
        //GC -> this only the piece that comes from the shift part in Z2 with f, that generates ONLY a \mu^2 and has an overall b1, after being then multiplied by the b1 in Z1...

        //GC: ORTHOGONAL -- start

        double *P12_2_dd_ortho;
        class_alloc(P12_2_dd_ortho, Nmax * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        double *P13_4_vv;
        double *P13UV_4_vv;
        double *P1loop_4_vv;
        class_alloc(P13_4_vv, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P13UV_4_vv, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P1loop_4_vv, Nmax * sizeof(double), pnlpt->error_message);
        double *P22_4_vv;
        class_alloc(P22_4_vv, Nmax * sizeof(double), pnlpt->error_message);

        double *P12_4_vv;
        class_alloc(P12_4_vv, Nmax * sizeof(double), pnlpt->error_message);
        //GC -> this one comes from the G2 piece, with the f\mu^2 from Z1. Then, the pain comes. I have that the shift part that has b1 cannot contribute to vv. Then, I need the part without b1, i.e. the f^2\mu^2 subset of the shift piece. Then, I need to multiply this by a f\mu^2 in the Z1. Misha here stops at \mu^4, so no other pieces are needed...

        //GC: ORTHOGONAL -- start

        double *P12_4_vv_ortho;
        class_alloc(P12_4_vv_ortho, Nmax * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        double *P13_4_vd;
        double *P13UV_4_vd;
        double *P1loop_4_vd;
        class_alloc(P13_4_vd, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P13UV_4_vd, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(P1loop_4_vd, Nmax * sizeof(double), pnlpt->error_message);
        double *P22_4_vd;
        class_alloc(P22_4_vd, Nmax * sizeof(double), pnlpt->error_message);

        double *P12_4_vd;
        class_alloc(P12_4_vd, Nmax * sizeof(double), pnlpt->error_message);
        //GC -> this requires the b1 f\mu^2 piece in the shift, multiplied by the f\mu^2 in Z1. Then, in the f^2 part of the shift I have a \mu^4 part that I multiply by the b1 in Z1, and then nothing else...

        //GC: ORTHOGONAL -- start

        double *P12_4_vd_ortho;
        class_alloc(P12_4_vd_ortho, Nmax * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        double *P22_4_dd;
        class_alloc(P22_4_dd, Nmax * sizeof(double), pnlpt->error_message);
        double *P1loop_4_dd;
        class_alloc(P1loop_4_dd, Nmax * sizeof(double), pnlpt->error_message);

        double *P12_4_dd;
        class_alloc(P12_4_dd, Nmax * sizeof(double), pnlpt->error_message);
        //GC -> this should be actually zero... Again, notice that Misha here stops at \mu^4... This is zero now... Notice that actually I need to allocate it for the AP below... All the others were essentially already allocated because I have a contribution even before AP...

        //GC: ORTHOGONAL -- start

        double *P12_4_dd_ortho;
        class_alloc(P12_4_dd_ortho, Nmax * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        double *P_CTR_0;
        class_alloc(P_CTR_0, Nmax * sizeof(double), pnlpt->error_message);
        double *P_CTR_2;
        class_alloc(P_CTR_2, Nmax * sizeof(double), pnlpt->error_message);
        double *P_CTR_4;
        class_alloc(P_CTR_4, Nmax * sizeof(double), pnlpt->error_message);

        //    double help1; //GC!
        //    double help2; //GC!
        //    double help3; //GC!

        //    double complex log1; //GC!
        //    double complex log2; //GC!
        //    double complex log11; //GC!
        //    double complex log12; //GC!
        //    double complex log13; //GC!
        //    double complex log21; //GC!
        //    double complex log22; //GC!
        //    double complex log23; //GC!

        //    double twodouble = 2.0; //GC!

        //  if (pnlpt->irres == irres_no){

        //  if (pnlpt->irres == irres_yes){
        if (irindex == 0)
        {

            //GC: here is redshift space (rsd == rsd_yes), without IR resummation -> I am in RSD, but with irindex == 0 which mean Tbin = Tnw... Ok, so never double counting, and I do the resummation with Tnw * (1 + ratio) * \int Tnw Tnw + 2 Tnw * (\int Tw Tnw) * exp, where ratio = (Tw/Tnw)*exp... So in this sense, Tbin is not really used ever I think, i.e. it is not used when I actually have the wiggly things resummation...
            //GC -> notice that also this bracket is screwed: ends at the line "// end of IR resummation expression"...

            //  printf("Computing RSD without IR resummation...\n");

            // Computing P_{vv} contribution

            //GC -> here it means no IR resummation in redshift space I think... The coefficients for the FFTlog are the same used for the matter above... Also, notice that maybe then 0, 2 and 4 did not correspond to the powers of \mu, but already to monopole, quadrupole and hexadecapole... I need to collect my matrices correctly on a notebook/note... Notice that for sure 0 is monopole, while mu0 is \mu^0... Confirmed by the fact that f appears in the matrices... Also, notice that f is extracted above... OK -> so here I have already all the matrices. There are a few entries that are zero in my case. Notice that I do not have FoG in Z1, and also I should have the counterterms as scale-dependent bias (which comes from the high-k limit of the 12 integral) and some correction to the noise, which comes from the low-k limit of the integral. This correction will not be exactly local... Notice -> that I was incorrect: of course, the counterterms come from the UV limit of the integral only, i.e. k -> 0...

            count = 0;
            index_l = 0;
            index_i = 0;
            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam[index_i];
                    nu2 = -0.5 * etam[index_l];
                    nu12 = nu1 + nu2;

                    pnlpt->M22_oneline_0_vv_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * (f * f * (14. * f * f * (24. - 8. * nu2 - 15. * nu2 * nu2 + 5. * nu2 * nu2 * nu2 + 5. * nu1 * nu1 * nu1 * (1. + 7. * nu2) + 5. * nu1 * nu1 * (-3. - 10. * nu2 + 14. * nu2 * nu2) + nu1 * (-8. - 24. * nu2 - 50. * nu2 * nu2 + 35. * nu2 * nu2 * nu2)) + 18. * f * (36. - 8. * nu2 + 70. * nu1 * nu1 * nu1 * nu2 - 23. * nu2 * nu2 + nu1 * nu1 * (-23. - 94. * nu2 + 140. * nu2 * nu2) + nu1 * (-8. - 42. * nu2 - 94. * nu2 * nu2 + 70. * nu2 * nu2 * nu2)) + 9. * (50. - 9. * nu2 + 98. * nu1 * nu1 * nu1 * nu2 - 35. * nu2 * nu2 + 7. * nu1 * nu1 * (-5. - 18. * nu2 + 28. * nu2 * nu2) + nu1 * (-9. - 66. * nu2 - 126. * nu2 * nu2 + 98. * nu2 * nu2 * nu2)))) / 8820.;

                    //           printf("%d %d %.16f",index_i,index_l,-0.5*etam[index_i]*_Complex_I); //GC -> good, it is because they are complex...
                    //            printf("\n");

                    if (SWITCH_index == 1) {
                        
                    //GC - SWITCH!

                    
                    nu1 = -0.5 * etam_transfer[index_i];
                    nu2 = -0.5 * etam_transfer[index_l]; //GC!

                    //GC! ALL THESE THINGS MUST BE ALLOCATED ABOVE!!! Notice that it is overkill to have etam_transfer. It is actually the same as the others? NO, THE BIAS CHANGES...

                    //pnlpt->M12_oneline_0_vv_complex[count] = (pnlpt->M12_oneline_complex[count])*(f*f*((-30.*(4. + 3.*f) + (145. + 99.*f)*nu2 + 64.*(7. + 5.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + (948. + 716.*f)*nu2*nu2 - 4.*(473. + 347.*f)*nu2*nu2*nu2 + 64.*(18. + 13.*f)*nu2*nu2*nu2*nu2 - 32.*(7. + 5.*f)*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(3.*(7. + 5.*f) + (438. + 322.*f)*nu2 - 16.*(51. + 37.*f)*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(-69. - 51.*f - 4.*(135. + 103.*f)*nu2 + 4.*(413. + 307.*f)*nu2*nu2 - 8.*(169. + 123.*f)*nu2*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2*nu2) + nu1*(321. + 243.*f + 8.*(99. + 79.*f)*nu2 - 8.*(605. + 459.*f)*nu2*nu2 + 8.*(751. + 557.*f)*nu2*nu2*nu2 - 256.*(11. + 8.*f)*nu2*nu2*nu2*nu2 + 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (32.*(7. + 5.*f)*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(18. + 13.*f - 4.*(11. + 8.*f)*nu2 + 3.*(7. + 5.*f)*nu2*nu2) + 3.*(-10.*(4. + 3.*f) + (107. + 81.*f)*nu2 - 4.*(23. + 17.*f)*nu2*nu2 + 4.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-473. - 347.*f + 2.*(751. + 557.*f)*nu2 - 8.*(169. + 123.*f)*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2) + nu1*(145. + 99.*f + 8.*(99. + 79.*f)*nu2 - 16.*(135. + 103.*f)*nu2*nu2 + 8.*(219. + 161.*f)*nu2*nu2*nu2 - 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(237. + 179.*f - 2.*(605. + 459.*f)*nu2 + 4.*(413. + 307.*f)*nu2*nu2 - 16.*(51. + 37.*f)*nu2*nu2*nu2 + 16.*(7. + 5.*f)*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (64.*(7. + 5.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 - 3.*(-10.*(4. + 3.*f) + (107. + 81.*f)*nu2 - 4.*(23. + 17.*f)*nu2*nu2 + 4.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-3.*(7. + 5.*f) + (466. + 342.*f)*nu2 - 8.*(81. + 59.*f)*nu2*nu2 + 32.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(69. + 51.*f - 4.*(179. + 135.*f)*nu2 + 4.*(289. + 215.*f)*nu2*nu2 - 8.*(81. + 59.*f)*nu2*nu2*nu2 + 16.*(7. + 5.*f)*nu2*nu2*nu2*nu2) - 1.*nu1*(321. + 243.*f - 16.*(116. + 89.*f)*nu2 + 16.*(179. + 135.*f)*nu2*nu2 - 8.*(233. + 171.*f)*nu2*nu2*nu2 + 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))/(5.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI))); //GC!

                    //pnlpt->M12_oneline_0_vv_complex[count] = csin(2.*(nu1 + nu2)*M_PI); //(pnlpt->M12_oneline_complex[count]);

                    //pnlpt->M12_oneline_0_vv_complex[count] = (pnlpt->M12_oneline_complex[count]);
                    //GC -> if I use only this, it is fine... And it is a problem ALREADY if I do csin(2.*(nu1 + nu2)*M_PI)...

                    //pnlpt->M12_oneline_0_vv_complex[count] = (pnlpt->M12_oneline_complex[count])*(nu1+nu2);
                    //GC -> this is perfectly fine...

                    //printf("%.16f",csin(1.3*M_PI));
                    //printf("%.16f",csin(1.+_Complex_I*M_PI/1.3));
                    //printf("%.16f",nu1); //GC -> this prints always the same number? But also for the "non-transfer" one: I cannot believe this is the issue... But it must be related to it, because the csin(1.+_Complex_I*M_PI/1.3) works... See above... They are comparable with the ones for the old code because the only thing that changes is the REAL PART, the bias... The imaginary part is set by kmin and kmax, and it is always the same... Great, this is why before it was spitting out always a single number...
                    //printf("%d %d %.16f",index_i,index_l,nu1*_Complex_I); //GC -> good, it is because they are complex...
                    //printf("\n");
                    //printf("%.16f",csin(2.*47.3279497208556066*(1.+_Complex_I)*M_PI));
                    //GC -> this is the issue... It cannot compute these numbers...
                    //printf("%d %d %e",index_i,index_l,pnlpt->M12_oneline_0_vv_complex[count]);

                    //help=/*(f*f*((-30.*(4. + 3.*f) + (145. + 99.*f)*nu2 + 64.*(7. + 5.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + (948. + 716.*f)*nu2*nu2 - 4.*(473. + 347.*f)*nu2*nu2*nu2 + 64.*(18. + 13.*f)*nu2*nu2*nu2*nu2 - 32.*(7. + 5.*f)*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(3.*(7. + 5.*f) + (438. + 322.*f)*nu2 - 16.*(51. + 37.*f)*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(-69. - 51.*f - 4.*(135. + 103.*f)*nu2 + 4.*(413. + 307.*f)*nu2*nu2 - 8.*(169. + 123.*f)*nu2*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2*nu2) + nu1*(321. + 243.*f + 8.*(99. + 79.*f)*nu2 - 8.*(605. + 459.*f)*nu2*nu2 + 8.*(751. + 557.*f)*nu2*nu2*nu2 - 256.*(11. + 8.*f)*nu2*nu2*nu2*nu2 + 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (32.*(7. + 5.*f)*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(18. + 13.*f - 4.*(11. + 8.*f)*nu2 + 3.*(7. + 5.*f)*nu2*nu2) + 3.*(-10.*(4. + 3.*f) + (107. + 81.*f)*nu2 - 4.*(23. + 17.*f)*nu2*nu2 + 4.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-473. - 347.*f + 2.*(751. + 557.*f)*nu2 - 8.*(169. + 123.*f)*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2) + nu1*(145. + 99.*f + 8.*(99. + 79.*f)*nu2 - 16.*(135. + 103.*f)*nu2*nu2 + 8.*(219. + 161.*f)*nu2*nu2*nu2 - 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(237. + 179.*f - 2.*(605. + 459.*f)*nu2 + 4.*(413. + 307.*f)*nu2*nu2 - 16.*(51. + 37.*f)*nu2*nu2*nu2 + 16.*(7. + 5.*f)*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (64.*(7. + 5.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 - 3.*(-10.*(4. + 3.*f) + (107. + 81.*f)*nu2 - 4.*(23. + 17.*f)*nu2*nu2 + 4.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-3.*(7. + 5.*f) + (466. + 342.*f)*nu2 - 8.*(81. + 59.*f)*nu2*nu2 + 32.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(69. + 51.*f - 4.*(179. + 135.*f)*nu2 + 4.*(289. + 215.*f)*nu2*nu2 - 8.*(81. + 59.*f)*nu2*nu2*nu2 + 16.*(7. + 5.*f)*nu2*nu2*nu2*nu2) - 1.*nu1*(321. + 243.*f - 16.*(116. + 89.*f)*nu2 + 16.*(179. + 135.*f)*nu2*nu2 - 8.*(233. + 171.*f)*nu2*nu2*nu2 + 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))*/1./(5.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))/**csin(2.*(nu1 + nu2)*M_PI)*/));

                    //help=nu1*_Complex_I;

                    //help = 1./csin(2.*(-0.4)*M_PI);

                    //help = -2.*nu2*M_PI*_Complex_I; //1./csin(2.*nu2*M_PI);

                    //help = csin(2.*(nu1+0.0001)*M_PI)/csin(2.*nu2*M_PI); *****!!! IT WORKS...

                    //help1=(224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI);

                    //help2=(f*f*((-30.*(4. + 3.*f) + (145. + 99.*f)*nu2 + 64.*(7. + 5.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + (948. + 716.*f)*nu2*nu2 - 4.*(473. + 347.*f)*nu2*nu2*nu2 + 64.*(18. + 13.*f)*nu2*nu2*nu2*nu2 - 32.*(7. + 5.*f)*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(3.*(7. + 5.*f) + (438. + 322.*f)*nu2 - 16.*(51. + 37.*f)*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(-69. - 51.*f - 4.*(135. + 103.*f)*nu2 + 4.*(413. + 307.*f)*nu2*nu2 - 8.*(169. + 123.*f)*nu2*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2*nu2) + nu1*(321. + 243.*f + 8.*(99. + 79.*f)*nu2 - 8.*(605. + 459.*f)*nu2*nu2 + 8.*(751. + 557.*f)*nu2*nu2*nu2 - 256.*(11. + 8.*f)*nu2*nu2*nu2*nu2 + 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (32.*(7. + 5.*f)*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(18. + 13.*f - 4.*(11. + 8.*f)*nu2 + 3.*(7. + 5.*f)*nu2*nu2) + 3.*(-10.*(4. + 3.*f) + (107. + 81.*f)*nu2 - 4.*(23. + 17.*f)*nu2*nu2 + 4.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-473. - 347.*f + 2.*(751. + 557.*f)*nu2 - 8.*(169. + 123.*f)*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2) + nu1*(145. + 99.*f + 8.*(99. + 79.*f)*nu2 - 16.*(135. + 103.*f)*nu2*nu2 + 8.*(219. + 161.*f)*nu2*nu2*nu2 - 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(237. + 179.*f - 2.*(605. + 459.*f)*nu2 + 4.*(413. + 307.*f)*nu2*nu2 - 16.*(51. + 37.*f)*nu2*nu2*nu2 + 16.*(7. + 5.*f)*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (64.*(7. + 5.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 - 3.*(-10.*(4. + 3.*f) + (107. + 81.*f)*nu2 - 4.*(23. + 17.*f)*nu2*nu2 + 4.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-3.*(7. + 5.*f) + (466. + 342.*f)*nu2 - 8.*(81. + 59.*f)*nu2*nu2 + 32.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(69. + 51.*f - 4.*(179. + 135.*f)*nu2 + 4.*(289. + 215.*f)*nu2*nu2 - 8.*(81. + 59.*f)*nu2*nu2*nu2 + 16.*(7. + 5.*f)*nu2*nu2*nu2*nu2) - 1.*nu1*(321. + 243.*f - 16.*(116. + 89.*f)*nu2 + 16.*(179. + 135.*f)*nu2*nu2 - 8.*(233. + 171.*f)*nu2*nu2*nu2 + 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2))//*csin(2.*(nu1 + nu2)*M_PI)
                    //));

                    //help3 = (5.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))//*csin(2.*(nu1 + nu2)*M_PI)
                    //));

                    /*
            if (index_l < 10 && index_i < 10) {
                //printf("%e",help2);
                printf("%e",help2);
                printf("\n");
                printf("%e",csin(2.*(nu1 + nu2)*M_PI));
                printf("\n");
                printf("*****************************");
                printf("\n");
                //printf("%e",(-5.641947e+164)/(-4.056806e+301)); //GC -> works...
                printf("%e",(help2)/(csin(2.*(nu1 + nu2)*M_PI)));
                printf("\n");
                printf("#############################");
                printf("\n");
            }
             */

                    /*
            
            if (index_l < 10 && index_i < 10) {
                //printf("%e",help2);
                printf("%e",help3);
                printf("\n");
                printf("%e",csin(2.*(nu1 + nu2)*M_PI));
                printf("\n");
                printf("*****************************");
                printf("\n");
                //printf("%e",(-2.819291e+165)/(-4.056806e+301)); //GC -> works...
                //printf("\n");
                //printf("+++++++++++++++++++++++++++++");
                //printf("\n");
                //printf("%e",(help3)/(csin(2.*(nu1 + nu2)*M_PI)));
                //printf("%e",(-2.819291e+165)/(csin(2.*(nu1 + nu2)*M_PI)));
                printf("%e",(help3)/(-4.056806e+301)); //GC -> this works...
                printf("\n");
                printf("#############################");
                printf("\n");
            }
             
             */

                    //help1=csin(twodouble*nu1*M_PI);//1.;//csin(2.*nu1*M_PI);
                    //help2=csin(twodouble*nu2*M_PI);//1.;//csin(2.*nu2*M_PI);
                    //help3=csin(twodouble*(nu1 + nu2)*M_PI);//1.;//csin(2.*(nu1 + nu2)*M_PI);

                    //pnlpt->M12_oneline_0_vv_complex[count] = (pnlpt->M12_oneline_complex[count])*(f*f*((-30.*(4. + 3.*f) + (145. + 99.*f)*nu2 + 64.*(7. + 5.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + (948. + 716.*f)*nu2*nu2 - 4.*(473. + 347.*f)*nu2*nu2*nu2 + 64.*(18. + 13.*f)*nu2*nu2*nu2*nu2 - 32.*(7. + 5.*f)*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(3.*(7. + 5.*f) + (438. + 322.*f)*nu2 - 16.*(51. + 37.*f)*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(-69. - 51.*f - 4.*(135. + 103.*f)*nu2 + 4.*(413. + 307.*f)*nu2*nu2 - 8.*(169. + 123.*f)*nu2*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2*nu2) + nu1*(321. + 243.*f + 8.*(99. + 79.*f)*nu2 - 8.*(605. + 459.*f)*nu2*nu2 + 8.*(751. + 557.*f)*nu2*nu2*nu2 - 256.*(11. + 8.*f)*nu2*nu2*nu2*nu2 + 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2*nu2))*help1 + (32.*(7. + 5.*f)*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(18. + 13.*f - 4.*(11. + 8.*f)*nu2 + 3.*(7. + 5.*f)*nu2*nu2) + 3.*(-10.*(4. + 3.*f) + (107. + 81.*f)*nu2 - 4.*(23. + 17.*f)*nu2*nu2 + 4.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-473. - 347.*f + 2.*(751. + 557.*f)*nu2 - 8.*(169. + 123.*f)*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2) + nu1*(145. + 99.*f + 8.*(99. + 79.*f)*nu2 - 16.*(135. + 103.*f)*nu2*nu2 + 8.*(219. + 161.*f)*nu2*nu2*nu2 - 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(237. + 179.*f - 2.*(605. + 459.*f)*nu2 + 4.*(413. + 307.*f)*nu2*nu2 - 16.*(51. + 37.*f)*nu2*nu2*nu2 + 16.*(7. + 5.*f)*nu2*nu2*nu2*nu2))*help2 + (64.*(7. + 5.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 - 3.*(-10.*(4. + 3.*f) + (107. + 81.*f)*nu2 - 4.*(23. + 17.*f)*nu2*nu2 + 4.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-3.*(7. + 5.*f) + (466. + 342.*f)*nu2 - 8.*(81. + 59.*f)*nu2*nu2 + 32.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(69. + 51.*f - 4.*(179. + 135.*f)*nu2 + 4.*(289. + 215.*f)*nu2*nu2 - 8.*(81. + 59.*f)*nu2*nu2*nu2 + 16.*(7. + 5.*f)*nu2*nu2*nu2*nu2) - 1.*nu1*(321. + 243.*f - 16.*(116. + 89.*f)*nu2 + 16.*(179. + 135.*f)*nu2*nu2 - 8.*(233. + 171.*f)*nu2*nu2*nu2 + 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2))*help3))/(5.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*help1 + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*help2 + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*help3)); //GC!

                    //log1=clog((f*f*((-30.*(4. + 3.*f) + (145. + 99.*f)*nu2 + 64.*(7. + 5.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + (948. + 716.*f)*nu2*nu2 - 4.*(473. + 347.*f)*nu2*nu2*nu2 + 64.*(18. + 13.*f)*nu2*nu2*nu2*nu2 - 32.*(7. + 5.*f)*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(3.*(7. + 5.*f) + (438. + 322.*f)*nu2 - 16.*(51. + 37.*f)*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(-69. - 51.*f - 4.*(135. + 103.*f)*nu2 + 4.*(413. + 307.*f)*nu2*nu2 - 8.*(169. + 123.*f)*nu2*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2*nu2) + nu1*(321. + 243.*f + 8.*(99. + 79.*f)*nu2 - 8.*(605. + 459.*f)*nu2*nu2 + 8.*(751. + 557.*f)*nu2*nu2*nu2 - 256.*(11. + 8.*f)*nu2*nu2*nu2*nu2 + 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2*nu2))*help1 + (32.*(7. + 5.*f)*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(18. + 13.*f - 4.*(11. + 8.*f)*nu2 + 3.*(7. + 5.*f)*nu2*nu2) + 3.*(-10.*(4. + 3.*f) + (107. + 81.*f)*nu2 - 4.*(23. + 17.*f)*nu2*nu2 + 4.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-473. - 347.*f + 2.*(751. + 557.*f)*nu2 - 8.*(169. + 123.*f)*nu2*nu2 + 48.*(7. + 5.*f)*nu2*nu2*nu2) + nu1*(145. + 99.*f + 8.*(99. + 79.*f)*nu2 - 16.*(135. + 103.*f)*nu2*nu2 + 8.*(219. + 161.*f)*nu2*nu2*nu2 - 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(237. + 179.*f - 2.*(605. + 459.*f)*nu2 + 4.*(413. + 307.*f)*nu2*nu2 - 16.*(51. + 37.*f)*nu2*nu2*nu2 + 16.*(7. + 5.*f)*nu2*nu2*nu2*nu2))*help2 + (64.*(7. + 5.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 - 3.*(-10.*(4. + 3.*f) + (107. + 81.*f)*nu2 - 4.*(23. + 17.*f)*nu2*nu2 + 4.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-3.*(7. + 5.*f) + (466. + 342.*f)*nu2 - 8.*(81. + 59.*f)*nu2*nu2 + 32.*(7. + 5.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(69. + 51.*f - 4.*(179. + 135.*f)*nu2 + 4.*(289. + 215.*f)*nu2*nu2 - 8.*(81. + 59.*f)*nu2*nu2*nu2 + 16.*(7. + 5.*f)*nu2*nu2*nu2*nu2) - 1.*nu1*(321. + 243.*f - 16.*(116. + 89.*f)*nu2 + 16.*(179. + 135.*f)*nu2*nu2 - 8.*(233. + 171.*f)*nu2*nu2*nu2 + 64.*(7. + 5.*f)*nu2*nu2*nu2*nu2))*help3)));

                    //log2=clog((5.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*help1 + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*help2 + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*help3)));

                    //log21=clog((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*help1);

                    //log22=clog((224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*help2);

                    //log23=clog((60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*help3);

                    
                    pnlpt->M12_oneline_0_vv_complex[count] = f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2[count]) + f * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f3[count]);

                    //GC -> 2 of 14...

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_oneline_0_vv_complex_ortho[count] = f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2_ortho[count]) + f * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f3_ortho[count]);

                    //GC: ORTHOGONAL -- finish

                    //if (index_l < 10 && index_i < 10) { //GC -> FORGOT TO DECOMMENT THIS!

                    //printf("%e",pnlpt->M12_oneline_0_vv_complex[count]); //GC!!!
                    //printf("%e",log21); //GC!!!
                    //printf("\n");
                    //printf("%e",log22); //GC!!!
                    //printf("\n");
                    //printf("%e",log23); //GC!!!
                    //printf("\n");
                    //printf("%e",log21+log22+log23); //GC!!!
                    //printf("%e",log21+log22+log23-500.-500.*_Complex_I); //GC!!!
                    //printf("\n");
                    //printf("[][][][][][][][][][][][][][][][][][][][]\n");
                    //printf("%e",cexp(log21+log22+log23-500.-500.*_Complex_I)); //GC!!!
                    //printf("%e",cexp(300.+200.*_Complex_I)); //GC -> this works... P**** T****!!!
                    //printf("%e",cexp(1400.+1200.*_Complex_I)); //GC -> this DOES NOT work...
                    //printf("\n");
                    //printf("[*][*][*][*][*][*][*][*][*][*][*][*][*][*]\n");

                    //GC -> ANY OPERATION I CAN TRY WOULD LOSE ACCURACY, MOST LIKELY... GO BACK TO PYTHON...

                    //} //GC!
                        
                    }

                    count++;
                }
            }

            double complex *f13_0_vv;
            class_alloc(f13_0_vv, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_0_vv;
            class_alloc(f22_0_vv, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_0_vv;
            class_alloc(f12_0_vv, Nmax * sizeof(complex double), pnlpt->error_message); //GC!

            /*
        
        //GC: ORTHOGONAL -- start

        double complex *f12_0_vv_ortho;
        class_alloc(f12_0_vv_ortho,Nmax*sizeof(complex double),pnlpt->error_message); //GC!

        //GC: ORTHOGONAL -- finish
         
        //GC: NO, WE DO NOT NEED THIS, as discussed above...
         
        */

            index_i = 0;
            for (index_i = 0; index_i < Nmax + 1; index_i++)
            {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_0_vv_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * (3. * (f * f) * (7. * (-5. + 3. * nu1) + 6. * f * (-7. + 5. * nu1))) / 3920.;
            }

            index_j = 0;
            count = 0;
            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                }
                f13_0_vv[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_0_vv_oneline_complex, &inc);
                P13UV_0_vv[index_j] = -1. * Pbin[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav * (f * f * (441. + 566. * f + 175. * f * f) / 1225.);
                P13_0_vv[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_0_vv[index_j] * Pbin[index_j]) + P13UV_0_vv[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
            }

            index_j = 0;
            count = 0;

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                    
                    if (SWITCH_index == 1) {

                    x_transfer[count] = cmsym_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]); //GC!
                    //printf("%.16f",x_transfer[count]); //GC!
                        
                    }
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_0_vv_complex, x, &inc, &beta, y, &inc);
                
                if (SWITCH_index == 1) {

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_0_vv_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC -> again, y is not used??? No, see above. I could get away by overwriting stuff but whatever. It is better to avoid problems...
                //zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_0_vv_complex, x_transfer, &inc, &beta, y_transfer, &inc);
                //GC -> ok, the issue is the matrix...
                //printf("y_transfer is "); //GC!
                //printf("***********************************************************\n");
                //printf("%.16f",y_transfer); //GC!
                //printf("y is "); //GC!
                //printf("%.16f",y); //GC!
                //for (count=0; count < Nmax+1; count++){
                //    printf("%.16f",y_transfer[count]);
                //}
                //printf("***********************************************************\n");
                    
                }
                
                f22_0_vv[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                
                if (SWITCH_index == 1) {

                f12_0_vv[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC -> here I will have y_transfer...
                    
                }
                
                P22_0_vv[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_0_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                
                if (SWITCH_index == 1) {
                
                P12_0_vv[index_j] = Tbin[index_j] *                                                                               //[FACTOR] *
                                    creal(cpow(kdisc[index_j], 3.) * f12_0_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC! RECHECK THE CALCULATION BY HAND THAT ALSO IN REDSHIFT SPACE I HAVE THIS Tbin OVERALL!!!
                }
                
                //GC: ORTHOGONAL -- start

                if (SWITCH_index == 1) { //GC -> I put a switch here... This could have been written more compactly...
                
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_0_vv_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);
                f12_0_vv[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);
                P12_0_vv_ortho[index_j] = Tbin[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_0_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                    
                }

                //GC: do I print here? No, see below...

                //GC: ORTHOGONAL -- finish

                P1loop_0_vv[index_j] = Pbin[index_j] * (f * f / 5.) * 0. + (P13_0_vv[index_j] + P22_0_vv[index_j]); //GC -> notice the overall *0. that sets to zero the tree-level part, which is instead computed below...
                P_CTR_0[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j];
                Ptree_0_vv[index_j] = Pbin[index_j] * (f * f / 5.);
                //    printf("%le %le\n",kdisc[j],P22[j]);

                //printf("%.16f %.16f %.16f",kdisc[index_j],Tbin[index_j],P12_0_vv[index_j]); //GC!!!
                //printf("%.16f",f12_0_vv[index_j]); //GC!
                //printf("%.16f",f22_0_vv[index_j]); //GC!
                //printf("%.20f %.20f %.20f",kdisc[index_j],Tbin[index_j],P12_0_vv[index_j]);
                //printf("\n"); //GC!

                //printf("%.20f %.20f %.20f",kdisc[index_j],Tbin[index_j],P12_0_vv_ortho[index_j]);
                //printf("\n"); //GC!
            }

            // Computing P_{vd} contribution

            count = 0;
            index_l = 0;
            index_i = 0;
            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu12 = -0.5 * etam[index_i] - 0.5 * etam[index_l];
                    nu1 = -0.5 * etam[index_i];
                    nu2 = -0.5 * etam[index_l];
                    pnlpt->M22_oneline_0_vd_complex[count] = pnlpt->M22_oneline_complex[count] * 196. / ((nu1 * nu2 * (98. * nu12 * nu12 - 14. * nu12 + 36.) - 91. * nu12 * nu12 + 3. * nu12 + 58.)) * (f * (21. * f * f * (6. + 3. * nu2 - 10. * nu2 * nu2 + 2. * nu2 * nu2 * nu2 + 2. * nu1 * nu1 * nu1 * (1. + 5. * nu2) + 2 * nu1 * nu1 * (-5. - 2. * nu2 + 10. * nu2 * nu2) + nu1 * (3. - 24. * nu2 - 4. * nu2 * nu2 + 10. * nu2 * nu2 * nu2)) + 14. * f * (18. + 11. * nu2 + 42. * nu1 * nu1 * nu1 * nu2 - 31. * nu2 * nu2 + nu1 * nu1 * (-31. - 22. * nu2 + 84. * nu2 * nu2) + nu1 * (11. - 74. * nu2 - 22. * nu2 * nu2 + 42. * nu2 * nu2 * nu2)) + 5. * (46. + 13. * nu2 + 98. * nu1 * nu1 * nu1 * nu2 - 63. * nu2 * nu2 + 7. * nu1 * nu1 * (-9. - 10. * nu2 + 28. * nu2 * nu2) + nu1 * (13. - 138. * nu2 - 70. * nu2 * nu2 + 98. * nu2 * nu2 * nu2)))) / 1470.;

                    
                    if (SWITCH_index == 1) {

                    nu1 = -0.5 * etam_transfer[index_i];
                    nu2 = -0.5 * etam_transfer[index_l]; //GC!

                    //pnlpt->M12_oneline_0_vd_complex[count] = pnlpt->M12_oneline_complex[count]*(2.*f*((-30.*(15. + 7.*f) + (885. + 595.*f)*nu2 + 448.*(5. + 3.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 4.*(845. + 371.*f)*nu2*nu2 - 4.*(2045. + 1099.*f)*nu2*nu2*nu2 + 64.*(85. + 49.*f)*nu2*nu2*nu2*nu2 - 224.*(5. + 3.*f)*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21.*(5. + 3.*f) + 2.*(935. + 497.*f)*nu2 - 16.*(235. + 133.*f)*nu2*nu2 + 336.*(5. + 3.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(-3.*(95. + 49.*f) - 20.*(89. + 35.*f)*nu2 + 28.*(235. + 117.*f)*nu2*nu2 - 8.*(765. + 427.*f)*nu2*nu2*nu2 + 336.*(5. + 3.*f)*nu2*nu2*nu2*nu2) + nu1*(1125. + 483.*f + 8.*(205. + 7.*f)*nu2 - 40.*(417. + 175.*f)*nu2*nu2 + 8.*(3035. + 1533.*f)*nu2*nu2*nu2 - 512.*(25. + 14.*f)*nu2*nu2*nu2*nu2 + 448.*(5. + 3.*f)*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*(5. + 3.*f)*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(85. + 49.*f - 8.*(25. + 14.*f)*nu2 + 21.*(5. + 3.*f)*nu2*nu2) + 3.*(-10.*(15. + 7.*f) + (375. + 161.*f)*nu2 - 4.*(95. + 49.*f)*nu2*nu2 + 28.*(5. + 3.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-2045. - 1099.*f + (6070. + 3066.*f)*nu2 - 8.*(765. + 427.*f)*nu2*nu2 + 336.*(5. + 3.*f)*nu2*nu2*nu2) + nu1*(885. + 595.*f + 8.*(205. + 7.*f)*nu2 - 80.*(89. + 35.*f)*nu2*nu2 + 8.*(935. + 497.*f)*nu2*nu2*nu2 - 448.*(5. + 3.*f)*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(845. + 371.*f - 10.*(417. + 175.*f)*nu2 + 28.*(235. + 117.*f)*nu2*nu2 - 16.*(235. + 133.*f)*nu2*nu2*nu2 + 112.*(5. + 3.*f)*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (448.*(5. + 3.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 - 3.*(-10.*(15. + 7.*f) + (375. + 161.*f)*nu2 - 4.*(95. + 49.*f)*nu2*nu2 + 28.*(5. + 3.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-21.*(5. + 3.*f) + 2.*(1005. + 539.*f)*nu2 - 8.*(365. + 203.*f)*nu2*nu2 + 224.*(5. + 3.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(3.*(95. + 49.*f) - 4.*(645. + 287.*f)*nu2 + 4.*(1145. + 567.*f)*nu2*nu2 - 8.*(365. + 203.*f)*nu2*nu2*nu2 + 112.*(5. + 3.*f)*nu2*nu2*nu2*nu2) - 1.*nu1*(1125. + 483.*f - 16.*(365. + 133.*f)*nu2 + 16.*(645. + 287.*f)*nu2*nu2 - 8.*(1005. + 539.*f)*nu2*nu2*nu2 + 448.*(5. + 3.*f)*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))/(15.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)));

                    
                    pnlpt->M12_oneline_0_vd_complex[count] = f * (pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1[count]) + f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd0_f2[count]);

                    //GC -> 4 of 14...

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_oneline_0_vd_complex_ortho[count] = f * (pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1_ortho[count]) + f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd0_f2_ortho[count]);

                    //GC: ORTHOGONAL -- finish
                    }
                    //GC!

                    count++;
                }
            }

            index_i = 0;

            for (index_i = 0; index_i < Nmax + 1; index_i++)
            {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_0_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * (f * (-35. - 18. * f + 45. * nu1 + 54. * f * nu1)) / 840.;
            }

            double complex *f13_0_vd;
            class_alloc(f13_0_vd, Nmax * sizeof(complex double), pnlpt->error_message);

            double complex *f22_0_vd;
            class_alloc(f22_0_vd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_0_vd;
            class_alloc(f12_0_vd, Nmax * sizeof(complex double), pnlpt->error_message); //GC!

            index_j = 0;
            count = 0;
            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                }
                f13_0_vd[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_0_vd_oneline_complex, &inc);
                P13UV_0_vd[index_j] = -1. * Pbin[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav * (2. * f * (625. + 558. * f + 315. * f * f) / 1575.);
                P13_0_vd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_0_vd[index_j] * Pbin[index_j]) + P13UV_0_vd[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                    
                    if (SWITCH_index == 1) {

                    x_transfer[count] = cmsym_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]); //GC!
                        
                    }
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_0_vd_complex, x, &inc, &beta, y, &inc);
                
                if (SWITCH_index == 1) {

                
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_0_vd_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC!
                    
                }
                
                f22_0_vd[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                
                if (SWITCH_index == 1) {

                f12_0_vd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                    
                }
                
                P22_0_vd[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_0_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                
                if (SWITCH_index == 1) {

                P12_0_vd[index_j] = Tbin[index_j] *                                                                               //[FACTOR] *
                                    creal(cpow(kdisc[index_j], 3.) * f12_0_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                    
                P1loop_0_vd[index_j] = Pbin[index_j] * (f * 2. / 3.) * 0. + (P13_0_vd[index_j] + P22_0_vd[index_j]);
                Ptree_0_vd[index_j] = Pbin[index_j] * (f * 2. / 3.);

                //printf("%.20f %.20f %.20f",kdisc[index_j],Tbin[index_j],P12_0_vd[index_j]);
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- start

                if (SWITCH_index == 1) {
                
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_0_vd_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);                        //GC!
                f12_0_vd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                                 //GC!
                P12_0_vd_ortho[index_j] = Tbin[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_0_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!

                }
                    
                //GC: ORTHOGONAL -- finish

                //printf("%.20f %.20f %.20f",kdisc[index_j],Tbin[index_j],P12_0_vd_ortho[index_j]);
                //printf("\n"); //GC!
            }

            // Computing P_{dd} contribution

            count = 0;
            index_l = 0;
            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu12 = -0.5 * etam[index_i] - 0.5 * etam[index_l];
                    nu1 = -0.5 * etam[index_i];
                    nu2 = -0.5 * etam[index_l];
                    pnlpt->M22_oneline_0_dd_complex[count] = pnlpt->M22_oneline_complex[count] * 196. / ((nu1 * nu2 * (98. * nu12 * nu12 - 14. * nu12 + 36.) - 91. * nu12 * nu12 + 3. * nu12 + 58.)) * (98. * f * f * (4. - 2. * nu2 - 5. * nu2 * nu2 + nu2 * nu2 * nu2 + nu1 * nu1 * nu1 * (1. + 3. * nu2) + nu1 * nu1 * (-5. + 2. * nu2 + 6. * nu2 * nu2) + nu1 * (-2. - 4. * nu2 + 2. * nu2 * nu2 + 3. * nu2 * nu2 * nu2)) + 70. * f * (10. - nu2 + 14. * nu1 * nu1 * nu1 * nu2 - 17. * nu2 * nu2 + nu1 * nu1 * (-17. + 6. * nu2 + 28. * nu2 * nu2) + nu1 * (-1. - 22. * nu2 + 6. * nu2 * nu2 + 14. * nu2 * nu2 * nu2)) + 15. * (58. + 3. * nu2 + 98. * nu1 * nu1 * nu1 * nu2 - 91. * nu2 * nu2 + 7. * nu1 * nu1 * (-13. - 2. * nu2 + 28. * nu2 * nu2) + nu1 * (3. - 146. * nu2 - 14. * nu2 * nu2 + 98. * nu2 * nu2 * nu2))) / 2940.;

                    if (SWITCH_index == 1) {

                    nu1 = -0.5 * etam_transfer[index_i];
                    nu2 = -0.5 * etam_transfer[index_l]; //GC!

                    //pnlpt->M12_oneline_0_dd_complex[count] = pnlpt->M12_oneline_complex[count]*((-180. + 3.*(209. + 91.*f)*nu2 + 448.*(3. + f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 - 4.*(-303. + 35.*f)*nu2*nu2 - 4.*(1035. + 217.*f)*nu2*nu2*nu2 + 128.*(24. + 7.*f)*nu2*nu2*nu2*nu2 - 224.*(3. + f)*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21.*(3. + f) + 2.*(465. + 91.*f)*nu2 - 16.*(129. + 35.*f)*nu2*nu2 + 336.*(3. + f)*nu2*nu2*nu2) + 4.*nu1*nu1*(-3.*(45. + 7.*f) + 4.*(-129. + 49.*f)*nu2 + 28.*(105. + 11.*f)*nu2*nu2 - 24.*(137. + 35.*f)*nu2*nu2*nu2 + 336.*(3. + f)*nu2*nu2*nu2*nu2) + nu1*(387. - 63.*f - 8.*(51. + 133.*f)*nu2 + 24.*(-229. + 49.*f)*nu2*nu2 + 8.*(1389. + 175.*f)*nu2*nu2*nu2 - 256.*(27. + 7.*f)*nu2*nu2*nu2*nu2 + 448.*(3. + f)*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*(3. + f)*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(48. + 14.*f - 4.*(27. + 7.*f)*nu2 + 21.*(3. + f)*nu2*nu2) + 3.*(-60. - 3.*(-43. + 7.*f)*nu2 - 4.*(45. + 7.*f)*nu2*nu2 + 28.*(3. + f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-1035. - 217.*f + (2778. + 350.*f)*nu2 - 24.*(137. + 35.*f)*nu2*nu2 + 336.*(3. + f)*nu2*nu2*nu2) + nu1*(627. + 273.*f - 8.*(51. + 133.*f)*nu2 + 16.*(-129. + 49.*f)*nu2*nu2 + 8.*(465. + 91.*f)*nu2*nu2*nu2 - 448.*(3. + f)*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(303. - 35.*f + 6.*(-229. + 49.*f)*nu2 + 28.*(105. + 11.*f)*nu2*nu2 - 16.*(129. + 35.*f)*nu2*nu2*nu2 + 112.*(3. + f)*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (180. + 9.*(-43. + 7.*f)*nu2 + 448.*(3. + f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 12.*(45. + 7.*f)*nu2*nu2 - 84.*(3. + f)*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21.*(3. + f) + 6.*(169. + 35.*f)*nu2 - 8.*(195. + 49.*f)*nu2*nu2 + 224.*(3. + f)*nu2*nu2*nu2) + nu1*(9.*(-43. + 7.*f) - 32.*(-45. + 28.*f)*nu2 + 48.*(-79. + 7.*f)*nu2*nu2 + 24.*(169. + 35.*f)*nu2*nu2*nu2 - 448.*(3. + f)*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(3.*(45. + 7.*f) + 12.*(-79. + 7.*f)*nu2 + 4.*(507. + 49.*f)*nu2*nu2 - 8.*(195. + 49.*f)*nu2*nu2*nu2 + 112.*(3. + f)*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI))/(3.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI))); //GC!

                    
                    pnlpt->M12_oneline_0_dd_complex[count] = (pnlpt->M12_oneline_complex_matter_multipoles_dd0_f0[count]) + f * (pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1[count]);

                    //GC -> 6 of 14...

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_oneline_0_dd_complex_ortho[count] = (pnlpt->M12_oneline_complex_matter_multipoles_dd0_f0_ortho[count]) + f * (pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1_ortho[count]);

                    //GC: ORTHOGONAL -- finish
                        
                    }

                    count++;
                }
            }
            index_i = 0;
            for (index_i = 0; index_i < Nmax + 1; index_i++)
            {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_0_dd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] / (1. + 9. * nu1) * (1. + 9. * nu1 + 6. * f * (1. + nu1));
            }

            double complex *f13_0_dd;
            class_alloc(f13_0_dd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_0_dd;
            class_alloc(f22_0_dd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_0_dd;
            class_alloc(f12_0_dd, Nmax * sizeof(complex double), pnlpt->error_message); //GC!
            index_j = 0;
            count = 0;
            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                }
                f13_0_dd[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_0_dd_oneline_complex, &inc);
                P13UV_0_dd[index_j] = -1. * Pbin[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav * ((61. - 2. * f + 35. * f * f) / 105.);
                P13_0_dd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_0_dd[index_j] * Pbin[index_j]) + P13UV_0_dd[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                    
                    if (SWITCH_index == 1) {
                    x_transfer[count] = cmsym_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]); //GC!
                    }
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_0_dd_complex, x, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_0_dd_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_0_dd[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_0_dd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_0_dd[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_0_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_0_dd[index_j] = Tbin[index_j] *                                                                               //[FACTOR] *
                                    creal(cpow(kdisc[index_j], 3.) * f12_0_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                P1loop_0_dd[index_j] = Pbin[index_j] * 0. + (P13_0_dd[index_j] + P22_0_dd[index_j]);
                Ptree_0_dd[index_j] = Pbin[index_j];

                //printf("%.24f %.24f %.24f",kdisc[index_j],Tbin[index_j],P12_0_dd[index_j]);
                //printf("\n"); //GC!
                if (SWITCH_index == 1) {
                //GC: ORTHOGONAL -- start

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_0_dd_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);                        //GC!
                f12_0_dd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                                 //GC!
                P12_0_dd_ortho[index_j] = Tbin[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_0_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!

                //GC: ORTHOGONAL -- finish
                }
                //        printf("%.20f %.20f %.20f",kdisc[index_j],Tbin[index_j],P12_0_dd_ortho[index_j]);
                //        printf("\n"); //GC!
            }

            // Computing P_{vv} contribution - Quadrupole

            count = 0;
            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    //  nu12 =-0.5*etam[index_i]-0.5*etam[index_l];
                    nu1 = -0.5 * etam[index_i];
                    nu2 = -0.5 * etam[index_l];

                    nu12 = nu1 + nu2;

                    pnlpt->M22_oneline_2_vv_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * (f * f * (396. * (50. - 9. * nu2 + 98. * nu1 * nu1 * nu1 * nu2 - 35. * nu2 * nu2 + 7. * nu1 * nu1 * (-5. - 18. * nu2 + 28. * nu2 * nu2) + nu1 * (-9. - 66. * nu2 - 126. * nu2 * nu2 + 98. * nu2 * nu2 * nu2)) + 231. * f * (142. - 21. * nu2 + 280. * nu1 * nu1 * nu1 * nu2 - 106. * nu2 * nu2 + 2. * nu1 * nu1 * (-53. - 174. * nu2 + 280. * nu2 * nu2) + nu1 * (-21. - 204. * nu2 - 348. * nu2 * nu2 + 280. * nu2 * nu2 * nu2)) + 49. * f * f * (336. - 62. * nu2 - 255. * nu2 * nu2 + 50. * nu2 * nu2 * nu2 + 10. * nu1 * nu1 * nu1 * (5. + 56. * nu2) + 5. * nu1 * nu1 * (-51. - 142. * nu2 + 224. * nu2 * nu2) + nu1 * (-62. - 486. * nu2 - 710. * nu2 * nu2 + 560. * nu2 * nu2 * nu2)))) / 135828.;
                
                    if (SWITCH_index == 1) {
                    nu1 = -0.5 * etam_transfer[index_i];
                    nu2 = -0.5 * etam_transfer[index_l]; //GC!

                    //pnlpt->M12_oneline_2_vv_complex[count] = (pnlpt->M12_oneline_complex[count])*(f*f*((-15.*(96. + 77.*f) + 2.*(870. + 749.*f)*nu2 + 896.*(6. + 5.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 144.*(79. + 63.*f)*nu2*nu2 - 8.*(2838. + 2317.*f)*nu2*nu2*nu2 + 96.*(144. + 119.*f)*nu2*nu2*nu2*nu2 - 448.*(6. + 5.*f)*nu2*nu2*nu2*nu2*nu2 + 24.*nu1*nu1*nu1*(7.*(6. + 5.*f) + (876. + 714.*f)*nu2 - 96.*(17. + 14.*f)*nu2*nu2 + 112.*(6. + 5.*f)*nu2*nu2*nu2) + 16.*nu1*nu1*(-3.*(69. + 56.*f) - 3.*(540. + 427.*f)*nu2 + 28.*(177. + 143.*f)*nu2*nu2 - 4.*(1014. + 833.*f)*nu2*nu2*nu2 + 168.*(6. + 5.*f)*nu2*nu2*nu2*nu2) + nu1*(6.*(642. + 511.*f) + 24.*(396. + 301.*f)*nu2 - 32.*(1815. + 1442.*f)*nu2*nu2 + 16.*(4506. + 3647.*f)*nu2*nu2*nu2 - 128.*(264. + 217.*f)*nu2*nu2*nu2*nu2 + 896.*(6. + 5.*f)*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (448.*(6. + 5.*f)*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 32.*nu1*nu1*nu1*nu1*(432. + 357.*f - 4.*(264. + 217.*f)*nu2 + 84.*(6. + 5.*f)*nu2*nu2) + 3.*(-5.*(96. + 77.*f) + 2.*(642. + 511.*f)*nu2 - 16.*(69. + 56.*f)*nu2*nu2 + 56.*(6. + 5.*f)*nu2*nu2*nu2) + 8.*nu1*nu1*nu1*(-2838. - 2317.*f + (9012. + 7294.*f)*nu2 - 8.*(1014. + 833.*f)*nu2*nu2 + 336.*(6. + 5.*f)*nu2*nu2*nu2) + 2.*nu1*(870. + 749.*f + 12.*(396. + 301.*f)*nu2 - 24.*(540. + 427.*f)*nu2*nu2 + 72.*(146. + 119.*f)*nu2*nu2*nu2 - 448.*(6. + 5.*f)*nu2*nu2*nu2*nu2) + 16.*nu1*nu1*(9.*(79. + 63.*f) - 2.*(1815. + 1442.*f)*nu2 + 28.*(177. + 143.*f)*nu2*nu2 - 144.*(17. + 14.*f)*nu2*nu2*nu2 + 56.*(6. + 5.*f)*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (896.*(6. + 5.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 - 3.*(-5.*(96. + 77.*f) + 2.*(642. + 511.*f)*nu2 - 16.*(69. + 56.*f)*nu2*nu2 + 56.*(6. + 5.*f)*nu2*nu2*nu2) + 8.*nu1*nu1*nu1*(-21.*(6. + 5.*f) + (2796. + 2282.*f)*nu2 - 24.*(162. + 133.*f)*nu2*nu2 + 224.*(6. + 5.*f)*nu2*nu2*nu2) + 16.*nu1*nu1*(3.*(69. + 56.*f) - 1.*(2148. + 1715.*f)*nu2 + 4.*(867. + 700.*f)*nu2*nu2 - 12.*(162. + 133.*f)*nu2*nu2*nu2 + 56.*(6. + 5.*f)*nu2*nu2*nu2*nu2) - 2.*nu1*(3.*(642. + 511.*f) - 4.*(2784. + 2191.*f)*nu2 + 8.*(2148. + 1715.*f)*nu2*nu2 - 8.*(1398. + 1141.*f)*nu2*nu2*nu2 + 448.*(6. + 5.*f)*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))/(21.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)));

                    pnlpt->M12_oneline_2_vv_complex[count] = 20. / 7. * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2[count]) + f * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv2_f3[count]);

                    //GC -> 8 of 14...

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_oneline_2_vv_complex_ortho[count] = 20. / 7. * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2_ortho[count]) + f * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv2_f3_ortho[count]);
                        
                        //printf("%d\n",SWITCH_index); //GC - SWITCH!


                    //GC: ORTHOGONAL -- finish
                    }
                    //GC!

                    count++;
                }
            }

            double complex *f13_2_vv;
            class_alloc(f13_2_vv, Nmax * sizeof(complex double), pnlpt->error_message);
            //    double *P_CTR0;
            //    class_alloc(P_CTR0,Nmax*sizeof(double),pnlpt->error_message);

            double complex *f22_2_vv;
            class_alloc(f22_2_vv, Nmax * sizeof(complex double), pnlpt->error_message);

            double complex *f12_2_vv;
            class_alloc(f12_2_vv, Nmax * sizeof(complex double), pnlpt->error_message);
            //GC!

            for (index_i = 0; index_i < Nmax + 1; index_i++)
            {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_2_vv_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * (3. * f * f * (-5. + 3. * nu1 + f * (-6. + 5. * nu1))) / 196.;
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                }
                f13_2_vv[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_2_vv_oneline_complex, &inc);
                P13UV_2_vv[index_j] = -1. * Pbin[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav * (2. * f * f * (54. + 74. * f + 25. * f * f) / 105.);
                P13_2_vv[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_2_vv[index_j] * Pbin[index_j]) + P13UV_2_vv[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                    if (SWITCH_index == 1) {
                        
                    x_transfer[count] = cmsym_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]); //GC!
                    }
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_2_vv_complex, x, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_2_vv_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_2_vv[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_2_vv[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_2_vv[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_2_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_2_vv[index_j] = Tbin[index_j] *                                                                               //[FACTOR] *
                                    creal(cpow(kdisc[index_j], 3.) * f12_2_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                P1loop_2_vv[index_j] = Pbin[index_j] * (f * f * 4. / 7.) * 0. + (P13_2_vv[index_j] + P22_2_vv[index_j]);
                P_CTR_2[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j] * f * 2. / 3.;
                Ptree_2_vv[index_j] = Pbin[index_j] * (f * f * 4. / 7.);
                //    printf("%le %le\n",kdisc[j],P22[j]);

                //printf("%.24f %.24f %.24f",kdisc[index_j],Tbin[index_j],P12_2_vv[index_j]);
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- start

                if (SWITCH_index == 1) {

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_2_vv_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);                        //GC!
                f12_2_vv[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                                 //GC!
                P12_2_vv_ortho[index_j] = Tbin[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_2_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!

                }
                    
                //GC: ORTHOGONAL -- finish

                //        printf("%.20f %.20f %.20f",kdisc[index_j],Tbin[index_j],P12_2_vv_ortho[index_j]);
                //        printf("\n"); //GC!
            }

            // Computing P_{vd} contribution - Quadrupole

            count = 0;
            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    //  nu12 =-0.5*etam[index_i]-0.5*etam[index_l];
                    nu1 = -0.5 * etam[index_i];
                    nu2 = -0.5 * etam[index_l];

                    nu12 = nu1 + nu2;

                    pnlpt->M22_oneline_2_vd_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * (f * (7. * f * f * (22. + 11. * nu2 - 40. * nu2 * nu2 + 4. * nu2 * nu2 * nu2 + nu1 * nu1 * nu1 * (4. + 40. * nu2) + 8. * nu1 * nu1 * (-5. - nu2 + 10. * nu2 * nu2) + nu1 * (11. - 88. * nu2 - 8. * nu2 * nu2 + 40. * nu2 * nu2 * nu2)) + 4. * (46. + 13. * nu2 + 98. * nu1 * nu1 * nu1 * nu2 - 63. * nu2 * nu2 + 7. * nu1 * nu1 * (-9. - 10. * nu2 + 28. * nu2 * nu2) + nu1 * (13. - 138. * nu2 - 70. * nu2 * nu2 + 98. * nu2 * nu2 * nu2)) + f * (306. + 161. * nu2 + 672. * nu1 * nu1 * nu1 * nu2 - 538. * nu2 * nu2 + 2. * nu1 * nu1 * (-269. - 134. * nu2 + 672. * nu2 * nu2) + nu1 * (161. - 1196. * nu2 - 268. * nu2 * nu2 + 672. * nu2 * nu2 * nu2)))) / 588.;

                    if (SWITCH_index == 1) {

                    nu1 = -0.5 * etam_transfer[index_i];
                    nu2 = -0.5 * etam_transfer[index_l]; //GC!

                    //pnlpt->M12_oneline_2_vd_complex[count] = (pnlpt->M12_oneline_complex[count])*(f*((-15.*(24. + 13.*f) + (708. + 728.*f)*nu2 + 256.*(7. + 6.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 8.*(338. + 161.*f)*nu2*nu2 - 16.*(409. + 290.*f)*nu2*nu2*nu2 + 32.*(136. + 109.*f)*nu2*nu2*nu2*nu2 - 128.*(7. + 6.*f)*nu2*nu2*nu2*nu2*nu2 + 16.*nu1*nu1*nu1*(3.*(7. + 6.*f) + (374. + 260.*f)*nu2 - 8.*(94. + 73.*f)*nu2*nu2 + 48.*(7. + 6.*f)*nu2*nu2*nu2) + 8.*nu1*nu1*(-3.*(38. + 25.*f) - 2.*(356. + 131.*f)*nu2 + 4.*(658. + 405.*f)*nu2*nu2 - 16.*(153. + 116.*f)*nu2*nu2*nu2 + 96.*(7. + 6.*f)*nu2*nu2*nu2*nu2) + 4.*nu1*(3.*(75. + 34.*f) + (328. - 158.*f)*nu2 - 4.*(834. + 359.*f)*nu2*nu2 + 8.*(607. + 384.*f)*nu2*nu2*nu2 - 32.*(80. + 61.*f)*nu2*nu2*nu2*nu2 + 64.*(7. + 6.*f)*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (128.*(7. + 6.*f)*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 32.*nu1*nu1*nu1*nu1*(136. + 109.*f - 4.*(80. + 61.*f)*nu2 + 24.*(7. + 6.*f)*nu2*nu2) + 3.*(-5.*(24. + 13.*f) + 4.*(75. + 34.*f)*nu2 - 8.*(38. + 25.*f)*nu2*nu2 + 16.*(7. + 6.*f)*nu2*nu2*nu2) + 16.*nu1*nu1*nu1*(-409. - 290.*f + 2.*(607. + 384.*f)*nu2 - 8.*(153. + 116.*f)*nu2*nu2 + 48.*(7. + 6.*f)*nu2*nu2*nu2) + 8.*nu1*nu1*(338. + 161.*f - 2.*(834. + 359.*f)*nu2 + 4.*(658. + 405.*f)*nu2*nu2 - 16.*(94. + 73.*f)*nu2*nu2*nu2 + 32.*(7. + 6.*f)*nu2*nu2*nu2*nu2) - 4.*nu1*(-177. - 182.*f + 2.*(-164. + 79.*f)*nu2 + 4.*(356. + 131.*f)*nu2*nu2 - 8.*(187. + 130.*f)*nu2*nu2*nu2 + 64.*(7. + 6.*f)*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (256.*(7. + 6.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 3.*(120. + 65.*f - 4.*(75. + 34.*f)*nu2 + 8.*(38. + 25.*f)*nu2*nu2 - 16.*(7. + 6.*f)*nu2*nu2*nu2) + 16.*nu1*nu1*nu1*(-3.*(7. + 6.*f) + (402. + 284.*f)*nu2 - 8.*(73. + 55.*f)*nu2*nu2 + 32.*(7. + 6.*f)*nu2*nu2*nu2) + 8.*nu1*nu1*(114. + 75.*f - 2.*(516. + 253.*f)*nu2 + 4.*(458. + 279.*f)*nu2*nu2 - 16.*(73. + 55.*f)*nu2*nu2*nu2 + 32.*(7. + 6.*f)*nu2*nu2*nu2*nu2) - 4.*nu1*(3.*(75. + 34.*f) - 2.*(584. + 175.*f)*nu2 + 4.*(516. + 253.*f)*nu2*nu2 - 8.*(201. + 142.*f)*nu2*nu2*nu2 + 64.*(7. + 6.*f)*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))/(3.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)));

                    pnlpt->M12_oneline_2_vd_complex[count] = 2. * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1[count]) + f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd2_f2[count]);

                    //GC -> 10 of 14...

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_oneline_2_vd_complex_ortho[count] = 2. * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1_ortho[count]) + f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd2_f2_ortho[count]);
                    }
                    //GC: ORTHOGONAL -- finish

                    //GC!

                    count++;
                }
            }

            double complex *f13_2_vd;
            class_alloc(f13_2_vd, Nmax * sizeof(complex double), pnlpt->error_message);
            //    double *P_CTR0;
            //    class_alloc(P_CTR0,Nmax*sizeof(double),pnlpt->error_message);

            double complex *f22_2_vd;
            class_alloc(f22_2_vd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_2_vd;
            class_alloc(f12_2_vd, Nmax * sizeof(complex double), pnlpt->error_message); //GC!

            for (index_i = 0; index_i < Nmax + 1; index_i++)
            {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_2_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * (f * (-49. - 9. * f + 63. * nu1 + 108. * f * nu1)) / 588.;
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                }
                f13_2_vd[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_2_vd_oneline_complex, &inc);
                P13UV_2_vd[index_j] = -1. * Pbin[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav * 4. * f * (175. + 180. * f + 126. * f * f) / 441.;
                P13_2_vd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_2_vd[index_j] * Pbin[index_j]) + P13UV_2_vd[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                    
                    if (SWITCH_index == 1) {

                    x_transfer[count] = cmsym_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]); //GC!
                        
                    }
                    
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_2_vd_complex, x, &inc, &beta, y, &inc);
                
                if (SWITCH_index == 1) {

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_2_vd_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC!
                    
                }
                f22_2_vd[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                
                if (SWITCH_index == 1) {

                f12_2_vd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                    
                }
                P22_2_vd[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_2_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                
                if (SWITCH_index == 1) {

                P12_2_vd[index_j] = Tbin[index_j] *                                                                               //[FACTOR] *
                                    creal(cpow(kdisc[index_j], 3.) * f12_2_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                P1loop_2_vd[index_j] = Pbin[index_j] * (f * 4. / 3.) * 0. + (P13_2_vd[index_j] + P22_2_vd[index_j]);
                Ptree_2_vd[index_j] = Pbin[index_j] * (f * 4. / 3.);
                //    P_CTR[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j];
                //printf("%le %le\n",kdisc[j],P22[j]);

                //printf("%.24f %.24f %.24f",kdisc[index_j],Tbin[index_j],P12_2_vd[index_j]);
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- start
                if (SWITCH_index == 1) {

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_2_vd_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);                        //GC!
                f12_2_vd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                                 //GC!
                P12_2_vd_ortho[index_j] = Tbin[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_2_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //GC: ORTHOGONAL -- finish

                //printf("%.20f %.20f %.20f",kdisc[index_j],Tbin[index_j],P12_2_vd_ortho[index_j]);
                //printf("\n"); //GC!
            }

            // Computing P_{dd} contribution - Quadrupole

            count = 0;
            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    //  nu12 =-0.5*etam[index_i]-0.5*etam[index_l];
                    nu1 = -0.5 * etam[index_i];
                    nu2 = -0.5 * etam[index_l];

                    nu12 = nu1 + nu2;

                    pnlpt->M22_oneline_2_dd_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * f * (4. * (10. - nu2 + 14. * nu1 * nu1 * nu1 * nu2 - 17. * nu2 * nu2 + nu1 * nu1 * (-17. + 6. * nu2 + 28. * nu2 * nu2) + nu1 * (-1. - 22. * nu2 + 6. * nu2 * nu2 + 14. * nu2 * nu2 * nu2)) + f * (26. - 13. * nu2 - 37. * nu2 * nu2 + 2. * nu2 * nu2 * nu2 + nu1 * nu1 * nu1 * (2. + 24. * nu2) + nu1 * nu1 * (-37. + 22. * nu2 + 48. * nu2 * nu2) + nu1 * (-13. - 26. * nu2 + 22. * nu2 * nu2 + 24. * nu2 * nu2 * nu2))) / 84.;

                    if (SWITCH_index == 1) {

                    nu1 = -0.5 * etam_transfer[index_i];
                    nu2 = -0.5 * etam_transfer[index_l]; //GC!

                    //pnlpt->M12_oneline_2_dd_complex[count] = (pnlpt->M12_oneline_complex[count])*(14.*f*((64.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 4.*nu1*nu1*nu1*(3. + 26.*nu2 - 80.*nu2*nu2 + 48.*nu2*nu2*nu2) + nu2*(39. - 20.*nu2 - 124.*nu2*nu2 + 128.*nu2*nu2*nu2 - 32.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(-3. + 28.*nu2 + 44.*nu2*nu2 - 120.*nu2*nu2*nu2 + 48.*nu2*nu2*nu2*nu2) + nu1*(-9. - 152.*nu2 + 168.*nu2*nu2 + 200.*nu2*nu2*nu2 - 256.*nu2*nu2*nu2*nu2 + 64.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (32.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(2. - 4.*nu2 + 3.*nu2*nu2) + 3.*nu2*(-3. - 4.*nu2 + 4.*nu2*nu2) + 4.*nu1*nu1*nu1*(-31. + 50.*nu2 - 120.*nu2*nu2 + 48.*nu2*nu2*nu2) + nu1*(39. - 152.*nu2 + 112.*nu2*nu2 + 104.*nu2*nu2*nu2 - 64.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(-5. + 42.*nu2 + 44.*nu2*nu2 - 80.*nu2*nu2*nu2 + 16.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (64.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 3.*nu2*(3. + 4.*nu2 - 4.*nu2*nu2) + 4.*nu1*nu1*nu1*(-3. + 30.*nu2 - 56.*nu2*nu2 + 32.*nu2*nu2*nu2) + nu1*(9. - 128.*nu2 + 48.*nu2*nu2 + 120.*nu2*nu2*nu2 - 64.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(3. + 12.*nu2 + 28.*nu2*nu2 - 56.*nu2*nu2*nu2 + 16.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))/(3.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)));

                    pnlpt->M12_oneline_2_dd_complex[count] = 2. * f * (pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1[count]);

                    //GC -> 11 of 14...

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_oneline_2_dd_complex_ortho[count] = 2. * f * (pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1_ortho[count]);

                    //GC: ORTHOGONAL -- finish

                    //GC!
                    }
                    count++;
                }
            }

            double complex *f13_2_dd;
            class_alloc(f13_2_dd, Nmax * sizeof(complex double), pnlpt->error_message);
            //    double *P_CTR0;
            //    class_alloc(P_CTR0,Nmax*sizeof(double),pnlpt->error_message);

            double complex *f22_2_dd;
            class_alloc(f22_2_dd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_2_dd;
            class_alloc(f12_2_dd, Nmax * sizeof(complex double), pnlpt->error_message); //GC!

            for (index_i = 0; index_i < Nmax + 1; index_i++)
            {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_2_dd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * (3. * f * (1. + nu1)) / 28.;
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                }
                f13_2_dd[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_2_dd_oneline_complex, &inc);
                P13UV_2_dd[index_j] = -1. * Pbin[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav * (2. * f * (35. * f - 2.) / 105.);
                P13_2_dd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_2_dd[index_j] * Pbin[index_j]) + P13UV_2_dd[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                    if (SWITCH_index == 1) {

                    x_transfer[count] = cmsym_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]); //GC!
                    }
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_2_dd_complex, x, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_2_dd_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_2_dd[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {

                f12_2_dd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_2_dd[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_2_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {

                P12_2_dd[index_j] = Tbin[index_j] *                                                                               //[FACTOR] *
                                    creal(cpow(kdisc[index_j], 3.) * f12_2_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                P1loop_2_dd[index_j] = P13_2_dd[index_j] + P22_2_dd[index_j];
                //    P_CTR[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j];
                //    printf("%le %le\n",kdisc[j],P22[j]);

                //printf("%.24f %.24f %.24f",kdisc[index_j],Tbin[index_j],P12_2_dd[index_j]);
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- start
                if (SWITCH_index == 1) {

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_2_dd_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);                        //GC!
                f12_2_dd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                                 //GC!
                P12_2_dd_ortho[index_j] = Tbin[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_2_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //GC: ORTHOGONAL -- finish

                //printf("%.20f %.20f %.20f",kdisc[index_j],Tbin[index_j],P12_2_dd_ortho[index_j]);
                //printf("\n"); //GC!
            }

            // Computing P_{vv} contribution - Hexadecapole

            //GC -> is it possible that the hexadecapole is large??? Here it is divided by one million, while I do NOT have any such division apparently...

            count = 0;
            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam[index_i];
                    nu2 = -0.5 * etam[index_l];
                    nu12 = nu1 + nu2;

                    pnlpt->M22_oneline_4_vv_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * (f * f * (1144. * (50. + 98. * nu1 * nu1 * nu1 * nu2 - nu2 * (9. + 35. * nu2) + 7. * nu1 * nu1 * (-5. + 2. * nu2 * (-9. + 14. * nu2)) + nu1 * (-9. + 2. * nu2 * (-33. + 7. * nu2 * (-9. + 7. * nu2)))) + 147. * f * f * (483. + 40. * nu1 * nu1 * nu1 * (-1. + 28. * nu2) - 2. * nu2 * (-57. + 10. * nu2 * (29. + 2. * nu2)) + 20. * nu1 * nu1 * (-29. + 2. * nu2 * (-25. + 56. * nu2)) + 2. * nu1 * (57. + 2. * nu2 * (-327. + 10. * nu2 * (-25. + 28. * nu2)))) + 728. * f * (206. + 420. * nu1 * nu1 * nu1 * nu2 + (7. - 208. * nu2) * nu2 + 8. * nu1 * nu1 * (-26. + nu2 * (-53. + 105. * nu2)) + nu1 * (7. + 4. * nu2 * (-108. + nu2 * (-106. + 105. * nu2)))))) / 980980.;

                    if (SWITCH_index == 1) {

                    nu1 = -0.5 * etam_transfer[index_i];
                    nu2 = -0.5 * etam_transfer[index_l]; //GC!

                    //pnlpt->M12_oneline_4_vv_complex[count] = (pnlpt->M12_oneline_complex[count])*(8.*f*f*((-15.*(88. + 91.*f) + 29.*(55. + 91.*f)*nu2 + 448.*(11. + 15.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 4.*(2607. + 2569.*f)*nu2*nu2 - 4.*(5203. + 6167.*f)*nu2*nu2*nu2 + 32.*(396. + 511.*f)*nu2*nu2*nu2*nu2 - 224.*(11. + 15.*f)*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(231. + 315.*f + (4818. + 5642.*f)*nu2 - 16.*(561. + 707.*f)*nu2*nu2 + 336.*(11. + 15.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(-759. - 861.*f - 4.*(1485. + 1358.*f)*nu2 + 28.*(649. + 711.*f)*nu2*nu2 - 8.*(1859. + 2303.*f)*nu2*nu2*nu2 + 336.*(11. + 15.*f)*nu2*nu2*nu2*nu2) + nu1*(3531. + 3423.*f + 8.*(1089. + 644.*f)*nu2 - 8.*(6655. + 6349.*f)*nu2*nu2 + 8.*(8261. + 9177.*f)*nu2*nu2*nu2 - 128.*(242. + 301.*f)*nu2*nu2*nu2*nu2 + 448.*(11. + 15.*f)*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*(11. + 15.*f)*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 32.*nu1*nu1*nu1*nu1*(396. + 511.*f - 4.*(242. + 301.*f)*nu2 + 42.*(11. + 15.*f)*nu2*nu2) + 3.*(-5.*(88. + 91.*f) + (1177. + 1141.*f)*nu2 - 4.*(253. + 287.*f)*nu2*nu2 + 28.*(11. + 15.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-5203. - 6167.*f + 2.*(8261. + 9177.*f)*nu2 - 8.*(1859. + 2303.*f)*nu2*nu2 + 336.*(11. + 15.*f)*nu2*nu2*nu2) + nu1*(29.*(55. + 91.*f) + 8.*(1089. + 644.*f)*nu2 - 16.*(1485. + 1358.*f)*nu2*nu2 + 8.*(2409. + 2821.*f)*nu2*nu2*nu2 - 448.*(11. + 15.*f)*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(2607. + 2569.*f - 2.*(6655. + 6349.*f)*nu2 + 28.*(649. + 711.*f)*nu2*nu2 - 16.*(561. + 707.*f)*nu2*nu2*nu2 + 112.*(11. + 15.*f)*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (448.*(11. + 15.*f)*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 3.*(5.*(88. + 91.*f) - 1.*(1177. + 1141.*f)*nu2 + 4.*(253. + 287.*f)*nu2*nu2 - 28.*(11. + 15.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-21.*(11. + 15.*f) + (5126. + 6062.*f)*nu2 - 8.*(891. + 1099.*f)*nu2*nu2 + 224.*(11. + 15.*f)*nu2*nu2*nu2) + 4.*nu1*nu1*(759. + 861.*f - 4.*(1969. + 1960.*f)*nu2 + 44.*(289. + 315.*f)*nu2*nu2 - 8.*(891. + 1099.*f)*nu2*nu2*nu2 + 112.*(11. + 15.*f)*nu2*nu2*nu2*nu2) - 1.*nu1*(3531. + 3423.*f - 2552.*(8. + 7.*f)*nu2 + 16.*(1969. + 1960.*f)*nu2*nu2 - 8.*(2563. + 3031.*f)*nu2*nu2*nu2 + 448.*(11. + 15.*f)*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))/(385.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)));

                    pnlpt->M12_oneline_4_vv_complex[count] = 8. / 7. * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2[count]) + f * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv4_f3[count]);

                    //GC -> 13 of 14...

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_oneline_4_vv_complex_ortho[count] = 8. / 7. * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2_ortho[count]) + f * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv4_f3_ortho[count]);

                    //GC: ORTHOGONAL -- finish
                    }
                    count++;
                }
            }

            double complex *f13_4_vv;
            class_alloc(f13_4_vv, Nmax * sizeof(complex double), pnlpt->error_message);
            //    double *P_CTR0;
            //    class_alloc(P_CTR0,Nmax*sizeof(double),pnlpt->error_message);

            double complex *f22_4_vv;
            class_alloc(f22_4_vv, Nmax * sizeof(complex double), pnlpt->error_message);

            double complex *f12_4_vv;
            class_alloc(f12_4_vv, Nmax * sizeof(complex double), pnlpt->error_message);
            //GC!

            for (index_i = 0; index_i < Nmax + 1; index_i++)
            {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_4_vv_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * (3. * f * f * (-55. + 33. * nu1 + f * (-66. + 90. * nu1))) / 5390.;
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                }
                f13_4_vv[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_4_vv_oneline_complex, &inc);
                P13UV_4_vv[index_j] = -1. * Pbin[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav * (24. * f * f * (33. + 58. * f + 25. * f * f) / 1925.);
                P13_4_vv[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_4_vv[index_j] * Pbin[index_j]) + P13UV_4_vv[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                    if (SWITCH_index == 1) {

                    x_transfer[count] = cmsym_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]); //GC!
                    }
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_4_vv_complex, x, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_4_vv_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_4_vv[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {

                f12_4_vv[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_4_vv[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_4_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {

                P12_4_vv[index_j] = Tbin[index_j] *                                                                               //[FACTOR] *
                                    creal(cpow(kdisc[index_j], 3.) * f12_4_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                P1loop_4_vv[index_j] = Pbin[index_j] * (f * f * 8. / 35.) * 0. + (P13_4_vv[index_j] + P22_4_vv[index_j]);
                P_CTR_4[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j] * (f * f * 8. / 35.);
                Ptree_4_vv[index_j] = Pbin[index_j] * (f * f * 8. / 35.);
                //    printf("%le %le\n",kdisc[j],P22[j]);

                //printf("%.24f %.24f %.24f",kdisc[index_j],Tbin[index_j],P12_4_vv[index_j]);
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- start
                if (SWITCH_index == 1) {

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_4_vv_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);                        //GC!
                f12_4_vv[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                                 //GC!
                P12_4_vv_ortho[index_j] = Tbin[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_4_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //GC: ORTHOGONAL -- finish

                //        printf("%.20f %.20f %.20f",kdisc[index_j],Tbin[index_j],P12_4_vv_ortho[index_j]);
                //        printf("\n"); //GC!
            }

            // Computing P_{vd} contribution - Hexadecapole

            count = 0;
            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    //  nu12 =-0.5*etam[index_i]-0.5*etam[index_l];
                    nu1 = -0.5 * etam[index_i];
                    nu2 = -0.5 * etam[index_l];

                    nu12 = nu1 + nu2;

                    pnlpt->M22_oneline_4_vd_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * f * f * (14. * f * (26. + 13. * nu2 - 60. * nu2 * nu2 - 8. * nu2 * nu2 * nu2 + nu1 * nu1 * nu1 * (-8. + 60. * nu2) + 4. * nu1 * nu1 * (-15. + 4. * nu2 + 30. * nu2 * nu2) + nu1 * (13. - 104. * nu2 + 16. * nu2 * nu2 + 60. * nu2 * nu2 * nu2)) + 11. * (58. + 21. * nu2 + 112. * nu1 * nu1 * nu1 * nu2 - 106. * nu2 * nu2 + 2. * nu1 * nu1 * (-53. - 6. * nu2 + 112. * nu2 * nu2) + nu1 * (21. - 204. * nu2 - 12. * nu2 * nu2 + 112. * nu2 * nu2 * nu2))) / 2695.;

                    if (SWITCH_index == 1) {

                    nu1 = -0.5 * etam_transfer[index_i];
                    nu2 = -0.5 * etam_transfer[index_l]; //GC!

                    //pnlpt->M12_oneline_4_vd_complex[count] = (pnlpt->M12_oneline_complex[count])*(4.*f*f*((-15. + 140.*nu2 + 256.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 56.*nu2*nu2 - 624.*nu2*nu2*nu2 + 544.*nu2*nu2*nu2*nu2 - 128.*nu2*nu2*nu2*nu2*nu2 + 16.*nu1*nu1*nu1*(3. + 34.*nu2 - 88.*nu2*nu2 + 48.*nu2*nu2*nu2) + 8.*nu1*nu1*(-9. + 10.*nu2 + 172.*nu2*nu2 - 272.*nu2*nu2*nu2 + 96.*nu2*nu2*nu2*nu2) + 4.*nu1*(3. - 94.*nu2 - 20.*nu2*nu2 + 344.*nu2*nu2*nu2 - 288.*nu2*nu2*nu2*nu2 + 64.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (128.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 32.*nu1*nu1*nu1*nu1*(17. - 36.*nu2 + 24.*nu2*nu2) + 3.*(-5. + 4.*nu2 - 24.*nu2*nu2 + 16.*nu2*nu2*nu2) + 16.*nu1*nu1*nu1*(-39. + 86.*nu2 - 136.*nu2*nu2 + 48.*nu2*nu2*nu2) + 8.*nu1*nu1*(7. - 10.*nu2 + 172.*nu2*nu2 - 176.*nu2*nu2*nu2 + 32.*nu2*nu2*nu2*nu2) - 4.*nu1*(-35. + 94.*nu2 - 20.*nu2*nu2 - 136.*nu2*nu2*nu2 + 64.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (15. - 12.*nu2 + 256.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 72.*nu2*nu2 - 48.*nu2*nu2*nu2 + 16.*nu1*nu1*nu1*(-3. + 38.*nu2 - 64.*nu2*nu2 + 32.*nu2*nu2*nu2) + 8.*nu1*nu1*(9. - 26.*nu2 + 116.*nu2*nu2 - 128.*nu2*nu2*nu2 + 32.*nu2*nu2*nu2*nu2) - 4.*nu1*(3. + 42.*nu2 + 52.*nu2*nu2 - 152.*nu2*nu2*nu2 + 64.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))/(5.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)));

                    pnlpt->M12_oneline_4_vd_complex[count] = f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd4_f2[count]);

                    //GC -> 14 of 14...

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_oneline_4_vd_complex_ortho[count] = f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd4_f2_ortho[count]);
                    }
                    //GC: ORTHOGONAL -- finish

                    //GC!

                    count++;
                }
            }

            double complex *f13_4_vd;
            class_alloc(f13_4_vd, Nmax * sizeof(complex double), pnlpt->error_message);
            //    double *P_CTR0;
            //    class_alloc(P_CTR0,Nmax*sizeof(double),pnlpt->error_message);

            double complex *f22_4_vd;
            class_alloc(f22_4_vd, Nmax * sizeof(complex double), pnlpt->error_message);

            double complex *f12_4_vd;
            class_alloc(f12_4_vd, Nmax * sizeof(complex double), pnlpt->error_message); //GC!

            for (index_i = 0; index_i < Nmax + 1; index_i++)
            {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_4_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * 9. * (f * f * (1. + 2. * nu1)) / 245.;
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                }
                f13_4_vd[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_4_vd_oneline_complex, &inc);
                P13UV_4_vd[index_j] = -1. * Pbin[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav * 16. * f * f * (22. + 35. * f) / 1225.;
                P13_4_vd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_4_vd[index_j] * Pbin[index_j]) + P13UV_4_vd[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                    
                    if (SWITCH_index == 1) {

                    x_transfer[count] = cmsym_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]); //GC!
                    }
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_4_vd_complex, x, &inc, &beta, y, &inc);

                if (SWITCH_index == 1) {

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_4_vd_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_4_vd[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                
                if (SWITCH_index == 1) {

                f12_4_vd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_4_vd[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_4_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                
                if (SWITCH_index == 1) {

                P12_4_vd[index_j] = Tbin[index_j] *                                                                               //[FACTOR] *
                                    creal(cpow(kdisc[index_j], 3.) * f12_4_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC! Careful to not forget "Tbin[index_j] * [FACTOR]" here...
                    
                }
                
                P1loop_4_vd[index_j] = (P13_4_vd[index_j] + P22_4_vd[index_j]);
                //    P_CTR[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j];
                //printf("%le %le\n",kdisc[j],P22[j]);

                //printf("%.24f %.24f %.24f",kdisc[index_j],Tbin[index_j],P12_4_vd[index_j]);
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- start

                
                if (SWITCH_index == 1) {

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_4_vd_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);                        //GC!
                f12_4_vd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                                 //GC!
                P12_4_vd_ortho[index_j] = Tbin[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_4_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC! Careful to not forget
                }
                //GC: ORTHOGONAL -- finish

                //printf("%.20f %.20f %.20f",kdisc[index_j],Tbin[index_j],P12_4_vd_ortho[index_j]);
                //printf("\n"); //GC!
            }

            // Computing P_{dd} contribution - Hexadecapole

            //GC -> I do not have anything here... I feel that I have to put something, otherwise I am screwed if I do not go into the AP case... I will see what happens when I check this...

            count = 0;
            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam[index_i];
                    nu2 = -0.5 * etam[index_l];

                    nu12 = nu1 + nu2;

                    pnlpt->M22_oneline_4_dd_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * f * f * (2. * nu1 - 1.) * (2. * nu2 - 1.) * (1. + nu12) * (2. + nu12) / 35.;

                    count++;
                }
            }

            double complex *f22_4_dd;
            class_alloc(f22_4_dd, Nmax * sizeof(complex double), pnlpt->error_message);

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym[count] * cpow(kdisc[index_j], etam[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_4_dd_complex, x, &inc, &beta, y, &inc);
                f22_4_dd[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                P22_4_dd[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_4_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                P1loop_4_dd[index_j] = P22_4_dd[index_j];

                //printf("%.24f %.24f %.24f",kdisc[index_j],Tbin[index_j],P12_4_dd[index_j]); //GC -> it already contains all zeros!!! Great...
                //printf("\n"); //GC!
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

            //GC! Notice that some of these are not generated so I do not need to allocate them or deallocate them...
            
            //GC - SWITCH ~> I allocated these, so it is ok to deallocate them! Notice that likely it would have been enough to allogate one f for everything... But the original code was written in this way...

            free(f12_0_vv);
            free(f12_0_vd);
            free(f12_0_dd);
            free(f12_2_vv);
            free(f12_2_vd);
            free(f12_2_dd);
            free(f12_4_vv);
            free(f12_4_vd);
            //free(f12_4_dd); //GC: here it is correct that it is NOT generated... But then when I have AP, the P12_4_dd is surely there... So I do not need to free this f12, but P12 yes, because in ANY CASE it is generated by AP. But ACTUALLY see that I was smart and I had allocated it, and discussed already, around line 3273...

            //GC -> I should have everything so far ~> need to deallocate some stuff but ok...

        } // end of IR resummation expression

        if (irindex == 1)
        {
            // printf("Computing IR resummed spectra in redshift space...\n");

            //Computing FFT for wiggly and non-wiggly parts

            double *input_real_nw;
            class_alloc(input_real_nw, (Nmax) * sizeof(double), pnlpt->error_message);
            double *input_imag_nw;
            class_alloc(input_imag_nw, (Nmax) * sizeof(double), pnlpt->error_message);
            double *output_real_nw;
            class_alloc(output_real_nw, (Nmax) * sizeof(double), pnlpt->error_message);
            double *output_imag_nw;
            class_alloc(output_imag_nw, (Nmax) * sizeof(double), pnlpt->error_message);

            double *input_real_w;
            class_alloc(input_real_w, (Nmax) * sizeof(double), pnlpt->error_message);
            double *input_imag_w;
            class_alloc(input_imag_w, (Nmax) * sizeof(double), pnlpt->error_message);
            double *output_real_w;
            class_alloc(output_real_w, (Nmax) * sizeof(double), pnlpt->error_message);
            double *output_imag_w;
            class_alloc(output_imag_w, (Nmax) * sizeof(double), pnlpt->error_message);

            double *input_real_nw_transfer;
            class_alloc(input_real_nw_transfer, (Nmax) * sizeof(double), pnlpt->error_message);
            double *input_imag_nw_transfer;
            class_alloc(input_imag_nw_transfer, (Nmax) * sizeof(double), pnlpt->error_message);
            double *output_real_nw_transfer;
            class_alloc(output_real_nw_transfer, (Nmax) * sizeof(double), pnlpt->error_message);
            double *output_imag_nw_transfer;
            class_alloc(output_imag_nw_transfer, (Nmax) * sizeof(double), pnlpt->error_message);

            double *input_real_w_transfer;
            class_alloc(input_real_w_transfer, (Nmax) * sizeof(double), pnlpt->error_message);
            double *input_imag_w_transfer;
            class_alloc(input_imag_w_transfer, (Nmax) * sizeof(double), pnlpt->error_message);
            double *output_real_w_transfer;
            class_alloc(output_real_w_transfer, (Nmax) * sizeof(double), pnlpt->error_message);
            double *output_imag_w_transfer;
            class_alloc(output_imag_w_transfer, (Nmax) * sizeof(double), pnlpt->error_message);

            //GC!

            index_kd = 0;
            for (index_kd = 0; index_kd < Nmax; index_kd++)
            {
                input_real_nw[index_kd] = Pnw[index_kd] * exp(-1. * index_kd * b * Delta);
                input_imag_nw[index_kd] = 0.;
                input_real_w[index_kd] = Pw[index_kd] * exp(-1. * index_kd * b * Delta);
                input_imag_w[index_kd] = 0.;
                //GC!
                input_real_nw_transfer[index_kd] = Tnw[index_kd] * exp(-1. * index_kd * b_transfer * Delta); //GC: CAREFUL about b_transfer!!!
                input_imag_nw_transfer[index_kd] = 0.;
                input_real_w_transfer[index_kd] = Tw[index_kd] * exp(-1. * index_kd * b_transfer * Delta); //GC: CAREFUL about b_transfer!!!
                input_imag_w_transfer[index_kd] = 0.;
            }

            FFT(input_real_nw, input_imag_nw, output_real_nw, output_imag_nw, Nmax, stepsize);

            FFT(input_real_nw_transfer, input_imag_nw_transfer, output_real_nw_transfer, output_imag_nw_transfer, Nmax, stepsize); //GC!

            double complex *cmsym_nw;
            class_alloc(cmsym_nw, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

            double complex *cmsym_nw_transfer;
            class_alloc(cmsym_nw_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);
            //GC!

            index_c = 0;
            for (index_c = 0; index_c < Nmax + 1; index_c++)
            {
                if (index_c < Nmax / 2)
                {
                    cmsym_nw[index_c] = cpow(kmin, -etam[index_c]) * (output_real_nw[Nmax / 2 - index_c] - _Complex_I * output_imag_nw[Nmax / 2 - index_c]) / Nmaxd;
                    cmsym_nw_transfer[index_c] = cpow(kmin, -etam_transfer[index_c]) * (output_real_nw_transfer[Nmax / 2 - index_c] - _Complex_I * output_imag_nw_transfer[Nmax / 2 - index_c]) / Nmaxd; //GC! CAREFUL about etam_transfer!!!
                }
                else
                {
                    cmsym_nw[index_c] = cpow(kmin, -etam[index_c]) * (output_real_nw[index_c - Nmax / 2] + _Complex_I * output_imag_nw[index_c - Nmax / 2]) / Nmaxd;
                    cmsym_nw_transfer[index_c] = cpow(kmin, -etam_transfer[index_c]) * (output_real_nw_transfer[index_c - Nmax / 2] + _Complex_I * output_imag_nw_transfer[index_c - Nmax / 2]) / Nmaxd; //GC! CAREFUL about etam_transfer!!!
                }
            }

            cmsym_nw[0] = cmsym_nw[0] / 2.;
            cmsym_nw[Nmax] = cmsym_nw[Nmax] / 2.;

            cmsym_nw_transfer[0] = cmsym_nw_transfer[0] / 2.;
            cmsym_nw_transfer[Nmax] = cmsym_nw_transfer[Nmax] / 2.; //GC!

            FFT(input_real_w, input_imag_w, output_real_w, output_imag_w, Nmax, stepsize);

            FFT(input_real_w_transfer, input_imag_w_transfer, output_real_w_transfer, output_imag_w_transfer, Nmax, stepsize); //GC!

            double complex *cmsym_w;
            class_alloc(cmsym_w, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

            double complex *cmsym_w_transfer;
            class_alloc(cmsym_w_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message); //GC!

            index_c = 0;
            for (index_c = 0; index_c < Nmax + 1; index_c++)
            {
                if (index_c < Nmax / 2)
                {
                    cmsym_w[index_c] = cpow(kmin, -etam[index_c]) * (output_real_w[Nmax / 2 - index_c] - _Complex_I * output_imag_w[Nmax / 2 - index_c]) / Nmaxd;
                    cmsym_w_transfer[index_c] = cpow(kmin, -etam_transfer[index_c]) * (output_real_w_transfer[Nmax / 2 - index_c] - _Complex_I * output_imag_w_transfer[Nmax / 2 - index_c]) / Nmaxd; //GC! NEED TO FIX etam_transfer HERE!!! They are defined around line 2760, and they benefit from the fact that for ALL my calculations I have a bias of -0.8, apart from the one for \delta^2, which is in the section below anyway... Be very careful to fix this etam_transfer around here...
                }
                else
                {
                    cmsym_w[index_c] = cpow(kmin, -etam[index_c]) * (output_real_w[index_c - Nmax / 2] + _Complex_I * output_imag_w[index_c - Nmax / 2]) / Nmaxd;
                    cmsym_w_transfer[index_c] = cpow(kmin, -etam_transfer[index_c]) * (output_real_w_transfer[index_c - Nmax / 2] + _Complex_I * output_imag_w_transfer[index_c - Nmax / 2]) / Nmaxd; //GC!
                }
            }

            cmsym_w[0] = cmsym_w[0] / 2.;
            cmsym_w[Nmax] = cmsym_w[Nmax] / 2.;

            cmsym_w_transfer[0] = cmsym_w_transfer[0] / 2.;
            cmsym_w_transfer[Nmax] = cmsym_w_transfer[Nmax] / 2.; //GC!

            free(input_real_nw);
            free(input_imag_nw);
            free(output_real_nw);
            free(output_imag_nw);

            free(input_real_w);
            free(input_imag_w);
            free(output_real_w);
            free(output_imag_w);

            //GC!

            free(input_real_nw_transfer);
            free(input_imag_nw_transfer);
            free(output_real_nw_transfer);
            free(output_imag_nw_transfer);

            free(input_real_w_transfer);
            free(input_imag_w_transfer);
            free(output_real_w_transfer);
            free(output_imag_w_transfer);

            //GC!

            //GC -> will need to do the same...
            
            //GC - SWITCH -> also here I keep the decomposition above, I just switch off the costly calculations...

            // Matrix multiplication

            double *P13_mu0_dd;
            double *P13UV_mu0_dd;
            class_alloc(P13_mu0_dd, Nmax * sizeof(double), pnlpt->error_message);
            class_alloc(P13UV_mu0_dd, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu0_dd;
            class_alloc(P22_mu0_dd, Nmax * sizeof(double), pnlpt->error_message);

            double *P13_mu2_dd;
            double *P13UV_mu2_dd;
            class_alloc(P13_mu2_dd, Nmax * sizeof(double), pnlpt->error_message);
            class_alloc(P13UV_mu2_dd, Nmax * sizeof(double), pnlpt->error_message);

            double *P13_mu2_vd;
            double *P13UV_mu2_vd;
            class_alloc(P13_mu2_vd, Nmax * sizeof(double), pnlpt->error_message);
            class_alloc(P13UV_mu2_vd, Nmax * sizeof(double), pnlpt->error_message);

            double *P13_mu4_vd;
            double *P13UV_mu4_vd;
            class_alloc(P13_mu4_vd, Nmax * sizeof(double), pnlpt->error_message);
            class_alloc(P13UV_mu4_vd, Nmax * sizeof(double), pnlpt->error_message);

            double *P13_mu4_vv;
            double *P13UV_mu4_vv;
            class_alloc(P13_mu4_vv, Nmax * sizeof(double), pnlpt->error_message);
            class_alloc(P13UV_mu4_vv, Nmax * sizeof(double), pnlpt->error_message);

            double *P13_mu6;
            double *P13UV_mu6;
            class_alloc(P13_mu6, Nmax * sizeof(double), pnlpt->error_message);
            class_alloc(P13UV_mu6, Nmax * sizeof(double), pnlpt->error_message);

            double *P22_mu2_vd;
            class_alloc(P22_mu2_vd, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu2_dd;
            class_alloc(P22_mu2_dd, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu4_vv;
            class_alloc(P22_mu4_vv, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu4_vd;
            class_alloc(P22_mu4_vd, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu4_dd;
            class_alloc(P22_mu4_dd, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu6_vv;
            class_alloc(P22_mu6_vv, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu6_vd;
            class_alloc(P22_mu6_vd, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu8;
            class_alloc(P22_mu8, Nmax * sizeof(double), pnlpt->error_message);

            //GC!

            //double *P12_mu0_dd;
            //class_alloc(P12_mu0_dd,Nmax*sizeof(double),pnlpt->error_message);

            //GC -> actually, this is NOT needed -> this is literally what I would compute via 2*F2*shape, and that I call pnlpt->M12_oneline_complex. This further confirms that M12 NOT BASIC is the pure matter piece (WITH THE FACTOR OF 2 FROM PERMUTATIONS), and that the (nu1,nu2) tricks are defined with respect to it... See below when I look at biased tracers for details on the "basic"...

            double *P12_mu0_dd;
            class_alloc(P12_mu0_dd, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu2_vd;
            class_alloc(P12_mu2_vd, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu2_dd;
            class_alloc(P12_mu2_dd, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu4_vv;
            class_alloc(P12_mu4_vv, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu4_vd;
            class_alloc(P12_mu4_vd, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu6_vv;
            class_alloc(P12_mu6_vv, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- start

            double *P12_mu0_dd_ortho;
            class_alloc(P12_mu0_dd_ortho, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu2_vd_ortho;
            class_alloc(P12_mu2_vd_ortho, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu2_dd_ortho;
            class_alloc(P12_mu2_dd_ortho, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu4_vv_ortho;
            class_alloc(P12_mu4_vv_ortho, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu4_vd_ortho;
            class_alloc(P12_mu4_vd_ortho, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu6_vv_ortho;
            class_alloc(P12_mu6_vv_ortho, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- finish

            //GC -> will need to allocate up to \mu^6 again after Z1 -> again for now no FoG... Here they are the powers of \mu... Notice that here again he adds f in the matrix. So I will have to do that as well. Nothing complicated... The confusing part is: why does he extract the coefficients of the wiggly and non-wiggly part? This must be related to 2.32 of their paper, and 2.33. 2.33 first line I get, but in the second line I need nw + w, not only w... And Pw is definitely ONLY the wiggly part... I understand that in 2.33 you don't want the Pw with the damping, so I cannot use Pbin. But then at least I should have some cross term... HE HAS IT!!! It is a factor of 2 that multiplies the coefficients of the wiggly part (see "x_w[count]= 2.*cmsym_w[count]* cpow(kdisc[index_j], etam[count]);") => then he has 2.33 expanded "functionally" at leading order in wiggly... So here I proceed again as a bullet... THERE IS STILL ONE TRICKY THING. First, for P22 this is perfect. But for P13 it is not. Indeed, he computes the actual integral with the wiggly PS, but on the outside he multiplies it by Pnw. Now -> this is NOT the full result. How does he fix this? He defines Pnw_ap_out and Pw_ap_out that both start at zero. Now, the idea should be, I THINK, that in P13 the suppression of the PS outside of the integral can be obtained by replacing, in the full non-wiggly P13, the power spectrum outside with the exponential and the wiggly... So if I have Pnw, I multiply it by 1 + (Pw/Pnw)*Exp. This is what he does below, but AFTER doing the AP effect... I do not know why it is fine to do it after the AP effect, AND if it is correct for me to do the same with P12 (since P12 has an integral as P22, but it also has a "power spectrum" outside as P13). But for now I will do the same... Essentially, the point is that he does the exponential modified by RSD stuff ONLY on the outside... So now the first thing is derive the matrices for the powers of \mu, not the multipoles...

            count = 0;
            index_l = 0;
            index_i = 0;
            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam[index_i];
                    nu2 = -0.5 * etam[index_l];
                    nu12 = nu1 + nu2;

                    pnlpt->M22_oneline_mu2_vd_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * (-1.) * f * (7. * f * (-1. + 2. * nu1) * (-1. + 2. * nu2) * (6. + 7. * nu12) - 4. * (46. + 13. * nu2 + 98. * nu1 * nu1 * nu1 * nu2 - 63. * nu2 * nu2 + 7. * nu1 * nu1 * (-9. - 10. * nu2 + 28. * nu2 * nu2) + nu1 * (13. - 138. * nu2 - 70. * nu2 * nu2 + 98. * nu2 * nu2 * nu2))) / 392.;

                    pnlpt->M22_oneline_mu2_dd_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * f * (7. * f * (2. + 2. * nu1 * nu1 * nu1 - nu2 - nu2 * nu2 + 2. * nu2 * nu2 * nu2 - nu1 * nu1 * (1. + 2. * nu2) - nu1 * (1. + 2. * nu2 + 2. * nu2 * nu2)) + 4. * (10. - nu2 + 14. * nu1 * nu1 * nu1 * nu2 - 17. * nu2 * nu2 + nu1 * nu1 * (-17. + 6. * nu2 + 28. * nu2 * nu2) + nu1 * (-1. - 22. * nu2 + 6. * nu2 * nu2 + 14. * nu2 * nu2 * nu2))) / 56.;

                    pnlpt->M22_oneline_mu4_vv_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * f * f * (147. * f * f * (-1. + 2. * nu1) * (-1. + 2. * nu2) - 28. * f * (-1. + 2. * nu1) * (-1. + 2. * nu2) * (-2. + 7. * nu12) + 8. * (50. - 9. * nu2 + 98 * nu1 * nu1 * nu1 * nu2 - 35. * nu2 * nu2 + 7. * nu1 * nu1 * (-5. - 18. * nu2 + 28 * nu2 * nu2) + nu1 * (-9. - 66. * nu2 - 126 * nu2 * nu2 + 98. * nu2 * nu2 * nu2))) / 1568.;

                    pnlpt->M22_oneline_mu4_vd_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * f * f * (58. + 21. * nu2 + 112. * nu1 * nu1 * nu1 * nu2 - 106. * nu2 * nu2 + 2. * nu1 * nu1 * (-53. - 6. * nu2 + 112. * nu2 * nu2) + 7. * f * (2. + nu1 + 4. * nu1 * nu1 * nu1 + nu2 - 8. * nu1 * nu2 - 8. * nu1 * nu1 * nu2 - 8. * nu1 * nu2 * nu2 + 4. * nu2 * nu2 * nu2) + nu1 * (21. - 204. * nu2 - 12. * nu2 * nu2 + 112. * nu2 * nu2 * nu2)) / 56.;

                    pnlpt->M22_oneline_mu4_dd_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * f * f * (2. * nu1 - 1.) * (2. * nu2 - 1.) * (2. + nu1 * nu1 + 3. * nu2 + nu2 * nu2 + nu1 * (3. + 2. * nu2)) / 8.;

                    pnlpt->M22_oneline_mu6_vv_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * f * f * f * (7. * f * (1. + 4. * nu1 * nu1 * nu1 + nu1 * nu1 * (2. - 12. * nu2) + 2. * nu2 + 2. * nu2 * nu2 + 4. * nu2 * nu2 * nu2 - 2. * nu1 * (-1. + 4. * nu2 + 6. * nu2 * nu2)) + 2. * (26. + 9. * nu2 + 56. * nu1 * nu1 * nu1 * nu2 - 38. * nu2 * nu2 + 2. * nu1 * nu1 * (-19. - 18. * nu2 + 56. * nu2 * nu2) + nu1 * (9. - 84. * nu2 - 36. * nu2 * nu2 + 56. * nu2 * nu2 * nu2))) / 112.;

                    pnlpt->M22_oneline_mu6_vd_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * f * f * f * (2. * nu1 - 1.) * (2. * nu2 - 1.) * (2. + 2. * nu1 * nu1 + 5. * nu2 + 2. * nu2 * nu2 + nu1 * (5. + 4. * nu2)) / 8.;

                    pnlpt->M22_oneline_mu8_complex[count] = (pnlpt->M22_oneline_complex[count] * 196. / (98. * nu1 * nu2 * nu12 * nu12 - 91. * nu12 * nu12 + 36. * nu1 * nu2 - 14. * nu1 * nu2 * nu12 + 3. * nu12 + 58.)) * f * f * f * f * (2. * nu1 - 1.) * (2. * nu2 - 1.) * (3. + 4. * nu1 * nu1 + 8. * nu2 + 4. * nu2 * nu2 + 8. * nu1 * (1. + nu2)) / 32.;

                    //GC!

                    
                    if (SWITCH_index == 1) {
                    
                    nu1 = -0.5 * etam_transfer[index_i];
                    nu2 = -0.5 * etam_transfer[index_l]; //GC!

                    //pnlpt->M12_oneline_mu2_vd_complex[count] = (pnlpt->M12_oneline_complex[count])*(-1.*(f*((7.*f*(15. + 16.*nu2 + 128.*nu1*nu1*nu1*(-1. + nu2)*nu2 - 136.*nu2*nu2 + 128.*nu2*nu2*nu2 - 32.*nu2*nu2*nu2*nu2 + 8.*nu1*nu1*(3. + 46.*nu2 - 84.*nu2*nu2 + 32.*nu2*nu2*nu2) + 8.*nu1*(-6. - 29.*nu2 + 94.*nu2*nu2 - 72.*nu2*nu2*nu2 + 16.*nu2*nu2*nu2*nu2)) - 4.*(-90. + 177.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 676.*nu2*nu2 - 1636.*nu2*nu2*nu2 + 1088.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 374.*nu2 - 752.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-57. - 356.*nu2 + 1316.*nu2*nu2 - 1224.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(225. + 328.*nu2 - 3336.*nu2*nu2 + 4856.*nu2*nu2*nu2 - 2560.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2)))*csin(2.*nu1*M_PI) + (7.*f*(32.*nu1*nu1*nu1*nu1*(-1. + 4.*nu2) + 64.*nu1*nu1*nu1*(2. - 9.*nu2 + 4.*nu2*nu2) + 3.*(5. - 16.*nu2 + 8.*nu2*nu2) + 8.*nu1*nu1*(-17. + 94.*nu2 - 84.*nu2*nu2 + 16.*nu2*nu2*nu2) - 8.*nu1*(-2. + 29.*nu2 - 46.*nu2*nu2 + 16.*nu2*nu2*nu2)) - 4.*(224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(17. - 40.*nu2 + 21.*nu2*nu2) + 3.*(-30. + 75.*nu2 - 76.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-409. + 1214.*nu2 - 1224.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(177. + 328.*nu2 - 1424.*nu2*nu2 + 1496.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(169. - 834.*nu2 + 1316.*nu2*nu2 - 752.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2)))*csin(2.*nu2*M_PI) + (7.*f*(128.*nu1*nu1*nu1*(-1. + nu2)*nu2 - 3.*(5. - 16.*nu2 + 8.*nu2*nu2) + 8.*nu1*nu1*(-3. + 50.*nu2 - 60.*nu2*nu2 + 16.*nu2*nu2*nu2) - 8.*nu1*(-6. + 43.*nu2 - 50.*nu2*nu2 + 16.*nu2*nu2*nu2)) - 4.*(90. - 225.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 228.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 402.*nu2 - 584.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-225. + 1168.*nu2 - 2064.*nu2*nu2 + 1608.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(57. - 516.*nu2 + 916.*nu2*nu2 - 584.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2)))*csin(2.*(nu1 + nu2)*M_PI)))/(2.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI))));

                    pnlpt->M12_oneline_mu2_vd_complex[count] = f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f1[count]) + f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2[count]);

                    //GC -> 2 of 7...

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_oneline_mu2_vd_complex_ortho[count] = f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f1_ortho[count]) + f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2_ortho[count]);

                    //GC: ORTHOGONAL -- finish

                    //pnlpt->M12_oneline_mu2_dd_complex[count] = (pnlpt->M12_oneline_complex[count])*(7.*f*((64.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 4.*nu1*nu1*nu1*(3. + 26.*nu2 - 80.*nu2*nu2 + 48.*nu2*nu2*nu2) + nu2*(39. - 20.*nu2 - 124.*nu2*nu2 + 128.*nu2*nu2*nu2 - 32.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(-3. + 28.*nu2 + 44.*nu2*nu2 - 120.*nu2*nu2*nu2 + 48.*nu2*nu2*nu2*nu2) + nu1*(-9. - 152.*nu2 + 168.*nu2*nu2 + 200.*nu2*nu2*nu2 - 256.*nu2*nu2*nu2*nu2 + 64.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (32.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(2. - 4.*nu2 + 3.*nu2*nu2) + 3.*nu2*(-3. - 4.*nu2 + 4.*nu2*nu2) + 4.*nu1*nu1*nu1*(-31. + 50.*nu2 - 120.*nu2*nu2 + 48.*nu2*nu2*nu2) + nu1*(39. - 152.*nu2 + 112.*nu2*nu2 + 104.*nu2*nu2*nu2 - 64.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(-5. + 42.*nu2 + 44.*nu2*nu2 - 80.*nu2*nu2*nu2 + 16.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (64.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 3.*nu2*(3. + 4.*nu2 - 4.*nu2*nu2) + 4.*nu1*nu1*nu1*(-3. + 30.*nu2 - 56.*nu2*nu2 + 32.*nu2*nu2*nu2) + nu1*(9. - 128.*nu2 + 48.*nu2*nu2 + 120.*nu2*nu2*nu2 - 64.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(3. + 12.*nu2 + 28.*nu2*nu2 - 56.*nu2*nu2*nu2 + 16.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))/((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI));

                    pnlpt->M12_oneline_mu2_dd_complex[count] = f * (pnlpt->M12_oneline_complex_matter_mu_powers_dd2_f1[count]);

                    //GC -> 3 of 7...

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_oneline_mu2_dd_complex_ortho[count] = f * (pnlpt->M12_oneline_complex_matter_mu_powers_dd2_f1_ortho[count]);

                    //GC: ORTHOGONAL -- finish

                    //pnlpt->M12_oneline_mu4_vv_complex[count] = (pnlpt->M12_oneline_complex[count])*(-1.*(f*f*((7.*f*(15. + 16.*nu2 + 128.*nu1*nu1*nu1*(-1. + nu2)*nu2 - 136.*nu2*nu2 + 128.*nu2*nu2*nu2 - 32.*nu2*nu2*nu2*nu2 + 8.*nu1*nu1*(3. + 46.*nu2 - 84.*nu2*nu2 + 32.*nu2*nu2*nu2) + 8.*nu1*(-6. - 29.*nu2 + 94.*nu2*nu2 - 72.*nu2*nu2*nu2 + 16.*nu2*nu2*nu2*nu2)) - 2.*(-120. + 145.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 948.*nu2*nu2 - 1892.*nu2*nu2*nu2 + 1152.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 12.*nu1*nu1*nu1*(7. + 146.*nu2 - 272.*nu2*nu2 + 112.*nu2*nu2*nu2) + 4.*nu1*nu1*(-69. - 540.*nu2 + 1652.*nu2*nu2 - 1352.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(321. + 792.*nu2 - 4840.*nu2*nu2 + 6008.*nu2*nu2*nu2 - 2816.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2)))*csin(2.*nu1*M_PI) + (7.*f*(32.*nu1*nu1*nu1*nu1*(-1. + 4.*nu2) + 64.*nu1*nu1*nu1*(2. - 9.*nu2 + 4.*nu2*nu2) + 3.*(5. - 16.*nu2 + 8.*nu2*nu2) + 8.*nu1*nu1*(-17. + 94.*nu2 - 84.*nu2*nu2 + 16.*nu2*nu2*nu2) - 8.*nu1*(-2. + 29.*nu2 - 46.*nu2*nu2 + 16.*nu2*nu2*nu2)) - 2.*(224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(18. - 44.*nu2 + 21.*nu2*nu2) + 3.*(-40. + 107.*nu2 - 92.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-473. + 1502.*nu2 - 1352.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(145. + 792.*nu2 - 2160.*nu2*nu2 + 1752.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(237. - 1210.*nu2 + 1652.*nu2*nu2 - 816.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2)))*csin(2.*nu2*M_PI) + (7.*f*(128.*nu1*nu1*nu1*(-1. + nu2)*nu2 - 3.*(5. - 16.*nu2 + 8.*nu2*nu2) + 8.*nu1*nu1*(-3. + 50.*nu2 - 60.*nu2*nu2 + 16.*nu2*nu2*nu2) - 8.*nu1*(-6. + 43.*nu2 - 50.*nu2*nu2 + 16.*nu2*nu2*nu2)) + 2.*(-448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 3.*(-40. + 107.*nu2 - 92.*nu2*nu2 + 28.*nu2*nu2*nu2) - 4.*nu1*nu1*nu1*(-21. + 466.*nu2 - 648.*nu2*nu2 + 224.*nu2*nu2*nu2) - 4.*nu1*nu1*(69. - 716.*nu2 + 1156.*nu2*nu2 - 648.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2) + nu1*(321. - 1856.*nu2 + 2864.*nu2*nu2 - 1864.*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2)))*csin(2.*(nu1 + nu2)*M_PI)))/(2.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI))));

                    pnlpt->M12_oneline_mu4_vv_complex[count] = f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vv4_f2[count]) + f * f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2[count]);

                    //GC -> 5 of 7...

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_oneline_mu4_vv_complex_ortho[count] = f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vv4_f2_ortho[count]) + f * f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2_ortho[count]);

                    //GC: ORTHOGONAL -- finish

                    //pnlpt->M12_oneline_mu4_vd_complex[count] = (pnlpt->M12_oneline_complex[count])*(7.*f*f*((-15. + 140.*nu2 + 256.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 56.*nu2*nu2 - 624.*nu2*nu2*nu2 + 544.*nu2*nu2*nu2*nu2 - 128.*nu2*nu2*nu2*nu2*nu2 + 16.*nu1*nu1*nu1*(3. + 34.*nu2 - 88.*nu2*nu2 + 48.*nu2*nu2*nu2) + 8.*nu1*nu1*(-9. + 10.*nu2 + 172.*nu2*nu2 - 272.*nu2*nu2*nu2 + 96.*nu2*nu2*nu2*nu2) + 4.*nu1*(3. - 94.*nu2 - 20.*nu2*nu2 + 344.*nu2*nu2*nu2 - 288.*nu2*nu2*nu2*nu2 + 64.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (128.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 32.*nu1*nu1*nu1*nu1*(17. - 36.*nu2 + 24.*nu2*nu2) + 3.*(-5. + 4.*nu2 - 24.*nu2*nu2 + 16.*nu2*nu2*nu2) + 16.*nu1*nu1*nu1*(-39. + 86.*nu2 - 136.*nu2*nu2 + 48.*nu2*nu2*nu2) + 8.*nu1*nu1*(7. - 10.*nu2 + 172.*nu2*nu2 - 176.*nu2*nu2*nu2 + 32.*nu2*nu2*nu2*nu2) - 4.*nu1*(-35. + 94.*nu2 - 20.*nu2*nu2 - 136.*nu2*nu2*nu2 + 64.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (15. - 12.*nu2 + 256.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 72.*nu2*nu2 - 48.*nu2*nu2*nu2 + 16.*nu1*nu1*nu1*(-3. + 38.*nu2 - 64.*nu2*nu2 + 32.*nu2*nu2*nu2) + 8.*nu1*nu1*(9. - 26.*nu2 + 116.*nu2*nu2 - 128.*nu2*nu2*nu2 + 32.*nu2*nu2*nu2*nu2) - 4.*nu1*(3. + 42.*nu2 + 52.*nu2*nu2 - 152.*nu2*nu2*nu2 + 64.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))/(2.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)));

                    pnlpt->M12_oneline_mu4_vd_complex[count] = f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd4_f2[count]);

                    //GC -> 6 of 7...

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_oneline_mu4_vd_complex_ortho[count] = f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd4_f2_ortho[count]);

                    //GC: ORTHOGONAL -- finish

                    //pnlpt->M12_oneline_mu6_vv_complex[count] = (pnlpt->M12_oneline_complex[count])*(7.*f*f*f*((-15. + 62.*nu2 + 128.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 96.*nu2*nu2 - 376.*nu2*nu2*nu2 + 288.*nu2*nu2*nu2*nu2 - 64.*nu2*nu2*nu2*nu2*nu2 + 24.*nu1*nu1*nu1*(1. + 14.*nu2 - 32.*nu2*nu2 + 16.*nu2*nu2*nu2) + 16.*nu1*nu1*(-3. - 9.*nu2 + 64.*nu2*nu2 - 76.*nu2*nu2*nu2 + 24.*nu2*nu2*nu2*nu2) + 2.*nu1*(15. - 36.*nu2 - 208.*nu2*nu2 + 488.*nu2*nu2*nu2 - 320.*nu2*nu2*nu2*nu2 + 64.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (64.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 32.*nu1*nu1*nu1*nu1*(9. - 20.*nu2 + 12.*nu2*nu2) + 3.*(-5. + 10.*nu2 - 16.*nu2*nu2 + 8.*nu2*nu2*nu2) + 8.*nu1*nu1*nu1*(-47. + 122.*nu2 - 152.*nu2*nu2 + 48.*nu2*nu2*nu2) + 32.*nu1*nu1*(3. - 13.*nu2 + 32.*nu2*nu2 - 24.*nu2*nu2*nu2 + 4.*nu2*nu2*nu2*nu2) - 2.*nu1*(-31. + 36.*nu2 + 72.*nu2*nu2 - 168.*nu2*nu2*nu2 + 64.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (15. - 30.*nu2 + 128.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 48.*nu2*nu2 - 24.*nu2*nu2*nu2 + 8.*nu1*nu1*nu1*(-3. + 46.*nu2 - 72.*nu2*nu2 + 32.*nu2*nu2*nu2) + 16.*nu1*nu1*(3. - 19.*nu2 + 44.*nu2*nu2 - 36.*nu2*nu2*nu2 + 8.*nu2*nu2*nu2*nu2) - 2.*nu1*(15. - 44.*nu2 + 152.*nu2*nu2 - 184.*nu2*nu2*nu2 + 64.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))/(2.*((-60. + 209.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 404.*nu2*nu2 - 1380.*nu2*nu2*nu2 + 1024.*nu2*nu2*nu2*nu2 - 224.*nu2*nu2*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(21. + 310.*nu2 - 688.*nu2*nu2 + 336.*nu2*nu2*nu2) + 4.*nu1*nu1*(-45. - 172.*nu2 + 980.*nu2*nu2 - 1096.*nu2*nu2*nu2 + 336.*nu2*nu2*nu2*nu2) + nu1*(129. - 136.*nu2 - 1832.*nu2*nu2 + 3704.*nu2*nu2*nu2 - 2304.*nu2*nu2*nu2*nu2 + 448.*nu2*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (224.*nu1*nu1*nu1*nu1*nu1*(-1. + 2.*nu2) + 64.*nu1*nu1*nu1*nu1*(16. - 36.*nu2 + 21.*nu2*nu2) + 3.*(-20. + 43.*nu2 - 60.*nu2*nu2 + 28.*nu2*nu2*nu2) + 4.*nu1*nu1*nu1*(-345. + 926.*nu2 - 1096.*nu2*nu2 + 336.*nu2*nu2*nu2) + nu1*(209. - 136.*nu2 - 688.*nu2*nu2 + 1240.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(101. - 458.*nu2 + 980.*nu2*nu2 - 688.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (60. - 129.*nu2 + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 + 180.*nu2*nu2 - 84.*nu2*nu2*nu2 + 4.*nu1*nu1*nu1*(-21. + 338.*nu2 - 520.*nu2*nu2 + 224.*nu2*nu2*nu2) + nu1*(-129. + 480.*nu2 - 1264.*nu2*nu2 + 1352.*nu2*nu2*nu2 - 448.*nu2*nu2*nu2*nu2) + 4.*nu1*nu1*(45. - 316.*nu2 + 676.*nu2*nu2 - 520.*nu2*nu2*nu2 + 112.*nu2*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)));

                    pnlpt->M12_oneline_mu6_vv_complex[count] = f * f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vv6_f3[count]);

                    //GC -> 7 of 7...

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_oneline_mu6_vv_complex_ortho[count] = f * f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vv6_f3_ortho[count]);
                        
                        //printf("%d\n",SWITCH_index); //GC - SWITCH!


                    //GC: ORTHOGONAL -- finish

                    }
                        
                    count++;
                }
            }

            for (index_i = 0; index_i < Nmax + 1; index_i++)
            {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_mu2_dd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 2. / (1. + 9. * nu1) * 9. * f * (1. + nu1);
                pnlpt->M13_mu2_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 2. / (1. + 9. * nu1) * (-1.) * (f * (7. + 9. * f - 9. * nu1));
                pnlpt->M13_mu4_vv_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] / (1. + 9. * nu1) * (-1.) * 3. * f * f * (5. + 6. * f - 3. * nu1);
                pnlpt->M13_mu4_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 2. / (1. + 9. * nu1) * 9. * f * f * (1. + 2. * nu1);
                pnlpt->M13_mu6_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 2. / (1. + 9. * nu1) * 9. * f * f * f * nu1;
            }

            double complex *f13_mu0_dd;
            class_alloc(f13_mu0_dd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu0_dd;
            class_alloc(f22_mu0_dd, Nmax * sizeof(complex double), pnlpt->error_message);

            double complex *f13_mu2_vd;
            class_alloc(f13_mu2_vd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f13_mu2_dd;
            class_alloc(f13_mu2_dd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f13_mu4_vv;
            class_alloc(f13_mu4_vv, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f13_mu4_vd;
            class_alloc(f13_mu4_vd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f13_mu6;
            class_alloc(f13_mu6, Nmax * sizeof(complex double), pnlpt->error_message);

            double complex *f22_mu2_vd;
            class_alloc(f22_mu2_vd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu2_dd;
            class_alloc(f22_mu2_dd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu4_vd;
            class_alloc(f22_mu4_vd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu4_dd;
            class_alloc(f22_mu4_dd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu4_vv;
            class_alloc(f22_mu4_vv, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu6_vv;
            class_alloc(f22_mu6_vv, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu6_vd;
            class_alloc(f22_mu6_vd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu8;
            class_alloc(f22_mu8, Nmax * sizeof(complex double), pnlpt->error_message);

            //GC!

            double complex *f12_mu0_dd;
            class_alloc(f12_mu0_dd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_mu2_vd;
            class_alloc(f12_mu2_vd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_mu2_dd;
            class_alloc(f12_mu2_dd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_mu4_vv;
            class_alloc(f12_mu4_vv, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_mu4_vd;
            class_alloc(f12_mu4_vd, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_mu6_vv;
            class_alloc(f12_mu6_vv, Nmax * sizeof(complex double), pnlpt->error_message);

            count = 0;
            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym_nw[count] * cpow(kdisc[index_j], etam[count]);
                }

                f13_mu0_dd[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_oneline_complex, &inc);
                P13UV_mu0_dd[index_j] = -61. * Pnw[index_j] * pow(kdisc[index_j], 2.) * sigmav / 105.;
                P13_mu0_dd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu0_dd[index_j] * Pnw[index_j]) + P13UV_mu0_dd[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));

                P13UV_mu2_dd[index_j] = -1. * Pnw[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav * f * (105. * f - 6.) / 105.;
                P13UV_mu2_vd[index_j] = -1. * Pnw[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav * f * (250. + 144. * f) / 105.;
                P13UV_mu4_vv[index_j] = -1. * Pnw[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav * f * f * (63. + 48. * f) / 35.;
                P13UV_mu4_vd[index_j] = -1. * Pnw[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav * f * f * (44. + 70. * f) / 35.;
                P13UV_mu6[index_j] = -1. * Pnw[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav * f * f * f * (46. + 35. * f) / 35.;

                f13_mu2_dd[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_mu2_dd_oneline_complex, &inc);
                f13_mu2_vd[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_mu2_vd_oneline_complex, &inc);
                f13_mu4_vv[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_mu4_vv_oneline_complex, &inc);
                f13_mu4_vd[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_mu4_vd_oneline_complex, &inc);
                f13_mu6[index_j] = zdotu_(&Nmaxf, x, &inc, pnlpt->M13_mu6_oneline_complex, &inc);

                P13_mu2_dd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu2_dd[index_j] * Pnw[index_j]) + P13UV_mu2_dd[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
                P13_mu2_vd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu2_vd[index_j] * Pnw[index_j]) + P13UV_mu2_vd[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
                P13_mu4_vv[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu4_vv[index_j] * Pnw[index_j]) + P13UV_mu4_vv[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
                P13_mu4_vd[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu4_vd[index_j] * Pnw[index_j]) + P13UV_mu4_vd[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
                P13_mu6[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu6[index_j] * Pnw[index_j]) + P13UV_mu6[index_j]) * exp(-pow(kdisc[index_j] / cutoff, 6.));
            }

            count = 0;
            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym_nw[count] * cpow(kdisc[index_j], etam[count]);
                    
                    if (SWITCH_index == 1) {
                    x_transfer[count] = cmsym_nw_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]); //GC -> need to be careful to use etam_transfer when computing the matrices!!!
                    }
                }

                //GC! Why doesn't he use a single loop also above? He uses zspmv_ many many times... Ah no, also here... The issue is that y get overwritten every time...

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_complex, x, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_mu0_dd[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_mu0_dd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_mu0_dd[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu0_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_mu0_dd[index_j] = Tnw[index_j] *                                                                                  //[FACTOR]*
                                      creal(cpow(kdisc[index_j], 3.) * f12_mu0_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC! This is the NON-WIGGLY part...
                }
                //        printf("%.16e %.16e %.16e",kdisc[index_j],Tnw[index_j],P12_mu0_dd[index_j]);
                //        printf("\n"); //GC!

                //GC: ORTHOGONAL -- start

                if (SWITCH_index == 1) {
                
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);                                //GC!
                f12_mu0_dd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                                  //GC!
                P12_mu0_dd_ortho[index_j] = Tnw[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_mu0_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC! This is the NON-WIGGLY part...
                }
                //printf("%.16e %.16e %.16e",kdisc[index_j],Tnw[index_j],P12_mu0_dd_ortho[index_j]);
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- finish

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu2_vd_complex, x, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu2_vd_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_mu2_vd[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_mu2_vd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_mu2_vd[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu2_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_mu2_vd[index_j] = Tnw[index_j] *                                                                                  //[FACTOR]*
                                      creal(cpow(kdisc[index_j], 3.) * f12_mu2_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //        printf("%.16e %.16e %.16e",kdisc[index_j],Tnw[index_j],P12_mu2_vd[index_j]);
                //        printf("\n"); //GC!

                //GC: ORTHOGONAL -- start
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu2_vd_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);                         //GC!
                f12_mu2_vd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                                  //GC!
                P12_mu2_vd_ortho[index_j] = Tnw[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_mu2_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //printf("%.16e %.16e %.16e",kdisc[index_j],Tnw[index_j],P12_mu2_vd_ortho[index_j]);
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- finish

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu2_dd_complex, x, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu2_dd_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_mu2_dd[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_mu2_dd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC -> as I said, I will need y_transfer here...
                }
                P22_mu2_dd[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu2_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_mu2_dd[index_j] = Tnw[index_j] *                                                                                  //[FACTOR]*
                                      creal(cpow(kdisc[index_j], 3.) * f12_mu2_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //        printf("%.16e %.16e %.16e",kdisc[index_j],Tnw[index_j],P12_mu2_dd[index_j]);
                //        printf("\n"); //GC!

                //GC: ORTHOGONAL -- start
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu2_dd_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);                         //GC!
                f12_mu2_dd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                                  //GC -> as I said, I will need y_transfer here...
                P12_mu2_dd_ortho[index_j] = Tnw[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_mu2_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //printf("%.16e %.16e %.16e",kdisc[index_j],Tnw[index_j],P12_mu2_dd_ortho[index_j]);
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- finish

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu4_vv_complex, x, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu4_vv_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_mu4_vv[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_mu4_vv[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_mu4_vv[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu4_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_mu4_vv[index_j] = Tnw[index_j] *                                                                                  //[FACTOR]*
                                      creal(cpow(kdisc[index_j], 3.) * f12_mu4_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //        printf("%.16e %.16e %.16e",kdisc[index_j],Tnw[index_j],P12_mu4_vv[index_j]);
                //        printf("\n"); //GC!

                //GC: ORTHOGONAL -- start
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu4_vv_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);                         //GC!
                f12_mu4_vv[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                                  //GC!
                P12_mu4_vv_ortho[index_j] = Tnw[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_mu4_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //printf("%.16e %.16e %.16e",kdisc[index_j],Tnw[index_j],P12_mu4_vv_ortho[index_j]);
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- finish

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu4_vd_complex, x, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu4_vd_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_mu4_vd[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_mu4_vd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_mu4_vd[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu4_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_mu4_vd[index_j] = Tnw[index_j] *                                                                                  //[FACTOR]*
                                      creal(cpow(kdisc[index_j], 3.) * f12_mu4_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //        printf("%.16e %.16e %.16e",kdisc[index_j],Tnw[index_j],P12_mu4_vd[index_j]);
                //        printf("\n"); //GC!

                //GC: ORTHOGONAL -- start
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu4_vd_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);                         //GC!
                f12_mu4_vd[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                                  //GC!
                P12_mu4_vd_ortho[index_j] = Tnw[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_mu4_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //printf("%.16e %.16e %.16e",kdisc[index_j],Tnw[index_j],P12_mu4_vd_ortho[index_j]);
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- finish

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu4_dd_complex, x, &inc, &beta, y, &inc);
                f22_mu4_dd[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                P22_mu4_dd[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu4_dd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu6_vv_complex, x, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu6_vv_complex, x_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_mu6_vv[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_mu6_vv[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_mu6_vv[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu6_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_mu6_vv[index_j] = Tnw[index_j] *                                                                                  //[FACTOR]*
                                      creal(cpow(kdisc[index_j], 3.) * f12_mu6_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //        printf("%.16e %.16e %.16e",kdisc[index_j],Tnw[index_j],P12_mu6_vv[index_j]);
                //        printf("\n"); //GC!

                //GC: ORTHOGONAL -- start
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu6_vv_complex_ortho, x_transfer, &inc, &beta, y_transfer, &inc);                         //GC!
                f12_mu6_vv[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                                  //GC!
                P12_mu6_vv_ortho[index_j] = Tnw[index_j] * creal(cpow(kdisc[index_j], 3.) * f12_mu6_vv[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //printf("%.16e %.16e %.16e",kdisc[index_j],Tnw[index_j],P12_mu6_vv_ortho[index_j]);
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- finish

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu6_vd_complex, x, &inc, &beta, y, &inc);
                f22_mu6_vd[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                P22_mu6_vd[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu6_vd[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu8_complex, x, &inc, &beta, y, &inc);
                f22_mu8[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                P22_mu8[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu8[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
            }

            // Matrix multiplication in case of Pw

            double *P13_mu0_dd_w;
            double *P22_mu0_dd_w;
            class_alloc(P13_mu0_dd_w, Nmax * sizeof(double), pnlpt->error_message);
            class_alloc(P22_mu0_dd_w, Nmax * sizeof(double), pnlpt->error_message);

            double *P13_mu2_dd_w;
            class_alloc(P13_mu2_dd_w, Nmax * sizeof(double), pnlpt->error_message);
            double *P13_mu2_vd_w;
            class_alloc(P13_mu2_vd_w, Nmax * sizeof(double), pnlpt->error_message);
            double *P13_mu4_vd_w;
            class_alloc(P13_mu4_vd_w, Nmax * sizeof(double), pnlpt->error_message);
            double *P13_mu4_vv_w;
            class_alloc(P13_mu4_vv_w, Nmax * sizeof(double), pnlpt->error_message);
            double *P13_mu6_w;
            class_alloc(P13_mu6_w, Nmax * sizeof(double), pnlpt->error_message);

            double *P22_mu2_vd_w;
            class_alloc(P22_mu2_vd_w, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu2_dd_w;
            class_alloc(P22_mu2_dd_w, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu4_vv_w;
            class_alloc(P22_mu4_vv_w, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu4_vd_w;
            class_alloc(P22_mu4_vd_w, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu4_dd_w;
            class_alloc(P22_mu4_dd_w, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu6_vv_w;
            class_alloc(P22_mu6_vv_w, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu6_vd_w;
            class_alloc(P22_mu6_vd_w, Nmax * sizeof(double), pnlpt->error_message);
            double *P22_mu8_w;
            class_alloc(P22_mu8_w, Nmax * sizeof(double), pnlpt->error_message);

            double complex *f13_mu0_dd_w;
            class_alloc(f13_mu0_dd_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu0_dd_w;
            class_alloc(f22_mu0_dd_w, Nmax * sizeof(complex double), pnlpt->error_message);

            double complex *f13_mu2_vd_w;
            class_alloc(f13_mu2_vd_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f13_mu2_dd_w;
            class_alloc(f13_mu2_dd_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f13_mu4_vv_w;
            class_alloc(f13_mu4_vv_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f13_mu4_vd_w;
            class_alloc(f13_mu4_vd_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f13_mu6_w;
            class_alloc(f13_mu6_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu2_vd_w;
            class_alloc(f22_mu2_vd_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu2_dd_w;
            class_alloc(f22_mu2_dd_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu4_vd_w;
            class_alloc(f22_mu4_vd_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu4_dd_w;
            class_alloc(f22_mu4_dd_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu4_vv_w;
            class_alloc(f22_mu4_vv_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu6_vv_w;
            class_alloc(f22_mu6_vv_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu6_vd_w;
            class_alloc(f22_mu6_vd_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f22_mu8_w;
            class_alloc(f22_mu8_w, Nmax * sizeof(complex double), pnlpt->error_message);

            //GC!

            double *P12_mu0_dd_w;
            class_alloc(P12_mu0_dd_w, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu2_vd_w;
            class_alloc(P12_mu2_vd_w, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu2_dd_w;
            class_alloc(P12_mu2_dd_w, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu4_vv_w;
            class_alloc(P12_mu4_vv_w, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu4_vd_w;
            class_alloc(P12_mu4_vd_w, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu6_vv_w;
            class_alloc(P12_mu6_vv_w, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- start

            double *P12_mu0_dd_w_ortho;
            class_alloc(P12_mu0_dd_w_ortho, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu2_vd_w_ortho;
            class_alloc(P12_mu2_vd_w_ortho, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu2_dd_w_ortho;
            class_alloc(P12_mu2_dd_w_ortho, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu4_vv_w_ortho;
            class_alloc(P12_mu4_vv_w_ortho, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu4_vd_w_ortho;
            class_alloc(P12_mu4_vd_w_ortho, Nmax * sizeof(double), pnlpt->error_message);

            double *P12_mu6_vv_w_ortho;
            class_alloc(P12_mu6_vv_w_ortho, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- finish

            double complex *f12_mu0_dd_w;
            class_alloc(f12_mu0_dd_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_mu2_vd_w;
            class_alloc(f12_mu2_vd_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_mu2_dd_w;
            class_alloc(f12_mu2_dd_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_mu4_vv_w;
            class_alloc(f12_mu4_vv_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_mu4_vd_w;
            class_alloc(f12_mu4_vd_w, Nmax * sizeof(complex double), pnlpt->error_message);
            double complex *f12_mu6_vv_w;
            class_alloc(f12_mu6_vv_w, Nmax * sizeof(complex double), pnlpt->error_message);

            //GC!

            count = 0;
            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    //   x[count]= cmsym_nw[count]* cpow(kdisc[index_j], etam[count]);
                    x_w[count] = cmsym_w[count] * cpow(kdisc[index_j], etam[count]);
                }

                f13_mu2_dd_w[index_j] = zdotu_(&Nmaxf, x_w, &inc, pnlpt->M13_mu2_dd_oneline_complex, &inc);
                f13_mu2_vd_w[index_j] = zdotu_(&Nmaxf, x_w, &inc, pnlpt->M13_mu2_vd_oneline_complex, &inc);
                f13_mu4_vv_w[index_j] = zdotu_(&Nmaxf, x_w, &inc, pnlpt->M13_mu4_vv_oneline_complex, &inc);
                f13_mu4_vd_w[index_j] = zdotu_(&Nmaxf, x_w, &inc, pnlpt->M13_mu4_vd_oneline_complex, &inc);
                f13_mu6_w[index_j] = zdotu_(&Nmaxf, x_w, &inc, pnlpt->M13_mu6_oneline_complex, &inc);
                f13_mu0_dd_w[index_j] = zdotu_(&Nmaxf, x_w, &inc, pnlpt->M13_oneline_complex, &inc);

                P13_mu0_dd_w[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu0_dd_w[index_j] * Pnw[index_j])) * exp(-pow(kdisc[index_j] / cutoff, 6.));
                P13_mu2_dd_w[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu2_dd_w[index_j] * Pnw[index_j])) * exp(-pow(kdisc[index_j] / cutoff, 6.));
                P13_mu2_vd_w[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu2_vd_w[index_j] * Pnw[index_j])) * exp(-pow(kdisc[index_j] / cutoff, 6.));
                P13_mu4_vv_w[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu4_vv_w[index_j] * Pnw[index_j])) * exp(-pow(kdisc[index_j] / cutoff, 6.));
                P13_mu4_vd_w[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu4_vd_w[index_j] * Pnw[index_j])) * exp(-pow(kdisc[index_j] / cutoff, 6.));
                P13_mu6_w[index_j] = (creal(cpow(kdisc[index_j], 3.) * f13_mu6_w[index_j] * Pnw[index_j])) * exp(-pow(kdisc[index_j] / cutoff, 6.));
            }

            count = 0;
            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x[count] = cmsym_nw[count] * cpow(kdisc[index_j], etam[count]);
                    x_w[count] = 2. * cmsym_w[count] * cpow(kdisc[index_j], etam[count]);
                    if (SWITCH_index == 1) {
                    x_transfer[count] = cmsym_nw_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]);
                    x_w_transfer[count] = 2. * cmsym_w_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]);
                    }
                    //GC -> recall that I still need to allocate and define this x_w_transfer (while y_transfer, that is not yet used here, doesn't need a y_w_transfer counterpart)... For now I just use it... Notice the factor of 2 => it is because of the functional expansion of a "P22-like" P_{1-loop}[Pnw+Pw]... I have to do the same. For the fact that I have to use etam_transfer also in other points above, see comments! There are many things that need to be fixed...
                }

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_complex, x_w, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_complex, x_w_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                //GC -> here the matrix does not have any names like M12_oneline_mu2_vd_complex because in the file matrix_summary-mu_powers-with_respect_to_matter.c you see that dd0 has ratio = 1! So it is JUST the matrix that I would use for 2F2, which is the PURE MATTER TERM IN REAL SPACE (here REAL means before redshift space)...
                f22_mu0_dd_w[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_mu0_dd_w[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC! NOT x_w_transfer here!!! It is like a P22. Notice that I will need y_transfer!!! It becomes like a P13 later after I multiply by the "ratio" 1+w/nw or whatever, see below... Will need to add some formulas to this code... The point here is that the wiggles OUTSIDE ONLY are damped by a \mu-dependent part... This is the reason for all this thing...
                }
                P22_mu0_dd_w[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu0_dd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_mu0_dd_w[index_j] = Tnw[index_j] *                                                                                    //[FACTOR]*
                                        creal(cpow(kdisc[index_j], 3.) * f12_mu0_dd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!

                //GC -> careful, I had put f12_mu0_dd...

                //GC: ORTHOGONAL -- start

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_complex_ortho, x_w_transfer, &inc, &beta, y_transfer, &inc); //GC!
                //GC -> here the matrix does not have any names like M12_oneline_mu2_vd_complex because in the file matrix_summary-mu_powers-with_respect_to_matter.c you see that dd0 has ratio = 1! So it is JUST the matrix that I would use for 2F2, which is the PURE MATTER TERM IN REAL SPACE (here REAL means before redshift space)...
                f12_mu0_dd_w[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                     //GC! NOT x_w_transfer here!!! It is like a P22. Notice that I will need y_transfer!!! It becomes like a P13 later after I multiply by the "ratio" 1+w/nw or whatever, see below... Will need to add some formulas to this code... The point here is that the wiggles OUTSIDE ONLY are damped by a \mu-dependent part... This is the reason for all this thing...
                P12_mu0_dd_w_ortho[index_j] = Tnw[index_j] *                                                                                    //[FACTOR]*
                                              creal(cpow(kdisc[index_j], 3.) * f12_mu0_dd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                    
                    //printf("%.16e\n",P12_mu0_dd_w_ortho[index_j]); //GC!
                    
                }
                //GC: ORTHOGONAL -- finish

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu2_vd_complex, x_w, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu2_vd_complex, x_w_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_mu2_vd_w[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_mu2_vd_w[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_mu2_vd_w[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu2_vd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_mu2_vd_w[index_j] = Tnw[index_j] *                                                                                    //[FACTOR]*
                                        creal(cpow(kdisc[index_j], 3.) * f12_mu2_vd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!

                //GC: ORTHOGONAL -- start

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu2_vd_complex_ortho, x_w_transfer, &inc, &beta, y_transfer, &inc);            //GC!
                f12_mu2_vd_w[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                     //GC!
                P12_mu2_vd_w_ortho[index_j] = Tnw[index_j] *                                                                                    //[FACTOR]*
                                              creal(cpow(kdisc[index_j], 3.) * f12_mu2_vd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //GC: ORTHOGONAL -- finish

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu2_dd_complex, x_w, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu2_dd_complex, x_w_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_mu2_dd_w[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_mu2_dd_w[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_mu2_dd_w[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu2_dd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_mu2_dd_w[index_j] = Tnw[index_j] *                                                                                    //[FACTOR]*
                                        creal(cpow(kdisc[index_j], 3.) * f12_mu2_dd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!

                //GC: ORTHOGONAL -- start

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu2_dd_complex_ortho, x_w_transfer, &inc, &beta, y_transfer, &inc);            //GC!
                f12_mu2_dd_w[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                     //GC!
                P12_mu2_dd_w_ortho[index_j] = Tnw[index_j] *                                                                                    //[FACTOR]*
                                              creal(cpow(kdisc[index_j], 3.) * f12_mu2_dd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //GC: ORTHOGONAL -- finish

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu4_vv_complex, x_w, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu4_vv_complex, x_w_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_mu4_vv_w[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_mu4_vv_w[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_mu4_vv_w[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu4_vv_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_mu4_vv_w[index_j] = Tnw[index_j] *                                                                                    //[FACTOR]*
                                        creal(cpow(kdisc[index_j], 3.) * f12_mu4_vv_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!

                //GC: ORTHOGONAL -- start

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu4_vv_complex_ortho, x_w_transfer, &inc, &beta, y_transfer, &inc);            //GC!
                f12_mu4_vv_w[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                     //GC!
                P12_mu4_vv_w_ortho[index_j] = Tnw[index_j] *                                                                                    //[FACTOR]*
                                              creal(cpow(kdisc[index_j], 3.) * f12_mu4_vv_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //GC: ORTHOGONAL -- finish

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu4_vd_complex, x_w, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu4_vd_complex, x_w_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_mu4_vd_w[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_mu4_vd_w[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_mu4_vd_w[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu4_vd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_mu4_vd_w[index_j] = Tnw[index_j] *                                                                                    //[FACTOR]*
                                        creal(cpow(kdisc[index_j], 3.) * f12_mu4_vd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!

                //GC: ORTHOGONAL -- start

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu4_vd_complex_ortho, x_w_transfer, &inc, &beta, y_transfer, &inc);            //GC!
                f12_mu4_vd_w[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                     //GC!
                P12_mu4_vd_w_ortho[index_j] = Tnw[index_j] *                                                                                    //[FACTOR]*
                                              creal(cpow(kdisc[index_j], 3.) * f12_mu4_vd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //GC: ORTHOGONAL -- finish

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu4_dd_complex, x_w, &inc, &beta, y, &inc);
                f22_mu4_dd_w[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                P22_mu4_dd_w[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu4_dd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu6_vv_complex, x_w, &inc, &beta, y, &inc);
                if (SWITCH_index == 1) {
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu6_vv_complex, x_w_transfer, &inc, &beta, y_transfer, &inc); //GC!
                }
                f22_mu6_vv_w[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                if (SWITCH_index == 1) {
                f12_mu6_vv_w[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc); //GC!
                }
                P22_mu6_vv_w[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu6_vv_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
                if (SWITCH_index == 1) {
                P12_mu6_vv_w[index_j] = Tnw[index_j] *                                                                                    //[FACTOR]*
                                        creal(cpow(kdisc[index_j], 3.) * f12_mu6_vv_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!

                //GC: ORTHOGONAL -- start

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_oneline_mu6_vv_complex_ortho, x_w_transfer, &inc, &beta, y_transfer, &inc);            //GC!
                f12_mu6_vv_w[index_j] = zdotu_(&Nmaxf, x_transfer, &inc, y_transfer, &inc);                                                     //GC!
                P12_mu6_vv_w_ortho[index_j] = Tnw[index_j] *                                                                                    //[FACTOR]*
                                              creal(cpow(kdisc[index_j], 3.) * f12_mu6_vv_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.))); //GC!
                }
                //GC: ORTHOGONAL -- finish

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu6_vd_complex, x_w, &inc, &beta, y, &inc);
                f22_mu6_vd_w[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                P22_mu6_vd_w[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu6_vd_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_oneline_mu8_complex, x_w, &inc, &beta, y, &inc);
                f22_mu8_w[index_j] = zdotu_(&Nmaxf, x, &inc, y, &inc);
                P22_mu8_w[index_j] = creal(cpow(kdisc[index_j], 3.) * f22_mu8_w[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.)));
            }

            //GC ~> why allocated in such strange order???

            double *dd_P13_mu4_vv;
            class_alloc(dd_P13_mu4_vv, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P13_mu4_vv, 1, dd_P13_mu4_vv, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P13_mu4_vv_ap_out = 0;

            double *dd_P22_mu4_vv;
            class_alloc(dd_P22_mu4_vv, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu4_vv, 1, dd_P22_mu4_vv, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu4_vv_ap_out = 0;

            double *dd_P22_mu4_vv_w;
            class_alloc(dd_P22_mu4_vv_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu4_vv_w, 1, dd_P22_mu4_vv_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu4_vv_w_ap_out = 0;

            double *dd_P13_mu4_vv_w;
            class_alloc(dd_P13_mu4_vv_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P13_mu4_vv_w, 1, dd_P13_mu4_vv_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P13_mu4_vv_w_ap_out = 0;

            double *dd_P13_mu6;
            class_alloc(dd_P13_mu6, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P13_mu6, 1, dd_P13_mu6, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P13_mu6_ap_out = 0;

            double *dd_P22_mu6_vv;
            class_alloc(dd_P22_mu6_vv, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu6_vv, 1, dd_P22_mu6_vv, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu6_vv_ap_out = 0;

            double *dd_P22_mu6_vv_w;
            class_alloc(dd_P22_mu6_vv_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu6_vv_w, 1, dd_P22_mu6_vv_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu6_vv_w_ap_out = 0;

            double *dd_P13_mu6_w;
            class_alloc(dd_P13_mu6_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P13_mu6_w, 1, dd_P13_mu6_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P13_mu6_w_ap_out = 0;

            double *dd_P22_mu8;
            class_alloc(dd_P22_mu8, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu8, 1, dd_P22_mu8, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu8_ap_out = 0;

            double *dd_P22_mu8_w;
            class_alloc(dd_P22_mu8_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu8_w, 1, dd_P22_mu8_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu8_w_ap_out = 0;

            double *dd_P13_mu0_dd;
            class_alloc(dd_P13_mu0_dd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P13_mu0_dd, 1, dd_P13_mu0_dd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P13_mu0_dd_ap_out = 0;

            double *dd_P22_mu0_dd;
            class_alloc(dd_P22_mu0_dd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu0_dd, 1, dd_P22_mu0_dd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu0_dd_ap_out = 0;

            double *dd_P13_mu0_dd_w;
            class_alloc(dd_P13_mu0_dd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P13_mu0_dd_w, 1, dd_P13_mu0_dd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P13_mu0_dd_w_ap_out = 0;

            double *dd_P22_mu0_dd_w;
            class_alloc(dd_P22_mu0_dd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu0_dd_w, 1, dd_P22_mu0_dd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu0_dd_w_ap_out = 0;

            double *dd_P22_mu2_dd;
            class_alloc(dd_P22_mu2_dd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu2_dd, 1, dd_P22_mu2_dd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu2_dd_ap_out = 0;

            double *dd_P13_mu2_dd;
            class_alloc(dd_P13_mu2_dd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P13_mu2_dd, 1, dd_P13_mu2_dd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P13_mu2_dd_ap_out = 0;

            double *dd_P22_mu2_dd_w;
            class_alloc(dd_P22_mu2_dd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu2_dd_w, 1, dd_P22_mu2_dd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu2_dd_w_ap_out = 0;

            double *dd_P13_mu2_dd_w;
            class_alloc(dd_P13_mu2_dd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P13_mu2_dd_w, 1, dd_P13_mu2_dd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P13_mu2_dd_w_ap_out = 0;

            double *dd_P22_mu4_dd;
            class_alloc(dd_P22_mu4_dd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu4_dd, 1, dd_P22_mu4_dd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu4_dd_ap_out = 0;

            double *dd_P22_mu4_dd_w;
            class_alloc(dd_P22_mu4_dd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu4_dd_w, 1, dd_P22_mu4_dd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu4_dd_w_ap_out = 0;

            double *dd_P13_mu2_vd;
            class_alloc(dd_P13_mu2_vd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P13_mu2_vd, 1, dd_P13_mu2_vd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P13_mu2_vd_ap_out = 0;

            double *dd_P22_mu2_vd;
            class_alloc(dd_P22_mu2_vd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu2_vd, 1, dd_P22_mu2_vd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu2_vd_ap_out = 0;

            double *dd_P22_mu2_vd_w;
            class_alloc(dd_P22_mu2_vd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu2_vd_w, 1, dd_P22_mu2_vd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu2_vd_w_ap_out = 0;

            double *dd_P13_mu2_vd_w;
            class_alloc(dd_P13_mu2_vd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P13_mu2_vd_w, 1, dd_P13_mu2_vd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P13_mu2_vd_w_ap_out = 0;

            double *dd_P13_mu4_vd;
            class_alloc(dd_P13_mu4_vd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P13_mu4_vd, 1, dd_P13_mu4_vd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P13_mu4_vd_ap_out = 0;

            double *dd_P22_mu4_vd;
            class_alloc(dd_P22_mu4_vd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu4_vd, 1, dd_P22_mu4_vd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu4_vd_ap_out = 0;

            double *dd_P22_mu4_vd_w;
            class_alloc(dd_P22_mu4_vd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu4_vd_w, 1, dd_P22_mu4_vd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu4_vd_w_ap_out = 0;

            double *dd_P13_mu4_vd_w;
            class_alloc(dd_P13_mu4_vd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P13_mu4_vd_w, 1, dd_P13_mu4_vd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P13_mu4_vd_w_ap_out = 0;

            double *dd_P22_mu6_vd;
            class_alloc(dd_P22_mu6_vd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu6_vd, 1, dd_P22_mu6_vd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu6_vd_ap_out = 0;

            double *dd_P22_mu6_vd_w;
            class_alloc(dd_P22_mu6_vd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P22_mu6_vd_w, 1, dd_P22_mu6_vd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P22_mu6_vd_w_ap_out = 0;

            //GC!

            double *dd_P12_mu6_vv;
            class_alloc(dd_P12_mu6_vv, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu6_vv, 1, dd_P12_mu6_vv, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu6_vv_ap_out = 0;

            double *dd_P12_mu4_vd;
            class_alloc(dd_P12_mu4_vd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu4_vd, 1, dd_P12_mu4_vd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu4_vd_ap_out = 0;

            double *dd_P12_mu4_vv;
            class_alloc(dd_P12_mu4_vv, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu4_vv, 1, dd_P12_mu4_vv, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu4_vv_ap_out = 0;

            double *dd_P12_mu2_vd;
            class_alloc(dd_P12_mu2_vd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu2_vd, 1, dd_P12_mu2_vd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu2_vd_ap_out = 0;

            double *dd_P12_mu2_dd;
            class_alloc(dd_P12_mu2_dd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu2_dd, 1, dd_P12_mu2_dd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu2_dd_ap_out = 0;

            double *dd_P12_mu0_dd;
            class_alloc(dd_P12_mu0_dd, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu0_dd, 1, dd_P12_mu0_dd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu0_dd_ap_out = 0;

            double *dd_P12_mu0_dd_w;
            class_alloc(dd_P12_mu0_dd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu0_dd_w, 1, dd_P12_mu0_dd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu0_dd_w_ap_out = 0;

            double *dd_P12_mu2_vd_w;
            class_alloc(dd_P12_mu2_vd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu2_vd_w, 1, dd_P12_mu2_vd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu2_vd_w_ap_out = 0;

            double *dd_P12_mu2_dd_w;
            class_alloc(dd_P12_mu2_dd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu2_dd_w, 1, dd_P12_mu2_dd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu2_dd_w_ap_out = 0;

            double *dd_P12_mu4_vv_w;
            class_alloc(dd_P12_mu4_vv_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu4_vv_w, 1, dd_P12_mu4_vv_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu4_vv_w_ap_out = 0;

            double *dd_P12_mu4_vd_w;
            class_alloc(dd_P12_mu4_vd_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu4_vd_w, 1, dd_P12_mu4_vd_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu4_vd_w_ap_out = 0;

            double *dd_P12_mu6_vv_w;
            class_alloc(dd_P12_mu6_vv_w, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu6_vv_w, 1, dd_P12_mu6_vv_w, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu6_vv_w_ap_out = 0;

            //GC!

            //GC: ORTHOGONAL -- start

            double *dd_P12_mu6_vv_ortho;
            class_alloc(dd_P12_mu6_vv_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu6_vv_ortho, 1, dd_P12_mu6_vv_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu6_vv_ap_out_ortho = 0;

            double *dd_P12_mu4_vd_ortho;
            class_alloc(dd_P12_mu4_vd_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu4_vd_ortho, 1, dd_P12_mu4_vd_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu4_vd_ap_out_ortho = 0;

            double *dd_P12_mu4_vv_ortho;
            class_alloc(dd_P12_mu4_vv_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu4_vv_ortho, 1, dd_P12_mu4_vv_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu4_vv_ap_out_ortho = 0;

            double *dd_P12_mu2_vd_ortho;
            class_alloc(dd_P12_mu2_vd_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu2_vd_ortho, 1, dd_P12_mu2_vd_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu2_vd_ap_out_ortho = 0;

            double *dd_P12_mu2_dd_ortho;
            class_alloc(dd_P12_mu2_dd_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu2_dd_ortho, 1, dd_P12_mu2_dd_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu2_dd_ap_out_ortho = 0;

            double *dd_P12_mu0_dd_ortho;
            class_alloc(dd_P12_mu0_dd_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu0_dd_ortho, 1, dd_P12_mu0_dd_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu0_dd_ap_out_ortho = 0;

            double *dd_P12_mu0_dd_w_ortho;
            class_alloc(dd_P12_mu0_dd_w_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu0_dd_w_ortho, 1, dd_P12_mu0_dd_w_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu0_dd_w_ap_out_ortho = 0;

            double *dd_P12_mu2_vd_w_ortho;
            class_alloc(dd_P12_mu2_vd_w_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu2_vd_w_ortho, 1, dd_P12_mu2_vd_w_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu2_vd_w_ap_out_ortho = 0;

            double *dd_P12_mu2_dd_w_ortho;
            class_alloc(dd_P12_mu2_dd_w_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu2_dd_w_ortho, 1, dd_P12_mu2_dd_w_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu2_dd_w_ap_out_ortho = 0;

            double *dd_P12_mu4_vv_w_ortho;
            class_alloc(dd_P12_mu4_vv_w_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu4_vv_w_ortho, 1, dd_P12_mu4_vv_w_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu4_vv_w_ap_out_ortho = 0;

            double *dd_P12_mu4_vd_w_ortho;
            class_alloc(dd_P12_mu4_vd_w_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu4_vd_w_ortho, 1, dd_P12_mu4_vd_w_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu4_vd_w_ap_out_ortho = 0;

            double *dd_P12_mu6_vv_w_ortho;
            class_alloc(dd_P12_mu6_vv_w_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_mu6_vv_w_ortho, 1, dd_P12_mu6_vv_w_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_mu6_vv_w_ap_out_ortho = 0;

            //GC: ORTHOGONAL -- finish

            //     printf("Until now all good\n");

            // Numerical integration over mu

            double Sigmatot = 0.;
            double LegendreP0 = 1.;
            double LegendreP2 = 0.;
            double LegendreP4 = 0.;
            double P1loopvv = 0.;
            double P1loopvd = 0.;
            double P1loopdd = 0.;
            double P12vv = 0.; //GC!
            double P12vd = 0.; //GC!
            double P12dd = 0.; //GC!
            //GC: ORTHOGONAL -- start
            double P12vv_ortho = 0.; //GC!
            double P12vd_ortho = 0.; //GC!
            double P12dd_ortho = 0.; //GC!
            //GC: ORTHOGONAL -- finish
            double Exp = 0.;
            double p_tree = 0.;
            double p_tree_vv = 0.;
            double p_tree_vd = 0.;
            double p_tree_dd = 0.;
            double P13ratio = 0.;
            double P12ratio = 0.; //GC!
            //GC: ORTHOGONAL ~> HERE I DO NOT NEED ITS OWN VERSION!!!
            double Pctr0 = 0.;
            double Pctr2 = 0.;
            double Pctr4 = 0.;
            double P1 = 0.;
            double P1b1 = 0.;
            double mu = 0.;
            double ktrue = 0.;
            double mutrue = 0.;
            last_index = 0;

            int index_gauss2 = 0;

            double P1loopvv_ap_ir = 0.;
            double P1loopvd_ap_ir = 0.;
            double P1loopdd_ap_ir = 0.;

            //GC!
            double P12vv_ap_ir = 0.;
            double P12vd_ap_ir = 0.;
            double P12dd_ap_ir = 0.;

            //GC: ORTHOGONAL -- start
            double P12vv_ap_ir_ortho = 0.;
            double P12vd_ap_ir_ortho = 0.;
            double P12dd_ap_ir_ortho = 0.;
            //GC: ORTHOGONAL -- finish

            /*
    //Dratio = 1.0; //GC!
    //hratio = 1.0; //GC!
    */

            //GC -> I cannot get them to 1 easily via the .ini???
            
            //GC - SWITCH ~> I keep these manipulations here. There is no segmentation fault, and they are not matrix manipulations. The brunt of the time losses was not at all here... Well, actually let me remove some calls...

            //    printf("Dratio=%lf\n",Dratio);
            //    printf("hratio=%lf\n",hratio);
            //printf("Dratio=%.16e\n",Dratio); //GC!
            //printf("hratio=%.16e\n",hratio); //GC!
            for (index_j = 0; index_j < Nside; index_j++)
            {

                P1loop_0_vv[index_j] = 0.;
                P1loop_0_dd[index_j] = 0.;
                P1loop_0_vd[index_j] = 0.;
                P1loop_2_vv[index_j] = 0.;
                P1loop_2_dd[index_j] = 0.;
                P1loop_2_vd[index_j] = 0.;
                P1loop_4_vv[index_j] = 0.;
                P1loop_4_dd[index_j] = 0.;
                P1loop_4_vd[index_j] = 0.;
                //GC!
                P12_0_vv[index_j] = 0.;
                P12_0_dd[index_j] = 0.;
                P12_0_vd[index_j] = 0.;
                P12_2_vv[index_j] = 0.;
                P12_2_dd[index_j] = 0.;
                P12_2_vd[index_j] = 0.;
                P12_4_vv[index_j] = 0.;
                P12_4_dd[index_j] = 0.;
                P12_4_vd[index_j] = 0.;
                //GC -> I didn't have a dd hexadecapole. Does the AP generate it??? Yes, but it is not here that... Ah no, it was the wrong order... Come on! One should be consistent... Then, I do not really get why I have P12vv, P12vd and P12dd. What's their function??? It is the sum of all the mu powers times actually mu to that power... So I need it... NOW -> notice that P12_4_dd, while not used above before AP, needs to be allocated! I allocate it above...
                //GC!
                //GC: ORTHOGONAL -- start
                P12_0_vv_ortho[index_j] = 0.;
                P12_0_dd_ortho[index_j] = 0.;
                P12_0_vd_ortho[index_j] = 0.;
                P12_2_vv_ortho[index_j] = 0.;
                P12_2_dd_ortho[index_j] = 0.;
                P12_2_vd_ortho[index_j] = 0.;
                P12_4_vv_ortho[index_j] = 0.;
                P12_4_dd_ortho[index_j] = 0.;
                P12_4_vd_ortho[index_j] = 0.;
                //GC: ORTHOGONAL -- finish
                P_CTR_0[index_j] = 0.;
                P_CTR_2[index_j] = 0.;
                P_CTR_4[index_j] = 0.;
                Ptree_0_vv[index_j] = 0.;
                Ptree_0_vd[index_j] = 0.;
                Ptree_0_dd[index_j] = 0.;
                Ptree_2_vv[index_j] = 0.;
                Ptree_2_vd[index_j] = 0.;
                Ptree_4_vv[index_j] = 0.;

                P1loop_0_vv[Nmax - 1 - index_j] = 0.;
                P1loop_0_dd[Nmax - 1 - index_j] = 0.;
                P1loop_0_vd[Nmax - 1 - index_j] = 0.;
                P1loop_2_vv[Nmax - 1 - index_j] = 0.;
                P1loop_2_dd[Nmax - 1 - index_j] = 0.;
                P1loop_2_vd[Nmax - 1 - index_j] = 0.;
                P1loop_4_vv[Nmax - 1 - index_j] = 0.;
                P1loop_4_dd[Nmax - 1 - index_j] = 0.;
                P1loop_4_vd[Nmax - 1 - index_j] = 0.;
                //GC!
                P12_0_vv[Nmax - 1 - index_j] = 0.;
                P12_0_dd[Nmax - 1 - index_j] = 0.;
                P12_0_vd[Nmax - 1 - index_j] = 0.;
                P12_2_vv[Nmax - 1 - index_j] = 0.;
                P12_2_dd[Nmax - 1 - index_j] = 0.;
                P12_2_vd[Nmax - 1 - index_j] = 0.;
                P12_4_vv[Nmax - 1 - index_j] = 0.;
                P12_4_dd[Nmax - 1 - index_j] = 0.;
                P12_4_vd[Nmax - 1 - index_j] = 0.;
                //GC!
                //GC: ORTHOGONAL -- start
                P12_0_vv_ortho[Nmax - 1 - index_j] = 0.;
                P12_0_dd_ortho[Nmax - 1 - index_j] = 0.;
                P12_0_vd_ortho[Nmax - 1 - index_j] = 0.;
                P12_2_vv_ortho[Nmax - 1 - index_j] = 0.;
                P12_2_dd_ortho[Nmax - 1 - index_j] = 0.;
                P12_2_vd_ortho[Nmax - 1 - index_j] = 0.;
                P12_4_vv_ortho[Nmax - 1 - index_j] = 0.;
                P12_4_dd_ortho[Nmax - 1 - index_j] = 0.;
                P12_4_vd_ortho[Nmax - 1 - index_j] = 0.;
                //GC: ORTHOGONAL -- finish
                P_CTR_0[Nmax - 1 - index_j] = 0.;
                P_CTR_2[Nmax - 1 - index_j] = 0.;
                P_CTR_4[Nmax - 1 - index_j] = 0.;
                Ptree_0_vv[Nmax - 1 - index_j] = 0.;
                Ptree_0_vd[Nmax - 1 - index_j] = 0.;
                Ptree_0_dd[Nmax - 1 - index_j] = 0.;
                Ptree_2_vv[Nmax - 1 - index_j] = 0.;
                Ptree_2_vd[Nmax - 1 - index_j] = 0.;
                Ptree_4_vv[Nmax - 1 - index_j] = 0.;
            }

            //    for (index_j=0; index_j < Nmax; index_j++){
            for (index_j = Nside; index_j < Nmax - Nside; index_j++)
            {

                P1loop_0_vv[index_j] = 0.;
                P1loop_0_dd[index_j] = 0.;
                P1loop_0_vd[index_j] = 0.;

                P1loop_2_vv[index_j] = 0.;
                P1loop_2_dd[index_j] = 0.;
                P1loop_2_vd[index_j] = 0.;

                P1loop_4_vv[index_j] = 0.;
                P1loop_4_dd[index_j] = 0.;
                P1loop_4_vd[index_j] = 0.;

                //GC!

                P12_0_vv[index_j] = 0.;
                P12_0_dd[index_j] = 0.;
                P12_0_vd[index_j] = 0.;

                P12_2_vv[index_j] = 0.;
                P12_2_dd[index_j] = 0.;
                P12_2_vd[index_j] = 0.;

                P12_4_vv[index_j] = 0.;
                P12_4_dd[index_j] = 0.;
                P12_4_vd[index_j] = 0.;

                //GC: ORTHOGONAL -- start

                P12_0_vv_ortho[index_j] = 0.;
                P12_0_dd_ortho[index_j] = 0.;
                P12_0_vd_ortho[index_j] = 0.;

                P12_2_vv_ortho[index_j] = 0.;
                P12_2_dd_ortho[index_j] = 0.;
                P12_2_vd_ortho[index_j] = 0.;

                P12_4_vv_ortho[index_j] = 0.;
                P12_4_dd_ortho[index_j] = 0.;
                P12_4_vd_ortho[index_j] = 0.;

                //GC: ORTHOGONAL -- finish

                //GC!

                P_CTR_0[index_j] = 0.;
                P_CTR_2[index_j] = 0.;
                P_CTR_4[index_j] = 0.;

                P10b1[index_j] = 0.;
                P10[index_j] = 0.;
                P12[index_j] = 0.;
                //GC -> I changed my added names to not confuse with this P12, which was there before...

                Ptree_0_vv[index_j] = 0.;
                Ptree_0_vd[index_j] = 0.;
                Ptree_0_dd[index_j] = 0.;
                Ptree_2_vv[index_j] = 0.;
                Ptree_2_vd[index_j] = 0.;
                Ptree_4_vv[index_j] = 0.;

                for (index_gauss2 = 0; index_gauss2 < 40; index_gauss2++)
                {

                    mu = pnlpt->gauss_x[index_gauss2];

                    if (pnlpt->AP_effect == AP_effect_yes)
                    {
                        mutrue = mu * hratio / pow((1. / Dratio / Dratio + (hratio * hratio - 1. / Dratio / Dratio) * mu * mu), 0.5);
                        ktrue = kdisc[index_j] * pow((1. / Dratio / Dratio + (hratio * hratio - 1. / Dratio / Dratio) * mu * mu), 0.5);
                    }

                    else
                    {
                        mutrue = mu;
                        ktrue = kdisc[index_j];
                    }

                    class_call(array_interpolate_spline(kdisc, Nmax, Pnw, dd_Pnw, 1, ktrue, &last_index, &Pnw_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, Pw, dd_Pw, 1, ktrue, &last_index, &Pw_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    if (SWITCH_index == 1) {
                    
                    //GC!

                    class_call(array_interpolate_spline(kdisc, Nmax, Tnw, dd_Tnw, 1, ktrue, &last_index, &Tnw_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, Tw, dd_Tw, 1, ktrue, &last_index, &Tw_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    //GC!
                        
                    }

                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu4_vv, dd_P22_mu4_vv, 1, ktrue, &last_index, &P22_mu4_vv_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P13_mu4_vv, dd_P13_mu4_vv, 1, ktrue, &last_index, &P13_mu4_vv_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu4_vv_w, dd_P22_mu4_vv_w, 1, ktrue, &last_index, &P22_mu4_vv_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P13_mu4_vv_w, dd_P13_mu4_vv_w, 1, ktrue, &last_index, &P13_mu4_vv_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P13_mu6, dd_P13_mu6, 1, ktrue, &last_index, &P13_mu6_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu6_vv, dd_P22_mu6_vv, 1, ktrue, &last_index, &P22_mu6_vv_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu6_vv_w, dd_P22_mu6_vv_w, 1, ktrue, &last_index, &P22_mu6_vv_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P13_mu6_w, dd_P13_mu6_w, 1, ktrue, &last_index, &P13_mu6_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu8, dd_P22_mu8, 1, ktrue, &last_index, &P22_mu8_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu8_w, dd_P22_mu8_w, 1, ktrue, &last_index, &P22_mu8_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu0_dd, dd_P22_mu0_dd, 1, ktrue, &last_index, &P22_mu0_dd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P13_mu0_dd, dd_P13_mu0_dd, 1, ktrue, &last_index, &P13_mu0_dd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P13_mu0_dd_w, dd_P13_mu0_dd_w, 1, ktrue, &last_index, &P13_mu0_dd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu0_dd_w, dd_P22_mu0_dd_w, 1, ktrue, &last_index, &P22_mu0_dd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu2_dd, dd_P22_mu2_dd, 1, ktrue, &last_index, &P22_mu2_dd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P13_mu2_dd, dd_P13_mu2_dd, 1, ktrue, &last_index, &P13_mu2_dd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu2_dd_w, dd_P22_mu2_dd_w, 1, ktrue, &last_index, &P22_mu2_dd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P13_mu2_dd_w, dd_P13_mu2_dd_w, 1, ktrue, &last_index, &P13_mu2_dd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu4_dd, dd_P22_mu4_dd, 1, ktrue, &last_index, &P22_mu4_dd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu4_dd_w, dd_P22_mu4_dd_w, 1, ktrue, &last_index, &P22_mu4_dd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P13_mu2_vd, dd_P13_mu2_vd, 1, ktrue, &last_index, &P13_mu2_vd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu2_vd, dd_P22_mu2_vd, 1, ktrue, &last_index, &P22_mu2_vd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu2_vd_w, dd_P22_mu2_vd_w, 1, ktrue, &last_index, &P22_mu2_vd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P13_mu2_vd_w, dd_P13_mu2_vd_w, 1, ktrue, &last_index, &P13_mu2_vd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P13_mu4_vd, dd_P13_mu4_vd, 1, ktrue, &last_index, &P13_mu4_vd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu4_vd, dd_P22_mu4_vd, 1, ktrue, &last_index, &P22_mu4_vd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu4_vd_w, dd_P22_mu4_vd_w, 1, ktrue, &last_index, &P22_mu4_vd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P13_mu4_vd_w, dd_P13_mu4_vd_w, 1, ktrue, &last_index, &P13_mu4_vd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu6_vd, dd_P22_mu6_vd, 1, ktrue, &last_index, &P22_mu6_vd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu6_vd_w, dd_P22_mu6_vd_w, 1, ktrue, &last_index, &P22_mu6_vd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P22_mu6_vd_w, dd_P22_mu6_vd_w, 1, ktrue, &last_index, &P22_mu6_vd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    if (SWITCH_index == 1) {
                    
                    //GC -> again in strange order...

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu0_dd, dd_P12_mu0_dd, 1, ktrue, &last_index, &P12_mu0_dd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu2_vd, dd_P12_mu2_vd, 1, ktrue, &last_index, &P12_mu2_vd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu2_dd, dd_P12_mu2_dd, 1, ktrue, &last_index, &P12_mu2_dd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu4_vv, dd_P12_mu4_vv, 1, ktrue, &last_index, &P12_mu4_vv_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu4_vd, dd_P12_mu4_vd, 1, ktrue, &last_index, &P12_mu4_vd_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu6_vv, dd_P12_mu6_vv, 1, ktrue, &last_index, &P12_mu6_vv_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    //GC: ORTHOGONAL -- start

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu0_dd_ortho, dd_P12_mu0_dd_ortho, 1, ktrue, &last_index, &P12_mu0_dd_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu2_vd_ortho, dd_P12_mu2_vd_ortho, 1, ktrue, &last_index, &P12_mu2_vd_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu2_dd_ortho, dd_P12_mu2_dd_ortho, 1, ktrue, &last_index, &P12_mu2_dd_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu4_vv_ortho, dd_P12_mu4_vv_ortho, 1, ktrue, &last_index, &P12_mu4_vv_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu4_vd_ortho, dd_P12_mu4_vd_ortho, 1, ktrue, &last_index, &P12_mu4_vd_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu6_vv_ortho, dd_P12_mu6_vv_ortho, 1, ktrue, &last_index, &P12_mu6_vv_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    //GC: ORTHOGONAL -- finish

                    //GC -> missing the wiggly ones... Added below! They are the same as the non-wiggly: dd_0, vd_2, dd_2, vv_4, vd_4, vv_6...

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu0_dd_w, dd_P12_mu0_dd_w, 1, ktrue, &last_index, &P12_mu0_dd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu2_vd_w, dd_P12_mu2_vd_w, 1, ktrue, &last_index, &P12_mu2_vd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu2_dd_w, dd_P12_mu2_dd_w, 1, ktrue, &last_index, &P12_mu2_dd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu4_vv_w, dd_P12_mu4_vv_w, 1, ktrue, &last_index, &P12_mu4_vv_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu4_vd_w, dd_P12_mu4_vd_w, 1, ktrue, &last_index, &P12_mu4_vd_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu6_vv_w, dd_P12_mu6_vv_w, 1, ktrue, &last_index, &P12_mu6_vv_w_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    //GC: ORTHOGONAL -- start

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu0_dd_w_ortho, dd_P12_mu0_dd_w_ortho, 1, ktrue, &last_index, &P12_mu0_dd_w_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu2_vd_w_ortho, dd_P12_mu2_vd_w_ortho, 1, ktrue, &last_index, &P12_mu2_vd_w_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu2_dd_w_ortho, dd_P12_mu2_dd_w_ortho, 1, ktrue, &last_index, &P12_mu2_dd_w_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu4_vv_w_ortho, dd_P12_mu4_vv_w_ortho, 1, ktrue, &last_index, &P12_mu4_vv_w_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu4_vd_w_ortho, dd_P12_mu4_vd_w_ortho, 1, ktrue, &last_index, &P12_mu4_vd_w_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_mu6_vv_w_ortho, dd_P12_mu6_vv_w_ortho, 1, ktrue, &last_index, &P12_mu6_vv_w_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    //GC: ORTHOGONAL -- finish

                    }
                        
                    //GC!

                    LegendreP2 = (3. * pow(mu, 2.) - 1.) / 2.;
                    LegendreP4 = (35. * pow(mu, 4.) - 30. * pow(mu, 2.) + 3.) / 8.;

                    Sigmatot = SigmaBAO * (1. + f * mutrue * mutrue * (2. + f)) + f * f * mutrue * mutrue * (mutrue * mutrue - 1.) * deltaSigmaBAO;
                    Exp = exp(-Sigmatot * pow(ktrue, 2.));

                    p_tree = (Pnw_ap_out + (1. + Sigmatot * pow(ktrue, 2.)) * Pw_ap_out * Exp);

                    P13ratio = 1. + (Pw_ap_out / Pnw_ap_out) * Exp;

                    if (SWITCH_index == 1) {
                    
                    P12ratio = 1. + (Tw_ap_out / Tnw_ap_out) * Exp; //GC!
                        
                    }

                    //GC: ORTHOGONAL -> this can remain unchanged...

                    P1b1 = (Pnw[index_j] + Pw[index_j] * Exp) * pnlpt->gauss_w[index_gauss2];
                    P1 = (Pnw[index_j] + Pw[index_j] * Exp) * f * pow(pnlpt->gauss_x[index_gauss2], 2.) * pnlpt->gauss_w[index_gauss2];

                    //GC -> I don't think I need to touch these two above...

                    p_tree_vv = (Pnw_ap_out + (1. + Sigmatot * pow(ktrue, 2.)) * Pw_ap_out * Exp) * pow(f * pow(mutrue, 2.), 2.) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    p_tree_vd = (Pnw_ap_out + (1. + Sigmatot * pow(ktrue, 2.)) * Pw_ap_out * Exp) * 2. * f * pow(mutrue, 2.) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    p_tree_dd = (Pnw_ap_out + (1. + Sigmatot * pow(ktrue, 2.)) * Pw_ap_out * Exp) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    Pctr0 = pow(ktrue, 2.) * (Pnw_ap_out + Pw_ap_out * Exp) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    Pctr2 = (Pnw_ap_out + Pw_ap_out * Exp) * pnlpt->gauss_w[index_gauss2] * f * pow(mutrue * ktrue, 2.) * hratio / Dratio / Dratio;
                    Pctr4 = pow(ktrue, 2.) * (Pnw_ap_out + Pw_ap_out * Exp) * pnlpt->gauss_w[index_gauss2] * f * f * pow(mutrue, 4.) * hratio / Dratio / Dratio;

                    //GC -> this is also not needed... Well, until I add the counterterms at least...

                    //GC -> dd_0, vd_2, dd_2, vv_4, vd_4, vv_6...

                    P1loopvv = (p_tree * 0. + (P13_mu4_vv_ap_out * P13ratio + P22_mu4_vv_ap_out + (P22_mu4_vv_w_ap_out + P13_mu4_vv_w_ap_out) * Exp) * pow(mutrue, 4.) + (P13_mu6_ap_out * P13ratio + P22_mu6_vv_ap_out + (P22_mu6_vv_w_ap_out + P13_mu6_w_ap_out) * Exp) * pow(mutrue, 6.) + (P22_mu8_ap_out + P22_mu8_w_ap_out * Exp) * pow(mutrue, 8.)) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio; //GC: Misha doesn't make a "vv", "vd", "dd" distinction for P13_mu6_ap_out, since it seems it doesn't appear beyond vv and \mu^6? I always make a distinction...

                    if (SWITCH_index == 1) {
                    
                    P12vv = ((P12_mu4_vv_ap_out * P12ratio + P12_mu4_vv_w_ap_out * Exp) * pow(mutrue, 4.) + (P12_mu6_vv_ap_out * P12ratio + P12_mu6_vv_w_ap_out * Exp) * pow(mutrue, 6.)) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio; //GC! ORIGINAL!!!

                    /*
            
            P12vv = ( (P12_mu4_vv_ap_out//*P12ratio
                       + 0.*P12_mu4_vv_w_ap_out*Exp)*pow(mutrue,4.) + (P12_mu6_vv_ap_out//*P12ratio
                                                                       + 0.*P12_mu6_vv_w_ap_out*Exp)*pow(mutrue,6.) )*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio; //GC -> ATTEMPT TO LOOK FOR ISSUES...
            
            */

                    //GC: ORTHOGONAL -- start

                    P12vv_ortho = ((P12_mu4_vv_ap_out_ortho * P12ratio + P12_mu4_vv_w_ap_out_ortho * Exp) * pow(mutrue, 4.) + (P12_mu6_vv_ap_out_ortho * P12ratio + P12_mu6_vv_w_ap_out_ortho * Exp) * pow(mutrue, 6.)) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio; //GC! ORIGINAL!!!

                    //GC: ORTHOGONAL -- finish
                        
                    }

                    P1loopdd = ((p_tree * 0. + P22_mu0_dd_ap_out + P13_mu0_dd_ap_out * P13ratio + (P13_mu0_dd_w_ap_out + P22_mu0_dd_w_ap_out) * Exp) + (P22_mu2_dd_ap_out + P13_mu2_dd_ap_out * P13ratio + (P22_mu2_dd_w_ap_out + P13_mu2_dd_w_ap_out) * Exp) * pow(mutrue, 2.) + (P22_mu4_dd_ap_out + P22_mu4_dd_w_ap_out * Exp) * pow(mutrue, 4.)) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                
                    if (SWITCH_index == 1) {
                        
                    P12dd = ((P12_mu0_dd_ap_out * P12ratio + P12_mu0_dd_w_ap_out * Exp) + (P12_mu2_dd_ap_out * P12ratio + P12_mu2_dd_w_ap_out * Exp) * pow(mutrue, 2.)) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio; //GC!

                    //GC: ORTHOGONAL -- start

                    P12dd_ortho = ((P12_mu0_dd_ap_out_ortho * P12ratio + P12_mu0_dd_w_ap_out_ortho * Exp) + (P12_mu2_dd_ap_out_ortho * P12ratio + P12_mu2_dd_w_ap_out_ortho * Exp) * pow(mutrue, 2.)) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio; //GC!

                    //GC: ORTHOGONAL -- finish
                        
                    }

                    //GC -> WHY IN THIS ORDER??? Come on...

                    P1loopvd = ((p_tree * 2. * f * 0. + P13_mu2_vd_ap_out * P13ratio + P22_mu2_vd_ap_out + (P22_mu2_vd_w_ap_out + P13_mu2_vd_w_ap_out) * Exp) * pow(mutrue, 2.) + (P13_mu4_vd_ap_out * P13ratio + P22_mu4_vd_ap_out + (P22_mu4_vd_w_ap_out + P13_mu4_vd_w_ap_out) * Exp) * pow(mutrue, 4.) + (P22_mu6_vd_ap_out + P22_mu6_vd_w_ap_out * Exp) * pow(mutrue, 6.)) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    if (SWITCH_index == 1) {
                    
                    P12vd = ((P12_mu2_vd_ap_out * P12ratio + P12_mu2_vd_w_ap_out * Exp) * pow(mutrue, 2.) + (P12_mu4_vd_ap_out * P12ratio + P12_mu4_vd_w_ap_out * Exp) * pow(mutrue, 4.)) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio; //GC!

                    //GC: ORTHOGONAL -- start

                    P12vd_ortho = ((P12_mu2_vd_ap_out_ortho * P12ratio + P12_mu2_vd_w_ap_out_ortho * Exp) * pow(mutrue, 2.) + (P12_mu4_vd_ap_out_ortho * P12ratio + P12_mu4_vd_w_ap_out_ortho * Exp) * pow(mutrue, 4.)) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio; //GC!

                    //GC: ORTHOGONAL -- finish
                        
                    }

                    P1loopdd_ap_ir = ((p_tree + P22_mu0_dd_ap_out + P13_mu0_dd_ap_out * P13ratio + (P13_mu0_dd_w_ap_out + P22_mu0_dd_w_ap_out) * Exp) + (P22_mu2_dd_ap_out + P13_mu2_dd_ap_out * P13ratio + (P22_mu2_dd_w_ap_out + P13_mu2_dd_w_ap_out) * Exp) * pow(mutrue, 2.) + (P22_mu4_dd_ap_out + P22_mu4_dd_w_ap_out * Exp) * pow(mutrue, 4.)) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    P1loopvd_ap_ir = ((p_tree * 2. * f + P13_mu2_vd_ap_out * P13ratio + P22_mu2_vd_ap_out + (P22_mu2_vd_w_ap_out + P13_mu2_vd_w_ap_out) * Exp) * pow(mutrue, 2.) + (P13_mu4_vd_ap_out * P13ratio + P22_mu4_vd_ap_out + (P22_mu4_vd_w_ap_out + P13_mu4_vd_w_ap_out) * Exp) * pow(mutrue, 4.) + (P22_mu6_vd_ap_out + P22_mu6_vd_w_ap_out * Exp) * pow(mutrue, 6.)) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    //GC -> what is this? It seems to be exactly as P1loopdd and P1loopvd, only with the tree-level contribution added...

                    P1loop_0_vv[index_j] += P1loopvv * LegendreP0 / 2.;
                    P1loop_2_vv[index_j] += P1loopvv * LegendreP2 * 2.5;
                    P1loop_4_vv[index_j] += P1loopvv * LegendreP4 * 4.5;

                    P1loop_0_dd[index_j] += P1loopdd * LegendreP0 / 2.;
                    P1loop_2_dd[index_j] += P1loopdd_ap_ir * LegendreP2 * 2.5;
                    P1loop_4_dd[index_j] += P1loopdd_ap_ir * LegendreP4 * 4.5;

                    P1loop_0_vd[index_j] += P1loopvd * LegendreP0 / 2.;
                    P1loop_2_vd[index_j] += P1loopvd * LegendreP2 * 2.5;
                    P1loop_4_vd[index_j] += P1loopvd_ap_ir * LegendreP4 * 4.5;

                    
                    if (SWITCH_index == 1) {
                        
                        
                    //GC!

                    P12_0_vv[index_j] += P12vv * LegendreP0 / 2.;
                    P12_2_vv[index_j] += P12vv * LegendreP2 * 2.5;
                    P12_4_vv[index_j] += P12vv * LegendreP4 * 4.5;

                    P12_0_dd[index_j] += P12dd * LegendreP0 / 2.;
                    P12_2_dd[index_j] += P12dd * LegendreP2 * 2.5; //GC: I do not have to worry about "ap_ir" statements...
                    P12_4_dd[index_j] += P12dd * LegendreP4 * 4.5;

                    P12_0_vd[index_j] += P12vd * LegendreP0 / 2.;
                    P12_2_vd[index_j] += P12vd * LegendreP2 * 2.5;
                    P12_4_vd[index_j] += P12vd * LegendreP4 * 4.5;

                    //GC: where are the integrals over \mu? Now, I have all the \mu powers multiplied by P_0, P_2 and P_4 (and 2\ell+1 divided by 2). So f[\mu] times P_\ell[\mu] times normalization... It seems that it is above, where the "for (index_gauss2=0; index_gauss2 < 40; index_gauss2++){..." code is, but where are the Gaussian weights? Yes, look that there is gauss_w multiplying everything...

                    //GC: ORTHOGONAL -- start

                    P12_0_vv_ortho[index_j] += P12vv_ortho * LegendreP0 / 2.;
                    P12_2_vv_ortho[index_j] += P12vv_ortho * LegendreP2 * 2.5;
                    P12_4_vv_ortho[index_j] += P12vv_ortho * LegendreP4 * 4.5;

                    P12_0_dd_ortho[index_j] += P12dd_ortho * LegendreP0 / 2.;
                    P12_2_dd_ortho[index_j] += P12dd_ortho * LegendreP2 * 2.5; //GC: I do not have to worry about "ap_ir" statements...
                    P12_4_dd_ortho[index_j] += P12dd_ortho * LegendreP4 * 4.5;

                    P12_0_vd_ortho[index_j] += P12vd_ortho * LegendreP0 / 2.;
                    P12_2_vd_ortho[index_j] += P12vd_ortho * LegendreP2 * 2.5;
                    P12_4_vd_ortho[index_j] += P12vd_ortho * LegendreP4 * 4.5;

                    //GC: ORTHOGONAL -- finish

                    //GC!
                        
                    }

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

                    P_CTR_0[index_j] += Pctr0 * LegendreP0 / 2.;
                    P_CTR_2[index_j] += Pctr2 * LegendreP2 * 2.5;
                    P_CTR_4[index_j] += Pctr4 * LegendreP4 * 4.5;

                    Ptree_0_vv[index_j] += p_tree_vv * LegendreP0 / 2.;
                    Ptree_0_vd[index_j] += p_tree_vd * LegendreP0 / 2.;
                    Ptree_0_dd[index_j] += p_tree_dd * LegendreP0 / 2.;
                    Ptree_2_vv[index_j] += p_tree_vv * LegendreP2 * 2.5;
                    Ptree_2_vd[index_j] += p_tree_vd * LegendreP2 * 2.5;
                    Ptree_4_vv[index_j] += p_tree_vv * LegendreP4 * 4.5;

                    P10b1[index_j] += P1b1 * LegendreP0 / 2.;
                    P10[index_j] += P1 * LegendreP0 / 2.;
                    P12[index_j] += P1 * LegendreP2 * 2.5;
                    //GC -> no change here. I changed my P12 -> P12_fNL to avoid confusion... The integral should be here since I have +=, and at each point I multiply by the corresponding Legendre polynomial. Notice that I still don't understand if he has the Gaussian weights or not...

                    /* P1loopvv = ((Pbin[index_j]*f*f + P13_mu4_vv[index_j]+P22_mu4_vv[index_j])*pow(pnlpt->gauss_x[index_gauss2],4.) + (P13_mu6[index_j]+P22_mu6_vv[index_j])*pow(pnlpt->gauss_x[index_gauss2],6.)+P22_mu8[index_j]*pow(pnlpt->gauss_x[index_gauss2],8.))*pnlpt->gauss_w[index_gauss2];
            P1loopdd = (Pbin[index_j] + P13[index_j] + P22[index_j] + (P13_mu2_dd[index_j] + P22_mu2_dd[index_j])*pow(pnlpt->gauss_x[index_gauss2],2.) + P22_mu4_dd[index_j]*pow(pnlpt->gauss_x[index_gauss2],4.))*pnlpt->gauss_w[index_gauss2];
            P1loopvd = (((Pbin[index_j])*2.*f + P13_mu2_vd[index_j] + P22_mu2_vd[index_j])*pow(pnlpt->gauss_x[index_gauss2],2.)+(P13_mu4_vd[index_j]+P22_mu4_vd[index_j])*pow(pnlpt->gauss_x[index_gauss2],4.)+P22_mu6_vd[index_j]*pow(pnlpt->gauss_x[index_gauss2],6.))*pnlpt->gauss_w[index_gauss2];*/
                    //  P1loopvv = (((Pnw[index_j] + (1. + Sigmatot*kdisc[index_j]*kdisc[index_j])* Pw[index_j] * exp(-Sigmatot * pow(kdisc[index_j],2.)))*f*f + P13_mu4_vv[index_j]+P22_mu4_vv[index_j])*pow(gauss_x[index_gauss],4.) + (P13_mu6[index_j]+P22_mu6_vv[index_j])*pow(gauss_x[index_gauss],6.)+P22_mu8[index_j]*pow(gauss_x[index_gauss],8.))*gauss_w[index_gauss];
                    //   P1loopdd = ((Pnw[index_j] + (1. + Sigmatot*kdisc[index_j]*kdisc[index_j])* Pw[index_j] * exp(-Sigmatot * pow(kdisc[index_j],2.))) + P13[index_j] + P22[index_j] + (P13_mu2_dd[index_j] + P22_mu2_dd[index_j])*pow(gauss_x[index_gauss],2.) + P22_mu4_dd[index_j]*pow(gauss_x[index_gauss],4.))*gauss_w[index_gauss];
                    //   P1loopvd = (((Pnw[index_j] + (1. + Sigmatot*kdisc[index_j]*kdisc[index_j])* Pw[index_j] * exp(-Sigmatot * pow(kdisc[index_j],2.)))*2.*f + P13_mu2_vd[index_j] + P22_mu2_vd[index_j])*pow(gauss_x[index_gauss],2.)+(P13_mu4_vd[index_j]+P22_mu4_vd[index_j])*pow(gauss_x[index_gauss],4.)+P22_mu6_vd[index_j]*pow(gauss_x[index_gauss],6.))*gauss_w[index_gauss];
                }

                //      printf("%lf %lf %lf %lf\n",kdisc[index_j],P1loop_4_vv[index_j],P1loop_4_vd[index_j],P1loop_4_dd[index_j]);

                //        printf("%lf %lf %lf %lf\n",kdisc[index_j],P_CTR_0[index_j],P_CTR_2[index_j],P_CTR_4[index_j]);

                //printf("%.16e %.16e\n",kdisc[index_j],P12_0_dd[index_j]); //GC -> CHECKED...

                P1loop_4_dd[index_j] = P1loop_4_dd[index_j];
                P1loop_4_vd[index_j] = P1loop_4_vd[index_j];
                P1loop_2_dd[index_j] = P1loop_2_dd[index_j];

                //GC??? I do not touch this... Why is it only for these specific ones? It is OUTSIDE of the loop that does the integral over \mu... Maybe there was some factor in front... Ok, let's leave this... Look that maybe these specific terms are the ones NOT generated by perturbation theory, only by AP...
            }

            //printf("{{--}}{{--}}{{--}}{{--}}{{--}}{{--}}{{--}}\n"); //GC!
            //printf("{{--}}{{--}}{{--}}{{--}}{{--}}{{--}}{{--}}\n"); //GC!

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

            //GC!

            free(dd_P12_mu6_vv);
            free(dd_P12_mu4_vd);
            free(dd_P12_mu4_vv);
            free(dd_P12_mu2_vd);
            free(dd_P12_mu2_dd);
            free(dd_P12_mu0_dd);

            free(dd_P12_mu6_vv_w);
            free(dd_P12_mu4_vd_w);
            free(dd_P12_mu4_vv_w);
            free(dd_P12_mu2_vd_w);
            free(dd_P12_mu2_dd_w);
            free(dd_P12_mu0_dd_w);

            free(f12_mu6_vv);
            free(f12_mu4_vd);
            free(f12_mu4_vv);
            free(f12_mu2_vd);
            free(f12_mu2_dd);
            free(f12_mu0_dd);

            free(f12_mu6_vv_w);
            free(f12_mu4_vd_w);
            free(f12_mu4_vv_w);
            free(f12_mu2_vd_w);
            free(f12_mu2_dd_w);
            free(f12_mu0_dd_w);

            free(P12_mu6_vv);
            free(P12_mu4_vd);
            free(P12_mu4_vv);
            free(P12_mu2_vd);
            free(P12_mu2_dd);
            free(P12_mu0_dd);

            free(P12_mu6_vv_w);
            free(P12_mu4_vd_w);
            free(P12_mu4_vv_w);
            free(P12_mu2_vd_w);
            free(P12_mu2_dd_w);
            free(P12_mu0_dd_w);

            free(cmsym_w_transfer);
            free(cmsym_nw_transfer);

            //GC: ORTHOGONAL -- start

            free(dd_P12_mu6_vv_ortho);
            free(dd_P12_mu4_vd_ortho);
            free(dd_P12_mu4_vv_ortho);
            free(dd_P12_mu2_vd_ortho);
            free(dd_P12_mu2_dd_ortho);
            free(dd_P12_mu0_dd_ortho);

            free(dd_P12_mu6_vv_w_ortho);
            free(dd_P12_mu4_vd_w_ortho);
            free(dd_P12_mu4_vv_w_ortho);
            free(dd_P12_mu2_vd_w_ortho);
            free(dd_P12_mu2_dd_w_ortho);
            free(dd_P12_mu0_dd_w_ortho);

            free(P12_mu6_vv_ortho);
            free(P12_mu4_vd_ortho);
            free(P12_mu4_vv_ortho);
            free(P12_mu2_vd_ortho);
            free(P12_mu2_dd_ortho);
            free(P12_mu0_dd_ortho);

            free(P12_mu6_vv_w_ortho);
            free(P12_mu4_vd_w_ortho);
            free(P12_mu4_vv_w_ortho);
            free(P12_mu2_vd_w_ortho);
            free(P12_mu2_dd_w_ortho);
            free(P12_mu0_dd_w_ortho);

            //GC: ORTHOGONAL -- finish

            //    printf("%lf\n",SigmaBAO);

        } // end of second IR resummation condition

        //printf("\n"); //GC!
        //printf("[*][*][*][*][*][*][*][*][*][*][*][*][*]\n"); //GC!
        //printf("[*][*][*][*][*][*][*][*][*][*][*][*][*]\n"); //GC!
        //printf("\n"); //GC!

        // Constructing the final output spectra
        
        //GC - SWITCH ~> THIS REMAINS, as it does for the matter-only part above... It just exports junk but no segmentation fault...

        double *ddpk_nl_0_vv;
        class_alloc(ddpk_nl_0_vv, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P1loop_0_vv,
                                              1,
                                              ddpk_nl_0_vv,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index = 0;
        pk_nl_out = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {
            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {
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
                pk_l_0_vv[index_k] = pk_nl_out + large_for_logs_big; //GC: notice that "l" here means LOOPS, not linear...
            }
            else
            {

                pk_l_0_vv[index_k] = -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * (f * f * (441. + 566. * f + 175. * f * f) / 1225.) + large_for_logs_big;
                //  pk_l_0_vv[index_k] = exp(lnpk_l[index_k])*(f*f/5.);

                //GC: minus makes sense -> having large vv term (and velocity dispersion) prevents from clustering... Careful: kmaxnew is the one of FFTlog, so 100. Same for the min. See around line 2160... He could have extrapolated by putting zero on either side. But now he puts the low-k limit also at high k, which is not optimal. I will need to find the low-k limit also for me... ir4dd is extrapolation. AP ignored but Misha knows it. It is under control... What did I mean here? The issue is that for me the low-k limit will not be proportional to P? It happens for the contribution that is ir4dd: it is just the fact that in the low-k limit I have k^4 times \int P^2... I do not really know how to fix it yet... Whatever I do, I will neglect the AP effect. For now I will put the else to 1.e7...

                //   -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*(f*f*(441.+566.*f+175.*f*f)/1225.)
            }
            //     printf("%le %le\n",pnlpt->k[index_k]/pba->h,pk_l_0_vv[index_k]*pow(pba->h,3));

            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_0_vv[index_k] - 1.e7); //GC -> check if also Misha's results change because of AP...

            //printf("%.24e %.24e\n",pnlpt->k[index_k],pk_l_0_vv[index_k] - 1.e7); //GC -> check if also Misha's results change because of AP...
        }

        double *ddpk_nl_0_vd;
        class_alloc(ddpk_nl_0_vd, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P1loop_0_vd,
                                              1,
                                              ddpk_nl_0_vd,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_out = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

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

                pk_l_0_vd[index_k] = pk_nl_out + large_for_logs_big;
            }
            else
            {
                //      pk_l_0_vd[index_k] = (f*2./3.)*exp(lnpk_l[index_k]);
                pk_l_0_vd[index_k] = -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * (2. * f * (625. + 558. * f + 315. * f * f) / 1575.) + large_for_logs_big;
            }
            //         printf("%le %le\n",pnlpt->k[index_k]/pba->h,pk_l_0_vd[index_k]*pow(pba->h,3));
        }

        double *ddpk_nl_0_dd;
        class_alloc(ddpk_nl_0_dd, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P1loop_0_dd,
                                              1,
                                              ddpk_nl_0_dd,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_out = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

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

                pk_l_0_dd[index_k] = pk_nl_out + large_for_logs_big;
            }
            else
            {
                //  pk_l_0_dd[index_k] = exp(lnpk_l[index_k]);
                pk_l_0_dd[index_k] = -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * ((61. - 2. * f + 35. * f * f) / 105.) + large_for_logs_big;
            }
            //  -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*((61. -2.*f + 35.*f*f)/105.);
            //     printf("%le %le\n",pnlpt->k[index_k]/pba->h,pk_l_0_dd[index_k]*pow(pba->h,3));
        }

        double *ddpk_nl_2_vv;
        class_alloc(ddpk_nl_2_vv, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P1loop_2_vv,
                                              1,
                                              ddpk_nl_2_vv,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_out = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

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

                pk_l_2_vv[index_k] = pk_nl_out + large_for_logs_big;
            }

            else
            {
                pk_l_2_vv[index_k] = -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * (2. * f * f * (54. + 74. * f + 25. * f * f) / 105.) + large_for_logs_big;
                //pk_l_2_vv[index_k] = exp(lnpk_l[index_k])*(f*f*4./7.);
                // -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*(2.*f*f*(54.+74.*f+25.*f*f)/105.);
            }
        }

        double *ddpk_nl_2_vd;
        class_alloc(ddpk_nl_2_vd, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P1loop_2_vd,
                                              1,
                                              ddpk_nl_2_vd,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_out = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

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

                pk_l_2_vd[index_k] = pk_nl_out + large_for_logs_big;
            }

            else
            {
                pk_l_2_vd[index_k] = -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * 4. * f * (175. + 180. * f + 126. * f * f) / 441. + large_for_logs_big;
            }
            //    -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*4.*f*(175.+180.*f+126.*f*f)/441.;
        }

        double *ddpk_nl_2_dd;
        class_alloc(ddpk_nl_2_dd, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P1loop_2_dd,
                                              1,
                                              ddpk_nl_2_dd,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_out = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

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

                pk_l_2_dd[index_k] = pk_nl_out + large_for_logs_big;
            }

            else
            {
                pk_l_2_dd[index_k] = -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * (2. * f * (35. * f - 2.) / 105.) + large_for_logs_big;
            }

            //    -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*(2.*f*(35.*f-2.)/105.);
        }

        double *ddpk_nl_4_vv;
        class_alloc(ddpk_nl_4_vv, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P1loop_4_vv,
                                              1,
                                              ddpk_nl_4_vv,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_out = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

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

                pk_l_4_vv[index_k] = pk_nl_out + large_for_logs_big;
            }

            else
            {
                pk_l_4_vv[index_k] = -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * (24. * f * f * (33. + 58. * f + 25. * f * f) / 1925.) + large_for_logs_big;
                //-1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*(24.*f*f*(33.+58.*f+25.*f*f)/1925.);
            }
        }

        double *ddpk_nl_4_vd;
        class_alloc(ddpk_nl_4_vd, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P1loop_4_vd,
                                              1,
                                              ddpk_nl_4_vd,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_out = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

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

                pk_l_4_vd[index_k] = pk_nl_out + large_for_logs_big;
            }

            else
            {
                pk_l_4_vd[index_k] = -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * 16. * f * f * (22. + 35. * f) / 1225. + large_for_logs_big;
            }
            //    -1.*Pbin[index_j] * kdisc[index_j]*kdisc[index_j]* sigmav*16.*f*f*(22.+35.*f)/1225.;
        }

        double *ddpk_nl_4_dd;
        class_alloc(ddpk_nl_4_dd, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P1loop_4_dd,
                                              1,
                                              ddpk_nl_4_dd,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_out = 0;
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

        last_index = 0;
        pk_nl_out = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

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

                pk_l_4_dd[index_k] = pk_nl_out * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_big;
            }

            else
            {
                pk_l_4_dd[index_k] = large_for_logs_big + ir4dd * pow((pnlpt->k[index_k] / kminnew), 4.) * exp(-pow(pnlpt->k[index_k] / 3., 4.));
            }
        }

        //GC: monopole dd...

        double *ddpk_nl_fNL_0_dd;
        class_alloc(ddpk_nl_fNL_0_dd, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_0_dd,
                                              1,
                                              ddpk_nl_fNL_0_dd,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out = 0;

        /*
    double irfNL0dd; //GC -> COULD BE USELESS! ESPECIALLY SINCE I MOST LIKELY DO NOT HAVE A SIMPLE POWER LAW BUT I GET THE SQUARE ROOT OF THE POWER SPECTRUM! I WILL NEED TO CHECK THIS!!!
    class_call(array_interpolate_spline(kdisc,
                                        Nmax,
                                        P12_0_dd,
                                        ddpk_nl_fNL_0_dd,
                                        1,
                                        kminnew,
                                        &last_index,
                                        &irfNL0dd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    last_index=0;
    pk_nl_fNL_out=0; //GC: redundant???
    */

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_0_dd,
                                                    ddpk_nl_fNL_0_dd,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_0_dd[index_k] = pk_nl_fNL_out * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_0_dd[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL0dd*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_0_dd[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_0_dd[index_k]-0.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- start

        double *ddpk_nl_fNL_0_dd_ortho;
        class_alloc(ddpk_nl_fNL_0_dd_ortho, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_0_dd_ortho,
                                              1,
                                              ddpk_nl_fNL_0_dd_ortho,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out_ortho = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_0_dd_ortho,
                                                    ddpk_nl_fNL_0_dd_ortho,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out_ortho,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_0_dd_ortho[index_k] = pk_nl_fNL_out_ortho * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_0_dd_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL0dd*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_0_dd[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_0_dd[index_k]-0.*large_for_logs_fNL); //GC!

            //printf("%.20f %.20f\n",pnlpt->k[index_k],pk_l_fNL_0_dd_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            //printf("%.20e %.20e\n",pnlpt->k[index_k],pk_l_fNL_0_dd_ortho[index_k]-1.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- finish

        //GC: monopole vd...

        double *ddpk_nl_fNL_0_vd;
        class_alloc(ddpk_nl_fNL_0_vd, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_0_vd,
                                              1,
                                              ddpk_nl_fNL_0_vd,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out = 0;

        /*

     double irfNL0vd; //GC -> COULD BE USELESS! ESPECIALLY SINCE I MOST LIKELY DO NOT HAVE A SIMPLE POWER LAW BUT I GET THE SQUARE ROOT OF THE POWER SPECTRUM! I WILL NEED TO CHECK THIS!!!
    class_call(array_interpolate_spline(kdisc,
                                        Nmax,
                                        P12_0_vd,
                                        ddpk_nl_fNL_0_vd,
                                        1,
                                        kminnew,
                                        &last_index,
                                        &irfNL0vd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    last_index=0;
    pk_nl_fNL_out=0; //GC: redundant??? Notice that BELOW I essentially repeated the same command twice... It is not the only place where this happens in this code, so whatever...
     
     */

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_0_vd,
                                                    ddpk_nl_fNL_0_vd,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_0_vd[index_k] = pk_nl_fNL_out * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;

                //GC -> notice that there is ANOTHER cutoff here, on top of the one that was already present when doing the matrix multiplication... I really don't get it... In any case -> could this be the reason for the problems of the double zero at low k? I cannot believe it...
            }

            else
            {
                pk_l_fNL_0_vd[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL0vd*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.20f %.20f\n",pnlpt->k[index_k],pk_l_fNL_0_vd[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_0_vd[index_k]-0.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- start

        double *ddpk_nl_fNL_0_vd_ortho;
        class_alloc(ddpk_nl_fNL_0_vd_ortho, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_0_vd_ortho,
                                              1,
                                              ddpk_nl_fNL_0_vd_ortho,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out_ortho = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_0_vd_ortho,
                                                    ddpk_nl_fNL_0_vd_ortho,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out_ortho,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_0_vd_ortho[index_k] = pk_nl_fNL_out_ortho * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
            }

            else
            {
                pk_l_fNL_0_vd_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL0vd*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.20f %.20f\n",pnlpt->k[index_k],pk_l_fNL_0_vd[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_0_vd[index_k]-0.*large_for_logs_fNL); //GC!

            //printf("%.20f %.20f\n",pnlpt->k[index_k],pk_l_fNL_0_vd_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            //printf("%.20e %.20e\n",pnlpt->k[index_k],pk_l_fNL_0_vd_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_0_vd[index_k]-1.*large_for_logs_fNL); //GC!

            //GC: could have been that the oscillation at extremely low k when printing the interpolation was due to using f instead of e and losing some precision? Whatever... I already investigated this enough, and it was a non-issue...
        }

        //GC: ORTHOGONAL -- finish

        //GC: monopole vv...

        double *ddpk_nl_fNL_0_vv;
        class_alloc(ddpk_nl_fNL_0_vv, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_0_vv,
                                              1,
                                              ddpk_nl_fNL_0_vv,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out = 0;

        /*
    
    double irfNL0vv; //GC -> COULD BE USELESS! ESPECIALLY SINCE I MOST LIKELY DO NOT HAVE A SIMPLE POWER LAW BUT I GET THE SQUARE ROOT OF THE POWER SPECTRUM! I WILL NEED TO CHECK THIS!!!
    class_call(array_interpolate_spline(kdisc,
                                        Nmax,
                                        P12_0_vv,
                                        ddpk_nl_fNL_0_vv,
                                        1,
                                        kminnew,
                                        &last_index,
                                        &irfNL0vv,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);
     
     //GC -> as I said, for now we do not use this thing about the low-k and high-k limits...
     
     */

        /*
    
    last_index=0;
    pk_nl_fNL_out=0; //GC: redundant???
    
     //GC: it is redundant -> literally the same command is repeated twice...
    
     */

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_0_vv,
                                                    ddpk_nl_fNL_0_vv,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_0_vv[index_k] = pk_nl_fNL_out * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; // + 1.*epsilon_for_logs_fNL; //+1.e7; //GC!!!
            }

            else
            {
                pk_l_fNL_0_vv[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL0vv*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.));
                //GC -> as discussed, we just put this to ZERO or whatever (see above) outside...
            }

            //printf("%.16f %.16f \n",pnlpt->k[index_k],pk_l_fNL_0_vv[index_k]-0.*large_for_logs_fNL); //GC!!!
            //printf("%.20f %.20f\n",pnlpt->k[index_k],pk_l_fNL_0_vv[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_0_vv[index_k]-0.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- start

        double *ddpk_nl_fNL_0_vv_ortho;
        class_alloc(ddpk_nl_fNL_0_vv_ortho, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_0_vv_ortho,
                                              1,
                                              ddpk_nl_fNL_0_vv_ortho,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out_ortho = 0;

        //last_index=0;
        //pk_nl_fNL_out=0; //GC: redundant???

        //GC: ORTHOGONAL -> yes, it is redundant...

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_0_vv_ortho,
                                                    ddpk_nl_fNL_0_vv_ortho,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out_ortho,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_0_vv_ortho[index_k] = pk_nl_fNL_out_ortho * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; // + 1.*epsilon_for_logs_fNL; //+1.e7; //GC!!!
            }

            else
            {
                pk_l_fNL_0_vv_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL0vv*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.));
                //GC -> as discussed, we just put this to ZERO or whatever (see above) outside...
            }

            //printf("%.16f %.16f \n",pnlpt->k[index_k],pk_l_fNL_0_vv[index_k]-0.*large_for_logs_fNL); //GC!!!
            //printf("%.20f %.20f\n",pnlpt->k[index_k],pk_l_fNL_0_vv[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_0_vv[index_k]-0.*large_for_logs_fNL); //GC!

            //printf("%.20f %.20f\n",pnlpt->k[index_k],pk_l_fNL_0_vv_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            //printf("%.20e %.20e\n",pnlpt->k[index_k],pk_l_fNL_0_vv_ortho[index_k]-1.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- finish

        //GC: quadrupole dd...

        double *ddpk_nl_fNL_2_dd;
        class_alloc(ddpk_nl_fNL_2_dd, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_2_dd,
                                              1,
                                              ddpk_nl_fNL_2_dd,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out = 0;

        /*
    
    double irfNL2dd; //GC -> COULD BE USELESS! ESPECIALLY SINCE I MOST LIKELY DO NOT HAVE A SIMPLE POWER LAW BUT I GET THE SQUARE ROOT OF THE POWER SPECTRUM! I WILL NEED TO CHECK THIS!!!
    class_call(array_interpolate_spline(kdisc,
                                        Nmax,
                                        P12_2_dd,
                                        ddpk_nl_fNL_2_dd,
                                        1,
                                        kminnew,
                                        &last_index,
                                        &irfNL2dd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    last_index=0;
    pk_nl_fNL_out=0; //GC: redundant???
    
    */

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_2_dd,
                                                    ddpk_nl_fNL_2_dd,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_2_dd[index_k] = pk_nl_fNL_out * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_2_dd[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!
                //1.e7 + irfNL2dd*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_l_fNL_2_dd[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_2_dd[index_k]-0.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- start

        double *ddpk_nl_fNL_2_dd_ortho;
        class_alloc(ddpk_nl_fNL_2_dd_ortho, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_2_dd_ortho,
                                              1,
                                              ddpk_nl_fNL_2_dd_ortho,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out_ortho = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_2_dd_ortho,
                                                    ddpk_nl_fNL_2_dd_ortho,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out_ortho,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_2_dd_ortho[index_k] = pk_nl_fNL_out_ortho * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_2_dd_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!
                //1.e7 + irfNL2dd*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_l_fNL_2_dd[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_2_dd[index_k]-0.*large_for_logs_fNL); //GC!

            //printf("%.20f %.20f\n",pnlpt->k[index_k],pk_l_fNL_2_dd_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            //printf("%.20e %.20e\n",pnlpt->k[index_k],pk_l_fNL_2_dd_ortho[index_k]-1.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- finish

        //GC: quadrupole vd...

        double *ddpk_nl_fNL_2_vd;
        class_alloc(ddpk_nl_fNL_2_vd, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_2_vd,
                                              1,
                                              ddpk_nl_fNL_2_vd,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out = 0;

        /*
    
    double irfNL2vd; //GC -> COULD BE USELESS! ESPECIALLY SINCE I MOST LIKELY DO NOT HAVE A SIMPLE POWER LAW BUT I GET THE SQUARE ROOT OF THE POWER SPECTRUM! I WILL NEED TO CHECK THIS!!!
    class_call(array_interpolate_spline(kdisc,
                                        Nmax,
                                        P12_2_vd,
                                        ddpk_nl_fNL_2_vd,
                                        1,
                                        kminnew,
                                        &last_index,
                                        &irfNL2vd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    last_index=0;
    pk_nl_fNL_out=0; //GC: redundant???
    
    */

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_2_vd,
                                                    ddpk_nl_fNL_2_vd,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_2_vd[index_k] = pk_nl_fNL_out * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_2_vd[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL2vd*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_l_fNL_2_vd[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_2_vd[index_k]-0.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- start

        double *ddpk_nl_fNL_2_vd_ortho;
        class_alloc(ddpk_nl_fNL_2_vd_ortho, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_2_vd_ortho,
                                              1,
                                              ddpk_nl_fNL_2_vd_ortho,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out_ortho = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_2_vd_ortho,
                                                    ddpk_nl_fNL_2_vd_ortho,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out_ortho,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_2_vd_ortho[index_k] = pk_nl_fNL_out_ortho * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_2_vd_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL2vd*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_l_fNL_2_vd[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_2_vd[index_k]-0.*large_for_logs_fNL); //GC!

            //printf("%.20f %.20f\n",pnlpt->k[index_k],pk_l_fNL_2_vd_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            //printf("%.20e %.20e\n",pnlpt->k[index_k],pk_l_fNL_2_vd_ortho[index_k]-1.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- finish

        //GC: quadrupole vv...

        double *ddpk_nl_fNL_2_vv;
        class_alloc(ddpk_nl_fNL_2_vv, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_2_vv,
                                              1,
                                              ddpk_nl_fNL_2_vv,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out = 0;

        /*
    
    double irfNL2vv; //GC -> COULD BE USELESS! ESPECIALLY SINCE I MOST LIKELY DO NOT HAVE A SIMPLE POWER LAW BUT I GET THE SQUARE ROOT OF THE POWER SPECTRUM! I WILL NEED TO CHECK THIS!!!
    class_call(array_interpolate_spline(kdisc,
                                        Nmax,
                                        P12_2_vv,
                                        ddpk_nl_fNL_2_vv,
                                        1,
                                        kminnew,
                                        &last_index,
                                        &irfNL2vv,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    last_index=0;
    pk_nl_fNL_out=0; //GC: redundant???
    
    */

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_2_vv,
                                                    ddpk_nl_fNL_2_vv,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_2_vv[index_k] = pk_nl_fNL_out * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_2_vv[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL2vv*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_l_fNL_2_vv[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_2_vv[index_k]-0.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- start

        double *ddpk_nl_fNL_2_vv_ortho;
        class_alloc(ddpk_nl_fNL_2_vv_ortho, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_2_vv_ortho,
                                              1,
                                              ddpk_nl_fNL_2_vv_ortho,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out_ortho = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_2_vv_ortho,
                                                    ddpk_nl_fNL_2_vv_ortho,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out_ortho,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_2_vv_ortho[index_k] = pk_nl_fNL_out_ortho * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_2_vv_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL2vv*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_l_fNL_2_vv[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_2_vv[index_k]-0.*large_for_logs_fNL); //GC!

            //printf("%.20f %.20f\n",pnlpt->k[index_k],pk_l_fNL_2_vv_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            //printf("%.20e %.20e\n",pnlpt->k[index_k],pk_l_fNL_2_vv_ortho[index_k]-1.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- finish

        //GC: hexadecapole dd...

        double *ddpk_nl_fNL_4_dd;
        class_alloc(ddpk_nl_fNL_4_dd, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_4_dd,
                                              1,
                                              ddpk_nl_fNL_4_dd,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out = 0;

        /*
    
    double irfNL4dd; //GC -> COULD BE USELESS! ESPECIALLY SINCE I MOST LIKELY DO NOT HAVE A SIMPLE POWER LAW BUT I GET THE SQUARE ROOT OF THE POWER SPECTRUM! I WILL NEED TO CHECK THIS!!! No, I checked and it already exports only zeros... Great -> notice that also there is the fact that I am now adding a 1.0*epsilon_for_logs_fNL likely affects the checks I was doing on C vs. C...
    class_call(array_interpolate_spline(kdisc,
                                        Nmax,
                                        P12_4_dd,
                                        ddpk_nl_fNL_4_dd,
                                        1,
                                        kminnew,
                                        &last_index,
                                        &irfNL4dd,
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    last_index=0;
    pk_nl_fNL_out=0; //GC: redundant???
    
    */

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_4_dd,
                                                    ddpk_nl_fNL_4_dd,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_4_dd[index_k] = pk_nl_fNL_out * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_4_dd[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL4dd*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_l_fNL_4_dd[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_4_dd[index_k]-0.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- start

        double *ddpk_nl_fNL_4_dd_ortho;
        class_alloc(ddpk_nl_fNL_4_dd_ortho, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_4_dd_ortho,
                                              1,
                                              ddpk_nl_fNL_4_dd_ortho,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out_ortho = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_4_dd_ortho,
                                                    ddpk_nl_fNL_4_dd_ortho,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out_ortho,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_4_dd_ortho[index_k] = pk_nl_fNL_out_ortho * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_4_dd_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL4dd*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_l_fNL_4_dd[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_4_dd[index_k]-0.*large_for_logs_fNL); //GC!

            //printf("%.20e %.20e\n",pnlpt->k[index_k],pk_l_fNL_4_dd_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            //GC: this is zero unless IR resummation... Notice that before I was summing the epsilon_for_logs_fNL also in the non-padded piece, when doing the checks. The final version sent to Misha of course does NOT do that. The version with orthogonal also does NOT since I copied the lines from the version with equilateral sent to Misha. You can see it here since indeed I output zeros... The checks against Mathematica were swamped by interpolation differences, and the C vs. C checks automatically had \epsilon - \epsilon... Now I just do not have them...
            //GC: same with the fact that now, correctly, I have to subtract the large_for_logs_fNL since I use code from the final version where it is needed...
        }

        //GC: ORTHOGONAL -- finish

        //GC: hexadecapole vd...

        double *ddpk_nl_fNL_4_vd;
        class_alloc(ddpk_nl_fNL_4_vd, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_4_vd,
                                              1,
                                              ddpk_nl_fNL_4_vd,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out = 0;

        /*
    
    double irfNL4vd; //GC -> COULD BE USELESS! ESPECIALLY SINCE I MOST LIKELY DO NOT HAVE A SIMPLE POWER LAW BUT I GET THE SQUARE ROOT OF THE POWER SPECTRUM! I WILL NEED TO CHECK THIS!!!
    class_call(array_interpolate_spline(kdisc,
                                        Nmax,
                                        P12_4_vd,
                                        ddpk_nl_fNL_4_vd,
                                        1,
                                        kminnew,
                                        &last_index,
                                        &irfNL4vd, //GC -> if I interpolate at high k, I will need a "uv" term as well...
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    last_index=0;
    pk_nl_fNL_out=0; //GC: redundant???
    
    */

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_4_vd,
                                                    ddpk_nl_fNL_4_vd,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_4_vd[index_k] = pk_nl_fNL_out * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_4_vd[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL4vd*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_l_fNL_4_vd[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_4_vd[index_k]-0.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- start

        double *ddpk_nl_fNL_4_vd_ortho;
        class_alloc(ddpk_nl_fNL_4_vd_ortho, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_4_vd_ortho,
                                              1,
                                              ddpk_nl_fNL_4_vd_ortho,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out_ortho = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_4_vd_ortho,
                                                    ddpk_nl_fNL_4_vd_ortho,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out_ortho,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_4_vd_ortho[index_k] = pk_nl_fNL_out_ortho * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_4_vd_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL4vd*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_l_fNL_4_vd[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_4_vd[index_k]-0.*large_for_logs_fNL); //GC!

            //printf("%.20f %.20f\n",pnlpt->k[index_k],pk_l_fNL_4_vd_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            //printf("%.20e %.20e\n",pnlpt->k[index_k],pk_l_fNL_4_vd_ortho[index_k]-1.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- finish

        //GC: hexadecapole vv...

        double *ddpk_nl_fNL_4_vv;
        class_alloc(ddpk_nl_fNL_4_vv, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_4_vv,
                                              1,
                                              ddpk_nl_fNL_4_vv,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out = 0;

        /*
    
    double irfNL4vv; //GC -> COULD BE USELESS! ESPECIALLY SINCE I MOST LIKELY DO NOT HAVE A SIMPLE POWER LAW BUT I GET THE SQUARE ROOT OF THE POWER SPECTRUM! I WILL NEED TO CHECK THIS!!!
    class_call(array_interpolate_spline(kdisc,
                                        Nmax,
                                        P12_4_vv,
                                        ddpk_nl_fNL_4_vv,
                                        1,
                                        kminnew,
                                        &last_index,
                                        &irfNL4vv, //GC -> if I interpolate at high k, I will need a "uv" term as well... And maybe like this I can remove the cutoff...
                                        1,
                                        pnlpt->error_message),
               pnlpt->error_message,
               pnlpt->error_message);

    last_index=0;
    pk_nl_fNL_out=0; //GC: redundant???
    
   */

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_4_vv,
                                                    ddpk_nl_fNL_4_vv,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_4_vv[index_k] = pk_nl_fNL_out * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_4_vv[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL4vv*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_l_fNL_4_vv[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_4_vv[index_k]-0.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- start

        double *ddpk_nl_fNL_4_vv_ortho;
        class_alloc(ddpk_nl_fNL_4_vv_ortho, sizeof(double) * Nmax, pnlpt->error_message);

        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P12_4_vv_ortho,
                                              1,
                                              ddpk_nl_fNL_4_vv_ortho,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_fNL_out_ortho = 0;

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P12_4_vv_ortho,
                                                    ddpk_nl_fNL_4_vv_ortho,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_nl_fNL_out_ortho,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);

                pk_l_fNL_4_vv_ortho[index_k] = pk_nl_fNL_out_ortho * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                //+1.e7; //GC: needs to be allocated!!!
            }

            else
            {
                pk_l_fNL_4_vv_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                //1.e7 + irfNL4vv*pow((pnlpt->k[index_k]/kminnew),4.)*exp(-pow(pnlpt->k[index_k]/3.,4.)); //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_l_fNL_4_vv[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_fNL_4_vv[index_k]-0.*large_for_logs_fNL); //GC!

            //printf("%.20f %.20f\n",pnlpt->k[index_k],pk_l_fNL_4_vv_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            //printf("%.20e %.20e\n",pnlpt->k[index_k],pk_l_fNL_4_vv_ortho[index_k]-1.*large_for_logs_fNL); //GC!
        }

        //GC: ORTHOGONAL -- finish

        //GC!

        double *ddpk_CTR_0;
        class_alloc(ddpk_CTR_0, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_CTR_0,
                                              1,
                                              ddpk_CTR_0,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_out = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

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

                if (pk_nl_out<=0) pk_nl_out = 1e-16;
                pk_CTR_0[index_k] = pk_nl_out;
            }
            else
            {
                pk_CTR_0[index_k] = exp(lnpk_l[index_k] + 2. * lnk_l[index_k]);
            }
        }

        double *ddpk_CTR_2;
        class_alloc(ddpk_CTR_2, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_CTR_2,
                                              1,
                                              ddpk_CTR_2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        pk_nl_out = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {
            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {
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
                if (pk_nl_out<=0) pk_nl_out = 1e-16;
                pk_CTR_2[index_k] = pk_nl_out;
            }
            else
            {
                pk_CTR_2[index_k] = exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * f * 2. / 3.;
                //    pk_CTR_2[index_k] = exp(lnpk_l[index_k]+2.*lnk_l[index_k])*(f*f*2./3.*2.2*2.2+8.*2.2*pow(f,3.)/7.+10.*pow(f,4.)/21.);
            }
        }

        double *ddpk_CTR_4;
        class_alloc(ddpk_CTR_4, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc, Nmax, P_CTR_4, 1, ddpk_CTR_4, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
        last_index = 0;
        pk_nl_out = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {
            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {
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
                if (pk_nl_out<=0) pk_nl_out = 1e-16;
                pk_CTR_4[index_k] = pk_nl_out;
            }

            else
            {
                pk_CTR_4[index_k] = exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * f * f * 8. / 35.;
            }
        }

        double *ddpk_tree_0_vv;
        class_alloc(ddpk_tree_0_vv, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc, Nmax, Ptree_0_vv, 1, ddpk_tree_0_vv, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
        last_index = 0;
        pk_nl_out = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {
            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {
                class_call(array_interpolate_spline(kdisc, Nmax, Ptree_0_vv, ddpk_tree_0_vv, 1, pnlpt->k[index_k], &last_index, &pk_nl_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                pk_Tree_0_vv[index_k] = pk_nl_out + large_for_logs_big;
            }
            else
            {
                pk_Tree_0_vv[index_k] = exp(lnpk_l[index_k]) * f * f / 5. + large_for_logs_big;
            }
        }
        free(ddpk_tree_0_vv);
        free(Ptree_0_vv);

        double *ddpk_tree_0_vd;
        class_alloc(ddpk_tree_0_vd, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc, Nmax, Ptree_0_vd, 1, ddpk_tree_0_vd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
        last_index = 0;
        pk_nl_out = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {
            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {
                class_call(array_interpolate_spline(kdisc, Nmax, Ptree_0_vd, ddpk_tree_0_vd, 1, pnlpt->k[index_k], &last_index, &pk_nl_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                pk_Tree_0_vd[index_k] = pk_nl_out + large_for_logs_big;
            }
            else
            {
                pk_Tree_0_vd[index_k] = exp(lnpk_l[index_k]) * 2. * f / 3. + large_for_logs_big;
            }
        }
        free(ddpk_tree_0_vd);
        free(Ptree_0_vd);

        double *ddpk_tree_0_dd;
        class_alloc(ddpk_tree_0_dd, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc, Nmax, Ptree_0_dd, 1, ddpk_tree_0_dd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
        last_index = 0;
        pk_nl_out = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {
            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {
                class_call(array_interpolate_spline(kdisc, Nmax, Ptree_0_dd, ddpk_tree_0_dd, 1, pnlpt->k[index_k], &last_index, &pk_nl_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                pk_Tree_0_dd[index_k] = pk_nl_out + large_for_logs_big;
            }
            else
            {
                pk_Tree_0_dd[index_k] = exp(lnpk_l[index_k]) + large_for_logs_big;
            }
        }
        free(ddpk_tree_0_dd);
        free(Ptree_0_dd);

        double *ddpk_tree_2_vv;
        class_alloc(ddpk_tree_2_vv, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc, Nmax, Ptree_2_vv, 1, ddpk_tree_2_vv, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
        last_index = 0;
        pk_nl_out = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {
            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {
                class_call(array_interpolate_spline(kdisc, Nmax, Ptree_2_vv, ddpk_tree_2_vv, 1, pnlpt->k[index_k], &last_index, &pk_nl_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                pk_Tree_2_vv[index_k] = pk_nl_out + large_for_logs_big;
            }
            else
            {
                pk_Tree_2_vv[index_k] = exp(lnpk_l[index_k]) * 4. * f * f / 7. + large_for_logs_big;
            }
        }
        free(ddpk_tree_2_vv);
        free(Ptree_2_vv);

        double *ddpk_tree_2_vd;
        class_alloc(ddpk_tree_2_vd, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc, Nmax, Ptree_2_vd, 1, ddpk_tree_2_vd, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
        last_index = 0;
        pk_nl_out = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {
            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {
                class_call(array_interpolate_spline(kdisc, Nmax, Ptree_2_vd, ddpk_tree_2_vd, 1, pnlpt->k[index_k], &last_index, &pk_nl_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                pk_Tree_2_vd[index_k] = pk_nl_out + large_for_logs_big;
            }
            else
            {
                pk_Tree_2_vd[index_k] = exp(lnpk_l[index_k]) * 4. * f / 3. + large_for_logs_big;
            }
        }
        free(ddpk_tree_2_vd);
        free(Ptree_2_vd);

        double *ddpk_tree_4_vv;
        class_alloc(ddpk_tree_4_vv, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc, Nmax, Ptree_4_vv, 1, ddpk_tree_4_vv, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
        last_index = 0;
        pk_nl_out = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {
            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {
                class_call(array_interpolate_spline(kdisc, Nmax, Ptree_4_vv, ddpk_tree_4_vv, 1, pnlpt->k[index_k], &last_index, &pk_nl_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                pk_Tree_4_vv[index_k] = pk_nl_out + large_for_logs_big;
            }
            else
            {
                pk_Tree_4_vv[index_k] = exp(lnpk_l[index_k]) * 8. * f * f / 35. + large_for_logs_big;
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

        //GC! Again strange order below :)

        free(ddpk_nl_fNL_2_vv);
        free(ddpk_nl_fNL_2_vd);
        free(ddpk_nl_fNL_2_dd);

        free(ddpk_nl_fNL_4_vv);
        free(ddpk_nl_fNL_4_vd);
        free(ddpk_nl_fNL_4_dd);

        free(ddpk_nl_fNL_0_vv);
        free(ddpk_nl_fNL_0_vd);
        free(ddpk_nl_fNL_0_dd);

        //GC: ORTHOGONAL -- start

        free(ddpk_nl_fNL_2_vv_ortho);
        free(ddpk_nl_fNL_2_vd_ortho);
        free(ddpk_nl_fNL_2_dd_ortho);

        free(ddpk_nl_fNL_4_vv_ortho);
        free(ddpk_nl_fNL_4_vd_ortho);
        free(ddpk_nl_fNL_4_dd_ortho);

        free(ddpk_nl_fNL_0_vv_ortho);
        free(ddpk_nl_fNL_0_vd_ortho);
        free(ddpk_nl_fNL_0_dd_ortho);

        //GC: ORTHOGONAL -- finish

        //GC: no dd hexadecapole??? WHY WOULD THAT BE, after AP? No, everything is there... WAIT: these are defined as e.g. pk_l_fNL_4_vv at the end. Where does... No problem, a quick search for P12_4_dd shows that it was defined already... It was because around line 3900 there was "//free(f12_4_dd);" -> the point is that BEFORE AP I do not have it... This must be checked: see above...
        free(P12_4_dd);
        free(P12_4_vd);
        free(P12_4_vv);
        free(P12_2_dd);
        free(P12_2_vd);
        free(P12_2_vv);
        free(P12_0_dd);
        free(P12_0_vd);
        free(P12_0_vv);

        //GC: ORTHOGONAL -- start

        free(P12_4_dd_ortho);
        free(P12_4_vd_ortho);
        free(P12_4_vv_ortho);
        free(P12_2_dd_ortho);
        free(P12_2_vd_ortho);
        free(P12_2_vv_ortho);
        free(P12_0_dd_ortho);
        free(P12_0_vd_ortho);
        free(P12_0_vv_ortho);

        //GC: ORTHOGONAL -- finish

        //GC!

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

    } // end of RSD conditional expression

    free(f13);
    free(x);
    free(x_transfer); //GC! Essentially everything that was allocated...
    free(x_w);
    free(x_w_transfer); //GC! RECALL THAT I STILL HAVE TO USE IT CORRECTLY!!! y is not freed, so I do not need to free y_transfer. It is freed at the end...
    free(f22);
    free(f12); //GC!
    free(Pdisc);
    free(PPRIMdisc); //GC!
    free(ddpk_nl);
    free(ddpk_nl_fNL); //GC!
    //GC: ORTHOGONAL -- start
    free(ddpk_nl_fNL_ortho); //GC!
    //GC: ORTHOGONAL -- finish
    free(ddpk_CTR);
    free(ddpk_Tree);
    free(etam);
    //free(etam_transfer); //GC: I can keep using them for G2!!! Since there I have the same bias...
    free(cmsym);
    //free(cmsym_transfer); //GC: I can keep using them for G2!!! Since there I have the same bias...
    free(P22);
    free(P13);
    //free(P12); //GC: CAREFUL!!! THIS HAS NOTHING TO DO WITH MY THINGS!!!
    free(P12_fNL); //GC -> recall name confusions...
    //GC: ORTHOGONAL -- start
    free(P12_fNL_ortho); //GC -> recall name confusions...
    //GC: ORTHOGONAL -- finish
    free(P13UV);
    free(P_CTR);
    free(P1loop);
    free(Ptree);
    free(myddlnpk);

    /* Computing the power spectra for biased tracers. For this reason we have to compute the FFTLog coefficients for a new 'bias' exponent b2 */

    if (pnlpt->bias == bias_yes)
    {

        if (pnlpt->nonlinear_pt_verbose > 0)
            printf("Computing the spectra for biased tracers...\n");

        double complex *etam2;
        class_alloc(etam2, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

        double complex *etam2_transfer;
        class_alloc(etam2_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

        double b2 = -1.6000001;
        double b2_transfer = -1.25; //GC: why the 0.00...001?
        int index_c2 = 0;
        for (index_c2 = 0; index_c2 < Nmax + 1; index_c2++)
        {
            etam2[index_c2] = b2 + 2. * M_PI * _Complex_I * js[index_c2] / Nmaxd / Delta;
            etam2_transfer[index_c2] = b2_transfer + 2. * M_PI * _Complex_I * js[index_c2] / Nmaxd / Delta; //GC: wait -> it seems that for biased tracers he takes -1.6 ALSO for e.g. G2G2? In that case, it was ok to deallocate etam and cmsym. BUT now I only have a new bias for \delta^2... So I can avoid deallocating the ones above and keep using them...
        }

        double *input_real_bias;
        class_alloc(input_real_bias, (Nmax) * sizeof(double), pnlpt->error_message);
        double *input_imag_bias;
        class_alloc(input_imag_bias, (Nmax) * sizeof(double), pnlpt->error_message);
        double *output_real_bias;
        class_alloc(output_real_bias, (Nmax) * sizeof(double), pnlpt->error_message);
        double *output_imag_bias;
        class_alloc(output_imag_bias, (Nmax) * sizeof(double), pnlpt->error_message);

        //GC!
        double *input_real_transfer2;
        class_alloc(input_real_transfer2, (Nmax) * sizeof(double), pnlpt->error_message);
        double *input_imag_transfer2;
        class_alloc(input_imag_transfer2, (Nmax) * sizeof(double), pnlpt->error_message);
        double *output_real_transfer2;
        class_alloc(output_real_transfer2, (Nmax) * sizeof(double), pnlpt->error_message);
        double *output_imag_transfer2;
        class_alloc(output_imag_transfer2, (Nmax) * sizeof(double), pnlpt->error_message);
        //GC!

        index_c2 = 0;
        for (index_c2 = 0; index_c2 < Nmax; index_c2++)
        {
            input_real_bias[index_c2] = Pbin[index_c2] * exp(-1. * index_c2 * b2 * Delta);
            input_imag_bias[index_c2] = 0.;
            //GC!
            input_real_transfer2[index_c2] = Tbin[index_c2] * exp(-1. * index_c2 * b2_transfer * Delta);
            input_imag_transfer2[index_c2] = 0.;
        }

        FFT(input_real_bias, input_imag_bias, output_real_bias, output_imag_bias, Nmax, stepsize);
        FFT(input_real_transfer2, input_imag_transfer2, output_real_transfer2, output_imag_transfer2, Nmax, stepsize); //GC!

        //GC - SWITCH ~> as discussed before, these are kept...
        
        double complex *cmsym2;
        class_alloc(cmsym2, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

        double complex *cmsym2_transfer;
        class_alloc(cmsym2_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

        //     double complex cmsym[Nmax+1];
        index_c2 = 0;
        for (index_c2 = 0; index_c2 < Nmax + 1; index_c2++)
        {
            if (index_c2 < Nmax / 2)
            {
                cmsym2[index_c2] = cpow(kmin, -etam2[index_c2]) * (output_real_bias[Nmax / 2 - index_c2] - _Complex_I * output_imag_bias[Nmax / 2 - index_c2]) / Nmaxd;

                cmsym2_transfer[index_c2] = cpow(kmin, -etam2_transfer[index_c2]) * (output_real_transfer2[Nmax / 2 - index_c2] - _Complex_I * output_imag_transfer2[Nmax / 2 - index_c2]) / Nmaxd; //GC!
            }
            else
            {

                cmsym2[index_c2] = cpow(kmin, -etam2[index_c2]) * (output_real_bias[index_c2 - Nmax / 2] + _Complex_I * output_imag_bias[index_c2 - Nmax / 2]) / Nmaxd; //GC -> careful, I had modified this without copying it...

                cmsym2_transfer[index_c2] = cpow(kmin, -etam2_transfer[index_c2]) * (output_real_transfer2[index_c2 - Nmax / 2] + _Complex_I * output_imag_transfer2[index_c2 - Nmax / 2]) / Nmaxd; //GC! CAREFUL, I had "cancelled" the one above...
            }
        }

        cmsym2[0] = cmsym2[0] / 2.;
        cmsym2[Nmax] = cmsym2[Nmax] / 2.;

        //GC!

        cmsym2_transfer[0] = cmsym2_transfer[0] / 2.;
        cmsym2_transfer[Nmax] = cmsym2_transfer[Nmax] / 2.;

        free(input_real_bias);
        free(input_imag_bias);
        free(output_real_bias);
        free(output_imag_bias);

        //GC!

        free(input_real_transfer2);
        free(input_imag_transfer2);
        free(output_real_transfer2);
        free(output_imag_transfer2);

        /* Computing Id2d2 */

        double complex *f22_Id2d2;
        double *P_Id2d2;
        class_alloc(f22_Id2d2, Nmax * sizeof(complex double), pnlpt->error_message);
        class_alloc(P_Id2d2, Nmax * sizeof(double), pnlpt->error_message);

        double epsilon_for_logs = 1. * pow(10., -6.);

        int count2 = 0;
        index_i = 0;
        index_l = 0;
        index_j = 0;
        count = 0;

        double complex *x2;
        class_alloc(x2, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

        for (index_j = 0; index_j < Nmax; index_j++)
        {

            for (count = 0; count < Nmax + 1; count++)
            {
                x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
            }
            zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22basic_oneline_complex, x2, &inc, &beta, y, &inc);
            f22_Id2d2[index_j] = 2. * zdotu_(&Nmaxf, x2, &inc, y, &inc);
            //GC -> exactly ~> d^2d^2, i.e. the P22 of this, has I(n1,n2) times 2 from two permutations! Perfect. So this is why in the python, the Md2d2 is the one extracting I(n1,n2) WITHOUT the overall 2... I will do somewhat differently. First, I have no simple factor of 2 since there is always the shape. So I just always include the 2 in the matrices, and then they are the correct matrices BY THEMSELVES to multiply by b2/2 and bG2. Then, since I have only two of them, I just export them from python by themselves as M12delta and M12G2 or something, and use them here directly in this loop, being careful about biases (which I can define from scratch), etc.

            P_Id2d2[index_j] = fabs(creal(cpow(kdisc[index_j], 3.) * f22_Id2d2[index_j]) - creal(cpow(kdisc[0], 3.) * f22_Id2d2[0]) + epsilon_for_logs);
        }

        //GC -> I cannot use simply M22basic_oneline_complex (and multiply by 2 because of the permutations). Of course the permutations are included in the matrices I will write now [WILL RECHECK, but I checked already a billion times]... Notice that however I will have to be careful that now I will use pnlpt->something[abc], but this could cause problems since these matrices are not "preallocated" in that structure?

        /*  Having eta2 we compute the rest of the PT matrices */

        count2 = 0;
        for (index_l = 0; index_l < Nmax + 1; index_l++)
        {
            for (index_i = index_l; index_i < Nmax + 1; index_i++)
            {
                pnlpt->M_IG2G2[count2] = pnlpt->M22basic_oneline_complex[count2] * ((3. + etam2[index_i] + etam2[index_l]) * (1. + etam2[index_i] + etam2[index_l]) / ((-0.5 * etam2[index_i]) * (-0.5 * etam2[index_l]) * (1. - 0.5 * etam2[index_i]) * (1. - 0.5 * etam2[index_l])));

                pnlpt->M_Id2[count2] = pnlpt->M22basic_oneline_complex[count2] * ((3. + etam2[index_i] + etam2[index_l]) * (4. + 3.5 * (etam2[index_i] + etam2[index_l])) / (14. * (-0.5 * etam2[index_l]) * (-0.5 * etam2[index_i])));

                pnlpt->M_Id2G2[count2] = pnlpt->M22basic_oneline_complex[count2] * ((3. + etam2[index_i] + etam2[index_l]) / ((-0.5 * etam2[index_i]) * (-0.5 * etam2[index_l])));

                pnlpt->M_IG2[count2] = pnlpt->M22basic_oneline_complex[count2] * (-1. * (3. + etam2[index_i] + etam2[index_l]) * (1. + etam2[index_i] + etam2[index_l]) * (6. - 3.5 * (etam2[index_i] + etam2[index_l])) / (28. * (1. - 0.5 * etam2[index_i]) * (1. - 0.5 * etam2[index_l]) * (-0.5 * etam2[index_i]) * (-0.5 * etam2[index_l])));

                if (SWITCH_index == 1) {
                
                nu1 = -0.5 * etam2_transfer[index_i];
                nu2 = -0.5 * etam2_transfer[index_l]; //GC! IMPORTANT TO USE THE CORRECT \eta, one has one bias and one another... I use the fact that the second also has -0.8, as the matter ones...

                //pnlpt->M_fNLd2[count2] = pnlpt->M22basic_oneline_complex[count2] * (1./(csin(nu1*M_PI))/(csin(nu2*M_PI))/(csin((nu1 + nu2)*M_PI))*((-15. - 8.*nu1*(-1. + nu2) - 8.*(-3. + nu2)*nu2)*csin(2.*nu1*M_PI) + (-15. + 8.*nu2 - 8.*nu1*(-3. + nu1 + nu2))*csin(2.*nu2*M_PI) + (-9. + 8.*nu1 + 8.*nu2 - 8.*nu1*nu2)*csin(2.*(nu1 + nu2)*M_PI)))/((-1. + 2.*nu1)*(-1. + 2.*nu2)*(-5. + 2.*nu1 + 2.*nu2));

                pnlpt->M_fNLd2[count2] = pnlpt->M12_oneline_complex_bias_real_space_b2[count2];

                //GC: ORTHOGONAL -- start

                pnlpt->M_fNLd2_ortho[count2] = pnlpt->M12_oneline_complex_bias_real_space_b2_ortho[count2];

                //GC: ORTHOGONAL -- finish

                nu1 = -0.5 * etam_transfer[index_i];
                nu2 = -0.5 * etam_transfer[index_l]; //GC!

                //pnlpt->M_fNLG2[count2] = pnlpt->M22basic_oneline_complex[count2] * (1./(csin(nu1*M_PI))/(csin(nu2*M_PI))*((-30. + 96.*nu2 + 16.*(-3.*(-2. + nu1)*nu1 + nu1*(-43. + 2.*(25. - 8.*nu1)*nu1)*nu2 + (-3. + 2.*nu1*(-5. + 2.*nu1)*(-5. + 4.*nu1))*nu2*nu2 + 16.*(-1. + nu1)*nu1*nu2*nu2*nu2))*ccos((nu1 + nu2)*M_PI) + 1./(csin((nu1 + nu2)*M_PI))*((15. + 128.*nu1*nu1*nu1*(-1. + nu2)*nu2 - 8.*(-2. + nu2)*nu2*(1. + 4.*(-2. + nu2)*nu2) + 8.*nu1*(-3. + 2.*nu2)*(2. + nu2*(11. + 8.*(-3. + nu2)*nu2)) + 8.*nu1*nu1*(3. + 2.*nu2*(23. + 2.*nu2*(-21. + 8.*nu2))))*csin(2.*nu1*M_PI) + 15.*csin(2.*nu2*M_PI) + 8.*((2. - 1.*nu1)*nu1*(1. + 4.*(-2. + nu1)*nu1) + (-3. + 2.*nu1)*(2. + nu1*(11. + 8.*(-3. + nu1)*nu1))*nu2 + (3. + 2.*nu1*(23. + 2.*nu1*(-21. + 8.*nu1)))*nu2*nu2 + 16.*(-1. + nu1)*nu1*nu2*nu2*nu2)*csin(2.*nu2*M_PI))))/(2.*(nu1)*(-1. + 4.*nu1*nu1)*(nu2)*(-5. + 2.*nu1 + 2.*nu2)*(-1. + 4.*nu2*nu2));

                pnlpt->M_fNLG2[count2] = pnlpt->M12_oneline_complex_bias_real_space_bG2[count2];

                //GC: ORTHOGONAL -- start

                pnlpt->M_fNLG2_ortho[count2] = pnlpt->M12_oneline_complex_bias_real_space_bG2_ortho[count2];

                //GC: ORTHOGONAL -- finish

                //GC!
                    
                }

                count2++;
            }
        }

        //if (pnlpt->rsd == rsd_yes && pnlpt->rsd_only == rsd_only_no || pnlpt->rsd == rsd_no) {

        /* //GC -> I put the P12 stuff before, this time... */

        //GC: contribution from \delta^2...

        double complex *f12_fNLd2;
        double *P_fNLd2;
        class_alloc(f12_fNLd2, Nmax * sizeof(complex double), pnlpt->error_message);
        class_alloc(P_fNLd2, Nmax * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        double *P_fNLd2_ortho;
        class_alloc(P_fNLd2_ortho, Nmax * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        double complex *x2_transfer;
        class_alloc(x2_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);
        //GC: since these were constantly overwritten it is fine to add another here, that will be used for BOTH BIAS TERMS!!! EVEN IF ONE HAS THE SAME BIAS AS THE OTHER TERMS...

        index_j = 0;
        count = 0;
        
        if (SWITCH_index == 1) {

        for (index_j = 0; index_j < Nmax; index_j++)
        {

            for (count = 0; count < Nmax + 1; count++)
            {
                x2_transfer[count] = cmsym2_transfer[count] * cpow(kdisc[index_j], etam2_transfer[count]); //GC: this is the one of \delta squared, so it is correct that it must have the "2"...
            }

            zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M_fNLd2, x2_transfer, &inc, &beta, y_transfer, &inc); //GC: notice the use of y_transfer!
            f12_fNLd2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
            P_fNLd2[index_j] = Tbin[index_j] * //[FACTOR] *
                               creal(cpow(kdisc[index_j], 3) * f12_fNLd2[index_j]);

            //printf("%.24f %.24f %.24f",kdisc[index_j],Tbin[index_j],P_fNLd2[index_j]); //GC!
            //printf("\n"); //GC!

            //GC: ORTHOGONAL -- start

            zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M_fNLd2_ortho, x2_transfer, &inc, &beta, y_transfer, &inc); //GC: notice the use of y_transfer!
            f12_fNLd2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
            P_fNLd2_ortho[index_j] = Tbin[index_j] * //[FACTOR] *
                                     creal(cpow(kdisc[index_j], 3) * f12_fNLd2[index_j]);

            //GC: ORTHOGONAL -- finish

            //printf("%.24f %.24f %.24f",kdisc[index_j],Tbin[index_j],P_fNLd2_ortho[index_j]); //GC!
            //printf("\n"); //GC!
        }
            
        }

        //GC: contribution from G^2...

        double complex *f12_fNLG2;
        double *P_fNLG2;
        class_alloc(f12_fNLG2, Nmax * sizeof(complex double), pnlpt->error_message);
        class_alloc(P_fNLG2, Nmax * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- start

        double *P_fNLG2_ortho;
        class_alloc(P_fNLG2_ortho, Nmax * sizeof(double), pnlpt->error_message);

        //GC: ORTHOGONAL -- finish

        //double complex *x2_transfer;
        //class_alloc(x2_transfer,(Nmax+1)*sizeof(complex double),pnlpt->error_message);
        //GC: already allocated above: now this just gets overwritten as it was before...

        index_j = 0;
        count = 0;

        if (SWITCH_index == 1) {
        
        for (index_j = 0; index_j < Nmax; index_j++)
        {

            for (count = 0; count < Nmax + 1; count++)
            {
                x2_transfer[count] = cmsym_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]); //GC: this is the one of G2, so it uses the old coefficients...
            }

            zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M_fNLG2, x2_transfer, &inc, &beta, y_transfer, &inc); //GC: notice the use of y_transfer!
            f12_fNLG2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
            P_fNLG2[index_j] = Tbin[index_j] * //[FACTOR] *
                               creal(cpow(kdisc[index_j], 3) * f12_fNLG2[index_j]);

            //printf("%.24f %.24f %.24f",kdisc[index_j],Tbin[index_j],P_fNLG2[index_j]); //GC!
            //printf("\n"); //GC!

            //GC: ORTHOGONAL -- start

            zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M_fNLG2_ortho, x2_transfer, &inc, &beta, y_transfer, &inc); //GC: notice the use of y_transfer!
            f12_fNLG2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
            P_fNLG2_ortho[index_j] = Tbin[index_j] * //[FACTOR] *
                                     creal(cpow(kdisc[index_j], 3) * f12_fNLG2[index_j]);

            //GC: ORTHOGONAL -- finish

            //printf("%.24f %.24f %.24f",kdisc[index_j],Tbin[index_j],P_fNLG2_ortho[index_j]); //GC!
            //printf("\n"); //GC!
        }
            
        }

        //NOTHING ELSE TO DO HERE!!!

        /* Computing Id2 */

        double complex *f22_Id2;
        double *P_Id2;
        class_alloc(f22_Id2, Nmax * sizeof(complex double), pnlpt->error_message);
        class_alloc(P_Id2, Nmax * sizeof(double), pnlpt->error_message);

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

        for (index_j = 0; index_j < Nmax; index_j++)
        {

            for (count = 0; count < Nmax + 1; count++)
            {
                x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
            }

            zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M_Id2, x2, &inc, &beta, y, &inc);
            f22_Id2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
            P_Id2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_Id2[index_j]);
        }

        /* Computing IG2 */

        double complex *f22_IG2;
        double *P_IG2;
        class_alloc(f22_IG2, Nmax * sizeof(complex double), pnlpt->error_message);
        class_alloc(P_IG2, Nmax * sizeof(double), pnlpt->error_message);

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

        for (index_j = 0; index_j < Nmax; index_j++)
        {
            for (count = 0; count < Nmax + 1; count++)
            {
                x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
            };
            zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M_IG2, x2, &inc, &beta, y, &inc);
            f22_IG2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
            //printf("f22_IG2_real=%.18le f22_IG2_imag=%.18le\n",creal(f22_IG2[index_j]),cimag(f22_IG2[index_j]));

            P_IG2[index_j] = fabs(creal(cpow(kdisc[index_j], 3) * f22_IG2[index_j]));
            //printf("%le %le\n",kdisc[j],P_Id2[j]);
        }
        
        //GC - SWITCH ~> as discussed, these are kept...

        double *ddpk_PfNLd2;
        double pk_PfNLd2_out = 0.;
        class_alloc(ddpk_PfNLd2, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_fNLd2,
                                              1,
                                              ddpk_PfNLd2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmax && pnlpt->k[index_k] >= kmin)
            {
                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P_fNLd2,
                                                    ddpk_PfNLd2,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_PfNLd2_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);
                pk_fNLd2[index_k] =                     //10. + //GC!
                    pk_PfNLd2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC -> NEED TO CHECK THE CONSTANTS HERE!!! The outputs here are allocated around line 2000...
            }

            else
            {
                pk_fNLd2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //10.; //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_fNLd2[index_k]-0.*large_for_logs_fNL); //GC!

            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_fNLd2[index_k]-0.*large_for_logs_fNL); //GC!
        }
        free(ddpk_PfNLd2);

        free(f12_fNLd2);
        free(P_fNLd2);

        //GC: ORTHOGONAL -- start

        double *ddpk_PfNLd2_ortho;
        double pk_PfNLd2_out_ortho = 0.;
        class_alloc(ddpk_PfNLd2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_fNLd2_ortho,
                                              1,
                                              ddpk_PfNLd2_ortho,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmax && pnlpt->k[index_k] >= kmin)
            {
                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P_fNLd2_ortho,
                                                    ddpk_PfNLd2_ortho,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_PfNLd2_out_ortho,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);
                pk_fNLd2_ortho[index_k] =                     //10. + //GC!
                    pk_PfNLd2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC -> NEED TO CHECK THE CONSTANTS HERE!!! The outputs here are allocated around line 2000...
            }

            else
            {
                pk_fNLd2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //10.; //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_fNLd2[index_k]-0.*large_for_logs_fNL); //GC!

            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_fNLd2[index_k]-0.*large_for_logs_fNL); //GC!

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_fNLd2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            //printf("%.24e %.24e\n",pnlpt->k[index_k],pk_fNLd2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
        }
        free(ddpk_PfNLd2_ortho);

        //free(f12_fNLd2);
        free(P_fNLd2_ortho);

        //GC: ORTHOGONAL -- finish

        double *ddpk_PfNLG2;
        double pk_PfNLG2_out = 0.;
        class_alloc(ddpk_PfNLG2, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_fNLG2,
                                              1,
                                              ddpk_PfNLG2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmax && pnlpt->k[index_k] >= kmin)
            {
                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P_fNLG2,
                                                    ddpk_PfNLG2,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_PfNLG2_out,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);
                pk_fNLG2[index_k] =                     //10. + //GC!
                    pk_PfNLG2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC -> NEED TO CHECK THE CONSTANTS HERE!!! The outputs here are allocated around line 2000...
            }

            else
            {
                pk_fNLG2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //10.; //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_fNLG2[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_fNLG2[index_k]-0.*large_for_logs_fNL); //GC!
        }
        free(ddpk_PfNLG2);

        free(f12_fNLG2);
        free(P_fNLG2);

        //GC: ORTHOGONAL -- start

        double *ddpk_PfNLG2_ortho;
        double pk_PfNLG2_out_ortho = 0.;
        class_alloc(ddpk_PfNLG2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_fNLG2_ortho,
                                              1,
                                              ddpk_PfNLG2_ortho,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmax && pnlpt->k[index_k] >= kmin)
            {
                class_call(array_interpolate_spline(kdisc,
                                                    Nmax,
                                                    P_fNLG2_ortho,
                                                    ddpk_PfNLG2_ortho,
                                                    1,
                                                    pnlpt->k[index_k],
                                                    &last_index,
                                                    &pk_PfNLG2_out_ortho,
                                                    1,
                                                    pnlpt->error_message),
                           pnlpt->error_message,
                           pnlpt->error_message);
                pk_fNLG2_ortho[index_k] =                     //10. + //GC!
                    pk_PfNLG2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC -> NEED TO CHECK THE CONSTANTS HERE!!! The outputs here are allocated around line 2000...
            }

            else
            {
                pk_fNLG2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //10.; //GC!
            }

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_fNLG2[index_k]-0.*large_for_logs_fNL); //GC!
            //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_fNLG2[index_k]-0.*large_for_logs_fNL); //GC!

            //printf("%.24f %.24f\n",pnlpt->k[index_k],pk_fNLG2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            //printf("%.24e %.24e\n",pnlpt->k[index_k],pk_fNLG2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
        }
        free(ddpk_PfNLG2_ortho);

        //free(f12_fNLG2);
        free(P_fNLG2_ortho);

        //GC: ORTHOGONAL -- finish

        //GC: DONE HERE!!! Now we go to RSDs...

        double *ddpk_PId2;
        double pk_Id2_out = 0.;
        class_alloc(ddpk_PId2, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_Id2,
                                              1,
                                              ddpk_PId2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmax && pnlpt->k[index_k] >= kmin)
            {
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
                pk_Id2[index_k] = large_for_logs_small + pk_Id2_out;
            }

            else
            {
                pk_Id2[index_k] = large_for_logs_small;
            }
        }
        free(ddpk_PId2);

        free(f22_Id2);
        free(P_Id2);

        double *ddpk_IG2;
        double pk_IG2_out = 0.;
        class_alloc(ddpk_IG2, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_IG2,
                                              1,
                                              ddpk_IG2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmax && pnlpt->k[index_k] >= kmin)
            {
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
                pk_IG2[index_k] = pk_IG2_out + large_for_logs_small;
            }
            else
            {
                pk_IG2[index_k] = epsilon_for_logs + large_for_logs_small;
            }
        }
        free(ddpk_IG2);

        free(f22_IG2);
        free(P_IG2);

        //  }// end of RSD_only condition

        /* Computing Id2G2 */

        double complex *f22_Id2G2;
        double *P_Id2G2;
        class_alloc(f22_Id2G2, Nmax * sizeof(complex double), pnlpt->error_message);
        class_alloc(P_Id2G2, Nmax * sizeof(double), pnlpt->error_message);

        index_j = 0;
        /*
        count2 = 0;
        for (index_l=0; index_l < Nmax+1; index_l++){
            for (index_i=index_l; index_i < Nmax+1; index_i++){
                pnlpt->M_Id2G2[count2] = pnlpt->M22basic_oneline_complex[count2] * ((3.+etam2[index_i]+etam2[index_l])/((-0.5*etam2[index_i])*(-0.5*etam2[index_l])));
                count2++;
            }
        }*/

        for (index_j = 0; index_j < Nmax; index_j++)
        {

            for (count = 0; count < Nmax + 1; count++)
            {
                x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
            }
            zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M_Id2G2, x2, &inc, &beta, y, &inc);
            f22_Id2G2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);

            P_Id2G2[index_j] = fabs(creal(cpow(kdisc[index_j], 3.) * f22_Id2G2[index_j]));
        }

        /* Computing IG2G2 */

        double complex *f22_IG2G2;
        double *P_IG2G2;
        class_alloc(f22_IG2G2, Nmax * sizeof(complex double), pnlpt->error_message);
        class_alloc(P_IG2G2, Nmax * sizeof(double), pnlpt->error_message);

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

        for (index_j = 0; index_j < Nmax; index_j++)
        {

            for (count = 0; count < Nmax + 1; count++)
            {
                x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
            }
            zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M_IG2G2, x2, &inc, &beta, y, &inc);
            f22_IG2G2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);

            P_IG2G2[index_j] = fabs(creal(cpow(kdisc[index_j], 3.) * f22_IG2G2[index_j]));
        }

        /* Computing IFG2 */

        double complex *f13_IFG2;
        double *P_IFG2;
        class_alloc(f13_IFG2, Nmax * sizeof(complex double), pnlpt->error_message);
        class_alloc(P_IFG2, Nmax * sizeof(double), pnlpt->error_message);

        for (index_j = 0; index_j < Nmax; index_j++)
        {
            /*f13_IFG2[index_j] = 0.;
         for (index_i=0; index_i < Nmax+1; index_i++){
             f13_IFG2[index_j] += cmsym2[index_i]* cpow(kdisc[index_j], etam2[index_i]) * (pnlpt->IFG2_oneline[index_i]+ _Complex_I * pnlpt->IFG2_oneline[index_i + Nmax+1]);
         }*/

            for (count = 0; count < Nmax + 1; count++)
            {
                x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
            };

            f13_IFG2[index_j] = zdotu_(&Nmaxf, x2, &inc, pnlpt->IFG2_oneline_complex, &inc);
            //         printf("f13_IFG2_real=%le f13_IFG2_imag=%le\n",creal(f13_IFG2[index_j]),cimag(f13_IFG2[index_j]));
            P_IFG2[index_j] = fabs(creal(cpow(kdisc[index_j], 3.) * f13_IFG2[index_j] * Pbin[index_j]));
        }

        double *ddpk_IFG2;
        double pk_IFG2_out = 0.;
        class_alloc(ddpk_IFG2, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_IFG2,
                                              1,
                                              ddpk_IFG2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {
            if (pnlpt->k[index_k] <= kmax && pnlpt->k[index_k] >= kmin)
            {
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
                pk_IFG2[index_k] = pk_IFG2_out + large_for_logs_small;
            }
            else
            {
                pk_IFG2[index_k] = epsilon_for_logs + large_for_logs_small;
            }
        }

        if (pnlpt->rsd == rsd_yes)
        {

            //GC!

            /* Computing vd monopole from \delta^2 in RSD -> of course as always this comes with b_2[~~~/2~~~]... */

            double complex *f12_0_b1b2;
            double *P12_0_b1b2;
            class_alloc(f12_0_b1b2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P12_0_b1b2, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- start

            double *P12_0_b1b2_ortho;
            class_alloc(P12_0_b1b2_ortho, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- finish

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;

            if (SWITCH_index == 1) {
            
            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam2_transfer[index_i];
                    nu2 = -0.5 * etam2_transfer[index_l];

                    //pnlpt->M12_0_b1b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-1./(csin(nu1*M_PI))/(csin(nu2*M_PI))/(csin((nu1 + nu2)*M_PI))*((15. + 8.*nu1*(-1. + nu2) - 24.*nu2 + 8.*nu2*nu2)*csin(2.*nu1*M_PI) + (15. + 8.*nu1*nu1 + 8.*nu1*(-3. + nu2) - 8.*nu2)*csin(2.*nu2*M_PI) + (9. + 8.*nu1*(-1. + nu2) - 8.*nu2)*csin(2.*(nu1 + nu2)*M_PI)))/((-1. + 2.*nu1)*(-1. + 2.*nu2)*(-5. + 2.*nu1 + 2.*nu2));

                    pnlpt->M12_0_b1b2_oneline_complex[count2] = 3. * (pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1[count2]);

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_0_b1b2_oneline_complex_ortho[count2] = 3. * (pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1_ortho[count2]);

                    //GC: ORTHOGONAL -- finish

                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2_transfer[count] = cmsym2_transfer[count] * cpow(kdisc[index_j], etam2_transfer[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_0_b1b2_oneline_complex, x2_transfer, &inc, &beta, y_transfer, &inc);
                f12_0_b1b2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
                P12_0_b1b2[index_j] = Tbin[index_j] *                                       //[FACTOR] *
                                      creal(cpow(kdisc[index_j], 3) * f12_0_b1b2[index_j]); //GC: recheck that indeed it is correct that I multiply by "Tbin[index_j] * [FACTOR]" and a k^3 only. 99.99999% sure...

                //printf("%.24f %.24f %.24f",kdisc[index_j],Tbin[index_j],P12_0_b1b2[index_j]); //GC!
                //printf("%.16e %.16e",kdisc[index_j],P12_0_b1b2[index_j]); //GC! USELESS: see README_about_bias_checks.txt...
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- start

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_0_b1b2_oneline_complex_ortho, x2_transfer, &inc, &beta, y_transfer, &inc);
                f12_0_b1b2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
                P12_0_b1b2_ortho[index_j] = Tbin[index_j] *                                       //[FACTOR] *
                                            creal(cpow(kdisc[index_j], 3) * f12_0_b1b2[index_j]); //GC: recheck that indeed it is correct that I multiply by "Tbin[index_j] * [FACTOR]" and a k^3 only. 99.99999% sure...

                //GC: ORTHOGONAL -- finish

                //printf("%.24e %.24e %.24e",kdisc[index_j],Tbin[index_j],P12_0_b1b2_ortho[index_j]); //GC!
                //printf("\n"); //GC!
            }
                
            }

            //GC!

            //GC: now comes the vv part...

            double complex *f12_0_b2;
            double *P12_0_b2;
            class_alloc(f12_0_b2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P12_0_b2, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- start

            double *P12_0_b2_ortho;
            class_alloc(P12_0_b2_ortho, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- finish

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;
            
            if (SWITCH_index == 1) {

            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam2_transfer[index_i];
                    nu2 = -0.5 * etam2_transfer[index_l];

                    //GC: here I still need the etam2_transfer!!!

                    //pnlpt->M12_0_b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-1.0*(f/(csin(nu1*M_PI))/(csin(nu2*M_PI))/(csin((nu1 + nu2)*M_PI))*((15. + 8.*nu1*(-1. + nu2) - 24.*nu2 + 8.*nu2*nu2)*csin(2.*nu1*M_PI) + (15. + 8.*nu1*nu1 + 8.*nu1*(-3. + nu2) - 8.*nu2)*csin(2.*nu2*M_PI) + (9. + 8.*nu1*(-1. + nu2) - 8.*nu2)*csin(2.*(nu1 + nu2)*M_PI)))/(3.*(-1. + 2.*nu1)*(-1. + 2.*nu2)*(-5. + 2.*nu1 + 2.*nu2)));

                    pnlpt->M12_0_b2_oneline_complex[count2] = f * (pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1[count2]);

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_0_b2_oneline_complex_ortho[count2] = f * (pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1_ortho[count2]);

                    //GC: ORTHOGONAL -- finish

                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2_transfer[count] = cmsym2_transfer[count] * cpow(kdisc[index_j], etam2_transfer[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_0_b2_oneline_complex, x2_transfer, &inc, &beta, y_transfer, &inc);
                f12_0_b2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
                P12_0_b2[index_j] = Tbin[index_j] *                                     //[FACTOR] *
                                    creal(cpow(kdisc[index_j], 3) * f12_0_b2[index_j]); //GC: recheck that indeed it is correct that I multiply by "Tbin[index_j] * [FACTOR]" and a k^3 only. 99.99999% sure...

                //printf("%.24f %.24f %.24f",kdisc[index_j],Tbin[index_j],P12_0_b2[index_j]); //GC!
                //printf("%.16e %.16e %.16e",kdisc[index_j],Tbin[index_j],P12_0_b2[index_j]); //GC!
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- start

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_0_b2_oneline_complex_ortho, x2_transfer, &inc, &beta, y_transfer, &inc);
                f12_0_b2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
                P12_0_b2_ortho[index_j] = Tbin[index_j] *                                     //[FACTOR] *
                                          creal(cpow(kdisc[index_j], 3) * f12_0_b2[index_j]); //GC: recheck that indeed it is correct that I multiply by "Tbin[index_j] * [FACTOR]" and a k^3 only. 99.99999% sure...

                //GC: ORTHOGONAL -- finish

                //printf("%.24e %.24e %.24e",kdisc[index_j],Tbin[index_j],P12_0_b2_ortho[index_j]); //GC!
                //printf("\n"); //GC!
            }
                
            }

            //GC: switch to G2, again monopole... Start from vd piece...

            double complex *f12_0_b1bG2;
            double *P12_0_b1bG2;
            class_alloc(f12_0_b1bG2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P12_0_b1bG2, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- start

            double *P12_0_b1bG2_ortho;
            class_alloc(P12_0_b1bG2_ortho, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- finish

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;
            
            if (SWITCH_index == 1) {

            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam_transfer[index_i];
                    nu2 = -0.5 * etam_transfer[index_l];
                    //pnlpt->M12_0_b1bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (1./(csin(nu1*M_PI))/(csin(nu2*M_PI))/(csin((nu1 + nu2)*M_PI))*((15. + 16.*nu2 + 128.*nu1*nu1*nu1*(-1. + nu2)*nu2 - 136.*nu2*nu2 + 128.*nu2*nu2*nu2 - 32.*nu2*nu2*nu2*nu2 + 8.*nu1*nu1*(3. + 46.*nu2 - 84.*nu2*nu2 + 32.*nu2*nu2*nu2) + 8.*nu1*(-6. - 29.*nu2 + 94.*nu2*nu2 - 72.*nu2*nu2*nu2 + 16.*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (32.*nu1*nu1*nu1*nu1*(-1. + 4.*nu2) + 64.*nu1*nu1*nu1*(2. - 9.*nu2 + 4.*nu2*nu2) + 3.*(5. - 16.*nu2 + 8.*nu2*nu2) + 8.*nu1*nu1*(-17. + 94.*nu2 - 84.*nu2*nu2 + 16.*nu2*nu2*nu2) - 8.*nu1*(-2. + 29.*nu2 - 46.*nu2*nu2 + 16.*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (128.*nu1*nu1*nu1*(-1. + nu2)*nu2 - 3.*(5. - 16.*nu2 + 8.*nu2*nu2) + 8.*nu1*nu1*(-3. + 50.*nu2 - 60.*nu2*nu2 + 16.*nu2*nu2*nu2) - 8.*nu1*(-6. + 43.*nu2 - 50.*nu2*nu2 + 16.*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))/(2.*(nu1)*(-1. + 4.*nu1*nu1)*(nu2)*(-5. + 2.*nu1 + 2.*nu2)*(-1. + 4.*nu2*nu2));

                    pnlpt->M12_0_b1bG2_oneline_complex[count2] = 3. * (pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1[count2]);

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_0_b1bG2_oneline_complex_ortho[count2] = 3. * (pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1_ortho[count2]);

                    //GC: ORTHOGONAL -- finish

                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2_transfer[count] = cmsym_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_0_b1bG2_oneline_complex, x2_transfer, &inc, &beta, y_transfer, &inc);
                f12_0_b1bG2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
                P12_0_b1bG2[index_j] = Tbin[index_j] * //[FACTOR] *
                                       creal(cpow(kdisc[index_j], 3) * f12_0_b1bG2[index_j]);

                //printf("%.16e %.16e %.16e",kdisc[index_j],Tbin[index_j],P12_0_b1bG2[index_j]); //GC!
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- start

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_0_b1bG2_oneline_complex_ortho, x2_transfer, &inc, &beta, y_transfer, &inc);
                f12_0_b1bG2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
                P12_0_b1bG2_ortho[index_j] = Tbin[index_j] * //[FACTOR] *
                                             creal(cpow(kdisc[index_j], 3) * f12_0_b1bG2[index_j]);

                //GC: ORTHOGONAL -- finish

                //printf("%.24e %.24e %.24e",kdisc[index_j],Tbin[index_j],P12_0_b1bG2_ortho[index_j]); //GC!
                //printf("\n"); //GC!
            }
                
            }

            //GC: vv piece...

            double complex *f12_0_bG2;
            double *P12_0_bG2;
            class_alloc(f12_0_bG2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P12_0_bG2, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- start

            double *P12_0_bG2_ortho;
            class_alloc(P12_0_bG2_ortho, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- finish

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;
            
            if (SWITCH_index == 1) {

            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam_transfer[index_i];
                    nu2 = -0.5 * etam_transfer[index_l];
                    //pnlpt->M12_0_bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (f/(csin(nu1*M_PI))/(csin(nu2*M_PI))/(csin((nu1 + nu2)*M_PI))*((15. + 16.*nu2 + 128.*nu1*nu1*nu1*(-1. + nu2)*nu2 - 136.*nu2*nu2 + 128.*nu2*nu2*nu2 - 32.*nu2*nu2*nu2*nu2 + 8.*nu1*nu1*(3. + 46.*nu2 - 84.*nu2*nu2 + 32.*nu2*nu2*nu2) + 8.*nu1*(-6. - 29.*nu2 + 94.*nu2*nu2 - 72.*nu2*nu2*nu2 + 16.*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (32.*nu1*nu1*nu1*nu1*(-1. + 4.*nu2) + 64.*nu1*nu1*nu1*(2. - 9.*nu2 + 4.*nu2*nu2) + 3.*(5. - 16.*nu2 + 8.*nu2*nu2) + 8.*nu1*nu1*(-17. + 94.*nu2 - 84.*nu2*nu2 + 16.*nu2*nu2*nu2) - 8.*nu1*(-2. + 29.*nu2 - 46.*nu2*nu2 + 16.*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (128.*nu1*nu1*nu1*(-1. + nu2)*nu2 - 3.*(5. - 16.*nu2 + 8.*nu2*nu2) + 8.*nu1*nu1*(-3. + 50.*nu2 - 60.*nu2*nu2 + 16.*nu2*nu2*nu2) - 8.*nu1*(-6. + 43.*nu2 - 50.*nu2*nu2 + 16.*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))/(6.*(nu1)*(-1. + 4.*nu1*nu1)*(nu2)*(-5. + 2.*nu1 + 2.*nu2)*(-1. + 4.*nu2*nu2));

                    pnlpt->M12_0_bG2_oneline_complex[count2] = f * (pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1[count2]);

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_0_bG2_oneline_complex_ortho[count2] = f * (pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1_ortho[count2]);

                    //GC: ORTHOGONAL -- finish

                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2_transfer[count] = cmsym_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_0_bG2_oneline_complex, x2_transfer, &inc, &beta, y_transfer, &inc);
                f12_0_bG2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
                P12_0_bG2[index_j] = Tbin[index_j] *                                      //[FACTOR] *
                                     creal(cpow(kdisc[index_j], 3) * f12_0_bG2[index_j]); //GC: notice that Misha here had P_0_bG2, so he was not using "f"... Indeed here f seems superfluous...

                //printf("%.16e %.16e %.16e",kdisc[index_j],Tbin[index_j],P12_0_bG2[index_j]); //GC!
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- start

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_0_bG2_oneline_complex_ortho, x2_transfer, &inc, &beta, y_transfer, &inc);
                f12_0_bG2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
                P12_0_bG2_ortho[index_j] = Tbin[index_j] *                                      //[FACTOR] *
                                           creal(cpow(kdisc[index_j], 3) * f12_0_bG2[index_j]); //GC: notice that Misha here had P_0_bG2, so he was not using "f"... Indeed here f seems superfluous...

                //GC: ORTHOGONAL -- finish

                //printf("%.24e %.24e %.24e",kdisc[index_j],Tbin[index_j],P12_0_bG2_ortho[index_j]); //GC!
                //printf("\n"); //GC!
            }
                
            }

            /* Computing Pb1b2 correction in RSD */

            //GC: here Misha means the MONOPOLE...

            double complex *f22_0_b1b2;
            double *P_0_b1b2;
            class_alloc(f22_0_b1b2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P_0_b1b2, Nmax * sizeof(double), pnlpt->error_message);

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;

            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam2[index_i];
                    nu2 = -0.5 * etam2[index_l];

                    pnlpt->M22_0_b1b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3. + 2. * nu1 + 2. * nu2) * (-12. + 7. * (3. + f) * nu1 + 7. * (3. + f) * nu2) / (42. * nu1 * nu2);

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

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_0_b1b2_oneline_complex, x2, &inc, &beta, y, &inc);
                f22_0_b1b2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
                P_0_b1b2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_0_b1b2[index_j]);
            }

            /* Computing b2 correction in RSD */

            double complex *f22_0_b2;
            double *P_0_b2;
            class_alloc(f22_0_b2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P_0_b2, Nmax * sizeof(double), pnlpt->error_message);

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;

            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam2[index_i];
                    nu2 = -0.5 * etam2[index_l];
                    pnlpt->M22_0_b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (7. * f * f * (12. + 6. * nu1 * nu1 - 17. * nu2 + 6. * nu2 * nu2 + nu1 * (-17. + 12. * nu2)) + 5. * f * (24. + 14. * nu1 * nu1 - 37. * nu2 + 14. * nu2 * nu2 + nu1 * (-37. + 28. * nu2))) / (210. * nu1 * nu2);
                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_0_b2_oneline_complex, x2, &inc, &beta, y, &inc);
                f22_0_b2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
                P_0_b2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_0_b2[index_j]);
            }

            /* Computing b1bG2 correction in RSD */

            double complex *f22_0_b1bG2;
            double *P_0_b1bG2;
            class_alloc(f22_0_b1bG2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P_0_b1bG2, Nmax * sizeof(double), pnlpt->error_message);

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;

            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam2[index_i];
                    nu2 = -0.5 * etam2[index_l];
                    pnlpt->M22_0_b1bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3. + 2. * nu1 + 2. * nu2) * (-1. + 2. * nu1 + 2. * nu2) * (7. * f * (2. + nu1 + nu2) + 3. * (6. + 7. * nu1 + 7. * nu2)) / (42. * nu1 * (1. + nu1) * nu2 * (1. + nu2));
                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_0_b1bG2_oneline_complex, x2, &inc, &beta, y, &inc);
                f22_0_b1bG2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
                P_0_b1bG2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_0_b1bG2[index_j]);
            }

            /* Computing bG2 correction in RSD */

            double complex *f22_0_bG2;
            double *P_0_bG2;
            class_alloc(f22_0_bG2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P_0_bG2, Nmax * sizeof(double), pnlpt->error_message);

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;

            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam2[index_i];
                    nu2 = -0.5 * etam2[index_l];
                    pnlpt->M22_0_bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3. + 2. * nu1 + 2. * nu2) * (-1. + 2. * nu1 + 2. * nu2) * (-10. * f + 7. * f * (5. * (nu1 + nu2) + f * (-2. + 3. * nu1 + 3. * nu2))) / (210. * nu1 * (1. + nu1) * nu2 * (1. + nu2));
                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_0_bG2_oneline_complex, x2, &inc, &beta, y, &inc);
                P_0_bG2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
                P_0_bG2[index_j] = creal(cpow(kdisc[index_j], 3) * P_0_bG2[index_j]);
            }

            //GC!
            /* GOING TO THE QUADRUPOLE!!! */
            //GC!

            //GC: start from b2, vd... There is nothing!!!

            //GC: go to b2, vv, quadrupole...

            double complex *f12_2_b2;
            double *P12_2_b2;
            class_alloc(f12_2_b2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P12_2_b2, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- start

            double *P12_2_b2_ortho;
            class_alloc(P12_2_b2_ortho, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- finish

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;
            
            
            if (SWITCH_index == 1) {

            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam2_transfer[index_i];
                    nu2 = -0.5 * etam2_transfer[index_l];
                    //pnlpt->M12_2_b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-2.*f/(csin(nu1*M_PI))/(csin(nu2*M_PI))/(csin((nu1 + nu2)*M_PI))*((15. + 8.*nu1*(-1. + nu2) - 24.*nu2 + 8.*nu2*nu2)*csin(2.*nu1*M_PI) + (15. + 8.*nu1*nu1 + 8.*nu1*(-3. + nu2) - 8.*nu2)*csin(2.*nu2*M_PI) + (9. + 8.*nu1*(-1. + nu2) - 8.*nu2)*csin(2.*(nu1 + nu2)*M_PI)))/(3.*(-1. + 2.*nu1)*(-1. + 2.*nu2)*(-5. + 2.*nu1 + 2.*nu2));

                    pnlpt->M12_2_b2_oneline_complex[count2] = 2. * f * (pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1[count2]);

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_2_b2_oneline_complex_ortho[count2] = 2. * f * (pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1_ortho[count2]);

                    //GC: ORTHOGONAL -- finish

                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2_transfer[count] = cmsym2_transfer[count] * cpow(kdisc[index_j], etam2_transfer[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_2_b2_oneline_complex, x2_transfer, &inc, &beta, y_transfer, &inc);
                f12_2_b2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
                P12_2_b2[index_j] = Tbin[index_j] * //[FACTOR] *
                                    creal(cpow(kdisc[index_j], 3) * f12_2_b2[index_j]);

                //printf("%.16e %.16e %.16e",kdisc[index_j],Tbin[index_j],P12_2_b2[index_j]); //GC!
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- start

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_2_b2_oneline_complex_ortho, x2_transfer, &inc, &beta, y_transfer, &inc);
                f12_2_b2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
                P12_2_b2_ortho[index_j] = Tbin[index_j] * //[FACTOR] *
                                          creal(cpow(kdisc[index_j], 3) * f12_2_b2[index_j]);

                //GC: ORTHOGONAL -- finish

                //printf("%.24e %.24e %.24e",kdisc[index_j],Tbin[index_j],P12_2_b2_ortho[index_j]); //GC!
                //printf("\n"); //GC!
            }
                
            }

            //GC: go to bG2. Also here, only vv...

            double complex *f12_2_bG2; //GC: again it was unused for some reason...
            double *P12_2_bG2;
            class_alloc(f12_2_bG2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P12_2_bG2, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- start

            double *P12_2_bG2_ortho;
            class_alloc(P12_2_bG2_ortho, Nmax * sizeof(double), pnlpt->error_message);

            //GC: ORTHOGONAL -- finish

            count2 = 0;
            
            if (SWITCH_index == 1) {
            
            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam_transfer[index_i];
                    nu2 = -0.5 * etam_transfer[index_l];
                    //pnlpt->M12_2_bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (f/(csin(nu1*M_PI))/(csin(nu2*M_PI))/(csin((nu1 + nu2)*M_PI))*((15. + 16.*nu2 + 128.*nu1*nu1*nu1*(-1. + nu2)*nu2 - 136.*nu2*nu2 + 128.*nu2*nu2*nu2 - 32.*nu2*nu2*nu2*nu2 + 8.*nu1*nu1*(3. + 46.*nu2 - 84.*nu2*nu2 + 32.*nu2*nu2*nu2) + 8.*nu1*(-6. - 29.*nu2 + 94.*nu2*nu2 - 72.*nu2*nu2*nu2 + 16.*nu2*nu2*nu2*nu2))*csin(2.*nu1*M_PI) + (32.*nu1*nu1*nu1*nu1*(-1. + 4.*nu2) + 64.*nu1*nu1*nu1*(2. - 9.*nu2 + 4.*nu2*nu2) + 3.*(5. - 16.*nu2 + 8.*nu2*nu2) + 8.*nu1*nu1*(-17. + 94.*nu2 - 84.*nu2*nu2 + 16.*nu2*nu2*nu2) - 8.*nu1*(-2. + 29.*nu2 - 46.*nu2*nu2 + 16.*nu2*nu2*nu2))*csin(2.*nu2*M_PI) + (128.*nu1*nu1*nu1*(-1. + nu2)*nu2 - 3.*(5. - 16.*nu2 + 8.*nu2*nu2) + 8.*nu1*nu1*(-3. + 50.*nu2 - 60.*nu2*nu2 + 16.*nu2*nu2*nu2) - 8.*nu1*(-6. + 43.*nu2 - 50.*nu2*nu2 + 16.*nu2*nu2*nu2))*csin(2.*(nu1 + nu2)*M_PI)))/(3.*(nu1)*(-1. + 4.*nu1*nu1)*(nu2)*(-5. + 2.*nu1 + 2.*nu2)*(-1. + 4.*nu2*nu2));

                    pnlpt->M12_2_bG2_oneline_complex[count2] = 2. * f * (pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1[count2]);

                    //GC: ORTHOGONAL -- start

                    pnlpt->M12_2_bG2_oneline_complex_ortho[count2] = 2. * f * (pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1_ortho[count2]);

                    //GC: ORTHOGONAL -- finish

                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2_transfer[count] = cmsym_transfer[count] * cpow(kdisc[index_j], etam_transfer[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_2_bG2_oneline_complex, x2_transfer, &inc, &beta, y_transfer, &inc);
                f12_2_bG2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
                P12_2_bG2[index_j] = Tbin[index_j] * //[FACTOR] *
                                     creal(cpow(kdisc[index_j], 3) * f12_2_bG2[index_j]);

                //printf("%.16e %.16e %.16e",kdisc[index_j],Tbin[index_j],P12_2_bG2[index_j]); //GC!
                //printf("\n"); //GC!

                //GC: ORTHOGONAL -- start

                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M12_2_bG2_oneline_complex_ortho, x2_transfer, &inc, &beta, y_transfer, &inc);
                f12_2_bG2[index_j] = zdotu_(&Nmaxf, x2_transfer, &inc, y_transfer, &inc);
                P12_2_bG2_ortho[index_j] = Tbin[index_j] * //[FACTOR] *
                                           creal(cpow(kdisc[index_j], 3) * f12_2_bG2[index_j]);

                //GC: ORTHOGONAL -- finish

                //printf("%.24e %.24e %.24e",kdisc[index_j],Tbin[index_j],P12_2_bG2_ortho[index_j]); //GC!
                //printf("\n"); //GC!
            }
                
            }

            /* Computing Pb1b2 correction for the Quadrupole */

            double complex *f22_2_b1b2;
            double *P_2_b1b2;
            class_alloc(f22_2_b1b2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P_2_b1b2, Nmax * sizeof(double), pnlpt->error_message);

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;

            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam2[index_i];
                    nu2 = -0.5 * etam2[index_l];
                    pnlpt->M22_2_b1b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3. + 2. * nu1 + 2. * nu2) * f * (nu1 + nu2) / (3. * nu1 * nu2);
                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_2_b1b2_oneline_complex, x2, &inc, &beta, y, &inc);
                f22_2_b1b2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
                P_2_b1b2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_2_b1b2[index_j]);
            }

            /* Computing Pb2 correction for the Quadrupole */

            double complex *f22_2_b2;
            double *P_2_b2;
            class_alloc(f22_2_b2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P_2_b2, Nmax * sizeof(double), pnlpt->error_message);

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;

            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam2[index_i];
                    nu2 = -0.5 * etam2[index_l];
                    pnlpt->M22_2_b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3. + 2. * nu1 + 2. * nu2) * f * (-16. + 14. * (nu1 + nu2) + f * (-13. + 12. * (nu1 + nu2))) / (42. * nu1 * nu2);
                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_2_b2_oneline_complex, x2, &inc, &beta, y, &inc);
                f22_2_b2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
                P_2_b2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_2_b2[index_j]);
            }

            /* Computing b1bG2 correction for Quadrupole */

            double complex *f22_2_b1bG2;
            double *P_2_b1bG2;
            class_alloc(f22_2_b1bG2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P_2_b1bG2, Nmax * sizeof(double), pnlpt->error_message);

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;

            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam2[index_i];
                    nu2 = -0.5 * etam2[index_l];
                    pnlpt->M22_2_b1bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3. + 2. * nu1 + 2. * nu2) * (-1. + 2. * nu1 + 2. * nu2) * f * (2. + nu1 + nu2) / (3. * nu1 * (1. + nu1) * nu2 * (1. + nu2));
                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_2_b1bG2_oneline_complex, x2, &inc, &beta, y, &inc);
                f22_2_b1bG2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
                P_2_b1bG2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_2_b1bG2[index_j]);
            }

            /* Computing bG2 correction for Quadrupole */

            double complex *f22_2_bG2;
            double *P_2_bG2;
            class_alloc(f22_2_bG2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P_2_bG2, Nmax * sizeof(double), pnlpt->error_message);

            count2 = 0;
            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam2[index_i];
                    nu2 = -0.5 * etam2[index_l];
                    pnlpt->M22_2_bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3. + 2. * nu1 + 2. * nu2) * (-1. + 2. * nu1 + 2. * nu2) * f * (-2. - f + 7. * (nu1 + nu2) + 6. * f * (nu1 + nu2)) / (21. * nu1 * (1. + nu1) * nu2 * (1. + nu2));
                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_2_bG2_oneline_complex, x2, &inc, &beta, y, &inc);
                P_2_bG2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
                P_2_bG2[index_j] = creal(cpow(kdisc[index_j], 3) * P_2_bG2[index_j]);
            }

            /* Computing Pb2 correction for Hexadecapole */

            double complex *f22_4_b2;
            double *P_4_b2;
            class_alloc(f22_4_b2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P_4_b2, Nmax * sizeof(double), pnlpt->error_message);

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;

            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam2[index_i];
                    nu2 = -0.5 * etam2[index_l];
                    pnlpt->M22_4_b2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3. + 2. * nu1 + 2. * nu2) * (-1. + 2. * nu1 + 2. * nu2) * 2. * f * f / (35. * nu1 * nu2);
                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
                }
                zspmv_(&uplo, &Nmaxf, &alpha, pnlpt->M22_4_b2_oneline_complex, x2, &inc, &beta, y, &inc);
                f22_4_b2[index_j] = zdotu_(&Nmaxf, x2, &inc, y, &inc);
                P_4_b2[index_j] = creal(cpow(kdisc[index_j], 3) * f22_4_b2[index_j]);
            }

            /* Computing bG2 correction for Hexadecapole */

            double complex *f22_4_bG2;
            double *P_4_bG2;
            class_alloc(f22_4_bG2, Nmax * sizeof(complex double), pnlpt->error_message);
            class_alloc(P_4_bG2, Nmax * sizeof(double), pnlpt->error_message);

            count2 = 0;
            index_i = 0;
            index_l = 0;
            index_j = 0;

            for (index_l = 0; index_l < Nmax + 1; index_l++)
            {
                for (index_i = index_l; index_i < Nmax + 1; index_i++)
                {
                    nu1 = -0.5 * etam2[index_i];
                    nu2 = -0.5 * etam2[index_l];
                    pnlpt->M22_4_bG2_oneline_complex[count2] = pnlpt->M22basic_oneline_complex[count2] * (-3. + 2. * nu1 + 2. * nu2) * (-1. + 2. * nu1 + 2. * nu2) * 4. * f * f * (1. + nu1 + nu2) / (35. * nu1 * (1. + nu1) * nu2 * (1. + nu2));
                    count2++;
                }
            }

            for (index_j = 0; index_j < Nmax; index_j++)
            {
                for (count = 0; count < Nmax + 1; count++)
                {
                    x2[count] = cmsym2[count] * cpow(kdisc[index_j], etam2[count]);
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
            class_alloc(P_Id2d2_2, Nmax * sizeof(double), pnlpt->error_message);
            double *P_Id2d2_4;
            class_alloc(P_Id2d2_4, Nmax * sizeof(double), pnlpt->error_message);
            double *P_Id2G2_2;
            class_alloc(P_Id2G2_2, Nmax * sizeof(double), pnlpt->error_message);
            double *P_Id2G2_4;
            class_alloc(P_Id2G2_4, Nmax * sizeof(double), pnlpt->error_message);
            double *P_IG2G2_2;
            class_alloc(P_IG2G2_2, Nmax * sizeof(double), pnlpt->error_message);
            double *P_IG2G2_4;
            class_alloc(P_IG2G2_4, Nmax * sizeof(double), pnlpt->error_message);
            double *P_4_b1b2;
            class_alloc(P_4_b1b2, Nmax * sizeof(double), pnlpt->error_message);
            double *P_4_b1bG2;
            class_alloc(P_4_b1bG2, Nmax * sizeof(double), pnlpt->error_message);

            double *P_IFG2_0b1_x;
            class_alloc(P_IFG2_0b1_x, Nmax * sizeof(double), pnlpt->error_message);
            double *P_IFG2_0;
            class_alloc(P_IFG2_0, Nmax * sizeof(double), pnlpt->error_message);
            double *P_IFG2_2;
            class_alloc(P_IFG2_2, Nmax * sizeof(double), pnlpt->error_message);

            double *P_Id2d2_new;
            class_alloc(P_Id2d2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P_Id2G2_new;
            class_alloc(P_Id2G2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P_IG2G2_new;
            class_alloc(P_IG2G2_new, sizeof(double) * Nmax, pnlpt->error_message);

            double *P_0_bG2_new;
            class_alloc(P_0_bG2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P_2_bG2_new;
            class_alloc(P_2_bG2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P_4_bG2_new;
            class_alloc(P_4_bG2_new, sizeof(double) * Nmax, pnlpt->error_message);

            double *P_IFG2_new;
            class_alloc(P_IFG2_new, sizeof(double) * Nmax, pnlpt->error_message);

            double *P_0_b1b2_new;
            class_alloc(P_0_b1b2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P_2_b1b2_new;
            class_alloc(P_2_b1b2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P_0_b1bG2_new;
            class_alloc(P_0_b1bG2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P_2_b1bG2_new;
            class_alloc(P_2_b1bG2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P_0_b2_new;
            class_alloc(P_0_b2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P_2_b2_new;
            class_alloc(P_2_b2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P_4_b2_new;
            class_alloc(P_4_b2_new, sizeof(double) * Nmax, pnlpt->error_message);

            //GC: I need to look at these? I think so...

            //GC: why does he use also Id2d2??? No, this is a name given for something that will be used later? It seems that Id2d2, Id2G2, and IG2G2 are for the pure bias 22 terms... These do NOT come with any power of \mu^2 also in redshift space, so he just takes them from the real space calculation -> they do not have any multiplication by Z1. Things like P_Id2d2_2 are defined first above here. Also things like P_4_b1b2. But not P_0_b1b2 or P_0_b1b2, which enters already before. What does this mean? Oh, it means that such term can only arise after the AP!!! Cool... Notice that in my case I never have pure bias terms. I always multiply by Z1. So I will just define the following objects -> the AP generates all the multipoles, so we will have 0, 2 and 4. But we need to be careful about the powers of b1.

            double *P12_0_b1b2_new;
            class_alloc(P12_0_b1b2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_0_b1bG2_new;
            class_alloc(P12_0_b1bG2_new, sizeof(double) * Nmax, pnlpt->error_message);

            double *P12_0_b2_new;
            class_alloc(P12_0_b2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_2_b2_new;
            class_alloc(P12_2_b2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_0_bG2_new;
            class_alloc(P12_0_bG2_new, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_2_bG2_new;
            class_alloc(P12_2_bG2_new, sizeof(double) * Nmax, pnlpt->error_message);

            //GC: ORTHOGONAL -- start

            double *P12_0_b1b2_new_ortho;
            class_alloc(P12_0_b1b2_new_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_0_b1bG2_new_ortho;
            class_alloc(P12_0_b1bG2_new_ortho, sizeof(double) * Nmax, pnlpt->error_message);

            double *P12_0_b2_new_ortho;
            class_alloc(P12_0_b2_new_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_2_b2_new_ortho;
            class_alloc(P12_2_b2_new_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_0_bG2_new_ortho;
            class_alloc(P12_0_bG2_new_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_2_bG2_new_ortho;
            class_alloc(P12_2_bG2_new_ortho, sizeof(double) * Nmax, pnlpt->error_message);

            //GC: ORTHOGONAL -- finish

            //GC: first, these are what I already had... 0 -> P12_0_b2, P12_0_b1b2, P12_0_bG2, P12_0_b1bG2; 2 -> P12_2_b2, P12_2_bG2...

            //GC: then, we need all these up to 4, via AP, with NO CHANGE IN b1 POWERS...

            double *P12_2_b1b2;
            class_alloc(P12_2_b1b2, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_2_b1bG2;
            class_alloc(P12_2_b1bG2, sizeof(double) * Nmax, pnlpt->error_message);

            double *P12_4_b1b2;
            class_alloc(P12_4_b1b2, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_4_b1bG2;
            class_alloc(P12_4_b1bG2, sizeof(double) * Nmax, pnlpt->error_message);

            double *P12_4_b2;
            class_alloc(P12_4_b2, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_4_bG2;
            class_alloc(P12_4_bG2, sizeof(double) * Nmax, pnlpt->error_message);

            //GC: ORTHOGONAL -- start

            double *P12_2_b1b2_ortho;
            class_alloc(P12_2_b1b2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_2_b1bG2_ortho;
            class_alloc(P12_2_b1bG2_ortho, sizeof(double) * Nmax, pnlpt->error_message);

            double *P12_4_b1b2_ortho;
            class_alloc(P12_4_b1b2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_4_b1bG2_ortho;
            class_alloc(P12_4_b1bG2_ortho, sizeof(double) * Nmax, pnlpt->error_message);

            double *P12_4_b2_ortho;
            class_alloc(P12_4_b2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            double *P12_4_bG2_ortho;
            class_alloc(P12_4_bG2_ortho, sizeof(double) * Nmax, pnlpt->error_message);

            //GC: ORTHOGONAL -- finish

            for (index_j = 0; index_j < Nmax; index_j++)
            {
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

                //GC...

                P12_0_b1b2_new[index_j] = P12_0_b1b2[index_j];
                P12_0_b1bG2_new[index_j] = P12_0_b1bG2[index_j];
                P12_0_b2_new[index_j] = P12_0_b2[index_j];
                P12_2_b2_new[index_j] = P12_2_b2[index_j];
                P12_0_bG2_new[index_j] = P12_0_bG2[index_j];
                P12_2_bG2_new[index_j] = P12_2_bG2[index_j];

                //GC: ORTHOGONAL -- start

                P12_0_b1b2_new_ortho[index_j] = P12_0_b1b2_ortho[index_j];
                P12_0_b1bG2_new_ortho[index_j] = P12_0_b1bG2_ortho[index_j];
                P12_0_b2_new_ortho[index_j] = P12_0_b2_ortho[index_j];
                P12_2_b2_new_ortho[index_j] = P12_2_b2_ortho[index_j];
                P12_0_bG2_new_ortho[index_j] = P12_0_bG2_ortho[index_j];
                P12_2_bG2_new_ortho[index_j] = P12_2_bG2_ortho[index_j];

                //GC: ORTHOGONAL -- finish
            }

            //GC: this seems to save the multipoles of before? Unclear: it does not use the "_new" variables? Some do, for some reason... Those above are multipoles, right?

            //GC: things to ask. Why in the lines directly below someone uses the "_new" variables and the others do not? Then, here he initializes the "_ap_out" variables to zero... The structure is very similar to what happens around 4530... But if I took this statement at face value, the things like P_4_bG2, for example, would NOT be hexadecapole contributions, but \mu^4 contributions? No, what he does is the following -> he simply takes MONOPOLE * P_0 + QUADRUPOLE * P_2 and then gets the full powers!!! OF COURSE!!! Ok. Now I just want to remember what happens when I do the IR resummation -> what happens is that Tbin and Pbin are the ones that, if you toggle in the IR resummation, have the wiggly part suppressed by the exponential. For matter he first does the calculation with Pbin both in the loop and outside, which is what I also did there, and now for biased tracers. Then, the only tricky part is that at some point he does the IR resummation also in redshift space, and by this he means that the IR resummation outside is done via a \mu-dependent damping factor. Only outside, though: before, the multipoles were computed with the IR-resummed power spectrum in real space everywhere (which was Pbin, or Tbin for me). Now, the difference is that only the outside piece is computed with a \mu-dependent damping... But let's be careful -> before, when just computing the multipoles, [everything] is computed in terms of Pbin. Which is damped. So, P13, P22 and P12 have all a damped PS. They are a functional of a linear damped PS. No, wait -> maybe the point is that there is a "IF" flag, that there is set to NO resummation, so Tbin and Pbin are untouched... Is that it? But then Tbin for me is essentially never used??? This I must still understand... The "IF" chains must be understood... Let's try to understand them AFTER I coded up everything... SEE .tex file on this!!!

            double *dd_P_Id2d2;
            class_alloc(dd_P_Id2d2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_Id2d2, 1, dd_P_Id2d2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_Id2d2_ap_out = 0;

            double *dd_P_Id2G2;
            class_alloc(dd_P_Id2G2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_Id2G2, 1, dd_P_Id2G2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_Id2G2_ap_out = 0;

            double *dd_P_IG2G2;
            class_alloc(dd_P_IG2G2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_IG2G2, 1, dd_P_IG2G2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_IG2G2_ap_out = 0;

            double *dd_P_0_b1b2;
            class_alloc(dd_P_0_b1b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_0_b1b2_new, 1, dd_P_0_b1b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_0_b1b2_ap_out = 0;

            double *dd_P_2_b1b2;
            class_alloc(dd_P_2_b1b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_2_b1b2_new, 1, dd_P_2_b1b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_2_b1b2_ap_out = 0;

            double *dd_P_0_b1bG2;
            class_alloc(dd_P_0_b1bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_0_b1bG2, 1, dd_P_0_b1bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_0_b1bG2_ap_out = 0;

            double *dd_P_2_b1bG2;
            class_alloc(dd_P_2_b1bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_2_b1bG2, 1, dd_P_2_b1bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_2_b1bG2_ap_out = 0;

            double *dd_P_0_b2;
            class_alloc(dd_P_0_b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_0_b2_new, 1, dd_P_0_b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_0_b2_ap_out = 0;
            double *dd_P_2_b2;
            class_alloc(dd_P_2_b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_2_b2_new, 1, dd_P_2_b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_2_b2_ap_out = 0;
            double *dd_P_4_b2;
            class_alloc(dd_P_4_b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_4_b2_new, 1, dd_P_4_b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_4_b2_ap_out = 0;

            double *dd_P_0_bG2;
            class_alloc(dd_P_0_bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_0_bG2_new, 1, dd_P_0_bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_0_bG2_ap_out = 0;
            double *dd_P_2_bG2;
            class_alloc(dd_P_2_bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_2_bG2_new, 1, dd_P_2_bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_2_bG2_ap_out = 0;
            double *dd_P_4_bG2;
            class_alloc(dd_P_4_bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_4_bG2_new, 1, dd_P_4_bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_4_bG2_ap_out = 0;

            double *dd_P_IFG2;
            class_alloc(dd_P_IFG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_IFG2_new, 1, dd_P_IFG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P_IFG2_ap_out = 0;

            double *dd_Pbin;
            class_alloc(dd_Pbin, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, Pbin, 1, dd_Pbin, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double Pbin_ap_out = 0;

            //GC -> here I interpolate the ones I already computed (indeed see that you call the "new" ones in the call)... I also get Tbin...

            //GC -> what I already had were: 0 -> P12_0_b2, P12_0_b1b2, P12_0_bG2, P12_0_b1bG2; 2 -> P12_2_b2, P12_2_bG2...

            double *dd_P12_0_b2;
            class_alloc(dd_P12_0_b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_b2_new, 1, dd_P12_0_b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_0_b2_ap_out = 0;

            //GC: ORTHOGONAL -- start

            double *dd_P12_0_b2_ortho;
            class_alloc(dd_P12_0_b2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_b2_new_ortho, 1, dd_P12_0_b2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_0_b2_ap_out_ortho = 0;

            //GC: ORTHOGONAL -- finish

            //GC -> the structure is that the "new" are interpolated, and then called at the AP values of momenta -> the output of this call is stored in the "out" variables...

            double *dd_P12_0_b1b2;
            class_alloc(dd_P12_0_b1b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_b1b2_new, 1, dd_P12_0_b1b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_0_b1b2_ap_out = 0;

            //GC: ORTHOGONAL -- start

            double *dd_P12_0_b1b2_ortho;
            class_alloc(dd_P12_0_b1b2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_b1b2_new_ortho, 1, dd_P12_0_b1b2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_0_b1b2_ap_out_ortho = 0;

            //GC: ORTHOGONAL -- finish

            double *dd_P12_0_bG2;
            class_alloc(dd_P12_0_bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_bG2_new, 1, dd_P12_0_bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_0_bG2_ap_out = 0;

            //GC: ORTHOGONAL -- start

            double *dd_P12_0_bG2_ortho;
            class_alloc(dd_P12_0_bG2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_bG2_new_ortho, 1, dd_P12_0_bG2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_0_bG2_ap_out_ortho = 0;

            //GC: ORTHOGONAL -- finish

            double *dd_P12_0_b1bG2;
            class_alloc(dd_P12_0_b1bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_b1bG2_new, 1, dd_P12_0_b1bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_0_b1bG2_ap_out = 0;

            //GC: ORTHOGONAL -- start

            double *dd_P12_0_b1bG2_ortho;
            class_alloc(dd_P12_0_b1bG2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_b1bG2_new_ortho, 1, dd_P12_0_b1bG2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_0_b1bG2_ap_out_ortho = 0;

            //GC: ORTHOGONAL -- finish

            double *dd_P12_2_b2;
            class_alloc(dd_P12_2_b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_2_b2_new, 1, dd_P12_2_b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_2_b2_ap_out = 0;

            //GC: ORTHOGONAL -- start

            double *dd_P12_2_b2_ortho;
            class_alloc(dd_P12_2_b2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_2_b2_new_ortho, 1, dd_P12_2_b2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_2_b2_ap_out_ortho = 0;

            //GC: ORTHOGONAL -- finish

            double *dd_P12_2_bG2;
            class_alloc(dd_P12_2_bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_2_bG2_new, 1, dd_P12_2_bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_2_bG2_ap_out = 0;

            //GC: ORTHOGONAL -- start

            double *dd_P12_2_bG2_ortho;
            class_alloc(dd_P12_2_bG2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_2_bG2_new_ortho, 1, dd_P12_2_bG2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double P12_2_bG2_ap_out_ortho = 0;

            //GC: ORTHOGONAL -- finish

            double *dd_Tbin;
            class_alloc(dd_Tbin, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, Tbin, 1, dd_Tbin, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            double Tbin_ap_out = 0;

            //GC -> maybe this is unneeded!!!

            double Pd2d2_in = 0.;
            double Pd2G2_in = 0.;
            double PG2G2_in = 0.;
            double Pb1b2_in = 0.;
            double Pb1bG2_in = 0.;
            double Pb2_in = 0.;
            double PbG2_in = 0.;
            double IFG2_in = 0.;

            //GC -> here for me I just need four things...

            double P12b1b2_in = 0.;
            double P12b1bG2_in = 0.;
            double P12b2_in = 0.;
            double P12bG2_in = 0.;

            //GC: ORTHOGONAL -- start

            double P12b1b2_in_ortho = 0.;
            double P12b1bG2_in_ortho = 0.;
            double P12b2_in_ortho = 0.;
            double P12bG2_in_ortho = 0.;

            //GC: ORTHOGONAL -- finish

            double LegendreP0 = 1.;
            double LegendreP2 = 0.;
            double LegendreP4 = 0.;

            double LegendreP2true = 0.;
            double LegendreP4true = 0.;

            int index_gauss2 = 0;
            double mu = 0.;
            double mutrue = 0.;
            double ktrue = 0.;
            last_index = 0;
            double Sigmatot = 0.;
            double p_lo = 0.;
            double Exp = 0.;

            //GC!

            double t_lo = 0.;

            //GC!

            /*
    hratio = 1.;
    Dratio = 1.;
    */

            for (index_j = 0; index_j < Nside; index_j++)
            {
                P_IFG2_0b1_x[index_j] = 0.;
                P_IFG2_0[index_j] = 0.;
                P_IFG2_2[index_j] = 0.;

                P_IFG2_0b1_x[Nmax - 1 - index_j] = 0.;
                P_IFG2_0[Nmax - 1 - index_j] = 0.;
                P_IFG2_2[Nmax - 1 - index_j] = 0.;

                P_Id2d2_2[index_j] = 0.;
                P_Id2d2_4[index_j] = 0.;
                P_Id2d2_2[Nmax - 1 - index_j] = 0.;
                P_Id2d2_4[Nmax - 1 - index_j] = 0.;

                P_Id2G2_2[index_j] = 0.;
                P_Id2G2_4[index_j] = 0.;
                P_Id2G2_2[Nmax - 1 - index_j] = 0.;
                P_Id2G2_4[Nmax - 1 - index_j] = 0.;

                P_IG2G2_2[index_j] = 0.;
                P_IG2G2_4[index_j] = 0.;
                P_IG2G2_2[Nmax - 1 - index_j] = 0.;
                P_IG2G2_4[Nmax - 1 - index_j] = 0.;

                P_4_b1b2[index_j] = 0.;
                P_4_b1bG2[index_j] = 0.;
                P_4_b1b2[Nmax - 1 - index_j] = 0.;
                P_4_b1bG2[Nmax - 1 - index_j] = 0.;

                //GC! It seems I only put the ones I do not already have, the ones that have not been "newed"...

                P12_2_b1b2[index_j] = 0.;
                P12_2_b1b2[Nmax - 1 - index_j] = 0.;

                P12_2_b1bG2[index_j] = 0.;
                P12_2_b1bG2[Nmax - 1 - index_j] = 0.;

                P12_4_b1b2[index_j] = 0.;
                P12_4_b1b2[Nmax - 1 - index_j] = 0.;

                P12_4_b1bG2[index_j] = 0.;
                P12_4_b1bG2[Nmax - 1 - index_j] = 0.;

                P12_4_b2[index_j] = 0.;
                P12_4_b2[Nmax - 1 - index_j] = 0.;

                P12_4_bG2[index_j] = 0.;
                P12_4_bG2[Nmax - 1 - index_j] = 0.;

                //GC: ORTHOGONAL -- start

                P12_2_b1b2_ortho[index_j] = 0.;
                P12_2_b1b2_ortho[Nmax - 1 - index_j] = 0.;

                P12_2_b1bG2_ortho[index_j] = 0.;
                P12_2_b1bG2_ortho[Nmax - 1 - index_j] = 0.;

                P12_4_b1b2_ortho[index_j] = 0.;
                P12_4_b1b2_ortho[Nmax - 1 - index_j] = 0.;

                P12_4_b1bG2_ortho[index_j] = 0.;
                P12_4_b1bG2_ortho[Nmax - 1 - index_j] = 0.;

                P12_4_b2_ortho[index_j] = 0.;
                P12_4_b2_ortho[Nmax - 1 - index_j] = 0.;

                P12_4_bG2_ortho[index_j] = 0.;
                P12_4_bG2_ortho[Nmax - 1 - index_j] = 0.;

                //GC: ORTHOGONAL -- finish
            }

            for (index_j = Nside; index_j < Nmax - Nside; index_j++)
            {

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

                //GC!

                //GC: here instead I have all of them it seems. Which should make sense since now I have the numerical integration... These are the things I sum over in the loop for the integration. The final output!!! Indeed, they are 24 terms, and 24 are the terms below that follow the structure e.g. "P_0_b2[index_j] +=  Pb2_in*LegendreP0/2.;". Notice that P_Id2d2 and P_IG2G2 are the MONOPOLES! Ok -> I don't know what happens for FG2 stuff, with the _x subscript, but ok...

                P12_0_b2[index_j] = 0.;
                P12_2_b2[index_j] = 0.;
                P12_4_b2[index_j] = 0.;

                P12_0_bG2[index_j] = 0.;
                P12_2_bG2[index_j] = 0.;
                P12_4_bG2[index_j] = 0.;

                P12_0_b1b2[index_j] = 0.;
                P12_2_b1b2[index_j] = 0.;
                P12_4_b1b2[index_j] = 0.;

                P12_0_b1bG2[index_j] = 0.;
                P12_2_b1bG2[index_j] = 0.;
                P12_4_b1bG2[index_j] = 0.;

                //GC: ORTHOGONAL -- start

                P12_0_b2_ortho[index_j] = 0.;
                P12_2_b2_ortho[index_j] = 0.;
                P12_4_b2_ortho[index_j] = 0.;

                P12_0_bG2_ortho[index_j] = 0.;
                P12_2_bG2_ortho[index_j] = 0.;
                P12_4_bG2_ortho[index_j] = 0.;

                P12_0_b1b2_ortho[index_j] = 0.;
                P12_2_b1b2_ortho[index_j] = 0.;
                P12_4_b1b2_ortho[index_j] = 0.;

                P12_0_b1bG2_ortho[index_j] = 0.;
                P12_2_b1bG2_ortho[index_j] = 0.;
                P12_4_b1bG2_ortho[index_j] = 0.;

                //GC: ORTHOGONAL -- finish

                Tnw_ap_out = 0.;
                Tw_ap_out = 0.;

                for (index_gauss2 = 0; index_gauss2 < 40; index_gauss2++)
                {

                    /*
            mu = pnlpt->gauss_x[index_gauss2];
            mutrue = mu*hratio/pow((1./Dratio/Dratio + (hratio*hratio - 1./Dratio/Dratio)*mu*mu),0.5);
            ktrue = kdisc[index_j]*pow((1./Dratio/Dratio + (hratio*hratio - 1./Dratio/Dratio)*mu*mu),0.5);


            mutrue = mu;
            ktrue = kdisc[index_j];
             */

                    mu = pnlpt->gauss_x[index_gauss2];

                    if (pnlpt->AP_effect == AP_effect_yes)
                    {
                        mutrue = mu * hratio / pow((1. / Dratio / Dratio + (hratio * hratio - 1. / Dratio / Dratio) * mu * mu), 0.5);
                        ktrue = kdisc[index_j] * pow((1. / Dratio / Dratio + (hratio * hratio - 1. / Dratio / Dratio) * mu * mu), 0.5);
                    }

                    else
                    {
                        mutrue = mu;
                        ktrue = kdisc[index_j];
                    }

                    LegendreP2 = (3. * pow(mu, 2.) - 1.) / 2.;
                    LegendreP4 = (35. * pow(mu, 4.) - 30. * pow(mu, 2.) + 3.) / 8.;
                    LegendreP2true = (3. * pow(mutrue, 2.) - 1.) / 2.;
                    LegendreP4true = (35. * pow(mutrue, 4.) - 30. * pow(mutrue, 2.) + 3.) / 8.;

                    class_call(array_interpolate_spline(kdisc, Nmax, Pnw, dd_Pnw, 1, ktrue, &last_index, &Pnw_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, Pw, dd_Pw, 1, ktrue, &last_index, &Pw_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, Pbin, dd_Pbin, 1, ktrue, &last_index, &Pbin_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    
                    
                    if (SWITCH_index == 1) {

                    //GC...
                    class_call(array_interpolate_spline(kdisc, Nmax, Tnw, dd_Tnw, 1, ktrue, &last_index, &Tnw_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, Tw, dd_Tw, 1, ktrue, &last_index, &Tw_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, Tbin, dd_Tbin, 1, ktrue, &last_index, &Tbin_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    //GC... Calls the transfer functions at the AP values of k...

                    //GC!

                    //GC -> start with evaluating the "new" at the modified k. But now such k depends on \mu. Does it make sense to evaluate the multipoles at a k that depends on \mu? Well, let's see -> I take e.g. P_2(k) * Leg_2(\mu) -> this is a P(k,\mu). Then, I evaluate at the true k on the right. It is then equivalent to evaluating P_2(k) at the true, \mu-dependent k...

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_b1b2_new, dd_P12_0_b1b2, 1, ktrue, &last_index, &P12_0_b1b2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_b1bG2_new, dd_P12_0_b1bG2, 1, ktrue, &last_index, &P12_0_b1bG2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_b2_new, dd_P12_0_b2, 1, ktrue, &last_index, &P12_0_b2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_2_b2_new, dd_P12_2_b2, 1, ktrue, &last_index, &P12_2_b2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_bG2_new, dd_P12_0_bG2, 1, ktrue, &last_index, &P12_0_bG2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_2_bG2_new, dd_P12_2_bG2, 1, ktrue, &last_index, &P12_2_bG2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    //GC: ORTHOGONAL -- start

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_b1b2_new_ortho, dd_P12_0_b1b2_ortho, 1, ktrue, &last_index, &P12_0_b1b2_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_b1bG2_new_ortho, dd_P12_0_b1bG2_ortho, 1, ktrue, &last_index, &P12_0_b1bG2_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_b2_new_ortho, dd_P12_0_b2_ortho, 1, ktrue, &last_index, &P12_0_b2_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_2_b2_new_ortho, dd_P12_2_b2_ortho, 1, ktrue, &last_index, &P12_2_b2_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_bG2_new_ortho, dd_P12_0_bG2_ortho, 1, ktrue, &last_index, &P12_0_bG2_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P12_2_bG2_new_ortho, dd_P12_2_bG2_ortho, 1, ktrue, &last_index, &P12_2_bG2_ap_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    //GC: ORTHOGONAL -- finish
                        
                    }
                    

                    //GC!

                    class_call(array_interpolate_spline(kdisc, Nmax, P_IFG2_new, dd_P_IFG2, 1, ktrue, &last_index, &P_IFG2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P_Id2d2_new, dd_P_Id2d2, 1, ktrue, &last_index, &P_Id2d2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P_Id2G2_new, dd_P_Id2G2, 1, ktrue, &last_index, &P_Id2G2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P_IG2G2_new, dd_P_IG2G2, 1, ktrue, &last_index, &P_IG2G2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P_0_b1b2_new, dd_P_0_b1b2, 1, ktrue, &last_index, &P_0_b1b2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P_2_b1b2_new, dd_P_2_b1b2, 1, ktrue, &last_index, &P_2_b1b2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P_0_b1bG2_new, dd_P_0_b1bG2, 1, ktrue, &last_index, &P_0_b1bG2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P_2_b1bG2_new, dd_P_2_b1bG2, 1, ktrue, &last_index, &P_2_b1bG2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P_0_b2_new, dd_P_0_b2, 1, ktrue, &last_index, &P_0_b2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P_2_b2_new, dd_P_2_b2, 1, ktrue, &last_index, &P_2_b2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P_4_b2_new, dd_P_4_b2, 1, ktrue, &last_index, &P_4_b2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    class_call(array_interpolate_spline(kdisc, Nmax, P_0_bG2_new, dd_P_0_bG2, 1, ktrue, &last_index, &P_0_bG2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P_2_bG2_new, dd_P_2_bG2, 1, ktrue, &last_index, &P_2_bG2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    class_call(array_interpolate_spline(kdisc, Nmax, P_4_bG2_new, dd_P_4_bG2, 1, ktrue, &last_index, &P_4_bG2_ap_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);

                    //          Sigmatot = SigmaBAO*(1.+f*mutrue*mutrue*(2.+f))+f*f*mutrue*mutrue*(mutrue*mutrue-1.)*deltaSigmaBAO;
                    //          Exp = exp(-Sigmatot * pow(ktrue,2.));
                    //           P13ratio = 1.+(Pw_ap_out/Pnw_ap_out)*Exp;
                    //           P1b1 = (Pnw[index_j]+Pw[index_j]*Exp)*pnlpt->gauss_w[index_gauss2];
                    //           P1 = (Pnw[index_j]+Pw[index_j]*Exp)*f*pow(pnlpt->gauss_x[index_gauss2],2.)*pnlpt->gauss_w[index_gauss2];

                    Sigmatot = SigmaBAO * (1. + f * mutrue * mutrue * (2. + f)) + f * f * mutrue * mutrue * (mutrue * mutrue - 1.) * deltaSigmaBAO;
                    Exp = exp(-Sigmatot * pow(ktrue, 2.));
                    p_lo = (Pnw_ap_out + Pw_ap_out * Exp) / Pbin_ap_out;

                    IFG2_in = p_lo * P_IFG2_ap_out * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    Pd2d2_in = P_Id2d2_ap_out * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    Pd2G2_in = P_Id2G2_ap_out * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    PG2G2_in = P_IG2G2_ap_out * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    Pb1b2_in = (P_0_b1b2_ap_out + LegendreP2true * P_2_b1b2_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    Pb1bG2_in = (P_0_b1bG2_ap_out + LegendreP2true * P_2_b1bG2_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    Pb2_in = (P_0_b2_ap_out + LegendreP2true * P_2_b2_ap_out + LegendreP4true * P_4_b2_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    PbG2_in = (P_0_bG2_ap_out + LegendreP2true * P_2_bG2_ap_out + LegendreP4true * P_4_bG2_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    //GC!
                    
                    if (SWITCH_index == 1) {

                    t_lo = (Tnw_ap_out + Tw_ap_out * Exp) / Tbin_ap_out;

                    //GC -> recall what I have: P12_0_b1b2_ap_out, P12_0_b2_ap_out, P12_2_b2_ap_out, P12_0_b1bG2_ap_out, P12_0_bG2_ap_out, P12_2_bG2_ap_out...

                    P12b1b2_in = t_lo * (P12_0_b1b2_ap_out)*pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    P12b1bG2_in = t_lo * (P12_0_b1bG2_ap_out)*pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    P12b2_in = t_lo * (P12_0_b2_ap_out + LegendreP2true * P12_2_b2_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    P12bG2_in = t_lo * (P12_0_bG2_ap_out + LegendreP2true * P12_2_bG2_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    //GC: ORTHOGONAL -- start

                    P12b1b2_in_ortho = t_lo * (P12_0_b1b2_ap_out_ortho)*pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    P12b1bG2_in_ortho = t_lo * (P12_0_b1bG2_ap_out_ortho)*pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    P12b2_in_ortho = t_lo * (P12_0_b2_ap_out_ortho + LegendreP2true * P12_2_b2_ap_out_ortho) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    P12bG2_in_ortho = t_lo * (P12_0_bG2_ap_out_ortho + LegendreP2true * P12_2_bG2_ap_out_ortho) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    //GC: ORTHOGONAL -- finish

                    //GC!
                        
                    }

                    /*
            P1loopdd = ((p_tree + P22_mu0_dd_ap_out + P13_mu0_dd_ap_out*P13ratio + (P13_mu0_dd_w_ap_out+P22_mu0_dd_w_ap_out)*Exp)+(P22_mu2_dd_ap_out + P13_mu2_dd_ap_out*P13ratio + (P22_mu2_dd_w_ap_out+P13_mu2_dd_w_ap_out)*Exp)*pow(mutrue,2.)+(P22_mu4_dd_ap_out+P22_mu4_dd_w_ap_out*Exp)*pow(mutrue,4.))*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;

            P1loopvd = ((p_tree*2.*f + P13_mu2_vd_ap_out*P13ratio + P22_mu2_vd_ap_out+ (P22_mu2_vd_w_ap_out+P13_mu2_vd_w_ap_out)*Exp)*pow(mutrue,2.)+(P13_mu4_vd_ap_out*P13ratio + P22_mu4_vd_ap_out + (P22_mu4_vd_w_ap_out+P13_mu4_vd_w_ap_out)*Exp)*pow(mutrue,4.)+(P22_mu6_vd_ap_out+P22_mu6_vd_w_ap_out*Exp)*pow(mutrue,6.))*pnlpt->gauss_w[index_gauss2]*hratio/Dratio/Dratio;
            */
                    P_IFG2_0b1_x[index_j] += IFG2_in * LegendreP0 / 2.;
                    P_IFG2_0[index_j] += IFG2_in * f * pow(mutrue, 2.) * LegendreP0 / 2.;
                    P_IFG2_2[index_j] += IFG2_in * f * pow(mutrue, 2.) * LegendreP2 * 2.5;
                    P_Id2d2[index_j] += Pd2d2_in * LegendreP0 / 2.;
                    P_Id2d2_2[index_j] += Pd2d2_in * LegendreP2 * 2.5;
                    P_Id2d2_4[index_j] += Pd2d2_in * LegendreP4 * 4.5;
                    P_Id2G2[index_j] += Pd2G2_in * LegendreP0 / 2.;
                    P_Id2G2_2[index_j] += Pd2G2_in * LegendreP2 * 2.5;
                    P_Id2G2_4[index_j] += Pd2G2_in * LegendreP4 * 4.5;
                    P_IG2G2[index_j] += PG2G2_in * LegendreP0 / 2.;
                    P_IG2G2_2[index_j] += PG2G2_in * LegendreP2 * 2.5;
                    P_IG2G2_4[index_j] += PG2G2_in * LegendreP4 * 4.5;
                    P_0_b1b2[index_j] += Pb1b2_in * LegendreP0 / 2.;
                    P_2_b1b2[index_j] += Pb1b2_in * LegendreP2 * 2.5;
                    P_0_b1bG2[index_j] += Pb1bG2_in * LegendreP0 / 2.;
                    P_2_b1bG2[index_j] += Pb1bG2_in * LegendreP2 * 2.5;
                    P_0_b2[index_j] += Pb2_in * LegendreP0 / 2.;
                    P_2_b2[index_j] += Pb2_in * LegendreP2 * 2.5;
                    P_4_b2[index_j] += Pb2_in * LegendreP4 * 4.5;
                    P_0_bG2[index_j] += PbG2_in * LegendreP0 / 2.;
                    P_2_bG2[index_j] += PbG2_in * LegendreP2 * 2.5;
                    P_4_bG2[index_j] += PbG2_in * LegendreP4 * 4.5;
                    P_4_b1b2[index_j] += Pb1b2_in * LegendreP4 * 4.5;
                    P_4_b1bG2[index_j] += Pb1bG2_in * LegendreP4 * 4.5; //; //GC: why two ;;???

                    //GC!
                    
                    if (SWITCH_index == 1) {

                    P12_0_b2[index_j] += P12b2_in * LegendreP0 / 2.;
                    P12_2_b2[index_j] += P12b2_in * LegendreP2 * 2.5;
                    P12_4_b2[index_j] += P12b2_in * LegendreP4 * 4.5;

                    P12_0_bG2[index_j] += P12bG2_in * LegendreP0 / 2.;
                    P12_2_bG2[index_j] += P12bG2_in * LegendreP2 * 2.5;
                    P12_4_bG2[index_j] += P12bG2_in * LegendreP4 * 4.5;

                    P12_0_b1b2[index_j] += P12b1b2_in * LegendreP0 / 2.;
                    P12_2_b1b2[index_j] += P12b1b2_in * LegendreP2 * 2.5;
                    P12_4_b1b2[index_j] += P12b1b2_in * LegendreP4 * 4.5;

                    P12_0_b1bG2[index_j] += P12b1bG2_in * LegendreP0 / 2.;
                    P12_2_b1bG2[index_j] += P12b1bG2_in * LegendreP2 * 2.5;
                    P12_4_b1bG2[index_j] += P12b1bG2_in * LegendreP4 * 4.5;

                    //GC: ORTHOGONAL -- start

                    P12_0_b2_ortho[index_j] += P12b2_in_ortho * LegendreP0 / 2.;
                    P12_2_b2_ortho[index_j] += P12b2_in_ortho * LegendreP2 * 2.5;
                    P12_4_b2_ortho[index_j] += P12b2_in_ortho * LegendreP4 * 4.5;

                    P12_0_bG2_ortho[index_j] += P12bG2_in_ortho * LegendreP0 / 2.;
                    P12_2_bG2_ortho[index_j] += P12bG2_in_ortho * LegendreP2 * 2.5;
                    P12_4_bG2_ortho[index_j] += P12bG2_in_ortho * LegendreP4 * 4.5;

                    P12_0_b1b2_ortho[index_j] += P12b1b2_in_ortho * LegendreP0 / 2.;
                    P12_2_b1b2_ortho[index_j] += P12b1b2_in_ortho * LegendreP2 * 2.5;
                    P12_4_b1b2_ortho[index_j] += P12b1b2_in_ortho * LegendreP4 * 4.5;

                    P12_0_b1bG2_ortho[index_j] += P12b1bG2_in_ortho * LegendreP0 / 2.;
                    P12_2_b1bG2_ortho[index_j] += P12b1bG2_in_ortho * LegendreP2 * 2.5;
                    P12_4_b1bG2_ortho[index_j] += P12b1bG2_in_ortho * LegendreP4 * 4.5;
                        
                    }

                    //GC: ORTHOGONAL -- finish

                    //GC!
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

            //GC!
            
            //GC - SWITCH ~> here as I said many times I leave the freeing... One could be smarter but for now we just remove the biggest part, which is the import of matrices and calculation of the contributions...

            free(dd_Tbin);

            free(P12_0_b1b2_new);
            free(dd_P12_0_b1b2);

            free(P12_0_b1bG2_new);
            free(dd_P12_0_b1bG2);

            free(P12_0_b2_new);
            free(dd_P12_0_b2);

            free(P12_2_b2_new);
            free(dd_P12_2_b2);

            free(P12_0_bG2_new);
            free(dd_P12_0_bG2);

            free(P12_2_bG2_new);
            free(dd_P12_2_bG2);

            //GC: ORTHOGONAL -- start

            free(P12_0_b1b2_new_ortho);
            free(dd_P12_0_b1b2_ortho);

            free(P12_0_b1bG2_new_ortho);
            free(dd_P12_0_b1bG2_ortho);

            free(P12_0_b2_new_ortho);
            free(dd_P12_0_b2_ortho);

            free(P12_2_b2_new_ortho);
            free(dd_P12_2_b2_ortho);

            free(P12_0_bG2_new_ortho);
            free(dd_P12_0_bG2_ortho);

            free(P12_2_bG2_new_ortho);
            free(dd_P12_2_bG2_ortho);

            //GC: ORTHOGONAL -- finish

            //GC -> should be all we need...

            //GC -> now he starts interpolating the result. Notice that he starts from the multipoles that ARE NOT PRESENT before the AP, which would be the ones that didn't need to be "newed"... I start from the ones I already have...

            //GC: was newed...
            
            //GC - SWITCH -> this remains, as before...

            double pk12_0_b1b2_out = 0;
            double *ddpk12_0_b1b2;
            class_alloc(ddpk12_0_b1b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_b1b2, 1, ddpk12_0_b1b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_b1b2, ddpk12_0_b1b2, 1, pnlpt->k[index_k], &last_index, &pk12_0_b1b2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_0_b1b2[index_k] = pk12_0_b1b2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                    //+ 1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!! Now they have been defined... Below there still are some comments like this...
                }
                else
                {
                    pk12_l_0_b1b2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.24f %.24f\n",pnlpt->k[index_k],pk12_l_0_b1b2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_b1b2[index_k]-0.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_0_b1b2);
            free(f12_0_b1b2);
            free(P12_0_b1b2);

            //GC: ORTHOGONAL -- start

            double pk12_0_b1b2_out_ortho = 0;
            double *ddpk12_0_b1b2_ortho;
            class_alloc(ddpk12_0_b1b2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_b1b2_ortho, 1, ddpk12_0_b1b2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_b1b2_ortho, ddpk12_0_b1b2_ortho, 1, pnlpt->k[index_k], &last_index, &pk12_0_b1b2_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_0_b1b2_ortho[index_k] = pk12_0_b1b2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                    //+ 1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!! Now they have been defined... Below there still are some comments like this...
                }
                else
                {
                    pk12_l_0_b1b2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.24f %.24f\n",pnlpt->k[index_k],pk12_l_0_b1b2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_b1b2[index_k]-0.*large_for_logs_fNL); //GC!

                //printf("%.24e %.24e\n",pnlpt->k[index_k],pk12_l_0_b1b2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_0_b1b2_ortho);
            //free(f12_0_b1b2);
            free(P12_0_b1b2_ortho);

            //GC: ORTHOGONAL -- finish

            //GC: was newed...

            double pk12_0_b1bG2_out = 0;
            double *ddpk12_0_b1bG2;
            class_alloc(ddpk12_0_b1bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_b1bG2, 1, ddpk12_0_b1bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_b1bG2, ddpk12_0_b1bG2, 1, pnlpt->k[index_k], &last_index, &pk12_0_b1bG2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_0_b1bG2[index_k] = pk12_0_b1bG2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                    //+ 1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_0_b1bG2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_b1bG2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_b1bG2[index_k]-0.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_0_b1bG2);
            free(f12_0_b1bG2);
            free(P12_0_b1bG2);

            //GC: ORTHOGONAL -- start

            double pk12_0_b1bG2_out_ortho = 0;
            double *ddpk12_0_b1bG2_ortho;
            class_alloc(ddpk12_0_b1bG2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_b1bG2_ortho, 1, ddpk12_0_b1bG2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_b1bG2_ortho, ddpk12_0_b1bG2_ortho, 1, pnlpt->k[index_k], &last_index, &pk12_0_b1bG2_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_0_b1bG2_ortho[index_k] = pk12_0_b1bG2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                    //+ 1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_0_b1bG2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_b1bG2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_b1bG2[index_k]-0.*large_for_logs_fNL); //GC!

                //printf("%.24e %.24e\n",pnlpt->k[index_k],pk12_l_0_b1bG2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_0_b1bG2_ortho);
            //free(f12_0_b1bG2);
            free(P12_0_b1bG2_ortho);

            //GC: ORTHOGONAL -- finish

            //GC: was newed...

            double pk12_0_b2_out = 0;
            double *ddpk12_0_b2;
            class_alloc(ddpk12_0_b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_b2, 1, ddpk12_0_b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_b2, ddpk12_0_b2, 1, pnlpt->k[index_k], &last_index, &pk12_0_b2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_0_b2[index_k] = pk12_0_b2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                    //+ 1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_0_b2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.24f %.24f\n",pnlpt->k[index_k],pk12_l_0_b2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_b2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_b2[index_k]-0.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_0_b2);
            free(f12_0_b2);
            free(P12_0_b2);

            //GC: ORTHOGONAL -- start

            double pk12_0_b2_out_ortho = 0;
            double *ddpk12_0_b2_ortho;
            class_alloc(ddpk12_0_b2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_b2_ortho, 1, ddpk12_0_b2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_b2_ortho, ddpk12_0_b2_ortho, 1, pnlpt->k[index_k], &last_index, &pk12_0_b2_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_0_b2_ortho[index_k] = pk12_0_b2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                    //+ 1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_0_b2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.24f %.24f\n",pnlpt->k[index_k],pk12_l_0_b2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_b2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_b2[index_k]-0.*large_for_logs_fNL); //GC!

                //printf("%.24e %.24e\n",pnlpt->k[index_k],pk12_l_0_b2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_0_b2_ortho);
            //free(f12_0_b2);
            free(P12_0_b2_ortho);

            //GC: ORTHOGONAL -- finish

            //GC: was newed...

            double pk12_0_bG2_out = 0;
            double *ddpk12_0_bG2;
            class_alloc(ddpk12_0_bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_bG2, 1, ddpk12_0_bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_bG2, ddpk12_0_bG2, 1, pnlpt->k[index_k], &last_index, &pk12_0_bG2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_0_bG2[index_k] = pk12_0_bG2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                    //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_0_bG2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_bG2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_bG2[index_k]-0.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_0_bG2);
            free(f12_0_bG2);
            free(P12_0_bG2);

            //GC: ORTHOGONAL -- start

            double pk12_0_bG2_out_ortho = 0;
            double *ddpk12_0_bG2_ortho;
            class_alloc(ddpk12_0_bG2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_0_bG2_ortho, 1, ddpk12_0_bG2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_0_bG2_ortho, ddpk12_0_bG2_ortho, 1, pnlpt->k[index_k], &last_index, &pk12_0_bG2_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_0_bG2_ortho[index_k] = pk12_0_bG2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                    //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_0_bG2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_bG2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_0_bG2[index_k]-0.*large_for_logs_fNL); //GC!

                //printf("%.24e %.24e\n",pnlpt->k[index_k],pk12_l_0_bG2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_0_bG2_ortho);
            //free(f12_0_bG2);
            free(P12_0_bG2_ortho);

            //GC: ORTHOGONAL -- finish

            //GC: was newed...

            double pk12_2_b2_out = 0;
            double *ddpk12_2_b2;
            class_alloc(ddpk12_2_b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_2_b2, 1, ddpk12_2_b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_2_b2, ddpk12_2_b2, 1, pnlpt->k[index_k], &last_index, &pk12_2_b2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_2_b2[index_k] = pk12_2_b2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                    //+ 0.*1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_2_b2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 0.*1.e7; //GC!
                }

                //printf("%.16f %.16f\n", pnlpt->k[index_k], pk12_l_2_b2[index_k]); //GC -> check with Misha was here... Before it was "printf("%lf\n", pk_l_0_b1b2[index_k]);", as it is around everywhere...

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_b2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_b2[index_k]-0.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_2_b2);
            free(f12_2_b2);
            free(P12_2_b2);

            //GC: ORTHOGONAL -- start

            double pk12_2_b2_out_ortho = 0;
            double *ddpk12_2_b2_ortho;
            class_alloc(ddpk12_2_b2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_2_b2_ortho, 1, ddpk12_2_b2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_2_b2_ortho, ddpk12_2_b2_ortho, 1, pnlpt->k[index_k], &last_index, &pk12_2_b2_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_2_b2_ortho[index_k] = pk12_2_b2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                    //+ 0.*1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_2_b2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 0.*1.e7; //GC!
                }

                //printf("%.16f %.16f\n", pnlpt->k[index_k], pk12_l_2_b2[index_k]); //GC -> check with Misha was here... Before it was "printf("%lf\n", pk_l_0_b1b2[index_k]);", as it is around everywhere...

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_b2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_b2[index_k]-0.*large_for_logs_fNL); //GC!

                //printf("%.24e %.24e\n",pnlpt->k[index_k],pk12_l_2_b2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_2_b2_ortho);
            //free(f12_2_b2);
            free(P12_2_b2_ortho);

            //GC: ORTHOGONAL -- finish

            //GC: was newed...

            double pk12_2_bG2_out = 0;
            double *ddpk12_2_bG2;
            class_alloc(ddpk12_2_bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_2_bG2, 1, ddpk12_2_bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_2_bG2, ddpk12_2_bG2, 1, pnlpt->k[index_k], &last_index, &pk12_2_bG2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_2_bG2[index_k] = pk12_2_bG2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                    //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_2_bG2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                    //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_bG2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_bG2[index_k]-0.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_2_bG2);
            free(f12_2_bG2);
            free(P12_2_bG2);

            //GC: ORTHOGONAL -- start

            double pk12_2_bG2_out_ortho = 0;
            double *ddpk12_2_bG2_ortho;
            class_alloc(ddpk12_2_bG2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_2_bG2_ortho, 1, ddpk12_2_bG2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_2_bG2_ortho, ddpk12_2_bG2_ortho, 1, pnlpt->k[index_k], &last_index, &pk12_2_bG2_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_2_bG2_ortho[index_k] = pk12_2_bG2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                    //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_2_bG2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL;
                    //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_bG2[index_k]-0.*large_for_logs_fNL); //GC!
                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_bG2[index_k]-0.*large_for_logs_fNL); //GC!

                //printf("%.24e %.24e\n",pnlpt->k[index_k],pk12_l_2_bG2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_2_bG2_ortho);
            //free(f12_2_bG2);
            free(P12_2_bG2_ortho);

            //GC: ORTHOGONAL -- finish

            //GC -> now I go to the ones that weren't there before. I need to free less stuff in this case... Why doesn't he free the things called e.g. P_n_bG2 or whatever? He frees them at the end for some reason... Anyway, these are -> P12_4_b2, P12_4_bG2, P12_2_b1b2, P12_4_b1b2, P12_2_b1bG2, P12_4_b1bG2... For these, I just need to NOT free the f12 matrices, since they do not exist...

            double pk12_4_b2_out = 0;
            double *ddpk12_4_b2;
            class_alloc(ddpk12_4_b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_4_b2, 1, ddpk12_4_b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_4_b2, ddpk12_4_b2, 1, pnlpt->k[index_k], &last_index, &pk12_4_b2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_4_b2[index_k] = pk12_4_b2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_4_b2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!
                    //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_4_b2[index_k]-0.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_4_b2);
            free(P12_4_b2);

            //GC: ORTHOGONAL -- start

            double pk12_4_b2_out_ortho = 0;
            double *ddpk12_4_b2_ortho;
            class_alloc(ddpk12_4_b2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_4_b2_ortho, 1, ddpk12_4_b2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_4_b2_ortho, ddpk12_4_b2_ortho, 1, pnlpt->k[index_k], &last_index, &pk12_4_b2_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_4_b2_ortho[index_k] = pk12_4_b2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_4_b2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!
                    //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_4_b2[index_k]-0.*large_for_logs_fNL); //GC!

                //printf("%.24e %.24e\n",pnlpt->k[index_k],pk12_l_4_b2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_4_b2_ortho);
            free(P12_4_b2_ortho);

            //GC: ORTHOGONAL -- finish

            double pk12_4_bG2_out = 0;
            double *ddpk12_4_bG2;
            class_alloc(ddpk12_4_bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_4_bG2, 1, ddpk12_4_bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_4_bG2, ddpk12_4_bG2, 1, pnlpt->k[index_k], &last_index, &pk12_4_bG2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_4_bG2[index_k] = pk12_4_bG2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                    //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_4_bG2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!
                    //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_4_bG2[index_k]-0.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_4_bG2);
            free(P12_4_bG2);

            //GC: ORTHOGONAL -- start

            double pk12_4_bG2_out_ortho = 0;
            double *ddpk12_4_bG2_ortho;
            class_alloc(ddpk12_4_bG2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_4_bG2_ortho, 1, ddpk12_4_bG2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_4_bG2_ortho, ddpk12_4_bG2_ortho, 1, pnlpt->k[index_k], &last_index, &pk12_4_bG2_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_4_bG2_ortho[index_k] = pk12_4_bG2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                    //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_4_bG2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!
                    //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_4_bG2[index_k]-0.*large_for_logs_fNL); //GC!

                //printf("%.24e %.24e\n",pnlpt->k[index_k],pk12_l_4_bG2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_4_bG2_ortho);
            free(P12_4_bG2_ortho);

            //GC: ORTHOGONAL -- finish

            double pk12_2_b1b2_out = 0;
            double *ddpk12_2_b1b2;
            class_alloc(ddpk12_2_b1b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_2_b1b2, 1, ddpk12_2_b1b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_2_b1b2, ddpk12_2_b1b2, 1, pnlpt->k[index_k], &last_index, &pk12_2_b1b2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_2_b1b2[index_k] = pk12_2_b1b2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                    //+ 1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_2_b1b2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_b1b2[index_k]-0.*large_for_logs_fNL); //GC! TEST!!! It most likely adds this to zero, if AP is set to zero? No, there could be some random junk in it... In any case, this is irrelevant since we will only care about the case where there IS the AP effect...

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_b1b2[index_k]-0.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_2_b1b2);
            free(P12_2_b1b2);

            //GC: ORTHOGONAL -- start

            double pk12_2_b1b2_out_ortho = 0;
            double *ddpk12_2_b1b2_ortho;
            class_alloc(ddpk12_2_b1b2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_2_b1b2_ortho, 1, ddpk12_2_b1b2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_2_b1b2_ortho, ddpk12_2_b1b2_ortho, 1, pnlpt->k[index_k], &last_index, &pk12_2_b1b2_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_2_b1b2_ortho[index_k] = pk12_2_b1b2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL;
                    //+ 1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_2_b1b2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_b1b2[index_k]-0.*large_for_logs_fNL); //GC! TEST!!! It most likely adds this to zero, if AP is set to zero? No, there could be some random junk in it... In any case, this is irrelevant since we will only care about the case where there IS the AP effect...

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_b1b2[index_k]-0.*large_for_logs_fNL); //GC!

                //printf("%.24e %.24e\n",pnlpt->k[index_k],pk12_l_2_b1b2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_2_b1b2_ortho);
            free(P12_2_b1b2_ortho);

            //GC: ORTHOGONAL -- finish

            double pk12_4_b1b2_out = 0;
            double *ddpk12_4_b1b2;
            class_alloc(ddpk12_4_b1b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_4_b1b2, 1, ddpk12_4_b1b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_4_b1b2, ddpk12_4_b1b2, 1, pnlpt->k[index_k], &last_index, &pk12_4_b1b2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_4_b1b2[index_k] = pk12_4_b1b2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                    //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_4_b1b2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!
                    //epsilon_for_logs + 1.e7; //GC!
                }

                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_4_b1b2[index_k]-0.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_4_b1b2);
            free(P12_4_b1b2);

            //GC: ORTHOGONAL -- start

            double pk12_4_b1b2_out_ortho = 0;
            double *ddpk12_4_b1b2_ortho;
            class_alloc(ddpk12_4_b1b2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_4_b1b2_ortho, 1, ddpk12_4_b1b2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_4_b1b2_ortho, ddpk12_4_b1b2_ortho, 1, pnlpt->k[index_k], &last_index, &pk12_4_b1b2_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_4_b1b2_ortho[index_k] = pk12_4_b1b2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                    //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_4_b1b2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!
                    //epsilon_for_logs + 1.e7; //GC!
                }

                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_4_b1b2[index_k]-0.*large_for_logs_fNL); //GC!

                //printf("%.24e %.24e\n",pnlpt->k[index_k],pk12_l_4_b1b2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_4_b1b2_ortho);
            free(P12_4_b1b2_ortho);

            //GC: ORTHOGONAL -- finish

            double pk12_2_b1bG2_out = 0;
            double *ddpk12_2_b1bG2;
            class_alloc(ddpk12_2_b1bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_2_b1bG2, 1, ddpk12_2_b1bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_2_b1bG2, ddpk12_2_b1bG2, 1, pnlpt->k[index_k], &last_index, &pk12_2_b1bG2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_2_b1bG2[index_k] = pk12_2_b1bG2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                    //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_2_b1bG2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!
                    //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_b1bG2[index_k]-0.*large_for_logs_fNL); //GC! TEST!!!

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_b1bG2[index_k]-0.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_2_b1bG2);
            free(P12_2_b1bG2);

            //GC: ORTHOGONAL -- start

            double pk12_2_b1bG2_out_ortho = 0;
            double *ddpk12_2_b1bG2_ortho;
            class_alloc(ddpk12_2_b1bG2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_2_b1bG2_ortho, 1, ddpk12_2_b1bG2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_2_b1bG2_ortho, ddpk12_2_b1bG2_ortho, 1, pnlpt->k[index_k], &last_index, &pk12_2_b1bG2_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_2_b1bG2_ortho[index_k] = pk12_2_b1bG2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                    //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_2_b1bG2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //GC!
                    //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_b1bG2[index_k]-0.*large_for_logs_fNL); //GC! TEST!!!

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_2_b1bG2[index_k]-0.*large_for_logs_fNL); //GC!

                //printf("%.24e %.24e\n",pnlpt->k[index_k],pk12_l_2_b1bG2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_2_b1bG2_ortho);
            free(P12_2_b1bG2_ortho);

            //GC: ORTHOGONAL -- finish

            double pk12_4_b1bG2_out = 0;
            double *ddpk12_4_b1bG2;
            class_alloc(ddpk12_4_b1bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_4_b1bG2, 1, ddpk12_4_b1bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_4_b1bG2, ddpk12_4_b1bG2, 1, pnlpt->k[index_k], &last_index, &pk12_4_b1bG2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_4_b1bG2[index_k] = pk12_4_b1bG2_out + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                    //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_4_b1bG2[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_4_b1bG2[index_k]-0.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_4_b1bG2);
            free(P12_4_b1bG2);

            //GC: ORTHOGONAL -- start

            double pk12_4_b1bG2_out_ortho = 0;
            double *ddpk12_4_b1bG2_ortho;
            class_alloc(ddpk12_4_b1bG2_ortho, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P12_4_b1bG2_ortho, 1, ddpk12_4_b1bG2_ortho, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P12_4_b1bG2_ortho, ddpk12_4_b1bG2_ortho, 1, pnlpt->k[index_k], &last_index, &pk12_4_b1bG2_out_ortho, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk12_l_4_b1bG2_ortho[index_k] = pk12_4_b1bG2_out_ortho + large_for_logs_fNL; //+ 1.*epsilon_for_logs_fNL; //GC!
                    //1.e7; //GC!
                    //GC -> pk12_l_0_b1b2 and co. will need to be defined around line 2000!!!
                }
                else
                {
                    pk12_l_4_b1bG2_ortho[index_k] = large_for_logs_fNL + 1. * epsilon_for_logs_fNL; //epsilon_for_logs + 1.e7; //GC!
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk12_l_4_b1bG2[index_k]-0.*large_for_logs_fNL); //GC!

                //printf("%.24e %.24e\n",pnlpt->k[index_k],pk12_l_4_b1bG2_ortho[index_k]-1.*large_for_logs_fNL); //GC!
            }
            free(ddpk12_4_b1bG2_ortho);
            free(P12_4_b1bG2_ortho);

            //GC: ORTHOGONAL -- finish

            //GC!

            double pk_Id2d2_2_out = 0.;
            double *ddpk_PId2d2_2;
            class_alloc(ddpk_PId2d2_2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_Id2d2_2, 1, ddpk_PId2d2_2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P_Id2d2_2, ddpk_PId2d2_2, 1, pnlpt->k[index_k], &last_index, &pk_Id2d2_2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk_Id2d2_2[index_k] = pk_Id2d2_2_out + large_for_logs_big;
                    // printf("%lf\n",pk_Id2d2_2[index_k]);
                }
                else
                {
                    pk_Id2d2_2[index_k] = large_for_logs_big;
                }

                //printf("%.24e %.24e\n",pnlpt->k[index_k],pk_Id2d2_2[index_k]-1.e7); //GC!

                //GC: not best as check...
            }
            free(ddpk_PId2d2_2);

            double pk_Id2d2_4_out = 0;
            double *ddpk_PId2d2_4;
            class_alloc(ddpk_PId2d2_4, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_Id2d2_4, 1, ddpk_PId2d2_4, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P_Id2d2_4, ddpk_PId2d2_4, 1, pnlpt->k[index_k], &last_index, &pk_Id2d2_4_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk_Id2d2_4[index_k] = pk_Id2d2_4_out + large_for_logs_big;
                }
                else
                {
                    pk_Id2d2_4[index_k] = large_for_logs_big;
                }
            }
            free(ddpk_PId2d2_4);

            double pk_Id2G2_2_out = 0.;
            double *ddpk_PId2G2_2;
            class_alloc(ddpk_PId2G2_2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_Id2G2_2, 1, ddpk_PId2G2_2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P_Id2G2_2, ddpk_PId2G2_2, 1, pnlpt->k[index_k], &last_index, &pk_Id2G2_2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk_Id2G2_2[index_k] = pk_Id2G2_2_out + large_for_logs_big;
                }
                else
                {
                    pk_Id2G2_2[index_k] = 0. * epsilon_for_logs + large_for_logs_big;
                }
            }
            free(ddpk_PId2G2_2);

            double pk_Id2G2_4_out = 0.;
            double *ddpk_PId2G2_4;
            class_alloc(ddpk_PId2G2_4, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_Id2G2_4, 1, ddpk_PId2G2_4, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P_Id2G2_4, ddpk_PId2G2_4, 1, pnlpt->k[index_k], &last_index, &pk_Id2G2_4_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk_Id2G2_4[index_k] = pk_Id2G2_4_out + large_for_logs_big;
                }
                else
                {
                    pk_Id2G2_4[index_k] = 0. * epsilon_for_logs + large_for_logs_big;
                }
            }
            free(ddpk_PId2G2_4);

            double pk_IG2G2_2_out = 0.;
            double *ddpk_PIG2G2_2;
            class_alloc(ddpk_PIG2G2_2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_IG2G2_2, 1, ddpk_PIG2G2_2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P_IG2G2_2, ddpk_PIG2G2_2, 1, pnlpt->k[index_k], &last_index, &pk_IG2G2_2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk_IG2G2_2[index_k] = pk_IG2G2_2_out + large_for_logs_big;
                }
                else
                {
                    pk_IG2G2_2[index_k] = 0. * epsilon_for_logs + large_for_logs_big;
                }
                //        printf("%le %le %le\n", pk_Id2d2_2[index_k], pk_Id2G2_2[index_k],pk_IG2G2_2[index_k]);
            }
            free(ddpk_PIG2G2_2);

            double pk_IG2G2_4_out = 0.;
            double *ddpk_PIG2G2_4;
            class_alloc(ddpk_PIG2G2_4, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_IG2G2_4, 1, ddpk_PIG2G2_4, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P_IG2G2_4, ddpk_PIG2G2_4, 1, pnlpt->k[index_k], &last_index, &pk_IG2G2_4_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk_IG2G2_4[index_k] = pk_IG2G2_4_out + large_for_logs_big;
                }
                else
                {
                    pk_IG2G2_4[index_k] = 0. * epsilon_for_logs + large_for_logs_big;
                }
                //       printf("%lf\n",pk_IG2G2_4[index_k]);
                //       printf("%le %le %le\n", pk_Id2d2_4[index_k], pk_Id2G2_4[index_k],pk_IG2G2_4[index_k]);
            }
            free(ddpk_PIG2G2_4);

            double pk_0_b1b2_out = 0;
            double *ddpk_0_b1b2;
            class_alloc(ddpk_0_b1b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_0_b1b2, 1, ddpk_0_b1b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P_0_b1b2, ddpk_0_b1b2, 1, pnlpt->k[index_k], &last_index, &pk_0_b1b2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk_l_0_b1b2[index_k] = pk_0_b1b2_out + large_for_logs_big;
                }
                else
                {
                    pk_l_0_b1b2[index_k] = epsilon_for_logs + large_for_logs_big;
                }
                //      printf("%lf\n", pk_l_0_b1b2[index_k]);
            }
            free(ddpk_0_b1b2);
            free(f22_0_b1b2);
            free(P_0_b1b2);

            double pk_0_b2_out = 0;
            double *ddpk_0_b2;
            class_alloc(ddpk_0_b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc,
                                                  Nmax,
                                                  P_0_b2,
                                                  1,
                                                  ddpk_0_b2,
                                                  _SPLINE_NATURAL_,
                                                  pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);

            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {

                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
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
                    pk_l_0_b2[index_k] = pk_0_b2_out + large_for_logs_big;
                }
                else
                {
                    pk_l_0_b2[index_k] = epsilon_for_logs + large_for_logs_big;
                }
            }
            free(ddpk_0_b2);
            free(f22_0_b2);
            free(P_0_b2);

            double pk_0_b1bG2_out;
            double *ddpk_0_b1bG2;
            class_alloc(ddpk_0_b1bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc,
                                                  Nmax,
                                                  P_0_b1bG2,
                                                  1,
                                                  ddpk_0_b1bG2,
                                                  _SPLINE_NATURAL_,
                                                  pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);

            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {

                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
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
                    pk_l_0_b1bG2[index_k] = pk_0_b1bG2_out + large_for_logs_big;
                }
                else
                {
                    pk_l_0_b1bG2[index_k] = epsilon_for_logs + large_for_logs_big;
                }
            }
            free(ddpk_0_b1bG2);
            free(f22_0_b1bG2);
            free(P_0_b1bG2);

            double pk_2_b1b2_out;
            double *ddpk_2_b1b2;
            class_alloc(ddpk_2_b1b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc,
                                                  Nmax,
                                                  P_2_b1b2,
                                                  1,
                                                  ddpk_2_b1b2,
                                                  _SPLINE_NATURAL_,
                                                  pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);

            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {

                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
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
                    pk_l_2_b1b2[index_k] = pk_2_b1b2_out + large_for_logs_big;
                }
                else
                {
                    pk_l_2_b1b2[index_k] = epsilon_for_logs + large_for_logs_big;
                }
            }
            free(ddpk_2_b1b2);
            free(f22_2_b1b2);
            free(P_2_b1b2);

            double pk_4_b1b2_out;
            double *ddpk_4_b1b2;
            class_alloc(ddpk_4_b1b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_4_b1b2, 1, ddpk_4_b1b2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P_4_b1b2, ddpk_4_b1b2, 1, pnlpt->k[index_k], &last_index, &pk_4_b1b2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk_l_4_b1b2[index_k] = pk_4_b1b2_out + large_for_logs_big;
                }
                else
                {
                    pk_l_4_b1b2[index_k] = epsilon_for_logs + large_for_logs_big;
                }
            }
            free(ddpk_4_b1b2);

            double pk_0_bG2_out;
            double *ddpk_0_bG2;
            class_alloc(ddpk_0_bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc,
                                                  Nmax,
                                                  P_0_bG2,
                                                  1,
                                                  ddpk_0_bG2,
                                                  _SPLINE_NATURAL_,
                                                  pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);

            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {

                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
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
                    pk_l_0_bG2[index_k] = pk_0_bG2_out + large_for_logs_big;
                }
                else
                {
                    pk_l_0_bG2[index_k] = epsilon_for_logs + large_for_logs_big;
                }
            }
            free(ddpk_0_bG2);
            free(f22_0_bG2);
            free(P_0_bG2);

            double pk_2_b2_out;
            double *ddpk_2_b2;
            class_alloc(ddpk_2_b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc,
                                                  Nmax,
                                                  P_2_b2,
                                                  1,
                                                  ddpk_2_b2,
                                                  _SPLINE_NATURAL_,
                                                  pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);

            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {

                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
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
                    pk_l_2_b2[index_k] = pk_2_b2_out + large_for_logs_big;
                }
                else
                {
                    pk_l_2_b2[index_k] = epsilon_for_logs + large_for_logs_big;
                }
            }
            free(ddpk_2_b2);
            free(f22_2_b2);
            free(P_2_b2);

            double *ddpk_2_b1bG2;
            double pk_2_b1bG2_out;
            class_alloc(ddpk_2_b1bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc,
                                                  Nmax,
                                                  P_2_b1bG2,
                                                  1,
                                                  ddpk_2_b1bG2,
                                                  _SPLINE_NATURAL_,
                                                  pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);

            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
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
                    pk_l_2_b1bG2[index_k] = pk_2_b1bG2_out + large_for_logs_big;
                }
                else
                {
                    pk_l_2_b1bG2[index_k] = epsilon_for_logs + large_for_logs_big;
                }
            }
            free(ddpk_2_b1bG2);
            free(f22_2_b1bG2);
            free(P_2_b1bG2);

            double pk_4_b1bG2_out;
            double *ddpk_4_b1bG2;
            class_alloc(ddpk_4_b1bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_4_b1bG2, 1, ddpk_4_b1bG2, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P_4_b1bG2, ddpk_4_b1bG2, 1, pnlpt->k[index_k], &last_index, &pk_4_b1bG2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk_l_4_b1bG2[index_k] = pk_4_b1bG2_out + large_for_logs_big;
                }
                else
                {
                    pk_l_4_b1bG2[index_k] = epsilon_for_logs + large_for_logs_big;
                }
                //      printf("%le %le\n", pk_l_4_b1b2[index_k],pk_l_4_b1bG2[index_k]);
            }
            free(ddpk_4_b1bG2);

            double *ddpk_2_bG2;
            double pk_2_bG2_out;
            class_alloc(ddpk_2_bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc,
                                                  Nmax,
                                                  P_2_bG2,
                                                  1,
                                                  ddpk_2_bG2,
                                                  _SPLINE_NATURAL_,
                                                  pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);

            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
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
                    pk_l_2_bG2[index_k] = pk_2_bG2_out + large_for_logs_big;
                }
                else
                {
                    pk_l_2_bG2[index_k] = epsilon_for_logs + large_for_logs_big;
                }
            }
            free(ddpk_2_bG2);
            free(f22_2_bG2);
            free(P_2_bG2);

            double pk_4_b2_out;
            double *ddpk_4_b2;
            class_alloc(ddpk_4_b2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc,
                                                  Nmax,
                                                  P_4_b2,
                                                  1,
                                                  ddpk_4_b2,
                                                  _SPLINE_NATURAL_,
                                                  pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);

            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
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
                    pk_l_4_b2[index_k] = pk_4_b2_out + large_for_logs_big;
                }
                else
                {
                    pk_l_4_b2[index_k] = epsilon_for_logs + large_for_logs_big;
                }
            }
            free(ddpk_4_b2);
            free(f22_4_b2);
            free(P_4_b2);

            double pk_4_bG2_out;
            double *ddpk_4_bG2;
            class_alloc(ddpk_4_bG2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc,
                                                  Nmax,
                                                  P_4_bG2,
                                                  1,
                                                  ddpk_4_bG2,
                                                  _SPLINE_NATURAL_,
                                                  pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);

            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
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
                    pk_l_4_bG2[index_k] = pk_4_bG2_out + large_for_logs_big; //GC -> ORIGINAL!!!
                    //pk_l_4_bG2[index_k] = pk_4_bG2_out; //GC!
                }
                else
                {
                    pk_l_4_bG2[index_k] = epsilon_for_logs + large_for_logs_big; //GC -> ORIGINAL!!!
                    //pk_l_4_bG2[index_k] = epsilon_for_logs; //GC!
                }

                //printf("%.16e %.16e\n",pnlpt->k[index_k],pk_l_4_bG2[index_k]); //GC -> compare with Misha...
            }
            free(ddpk_4_bG2);
            free(f22_4_bG2);
            free(P_4_bG2);

            double fg2_out = 0;
            double *ddpk_IFG2_0b1_mmm;
            class_alloc(ddpk_IFG2_0b1_mmm, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_IFG2_0b1_x, 1, ddpk_IFG2_0b1_mmm, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            int last_index_2 = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P_IFG2_0b1_x, ddpk_IFG2_0b1_mmm, 1, pnlpt->k[index_k], &last_index_2, &fg2_out, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk_IFG2_0b1[index_k] = fg2_out + large_for_logs_big;
                    //    printf("%lf\n", fg2_out);
                }
                else
                {
                    pk_IFG2_0b1[index_k] = exp(lnpk_l[index_k]) * 0. + large_for_logs_big;
                }
                //   printf("%lf\n",pk_Id2d2_2[index_k]);
                //        printf("%le %le %le\n", pk_IFG2_0b1[index_k], pk_Id2d2_2[index_k],pk_Id2d2_4[index_k]);
            }
            free(ddpk_IFG2_0b1_mmm);
            free(P_IFG2_new);

            double fg2_out_2 = 0.;
            double *ddpk_IFG2_0;
            class_alloc(ddpk_IFG2_0, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc, Nmax, P_IFG2_0, 1, ddpk_IFG2_0, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
                    class_call(array_interpolate_spline(kdisc, Nmax, P_IFG2_0, ddpk_IFG2_0, 1, pnlpt->k[index_k], &last_index, &fg2_out_2, 1, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
                    pk_IFG2_0[index_k] = fg2_out_2 + large_for_logs_big;
                }
                else
                {
                    pk_IFG2_0[index_k] = exp(lnpk_l[index_k]) * f / 3. * 0. + large_for_logs_big;
                }
                //      printf("%lf\n", pk_IFG2_0[index_k]);
            }
            free(ddpk_IFG2_0);

            double fg2_out_3 = 0.;
            double *ddpk_IFG2_2;
            class_alloc(ddpk_IFG2_2, sizeof(double) * Nmax, pnlpt->error_message);
            class_call(array_spline_table_columns(kdisc,
                                                  Nmax,
                                                  P_IFG2_2,
                                                  1,
                                                  ddpk_IFG2_2,
                                                  _SPLINE_NATURAL_,
                                                  pnlpt->error_message),
                       pnlpt->error_message,
                       pnlpt->error_message);
            last_index = 0;
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
            {
                if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
                {
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
                    pk_IFG2_2[index_k] = fg2_out_3 + large_for_logs_big;
                }
                else
                {
                    pk_IFG2_2[index_k] = exp(lnpk_l[index_k]) * 2 * f / 3. * 0. + large_for_logs_big;
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

        } //end of RSD conditional expression

        //GC -> here he interpolates in the case of NO RSDs... So for me only TWO terms... Notice that I made a mistake... In the sense that the reason why he does not deallocate... No, wait a second... Why does he put these outside of "if (pnlpt->rsd == rsd_yes) {"? I think just to save calls of allocation. I already splined my real-space results above!!! So I can skip this part... It is put here but it should have been above... The "if" is put in a somehow awkward place...

        /* double kmaxnew = kdisc[Nmax-2];
        double kminnew = kdisc[1];*/

        double pk_Id2d2_out;
        double *ddpk_PId2d2;
        class_alloc(ddpk_PId2d2, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_Id2d2,
                                              1,
                                              ddpk_PId2d2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);

        last_index = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {
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
                pk_Id2d2[index_k] = pk_Id2d2_out + large_for_logs_big;
            }
            else
            {
                pk_Id2d2[index_k] = epsilon_for_logs + large_for_logs_big;
            }
        }
        free(ddpk_PId2d2);
        free(f22_Id2d2);
        free(P_Id2d2);

        double pk_Id2G2_out;
        double *ddpk_Id2G2;
        class_alloc(ddpk_Id2G2, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_Id2G2,
                                              1,
                                              ddpk_Id2G2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {

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
                pk_Id2G2[index_k] = pk_Id2G2_out + large_for_logs_big;
            }

            else
            {
                pk_Id2G2[index_k] = epsilon_for_logs + large_for_logs_big;
            }

            //printf("%.24e %.24e\n",pnlpt->k[index_k],pk_Id2G2[index_k]-1.e7); //GC!

            //GC: not the best of checks. Let's check matter...
        }
        free(ddpk_Id2G2);

        free(f22_Id2G2);
        free(P_Id2G2);

        double pk_IG2G2_out;
        double *ddpk_IG2G2;
        class_alloc(ddpk_IG2G2, sizeof(double) * Nmax, pnlpt->error_message);
        class_call(array_spline_table_columns(kdisc,
                                              Nmax,
                                              P_IG2G2,
                                              1,
                                              ddpk_IG2G2,
                                              _SPLINE_NATURAL_,
                                              pnlpt->error_message),
                   pnlpt->error_message,
                   pnlpt->error_message);
        last_index = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {

            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew)
            {
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
                pk_IG2G2[index_k] = pk_IG2G2_out + large_for_logs_big;
            }
            else
            {
                pk_IG2G2[index_k] = epsilon_for_logs + large_for_logs_big;
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

        //GC -> now I need to free my stuff...
        
        //GC -> SWITCH -> all allocations remained, so they must be freed...

        free(etam2_transfer);
        free(cmsym2_transfer);
        free(x2_transfer);

        /*
        int* addressOfX = &index_k;
        printf("%p\n", &index_k);
         */

    } // end of bias conditional expression

    else
    {

        if (pnlpt->nonlinear_pt_verbose > 0)
            printf("No bias tracers requested.\n");

        double epsilon_for_logs = 1. * pow(10., -6.);

        for (index_k = 0; index_k < pnlpt->k_size; index_k++)
        {
            pk_Id2d2[index_k] = epsilon_for_logs;
            pk_Id2[index_k] = epsilon_for_logs;
            pk_IG2[index_k] = epsilon_for_logs;
            pk_Id2G2[index_k] = epsilon_for_logs;
            pk_IG2G2[index_k] = epsilon_for_logs;
            pk_IFG2[index_k] = epsilon_for_logs;
        }
    }
    //int end1=clock();

    //if (pnlpt->nonlinear_pt_verbose > 0)
    //printf("All matrices are calculated in %d mus\n",end1-start1);
    
    //GC -> SWITCH -> all allocations remained, so they must be freed...

    free(js);
    free(kdisc);
    free(Pbin);
    free(y);
    free(Tbin); //GC!

    free(y_transfer); //GC!

    free(Pw);
    free(Pnw);
    free(dd_Pnw);
    free(dd_Pw);
    //GC!

    free(Tw);
    free(Tnw);
    free(dd_Tnw);
    free(dd_Tw);

    //GC -> for now let's keep this...

    free(P10b1);
    free(P10);
    free(P12);

    return _SUCCESS_;
}
