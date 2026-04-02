/**
 * @file nonlinear_pt.c
 * @brief One-loop perturbation theory engine for CLASS-PT.
 *
 * @author Mikhail M. Ivanov (original), with extensions by
 *         Chudaykin, Philcox, Simonovic, Cabass et al.
 * @date 17.06.2018
 *
 * Computes one-loop SPT power spectra for dark matter and biased tracers
 * in real space and redshift space, with:
 *  - FFTLog decomposition of loop integrals (EdS approximation)
 *  - IR resummation of BAO (Senatore & Zaldarriaga method)
 *  - Alcock-Paczynski projection
 *  - Non-Gaussian initial conditions (local/equilateral/orthogonal fNL)
 *
 * Pipeline: input -> background -> thermodynamics -> perturbations ->
 *           primordial -> nonlinear -> [this module] -> transfer -> spectra
 *
 * Reference: Chudaykin, Ivanov, Philcox, Simonovic (arXiv:2004.10607)
 */

#define _GNU_SOURCE
#include "nonlinear_pt.h"
#include "fft.h"
#include <time.h>
#include <omp.h>

/* OpenBLAS thread control — allows setting BLAS to single-threaded
 * so we can run many independent BLAS calls in parallel via OpenMP */
extern void openblas_set_num_threads(int num_threads);
extern int openblas_get_num_threads(void);

#define TIMER_DECL(name) double _t_##name = 0.; struct timespec _ts0_##name, _ts1_##name
#define TIMER_START(name) clock_gettime(CLOCK_MONOTONIC, &_ts0_##name)
#define TIMER_ADD(name) do { \
    clock_gettime(CLOCK_MONOTONIC, &_ts1_##name); \
    _t_##name += (_ts1_##name.tv_sec - _ts0_##name.tv_sec) \
               + 1e-9 * (_ts1_##name.tv_nsec - _ts0_##name.tv_nsec); \
} while (0)

/* === K-GRID CONSTRUCTION === */

/**
 * Build the wavenumber grid for the PT module.
 *
 * Constructs a non-uniform k-grid with dense sampling around the BAO
 * scale (k_rec = 2pi/r_s) and logarithmic spacing at high k.
 * Three regimes: k < k_max_cmb (linear steps), k_max_cmb..k_max_cl
 * (log steps with BAO refinement), k_max_cl..k_max (log steps).
 */
int perturb_get_k_list_nl(
    struct precision *ppr,
    struct background *pba,
    struct thermodynamics *pth,
    struct perturbations *ppt,
    struct nonlinear_pt *pnlpt) {
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

    if (ppt->has_scalars == _TRUE_) {
        k_min = ppr->k_min_tau0 / pba->conformal_age;

        /* first value */
        if (pba->sgnK == 0) {
            /* K<0 (flat)  : start close to zero */
            k_min = ppr->k_min_tau0 / pba->conformal_age;
        }
        else if (pba->sgnK == -1) {
            /* K<0 (open): start close to sqrt(-K)
             * (in transfer modules, for scalars, this will correspond to q close to zero;
             * for vectors and tensors, this value is even smaller than the minimum necessary value) */
            k_min = sqrt(-pba->K + pow(ppr->k_min_tau0 / pba->conformal_age / pth->angular_rescaling, 2));
        }
        else if (pba->sgnK == 1) {
            /* K>0 (closed): start from q=sqrt(k2+(1+m)K) equal to 3sqrt(K), i.e. k=sqrt((8-m)K) */
            k_min = sqrt((8. - 1.e-4) * pba->K);
        }

        /** - --> find k_max (as well as k_max_cmb[ppt->index_md_scalars], k_max_cl[ppt->index_md_scalars]) */

        k_rec = 2. * _PI_ / pth->rs_rec; /* comoving scale corresponding to sound horizon at recombination */

        k_max_cmb = k_min;
        k_max_cl = k_min;
        k_max = k_min;

        /* find k_max: */

        if ((ppt->has_pk_matter == _TRUE_) || (ppt->has_density_transfers == _TRUE_) || (ppt->has_velocity_transfers == _TRUE_))
            k_max = MAX(k_max, ppt->k_max_for_pk);

        if (ppt->has_nl_corrections_based_on_delta_m == _TRUE_)
            k_max = MAX(k_max, ppr->nonlinear_min_k_max);

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

        while (k < k_max_cmb) {
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

            scale2 = pow(pba->H0, 2) + fabs(pba->K); /* a_today=1 in v3 */

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

        while (k < k_max_cl) {
            k *= pow(10., 1. / (ppr->k_per_decade_for_pk + (ppr->k_per_decade_for_bao - ppr->k_per_decade_for_pk) * (1. - tanh(pow((log(k) - log(ppr->k_bao_center * k_rec)) / log(ppr->k_bao_width), 4)))));

            pnlpt->k[index_k] = k;
            index_k++;
        }

        pnlpt->k_size_cl = index_k;

        /* values until k_max */

        while (k < k_max) {
            k *= pow(10., 1. / (ppr->k_per_decade_for_pk + (ppr->k_per_decade_for_bao - ppr->k_per_decade_for_pk) * (1. - tanh(pow((log(k) - log(ppr->k_bao_center * k_rec)) / log(ppr->k_bao_width), 4)))));

            pnlpt->k[index_k] = k;
            index_k++;
        }

        pnlpt->k_size = index_k;

        class_realloc(pnlpt->k,
                      pnlpt->k_size * sizeof(double),
                      pnlpt->error_message);
    }

    /** - vector modes skipped */

    /** - tensor modes skipped */

    /** - If user asked for k_output_values, add those to all k lists: */
    if (ppt->k_output_values_num > 0) {
        /* Allocate storage */
        class_alloc(ppt->index_k_output_values, sizeof(double) * 1 * ppt->k_output_values_num, pnlpt->error_message);

        /** - --> Find indices in ppt->k[index_md] corresponding to 'k_output_values'.
        We are assuming that ppt->k is sorted and growing, and we have made sure
        that ppt->k_output_values is also sorted and growing.*/

        newk_size = pnlpt->k_size + ppt->k_output_values_num;

        class_alloc(tmp_k_list, sizeof(double) * newk_size, pnlpt->error_message);

        index_k = 0;
        index_k_output = 0;
        for (index_newk = 0; index_newk < newk_size; index_newk++) {
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

            if (add_k_output_value == _TRUE_) {
                tmp_k_list[index_newk] = ppt->k_output_values[index_k_output];
                ppt->index_k_output_values[0 * ppt->k_output_values_num + index_k_output] = index_newk;
                index_k_output++;
            } else {
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
    }

    /** - finally, find the global k_min and k_max for the ensemble of all modes scalars, vectors, tensors) */

    pnlpt->k_min = _HUGE_;
    pnlpt->k_max = 0.;
    if (ppt->has_scalars == _TRUE_) {
        pnlpt->k_min = MIN(pnlpt->k_min, pnlpt->k[0]);                 /* first value, inferred from perturbations structure */
        pnlpt->k_max = MAX(pnlpt->k_max, pnlpt->k[pnlpt->k_size - 1]); /* last value, inferred from perturbations structure */
    }

    return _SUCCESS_;
}

/* === MODULE INITIALIZATION === */

/**
 * Initialize the non-linear PT module.
 *
 * Reads precomputed kernel matrices (M22, M13, IFG2, fNL),
 * allocates output arrays, extracts the linear P(k) at each
 * requested redshift, and calls nonlinear_pt_loop() to compute
 * the one-loop spectra.
 */
/* Offset tables for data-driven allocation and storage of P(k) arrays
 * and M22/M13 complex kernel matrices. These use offsetof() to enable
 * a single loop to allocate/initialize all struct members, replacing
 * ~100 individual allocation calls. */
#include <stddef.h>

/** Offsets of all ln_pk_* output arrays in struct nonlinear_pt (87 total) */
static const size_t ln_pk_offsets[] = {
    offsetof(struct nonlinear_pt, ln_pk_nl),
    offsetof(struct nonlinear_pt, ln_pk_Id2d2),
    offsetof(struct nonlinear_pt, ln_pk_Id2d2_2),
    offsetof(struct nonlinear_pt, ln_pk_Id2d2_4),
    offsetof(struct nonlinear_pt, ln_pk_Id2),
    offsetof(struct nonlinear_pt, ln_pk_IG2),
    offsetof(struct nonlinear_pt, ln_pk_Id2G2),
    offsetof(struct nonlinear_pt, ln_pk_Id2G2_2),
    offsetof(struct nonlinear_pt, ln_pk_Id2G2_4),
    offsetof(struct nonlinear_pt, ln_pk_IG2G2),
    offsetof(struct nonlinear_pt, ln_pk_IG2G2_2),
    offsetof(struct nonlinear_pt, ln_pk_IG2G2_4),
    offsetof(struct nonlinear_pt, ln_pk_IFG2),
    offsetof(struct nonlinear_pt, ln_pk_IFG2_0b1),
    offsetof(struct nonlinear_pt, ln_pk_IFG2_0),
    offsetof(struct nonlinear_pt, ln_pk_IFG2_2),
    offsetof(struct nonlinear_pt, ln_pk_CTR),
    offsetof(struct nonlinear_pt, ln_pk_CTR_0),
    offsetof(struct nonlinear_pt, ln_pk_CTR_2),
    offsetof(struct nonlinear_pt, ln_pk_CTR_4),
    offsetof(struct nonlinear_pt, ln_pk_Tree),
    offsetof(struct nonlinear_pt, ln_pk_Tree_0_vv),
    offsetof(struct nonlinear_pt, ln_pk_Tree_0_vd),
    offsetof(struct nonlinear_pt, ln_pk_Tree_0_dd),
    offsetof(struct nonlinear_pt, ln_pk_Tree_2_vv),
    offsetof(struct nonlinear_pt, ln_pk_Tree_2_vd),
    offsetof(struct nonlinear_pt, ln_pk_Tree_4_vv),
    offsetof(struct nonlinear_pt, ln_pk_0_vv),
    offsetof(struct nonlinear_pt, ln_pk_0_vd),
    offsetof(struct nonlinear_pt, ln_pk_0_dd),
    offsetof(struct nonlinear_pt, ln_pk_2_vv),
    offsetof(struct nonlinear_pt, ln_pk_2_vd),
    offsetof(struct nonlinear_pt, ln_pk_2_dd),
    offsetof(struct nonlinear_pt, ln_pk_4_vv),
    offsetof(struct nonlinear_pt, ln_pk_4_vd),
    offsetof(struct nonlinear_pt, ln_pk_4_dd),
    offsetof(struct nonlinear_pt, ln_pk_0_b1b2),
    offsetof(struct nonlinear_pt, ln_pk_0_b2),
    offsetof(struct nonlinear_pt, ln_pk_0_b1bG2),
    offsetof(struct nonlinear_pt, ln_pk_0_bG2),
    offsetof(struct nonlinear_pt, ln_pk_2_b1b2),
    offsetof(struct nonlinear_pt, ln_pk_2_b2),
    offsetof(struct nonlinear_pt, ln_pk_2_b1bG2),
    offsetof(struct nonlinear_pt, ln_pk_2_bG2),
    offsetof(struct nonlinear_pt, ln_pk_4_b2),
    offsetof(struct nonlinear_pt, ln_pk_4_bG2),
    offsetof(struct nonlinear_pt, ln_pk_4_b1b2),
    offsetof(struct nonlinear_pt, ln_pk_4_b1bG2),
    offsetof(struct nonlinear_pt, ln_pk_fNL_0_vv),
    offsetof(struct nonlinear_pt, ln_pk_fNL_0_vd),
    offsetof(struct nonlinear_pt, ln_pk_fNL_0_dd),
    offsetof(struct nonlinear_pt, ln_pk_fNL_2_vv),
    offsetof(struct nonlinear_pt, ln_pk_fNL_2_vd),
    offsetof(struct nonlinear_pt, ln_pk_fNL_2_dd),
    offsetof(struct nonlinear_pt, ln_pk_fNL_4_vv),
    offsetof(struct nonlinear_pt, ln_pk_fNL_4_vd),
    offsetof(struct nonlinear_pt, ln_pk_fNL_4_dd),
    offsetof(struct nonlinear_pt, ln_pk_fNL_0_vv_ortho),
    offsetof(struct nonlinear_pt, ln_pk_fNL_0_vd_ortho),
    offsetof(struct nonlinear_pt, ln_pk_fNL_0_dd_ortho),
    offsetof(struct nonlinear_pt, ln_pk_fNL_2_vv_ortho),
    offsetof(struct nonlinear_pt, ln_pk_fNL_2_vd_ortho),
    offsetof(struct nonlinear_pt, ln_pk_fNL_2_dd_ortho),
    offsetof(struct nonlinear_pt, ln_pk_fNL_4_vv_ortho),
    offsetof(struct nonlinear_pt, ln_pk_fNL_4_vd_ortho),
    offsetof(struct nonlinear_pt, ln_pk_fNL_4_dd_ortho),
    offsetof(struct nonlinear_pt, ln_pk12_0_b1b2),
    offsetof(struct nonlinear_pt, ln_pk12_0_b2),
    offsetof(struct nonlinear_pt, ln_pk12_0_b1bG2),
    offsetof(struct nonlinear_pt, ln_pk12_0_bG2),
    offsetof(struct nonlinear_pt, ln_pk12_2_b1b2),
    offsetof(struct nonlinear_pt, ln_pk12_2_b2),
    offsetof(struct nonlinear_pt, ln_pk12_2_b1bG2),
    offsetof(struct nonlinear_pt, ln_pk12_2_bG2),
    offsetof(struct nonlinear_pt, ln_pk12_4_b1b2),
    offsetof(struct nonlinear_pt, ln_pk12_4_b2),
    offsetof(struct nonlinear_pt, ln_pk12_4_b1bG2),
    offsetof(struct nonlinear_pt, ln_pk12_4_bG2),
    offsetof(struct nonlinear_pt, ln_pk12_0_b1b2_ortho),
    offsetof(struct nonlinear_pt, ln_pk12_0_b2_ortho),
    offsetof(struct nonlinear_pt, ln_pk12_0_b1bG2_ortho),
    offsetof(struct nonlinear_pt, ln_pk12_0_bG2_ortho),
    offsetof(struct nonlinear_pt, ln_pk12_2_b1b2_ortho),
    offsetof(struct nonlinear_pt, ln_pk12_2_b2_ortho),
    offsetof(struct nonlinear_pt, ln_pk12_2_b1bG2_ortho),
    offsetof(struct nonlinear_pt, ln_pk12_2_bG2_ortho),
    offsetof(struct nonlinear_pt, ln_pk12_4_b1b2_ortho),
    offsetof(struct nonlinear_pt, ln_pk12_4_b2_ortho),
    offsetof(struct nonlinear_pt, ln_pk12_4_b1bG2_ortho),
    offsetof(struct nonlinear_pt, ln_pk12_4_bG2_ortho),
    offsetof(struct nonlinear_pt, ln_pk_nl_fNL),
    offsetof(struct nonlinear_pt, ln_pk_fNLd2),
    offsetof(struct nonlinear_pt, ln_pk_fNLG2),
    offsetof(struct nonlinear_pt, ln_pk_nl_fNL_ortho),
    offsetof(struct nonlinear_pt, ln_pk_fNLd2_ortho),
    offsetof(struct nonlinear_pt, ln_pk_fNLG2_ortho),
};
#define N_LN_PK (sizeof(ln_pk_offsets) / sizeof(ln_pk_offsets[0]))


static const size_t m22_complex_offsets[] = {
    offsetof(struct nonlinear_pt, M22_oneline_complex),
    offsetof(struct nonlinear_pt, M22basic_oneline_complex),
    offsetof(struct nonlinear_pt, M22_0_b1b2_oneline_complex),
    offsetof(struct nonlinear_pt, M22_0_b2_oneline_complex),
    offsetof(struct nonlinear_pt, M22_0_b1bG2_oneline_complex),
    offsetof(struct nonlinear_pt, M22_0_bG2_oneline_complex),
    offsetof(struct nonlinear_pt, M22_2_b1b2_oneline_complex),
    offsetof(struct nonlinear_pt, M22_2_b2_oneline_complex),
    offsetof(struct nonlinear_pt, M22_2_b1bG2_oneline_complex),
    offsetof(struct nonlinear_pt, M22_2_bG2_oneline_complex),
    offsetof(struct nonlinear_pt, M22_4_b2_oneline_complex),
    offsetof(struct nonlinear_pt, M22_4_bG2_oneline_complex),
    offsetof(struct nonlinear_pt, M_IG2G2),
    offsetof(struct nonlinear_pt, M_Id2),
    offsetof(struct nonlinear_pt, M_IG2),
    offsetof(struct nonlinear_pt, M_Id2G2),
    offsetof(struct nonlinear_pt, M22_oneline_0_vv_complex),
    offsetof(struct nonlinear_pt, M22_oneline_0_vd_complex),
    offsetof(struct nonlinear_pt, M22_oneline_0_dd_complex),
    offsetof(struct nonlinear_pt, M22_oneline_2_vv_complex),
    offsetof(struct nonlinear_pt, M22_oneline_2_vd_complex),
    offsetof(struct nonlinear_pt, M22_oneline_2_dd_complex),
    offsetof(struct nonlinear_pt, M22_oneline_4_vv_complex),
    offsetof(struct nonlinear_pt, M22_oneline_4_vd_complex),
    offsetof(struct nonlinear_pt, M22_oneline_4_dd_complex),
    offsetof(struct nonlinear_pt, M22_oneline_mu2_vd_complex),
    offsetof(struct nonlinear_pt, M22_oneline_mu2_dd_complex),
    offsetof(struct nonlinear_pt, M22_oneline_mu4_vd_complex),
    offsetof(struct nonlinear_pt, M22_oneline_mu4_vv_complex),
    offsetof(struct nonlinear_pt, M22_oneline_mu4_dd_complex),
    offsetof(struct nonlinear_pt, M22_oneline_mu6_vv_complex),
    offsetof(struct nonlinear_pt, M22_oneline_mu6_vd_complex),
    offsetof(struct nonlinear_pt, M22_oneline_mu8_complex),
    offsetof(struct nonlinear_pt, M12_oneline_complex),
    offsetof(struct nonlinear_pt, M12_oneline_complex_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vv0_f2),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vv0_f3),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vd0_f1),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vd0_f2),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_dd0_f0),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_dd0_f1),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vv2_f3),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vd2_f2),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vv4_f3),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vd4_f2),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vv0_f2_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vv0_f3_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vd0_f1_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vd0_f2_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_dd0_f0_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_dd0_f1_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vv2_f3_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vd2_f2_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vv4_f3_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_multipoles_vd4_f2_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_mu_powers_vd2_f1),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_mu_powers_vd2_f2),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_mu_powers_dd2_f1),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_mu_powers_vv4_f2),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_mu_powers_vd4_f2),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_mu_powers_vv6_f3),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_mu_powers_vd2_f1_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_mu_powers_vd2_f2_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_mu_powers_dd2_f1_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_mu_powers_vv4_f2_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_mu_powers_vd4_f2_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_matter_mu_powers_vv6_f3_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_bias_real_space_b2),
    offsetof(struct nonlinear_pt, M12_oneline_complex_bias_real_space_bG2),
    offsetof(struct nonlinear_pt, M12_oneline_complex_bias_real_space_b2_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_bias_real_space_bG2_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_bias_multipoles_b2_vv0_f1),
    offsetof(struct nonlinear_pt, M12_oneline_complex_bias_multipoles_bG2_vv0_f1),
    offsetof(struct nonlinear_pt, M12_oneline_complex_bias_multipoles_b2_vv0_f1_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_complex_bias_multipoles_bG2_vv0_f1_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_0_vv_complex),
    offsetof(struct nonlinear_pt, M12_oneline_0_vd_complex),
    offsetof(struct nonlinear_pt, M12_oneline_0_dd_complex),
    offsetof(struct nonlinear_pt, M12_oneline_2_vv_complex),
    offsetof(struct nonlinear_pt, M12_oneline_2_vd_complex),
    offsetof(struct nonlinear_pt, M12_oneline_2_dd_complex),
    offsetof(struct nonlinear_pt, M12_oneline_4_vv_complex),
    offsetof(struct nonlinear_pt, M12_oneline_4_vd_complex),
    offsetof(struct nonlinear_pt, M12_oneline_mu2_vd_complex),
    offsetof(struct nonlinear_pt, M12_oneline_mu2_dd_complex),
    offsetof(struct nonlinear_pt, M12_oneline_mu4_vv_complex),
    offsetof(struct nonlinear_pt, M12_oneline_mu4_vd_complex),
    offsetof(struct nonlinear_pt, M12_oneline_mu6_vv_complex),
    offsetof(struct nonlinear_pt, M12_oneline_0_vv_complex_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_0_vd_complex_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_0_dd_complex_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_2_vv_complex_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_2_vd_complex_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_2_dd_complex_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_4_vv_complex_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_4_vd_complex_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_mu2_vd_complex_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_mu2_dd_complex_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_mu4_vv_complex_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_mu4_vd_complex_ortho),
    offsetof(struct nonlinear_pt, M12_oneline_mu6_vv_complex_ortho),
    offsetof(struct nonlinear_pt, M12_2_bG2_oneline_complex),
    offsetof(struct nonlinear_pt, M12_2_b2_oneline_complex),
    offsetof(struct nonlinear_pt, M12_0_bG2_oneline_complex),
    offsetof(struct nonlinear_pt, M12_0_b1bG2_oneline_complex),
    offsetof(struct nonlinear_pt, M12_0_b2_oneline_complex),
    offsetof(struct nonlinear_pt, M12_0_b1b2_oneline_complex),
    offsetof(struct nonlinear_pt, M12_2_bG2_oneline_complex_ortho),
    offsetof(struct nonlinear_pt, M12_2_b2_oneline_complex_ortho),
    offsetof(struct nonlinear_pt, M12_0_bG2_oneline_complex_ortho),
    offsetof(struct nonlinear_pt, M12_0_b1bG2_oneline_complex_ortho),
    offsetof(struct nonlinear_pt, M12_0_b2_oneline_complex_ortho),
    offsetof(struct nonlinear_pt, M12_0_b1b2_oneline_complex_ortho),
    offsetof(struct nonlinear_pt, M_fNLd2),
    offsetof(struct nonlinear_pt, M_fNLG2),
    offsetof(struct nonlinear_pt, M_fNLd2_ortho),
    offsetof(struct nonlinear_pt, M_fNLG2_ortho),
};
#define N_M22_COMPLEX (sizeof(m22_complex_offsets) / sizeof(m22_complex_offsets[0]))

static const size_t m13_complex_offsets[] = {
    offsetof(struct nonlinear_pt, M13_oneline_complex),
    offsetof(struct nonlinear_pt, IFG2_oneline_complex),
    offsetof(struct nonlinear_pt, M13_0_vv_oneline_complex),
    offsetof(struct nonlinear_pt, M13_0_vd_oneline_complex),
    offsetof(struct nonlinear_pt, M13_0_dd_oneline_complex),
    offsetof(struct nonlinear_pt, M13_2_vv_oneline_complex),
    offsetof(struct nonlinear_pt, M13_2_vd_oneline_complex),
    offsetof(struct nonlinear_pt, M13_2_dd_oneline_complex),
    offsetof(struct nonlinear_pt, M13_4_vv_oneline_complex),
    offsetof(struct nonlinear_pt, M13_4_vd_oneline_complex),
    offsetof(struct nonlinear_pt, M13_mu2_dd_oneline_complex),
    offsetof(struct nonlinear_pt, M13_mu2_vd_oneline_complex),
    offsetof(struct nonlinear_pt, M13_mu4_vv_oneline_complex),
    offsetof(struct nonlinear_pt, M13_mu4_vd_oneline_complex),
    offsetof(struct nonlinear_pt, M13_mu6_oneline_complex),
};
#define N_M13_COMPLEX (sizeof(m13_complex_offsets) / sizeof(m13_complex_offsets[0]))

/* Compact spline setup: compute spline derivative table */
#define SPLINE_SETUP(x, n, data, dd, mode)                                    \
    class_call(array_spline_table_columns(x, n, data, 1, dd, mode,            \
                                          pnlpt->error_message),              \
               pnlpt->error_message, pnlpt->error_message)

/* Compact spline interpolation: evaluate spline at a single point */
#define SPLINE_EVAL(x, n, data, dd, xval, last_idx, out)                      \
    class_call(array_interpolate_spline(x, n, data, dd, 1, xval,              \
                                        last_idx, out, 1,                     \
                                        pnlpt->error_message),               \
               pnlpt->error_message, pnlpt->error_message)

int nonlinear_pt_init(
    struct precision *ppr,
    struct background *pba,
    struct thermodynamics *pth,
    struct perturbations *ppt,
    struct primordial *ppm,
    struct nonlinear_pt *pnlpt) {
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

    /* Primordial P(k) arrays, used to extract the transfer function T(k) for fNL */
    double *pPRIMk_l;
    double *lnpPRIMk_l;
    double *ddlnpPRIMk_l;
    double *pk_l_fNL_0_vv;
    double *pk_l_fNL_0_vd;
    double *pk_l_fNL_0_dd;
    double *pk_l_fNL_2_vv;
    double *pk_l_fNL_2_vd;
    double *pk_l_fNL_2_dd;
    double *pk_l_fNL_4_vv;
    double *pk_l_fNL_4_vd;
    double *pk_l_fNL_4_dd;
    double *pk_l_fNL_0_vv_ortho;
    double *pk_l_fNL_0_vd_ortho;
    double *pk_l_fNL_0_dd_ortho;
    double *pk_l_fNL_2_vv_ortho;
    double *pk_l_fNL_2_vd_ortho;
    double *pk_l_fNL_2_dd_ortho;
    double *pk_l_fNL_4_vv_ortho;
    double *pk_l_fNL_4_vd_ortho;
    double *pk_l_fNL_4_dd_ortho;
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
    double *pk_nl_fNL;
    double *pk_fNLd2;
    double *pk_fNLG2;
    double *pk_nl_fNL_ortho;
    double *pk_fNLd2_ortho;
    double *pk_fNLG2_ortho;

    short print_warning = _FALSE_;
    int last_index;
    double a, z;

    last_index = 0;

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

        pnlpt->index_md_scalars = ppt->index_md_scalars;
        pnlpt->ic_size = ppt->ic_size[ppt->index_md_scalars];
        pnlpt->tp_size = ppt->tp_size[ppt->index_md_scalars];

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
            class_alloc(pnlpt->k, pnlpt->k_size * sizeof(double), pnlpt->error_message);
            for (index_k = 0; index_k < pnlpt->k_size; index_k++) {
                pnlpt->k[index_k] = ppt->k[ppt->index_md_scalars][index_k];
            }
        }

        pnlpt->ln_k_size = pnlpt->k_size;
        class_alloc(pnlpt->ln_k, pnlpt->ln_k_size * sizeof(double), pnlpt->error_message);
        for (index_k = 0; index_k < pnlpt->k_size; index_k++) {
            pnlpt->ln_k[index_k] = log(pnlpt->k[index_k]);
        }

        pnlpt->tau_size = ppt->tau_size;
        class_alloc(pnlpt->tau, pnlpt->tau_size * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->ln_tau, pnlpt->tau_size * sizeof(double), pnlpt->error_message);
        for (index_tau = 0; index_tau < pnlpt->tau_size; index_tau++) {
            pnlpt->tau[index_tau] = ppt->tau_sampling[index_tau];
            pnlpt->ln_tau[index_tau] = log(ppt->tau_sampling[index_tau]);
        }

        /* --- Allocate output arrays for all one-loop spectra --- */
        class_alloc(pnlpt->nl_corr_density, pnlpt->tau_size * ppt->k_size[pnlpt->index_md_scalars] * sizeof(double), pnlpt->error_message);
        {
            double **_pk_ptrs[] = {
                &pk_l, &pk_l_0_vv, &pk_l_0_vd, &pk_l_0_dd,
                &pk_l_2_vv, &pk_l_2_vd, &pk_l_2_dd, &pk_l_4_vv, &pk_l_4_vd, &pk_l_4_dd,
                &pk_l_0_b1b2, &pk_l_0_b2, &pk_l_0_b1bG2, &pk_l_0_bG2,
                &pk_l_2_b1b2, &pk_l_2_b2, &pk_l_2_b1bG2, &pk_l_2_bG2,
                &pk_l_4_b2, &pk_l_4_bG2, &pk_l_4_b1b2, &pk_l_4_b1bG2,
                &pk_nl, &pk_Id2d2, &pk_Id2d2_2, &pk_Id2d2_4,
                &pk_Id2, &pk_IG2, &pk_Id2G2, &pk_Id2G2_2, &pk_Id2G2_4,
                &pk_IG2G2, &pk_IG2G2_2, &pk_IG2G2_4,
                &pk_IFG2, &pk_IFG2_0b1, &pk_IFG2_0, &pk_IFG2_2,
                &pk_CTR, &pk_CTR_0, &pk_CTR_2, &pk_CTR_4,
                &pk_Tree, &pk_Tree_0_vv, &pk_Tree_0_vd, &pk_Tree_0_dd,
                &pk_Tree_2_vv, &pk_Tree_2_vd, &pk_Tree_4_vv,
                &lnk_l, &lnpk_l, &ddlnpk_l, &pPRIMk_l, &lnpPRIMk_l, &ddlnpPRIMk_l,
                &pk_l_fNL_0_vv, &pk_l_fNL_0_vd, &pk_l_fNL_0_dd,
                &pk_l_fNL_2_vv, &pk_l_fNL_2_vd, &pk_l_fNL_2_dd,
                &pk_l_fNL_4_vv, &pk_l_fNL_4_vd, &pk_l_fNL_4_dd,
                &pk_l_fNL_0_vv_ortho, &pk_l_fNL_0_vd_ortho, &pk_l_fNL_0_dd_ortho,
                &pk_l_fNL_2_vv_ortho, &pk_l_fNL_2_vd_ortho, &pk_l_fNL_2_dd_ortho,
                &pk_l_fNL_4_vv_ortho, &pk_l_fNL_4_vd_ortho, &pk_l_fNL_4_dd_ortho,
                &pk12_l_0_b1b2, &pk12_l_0_b2, &pk12_l_0_b1bG2, &pk12_l_0_bG2,
                &pk12_l_2_b1b2, &pk12_l_2_b2, &pk12_l_2_b1bG2, &pk12_l_2_bG2,
                &pk12_l_4_b1b2, &pk12_l_4_b2, &pk12_l_4_b1bG2, &pk12_l_4_bG2,
                &pk12_l_0_b1b2_ortho, &pk12_l_0_b2_ortho, &pk12_l_0_b1bG2_ortho, &pk12_l_0_bG2_ortho,
                &pk12_l_2_b1b2_ortho, &pk12_l_2_b2_ortho, &pk12_l_2_b1bG2_ortho, &pk12_l_2_bG2_ortho,
                &pk12_l_4_b1b2_ortho, &pk12_l_4_b2_ortho, &pk12_l_4_b1bG2_ortho, &pk12_l_4_bG2_ortho,
                &pk_nl_fNL, &pk_fNLd2, &pk_fNLG2,
                &pk_nl_fNL_ortho, &pk_fNLd2_ortho, &pk_fNLG2_ortho
            };
            const size_t _pk_sz = pnlpt->k_size * sizeof(double);
            for (int ai_ = 0; ai_ < 103; ai_++)
                class_alloc(*_pk_ptrs[ai_], _pk_sz, pnlpt->error_message);
        }
        class_alloc(tau_req, pnlpt->z_pk_num * sizeof(double), pnlpt->error_message);

        int i_z = 0;
        for (i_z = 0; i_z < pnlpt->z_pk_num; i_z++) {
            class_call(background_tau_of_z(pba, pnlpt->z_pk[i_z], &tau_req[i_z]),
                       pba->error_message,
                       pnlpt->error_message);
        }

        /* --- Read precomputed PT kernel matrices from disk --- */
        /* M22 (P22 kernels), M13 (P13 kernels), IFG2 (bias),
         * M12 (fNL transfer), and associated RSD/bias/multipole matrices.
         * Matrix size is set by ppr->nmax_nlpt (128, 256, 512, or 1024). */

        char file2openM22[256];
        char file2openM13[256];
        char file2openM22basic[256];
        char file2openM13basic[256];

        int has_fNL = 0;

        if (pnlpt->fNL_equil_ortho_switch == fNL_equil_ortho_yes) {
            has_fNL = 1;
        }

        /* Determine resolution string for matrix file paths */
        const char *Nstr;
        const char *Nstr_M22;
        if (ppr->nmax_nlpt == 128) {
            Nstr = "N128";
            Nstr_M22 = "N128";
        }
        else if (ppr->nmax_nlpt == 512) {
            Nstr = "N512";
            Nstr_M22 = "N512";
        } else {
            Nstr = "N256";
            Nstr_M22 = "N256_packed";
        }

        sprintf(file2openM22, "%s/pt_matrices/M22oneline_%s.dat", __CLASSDIR__, Nstr_M22);
        sprintf(file2openM13, "%s/pt_matrices/M13oneline_%s.dat", __CLASSDIR__, Nstr);
        sprintf(file2openM22basic, "%s/pt_matrices/M22basiconeline_%s.dat", __CLASSDIR__, Nstr_M22);
        sprintf(file2openM13basic, "%s/pt_matrices/IFG2oneline_%s.dat", __CLASSDIR__, Nstr);


        /* --- Macros for matrix I/O and real-to-complex conversion --- */

#define LOAD_MATRIX(array, filepath, size)                                   \
    do {                                                                     \
        class_alloc(array, (size) * sizeof(double), pnlpt->error_message);   \
        char binpath_[512];                                                  \
        snprintf(binpath_, sizeof(binpath_), "%s.bin", filepath);            \
        FILE *fptr_lm_ = fopen(binpath_, "rb");                              \
        if (fptr_lm_ != NULL) {                                              \
            size_t nread_ = fread((array), sizeof(double), (size), fptr_lm_);\
            fclose(fptr_lm_);                                                \
            if ((int)nread_ != (size)) {                                     \
                printf("Error: binary read %s: got %zu of %d\n",             \
                       binpath_, nread_, (size));                             \
                exit(1);                                                     \
            }                                                                \
        } else {                                                             \
            fptr_lm_ = fopen(filepath, "r");                                 \
            if (fptr_lm_ == NULL) {                                          \
                printf("Error: cannot open %s\n", filepath);                 \
                exit(1);                                                     \
            }                                                                \
            for (int idx_lm_ = 0; idx_lm_ < (size); idx_lm_++)              \
                fscanf(fptr_lm_, "%lf", &(array)[idx_lm_]);                  \
            fclose(fptr_lm_);                                                \
        }                                                                    \
    } while (0)

#define CONVERT_REAL_TO_COMPLEX_M22(dest, src, Np1)                                         \
    for (count = 0; count < (Np1) * ((Np1) + 1) / 2; count++)                               \
    {                                                                                       \
        (dest)[count] = (src)[count] + _Complex_I * (src)[count + (Np1) * ((Np1) + 1) / 2]; \
    }

#define CONVERT_REAL_TO_COMPLEX_M13(dest, src, Np1)                       \
    for (count = 0; count < (Np1); count++)                               \
    {                                                                     \
        (dest)[count] = (src)[count] + _Complex_I * (src)[count + (Np1)]; \
    }

/* Store log(pk) into pnlpt struct array at redshift slice i_z */
#define STORE_LN_PK(member, pk_arr) \
    pnlpt->member[i_z * pnlpt->k_size + index_k] = log(pk_arr[index_k])

        LOAD_MATRIX(pnlpt->M22_oneline, file2openM22, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)));

        LOAD_MATRIX(pnlpt->M13_oneline, file2openM13, ((ppr->nmax_nlpt + 1) * 2));

        LOAD_MATRIX(pnlpt->M22basic_oneline, file2openM22basic, ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2)));

        LOAD_MATRIX(pnlpt->IFG2_oneline, file2openM13basic, ((ppr->nmax_nlpt + 1) * 2));

        /* M12 matrix loading — bundled binary for fNL, zero-alloc otherwise */
        {
            const int m12_size = (ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2);
            const int n_m12 = 21;  /* must match M12_PAIRS in bundle_m12.py */
            double **m12_regular[21] = {
                &pnlpt->M12_oneline, &pnlpt->M12_oneline_matter_multipoles_vv0_f2,
                &pnlpt->M12_oneline_matter_multipoles_vv0_f3, &pnlpt->M12_oneline_matter_multipoles_vd0_f1,
                &pnlpt->M12_oneline_matter_multipoles_vd0_f2, &pnlpt->M12_oneline_matter_multipoles_dd0_f0,
                &pnlpt->M12_oneline_matter_multipoles_dd0_f1, &pnlpt->M12_oneline_matter_multipoles_vv2_f3,
                &pnlpt->M12_oneline_matter_multipoles_vd2_f2, &pnlpt->M12_oneline_matter_multipoles_vv4_f3,
                &pnlpt->M12_oneline_matter_multipoles_vd4_f2, &pnlpt->M12_oneline_matter_mu_powers_vd2_f1,
                &pnlpt->M12_oneline_matter_mu_powers_vd2_f2, &pnlpt->M12_oneline_matter_mu_powers_dd2_f1,
                &pnlpt->M12_oneline_matter_mu_powers_vv4_f2, &pnlpt->M12_oneline_matter_mu_powers_vd4_f2,
                &pnlpt->M12_oneline_matter_mu_powers_vv6_f3, &pnlpt->M12_oneline_bias_real_space_b2,
                &pnlpt->M12_oneline_bias_real_space_bG2, &pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1,
                &pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1,
            };
            double **m12_ortho[21] = {
                &pnlpt->M12_oneline_ortho, &pnlpt->M12_oneline_matter_multipoles_vv0_f2_ortho,
                &pnlpt->M12_oneline_matter_multipoles_vv0_f3_ortho, &pnlpt->M12_oneline_matter_multipoles_vd0_f1_ortho,
                &pnlpt->M12_oneline_matter_multipoles_vd0_f2_ortho, &pnlpt->M12_oneline_matter_multipoles_dd0_f0_ortho,
                &pnlpt->M12_oneline_matter_multipoles_dd0_f1_ortho, &pnlpt->M12_oneline_matter_multipoles_vv2_f3_ortho,
                &pnlpt->M12_oneline_matter_multipoles_vd2_f2_ortho, &pnlpt->M12_oneline_matter_multipoles_vv4_f3_ortho,
                &pnlpt->M12_oneline_matter_multipoles_vd4_f2_ortho, &pnlpt->M12_oneline_matter_mu_powers_vd2_f1_ortho,
                &pnlpt->M12_oneline_matter_mu_powers_vd2_f2_ortho, &pnlpt->M12_oneline_matter_mu_powers_dd2_f1_ortho,
                &pnlpt->M12_oneline_matter_mu_powers_vv4_f2_ortho, &pnlpt->M12_oneline_matter_mu_powers_vd4_f2_ortho,
                &pnlpt->M12_oneline_matter_mu_powers_vv6_f3_ortho, &pnlpt->M12_oneline_bias_real_space_b2_ortho,
                &pnlpt->M12_oneline_bias_real_space_bG2_ortho, &pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1_ortho,
                &pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1_ortho,
            };

            if (has_fNL) {
                /* Load from bundled binary files (2 fopen instead of 42) */
                char bpath_[256];
                sprintf(bpath_, "%s/pt_matrices/M12_bundle_%s.bin", __CLASSDIR__, Nstr);
                FILE *fb_ = fopen(bpath_, "rb");
                if (fb_ == NULL) { printf("Error: cannot open %s\n", bpath_); exit(1); }
                for (int m12i_ = 0; m12i_ < n_m12; m12i_++) {
                    class_alloc(*m12_regular[m12i_], m12_size * sizeof(double), pnlpt->error_message);
                    if (fread(*m12_regular[m12i_], sizeof(double), m12_size, fb_) != (size_t)m12_size) {
                        printf("Error: short read from %s at matrix %d\n", bpath_, m12i_); exit(1);
                    }
                }
                fclose(fb_);

                sprintf(bpath_, "%s/pt_matrices/M12_bundle_%s_ortho.bin", __CLASSDIR__, Nstr);
                fb_ = fopen(bpath_, "rb");
                if (fb_ == NULL) { printf("Error: cannot open %s\n", bpath_); exit(1); }
                for (int m12i_ = 0; m12i_ < n_m12; m12i_++) {
                    class_alloc(*m12_ortho[m12i_], m12_size * sizeof(double), pnlpt->error_message);
                    if (fread(*m12_ortho[m12i_], sizeof(double), m12_size, fb_) != (size_t)m12_size) {
                        printf("Error: short read from %s at matrix %d\n", bpath_, m12i_); exit(1);
                    }
                }
                fclose(fb_);
            } else {
                for (int m12i_ = 0; m12i_ < n_m12; m12i_++) {
                    class_alloc(*m12_regular[m12i_], m12_size * sizeof(double), pnlpt->error_message);
                    class_alloc(*m12_ortho[m12i_], m12_size * sizeof(double), pnlpt->error_message);
                }
            }
        }

        char file2openGauss[256];
        sprintf(file2openGauss, "%s/pt_matrices/gauss_tab.dat", __CLASSDIR__);
        int index_gauss = 0;

        class_alloc(pnlpt->gauss_w, 40 * sizeof(double), pnlpt->error_message);
        class_alloc(pnlpt->gauss_x, 40 * sizeof(double), pnlpt->error_message);
        LOAD_MATRIX(pnlpt->gauss, file2openGauss, 80);

        for (index_gauss = 0; index_gauss < 40; index_gauss++) {
            pnlpt->gauss_x[index_gauss] = pnlpt->gauss[index_gauss];
            pnlpt->gauss_w[index_gauss] = pnlpt->gauss[40 + index_gauss];
        }

        /* --- Convert real-valued kernel matrices to packed complex (Hermitian) form --- */
        /* The complex matrices are stored in LAPACK packed format for use
         * with zspmv (symmetric packed matrix-vector multiply). */

        /* Allocate all M22-type complex matrices */
        {
            const int m22_csize = ((ppr->nmax_nlpt + 1) * (ppr->nmax_nlpt + 2) / 2) * sizeof(complex double);
            for (int ci_ = 0; ci_ < N_M22_COMPLEX; ci_++) {
                complex double **ptr = (complex double **)((char *)pnlpt + m22_complex_offsets[ci_]);
                class_alloc(*ptr, m22_csize, pnlpt->error_message);
            }
            const int m13_csize = (ppr->nmax_nlpt + 1) * sizeof(complex double);
            for (int ci_ = 0; ci_ < N_M13_COMPLEX; ci_++) {
                complex double **ptr = (complex double **)((char *)pnlpt + m13_complex_offsets[ci_]);
                class_alloc(*ptr, m13_csize, pnlpt->error_message);
            }
        }

        int count = 0;

        CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M22_oneline_complex, pnlpt->M22_oneline, ppr->nmax_nlpt + 1);
        CONVERT_REAL_TO_COMPLEX_M13(pnlpt->M13_oneline_complex, pnlpt->M13_oneline, ppr->nmax_nlpt + 1);
        CONVERT_REAL_TO_COMPLEX_M13(pnlpt->IFG2_oneline_complex, pnlpt->IFG2_oneline, ppr->nmax_nlpt + 1);
        CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M22basic_oneline_complex, pnlpt->M22basic_oneline, ppr->nmax_nlpt + 1);
        if (has_fNL) {
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex, pnlpt->M12_oneline, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_ortho, pnlpt->M12_oneline_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2, pnlpt->M12_oneline_matter_multipoles_vv0_f2, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vv0_f3, pnlpt->M12_oneline_matter_multipoles_vv0_f3, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1, pnlpt->M12_oneline_matter_multipoles_vd0_f1, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vd0_f2, pnlpt->M12_oneline_matter_multipoles_vd0_f2, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_dd0_f0, pnlpt->M12_oneline_matter_multipoles_dd0_f0, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1, pnlpt->M12_oneline_matter_multipoles_dd0_f1, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vv2_f3, pnlpt->M12_oneline_matter_multipoles_vv2_f3, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vd2_f2, pnlpt->M12_oneline_matter_multipoles_vd2_f2, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vv4_f3, pnlpt->M12_oneline_matter_multipoles_vv4_f3, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vd4_f2, pnlpt->M12_oneline_matter_multipoles_vd4_f2, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2_ortho, pnlpt->M12_oneline_matter_multipoles_vv0_f2_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vv0_f3_ortho, pnlpt->M12_oneline_matter_multipoles_vv0_f3_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1_ortho, pnlpt->M12_oneline_matter_multipoles_vd0_f1_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vd0_f2_ortho, pnlpt->M12_oneline_matter_multipoles_vd0_f2_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_dd0_f0_ortho, pnlpt->M12_oneline_matter_multipoles_dd0_f0_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1_ortho, pnlpt->M12_oneline_matter_multipoles_dd0_f1_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vv2_f3_ortho, pnlpt->M12_oneline_matter_multipoles_vv2_f3_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vd2_f2_ortho, pnlpt->M12_oneline_matter_multipoles_vd2_f2_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vv4_f3_ortho, pnlpt->M12_oneline_matter_multipoles_vv4_f3_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_multipoles_vd4_f2_ortho, pnlpt->M12_oneline_matter_multipoles_vd4_f2_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f1, pnlpt->M12_oneline_matter_mu_powers_vd2_f1, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2, pnlpt->M12_oneline_matter_mu_powers_vd2_f2, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_mu_powers_dd2_f1, pnlpt->M12_oneline_matter_mu_powers_dd2_f1, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_mu_powers_vv4_f2, pnlpt->M12_oneline_matter_mu_powers_vv4_f2, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_mu_powers_vd4_f2, pnlpt->M12_oneline_matter_mu_powers_vd4_f2, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_mu_powers_vv6_f3, pnlpt->M12_oneline_matter_mu_powers_vv6_f3, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f1_ortho, pnlpt->M12_oneline_matter_mu_powers_vd2_f1_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2_ortho, pnlpt->M12_oneline_matter_mu_powers_vd2_f2_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_mu_powers_dd2_f1_ortho, pnlpt->M12_oneline_matter_mu_powers_dd2_f1_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_mu_powers_vv4_f2_ortho, pnlpt->M12_oneline_matter_mu_powers_vv4_f2_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_mu_powers_vd4_f2_ortho, pnlpt->M12_oneline_matter_mu_powers_vd4_f2_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_matter_mu_powers_vv6_f3_ortho, pnlpt->M12_oneline_matter_mu_powers_vv6_f3_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_bias_real_space_b2, pnlpt->M12_oneline_bias_real_space_b2, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_bias_real_space_bG2, pnlpt->M12_oneline_bias_real_space_bG2, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_bias_real_space_b2_ortho, pnlpt->M12_oneline_bias_real_space_b2_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_bias_real_space_bG2_ortho, pnlpt->M12_oneline_bias_real_space_bG2_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1, pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1, pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1_ortho, pnlpt->M12_oneline_bias_multipoles_b2_vv0_f1_ortho, ppr->nmax_nlpt + 1);
            CONVERT_REAL_TO_COMPLEX_M22(pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1_ortho, pnlpt->M12_oneline_bias_multipoles_bG2_vv0_f1_ortho, ppr->nmax_nlpt + 1);
        }

        /* --- All PT kernel matrices loaded; allocate output storage arrays --- */
        /* Output arrays for spectra at each redshift (accessed by classy via spectra_pk_nl_bias_at_z) */
        /* Allocate all ln_pk output arrays */
        {
            const size_t ln_pk_size = sizeof(double) * pnlpt->z_pk_num * pnlpt->k_size;
            for (int lpi_ = 0; lpi_ < N_LN_PK; lpi_++) {
                double **ptr = (double **)((char *)pnlpt + ln_pk_offsets[lpi_]);
                class_alloc(*ptr, ln_pk_size, pnlpt->error_message);
            }
        }

        int index_md;
        index_md = pnlpt->index_md_scalars;
        int index_ic = 0;

        /* --- Fast output: spline sources for interpolation to arbitrary tau --- */
        if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {
            if (pnlpt->cb == _TRUE_) {
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

                for (index_tau = pnlpt->tau_size - 1; index_tau >= 0; index_tau--) {
                    for (index_k = 0; index_k < pnlpt->k_size; index_k++) {
                        if (pnlpt->k[index_k] <= ppt->k[index_md][ppt->k_size[index_md] - 1]) {
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

                        } /* end k <= k_max condition */
                        /* Crude padding for k beyond perturbation module range (only relevant for non-zero curvature) */
                        /* This is done only for the non-zero spatial curvature, for which the kmax from perturbations does not coincude with the kmax from the spectra */
                        else
                        {
                            pnlpt->sources_tp_delta_cb[index_tau * pnlpt->k_size + index_k] = pnlpt->sources_tp_delta_cb[index_tau * pnlpt->k_size + index_k - 1];
                        }
                    }
                }
            } else {
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
                for (index_tau = pnlpt->tau_size - 1; index_tau >= 0; index_tau--) {
                    for (index_k = 0; index_k < pnlpt->k_size; index_k++) {
                        if (pnlpt->k[index_k] <= ppt->k[index_md][ppt->k_size[index_md] - 1]) {
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
                        } /* end k <= k_max condition */
                        /* Crude padding for k beyond perturbation module range (only relevant for non-zero curvature) */
                        /* This is done only for curvature, for which the kmax from perturbations does not coincide with the kmax from the spectra */
                        else
                        {
                            pnlpt->sources_tp_delta_m[index_tau * pnlpt->k_size + index_k] = pnlpt->sources_tp_delta_m[index_tau * pnlpt->k_size + index_k - 1];
                        }
                    }
                }
            }
            /* End cb (cold baryon) transfer function interpolation */
        }

        /* End fast_output source interpolation */

        /* Build full P_L(k, tau), ln P_L, and spline tables over all tau */
        double *lnpk_l_full;
        double *pk_l_full;
        double *ddlnpk_l_full;

        class_alloc(lnpk_l_full, sizeof(double) * pnlpt->tau_size * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pk_l_full, sizeof(double) * pnlpt->tau_size * pnlpt->k_size, pnlpt->error_message);
        class_alloc(ddlnpk_l_full, sizeof(double) * pnlpt->tau_size * pnlpt->k_size, pnlpt->error_message);

        /* Loop over conformal time to fill P_L(k,tau) arrays */
        for (index_tau = pnlpt->tau_size - 1; index_tau >= 0; index_tau--) {
            /* get P_L(k) at this time */

            class_call(nonlinear_pt_pk_l(pba, ppt, ppm, pnlpt, index_tau, pk_l, lnk_l, lnpk_l, ddlnpk_l),
                       pnlpt->error_message,
                       pnlpt->error_message);

            /* get P_L(k,tau) lnP_L(k,tau) and ddP_L(k,tau) */

            for (index_k = 0; index_k < pnlpt->k_size; index_k++) {
                lnpk_l_full[index_tau * pnlpt->k_size + index_k] = lnpk_l[index_k];
                pk_l_full[index_tau * pnlpt->k_size + index_k] = pk_l[index_k];
                ddlnpk_l_full[index_tau * pnlpt->k_size + index_k] = ddlnpk_l[index_k];
            }
        }
        /*end cycle over tau*/

        class_call(nonlinear_pt_pPRIMk_l(
                       ppt, ppm, pnlpt, lnk_l,
                       pPRIMk_l, lnpPRIMk_l, ddlnpPRIMk_l),
                   pnlpt->error_message,
                   pnlpt->error_message);

        /* --- Interpolate P_L(k) and P_PRIM(k) to each requested redshift --- */
        double *ln_pk_l_at_z_req;
        double *pk_l_at_z_req;
        double *ln_pPRIMk_l_req;
        double *pPRIMk_l_req;

        class_alloc(ln_pPRIMk_l_req, sizeof(double) * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pPRIMk_l_req, sizeof(double) * pnlpt->k_size, pnlpt->error_message);

        double Dref;

        class_alloc(ln_pk_l_at_z_req, sizeof(double) * pnlpt->k_size, pnlpt->error_message);
        class_alloc(pk_l_at_z_req, sizeof(double) * pnlpt->k_size, pnlpt->error_message);

        /* Retrieve the logarithmic growth rate f for RSD */

        class_alloc(pnlpt->growthf, sizeof(double) * pnlpt->z_pk_num, pnlpt->error_message);
        class_alloc(pnlpt->hratio_array, sizeof(double) * pnlpt->z_pk_num, pnlpt->error_message);
        class_alloc(pnlpt->Dratio_array, sizeof(double) * pnlpt->z_pk_num, pnlpt->error_message);

        if (pnlpt->rsd == rsd_yes) {
            double *pvecbackf;
            int last_indexf;
            class_alloc(pvecbackf, pba->bg_size * sizeof(double), pnlpt->error_message);

            /* Variables for Alcock-Paczynski projection */
            int j;
            int Nz = 2000;
            double Omfid = pnlpt->OmfidAP;
            double dz;
            double kmsMpc = 3.33564095198145e-6;

            if (pnlpt->nonlinear_pt_verbose > 0)
                if (pnlpt->AP_effect == AP_effect_yes) {
                    printf("Computing the Alcock-Paczynski effect for fiducial cosmology with Om=%lf\n", Omfid);
                } else {
                    printf("No Alcock-Paczynski effect.\n");
                }

            /* --- Main redshift loop: compute 1-loop spectra at each z_pk --- */
            for (i_z = 0; i_z < pnlpt->z_pk_num; i_z++) {
                /* Extract background quantities (D, f, H, D_A) at this redshift */
                double Da = 0;
                double Dfid = 0;
                double hnew = 0;
                double hfid = 0;

                class_call(background_at_tau(pba, tau_req[i_z], long_info, inter_normal, &last_indexf, pvecbackf),
                           pnlpt->error_message,
                           pnlpt->error_message);

                /* Override f(z) and D(z) with user-supplied values if replace_background is set */
                if (pnlpt->replace_background) {
                    pnlpt->growthf[i_z] = pnlpt->replace_fz_value;
                    Dref = pnlpt->replace_Dz_value;
                } else {
                    pnlpt->growthf[i_z] = pvecbackf[pba->index_bg_f];
                    Dref = pvecbackf[pba->index_bg_D];
                }

                if (pnlpt->AP_effect == AP_effect_yes) {
                    if (pnlpt->z_pk[i_z] == 0.) {
                        pnlpt->hratio_array[i_z] = 1.;
                        pnlpt->Dratio_array[i_z] = 1.;
                        hnew = 1.;
                    } else {
                        hfid = pow((Omfid * pow((1. + pnlpt->z_pk[i_z]), 3.) + (1. - Omfid) + (pba->Omega0_g) * pow((1. + pnlpt->z_pk[i_z]), 4.)), 0.5);
                        if (pnlpt->replace_background) {
                            hnew = pnlpt->replace_Hz_value / kmsMpc / 100.0 / pba->h;
                        } else {
                            hnew = pvecbackf[pba->index_bg_H] / kmsMpc / 100.0 / pba->h;
                        }
                        pnlpt->hratio_array[i_z] = hnew / hfid;

                        dz = pnlpt->z_pk[i_z] / (double)(Nz - 1);

                        for (j = 1; j < Nz; j++) {
                            Dfid = Dfid + dz * (1. / pow((Omfid * pow((1. + dz * j), 3.) + (1. - Omfid) + (pba->Omega0_g) * pow((1. + pnlpt->z_pk[i_z]), 4.)), 0.5) + 1. / pow((Omfid * pow((1. + dz * (j - 1)), 3.) + (1. - Omfid) + (pba->Omega0_g) * pow((1. + pnlpt->z_pk[i_z]), 4.)), 0.5)) / 2.;
                        }
                        if (pnlpt->replace_background) {
                            Da = pnlpt->replace_DAz_value * kmsMpc * 100 * pba->h * (1. + pnlpt->z_pk[i_z]);
                        } else {
                            Da = pvecbackf[pba->index_bg_ang_distance] * kmsMpc * 100 * pba->h * (1. + pnlpt->z_pk[i_z]);
                        }

                        pnlpt->Dratio_array[i_z] = Da / Dfid;
                    }
                } else {
                    pnlpt->hratio_array[i_z] = 1.;
                    pnlpt->Dratio_array[i_z] = 1.;
                }
            }

            free(pvecbackf);
        } else {
            /* No RSD: set growth/AP ratios to unity, retrieve D(z) from background */
            double *pvecbackf;
            int last_indexf = 0;
            class_alloc(pvecbackf, pba->bg_size * sizeof(double), pnlpt->error_message);

            for (i_z = 0; i_z < pnlpt->z_pk_num; i_z++) {
                pnlpt->growthf[i_z] = 1.;
                pnlpt->hratio_array[i_z] = 1.;
                pnlpt->Dratio_array[i_z] = 1.;
                class_call(background_at_tau(pba, tau_req[i_z], long_info, inter_normal, &last_indexf, pvecbackf),
                           pnlpt->error_message,
                           pnlpt->error_message);
                Dref = pvecbackf[pba->index_bg_D];
            }
            free(pvecbackf);
        }

        /* End RSD/AP parameter setup */

        /* Interpolate P_L(k) at each requested redshift */

        last_index = 0;
        for (i_z = 0; i_z < pnlpt->z_pk_num; i_z++) {
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

            for (index_k = 0; index_k < pnlpt->k_size; index_k++) {
                ln_pk_l_at_z_req[index_k] = ln_pk_l_at_z_req[index_k];
                pk_l_at_z_req[index_k] = exp(ln_pk_l_at_z_req[index_k]);

                ln_pPRIMk_l_req[index_k] = lnpPRIMk_l[index_k];
                pPRIMk_l_req[index_k] = exp(ln_pPRIMk_l_req[index_k]);
            }

            /* get P_NL(k) at tau_req */

            if (print_warning == _FALSE_) {
                /* Call the main 1-loop computation engine */
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
                                             pPRIMk_l_req,
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
                                             pk_nl_fNL,
                                             pk_fNLd2,
                                             pk_fNLG2,
                                             pk_nl_fNL_ortho,
                                             pk_fNLd2_ortho,
                                             pk_fNLG2_ortho,
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
                                             lnk_l,
                                             ln_pk_l_at_z_req,
                                             ln_pPRIMk_l_req),
                           pnlpt->error_message,
                           pnlpt->error_message);
                /* --- Store results: copy 1-loop spectra into pnlpt arrays --- */
                {
                    double *store_pk_arrs_[N_LN_PK] = {
                        pk_nl, pk_Id2d2, pk_Id2d2_2, pk_Id2d2_4,
                        pk_Id2, pk_IG2, pk_Id2G2, pk_Id2G2_2, pk_Id2G2_4,
                        pk_IG2G2, pk_IG2G2_2, pk_IG2G2_4,
                        pk_IFG2, pk_IFG2_0b1, pk_IFG2_0, pk_IFG2_2,
                        pk_CTR, pk_CTR_0, pk_CTR_2, pk_CTR_4,
                        pk_Tree, pk_Tree_0_vv, pk_Tree_0_vd, pk_Tree_0_dd,
                        pk_Tree_2_vv, pk_Tree_2_vd, pk_Tree_4_vv,
                        pk_l_0_vv, pk_l_0_vd, pk_l_0_dd,
                        pk_l_2_vv, pk_l_2_vd, pk_l_2_dd,
                        pk_l_4_vv, pk_l_4_vd, pk_l_4_dd,
                        pk_l_0_b1b2, pk_l_0_b2, pk_l_0_b1bG2, pk_l_0_bG2,
                        pk_l_2_b1b2, pk_l_2_b2, pk_l_2_b1bG2, pk_l_2_bG2,
                        pk_l_4_b2, pk_l_4_bG2, pk_l_4_b1b2, pk_l_4_b1bG2,
                        pk_l_fNL_0_vv, pk_l_fNL_0_vd, pk_l_fNL_0_dd,
                        pk_l_fNL_2_vv, pk_l_fNL_2_vd, pk_l_fNL_2_dd,
                        pk_l_fNL_4_vv, pk_l_fNL_4_vd, pk_l_fNL_4_dd,
                        pk_l_fNL_0_vv_ortho, pk_l_fNL_0_vd_ortho, pk_l_fNL_0_dd_ortho,
                        pk_l_fNL_2_vv_ortho, pk_l_fNL_2_vd_ortho, pk_l_fNL_2_dd_ortho,
                        pk_l_fNL_4_vv_ortho, pk_l_fNL_4_vd_ortho, pk_l_fNL_4_dd_ortho,
                        pk12_l_0_b1b2, pk12_l_0_b2, pk12_l_0_b1bG2, pk12_l_0_bG2,
                        pk12_l_2_b1b2, pk12_l_2_b2, pk12_l_2_b1bG2, pk12_l_2_bG2,
                        pk12_l_4_b1b2, pk12_l_4_b2, pk12_l_4_b1bG2, pk12_l_4_bG2,
                        pk12_l_0_b1b2_ortho, pk12_l_0_b2_ortho, pk12_l_0_b1bG2_ortho, pk12_l_0_bG2_ortho,
                        pk12_l_2_b1b2_ortho, pk12_l_2_b2_ortho, pk12_l_2_b1bG2_ortho, pk12_l_2_bG2_ortho,
                        pk12_l_4_b1b2_ortho, pk12_l_4_b2_ortho, pk12_l_4_b1bG2_ortho, pk12_l_4_bG2_ortho,
                        pk_nl_fNL, pk_fNLd2, pk_fNLG2,
                        pk_nl_fNL_ortho, pk_fNLd2_ortho, pk_fNLG2_ortho
                    };
                    for (index_k = 0; index_k < pnlpt->k_size; index_k++)
                        for (int si_ = 0; si_ < N_LN_PK; si_++) {
                            double *dest_ = *(double **)((char *)pnlpt + ln_pk_offsets[si_]);
                            dest_[i_z * pnlpt->k_size + index_k] = log(store_pk_arrs_[si_][index_k]);
                        }
                }
            }
        }

        /* This is required for lensing with one-loop ! */

        /* --- Lensing module: compute nl_corr_density ratio for CMB lensing --- */
        if (ppt->has_cl_cmb_lensing_potential == _TRUE_) {
            /* Lensing corrections active */
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

            last_index = 0;

            double large_for_logs_matter = 5.e4;

            if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {
                class_alloc(dd_pk_Tree, pnlpt->k_size * sizeof(double), pnlpt->error_message);
                class_alloc(dd_pk_nl, pnlpt->k_size * sizeof(double), pnlpt->error_message);
                class_alloc(dd_pk_l_at_z_req, pnlpt->k_size * sizeof(double), pnlpt->error_message);
                class_alloc(pk_Tree_int, ppt->k_size[pnlpt->index_md_scalars] * sizeof(double), pnlpt->error_message);
                class_alloc(pk_nl_int, ppt->k_size[pnlpt->index_md_scalars] * sizeof(double), pnlpt->error_message);
                class_alloc(pk_l_at_z_req_int, ppt->k_size[pnlpt->index_md_scalars] * sizeof(double), pnlpt->error_message);

                SPLINE_SETUP(pnlpt->k, pnlpt->k_size, pk_Tree, dd_pk_Tree, _SPLINE_EST_DERIV_);
                SPLINE_SETUP(pnlpt->k, pnlpt->k_size, pk_nl, dd_pk_nl, _SPLINE_EST_DERIV_);
                SPLINE_SETUP(pnlpt->k, pnlpt->k_size, pk_l_at_z_req, dd_pk_l_at_z_req, _SPLINE_EST_DERIV_);
                for (index_k = 0; index_k < ppt->k_size[pnlpt->index_md_scalars]; index_k++) {
                    if (ppt->k[pnlpt->index_md_scalars][index_k] < pnlpt->k[pnlpt->k_size - 1]) {
                        SPLINE_EVAL(pnlpt->k, pnlpt->k_size, pk_Tree, dd_pk_Tree, ppt->k[pnlpt->index_md_scalars][index_k], &last_index, &pk_Tree_int[index_k]);
                        SPLINE_EVAL(pnlpt->k, pnlpt->k_size, pk_nl, dd_pk_nl, ppt->k[pnlpt->index_md_scalars][index_k], &last_index, &pk_nl_int[index_k]);
                        SPLINE_EVAL(pnlpt->k, pnlpt->k_size, pk_l_at_z_req, dd_pk_l_at_z_req, ppt->k[pnlpt->index_md_scalars][index_k], &last_index, &pk_l_at_z_req_int[index_k]);
                    } else {
                        pk_Tree_int[index_k] = pk_Tree_int[index_k - 1];
                        pk_nl_int[index_k] = pk_nl_int[index_k - 1];
                        pk_l_at_z_req_int[index_k] = pk_l_at_z_req_int[index_k - 1];
                    }
                }
            }

            for (index_tau = pnlpt->tau_size - 1; index_tau >= 0; index_tau--) {
                class_call(background_at_tau(pba, pnlpt->tau[index_tau], long_info, inter_normal, &last_indexD, pvecbackD), pnlpt->error_message, pnlpt->error_message);

                Dplus = pvecbackD[pba->index_bg_D];
                if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {
                    for (index_k = 0; index_k < ppt->k_size[pnlpt->index_md_scalars]; index_k++)
                        pnlpt->nl_corr_density[index_tau * ppt->k_size[pnlpt->index_md_scalars] + index_k] = sqrt(fabs((pk_Tree_int[index_k] + Dplus * Dplus * (pk_nl_int[index_k] - large_for_logs_matter) /
                                                                                                                                                   Dref / Dref)) /
                                                                                                                  pk_l_at_z_req_int[index_k]);
                } else {
                    for (index_k = 0; index_k < ppt->k_size[pnlpt->index_md_scalars]; index_k++)
                        pnlpt->nl_corr_density[index_tau * ppt->k_size[pnlpt->index_md_scalars] + index_k] = sqrt(fabs((pk_Tree[index_k] + Dplus * Dplus * (pk_nl[index_k] - large_for_logs_matter) /
                                                                                                                                               Dref / Dref)) /
                                                                                                                  pk_l_at_z_req[index_k]);
                }
            }
            free(pvecbackD);

            if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {
                free(dd_pk_Tree);
                free(dd_pk_nl);
                free(dd_pk_l_at_z_req);
                free(pk_Tree_int);
                free(pk_nl_int);
                free(pk_l_at_z_req_int);
            }
        } else {
            /* No CMB lensing: set nl_corr_density = 1 */
            for (index_tau = pnlpt->tau_size - 1; index_tau >= 0; index_tau--) {
                for (index_k = 0; index_k < ppt->k_size[pnlpt->index_md_scalars]; index_k++) {
                    pnlpt->nl_corr_density[index_tau * ppt->k_size[pnlpt->index_md_scalars] + index_k] = 1.;
                }
            }
        }

        /* --- End lensing module --- */

        { void *_f[] = {
            pk_l, pk_l_0_vv, pk_l_0_vd, pk_l_0_dd, pk_l_2_vv, pk_l_2_vd, pk_l_2_dd,
            pk_l_4_vv, pk_l_4_vd, pk_l_4_dd,
            pk_l_0_b1b2, pk_l_0_b2, pk_l_0_b1bG2, pk_l_0_bG2,
            pk_l_2_b1b2, pk_l_2_b2, pk_l_2_b1bG2, pk_l_2_bG2,
            pk_l_4_b2, pk_l_4_bG2, pk_l_4_b1b2, pk_l_4_b1bG2,
            pk_nl, pk_Id2d2, pk_Id2d2_2, pk_Id2d2_4,
            pk_IG2, pk_Id2G2, pk_IG2G2, pk_Id2G2_2, pk_IG2G2_2, pk_Id2G2_4, pk_IG2G2_4,
            pk_IFG2, pk_IFG2_0, pk_IFG2_0b1, pk_IFG2_2, pk_Id2,
            pk_CTR, pk_CTR_0, pk_CTR_2, pk_CTR_4,
            pk_Tree, pk_Tree_0_vv, pk_Tree_0_vd, pk_Tree_0_dd, pk_Tree_2_vv, pk_Tree_2_vd, pk_Tree_4_vv,
            lnk_l, lnpk_l, ddlnpk_l, tau_req,
            lnpk_l_full, pk_l_full, ddlnpk_l_full,
            ln_pk_l_at_z_req, pk_l_at_z_req, ln_pPRIMk_l_req, pPRIMk_l_req,
            pk_l_fNL_0_vv, pk_l_fNL_0_vd, pk_l_fNL_0_dd,
            pk_l_fNL_2_vv, pk_l_fNL_2_vd, pk_l_fNL_2_dd,
            pk_l_fNL_4_vv, pk_l_fNL_4_vd, pk_l_fNL_4_dd,
            pk12_l_0_b1b2, pk12_l_0_b2, pk12_l_0_b1bG2, pk12_l_0_bG2,
            pk12_l_2_b1b2, pk12_l_2_b2, pk12_l_2_b1bG2, pk12_l_2_bG2,
            pk12_l_4_b1b2, pk12_l_4_b2, pk12_l_4_b1bG2, pk12_l_4_bG2,
            pk_nl_fNL, pk_fNLd2, pk_fNLG2,
            pk_l_fNL_0_vv_ortho, pk_l_fNL_0_vd_ortho, pk_l_fNL_0_dd_ortho,
            pk_l_fNL_2_vv_ortho, pk_l_fNL_2_vd_ortho, pk_l_fNL_2_dd_ortho,
            pk_l_fNL_4_vv_ortho, pk_l_fNL_4_vd_ortho, pk_l_fNL_4_dd_ortho,
            pk12_l_0_b1b2_ortho, pk12_l_0_b2_ortho, pk12_l_0_b1bG2_ortho, pk12_l_0_bG2_ortho,
            pk12_l_2_b1b2_ortho, pk12_l_2_b2_ortho, pk12_l_2_b1bG2_ortho, pk12_l_2_bG2_ortho,
            pk12_l_4_b1b2_ortho, pk12_l_4_b2_ortho, pk12_l_4_b1bG2_ortho, pk12_l_4_bG2_ortho,
            pk_nl_fNL_ortho, pk_fNLd2_ortho, pk_fNLG2_ortho,
            pPRIMk_l, lnpPRIMk_l, ddlnpPRIMk_l
        };
        for (int fi_ = 0; fi_ < 111; fi_++) free(_f[fi_]); }

        if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {
            if (pnlpt->cb == _TRUE_) {
                free(pnlpt->dd_sources_tp_delta_cb);
                free(pnlpt->sources_tp_delta_cb);
            } else {
                free(pnlpt->dd_sources_tp_delta_m);
                free(pnlpt->sources_tp_delta_m);
            }
        }

        if (pnlpt->nonlinear_pt_verbose > 0)
            printf(" 'nonlinear_pt_init' module executed successfully\n");
    } else {
        class_stop(pnlpt->error_message,
                   "Your non-linear method variable is set to %d, out of the range defined in nonlinear_pt.h", pnlpt->method);
    }

    return _SUCCESS_;
}

/* === MODULE CLEANUP === */

/**
 * Free all memory allocated by the PT module.
 */
int nonlinear_pt_free(
    struct nonlinear_pt *pnlpt) {
    if (pnlpt->method > nlpt_none) {
        if (pnlpt->method == nlpt_spt) {
            /* Free kernel matrices */

            free(pnlpt->M13_oneline);
            free(pnlpt->M22_oneline);
            free(pnlpt->M22basic_oneline);
            free(pnlpt->IFG2_oneline);
            free(pnlpt->k);
            free(pnlpt->ln_k);
            free(pnlpt->tau);
            free(pnlpt->ln_tau);
            /* Free all complex matrices */
            for (int ci_ = 0; ci_ < N_M22_COMPLEX; ci_++) {
                complex double **ptr = (complex double **)((char *)pnlpt + m22_complex_offsets[ci_]);
                free(*ptr);
            }
            for (int ci_ = 0; ci_ < N_M13_COMPLEX; ci_++) {
                complex double **ptr = (complex double **)((char *)pnlpt + m13_complex_offsets[ci_]);
                free(*ptr);
            }
            free(pnlpt->growthf);
            free(pnlpt->hratio_array);
            free(pnlpt->Dratio_array);
            free(pnlpt->nl_corr_density);
            /* Free all ln_pk output arrays */
            for (int lpi_ = 0; lpi_ < N_LN_PK; lpi_++) {
                double **ptr = (double **)((char *)pnlpt + ln_pk_offsets[lpi_]);
                free(*ptr);
            }
            free(pnlpt->gauss_x);
            free(pnlpt->gauss);
            free(pnlpt->gauss_w);
        }
    }
    return _SUCCESS_;
}

/* === LINEAR POWER SPECTRUM EXTRACTION === */

/**
 * Extract the linear matter power spectrum P_L(k) at a given time.
 *
 * Combines perturbation source functions and primordial spectrum
 * to produce P_L(k), stored as ln(P) with spline derivatives
 * for interpolation.
 */
int nonlinear_pt_pk_l(
    struct background *pba,
    struct perturbations *ppt,
    struct primordial *ppm,
    struct nonlinear_pt *pnlpt,
    int index_tau,
    double *pk_l,
    double *lnk,
    double *lnpk,
    double *ddlnpk) {
    int index_md;
    int index_ic = 0;
    int index_type;
    int index_k;
    int index_ic1, index_ic2, index_ic1_ic2;
    double *primordial_pk;
    double source_ic1, source_ic2;

    /* Replace linear spectrum with a given input file */
    if (pnlpt->replace_pk) {
        /* open input interpolation file */
        char line[100000];
        FILE *fp;
        errno = 0;
        fp = fopen(pnlpt->input_pk, "r");

        if (fp == NULL) {
            fprintf(stderr, "Interpolation file %s not found\n", pnlpt->input_pk);
            fprintf(stderr, "Trying again in 3sec\n");
            sleep(3);
            fp = fopen(pnlpt->input_pk, "r");
            if (fp == NULL) {
                fprintf(stderr, "Interpolation file %s still not found\n Aborting!\n", pnlpt->input_pk);
                fprintf(stderr, "The error code is %d.\n", errno);
                fprintf(stderr, "Error opening file: %s\n", strerror(errno));
                exit(0);
            }
        }
        /* Count lines to construct the correct size */
        int nline = 0;
        while (fgets(line, 10000, fp) != NULL) {
            if (line[0] == '#')
                continue; /* skip comment lines starting with '#' */
            if (line[0] == '\n')
                continue;
            nline++;
        }
        rewind(fp);

        /* Allocate memory for k and P(k) arrays */
        double kint[nline], pkint[nline], kmin, kmax;

        int line_count = 0;
        int counter = 0;

        /* Read in values from file */
        while (fgets(line, 100000, fp) != NULL) {
            /* Skip comment and blank lines */
            if (line[0] == '#')
                continue;
            if (line[0] == '\n')
                continue;

            /* Split into tab-separated variables */
            char *split_string;
            split_string = strtok(line, "\t");
            counter = 0;

            /* Iterate over columns in line */
            while (split_string != NULL) {
                if (counter == 0) {
                    kint[line_count] = atof(split_string);
                }
                if (counter == 1) {
                    pkint[line_count] = atof(split_string);
                }
                if (counter > 1) {
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

        for (index_k = 0; index_k < pnlpt->k_size; index_k++) {
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
                for (segment = 0; segment < nline - 1; segment++) {
                    if ((kint[segment + 1] >= this_k) && (kint[segment] <= this_k)) {
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
    /* Compute linear power from CLASS as usual */
    else
    {
        index_md = ppt->index_md_scalars;

        class_alloc(primordial_pk, ppm->ic_ic_size[index_md] * sizeof(double), pnlpt->error_message);

        for (index_k = 0; index_k < pnlpt->k_size; index_k++) {
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
                index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic1, ppm->ic_size[index_md]);

                if (ppt->has_cls == _TRUE_ && pnlpt->fast_output == _TRUE_) {
                    if (pnlpt->cb == _TRUE_) {
                        source_ic1 = pnlpt->sources_tp_delta_cb[index_tau * pnlpt->k_size + index_k];
                    } else {
                        source_ic1 = pnlpt->sources_tp_delta_m[index_tau * pnlpt->k_size + index_k];
                    }
                } else {
                    if (pnlpt->cb == _TRUE_) {
                        source_ic1 = ppt->sources[index_md][index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_cb][index_tau * ppt->k_size[index_md] + index_k];
                    } else {
                        source_ic1 = ppt->sources[index_md][index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m][index_tau * ppt->k_size[index_md] + index_k];
                    }
                }

                /* source_ic1 and source_ic2 are the transfer functions */

                pk_l[index_k] += 2. * _PI_ * _PI_ / pow(pnlpt->k[index_k], 3) * source_ic1 * source_ic1 * primordial_pk[index_ic1_ic2];
            }

            /* part non-diagonal in initial conditions */
            for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
                for (index_ic2 = index_ic1 + 1; index_ic2 < ppm->ic_size[index_md]; index_ic2++) {
                    index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic2, ppm->ic_size[index_md]);
                    if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {
                        source_ic1 = ppt->sources[index_md]
                                                 [index_ic1 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
                                                 [index_tau * ppt->k_size[index_md] + index_k];
                        source_ic2 = ppt->sources[index_md]
                                                 [index_ic2 * ppt->tp_size[index_md] + ppt->index_tp_delta_m]
                                                 [index_tau * ppt->k_size[index_md] + index_k];
                        pk_l[index_k] += 2. * 2. * _PI_ * _PI_ / pow(pnlpt->k[index_k], 3) * source_ic1 * source_ic2 * primordial_pk[index_ic1_ic2]; /* extra factor of 2 for the symmetric term (ic2,ic1) */
                    }
                }
            }

            lnk[index_k] = log(pnlpt->k[index_k]);
            lnpk[index_k] = log(pk_l[index_k]);
        }
        free(primordial_pk);
    }

    SPLINE_SETUP(lnk, pnlpt->k_size, lnpk, ddlnpk, _SPLINE_NATURAL_);

    return _SUCCESS_;

} /* end nonlinear_pt_pk_l */

/* === PRIMORDIAL POWER SPECTRUM EXTRACTION === */

/**
 * Extract the primordial power spectrum P_prim(k).
 *
 * Used for fNL calculations: the transfer function T(k) is obtained
 * by dividing P_L(k) by P_prim(k), since P_L = T^2 * P_prim.
 */
int nonlinear_pt_pPRIMk_l(
    struct perturbations *ppt,
    struct primordial *ppm,
    struct nonlinear_pt *pnlpt,
    double *lnk,
    double *pPRIMk_l,
    double *lnpPRIMk,
    double *ddlnpPRIMk) {
    int index_md;
    int index_ic = 0;
    int index_type;
    int index_k;
    int index_ic1, index_ic2, index_ic1_ic2;
    double *primordial_pk;

    index_md = ppt->index_md_scalars;

    class_alloc(primordial_pk, ppm->ic_ic_size[index_md] * sizeof(double), pnlpt->error_message);

    for (index_k = 0; index_k < pnlpt->k_size; index_k++) {
        class_call(primordial_spectrum_at_k(ppm,
                                            index_md,
                                            linear,
                                            pnlpt->k[index_k],
                                            primordial_pk),
                   ppm->error_message,
                   pnlpt->error_message);

        pPRIMk_l[index_k] = 0;

        /* part diagonal in initial conditions */
        for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
            index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic1, ppm->ic_size[index_md]);

            /* source_ic1 and source_ic2 are the transfer functions */

            pPRIMk_l[index_k] += 2. * _PI_ * _PI_ * (pnlpt->k[index_k]) * primordial_pk[index_ic1_ic2] / ppm->A_s;
        }

        /* part non-diagonal in initial conditions */
        for (index_ic1 = 0; index_ic1 < ppm->ic_size[index_md]; index_ic1++) {
            for (index_ic2 = index_ic1 + 1; index_ic2 < ppm->ic_size[index_md]; index_ic2++) {
                index_ic1_ic2 = index_symmetric_matrix(index_ic1, index_ic2, ppm->ic_size[index_md]);
                if (ppm->is_non_zero[index_md][index_ic1_ic2] == _TRUE_) {
                    pPRIMk_l[index_k] += 2. * 2. * _PI_ * _PI_ * (pnlpt->k[index_k]) * primordial_pk[index_ic1_ic2] / ppm->A_s;
                }
            }
        }

        lnk[index_k] = log(pnlpt->k[index_k]);

        lnpPRIMk[index_k] = log(pPRIMk_l[index_k]);
    }

    free(primordial_pk);

    SPLINE_SETUP(lnk, pnlpt->k_size, lnpPRIMk, ddlnpPRIMk, _SPLINE_NATURAL_);

    return _SUCCESS_;

} /* end nonlinear_pt_pk_l */

/* ================================================================== */
/* === MAIN ONE-LOOP COMPUTATION ENGINE === */
/* ================================================================== */

/**
 * Compute all one-loop power spectra at a single redshift.
 *
 * This is the core computational routine. Given the linear P(k) and
 * growth rate f at one redshift, it computes:
 *   1. FFTLog decomposition of P_L(k) into power-law basis
 *   2. P22 and P13 one-loop integrals via precomputed kernel matrices
 *   3. IR resummation of BAO wiggles (Senatore & Zaldarriaga)
 *   4. RSD multipoles (ell=0,2,4) for velocity-velocity, velocity-density,
 *      and density-density components
 *   5. Alcock-Paczynski projection via Gauss-Legendre quadrature
 *   6. Bias operator integrals (Id2d2, Id2G2, IG2G2, IFG2)
 *   7. fNL contributions (local, equilateral, orthogonal)
 */

/* ==================================================================
 * BATCH BLAS-3 HELPERS FOR REAL-VALUED PT COMPUTATIONS
 *
 * The FFTLog decomposition represents P(k) as a sum of power laws:
 *   P(k) = sum_m c_m * k^{eta_m}
 * where c_m are complex coefficients and eta_m = b + 2*pi*i*m/(N*Delta).
 *
 * The one-loop integrals become matrix operations on X[c,j]:
 *   X[c,j] = c_m * k_j^{eta_m}   (the "X matrix")
 *   P13[j] = X[j,:]^T * M13      (dot product, dgemv)
 *   P22[j] = X[j,:]^T * M22 * X[j,:]  (quadratic form, dsymm+ddot)
 *   P12[j] = Xl[j,:]^T * M12 * Xr[j,:]  (bilinear form)
 *
 * All operations use real BLAS (dgemv, dsymm, ddot) by separating
 * real and imaginary parts, avoiding complex BLAS dependencies.
 * ================================================================== */

/** Build X matrix: X[c,j] = coeffs[c] * k[j]^exponents[c] for all j.
 *  Layout: column-major, X[c + j*ncoeff] for coefficient c at wavenumber j.
 *
 *  Optimization: exponents have the form b + i*omega[c] where b is a real
 *  constant (the FFTLog bias). So k^(b+i*omega) = k^b * exp(i*omega*ln(k)).
 *  kb_cache and lnk_cache are pre-computed arrays of k^b and ln(k);
 *  if NULL, they are computed internally (fallback). */
static void build_X_matrix(
    int nk, int ncoeff,
    const double complex *coeffs,
    const double complex *exponents,
    const double *k,
    double *Xr, double *Xi) {
    double b_real = creal(exponents[0]);

    for (int j = 0; j < nk; j++) {
        double lnk = log(k[j]);
        double kb = pow(k[j], b_real);
        int offset = j * ncoeff;
        for (int c = 0; c < ncoeff; c++) {
            double omega_c = cimag(exponents[c]);
            double s, co;
            sincos(omega_c * lnk, &s, &co);
            double cr_c = creal(coeffs[c]);
            double ci_c = cimag(coeffs[c]);
            Xr[c + offset] = kb * (cr_c * co - ci_c * s);
            Xi[c + offset] = kb * (cr_c * s + ci_c * co);
        }
    }
}

/** Build X matrix with cached k^b and ln(k) arrays to avoid recomputing pow/log. */
static void build_X_matrix_cached(
    int nk, int ncoeff,
    const double complex *coeffs,
    const double complex *exponents,
    const double *kb_cache,
    const double *lnk_cache,
    double *Xr, double *Xi) {
    for (int j = 0; j < nk; j++) {
        double lnk = lnk_cache[j];
        double kb = kb_cache[j];
        int offset = j * ncoeff;
        for (int c = 0; c < ncoeff; c++) {
            double omega_c = cimag(exponents[c]);
            double s, co;
            sincos(omega_c * lnk, &s, &co);
            double cr_c = creal(coeffs[c]);
            double ci_c = cimag(coeffs[c]);
            Xr[c + offset] = kb * (cr_c * co - ci_c * s);
            Xi[c + offset] = kb * (cr_c * s + ci_c * co);
        }
    }
}

/** Precompute sincos table for given exponents and ln(k) values.
 *  Layout: sc_cos[c + j*ncoeff], sc_sin[c + j*ncoeff].
 *  Caller must allocate sc_cos, sc_sin of size nk*ncoeff. */
static void precompute_sincos_table(
    int nk, int ncoeff,
    const double complex *exponents,
    const double *lnk_cache,
    double *sc_cos, double *sc_sin) {
    for (int j = 0; j < nk; j++) {
        double lnk = lnk_cache[j];
        int offset = j * ncoeff;
        for (int c = 0; c < ncoeff; c++) {
            double omega_c = cimag(exponents[c]);
            sincos(omega_c * lnk, &sc_sin[c + offset], &sc_cos[c + offset]);
        }
    }
}

/** Build X matrix from precomputed sincos table (no trig calls). */
static void build_X_from_tables(
    int nk, int ncoeff,
    const double complex *coeffs,
    const double *kb_cache,
    const double *sc_cos, const double *sc_sin,
    double *Xr, double *Xi) {
    /* Only parallelize if we have enough work per thread (nk=128, reasonable up to 8 threads) */
    #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
    for (int j = 0; j < nk; j++) {
        double kb = kb_cache[j];
        int offset = j * ncoeff;
        for (int c = 0; c < ncoeff; c++) {
            double co = sc_cos[c + offset];
            double s = sc_sin[c + offset];
            double cr_c = creal(coeffs[c]);
            double ci_c = cimag(coeffs[c]);
            Xr[c + offset] = kb * (cr_c * co - ci_c * s);
            Xi[c + offset] = kb * (cr_c * s + ci_c * co);
        }
    }
}

/** Batch dot product: f[j] = X[j,:]^T · m  for all j (using real BLAS)
 *  mr_ws, mi_ws: workspace arrays of size ncoeff (caller-allocated) */
static void batch_dot_product(
    int ncoeff, int nk,
    const double *Xr, const double *Xi,
    const double complex *m_complex,
    double *f,
    double *mr_ws, double *mi_ws) {
    for (int c = 0; c < ncoeff; c++) {
        mr_ws[c] = creal(m_complex[c]);
        mi_ws[c] = cimag(m_complex[c]);
    }

    double one = 1., zero = 0., neg = -1.;
    int inc = 1;
    dgemv_("T", &ncoeff, &nk, &one, (double *)Xr, &ncoeff, mr_ws, &inc, &zero, f, &inc);
    dgemv_("T", &ncoeff, &nk, &neg, (double *)Xi, &ncoeff, mi_ws, &inc, &one, f, &inc);
}

/** Batch quadratic form: f[j] = X[j,:]^T · M · X[j,:]  for all j */
static void batch_quadratic_form(
    int ncoeff, int nk,
    const double complex *M_packed,
    const double *Xr, const double *Xi,
    double *f,
    double *Mr_ws, double *Mi_ws,
    double *Y1_ws, double *Y2_ws) {
    /* Unpack complex symmetric M (upper-triangular packed order) into real Mr and Mi */
    int k = 0;
    for (int j = 0; j < ncoeff; j++) {
        for (int i = j; i < ncoeff; i++) {
            Mr_ws[j + i * ncoeff] = Mr_ws[i + j * ncoeff] = creal(M_packed[k]);
            Mi_ws[j + i * ncoeff] = Mi_ws[i + j * ncoeff] = cimag(M_packed[k]);
            k++;
        }
    }

    /* f = diag(Xr^T Mr Xr - Xi^T Mr Xi - 2 Xr^T Mi Xi) via dsymm */
    double one = 1., zero = 0.;

    /* Y1 = Mr * Xr, Y2 = Mr * Xi */
    dsymm_("L", "U", &ncoeff, &nk, &one, Mr_ws, &ncoeff, (double *)Xr, &ncoeff, &zero, Y1_ws, &ncoeff);
    dsymm_("L", "U", &ncoeff, &nk, &one, Mr_ws, &ncoeff, (double *)Xi, &ncoeff, &zero, Y2_ws, &ncoeff);

    /* Extract diagonal: f[j] = Xr[:,j].Y1[:,j] - Xi[:,j].Y2[:,j] */
    for (int j = 0; j < nk; j++) {
        const double *xr_j = &Xr[j * ncoeff];
        const double *xi_j = &Xi[j * ncoeff];
        const double *y1_j = &Y1_ws[j * ncoeff];
        const double *y2_j = &Y2_ws[j * ncoeff];
        double s = 0.;
        for (int c = 0; c < ncoeff; c++)
            s += xr_j[c] * y1_j[c] - xi_j[c] * y2_j[c];
        f[j] = s;
    }

    /* Y1 = Mi * Xi */
    dsymm_("L", "U", &ncoeff, &nk, &one, Mi_ws, &ncoeff, (double *)Xi, &ncoeff, &zero, Y1_ws, &ncoeff);

    /* f[j] -= 2 * Xr[:,j].Y1[:,j] */
    for (int j = 0; j < nk; j++) {
        const double *xr_j = &Xr[j * ncoeff];
        const double *y1_j = &Y1_ws[j * ncoeff];
        double s = 0.;
        for (int c = 0; c < ncoeff; c++)
            s += xr_j[c] * y1_j[c];
        f[j] -= 2. * s;
    }
}

/** Batch bilinear form: f[j] = Re[ Xl[j,:]^H · M · Xr[j,:] ] for all j
 *  where Xl = Xlr + i*Xli, Xr = Xrr + i*Xri, M = Mr + i*Mi (symmetric).
 *  Uses 4 dsymm calls: Mr*Xrr, Mi*Xrr, Mr*Xri, Mi*Xri. */
static void batch_bilinear_form(
    int ncoeff, int nk,
    const double complex *M_packed,
    const double *Xlr, const double *Xli,
    const double *Xrr, const double *Xri,
    double *f,
    double *Mr_ws, double *Mi_ws, double *Y_ws) {
    /* Unpack M (upper-triangular packed order) */
    int k = 0;
    for (int j = 0; j < ncoeff; j++) {
        for (int i = j; i < ncoeff; i++) {
            Mr_ws[j + i * ncoeff] = Mr_ws[i + j * ncoeff] = creal(M_packed[k]);
            Mi_ws[j + i * ncoeff] = Mi_ws[i + j * ncoeff] = cimag(M_packed[k]);
            k++;
        }
    }

    double one = 1., zero = 0.;

    /* Y = Mr * Xrr;  f[j] = Xlr[:,j] . Y[:,j] */
    dsymm_("L", "U", &ncoeff, &nk, &one, Mr_ws, &ncoeff, (double *)Xrr, &ncoeff, &zero, Y_ws, &ncoeff);
    for (int j = 0; j < nk; j++) {
        double s = 0.;
        for (int c = 0; c < ncoeff; c++)
            s += Xlr[j * ncoeff + c] * Y_ws[j * ncoeff + c];
        f[j] = s;
    }

    /* Y = Mi * Xrr;  f[j] -= Xli[:,j] . Y[:,j] */
    dsymm_("L", "U", &ncoeff, &nk, &one, Mi_ws, &ncoeff, (double *)Xrr, &ncoeff, &zero, Y_ws, &ncoeff);
    for (int j = 0; j < nk; j++) {
        double s = 0.;
        for (int c = 0; c < ncoeff; c++)
            s += Xli[j * ncoeff + c] * Y_ws[j * ncoeff + c];
        f[j] -= s;
    }

    /* Y = Mi * Xri;  f[j] -= Xlr[:,j] . Y[:,j] */
    dsymm_("L", "U", &ncoeff, &nk, &one, Mi_ws, &ncoeff, (double *)Xri, &ncoeff, &zero, Y_ws, &ncoeff);
    for (int j = 0; j < nk; j++) {
        double s = 0.;
        for (int c = 0; c < ncoeff; c++)
            s += Xlr[j * ncoeff + c] * Y_ws[j * ncoeff + c];
        f[j] -= s;
    }

    /* Y = Mr * Xri;  f[j] -= Xli[:,j] . Y[:,j] */
    dsymm_("L", "U", &ncoeff, &nk, &one, Mr_ws, &ncoeff, (double *)Xri, &ncoeff, &zero, Y_ws, &ncoeff);
    for (int j = 0; j < nk; j++) {
        double s = 0.;
        for (int c = 0; c < ncoeff; c++)
            s += Xli[j * ncoeff + c] * Y_ws[j * ncoeff + c];
        f[j] -= s;
    }
}

/* Assembly helpers: combine FFTLog coefficients with UV counterterms and
 * apply high-k exponential cutoff exp(-(k/cutoff)^6) to regularize. */

/** Assemble P13: P13[j] = (k^3 * f * P_L + P13UV) * exp(-cutoff) */
static void assemble_P13(int nk, const double *k, const double *f,
                         const double *Pbin, const double *P13UV,
                         double cutoff, double *P13) {
    double inv_cutoff = 1.0 / cutoff;
    for (int j = 0; j < nk; j++) {
        double kc = k[j] * inv_cutoff; double kc2 = kc * kc; double kc6 = kc2 * kc2 * kc2;
        P13[j] = (k[j] * k[j] * k[j] * f[j] * Pbin[j] + P13UV[j]) * exp(-kc6);
    }
}

/** Assemble P22: P22[j] = k^3 * f * exp(-cutoff) */
static void assemble_P22(int nk, const double *k, const double *f,
                         double cutoff, double *P22) {
    double inv_cutoff = 1.0 / cutoff;
    for (int j = 0; j < nk; j++) {
        double kc = k[j] * inv_cutoff; double kc2 = kc * kc; double kc6 = kc2 * kc2 * kc2;
        P22[j] = k[j] * k[j] * k[j] * f[j] * exp(-kc6);
    }
}

/** Assemble P12: P12[j] = T(k) * k^3 * f * exp(-cutoff) (fNL transfer) */
static void assemble_P12(int nk, const double *k, const double *f,
                         const double *Tbin, double cutoff, double *P12) {
    double inv_cutoff = 1.0 / cutoff;
    for (int j = 0; j < nk; j++) {
        double kc = k[j] * inv_cutoff; double kc2 = kc * kc; double kc6 = kc2 * kc2 * kc2;
        P12[j] = Tbin[j] * k[j] * k[j] * k[j] * f[j] * exp(-kc6);
    }
}

/* ==================================================================
 * MACROS FOR REPETITIVE PT COMPUTATION PATTERNS
 *
 * These macros eliminate ~200+ near-identical code blocks that appear
 * in the RSD, IR resummation, and bias computation sections. Each macro
 * encapsulates a common pattern involving batch BLAS operations, spline
 * interpolation, or Legendre projection.
 *
 * MACRO REFERENCE:
 *   COMPUTE_P13(suffix, M13, uv)     - P13 one-loop via batch dot product
 *   COMPUTE_P22(suffix, M22)         - P22 one-loop via batch quadratic form
 *   COMPUTE_P12(suffix, M12)         - P12 fNL correction (if has_fNL==1)
 *   SPLINE_SETUP(x,n,data,dd,mode)   - Compute spline derivative table
 *   SPLINE_EVAL(x,n,data,dd,xv,li,o) - Evaluate spline at single point
 *   AP_SPLINE_SETUP(name)             - Alloc + spline for AP projection
 *   AP_INTERP(name)                   - Interpolate AP-rescaled k
 *   AP_SPLINE_FREE(name)              - Free AP spline derivative array
 *   AP_SPLINE_SETUP_EX(label, arr)    - Extended AP setup (decoupled names)
 *   AP_INTERP_EX(label, arr)          - Extended AP interpolation
 *   COMPUTE_NW_P22_P12(...)           - IR no-wiggle P22+P12 assembly
 *   COMPUTE_NW_P22_ONLY(...)          - IR no-wiggle P22 only
 *   COMPUTE_W_P22_P12(...)            - IR wiggle P22+P12 assembly
 *   COMPUTE_W_P22_ONLY(...)           - IR wiggle P22 only
 *   SPLINE_INTERP_OUTPUT(...)         - Spline + interpolate to output grid
 *   SPLINE_INTERP_BATCH(...)          - Batch version of SPLINE_INTERP_OUTPUT
 *   COMPUTE_P12_BIAS_PAIR(...)        - P12 pair (base+ortho) for bias
 *   LEGENDRE_PROJECT(val, P0,P2,P4)   - Project into ell=0,2,4 multipoles
 *   ASSEMBLE_P22_BIAS(M22, P_out)     - k^3-scaled P22 for bias operators
 *   TIDAL_P22(M, Xr, Xi, P, post)     - Self-contained tidal P22 integral
 *   COMPUTE_FNL_P12_PAIR(...)         - fNL P12 pair (base+ortho)
 *   FILL_M12_BIAS(...)                - Fill M12 matrix + compute P12 pair
 *   FILL_M22_BIAS(...)                - Fill M22 matrix + assemble P22
 *   FFT_TO_CMSYM(P, b, eta, cmsym)    - Full FFTLog decomposition pipeline
 *
 * Free variables assumed in scope (set up in nonlinear_pt_loop):
 *   Nmax, Np1, kdisc, Pbin, Tbin, cutoff, has_fNL
 *   Xr, Xi, Xr_transfer, Xi_transfer (FFTLog X matrices)
 *   Mr_ws, Mi_ws, Y1_ws, Y2_ws (BLAS workspace)
 *   pnlpt->error_message (for class_call/class_alloc)
 * ================================================================== */

/** Compute P13 for a given matrix and UV counterterm formula.
 *  @param suffix   name suffix (e.g. 0_vv) -> P13_0_vv, P13UV_0_vv
 *  @param M13_matrix  complex M13 kernel vector
 *  @param uv_formula  C expression for P13UV[j] (may reference j, Pbin, kdisc, sigmav, f) */
#define COMPUTE_P13(suffix, M13_matrix, uv_formula)                                     \
    do                                                                                  \
    {                                                                                   \
        double *f13_tmp = malloc(Nmax * sizeof(double));                                \
        TIMER_START(blas);                                                               \
        batch_dot_product(Np1, Nmax, Xr, Xi, M13_matrix, f13_tmp, dot_mr_ws, dot_mi_ws);                      \
        TIMER_ADD(blas);                                                                 \
        for (int j = 0; j < Nmax; j++)                                                  \
            P13UV_##suffix[j] = uv_formula;                                             \
        assemble_P13(Nmax, kdisc, f13_tmp, Pbin, P13UV_##suffix, cutoff, P13_##suffix); \
        free(f13_tmp);                                                                  \
    } while (0)

/** Compute P22 one-loop integral via batch quadratic form.
 *  @param suffix      name suffix -> P22_{suffix}
 *  @param M22_matrix  packed complex symmetric M22 kernel matrix */
#define COMPUTE_P22(suffix, M22_matrix)                            \
    do                                                             \
    {                                                              \
        double *f22_tmp = malloc(Nmax * sizeof(double));           \
        TIMER_START(blas);                                          \
        batch_quadratic_form(Np1, Nmax, M22_matrix, Xr, Xi,        \
                             f22_tmp, Mr_ws, Mi_ws, Y1_ws, Y2_ws); \
        TIMER_ADD(blas);                                            \
        assemble_P22(Nmax, kdisc, f22_tmp, cutoff, P22_##suffix);  \
        free(f22_tmp);                                             \
    } while (0)

/** Compute P12 fNL correction via batch quadratic form with transfer X.
 *  Original: f12[j] = x_transfer^T * M12 * x_transfer (quadratic form).
 *  Only executes when has_fNL (fNL enabled).
 *  @param suffix      name suffix -> P12_{suffix}
 *  @param M12_matrix  complex M12 kernel matrix */
#define COMPUTE_P12(suffix, M12_matrix)                                     \
    do                                                                      \
    {                                                                       \
        if (has_fNL)                                              \
        {                                                                   \
            double *f12_tmp = malloc(Nmax * sizeof(double));                \
            TIMER_START(blas);                                               \
            batch_quadratic_form(Np1, Nmax, M12_matrix,                     \
                                 Xr_transfer, Xi_transfer,                  \
                                 f12_tmp, Mr_ws, Mi_ws, Y1_ws, Y2_ws);     \
            TIMER_ADD(blas);                                                 \
            assemble_P12(Nmax, kdisc, f12_tmp, Tbin, cutoff, P12_##suffix); \
            free(f12_tmp);                                                  \
        }                                                                   \
    } while (0)

/* AP spline setup: allocate dd array, spline, declare output variable */
#define AP_SPLINE_SETUP(name)                                                      \
    double *dd_##name;                                                             \
    class_alloc(dd_##name, sizeof(double) * Nmax, pnlpt->error_message);           \
    class_call(array_spline_table_columns(kdisc, Nmax, name, 1, dd_##name,         \
                                          _SPLINE_NATURAL_, pnlpt->error_message), \
               pnlpt->error_message, pnlpt->error_message);                        \
    double name##_ap_out = 0

/* AP spline interpolation at ktrue (slow version — individual binary search) */
#define AP_INTERP(name)                                                                               \
    class_call(array_interpolate_spline(kdisc, Nmax, name, dd_##name, 1,                              \
                                        ktrue, &last_index, &name##_ap_out, 1, pnlpt->error_message), \
               pnlpt->error_message, pnlpt->error_message)

/* Fast AP interpolation: reuses precomputed binary search result.
 * Before using, compute ap_inf/ap_a/ap_b/ap_a3a/ap_b3b/ap_h2_6 once per ktrue. */
#define AP_INTERP_FAST(name) \
    name##_ap_out = ap_a_ * name[ap_inf_] + ap_b_ * name[ap_sup_] \
                  + (ap_a3a_ * dd_##name[ap_inf_] + ap_b3b_ * dd_##name[ap_sup_]) * ap_h2_6_

/* Fast version for extended labels */
#define AP_INTERP_EX_FAST(label, arr) \
    label##_ap_out = ap_a_ * arr[ap_inf_] + ap_b_ * arr[ap_sup_] \
                   + (ap_a3a_ * dd_##label[ap_inf_] + ap_b3b_ * dd_##label[ap_sup_]) * ap_h2_6_

/* Binary search on kdisc to find interval for ktrue, then precompute cubic spline coefficients.
 * Must be called once before AP_INTERP_FAST calls. */
#define AP_BSEARCH_SETUP()                                                   \
    int ap_inf_ = 0, ap_sup_ = Nmax - 1;                                    \
    while (ap_sup_ - ap_inf_ > 1) {                                         \
        int ap_mid_ = (ap_inf_ + ap_sup_) >> 1;                             \
        if (ktrue < kdisc[ap_mid_]) ap_sup_ = ap_mid_; else ap_inf_ = ap_mid_; \
    }                                                                        \
    double ap_h_ = kdisc[ap_sup_] - kdisc[ap_inf_];                         \
    double ap_b_ = (ktrue - kdisc[ap_inf_]) / ap_h_;                        \
    double ap_a_ = 1. - ap_b_;                                              \
    double ap_a3a_ = ap_a_ * ap_a_ * ap_a_ - ap_a_;                        \
    double ap_b3b_ = ap_b_ * ap_b_ * ap_b_ - ap_b_;                        \
    double ap_h2_6_ = ap_h_ * ap_h_ / 6.

/* AP spline cleanup */
#define AP_SPLINE_FREE(name) free(dd_##name)

/* Extended AP macros: label controls dd_/ap_out names, arr is the data array */
#define AP_SPLINE_SETUP_EX(label, arr)                                             \
    double *dd_##label;                                                            \
    class_alloc(dd_##label, sizeof(double) * Nmax, pnlpt->error_message);          \
    class_call(array_spline_table_columns(kdisc, Nmax, arr, 1, dd_##label,         \
                                          _SPLINE_NATURAL_, pnlpt->error_message), \
               pnlpt->error_message, pnlpt->error_message);                        \
    double label##_ap_out = 0

#define AP_INTERP_EX(label, arr)                                                                       \
    class_call(array_interpolate_spline(kdisc, Nmax, arr, dd_##label, 1,                               \
                                        ktrue, &last_index, &label##_ap_out, 1, pnlpt->error_message), \
               pnlpt->error_message, pnlpt->error_message)


/* IR resummation: NW P22 + P12 (with ortho) for a given mu-power */

#define COMPUTE_NW_P22_P12(suffix, M22_mat, M12_mat, M12_ortho_mat)                    \
    batch_quadratic_form(Np1, Nmax, M22_mat, Xr_nw, Xi_nw,                             \
                         f22_tmp, Mr_ws, Mi_ws, Y1_ws, Y2_ws);                         \
    assemble_P22(Nmax, kdisc, f22_tmp, cutoff, P22_##suffix);                          \
    if (has_fNL)                                                             \
    {                                                                                  \
        batch_quadratic_form(Np1, Nmax, M12_mat, Xr_nw_transfer, Xi_nw_transfer,       \
                             f12_tmp, Mr_ws, Mi_ws, Y1_ws, Y2_ws);                     \
        assemble_P12(Nmax, kdisc, f12_tmp, Tnw, cutoff, P12_##suffix);                 \
        batch_quadratic_form(Np1, Nmax, M12_ortho_mat, Xr_nw_transfer, Xi_nw_transfer, \
                             f12_tmp, Mr_ws, Mi_ws, Y1_ws, Y2_ws);                     \
        assemble_P12(Nmax, kdisc, f12_tmp, Tnw, cutoff, P12_##suffix##_ortho);         \
    }

/* IR resummation: NW P22 only (no P12) */
#define COMPUTE_NW_P22_ONLY(suffix, M22_mat)                   \
    batch_quadratic_form(Np1, Nmax, M22_mat, Xr_nw, Xi_nw,     \
                         f22_tmp, Mr_ws, Mi_ws, Y1_ws, Y2_ws); \
    assemble_P22(Nmax, kdisc, f22_tmp, cutoff, P22_##suffix)

/* IR resummation: wiggle P22 + P12 (with ortho) for a given mu-power */
#define COMPUTE_W_P22_P12(suffix, M22_mat, M12_mat, M12_ortho_mat)                           \
    batch_bilinear_form(Np1, Nmax, M22_mat, Xr_nw, Xi_nw, Xr_2w, Xi_2w,                      \
                        f22_w_tmp, Mr_ws, Mi_ws, Y1_ws);                                     \
    assemble_P22(Nmax, kdisc, f22_w_tmp, cutoff, P22_##suffix##_w);                          \
    if (has_fNL)                                                                   \
    {                                                                                        \
        batch_bilinear_form(Np1, Nmax, M12_mat, Xr_nw_transfer, Xi_nw_transfer,              \
                            Xr_2w_transfer, Xi_2w_transfer, f12_w_tmp, Mr_ws, Mi_ws, Y1_ws); \
        assemble_P12(Nmax, kdisc, f12_w_tmp, Tnw, cutoff, P12_##suffix##_w);                 \
        batch_bilinear_form(Np1, Nmax, M12_ortho_mat, Xr_nw_transfer, Xi_nw_transfer,        \
                            Xr_2w_transfer, Xi_2w_transfer, f12_w_tmp, Mr_ws, Mi_ws, Y1_ws); \
        assemble_P12(Nmax, kdisc, f12_w_tmp, Tnw, cutoff, P12_##suffix##_w_ortho);           \
    }

/* IR resummation: wiggle P22 only (no P12) */
#define COMPUTE_W_P22_ONLY(suffix, M22_mat)                             \
    batch_bilinear_form(Np1, Nmax, M22_mat, Xr_nw, Xi_nw, Xr_2w, Xi_2w, \
                        f22_w_tmp, Mr_ws, Mi_ws, Y1_ws);                \
    assemble_P22(Nmax, kdisc, f22_w_tmp, cutoff, P22_##suffix##_w)

/* Output spline interpolation: transfer results from the FFTLog grid (kdisc[])
 * to the output k-grid (pnlpt->k[]). Inside [kmin_b, kmax_b], values are
 * spline-interpolated; outside, the extrapolation expression is used.
 *
 * @param P_in        input array on kdisc grid (size Nmax)
 * @param pk_out      output array on pnlpt->k grid (size pnlpt->k_size)
 * @param kmin_b      lower bound of valid interpolation range
 * @param kmax_b      upper bound of valid interpolation range
 * @param interp_expr C expression for interpolated value (may use out_tmp_)
 * @param extrap_expr C expression for extrapolated value
 * @param do_free     if nonzero, free(P_in) after interpolation */
#define SPLINE_INTERP_OUTPUT(P_in, pk_out, kmin_b, kmax_b, interp_expr, extrap_expr, do_free)     \
    do                                                                                            \
    {                                                                                             \
        double *dd_sio_;                                                                          \
        double out_tmp_ = 0.;                                                                     \
        class_alloc(dd_sio_, sizeof(double) * Nmax, pnlpt->error_message);                        \
        class_call(array_spline_table_columns(kdisc, Nmax, P_in, 1, dd_sio_,                      \
                                              _SPLINE_NATURAL_, pnlpt->error_message),            \
                   pnlpt->error_message, pnlpt->error_message);                                   \
        last_index = 0;                                                                           \
        for (index_k = 0; index_k < pnlpt->k_size; index_k++)                                     \
        {                                                                                         \
            if (pnlpt->k[index_k] <= (kmax_b) && pnlpt->k[index_k] >= (kmin_b))                   \
            {                                                                                     \
                class_call(array_interpolate_spline(kdisc, Nmax, P_in, dd_sio_, 1,                \
                                                    pnlpt->k[index_k], &last_index, &out_tmp_, 1, \
                                                    pnlpt->error_message),                        \
                           pnlpt->error_message, pnlpt->error_message);                           \
                pk_out[index_k] = (interp_expr);                                                  \
            }                                                                                     \
            else                                                                                  \
            {                                                                                     \
                pk_out[index_k] = (extrap_expr);                                                  \
            }                                                                                     \
        }                                                                                         \
        free(dd_sio_);                                                                            \
        if (do_free)                                                                              \
            free(P_in);                                                                           \
    } while (0)

/* Batch spline-interpolate: applies SPLINE_INTERP_OUTPUT to n_pairs arrays at once.
 * All pairs share the same k-range and interpolation/extrapolation expressions.
 *
 * @param p_in_arr   array of double* pointers (inputs on kdisc grid)
 * @param p_out_arr  array of double* pointers (outputs on pnlpt->k grid)
 * @param n_pairs    number of (input, output) pairs to process
 * @param do_free    if nonzero, free each input array after interpolation */
#define SPLINE_INTERP_BATCH(p_in_arr, p_out_arr, n_pairs, kmin_b, kmax_b,     \
                            interp_expr, extrap_expr, do_free)                 \
    do {                                                                       \
        for (int sib_i_ = 0; sib_i_ < (n_pairs); sib_i_++) {                  \
            double *sib_pin_ = (p_in_arr)[sib_i_];                            \
            double *sib_pout_ = (p_out_arr)[sib_i_];                          \
            double *dd_sio_;                                                   \
            double out_tmp_ = 0.;                                              \
            class_alloc(dd_sio_, sizeof(double) * Nmax, pnlpt->error_message); \
            class_call(array_spline_table_columns(kdisc, Nmax, sib_pin_, 1,    \
                       dd_sio_, _SPLINE_NATURAL_, pnlpt->error_message),       \
                       pnlpt->error_message, pnlpt->error_message);           \
            last_index = 0;                                                    \
            for (index_k = 0; index_k < pnlpt->k_size; index_k++) {           \
                if (pnlpt->k[index_k] <= (kmax_b) &&                          \
                    pnlpt->k[index_k] >= (kmin_b)) {                          \
                    class_call(array_interpolate_spline(kdisc, Nmax,           \
                               sib_pin_, dd_sio_, 1, pnlpt->k[index_k],       \
                               &last_index, &out_tmp_, 1,                      \
                               pnlpt->error_message),                         \
                               pnlpt->error_message, pnlpt->error_message);   \
                    sib_pout_[index_k] = (interp_expr);                        \
                } else {                                                       \
                    sib_pout_[index_k] = (extrap_expr);                        \
                }                                                              \
            }                                                                  \
            free(dd_sio_);                                                     \
            if (do_free) free(sib_pin_);                                       \
        }                                                                      \
    } while (0)

/* Compute P12 pair (base + ortho) from transfer matrices, with k^3*T(k) scaling.
 * Allocates/frees X transfer matrix and f12 workspace internally. */
#define COMPUTE_P12_BIAS_PAIR(M12_base, M12_ortho, P12_out, P12_out_ortho,     \
                              cmsym_t, kb_t_cache, sc_c_tab, sc_s_tab)          \
    do {                                                                        \
        double *Xr_t_, *Xi_t_;                                                 \
        class_alloc(Xr_t_, Np1 * Nmax * sizeof(double), pnlpt->error_message); \
        class_alloc(Xi_t_, Np1 * Nmax * sizeof(double), pnlpt->error_message); \
        TIMER_START(buildx);                                                    \
        build_X_from_tables(Nmax, Np1, cmsym_t, kb_t_cache,                    \
                            sc_c_tab, sc_s_tab, Xr_t_, Xi_t_);                \
        TIMER_ADD(buildx);                                                      \
        double *f12_t_ = malloc(Nmax * sizeof(double));                        \
        TIMER_START(blas);                                                      \
        batch_quadratic_form(Np1, Nmax, M12_base, Xr_t_, Xi_t_,               \
                             f12_t_, Mr_ws, Mi_ws, Y1_ws, Y2_ws);             \
        TIMER_ADD(blas);                                                        \
        for (int j_ = 0; j_ < Nmax; j_++)                                     \
            P12_out[j_] = Tbin[j_] * kdisc[j_] * kdisc[j_] * kdisc[j_] * f12_t_[j_]; \
        TIMER_START(blas);                                                      \
        batch_quadratic_form(Np1, Nmax, M12_ortho, Xr_t_, Xi_t_,              \
                             f12_t_, Mr_ws, Mi_ws, Y1_ws, Y2_ws);             \
        TIMER_ADD(blas);                                                        \
        for (int j_ = 0; j_ < Nmax; j_++)                                     \
            P12_out_ortho[j_] = Tbin[j_] * kdisc[j_] * kdisc[j_] * kdisc[j_] * f12_t_[j_]; \
        free(f12_t_); free(Xr_t_); free(Xi_t_);                               \
    } while (0)

/* Project scalar into Legendre multipoles: ell=0,2,4 */
#define LEGENDRE_PROJECT(val, P0, P2, P4)                                      \
    P0[index_j] += (val) * LegendreP0 / 2.;                                   \
    P2[index_j] += (val) * LegendreP2 * 2.5;                                  \
    P4[index_j] += (val) * LegendreP4 * 4.5

/* Compute P22 from bias matrix and apply k^3 scaling.
 * f22_tmp and Xr2/Xi2 must already be in scope. */
#define ASSEMBLE_P22_BIAS(M22_mat, P22_out)                                    \
    TIMER_START(blas);                                                          \
    batch_quadratic_form(Np1, Nmax, M22_mat, Xr2, Xi2,                        \
                         f22_tmp, Mr_ws, Mi_ws, Y1_ws, Y2_ws);                \
    TIMER_ADD(blas);                                                            \
    for (int j_ = 0; j_ < Nmax; j_++)                                         \
        P22_out[j_] = kdisc[j_] * kdisc[j_] * kdisc[j_] * f22_tmp[j_]

/* Self-contained tidal P22: alloc workspace, quadratic form, k^3 scale, free.
 * post_fn is applied to k^3*f result: use IDENTITY for plain, fabs for abs. */
#define TIDAL_P22(M_mat, Xr_, Xi_, P_out, post_fn)                            \
    do {                                                                       \
        double *_f22_ = malloc(Nmax * sizeof(double));                         \
        batch_quadratic_form(Np1, Nmax, M_mat, Xr_, Xi_,                      \
                             _f22_, Mr_ws, Mi_ws, Y1_ws, Y2_ws);              \
        for (int j_ = 0; j_ < Nmax; j_++)                                     \
            P_out[j_] = post_fn(kdisc[j_] * kdisc[j_] * kdisc[j_] * _f22_[j_]); \
        free(_f22_);                                                           \
    } while (0)
#define IDENTITY(x) (x)

/* Compute fNL P12 pair (base + ortho) from a single X matrix.
 * Allocates workspace, runs two quadratic forms, applies Tbin*k^3, frees. */
#define COMPUTE_FNL_P12_PAIR(M_base, M_ortho, Xr_, Xi_, P_out, P_out_ortho)   \
    do {                                                                       \
        double *_f12_ = malloc(Nmax * sizeof(double));                         \
        batch_quadratic_form(Np1, Nmax, M_base, Xr_, Xi_,                     \
                             _f12_, Mr_ws, Mi_ws, Y1_ws, Y2_ws);              \
        for (int j_ = 0; j_ < Nmax; j_++)                                     \
            P_out[j_] = Tbin[j_] * kdisc[j_] * kdisc[j_] * kdisc[j_] * _f12_[j_]; \
        batch_quadratic_form(Np1, Nmax, M_ortho, Xr_, Xi_,                    \
                             _f12_, Mr_ws, Mi_ws, Y1_ws, Y2_ws);              \
        for (int j_ = 0; j_ < Nmax; j_++)                                     \
            P_out_ortho[j_] = Tbin[j_] * kdisc[j_] * kdisc[j_] * kdisc[j_] * _f12_[j_]; \
        free(_f12_);                                                           \
    } while (0)

/* Fill M12 bias matrix from source with coefficient, then compute P12 pair.
 * Only runs when has_fNL (fNL enabled). */
#define FILL_M12_BIAS(M_dst, M_dst_o, coeff, M_src, M_src_o, P_out, P_out_o, cs, kb_t_c, sc_c_tab, sc_s_tab) \
    count2 = 0;                                                                \
    if (has_fNL) {                                                   \
        for (index_l = 0; index_l < Nmax + 1; index_l++)                       \
            for (index_i = index_l; index_i < Nmax + 1; index_i++) {           \
                (M_dst)[count2] = (coeff) * (M_src)[count2];                   \
                (M_dst_o)[count2] = (coeff) * (M_src_o)[count2];              \
                count2++;                                                      \
            }                                                                  \
        COMPUTE_P12_BIAS_PAIR(M_dst, M_dst_o, P_out, P_out_o, cs, kb_t_c, sc_c_tab, sc_s_tab);\
    }

/* Fill M22 bias matrix with kernel expression (in nu1, nu2), then assemble P22.
 * nu1 = -0.5*etam2[index_i], nu2 = -0.5*etam2[index_l] are set automatically. */
#define FILL_M22_BIAS(M_out, kernel_expr, P_out)                               \
    {                                                                          \
        int _Nb_ = Nmax + 1;                                                   \
        _Pragma("omp parallel for schedule(static) if(omp_get_max_threads() <= 8)") \
        for (int _lb_ = 0; _lb_ < _Nb_; _lb_++) {                             \
            int _cbb_ = _lb_ * _Nb_ - _lb_ * (_lb_ - 1) / 2;                 \
            for (int _ib_ = _lb_; _ib_ < _Nb_; _ib_++) {                      \
                int _idx_ = _cbb_ + (_ib_ - _lb_);                            \
                double complex nu1 = -0.5 * etam2[_ib_];                       \
                double complex nu2 = -0.5 * etam2[_lb_];                       \
                (M_out)[_idx_] = pnlpt->M22basic_oneline_complex[_idx_] * (kernel_expr); \
            }                                                                  \
        }                                                                      \
    }                                                                          \
    ASSEMBLE_P22_BIAS(M_out, P_out)

/* FFTLog decomposition: the full pipeline from spectrum P(k) to symmetric
 * FFTLog coefficients c_m, implementing the algorithm of Hamilton (2000):
 *
 *   1. Bias subtraction: input_j = P(k_j) * exp(-j * b * Delta)
 *   2. FFT of biased input to get Fourier coefficients
 *   3. Reorder and apply k_min^{-eta} normalization
 *   4. Halve endpoints (trapezoidal rule correction)
 *
 * Output: cmsym_out[0..Nmax] with symmetrized coefficients.
 * Allocates cmsym_out internally; temporary FFT arrays are freed.
 *
 * @param P_src      power spectrum on kdisc grid (size Nmax)
 * @param b_exp      FFTLog bias exponent (e.g. -0.3 for matter, -1.6 for bias)
 * @param eta_arr    complex frequency array eta_m (size Nmax+1)
 * @param cmsym_out  output: allocated array of Nmax+1 complex coefficients
 *
 * Free variables: Nmax, Nmaxd, kmin, Delta, stepsize, pnlpt->error_message */
#define FFT_TO_CMSYM(P_src, b_exp, eta_arr, cmsym_out)                             \
    do {                                                                             \
        double *_ir, *_ii, *_or, *_oi;                                               \
        class_alloc(_ir, Nmax * sizeof(double), pnlpt->error_message);               \
        class_alloc(_ii, Nmax * sizeof(double), pnlpt->error_message);               \
        class_alloc(_or, Nmax * sizeof(double), pnlpt->error_message);               \
        class_alloc(_oi, Nmax * sizeof(double), pnlpt->error_message);               \
        for (int _k = 0; _k < Nmax; _k++) {                                         \
            _ir[_k] = (P_src)[_k] * exp(-1. * _k * (b_exp) * Delta);                \
            _ii[_k] = 0.;                                                            \
        }                                                                            \
        FFT(_ir, _ii, _or, _oi, Nmax, stepsize);                                    \
        class_alloc(cmsym_out, (Nmax + 1) * sizeof(complex double),                 \
                    pnlpt->error_message);                                           \
        /* kmin^(-eta) = kmin^(-b) * exp(-i*imag(eta)*ln(kmin))                      \
         * Use sincos to avoid expensive cpow */                                     \
        double _km_neg_b = pow(kmin, -(b_exp));                                      \
        double _lnkmin = log(kmin);                                                  \
        for (int _k = 0; _k < Nmax + 1; _k++) {                                     \
            int _idx = (_k < Nmax / 2) ? (Nmax / 2 - _k) : (_k - Nmax / 2);        \
            double _sign = (_k < Nmax / 2) ? -1. : 1.;                              \
            double _omega = -cimag((eta_arr)[_k]) * _lnkmin;                         \
            double _s, _c;                                                           \
            sincos(_omega, &_s, &_c);                                                \
            double complex _kpow = _km_neg_b * (_c + _Complex_I * _s);              \
            (cmsym_out)[_k] = _kpow *                                                \
                (_or[_idx] + _sign * _Complex_I * _oi[_idx]) / Nmaxd;               \
        }                                                                            \
        (cmsym_out)[0] /= 2.;                                                       \
        (cmsym_out)[Nmax] /= 2.;                                                    \
        free(_ir); free(_ii); free(_or); free(_oi);                                  \
    } while (0)

/* ====================================================================== */
/* IR RESUMMATION (Senatore & Zaldarriaga 2015, arXiv:1404.5954)           */
/*                                                                          */
/* Extracts the smooth (no-wiggle) power spectrum P_nw(k) from P_L(k)     */
/* using a DST-II transform of log(k*P), then computes the BAO damping    */
/* scales Sigma^2 and deltaSigma^2, and assembles the IR-resummed spectra: */
/*                                                                          */
/*   P_tree(k) = P_nw(k) + (1 + Sigma^2*k^2) * P_w(k) * exp(-Sigma^2*k^2)*/
/*   P_bin(k)  = P_nw(k) + P_w(k) * exp(-Sigma^2*k^2)                     */
/*                                                                          */
/* where P_w(k) = P_L(k) - P_nw(k) is the oscillatory (wiggle) component.  */
/*                                                                          */
/* The algorithm:                                                           */
/*   1. DST-II of log(k*P) on a fine grid (Nir=65536)                     */
/*   2. Zero Fourier harmonics between Nleft=120 and Nright=240            */
/*      (these contain the BAO bump in Fourier space)                      */
/*   3. Inverse DST to recover P_nw(k)                                    */
/*   4. Integrate (P_L - P_nw) to get Sigma^2 = int P_nw j0(kr)/(6pi^2)  */
/*   5. Interpolate P_nw onto the FFTLog grid and build resummed spectra   */
/* ====================================================================== */
static int nonlinear_pt_ir_resummation(
    struct nonlinear_pt *pnlpt,
    struct background *pba,
    struct primordial *ppm,
    struct thermodynamics *pth,
    int Nmax,
    double *lnk_l,
    double *lnpk_l,
    double *myddlnpk,
    double *kdisc,
    double *Pdisc,
    double *PPRIMdisc,
    /* outputs: pre-allocated arrays of size Nmax */
    double *Pnw,
    double *Pw,
    double *Pbin,
    double *Ptree,
    double *Tnw,
    double *Tw,
    double *Tbin,
    /* output scalars */
    double *SigmaBAO_out,
    double *deltaSigmaBAO_out) {

    int last_index = 0, last_index2 = 0;
    double lnpk_out;

    if (pnlpt->nonlinear_pt_verbose > 0)
        printf("Performing IR resummation...\n");

    /* Step 1: Extract smooth (no-wiggle) P(k) via DST-II of log(k*P) */

    /* DST grid for BAO extraction (Nir must match PT matrix resolution) */
    int Nir = 65536;
    int Nirover2 = Nir / 2;      /* half-grid for odd/even decomposition */
    int Nirby4 = 4 * Nir;        /* zero-padded FFT size */
    double Nird = (double)Nir;
    double kmin2 = 0.00007;       /* k range for DST interpolation [1/Mpc] */
    double kmax2 = 7;
    /* Allocate FFT buffers once, reuse for both forward and inverse transforms */
    double *logkPdiscr, *fft_ir, *fft_ii, *fft_or, *fft_oi;
    double logPbin2, kbin2;
    class_alloc(logkPdiscr, Nir * sizeof(double), pnlpt->error_message);
    class_alloc(fft_ir, Nirby4 * sizeof(double), pnlpt->error_message);
    class_alloc(fft_ii, Nirby4 * sizeof(double), pnlpt->error_message);
    class_alloc(fft_or, Nirby4 * sizeof(double), pnlpt->error_message);
    class_alloc(fft_oi, Nirby4 * sizeof(double), pnlpt->error_message);

    last_index = 0;
    int index_ir = 0;
    /* Zero the imaginary input and real even-indexed slots upfront */
    memset(fft_ii, 0, Nirby4 * sizeof(double));
    memset(fft_ir, 0, Nirby4 * sizeof(double));
    double kstep = (kmax2 - kmin2) / (Nird - 1.);
    double lnk0 = lnk_l[0];
    double exp_lnk0 = exp(lnk0);
    #pragma omp parallel
    {
        int _last_idx = 0;
        #pragma omp for schedule(static)
        for (index_ir = 0; index_ir < Nir; index_ir++) {
            double _kbin2 = kmin2 + index_ir * kstep;

            if (_kbin2 >= exp_lnk0) {
                double _logPbin2;
                array_interpolate_spline(lnk_l, pnlpt->k_size, lnpk_l, myddlnpk, 1, log(_kbin2),
                                         &_last_idx, &_logPbin2, 1, pnlpt->error_message);
                logkPdiscr[index_ir] = log(_kbin2) + _logPbin2;
            } else {
                logkPdiscr[index_ir] = log(_kbin2) + lnpk_l[0] + (ppm->n_s) * (log(_kbin2) - lnk0);
            }

            /* DST-II symmetry: alternating sign for odd-indexed real FFT input */
            double val = ((index_ir % 2 == 0) ? 1. : -1.) * logkPdiscr[index_ir];
            fft_ir[2 * index_ir + 1] = val;
            fft_ir[Nirby4 - 2 * index_ir - 1] = val;
        }
    }

    int stepsize = 1;
    FFT(fft_ir, fft_ii, fft_or, fft_oi, Nirby4, stepsize);

    double *out_ir;
    class_alloc(out_ir, Nir * sizeof(double), pnlpt->error_message);
    int index_ir2 = 0;
    for (index_ir2 = 0; index_ir2 < Nir; index_ir2++)
        out_ir[index_ir2] = fft_or[Nir - index_ir2 - 1];

    /* Step 2: Erase BAO bump from Fourier harmonics and spline-interpolate */

    double *cmodd, *cmeven;
    int *ivar;
    class_alloc(cmodd, Nirover2 * sizeof(double), pnlpt->error_message);
    class_alloc(cmeven, Nirover2 * sizeof(double), pnlpt->error_message);
    class_alloc(ivar, Nirover2 * sizeof(int), pnlpt->error_message);

    for (index_ir = 0; index_ir < Nirover2; index_ir++) {
        ivar[index_ir] = index_ir + 1;
        cmodd[index_ir] = out_ir[2 * index_ir];
        cmeven[index_ir] = out_ir[2 * index_ir + 1];
    }

    /* Fourier harmonic window containing BAO feature: modes [Nleft, Nright)
     * are removed to produce the no-wiggle spectrum P_nw */
    int Nleft = 120;
    int Nright = 240;
    int Nthrow = Nright - Nleft;
    int Nnew = Nirover2 - Nthrow;
    double *cmoddnw, *cmevennw, *inew;
    class_alloc(cmoddnw, Nnew * sizeof(double), pnlpt->error_message);
    class_alloc(cmevennw, Nnew * sizeof(double), pnlpt->error_message);
    class_alloc(inew, Nnew * sizeof(double), pnlpt->error_message);

    for (index_ir = 0; index_ir < Nleft; index_ir++) {
        cmoddnw[index_ir] = cmodd[index_ir];
        cmevennw[index_ir] = cmeven[index_ir];
        inew[index_ir] = index_ir + 1.;
    }
    for (index_ir = Nleft; index_ir < Nnew; index_ir++) {
        cmoddnw[index_ir] = cmodd[index_ir + Nthrow];
        cmevennw[index_ir] = cmeven[index_ir + Nthrow];
        inew[index_ir] = index_ir + 1. + Nthrow;
    }

    last_index = 0;
    double *dd_cmoddnw, *dd_cmevennw;
    double cmodd_newval, cmeven_newval;
    class_alloc(dd_cmoddnw, sizeof(double) * Nnew, pnlpt->error_message);
    SPLINE_SETUP(inew, Nnew, cmoddnw, dd_cmoddnw, _SPLINE_NATURAL_);
    class_alloc(dd_cmevennw, sizeof(double) * Nnew, pnlpt->error_message);
    SPLINE_SETUP(inew, Nnew, cmevennw, dd_cmevennw, _SPLINE_NATURAL_);

    double *cmnew;
    class_alloc(cmnew, Nir * sizeof(double), pnlpt->error_message);
    last_index = 0;
    last_index2 = 0;
    for (index_ir = 0; index_ir < Nirover2; index_ir++) {
        SPLINE_EVAL(inew, Nnew, cmoddnw, dd_cmoddnw, (double)(index_ir + 1), &last_index, &cmodd_newval);
        SPLINE_EVAL(inew, Nnew, cmevennw, dd_cmevennw, (double)(index_ir + 1), &last_index2, &cmeven_newval);
        cmnew[index_ir * 2] = cmodd_newval;
        cmnew[index_ir * 2 + 1] = cmeven_newval;
    }

    /* Step 3: Inverse DST to get smooth P_nw(k) and spline-interpolate */

    /* Reuse FFT buffers for the inverse transform (zero imaginary input) */
    double *out_2;
    class_alloc(out_2, Nir * sizeof(double), pnlpt->error_message);
    fft_ir[0] = cmnew[Nir - 1] * 0.5;
    fft_ir[Nir] = 0.;
    fft_ir[2 * Nir] = -1. * cmnew[Nir - 1] * 0.5;
    fft_ir[3 * Nir] = 0.;

    fft_ii[0] = 0.;
    fft_ii[Nir] = 0.;
    fft_ii[2 * Nir] = 0.;
    fft_ii[3 * Nir] = 0.;

    for (index_ir = 1; index_ir < Nir; index_ir++) {
        double half_cm = 0.5 * cmnew[Nir - 1 - index_ir];
        fft_ir[index_ir] = half_cm;
        fft_ir[4 * Nir - index_ir] = half_cm;
        fft_ir[2 * Nir - index_ir] = -half_cm;
        fft_ir[2 * Nir + index_ir] = -half_cm;

        fft_ii[index_ir] = 0.;
        fft_ii[4 * Nir - index_ir] = 0.;
        fft_ii[2 * Nir - index_ir] = 0.;
        fft_ii[2 * Nir + index_ir] = 0.;
    }

    FFT(fft_ir, fft_ii, fft_or, fft_oi, Nirby4, stepsize);

    for (index_ir = 0; index_ir < Nir; index_ir++)
        out_2[index_ir] = ((index_ir % 2 == 0) ? 1. : -1.) * fft_or[2 * index_ir + 1];

    /* Free FFT buffers (no longer needed) */
    free(fft_ir);
    free(fft_ii);
    free(fft_or);
    free(fft_oi);

    double *Pnw_ir, *knw_ir;
    class_alloc(Pnw_ir, Nir * sizeof(double), pnlpt->error_message);
    class_alloc(knw_ir, Nir * sizeof(double), pnlpt->error_message);
    for (index_ir = 0; index_ir < Nir; index_ir++) {
        knw_ir[index_ir] = kmin2 + index_ir * (kmax2 - kmin2) / (Nird - 1.);
        Pnw_ir[index_ir] = exp(out_2[index_ir] / (2 * Nird)) / (knw_ir[index_ir]);
    }

    free(dd_cmoddnw);
    free(dd_cmevennw);

    double *ddPnw;
    double Pnwval;
    class_alloc(ddPnw, sizeof(double) * Nir, pnlpt->error_message);
    SPLINE_SETUP(knw_ir, Nir, Pnw_ir, ddPnw, _SPLINE_NATURAL_);

    /* Step 4: Compute BAO damping factor Sigma^2 from P_nw */

    double rbao = pth->rs_d;
    int Nint2 = 500;

    double *qint2, *IntegrandBAO, *IntegrandBAO2;
    class_alloc(qint2, (Nint2 + 1) * sizeof(double), pnlpt->error_message);
    class_alloc(IntegrandBAO, (Nint2 + 1) * sizeof(double), pnlpt->error_message);
    class_alloc(IntegrandBAO2, (Nint2 + 1) * sizeof(double), pnlpt->error_message);
    double ks = 0.2 * pba->h;
    #pragma omp parallel
    {
        int _last_idx = 0;
        #pragma omp for schedule(static)
        for (index_ir = 0; index_ir < Nint2 + 1; index_ir++) {
            double _Pnwval;
            qint2[index_ir] = kmin2 * exp(index_ir * log(ks / kmin2) / (Nint2));
            array_interpolate_spline(knw_ir, Nir, Pnw_ir, ddPnw, 1, qint2[index_ir],
                                     &_last_idx, &_Pnwval, 1, pnlpt->error_message);
            { double qr = qint2[index_ir] * rbao;
              double qr2 = qr * qr, qr3 = qr2 * qr;
              double sinqr = sin(qr), cosqr = cos(qr);
            IntegrandBAO[index_ir] = _Pnwval * (1. - 3. * sinqr / qr + 6. * (sinqr / qr3 - cosqr / qr2));
            IntegrandBAO2[index_ir] = -1. * _Pnwval * (3. * cosqr * qr + (-3. + qr2) * sinqr) / qr3; }
        }
    }

    double SigmaBAO = 0., deltaSigmaBAO = 0.;
    for (index_ir = 0; index_ir < Nint2; index_ir++) {
        SigmaBAO += (log(qint2[index_ir + 1]) - log(qint2[index_ir])) * (qint2[index_ir + 1] * IntegrandBAO[index_ir + 1] + qint2[index_ir] * IntegrandBAO[index_ir]) / 2. / (6. * (M_PI * M_PI));
        deltaSigmaBAO += (log(qint2[index_ir + 1]) - log(qint2[index_ir])) * (qint2[index_ir + 1] * IntegrandBAO2[index_ir + 1] + qint2[index_ir] * IntegrandBAO2[index_ir]) / 2. / (2. * (M_PI * M_PI));
    }

    /* Step 5: Compute LO IR-resummed spectrum: P_nw + P_w * exp(-k^2 * Sigma^2) */

    #pragma omp parallel
    {
    int _last_idx = 0;
    #pragma omp for schedule(static)
    for (int i = 0; i < Nmax; i++) {
        double _Pnwval2, _Pnwval_rescaled, _Pwval_rescaled;
        if (kdisc[i] <= kmax2 && kdisc[i] >= kmin2 && pnlpt->alpha_rs * kdisc[i] <= kmax2 && pnlpt->alpha_rs * kdisc[i] >= kmin2) {
            array_interpolate_spline(knw_ir, Nir, Pnw_ir, ddPnw, 1, kdisc[i],
                                     &_last_idx, &_Pnwval2, 1, pnlpt->error_message);
            array_interpolate_spline(knw_ir, Nir, Pnw_ir, ddPnw, 1, pnlpt->alpha_rs * kdisc[i],
                                     &_last_idx, &_Pnwval_rescaled, 1, pnlpt->error_message);

            if (kdisc[i] >= exp(lnk_l[0])) {
                double _lnpk_out;
                array_interpolate_spline(lnk_l, pnlpt->k_size, lnpk_l, myddlnpk, 1,
                                         log(pnlpt->alpha_rs * kdisc[i]),
                                         &_last_idx, &_lnpk_out, 1, pnlpt->error_message);
                _Pwval_rescaled = exp(_lnpk_out);
            } else {
                _Pwval_rescaled = exp(lnpk_l[0]) * pow(pnlpt->alpha_rs * kdisc[i] / exp(lnk_l[0]), ppm->n_s);
            }

            Pnw[i] = _Pnwval2;
            Tnw[i] = sqrt(_Pnwval2 / PPRIMdisc[i]) * (5. / 3.);
            Pw[i] = _Pwval_rescaled - _Pnwval_rescaled;
            Tw[i] = Tnw[i] * (Pw[i] / _Pnwval2 / 2. - Pw[i] * Pw[i] / _Pnwval2 / _Pnwval2 / 8. + Pw[i] * Pw[i] * Pw[i] / _Pnwval2 / _Pnwval2 / _Pnwval2 / 16.);

            if (pnlpt->no_wiggle)
                Pw[i] = 0.;
            if (pnlpt->wiggle_only)
                Pnw[i] = 0.;

            Pbin[i] = Pnw[i] + Pw[i] * exp(-SigmaBAO * kdisc[i] * kdisc[i]);
            Tbin[i] = Tnw[i] + Tw[i] * exp(-SigmaBAO * kdisc[i] * kdisc[i]);
        } else {
            Pnw[i] = Pdisc[i];
            Pw[i] = 0.;
            Pbin[i] = Pdisc[i];
            Tnw[i] = sqrt(Pdisc[i] / PPRIMdisc[i]) * (5. / 3.);
            Tw[i] = 0.;
            Tbin[i] = Tnw[i];
        }
        Ptree[i] = Pnw[i] + Pw[i] * exp(-SigmaBAO * kdisc[i] * kdisc[i]) * (1. + SigmaBAO * kdisc[i] * kdisc[i]);
    }
    } /* end omp parallel */

    *SigmaBAO_out = SigmaBAO;
    *deltaSigmaBAO_out = deltaSigmaBAO;

    { void *_f[] = {
        ddPnw, out_ir, out_2, cmnew, logkPdiscr,
        cmeven, cmodd, ivar, inew, cmoddnw,
        cmevennw, knw_ir, Pnw_ir, qint2, IntegrandBAO,
        IntegrandBAO2,
    };
    for (int fi_ = 0; fi_ < 16; fi_++) free(_f[fi_]); }

    return _SUCCESS_;
}

/**
 * Compute all one-loop power spectra at a single redshift.
 *
 * ALGORITHM OUTLINE:
 *   1. Set up FFTLog grid: k_j = kmin * exp(j*Delta), j=0..Nmax-1
 *   2. Interpolate linear P_L(k) and primordial P_prim(k) onto the grid
 *   3. IR resummation: split P_L = P_nw + P_w, compute BAO damping Sigma^2
 *   4. FFTLog decomposition: P_L(k) -> c_m coefficients via FFT
 *   5. One-loop integrals via BLAS:
 *      - P13[k] = sum_m c_m * M13_m * k^{eta_m}  (dot product)
 *      - P22[k] = c^T * M22 * c                   (quadratic form)
 *      - P12[k] = c_transfer^T * M12 * c           (bilinear, fNL only)
 *   6. Assemble EFT counterterm P_CTR = k^2 * P_L
 *   7. Interpolate 1-loop results to output k-grid
 *   8. RSD multipoles (ell=0,2,4): compute P_{ell}^{vv,vd,dd}
 *      - Fill M22/M13 kernels with RSD factors (f-dependent polynomials)
 *      - IR resummation: separate NW and wiggle P22/P13/P12
 *      - Alcock-Paczynski: Gauss-Legendre integration over mu
 *   9. Biased tracers: I[d2d2], I[d2G2], I[G2G2], I[FG2], fNL bias
 *
 * @param tau     conformal time at redshift z
 * @param f       logarithmic growth rate d(ln D)/d(ln a)
 * @param hratio  H_fid/H_true for Alcock-Paczynski rescaling
 * @param Dratio  D_A,true/D_A,fid for Alcock-Paczynski rescaling
 * @param pk_l    linear power spectrum P_L(k) on output grid (input)
 * @param pk_l_X_Y  output arrays for 1-loop multipoles (see below)
 */
int nonlinear_pt_loop(
    struct precision *ppr,
    struct background *pba,
    struct primordial *ppm,
    struct thermodynamics *pth,
    struct nonlinear_pt *pnlpt,
    double tau,
    double f,
    double hratio,
    double Dratio,
    double *pk_l,
    double *pPRIMk_l,
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
    double *pk_nl_fNL,
    double *pk_fNLd2,
    double *pk_fNLG2,
    double *pk_nl_fNL_ortho,
    double *pk_fNLd2_ortho,
    double *pk_fNLG2_ortho,
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
    double *lnk_l,
    double *lnpk_l,
    double *lnpPRIMk_l) {
    int index_k = 0;
    int index_j = 0;
    int index_i = 0;
    int index_l = 0;

    double sigmav = 0.;
    double *pvecback;

    int last_index = 0;
    int last_index2 = 0;

    TIMER_DECL(blas); TIMER_DECL(fft); TIMER_DECL(buildx);
    TIMER_DECL(ap_spline); TIMER_DECL(rsd_ap); TIMER_DECL(bias_ap);
    TIMER_DECL(spline_out); TIMER_DECL(m22fill);
    TIMER_DECL(alloc); TIMER_DECL(ir); TIMER_DECL(macro);
    struct timespec _ts_total0, _ts_total1;
    clock_gettime(CLOCK_MONOTONIC, &_ts_total0);

    class_alloc(pvecback, pba->bg_size * sizeof(double), pnlpt->error_message);

    class_call(background_at_tau(pba, tau, long_info, inter_normal, &last_index, pvecback),
               pba->error_message,
               pnlpt->error_message);

    if (pnlpt->nonlinear_pt_verbose > 0)
        printf("Computing one-loop power spectra at z=%e\n", 1. / pvecback[pba->index_bg_a] - 1.); /* a_today=1 in v3 */

    free(pvecback);

    /* --- FFTLog discretization parameters --- */
    int Nmax = ppr->nmax_nlpt; /* Number of FFTLog sampling points */
    int Nmaxf = ppr->nmax_nlpt + 1;
    double Nmaxd = Nmax * 1.;
    double kmin = 0.00005 * pba->h; /* FFTLog k-range (must match PT matrices) */
    double kmax = 100. * pba->h;
    double *js, *kdisc, *Pbin, *Pdisc, *Ptree, *PPRIMdisc, *Tbin;
    double *P10, *P10b1, *P12;
    double *Pnw, *Pw, *dd_Pnw, *dd_Pw, *Tnw, *Tw, *dd_Tnw, *dd_Tw;
    class_alloc(js, (Nmax + 1) * sizeof(double), pnlpt->error_message);
    { /* Allocate FFTLog grid arrays */
        double **_ptrs[] = {&kdisc, &Pbin, &Pdisc, &Ptree, &PPRIMdisc, &Tbin,
                            &P10, &P10b1, &P12, &Pnw, &Pw, &dd_Pnw, &dd_Pw,
                            &Tnw, &Tw, &dd_Tnw, &dd_Tw};
        for (int ai_ = 0; ai_ < 17; ai_++)
            class_alloc(*_ptrs[ai_], Nmax * sizeof(double), pnlpt->error_message);
    }
    double SigmaBAO = 0.;
    double deltaSigmaBAO = 0.;
    double Delta = log(kmax / kmin) / (Nmaxd - 1);

    /* Interpolate linear P(k) onto FFTLog grid */

    double lnpk_out = 0.;
    double *myddlnpk;
    double lnpPRIMk_out = 0.;
    double *myddlnpPRIMk;
    class_alloc(myddlnpk, sizeof(double) * pnlpt->k_size, pnlpt->error_message);
    SPLINE_SETUP(lnk_l, pnlpt->k_size, lnpk_l, myddlnpk, _SPLINE_NATURAL_);
    class_alloc(myddlnpPRIMk, sizeof(double) * pnlpt->k_size, pnlpt->error_message);
    SPLINE_SETUP(lnk_l, pnlpt->k_size, lnpPRIMk_l, myddlnpPRIMk, _SPLINE_NATURAL_);

    last_index = 0;
    int i_kdisc = 0;

    for (i_kdisc = 0; i_kdisc < Nmax; i_kdisc++) {
        kdisc[i_kdisc] = kmin * exp(i_kdisc * Delta);

        if (kdisc[i_kdisc] >= exp(lnk_l[0])) {
            SPLINE_EVAL(lnk_l, pnlpt->k_size, lnpk_l, myddlnpk, log(kdisc[i_kdisc]), &last_index, &lnpk_out);
            SPLINE_EVAL(lnk_l, pnlpt->k_size, lnpPRIMk_l, myddlnpPRIMk, log(kdisc[i_kdisc]), &last_index, &lnpPRIMk_out);

            Pdisc[i_kdisc] = exp(lnpk_out);
            PPRIMdisc[i_kdisc] = exp(lnpPRIMk_out);

        } /* end k in CLASS range */
        else
        {
            Pdisc[i_kdisc] = exp(lnpk_l[0]) * pow(kdisc[i_kdisc] / exp(lnk_l[0]), ppm->n_s);
            PPRIMdisc[i_kdisc] = exp(lnpPRIMk_l[0]) * pow(kdisc[i_kdisc] / exp(lnk_l[0]), ppm->n_s);
        }
    }

    int Nside;

    /*
     To be safe we're shifting the range of kmaxes and kmins used in the analysis in order to be able to compute the AP,
     this should be increased if needed
     */

    if (pnlpt->AP_effect == AP_effect_yes) {
        Nside = 10;
    } else {
        Nside = 1;
    }
    double kmaxnew = kdisc[Nmax - 1 - Nside];
    double kminnew = kdisc[Nside];

    /* --- Compute velocity dispersion sigma_v^2 = int P_L(k)/(6pi^2) dk --- */
    /* Used for Fingers-of-God damping in IR resummation */
    sigmav = 0.;
    for (index_k = 0; index_k < pnlpt->k_size - 1; index_k++) {
        sigmav += (lnk_l[index_k + 1] - lnk_l[index_k]) * (exp(lnk_l[index_k + 1]) * exp(lnpk_l[index_k + 1]) + exp(lnpk_l[index_k]) * exp(lnk_l[index_k])) / 2. / (6. * (M_PI * M_PI));
    }

    /* --- IR Resummation (Senatore & Zaldarriaga method) --- */
    /* irindex=0: no-wiggle contributions only
     * irindex=1: wiggly (BAO) contributions included via damped oscillations */
    int irindex = 0;

    if (pnlpt->irres == irres_yes) {
        irindex = 1;
        TIMER_START(ir);
        class_call(nonlinear_pt_ir_resummation(
                       pnlpt, pba, ppm, pth, Nmax,
                       lnk_l, lnpk_l, myddlnpk, kdisc, Pdisc, PPRIMdisc,
                       Pnw, Pw, Pbin, Ptree, Tnw, Tw, Tbin,
                       &SigmaBAO, &deltaSigmaBAO),
                   pnlpt->error_message, pnlpt->error_message);
        TIMER_ADD(ir);
    } /* End of IR resummation conditional expression */

    else
    {
        if (pnlpt->nonlinear_pt_verbose > 0)
            printf("IR resummation skipped.\n");
        int index_kd = 0;
        for (index_kd = 0; index_kd < Nmax; index_kd++) {
            Pbin[index_kd] = Pdisc[index_kd];
            Ptree[index_kd] = Pbin[index_kd];
            Pnw[index_kd] = Pdisc[index_kd];
            Pw[index_kd] = 0.;
            Tbin[index_kd] = pow(Pdisc[index_kd] / PPRIMdisc[index_kd], 0.5) * (5. / 3.);
            Tnw[index_kd] = Tbin[index_kd];
            Tw[index_kd] = 0.;
        }
    }

    class_call(array_spline_table_columns(kdisc, Nmax, Pnw, 1, dd_Pnw, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
    last_index = 0;
    double Pnw_ap_out = 0;

    class_call(array_spline_table_columns(kdisc, Nmax, Pw, 1, dd_Pw, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
    double Pw_ap_out = 0;

    class_call(array_spline_table_columns(kdisc, Nmax, Tnw, 1, dd_Tnw, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
    last_index = 0;
    double Tnw_ap_out = 0;

    class_call(array_spline_table_columns(kdisc, Nmax, Tw, 1, dd_Tw, _SPLINE_NATURAL_, pnlpt->error_message), pnlpt->error_message, pnlpt->error_message);
    double Tw_ap_out = 0;

    /* --- FFTLog basis: decompose P(k) = sum_j c_j * k^{b + 2*pi*i*j/(N*Delta)} --- */
    /* etam[j] = b + 2*pi*i*j/(N*Delta) are the complex FFTLog exponents;
     * b is the bias parameter controlling convergence */

    double complex *etam;
    double complex *etam_transfer;
    class_alloc(etam, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);
    class_alloc(etam_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

    int index_c = 0;
    double b = -0.3;
    double b_transfer = -0.8;
    for (index_c = 0; index_c < Nmax + 1; index_c++) {
        js[index_c] = index_c - Nmaxd / 2;
        etam[index_c] = b + 2. * M_PI * _Complex_I * js[index_c] / Nmaxd / Delta;
        etam_transfer[index_c] = b_transfer + 2. * M_PI * _Complex_I * js[index_c] / Nmaxd / Delta;
    }

    int stepsize = 1;

    /* --- FFT of P_L(k) and T(k) to get FFTLog coefficients c_m --- */
    TIMER_START(fft);
    double complex *cmsym;
    FFT_TO_CMSYM(Pbin, b, etam, cmsym);
    double complex *cmsym_transfer;
    FFT_TO_CMSYM(Tbin, b_transfer, etam_transfer, cmsym_transfer);
    TIMER_ADD(fft);

    /* Load precomputed PT kernel matrices into complex arrays */

    double complex nu1, nu2, nu12;

    /* Offset constants for log-safe storage: output arrays store
     * P(k) + offset so that values remain positive for log interpolation
     * in the spectra module. Different offsets suit different magnitude ranges:
     *   large_for_logs_big    = 1e6  (1-loop RSD multipoles, bias P22)
     *   large_for_logs_matter = 5e4  (matter 1-loop)
     *   large_for_logs_fNL    = 5e4  (fNL P12 corrections)
     *   large_for_logs_small  = 100  (tidal field integrals Id2, IG2, IFG2)
     *   epsilon_for_logs_fNL  = 1e-6 (floor for extrapolation regions) */
    double epsilon_for_logs_fNL = 1.e-6;
    double large_for_logs_fNL = 5.e4;
    double large_for_logs_matter = 5.e4;
    double large_for_logs_big = 1000000.;
    double large_for_logs_small = 100.;

    /* ================================================================== */
    /* === MATTER ONE-LOOP INTEGRALS (P13 + P22) === */
    /* ================================================================== */
    /* The one-loop SPT integrals are computed via matrix multiplication
     * on the FFTLog basis coefficients c_m:
     *   f13[k] = sum_m c_m * M13_m * k^{eta_m}  (vector dot product)
     *   f22[k] = c^T * M22 * c                   (symmetric quadratic form)
     * where M13, M22 are precomputed kernel matrices encoding the F_2 and
     * F_3 symmetrized perturbation theory kernels (EdS approximation).
     *
     * P13UV is the UV counterterm that ensures UV convergence.
     * P_CTR = k^2 * P_L is the leading EFT counterterm (speed of sound).
     * cutoff = 3*h applies exp(-(k/cutoff)^6) to suppress numerical noise. */
    double cutoff = 3. * pba->h;
    double *P13;
    double *P13UV;
    double *P1loop;
    double *f13; /* Now real-valued */
    class_alloc(P13, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(P13UV, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(P1loop, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(f13, Nmax * sizeof(double), pnlpt->error_message);

    double *f22; /* Now real-valued */
    double *P22;
    double *P_CTR;
    class_alloc(f22, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(P22, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(P_CTR, Nmax * sizeof(double), pnlpt->error_message);

    double *P12_fNL;
    double *f12; /* Now real-valued */
    class_alloc(P12_fNL, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(f12, Nmax * sizeof(double), pnlpt->error_message);

    double *P12_fNL_ortho;
    class_alloc(P12_fNL_ortho, Nmax * sizeof(double), pnlpt->error_message);

    int count = 0;

    /* has_fNL: controls whether fNL (equilateral/orthogonal) P12 corrections
     * and their orthogonal variants are computed. Set to 1 when fNL is enabled. */
    int has_fNL = 0;
    if (pnlpt->fNL_equil_ortho_switch == fNL_equil_ortho_yes)
        has_fNL = 1;

    int Np1 = Nmax + 1; /* Number of FFTLog coefficients */

    double complex *x;
    double complex *x_w;
    class_alloc(x, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

    /* --- OpenMP + BLAS thread control ---
     * Adaptive threading strategy:
     *   OMP ≤ 8: Use task-level parallelism with single-threaded BLAS
     *   OMP > 8: Use multi-threaded BLAS for coarse-grained work
     * The fine-grained loops (M22 fills, build_X) have if(omp <= 8) guards,
     * so at high thread counts they run serially while BLAS can use threads. */
    int _blas_threads_saved = openblas_get_num_threads();
    int _n_threads = omp_get_max_threads();

    if (_n_threads <= 8) {
        /* Low thread count: single-threaded BLAS, parallelize tasks */
        openblas_set_num_threads(1);
    } else if (_n_threads <= 16) {
        /* Medium: 2-thread BLAS for each of the ~9 parallel BLAS dispatch tasks */
        openblas_set_num_threads(2);
    } else if (_n_threads <= 32) {
        /* High: 4-thread BLAS, fewer parallel tasks */
        openblas_set_num_threads(4);
    } else {
        /* Very high (64+): 8-thread BLAS, minimal parallel tasks */
        openblas_set_num_threads(8);
    }

    /* Workspace for batch real BLAS operations — per-thread arrays */
    double *Xr, *Xi; /* X matrices: Np1 x Nmax */
    double *Xr_transfer, *Xi_transfer;

    /* Per-thread workspace: Mr, Mi (Np1*Np1), Y1, Y2 (Np1*Nmax), dot_mr, dot_mi (Np1) */
    double **ws_Mr, **ws_Mi, **ws_Y1, **ws_Y2, **ws_dot_mr, **ws_dot_mi;
    ws_Mr = (double **)malloc(_n_threads * sizeof(double *));
    ws_Mi = (double **)malloc(_n_threads * sizeof(double *));
    ws_Y1 = (double **)malloc(_n_threads * sizeof(double *));
    ws_Y2 = (double **)malloc(_n_threads * sizeof(double *));
    ws_dot_mr = (double **)malloc(_n_threads * sizeof(double *));
    ws_dot_mi = (double **)malloc(_n_threads * sizeof(double *));
    for (int _t = 0; _t < _n_threads; _t++) {
        ws_Mr[_t]     = (double *)calloc(Np1 * Np1, sizeof(double));
        ws_Mi[_t]     = (double *)calloc(Np1 * Np1, sizeof(double));
        ws_Y1[_t]     = (double *)calloc(Np1 * Nmax, sizeof(double));
        ws_Y2[_t]     = (double *)calloc(Np1 * Nmax, sizeof(double));
        ws_dot_mr[_t] = (double *)calloc(Np1, sizeof(double));
        ws_dot_mi[_t] = (double *)calloc(Np1, sizeof(double));
    }
    /* Thread-0 aliases for serial code paths (backward compatibility) */
    double *Mr_ws = ws_Mr[0], *Mi_ws = ws_Mi[0];
    double *Y1_ws = ws_Y1[0], *Y2_ws = ws_Y2[0];
    double *dot_mr_ws = ws_dot_mr[0], *dot_mi_ws = ws_dot_mi[0];

    class_alloc(Xr, Np1 * Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(Xi, Np1 * Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(Xr_transfer, Np1 * Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(Xi_transfer, Np1 * Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(x_w, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);
    double complex *x_transfer;
    class_alloc(x_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

    double complex *x_w_transfer;

    class_alloc(x_w_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

    /* Precompute k^b and ln(k) arrays for reuse across build_X_matrix calls.
     * Four distinct b values: b=-0.3, b_transfer=-0.8, b2=-1.6000001, b2_transfer=-1.25 */
    double *lnk_cache, *kb_cache, *kb_transfer_cache;
    class_alloc(lnk_cache, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(kb_cache, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(kb_transfer_cache, Nmax * sizeof(double), pnlpt->error_message);
    for (int j = 0; j < Nmax; j++) {
        lnk_cache[j] = log(kdisc[j]);
        kb_cache[j] = pow(kdisc[j], b);
        kb_transfer_cache[j] = pow(kdisc[j], b_transfer);
    }

    /* Precompute sincos tables for each unique etam, to share across
     * all build_X calls with the same exponents (only coeffs differ). */
    double *sc_cos, *sc_sin;           /* sincos table for etam */
    double *sc_cos_t, *sc_sin_t;       /* sincos table for etam_transfer */
    class_alloc(sc_cos, Np1 * Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(sc_sin, Np1 * Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(sc_cos_t, Np1 * Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(sc_sin_t, Np1 * Nmax * sizeof(double), pnlpt->error_message);

    TIMER_START(buildx);
    precompute_sincos_table(Nmax, Np1, etam, lnk_cache, sc_cos, sc_sin);
    precompute_sincos_table(Nmax, Np1, etam_transfer, lnk_cache, sc_cos_t, sc_sin_t);

    /* Build X matrices for the primary FFTLog basis (matter) */
    build_X_from_tables(Nmax, Np1, cmsym, kb_cache, sc_cos, sc_sin, Xr, Xi);
    build_X_from_tables(Nmax, Np1, cmsym_transfer, kb_transfer_cache, sc_cos_t, sc_sin_t, Xr_transfer, Xi_transfer);
    TIMER_ADD(buildx);

    /* Compute basic P13 using batch dot product */
    {
        double *f13_real = malloc(Nmax * sizeof(double));
        TIMER_START(blas);
        batch_dot_product(Np1, Nmax, Xr, Xi, pnlpt->M13_oneline_complex, f13_real, dot_mr_ws, dot_mi_ws);
        TIMER_ADD(blas);
        for (index_j = 0; index_j < Nmax; index_j++) {
            P13UV[index_j] = -61. * Pbin[index_j] * kdisc[index_j] * kdisc[index_j] * sigmav / 105.;
            { double kc_ = kdisc[index_j] / cutoff; double kc2_ = kc_ * kc_; double kc6_ = kc2_ * kc2_ * kc2_;
            P13[index_j] = (kdisc[index_j] * kdisc[index_j] * kdisc[index_j] * f13_real[index_j] * Pbin[index_j] + P13UV[index_j]) * exp(-kc6_); }
        }
        free(f13_real);
    }
    index_j = 0;
    count = 0;
    /* Compute basic P22 using batch quadratic form */
    {
        double *f22_real = malloc(Nmax * sizeof(double));
        TIMER_START(blas);
        batch_quadratic_form(Np1, Nmax, pnlpt->M22_oneline_complex, Xr, Xi, f22_real, Mr_ws, Mi_ws, Y1_ws, Y2_ws);
        TIMER_ADD(blas);
        for (index_j = 0; index_j < Nmax; index_j++) {
            P22[index_j] = kdisc[index_j] * kdisc[index_j] * kdisc[index_j] * f22_real[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.));
            P1loop[index_j] = P13[index_j] + P22[index_j];
        }
        free(f22_real);
    }

    if (has_fNL) {
        double *f12_real = malloc(Nmax * sizeof(double));
        TIMER_START(blas);
        batch_quadratic_form(Np1, Nmax, pnlpt->M12_oneline_complex, Xr_transfer, Xi_transfer, f12_real, Mr_ws, Mi_ws, Y1_ws, Y2_ws);
        TIMER_ADD(blas);
        for (index_j = 0; index_j < Nmax; index_j++) {
            P12_fNL[index_j] = Tbin[index_j] * kdisc[index_j] * kdisc[index_j] * kdisc[index_j] * f12_real[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.));
        }

        TIMER_START(blas);
        batch_quadratic_form(Np1, Nmax, pnlpt->M12_oneline_complex_ortho, Xr_transfer, Xi_transfer, f12_real, Mr_ws, Mi_ws, Y1_ws, Y2_ws);
        TIMER_ADD(blas);
        for (index_j = 0; index_j < Nmax; index_j++) {
            P12_fNL_ortho[index_j] = Tbin[index_j] * kdisc[index_j] * kdisc[index_j] * kdisc[index_j] * f12_real[index_j] * exp(-pow(kdisc[index_j] / cutoff, 6.));
        }
        free(f12_real);
    }

    /* --- Assemble EFT counterterm and spline-interpolate all 1-loop results --- */
    for (index_j = 0; index_j < Nmax; index_j++) {
        P_CTR[index_j] = kdisc[index_j] * kdisc[index_j] * Pbin[index_j];
    }

    double *ddpk_nl;
    double *ddpk_nl_fNL;
    double *ddpk_nl_fNL_ortho;
    double *ddpk_CTR;
    double *ddpk_Tree;
    class_alloc(ddpk_nl, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(ddpk_nl_fNL, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(ddpk_nl_fNL_ortho, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(ddpk_CTR, Nmax * sizeof(double), pnlpt->error_message);
    class_alloc(ddpk_Tree, Nmax * sizeof(double), pnlpt->error_message);

    SPLINE_SETUP(kdisc, Nmax, P1loop, ddpk_nl, _SPLINE_NATURAL_);
    if (has_fNL) {
        SPLINE_SETUP(kdisc, Nmax, P12_fNL, ddpk_nl_fNL, _SPLINE_NATURAL_);
        SPLINE_SETUP(kdisc, Nmax, P12_fNL_ortho, ddpk_nl_fNL_ortho, _SPLINE_NATURAL_);
    }
    SPLINE_SETUP(kdisc, Nmax, P_CTR, ddpk_CTR, _SPLINE_NATURAL_);
    SPLINE_SETUP(kdisc, Nmax, Ptree, ddpk_Tree, _SPLINE_NATURAL_);

    double pk_nl_out;
    double pk_nl_fNL_out;
    double pk_nl_fNL_out_ortho;
    double pk_CTR_out;
    double pk_Tree_out;

    last_index = 0;
    last_index2 = 0;
    int last_index3 = 0;
    int last_index4 = 0;

    int last_index4_ortho = 0;

    /* --- Interpolate matter 1-loop results to the output k-grid --- */
    for (index_k = 0; index_k < pnlpt->k_size; index_k++) {
        if (pnlpt->k[index_k] >= kmin && pnlpt->k[index_k] <= kmax) {
            SPLINE_EVAL(kdisc, Nmax, P1loop, ddpk_nl, pnlpt->k[index_k], &last_index2, &pk_nl_out);
            SPLINE_EVAL(kdisc, Nmax, P_CTR, ddpk_CTR, pnlpt->k[index_k], &last_index, &pk_CTR_out);
            SPLINE_EVAL(kdisc, Nmax, Ptree, ddpk_Tree, pnlpt->k[index_k], &last_index3, &pk_Tree_out);

            pk_nl[index_k] = pk_nl_out + large_for_logs_matter;

            if (has_fNL) {
                SPLINE_EVAL(kdisc, Nmax, P12_fNL, ddpk_nl_fNL, pnlpt->k[index_k], &last_index4, &pk_nl_fNL_out);
                SPLINE_EVAL(kdisc, Nmax, P12_fNL_ortho, ddpk_nl_fNL_ortho, pnlpt->k[index_k], &last_index4_ortho, &pk_nl_fNL_out_ortho);
                pk_nl_fNL[index_k] = pk_nl_fNL_out + large_for_logs_fNL;
                pk_nl_fNL_ortho[index_k] = pk_nl_fNL_out_ortho + large_for_logs_fNL;
            } else {
                pk_nl_fNL[index_k] = large_for_logs_fNL;
                pk_nl_fNL_ortho[index_k] = large_for_logs_fNL;
            }

            if (pk_CTR_out <= 0)
                pk_CTR_out = 1e-16;
            if (pk_Tree_out <= 0)
                pk_Tree_out = 1e-16;
            pk_CTR[index_k] = pk_CTR_out;
            pk_Tree[index_k] = pk_Tree_out;
        } else {
            if (pnlpt->k[index_k] < kmin) {
                pk_nl[index_k] = large_for_logs_matter - 61. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav / 105.;
            }
            if (pnlpt->k[index_k] > kmax) {
                pk_nl[index_k] = large_for_logs_matter;
            }
            pk_nl_fNL[index_k] = large_for_logs_fNL + epsilon_for_logs_fNL;
            pk_nl_fNL_ortho[index_k] = large_for_logs_fNL + epsilon_for_logs_fNL;
            pk_CTR[index_k] = exp(lnpk_l[index_k] + 2. * lnk_l[index_k]);
            pk_Tree[index_k] = exp(lnpk_l[index_k]);
        }
    }

    /* ================================================================== */
    /* === REDSHIFT-SPACE DISTORTIONS (RSD) MULTIPOLES === */
    /* ================================================================== */
    /* Compute P_ell(k) for ell=0,2,4 (monopole, quadrupole, hexadecapole).
     * Each multipole is decomposed into velocity-velocity (vv),
     * velocity-density (vd), and density-density (dd) components,
     * with powers of the growth rate f from the RSD kernel Z_n.
     * Two IR passes: irindex=0 uses P_nw, irindex=1 adds P_w corrections. */

    if (pnlpt->rsd == rsd_yes) {
        if (pnlpt->nonlinear_pt_verbose > 0)
            printf("Computing RSD...\n");
        if (pnlpt->nonlinear_pt_verbose > 0)
            printf("Logarithmic growth factor f=%.16e\n", f);

        /* RSD 1-loop power spectrum arrays: tree, P13, P13UV, P1loop, P22, P12 per (ell, component) */
        double *Ptree_0_vv, *Ptree_0_vd, *Ptree_0_dd, *Ptree_2_vv, *Ptree_2_vd, *Ptree_2_dd;
        double *Ptree_4_vv, *Ptree_4_vd;
        double *P13_0_vv, *P13UV_0_vv, *P1loop_0_vv, *P22_0_vv, *P12_0_vv, *P12_0_vv_ortho;
        double *P13_0_vd, *P13UV_0_vd, *P1loop_0_vd, *P22_0_vd, *P12_0_vd, *P12_0_vd_ortho;
        double *P13_0_dd, *P13UV_0_dd, *P1loop_0_dd, *P22_0_dd, *P12_0_dd, *P12_0_dd_ortho;
        double *P13_2_vv, *P13UV_2_vv, *P1loop_2_vv, *P22_2_vv, *P12_2_vv, *P12_2_vv_ortho;
        double *P13_2_vd, *P13UV_2_vd, *P1loop_2_vd, *P22_2_vd, *P12_2_vd, *P12_2_vd_ortho;
        double *P13_2_dd, *P13UV_2_dd, *P1loop_2_dd, *P22_2_dd, *P12_2_dd, *P12_2_dd_ortho;
        double *P13_4_vv, *P13UV_4_vv, *P1loop_4_vv, *P22_4_vv, *P12_4_vv, *P12_4_vv_ortho;
        double *P13_4_vd, *P13UV_4_vd, *P1loop_4_vd, *P22_4_vd, *P12_4_vd, *P12_4_vd_ortho;
        double *P22_4_dd, *P1loop_4_dd, *P12_4_dd, *P12_4_dd_ortho;
        double *P_CTR_0, *P_CTR_2, *P_CTR_4;
        { /* Batch-allocate all RSD multipole arrays (each Nmax doubles) */
            double **_rsd_ptrs[] = {
                /* P13, P13UV, P1loop */
                &P13_0_vv, &P13UV_0_vv, &P1loop_0_vv, &P13_0_vd, &P13UV_0_vd, &P1loop_0_vd,
                &P13_0_dd, &P13UV_0_dd, &P1loop_0_dd, &P13_2_vv, &P13UV_2_vv, &P1loop_2_vv,
                &P13_2_vd, &P13UV_2_vd, &P1loop_2_vd, &P13_2_dd, &P13UV_2_dd, &P1loop_2_dd,
                &P13_4_vv, &P13UV_4_vv, &P1loop_4_vv, &P13_4_vd, &P13UV_4_vd, &P1loop_4_vd,
                &P22_4_dd, &P1loop_4_dd, &P12_4_dd, &P12_4_dd_ortho,
                &P_CTR_0, &P_CTR_2, &P_CTR_4,
                /* Ptree */
                &Ptree_0_vv, &Ptree_0_vd, &Ptree_0_dd, &Ptree_2_vv, &Ptree_2_vd, &Ptree_2_dd,
                &Ptree_4_vv, &Ptree_4_vd,
                /* P22 */
                &P22_0_vv, &P22_0_vd, &P22_0_dd, &P22_2_vv, &P22_2_vd, &P22_2_dd,
                &P22_4_vv, &P22_4_vd,
                /* P12 */
                &P12_0_vv, &P12_0_vd, &P12_0_dd, &P12_2_vv, &P12_2_vd, &P12_2_dd,
                &P12_4_vv, &P12_4_vd,
                /* P12_ortho */
                &P12_0_vv_ortho, &P12_0_vd_ortho, &P12_0_dd_ortho,
                &P12_2_vv_ortho, &P12_2_vd_ortho, &P12_2_dd_ortho,
                &P12_4_vv_ortho, &P12_4_vd_ortho};
            for (int ai_ = 0; ai_ < 63; ai_++)
                class_alloc(*_rsd_ptrs[ai_], Nmax * sizeof(double), pnlpt->error_message);
        }

        if (irindex == 0) {
            /* --- Fill RSD M22 kernels with f-dependent polynomials ---
             *
             * The M22 kernel matrices for each RSD mu-power are derived from the
             * general M22 kernel by multiplying with rational polynomials in
             * nu1 = -eta_i/2, nu2 = -eta_l/2 (the FFTLog exponents) and powers
             * of the growth rate f. These expressions encode the angular
             * dependence of the redshift-space Z_2 kernel:
             *
             *   Z_2(k1,k2) = b1*F_2(k1,k2) + f*mu^2*k/k1/k2*[...] + ...
             *
             * Each M22 mu-power captures the coefficient of mu^{2n} in the
             * one-loop integrand (see Eq. 2.12-2.15 of arXiv:2004.10607).
             *
             * M12 kernels (fNL) are precomputed and simply combined with
             * powers of f here. */

            TIMER_START(m22fill);
            {
                int _N = Nmax + 1;
                #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
                for (int _l = 0; _l < _N; _l++) {
                    int _count_base = _l * _N - _l * (_l - 1) / 2;
                    for (int _i = _l; _i < _N; _i++) {
                        int _c = _count_base + (_i - _l);
                        double complex _nu1 = -0.5 * etam[_i];
                        double complex _nu2 = -0.5 * etam[_l];
                        double complex _nu12 = _nu1 + _nu2;

                        pnlpt->M22_oneline_0_vv_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * (f * f * (14. * f * f * (24. - 8. * _nu2 - 15. * _nu2 * _nu2 + 5. * _nu2 * _nu2 * _nu2 + 5. * _nu1 * _nu1 * _nu1 * (1. + 7. * _nu2) + 5. * _nu1 * _nu1 * (-3. - 10. * _nu2 + 14. * _nu2 * _nu2) + _nu1 * (-8. - 24. * _nu2 - 50. * _nu2 * _nu2 + 35. * _nu2 * _nu2 * _nu2)) + 18. * f * (36. - 8. * _nu2 + 70. * _nu1 * _nu1 * _nu1 * _nu2 - 23. * _nu2 * _nu2 + _nu1 * _nu1 * (-23. - 94. * _nu2 + 140. * _nu2 * _nu2) + _nu1 * (-8. - 42. * _nu2 - 94. * _nu2 * _nu2 + 70. * _nu2 * _nu2 * _nu2)) + 9. * (50. - 9. * _nu2 + 98. * _nu1 * _nu1 * _nu1 * _nu2 - 35. * _nu2 * _nu2 + 7. * _nu1 * _nu1 * (-5. - 18. * _nu2 + 28. * _nu2 * _nu2) + _nu1 * (-9. - 66. * _nu2 - 126. * _nu2 * _nu2 + 98. * _nu2 * _nu2 * _nu2)))) / 8820.;

                        if (has_fNL) {
                            pnlpt->M12_oneline_0_vv_complex[_c] = f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2[_c]) + f * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f3[_c]);

                            pnlpt->M12_oneline_0_vv_complex_ortho[_c] = f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2_ortho[_c]) + f * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f3_ortho[_c]);
                        }
                    }
                }
            }
            TIMER_ADD(m22fill);

            /* --- Monopole (ell=0), velocity-velocity component: P_{0,vv} ~ f^2 --- */
            /* f13/f22/f12 computed inline by batch BLAS macros */

            for (index_i = 0; index_i < Nmax + 1; index_i++) {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_0_vv_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * (3. * (f * f) * (7. * (-5. + 3. * nu1) + 6. * f * (-7. + 5. * nu1))) / 3920.;
            }

            /* Batch compute P13_0_vv for all k values at once */
            COMPUTE_P13(0_vv, pnlpt->M13_0_vv_oneline_complex,
                        -1. * Pbin[j] * kdisc[j] * kdisc[j] * sigmav * (f * f * (441. + 566. * f + 175. * f * f) / 1225.));

            /* Batch compute P22_0_vv and P12_0_vv for all k values */
            TIMER_START(blas);
            COMPUTE_P22(0_vv, pnlpt->M22_oneline_0_vv_complex);
            COMPUTE_P12(0_vv, pnlpt->M12_oneline_0_vv_complex);
            TIMER_ADD(blas);

            /* Compute orthogonal fNL variant if needed */
            if (has_fNL) {
                TIMER_START(blas);
                COMPUTE_P12(0_vv_ortho, pnlpt->M12_oneline_0_vv_complex_ortho);
                TIMER_ADD(blas);
            }

            /* Assemble P1loop and related quantities */
            for (int j = 0; j < Nmax; j++) {
                P1loop_0_vv[j] = P13_0_vv[j] + P22_0_vv[j];
                P_CTR_0[j] = kdisc[j] * kdisc[j] * Pbin[j];
                Ptree_0_vv[j] = Pbin[j] * (f * f / 5.);
            }

            /* Computing P_{vd} contribution */

            TIMER_START(m22fill);
            {
                int _N = Nmax + 1;
                #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
                for (int _l = 0; _l < _N; _l++) {
                    int _count_base = _l * _N - _l * (_l - 1) / 2;
                    for (int _i = _l; _i < _N; _i++) {
                        int _c = _count_base + (_i - _l);
                        double complex _nu1 = -0.5 * etam[_i];
                        double complex _nu2 = -0.5 * etam[_l];
                        double complex _nu12 = _nu1 + _nu2;
                        pnlpt->M22_oneline_0_vd_complex[_c] = pnlpt->M22_oneline_complex[_c] * 196. / ((_nu1 * _nu2 * (98. * _nu12 * _nu12 - 14. * _nu12 + 36.) - 91. * _nu12 * _nu12 + 3. * _nu12 + 58.)) * (f * (21. * f * f * (6. + 3. * _nu2 - 10. * _nu2 * _nu2 + 2. * _nu2 * _nu2 * _nu2 + 2. * _nu1 * _nu1 * _nu1 * (1. + 5. * _nu2) + 2 * _nu1 * _nu1 * (-5. - 2. * _nu2 + 10. * _nu2 * _nu2) + _nu1 * (3. - 24. * _nu2 - 4. * _nu2 * _nu2 + 10. * _nu2 * _nu2 * _nu2)) + 14. * f * (18. + 11. * _nu2 + 42. * _nu1 * _nu1 * _nu1 * _nu2 - 31. * _nu2 * _nu2 + _nu1 * _nu1 * (-31. - 22. * _nu2 + 84. * _nu2 * _nu2) + _nu1 * (11. - 74. * _nu2 - 22. * _nu2 * _nu2 + 42. * _nu2 * _nu2 * _nu2)) + 5. * (46. + 13. * _nu2 + 98. * _nu1 * _nu1 * _nu1 * _nu2 - 63. * _nu2 * _nu2 + 7. * _nu1 * _nu1 * (-9. - 10. * _nu2 + 28. * _nu2 * _nu2) + _nu1 * (13. - 138. * _nu2 - 70. * _nu2 * _nu2 + 98. * _nu2 * _nu2 * _nu2)))) / 1470.;

                        if (has_fNL) {
                            pnlpt->M12_oneline_0_vd_complex[_c] = f * (pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1[_c]) + f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd0_f2[_c]);

                            pnlpt->M12_oneline_0_vd_complex_ortho[_c] = f * (pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1_ortho[_c]) + f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd0_f2_ortho[_c]);
                        }
                    }
                }
            }
            TIMER_ADD(m22fill);

            for (index_i = 0; index_i < Nmax + 1; index_i++) {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_0_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * (f * (-35. - 18. * f + 45. * nu1 + 54. * f * nu1)) / 840.;
            }

            /* --- Monopole (ell=0), velocity-density component: P_{0,vd} ~ f --- */
            /* f13/f22/f12 computed inline by batch BLAS macros */

            COMPUTE_P13(0_vd, pnlpt->M13_0_vd_oneline_complex,
                        -1. * Pbin[j] * kdisc[j] * kdisc[j] * sigmav * (2. * f * (625. + 558. * f + 315. * f * f) / 1575.));

            TIMER_START(blas);
            COMPUTE_P22(0_vd, pnlpt->M22_oneline_0_vd_complex);
            COMPUTE_P12(0_vd, pnlpt->M12_oneline_0_vd_complex);
            TIMER_ADD(blas);

            if (has_fNL) {
                TIMER_START(blas);
                COMPUTE_P12(0_vd_ortho, pnlpt->M12_oneline_0_vd_complex_ortho);
                TIMER_ADD(blas);
            }

            for (int j = 0; j < Nmax; j++) {
                P1loop_0_vd[j] = P13_0_vd[j] + P22_0_vd[j];
                Ptree_0_vd[j] = Pbin[j] * (f * 2. / 3.);
            }

            /* Computing P_{dd} contribution */

            TIMER_START(m22fill);
            {
                int _N = Nmax + 1;
                #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
                for (int _l = 0; _l < _N; _l++) {
                    int _count_base = _l * _N - _l * (_l - 1) / 2;
                    for (int _i = _l; _i < _N; _i++) {
                        int _c = _count_base + (_i - _l);
                        double complex _nu1 = -0.5 * etam[_i];
                        double complex _nu2 = -0.5 * etam[_l];
                        double complex _nu12 = _nu1 + _nu2;
                        pnlpt->M22_oneline_0_dd_complex[_c] = pnlpt->M22_oneline_complex[_c] * 196. / ((_nu1 * _nu2 * (98. * _nu12 * _nu12 - 14. * _nu12 + 36.) - 91. * _nu12 * _nu12 + 3. * _nu12 + 58.)) * (98. * f * f * (4. - 2. * _nu2 - 5. * _nu2 * _nu2 + _nu2 * _nu2 * _nu2 + _nu1 * _nu1 * _nu1 * (1. + 3. * _nu2) + _nu1 * _nu1 * (-5. + 2. * _nu2 + 6. * _nu2 * _nu2) + _nu1 * (-2. - 4. * _nu2 + 2. * _nu2 * _nu2 + 3. * _nu2 * _nu2 * _nu2)) + 70. * f * (10. - _nu2 + 14. * _nu1 * _nu1 * _nu1 * _nu2 - 17. * _nu2 * _nu2 + _nu1 * _nu1 * (-17. + 6. * _nu2 + 28. * _nu2 * _nu2) + _nu1 * (-1. - 22. * _nu2 + 6. * _nu2 * _nu2 + 14. * _nu2 * _nu2 * _nu2)) + 15. * (58. + 3. * _nu2 + 98. * _nu1 * _nu1 * _nu1 * _nu2 - 91. * _nu2 * _nu2 + 7. * _nu1 * _nu1 * (-13. - 2. * _nu2 + 28. * _nu2 * _nu2) + _nu1 * (3. - 146. * _nu2 - 14. * _nu2 * _nu2 + 98. * _nu2 * _nu2 * _nu2))) / 2940.;

                        if (has_fNL) {
                            pnlpt->M12_oneline_0_dd_complex[_c] = (pnlpt->M12_oneline_complex_matter_multipoles_dd0_f0[_c]) + f * (pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1[_c]);

                            pnlpt->M12_oneline_0_dd_complex_ortho[_c] = (pnlpt->M12_oneline_complex_matter_multipoles_dd0_f0_ortho[_c]) + f * (pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1_ortho[_c]);
                        }
                    }
                }
            }
            TIMER_ADD(m22fill);
            for (index_i = 0; index_i < Nmax + 1; index_i++) {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_0_dd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] / (1. + 9. * nu1) * (1. + 9. * nu1 + 6. * f * (1. + nu1));
            }

            /* --- Monopole (ell=0), density-density component: P_{0,dd} ~ f^0 --- */
            /* f13/f22/f12 computed inline by batch BLAS macros */

            COMPUTE_P13(0_dd, pnlpt->M13_0_dd_oneline_complex,
                        -1. * Pbin[j] * kdisc[j] * kdisc[j] * sigmav * ((61. - 2. * f + 35. * f * f) / 105.));

            TIMER_START(blas);
            COMPUTE_P22(0_dd, pnlpt->M22_oneline_0_dd_complex);
            COMPUTE_P12(0_dd, pnlpt->M12_oneline_0_dd_complex);
            TIMER_ADD(blas);

            if (has_fNL) {
                TIMER_START(blas);
                COMPUTE_P12(0_dd_ortho, pnlpt->M12_oneline_0_dd_complex_ortho);
                TIMER_ADD(blas);
            }

            for (int j = 0; j < Nmax; j++) {
                P1loop_0_dd[j] = P13_0_dd[j] + P22_0_dd[j];
                Ptree_0_dd[j] = Pbin[j];
            }

            /* Computing P_{vv} contribution - Quadrupole */

            TIMER_START(m22fill);
            {
                int _N = Nmax + 1;
                #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
                for (int _l = 0; _l < _N; _l++) {
                    int _count_base = _l * _N - _l * (_l - 1) / 2;
                    for (int _i = _l; _i < _N; _i++) {
                        int _c = _count_base + (_i - _l);
                        double complex _nu1 = -0.5 * etam[_i];
                        double complex _nu2 = -0.5 * etam[_l];
                        double complex _nu12 = _nu1 + _nu2;

                        pnlpt->M22_oneline_2_vv_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * (f * f * (396. * (50. - 9. * _nu2 + 98. * _nu1 * _nu1 * _nu1 * _nu2 - 35. * _nu2 * _nu2 + 7. * _nu1 * _nu1 * (-5. - 18. * _nu2 + 28. * _nu2 * _nu2) + _nu1 * (-9. - 66. * _nu2 - 126. * _nu2 * _nu2 + 98. * _nu2 * _nu2 * _nu2)) + 231. * f * (142. - 21. * _nu2 + 280. * _nu1 * _nu1 * _nu1 * _nu2 - 106. * _nu2 * _nu2 + 2. * _nu1 * _nu1 * (-53. - 174. * _nu2 + 280. * _nu2 * _nu2) + _nu1 * (-21. - 204. * _nu2 - 348. * _nu2 * _nu2 + 280. * _nu2 * _nu2 * _nu2)) + 49. * f * f * (336. - 62. * _nu2 - 255. * _nu2 * _nu2 + 50. * _nu2 * _nu2 * _nu2 + 10. * _nu1 * _nu1 * _nu1 * (5. + 56. * _nu2) + 5. * _nu1 * _nu1 * (-51. - 142. * _nu2 + 224. * _nu2 * _nu2) + _nu1 * (-62. - 486. * _nu2 - 710. * _nu2 * _nu2 + 560. * _nu2 * _nu2 * _nu2)))) / 135828.;

                        if (has_fNL) {
                            pnlpt->M12_oneline_2_vv_complex[_c] = 20. / 7. * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2[_c]) + f * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv2_f3[_c]);

                            pnlpt->M12_oneline_2_vv_complex_ortho[_c] = 20. / 7. * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2_ortho[_c]) + f * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv2_f3_ortho[_c]);
                        }
                    }
                }
            }
            TIMER_ADD(m22fill);

            /* --- Quadrupole (ell=2), velocity-velocity component: P_{2,vv} ~ f^2 --- */
            /* f13/f22/f12 computed inline by batch BLAS macros */

            for (index_i = 0; index_i < Nmax + 1; index_i++) {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_2_vv_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * (3. * f * f * (-5. + 3. * nu1 + f * (-6. + 5. * nu1))) / 196.;
            }

            COMPUTE_P13(2_vv, pnlpt->M13_2_vv_oneline_complex,
                        -1. * Pbin[j] * kdisc[j] * kdisc[j] * sigmav * (2. * f * f * (54. + 74. * f + 25. * f * f) / 105.));

            TIMER_START(blas);
            COMPUTE_P22(2_vv, pnlpt->M22_oneline_2_vv_complex);
            COMPUTE_P12(2_vv, pnlpt->M12_oneline_2_vv_complex);
            TIMER_ADD(blas);

            if (has_fNL) {
                TIMER_START(blas);
                COMPUTE_P12(2_vv_ortho, pnlpt->M12_oneline_2_vv_complex_ortho);
                TIMER_ADD(blas);
            }

            for (int j = 0; j < Nmax; j++) {
                P1loop_2_vv[j] = P13_2_vv[j] + P22_2_vv[j];
                P_CTR_2[j] = kdisc[j] * kdisc[j] * Pbin[j] * f * 2. / 3.;
                Ptree_2_vv[j] = Pbin[j] * (f * f * 4. / 7.);
            }

            /* Computing P_{vd} contribution - Quadrupole */

            TIMER_START(m22fill);
            {
                int _N = Nmax + 1;
                #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
                for (int _l = 0; _l < _N; _l++) {
                    int _count_base = _l * _N - _l * (_l - 1) / 2;
                    for (int _i = _l; _i < _N; _i++) {
                        int _c = _count_base + (_i - _l);
                        double complex _nu1 = -0.5 * etam[_i];
                        double complex _nu2 = -0.5 * etam[_l];
                        double complex _nu12 = _nu1 + _nu2;

                        pnlpt->M22_oneline_2_vd_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * (f * (7. * f * f * (22. + 11. * _nu2 - 40. * _nu2 * _nu2 + 4. * _nu2 * _nu2 * _nu2 + _nu1 * _nu1 * _nu1 * (4. + 40. * _nu2) + 8. * _nu1 * _nu1 * (-5. - _nu2 + 10. * _nu2 * _nu2) + _nu1 * (11. - 88. * _nu2 - 8. * _nu2 * _nu2 + 40. * _nu2 * _nu2 * _nu2)) + 4. * (46. + 13. * _nu2 + 98. * _nu1 * _nu1 * _nu1 * _nu2 - 63. * _nu2 * _nu2 + 7. * _nu1 * _nu1 * (-9. - 10. * _nu2 + 28. * _nu2 * _nu2) + _nu1 * (13. - 138. * _nu2 - 70. * _nu2 * _nu2 + 98. * _nu2 * _nu2 * _nu2)) + f * (306. + 161. * _nu2 + 672. * _nu1 * _nu1 * _nu1 * _nu2 - 538. * _nu2 * _nu2 + 2. * _nu1 * _nu1 * (-269. - 134. * _nu2 + 672. * _nu2 * _nu2) + _nu1 * (161. - 1196. * _nu2 - 268. * _nu2 * _nu2 + 672. * _nu2 * _nu2 * _nu2)))) / 588.;

                        if (has_fNL) {
                            pnlpt->M12_oneline_2_vd_complex[_c] = 2. * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1[_c]) + f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd2_f2[_c]);

                            pnlpt->M12_oneline_2_vd_complex_ortho[_c] = 2. * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd0_f1_ortho[_c]) + f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd2_f2_ortho[_c]);
                        }
                    }
                }
            }
            TIMER_ADD(m22fill);

            /* --- Quadrupole (ell=2), velocity-density component: P_{2,vd} ~ f --- */
            /* f13/f22/f12 computed inline by batch BLAS macros */

            for (index_i = 0; index_i < Nmax + 1; index_i++) {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_2_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * (f * (-49. - 9. * f + 63. * nu1 + 108. * f * nu1)) / 588.;
            }

            COMPUTE_P13(2_vd, pnlpt->M13_2_vd_oneline_complex,
                        -1. * Pbin[j] * kdisc[j] * kdisc[j] * sigmav * 4. * f * (175. + 180. * f + 126. * f * f) / 441.);
            TIMER_START(blas);
            COMPUTE_P22(2_vd, pnlpt->M22_oneline_2_vd_complex);
            COMPUTE_P12(2_vd, pnlpt->M12_oneline_2_vd_complex);
            TIMER_ADD(blas);
            if (has_fNL) {
                TIMER_START(blas);
                COMPUTE_P12(2_vd_ortho, pnlpt->M12_oneline_2_vd_complex_ortho);
                TIMER_ADD(blas);
            }

            for (int j = 0; j < Nmax; j++) {
                P1loop_2_vd[j] = P13_2_vd[j] + P22_2_vd[j];
                Ptree_2_vd[j] = Pbin[j] * (f * 4. / 3.);
            }

            /* Computing P_{dd} contribution - Quadrupole */

            TIMER_START(m22fill);
            {
                int _N = Nmax + 1;
                #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
                for (int _l = 0; _l < _N; _l++) {
                    int _count_base = _l * _N - _l * (_l - 1) / 2;
                    for (int _i = _l; _i < _N; _i++) {
                        int _c = _count_base + (_i - _l);
                        double complex _nu1 = -0.5 * etam[_i];
                        double complex _nu2 = -0.5 * etam[_l];
                        double complex _nu12 = _nu1 + _nu2;

                        pnlpt->M22_oneline_2_dd_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * f * (4. * (10. - _nu2 + 14. * _nu1 * _nu1 * _nu1 * _nu2 - 17. * _nu2 * _nu2 + _nu1 * _nu1 * (-17. + 6. * _nu2 + 28. * _nu2 * _nu2) + _nu1 * (-1. - 22. * _nu2 + 6. * _nu2 * _nu2 + 14. * _nu2 * _nu2 * _nu2)) + f * (26. - 13. * _nu2 - 37. * _nu2 * _nu2 + 2. * _nu2 * _nu2 * _nu2 + _nu1 * _nu1 * _nu1 * (2. + 24. * _nu2) + _nu1 * _nu1 * (-37. + 22. * _nu2 + 48. * _nu2 * _nu2) + _nu1 * (-13. - 26. * _nu2 + 22. * _nu2 * _nu2 + 24. * _nu2 * _nu2 * _nu2))) / 84.;

                        if (has_fNL) {
                            pnlpt->M12_oneline_2_dd_complex[_c] = 2. * f * (pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1[_c]);

                            pnlpt->M12_oneline_2_dd_complex_ortho[_c] = 2. * f * (pnlpt->M12_oneline_complex_matter_multipoles_dd0_f1_ortho[_c]);
                        }
                    }
                }
            }
            TIMER_ADD(m22fill);

            /* f13/f22/f12 computed inline by batch BLAS macros */

            for (index_i = 0; index_i < Nmax + 1; index_i++) {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_2_dd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * (3. * f * (1. + nu1)) / 28.;
            }

            COMPUTE_P13(2_dd, pnlpt->M13_2_dd_oneline_complex,
                        -1. * Pbin[j] * kdisc[j] * kdisc[j] * sigmav * (2. * f * (35. * f - 2.) / 105.));
            TIMER_START(blas);
            COMPUTE_P22(2_dd, pnlpt->M22_oneline_2_dd_complex);
            COMPUTE_P12(2_dd, pnlpt->M12_oneline_2_dd_complex);
            TIMER_ADD(blas);
            if (has_fNL) {
                TIMER_START(blas);
                COMPUTE_P12(2_dd_ortho, pnlpt->M12_oneline_2_dd_complex_ortho);
                TIMER_ADD(blas);
            }

            for (int j = 0; j < Nmax; j++) {
                P1loop_2_dd[j] = P13_2_dd[j] + P22_2_dd[j];
                Ptree_2_dd[j] = Pbin[j];
            }

            /* Computing P_{vv} contribution - Hexadecapole */

            TIMER_START(m22fill);
            {
                int _N = Nmax + 1;
                #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
                for (int _l = 0; _l < _N; _l++) {
                    int _count_base = _l * _N - _l * (_l - 1) / 2;
                    for (int _i = _l; _i < _N; _i++) {
                        int _c = _count_base + (_i - _l);
                        double complex _nu1 = -0.5 * etam[_i];
                        double complex _nu2 = -0.5 * etam[_l];
                        double complex _nu12 = _nu1 + _nu2;

                        pnlpt->M22_oneline_4_vv_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * (f * f * (1144. * (50. + 98. * _nu1 * _nu1 * _nu1 * _nu2 - _nu2 * (9. + 35. * _nu2) + 7. * _nu1 * _nu1 * (-5. + 2. * _nu2 * (-9. + 14. * _nu2)) + _nu1 * (-9. + 2. * _nu2 * (-33. + 7. * _nu2 * (-9. + 7. * _nu2)))) + 147. * f * f * (483. + 40. * _nu1 * _nu1 * _nu1 * (-1. + 28. * _nu2) - 2. * _nu2 * (-57. + 10. * _nu2 * (29. + 2. * _nu2)) + 20. * _nu1 * _nu1 * (-29. + 2. * _nu2 * (-25. + 56. * _nu2)) + 2. * _nu1 * (57. + 2. * _nu2 * (-327. + 10. * _nu2 * (-25. + 28. * _nu2)))) + 728. * f * (206. + 420. * _nu1 * _nu1 * _nu1 * _nu2 + (7. - 208. * _nu2) * _nu2 + 8. * _nu1 * _nu1 * (-26. + _nu2 * (-53. + 105. * _nu2)) + _nu1 * (7. + 4. * _nu2 * (-108. + _nu2 * (-106. + 105. * _nu2)))))) / 980980.;

                        if (has_fNL) {
                            pnlpt->M12_oneline_4_vv_complex[_c] = 8. / 7. * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2[_c]) + f * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv4_f3[_c]);

                            pnlpt->M12_oneline_4_vv_complex_ortho[_c] = 8. / 7. * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv0_f2_ortho[_c]) + f * f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vv4_f3_ortho[_c]);
                        }
                    }
                }
            }
            TIMER_ADD(m22fill);

            /* --- Hexadecapole (ell=4), velocity-velocity component: P_{4,vv} ~ f^2 --- */
            /* f13/f22/f12 computed inline by batch BLAS macros */

            for (index_i = 0; index_i < Nmax + 1; index_i++) {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_4_vv_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * (3. * f * f * (-55. + 33. * nu1 + f * (-66. + 90. * nu1))) / 5390.;
            }

            COMPUTE_P13(4_vv, pnlpt->M13_4_vv_oneline_complex,
                        -1. * Pbin[j] * kdisc[j] * kdisc[j] * sigmav * (24. * f * f * (33. + 58. * f + 25. * f * f) / 1925.));
            TIMER_START(blas);
            COMPUTE_P22(4_vv, pnlpt->M22_oneline_4_vv_complex);
            COMPUTE_P12(4_vv, pnlpt->M12_oneline_4_vv_complex);
            TIMER_ADD(blas);
            if (has_fNL) {
                TIMER_START(blas);
                COMPUTE_P12(4_vv_ortho, pnlpt->M12_oneline_4_vv_complex_ortho);
                TIMER_ADD(blas);
            }

            for (int j = 0; j < Nmax; j++) {
                P1loop_4_vv[j] = P13_4_vv[j] + P22_4_vv[j];
                P_CTR_4[j] = kdisc[j] * kdisc[j] * Pbin[j] * (f * f * 8. / 35.);
                Ptree_4_vv[j] = Pbin[j] * (f * f * 8. / 35.);
            }

            /* Computing P_{vd} contribution - Hexadecapole */

            TIMER_START(m22fill);
            {
                int _N = Nmax + 1;
                #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
                for (int _l = 0; _l < _N; _l++) {
                    int _count_base = _l * _N - _l * (_l - 1) / 2;
                    for (int _i = _l; _i < _N; _i++) {
                        int _c = _count_base + (_i - _l);
                        double complex _nu1 = -0.5 * etam[_i];
                        double complex _nu2 = -0.5 * etam[_l];
                        double complex _nu12 = _nu1 + _nu2;

                        pnlpt->M22_oneline_4_vd_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * f * f * (14. * f * (26. + 13. * _nu2 - 60. * _nu2 * _nu2 - 8. * _nu2 * _nu2 * _nu2 + _nu1 * _nu1 * _nu1 * (-8. + 60. * _nu2) + 4. * _nu1 * _nu1 * (-15. + 4. * _nu2 + 30. * _nu2 * _nu2) + _nu1 * (13. - 104. * _nu2 + 16. * _nu2 * _nu2 + 60. * _nu2 * _nu2 * _nu2)) + 11. * (58. + 21. * _nu2 + 112. * _nu1 * _nu1 * _nu1 * _nu2 - 106. * _nu2 * _nu2 + 2. * _nu1 * _nu1 * (-53. - 6. * _nu2 + 112. * _nu2 * _nu2) + _nu1 * (21. - 204. * _nu2 - 12. * _nu2 * _nu2 + 112. * _nu2 * _nu2 * _nu2))) / 2695.;

                        if (has_fNL) {
                            pnlpt->M12_oneline_4_vd_complex[_c] = f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd4_f2[_c]);

                            pnlpt->M12_oneline_4_vd_complex_ortho[_c] = f * f * (pnlpt->M12_oneline_complex_matter_multipoles_vd4_f2_ortho[_c]);
                        }
                    }
                }
            }
            TIMER_ADD(m22fill);

            /* f13/f22/f12 computed inline by batch BLAS macros */

            for (index_i = 0; index_i < Nmax + 1; index_i++) {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_4_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 112. / (1. + 9. * nu1) * 9. * (f * f * (1. + 2. * nu1)) / 245.;
            }

            COMPUTE_P13(4_vd, pnlpt->M13_4_vd_oneline_complex,
                        -1. * Pbin[j] * kdisc[j] * kdisc[j] * sigmav * 16. * f * f * (22. + 35. * f) / 1225.);
            TIMER_START(blas);
            COMPUTE_P22(4_vd, pnlpt->M22_oneline_4_vd_complex);
            COMPUTE_P12(4_vd, pnlpt->M12_oneline_4_vd_complex);
            TIMER_ADD(blas);
            if (has_fNL) {
                TIMER_START(blas);
                COMPUTE_P12(4_vd_ortho, pnlpt->M12_oneline_4_vd_complex_ortho);
                TIMER_ADD(blas);
            }

            for (int j = 0; j < Nmax; j++) {
                P1loop_4_vd[j] = P13_4_vd[j] + P22_4_vd[j];
                Ptree_4_vd[j] = Pbin[j] * (f * f * 4. / 7.);
            }

            /* Computing P_{dd} contribution - Hexadecapole */

            TIMER_START(m22fill);
            {
                int _N = Nmax + 1;
                #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
                for (int _l = 0; _l < _N; _l++) {
                    int _count_base = _l * _N - _l * (_l - 1) / 2;
                    for (int _i = _l; _i < _N; _i++) {
                        int _c = _count_base + (_i - _l);
                        double complex _nu1 = -0.5 * etam[_i];
                        double complex _nu2 = -0.5 * etam[_l];
                        double complex _nu12 = _nu1 + _nu2;

                        pnlpt->M22_oneline_4_dd_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * f * f * (2. * _nu1 - 1.) * (2. * _nu2 - 1.) * (1. + _nu12) * (2. + _nu12) / 35.;
                    }
                }
            }
            TIMER_ADD(m22fill);

            double complex *f22_4_dd;
            class_alloc(f22_4_dd, Nmax * sizeof(complex double), pnlpt->error_message);

            /* P22_4_dd - hexadecapole density-density (no P13, no P12) */
            double *f22_tmp = malloc(Nmax * sizeof(double));
            TIMER_START(blas);
            batch_quadratic_form(Np1, Nmax, pnlpt->M22_oneline_4_dd_complex, Xr, Xi,
                                 f22_tmp, Mr_ws, Mi_ws, Y1_ws, Y2_ws);
            TIMER_ADD(blas);
            assemble_P22(Nmax, kdisc, f22_tmp, cutoff, P22_4_dd);
            for (int j = 0; j < Nmax; j++)
                P1loop_4_dd[j] = P22_4_dd[j];
            free(f22_tmp);

        } /* end of IR resummation expression */

        if (irindex == 1) {
            /* Computing FFT for wiggly and non-wiggly parts */

            /* FFTLog decomposition of NW and wiggle spectra */
            TIMER_START(fft);
            double complex *cmsym_nw, *cmsym_nw_transfer, *cmsym_w, *cmsym_w_transfer;
            FFT_TO_CMSYM(Pnw, b, etam, cmsym_nw);
            FFT_TO_CMSYM(Tnw, b_transfer, etam_transfer, cmsym_nw_transfer);
            FFT_TO_CMSYM(Pw, b, etam, cmsym_w);
            FFT_TO_CMSYM(Tw, b_transfer, etam_transfer, cmsym_w_transfer);
            TIMER_ADD(fft);

            /* --- IR resummation: NW and wiggle one-loop integrals ---
             *
             * With IR resummation (irindex=1), we split each one-loop integral
             * into no-wiggle (NW) and wiggle (W) parts:
             *   P_1loop = P_1loop^{NW} + P_1loop^{W} * exp(-Sigma^2 * k^2)
             *
             * The NW part uses P_nw(k) in the FFTLog basis (quadratic form),
             * while the W part uses a bilinear form between P_nw and P_w.
             *
             * mu-power decomposition for Alcock-Paczynski projection:
             * P(k,mu) = sum_{n=0,2,4,6,8} P_mu^{2n}(k) * mu^{2n}
             * These are later integrated over mu with Legendre polynomials
             * and the AP rescaling k -> k(mu), to get P_ell(k). */

            double *P13_mu0_dd, *P13UV_mu0_dd, *P13_mu2_dd, *P13UV_mu2_dd;
            double *P13_mu2_vd, *P13UV_mu2_vd, *P13_mu4_vd, *P13UV_mu4_vd;
            double *P13_mu4_vv, *P13UV_mu4_vv, *P13_mu6, *P13UV_mu6;
            { /* Batch-allocate P13/P13UV IR multipole arrays */
                double **_p13ir[] = {
                    &P13_mu0_dd, &P13UV_mu0_dd, &P13_mu2_dd, &P13UV_mu2_dd,
                    &P13_mu2_vd, &P13UV_mu2_vd, &P13_mu4_vd, &P13UV_mu4_vd,
                    &P13_mu4_vv, &P13UV_mu4_vv, &P13_mu6, &P13UV_mu6};
                for (int ai_ = 0; ai_ < 12; ai_++)
                    class_alloc(*_p13ir[ai_], Nmax * sizeof(double), pnlpt->error_message);
            }
            double *P22_mu0_dd, *P22_mu2_vd, *P22_mu2_dd, *P22_mu4_vv, *P22_mu4_vd;
            double *P22_mu4_dd, *P22_mu6_vv, *P22_mu6_vd, *P22_mu8;
            double *P12_mu0_dd, *P12_mu2_vd, *P12_mu2_dd;
            double *P12_mu4_vv, *P12_mu4_vd, *P12_mu6_vv;
            double *P12_mu0_dd_ortho, *P12_mu2_vd_ortho, *P12_mu2_dd_ortho;
            double *P12_mu4_vv_ortho, *P12_mu4_vd_ortho, *P12_mu6_vv_ortho;
            { /* Batch-allocate P22/P12 NW arrays */
                double **_p22p12[] = {
                    &P22_mu0_dd, &P22_mu2_vd, &P22_mu2_dd, &P22_mu4_vv, &P22_mu4_vd,
                    &P22_mu4_dd, &P22_mu6_vv, &P22_mu6_vd, &P22_mu8,
                    &P12_mu0_dd, &P12_mu2_vd, &P12_mu2_dd,
                    &P12_mu4_vv, &P12_mu4_vd, &P12_mu6_vv,
                    &P12_mu0_dd_ortho, &P12_mu2_vd_ortho, &P12_mu2_dd_ortho,
                    &P12_mu4_vv_ortho, &P12_mu4_vd_ortho, &P12_mu6_vv_ortho};
                for (int ai_ = 0; ai_ < 21; ai_++)
                    *_p22p12[ai_] = malloc(Nmax * sizeof(double));
            }

            /* M22/M12 kernel fill: parallelize over outer index of triangular loop.
             * count = l*(Nmax+1) - l*(l-1)/2 + (i - l) for (l, i) with i >= l. */
            TIMER_START(m22fill);
            {
                int _N = Nmax + 1;
                #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
                for (int _l = 0; _l < _N; _l++) {
                    int _count_base = _l * _N - _l * (_l - 1) / 2;
                    for (int _i = _l; _i < _N; _i++) {
                        int _c = _count_base + (_i - _l);
                        double complex _nu1 = -0.5 * etam[_i];
                        double complex _nu2 = -0.5 * etam[_l];
                        double complex _nu12 = _nu1 + _nu2;

                        pnlpt->M22_oneline_mu2_vd_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * (-1.) * f * (7. * f * (-1. + 2. * _nu1) * (-1. + 2. * _nu2) * (6. + 7. * _nu12) - 4. * (46. + 13. * _nu2 + 98. * _nu1 * _nu1 * _nu1 * _nu2 - 63. * _nu2 * _nu2 + 7. * _nu1 * _nu1 * (-9. - 10. * _nu2 + 28. * _nu2 * _nu2) + _nu1 * (13. - 138. * _nu2 - 70. * _nu2 * _nu2 + 98. * _nu2 * _nu2 * _nu2))) / 392.;

                        pnlpt->M22_oneline_mu2_dd_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * f * (7. * f * (2. + 2. * _nu1 * _nu1 * _nu1 - _nu2 - _nu2 * _nu2 + 2. * _nu2 * _nu2 * _nu2 - _nu1 * _nu1 * (1. + 2. * _nu2) - _nu1 * (1. + 2. * _nu2 + 2. * _nu2 * _nu2)) + 4. * (10. - _nu2 + 14. * _nu1 * _nu1 * _nu1 * _nu2 - 17. * _nu2 * _nu2 + _nu1 * _nu1 * (-17. + 6. * _nu2 + 28. * _nu2 * _nu2) + _nu1 * (-1. - 22. * _nu2 + 6. * _nu2 * _nu2 + 14. * _nu2 * _nu2 * _nu2))) / 56.;

                        pnlpt->M22_oneline_mu4_vv_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * f * f * (147. * f * f * (-1. + 2. * _nu1) * (-1. + 2. * _nu2) - 28. * f * (-1. + 2. * _nu1) * (-1. + 2. * _nu2) * (-2. + 7. * _nu12) + 8. * (50. - 9. * _nu2 + 98 * _nu1 * _nu1 * _nu1 * _nu2 - 35. * _nu2 * _nu2 + 7. * _nu1 * _nu1 * (-5. - 18. * _nu2 + 28 * _nu2 * _nu2) + _nu1 * (-9. - 66. * _nu2 - 126 * _nu2 * _nu2 + 98. * _nu2 * _nu2 * _nu2))) / 1568.;

                        pnlpt->M22_oneline_mu4_vd_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * f * f * (58. + 21. * _nu2 + 112. * _nu1 * _nu1 * _nu1 * _nu2 - 106. * _nu2 * _nu2 + 2. * _nu1 * _nu1 * (-53. - 6. * _nu2 + 112. * _nu2 * _nu2) + 7. * f * (2. + _nu1 + 4. * _nu1 * _nu1 * _nu1 + _nu2 - 8. * _nu1 * _nu2 - 8. * _nu1 * _nu1 * _nu2 - 8. * _nu1 * _nu2 * _nu2 + 4. * _nu2 * _nu2 * _nu2) + _nu1 * (21. - 204. * _nu2 - 12. * _nu2 * _nu2 + 112. * _nu2 * _nu2 * _nu2)) / 56.;

                        pnlpt->M22_oneline_mu4_dd_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * f * f * (2. * _nu1 - 1.) * (2. * _nu2 - 1.) * (2. + _nu1 * _nu1 + 3. * _nu2 + _nu2 * _nu2 + _nu1 * (3. + 2. * _nu2)) / 8.;

                        pnlpt->M22_oneline_mu6_vv_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * f * f * f * (7. * f * (1. + 4. * _nu1 * _nu1 * _nu1 + _nu1 * _nu1 * (2. - 12. * _nu2) + 2. * _nu2 + 2. * _nu2 * _nu2 + 4. * _nu2 * _nu2 * _nu2 - 2. * _nu1 * (-1. + 4. * _nu2 + 6. * _nu2 * _nu2)) + 2. * (26. + 9. * _nu2 + 56. * _nu1 * _nu1 * _nu1 * _nu2 - 38. * _nu2 * _nu2 + 2. * _nu1 * _nu1 * (-19. - 18. * _nu2 + 56. * _nu2 * _nu2) + _nu1 * (9. - 84. * _nu2 - 36. * _nu2 * _nu2 + 56. * _nu2 * _nu2 * _nu2))) / 112.;

                        pnlpt->M22_oneline_mu6_vd_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * f * f * f * (2. * _nu1 - 1.) * (2. * _nu2 - 1.) * (2. + 2. * _nu1 * _nu1 + 5. * _nu2 + 2. * _nu2 * _nu2 + _nu1 * (5. + 4. * _nu2)) / 8.;

                        pnlpt->M22_oneline_mu8_complex[_c] = (pnlpt->M22_oneline_complex[_c] * 196. / (98. * _nu1 * _nu2 * _nu12 * _nu12 - 91. * _nu12 * _nu12 + 36. * _nu1 * _nu2 - 14. * _nu1 * _nu2 * _nu12 + 3. * _nu12 + 58.)) * f * f * f * f * (2. * _nu1 - 1.) * (2. * _nu2 - 1.) * (3. + 4. * _nu1 * _nu1 + 8. * _nu2 + 4. * _nu2 * _nu2 + 8. * _nu1 * (1. + _nu2)) / 32.;

                        if (has_fNL) {
                            pnlpt->M12_oneline_mu2_vd_complex[_c] = f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f1[_c]) + f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2[_c]);
                            pnlpt->M12_oneline_mu2_vd_complex_ortho[_c] = f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f1_ortho[_c]) + f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2_ortho[_c]);
                            pnlpt->M12_oneline_mu2_dd_complex[_c] = f * (pnlpt->M12_oneline_complex_matter_mu_powers_dd2_f1[_c]);
                            pnlpt->M12_oneline_mu2_dd_complex_ortho[_c] = f * (pnlpt->M12_oneline_complex_matter_mu_powers_dd2_f1_ortho[_c]);
                            pnlpt->M12_oneline_mu4_vv_complex[_c] = f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vv4_f2[_c]) + f * f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2[_c]);
                            pnlpt->M12_oneline_mu4_vv_complex_ortho[_c] = f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vv4_f2_ortho[_c]) + f * f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd2_f2_ortho[_c]);
                            pnlpt->M12_oneline_mu4_vd_complex[_c] = f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd4_f2[_c]);
                            pnlpt->M12_oneline_mu4_vd_complex_ortho[_c] = f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vd4_f2_ortho[_c]);
                            pnlpt->M12_oneline_mu6_vv_complex[_c] = f * f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vv6_f3[_c]);
                            pnlpt->M12_oneline_mu6_vv_complex_ortho[_c] = f * f * f * (pnlpt->M12_oneline_complex_matter_mu_powers_vv6_f3_ortho[_c]);
                        }
                    }
                }
            }
            TIMER_ADD(m22fill);

            for (index_i = 0; index_i < Nmax + 1; index_i++) {
                nu1 = -0.5 * etam[index_i];
                pnlpt->M13_mu2_dd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 2. / (1. + 9. * nu1) * 9. * f * (1. + nu1);
                pnlpt->M13_mu2_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 2. / (1. + 9. * nu1) * (-1.) * (f * (7. + 9. * f - 9. * nu1));
                pnlpt->M13_mu4_vv_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] / (1. + 9. * nu1) * (-1.) * 3. * f * f * (5. + 6. * f - 3. * nu1);
                pnlpt->M13_mu4_vd_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 2. / (1. + 9. * nu1) * 9. * f * f * (1. + 2. * nu1);
                pnlpt->M13_mu6_oneline_complex[index_i] = pnlpt->M13_oneline_complex[index_i] * 2. / (1. + 9. * nu1) * 9. * f * f * f * nu1;
            }

            /* f13/f22/f12 for all mu-powers computed inline by batch BLAS macros */

            /* Build X matrix for NW once and compute all P13 mu-powers */
            double *Xr_nw, *Xi_nw;
            class_alloc(Xr_nw, Np1 * Nmax * sizeof(double), pnlpt->error_message);
            class_alloc(Xi_nw, Np1 * Nmax * sizeof(double), pnlpt->error_message);
            TIMER_START(buildx);
            build_X_from_tables(Nmax, Np1, cmsym_nw, kb_cache, sc_cos, sc_sin, Xr_nw, Xi_nw);
            TIMER_ADD(buildx);

            /* P13 NW components: batch_dot_product + UV formula + assemble (parallel) */
            {
                double complex *_m13nw[] = {
                    pnlpt->M13_oneline_complex, pnlpt->M13_mu2_dd_oneline_complex,
                    pnlpt->M13_mu2_vd_oneline_complex, pnlpt->M13_mu4_vv_oneline_complex,
                    pnlpt->M13_mu4_vd_oneline_complex, pnlpt->M13_mu6_oneline_complex};
                double *_p13_out[] = {P13_mu0_dd, P13_mu2_dd, P13_mu2_vd,
                                      P13_mu4_vv, P13_mu4_vd, P13_mu6};
                double *_p13uv[] = {P13UV_mu0_dd, P13UV_mu2_dd, P13UV_mu2_vd,
                                    P13UV_mu4_vv, P13UV_mu4_vd, P13UV_mu6};
                double _uv_coeff[] = {
                    61. / 105.,
                    f * (105. * f - 6.) / 105.,
                    f * (250. + 144. * f) / 105.,
                    f * f * (63. + 48. * f) / 35.,
                    f * f * (44. + 70. * f) / 35.,
                    f * f * f * (46. + 35. * f) / 35.};
                #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
                for (int pi_ = 0; pi_ < 6; pi_++) {
                    int _tid = omp_get_thread_num();
                    double *_f13 = (double *)malloc(Nmax * sizeof(double));
                    batch_dot_product(Np1, Nmax, Xr_nw, Xi_nw, _m13nw[pi_], _f13, ws_dot_mr[_tid], ws_dot_mi[_tid]);
                    for (int j = 0; j < Nmax; j++)
                        _p13uv[pi_][j] = -1. * Pnw[j] * kdisc[j] * kdisc[j] * sigmav * _uv_coeff[pi_];
                    assemble_P13(Nmax, kdisc, _f13, Pnw, _p13uv[pi_], cutoff, _p13_out[pi_]);
                    free(_f13);
                }
            }

            /* Build X matrices for NW transfer once */
            double *Xr_nw_transfer, *Xi_nw_transfer;
            if (has_fNL) {
                class_alloc(Xr_nw_transfer, Np1 * Nmax * sizeof(double), pnlpt->error_message);
                class_alloc(Xi_nw_transfer, Np1 * Nmax * sizeof(double), pnlpt->error_message);
                TIMER_START(buildx);
                build_X_from_tables(Nmax, Np1, cmsym_nw_transfer, kb_transfer_cache, sc_cos_t, sc_sin_t, Xr_nw_transfer, Xi_nw_transfer);
                TIMER_ADD(buildx);
            }

            /* Compute all P22 and P12 NW components in parallel.
             * 9 independent quadratic form computations, each with its own workspace. */
            {
                const double complex *_nw_M22[] = {
                    pnlpt->M22_oneline_complex, pnlpt->M22_oneline_mu2_vd_complex,
                    pnlpt->M22_oneline_mu2_dd_complex, pnlpt->M22_oneline_mu4_vv_complex,
                    pnlpt->M22_oneline_mu4_vd_complex, pnlpt->M22_oneline_mu4_dd_complex,
                    pnlpt->M22_oneline_mu6_vv_complex, pnlpt->M22_oneline_mu6_vd_complex,
                    pnlpt->M22_oneline_mu8_complex};
                double *_nw_P22_out[] = {
                    P22_mu0_dd, P22_mu2_vd, P22_mu2_dd, P22_mu4_vv, P22_mu4_vd,
                    P22_mu4_dd, P22_mu6_vv, P22_mu6_vd, P22_mu8};
                /* M12 matrices: NULL for P22-only entries (indices 5,7,8) */
                const double complex *_nw_M12[] = {
                    pnlpt->M12_oneline_complex, pnlpt->M12_oneline_mu2_vd_complex,
                    pnlpt->M12_oneline_mu2_dd_complex, pnlpt->M12_oneline_mu4_vv_complex,
                    pnlpt->M12_oneline_mu4_vd_complex, NULL,
                    pnlpt->M12_oneline_mu6_vv_complex, NULL, NULL};
                const double complex *_nw_M12_ortho[] = {
                    pnlpt->M12_oneline_complex_ortho, pnlpt->M12_oneline_mu2_vd_complex_ortho,
                    pnlpt->M12_oneline_mu2_dd_complex_ortho, pnlpt->M12_oneline_mu4_vv_complex_ortho,
                    pnlpt->M12_oneline_mu4_vd_complex_ortho, NULL,
                    pnlpt->M12_oneline_mu6_vv_complex_ortho, NULL, NULL};
                double *_nw_P12_out[] = {
                    P12_mu0_dd, P12_mu2_vd, P12_mu2_dd, P12_mu4_vv, P12_mu4_vd,
                    NULL, P12_mu6_vv, NULL, NULL};
                double *_nw_P12_ortho_out[] = {
                    P12_mu0_dd_ortho, P12_mu2_vd_ortho, P12_mu2_dd_ortho,
                    P12_mu4_vv_ortho, P12_mu4_vd_ortho, NULL,
                    P12_mu6_vv_ortho, NULL, NULL};

                #pragma omp parallel for schedule(dynamic)
                for (int _i = 0; _i < 9; _i++) {
                    int _tid = omp_get_thread_num();
                    double *_f22 = (double *)malloc(Nmax * sizeof(double));
                    /* P22: quadratic form with NW X matrices */
                    batch_quadratic_form(Np1, Nmax, _nw_M22[_i], Xr_nw, Xi_nw,
                                         _f22, ws_Mr[_tid], ws_Mi[_tid], ws_Y1[_tid], ws_Y2[_tid]);
                    assemble_P22(Nmax, kdisc, _f22, cutoff, _nw_P22_out[_i]);
                    /* P12: only for entries with fNL and non-NULL M12 */
                    if (has_fNL && _nw_M12[_i] != NULL) {
                        double *_f12 = (double *)malloc(Nmax * sizeof(double));
                        batch_quadratic_form(Np1, Nmax, _nw_M12[_i], Xr_nw_transfer, Xi_nw_transfer,
                                             _f12, ws_Mr[_tid], ws_Mi[_tid], ws_Y1[_tid], ws_Y2[_tid]);
                        assemble_P12(Nmax, kdisc, _f12, Tnw, cutoff, _nw_P12_out[_i]);
                        batch_quadratic_form(Np1, Nmax, _nw_M12_ortho[_i], Xr_nw_transfer, Xi_nw_transfer,
                                             _f12, ws_Mr[_tid], ws_Mi[_tid], ws_Y1[_tid], ws_Y2[_tid]);
                        assemble_P12(Nmax, kdisc, _f12, Tnw, cutoff, _nw_P12_ortho_out[_i]);
                        free(_f12);
                    }
                    free(_f22);
                }
            }
            /* Xr_nw_transfer/Xi_nw_transfer freed later — still needed by W dispatch */

            /* Wiggle (BAO) P13/P22/P12 arrays */
            double *P13_mu0_dd_w, *P22_mu0_dd_w;
            double *P13_mu2_dd_w, *P13_mu2_vd_w, *P13_mu4_vd_w, *P13_mu4_vv_w, *P13_mu6_w;
            double *P22_mu2_vd_w, *P22_mu2_dd_w, *P22_mu4_vv_w, *P22_mu4_vd_w, *P22_mu4_dd_w;
            double *P22_mu6_vv_w, *P22_mu6_vd_w, *P22_mu8_w;
            double *P12_mu0_dd_w, *P12_mu2_vd_w, *P12_mu2_dd_w;
            double *P12_mu4_vv_w, *P12_mu4_vd_w, *P12_mu6_vv_w;
            double *P12_mu0_dd_w_ortho, *P12_mu2_vd_w_ortho, *P12_mu2_dd_w_ortho;
            double *P12_mu4_vv_w_ortho, *P12_mu4_vd_w_ortho, *P12_mu6_vv_w_ortho;
            { /* Batch-allocate wiggle P13/P22/P12 arrays */
                double **_pwiggle[] = {
                    &P13_mu0_dd_w, &P22_mu0_dd_w,
                    &P13_mu2_dd_w, &P13_mu2_vd_w, &P13_mu4_vd_w, &P13_mu4_vv_w, &P13_mu6_w,
                    &P22_mu2_vd_w, &P22_mu2_dd_w, &P22_mu4_vv_w, &P22_mu4_vd_w, &P22_mu4_dd_w,
                    &P22_mu6_vv_w, &P22_mu6_vd_w, &P22_mu8_w,
                    &P12_mu0_dd_w, &P12_mu2_vd_w, &P12_mu2_dd_w,
                    &P12_mu4_vv_w, &P12_mu4_vd_w, &P12_mu6_vv_w,
                    &P12_mu0_dd_w_ortho, &P12_mu2_vd_w_ortho, &P12_mu2_dd_w_ortho,
                    &P12_mu4_vv_w_ortho, &P12_mu4_vd_w_ortho, &P12_mu6_vv_w_ortho};
                for (int ai_ = 0; ai_ < 27; ai_++)
                    *_pwiggle[ai_] = malloc(Nmax * sizeof(double));
            }

            /* Build X matrix for wiggle P13 components */
            double *Xr_w, *Xi_w;
            class_alloc(Xr_w, Np1 * Nmax * sizeof(double), pnlpt->error_message);
            class_alloc(Xi_w, Np1 * Nmax * sizeof(double), pnlpt->error_message);
            TIMER_START(buildx);
            build_X_from_tables(Nmax, Np1, cmsym_w, kb_cache, sc_cos, sc_sin, Xr_w, Xi_w);
            TIMER_ADD(buildx);

            /* Compute P13 wiggle components in parallel */
            {
                double complex *_m13w[] = {
                    pnlpt->M13_oneline_complex, pnlpt->M13_mu2_dd_oneline_complex,
                    pnlpt->M13_mu2_vd_oneline_complex, pnlpt->M13_mu4_vv_oneline_complex,
                    pnlpt->M13_mu4_vd_oneline_complex, pnlpt->M13_mu6_oneline_complex};
                double *_p13w_out[] = {
                    P13_mu0_dd_w, P13_mu2_dd_w, P13_mu2_vd_w,
                    P13_mu4_vv_w, P13_mu4_vd_w, P13_mu6_w};
                #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
                for (int wi_ = 0; wi_ < 6; wi_++) {
                    int _tid = omp_get_thread_num();
                    double *_f13w = (double *)malloc(Nmax * sizeof(double));
                    batch_dot_product(Np1, Nmax, Xr_w, Xi_w, _m13w[wi_], _f13w, ws_dot_mr[_tid], ws_dot_mi[_tid]);
                    for (int j = 0; j < Nmax; j++)
                        _p13w_out[wi_][j] = kdisc[j] * kdisc[j] * kdisc[j] * _f13w[j] * Pnw[j] * exp(-pow(kdisc[j] / cutoff, 6.));
                    free(_f13w);
                }
            }

            /* Build X matrices for wiggle (scaled by 2 for P22/P12) */
            double *Xr_2w, *Xi_2w;
            class_alloc(Xr_2w, Np1 * Nmax * sizeof(double), pnlpt->error_message);
            class_alloc(Xi_2w, Np1 * Nmax * sizeof(double), pnlpt->error_message);
            for (int k = 0; k < Np1 * Nmax; k++) {
                Xr_2w[k] = 2. * Xr_w[k];
                Xi_2w[k] = 2. * Xi_w[k];
            }

            double *Xr_w_transfer, *Xi_w_transfer, *Xr_2w_transfer, *Xi_2w_transfer;
            if (has_fNL) {
                class_alloc(Xr_w_transfer, Np1 * Nmax * sizeof(double), pnlpt->error_message);
                class_alloc(Xi_w_transfer, Np1 * Nmax * sizeof(double), pnlpt->error_message);
                TIMER_START(buildx);
                build_X_from_tables(Nmax, Np1, cmsym_w_transfer, kb_transfer_cache, sc_cos_t, sc_sin_t, Xr_w_transfer, Xi_w_transfer);
                TIMER_ADD(buildx);
                class_alloc(Xr_2w_transfer, Np1 * Nmax * sizeof(double), pnlpt->error_message);
                class_alloc(Xi_2w_transfer, Np1 * Nmax * sizeof(double), pnlpt->error_message);
                for (int k = 0; k < Np1 * Nmax; k++) {
                    Xr_2w_transfer[k] = 2. * Xr_w_transfer[k];
                    Xi_2w_transfer[k] = 2. * Xi_w_transfer[k];
                }
            }

            /* Compute all P22/P12 wiggle components in parallel (bilinear forms: nw x w) */
            {
                const double complex *_w_M22[] = {
                    pnlpt->M22_oneline_complex, pnlpt->M22_oneline_mu2_vd_complex,
                    pnlpt->M22_oneline_mu2_dd_complex, pnlpt->M22_oneline_mu4_vv_complex,
                    pnlpt->M22_oneline_mu4_vd_complex, pnlpt->M22_oneline_mu4_dd_complex,
                    pnlpt->M22_oneline_mu6_vv_complex, pnlpt->M22_oneline_mu6_vd_complex,
                    pnlpt->M22_oneline_mu8_complex};
                double *_w_P22_out[] = {
                    P22_mu0_dd_w, P22_mu2_vd_w, P22_mu2_dd_w, P22_mu4_vv_w, P22_mu4_vd_w,
                    P22_mu4_dd_w, P22_mu6_vv_w, P22_mu6_vd_w, P22_mu8_w};
                const double complex *_w_M12[] = {
                    pnlpt->M12_oneline_complex, pnlpt->M12_oneline_mu2_vd_complex,
                    pnlpt->M12_oneline_mu2_dd_complex, pnlpt->M12_oneline_mu4_vv_complex,
                    pnlpt->M12_oneline_mu4_vd_complex, NULL,
                    pnlpt->M12_oneline_mu6_vv_complex, NULL, NULL};
                const double complex *_w_M12_ortho[] = {
                    pnlpt->M12_oneline_complex_ortho, pnlpt->M12_oneline_mu2_vd_complex_ortho,
                    pnlpt->M12_oneline_mu2_dd_complex_ortho, pnlpt->M12_oneline_mu4_vv_complex_ortho,
                    pnlpt->M12_oneline_mu4_vd_complex_ortho, NULL,
                    pnlpt->M12_oneline_mu6_vv_complex_ortho, NULL, NULL};
                double *_w_P12_out[] = {
                    P12_mu0_dd_w, P12_mu2_vd_w, P12_mu2_dd_w, P12_mu4_vv_w, P12_mu4_vd_w,
                    NULL, P12_mu6_vv_w, NULL, NULL};
                double *_w_P12_ortho_out[] = {
                    P12_mu0_dd_w_ortho, P12_mu2_vd_w_ortho, P12_mu2_dd_w_ortho,
                    P12_mu4_vv_w_ortho, P12_mu4_vd_w_ortho, NULL,
                    P12_mu6_vv_w_ortho, NULL, NULL};

                #pragma omp parallel for schedule(dynamic)
                for (int _i = 0; _i < 9; _i++) {
                    int _tid = omp_get_thread_num();
                    double *_f22w = (double *)malloc(Nmax * sizeof(double));
                    /* P22 wiggle: bilinear form with nw x 2w */
                    batch_bilinear_form(Np1, Nmax, _w_M22[_i], Xr_nw, Xi_nw, Xr_2w, Xi_2w,
                                        _f22w, ws_Mr[_tid], ws_Mi[_tid], ws_Y1[_tid]);
                    assemble_P22(Nmax, kdisc, _f22w, cutoff, _w_P22_out[_i]);
                    /* P12 wiggle: only for entries with fNL and non-NULL M12 */
                    if (has_fNL && _w_M12[_i] != NULL) {
                        double *_f12w = (double *)malloc(Nmax * sizeof(double));
                        batch_bilinear_form(Np1, Nmax, _w_M12[_i], Xr_nw_transfer, Xi_nw_transfer,
                                            Xr_2w_transfer, Xi_2w_transfer, _f12w, ws_Mr[_tid], ws_Mi[_tid], ws_Y1[_tid]);
                        assemble_P12(Nmax, kdisc, _f12w, Tnw, cutoff, _w_P12_out[_i]);
                        batch_bilinear_form(Np1, Nmax, _w_M12_ortho[_i], Xr_nw_transfer, Xi_nw_transfer,
                                            Xr_2w_transfer, Xi_2w_transfer, _f12w, ws_Mr[_tid], ws_Mi[_tid], ws_Y1[_tid]);
                        assemble_P12(Nmax, kdisc, _f12w, Tnw, cutoff, _w_P12_ortho_out[_i]);
                        free(_f12w);
                    }
                    free(_f22w);
                }
            }
            free(Xr_w);
            free(Xi_w);
            free(Xr_2w);
            free(Xi_2w);
            free(Xr_nw);
            free(Xi_nw);
            if (has_fNL) {
                free(Xr_nw_transfer);
                free(Xi_nw_transfer);
                free(Xr_w_transfer);
                free(Xi_w_transfer);
                free(Xr_2w_transfer);
                free(Xi_2w_transfer);
            }

            /* AP spline setup for all NW/wiggle components */
            TIMER_START(ap_spline);
            AP_SPLINE_SETUP(P13_mu4_vv);
            AP_SPLINE_SETUP(P22_mu4_vv);
            AP_SPLINE_SETUP(P22_mu4_vv_w);
            AP_SPLINE_SETUP(P13_mu4_vv_w);
            AP_SPLINE_SETUP(P13_mu6);
            AP_SPLINE_SETUP(P22_mu6_vv);
            AP_SPLINE_SETUP(P22_mu6_vv_w);
            AP_SPLINE_SETUP(P13_mu6_w);
            AP_SPLINE_SETUP(P22_mu8);
            AP_SPLINE_SETUP(P22_mu8_w);
            AP_SPLINE_SETUP(P13_mu0_dd);
            AP_SPLINE_SETUP(P22_mu0_dd);
            AP_SPLINE_SETUP(P13_mu0_dd_w);
            AP_SPLINE_SETUP(P22_mu0_dd_w);
            AP_SPLINE_SETUP(P22_mu2_dd);
            AP_SPLINE_SETUP(P13_mu2_dd);
            AP_SPLINE_SETUP(P22_mu2_dd_w);
            AP_SPLINE_SETUP(P13_mu2_dd_w);
            AP_SPLINE_SETUP(P22_mu4_dd);
            AP_SPLINE_SETUP(P22_mu4_dd_w);
            AP_SPLINE_SETUP(P13_mu2_vd);
            AP_SPLINE_SETUP(P22_mu2_vd);
            AP_SPLINE_SETUP(P22_mu2_vd_w);
            AP_SPLINE_SETUP(P13_mu2_vd_w);
            AP_SPLINE_SETUP(P13_mu4_vd);
            AP_SPLINE_SETUP(P22_mu4_vd);
            AP_SPLINE_SETUP(P22_mu4_vd_w);
            AP_SPLINE_SETUP(P13_mu4_vd_w);
            AP_SPLINE_SETUP(P22_mu6_vd);
            AP_SPLINE_SETUP(P22_mu6_vd_w);
            AP_SPLINE_SETUP(P12_mu6_vv);
            AP_SPLINE_SETUP(P12_mu4_vd);
            AP_SPLINE_SETUP(P12_mu4_vv);
            AP_SPLINE_SETUP(P12_mu2_vd);
            AP_SPLINE_SETUP(P12_mu2_dd);
            AP_SPLINE_SETUP(P12_mu0_dd);
            AP_SPLINE_SETUP(P12_mu0_dd_w);
            AP_SPLINE_SETUP(P12_mu2_vd_w);
            AP_SPLINE_SETUP(P12_mu2_dd_w);
            AP_SPLINE_SETUP(P12_mu4_vv_w);
            AP_SPLINE_SETUP(P12_mu4_vd_w);
            AP_SPLINE_SETUP(P12_mu6_vv_w);
            AP_SPLINE_SETUP(P12_mu6_vv_ortho);
            AP_SPLINE_SETUP(P12_mu4_vd_ortho);
            AP_SPLINE_SETUP(P12_mu4_vv_ortho);
            AP_SPLINE_SETUP(P12_mu2_vd_ortho);
            AP_SPLINE_SETUP(P12_mu2_dd_ortho);
            AP_SPLINE_SETUP(P12_mu0_dd_ortho);
            AP_SPLINE_SETUP(P12_mu0_dd_w_ortho);
            AP_SPLINE_SETUP(P12_mu2_vd_w_ortho);
            AP_SPLINE_SETUP(P12_mu2_dd_w_ortho);
            AP_SPLINE_SETUP(P12_mu4_vv_w_ortho);
            AP_SPLINE_SETUP(P12_mu4_vd_w_ortho);
            AP_SPLINE_SETUP(P12_mu6_vv_w_ortho);
            TIMER_ADD(ap_spline);

            /* Numerical integration over mu */

            double Sigmatot = 0., Exp = 0., mu = 0., ktrue = 0., mutrue = 0.;
            double LegendreP0 = 1., LegendreP2 = 0., LegendreP4 = 0.;
            double P1loopvv = 0., P1loopvd = 0., P1loopdd = 0.;
            double P12vv = 0., P12vd = 0., P12dd = 0.;
            double P12vv_ortho = 0., P12vd_ortho = 0., P12dd_ortho = 0.;
            double p_tree = 0., p_tree_vv = 0., p_tree_vd = 0., p_tree_dd = 0.;
            double P13ratio = 0., P12ratio = 0.;
            double Pctr0 = 0., Pctr2 = 0., Pctr4 = 0.;
            double P1 = 0., P1b1 = 0.;
            double P1loopvv_ap_ir = 0., P1loopvd_ap_ir = 0., P1loopdd_ap_ir = 0.;
            double P12vv_ap_ir = 0., P12vd_ap_ir = 0., P12dd_ap_ir = 0.;
            double P12vv_ap_ir_ortho = 0., P12vd_ap_ir_ortho = 0., P12dd_ap_ir_ortho = 0.;
            int index_gauss2 = 0;
            last_index = 0;

            /* ============================================================
             * ALCOCK-PACZYNSKI PROJECTION VIA GAUSS-LEGENDRE QUADRATURE
             *
             * The AP effect distorts observed (k, mu) when the fiducial
             * cosmology differs from the true one. The mapping is:
             *   k_true = k_fid * sqrt(1/D^2 + (h^2 - 1/D^2) * mu^2)
             *   mu_true = mu * h / sqrt(...)
             * where h = H_fid/H_true, D = D_A,true/D_A,fid.
             *
             * For each k on the FFTLog grid, we integrate over mu using
             * 40-point Gauss-Legendre quadrature. At each mu, we:
             *   1. Compute k_true(mu) and mu_true(mu)
             *   2. Interpolate all P(k) arrays at k_true via splines
             *   3. Build P(k,mu) from mu-power decomposition
             *   4. Project onto Legendre multipoles P_ell(k)
             *
             * The gauss_x[] and gauss_w[] arrays are precomputed and stored
             * in pnlpt during initialization.
             * ============================================================ */

            /* Zero all RSD AP projection arrays */
            { void *_z[] = {
                P1loop_0_vv, P1loop_0_dd, P1loop_0_vd, P1loop_2_vv, P1loop_2_dd, P1loop_2_vd,
                P1loop_4_vv, P1loop_4_dd, P1loop_4_vd,
                P12_0_vv, P12_0_dd, P12_0_vd, P12_2_vv, P12_2_dd, P12_2_vd,
                P12_4_vv, P12_4_dd, P12_4_vd,
                P12_0_vv_ortho, P12_0_dd_ortho, P12_0_vd_ortho,
                P12_2_vv_ortho, P12_2_dd_ortho, P12_2_vd_ortho,
                P12_4_vv_ortho, P12_4_dd_ortho, P12_4_vd_ortho,
                P_CTR_0, P_CTR_2, P_CTR_4,
                Ptree_0_vv, Ptree_0_vd, Ptree_0_dd, Ptree_2_vv, Ptree_2_vd, Ptree_4_vv,
                P10b1, P10, P12
            };
            for (int zi_ = 0; zi_ < 39; zi_++) memset(_z[zi_], 0, sizeof(double) * Nmax); }

            TIMER_START(rsd_ap);
            /* Precompute AP constants outside loops */
            double ap_inv_Dr2_ = 1. / (Dratio * Dratio);
            double ap_hr2_minus_inv_Dr2_ = hratio * hratio - ap_inv_Dr2_;
            double ap_w_hr_Dr2_ = hratio / (Dratio * Dratio);

            for (index_j = Nside; index_j < Nmax - Nside; index_j++) {
                for (index_gauss2 = 0; index_gauss2 < 40; index_gauss2++) {
                    mu = pnlpt->gauss_x[index_gauss2];

                    if (pnlpt->AP_effect == AP_effect_yes) {
                        /* AP rescaling: k_true(mu) and mu_true(mu) from fiducial -> true cosmology */
                        double ap_fac_ = sqrt(ap_inv_Dr2_ + ap_hr2_minus_inv_Dr2_ * mu * mu);
                        mutrue = mu * hratio / ap_fac_;
                        ktrue = kdisc[index_j] * ap_fac_;
                    } else {
                        mutrue = mu;
                        ktrue = kdisc[index_j];
                    }

                    /* Single binary search for all spline interpolations at this ktrue */
                    AP_BSEARCH_SETUP();

                    AP_INTERP_FAST(Pnw);
                    AP_INTERP_FAST(Pw);

                    if (has_fNL) {
                        AP_INTERP_FAST(Tnw);
                        AP_INTERP_FAST(Tw);
                    }

                    AP_INTERP_FAST(P22_mu4_vv);
                    AP_INTERP_FAST(P13_mu4_vv);
                    AP_INTERP_FAST(P22_mu4_vv_w);
                    AP_INTERP_FAST(P13_mu4_vv_w);
                    AP_INTERP_FAST(P13_mu6);
                    AP_INTERP_FAST(P22_mu6_vv);
                    AP_INTERP_FAST(P22_mu6_vv_w);
                    AP_INTERP_FAST(P13_mu6_w);
                    AP_INTERP_FAST(P22_mu8);
                    AP_INTERP_FAST(P22_mu8_w);
                    AP_INTERP_FAST(P22_mu0_dd);
                    AP_INTERP_FAST(P13_mu0_dd);
                    AP_INTERP_FAST(P13_mu0_dd_w);
                    AP_INTERP_FAST(P22_mu0_dd_w);
                    AP_INTERP_FAST(P22_mu2_dd);
                    AP_INTERP_FAST(P13_mu2_dd);
                    AP_INTERP_FAST(P22_mu2_dd_w);
                    AP_INTERP_FAST(P13_mu2_dd_w);
                    AP_INTERP_FAST(P22_mu4_dd);
                    AP_INTERP_FAST(P22_mu4_dd_w);
                    AP_INTERP_FAST(P13_mu2_vd);
                    AP_INTERP_FAST(P22_mu2_vd);
                    AP_INTERP_FAST(P22_mu2_vd_w);
                    AP_INTERP_FAST(P13_mu2_vd_w);
                    AP_INTERP_FAST(P13_mu4_vd);
                    AP_INTERP_FAST(P22_mu4_vd);
                    AP_INTERP_FAST(P22_mu4_vd_w);
                    AP_INTERP_FAST(P13_mu4_vd_w);
                    AP_INTERP_FAST(P22_mu6_vd);
                    AP_INTERP_FAST(P22_mu6_vd_w);

                    if (has_fNL) {
                        AP_INTERP_FAST(P12_mu0_dd);
                        AP_INTERP_FAST(P12_mu2_vd);
                        AP_INTERP_FAST(P12_mu2_dd);
                        AP_INTERP_FAST(P12_mu4_vv);
                        AP_INTERP_FAST(P12_mu4_vd);
                        AP_INTERP_FAST(P12_mu6_vv);
                        AP_INTERP_FAST(P12_mu0_dd_ortho);
                        AP_INTERP_FAST(P12_mu2_vd_ortho);
                        AP_INTERP_FAST(P12_mu2_dd_ortho);
                        AP_INTERP_FAST(P12_mu4_vv_ortho);
                        AP_INTERP_FAST(P12_mu4_vd_ortho);
                        AP_INTERP_FAST(P12_mu6_vv_ortho);
                        AP_INTERP_FAST(P12_mu0_dd_w);
                        AP_INTERP_FAST(P12_mu2_vd_w);
                        AP_INTERP_FAST(P12_mu2_dd_w);
                        AP_INTERP_FAST(P12_mu4_vv_w);
                        AP_INTERP_FAST(P12_mu4_vd_w);
                        AP_INTERP_FAST(P12_mu6_vv_w);
                        AP_INTERP_FAST(P12_mu0_dd_w_ortho);
                        AP_INTERP_FAST(P12_mu2_vd_w_ortho);
                        AP_INTERP_FAST(P12_mu2_dd_w_ortho);
                        AP_INTERP_FAST(P12_mu4_vv_w_ortho);
                        AP_INTERP_FAST(P12_mu4_vd_w_ortho);
                        AP_INTERP_FAST(P12_mu6_vv_w_ortho);
                    }

                    /* Legendre polynomials: P_2(mu) and P_4(mu) */
                    LegendreP2 = (3. * mu * mu - 1.) / 2.;
                    LegendreP4 = (35. * mu * mu * mu * mu - 30. * mu * mu + 3.) / 8.;

                    /* Pre-compute powers of mutrue for RSD mu-expansion */
                    double mu2t = mutrue * mutrue;
                    double mu4t = mu2t * mu2t;
                    double mu6t = mu4t * mu2t;
                    double mu8t = mu4t * mu4t;

                    /* IR resummation damping: Sigma_tot(k, mu) from Senatore & Zaldarriaga */
                    Sigmatot = SigmaBAO * (1. + f * mu2t * (2. + f)) + f * f * mu2t * (mu2t - 1.) * deltaSigmaBAO;
                    Exp = exp(-Sigmatot * ktrue * ktrue);

                    p_tree = (Pnw_ap_out + (1. + Sigmatot * ktrue * ktrue) * Pw_ap_out * Exp);

                    P13ratio = 1. + (Pw_ap_out / Pnw_ap_out) * Exp;

                    if (has_fNL) {
                        P12ratio = 1. + (Tw_ap_out / Tnw_ap_out) * Exp;
                    }

                    P1b1 = (Pnw[index_j] + Pw[index_j] * Exp) * pnlpt->gauss_w[index_gauss2];
                    P1 = (Pnw[index_j] + Pw[index_j] * Exp) * f * pnlpt->gauss_x[index_gauss2] * pnlpt->gauss_x[index_gauss2] * pnlpt->gauss_w[index_gauss2];

                    p_tree_vv = (Pnw_ap_out + (1. + Sigmatot * ktrue * ktrue) * Pw_ap_out * Exp) * f * f * mu4t * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    p_tree_vd = (Pnw_ap_out + (1. + Sigmatot * ktrue * ktrue) * Pw_ap_out * Exp) * 2. * f * mu2t * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    p_tree_dd = (Pnw_ap_out + (1. + Sigmatot * ktrue * ktrue) * Pw_ap_out * Exp) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    Pctr0 = ktrue * ktrue * (Pnw_ap_out + Pw_ap_out * Exp) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    Pctr2 = (Pnw_ap_out + Pw_ap_out * Exp) * pnlpt->gauss_w[index_gauss2] * f * mu2t * ktrue * ktrue * hratio / Dratio / Dratio;
                    Pctr4 = ktrue * ktrue * (Pnw_ap_out + Pw_ap_out * Exp) * pnlpt->gauss_w[index_gauss2] * f * f * mu4t * hratio / Dratio / Dratio;

                    /* P_{1-loop,vv}: mu^4 + mu^6 + mu^8 terms (tree-level is in p_tree_vv) */
                    P1loopvv = ((P13_mu4_vv_ap_out * P13ratio + P22_mu4_vv_ap_out + (P22_mu4_vv_w_ap_out + P13_mu4_vv_w_ap_out) * Exp) * mu4t + (P13_mu6_ap_out * P13ratio + P22_mu6_vv_ap_out + (P22_mu6_vv_w_ap_out + P13_mu6_w_ap_out) * Exp) * mu6t + (P22_mu8_ap_out + P22_mu8_w_ap_out * Exp) * mu8t) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    if (has_fNL) {
                        P12vv = ((P12_mu4_vv_ap_out * P12ratio + P12_mu4_vv_w_ap_out * Exp) * mu4t + (P12_mu6_vv_ap_out * P12ratio + P12_mu6_vv_w_ap_out * Exp) * mu6t) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                        P12vv_ortho = ((P12_mu4_vv_ortho_ap_out * P12ratio + P12_mu4_vv_w_ortho_ap_out * Exp) * mu4t + (P12_mu6_vv_ortho_ap_out * P12ratio + P12_mu6_vv_w_ortho_ap_out * Exp) * mu6t) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    }

                    /* P_{1-loop,dd}: mu^0 + mu^2 + mu^4 terms (tree-level is in p_tree_dd) */
                    P1loopdd = ((P22_mu0_dd_ap_out + P13_mu0_dd_ap_out * P13ratio + (P13_mu0_dd_w_ap_out + P22_mu0_dd_w_ap_out) * Exp) + (P22_mu2_dd_ap_out + P13_mu2_dd_ap_out * P13ratio + (P22_mu2_dd_w_ap_out + P13_mu2_dd_w_ap_out) * Exp) * mu2t + (P22_mu4_dd_ap_out + P22_mu4_dd_w_ap_out * Exp) * mu4t) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    if (has_fNL) {
                        P12dd = ((P12_mu0_dd_ap_out * P12ratio + P12_mu0_dd_w_ap_out * Exp) + (P12_mu2_dd_ap_out * P12ratio + P12_mu2_dd_w_ap_out * Exp) * mu2t) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                        P12dd_ortho = ((P12_mu0_dd_ortho_ap_out * P12ratio + P12_mu0_dd_w_ortho_ap_out * Exp) + (P12_mu2_dd_ortho_ap_out * P12ratio + P12_mu2_dd_w_ortho_ap_out * Exp) * mu2t) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    }

                    /* P_{1-loop,vd}: mu^2 + mu^4 + mu^6 terms (tree-level is in p_tree_vd) */
                    P1loopvd = ((P13_mu2_vd_ap_out * P13ratio + P22_mu2_vd_ap_out + (P22_mu2_vd_w_ap_out + P13_mu2_vd_w_ap_out) * Exp) * mu2t + (P13_mu4_vd_ap_out * P13ratio + P22_mu4_vd_ap_out + (P22_mu4_vd_w_ap_out + P13_mu4_vd_w_ap_out) * Exp) * mu4t + (P22_mu6_vd_ap_out + P22_mu6_vd_w_ap_out * Exp) * mu6t) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    if (has_fNL) {
                        P12vd = ((P12_mu2_vd_ap_out * P12ratio + P12_mu2_vd_w_ap_out * Exp) * mu2t + (P12_mu4_vd_ap_out * P12ratio + P12_mu4_vd_w_ap_out * Exp) * mu4t) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                        P12vd_ortho = ((P12_mu2_vd_ortho_ap_out * P12ratio + P12_mu2_vd_w_ortho_ap_out * Exp) * mu2t + (P12_mu4_vd_ortho_ap_out * P12ratio + P12_mu4_vd_w_ortho_ap_out * Exp) * mu4t) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    }

                    P1loopdd_ap_ir = ((p_tree + P22_mu0_dd_ap_out + P13_mu0_dd_ap_out * P13ratio + (P13_mu0_dd_w_ap_out + P22_mu0_dd_w_ap_out) * Exp) + (P22_mu2_dd_ap_out + P13_mu2_dd_ap_out * P13ratio + (P22_mu2_dd_w_ap_out + P13_mu2_dd_w_ap_out) * Exp) * mu2t + (P22_mu4_dd_ap_out + P22_mu4_dd_w_ap_out * Exp) * mu4t) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    P1loopvd_ap_ir = ((p_tree * 2. * f + P13_mu2_vd_ap_out * P13ratio + P22_mu2_vd_ap_out + (P22_mu2_vd_w_ap_out + P13_mu2_vd_w_ap_out) * Exp) * mu2t + (P13_mu4_vd_ap_out * P13ratio + P22_mu4_vd_ap_out + (P22_mu4_vd_w_ap_out + P13_mu4_vd_w_ap_out) * Exp) * mu4t + (P22_mu6_vd_ap_out + P22_mu6_vd_w_ap_out * Exp) * mu6t) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    LEGENDRE_PROJECT(P1loopvv, P1loop_0_vv, P1loop_2_vv, P1loop_4_vv);
                    P1loop_0_dd[index_j] += P1loopdd * LegendreP0 / 2.;
                    P1loop_2_dd[index_j] += P1loopdd_ap_ir * LegendreP2 * 2.5;
                    P1loop_4_dd[index_j] += P1loopdd_ap_ir * LegendreP4 * 4.5;
                    P1loop_0_vd[index_j] += P1loopvd * LegendreP0 / 2.;
                    P1loop_2_vd[index_j] += P1loopvd * LegendreP2 * 2.5;
                    P1loop_4_vd[index_j] += P1loopvd_ap_ir * LegendreP4 * 4.5;

                    if (has_fNL) {
                        LEGENDRE_PROJECT(P12vv, P12_0_vv, P12_2_vv, P12_4_vv);
                        LEGENDRE_PROJECT(P12dd, P12_0_dd, P12_2_dd, P12_4_dd);
                        LEGENDRE_PROJECT(P12vd, P12_0_vd, P12_2_vd, P12_4_vd);
                        LEGENDRE_PROJECT(P12vv_ortho, P12_0_vv_ortho, P12_2_vv_ortho, P12_4_vv_ortho);
                        LEGENDRE_PROJECT(P12dd_ortho, P12_0_dd_ortho, P12_2_dd_ortho, P12_4_dd_ortho);
                        LEGENDRE_PROJECT(P12vd_ortho, P12_0_vd_ortho, P12_2_vd_ortho, P12_4_vd_ortho);
                    }

                    P_CTR_0[index_j] += Pctr0 * LegendreP0 / 2.;
                    P_CTR_2[index_j] += Pctr2 * LegendreP2 * 2.5;
                    P_CTR_4[index_j] += Pctr4 * LegendreP4 * 4.5;

                    LEGENDRE_PROJECT(p_tree_vv, Ptree_0_vv, Ptree_2_vv, Ptree_4_vv);
                    Ptree_0_vd[index_j] += p_tree_vd * LegendreP0 / 2.;
                    Ptree_0_dd[index_j] += p_tree_dd * LegendreP0 / 2.;
                    Ptree_2_vd[index_j] += p_tree_vd * LegendreP2 * 2.5;

                    P10b1[index_j] += P1b1 * LegendreP0 / 2.;
                    P10[index_j] += P1 * LegendreP0 / 2.;
                    P12[index_j] += P1 * LegendreP2 * 2.5;
                }
            }
            TIMER_ADD(rsd_ap);

            { void *_f[] = {
                dd_P13_mu4_vv, dd_P22_mu4_vv, dd_P22_mu4_vv_w, dd_P13_mu4_vv_w,
                dd_P13_mu6, dd_P22_mu6_vv, dd_P22_mu6_vv_w, dd_P13_mu6_w,
                dd_P22_mu8, dd_P22_mu8_w, dd_P13_mu0_dd, dd_P22_mu0_dd,
                dd_P13_mu0_dd_w, dd_P22_mu0_dd_w, dd_P22_mu2_dd, dd_P13_mu2_dd,
                dd_P22_mu2_dd_w, dd_P13_mu2_dd_w, dd_P22_mu4_dd, dd_P22_mu4_dd_w,
                dd_P13_mu2_vd, dd_P22_mu2_vd, dd_P22_mu2_vd_w, dd_P13_mu2_vd_w,
                dd_P13_mu4_vd, dd_P22_mu4_vd, dd_P22_mu4_vd_w, dd_P13_mu4_vd_w,
                dd_P22_mu6_vd, dd_P22_mu6_vd_w
            };
            for (int fi_ = 0; fi_ < 30; fi_++) free(_f[fi_]); }

            { void *_f[] = {
                P13_mu0_dd_w, P13_mu0_dd, P13UV_mu0_dd, P22_mu0_dd, P22_mu0_dd_w,
                P22_mu8, P22_mu6_vd, P22_mu6_vv, P22_mu4_vv, P22_mu4_vd,
                P22_mu4_dd, P22_mu2_vd, P22_mu2_dd, P13_mu4_vd, P13_mu4_vv,
                P13_mu2_vd, P13_mu2_dd, P13_mu6, P13UV_mu4_vd, P13UV_mu4_vv,
                P13UV_mu2_vd, P13UV_mu2_dd, P13UV_mu6, P22_mu8_w, P22_mu6_vd_w,
                P22_mu6_vv_w, P22_mu4_vv_w, P22_mu4_vd_w, P22_mu4_dd_w, P22_mu2_vd_w,
                P22_mu2_dd_w, P13_mu4_vd_w, P13_mu4_vv_w, P13_mu2_vd_w, P13_mu2_dd_w,
                P13_mu6_w, cmsym_w, cmsym_nw,
            };
            for (int fi_ = 0; fi_ < 38; fi_++) free(_f[fi_]); }

            { void *_f[] = {
                dd_P12_mu6_vv, dd_P12_mu4_vd, dd_P12_mu4_vv, dd_P12_mu2_vd, dd_P12_mu2_dd,
                dd_P12_mu0_dd, dd_P12_mu6_vv_w, dd_P12_mu4_vd_w, dd_P12_mu4_vv_w, dd_P12_mu2_vd_w,
                dd_P12_mu2_dd_w, dd_P12_mu0_dd_w,
                P12_mu6_vv, P12_mu4_vd, P12_mu4_vv, P12_mu2_vd, P12_mu2_dd,
                P12_mu0_dd, P12_mu6_vv_w, P12_mu4_vd_w, P12_mu4_vv_w, P12_mu2_vd_w,
                P12_mu2_dd_w, P12_mu0_dd_w, cmsym_w_transfer, cmsym_nw_transfer,
            };
            for (int fi_ = 0; fi_ < 26; fi_++) free(_f[fi_]); }
            { void *_f[] = {
                dd_P12_mu6_vv_ortho, dd_P12_mu4_vd_ortho, dd_P12_mu4_vv_ortho,
                dd_P12_mu2_vd_ortho, dd_P12_mu2_dd_ortho, dd_P12_mu0_dd_ortho,
                dd_P12_mu6_vv_w_ortho, dd_P12_mu4_vd_w_ortho, dd_P12_mu4_vv_w_ortho,
                dd_P12_mu2_vd_w_ortho, dd_P12_mu2_dd_w_ortho, dd_P12_mu0_dd_w_ortho,
                P12_mu6_vv_ortho, P12_mu4_vd_ortho, P12_mu4_vv_ortho, P12_mu2_vd_ortho, P12_mu2_dd_ortho,
                P12_mu0_dd_ortho, P12_mu6_vv_w_ortho, P12_mu4_vd_w_ortho, P12_mu4_vv_w_ortho, P12_mu2_vd_w_ortho,
                P12_mu2_dd_w_ortho, P12_mu0_dd_w_ortho,
            };
            for (int fi_ = 0; fi_ < 24; fi_++) free(_f[fi_]); }

        } /* end of second IR resummation condition */

        /* Constructing the final output spectra */

        TIMER_START(spline_out);
        SPLINE_INTERP_OUTPUT(P1loop_0_vv, pk_l_0_vv, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * (f * f * (441. + 566. * f + 175. * f * f) / 1225.) + large_for_logs_big, 0);

        SPLINE_INTERP_OUTPUT(P1loop_0_vd, pk_l_0_vd, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * (2. * f * (625. + 558. * f + 315. * f * f) / 1575.) + large_for_logs_big, 0);

        SPLINE_INTERP_OUTPUT(P1loop_0_dd, pk_l_0_dd, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * ((61. - 2. * f + 35. * f * f) / 105.) + large_for_logs_big, 0);

        SPLINE_INTERP_OUTPUT(P1loop_2_vv, pk_l_2_vv, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * (2. * f * f * (54. + 74. * f + 25. * f * f) / 105.) + large_for_logs_big, 0);

        SPLINE_INTERP_OUTPUT(P1loop_2_vd, pk_l_2_vd, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * 4. * f * (175. + 180. * f + 126. * f * f) / 441. + large_for_logs_big, 0);

        SPLINE_INTERP_OUTPUT(P1loop_2_dd, pk_l_2_dd, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * (2. * f * (35. * f - 2.) / 105.) + large_for_logs_big, 0);

        SPLINE_INTERP_OUTPUT(P1loop_4_vv, pk_l_4_vv, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * (24. * f * f * (33. + 58. * f + 25. * f * f) / 1925.) + large_for_logs_big, 0);

        SPLINE_INTERP_OUTPUT(P1loop_4_vd, pk_l_4_vd, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             -1. * exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * sigmav * 16. * f * f * (22. + 35. * f) / 1225. + large_for_logs_big, 0);

        double *ddpk_nl_4_dd;
        class_alloc(ddpk_nl_4_dd, Nmax * sizeof(double), pnlpt->error_message);

        SPLINE_SETUP(kdisc, Nmax, P1loop_4_dd, ddpk_nl_4_dd, _SPLINE_NATURAL_);
        last_index = 0;
        pk_nl_out = 0;
        double ir4dd;
        SPLINE_EVAL(kdisc, Nmax, P1loop_4_dd, ddpk_nl_4_dd, kminnew, &last_index, &ir4dd);
        last_index = 0;
        pk_nl_out = 0;
        for (index_k = 0; index_k < pnlpt->k_size; index_k++) {
            if (pnlpt->k[index_k] <= kmaxnew && pnlpt->k[index_k] >= kminnew) {
                SPLINE_EVAL(kdisc, Nmax, P1loop_4_dd, ddpk_nl_4_dd, pnlpt->k[index_k], &last_index, &pk_nl_out);

                pk_l_4_dd[index_k] = pk_nl_out * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_big;
            } else {
                pk_l_4_dd[index_k] = large_for_logs_big + ir4dd * pow((pnlpt->k[index_k] / kminnew), 4.) * exp(-pow(pnlpt->k[index_k] / 3., 4.));
            }
        }

        { /* fNL P12 RSD multipoles with exponential damping */
            double *_p12_in[] = {
                P12_0_dd, P12_0_dd_ortho, P12_0_vd, P12_0_vd_ortho,
                P12_0_vv, P12_0_vv_ortho, P12_2_dd, P12_2_dd_ortho,
                P12_2_vd, P12_2_vd_ortho, P12_2_vv, P12_2_vv_ortho,
                P12_4_dd, P12_4_dd_ortho, P12_4_vd, P12_4_vd_ortho,
                P12_4_vv, P12_4_vv_ortho};
            double *_p12_out[] = {
                pk_l_fNL_0_dd, pk_l_fNL_0_dd_ortho, pk_l_fNL_0_vd, pk_l_fNL_0_vd_ortho,
                pk_l_fNL_0_vv, pk_l_fNL_0_vv_ortho, pk_l_fNL_2_dd, pk_l_fNL_2_dd_ortho,
                pk_l_fNL_2_vd, pk_l_fNL_2_vd_ortho, pk_l_fNL_2_vv, pk_l_fNL_2_vv_ortho,
                pk_l_fNL_4_dd, pk_l_fNL_4_dd_ortho, pk_l_fNL_4_vd, pk_l_fNL_4_vd_ortho,
                pk_l_fNL_4_vv, pk_l_fNL_4_vv_ortho};
            SPLINE_INTERP_BATCH(_p12_in, _p12_out, 18, kminnew, kmaxnew,
                out_tmp_ * exp(-pow(pnlpt->k[index_k] / 3., 4.)) + large_for_logs_fNL,
                large_for_logs_fNL + epsilon_for_logs_fNL, 0);
        }

        SPLINE_INTERP_OUTPUT(P_CTR_0, pk_CTR_0, kminnew, kmaxnew,
                             (out_tmp_ <= 0 ? 1e-16 : out_tmp_),
                             exp(lnpk_l[index_k] + 2. * lnk_l[index_k]), 0);

        SPLINE_INTERP_OUTPUT(P_CTR_2, pk_CTR_2, kminnew, kmaxnew,
                             (out_tmp_ <= 0 ? 1e-16 : out_tmp_),
                             exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * f * 2. / 3., 0);

        SPLINE_INTERP_OUTPUT(P_CTR_4, pk_CTR_4, kminnew, kmaxnew,
                             (out_tmp_ <= 0 ? 1e-16 : out_tmp_),
                             exp(lnpk_l[index_k] + 2. * lnk_l[index_k]) * f * f * 8. / 35., 0);

        SPLINE_INTERP_OUTPUT(Ptree_0_vv, pk_Tree_0_vv, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             exp(lnpk_l[index_k]) * f * f / 5. + large_for_logs_big, 1);

        SPLINE_INTERP_OUTPUT(Ptree_0_vd, pk_Tree_0_vd, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             exp(lnpk_l[index_k]) * 2. * f / 3. + large_for_logs_big, 1);

        SPLINE_INTERP_OUTPUT(Ptree_0_dd, pk_Tree_0_dd, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             exp(lnpk_l[index_k]) + large_for_logs_big, 1);

        SPLINE_INTERP_OUTPUT(Ptree_2_vv, pk_Tree_2_vv, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             exp(lnpk_l[index_k]) * 4. * f * f / 7. + large_for_logs_big, 1);

        SPLINE_INTERP_OUTPUT(Ptree_2_vd, pk_Tree_2_vd, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             exp(lnpk_l[index_k]) * 4. * f / 3. + large_for_logs_big, 1);

        SPLINE_INTERP_OUTPUT(Ptree_4_vv, pk_Tree_4_vv, kminnew, kmaxnew,
                             out_tmp_ + large_for_logs_big,
                             exp(lnpk_l[index_k]) * 8. * f * f / 35. + large_for_logs_big, 1);
        TIMER_ADD(spline_out);

        { void *_f[] = {
            ddpk_nl_4_dd, P13UV_0_vv, P13_0_vv, P1loop_0_vv, P22_0_vv,
            P13UV_0_vd, P13_0_vd, P1loop_0_vd, P22_0_vd, P13UV_0_dd,
            P13_0_dd, P1loop_0_dd, P22_0_dd, P13UV_2_vv, P13_2_vv,
            P1loop_2_vv, P22_2_vv, P13UV_2_vd, P13_2_vd, P1loop_2_vd,
            P22_2_vd, P13UV_2_dd, P13_2_dd, P1loop_2_dd, P22_2_dd,
            P13UV_4_vv, P13_4_vv, P1loop_4_vv, P22_4_vv, P13UV_4_vd,
            P13_4_vd, P1loop_4_vd, P22_4_vd, P1loop_4_dd, P22_4_dd,
            P12_4_dd, P12_4_vd, P12_4_vv, P12_2_dd, P12_2_vd,
            P12_2_vv, P12_0_dd, P12_0_vd, P12_0_vv, P12_4_dd_ortho,
            P12_4_vd_ortho, P12_4_vv_ortho, P12_2_dd_ortho, P12_2_vd_ortho, P12_2_vv_ortho,
            P12_0_dd_ortho, P12_0_vd_ortho, P12_0_vv_ortho, P_CTR_0, P_CTR_2,
            P_CTR_4,
        };
        for (int fi_ = 0; fi_ < 56; fi_++) free(_f[fi_]); }

    } /* end of RSD conditional expression */

    { void *_f[] = {
        f13, x, x_transfer, x_w, x_w_transfer,
        f22, f12, Pdisc, PPRIMdisc, ddpk_nl,
        ddpk_nl_fNL, ddpk_nl_fNL_ortho, ddpk_CTR, ddpk_Tree, etam,
        cmsym, P22, P13, P12_fNL, P12_fNL_ortho,
        P13UV, P_CTR, P1loop, Ptree, myddlnpk,
    };
    for (int fi_ = 0; fi_ < 25; fi_++) free(_f[fi_]); }

    /* ================================================================== */
    /* === BIASED TRACER POWER SPECTRA === */
    /* ================================================================== */
    /* Compute FFTLog coefficients with bias exponents (b=-1.6, b2=-1.25)
     * for the quadratic bias operator integrals. */

    if (pnlpt->bias == bias_yes) {
        if (pnlpt->nonlinear_pt_verbose > 0)
            printf("Computing the spectra for biased tracers...\n");

        double complex *etam2;
        class_alloc(etam2, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

        double complex *etam2_transfer;
        class_alloc(etam2_transfer, (Nmax + 1) * sizeof(complex double), pnlpt->error_message);

        double b2 = -1.6000001;
        double b2_transfer = -1.25;
        int index_c2 = 0;
        for (index_c2 = 0; index_c2 < Nmax + 1; index_c2++) {
            etam2[index_c2] = b2 + 2. * M_PI * _Complex_I * js[index_c2] / Nmaxd / Delta;
            etam2_transfer[index_c2] = b2_transfer + 2. * M_PI * _Complex_I * js[index_c2] / Nmaxd / Delta;
        }

        /* Precompute k^b2 and k^b2_transfer caches for build_X_matrix_cached */
        double *kb2_cache, *kb2_transfer_cache;
        class_alloc(kb2_cache, Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(kb2_transfer_cache, Nmax * sizeof(double), pnlpt->error_message);
        for (int j = 0; j < Nmax; j++) {
            kb2_cache[j] = pow(kdisc[j], b2);
            kb2_transfer_cache[j] = pow(kdisc[j], b2_transfer);
        }

        /* FFTLog decomposition with bias exponents */
        TIMER_START(fft);
        double complex *cmsym2;
        FFT_TO_CMSYM(Pbin, b2, etam2, cmsym2);
        double complex *cmsym2_transfer;
        FFT_TO_CMSYM(Tbin, b2_transfer, etam2_transfer, cmsym2_transfer);
        TIMER_ADD(fft);

        /* Precompute sincos tables for bias exponents */
        double *sc_cos2, *sc_sin2;     /* sincos table for etam2 */
        double *sc_cos2_t, *sc_sin2_t; /* sincos table for etam2_transfer */
        class_alloc(sc_cos2, Np1 * Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(sc_sin2, Np1 * Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(sc_cos2_t, Np1 * Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(sc_sin2_t, Np1 * Nmax * sizeof(double), pnlpt->error_message);
        TIMER_START(buildx);
        precompute_sincos_table(Nmax, Np1, etam2, lnk_cache, sc_cos2, sc_sin2);
        precompute_sincos_table(Nmax, Np1, etam2_transfer, lnk_cache, sc_cos2_t, sc_sin2_t);
        TIMER_ADD(buildx);

        /* ================================================================== */
        /* === BIAS OPERATOR LOOP INTEGRALS === */
        /* ================================================================== */
        /* Compute one-loop integrals of quadratic bias operators for the
         * galaxy power spectrum model:
         *
         *   P_gg(k) = b1^2 * P_mm(k) + b1*b2 * P_{b1b2}(k)
         *           + b2^2 * I[d2d2](k) + b1*bG2 * P_{b1bG2}(k)
         *           + b2*bG2 * I[d2G2](k) + bG2^2 * I[G2G2](k) + ...
         *
         * The bias operators are (Desjacques, Jeong & Schmidt 2018):
         *   d2 = delta^2 - <delta^2>     (local density squared)
         *   G2 = (nabla_ij Phi)^2 - (nabla^2 Phi)^2  (tidal field squared)
         *   FG2 = F_2[G2]                 (non-local gravity-tidal coupling)
         *
         * Uses a separate FFTLog basis with bias exponents b2=-1.6,
         * b2_transfer=-1.25 (different from the matter b=-0.3) for
         * better convergence of the bias integrals.
         *
         * For fNL (has_fNL==1), P12 corrections couple the primordial
         * transfer function T(k) to the bias operators.
         *
         * See Section 4 of arXiv:2004.10607 for the full model. */

        /* --- I[delta^2 delta^2] --- */

        double complex *f22_Id2d2;
        double *P_Id2d2;
        class_alloc(f22_Id2d2, Nmax * sizeof(complex double), pnlpt->error_message);
        class_alloc(P_Id2d2, Nmax * sizeof(double), pnlpt->error_message);

        double epsilon_for_logs = 1.e-6;

        int count2 = 0;

        /* Build X matrix for eta2 basis (used by Id2d2 and other bias operators) */
        double *Xr2_bias, *Xi2_bias;
        class_alloc(Xr2_bias, Np1 * Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(Xi2_bias, Np1 * Nmax * sizeof(double), pnlpt->error_message);
        TIMER_START(buildx);
        build_X_from_tables(Nmax, Np1, cmsym2, kb2_cache, sc_cos2, sc_sin2, Xr2_bias, Xi2_bias);
        TIMER_ADD(buildx);

        /* Compute Id2d2 using batch quadratic form */
        {
            double *f22_Id2d2_real = malloc(Nmax * sizeof(double));
            TIMER_START(blas);
            batch_quadratic_form(Np1, Nmax, pnlpt->M22basic_oneline_complex, Xr2_bias, Xi2_bias, f22_Id2d2_real, Mr_ws, Mi_ws, Y1_ws, Y2_ws);
            TIMER_ADD(blas);
            for (index_j = 0; index_j < Nmax; index_j++) {
                f22_Id2d2[index_j] = 2. * f22_Id2d2_real[index_j];
                P_Id2d2[index_j] = fabs(kdisc[index_j] * kdisc[index_j] * kdisc[index_j] * f22_Id2d2[index_j] - kdisc[0] * kdisc[0] * kdisc[0] * 2. * f22_Id2d2_real[0] + epsilon_for_logs);
            }
            free(f22_Id2d2_real);
        }

        /*  Having eta2 we compute the rest of the PT matrices */

        TIMER_START(m22fill);
        {
            int _N = Nmax + 1;
            #pragma omp parallel for schedule(static) if(omp_get_max_threads() <= 8)
            for (int _l = 0; _l < _N; _l++) {
                int _count_base = _l * _N - _l * (_l - 1) / 2;
                for (int _i = _l; _i < _N; _i++) {
                    int _c = _count_base + (_i - _l);
                    double complex _ei = etam2[_i];
                    double complex _el = etam2[_l];

                    pnlpt->M_IG2G2[_c] = pnlpt->M22basic_oneline_complex[_c] * ((3. + _ei + _el) * (1. + _ei + _el) / ((-0.5 * _ei) * (-0.5 * _el) * (1. - 0.5 * _ei) * (1. - 0.5 * _el)));

                    pnlpt->M_Id2[_c] = pnlpt->M22basic_oneline_complex[_c] * ((3. + _ei + _el) * (4. + 3.5 * (_ei + _el)) / (14. * (-0.5 * _el) * (-0.5 * _ei)));

                    pnlpt->M_Id2G2[_c] = pnlpt->M22basic_oneline_complex[_c] * ((3. + _ei + _el) / ((-0.5 * _ei) * (-0.5 * _el)));

                    pnlpt->M_IG2[_c] = pnlpt->M22basic_oneline_complex[_c] * (-1. * (3. + _ei + _el) * (1. + _ei + _el) * (6. - 3.5 * (_ei + _el)) / (28. * (1. - 0.5 * _ei) * (1. - 0.5 * _el) * (-0.5 * _ei) * (-0.5 * _el)));

                    if (has_fNL) {
                        pnlpt->M_fNLd2[_c] = pnlpt->M12_oneline_complex_bias_real_space_b2[_c];
                        pnlpt->M_fNLd2_ortho[_c] = pnlpt->M12_oneline_complex_bias_real_space_b2_ortho[_c];
                        pnlpt->M_fNLG2[_c] = pnlpt->M12_oneline_complex_bias_real_space_bG2[_c];
                        pnlpt->M_fNLG2_ortho[_c] = pnlpt->M12_oneline_complex_bias_real_space_bG2_ortho[_c];
                    }
                }
            }
        }
        TIMER_ADD(m22fill);

        /* --- fNL contributions: P12 integrals for delta^2 bias --- */
        double *P_fNLd2 = calloc(Nmax, sizeof(double));
        double *P_fNLd2_ortho = calloc(Nmax, sizeof(double));

        /* Build X matrix for transfer with eta2 basis */
        double *Xr2_transfer, *Xi2_transfer;
        class_alloc(Xr2_transfer, Np1 * Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(Xi2_transfer, Np1 * Nmax * sizeof(double), pnlpt->error_message);
        if (has_fNL) {
            TIMER_START(buildx);
            build_X_from_tables(Nmax, Np1, cmsym2_transfer, kb2_transfer_cache, sc_cos2_t, sc_sin2_t, Xr2_transfer, Xi2_transfer);
            TIMER_ADD(buildx);
        }

        if (has_fNL) {
            TIMER_START(blas);
            COMPUTE_FNL_P12_PAIR(pnlpt->M_fNLd2, pnlpt->M_fNLd2_ortho,
                                 Xr2_transfer, Xi2_transfer, P_fNLd2, P_fNLd2_ortho);
            TIMER_ADD(blas);
        }

        /* fNL contribution from G2 bias */
        double *P_fNLG2 = calloc(Nmax, sizeof(double));
        double *P_fNLG2_ortho = calloc(Nmax, sizeof(double));

        if (has_fNL) {
            TIMER_START(blas);
            COMPUTE_FNL_P12_PAIR(pnlpt->M_fNLG2, pnlpt->M_fNLG2_ortho,
                                 Xr_transfer, Xi_transfer, P_fNLG2, P_fNLG2_ortho);
            TIMER_ADD(blas);
        }

        /* Tidal field P22 integrals: I[d2], I[G2] */
        double *P_Id2 = malloc(Nmax * sizeof(double));
        TIMER_START(blas);
        TIDAL_P22(pnlpt->M_Id2, Xr2_bias, Xi2_bias, P_Id2, IDENTITY);
        TIMER_ADD(blas);

        double *P_IG2 = malloc(Nmax * sizeof(double));
        TIMER_START(blas);
        TIDAL_P22(pnlpt->M_IG2, Xr2_bias, Xi2_bias, P_IG2, fabs);
        TIMER_ADD(blas);

        TIMER_START(spline_out);
        if (has_fNL) { /* fNL bias contributions */
            double *_fnl_in[] = {P_fNLd2, P_fNLd2_ortho, P_fNLG2, P_fNLG2_ortho};
            double *_fnl_out[] = {pk_fNLd2, pk_fNLd2_ortho, pk_fNLG2, pk_fNLG2_ortho};
            SPLINE_INTERP_BATCH(_fnl_in, _fnl_out, 4, kmin, kmax,
                out_tmp_ + large_for_logs_fNL, large_for_logs_fNL + epsilon_for_logs_fNL, 1);
        } else {
            double *_fnl_in[] = {P_fNLd2, P_fNLd2_ortho, P_fNLG2, P_fNLG2_ortho};
            double *_fnl_out[] = {pk_fNLd2, pk_fNLd2_ortho, pk_fNLG2, pk_fNLG2_ortho};
            for (int fi_ = 0; fi_ < 4; fi_++) {
                for (index_k = 0; index_k < pnlpt->k_size; index_k++)
                    _fnl_out[fi_][index_k] = large_for_logs_fNL;
                free(_fnl_in[fi_]);
            }
        }

        SPLINE_INTERP_OUTPUT(P_Id2, pk_Id2, kmin, kmax,
                             large_for_logs_small + out_tmp_,
                             large_for_logs_small, 1);

        SPLINE_INTERP_OUTPUT(P_IG2, pk_IG2, kmin, kmax,
                             out_tmp_ + large_for_logs_small,
                             epsilon_for_logs + large_for_logs_small, 1);
        TIMER_ADD(spline_out);

        /* Tidal field P22 integrals: I[d2 G2], I[G2 G2] */
        double *P_Id2G2 = malloc(Nmax * sizeof(double));
        TIMER_START(blas);
        TIDAL_P22(pnlpt->M_Id2G2, Xr2_bias, Xi2_bias, P_Id2G2, fabs);
        TIMER_ADD(blas);

        double *P_IG2G2 = malloc(Nmax * sizeof(double));
        TIMER_START(blas);
        TIDAL_P22(pnlpt->M_IG2G2, Xr2_bias, Xi2_bias, P_IG2G2, fabs);
        TIMER_ADD(blas);

        /* --- I[F2 G2] (non-local bias from gravity-tidal coupling) --- */

        double *P_IFG2 = malloc(Nmax * sizeof(double));

        /* Build X matrix for IFG2 (uses cmsym2, etam2) */
        double *Xr2, *Xi2;
        class_alloc(Xr2, Np1 * Nmax * sizeof(double), pnlpt->error_message);
        class_alloc(Xi2, Np1 * Nmax * sizeof(double), pnlpt->error_message);
        TIMER_START(buildx);
        build_X_from_tables(Nmax, Np1, cmsym2, kb2_cache, sc_cos2, sc_sin2, Xr2, Xi2);
        TIMER_ADD(buildx);

        /* Compute P_IFG2 using batch operations */
        double *f13_IFG2_tmp = malloc(Nmax * sizeof(double));
        TIMER_START(blas);
        batch_dot_product(Np1, Nmax, Xr2, Xi2, pnlpt->IFG2_oneline_complex, f13_IFG2_tmp, dot_mr_ws, dot_mi_ws);
        TIMER_ADD(blas);
        for (int j = 0; j < Nmax; j++) {
            P_IFG2[j] = fabs(kdisc[j] * kdisc[j] * kdisc[j] * f13_IFG2_tmp[j] * Pbin[j]);
        }
        free(f13_IFG2_tmp);
        free(Xr2);
        free(Xi2);

        /* Note: P_IFG2 is NOT freed here (do_free=0) because the RSD+AP
         * bias section uses it later for AP spline setup via memcpy. */
        TIMER_START(spline_out);
        SPLINE_INTERP_OUTPUT(P_IFG2, pk_IFG2, kmin, kmax,
                             out_tmp_ + large_for_logs_small,
                             epsilon_for_logs + large_for_logs_small, 0);
        TIMER_ADD(spline_out);

        /* --- fNL bias operator integrals in RSD --- */
        if (pnlpt->rsd == rsd_yes) {
            /* Monopole fNL P12 bias contributions (b1b2, b2, b1bG2, bG2) */
            double *P12_0_b1b2 = malloc(Nmax * sizeof(double));
            double *P12_0_b1b2_ortho = malloc(Nmax * sizeof(double));
            FILL_M12_BIAS(pnlpt->M12_0_b1b2_oneline_complex, pnlpt->M12_0_b1b2_oneline_complex_ortho,
                3., pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1, pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1_ortho,
                P12_0_b1b2, P12_0_b1b2_ortho, cmsym2_transfer, kb2_transfer_cache, sc_cos2_t, sc_sin2_t);

            double *P12_0_b2 = malloc(Nmax * sizeof(double));
            double *P12_0_b2_ortho = malloc(Nmax * sizeof(double));
            FILL_M12_BIAS(pnlpt->M12_0_b2_oneline_complex, pnlpt->M12_0_b2_oneline_complex_ortho,
                f, pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1, pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1_ortho,
                P12_0_b2, P12_0_b2_ortho, cmsym2_transfer, kb2_transfer_cache, sc_cos2_t, sc_sin2_t);

            double *P12_0_b1bG2 = malloc(Nmax * sizeof(double));
            double *P12_0_b1bG2_ortho = malloc(Nmax * sizeof(double));
            FILL_M12_BIAS(pnlpt->M12_0_b1bG2_oneline_complex, pnlpt->M12_0_b1bG2_oneline_complex_ortho,
                3., pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1, pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1_ortho,
                P12_0_b1bG2, P12_0_b1bG2_ortho, cmsym_transfer, kb_transfer_cache, sc_cos_t, sc_sin_t);

            double *P12_0_bG2 = malloc(Nmax * sizeof(double));
            double *P12_0_bG2_ortho = malloc(Nmax * sizeof(double));
            FILL_M12_BIAS(pnlpt->M12_0_bG2_oneline_complex, pnlpt->M12_0_bG2_oneline_complex_ortho,
                f, pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1, pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1_ortho,
                P12_0_bG2, P12_0_bG2_ortho, cmsym_transfer, kb_transfer_cache, sc_cos_t, sc_sin_t);

            /* Computing bias P22 corrections in RSD */

            double *P_0_b1b2, *P_0_b2, *P_0_b1bG2, *P_0_bG2;
            double *P_2_b1b2, *P_2_b2, *P_2_b1bG2, *P_2_bG2;
            double *P_4_b2, *P_4_bG2;
            { /* Allocate all P22 bias arrays */
                double **_p22b[] = {
                    &P_0_b1b2, &P_0_b2, &P_0_b1bG2, &P_0_bG2,
                    &P_2_b1b2, &P_2_b2, &P_2_b1bG2, &P_2_bG2,
                    &P_4_b2, &P_4_bG2};
                for (int ai_ = 0; ai_ < 10; ai_++)
                    class_alloc(*_p22b[ai_], Nmax * sizeof(double), pnlpt->error_message);
            }

            double *Xr2, *Xi2;
            class_alloc(Xr2, Np1 * Nmax * sizeof(double), pnlpt->error_message);
            class_alloc(Xi2, Np1 * Nmax * sizeof(double), pnlpt->error_message);
            TIMER_START(buildx);
            build_X_from_tables(Nmax, Np1, cmsym2, kb2_cache, sc_cos2, sc_sin2, Xr2, Xi2);
            TIMER_ADD(buildx);
            double *f22_tmp = malloc(Nmax * sizeof(double));

            /* Monopole (ell=0) M22 bias corrections */
            TIMER_START(m22fill);
            FILL_M22_BIAS(pnlpt->M22_0_b1b2_oneline_complex,
                (-3. + 2. * nu1 + 2. * nu2) * (-12. + 7. * (3. + f) * nu1 + 7. * (3. + f) * nu2) / (42. * nu1 * nu2), P_0_b1b2);
            FILL_M22_BIAS(pnlpt->M22_0_b2_oneline_complex,
                (7. * f * f * (12. + 6. * nu1 * nu1 - 17. * nu2 + 6. * nu2 * nu2 + nu1 * (-17. + 12. * nu2)) + 5. * f * (24. + 14. * nu1 * nu1 - 37. * nu2 + 14. * nu2 * nu2 + nu1 * (-37. + 28. * nu2))) / (210. * nu1 * nu2), P_0_b2);
            FILL_M22_BIAS(pnlpt->M22_0_b1bG2_oneline_complex,
                (-3. + 2. * nu1 + 2. * nu2) * (-1. + 2. * nu1 + 2. * nu2) * (7. * f * (2. + nu1 + nu2) + 3. * (6. + 7. * nu1 + 7. * nu2)) / (42. * nu1 * (1. + nu1) * nu2 * (1. + nu2)), P_0_b1bG2);
            FILL_M22_BIAS(pnlpt->M22_0_bG2_oneline_complex,
                (-3. + 2. * nu1 + 2. * nu2) * (-1. + 2. * nu1 + 2. * nu2) * (-10. * f + 7. * f * (5. * (nu1 + nu2) + f * (-2. + 3. * nu1 + 3. * nu2))) / (210. * nu1 * (1. + nu1) * nu2 * (1. + nu2)), P_0_bG2);
            TIMER_ADD(m22fill);

            /* Quadrupole (ell=2) fNL P12 bias contributions */
            double *P12_2_b2 = malloc(Nmax * sizeof(double));
            double *P12_2_b2_ortho = malloc(Nmax * sizeof(double));
            FILL_M12_BIAS(pnlpt->M12_2_b2_oneline_complex, pnlpt->M12_2_b2_oneline_complex_ortho,
                2. * f, pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1, pnlpt->M12_oneline_complex_bias_multipoles_b2_vv0_f1_ortho,
                P12_2_b2, P12_2_b2_ortho, cmsym2_transfer, kb2_transfer_cache, sc_cos2_t, sc_sin2_t);
            double *P12_2_bG2 = malloc(Nmax * sizeof(double));
            double *P12_2_bG2_ortho = malloc(Nmax * sizeof(double));
            FILL_M12_BIAS(pnlpt->M12_2_bG2_oneline_complex, pnlpt->M12_2_bG2_oneline_complex_ortho,
                2. * f, pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1, pnlpt->M12_oneline_complex_bias_multipoles_bG2_vv0_f1_ortho,
                P12_2_bG2, P12_2_bG2_ortho, cmsym_transfer, kb_transfer_cache, sc_cos_t, sc_sin_t);

            /* Quadrupole (ell=2) M22 bias corrections */
            TIMER_START(m22fill);
            FILL_M22_BIAS(pnlpt->M22_2_b1b2_oneline_complex,
                (-3. + 2. * nu1 + 2. * nu2) * f * (nu1 + nu2) / (3. * nu1 * nu2), P_2_b1b2);
            FILL_M22_BIAS(pnlpt->M22_2_b2_oneline_complex,
                (-3. + 2. * nu1 + 2. * nu2) * f * (-16. + 14. * (nu1 + nu2) + f * (-13. + 12. * (nu1 + nu2))) / (42. * nu1 * nu2), P_2_b2);
            FILL_M22_BIAS(pnlpt->M22_2_b1bG2_oneline_complex,
                (-3. + 2. * nu1 + 2. * nu2) * (-1. + 2. * nu1 + 2. * nu2) * f * (2. + nu1 + nu2) / (3. * nu1 * (1. + nu1) * nu2 * (1. + nu2)), P_2_b1bG2);

            FILL_M22_BIAS(pnlpt->M22_2_bG2_oneline_complex,
                (-3. + 2. * nu1 + 2. * nu2) * (-1. + 2. * nu1 + 2. * nu2) * f * (-2. - f + 7. * (nu1 + nu2) + 6. * f * (nu1 + nu2)) / (21. * nu1 * (1. + nu1) * nu2 * (1. + nu2)), P_2_bG2);
            TIMER_ADD(m22fill);

            /* Hexadecapole (ell=4) M22 bias corrections */
            TIMER_START(m22fill);
            FILL_M22_BIAS(pnlpt->M22_4_b2_oneline_complex,
                (-3. + 2. * nu1 + 2. * nu2) * (-1. + 2. * nu1 + 2. * nu2) * 2. * f * f / (35. * nu1 * nu2), P_4_b2);
            FILL_M22_BIAS(pnlpt->M22_4_bG2_oneline_complex,
                (-3. + 2. * nu1 + 2. * nu2) * (-1. + 2. * nu1 + 2. * nu2) * 4. * f * f * (1. + nu1 + nu2) / (35. * nu1 * (1. + nu1) * nu2 * (1. + nu2)), P_4_bG2);
            TIMER_ADD(m22fill);

            free(f22_tmp);
            free(Xr2);
            free(Xi2);

            /* Numerical integration for bias tracers */
            double *P_Id2d2_2, *P_Id2d2_4, *P_Id2G2_2, *P_Id2G2_4;
            double *P_IG2G2_2, *P_IG2G2_4, *P_4_b1b2, *P_4_b1bG2;
            double *P_IFG2_0b1_x, *P_IFG2_0, *P_IFG2_2;
            double *P_Id2d2_new, *P_Id2G2_new, *P_IG2G2_new, *P_IFG2_new;
            double *P_0_bG2_new, *P_2_bG2_new, *P_4_bG2_new;
            double *P_0_b1b2_new, *P_2_b1b2_new, *P_0_b1bG2_new, *P_2_b1bG2_new;
            double *P_0_b2_new, *P_2_b2_new, *P_4_b2_new;
            double *P12_0_b1b2_new, *P12_0_b1bG2_new, *P12_0_b2_new;
            double *P12_2_b2_new, *P12_0_bG2_new, *P12_2_bG2_new;
            double *P12_0_b1b2_new_ortho, *P12_0_b1bG2_new_ortho, *P12_0_b2_new_ortho;
            double *P12_2_b2_new_ortho, *P12_0_bG2_new_ortho, *P12_2_bG2_new_ortho;
            double *P12_2_b1b2, *P12_2_b1bG2, *P12_4_b1b2, *P12_4_b1bG2;
            double *P12_4_b2, *P12_4_bG2;
            double *P12_2_b1b2_ortho, *P12_2_b1bG2_ortho;
            double *P12_4_b1b2_ortho, *P12_4_b1bG2_ortho;
            double *P12_4_b2_ortho, *P12_4_bG2_ortho;

            { /* Allocate all bias AP arrays */
                double **_bap[] = {
                    &P_Id2d2_2, &P_Id2d2_4, &P_Id2G2_2, &P_Id2G2_4,
                    &P_IG2G2_2, &P_IG2G2_4, &P_4_b1b2, &P_4_b1bG2,
                    &P_IFG2_0b1_x, &P_IFG2_0, &P_IFG2_2,
                    &P_Id2d2_new, &P_Id2G2_new, &P_IG2G2_new, &P_IFG2_new,
                    &P_0_bG2_new, &P_2_bG2_new, &P_4_bG2_new,
                    &P_0_b1b2_new, &P_2_b1b2_new, &P_0_b1bG2_new, &P_2_b1bG2_new,
                    &P_0_b2_new, &P_2_b2_new, &P_4_b2_new,
                    &P12_0_b1b2_new, &P12_0_b1bG2_new, &P12_0_b2_new,
                    &P12_2_b2_new, &P12_0_bG2_new, &P12_2_bG2_new,
                    &P12_0_b1b2_new_ortho, &P12_0_b1bG2_new_ortho, &P12_0_b2_new_ortho,
                    &P12_2_b2_new_ortho, &P12_0_bG2_new_ortho, &P12_2_bG2_new_ortho,
                    &P12_2_b1b2, &P12_2_b1bG2, &P12_4_b1b2, &P12_4_b1bG2,
                    &P12_4_b2, &P12_4_bG2,
                    &P12_2_b1b2_ortho, &P12_2_b1bG2_ortho,
                    &P12_4_b1b2_ortho, &P12_4_b1bG2_ortho,
                    &P12_4_b2_ortho, &P12_4_bG2_ortho};
                for (int ai_ = 0; ai_ < 49; ai_++)
                    class_alloc(*_bap[ai_], Nmax * sizeof(double), pnlpt->error_message);
            }

            { /* Copy bias arrays to _new variants for AP spline setup */
                size_t sz = Nmax * sizeof(double);
                double *_src[] = {
                    P_Id2d2, P_Id2G2, P_IG2G2, P_0_b1b2, P_2_b1b2,
                    P_0_b1bG2, P_2_b1bG2, P_4_b2, P_2_b2, P_0_b2,
                    P_4_bG2, P_2_bG2, P_0_bG2, P_IFG2,
                    P12_0_b1b2, P12_0_b1bG2, P12_0_b2, P12_2_b2, P12_0_bG2, P12_2_bG2,
                    P12_0_b1b2_ortho, P12_0_b1bG2_ortho, P12_0_b2_ortho,
                    P12_2_b2_ortho, P12_0_bG2_ortho, P12_2_bG2_ortho};
                double *_dst[] = {
                    P_Id2d2_new, P_Id2G2_new, P_IG2G2_new, P_0_b1b2_new, P_2_b1b2_new,
                    P_0_b1bG2_new, P_2_b1bG2_new, P_4_b2_new, P_2_b2_new, P_0_b2_new,
                    P_4_bG2_new, P_2_bG2_new, P_0_bG2_new, P_IFG2_new,
                    P12_0_b1b2_new, P12_0_b1bG2_new, P12_0_b2_new, P12_2_b2_new, P12_0_bG2_new, P12_2_bG2_new,
                    P12_0_b1b2_new_ortho, P12_0_b1bG2_new_ortho, P12_0_b2_new_ortho,
                    P12_2_b2_new_ortho, P12_0_bG2_new_ortho, P12_2_bG2_new_ortho};
                for (int ci_ = 0; ci_ < 26; ci_++) memcpy(_dst[ci_], _src[ci_], sz);
            }

            TIMER_START(ap_spline);
            AP_SPLINE_SETUP(P_Id2d2);
            AP_SPLINE_SETUP(P_Id2G2);
            AP_SPLINE_SETUP(P_IG2G2);
            AP_SPLINE_SETUP_EX(P_0_b1b2, P_0_b1b2_new);
            AP_SPLINE_SETUP_EX(P_2_b1b2, P_2_b1b2_new);
            AP_SPLINE_SETUP(P_0_b1bG2);
            AP_SPLINE_SETUP(P_2_b1bG2);
            AP_SPLINE_SETUP_EX(P_0_b2, P_0_b2_new);
            AP_SPLINE_SETUP_EX(P_2_b2, P_2_b2_new);
            AP_SPLINE_SETUP_EX(P_4_b2, P_4_b2_new);
            AP_SPLINE_SETUP_EX(P_0_bG2, P_0_bG2_new);
            AP_SPLINE_SETUP_EX(P_2_bG2, P_2_bG2_new);
            AP_SPLINE_SETUP_EX(P_4_bG2, P_4_bG2_new);
            AP_SPLINE_SETUP_EX(P_IFG2, P_IFG2_new);
            AP_SPLINE_SETUP(Pbin);
            /* P12 bias AP splines + Tbin: only needed for fNL.
             * We declare dd_ and _ap_out variables manually so they exist
             * in function scope for the free block and AP_INTERP macros. */
            double *dd_P12_0_b2 = NULL, *dd_P12_0_b2_ortho = NULL;
            double *dd_P12_0_b1b2 = NULL, *dd_P12_0_b1b2_ortho = NULL;
            double *dd_P12_0_bG2 = NULL, *dd_P12_0_bG2_ortho = NULL;
            double *dd_P12_0_b1bG2 = NULL, *dd_P12_0_b1bG2_ortho = NULL;
            double *dd_P12_2_b2 = NULL, *dd_P12_2_b2_ortho = NULL;
            double *dd_P12_2_bG2 = NULL, *dd_P12_2_bG2_ortho = NULL;
            double *dd_Tbin = NULL;
            double P12_0_b2_ap_out = 0, P12_0_b2_ortho_ap_out = 0;
            double P12_0_b1b2_ap_out = 0, P12_0_b1b2_ortho_ap_out = 0;
            double P12_0_bG2_ap_out = 0, P12_0_bG2_ortho_ap_out = 0;
            double P12_0_b1bG2_ap_out = 0, P12_0_b1bG2_ortho_ap_out = 0;
            double P12_2_b2_ap_out = 0, P12_2_b2_ortho_ap_out = 0;
            double P12_2_bG2_ap_out = 0, P12_2_bG2_ortho_ap_out = 0;
            double Tbin_ap_out = 0;
            if (has_fNL) {
                /* Allocate and compute spline tables for P12 bias arrays */
                double **_dd_p12[] = {
                    &dd_P12_0_b2, &dd_P12_0_b2_ortho, &dd_P12_0_b1b2, &dd_P12_0_b1b2_ortho,
                    &dd_P12_0_bG2, &dd_P12_0_bG2_ortho, &dd_P12_0_b1bG2, &dd_P12_0_b1bG2_ortho,
                    &dd_P12_2_b2, &dd_P12_2_b2_ortho, &dd_P12_2_bG2, &dd_P12_2_bG2_ortho, &dd_Tbin};
                double *_p12_data[] = {
                    P12_0_b2_new, P12_0_b2_new_ortho, P12_0_b1b2_new, P12_0_b1b2_new_ortho,
                    P12_0_bG2_new, P12_0_bG2_new_ortho, P12_0_b1bG2_new, P12_0_b1bG2_new_ortho,
                    P12_2_b2_new, P12_2_b2_new_ortho, P12_2_bG2_new, P12_2_bG2_new_ortho, Tbin};
                for (int si_ = 0; si_ < 13; si_++) {
                    class_alloc(*_dd_p12[si_], sizeof(double) * Nmax, pnlpt->error_message);
                    class_call(array_spline_table_columns(kdisc, Nmax, _p12_data[si_], 1, *_dd_p12[si_],
                                                          _SPLINE_NATURAL_, pnlpt->error_message),
                               pnlpt->error_message, pnlpt->error_message);
                }
            }
            TIMER_ADD(ap_spline);

            double Pd2d2_in = 0., Pd2G2_in = 0., PG2G2_in = 0.;
            double Pb1b2_in = 0., Pb1bG2_in = 0., Pb2_in = 0., PbG2_in = 0.;
            double IFG2_in = 0.;
            double P12b1b2_in = 0., P12b1bG2_in = 0., P12b2_in = 0., P12bG2_in = 0.;
            double P12b1b2_in_ortho = 0., P12b1bG2_in_ortho = 0.;
            double P12b2_in_ortho = 0., P12bG2_in_ortho = 0.;
            double LegendreP0 = 1., LegendreP2 = 0., LegendreP4 = 0.;
            double LegendreP2true = 0., LegendreP4true = 0.;

            int index_gauss2 = 0;
            double mu = 0.;
            double mutrue = 0.;
            double ktrue = 0.;
            last_index = 0;
            double Sigmatot = 0.;
            double p_lo = 0.;
            double Exp = 0.;
            double t_lo = 0.;

            /* Zero all bias AP projection arrays */
            { void *_z[] = {
                P_IFG2_0b1_x, P_IFG2_0, P_IFG2_2, P_Id2d2_2, P_Id2d2_4,
                P_Id2G2_2, P_Id2G2_4, P_IG2G2_2, P_IG2G2_4,
                P_4_b1b2, P_4_b1bG2,
                P12_2_b1b2, P12_2_b1bG2, P12_4_b1b2, P12_4_b1bG2, P12_4_b2, P12_4_bG2,
                P12_2_b1b2_ortho, P12_2_b1bG2_ortho,
                P12_4_b1b2_ortho, P12_4_b1bG2_ortho, P12_4_b2_ortho, P12_4_bG2_ortho,
                P_Id2d2, P_Id2G2, P_IG2G2,
                P_0_b1b2, P_2_b1b2, P_0_b1bG2, P_2_b1bG2,
                P_0_b2, P_2_b2, P_4_b2, P_0_bG2, P_2_bG2, P_4_bG2,
                P12_0_b2, P12_2_b2, P12_0_bG2, P12_2_bG2, P12_0_b1b2, P12_0_b1bG2,
                P12_0_b2_ortho, P12_2_b2_ortho, P12_0_bG2_ortho, P12_2_bG2_ortho,
                P12_0_b1b2_ortho, P12_0_b1bG2_ortho
            };
            for (int zi_ = 0; zi_ < 48; zi_++) memset(_z[zi_], 0, sizeof(double) * Nmax); }

            TIMER_START(bias_ap);
            /* Precompute AP constants (reuse from RSD loop if already set) */
            double bias_ap_inv_Dr2_ = 1. / (Dratio * Dratio);
            double bias_ap_hr2_minus_inv_Dr2_ = hratio * hratio - bias_ap_inv_Dr2_;

            for (index_j = Nside; index_j < Nmax - Nside; index_j++) {
                Pnw_ap_out = 0.;
                Pw_ap_out = 0.;
                Tnw_ap_out = 0.;
                Tw_ap_out = 0.;

                for (index_gauss2 = 0; index_gauss2 < 40; index_gauss2++) {
                    mu = pnlpt->gauss_x[index_gauss2];

                    if (pnlpt->AP_effect == AP_effect_yes) {
                        double ap_fac_ = sqrt(bias_ap_inv_Dr2_ + bias_ap_hr2_minus_inv_Dr2_ * mu * mu);
                        mutrue = mu * hratio / ap_fac_;
                        ktrue = kdisc[index_j] * ap_fac_;
                    } else {
                        mutrue = mu;
                        ktrue = kdisc[index_j];
                    }

                    LegendreP2 = (3. * mu * mu - 1.) / 2.;
                    LegendreP4 = (35. * mu * mu * mu * mu - 30. * mu * mu + 3.) / 8.;
                    double mu2t = mutrue * mutrue;
                    double mu4t = mu2t * mu2t;
                    LegendreP2true = (3. * mu2t - 1.) / 2.;
                    LegendreP4true = (35. * mu4t - 30. * mu2t + 3.) / 8.;

                    /* Single binary search for all spline interpolations at this ktrue */
                    AP_BSEARCH_SETUP();

                    AP_INTERP_FAST(Pnw);
                    AP_INTERP_FAST(Pw);
                    AP_INTERP_FAST(Pbin);

                    if (has_fNL) {
                        AP_INTERP_FAST(Tnw);
                        AP_INTERP_FAST(Tw);
                        AP_INTERP_FAST(Tbin);
                        AP_INTERP_EX_FAST(P12_0_b1b2, P12_0_b1b2_new);
                        AP_INTERP_EX_FAST(P12_0_b1bG2, P12_0_b1bG2_new);
                        AP_INTERP_EX_FAST(P12_0_b2, P12_0_b2_new);
                        AP_INTERP_EX_FAST(P12_2_b2, P12_2_b2_new);
                        AP_INTERP_EX_FAST(P12_0_bG2, P12_0_bG2_new);
                        AP_INTERP_EX_FAST(P12_2_bG2, P12_2_bG2_new);
                        AP_INTERP_EX_FAST(P12_0_b1b2_ortho, P12_0_b1b2_new_ortho);
                        AP_INTERP_EX_FAST(P12_0_b1bG2_ortho, P12_0_b1bG2_new_ortho);
                        AP_INTERP_EX_FAST(P12_0_b2_ortho, P12_0_b2_new_ortho);
                        AP_INTERP_EX_FAST(P12_2_b2_ortho, P12_2_b2_new_ortho);
                        AP_INTERP_EX_FAST(P12_0_bG2_ortho, P12_0_bG2_new_ortho);
                        AP_INTERP_EX_FAST(P12_2_bG2_ortho, P12_2_bG2_new_ortho);
                    }

                    AP_INTERP_EX_FAST(P_IFG2, P_IFG2_new);
                    AP_INTERP_EX_FAST(P_Id2d2, P_Id2d2_new);
                    AP_INTERP_EX_FAST(P_Id2G2, P_Id2G2_new);
                    AP_INTERP_EX_FAST(P_IG2G2, P_IG2G2_new);
                    AP_INTERP_EX_FAST(P_0_b1b2, P_0_b1b2_new);
                    AP_INTERP_EX_FAST(P_2_b1b2, P_2_b1b2_new);
                    AP_INTERP_EX_FAST(P_0_b1bG2, P_0_b1bG2_new);
                    AP_INTERP_EX_FAST(P_2_b1bG2, P_2_b1bG2_new);
                    AP_INTERP_EX_FAST(P_0_b2, P_0_b2_new);
                    AP_INTERP_EX_FAST(P_2_b2, P_2_b2_new);
                    AP_INTERP_EX_FAST(P_4_b2, P_4_b2_new);
                    AP_INTERP_EX_FAST(P_0_bG2, P_0_bG2_new);
                    AP_INTERP_EX_FAST(P_2_bG2, P_2_bG2_new);
                    AP_INTERP_EX_FAST(P_4_bG2, P_4_bG2_new);

                    Sigmatot = SigmaBAO * (1. + f * mu2t * (2. + f)) + f * f * mu2t * (mu2t - 1.) * deltaSigmaBAO;
                    Exp = exp(-Sigmatot * ktrue * ktrue);
                    p_lo = (Pnw_ap_out + Pw_ap_out * Exp) / Pbin_ap_out;

                    IFG2_in = p_lo * P_IFG2_ap_out * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    Pd2d2_in = P_Id2d2_ap_out * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    Pd2G2_in = P_Id2G2_ap_out * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    PG2G2_in = P_IG2G2_ap_out * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    Pb1b2_in = (P_0_b1b2_ap_out + LegendreP2true * P_2_b1b2_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    Pb1bG2_in = (P_0_b1bG2_ap_out + LegendreP2true * P_2_b1bG2_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    Pb2_in = (P_0_b2_ap_out + LegendreP2true * P_2_b2_ap_out + LegendreP4true * P_4_b2_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    PbG2_in = (P_0_bG2_ap_out + LegendreP2true * P_2_bG2_ap_out + LegendreP4true * P_4_bG2_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                    if (has_fNL) {
                        t_lo = (Tnw_ap_out + Tw_ap_out * Exp) / Tbin_ap_out;

                        P12b1b2_in = t_lo * (P12_0_b1b2_ap_out)*pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                        P12b1bG2_in = t_lo * (P12_0_b1bG2_ap_out)*pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                        P12b2_in = t_lo * (P12_0_b2_ap_out + LegendreP2true * P12_2_b2_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                        P12bG2_in = t_lo * (P12_0_bG2_ap_out + LegendreP2true * P12_2_bG2_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;

                        P12b1b2_in_ortho = t_lo * (P12_0_b1b2_ortho_ap_out)*pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                        P12b1bG2_in_ortho = t_lo * (P12_0_b1bG2_ortho_ap_out)*pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                        P12b2_in_ortho = t_lo * (P12_0_b2_ortho_ap_out + LegendreP2true * P12_2_b2_ortho_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                        P12bG2_in_ortho = t_lo * (P12_0_bG2_ortho_ap_out + LegendreP2true * P12_2_bG2_ortho_ap_out) * pnlpt->gauss_w[index_gauss2] * hratio / Dratio / Dratio;
                    }

                    P_IFG2_0b1_x[index_j] += IFG2_in * LegendreP0 / 2.;
                    P_IFG2_0[index_j] += IFG2_in * f * mu2t * LegendreP0 / 2.;
                    P_IFG2_2[index_j] += IFG2_in * f * mu2t * LegendreP2 * 2.5;
                    LEGENDRE_PROJECT(Pd2d2_in, P_Id2d2, P_Id2d2_2, P_Id2d2_4);
                    LEGENDRE_PROJECT(Pd2G2_in, P_Id2G2, P_Id2G2_2, P_Id2G2_4);
                    LEGENDRE_PROJECT(PG2G2_in, P_IG2G2, P_IG2G2_2, P_IG2G2_4);
                    LEGENDRE_PROJECT(Pb1b2_in, P_0_b1b2, P_2_b1b2, P_4_b1b2);
                    LEGENDRE_PROJECT(Pb1bG2_in, P_0_b1bG2, P_2_b1bG2, P_4_b1bG2);
                    LEGENDRE_PROJECT(Pb2_in, P_0_b2, P_2_b2, P_4_b2);
                    LEGENDRE_PROJECT(PbG2_in, P_0_bG2, P_2_bG2, P_4_bG2);

                    if (has_fNL) {
                        LEGENDRE_PROJECT(P12b2_in, P12_0_b2, P12_2_b2, P12_4_b2);
                        LEGENDRE_PROJECT(P12bG2_in, P12_0_bG2, P12_2_bG2, P12_4_bG2);
                        LEGENDRE_PROJECT(P12b1b2_in, P12_0_b1b2, P12_2_b1b2, P12_4_b1b2);
                        LEGENDRE_PROJECT(P12b1bG2_in, P12_0_b1bG2, P12_2_b1bG2, P12_4_b1bG2);
                        LEGENDRE_PROJECT(P12b2_in_ortho, P12_0_b2_ortho, P12_2_b2_ortho, P12_4_b2_ortho);
                        LEGENDRE_PROJECT(P12bG2_in_ortho, P12_0_bG2_ortho, P12_2_bG2_ortho, P12_4_bG2_ortho);
                        LEGENDRE_PROJECT(P12b1b2_in_ortho, P12_0_b1b2_ortho, P12_2_b1b2_ortho, P12_4_b1b2_ortho);
                        LEGENDRE_PROJECT(P12b1bG2_in_ortho, P12_0_b1bG2_ortho, P12_2_b1bG2_ortho, P12_4_b1bG2_ortho);
                    }
                }
            }
            TIMER_ADD(bias_ap);

            { void *_f[] = {
                P_Id2d2_new, P_Id2G2_new, P_IG2G2_new, P_0_b1b2_new, P_2_b1b2_new,
                P_0_b1bG2_new, P_2_b1bG2_new, P_0_b2_new, P_2_b2_new, P_4_b2_new,
                P_0_bG2_new, P_2_bG2_new, P_4_bG2_new, dd_P_Id2d2, dd_P_Id2G2,
                dd_P_IG2G2, dd_P_0_b1b2, dd_P_2_b1b2, dd_P_0_b1bG2, dd_P_2_b1bG2,
                dd_P_0_b2, dd_P_2_b2, dd_P_4_b2, dd_P_0_bG2, dd_P_2_bG2,
                dd_P_4_bG2, dd_P_IFG2, dd_Pbin, dd_Tbin, P12_0_b1b2_new,
                dd_P12_0_b1b2, P12_0_b1bG2_new, dd_P12_0_b1bG2, P12_0_b2_new, dd_P12_0_b2,
                P12_2_b2_new, dd_P12_2_b2, P12_0_bG2_new, dd_P12_0_bG2, P12_2_bG2_new,
                dd_P12_2_bG2, P12_0_b1b2_new_ortho, dd_P12_0_b1b2_ortho, P12_0_b1bG2_new_ortho, dd_P12_0_b1bG2_ortho,
                P12_0_b2_new_ortho, dd_P12_0_b2_ortho, P12_2_b2_new_ortho, dd_P12_2_b2_ortho, P12_0_bG2_new_ortho,
                dd_P12_0_bG2_ortho, P12_2_bG2_new_ortho, dd_P12_2_bG2_ortho,
            };
            for (int fi_ = 0; fi_ < 53; fi_++) free(_f[fi_]); }

            TIMER_START(spline_out);
            if (has_fNL) { /* fNL P12 biased tracer multipoles */
                double *_p12b_in[] = {
                    P12_0_b1b2, P12_0_b1b2_ortho, P12_0_b1bG2, P12_0_b1bG2_ortho,
                    P12_0_b2, P12_0_b2_ortho, P12_0_bG2, P12_0_bG2_ortho,
                    P12_2_b2, P12_2_b2_ortho, P12_2_bG2, P12_2_bG2_ortho,
                    P12_4_b2, P12_4_b2_ortho, P12_4_bG2, P12_4_bG2_ortho,
                    P12_2_b1b2, P12_2_b1b2_ortho, P12_4_b1b2, P12_4_b1b2_ortho,
                    P12_2_b1bG2, P12_2_b1bG2_ortho, P12_4_b1bG2, P12_4_b1bG2_ortho};
                double *_p12b_out[] = {
                    pk12_l_0_b1b2, pk12_l_0_b1b2_ortho, pk12_l_0_b1bG2, pk12_l_0_b1bG2_ortho,
                    pk12_l_0_b2, pk12_l_0_b2_ortho, pk12_l_0_bG2, pk12_l_0_bG2_ortho,
                    pk12_l_2_b2, pk12_l_2_b2_ortho, pk12_l_2_bG2, pk12_l_2_bG2_ortho,
                    pk12_l_4_b2, pk12_l_4_b2_ortho, pk12_l_4_bG2, pk12_l_4_bG2_ortho,
                    pk12_l_2_b1b2, pk12_l_2_b1b2_ortho, pk12_l_4_b1b2, pk12_l_4_b1b2_ortho,
                    pk12_l_2_b1bG2, pk12_l_2_b1bG2_ortho, pk12_l_4_b1bG2, pk12_l_4_b1bG2_ortho};
                SPLINE_INTERP_BATCH(_p12b_in, _p12b_out, 24, kminnew, kmaxnew,
                    out_tmp_ + large_for_logs_fNL, large_for_logs_fNL + epsilon_for_logs_fNL, 1);
            } else { /* No fNL: fill with offset constant and free input arrays */
                double *_p12b_in[] = {
                    P12_0_b1b2, P12_0_b1b2_ortho, P12_0_b1bG2, P12_0_b1bG2_ortho,
                    P12_0_b2, P12_0_b2_ortho, P12_0_bG2, P12_0_bG2_ortho,
                    P12_2_b2, P12_2_b2_ortho, P12_2_bG2, P12_2_bG2_ortho,
                    P12_4_b2, P12_4_b2_ortho, P12_4_bG2, P12_4_bG2_ortho,
                    P12_2_b1b2, P12_2_b1b2_ortho, P12_4_b1b2, P12_4_b1b2_ortho,
                    P12_2_b1bG2, P12_2_b1bG2_ortho, P12_4_b1bG2, P12_4_b1bG2_ortho};
                double *_p12b_out[] = {
                    pk12_l_0_b1b2, pk12_l_0_b1b2_ortho, pk12_l_0_b1bG2, pk12_l_0_b1bG2_ortho,
                    pk12_l_0_b2, pk12_l_0_b2_ortho, pk12_l_0_bG2, pk12_l_0_bG2_ortho,
                    pk12_l_2_b2, pk12_l_2_b2_ortho, pk12_l_2_bG2, pk12_l_2_bG2_ortho,
                    pk12_l_4_b2, pk12_l_4_b2_ortho, pk12_l_4_bG2, pk12_l_4_bG2_ortho,
                    pk12_l_2_b1b2, pk12_l_2_b1b2_ortho, pk12_l_4_b1b2, pk12_l_4_b1b2_ortho,
                    pk12_l_2_b1bG2, pk12_l_2_b1bG2_ortho, pk12_l_4_b1bG2, pk12_l_4_b1bG2_ortho};
                for (int fi_ = 0; fi_ < 24; fi_++) {
                    for (index_k = 0; index_k < pnlpt->k_size; index_k++)
                        _p12b_out[fi_][index_k] = large_for_logs_fNL;
                    free(_p12b_in[fi_]); /* match do_free=1 behavior */
                }
            }
            { /* Id2d2/Id2G2/IG2G2 higher multipoles */
                double *_id_in[] = {P_Id2d2_2, P_Id2d2_4, P_Id2G2_2, P_Id2G2_4, P_IG2G2_2, P_IG2G2_4};
                double *_id_out[] = {pk_Id2d2_2, pk_Id2d2_4, pk_Id2G2_2, pk_Id2G2_4, pk_IG2G2_2, pk_IG2G2_4};
                SPLINE_INTERP_BATCH(_id_in, _id_out, 6, kminnew, kmaxnew,
                    out_tmp_ + large_for_logs_big, large_for_logs_big, 1);
            }
            { /* Biased tracer P22 multipoles */
                double *_bt_in[] = {
                    P_0_b1b2, P_0_b2, P_0_b1bG2, P_2_b1b2, P_4_b1b2, P_0_bG2,
                    P_2_b2, P_2_b1bG2, P_4_b1bG2, P_2_bG2, P_4_b2, P_4_bG2};
                double *_bt_out[] = {
                    pk_l_0_b1b2, pk_l_0_b2, pk_l_0_b1bG2, pk_l_2_b1b2, pk_l_4_b1b2, pk_l_0_bG2,
                    pk_l_2_b2, pk_l_2_b1bG2, pk_l_4_b1bG2, pk_l_2_bG2, pk_l_4_b2, pk_l_4_bG2};
                SPLINE_INTERP_BATCH(_bt_in, _bt_out, 12, kminnew, kmaxnew,
                    out_tmp_ + large_for_logs_big, epsilon_for_logs + large_for_logs_big, 1);
            }
            free(P_IFG2_new);
            { /* IFG2 multipoles */
                double *_ifg_in[] = {P_IFG2_0b1_x, P_IFG2_0, P_IFG2_2};
                double *_ifg_out[] = {pk_IFG2_0b1, pk_IFG2_0, pk_IFG2_2};
                SPLINE_INTERP_BATCH(_ifg_in, _ifg_out, 3, kminnew, kmaxnew,
                    out_tmp_ + large_for_logs_big, large_for_logs_big, 0);
            }
            TIMER_ADD(spline_out);
            free(P_IFG2_0b1_x);
            free(P_IFG2_0);
            free(P_IFG2_2);

            /* P22 bias coefficient arrays (P_0_b1b2, ..., P_4_bG2) are freed
             * by SPLINE_INTERP_BATCH above (do_free=1), so no explicit free here. */

        } /* end of RSD conditional expression */

        TIMER_START(spline_out);
        { /* Id2d2/Id2G2/IG2G2 monopole */
            double *_id0_in[] = {P_Id2d2, P_Id2G2, P_IG2G2};
            double *_id0_out[] = {pk_Id2d2, pk_Id2G2, pk_IG2G2};
            SPLINE_INTERP_BATCH(_id0_in, _id0_out, 3, kminnew, kmaxnew,
                out_tmp_ + large_for_logs_big, epsilon_for_logs + large_for_logs_big, 1);
        }
        TIMER_ADD(spline_out);
        free(Xr2_bias);
        free(Xi2_bias);
        free(Xr2_transfer);
        free(Xi2_transfer);
        free(etam2);
        free(cmsym2);
        free(etam2_transfer);
        free(cmsym2_transfer);
        free(P_IFG2);

    } /* end of bias conditional expression */

    else
    {
        if (pnlpt->nonlinear_pt_verbose > 0)
            printf("No bias tracers requested.\n");

        double epsilon_for_logs = 1.e-6;
        double *_eps_arr[] = {pk_Id2d2, pk_Id2, pk_IG2, pk_Id2G2, pk_IG2G2, pk_IFG2};
        for (int ai_ = 0; ai_ < 6; ai_++)
            for (index_k = 0; index_k < pnlpt->k_size; index_k++)
                _eps_arr[ai_][index_k] = epsilon_for_logs;
    }
    /* --- Free temporary arrays --- */
    { void *_f[] = {
        js, kdisc, Pbin, Tbin,
        Pw, Pnw, dd_Pnw, dd_Pw,
        Tw, Tnw, dd_Tnw, dd_Tw, P10b1,
        P10, P12,
    };
    for (int fi_ = 0; fi_ < 15; fi_++) free(_f[fi_]); }

    /* Free per-thread workspace and restore BLAS threads */
    for (int _t = 0; _t < _n_threads; _t++) {
        free(ws_Mr[_t]); free(ws_Mi[_t]);
        free(ws_Y1[_t]); free(ws_Y2[_t]);
        free(ws_dot_mr[_t]); free(ws_dot_mi[_t]);
    }
    free(ws_Mr); free(ws_Mi); free(ws_Y1); free(ws_Y2);
    free(ws_dot_mr); free(ws_dot_mi);
    openblas_set_num_threads(_blas_threads_saved);

    clock_gettime(CLOCK_MONOTONIC, &_ts_total1);
    double _t_total = (_ts_total1.tv_sec - _ts_total0.tv_sec)
                    + 1e-9 * (_ts_total1.tv_nsec - _ts_total0.tv_nsec);
    {
        FILE *_tf = fopen("/home/groups/ophilcox/CLASS-PT/pt_timing.log", "w");
        if (_tf) {
            fprintf(_tf, "=== nonlinear_pt_loop timing (%.3fs total) ===\n", _t_total);
            fprintf(_tf, "  BLAS:          %.3fs (%.0f%%)\n", _t_blas, 100.*_t_blas/_t_total);
            fprintf(_tf, "  FFT:           %.3fs (%.0f%%)\n", _t_fft, 100.*_t_fft/_t_total);
            fprintf(_tf, "  build_X:       %.3fs (%.0f%%)\n", _t_buildx, 100.*_t_buildx/_t_total);
            fprintf(_tf, "  IR resum:      %.3fs (%.0f%%)\n", _t_ir, 100.*_t_ir/_t_total);
            fprintf(_tf, "  AP spline:     %.3fs (%.0f%%)\n", _t_ap_spline, 100.*_t_ap_spline/_t_total);
            fprintf(_tf, "  RSD AP loop:   %.3fs (%.0f%%)\n", _t_rsd_ap, 100.*_t_rsd_ap/_t_total);
            fprintf(_tf, "  Bias AP loop:  %.3fs (%.0f%%)\n", _t_bias_ap, 100.*_t_bias_ap/_t_total);
            fprintf(_tf, "  Spline out:    %.3fs (%.0f%%)\n", _t_spline_out, 100.*_t_spline_out/_t_total);
            fprintf(_tf, "  M22 fill:      %.3fs (%.0f%%)\n", _t_m22fill, 100.*_t_m22fill/_t_total);
            double _t_tagged = _t_blas + _t_fft + _t_buildx + _t_ir + _t_ap_spline
                             + _t_rsd_ap + _t_bias_ap + _t_spline_out + _t_m22fill;
            double _t_other = _t_total - _t_tagged;
            fprintf(_tf, "  Other:         %.3fs (%.0f%%)\n", _t_other, 100.*_t_other/_t_total);
            fclose(_tf);
        }
    }

    return _SUCCESS_;
}
