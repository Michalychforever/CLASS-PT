/** @file fft.c Documented fast fourier transform module
 *
 * Nils Schoenberg, 12.10.2017
 *
 * This module computes the fast fourier transform (FFT) of any function
 */
 #include "fft.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
//Mathematical constants used in this file
 #define MATH_PI_2 6.2831853071795864769252867665590057683943387987502
 #define INV_SQRT_2 0.70710678118654752440084436210484903928483593768847
//Function implementations
 /**
 * Computes the fast fourier transform of some arbitrary input array input_real and input_imag of sizes N
 * Returns the output by writing into output_real and output_imag
 *
 * It is assumed that all the arrays are allocated and of size N
 * N must be a power of 2
 *
 * For recursion there is a stepsize parameter
 * If the full FFT should be calculated, set
 *
 *          ** stepsize = 1 **
 *
 * Uses iterative radix-2 DIT with precomputed twiddle table for
 * better cache locality and fewer trig calls (N/2 vs N*log2N).
 * */
void FFT(double* input_real, double* input_imag, double* output_real, double* output_imag, int N, int stepsize){
	if (N <= 1) {
		output_real[0] = input_real[0];
		output_imag[0] = input_imag[0];
		return;
	}
	if (N == 2) {
		output_real[0] = input_real[0] + input_real[stepsize];
		output_real[1] = input_real[0] - input_real[stepsize];
		output_imag[0] = input_imag[0] + input_imag[stepsize];
		output_imag[1] = input_imag[0] - input_imag[stepsize];
		return;
	}

	/* Compute log2(N) */
	int log2N = 0;
	{ int tmp = N; while (tmp > 1) { tmp >>= 1; log2N++; } }

	/* Bit-reversal permutation: copy strided input to output in bit-reversed order */
	for (int i = 0; i < N; i++) {
		int rev = 0, x = i;
		for (int b = 0; b < log2N; b++) {
			rev = (rev << 1) | (x & 1);
			x >>= 1;
		}
		output_real[rev] = input_real[i * stepsize];
		output_imag[rev] = input_imag[i * stepsize];
	}

	/* Stages 1-3 (butterfly sizes 2, 4, 8): hardcoded twiddle factors */

	/* Stage 1: butterflies of size 2, twiddle = 1 */
	for (int k = 0; k < N; k += 2) {
		double t0r = output_real[k], t0i = output_imag[k];
		double t1r = output_real[k+1], t1i = output_imag[k+1];
		output_real[k]   = t0r + t1r;
		output_imag[k]   = t0i + t1i;
		output_real[k+1] = t0r - t1r;
		output_imag[k+1] = t0i - t1i;
	}

	/* Stage 2: butterflies of size 4, twiddles = {1, -i} */
	for (int k = 0; k < N; k += 4) {
		/* j=0: twiddle = 1 */
		double er = output_real[k], ei = output_imag[k];
		double odr = output_real[k+2], odi = output_imag[k+2];
		output_real[k]   = er + odr;
		output_imag[k]   = ei + odi;
		output_real[k+2] = er - odr;
		output_imag[k+2] = ei - odi;
		/* j=1: twiddle = cos(pi/2) + i*sin(pi/2) → (0, 1) → tr = 0*or + 1*oi, ti = 0*oi - 1*or */
		er = output_real[k+1]; ei = output_imag[k+1];
		odr = output_real[k+3]; odi = output_imag[k+3];
		output_real[k+1] = er + odi;
		output_imag[k+1] = ei - odr;
		output_real[k+3] = er - odi;
		output_imag[k+3] = ei + odr;
	}

	if (log2N >= 3) {
		/* Stage 3: butterflies of size 8, twiddles at multiples of pi/4 */
		for (int k = 0; k < N; k += 8) {
			/* j=0: twiddle = 1 */
			double er, ei, odr, odi, tr, ti;
			er = output_real[k]; ei = output_imag[k];
			odr = output_real[k+4]; odi = output_imag[k+4];
			output_real[k]   = er + odr;
			output_imag[k]   = ei + odi;
			output_real[k+4] = er - odr;
			output_imag[k+4] = ei - odi;
			/* j=1: twiddle = (INV_SQRT_2, INV_SQRT_2) */
			er = output_real[k+1]; ei = output_imag[k+1];
			odr = output_real[k+5]; odi = output_imag[k+5];
			tr = INV_SQRT_2*odr + INV_SQRT_2*odi;
			ti = INV_SQRT_2*odi - INV_SQRT_2*odr;
			output_real[k+1] = er + tr;
			output_imag[k+1] = ei + ti;
			output_real[k+5] = er - tr;
			output_imag[k+5] = ei - ti;
			/* j=2: twiddle = (0, 1) */
			er = output_real[k+2]; ei = output_imag[k+2];
			odr = output_real[k+6]; odi = output_imag[k+6];
			output_real[k+2] = er + odi;
			output_imag[k+2] = ei - odr;
			output_real[k+6] = er - odi;
			output_imag[k+6] = ei + odr;
			/* j=3: twiddle = (-INV_SQRT_2, INV_SQRT_2) */
			er = output_real[k+3]; ei = output_imag[k+3];
			odr = output_real[k+7]; odi = output_imag[k+7];
			tr = -INV_SQRT_2*odr + INV_SQRT_2*odi;
			ti = -INV_SQRT_2*odr - INV_SQRT_2*odi;
			output_real[k+3] = er + tr;
			output_imag[k+3] = ei + ti;
			output_real[k+7] = er - tr;
			output_imag[k+7] = ei - ti;
		}
	}

	/* Remaining stages (s >= 4): use precomputed twiddle table */
	if (log2N >= 4) {
		int half_N = N >> 1;
		/* Precompute twiddle factors for largest butterfly.
		 * For stage s with m=2^s, twiddle[j] = tw[j * (N/m)] */
		double *tw_cos = (double *)malloc(half_N * sizeof(double));
		double *tw_sin = (double *)malloc(half_N * sizeof(double));
		double base_angle = MATH_PI_2 / (double)N;
		for (int i = 0; i < half_N; i++) {
			sincos(i * base_angle, &tw_sin[i], &tw_cos[i]);
		}

		for (int s = 4; s <= log2N; s++) {
			int m = 1 << s;
			int m2 = m >> 1;
			int tw_stride = half_N >> (s - 1);

			for (int k = 0; k < N; k += m) {
				for (int j = 0; j < m2; j++) {
					int tw_idx = j * tw_stride;
					double cos_val = tw_cos[tw_idx];
					double sin_val = tw_sin[tw_idx];

					int idx_e = k + j;
					int idx_o = k + j + m2;
					double er = output_real[idx_e], ei = output_imag[idx_e];
					double odr = output_real[idx_o], odi = output_imag[idx_o];

					double tr = cos_val * odr + sin_val * odi;
					double ti = cos_val * odi - sin_val * odr;

					output_real[idx_e] = er + tr;
					output_imag[idx_e] = ei + ti;
					output_real[idx_o] = er - tr;
					output_imag[idx_o] = ei - ti;
				}
			}
		}

		free(tw_cos);
		free(tw_sin);
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
  double dct_step = -0.25*MATH_PI_2/N;
  for(i=0;i<N;++i){
    sincos(i*dct_step, &s, &c);
    temp = output_real[i];
    output_real[i] = c*temp-s*output_imag[i];
    output_imag[i] = c*output_imag[i]+s*temp;
  }
}
