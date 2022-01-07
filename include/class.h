#ifndef __CLASS__
#define __CLASS__

/* standard libraries */
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "float.h"
#include <complex.h>
#undef I

#include <unistd.h>
#include <stdio.h>
#include <limits.h>

#ifdef _OPENMP
#include "omp.h"
#endif

/* tools for class */
#include "quadrature.h"
#include "growTable.h"
#include "arrays.h"
#include "dei_rkck.h"
#include "parser.h"

/* class modules */
#include "common.h"
#include "input.h"
#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "primordial.h"
#include "fourier.h"
#include "nonlinear_pt.h"
#include "transfer.h"
#include "harmonic.h"
#include "distortions.h"
#include "lensing.h"
#include "output.h"

#endif
