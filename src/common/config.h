#pragma once

#ifdef MKL_FOUND
#include <mkl_cblas.h>
#define ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS
#define ARMA_USE_MKL_TYPES
#else

#if defined(__APPLE__)
#include "vecLib/cblas.h"
#elif defined(HAVE_FLEXIBLAS_CBLAS_H)
#include "flexiblas/cblas.h"
#elif defined(HAVE_OPENBLAS_CBLAS_H)
#include "openblas/cblas.h"
#else
#include "cblas.h"
#endif
#endif

// #if !defined(ARMA_64BIT_WORD)
// #define ARMA_64BIT_WORD
#ifndef ARMA_DONT_PRINT_FAST_MATH_WARNING
#define ARMA_DONT_PRINT_FAST_MATH_WARNING
#endif
#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_BLAS
#define ARMA_USE_LAPACK
// #endif
