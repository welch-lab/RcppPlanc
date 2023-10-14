#pragma once
#include "config.h"
/* Copyright 2016 Ramakrishnan Kannan */
// utility functions

// #ifndef _VERBOSE
// #define _VERBOSE 1
// #endif

enum algotype { MU, HALS, ANLSBPP, NAIVEANLSBPP, AOADMM,
        NESTEROV, CPALS, GNSYM, R2, PGD, PGNCG };

enum normtype { NONE, L2NORM, MAXNORM };

enum helptype { NMF, DISTNMF, NTF, DISTNTF, JOINTNMF, DISTJOINTNMF, HIERNMF };


#include <armadillo>
#include <cmath>
#include <iostream>
#include <vector>

// using namespace std;

#ifndef ERR
#define ERR std::cerr
#endif

#ifndef WARN
#define WARN std::cerr
#endif

#ifndef INFO
#define INFO std::cout
#endif

#ifndef OUTPUT
#define OUTPUT std::cout
#endif

constexpr auto EPSILON_1EMINUS16 = 0.00000000000000001;
constexpr auto EPSILON_1EMINUS8=0.00000001;
constexpr auto EPSILON = 0.000001;
constexpr auto EPSILON_1EMINUS12 = 1e-12;
constexpr auto NUMBEROF_DECIMAL_PLACES = 12;
constexpr auto RAND_SEED = 100;
constexpr auto RAND_SEED_SPARSE = 100;
constexpr auto WTRUE_SEED=1196089;
constexpr auto HTRUE_SEED=1230587;


#define PRINTMATINFO(A) "::" #A "::" << (A).n_rows << "x" << (A).n_cols

#define PRINTMAT(A) PRINTMATINFO((A)) << std::endl << (A)

typedef std::vector<int> STDVEC;
typedef uint64_t ULONG;

void absmat(const arma::fmat *X);

inline void tic();
inline double toc();

int random_sieve(const int);

template <typename FVT>
inline void fillVector(const FVT value, std::vector<FVT> *a) {
  for (unsigned int ii = 0; ii < a->size(); ii++) {
    (*a)[ii] = value;
  }
}
