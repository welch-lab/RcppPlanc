#pragma once
/* Copyright 2016 Ramakrishnan Kannan */
// utility functions

#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif //_OPENMP
#include <cstdio>
#include <chrono>
#include <ctime>
#include <stack>
#include <typeinfo>
#include <utility>
#include <vector>
#include <random>
#include "data.hpp"
#include "config.h"

#include "utils.h"

extern "C" {
  #include "hw_detect.h"
}

static uint64_t powersof10[16] = {1,
                                  10,
                                  100,
                                  1000,
                                  10000,
                                  100000,
                                  1000000,
                                  10000000,
                                  100000000,
                                  1000000000,
                                  10000000000,
                                  100000000000,
                                  1000000000000,
                                  10000000000000,
                                  100000000000000,
                                  1000000000000000};

static std::stack<std::chrono::steady_clock::time_point> tictoc_stack;
static std::stack<double> tictoc_stack_omp_clock;


/// start the timer. easy to call as tic(); some code; double t=toc();
inline void tic() { tictoc_stack.push(std::chrono::steady_clock::now()); }

/***
 * Returns the time taken between the most recent tic() to itself.
 * @return time in seconds.
*/
inline double toc() {
  std::chrono::duration<double> time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - tictoc_stack.top());
  double rc = time_span.count();
  tictoc_stack.pop();
  return rc;
}

template <class T>
void fixNumericalError(T *X, const double prec = EPSILON_1EMINUS16,
                       const double repl = 0.0) {
  (*X).for_each(
      [&](typename T::elem_type &val) { val = (val < prec) ? repl : val; });
}

template <class T>
void fixAbsNumericalError(T *X, const double prec = EPSILON_1EMINUS16,
                       const double repl = 0.0) {
  (*X).for_each([&](typename T::elem_type &val) {
        val = (std::abs(val) < prec) ? repl : val;
      });
}

template <class T>
void fixDecimalPlaces(T *X, const int places = NUMBEROF_DECIMAL_PLACES) {
  (*X).for_each([&](typename T::elem_type &val) {
    val = floorf(val * powersof10[places]) / powersof10[places];
  });
}

/*
 * Returns the nth prime number.
 * There are totally 10000 prime numbers within 104000;
 */
int random_sieve(const int nthprime) {
  int i, m, k;
  int klimit, nlimit;
  int *mark;

  nlimit = 104000;

  mark = reinterpret_cast<int *>(calloc(nlimit, sizeof(int)));

  /* Calculate limit for k */
  klimit = static_cast<int>(sqrt(static_cast<double>(nlimit) + 1));

  /* Mark the composites */
  /* Special case */
  mark[1] = -1;

  /* Set k=1. Loop until k >= sqrt(n) */
  for (k = 3; k <= klimit; k = m) {
    /* Find first non-composite in list > k */
    for (m = k + 1; m < nlimit; m++)
      if (!mark[m]) break;

    /* Mark the numbers 2m, 3m, 4m, ... */
    for (i = m * 2; i < nlimit; i += m) mark[i] = -1;
  }

  /* Now display results - all unmarked numbers are prime */
  int rcprime = -1;
  for (k = 0, i = 1; i < nlimit; i++) {
    if (!mark[i]) {
      k++;
      if (k == nthprime + 1) {
        rcprime = i;
        break;
      }
    }
  }
  free(mark);
  return rcprime;
}

template <class T>
void absmat(T *X) {
  arma::uvec negativeIdx = find((*X) < 0);
  (*X)(negativeIdx) = (*X)(negativeIdx) * -1;
}

template <class T>
void makeSparse(const double sparsity, T(*X)) {
  // make a matrix sparse
#ifndef USING_R
  srand(RAND_SEED_SPARSE);
#endif
#pragma omp parallel for default(none) shared(sparsity, X)
  for (arma::uword j = 0; j < X->n_cols; j++) {
    for (arma::uword i = 0; i < X->n_rows; i++) {
      if (arma::randu() > sparsity) (*X)(i, j) = 0;
    }
  }
}

void randNMF(const arma::uword m, const arma::uword n, const arma::uword k, const double sparsity,
             arma::mat *A) {
#ifndef USING_R
  srand(RAND_SEED);
#endif
  arma::mat W = 10 * arma::randu<arma::mat>(m, k);
  arma::mat H = 10 * arma::randu<arma::mat>(n, k);
  if (sparsity < 1) {
    makeSparse<arma::mat>(sparsity, &W);
    makeSparse<arma::mat>(sparsity, &H);
  }
  arma::mat temp = ceil(W * trans(H));
  A = &temp;
}
void randNMF(const arma::uword m, const arma::uword n, const arma::uword k, const double sparsity,
             arma::sp_mat *A) {
  auto temp = arma::sprandu<arma::sp_mat>(m, n, sparsity);
  A = &temp;
}

template <class T>
void printVector(const std::vector<T> &x) {
  for (int i = 0; i < x.size(); i++) {
    INFO << x[i] << ' ';
  }
  INFO << std::endl;
}

std::vector<std::vector<size_t>> cartesian_product(
    const std::vector<std::vector<size_t>> &v) {
  std::vector<std::vector<size_t>> s = {{}};
  for (auto &u : v) {
    std::vector<std::vector<size_t>> r;
    for (auto y : u) {
      for (auto &x : s) {
        r.push_back(x);
        r.back().push_back(y);
      }
    }
    s.swap(r);
  }
  return s;
}

/*
 * can be called by external people for sparse input matrix.
 */
template <class INPUTTYPE, class LRTYPE>
double computeObjectiveError(const INPUTTYPE &A, const LRTYPE &W,
                             const LRTYPE &H) {
  // 1. over all nnz (a_ij - w_i h_j)^2
  // 2. over all nnz (w_i h_j)^2
  // 3. Compute R of W ahd L of H through QR
  // 4. use sgemm to compute RL
  // 5. use slange to compute ||RL||_F^2
  // 6. return nnzsse+nnzwh-||RL||_F^2
  arma::uword k = W.n_cols;
  arma::uword m = A.n_rows;
  arma::uword n = A.n_cols;
  tic();
  double nnzsse = 0;
  double nnzwh = 0;
  LRTYPE Rw(k, k);
  LRTYPE Rh(k, k);
  LRTYPE Qw(m, k);
  LRTYPE Qh(n, k);
  LRTYPE RwRh(k, k);
#pragma omp parallel for reduction(+ : nnzsse, nnzwh) default(none)
  for (arma::uword jj = 1; jj <= A.n_cols; jj++) {
    arma::uword startIdx = A.col_ptrs[jj - 1];
    arma::uword endIdx = A.col_ptrs[jj];
    arma::uword col = jj - 1;
    double nnzssecol = 0;
    double nnzwhcol = 0;
    for (arma::uword ii = startIdx; ii < endIdx; ii++) {
      arma::uword row = A.row_indices[ii];
      double tempsum = 0;
      for (arma::uword kk = 0; kk < k; kk++) {
        tempsum += (W(row, kk) * H(col, kk));
      }
      nnzwhcol += tempsum * tempsum;
      nnzssecol += (A.values[ii] - tempsum) * (A.values[ii] - tempsum);
    }
    nnzsse += nnzssecol;
    nnzwh += nnzwhcol;
  }
  qr_econ(Qw, Rw, W);
  qr_econ(Qh, Rh, H);
  RwRh = Rw * Rh.t();
  double normWH = arma::norm(RwRh, "fro");
  Rw.clear();
  Rh.clear();
  Qw.clear();
  Qh.clear();
  RwRh.clear();
#ifdef _VERBOSE
  INFO << "error compute time " << toc() << std::endl;
#endif // _VERBOSE
  double fastErr = sqrt(nnzsse + (normWH * normWH - nnzwh));
  return (fastErr);
}

/*
 * This is an sgemm wrapper for armadillo matrices
 * Something is going crazy with armadillo
 */

void cblas_sgemm(const arma::mat &A, const arma::mat &B, double *C) {
  arma::uword m = A.n_rows;
  arma::uword n = B.n_cols;
  arma::uword k = A.n_cols;
  double alpha = 1.0;
  double beta = 0.0;
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (int32_t)m, (int32_t)n, (int32_t)k, alpha,
              A.memptr(), (int32_t)m, B.memptr(), (int32_t)k, beta, C, (int32_t)m);
}

/**
 * Generates the same low rank matrix. Matrix columns are seeded to
 * easily create different row slices of the matrix. The lowrank matrix
 * X is expected to be of size n x k where k << n.
 *
 * @param[in] row_start is the starting row of the block to be generated
 * @param[in] nrows is the number of rows to be generated
 * @param[in] k is the lowrank (length of each row vector)
 * @param[in] X is the reference to the matrix to populate
 * @param[in] trans is a flag to indicate that Xt is sent in.
 * @param[in] mseed is the seed for the first column of the matrix
 */
void gen_discard(arma::uword row_start, arma::uword nrows, arma::uword k,
        arma::mat &X, bool trans, int mseed=7907) {
  for(unsigned int j = 0; j < k; ++j) {
    std::mt19937 gen(mseed + j);
    gen.discard(row_start);
    for(unsigned int i = 0; i < nrows; ++i) {
      if (trans) {
          X(j, i) =  ((double)gen()) / std::mt19937::max();
      } else {
          X(i, j) =  ((double)gen()) / std::mt19937::max();
      }
    }
  }
}

/*
 * Read in a dense matrix
 */
void read_input_matrix(arma::mat &A, std::string fname) {
  A.load(std::move(fname));
}

/*
 * Read in a sparse matrix
 */
void read_input_matrix(arma::sp_mat &A, std::string fname) {
  A.load(std::move(fname), arma::coord_ascii);
}

/*
 * Generate random dense matrix
 */
void generate_rand_matrix(arma::mat &A, const std::string& rtype,
        arma::uword m, arma::uword n, arma::uword k, double density, bool symm_flag = false,
        bool adjrand = false, int kalpha = 1, int kbeta = 0) {
  if (rtype == "uniform") {
    if (symm_flag) {
      A = arma::randu<arma::mat>(m, n);
      A = 0.5 * (A + A.t());
    } else {
      A = arma::randu<arma::mat>(m, n);
    }
  } else if (rtype == "normal") {
    if (symm_flag) {
      A = arma::randn<arma::mat>(m, n);
      A = 0.5 * (A + A.t());
    } else {
      A = arma::randn<arma::mat>(m, n);
    }
    A.elem(find(A < 0)).zeros();
  } else {
    if (symm_flag) {
      arma::mat Htrue = arma::zeros<arma::mat>(n, k);
      gen_discard(0, n, k, Htrue, false, HTRUE_SEED);
      A = Htrue * Htrue.t();

      // Free auxiliary variables
      Htrue.clear();
    } else {
      arma::mat Wtrue = arma::zeros<arma::mat>(m, k);
      gen_discard(0, m, k, Wtrue, false, WTRUE_SEED);
      arma::mat Htrue = arma::zeros<arma::mat>(k, n);
      gen_discard(0, n, k, Htrue, true, HTRUE_SEED);
      A = Wtrue * Htrue;

      // Free auxiliary variables
      Wtrue.clear();
      Htrue.clear();
    }
  }
  if (adjrand) {
    A = kalpha * (A) + kbeta;
    A = ceil(A);
  }
}

/*
 * Generate random sparse matrix
 */
void generate_rand_matrix(arma::sp_mat &A, const std::string& rtype,
        arma::uword m, arma::uword n, arma::uword k, double density, bool symm_flag = false,
        bool adjrand = false, int kalpha = 5, int kbeta = 10) {
  if (rtype == "uniform") {
    if (symm_flag) {
      double dens = 0.5 * density;
      A = arma::sprandu<arma::sp_mat>(m, n, dens);
      A = 0.5 * (A + A.t());
    } else {
      A = arma::sprandu<arma::sp_mat>(m, n, density);
      INFO << size(nonzeros(A)) << std::endl;
    }
  } else if (rtype == "normal") {
    if (symm_flag) {
      double dens = 0.5 * density;
      A = arma::sprandn<arma::sp_mat>(m, n, dens);
      A = 0.5 * (A + A.t());
    } else {
      A = arma::sprandn<arma::sp_mat>(m, n, density);
    }
  } else if (rtype == "lowrank") {
    if (symm_flag) {
      double dens = 0.5 * density;
      auto mask = arma::sprandu<arma::sp_mat>(m, n, dens);
      mask = 0.5 * (mask + mask.t());
      mask = arma::spones(mask);
      arma::mat Htrue = arma::zeros(n, k);
      gen_discard(0, n, k, Htrue, false, HTRUE_SEED);
      A = arma::sp_mat(mask % (Htrue * Htrue.t()));

      // Free auxiliary space
      Htrue.clear();
      mask.clear();
    } else {
      auto mask = arma::sprandu<arma::sp_mat>(m, n, density);
      mask = arma::spones(mask);
      arma::mat Wtrue = arma::zeros(m, k);
      gen_discard(0, m, k, Wtrue, false, WTRUE_SEED);
      arma::mat Htrue = arma::zeros(k, n);
      gen_discard(0, n, k, Htrue, true, HTRUE_SEED);
      A = arma::sp_mat(mask % (Wtrue * Htrue));

      // Free auxiliary space
      Wtrue.clear();
      Htrue.clear();
      mask.clear();
    }
  }
  // Adjust and project non-zeros
  arma::sp_mat::iterator start_it = A.begin();
  arma::sp_mat::iterator end_it = A.end();
  for (arma::sp_mat::iterator it = start_it; it != end_it; ++it) {
    double curVal = (*it);
    if (adjrand) {
      (*it) = ceil(kalpha * curVal + kbeta);
    }
    if ((*it) < 0) (*it) = kbeta;
  }
}

int debug_hook(){
  int i = 0;
  while(i < 1){}
  return 0;
}

template<typename T>
arma::uword chunk_size_dense(arma::uword rank) {
#ifdef _OPENMP
return (get_l1_data_cache() / (rank * sizeof(T)));
#else
return (get_l2_data_cache() / (rank * sizeof(T)));
#endif
}

