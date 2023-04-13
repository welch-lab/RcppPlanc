// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "bppnmf.hpp"
#include "nmf.hpp"
#include "utils.hpp"
#include "bppnnls.hpp"
#include "nnls.hpp"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

normtype m_input_normalization;
int m_initseed;
arma::uword m_m, m_n;
int m_symm_flag;
double m_symm_reg;
arma::fvec m_regW;
arma::fvec m_regH;

//template <class NMFTYPE>
//void runNMF() {
//
//}


// Use the ANLS-BPP algorithm to factor a given matrix
// at the given rank. TODO FIX set.seed

// [[Rcpp::export]]
Rcpp::List rcppplanc_bppnmf(const arma::sp_mat & x, const int & k, const int & niter) {
  // Set parameters and call NMF
  arma::sp_mat A = x.t();
  // arma::arma_rng::set_seed(m_initseed);
  m_m = x.n_rows;
  m_n = x.n_cols;
  m_symm_flag = 0;
  m_symm_reg = -1;
  arma::mat W = arma::randu<arma::mat>(m_m, k);
  arma::mat H = arma::randu<arma::mat>(m_n, k);
  // symmreg here
  // double meanA = arma::mean(arma::mean(x));
  //  H = 2 * std::sqrt(meanA / k) * H;
  //  W = H;
  //  if (m_symm_reg == 0.0)
  //  {
  //    double symreg = x.max();
  //    m_symm_reg = symreg * symreg;
  //  }
  planc::BPPNMF<arma::sp_mat> MyNMF = planc::BPPNMF<arma::sp_mat>(x, W, H);
  MyNMF.num_iterations(niter);
  MyNMF.symm_reg(m_symm_reg);
  // MyNMF.compute_error(this->m_compute_error);
  // MyNMF.algorithm(this->m_nmfalgo);

  if (!m_regW.empty())
  {
    MyNMF.regW(m_regW);
  }
  if (!m_regH.empty())
  {
    MyNMF.regH(m_regH);
  }
  MyNMF.computeNMF();
  return Rcpp::List::create(
      Rcpp::Named("W") = MyNMF.getLeftLowRankFactor(),
      Rcpp::Named("H") = MyNMF.getRightLowRankFactor()
    );
}

// [[Rcpp::export]]
arma::mat rcppplanc_bppnnls(const arma::sp_mat &A, const arma::mat &B)
{
  arma::uword m_n = B.n_cols;
  arma::uword m_m = A.n_cols;
  int m_k = B.n_rows;
  arma::mat outmat = arma::randu<arma::mat>(m_m, m_k);
  arma::mat *outmatptr;
  UINT numChunks = m_n / ONE_THREAD_MATRIX_SIZE;
#pragma omp parallel for schedule(auto)
  for (UINT i = 0; i < numChunks; i++)
  {
    UINT spanStart = i * ONE_THREAD_MATRIX_SIZE;
    UINT spanEnd = (i + 1) * ONE_THREAD_MATRIX_SIZE - 1;
    if (spanEnd > m_n - 1)
    {
      spanEnd = m_n - 1;
    }
    // double start = omp_get_wtime();
    BPPNNLS<AMAT, VEC> solveProblem(B.t(), (AMAT)A.cols(spanStart, spanEnd));
    solveProblem.solveNNLS();
    // double end = omp_get_wtime();
    // titer = end - start;
    // #ifdef _VERBOSE
    // INFO << " start=" << spanStart
    //     << ", end=" << spanEnd
    //     << ", tid=" << omp_get_thread_num() << " cpu=" << sched_getcpu() << std::endl;
    //     // << " time taken=" << titer << std::endl;
    // #endif
    (*outmatptr).rows(spanStart, spanEnd) = solveProblem.getSolutionMatrix().t();
  };
  return outmat;
}
