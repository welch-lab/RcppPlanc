// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "bppnmf.hpp"
#include "nmf.hpp"
#include "utils.hpp"
#include "bppnnls.hpp"
#include "aoadmm.hpp"
#include "gnsym.hpp"
#include "hals.hpp"
#include "mu.hpp"
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

template <class T>
Rcpp::List RcallNMF(arma::sp_mat x, int k, int niter)
{
  arma::sp_mat A;
  m_m = x.n_rows;
  m_n = x.n_cols;
  m_symm_flag = 0;
  m_symm_reg = -1;
  arma::mat W = arma::randu<arma::mat>(m_m, k);
  arma::mat H = arma::randu<arma::mat>(m_n, k);
  T MyNMF(x, W, H);
  // symmreg here
  // double meanA = arma::mean(arma::mean(x));
  //  H = 2 * std::sqrt(meanA / k) * H;
  //  W = H;
  //  if (m_symm_reg == 0.0)
  //  {
  //    double symreg = x.max();
  //    m_symm_reg = symreg * symreg;
  //  }
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
      Rcpp::Named("H") = MyNMF.getRightLowRankFactor());
}

//' Alternating Direction Method of Multipliers NMF
//'
//' Use the AOADMM algorithm to factor a given matrix at the given rank.
//'
//' @param x Input matrix for factorization
//' @param k Factor matrix rank
//' @param niter Maximum number of nmf iterations
//' @export
//' @returns The calculated factor matrices as an Rcpp::list
//' @examplesIf require("Matrix")
//' aoadmmnmf(rsparsematrix(nrow = 100, ncol = 100, nnz = 10, symmetric = TRUE), 10, 10)
// [[Rcpp::export]]
Rcpp::List aoadmmnmf(const arma::sp_mat & x, const int & k, const int & niter) {
  return RcallNMF<planc::AOADMMNMF<arma::sp_mat>>(x, k, niter);
}

//' Gauss-Newton using Conjugate Gradients NMF
//'
//' Use the Gauss-Newton algorithm to factor a given matrix at the given rank.
//'
//' @param x Input matrix for factorization
//' @param k Factor matrix rank
//' @param niter Maximum number of nmf iterations
//' @export
//' @returns The calculated factor matrices as an Rcpp::list
//' @examplesIf require("Matrix")
//' gnsymnmf(rsparsematrix(nrow = 100, ncol = 100, nnz = 10, symmetric = TRUE), 10, 10)
// [[Rcpp::export]]
Rcpp::List gnsymnmf(const arma::sp_mat &x, const int &k, const int &niter)
{
  return RcallNMF<planc::GNSYMNMF<arma::sp_mat>>(x, k, niter);
}

//' Hierarchical Alternating Least Squares NMF
//'
//' Use the HALS algorithm to factor a given matrix at the given rank.
//'
//' @param x Input matrix for factorization
//' @param k Factor matrix rank
//' @param niter Maximum number of nmf iterations
//' @export
//' @returns The calculated factor matrices as an Rcpp::list
//' @examplesIf require("Matrix")
//' halsmnmf(rsparsematrix(nrow = 100, ncol = 100, nnz = 10, symmetric = TRUE), 10, 10)
// [[Rcpp::export]]
Rcpp::List halsnmf(const arma::sp_mat &x, const int &k, const int &niter)
{
  return RcallNMF<planc::HALSNMF<arma::sp_mat>>(x, k, niter);
}
//' Multiplicative Update NMF
//'
//' Use the MU algorithm to factor a given matrix
//' at the given rank.
//' @param x Input matrix for factorization
//' @param k Factor matrix rank
//' @param niter Maximum number of nmf iterations
//' @export
//' @returns The calculated factor matrices as an Rcpp::list
//' @examplesIf require("Matrix")
//' halsnmf(rsparsematrix(nrow = 100, ncol = 100, nnz = 10, symmetric = TRUE), 10, 10)
// [[Rcpp::export]]
Rcpp::List munmf(const arma::sp_mat &x, const int &k, const int &niter)
{
  return RcallNMF<planc::MUNMF<arma::sp_mat>>(x, k, niter);
}
//' Alternating  Nonnegative Least Squares with Block Principal Pivoting NMF
//'
//' Use the ANLS-BPP algorithm to factor a given matrix at the given rank.
//'
//' @param x Input matrix for factorization
//' @param k Factor matrix rank
//' @param niter Maximum number of nmf iterations
//' @export
//' @returns The calculated factor matrices as an Rcpp::list
//' @examplesIf require("Matrix")
//' munmf(rsparsematrix(nrow = 100, ncol = 100, nnz = 10, symmetric = TRUE), 10, 10)
// [[Rcpp::export]]
Rcpp::List bppnmf(const arma::sp_mat & x, const int & k, const int & niter) {
  return RcallNMF<planc::BPPNMF<arma::sp_mat>>(x, k, niter);
}
//' Block Principal Pivoted Non-Negative Least Squares
//'
//' Use the BPP algorithm to get the nonnegative least squares solution for the given matrices.
//'
//' @param A Input sparse matrix
//' @param B Input factor dense matrix
//' @export
//' @returns The calculated solution matrix in dense form.
//' @examplesIf require("Matrix")
//' bppnnls(rsparsematrix(nrow=20,ncol=20,nnz=10), Matrix(runif(n=200,min=0,max=2),20,10))
// [[Rcpp::export]]
arma::mat bppnnls(const arma::sp_mat &A, const arma::mat &B)
{
  arma::uword m_n = B.n_cols;
  arma::uword m_m = A.n_cols;
  int m_k = B.n_rows;
  arma::mat outmat = arma::randu<arma::mat>(m_m, m_k);
  arma::mat *outmatptr;
  unsigned int numChunks = m_n / ONE_THREAD_MATRIX_SIZE;
#pragma omp parallel for schedule(auto)
  for (unsigned int i = 0; i < numChunks; i++)
  {
    unsigned int spanStart = i * ONE_THREAD_MATRIX_SIZE;
    unsigned int spanEnd = (i + 1) * ONE_THREAD_MATRIX_SIZE - 1;
    if (spanEnd > m_n - 1)
    {
      spanEnd = m_n - 1;
    }
    // double start = omp_get_wtime();
    BPPNNLS<arma::mat, arma::vec> solveProblem(B.t(), (arma::mat)A.cols(spanStart, spanEnd));
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
