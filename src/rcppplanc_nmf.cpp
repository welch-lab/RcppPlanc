// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "config.h"
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <functional>
#include <utility>
#include "bppnmf.hpp"
#include "bppinmf.hpp"
#include "bppnnls.hpp"
#include "aoadmm.hpp"
#include "gnsym.hpp"
#include "hals.hpp"
#include "mu.hpp"
#include "data.hpp"
#include "onlineinmf.hpp"
#include "uinmf.hpp"
extern "C" {
#include "detect_blas.h"
}

// [[Rcpp::plugins(openmp)]]
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

normtype m_input_normalization;
int m_initseed;
arma::fvec m_regW;
arma::fvec m_regH;

// T2 e.g. arma::sp_mat
template <typename T2>
Rcpp::List runNMF(T2 x, arma::uword k, const std::string& algo, const arma::uword& niter, const int& nCores,
                  Rcpp::Nullable<Rcpp::NumericMatrix> Winit,
                  Rcpp::Nullable<Rcpp::NumericMatrix> Hinit) {
    libcall = planc::nmf(x, k, niter, algo, Winit, Hinit);
    return Rcpp::List::create(
        Rcpp::Named("W") = MyNMF.getLeftLowRankFactor(),
        Rcpp::Named("H") = MyNMF.getRightLowRankFactor(),
        Rcpp::Named("objErr") = MyNMF.objErr()
    );
}

//' @title Perform Non-negative Matrix Factorization
//' @description
//' Regularly, Non-negative Matrix Factorization (NMF) is factorizes input
//' matrix \eqn{X} into low rank matrices \eqn{W} and \eqn{H}, so that
//' \eqn{X \approx WH}. The objective function can be stated as
//' \eqn{\arg\min_{W\ge0,H\ge0}||X-WH||_F^2}. In practice, \eqn{X} is usually
//' regarded as a matrix of \eqn{m} features by \eqn{n} sample points. And the
//' result matrix \eqn{W} should have the dimensionality of \eqn{m \times k} and
//' H with \eqn{n \times k} (transposed).
//' This function wraps the algorithms implemented in PLANC library to solve
//' NMF problems. Algorithms includes Alternating Non-negative Least Squares
//' with Block Principal Pivoting (ANLS-BPP), Alternating Direction Method of
//' Multipliers (ADMM), Hierarchical Alternating Least Squares (HALS), and
//' Multiplicative Update (MU).
//' @param x Input matrix for factorization. Can be either dense or sparse.
//' @param k Integer. Factor matrix rank.
//' @param niter Integer. Maximum number of NMF interations.
//' @param algo Algorithm to perform the factorization, choose from "anlsbpp",
//' "admm", "hals" or "mu". See detailed sections.
//' @param nCores The number of parallel tasks that will be spawned. Only applies to anlsbpp.
//' Default \code{2}
//' @param Winit Initial left-hand factor matrix, must be of size m x k.
//' @param Hinit Initial right-hand factor matrix, must be of size n x k.
//' @returns A list with the following elements:
//' \itemize{
//'  \item{\code{W} - the result left-hand factor matrix}
//'  \item{\code{H} - the result right hand matrix.}
//'  \item{\code{objErr} - the objective error of the factorization.}
//' }
//' @references
//' Ramakrishnan Kannan and et al., A High-Performance Parallel Algorithm for
//' Nonnegative Matrix Factorization, PPoPP '16, 2016, 10.1145/2851141.2851152
// [[Rcpp::export]]
Rcpp::List nmf(const SEXP& x, const arma::uword &k, const arma::uword &niter = 30,
               const std::string &algo = "anlsbpp",
               const int& nCores = 2,
               const Rcpp::Nullable<Rcpp::NumericMatrix> &Winit = R_NilValue,
               const Rcpp::Nullable<Rcpp::NumericMatrix> &Hinit = R_NilValue) {
    Rcpp::List outlist;
    if (Rf_isS4(x)) {
        // Assume using dgCMatrix
        outlist = runNMF<arma::sp_mat>(Rcpp::as<arma::sp_mat>(x), k, algo, niter, nCores, Winit, Hinit);
    } else {
        // Assume regular dense matrix
        outlist = runNMF<arma::mat>(Rcpp::as<arma::mat>(x), k, algo, niter, nCores, Winit, Hinit);
    }
    return outlist;
}

// T1 e.g. BPPNMF<arma::sp_mat>
// T2 e.g. arma::sp_mat
template <class T1, class T2>
Rcpp::List runSymNMF(T2 x, arma::uword k, const int& nCores, arma::uword niter, double symm_reg,
                     Rcpp::Nullable<Rcpp::NumericMatrix> Hinit) {
  arma::uword m = x.n_rows;
  arma::uword n = x.n_cols;
  if (m != n) {
    Rcpp::stop("Input `x` is not square.");
  }
  if (k >= m) {
    Rcpp::stop("`k` must be less than `nrow(x)`");
  }
  arma::mat W(n, k); // W and H are the same for symNMF!
  arma::mat H(n, k);

  if (Hinit.isNotNull()) {
    H = Rcpp::as<arma::mat>(Hinit);
    if (H.n_rows != n || H.n_cols != k) {
      Rcpp::stop("Hinit must be of size " +
        std::to_string(n) + " x " + std::to_string(k));
    }
  } else {
    H = arma::randu<arma::mat>(n, k);
    double meanX = arma::mean(arma::mean(x));
    H = 2 * std::sqrt(meanX / k) * H;
    W = H;
    if (symm_reg == 0.0) {
      double symreg = x.max();
      symm_reg = symreg * symreg;
    }
  }
  T1 MyNMF(x, W, H, nCores);
  MyNMF.num_iterations(niter);
  MyNMF.symm_reg(symm_reg);
  if (!m_regW.empty()) MyNMF.regW(m_regW);
  if (!m_regH.empty()) MyNMF.regH(m_regH);
  MyNMF.computeNMF();
  return Rcpp::List::create(
    Rcpp::Named("W") = MyNMF.getLeftLowRankFactor(),
    Rcpp::Named("H") = MyNMF.getRightLowRankFactor(),
    Rcpp::Named("objErr") = MyNMF.objErr()
  );
}

//' Perform Symmetric Non-negative Matrix Factorization
//'
//' Symmetric input matrix \eqn{X} of size \eqn{n \times n} is required. Two
//' approaches are provided. Alternating Non-negative Least Squares Block
//' Principal Pivoting algorithm (ANLSBPP) with symmetric regularization, where
//' the objective function is set to be \eqn{\arg\min_{H\ge0,W\ge0}||X-WH||_F^2+
//' \lambda||W-H||_F^2}, can be run with \code{algo = "anlsbpp"}.
//' Gaussian-Newton algorithm, where the objective function is set to be
//' \eqn{\arg\min_{H\ge0}||X-H^\mathsf{T}H||_F^2}, can be run with \code{algo =
//' "gnsym"}. In the objectives, \eqn{W} is of size \eqn{n \times k} and \eqn{H}
//' is of size \eqn{k \times n}. The returned results will all be
//' \eqn{n \times k}.
//'
//' @param x Input matrix for factorization. Must be symmetric. Can be either
//' dense or sparse.
//' @param k Integer. Factor matrix rank.
//' @param niter Integer. Maximum number of symNMF interations.
//' Default \code{30}
//' @param lambda Symmetric regularization parameter. Must be
//' non-negative. Default \code{0.0} uses the square of the maximum value in
//' \code{x}.
//' @param algo Algorithm to perform the factorization, choose from "gnsym" or
//' "anlsbpp". Default \code{"gnsym"}
//' @param nCores The number of parallel tasks that will be spawned. Only applies to anlsbpp.
//' Default \code{2}
//' @param Hinit Initial right-hand factor matrix, must be of size n x k.
//' Default \code{NULL}.
//' @returns A list with the following elements:
//' \itemize{
//'  \item{\code{W} - the result left-hand factor matrix, non-empty when using
//'  \code{"anlsbpp"}}
//'  \item{\code{H} - the result right hand matrix.}
//'  \item{\code{objErr} - the objective error of the factorization.}
//' }
//' @references
//' Srinivas Eswar and et al., Distributed-Memory Parallel Symmetric Nonnegative
//' Matrix Factorization, SC '20, 2020, 10.5555/3433701.3433799
// [[Rcpp::export()]]
Rcpp::List symNMF(const SEXP& x, const arma::uword& k, const arma::uword& niter = 30,
                 const double& lambda = 0.0, const std::string& algo = "gnsym", const int& nCores = 2,
                 const Rcpp::Nullable<Rcpp::NumericMatrix> &Hinit = R_NilValue) {
//   arma::mat out;
  Rcpp::List out;
  if (Rf_isS4(x)) {
    // Assume using dgCMatrix
    if (algo == "gnsym") {
      out = runSymNMF<planc::GNSYMNMF<arma::sp_mat>, arma::sp_mat>(
        Rcpp::as<arma::sp_mat>(x), k, nCores, niter, lambda, Hinit
      );
    } else if (algo == "anlsbpp") {
      out = runSymNMF<planc::BPPNMF<arma::sp_mat>, arma::sp_mat>(
        Rcpp::as<arma::sp_mat>(x), k, nCores, niter, lambda, Hinit
      );
    } else {
      Rcpp::stop(R"(Please choose `algo` from "gnsym" or "anlsbpp".)");
    }
  } else {
    // Assume using default matrix
    if (algo == "gnsym") {
      out = runSymNMF<planc::GNSYMNMF<arma::mat>, arma::mat>(
        Rcpp::as<arma::mat>(x), k, nCores, niter, lambda, Hinit
      );
    } else if (algo == "anlsbpp") {
      out = runSymNMF<planc::BPPNMF<arma::mat>, arma::mat>(
        Rcpp::as<arma::mat>(x), k, nCores, niter, lambda, Hinit
      );
    } else {
      Rcpp::stop(R"(Please choose `algo` from "gnsym" or "anlsbpp".)");
    }
  }
  return out;
}


// template <class T>
// Rcpp::List RcallNMF(arma::sp_mat x, int k, int niter)
// {
//   arma::uword m_m = x.n_rows;
//   arma::uword m_n = x.n_cols;
//   arma::mat W = arma::randu<arma::mat>(m_m, k);
//   arma::mat H = arma::randu<arma::mat>(m_n, k);
//   T MyNMF(x, W, H);
//   // symmreg here
//   // double meanA = arma::mean(arma::mean(x));
//   //  H = 2 * std::sqrt(meanA / k) * H;
//   //  W = H;
//   //  if (m_symm_reg == 0.0)
//   //  {
//   //    double symreg = x.max();
//   //    m_symm_reg = symreg * symreg;
//   //  }
//   MyNMF.num_iterations(niter);
//   MyNMF.symm_reg(-1);
//   // MyNMF.compute_error(this->m_compute_error);

//   if (!m_regW.empty())
//   {
//     MyNMF.regW(m_regW);
//   }
//   if (!m_regH.empty())
//   {
//     MyNMF.regH(m_regH);
//   }
//   MyNMF.computeNMF();
//   return Rcpp::List::create(
//       Rcpp::Named("W") = MyNMF.getLeftLowRankFactor(),
//       Rcpp::Named("H") = MyNMF.getRightLowRankFactor());
// }

// template <class T>
// Rcpp::List RcallNMF(arma::sp_mat x, int k, int niter,
//                     arma::mat Winit, arma::mat Hinit) {
//   T MyNMF(x, Winit, Hinit);
//   MyNMF.num_iterations(niter);
//   MyNMF.symm_reg(-1);
//   MyNMF.computeNMF();
//   return Rcpp::List::create(
//       Rcpp::Named("W") = MyNMF.getLeftLowRankFactor(),
//       Rcpp::Named("H") = MyNMF.getRightLowRankFactor());
// }





// //' Alternating Direction Method of Multipliers NMF
// //'
// //' Use the AOADMM algorithm to factor a given matrix at the given rank.
// //'
// //' @param x Input matrix for factorization
// //' @param k Factor matrix rank
// //' @param niter Maximum number of nmf iterations
// //' @param H_init Initial right-hand factor matrix (Optional)
// //' @param W_init Initial left-hand factor matrix (Optional)
// //' @returns The calculated factor matrices as an Rcpp::List
// //' @examplesIf require("Matrix")
// //' aoadmmnmf(rsparsematrix(nrow = 25, ncol = 25, nnz = 5, symmetric = TRUE), 10, 10)
// // [[Rcpp::export]]
// Rcpp::List aoadmmnmf(const arma::sp_mat &x, const int &k, const int &niter,
//                      const Rcpp::Nullable<Rcpp::NumericMatrix> &W_init = R_NilValue,
//                      const Rcpp::Nullable<Rcpp::NumericMatrix> &H_init = R_NilValue) {
//   Rcpp::List outlist;
//   if (W_init.isNotNull() && W_init.isNotNull()) {
//     outlist = RcallNMF<planc::AOADMMNMF<arma::sp_mat>>(x, k, niter,
//                                                        Rcpp::as<arma::mat>(W_init),
//                                                        Rcpp::as<arma::mat>(H_init));
//   }
//   else {
//     outlist = RcallNMF<planc::AOADMMNMF<arma::sp_mat>>(x, k, niter);
//   }
//   return outlist;
// }

// //' Gauss-Newton using Conjugate Gradients NMF
// //'
// //' Use the Gauss-Newton algorithm to factor a given matrix at the given rank.
// //'
// //' @param x Input matrix for factorization
// //' @param k Factor matrix rank
// //' @param niter Maximum number of nmf iterations
// //' @param H_init Initial right-hand factor matrix (Optional)
// //' @param W_init Initial left-hand factor matrix (Optional)
// //' @returns The calculated factor matrices as an Rcpp::List
// //' @examplesIf require("Matrix")
// //' gnsymnmf(rsparsematrix(nrow = 25, ncol = 25, nnz = 5, symmetric = TRUE), 10, 10)
// // [[Rcpp::export]]
// Rcpp::List gnsymnmf(const arma::sp_mat &x, const int &k, const int &niter,
//                     const Rcpp::Nullable<Rcpp::NumericMatrix> &W_init = R_NilValue,
//                     const Rcpp::Nullable<Rcpp::NumericMatrix> &H_init = R_NilValue)
// {
//   Rcpp::List outlist;
//   if (W_init.isNotNull() && W_init.isNotNull()) {
//     outlist = RcallNMF<planc::GNSYMNMF<arma::sp_mat>>(x, k, niter,
//                                                       Rcpp::as<arma::mat>(W_init),
//                                                       Rcpp::as<arma::mat>(H_init));
//   }
//   else {
//     outlist = RcallNMF<planc::GNSYMNMF<arma::sp_mat>>(x, k, niter);
//   }
//   return outlist;
// }

// //' Hierarchical Alternating Least Squares NMF
// //'
// //' Use the HALS algorithm to factor a given matrix at the given rank.
// //'
// //' @param x Input matrix for factorization
// //' @param k Factor matrix rank
// //' @param niter Maximum number of nmf iterations
// //' @param H_init Initial right-hand factor matrix (Optional)
// //' @param W_init Initial left-hand factor matrix (Optional)
// //' @returns The calculated factor matrices as an Rcpp::List
// //' @examplesIf require("Matrix")
// //' halsnmf(rsparsematrix(nrow = 25, ncol = 25, nnz = 5, symmetric = TRUE), 10, 10)
// // [[Rcpp::export]]
// Rcpp::List halsnmf(const arma::sp_mat &x, const int &k, const int &niter,
//                    const Rcpp::Nullable<Rcpp::NumericMatrix> &W_init = R_NilValue,
//                    const Rcpp::Nullable<Rcpp::NumericMatrix> &H_init = R_NilValue)
// {
//   Rcpp::List outlist;
//   if (W_init.isNotNull() && W_init.isNotNull()) {
//     outlist = RcallNMF<planc::HALSNMF<arma::sp_mat>>(x, k, niter,
//                                                      Rcpp::as<arma::mat>(W_init),
//                                                      Rcpp::as<arma::mat>(H_init));
//   }
//   else {
//     outlist = RcallNMF<planc::HALSNMF<arma::sp_mat>>(x, k, niter);
//   }
//   return outlist;
// }
// //' Multiplicative Update NMF
// //'
// //' Use the MU algorithm to factor a given matrix
// //' at the given rank.
// //' @param x Input matrix for factorization
// //' @param k Factor matrix rank
// //' @param niter Maximum number of nmf iterations
// //' @param H_init Initial right-hand factor matrix (Optional)
// //' @param W_init Initial left-hand factor matrix (Optional)
// //' @returns The calculated factor matrices as an Rcpp::List
// //' @examplesIf require("Matrix")
// //' halsnmf(rsparsematrix(nrow = 25, ncol = 25, nnz = 5, symmetric = TRUE), 10, 10)
// // [[Rcpp::export]]
// Rcpp::List munmf(const arma::sp_mat &x, const int &k, const int &niter,
//                  const Rcpp::Nullable<Rcpp::NumericMatrix> &W_init = R_NilValue,
//                  const Rcpp::Nullable<Rcpp::NumericMatrix> &H_init = R_NilValue)
// {
//   Rcpp::List outlist;
//   if (W_init.isNotNull() && W_init.isNotNull()) {
//     outlist = RcallNMF<planc::MUNMF<arma::sp_mat>>(x, k, niter,
//                                                    Rcpp::as<arma::mat>(W_init),
//                                                    Rcpp::as<arma::mat>(H_init));
//   }
//   else {
//     outlist = RcallNMF<planc::MUNMF<arma::sp_mat>>(x, k, niter);
//   }
//   return outlist;
// }
// //' Alternating  Nonnegative Least Squares with Block Principal Pivoting NMF
// //'
// //' Use the ANLS-BPP algorithm to factor a given matrix at the given rank.
// //'
// //' @param x Input matrix for factorization
// //' @param k Factor matrix rank
// //' @param niter Maximum number of nmf iterations
// //' @param H_init Initial right-hand factor matrix (Optional)
// //' @param W_init Initial left-hand factor matrix (Optional)
// //' @returns The calculated factor matrices as an Rcpp::List
// //' @examplesIf require("Matrix")
// //' munmf(rsparsematrix(nrow = 25, ncol = 25, nnz = 5, symmetric = TRUE), 10, 10)
// // [[Rcpp::export]]
// Rcpp::List bppnmf(const arma::sp_mat & x, const int & k, const int & niter,
//                   const Rcpp::Nullable<Rcpp::NumericMatrix> &W_init = R_NilValue,
//                   const Rcpp::Nullable<Rcpp::NumericMatrix> &H_init = R_NilValue) {
//   Rcpp::List outlist;
//   if (W_init.isNotNull() && W_init.isNotNull()) {
//     outlist = RcallNMF<planc::BPPNMF<arma::sp_mat>>(x, k, niter,
//                                                     Rcpp::as<arma::mat>(W_init),
//                                                     Rcpp::as<arma::mat>(H_init));
//   }
//   else {
//     outlist = RcallNMF<planc::BPPNMF<arma::sp_mat>>(x, k, niter);
//   }
//   return outlist;
// }

template <typename T>
arma::mat runbppnnls(const arma::mat &C, const T &B, int ncores) {
    arma::uword m_n = B.n_cols;
    arma::uword m_k = C.n_cols;
    arma::mat CtC = C.t() * C;
    arma::mat outmat = arma::zeros<arma::mat>(m_k, m_n);
    arma::mat *outmatptr;
    outmatptr = &outmat;
    arma::uword ONE_THREAD_MATRIX_SIZE = chunk_size_dense<double>(m_k);
    unsigned int numChunks = m_n / ONE_THREAD_MATRIX_SIZE;
    if (numChunks*ONE_THREAD_MATRIX_SIZE < m_n) numChunks++;
#pragma omp parallel for schedule(dynamic) default(none) shared(numChunks, ONE_THREAD_MATRIX_SIZE, m_n, outmatptr, C, B, CtC) num_threads(ncores)
    for (unsigned int i = 0; i < numChunks; i++) {
        unsigned int spanStart = i * ONE_THREAD_MATRIX_SIZE;
        unsigned int spanEnd = (i + 1) * ONE_THREAD_MATRIX_SIZE - 1;
        if (spanEnd > m_n - 1) spanEnd = m_n - 1;
        // double start = omp_get_wtime();
        arma::mat CtBChunk = C.t() * B.cols(spanStart, spanEnd);
        BPPNNLS<arma::mat, arma::vec> solveProblem(CtC, CtBChunk, true);
        solveProblem.solveNNLS();
        (*outmatptr).cols(spanStart, spanEnd) = solveProblem.getSolutionMatrix();
    }
    return outmat;
}

//' Block Principal Pivoted Non-Negative Least Squares
//'
//' Use the BPP algorithm to get the nonnegative least squares solution. Regular
//' NNLS problem is described as optimizing \eqn{\min_{x\ge0}||CX - B||_F^2}
//' where \eqn{C} and \eqn{B} are given and \eqn{X} is to be solved.
//' \code{bppnnls} takes \eqn{C} and \eqn{B} as input. \code{bppnnls_prod} takes
//' \eqn{C^\mathsf{T}C} and \eqn{C^\mathsf{T}B} as
//' input to directly go for the intermediate step of BPP algorithm. This can be
//' useful when the dimensionality of \eqn{C} and \eqn{B} is large while
//' pre-calculating \eqn{C^\mathsf{T}C} and \eqn{C^\mathsf{T}B} is cheap.
//'
//' @param C Input dense \eqn{C} matrix
//' @param B Input \eqn{B} matrix of either dense or sparse form
//' @param nCores The number of parallel tasks that will be spawned.
//' Default \code{2}
//' @returns The calculated solution matrix in dense form.
//' @rdname bppnnls
//' @examples
//' set.seed(1)
//' C <- matrix(rnorm(250), nrow = 25)
//' B <- matrix(rnorm(375), nrow = 25)
//' res1 <- bppnnls(C, B)
//' dim(res1)
//' res2 <- bppnnls_prod(t(C) %*% C, t(C) %*% B)
//' all.equal(res1, res2)
// [[Rcpp::export]]
arma::mat bppnnls(const arma::mat &C, const SEXP &B, const int& nCores = 2) {
    if (Rf_isS4(B)) {
        return runbppnnls<arma::sp_mat>(C, Rcpp::as<arma::sp_mat>(B), nCores);
    } else {
        return runbppnnls<arma::mat>(C, Rcpp::as<arma::mat>(B), nCores);
    }
    return {};
}

//' @param CtC The \eqn{C^\mathsf{T}C} matrix, see description.
//' @param CtB The \eqn{C^\mathsf{T}B} matrix, see description.
//' @param nCores The number of parallel tasks that will be spawned.
//' Default \code{2}
//' @rdname bppnnls
// [[Rcpp::export]]
arma::mat bppnnls_prod(const arma::mat &CtC, const arma::mat &CtB, const int& nCores = 2) {
    arma::uword n = CtB.n_cols;
    arma::uword k = CtC.n_cols;
    arma::uword ONE_THREAD_MATRIX_SIZE = chunk_size_dense<double>(k);
    arma::mat outmat = arma::zeros<arma::mat>(k, n);
    arma::mat *outmatptr;
    outmatptr = &outmat;
    unsigned int numChunks = n / ONE_THREAD_MATRIX_SIZE;
    if (numChunks*ONE_THREAD_MATRIX_SIZE < n) numChunks++;
#pragma omp parallel for schedule(dynamic) default(none) shared(numChunks, CtB, ONE_THREAD_MATRIX_SIZE, outmatptr, CtC, n) num_threads(nCores)
    for (unsigned int i = 0; i < numChunks; i++) {
        unsigned int spanStart = i * ONE_THREAD_MATRIX_SIZE;
        unsigned int spanEnd = (i + 1) * ONE_THREAD_MATRIX_SIZE - 1;
        if (spanEnd > n - 1) spanEnd = n - 1;
        arma::mat CtBChunk = CtB.cols(spanStart, spanEnd);
        BPPNNLS<arma::mat, arma::vec> solveProblem(CtC, CtBChunk, true);
        solveProblem.solveNNLS();
        (*outmatptr).cols(spanStart, spanEnd) = solveProblem.getSolutionMatrix();
    }
    return outmat;
}


// %%%%%%%%%%%%%%%%%% BPPINMF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template <typename T>
std::vector<std::unique_ptr<T>> initMemMatPtr(std::vector<T> objectList)
{
    std::vector<std::unique_ptr<T>> matPtrVec;
    for (arma::uword i = 0; i < objectList.size(); ++i)
    {
        T E = T(objectList[i]);
        std::unique_ptr<T> ptr = std::make_unique<T>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    return matPtrVec;
}

template <typename T>
Rcpp::List runINMF(std::vector<T> objectList, arma::uword k, double lambda,
                   arma::uword niter, bool verbose, const int& ncores)
{
    std::vector<std::unique_ptr<T>> matPtrVec;
    matPtrVec = initMemMatPtr<T>(objectList);
    planc::BPPINMF<T> solver(matPtrVec, k, lambda);
    solver.initW();
    solver.initV();
    solver.initH();
    solver.optimizeALS(niter, verbose, ncores);
    std::vector<Rcpp::NumericMatrix> HList;
    std::vector<Rcpp::NumericMatrix> VList;
    for (arma::uword i = 0; i < objectList.size(); ++i)
    {
        HList.push_back(Rcpp::wrap(solver.getHi(i)));
        VList.push_back(Rcpp::wrap(solver.getVi(i)));
    }
    return Rcpp::List::create(
        Rcpp::Named("H") = Rcpp::wrap(HList),
        Rcpp::Named("V") = Rcpp::wrap(VList),
        Rcpp::Named("W") = solver.getW(),
        Rcpp::Named("objErr") = solver.objErr());
}
template <typename T>
Rcpp::List runINMF(std::vector<T> objectList, arma::uword k, double lambda,
                   arma::uword niter, bool verbose,
                   std::vector<arma::mat> HinitList, std::vector<arma::mat> VinitList, arma::mat Winit,
                   const int& ncores)
{
    std::vector<std::unique_ptr<T>> matPtrVec;
    matPtrVec = initMemMatPtr<T>(objectList);
    planc::BPPINMF<T> solver(matPtrVec, k, lambda);
    solver.initW(Winit);
    solver.initV(VinitList);
    solver.initH(HinitList);
    solver.optimizeALS(niter, verbose, ncores);
    std::vector<Rcpp::NumericMatrix> HList;
    std::vector<Rcpp::NumericMatrix> VList;
    for (arma::uword i = 0; i < objectList.size(); ++i)
    {
        HList.push_back(Rcpp::wrap(solver.getHi(i)));
        VList.push_back(Rcpp::wrap(solver.getVi(i)));
    }
    return Rcpp::List::create(
        Rcpp::Named("H") = Rcpp::wrap(HList),
        Rcpp::Named("V") = Rcpp::wrap(VList),
        Rcpp::Named("W") = solver.getW(),
        Rcpp::Named("objErr") = solver.objErr());
}

Rcpp::List bppinmf_dense(const std::vector<arma::mat>& objectList, arma::uword k,
                         double lambda, arma::uword niter, const int& nCores, bool verbose = true,
                         Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
                         Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
                         Rcpp::Nullable<arma::mat> Winit = R_NilValue)
{
    if (Hinit.isNotNull() && Vinit.isNotNull() && Winit.isNotNull())
    {
        return runINMF<arma::mat>(objectList, k, lambda,
                                     niter, verbose,
                                     Rcpp::as<std::vector<arma::mat>>(Hinit),
                                     Rcpp::as<std::vector<arma::mat>>(Vinit),
                                     Rcpp::as<arma::mat>(Winit), nCores);
    }
    else
    {
        return runINMF<arma::mat>(objectList, k, lambda,
                                     niter, verbose, nCores);
    }
}

Rcpp::List bppinmf_sparse(const std::vector<arma::sp_mat> &objectList, arma::uword k, double lambda,
                          arma::uword niter, const int &nCores, bool verbose = true,
                          Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
                          Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
                          Rcpp::Nullable<arma::mat> Winit = R_NilValue)
{
    if (Hinit.isNotNull() && Vinit.isNotNull() && Winit.isNotNull()) {
        return runINMF<arma::sp_mat>(objectList, k, lambda,
                                     niter, verbose,
                                     Rcpp::as<std::vector<arma::mat>>(Hinit),
                                     Rcpp::as<std::vector<arma::mat>>(Vinit),
                                     Rcpp::as<arma::mat>(Winit), nCores);
    }
    else {
        return runINMF<arma::sp_mat>(objectList, k, lambda,
                niter, verbose, nCores);
    }
}

// [[Rcpp::export(.bppinmf)]]
Rcpp::List bppinmf(Rcpp::List objectList, const arma::uword k, const int& nCores,
                   const double lambda = 5, const arma::uword niter = 30,
                   const bool verbose = true,
                   Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
                   Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
                   Rcpp::Nullable<arma::mat> Winit = R_NilValue) {
    if (Rf_isS4(objectList[0])) {
        return bppinmf_sparse(Rcpp::as<std::vector<arma::sp_mat>>(objectList), k, lambda,
                        niter, nCores, verbose, Hinit, Vinit, Winit);
    } else {
        return bppinmf_dense(Rcpp::as<std::vector<arma::mat>>(objectList), k, lambda,
                            niter, nCores, verbose, Hinit, Vinit, Winit);
    }
    return Rcpp::List::create();
}

// [[Rcpp::export(.bppinmf_h5dense)]]
Rcpp::List bppinmf_h5dense(std::vector<std::string> filenames, std::vector<std::string> dataPath,
    arma::uword k, const int& nCores, double lambda, arma::uword niter, bool verbose = true,
    Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
    Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
    Rcpp::Nullable<arma::mat> Winit  = R_NilValue) {
    std::vector<std::unique_ptr<planc::H5Mat>> matPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        planc::H5Mat h5m(filenames[i], dataPath[i]);
        std::unique_ptr<planc::H5Mat> ptr = std::make_unique<planc::H5Mat>(h5m);
        matPtrVec.push_back(std::move(ptr));
    }
    planc::BPPINMF<planc::H5Mat> solver(matPtrVec, k, lambda);

    if (Winit.isNotNull()) {
        auto W = Rcpp::as<arma::mat>(Winit);
        solver.initW(W);
    } else {
        solver.initW();
    }

    if (Vinit.isNotNull()) {
        auto VinitList = Rcpp::as<std::vector<arma::mat>>(Vinit);
        solver.initV(VinitList);
    } else {
        solver.initV();
    }

    if (Hinit.isNotNull()) {
        auto HinitList = Rcpp::as<std::vector<arma::mat>>(Hinit);
        solver.initH(HinitList);
    } else {
        solver.initH();
    }

    solver.optimizeALS(niter, verbose, nCores);
    Rcpp::List HList(filenames.size());
    Rcpp::List VList(filenames.size());
    for (arma::uword i = 0; i < filenames.size(); ++i)
    {
        HList[i] = solver.getHi(i);
        VList[i] = solver.getVi(i);
    }
    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("H") = Rcpp::wrap(HList),
        Rcpp::Named("V") = Rcpp::wrap(VList),
        Rcpp::Named("W") = Rcpp::wrap(solver.getW()),
        Rcpp::Named("objErr") = Rcpp::wrap(solver.objErr())
    );
    return output;
    }

// [[Rcpp::export(.bppinmf_h5sparse)]]
Rcpp::List bppinmf_h5sparse(
    std::vector<std::string> filenames,
    std::vector<std::string> valuePath,
    std::vector<std::string> rowindPath,
    std::vector<std::string> colptrPath,
    arma::uvec nrow, arma::uvec ncol,
    arma::uword k, const int& nCores, double lambda, arma::uword niter,
    bool verbose = true,
    Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
    Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
    Rcpp::Nullable<arma::mat> Winit  = R_NilValue) {
    std::vector<std::unique_ptr<planc::H5SpMat>> matPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        planc::H5SpMat h5spm(filenames[i], rowindPath[i], colptrPath[i], valuePath[i], nrow[i], ncol[i]);
        std::unique_ptr<planc::H5SpMat> ptr = std::make_unique<planc::H5SpMat>(h5spm);
        matPtrVec.push_back(std::move(ptr));
    }
    planc::BPPINMF<planc::H5SpMat> solver(matPtrVec, k, lambda);

    if (Winit.isNotNull()) {
        auto W = Rcpp::as<arma::mat>(Winit);
        solver.initW(W);
    } else {
        solver.initW();
    }

    if (Vinit.isNotNull()) {
        auto VinitList = Rcpp::as<std::vector<arma::mat>>(Vinit);
        solver.initV(VinitList);
    } else {
        solver.initV();
    }

    if (Hinit.isNotNull()) {
        auto HinitList = Rcpp::as<std::vector<arma::mat>>(Hinit);
        solver.initH(HinitList);
    } else {
        solver.initH();
    }

    solver.optimizeALS(niter, verbose, nCores);
    Rcpp::List HList(filenames.size());
    Rcpp::List VList(filenames.size());
    for (arma::uword i = 0; i < filenames.size(); ++i)
    {
        HList[i] = solver.getHi(i);
        VList[i] = solver.getVi(i);
    }
    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("H") = Rcpp::wrap(HList),
        Rcpp::Named("V") = Rcpp::wrap(VList),
        Rcpp::Named("W") = Rcpp::wrap(solver.getW()),
        Rcpp::Named("objErr") = Rcpp::wrap(solver.objErr())
    );
    return output;
}


// %%%%%%%%%%%%%%%%%% online INMF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template <typename T>
Rcpp::List onlineINMF_S1_mem(std::vector<T> objectList, arma::uword k, const int& nCores,
    double lambda, arma::uword maxEpoch = 5, arma::uword minibatchSize = 5000,
    arma::uword maxHALSIter = 1, bool verbose = true) {
    std::vector<std::unique_ptr<T>> matPtrVec = initMemMatPtr<T>(objectList);
    planc::ONLINEINMF<T, T> solver(matPtrVec, k, lambda);
    solver.runOnlineINMF(minibatchSize, maxEpoch, maxHALSIter, verbose, nCores);

    Rcpp::List HList(objectList.size());
    Rcpp::List VList(objectList.size());
    Rcpp::List AList(objectList.size());
    Rcpp::List BList(objectList.size());
    for (arma::uword i = 0; i < objectList.size(); ++i)
    {
        HList[i] = solver.getHi(i);
        VList[i] = solver.getVi(i);
        AList[i] = solver.getAi(i);
        BList[i] = solver.getBi(i);
    }
    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("H") = Rcpp::wrap(HList),
        Rcpp::Named("V") = Rcpp::wrap(VList),
        Rcpp::Named("W") = Rcpp::wrap(solver.getW()),
        Rcpp::Named("A") = Rcpp::wrap(AList),
        Rcpp::Named("B") = Rcpp::wrap(BList),
        Rcpp::Named("objErr") = Rcpp::wrap(solver.objErr())
    );
    return output;
    }

// [[Rcpp::export(.onlineINMF_S1)]]
Rcpp::List onlineINMF_S1(Rcpp::List objectList, arma::uword k, const int& nCores,
    double lambda, arma::uword maxEpoch = 5, arma::uword minibatchSize = 5000,
    arma::uword maxHALSIter = 1, bool verbose = true) {
    if (Rf_isS4(objectList[0])) {
        return onlineINMF_S1_mem<arma::sp_mat>(Rcpp::as<std::vector<arma::sp_mat>>(objectList),
            k, nCores, lambda, maxEpoch, minibatchSize, maxHALSIter, verbose);
    } else {
        return onlineINMF_S1_mem<arma::mat>(Rcpp::as<std::vector<arma::mat>>(objectList),
            k, nCores, lambda, maxEpoch, minibatchSize, maxHALSIter, verbose);
    }
    return Rcpp::List::create();
}

// [[Rcpp::export(.onlineINMF_S1_h5dense)]]
Rcpp::List onlineINMF_S1_h5dense(std::vector<std::string> filenames,
    std::vector<std::string> dataPaths, arma::uword k, int nCores,
    double lambda, arma::uword maxEpoch = 5, arma::uword minibatchSize = 5000,
    arma::uword maxHALSIter = 1, bool verbose = true) {
    std::vector<std::unique_ptr<planc::H5Mat>> matPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i)
    {
        planc::H5Mat E(filenames[i], dataPaths[i]);
        std::unique_ptr<planc::H5Mat> ptr = std::make_unique<planc::H5Mat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    planc::ONLINEINMF<planc::H5Mat, arma::mat> solver(matPtrVec, k, lambda);
    solver.runOnlineINMF(minibatchSize, maxEpoch, maxHALSIter, verbose, nCores);

    Rcpp::List HList(filenames.size());
    Rcpp::List VList(filenames.size());
    Rcpp::List AList(filenames.size());
    Rcpp::List BList(filenames.size());
    for (arma::uword i = 0; i < filenames.size(); ++i)
    {
        HList[i] = solver.getHi(i);
        VList[i] = solver.getVi(i);
        AList[i] = solver.getAi(i);
        BList[i] = solver.getBi(i);
    }
    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("H") = Rcpp::wrap(HList),
        Rcpp::Named("V") = Rcpp::wrap(VList),
        Rcpp::Named("W") = Rcpp::wrap(solver.getW()),
        Rcpp::Named("A") = Rcpp::wrap(AList),
        Rcpp::Named("B") = Rcpp::wrap(BList),
        Rcpp::Named("objErr") = Rcpp::wrap(solver.objErr())
    );
    return output;
}

// [[Rcpp::export(.onlineINMF_S1_h5sparse)]]
Rcpp::List onlineINMF_S1_h5sparse(
    std::vector<std::string> filenames,
    std::vector<std::string> valuePaths,
    std::vector<std::string> rowindPaths,
    std::vector<std::string> colptrPaths,
    arma::uvec nrows, arma::uvec ncols,
    arma::uword k, const int& nCores, double lambda,
    arma::uword maxEpoch = 5, arma::uword minibatchSize = 5000,
    arma::uword maxHALSIter = 1, bool verbose = true
) {
    std::vector<std::unique_ptr<planc::H5SpMat>> matPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i)
    {
        planc::H5SpMat E(filenames[i], rowindPaths[i], colptrPaths[i], valuePaths[i], nrows[i], ncols[i]);
        std::unique_ptr<planc::H5SpMat> ptr = std::make_unique<planc::H5SpMat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    planc::ONLINEINMF<planc::H5SpMat, arma::sp_mat> solver(matPtrVec, k, lambda);
    solver.runOnlineINMF(minibatchSize, maxEpoch, maxHALSIter, verbose, nCores);

    Rcpp::List HList(filenames.size());
    Rcpp::List VList(filenames.size());
    Rcpp::List AList(filenames.size());
    Rcpp::List BList(filenames.size());
    for (arma::uword i = 0; i < filenames.size(); ++i)
    {
        HList[i] = solver.getHi(i);
        VList[i] = solver.getVi(i);
        AList[i] = solver.getAi(i);
        BList[i] = solver.getBi(i);
    }
    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("H") = Rcpp::wrap(HList),
        Rcpp::Named("V") = Rcpp::wrap(VList),
        Rcpp::Named("W") = Rcpp::wrap(solver.getW()),
        Rcpp::Named("A") = Rcpp::wrap(AList),
        Rcpp::Named("B") = Rcpp::wrap(BList),
        Rcpp::Named("objErr") = Rcpp::wrap(solver.objErr())
    );
    return output;
}

template <typename T>
Rcpp::List onlineINMF_S23_mem(std::vector<T> objectList,
    std::vector<arma::mat> Vinit, arma::mat Winit,
    std::vector<arma::mat> Ainit, std::vector<arma::mat> Binit,
    std::vector<T> objectListNew,
    arma::uword k, const int& nCores, double lambda, bool project = false, arma::uword maxEpoch = 5,
    arma::uword minibatchSize = 5000, arma::uword maxHALSIter = 1, bool verbose = true) {
    std::vector<std::unique_ptr<T>> matPtrVec = initMemMatPtr<T>(objectList);
    std::vector<std::unique_ptr<T>> matPtrVecNew = initMemMatPtr<T>(objectListNew);
    planc::ONLINEINMF<T, T> solver(matPtrVec, k, lambda);
    if (!project) solver.initV(Vinit, false);
    solver.initW(Winit, false);
    if (!project) solver.initA(Ainit);
    if (!project) solver.initB(Binit);
    solver.runOnlineINMF(matPtrVecNew, project, minibatchSize, maxEpoch, maxHALSIter, verbose, nCores);

    if (!project) {
        // Scenario 2
        arma::uword nDatasets = objectList.size() + objectListNew.size();
        Rcpp::List HList(nDatasets);
        Rcpp::List VList(nDatasets);
        Rcpp::List AList(nDatasets);
        Rcpp::List BList(nDatasets);
        for (arma::uword i = 0; i < nDatasets; ++i) {
            HList[i] = solver.getHi(i);
            VList[i] = solver.getVi(i);
            AList[i] = solver.getAi(i);
            BList[i] = solver.getBi(i);
        }
        Rcpp::List output = Rcpp::List::create(
            Rcpp::Named("H") = Rcpp::wrap(HList),
            Rcpp::Named("V") = Rcpp::wrap(VList),
            Rcpp::Named("W") = Rcpp::wrap(solver.getW()),
            Rcpp::Named("A") = Rcpp::wrap(AList),
            Rcpp::Named("B") = Rcpp::wrap(BList),
            Rcpp::Named("objErr") = Rcpp::wrap(solver.objErr())
        );
        return output;
    } else {
        // Scenario 3
        Rcpp::List HList(objectList.size());
        for (arma::uword i = 0; i < objectListNew.size(); ++i) {
            HList[i] = solver.getHi(i);
        }
        return Rcpp::List::create(
            Rcpp::Named("H") = HList
        );
    }
}

// [[Rcpp::export(.onlineINMF_S23)]]
Rcpp::List onlineINMF_S23(
    Rcpp::List objectList,
    const std::vector<arma::mat>& Vinit, const arma::mat& Winit,
    const std::vector<arma::mat>& Ainit, const std::vector<arma::mat>& Binit,
    const Rcpp::List& objectListNew,
    arma::uword k, const int& nCores, double lambda, bool project = false, arma::uword maxEpoch = 5,
    arma::uword minibatchSize = 5000, arma::uword maxHALSIter = 1, bool verbose = true) {
    if (Rf_isS4(objectList[0])) {
        return onlineINMF_S23_mem<arma::sp_mat>(
            Rcpp::as<std::vector<arma::sp_mat>>(objectList),
            Vinit, Winit, Ainit, Binit,
            Rcpp::as<std::vector<arma::sp_mat>>(objectListNew),
            k, nCores,lambda, project, maxEpoch, minibatchSize, maxHALSIter, verbose);
    } else {
        return onlineINMF_S23_mem<arma::mat>(
            Rcpp::as<std::vector<arma::mat>>(objectList),
            Vinit, Winit, Ainit, Binit,
            Rcpp::as<std::vector<arma::mat>>(objectListNew),
            k, nCores, lambda, project, maxEpoch, minibatchSize, maxHALSIter, verbose);
    }
    return Rcpp::List::create();
}

// [[Rcpp::export(.onlineINMF_S23_h5dense)]]
Rcpp::List onlineINMF_S23_h5dense(
    std::vector<std::string> filenames, std::vector<std::string> dataPaths,
    std::vector<std::string> filenamesNew, std::vector<std::string> dataPathsNew,
    std::vector<arma::mat> Vinit, const arma::mat& Winit,
    std::vector<arma::mat> Ainit, std::vector<arma::mat> Binit,
    arma::uword k, const int& nCores, double lambda, bool project = false, arma::uword maxEpoch = 5,
    arma::uword minibatchSize = 5000, arma::uword maxHALSIter = 1, bool verbose = true) {
    std::vector<std::unique_ptr<planc::H5Mat>> matPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i)
    {
        planc::H5Mat E(filenames[i], dataPaths[i]);
        std::unique_ptr<planc::H5Mat> ptr = std::make_unique<planc::H5Mat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    std::vector<std::unique_ptr<planc::H5Mat>> matPtrVecNew;
    for (arma::uword i = 0; i < filenamesNew.size(); ++i)
    {
        planc::H5Mat E(filenamesNew[i], dataPathsNew[i]);
        std::unique_ptr<planc::H5Mat> ptr = std::make_unique<planc::H5Mat>(E);
        matPtrVecNew.push_back(std::move(ptr));
    }
    planc::ONLINEINMF<planc::H5Mat, arma::mat> solver(matPtrVec, k, lambda);
    solver.initV(Vinit, false);
    solver.initW(Winit, false);
    solver.initA(Ainit);
    solver.initB(Binit);
    solver.runOnlineINMF(matPtrVecNew, project, minibatchSize, maxEpoch, maxHALSIter, verbose, nCores);
    if (!project) {
        // Scenario 2
        arma::uword nDatasets = filenames.size() + filenamesNew.size();
        Rcpp::List HList(nDatasets);
        Rcpp::List VList(nDatasets);
        Rcpp::List AList(nDatasets);
        Rcpp::List BList(nDatasets);
        for (arma::uword i = 0; i < nDatasets; ++i) {
            HList[i] = solver.getHi(i);
            VList[i] = solver.getVi(i);
            AList[i] = solver.getAi(i);
            BList[i] = solver.getBi(i);
        }
        Rcpp::List output = Rcpp::List::create(
            Rcpp::Named("H") = Rcpp::wrap(HList),
            Rcpp::Named("V") = Rcpp::wrap(VList),
            Rcpp::Named("W") = Rcpp::wrap(solver.getW()),
            Rcpp::Named("A") = Rcpp::wrap(AList),
            Rcpp::Named("B") = Rcpp::wrap(BList),
            Rcpp::Named("objErr") = Rcpp::wrap(solver.objErr())
        );
        return output;
    } else {
        // Scenario 3
        Rcpp::List HList(filenamesNew.size());
        for (arma::uword i = 0; i < filenamesNew.size(); ++i) {
            HList[i] = solver.getHi(i);
        }
        return Rcpp::List::create(
            Rcpp::Named("H") = HList
        );
    }
}

// [[Rcpp::export(.onlineINMF_S23_h5sparse)]]
Rcpp::List onlineINMF_S23_h5sparse(
    std::vector<std::string> filenames, std::vector<std::string> valuePaths,
    std::vector<std::string> rowindPaths, std::vector<std::string> colptrPaths,
    arma::uvec nrows, arma::uvec ncols,
    std::vector<std::string> filenamesNew, std::vector<std::string> valuePathsNew,
    std::vector<std::string> rowindPathsNew, std::vector<std::string> colptrPathsNew,
    arma::uvec nrowsNew, arma::uvec ncolsNew,
    std::vector<arma::mat> Vinit, const arma::mat& Winit,
    std::vector<arma::mat> Ainit, std::vector<arma::mat> Binit,
    arma::uword k, const int& nCores, double lambda, bool project = false, arma::uword maxEpoch = 5,
    arma::uword minibatchSize = 5000, arma::uword maxHALSIter = 1, bool verbose = true) {
    std::vector<std::unique_ptr<planc::H5SpMat>> matPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i)
    {
        planc::H5SpMat E(filenames[i], rowindPaths[i], colptrPaths[i], valuePaths[i], nrows[i], ncols[i]);
        std::unique_ptr<planc::H5SpMat> ptr = std::make_unique<planc::H5SpMat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    std::vector<std::unique_ptr<planc::H5SpMat>> matPtrVecNew;
    for (arma::uword i = 0; i < filenamesNew.size(); ++i)
    {
        planc::H5SpMat E(filenamesNew[i], rowindPathsNew[i], colptrPathsNew[i], valuePathsNew[i], nrowsNew[i], ncolsNew[i]);
        std::unique_ptr<planc::H5SpMat> ptr = std::make_unique<planc::H5SpMat>(E);
        matPtrVecNew.push_back(std::move(ptr));
    }
    planc::ONLINEINMF<planc::H5SpMat, arma::sp_mat> solver(matPtrVec, k, lambda);
    solver.initV(Vinit, false);
    solver.initW(Winit, false);
    solver.initA(Ainit);
    solver.initB(Binit);
    solver.runOnlineINMF(matPtrVecNew, project, minibatchSize, maxEpoch, maxHALSIter, verbose, nCores);
    if (!project) {
        // Scenario 2
        arma::uword nDatasets = filenames.size() + filenamesNew.size();
        Rcpp::List HList(nDatasets);
        Rcpp::List VList(nDatasets);
        Rcpp::List AList(nDatasets);
        Rcpp::List BList(nDatasets);
        for (arma::uword i = 0; i < nDatasets; ++i) {
            HList[i] = solver.getHi(i);
            VList[i] = solver.getVi(i);
            AList[i] = solver.getAi(i);
            BList[i] = solver.getBi(i);
        }
        Rcpp::List output = Rcpp::List::create(
            Rcpp::Named("H") = Rcpp::wrap(HList),
            Rcpp::Named("V") = Rcpp::wrap(VList),
            Rcpp::Named("W") = Rcpp::wrap(solver.getW()),
            Rcpp::Named("A") = Rcpp::wrap(AList),
            Rcpp::Named("B") = Rcpp::wrap(BList),
            Rcpp::Named("objErr") = Rcpp::wrap(solver.objErr())
        );
        return output;
    } else {
        // Scenario 3
        Rcpp::List HList(filenamesNew.size());
        for (arma::uword i = 0; i < filenamesNew.size(); ++i) {
            HList[i] = solver.getHi(i);
        }
        return Rcpp::List::create(
            Rcpp::Named("H") = HList
        );
    }
}
// %%%%%%%%%%%%%% UINMF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template <typename T>
Rcpp::List uinmf_mem(std::vector<T> objectList,
                    std::vector<T> unsharedList,
                    std::vector<int> whichUnshared,
                    arma::uword k, const int& nCores, arma::vec lambda,
                    arma::uword niter, bool verbose)
{
    std::vector<std::unique_ptr<T>> matPtrVec;
    std::vector<std::unique_ptr<T>> unsharedPtrVec;
    matPtrVec = initMemMatPtr<T>(objectList);
    unsharedPtrVec = initMemMatPtr<T>(unsharedList);
    planc::UINMF<T> solver(matPtrVec, unsharedPtrVec, whichUnshared, k, lambda);
    solver.optimizeUANLS(niter, verbose, nCores);

    Rcpp::List HList(objectList.size());
    Rcpp::List VList(objectList.size());
    Rcpp::List UList(unsharedList.size());
    for (arma::uword i = 0; i < objectList.size(); ++i)
    {
        HList[i] = solver.getHi(i);
        VList[i] = solver.getVi(i);
        int uidx = whichUnshared[i];
        if (uidx >= 0) {
            UList[uidx] = solver.getUi(uidx);
        }
    }
    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("H") = Rcpp::wrap(HList),
        Rcpp::Named("V") = Rcpp::wrap(VList),
        Rcpp::Named("W") = Rcpp::wrap(solver.getW()),
        Rcpp::Named("U") = Rcpp::wrap(UList),
        Rcpp::Named("objErr") = Rcpp::wrap(solver.objErr())
    );
    return output;
}

// [[Rcpp::export(.uinmf_rcpp)]]
Rcpp::List uinmf_rcpp(Rcpp::List objectList, const Rcpp::List& unsharedList,
                      const std::vector<int>& whichUnshared, arma::uword k, const int& nCores,
                      const arma::vec& lambda, arma::uword niter, bool verbose) {
    if (Rf_isS4(objectList[0])) {
        return uinmf_mem<arma::sp_mat>(
            Rcpp::as<std::vector<arma::sp_mat>>(objectList),
            Rcpp::as<std::vector<arma::sp_mat>>(unsharedList),
            whichUnshared, k, nCores, lambda, niter, verbose
        );
    } else {
        return uinmf_mem<arma::mat>(
            Rcpp::as<std::vector<arma::mat>>(objectList),
            Rcpp::as<std::vector<arma::mat>>(unsharedList),
            whichUnshared, k, nCores, lambda, niter, verbose
        );
    }
    return Rcpp::List::create();
}

// [[Rcpp::export(.uinmf_h5dense)]]
Rcpp::List uinmf_h5dense(std::vector<std::string> filenames,
                         std::vector<std::string> dataPaths,
                         std::vector<std::string> unsharedFilenames,
                         std::vector<std::string> unsharedDataPaths,
                         std::vector<int> whichUnshared,
                         arma::uword k, const int& nCores, const arma::vec& lambda,
                         arma::uword niter, bool verbose) {
    std::vector<std::unique_ptr<planc::H5Mat>> matPtrVec;
    std::vector<std::unique_ptr<planc::H5Mat>> unsharedPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        planc::H5Mat E(filenames[i], dataPaths[i]);
        std::unique_ptr<planc::H5Mat> ptr = std::make_unique<planc::H5Mat>(E);
        matPtrVec.push_back(std::move(ptr));

        planc::H5Mat E_unshared(unsharedFilenames[i], unsharedDataPaths[i]);
        std::unique_ptr<planc::H5Mat> ptr_unshared = std::make_unique<planc::H5Mat>(E_unshared);
        unsharedPtrVec.push_back(std::move(ptr_unshared));
    }
    planc::UINMF<planc::H5Mat> solver(matPtrVec, unsharedPtrVec, whichUnshared, k, lambda);

    solver.optimizeUANLS(niter, verbose, nCores);

    Rcpp::List HList(matPtrVec.size());
    Rcpp::List VList(matPtrVec.size());
    Rcpp::List UList(matPtrVec.size());
    for (arma::uword i = 0; i < matPtrVec.size(); ++i)
    {
        HList[i] = solver.getHi(i);
        VList[i] = solver.getVi(i);
        int uidx = whichUnshared[i];
        if (uidx >= 0) {
            UList[uidx] = solver.getUi(uidx);
        }
    }
    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("H") = Rcpp::wrap(HList),
        Rcpp::Named("V") = Rcpp::wrap(VList),
        Rcpp::Named("W") = Rcpp::wrap(solver.getW()),
        Rcpp::Named("U") = Rcpp::wrap(UList),
        Rcpp::Named("objErr") = Rcpp::wrap(solver.objErr())
    );
    return output;
}

// [[Rcpp::export(.uinmf_h5sparse)]]
Rcpp::List uinmf_h5sparse(std::vector<std::string> filenames,
                          std::vector<std::string> rowindPaths,
                          std::vector<std::string> colptrPaths,
                          std::vector<std::string> valuePaths,
                          arma::uvec nrows, arma::uvec ncols,
                          std::vector<std::string> unsharedFilenames,
                          std::vector<std::string> unsharedRowindPaths,
                          std::vector<std::string> unsharedColptrPaths,
                          std::vector<std::string> unsharedValuePaths,
                          arma::uvec unsharedNrows, arma::uvec unsharedNcols,
                          std::vector<int> whichUnshared,
                          arma::uword k, const int& nCores, const arma::vec& lambda,
                          arma::uword niter, bool verbose) {
    std::vector<std::unique_ptr<planc::H5SpMat>> matPtrVec;
    std::vector<std::unique_ptr<planc::H5SpMat>> unsharedPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        planc::H5SpMat E(filenames[i], rowindPaths[i], colptrPaths[i], valuePaths[i], nrows[i], ncols[i]);
        std::unique_ptr<planc::H5SpMat> ptr = std::make_unique<planc::H5SpMat>(E);
        matPtrVec.push_back(std::move(ptr));

        planc::H5SpMat E_unshared(unsharedFilenames[i], unsharedRowindPaths[i], unsharedColptrPaths[i], unsharedValuePaths[i], unsharedNrows[i], unsharedNcols[i]);
        std::unique_ptr<planc::H5SpMat> ptr_unshared = std::make_unique<planc::H5SpMat>(E_unshared);
        unsharedPtrVec.push_back(std::move(ptr_unshared));
    }
    planc::UINMF<planc::H5SpMat> solver(matPtrVec, unsharedPtrVec, whichUnshared, k, lambda);
    solver.optimizeUANLS(niter, verbose, nCores);
    Rcpp::List HList(matPtrVec.size());
    Rcpp::List VList(matPtrVec.size());
    Rcpp::List UList(matPtrVec.size());
    for (arma::uword i = 0; i < matPtrVec.size(); ++i)
    {
        HList[i] = solver.getHi(i);
        VList[i] = solver.getVi(i);
        int uidx = whichUnshared[i];
        if (uidx >= 0) {
            UList[uidx] = solver.getUi(uidx);
        }
    }
    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("H") = Rcpp::wrap(HList),
        Rcpp::Named("V") = Rcpp::wrap(VList),
        Rcpp::Named("W") = Rcpp::wrap(solver.getW()),
        Rcpp::Named("U") = Rcpp::wrap(UList),
        Rcpp::Named("objErr") = Rcpp::wrap(solver.objErr())
    );
    return output;
}

// [[Rcpp::export(.testCacheCalc)]]
arma::uword testcacheCalc(int rank) {
    return chunk_size_dense<double>(rank);
}

// [[Rcpp::export(.getBoundThreads)]]
arma::uword getBoundThreadCount() {
    return get_num_bound_threads();
}

// [[Rcpp::export(.openblaspthreadoff)]]
void openblas_pthread_off(Rcpp::XPtr<void*> libloc) {
    if (is_openmp()) {
        if (const std::function<int()> openblas_parallel = get_openblas_parallel(libloc))
        {
            if (openblas_parallel() == 1) {
                const std::function<void(int)> openblas_set = get_openblas_set(libloc);
                openblas_set(1);
            }
        }
    }
}

// [[Rcpp::export(.openblaspthreadon)]]
void openblas_pthread_on(Rcpp::XPtr<void*> libloc) {
    if (is_openmp()) {
        if (const std::function<int()> openblas_parallel = get_openblas_parallel(libloc))
        {
            if (openblas_parallel() == 1) {
                const std::function<void(int)> openblas_set = get_openblas_set(libloc);
                openblas_set(0);
            }
        }
    }
}