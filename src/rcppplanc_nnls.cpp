//
// Created by andrew on 3/21/2025.
//
#ifndef ARMA_64BIT_WORD
#define ARMA_64BIT_WORD
#endif
#include <RcppArmadillo.h>
#include "config.h"
#include <progress.hpp>
#include <nnls_lib.hpp>


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
         return planc::nnlslib<arma::sp_mat>::runbppnnls(C, Rcpp::as<arma::sp_mat>(B), nCores);
     } else {
         return planc::nnlslib<arma::mat>::runbppnnls(C, Rcpp::as<arma::mat>(B), nCores);
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
     return planc::nnlslib<arma::mat>::bppnnls_prod(CtC, CtB, nCores);
 }