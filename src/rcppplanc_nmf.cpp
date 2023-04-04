// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "bppnmf.hpp"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//
//
//
// Use the ANLS-BPP algorithm to factor a given matrix
// at the given rank.
//
// [[Rcpp::export]]
Rcpp::List rcppplanc_bppnmf(const arma::sp_mat & x, int & k) {
  planc::BPPNMF<arma::sp_mat> MyNMF = planc::BPPNMF<arma::sp_mat>(x, k);
  MyNMF.computeNMF();
  return Rcpp::List::create(
      Rcpp::Named("W") = MyNMF.getLeftLowRankFactor(),
      Rcpp::Named("H") = MyNMF.getRightLowRankFactor()
    );
}
