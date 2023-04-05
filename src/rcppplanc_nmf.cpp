// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "bppnmf.hpp"

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

// Use the ANLS-BPP algorithm to factor a given matrix
// at the given rank. TODO FIX set.seed
//
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
