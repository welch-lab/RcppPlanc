// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "bppnmf.hpp"
#include "bppinmf.hpp"
#include "bppnnls.hpp"
#include "aoadmm.hpp"
#include "gnsym.hpp"
#include "hals.hpp"
#include "mu.hpp"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

normtype m_input_normalization;
int m_initseed;
arma::fvec m_regW;
arma::fvec m_regH;

template <class T>
Rcpp::List RcallNMF(arma::sp_mat x, int k, int niter)
{
  arma::uword m_m = x.n_rows;
  arma::uword m_n = x.n_cols;
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
  MyNMF.symm_reg(-1);
  // MyNMF.compute_error(this->m_compute_error);

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

template <class T>
Rcpp::List RcallNMF(arma::sp_mat x, int k, int niter,
                    arma::mat Winit, arma::mat Hinit) {
  T MyNMF(x, Winit, Hinit);
  MyNMF.num_iterations(niter);
  MyNMF.symm_reg(-1);
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
//' @param H_init Initial right-hand factor matrix (Optional)
//' @param W_init Initial left-hand factor matrix (Optional)
//' @export
//' @returns The calculated factor matrices as an Rcpp::List
//' @examplesIf require("Matrix")
//' aoadmmnmf(rsparsematrix(nrow = 100, ncol = 100, nnz = 10, symmetric = TRUE), 10, 10)
// [[Rcpp::export]]
Rcpp::List aoadmmnmf(const arma::sp_mat &x, const int &k, const int &niter,
                     const Rcpp::Nullable<Rcpp::NumericMatrix> &W_init = R_NilValue,
                     const Rcpp::Nullable<Rcpp::NumericMatrix> &H_init = R_NilValue) {
  Rcpp::List outlist;
  if (W_init.isNotNull() && W_init.isNotNull()) {
    outlist = RcallNMF<planc::AOADMMNMF<arma::sp_mat>>(x, k, niter,
                                                       Rcpp::as<arma::mat>(W_init),
                                                       Rcpp::as<arma::mat>(H_init));
  }
  else {
    outlist = RcallNMF<planc::AOADMMNMF<arma::sp_mat>>(x, k, niter);
  }
  return outlist;
}

//' Gauss-Newton using Conjugate Gradients NMF
//'
//' Use the Gauss-Newton algorithm to factor a given matrix at the given rank.
//'
//' @param x Input matrix for factorization
//' @param k Factor matrix rank
//' @param niter Maximum number of nmf iterations
//' @param H_init Initial right-hand factor matrix (Optional)
//' @param W_init Initial left-hand factor matrix (Optional)
//' @export
//' @returns The calculated factor matrices as an Rcpp::List
//' @examplesIf require("Matrix")
//' gnsymnmf(rsparsematrix(nrow = 100, ncol = 100, nnz = 10, symmetric = TRUE), 10, 10)
// [[Rcpp::export]]
Rcpp::List gnsymnmf(const arma::sp_mat &x, const int &k, const int &niter,
                    const Rcpp::Nullable<Rcpp::NumericMatrix> &W_init = R_NilValue,
                    const Rcpp::Nullable<Rcpp::NumericMatrix> &H_init = R_NilValue)
{
  Rcpp::List outlist;
  if (W_init.isNotNull() && W_init.isNotNull()) {
    outlist = RcallNMF<planc::GNSYMNMF<arma::sp_mat>>(x, k, niter,
                                                      Rcpp::as<arma::mat>(W_init),
                                                      Rcpp::as<arma::mat>(H_init));
  }
  else {
    outlist = RcallNMF<planc::GNSYMNMF<arma::sp_mat>>(x, k, niter);
  }
  return outlist;
}

//' Hierarchical Alternating Least Squares NMF
//'
//' Use the HALS algorithm to factor a given matrix at the given rank.
//'
//' @param x Input matrix for factorization
//' @param k Factor matrix rank
//' @param niter Maximum number of nmf iterations
//' @param H_init Initial right-hand factor matrix (Optional)
//' @param W_init Initial left-hand factor matrix (Optional)
//' @export
//' @returns The calculated factor matrices as an Rcpp::List
//' @examplesIf require("Matrix")
//' halsmnmf(rsparsematrix(nrow = 100, ncol = 100, nnz = 10, symmetric = TRUE), 10, 10)
// [[Rcpp::export]]
Rcpp::List halsnmf(const arma::sp_mat &x, const int &k, const int &niter,
                   const Rcpp::Nullable<Rcpp::NumericMatrix> &W_init = R_NilValue,
                   const Rcpp::Nullable<Rcpp::NumericMatrix> &H_init = R_NilValue)
{
  Rcpp::List outlist;
  if (W_init.isNotNull() && W_init.isNotNull()) {
    outlist = RcallNMF<planc::HALSNMF<arma::sp_mat>>(x, k, niter,
                                                     Rcpp::as<arma::mat>(W_init),
                                                     Rcpp::as<arma::mat>(H_init));
  }
  else {
    outlist = RcallNMF<planc::HALSNMF<arma::sp_mat>>(x, k, niter);
  }
  return outlist;
}
//' Multiplicative Update NMF
//'
//' Use the MU algorithm to factor a given matrix
//' at the given rank.
//' @param x Input matrix for factorization
//' @param k Factor matrix rank
//' @param niter Maximum number of nmf iterations
//' @param H_init Initial right-hand factor matrix (Optional)
//' @param W_init Initial left-hand factor matrix (Optional)
//' @export
//' @returns The calculated factor matrices as an Rcpp::List
//' @examplesIf require("Matrix")
//' halsnmf(rsparsematrix(nrow = 100, ncol = 100, nnz = 10, symmetric = TRUE), 10, 10)
// [[Rcpp::export]]
Rcpp::List munmf(const arma::sp_mat &x, const int &k, const int &niter,
                 const Rcpp::Nullable<Rcpp::NumericMatrix> &W_init = R_NilValue,
                 const Rcpp::Nullable<Rcpp::NumericMatrix> &H_init = R_NilValue)
{
  Rcpp::List outlist;
  if (W_init.isNotNull() && W_init.isNotNull()) {
    outlist = RcallNMF<planc::MUNMF<arma::sp_mat>>(x, k, niter,
                                                   Rcpp::as<arma::mat>(W_init),
                                                   Rcpp::as<arma::mat>(H_init));
  }
  else {
    outlist = RcallNMF<planc::MUNMF<arma::sp_mat>>(x, k, niter);
  }
  return outlist;
}
//' Alternating  Nonnegative Least Squares with Block Principal Pivoting NMF
//'
//' Use the ANLS-BPP algorithm to factor a given matrix at the given rank.
//'
//' @param x Input matrix for factorization
//' @param k Factor matrix rank
//' @param niter Maximum number of nmf iterations
//' @param H_init Initial right-hand factor matrix (Optional)
//' @param W_init Initial left-hand factor matrix (Optional)
//' @export
//' @returns The calculated factor matrices as an Rcpp::List
//' @examplesIf require("Matrix")
//' munmf(rsparsematrix(nrow = 100, ncol = 100, nnz = 10, symmetric = TRUE), 10, 10)
// [[Rcpp::export]]
Rcpp::List bppnmf(const arma::sp_mat & x, const int & k, const int & niter,
                  const Rcpp::Nullable<Rcpp::NumericMatrix> &W_init = R_NilValue,
                  const Rcpp::Nullable<Rcpp::NumericMatrix> &H_init = R_NilValue) {
  Rcpp::List outlist;
  if (W_init.isNotNull() && W_init.isNotNull()) {
    outlist = RcallNMF<planc::BPPNMF<arma::sp_mat>>(x, k, niter,
                                                    Rcpp::as<arma::mat>(W_init),
                                                    Rcpp::as<arma::mat>(H_init));
  }
  else {
    outlist = RcallNMF<planc::BPPNMF<arma::sp_mat>>(x, k, niter);
  }
  return outlist;
}

//' Block Principal Pivoted Non-Negative Least Squares
//'
//' Use the BPP algorithm to get the nonnegative least squares solution for the given matrices.
//'
//' @param C Input factor dense matrix
//' @param B Input sparse matrix
//' @export
//' @returns The calculated solution matrix in dense form.
// [[Rcpp::export]]
arma::mat bppnnls(const arma::mat &C, const arma::sp_mat &B) {
    arma::uword m_n = B.n_cols;
    arma::uword m_k = C.n_cols;
    arma::mat outmat = arma::zeros<arma::mat>(m_k, m_n);
    arma::mat *outmatptr;
    outmatptr = &outmat;
    unsigned int numChunks = m_n / ONE_THREAD_MATRIX_SIZE;
    if (numChunks*ONE_THREAD_MATRIX_SIZE < m_n) numChunks++;
#pragma omp parallel for schedule(auto)
    for (unsigned int i = 0; i < numChunks; i++) {
        unsigned int spanStart = i * ONE_THREAD_MATRIX_SIZE;
        unsigned int spanEnd = (i + 1) * ONE_THREAD_MATRIX_SIZE - 1;
        if (spanEnd > m_n - 1) spanEnd = m_n - 1;
        // double start = omp_get_wtime();
        BPPNNLS<arma::mat, arma::vec> solveProblem(C, (arma::mat)B.cols(spanStart, spanEnd));
        solveProblem.solveNNLS();
        // double end = omp_get_wtime();
      // titer = end - start;
      // #ifdef _VERBOSE
      // INFO << " start=" << spanStart
      //     << ", end=" << spanEnd
      //     << ", tid=" << omp_get_thread_num() << " cpu=" << sched_getcpu() << std::endl;
      //     // << " time taken=" << titer << std::endl;
      // #endif
        (*outmatptr).cols(spanStart, spanEnd) = solveProblem.getSolutionMatrix();
    };

    return outmat;
}

//' Block Principal Pivoted Iterative Non-Negative Matrix Factorization
//'
//' Use the BPP algorithm to iteratively factor the given datasets.
//'
//' @param Ei List of datasets in dense matrix form.
//' @export
//' @returns The calculated solution matrix in dense form.
//' @examplesIf require("Matrix")
//' bppinmf(rsparsematrix(nrow=20,ncol=20,nnz=10), Matrix(runif(n=200,min=0,max=2),20,10))
// [[Rcpp::export]]
Rcpp::List bppinmf(std::vector<arma::mat> objectList, arma::uword k, double lambda, arma::uword maxIter, double thresh) {
    // std::vector<arma::mat> matVec;
    // std::vector<std::unique_ptr<arma::mat>> matPtrVec;
    // for (arma::uword i = 0; i < objectList.size(); ++i) {
    //     matVec.push_back(arma::mat(objectList[i].begin(), objectList[i].nrow(), objectList[i].ncol(), false));
    // }
    // for (arma::uword i = 0; i < objectList.size(); ++i)
    // {
    //     matPtrVec.push_back(std::unique_ptr<arma::mat>(&matVec[i]));
    // }


    std::vector<std::unique_ptr<arma::mat>> matPtrVec;
    for (arma::uword i = 0; i < objectList.size(); ++i)
    {
        arma::mat E = arma::mat(objectList[i].begin(), objectList[i].n_rows, objectList[i].n_cols, false, true);
        std::unique_ptr<arma::mat> ptr = std::make_unique<arma::mat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    planc::BPPINMF<arma::mat> solver(matPtrVec, k, lambda);
    solver.optimizeALS(maxIter, thresh);
    std::cout << "iNMF finished" << std::endl;
    Rcpp::List HList = Rcpp::List::create();
    Rcpp::List VList = Rcpp::List::create();
    for (arma::uword i = 0; i < objectList.size(); ++i) {
        HList.push_back(solver.getHi(i));
        VList.push_back(solver.getVi(i));
    }
    // return HList;
    return Rcpp::List::create(
        Rcpp::Named("H") = HList,
        Rcpp::Named("V") = VList,
        Rcpp::Named("W") = solver.getW()
    );
}

//' @export
// [[Rcpp::export]]
Rcpp::List bppinmf_sparse(Rcpp::List objectList, arma::uword k, double lambda) {
    std::cout << "Testing bppinmf" << std::endl;
    std::vector<std::unique_ptr<arma::sp_mat>> matPtrVec;
    for (arma::uword i = 0; i < objectList.size(); ++i)
    {
        std::cout << "i=" << i << std::endl;
        arma::sp_mat E = arma::sp_mat(objectList[i]);
        std::unique_ptr<arma::sp_mat> ptr = std::make_unique<arma::sp_mat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    std::cout << "matPtrVec size=" << matPtrVec.size() << std::endl;
    planc::BPPINMF<arma::sp_mat> solver(matPtrVec, k, lambda);
    std::cout << "solver created" << std::endl;
    solver.optimizeALS(5u, 0.000001);
    Rcpp::List HList = Rcpp::List::create();
    Rcpp::List VList = Rcpp::List::create();
    for (arma::uword i = 0; i < objectList.size(); ++i) {
        HList.push_back(solver.getHi(i));
        VList.push_back(solver.getVi(i));
    }
    return Rcpp::List::create(
        Rcpp::Named("H") = HList,
        Rcpp::Named("V") = VList,
        Rcpp::Named("W") = solver.getW()
    );
}
