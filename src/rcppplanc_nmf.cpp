// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#define ARMA_DONT_PRINT_FAST_MATH_WARNING
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <progress.hpp>
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
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

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
//' @returns The calculated factor matrices as an Rcpp::List
//' @examplesIf require("Matrix")
//' halsnmf(rsparsematrix(nrow = 100, ncol = 100, nnz = 10, symmetric = TRUE), 10, 10)
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

template <typename T>
arma::mat runbppnnls(const arma::mat &C, const T &B) {
    arma::uword m_n = B.n_cols;
    arma::uword m_k = C.n_cols;
    arma::mat CtC = C.t() * C;
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
        arma::mat CtBChunk = C.t() * B.cols(spanStart, spanEnd);
        BPPNNLS<arma::mat, arma::vec> solveProblem(CtC, CtBChunk, true);
        solveProblem.solveNNLS();
        (*outmatptr).cols(spanStart, spanEnd) = solveProblem.getSolutionMatrix();
    };
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
//' @returns The calculated solution matrix in dense form.
//' @rdname bppnnls
//' @examples
//' set.seed(1)
//' C <- matrix(rnorm(1000), nrow = 100)
//' B <- matrix(rnorm(1500), nrow = 100)
//' res1 <- bppnnls(C, B)
//' dim(res1)
//' res2 <- bppnnls_prod(t(C) %*% C, t(C) %*% B)
//' all.equal(res1, res2)
// [[Rcpp::export]]
arma::mat bppnnls(const arma::mat &C, const SEXP &B) {
    if (Rf_isS4(B)) {
        return runbppnnls<arma::sp_mat>(C, Rcpp::as<arma::sp_mat>(B));
    } else {
        return runbppnnls<arma::mat>(C, Rcpp::as<arma::mat>(B));
    }
    return arma::mat();
}

//' @param CtC The \eqn{C^\mathsf{T}C} matrix, see description.
//' @param CtB The \eqn{C^\mathsf{T}B} matrix, see description.
//' @rdname bppnnls
// [[Rcpp::export]]
arma::mat bppnnls_prod(const arma::mat &CtC, const arma::mat &CtB) {
    arma::uword n = CtB.n_cols;
    arma::uword k = CtC.n_cols;
    arma::mat outmat = arma::zeros<arma::mat>(k, n);
    arma::mat *outmatptr;
    outmatptr = &outmat;
    unsigned int numChunks = n / ONE_THREAD_MATRIX_SIZE;
    if (numChunks*ONE_THREAD_MATRIX_SIZE < n) numChunks++;
#pragma omp parallel for schedule(auto)
    for (unsigned int i = 0; i < numChunks; i++) {
        unsigned int spanStart = i * ONE_THREAD_MATRIX_SIZE;
        unsigned int spanEnd = (i + 1) * ONE_THREAD_MATRIX_SIZE - 1;
        if (spanEnd > n - 1) spanEnd = n - 1;
        arma::mat CtBChunk = CtB.cols(spanStart, spanEnd);
        BPPNNLS<arma::mat, arma::vec> solveProblem(CtC, CtBChunk, true);
        solveProblem.solveNNLS();
        (*outmatptr).cols(spanStart, spanEnd) = solveProblem.getSolutionMatrix();
    };
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
    return std::move(matPtrVec);
}

template <typename T>
Rcpp::List runINMF(std::vector<T> objectList, arma::uword k, double lambda,
                   arma::uword niter, bool verbose)
{
    std::vector<std::unique_ptr<T>> matPtrVec;
    matPtrVec = initMemMatPtr<T>(objectList);
    planc::BPPINMF<T> solver(matPtrVec, k, lambda);
    solver.initW();
    solver.initV();
    solver.initH();
    solver.optimizeALS(niter, verbose);
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
                   std::vector<arma::mat> HinitList, std::vector<arma::mat> VinitList, arma::mat Winit)
{
    std::vector<std::unique_ptr<T>> matPtrVec;
    matPtrVec = initMemMatPtr<T>(objectList);
    planc::BPPINMF<T> solver(matPtrVec, k, lambda);
    solver.initW(Winit);
    solver.initV(VinitList);
    solver.initH(HinitList);
    solver.optimizeALS(niter, verbose);
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

Rcpp::List bppinmf_dense(std::vector<arma::mat> objectList, arma::uword k,
                         double lambda, arma::uword niter, bool verbose = true,
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
                                     Rcpp::as<arma::mat>(Winit));
    }
    else
    {
        return runINMF<arma::mat>(objectList, k, lambda,
                                     niter, verbose);
    }
}

Rcpp::List bppinmf_sparse(std::vector<arma::sp_mat> objectList, arma::uword k, double lambda,
    arma::uword niter, bool verbose = true,
    Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
    Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
    Rcpp::Nullable<arma::mat> Winit  = R_NilValue) {
    if (Hinit.isNotNull() && Vinit.isNotNull() && Winit.isNotNull()) {
        return runINMF<arma::sp_mat>(objectList, k, lambda,
                                     niter, verbose,
                                     Rcpp::as<std::vector<arma::mat>>(Hinit),
                                     Rcpp::as<std::vector<arma::mat>>(Vinit),
                                     Rcpp::as<arma::mat>(Winit));
    }
    else {
        return runINMF<arma::sp_mat>(objectList, k, lambda,
                niter, verbose);
    }
}

// [[Rcpp::export(.bppinmf)]]
Rcpp::List bppinmf(Rcpp::List objectList, const arma::uword k,
                   const double lambda = 5, const arma::uword niter = 30,
                   const bool verbose = true,
                   Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
                   Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
                   Rcpp::Nullable<arma::mat> Winit = R_NilValue) {
    if (Rf_isS4(objectList[0])) {
        return bppinmf_sparse(Rcpp::as<std::vector<arma::sp_mat>>(objectList), k, lambda,
                        niter, verbose, Hinit, Vinit, Winit);
    } else {
        return bppinmf_dense(Rcpp::as<std::vector<arma::mat>>(objectList), k, lambda,
                            niter, verbose, Hinit, Vinit, Winit);
    }
    return Rcpp::List::create();
}

// [[Rcpp::export(.bppinmf_h5dense)]]
Rcpp::List bppinmf_h5dense(std::vector<std::string> filenames, std::vector<std::string> dataPath,
    arma::uword k, double lambda, arma::uword niter, bool verbose = true,
    Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
    Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
    Rcpp::Nullable<arma::mat> Winit  = R_NilValue) {
    std::vector<std::unique_ptr<planc::H5Mat>> matPtrVec;
    for (int i = 0; i < filenames.size(); ++i) {
        planc::H5Mat h5m(filenames[i], dataPath[i]);
        std::unique_ptr<planc::H5Mat> ptr = std::make_unique<planc::H5Mat>(h5m);
        matPtrVec.push_back(std::move(ptr));
    }
    planc::BPPINMF<planc::H5Mat> solver(matPtrVec, k, lambda);

    if (Winit.isNotNull()) {
        arma::mat W = Rcpp::as<arma::mat>(Winit);
        solver.initW(W);
    } else {
        solver.initW();
    }

    if (Vinit.isNotNull()) {
        std::vector<arma::mat> VinitList = Rcpp::as<std::vector<arma::mat>>(Vinit);
        solver.initV(VinitList);
    } else {
        solver.initV();
    }

    if (Hinit.isNotNull()) {
        std::vector<arma::mat> HinitList = Rcpp::as<std::vector<arma::mat>>(Hinit);
        solver.initH(HinitList);
    } else {
        solver.initH();
    }

    solver.optimizeALS(niter, verbose);
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
    arma::uword k, double lambda, arma::uword niter,
    bool verbose = true,
    Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
    Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
    Rcpp::Nullable<arma::mat> Winit  = R_NilValue) {
    std::vector<std::unique_ptr<planc::H5SpMat>> matPtrVec;
    for (int i = 0; i < filenames.size(); ++i) {
        planc::H5SpMat h5spm(filenames[i], rowindPath[i], colptrPath[i], valuePath[i], nrow[i], ncol[i]);
        std::unique_ptr<planc::H5SpMat> ptr = std::make_unique<planc::H5SpMat>(h5spm);
        matPtrVec.push_back(std::move(ptr));
    }
    planc::BPPINMF<planc::H5SpMat> solver(matPtrVec, k, lambda);

    if (Winit.isNotNull()) {
        arma::mat W = Rcpp::as<arma::mat>(Winit);
        solver.initW(W);
    } else {
        solver.initW();
    }

    if (Vinit.isNotNull()) {
        std::vector<arma::mat> VinitList = Rcpp::as<std::vector<arma::mat>>(Vinit);
        solver.initV(VinitList);
    } else {
        solver.initV();
    }

    if (Hinit.isNotNull()) {
        std::vector<arma::mat> HinitList = Rcpp::as<std::vector<arma::mat>>(Hinit);
        solver.initH(HinitList);
    } else {
        solver.initH();
    }

    solver.optimizeALS(niter, verbose);
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
Rcpp::List onlineINMF_S1_mem(std::vector<T> objectList, arma::uword k,
    double lambda, arma::uword maxEpoch = 5, arma::uword minibatchSize = 5000,
    arma::uword maxHALSIter = 1, bool verbose = true) {
    std::vector<std::unique_ptr<T>> matPtrVec = initMemMatPtr<T>(objectList);
    planc::ONLINEINMF<T, T> solver(matPtrVec, k, lambda);
    solver.runOnlineINMF(minibatchSize, maxEpoch, maxHALSIter, verbose);

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
Rcpp::List onlineINMF_S1(Rcpp::List objectList, arma::uword k,
    double lambda, arma::uword maxEpoch = 5, arma::uword minibatchSize = 5000,
    arma::uword maxHALSIter = 1, bool verbose = true) {
    if (Rf_isS4(objectList[0])) {
        return onlineINMF_S1_mem<arma::sp_mat>(Rcpp::as<std::vector<arma::sp_mat>>(objectList),
            k, lambda, maxEpoch, minibatchSize, maxHALSIter, verbose);
    } else {
        return onlineINMF_S1_mem<arma::mat>(Rcpp::as<std::vector<arma::mat>>(objectList),
            k, lambda, maxEpoch, minibatchSize, maxHALSIter, verbose);
    }
    return Rcpp::List::create();
}

// [[Rcpp::export(.onlineINMF_S1_h5dense)]]
Rcpp::List onlineINMF_S1_h5dense(std::vector<std::string> filenames,
    std::vector<std::string> dataPaths, arma::uword k,
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
    solver.runOnlineINMF(minibatchSize, maxEpoch, maxHALSIter, verbose);

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
    arma::uword k, double lambda,
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
    solver.runOnlineINMF(minibatchSize, maxEpoch, maxHALSIter, verbose);

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
    arma::uword k, double lambda, bool project = false, arma::uword maxEpoch = 5,
    arma::uword minibatchSize = 5000, arma::uword maxHALSIter = 1, bool verbose = true) {
    std::vector<std::unique_ptr<T>> matPtrVec = initMemMatPtr<T>(objectList);
    std::vector<std::unique_ptr<T>> matPtrVecNew = initMemMatPtr<T>(objectListNew);
    planc::ONLINEINMF<T, T> solver(matPtrVec, k, lambda);
    solver.initV(Vinit, false);
    solver.initW(Winit, false);
    solver.initA(Ainit);
    solver.initB(Binit);
    solver.runOnlineINMF(matPtrVecNew, project, minibatchSize, maxEpoch, maxHALSIter, verbose);

    if (!project) {
        // Scenario 2
        int nDatasets = objectList.size() + objectListNew.size();
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
    std::vector<arma::mat> Vinit, arma::mat Winit,
    std::vector<arma::mat> Ainit, std::vector<arma::mat> Binit,
    Rcpp::List objectListNew,
    arma::uword k, double lambda, bool project = false, arma::uword maxEpoch = 5,
    arma::uword minibatchSize = 5000, arma::uword maxHALSIter = 1, bool verbose = true) {
    if (Rf_isS4(objectList[0])) {
        return onlineINMF_S23_mem<arma::sp_mat>(
            Rcpp::as<std::vector<arma::sp_mat>>(objectList),
            Vinit, Winit, Ainit, Binit,
            Rcpp::as<std::vector<arma::sp_mat>>(objectListNew),
            k, lambda, project, maxEpoch, minibatchSize, maxHALSIter, verbose);
    } else {
        return onlineINMF_S23_mem<arma::mat>(
            Rcpp::as<std::vector<arma::mat>>(objectList),
            Vinit, Winit, Ainit, Binit,
            Rcpp::as<std::vector<arma::mat>>(objectListNew),
            k, lambda, project, maxEpoch, minibatchSize, maxHALSIter, verbose);
    }
    return Rcpp::List::create();
}

// [[Rcpp::export(.onlineINMF_S23_h5dense)]]
Rcpp::List onlineINMF_S23_h5dense(
    std::vector<std::string> filenames, std::vector<std::string> dataPaths,
    std::vector<std::string> filenamesNew, std::vector<std::string> dataPathsNew,
    std::vector<arma::mat> Vinit, arma::mat Winit,
    std::vector<arma::mat> Ainit, std::vector<arma::mat> Binit,
    arma::uword k, double lambda, bool project = false, arma::uword maxEpoch = 5,
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
    solver.runOnlineINMF(matPtrVecNew, project, minibatchSize, maxEpoch, maxHALSIter, verbose);
    if (!project) {
        // Scenario 2
        int nDatasets = filenames.size() + filenamesNew.size();
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
    std::vector<arma::mat> Vinit, arma::mat Winit,
    std::vector<arma::mat> Ainit, std::vector<arma::mat> Binit,
    arma::uword k, double lambda, bool project = false, arma::uword maxEpoch = 5,
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
    solver.runOnlineINMF(matPtrVecNew, project, minibatchSize, maxEpoch, maxHALSIter, verbose);
    if (!project) {
        // Scenario 2
        int nDatasets = filenames.size() + filenamesNew.size();
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
                    arma::uword k, arma::vec lambda,
                    arma::uword niter, bool verbose)
{
    std::vector<std::unique_ptr<T>> matPtrVec;
    std::vector<std::unique_ptr<T>> unsharedPtrVec;
    matPtrVec = initMemMatPtr<T>(objectList);
    unsharedPtrVec = initMemMatPtr<T>(unsharedList);
    planc::UINMF<T> solver(matPtrVec, unsharedPtrVec, k, lambda);
    solver.optimizeUANLS(niter, verbose);

    Rcpp::List HList(objectList.size());
    Rcpp::List UList(objectList.size());
    Rcpp::List VList(objectList.size());
    for (arma::uword i = 0; i < objectList.size(); ++i)
    {
        HList[i] = solver.getHi(i);
        VList[i] = solver.getVi(i);
        UList[i] = solver.getUi(i);
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
Rcpp::List uinmf_rcpp(Rcpp::List objectList, Rcpp::List unsharedList,
                 arma::uword k, arma::vec lambda,
                 arma::uword niter, bool verbose) {
    if (Rf_isS4(objectList[0])) {
        return uinmf_mem<arma::sp_mat>(Rcpp::as<std::vector<arma::sp_mat>>(objectList),
                                Rcpp::as<std::vector<arma::sp_mat>>(unsharedList),
                                k, lambda, niter, verbose);
    } else {
        return uinmf_mem<arma::mat>(Rcpp::as<std::vector<arma::mat>>(objectList),
                                Rcpp::as<std::vector<arma::mat>>(unsharedList),
                                k, lambda, niter, verbose);
    }
    return Rcpp::List::create();
}

// [[Rcpp::export(.uinmf_h5dense)]]
Rcpp::List uinmf_h5dense(std::vector<std::string> filenames,
                         std::vector<std::string> dataPaths,
                         std::vector<std::string> unsharedFilenames,
                         std::vector<std::string> unsharedDataPaths,
                         arma::uword k, arma::vec lambda,
                         arma::uword niter, bool verbose) {
    std::vector<std::unique_ptr<planc::H5Mat>> matPtrVec;
    std::vector<std::unique_ptr<planc::H5Mat>> unsharedPtrVec;
    for (int i = 0; i < filenames.size(); ++i) {
        planc::H5Mat E(filenames[i], dataPaths[i]);
        std::unique_ptr<planc::H5Mat> ptr = std::make_unique<planc::H5Mat>(E);
        matPtrVec.push_back(std::move(ptr));

        planc::H5Mat E_unshared(unsharedFilenames[i], unsharedDataPaths[i]);
        std::unique_ptr<planc::H5Mat> ptr_unshared = std::make_unique<planc::H5Mat>(E_unshared);
        unsharedPtrVec.push_back(std::move(ptr_unshared));
    }
    planc::UINMF<planc::H5Mat> solver(matPtrVec, unsharedPtrVec, k, lambda);

    solver.optimizeUANLS(niter, verbose);

    Rcpp::List HList(matPtrVec.size());
    Rcpp::List UList(matPtrVec.size());
    Rcpp::List VList(matPtrVec.size());
    for (arma::uword i = 0; i < matPtrVec.size(); ++i)
    {
        HList[i] = solver.getHi(i);
        VList[i] = solver.getVi(i);
        UList[i] = solver.getUi(i);
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
                          arma::uword k, arma::vec lambda,
                          arma::uword niter, bool verbose) {
    std::vector<std::unique_ptr<planc::H5SpMat>> matPtrVec;
    std::vector<std::unique_ptr<planc::H5SpMat>> unsharedPtrVec;
    for (int i = 0; i < filenames.size(); ++i) {
        planc::H5SpMat E(filenames[i], rowindPaths[i], colptrPaths[i], valuePaths[i], nrows[i], ncols[i]);
        std::unique_ptr<planc::H5SpMat> ptr = std::make_unique<planc::H5SpMat>(E);
        matPtrVec.push_back(std::move(ptr));

        planc::H5SpMat E_unshared(unsharedFilenames[i], unsharedRowindPaths[i], unsharedColptrPaths[i], unsharedValuePaths[i], unsharedNrows[i], unsharedNcols[i]);
        std::unique_ptr<planc::H5SpMat> ptr_unshared = std::make_unique<planc::H5SpMat>(E_unshared);
        unsharedPtrVec.push_back(std::move(ptr_unshared));
    }
    planc::UINMF<planc::H5SpMat> solver(matPtrVec, unsharedPtrVec, k, lambda);
    solver.optimizeUANLS(niter, verbose);
    Rcpp::List HList(matPtrVec.size());
    Rcpp::List UList(matPtrVec.size());
    Rcpp::List VList(matPtrVec.size());
    for (arma::uword i = 0; i < matPtrVec.size(); ++i)
    {
        HList[i] = solver.getHi(i);
        VList[i] = solver.getVi(i);
        UList[i] = solver.getUi(i);
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
