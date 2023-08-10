// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

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
        // BPPNNLS<arma::mat, arma::vec> solveProblem(C, (arma::mat)B.cols(spanStart, spanEnd));
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

template <typename T>
std::vector<std::unique_ptr<T>> bppinmfinit(std::vector<T> objectList)
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
                   arma::uword maxIter, double thresh, bool verbose)
{
    std::vector<std::unique_ptr<T>> matPtrVec;
    matPtrVec = bppinmfinit<T>(objectList);
    planc::BPPINMF<T> solver(matPtrVec, k, lambda);
    solver.initH();
    solver.initV();
    solver.initW();
    solver.optimizeALS(maxIter, thresh, verbose);
    std::vector<arma::mat> HList;
    std::vector<arma::mat> VList;
    for (arma::uword i = 0; i < objectList.size(); ++i)
    {
        HList.push_back(solver.getHi(i));
        VList.push_back(solver.getVi(i));
    }
    return Rcpp::List::create(
        Rcpp::Named("H") = HList,
        Rcpp::Named("V") = VList,
        Rcpp::Named("W") = solver.getW(),
        Rcpp::Named("objErr") = solver.objErr());
}
template <typename T>
Rcpp::List runINMF(std::vector<T> objectList, arma::uword k, double lambda,
                   arma::uword maxIter, double thresh, bool verbose,
                   std::vector<arma::mat> HinitList, std::vector<arma::mat> VinitList, arma::mat Winit)
{
    std::vector<std::unique_ptr<T>> matPtrVec;
    matPtrVec = bppinmfinit<T>(objectList);
    planc::BPPINMF<T> solver(matPtrVec, k, lambda);
    solver.initW(Winit);
    solver.initH(HinitList);
    solver.initV(VinitList);
    solver.optimizeALS(maxIter, thresh, verbose);
    std::vector<arma::mat> HList;
    std::vector<arma::mat> VList;
    for (arma::uword i = 0; i < objectList.size(); ++i)
    {
        HList.push_back(solver.getHi(i));
        VList.push_back(solver.getVi(i));
    }
    return Rcpp::List::create(
        Rcpp::Named("H") = HList,
        Rcpp::Named("V") = VList,
        Rcpp::Named("W") = solver.getW(),
        Rcpp::Named("objErr") = solver.objErr());
}

Rcpp::List bppinmf_dense(std::vector<arma::mat> objectList, arma::uword k,
                         double lambda, arma::uword maxIter, double thresh, bool verbose = true,
                         Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
                         Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
                         Rcpp::Nullable<arma::mat> Winit = R_NilValue)
{
    if (Hinit.isNotNull() && Vinit.isNotNull() && Winit.isNotNull())
    {
        return runINMF<arma::mat>(objectList, k, lambda,
                                     maxIter, thresh, verbose,
                                     Rcpp::as<std::vector<arma::mat>>(Hinit),
                                     Rcpp::as<std::vector<arma::mat>>(Vinit),
                                     Rcpp::as<arma::mat>(Winit));
    }
    else
    {
        return runINMF<arma::mat>(objectList, k, lambda,
                                     maxIter, thresh, verbose);
    }
}

Rcpp::List bppinmf_sparse(std::vector<arma::sp_mat> objectList, arma::uword k, double lambda,
    arma::uword maxIter, double thresh, bool verbose = true,
    Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
    Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
    Rcpp::Nullable<arma::mat> Winit  = R_NilValue) {
    if (Hinit.isNotNull() && Vinit.isNotNull() && Winit.isNotNull()) {
        return runINMF<arma::sp_mat>(objectList, k, lambda,
                                     maxIter, thresh, verbose,
                                     Rcpp::as<std::vector<arma::mat>>(Hinit),
                                     Rcpp::as<std::vector<arma::mat>>(Vinit),
                                     Rcpp::as<arma::mat>(Winit));
    }
    else {
        return runINMF<arma::sp_mat>(objectList, k, lambda,
                maxIter, thresh, verbose);
    }
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
Rcpp::List bppinmf(Rcpp::List objectList, const arma::uword k,
                   const double lambda, const arma::uword maxIter,
                   const double thresh, const bool verbose = true,
                   Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
                   Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
                   Rcpp::Nullable<arma::mat> Winit = R_NilValue) {
                    if (Rf_isS4(objectList[0])) {
                      // warning: non-void function does not return a value in all control paths
                      if (Rf_inherits(objectList[0], "dgCMatrix")) {
                        bppinmf_sparse(Rcpp::as<std::vector<arma::sp_mat>>(objectList), k, lambda,
                                       maxIter, thresh, verbose,
                                       Hinit,Vinit,Winit);
                      }
                      // warning: non-void function does not return a value in all control paths
                    } else {
                      return bppinmf_dense(Rcpp::as<std::vector<arma::mat>>(objectList), k, lambda,
                                            maxIter, thresh, verbose,
                                            Hinit, Vinit, Winit);
                    }
                   }

//' @export
// [[Rcpp::export]]
Rcpp::List bppinmf_h5dense(std::vector<std::string> filenames, std::vector<std::string> dataPath,
    arma::uword k, double lambda, arma::uword maxIter, double thresh,
    bool verbose = true,
    Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
    Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
    Rcpp::Nullable<arma::mat> Winit  = R_NilValue) {
        if (dataPath.size() == 1) {
            for (int i = 0; i < filenames.size() - 1; ++i) {
                dataPath.push_back(dataPath[0]);
            }
        }

        std::vector<std::unique_ptr<planc::H5Mat>> matPtrVec;
        for (int i = 0; i < filenames.size(); ++i) {
            std::cout << filenames[i] << " " << dataPath[i] << std::endl;
            planc::H5Mat h5m(filenames[i], dataPath[i]);
            std::unique_ptr<planc::H5Mat> ptr = std::make_unique<planc::H5Mat>(h5m);
            matPtrVec.push_back(std::move(ptr));
        }
        planc::BPPINMF<planc::H5Mat> solver(matPtrVec, k, lambda);

        if (Hinit.isNotNull()) {
            std::vector<arma::mat> HinitList = Rcpp::as<std::vector<arma::mat>>(Hinit);
            solver.initH(HinitList);
        } else {
            solver.initH();
        }

        if (Vinit.isNotNull()) {
            std::vector<arma::mat> VinitList = Rcpp::as<std::vector<arma::mat>>(Vinit);
            solver.initV(VinitList);
        } else {
            solver.initV();
        }

        if (Winit.isNotNull()) {
            arma::mat W = Rcpp::as<arma::mat>(Winit);
            solver.initW(W);
        } else {
            solver.initW();
        }

        solver.optimizeALS(maxIter, thresh, verbose);
        Rcpp::List HList = Rcpp::List::create();
        Rcpp::List VList = Rcpp::List::create();
        for (arma::uword i = 0; i < filenames.size(); ++i) {
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
void bppinmf_h5sparse(
    std::vector<std::string> filenames,
    std::vector<std::string> rowindPath,
    std::vector<std::string> colptrPath,
    std::vector<std::string> valuePath,
    arma::uword nrow, arma::uvec ncol,
    arma::uword k, double lambda, arma::uword maxIter, double thresh,
    bool verbose = true,
    Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
    Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
    Rcpp::Nullable<arma::mat> Winit  = R_NilValue) {
        if (rowindPath.size() == 1) {
            for (int i = 0; i < filenames.size() - 1; ++i) {
                rowindPath.push_back(rowindPath[0]);
            }
        }
        if (colptrPath.size() == 1) {
            for (int i = 0; i < filenames.size() - 1; ++i) {
                colptrPath.push_back(colptrPath[0]);
            }
        }
        if (valuePath.size() == 1) {
            for (int i = 0; i < filenames.size() - 1; ++i) {
                valuePath.push_back(valuePath[0]);
            }
        }

        std::vector<std::unique_ptr<planc::H5SpMat>> matPtrVec;
        for (int i = 0; i < filenames.size(); ++i) {
            planc::H5SpMat h5spm(filenames[i], rowindPath[i], colptrPath[i], valuePath[i], nrow, ncol[i]);
            std::unique_ptr<planc::H5SpMat> ptr = std::make_unique<planc::H5SpMat>(h5spm);
            matPtrVec.push_back(std::move(ptr));
        }
        // planc::BPPINMF<planc::H5SpMat> solver(matPtrVec, k, lambda);

        // if (Hinit.isNotNull()) {
        //     std::vector<arma::mat> HinitList = Rcpp::as<std::vector<arma::mat>>(Hinit);
        //     solver.initH(HinitList);
        // } else {
        //     solver.initH();
        // }

        // if (Vinit.isNotNull()) {
        //     std::vector<arma::mat> VinitList = Rcpp::as<std::vector<arma::mat>>(Vinit);
        //     solver.initV(VinitList);
        // } else {
        //     solver.initV();
        // }

        // if (Winit.isNotNull()) {
        //     arma::mat W = Rcpp::as<arma::mat>(Winit);
        //     solver.initW(W);
        // } else {
        //     solver.initW();
        // }

        // solver.optimizeALS(maxIter, thresh, verbose);
        // Rcpp::List HList = Rcpp::List::create();
        // Rcpp::List VList = Rcpp::List::create();
        // for (arma::uword i = 0; i < filenames.size(); ++i) {
        //     HList.push_back(solver.getHi(i));
        //     VList.push_back(solver.getVi(i));
        // }
        // // return HList;
        // return Rcpp::List::create(
        //     Rcpp::Named("H") = HList,
        //     Rcpp::Named("V") = VList,
        //     Rcpp::Named("W") = solver.getW()
        // );
    }

//' @export
// [[Rcpp::export]]
Rcpp::List onlineINMF_dense(std::vector<arma::mat> objectList, arma::uword k,
    double lambda, arma::uword maxEpoch = 5, arma::uword minibatchSize = 5000,
    arma::uword maxHALSIter = 1) {
    std::vector<std::unique_ptr<arma::mat>> matPtrVec;
    for (arma::uword i = 0; i < objectList.size(); ++i)
    {
        arma::mat E = arma::mat(objectList[i].begin(), objectList[i].n_rows, objectList[i].n_cols, false, true);
        std::unique_ptr<arma::mat> ptr = std::make_unique<arma::mat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    planc::ONLINEINMF<arma::mat, arma::mat> solver(matPtrVec, k, lambda);
    solver.runOnlineINMF(minibatchSize, maxEpoch, maxHALSIter);
    Rcpp::List HList = Rcpp::List::create();
    Rcpp::List VList = Rcpp::List::create();
    Rcpp::List AList = Rcpp::List::create();
    Rcpp::List BList = Rcpp::List::create();
    for (arma::uword i = 0; i < objectList.size(); ++i) {
        HList.push_back(solver.getHi(i));
        VList.push_back(solver.getVi(i));
        AList.push_back(solver.getAi(i));
        BList.push_back(solver.getBi(i));
    }
    // return HList;
    return Rcpp::List::create(
        Rcpp::Named("H") = HList,
        Rcpp::Named("V") = VList,
        Rcpp::Named("W") = solver.getW(),
        Rcpp::Named("A") = AList,
        Rcpp::Named("B") = BList
    );
}


//' @export
// [[Rcpp::export]]
Rcpp::List onlineINMF_sparse(Rcpp::List objectList, arma::uword k,
    double lambda, arma::uword maxEpoch = 5, arma::uword minibatchSize = 5000,
    arma::uword maxHALSIter = 1) {
    std::vector<std::unique_ptr<arma::sp_mat>> matPtrVec;
    for (arma::uword i = 0; i < objectList.size(); ++i)
    {
        arma::sp_mat E = objectList[i];
        std::unique_ptr<arma::sp_mat> ptr = std::make_unique<arma::sp_mat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    planc::ONLINEINMF<arma::sp_mat, arma::sp_mat> solver(matPtrVec, k, lambda);
    solver.runOnlineINMF(minibatchSize, maxEpoch, maxHALSIter);

    Rcpp::List HList = Rcpp::List::create();
    Rcpp::List VList = Rcpp::List::create();
    Rcpp::List AList = Rcpp::List::create();
    Rcpp::List BList = Rcpp::List::create();
    for (arma::uword i = 0; i < objectList.size(); ++i) {
        HList.push_back(solver.getHi(i));
        VList.push_back(solver.getVi(i));
        AList.push_back(solver.getAi(i));
        BList.push_back(solver.getBi(i));
    }
    // return HList;
    return Rcpp::List::create(
        Rcpp::Named("H") = HList,
        Rcpp::Named("V") = VList,
        Rcpp::Named("W") = solver.getW(),
        Rcpp::Named("A") = AList,
        Rcpp::Named("B") = BList
    );
}

//' @export
// [[Rcpp::export]]
Rcpp::List onlineINMF_H5Dense(std::vector<std::string> filenames,
    std::vector<std::string> dataPaths, arma::uword k,
    double lambda, arma::uword maxEpoch = 5, arma::uword minibatchSize = 5000,
    arma::uword maxHALSIter = 1) {
    std::vector<std::unique_ptr<planc::H5Mat>> matPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i)
    {
        planc::H5Mat E(filenames[i], dataPaths[i]);
        std::unique_ptr<planc::H5Mat> ptr = std::make_unique<planc::H5Mat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    planc::ONLINEINMF<planc::H5Mat, arma::mat> solver(matPtrVec, k, lambda);
    solver.runOnlineINMF(minibatchSize, maxEpoch, maxHALSIter);

    Rcpp::List HList = Rcpp::List::create();
    Rcpp::List VList = Rcpp::List::create();
    Rcpp::List AList = Rcpp::List::create();
    Rcpp::List BList = Rcpp::List::create();
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        HList.push_back(solver.getHi(i));
        VList.push_back(solver.getVi(i));
        AList.push_back(solver.getAi(i));
        BList.push_back(solver.getBi(i));
    }
    // return HList;
    return Rcpp::List::create(
        Rcpp::Named("H") = HList,
        Rcpp::Named("V") = VList,
        Rcpp::Named("W") = solver.getW(),
        Rcpp::Named("A") = AList,
        Rcpp::Named("B") = BList
    );
}

//' @export
// [[Rcpp::export]]
Rcpp::List onlineINMF_H5Sparse(
    std::vector<std::string> filenames,
    std::vector<std::string> valuePaths,
    std::vector<std::string> rowindPaths,
    std::vector<std::string> colptrPaths,
    arma::uvec nrows, arma::uvec ncols,
    arma::uword k, double lambda,
    arma::uword maxEpoch = 5, arma::uword minibatchSize = 5000,
    arma::uword maxHALSIter = 1
) {
    std::vector<std::unique_ptr<planc::H5SpMat>> matPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i)
    {
        planc::H5SpMat E(filenames[i], rowindPaths[i], colptrPaths[i], valuePaths[i], nrows[i], ncols[i]);
        std::unique_ptr<planc::H5SpMat> ptr = std::make_unique<planc::H5SpMat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    planc::ONLINEINMF<planc::H5SpMat, arma::sp_mat> solver(matPtrVec, k, lambda);
    solver.runOnlineINMF(minibatchSize, maxEpoch, maxHALSIter);

    Rcpp::List HList = Rcpp::List::create();
    Rcpp::List VList = Rcpp::List::create();
    Rcpp::List AList = Rcpp::List::create();
    Rcpp::List BList = Rcpp::List::create();
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        HList.push_back(solver.getHi(i));
        VList.push_back(solver.getVi(i));
        AList.push_back(solver.getAi(i));
        BList.push_back(solver.getBi(i));
    }
    // return HList;
    return Rcpp::List::create(
        Rcpp::Named("H") = HList,
        Rcpp::Named("V") = VList,
        Rcpp::Named("W") = solver.getW(),
        Rcpp::Named("A") = AList,
        Rcpp::Named("B") = BList
    );
}

//' @export
// [[Rcpp::export]]
Rcpp::List onlineINMF_Xnew_sparse(
    std::vector<arma::sp_mat> objectList,
    std::vector<arma::mat>& Vinit, arma::mat& Winit,
    std::vector<arma::mat>& Ainit, std::vector<arma::mat>& Binit,
    std::vector<arma::sp_mat> objectListNew,
    arma::uword k, double lambda, bool project = false, arma::uword maxEpoch = 5,
    arma::uword minibatchSize = 5000, arma::uword maxHALSIter = 1, bool verbose = true) {
    std::vector<std::unique_ptr<arma::sp_mat>> matPtrVec;
    for (arma::uword i = 0; i < objectList.size(); ++i)
    {
        arma::sp_mat E = objectList[i];
        std::unique_ptr<arma::sp_mat> ptr = std::make_unique<arma::sp_mat>(E);
        matPtrVec.push_back(std::move(ptr));
    }

    std::vector<std::unique_ptr<arma::sp_mat>> matPtrVecNew;
    for (arma::uword i = 0; i < objectListNew.size(); ++i)
    {
        arma::sp_mat E_new = objectListNew[i];
        std::unique_ptr<arma::sp_mat> ptr = std::make_unique<arma::sp_mat>(E_new);
        matPtrVecNew.push_back(std::move(ptr));
    }
    planc::ONLINEINMF<arma::sp_mat, arma::sp_mat> solver(matPtrVec, k, lambda);
    solver.initV(Vinit, false);
    solver.initW(Winit, false);
    solver.initA(Ainit);
    solver.initB(Binit);
    solver.runOnlineINMF(matPtrVecNew, project, minibatchSize, maxEpoch, maxHALSIter, verbose);

    if (!project) {
        // Scenario 2
        Rcpp::List HList = Rcpp::List::create();
        Rcpp::List VList = Rcpp::List::create();
        Rcpp::List AList = Rcpp::List::create();
        Rcpp::List BList = Rcpp::List::create();
        for (arma::uword i = 0; i < objectList.size() + objectListNew.size(); ++i) {
            HList.push_back(solver.getHi(i));
            VList.push_back(solver.getVi(i));
            AList.push_back(solver.getAi(i));
            BList.push_back(solver.getBi(i));
        }
        return Rcpp::List::create(
            Rcpp::Named("H") = HList,
            Rcpp::Named("V") = VList,
            Rcpp::Named("W") = solver.getW(),
            Rcpp::Named("A") = AList,
            Rcpp::Named("B") = BList
        );
    } else {
        // Scenario 3
        Rcpp::List HList = Rcpp::List::create();
        for (arma::uword i = 0; i < objectListNew.size(); ++i) {
            HList.push_back(solver.getHi(i));
        }
        return Rcpp::List::create(
            Rcpp::Named("H") = HList
        );
    }

}
