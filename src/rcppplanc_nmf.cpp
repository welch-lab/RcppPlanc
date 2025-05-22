// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#ifndef ARMA_64BIT_WORD
#define ARMA_64BIT_WORD
#endif
#include <RcppArmadillo.h>
#include "config.h"
#include <progress.hpp>
#include <functional>
#include <variant>
#include <data.hpp>
#include <nmf_lib.hpp>

// [[Rcpp::plugins(openmp)]]
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

// normtype m_input_normalization;
// int m_initseed;
// arma::fvec m_regW;
// arma::fvec m_regH;

// T2 e.g. arma::sp_mat
template<typename T2, typename eT = typename T2::elem_type>
Rcpp::List runNMF(T2 x, arma::uword k, const std::string&algo, const arma::uword&niter, const int&nCores,
                  const Rcpp::Nullable<arma::mat> Winit,
                  const Rcpp::Nullable<arma::mat> Hinit) {
    typedef planc::nmfOutput<eT> (*nmfCallType)(const T2&, const arma::uword&, const arma::uword&, const std::string&,
                                                const int&, const arma::mat&, const arma::mat&);
    planc::nmfOutput<eT> libcall{};
    std::function nmfcall = static_cast<nmfCallType>(planc::nmflib<T2>::nmf);
    if (!Winit.isNotNull() and !Hinit.isNotNull()) libcall = planc::nmflib<T2>::nmf(x, k, niter, algo, nCores);
    else {
        using namespace std::placeholders;
        auto nmfcallBound = [nmfcall, capture0 = Rcpp::as<arma::mat>(Winit), capture1 = Rcpp::as<arma::mat>(Hinit)
                ](auto&&PH1, auto&&PH2, auto&&PH3, auto&&PH4, auto&&PH5) {
            return nmfcall(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2),
                           std::forward<decltype(PH3)>(PH3), std::forward<decltype(PH4)>(PH4),
                           std::forward<decltype(PH5)>(PH5), capture0, capture1);
        };
        libcall = nmfcallBound(x, k, niter, algo, nCores);
    }
    return Rcpp::List::create(
        Rcpp::Named("W") = libcall.outW,
        Rcpp::Named("H") = libcall.outH,
        Rcpp::Named("objErr") = libcall.objErr
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
Rcpp::List nmf(const SEXP&x, const arma::uword&k, const arma::uword&niter = 30,
               const std::string&algo = "anlsbpp",
               const int&nCores = 2,
               const Rcpp::Nullable<arma::mat>&Winit = R_NilValue,
               const Rcpp::Nullable<arma::mat>&Hinit = R_NilValue) {
    Rcpp::List outlist;
    try {
        if (Rf_isS4(x)) {
            // Assume using dgCMatrix
            outlist = runNMF<arma::sp_mat>(Rcpp::as<arma::sp_mat>(x), k, algo, niter, nCores, Winit, Hinit);
        }
        else {
            // Assume regular dense matrix
            outlist = runNMF<arma::mat>(Rcpp::as<arma::mat>(x), k, algo, niter, nCores, Winit, Hinit);
        }
    }
    catch (const std::exception&e) {
        throw Rcpp::exception(e.what());
    }
    return outlist;
}

// // T1 e.g. BPPNMF<arma::sp_mat>
// // T2 e.g. arma::sp_mat
template<class T2, typename eT = typename T2::elem_type>
Rcpp::List runSymNMF(const T2&x, const arma::uword&k, const int&nCores, const arma::uword&niter, const double&symm_reg,
                     const std::string&algo,
                     const Rcpp::Nullable<arma::mat> Hinit) {
    typedef planc::nmfOutput<eT> (*symnmfCallType)(const T2&, const arma::uword&, const arma::uword&, const double&,
                                                   const std::string&, const int&, const arma::mat&);
    planc::nmfOutput<eT> libcall{};
    std::function nmfcall = static_cast<symnmfCallType>(planc::nmflib<T2>::symNMF);
    if (!Hinit.isNotNull()) {
        arma::mat nullMat{};
        libcall = planc::nmflib<T2>::symNMF(x, k, niter, symm_reg, algo, nCores, nullMat);
    }
    else {
        using namespace std::placeholders;
        auto nmfcallBoundH = std::bind(nmfcall, _1, _2, _3, _4, _5, _6, Rcpp::as<arma::mat>(Hinit));
        libcall = nmfcallBoundH(x, k, niter, symm_reg, algo, nCores);
    }
    return Rcpp::List::create(
        Rcpp::Named("W") = libcall.outW,
        Rcpp::Named("H") = libcall.outH,
        Rcpp::Named("objErr") = libcall.objErr
    );
}

//
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
Rcpp::List symNMF(const SEXP&x, const arma::uword&k, const arma::uword&niter = 30,
                  const double&lambda = 0.0, const std::string&algo = "gnsym", const int&nCores = 2,
                  const Rcpp::Nullable<arma::mat>&Hinit = R_NilValue) {
    //   arma::mat out;
    Rcpp::List out;
    try {
        if (Rf_isS4(x)) {
            // Assume using dgCMatrix
            out = runSymNMF<arma::sp_mat>(
                Rcpp::as<arma::sp_mat>(x), k, nCores, niter, lambda, algo, Hinit
            );
        }
        else {
            // Assume using default matrix
            out = runSymNMF<arma::mat>(
                Rcpp::as<arma::mat>(x), k, nCores, niter, lambda, algo, Hinit
            );
        }
    }
    catch (const std::exception&e) {
        throw Rcpp::exception(e.what());
    }
    return out;
}



//
//
// // %%%%%%%%%%%%%%%%%% BPPINMF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//
template<typename T, typename eT = typename T::elem_type>
Rcpp::List runINMF(std::vector<T> objectList, arma::uword k, double lambda,
                   arma::uword niter, bool verbose, const int&ncores) {
    planc::inmfOutput<eT> libcall{};
    libcall = planc::nmflib<T>::bppinmf(objectList, k, lambda, niter, verbose, ncores);
    return Rcpp::List::create(
        Rcpp::Named("H") = Rcpp::wrap(libcall.outHList),
        Rcpp::Named("V") = Rcpp::wrap(libcall.outVList),
        Rcpp::Named("W") = libcall.outW,
        Rcpp::Named("objErr") = libcall.objErr);
}

arma::mat deRm(const Rcpp::Nullable<arma::mat> wrapped) {
    if (wrapped.isNotNull()) {
        return Rcpp::as<arma::mat>(wrapped);
    }
    return arma::mat{};
}

// this is gross and i don't like it
std::vector<arma::mat> deRm(const Rcpp::Nullable<std::vector<arma::mat>> wrapped, const size_t len) {
    if (wrapped.isNotNull()) {
        return Rcpp::as<std::vector<arma::mat>>(wrapped);
    }
    std::vector<arma::mat> null{};
    for (arma::uword i = 0; i < len; ++i) {
        null.push_back(arma::mat{});
    }
    return null;
}

// the correct way TODO this is to implement an Rcpp exporter and use the templated function below.
// [[Rcpp::export(.bppinmf_h5dense)]]
Rcpp::List bppinmf_h5(const std::vector<std::string>&filenames, const std::vector<std::string>&dataPath,
                      const arma::uword k, const int&nCores, const double lambda, const arma::uword niter,
                      const bool verbose = true,
                      const Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
                      const Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
                      const Rcpp::Nullable<arma::mat> Winit = R_NilValue) {
    std::vector<std::shared_ptr<planc::H5Mat>> matVec;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        matVec.push_back(std::make_shared<planc::H5Mat>(planc::H5Mat(filenames[i], dataPath[i])));
    }
    planc::inmfOutput<double> libcall{};
    if (Hinit.isNull() && Vinit.isNull() && Winit.isNull()) {
        libcall = planc::nmflib<planc::H5Mat>::bppinmf(matVec, k, lambda, niter, verbose, nCores);
    }
    else {
        const arma::mat WinitPass = deRm(Winit);
        const std::vector<arma::mat> VinitPass = deRm(Vinit, filenames.size());
        const std::vector<arma::mat> HinitPass = deRm(Hinit, filenames.size());
        libcall = planc::nmflib<planc::H5Mat>::bppinmf(matVec, k, lambda, niter, verbose, HinitPass, VinitPass,
                                                       WinitPass, nCores);
    }
    std::vector<Rcpp::NumericMatrix> HList;
    std::vector<Rcpp::NumericMatrix> VList;
    for (arma::uword i = 0; i < matVec.size(); ++i) {
        HList.push_back(Rcpp::wrap(libcall.outHList[i]));
        VList.push_back(Rcpp::wrap(libcall.outVList[i]));
    }
    return Rcpp::List::create(
        Rcpp::Named("H") = wrap(HList),
        Rcpp::Named("V") = wrap(VList),
        Rcpp::Named("W") = libcall.outW,
        Rcpp::Named("objErr") = libcall.objErr);
}

// the correct way TODO this is to implement an Rcpp exporter and use the templated function below.
// [[Rcpp::export(.bppinmf_h5sparse)]]
Rcpp::List bppinmf_h5sp(
    std::vector<std::string> filenames,
    std::vector<std::string> valuePath,
    std::vector<std::string> rowindPath,
    std::vector<std::string> colptrPath,
    arma::uvec nrow, arma::uvec ncol,
    arma::uword k, const int&nCores, double lambda, arma::uword niter,
    bool verbose = true,
    Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
    Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
    Rcpp::Nullable<arma::mat> Winit = R_NilValue) {
    std::vector<std::shared_ptr<planc::H5SpMat>> matVec;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        auto h5spm = std::make_shared<planc::H5SpMat>(filenames[i], rowindPath[i], colptrPath[i], valuePath[i], nrow[i],
                                                      ncol[i]);
        matVec.push_back(h5spm);
    }
    planc::inmfOutput<double> libcall{};
    if (Hinit.isNull() && Vinit.isNull() && Winit.isNull()) {
        libcall = planc::nmflib<planc::H5SpMat>::bppinmf(matVec, k, lambda, niter, verbose, nCores);
    }
    else {
        arma::mat WinitPass = deRm(Winit);
        std::vector<arma::mat> VinitPass = deRm(Vinit, filenames.size());
        std::vector<arma::mat> HinitPass = deRm(Hinit, filenames.size());
        libcall = planc::nmflib<planc::H5SpMat>::bppinmf(matVec, k, lambda, niter, verbose, HinitPass, VinitPass,
                                                         WinitPass, nCores);
    }
    std::vector<Rcpp::NumericMatrix> HList;
    std::vector<Rcpp::NumericMatrix> VList;
    for (arma::uword i = 0; i < matVec.size(); ++i) {
        HList.push_back(Rcpp::wrap(libcall.outHList[i]));
        VList.push_back(Rcpp::wrap(libcall.outVList[i]));
    }
    return Rcpp::List::create(
        Rcpp::Named("H") = wrap(HList),
        Rcpp::Named("V") = wrap(VList),
        Rcpp::Named("W") = libcall.outW,
        Rcpp::Named("objErr") = libcall.objErr);
}


template<typename T>
Rcpp::List bppinmf(const std::vector<std::shared_ptr<T>>&objectList, arma::uword k, double lambda,
                   arma::uword niter, bool verbose,
                   const Rcpp::Nullable<std::vector<arma::mat>> HinitList,
                   const Rcpp::Nullable<std::vector<arma::mat>> VinitList,
                   const Rcpp::Nullable<arma::mat> Winit,
                   const int&ncores) {
    planc::inmfOutput<double> libcall{};
    if (HinitList.isNull() && VinitList.isNull() && Winit.isNull()) {
        libcall = planc::nmflib<T>::bppinmf(objectList, k, lambda, niter, verbose, ncores);
    }
    else {
        libcall = planc::nmflib<T>::bppinmf(objectList, k, lambda, niter, verbose,
                                            Rcpp::as<std::vector<arma::mat>>(HinitList),
                                            Rcpp::as<std::vector<arma::mat>>(VinitList), Rcpp::as<arma::mat>(Winit),
                                            ncores);
    }
    std::vector<Rcpp::NumericMatrix> HList;
    std::vector<Rcpp::NumericMatrix> VList;
    for (arma::uword i = 0; i < objectList.size(); ++i) {
        HList.push_back(Rcpp::wrap(libcall.outHList[i]));
        VList.push_back(Rcpp::wrap(libcall.outVList[i]));
    }
    return Rcpp::List::create(
        Rcpp::Named("H") = wrap(HList),
        Rcpp::Named("V") = wrap(VList),
        Rcpp::Named("W") = libcall.outW,
        Rcpp::Named("objErr") = libcall.objErr);
}

// [[Rcpp::export(.bppinmf)]]
Rcpp::List bppinmf(const Rcpp::List&objectList, const arma::uword k, const int&nCores,
                   const double lambda = 5, const arma::uword niter = 30,
                   const bool verbose = true,
                   Rcpp::Nullable<std::vector<arma::mat>> Hinit = R_NilValue,
                   Rcpp::Nullable<std::vector<arma::mat>> Vinit = R_NilValue,
                   Rcpp::Nullable<arma::mat> Winit = R_NilValue) {
    std::variant<std::vector<std::shared_ptr<arma::mat>>, std::vector<std::shared_ptr<arma::sp_mat>>> unwrap;
    if (Rf_isS4(objectList[0])) {
        unwrap = planc::nmflib<arma::sp_mat>::initMemSharedPtr(Rcpp::as<std::vector<arma::sp_mat>>(objectList));
    }
    else {
        unwrap = planc::nmflib<arma::mat>::initMemSharedPtr(Rcpp::as<std::vector<arma::mat>>(objectList));
    }

    return std::visit([k, lambda, niter, verbose, Hinit, Vinit, Winit, nCores](auto&&arg) {
        return bppinmf(arg, k, lambda, niter, verbose, Hinit, Vinit, Winit, nCores);
    }, unwrap);
}


// %%%%%%%%%%%%%%%%%% online INMF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template<typename T>
Rcpp::List onlineINMF(std::vector<std::shared_ptr<T>> matPtrVec, arma::uword k, const int&nCores,
                      double lambda, arma::uword maxEpoch = 5, arma::uword minibatchSize = 5000,
                      arma::uword maxHALSIter = 1, arma::uword permuteChunkSize = 1000, bool verbose = true) {
    const planc::oinmfOutput<double> libcall = planc::nmflib<T, double>::oinmf(
        matPtrVec, k, nCores, lambda, maxEpoch, minibatchSize, maxHALSIter, permuteChunkSize, verbose);

    const arma::uword nDatasets = matPtrVec.size();
    std::vector<Rcpp::NumericMatrix> HList;
    std::vector<Rcpp::NumericMatrix> VList;
    std::vector<Rcpp::NumericMatrix> AList;
    std::vector<Rcpp::NumericMatrix> BList;
    for (arma::uword i = 0; i < nDatasets; ++i) {
        HList.push_back(Rcpp::wrap(libcall.outHList[i]));
        VList.push_back(Rcpp::wrap(libcall.outVList[i]));
        AList.push_back(Rcpp::wrap(libcall.outAList[i]));
        BList.push_back(Rcpp::wrap(libcall.outBList[i]));
    }
    return Rcpp::List::create(
        Rcpp::Named("H") = wrap(HList),
        Rcpp::Named("V") = wrap(VList),
        Rcpp::Named("W") = Rcpp::wrap(libcall.outW),
        Rcpp::Named("A") = wrap(AList),
        Rcpp::Named("B") = wrap(BList),
        Rcpp::Named("objErr") = Rcpp::wrap(libcall.objErr)
    );
}

// [[Rcpp::export(.onlineINMF_h5dense)]]
Rcpp::List onlineINMF_h5dense(const std::vector<std::string>&filenames,
                              const std::vector<std::string>&dataPaths,
                              const arma::uword k,
                              const int nCores,
                              const double lambda,
                              const arma::uword maxEpoch = 5,
                              const arma::uword minibatchSize = 5000,
                              const arma::uword maxHALSIter = 1,
                              const arma::uword permuteChunkSize = 1000,
                              const bool verbose = true) {
    std::vector<std::shared_ptr<planc::H5Mat>> matPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        planc::H5Mat E(filenames[i], dataPaths[i]);
        auto ptr = std::make_shared<planc::H5Mat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    return onlineINMF<planc::H5Mat>(matPtrVec, k, nCores, lambda, maxEpoch, minibatchSize, maxHALSIter, permuteChunkSize, verbose);
}

// [[Rcpp::export(.onlineINMF_h5sparse)]]
Rcpp::List onlineINMF_h5sparse(
    const std::vector<std::string>&filenames,
    const std::vector<std::string>&valuePaths,
    const std::vector<std::string>&rowindPaths,
    const std::vector<std::string>&colptrPaths,
    arma::uvec nrows, arma::uvec ncols,
    const arma::uword k, const int&nCores, const double lambda,
    const arma::uword maxEpoch = 5, const arma::uword minibatchSize = 5000,
    const arma::uword maxHALSIter = 1, const arma::uword permuteChunkSize = 1000,
    const bool verbose = true
) {
    std::vector<std::shared_ptr<planc::H5SpMat>> matPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        planc::H5SpMat E(filenames[i], rowindPaths[i], colptrPaths[i], valuePaths[i], nrows[i], ncols[i]);
        auto ptr = std::make_shared<planc::H5SpMat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    return onlineINMF<planc::H5SpMat>(matPtrVec, k, nCores, lambda, maxEpoch, minibatchSize, maxHALSIter, permuteChunkSize, verbose);
}

template<typename T>
Rcpp::List onlineINMF_initial(std::vector<std::shared_ptr<T>> matPtrVec,
                              std::vector<arma::mat> Hinit,
                              std::vector<arma::mat> Vinit, arma::mat Winit,
                              std::vector<arma::mat> Ainit, std::vector<arma::mat> Binit,
                              std::vector<std::shared_ptr<T>> matPtrVecNew,
                              arma::uword k, const int&nCores, double lambda, arma::uword maxEpoch = 5,
                              arma::uword minibatchSize = 5000, arma::uword maxHALSIter = 1,
                              arma::uword permuteChunkSize = 1000, bool verbose = true) {
    const planc::oinmfOutput<double> libcall = planc::nmflib<T, double>::oinmf(
        matPtrVec, Hinit, Vinit, Winit, Ainit, Binit, matPtrVecNew, k, nCores,
        lambda, maxEpoch, minibatchSize, maxHALSIter, permuteChunkSize, verbose
    );

    // Scenario 2
    const arma::uword nDatasets = matPtrVec.size() + matPtrVecNew.size();
    std::vector<Rcpp::NumericMatrix> HList;
    std::vector<Rcpp::NumericMatrix> VList;
    std::vector<Rcpp::NumericMatrix> AList;
    std::vector<Rcpp::NumericMatrix> BList;
    for (arma::uword i = 0; i < nDatasets; ++i) {
        HList.push_back(Rcpp::wrap(libcall.outHList[i]));
        VList.push_back(Rcpp::wrap(libcall.outVList[i]));
        AList.push_back(Rcpp::wrap(libcall.outAList[i]));
        BList.push_back(Rcpp::wrap(libcall.outBList[i]));
    }
    return Rcpp::List::create(
        Rcpp::Named("H") = wrap(HList),
        Rcpp::Named("V") = wrap(VList),
        Rcpp::Named("W") = Rcpp::wrap(libcall.outW),
        Rcpp::Named("A") = wrap(AList),
        Rcpp::Named("B") = wrap(BList),
        Rcpp::Named("objErr") = Rcpp::wrap(libcall.objErr)
    );
}

// [[Rcpp::export(.onlineINMF)]]
Rcpp::List onlineINMF(Rcpp::List objectList, arma::uword k, const int&nCores,
                      double lambda, arma::uword maxEpoch = 5, arma::uword minibatchSize = 5000,
                      arma::uword maxHALSIter = 1, arma::uword permuteChunkSize = 1000, bool verbose = true) {
    std::variant<std::vector<std::shared_ptr<arma::mat>>, std::vector<std::shared_ptr<arma::sp_mat>>> unwrap;
    if (Rf_isS4(objectList[0])) {
        unwrap = planc::nmflib<arma::sp_mat>::initMemSharedPtr(Rcpp::as<std::vector<arma::sp_mat>>(objectList));
    }
    else {
        unwrap = planc::nmflib<arma::mat>::initMemSharedPtr(Rcpp::as<std::vector<arma::mat>>(objectList));
    }
    return std::visit([k, nCores, lambda, maxEpoch, minibatchSize, maxHALSIter, permuteChunkSize, verbose](auto&&arg) {
        return onlineINMF(arg, k, nCores, lambda, maxEpoch, minibatchSize, maxHALSIter, permuteChunkSize, verbose);
    }, unwrap);
}

template<typename T>
Rcpp::List onlineINMF_project(
    std::vector<std::shared_ptr<T>> matPtrVec, arma::mat Winit,
    std::vector<std::shared_ptr<T>> matPtrVecNew,
    arma::uword k, const int&nCores, double lambda) {
    const std::vector<arma::mat> libcall = planc::nmflib<T>::oinmf_project(
        matPtrVec, Winit, matPtrVecNew, k, nCores, lambda);

    // Scenario 3
    std::vector<Rcpp::NumericMatrix> HList;
    for (arma::uword i = 0; i < matPtrVecNew.size(); ++i) {
        HList.push_back(Rcpp::wrap(libcall[i]));
    }

    return Rcpp::List::create(
        Rcpp::Named("H") = Rcpp::wrap(HList)
    );
}

// [[Rcpp::export(.onlineINMF_withInitial)]]
Rcpp::List onlineINMF_withInitial(
    Rcpp::List objectList,
    const std::vector<arma::mat>&Hinit,
    const std::vector<arma::mat>&Vinit, const arma::mat&Winit,
    const std::vector<arma::mat>&Ainit, const std::vector<arma::mat>&Binit,
    const Rcpp::List&objectListNew,
    const arma::uword k, const int&nCores, const double lambda,
    const arma::uword maxEpoch = 5,
    const arma::uword minibatchSize = 5000,
    const arma::uword maxHALSIter = 1,
    const arma::uword permuteChunkSize = 1000,
    const bool verbose = true) {
    //std::variant<std::vector<std::shared_ptr<arma::mat>>, std::vector<std::shared_ptr<arma::sp_mat>>> unwrap;
    //std::variant<std::vector<std::shared_ptr<arma::mat>>, std::vector<std::shared_ptr<arma::sp_mat>>> unwrap2;
    if (Rf_isS4(objectList[0])) {
        const std::vector<std::shared_ptr<arma::sp_mat>> matvec = planc::nmflib<arma::sp_mat>::initMemSharedPtr(
            Rcpp::as<std::vector<arma::sp_mat>>(objectList));
        const std::vector<std::shared_ptr<arma::sp_mat>> matvecnew = planc::nmflib<arma::sp_mat>::initMemSharedPtr(
            Rcpp::as<std::vector<arma::sp_mat>>(objectListNew));
        return onlineINMF_initial(matvec, Hinit, Vinit, Winit, Ainit, Binit, matvecnew, k, nCores, lambda, maxEpoch,
                                  minibatchSize, maxHALSIter, permuteChunkSize, verbose);
    }
    const std::vector<std::shared_ptr<arma::mat>> matvec = planc::nmflib<arma::mat>::initMemSharedPtr(
        Rcpp::as<std::vector<arma::mat>>(objectList));
    const std::vector<std::shared_ptr<arma::mat>> matvecnew = planc::nmflib<arma::mat>::initMemSharedPtr(
        Rcpp::as<std::vector<arma::mat>>(objectListNew));
    return onlineINMF_initial(matvec, Hinit, Vinit, Winit, Ainit, Binit, matvecnew, k, nCores, lambda, maxEpoch,
                              minibatchSize, maxHALSIter, permuteChunkSize, verbose);
    // return std::visit([Vinit, Winit, Ainit, Binit, k, nCores, lambda, maxEpoch, minibatchSize, maxHALSIter, verbose](auto&& arg, auto&& arg2) {
    //     return onlineINMF_initial(arg, Vinit, Winit, Ainit, Binit, arg2, k, nCores, lambda, maxEpoch, minibatchSize, maxHALSIter, verbose);
    // }, unwrap, unwrap2);
}

// [[Rcpp::export(.onlineINMF_project)]]
Rcpp::List onlineINMF_project(
    Rcpp::List objectList, const arma::mat&Winit,
    const Rcpp::List&objectListNew, const arma::uword k, const int&nCores, const double lambda) {
    //std::variant<std::vector<std::shared_ptr<arma::mat>>, std::vector<std::shared_ptr<arma::sp_mat>>> unwrap;
    //std::variant<std::vector<std::shared_ptr<arma::mat>>, std::vector<std::shared_ptr<arma::sp_mat>>> unwrap2;
    if (Rf_isS4(objectList[0])) {
        const std::vector<std::shared_ptr<arma::sp_mat>> matvec = planc::nmflib<arma::sp_mat>::initMemSharedPtr(
            Rcpp::as<std::vector<arma::sp_mat>>(objectList));
        const std::vector<std::shared_ptr<arma::sp_mat>> matvecnew = planc::nmflib<arma::sp_mat>::initMemSharedPtr(
            Rcpp::as<std::vector<arma::sp_mat>>(objectListNew));
        return onlineINMF_project(matvec, Winit, matvecnew, k, nCores, lambda);
    }
    const std::vector<std::shared_ptr<arma::mat>> matvec = planc::nmflib<arma::mat>::initMemSharedPtr(
        Rcpp::as<std::vector<arma::mat>>(objectList));
    const std::vector<std::shared_ptr<arma::mat>> matvecnew = planc::nmflib<arma::mat>::initMemSharedPtr(
        Rcpp::as<std::vector<arma::mat>>(objectListNew));
    return onlineINMF_project(matvec, Winit, matvecnew, k, nCores, lambda);
    // return std::visit([Winit, k, nCores, lambda](auto&& arg, auto&& arg2) {
    //     return onlineINMF_project(arg, Winit, arg2, k, nCores, lambda);
    // }, unwrap, unwrap2);
}

// [[Rcpp::export(.onlineINMF_h5dense_withInitial)]]
Rcpp::List onlineINMF_h5dense_withInitial(
    const std::vector<std::string>&filenames, const std::vector<std::string>&dataPaths,
    const std::vector<std::string>&filenamesNew, const std::vector<std::string>&dataPathsNew,
    const std::vector<arma::mat>&Hinit,
    const std::vector<arma::mat>&Vinit, const arma::mat&Winit,
    const std::vector<arma::mat>&Ainit, const std::vector<arma::mat>&Binit,
    const arma::uword k, const int&nCores, const double lambda, const arma::uword maxEpoch = 5,
    const arma::uword minibatchSize = 5000, const arma::uword maxHALSIter = 1,
    const arma::uword permuteChunkSize = 1000, const bool verbose = true) {
    std::vector<std::shared_ptr<planc::H5Mat>> matPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        planc::H5Mat E(filenames[i], dataPaths[i]);
        auto ptr = std::make_shared<planc::H5Mat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    std::vector<std::shared_ptr<planc::H5Mat>> matPtrVecNew;
    for (arma::uword i = 0; i < filenamesNew.size(); ++i) {
        planc::H5Mat E(filenamesNew[i], dataPathsNew[i]);
        auto ptr = std::make_shared<planc::H5Mat>(E);
        matPtrVecNew.push_back(std::move(ptr));
    } // Scenario 2
    return onlineINMF_initial<planc::H5Mat>(matPtrVec, Hinit, Vinit, Winit, Ainit, Binit, matPtrVecNew, k, nCores, lambda,
                                            maxEpoch, minibatchSize, maxHALSIter, permuteChunkSize, verbose);
}

// [[Rcpp::export(.onlineINMF_project_h5dense)]]
Rcpp::List onlineINMF_project_h5dense(
    const std::vector<std::string>&filenames, const std::vector<std::string>&dataPaths,
    const std::vector<std::string>&filenamesNew, const std::vector<std::string>&dataPathsNew,
    const arma::mat&Winit, const arma::uword k, const int&nCores, const double lambda) {
    std::vector<std::shared_ptr<planc::H5Mat>> matPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        auto E = planc::H5Mat(filenames[i], dataPaths[i]);
        auto ptr = std::make_shared<planc::H5Mat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    std::vector<std::shared_ptr<planc::H5Mat>> matPtrVecNew;
    for (arma::uword i = 0; i < filenamesNew.size(); ++i) {
        auto E = planc::H5Mat(filenamesNew[i], dataPathsNew[i]);
        auto ptr = std::make_shared<planc::H5Mat>(E);
        matPtrVecNew.push_back(std::move(ptr));
    }
    return onlineINMF_project<planc::H5Mat>(matPtrVec, Winit, matPtrVecNew,
                                            k, nCores, lambda);
}

// [[Rcpp::export(.onlineINMF_project_h5sparse)]]
Rcpp::List onlineINMF_project_h5sparse(
    std::vector<std::string> filenames, std::vector<std::string> valuePaths,
    std::vector<std::string> rowindPaths, std::vector<std::string> colptrPaths, 
    arma::uvec nrows, arma::uvec ncols, std::vector<std::string> filenamesNew, 
    std::vector<std::string> valuePathsNew, std::vector<std::string> rowindPathsNew, 
    std::vector<std::string> colptrPathsNew, arma::uvec nrowsNew, 
    arma::uvec ncolsNew, const arma::mat&Winit, arma::uword k, 
    const int&nCores, double lambda) {
    std::vector<std::shared_ptr<planc::H5SpMat>> matPtrVec;
    std::vector<std::shared_ptr<planc::H5SpMat>> matPtrVecNew;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        auto E = planc::H5SpMat(filenames[i], rowindPaths[i], colptrPaths[i], valuePaths[i], nrows[i], ncols[i]);
        auto ptr = std::make_shared<planc::H5SpMat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    for (arma::uword i = 0; i < filenamesNew.size(); ++i) {
        auto E = planc::H5SpMat(filenamesNew[i], rowindPathsNew[i], colptrPathsNew[i], valuePathsNew[i], nrowsNew[i],
                                ncolsNew[i]);
        auto ptr = std::make_shared<planc::H5SpMat>(E);
        matPtrVecNew.push_back(std::move(ptr));
    }
    return onlineINMF_project<planc::H5SpMat>(matPtrVec, Winit, matPtrVecNew,
                                              k, nCores, lambda);
}

// [[Rcpp::export(.onlineINMF_h5sparse_withInitial)]]
Rcpp::List onlineINMF_h5sparse_withInitial(
    std::vector<std::string> filenames, std::vector<std::string> valuePaths,
    std::vector<std::string> rowindPaths, std::vector<std::string> colptrPaths,
    arma::uvec nrows, arma::uvec ncols,
    std::vector<std::string> filenamesNew, std::vector<std::string> valuePathsNew,
    std::vector<std::string> rowindPathsNew, std::vector<std::string> colptrPathsNew,
    arma::uvec nrowsNew, arma::uvec ncolsNew,
    const std::vector<arma::mat>&Hinit,
    std::vector<arma::mat> Vinit, const arma::mat&Winit,
    std::vector<arma::mat> Ainit, std::vector<arma::mat> Binit,
    arma::uword k, const int&nCores, double lambda, arma::uword maxEpoch = 5,
    arma::uword minibatchSize = 5000, arma::uword maxHALSIter = 1,
    arma::uword permuteChunkSize = 1000, bool verbose = true
) {
    std::vector<std::shared_ptr<planc::H5SpMat>> matPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        planc::H5SpMat E(filenames[i], rowindPaths[i], colptrPaths[i], valuePaths[i], nrows[i], ncols[i]);
        auto ptr = std::make_shared<planc::H5SpMat>(E);
        matPtrVec.push_back(std::move(ptr));
    }
    std::vector<std::shared_ptr<planc::H5SpMat>> matPtrVecNew;
    for (arma::uword i = 0; i < filenamesNew.size(); ++i) {
        planc::H5SpMat E(filenamesNew[i], rowindPathsNew[i], colptrPathsNew[i], valuePathsNew[i], nrowsNew[i],
                         ncolsNew[i]);
        auto ptr = std::make_shared<planc::H5SpMat>(E);
        matPtrVecNew.push_back(std::move(ptr));
    }
    return onlineINMF_initial(matPtrVec, Hinit, Vinit, Winit, Ainit, Binit, matPtrVecNew, k, nCores, lambda, maxEpoch,
                              minibatchSize, maxHALSIter, permuteChunkSize, verbose);
}

// %%%%%%%%%%%%%% UINMF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template<typename T>
Rcpp::List uinmf_mem(std::vector<T> objectList,
                     std::vector<T> unsharedList,
                     std::vector<int> whichUnshared,
                     arma::uword k, const int&nCores, arma::vec lambda,
                     arma::uword niter, bool verbose) {
    const std::vector<std::shared_ptr<T>> matPtrVec = planc::nmflib<T>::initMemSharedPtr(
        objectList);
    const std::vector<std::shared_ptr<T>> unsharedPtrVec = planc::nmflib<T>::initMemSharedPtr(unsharedList);
    planc::uinmfOutput libcall = planc::nmflib<T>::uinmf(matPtrVec, unsharedPtrVec, whichUnshared, k, nCores, lambda,
                                                         niter, verbose);
    std::vector<Rcpp::NumericMatrix> HList;
    std::vector<Rcpp::NumericMatrix> VList;
    std::vector<Rcpp::NumericMatrix> UList;
    for (arma::uword i = 0; i < objectList.size(); ++i) {
        HList.push_back(Rcpp::wrap(libcall.outHList[i]));
        VList.push_back(Rcpp::wrap(libcall.outVList[i]));
    }
    for (arma::uword i = 0; i < unsharedList.size(); ++i) {
        UList.push_back(Rcpp::wrap(libcall.outUList[i]));
    }
    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("H") = wrap(HList),
        Rcpp::Named("V") = wrap(VList),
        Rcpp::Named("W") = Rcpp::wrap(libcall.outW),
        Rcpp::Named("U") = wrap(UList),
        Rcpp::Named("objErr") = Rcpp::wrap(libcall.objErr)
    );
    return output;
}

// [[Rcpp::export(.uinmf_rcpp)]]
Rcpp::List uinmf_rcpp(Rcpp::List objectList, const Rcpp::List&unsharedList,
                      const std::vector<int>&whichUnshared, arma::uword k, const int&nCores,
                      const arma::vec&lambda, arma::uword niter, bool verbose) {
    if (Rf_isS4(objectList[0])) {
        return uinmf_mem<arma::sp_mat>(
            Rcpp::as<std::vector<arma::sp_mat>>(objectList),
            Rcpp::as<std::vector<arma::sp_mat>>(unsharedList),
            whichUnshared, k, nCores, lambda, niter, verbose
        );
    }
    else {
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
                         arma::uword k, const int&nCores, const arma::vec&lambda,
                         arma::uword niter, bool verbose) {
    std::vector<std::shared_ptr<planc::H5Mat>> matPtrVec;
    std::vector<std::shared_ptr<planc::H5Mat>> unsharedPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        planc::H5Mat E(filenames[i], dataPaths[i]);
        std::shared_ptr<planc::H5Mat> ptr = std::make_shared<planc::H5Mat>(E);
        matPtrVec.push_back(std::move(ptr));

        planc::H5Mat E_unshared(unsharedFilenames[i], unsharedDataPaths[i]);
        std::shared_ptr<planc::H5Mat> ptr_unshared = std::make_shared<planc::H5Mat>(E_unshared);
        unsharedPtrVec.push_back(std::move(ptr_unshared));
    }
    planc::uinmfOutput libcall = planc::nmflib<planc::H5Mat>::uinmf(matPtrVec, unsharedPtrVec, whichUnshared, k, nCores,
                                                                    lambda, niter, verbose);
    std::vector<Rcpp::NumericMatrix> HList;
    std::vector<Rcpp::NumericMatrix> VList;
    std::vector<Rcpp::NumericMatrix> UList;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        HList.push_back(Rcpp::wrap(libcall.outHList[i]));
        VList.push_back(Rcpp::wrap(libcall.outVList[i]));
    }
    for (arma::uword i = 0; i < unsharedFilenames.size(); ++i) {
        UList.push_back(Rcpp::wrap(libcall.outUList[i]));
    }
    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("H") = wrap(HList),
        Rcpp::Named("V") = wrap(VList),
        Rcpp::Named("W") = Rcpp::wrap(libcall.outW),
        Rcpp::Named("U") = wrap(UList),
        Rcpp::Named("objErr") = Rcpp::wrap(libcall.objErr)
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
                          arma::uword k, const int&nCores, const arma::vec&lambda,
                          arma::uword niter, bool verbose) {
    std::vector<std::shared_ptr<planc::H5SpMat>> matPtrVec;
    std::vector<std::shared_ptr<planc::H5SpMat>> unsharedPtrVec;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        planc::H5SpMat E(filenames[i], rowindPaths[i], colptrPaths[i], valuePaths[i], nrows[i], ncols[i]);
        std::shared_ptr<planc::H5SpMat> ptr = std::make_shared<planc::H5SpMat>(E);
        matPtrVec.push_back(std::move(ptr));

        planc::H5SpMat E_unshared(unsharedFilenames[i], unsharedRowindPaths[i], unsharedColptrPaths[i],
                                  unsharedValuePaths[i], unsharedNrows[i], unsharedNcols[i]);
        std::shared_ptr<planc::H5SpMat> ptr_unshared = std::make_shared<planc::H5SpMat>(E_unshared);
        unsharedPtrVec.push_back(std::move(ptr_unshared));
    }
    planc::uinmfOutput libcall = planc::nmflib<planc::H5SpMat>::uinmf(matPtrVec, unsharedPtrVec, whichUnshared, k,
                                                                      nCores, lambda, niter, verbose);
    std::vector<Rcpp::NumericMatrix> HList;
    std::vector<Rcpp::NumericMatrix> VList;
    std::vector<Rcpp::NumericMatrix> UList;
    for (arma::uword i = 0; i < filenames.size(); ++i) {
        HList.push_back(Rcpp::wrap(libcall.outHList[i]));
        VList.push_back(Rcpp::wrap(libcall.outVList[i]));
    }
    for (arma::uword i = 0; i < unsharedFilenames.size(); ++i) {
        UList.push_back(Rcpp::wrap(libcall.outUList[i]));
    }
    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("H") = wrap(HList),
        Rcpp::Named("V") = wrap(VList),
        Rcpp::Named("W") = Rcpp::wrap(libcall.outW),
        Rcpp::Named("U") = wrap(UList),
        Rcpp::Named("objErr") = Rcpp::wrap(libcall.objErr)
    );
    return output;
}


// [[Rcpp::export(.openblaspthreadoff)]]
void openblas_pthread_off(Rcpp::XPtr<void *> libloc) {
    planc::openblas_pthread_off(libloc);
}

// [[Rcpp::export(.openblaspthreadon)]]
void openblas_pthread_on(Rcpp::XPtr<void *> libloc) {
    planc::openblas_pthread_on(libloc);
}
