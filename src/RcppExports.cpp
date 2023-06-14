// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// aoadmmnmf
Rcpp::List aoadmmnmf(const arma::sp_mat& x, const int& k, const int& niter, const Rcpp::Nullable<Rcpp::NumericMatrix>& W_init, const Rcpp::Nullable<Rcpp::NumericMatrix>& H_init);
RcppExport SEXP _RcppPlanc_aoadmmnmf(SEXP xSEXP, SEXP kSEXP, SEXP niterSEXP, SEXP W_initSEXP, SEXP H_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int& >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type W_init(W_initSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type H_init(H_initSEXP);
    rcpp_result_gen = Rcpp::wrap(aoadmmnmf(x, k, niter, W_init, H_init));
    return rcpp_result_gen;
END_RCPP
}
// gnsymnmf
Rcpp::List gnsymnmf(const arma::sp_mat& x, const int& k, const int& niter, const Rcpp::Nullable<Rcpp::NumericMatrix>& W_init, const Rcpp::Nullable<Rcpp::NumericMatrix>& H_init);
RcppExport SEXP _RcppPlanc_gnsymnmf(SEXP xSEXP, SEXP kSEXP, SEXP niterSEXP, SEXP W_initSEXP, SEXP H_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int& >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type W_init(W_initSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type H_init(H_initSEXP);
    rcpp_result_gen = Rcpp::wrap(gnsymnmf(x, k, niter, W_init, H_init));
    return rcpp_result_gen;
END_RCPP
}
// halsnmf
Rcpp::List halsnmf(const arma::sp_mat& x, const int& k, const int& niter, const Rcpp::Nullable<Rcpp::NumericMatrix>& W_init, const Rcpp::Nullable<Rcpp::NumericMatrix>& H_init);
RcppExport SEXP _RcppPlanc_halsnmf(SEXP xSEXP, SEXP kSEXP, SEXP niterSEXP, SEXP W_initSEXP, SEXP H_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int& >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type W_init(W_initSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type H_init(H_initSEXP);
    rcpp_result_gen = Rcpp::wrap(halsnmf(x, k, niter, W_init, H_init));
    return rcpp_result_gen;
END_RCPP
}
// munmf
Rcpp::List munmf(const arma::sp_mat& x, const int& k, const int& niter, const Rcpp::Nullable<Rcpp::NumericMatrix>& W_init, const Rcpp::Nullable<Rcpp::NumericMatrix>& H_init);
RcppExport SEXP _RcppPlanc_munmf(SEXP xSEXP, SEXP kSEXP, SEXP niterSEXP, SEXP W_initSEXP, SEXP H_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int& >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type W_init(W_initSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type H_init(H_initSEXP);
    rcpp_result_gen = Rcpp::wrap(munmf(x, k, niter, W_init, H_init));
    return rcpp_result_gen;
END_RCPP
}
// bppnmf
Rcpp::List bppnmf(const arma::sp_mat& x, const int& k, const int& niter, const Rcpp::Nullable<Rcpp::NumericMatrix>& W_init, const Rcpp::Nullable<Rcpp::NumericMatrix>& H_init);
RcppExport SEXP _RcppPlanc_bppnmf(SEXP xSEXP, SEXP kSEXP, SEXP niterSEXP, SEXP W_initSEXP, SEXP H_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int& >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type W_init(W_initSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::NumericMatrix>& >::type H_init(H_initSEXP);
    rcpp_result_gen = Rcpp::wrap(bppnmf(x, k, niter, W_init, H_init));
    return rcpp_result_gen;
END_RCPP
}
// bppnnls
arma::mat bppnnls(const arma::mat& C, const arma::sp_mat& B);
RcppExport SEXP _RcppPlanc_bppnnls(SEXP CSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(bppnnls(C, B));
    return rcpp_result_gen;
END_RCPP
}
// bppinmf
Rcpp::List bppinmf(std::vector<Rcpp::NumericMatrix> objectList);
RcppExport SEXP _RcppPlanc_bppinmf(SEXP objectListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<Rcpp::NumericMatrix> >::type objectList(objectListSEXP);
    rcpp_result_gen = Rcpp::wrap(bppinmf(objectList));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RcppPlanc_aoadmmnmf", (DL_FUNC) &_RcppPlanc_aoadmmnmf, 5},
    {"_RcppPlanc_gnsymnmf", (DL_FUNC) &_RcppPlanc_gnsymnmf, 5},
    {"_RcppPlanc_halsnmf", (DL_FUNC) &_RcppPlanc_halsnmf, 5},
    {"_RcppPlanc_munmf", (DL_FUNC) &_RcppPlanc_munmf, 5},
    {"_RcppPlanc_bppnmf", (DL_FUNC) &_RcppPlanc_bppnmf, 5},
    {"_RcppPlanc_bppnnls", (DL_FUNC) &_RcppPlanc_bppnnls, 2},
    {"_RcppPlanc_bppinmf", (DL_FUNC) &_RcppPlanc_bppinmf, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_RcppPlanc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
