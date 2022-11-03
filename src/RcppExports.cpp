// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// GWAS
SEXP GWAS(Eigen::MatrixXd Y, Eigen::MatrixXd GEN, Eigen::MatrixXd M, int maxit, double logtol, int cores);
RcppExport SEXP _mas_GWAS(SEXP YSEXP, SEXP GENSEXP, SEXP MSEXP, SEXP maxitSEXP, SEXP logtolSEXP, SEXP coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type GEN(GENSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type logtol(logtolSEXP);
    Rcpp::traits::input_parameter< int >::type cores(coresSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS(Y, GEN, M, maxit, logtol, cores));
    return rcpp_result_gen;
END_RCPP
}
// MLM
SEXP MLM(Eigen::MatrixXd Y, Eigen::MatrixXd X, Eigen::MatrixXd Z, int maxit, double logtol, int cores);
RcppExport SEXP _mas_MLM(SEXP YSEXP, SEXP XSEXP, SEXP ZSEXP, SEXP maxitSEXP, SEXP logtolSEXP, SEXP coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type logtol(logtolSEXP);
    Rcpp::traits::input_parameter< int >::type cores(coresSEXP);
    rcpp_result_gen = Rcpp::wrap(MLM(Y, X, Z, maxit, logtol, cores));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mas_GWAS", (DL_FUNC) &_mas_GWAS, 6},
    {"_mas_MLM", (DL_FUNC) &_mas_MLM, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_mas(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
