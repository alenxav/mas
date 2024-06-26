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
SEXP GWAS(Eigen::MatrixXf Y, Eigen::MatrixXf GEN, Eigen::MatrixXf M, int NumPCs, int maxit, float logtol, int cores, bool verb);
RcppExport SEXP _mas_GWAS(SEXP YSEXP, SEXP GENSEXP, SEXP MSEXP, SEXP NumPCsSEXP, SEXP maxitSEXP, SEXP logtolSEXP, SEXP coresSEXP, SEXP verbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXf >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXf >::type GEN(GENSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXf >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type NumPCs(NumPCsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< float >::type logtol(logtolSEXP);
    Rcpp::traits::input_parameter< int >::type cores(coresSEXP);
    Rcpp::traits::input_parameter< bool >::type verb(verbSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS(Y, GEN, M, NumPCs, maxit, logtol, cores, verb));
    return rcpp_result_gen;
END_RCPP
}
// MLM
SEXP MLM(Eigen::MatrixXf Y, Eigen::MatrixXf X, Eigen::MatrixXf Z, int maxit, float logtol, int cores);
RcppExport SEXP _mas_MLM(SEXP YSEXP, SEXP XSEXP, SEXP ZSEXP, SEXP maxitSEXP, SEXP logtolSEXP, SEXP coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXf >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXf >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXf >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< float >::type logtol(logtolSEXP);
    Rcpp::traits::input_parameter< int >::type cores(coresSEXP);
    rcpp_result_gen = Rcpp::wrap(MLM(Y, X, Z, maxit, logtol, cores));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mas_GWAS", (DL_FUNC) &_mas_GWAS, 8},
    {"_mas_MLM", (DL_FUNC) &_mas_MLM, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_mas(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
