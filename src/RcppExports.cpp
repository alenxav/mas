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
SEXP GWAS(Eigen::MatrixXf Y, Eigen::MatrixXf GEN, Eigen::MatrixXf M, int maxit, float logtol, int cores, bool verb);
RcppExport SEXP _mas_GWAS(SEXP YSEXP, SEXP GENSEXP, SEXP MSEXP, SEXP maxitSEXP, SEXP logtolSEXP, SEXP coresSEXP, SEXP verbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXf >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXf >::type GEN(GENSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXf >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< float >::type logtol(logtolSEXP);
    Rcpp::traits::input_parameter< int >::type cores(coresSEXP);
    Rcpp::traits::input_parameter< bool >::type verb(verbSEXP);
    rcpp_result_gen = Rcpp::wrap(GWAS(Y, GEN, M, maxit, logtol, cores, verb));
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
// MRR
SEXP MRR(Eigen::MatrixXf Y, Eigen::MatrixXf X, int maxit, float tol, int cores, bool TH, float NonLinearFactor, bool InnerGS, bool NoInversion, bool HCS, bool XFA, int NumXFA, float prior_R2, float gc_prior_df, float var_prior_df, float weight_prior_h2, float weight_prior_gc, float PenCor, float MinCor, float uncorH2below, float roundGCupFrom, float roundGCupTo, float roundGCdownFrom, float roundGCdownTo, float bucketGCfrom, float bucketGCto, float DeflateMax, float DeflateBy, bool OneVarB, bool OneVarE, bool verbose);
RcppExport SEXP _mas_MRR(SEXP YSEXP, SEXP XSEXP, SEXP maxitSEXP, SEXP tolSEXP, SEXP coresSEXP, SEXP THSEXP, SEXP NonLinearFactorSEXP, SEXP InnerGSSEXP, SEXP NoInversionSEXP, SEXP HCSSEXP, SEXP XFASEXP, SEXP NumXFASEXP, SEXP prior_R2SEXP, SEXP gc_prior_dfSEXP, SEXP var_prior_dfSEXP, SEXP weight_prior_h2SEXP, SEXP weight_prior_gcSEXP, SEXP PenCorSEXP, SEXP MinCorSEXP, SEXP uncorH2belowSEXP, SEXP roundGCupFromSEXP, SEXP roundGCupToSEXP, SEXP roundGCdownFromSEXP, SEXP roundGCdownToSEXP, SEXP bucketGCfromSEXP, SEXP bucketGCtoSEXP, SEXP DeflateMaxSEXP, SEXP DeflateBySEXP, SEXP OneVarBSEXP, SEXP OneVarESEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXf >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXf >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< float >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type cores(coresSEXP);
    Rcpp::traits::input_parameter< bool >::type TH(THSEXP);
    Rcpp::traits::input_parameter< float >::type NonLinearFactor(NonLinearFactorSEXP);
    Rcpp::traits::input_parameter< bool >::type InnerGS(InnerGSSEXP);
    Rcpp::traits::input_parameter< bool >::type NoInversion(NoInversionSEXP);
    Rcpp::traits::input_parameter< bool >::type HCS(HCSSEXP);
    Rcpp::traits::input_parameter< bool >::type XFA(XFASEXP);
    Rcpp::traits::input_parameter< int >::type NumXFA(NumXFASEXP);
    Rcpp::traits::input_parameter< float >::type prior_R2(prior_R2SEXP);
    Rcpp::traits::input_parameter< float >::type gc_prior_df(gc_prior_dfSEXP);
    Rcpp::traits::input_parameter< float >::type var_prior_df(var_prior_dfSEXP);
    Rcpp::traits::input_parameter< float >::type weight_prior_h2(weight_prior_h2SEXP);
    Rcpp::traits::input_parameter< float >::type weight_prior_gc(weight_prior_gcSEXP);
    Rcpp::traits::input_parameter< float >::type PenCor(PenCorSEXP);
    Rcpp::traits::input_parameter< float >::type MinCor(MinCorSEXP);
    Rcpp::traits::input_parameter< float >::type uncorH2below(uncorH2belowSEXP);
    Rcpp::traits::input_parameter< float >::type roundGCupFrom(roundGCupFromSEXP);
    Rcpp::traits::input_parameter< float >::type roundGCupTo(roundGCupToSEXP);
    Rcpp::traits::input_parameter< float >::type roundGCdownFrom(roundGCdownFromSEXP);
    Rcpp::traits::input_parameter< float >::type roundGCdownTo(roundGCdownToSEXP);
    Rcpp::traits::input_parameter< float >::type bucketGCfrom(bucketGCfromSEXP);
    Rcpp::traits::input_parameter< float >::type bucketGCto(bucketGCtoSEXP);
    Rcpp::traits::input_parameter< float >::type DeflateMax(DeflateMaxSEXP);
    Rcpp::traits::input_parameter< float >::type DeflateBy(DeflateBySEXP);
    Rcpp::traits::input_parameter< bool >::type OneVarB(OneVarBSEXP);
    Rcpp::traits::input_parameter< bool >::type OneVarE(OneVarESEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(MRR(Y, X, maxit, tol, cores, TH, NonLinearFactor, InnerGS, NoInversion, HCS, XFA, NumXFA, prior_R2, gc_prior_df, var_prior_df, weight_prior_h2, weight_prior_gc, PenCor, MinCor, uncorH2below, roundGCupFrom, roundGCupTo, roundGCdownFrom, roundGCdownTo, bucketGCfrom, bucketGCto, DeflateMax, DeflateBy, OneVarB, OneVarE, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mas_GWAS", (DL_FUNC) &_mas_GWAS, 7},
    {"_mas_MLM", (DL_FUNC) &_mas_MLM, 6},
    {"_mas_MRR", (DL_FUNC) &_mas_MRR, 31},
    {NULL, NULL, 0}
};

RcppExport void R_init_mas(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
