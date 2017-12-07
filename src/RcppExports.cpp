// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mvtSampler
NumericVector mvtSampler(NumericVector y, NumericVector mu, IntegerVector selected, NumericMatrix threshold, NumericMatrix precision, int nsamp, int burnin, int trim, bool verbose);
RcppExport SEXP _variationalReg_mvtSampler(SEXP ySEXP, SEXP muSEXP, SEXP selectedSEXP, SEXP thresholdSEXP, SEXP precisionSEXP, SEXP nsampSEXP, SEXP burninSEXP, SEXP trimSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type selected(selectedSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< int >::type nsamp(nsampSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type trim(trimSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(mvtSampler(y, mu, selected, threshold, precision, nsamp, burnin, trim, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_variationalReg_mvtSampler", (DL_FUNC) &_variationalReg_mvtSampler, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_variationalReg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}