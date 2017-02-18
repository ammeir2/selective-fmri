// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sampleTruncNorm
NumericVector sampleTruncNorm(NumericVector sample, NumericVector lower, NumericVector upper, NumericVector mean, NumericMatrix precision, int cycles);
RcppExport SEXP selectivefmri_sampleTruncNorm(SEXP sampleSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP meanSEXP, SEXP precisionSEXP, SEXP cyclesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< int >::type cycles(cyclesSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleTruncNorm(sample, lower, upper, mean, precision, cycles));
    return rcpp_result_gen;
END_RCPP
}
