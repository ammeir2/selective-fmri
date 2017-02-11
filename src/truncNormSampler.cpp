#include <Rcpp.h>
using namespace Rcpp;

double sampleUnivTruncNorm(double mu, double sd, double lower, double upper) {
  double u = runif(1)[0] ;
  double phiB = R::pnorm5(upper, mu, sd, 1, 0) ;
  double phiA = R::pnorm5(lower, mu, sd, 1, 0) ;
  double quantile = u * phiB + phiA * (1 - u) ;
  double sample = R::qnorm(quantile, mu, sd, 1, 0) ;

  int tries = 0;

  // Sometimes the sampler doesn't work well, in those cases we use a rejection
  // sampler. They do something similar in the truncnorm package.
  while(sample < lower | sample > upper | isinf(std::abs(sample))) {
    sample = rnorm(1, mu, sd)[0] ;
    if(++tries > 1000) break ;
  }

  return sample ;
}

double computeConditionalMean(NumericVector mu, NumericVector samp,
                              NumericMatrix precision, int index) {
  double result = 0 ;

  for(int j = 0; j < mu.length() ; j ++) {
    if(j != index) {
      result += precision(index, j) * (samp[j] - mu[j]) ;
    }
  }

  result = result / precision(index, index) ;
  result = mu[index] - result ;
  return result ;
}

// [[Rcpp::export]]
NumericVector sampleTruncNorm(NumericVector sample,
                     NumericVector lower, NumericVector upper,
                     NumericVector mean, NumericMatrix precision,
                     int cycles) {
  double condMean ;
  double condVar ;
  int i, j ;
  NumericVector samp = clone(sample) ;

  for(i = 0; i < cycles ; i++) {
    for(j = 0; j < samp.length() ; j ++) {
      condMean = computeConditionalMean(mean, samp, precision, j) ;
      condVar = 1 / precision(j, j) ;
      samp[j] = sampleUnivTruncNorm(condMean, std::sqrt(condVar), lower[j], upper[j]) ;
    }
  }

  return samp ;
}


