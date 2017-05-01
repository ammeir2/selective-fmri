#include <Rcpp.h>
using namespace Rcpp;

double sampleExtreme(double mu, double sd, double lower, double upper) {
  double sign = 1 ;
  double threshold ;
  double proposal ;
  double alpha ;
  double phi ;

  if(isinf(lower)) {
    sign = -1 ;
    mu *= sign ;
    threshold = upper * sign ;
  } else {
    threshold = lower ;
  }

  // rescaling
  threshold = (threshold - mu) / sd ;
  alpha = (threshold + std::sqrt(std::pow(threshold, 2) + 4)) / 2 ;

  bool reject = true ;
  while(reject) {
    proposal = threshold + R::rexp(alpha) ;
    phi = std::exp(-std::pow(proposal - alpha, 2) / 2) ;
    if(runif(1)[0] < phi) {
      reject = false ;
    }
  }

  proposal = proposal * sd + mu ;
  return proposal * sign;
}

double sampleUnivTruncNorm(double mu, double sd, double lower, double upper) {
  double u = runif(1)[0] ;
  double phiB = R::pnorm5(upper, mu, sd, 1, 0) ;
  double phiA = R::pnorm5(lower, mu, sd, 1, 0) ;
  double quantile = u * phiB + phiA * (1 - u) ;
  double sample = R::qnorm(quantile, mu, sd, 1, 0) ;

  // Sometimes the sampler doesn't work well, in those cases we use a rejection
  // sampler. They do something similar in the truncnorm package.
  while(sample < lower | sample > upper | isinf(std::abs(sample))) {
    sample = sampleExtreme(mu, sd, lower, upper) ;
  }

  return sample ;
}

double computeConditionalMean(NumericVector mu,
                              NumericVector samp,
                              const NumericMatrix sigma,
                              int index) {
  double result = 0 ;

  for(int j = 0; j < mu.length() ; j ++) {
    if(j != index) {
      result += sigma(index, j) * (samp[j] - mu[j]) ;
    }
  }

  result = result / sigma(index, index) ;
  result = mu[index] - result ;
  return result ;
}

// [[Rcpp::export]]
NumericVector sampleTruncNorm(NumericVector sample,
                     NumericVector lower, NumericVector upper,
                     NumericVector mean, NumericMatrix sigma,
                     NumericVector condSigma,
                     int cycles) {
  double condMean ;
  double condSD ;
  int i, j ;
  NumericVector samp = clone(sample) ;

  for(i = 0; i < cycles ; i++) {
    for(j = 0; j < samp.length() ; j ++) {
      condMean = computeConditionalMean(mean, samp, sigma, j) ;
      condSD = std::sqrt(condSigma[j]) ;
      if((isinf(std::abs(lower[j])) & ((condMean - upper[j]) / condSD) > 3.5) |
         (isinf(std::abs(upper[j])) & ((lower[j] - condMean) / condSD) > 3.5)) {
        samp[j] = sampleExtreme(condMean, condSD, lower[j], upper[j]) ;
      } else {
        samp[j] = sampleUnivTruncNorm(condMean, condSD, lower[j], upper[j]) ;
      }
    }
  }

  return samp ;
}


