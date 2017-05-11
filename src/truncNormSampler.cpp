#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sampleExtreme(double mu, double sd, double lower, double upper) {
  double sign = 1 ;
  double threshold ;
  double proposal ;
  double alpha ;
  double phi ;

  if(std::isinf(lower)) {
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
  int iter = 0;
  while(reject & (iter++ < 10000)) {
    proposal = threshold + R::rexp(1 / alpha) ;
    phi = std::exp(-std::pow(proposal - alpha, 2) / 2) ;
    if(runif(1)[0] < phi) {
      reject = false ;
    }
  }
  //Rcpp::Rcout<<iter<<"\n" ;

  proposal = proposal * sd + mu ;
  return proposal * sign;
}

// [[Rcpp::export]]
double sampleUnivTruncNorm(double mu, double sd, double lower, double upper) {
  double u = runif(1)[0] ;
  double phiB = R::pnorm5(upper, mu, sd, 1, 0) ;
  double phiA = R::pnorm5(lower, mu, sd, 1, 0) ;
  double quantile = u * phiB + phiA * (1 - u) ;
  double sample = R::qnorm(quantile, mu, sd, 1, 0) ;

  // Sometimes the sampler doesn't work well, in those cases we use a rejection
  // sampler. They do something similar in the truncnorm package.
  while(sample < lower | sample > upper | std::isinf(std::abs(sample))) {
    sample = sampleExtreme(mu, sd, lower, upper) ;
  }

  return sample ;
}

double computeConditionalMean(NumericVector mu,
                              NumericVector samp,
                              const NumericMatrix precision,
                              int index) {
  double result = 0 ;
  double cor = 0;

  for(int j = 0; j < mu.length() ; j ++) {
    if(j != index) {
      result += precision(index, j) * (samp[j] - mu[j]) ;
      cor += precision(index, j) ;
    }
  }

  result = result / precision(index, index) ;
  //Rcpp::Rcout<<cor<<" "<<result<<" "<< samp[index]<<" ";
  result = mu[index] - result ;

  return result ;
}

// [[Rcpp::export]]
NumericVector sampleTruncNorm(NumericVector sample,
                     NumericVector lower, NumericVector upper,
                     NumericVector mean, NumericMatrix precision,
                     int cycles) {
  double condMean ;
  double condSD ;
  int i, j ;
  int EXTREME_THRESHOLD = 4 ;
  NumericVector samp = clone(sample) ;

  for(i = 0; i < cycles ; i++) {
    for(j = 0; j < samp.length() ; j ++) {
      condMean = computeConditionalMean(mean, samp, precision, j) ;
      //Rcpp::Rcout<<mean[j]<<" "<<condMean<<"\n" ;
      condSD = std::sqrt(1 / precision(j, j)) ;
      if((std::isinf(std::abs(lower[j])) & ((condMean - upper[j]) / condSD) > EXTREME_THRESHOLD) |
         (std::isinf(std::abs(upper[j])) & ((lower[j] - condMean) / condSD) > EXTREME_THRESHOLD)) {
        samp[j] = sampleExtreme(condMean, condSD, lower[j], upper[j]) ;
      } else {
        samp[j] = sampleUnivTruncNorm(condMean, condSD, lower[j], upper[j]) ;
      }
    }
   //Rcpp::Rcout<<samp<<"\n" ;
  }

  return samp ;
}


