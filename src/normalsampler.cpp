#include <Rcpp.h>
using namespace Rcpp;

double computeConditionalMean(NumericVector mu,
                              NumericVector samp,
                              NumericMatrix invcov,
                              int index) {
  double result = 0 ;

  for(int j = 0; j < mu.length() ; j ++) {
    if(j != index) {
      result += invcov(index, j) * (samp[j] - mu[j]) ;
    }
  }

  result = mu[index] - result / invcov(index, index) ;
  return result ;
}

double sampleExtreme(double mu, double sd, double threshold) {
  double sign = 1 ;
  double proposal ;
  double alpha ;
  double phi ;

  sign = -1 ;
  mu *= sign ;
  threshold = threshold * sign ;

  // rescaling
  threshold = (threshold - mu) / sd ;
  alpha = (threshold + sqrt(std::pow(threshold, 2) + 4)) / 2 ;

  bool reject = true ;
  int iter = 0;
  while(reject & (iter++ < 10000)) {
    proposal = threshold + R::rexp(1 / alpha) ;
    phi = exp(-std::pow(proposal - alpha, 2) / 2) ;
    if(runif(1)[0] < phi) {
      reject = false ;
    }
  }

  proposal = proposal * sd + mu ;
  return proposal * sign;
}

double sampleUnivTruncNorm(double mu, double sd, double threshold) {
  double u = runif(1)[0] ;
  double phiThreshold, sample ;

  // if(isnan(mu)) {
  //   Rcpp::Rcout<<"mu is nan \n" ;
  //   return 0 ;
  // }

  if((std::abs(mu - threshold) / sd) > 2.5) {
    return sampleExtreme(mu, sd, threshold) ;
  }

  phiThreshold = R::pnorm5(threshold, mu, sd, 1, 0) ;
  sample = R::qnorm5(u * phiThreshold, mu, sd, 1, 0) ;

  int tries = 0 ;
  while(isnan(sample) & tries ++ < 10) {
    sample =  sampleExtreme(mu, sd, threshold) ;
  }

  return sample ;
}

double sampleBoundedTruncNorm(double mu, double sd, double lower, double upper) {
  double u = runif(1)[0] ;
  double phiB = R::pnorm5(upper, mu, sd, 1, 0) ;
  double phiA = R::pnorm5(lower, mu, sd, 1, 0) ;
  double quantile = u * phiB + phiA * (1 - u) ;
  double sample = R::qnorm(quantile, mu, sd, 1, 0) ;
  return sample ;
}

// [[Rcpp::export]]
NumericVector mvtSampler(NumericVector y,
                         NumericVector mu,
                         IntegerVector selected,
                         NumericMatrix threshold,
                         NumericMatrix precision,
                         int nsamp, int burnin, int trim,
                         IntegerVector samporder,
                         bool verbose) {
  int totalIter = burnin + (nsamp - 1) * trim ;
  int frac = round(totalIter / 5) ;
  int p = y.length() ;
  NumericMatrix samples(nsamp, p) ;
  NumericVector samp = clone(y) ;
  double condmean, condsd, pprob, nprob, u ;
  int row = 0;

  for(int i = 0 ; i < (burnin + trim * nsamp) ; i ++) {
    //Rcpp::Rcout<<"\n"<<i<<" " ;
    for(int k = 0 ; k < p ; k ++) {
      int j = samporder[k] - 1 ;
      //Rcpp::Rcout<<j<<" " ;
      condmean = computeConditionalMean(mu, samp, precision, j) ;
      //Rcpp::Rcout<<condmean<<" "<<samp[j]<<"\n" ;
      condsd = 1 / std::sqrt(precision(j, j)) ;
      if(selected[j] == 1) {
        nprob = R::pnorm(threshold(j, 0), condmean, condsd, 1, 1) ;
        pprob = R::pnorm(threshold(j, 1), condmean, condsd, 0, 1) ;
        nprob = 1 / (1 + std::exp(pprob - nprob)) ;
        u = runif(1)[0] ;
        if(u < nprob) {
          samp[j] = sampleUnivTruncNorm(condmean, condsd, threshold(j, 0)) ;
        } else {
          samp[j] = sampleUnivTruncNorm(-condmean, condsd, -threshold(j, 1)) ;
          samp[j] = -samp[j] ;
        }
      } else {
        samp[j] = sampleBoundedTruncNorm(condmean, condsd, threshold(j, 0), threshold(j, 1)) ;
      }
    }

    if(verbose & ((((i + 1) % frac) == 0) | (i == totalIter - 1))) {
      int out = (i + 1.0) / (1.0 * totalIter) * 100 ;
      Rcpp::Rcout<<out<<"% " ;
    }

    if(i >= (burnin - 1) & (i - burnin + 1) % trim == 0) {
      for(int j = 0 ; j < samples.ncol() ; j++) {
        samples(row, j) = samp[j] ;
      }
      if(++row == nsamp) {
        break ;
      }
    }
  }

  if(verbose) Rcpp::Rcout<<"\n" ;

  return samples ;
}

bool checkSignEqual(double a, double b) {
  if(a < 0 & b > 0) {
    return false ;
  } else {
    return true ;
  }
}

// [[Rcpp::export]]
NumericVector modmvtSampler(NumericVector y,
                         NumericVector mu,
                         IntegerVector selected,
                         NumericMatrix threshold,
                         NumericMatrix precision,
                         int nsamp, int burnin, int trim,
                         bool verbose, int allowedNonModel) {
  int totalIter = burnin + (nsamp - 1) * trim ;
  int frac = round(nsamp / 5) ;
  int p = y.length() ;
  NumericMatrix samples(nsamp, p) ;
  NumericVector samp = clone(y) ;
  double condmean, condsd, pprob, nprob, u ;
  int row = 0;
  double proposal ;
  IntegerVector nonModel(p, 0) ;
  double tryCount = 0 ;
  double aprop = 0;
  int i = 0 ;
  int fromLast = 0 ;

  while(row <= nsamp) {
    i++ ;
    fromLast++ ;
    //Rcpp::Rcout<<"\n"<<i<<" " ;
    for(int j = 0 ; j < p ; j ++) {
      condsd = 1 / std::sqrt(precision(j, j)) ;
      condmean = computeConditionalMean(mu, samp, precision, j) ;
      if(selected[j] == 1) {
        bool accepted = false ;
        //Rcpp::Rcout<<condmean<<" " ;
        proposal = rnorm(1, condmean, condsd)[0] ;
        if(proposal < threshold(j, 0) || proposal > threshold(j, 1)) {
          samp[j] = proposal ;
          nonModel[j] = 0 ;
          accepted = true ;
        } else if(sum(nonModel) < allowedNonModel) {
          samp[j] = proposal ;
          nonModel[j] = 1;
          accepted = true ;
        } else if(nonModel[j] == 1) {
          samp[j] = proposal ;
          accepted = true ;
        }

        //Rcpp::Rcout<<samp[j]<<" " ;

        tryCount++ ;
        aprop *= (tryCount - 1.0) / tryCount ;
        if(accepted) {
           aprop += 1.0 / tryCount ;
        }

        if(sum(nonModel) == 0 & fromLast >= trim & i > burnin) {
          break ;
        }
      } else {
        samp[j] = sampleBoundedTruncNorm(condmean, condsd, threshold(j, 0), threshold(j, 1)) ;
      }
    }

    //Rcpp::Rcout<<sum(nonModel)<<" "<<fromLast<<"\n" ;
    if(verbose & ((((row + 1) % frac) == 0) | (row == nsamp))) {
      int out = (i + 1.0) / (1.0 * totalIter) * 100 ;
      //Rcpp::Rcout<<out<<"% " ;
    }

    if(i >= burnin & fromLast >= trim & sum(nonModel) == 0) {
      Rcpp::Rcout<<row<<" "<<fromLast<<" " ;
      for(int j = 0 ; j < samples.ncol() ; j++) {
        samples(row, j) = samp[j] ;
      }
      row++ ;
      fromLast = 0 ;
    }
    //Rcpp::Rcout<<fromLast<<" "<<sum(nonModel)<<" ";
  }

  if(verbose) Rcpp::Rcout<<"\n" ;

  return samples ;
}

