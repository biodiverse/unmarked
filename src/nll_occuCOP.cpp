#include <RcppArmadillo.h>
#include <float.h>

using namespace Rcpp ;

// [[Rcpp::export]]
double nll_occuCOP(arma::icolvec y, arma::icolvec L, 
    arma::mat X, arma::mat V,
    arma::colvec beta_psi, arma::colvec beta_lambda,
    Rcpp::LogicalVector removed) {
  
  // Number of sites M and obs J
  int M = X.n_rows;
  int J = y.n_elem / M;

  //Calculate psi back-transformed from logit
  arma::colvec logit_psi = X*beta_psi;
  arma::colvec psi(M);
  psi = 1.0/(1.0+exp(-logit_psi));

  //Calculate lambda back-transformed from log
  arma::colvec log_lambda = V*beta_lambda;
  arma::colvec lambda = exp(log_lambda);
  

  double ll=0.0;
  int k=0; // counter
  // for each site i in 1:M
  for(int i=0; i<M; i++) {
    double iLambdaL=0.0; // init sum(lambda_ij * L_ij)
    double iN=0.0; // init sum(y) = total count of detec at site i
    for(int j=0; j<J; j++) {
      if(!removed(k)) {
        iLambdaL += lambda(k)*L(k);
        iN += y(k);
      }
      k++;
    }
    if(iN==0) {
      ll += log(psi(i) * pow(iLambdaL, iN) / tgamma(iN + 1) * exp(-iLambdaL));
    } else {
      ll += log(psi(i) * exp(-iLambdaL) + 1-psi(i));
    }
  }
  return -ll;
}
