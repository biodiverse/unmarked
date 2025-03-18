#include <RcppArmadillo.h>
#include <float.h>
#include "utils.h"

using namespace Rcpp;
using namespace arma;

arma::vec inv_logit( arma::vec inp ){
  return(1 / (1 + exp(-1 * inp)));
}

// [[Rcpp::export]]
double nll_occu_Cpp(arma::vec params, Rcpp::List inputs) {

  mat y = as<mat>(inputs["y"]);
  ivec no_detect = 1 - as<ivec>(inputs["Kmin"]);

  SUBMODEL_INPUTS(state);
  int invlink_state = inputs["invlink_state"];
  ivec known_occ = as<ivec>(inputs["known_occ_state"]);

  SUBMODEL_INPUTS(det);

  int M = y.n_rows;
  int J = y.n_cols;

  //Calculate psi
  vec psi = X_state * beta_state + offset_state;
  if(invlink_state == 2){ //inverse cloglog
    psi = 1 - exp(-exp(psi));
  } else {
    psi = inv_logit(psi);
  }

  //Calculate p
  vec p = X_det * beta_det + offset_det;
  p = inv_logit(p);

  double nll = 0;
  //colvec nll = zeros(M);
  bool all_is_na;
  for (int m=0; m<M; m++){
    int pind = m * J;
    double cp = 1.0;

    all_is_na = true;
    for (int j=0; j<J; j++){
      if(!is_finite(y(m,j))){
        pind += 1;
        continue;
      }
      cp *= pow(p(pind), y(m,j)) * pow(1-p(pind), 1-y(m,j));
      pind += 1;
      all_is_na = false;
    }
    if(all_is_na) continue;

    if(known_occ(m)) psi(m) = 1;
    nll -= log(psi(m) * cp + (1-psi(m)) * no_detect(m));
  }

  return nll;
}


// [[Rcpp::export]]
double nll_occu(arma::icolvec y, arma::mat X, arma::mat V,
    arma::colvec beta_psi, arma::colvec beta_p,
    Rcpp::IntegerVector nd, Rcpp::LogicalVector knownOcc, Rcpp::LogicalVector navec,
    arma::colvec X_offset, arma::colvec V_offset, std::string link_psi) {

  int R = X.n_rows;
  int J = y.n_elem / R;

  //Calculate psi
  arma::colvec psi_lp = X*beta_psi + X_offset;
  arma::colvec psi(R);
  if(link_psi == "cloglog"){
    psi = 1 - exp(-exp(psi_lp));
  } else {
    psi = 1.0/(1.0+exp(-psi_lp));
  }

  //Calculate p
  arma::colvec logit_p = V*beta_p + V_offset;
  arma::colvec p = 1.0/(1.0+exp(-logit_p));

  double ll=0.0;
  int k=0; // counter
  for(int i=0; i<R; i++) {
    double cp=1.0;
    for(int j=0; j<J; j++) {
      if(!navec(k))
	cp *= pow(p(k),y(k)) * pow(1-p(k), 1-y(k));
      k++;
    }
    if(knownOcc(i))
      psi(i) = 1.0;
    if(nd(i)==0)
      ll += log(cp * psi(i) + DBL_MIN);
    else if(nd(i)==1)
      ll += log(cp * psi(i) + (1-psi(i)) + DBL_MIN);
  }
  return -ll;
}
