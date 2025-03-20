#include <RcppArmadillo.h>
#include <float.h>
#include "utils.h"

using namespace Rcpp;
using namespace arma;

arma::vec inv_logit( arma::vec inp ){
  return(1 / (1 + exp(-1 * inp)));
}

double calculate_occu_penalty(arma::vec beta_state, arma::vec beta_det,
    Rcpp::List inputs){

  int type = inputs["pen_type"];
  if(type == 0) return(0);

  double lambda = inputs["pen_lambda"];

  double penalty = 0;

  if(type == 1){ // Bayes
    penalty += accu(pow(beta_state, 2)) + accu(pow(beta_det, 2));
    penalty = penalty * lambda * 0.5;

  } else if(type == 2){ // Ridge
    bool has_int_state = inputs["has_int_state"];
    bool has_int_det = inputs["has_int_det"];
    if(has_int_state & (beta_state.size() > 1)){
      penalty += accu(pow(beta_state.subvec(1, beta_state.size()-1), 2));
    } else if(!has_int_state){
      penalty += accu(pow(beta_state, 2));
    }
    if(has_int_det & (beta_det.size() > 1)){
      penalty += accu(pow(beta_det.subvec(1, beta_det.size()-1), 2));
    } else if(!has_int_det){
      penalty += accu(pow(beta_det, 2));
    }
    penalty = penalty * lambda * 0.5;

  } else if(type ==3){ //MPLE
    vec LR_est = as<vec>(inputs["LR_est"]);
    penalty = accu(abs(beta_state - LR_est)) * lambda;
  }

  return penalty;
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

  double penalty = calculate_occu_penalty(beta_state, beta_det, inputs);
  nll += penalty;

  return nll;
}
