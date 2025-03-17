#include <RcppArmadillo.h>
#include <float.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nll_occu_Cpp(arma::colvec params, Rcpp::List inputs) {

  mat y = as<mat>(inputs["y"]);
  icolvec no_detect = 1 - as<icolvec>(inputs["Kmin"]);

  colvec idx_state = as<colvec>(inputs["idx_state"]) - 1;
  colvec beta_state = params.subvec(idx_state(0), idx_state(1));
  mat X_state = as<mat>(inputs["X_state"]);
  colvec offset_state = as<colvec>(inputs["offset_state"]);
  int invlink_state = inputs["invlink_state"];
  icolvec known_occ = as<icolvec>(inputs["known_occ_state"]);

  colvec idx_det = as<colvec>(inputs["idx_det"]) - 1;
  colvec beta_det = params.subvec(idx_det(0), idx_det(1));
  mat X_det = as<mat>(inputs["X_det"]);
  colvec offset_det = as<colvec>(inputs["offset_det"]);

  int M = y.n_rows;
  int J = y.n_cols;

  //Calculate psi
  colvec psi = X_state * beta_state + offset_state;
  if(invlink_state == 2){ //inverse cloglog
    psi = 1 - exp(-exp(psi));
  } else {
    psi = 1.0/(1.0+exp(-psi));
  }

  //Calculate p
  colvec p = X_det * beta_det + offset_det;
  p = 1.0/(1.0+exp(-p));

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
