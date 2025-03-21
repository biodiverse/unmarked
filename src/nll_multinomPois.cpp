#include <RcppArmadillo.h>
#include "pifun.h"
#include "utils.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nll_multinomPois_Cpp(arma::vec params, Rcpp::List inputs){

  mat y = as<mat>(inputs["y"]);
  int M = y.n_rows;
  int J = y.n_cols;

  SUBMODEL_INPUTS(state);

  SUBMODEL_INPUTS(det);
  int R = inputs["R_det"];
  std::string pi_name = inputs["pi_name_det"];

  vec lambda = exp(X_state * beta_state + offset_state);

  vec p = inv_logit(X_det * beta_det + offset_det);

  double loglik = 0;
  int p_start, p_stop;
  uvec fin;
  vec pi(J);

  for (int m=0; m<M; m++){
    fin = find_finite(y.row(m));
    if(fin.size() == 0) continue;

    p_start = m * R;
    p_stop = p_start + R - 1;
    pi = piFun(p.subvec(p_start, p_stop), pi_name);

    for (unsigned j=0; j<fin.size(); j++){
      loglik += Rf_dpois(y(m, fin(j)), lambda(m) * pi(fin(j)), true);
    }
  }

  return -loglik;
}
