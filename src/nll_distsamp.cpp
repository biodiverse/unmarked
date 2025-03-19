#include <RcppArmadillo.h>
#include <float.h>
#include "distprob.h"
#include "utils.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nll_distsamp_Cpp(arma::vec params, Rcpp::List inputs) {

  mat y = as<mat>(inputs["y"]);
  int M = y.n_rows;
  int J = y.n_cols;

  SUBMODEL_INPUTS(state);
  vec A_state = inputs["A_state"];

  vec lambda = exp(X_state * beta_state + offset_state) % A_state;

  std::string keyfun_det = inputs["keyfun_det"];
  std::string survey_det = inputs["survey_det"];
  vec db_det = inputs["db_det"];
  vec w_det = inputs["w_det"];
  mat a_det = inputs["a_det"];
  mat u_det = inputs["u_det"];
  u_det = u_det.t(); // for compatability later

  vec sigma(M);
  if(keyfun_det != "uniform"){
    SUBMODEL_INPUTS(det);
    sigma = exp(X_det * beta_det + offset_det);
  }

  double scale;
  if(keyfun_det == "hazard"){
    ivec idx_scale = as<ivec>(inputs["idx_scale"]) - 1;
    scale = exp(params(idx_scale(0)));
  }

  double loglik = 0;
  uvec fin;
  vec p(J);
  for (int m=0; m<M; m++){
    fin = find_finite(y.row(m));
    if(fin.size() == 0) continue;

    p = distprob(keyfun_det, sigma(m), scale, survey_det, db_det, w_det, a_det.row(m));
    p = p % u_det.col(m);

    for (unsigned j=0; j<fin.size(); j++){
      loglik += Rf_dpois(y(m, fin(j)), lambda(m) * p(fin(j)), true);
    }
  }

  return -loglik;

}
