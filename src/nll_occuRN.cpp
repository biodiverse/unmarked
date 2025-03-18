#include <RcppArmadillo.h>
#include <float.h>
#include "utils.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

double lp_site_occuRN(const rowvec y, double lam, const vec q, int K, int Kmin){

  //y must be arma::vec for this to work (not uvec/ivec)
  const uvec fin = find_finite(y);
  if(fin.size() == 0) return 0.0;

  double f, g, p, out = 0.0;

  for (int k=Kmin; k<(K+1); k++){
    f = Rf_dpois(k, lam, 0);
    g = 0.0;
    for (unsigned j=0; j<fin.size(); j++){
      p = 1 - pow(q(fin(j)), k);
      g += Rf_dbinom(y(fin(j)), 1, p, true);
    }
    out += f * exp(g);
  }
  return log(out + DBL_MIN);
}

// [[Rcpp::export]]
double nll_occuRN_Cpp(arma::vec params, Rcpp::List inputs){

  mat y = as<mat>(inputs["y"]);
  uvec Kmin = as<uvec>(inputs["Kmin"]);
  int Kmax = inputs["Kmax"];

  SUBMODEL_INPUTS(state);
  SUBMODEL_INPUTS(det);

  int threads = inputs["threads"];

  int M = y.n_rows;
  int J = y.n_cols;

  const vec lam = exp(X_state * beta_state + offset_state);
  const vec q = 1 - inv_logit(X_det * beta_det + offset_det);

  #ifdef _OPENMP
    omp_set_num_threads(threads);
  #endif

  double loglik = 0.0;

  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for (int i=0; i<M; i++){
    int pstart = i * J;
    int pstop = i * J + J - 1;
    loglik += lp_site_occuRN(y.row(i), lam(i), q.subvec(pstart, pstop),
                             Kmax, Kmin(i));
  }

  return -loglik;

}
