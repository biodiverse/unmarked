#include <RcppArmadillo.h>
#include <float.h>
#include "distr.h"
#include "utils.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

double lp_site_pcount(const rowvec y, int mixture, double lam, double log_alpha,
                      const vec p, int K, int Kmin){

  //y must be arma::vec for this to work (not uvec/ivec)
  uvec fin = find_finite(y);
  if(fin.size() == 0) return 0.0;

  double f, g, out = 0.0;
  for (int k=Kmin; k<(K+1); k++){
    f = N_density(mixture, k, lam, log_alpha);
    g = 0.0;
    for (unsigned j=0; j<fin.size(); j++){
      g += Rf_dbinom(y(fin(j)), k, p(fin(j)), true);
    }
    out += f * exp(g);
  }
  return log(out + DBL_MIN);
}


// [[Rcpp::export]]
double nll_pcount_Cpp(arma::vec params, Rcpp::List inputs){

  mat y = as<mat>(inputs["y"]);
  uvec Kmin = as<uvec>(inputs["Kmin"]);
  int Kmax = inputs["Kmax"];

  SUBMODEL_INPUTS(state);
  int family_state = inputs["family_state"];
  SUBMODEL_INPUTS(det);

  int threads = inputs["threads"];

  int M = y.n_rows;
  int J = y.n_cols;

  vec lam = exp(X_state * beta_state + offset_state);
  vec p = inv_logit(X_det * beta_det + offset_det);

  double par2 = 0;
  if(family_state == 2){
    ivec idx = as<ivec>(inputs["idx_alpha"]) - 1;
    par2 = params(idx(0));
  } else if(family_state == 3){
    ivec idx = as<ivec>(inputs["idx_psi"]) - 1;
    par2 = params(idx(0));
  }

  #ifdef _OPENMP
    omp_set_num_threads(threads);
  #endif

  double loglik = 0.0;

  //This will compile but throw an unknown pragma warning if no openMP
  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for (int i=0; i<M; i++){
    int pstart = i * J;
    int pstop = i * J + J - 1;
    loglik += lp_site_pcount(y.row(i), family_state, lam(i), par2,
                             p.subvec(pstart, pstop), Kmax, Kmin(i));
  }

  return -loglik;

}
