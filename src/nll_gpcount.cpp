#include <RcppArmadillo.h>
#include "distr.h"
#include "utils.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nll_gpcount_Cpp(arma::vec params, Rcpp::List inputs){

  #ifdef _OPENMP
    int threads = inputs["threads"];
    omp_set_num_threads(threads);
  #endif

  mat y = as<mat>(inputs["y"]);
  int Kmax = inputs["Kmax"];
  uvec Kmin = as<uvec>(inputs["Kmin"]);

  int M = y.n_rows;
  int T = inputs["T_phi"];
  int J = inputs["J_phi"];

  SUBMODEL_INPUTS(lambda);
  int family_lambda = inputs["family_lambda"];
  vec lam = exp(X_lambda * beta_lambda + offset_lambda);

  double par2 = 0;
  if(family_lambda == 2){
    ivec idx = as<ivec>(inputs["idx_alpha"]) - 1;
    par2 = params(idx(0));
  } else if(family_lambda == 3){
    ivec idx = as<ivec>(inputs["idx_psi"]) - 1;
    par2 = params(idx(0));
  }

  SUBMODEL_INPUTS(phi);
  vec phi = inv_logit(X_phi * beta_phi + offset_phi);

  SUBMODEL_INPUTS(det);
  vec p = inv_logit(X_det * beta_det + offset_det);

  //Indices etc.
  int phi_idx;
  double phi_t;
  int y_start;
  int y_end;
  int p_start;
  int p_end;
  rowvec y_m(T*J);
  rowvec y_t(J);
  vec p_t(J);
  uvec fin;

  double loglik = 0;
  //Intermediate sums etc.
  double f;
  double lik_site;
  double ll_k;
  double Nt_lik;

  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for (int m=0; m<M; m++){
    lik_site = 0;
    y_m = y.row(m);

    fin = find_finite(y_m);
    if(fin.size() == 0) continue;

    //Iterate over possible true population size K
    for (int k=Kmin(m); k<(Kmax+1); k++){

      y_start = 0;
      p_start = m*T*J;
      phi_idx = m*T;

      ll_k = log(N_density(family_lambda, k, lam(m), par2));

      for (int t=0; t<T; t++){

        //Extract data for period T and increment
        phi_t = phi(phi_idx);
        phi_idx += 1;

        y_end = y_start + J - 1;
        p_end = p_start + J - 1;
        y_t = y_m.subvec(y_start, y_end);
        p_t = p.subvec(p_start, p_end);
        y_start += J;
        p_start += J;

        //If no obs for this period, skip
        fin = find_finite(y_t);
        if(fin.size() == 0) continue;

        //Iterate over possible # of available animals N(t)
        //Cannot be smaller than max observed in period t
        Nt_lik = 0;
        for (int n=max(y_t); n<(k+1); n++){
          //Availability log-likelihood
          f = 0;
          if(is_finite(phi_t)){
            f = Rf_dbinom(n, k, phi_t, true);
          }
          //Detection log-likelihood
          for (unsigned j=0; j<fin.size(); j++){
            f += Rf_dbinom(y_t(fin(j)), n, p_t(fin(j)), true);
          }
          Nt_lik += exp(f);
        }
        ll_k += log(Nt_lik);
      }
      lik_site += exp(ll_k);
    }
    loglik += log(lik_site);
  }

  return -loglik;
}
