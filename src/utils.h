#ifndef _unmarked_UTILS_H
#define _unmarked_UTILS_H

#include <RcppArmadillo.h>

arma::mat inv_logit( arma::mat inp );

arma::vec inv_logit( arma::vec inp );

double inv_logit(double x);

arma::vec beta_sub(arma::vec beta, arma::uvec n_param, unsigned idx);

double dmultinom(arma::vec x, arma::vec prob);

// Macro for submodel inputs
#define STR(x) #x

#define SUBMODEL_INPUTS(type) \
  vec idx_##type = as<vec>(inputs["idx_" STR(type)]) - 1; \
  vec beta_##type = params.subvec(idx_##type(0), idx_##type(1)); \
  mat X_##type = as<mat>(inputs["X_" STR(type)]); \
  vec offset_##type = as<vec>(inputs["offset_" STR(type)]);

#endif
