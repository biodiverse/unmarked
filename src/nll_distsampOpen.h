#ifndef _unmarked_NLL_DISTSAMPOPEN_H
#define _unmarked_NLL_DISTSAMPOPEN_H

#include <RcppArmadillo.h>

RcppExport SEXP nll_distsampOpen( SEXP y_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, SEXP Xsig_, SEXP beta_lam_, SEXP beta_gam_, SEXP beta_om_, SEXP beta_sig_, SEXP log_alpha_, SEXP Xlam_offset_, SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xsig_offset_, SEXP ytr_, SEXP yr_, SEXP lk_, SEXP mixture_, SEXP first_, SEXP last_, SEXP M_, SEXP J_, SEXP T_, SEXP delta_, SEXP dynamics_, SEXP fix_, SEXP go_dims_, SEXP I_, SEXP I1_, SEXP Ib_, SEXP Ip_, SEXP a_, SEXP u_, SEXP db_ ) ;



// constant model
//void tp1(arma::mat& g3, int lk, double gam, double om);
void tp1(arma::mat& g3, int nrI, int nrI1, Rcpp::IntegerVector N, arma::imat I, arma::imat I1, Rcpp::List Ib, Rcpp::List Ip, double gam, double om);



// autoregressive model
void tp2(arma::mat& g3, int lk, double gam, double om);



// trend model
void tp3(arma::mat& g3, int lk, double gam);


#endif
