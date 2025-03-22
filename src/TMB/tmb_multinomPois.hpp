#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
template <class Type>
Type tmb_multinomPois(objective_function<Type>* obj) {
  //Describe input data
  DATA_MATRIX(y); //observations
  int M = y.rows(); // # of sites
  int J = y.cols(); // # of obs per site
  
  SUBMODEL_INPUTS(state);

  SUBMODEL_INPUTS(det);
  DATA_INTEGER(R_det);
  DATA_INTEGER(pi_code_det);

  //Define the log likelihood so that it can be calculated in parallel over sites
  Type nll = 0.0;

  //Construct lambda vector
  vector<Type> lam = X_state * beta_state + offset_state;
  lam = add_ranef(lam, nll, b_state, Z_state, lsigma_state, 
                  n_group_vars_state, n_grouplevels_state);
  lam = exp(lam);

  //Construct p vector
  vector<Type> p = X_det * beta_det + offset_det;
  p = add_ranef(p, nll, b_det, Z_det, lsigma_det, 
                n_group_vars_det, n_grouplevels_det);
  p = invlogit(p);
  
  //Likelihood
  for (int m=0; m<M; m++){
    vector<Type> ysub = y.row(m);
    if(all_na(ysub)) continue;

    int pstart = m * R_det;
    vector<Type> psub = p.segment(pstart, R_det);
    vector<Type> pi = pifun(psub, pi_code_det);
    
    for (int j=0; j<J; j++){
      if(is_na(ysub(j))) continue;
      nll -= dpois(ysub(j), lam(m) * pi(j), true);
    }
  }

  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
