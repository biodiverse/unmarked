#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
template <class Type>
Type tmb_distsamp(objective_function<Type>* obj) {

  DATA_MATRIX(y); //observations
  int M = y.rows(); // # of sites
  int J = y.cols(); // # of distance categories per site

  SUBMODEL_INPUTS(state);
  DATA_VECTOR(A_state); // Area

  // Distance sampling info
  DATA_INTEGER(keyfun_code_det);
  DATA_INTEGER(survey_code_det);
  DATA_VECTOR(db_det);
  DATA_VECTOR(w_det);
  DATA_MATRIX(a_det);
  DATA_MATRIX(u_det);

  Type loglik = 0;

  //Construct lambda vector
  vector<Type> lambda = X_state * beta_state + offset_state;
  lambda = add_ranef(lambda, loglik, b_state, Z_state, lsigma_state, 
                     n_group_vars_state, n_grouplevels_state);
  lambda = exp(lambda);
  lambda = lambda.array() * A_state.array();

  //Construct sigma vector if needed
  vector<Type> sigma(M);
  if(keyfun_code_det != 0){ // If not uniform
    SUBMODEL_INPUTS(det);
    sigma = X_det * beta_det + offset_det;
    sigma = add_ranef(sigma, loglik, b_det, Z_det, lsigma_det, 
                      n_group_vars_det, n_grouplevels_det);
    sigma = exp(sigma);
  }

  Type scale; 
  if(keyfun_code_det == 3){ // If hazard
    PARAMETER_VECTOR(beta_scale);
    scale = exp(beta_scale(0));
  }

  //Likelihood
  for (int m=0; m<M; m++){
    vector<Type> ysub = y.row(m);
    if(all_na(ysub)) continue;
    
    //Not sure if defining this inside loop is necessary for parallel
    vector<Type> asub = a_det.row(m);
    vector<Type> usub = u_det.row(m);
    vector<Type> cp = distance_prob(keyfun_code_det, sigma(m), scale, 
                                    survey_code_det, db_det, w_det, asub, usub); 

    Type site_lp = 0;
    
    for (int j=0; j<J; j++){
      if(is_na(ysub(j))) continue; //If NA found skip
      site_lp += dpois(ysub(j), lambda(m) * cp(j), true);
    }
    loglik += site_lp;
  }

  return -loglik;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
