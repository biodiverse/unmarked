#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
template <class Type>
Type tmb_occu(objective_function<Type>* obj) {
  //Describe input data
  DATA_MATRIX(y); //observations
  DATA_VECTOR(Kmin); //Minimum value recorded at site
  vector<Type> no_detect = 1 - Kmin;

  SUBMODEL_INPUTS(state);
  DATA_INTEGER(invlink_state);
  DATA_IVECTOR(known_occ_state);

  SUBMODEL_INPUTS(det);

  Type loglik = 0.0;

  int M = y.rows(); //# of sites
  int J = y.cols(); //# of observations per site

  //Construct psi vector
  vector<Type> psi = X_state * beta_state + offset_state;
  psi = add_ranef(psi, loglik, b_state, Z_state, lsigma_state,
                  n_group_vars_state, n_grouplevels_state);
  if(invlink_state == 2){
    psi = invcloglog(psi);
  } else {
    psi = invlogit(psi);
  }

  //Construct p vector
  vector<Type> p = X_det * beta_det + offset_det;
  p = add_ranef(p, loglik, b_det, Z_det, lsigma_det,
                n_group_vars_det, n_grouplevels_det);
  p = invlogit(p);

  //Standard occupancy likelihood calculation
  for (int i=0; i<M; i++){
    if(is_na(psi(i))) continue;
    int pind = i * J;
    int all_na = true;
    Type cp = 1.0;
    for (int j=0; j<J; j++){
      if(is_na(y(i,j)) | is_na(p(pind))){
        pind += 1;
        continue;
      }
      cp *= pow(p(pind), y(i,j)) * pow(1-p(pind), 1-y(i,j));
      pind += 1;
      all_na = false;
    }
    if(all_na) continue;
    if(known_occ_state(i)) psi(i) = 1.0;
    loglik += log(psi(i) * cp + (1-psi(i)) * no_detect(i));
  }

  return -loglik;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
