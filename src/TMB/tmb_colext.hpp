#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
template <class Type>
Type tmb_colext(objective_function<Type>* obj) {

  DATA_MATRIX(y);
  DATA_MATRIX(Kmin);
  matrix<Type> no_detect = 1 - Kmin.array();
  DATA_MATRIX(site_sampled);

  Type loglik = 0.0;

  SUBMODEL_INPUTS(psi);
  int M = X_psi.rows();
  vector<Type> psi = X_psi * beta_psi + offset_psi;
  psi = add_ranef(psi, loglik, b_psi, Z_psi, lsigma_psi,
                  n_group_vars_psi, n_grouplevels_psi);
  psi = invlogit(psi);
 
  SUBMODEL_INPUTS(col);
  DATA_INTEGER(T_col);
  DATA_INTEGER(J_col);
  int T = T_col;
  int J = J_col;
  vector<Type> col = X_col * beta_col + offset_col;
  col = add_ranef(col, loglik, b_col, Z_col, lsigma_col,
                  n_group_vars_col, n_grouplevels_col);
  col = invlogit(col);

  SUBMODEL_INPUTS(ext);
  vector<Type> ext = X_ext * beta_ext + offset_ext;
  ext = add_ranef(ext, loglik, b_ext, Z_ext, lsigma_ext,
                  n_group_vars_ext, n_grouplevels_ext);
  ext = invlogit(ext);

  SUBMODEL_INPUTS(det);
  vector<Type> p = X_det * beta_det + offset_det;
  p = add_ranef(p, loglik, b_det, Z_det, lsigma_det,
                n_group_vars_det, n_grouplevels_det);
  p = invlogit(p);

  int Tidx = 0;
  int yind;
  int pind = 0;
  vector<Type> ys(y.cols());

  vector<Type> Dy(2);
  matrix<Type> phi_total(2, 2);
  matrix<Type> phi(2, 2);
  vector<Type> phi_end(2);
  vector<Type> psi_vec(2);
  matrix<Type> Dy_diag(2,2);
  Dy_diag = Dy_diag.setZero();

  for (int i=0; i<M; i++){
    yind = 0;
    ys = y.row(i);
    if(all_na(ys)){
      pind += J*T;
      Tidx += (T-1);
      continue;
    }
    phi_total = phi.setIdentity();

    //All but the last time period
    for (int t=0; t<(T-1); t++){

      phi = phi.setZero();
      phi(0, 0) = 1 - ext(Tidx);
      phi(0, 1) = ext(Tidx);
      phi(1, 0) = col(Tidx);
      phi(1, 1) = 1 - col(Tidx);
      Tidx += 1;

      if(site_sampled(i, t) == 1){
        Dy = Dy.setOnes();

        if(no_detect(i, t) == 0){ // if at least one detection
          Dy(1) = 0;
        }

        for (int j=0; j<J; j++){
          if(is_na(ys(yind))){
            yind += 1;
            pind += 1;
            continue;
          }

          Dy(0) *= pow(p(pind), ys(yind)) * pow(1-p(pind), 1-ys(yind));
          yind += 1;
          pind += 1;
        }
        Dy_diag.diagonal() = Dy;
        //phi_total = phi_total * phi * Dy_diag;
        phi_total = phi_total * Dy_diag * phi;
      } else {
        phi_total = phi_total * phi;
        yind += J;
        pind += J;
      }
    }
    
    //Final period; no psi calculation
    if(site_sampled(i, (T-1)) == 1){
      Dy = Dy.setOnes();

      if(no_detect(i, (T-1)) == 0){ // if at least one detection
        Dy(1) = 0;
      }

      for (int j=0; j<J; j++){
        if(is_na(ys(yind))){
          yind += 1;
          pind += 1;
          continue;
        }

        Dy(0) *= pow(p(pind), ys(yind)) * pow(1-p(pind), 1-ys(yind));
        yind += 1;
        pind += 1;
      }

      phi_end = phi_total * Dy;
    } else {
      Dy = Dy.setOnes();
      phi_end = phi_total * Dy;
      yind += J;
      pind += J;
    }

    psi_vec(0) = psi(i);
    psi_vec(1) = 1 - psi(i);
    // Dot product
    vector<Type> site_lik = psi_vec.array() * phi_end.array();
    loglik += log(sum(site_lik));
  }

  return(-loglik);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
