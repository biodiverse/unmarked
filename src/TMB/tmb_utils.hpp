template<class Type>
vector<Type> invcloglog(vector<Type> inp) {
  int sz = inp.size();
  vector<Type> out(sz);
  for (int i=0; i<sz; i++){
    out(i) = 1 - exp(-exp(inp(i)));
  }
  return out;
}

template<class Type>
vector<Type> add_ranef(vector<Type> par, Type& loglik, 
                vector<Type> b, Eigen::SparseMatrix<Type> Z,
                vector<Type> lsigma, int n_group_vars, vector<int> n_grouplevels) {
  
  if(n_group_vars == 0) return par;
  vector<Type> sigma = exp(lsigma);
  int idx = 0;
  for (int i=0; i<n_group_vars; i++){
    for (int j=0; j<n_grouplevels(i); j++){
      loglik += dnorm(b(idx), Type(0.0), sigma(i), true);
      idx += 1;
    }
  }
  par += Z * b;
  return par;
}

template<class Type>
bool is_na(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
bool all_na(vector<Type> x){
  for (int i = 0; i < x.size(); i++){
    if(!is_na(x(i))){
      return false;
    }
  }
  return true;
}

// Macro for submodel inputs
#define SUBMODEL_INPUTS(type) \
  DATA_MATRIX(X_##type) \
  DATA_SPARSE_MATRIX(Z_##type) \
  DATA_VECTOR(offset_##type) \
  DATA_INTEGER(n_group_vars_##type) \
  DATA_IVECTOR(n_grouplevels_##type) \
  PARAMETER_VECTOR(beta_##type) \
  PARAMETER_VECTOR(b_##type) \
  PARAMETER_VECTOR(lsigma_##type)
