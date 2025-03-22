
# Fit the Occupancy model of Royle and Nichols

occuRN <- function(formula, data, K = 25, starts = NULL, method = "BFGS",
                   se = TRUE, engine=c("C","R"), threads = 1, ...){

  # Check arguments
  if(!is(data, "unmarkedFrameOccu")) stop("Data is not an unmarkedFrameOccu object.")
  engine <- match.arg(engine)
  forms <- split_formula(formula)
  check_no_support(forms)
    
  # Build submodels
  submodels <- unmarkedSubmodelList(
    state = unmarkedSubmodelState(name = "Abundance", short_name = "lam", 
                                  type = "state", formula = forms[[2]], 
                                  data = data, family = "poisson", link = "log"),

    det = unmarkedSubmodelDet(name = "Detection", short_name = "p", 
                              type = "det", formula = forms[[1]], data = data, 
                              family = "binomial", link = "logit")
  )

  # Build response object
  response <- unmarkedResponseBinary(data, submodels, Kmax = K)

  # Get nll inputs
  inputs <- nll_inputs(response, submodels, engine)
  inputs$threads <- threads
  nll_fun <- switch(engine, R = nll_occuRN_R, C = nll_occuRN_Cpp)

  # Fit model
  fit <- fit_model(nll_fun, inputs = inputs, submodels = submodels,
                   starts = starts, method = method, se = se, ...)

  new("unmarkedFitOccuRN", fitType = "occuRN",
      call = match.call(), formula = formula, data = data,
      sitesRemoved = removed_sites(response), estimates = fit$submodels,
      AIC = fit$AIC, opt = fit$opt, negLogLike = fit$opt$value, 
      nllFun = fit$nll, K = response@Kmax)
}


nll_occuRN_R <- function(params, inputs){
  with(inputs, {

  beta_state <- params[idx_state[1]:idx_state[2]]
  beta_det <- params[idx_det[1]:idx_det[2]]
  
  M <- nrow(inputs$y)
  J <- ncol(inputs$y)
  
  # compute individual level detection probabilities
  r_ij <- matrix(plogis(X_det %*% beta_det + offset_det), M, J,
                 byrow = TRUE)
  
  # compute list of detection probabilities along N
  n <- 0:Kmax
  p_ij_list <- lapply(n, function(k) 1 - (1 - r_ij)^k)

  # compute P(y_{ij} | N) (cell probabilities) along N
  cp_ij_list <- lapply(p_ij_list, function(pmat) pmat^y * (1-pmat)^(1-y))

  # multiply across J to get P(y_i | N) along N
  cp_in <- sapply(cp_ij_list, rowProds, na.rm = TRUE)

  # compute P(N = n | lambda_i) along i
  lambda_i <- exp(X_state %*% beta_state + offset_state)
  lambda_in <- sapply(n, function(x) dpois(x, lambda_i))

  # integrate over P(y_i | N = n) * P(N = n | lambda_i) wrt n
  like_i <- rowSums(cp_in * lambda_in)
  -sum(log(like_i), na.rm = TRUE)
  })
}
