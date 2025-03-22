
#' Fit the N-mixture point count model

pcount <- function(formula, data, K = NULL, mixture = c("P", "NB", "ZIP"), 
                   starts = NULL, method = "BFGS", se = TRUE, 
                   engine = c("C", "R", "TMB"), threads = 1, ...){
  
  # Check arguments
  if(!is(data, "unmarkedFramePCount")){
    stop("Data is not an unmarkedFramePCount object.")
  }
  engine <- match.arg(engine)
  forms <- split_formula(formula)
  if(any(sapply(forms, has_random))) engine <- "TMB"
  mixture <- match.arg(mixture)
  if(identical(mixture, "ZIP") & engine == "R"){
    stop("ZIP mixture not available for R engine")
  }

  # Build submodels
  fam <- switch(mixture, P = "poisson", NB = "negative_binomial", ZIP = "ZIP")
  submodels <- unmarkedSubmodelList(
    state = unmarkedSubmodelState(name = "Abundance", short_name = "lam", 
                                  type = "state", formula = forms[[2]], data = data, 
                                  family = fam, link = "log"),

    det = unmarkedSubmodelDet(name = "Detection", short_name = "p", 
                              type = "det", formula = forms[[1]], data = data, 
                              family = "binomial", link = "logit")
  )

  if(mixture == "NB"){
    submodels['alpha'] <- unmarkedSubmodelScalar(name = "Dispersion", 
                                                 short_name = "alpha", 
                                                 type = "alpha", link = "log")
  } else if(mixture == "ZIP"){
    submodels['psi'] <- unmarkedSubmodelScalar(name = "Zero-inflation", 
                                               short_name = "psi", 
                                               type = "psi", link = "logit")
  }

  # Build response object
  response <- unmarkedResponseCount(data, submodels, Kmax = K)

  # Get nll inputs
  inputs <- nll_inputs(response, submodels, engine)
  inputs$threads <- threads
  nll_fun <- switch(engine, R = nll_pcount_R, C = nll_pcount_Cpp, TMB = "tmb_pcount")

  # Fit model
  fit <- fit_model(nll_fun, inputs = inputs, submodels = submodels,
                   starts = starts, method = method, se = se, ...)
   
  # Create unmarkedFit object
  new("unmarkedFitPCount", fitType="pcount", call=match.call(),
      formula = formula, data = data, sitesRemoved = removed_sites(response),
      estimates = fit$submodels, AIC = fit$AIC, opt = fit$opt,
      negLogLike = fit$opt$value, nllFun = fit$nll, K = response@Kmax, 
      mixture = mixture, TMB=fit$TMB)
}


nll_pcount_R <- function(params, inputs){
  with(inputs, {

  beta_state <- params[idx_state[1]:idx_state[2]]
  beta_det <- params[idx_det[1]:idx_det[2]]

  M <- nrow(y)
  J <- ncol(y)

  k <- 0:Kmax
  k_ik <- rep(k, M)
  k_ijk <- rep(k, M*J)
  y_ij <- as.numeric(t(y))
  y_ijk <- rep(y_ij, each = Kmax + 1)

  ijk <- expand.grid(k = 0:Kmax, j = 1:J, i = 1:M)
  ijk_to_ikj <- with(ijk, order(i, k, j))

  theta_i <- exp(X_state %*% beta_state + offset_state)
  p_ij <- plogis(X_det %*% beta_det + offset_det)
  theta_ik <- rep(theta_i, each = Kmax + 1)
  p_ijk <- rep(p_ij, each = Kmax + 1)

  bin_ijk <- dbinom(y_ijk, k_ijk, p_ijk)
  bin_ijk[which(is.na(bin_ijk))] <- 1
  
  bin_ik_mat <- matrix(bin_ijk[ijk_to_ikj], M * (Kmax + 1), J,
                       byrow = TRUE)
  g_ik <- rowProds(bin_ik_mat)

  if(family_state == 1) {
    f_ik <- dpois(k_ik, theta_ik)
  } else if(family_state == 2){
    alpha <- exp(params[idx_alpha[1]])
    f_ik <- dnbinom(k_ik, mu = theta_ik, size = alpha)
  }
  
  dens_i_mat <- matrix(f_ik * g_ik, M, Kmax + 1, byrow = TRUE)
  dens_i <- rowSums(dens_i_mat)  # sum over the K

  -sum(log(dens_i), na.rm = TRUE)

  })
}
