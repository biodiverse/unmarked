
distsamp <- function(formula, data,
    keyfun=c("halfnorm", "exp", "hazard", "uniform"),
    output=c("density", "abund"), unitsOut=c("ha", "kmsq"), starts = NULL,
    method="BFGS", se = TRUE, engine = c("C", "R", "TMB"), ...){

  # Check arguments
  if(!is(data, "unmarkedFrameDS")){
    stop("Data is not an unmarkedFrameDS object.")
  }  
  engine <- match.arg(engine)
  forms <- split_formula(formula)
  if(any(sapply(split_formula(formula), has_random))) engine <- "TMB"
  keyfun <- match.arg(keyfun)
  output <- match.arg(output)
  unitsOut <- match.arg(unitsOut)

  # Build submodels
  state_name <- switch(output, abund = "Abundance", density = "Density")
  A <- get_ds_area(data, unitsOut, output)
  submodels <- unmarkedSubmodelList(
    state = unmarkedSubmodelState(name = state_name, short_name = "lam", 
              type = "state", formula = forms[[2]], data = data, 
              family = "poisson", link = "log", 
              auxiliary = list(output = output, unitsOut = unitsOut, A = A)),

    det = unmarkedSubmodelDistance(name = "Detection", short_name = "p", 
                                   type = "det", formula = forms[[1]], 
                                   data = data, keyfun = keyfun, link = "log")
  )

  if(keyfun == "hazard"){
    submodels["scale"] <- unmarkedSubmodelScalar(name = "Hazard-rate (scale)", 
                            short_name = "p", type = "scale", link = "log")
    # for backwards compatability, initialize log(scale) at 1
  }
  
  # Build response object
  response <- unmarkedResponseCount(data, submodels, Kmax = NULL)

  # Get nll inputs
  inputs <- nll_inputs(response, submodels, engine)
  nll_fun <- switch(engine, R = nll_distsamp_R, C = nll_distsamp_Cpp, 
                    TMB = "tmb_distsamp")

  # Fit model
  fit <- fit_model(nll_fun, inputs = inputs, submodels = submodels,
                   starts = starts, method = method, se = se, ...)

  new("unmarkedFitDS", fitType = "distsamp", call = match.call(),
      opt = fit$opt, formula = formula, data = data, keyfun=keyfun,
      sitesRemoved = removed_sites(response), unitsOut=unitsOut,
      estimates = fit$submodels, AIC = fit$AIC, negLogLike = fit$opt$value,
      nllFun = fit$nll, output=output, TMB=fit$TMB)
}


nll_distsamp_R <- function(params, inputs){
  with(inputs, {

  M <- nrow(y)
  J <- ncol(y)
  
  beta_state <- params[idx_state[1]:idx_state[2]]
  lambda <- drop(exp(X_state %*% beta_state + offset_state)) * A_state

  sigma <- rep(NA, M)
  if(keyfun_det != "uniform"){
    beta_det <- params[idx_det[1]:idx_det[2]]
    sigma <- exp(X_det %*% beta_det + offset_det)
  }

  scale <- NA
  if(keyfun_det == "hazard"){
    scale <- exp(params[idx_scale[1]])
  }

  cp <- matrix(NA, M, J) 
  for (i in 1:M){
    cp[i,] <- getDistCP(keyfun_det, sigma[i], scale, survey_det, db_det,
                        w_det, a_det[i,], u_det[i,])
  }

  ll <- dpois(y, lambda * cp, log=TRUE)
  -sum(ll, na.rm = TRUE)
  })
}
