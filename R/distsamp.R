
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
  A <- rep(1, numSites(data))
  if(output == "density"){
    A <- rep(get_ds_area(data, unitsOut), numSites(data))
  }
  submodels <- unmarkedSubmodelList(
    state = unmarkedSubmodelState(name = state_name, short_name = "lam", 
              type = "state", formula = forms[[2]], data = data, 
              family = "poisson", link = "log", 
              auxiliary = list(output = output, unitsOut = unitsOut, A = A))
  )

  if(keyfun != "uniform"){
    submodels["det"] <- unmarkedSubmodelDistance(name = "Detection", 
                          short_name = "p", type = "det", 
                          formula = forms[[1]], data = data, 
                          keyfun = keyfun, link = "log")
  }
  if(keyfun == "hazard"){
    submodels['scale'] <- unmarkedSubmodelScalar(name = "Hazard-rate (scale)", 
                            short_name = "p", type = "alpha", link = "log")
    #if(is.null(starts)){
    #  starts <- default_starts(submodels)
      # for backwards compatability, initialize log(scale) at 1
      # it's not clear this is actually better than 0
      # note also that exp intercept used to be initialized at 0 
      # instead of log(median(db))
    #  starts[length(starts)] <- 1
    #}
  }
  
  # Build response object
  response <- unmarkedResponse(data, Kmax = NULL)
  # Handle missing values in covariates
  response <- add_missing(response, submodels)

  # Generate engine inputs
  if(engine == "TMB"){
    inputs <- engine_inputs_TMB(response, submodels)
  } else {
    inputs <- engine_inputs_CR(response, submodels)
  }
  # There's no submodel for a uniform key function, but we still need this info
  if(keyfun == "uniform"){
    ua <- getUA(data)
    inputs <- c(inputs, list(keyfun_det = "uniform", survey_det = data@survey,
                             db_det = data@dist.breaks, w_det = diff(data@dist.breaks),
                             u_det=ua$u, a_det=ua$a))
  }

  # Fit model
  nll_fun <- switch(engine, R = nll_distsamp_R, C = nll_distsamp_Cpp, 
                    TMB = "tmb_distsamp")
  fit <- fit_model(nll_fun, inputs = inputs, submodels = submodels,
                   starts = starts, method = method, se = se, ...)

  new("unmarkedFitDS", fitType = "distsamp", call = match.call(),
      opt = fit$opt, formula = formula, data = data, keyfun=keyfun,
      sitesRemoved = removed_sites(response), unitsOut=unitsOut,
      estimates = fit$submodels, AIC = fit$AIC, negLogLike = fit$opt$value,
      nllFun = fit$nll, output=output, TMB=fit$opt$TMB)
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


# Detection functions
gxhn <- function(x, sigma) exp(-x^2/(2 * sigma^2))
gxexp <- function(x, rate) exp(-x / rate)
gxhaz <- function(x, shape, scale)  1 - exp(-(x/shape)^-scale)
grhn <- function(r, sigma) exp(-r^2/(2 * sigma^2)) * r
grexp <- function(r, rate) exp(-r / rate) * r
grhaz <- function(r, shape, scale)  (1 - exp(-(r/shape)^-scale)) * r

dxhn <- function(x, sigma)
	gxhn(x=x, sigma=sigma) / integrate(gxhn, 0, Inf, sigma=sigma)$value
drhn <- function(r, sigma)
	grhn(r=r, sigma=sigma) / integrate(grhn, 0, Inf, sigma=sigma)$value
dxexp <- function(x, rate)
	gxexp(x=x, rate=rate) / integrate(gxexp, 0, Inf, rate=rate)$value
drexp <- function(r, rate)
	grexp(r=r, rate=rate) / integrate(grexp, 0, Inf, rate=rate)$value
dxhaz <- function(x, shape, scale)
	gxhaz(x=x, shape=shape, scale=scale) / integrate(gxhaz, 0, Inf,
		shape=shape, scale=scale)$value
drhaz <- function(r, shape, scale)
	grhaz(r=r, shape=shape, scale=scale) / integrate(grhaz, 0, Inf,
		shape=shape, scale=scale)$value


