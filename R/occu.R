
#  Fit the occupancy model of MacKenzie et al (2002).

occu <- function(formula, data, knownOcc = numeric(0),
                 linkPsi = c("logit", "cloglog"), starts = NULL, method = "BFGS",
                 se = TRUE, engine = c("C", "R", "TMB"), threads=1, ...) {

  # Check arguments
  if(!is(data, "unmarkedFrameOccu")) stop("Data is not an unmarkedFrameOccu object.")

  engine <- match.arg(engine)
  forms <- split_formula(formula)
  if(any(sapply(forms, has_random))) engine <- "TMB"

  linkPsi <- match.arg(linkPsi)

  known_occ <- rep(0, numSites(data))
  known_occ[knownOcc] <- 1

  # Build submodels
  submodels <- unmarkedSubmodelList(
    state = unmarkedSubmodelState(name = "Occupancy", short_name = "psi", 
                                  type = "state", formula = forms[[2]], data = data, 
                                  family = "binomial", link = linkPsi,
                                  auxiliary = list(known_occ = known_occ)),

    det = unmarkedSubmodelDet(name = "Detection", short_name = "p", 
                              type = "det", formula = forms[[1]], data = data, 
                              family = "binomial", link = "logit")
  )

  # Build response object
  response <- unmarkedResponse(data)
  # Handle missing values in covariates
  response <- add_missing(response, submodels)

  # Fit model
  nll_fun <- switch(engine, R = nll_occu_R, C = nll_occu_Cpp, TMB = "tmb_occu")
  fit <- fit_model(nll_fun, response = response, submodels = submodels,
                   starts = starts, method = method, se = se, ...)

  # Create unmarkedFit object
  new("unmarkedFitOccu", fitType = "occu", call = match.call(),
      formula = formula, data = data, sitesRemoved = removed_sites(response),
      estimates = fit$submodels, AIC = fit$AIC, opt = fit$opt,
      negLogLike = fit$opt$value, nllFun = fit$nll, 
      knownOcc = as.logical(known_occ), TMB=fit$TMB)
}

nll_occu_R <- function(params, inputs){
  with(inputs, {
  
  beta_state <- params[idx_state[1]:idx_state[2]]
  beta_det <- params[idx_det[1]:idx_det[2]]

  J <- ncol(inputs$y)
  M <- nrow(inputs$y)
  y <- as.vector(t(inputs$y))
  nd <- 1 - Kmin

  invlink <- ifelse(invlink_state == 0, plogis, cloglog)
  psi <- invlink(X_state %*% beta_state + offset_state)
  psi[known_occ_state == 1] <- 1
  pvec <- plogis(X_det %*% beta_det + offset_det)
  cp <- (pvec^y) * ((1 - pvec)^(1 - y))
  cpmat <- matrix(cp, M, J, byrow = TRUE)
  loglik <- log(rowProds(cpmat, na.rm = TRUE) * psi + nd * (1 - psi))
  -sum(loglik, na.rm = TRUE)
  
  })
}
