
#  Fit the occupancy model of MacKenzie et al (2002).

occu <- function(formula, data, knownOcc = numeric(0),
                 linkPsi = c("logit", "cloglog"), starts = NULL, method = "BFGS",
                 se = TRUE, engine = c("C", "R", "TMB"),
                 lambda = NULL, pen.type = c("Bayes", "Ridge", "MPLE"), ...){

  # Check arguments
  if(!is(data, "unmarkedFrameOccu")){
    stop("Data is not an unmarkedFrameOccu object.")
  }
  engine <- match.arg(engine)
  forms <- split_formula(formula)
  if(any(sapply(forms, has_random))) engine <- "TMB"
  linkPsi <- match.arg(linkPsi)
  pen.type <- match.arg(pen.type)

  # Build submodels
  known_occ <- rep(0, numSites(data))
  known_occ[knownOcc] <- 1
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
  response <- unmarkedResponseBinary(data, submodels)

  # Get nll inputs
  inputs <- nll_inputs(response, submodels, engine)
  nll_fun <- switch(engine, R = nll_occu_R, C = nll_occu_Cpp, TMB = "tmb_occu")

  # Handle penalized likelihood
  inputs$pen_type <- 0 # Skip penalty by default
  use_penalty <- !is.null(lambda)
  if(use_penalty){
    check_occu_penalty_args(formula, data, knownOcc, linkPsi, starts,
                            method, engine, lambda, pen.type)
    se <- FALSE
    pen_info <- penalty_info(submodels, response, pen.type, lambda, starts)
    inputs <- utils::modifyList(inputs, pen_info)
  }

  # Fit model
  fit <- fit_model(nll_fun, inputs = inputs, submodels = submodels,
                   starts = starts, method = method, se = se, ...)

  # Create unmarkedFit object
  umfit <- new("unmarkedFitOccu", fitType = "occu", call = match.call(),
    formula = formula, data = data, sitesRemoved = removed_sites(response),
    estimates = fit$submodels, AIC = fit$AIC, opt = fit$opt,
    negLogLike = fit$opt$value, nllFun = fit$nll, 
    knownOcc = as.logical(known_occ), TMB=fit$TMB)

  # If penalized likelihood was used, convert to unmarkedFitOccuPEN
  if(use_penalty){
    umfit <- occu_to_occuPEN(umfit, lambda = lambda, pen.type = pen.type)
  }

  umfit
}

nll_occu_R <- function(params, inputs){

  penalty <- calculate_occu_penalty(params, inputs)
  
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
  -(sum(loglik, na.rm = TRUE) - penalty)
  
  })
}

check_occu_penalty_args <- function(formula, data, knownOcc, linkPsi, starts,
                                    method, engine, lambda, 
                                    pen.type = c("Bayes","Ridge","MPLE")){
 
  if(engine == "TMB"){
    stop("Penalized likelihood not supported by TMB engine", call.=FALSE)
  }
  if(pen.type == "MPLE"){
    mple <- computeMPLElambda(formula, data, knownOcc=knownOcc, starts = starts,
                              method = method, engine = engine)
    if(round(mple, 4) != round(lambda, 4)){
      warning("Supplied lambda does not match the computed value.", call.=FALSE)
    }
  }
  if(linkPsi == "cloglog"){
    stop("Penalties not supported with cloglog", call.=FALSE)
  }

  invisible()
}

penalty_info <- function(submodels, response, type, lambda, starts){
  if(is.null(lambda)) return(list(pen_type = 0))
  out <- list(pen_lambda = lambda,
              pen_type = switch(type, Bayes=1, Ridge=2, MPLE=3))

  if(type == "Ridge"){
    has_int <- lapply(submodels(submodels), function(x){
      attr(stats::terms(x@formula), "intercept")
    })
    has_cov <- sapply(submodels(submodels), function(x) ncol(model.matrix(x)) > 1)
    if(all(!has_cov)) stop("Ridge requires covariates", call.=FALSE)
    names(has_int) <- paste0("has_int_", names(submodels(submodels)))
    out <- c(out, has_int)
  }

  if(type == "MPLE"){
    X <- model.matrix(submodels["state"])
    if(ncol(X) == 1) stop("MPLE requires occupancy covariates", call.=FALSE)
    y <- response@Kmin
    starts <- if(is.null(starts)) default_starts(submodels["state"])
    fit <- glm.fit(x=X, y=y, family=binomial(), intercept=FALSE,
                       start=starts[1:ncol(X)])
    out$LR_est <- fit$coefficients
  }

  out
}

calculate_occu_penalty <- function(params, inputs){
  with(inputs, {
  
  if(pen_type == 0) return(0)

  beta_state <- params[idx_state[1]:idx_state[2]]
  beta_det <- params[idx_det[1]:idx_det[2]]

  if(pen_type == 1){        # Bayes
    penalty <- sum(params^2) * pen_lambda * 0.5
  } else if(pen_type == 2){ # Ridge
    pars_ridge <- c()
    if(has_int_state & length(beta_state) > 1){
      pars_ridge <- c(pars_ridge, beta_state[2:length(beta_state)])
    } else if(!has_int_state){
      pars_ridge <- c(pars_ridge, beta_state)
    }
    if(has_int_det & length(beta_det) > 1){
      pars_ridge <- c(pars_ridge, beta_det[2:length(beta_det)])
    } else if(!has_int_det){
      pars_ridge <- c(pars_ridge, beta_det)
    }
    penalty <- sum(pars_ridge^2) * pen_lambda * 0.5
  } else if(pen_type == 3){ # MPLE
    penalty <- sum(abs(beta_state - LR_est)) * pen_lambda
  }

  penalty
  })
}

computeMPLElambda = function(formula, data, knownOcc = numeric(0), 
                             starts = NULL, method = "BFGS", engine = c("C", "R")){
  
  if(split_formula(formula)[[2]] == ~1){
    stop("MPLE requires occupancy covariates", call.=FALSE)
  }
  engine <- match.arg(engine)
  mle <- occu(formula = formula, data = data, knownOcc = knownOcc,
              starts = starts, method = method, engine = engine)

  X_state <- model.matrix(mle["state"])
  X_det <- model.matrix(mle["det"])

  response <- unmarkedResponseBinary(data, mle@estimates)
  ymax <- response@Kmin

  starts <- if(is.null(starts)) starts <- default_starts(mle["state"])
  starts <- starts[1:ncol(X_state)]

  LR_fit <- stats::glm.fit(x=X_state, y=ymax, family=binomial(), intercept=FALSE,
                           start=starts)

  naive_occ <- mean(LR_fit$fitted.values)  
  mean_det <- mean((1 + exp(-coef(mle, "det") %*% t(X_det)))^-1)
  J <- obsNum(data)
  mple_lambda <- sqrt(sum(diag(vcov(mle, "det")))) *
                 (1- (1 - mean_det)^J) * (1 - naive_occ)

  mple_lambda
}

occu_to_occuPEN <- function(umfit, lambda, pen.type){
  new("unmarkedFitOccuPEN", fitType = "occu", call = umfit@call,
      formula = umfit@formula, data = umfit@data, sitesRemoved = umfit@sitesRemoved,
      estimates = umfit@estimates, AIC = umfit@AIC, opt = umfit@opt,
      negLogLike = umfit@negLogLike, nllFun = umfit@nllFun, 
      knownOcc = umfit@knownOcc, TMB=umfit@TMB, 
      lambda = lambda, pen.type = pen.type)
}
