# Fit the multinomial-Poisson abundance mixture model.

multinomPois <- function(formula, data, starts = NULL, method = "BFGS",
                         se = TRUE, engine = c("C","R","TMB"), ...){

  # Check arguments
  if(!is(data, "unmarkedFrameMPois")){
    stop("Data is not a data frame or unmarkedFrame.")
  }
  engine <- match.arg(engine, c("C", "R", "TMB"))
  forms <- split_formula(formula)
  if(any(sapply(forms, has_random))) engine <- "TMB"

  c_support <- c('doublePiFun','removalPiFun','depDoublePiFun')
  if(engine %in% c("C", "TMB") & !data@piFun %in% c_support){
    warning("Custom pi functions are not supported by C engine. Using R engine instead.",
            call.=FALSE)
    engine <- "R"
  }

  # Build submodels
  submodels <- unmarkedSubmodelList(
    state = unmarkedSubmodelState(name = "Abundance", short_name = "lambda", 
                                  type = "state", formula = forms[[2]], 
                                  data = data, family = "poisson", link = "log"),

    det = unmarkedSubmodelMultinom(name = "Detection", short_name = "p", 
                                   type = "det", formula = forms[[1]], 
                                   data = data, link = "logit")
  )

  # Build response object
  response <- unmarkedResponseCount(data, submodels, Kmax = NULL)

  # Generate engine inputs
  if(engine == "TMB"){
    inputs <- engine_inputs_TMB(response, submodels)
  } else {
    inputs <- engine_inputs_CR(response, submodels)
  }

  # Fit model
  nll_fun <- switch(engine, R = nll_multinomPois_R, C = nll_multinomPois_Cpp, 
                    TMB = "tmb_multinomPois")
  fit <- fit_model(nll_fun, inputs = inputs, submodels = submodels,
                   starts = starts, method = method, se = se)#, ...)
   
  # Create unmarkedFit object
  new("unmarkedFitMPois", fitType="multinomPois", call=match.call(),
      formula = formula, data = data, sitesRemoved = removed_sites(response),
      estimates = fit$submodels, AIC = fit$AIC, opt = fit$opt,
      negLogLike = fit$opt$value, nllFun = fit$nll, TMB=fit$TMB)
}


nll_multinomPois_R <- function(params, inputs){
  with(inputs, {
  beta_state <- params[idx_state[1]:idx_state[2]]
  beta_det <- params[idx_det[1]:idx_det[2]]

  M <- nrow(y)
  J <- ncol(y)

  lambda <- exp(X_state %*% beta_state + offset_state)
  lambda <- matrix(rep(lambda, J), nrow = M, ncol = J) 
  p <- plogis(X_det %*% beta_det + offset_det)
  p <- matrix(p, M, R_det, byrow = TRUE)
  pi <- do.call(pi_fun_det, list(p = p))

  loglik <- dpois(y, lambda * pi, log = TRUE)
   
  -sum(loglik, na.rm = TRUE)
  })
}
