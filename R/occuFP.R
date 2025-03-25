#  Fit occupancy data where false positives oc cur
#  three data types are allowed
#  1) standard occupancy data as in MacKenzie et al (2002).
#  2) Royle-Link (2006) type
#  3) certain-uncertain data like Miller (2011)

occuFP <- function(detformula = ~ 1, FPformula = ~ 1, Bformula = ~ 1, 
                   stateformula = ~ 1, data, starts = NULL, method = "BFGS", 
                   se = TRUE, engine = "R", ...) {

  # Check arguments
  if(!is(data, "unmarkedFrameOccuFP")){
    stop("Data is not an unmarkedFrameOccuFP object.")
  }
  # these checks should be in unmarkedFrameOccuFP
  if(sum(data@type[2:3])==0){
    stop("Only type 1 data. No data types with false positives. Use occu instead.")
  }
  if(any(data@type[1:2] > 0)){
    if(any(data@y[,1:sum(data@type[1:2])] > 1, na.rm = TRUE)){
      stop("Values of y for type 1 and type 2 data must be 0 or 1.", call.=FALSE)
    }
  }
  if(data@type[3] > 0){
    if(any(data@y > 2, na.rm = TRUE)){
      stop("Values of y for type 3 data must be 0, 1, or 2.", call.=FALSE)
    }
  }
  check_no_support(list(detformula,FPformula,Bformula,stateformula))
  engine <- match.arg(engine)

  # Build submodels
  submodels <- unmarkedSubmodelList(
    state = unmarkedSubmodelState(name = "Occupancy", short_name = "psi", 
                                  type = "state", formula = stateformula, data = data, 
                                  family = "binomial", link = "logit"),

    det = unmarkedSubmodelDet(name = "Detection", short_name = "p", 
                              type = "det", formula = detformula, data = data, 
                              family = "binomial", link = "logit"),

    fp = unmarkedSubmodelDet(name = "false positive", short_name = "fp",
                             type = "fp", formula = FPformula, data = data,
                             family = "binomial", link = "logit")
  )

  if(data@type[3] != 0){
     submodels["b"] = unmarkedSubmodelDet(name = "Pcertain", short_name = "b",
                                          type = "b", formula = Bformula, 
                                          data = data, family = "binomial", 
                                          link = "logit") 
  }

  response <- unmarkedResponse(data, submodels, Kmax = NULL,
                               auxiliary = list(data_type = data@type))

  # Get nll inputs
  inputs <- nll_inputs(response, submodels, engine)

  fit <- fit_model(nll_occuFP_R, inputs = inputs, submodels = submodels,
                   starts = starts, method = method, se = se, ...)

  new("unmarkedFitOccuFP", fitType = "occuFP", call = match.call(),
      detformula = detformula, FPformula = FPformula, Bformula = Bformula,
      stateformula = stateformula, formula = ~1, type = data@type, data = data,
      sitesRemoved = removed_sites(response), estimates = fit$submodels, 
      AIC = fit$AIC, opt = fit$opt, negLogLike = fit$opt$value, nllFun = fit$nll)
}

nll_occuFP_R <- function(params, inputs){

  with(inputs, {

  M <- nrow(y)
  J <- ncol(y)

  yvec0 <- as.numeric(t(y==0))
  yvec1 <- as.numeric(t(y==1))
  yvec2 <- as.numeric(t(y==2))

  beta_state <- params[idx_state[1]:idx_state[2]]
  beta_det <- params[idx_det[1]:idx_det[2]]
  beta_fp <- params[idx_fp[1]:idx_fp[2]]

  psi <- plogis(X_state %*% beta_state + offset_state)
  pvec <- plogis(X_det %*% beta_det + offset_det)
  fvec <- plogis(X_fp %*% beta_fp + offset_fp)
    
  if(data_type[1] != 0){
    fvec[rep(c(rep(TRUE, data_type[1]),rep(FALSE, sum(data_type[2:3]))),M)] <- 0
  }

  if(data_type[3]!=0){
    beta_b <- params[idx_b[1]:idx_b[2]]
    bvec <- plogis(X_b %*% beta_b + offset_b)
    if (data_type[1]!=0 | data_type[2]!=0){
      bvec[rep(c(rep(TRUE,sum(data_type[1:2])),rep(FALSE,data_type[3])),M)] <- 0
    }
  } else if (data_type[3]==0){
    bvec <- matrix(0,M*J,1)
  }
  cp0 <- (1-fvec)^(yvec0) * (fvec)^(yvec1) * (1-yvec2)
  cp1 <- (1 - pvec)^(yvec0)*(pvec*(1-bvec))^(yvec1)*(pvec*bvec)^(yvec2)

  na_vec <- is.na(t(y))
  cp0[na_vec] <- 1
  cp1[na_vec] <- 1

  cp0mat <- matrix(cp0, M, J, byrow = TRUE) #
  cp1mat <- matrix(cp1, M, J, byrow = TRUE) #
  loglik <- log(rowProds(cp0mat) *(1-psi) + rowProds(cp1mat) * psi)
  -sum(loglik, na.rm = TRUE)
  })
}
