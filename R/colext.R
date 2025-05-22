colext <- function(psiformula = ~ 1, gammaformula = ~ 1,
                    epsilonformula = ~ 1, pformula = ~ 1,
                    data, starts = NULL, method = "BFGS", se = TRUE, ...){

  # Check inputs
  formula <- list(psiformula = psiformula, gammaformula = gammaformula,
                  epsilonformula = epsilonformula, pformula = pformula)
  check_no_support(formula)

  # Build submodels
  submodels <- unmarkedSubmodelList(    
    psi = unmarkedSubmodelState(name = "Initial", short_name = "psi", 
                                  type = "psi", formula = psiformula, data = data, 
                                  family = "binomial", link = "logit"),

    col = unmarkedSubmodelTransition(name = "Colonization", short_name = "col", 
                                     type = "col", formula = gammaformula, 
                                     data = data, family = "binomial", link = "logit"),

    ext = unmarkedSubmodelTransition(name = "Extinction", short_name = "ext", 
                                     type = "ext", formula = epsilonformula, 
                                     data = data, family = "binomial", link = "logit"),

    det = unmarkedSubmodelDet(name = "Detection", short_name = "p", 
                              type = "det", formula = pformula, data = data, 
                              family = "binomial", link = "logit")
  )

  # Build response object
  response <- unmarkedResponseBinary(data, submodels)

  # Get nll inputs
  inputs <- nll_inputs(response, submodels, engine="TMB")

  # Adjust Kmin to be by period
  # TODO: do this automatically
  inputs$Kmin <- Kmin_by_T(inputs$y, data@numPrimary)

  # Determine which periods were sampled at each site
  M <- nrow(inputs$X_psi)
  T <- inputs$T_col
  J <- inputs$J_col
  inputs$site_sampled <- matrix(1, M, T)
  Trep <- rep(1:T, each = J)
  for (t in 1:T){
    ind <- which(Trep == t)
    for (i in 1:M){
      if(all(is.na(inputs$y[i, ind]))) inputs$site_sampled[i,t] <- 0
    }
  }

  # Fit model
  fit <- fit_model("tmb_colext", inputs = inputs, submodels = submodels,
                   starts = starts, method = method, se = se, ...)

  # Compute projected estimates
  M <- nrow(inputs$X_psi)
  T <- data@numPrimary
  J <- inputs$J_col
  psis <- plogis(inputs$X_psi %*% coef(fit$submodels["psi"]))
  phis <- array(NA,c(2,2,T-1,M))
  phis[,1,,] <- plogis(inputs$X_col %x% c(-1,1) %*% coef(fit$submodels["col"]))
  phis[,2,,] <- plogis(inputs$X_ext %x% c(-1,1) %*% -coef(fit$submodels["ext"]))

  projected <- array(NA, c(2, T, M))
  projected[1,1,] <- 1 - psis
  projected[2,1,] <- psis
  for(i in 1:M) {
    for(t in 2:T) {
      projected[,t,i] <- phis[,,t-1,i] %*% projected[,t-1,i]
    }
  }
  projected.mean <- apply(projected, 1:2, mean)
  rownames(projected.mean) <- c("unoccupied","occupied")
  colnames(projected.mean) <- 1:T

  # Compute smoothed estimates
  smoothed <- calculate_smooth(y = inputs$y, psi = psis,
                               gam = plogis(inputs$X_col %*% coef(fit$submodels["col"])),
                               eps = plogis(inputs$X_ext %*% coef(fit$submodels["ext"])),
                               p = plogis(inputs$X_det %*% coef(fit$submodels["det"])),
                               M = M, T = T, J = J)
  smoothed.mean <- apply(smoothed, 1:2, mean)
  rownames(smoothed.mean) <- c("unoccupied","occupied")
  colnames(smoothed.mean) <- 1:T

  new("unmarkedFitColExt", fitType = "colext", call = match.call(),
      formula = as.formula(paste(unlist(formula),collapse=" ")),
      psiformula = psiformula, gamformula = gammaformula, 
      epsformula = epsilonformula, detformula = pformula,
      data = data, sitesRemoved = removed_sites(response),
      estimates = fit$submodels, AIC = fit$AIC, opt = fit$opt,
      negLogLike = fit$opt$value, nllFun = fit$nll,
      projected = projected, projected.mean = projected.mean,
      smoothed = smoothed, smoothed.mean = smoothed.mean)
}

# Based on Weir, Fiske, Royle 2009 "TRENDS IN ANURAN OCCUPANCY"
# Appendix 1
calculate_smooth <- function(y, psi, gam, eps, p, M, T, J){

  smoothed <- array(NA, c(2, T, M))

  # Turn parameters into matrices
  p <- matrix(p, M, T*J, byrow=TRUE)
  gam <- matrix(gam, M, (T-1), byrow=TRUE)
  eps <- matrix(eps, M, (T-1), byrow = TRUE)

  tind <- rep(1:T, each = J)

  for (i in 1:M){

    # Forward pass
    alpha1 <- matrix(NA, M, T)
    alpha0 <- matrix(NA, M, T)

    # Initialize at t=1
    ysub <- y[i, tind == 1]
    #no_detects <- all(na.omit(ysub) == 0) * 1 
    psub <- p[i, tind == 1]

    if(all(is.na(ysub))){
      # Don't include detection likelihood in calculation if no data
      alpha1[i,1] <- psi[i]
      alpha0[i,1] <- (1-psi[i])
    } else {
      # Case when z = 1
      cp <- prod(na.omit(dbinom(ysub, 1, psub)))
      alpha1[i,1] <- psi[i] * cp 

      # Case when z = 0
      cp <- prod(na.omit(dbinom(ysub, 1, 0)))
      alpha0[i,1] <- (1-psi[i]) * cp
    }

    for (t in 2:T){
      ysub <- y[i, tind == t]
      psub <- p[i, tind == t]

      if(all(is.na(ysub))){
        alpha1[i,t] <- alpha0[i,t-1] * gam[i,t-1] + alpha1[i,t-1] * (1 - eps[i,t-1])
        alpha0[i,t] <- alpha0[i,t-1] * (1-gam[i,t-1]) + alpha1[i,t-1] * eps[i,t-1]
      } else {
        # Case when z = 1
        cp <- prod(na.omit(dbinom(ysub, 1, psub)))
        alpha1[i,t] <- (alpha0[i,t-1] * gam[i,t-1] + alpha1[i,t-1] * (1 - eps[i,t-1])) * cp  

        # Case when z = 0
        cp <- prod(na.omit(dbinom(ysub, 1, 0)))
        alpha0[i,t] <- (alpha0[i,t-1] * (1-gam[i,t-1]) + alpha1[i,t-1] * eps[i,t-1]) * cp
      }
    }

    # Backwards pass
    beta1 <- matrix(NA, M, T)
    beta0 <- matrix(NA, M, T)

    # Initialize
    beta1[i, T] <- 1
    beta0[i, T] <- 1
  
    for (t in (T-1):1){
      ysub <- y[i, tind == t+1]
      psub <- p[i, tind == t+1]

      if(all(is.na(ysub))){
        beta1[i, t] <- eps[i,t] * beta0[i, t+1] + (1-eps[i,t]) * beta1[i, t+1]
        beta0[i, t] <- (1-gam[i,t]) * beta0[i, t+1] + gam[i,t] * beta1[i, t+1]
      } else {
        cp1 <- prod(na.omit(dbinom(ysub, 1, psub)))
        cp0 <- prod(na.omit(dbinom(ysub, 1, 0)))

        # Case when z = 1
        beta1[i, t] <- eps[i,t] * cp0 * beta0[i, t+1] + (1-eps[i,t]) * cp1 * beta1[i, t+1]
        # Case when z = 0
        beta0[i, t] <- (1-gam[i,t]) * cp0 * beta0[i, t+1] + gam[i,t] * cp1 * beta1[i, t+1]
      }
    }

    out <- rep(0, T)
    for (t in 1:T){
      out[t] <- (alpha1[i,t] * beta1[i,t]) / (alpha0[i,t]*beta0[i,t] + alpha1[i,t] * beta1[i,t])
    }

    smoothed[1, 1:T, i] <- 1 - out
    smoothed[2, 1:T, i] <- out
  }

  smoothed
}
