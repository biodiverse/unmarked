colext2 <- function(psiformula = ~ 1, gammaformula = ~ 1,
                    epsilonformula = ~ 1, pformula = ~ 1,
                    data, starts, method = "BFGS", se = TRUE, ...){

  ## truncate to 1
  data@y <- truncateToBinary(data@y)
  M <- numSites(data)
  T <- data@numPrimary
  J <- ncol(data@y) / T
    
  formula <- list(psiformula = psiformula, gammaformula = gammaformula,
                  epsilonformula = epsilonformula, pformula = pformula)
  check_no_support(formula)
  designMats <- getDesign(data, formula = as.formula(paste(unlist(formula), collapse=" ")))
  X_psi <- designMats$W
  X_gam <- designMats$X.gam
  X_eps <- designMats$X.eps
  X_det <- designMats$V
  y <- designMats$y
  M <- nrow(y)
  T <- ncol(y) / J
  psiParms <- colnames(X_psi)
  gamParms <- colnames(X_gam)
  epsParms <- colnames(X_eps)
  detParms <- colnames(X_det)

  ## remove final year from transition prob design matrices
  X_gam <- as.matrix(X_gam[-seq(T,M*T,by=T),])
  X_eps <- as.matrix(X_eps[-seq(T,M*T,by=T),])

  # Determine which periods were sampled at each site
  pers_sampled <- matrix(1, M, T)
  Trep <- rep(1:T, each = J)
  for (t in 1:T){
    ind <- which(Trep == t)
    for (i in 1:M){
      if(all(is.na(y[i, ind]))) pers_sampled[i,t] <- 0
    }
  }
  T_site <- apply(pers_sampled, 1, sum)

  Tsamp <- matrix(NA, M, T)
  for (i in 1:M){
    which_sampled <- which(pers_sampled[i,] == 1)
    Tsamp[i,1:length(which_sampled)] <- which_sampled
  }

  no_detects <- matrix(1, M, T)
  for (i in 1:M){
    for (t in 1:T){
      ind <- which(Trep == t)
      ysub <- y[i, ind]
      if(all(is.na(ysub))) next
      ysub <- na.omit(ysub)
      if(all(ysub == 0)){
        no_detects[i,t] <- 1
      } else {
        no_detects[i,t] <- 0
      }
    }
  }

  nDP <- length(detParms)
  nGP <- length(gamParms)
  nEP <- length(epsParms)
  nSP <- length(psiParms)
  pind_mat <- matrix(0, 4, 2)
  pind_mat[1,] <- c(1, length(psiParms))
  pind_mat[2,] <- max(pind_mat[1,]) + c(1, length(gamParms))
  pind_mat[3,] <- max(pind_mat[2,]) + c(1, length(epsParms))
  pind_mat[4,] <- max(pind_mat[3,]) + c(1, length(detParms))
  #nDMP <-  1
  
  tmb_dat <- list(y = as.vector(t(y)), X_psi = X_psi, X_gam = X_gam, X_eps = X_eps, X_det = X_det,
                  M = M, T = T, J = J, site_sampled = pers_sampled, nd = no_detects)

  tmb_pars <- list(beta_psi = rep(0, nSP), beta_gam = rep(0, nGP),
                   beta_eps = rep(0, nEP), beta_det = rep(0, nDP))



  #library(TMB)  
  #compile("~/downloads/colext.cpp")
  #dyn.load(dynlib("~/downloads/colext"))

  #tmb_obj <- MakeADFun(data = tmb_dat, parameters=tmb_pars,
  #                     DLL = "colext", silent=TRUE)

  tmb_obj <- TMB::MakeADFun(data = c(model = "tmb_colext", tmb_dat), 
                            parameters = tmb_pars,
                            DLL = "unmarked_TMBExports", silent=TRUE)
  
  opt <- optim(unlist(tmb_pars), fn=tmb_obj$fn, gr=tmb_obj$gr, 
               method=method, hessian = se)#, ...)

  fmAIC <- 2 * opt$value + 2 * length(unlist(tmb_pars))

  sdr <- TMB::sdreport(tmb_obj)

  psi_coef <- get_coef_info(sdr, "psi", psiParms, pind_mat[1,1]:pind_mat[1,2])
  gam_coef <- get_coef_info(sdr, "gam", gamParms, pind_mat[2,1]:pind_mat[2,2])
  eps_coef <- get_coef_info(sdr, "eps", epsParms, pind_mat[3,1]:pind_mat[3,2])
  det_coef <- get_coef_info(sdr, "det", detParms, pind_mat[4,1]:pind_mat[4,2])
  
  psi <- unmarkedEstimate(name = "Initial", short.name = "psi",
                          estimates = psi_coef$ests,
                          covMat = psi_coef$cov,
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

  col <- unmarkedEstimate(name = "Colonization", short.name = "col",
                          estimates = gam_coef$ests,
                          covMat = gam_coef$cov,
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

  ext <- unmarkedEstimate(name = "Extinction", short.name = "ext",
                          estimates = eps_coef$ests,
                          covMat = eps_coef$cov,
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

  det <- unmarkedEstimate(name = "Detection", short.name = "p",
                          estimates = det_coef$ests,
                          covMat = det_coef$cov,
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

  estimateList <- unmarkedEstimateList(list(psi = psi, col = col,
                                            ext = ext, det=det))
  
  psis <- plogis(X_psi %*% psi_coef$ests)

  ## computed projected estimates
  phis <- array(NA,c(2,2,T-1,M))
  phis[,1,,] <- plogis(X_gam %x% c(-1,1) %*% gam_coef$ests)
  phis[,2,,] <- plogis(X_eps %x% c(-1,1) %*% -eps_coef$ests)

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

  ## smoothing (not supported right now)
  smoothed <- calculate_smooth(y = y, psi = psis,
                               gam = plogis(X_gam %*% gam_coef$ests),
                               eps = plogis(X_eps %*% eps_coef$ests),
                               p = plogis(X_det %*% det_coef$ests),
                               M = M, T = T, J = J)
  smoothed.mean <- apply(smoothed, 1:2, mean)
  rownames(smoothed.mean) <- c("unoccupied","occupied")
  colnames(smoothed.mean) <- 1:T

  #smoothed <- array(NA, c(2, T, M))
  #smoothed.mean <- projected.mean
  #smoothed.mean[] <- NA
  #forward(detParams, phis, psis, storeAlpha = TRUE)
  #  beta <- backward(detParams, phis)
  #  gamma <- array(NA, c(K + 1, nY, M))
  #  for(i in 1:M) {
  #  for(t in 1:nY) {
  #      gamma[,t,i] <- alpha[,t,i]*beta[,t,i] / sum(alpha[,t,i]*beta[,t,i])
  #  }}
  #  smoothed.mean <- apply(gamma, 1:2, mean)
  #  rownames(smoothed.mean) <- c("unoccupied","occupied")
  #  colnames(smoothed.mean) <- 1:nY

  umfit <- new("unmarkedFitColExt", fitType = "colext",
                call = match.call(),
                formula = as.formula(paste(unlist(formula),collapse=" ")),
                psiformula = psiformula,
                gamformula = gammaformula,
                epsformula = epsilonformula,
                detformula = pformula,
                data = data, sitesRemoved = designMats$removed.sites,
                estimates = estimateList,
                AIC = fmAIC, opt = opt, negLogLike = opt$value,
                nllFun = tmb_obj$fn,
                projected = projected,
                projected.mean = projected.mean,
                smoothed = smoothed, smoothed.mean = smoothed.mean)

  return(umfit)
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
