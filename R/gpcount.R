
# data will need to be an unmarkedMultFrame
gpcount <- function(lambdaformula, phiformula, pformula, data,
    mixture=c('P', 'NB', 'ZIP'), K = NULL, starts = NULL, method = "BFGS", 
    se = TRUE, engine=c('C', 'R'), threads = 1, ...){

  if(!is(data, "unmarkedFrameGPC")){
    stop("Data is not of class unmarkedFrameGPC.")
  }
  if(data@numPrimary < 2) stop("Use pcount instead", call.=FALSE)
  mixture <- match.arg(mixture)
  engine <- match.arg(engine)
  formlist <- list(lambdaformula = lambdaformula, phiformula = phiformula,
                   pformula = pformula)
  check_no_support(formlist)
  form <- as.formula(paste(unlist(formlist), collapse=" "))

  # Build submodels
  fam <- switch(mixture, P = "poisson", NB = "negative_binomial", ZIP = "ZIP")
  submodels <- unmarkedSubmodelList(
    lambda = unmarkedSubmodelState(name = "Abundance", short_name = "lam", 
                                  type = "lambda", formula = lambdaformula,
                                  data = data, family = fam, link = "log"),

    phi = unmarkedSubmodelAvail(name = "Availability", short_name = "phi", 
                                type = "phi", formula = phiformula, data = data, 
                                family = "binomial", link = "logit"),

    det = unmarkedSubmodelDet(name = "Detection", short_name = "p", 
                              type = "det", formula = pformula, data = data, 
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
  nll_fun <- switch(engine, R = nll_gpcount_R, C = nll_gpcount_Cpp)

  # Fit model
  fit <- fit_model(nll_fun, inputs = inputs, submodels = submodels,
                   starts = starts, method = method, se = se, ...)

  new("unmarkedFitGPC", fitType = "gpcount",
      call = match.call(), formula = form, formlist = formlist,
      data = data, estimates = fit$submodels, 
      sitesRemoved = removed_sites(response), AIC = fit$AIC, opt = fit$opt,
      negLogLike = fit$opt$value, nllFun = fit$nll, mixture=mixture, 
      K=response@Kmax)
}


nll_gpcount_R <- function(params, inputs){

  with(inputs, {

  M <- N <- 0:Kmax
  lM <- length(M)

  I  <- nrow(y)
  T <- T_phi
  J <- J_phi
  y <- array(y, c(I, J, T))

  beta_lambda <- params[idx_lambda[1]:idx_lambda[2]]
  beta_phi <- params[idx_phi[1]:idx_phi[2]]
  beta_det <- params[idx_det[1]:idx_det[2]]

  lam <- exp(X_lambda %*% beta_lambda + offset_lambda)
  phi <- plogis(X_phi %*% beta_phi + offset_phi)
  phi <- matrix(phi, I, T, byrow=TRUE)
  p <- plogis(X_det %*% beta_det + offset_det)
  p <- matrix(p, I, byrow=TRUE)
  p <- array(p, c(I, J, T))  # byrow?

  par2 <- 0
  if(family_lambda == 2){ # negbin
    par2 <- exp(params[idx_alpha[1]])
  } else if(family_lambda == 3){
    par2 <- plogis(params[idx_psi[1]])
  }

  L <- rep(NA, I)

    for(i in 1:I) {
        if(all(is.na(y[i,,]))) next
        f <- switch(family_lambda,
            "1" = dpois(M, lam[i], log=TRUE),
            "2" = dnbinom(M, mu=lam[i], size=par2, log=TRUE),
            "3" = log(dzip(M, lambda=lam[i], psi=par2))
        )
        ghi <- rep(0, lM)
        for(t in 1:T) {
          if(is.na(phi[i,t])) next
            gh <- matrix(-Inf, lM, lM)
            for(m in M) {
                if(m < max(y[i,,], na.rm=TRUE)) {
                    gh[,m+1] <- -Inf
                    next
                }
                # Changing this 5/21/25, which breaks backward compatibility
                # when phi[i,t] is missing. It should be equivalent to 
                # simply having all ys for this T missing (see ~L106)
                #if(is.na(phi[i,t])) {
                #    g <- rep(0, lM)
                #    g[N>m] <- -Inf
                #}
                #else
                g <- dbinom(N, m, phi[i,t], log=TRUE)
                h <- rep(0, lM)
                for(j in 1:J) {
                    if(is.na(y[i,j,t]))
                        next
                    h <- h + dbinom(y[i,j,t], N, p[i,j,t], log=TRUE)
                }
                gh[,m+1] <- g + h
            }
            ghi <- ghi + log(colSums(exp(gh))) # sum over N(t)
        }
        fgh <- f + ghi
        L[i] <- sum(exp(fgh)) # sum over M
    }
    return(-sum(log(L), na.rm = TRUE))
  })
}
