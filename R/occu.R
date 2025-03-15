#  Fit the occupancy model of MacKenzie et al (2002).

setClass("unmarkedFit2Occu", contains = "unmarkedFit2")

occu <- function(formula, data, knownOcc = numeric(0),
                  linkPsi = c("logit", "cloglog"), starts = NULL, method = "BFGS",
                  se = TRUE, engine = c("TMB", "R", "C"), ...) {

  # Check arguments------------------------------------------------------------
  if(!is(data, "unmarkedFrameOccu"))
    stop("Data is not an unmarkedFrameOccu object.")

  engine <- match.arg(engine)
  if(engine == "C") engine <- "TMB"
  forms <- split_formula(formula)
  if(engine == "R" & any(sapply(forms, has_random))){
    message("Switching to TMB engine to fit random effects")
    engine <- "TMB"
  }

  known_occ <- rep(0, numSites(data))
  known_occ[knownOcc] <- 1

  linkPsi <- match.arg(linkPsi)

  submods <- unmarkedSubmodelList(
    state = unmarkedSubmodelState(long_name = "Occupancy", short_name = "psi", 
                                  type = "state", formula = forms[[2]], data = data, 
                                  family = "binomial", link = linkPsi,
                                  auxiliary = list(known_occ = known_occ)),

    det = unmarkedSubmodelDet(long_name = "Detection", short_name = "p", 
                              type = "det", formula = forms[[1]], data = data, 
                              family = "binomial", link = "logit")
  )

  resp <- unmarkedResponse(data)

  comps <- unmarkedModelComponents(resp, submods, auxiliary = list())

  fit <- fit_TMB2("tmb_occu", components = comps, starts = starts, se = se, 
                  method = method, ...)
  
  out <- new("unmarkedFit2Occu", call = match.call(),  data = data,
             components = fit$components, opt = fit$opt, TMB = fit$TMB)

  out
}

knownOcc <- function(object){
  object['state']@auxiliary$known_occ
}

setMethod("simulate_internal", "unmarkedFit2Occu",
  function(object, nsim){
  p <- getP(object, na.rm = FALSE)
  M <- nrow(p)
  J <- ncol(p)
  p <- as.vector(t(p))
  psi <- predict(object, type = "state", level = NULL)$Predicted
  
  lapply(1:nsim, function(i){
    z <- stats::rbinom(M, 1, psi)
    z[knownOcc(object) == 1] <- 1
    z <- rep(z, each = J)
    yvec <- z * stats::rbinom(M*J, 1, p)
    matrix(yvec, M, J, byrow = TRUE)
  })
})

setMethod("ranef_internal", "unmarkedFit2Occu", function(object, ...){
    psi <- predict(object, type="state", level=NULL)$Predicted
    R <- length(psi)
    p <- getP(object)
    z <- 0:1
    y <- getY(getData(object))
    y[y>1] <- 1
    post <- array(NA, c(R,2,1))
    colnames(post) <- z
    for(i in 1:R) {
        if(all(is.na(y[i,]))) next
        f <- dbinom(z, 1, psi[i])
        g <- rep(1, 2)
        for(j in 1:ncol(y)) {
            if(is.na(y[i,j]) | is.na(p[i,j]))
                next
            g <- g * dbinom(y[i,j], 1, z*p[i,j])
        }
        fudge <- f*g
        post[i,,1] <- fudge / sum(fudge)
    }
    new("unmarkedRanef", post=post)
})






occu_old <- function(formula, data, knownOcc = numeric(0),
                 linkPsi = c("logit", "cloglog"), starts, method = "BFGS",
                 se = TRUE, engine = c("C", "R", "TMB"), threads=1, ...) {

  # Check arguments------------------------------------------------------------
  if(!is(data, "unmarkedFrameOccu"))
    stop("Data is not an unmarkedFrameOccu object.")

  engine <- match.arg(engine, c("C", "R", "TMB"))
  if(any(sapply(split_formula(formula), has_random))) engine <- "TMB"
  if(length(knownOcc)>0 & engine == "TMB"){
    stop("TMB engine does not support knownOcc argument", call.=FALSE)
  }

  linkPsi <- match.arg(linkPsi, c("logit","cloglog"))
  psiLinkFunc <- ifelse(linkPsi=="cloglog", cloglog, plogis)
  psiInvLink <- ifelse(linkPsi=="cloglog", "cloglog", "logistic")
  psiLinkGrad <- ifelse(linkPsi=="cloglog", "cloglog.grad", "logistic.grad")

  # Format input data----------------------------------------------------------
  designMats <- getDesign(data, formula)
  y <- truncateToBinary(designMats$y)
  X <- designMats$X; V <- designMats$V;
  Z_state <- designMats$Z_state; Z_det <- designMats$Z_det
  removed <- designMats$removed.sites
  X.offset <- designMats$X.offset; V.offset <- designMats$V.offset

  # Re-format some variables for C and R engines
  yvec <- as.numeric(t(y))
  navec <- is.na(yvec)
  nd <- ifelse(rowSums(y,na.rm=TRUE) == 0, 1, 0) # no det at site i

  # convert knownOcc to logical so we can correctly to handle NAs.
  knownOccLog <- rep(FALSE, numSites(data))
  knownOccLog[knownOcc] <- TRUE
  if(length(removed)>0) knownOccLog <- knownOccLog[-removed]

  # Set up parameter names and indices-----------------------------------------
  occParms <- colnames(X)
  detParms <- colnames(V)
  nDP <- ncol(V)
  nOP <- ncol(X)
  nP <- nDP + nOP
  psiIdx <- 1:nOP
  pIdx <- (nOP+1):nP

  # Set up negative log likelihood functions for C++ and R engines-----_-------
  if(identical(engine, "C")) {
    nll <- function(params) {
      beta.psi <- params[1:nOP]
      beta.p <- params[(nOP+1):nP]
      nll_occu(
        yvec, X, V, beta.psi, beta.p, nd, knownOccLog, navec,
        X.offset, V.offset, linkPsi
      )
    }
  } else if (identical(engine, "R")){

    J <- ncol(y)
    M <- nrow(y)

    nll <- function(params) {
      psi <- psiLinkFunc(X %*% params[1 : nOP] + X.offset)
      psi[knownOccLog] <- 1
      pvec <- plogis(V %*% params[(nOP + 1) : nP] + V.offset)
      cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
      cp[navec] <- 1 # so that NA's don't modify likelihood
      cpmat <- matrix(cp, M, J, byrow = TRUE) #
      loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
      -sum(loglik)
    }
  }

  # Fit model with C++ and R engines-------------------------------------------
  if(engine %in% c("C", "R")){
    if(missing(starts)) starts <- rep(0, nP)
    if(length(starts) != nP){
      stop(paste("The number of starting values should be", nP))
    }
    fm <- optim(starts, nll, method = method, hessian = se, ...)
    covMat <- invertHessian(fm, nP, se)
    ests <- fm$par
    tmb_mod <- NULL
    fmAIC <- 2 * fm$value + 2 * nP #+ 2*nP*(nP + 1)/(M - nP - 1)

    # Organize effect estimates
    names(ests) <- c(occParms, detParms)
    state_coef <- list(ests=ests[1:nOP], cov=as.matrix(covMat[1:nOP,1:nOP]))
    det_coef <- list(ests=ests[(nOP+1):nP],
                     cov=as.matrix(covMat[(nOP+1):nP, (nOP+1):nP]))

    # No random effects for C++ and R engines
    state_rand_info <- det_rand_info <- list()

  # Fit model with TMB engine--------------------------------------------------
  } else if(identical(engine, "TMB")){

    # Set up TMB input data
    forms <- split_formula(formula)
    obs_all <- add_covariates(obsCovs(data), siteCovs(data), length(getY(data)))
    inps <- get_ranef_inputs(forms, list(det=obs_all, state=siteCovs(data)),
                             list(V, X), designMats[c("Z_det","Z_state")])

    tmb_dat <- c(list(y=y, no_detect=nd, link=ifelse(linkPsi=="cloglog",1,0),
                      offset_state=X.offset, offset_det=V.offset), inps$data)

    # Fit model with TMB
    if(missing(starts)) starts <- NULL
    tmb_out <- fit_TMB("tmb_occu", tmb_dat, inps$pars, inps$rand_ef,
                        starts=starts, method, ...)
    tmb_mod <- tmb_out$TMB
    fm <- tmb_out$opt
    fmAIC <- tmb_out$AIC
    nll <- tmb_mod$fn

    # Organize fixed effect estimates
    state_coef <- get_coef_info(tmb_out$sdr, "state", occParms, 1:nOP)
    det_coef <- get_coef_info(tmb_out$sdr, "det", detParms, (nOP+1):nP)

    # Organize random effect estimates
    state_rand_info <- get_randvar_info(tmb_out$sdr, "state", forms[[2]],
                                        siteCovs(data))
    det_rand_info <- get_randvar_info(tmb_out$sdr, "det", forms[[1]],
                                      obs_all)

    }

  # Create unmarkedEstimates---------------------------------------------------
  state <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                            estimates = state_coef$est,
                            covMat = state_coef$cov,
                            fixed = 1:nOP,
                            invlink = psiInvLink,
                            invlinkGrad = psiLinkGrad,
                            randomVarInfo=state_rand_info
                            )

  det <- unmarkedEstimate(name = "Detection", short.name = "p",
                          estimates =det_coef$est,
                          covMat = det_coef$cov,
                          fixed = 1:nDP,
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad",
                          randomVarInfo=det_rand_info
                          )

  estimateList <- unmarkedEstimateList(list(state=state, det=det))

  # Create unmarkedFit object--------------------------------------------------
  umfit <- new("unmarkedFitOccu", fitType = "occu", call = match.call(),
                 formula = formula, data = data,
                 sitesRemoved = designMats$removed.sites,
                 estimates = estimateList, AIC = fmAIC, opt = fm,
                 negLogLike = fm$value,
                 nllFun = nll, knownOcc = knownOccLog, TMB=tmb_mod)

  return(umfit)
}
