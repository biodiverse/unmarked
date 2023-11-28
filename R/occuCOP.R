
# Prerequisites (while its not integrated within the package)
library(unmarked)
source("R/unmarkedFrame.R")
source("R/utils.R")
source("R/unmarkedEstimate.R")
source("R/unmarkedFitList.R")
source("R/predict.R")
source("R/getDesign.R")
source("R/mixedModelTools.R")
# sapply(paste0("./R/", list.files("./R/")), source)

library(lattice)
library(testthat)

# Fit the occupancy model COP
# (Counting Occurrences Process)

# Occupancy
# z_i ~ Bernoulli(psi)
# 
# with:
#   z_i = Occupancy state of site i
#       = 1 if the site i is occupied
#       = 0 else

# Detection
# N_ijs | z_i = 1 ~ Poisson(lambda*T_s)
# N_ijs | z_i = 0 ~ 0
# 
# with:
#   N_is = Number of detections of site i during session s
#   z_i = Occupancy state of site i
#   lambda = Detection rate
#   T_s = Duration of the session s

# CLASSES ----------------------------------------------------------------------

## unmarkedFrameCOP class ----
setClass(
  "unmarkedFrameCOP",
  representation(obsLength = "matrix"),
  contains = "unmarkedFrame",
  validity = function(object) {
    errors <- character(0)
    M <- nrow(object@y)
    J <- ncol(object@y)
    y_integers = sapply(object@y, check.integer, na.ignore = T)
    if (!all(y_integers)) {
      errors <- c(errors,
                  paste(
                    "Count detection should be integers. Non-integer values:",
                    paste(object@y[which(!y_integers)], collapse = ', ')
                  ))
    }
    if (!all(all(dim(object@obsLength) == dim(object@y)))){
      errors <- c( errors, paste(
          "obsLength should be a matrix of the same dimension as y, with M =", M,
          "rows (sites) and J =", J, "columns (sampling occasions)."
    ))}
    if (length(errors) == 0) TRUE
    else errors
  }
)

## unmarkedFitCOP class ----
setClass("unmarkedFitCOP",
         representation(removed_obs = "matrix",
                        formlist = "list",
                        nll = "optionalNumeric",
                        convergence="optionalNumeric"),
         contains = "unmarkedFit")


# METHODS ----------------------------------------------------------------------

## getDesign method ----
# Example for occu: https://github.com/rbchan/unmarked/blob/536e32ad7b2526f8fac1b315ed2e99accac0d50f/R/getDesign.R#L10-L97
# Example for GDR: https://github.com/rbchan/unmarked/blob/c82e63947d7df7dfc896066e51dbf63bda3babf4/R/gdistremoval.R#L177-L253
setMethod(
  "getDesign", "unmarkedFrameCOP",
  function(umf, formlist, na.rm = TRUE) {
    
    # Retrieve useful informations from umf
    M <- numSites(umf)
    J <- obsNum(umf)
    y <- getY(umf)
    obsLength <- umf@obsLength

    # Occupancy submodel -------------------------------------------------------

    # Retrieve the fixed-effects part of the formula
    psiformula <- lme4::nobars(as.formula(formlist$psiformula))
    psiVars <- all.vars(psiformula)
    
    # Retrieve the site covariates
    sc <- siteCovs(umf)
    if(is.null(sc)) {
      sc <- data.frame(.dummy = rep(0, M))
    }
    
    # Check for missing variables
    psiMissingVars <- psiVars[!(psiVars %in% names(sc))]
    if (length(psiMissingVars) > 0) {
      stop(paste(
        "Variable(s)",
        paste(psiMissingVars, collapse = ", "),
        "not found in siteCovs"
      ), call. = FALSE)
    }

    # State model matrix
    Xpsi <- model.matrix(
      psiformula,
      model.frame(psiformula, sc, na.action = NULL)
    )
    Zpsi <- get_Z(formlist$psiformula, sc)
    
    # Detection submodel -------------------------------------------------------
    
    # Retrieve the fixed-effects part of the formula
    lambdaformula <- lme4::nobars(as.formula(formlist$lambdaformula))
    lambdaVars <- all.vars(lambdaformula)
    
    # Retrieve the observation covariates
    oc <- obsCovs(umf)
    if(is.null(oc)) {
      oc <- data.frame(.dummy = rep(0, M*J))
    }
    
    # Check for missing variables
    lambdaMissingVars <- lambdaVars[!(lambdaVars %in% names(oc))]
    if (length(lambdaMissingVars) > 0) {
      stop(paste(
        "Variable(s)",
        paste(lambdaMissingVars, collapse = ", "),
        "not found in obsCovs"
      ), call. = FALSE)
    }
    
    # Detection model matrix
    Xlambda <- model.matrix(
      lambdaformula,
      model.frame(lambdaformula, oc, na.action = NULL)
    )
    Zlambda <- get_Z(formlist$lambdaformula, oc)
    
    # Missing data -------------------------------------------------------------
    
    # Missing count data
    missing_y <- is.na(y)
    
    # Missing site covariates
    # (TRUE if at least one site covariate is missing in a site)
    missing_sc <- apply(Xpsi, 1, function(x) any(is.na(x)))
    
    # Missing observation covariates
    # (TRUE if at least one observation covariate is missing in a sampling occasion in a site)
    missing_oc <- apply(Xlambda, 1, function(x) any(is.na(x)))
    
    # Matrix MxJ of values to not use because there is some data missing
    removed_obs <- 
      # If there is count data missing in site i and obs j
      missing_y | 
      # If there is any site covariate missing in site i
      replicate(n = J, missing_sc) |
      # If there is any observation covariate missing in site i and obs j
      matrix(missing_oc, M, J, byrow = T)
    
    if (any(removed_obs)) {
      if (na.rm) {
        nb_missing_sites <- sum(apply(removed_obs, 1, function(x) all(is.na(x))))
        nb_missing_observations <- sum(is.na(removed_obs))
        warning("There are missing data (count data, site covariate, observation covariate):\n\t",
                "Only data from ", M-nb_missing_sites, " sites out of ", M, " are used.\n\t",
                "Only data from ", (M*J)-nb_missing_sites, " observations out of ", (M*J), " are used.")
      } else {
        stop("There is missing data:\n\t",
             sum(missing_y), " missing count data (y)\n\t",
             sum(missing_sc), " missing site covariates (siteCovs)\n\t",
             sum(missing_oc), " missing observation covariates (obsCovs)")
      }
    }
    
    # Output -------------------------------------------------------------------
    return(list(
        y = y,
        Xpsi = Xpsi,
        Zpsi = Zpsi,
        Xlambda = Xlambda,
        Zlambda = Zlambda,
        removed_obs = removed_obs
    ))
  })


## show method ----
setMethod("show", "unmarkedFrameCOP", function(object) {
  J <- ncol(object@obsLength)
  df_unmarkedFrame <- as(object, "data.frame")
  df_obsLength <- data.frame(object@obsLength)
  colnames(df_obsLength) <- paste0("obsLength.", 1:J)
  if (ncol(df_unmarkedFrame) > J) {
    df <-
      cbind(df_unmarkedFrame[, 1:J], df_obsLength, df_unmarkedFrame[, (J + 1):ncol(df_unmarkedFrame)])
  } else{
    df <-
      cbind(df_unmarkedFrame[, 1:J], df_obsLength)
  }
  cat("Data frame representation of unmarkedFrame object.\n")
  print(df)
})
# LP note: as is defined in unmarkedFrame.R part "COERCION"


## summary method ----
setMethod("summary", "unmarkedFrameCOP", function(object,...) {
  cat("unmarkedFrameCOP Object\n\n")
  
  cat(nrow(object@y), "sites\n")
  cat("Maximum number of sampling occasions per site:",obsNum(object),"\n")
  mean.obs <- mean(rowSums(!is.na(getY(object))))
  cat("Mean number of sampling occasions per site:",round(mean.obs,2),"\n")
  cat("Sites with at least one detection:",
      sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  
  cat("Tabulation of y observations:")
  print(table(object@y, exclude=NULL))
  
  if(!is.null(object@obsLength)) {
    cat("\nTabulation of sampling occasions length:")
    print(table(object@obsLength))
  }

  if(!is.null(object@siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object@siteCovs))
  }
  
  if(!is.null(object@obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object@obsCovs))
  }

})

## umf[i, j] ----
setMethod("[", c("unmarkedFrameCOP", "numeric", "numeric", "missing"),
          function(x, i, j) {
            # Gey dimensions of x
            M <- numSites(x)
            J <- obsNum(x)
            
            if (length(i) == 0 & length(j) == 0) {
              return(x)
            }
            
            # Check i
            if (any(i < 0) &&
                any(i > 0)) {
              stop("i must be all positive or all negative indices.")
            }
            if (all(j < 0)) {
              i <- (1:M)[i]
            }
            
            # Check j
            if (any(j < 0) &&
                any(j > 0)) {
              stop("j must be all positive or all negative indices.")
            }
            if (all(j < 0)) {
              j <- (1:J)[j]
            }
            
            # y observation count data subset
            y <- getY(x)[i, j]
            if (min(length(i), length(j)) == 1) {
              y <- t(y)
            }
            
            # obsLength subset
            obsLength <- x@obsLength[i, j]
            if (min(length(i), length(j)) == 1) {
              obsLength <- t(obsLength)
            }
            
            # siteCovs subset
            siteCovs <- siteCovs(x)
            if (!is.null(siteCovs)) {
              siteCovs <- siteCovs(x)[i, , drop = FALSE]
            }
            
            # obsCovs subset
            obsCovs <- obsCovs(x)
            if (!is.null(obsCovs)) {
              MJ_site <- rep(1:M, each = J)
              MJ_obs <- rep(1:J, times = M)
              obsCovs <- obsCovs[((MJ_obs %in% j) & (MJ_site %in% i)), ]
              rownames(obsCovs) <- NULL
            }
            
            # Recreate umf
            new(
              Class = "unmarkedFrameCOP",
              y = y,
              obsLength = obsLength,
              siteCovs = siteCovs,
              obsCovs = obsCovs,
              obsToY = diag(length(j)),
              mapInfo = x@mapInfo
            )
          })


## umf[i, ] ----
setMethod("[", c("unmarkedFrameCOP", "numeric", "missing", "missing"),
          function(x, i) {
            x[i, 1:obsNum(x)]
          })

## umf[, j] ----
setMethod("[", c("unmarkedFrameCOP", "missing", "numeric", "missing"),
          function(x, j) {
            x[1:numSites(x), j]
          })


## fl_getY ----
setMethod("fl_getY", "unmarkedFitCOP", function(fit, ...){
  getDesign(getData(fit), fit@formlist)$y
})

# PREDICT METHODS ----
# Fit type-specific methods to generate different components of prediction
# 5. predict_by_chunk(): Take inputs and generate predictions
# Basic methods are shown below; fit type-specific methods in their own sections

## predict_inputs_from_umf ----
setMethod("predict_inputs_from_umf", "unmarkedFitCOP",
          function(object, type, newdata, na.rm, re.form = NULL) {
            designMats = getDesign(umf = newdata,
                                   formlist = object@formlist,
                                   na.rm = na.rm)
            if (type == "psi") list_els <- c("Xpsi", "Zpsi")
            if (type == "lambda") list_els <- c("Xlambda", "Zlambda")
            X <- designMats[[list_els[1]]]
            if (is.null(re.form)) X <- cbind(X, designMats[[list_els[2]]])
            return(list(X = X, offset = NULL))
          })

## get_formula ----
setMethod("get_formula", "unmarkedFitCOP", function(object, type, ...) {
  fl <- object@formlist
  switch(type, psi = fl$psiformula, lambda = fl$lambdaformula)
})

## get_orig_data ----
setMethod("get_orig_data", "unmarkedFitCOP", function(object, type, ...){
  clean_covs <- clean_up_covs(object@data, drop_final=FALSE)
  datatype <- switch(type, psi = 'site_covs', lambda = 'obs_covs')
  clean_covs[[datatype]]
})



# Useful functions -------------------------------------------------------------

replaceinf = function(x) {
  max(.Machine$double.xmin, min(.Machine$double.xmax, x))
}

check.integer <- function(x, na.ignore = F) {
  if (is.na(x) & na.ignore) {
    return(T)
  }
  x %% 1 == 0
}

# unmarkedFrame ----------------------------------------------------------------

unmarkedFrameCOP <- function(y, obsLength, siteCovs = NULL, obsCovs = NULL, mapInfo = NULL) {
  
  # Verification that they are non-NA data in y
  if (all(is.na(y))) {
    stop("y only contains NA. y needs to contain non-NA integers.")
  }
  
  # Verification that these are count data (and not detection/non-detection)
  if (max(y, na.rm = T) == 1) {
    warning("unmarkedFrameCOP is for count data. ",
            "The data furnished appear to be detection/non-detection.")
  }
  
  # Number of sampling occasions
  J <- ncol(y)
  
  # If missing obsLength: replace by matrix of 1
  # and the lambda will be the detection rate per observation length
  if (missing(obsLength)) {
    obsLength <- y * 0 + 1
    warning("obsLength is missing, replacing it by a matrix of 1.")
  } else if (is.null(obsLength)) {
    obsLength <- y * 0 + 1
    warning("obsLength is missing, replacing it by a matrix of 1.")
  }
  
  # Transformation observation covariates
  obsCovs <- covsToDF(
    covs = obsCovs,
    name = "obsCovs",
    obsNum = J,
    numSites = nrow(y)
  )
  
  # Create S4 object of class unmarkedFrameCOP
  umf <- new(
    Class = "unmarkedFrameCOP",
    y = y,
    obsLength = obsLength,
    siteCovs = siteCovs,
    obsCovs = obsCovs,
    obsToY = diag(J),
    mapInfo = mapInfo
  )
  
  return(umf)
}


# occuCOP ----------------------------------------------------------------------

occuCOP <- function(data,
                    psiformula = ~1,
                    lambdaformula = ~1,
                    psistarts = rep(0, length(attr(terms(psiformula), "term.labels")) + 1),
                    lambdastarts = rep(0, length(attr(terms(lambdaformula), "term.labels")) + 1),
                    method = "BFGS",
                    se = TRUE,
                    engine = "R",
                    threads = 1L,
                    na.rm = TRUE,
                    get.NLL.params = NULL,
                    ...) {
  #TODO: engines
  #TODO: NA
  
  # Neg loglikelihood COP ------------------------------------------------------
  nll_COP <- function(params) {
    
    # Reading and transforming params
    # Taking into account the covariates
    
    # psi as a function of covariates
    #   psi in params are transformed with a logit transformation (qlogis)
    #       so they're back-transformed to a proba with inverse logit (plogis)
    #   Xpsi is the matrix with occupancy covariates
    #   params is the vector with all the parameters
    #   psiIdx is the index of Occupancy Parameters in params
    psi <- plogis(Xpsi %*% params[psiIdx])
    
    # lambda as a function of covariates
    #   lambda in params are transformed with a log-transformation
    #          so they're back-transformed to a rate with exp here
    #   Xlambda is the matrix with detection covariates
    #   params is the vector with all the parameters
    #   lambdaIdx is the index of Occupancy Parameters in params
    lambda <- exp(Xlambda %*% params[lambdaIdx])
    
    # Probability for each site (i)
    iProb <- rep(NA, M)
    
    for (i in 1:M) {
      # Probability for the site i for each sampling occasion (j)
      ijProb <- rep(NA, J)
      
      # Total obs length in site i
      obsLengthSitei <- sum(data@obsLength[i,])
      
      # iIdx is the index to access the vectorised vectors of all obs in site i
      iIdx <- ((i - 1) * J + 1):(i * J)
      
      if (SitesWithDetec[i]) {
        # If there is at least one detection in site i
        # factNij = sapply(factorial(yvec[iIdx]),replaceinf)
        # iProb[i] = psi[i] * sum(lambda[iIdx] * obsLengthvec[iIdx] / factNij * exp(-lambda[iIdx] * obsLengthvec[iIdx]))
        
        
        iProb[i] = psi[i] * (
          (sum(lambda[iIdx] * obsLengthvec[iIdx])) ^ sum(yvec[iIdx]) / 
            factorial(sum(yvec[iIdx])) *
            exp(-sum(lambda[iIdx] * obsLengthvec[iIdx]))
        )
        
      } else {
        # If there is zero detection in site i
        # iProb[i] = psi[i] * sum(exp(-lambda[iIdx] * obsLengthvec[iIdx])) + (1 - psi[i])
        
        iProb[i] = psi[i] * exp(-sum(lambda[iIdx] * obsLengthvec[iIdx])) + (1 - psi[i]) 
        
      }
      
    }
    # Note: Why is there "replaceinf(factorial(data@y[i,]))" 
    #       instead of just "factorial(data@y[i,])"?
    # Because if there are a lot of detections (n_ij >= 171 on my machine),
    #   then factorial(n_ij) = Inf, 
    #   then 1/factorial(n_ij) = 0, 
    #   then log(0) = -Inf,
    #   then loglikelihood = -Inf
    #   then optim does not work.
    # 
    # We want to maximise likelihood, 
    # (technically, minimise negative loglikelihood).
    # We don't do anything else with this function.
    # So as long as we have a tiny value for likelihood, 
    # (ie a huge value for negative loglikelihood)
    # it doesn't matter if its an approximation.
    
    # log-likelihood
    ll = sum(log(iProb))
    return(-ll)
  }
  
  # Check arguments ------------------------------------------------------------
  if (!is(data, "unmarkedFrameCOP")) {
    stop("Data is not an unmarkedFrameCOP object. See ?unmarkedFrameCOP if necessary.")
  }
  stopifnot(class(psiformula) == "formula")
  stopifnot(class(lambdaformula) == "formula")
  stopifnot(class(psistarts) %in%  c("numeric", "double", "integer"))
  stopifnot(class(lambdastarts) %in%  c("numeric", "double", "integer"))
  stopifnot(class(threads) %in%  c("numeric", "double", "integer"))
  se = as.logical(match.arg(
    arg = as.character(se),
    choices = c("TRUE", "FALSE", "0", "1")
  ))
  na.rm = as.logical(match.arg(
    arg = as.character(na.rm),
    choices = c("TRUE", "FALSE", "0", "1")
  ))
  engine <- match.arg(engine, c("R"))
  
  # Format input data ----------------------------------------------------------
  
  # Retrieve formlist
  formlist <- mget(c("psiformula", "lambdaformula"))
  
  # Get the design matrix (calling the getDesign method for unmarkedFrame)
  # For more informations, see: getMethod("getDesign", "unmarkedFrameCOP")
  designMats <- getDesign(umf = data, formlist = formlist, na.rm = na.rm)
  
  # y is the count detection data (matrix of size M sites x J sessions)
  y <- getY(data)
  
  # obsLength is the length of sessions (matrix of size M sites x J sessions)
  obsLength = data@obsLength
  
  # Xpsi is the matrix with the occupancy covariates
  Xpsi <- designMats$Xpsi
  
  # Xlambda is the matrix with the detection covariates
  Xlambda <- designMats$Xlambda
  
  # Zpsi is ???
  Zpsi <- designMats$Zpsi
  
  # Zlambda is ???
  Zlambda <- designMats$Zlambda
  
  # removed_obs is a M x J matrix of the observations removed from the analysis
  removed_obs <- designMats$removed_obs
  sitesRemoved <- unname(which(apply(removed_obs, 1, function(x) all(x))))
  
  # Number of sites
  M <- nrow(y)
  
  # Number of sampling occasions
  J <- ncol(y)
  
  # Total number of detection per site
  NbDetecPerSite = rowSums(y)
  
  # Sites where there was at least one detection
  SitesWithDetec = NbDetecPerSite > 0
  
  # Number of sites where there was at least one detection
  NbSitesWithDetec = sum(SitesWithDetec)
  
  # Set up parameter names and indices-----------------------------------------
  
  # ParamPsi Occupancy parameter names
  ParamPsi <- colnames(Xpsi)
  
  # ParamLambda Detection parameter names
  ParamLambda <- colnames(Xlambda)
  
  # NbParamPsi Number of occupancy parameters
  NbParamPsi <- ncol(Xpsi)
  
  # NbParamLambda Number of detection parameters
  NbParamLambda <- ncol(Xlambda)
  
  # nP Total number of parameters
  nP <- NbParamPsi + NbParamLambda
  
  # psiIdx Index of the occupancy parameters
  psiIdx <- 1:NbParamPsi
  
  # lambdaIdx Index of the detection parameters
  lambdaIdx <- (NbParamPsi+1):nP
  
  # Re-format some variables for C and R engines
  #   From Matrix of dim MxJ to vector of length MxJ:
  #   c(ySite1Obs1, ySite1Obs2, ..., ySite1ObsJ, ysite2Obs1, ...)
  yvec <- as.numeric(t(y))
  obsLengthvec <- as.numeric(t(obsLength))
  
  # get.NLL.params -------------------------------------------------------------
  if (!is.null(get.NLL.params)) {
    df_NLL = data.frame(t(as.data.frame(get.NLL.params)))
    rownames(df_NLL) = NULL
    colnames(df_NLL) = c(paste0("logit(psi).", ParamPsi),
                         paste0("log(lambda).", ParamLambda))
    df_NLL$nll = NA
    for (i in 1:nrow(df_NLL)) {
      df_NLL$nll[i] = nll_COP(params = as.numeric(as.vector(df_NLL[i, -ncol(df_NLL)])))
    }
    return(df_NLL)
  }
  
  # Optimisation ---------------------------------------------------------------
  
  ## Checking the starting point for optim
  if (missing(lambdastarts)) {
    # If lambda starts argument was not given:
    # 0 by default 
    # so lambda = exp(0) = 1 by default
    lambdastarts = rep(0, NbParamLambda)
    message(
      "No lambda initial values provided for optim. Using lambdastarts = c(",
      paste(lambdastarts, collapse = ", "),
      "), equivalent to detection rate",
      ifelse(length(lambdastarts) == 1, "", "s"),
      " of 1."
    )
  } else if (length(lambdastarts) != NbParamLambda) {
    stop("lambdastarts (", paste(lambdastarts, collapse = ", "), ") ",
         "should be of length ", NbParamLambda, " with lambdaformula ", lambdaformula)
  }
  
  if (missing(psistarts)) {
    # If psi starts argument was not given
    # 0 by default 
    # so psi = plogis(0) = 0.5 by default
    psistarts = rep(0, NbParamPsi)
    message(
      "No psi initial values provided for optim. Using psistarts = c(",
      paste(psistarts, collapse = ", "),
      "), equivalent to occupancy probabilit",
      ifelse(length(psistarts) == 1, "y", "ies"),
      " of 0.5."
    )
  } else if (length(psistarts) != NbParamPsi) {
    stop("psistarts (", paste(psistarts, collapse = ", "), ") ",
         "should be of length ", NbParamPsi, " with psiformula ", psiformula)
  }

  starts <- c(psistarts, lambdastarts)
  
  ## Run optim
  opt <- optim(
    starts,
    nll_COP,
    method = method,
    hessian = se,
    ...
  )
  
  # Get output -----------------------------------------------------------------
  covMat <- invertHessian(opt, nP, se)
  ests <- opt$par
  tmb_mod <- NULL
  nll <- opt$value
  fmAIC <- 2 * opt$value + 2 * nP
  
  # Organize effect estimates
  names(ests) <- c(ParamPsi, ParamLambda)
  psi_coef <- list(ests = ests[psiIdx], cov = as.matrix(covMat[psiIdx, psiIdx]))
  lambda_coef <- list(ests = ests[lambdaIdx],
                      cov = as.matrix(covMat[lambdaIdx, lambdaIdx]))
  
  # No random effects
  psi_rand_info <- lambda_rand_info <- list()
  
  # Create unmarkedEstimates ---------------------------------------------------
  psi_uE <- unmarkedEstimate(
    name = "Occupancy probability",
    short.name = "psi",
    estimates = psi_coef$ests,
    covMat = psi_coef$cov,
    fixed = 1:NbParamPsi,
    invlink = "logistic",
    invlinkGrad = "logistic.grad",
    randomVarInfo = psi_rand_info
  )
  
  lambda_uE <- unmarkedEstimate(
    name = "Detection rate",
    short.name = "lambda",
    estimates = lambda_coef$ests,
    covMat = lambda_coef$cov,
    fixed = 1:NbParamLambda,
    invlink = "exp",
    invlinkGrad = "exp",
    randomVarInfo = lambda_rand_info
  )
  
  estimateList <- unmarkedEstimateList(list(psi = psi_uE, lambda = lambda_uE))
  
  # Create unmarkedFit object--------------------------------------------------
  umfit <- new(
    "unmarkedFitCOP",
    fitType = "occuCOP",
    call = match.call(),
    # formula = as.formula(paste(formlist, collapse = "")),
    formula = as.formula(paste(
      formlist["lambdaformula"], formlist["psiformula"], collapse = ""
    )),
    formlist = formlist,
    data = data,
    estimates = estimateList,
    sitesRemoved = sitesRemoved,
    removed_obs = removed_obs,
    AIC = fmAIC,
    opt = opt,
    negLogLike = opt$value,
    nllFun = nll_COP,
    nll = nll,
    convergence = opt$convergence,
    TMB = tmb_mod
  )
  
  return(umfit)
}

