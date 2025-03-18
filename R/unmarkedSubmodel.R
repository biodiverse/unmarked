setClass("unmarkedSubmodel",
  slots = c(
    type = "character",
    data = "data.frame",
    formula = "formula",
    family = "character",
    auxiliary = "list"
  ),
  prototype = list(
    type = NA_character_,
    data = data.frame(),
    formula = ~1,
    family = NA_character_,
    auxiliary = list()
  ),
  contains = "unmarkedEstimate"
)

setClass("unmarkedSubmodelState",
  contains = "unmarkedSubmodel")

setClass("unmarkedSubmodelDet",
  contains = "unmarkedSubmodel")

unmarkedSubmodelDet <- function(name, short_name, type, formula, data, 
                                family, link, auxiliary = list()){
  data <- clean_up_covs(data, drop_final = FALSE)$obs_covs
  data <- subset_covs(data, formula)
  out <- new("unmarkedSubmodelDet",
    name = name, short.name = short_name, type = type,
    formula = formula, data = data, family = family,
    invlink = get_invlink(link), invlinkGrad = get_grad(link),
    auxiliary = auxiliary)
  out@fixed <- 1:ncol(model.matrix(out))
  out
}

unmarkedSubmodelState <- function(name, short_name, type, formula, data, 
                                  family, link, auxiliary = list()){
  data <- clean_up_covs(data, drop_final = FALSE)$site_covs
  data <- subset_covs(data, formula)
  out <- new("unmarkedSubmodelState",
    name = name, short.name = short_name, type = type,
    formula = formula, data = data, family = family,
    invlink = get_invlink(link), invlinkGrad = get_grad(link),
    auxiliary = auxiliary)
  out@fixed <- 1:ncol(model.matrix(out))
  out
}


subset_covs <- function(covs, formula){
  vars <- all.vars(formula)
  cov_nms <- names(covs)
  miss <- vars[! vars %in% cov_nms]
  if(length(miss) > 0){
    stop(paste("Variable(s)", paste(miss, collapse = ", "),
               "not in covariates"), call.=FALSE)
  }
  covs[,vars,drop=FALSE]
}

get_invlink <- function(link){
  switch(link, 
         logit = "logistic",
         log = "exp",
         cloglog = "cloglog"
        )
}

get_grad <- function(link){
  switch(link, 
         logit = "logistic.grad",
         log = "exp",
         cloglog = "cloglog.grad"
        )
}


setMethod("model.matrix", "unmarkedSubmodel",
  function(object, newdata = NULL, ...){
  form <- reformulas::nobars(object@formula)
  mf <- model.frame(form, object@data, na.action = stats::na.pass)
  
  if(!is.null(newdata)){
    X_terms <- stats::terms(mf)
    xlevs <- list()
    if(ncol(object@data) > 0){
      fac_cols <- sapply(object@data, is.factor)
      fac_cols <- object@data[, sapply(object@data, is.factor), drop=FALSE]
      xlevs <- lapply(fac_cols, levels)
      xlevs <- xlevs[names(xlevs) %in% names(mf)]
    }
    mf <- model.frame(X_terms, newdata, na.action=stats::na.pass, xlev=xlevs)
  }

  model.matrix(form, mf)
})

setGeneric("has_random", function(object) standardGeneric("has_random"))

setMethod("has_random", "formula", function(object){
  length(reformulas::findbars(object)) > 0
})

setMethod("has_random", "unmarkedSubmodel", function(object){
  has_random(object@formula) 
})

setGeneric("Z_matrix", function(object, ...){
  standardGeneric("Z_matrix")
})

setMethod("Z_matrix", "unmarkedSubmodel",
  function(object, newdata = NULL, ...){
  if(!has_random(object)){
    return(as(as(Matrix::Matrix(0,0,0), "generalMatrix"), "CsparseMatrix"))
  }
  get_Z(object@formula, object@data, newdata)
})

setGeneric("get_offset", function(object, ...){
  standardGeneric("get_offset")
})

setMethod("get_offset", "unmarkedSubmodel",
  function(object, newdata = NULL, ...){
  mf <- model.frame(reformulas::nobars(object@formula), object@data, na.action=stats::na.pass)
  if(!is.null(newdata)){
    X_terms <- stats::terms(mf)
    mf <- model.frame(X_terms, newdata, na.action=stats::na.pass)
  }
  off <- as.vector(stats::model.offset(mf))
  if(is.null(off)) off <- rep(0, nrow(mf))
  off
})

setGeneric("find_missing", function(object, ...){ 
  standardGeneric("find_missing")
})

setMethod("find_missing", "unmarkedSubmodel",
  function(object, ...){
  mm <- model.matrix(object)
  apply(mm, 1, function(x) any(is.na(x)))
})

setGeneric("default_starts", function(object, ...) standardGeneric("default_starts"))

setMethod("default_starts", "unmarkedSubmodel", function(object, ...){
  nms <- paste0(object@short.name, "(", colnames(model.matrix(object)), ")")
  nms <- gsub("(Intercept)", "Int", nms, fixed = TRUE)
  stats::setNames(rep(0, length(nms)), nms)
})

setGeneric("engine_inputs", function(object, object2) standardGeneric("engine_inputs"))

setMethod("engine_inputs", c("unmarkedSubmodel", "missing"), function(object, object2){
  out <- list(family = get_family_code(object@family),
              invlink = get_invlink_code(object@invlink),
              X = model.matrix(object), 
              offset = get_offset(object))
  out <- c(out, object@auxiliary)
  names(out) <- paste0(names(out), "_", object@type)
  out 
})

setGeneric("engine_inputs_TMB", function(object, object2) standardGeneric("engine_inputs_TMB"))

setMethod("engine_inputs_TMB", c("unmarkedSubmodel", "missing"), function(object, object2){
  out <- engine_inputs(object)
  n_grouplevels <- get_nrandom(object@formula, object@data)
  tmb <- list(Z = Z_matrix(object),
              n_group_vars = length(n_grouplevels[n_grouplevels > 0]),
              n_grouplevels = n_grouplevels)
  names(tmb) <- paste0(names(tmb), "_", object@type)
  c(out, tmb)
})

get_family_code <- function(family){
  switch(family, binomial = 0, poisson = 1)
}

get_invlink_code <- function(link){
  switch(link, logistic = 0, exp = 1, cloglog = 2)
}

setGeneric("get_TMB_pars", function(object, starts){
  standardGeneric("get_TMB_pars")
})

setMethod("get_TMB_pars", "unmarkedSubmodel", function(object, starts){

  out <- list()

  # Fixed pars
  fixed <- colnames(model.matrix(object))
  out$beta <- stats::setNames(rep(0, length(fixed)), fixed)

  # Random effect names
  if(has_random(object)){
    re_info <- random_effect_names(object)

    # Sigma pars
    sig_names <- re_info$sigma
    sig_names <- apply(sig_names, 1, paste, collapse = "_")
    sig_names <- gsub("(Intercept)", "Int", sig_names, fixed = TRUE)
    out$lsigma <- stats::setNames(rep(0, length(sig_names)), sig_names)

    # Random effect pars
    b_names <- apply(re_info$random_effects, 1, paste, collapse = "_")
    b_names <- gsub("(Intercept)", "Int", b_names, fixed = TRUE)
    out$b <- stats::setNames(rep(0, length(b_names)), b_names)
  } else {
    out$lsigma <- numeric(0)
    out$b <- numeric(0)
  }
  
  names(out) <- paste0(names(out), "_", object@type)
  
  out
})

random_effect_names <- function(object){

  re_info <- get_reTrms(object@formula, object@data)[c("cnms", "flist")]

  Groups <- lapply(1:length(re_info$cnms), function(x){
                  gn <- names(re_info$cnms)[x]
                  rep(gn, length(levels(re_info$flist[[gn]])))
            })
  Groups <- do.call(c, Groups)

  Name <-  lapply(1:length(re_info$cnms), function(x){
                  gn <- names(re_info$cnms)[x]
                  var <- re_info$cnms[[x]]
                  rep(var, length(levels(re_info$flist[[gn]])))
            })
  Name <- do.call(c, Name)

  Levels <- lapply(1:length(re_info$cnms), function(x){
              gn <- names(re_info$cnms)[x]
              levels(re_info$flist[[gn]])
            })
  Levels <- do.call(c, Levels)

  re <- data.frame(Groups = Groups, Name = Name, Levels = Levels)
 
  list(sigma = unique(re[,c("Groups", "Name")]),
       random_effects = re)
}

setGeneric("get_TMB_random", function(object){
  standardGeneric("get_TMB_random")
})

setMethod("get_TMB_random", "unmarkedSubmodel", 
  function(object){
  if(!has_random(object)) return(NULL)
  paste0("b_", object@type)
})

# unmarkedSubmodelList---------------------------------------------------------

setClass("unmarkedSubmodelList",
  contains = "unmarkedEstimateList"
)

unmarkedSubmodelList <- function(...){
  submodels <- list(...)
  names(submodels) <- sapply(submodels, function(x) x@type)
  new("unmarkedSubmodelList", estimates=submodels)
}

setGeneric("submodels", function(object) standardGeneric("submodels"))

setMethod("submodels", "unmarkedSubmodelList", function(object){
  object@estimates
})

setGeneric("submodels<-", function(object, value) standardGeneric("submodels<-"))

setMethod("submodels<-", "unmarkedSubmodelList", function(object, value){
  object@estimates <- value
  object
})

setMethod("has_random", "unmarkedSubmodelList", function(object){
  sapply(submodels(object), has_random)
})

setGeneric("get_parameter_idx", function(object, ...){
  standardGeneric("get_parameter_idx")
})

setMethod("get_parameter_idx", "unmarkedSubmodelList",
  function(object, ...){
  idx <- lapply(submodels(object), function(x){
    1:ncol(model.matrix(x))
  })
  for (i in 2:length(idx)){
    idx[[i]] <- idx[[i]] + max(idx[[i-1]])
  }
  lapply(idx, range)
})

setMethod("engine_inputs", c("unmarkedSubmodelList", "missing"), function(object, object2){
  out <- lapply(submodels(object), engine_inputs)
  do.call(c, unname(out))
})

setMethod("engine_inputs_TMB", c("unmarkedSubmodelList", "missing"), function(object, object2){
  out <- lapply(submodels(object), engine_inputs_TMB)
  do.call(c, unname(out))
})

setMethod("default_starts", "unmarkedSubmodelList", function(object, ...){
  unlist(unname(lapply(submodels(object), default_starts)))
})

setMethod("get_TMB_random", "unmarkedSubmodelList",
  function(object){
  unlist(sapply(submodels(object), get_TMB_random))
})

setMethod("get_TMB_pars", "unmarkedSubmodelList", function(object, starts){
  unlist(lapply(unname(submodels(object)), get_TMB_pars, starts = starts), 
         recursive = FALSE)
})

# unmarkedResponse-------------------------------------------------------------

setClass("unmarkedResponse",
  slots = c(
    y = "ANY",
    Kmin = "numeric",
    Kmax = "numeric"
  ), 
  prototype = list(
    y = numeric(0),
    Kmin = numeric(0),
    Kmax = numeric(0)
  )
)

unmarkedResponse <- function(data, Kmax = 1){
  response <- new("unmarkedResponse", y = data@y, Kmax = Kmax)
  #response <- add_missing(response, submodels)
  Kmin <- apply(data@y, 1, function(x){
    ifelse(all(is.na(x)), NA, max(x, na.rm=TRUE))
  })
  response@Kmin <- Kmin
  response
}

setGeneric("add_missing", function(object, submodel, ...){
  standardGeneric("add_missing")
})

setMethod("add_missing", c("unmarkedResponse", "unmarkedSubmodelState"),
  function(object, submodel, ...){
  
  object@y[find_missing(submodel),] <- NA
  object
})

setMethod("add_missing", c("unmarkedResponse", "unmarkedSubmodelDet"),
  function(object, submodel, ...){
  yt <- t(object@y)
  yt[find_missing(submodel)] <- NA
  object@y <- t(yt)
  Kmin <- apply(object@y, 1, function(x){
    ifelse(all(is.na(x)), NA, max(x, na.rm=TRUE))
  })
  object@Kmin <- Kmin
  object
})

setMethod("add_missing", c("unmarkedResponse", "unmarkedSubmodelList"),
  function(object, submodel, ...){  
  for (i in 1:length(submodels(submodel))){
    object <- add_missing(object, submodels(submodel)[[i]])
  }
  object
})

setGeneric("removed_sites", function(object, ...){
  standardGeneric("removed_sites")
})

setMethod("removed_sites", "unmarkedResponse", function(object, ...){
  all_na <- apply(object@y, 1, function(x) all(is.na(x)))
  which(all_na)
})

setMethod("engine_inputs", c("unmarkedResponse", "missing"), function(object, object2){
  list(y = object@y, Kmin = object@Kmin, Kmax = object@Kmax)
})


setMethod("engine_inputs", c("unmarkedResponse", "unmarkedSubmodelList"),
  function(object, object2){
  c(engine_inputs(object), engine_inputs(object2))
})


setGeneric("engine_inputs_CR", function(object, object2) standardGeneric("engine_inputs_CR"))

setMethod("engine_inputs_CR", c("unmarkedResponse", "unmarkedSubmodelList"), 
  function(object, object2){
  out <- engine_inputs(object, object2)
  idx <- get_parameter_idx(object2)
  names(idx) <- paste0("idx_", names(idx))
  c(out, idx)
})

setMethod("engine_inputs_TMB", c("unmarkedResponse", "unmarkedSubmodelList"), 
  function(object, object2){
  c(engine_inputs(object), engine_inputs_TMB(object2))
})


# Fitting models---------------------------------------------------------------

fit_model <- function(nll_fun, response, submodels, starts, method, se, ...){
  # C and R engines
  if(is.function(nll_fun)){ 
    fit <- fit_optim(nll_fun, response = response, submodels = submodels,
                     starts = starts, method = method, se = se, ...) 
  } else if(grepl("tmb_", nll_fun)){
    fit <- fit_TMB2(nll_fun, response = response, submodels = submodels, 
                    starts = starts, method = method, se = se, ...)
  }
  fit
}

fit_optim <- function(nll_fun, response, submodels, starts, method, se, ...){
  inputs <- engine_inputs_CR(response, submodels)

  npars <- max(unlist(get_parameter_idx(submodels)))
  if(is.null(starts)) starts <- default_starts(submodels)
  if(length(starts) != npars){
    stop(paste("The number of starting values should be", npars))
  }

  # Capture the environment for profile likelihood
  nll <- function(params) nll_fun(params, inputs)

  opt <- optim(starts, nll, method = method, hessian = se, ...)            
  cov_mat <- invertHessian(opt, length(opt$par), se)

  fit <- list(opt = opt, cov_mat = cov_mat, TMB = NULL, nll = nll)
  submodels <- add_estimates(submodels, fit)

  fit$submodels <- submodels
  fit$AIC <- 2 * opt$value + 2 * length(opt$par)
  fit
}

fit_TMB2 <- function(model, response, submodels, starts, method, se, ...){
  data <- engine_inputs_TMB(response, submodels)
  params <- get_TMB_pars(submodels, starts = starts)
  random <- get_TMB_random(submodels)

  fixed_sub <- names(params)[!names(params) %in% random]
  nfixed <- length(unlist(params[fixed_sub]))
  list_fixed_only <- params[fixed_sub]
  plengths <- sapply(list_fixed_only, length)
  starts_order <- rep(fixed_sub, plengths)

  if(!is.null(starts)){
    if(length(starts) != nfixed){
      stop(paste("The number of starting values should be", nfixed))
    }
    list_fixed_only <- params[fixed_sub]
    list_fixed_only <- utils::relist(starts, list_fixed_only)
    params <- replace(params, names(list_fixed_only), list_fixed_only)
  }

  tmb_mod <- TMB::MakeADFun(data = c(model = model, data),
                            parameters = params,
                            random = random,
                            silent=TRUE,
                            DLL = "unmarked_TMBExports")
  tmb_mod$starts_order <- starts_order
  class(tmb_mod) <- "TMB"

  opt <- optim(tmb_mod$par, fn=tmb_mod$fn, gr=tmb_mod$gr, method=method, ...)

  submodels <- add_estimates(submodels, tmb_mod, se)
  
  AIC = 2 * opt$value + 2 * length(opt$par)
  
  list(opt=opt, TMB=tmb_mod, nll = tmb_mod$fn, submodels = submodels,
       AIC = AIC)
}


# Add estimates to submodels---------------------------------------------------

setGeneric("add_estimates", function(object, fit, ...) standardGeneric("add_estimates"))

# optim method (R/C engines)
setMethod("add_estimates", c("unmarkedSubmodelList", "list"),
  function(object, fit, ...){
  idx <- get_parameter_idx(object)
  for (i in names(submodels(object))){
    idx_rng <- idx[[i]][1]:idx[[i]][2] 
    est <- fit$opt$par[idx_rng]
    names(est) <- colnames(model.matrix(object@estimates[[i]]))
    object@estimates[[i]]@estimates <- est
    v <- fit$cov_mat[idx_rng, idx_rng, drop = FALSE]
    rownames(v) <- colnames(v) <- names(est) 
    object@estimates[[i]]@covMat <- v
  }
  object
})

setOldClass("TMB")

setMethod("add_estimates", c("unmarkedSubmodel", "TMB"), 
  function(object, fit, sdr, ...){   
  
  # Add estimates
  pars <- fit$env$last.par.best
  fixed_type <- grepl(paste0("beta_", object@type), names(pars))
  names(pars)[fixed_type] <- colnames(model.matrix(object))

  rand_type <- grepl(paste0("b_", object@type), names(pars))

  beta_or_b <- fixed_type | rand_type
  object@estimates <- pars[beta_or_b]

  # Add variance-covariance matrix
  # default is blank
  covMat <- matrix(NA, sum(beta_or_b), sum(beta_or_b))
  if(!is.null(sdr)){ # If there's an SD report
    if(!is.null(sdr$cov_all)){
      # If there's a joint precision, use it
      covMat <- sdr$cov_all[beta_or_b, beta_or_b, drop = FALSE]
    } else {
      # Otherwise we can only get the fixed effects covmat
      fixed_type <- grepl(paste0("beta_", object@type), colnames(sdr$cov.fixed))
      covMat <- sdr$cov.fixed[fixed_type, fixed_type, drop = FALSE]
    }
  }
  object@covMat <- covMat

  # Random var info - this should be done more efficiently
  if(has_random(object)){
    sig_type <- grepl(paste0("lsigma_", object@type), names(pars))
    sig_est <- pars[sig_type]
    
    sig_cov <- matrix(NA, sum(sig_type), sum(sig_type))
    if(!is.null(sdr)){
      sig_type <- grepl(paste0("lsigma_", object@type), colnames(sdr$cov.fixed))
      sig_cov <- sdr$cov.fixed[sig_type, sig_type, drop = FALSE]
    }
    
    re <- get_reTrms(object@formula, object@data)
    Z <- get_Z(object@formula, object@data) # inefficient!!
    
    rv_info <- list(names=sigma_names(object@formula, object@data), estimates=sig_est, 
                    covMat=sig_cov, invlink="exp", invlinkGrad="exp",
                    n_levels=lapply(re$flist, function(x) length(levels(x))), 
                    cnms=re$cnms, levels=colnames(Z))
    object@randomVarInfo <- rv_info
  }

  object
})

setMethod("add_estimates", c("unmarkedSubmodelList", "TMB"), 
  function(object, fit, se, ...){   
  
  sdr <- NULL
  if(se){
    if(any(has_random(object))){
      sdr <- TMB::sdreport(fit, getJointPrecision = TRUE)
      sdr$cov_all <- solve(sdr$jointPrecision)
    } else {
      sdr <- TMB::sdreport(fit)
    }
  }

  object@estimates <- lapply(object@estimates, add_estimates, fit = fit, sdr = sdr)

  object
})
